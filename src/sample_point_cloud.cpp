#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

#include <fstream>
#include <iostream>
#include <functional>

#include "common/common.h"
#include "common/vcg_utils.h"


namespace {

using namespace vcg;
class VCGMeshEdge;
class VCGMeshFace;
class VCGMeshVertex;
struct VCGMeshUsedTypes : public UsedTypes<	Use<VCGMeshVertex>   ::AsVertexType,
                                            Use<VCGMeshEdge>     ::AsEdgeType,
                                            Use<VCGMeshFace>     ::AsFaceType>{};
class VCGMeshVertex  : public Vertex<VCGMeshUsedTypes, vertex::Coord3d, vertex::Normal3d, vertex::BitFlags> {};
class VCGMeshFace    : public Face<VCGMeshUsedTypes, face::FFAdj,  face::Normal3d, face::VertexRef, face::BitFlags> {};
class VCGMeshEdge    : public Edge<VCGMeshUsedTypes>{};
class VCGMesh : public tri::TriMesh<std::vector<VCGMeshVertex>, std::vector<VCGMeshFace>, std::vector<VCGMeshEdge>> {};


template <typename PointScalar, typename AttribScalar>
struct AccumulatedPoint {
    using AttribVecT = Eigen::Matrix<AttribScalar, 1, Eigen::Dynamic>;
    using PointT = Eigen::Matrix<PointScalar, 3, 1>;

    AccumulatedPoint()
            : num_of_points_(0),
              point_(PointScalar(0.0), PointScalar(0.0), PointScalar(0.0)) {}

    void AddPoint(const PointT& point, const AttribVecT& attrib) {
        point_ += point;
        if (attrib.cols() > 0) {
            attrib_ += attrib;
        }
        num_of_points_ += 1;
    }

    PointT GetAveragePoint() const {
        return point_ / double(num_of_points_);
    }

    AttribVecT GetAverageAttrib() const {
        // Call NormalizeNormals() afterwards if necessary
        return attrib_ / AttribScalar(num_of_points_);
    }

    int GetNumPoints() const {
        return num_of_points_;
    }

    int num_of_points_;
    PointT point_;
    AttribVecT attrib_;
};


template <typename T>
struct hash_eigen {
    std::size_t operator()(T const& matrix) const {
        size_t seed = 0;
        for (int i = 0; i < (int)matrix.size(); i++) {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


template <typename DerivedV, typename DerivedAttrib, typename DerivedOutV, typename DerivedOutA>
void downsample_point_cloud_to_voxels(const DerivedV& v,
                                      const DerivedAttrib& v_attrib,
                                      const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& voxel_size3,
                                      const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& voxel_min_bound,
                                      const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& voxel_max_bound,
                                      DerivedOutV& out_v,
                                      DerivedOutA& out_attrib,
                                      int min_pts_per_bin) {
    // The way this is templated is dangerous since you could pass in an Eigen expression template. However,
    // we always pass either an Eigen::Map or Eigen::Matrix so it should be okay.

    for (int i = 0; i < 3; i++) {
        if (voxel_size3[i] <= 0.0) {
            throw pybind11::value_error("Voxel size is negative");
        }
        if (voxel_size3[i] * std::numeric_limits<int>::max() < (voxel_max_bound - voxel_min_bound)[i]) {
            throw pybind11::value_error("Voxel size is too small");
        }
    }

    bool has_attribute = (v_attrib.rows() != 0 && v_attrib.cols() != 0);
    if (has_attribute && v_attrib.rows() != v.rows()) {
        throw pybind11::value_error("Invalid number of attributes (" + std::to_string(v_attrib.rows()) +
                                    "). Must match number of input vertices (" +
                                    std::to_string(v.rows()) + ") or be 0.");
    }

    using AccPointT = AccumulatedPoint<typename DerivedV::Scalar, typename DerivedAttrib::Scalar>;
    using AttribVecT = typename AccPointT::AttribVecT;
    using PointT = typename AccPointT::PointT;
    std::unordered_map<Eigen::Vector3i, AccPointT, hash_eigen<Eigen::Vector3i>> voxelindex_to_accpoint;

    PointT ref_coord;
    PointT ref_point;
    AttribVecT ref_attrib;
    Eigen::Vector3i voxel_index;
    for (int i = 0; i < v.rows(); i++) {
        ref_point = PointT(v(i, 0), v(i, 1), v(i, 2));
        ref_coord = (ref_point - voxel_min_bound).array() / voxel_size3.array();
        voxel_index << int(floor(ref_coord(0))), int(floor(ref_coord(1))), int(floor(ref_coord(2)));
        if (has_attribute) {
            ref_attrib = v_attrib.row(i).eval();
        }
        voxelindex_to_accpoint[voxel_index].AddPoint(ref_point, ref_attrib);
    }

    int num_output_points = voxelindex_to_accpoint.size();
    out_v.resize(num_output_points, 3);
    if (has_attribute) {
        out_attrib.resize(num_output_points, v_attrib.cols());
    }

    int count = 0;
    for (auto accpoint : voxelindex_to_accpoint) {
        if (accpoint.second.GetNumPoints() < min_pts_per_bin) {
            continue;
        }
        PointT accpoint_pos = accpoint.second.GetAveragePoint();
        for (int i = 0; i < 3; i++) {out_v(count, i) = accpoint_pos[i];}

        if (has_attribute) {
            AttribVecT accpoint_attrib = accpoint.second.GetAverageAttrib();
            for (int i = 0; i < accpoint_attrib.cols(); i++) {out_attrib(count, i) = accpoint_attrib[i];}
        }
        count += 1;
    }
    if (count < num_output_points) {
        out_v.conservativeResize(count, v.cols());
        if (has_attribute) {
            out_attrib.conservativeResize(count, v_attrib.cols());
        }
    }
}
}



const char* downsample_point_cloud_poisson_disk_doc = R"Qu8mg5v7(
Downsample a point set so that samples are approximately evenly spaced.
This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

Parameters
----------
v : #v by 3 array of vertex positions
n : #v by 3 array of vertex normals
num_samples: desired number of Poisson Disk samples. Note that the actual number of returned samples
             will not be exactly this value (see sample_num_tolerance) to control the range of possible
             returned samples.
             Note: If this value <= 0, then the parameter radius is used to decide the number of samples
radius : desired separation between points, if num_samples <= 0, then this value is used to determine the
         sampling (-1.0, by default).
best_choice_sampling : When downsampling, always keep the sample that will remove the
                       fewest number of samples. (True by default).
random_seed : A random seed used to generate the samples.
              Passing in 0 will use the current time. (0 by default).
sample_num_tolerance: If you requested a target number of samples, by passsing num_samples > 0, then this function will return
                      between (1 - sample_num_tolerance) * num_samples and (1 + sample_num_tolerance) * num_samples.
                      Setting a very small value for this parameter will increase convergence time. (0.04 by default).

Returns
-------
A (m,) shaped array of indices into v where m is the number of Poisson-disk samples

)Qu8mg5v7";
npe_function(downsample_point_cloud_poisson_disk)
    npe_arg(v, dense_float, dense_double)
    npe_arg(num_samples, int)
    npe_default_arg(radius, double, 0.0)
    npe_default_arg(best_choice_sampling, bool, true)
    npe_default_arg(random_seed, unsigned int, 0)
    npe_default_arg(sample_num_tolerance, float, 0.04)
    npe_doc(downsample_point_cloud_poisson_disk_doc)
    npe_begin_code()
    {
        VCGMesh m;
        vcg_mesh_from_v(v, m);

        if (num_samples <= 0 && radius <= 0.0) {
            throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
        }
        if (sample_num_tolerance > 1.0 || sample_num_tolerance <= 0.0) {
            throw pybind11::value_error("sample_num_tolerance must be in (0, 1]");
        }

        typedef EigenVertexIndexSampler<VCGMesh> PoissonDiskSampler;
        typedef PoissonDiskSampler::IndexArray EigenRetI;

        EigenRetI ret_i;
        PoissonDiskSampler mps(m, ret_i);

        typename tri::SurfaceSampling<VCGMesh, PoissonDiskSampler>::PoissonDiskParam pp;
        typename tri::SurfaceSampling<VCGMesh, PoissonDiskSampler>::PoissonDiskParam::Stat pds;
        pp.pds = pds;
        pp.bestSampleChoiceFlag = best_choice_sampling;
        pp.geodesicDistanceFlag = false;
        pp.randomSeed = random_seed;

        if(random_seed) {
            tri::SurfaceSampling<VCGMesh, PoissonDiskSampler>::SamplingRandomGenerator().initialize(random_seed);
        }

        if (radius <= 0.0 && num_samples > 0) {
            num_samples = std::min(num_samples, (int)v.rows());
            tri::SurfaceSampling<VCGMesh, PoissonDiskSampler>::PoissonDiskPruningByNumber(mps, m, num_samples, radius, pp, sample_num_tolerance);
        } else if (radius > 0.0 && num_samples <= 0) {
            tri::SurfaceSampling<VCGMesh, PoissonDiskSampler>::PoissonDiskPruning(mps, m, radius, pp);
        }
        mps.trim();
        return npe::move(ret_i);
    }
npe_end_code()


npe_function(downsample_point_cloud_voxel_grid_internal)
    npe_arg(v, dense_float, dense_double)
    npe_arg(attrib, dense_float, dense_double)
    npe_arg(voxel_size_x, double)
    npe_arg(voxel_size_y, double)
    npe_arg(voxel_size_z, double)
    npe_arg(voxel_min_x, double)
    npe_arg(voxel_min_y, double)
    npe_arg(voxel_min_z, double)
    npe_arg(voxel_max_x, double)
    npe_arg(voxel_max_y, double)
    npe_arg(voxel_max_z, double)
    npe_arg(min_points_per_voxel, int)
    npe_begin_code()
    {
        EigenDense<npe_Scalar_v> out_v;
        EigenDense<npe_Scalar_attrib> out_attrib;
        Eigen::Matrix<npe_Scalar_v, 3, 1> voxel_size((npe_Scalar_v) voxel_size_x,
                                                     (npe_Scalar_v) voxel_size_y,
                                                     (npe_Scalar_v) voxel_size_z);
        Eigen::Matrix<npe_Scalar_v, 3, 1> voxel_min_bound((npe_Scalar_v) voxel_min_x,
                                                          (npe_Scalar_v) voxel_min_y,
                                                          (npe_Scalar_v) voxel_min_z);
        Eigen::Matrix<npe_Scalar_v, 3, 1> voxel_max_bound((npe_Scalar_v) voxel_max_x,
                                                          (npe_Scalar_v) voxel_max_y,
                                                          (npe_Scalar_v) voxel_max_z);

        downsample_point_cloud_to_voxels(v, attrib, voxel_size, voxel_min_bound, voxel_max_bound,
                                         out_v, out_attrib, min_points_per_voxel);

        return std::make_tuple(npe::move(out_v), npe::move(out_attrib));
    }
npe_end_code()