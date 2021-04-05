#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

#include <fstream>
#include <iostream>
#include <functional>

#include "common.h"
#include "vcg_utils.h"


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


class AccumulatedPoint {
public:
    AccumulatedPoint()
            : num_of_points_(0),
              point_(0.0, 0.0, 0.0),
              normal_(0.0, 0.0, 0.0),
              color_(0.0, 0.0, 0.0),
              alpha_(0.0) {}

public:
    void AddPoint(Eigen::Vector3d point, Eigen::Vector3d normal, Eigen::Vector3d color, double alpha) {
        point_ += point;
        normal_ += normal;
        color_ += color;
        alpha_ += alpha;
        num_of_points_ += 1;
    }

    Eigen::Vector3d GetAveragePoint() const {
        return point_ / double(num_of_points_);
    }

    Eigen::Vector3d GetAverageNormal() const {
        // Call NormalizeNormals() afterwards if necessary
        return normal_ / double(num_of_points_);
    }

    Eigen::Vector3d GetAverageColor() const {
        return color_ / double(num_of_points_);
    }

    double GetAverageAlpha() const {
        return alpha_ / double(num_of_points_);
    }

    int GetNumPoints() const {
        return num_of_points_;
    }

public:
    int num_of_points_;
    Eigen::Vector3d point_;
    Eigen::Vector3d normal_;
    Eigen::Vector3d color_;
    double alpha_;
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


template <typename DerivedV, typename DerivedN, typename DerivedC,
        typename DerivedOutV, typename DerivedOutN, typename DerivedOutC>
void downsample_point_cloud_to_voxels(const DerivedV& V,
                                      const DerivedN& N,
                                      const DerivedC& C,
                                      Eigen::Vector3d voxel_size3,
                                      Eigen::Vector3d voxel_min_bound,
                                      Eigen::Vector3d voxel_max_bound,
                                      DerivedOutV& outV,
                                      DerivedOutN& outN,
                                      DerivedOutC& outC,
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

    bool has_normals = (N.rows() != 0 && N.cols() != 0);
    bool has_colors = (C.rows() != 0 && C.cols() != 0);
    if (has_normals && N.rows() != V.rows()) {
        throw pybind11::value_error("Invalid number of normals (" + std::to_string(N.rows()) +
                                    "). Must match number of input vertices (" +
                                    std::to_string(V.rows()) + ") or b 0.");
    }
    if (has_normals && N.cols() != 3) {
        throw pybind11::value_error("Invalid dimension for normals. Got shape (" + std::to_string(N.rows()) +
                                    ", " + std::to_string(N.cols()) + "). Must have shape (N, 3)");
    }
    if (has_colors && C.rows() != V.rows()) {
        throw pybind11::value_error("Invalid number of colors (" + std::to_string(C.rows()) +
                                    "). Must match number of input vertices (" +
                                    std::to_string(V.rows()) + ") or b 0.");
    }
    if (has_colors && (C.cols() != 3 && C.cols() != 4)) {
        throw pybind11::value_error("Invalid dimension for colors. Got shape (" + std::to_string(C.rows()) +
                                    ", " + std::to_string(C.cols()) + "). Must have shape (N, 3) or (N, 4)");
    }
    bool color_has_alpha_channel = C.cols() == 4;

    std::unordered_map<Eigen::Vector3i, AccumulatedPoint, hash_eigen<Eigen::Vector3i>> voxelindex_to_accpoint;

    Eigen::Vector3d ref_coord;
    Eigen::Vector3d ref_point;
    Eigen::Vector3d ref_normal;
    Eigen::Vector3d ref_color;
    double ref_alpha;
    Eigen::Vector3i voxel_index;
    for (int i = 0; i < V.rows(); i++) {
        ref_point = Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2));
        ref_coord = (ref_point - voxel_min_bound).array() / voxel_size3.array();
        voxel_index << int(floor(ref_coord(0))),
                int(floor(ref_coord(1))),
                int(floor(ref_coord(2)));
        if (has_normals) {
            ref_normal = Eigen::Vector3d(N(i, 0), N(i, 1), N(i, 2));
        }
        if (has_colors) {
            ref_color = Eigen::Vector3d(C(i, 0), C(i, 1), C(i, 2));
            ref_alpha = color_has_alpha_channel ? C(i, 3) : 1.0;
        }
        voxelindex_to_accpoint[voxel_index].AddPoint(ref_point, ref_normal, ref_color, ref_alpha);
    }

    int num_output_points = voxelindex_to_accpoint.size();
    outV.resize(num_output_points, 3);
    if (has_normals) {
        outN.resize(num_output_points, 3);
    }
    if (has_colors) {
        int dim = color_has_alpha_channel ? 4 : 3;
        outC.resize(num_output_points, dim);
    }

    int count = 0;
    for (auto accpoint : voxelindex_to_accpoint) {
        if (accpoint.second.GetNumPoints() < min_pts_per_bin) {
            continue;
        }
        Eigen::RowVector3d accpoint_pos = accpoint.second.GetAveragePoint();
        for (int i = 0; i < 3; i++) {outV(count, i) = accpoint_pos[i];}

        if (has_normals) {
            Eigen::RowVector3d accpoint_nml = accpoint.second.GetAverageNormal();
            for (int i = 0; i < 3; i++) {outN(count, i) = accpoint_nml[i];}
        }
        if (has_colors) {
            Eigen::Vector3d accpoint_rgb = accpoint.second.GetAverageColor();
            for (int i = 0; i < 3; i++) {outC(count, i) = accpoint_rgb[i];}
            if (color_has_alpha_channel) {
                outC(count, 3) = accpoint.second.GetAverageAlpha();
            }
        }
        count += 1;
    }
    if (count < num_output_points) {
        outV.conservativeResize(count, 3);
        if (has_normals) {
            outN.conservativeResize(count, 3);
        }
        if (has_colors) {
            int dim = color_has_alpha_channel ? 4 : 3;
            outC.conservativeResize(count, dim);
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
    npe_arg(n, npe_matches(v))
    npe_arg(c, npe_matches(v))
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
        EigenDense<npe_Scalar_v> outV;
        EigenDense<npe_Scalar_n> outN;
        EigenDense<npe_Scalar_c> outC;
        Eigen::Vector3d voxel_size(voxel_size_x, voxel_size_y, voxel_size_z);
        Eigen::Vector3d voxel_min_bound(voxel_min_x, voxel_min_y, voxel_min_z);
        Eigen::Vector3d voxel_max_bound(voxel_max_x, voxel_max_y, voxel_max_z);

        downsample_point_cloud_to_voxels(v, n, c,
                                         voxel_size, voxel_min_bound, voxel_max_bound,
                                         outV, outN, outC, min_points_per_voxel);

        return std::make_tuple(npe::move(outV), npe::move(outN), npe::move(outC));
    }
npe_end_code()