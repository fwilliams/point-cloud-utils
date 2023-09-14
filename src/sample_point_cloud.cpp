#include <npe.h>

#include <fstream>
#include <iostream>
#include <functional>
#include <algorithm>

#include <igl/blue_noise.h>

#include "common/common.h"

namespace {

template <typename DerivedV, typename DerivedXI, typename URBG>
IGL_INLINE void blue_noise_downsample(const Eigen::MatrixBase<DerivedV> & X,
                                      const typename DerivedV::Scalar r,
                                      Eigen::PlainObjectBase<DerivedXI> & XI,
                                      URBG && urbg = igl::generate_default_urbg()) {
    typedef int64_t BlueNoiseKeyType;
    typedef typename DerivedV::Scalar Scalar;
    // float+RowMajor is faster...
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,3,Eigen::RowMajor> MatrixX3S;
    assert(X.cols() == 3 && "Only 3D embeddings allowed");
    // minimum radius
    const Scalar min_r = r;
    // cell size based on 3D distance
    // It works reasonably well (but is probably biased to use s=2*r/√3 here and
    // g=1 in the outer loop below.
    //
    // One thing to try would be to store a list in S (rather than a single point)
    // or equivalently a mask over M and just use M as a generic spatial hash
    // (with arbitrary size) and then tune its size (being careful to make g a
    // function of r and s; and removing the `if S=-1 checks`)
    const Scalar s = r/sqrt(3.0);
  
    // Make a uniform random sampling with 30*expected_number_of_points.
    const int nx = X.rows();
  
    // Rescale so that s = 1
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> Xs = ((X.rowwise()-X.colwise().minCoeff())/s).template cast<int>();
    const int w = Xs.maxCoeff() + 1;
    Eigen::VectorXi SortIdx;
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic, DerivedXI::Options> Xsorted;
    {
        igl::sortrows(decltype(Xs)(Xs), true, Xs, SortIdx);
        Xsorted = X(SortIdx, Eigen::all).eval();
    }
    // Initialization
    std::unordered_map<BlueNoiseKeyType,std::vector<int> > M;
    std::unordered_map<BlueNoiseKeyType, int > S;
    // attempted to seed
    std::unordered_map<BlueNoiseKeyType, int > A;
    // Q: Too many?
    // A: Seems to help though.
    M.reserve(Xs.rows());
    S.reserve(Xs.rows());
    for(int i = 0;i<Xs.rows();i++) {
        BlueNoiseKeyType k = igl::blue_noise_key(w,Xs(i,0),Xs(i,1),Xs(i,2));
        const auto Miter = M.find(k);
        if(Miter  == M.end()) {
            M.insert({k,{i}});
        } else {
            Miter->second.push_back(i);
        }
        S.emplace(k,-1);
        A.emplace(k,false);
    }

    std::vector<int> active;
    // precompute r²
    // Q: is this necessary?
    const double rr = r*r;
    std::vector<int> collected;
    collected.reserve(nx);

    auto Mouter = M.begin();
    // Just take the first point as the initial seed
    const auto initialize = [&]()->bool {
        while(true) {
            if(Mouter == M.end()) {
                return false;
            }
            const BlueNoiseKeyType k = Mouter->first;
            // Haven't placed in this cell yet
            if(S[k]<0) {
                if(igl::activate(Xsorted,Xs,rr,-1,w,k,M,S,active)) { return true; }
            }
            Mouter++;
        }
        assert(false && "should not be reachable.");
    };

    // important if mesh contains many connected components
    while(initialize()) {
        while(active.size()>0) {
            igl::step(Xsorted,Xs,rr,w,urbg,M,S,active,collected);
        }
    }

    {
      const int n = collected.size();
      XI.resize(n);
      for(int i = 0;i<n;i++) {
          const int c = collected[i];
          XI(i) = SortIdx[c];
          //   P.row(i) = X.row(c).template cast<typename DerivedP::Scalar>();
      }
    }
}

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
            if (attrib_.cols() == 0) {
                attrib_.resize(1, attrib.cols());
            }
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

Args:
    v: \#v by 3 array of vertex positions
    radius: desired separation between points.
    target_num_samples: If set to a positive value, iterate to generate points as close to this target as possible (determined by sample_num_tolerance). 
    random_seed: A random seed used to generate the samples. Passing in 0 will use the current time. (0 by default).
    sample_num_tolerance: If you requested a target number of samples, by passsing num_samples > 0, then this function will return between (1 - sample_num_tolerance) * num_samples and (1 + sample_num_tolerance) * num_samples. Setting a very small value for this parameter will increase convergence time. (0.04 by default).

Returns:
    p_idx : A (m,) shaped array of indices into v where m is the number of Poisson-disk samples
)Qu8mg5v7";
npe_function(downsample_point_cloud_poisson_disk)
npe_arg(v, dense_float, dense_double)
npe_arg(radius, double)
npe_default_arg(target_num_samples, int, -1)
npe_default_arg(random_seed, unsigned int, 0)
npe_default_arg(sample_num_tolerance, float, 0.04)
npe_doc(downsample_point_cloud_poisson_disk_doc)
npe_begin_code()
{
    // VCGMesh m;
    // vcg_mesh_from_v(v, m);

    if (target_num_samples <= 0 && radius <= 0.0) {
        throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
    }
    if (sample_num_tolerance > 1.0 || sample_num_tolerance <= 0.0) {
        throw pybind11::value_error("sample_num_tolerance must be in (0, 1]");
    }
    if (random_seed != 0) {
        srand(random_seed);
    }


    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> ret_i;

    if (target_num_samples > 0) {
        if (target_num_samples >= v.rows()) {
            ret_i.resize(v.rows(), 1);
            for (int i = 0; i < ret_i.rows(); i += 1) { ret_i(i) = i; }
            return npe::move(ret_i);
        }

        size_t num_samples_min = int(npe_Scalar_v(target_num_samples)*(1.0f-sample_num_tolerance));
        size_t num_samples_max = int(npe_Scalar_v(target_num_samples)*(1.0f+sample_num_tolerance));

        npe_Matrix_v bmin = v.colwise().minCoeff();
        npe_Matrix_v bmax = v.colwise().maxCoeff();
        npe_Scalar_v bbsize = (bmax - bmin).norm();
        npe_Scalar_v range_min_rad = bbsize/50.0;
        npe_Scalar_v range_max_rad = bbsize/50.0;
        size_t range_min_rad_num = -1;
        size_t range_max_rad_num = -1;
        target_num_samples = std::min(target_num_samples, (int) v.rows());
        do {
            ret_i.conservativeResize(0, 1);
            range_min_rad /= 2.0;
            blue_noise_downsample(v, (npe_Scalar_v) range_min_rad, ret_i, igl::generate_default_urbg());
            range_min_rad_num = ret_i.rows();
        } while(range_min_rad_num < target_num_samples);

        do {
            ret_i.conservativeResize(0, 1);
            range_max_rad *= 2.0;
            blue_noise_downsample(v, (npe_Scalar_v) range_max_rad, ret_i, igl::generate_default_urbg());
            range_max_rad_num = ret_i.rows();
        } while(range_max_rad_num > target_num_samples);


        npe_Scalar_v current_rad = range_max_rad;
        int iter_count = 0;
        while (iter_count < 20 && (ret_i.rows() < num_samples_min || ret_i.rows() > num_samples_max)) {
            iter_count += 1;
            ret_i.conservativeResize(0, 1);
            current_rad = (range_min_rad + range_max_rad) / 2.0;
            blue_noise_downsample(v, (npe_Scalar_v) current_rad, ret_i, igl::generate_default_urbg());
            if (ret_i.rows() > target_num_samples) {
                range_min_rad = current_rad;
                range_min_rad_num = ret_i.rows();
            } 
            if (ret_i.rows() < target_num_samples) {
                range_max_rad = current_rad;
                range_max_rad_num = ret_i.rows();
            }
        }
    } else {
        blue_noise_downsample(v, (npe_Scalar_v) radius, ret_i, igl::generate_default_urbg());
    }

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