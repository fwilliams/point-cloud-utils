#include <pybind11/pybind11.h>
#include <Eigen/Core>

#if defined(_OPENMP)
#include <thread>
#include <omp.h>
#endif

#ifndef COMMON_H
#define COMMON_H

const int IglDefaultOptions = Eigen::RowMajor;

//constexpr int extract_options(int options) {
//    if (options == (Eigen::ColMajor | Eigen::DontAlign) || options == (Eigen::RowMajor | Eigen::DontAlign)) {
//        return Eigen::RowMajor;
//    } else {
//        return options;
//    }
//}


constexpr bool opts_dontalign(int options) {
    return (options == (Eigen::ColMajor | Eigen::DontAlign) || options == (Eigen::RowMajor | Eigen::DontAlign));
}

template <bool DontAlign, int Opts>
struct OptExtractor {
    enum {
        Options = Opts
    };
};

template <int Opts>
struct OptExtractor<true, Opts> {
    enum Options {
        Options = Eigen::RowMajor
    };
};

template <typename LikeT>
using EigenSparseLike = Eigen::SparseMatrix<typename LikeT::Scalar, Eigen::ColMajor>; // FIXME: Maybe we should output CSR if LikeT is row major

template <typename LikeT>
using EigenDenseLike = Eigen::Matrix<typename LikeT::Scalar, Eigen::Dynamic, Eigen::Dynamic, OptExtractor<opts_dontalign(LikeT::Options), LikeT::Options>::Options, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using EigenDense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic>;

typedef Eigen::Matrix<std::float_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseF32;
typedef Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseF64;
typedef Eigen::Matrix<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseI32;
typedef Eigen::Matrix<std::int64_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseI64;


template <typename TV>
//void validate_mesh(const Eigen::MatrixBase<TV>& v, const Eigen::MatrixBase<TF>& f) {
void validate_point_cloud(const TV& v, bool allow_0=true) {
    if (!allow_0) {
        if (v.rows() == 0) {
            std::stringstream ss;
            ss << "Invalid input point cloud with zero points: points must have shape (n, 3) (n > 0). "
               << "Got points.shape =(" << v.rows() << ", " << v.cols() << ").";
            throw pybind11::value_error(ss.str());
        }
    }

    if (v.cols() != 3) {
        std::stringstream ss;
        ss << "Only 3D inputs are supported: v must have shape (n, 3) (n > 0). "
           << "Got points.shape =(" << v.rows() << ", " << v.cols() << ").";
        throw pybind11::value_error(ss.str());
    }
}

template <typename TV, typename TN>
//void validate_mesh(const Eigen::MatrixBase<TV>& v, const Eigen::MatrixBase<TF>& f) {
void validate_point_cloud_normals(const TV& v, const TN& n, bool allow_0=true) {
    if (!allow_0) {
        if (v.rows() == 0) {
            std::stringstream ss;
            ss << "Invalid input point cloud with zero points: points must have shape (n, 3) (n > 0). "
               << "Got points.shape =(" << v.rows() << ", " << v.cols() << ").";
            throw pybind11::value_error(ss.str());
        }
    }

    if (v.cols() != 3) {
        std::stringstream ss;
        ss << "Only 3D inputs are supported: v must have shape (n, 3) (n > 0). "
           << "Got points.shape =(" << v.rows() << ", " << v.cols() << ").";
        throw pybind11::value_error(ss.str());
    }

    if (n.cols() != 3) {
        std::stringstream ss;
        ss << "Invalid shape for normals: must have shape (n, 3) (n > 0). "
           << "Got normals.shape =(" << n.rows() << ", " << n.cols() << ").";
        throw pybind11::value_error(ss.str());
    }

    if (n.rows() != v.rows()) {
        std::stringstream ss;
        ss << "Invalid input point cloud. Number of normals must match number of points. "
           << "Got points.shape =(" << v.rows() << ", " << v.cols() << ") and normals.shape = "
           << n.rows() << ", " << n.cols();
        throw pybind11::value_error(ss.str());
    }
}


template <typename TV, typename TF>
//void validate_mesh(const Eigen::MatrixBase<TV>& v, const Eigen::MatrixBase<TF>& f) {
void validate_mesh(const TV& v, const TF& f) {
  if (v.rows() == 0 || f.rows() == 0) {
    std::stringstream ss;
    ss << "Invalid input mesh with zero elements: v and f must have shape (n, 3) and (m, 3) (n, m > 0). Got v.shape =("
       << v.rows() << ", " << v.cols() << "), f.shape = ("  << f.rows() << ", " << f.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  if (v.cols() != 3 || f.cols() != 3) {
    std::stringstream ss;
    ss << "Only 3D inputs are supported: v and f must have shape (n, 3) and (m, 3) (n, m > 0). Got v.shape =("
       << v.rows() << ", " << v.cols() << "), f.shape = ("  << f.rows() << ", " << f.cols() << ").";
    throw pybind11::value_error(ss.str());
  }
}


template <typename TV, typename TF, typename TN>
//void validate_mesh(const Eigen::MatrixBase<TV>& v, const Eigen::MatrixBase<TF>& f) {
void validate_mesh(const TV& v, const TF& f, const TN& n) {
  if (v.rows() == 0 || f.rows() == 0) {
    std::stringstream ss;
    ss << "Invalid input mesh with zero elements: v and f must have shape (n, 3) and (m, 3) (n, m > 0). Got v.shape =("
       << v.rows() << ", " << v.cols() << "), f.shape = ("  << f.rows() << ", " << f.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  if (v.cols() != 3 || f.cols() != 3) {
    std::stringstream ss;
    ss << "Only 3D inputs are supported: v and f must have shape (n, 3) and (m, 3) (n, m > 0). Got v.shape =("
       << v.rows() << ", " << v.cols() << "), f.shape = ("  << f.rows() << ", " << f.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  if (n.rows() > 0 && n.cols() != 3) {
      std::stringstream ss;
      ss << "Invalid normals must have shape (#v, 3). Got n.shape = ("
         << n.rows() << ", " << n.cols() << ").";
      throw pybind11::value_error(ss.str());
  }
  if (n.rows() > 0 && v.rows() != n.rows()) {
      std::stringstream ss;
      ss << "Invalid normals must have the same shape as vertices. Got v.shape = ("
         << v.rows() << ", " << v.cols() << "), n.shape = ("  << n.rows() << ", " << n.cols() << ").";
      throw pybind11::value_error(ss.str());
  }
}


class OmpSetParallelism {
private:
    int old_num_threads;
    bool is_nop;
public:
    OmpSetParallelism(int num_threads, bool use_parallel=true) {
        #if defined(_OPENMP)
        is_nop = !use_parallel;
        if (is_nop) {
            return;
        }
        old_num_threads = omp_get_num_threads();
        if (num_threads < 0) {  // Use all processors if you pass in -1 for number of threads
            const auto processor_count = std::thread::hardware_concurrency();
            omp_set_num_threads(processor_count);
            // std::cout << "PROCESSOR COUNT " << processor_count << std::endl;
        } else { // Otherwise, set the number of threads explicitly
            omp_set_num_threads(num_threads);
        }
        #endif
    }

    ~OmpSetParallelism() {
        #if defined(_OPENMP)
        if (is_nop) {
            return;
        }
        omp_set_num_threads(old_num_threads);
        #endif
    }
};

#endif // COMMON_H
