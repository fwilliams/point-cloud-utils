#include <Eigen/Core>

#ifndef COMMON_H
#define COMMON_H

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

#endif // COMMON_H
