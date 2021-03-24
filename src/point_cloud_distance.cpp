#include <npe.h>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>

#include "nanoflann.hpp"
#include "common.h"


namespace {

/*
 * Compute shortest distance from one point cloud to another using nanoflann
 */
template <typename Source, typename Target, typename DerivedCorrs, typename DerivedDists>
void compute_distance(
    Source& query_mat,
    Target& dataset_mat,
    Eigen::PlainObjectBase<DerivedCorrs> &corrs,
    Eigen::PlainObjectBase<DerivedDists> &distances) {
  assert(query_mat.cols() == 3);
  assert(dataset_mat.cols() == 3);

  using namespace nanoflann;
  using kd_tree = nanoflann::KDTreeEigenMatrixAdaptor<Source, 3, nanoflann::metric_L2>;
  using IndexType = typename kd_tree::IndexType;

  kd_tree mat_index(3, std::cref(dataset_mat), 10 /* max leaf */);
  mat_index.index->buildIndex();

  std::array<typename Source::Scalar, 3> query_point;
  std::array<IndexType, 1> out_indices;
  std::array<typename Source::Scalar, 1> out_dists_sqr;

  corrs.resize(query_mat.rows(), 1);
  distances.resize(query_mat.rows(), 1);

  for(int i = 0; i < query_mat.rows(); ++i) {
    for (int j = 0; j < query_mat.cols(); ++j) { query_point[j] = query_mat(i, j); }

    size_t founds = mat_index.index->knnSearch(query_point.data(), 1, &out_indices[0], &out_dists_sqr[0]);
    assert(founds == 1);

    corrs(i) = out_indices[0];
    distances(i) = out_dists_sqr[0];
  }
}

} // namespace









const char* shortest_squared_distances_doc = R"Qu8mg5v7(
Compute the shortest squared L2 distances from each point in the source point cloud to the target point cloud.

Parameters
----------
source : n by 3 array of representing a set of n points (each row is a point of dimension 3)
target : m by 3 array of representing a set of m points (each row is a point of dimension 3)

Returns
-------
A pair `(dists, corrs)` where dists and corrs have shape (n,). `dists[i]` contains the shortest
squared L2 distance from the point `source[i, :]` to `target`. `corrs[i]` contains the index into
`target` of the nearest point to `source[i, :]`.
)Qu8mg5v7";

npe_function(shortest_squared_distances)
npe_arg(source, dense_float, dense_double)
npe_arg(target, npe_matches(source))
npe_doc(shortest_squared_distances_doc)
npe_begin_code()
{
  if (source.rows() == 0 || target.rows() == 0) {
    std::stringstream ss;
    ss << "Invalid input set with zero elements: source and targets must have shape (n, 3) and (m, 3). Got source.shape = ("
       << source.rows() << ", " << source.cols() << "), " << "target.shape = (" << target.rows() << ", " << target.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  if (source.cols() != 3 || target.cols() != 3) {
    std::stringstream ss;
    ss << "Only 3D inputs are supported: source and targets must have shape (n, 3) and (m, 3). Got source.shape = ("
       << source.rows() << ", " << source.cols() << "), "
       << "target.shape = (" << target.rows() << ", " << target.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  // FIXME: nanoflann does not work with Eigen::Maps so we have to do a copy here :(
  EigenDenseLike<npe_Matrix_source> src = source;
  EigenDenseLike<npe_Matrix_target> dst = target;
  EigenDenseLike<npe_Matrix_source> dists;

  using kd_tree = nanoflann::KDTreeEigenMatrixAdaptor<EigenDenseLike<npe_Matrix_source>, 3, nanoflann::metric_L2>;
  using IndexType = typename kd_tree::IndexType;
  Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> corrs;

  compute_distance(src, dst, corrs, dists);

  return std::make_tuple(npe::move(dists), npe::move(corrs));
}
npe_end_code()




const char* squared_hausdorff_distance_doc = R"Qu8mg5v7(
Compute the one sided squared Hausdorff distance from source to target

Parameters
----------
source : n by 3 array of representing a set of n points (each row is a point of dimension 3)
target : m by 3 array of representing a set of m points (each row is a point of dimension 3)
return_index : Optionally return the index pair `(i, j)` into source and target such that
               `source[i, :]` and `target[j, :]` are the two points with maximum shortest distance.
Returns
-------
The largest shortest distance, `d` between each point in `source` and the points in `target`.
If `return_index` is set, then this function returns a tuple (d, i, j) where `d` is as described above
and `(i, j)` are such that `source[i, :]` and `target[j, :]` are the two points with maximum shortest
distance.

)Qu8mg5v7";

npe_function(squared_hausdorff_distance)
npe_arg(source, dense_float, dense_double)
npe_arg(target, npe_matches(source))
npe_default_arg(return_index, bool, false)
npe_doc(squared_hausdorff_distance_doc)
npe_begin_code()
{
  if (source.rows() == 0 || target.rows() == 0) {
    std::stringstream ss;
    ss << "Invalid input set with zero elements: source and targets must have shape (n, 3) and (m, 3). Got source.shape = ("
       << source.rows() << ", " << source.cols() << "), " << "target.shape = (" << target.rows() << ", " << target.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  if (source.cols() != 3 || target.cols() != 3) {
    std::stringstream ss;
    ss << "Only 3D inputs are supported: source and targets must have shape (n, 3) and (m, 3). Got source.shape = ("
       << source.rows() << ", " << source.cols() << "), "
       << "target.shape = (" << target.rows() << ", " << target.cols() << ").";
    throw pybind11::value_error(ss.str());
  }

  // FIXME: nanoflann does not work with Eigen::Maps so we have to do a copy here :(
  EigenDenseLike<npe_Matrix_source> src = source;
  EigenDenseLike<npe_Matrix_target> dst = target;
  EigenDenseLike<npe_Matrix_source> dists;

  using kd_tree = nanoflann::KDTreeEigenMatrixAdaptor<EigenDenseLike<npe_Matrix_source>, 3, nanoflann::metric_L2>;
  using IndexType = typename kd_tree::IndexType;
  Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> corrs;

  compute_distance(src, dst, corrs, dists);

  size_t max_index_source = -1;
  size_t dummy = -1;
  npe_Scalar_source max_dist = dists.maxCoeff(&max_index_source, &dummy);
  assert(dummy == 0);
  size_t max_index_target = corrs(max_index_source, 0);

  if (return_index) {
    auto ret = std::make_tuple(max_dist, max_index_source, max_index_target);
    return pybind11::cast(ret);
  } else {
    return pybind11::cast(max_dist);
  }
}
npe_end_code()
