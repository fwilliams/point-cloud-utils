#include <npe.h>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>

#include "nanoflann.hpp"
#include "common.h"


namespace {

/*
 * Compute shortest distance from one point cloud to another using nanoflann
 */
template <typename Source, typename Target, typename DerivedCorrs, typename DerivedDists>
void shortest_distances_nanoflann(
    Source& query_mat,
    Target& dataset_mat,
    Eigen::PlainObjectBase<DerivedCorrs> &corrs,
    Eigen::PlainObjectBase<DerivedDists> &distances,
    int num_nbrs=1,
    bool squared_dist=false) {
  assert(query_mat.cols() == 3);
  assert(dataset_mat.cols() == 3);

  using namespace nanoflann;
  using kd_tree = nanoflann::KDTreeEigenMatrixAdaptor<Source, 3, nanoflann::metric_L2_Simple>;
  using IndexType = typename kd_tree::IndexType;

  kd_tree mat_index(3, std::cref(dataset_mat), 10 /* max leaf */);
  mat_index.index->buildIndex();

  std::array<typename Source::Scalar, 3> query_point;
  std::vector<IndexType> out_indices(num_nbrs);
  std::vector<typename Source::Scalar> out_dists_sqr(num_nbrs);

  corrs.resize(query_mat.rows(), num_nbrs);
  distances.resize(query_mat.rows(), num_nbrs);

  for(int i = 0; i < query_mat.rows(); ++i) {
    for (int j = 0; j < query_mat.cols(); ++j) { query_point[j] = query_mat(i, j); }

    size_t founds = mat_index.index->knnSearch(query_point.data(), num_nbrs, &out_indices[0], &out_dists_sqr[0]);
    assert(founds >= 1);

    for (int k = 0; k < founds; k++) {
        corrs(i, k) = out_indices[k];
        if (squared_dist) {
            distances(i, k) = out_dists_sqr[0];
        } else {
            distances(i, k) = sqrt(out_dists_sqr[0]);
        }
    }
    for (int k = founds; k < num_nbrs; k++) {
        corrs(i, k) = -1;
        distances(i, k) = -1.0;
    }
  }
}

} // namespace




const char* k_nearest_neighbors_doc = R"Qu8mg5v7(
Compute the k nearest neighbors (L2 distance) from each point in the source point cloud to the target point cloud.

Parameters
----------
source : n by 3 array of representing a set of n points (each row is a point of dimension 3)
target : m by 3 array of representing a set of m points (each row is a point of dimension 3)

Returns
-------
A pair `(dists, corrs)` where dists and corrs have shape (n, k). `dists[i, k]` contains the k^th shortest
 L2 distance from the point `source[i, :]` to `target`. `corrs[i, k]` contains the index into
`target` of the k^th nearest point to `source[i, :]`.
)Qu8mg5v7";

npe_function(k_nearest_neighbors)
npe_arg(source, dense_float, dense_double)
npe_arg(target, npe_matches(source))
npe_arg(k, int)
npe_doc(k_nearest_neighbors_doc)
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

  shortest_distances_nanoflann(src, dst, corrs, dists, k);

  return std::make_tuple(npe::move(dists), npe::move(corrs));
}
npe_end_code()




const char* one_sided_hausdorff_distance_doc = R"Qu8mg5v7(
Compute the one sided Hausdorff distance from source to target

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

npe_function(one_sided_hausdorff_distance)
npe_arg(source, dense_float, dense_double)
npe_arg(target, npe_matches(source))
npe_default_arg(return_index, bool, false)
npe_doc(one_sided_hausdorff_distance_doc)
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

  shortest_distances_nanoflann(src, dst, corrs, dists);

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
