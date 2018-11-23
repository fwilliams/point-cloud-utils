#include <npe.h>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>

#include "nanoflann.hpp"


template <typename Source, typename Target, typename DerivedCorrs, typename DerivedDists>
void compute_distance(
    Source& query_mat,
    Target& dataset_mat,
    Eigen::PlainObjectBase<DerivedCorrs> &corrs,
    Eigen::PlainObjectBase<DerivedDists> &distances)
{
  assert(query_mat.cols() == 3);
  assert(dataset_mat.cols() == 3);

  using namespace nanoflann;
  typedef KDTreeEigenMatrixAdaptor<Source, 3, nanoflann::metric_L2> kd_tree;

  kd_tree mat_index(3, std::cref(dataset_mat), 10 /* max leaf */);
  mat_index.index->buildIndex();

  std::array<typename Source::Scalar, 3> query_point;
  std::array<long int, 1> out_indices;
  std::array<typename Source::Scalar, 1> out_dists_sqr;

  corrs.resize(query_mat.rows(), 1);
  distances.resize(query_mat.rows(), 1);

  for(int i = 0; i < query_mat.rows(); ++i)
  {
    query_point[0] = query_mat(i, 0);
    query_point[1] = query_mat(i, 1);
    query_point[2] = query_mat(i, 2);

    size_t founds = mat_index.index->knnSearch(&query_point[0], 1, &out_indices[0], &out_dists_sqr[0]);
    assert(founds == 1);

    corrs(i) = out_indices[0];
    distances(i) = out_dists_sqr[0];
  }
}

const char* point_cloud_distance_ds = R"Qu8mg5v7(
Compute the shortest distances from each point in the source point cloud to the target point cloud.

Parameters
----------
source : n by d array of representing a set of n points (each row is a point of dimension d)
target : m by d array of representing a set of m points (each row is a point of dimension d)

Returns
-------
A pair `(dists, corrs)` where dists and corrs have shape (n,). `dists[i]` contains the shortest
distance from the point `source[i, :]` to `target`. `corrs[i]` contains the index into `target`
of the nearest point to `source[i, :]`.
)Qu8mg5v7";


npe_function(point_cloud_distance)
npe_arg(source, dense_f32, dense_f64)
npe_arg(target, npe_matches(source))
npe_doc(point_cloud_distance_ds)
npe_begin_code()

npe_Matrix_source src = source;
npe_Matrix_target dst = target;
npe_Matrix_source dists;
Eigen::Matrix<long int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> corrs;

compute_distance(src, dst, corrs, dists);

return std::make_tuple(npe::move(dists), npe::move(corrs));

npe_end_code()
