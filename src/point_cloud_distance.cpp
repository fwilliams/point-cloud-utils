#include <npe.h>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "nanoflann.hpp"
#include "common.h"


namespace {

/*
 * Compute shortest distance from one point cloud to another using nanoflann
 */
template <typename Source, typename DerivedCorrs, typename DerivedDists>
void shortest_distances_nanoflann(Source& query_mat,
                                  Source& dataset_mat,
                                  Eigen::PlainObjectBase<DerivedCorrs> &corrs,
                                  Eigen::PlainObjectBase<DerivedDists> &distances,
                                  int num_nbrs=1, bool squared_dist=false,
                                  int max_points_per_leaf=10,
                                  int num_threads=0) {
    const int MIN_PARALLEL_INPUT_SIZE = 100000;
    const bool run_parallel = query_mat.rows() >= MIN_PARALLEL_INPUT_SIZE && num_threads != 0;
    auto set_parallel = OmpSetParallelism(num_threads, run_parallel);

    assert(query_mat.cols() == 3);
    assert(dataset_mat.cols() == 3);

    using namespace nanoflann;
    using KdTreeType = nanoflann::KDTreeEigenMatrixAdaptor<Source, 3, nanoflann::metric_L2_Simple>;
    using IndexType = typename KdTreeType::IndexType;
    using ScalarType = typename KdTreeType::num_t;

    KdTreeType mat_index(3, std::cref(dataset_mat), max_points_per_leaf /* max leaf */);
    mat_index.index->buildIndex();

    corrs.resize(query_mat.rows(), num_nbrs);
    distances.resize(query_mat.rows(), num_nbrs);

    bool threw_exception = false;

    #if defined(_OPENMP)
    #pragma omp parallel if (run_parallel)
    #endif
    {
        std::array<ScalarType, 3> query_point;
        std::vector<IndexType> out_indices(num_nbrs);
        std::vector<ScalarType> out_dists_sqr(num_nbrs);
        #if defined(_OPENMP)
        #pragma omp for
        #endif
        for(int i = 0; i < query_mat.rows(); ++i) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                    // This doesn't work on Windows :(
                    //#pragma omp cancel for
                #else
                    threw_exception = true;
                    break;
                #endif
            }
            for (int j = 0; j < query_mat.cols(); ++j) { query_point[j] = query_mat(i, j); }

            const size_t founds = mat_index.index->knnSearch(query_point.data(), num_nbrs,
                                                             out_indices.data(), out_dists_sqr.data());
            assert(founds >= 1);

            for (int k = 0; k < founds; k++) {
                corrs(i, k) = out_indices[k];
                if (squared_dist) {
                    distances(i, k) = out_dists_sqr[k];
                } else {
                    distances(i, k) = sqrt(out_dists_sqr[k]);
                }
            }
            for (int k = founds; k < num_nbrs; k++) {
                corrs(i, k) = -1;
                distances(i, k) = -1.0;
            }
        }
    }
    if (threw_exception) {
        throw pybind11::error_already_set();
    }
}

} // namespace




const char* k_nearest_neighbors_doc = R"Qu8mg5v7(
Compute the k nearest neighbors (L2 distance) from each point in the query point cloud to the dataset point cloud.

Parameters
----------
query_points : n by 3 array of representing a set of n points (each row is a point of dimension 3).
dataset_points : m by 3 array of representing a set of m points (each row is a point of dimension 3).
k      : the number of nearest neighbors to query per point.
squared_distances : If set to True, then return squared L2 distances. Default is False.
max_points_per_leaf : The maximum number of points per leaf node in the KD tree used by this function.
                      Default is 10.
num_threads : Number of threads to use. If set to -1, will use all available CPUs. If set to 0, will run in serial. Default is -1

Returns
-------
A pair `(dists, corrs)` where dists and corrs have shape (n, k). `dists[i, k]` contains the k^th shortest
 L2 distance from the point `query_points[i, :]` to `dataset_points`. `corrs[i, k]` contains the index into
`dataset_points` of the k^th nearest point to `query_points[i, :]`.
)Qu8mg5v7";

npe_function(k_nearest_neighbors)
npe_arg(query_points, dense_float, dense_double)
npe_arg(dataset_points, npe_matches(query_points))
npe_arg(k, int)
npe_default_arg(squared_distances, bool, false)
npe_default_arg(max_points_per_leaf, int, 10)
npe_default_arg(num_threads, int, -1)
npe_doc(k_nearest_neighbors_doc)
npe_begin_code()
{
    if (k <= 0) {
        throw pybind11::value_error("Invalid value for k (" + std::to_string(k) + ") must be greater than 0.");
    }
    if (query_points.rows() == 0 || dataset_points.rows() == 0) {
        std::stringstream ss;
        ss << "Invalid input set with zero elements: query_points and dataset_points must have shape (n, 3) and (m, 3). Got query_points.shape = ("
           << query_points.rows() << ", " << query_points.cols() << "), " << "dataset_points.shape = (" << dataset_points.rows() << ", " << dataset_points.cols() << ").";
        throw pybind11::value_error(ss.str());
    }

    if (query_points.cols() != 3 || dataset_points.cols() != 3) {
        std::stringstream ss;
        ss << "Only 3D inputs are supported: query_points and dataset_points must have shape (n, 3) and (m, 3). Got query_points.shape = ("
           << query_points.rows() << ", " << query_points.cols() << "), "
           << "dataset_points.shape = (" << dataset_points.rows() << ", " << dataset_points.cols() << ").";
        throw pybind11::value_error(ss.str());
    }

    // FIXME: nanoflann does not work with Eigen::Maps so we have to do a copy here :(
    EigenDenseLike<npe_Matrix_query_points> src = query_points;
    EigenDenseLike<npe_Matrix_dataset_points> dst = dataset_points;
    EigenDenseLike<npe_Matrix_query_points> dists;

    //  using kd_tree = nanoflann::KDTreeEigenMatrixAdaptor<EigenDenseLike<npe_Matrix_query_points>, 3, nanoflann::metric_L2>;
    using IndexType = typename npe_Matrix_query_points::Index;
    Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> corrs;

    shortest_distances_nanoflann(src, dst, corrs, dists, k, squared_distances, max_points_per_leaf, num_threads);

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
squared_distances : If set to True, then return squared L2 distances.
max_points_per_leaf : the maximum number of points per leaf node in the KD tree used by this function. Default is 10.

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
npe_default_arg(squared_distances, bool, false)
npe_default_arg(max_points_per_leaf, int, 10)
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

    shortest_distances_nanoflann(src, dst, corrs, dists, 1, squared_distances, max_points_per_leaf);

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
