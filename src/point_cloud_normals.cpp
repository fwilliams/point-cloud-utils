#include <npe.h>
#include <Eigen/SVD>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <thread>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "nanoflann.hpp"
#include "common/common.h"


namespace {

// So annoying that this isn't in the standard library
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename PType, typename DType>
void validate_input(const PType& points, const DType& view_dirs) {
    if (points.rows() == 0 || points.cols() != 3) {
        std::stringstream ss;
        ss << "Invalid point set with zero elements: points must have shape (n, 3), but got ot points.shape = ("
           << points.rows() << ", " << points.cols() << ").";
        throw pybind11::value_error(ss.str());
    }
    if (view_dirs.rows() != 0 && (view_dirs.rows() != points.rows() || view_dirs.cols() != 3)) {
        std::stringstream ss;
        ss << "Invalid view directions does not match the number of points. If view directions are passed in, "
           << "they must have the same shape as points. Got points.shape = ("
           << points.rows() << ", " << points.cols() << "), and view_dirs.shape = ("
           << view_dirs.rows() << ", " << view_dirs.cols() << ").";
        throw pybind11::value_error(ss.str());
    }
}


/*
 * Estimate the normal of a point in a local neighborhood weighted by a radial basis function
 */
template <typename KdTreeType, typename PointType, typename DirType>
std::tuple<bool, Eigen::RowVector3d> estimate_local_normal_rbf(const KdTreeType& tree,
                   const PointType& points,
                   const DirType& view_dirs,
                   int point_index,
                   double ball_radius,
                   int min_pts_per_ball,
                   int max_pts_per_ball,
                   double drop_angle_threshold,
                   const std::function<double(double, double)>& weighting_function)
{
    const bool has_view_directions = view_dirs.rows() > 0;

    using namespace nanoflann;
    using IndexType = typename KdTreeType::IndexType;
    using ScalarType = typename KdTreeType::num_t;

    // If we have view directions, get the view direction for this point
    Eigen::RowVector3d dir;
    if (has_view_directions) {
        dir = Eigen::RowVector3d(view_dirs(point_index, 0), view_dirs(point_index, 1), view_dirs(point_index, 2));
    }
    const std::array<ScalarType, 3> query_point =
        { points(point_index, 0), points(point_index, 1), points(point_index, 2) };

    std::vector<std::pair<IndexType, ScalarType>> out_nbrs;
    size_t founds = tree.index->radiusSearch(query_point.data(), ball_radius, out_nbrs, SearchParams());
    if (founds < min_pts_per_ball) {
        return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
    }

    if (max_pts_per_ball > 0) {
        std::random_shuffle(out_nbrs.begin(), out_nbrs.end());
        founds = std::min((size_t) max_pts_per_ball, founds);
    }

    Eigen::MatrixXd neighbors(founds, 3);
    for (int i = 0; i < founds; i++) {
        const int nbr_idx = out_nbrs[i].first;
        const double nbr_dist = out_nbrs[i].second;
        const double weight = weighting_function(sqrt(nbr_dist), ball_radius);

        for (int j = 0; j < 3; j++) { neighbors(i, j) = (points(nbr_idx, j) - query_point[j]) * weight; }
    }

    // Estimate the normal at this point by fitting a plane to its neighborhood
    Eigen::RowVector3d normal;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(neighbors, Eigen::ComputeThinV);
    Eigen::MatrixXd V = svd.matrixV();
    for (int l = 0; l < 3; l++) {
        normal[l] = V(l, 2);
    }

    if (has_view_directions) {
        // Flip normal to align with view direction
        normal *= sign(normal.dot(dir));
        // Filter out point/normal if the angle with the view direction is too large
        double normal_to_view_angle = acos(normal.dot(dir));
        if (normal_to_view_angle > drop_angle_threshold) {
            return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
        }
    }
    return std::make_tuple(true, normal);

}


template <typename KdTreeType, typename PointType, typename DirType>
std::tuple<bool, Eigen::RowVector3d> estimate_local_normal_knn(const KdTreeType& tree,
                   const PointType& points,
                   const DirType& view_dirs,
                   int point_index,
                   int num_neighbors,
                   double drop_angle_threshold)
{
    const bool has_view_directions = view_dirs.rows() > 0;

    using namespace nanoflann;
    using IndexType = typename KdTreeType::IndexType;
    using ScalarType = typename KdTreeType::num_t;

    // If we have view directions, get the view direction for this point
    Eigen::RowVector3d dir;
    if (has_view_directions) {
        dir = Eigen::RowVector3d(view_dirs(point_index, 0), view_dirs(point_index, 1), view_dirs(point_index, 2));
    }
    const std::array<ScalarType, 3> query_point =
        { points(point_index, 0), points(point_index, 1), points(point_index, 2) };

    std::vector<IndexType> out_indices(num_neighbors);
    std::vector<ScalarType> out_dists_sqr(num_neighbors);
    const size_t founds = tree.index->knnSearch(query_point.data(), num_neighbors,
                                                out_indices.data(), out_dists_sqr.data());

    // If we don't have enough neighbors, then discard this point
    if (founds < num_neighbors) {
        return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
    }

    Eigen::MatrixXd neighbors(founds, 3);
    for (int i = 0; i < founds; i++) {
        const int nbr_idx = out_indices[i];

        for (int j = 0; j < 3; j++) { neighbors(i, j) = points(nbr_idx, j) - query_point[j]; }
    }

    // Estimate the normal at this point by fitting a plane to its neighborhood
    Eigen::RowVector3d normal;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(neighbors, Eigen::ComputeThinV);
    Eigen::MatrixXd V = svd.matrixV();
    for (int l = 0; l < 3; l++) {
        normal[l] = V(l, 2);
    }

    if (has_view_directions) {
        // Flip normal to align with view direction
        normal *= sign(normal.dot(dir));
        // Filter out point/normal if the angle with the view direction is too large
        double normal_to_view_angle = acos(normal.dot(dir));
        if (normal_to_view_angle > drop_angle_threshold) {
            return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
        }
    }

    return std::make_tuple(true, normal);
}

template <typename InPointsType, typename InViewDirsType, typename OutIndexType, typename OutNormalsType>
void estimate_normals_parallel(const InPointsType& points, const InViewDirsType& view_dirs,
                               OutIndexType& filtered_indices, OutNormalsType& filtered_normals,
                               std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator,
                               int num_threads, int random_seed)
{

    auto set_parallel = OmpSetParallelism(num_threads);

    filtered_indices.resize(points.rows(), 1);
    filtered_normals.resize(points.rows(), 3);
    int num_inserted = 0;

    bool threw_exception = false;

    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        using index_t = typename OutIndexType::Scalar;
        std::vector<index_t> selected_indices;
        std::vector<Eigen::RowVector3d> estimated_normals;

        if (random_seed > 0) {
            #if defined(_OPENMP)
            std::srand(random_seed ^ omp_get_thread_num());  // Hopefully everything uses rand()
            #endif
        }

        #if defined(_OPENMP)
        #pragma omp for nowait
        #endif
        for (int i = 0; i < points.rows(); i++) {
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

            std::tuple<bool, Eigen::RowVector3d> res = estimator(i);
            if (std::get<0>(res)) {
                selected_indices.push_back(i);
                estimated_normals.push_back(std::get<1>(res));
            }
        }

        #if defined(_OPENMP)
        #pragma omp critical
        #endif
        {
            if (!threw_exception) {
                for (int i = 0; i < selected_indices.size(); i++) {
                    filtered_indices(num_inserted, 0) = selected_indices[i];
                    for (int j = 0; j < 3; j++) {
                        filtered_normals(num_inserted, j) = estimated_normals[i][j];
                    }
                    num_inserted += 1;
                }
            }
        }
    }

    if (threw_exception) {
        throw pybind11::error_already_set();
    }
    filtered_indices.conservativeResize(num_inserted, 1);
    filtered_normals.conservativeResize(num_inserted, 3);
}

template <typename InPointsType, typename InViewDirsType, typename OutIndexType, typename OutNormalsType>
void estimate_normals_single_thread(const InPointsType& points, const InViewDirsType& view_dirs,
                                    OutIndexType& filtered_indices, OutNormalsType& filtered_normals,
                                    std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator,
                                    int random_seed)
{
    if (random_seed > 0) {
        std::srand(random_seed);  // Hopefully everything uses rand()
    }

    filtered_indices.resize(points.rows(), 1);
    filtered_normals.resize(points.rows(), 3);
    int num_filtered_points = 0;
    for (int i = 0; i < points.rows(); i++) {
        if (PyErr_CheckSignals() != 0) { throw pybind11::error_already_set(); }
        std::tuple<bool, Eigen::RowVector3d> res = estimator(i);
        if (std::get<0>(res)) {
            Eigen::RowVector3d normal = std::get<1>(res);
            filtered_indices(num_filtered_points, 0) = i;
            for (int j = 0; j < 3; j++) {
                filtered_normals(num_filtered_points, j) = normal[j];
            }
            num_filtered_points += 1;
        }
    }

    filtered_indices.conservativeResize(num_filtered_points, 1);
    filtered_normals.conservativeResize(num_filtered_points, 3);
}


template <typename InPointsType, typename InViewDirsType, typename OutIndexType, typename OutNormalsType>
void estimate_normals(const InPointsType& points, const InViewDirsType& view_dirs,
                      OutIndexType& filtered_indices, OutNormalsType& filtered_normals,
                      std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator,
                      int num_threads, int random_seed)
{
    #ifdef _OPENMP
    const int MIN_PARALLEL_INPUT_SIZE = 1000000;
    if (num_threads != 0 && num_threads != 1 && points.rows() > MIN_PARALLEL_INPUT_SIZE) {
        estimate_normals_parallel(points, view_dirs, filtered_indices, filtered_normals, estimator, num_threads, random_seed);
    } else {
        estimate_normals_single_thread(points, view_dirs, filtered_indices, filtered_normals, estimator, random_seed);
    }
    #else
        estimate_normals_single_thread(points, view_dirs, filtered_indices, filtered_normals, estimator, random_seed);
    #endif
}
} // namespace


npe_function(estimate_point_cloud_normals_ball_internal)
npe_arg(points, dense_float, dense_double)
npe_arg(view_dirs, npe_matches(points))
npe_arg(radius, double)
npe_arg(min_pts_per_ball, int)
npe_default_arg(max_pts_per_ball, int, -1)
npe_default_arg(drop_angle_threshold, double, M_PI / 2.0)
npe_default_arg(max_points_per_leaf, int, 10)
npe_default_arg(num_threads, int, 0)
npe_default_arg(weight_function, std::string, "constant")
npe_default_arg(random_seed, int, -1)
npe_begin_code()
{
    // Check that the inputs are good
    if (radius <= 0.0) {
        throw pybind11::value_error("Invalid radius (" + std::to_string(radius) + ") must be greater than 0.");
    }
    if (min_pts_per_ball < 3) {
        throw pybind11::value_error("Invalid min_pts_per_ball (" + std::to_string(min_pts_per_ball) +
                                    ") must be greater than 3.");
    }
    if (max_pts_per_ball > 0 && max_pts_per_ball < 3) {
        throw pybind11::value_error("Invalid max_pts_per_ball (" + std::to_string(max_pts_per_ball) +
                            ") must either be negative (no max) or a number greater than 3.");
    }
    validate_input(points, view_dirs);

    // Build kdtree
    using namespace nanoflann;
    using MatrixType = EigenDenseLike<npe_Matrix_points>;
    using KdTreeType = nanoflann::KDTreeEigenMatrixAdaptor<decltype(points), 3, nanoflann::metric_L2_Simple>;
    KdTreeType tree(3, std::cref(points), max_points_per_leaf /* max leaf */);

    // Create a callback function which estimates the normal of the i^th point in the input point cloud
    auto wendland_rbf = [](double d, double rad) {
        double r = d / rad;
        double val1 = 1.0 - r;
        double val2 = 4 * r + 1.0;
        return val1 * val1 * val1 * val1 * val2;
    };
    auto constant_rbf = [](double d, double rad) {
        return 1.0;
    };
    std::function<double(double, double)> weight_callback;
    if (weight_function == "constant") {
        weight_callback  = constant_rbf;
    } else if (weight_function == "rbf") {
        weight_callback = wendland_rbf;
    } else {
        throw pybind11::value_error("Invalid weight_function, must be one of 'constant' or 'rbf'.");
    }

    auto normal_estimator = [&](int pt_index) {
        return estimate_local_normal_rbf(tree, points, view_dirs, pt_index, radius,
                                         min_pts_per_ball, max_pts_per_ball, drop_angle_threshold,
                                         weight_callback);
    };

    // Estimate the normals at every point (possibly in parallel)
    Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> filtered_indices;
    EigenDenseLike<npe_Matrix_view_dirs> filtered_normals;
    estimate_normals(points, view_dirs, filtered_indices, filtered_normals, normal_estimator,
                     num_threads, random_seed);

    return std::make_tuple(npe::move(filtered_indices), npe::move(filtered_normals));
}
npe_end_code()


npe_function(estimate_point_cloud_normals_knn_internal)
npe_arg(points, dense_float, dense_double)
npe_arg(view_dirs, npe_matches(points))
npe_arg(num_neighbors, int)
npe_default_arg(max_points_per_leaf, int, 10)
npe_default_arg(drop_angle_threshold, double, M_PI / 2.0)
npe_default_arg(num_threads, int, 0)
npe_default_arg(random_seed, int, -1)
npe_begin_code()
{
    if (num_neighbors <= 0) {
        throw pybind11::value_error("Invalid number of neighbors (" +
                                    std::to_string(num_neighbors) +
                                    ") must be greater than 0.");
    }
    validate_input(points, view_dirs);

    // Build kdtree
    using namespace nanoflann;
    using MatrixType = EigenDenseLike<npe_Matrix_points>;
    using KdTreeType = nanoflann::KDTreeEigenMatrixAdaptor<decltype(points), 3, nanoflann::metric_L2_Simple>;
    KdTreeType tree(3, std::cref(points), max_points_per_leaf /* max leaf */);

    // Create a callback function which estimates the normal of the i^th point in the input point cloud
    auto normal_estimator = [&](int pt_index) {
        return estimate_local_normal_knn(tree, points, view_dirs, pt_index, num_neighbors, drop_angle_threshold);
    };

    Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> filtered_indices;
    EigenDenseLike<npe_Matrix_view_dirs> filtered_normals;
    estimate_normals(points, view_dirs, filtered_indices, filtered_normals, normal_estimator,
                     num_threads, random_seed);

    // Estimate the normals at every point (possibly in parallel)
    return std::make_tuple(npe::move(filtered_indices), npe::move(filtered_normals));
}
npe_end_code()
