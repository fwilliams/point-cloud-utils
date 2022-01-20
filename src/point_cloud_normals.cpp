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
#include "common.h"


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
    const size_t founds = tree.index->radiusSearch(query_point.data(), ball_radius, out_nbrs, SearchParams());
    if (founds < min_pts_per_ball) {
        return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
    }

    int num_valid = 0;
    Eigen::MatrixXd neighbors(founds, 3);
    for (int i = 0; i < founds; i++) {
        const int nbr_idx = out_nbrs[i].first;
        const double nbr_dist = out_nbrs[i].second;
        const double weight = weighting_function(sqrt(nbr_dist), ball_radius);

        if (has_view_directions) {
            // Add neighboring point if its view direction has the same sign as this point
            Eigen::RowVector3d d_i(view_dirs(nbr_idx, 0), view_dirs(nbr_idx, 1), view_dirs(nbr_idx, 2));
            if (d_i.dot(dir) > 0.0) {
                for (int j = 0; j < 3; j++) { neighbors(num_valid, j) = points(nbr_idx, j) * weight; }
                num_valid += 1;
            }
        } else {
            // If we don't have view directions, just add the whole neighborhood
            for (int j = 0; j < 3; j++) { neighbors(num_valid, j) = points(nbr_idx, j) * weight; }
            num_valid += 1;
        }
    }

    // Filter point/normal if it does not have enought neighbors
    if (num_valid < min_pts_per_ball) {
        return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
    }

    neighbors.conservativeResize(num_valid, 3);

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

    std::vector<IndexType> out_indices(2 * num_neighbors);
    std::vector<ScalarType> out_dists_sqr(2 * num_neighbors);
    const size_t founds = tree.index->knnSearch(query_point.data(), 2 * num_neighbors,
                                                     out_indices.data(), out_dists_sqr.data());
    int num_valid = 0;
    Eigen::MatrixXd neighbors(founds, 3);
    for (int i = 0; i < founds; i++) {
        const int nbr_idx = out_indices[i];
        const double nbr_dist = out_dists_sqr[i];

        if (has_view_directions) {
            // Add neighboring point if its view direction has the same sign as this point
            Eigen::RowVector3d d_i(view_dirs(nbr_idx, 0), view_dirs(nbr_idx, 1), view_dirs(nbr_idx, 2));
            if (d_i.dot(dir) > 0.0) {
                for (int j = 0; j < 3; j++) { neighbors(num_valid, j) = points(nbr_idx, j); }
                num_valid += 1;
            }
        } else {
            // If we don't have view directions, just add the whole neighborhood
            for (int j = 0; j < 3; j++) { neighbors(num_valid, j) = points(nbr_idx, j); }
            num_valid += 1;
        }
    }

    // If we don't have enough neighbors facing in the same direction, then discard this point
    if (num_valid < num_neighbors) {
        return std::make_tuple(false, Eigen::RowVector3d(0., 0., 0.));
    }

    neighbors.conservativeResize(num_neighbors, 3);

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

template <typename InPointsType, typename InViewDirsType, typename OutPointsType, typename OutNormalsType>
void estimate_normals_parallel(const InPointsType& points, const InViewDirsType& view_dirs,
                               OutPointsType& filtered_points, OutNormalsType& filtered_normals,
                               std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator, int num_threads)
{

    #if defined(_OPENMP)
    if (num_threads < 0) {
        const auto processor_count = std::thread::hardware_concurrency();
        omp_set_num_threads(processor_count / 2);
        // std::cout << "PROCESSOR COUNT " << processor_count << std::endl;
    } else {
        omp_set_num_threads(num_threads);
    }
    // std::cout << "OPENMP NUMBER OF THREADS: " << omp_get_num_threads() << " / " << omp_get_max_threads() << std::endl;
    #endif

    std::vector<bool> mask(points.rows());
    std::vector<Eigen::RowVector3d> estimated_normals(points.rows());

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for (int i = 0; i < points.rows(); i++) {
        std::tuple<bool, Eigen::RowVector3d> res = estimator(i);
        mask[i] = std::get<0>(res);
        estimated_normals[i] = std::get<1>(res);
    }

    filtered_points.resize(points.rows(), 3);
    filtered_normals.resize(view_dirs.rows(), 3);
    int num_filtered_points = 0;
    for (int i = 0; i < points.rows(); i++) {
        if (mask[i]) {
            const Eigen::RowVector3d n_i = estimated_normals[i];
            for (int j = 0; j < 3; j++) {
                filtered_normals(num_filtered_points, j) = n_i[j];
                filtered_points(num_filtered_points, j) = points(i, j);
            }
            num_filtered_points += 1;
        }
    }
    filtered_points.conservativeResize(num_filtered_points, 3);
    filtered_normals.conservativeResize(num_filtered_points, 3);
}

template <typename InPointsType, typename InViewDirsType, typename OutPointsType, typename OutNormalsType>
void estimate_normals_single_thread(const InPointsType& points, const InViewDirsType& view_dirs,
                                    OutPointsType& filtered_points, OutNormalsType& filtered_normals,
                                    std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator)
{
    filtered_points.resize(view_dirs.rows(), 3);
    filtered_normals.resize(view_dirs.rows(), 3);
    int num_filtered_points = 0;
    for (int i = 0; i < points.rows(); i++) {
        std::tuple<bool, Eigen::RowVector3d> res = estimator(i);
        if (std::get<0>(res)) {
            Eigen::RowVector3d normal = std::get<1>(res);
            for (int j = 0; j < 3; j++) {
                filtered_normals(num_filtered_points, j) = normal[j];
                filtered_points(num_filtered_points, j) = points(i, j);
            }
            num_filtered_points += 1;
        }
    }

    filtered_points.conservativeResize(num_filtered_points, 3);
    filtered_normals.conservativeResize(num_filtered_points, 3);
}


template <typename InPointsType, typename InViewDirsType, typename OutPointsType, typename OutNormalsType>
void estimate_normals(const InPointsType& points, const InViewDirsType& view_dirs,
                      OutPointsType& filtered_points, OutNormalsType& filtered_normals,
                      std::function<std::tuple<bool, Eigen::RowVector3d>(int)> estimator, int num_threads)
{
    #ifdef _OPENMP
    const int MIN_PARALLEL_INPUT_SIZE = 1000000;
    if (num_threads != 0 && num_threads != 1 && points.rows() > MIN_PARALLEL_INPUT_SIZE) {
        estimate_normals_parallel(points, view_dirs, filtered_points, filtered_normals, estimator, num_threads);
    } else {
        estimate_normals_single_thread(points, view_dirs, filtered_points, filtered_normals, estimator);
    }
    #else
        estimate_normals_single_thread(points, view_dirs, filtered_points, filtered_normals, estimator);
    #endif
}
} // namespace


npe_function(estimate_point_cloud_normals_ball_internal)
npe_arg(points, dense_float, dense_double)
npe_arg(view_dirs, npe_matches(points))
npe_arg(radius, double)
npe_arg(min_pts_per_ball, int)
npe_default_arg(drop_angle_threshold, double, M_PI / 2.0)
npe_default_arg(max_points_per_leaf, int, 10)
npe_default_arg(num_threads, int, 0)
npe_default_arg(weight_function, std::string, "rbf")
npe_begin_code()
{
    // Check that the inputs are good
    if (radius <= 0.0) {
        throw pybind11::value_error("Invalid radius (" + std::to_string(radius) + ") must be greater than 0.");
    }
    validate_input(points, view_dirs);

    // Build kdtree
    using namespace nanoflann;
    using MatrixType = EigenDenseLike<npe_Matrix_points>;
    using KdTreeType = nanoflann::KDTreeEigenMatrixAdaptor<MatrixType, 3, nanoflann::metric_L2_Simple>;
    EigenDenseLike<npe_Matrix_points> pts = points;  // nanoflann doesn't work with Eigen::map so copy :(
    KdTreeType tree(3, std::cref(pts), max_points_per_leaf /* max leaf */);
    tree.index->buildIndex();

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
                                         min_pts_per_ball, drop_angle_threshold, weight_callback);
    };

    // Estimate the normals at every point (possibly in parallel)
    EigenDenseLike<npe_Matrix_view_dirs> filtered_points;
    EigenDenseLike<npe_Matrix_view_dirs> filtered_normals;
    estimate_normals(points, view_dirs, filtered_points, filtered_normals, normal_estimator, num_threads);

    return std::make_tuple(npe::move(filtered_points), npe::move(filtered_normals));
}
npe_end_code()


npe_function(estimate_point_cloud_normals_knn_internal)
npe_arg(points, dense_float, dense_double)
npe_arg(view_dirs, npe_matches(points))
npe_arg(num_neighbors, int)
npe_default_arg(max_points_per_leaf, int, 10)
npe_default_arg(drop_angle_threshold, double, M_PI / 2.0)
npe_default_arg(num_threads, int, 0)
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
    using KdTreeType = nanoflann::KDTreeEigenMatrixAdaptor<MatrixType, 3, nanoflann::metric_L2_Simple>;
    EigenDenseLike<npe_Matrix_points> pts = points;  // nanoflann doesn't work with Eigen::map so copy :(
    KdTreeType tree(3, std::cref(pts), max_points_per_leaf /* max leaf */);
    tree.index->buildIndex();

    // Create a callback function which estimates the normal of the i^th point in the input point cloud
    auto normal_estimator = [&](int pt_index) {
        return estimate_local_normal_knn(tree, points, view_dirs, pt_index, num_neighbors, drop_angle_threshold);
    };

    EigenDenseLike<npe_Matrix_view_dirs> filtered_points;
    EigenDenseLike<npe_Matrix_view_dirs> filtered_normals;
    estimate_normals(points, view_dirs, filtered_points, filtered_normals, normal_estimator, num_threads);

    // Estimate the normals at every point (possibly in parallel)
    return std::make_tuple(npe::move(filtered_points), npe::move(filtered_normals));
}
npe_end_code()