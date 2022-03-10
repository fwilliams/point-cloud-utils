import numpy as np


def estimate_point_cloud_normals_knn(points, num_neighbors, view_directions=None,
                                     drop_angle_threshold=np.deg2rad(90.0),
                                     max_points_per_leaf=10, num_threads=-1):
    """
    Estimate normals for a point cloud by locally fitting a plane to the k nearest neighbors of each point.

    This function can optionally consider directions to the sensor for each point to compute neighborhoods of points
    which are all facing the same direction, and align the final normal directions.
    
    Parameters
    ----------
    points: (n, 3)-shaped NumPy array of point positions (each row is a point)
    num_neighbors: Integer number of neighbors to use in each neigghborhood.
    view_directions: (n, 3)-shaped NumPy array or None, representing the unit direction to the sensor for each point.
                      This parameter is used to align the normals and compute neighborhoods of similar facing points.
    drop_angle_threshold: If view_directions is passed in, drop points whose angle between the normal and view direction
                          exceeds drop_angle_threshold (in radians). Useful for filtering out low quality points.
    max_points_per_leaf: Maximum number of points in each leaf node of the KD-tree used for nearest neighbor queries.
                          Tuning this can potentially improve performance on large point clouds.
    num_threads: Number of threads used to parallelize computation. If set to 0 ir 1, will run in single threaded mode.
                 If set to a positive number t > 1, will use t threads.
                 If set to -1, will use #processors threads for inputs greater than 1 million points.

    Returns
    -------
    A tuple (p, n) of filtered points and normals (filtered because some points may be discarded), where:
      - p is an (m, 3)-shaped Numpy array of 3d points
      - n is an (m, 3)-shaped Numpy array of unit normals for each point in p

    See Also
    --------
    estimate_pointcloud_normals_ball
    """
    from ._pcu_internal import estimate_point_cloud_normals_knn_internal

    if view_directions is None:
        view_directions = np.zeros([0, 3], dtype=points.dtype)

    if type(view_directions) != np.ndarray:
        raise ValueError("Invalid type for view_directions, must be None or a NumPy array, "
                         "but got " + type(view_directions) + ".")
    if type(points) != np.ndarray:
        raise ValueError("Invalid type for points, must be None or a NumPy array, but got " + type(points) + ".")

    if len(points.shape) != 2:
        raise ValueError("Invalid shape for points, must be (n, 3) but got " + str(points.shape))
    if points.shape[-1] != 3:
        raise ValueError("Invalid shape for points, must be (n, 3) but got " + str(points.shape))

    if len(view_directions.shape) != 2:
        raise ValueError("Invalid shape for view_directions, must be (n, 3) but got " + str(view_directions.shape))
    if points.shape[-1] != 3:
        raise ValueError("Invalid shape for view_directions, must be (n, 3) but got " + str(points.shape))

    points, normals = estimate_point_cloud_normals_knn_internal(points, view_directions, num_neighbors,
                                                                max_points_per_leaf, drop_angle_threshold,
                                                                num_threads)

    return points, normals


def estimate_point_cloud_normals_ball(points, ball_radius, view_directions=None,
                                      drop_angle_threshold=np.deg2rad(90.0),
                                      min_pts_per_ball=3,
                                      weight_function="rbf",
                                      max_points_per_leaf=10, num_threads=-1):
    """
    Estimate normals for a point cloud by locally fitting a plane to all points within a radius of each point
    (possibly weighted by a radial basis function).

    This function can optionally consider directions to the sensor for each point to compute neighborhoods of points
    which are all facing the same direction, and align the final normal directions.

    Parameters
    ----------
    points: (n, 3)-shaped NumPy array of point positions (each row is a point)
    ball_radius: The radius of each neighborhood used to estimate normals
    view_directions: (n, 3)-shaped NumPy array or None, representing the unit direction to the sensor for each point.
                      This parameter is used to align the normals and compute neighborhoods of similar facing points.
    drop_angle_threshold: If view_directions is passed in, drop points whose angle between the normal and view direction
                          exceeds drop_angle_threshold (in radians). Useful for filtering out low quality points.
    min_pts_per_ball: Discard points whose neighborhood contains fewer than min_pts_per_ball points.
    weight_function: Weighting function for points in a neighborhood. Must be one of 'constant' or 'rbf' where:
      - 'rbf' weights points as (1 - d/r)^4 * (4* d/r + 1) where d = distance to the center point and r = ball_radius
      - 'constant' weights points as 1.0 for every point
    max_points_per_leaf: Maximum number of points in each leaf node of the KD-tree used for nearest neighbor queries.
                          Tuning this can potentially improve performance on large point clouds.
    num_threads: Number of threads used to parallelize computation. If set to 0 ir 1, will run in single threaded mode.
                 If set to a positive number t > 1, will use t threads.
                 If set to -1, will use #processors threads for inputs greater than 1 million points.

    Returns
    -------
    A tuple (p, n) of filtered points and normals (filtered because some points may be discarded), where:
      - p is an (m, 3)-shaped Numpy array of 3d points
      - n is an (m, 3)-shaped Numpy array of unit normals for each point in p

    See Also
    --------
    estimate_pointcloud_normals_ball
    """
    from ._pcu_internal import estimate_point_cloud_normals_ball_internal

    if view_directions is None:
        view_directions = np.zeros([0, 3], dtype=points.dtype)

    if type(view_directions) != np.ndarray:
        raise ValueError("Invalid type for view_directions, must be None or a NumPy array, "
                         "but got " + type(view_directions) + ".")
    if type(points) != np.ndarray:
        raise ValueError("Invalid type for points, must be None or a NumPy array, but got " + type(points) + ".")

    if len(points.shape) != 2:
        raise ValueError("Invalid shape for points, must be (n, 3) but got " + str(points.shape))
    if points.shape[-1] != 3:
        raise ValueError("Invalid shape for points, must be (n, 3) but got " + str(points.shape))

    if len(view_directions.shape) != 2:
        raise ValueError("Invalid shape for view_directions, must be (n, 3) but got " + str(view_directions.shape))
    if points.shape[-1] != 3:
        raise ValueError("Invalid shape for view_directions, must be (n, 3) but got " + str(points.shape))

    points, normals = estimate_point_cloud_normals_ball_internal(points, view_directions, ball_radius, min_pts_per_ball,
                                                                 drop_angle_threshold, max_points_per_leaf,
                                                                 num_threads, weight_function)

    return points, normals
