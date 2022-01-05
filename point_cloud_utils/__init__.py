from ._pcu_internal import sample_mesh_poisson_disk, sample_mesh_random, \
    downsample_point_cloud_poisson_disk, estimate_point_cloud_normals, estimate_mesh_normals, \
    k_nearest_neighbors, one_sided_hausdorff_distance, \
    morton_encode, morton_decode, morton_knn, \
    lloyd_2d, lloyd_3d, voronoi_centroids_unit_cube, sample_mesh_lloyd, \
    deduplicate_point_cloud, deduplicate_mesh_vertices, signed_distance_to_mesh, \
    closest_points_on_mesh, connected_components, ray_mesh_intersection, laplacian_smooth_mesh, \
    make_mesh_watertight, mesh_principal_curvatures, mesh_mean_and_gaussian_curvatures

from ._sinkhorn import *
from ._mesh_io import *
import numpy as np
from ._octree import *


def hausdorff_distance(x, y, return_index=False, squared_distances=False, max_points_per_leaf=10):
    """
    Compute the Hausdorff distance between x and y

    Parameters
    ----------
    x : n by 3 array of representing a set of n points (each row is a point of dimension 3)
    y : m by 3 array of representing a set of m points (each row is a point of dimension 3)
    return_index : Optionally return the index pair `(i, j)` into x and y such that
                   `x[i, :]` and `y[j, :]` are the two points with maximum shortest distance.
    squared_distances : If set to True, then return squared L2 distances. Default is False.
    max_points_per_leaf : The maximum number of points per leaf node in the KD tree used by this function.
                          Default is 10.

    Returns
    -------
    The largest shortest distance, `d` between each point in `source` and the points in `target`.
    If `return_index` is set, then this function returns a tuple (d, i, j) where `d` is as described above
    and `(i, j)` are such that `source[i, :]` and `target[j, :]` are the two points with maximum shortest
    distance.
    """
    hausdorff_x_to_y, idx_x1, idx_y1 = one_sided_hausdorff_distance(x, y, return_index=True,
                                                                    squared_distances=squared_distances,
                                                                    max_points_per_leaf=max_points_per_leaf)
    hausdorff_y_to_x, idx_y2, idx_x2 = one_sided_hausdorff_distance(y, x, return_index=True,
                                                                    squared_distances=squared_distances,
                                                                    max_points_per_leaf=max_points_per_leaf)

    hausdorff = max(hausdorff_x_to_y, hausdorff_y_to_x)
    if return_index and hausdorff_x_to_y > hausdorff_y_to_x:
        return hausdorff, idx_x1, idx_y1
    elif return_index and hausdorff_x_to_y <= hausdorff_y_to_x:
        return hausdorff, idx_x2, idx_y2
    return hausdorff


def chamfer_distance(x, y, return_index=False, p_norm=2, max_points_per_leaf=10):
    """
    Compute the chamfer distance between two point clouds x, and y

    Parameters
    ----------
    x : A m-sized minibatch of point sets in R^d. i.e. shape [m, n_a, d]
    y : A m-sized minibatch of point sets in R^d. i.e. shape [m, n_b, d]
    return_index: If set to True, will return a pair (corrs_x_to_y, corrs_y_to_x) where
                  corrs_x_to_y[i] stores the index into y of the closest point to x[i]
                  (i.e. y[corrs_x_to_y[i]] is the nearest neighbor to x[i] in y).
                  corrs_y_to_x is similar to corrs_x_to_y but with x and y reversed.
    max_points_per_leaf : The maximum number of points per leaf node in the KD tree used by this function.
                          Default is 10.
    p_norm : Which norm to use. p_norm can be any real number, inf (for the max norm) -inf (for the min norm),
             0 (for sum(x != 0))
    Returns
    -------
    The chamfer distance between x an dy.
    If return_index is set, then this function returns a tuple (chamfer_dist, corrs_x_to_y, corrs_y_to_x) where
    corrs_x_to_y and corrs_y_to_x are described above.
    """

    dists_x_to_y, corrs_x_to_y = k_nearest_neighbors(x, y, k=1,
                                                     squared_distances=False,
                                                     max_points_per_leaf=max_points_per_leaf)
    dists_y_to_x, corrs_y_to_x = k_nearest_neighbors(y, x, k=1,
                                                     squared_distances=False,
                                                     max_points_per_leaf=max_points_per_leaf)

    dists_x_to_y = np.linalg.norm(x[corrs_y_to_x] - y, axis=-1, ord=p_norm).mean()
    dists_y_to_x = np.linalg.norm(y[corrs_x_to_y] - x, axis=-1, ord=p_norm).mean()

    cham_dist = np.mean(dists_x_to_y) + np.mean(dists_y_to_x)

    if return_index:
        return cham_dist, corrs_x_to_y, corrs_y_to_x

    return cham_dist


def downsample_point_cloud_voxel_grid(voxel_size, points, normals=None, colors=None, min_bound=None, max_bound=None,
                                      min_points_per_voxel=1):
    """
    Downsample a point set to conform with a voxel grid by taking the average of points within each voxel.

    Parameters
    ----------
    voxel_size : a scalar representing the size of each voxel or a 3 tuple representing the size per axis of each voxel.
    points: a #v x 3 array of 3d points.
    normals: a #v x 3 array of 3d normals per point or None for no normals.
    colors: a #v x 3 or #v x 4 array of colors per point or None for no colors.
    min_bound: a 3 tuple representing the minimum coordinate of the voxel grid or None to use the bounding box of the
               input point cloud.
    max_bound: a 3 tuple representing the maximum coordinate of the voxel grid or None to use the bounding box of the
               input point cloud.
    min_points_per_voxel: If a voxel contains fewer than this many points, then don't include the points in that voxel
                          in the output.

    Returns
    -------
    A triple (v, n, c) of downsampled vertices, normals and colors. If no vertices or colors are passed in, then
    n and c are None.
    """
    from ._pcu_internal import downsample_point_cloud_voxel_grid_internal

    if np.isscalar(voxel_size):
        voxel_size = np.array([voxel_size] * 3)
    else:
        voxel_size = np.array(voxel_size)
        if len(voxel_size) != 3:
            raise ValueError("Invalid voxel size must be a 3-tuple or a single float")
    has_normals = True
    has_colors = True
    if normals is None:
        normals = np.zeros([0, 0]).astype(points.dtype)
        has_normals = False
    if colors is None:
        colors = np.zeros([0, 0]).astype(points.dtype)
        has_colors = False
    if min_bound is None:
        min_bound = np.min(points, axis=0) - voxel_size * 0.5
    if max_bound is None:
        max_bound = np.max(points, axis=0) + voxel_size * 0.5

    min_bound = np.array(min_bound)
    max_bound = np.array(max_bound)
    if len(min_bound) != 3:
        raise ValueError("min_bound must be a 3 tuple")
    if len(max_bound) != 3:
        raise ValueError("max_bound must be a 3 tuple")

    if np.any(max_bound - min_bound <= 0.0):
        raise ValueError("Invalid min_bound and max_bound. max_bound must be greater than min_bound in all dimensions")

    ret_v, ret_n, ret_c = downsample_point_cloud_voxel_grid_internal(points,
                                                                     normals.astype(points.dtype),
                                                                     colors.astype(points.dtype),
                                                                     voxel_size[0], voxel_size[1], voxel_size[2],
                                                                     min_bound[0], min_bound[1], min_bound[2],
                                                                     max_bound[0], max_bound[1], max_bound[2],
                                                                     min_points_per_voxel)
    ret = [ret_v, None, None]
    if has_normals:
        ret[1] = ret_n
    if has_colors:
        ret[2] = ret_c
    return tuple(ret)


def interpolate_barycentric_coords(f, fi, bc, attribute):
    """
    Interpolate an attribute stored at each vertex of a mesh across the faces of a triangle mesh using
    barycentric coordinates

    Parameters
    ----------
    f : a (#faces, 3)-shaped NumPy array of mesh faces (indexing into some vertex array).
    fi: a (#attribs,)-shaped NumPy array of indexes into f indicating which face each attribute lies within.
    bc: a (#attribs, 3)-shaped NumPy array of barycentric coordinates for each attribute
    attribute: a (#vertices, dim)-shaped NumPy array of attributes at each of the mesh vertices

    Returns
    -------
    A (#attribs, dim)-shaped array of interpolated attributes.
    """
    return (attribute[f[fi]] * bc[:, :, np.newaxis]).sum(1)
