from .pcu_internal import *
from .sinkhorn import *
import numpy as np


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
    from .pcu_internal import __internal_downsample_point_cloud_voxel_grid

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

    ret_v, ret_n, ret_c = __internal_downsample_point_cloud_voxel_grid(points, normals, colors,
                                                                       *voxel_size, *min_bound, *max_bound,
                                                                       min_points_per_voxel)
    ret = [ret_v, None, None]
    if has_normals:
        ret[1] = ret_n
    if has_colors:
        ret[2] = ret_c
    return tuple(ret)
