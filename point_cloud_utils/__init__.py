from warnings import warn


from ._pcu_internal import sample_mesh_poisson_disk, sample_mesh_random, \
    downsample_point_cloud_poisson_disk, estimate_mesh_vertex_normals, \
    estimate_mesh_face_normals, orient_mesh_faces, \
    k_nearest_neighbors, one_sided_hausdorff_distance, \
    morton_encode, morton_decode, morton_knn, \
    lloyd_2d, lloyd_3d, voronoi_centroids_unit_cube, sample_mesh_lloyd, \
    deduplicate_point_cloud, deduplicate_mesh_vertices, signed_distance_to_mesh, \
    closest_points_on_mesh, connected_components, ray_mesh_intersection, laplacian_smooth_mesh, \
    make_mesh_watertight, mesh_principal_curvatures, \
    morton_add, morton_subtract, point_cloud_fast_winding_number, \
    sparse_voxel_grid_boundary, marching_cubes_sparse_voxel_grid, decimate_triangle_mesh, \
    remove_unreferenced_mesh_vertices, mesh_face_areas, triangle_soup_fast_winding_number, \
    _voxel_mesh_internal

from ._sinkhorn import *
from ._mesh_io import *
from ._pointcloud_normals import estimate_point_cloud_normals_knn, estimate_point_cloud_normals_ball
from ._ray_mesh_intersector import RayMeshIntersector
from ._ray_point_cloud_intersector import ray_surfel_intersection, RaySurfelIntersector
from ._point_cloud_geometry import voxel_grid_geometry, pointcloud_sphere_geometry, pointcloud_surfel_geometry
from ._voxels import flood_fill_3d, voxelize_triangle_mesh
from ._mesh_primitives import sphere_mesh, cube_mesh, cylinder_mesh, cone_mesh

MORTON_MIN = -1048576
MORTON_MAX = 1048576


def mesh_mean_and_gaussian_curvatures(v, f, r=-1.0):
    """
    Estimate mean and Gaussian curvatures for a mesh

    Args:
        v : \#v by 3 Matrix of mesh vertex 3D positions
        f : \#f by 3 Matrix of face (triangle) indices
        r : optional floating point radius of neighborhood to consider when estimating curvature
            If set to a positive value, will use a more robust curvature estimation method (but may require some tuning)

    Returns:
        kh : an array of shape (#v,) of per-vertex mean curvatures
        kg : an array of shape (#v,) of per-vertex Gaussian curvatures

    """
    k1, k2, _, _ = mesh_principal_curvatures(v, f, r)

    return 0.5 * (k1 + k2), (k1 * k2)


def hausdorff_distance(x, y, return_index=False, squared_distances=False, max_points_per_leaf=10):
    """
    Compute the Hausdorff distance between x and y

    Args:
        x : n by 3 array of representing a set of n points (each row is a point of dimension 3)
        y : m by 3 array of representing a set of m points (each row is a point of dimension 3)
        return_index : Optionally return the index pair `(i, j)` into x and y such that `x[i, :]` and `y[j, :]` are the two points with maximum shortest distance.
        squared_distances : If set to True, then return squared L2 distances. Default is False.
        max_points_per_leaf : The maximum number of points per leaf node in the KD tree used by this function. Default is 10.

    Returns:
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

    Args:
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
    Returns:
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


def downsample_point_cloud_on_voxel_grid(voxel_size, points, *args, min_bound=None, max_bound=None, min_points_per_voxel=1):
    """
    Downsample a point set to conform with a voxel grid by taking the average of points within each voxel.

    Args:
        voxel_size : a scalar representing the size of each voxel or a 3 tuple representing the size per axis of each voxel.
        points: a [#v, 3]-shaped array of 3d points.
        *args: Any additional arguments of shape [#v, *] are treated as attributes and will averaged into each voxel along with the points
            These will be returns
        min_bound: a 3 tuple representing the minimum coordinate of the voxel grid or None to use the bounding box of the
                input point cloud.
        max_bound: a 3 tuple representing the maximum coordinate of the voxel grid or None to use the bounding box of the
                input point cloud.
        min_points_per_voxel: If a voxel contains fewer than this many points, then don't include the points in that voxel
                            in the output.

    Returns:
        A tuple (v, attrib0, attrib1, ....) of downsampled points, and point attributes.
        Attributes are returned in the same order they are passed in.
        If no attributes are passed in, then this function simply returns vertices.
    """
    from ._pcu_internal import downsample_point_cloud_voxel_grid_internal

    if np.isscalar(voxel_size):
        voxel_size = np.array([voxel_size] * 3)
    else:
        voxel_size = np.array(voxel_size)
        if len(voxel_size) != 3:
            raise ValueError("Invalid voxel size must be a 3-tuple or a single float")

    if type(points) != np.ndarray:
        raise ValueError("points must be a numpy array but got type " + str(type(points)))
    attribs = []
    for i, arg in enumerate(args):
        if type(arg) != np.ndarray:
            raise ValueError("Additional arguments after points and before keyword arguments must be numpy arrays")
        if arg.shape[0] != points.shape[0]:
            raise ValueError("Attribute " + str(i) + " must have same first dimension as number of points (" + str(points.shape) + " but got attrib.shape = " + str(arg.shape))
        attribs.append(arg)

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


    a0 = attribs[0] if len(attribs) > 0 else np.zeros([0, 0])
    ret_v, ret_a0 = downsample_point_cloud_voxel_grid_internal(points, a0,
                                                               voxel_size[0], voxel_size[1], voxel_size[2],
                                                               min_bound[0], min_bound[1], min_bound[2],
                                                               max_bound[0], max_bound[1], max_bound[2],
                                                               min_points_per_voxel)
    if ret_a0.size > 0:
        ret = [ret_v, ret_a0]
    else:
        ret = [ret_v]
    for i in range(1, len(attribs)):
        _, ret_ai = downsample_point_cloud_voxel_grid_internal(points, attribs[i],
                                                               voxel_size[0], voxel_size[1], voxel_size[2],
                                                               min_bound[0], min_bound[1], min_bound[2],
                                                               max_bound[0], max_bound[1], max_bound[2],
                                                               min_points_per_voxel)
        ret.append(ret_ai)

    if len(ret) > 1:
        return tuple(ret)
    else:
        return ret_v


def interpolate_barycentric_coords(f, fi, bc, attribute):
    """
    Interpolate an attribute stored at each vertex of a mesh across the faces of a triangle mesh using
    barycentric coordinates

    Args:
        f : a (#faces, 3)-shaped NumPy array of mesh faces (indexing into some vertex array).
        fi: a (#attribs,)-shaped NumPy array of indexes into f indicating which face each attribute lies within.
        bc: a (#attribs, 3)-shaped NumPy array of barycentric coordinates for each attribute
        attribute: a (#vertices, dim)-shaped NumPy array of attributes at each of the mesh vertices

    Returns:
        A (#attribs, dim)-shaped array of interpolated attributes.
    """
    return (attribute[f[fi]] * bc[:, :, np.newaxis]).sum(1)
