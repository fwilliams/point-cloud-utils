from ._pcu_internal import sample_mesh_poisson_disk, sample_mesh_random, \
    downsample_point_cloud_poisson_disk, estimate_mesh_normals, \
    k_nearest_neighbors, one_sided_hausdorff_distance, \
    morton_encode, morton_decode, morton_knn, \
    lloyd_2d, lloyd_3d, voronoi_centroids_unit_cube, sample_mesh_lloyd, \
    deduplicate_point_cloud, deduplicate_mesh_vertices, signed_distance_to_mesh, \
    closest_points_on_mesh, connected_components, ray_mesh_intersection, laplacian_smooth_mesh, \
    make_mesh_watertight, mesh_principal_curvatures, mesh_mean_and_gaussian_curvatures, \
    morton_add, morton_subtract, \
    sparse_voxel_grid_boundary, marching_cubes_sparse_voxel_grid

from ._sinkhorn import *
from ._mesh_io import *
from ._octree import *
from ._pointcloud_normals import estimate_point_cloud_normals_knn, estimate_point_cloud_normals_ball
from ._ray_mesh_intersector import RayMeshIntersector
from ._ray_point_cloud_intersector import ray_surfel_intersection, surfel_geometry, RaySurfelIntersector


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


def sparse_voxel_grid_from_pointcloud(points, voxel_size, origin=(0.0, 0.0, 0.0)):
    """
    Construct a sparse voxel grid containing a point cloud, optionally returning the voxel index for each point.

    Parameters
    ----------
    points: An (n, 3) shaped array of points (one per row)
    voxel_size: Either a scalar representing the sidelength of a voxel or a triple representing the sidelength of a
                voxel along each axis.
    origin: The origin of the point cloud corresponding to the bottom, back, left corner of the (0, 0, 0) voxel.

    Returns
    -------
    grid_coordinates: An (n, 3)-shaped array of integer grid coordinates corresponding to voxels in the sparse grid
    point_to_vox_idx: An (#points,)-shaped array where the i^th entry point_to_vox_idx[i] is the index into
                      grid_coordinates of of the voxel containing points[i]
    """
    from ._pcu_internal import sparse_voxel_grid_from_pointcloud_internal

    if np.isscalar(voxel_size):
        voxel_size = np.array([voxel_size] * 3)
    elif len(voxel_size) == 1:
        voxel_size = np.array([np.squeeze(voxel_size).reshape((1,))[0]] * 3)
    elif len(voxel_size) == 3:
        voxel_size = np.squeeze(voxel_size).reshape(3,)
    else:
        raise ValueError("voxel_size must be a scalar or an array/tuple/list of length 3")

    if not np.all(np.isreal(voxel_size)):
        raise ValueError("voxel_size must have a real numeric type")

    if len(origin) == 3:
        origin = np.squeeze(origin).reshape(3,)
    else:
        raise ValueError("origin must be an array/tuple/list of length 3")

    if not np.all(np.isreal(origin)):
        raise ValueError("origin must have a real numeric type")

    grid_coords, point_to_vox = sparse_voxel_grid_from_pointcloud_internal(
        points, voxel_size[0], voxel_size[1], voxel_size[2],
        origin[0], origin[1], origin[2], return_point_to_vox=True
    )

    return grid_coords, point_to_vox


def subdivide_sparse_voxel_grid(grid_coordinates, num_subdivs=1):
    """
    Subdivide a sparse voxel grid such that every subdivision round splits each voxel into 8 subvoxels.

    Parameters
    ----------
    grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
    num_subdivs: The number of subdivision iterations to perform

    Returns
    -------
    new_grid_coordinates: An (m, 3) shaped array containing (fractional) coordinates representing the subdivided grid
                          voxel grid. To convert to integer coordinates you can run:
                          integer_grid_coords = quantize_subdivided_sparse_voxel_grid(num_subdivs)

    """

    from ._pcu_internal import subdivide_sparse_voxel_grid_internal

    gc = grid_coordinates
    for _ in range(num_subdivs):
        gc = subdivide_sparse_voxel_grid_internal(gc)

    return gc / 2 ** num_subdivs


def quantize_subdivided_sparse_voxel_grid(fractional_grid_coordinates, num_subdivs):
    """
    subdivide_sparse_voxel_grid returns fractional voxel coordinates (since we are subdividing an integer grid).
    This function transforms these fractional coordinates into integer coordinates. This operation can be undone
    with unquantize_subdivided_sparse_voxel_grid

    Parameters
    ----------
    fractional_grid_coordinates: An (n, 3) shaped array of fractional voxels coordinates returned by
                                 subdivide_sparse_voxel_grid
    num_subdivs: The number of subdivision iterations used to compute fractional_grid_coordinates

    Returns
    -------
    quantized_grid_coordinates: An (n, 3) shaped array of quantized grid coordinates

    See Also
    --------
    subdivide_sparse_voxel_grid
    unquantize_subdivided_sparse_voxel_grid
    """
    return (fractional_grid_coordinates * (2 ** num_subdivs) - 1).astype(np.int32) // 2


def unquantize_subdivided_sparse_voxel_grid(quantized_grid_coordinates, num_subdivs):
    """
    Undo quantization done by quantize_subdivided_sparse_voxel_grid.

    Parameters
    ----------
    quantized_grid_coordinates: An (n, 3) shaped array of quantized voxels coordinates returned by
                                quantize_subdivided_sparse_voxel_grid
    num_subdivs: The number of subdivision iterations used to compute the original fractional coordinates

    Returns
    -------
    fractional_grid_coordinates: An (n, 3) shaped array of fractional grid coordinates

    See Also
    --------
    subdivide_sparse_voxel_grid
    quantize_subdivided_sparse_voxel_grid
    """
    return 2.0 * quantized_grid_coordinates / (2 ** num_subdivs) + 1.0 / (2 ** num_subdivs)


def dilate_sparse_voxel_grid(grid_coordinates, count=1):
    """
    Perform binary dilation on a sparse voxel grid, using a 3x3 filter of ones applied count times.
    i.e. convolve the filter:
        [[1, 1, 1],
         [1, 1, 1],
         [1, 1, 1]]
    with the sparse voxel grid.

    Parameters
    ----------
    grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
    count: The number of iterations of binary dilation to run

    Returns
    -------
    dilated_grid_coordinates: An (m, 3) integer array encoding the coordinates of the sparse voxel grid after dilation
    """
    from ._pcu_internal import dilate_sparse_voxel_grid_internal
    for _ in range(count):
        grid_coordinates = dilate_sparse_voxel_grid_internal(grid_coordinates)

    return grid_coordinates


def erode_sparse_voxel_grid(grid_coordinates, count=1):
    """
    Perform binary erosion on a sparse voxel grid, using a 3x3 filter of ones applied count times.
    i.e. take the negation of the convolution of the filter:
        [[1, 1, 1],
         [1, 1, 1],
         [1, 1, 1]]
    with the sparse voxel grid.

    Parameters
    ----------
    grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
    count: The number of iterations of binary erosion to run

    Returns
    -------
    eroded_grid_coordinates: An (m, 3) integer array encoding the coordinates of the sparse voxel grid after erosion
    """
    from ._pcu_internal import erode_sparse_voxel_grid_internal
    for _ in range(count):
        grid_coordinates = erode_sparse_voxel_grid_internal(grid_coordinates)

    return grid_coordinates


def sparse_voxel_grid_to_hex_mesh(grid_coordinates, voxel_size=1.0, origin=(0.0, 0.0, 0.0), dtype=np.float64):
    """
    Convert a sparse voxel grid where the voxel indices represent the centers of each voxel coordinate into a
    hexahedral mesh.

    i.e. a pair of arrays (v, c) where v has shape (n, 3) and each row is a vertex position
         and c has shape (m, 8) where the i^th row c[i, :] are the indices into v of the 8 corners of the i^th cube

    Parameters
    ----------
    grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
    voxel_size: Either a scalar encoding the length of one voxel along an axis or a triple encoding the length of a
                voxel on each of 3 axes
    origin: A triple encoding a global offset for the hex mesh (i.e. all vertices will have origin added to them).
            Default is (0.0, 0.0, 0.0).
    dtype: The scalar type of vertices to return.
           Default is np.float64.

    Returns
    -------
    vertices: An (n, 3) shaped array where each row vertices[i, :] is a vertex position in the hex mesh
    cubes: An (m, 8) shaped array of indices into v where each row cubes[j, :] are indices of the 8 points forming
           the j^th cube
    """
    from ._pcu_internal import sparse_voxel_grid_to_hex_mesh_internal
    if np.isscalar(voxel_size):
        voxel_size = np.array([voxel_size] * 3)
    elif len(voxel_size) == 1:
        voxel_size = np.array([np.squeeze(voxel_size).reshape((1,))[0]] * 3)
    elif len(voxel_size) == 3:
        voxel_size = np.squeeze(voxel_size).reshape(3,)
    else:
        raise ValueError("voxel_size must be a scalar or an array/tuple/list of length 3")
    if not np.all(np.isreal(voxel_size)):
        raise ValueError("voxel_size must have a real numeric type")

    if len(origin) == 3:
        origin = np.squeeze(origin).reshape(3,)
    else:
        raise ValueError("origin must be an array/tuple/list of length 3")
    if not np.all(np.isreal(origin)):
        raise ValueError("origin must have a real numeric type")

    ret_v, ret_c = sparse_voxel_grid_to_hex_mesh_internal(
        grid_coordinates, voxel_size[0], voxel_size[1], voxel_size[2],
        origin[0], origin[1], origin[2])

    return ret_v.astype(dtype), ret_c

