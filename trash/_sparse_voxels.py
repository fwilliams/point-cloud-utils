# def sparse_voxel_grid_from_pointcloud(points, voxel_size, origin=(0.0, 0.0, 0.0)):
#     """
#     Construct a sparse voxel grid containing a point cloud, optionally returning the voxel index for each point.

#     Args:
#         points: An (n, 3) shaped array of points (one per row)
#         voxel_size: Either a scalar representing the sidelength of a voxel or a triple representing the sidelength of a
#                     voxel along each axis.
#         origin: The origin of the point cloud corresponding to the bottom, back, left corner of the (0, 0, 0) voxel.

#     Returns:
#         grid_coordinates: An (n, 3)-shaped array of integer grid coordinates corresponding to voxels in the sparse grid
#         point_to_vox_idx: An (#points,)-shaped array where the i^th entry point_to_vox_idx[i] is the index into
#                         grid_coordinates of of the voxel containing points[i]
#     """
#     from ._pcu_internal import sparse_voxel_grid_from_pointcloud_internal

#     if np.isscalar(voxel_size):
#         voxel_size = np.array([voxel_size] * 3)
#     elif len(voxel_size) == 1:
#         voxel_size = np.array([np.squeeze(voxel_size).reshape((1,))[0]] * 3)
#     elif len(voxel_size) == 3:
#         voxel_size = np.squeeze(voxel_size).reshape(3,)
#     else:
#         raise ValueError("voxel_size must be a scalar or an array/tuple/list of length 3")

#     if not np.all(np.isreal(voxel_size)):
#         raise ValueError("voxel_size must have a real numeric type")

#     if len(origin) == 3:
#         origin = np.squeeze(origin).reshape(3,)
#     else:
#         raise ValueError("origin must be an array/tuple/list of length 3")

#     if not np.all(np.isreal(origin)):
#         raise ValueError("origin must have a real numeric type")

#     grid_coords, point_to_vox = sparse_voxel_grid_from_pointcloud_internal(
#         points, voxel_size[0], voxel_size[1], voxel_size[2],
#         origin[0], origin[1], origin[2], return_point_to_vox=True
#     )

#     return grid_coords, point_to_vox


# def subdivide_sparse_voxel_grid(grid_coordinates, num_subdivs=1):
#     """
#     Subdivide a sparse voxel grid such that every subdivision round splits each voxel into 8 subvoxels.

#     Args:
#         grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
#         num_subdivs: The number of subdivision iterations to perform

#     Returns:
#     new_grid_coordinates: An (m, 3) shaped array containing (fractional) coordinates representing the subdivided grid
#                           voxel grid. To convert to integer coordinates you can run:
#                           integer_grid_coords = quantize_subdivided_sparse_voxel_grid(num_subdivs)

#     """

#     from ._pcu_internal import subdivide_sparse_voxel_grid_internal

#     gc = grid_coordinates
#     for _ in range(num_subdivs):
#         gc = subdivide_sparse_voxel_grid_internal(gc)

#     return gc / 2 ** num_subdivs


# def quantize_subdivided_sparse_voxel_grid(fractional_grid_coordinates, num_subdivs):
#     """
#     subdivide_sparse_voxel_grid returns fractional voxel coordinates (since we are subdividing an integer grid).
#     This function transforms these fractional coordinates into integer coordinates. This operation can be undone
#     with unquantize_subdivided_sparse_voxel_grid

#     Args:
#         fractional_grid_coordinates: An (n, 3) shaped array of fractional voxels coordinates returned by
#                                     subdivide_sparse_voxel_grid
#         num_subdivs: The number of subdivision iterations used to compute fractional_grid_coordinates

#     Returns:
#         quantized_grid_coordinates: An (n, 3) shaped array of quantized grid coordinates

#     See Also:
#         subdivide_sparse_voxel_grid
#         unquantize_subdivided_sparse_voxel_grid
#     """
#     return (fractional_grid_coordinates * (2 ** num_subdivs) - 1).astype(np.int32) // 2


# def unquantize_subdivided_sparse_voxel_grid(quantized_grid_coordinates, num_subdivs):
#     """
#     Undo quantization done by quantize_subdivided_sparse_voxel_grid.

#     Args:
#         quantized_grid_coordinates: An (n, 3) shaped array of quantized voxels coordinates returned by
#                                     quantize_subdivided_sparse_voxel_grid
#         num_subdivs: The number of subdivision iterations used to compute the original fractional coordinates

#     Returns:
#         fractional_grid_coordinates: An (n, 3) shaped array of fractional grid coordinates

#     See Also:
#         subdivide_sparse_voxel_grid
#         quantize_subdivided_sparse_voxel_grid
#     """
#     return 2.0 * quantized_grid_coordinates / (2 ** num_subdivs) + 1.0 / (2 ** num_subdivs)


# def dilate_sparse_voxel_grid(grid_coordinates, count=1):
#     """
#     Perform binary dilation on a sparse voxel grid, using a 3x3 filter of ones applied count times.
#     i.e. convolve the filter:
#         [[1, 1, 1],
#          [1, 1, 1],
#          [1, 1, 1]]
#     with the sparse voxel grid.

#     Args:
#         grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
#         count: The number of iterations of binary dilation to run

#     Returns:
#         dilated_grid_coordinates: An (m, 3) integer array encoding the coordinates of the sparse voxel grid after dilation
#     """
#     from ._pcu_internal import dilate_sparse_voxel_grid_internal
#     for _ in range(count):
#         grid_coordinates = dilate_sparse_voxel_grid_internal(grid_coordinates)

#     return grid_coordinates


# def erode_sparse_voxel_grid(grid_coordinates, count=1):
#     """
#     Perform binary erosion on a sparse voxel grid, using a 3x3 filter of ones applied count times.
#     i.e. take the negation of the convolution of the filter:
#         [[1, 1, 1],
#          [1, 1, 1],
#          [1, 1, 1]]
#     with the sparse voxel grid.

#     Args:
#         grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
#         count: The number of iterations of binary erosion to run

#     Returns:
#         eroded_grid_coordinates: An (m, 3) integer array encoding the coordinates of the sparse voxel grid after erosion
#     """
#     from ._pcu_internal import erode_sparse_voxel_grid_internal
#     for _ in range(count):
#         grid_coordinates = erode_sparse_voxel_grid_internal(grid_coordinates)

#     return grid_coordinates


# def sparse_voxel_grid_to_hex_mesh(grid_coordinates, voxel_size=1.0, origin=(0.0, 0.0, 0.0), dtype=np.float64):
#     """
#     Convert a sparse voxel grid where the voxel indices represent the centers of each voxel coordinate into a
#     hexahedral mesh.

#     i.e. a pair of arrays (v, c) where v has shape (n, 3) and each row is a vertex position
#          and c has shape (m, 8) where the i^th row c[i, :] are the indices into v of the 8 corners of the i^th cube

#     Args:
#         grid_coordinates: An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid
#         voxel_size: Either a scalar encoding the length of one voxel along an axis or a triple encoding the length of a
#                     voxel on each of 3 axes
#         origin: A triple encoding a global offset for the hex mesh (i.e. all vertices will have origin added to them).
#                 Default is (0.0, 0.0, 0.0).
#         dtype: The scalar type of vertices to return.
#             Default is np.float64.

#     Returns:
#         vertices: An (n, 3) shaped array where each row vertices[i, :] is a vertex position in the hex mesh
#         cubes: An (m, 8) shaped array of indices into v where each row cubes[j, :] are indices of the 8 points forming
#             the j^th cube
#     """
#     from ._pcu_internal import sparse_voxel_grid_to_hex_mesh_internal
#     if np.isscalar(voxel_size):
#         voxel_size = np.array([voxel_size] * 3)
#     elif len(voxel_size) == 1:
#         voxel_size = np.array([np.squeeze(voxel_size).reshape((1,))[0]] * 3)
#     elif len(voxel_size) == 3:
#         voxel_size = np.squeeze(voxel_size).reshape(3,)
#     else:
#         raise ValueError("voxel_size must be a scalar or an array/tuple/list of length 3")
#     if not np.all(np.isreal(voxel_size)):
#         raise ValueError("voxel_size must have a real numeric type")

#     if len(origin) == 3:
#         origin = np.squeeze(origin).reshape(3,)
#     else:
#         raise ValueError("origin must be an array/tuple/list of length 3")
#     if not np.all(np.isreal(origin)):
#         raise ValueError("origin must have a real numeric type")

#     ret_v, ret_c = sparse_voxel_grid_to_hex_mesh_internal(
#         grid_coordinates, voxel_size[0], voxel_size[1], voxel_size[2],
#         origin[0], origin[1], origin[2])

#     return ret_v.astype(dtype), ret_c

