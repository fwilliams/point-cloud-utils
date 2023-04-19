import numpy as np

from ._pcu_internal import _flood_fill_3d_internal, _voxelize_triangle_mesh_internal
from ._point_cloud_geometry import _coord3d_to_array, _number_or_coord3d_to_array
from . import morton_encode, morton_decode

def flood_fill_3d(grid, start_coord, fill_value):
    """
    Flood fill a 3D grid starting from_start_coord with fill_value. This will return a
    copy of grid where the region in the input grid which is connected to start_coord
    and shares the value of start_coord is set to fill_value.

    Args:
        grid (np.ndarray): [w, h, d] array of scalars as input to the flood fill
        start_coord: (i, j, k) integer coordinate to start the flood fill
        fill_value: scalar value to flood fill

    Returns:
        _type_: A flood filled copy of grid where all voxels which are connected to start_coord
                are set to fill_value
    """
    start_coord = _coord3d_to_array(start_coord, np.int32)
    sizes = grid.shape
    if len(sizes) != 3:
        raise ValueError("grid must have shape [w, h, d]")

    ret = np.ascontiguousarray(grid).ravel().copy()
    _flood_fill_3d_internal(ret, start_coord[0], start_coord[1], start_coord[2],
                            sizes[0], sizes[1], sizes[2], float(fill_value))
    return ret.reshape(sizes)


def voxelize_triangle_mesh(v, f, voxel_size, voxel_origin):
    """
    Return ijk coordinates of voxels which intersect the given mesh.
    Each voxel is assumed to have size voxel_size (scalar or triple of floats) and
    the (0, 0, 0) voxel has its bottom-back-left corner at voxel_origin

    Args:
        v (np.ndarray): [num_vertices, 3] array of triangle mesh vertices
        f (np.ndarray): [num_faces, 3] array of face indexes into v
        voxel_size: A float or triple specifying the size of each voxel
        voxel_origin: A float or triple specifying the the position of the
                      bottom-back-left corner of the (0, 0, 0) voxel

    Returns:
        ijk: [num_vox, 3] array of integer ijk coordinates for each voxel intersecting the mesh
    """
    voxel_origin = _coord3d_to_array(voxel_origin, np.float64)
    voxel_size = _number_or_coord3d_to_array(voxel_size, np.float64)

    ijk = _voxelize_triangle_mesh_internal(v, f, voxel_origin, voxel_size)
    ijk = morton_decode(np.unique(morton_encode(ijk)))
    return ijk

