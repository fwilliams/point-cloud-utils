from ._pcu_internal import point_cloud_splatting_geometry_internal_, \
    _voxel_mesh_internal

import numpy as np


def _coord3d_to_array(coord, dtype = np.float64):
    if not hasattr(coord, "__len__") or len(coord) != 3:
        raise ValueError("expected 3D coordinate")

    if isinstance(coord, np.ndarray):
        return np.array([c.item() for c in coord], dtype=dtype)
    else:
        return np.array([c for c in coord], dtype=dtype)


def _number_or_coord3d_to_array(coord_or_number, dtype=np.float64):
    if isinstance(coord_or_number, (float, int)):
        return np.array([coord_or_number] * 3, dtype=dtype)

    return _coord3d_to_array(coord_or_number, dtype=dtype)


def _validate_point_radius_internal(p, r):
    if not isinstance(r, np.ndarray):
        if np.isscalar(r):
            r = r * np.ones(p.shape[0])
        elif isinstance(r, list) or isinstance(r, tuple):
            r = np.array(r).astype(p.dtype)
        else:
            raise ValueError("Argument r must be a scalar or numpy array with the same number of rows as p")

    if r.shape[0] != p.shape[0]:
        raise ValueError("Argument r have the same number of rows as p")

    if len(r.shape) != 1:
        r = np.squeeze(r)

    if len(r.shape) != 1:
        raise ValueError("Invalid shape for argument r, must have shape (N,) or (N, 1)")

    return r


def voxel_grid_geometry(ijk, voxel_size=np.array((1., 1., 1.)), voxel_origin=np.array((0., 0., 0.)), gap_fraction=0.0):
    """
    Generate a triangle mesh of cubes for voxel coordinates ijk. The [0, 0, 0] voxel has its
    center at voxel_origin and each voxel has voxel_size.

    Args:
        ijk np.ndarray: [num_voxels, 3] array of integer voxel coordinates
        voxel_size: Float or triple representing the size of each voxel.
                    Defaults to np.array((1., 1., 1.)).
        voxel_origin: Center coordinate of the [0, 0, 0] voxel. Defaults to np.array((0., 0., 0.)).
        gap_fraction: Fraction of a voxel to leave as a gap between voxels (default 0.0)

    Returns:
        v: Numpy array of vertices for the cube mesh
        f: Numpy array of faces for the cube mesh
    """
    return _voxel_mesh_internal(ijk, gap_fraction,
                                _coord3d_to_array(voxel_origin, dtype=np.float64),
                                _number_or_coord3d_to_array(voxel_size, dtype=np.float64))


def pointcloud_sphere_geometry(p, r, num_stacks, num_slices):
    """
    Generate sphere geometry for a point cloud (i.e. one sphere per point)

    Args:
        p : \#p by 3 array of vertex positions (each row is a vertex)
        r : Array or Scalar describing the radius along each axis (Either one radius per vertex, or a global size for the whole point cloud)
        num_stacks : Number of latitudal subdivisions
        num_slices : Number of longitudal subdivisions

    Returns:
        verts : an array of shape (#output_vertices, 3)
        faces : an array of shape (#output_faces, 3) indexing into verts

    """
    r = _validate_point_radius_internal(p, r)
    return point_cloud_splatting_geometry_internal_(p, p, "sphere", r.astype(p.dtype), num_stacks, num_slices)


def pointcloud_surfel_geometry(p, n, r=0.1, subdivs=7):
    """
    Generate geometry for a point cloud encoded as surfels (i.e. circular patches centered at each point and oriented
    perpendicularly to each normal)

    Args:
        p : \#p by 3 array of vertex positions (each row is a vertex)
        n : \#p by 3 array of vertex normals (each row is a vertex)
        r : Array or Scalar describing the size of each geometry element (Either one radius per vertex, or a global size for the whole point cloud)
        subdivs : Number of times to subdivide the patch geometry for each point (i.e. # tris per cicle)

    Returns:
        verts : an array of shape (#output_vertices, 3)
        faces : an array of shape (#output_faces, 3) indexing into verts

    """
    r = _validate_point_radius_internal(p, r)
    return point_cloud_splatting_geometry_internal_(p, n, "circle", r.astype(p.dtype), subdivs)