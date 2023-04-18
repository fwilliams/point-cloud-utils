import numpy as np

from ._point_cloud_geometry import _validate_point_radius_internal


def ray_surfel_intersection(p, n, ray_o, ray_d, r=0.1, subdivs=4, ray_near=0.0, ray_far=np.inf):
    """
    Compute intersection between a set of rays and a point cloud converted to surfels (i.e. circular patches oriented
    with the point normals)

    Args:
        p : \\#p by 3 array of vertex positions (each row is a vertex)
        n : \\#p by 3 Matrix of vertex normals (each row is a vertex)
        ray_o : array of shape (#rays, 3) of ray origins (one per row) or a single array of shape (3,) to use
        ray_d : array of shape (#rays, 3) of ray directions (one per row)
        r : Array or Scalar describing the size of each geometry element (Either one radius per vertex, or a global size for the whole point cloud)
        subdivs : Number of times to subdivide the patch geometry for each point (i.e. # tris per cicle)
        ray_near : an optional floating point value indicating the distance along each ray to start searching (default 0.0)
        ray_far : an optional floating point value indicating the maximum distance along each ray to search (default inf)

    Returns:
        t : a (#rays,) shaped array encoding the distance between the ray origin and intersection point for each ray (inf for missed rays)
        pid : a (#rays,) shaped array of integer indices corresponding to which points were hit (-1 for a ray miss)
    """
    from ._pcu_internal import ray_point_cloud_intersection_internal_
    r = _validate_point_radius_internal(p, r)
    return ray_point_cloud_intersection_internal_(p, n, ray_o, ray_d, "circle", r.astype(p.dtype), subdivs, subdivs, ray_near, ray_far)


class RaySurfelIntersector:
    """
    Class used to find the intersection between rays and a point cloud converted to surfels (i.e. a point cloud
    represented as circles centered at each point and oriented perpendicularly to each normal).
    """
    def __init__(self, p, n, r=0.1, subdivs=7):
        """
        Construct a RayPointCloudIntersector which can be used to compute the intersection between a set of rays and a
        point cloud converted to circular patches oriented with the point normals

        Args:
            p : \\#p by 3 array of vertex positions (each row is a vertex)
            n : \\#p by 3 Matrix of vertex normals (each row is a vertex)
            r : Array or Scalar describing the size of each geometry element (Either one radius per vertex, or a global size for the whole point cloud)
            subdivs : Number of times to subdivide the patch geometry for each point (i.e. # tris per cicle)
        """
        from ._pcu_internal import _RayMeshIntersectorInternal, _populate_ray_point_cloud_intersector_internal
        r = _validate_point_radius_internal(p, r)

        self.p = p
        self.n = n
        self.r = r
        self.num_subdivs = subdivs
        self.__internal_intersector = _RayMeshIntersectorInternal()
        self.__num_faces_per_geom = _populate_ray_point_cloud_intersector_internal(p, n, "circle", r, subdivs, subdivs,
                                                                                   self.__internal_intersector)

    def intersect_rays(self, ray_o, ray_d, ray_near=0.0, ray_far=np.inf):
        """
        Compute intersection between a set of rays and the point cloud converted to circular patches oriented with the
        point normals

        Args:
            ray_o : array of shape (#rays, 3) of ray origins (one per row) or a single array of shape (3,) to use
            ray_d : array of shape (#rays, 3) of ray directions (one per row)
            ray_near : an optional floating point value indicating the distance along each ray to start searching (default 0.0)
            ray_far : an optional floating point value indicating the maximum distance along each ray to search (default inf)

        Returns:
            t : a (#rays,) shaped array encoding the distance between the ray origin and intersection point for each ray (inf for missed rays)
            pid : a (#rays,) shaped array of integer indices corresponding to which points were hit (-1 for a ray miss)
        """
        from ._pcu_internal import _intersect_ray_point_cloud_intersector_internal
        return _intersect_ray_point_cloud_intersector_internal(ray_o, ray_d, self.__internal_intersector,
                                                               self.__num_faces_per_geom, ray_near, ray_far)
