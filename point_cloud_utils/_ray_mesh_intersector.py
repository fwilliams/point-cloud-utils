import numpy as np


class RayMeshIntersector:
    """
    Class used to find the intersection between rays and a triangle mesh.
    """
    def __init__(self, mesh_v, mesh_f):
        """
        Create a RayMeshIntersector object which can be used to do ray/mesh queries with a triangle mesh.

        Parameters
        ----------
        mesh_v : #v by 3 array of vertex positions (each row is a vertex)
        mesh_f: #f by 3 Matrix of face (triangle) indices
        """
        from ._pcu_internal import _RayMeshIntersectorInternal, _populate_ray_intersector_internal
        self.__internal_intersector = _RayMeshIntersectorInternal()
        self.v = mesh_v
        self.f = mesh_f
        _populate_ray_intersector_internal(mesh_v, mesh_f, self.__internal_intersector)

    def intersect_rays(self, ray_o, ray_d, ray_near=0.0, ray_far=np.inf):
        """
        Compute intersection between a set of rays and the triangle mesh enclosed in this class

        Parameters
        ----------
        ray_o : array of shape (#rays, 3) of ray origins (one per row) or a single array of shape (3,) to use
        ray_d : array of shape (#rays, 3) of ray directions (one per row)
        ray_near : an optional floating point value indicating the distance along each ray to start searching (default 0.0)
        ray_far : an optional floating point value indicating the maximum distance along each ray to search (default inf)

        Returns
        -------
        A tuple (f_id, bc, t) where:
          - f_id is an array of shape (#rays,) representing the face id hit by each ray
          - bc is an array of shape (#rays, 3) where each row is the barycentric coordinates
            within each face of the ray intersection
          - t is the distance along each ray to the intersection
        """
        from ._pcu_internal import _intersect_ray_intersector_internal
        return _intersect_ray_intersector_internal(ray_o, ray_d, self.__internal_intersector, ray_near, ray_far)