import numpy as np


class RayMeshIntersector:
    def __init__(self, mesh_v, mesh_f):
        from ._pcu_internal import _RayMeshIntersectorInternal, _populate_ray_intersector_internal
        self.__internal_intersector = _RayMeshIntersectorInternal()
        self.v = mesh_v
        self.f = mesh_f
        _populate_ray_intersector_internal(mesh_v, mesh_f, self.__internal_intersector)

    def intersect_rays(self, ray_o, ray_d, ray_near=0.0, ray_far=np.inf):
        from ._pcu_internal import _intersect_ray_intersector_internal
        return _intersect_ray_intersector_internal(ray_o, ray_d, self.__internal_intersector, ray_near, ray_far)