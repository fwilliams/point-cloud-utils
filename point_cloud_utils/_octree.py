import numpy as np


class Octree:
    def __init__(self, max_depth, origin=np.array((0.0, 0.0, 0.0)), size=0.0):
        from ._pcu_internal import Octree
        self.__internal_octree = Octree(max_depth, origin[0], origin[1], origin[2], size)

    def clear(self):
        self.__internal_octree.clear()

    def is_empty(self):
        return self.__internal_octree.is_empty()

    @property
    def min_bound(self):
        return np.array(self.__internal_octree.get_min_bound())

    @property
    def max_bound(self):
        return np.array(self.__internal_octree.get_max_bound())

    @property
    def center(self):
        return np.array(self.__internal_octree.get_center())

    def build_from_point_cloud(self, points, pad_amount=0.01):
        points = self._check_shape(points)
        from ._pcu_internal import build_octree_from_pointcloud_internal
        build_octree_from_pointcloud_internal(self.__internal_octree, points, pad_amount)

    def find(self, points):
        points = self._check_shape(points)
        from ._pcu_internal import get_octree_point_leaves_internal
        return get_octree_point_leaves_internal(self.__internal_octree, points)

    def point_depths(self, points):
        points = self._check_shape(points)
        from ._pcu_internal import get_octree_point_depths_internal
        return get_octree_point_depths_internal(self.__internal_octree, points)

    @staticmethod
    def _check_shape(points):
        if not isinstance(points, np.ndarray):
            raise ValueError("points must be a numpy array of shape (N, 3) but got points of "
                             "type %s" % str(type(points)))
        if len(points.shape) == 1:
            points = points[np.newaxis, :]
        if len(points.shape) != 2:
            raise ValueError("Invalid input points must have shape (N, 3), but got %s" % str(points.shape))
        if points.shape[0] <= 0:
            raise ValueError("Invalid input points must have greater than zero points")
        if points.shape[1] != 3:
            raise ValueError("Invalid input points must have shape (N, 3), but got %s" % str(points.shape))
        return points
