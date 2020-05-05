import point_cloud_utils as pcu
import numpy as np
from mayavi import mlab
from scipy.spatial import cKDTree
k = 50

div = 1.0 / 1000.0
pts = (np.random.rand(10_000, 3) / div).astype(np.int32)
qpts = (np.random.rand(100, 3) / div).astype(np.int32)
pts[:, 2] = 0.0
qpts[:, 2] = 0.0

codes = pcu.morton_encode(pts)
codes_sorted_idx = np.argsort(codes)
codes_sorted = codes[codes_sorted_idx]

qcodes = pcu.morton_encode(qpts)

nn_idx = pcu.morton_knn(codes_sorted, qcodes, k)
nn_pts = pcu.morton_decode(codes_sorted[nn_idx[0]])

kdt = cKDTree(pts)
_, nn_gt_idx = kdt.query(np.array([qpts[0]]), k=k)
print(nn_gt_idx.shape)
nn_gt_pts = pts[nn_gt_idx[0]]

mlab.points3d(pts[:, 0], pts[:, 1], pts[:, 2], scale_factor=5.0)
mlab.points3d([qpts[0, 0]], [qpts[0, 1]], [qpts[0, 2]], scale_factor=7.0, color=(1.0, 0.0, 0.0))
mlab.points3d(nn_pts[:, 0], nn_pts[:, 1], nn_pts[:, 2], scale_factor=9.0, color=(0.0, 0.0, 1.0), opacity=0.5)
mlab.points3d(nn_gt_pts[:, 0], nn_gt_pts[:, 1], nn_gt_pts[:, 2], scale_factor=9.0, color=(0.0, 1.0, 0.0), opacity=0.5)
mlab.show()
