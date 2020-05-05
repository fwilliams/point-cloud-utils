import point_cloud_utils as pcu
import numpy as np

# Lots of data
div = 1.0 / 1000.0
pts = np.random.rand(1_000_000, 3) / div
pts_int = pts.astype(np.int32)

qpts = np.random.rand(10_000, 3) / div
qpts_int = qpts.astype(np.int32)

codes = pcu.morton_encode(pts_int)
codes_sort_idx = np.argsort(codes)
codes_sorted = codes[codes_sort_idx]

qcodes = pcu.morton_encode(qpts_int)

nn_idx = pcu.morton_knn(codes_sorted, qcodes, 7)
codes_sorted[nn_idx]
print("0 --------------------------------")
print(nn_idx.shape)
print(nn_idx)
print("\n")

# 1nn
div = 1.0 / 1000.0
pts = np.random.rand(1_000_000, 3) / div
pts_int = pts.astype(np.int32)

qpts = np.random.rand(10_000, 3) / div
qpts_int = qpts.astype(np.int32)

codes = pcu.morton_encode(pts_int)
codes_sort_idx = np.argsort(codes)
codes_sorted = codes[codes_sort_idx]

qcodes = pcu.morton_encode(qpts_int)

nn_idx = pcu.morton_knn(codes_sorted, qcodes, 1)
print("1 --------------------------------")
print(nn_idx.shape)
print(nn_idx)
print("\n")


# Tiny amount of data to hit the boundary cases
div = 1.0 / 1000.0
pts = np.random.rand(10, 3) / div
pts_int = pts.astype(np.int32)

qpts = np.random.rand(10_000, 3) / div
qpts_int = qpts.astype(np.int32)

codes = pcu.morton_encode(pts_int)
codes_sort_idx = np.argsort(codes)
codes_sorted = codes[codes_sort_idx]

qcodes = pcu.morton_encode(qpts_int)

nn_idx = pcu.morton_knn(codes_sorted, qcodes, 7)
print("2 --------------------------------")
print(nn_idx.shape)
print(nn_idx)
print("\n")


# k > data
div = 1.0 / 1000.0
pts = np.random.rand(10, 3) / div
pts_int = pts.astype(np.int32)

qpts = np.random.rand(10_000, 3) / div
qpts_int = qpts.astype(np.int32)

codes = pcu.morton_encode(pts_int)
codes_sort_idx = np.argsort(codes)
codes_sorted = codes[codes_sort_idx]

qcodes = pcu.morton_encode(qpts_int)

nn_idx = pcu.morton_knn(codes_sorted, qcodes, 15)
print("3 --------------------------------")
print(nn_idx.shape)
print(nn_idx)
