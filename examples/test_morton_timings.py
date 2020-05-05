import point_cloud_utils as pcu
import numpy as np
import time

start_time = time.time()

div = 1.0 / 1000.0

pts = np.random.rand(1_000_000, 3) / div
pts_int = pts.astype(np.int32)

qpts = np.random.rand(10_000, 3) / div
qpts_int = qpts.astype(np.int32)

init_time = time.time()


codes = pcu.morton_encode(pts_int, sort=False)
codes_time = time.time()

idxs = np.argsort(codes)
sort_time1 = time.time()
codes = codes[idxs]
sort_time2 = time.time()

qcodes = pcu.morton_encode(qpts_int)
qcodes_time = time.time()

nn_idx = pcu.morton_knn(codes, qcodes, k=20)
knn_time = time.time()

print(nn_idx.shape)

print("Timing results: ")
print("  Init time %f ms" % (1000 * (init_time - start_time)))
print("  Codes time %f ms" % (1000 * (codes_time - init_time)))
print("  Sort time 1 %f ms" % (1000 * (sort_time1 - codes_time)))
print("  Sort time 2 %f ms" % (1000 * (sort_time2 - sort_time1)))
print("  QCodes time %f ms" % (1000 * (qcodes_time - sort_time2)))
print("  KNN time %f ms" % (1000 * (knn_time - qcodes_time)))
