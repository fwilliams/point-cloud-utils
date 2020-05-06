from __future__ import print_function
import unittest
import os


class TestDenseBindings(unittest.TestCase):
    def setUp(self):
        self.test_path = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), "..", "data")

    def test_poisson_disk_sampling(self):
        import point_cloud_utils as pcu
        import numpy as np

        # v is a nv by 3 NumPy array of vertices
        # f is an nf by 3 NumPy array of face indexes into v
        # n is a nv by 3 NumPy array of vertex normals if they are specified, otherwise an empty array
        v, f, n = pcu.read_obj(os.path.join(self.test_path, "cube_twist.obj"))
        bbox = np.max(v, axis=0) - np.min(v, axis=0)
        bbox_diag = np.linalg.norm(bbox)

        # Generate very dense  random samples on the mesh (v, f, n)
        # Note that this function works with no normals, just pass in an empty array np.array([], dtype=v.dtype)
        # v_dense is an array with shape (100*v.shape[0], 3) where each row is a point on the mesh (v, f)
        # n_dense is an array with shape (100*v.shape[0], 3) where each row is a the normal of a point in v_dense
        v_dense, n_dense = pcu.sample_mesh_random(v, f, n, num_samples=v.shape[0] * 100)

        # Downsample v_dense to be from a blue noise distribution:
        #
        # v_poisson is a downsampled version of v where points are separated by approximately
        # `radius` distance, use_geodesic_distance indicates that the distance should be measured on the mesh.
        #
        # n_poisson are the corresponding normals of v_poisson
        v_poisson, n_poisson = pcu.prune_point_cloud_poisson_disk(v_dense, n_dense, 0.1*bbox_diag)

        v_poisson, n_poisson = pcu.sample_mesh_poisson_disk(
            v, f, n, num_samples=7777, use_geodesic_distance=True)

        v_poisson, n_poisson = pcu.sample_mesh_poisson_disk(
            v, f, n, num_samples=-1, radius=0.01*bbox_diag, use_geodesic_distance=True)



    def test_lloyd_relaxation(self):
        import point_cloud_utils as pcu

        # v is a nv by 3 NumPy array of vertices
        # f is an nf by 3 NumPy array of face indexes into v
        v, f, n = pcu.read_obj(os.path.join(self.test_path, "cube_twist.obj"))

        # Generate 1000 points on the mesh with Lloyd's algorithm
        samples = pcu.sample_mesh_lloyd(v, f, 1000)

        # Generate 100 points on the unit square with Lloyd's algorithm
        samples_2d = pcu.lloyd_2d(100)

        # Generate 100 points on the unit cube with Lloyd's algorithm
        samples_3d = pcu.lloyd_3d(100)

    def test_sinkhorn(self):
        import point_cloud_utils as pcu
        import numpy as np

        # a and b are arrays where each row contains a point
        # Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
        a = np.random.rand(100, 3)
        b = np.random.rand(100, 3)

        # M is a 100x100 array where each entry  (i, j) is the squared distance between point a[i, :] and b[j, :]
        M = pcu.pairwise_distances(a, b)

        # w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
        w_a = np.ones(a.shape[0])
        w_b = np.ones(b.shape[0])

        # P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to
        # better approximation of the true Wasserstein distance at the expense of slower convergence
        P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

        # To get the distance as a number just compute the frobenius inner product <M, P>
        sinkhorn_dist = (M * P).sum()

    def test_batched_sinkhorn(self):
        import point_cloud_utils as pcu
        import numpy as np

        # a and b are each contain 10 batches each of which contain 100 points  of dimension 3
        # i.e. a[i, :, :] is the i^th point set which contains 100 points
        # Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
        a = np.random.rand(10, 100, 3)
        b = np.random.rand(10, 100, 3)

        # M is a 10x100x100 array where each entry (k, i, j) is the squared distance between point a[k, i, :]
        # and b[k, j, :]
        M = pcu.pairwise_distances(a, b)

        # w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
        w_a = np.ones(a.shape[:2])
        w_b = np.ones(b.shape[:2])

        # P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to
        # better approximation of the true Wasserstein distance at the expense of slower convergence
        P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

        # To get the distance as a number just compute the frobenius inner product <M, P>
        sinkhorn_dist = (M * P).sum()

    def test_chamfer(self):
        import point_cloud_utils as pcu
        import numpy as np

        # a and b are arrays where each row contains a point
        # Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
        a = np.random.rand(100, 3)
        b = np.random.rand(100, 3)

        chamfer_dist = pcu.chamfer(a, b)

    def test_batched_chamfer(self):
        import point_cloud_utils as pcu
        import numpy as np

        # a and b are each contain 10 batches each of which contain 100 points  of dimension 3
        # i.e. a[i, :, :] is the i^th point set which contains 100 points
        # Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
        a = np.random.rand(10, 100, 3)
        b = np.random.rand(10, 100, 3)

        chamfer_dist = pcu.chamfer(a, b)

    def test_hausdorff_and_nearest_neighbor(self):
        import point_cloud_utils as pcu
        import numpy as np

        # Generate two random point sets
        a = np.random.rand(1000, 3)
        b = np.random.rand(500, 3)

        # dists_a_to_b is of shape (a.shape[0],) and contains the shortest squared distance
        # between each point in a and the points in b
        # corrs_a_to_b is of shape (a.shape[0],) and contains the index into b of the
        # closest point for each point in a
        dists_a_to_b, corrs_a_to_b = pcu.point_cloud_distance(a, b)

        # Compute each one sided squared Hausdorff distances
        hausdorff_a_to_b = pcu.hausdorff(a, b)
        hausdorff_b_to_a = pcu.hausdorff(b, a)

        # Take a max of the one sided squared  distances to get the two sided Hausdorff distance
        hausdorff_dist = max(hausdorff_a_to_b, hausdorff_b_to_a)

        # Find the index pairs of the two points with maximum shortest distancce
        hausdorff_b_to_a, idx_b, idx_a = pcu.hausdorff(b, a, return_index=True)
        self.assertAlmostEqual(np.sum((a[idx_a] - b[idx_b])**2), hausdorff_b_to_a)

    def test_estimate_normals(self):
        import point_cloud_utils as pcu
        import numpy as np

        # v is a nv by 3 NumPy array of vertices
        # f is an nf by 3 NumPy array of face indexes into v
        # n is a nv by 3 NumPy array of vertex normals if they are specified, otherwise an empty array
        v, f, n = pcu.read_obj(os.path.join(self.test_path, "cube_twist.obj"))

        # Estimate normals for the point set, v using 12 nearest neighbors per point
        n = pcu.estimate_normals(v, k=12)
        self.assertEqual(n.shape, v.shape)

    def test_morton_coding_big_data(self):
        import point_cloud_utils as pcu
        import numpy as np
        import os
        if os.name == 'nt':
            num_pts = 10000
            num_qpts = 100
        div = 1.0 / 1000.0
        pts = np.random.rand(num_pts, 3) / div
        pts_int = pts.astype(np.int32)

        qpts = np.random.rand(num_qpts, 3) / div
        qpts_int = qpts.astype(np.int32)

        codes = pcu.morton_encode(pts_int)
        codes_sort_idx = np.argsort(codes)
        codes_sorted = codes[codes_sort_idx]

        qcodes = pcu.morton_encode(qpts_int)

        nn_idx = pcu.morton_knn(codes_sorted, qcodes, 7)
        codes_sorted[nn_idx]
        self.assertEqual(nn_idx.shape, (num_qpts, 7))

    def test_morton_coding_small_data(self):
        import point_cloud_utils as pcu
        import numpy as np
        div = 1.0 / 1000.0
        pts = np.random.rand(10, 3) / div
        pts_int = pts.astype(np.int32)

        qpts = np.random.rand(10000, 3) / div
        qpts_int = qpts.astype(np.int32)

        codes = pcu.morton_encode(pts_int)
        codes_sort_idx = np.argsort(codes)
        codes_sorted = codes[codes_sort_idx]

        qcodes = pcu.morton_encode(qpts_int)

        nn_idx = pcu.morton_knn(codes_sorted, qcodes, 7)
        codes_sorted[nn_idx]
        self.assertEqual(nn_idx.shape, (10000, 7))

    def test_morton_coding_tiny_data(self):
        import point_cloud_utils as pcu
        import numpy as np
        div = 1.0 / 1000.0
        pts = np.random.rand(10, 3) / div
        pts_int = pts.astype(np.int32)

        qpts = np.random.rand(10000, 3) / div
        qpts_int = qpts.astype(np.int32)

        codes = pcu.morton_encode(pts_int)
        codes_sort_idx = np.argsort(codes)
        codes_sorted = codes[codes_sort_idx]

        qcodes = pcu.morton_encode(qpts_int)

        nn_idx = pcu.morton_knn(codes_sorted, qcodes, 15)
        codes_sorted[nn_idx]
        self.assertEqual(nn_idx.shape, (10000, 10))

if __name__ == '__main__':
    unittest.main()
