from __future__ import print_function
import unittest
import os
class TestDenseBindings(unittest.TestCase):
    def setUp(self):
        self.test_path = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), "..", "data")

    def load_mesh(self):
        import numpy as np
        # v is a nv by 3 NumPy array of vertices
        v_path = os.path.join(self.test_path, "1_vertices.npy")
        # f is an nf by 3 NumPy array of face indexes into v
        f_path = os.path.join(self.test_path, "1_faces.npy")
        # n is a nv by 3 NumPy array of vertex normals if they are specified, otherwise an empty array
        n_path = os.path.join(self.test_path, "1_vertex_normals.npy")
        nf_path = os.path.join(self.test_path, "1_face_normals.npy")
        return np.load(v_path), np.load(f_path), np.load(n_path), np.load(nf_path)
        
    def test_trianglemesh_from_arrays(self):
        import point_cloud_utils as pcu
        v,f,n,fn = self.load_mesh()
        face_data = pcu.TriangleMesh.FaceData()
        face_data.normals = fn
        face_data.vertex_ids = f
        vertex_data = pcu.TriangleMesh.VertexData()
        vertex_data.positions = v
        vertex_data.normals = n
        mesh = pcu.TriangleMesh()
        mesh.face_data = face_data
        mesh.vertex_data = vertex_data
        

    def test_lloyd_relaxation(self):
        import point_cloud_utils as pcu

        # v is a nv by 3 NumPy array of vertices
        # f is an nf by 3 NumPy array of face indexes into v
        v, f, n, fn = self.load_mesh()

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

        chamfer_dist = pcu.chamfer_distance(a, b)
        chamfer_dist, c_a_to_b, c_b_to_a = pcu.chamfer_distance(a, b, return_index=True)

    def test_knn(self):
        import point_cloud_utils as pcu
        import numpy as np

        for i in range(3):
            # Generate two random point sets
            a = np.random.rand(1000, 3)
            b = np.random.rand(500, 3)
            # dists_a_to_b is of shape (a.shape[0],) and contains the shortest squared distance
            # between each point in a and the points in b
            # corrs_a_to_b is of shape (a.shape[0],) and contains the index into b of the
            # closest point for each point in a
            k = np.random.randint(10) + 1
            dists_a_to_b, corrs_a_to_b = pcu.k_nearest_neighbors(a, b, k)
            if k > 1:
                self.assertEqual(dists_a_to_b.shape, (a.shape[0], k))
                self.assertEqual(corrs_a_to_b.shape, (a.shape[0], k))
            else:
                self.assertEqual(dists_a_to_b.shape, (a.shape[0],))
                self.assertEqual(corrs_a_to_b.shape, (a.shape[0],))

            if k == 1:
                dists_a_to_b = dists_a_to_b[:, np.newaxis]
                corrs_a_to_b = corrs_a_to_b[:, np.newaxis]

            for i in range(dists_a_to_b.shape[1]):
                b_map = b[corrs_a_to_b[:, i]]
                dists = np.linalg.norm(a - b_map, axis=-1)
                diff_dists = dists - dists_a_to_b[:, i]
                self.assertTrue(np.all(np.abs(diff_dists) < 1e-5))

            b_map = b[corrs_a_to_b]
            dists = np.linalg.norm(a[:, np.newaxis, :] - b_map, axis=-1)
            diff_dists = dists - dists_a_to_b
            self.assertTrue(np.all(np.abs(diff_dists) < 1e-5))

        with self.assertRaises(ValueError):
            a = np.random.rand(1000, 3)
            b = np.random.rand(500, 3)
            dists_a_to_b, corrs_a_to_b = pcu.k_nearest_neighbors(a, b, 0)

        a = np.random.rand(100, 3)
        b = np.random.rand(50, 3)
        dists_a_to_b, corrs_a_to_b = pcu.k_nearest_neighbors(a, b, 3)
        dists_a_to_b2, corrs_a_to_b2 = pcu.k_nearest_neighbors(a, b, 3, squared_distances=True)

        self.assertTrue(np.all(corrs_a_to_b == corrs_a_to_b2))
        self.assertTrue(np.all(np.abs(dists_a_to_b ** 2.0 - dists_a_to_b2) < 1e-5))

    def test_hausdorff(self):
        import point_cloud_utils as pcu
        import numpy as np

        # Generate two random point sets
        a = np.random.rand(1000, 3)
        b = np.random.rand(500, 3)

        # Compute each one sided squared Hausdorff distances
        hausdorff_a_to_b, idx_a1, idx_b1 = pcu.one_sided_hausdorff_distance(a, b, return_index=True)
        hausdorff_b_to_a, idx_b2, idx_a2 = pcu.one_sided_hausdorff_distance(b, a, return_index=True)

        # Take a max of the one sided squared  distances to get the two sided Hausdorff distance
        hausdorff_a_b = max(hausdorff_a_to_b, hausdorff_b_to_a)
        hausdorff_a_b_pcu, i1, i2 = pcu.hausdorff_distance(a, b, return_index=True)
        self.assertAlmostEqual(hausdorff_a_b, hausdorff_a_b_pcu)
        self.assertAlmostEqual(hausdorff_a_b, np.linalg.norm(a[i1] - b[i2]))
        self.assertAlmostEqual(hausdorff_a_b_pcu, np.linalg.norm(a[i1] - b[i2]))
        if hausdorff_a_to_b > hausdorff_b_to_a:
            self.assertEqual(i1, idx_a1)
            self.assertEqual(i2, idx_b1)
        else:
            self.assertEqual(i1, idx_a2)
            self.assertEqual(i2, idx_b2)

        # Find the index pairs of the two points with maximum shortest distancce
        hausdorff_b_to_a, idx_b, idx_a = pcu.one_sided_hausdorff_distance(b, a, return_index=True)
        self.assertAlmostEqual(np.linalg.norm(a[idx_a] - b[idx_b]), hausdorff_b_to_a)

    def test_shortest_dist_two_point_clouts(self):
        import point_cloud_utils as pcu
        import numpy as np

        # Generate two random point sets
        a = np.random.rand(1000, 3)
        # Comparing identical pointset, expecting all 0
        b = a.copy()

        # dists_a_to_b is of shape (a.shape[0],) and contains the shortest squared distance
        # between each point in a and the points in b
        # corrs_a_to_b is of shape (a.shape[0],) and contains the index into b of the
        # closest point for each point in a
        dists_a_to_b, corrs_a_to_b = pcu.k_nearest_neighbors(a, b, k=1)
        self.assertEqual(np.sum(dists_a_to_b), 0.)

    def skip_test_morton_coding_big_data(self):
        import point_cloud_utils as pcu
        import numpy as np
        import os
        if os.name == 'nt':
            num_pts = 1000
            num_qpts = 10
        else:
            num_pts = 1000000
            num_qpts = 10000
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

    def skip_test_morton_coding_small_data(self):
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

    def skip_test_morton_coding_tiny_data(self):
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

    def test_remove_duplicate_points(self):
        import point_cloud_utils as pcu
        import numpy as np

        # v is a nv by 3 NumPy array of vertices
        v = np.load(os.path.join(self.test_path, "duplicate_vertices.npy"))

        v2, idx_v_to_v2, idx_v2_to_v = pcu.deduplicate_point_cloud(v, 1e-11, return_index=True)
        self.assertLess(v2.shape[0], v.shape[0])
        self.assertTrue(np.all(np.equal(v2[idx_v2_to_v], v)))
        self.assertTrue(np.all(np.equal(v[idx_v_to_v2], v2)))

    def test_closest_points_on_mesh(self):
        import point_cloud_utils as pcu
        import numpy as np

        v, f, _, _ = self.load_mesh()

        v -= v.min(0)
        v /= v.max(0)
        v -= 0.5
        v *= 2.0

        p = np.random.rand(1000, 3)

        d, fi, bc = pcu.closest_points_on_mesh(p, v, f)

        closest_points = (v[f[fi]] * bc[:, :, np.newaxis]).sum(1)
        d2 = np.linalg.norm(p - closest_points, axis=-1)

        self.assertAlmostEqual(np.max(d - d2), 0.0)

    def test_signed_distance_to_mesh(self):
        import point_cloud_utils as pcu
        import numpy as np

        v, f, _, _ = self.load_mesh()
        # compare against itself
        d,fi,b = pcu.signed_distance_to_mesh(v,v,f)
        self.assertEqual(np.sum(d),0.)

    def test_manifold(self):
        import point_cloud_utils as pcu
        v, f, _, _ = self.load_mesh()
        resolution = 20_000  
        v_watertight, f_watertight = pcu.make_mesh_watertight(v, f, resolution=resolution)
        self.assertEqual(v_watertight.shape[0], 64748)
        self.assertEqual(f_watertight.shape[0], 129494)

    def test_ray_mesh_intersection(self):
        import point_cloud_utils as pcu
        import numpy as np

        v, f, _, _ = self.load_mesh()

        v -= v.min(0)
        v /= v.max(0)
        v -= 0.5
        v *= 2.0

        d = np.concatenate([
                np.stack([a.ravel() for a in np.mgrid[-0.1:0.1:64j, -0.1:0.1:64j]], axis=-1),
                0.1 * np.ones([64 ** 2, 1])],
            axis=-1)
        d /= np.linalg.norm(d, axis=-1, keepdims=True)

        o1 = np.array([0., 0., -2.])
        fid1, bc1, t1 = pcu.ray_mesh_intersection(v, f, o1, d)
        mask1 = np.isfinite(t1)
        self.assertTrue(mask1.sum() > 0)

        p11 = pcu.interpolate_barycentric_coords(f, fid1[mask1], bc1[mask1], v)
        p12 = o1 + t1[mask1, np.newaxis] * d[mask1]
        self.assertTrue(np.allclose(p11, p12, atol=1e-5))

        o2 = np.stack([o1]*d.shape[0])
        fid2, bc2, t2 = pcu.ray_mesh_intersection(v, f, o2, d)
        mask2 = np.isfinite(t2)
        self.assertTrue(mask2.sum() > 0)

        p21 = pcu.interpolate_barycentric_coords(f, fid2[mask2], bc2[mask2], v)
        p22 = o2[mask2] + t2[mask2, np.newaxis] * d[mask2]
        self.assertTrue(np.allclose(p21, p22, atol=1e-5))

        self.assertTrue(np.alltrue(mask1 == mask2))
        self.assertTrue(np.alltrue(fid2 == fid1))
        self.assertTrue(np.allclose(bc2, bc1))
        self.assertTrue(np.allclose(t1, t2))


if __name__ == '__main__':
    unittest.main()
