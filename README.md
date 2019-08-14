# Point Cloud Utils (pcu) - A Python library for common tasks on 3D point clouds

**Point Cloud Utils (pcu)** is a utility library providing the following functionality:
 - A series of algorithms for generating point samples on meshes:
   - Poisson-Disk-Sampling of a mesh based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface](http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
   - Sampling a mesh with [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
   - Sampling a mesh uniformly
 - Clustering point-cloud vertices into bins
 - Very fast pairwise nearest neighbor between point clouds (based on [nanoflann](https://github.com/jlblancoc/nanoflann))
 - Hausdorff distances between point-clouds.
 - Chamfer distnaces between point-clouds.
 - Approximate Wasserstein distances between point-clouds using the [Sinkhorn](https://arxiv.org/abs/1306.0895) method.
 - Pairwise distances between point clouds
 - Utility functions for reading and writing common mesh formats (OBJ, OFF, PLY)
 
![Example of Poisson Disk Sampling](/img/blue_noise.png?raw=true "Example of Poisson Disk Sampling")

# Installation Instructions
### With `conda` (recommended)
Simply run:
```
conda install -c conda-forge point_cloud_utils
```

### With `pip`
```
pip install git+git://github.com/fwilliams/point-cloud-utils
```
The following dependencies are required to install with `pip`:
* A C++ compiler supporting C++14 or later
* CMake 3.2 or later.
* git

# Examples

### Poisson-Disk-Sampling
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
# n is a nv by 3 NumPy array of vertex normals
v, f, n, _ = pcu.read_ply("my_model.ply")
bbox = np.max(v, axis=0) - np.min(v, axis=0)
bbox_diag = np.linalg.norm(bbox)

# Generate very dense  random samples on the mesh (v, f, n)
# Note that this function works with no normals, just pass in an empty array np.array([], dtype=v.dtype)
# v_dense is an array with shape (100*v.shape[0], 3) where each row is a point on the mesh (v, f)
# n_dense is an array with shape (100*v.shape[0], 3) where each row is a the normal of a point in v_dense
v_dense, n_dense = pcu.sample_mesh_random(v, f, n, num_samples=v.shape[0]*100)

# Downsample v_dense to be from a blue noise distribution: 
#
# v_poisson is a downsampled version of v where points are separated by approximately 
# `radius` distance, use_geodesic_distance indicates that the distance should be measured on the mesh.
#
# n_poisson are the corresponding normals of v_poisson
v_poisson, n_poisson = pcu.sample_mesh_poisson_disk(
    v_dense, f, n_dense, radius=0.01*bbox_diag, use_geodesic_distance=True)
```

### Lloyd Relaxation
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
v, f, _, _ = pcu.read_ply("my_model.ply")

# Generate 1000 points on the mesh with Lloyd's algorithm
samples = pcu.sample_mesh_lloyd(v, f, 1000)

# Generate 100 points on the unit square with Lloyd's algorithm
samples_2d = pcu.lloyd_2d(100)
```

### Approximate Wasserstein (Sinkhorn) Distance Between Point-Clouds

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are arrays where each row contains a point 
# Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
a = np.random.rand(100, 3)
b = np.random.rand(100, 3)

# M is a 100x100 array where each entry  (i, j) is the squared distance between point a[i, :] and b[j, :]
M = pcu.pairwise_distances(a, b)

# w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
w_a = np.ones_like(a)
w_b = np.ones_like(b)

# P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to 
# better approximation of the true Wasserstein distance at the expense of slower convergence
P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

# To get the distance as a number just compute the frobenius inner product <M, P>
sinkhorn_dist = (M*P).sum() 
```


##### Batched Version:

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are each contain 10 batches each of which contain 100 points  of dimension 3
# i.e. a[i, :, :] is the i^th point set which contains 100 points 
# Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
a = np.random.rand(10, 100, 3)
b = np.random.rand(10, 100, 3)

# M is a 10x100x100 array where each entry (k, i, j) is the squared distance between point a[k, i, :] and b[k, j, :]
M = pcu.pairwise_distances(a, b)

# w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
w_a = np.ones_like(a)
w_b = np.ones_like(b)

# P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to 
# better approximation of the true Wasserstein distance at the expense of slower convergence
P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

# To get the distance as a number just compute the frobenius inner product <M, P>
sinkhorn_dist = (M*P).sum() 
```


### Chamfer Distance Between Point-Clouds
```python
import point_cloud_utils as pcu
import numpy as np

# a and b are arrays where each row contains a point 
# Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
a = np.random.rand(100, 3)
b = np.random.rand(100, 3)

chamfer_dist = pcu.chamfer(a, b)
```
##### Batched Version:

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are each contain 10 batches each of which contain 100 points  of dimension 3
# i.e. a[i, :, :] is the i^th point set which contains 100 points 
# Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
a = np.random.rand(10, 100, 3)
b = np.random.rand(10, 100, 3)

chamfer_dist = pcu.chamfer(a, b)
```


### Nearest-Neighbors and Hausdorff Distances Between Point-Clouds
```python
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
hausdorff_dist = np.max(hausdorff_a_to_b, hausdorff_b_to_a)

# Find the index pairs of the two points with maximum shortest distancce
hausdorff_b_to_a, idx_b, idx_a = pcu.hausdorff(b, a, return_index=True)
assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_b_to_a) < 1e-5, "These values should be almost equal"
```

