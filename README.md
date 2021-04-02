# Point Cloud Utils (pcu) - A Python library for common tasks on 3D point clouds

[![Build Status](https://travis-ci.org/fwilliams/point-cloud-utils.svg?branch=master)](https://travis-ci.org/fwilliams/point-cloud-utils)
[![Build status](https://ci.appveyor.com/api/projects/status/ujv44lqbeosgl9ij/branch/master?svg=true)](https://ci.appveyor.com/project/fwilliams/point-cloud-utils/branch/master)

**Point Cloud Utils (pcu)** is a utility library providing the following functionality:
 - A series of algorithms for generating point samples on meshes:
   - Poisson-Disk-Sampling of a mesh based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface](http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
   - Sampling a mesh with [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
   - Monte-Carlo sampling on a mesh
 - Normal estimation from point clouds
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

### With `pip` (not recommended)
```
pip install git+git://github.com/fwilliams/point-cloud-utils
```
The following dependencies are required to install with `pip`:
* A C++ compiler supporting C++14 or later
* CMake 3.2 or later.
* git

# Examples

### Poisson-Disk-Sampling
Generate 10000 samples on a mesh with poisson disk samples
```python
import point_cloud_utils as pcu
import numpy as np

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
# n is a nv by 3 NumPy array of vertex normals
v, f, n, _ = pcu.read_ply("my_model.ply")

# Generate 10000 samples on a mesh with poisson disk samples
# f_i are the face indices of each sample and bc are barycentric coordinates of the sample within a face
f_i, bc = pcu.sample_mesh_poisson_disk(v, f, n, 10000)

# Use the face indices and barycentric coordinate to compute sample positions and normals
v_poisson = (v[f[f_i]] * bc[:, np.newaxis]).sum(1)
n_poisson = (n[f[f_i]] * bc[:, np.newaxis]).sum(1)
```

Generate samples on a mesh with poisson disk samples seperated by approximately 0.01 times the boundinb box diagonal
```python
import point_cloud_utils as pcu
import numpy as np
# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
# n is a nv by 3 NumPy array of vertex normals
v, f, n, _ = pcu.read_ply("my_model.ply")


# Generate samples on a mesh with poisson disk samples seperated by approximately 0.01 times 
# the length of the bounding box diagonal
bbox = np.max(v, axis=0) - np.min(v, axis=0)
bbox_diag = np.linalg.norm(bbox)

# f_i are the face indices of each sample and bc are barycentric coordinates of the sample within a face
f_i, bc = pcu.sample_mesh_poisson_disk(v, f, n, 10000)

# Use the face indices and barycentric coordinate to compute sample positions and normals
v_poisson = (v[f[f_i]] * bc[:, np.newaxis]).sum(1)
n_poisson = (n[f[f_i]] * bc[:, np.newaxis]).sum(1)
    
```

### Monte-Carlo Sampling on a mesh
```python
import point_cloud_utils as pcu
import numpy as np

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
# n is a nv by 3 NumPy array of vertex normals
v, f, n, _ = pcu.read_ply("my_model.ply")

# Generate very dense random samples on the mesh (v, f, n)
# f_i are the face indices of each sample and bc are barycentric coordinates of the sample within a face
f_idx, bc = pcu.sample_mesh_random(v, f, num_samples=v.shape[0] * 40)

# Use the face indices and barycentric coordinate to compute sample positions and normals
v_dense = (v[f[f_idx]] * bc[:, np.newaxis]).sum(1)
n_dense = (n[f[f_idx]] * bc[:, np.newaxis]).sum(1)
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

# Generate 100 points on the unit cube with Lloyd's algorithm
samples_3d = pcu.lloyd_3d(100)
```

### Estimating Normals From a Point Cloud
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
v, _, _, _ = pcu.read_ply("my_model.ply")

# Estimate a normal at each point (row of v) using its 16 nearest neighbors
n = pcu.estimate_point_cloud_normals(n, k=16)
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
w_a = np.ones(a.shape[0])
w_b = np.ones(b.shape[0])

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

chamfer_dist = pcu.chamfer_distance(a, b)
```

### Hausdorff Distances Between Point-Clouds
```python
import point_cloud_utils as pcu
import numpy as np

# Generate two random point sets
a = np.random.rand(1000, 3)
b = np.random.rand(500, 3)

# Compute one-sided squared Hausdorff distances
hausdorff_a_to_b = pcu.one_sided_hausdorff_distance(a, b)
hausdorff_b_to_a = pcu.one_sided_hausdorff_distance(b, a)

# Take a max of the one sided squared  distances to get the two sided Hausdorff distance
hausdorff_dist = pcu.hausdorff_distance(a, b)

# Find the index pairs of the two points with maximum shortest distancce
hausdorff_b_to_a, idx_b, idx_a = pcu.one_sided_hausdorff_distance(b, a, return_index=True)
assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_b_to_a) < 1e-5, "These values should be almost equal"

# Find the index pairs of the two points with maximum shortest distancce
hausdorff_dist, idx_b, idx_a = pcu.hausdorff_distance(b, a, return_index=True)
assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_dist) < 1e-5, "These values should be almost equal"

```

### K-Nearest-Neighbors Between Point Clouds
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
dists_a_to_b, corrs_a_to_b = pcu.shortest_distance_pairs(a, b)
```
