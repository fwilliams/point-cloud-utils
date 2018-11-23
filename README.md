# Point Cloud Utilities (pcu) - A Python library for common task on 3D point clouds

**pcu** is a utility library for common tasks on point clouds. It provides the following operations:
 - A series of algorithms for generating point samples on meshes:
   - Poisson-Disk-Sampling of a mesh based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface](http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
   - Sampling a mesh with [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
   - Sampling a mesh uniformly
 - Clustering point-cloud vertices into bins
 - Very fast pairwise nearest neighbor between point clouds (based on [nanoflann](https://github.com/jlblancoc/nanoflann))
 - Hausdorff distances between point-clouds.
 
![Example of Poisson Disk Sampling](/img/blue_noise.png?raw=true "Example of Poisson Disk Sampling")

# Installation Instructions
Simply run:
```
pip install git+git://github.com/fwilliams/point-cloud-utils
```
The only dependencies required are a valid Python installation with SciPy, a C++ compiler supporting C++14 or later, and CMake 3.2 or later.

# Examples

### Poisson-Disk-Sampling
```python
import point_cloud_utils as pcu

# Assume v is a nv by 3 NumPy array of vertices
# Assume f is an nf by 3 NumPy array of face indexes into v 
bbox = np.max(v, axis=0) - np.min(v, axis=0)
bbox_diag = np.linalg.norm(bbox)

# Generate very dense  random samples on the mesh (v, f)
# v_dense is an array with shape (100*v.shape[0], 3) where each row is a point on the mesh (v, f)
v_dense = pcu.sample_mesh_uniform(v, f, num_samples=v.shape[0]*100)

# Downsample v_dense to be from a blue noise distribution 
# v_poisson is a downsampled version of v where points are separated by approximately 
# `radius` distance, use_geodesic_distance indicates that the distance should be measured on the mesh.
v_poisson = sm.sample_mesh_poisson_disk(v, f, radius=0.01*bbox_diag, use_geodesic_distance=True)
```

### Lloyd Relaxation
```python
import point_cloud_utils as pcu

# Assume v is a nv by 3 NumPy array of vertices
# Assume f is an nf by 3 NumPy array of face indexes into v 

# Generate 1000 points on the mesh with Lloyd's algorithm
samples = pcu.sample_mesh_lloyd(v, f, 1000)

# Generate 100 points on the unit square with Lloyd's algorithm
samples_2d = pcu.lloyd_2d(100)
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
hausdorff_dist = max(hausdorff_a_to_b, hausdorff_b_to_a)

# Find the index pairs of the two points with maximum shortest distancce
hausdorff_b_to_a, idx_a, idx_b = pcu.hausdorff(b, a, return_index=True)
assert np.abs(np.linalg.norm(a[idx_a] - b[idx_b])**2 - hausdorff_b_to_a) < 1e-5, \
       "These values should be close"
```

