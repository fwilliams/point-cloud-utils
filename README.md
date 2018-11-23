# Point Cloud Utilities (pcu)
A Python Library for 3D Point Clouds

**pcu** is a utility library for common tasks on point clouds. It provides the following operations:
 - A series of algorithms for generating point samples on meshes:
   - Poisson-Disk-Sampling of a mesh based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface (http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
   - Sampling a mesh with [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
   - Sampling a mesh uniformly
 - Clustering point-cloud vertices into bins
 - Very fast pairwise nearest neighbor between point clouds (based on [nanoflann](https://github.com/jlblancoc/nanoflann))

![Example of Poisson Disk Sampling](/img/blue_noise.png?raw=true "Example of Poisson Disk Sampling")

# Installation Instructions
Simply run:
```
pip install git+git://github.com/fwilliams/point-cloud-utils
```
The only dependencies required are a valid Python installation with SciPy, a C++ compiler supporting C++14 or later, and CMake 3.2 or later.

# Example Poisson Disk Sampling
```python
import point_cloud_utils as pcu

# Assume v is a nv by 3 Numpy array of vertices
# Assume f is an nf by 3 Numpy array of face indexes into v 
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
