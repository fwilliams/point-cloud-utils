# py-sample-mesh
A Python Library for Generating Point Samples on Meshes

py-sample mesh is a very simple library to add random samples to a triangle mesh. It provides the following operations:
 - Poisson-Disk-Sampling based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface](http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
 - Clustering point-cloud vertices into bins
 - Uniform random sampling on meshes

![Example of Poisson Disk Sampling](https://github.com/fwilliams/py-sample-mesh/edit/master/img/blue_noise.png)

# Installation Instructions
Simply run:
```
pip install git+git://github.com/fwilliams/py-sample-mesh
```
The only dependency required is a valid Python installation with SciPy and a C++ compiler supporting C++14 or later.

# Example Usage
```python
import sample_mesh as sm

# Assume v is a nv by 3 Numpy array of vertices
# Assum f is an nf by 3 Numpy array of face indexes into v 
bbox = np.max(v, axis=0) - np.min(v, axis=0)
bbox_diag = np.linalg.norm(bbox)

# v_poisson is a downsampled version of v where points are separated
# by approximately `radius` distance, use_geodesic_distance indicates
# that the distance should be measured on the mesh.
v_poisson = sm.poisson_disk_sample(v, f, radius=0.01*bbox_diag, use_geodesic_distance=True)

# v_poisson is a downsampled version of v where points are separated
# by approximately `radius` distance
v_cluster = sm.cluster_vertices(v, cell_size=0.01*bbox_diag)

# v_random_samples contains `num_samples` samples placed uniformly at random on
# the mesh
ns = 1024
v_random_samples = sm.random_sample(v, f, num_samples=ns)
```
