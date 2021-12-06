# Point Cloud Utils (pcu) - A Python library for common tasks on 3D point clouds (MPLv2)

**Modified license version (MPLv2): Utils and libraries from the original point-cloud-utils that were licensed under GNU/GPL have been removed from this fork.**

**Point Cloud Utils (pcu)** is a utility library providing the following functionality. See the [Examples section](#Examples) for documentation on how to use these:
 - With this fork, data needs to be imported via numpy. You can use several external libraries (e.g. trimesh or open3d) to read data from
   file and use arrays from those directly with functionality in this library.
 - Closest points between a point cloud and a mesh
 - Fast k-nearest-neighbor search between point clouds (based on [nanoflann](https://github.com/jlblancoc/nanoflann)).
 - Hausdorff distances between point-clouds.
 - Chamfer distances between point-clouds.
 - Approximate Wasserstein distances between point-clouds using the [Sinkhorn](https://arxiv.org/abs/1306.0895) method.
 - Compute signed distances between a point cloud and a mesh using [Fast Winding Numbers](https://www.dgp.toronto.edu/projects/fast-winding-numbers/)
 - Compute closest points on a mesh to a point cloud
 - Deduplicating point clouds and mesh vertices
 - Making a mesh watertight (based on the [Watertight Manifold](https://github.com/hjwdzh/Manifold) algorithm)
![Example of Poisson Disk Sampling](/img/blue_noise.png?raw=true "Example of Poisson Disk Sampling")

# Installation Instructions
### With `pip`
```
pip install git+git://github.com/{path to this repository}
```
The following dependencies are required to install with `pip`:
* A C++ compiler supporting C++14 or later
* git

### Install from wheel
A pre-built wheel for Python 3.7 has been included in ./wheel/point_cloud_utils-0.21.0-cp37-cp37m-win_amd64.whl
```
pip install path/to/point_cloud_utils-0.21.0-cp37-cp37m-win_amd64.whl
```

# Build Instructions
### Clone the repository and run (outputs into ./build):
```
python setup.py build
```
* Requires a C++ compiler supporting C++14 or later

### Create a wheel from the build (outputs into ./dist):

```
python setup.py bdist_wheel
```

# Examples

### List of examples
- [Point Cloud Utils (pcu) - A Python library for common tasks on 3D point clouds (MPLv2)](#point-cloud-utils-pcu---a-python-library-for-common-tasks-on-3d-point-clouds-mplv2)
- [Installation Instructions](#installation-instructions)
    - [With `pip`](#with-pip)
    - [Install from wheel](#install-from-wheel)
- [Build Instructions](#build-instructions)
    - [Clone the repository and run (outputs into ./build):](#clone-the-repository-and-run-outputs-into-build)
    - [Create a wheel from the build (outputs into ./dist):](#create-a-wheel-from-the-build-outputs-into-dist)
- [Examples](#examples)
    - [List of examples](#list-of-examples)
    - [Loading meshes and point clouds](#loading-meshes-and-point-clouds)
    - [Compute closest points on a mesh](#compute-closest-points-on-a-mesh)
    - [Approximate Wasserstein (Sinkhorn) distance between two point clouds](#approximate-wasserstein-sinkhorn-distance-between-two-point-clouds)
    - [Chamfer distance between two point clouds](#chamfer-distance-between-two-point-clouds)
    - [Hausdorff distance between two point clouds](#hausdorff-distance-between-two-point-clouds)
    - [K-nearest-neighbors between two point clouds](#k-nearest-neighbors-between-two-point-clouds)
    - [Compute shortest signed distances to a triangle mesh with fast winding numbers](#compute-shortest-signed-distances-to-a-triangle-mesh-with-fast-winding-numbers)
    - [Deduplicating Point Clouds and Meshes](#deduplicating-point-clouds-and-meshes)
      - [Point Clouds:](#point-clouds)
      - [Meshes:](#meshes)
    - [Making a Mesh Watertight](#making-a-mesh-watertight)


### Loading meshes and point clouds
With this fork, data needs to be imported via numpy. You can use several external libraries (e.g. trimesh or open3d) to read data from file and use arrays from those directly with functionality in this library.
```python
# Import from raw numpy array files
import numpy as np
# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v
# n is a nv by 3 NumPy array of vertex normals if they are specified, otherwise an empty array
v,f,n = np.load("../data/1_vertices.npy"), np.load("../data/1_faces.npy"), np.load("../data/1_vertex_normals.npy")

# Import from trimesh
import trimesh
mesh = trimesh.load_mesh("path/to/mesh")
v, f, n = mesh.vertices, mesh.faces, mesh.vertex_normals

```

For meshes and point clouds with more complex attributes, use a `TriangleMesh` object.

```python
import point_cloud_utils as pcu

# mesh is a lightweight TriangleMesh container object holding mesh vertices, faces, and their attributes.
# Any attributes which aren't loaded (because they aren't present in the file) are set to None.
# The data in TriangleMesh is layed out as follows (run help(pcu.TriangleMesh) for more details):
# TriangleMesh:
#   vertex_data:
#       positions: [V, 3]-shaped numpy array of per-vertex positions
#       normals: [V, 3]-shaped numpy array of per-vertex normals (or None)
#       texcoords: [V, 2]-shaped numpy array of per-vertex uv coordinates (or None)
#       tex_ids: [V,]-shaped numpy array of integer indices into TriangleMesh.textures indicating which texture to
#                use at this vertex (or None)
#       colors: [V, 4]-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] (or None)
#       radius: [V,]-shaped numpy array of per-vertex curvature radii (or None)
#       quality: [V,]-shaped numpy array of per-vertex quality measures (or None)
#       flags: [V,]-shaped numpy array of 32-bit integer flags per vertex (or None)
#   face_data:
#       vertex_ids: [F, 3]-shaped numpy array of integer face indices into TrianglMesh.vertex_data.positions
#       normals: [F, 3]-shaped numpy array of per-face normals (or None)
#       colors: [F, 4]-shaped numpy array of per-face RBGA colors in [0.0, 1.0] (or None)
#       quality: [F,]-shaped numpy array of per-face quality measures (or None)
#       flags: [F,]-shaped numpy array of 32-bit integer flags per face (or None)
#
#       wedge_colors: [F, 3, 4]-shaped numpy array of per-wedge RBGA colors in [0.0, 1.0] (or None)
#       wedge_normals: [F, 3, 3]-shaped numpy array of per-wedge normals (or None)
#       wedge_texcoords: [F, 3, 2]-shaped numpy array of per-wedge] uv coordinates (or None)
#       wedge_tex_ids: [F, 3]-shaped numpy array of integer indices into TriangleMesh.textures indicating which
#                      texture to use at this wedge (or None)
#   textures: A list of paths to texture image files for this mesh
#   normal_maps: A list of paths to texture image files for this mesh

#   Initialise a TriangleMesh from arrays

v,f,n,fn= np.load("../data/1_vertices.npy"), np.load("../data/1_faces.npy"), np.load("../data/1_vertex_normals.npy"), np.load("../data/_face_normals.npy")
face_data = pcu.TriangleMesh.FaceData()
face_data.normals = fn
face_data.vertex_ids = f
vertex_data = pcu.TriangleMesh.VertexData()
vertex_data.positions = v
vertex_data.normals = n
mesh = pcu.TriangleMesh()
mesh.face_data = face_data
mesh.vertex_data = vertex_data

```

### Compute closest points on a mesh
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
v,f= np.load("../data/1_vertices.npy"), np.load("../data/1_faces.npy")

# Generate 1000 random query points. We will find the closest point on the mesh for each of these
p = np.random.rand(1000, 3)

# For each query point, find the closest point on the mesh.
# Here:
#  - d is an array of closest distances for each query point with shape (1000,)
#  - fi is an array of closest face indices for each point with shape (1000,)
#  - bc is an array of barycentric coordinates within each face (shape (1000, 3)
#    of the closest point for each query point
d, fi, bc = pcu.closest_points_on_mesh(p, v, f)

# Convert barycentric coordinates to 3D positions
closest_points = pcu.interpolate_barycentric_coords(f, fi, bc, v)

# Alternatively broadcast the distances
filter_dist = 10.
ids_within_dist = np.where((d < filter_dist))
closest_points = p[ids_within_dist]
```

### Approximate Wasserstein (Sinkhorn) distance between two point clouds
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

### Chamfer distance between two point clouds
```python
import point_cloud_utils as pcu
import numpy as np

# a and b are arrays where each row contains a point
# Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
a = np.random.rand(100, 3)
b = np.random.rand(100, 3)

chamfer_dist = pcu.chamfer_distance(a, b)
```

### Hausdorff distance between two point clouds
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

### K-nearest-neighbors between two point clouds
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
# k is number of neighbors to query per point.
dists_a_to_b, corrs_a_to_b = pcu.k_nearest_neighbors(a, b, k=1)
```

### Compute shortest signed distances to a triangle mesh with [fast winding numbers](https://www.dgp.toronto.edu/projects/fast-winding-numbers/)
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v
v,f= np.load("../data/1_vertices.npy"), np.load("../data/1_faces.npy")

# Generate 1000 points in the volume around the mesh. We'll compute the signed distance to the
# mesh at each of these points
pts = np.random.rand(1000, 3) * (v.max(0) - v.min(0)) + v.min(0)

# Compute the sdf, the index of the closest face in the mesh, and the barycentric coordinates of
# closest point on the mesh, for each point in pts
sdfs, face_ids, barycentric_coords = pcu.signed_distance_to_mesh(pts, v, f)
```


### Deduplicating Point Clouds and Meshes
#### Point Clouds:
```python
import point_cloud_utils as pcu

# p is a (n, 3)-shaped array of points (one per row)
# p is a (n, 3)-shaped array of normals at each point
v,n= np.load("../data/1_vertices.npy"), np.load("../data/1_vertex_normals.npy"))

# Treat any points closer than 1e-7 apart as the same point
# idx_i is an array of indices such that p_dedup = p[idx_i]
# idx_j is an array of indices such that p = p_dedup[idx_j]
p_dedup, idx_i, idx_j  = deduplicate_point_cloud(p, 1e-7)

# Use idx_i to deduplicate the normals
n_dedup = n[idx_i]
```

#### Meshes:
```python
# v is a (nv, 3)-shaped NumPy array of vertices
# f is an (nf, 3)-shaped NumPy array of face indexes into v
# c is a (nv, 4)-shaped numpy array of per-vertex colors
v,f,c= np.load("../data/1_vertices.npy"), np.load("../data/1_vertex_normals.npy"), np.load("../data/1_vertex_colors.npy"))

# Treat any points closer than 1e-7 apart as the same point
# idx_i is an array of indices such that v_dedup = v[idx_i]
# idx_j is an array of indices such that v = v_dedup[idx_j]
v_dedup, f_dedup, idx_i, idx_j = pcu.deduplicate_mesh_vertices(v, f, 1e-7)

# Use idx_i to deduplicate the colors
c_dedup = c[idx_i]
```


### Making a Mesh Watertight
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v
v,f= np.load("../data/1_vertices.npy"), np.load("../data/1_vertex_normals.npy"))

# Optional resolution parameter (default is 20_000).
# See https://github.com/hjwdzh/Manifold for details
resolution = 20_000  
v_watertight, f_watertight = pcu.make_mesh_watertight(v, f, resolution=resolution)
```
