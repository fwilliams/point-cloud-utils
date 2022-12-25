<p align="center">
  <img src="imgs/logo_solarized_tint.png">
</p>

**Point Cloud Utils** is an _easy-to-use_ Python library for processing
and manipulating 3D point clouds and meshes.

## Installation
``` shell
pip install point-cloud-utils
```

## A very simple example
Point Cloud Utils uses NumPy arrays as fundamental data structure, making it very easy to integrate with existing numerical code.
For example, here's how to remove all points in a point cloud which are greater than some distance from a mesh.

``` py
import point_cloud_utils as pcu
# Load a mesh stored in my_mesh.ply:
#   v is a NumPy array of coordinates with shape (V, 3)
#   f is a NumPy array of face indices with shape (F, 3)
v, f = pcu.load_mesh_vf("my_mesh.ply")

# Load a point cloud stored in my_point_cloud.ply:
#   p is a NumPy array of point coordinates with shape (P, 3)
p = pcu.load_mesh_v("my_point_cloud.ply")

# Compute the shortest distance between each point in p and the mesh:
#   dists is a NumPy array of shape (P,) where dists[i] is the
#   shortest distnace between the point p[i, :] and the mesh (v, f)
dists, _, _ = pcu.closest_points_on_mesh(p, v, f)

# Delete all points which are farther than some distance away from the mesh
dist_thresh = 0.1
keep_points = p[dists < dist_thresh]

# Save the filtered point cloud to my_point_cloud_trimmed.ply
pcu.save_mesh_v("my_point_cloud_trimmed.ply", keep_points)
```



## Core Features
Point Cloud Utils includes utilities to perform the following tasks:

* [Mesh and Point Cloud I/O](sections/mesh_io). Point Cloud Utils can handle any file that can be opened in MeshLab.
* [Resampling Point Clouds](sections/point_cloud_resampling) to have different distributions.
* [Generating Point Samples on a Mesh](sections/mesh_sampling).
* [Computing Metrics Between Point Clouds](sections/shape_metrics) (e.g. Chamfer Distance, Hausdorff Distance, etc...).
* [Making a Mesh a Watertight Manifold](sections/watertight_manifold).
* [Computing Signed Distances (SDFs) to Meshes](sections/mesh_sdf).
* [Estimating Normals for Point Clouds](sections/point_cloud_normal_estimation).
* [Computing Mesh Normals](sections/mesh_normal_estimation).
* [Ray/Mesh Intersection](sections/ray_mesh_intersection).
* [Generating Surfels from Point Clouds](sections/surfels).
* [Smooting a Mesh](sections/mesh_smoothing).
* [Decimating a Mesh](sections/mesh_decimation).
* [Computing Connected Components of a Mesh](sections/computing_connected_components_of_a_mesh.md).
* [Deduplicating Point Clouds and Meshes](sections/deduplicating_point_clouds_and_meshes).
* [Calculating Mesh Face Areas](sections/calculating_mesh_face_areas).
* [Consistently Orienting Mesh Faces](sections/consistently_orienting_mesh_faces.md)


## [API Reference](sections/api_reference.md)
All available functions and classes are listed in the [API reference](sections/api_reference.md)