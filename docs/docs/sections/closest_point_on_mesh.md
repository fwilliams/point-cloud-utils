# Finding Closest Points Between Point Clouds and Meshes
You can easily find the closest point on a mesh to any 3D point in space with Point-Cloud-Utils.

```python
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
dists, fid, bc = pcu.closest_points_on_mesh(p, v, f)

# Interpolate the barycentric coordinates to get the coordinates of 
# the closest points on the mesh to each point in p
closest_pts = pcu.interpolate_barycentric_coords(f, fid, bc, f)
```

