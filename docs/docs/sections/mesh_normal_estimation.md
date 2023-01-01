# Computing Mesh Normals
Point-Cloud-Utils lets you easily compute both vertex and face normals for a mesh

## Per-Vertex Normals
```python
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("bunny.ply")

# n is a NumPy array with the same shape as v containing vertex normals
n = pcu.estimate_mesh_normals(v, f)
```
<p align="center">
    <img src="../../imgs/mesh_vertex_normals_2.png" alt="Estimated vertex normals for a triangle mesh" style="width:70%">
    <figcaption style="text-align: center; font-style: italic;">Estimating vertex normals for a triangle mesh.</figcaption>
</p>


## Per-Face Normals for a Mesh
```python
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("bunny.ply")

# n is a NumPy array with the same shape as f containing face normals
n = pcu.estimate_mesh_face_normals(v, f)
```
<p align="center">
    <img src="../../imgs/mesh_face_normals_2.png" alt="Estimated face normals for a triangle mesh" style="width:70%">
    <figcaption style="text-align: center; font-style: italic;">Estimating face normals for a triangle mesh</figcaption>
</p>