# Computing Connected Components of a Mesh
Point Cloud Utils lets you easily find connected components of a mesh using the `connected_components` function:
```python
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("tree.ply")

# * cv a [num_vertices,]-shaped array where cv[i] is the integer 
#   id of the connected component for vertex v[i]
# * nv is a [num_connected_components]-shaped array where nv[j] is
#   the number of vertices in connected component j
# * cf a [num_faces,]-shaped array where cf[i] is the integer 
#   id of the connected component for face f[i]
# * nv is a [num_connected_components]-shaped array where nf[j] is
#   the number of faces in connected component j
cv, nv, cf, nf = pcu.connected_components(v, f)
```

<p align="center">
    <img src="../../imgs/connected_components_3.png" alt="Mesh colored by connected components" style="width:80%">
    <figcaption style="text-align: center; font-style: italic;">A mesh where each face is colored according to the connected component that face belongs to.</figcaption>
</p>