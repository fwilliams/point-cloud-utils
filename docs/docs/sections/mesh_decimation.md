# Decimating a Mesh
Point Cloud Utils supports mesh decimation using the [Quadric Edge Collapse](https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf) algorithm.
```python
import numpy as np
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("bunny.ply")

# v_quarter, f_quarter are the mesh vertices and faces of the decimated mesh
# corr_qv and corr_qf are the correspondences between vertices and faces in the 
# original mesh and the decimated mesh. 
# i.e. corr_qv[i] is the index of the vertex in v which generated v_quarter[i] 
#      corr_qf[i] is the index of the face in f which generated f_quarter[i]
v_quarter, f_quarter, corr_qv, corr_qf = pcu.decimate_triangle_mesh(v, f, max_faces=f.shape[0]//4)

v_eighth, f_eighth, _, _ = pcu.decimate_triangle_mesh(v, f, max_faces=f.shape[0]//8)
```

<p align="center">
    <img src="../../imgs/decimation.png" alt="Decimation hierarchy of meshes" style="width:100%">
    <figcaption style="text-align: center; font-style: italic;"><b>From Left to Right:</b> full mesh, a quarter of the number of faces, an eighth the number of faces.</figcaption>
</p>
