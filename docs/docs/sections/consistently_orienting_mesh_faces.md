# Consistently Orienting Mesh Faces
Often, meshes in the wild have inconsistently oriented faces, which can cause issues (e.g. due to flipped normals). Point Cloud Utils lets you consistently orient normals within each component
```python
import numpy as np
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("truck.ply")

# * f_oriented is a new face array where faces within each connected
#   component are consistently oriented
# * f_comp is a [num_faces,]-shaped array where f_comp[i] is the 
#   connected component of the i^th face
f_oriented, f_comp = pcu.orient_mesh_faces(f)
```

<p align="center">
  <div class="row">
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/oriented_bad.png" alt="Input ShapeNet Mesh" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/oriented_good.png" alt="Mesh with consistently oriented faces" style="width:100%">
    </div>
  </div>
  <figcaption style="text-align: center; font-style: italic;"><b>Left:</b> Mesh with inconsistently oriented faces. <b>Right:</b> Cleaned mesh with consistently oriented faces.</figcaption>
</p>