# Voxelizing a Triangle Mesh
Point-Cloud-Utils supports rasterizing triangle meshes to voxels as well as generating cube meshes from those voxels.

```python
import numpy as np
import point_cloud_utils as pcu

vox_size = 1.0 / 256.0  # Size of each voxel
vox_origin = [0, 0, 0]  # location of the bottom-left-back corner of the [0, 0, 0] voxel
v, f = pcu.load_mesh_vf("truck.ply")

# Voxelize the input mesh, vox_ijk is an array of integer voxel coordinates
vox_ijk = pcu.voxelize_triangle_mesh(vox_ijk, vox_size, vox_origin)

# Generate a cube mesh of voxels with a spacing of 0.02 voxels between each cube
cube_v, cube_f = pcu.voxel_grid_geometry(vox_ijk, vox_size, vox_origin, gap_fraction=0.02)
```

<p align="center">
  <div class="row">
    <div class="column" style="float: left; width: 47%; padding: 5px;">
      <img src="../../imgs/voxelize_2.png" alt="Mesh to voxelize" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 49%; padding: 5px;">
      <img src="../../imgs/voxelize_1.png" alt="Vizualized voxels with gaps" style="width:100%">
    </div>
  </div>
  <figcaption style="text-align: center; font-style: italic;"><b>Left:</b> Input mesh to voxelize. <b>Right:</b> Visualization of voxelized mesh.</figcaption>
</p>