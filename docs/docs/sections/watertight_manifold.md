# Making a Mesh a Watertight Manifold
Many algorithms in geometry processing require a mesh to be *manifold*, and *watertight*. For example, [computing signed distances from meshes](sections/mesh_sdf.md) requires watertight shapes.

Point-Cloud-Utils implements the [Robust Watertight Manifold Surface Generation Method for ShapeNet Models](https://arxiv.org/abs/1802.01698) algorithm by Huang et.al. for converting meshes to watertight manifolds. 

```python
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("chair.ply")


# The resolution parameter controls the density of the output mesh
# It is linearly proportional to the number of faces in the output
# mesh. A higher value corresponds to a denser mesh.
resolution = 50_000
vw, fw = pcu.make_mesh_watertight(v, f, resolution)
```
<p align="center">
  <div class="row">
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/non_manifold.png" alt="Non-manifold, non-watertight mesh" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/manifold.png" alt="Watertight, Manifold version of the input mesh" style="width:100%">
    </div>
    <figcaption style="text-align: center; font-style: italic;">The left mesh is non-manifold and non-watertight. The Manifold algorithm converts it to a watertight manifold on the left.</figcaption>
  </div>
</p>