# Smoothing a Mesh
Point Cloud Utils supports simple [Laplacian Smoothing](https://en.wikipedia.org/wiki/Laplacian_smoothing) of a triangle mesh.

```python
import numpy as np
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("bunny.ply")
n = pcu.estimate_mesh_vertex_normals(v, f)

# Add Gaussian noise to each vertex along its normal
v += np.random.randn(v.shape[0], 1) * n * 0.02

# Run 4 iterations of Laplacian smoothing. use_cotan_weights is an 
# optional parameter specifying whether to use uniform or cotangent
# weights.
# v_smooth has the same shape as v and contains the smooth vertices
v_smooth = pcu.laplacian_smooth_mesh(v, f, num_iters=4, use_cotan_weights=True)
```
<p align="center">
  <div class="row">
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/mesh_noisy.png" alt="A noisy triangle mesh" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 50%; padding: 5px;">
      <img src="../../imgs/mesh_smooth_4_iters.png" alt="A smoothed version of the noisy mesh" style="width:100%">
    </div>
    <figcaption style="text-align: center; font-style: italic;"><b>Left:</b> A noisy triangle mesh. <b>Right:</b> The noisy mesh after 4 iterations of Laplacian smoothing.</figcaption>
  </div>
</p>