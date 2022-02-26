# Point Cloud Re-Sampling
Point Cloud Utils provides utilities for re-sampling a point cloud in several ways.

## Binning points in a voxel grid
In many 3D tasks, it is often advantageous to quantize a point cloud to a regular voxel grid in space.
Such quantization can reduce the number of input points to an algorithm to match the resolution 
at which that algorithm operates.

Point Cloud Utils provides a simple function for quantizing a point cloud with attributes into a voxel grid.

```python
import point_cloud_utils as pcu

pcu.downsample_point_cloud_voxel_grid
```

## Downsampling a point cloud to have a blue noise distribution
Often times, we want to resample a dense point cloud to have approximately even spacing between points. 
This can be achieved with Poisson Disk Sampling. Point Cloud utils provides a r
## Deduplicating a point cloud

```python
import point_cloud_utils as pcu
pcu.deduplicate_point_cloud
pcu.downsample_point_cloud_voxel_grid
pcu.downsample_point_cloud_poisson_disk
```