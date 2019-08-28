# Point Cloud Utils (pcu) - A Python library for common tasks on 3D point clouds

[![Build Status](https://travis-ci.org/fwilliams/point-cloud-utils.svg?branch=master)](https://travis-ci.org/fwilliams/point-cloud-utils)
[![Build status](https://ci.appveyor.com/api/projects/status/ujv44lqbeosgl9ij/branch/master?svg=true)](https://ci.appveyor.com/project/fwilliams/point-cloud-utils/branch/master)

**Point Cloud Utils (pcu)** is a utility library providing the following functionality:
 - A series of algorithms for generating point samples on meshes:
   - Poisson-Disk-Sampling of a mesh based on "[Parallel Poisson Disk Sampling with Spectrum Analysis on Surface](http://graphics.cs.umass.edu/pubs/sa_2010.pdf)".
   - Sampling a mesh with [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
   - Sampling a mesh uniformly
 - Clustering point-cloud vertices into bins
 - Very fast pairwise nearest neighbor between point clouds (based on [nanoflann](https://github.com/jlblancoc/nanoflann))
 - Hausdorff distances between point-clouds.
 - Chamfer distnaces between point-clouds.
 - Approximate Wasserstein distances between point-clouds using the [Sinkhorn](https://arxiv.org/abs/1306.0895) method.
 - Pairwise distances between point clouds
 - Utility functions for reading and writing common mesh formats (OBJ, OFF, PLY)



<script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js" integrity="sha256-Ae2Vz/4ePdIu6ZyI/5ZGsYnb+m0JlOmKPjt6XZ9JJkA=" crossorigin="anonymous"></script>
<script src="https://unpkg.com/@jupyter-widgets/html-manager@^0.18.0/dist/embed-amd.js" crossorigin="anonymous"></script>
<script type="application/vnd.jupyter.widget-state+json">
{
  "version_major": 2,
  "version_minor": 0,
  "state": {
    "c4da9a73f0174e10a87b2669312571a0": {
      "model_name": "RendererModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "_antialias": true,
        "_dom_classes": [],
        "_height": 600,
        "_width": 600,
        "camera": "IPY_MODEL_70aed017101a48598d911cce099cbb53",
        "clippingPlanes": [],
        "controls": [
          "IPY_MODEL_d65670096fea4f21932f240499fb3b52"
        ],
        "layout": "IPY_MODEL_8665d11883f4418a9a403aa07bf93365",
        "scene": "IPY_MODEL_f7432fe85a1e4b709370d44575d962ac",
        "shadowMap": "IPY_MODEL_2494ff4a3ba64a28a78d919d7ab07b8b"
      }
    },
    "70aed017101a48598d911cce099cbb53": {
      "model_name": "PerspectiveCameraModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "children": [
          "IPY_MODEL_fcd523faa9b840d6815e2403924036d6"
        ],
        "fov": 30.0,
        "position": [
          0.0,
          0.0,
          0.3695761236380554
        ],
        "quaternion": [
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "scale": [
          1.0,
          1.0,
          1.0
        ],
        "up": [
          0.0,
          1.0,
          0.0
        ]
      }
    },
    "fcd523faa9b840d6815e2403924036d6": {
      "model_name": "DirectionalLightModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "children": [],
        "color": "white",
        "intensity": 0.6,
        "position": [
          0.0,
          0.0,
          0.3695761236380554
        ],
        "quaternion": [
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "scale": [
          1.0,
          1.0,
          1.0
        ],
        "shadow": "uninitialized",
        "target": "uninitialized",
        "up": [
          0.0,
          1.0,
          0.0
        ]
      }
    },
    "d65670096fea4f21932f240499fb3b52": {
      "model_name": "OrbitControlsModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "controlling": "IPY_MODEL_70aed017101a48598d911cce099cbb53",
        "maxAzimuthAngle": "inf",
        "maxDistance": "inf",
        "maxZoom": "inf",
        "minAzimuthAngle": "-inf"
      }
    },
    "8665d11883f4418a9a403aa07bf93365": {
      "model_name": "LayoutModel",
      "model_module": "@jupyter-widgets/base",
      "model_module_version": "1.2.0",
      "state": {}
    },
    "f7432fe85a1e4b709370d44575d962ac": {
      "model_name": "SceneModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "children": [
          "IPY_MODEL_70aed017101a48598d911cce099cbb53",
          "IPY_MODEL_a3ced971b3604c7d914d04edc12d95f2",
          "IPY_MODEL_1ab65b433eef49b1b99c841aae78c03f"
        ],
        "fog": null,
        "overrideMaterial": null,
        "quaternion": [
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "scale": [
          1.0,
          1.0,
          1.0
        ],
        "up": [
          0.0,
          1.0,
          0.0
        ]
      }
    },
    "a3ced971b3604c7d914d04edc12d95f2": {
      "model_name": "AmbientLightModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "children": [],
        "intensity": 0.5,
        "quaternion": [
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "scale": [
          1.0,
          1.0,
          1.0
        ],
        "up": [
          0.0,
          1.0,
          0.0
        ]
      }
    },
    "1ab65b433eef49b1b99c841aae78c03f": {
      "model_name": "PointsModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "children": [],
        "geometry": "IPY_MODEL_958f938a6ce54e2e95a352080b0d4ec5",
        "material": "IPY_MODEL_f9db48a1e8b84b7381778f032b5a0f82",
        "quaternion": [
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "scale": [
          1.0,
          1.0,
          1.0
        ],
        "up": [
          0.0,
          1.0,
          0.0
        ]
      }
    },
    "958f938a6ce54e2e95a352080b0d4ec5": {
      "model_name": "BufferGeometryModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "_ref_geometry": null,
        "attributes": {
          "position": "IPY_MODEL_b84ae18401b44da092b9f55514f057c5",
          "color": "IPY_MODEL_da4833cc15d6471fa4981adc88e096fd"
        },
        "index": null,
        "morphAttributes": {},
        "userData": {}
      }
    },
    "b84ae18401b44da092b9f55514f057c5": {
      "model_name": "BufferAttributeModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "array": {
          "shape": [
            1000,
            3
          ],
          "dtype": "float32"
        }
      },
      "buffers": [
        {
          "encoding": "base64",
          "path": [
            "array",
            "buffer"
          ],
          "data": "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"
        }
      ]
    },
    "da4833cc15d6471fa4981adc88e096fd": {
      "model_name": "BufferAttributeModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "array": {
          "shape": [
            1000,
            3
          ],
          "dtype": "float32"
        }
      },
      "buffers": [
        {
          "encoding": "base64",
          "path": [
            "array",
            "buffer"
          ],
          "data": "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"
        }
      ]
    },
    "f9db48a1e8b84b7381778f032b5a0f82": {
      "model_name": "PointsMaterialModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "alphaTest": 0.5,
        "clippingPlanes": [],
        "defines": null,
        "map": "IPY_MODEL_0b2d913b83c4419b92b58574ae4d0bb1",
        "size": 0.022,
        "vertexColors": "VertexColors"
      }
    },
    "0b2d913b83c4419b92b58574ae4d0bb1": {
      "model_name": "DataTextureModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {
        "data": {
          "shape": [
            16,
            16,
            4
          ],
          "dtype": "float32"
        },
        "repeat": [
          1.0,
          1.0
        ],
        "type": "FloatType"
      },
      "buffers": [
        {
          "encoding": "base64",
          "path": [
            "data",
            "buffer"
          ],
          "data": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=="
        }
      ]
    },
    "2494ff4a3ba64a28a78d919d7ab07b8b": {
      "model_name": "WebGLShadowMapModel",
      "model_module": "jupyter-threejs",
      "model_module_version": "^2.1.0",
      "state": {}
    }
  }
}
</script>
<script type="application/vnd.jupyter.widget-view+json">
{"version_major": 2, "version_minor": 0, "model_id": "c4da9a73f0174e10a87b2669312571a0"}
</script>



# Installation Instructions
### With `conda` (recommended)
Simply run:
```
conda install -c conda-forge point_cloud_utils
```

### With `pip` (not recommended)
```
pip install git+git://github.com/fwilliams/point-cloud-utils
```
The following dependencies are required to install with `pip`:
* A C++ compiler supporting C++14 or later
* CMake 3.2 or later.
* git

# Examples

### Poisson-Disk-Sampling
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
# n is a nv by 3 NumPy array of vertex normals
v, f, n, _ = pcu.read_ply("my_model.ply")
bbox = np.max(v, axis=0) - np.min(v, axis=0)
bbox_diag = np.linalg.norm(bbox)

# Generate very dense  random samples on the mesh (v, f, n)
# Note that this function works with no normals, just pass in an empty array np.array([], dtype=v.dtype)
# v_dense is an array with shape (100*v.shape[0], 3) where each row is a point on the mesh (v, f)
# n_dense is an array with shape (100*v.shape[0], 3) where each row is a the normal of a point in v_dense
v_dense, n_dense = pcu.sample_mesh_random(v, f, n, num_samples=v.shape[0]*100)

# Downsample v_dense to be from a blue noise distribution: 
#
# v_poisson is a downsampled version of v where points are separated by approximately 
# `radius` distance, use_geodesic_distance indicates that the distance should be measured on the mesh.
#
# n_poisson are the corresponding normals of v_poisson
v_poisson, n_poisson = pcu.sample_mesh_poisson_disk(
    v_dense, f, n_dense, radius=0.01*bbox_diag, use_geodesic_distance=True)
```

### Lloyd Relaxation
```python
import point_cloud_utils as pcu

# v is a nv by 3 NumPy array of vertices
# f is an nf by 3 NumPy array of face indexes into v 
v, f, _, _ = pcu.read_ply("my_model.ply")

# Generate 1000 points on the mesh with Lloyd's algorithm
samples = pcu.sample_mesh_lloyd(v, f, 1000)

# Generate 100 points on the unit square with Lloyd's algorithm
samples_2d = pcu.lloyd_2d(100)
```

### Approximate Wasserstein (Sinkhorn) Distance Between Point-Clouds

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are arrays where each row contains a point 
# Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
a = np.random.rand(100, 3)
b = np.random.rand(100, 3)

# M is a 100x100 array where each entry  (i, j) is the squared distance between point a[i, :] and b[j, :]
M = pcu.pairwise_distances(a, b)

# w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
w_a = np.ones(a.shape[0])
w_b = np.ones(b.shape[0])

# P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to 
# better approximation of the true Wasserstein distance at the expense of slower convergence
P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

# To get the distance as a number just compute the frobenius inner product <M, P>
sinkhorn_dist = (M*P).sum() 
```


##### Batched Version:

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are each contain 10 batches each of which contain 100 points  of dimension 3
# i.e. a[i, :, :] is the i^th point set which contains 100 points 
# Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
a = np.random.rand(10, 100, 3)
b = np.random.rand(10, 100, 3)

# M is a 10x100x100 array where each entry (k, i, j) is the squared distance between point a[k, i, :] and b[k, j, :]
M = pcu.pairwise_distances(a, b)

# w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
w_a = np.ones(a.shape[:2])
w_b = np.ones(b.shape[:2])

# P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to 
# better approximation of the true Wasserstein distance at the expense of slower convergence
P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

# To get the distance as a number just compute the frobenius inner product <M, P>
sinkhorn_dist = (M*P).sum() 
```


### Chamfer Distance Between Point-Clouds
```python
import point_cloud_utils as pcu
import numpy as np

# a and b are arrays where each row contains a point 
# Note that the point sets can have different sizes (e.g [100, 3], [111, 3])
a = np.random.rand(100, 3)
b = np.random.rand(100, 3)

chamfer_dist = pcu.chamfer(a, b)
```
##### Batched Version:

```python
import point_cloud_utils as pcu
import numpy as np

# a and b are each contain 10 batches each of which contain 100 points  of dimension 3
# i.e. a[i, :, :] is the i^th point set which contains 100 points 
# Note that the point sets can have different sizes (e.g [10, 100, 3], [10, 111, 3])
a = np.random.rand(10, 100, 3)
b = np.random.rand(10, 100, 3)

chamfer_dist = pcu.chamfer(a, b)
```


### Nearest-Neighbors and Hausdorff Distances Between Point-Clouds
```python
import point_cloud_utils as pcu
import numpy as np

# Generate two random point sets
a = np.random.rand(1000, 3)
b = np.random.rand(500, 3)

# dists_a_to_b is of shape (a.shape[0],) and contains the shortest squared distance 
# between each point in a and the points in b
# corrs_a_to_b is of shape (a.shape[0],) and contains the index into b of the 
# closest point for each point in a
dists_a_to_b, corrs_a_to_b = pcu.point_cloud_distance(a, b)

# Compute each one sided squared Hausdorff distances
hausdorff_a_to_b = pcu.hausdorff(a, b)
hausdorff_b_to_a = pcu.hausdorff(b, a)

# Take a max of the one sided squared  distances to get the two sided Hausdorff distance
hausdorff_dist = np.max(hausdorff_a_to_b, hausdorff_b_to_a)

# Find the index pairs of the two points with maximum shortest distancce
hausdorff_b_to_a, idx_b, idx_a = pcu.hausdorff(b, a, return_index=True)
assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_b_to_a) < 1e-5, "These values should be almost equal"
```

