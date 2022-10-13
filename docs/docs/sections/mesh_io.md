# Mesh and Point Cloud I/O
**If your file can be imported into MeshLab, we can read it!**

Point Cloud Utils supports reading and writing many common 3D formats (PLY, STL, OFF, OBJ, 3DS, VRML 2.0,
X3D, COLLADA).
In all the functions in this section, the type of file is inferred from its file extension.

## 3D Data Representation in Point Cloud Utils
Point Cloud Utils uses NumPy arrays as a fundamental data structure for storing 3D data.

### Representing Point Clouds
A point cloud with `#p` points is represented as a simple Numpy array with shape `(#p, 3)` (one point per row). The example below illustrates a point cloud with six points.
<p align="center">
  <img src="/imgs/mesh_format.svg">
</p>


### Representing Triangle Meshes
 A triangle mesh with `#v` vertices and `#f` faces is stored using a pair of NumPy arrays `v, f` with
shape `(#v, 3)` and `(#f, 3)` respectively. Each row of `v` stores a vertex coordinate and each row of `f` stores
three integer indexes into `v` denoting the three vertices forming that face. Per-vertex attributes are stored as
separate NumPy arrays with `#v` rows (one row per vertex).
<p align="center">
  <img src="/imgs/mesh_format_transparent.png">
</p>

For example, consider the mesh with 6 vertices and 5 faces illustrated above.
We would store this as two numpy arrays  `v, f` with `v.shape = (6, 3)` and `f.shape = (5, 3)`. The i<sup>th</sup> row
of `v` is a 3D vector encoding the position of the vertex _v<sub>i</sub>_, and j<sup>th</sup> row of
`f` holds three indices to each of the vertices of the j<sup>th</sup> face. In the figure, the second face is formed by
the vertices _v<sub>1</sub>, v<sub>2</sub>, v<sub>5</sub>_, so the second row of `f` are the integers `(1, 2, 5)`.


### Representing Point and Face Attributes
In addition to vertices and face data, we often want to store attributes alongside a point cloud or mesh. Attributes can be stored per-point as well as per-face. Attributes are stored as NumPy arrays with one attribute per row. Examples of attributes include point cloud normals, point colors, per-face colors, or per-vertex texture coordinates.

**Example - Storing face colors and vertex normals in a mesh:**
For example, consider the mesh described above with 6 vertices and 5 faces, stored as two arrays `v` and `f` with shapes `(6, 3)` and `(5, 3)` respectively. We could store per-face colors as
a NumPy array with shape `(5, 3)` each row is the color of the correponding face:

<p align="center">
  <img src="/imgs/mesh_format_attrib.svg" alt="Per face colors are stored in an array with the same number of rows as the face array">
  <figcaption style="text-align: center; font-style: italic;">Per-face attributes are stored in an array with the same number of rows as the face array</figcaption>
</p>


We could store per-vertex normals as a NumPy array of shape `(6, 3)` where each row contains the normal at the corresponding vertex:
<p align="center">
  <img src="/imgs/mesh_format_attrib_vertex.svg" alt="Per-vertex normals are stored in an array with the same number of rows as the vertex array">
  <figcaption style="text-align: center; font-style: italic">Per-vertex attributes (e.g. normals) are stored in an array with the same number of rows as the vertex array</figcaption>
</p>



## Shorthand functions for loading and saving 3D data
Many times, we only need to load a few attributes from a mesh or point cloud. For example, we may only want the vertices, faces, and vertex colors for mesh. Point Cloud Utils provides a number of
shorthand helper functions which can load these attributes directly into NumPy arrays. These functions have the form `load_mesh_*` and `save_mesh_*` where `*` indicates which data to load. If some attribute is not present in the file being loaded, these functions will simply return an empty array.


#### Load/Save only points
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
v = pcu.load_mesh_v("path/to/mesh")
pcu.save_mesh_v("path/to/mash", v)
```

#### Load/Save mesh vertices and faces
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
# f is a NumPy array of shape (#f, 3) of face indices
v, f = pcu.load_mesh_vf("path/to/mesh")
pcu.save_mesh_vf("path/to/mash", v, f)
```

#### Load/Save points and per-point normals
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
# n is a NumPy array of shape (#v, 3) of per-vertex normals
v, n = pcu.load_mesh_vn("path/to/mesh")
pcu.save_mesh_vn("path/to/mash", v, n)
```

#### Load/Save mesh vertices, faces, and vertex colors
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
# f is a NumPy array of shape (#f, 3) of face indices
# c is a NumPy array of shape (#v, 4) of RBGA per-vertex colors
v, f, c = pcu.load_mesh_vfc("path/to/mesh")
pcu.save_mesh_vfc("path/to/mash", v, f, c)
```

#### Load/Save mesh vertices, faces, and vertex normals
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
# f is a NumPy array of shape (#f, 3) of face indices
# n is a NumPy array of shape (#v, 3) of per-vertex normals
v, f, n = pcu.load_mesh_vfn("path/to/mesh")
pcu.save_mesh_vfn("path/to/mash", v, f, n)
```

#### Load/Save mesh vertices, faces, vertex normals, and vertex colors
```python
import point_cloud_utils as pcu

# v is a NumPy array of shape (#v, 3) of points
# f is a NumPy array of shape (#f, 3) of face indices
# n is a NumPy array of shape (#v, 3) of per-vertex normals
# c is a NumPy array of shape (#v, 4) of RBGA per-vertex colors
v, f, n, c = pcu.load_mesh_vfn("path/to/mesh")
pcu.save_mesh_vfnc("path/to/mash", v, f, n, c)
```

## Loading 3D data with all attributes
Some meshes or point clouds may have more complex attribute structures, or we may not know which attributes are stored
in the mesh file before loading it. In this case, Point Cloud Utils provides utilities to load a mesh into
a `TriangleMesh` class. `TriangleMesh` is a lightweight container containing numpy arrays for
vertices, faces, textures, and other attributes.

To load a mesh into a `TriangleMesh` class simply construct it with the path to the mesh:

```python
import point_cloud_utils as pcu

mesh = pcu.TriangleMesh("path/to/mesh")
```
To save a `TriangleMesh` object, simply call the `save` method with the desired path:
```python
mesh.save("path/to/mesh")
```

The `TriangleMesh` class contains attributes encoding vertex and face data as well as texture information. It is structured as follows:

* **`TriangleMesh`**:
    * **`vertex_data`**:
        * **`positions`**: (#v, 3)-shaped NumPy array of per-vertex positions
        * **`normals`**: (#v, 3)-shaped NumPy array of per-vertex normals (or `None`)
        * **`texcoords`**: (#v, 2)-shaped NumPy array of per-vertex uv coordinates (or `None`)
        * **`tex_ids`**: (#v,)-shaped NumPy array of integer indices into TriangleMesh.textures indicating which texture to use at this vertex (or `None`)
        * **`colors`**: (#v, 4)-shaped NumPy array of per-vertex RBGA colors in `[0.0, 1.0]` (or `None`)
        * **`radius`**: (#v,)-shaped NumPy array of per-vertex curvature radii (or `None`)
        * **`quality`**: (#v,)-shaped NumPy array of per-vertex quality measures (or `None`)
        * **`flags`**: (#v,)-shaped NumPy array of 32-bit integer flags per vertex (or `None`)
    * **`face_data`**:
        * **`vertex_ids`**: (#f, 3)-shaped NumPy array of integer face indices into `TrianglMesh.vertex_data.positions`
        * **`normals`**: (#f, 3)-shaped NumPy array of per-face normals (or `None`)
        * **`colors`**: (#f, 4)-shaped NumPy array of per-face RBGA colors in `[0.0, 1.0]` (or `None`)
        * **`quality`**: (#f,)-shaped NumPy array of per-face quality measures (or `None`)
        * **`flags`**: (#f,(-shaped NumPy array of 32-bit integer flags per face (or `None`)
        * **`wedge_colors`**: (#f, 3, 4)-shaped NumPy array of per-wedge RBGA colors in `[0.0, 1.0]` (or `None`)
        * **`wedge_normals`**: (#f, 3, 3)-shaped NumPy array of per-wedge normals (or `None`)
        * **`wedge_texcoords`**: (#f, 3, 2)-shaped NumPy array of per-wedge] uv coordinates (or `None`)
        * **`wedge_tex_ids`**: (#f, 3)-shaped NumPy array of integer indices into `TriangleMesh.textures` indicating which texture to use at this wedge (or `None`)
    * **`textures`**: A list of paths to texture image files for this mesh
    * **`normal_maps`**: A list of paths to texture image files for this mesh

The hierarchy of the list above denotes composition. For example, to access the vertex colors of a `TriangleMesh`, you would
read the `TriangleMesh.vertex_data.colors` property.

!!! note "Remark on wedge face attributes"
    The `face_data` member of the `TriangleMesh` class contains a number of `wedge` attributes. These all start with
    `wedge_`. In this context a _wedge_ is a corner of a triangle face. Each face contains three wedges corresponding
     to each corner. Wedge attributes have shape `(#f, 3, d)` where `d` is the dimension of the attribute.
     Wedges attributes are useful for storing which are discontinous across face boundaries.


<!--
```text
TriangleMesh:
  vertex_data:
      positions: (#v, 3)-shaped numpy array of per-vertex positions
      normals: (#v, 3)-shaped numpy array of per-vertex normals (or None)
      texcoords: (#v, 2)-shaped numpy array of per-vertex uv coordinates (or None)
      tex_ids: (#v,)-shaped numpy array of integer indices into TriangleMesh.textures indicating which texture to
               use at this vertex (or None)
      colors: (#v, 4)-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] (or None)
      radius: (#v,)-shaped numpy array of per-vertex curvature radii (or None)
      quality: (#v,)-shaped numpy array of per-vertex quality measures (or None)
      flags: (#v,)-shaped numpy array of 32-bit integer flags per vertex (or None)
  face_data:
      vertex_ids: (#f, 3)-shaped numpy array of integer face indices into TrianglMesh.vertex_data.positions
      normals: (#f, 3)-shaped numpy array of per-face normals (or None)
      colors: (#f, 4)-shaped numpy array of per-face RBGA colors in [0.0, 1.0] (or None)
      quality: (#f,)-shaped numpy array of per-face quality measures (or None)
      flags: (#f,(-shaped numpy array of 32-bit integer flags per face (or None)

      wedge_colors: (#f, 3, 4)-shaped numpy array of per-wedge RBGA colors in [0.0, 1.0] (or None)
      wedge_normals: (#f, 3, 3)-shaped numpy array of per-wedge normals (or None)
      wedge_texcoords: (#f, 3, 2)-shaped numpy array of per-wedge] uv coordinates (or None)
      wedge_tex_ids: (#f, 3)-shaped numpy array of integer indices into TriangleMesh.textures indicating which
                     texture to use at this wedge (or None)
  textures: A list of paths to texture image files for this mesh
  normal_maps: A list of paths to texture image files for this mesh
``` --->