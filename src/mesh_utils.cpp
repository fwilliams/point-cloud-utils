#include <igl/per_vertex_normals.h>

#include <npe.h>
#include <npe_typedefs.h>

const char* ds_per_vertex_normals = R"igl_Qu8mg5v7(
Compute vertex normals of a mesh from its vertices and faces

Parameters
----------
v : #v by 3 Matrix of mesh vertex 3D positions
f : #f by 3 Matrix of face (triangle) indices
weighting : Weighting type (0 = uniform, 1 = area, 2 = angle, 3 = default)

Returns
-------
n : list of vertex normals of shape #v by 3

See also
--------

Notes
-----

Examples
--------
)igl_Qu8mg5v7";

npe_function(per_vertex_normals)
npe_doc(ds_per_vertex_normals)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_default_arg(weighting, int, 3)
npe_begin_code()
{
  npe_Matrix_v n;
  igl::PerVertexNormalsWeightingType wtype = igl::PerVertexNormalsWeightingType(weighting);
  igl::per_vertex_normals(v, f, wtype, n);

  return npe::move(n);
}
npe_end_code()
