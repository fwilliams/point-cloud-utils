#include <npe.h>

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <string>

#include "common/common.h"



const char* estimate_mesh_vertex_normals_doc = R"igl_Qu8mg5v7(
Compute vertex normals of a mesh from its vertices and faces using face area weighting

Args:
    v : (#v, 3)-shaped NumPy array of mesh vertex 3D positions
    f : (#f, 3)-shaped NumPy array of face (triangle) indices
    weighting_type : Weighting type must be one of 'uniform', 'angle', or 'area' (default is 'uniform')

Returns:
    n : (#v, 3)-shaped NumPy array of vertex normals (i.e. n[i] is the normal at vertex v[i]) estimate_mesh_face_normals

)igl_Qu8mg5v7";
npe_function(estimate_mesh_vertex_normals)
npe_doc(estimate_mesh_vertex_normals_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_default_arg(weighting_type, std::string, std::string("uniform"))
npe_begin_code()
{
    validate_mesh(v, f);
    int weighting = -1;
    if (weighting_type == "uniform") {
        weighting  = 0;
    } else if (weighting_type == "area") {
        weighting = 1;
    } else if (weighting_type == "angle") {
        weighting = 2;
    } else {
        throw pybind11::value_error("Invalid weighting_type must be one of 'uniform', 'angle', or 'area', but got '" +
                                    weighting_type + "'");
    }

    npe_Matrix_v n;
    igl::PerVertexNormalsWeightingType wtype = igl::PerVertexNormalsWeightingType(weighting);
    igl::per_vertex_normals(v, f, wtype, n);

    return npe::move(n);
}
npe_end_code()



const char* estimate_mesh_face_normals_doc = R"igl_Qu8mg5v7(
Compute vertex normals of a mesh from its vertices and faces using face area weighting

Args:
    v : (#v, 3)-shaped NumPy array of mesh vertex 3D positions
    f : (#f, 3)-shaped NumPy array of face (triangle) indices

Returns:
    n : (#f, 3)-shaped NumPy array of face normals (i.e. n[i] is the normal at face f[i]). Note that any degenerate faces will have a zero normal.

)igl_Qu8mg5v7";
npe_function(estimate_mesh_face_normals)
npe_doc(estimate_mesh_face_normals_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_begin_code()
{
    validate_mesh(v, f);

    npe_Matrix_v n;
    EigenDenseLike<npe_Matrix_v> bg(3, 1);
    bg.setZero();
    igl::per_face_normals(v, f, bg, n);

    return npe::move(n);
}
npe_end_code()
