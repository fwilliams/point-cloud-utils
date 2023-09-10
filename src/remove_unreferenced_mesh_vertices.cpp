#include <npe.h>

#include <igl/remove_unreferenced.h>

#include "common/common.h"


const char* remove_unreferenced_mesh_vertices_doc = R"igl_Qu8mg5v7(
Removes duplicated vertices from a triangle mesh two vertices are considered the same if their distance is below
some threshold

Args:
    v : (\#v, 3)-shaped array of mesh vertex 3D positions
    f : (\#f, 3)-shaped array of face (triangle) indices

Returns:
    v_new : (\#v_new, 3)-shaped array of mesh vertex positions without unreferenced vertices
    f_new : (\#f_new, 3)-shaped array of mesh faces indexing into v_new
    correspondences_v : (\#v, 1)-shaped array of indices so that v_new = correspondences_v[svi]
    correspondences_f : (\#f, 1)-shaped array of indices so that f_new = correspondences_f[svj]

)igl_Qu8mg5v7";
npe_function(remove_unreferenced_mesh_vertices)
npe_doc(remove_unreferenced_mesh_vertices_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_begin_code()
{
    validate_mesh(v, f);
    npe_Matrix_v v_out;
    npe_Matrix_f f_out, corr_f, corr_v;
    igl::remove_unreferenced(v, f, v_out, f_out, corr_f, corr_v);

    return std::make_tuple(npe::move(v_out), npe::move(f_out), npe::move(corr_v), npe::move(corr_f));
}
npe_end_code()