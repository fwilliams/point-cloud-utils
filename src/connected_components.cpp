#include <npe.h>
#include <igl/adjacency_matrix.h>
#include <igl/connected_components.h>
#include <numeric>

#include "common.h"

const char* connected_components_doc = R"igl_Qu8mg5v7(
Determine the connected components of a mesh

Parameters
----------
v: (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
f: (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)

Returns
-------
A tuple (cv, nv, cf, nf) where:
  - cv is a (#vertices,)-shaped array of integer indexes (starting from 0) indicating which
    component each vertex belongs to. i.e. cv[i] is the component of the vertex v[i].
  - nv is the number of vertices in each connected component. i.e. nv[j] is the number of vertices
    in component j
  - cf is a (#faces,)-shaped array of integer indexes (starting from 0) indicating which
    component each face belongs to. i.e. cf[i] is the component of the face f[i].
  - nf is the number of faces in each connected component. i.e. nf[j] is the number of faces in
    component j
)igl_Qu8mg5v7";
npe_function(connected_components)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_begin_code()
{
    Eigen::SparseMatrix<npe_Scalar_f> A;
    igl::adjacency_matrix(f, A);

    npe_Matrix_f c_v, c_f, v_counts, f_counts;
    igl::connected_components(A, c_v, v_counts);
    c_f.resize(f.rows(), 1);
    f_counts = npe_Matrix_f::Zero(v_counts.rows(), 1);
    for (int i = 0; i < c_f.rows(); i++) {
        npe_Scalar_f face_component = c_v(f(i, 0), 0);
        c_f(i, 0) = face_component;
        f_counts(face_component, 0) += 1;
    }
    pybind11::print("f.shape " + std::to_string(f.rows()) + " " + std::to_string(f.cols()));
    pybind11::print("A.shape " + std::to_string(A.rows()) + " " + std::to_string(A.cols()));
    return std::make_tuple(npe::move(c_v), npe::move(v_counts), npe::move(c_f), npe::move(f_counts));
}
npe_end_code()
