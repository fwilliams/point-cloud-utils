#include <npe.h>

#include <igl/remove_duplicate_vertices.h>

#include "common.h"

const char* remove_duplicate_points_doc = R"igl_Qu8mg5v7(
Removes duplicated points from a point cloud where two points are considered the same if their distance is below
some threshold

Parameters
----------
x : #x by 3 Matrix of 3D positions
epsilon: threshold below which two points are considered equal
return_index: If true, return indices to map between input and output

Returns
-------
x_new : #x_new x 3 Point cloud with duplicates removed
if return indices is set, this function also returns:
    svi : #x x 1 indices so that x_new = x[svi]
    svj : #x_new x 1 indices so that x = x_new[svj]

See also
--------
remove_duplicate_mesh_vertices
)igl_Qu8mg5v7";
npe_function(remove_duplicate_points)
npe_doc(remove_duplicate_points_doc)
npe_arg(points, dense_float, dense_double)
npe_arg(epsilon, double)
npe_default_arg(return_index, bool, false)
npe_begin_code()
{
    validate_point_cloud(points);
    Eigen::Matrix<npe_Scalar_points, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x_copy = points;
    Eigen::Matrix<npe_Scalar_points, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x_out;
    Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> svi;
    Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> svj;

    igl::remove_duplicate_vertices(x_copy, epsilon, x_out, svi, svj);

    if (return_index) {
        return pybind11::cast(std::make_tuple(npe::move(x_out), npe::move(svi), npe::move(svj)));
    } else {
        return npe::move(x_out);
    }
}
npe_end_code()



const char* remove_duplicate_mesh_vertices_doc = R"igl_Qu8mg5v7(
Removes duplicated vertices from a triangle mesh two vertices are considered the same if their distance is below
some threshold

Parameters
----------
v : #v by 3 Matrix of mesh vertex 3D positions
f : #f by 3 Matrix of face (triangle) indices
epsilon: threshold below which two points are considered equal
return_index: If true, return indices to map between input and output

Returns
-------
x_new : #x_new x 3 Point cloud with duplicates removed
if return indices is set, this function also returns:
    svi : #x x 1 indices so that x_new = x[svi]
    svj : #x_new x 1 indices so that x = x_new[svj]
)igl_Qu8mg5v7";
npe_function(remove_duplicate_mesh_vertices)
npe_doc(remove_duplicate_mesh_vertices_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_arg(epsilon, double)
npe_default_arg(return_index, bool, false)
npe_begin_code()
{
    validate_mesh(v, f);
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v_copy = v;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f_copy = f;
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v_out;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f_out;
    Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> svi;
    Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> svj;

    igl::remove_duplicate_vertices(v_copy, f_copy, epsilon, v_out, f_out, svi, svj);

    if (return_index) {
        return pybind11::cast(std::make_tuple(npe::move(v_out), npe::move(f_out), npe::move(svi), npe::move(svj)));
    } else {
        return pybind11::cast(std::make_tuple(npe::move(v_out), npe::move(f_out)));
    }
}
npe_end_code()