#include <npe.h>

#include <igl/remove_duplicate_vertices.h>

#include "common/common.h"



const char* remove_mesh_vertices_doc = R"igl_Qu8mg5v7(
Removes vertices specified by a mask from a triangle mesh and updates the face indices accordingly.

Args:
    v : \#v by 3 Matrix of mesh vertex 3D positions
    f : \#f by 3 Matrix of face (triangle) indices
    mask: A boolean mask of shape (\#v,) indicating which vertices to keep (True) or remove (False)

Returns:
    v_out : (\#v, 3)-shaped array of mesh vertices with the mask applied
    f_out : (\#f, 3)-shaped array of mesh faces corresponding to mesh with removed vertices
See also:
    deduplicate_point_cloud
)igl_Qu8mg5v7";
npe_function(remove_mesh_vertices)
npe_doc(remove_mesh_vertices_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64)
npe_arg(mask, dense_bool)
npe_begin_code()
{
    validate_mesh(v, f);
    if (mask.rows() != v.rows()) {
        throw std::invalid_argument("mask should have the same number of rows as v");
    }
    if (mask.cols() != 1) {
        throw std::invalid_argument("mask should have only one column");
    }

    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v_out(v.rows(), v.cols());
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f_out(f.rows(), v.cols());
    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> unmasked_to_masked(v.rows());

    uint32_t num_unmasked_vertices = 0;
    for (int i = 0; i < mask.size(); ++i) {
        if (mask(i, 0)) {
            unmasked_to_masked(i, 0) = num_unmasked_vertices;
            v_out.row(unmasked_to_masked(i, 0)) = v.row(i);
            num_unmasked_vertices += 1;
        } else {
            unmasked_to_masked(i, 0) = -1;
        }
    }

    uint32_t num_unmasked_faces = 0;
    for (int i = 0; i < f.rows(); i += 1) {
        bool keep = true;
        for (int j = 0; j < f.cols(); j += 1) {
            if (unmasked_to_masked(f(i, j), 0) == -1) {
                keep = false;
                break;
            }
        }
        if (keep) {
            for (int j = 0; j < f.cols(); j += 1) {
                f_out(num_unmasked_faces, j) = unmasked_to_masked(f(i, j), 0);
            }
            num_unmasked_faces += 1;
        }
    }

    v_out.conservativeResize(num_unmasked_vertices, Eigen::NoChange);
    f_out.conservativeResize(num_unmasked_faces, Eigen::NoChange);
    return std::make_tuple(npe::move(v_out), npe::move(f_out));
}
npe_end_code()