#include <npe.h>
#include <igl/decimate.h>


#include "common/common.h"

const char* decimate_triangle_mesh_doc = R"Qu8mg5v7(
Decimate a (manifold) triangle mesh by collapsing edges

Args:
    v : (#v, 3)-shaped array of mesh vertex positions
    f : (#f, 3)-shaped array of triangle face indices
    max_faces : The maximum number of faces in the decimated output mesh (must be between 0 and #f)
    decimation_heuristic : Which decimation heuristic to use. Currently only supports "shortest_edge".

Returns:
    v_out : (#v_out, 3)-shaped array of vertex positions for the decimated mesh
    f_out : (#f_out, 3)-shaped array of triangle face indices for the decimated mesh
    v_correspondences : (#v_out, 1)-shaped array where v_correspondences[i] is the index of the vertex in v which generated v_out[i]
    f_correspondences : (#f_out, 1)-shaped array where f_correspondences[i] is the index of the face in f which generated f_out[i]

)Qu8mg5v7";
npe_function(decimate_triangle_mesh)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(max_faces, int)
npe_default_arg(decimation_heuristic, std::string, "shortest_edge")
npe_doc(decimate_triangle_mesh_doc)
npe_begin_code()
{
    validate_mesh(v, f);
    if (max_faces <= 0) {
        throw pybind11::value_error("max_faces must be greater than 0, got " + std::to_string(max_faces) + ".");
    }
    if (max_faces > f.rows()) {
        throw pybind11::value_error("max_faces must be less than or equal to the number of input faces (" +
                                    std::to_string(f.rows()) + "), got " + std::to_string(max_faces) + ".");
    }
    npe_Matrix_v out_v;
    npe_Matrix_f out_f;
    npe_Matrix_f out_face_correspondences, out_vertex_correspondences;

    if (decimation_heuristic == "shortest_edge") {
        Eigen::MatrixXd v_copy = v.template cast<double>();
        Eigen::MatrixXi f_copy = f.template cast<int>();
        Eigen::MatrixXd out_v_copy;
        Eigen::MatrixXi out_f_copy;
        Eigen::VectorXi out_f_corr_copy, out_v_corr_copy;
        igl::decimate(v_copy, f_copy, (size_t) max_faces, out_v_copy, out_f_copy, out_f_corr_copy, out_v_corr_copy);

        if (out_v_copy.rows() == 0) {
            throw pybind11::value_error("decimate_triangle_mesh_doc produced an empty mesh. "
                                        "This can happen if the input is non-manifold.");
        }
        out_v = out_v_copy.template cast<npe_Scalar_v>();
        out_f = out_f_copy.template cast<npe_Scalar_f>();
        out_face_correspondences = out_f_corr_copy.template cast<npe_Scalar_f>();
        out_vertex_correspondences = out_v_corr_copy.template cast<npe_Scalar_f>();

        return std::make_tuple(npe::move(out_v), npe::move(out_f), npe::move(out_vertex_correspondences), npe::move(out_face_correspondences));
    } else {
        throw pybind11::value_error("Invalid decimation heuristic, must be 'shortest_edge'. Others will be implemented"
                                    " in future versions of Point Cloud Utils.");
    }
}
npe_end_code()