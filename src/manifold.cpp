#include <npe.h>
#include "Model_OBJ.h"
#include "glm/glm.hpp"
#include "common/common.h"


const char* make_mesh_watertight_doc = R"igl_Qu8mg5v7(
Convert a mesh into a watertight manifold

Args:
    v : (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
    f : (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)
    resolution : A resolution parameter for the algorithm

Returns:
    vw : a (#vw, 3)-shaped array of vertices (one per row)
    fw : a (#fw, 3)-shaped array of face indices into vw (one face per row)
)igl_Qu8mg5v7";
npe_function(make_mesh_watertight)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_default_arg(resolution, double, 20000)
npe_default_arg(seed, int, -1)
npe_doc(make_mesh_watertight_doc)
npe_begin_code()
{
    validate_mesh(v, f);

    // Make the example deterministic
    if (seed > 0) {
        srand((unsigned int) seed);
    }
    Model_OBJ model;
    model.vertices.resize(v.rows());
    model.face_indices.resize(f.rows());
    for (int i = 0; i < v.rows(); i += 1) {
        model.vertices[i] = glm::dvec3(v(i, 0), v(i, 1), v(i, 2));
    }
    for (int i = 0; i < f.rows(); i += 1) {
        model.face_indices[i] = glm::ivec3(f(i, 0), f(i, 1), f(i, 2));
    }
    model.Process_Manifold(resolution);
    EigenDense<npe_Scalar_v> v_watertight(model.vertices.size(), 3);
    EigenDense<npe_Scalar_f> f_watertight(model.face_indices.size(), 3);
    for (int i = 0; i < model.vertices.size(); i += 1) {
        for (int j = 0; j < 3; j += 1) {
            v_watertight(i, j) = (npe_Scalar_v) model.vertices[i][j];
        }
    }
    for (int i = 0; i < model.face_indices.size(); i += 1) {
        for (int j = 0; j < 3; j += 1) {
            f_watertight(i, j) = (npe_Scalar_f) model.face_indices[i][j];
        }
    }

    return std::make_tuple(npe::move(v_watertight), npe::move(f_watertight));
}
npe_end_code()
