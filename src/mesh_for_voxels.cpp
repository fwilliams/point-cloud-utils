
#include <npe.h>

#include "common/common.h"

#include <Eigen/Core>


namespace {

    template <typename MatrixIJK>
    void generate_cube_mesh(Eigen::Vector3d vox_origin, Eigen::Vector3d vox_size,
                            const MatrixIJK& in_ijk,
                            Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_v,
                            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_f) {

        std::array<int, 36> indices = {
            //Top
            2, 6, 7,
            2, 7, 3,

            //Bottom
            0, 4, 5,
            0, 5, 1,

            //Left
            0, 2, 6,
            0, 6, 4,

            //Right
            1, 3, 7,
            1, 7, 5,

            //Front
            0, 2, 3,
            0, 3, 1,

            //Back
            4, 6, 7,
            4, 7, 5
        };


        std::array<float, 24> vertices = {
            -1.0, -1.0,  1.0, //0
             1.0, -1.0,  1.0, //1
            -1.0,  1.0,  1.0, //2
             1.0,  1.0,  1.0, //3
            -1.0, -1.0, -1.0, //4
             1.0, -1.0, -1.0, //5
            -1.0,  1.0, -1.0, //6
             1.0,  1.0, -1.0  //7
        };

        out_v.resize(8 * in_ijk.rows(), 3);
        out_f.resize(12 * in_ijk.rows(), 3);

        for (int i = 0; i < in_ijk.rows(); i += 1) {
            for (int vi = 0; vi < 8; vi += 1) {
                Eigen::Vector3d v_world((vertices[vi*3+0] * 0.5 + (float) in_ijk(i, 0) + 0.5)* vox_size[0],
                                        (vertices[vi*3+1] * 0.5 + (float) in_ijk(i, 1) + 0.5)* vox_size[1],
                                        (vertices[vi*3+2] * 0.5 + (float) in_ijk(i, 2) + 0.5)* vox_size[0]);
                v_world += vox_origin;
                v_world = v_world.eval();
                out_v(i * 8 + vi, 0) = v_world[0];
                out_v(i * 8 + vi, 1) = v_world[1];
                out_v(i * 8 + vi, 2) = v_world[2];
            }

            for (int fi = 0; fi < 12; fi += 1) {
                out_f(i * 12 + fi, 0) = indices[fi*3+0] + i * 8;
                out_f(i * 12 + fi, 1) = indices[fi*3+1] + i * 8;
                out_f(i * 12 + fi, 2) = indices[fi*3+2] + i * 8;
            }
        }
    }
}


npe_function(_voxel_mesh_internal)
npe_arg(ijk, dense_int, dense_long, dense_longlong)
npe_arg(vox_origin, dense_double)
npe_arg(vox_size, dense_double)
npe_begin_code()
    using MatrixI = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MatrixF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    validate_point_cloud(ijk, false /*allow_0*/);
    if (vox_origin.size() != 3) {
        throw pybind11::value_error("Invalid shape");
    }

    if (vox_size.size() != 3) {
        throw pybind11::value_error("Invalid shape");
    }

    Eigen::Vector3d vsize = vox_size;
    Eigen::Vector3d vorgn = vox_origin;

    MatrixF geom_vertices;
    MatrixI geom_faces;

    generate_cube_mesh(vorgn, vsize, ijk, geom_vertices, geom_faces);

    return std::make_tuple(npe::move(geom_vertices), npe::move(geom_faces));
npe_end_code()