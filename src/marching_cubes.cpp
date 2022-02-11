#include <npe.h>
#include <igl/marching_cubes.h>
#include <tuple>
#include <numeric>

#include "common/common.h"

//npe_function(subdivide_sparse_voxel_grid)
//npe_arg(grid_coordinates, dense_int, dense_long_dense_longlong)
//npe_arg(cube_indices, dense_int, dense_long, dense_longlong)
//npe_arg(num_subdivs, int)
//npe_begin_code()
//{
//
//    if (cube_indices.rows() <= 0) {
//        throw pybind11::value_error("Invalid cube_indices has zero rows!");
//    }
//    if (cube_vertices.rows() <= 0) {
//        throw pybind11::value_error("Invalid cube_vertices has zero rows!");
//    }
//    if(cube_vertices.rows() <= 0 || cube_vertices.cols() != 3) {
//        throw pybind11::value_error(std::string("Invalid shape for cube_vertices must have shape (N, 3) but got (") +
//                                    std::to_string(grid_vertices.rows()) + ", " +
//                                    std::to_string(grid_vertices.cols()) + ")");
//    }
//
//    using index_t = npe_Scalar_cube_indices;
//    for (int c_i = 0; c_i < cube_indices.rows(); c_i++) {
//        index_t c000 = cube_indices(c_i, 0),
//                c100 = cube_indices(c_i, 1),
//                c110 = cube_indices(c_i, 2),
//                c010 = cube_indices(c_i, 3),
//                c001 = cube_indices(c_i, 4),
//                c101 = cube_indices(c_i, 5),
//                c111 = cube_indices(c_i, 6),
//                c011 = cube_indices(c_i, 7);
//
//    }
//}
//npe_end_code()


npe_function(sparse_marching_cubes)
npe_arg(grid_vertices, dense_float, dense_double)
npe_arg(grid_scalars, npe_matches(grid_vertices))
npe_arg(cube_indices, dense_int, dense_longlong)
npe_arg(isovalue, double)
npe_begin_code()
{
    if (grid_vertices.rows() <= 0) {
        throw pybind11::value_error("Invalid grid_vertices has zero rows!");
    }
    if(grid_vertices.rows() <= 0 || grid_vertices.cols() != 3) {
        throw pybind11::value_error(std::string("Invalid shape for grid_vertices must have shape (N, 3) but got (") +
                                    std::to_string(grid_vertices.rows()) + ", " +
                                    std::to_string(grid_vertices.cols()) + ")");
    }
    if (grid_scalars.cols() != 1) {
        throw pybind11::value_error(std::string("Invalid shape for grid_scalars must have shape (N,) or (N, 1) but got (") +
                                std::to_string(grid_scalars.rows()) + ", " +
                                std::to_string(grid_scalars.cols()) + ")");
    }
    if (grid_scalars.rows() != grid_vertices.rows()) {
        throw pybind11::value_error(std::string("grid_vertices and grid_scalars must have the same number of rows but got ") +
                                    std::string("grid_vertices.shape = (") +
                                    std::to_string(grid_vertices.rows()) + ", " +
                                    std::to_string(grid_vertices.cols()) + std::string("), and grid_scalars.shape = (") +
                                    std::to_string(grid_vertices.rows()) + ", " +
                                    std::to_string(grid_vertices.cols()) + ")");
    }
    if (cube_indices.rows() == 0) {
        throw pybind11::value_error("Invalid cube_indices has zero rows!");
    }
    if (cube_indices.cols() != 8) {
            throw pybind11::value_error(std::string("Invalid shape for cube_indices must have shape (N, 8) but got (") +
                                std::to_string(cube_indices.rows()) + ", " +
                                std::to_string(cube_indices.cols()) + ")");
    }

    EigenDenseLike<npe_Matrix_grid_vertices> v;
    EigenDenseLike<npe_Matrix_cube_indices> f;
    npe_Matrix_grid_scalars gs_copy = grid_scalars;
    igl::marching_cubes(gs_copy, grid_vertices, cube_indices, (npe_Scalar_grid_scalars) isovalue, v, f);

    return std::make_tuple(npe::move(v), npe::move(f));

}
npe_end_code()