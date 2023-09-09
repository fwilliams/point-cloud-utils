#include <npe.h>
#include <igl/marching_cubes.h>
#include <tuple>
#include <numeric>

#include "common/common.h"

const char* marching_cubes_sparse_voxel_grid_doc = R"Qu8mg5v7(
Run marching cubes on a hex mesh representing a sparse voxel grid (see also sparse_voxel_grid_to_hex_mesh)

Args:
    grid_scalars : An (n,) shaped array of scalar values at each hex mesh vertex
    grid_coordinates : An (n, 3) shaped array of hex mesh vertices
    cube_indices : An (m, 8) shaped array of indices into grid_coordinates where cube_indices[i, :] are the indices of the 8 points forming the i^th cube
                  Note the cube indices must be ordered as:
                  [[0, 0, 0],
                  [1, 0, 0],
                  [1, 1, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 1],
                  [1, 1, 1],
                  [0, 1, 1]]
                  where [i, j, k] indicates the offset along the (x, y, z) axes from the bottom, back, left
                  corner of the cube
    isovalue : Which level set to extract

Returns:
    v : A (#v, 3) array of triangle mesh vertices
    f : A (#f, 3) array of indices into v where f[i, :] are the indices into vof the 3 points forming the i^th triangle

See Also:
    sparse_voxel_grid_to_hex_mesh
)Qu8mg5v7";
npe_function(marching_cubes_sparse_voxel_grid)
npe_arg(grid_scalars, dense_float, dense_double)
npe_arg(grid_coordinates, dense_float, dense_double)
npe_arg(cube_indices, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(isovalue, double)
npe_doc(marching_cubes_sparse_voxel_grid_doc)
npe_begin_code()
{
    if (grid_coordinates.rows() <= 0) {
        throw pybind11::value_error("Invalid grid_coordinates has zero rows!");
    }
    if(grid_coordinates.rows() <= 0 || grid_coordinates.cols() != 3) {
        throw pybind11::value_error(std::string("Invalid shape for grid_coordinates must have shape (N, 3) but got (") +
                                    std::to_string(grid_coordinates.rows()) + ", " +
                                    std::to_string(grid_coordinates.cols()) + ")");
    }
    if (grid_scalars.cols() != 1) {
        throw pybind11::value_error(std::string("Invalid shape for grid_scalars must have shape (N,) or (N, 1) but got (") +
                                std::to_string(grid_scalars.rows()) + ", " +
                                std::to_string(grid_scalars.cols()) + ")");
    }
    if (grid_scalars.rows() != grid_coordinates.rows()) {
        throw pybind11::value_error(std::string("grid_coordinates and grid_scalars must have the same number of rows but got ") +
                                    std::string("grid_coordinates.shape = (") +
                                    std::to_string(grid_coordinates.rows()) + ", " +
                                    std::to_string(grid_coordinates.cols()) + std::string("), and grid_scalars.shape = (") +
                                    std::to_string(grid_coordinates.rows()) + ", " +
                                    std::to_string(grid_coordinates.cols()) + ")");
    }
    if (cube_indices.rows() == 0) {
        throw pybind11::value_error("Invalid cube_indices has zero rows!");
    }
    if (cube_indices.cols() != 8) {
            throw pybind11::value_error(std::string("Invalid shape for cube_indices must have shape (N, 8) but got (") +
                                std::to_string(cube_indices.rows()) + ", " +
                                std::to_string(cube_indices.cols()) + ")");
    }

    EigenDenseLike<npe_Matrix_grid_coordinates> v;
    EigenDenseLike<npe_Matrix_cube_indices> f;
    npe_Matrix_grid_coordinates gs_copy = grid_scalars.template cast<npe_Scalar_grid_coordinates>();
    npe_Matrix_grid_coordinates gc_copy = grid_coordinates.template cast<npe_Scalar_grid_coordinates>();
    igl::marching_cubes(gs_copy, grid_coordinates, cube_indices, (npe_Scalar_grid_scalars) isovalue, v, f);

    if (v.rows() == 0) {
        throw pybind11::value_error(std::string("No level set found for isovalue ") + std::to_string(isovalue));
    }
    return std::make_tuple(npe::move(v), npe::move(f));

}
npe_end_code()