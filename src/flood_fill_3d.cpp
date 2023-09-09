#include <npe.h>

#include <vector>
#include <queue>
#include <tuple>

namespace {

template <typename EigenMat, typename ScalarType>
void flood_fill(EigenMat& grid,
                int start_x, int start_y, int start_z,
                int size_x, int size_y, int size_z,
                ScalarType fill_val) {
    int start_offset = start_x * size_y * size_z + start_y * size_z + start_z;
    if (start_offset < 0 || start_offset >= grid.rows()) {
        throw pybind11::value_error("Invalid start coordinate must be in grid");
    }
    ScalarType start_val = grid(start_offset, 0);

    std::array<std::array<int, 3>, 8> directions = {{
        {1, 0, 0},
        {-1, 0, 0},
        {0, 1, 0},
        {0, -1, 0},
        {0, 0, 1},
        {0, 0, -1}
    }};
    std::queue<int> queue;
    queue.push({start_offset});

    while (!queue.empty()) {
        int offset = queue.front();

        int x = offset / (size_y * size_z);
        int y = (offset - (x * size_y * size_z)) / size_z;
        int z = offset - (x * size_y * size_z) - y * size_z;

        queue.pop();

        if (offset >= 0 && offset < grid.rows() && grid(offset, 0) == start_val) {
            grid(offset, 0) = fill_val;

            for (int i = 0; i < directions.size(); i += 1) {
                int dx = directions[i][0], dy = directions[i][1], dz = directions[i][2];
                int doffset = (x + dx) * size_y * size_z + (y + dy) * size_z + (z + dz);
                queue.push(doffset);
            }
        }
    }
}

}


npe_function(_flood_fill_3d_internal)
npe_arg(grid, dense_int32, dense_int64, dense_float, dense_double)
npe_arg(seed_x, int)
npe_arg(seed_y, int)
npe_arg(seed_z, int)
npe_arg(size_x, int)
npe_arg(size_y, int)
npe_arg(size_z, int)
npe_arg(flood_value, double)
npe_begin_code()

    if ((seed_x < 0 || seed_x >= size_x) || (seed_y < 0 || seed_y >= size_y) || (seed_z < 0 || seed_z >= size_z)) {
        throw pybind11::value_error("seed point must be inside grid");
    }
    if (grid.rows() != size_x * size_y * size_z) {
        throw pybind11::value_error("invalid size");
    }

    flood_fill(grid, seed_x, seed_y, seed_z, size_x, size_y, size_z, (npe_Scalar_grid) flood_value);

npe_end_code()