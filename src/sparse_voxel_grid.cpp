#include <npe.h>
#include <tuple>
#include <algorithm>
#include <igl/marching_cubes.h>

#include "common/common.h"
#include "common/morton_code.h"


namespace {

template <typename Tfrom, typename Tto>
const Tto safe_cast(Tfrom val) {
    // WARNING: Only use this function to cast from more expressive types to less expressive types
    // e.g. double to int32_t is okay since all int32_ts can be encoded as doubles
    //      but float to int32_t is not okay since all int32_ts cannot be encoded as floats
    if (val > (Tfrom)(std::numeric_limits<Tto>::max() - 1) || val < (Tfrom)(std::numeric_limits<Tto>::min() + 1)) {
        throw pybind11::value_error("Invalid vertex leads to an overflow integer. Perhaps grid_size is too small.");
    }
    return (Tto) val;
}

template <typename T>
constexpr T safe_cast(T val) {
    return val;
}


template <typename T>
std::ptrdiff_t vector_binsearch(const T& needle, const std::vector<T>& haystack) {
    const auto index_it = std::lower_bound(haystack.begin(), haystack.end(), needle);
    if (index_it == haystack.end() || (needle < *index_it)) {
        return -1;
    }
    const std::ptrdiff_t idx = index_it - haystack.begin();
    return idx;
}


template <typename T>
void sort_deduplicate(T& dup_array) {
    std::sort(dup_array.begin(), dup_array.end(),
              [](typename T::value_type a, typename T::value_type b) { if (PyErr_CheckSignals() != 0) { throw pybind11::error_already_set(); } return a < b;});
    const auto dedup_end = std::unique(dup_array.begin(), dup_array.end());
    dup_array.resize(dedup_end - dup_array.begin());
}


template <typename T>
void morton_encode_eigen_matrix_and_sort(T& mat, std::vector<MortonCode64>& codes) {
    using scalar_t = typename T::Scalar;
    for (int i = 0; i < mat.rows(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t sx = safe_cast<scalar_t, int32_t>(mat(i, 0));
        const int32_t sy = safe_cast<scalar_t, int32_t>(mat(i, 1));
        const int32_t sz = safe_cast<scalar_t, int32_t>(mat(i, 2));
        codes.push_back(MortonCode64(sx, sy, sz));
    }
    std::sort(codes.begin(), codes.end());
}

template <typename GV>
void validate_sparse_voxel_grid(const GV& grid_coordinates) {
    if (grid_coordinates.rows() <= 0) {
        throw pybind11::value_error("Invalid grid_coordinates has zero rows!");
    }
    if(grid_coordinates.cols() != 3) {
        throw pybind11::value_error(std::string("Invalid shape for grid_coordinates must have shape (N, 3) but got (") +
                                    std::to_string(grid_coordinates.rows()) + ", " +
                                    std::to_string(grid_coordinates.cols()) + ")");
    }
}

}

npe_function(sparse_voxel_grid_from_pointcloud_internal)
npe_arg(points, dense_float, dense_double)
npe_arg(gx, double)
npe_arg(gy, double)
npe_arg(gz, double)
npe_arg(ox, double)
npe_arg(oy, double)
npe_arg(oz, double)
npe_arg(return_point_to_vox, bool)
npe_begin_code()
{
    if (points.rows() <= 0) {
        throw pybind11::value_error("Invalid points has zero rows!");
    }
    if (points.cols() != 3) {
        throw pybind11::value_error(std::string("Invalid shape for points must have shape (N, 3) but got (") +
                                    std::to_string(points.rows()) + ", " +
                                    std::to_string(points.cols()) + ")");
    }

    // Add morton codes for new vertices
    std::vector<MortonCode64> cube_bottom_corners;
    cube_bottom_corners.reserve(points.rows());

    std::vector<MortonCode64> point_to_voxel_codes;
    if (return_point_to_vox) {
        point_to_voxel_codes.reserve(points.rows());
    }

    std::vector<uint64_t> new_vertex_codes;
    new_vertex_codes.reserve(points.rows());
    for (int i = 0; i < points.rows(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t px = safe_cast<double, int32_t>(round((points(i, 0) - ox) / gx));
        const int32_t py = safe_cast<double, int32_t>(round((points(i, 1) - oy) / gy));
        const int32_t pz = safe_cast<double, int32_t>(round((points(i, 2) - oz) / gz));

        const MortonCode64 point_code(px, py, pz);
        cube_bottom_corners.push_back(point_code);

        if (return_point_to_vox) {
            point_to_voxel_codes.push_back(point_code);
        }

    }

    // Sort and deduplicate morton codes for cube bottom corners
    sort_deduplicate(cube_bottom_corners);

    using grid_coord_matrix_t = Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    grid_coord_matrix_t ret_grid_coords(cube_bottom_corners.size(), 3);
    // Generate grid positions
    for (int i = 0; i < cube_bottom_corners.size(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        int32_t gx, gy, gz;
        cube_bottom_corners[i].decode(gx, gy, gz);
        ret_grid_coords(i, 0) = gx;
        ret_grid_coords(i, 1) = gy;
        ret_grid_coords(i, 2) = gz;
    }

    using point_vox_index_t = Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    point_vox_index_t ret_point_vox_index;
    if (return_point_to_vox) {
        ret_point_vox_index.resize(points.rows(), 1);
        for (int i = 0; i < point_to_voxel_codes.size(); i++) {
            if (PyErr_CheckSignals() != 0) {
                throw pybind11::error_already_set();
            }
            std::ptrdiff_t grid_index = vector_binsearch(point_to_voxel_codes[i], cube_bottom_corners);
            ret_point_vox_index(i, 0) = grid_index;
        }
    }

    return std::make_tuple(npe::move(ret_grid_coords), npe::move(ret_point_vox_index));
}
npe_end_code()


npe_function(dilate_sparse_voxel_grid_internal)
npe_arg(grid_coordinates, dense_int, dense_long, dense_longlong)
npe_begin_code()
{
    validate_sparse_voxel_grid(grid_coordinates);

    std::vector<MortonCode64> voxel_codes;
    voxel_codes.reserve(grid_coordinates.rows());
    morton_encode_eigen_matrix_and_sort(grid_coordinates, voxel_codes);

    const std::array<MortonCode64, 26> neighbor_offsets = {
        MortonCode64(-1, -1, -1),
        MortonCode64(-1, -1, 0),
        MortonCode64(-1, -1, 1),
        MortonCode64(-1, 0, -1),
        MortonCode64(-1, 0, 0),
        MortonCode64(-1, 0, 1),
        MortonCode64(-1, 1, -1),
        MortonCode64(-1, 1, 0),
        MortonCode64(-1, 0, 1),

        MortonCode64(0, -1, -1),
        MortonCode64(0, -1, 0),
        MortonCode64(0, -1, 1),
        MortonCode64(0, 0, -1),
        MortonCode64(0, 0, 1),
        MortonCode64(0, 1, -1),
        MortonCode64(0, 1, 0),
        MortonCode64(0, 0, 1),

        MortonCode64(1, -1, -1),
        MortonCode64(1, -1, 0),
        MortonCode64(1, -1, 1),
        MortonCode64(1, 0, -1),
        MortonCode64(1, 0, 0),
        MortonCode64(1, 0, 1),
        MortonCode64(1, 1, -1),
        MortonCode64(1, 1, 0),
        MortonCode64(1, 0, 1),
    };

    std::vector<MortonCode64> new_cube_codes;
    new_cube_codes.reserve(grid_coordinates.rows());

    for (int i = 0; i < grid_coordinates.rows(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t vx = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 0));
        const int32_t vy = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 1));
        const int32_t vz = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 2));
        const MortonCode64 vox_code(vx, vy, vz);

        for (int o = 0; o < neighbor_offsets.size(); o += 1) {
            MortonCode64 test_code = vox_code + neighbor_offsets[o];
            std::ptrdiff_t idx = vector_binsearch(test_code, voxel_codes);
            if (idx < 0) {
                new_cube_codes.push_back(test_code);
            }
        }
    }

    sort_deduplicate(new_cube_codes);

    using grid_coord_matrix_t = Eigen::Matrix<npe_Scalar_grid_coordinates, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    grid_coord_matrix_t ret(new_cube_codes.size() + grid_coordinates.rows(), 3);
    ret.block(0, 0, grid_coordinates.rows(), grid_coordinates.cols()) = grid_coordinates;
    for (int i = 0; i < new_cube_codes.size(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        int32_t vx, vy, vz;
        new_cube_codes[i].decode(vx, vy, vz);

        std::ptrdiff_t offset = i + grid_coordinates.rows();
        ret(offset, 0) = (npe_Scalar_grid_coordinates) vx;
        ret(offset, 1) = (npe_Scalar_grid_coordinates) vy;
        ret(offset, 2) = (npe_Scalar_grid_coordinates) vz;
    }

    return npe::move(ret);
}
npe_end_code()


npe_function(erode_sparse_voxel_grid_internal)
npe_arg(grid_coordinates, dense_int, dense_long, dense_longlong)
npe_begin_code()
{
    validate_sparse_voxel_grid(grid_coordinates);

    std::vector<MortonCode64> voxel_codes;
    voxel_codes.reserve(grid_coordinates.rows());
    morton_encode_eigen_matrix_and_sort(grid_coordinates, voxel_codes);

    const std::array<MortonCode64, 26> neighbor_offsets = {
        MortonCode64(-1, -1, -1),
        MortonCode64(-1, -1, 0),
        MortonCode64(-1, -1, 1),
        MortonCode64(-1, 0, -1),
        MortonCode64(-1, 0, 0),
        MortonCode64(-1, 0, 1),
        MortonCode64(-1, 1, -1),
        MortonCode64(-1, 1, 0),
        MortonCode64(-1, 0, 1),

        MortonCode64(0, -1, -1),
        MortonCode64(0, -1, 0),
        MortonCode64(0, -1, 1),
        MortonCode64(0, 0, -1),
        MortonCode64(0, 0, 1),
        MortonCode64(0, 1, -1),
        MortonCode64(0, 1, 0),
        MortonCode64(0, 0, 1),

        MortonCode64(1, -1, -1),
        MortonCode64(1, -1, 0),
        MortonCode64(1, -1, 1),
        MortonCode64(1, 0, -1),
        MortonCode64(1, 0, 0),
        MortonCode64(1, 0, 1),
        MortonCode64(1, 1, -1),
        MortonCode64(1, 1, 0),
        MortonCode64(1, 0, 1),
    };

    std::vector<std::ptrdiff_t> keep_cubes;
    keep_cubes.reserve(grid_coordinates.rows());

    for (int i = 0; i < grid_coordinates.rows(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t vx = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 0));
        const int32_t vy = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 1));
        const int32_t vz = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 2));
        const MortonCode64 vox_code(vx, vy, vz);

        bool delete_cube = false;
        for (int o = 0; o < neighbor_offsets.size(); o += 1) {
            MortonCode64 test_code = vox_code + neighbor_offsets[o];
            std::ptrdiff_t idx = vector_binsearch(test_code, voxel_codes);
            if (idx < 0) {
                delete_cube = true;
                break;
            }
        }
        if (!delete_cube) {
            keep_cubes.push_back(i);
        }
    }

    using grid_coord_matrix_t = Eigen::Matrix<npe_Scalar_grid_coordinates, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    grid_coord_matrix_t ret(keep_cubes.size(), 3);
    for (int i = 0; i < keep_cubes.size(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        ret(i, 0) = grid_coordinates(keep_cubes[i], 0);
        ret(i, 1) = grid_coordinates(keep_cubes[i], 1);
        ret(i, 2) = grid_coordinates(keep_cubes[i], 2);
    }

    return npe::move(ret);
}
npe_end_code()


npe_function(subdivide_sparse_voxel_grid_internal)
npe_arg(grid_coordinates, dense_int, dense_long, dense_longlong)
npe_begin_code()
{
    validate_sparse_voxel_grid(grid_coordinates);

    std::vector<MortonCode64> new_vox_codes;
    new_vox_codes.reserve(8*grid_coordinates.rows());

    const std::array<MortonCode64, 8> neighbor_offsets = {
        MortonCode64(1, 1, 1),
        MortonCode64(1, 1, -1),
        MortonCode64(1, -1, 1),
        MortonCode64(1, -1, -1),
        MortonCode64(-1, 1, 1),
        MortonCode64(-1, 1, -1),
        MortonCode64(-1, -1, 1),
        MortonCode64(-1, -1, -1),
    };

    for (int i = 0; i < grid_coordinates.rows(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }

        const int32_t sx = safe_cast<npe_Scalar_grid_coordinates, int32_t>(2 * grid_coordinates(i, 0));
        const int32_t sy = safe_cast<npe_Scalar_grid_coordinates, int32_t>(2 * grid_coordinates(i, 1));
        const int32_t sz = safe_cast<npe_Scalar_grid_coordinates, int32_t>(2 * grid_coordinates(i, 2));
        const MortonCode64 base_code(sx, sy, sz);

        for (int o = 0; o < neighbor_offsets.size(); o += 1) {
            const MortonCode64 new_code = base_code + neighbor_offsets[o];
            new_vox_codes.push_back(new_code);
        }
    }

    sort_deduplicate(new_vox_codes);

    EigenDenseLike<npe_Matrix_grid_coordinates> ret(new_vox_codes.size(), 3);
    for (int i = 0; i < ret.rows(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        int32_t vx, vy, vz;
        new_vox_codes[i].decode(vx, vy, vz);
        ret(i, 0) = vx;
        ret(i, 1) = vy;
        ret(i, 2) = vz;
    }

    return npe::move(ret);
    //    print("i love you")
}
npe_end_code()


npe_function(sparse_voxel_grid_to_hex_mesh_internal)
npe_arg(grid_coordinates, dense_int, dense_long, dense_longlong)
npe_arg(gx, double)  // voxel size along x
npe_arg(gy, double)  // voxel size along y
npe_arg(gz, double)  // voxel size along z
npe_arg(ox, double)  // origin x coordinate
npe_arg(oy, double)  // origin y coordinate
npe_arg(oz, double)  // origin z coordinate
npe_begin_code()
{
    const std::array<MortonCode64, 8> vertex_offsets = {
        MortonCode64(0, 0, 0),
        MortonCode64(1, 0, 0),
        MortonCode64(1, 1, 0),
        MortonCode64(0, 1, 0),
        MortonCode64(0, 0, 1),
        MortonCode64(1, 0, 1),
        MortonCode64(1, 1, 1),
        MortonCode64(0, 1, 1),
    };

    std::vector<MortonCode64> vertex_codes;
    vertex_codes.reserve(8*grid_coordinates.rows());

    for (int i = 0; i < grid_coordinates.rows(); i++) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t vx = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 0));
        const int32_t vy = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 1));
        const int32_t vz = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 2));
        const MortonCode64 corner_code(vx, vy, vz);

        for (int o = 0; o < vertex_offsets.size(); o += 1) {
            vertex_codes.push_back(corner_code + vertex_offsets[o]);
        }
    }

    sort_deduplicate(vertex_codes);

    Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ret_cubes(grid_coordinates.rows(), 8);
    for (int i = 0; i < grid_coordinates.rows(); i+= 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        // No need to safe_cast here since safe_cast succeeded above for the same inputs
        const int32_t vx = grid_coordinates(i, 0);
        const int32_t vy = grid_coordinates(i, 1);
        const int32_t vz = grid_coordinates(i, 2);
        const MortonCode64 corner_code(vx, vy, vz);

        for (int o = 0; o < vertex_offsets.size(); o += 1) {
            const std::ptrdiff_t idx = vector_binsearch(corner_code + vertex_offsets[o], vertex_codes);
            if (idx < 0) {
                throw std::runtime_error("Internal error. Neighbor lookup failed. This shouldn't happen! Please file an issue.");
            }
            ret_cubes(i, o) = idx;
        }
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ret_vertices(vertex_codes.size(), 3);
    for (int i = 0; i < vertex_codes.size(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        int32_t vx, vy, vz;
        vertex_codes[i].decode(vx, vy, vz);

        ret_vertices(i, 0) = ((double) vx) * gx + ox - 0.5 * gx;
        ret_vertices(i, 1) = ((double) vy) * gy + oy - 0.5 * gx;
        ret_vertices(i, 2) = ((double) vz) * gz + oz - 0.5 * gx;
    }

    return std::make_tuple(npe::move(ret_vertices), npe::move(ret_cubes));
}
npe_end_code()


const char* sparse_voxel_grid_boundary_doc = R"Qu8mg5v7(
Find the indices of the voxels which lie on the boundary of a sparse voxel grid

Args:
    grid_coordinates : An (n, 3) shaped integer array of voxels coordinates in the sparse voxel grid

Returns:
    boundary_voxels : An (m,) shaped array of indices into grid_coordinates encoding which voxels lie on the boundary

)Qu8mg5v7";
npe_function(sparse_voxel_grid_boundary)
npe_arg(grid_coordinates, dense_int, dense_long, dense_longlong)
// npe_doc(sparse_voxel_grid_boundary_doc)
npe_begin_code()
{
    validate_sparse_voxel_grid(grid_coordinates);

    std::vector<MortonCode64> voxel_codes;
    voxel_codes.reserve(grid_coordinates.rows());
    morton_encode_eigen_matrix_and_sort(grid_coordinates, voxel_codes);

    const std::array<MortonCode64, 6> neighbor_offsets = {
        MortonCode64(1, 0, 0),
        MortonCode64(-1, 0, 0),
        MortonCode64(0, 1, 0),
        MortonCode64(0, -1, 0),
        MortonCode64(0, 0, 1),
        MortonCode64(0, 0, -1)
    };

    Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ret_indices(grid_coordinates.rows(), 1);
    int num_boundary_voxels = 0;

    for (int i = 0; i < grid_coordinates.rows(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        const int32_t vx = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 0));
        const int32_t vy = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 1));
        const int32_t vz = safe_cast<npe_Scalar_grid_coordinates, int32_t>(grid_coordinates(i, 2));
        const MortonCode64 vox_code(vx, vy, vz);

        for (int o = 0; o < neighbor_offsets.size(); o += 1) {
            const MortonCode64 test_code = vox_code + neighbor_offsets[o];
            int32_t sx, sy, sz;
            test_code.decode(sx, sy, sz);
            std::ptrdiff_t idx = vector_binsearch(test_code, voxel_codes);
            if (idx < 0) { // Missing a neighbor, is a boundary voxel
                ret_indices(num_boundary_voxels, 0) = i;
                num_boundary_voxels += 1;
                break;
            }
        }
    }

    ret_indices.conservativeResize(num_boundary_voxels, 1);

    return npe::move(ret_indices);
}
npe_end_code()


