#include <cinttypes>
#include <cassert>
#include <algorithm>
#include <iostream>

#include <npe.h>

#include "common/common.h"
#include "common/morton_code.h"




const char* morton_add = R"Qu8mg5v7(
Add morton codes together (corresponding to adding the vectors they encode)

Parameters
----------
codes_1: an [n,] array of morton codes
codes_2 : an [n,] array of morton codes
num_threads : Number of threads to use. If set to -1, will use all available CPUs. If set to 0, will run in serial. Default is -1.
Returns
-------
an [n,] shaped array of added morton codes

)Qu8mg5v7";
npe_function(morton_add)
npe_arg(codes_1, dense_ulonglong)
npe_arg(codes_2, dense_ulonglong)
npe_default_arg(num_threads, int, -1)
npe_doc(morton_add)
npe_begin_code()
{
    if (codes_1.rows() <= 0) {
        throw pybind11::value_error("codes_1 must be an array of shape [n,] but got an empty array");
    }
    if (codes_1.cols() != 1) {
        throw pybind11::value_error("codes_1 must be an array of shape [n,] but got an invalid number of columns");
    }
    if (codes_2.rows() <= 0) {
        throw pybind11::value_error("codes_2 must be an array of shape [n,] but got an empty array");
    }
    if (codes_2.cols() != 1) {
        throw pybind11::value_error("codes_2 must be an array of shape [n,] but got an invalid number of columns");
    }
    if (codes_1.rows() != codes_1.rows()) {
        throw pybind11::value_error("codes_1 and codes_2  must have the same number of entries.");
    }


    Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> codes(codes_1.rows(), 1);

    const int MIN_PARALLEL_INPUT_SIZE = 10000;
    const bool run_parallel = codes_1.rows() >= MIN_PARALLEL_INPUT_SIZE && num_threads != 0;
    auto set_parallel = OmpSetParallelism(num_threads, run_parallel);

    bool threw_exception = false;
    #if defined(_OPENMP)
    #pragma omp parallel if (run_parallel)
    #endif
    {
        #if defined(_OPENMP)
        #pragma omp for
        #endif
        for(int i = 0; i < codes_1.rows(); i += 1) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                #else
                    threw_exception = true;
                    break;
                #endif
            }

            MortonCode64 c1 = MortonCode64(codes_1(i, 0));
            MortonCode64 c2 = MortonCode64(codes_2(i, 0));
            codes(i, 0) = (c1 + c2).get_data();
        }
    }
    if (threw_exception) {
        throw pybind11::error_already_set();
    }

    return npe::move(codes);
}
npe_end_code()


const char* morton_subtract = R"Qu8mg5v7(
Subtract morton codes from each other (corresponding to adding the vectors they encode)

Parameters
----------
codes_1: an [n,] array of morton codes
codes_2 : an [n,] array of morton codes
num_threads : Number of threads to use. If set to -1, will use all available CPUs. If set to 0, will run in serial. Default is -1.
Returns
-------
an [n,] shaped array of added morton codes

)Qu8mg5v7";
npe_function(morton_subtract)
npe_arg(codes_1, dense_ulonglong)
npe_arg(codes_2, dense_ulonglong)
npe_default_arg(num_threads, int, -1)
npe_doc(morton_subtract)
npe_begin_code()
{
    if (codes_1.rows() <= 0) {
        throw pybind11::value_error("codes_1 must be an array of shape [n,] but got an empty array");
    }
    if (codes_1.cols() != 1) {
        throw pybind11::value_error("codes_1 must be an array of shape [n,] but got an invalid number of columns");
    }
    if (codes_2.rows() <= 0) {
        throw pybind11::value_error("codes_2 must be an array of shape [n,] but got an empty array");
    }
    if (codes_2.cols() != 1) {
        throw pybind11::value_error("codes_2 must be an array of shape [n,] but got an invalid number of columns");
    }
    if (codes_1.rows() != codes_1.rows()) {
        throw pybind11::value_error("codes_1 and codes_2  must have the same number of entries.");
    }


    Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> codes(codes_1.rows(), 1);

    const int MIN_PARALLEL_INPUT_SIZE = 10000;
    const bool run_parallel = codes_1.rows() >= MIN_PARALLEL_INPUT_SIZE && num_threads != 0;
    auto set_parallel = OmpSetParallelism(num_threads, run_parallel);

    bool threw_exception = false;
    #if defined(_OPENMP)
    #pragma omp parallel if (run_parallel)
    #endif
    {
        #if defined(_OPENMP)
        #pragma omp for
        #endif
        for(int i = 0; i < codes_1.rows(); i += 1) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                #else
                    threw_exception = true;
                    break;
                #endif
            }

            MortonCode64 c1 = MortonCode64(codes_1(i, 0));
            MortonCode64 c2 = MortonCode64(codes_2(i, 0));
            codes(i, 0) = (c1 - c2).get_data();
        }
    }
    if (threw_exception) {
        throw pybind11::error_already_set();
    }

    return npe::move(codes);
}
npe_end_code()


const char* morton_encode = R"Qu8mg5v7(
Encode n 3D points using Morton coding, possibly sorting them

Parameters
----------
pts: an [n, 3] array of 3D points
num_threads : Number of threads to use. If set to -1, will use all available CPUs. If set to 0, will run in serial. Default is -1.
Returns
-------
an [n] shaped array of morton encoded points

)Qu8mg5v7";
npe_function(morton_encode)
npe_arg(pts, dense_int, dense_long, dense_longlong)
npe_default_arg(num_threads, int, -1)
npe_doc(morton_encode)
npe_begin_code()
{
    if (pts.rows() <= 0) {
        throw pybind11::value_error("pts must be an array of shape [n, 3] but got an empty array");
    }
    if (pts.cols() != 3) {
        throw pybind11::value_error("pts must be an array of shape [n, 3] but got an invalid number of columns");
    }

    Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> codes(pts.rows(), 1);

    const int MIN_PARALLEL_INPUT_SIZE = 10000;
    const bool run_parallel = pts.rows() >= MIN_PARALLEL_INPUT_SIZE && num_threads != 0;
    auto set_parallel = OmpSetParallelism(num_threads, run_parallel);

    bool threw_exception = false;
    #if defined(_OPENMP)
    #pragma omp parallel if (run_parallel)
    #endif
    {
        #if defined(_OPENMP)
        #pragma omp for
        #endif
        for(int i = 0; i < pts.rows(); i += 1) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                #else
                    threw_exception = true;
                    break;
                #endif
            }

            int32_t px = pts(i, 0), py = pts(i, 1), pz = pts(i, 2);
            MortonCode64 code(px, py, pz);
            codes[i] = code.get_data();
        }
    }
    if (threw_exception) {
        throw pybind11::error_already_set();
    }

    return npe::move(codes);
}
npe_end_code()



const char* morton_decode = R"Qu8mg5v7(
Decode n points along a Morton curve into 3D points

Parameters
----------
codes: an [n] shaped array of Morton codes

Returns
-------
an [n, 3] shaped array of 3D points

)Qu8mg5v7";
npe_function(morton_decode)
npe_arg(codes, dense_uint, dense_ulong, dense_ulonglong)
npe_doc(morton_decode)
npe_default_arg(num_threads, int, -1)
npe_begin_code()
{
    if (codes.rows() <= 0) {
        throw pybind11::value_error("codes must be an array of shape [n] but got an empty array");
    }
    if (codes.cols() != 1) {
        throw pybind11::value_error("pts must be an array of shape [n] but got an invalid number of columns");
    }

    const int MIN_PARALLEL_INPUT_SIZE = 10000;
    const bool run_parallel = codes.rows() >= MIN_PARALLEL_INPUT_SIZE && num_threads != 0;
    auto set_parallel = OmpSetParallelism(num_threads, run_parallel);

    Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor> pts(codes.rows(), 3);

    bool threw_exception = false;
    #if defined(_OPENMP)
    #pragma omp parallel if (run_parallel)
    #endif
    {
        #if defined(_OPENMP)
        #pragma omp for
        #endif
        for(int i = 0; i < codes.rows(); i += 1) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                #else
                    threw_exception = true;
                    break;
                #endif
            }
            int32_t px, py, pz;
            MortonCode64(codes(i, 0)).decode(px, py, pz);
            pts(i, 0) = px;
            pts(i, 1) = py;
            pts(i, 2) = pz;
        }
    }
    if (threw_exception) {
        throw pybind11::error_already_set();
    }
    return npe::move(pts);
}
npe_end_code()



const char* morton_knn = R"Qu8mg5v7(
Queries a sorted array of morton encoded points to find the (approximate) k nearest neighbors

Parameters
----------
codes: an [n] shaped array of morton codes
qcodes: an [m] shaped array of query codes
k: an integer representing the number of nearest neighbors
sort_dist: (optional, defaults to True) whether to return the nearest neigbors in distance sorted order
Returns
-------
an (m, k) shaped array of indices into codes

)Qu8mg5v7";
npe_function(morton_knn)
npe_arg(codes, dense_uint, dense_ulong, dense_ulonglong)
npe_arg(qcodes, npe_matches(codes))
npe_arg(k, int)
npe_default_arg(sort_dist, bool, true)
npe_doc(morton_knn)
npe_begin_code()
{
    if (k <= 0) {
        throw pybind11::value_error("k must be greater than 0");
    }
    if (codes.rows() <= 0) {
        throw pybind11::value_error("codes must be an array of shape [n] but got an empty array");
    }
    if (qcodes.rows() <= 0) {
        throw pybind11::value_error("codes must be an array of shape [n] but got an empty array");
    }
    if (codes.cols() != 1) {
        throw pybind11::value_error("codes must be an array of shape [n] but got an invalid shape");
    }
    if (qcodes.cols() != 1) {
        throw pybind11::value_error("qcodes must be an array of shape [n] but got an invalid shape");
    }

    k = std::min(k, (int)codes.rows());

    Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic> nn_idx(qcodes.rows(), k);

    //#pragma omp parallel for
    for(int i = 0; i < qcodes.rows(); i += 1) {
        if (PyErr_CheckSignals() != 0) {
            throw pybind11::error_already_set();
        }
        npe_Scalar_codes* code_ptr = std::lower_bound(codes.data(), codes.data() + codes.rows(), qcodes(i, 0));
        std::ptrdiff_t idx = code_ptr - codes.data();

        const int half_k_up = k / 2;
        const int half_k_down = k - half_k_up;

        std::ptrdiff_t upper_bound = idx + half_k_up;
        std::ptrdiff_t lower_bound = idx - half_k_down;

        if (upper_bound >= codes.rows()) {
            lower_bound -= (upper_bound - codes.rows());
            upper_bound = codes.rows();
        }
        if (lower_bound < 0) {
            upper_bound += -lower_bound;
            lower_bound = 0;
        }

        std::vector<std::ptrdiff_t> indices;
        indices.reserve(k);
        for (int j = 0; j < (upper_bound - lower_bound); j += 1) {
            indices.push_back(lower_bound + j);
        }

        MortonCode64 code_i(*code_ptr);
        if (sort_dist) {
            auto cmp_codes = [&](const std::ptrdiff_t& lhs, const std::ptrdiff_t& rhs) {
                int32_t q_x, q_y, q_z;
                int32_t lhs_x, lhs_y, lhs_z, rhs_x, rhs_y, rhs_z;
                MortonCode64(codes(lhs, 0)).decode(lhs_x, lhs_y, lhs_z);
                MortonCode64(codes(rhs, 0)).decode(rhs_x, rhs_y, rhs_z);

                double diff_lhs_x = q_x - lhs_x;
                double diff_lhs_y = q_y - lhs_y;
                double diff_lhs_z = q_z - lhs_z;

                double diff_rhs_x = q_x - rhs_x;
                double diff_rhs_y = q_y - rhs_y;
                double diff_rhs_z = q_z - rhs_z;

                double dist_lhs = diff_lhs_x * diff_lhs_x + diff_lhs_y * diff_lhs_y + diff_lhs_z * diff_lhs_z;
                double dist_rhs = diff_rhs_x * diff_rhs_x + diff_rhs_y * diff_rhs_y + diff_rhs_z * diff_rhs_z;

                return dist_lhs < dist_rhs;
            };

            std::sort(indices.begin(), indices.end(), cmp_codes);
        }
        for (int j = 0; j < (upper_bound - lower_bound); j += 1) {
            nn_idx(i, j) = indices[j];
        }

    }

    return npe::move(nn_idx);

}
npe_end_code()



