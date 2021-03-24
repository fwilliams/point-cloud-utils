#include <cinttypes>
#include <cassert>
#include <algorithm>
#include <iostream>

#include <npe.h>

#include "common.h"

namespace {
//Represents a three-dimensional 64-bit Morton Code.
class MortonCode64
{
public:
    MortonCode64();

    MortonCode64(int32_t x, int32_t y, int32_t z);

    MortonCode64(uint32_t x, uint32_t y, uint32_t z);

    MortonCode64(uint64_t);

    //Decodes the code into its three coordinates.
    void decode(int32_t& x, int32_t& y, int32_t& z) const;

    //Negates the specified coordinate.
    template <int DIM>
    MortonCode64 InvertDimension() const;

    //Under the assumption that all entries are positive
    MortonCode64 DivideDimensionBy2(int dim) const;
    MortonCode64 Negate() const;

    MortonCode64 operator+(const MortonCode64 rhs) const;
    MortonCode64 operator+(const int64_t rhs) const;
    MortonCode64& operator+=(const MortonCode64 rhs);
    MortonCode64 operator-(const MortonCode64 rhs) const;

    //Applies the bitshift to every dimension, assumes all entries to be positive
    MortonCode64 operator >> (int shift) const;
    MortonCode64 operator << (int shift) const;

    bool operator<(const MortonCode64 rhs) const { return data < rhs.data; }
    bool operator>(const MortonCode64 rhs) const { return data > rhs.data; }
    bool operator<=(const MortonCode64 rhs) const { return data <= rhs.data; }
    bool operator>=(const MortonCode64 rhs) const { return data >= rhs.data; }
    bool operator==(const MortonCode64 rhs) const { return data == rhs.data; }
    bool operator!=(const MortonCode64 rhs) const { return data != rhs.data; }

    explicit operator uint64_t() const { return data; }

    uint64_t get_data() const { return data; }

    static const MortonCode64 Zero;
    static const MortonCode64 UnitX;
    static const MortonCode64 UnitY;
    static const MortonCode64 UnitZ;

private:
    uint64_t data;
};


const MortonCode64 MortonCode64::Zero(0u, 0u, 0u);
const MortonCode64 MortonCode64::UnitX(1u, 0u, 0u);
const MortonCode64 MortonCode64::UnitY(0u, 1u, 0u);
const MortonCode64 MortonCode64::UnitZ(0u, 0u, 1u);


uint64_t SplitBy3Bits21(int32_t x)
{
    //int64_t r = r & 0x1fffff;
    uint64_t r = x;

    r = (r | r << 32) & 0x1f00000000ffff;
    r = (r | r << 16) & 0x1f0000ff0000ff;
    r = (r | r << 8) & 0x100f00f00f00f00f;
    r = (r | r << 4) & 0x10c30c30c30c30c3;
    r = (r | r << 2) & 0x1249249249249249;

    return r;
}

int32_t CompactBy3Bits21(uint64_t x)
{
    uint64_t d = x & 0x1249249249249249;
    d = (d | d >> 2) & 0x10c30c30c30c30c3;
    d = (d | d >> 4) & 0x100f00f00f00f00f;
    d = (d | d >> 8) & 0x1f0000ff0000ff;
    d = (d | d >> 16) & 0x1f00000000ffff;
    d = (d | d >> 32);

    //sign extension
    d = (d & 0x100000 ? d | 0xffe00000 : d);

    return (int32_t)d;
}

MortonCode64::MortonCode64()
{ }

MortonCode64::MortonCode64(int32_t x, int32_t y, int32_t z)
{
    assert(-(1 << 20) <= x && x < (1 << 20));
    assert(-(1 << 20) <= y && y < (1 << 20));
    assert(-(1 << 20) <= z && z < (1 << 20));
    //move sign bit to bit 20
    x = (x & 0x80000000) >> 11 | (x & 0x0fffff);
    y = (y & 0x80000000) >> 11 | (y & 0x0fffff);
    z = (z & 0x80000000) >> 11 | (z & 0x0fffff);

    data = SplitBy3Bits21(x) | SplitBy3Bits21(y) << 1 | SplitBy3Bits21(z) << 2;
    //invert sign bits
    data = data ^ 0x7000000000000000;

    /*int32_t dx, dy, dz;
    decode(dx, dy, dz);
    assert(x == dx);
    assert(y == dy);
    assert(z == dz);*/
}

MortonCode64::MortonCode64(uint32_t x, uint32_t y, uint32_t z)
{
    data = (SplitBy3Bits21(x) | SplitBy3Bits21(y) << 1 | SplitBy3Bits21(z) << 2) | 0x7000000000000000;
}

MortonCode64::MortonCode64(uint64_t d)
    : data(d)
{ }

void MortonCode64::decode(int32_t & x, int32_t & y, int32_t & z) const
{
    //invert sign bits
    uint64_t d = data ^ 0x7000000000000000;

    x = CompactBy3Bits21(d);
    y = CompactBy3Bits21(d >> 1);
    z = CompactBy3Bits21(d >> 2);
}

int64_t xMask = 0x1249249249249249; // ...001001001001001001001

template <int DIM>
MortonCode64 MortonCode64::InvertDimension() const
{
    static_assert(0 <= DIM && DIM < 3, "The dimension must be between 0 and 2.");
    uint64_t mask = xMask << DIM;

    //(2 - complement)
    //invert bits
    uint64_t c = data ^ mask;
    //add 1
    uint64_t sum = (c | ~mask) + 1; //fill the bits that are not in the mask with 1 to make sure carry gets carried on
    return MortonCode64((sum & mask) | (c & ~mask));
}

template MortonCode64 MortonCode64::InvertDimension<0>() const;
template MortonCode64 MortonCode64::InvertDimension<1>() const;
template MortonCode64 MortonCode64::InvertDimension<2>() const;

MortonCode64 MortonCode64::DivideDimensionBy2(int dim) const
{
    assert(data >> 60 == 7); //all signs must be positive
    assert(0 <= dim && dim < 3);
    uint64_t mask = (xMask << dim) & (0x0fffffffffffffff); //exclude signs from mask

    return MortonCode64((data & mask) >> 3 | (data & ~mask));
}

MortonCode64 MortonCode64::Negate() const
{
    uint64_t yMask = xMask << 1;
    uint64_t zMask = xMask << 2;

    //invert
    uint64_t d = ~data;
    uint64_t xSum = (d | ~xMask) + 1;
    uint64_t ySum = (d | ~yMask) + 1;
    uint64_t zSum = (d | ~zMask) + 1;

    return MortonCode64((xSum & xMask) | (ySum & yMask) | (zSum & zMask));
}

MortonCode64 MortonCode64::operator+(const MortonCode64 rhs) const
{
    //invert sign bits
    uint64_t code1 = data ^ 0x7000000000000000;
    uint64_t code2 = rhs.data ^ 0x7000000000000000;

    uint64_t yMask = xMask << 1;
    uint64_t zMask = xMask << 2;

    uint64_t xSum = (code1 | ~xMask) + (code2 & xMask);
    uint64_t ySum = (code1 | ~yMask) + (code2 & yMask);
    uint64_t zSum = (code1 | ~zMask) + (code2 & zMask);

    uint64_t result = (xSum & xMask) | (ySum & yMask) | (zSum & zMask);
    return MortonCode64(result ^ 0x7000000000000000);
}

MortonCode64 MortonCode64::operator+(const int64_t rhs) const
{
    return MortonCode64(data + rhs);
}

MortonCode64 & MortonCode64::operator+=(const MortonCode64 rhs)
{
    *this = *this + rhs;
    return *this;
}

MortonCode64 MortonCode64::operator-(const MortonCode64 rhs) const
{
    return *this + rhs.Negate();
}

MortonCode64 MortonCode64::operator>>(int shift) const
{
    assert(data >> 60 == 7); //all signs must be positive
    return MortonCode64((data & 0x0fffffffffffffff) >> (3 * shift) | 0x7000000000000000);
}

MortonCode64 MortonCode64::operator<<(int shift) const
{
    return MortonCode64(((data << (3 * shift)) & 0x0fffffffffffffff) | (data & 0x7000000000000000));
}

}




const char* morton_encode = R"Qu8mg5v7(
Encode n 3D points using Morton coding, possibly sorting them

Parameters
----------
pts: an [n, 3] array of 3D points
sort: (optional, default to false) sort the points

Returns
-------
an [n] shaped array of morton encoded points

)Qu8mg5v7";
npe_function(morton_encode)
npe_arg(pts, dense_int, dense_long, dense_longlong)
npe_default_arg(sort, bool, false)
npe_default_arg(parallel, bool, true)
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

    if (parallel) {
        #pragma omp parallel for
        for(int i = 0; i < pts.rows(); i += 1) {
            int32_t px = pts(i, 0), py = pts(i, 1), pz = pts(i, 2);
            MortonCode64 code(px, py, pz);
            codes[i] = code.get_data();
        }
    } else {
        for(int i = 0; i < pts.rows(); i += 1) {
            int32_t px = pts(i, 0), py = pts(i, 1), pz = pts(i, 2);
            MortonCode64 code(px, py, pz);
            codes[i] = code.get_data();
        }
    }

    if (sort) {
        std::sort(codes.data(), codes.data() + codes.rows());
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
npe_begin_code()
{
    if (codes.rows() <= 0) {
        throw pybind11::value_error("codes must be an array of shape [n] but got an empty array");
    }
    if (codes.cols() != 1) {
        throw pybind11::value_error("pts must be an array of shape [n] but got an invalid number of columns");
    }

    Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor> pts(codes.rows(), 3);

#pragma omp parallel for
    for(int i = 0; i < codes.rows(); i += 1) {
        int32_t px, py, pz;
        MortonCode64(codes(i, 0)).decode(px, py, pz);
        pts(i, 0) = px;
        pts(i, 1) = py;
        pts(i, 2) = pz;
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



