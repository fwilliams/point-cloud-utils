#include "morton_code.h"

#include <cinttypes>
#include <cassert>

//const MortonCode64 MortonCode64::Zero(0u, 0u, 0u);
//const MortonCode64 MortonCode64::UnitX(1u, 0u, 0u);
//const MortonCode64 MortonCode64::UnitY(0u, 1u, 0u);
//const MortonCode64 MortonCode64::UnitZ(0u, 0u, 1u);


namespace {
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