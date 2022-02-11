#include <cinttypes>


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