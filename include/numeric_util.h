// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

namespace slib {

template <typename T>
inline
T clamp(T value, T low, T high) {
    return (value < low) ? low : (value > high ? high : value);
}

/// check if an integer is a power of two.
/// @see http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
inline
bool IsPowerOfTwo(int n) {
    static_assert(sizeof(int) == 4, "non 32-bit integer");
    return n && !(n & (n - 1));
}

/// compute the smallest power of two larger than or equal to an integer.
/// @see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
inline
unsigned int FineNextPowerOfTwo(unsigned int n) {
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    if (sizeof(unsigned int) > 1) {
        n |= n >> 8;    // >= 16bit
    }
    if (sizeof(size_t) > 2) {
        n |= n >> 16;    // >= 32bit
    }
#if 0
    if (sizeof(size_t) > 4) {
        n |= n >> 32;    // >= 64bit
    }
#endif
    return ++n;
}

}