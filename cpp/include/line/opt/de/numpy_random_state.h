/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_DE_NUMPY_RANDOM_STATE_H
#define LINE_OPT_DE_NUMPY_RANDOM_STATE_H

/**
 * Bit-exact port of `matlab/src/opt/+opt/+de/NumpyRandomState.m`.
 *
 * Only the legacy NumPy RandomState operations consumed by line-opt are
 * exposed: random_sample, uniform, default-dtype randint, shuffle and
 * permutation.  Their word consumption is part of the contract.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "line/opt/de/mt19937.h"
#include "line/util/error.h"

namespace line {
namespace opt {
namespace de {

class NumpyRandomState {
public:
    explicit NumpyRandomState(std::uint64_t seed) : gen_(seed) {}

    MT19937& generator() { return gen_; }
    const MT19937& generator() const { return gen_; }

    /** NumPy `random_sample()`: one 53-bit double from two generator words. */
    double random_sample() {
        const std::uint64_t a = static_cast<std::uint64_t>(gen_.next_uint32() >> 5);
        const std::uint64_t b = static_cast<std::uint64_t>(gen_.next_uint32() >> 6);
        return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) /
               9007199254740992.0;
    }

    std::vector<double> random_sample(std::size_t n) {
        std::vector<double> out(n);
        for (double& value : out) value = random_sample();
        return out;
    }

    double uniform(double low = 0.0, double high = 1.0) {
        return low + (high - low) * random_sample();
    }

    std::vector<double> uniform(double low, double high, std::size_t n) {
        std::vector<double> out(n);
        for (double& value : out) value = low + (high - low) * random_sample();
        return out;
    }

    std::uint64_t randint(std::uint64_t low, std::uint64_t high) {
        if (high <= low) throw InputError("NumpyRandomState::randint: high must exceed low");
        const std::uint64_t range = high - 1 - low;
        if (range > UINT64_C(0xffffffff))
            throw InputError(
                "NumpyRandomState::randint: the line-opt RandomState subset uses 32-bit ranges");
        if (range == 0) return low;
        const std::uint64_t mask = fill_mask(range);
        for (;;) {
            const std::uint64_t value = static_cast<std::uint64_t>(gen_.next_uint32()) & mask;
            if (value <= range) return low + value;
        }
    }

    std::uint64_t randint(std::uint64_t high) { return randint(0, high); }

    std::uint64_t random_interval(std::uint64_t max_value) {
        if (max_value > UINT64_C(0xffffffff))
            throw InputError("NumpyRandomState::random_interval: max exceeds 32 bits");
        if (max_value == 0) return 0;
        const std::uint64_t mask = fill_mask(max_value);
        for (;;) {
            const std::uint64_t value = static_cast<std::uint64_t>(gen_.next_uint32()) & mask;
            if (value <= max_value) return value;
        }
    }

    template <class T>
    void shuffle(std::vector<T>& values) {
        for (std::size_t p = values.size(); p > 1; --p) {
            const std::size_t j = static_cast<std::size_t>(random_interval(p - 1));
            std::swap(values[p - 1], values[j]);
        }
    }

    std::vector<std::size_t> permutation(std::size_t n) {
        std::vector<std::size_t> out(n);
        std::iota(out.begin(), out.end(), std::size_t(0));
        shuffle(out);
        return out;
    }

private:
    static std::uint64_t fill_mask(std::uint64_t value) {
        value |= value >> 1;
        value |= value >> 2;
        value |= value >> 4;
        value |= value >> 8;
        value |= value >> 16;
        return value;
    }

    MT19937 gen_;
};

}  // namespace de
}  // namespace opt
}  // namespace line

#endif
