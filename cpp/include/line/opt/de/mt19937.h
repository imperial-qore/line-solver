/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_DE_MT19937_H
#define LINE_OPT_DE_MT19937_H

/**
 * Bit-exact port of `matlab/src/opt/+opt/+de/MT19937.m`.
 *
 * This is NumPy's legacy MT19937 core, including its scalar seeding rule.  It
 * intentionally does not use `std::mt19937`: the standard fixes the generator
 * recurrence but not NumPy RandomState's seed expansion and draw adapters, and
 * line-opt requires the same word stream in every language so Differential
 * Evolution follows the same trajectory for a fixed seed.
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace opt {
namespace de {

class MT19937 {
public:
    static constexpr std::size_t state_size = 624;

    explicit MT19937(std::uint64_t seed_value) { seed(seed_value); }

    /** NumPy legacy scalar seeding (`init_genrand`, or little-endian words). */
    void seed(std::uint64_t value) {
        if (value <= UINT64_C(0xffffffff)) {
            init_genrand(static_cast<std::uint32_t>(value));
            return;
        }
        std::vector<std::uint32_t> key;
        while (value != 0) {
            key.push_back(static_cast<std::uint32_t>(value & UINT64_C(0xffffffff)));
            value >>= 32;
        }
        init_by_array(key);
    }

    void init_genrand(std::uint32_t value) {
        mt_[0] = value;
        for (std::size_t i = 1; i < state_size; ++i) {
            const std::uint32_t prev = mt_[i - 1];
            mt_[i] = UINT32_C(1812433253) * (prev ^ (prev >> 30)) +
                     static_cast<std::uint32_t>(i);
        }
        pos_ = state_size;
    }

    void init_by_array(const std::vector<std::uint32_t>& key) {
        if (key.empty()) throw InputError("MT19937::init_by_array: the seed key is empty");
        init_genrand(UINT32_C(19650218));
        std::size_t i = 1, j = 0;
        std::size_t k = state_size > key.size() ? state_size : key.size();
        for (; k != 0; --k) {
            const std::uint32_t prev = mt_[i - 1];
            mt_[i] = (mt_[i] ^ ((prev ^ (prev >> 30)) * UINT32_C(1664525))) + key[j] +
                     static_cast<std::uint32_t>(j);
            ++i;
            ++j;
            if (i >= state_size) {
                mt_[0] = mt_[state_size - 1];
                i = 1;
            }
            if (j >= key.size()) j = 0;
        }
        for (k = state_size - 1; k != 0; --k) {
            const std::uint32_t prev = mt_[i - 1];
            mt_[i] = (mt_[i] ^ ((prev ^ (prev >> 30)) * UINT32_C(1566083941))) -
                     static_cast<std::uint32_t>(i);
            ++i;
            if (i >= state_size) {
                mt_[0] = mt_[state_size - 1];
                i = 1;
            }
        }
        mt_[0] = UINT32_C(0x80000000);
        pos_ = state_size;
    }

    /** Next tempered word, exactly `MT19937.nextUint32` in the MATLAB port. */
    std::uint32_t next_uint32() {
        constexpr std::size_t M = 397;
        constexpr std::uint32_t matrix_a = UINT32_C(0x9908b0df);
        constexpr std::uint32_t upper = UINT32_C(0x80000000);
        constexpr std::uint32_t lower = UINT32_C(0x7fffffff);
        if (pos_ >= state_size) {
            std::size_t k = 0;
            for (; k < state_size - M; ++k) {
                const std::uint32_t y = (mt_[k] & upper) | (mt_[k + 1] & lower);
                mt_[k] = mt_[k + M] ^ (y >> 1) ^ ((y & 1U) ? matrix_a : 0U);
            }
            for (; k < state_size - 1; ++k) {
                const std::uint32_t y = (mt_[k] & upper) | (mt_[k + 1] & lower);
                mt_[k] = mt_[k + (M - state_size)] ^ (y >> 1) ^
                         ((y & 1U) ? matrix_a : 0U);
            }
            const std::uint32_t y = (mt_[state_size - 1] & upper) | (mt_[0] & lower);
            mt_[state_size - 1] = mt_[M - 1] ^ (y >> 1) ^ ((y & 1U) ? matrix_a : 0U);
            pos_ = 0;
        }

        std::uint32_t y = mt_[pos_++];
        y ^= y >> 11;
        y ^= (y << 7) & UINT32_C(0x9d2c5680);
        y ^= (y << 15) & UINT32_C(0xefc60000);
        y ^= y >> 18;
        return y;
    }

    const std::array<std::uint32_t, state_size>& state_key() const { return mt_; }
    std::size_t position() const { return pos_; }

private:
    std::array<std::uint32_t, state_size> mt_{};
    std::size_t pos_ = state_size;
};

}  // namespace de
}  // namespace opt
}  // namespace line

#endif
