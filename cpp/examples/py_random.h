/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_EXAMPLES_PY_RANDOM_H
#define LINE_EXAMPLES_PY_RANDOM_H

/**
 * `random.seed(n)` / `random.random()`, bit for bit.
 *
 * Several gallery and example models DRAW their parameters -- a population, a
 * service mean -- from Python's `random` after a fixed seed, so the model that
 * `gallery_cqn(M)` builds is a function of that stream. Writing the realized
 * numbers here as literals would fix only the default M: the reference draws
 * M + 1 of them, so any other M silently gets the wrong model.
 *
 * IT IS THE SAME GENERATOR. CPython's `random` module is MT19937, which is
 * `std::mt19937`; what is NOT shared is the seeding and the double conversion,
 * and those are the whole of this header:
 *
 *   - `random.seed(n)` for a non-negative int calls `init_by_array` over the
 *     32-bit little-endian words of |n| (`_randommodule.c:random_seed`), NOT
 *     `init_genrand(n)` and not `std::seed_seq`.
 *   - `random.random()` is `genrand_res53`: two draws, 27 and 26 bits, combined
 *     as (a * 2^26 + b) / 2^53.
 *
 * The core is written out rather than delegated to `std::mt19937` because the
 * tempering and the twist have to interleave with `init_by_array`'s state
 * exactly; a `std::mt19937` seeded by any public route lands on a different
 * state vector, which is the whole failure mode this header exists to avoid.
 *
 * NOT numpy. `numpy.random.seed` scalar-seeds with `init_genrand` instead, so a
 * model drawn from numpy is a different stream and this class is not it.
 *
 * MATLAB DRAWS A THIRD STREAM. `rand` after `rng(seed)` is MT19937 with yet
 * another seeding, so `gallery_cqn` in MATLAB and in Python are ALREADY
 * different models; this header reproduces the Python one, which is the
 * reference `cpp/examples` is transcribed from, and that choice is stated here
 * rather than left to be inferred from the numbers.
 */

#include <cstddef>
#include <cstdint>
#include <vector>

namespace line {
namespace examples {

/** CPython's `random.Random`, restricted to `seed(int)` and `random()`. */
class PyRandom {
  public:
    explicit PyRandom(unsigned long long seed) { this->seed(seed); }

    /** `random.seed(n)`: init_by_array over the 32-bit words of n, low word first. */
    void seed(unsigned long long n) {
        std::vector<std::uint32_t> key;
        if (n == 0) {
            key.push_back(0);
        } else {
            for (unsigned long long v = n; v != 0; v >>= 32)
                key.push_back(static_cast<std::uint32_t>(v & 0xffffffffULL));
        }
        init_by_array(key);
    }

    /** `random.random()`: genrand_res53, the 53-bit double of two draws. */
    double random() {
        const std::uint32_t a = genrand_uint32() >> 5, b = genrand_uint32() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

  private:
    static const std::size_t N = 624, M = 397;
    static const std::uint32_t MATRIX_A = 0x9908b0dfUL;
    static const std::uint32_t UPPER_MASK = 0x80000000UL;
    static const std::uint32_t LOWER_MASK = 0x7fffffffUL;

    std::uint32_t mt_[N];
    std::size_t mti_ = N + 1;

    void init_genrand(std::uint32_t s) {
        mt_[0] = s;
        for (mti_ = 1; mti_ < N; ++mti_)
            mt_[mti_] = static_cast<std::uint32_t>(
                1812433253UL * (mt_[mti_ - 1] ^ (mt_[mti_ - 1] >> 30)) +
                static_cast<std::uint32_t>(mti_));
    }

    void init_by_array(const std::vector<std::uint32_t>& key) {
        init_genrand(19650218UL);
        std::size_t i = 1, j = 0;
        std::size_t k = N > key.size() ? N : key.size();
        for (; k; --k) {
            mt_[i] = static_cast<std::uint32_t>(
                         (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1664525UL)) + key[j]) +
                     static_cast<std::uint32_t>(j);
            ++i;
            ++j;
            if (i >= N) {
                mt_[0] = mt_[N - 1];
                i = 1;
            }
            if (j >= key.size()) j = 0;
        }
        for (k = N - 1; k; --k) {
            mt_[i] = static_cast<std::uint32_t>(
                         (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1566083941UL)) -
                         static_cast<std::uint32_t>(i));
            ++i;
            if (i >= N) {
                mt_[0] = mt_[N - 1];
                i = 1;
            }
        }
        mt_[0] = 0x80000000UL;
    }

    std::uint32_t genrand_uint32() {
        std::uint32_t y = 0;
        if (mti_ >= N) {
            static const std::uint32_t mag01[2] = {0x0UL, MATRIX_A};
            std::size_t kk = 0;
            if (mti_ == N + 1) init_genrand(5489UL);
            for (; kk < N - M; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk + 1] & LOWER_MASK);
                mt_[kk] = mt_[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            for (; kk < N - 1; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk + 1] & LOWER_MASK);
                mt_[kk] = mt_[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            y = (mt_[N - 1] & UPPER_MASK) | (mt_[0] & LOWER_MASK);
            mt_[N - 1] = mt_[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mti_ = 0;
        }
        y = mt_[mti_++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);
        return y;
    }
};

/**
 * Python's `round`, which is BANKER'S ROUNDING: a value exactly halfway goes to
 * the even neighbour, where C's `round` goes away from zero. The populations of
 * `gallery_cqn` and `gallery_repairmen` are `round(random() * ...)`, so the two
 * rules differ on a half-integer draw and the model would differ with them.
 */
inline double py_round(double v) {
    const double f = std::floor(v), diff = v - f;
    if (diff > 0.5) return f + 1.0;
    if (diff < 0.5) return f;
    return (std::fmod(f, 2.0) == 0.0) ? f : f + 1.0;
}

}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_PY_RANDOM_H
