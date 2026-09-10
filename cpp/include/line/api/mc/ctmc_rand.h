/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_RAND_H
#define LINE_API_MC_CTMC_RAND_H

/**
 * Random infinitesimal generator of a CTMC.
 *
 * Templated port of matlab/src/api/mc/ctmc_rand.m (identical to the kpctoolbox
 * copy) and jar/src/main/java/jline/api/mc/Ctmc_rand.java: an n x n matrix of
 * uniform [0,1) rates whose diagonal is then set by ctmc_makeinfgen.
 *
 * The MATLAB and Java versions draw from a global stream (rand / randMatrix),
 * so the caller has no way to reproduce a generator except by seeding that
 * stream. Here the source of randomness is an explicit parameter: any callable
 * returning a double in [0,1) with no global state. COMPARABILITY ACROSS
 * IMPLEMENTATIONS THEREFORE REQUIRES THE SAME GENERATOR, and a matrix built
 * here will not match one MATLAB built unless the same stream of variates is
 * fed in; what is guaranteed is that two calls with equal generator states
 * return the same matrix.
 *
 * The variates are consumed in row-major order, one per off-diagonal and
 * diagonal position alike (n*n draws), matching MATLAB's rand(n) column count
 * only in total, not in order: MATLAB fills column-major. Pass the transpose of
 * a MATLAB-filled stream if that ordering matters.
 *
 * Only the diagonal is arithmetic, so this is exact at Rational whenever the
 * generator's variates are (LcgUniform yields dyadic rationals, which are).
 */

#include <cstddef>
#include <cstdint>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * Deterministic uniform [0,1) source, so that a caller who does not have one
 * still has no reason to reach for a global stream. Numerical Recipes' 64-bit
 * linear congruential recurrence, whose high bits are the ones used.
 */
class LcgUniform {
public:
    explicit LcgUniform(std::uint64_t seed = 20260721ull) : s_(seed ? seed : 1ull) {}
    double operator()() {
        s_ = s_ * 6364136223846793005ull + 1442695040888963407ull;
        // Top 53 bits scaled into [0,1), so the value is an exact dyadic.
        return static_cast<double>(s_ >> 11) / 9007199254740992.0;
    }

private:
    std::uint64_t s_;
};

/**
 * @param n   order of the generator
 * @param gen callable returning a uniform variate in [0,1); n*n draws are made
 */
template <class T, class Gen>
Matrix<T> ctmc_rand(std::size_t n, Gen& gen) {
    if (n == 0) throw InputError("ctmc_rand: order must be positive");
    Matrix<T> R(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) R(i, j) = num_traits<T>::from_double(gen());
    return ctmc_makeinfgen(R);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_RAND_H
