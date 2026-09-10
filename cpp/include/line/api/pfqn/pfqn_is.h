/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_IS_H
#define LINE_API_PFQN_IS_H

/**
 * Importance-sampling estimate of the normalizing constant of a closed
 * LOAD-INDEPENDENT product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_is.m, which is one line: the
 * load-independent case mu_i(k) = 1 of pfqn_ld_is. This header keeps the same
 * shape so that the two entry points stay in step, rather than duplicating the
 * estimator.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION, gated on has_transcendental for the
 * same reason as pfqn_ld_is (a random output, plus lG = log(G)).
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_ld_is.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L       (M x R) per-class demands at the M single-server queues
 * @param N       (R) closed population vector
 * @param Z       (R) aggregated think times; empty or all zero for no delay
 * @param samples number of importance samples
 * @param rng     explicit generator, advanced by the call
 */
template <class T>
NcResult<T> pfqn_is(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                    std::size_t samples, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_is requires transcendental arithmetic: it is a Monte Carlo estimator, "
                  "inexact by construction, and reports the log of its own estimate");
    return pfqn_ld_is(L, N, Z, Matrix<T>(), samples, rng);
}

/** Reference default of 1e4 samples. */
template <class T>
NcResult<T> pfqn_is(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                    McRng& rng) {
    return pfqn_is(L, N, Z, static_cast<std::size_t>(10000), rng);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_IS_H
