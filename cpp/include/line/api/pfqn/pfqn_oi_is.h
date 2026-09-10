/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_OI_IS_H
#define LINE_API_PFQN_OI_IS_H

/**
 * Importance-sampling estimate of the normalizing constant of a closed
 * two-station order-independent (OI) tandem.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_oi_is.m. That file is
 * pfqn_pas_is.m with an EMPTY swap graph: the sampler places any present class
 * (no placement constraint), so the communicating class is the full set of
 * orderings and the estimand is the plain OI constant of pfqn_ncoi rather
 * than a per-communicating-class one. The two MATLAB files carry the same
 * estimator body verbatim, so this header delegates instead of duplicating it;
 * pas_placement on an empty H returns the empty closure, which makes the
 * placement test in pfqn_pas_is unconditionally true and reproduces
 * pfqn_oi_is.m line for line.
 *
 * There is no Java counterpart of pfqn_oi_is in jar/src/main/java/jline/api/
 * pfqn/nc/, so MATLAB is the only reference here.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION, gated on has_transcendental for the
 * same reason as pfqn_pas_is (a random output, plus lG = log(G)).
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/api/pfqn/pfqn_pas_is.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param N       (R) closed population vector
 * @param mu      the two OI rank-rate functions, station 1 then station 2
 * @param samples number of importance samples
 * @param rng     explicit generator, advanced by the call
 */
template <class T>
PasIsResult<T> pfqn_oi_is(const std::vector<int>& N, const std::vector<OiRateFun<T>>& mu,
                          std::size_t samples, McRng& rng, bool want_qlen = true) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_oi_is requires transcendental arithmetic: it is a Monte Carlo estimator, "
                  "inexact by construction, and reports the log of its own estimate");
    if (mu.size() != 2)
        throw InputError(
            "pfqn_oi_is models a two-station OI tandem: mu must have exactly two rate functions");
    return pfqn_pas_is(N, mu, Matrix<int>(), samples, rng, want_qlen);
}

/** Reference default of 1e4 samples. */
template <class T>
PasIsResult<T> pfqn_oi_is(const std::vector<int>& N, const std::vector<OiRateFun<T>>& mu,
                          McRng& rng, bool want_qlen = true) {
    return pfqn_oi_is(N, mu, static_cast<std::size_t>(10000), rng, want_qlen);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_OI_IS_H
