/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_MULTI_H
#define LINE_API_MC_CTMC_MULTI_H

/**
 * Two-level multigrid aggregation-disaggregation for a nearly completely
 * decomposable CTMC.
 *
 * Templated port of matlab/src/api/mc/ctmc_multi.m and
 * jar/src/main/java/jline/api/mc/Ctmc_multi.java. The construction is exactly
 * Courtois's: permute into macro-state order, decouple, solve each diagonal
 * block for its conditional distribution, and assemble the macro-state chain G.
 * The one difference, and the whole point of the method, is that G is not
 * solved directly but decomposed AGAIN, by a second Courtois step over the
 * macro-macro-states MSS, so the coarse problem is itself solved by aggregation.
 * This is the one-step, two-level instance of multigrid; a full multi-level
 * implementation is repeated coarsening with the same base method.
 *
 * The fine level is shared verbatim with ctmc_courtois rather than duplicated,
 * so the two cannot drift apart.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC: both levels are Courtois steps, and
 * epsMAX is an eigenvalue modulus.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_courtois.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct MultiResult {
    std::vector<T> p;       ///< approximate stationary vector, ORIGINAL ordering
    std::vector<T> pcourt;  ///< the plain Courtois estimate, for comparison
    Matrix<T> Qperm;        ///< Q reordered by macro-state
    T eps;                  ///< NCD index of the fine level
    T epsMAX;               ///< maximum admissible NCD index of the fine level
};

/**
 * @param Q   generator
 * @param MS  macro-states partitioning 0..n-1
 * @param MSS macro-macro-states partitioning 0..|MS|-1, the coarse partition
 * @param q   uniformization rate for the fine level
 */
template <class T>
MultiResult<T> ctmc_multi(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                          const std::vector<std::vector<std::size_t>>& MSS, const T& q) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_multi requires transcendental arithmetic: both levels are Courtois "
                  "decompositions, whose epsMAX is an iteratively computed eigenvalue modulus");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_multi: generator is not square");
    if (MSS.empty()) throw InputError("ctmc_multi: no macro-macro-states given");

    const detail::CourtoisCore<T> c = detail::courtois_core(Q, MS, q);
    // The coarse solve: a second Courtois decomposition of the macro chain, in
    // place of the direct dtmc_solve that ctmc_courtois performs.
    const std::vector<T> pMacro = ctmc_courtois(c.G, MSS).p;

    std::vector<T> pperm(n, num_traits<T>::from_int(0));
    std::size_t proc = 0;
    for (std::size_t i = 0; i < MS.size(); ++i) {
        for (std::size_t a = 0; a < MS[i].size(); ++a) pperm[proc + a] = pMacro[i] * c.pmicro[proc + a];
        proc += MS[i].size();
    }

    MultiResult<T> r;
    r.p = detail::unpermute_states(pperm, c.v);
    r.Qperm = c.Qperm;
    r.eps = c.eps;
    r.epsMAX = c.epsMAX;
    r.pcourt = ctmc_courtois(Q, MS, q).p;
    return r;
}

/** Overload deriving the rate as MATLAB does, q = (21/20) max|Qperm|. */
template <class T>
MultiResult<T> ctmc_multi(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                          const std::vector<std::vector<std::size_t>>& MSS) {
    return ctmc_multi(Q, MS, MSS, detail::courtois_default_rate(Q, MS));
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_MULTI_H
