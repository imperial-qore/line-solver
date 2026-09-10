/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `matlab/src/solvers/CTMC/ctmc_stationary.m`: the single entry point
 * for the stationary distribution of a CTMC generated from a NetworkStruct.
 *
 * All the stationary mass of a reducible chain lives in its bottom strongly
 * connected components, each weighted by the probability of being absorbed in
 * it from the declared initial state; every other state is transient and
 * carries zero. The block decomposition handles the irreducible case as the
 * degenerate one BSCC / no transient states, so every solve goes through it and
 * no dispatch can disagree with the algorithm about whether a chain is
 * reducible.
 *
 * The reference's local `ctmc_initial_distribution` turns (StateSpace, sn) into
 * the row of the initial state; here that job already belongs to
 * `analyzer_detail::init_state_index`, so this file takes the INDEX and stays
 * free of any `sn` dependency. `npos` means the initial state is absent from
 * the enumerated space -- stochastic complementation may have removed it, an
 * SPN whose immediate ENABLE states were eliminated being the usual case -- and
 * makes the block decomposition start in the SCCs with no incoming transition.
 */
#ifndef LINE_SOLVERS_CTMC_CTMC_STATIONARY_H
#define LINE_SOLVERS_CTMC_CTMC_STATIONARY_H

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_solve_reducible_blkdecomp.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/num/number.h"

namespace line {
namespace ctmc {

/** Magnitude above which an off-diagonal generator entry counts as an arc. */
static const double kArcTol = 1e-12;

template <class T>
struct CtmcStationaryResult {
    std::vector<T> pi;    ///< stationary distribution, length N
    bool seeded = false;  ///< true when a declared initial state selected the answer
    std::size_t nbscc = 0;  ///< number of closed communicating classes
    /**
     * Empty unless the chain is an UNSEEDED reducible mixture. The library
     * never writes to stderr (see `util/error.h`), so the caller decides
     * whether to surface it; the CLI prints it.
     */
    std::string warning;
};

namespace stationary_detail {

/**
 * Port of the reference's `warn_if_unseeded_mixture`.
 *
 * Without a seed the block decomposition invents a start distribution -- here a
 * uniform one over the SCCs with no incoming transition -- and no property of
 * the model implies it: on a reducible chain the stationary distribution is
 * fixed only by the initial state. Nor is it the product-form weighting, which
 * weights the recurrent classes by their unnormalized Kelly mass. The invented
 * start is order-independent and therefore looks more reproducible than the
 * answer the declared initial state selects, but that reproducibility is bought
 * by discarding the one input that makes the problem well posed.
 *
 * The fallback itself differs across codebases (python weights ALL the SCCs
 * equally, MATLAB and this port the source SCCs), which is a second reason not
 * to read the number as the model's answer.
 */
template <class T>
std::size_t count_bscc(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    if (n < 2) return 0;
    const T zero = num_traits<T>::from_int(0);
    // SIGN IS NOT A CRITERION: an ME generator embeds genuinely negative
    // off-diagonal entries, so the adjacency is taken on the MAGNITUDE.
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const double a = num_traits<T>::to_double(Q(i, j));
            if ((a < 0 ? -a : a) > kArcTol) A(i, j) = num_traits<T>::from_int(1);
        }
    const mc::SccResult s = mc::stronglyconncomp(A);
    std::size_t nbscc = 0;
    for (std::size_t c = 0; c < s.recurrent.size(); ++c)
        if (s.recurrent[c]) ++nbscc;
    return nbscc;
}

}  // namespace stationary_detail

/**
 * @param Q generator; the diagonal is recomputed by the block decomposition
 * @param init_index row of the declared initial state, or npos when absent
 */
template <class T>
CtmcStationaryResult<T> ctmc_stationary(const Matrix<T>& Q,
                                        std::size_t init_index = static_cast<std::size_t>(-1)) {
    const std::size_t npos = static_cast<std::size_t>(-1);
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_stationary: generator is not square");

    CtmcStationaryResult<T> out;
    std::vector<T> pi0;
    if (init_index != npos) {
        if (init_index >= n) throw InputError("ctmc_stationary: initial state index out of range");
        pi0.assign(n, num_traits<T>::from_int(0));
        pi0[init_index] = num_traits<T>::from_int(1);
        out.seeded = true;
    }

    out.nbscc = stationary_detail::count_bscc(Q);
    if (!out.seeded && out.nbscc > 1) {
        out.warning =
            "SolverCTMC: the generator has " + std::to_string(out.nbscc) +
            " closed communicating classes and the declared initial state could not be located in "
            "the enumerated state space, so the solve starts from a distribution the model never "
            "stated (uniform over the SCCs with no incoming transition). On a reducible chain the "
            "stationary distribution is determined only by the initial state, so this answer is "
            "not the model's. Call setState on the stations so the class the model actually starts "
            "in is the one solved.";
    }

    out.pi = mc::ctmc_solve_reducible_blkdecomp(Q, pi0).pi;
    return out;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_CTMC_STATIONARY_H
