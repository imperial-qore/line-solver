/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SIMULATE_H
#define LINE_API_MC_CTMC_SIMULATE_H

/**
 * Sample path of a continuous-time Markov chain given its generator.
 *
 * Templated port of matlab/src/api/mc/ctmc_simulate.m. The chain is simulated
 * by the standard jump-chain construction: from state i the holding time is
 * exponential with mean -1/Q(i,i), and the next state is drawn from the
 * embedded jump chain P(i,j) = Q(i,j) / sum_{k != i} Q(i,k).
 *
 * This is the one entry point in this batch that takes MATRICES and not an
 * `sn`: its arguments are the generator, an initial distribution and a step
 * count, so it is portable independently of the NetworkStruct layer.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental. The holding times
 * are exponential deviates, drawn by inverse transform as -mean * log(u), so
 * there is no exact instantiation: the sample path is a realization, not a
 * number that a rational field could represent.
 *
 * RANDOMNESS. The generator is `line::pfqn::McRng` (a std::mt19937_64) passed
 * by reference and advanced by the call, the same convention every Monte Carlo
 * entry point in this tree uses, so a caller controls reproducibility by
 * controlling the generator state. The stream is NOT comparable with MATLAB's
 * -- different generator, different mapping from bits to deviates -- so the
 * oracle for this function is distributional, never path-for-path. Two
 * generator steps are consumed per simulated step: one for the holding time
 * and one for the jump.
 *
 * REFERENCE DEFECT: the initial state is not drawn from pi0.
 *
 * ctmc_simulate.m selects the starting state with
 *
 *     [~, st] = min(abs(rand - cumsum(pi0)));
 *
 * which returns the state whose CUMULATIVE probability is nearest to the
 * uniform deviate, ties going to the lowest index because min returns the
 * first minimizer. That partitions [0,1] at the MIDPOINTS between consecutive
 * DISTINCT cumulative values c_k = sum_{j<=k} p_j, instead of at the
 * cumulative values themselves. Where the c_k are distinct this reduces to
 *
 *     P(1) = p_1 + p_2/2,   P(k) = (p_k + p_{k+1})/2,   P(n) = p_n/2,
 *
 * i.e. the LAST state always receives about half its intended mass and the
 * first receives an excess. Where two consecutive c_k coincide -- which is
 * exactly what a zero-probability entry produces -- the first of them takes
 * the whole window and the rest get nothing, so the naive midpoint reading
 * above does NOT apply there. Measured in MATLAB over 400000 draws:
 *
 *   pi0 = [0.5 0.5]        -> [0.7498 0.2502]   (want [0.5 0.5])
 *   pi0 = [1/3 1/3 1/3]    -> [0.4997 0.3335 0.1668]  (want [1/3 1/3 1/3])
 *   pi0 = [0.9 0 0.1]      -> [0.9497 0.0000 0.0503]  (want [0.9 0 0.1])
 *
 * the last of which is also reproducible through the entry point itself,
 * ctmc_simulate(Q, [0.9;0;0.1], 1) over 20000 calls giving
 * [0.9488 0.0000 0.0512]. The error washes out of a long ergodic run -- the
 * time-average occupancy of a two-state chain still converges to the exact
 * stationary law -- so it is invisible in steady-state use and corrupts
 * exactly the transient and short-run use that passing pi0 is for.
 *
 * The port draws the initial state by correct inverse transform. Reproducing
 * the defect was rejected: a sampler that does not sample from the
 * distribution it is handed has no contract left to preserve, and unlike a
 * closed-form value there is nothing downstream that could be calibrated
 * against the wrong answer. `ctmc_simulate_reference_initial_law` below
 * returns the law the reference actually realizes, so the discrepancy is
 * available to a caller as a value rather than only as prose.
 *
 * A second, smaller divergence: an ABSORBING state (a row whose off-diagonal
 * entries are all zero) gives MATLAB `F = 0/0 = NaN` on that row, after which
 * `find(rand - NaN > 0)` is empty and the chain silently jumps to state 1,
 * while the holding time is exprnd(Inf) = Inf. The port raises NumericError
 * naming the state instead: there is no Inf in an exact field, and a silent
 * jump to state 1 is not a property of the chain.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** One simulated sample path: the state visited at each step and its holding time. */
template <class T>
struct CtmcPath {
    std::vector<std::size_t> states;    ///< 0-based state index at each step
    std::vector<T> sojourn;             ///< holding time spent in that state
};

/**
 * The initial-state law that ctmc_simulate.m actually realizes for a given
 * pi0, as opposed to pi0 itself. Provided so a caller can measure the
 * reference defect documented above rather than take it on trust; it is not
 * used by the simulation.
 */
template <class T>
std::vector<T> ctmc_simulate_reference_initial_law(const std::vector<T>& pi0) {
    const std::size_t n = pi0.size();
    if (n == 0) throw InputError("ctmc_simulate_reference_initial_law: empty distribution");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T half = num_traits<T>::from_rational(1, 2);

    T tot = zero;
    for (std::size_t i = 0; i < n; ++i) {
        if (pi0[i] < zero)
            throw InputError("ctmc_simulate_reference_initial_law: negative entry");
        tot += pi0[i];
    }
    if (tot == zero) throw InputError("ctmc_simulate_reference_initial_law: zero total mass");

    std::vector<T> c(n);
    {
        T run = zero;
        for (std::size_t i = 0; i < n; ++i) {
            run += pi0[i] / tot;
            c[i] = run;
        }
    }

    // tie-breaking rationale: see _kb/03-api-layer.md (cpp port notes: mc)
    std::vector<std::size_t> rep;
    for (std::size_t k = 0; k < n; ++k)
        if (rep.empty() || c[k] != c[rep.back()]) rep.push_back(k);

    std::vector<T> q(n, zero);
    const std::size_t M = rep.size();
    for (std::size_t m = 0; m < M; ++m) {
        const T lo = (m == 0) ? zero : T(half * (c[rep[m - 1]] + c[rep[m]]));
        const T hi = (m + 1 == M) ? one : T(half * (c[rep[m]] + c[rep[m + 1]]));
        q[rep[m]] = hi - lo;
    }
    return q;
}

/**
 * Simulate n steps of the CTMC with generator Q.
 *
 * @param Q    (m x m) generator, negative diagonal and zero row sums
 * @param pi0  initial distribution; empty draws it uniformly at random and
 *             normalizes, as the reference does
 * @param n    number of steps
 * @param rng  generator, advanced by the call
 */
template <class T>
CtmcPath<T> ctmc_simulate(const Matrix<T>& Q, const std::vector<T>& pi0, std::size_t n,
                          pfqn::McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_simulate requires transcendental arithmetic: the holding times are "
                  "exponential deviates drawn as -mean * log(u)");
    const std::size_t m = Q.rows();
    if (Q.cols() != m) throw InputError("ctmc_simulate: the generator is not square");
    if (m == 0) throw InputError("ctmc_simulate: empty generator");
    const T zero = num_traits<T>::from_int(0);

    // Initial distribution: uniform random and normalized when not supplied,
    // which is what `r = rand(length(Q),1); pi0 = r/sum(r)` does.
    std::vector<T> p0;
    if (pi0.empty()) {
        p0.resize(m);
        T tot = zero;
        for (std::size_t i = 0; i < m; ++i) {
            p0[i] = pfqn::mc_uniform<T>(rng);
            tot += p0[i];
        }
        if (tot == zero) throw NumericError("ctmc_simulate: degenerate random initial distribution");
        for (std::size_t i = 0; i < m; ++i) p0[i] /= tot;
    } else {
        if (pi0.size() != m)
            throw InputError("ctmc_simulate: pi0 has the wrong length for the generator");
        T tot = zero;
        for (std::size_t i = 0; i < m; ++i) {
            if (pi0[i] < zero) throw InputError("ctmc_simulate: pi0 has a negative entry");
            tot += pi0[i];
        }
        if (tot == zero) throw InputError("ctmc_simulate: pi0 has zero total mass");
        p0 = pi0;
        for (std::size_t i = 0; i < m; ++i) p0[i] /= tot;
    }

    // Row-normalized cumulative jump probabilities over the off-diagonal.
    Matrix<T> F(m, m, zero);
    for (std::size_t i = 0; i < m; ++i) {
        T run = zero;
        for (std::size_t j = 0; j < m; ++j) {
            if (j != i) {
                if (Q(i, j) < zero)
                    throw InputError("ctmc_simulate: negative off-diagonal rate in the generator");
                run += Q(i, j);
            }
            F(i, j) = run;
        }
        if (run == zero)
            throw NumericError(
                "ctmc_simulate: state " + std::to_string(i) +
                " is absorbing (no outgoing rate), so the sample path cannot be continued; "
                "the reference silently jumps to state 1 here with an infinite holding time");
        for (std::size_t j = 0; j < m; ++j) F(i, j) /= run;
    }

    // Correct inverse transform for the initial state; see the header note.
    std::size_t st = m - 1;
    {
        const T u = pfqn::mc_uniform<T>(rng);
        T run = zero;
        for (std::size_t i = 0; i < m; ++i) {
            run += p0[i];
            if (u < run) {
                st = i;
                break;
            }
        }
    }

    CtmcPath<T> path;
    path.states.reserve(n);
    path.sojourn.reserve(n);
    using std::log;
    for (std::size_t k = 0; k < n; ++k) {
        path.states.push_back(st);
        const T rate = -Q(st, st);
        if (!(rate > zero))
            throw NumericError("ctmc_simulate: state " + std::to_string(st) +
                               " has a non-negative diagonal, the holding time is undefined");
        // Exponential of mean 1/rate by inverse transform. u is drawn on
        // [0,1), so 1-u is in (0,1] and the logarithm is always finite.
        const T u = pfqn::mc_uniform<T>(rng);
        path.sojourn.push_back(T(-log(T(num_traits<T>::from_int(1) - u)) / rate));

        const T v = pfqn::mc_uniform<T>(rng);
        std::size_t nxt = m - 1;
        for (std::size_t j = 0; j < m; ++j) {
            if (v < F(st, j)) {
                nxt = j;
                break;
            }
        }
        st = nxt;
    }
    return path;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SIMULATE_H
