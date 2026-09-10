/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_PROCESSES_MARKOV_CHAIN_H
#define LINE_LANG_PROCESSES_MARKOV_CHAIN_H

/**
 * The MarkovProcess / MarkovChain object surface.
 *
 * Port of `matlab/src/lang/processes/MarkovProcess.m` (a CTMC, carrying a
 * generator) and `MarkovChain.m` (a DTMC, carrying a transition matrix). Both
 * are thin objects over `api/mc`: the algorithms already lived there before
 * this header did, and what was missing was the surface that names them and
 * fixes which one each method calls. That distinction is the whole content of
 * the two classes and it is not cosmetic -- `solve` and `transient` each pick a
 * DIFFERENT primitive for the two chain kinds, and `toDTMC` and `toEmbedded`
 * both turn a CTMC into a DTMC while disagreeing about what the result means.
 *
 * ONE MODEL TYPE FOR BOTH KINDS. `MarkovChainModel<T>::discrete` says which
 * class the object would be in MATLAB. The reference dispatches on the class,
 * so every function here dispatches on that flag, and the constructors'
 * normalization (`ctmc_makeinfgen` / `dtmc_makestochastic`) is applied by the
 * factories exactly as the two constructors apply it. This type was previously
 * declared inside `solvers/ctmc/solver_ctmc_chain.h`, which aliases it now.
 *
 * WHAT IS DELIBERATELY ABSENT:
 *  - `plot` and `plot3`. They call graphViz4Matlab and MATLAB's `digraph`, a
 *    rendering layer this tree has no counterpart for. A caller wanting the
 *    graph has `mat` and `state_space` and can emit whatever format it wants.
 *  - `getGenerator` / `getTransMat` / `setStateSpace`, which are field access:
 *    `mat` and `state_space` are public members.
 *  - `isfinite`. It is a flag the reference stores and never reads.
 *
 * TWO PLACES WHERE THE REFERENCE IS NOT REPRODUCIBLE, and how this port stands:
 *  - `toMarkovChain` with no argument picks `q = max|Q| + rand`, so the SAME
 *    chain uniformizes to a different P on every call. Any `q > max|Q|` is a
 *    valid uniformization rate and all of them carry the same stationary law,
 *    but the returned matrix is not the same object, so a golden recorded
 *    against it would be noise. The default here is deterministic,
 *    `q = max|Q| * (1 + 1/16)`, and the rate is an explicit parameter for a
 *    caller that wants the reference's own draw.
 *  - `MarkovChain.sample` draws its own uniform initial law before simulating.
 *    Here the initial law is a parameter, empty meaning that same random draw,
 *    so the caller can make the path reproducible without a fork of the
 *    primitive.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_foxglynn.h"
#include "line/api/mc/ctmc_isfeasible.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_multi.h"
#include "line/api/mc/ctmc_passage.h"
#include "line/api/mc/ctmc_rand.h"
#include "line/api/mc/ctmc_relsolve.h"
#include "line/api/mc/ctmc_sens.h"
#include "line/api/mc/ctmc_simulate.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/api/mc/ctmc_timereverse.h"
#include "line/api/mc/ctmc_uniformization.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_rand.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/api/mc/dtmc_transient.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace lang {
namespace processes {

/** A user-supplied chain: a MarkovProcess when `discrete` is false, else a MarkovChain. */
template <class T>
struct MarkovChainModel {
    Matrix<T> mat;          ///< generator Q (CTMC) or transition matrix P (DTMC)
    bool discrete = false;  ///< true for a MarkovChain
    Matrix<T> state_space;  ///< optional; empty means the chain carries none

    MarkovChainModel() {}

    /** Mirrors the constructors: a generator is closed, a transition matrix normalized. */
    static MarkovChainModel<T> process(const Matrix<T>& infgen,
                                       const Matrix<T>& space = Matrix<T>()) {
        MarkovChainModel<T> m;
        m.mat = mc::ctmc_makeinfgen(infgen);
        m.discrete = false;
        m.state_space = space;
        return m;
    }

    static MarkovChainModel<T> chain(const Matrix<T>& transmat,
                                     const Matrix<T>& space = Matrix<T>()) {
        MarkovChainModel<T> m;
        m.mat = mc::dtmc_makestochastic(transmat);
        m.discrete = true;
        m.state_space = space;
        return m;
    }

    std::size_t order() const { return mat.rows(); }
};

namespace detail {

/** Every entry point below reads a square matrix; say so once, by name. */
template <class T>
void require_square(const MarkovChainModel<T>& m, const char* who) {
    if (m.mat.rows() == 0 || m.mat.cols() != m.mat.rows())
        throw InputError(std::string(who) + ": the chain matrix is empty or not square");
}

template <class T>
void require_kind(const MarkovChainModel<T>& m, bool discrete, const char* who) {
    if (m.discrete != discrete)
        throw InputError(std::string(who) + ": this is a method of " +
                         (discrete ? "MarkovChain (a DTMC); the object is a MarkovProcess"
                                   : "MarkovProcess (a CTMC); the object is a MarkovChain"));
}

/** The uniform law over n states, the default every `pi0` argument falls back to. */
template <class T>
std::vector<T> uniform_law(std::size_t n) {
    return std::vector<T>(n, num_traits<T>::from_int(1) / num_traits<T>::from_int(
                                                             static_cast<int>(n)));
}

template <class T>
std::vector<T> law_or_uniform(const std::vector<T>& pi0, std::size_t n, const char* who) {
    if (pi0.empty()) return uniform_law<T>(n);
    if (pi0.size() != n)
        throw InputError(std::string(who) + ": the initial distribution has the wrong length");
    return pi0;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Conversions between the two chain kinds
// ---------------------------------------------------------------------------

/**
 * `MarkovChain.toMarkovProcess` / `toCTMC`: read the DTMC as a CTMC with unit
 * exit rates, Q = P - I. The stationary law is preserved, since P and P - I
 * have the same left null structure up to the shift.
 */
template <class T>
MarkovChainModel<T> to_markov_process(const MarkovChainModel<T>& m) {
    detail::require_square(m, "to_markov_process");
    detail::require_kind(m, true, "to_markov_process");
    Matrix<T> Q = m.mat;
    const std::size_t n = Q.rows();
    for (std::size_t i = 0; i < n; ++i) Q(i, i) -= num_traits<T>::from_int(1);
    MarkovChainModel<T> out = MarkovChainModel<T>::process(Q, m.state_space);
    return out;
}

/** The uniformization rate this port uses when the caller names none. */
template <class T>
T default_uniformization_rate(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    T qmax = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const T a = num_abs(Q(i, j));
            if (a > qmax) qmax = a;
        }
    // Strictly above max|Q|, as the reference's `+rand` guarantees, but fixed.
    return T(qmax * num_traits<T>::from_double(17.0 / 16.0));
}

/**
 * `MarkovProcess.toMarkovChain` / `toDTMC`: the UNIFORMIZED chain, P = Q/q + I.
 *
 * This is the conversion that PRESERVES the stationary distribution, and it is
 * the one to reach for when the question is about long-run behaviour. Contrast
 * `to_embedded`, which does not.
 */
template <class T>
MarkovChainModel<T> to_markov_chain(const MarkovChainModel<T>& m, const T& q) {
    detail::require_square(m, "to_markov_chain");
    detail::require_kind(m, false, "to_markov_chain");
    if (!(q > num_traits<T>::from_int(0)))
        throw InputError("to_markov_chain: the uniformization rate must be positive");
    const std::size_t n = m.mat.rows();
    Matrix<T> P = m.mat;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) P(i, j) = T(P(i, j) / q);
    for (std::size_t i = 0; i < n; ++i) P(i, i) += num_traits<T>::from_int(1);
    return MarkovChainModel<T>::chain(P, m.state_space);
}

template <class T>
MarkovChainModel<T> to_markov_chain(const MarkovChainModel<T>& m) {
    detail::require_square(m, "to_markov_chain");
    detail::require_kind(m, false, "to_markov_chain");
    return to_markov_chain(m, default_uniformization_rate(m.mat));
}

/** `toDTMC`, the backwards-compatible alias of `toMarkovChain`. */
template <class T>
MarkovChainModel<T> to_dtmc(const MarkovChainModel<T>& m) {
    return to_markov_chain(m);
}

template <class T>
MarkovChainModel<T> to_dtmc(const MarkovChainModel<T>& m, const T& q) {
    return to_markov_chain(m, q);
}

/**
 * `MarkovProcess.toEmbedded`: the JUMP CHAIN, the DTMC of the states visited at
 * transition epochs.
 *
 * IT DOES NOT PRESERVE THE STATIONARY LAW, and that is the point of having it
 * separate from `to_markov_chain`. Dividing each off-diagonal row by the exit
 * rate throws away how long the chain lingers, so a state with a fast exit rate
 * is visited as often as a slow one and weighs the same here while weighing far
 * less in the CTMC. An absorbing state (exit rate zero) has no next jump, and
 * the reference makes it absorbing in the jump chain too rather than leaving an
 * all-zero row that `dtmc_makestochastic` would have to invent a law for.
 */
template <class T>
MarkovChainModel<T> to_embedded(const MarkovChainModel<T>& m) {
    detail::require_square(m, "to_embedded");
    detail::require_kind(m, false, "to_embedded");
    const std::size_t n = m.mat.rows();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> P = m.mat;
    for (std::size_t i = 0; i < n; ++i) {
        const T exit_rate = T(zero - m.mat(i, i));
        P(i, i) = zero;
        if (exit_rate > zero) {
            for (std::size_t j = 0; j < n; ++j) P(i, j) = T(P(i, j) / exit_rate);
        } else {
            P(i, i) = num_traits<T>::from_int(1);
        }
    }
    return MarkovChainModel<T>::chain(P, m.state_space);
}

/** `toTimeReversed` for either kind: the chain run backwards in time. */
template <class T>
MarkovChainModel<T> to_time_reversed(const MarkovChainModel<T>& m) {
    detail::require_square(m, "to_time_reversed");
    if (m.discrete)
        return MarkovChainModel<T>::chain(mc::dtmc_timereverse(m.mat), m.state_space);
    return MarkovChainModel<T>::process(mc::ctmc_timereverse(m.mat), m.state_space);
}

// ---------------------------------------------------------------------------
// Stationary analysis
// ---------------------------------------------------------------------------

/**
 * `MarkovProcess.solve` / `MarkovChain.solve`.
 *
 * THE REDUCIBLE SOLVER IS THE NUMERIC PATH IN BOTH CLASSES, not a fallback:
 * the reference reserves the plain `ctmc_solve` / `dtmc_solve` for a SYMBOLIC
 * matrix, where the reducible variant's component decomposition has nothing to
 * decide. This port has no symbolic element type, so the reducible one is what
 * every call takes. `solver_ctmc_chain` in `solvers/ctmc` deliberately does the
 * other thing -- primary first, reducible on a failed validity test -- because
 * it ports `solver_ctmc_chain.m`, not the class method, and those two disagree.
 */
template <class T>
std::vector<T> chain_solve(const MarkovChainModel<T>& m) {
    detail::require_square(m, "chain_solve");
    if (m.discrete) return mc::dtmc_solve_reducible(m.mat).pi;
    return mc::ctmc_solve_reducible(m.mat).pi;
}

/**
 * `MarkovProcess.solveRelative`: the equilibrium vector normalized so that
 * `refstate` carries one, which exists even where the normalizing constant does
 * not. `refstate` is 0-based here and 1-based in the reference.
 */
template <class T>
std::vector<T> solve_relative(const MarkovChainModel<T>& m, std::size_t refstate = 0) {
    detail::require_square(m, "solve_relative");
    detail::require_kind(m, false, "solve_relative");
    return mc::ctmc_relsolve(m.mat, refstate);
}

/** What `getProbState` returns: the probability and the two determinants behind it. */
template <class T>
struct ProbStateResult {
    T pi_i;  ///< probability of the state
    T num;   ///< determinant of the numerator matrix
    T den;   ///< determinant of the denominator matrix
};

/**
 * `MarkovProcess.getProbState`: the probability of ONE state by Cramer's rule.
 *
 * Column 0 of the generator is replaced by ones, which imposes the
 * normalization in place of the column the balance equations make redundant;
 * the numerator matrix additionally zeroes row `i` and puts a one back in its
 * first entry. The quotient of the two determinants is the probability.
 *
 * WHY A DETERMINANT AND NOT A SOLVE. This exists so that the probability of one
 * state can be written as a RATIO OF POLYNOMIALS in the generator's entries,
 * which is what makes it useful symbolically in the reference. Numerically a
 * full solve is cheaper and better conditioned, and `chain_solve` is that; this
 * one is kept faithful because a caller reaching for it wants `num` and `den`
 * separately.
 */
template <class T>
ProbStateResult<T> get_prob_state(const MarkovChainModel<T>& m, std::size_t i) {
    detail::require_square(m, "get_prob_state");
    detail::require_kind(m, false, "get_prob_state");
    const std::size_t n = m.mat.rows();
    if (i >= n) throw InputError("get_prob_state: state index is out of range");
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);

    Matrix<T> Q = m.mat;
    for (std::size_t r = 0; r < n; ++r) Q(r, 0) = one;
    Matrix<T> Qi = Q;
    for (std::size_t j = 0; j < n; ++j) Qi(i, j) = zero;
    Qi(i, 0) = one;

    ProbStateResult<T> out;
    out.num = lu_det(Qi);
    out.den = lu_det(Q);
    if (out.den == zero)
        throw NumericError(
            "get_prob_state: the normalized generator is singular, so Cramer's rule has no "
            "quotient; the chain is reducible -- use chain_solve, which decomposes it");
    out.pi_i = T(out.num / out.den);
    return out;
}

/** Row index of `state` in the chain's state space, or `n` when it carries none. */
template <class T>
std::size_t match_state(const MarkovChainModel<T>& m, const std::vector<T>& state) {
    const std::size_t rows = m.state_space.rows(), cols = m.state_space.cols();
    if (rows == 0)
        throw InputError("match_state: the chain carries no state space to match against");
    if (state.size() != cols)
        throw InputError("match_state: the state has the wrong number of columns");
    for (std::size_t r = 0; r < rows; ++r) {
        bool hit = true;
        for (std::size_t c = 0; c < cols && hit; ++c) hit = (m.state_space(r, c) == state[c]);
        if (hit) return r;
    }
    return rows;
}

/** `getProbState` addressed by the state itself rather than by its index. */
template <class T>
ProbStateResult<T> get_prob_state(const MarkovChainModel<T>& m, const std::vector<T>& state) {
    const std::size_t i = match_state(m, state);
    if (i >= m.state_space.rows())
        throw InputError("get_prob_state: the state is not in the chain's state space");
    return get_prob_state(m, i);
}

/** `isFeasible`: a valid generator, or a stochastic transition matrix. */
template <class T>
bool is_feasible(const MarkovChainModel<T>& m) {
    detail::require_square(m, "is_feasible");
    if (m.discrete) return mc::dtmc_isfeasible(m.mat) > 0;
    return mc::ctmc_isfeasible(m.mat);
}

// ---------------------------------------------------------------------------
// Transient analysis
// ---------------------------------------------------------------------------

/** What the CTMC `transient` returns: the law at t and the truncation it used. */
template <class T>
struct TransientAtResult {
    std::vector<T> pi;  ///< distribution at time t
    std::size_t kmax;   ///< Poisson terms used (the right truncation point for 'foxglynn')
};

/**
 * `MarkovProcess.transient`: the law at ONE time t, by uniformization.
 *
 * `method` is `"unif"` (Jensen, the default) or `"foxglynn"`, whose weights are
 * built by the Fox-Glynn recursion instead of by evaluating Poisson terms, so
 * it survives a `q t` large enough to underflow them. Note that `t` is a TIME
 * here; the DTMC counterpart takes a step count, which is why the two are
 * separate functions rather than one dispatching on `discrete`.
 */
template <class T>
TransientAtResult<T> chain_transient_at(const MarkovChainModel<T>& m, const std::vector<T>& pi0in,
                                        const T& t, const std::string& method = "unif") {
    detail::require_square(m, "chain_transient_at");
    detail::require_kind(m, false, "chain_transient_at");
    const std::size_t n = m.mat.rows();
    const std::vector<T> pi0 = detail::law_or_uniform(pi0in, n, "chain_transient_at");

    TransientAtResult<T> out;
    if (method == "foxglynn") {
        const mc::FoxGlynnResult<T> r = mc::ctmc_foxglynn(pi0, m.mat, t);
        out.pi = r.pi;
        out.kmax = r.right < 0 ? 0 : static_cast<std::size_t>(r.right);
    } else if (method == "unif" || method.empty()) {
        const mc::UniformizationResult<T> r = mc::ctmc_uniformization(pi0, m.mat, t);
        out.pi = r.pi;
        out.kmax = r.kmax;
    } else {
        throw InputError("chain_transient_at: unknown transient method '" + method +
                         "'; the reference offers 'unif' and 'foxglynn'");
    }
    return out;
}

/**
 * `MarkovChain.transient`: the law at every step 0..steps, one row per step.
 *
 * The DTMC's clock is the step count, so unlike the CTMC method this returns
 * the whole trajectory and takes no time argument.
 */
template <class T>
Matrix<T> chain_transient_steps(const MarkovChainModel<T>& m, const std::vector<T>& pi0in,
                                std::size_t steps = 1) {
    detail::require_square(m, "chain_transient_steps");
    detail::require_kind(m, true, "chain_transient_steps");
    const std::size_t n = m.mat.rows();
    const std::vector<T> pi0 = detail::law_or_uniform(pi0in, n, "chain_transient_steps");
    return mc::dtmc_transient(m.mat, pi0, steps);
}

/**
 * `MarkovChain.transientUnif`: the DTMC read as the randomized image of a CTMC,
 * so `t` is CONTINUOUS here where `chain_transient_steps` counts steps. The two
 * answer different questions about the same matrix and the reference keeps both.
 */
template <class T>
TransientAtResult<T> chain_transient_unif(const MarkovChainModel<T>& m,
                                          const std::vector<T>& pi0in, const T& t) {
    detail::require_square(m, "chain_transient_unif");
    detail::require_kind(m, true, "chain_transient_unif");
    const std::size_t n = m.mat.rows();
    const std::vector<T> pi0 = detail::law_or_uniform(pi0in, n, "chain_transient_unif");
    const mc::UniformizationResult<T> r = mc::dtmc_uniformization(pi0, m.mat, t);
    TransientAtResult<T> out;
    out.pi = r.pi;
    out.kmax = r.kmax;
    return out;
}

/** What `timeAverage` returns. */
template <class T>
struct TimeAverageOut {
    std::vector<T> pi_time_avg;  ///< time-averaged law over [0, t]
    std::vector<T> pi_exit;      ///< law at t
    std::size_t kmax;
};

/** `MarkovProcess.timeAverage`: the law averaged over [0,t], and its endpoint. */
template <class T>
TimeAverageOut<T> time_average(const MarkovChainModel<T>& m, const std::vector<T>& pi0in,
                               const T& t) {
    detail::require_square(m, "time_average");
    detail::require_kind(m, false, "time_average");
    const std::size_t n = m.mat.rows();
    const std::vector<T> pi0 = detail::law_or_uniform(pi0in, n, "time_average");
    const mc::TimeAverageResult<T> r = mc::ctmc_timeaverage(pi0, m.mat, t);
    TimeAverageOut<T> out;
    out.pi_time_avg = r.piTimeAvg;
    out.pi_exit = r.piExit;
    out.kmax = r.kmax;
    return out;
}

/**
 * `MarkovProcess.sens`: the derivative of the stationary law with respect to a
 * scalar parameter, given the derivative `dQ` of the generator. The reference
 * feeds it its own `solve()`, so this one does too.
 */
template <class T>
std::vector<T> chain_sens(const MarkovChainModel<T>& m, const Matrix<T>& dQ) {
    detail::require_square(m, "chain_sens");
    detail::require_kind(m, false, "chain_sens");
    return mc::ctmc_sens(m.mat, dQ, chain_solve(m));
}

// ---------------------------------------------------------------------------
// Aggregation-disaggregation
// ---------------------------------------------------------------------------

/** What `aggregate` returns: the approximate law and the two NCD indices. */
template <class T>
struct AggregateResult {
    std::vector<T> p;  ///< approximate stationary vector, ORIGINAL state ordering
    T eps;             ///< nearly-complete-decomposability index of the partition
    T epsMAX;          ///< the largest index for which the approximation is meant to hold
};

/**
 * `MarkovProcess.aggregate`: aggregation-disaggregation over a macrostate
 * partition.
 *
 * `method` is `"courtois"` (the default; `param` is the randomization rate q),
 * `"kms"` or `"takahashi"` (`param` is the sweep count, default 10), or
 * `"multi"`, which needs the SECOND-LEVEL partition and therefore does not go
 * through `param` at all -- it is refused here without `MSS`, as the reference
 * refuses it.
 *
 * READ `eps` BEFORE THE ANSWER. These are approximations whose error is
 * governed by how nearly decomposable the partition is; `eps > epsMAX` means
 * the partition does not justify the method, and the vector is returned anyway
 * because the reference returns it. It is a diagnostic, not a gate.
 *
 * `MS` is 0-based here and 1-based in the reference.
 */
template <class T>
AggregateResult<T> aggregate(const MarkovChainModel<T>& m,
                             const std::vector<std::vector<std::size_t>>& MS,
                             const std::string& method = "courtois",
                             const T* param = nullptr) {
    detail::require_square(m, "aggregate");
    detail::require_kind(m, false, "aggregate");
    if (MS.empty()) throw InputError("aggregate: the macrostate partition is empty");

    AggregateResult<T> out;
    if (method == "courtois" || method.empty()) {
        const mc::CourtoisResult<T> r =
            param ? mc::ctmc_courtois(m.mat, MS, *param) : mc::ctmc_courtois(m.mat, MS);
        out.p = r.p;
        out.eps = r.eps;
        out.epsMAX = r.epsMAX;
    } else if (method == "kms") {
        const std::size_t steps =
            param ? static_cast<std::size_t>(num_traits<T>::to_double(*param)) : 10;
        const mc::KmsResult<T> r = mc::ctmc_kms(m.mat, MS, steps);
        out.p = r.p;
        out.eps = r.eps;
        out.epsMAX = r.epsMAX;
    } else if (method == "takahashi") {
        const std::size_t steps =
            param ? static_cast<std::size_t>(num_traits<T>::to_double(*param)) : 10;
        const mc::TakahashiResult<T> r = mc::ctmc_takahashi(m.mat, MS, steps);
        out.p = r.p;
        out.eps = r.eps;
        out.epsMAX = r.epsMAX;
    } else if (method == "multi") {
        throw InputError(
            "aggregate: the 'multi' method requires the second-level partition MSS; call the "
            "aggregate_multi overload, which takes it");
    } else {
        throw InputError("aggregate: unknown aggregation method '" + method + "'");
    }
    return out;
}

/**
 * The `"multi"` arm of `aggregate`, separated because its parameter is a
 * PARTITION OF THE PARTITION and not a scalar. `MSS` partitions the macrostate
 * indices 0..|MS|-1.
 */
template <class T>
AggregateResult<T> aggregate_multi(const MarkovChainModel<T>& m,
                                   const std::vector<std::vector<std::size_t>>& MS,
                                   const std::vector<std::vector<std::size_t>>& MSS) {
    detail::require_square(m, "aggregate_multi");
    detail::require_kind(m, false, "aggregate_multi");
    if (MS.empty() || MSS.empty())
        throw InputError("aggregate_multi: both partitions must be non-empty");
    const mc::MultiResult<T> r = mc::ctmc_multi(m.mat, MS, MSS);
    AggregateResult<T> out;
    out.p = r.p;
    out.eps = r.eps;
    out.epsMAX = r.epsMAX;
    return out;
}

// ---------------------------------------------------------------------------
// Stochastic complementation
// ---------------------------------------------------------------------------

/** What `stochCompFull` returns for a CTMC; a DTMC fills the same blocks from P. */
template <class T>
struct StochCompOut {
    Matrix<T> S;    ///< the complement on the selected states
    Matrix<T> A11;  ///< the four blocks of the partitioned matrix
    Matrix<T> A12;
    Matrix<T> A21;
    Matrix<T> A22;
    Matrix<T> T12;  ///< the return-path term, so that S = A11 + T12
};

/**
 * `stochComp` / `stochCompFull` for either kind.
 *
 * `I` is 0-based here and 1-based in the reference. An empty `I` takes the
 * reference's own default, the first half of the state space.
 *
 * THE DTMC ARM FILLS ONLY `S`. `dtmc_stochcomp` returns the complement alone,
 * as `MarkovChain.stochCompFull` reports blocks its own primitive does not
 * separate; the four blocks are left empty rather than reconstructed here,
 * which would be a different function under the same name.
 */
template <class T>
StochCompOut<T> stoch_comp_full(const MarkovChainModel<T>& m,
                                const std::vector<std::size_t>& I = std::vector<std::size_t>()) {
    detail::require_square(m, "stoch_comp_full");
    StochCompOut<T> out;
    if (m.discrete) {
        std::vector<std::size_t> keep = I;
        if (keep.empty()) {
            const std::size_t half = m.mat.rows() / 2;
            if (half == 0)
                throw InputError(
                    "stoch_comp_full: the default partition takes the first half of the state "
                    "space, which is empty on a chain of order one; pass I explicitly");
            for (std::size_t i = 0; i < half; ++i) keep.push_back(i);
        }
        out.S = mc::dtmc_stochcomp(m.mat, keep);
        return out;
    }
    const mc::StochCompResult<T> r =
        I.empty() ? mc::ctmc_stochcomp(m.mat) : mc::ctmc_stochcomp(m.mat, I);
    out.S = r.S;
    out.A11 = r.Q11;
    out.A12 = r.Q12;
    out.A21 = r.Q21;
    out.A22 = r.Q22;
    out.T12 = r.T12;
    return out;
}

/** `stochComp`: the complement alone. */
template <class T>
Matrix<T> stoch_comp(const MarkovChainModel<T>& m,
                     const std::vector<std::size_t>& I = std::vector<std::size_t>()) {
    return stoch_comp_full(m, I).S;
}

// ---------------------------------------------------------------------------
// Hitting times and sampling
// ---------------------------------------------------------------------------

/**
 * `hittingTime`: the mean time (CTMC) or step count (DTMC) to reach any state
 * in `target`, zero on the target set itself and infinite from a state that
 * cannot reach it. `target` is 0-based here and 1-based in the reference.
 */
template <class T>
std::vector<T> hitting_time(const MarkovChainModel<T>& m,
                            const std::vector<std::size_t>& target) {
    detail::require_square(m, "hitting_time");
    if (target.empty()) throw InputError("hitting_time: the target set is empty");
    if (m.discrete) return mc::dtmc_hitting_time(m.mat, target);
    return mc::ctmc_hitting_time(m.mat, target);
}

/** A sampled path: the states visited, with their holding times for a CTMC. */
template <class T>
struct ChainPath {
    std::vector<std::size_t> states;  ///< 0-based state index at each step
    std::vector<T> sojourn;           ///< holding times; empty for a DTMC, which has none
};

/**
 * `sample`: simulate `n` steps.
 *
 * The reference draws the initial state from a uniform law it randomizes itself
 * (`MarkovChain.sample`) or from the primitive's default (`MarkovProcess`).
 * Here `pi0` is explicit and empty reproduces that default, so a caller can
 * make the path reproducible by naming the law and seeding `gen`.
 */
template <class T, class Gen>
ChainPath<T> chain_sample(const MarkovChainModel<T>& m, const std::vector<T>& pi0in,
                          std::size_t n, Gen& gen) {
    detail::require_square(m, "chain_sample");
    ChainPath<T> out;
    const std::size_t order = m.mat.rows();
    if (m.discrete) {
        std::vector<T> pi0 = pi0in;
        if (pi0.empty()) {
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            pi0.resize(order);
            double total = 0.0;
            for (std::size_t i = 0; i < order; ++i) {
                const double u = unif(gen);
                pi0[i] = num_traits<T>::from_double(u);
                total += u;
            }
            for (std::size_t i = 0; i < order; ++i)
                pi0[i] = T(pi0[i] / num_traits<T>::from_double(total));
        } else if (pi0.size() != order) {
            throw InputError("chain_sample: the initial distribution has the wrong length");
        }
        out.states = mc::dtmc_simulate(m.mat, pi0, n, gen);
        return out;
    }
    if (!pi0in.empty() && pi0in.size() != order)
        throw InputError("chain_sample: the initial distribution has the wrong length");
    pfqn::McRng rng(gen());
    const mc::CtmcPath<T> p = mc::ctmc_simulate(m.mat, pi0in, n, rng);
    out.states = p.states;
    out.sojourn = p.sojourn;
    return out;
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

/** `MarkovProcess.rand`: a random generator of the given order. */
template <class T, class Gen>
MarkovChainModel<T> rand_process(std::size_t n, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    auto draw = [&]() { return unif(gen); };
    return MarkovChainModel<T>::process(mc::ctmc_rand<T>(n, draw));
}

/** `MarkovChain.rand`: a random transition matrix of the given order. */
template <class T, class Gen>
MarkovChainModel<T> rand_chain(std::size_t n, Gen& gen) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    auto draw = [&]() { return unif(gen); };
    return MarkovChainModel<T>::chain(mc::dtmc_rand<T>(n, draw));
}

/**
 * `MarkovChain.fromSampleSysAggr`: estimate a DTMC from an observed trajectory.
 *
 * `sample_state` holds one row per observation, the columns being the aggregate
 * state; the reference joins the per-node trajectories COLUMN-WISE first, which
 * is time alignment, so a caller must pass the joined matrix. Distinct rows
 * become the state space in first-appearance order, transition counts between
 * consecutive observations become the matrix, and `dtmc_makestochastic`
 * normalizes it.
 *
 * FIRST-APPEARANCE ORDER, NOT SORTED ORDER, is a deliberate departure: MATLAB's
 * `unique(...,'rows')` sorts, and reproducing a lexicographic sort over rows of
 * an arbitrary element type would be a second, unstated, definition of order.
 * The estimated chain is the same up to the permutation, and `state_space`
 * carries the labelling, so a caller reading the two together is unaffected.
 * A caller comparing raw matrix entries against MATLAB is, and should permute.
 */
template <class T>
MarkovChainModel<T> from_sample_sys_aggr(const Matrix<T>& sample_state) {
    const std::size_t obs = sample_state.rows(), cols = sample_state.cols();
    if (obs < 2)
        throw InputError(
            "from_sample_sys_aggr: at least two observations are needed, since the estimate is "
            "built from the transitions between consecutive ones");

    std::vector<std::size_t> hash(obs);
    std::vector<std::vector<T>> space;
    for (std::size_t r = 0; r < obs; ++r) {
        std::vector<T> row(cols);
        for (std::size_t c = 0; c < cols; ++c) row[c] = sample_state(r, c);
        std::size_t found = space.size();
        for (std::size_t s = 0; s < space.size(); ++s) {
            bool hit = true;
            for (std::size_t c = 0; c < cols && hit; ++c) hit = (space[s][c] == row[c]);
            if (hit) {
                found = s;
                break;
            }
        }
        if (found == space.size()) space.push_back(row);
        hash[r] = found;
    }

    const std::size_t n = space.size();
    Matrix<T> counts(n, n, num_traits<T>::from_int(0));
    for (std::size_t r = 1; r < obs; ++r)
        counts(hash[r - 1], hash[r]) += num_traits<T>::from_int(1);

    Matrix<T> ss(n, cols);
    for (std::size_t s = 0; s < n; ++s)
        for (std::size_t c = 0; c < cols; ++c) ss(s, c) = space[s][c];
    return MarkovChainModel<T>::chain(counts, ss);
}

}  // namespace processes
}  // namespace lang
}  // namespace line

#endif  // LINE_LANG_PROCESSES_MARKOV_CHAIN_H
