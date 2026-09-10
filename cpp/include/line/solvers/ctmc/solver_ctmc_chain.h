/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `matlab/src/solvers/CTMC/solver_ctmc_chain.m` and
 * `solver_ctmc_chain_transient.m`: SolverCTMC applied to a USER-SUPPLIED Markov
 * chain rather than to a queueing network. The reference dispatches on the
 * class of its argument, a MarkovProcess (CTMC, generator Q) or a MarkovChain
 * (DTMC, transition matrix P).
 *
 * These two functions read only `getGenerator` / `getTransMat` and
 * `stateSpace`, so `MarkovChainModel` carries exactly that. The rest of the
 * MarkovProcess / MarkovChain class surface -- `toEmbedded`, `toDTMC`,
 * `aggregate`, `stochComp`, `timeAverage`, `sens`, `hittingTime`, `sample` and
 * the two constructors from a random draw or an observed trajectory -- lives in
 * `lang/processes/markov_chain.h`, which is where the model type is declared
 * now; the alias below keeps the name this header introduced.
 *
 * NOTE that `solver_ctmc_chain` and `processes::chain_solve` DELIBERATELY
 * disagree on how to solve. This one ports `solver_ctmc_chain.m`, which tries
 * the primary solver and falls back to the reducible one on a failed validity
 * test; the class method ports `MarkovProcess.solve`, which takes the reducible
 * one outright for a numeric matrix. Do not "align" them.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_CHAIN_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_CHAIN_H

#include <chrono>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/lang/lang_types.h"
#include "line/lang/processes/markov_chain.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/num/number.h"

namespace line {
namespace ctmc {

using lang::GlobalConstants;

/** A user-supplied chain: a MarkovProcess when `discrete` is false, else a MarkovChain. */
template <class T>
using MarkovChainModel = lang::processes::MarkovChainModel<T>;

template <class T>
struct CtmcChainSolution {
    std::vector<T> pi;      ///< stationary distribution, length n
    Matrix<T> infgen;       ///< Q for a CTMC, the uniformized P-I for a DTMC
    Matrix<T> state_space;  ///< the chain's space, or the state indices when it carries none
    double runtime = 0.0;   ///< seconds
};

namespace chain_detail {

/**
 * Port of the reference's `solver_ctmc_chain_isvalid`.
 *
 * Reject a solution the primary solver could not produce on a reducible chain,
 * so that the reducible fallback is used instead. A NaN or a negative mass is
 * not a near miss to be patched: it says the primary solve did not apply.
 */
template <class T>
bool chain_isvalid(const std::vector<T>& pi, std::size_t n) {
    if (pi.size() != n) return false;
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double v = num_traits<T>::to_double(pi[i]);
        if (!std::isfinite(v)) return false;
        if (v < -GlobalConstants::FineTol) return false;
        total += v;
    }
    return std::abs(total - 1.0) <= std::sqrt(GlobalConstants::FineTol);
}

/** The chain's own space, or the 1-based state indices as a column. */
template <class T>
Matrix<T> chain_space(const MarkovChainModel<T>& chain, std::size_t n) {
    if (chain.state_space.rows() > 0) return chain.state_space;
    Matrix<T> idx(n, 1);
    for (std::size_t i = 0; i < n; ++i) idx(i, 0) = num_traits<T>::from_int(static_cast<long long>(i + 1));
    return idx;
}

}  // namespace chain_detail

/**
 * Steady-state analysis of a user-supplied Markov chain.
 *
 * `infgen` is Q for a CTMC and the uniformized generator P-I for a DTMC, which
 * carries the same stationary vector; returning it under one name is what lets
 * a caller treat the two chain kinds alike.
 */
template <class T>
CtmcChainSolution<T> solver_ctmc_chain(const MarkovChainModel<T>& chain) {
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    CtmcChainSolution<T> out;
    const std::size_t n = chain.mat.rows();
    if (n == 0 || chain.mat.cols() != n)
        throw InputError("solver_ctmc_chain: the chain matrix is empty or not square");

    if (chain.discrete) {
        const Matrix<T>& P = chain.mat;
        out.pi = mc::dtmc_solve(P);
        if (!chain_detail::chain_isvalid(out.pi, n)) out.pi = mc::dtmc_solve_reducible(P).pi;
        out.infgen = P;
        for (std::size_t i = 0; i < n; ++i) out.infgen(i, i) -= num_traits<T>::from_int(1);
    } else {
        const Matrix<T>& Q = chain.mat;
        out.pi = mc::ctmc_solve(Q);
        if (!chain_detail::chain_isvalid(out.pi, n)) out.pi = mc::ctmc_solve_reducible(Q).pi;
        out.infgen = Q;
    }

    out.state_space = chain_detail::chain_space(chain, n);
    out.runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

template <class T>
struct CtmcChainTransientSolution {
    Matrix<T> pi_t;     ///< one row per time point
    std::vector<T> t;   ///< the time points; integer STEPS for a DTMC
};

/**
 * Transient distribution of a user-supplied Markov chain over [t0,t1].
 *
 * For a MarkovProcess the Kolmogorov forward equations are integrated from
 * `pi0`; for a MarkovChain the distribution is advanced one step per unit of
 * time, so the returned `t` holds the integer steps within the timespan.
 *
 * @param pi0 initial distribution; empty for the uniform one
 */
template <class T>
CtmcChainTransientSolution<T> solver_ctmc_chain_transient(const MarkovChainModel<T>& chain,
                                                          const std::vector<T>& pi0in, const T& t0in,
                                                          const T& t1) {
    const std::size_t n = chain.mat.rows();
    if (n == 0 || chain.mat.cols() != n)
        throw InputError("solver_ctmc_chain_transient: the chain matrix is empty or not square");

    std::vector<T> pi0 = pi0in;
    if (pi0.empty()) {
        pi0.assign(n, num_traits<T>::from_double(1.0 / static_cast<double>(n)));
    } else if (pi0.size() != n) {
        throw InputError("solver_ctmc_chain_transient: the initial distribution has the wrong length");
    } else {
        double total = 0.0;
        for (std::size_t i = 0; i < n; ++i) total += num_traits<T>::to_double(pi0[i]);
        if (std::abs(total - 1.0) > GlobalConstants::FineTol)
            throw InputError("solver_ctmc_chain_transient: the initial distribution must sum to one");
    }

    const double d1 = num_traits<T>::to_double(t1);
    if (!std::isfinite(d1))
        throw InputError(
            "solver_ctmc_chain_transient: a finite timespan is required, e.g. --timespan 0,T");
    double d0 = num_traits<T>::to_double(t0in);
    // An infinite LOWER end is the "no start declared" spelling, not an error.
    if (!std::isfinite(d0)) d0 = 0.0;
    const T t0 = num_traits<T>::from_double(d0);

    CtmcChainTransientSolution<T> out;
    if (chain.discrete) {
        const Matrix<T>& P = chain.mat;
        const long long k0 = static_cast<long long>(std::ceil(d0));
        const long long k1 = static_cast<long long>(std::floor(d1));
        if (k1 < k0)
            throw InputError("solver_ctmc_chain_transient: the timespan contains no integer step of "
                             "the DTMC");
        const std::size_t steps = static_cast<std::size_t>(k1 - k0 + 1);

        // pi0 * P^k0, accumulated one step at a time: the reference forms the
        // matrix power, but the row vector needs only the vector-matrix product
        // and that is what keeps the cost linear in k0 rather than cubic.
        std::vector<T> pik = pi0;
        const T zero = num_traits<T>::from_int(0);
        std::vector<T> next(n);
        for (long long k = 0; k < k0; ++k) {
            for (std::size_t j = 0; j < n; ++j) {
                T acc = zero;
                for (std::size_t i = 0; i < n; ++i) acc += pik[i] * P(i, j);
                next[j] = acc;
            }
            pik = next;
        }

        out.pi_t = Matrix<T>(steps, n);
        out.t.resize(steps);
        for (std::size_t s = 0; s < steps; ++s) {
            out.t[s] = num_traits<T>::from_int(k0 + static_cast<long long>(s));
            for (std::size_t j = 0; j < n; ++j) out.pi_t(s, j) = pik[j];
            for (std::size_t j = 0; j < n; ++j) {
                T acc = zero;
                for (std::size_t i = 0; i < n; ++i) acc += pik[i] * P(i, j);
                next[j] = acc;
            }
            pik = next;
        }
    } else {
        const mc::TransientResult<T> r = mc::ctmc_transient(chain.mat, pi0, t0, t1);
        out.pi_t = r.pi;
        out.t = r.t;
    }
    return out;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_CHAIN_H
