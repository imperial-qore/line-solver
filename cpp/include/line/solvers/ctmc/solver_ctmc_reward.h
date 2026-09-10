/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `solver_ctmc_reward.m` and the `@@SolverCTMC` reward surface
 * (`runRewardAnalyzer`, `getAvgReward`, `getTranReward`).
 *
 * TWO DIFFERENT QUANTITIES SHARE THE WORD "REWARD", and conflating them is the
 * trap this file exists to avoid:
 *
 *   steady state  E[r] = sum_s pi(s) r(s), a RATE -- the long-run average
 *                 reward earned per unit time.
 *   value function V^k(s), the reward ACCUMULATED over k uniformized steps
 *                 starting from s. It grows without bound in a recurrent chain,
 *                 because it is a total and not an average.
 *
 * `V` is therefore not "the transient version of E[r]" and does not converge to
 * it; its SLOPE does. The reference reports both and so does this port.
 *
 * UNIFORMIZATION IS WHAT MAKES THE VALUE ITERATION A CTMC ANSWER. The embedded
 * chain P = Q/q + I with q = max|diag(Q)| has the same stationary law as Q and
 * a uniform step of mean duration 1/q, so iteration index k maps to time k/q.
 * Any q at least as large as the maximum exit rate is valid; taking the maximum
 * is the tightest, hence the fastest-mixing, choice.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_REWARD_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_REWARD_H

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/** What the reward analyzer produces, per declared reward. */
template <class T>
struct CtmcReward {
    std::vector<std::string> names;
    std::vector<T> steady_state;            ///< E[r] per reward
    std::vector<Matrix<T>> V;               ///< V[r] is (Tmax+1 x nstates)
    std::vector<T> t;                       ///< iteration index / q
    Matrix<T> state_space_aggr;             ///< the rows the reward saw
    CtmcSolution<T> chain;
};

/**
 * Port of `solver_ctmc_reward.m`.
 *
 * @param tmax number of value-iteration steps; the reference's default is 1000
 * @param sn the refreshed network struct, carrying the reward definitions
 * @param opt CTMC options (state-space cutoff, tolerances, method)
 */
template <class T>
CtmcReward<T> solver_ctmc_reward(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                 std::size_t tmax = 1000) {
    if (sn.reward.empty())
        throw InputError(
            "solver_ctmc_reward: no rewards are defined; declare one with set_reward(name, fn) "
            "before asking for a reward analysis");

    CtmcReward<T> out;
    out.chain = solver_ctmc_analyzer(sn, opt);
    const std::size_t n = out.chain.chain.space.size();
    const std::size_t nr = sn.reward.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    out.state_space_aggr = ctmc_state_space_aggr(sn, out.chain.chain.space);

    // The reward vector of each declaration, evaluated once per state on the
    // aggregate row -- the same row `RewardState` wraps in the reference.
    Matrix<T> R(nr, n, zero);
    out.names.resize(nr);
    for (std::size_t r = 0; r < nr; ++r) {
        out.names[r] = sn.reward[r].name;
        for (std::size_t s = 0; s < n; ++s) {
            std::vector<T> row(out.state_space_aggr.cols());
            for (std::size_t c = 0; c < row.size(); ++c) row[c] = out.state_space_aggr(s, c);
            R(r, s) = sn.reward[r].fn(row);
        }
    }

    out.steady_state.assign(nr, zero);
    for (std::size_t r = 0; r < nr; ++r)
        for (std::size_t s = 0; s < n; ++s)
            out.steady_state[r] += T(out.chain.pi[s] * R(r, s));

    // Uniformization rate: the largest exit rate. A generator with none -- every
    // state absorbing -- would divide by zero, so it falls back to 1, which
    // leaves P = I and a value function that simply accumulates r(s).
    double q = 0;
    for (std::size_t s = 0; s < n; ++s)
        q = std::max(q, std::fabs(num_traits<T>::to_double(out.chain.chain.Q(s, s))));
    if (q == 0) q = 1.0;
    const T qq = num_traits<T>::from_double(q);

    Matrix<T> P(n, n, zero);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b)
            P(a, b) = a == b ? T(out.chain.chain.Q(a, b) / qq + one)
                             : T(out.chain.chain.Q(a, b) / qq);

    out.V.assign(nr, Matrix<T>(tmax + 1, n, zero));
    for (std::size_t r = 0; r < nr; ++r) {
        std::vector<T> prev(n, zero);
        for (std::size_t k = 1; k <= tmax; ++k) {
            std::vector<T> next(n, zero);
            // V^{k+1}(s) = r(s) + sum_s' P(s,s') V^k(s'), the reference's
            // `R(r,:) + v_prev*P'` -- note the TRANSPOSE, which makes it a
            // forward expectation over the successor rather than a backward
            // one over the predecessor.
            for (std::size_t s = 0; s < n; ++s) {
                T acc = R(r, s);
                for (std::size_t sp = 0; sp < n; ++sp) acc += T(P(s, sp) * prev[sp]);
                next[s] = acc;
            }
            for (std::size_t s = 0; s < n; ++s) out.V[r](k, s) = next[s];
            prev.swap(next);
        }
    }

    out.t.resize(tmax + 1);
    for (std::size_t k = 0; k <= tmax; ++k)
        out.t[k] = T(num_traits<T>::from_double(static_cast<double>(k)) / qq);
    return out;
}

/**
 * Port of `@@SolverCTMC/getTranReward`: E[r(X(t))] = sum_s pi_t(s) r(s).
 *
 * NOT THE VALUE FUNCTION `V` above. This is the expected reward RATE at time t,
 * which converges to the steady-state E[r]; `V` is the reward accumulated over
 * k uniformized steps and diverges. The two are related by V being roughly the
 * integral of this, and confusing them is the easiest mistake to make here.
 *
 * @param t0,t1 the timespan; an infinite one has no transient to report
 * @param sn the refreshed network struct, carrying the reward definitions
 * @param opt CTMC options (state-space cutoff, tolerances, method)
 * @param tout optional out-parameter receiving the integration time points
 * @param names optional out-parameter receiving the reward names, in result order
 * @return `out[r][i]` is reward r at time `t[i]`
 */
template <class T>
std::vector<std::vector<T>> solver_ctmc_tran_reward(const NetworkStruct<T>& sn,
                                                    const CtmcOptions& opt, const T& t0,
                                                    const T& t1, std::vector<T>* tout = nullptr,
                                                    std::vector<std::string>* names = nullptr) {
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ctmc_tran_reward integrates the forward equation, which needs "
                  "transcendental arithmetic; use --arith double or real");
    if (sn.reward.empty())
        throw InputError(
            "solver_ctmc_tran_reward: no rewards are defined; declare one with "
            "set_reward(name, fn) before asking for a transient reward");
    if (!std::isfinite(num_traits<T>::to_double(t1)))
        throw InputError(
            "solver_ctmc_tran_reward: a finite timespan is required; an unbounded one has no "
            "transient to report, so ask for the steady-state reward instead");

    const CtmcTransient<T> tr = solver_ctmc_transient_analyzer(sn, opt, t0, t1);
    const Matrix<T> A = ctmc_state_space_aggr(sn, tr.chain.chain.space);
    const std::size_t n = tr.chain.chain.space.size(), nr = sn.reward.size(), nt = tr.t.size();
    const T zero = num_traits<T>::from_int(0);

    // The reward vector is state-dependent only, so it is evaluated ONCE per
    // state rather than once per (state, time): r does not depend on t, and
    // re-evaluating a user callback nt times would be the dominant cost.
    Matrix<T> R(nr, n, zero);
    if (names) names->clear();
    for (std::size_t r = 0; r < nr; ++r) {
        if (names) names->push_back(sn.reward[r].name);
        for (std::size_t s = 0; s < n; ++s) {
            std::vector<T> row(A.cols());
            for (std::size_t c = 0; c < row.size(); ++c) row[c] = A(s, c);
            R(r, s) = sn.reward[r].fn(row);
        }
    }

    std::vector<std::vector<T>> out(nr, std::vector<T>(nt, zero));
    for (std::size_t r = 0; r < nr; ++r)
        for (std::size_t i = 0; i < nt; ++i)
            for (std::size_t s = 0; s < n; ++s) out[r][i] += T(tr.pit(i, s) * R(r, s));
    if (tout) *tout = tr.t;
    return out;
}

/** Port of `@@SolverCTMC/getAvgReward`: the steady-state expected rewards. */
template <class T>
std::vector<T> solver_ctmc_avg_reward(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                      std::vector<std::string>* names = nullptr) {
    const CtmcReward<T> r = solver_ctmc_reward(sn, opt);
    if (names) *names = r.names;
    return r.steady_state;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_REWARD_H
