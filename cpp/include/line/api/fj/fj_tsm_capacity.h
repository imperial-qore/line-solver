/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_TSM_CAPACITY_H
#define LINE_API_FJ_TSM_CAPACITY_H

/**
 * Saturation throughput of the team service model.
 *
 * Templated port of matlab/src/api/fj/fj_tsm_capacity.m.
 *
 * A class-k job seizes r(k) of the s servers at once, holds them for a mean
 * x(k), and releases them all together. The apparent saturation rate is
 *
 *   Lambda_max = s / sum_k f(k) r(k) x(k),
 *
 * attainable only when the scheduler can pack jobs into execution states that
 * leave no server idle. The attainable capacity is the largest arrival rate for
 * which some mixture p over the feasible execution states balances every class,
 *
 *   maximise Lambda  s.t.  sum_j p_j n(j,k)/x(k) = Lambda f(k),
 *                          sum_j p_j = 1,  p >= 0,
 *
 * over the multisets of jobs whose total server demand is at most s. The
 * optimum equals Lambda_max exactly when every state carrying positive
 * probability is full capacity.
 *
 * For the two-server two-class case with r = (1,2), strict first come first
 * served cannot pack at all and reaches only
 *
 *   lambda_FCFS = 2 mu1 mu2 / (f1^2 mu2 + 2 f2^2 mu1 + 2 f1 f2 (mu1+mu2)).
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/simplex.h"

namespace line {
namespace fj {

/** [Lmax, Llp, Lfcfs, states, prob] of fj_tsm_capacity. */
template <class T>
struct FJTsmCapacityResult {
    T Lmax;
    T Llp;
    T Lfcfs;
    bool fcfs_available;
    std::vector<std::vector<unsigned> > states;
    std::vector<T> prob;
};

namespace detail {

/** Every multiset of jobs whose total server demand is at most s, minus the empty one. */
inline void tsm_states(const std::vector<unsigned>& r, unsigned left, std::size_t k,
                       std::vector<unsigned>& stack,
                       std::vector<std::vector<unsigned> >& out) {
    if (k == r.size()) {
        for (std::size_t i = 0; i < stack.size(); ++i)
            if (stack[i] > 0) { out.push_back(stack); return; }
        return;
    }
    const unsigned nmax = left / r[k];
    for (unsigned n = 0; n <= nmax; ++n) {
        stack[k] = n;
        tsm_states(r, left - n * r[k], k + 1, stack, out);
    }
    stack[k] = 0;
}

}  // namespace detail

/**
 * @param s number of servers
 * @param f class frequencies in the arrival stream, summing to one
 * @param r per-class server requirements, integers in 1..s
 * @param x per-class mean service times, all positive
 * @return  the apparent and attainable capacities with the optimal state mixture
 */
template <class T>
FJTsmCapacityResult<T> fj_tsm_capacity(unsigned s, const std::vector<T>& f,
                                       const std::vector<unsigned>& r, const std::vector<T>& x) {
    const std::size_t K = f.size();
    if (s < 1) throw InputError("fj_tsm_capacity: s must be a positive integer");
    if (r.size() != K || x.size() != K)
        throw InputError("fj_tsm_capacity: f, r and x must have the same length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T tot = zero;
    for (std::size_t k = 0; k < K; ++k) {
        if (f[k] < zero) throw InputError("fj_tsm_capacity: the class frequencies must be non-negative");
        if (r[k] < 1 || r[k] > s)
            throw InputError("fj_tsm_capacity: the server requirements must lie in 1..s");
        if (!(x[k] > zero)) throw InputError("fj_tsm_capacity: the mean service times must be positive");
        tot += f[k];
    }
    if (!(tot - one < num_traits<T>::from_double(1e-9)) ||
        !(one - tot < num_traits<T>::from_double(1e-9)))
        throw InputError("fj_tsm_capacity: the class frequencies must sum to one");

    FJTsmCapacityResult<T> out;
    T den = zero;
    for (std::size_t k = 0; k < K; ++k)
        den += f[k] * num_traits<T>::from_int(static_cast<long>(r[k])) * x[k];
    out.Lmax = num_traits<T>::from_int(static_cast<long>(s)) / den;

    std::vector<unsigned> stack(K, 0);
    detail::tsm_states(r, s, 0, stack, out.states);
    const std::size_t ns = out.states.size();
    if (ns == 0) throw NumericError("fj_tsm_capacity: no feasible execution state");

    // Variables [p_1..p_ns, Lambda], maximise Lambda
    lp::LpModel<T> lpm(ns + 1);
    lpm.set_maximize(true);
    lpm.set_cost(ns, one);
    for (std::size_t k = 0; k < K; ++k) {
        if (!(f[k] > zero)) continue;
        lpm.row_clear();
        for (std::size_t j = 0; j < ns; ++j)
            lpm.row_add(j, num_traits<T>::from_int(static_cast<long>(out.states[j][k])) / x[k]);
        lpm.row_add(ns, -f[k]);
        lpm.emit(lp::LpSense::EQ, zero);
    }
    lpm.row_clear();
    for (std::size_t j = 0; j < ns; ++j) lpm.row_add(j, one);
    lpm.emit(lp::LpSense::EQ, one);

    const lp::LpSolution<T> sol = lp::simplex_solve(lpm);
    if (!sol.ok())
        throw NumericError(std::string("fj_tsm_capacity: the capacity linear program returned ") +
                           lp::lp_status_name(sol.status));
    out.Llp = sol.x[ns];
    out.prob.assign(sol.x.begin(), sol.x.begin() + static_cast<long>(ns));

    // Strict first come first served capacity of the two-server two-class case
    out.fcfs_available = false;
    out.Lfcfs = zero;
    if (s == 2 && K == 2 && ((r[0] == 1 && r[1] == 2) || (r[0] == 2 && r[1] == 1))) {
        const std::size_t a = (r[0] == 1) ? 0 : 1, b = 1 - a;
        const T f1 = f[a], f2 = f[b], mu1 = one / x[a], mu2 = one / x[b];
        const T two = num_traits<T>::from_int(2);
        out.Lfcfs = two * mu1 * mu2 /
                    (f1 * f1 * mu2 + two * f2 * f2 * mu1 + two * f1 * f2 * (mu1 + mu2));
        out.fcfs_available = true;
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_TSM_CAPACITY_H
