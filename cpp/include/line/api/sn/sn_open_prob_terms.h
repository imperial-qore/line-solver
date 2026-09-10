/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_OPEN_PROB_TERMS_H
#define LINE_API_SN_SN_OPEN_PROB_TERMS_H

/**
 * Open-class contribution to an aggregate state probability at one station.
 *
 * Templated port of jar/src/main/java/jline/api/sn/SnOpenProbTerms.java, which
 * factors the mixed branch shared by `@@SolverMVA/getProbAggr.m` and
 * `@@SolverMVA/getProbSysAggr.m` (and the native Python `_open_prob_aggr_terms`)
 * out of both getters. The three station shapes are:
 *
 *  - EXT (a Source): no contribution, its population belongs to the environment;
 *  - INF (a Delay):  an independent Poisson per open class, mean Q(i,r);
 *  - anything else:  the multinomial-geometric BCMP form
 *                    (1 - sum_r rho_r) (sum_r n_r)! prod_r rho_r^n_r / n_r!.
 *
 * The whole term is returned in logs, since the callers accumulate a log
 * probability and only exponentiate at the end.
 *
 * INFEASIBILITY IS A FLAG, NOT -Inf. The reference returns negative infinity for
 * a state the law gives zero mass to (a class present where its utilization is
 * zero, or a saturated station). Rational has no infinity, so the verdict rides
 * in `feasible` and the caller decides what sentinel to emit; every caller in
 * this port short-circuits, so the log value is never read when it is false.
 *
 * ARITHMETIC: transcendental. Needs log and lgamma, so it does not instantiate
 * under Rational.
 */

#include <cmath>
#include <cstddef>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/** The log contribution, and whether the state carries any mass at all. */
template <class T>
struct OpenProbTerm {
    T logp;         ///< log of the open-class factor at this station
    bool feasible;  ///< false when the law gives the state zero probability
};

/**
 * @param sn  the network struct
 * @param Q   mean queue lengths, stations by classes
 * @param U   utilizations, stations by classes
 * @param nir per-class job counts, stations by classes (row `ist` is read)
 * @param ist 0-based station index
 */
template <class T>
OpenProbTerm<T> sn_open_prob_terms(const qn::NetworkStruct<T>& sn, const Matrix<T>& Q,
                                   const Matrix<T>& U, const Matrix<T>& nir, std::size_t ist) {
    static_assert(num_traits<T>::has_transcendental,
                  "sn_open_prob_terms evaluates a product-form law and needs logarithms");
    using std::log;
    if (ist >= sn.stations.size())
        throw InputError("sn_open_prob_terms: station index out of range");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const qn::SchedStrategy sched = sn.stations[ist].sched;

    OpenProbTerm<T> out;
    out.logp = zero;
    out.feasible = true;
    if (sched == qn::SchedStrategy::EXT) return out;

    if (sched == qn::SchedStrategy::INF) {
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (std::isfinite(sn.classes[r].population)) continue;
            const T q = Q(ist, r);
            if (q > zero) {
                out.logp = T(out.logp + nir(ist, r) * log(q) - q -
                             pfqn::detail::num_factln<T>(nir(ist, r)));
            } else if (nir(ist, r) > zero) {
                out.feasible = false;
                return out;
            }
        }
        return out;
    }

    T rho_total = zero, n_total = zero;
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        if (std::isfinite(sn.classes[r].population)) continue;
        rho_total = T(rho_total + U(ist, r));
        n_total = T(n_total + nir(ist, r));
    }
    if (!(rho_total < one)) {  // a saturated station has no stationary law
        out.feasible = false;
        return out;
    }
    out.logp = T(out.logp + log(T(one - rho_total)) + pfqn::detail::num_factln<T>(n_total));
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        if (std::isfinite(sn.classes[r].population)) continue;
        if (!(nir(ist, r) > zero)) continue;
        const T rho_r = U(ist, r);
        if (!(rho_r > zero)) {
            out.feasible = false;
            return out;
        }
        out.logp =
            T(out.logp + nir(ist, r) * log(rho_r) - pfqn::detail::num_factln<T>(nir(ist, r)));
    }
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_OPEN_PROB_TERMS_H
