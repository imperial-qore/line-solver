/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_PROB_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_PROB_H

/**
 * The state-probability half of the SolverMVA class surface.
 *
 * MVA computes means, not distributions, so every member of this family is an
 * approximation FITTED to the means the analyzer returned, and the reference
 * says which one:
 *
 *   closed classes  a binomial with mean Q(i,r) over N(r) trials, from
 *                   R. Schmidt, "An approximate MVA algorithm for exponential,
 *                   class-dependent multiple servers", PEVA 29:245-254, 1997
 *   open  classes   the exact BCMP product form of the station: independent
 *                   Poisson at an infinite server, multinomial-geometric at a
 *                   queue
 *
 * Ported here: `getProbMarg` (a single class's queue-length distribution),
 * `getProbNormConstAggr` (log G, by re-entering the analyzer at method='exact'),
 * and `getProbAggr` / `getProbSysAggr` (the per-class joint at one station, and
 * the whole-system joint). The last two read the per-class occupancy `nir` of
 * the model's state through `State.toMarginal`; this port carries the DEFAULT
 * initial state -- every closed class's population sits at its reference station,
 * open classes hold no jobs -- for which `toMarginal` reduces to that same
 * reference-station allocation, so `nir` is computed directly. A custom initial
 * state (there is no `setState` in the builder) would need the full state
 * encoding; when one exists it must be refused rather than silently reported
 * against the default.
 *
 * Everything is computed in logs and exponentiated once, as the reference does,
 * so a large population does not underflow before it is normalized.
 */

#include <cmath>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/sn/sn_open_prob_terms.h"
#include "line/api/sn/sn_state.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace mva {

/** A marginal distribution and its logarithm, over the states asked for. */
template <class T>
struct MargResult {
    std::vector<T> P;
    std::vector<T> logP;
};

namespace detail {

/** MATLAB `nchoosekln(n,k)` = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1). */
template <class T>
T num_nchoosekln(const T& n, const T& k) {
    return T(pfqn::detail::num_factln<T>(n) - pfqn::detail::num_factln<T>(k) -
             pfqn::detail::num_factln<T>(T(n - k)));
}

/**
 * The Schmidt binomial term for one closed class: log C(N,n) + n log(Q/N) +
 * (N-n) log(1 - Q/N).
 *
 * Written exactly as the reference writes it, INCLUDING the behaviour at the
 * endpoints: Q = 0 makes the second term -inf whenever n > 0 and 0 * log(0) = 0
 * (NaN in IEEE) when n = 0, and the reference takes `real(exp(.))` of whatever
 * comes out. Guarding the n = 0 case is not a defensive workaround but the
 * documented convention that an empty class contributes nothing.
 */
template <class T>
T binom_logterm(const T& N, const T& n, const T& Q) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    // C(0,0) = 1 and there is nothing to weight: a class with no population
    // contributes zero in logs. Forming Q/0 first would give p = +Inf, and the
    // guards below only happen to keep it out of the result because n and N-n
    // are both zero there -- which is luck, not a rule the next edit will keep.
    if (N == num_traits<T>::from_int(0)) return num_traits<T>::from_int(0);
    T lp = num_nchoosekln<T>(N, n);
    const T p = T(Q / N);
    if (n > zero) lp += T(n * log(p));
    if (T(N - n) > zero) lp += T(T(N - n) * log(T(one - p)));
    return lp;
}

}  // namespace detail

/**
 * Port of `@@SolverMVA/getProbMarg.m`: P(n jobs of class `r` at station `i`) for
 * the states in `states` (or the reference's own default range when empty).
 *
 * The three cases are the reference's, keyed on the CLASS being open or closed
 * and, when open, on the station's discipline:
 *
 *   closed          binomial(N_r, Q(i,r)/N_r) over n = 0..N_r; `states` selects
 *                   from that vector and a state above N_r is an error
 *   open  at INF    Poisson with mean Q(i,r)
 *   open  at a queue  P(n) = (1-rho) rho_r^n / (1-rho+rho_r)^(n+1), the exact
 *                   multiclass open product-form marginal
 *   open  at EXT    a Source has no queue-length distribution; P = 1 at n = 0
 *
 * When `states` is empty the open cases pick their own range, as the reference
 * does: mean + 5 sigma for the Poisson, and for the geometric the n at which the
 * tail falls below 1e-10, capped at 1000.
 */
template <class T>
MargResult<T> solver_mva_get_prob_marg(const qn::NetworkStruct<T>& L, const AvgResult<T>& avg,
                                       std::size_t ist, std::size_t r,
                                       const std::vector<long>& states,
                                       const std::string& method = "default") {
    static_assert(num_traits<T>::has_transcendental,
                  "getProbMarg fits a binomial / Poisson / geometric law and needs logarithms");
    using std::ceil;
    using std::exp;
    using std::log;
    using std::sqrt;
    if (ist == 0 || ist > L.nstations)
        throw InputError("getProbMarg: station index exceeds the number of stations in the model");
    if (r == 0 || r > L.nclasses)
        throw InputError("getProbMarg: job class index exceeds the number of classes in the model");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const double Nr = L.classes[r - 1].population;
    MargResult<T> out;

    if (std::isfinite(Nr)) {
        // MATLAB getProbMarg exact branch requires all single-server stations
        // (the pfqn_mvaldmx result is discarded and the binomial below is
        // returned regardless, so only the guard is observable).
        if (method == "exact")
            for (std::size_t s = 0; s < L.nstations; ++s)
                if (L.stations[s].nservers != 1.0)
                    throw UnsupportedError(
                        "getProbMarg: exact marginalized probabilities require single-server "
                        "stations");
        // ---- closed class: the Schmidt binomial over 0..N_r ------------------
        const long n_max = static_cast<long>(Nr);
        std::vector<T> all(static_cast<std::size_t>(n_max) + 1, zero), alllog(all.size(), zero);
        const T N = num_traits<T>::from_double(Nr);
        for (long k = 0; k <= n_max; ++k) {
            const T lp =
                detail::binom_logterm<T>(N, num_traits<T>::from_int(k), avg.QN(ist - 1, r - 1));
            alllog[static_cast<std::size_t>(k)] = lp;
            all[static_cast<std::size_t>(k)] = exp(lp);
        }
        if (states.empty()) {
            out.P = all;
            out.logP = alllog;
            return out;
        }
        for (long s : states) {
            if (s < 0 || s > n_max)
                throw InputError(
                    "getProbMarg: the requested state exceeds the maximum population for this "
                    "class");
            out.P.push_back(all[static_cast<std::size_t>(s)]);
            out.logP.push_back(alllog[static_cast<std::size_t>(s)]);
        }
        return out;
    }

    // ---- open class -------------------------------------------------------
    const lang::SchedStrategy sched = L.stations[ist - 1].sched;
    const T neg_inf = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());

    if (sched == lang::SchedStrategy::EXT) {
        // a Source is not a queue; the reference returns the degenerate law
        out.P.push_back(one);
        out.logP.push_back(zero);
        return out;
    }

    std::vector<long> sm = states;
    if (sched == lang::SchedStrategy::INF) {
        const T lam = avg.QN(ist - 1, r - 1);
        const double lamd = num_traits<T>::to_double(lam);
        if (sm.empty()) {
            const long nm = std::max<long>(
                1, static_cast<long>(std::ceil(lamd + 5.0 * std::sqrt(std::max(lamd, 1.0)))));
            for (long n = 0; n <= nm; ++n) sm.push_back(n);
        }
        for (long n : sm) {
            T lp;
            if (lam > zero) {
                lp = T(num_traits<T>::from_int(n) * log(lam) - lam -
                       pfqn::detail::num_factln<T>(num_traits<T>::from_int(n)));
            } else {
                lp = (n == 0) ? zero : neg_inf;
            }
            out.logP.push_back(lp);
            out.P.push_back(n == 0 && !(lam > zero) ? one : exp(lp));
        }
        return out;
    }

    // a queueing station: the multiclass open product-form marginal. rho is
    // capped just below 1 exactly as the reference caps it, so a saturated
    // station yields the all-zero law rather than a negative logarithm.
    const T rho_r = avg.UN(ist - 1, r - 1);
    T rho_tot = zero;
    for (std::size_t k = 0; k < L.nclasses; ++k) {
        const T u = avg.UN(ist - 1, k);
        if (num_traits<T>::to_double(u) == num_traits<T>::to_double(u)) rho_tot += u;  // skip NaN
    }
    const T rho_cap = T(one - num_traits<T>::from_double(GlobalConstants::FineTol));
    if (rho_tot > rho_cap) rho_tot = rho_cap;

    if (sm.empty()) {
        long nm = 0;
        if (rho_r > zero && rho_tot < one) {
            const T denom = T(one - rho_tot + rho_r);
            const double ratio = num_traits<T>::to_double(T(rho_r / denom));
            if (ratio > 0.0 && ratio < 1.0)
                nm = std::max<long>(1, static_cast<long>(std::ceil(-std::log(1e-10) / -std::log(ratio))));
            nm = std::min<long>(nm, 1000);
        }
        for (long n = 0; n <= nm; ++n) sm.push_back(n);
    }
    if (!(rho_tot < rho_cap)) {
        out.P.assign(sm.size(), zero);
        out.logP.assign(sm.size(), neg_inf);
        return out;
    }
    const T denom = T(one - rho_tot + rho_r);
    for (long n : sm) {
        const T nn = num_traits<T>::from_int(n);
        T lp;
        if (rho_r > zero) {
            lp = T(log(T(one - rho_tot)) + nn * log(rho_r) - T(nn + one) * log(denom));
        } else if (n == 0) {
            lp = T(log(T(one - rho_tot)) - log(denom));
        } else {
            lp = neg_inf;
        }
        out.logP.push_back(lp);
        out.P.push_back(exp(lp));
    }
    return out;
}

/**
 * Port of `@@SolverMVA/getProbNormConstAggr.m`: log G.
 *
 * The reference re-enters the analyzer with `method='exact'` rather than reusing
 * whatever the last solve produced, because only the exact MVA recursion carries
 * a normalizing constant; an AMVA result has no G to report. A cached value from
 * a previous solve is returned unchanged, which is what `self.result.Prob` does.
 */
template <class T>
T solver_mva_get_prob_norm_const_aggr(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    MvaOptions o = opt;
    o.method = "exact";
    const DispatchResult<T> dr = mva_dispatch(L, o, Matrix<T>());
    return dr.sol.lG;
}

namespace detail {

/**
 * The per-class occupancy `nir` of the model's default initial state, station by
 * station: every closed class's whole population sits at its reference station,
 * open classes hold no jobs. This is what `State.toMarginal(sn.state{i})`
 * returns for that state, so `getProbAggr`/`getProbSysAggr` read it directly.
 */
template <class T>
Matrix<T> initial_marginal(const qn::NetworkStruct<T>& L) {
    // THE MODEL'S DECLARED STATE, not a rebuilt default marking. The reference
    // reads `State.toMarginal(sn, ist, state{isf})`, so a `setState` moves the
    // question these getters answer; rebuilding the default here answered about
    // the reference-station placement under the caller's name. `sn_declared_marginal`
    // falls back to that placement per station where nothing was declared.
    return api::sn_declared_marginal<T>(L);
}

}  // namespace detail

/**
 * Port of `@@SolverMVA/getProbAggr.m`: P(n1 jobs of class 1, n2 of class 2, ...)
 * at station `ist` for the model's state, a scalar in [0,1] with its log.
 *
 * Closed classes take the Schmidt binomial fitted to Q(i,r); open classes take
 * the BCMP product form of the station (independent Poisson at an INF server,
 * multinomial-geometric at a queue, nothing at an EXT source), exactly as the
 * reference splits them.
 */
template <class T>
struct AggrResult {
    T P;
    T logP;
};

template <class T>
AggrResult<T> solver_mva_get_prob_aggr(const qn::NetworkStruct<T>& L, const AvgResult<T>& avg,
                                       std::size_t ist, const std::string& method = "default") {
    static_assert(num_traits<T>::has_transcendental,
                  "getProbAggr fits a binomial / product-form law and needs logarithms");
    using std::exp;
    using std::log;
    if (ist == 0 || ist > L.nstations)
        throw InputError("getProbAggr: station number exceeds the number of stations in the model");
    const std::size_t i = ist - 1;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const Matrix<T> nir = detail::initial_marginal(L);

    bool all_closed = true;
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (!std::isfinite(L.classes[r].population)) all_closed = false;

    // MATLAB getProbAggr: a closed model with method='exact' errors -- the exact
    // marginal is not implemented; only the Schmidt binomial approximation is.
    if (all_closed && method == "exact")
        throw UnsupportedError(
            "getProbAggr: exact marginal state probabilities are not available yet in SolverMVA");

    T logP = zero;
    if (!all_closed) {
        // Open classes: BCMP product form of the station, shared with getProbSysAggr.
        const api::OpenProbTerm<T> term = api::sn_open_prob_terms(L, avg.QN, avg.UN, nir, i);
        if (!term.feasible) return {zero, num_traits<T>::from_double(-1e308)};
        logP = T(logP + term.logp);
    }
    // Closed classes: Schmidt binomial.
    for (std::size_t r = 0; r < L.nclasses; ++r) {
        if (!std::isfinite(L.classes[r].population)) continue;
        const T N = num_traits<T>::from_double(L.classes[r].population);
        logP = T(logP + detail::binom_logterm<T>(N, nir(i, r), avg.QN(i, r)));
    }
    return {T(exp(logP)), logP};
}

/**
 * Port of `@@SolverMVA/getProbSysAggr.m`: the joint probability of the model's
 * whole state across all stations, a scalar in [0,1] with its log. Closed
 * classes take the multinomial-binomial normalization sum(factln(N)) - the
 * per-station product; open classes take the same per-station BCMP form as
 * getProbAggr.
 */
template <class T>
AggrResult<T> solver_mva_get_prob_sys_aggr(const qn::NetworkStruct<T>& L,
                                           const AvgResult<T>& avg,
                                           const std::string& method = "default") {
    static_assert(num_traits<T>::has_transcendental,
                  "getProbSysAggr fits a binomial / product-form law and needs logarithms");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const Matrix<T> nir = detail::initial_marginal(L);

    bool all_closed = true;
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (!std::isfinite(L.classes[r].population)) all_closed = false;

    // MATLAB getProbSysAggr: a closed model with method='exact' errors -- only
    // the Schmidt product-form approximation is implemented.
    if (all_closed && method == "exact")
        throw UnsupportedError(
            "getProbSysAggr: exact joint state probabilities are not available yet in SolverMVA");

    T logP = zero;
    if (all_closed) {
        for (std::size_t r = 0; r < L.nclasses; ++r)
            logP = T(logP + pfqn::detail::num_factln<T>(num_traits<T>::from_double(
                                L.classes[r].population)));
        for (std::size_t i = 0; i < L.nstations; ++i)
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                // A CLASS WITH NO POPULATION CONTRIBUTES NOTHING. Under class
                // switching a class can carry N=0 and still show a positive
                // QN(i,r), because the jobs in it arrived by switching; then
                // log(QN/0) is +Inf while nir(i,r) is 0, and 0*Inf is NaN,
                // which propagates through exp() and makes the whole joint
                // probability nan. Its binomial factor is C(0,0)=1, which is
                // zero in logs -- exactly what skipping it records.
                if (L.classes[r].population == 0.0) continue;
                const T N = num_traits<T>::from_double(L.classes[r].population);
                logP = T(logP - pfqn::detail::num_factln<T>(nir(i, r)));
                if (avg.QN(i, r) > zero)
                    logP = T(logP + nir(i, r) * log(T(avg.QN(i, r) / N)));
            }
        return {T(exp(logP)), logP};
    }

    // Mixed / open: closed-class multinomial normalization, then per station.
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (std::isfinite(L.classes[r].population))
            logP = T(logP + pfqn::detail::num_factln<T>(num_traits<T>::from_double(
                                L.classes[r].population)));
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const api::OpenProbTerm<T> term = api::sn_open_prob_terms(L, avg.QN, avg.UN, nir, i);
        if (!term.feasible) return {zero, num_traits<T>::from_double(-1e308)};
        logP = T(logP + term.logp);
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            if (!std::isfinite(L.classes[r].population)) continue;
            // See the closed branch: N=0 must not reach log(QN/N).
            if (L.classes[r].population == 0.0) continue;
            const T N = num_traits<T>::from_double(L.classes[r].population);
            logP = T(logP - pfqn::detail::num_factln<T>(nir(i, r)));
            if (avg.QN(i, r) > zero) logP = T(logP + nir(i, r) * log(T(avg.QN(i, r) / N)));
        }
    }
    return {T(exp(logP)), logP};
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_PROB_H
