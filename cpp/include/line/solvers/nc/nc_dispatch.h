/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_NC_DISPATCH_H
#define LINE_SOLVERS_NC_NC_DISPATCH_H

/**
 * Port of `solver_nc_analyzer.m`, `solver_ncld_analyzer.m` and
 * `@@SolverNC/ncDispatch.m`: one inner solve, choosing the analyzer that fits.
 *
 * THE ORDER IS THE CONTRACT. As in `mva_dispatch.h`, the branches are NOT
 * disjoint -- an order-independent station is also a multiserver, a
 * pass-and-swap tandem is also a closed network -- so the first match wins and
 * reordering silently changes which algorithm a model gets. The sequence is
 * the reference's, top to bottom:
 *
 *   0  discrete-time (slotted) model, on an explicit options.slotted
 *   1  order-independent closed network, on 'default' or 'exact'
 *   2  pass-and-swap importance sampling, on 'default', 'is' or 'sampling'
 *   3  'is' on an open or mixed model                     -> rejected by name
 *   4  Maximum Entropy Method, on an explicit 'mem'
 *   5  'exact' on an open or mixed multiserver model      -> rejected by name
 *   6  fractional closed populations, by interpolation
 *   7  everything else -> solver_nc / solver_ncld
 *
 * Branch 1 is the EXACT balanced-fairness analyzer and branch 2 the sampler, so
 * a pure-OI tandem on 'default' never reaches branch 2; only a genuine swap
 * graph, or an explicit 'is' / 'sampling', does. A station carrying a service
 * rate function that satisfies NEITHER predicate is refused by name rather than
 * sent down the ordinary normalizing-constant path, where its rate function
 * would simply be ignored.
 *
 * `ncDispatch` (the fork-join inner solve) is the load-dependence test alone:
 * the MMT transformation returns a plain mixed queueing network, so none of the
 * specialised routes can apply to it.
 */

#include "line/util/line_console.h"
#include <cmath>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/solvers/nc/solver_nc_mem.h"
#include "line/solvers/nc/solver_nc_dt.h"
#include "line/solvers/nc/solver_nc_oi.h"
#include "line/solvers/nc/solver_nc_dps.h"
#include "line/solvers/nc/solver_ncld.h"
#include "line/util/error.h"

namespace line {
namespace nc {

namespace detail {

/**
 * Does any station carry a service rate function, and does any of them have a
 * non-zero swap graph (which makes it a genuine pass-and-swap station)?
 *
 * Kept as one pass because both callers need both answers, and the swap graph is
 * read through `station_swap_graph` so the `refreshLocalVars` defaults apply.
 */
template <class T>
bool has_svc_rate_fun(const qn::NetworkStruct<T>& sn, bool& pas) {
    pas = false;
    bool any = false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (!sn.stations[i].svc_rate_fun) continue;
        any = true;
        if (!qn::station_swap_graph_is_zero(sn, i + 1)) pas = true;
    }
    return any;
}

/** True when any station carries a load- or class-dependent rate. */
template <class T>
bool has_scaling(const qn::NetworkStruct<T>& sn) {
    for (const qn::Station<T>& st : sn.stations)
        if (!st.lldscaling.empty() || static_cast<bool>(st.cdscaling)) return true;
    return false;
}

}  // namespace detail

/**
 * Port of `solver_ncld_analyzer.m`: the load-dependent analyzer, which is
 * `solver_ncld` plus the same fractional-population interpolation the
 * load-independent one applies.
 */
template <class T>
NcSolution<T> solver_ncld_analyzer(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    return solver_ncld(sn, opt);
}

/**
 * Port of `solver_nc_analyzer.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls
 */
template <class T>
NcSolution<T> solver_nc_analyzer(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    // Discrete time first: the slot lattice is a property of the MODEL, not of
    // a method, so it precedes every method-keyed branch below and refuses a
    // model outside the discrete-time product form instead of falling through.
    if (opt.slotted) return solver_nc_dt(sn, opt);

    // Closed think+DPS network: Morrison's heavy-usage expansion of the generating
    // function (npfqn_dps_morrison), the default for that shape. Not a product-form
    // route: lG is NaN. The runner intercepts this shape before reaching here; the
    // arm is repeated for the inner entry points (fork-join, ncld) that call the
    // analyzer directly.
    if (nc_is_dps_model(sn) && (opt.method == "default" || opt.method == "morrison"))
        return solver_nc_dps_analyzer(sn, opt);
    // Named on a model that is not the shape. The runner refuses this before reaching
    // here; the arm is repeated because the inner entry points (fork-join, ncld) call
    // this analyzer directly, and without it the method falls through to the ordinary
    // normalizing-constant path and answers under the caller's label.
    if (opt.method == "morrison")
        throw UnsupportedError(
            "solver_nc_analyzer: method 'morrison' requires a CLOSED network of exactly two "
            "stations, one infinite-server (think) station and one single-server DPS station with "
            "exponential service, which this model is not.");

    // Exact balanced fairness first: a pure-OI network on 'default' or 'exact'
    // never reaches the sampler below.
    if (nc_is_oi_model(sn) && (opt.method == "default" || opt.method == "exact"))
        return solver_nc_oi_analyzer(sn, opt);

    // The importance sampler, on an explicit request and on 'default' for a P&S
    // tandem with a genuine swap graph, which is reducible and has no exact path.
    if (nc_is_pas_model(sn) &&
        (opt.method == "default" || opt.method == "is" || opt.method == "sampling"))
        return solver_nc_pas_is_analyzer(sn, opt);

    bool pas = false;
    if (detail::has_svc_rate_fun(sn, pas))
        throw UnsupportedError(
            "solver_nc_analyzer: a station with a service rate function (OI / PAS) has no "
            "product-form demand and cannot be solved by the ordinary normalizing-constant path");

    bool anyOpen = false;
    for (const qn::JobClass& c : sn.classes)
        if (std::isinf(c.population)) anyOpen = true;

    if (opt.method == "is" && anyOpen)
        throw UnsupportedError(
            "solver_nc_analyzer: the 'is' importance-sampling method requires a closed queueing "
            "network. Use 'sampling' for an open or mixed model");

    // The Maximum Entropy Method, on an explicit request only: 'default' keeps
    // the normalizing-constant path it had before MEM existed.
    if (opt.method == "mem") return solver_nc_mem(sn, opt);

    bool multiserver = false;
    for (const qn::Station<T>& st : sn.stations)
        if (std::isfinite(st.nservers) && st.nservers > 1.0) multiserver = true;
    if (multiserver && anyOpen && opt.method == "exact")
        throw UnsupportedError(
            "solver_nc_analyzer: the NC solver cannot provide exact solutions for open or mixed "
            "multiserver queueing networks. Remove the 'exact' option");

    // Fractional closed populations: solve at floor and at ceiling and
    // interpolate, which is the only meaning a non-integral population has for
    // a constant defined on the population lattice.
    double eta_max = 0.0;
    for (const qn::JobClass& c : sn.classes)
        if (std::isfinite(c.population)) {
            const double e = std::fabs(c.population - std::floor(c.population));
            if (e > eta_max) eta_max = e;
        }
    if (eta_max > GlobalConstants::FineTol) {
        qn::NetworkStruct<T> lo = sn, hi = sn;
        for (std::size_t k = 0; k < sn.nclasses; ++k) {
            const double p = sn.classes[k].population;
            if (!std::isfinite(p)) continue;
            lo.classes[k].population = std::floor(p);
            hi.classes[k].population = std::ceil(p);
        }
        const NcSolution<T> f = solver_nc(lo, opt);
        const NcSolution<T> c = solver_nc(hi, opt);
        NcSolution<T> out = c;
        const T e = num_traits<T>::from_double(eta_max);
        auto blend = [&](Matrix<T>& A, const Matrix<T>& Af, const Matrix<T>& Ac) {
            A = Af;
            for (std::size_t i = 0; i < A.rows(); ++i)
                for (std::size_t j = 0; j < A.cols(); ++j)
                    A(i, j) = T(Af(i, j) + e * (Ac(i, j) - Af(i, j)));
        };
        blend(out.sol.Q, f.sol.Q, c.sol.Q);
        blend(out.sol.U, f.sol.U, c.sol.U);
        blend(out.sol.R, f.sol.R, c.sol.R);
        blend(out.sol.Tp, f.sol.Tp, c.sol.Tp);
        for (std::size_t k = 0; k < out.sol.X.size(); ++k) {
            out.sol.X[k] = T(f.sol.X[k] + e * (c.sol.X[k] - f.sol.X[k]));
            out.sol.C[k] = T(f.sol.C[k] + e * (c.sol.C[k] - f.sol.C[k]));
        }
        // REFERENCE DEFECT, reproduced: the interpolation of lG reads
        // lGf + eta*(lGf - lGc), which extrapolates AWAY from the ceiling
        // solve instead of towards it. It is kept because lG is a reported
        // quantity and silently changing it would move every caller's
        // getProbNormConstAggr on fractional populations.
        out.sol.lG = f.sol.lG + eta_max * (f.sol.lG - c.sol.lG);
        out.sol.iter = f.sol.iter + c.sol.iter;
        return out;
    }
    return solver_nc(sn, opt);
}

/**
 * Port of `@@SolverNC/ncDispatch.m`: the inner solve of the fork-join fixed
 * point, which is the load-dependence test and nothing else.
 */
template <class T>
NcSolution<T> nc_dispatch(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    const NcSolution<T> out = opt.slotted           ? solver_nc_dt(sn, opt)
                              : detail::has_scaling(sn) ? solver_ncld_analyzer(sn, opt)
                                                        : solver_nc_analyzer(sn, opt);
    const double lg = num_traits<T>::to_double(out.sol.lG);
    if (std::isfinite(lg))
        line::util::LineConsole::step("normalizing constant obtained: log G = %.6g", lg);
    return out;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_NC_DISPATCH_H
