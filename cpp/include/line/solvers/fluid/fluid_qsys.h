/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_QSYS_H
#define LINE_SOLVERS_FLUID_FLUID_QSYS_H

/**
 * Port of matlab/src/solvers/FLD/solver_fluid_qsys_analyzer.m: the
 * single-station fluid limits.
 *
 * A Source -> Queue -> Sink model with one class, answered by a closed-form
 * fluid or Gaussian limit rather than by integrating the network drift.
 *
 * WHY THESE ARE FLUID METHODS AND NOT MVA ONES. Each depends on the service or
 * patience law BEYOND ITS MEAN -- the stationary point of the Liu-Whitt model is
 * where the patience ccdf crosses 1/rho, the Mt/G/inf mean is a convolution with
 * the service ccdf -- and each is the limit of a sequence of systems, not an
 * approximation to a fixed one. That is the fluid solver's contract.
 *
 * METHODS
 *   `ggisgi.fluid` stationary point of the G/GI/s+GI fluid model (Liu and Whitt,
 *                  Operations Research 60(5), 2012)
 *   `ggingi.tga`   truncated Gaussian approximation, the O(sqrt(n)) fluctuation
 *                  around that point (Liu, Whitt and Yu, NRL 63(3), 2016)
 *   `tvms`         the Gt/Mt/st+GI many-server fluid queue at CONSTANT staffing
 *                  (Liu and Whitt, INFORMS J. Computing 26(1), 2014)
 *   `mtginf`       the exact Mt/G/inf mean (Eick, Massey and Whitt, Management
 *                  Science 39(2), 1993)
 *   `mol`          the modified-offered-load approximation for a finite server
 *                  count (Massey and Whitt, Ann. Appl. Prob. 4(4), 1994)
 *
 * ARITHMETIC: transcendental. Every one of them integrates or bisects.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/qsys/qsys_ggingi_tga.h"
#include "line/api/qsys/qsys_ggisgi_fluid.h"
#include "line/api/qsys/qsys_gtmtst_fluid.h"
#include "line/api/qsys/qsys_mtginf.h"
#include "line/api/qsys/qsys_mtgs0_mol.h"
#include "line/api/sn/sn_arrival_rate_fun.h"
#include "line/api/sn/sn_patience_handles.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/fluid/solver_fluid.h"

namespace line {
namespace fluid {

namespace detail {

/**
 * The method names the single-station limits answer, AFTER `fluid_unqualify`.
 *
 * That strips a leading `fluid.`, so the aliases `fluid.ggisgi` and `fluid.tga`
 * arrive here as the bare `ggisgi` and `tga` while the primary names
 * `ggisgi.fluid` and `ggingi.tga` pass through untouched. All four are listed;
 * `fluid_qsys_canonical` maps them onto the two primary names, which is what
 * the analyzer switches on and what `FluidSolution::method` reports.
 */
inline bool fluid_qsys_handles(const std::string& m) {
    return m == "ggisgi.fluid" || m == "ggisgi" || m == "ggingi.tga" || m == "tga" ||
           m == "tvms" || m == "mtginf" || m == "mol";
}

/** The primary name of an alias. */
inline std::string fluid_qsys_canonical(const std::string& m) {
    if (m == "ggisgi") return "ggisgi.fluid";
    if (m == "tga") return "ggingi.tga";
    return m;
}

/**
 * The integration window of a time-varying single-station limit, asked as a
 * predicate rather than thrown.
 *
 * A horizon is a solver OPTION and not a model feature, so the feature registry
 * has no name for it and a report has to ask this predicate directly.
 * `fluid_qsys_horizon` asks the same one on the solve path, which is what keeps
 * the report and the run from disagreeing about whether a method can be asked
 * for.
 *
 * @param opt the fluid knobs, read for the end of the horizon
 * @return an empty string when the window is a finite non-empty interval
 */
inline std::string fluid_qsys_horizon_reason(const FluidOptions& opt) {
    const double t1 = opt.timespan_end;
    if (!std::isfinite(t1) || t1 <= 0.0)
        return "solver_fluid_qsys: a time-varying fluid method needs a finite horizon; set "
               "options.timespan_end";
    return std::string();
}

/**
 * The integration window, [0, timespan_end].
 *
 * The C++ FluidOptions carries only the END of the horizon, where MATLAB and
 * Python carry a pair; the start is 0 for every fluid route in this port, which
 * is also what `runAnalyzer.m` forces when the start is not finite. A
 * non-positive or infinite end has no trajectory to report, and rather than
 * substituting `fluid_default_horizon`'s 30/min_rate -- a stationary-model rule
 * of thumb that says nothing about the PERIOD of a time-varying arrival -- the
 * caller is asked for one.
 */
template <class T>
void fluid_qsys_horizon(const qn::NetworkStruct<T>&, const FluidOptions& opt, double& t0,
                        double& t1) {
    const std::string reason = fluid_qsys_horizon_reason(opt);
    if (!reason.empty()) throw UnsupportedError(reason);
    t0 = 0.0;
    t1 = opt.timespan_end;
}

}  // namespace detail

/**
 * The three single-station limits that report a TRAJECTORY rather than a
 * stationary point, and so need a finite horizon; `ggisgi` and `tga` are
 * stationary and are not among them.
 *
 * @param m the method name, already unqualified
 * @return true when the method integrates over a finite horizon
 */
inline bool fluid_is_time_varying_limit(const std::string& m) {
    return m == "tvms" || m == "mtginf" || m == "mol";
}

/**
 * The horizon rule the time-varying limits impose, as a public predicate a
 * REPORT can ask: empty when `method` may be asked for with these options.
 *
 * `solver_fluid_qsys` refuses through the same body, so a pair the report
 * offers is a pair the limit runs.
 *
 * @param method the method name, qualified or not
 * @param opt the fluid knobs, read for the end of the horizon
 * @return an empty string when the method may run, else the refusal
 */
inline std::string fluid_qsys_horizon_supports(const std::string& method,
                                               const FluidOptions& opt) {
    // The `fluid.` prefix is stripped here rather than through
    // `detail::fluid_unqualify`, which lives in fluid_runner.h: that header
    // includes this one, so the dependency only goes the one way.
    const std::string m = (method.size() > 6 && method.compare(0, 6, "fluid.") == 0)
                              ? method.substr(6)
                              : method;
    if (!fluid_is_time_varying_limit(m)) return std::string();
    const std::string reason = detail::fluid_qsys_horizon_reason(opt);
    if (reason.empty()) return std::string();
    return "the '" + method + "' method reports a trajectory. " + reason;
}

/**
 * Solve a single-station model with one of the closed-form fluid limits.
 *
 * The steady-state row of a TIME-VARYING model is the TIME AVERAGE over the
 * horizon, which is what a stationary reader of a periodic system measures; the
 * trajectory itself is available from `solver_fluid_qsys_transient`.
 */
template <class T>
FluidSolution solver_fluid_qsys(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                                std::vector<FluidTranPoint>* traj = nullptr) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_fluid_qsys: the single-station fluid limits integrate and bisect, so they need "
            "transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
    const std::size_t M = sn.nstations;
    const std::size_t K = sn.nclasses;
    FluidSolution out;
    out.QN = Matrix<double>(M, K, 0.0);
    out.UN = Matrix<double>(M, K, 0.0);
    out.RN = Matrix<double>(M, K, 0.0);
    out.TN = Matrix<double>(M, K, 0.0);
    out.CN.assign(K, 0.0);
    out.XN.assign(K, 0.0);
    out.iters = 1;

    std::size_t src = 0;
    std::size_t qi = 0;
    bool haveSrc = false;
    bool haveQ = false;
    for (std::size_t i = 0; i < sn.nof_nodes(); ++i) {
        if (sn.nodes[i].nodetype == qn::NodeType::Source) {
            src = sn.nodes[i].station - 1;
            haveSrc = true;
        } else if (sn.nodes[i].nodetype == qn::NodeType::Queue ||
                   sn.nodes[i].nodetype == qn::NodeType::Delay) {
            qi = sn.nodes[i].station - 1;
            haveQ = true;
        }
    }
    // THE SHAPE THESE LIMITS ARE STATED FOR, refused by name rather than
    // answered on a model they do not describe: one open class through one
    // queueing station. The MVA qsys analyzer is reached by a structural
    // dispatch that guarantees it; these methods are selected by NAME, so the
    // check has to live here.
    if (!haveSrc || !haveQ)
        throw UnsupportedError(
            "solver_fluid_qsys: the single-station fluid limits need a Source and a queueing "
            "station");
    if (K != 1 || sn.nclosedjobs() > 0)
        throw UnsupportedError("solver_fluid_qsys: the '" + opt.method +
                               "' method is a single-station limit: it needs one open class "
                               "through one Source and one queueing station");

    const std::size_t qstateful = sn.stateful_of_station(qi + 1);
    const T Vq = sn.visits[0](qstateful - 1, 0);
    const T lambda = T(sn.rates(src, 0) * Vq);
    const T mu = sn.rates(qi, 0);
    const double nserv = sn.stations[qi].nservers;
    const T scvS = sn.scv(qi, 0);
    const T ca = num_traits<T>::from_double(std::sqrt(num_traits<T>::to_double(sn.scv(src, 0))));
    const T cs = num_traits<T>::from_double(std::sqrt(num_traits<T>::to_double(scvS)));
    const api::PatienceHandles<T> h = api::sn_patience_handles(sn, qi, 0);

    // The service ccdf, needed by the two Mt/G methods: they are exact in the
    // service DISTRIBUTION, not in its mean, which is the whole point of the
    // Eick-Massey-Whitt lag.
    const lang::Distrib<T>& svc = sn.service[qi][0];
    std::function<T(const T&)> serviceCcdf;
    if (svc.has_map()) {
        const mam::Map<T> sm = lang::dist_to_map(svc);
        const T one = num_traits<T>::from_int(1);
        serviceCcdf = [sm, one](const T& x) {
            std::vector<T> pts(1, x);
            return T(one - mam::map_cdf(sm, pts)[0]);
        };
    } else {
        serviceCcdf = [mu](const T& x) { return qsys::detail::num_exp(T(-mu * x)); };
    }
    const T ES = T(num_traits<T>::from_int(1) / mu);
    const double ES2 =
        (1.0 + num_traits<T>::to_double(scvS)) * num_traits<T>::to_double(ES) *
        num_traits<T>::to_double(ES);

    const std::string m = detail::fluid_qsys_canonical(opt.method);
    const double VqD = num_traits<T>::to_double(Vq);
    const double lamD = num_traits<T>::to_double(lambda);
    const double muD = num_traits<T>::to_double(mu);

    // Little's law on the CARRIED rate, as every LINE solver reports a station
    // that loses work.
    auto stationary = [&](double Lsys, double Tq, double Uq) {
        const double R = Tq > 0 ? Lsys / Tq : 0.0;
        out.RN(qi, 0) = R;
        out.QN(qi, 0) = Lsys;
        out.UN(qi, 0) = Uq;
        out.TN(qi, 0) = Tq;
        out.TN(src, 0) = lamD / VqD;
        out.XN[0] = Tq;
        out.CN[0] = R * VqD;
    };

    auto transient = [&](const std::vector<T>& t, const std::vector<T>& Lt,
                         const std::vector<T>& Ut, const std::vector<T>& Tt,
                         const std::vector<T>& arrival) {
        const std::size_t n = t.size();
        std::vector<double> td(n), Ld(n), Ud(n), Td(n), Ad(n);
        for (std::size_t i = 0; i < n; ++i) {
            td[i] = num_traits<T>::to_double(t[i]);
            Ld[i] = num_traits<T>::to_double(Lt[i]);
            Ud[i] = num_traits<T>::to_double(Ut[i]);
            Td[i] = num_traits<T>::to_double(Tt[i]);
            Ad[i] = num_traits<T>::to_double(arrival[i]);
        }
        auto trapz = [&](const std::vector<double>& y) {
            double s = 0;
            for (std::size_t i = 1; i < n; ++i) s += 0.5 * (y[i] + y[i - 1]) * (td[i] - td[i - 1]);
            return s;
        };
        const double span = td[n - 1] - td[0];
        const double Lbar = span > 0 ? trapz(Ld) / span : Ld[0];
        const double Ubar = span > 0 ? trapz(Ud) / span : Ud[0];
        const double Tbar = span > 0 ? trapz(Td) / span : Td[0];
        const double Abar = span > 0 ? trapz(Ad) / span : Ad[0];
        out.QN(qi, 0) = Lbar;
        out.UN(qi, 0) = Ubar;
        out.TN(qi, 0) = Tbar;
        out.TN(src, 0) = Abar / VqD;
        out.RN(qi, 0) = Tbar > 0 ? Lbar / Tbar : 0.0;
        out.XN[0] = Tbar;
        out.CN[0] = out.RN(qi, 0) * VqD;
        if (traj != nullptr) {
            traj->clear();
            traj->reserve(n);
            for (std::size_t i = 0; i < n; ++i) {
                FluidTranPoint p;
                p.t = td[i];
                p.QN = Matrix<double>(M, K, 0.0);
                p.UN = Matrix<double>(M, K, 0.0);
                p.TN = Matrix<double>(M, K, 0.0);
                p.QN(qi, 0) = Ld[i];
                p.UN(qi, 0) = Ud[i];
                p.TN(qi, 0) = Td[i];
                p.TN(src, 0) = Ad[i] / VqD;
                traj->push_back(p);
            }
        }
    };

    auto requirePatience = [&]() {
        if (!h.present)
            throw UnsupportedError("solver_fluid_qsys: the '" + m +
                                   "' method needs a reneging patience law on the queue "
                                   "(Queue.setPatience)");
    };
    auto requireFiniteServers = [&]() {
        if (!std::isfinite(nserv) || nserv < 1)
            throw UnsupportedError("solver_fluid_qsys: the '" + m +
                                   "' method needs a finite number of servers");
    };
    auto linspaceT = [](double a, double b, std::size_t n) {
        std::vector<T> v(n);
        for (std::size_t i = 0; i < n; ++i)
            v[i] = num_traits<T>::from_double(a + (b - a) * static_cast<double>(i) /
                                                      static_cast<double>(n - 1));
        return v;
    };

    if (m == "ggisgi.fluid") {
        requirePatience();
        const qsys::QsysFluidAbandonResult<T> r = qsys::qsys_ggisgi_fluid<T>(
            lambda, mu, static_cast<unsigned>(std::llround(nserv)), h.ccdf);
        stationary(num_traits<T>::to_double(r.meanNumber), num_traits<T>::to_double(r.throughput),
                   num_traits<T>::to_double(r.utilization));
    } else if (m == "ggingi.tga") {
        requirePatience();
        requireFiniteServers();
        const qsys::QsysTgaResult<T> r = qsys::qsys_ggingi_tga<T>(
            lambda, mu, static_cast<unsigned>(std::llround(nserv)), ca, cs, h.ccdf, h.pdf,
            serviceCcdf);
        const double pa = num_traits<T>::to_double(r.probAbandon);
        stationary(num_traits<T>::to_double(r.meanNumber), lamD * (1.0 - pa),
                   std::min(num_traits<T>::to_double(r.meanNumberInService) / nserv, 1.0));
    } else if (m == "tvms") {
        requirePatience();
        requireFiniteServers();
        const api::ArrivalRateFun<T> rf = api::sn_arrival_rate_fun(sn, src, 0);
        double t0 = 0, t1 = 0;
        detail::fluid_qsys_horizon(sn, opt, t0, t1);
        // CONSTANT STAFFING. Nothing in a Network declares a time-varying server
        // count, so s(t) is the station's own s; the time variation the method
        // is for enters through lambda(t) alone. A staffing schedule would need
        // a model feature that does not exist, and inventing one here would make
        // the solver answer a model the user did not build.
        const T sT = num_traits<T>::from_double(nserv);
        qsys::TvFluidOptions<T> tvopt;
        tvopt.pdf = h.pdf;
        const qsys::QsysTvFluidResult<T> r = qsys::qsys_gtmtst_fluid<T>(
            rf.lambda, [sT](const T&) { return sT; }, [mu](const T&) { return mu; }, h.ccdf,
            num_traits<T>::from_double(t1 - t0), tvopt);
        std::vector<T> times = r.times;
        for (std::size_t i = 0; i < times.size(); ++i)
            times[i] = T(times[i] + num_traits<T>::from_double(t0));
        std::vector<T> served(r.B.size());
        for (std::size_t i = 0; i < r.B.size(); ++i) served[i] = T(mu * r.B[i]);
        transient(times, r.X, r.utilization, served, r.arrivalRate);
    } else if (m == "mtginf") {
        const api::ArrivalRateFun<T> rf = api::sn_arrival_rate_fun(sn, src, 0);
        double t0 = 0, t1 = 0;
        detail::fluid_qsys_horizon(sn, opt, t0, t1);
        const qsys::QsysMtginfResult<T> r = qsys::qsys_mtginf<T>(
            rf.lambda, serviceCcdf, ES, linspaceT(t0, t1, 200),
            -std::numeric_limits<double>::infinity(), ES2);
        // An infinite-server station serves everything that arrives, so the
        // throughput is the arrival rate and the busy-server count is what a
        // utilization column can carry.
        transient(r.times, r.meanNumber, r.meanNumber, r.arrivalRate, r.arrivalRate);
    } else if (m == "mol") {
        requireFiniteServers();
        const api::ArrivalRateFun<T> rf = api::sn_arrival_rate_fun(sn, src, 0);
        double t0 = 0, t1 = 0;
        detail::fluid_qsys_horizon(sn, opt, t0, t1);
        const double cap = sn.cap[qi];
        // A finite buffer beyond the servers is not part of the loss model the
        // approximation is for; only s servers and no waiting room is.
        const bool useDelay = std::isfinite(cap) && cap > nserv;
        const qsys::QsysMolResult<T> r = qsys::qsys_mtgs0_mol<T>(
            rf.lambda, serviceCcdf, ES, static_cast<unsigned>(std::llround(nserv)),
            linspaceT(t0, t1, 200), -std::numeric_limits<double>::infinity(), useDelay);
        std::vector<T> util(r.meanBusyMOL.size());
        std::vector<T> served(r.meanBusyMOL.size());
        const T sT = num_traits<T>::from_double(nserv);
        for (std::size_t i = 0; i < r.meanBusyMOL.size(); ++i) {
            util[i] = T(r.meanBusyMOL[i] / sT);
            served[i] = T(mu * r.meanBusyMOL[i]);
        }
        transient(r.times, r.meanBusyMOL, util, served, r.arrivalRate);
    } else {
        throw UnsupportedError("solver_fluid_qsys: the '" + m +
                               "' method is not a single-station fluid limit");
    }
    out.method = m;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_QSYS_H
