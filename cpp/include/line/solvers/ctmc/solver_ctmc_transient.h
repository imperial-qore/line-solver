/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `solver_ctmc_transient_analyzer.m`: the time-dependent counterpart of
 * `solver_ctmc_analyzer`, integrating dpi/dt = pi Q from a point mass on the
 * initial state instead of solving pi Q = 0.
 *
 * WHAT IS AND IS NOT A TIME AVERAGE. The reference deliberately reports the
 * INSTANTANEOUS occupancy pi(t), not its running mean -- the commented-out
 * `cumsum(...)/t` lines in the reference are the time-average it decided
 * against. Q(t) and U(t) are therefore the state of the chain at t, and their
 * limits as t grows are the stationary values, not their averages over [0,t].
 *
 * THE UTILIZATION SWITCH IS NOT THE STATIONARY ONE. In steady state the carried
 * rate is available and `T*E[S]/c` is used; here there is no stationary
 * throughput, so utilization is read off the occupancy directly as
 * min(n_k, c)/c, and the PS and DPS families take their capacity share. The
 * reference WARNS for every other discipline and returns the FCFS form as an
 * approximation, which this port reproduces rather than silently improving:
 * a caller comparing against MATLAB must get MATLAB's number.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_TRANSIENT_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_TRANSIENT_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_fau.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/tr/fj_tag_transform.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

namespace detail {

/**
 * Transient trajectory by FAST ADAPTIVE UNIFORMIZATION, selected with
 * `options.config.transient_method = "fau"` (see `mc::ctmc_fau`).
 *
 * MARCHED, NOT RESTARTED. pi(t_{k+1}) comes from pi(t_k) over the step rather
 * than from pi(0) over the whole horizon, which keeps the cost proportional to
 * the grid instead of quadratic in it. Every step removes a little mass and
 * none puts any back, so the per-step tolerance is `fau_epsilon` divided by the
 * number of steps and the total defect stays below it. Nothing is renormalized:
 * the point of the method is that its error is a measured quantity.
 */
template <class T>
mc::TransientResult<T> ctmc_fau_transient(const Matrix<T>& Q, const std::vector<T>& pi0,
                                          const T& t0, const T& t1, const CtmcOptions& opt,
                                          const std::vector<T>& grid) {
    const double t0d = num_traits<T>::to_double(t0);
    const double t1d = num_traits<T>::to_double(t1);
    if (!std::isfinite(t1d))
        throw InputError(
            "solver_ctmc_transient_analyzer: transient_method 'fau' needs a finite horizon");

    std::vector<T> ts;
    if (!grid.empty()) {
        ts = grid;
    } else if (opt.timestep > 0.0) {
        const std::size_t nstep =
            static_cast<std::size_t>(std::floor((t1d - t0d) / opt.timestep));
        for (std::size_t i = 0; i <= nstep; ++i)
            ts.push_back(num_traits<T>::from_double(t0d + i * opt.timestep));
        if (num_traits<T>::to_double(ts.back()) < t1d) ts.push_back(t1);
    } else {
        const std::size_t ngrid = (opt.fau_ngrid > 1) ? opt.fau_ngrid : 100;
        for (std::size_t i = 0; i < ngrid; ++i)
            ts.push_back(num_traits<T>::from_double(
                t0d + (t1d - t0d) * static_cast<double>(i) / static_cast<double>(ngrid - 1)));
    }
    const std::size_t nt = ts.size();
    const double eps = (opt.fau_epsilon > 0.0) ? opt.fau_epsilon : 1e-6;
    const double epsStep = eps / static_cast<double>((nt > 1) ? (nt - 1) : 1);

    mc::TransientResult<T> out;
    out.t = ts;
    out.pi = Matrix<T>(nt, pi0.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < pi0.size(); ++j) out.pi(0, j) = pi0[j];
    std::vector<T> cur = pi0;
    for (std::size_t k = 1; k < nt; ++k) {
        const T dt = ts[k] - ts[k - 1];
        const mc::FauResult<T> r = mc::ctmc_fau(cur, Q, dt, epsStep, opt.fau_delta, -1);
        cur = r.pit;
        for (std::size_t j = 0; j < cur.size(); ++j) out.pi(k, j) = cur[j];
    }
    return out;
}

}  // namespace detail

/** What one transient CTMC solve produces. */
template <class T>
struct CtmcTransient {
    std::vector<T> t;                       ///< the time grid the solver chose
    Matrix<T> pit;                          ///< (ntimes x nstates) occupancy
    std::vector<std::vector<std::vector<T>>> QNt, UNt, TNt;  ///< [station][class][time]
    CtmcSolution<T> chain;                  ///< the generator and its state space
};

/**
 * Port of `solver_ctmc_transient_analyzer.m`.
 *
 * @param t0,t1 the timespan; the initial state is the model's default one
 *
 * `options.config.rate_sched` -- the time-INHOMOGENEOUS generator, where a
 * per-(station, class) rate follows a schedule and Q(t) is rebuilt by probing
 * its linear dependence on `sn.rates` -- is not ported, and there is NO refusal
 * for it here because there is nothing to refuse: `CtmcOptions` carries no such
 * field, so the schedule cannot be requested at this entry at all. What this
 * function solves is always the constant-rate generator. Adding the option
 * means adding the refusal with it, or the schedule would be accepted and
 * silently ignored, which would answer a different question.
 * @param sn the refreshed network struct
 * @param opt CTMC options (state-space cutoff, tolerances, method)
 */
template <class T>
CtmcTransient<T> solver_ctmc_transient_analyzer(const NetworkStruct<T>& sn,
                                                const CtmcOptions& opt, const T& t0, const T& t1,
                                                const std::vector<T>& grid = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ctmc_transient_analyzer integrates the forward equation with an "
                  "adaptive Runge-Kutta step controller, whose error norm is transcendental; "
                  "use --arith double or real");
    check_method(opt.method);
    // The reference refuses this by name too (`@@SolverCTMC/runAnalyzer.m:113`).
    // A fork-join model is solved on the TAG-AUGMENTED copy, whose per-class
    // trajectories are indexed by the auxiliary classes; folding them back is a
    // steady-state aggregate (`sn_fj_foldback` recomputes response time by
    // Little's law), and Little's law does not hold pointwise in time.
    if (tr::has_fork_join(sn))
        throw UnsupportedError(
            "SolverCTMC: transient analysis of a fork-join model is not supported. The chain is "
            "the tag-augmented one, and folding the sibling classes back onto the original ones "
            "is a steady-state aggregate; use the stationary solve, or SolverLDES for transients");

    CtmcTransient<T> out;
    out.chain = solver_ctmc_analyzer(sn, opt);
    const std::size_t n = out.chain.chain.space.size();
    const std::size_t M = sn.stations.size(), K = sn.nclasses;

    // The model's initial DISTRIBUTION: a point mass on the default state where
    // no prior is declared, and the product of the declared per-node priors
    // where one is. Unlike the stationary solve, where the component alone
    // matters, the transient answer depends on WHERE the chain starts, so a
    // state that is not in the space is an error here.
    std::vector<T> pi0;
    if (!analyzer_detail::init_state_distribution(sn, out.chain.chain.space, pi0))
        throw InputError(
            "solver_ctmc_transient_analyzer: the initial state is not contained in the state "
            "space, so there is no distribution to start the integration from");

    mc::TransientResult<T> tr;
    if (opt.transient_method == "fau") {
        tr = detail::ctmc_fau_transient(out.chain.chain.Q, pi0, t0, t1, opt, grid);
    } else if (opt.transient_method == "ode") {
        tr = mc::ctmc_transient(out.chain.chain.Q, pi0, t0, t1);
        // `options.timestep`, applied where the reference applies it: on the way
        // out of the integrator, not inside it. An explicit `grid` is the same
        // resampling on points a uniform step cannot express, and it WINS: a caller
        // that names the abscissae is integrating something against this trajectory
        // and needs its own, which is what the environment coupling does.
        if (!grid.empty())
            tr = mc::ctmc_transient_on_grid(out.chain.chain.Q, tr, grid);
        else if (opt.timestep > 0.0)
            tr = mc::ctmc_transient_on_grid(out.chain.chain.Q, tr, t0, t1,
                                            num_traits<T>::from_double(opt.timestep));
    } else {
        throw InputError("solver_ctmc_transient_analyzer: unknown transient_method '" +
                         opt.transient_method + "'; use 'ode' or 'fau'");
    }
    out.t = tr.t;
    out.pit = tr.pi;
    const std::size_t nt = out.t.size();

    // Clamp the numerical dust the integrator leaves below the zero threshold,
    // as the reference does before reading any measure off pi(t).
    for (std::size_t i = 0; i < nt; ++i)
        for (std::size_t s = 0; s < n; ++s)
            if (num_traits<T>::to_double(out.pit(i, s)) < GlobalConstants::Zero)
                out.pit(i, s) = num_traits<T>::from_int(0);

    const Matrix<T> A = ctmc_state_space_aggr(sn, out.chain.chain.space);
    const T zero = num_traits<T>::from_int(0);
    out.QNt.assign(M, std::vector<std::vector<T>>(K, std::vector<T>(nt, zero)));
    out.UNt.assign(M, std::vector<std::vector<T>>(K, std::vector<T>(nt, zero)));
    out.TNt.assign(M, std::vector<std::vector<T>>(K, std::vector<T>(nt, zero)));

    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        if (isf == 0) continue;
        const bool is_source = sn.stations[ist - 1].nodetype == NodeType::Source;
        const double S = sn.stations[ist - 1].nservers;
        const SchedStrategy sched = sn.stations[ist - 1].sched;

        for (std::size_t k = 1; k <= K; ++k)
            for (std::size_t i = 0; i < nt; ++i) {
                T acc = zero;
                for (std::size_t s = 0; s < n; ++s)
                    acc += T(out.pit(i, s) * out.chain.chain.dep_rates[s][isf - 1][k - 1]);
                out.TNt[ist - 1][k - 1][i] = acc;
            }
        // A Source's marginal is an encoding sentinel, so its queue length and
        // utilization are reported as zero rather than read off it.
        if (is_source) continue;

        for (std::size_t k = 1; k <= K; ++k)
            for (std::size_t i = 0; i < nt; ++i) {
                T q = zero;
                for (std::size_t s = 0; s < n; ++s) q += T(out.pit(i, s) * A(s, (ist - 1) * K + k - 1));
                out.QNt[ist - 1][k - 1][i] = q;
            }

        if (sched == SchedStrategy::INF) {
            for (std::size_t k = 0; k < K; ++k) out.UNt[ist - 1][k] = out.QNt[ist - 1][k];
            continue;
        }
        if (sched == SchedStrategy::PS) {
            // The capacity share of a processor-sharing station: class k takes
            // n_k / sum_j n_j of the min(total, c) busy servers.
            for (std::size_t k = 1; k <= K; ++k)
                for (std::size_t i = 0; i < nt; ++i) {
                    T u = zero;
                    for (std::size_t s = 0; s < n; ++s) {
                        double tot = 0;
                        for (std::size_t j = 0; j < K; ++j)
                            tot += num_traits<T>::to_double(A(s, (ist - 1) * K + j));
                        if (tot <= 0) continue;
                        const double nk = num_traits<T>::to_double(A(s, (ist - 1) * K + k - 1));
                        u += T(out.pit(i, s) *
                               num_traits<T>::from_double(std::min(nk, S) * nk / tot / S));
                    }
                    out.UNt[ist - 1][k - 1][i] = u;
                }
            continue;
        }
        if (sched == SchedStrategy::DPS) {
            const std::vector<T>& w = sn.stations[ist - 1].schedparam;
            for (std::size_t k = 1; k <= K; ++k)
                for (std::size_t i = 0; i < nt; ++i) {
                    T u = zero;
                    for (std::size_t s = 0; s < n; ++s) {
                        double wtot = 0;
                        for (std::size_t j = 0; j < K; ++j)
                            wtot += num_traits<T>::to_double(w[j]) *
                                    num_traits<T>::to_double(A(s, (ist - 1) * K + j));
                        if (wtot <= 0) continue;
                        const double nk = num_traits<T>::to_double(A(s, (ist - 1) * K + k - 1));
                        u += T(out.pit(i, s) * num_traits<T>::from_double(
                                                   S * num_traits<T>::to_double(w[k - 1]) * nk /
                                                   wtot));
                    }
                    out.UNt[ist - 1][k - 1][i] = u;
                }
            continue;
        }
        // FCFS, HOL, SIRO, SEPT, LEPT, SJF -- and, as an APPROXIMATION the
        // reference warns about, every remaining discipline.
        for (std::size_t k = 1; k <= K; ++k) {
            const lang::Distrib<T>& d = sn.service[ist - 1][k - 1];
            if (d.disabled || d.D0.rows() == 0) continue;
            for (std::size_t i = 0; i < nt; ++i) {
                T u = zero;
                for (std::size_t s = 0; s < n; ++s) {
                    const double nk = num_traits<T>::to_double(A(s, (ist - 1) * K + k - 1));
                    u += T(out.pit(i, s) * num_traits<T>::from_double(std::min(nk, S) / S));
                }
                out.UNt[ist - 1][k - 1][i] = u;
            }
        }
    }
    return out;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_TRANSIENT_H
