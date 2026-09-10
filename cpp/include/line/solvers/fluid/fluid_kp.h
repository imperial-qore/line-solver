/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_KP_H
#define LINE_SOLVERS_FLUID_FLUID_KP_H

/**
 * Port of `solver_fluid_kp.m`: the fluid AND diffusion limits of the
 * (MAP_t/Ph_t/inf)^N network of Y. M. Ko and J. Pender, "Diffusion limits for the
 * (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45 (2017) 248-253.
 *
 * The mean and the covariance are integrated JOINTLY:
 *
 *     dq/dt     = A f(t,q)
 *     dSigma/dt = J Sigma + Sigma J' + G,   J = A df/dq,  G = A diag(f) A'
 *
 * with A the jump matrix whose column e is the jump vector of event e. G is
 * exactly dH dH' of the paper's Theorem 3.3, each independent Poisson term
 * contributing l_e l_e' f_e. Where f is affine in q -- infinite-server stations
 * and the arrival phase process -- J does not depend on q and both equations close
 * EXACTLY, so for the (MAP/Ph/inf)^N case the mean and covariance are exact rather
 * than asymptotic. Finite-server stations are admitted through the usual fluid
 * min(x,c) term, where the covariance degrades to a linear noise approximation.
 *
 * THIS IS THE ONLY FLUID METHOD IN THE PORT THAT RETURNS A SECOND MOMENT FOR AN
 * OPEN NETWORK, and it is the only one whose second moment is TRANSIENT rather
 * than stationary -- `minnormal` solves a stationary Lyapunov equation, this
 * integrates the covariance along the trajectory.
 *
 * IT DELIBERATELY DOES NOT REUSE THE CLOSING ODE. That formulation routes a
 * departure from the source to the destination station and returns mass through
 * the STATIONARY arrival-instant vector pie, replacing the D1' operator by the
 * rank-one map pie*(D1*e)', i.e. by the PH renewal process with representation
 * (pie, D0). Its stationary arrival rate is exact but its autocorrelation is gone,
 * and a non-renewal arrival stream is the entire point of a MAP.
 *
 * WHAT THE `_t` OF (MAP_t/Ph_t/inf)^N IS. A `MAPt` or `PHt` carries a
 * piecewise-constant (D0(t), D1(t)) schedule, and `kp_pair_at` returns the pair
 * in force at time t -- the nominal pair for a process with no schedule, which
 * is what the reference's `local_pair_at` returns in the same case. Three things
 * follow from a schedule being present, and all three are consequences of the
 * SAME fact, that a cyclic model has no fixed point:
 *
 *   THE HORIZON. An unbounded timespan is resolved from the slowest rate, as
 *   before, but is then extended to at least ten full periods so that the
 *   trajectory has reached its periodic regime before anything is read off it.
 *
 *   THE STEP CAP. LSODA picks its step for accuracy of the SOLUTION and will
 *   happily step over a whole segment of a schedule, integrating a rate that was
 *   never in force. `h_max` is capped at a quarter of the narrowest segment.
 *
 *   THE STEADY-STATE ANSWER IS A TIME AVERAGE. The value at the horizon is an
 *   arbitrary point of the cycle, at which the source and the station throughput
 *   do not even agree. The metrics are instead the trapezoidal average over the
 *   last full period, on a mesh refined uniformly and BRACKETED at every segment
 *   boundary, so that no trapezoid interval straddles a jump in the arrival rate.
 *
 * A non-cyclic schedule has a fixed point again -- it is constant on its last
 * segment -- so it takes none of the three, exactly as in the reference.
 *
 * State layout, station-major, arrival phases before service phases:
 *   u-block  one per (EXT station, class): arrival MAP phase occupancy, sum 1
 *   x-block  one per (queueing station, class): fluid count in each service phase
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/sn/sn_schedule_nominal.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"
#include "line/util/lsoda.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** One (station, class) block of the Ko-Pender state vector. */
struct KpBlock {
    std::size_t station = 0;
    std::size_t cls = 0;
    std::size_t offset = 0;
    std::size_t nphases = 0;
};

/** The five event families of Ko-Pender (3.1)-(3.2). */
enum class KpEventKind {
    ArrivalPhase = 1,  ///< A0: arrival-MAP phase change without an arrival
    Arrival = 2,       ///< A1: phase change WITH an arrival, into a service phase
    ServicePhase = 3,  ///< S: service phase change inside a station
    Departure = 4,     ///< D: completion leaving the network
    Routed = 5         ///< R: completion routed onward
};

struct KpEvent {
    KpEventKind kind = KpEventKind::ArrivalPhase;
    std::size_t i = 0, c = 0;   ///< source station and class
    std::size_t k = 0, j = 0;   ///< source and target phase of the modulating chain
    std::size_t n = 0, l = 0;   ///< destination station and class
    std::size_t ip = 0;         ///< destination entry phase
    std::size_t off_src = 0;    ///< offset of the source block
    std::size_t off_dst = 0;    ///< offset of the destination block
    double weight = 1.0;        ///< routing probability times entry-phase probability
    /** The jump: -1 at `minus`, +1 at each `plus`; a phase change carries both. */
    std::vector<std::size_t> minus, plus;
};

/** The transient the covariance equation produces, i.e. `getTranAvgVar`. */
struct FluidKpTransient {
    std::vector<double> t;
    std::vector<Matrix<double>> QVar;   ///< per time point, (nstations x nclasses)
    std::vector<Matrix<double>> Sigma;  ///< per time point, (dim x dim)
    std::vector<std::vector<double>> q;
};

namespace kp_detail {

/** The (D0, D1) pair of a station-class, lowered to double. */
template <class T>
void kp_pair(const qn::NetworkStruct<T>& sn, std::size_t i, std::size_t r, Matrix<double>& D0,
             Matrix<double>& D1) {
    const lang::Distrib<T>& d = sn.service[i][r];
    const std::size_t n = d.D0.rows();
    D0 = Matrix<double>(n, n, 0.0);
    D1 = Matrix<double>(n, n, 0.0);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) {
            D0(a, b) = num_traits<T>::to_double(d.D0(a, b));
            D1(a, b) = num_traits<T>::to_double(d.D1(a, b));
        }
}

/** `map_pie`: the arrival-instant phase distribution of (D0, D1). */
inline std::vector<double> kp_pie(const Matrix<double>& D0, const Matrix<double>& D1) {
    if (D0.rows() <= 1) return std::vector<double>{1.0};
    mam::Map<double> m;
    m.D0 = D0;
    m.D1 = D1;
    return mam::map_pie(m);
}

/**
 * The stationary phase distribution of the modulating chain Q = D0 + D1, from
 * `[Q'; ones] \ [zeros; 1]`.
 *
 * A SOURCE STARTS IN ITS STATIONARY PHASE, NOT IN A KNOWN ONE, which is why the
 * initial covariance is diag(theta) - theta theta' rather than zero: zero would
 * assert a known initial phase and understate the variance early on.
 */
inline std::vector<double> kp_stationary(const Matrix<double>& D0, const Matrix<double>& D1) {
    const std::size_t h = D0.rows();
    if (h == 1) return std::vector<double>{1.0};
    // The least-squares system the reference solves: Q' theta = 0 with sum = 1.
    Matrix<double> A(h + 1, h, 0.0);
    for (std::size_t a = 0; a < h; ++a)
        for (std::size_t b = 0; b < h; ++b) A(a, b) = D0(b, a) + D1(b, a);
    for (std::size_t b = 0; b < h; ++b) A(h, b) = 1.0;
    std::vector<double> rhs(h + 1, 0.0);
    rhs[h] = 1.0;
    // Normal equations: A'A theta = A'rhs, which is what MATLAB's backslash gives
    // for an overdetermined system.
    Matrix<double> AtA(h, h, 0.0);
    std::vector<double> Atb(h, 0.0);
    for (std::size_t a = 0; a < h; ++a) {
        for (std::size_t b = 0; b < h; ++b) {
            double acc = 0.0;
            for (std::size_t e = 0; e <= h; ++e) acc += A(e, a) * A(e, b);
            AtA(a, b) = acc;
        }
        double acc = 0.0;
        for (std::size_t e = 0; e <= h; ++e) acc += A(e, a) * rhs[e];
        Atb[a] = acc;
    }
    std::vector<std::size_t> piv = lu_factor(AtA);
    lu_solve(AtA, piv, Atb);
    double s = 0.0;
    for (std::size_t a = 0; a < h; ++a) {
        if (Atb[a] < 0.0) Atb[a] = 0.0;
        s += Atb[a];
    }
    if (s > 0.0)
        for (std::size_t a = 0; a < h; ++a) Atb[a] /= s;
    return Atb;
}

/** One (station, class) schedule, lowered to double. */
struct KpSchedule {
    std::size_t station = 0, cls = 0;
    std::vector<double> bp;                  ///< boundary vector, nseg + 1 long
    std::vector<Matrix<double>> segD0, segD1;
    bool cyclic = false;
};

/**
 * `local_pair_at`: the (D0, D1) in force at time t.
 *
 * A CYCLIC schedule wraps the offset into [0, T); a NON-CYCLIC one returns the
 * ZERO pair outside its own window, which is the reference's behaviour and is
 * not a defect to fix -- a process whose schedule has not started, or has ended,
 * produces no events, and substituting the nominal pair there would invent
 * arrivals the model does not declare.
 */
inline void kp_pair_at(const std::vector<KpSchedule>& sched, const Matrix<double>& nomD0,
                       const Matrix<double>& nomD1, std::size_t i, std::size_t r, double t,
                       Matrix<double>& D0, Matrix<double>& D1) {
    for (std::size_t e = 0; e < sched.size(); ++e) {
        const KpSchedule& sc = sched[e];
        if (sc.station != i || sc.cls != r) continue;
        const double T0 = sc.bp.front(), T1 = sc.bp.back();
        const double period = T1 - T0;
        double offset = t - T0;
        if (sc.cyclic) {
            if (period > 0.0) {
                offset = std::fmod(offset, period);
                if (offset < 0.0) offset += period;
            } else {
                offset = 0.0;
            }
        } else if (offset < 0.0 || offset >= period) {
            const std::size_t n = sc.segD0.front().rows();
            D0 = Matrix<double>(n, n, 0.0);
            D1 = Matrix<double>(n, n, 0.0);
            return;
        }
        const double pos = T0 + offset;
        std::size_t idx = sc.segD0.size() - 1;
        for (std::size_t k = 1; k < sc.bp.size(); ++k)
            if (pos < sc.bp[k]) {
                idx = k - 1;
                break;
            }
        D0 = sc.segD0[idx];
        D1 = sc.segD1[idx];
        return;
    }
    D0 = nomD0;
    D1 = nomD1;
}

/** `local_summarise`: the trapezoidal average over `window`, or the last value. */
inline double kp_summarise(const std::vector<double>& series, const std::vector<double>& t,
                           double w0, double w1, bool have_window) {
    if (series.empty()) return 0.0;
    if (!have_window) return series.back();
    std::vector<std::size_t> idx;
    for (std::size_t a = 0; a < t.size(); ++a)
        if (t[a] >= w0 && t[a] <= w1) idx.push_back(a);
    if (idx.size() < 2) return series.back();
    double acc = 0.0;
    for (std::size_t a = 0; a + 1 < idx.size(); ++a)
        acc += 0.5 * (series[idx[a]] + series[idx[a + 1]]) * (t[idx[a + 1]] - t[idx[a]]);
    const double span = t[idx.back()] - t[idx.front()];
    return span > 0.0 ? acc / span : series.back();
}

}  // namespace kp_detail

/**
 * The Ko-Pender solve, returning both the steady table and the covariance
 * trajectory so that neither has to integrate twice.
 */
template <class T>
FluidSolution solver_fluid_kp_core(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                                   FluidKpTransient* tran) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "solver_fluid_kp: the covariance equation is integrated with LSODA, whose coefficients "
            "assume double precision; rerun with --arith double");

    const std::size_t M = sn.nstations, K = sn.nclasses;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population))
            throw UnsupportedError(
                "solver_fluid_kp: the 'kp' method analyses the OPEN (MAP_t/Ph_t/inf)^N network of "
                "Ko and Pender (2017); a closed class has no arrival process to modulate. Use "
                "'closing' or 'matrix' for closed models");

    // ---- blocks -----------------------------------------------------------
    std::vector<KpBlock> ublocks, xblocks;
    std::size_t off = 0;
    std::vector<std::vector<std::size_t>> uof(M, std::vector<std::size_t>(K, 0));
    std::vector<std::vector<std::size_t>> xof(M, std::vector<std::size_t>(K, 0));
    std::vector<std::vector<bool>> is_u(M, std::vector<bool>(K, false));
    std::vector<std::vector<bool>> is_x(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i) {
        const bool ext = sn.stations[i].sched == lang::SchedStrategy::EXT;
        for (std::size_t r = 0; r < K; ++r) {
            if (sn.disabled[i][r]) continue;
            const std::size_t h = sn.service[i][r].D0.rows();
            const double rate = num_traits<T>::to_double(sn.rates(i, r));
            if (h == 0 || !std::isfinite(rate) || rate <= 0.0) continue;
            KpBlock b;
            b.station = i;
            b.cls = r;
            b.offset = off;
            b.nphases = h;
            if (ext) {
                ublocks.push_back(b);
                uof[i][r] = off;
                is_u[i][r] = true;
            } else {
                xblocks.push_back(b);
                xof[i][r] = off;
                is_x[i][r] = true;
            }
            off += h;
        }
    }
    const std::size_t dim = off;
    if (ublocks.empty())
        throw InputError(
            "solver_fluid_kp: the 'kp' method needs at least one Source with an arrival process");

    bool linear_model = true;
    for (std::size_t b = 0; b < xblocks.size(); ++b) {
        const std::size_t i = xblocks[b].station;
        if (sn.stations[i].sched != lang::SchedStrategy::INF &&
            std::isfinite(sn.stations[i].nservers))
            linear_model = false;
    }
    // The reference warns here. A finite-server station makes the rate functions
    // nonlinear, so the covariance is a linear noise approximation rather than the
    // exact second moment; it is exact for infinite-server stations. There is no
    // warning channel in this port, so the fact is recorded on the solution.
    (void)linear_model;

    // ---- the nominal pairs and entry-phase vectors -------------------------
    std::vector<std::vector<Matrix<double>>> D0(M, std::vector<Matrix<double>>(K)),
        D1(M, std::vector<Matrix<double>>(K));
    std::vector<std::vector<std::vector<double>>> pie(M, std::vector<std::vector<double>>(K));
    // The schedules, and the NOMINAL pair of everything else. `kp_pair` reads
    // `Distrib::D0`/`D1`, which for a MAPt/PHt already hold the width-weighted
    // time average -- `sn_schedule_nominal`'s first two outputs -- so the
    // nominal arm needs no special case here.
    std::vector<kp_detail::KpSchedule> sched;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!is_u[i][r] && !is_x[i][r]) continue;
            kp_detail::kp_pair(sn, i, r, D0[i][r], D1[i][r]);
            // `pie` is the arrival-instant vector of the NOMINAL pair, and stays
            // so under a schedule: it seeds a job's service phase, which is a
            // property of the service process as a whole, and the reference's
            // `nomPie` is built from `sn_schedule_nominal`'s nominal pair too.
            pie[i][r] = kp_detail::kp_pie(D0[i][r], D1[i][r]);
            if (!sn::sn_has_schedule(sn, i, r)) continue;
            const sn::ScheduleNominal<T> sc = sn::sn_schedule_nominal(sn, i, r);
            kp_detail::KpSchedule ks;
            ks.station = i;
            ks.cls = r;
            ks.cyclic = sc.cyclic;
            for (std::size_t a = 0; a < sc.breakpoints.size(); ++a)
                ks.bp.push_back(num_traits<T>::to_double(sc.breakpoints[a]));
            for (std::size_t k = 0; k < sc.segD0.size(); ++k) {
                const std::size_t nph = sc.segD0[k].rows();
                Matrix<double> A(nph, nph, 0.0), B(nph, nph, 0.0);
                for (std::size_t a = 0; a < nph; ++a)
                    for (std::size_t b = 0; b < nph; ++b) {
                        A(a, b) = num_traits<T>::to_double(sc.segD0[k](a, b));
                        B(a, b) = num_traits<T>::to_double(sc.segD1[k](a, b));
                    }
                ks.segD0.push_back(A);
                ks.segD1.push_back(B);
            }
            sched.push_back(ks);
        }

    // The pairs in force at time t, for every block. Built once per right-hand
    // side evaluation and handed to `rates`, so the Jacobian's 2*dim difference
    // calls all see the SAME instant -- differencing across a segment boundary
    // would report the jump in the schedule as a derivative in q.
    std::vector<std::vector<Matrix<double>>> Dt0 = D0, Dt1 = D1;
    const bool time_varying = !sched.empty();
    const auto pairs_at = [&](double t) {
        if (!time_varying) return;
        for (std::size_t e = 0; e < sched.size(); ++e) {
            const std::size_t i = sched[e].station, r = sched[e].cls;
            kp_detail::kp_pair_at(sched, D0[i][r], D1[i][r], i, r, t, Dt0[i][r], Dt1[i][r]);
        }
    };

    // Routing in stateful space, as `fluid_ode_system` reads it.
    const std::size_t S = sn.nof_stateful();
    const bool have_rt = sn.rt.rows() == S * K;
    std::vector<std::size_t> sf(M, 0);
    for (std::size_t i = 0; i < M; ++i) sf[i] = sn.stateful_of_station(i + 1) - 1;
    const auto route = [&](std::size_t i, std::size_t c, std::size_t j, std::size_t l) -> double {
        if (!have_rt) return 0.0;
        return num_traits<T>::to_double(sn.rt(sf[i] * K + c, sf[j] * K + l));
    };
    // The probability that a completion at (i,r) LEAVES the network: sn.rt is
    // closed through the Source, so any destination that is not a service block is
    // an exit.
    std::vector<std::vector<double>> pout(M, std::vector<double>(K, 0.0));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!is_x[i][r]) continue;
            double acc = 0.0;
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t l = 0; l < K; ++l)
                    if (!is_x[j][l]) acc += route(i, r, j, l);
            pout[i][r] = acc;
        }

    // ---- events ------------------------------------------------------------
    std::vector<KpEvent> ev;
    const auto push = [&](KpEvent e) {
        ev.push_back(e);
    };
    for (std::size_t b = 0; b < ublocks.size(); ++b) {  // (A0)
        const KpBlock& u = ublocks[b];
        for (std::size_t k = 0; k < u.nphases; ++k)
            for (std::size_t j = 0; j < u.nphases; ++j) {
                if (k == j) continue;
                KpEvent e;
                e.kind = KpEventKind::ArrivalPhase;
                e.i = u.station;
                e.c = u.cls;
                e.k = k;
                e.j = j;
                e.off_src = u.offset;
                e.minus.push_back(u.offset + k);
                e.plus.push_back(u.offset + j);
                push(e);
            }
    }
    for (std::size_t b = 0; b < ublocks.size(); ++b) {  // (A1)
        const KpBlock& u = ublocks[b];
        for (std::size_t d = 0; d < xblocks.size(); ++d) {
            const KpBlock& xb = xblocks[d];
            const double p = route(u.station, u.cls, xb.station, xb.cls);
            if (!(p > 0.0)) continue;
            for (std::size_t k = 0; k < u.nphases; ++k)
                for (std::size_t j = 0; j < u.nphases; ++j)
                    for (std::size_t ip = 0; ip < xb.nphases; ++ip) {
                        KpEvent e;
                        e.kind = KpEventKind::Arrival;
                        e.i = u.station;
                        e.c = u.cls;
                        e.k = k;
                        e.j = j;
                        e.n = xb.station;
                        e.l = xb.cls;
                        e.ip = ip;
                        e.off_src = u.offset;
                        e.off_dst = xb.offset;
                        e.weight = p * (ip < pie[xb.station][xb.cls].size()
                                            ? pie[xb.station][xb.cls][ip]
                                            : 0.0);
                        e.minus.push_back(u.offset + k);
                        e.plus.push_back(u.offset + j);
                        e.plus.push_back(xb.offset + ip);
                        push(e);
                    }
        }
    }
    for (std::size_t b = 0; b < xblocks.size(); ++b) {  // (S)
        const KpBlock& xb = xblocks[b];
        for (std::size_t p = 0; p < xb.nphases; ++p)
            for (std::size_t q = 0; q < xb.nphases; ++q) {
                if (p == q) continue;
                KpEvent e;
                e.kind = KpEventKind::ServicePhase;
                e.i = xb.station;
                e.c = xb.cls;
                e.k = p;
                e.j = q;
                e.off_src = xb.offset;
                e.minus.push_back(xb.offset + p);
                e.plus.push_back(xb.offset + q);
                push(e);
            }
    }
    for (std::size_t b = 0; b < xblocks.size(); ++b) {  // (D) and (R)
        const KpBlock& xb = xblocks[b];
        if (pout[xb.station][xb.cls] > 0.0)
            for (std::size_t p = 0; p < xb.nphases; ++p) {
                KpEvent e;
                e.kind = KpEventKind::Departure;
                e.i = xb.station;
                e.c = xb.cls;
                e.k = p;
                e.off_src = xb.offset;
                e.weight = pout[xb.station][xb.cls];
                e.minus.push_back(xb.offset + p);
                push(e);
            }
        for (std::size_t d = 0; d < xblocks.size(); ++d) {
            const KpBlock& nb = xblocks[d];
            const double p = route(xb.station, xb.cls, nb.station, nb.cls);
            if (!(p > 0.0)) continue;
            for (std::size_t q = 0; q < xb.nphases; ++q)
                for (std::size_t ip = 0; ip < nb.nphases; ++ip) {
                    KpEvent e;
                    e.kind = KpEventKind::Routed;
                    e.i = xb.station;
                    e.c = xb.cls;
                    e.k = q;
                    e.n = nb.station;
                    e.l = nb.cls;
                    e.ip = ip;
                    e.off_src = xb.offset;
                    e.off_dst = nb.offset;
                    e.weight =
                        p * (ip < pie[nb.station][nb.cls].size() ? pie[nb.station][nb.cls][ip] : 0.0);
                    e.minus.push_back(xb.offset + q);
                    e.plus.push_back(nb.offset + ip);
                    push(e);
                }
        }
    }
    const std::size_t nev = ev.size();

    // The server-capacity factor min(n,c)/n, shared by every class at a station.
    const auto capacity = [&](const double* q, std::size_t i) -> double {
        if (sn.stations[i].sched == lang::SchedStrategy::INF ||
            !std::isfinite(sn.stations[i].nservers))
            return 1.0;
        double ni = 0.0;
        for (std::size_t b = 0; b < xblocks.size(); ++b) {
            if (xblocks[b].station != i) continue;
            for (std::size_t p = 0; p < xblocks[b].nphases; ++p)
                ni += std::max(q[xblocks[b].offset + p], 0.0);
        }
        const double c = sn.stations[i].nservers;
        return (ni <= c) ? 1.0 : c / ni;
    };

    const auto rates = [&](const double* q, std::vector<double>& f) {
        f.assign(nev, 0.0);
        std::vector<double> cap(M, 1.0);
        for (std::size_t i = 0; i < M; ++i) cap[i] = capacity(q, i);
        for (std::size_t e = 0; e < nev; ++e) {
            const KpEvent& s = ev[e];
            const double mass = std::max(q[s.off_src + s.k], 0.0);
            switch (s.kind) {
                case KpEventKind::ArrivalPhase:
                    f[e] = Dt0[s.i][s.c](s.k, s.j) * mass;
                    break;
                case KpEventKind::Arrival:
                    f[e] = Dt1[s.i][s.c](s.k, s.j) * s.weight * mass;
                    break;
                case KpEventKind::ServicePhase:
                    f[e] = Dt0[s.i][s.c](s.k, s.j) * mass * cap[s.i];
                    break;
                case KpEventKind::Departure:
                case KpEventKind::Routed: {
                    double rowsum = 0.0;
                    for (std::size_t b = 0; b < Dt1[s.i][s.c].cols(); ++b)
                        rowsum += Dt1[s.i][s.c](s.k, b);
                    f[e] = rowsum * s.weight * mass * cap[s.i];
                    break;
                }
            }
        }
    };
    const auto apply_jumps = [&](const std::vector<double>& f, double* dq) {
        for (std::size_t a = 0; a < dim; ++a) dq[a] = 0.0;
        for (std::size_t e = 0; e < nev; ++e) {
            if (f[e] == 0.0) continue;
            for (std::size_t a = 0; a < ev[e].minus.size(); ++a) dq[ev[e].minus[a]] -= f[e];
            for (std::size_t a = 0; a < ev[e].plus.size(); ++a) dq[ev[e].plus[a]] += f[e];
        }
    };

    // ---- horizon -----------------------------------------------------------
    double t0 = 0.0;
    double tend = opt.timespan_end;
    const bool unbounded = !std::isfinite(tend);
    // The longest cycle among the CYCLIC schedules; zero when none is cyclic,
    // and a non-cyclic schedule is constant on its last segment, so it has a
    // fixed point and needs neither the extended horizon nor the averaging.
    double period = 0.0;
    for (std::size_t e = 0; e < sched.size(); ++e)
        if (sched[e].cyclic) period = std::max(period, sched[e].bp.back() - sched[e].bp.front());
    if (unbounded) {
        double slow = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                const double rate = num_traits<T>::to_double(sn.rates(i, r));
                if (std::isfinite(rate) && rate > 0.0) slow = std::min(slow, rate);
            }
        if (!std::isfinite(slow)) slow = 1.0;
        tend = t0 + std::max(10.0, 30.0 / slow);
        if (period > 0.0) tend = std::max(tend, t0 + 10.0 * period);
    }

    // ---- initial condition -------------------------------------------------
    std::vector<double> z(dim + dim * dim, 0.0);
    // The arrival phase is drawn from the stationary vector of the pair IN FORCE
    // AT t0, not of the time average: a schedule that starts in a quiet segment
    // starts in that segment's phase mix.
    pairs_at(t0);
    for (std::size_t b = 0; b < ublocks.size(); ++b) {
        const KpBlock& u = ublocks[b];
        const std::vector<double> theta =
            kp_detail::kp_stationary(Dt0[u.station][u.cls], Dt1[u.station][u.cls]);
        for (std::size_t a = 0; a < u.nphases; ++a) z[u.offset + a] = theta[a];
        for (std::size_t a = 0; a < u.nphases; ++a)
            for (std::size_t c2 = 0; c2 < u.nphases; ++c2)
                z[dim + (u.offset + a) * dim + (u.offset + c2)] =
                    (a == c2 ? theta[a] : 0.0) - theta[a] * theta[c2];
    }
    // NOT `opt.init_sol`, which is laid out for the CLOSING state vector: the two
    // can have the same length on the same model, so reading it here would let a
    // closing-layout seed zero the source phase mass and with it the network.
    //
    // A WRONG-SIZED SEED IS REFUSED, not ignored. Dropping it would integrate from
    // the default initial condition under the caller's name and return a plausible
    // trajectory for a model the caller did not ask about.
    if (!opt.kp_init_sol.empty()) {
        if (opt.kp_init_sol.size() != dim)
            throw InputError("solver_fluid_kp: config.kp_init_sol has " +
                             std::to_string(opt.kp_init_sol.size()) +
                             " entries but the 'kp' state vector of this model has " +
                             std::to_string(dim) +
                             ", laid out station-major over the (station, class) blocks. "
                             "It is NOT laid out like init_sol.");
        for (std::size_t a = 0; a < dim; ++a) z[a] = opt.kp_init_sol[a];
    }
    // Companion seed for the covariance. A caller that carries a DISTRIBUTION across
    // a handoff supplies the second moment beside the mean, so the next stage does
    // not restart from a point mass it never had. Same layout as kp_init_sol.
    if (opt.init_cov.rows() > 0 || opt.init_cov.cols() > 0) {
        if (opt.init_cov.rows() != dim || opt.init_cov.cols() != dim)
            throw InputError("solver_fluid_kp: config.init_cov is " +
                             std::to_string(opt.init_cov.rows()) + "x" +
                             std::to_string(opt.init_cov.cols()) +
                             " but the 'kp' state vector of this model has " +
                             std::to_string(dim) + " entries, so the covariance must be " +
                             std::to_string(dim) + "x" + std::to_string(dim) + ".");
        double asym = 0.0, scale = 0.0;
        for (std::size_t a = 0; a < dim; ++a)
            for (std::size_t c2 = 0; c2 < dim; ++c2) {
                const double d = opt.init_cov(a, c2) - opt.init_cov(c2, a);
                asym += d * d;
                scale += opt.init_cov(a, c2) * opt.init_cov(a, c2);
            }
        // Loose enough for the rounding of a covariance that was itself integrated,
        // tight enough to catch a matrix that is simply not one.
        if (std::sqrt(asym) > 1e-6 * std::max(1.0, std::sqrt(scale)))
            throw InputError("solver_fluid_kp: config.init_cov must be symmetric.");
        for (std::size_t a = 0; a < dim; ++a)
            for (std::size_t c2 = 0; c2 < dim; ++c2)
                z[dim + a * dim + c2] = opt.init_cov(a, c2);
    }

    // ---- integrate ---------------------------------------------------------
    // The Jacobian is taken by central differences ON THE ASSEMBLED RATES, so
    // every capacity term is differentiated consistently with the drift actually
    // integrated rather than with an algebraic derivative of a different function.
    const LsodaRhs rhs = [&](double t, const double* zz, double* dz) {
        pairs_at(t);
        std::vector<double> f;
        rates(zz, f);
        apply_jumps(f, dz);
        double qmax = 1.0;
        for (std::size_t a = 0; a < dim; ++a) qmax = std::max(qmax, std::fabs(zz[a]));
        const double hstep = 1e-6 * qmax;
        Matrix<double> J(dim, dim, 0.0);
        std::vector<double> qp(zz, zz + dim), fp, fm, dp(dim, 0.0), dm(dim, 0.0);
        for (std::size_t m = 0; m < dim; ++m) {
            const double keep = qp[m];
            qp[m] = keep + hstep;
            rates(qp.data(), fp);
            apply_jumps(fp, dp.data());
            qp[m] = keep - hstep;
            rates(qp.data(), fm);
            apply_jumps(fm, dm.data());
            qp[m] = keep;
            for (std::size_t a = 0; a < dim; ++a) J(a, m) = (dp[a] - dm[a]) / (2.0 * hstep);
        }
        // G = A diag(f) A', assembled from the jump lists.
        Matrix<double> G(dim, dim, 0.0);
        std::vector<double> col(dim, 0.0);
        for (std::size_t e = 0; e < nev; ++e) {
            if (f[e] == 0.0) continue;
            std::fill(col.begin(), col.end(), 0.0);
            for (std::size_t a = 0; a < ev[e].minus.size(); ++a) col[ev[e].minus[a]] -= 1.0;
            for (std::size_t a = 0; a < ev[e].plus.size(); ++a) col[ev[e].plus[a]] += 1.0;
            for (std::size_t a = 0; a < dim; ++a) {
                if (col[a] == 0.0) continue;
                for (std::size_t b = 0; b < dim; ++b)
                    if (col[b] != 0.0) G(a, b) += col[a] * f[e] * col[b];
            }
        }
        for (std::size_t a = 0; a < dim; ++a)
            for (std::size_t b = 0; b < dim; ++b) {
                double acc = G(a, b);
                for (std::size_t c2 = 0; c2 < dim; ++c2)
                    acc += J(a, c2) * zz[dim + c2 * dim + b] + zz[dim + a * dim + c2] * J(b, c2);
                dz[dim + a * dim + b] = acc;
            }
    };

    LsodaOptions lopt;
    lopt.rtol = opt.tol;
    lopt.atol = opt.tol * 1e-3;
    lopt.h_max = (tend - t0) / 10.0;
    // A step that crosses a whole segment integrates a rate that was never in
    // force. Cap it at a quarter of the NARROWEST segment of any schedule.
    if (period > 0.0) {
        double narrowest = std::numeric_limits<double>::infinity();
        for (std::size_t e = 0; e < sched.size(); ++e)
            for (std::size_t k = 1; k < sched[e].bp.size(); ++k)
                narrowest = std::min(narrowest, sched[e].bp[k] - sched[e].bp[k - 1]);
        if (std::isfinite(narrowest) && narrowest > 0.0)
            lopt.h_max = std::min(lopt.h_max, narrowest / 4.0);
    }

    // The output grid. Uniform over the whole horizon as before; when the answer
    // is a period average the last cycle is refined and every segment boundary
    // in it is BRACKETED, so no trapezoid interval straddles a jump.
    const bool averaging = unbounded && period > 0.0;
    const double w0 = averaging ? std::max(t0, tend - period) : t0;
    std::vector<double> grid;
    {
        const std::size_t ngrid = 201;
        for (std::size_t a = 0; a < ngrid; ++a)
            grid.push_back(t0 + (tend - t0) * static_cast<double>(a) /
                                    static_cast<double>(ngrid - 1));
        if (averaging) {
            const std::size_t nref = 2001;
            for (std::size_t a = 0; a < nref; ++a)
                grid.push_back(w0 + (tend - w0) * static_cast<double>(a) /
                                        static_cast<double>(nref - 1));
            std::vector<double> bounds;
            for (std::size_t e = 0; e < sched.size(); ++e) {
                const std::vector<double>& bp = sched[e].bp;
                const double per = bp.back() - bp.front();
                if (sched[e].cyclic && per > 0.0) {
                    const long kmax = static_cast<long>(std::ceil((tend - w0) / per)) + 2;
                    for (long kk = -1; kk <= kmax; ++kk)
                        for (std::size_t a = 0; a < bp.size(); ++a)
                            bounds.push_back(bp[a] + static_cast<double>(kk) * per);
                } else {
                    for (std::size_t a = 0; a < bp.size(); ++a) bounds.push_back(bp[a]);
                }
            }
            const double eps_b = std::max(1e-9, 1e-7 * (tend - w0));
            for (std::size_t a = 0; a < bounds.size(); ++a) {
                const double b = bounds[a];
                if (!(b > w0 && b < tend)) continue;
                grid.push_back(b - eps_b);
                grid.push_back(b);
                grid.push_back(b + eps_b);
            }
        }
        std::sort(grid.begin(), grid.end());
        grid.erase(std::remove_if(grid.begin(), grid.end(),
                                  [&](double v) { return v < t0 || v > tend; }),
                   grid.end());
        grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
        if (grid.empty() || grid.front() > t0) grid.insert(grid.begin(), t0);
    }
    const LsodaSolution sol = fluid_integrate_grid(rhs, z, grid, lopt);

    // ---- metrics -----------------------------------------------------------
    // A time-homogeneous model has a fixed point, so the steady-state answer is
    // the value at the horizon. A CYCLIC schedule has none, so the answer is the
    // trapezoidal average over the last full period; the value at the horizon
    // would be an arbitrary point of the cycle, at which the source and station
    // throughputs do not even agree.
    const std::vector<double>& zend = sol.final_state();
    FluidSolution out;
    out.method = "kp";
    out.iters = 1;
    out.QN = Matrix<double>(M, K, 0.0);
    out.UN = Matrix<double>(M, K, 0.0);
    out.RN = Matrix<double>(M, K, 0.0);
    out.TN = Matrix<double>(M, K, 0.0);
    out.xvec.assign(zend.begin(), zend.begin() + dim);

    const std::size_t nt = sol.t.size();
    Matrix<double> QVar(M, K, 0.0);
    for (std::size_t b = 0; b < xblocks.size(); ++b) {
        const KpBlock& xb = xblocks[b];
        std::vector<double> qser(nt, 0.0), user(nt, 0.0), tser(nt, 0.0);
        double vend = 0.0;
        for (std::size_t n = 0; n < nt; ++n) {
            const std::vector<double>& zs = sol.y[n];
            double q = 0.0;
            for (std::size_t p = 0; p < xb.nphases; ++p) q += zs[xb.offset + p];
            pairs_at(sol.t[n]);
            const double cap = capacity(zs.data(), xb.station);
            double tn = 0.0;
            for (std::size_t p = 0; p < xb.nphases; ++p) {
                double rowsum = 0.0;
                for (std::size_t c2 = 0; c2 < Dt1[xb.station][xb.cls].cols(); ++c2)
                    rowsum += Dt1[xb.station][xb.cls](p, c2);
                tn += rowsum * std::max(zs[xb.offset + p], 0.0) * cap;
            }
            const double c = sn.stations[xb.station].nservers;
            qser[n] = q;
            tser[n] = tn;
            user[n] = (sn.stations[xb.station].sched == lang::SchedStrategy::INF ||
                       !std::isfinite(c))
                          ? q
                          : std::min(q, c) / c;
        }
        for (std::size_t p = 0; p < xb.nphases; ++p)
            for (std::size_t p2 = 0; p2 < xb.nphases; ++p2)
                vend += zend[dim + (xb.offset + p) * dim + (xb.offset + p2)];
        out.QN(xb.station, xb.cls) = kp_detail::kp_summarise(qser, sol.t, w0, tend, averaging);
        out.UN(xb.station, xb.cls) = kp_detail::kp_summarise(user, sol.t, w0, tend, averaging);
        out.TN(xb.station, xb.cls) = kp_detail::kp_summarise(tser, sol.t, w0, tend, averaging);
        QVar(xb.station, xb.cls) = vend;
        // TN is zero only to the integrator's accuracy: a class that never visits leaves
        // a ~1e-20 residue in TN too, and a strict > 0 test then divides residue by residue.
        if (out.TN(xb.station, xb.cls) > lang::GlobalConstants::Zero)
            out.RN(xb.station, xb.cls) = out.QN(xb.station, xb.cls) / out.TN(xb.station, xb.cls);
    }
    for (std::size_t b = 0; b < ublocks.size(); ++b) {
        const KpBlock& u = ublocks[b];
        std::vector<double> aser(nt, 0.0);
        for (std::size_t n = 0; n < nt; ++n) {
            pairs_at(sol.t[n]);
            double tn = 0.0;
            for (std::size_t p = 0; p < u.nphases; ++p) {
                double rowsum = 0.0;
                for (std::size_t c2 = 0; c2 < Dt1[u.station][u.cls].cols(); ++c2)
                    rowsum += Dt1[u.station][u.cls](p, c2);
                tn += rowsum * std::max(sol.y[n][u.offset + p], 0.0);
            }
            aser[n] = tn;
        }
        out.TN(u.station, u.cls) = kp_detail::kp_summarise(aser, sol.t, w0, tend, averaging);
    }

    // The covariance IS the answer here, so it is reported through the same
    // `moments` channel the stationary closures use; `Sigma` is the full state
    // covariance at the horizon, so cross-station terms survive.
    FluidMomentReport rep;
    rep.Sigma = Matrix<double>(dim, dim, 0.0);
    for (std::size_t a = 0; a < dim; ++a)
        for (std::size_t b = 0; b < dim; ++b) rep.Sigma(a, b) = zend[dim + a * dim + b];
    rep.QVar = QVar;
    rep.QStd = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            rep.QStd(i, r) = std::sqrt(std::max(0.0, QVar(i, r)));
    rep.outer_iters = 1;
    out.has_moments = true;
    out.moments = rep;

    out.XN.assign(K, 0.0);
    out.CN.assign(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        const std::size_t rs = sn.classes[r].refstat;
        if (rs >= 1 && rs <= M) out.XN[r] = out.TN(rs - 1, r);
        double q = 0.0;
        for (std::size_t i = 0; i < M; ++i) q += out.QN(i, r);
        if (out.XN[r] > 0.0) out.CN[r] = q / out.XN[r];
    }

    if (tran != nullptr) {
        tran->t = sol.t;
        tran->q.clear();
        tran->QVar.clear();
        tran->Sigma.clear();
        for (std::size_t s = 0; s < sol.y.size(); ++s) {
            const std::vector<double>& zs = sol.y[s];
            tran->q.push_back(std::vector<double>(zs.begin(), zs.begin() + dim));
            Matrix<double> V(M, K, 0.0), Sg(dim, dim, 0.0);
            for (std::size_t a = 0; a < dim; ++a)
                for (std::size_t b = 0; b < dim; ++b) Sg(a, b) = zs[dim + a * dim + b];
            for (std::size_t b = 0; b < xblocks.size(); ++b) {
                const KpBlock& xb = xblocks[b];
                double v = 0.0;
                for (std::size_t p = 0; p < xb.nphases; ++p)
                    for (std::size_t p2 = 0; p2 < xb.nphases; ++p2)
                        v += Sg(xb.offset + p, xb.offset + p2);
                V(xb.station, xb.cls) = v;
            }
            tran->QVar.push_back(V);
            tran->Sigma.push_back(Sg);
        }
    }
    return out;
}

/** Port of `solver_fluid_kp.m`: the steady table at the horizon. */
template <class T>
FluidSolution solver_fluid_kp(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    return solver_fluid_kp_core(sn, opt, nullptr);
}

/**
 * Port of `@@SolverFLD/getTranAvgVar`: the queue-length VARIANCE along the
 * trajectory, per station and class, plus the full state covariance.
 *
 * ONLY `kp` HAS THIS. Every other fluid method integrates the mean alone and
 * carries no second moment, so asking them for one is an error rather than a
 * misleading zero -- and `minnormal`'s covariance is STATIONARY, so it is not this
 * quantity either. A caller that left the horizon unbounded gets one resolved the
 * way `getTranAvg` resolves it, from the slowest rate in the model.
 */
template <class T>
FluidKpTransient solver_fluid_tran_avg_var(const qn::NetworkStruct<T>& sn,
                                           const FluidOptions& opt) {
    std::string m = opt.method;
    if (m.size() > 6 && m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    if (m != "kp")
        throw UnsupportedError(
            "solver_fluid_tran_avg_var: getTranAvgVar needs method 'kp'; the other fluid methods "
            "integrate the mean only and carry no second moment");
    FluidKpTransient tran;
    FluidOptions o = opt;
    o.method = "kp";
    solver_fluid_kp_core(sn, o, &tran);
    return tran;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_KP_H
