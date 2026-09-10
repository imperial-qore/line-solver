/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_MATRIX_H
#define LINE_SOLVERS_FLUID_FLUID_MATRIX_H

/**
 * The `matrix` fluid method: a port of `solver_fluid_matrix.m`, the
 * formulation of Ruuskanen, Berg, Lehtinen et al., PEVA 151 (2021).
 *
 * WHY THIS METHOD EXISTS ALONGSIDE `closing`. Both integrate the same fluid
 * limit; they differ in how the drift is written. `closing` enumerates events
 * and sums their jumps. This one writes the whole drift as one matrix,
 *
 *     dx/dt = W' theta(x) + A_lambda,   W = Psi + B P A'
 *
 * where Psi is the block diagonal of the phase generators D0, B the column of
 * completion rates (the row sums of D1), P the station-to-station routing and
 * A the block diagonal of entry-phase vectors. `theta(x)` is the only
 * nonlinear part: it is the mass actually IN SERVICE,
 *
 *     theta = x / (sum of x at the station) * min(servers, sum of x)
 *
 * so it equals x while the station is under-loaded and saturates at the server
 * count above it. Writing the drift this way is what lets the p-norm smoothing
 * be applied in one place, and it is why MATLAB and native Python make this
 * the DEFAULT method.
 *
 * IT INTEGRATES ONCE, NOT ITERATIVELY. `closing` restarts from its own end
 * state repeatedly; this method picks a single horizon,
 * 10 * iter_max / min|W|, and integrates straight to it. Same fixed point,
 * different route.
 *
 * THE SOURCE IS NOT IN THE DYNAMICS. An EXT station's states are held at zero
 * with `theta = 0`, and its arrivals are injected directly into the phases of
 * the queues they route to, through `A_lambda`. The reference also strips the
 * routing back INTO a Source first, so an open class cannot cycle through it.
 * That is why the fluid state of an open model carries no Source mass here,
 * unlike the `closing` drift, which keeps a unit pool.
 *
 * DISABLED PAIRS. A (station, class) with no service still contributes one
 * placeholder column so the block structure lines up; the reference marks it
 * with NaN in A and then drops every NaN column from W. Reproduced, because
 * the surviving index order is what every result map is aligned to.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** The assembled matrix-form drift and the maps that read metrics off it. */
struct FluidMatrixSystem {
    Matrix<double> W;                 ///< (n x n) drift generator
    std::vector<double> alambda;      ///< (n) external arrivals per state
    std::vector<double> x0;           ///< (n) initial state
    std::vector<std::size_t> qa;      ///< (n) 0-based station of each state
    std::vector<double> sa;           ///< (n) server count of that station
    std::vector<bool> is_source;      ///< (n) state belongs to an EXT station
    std::vector<bool> is_inf;         ///< (n) state belongs to an INF station
    Matrix<double> sqc, suc, stc;     ///< (M*K x n) queue, utilization, throughput maps
    Matrix<double> src_arrival;       ///< (M x K) per-class arrival rate of each EXT station
    std::size_t nstates = 0;
    double min_rate = 1.0;            ///< smallest nonzero |W|, sets the horizon
    double pstar = 0.0;               ///< >0 selects the p-norm smoothing
    /**
     * When true, min(E[n], c) is replaced by E[min(n, c)] under the station's
     * equilibrium geometric marginal. Set ONLY by the degeneracy repair in
     * `solver_fluid`: min() is FLAT above the server count, so a network of
     * saturated stations has a CONTINUUM of fixed points and the integrator
     * returns whichever one it stopped at. See BUGS.md, _kb/06-solver-catalog.md.
     */
    bool var_closure = false;
};

namespace detail {

/**
 * Port of `compute_theta`: the mass in service.
 *
 * The FineTol in the denominator is the reference's, and is load bearing -- an
 * empty station would otherwise divide zero by zero and poison the drift.
 */
inline void fluid_matrix_theta(const FluidMatrixSystem& s, const double* x,
                               std::vector<double>& theta) {
    const std::size_t n = s.nstates;
    // The total mass at each station, which every state of it shares.
    std::vector<double> station_sum(s.sa.size(), 0.0);
    std::vector<double> sums(n, 0.0);
    std::size_t nst = 0;
    for (std::size_t i = 0; i < n; ++i) nst = std::max(nst, s.qa[i] + 1);
    station_sum.assign(nst, 0.0);
    for (std::size_t i = 0; i < n; ++i) station_sum[s.qa[i]] += x[i];
    for (std::size_t i = 0; i < n; ++i) {
        const double tot = lang::GlobalConstants::FineTol + station_sum[s.qa[i]];
        if (s.pstar > 0.0) {
            // p-norm smoothing: ghat stands in for min(1, c/tot) smoothly. An
            // INF station has k = infinity and so no min() to smooth; sa holds
            // the population there, which would read it as a k = N queue.
            const double c = s.sa[i];
            double g = 1.0;
            if (!s.is_inf[i] && tot > 0.0 && c > 0.0)
                g = 1.0 / std::pow(1.0 + std::pow(tot / c, s.pstar), 1.0 / s.pstar);
            if (std::isnan(g)) g = 0.0;
            theta[i] = x[i] * g;
        } else if (s.var_closure) {
            // E[min(n, c)] UNDER A GEOMETRIC MARGINAL, not min(E[n], c):
            //   n ~ Geometric(mean m) => E[min(n,c)] = sum_{k=1..c} p^k
            //                          = m * (1 - p^c),  p = m/(1+m).
            // Strictly increasing in m (slope 1/(1+m)^2 at c = 1, still 1e-2 at
            // m = 9, a restoring force the integrator can follow inside its
            // horizon) and with the same asymptotes, -> c as m -> inf and -> m
            // as m -> 0. An INF station has a server per job, so there is no
            // min() to close and sa holds the whole population there.
            double emin;
            if (s.is_inf[i]) {
                emin = tot;
            } else {
                const double p = tot > 0.0 ? tot / (1.0 + tot) : 0.0;
                emin = tot * (1.0 - std::pow(p, std::max(s.sa[i], 0.0)));
                if (!std::isfinite(emin)) emin = 0.0;
                emin = std::min(emin, tot);
            }
            theta[i] = x[i] / tot * emin;
        } else {
            theta[i] = x[i] / tot * std::min(s.sa[i], tot);
        }
        if (s.is_source[i]) theta[i] = 0.0;  // the Source is out of the dynamics
    }
}

}  // namespace detail

/** Assemble the matrix-form drift of `sn`. */
template <class T>
FluidMatrixSystem fluid_matrix_system(const qn::NetworkStruct<T>& sn,
                                      const std::vector<double>& init_sol, double pstar) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const FluidLayout L = fluid_layout(sn);
    FluidMatrixSystem s;
    s.pstar = pstar;

    // ---- station-to-station routing, via the stochastic complement ---------
    // sn.rt is over stateful nodes; a Router is stateful but not a station, so
    // the reference complements it out rather than indexing around it.
    const std::size_t S = sn.nof_stateful();
    std::vector<std::size_t> keep_idx;
    keep_idx.reserve(M * K);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
        for (std::size_t r = 0; r < K; ++r) keep_idx.push_back(isf * K + r);
    }
    Matrix<double> rt_full(S * K, S * K, 0.0);
    if (sn.rt.rows() == S * K)
        for (std::size_t a = 0; a < S * K; ++a)
            for (std::size_t b = 0; b < S * K; ++b)
                rt_full(a, b) = num_traits<T>::to_double(sn.rt(a, b));
    Matrix<double> P = mc::dtmc_stochcomp(rt_full, keep_idx);

    // Per-class arrival rate of each EXT station, and the removal of any
    // routing back INTO it: an open class leaves the Source once.
    Matrix<double> src_arrival(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched != lang::SchedStrategy::EXT) continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) continue;
            mam::Map<double> mp;
            const std::size_t nn = sn.service[i][r].D0.rows();
            mp.D0 = Matrix<double>(nn, nn, 0.0);
            mp.D1 = Matrix<double>(nn, nn, 0.0);
            for (std::size_t a = 0; a < nn; ++a)
                for (std::size_t b = 0; b < nn; ++b) {
                    mp.D0(a, b) = num_traits<T>::to_double(sn.service[i][r].D0(a, b));
                    mp.D1(a, b) = num_traits<T>::to_double(sn.service[i][r].D1(a, b));
                }
            const double mean = mam::map_mean(mp);
            if (mean > 0.0) src_arrival(i, r) = 1.0 / mean;
            if (src_arrival(i, r) > 0.0)
                for (std::size_t j = 0; j < M; ++j) {
                    if (j == i) continue;
                    for (std::size_t q = 0; q < K; ++q) P(j * K + q, i * K + r) = 0.0;
                }
        }
    }

    // ---- the block-structured state space ---------------------------------
    // One column per phase, or a single placeholder for a disabled pair. The
    // placeholder is dropped below; it exists so the blocks line up first.
    std::vector<std::size_t> blk_state;   // first state index of block (i,r)
    std::vector<std::size_t> blk_len;     // phases in block (i,r), 0 when disabled
    std::size_t nfull = 0;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            blk_state.push_back(nfull);
            const std::size_t p = L.kic[i][r];
            blk_len.push_back(p);
            nfull += (p == 0) ? 1 : p;
        }

    Matrix<double> Wfull(nfull, nfull, 0.0);
    std::vector<double> pie_full(nfull, 0.0), brate(nfull, 0.0);
    std::vector<bool> disabled_state(nfull, false);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t b = blk_state[i * K + r], p = blk_len[i * K + r];
            if (p == 0) {
                disabled_state[b] = true;  // the reference's NaN column
                continue;
            }
            const lang::Distrib<T>& d = sn.service[i][r];
            const std::vector<double> pie = detail::fluid_pie(d);
            for (std::size_t a = 0; a < p; ++a) {
                pie_full[b + a] = a < pie.size() ? pie[a] : 0.0;
                for (std::size_t c = 0; c < p; ++c)
                    Wfull(b + a, b + c) = num_traits<T>::to_double(d.D0(a, c));  // Psi
                double row = 0.0;
                for (std::size_t c = 0; c < d.D1.cols(); ++c)
                    row += num_traits<T>::to_double(d.D1(a, c));
                brate[b + a] = row;  // B
            }
        }
    // W += B P A': a completion at (i,r,a) routes to (j,l) and enters phase c.
    for (std::size_t ir = 0; ir < M * K; ++ir) {
        const std::size_t bi = blk_state[ir], pi = blk_len[ir];
        if (pi == 0) continue;
        for (std::size_t jl = 0; jl < M * K; ++jl) {
            const double p = P(ir, jl);
            if (!(p > 0.0)) continue;
            const std::size_t bj = blk_state[jl], pj = blk_len[jl];
            if (pj == 0) continue;
            for (std::size_t a = 0; a < pi; ++a)
                for (std::size_t c = 0; c < pj; ++c)
                    Wfull(bi + a, bj + c) += brate[bi + a] * p * pie_full[bj + c];
        }
    }

    // ---- external arrivals, injected at the queues, not the Source --------
    std::vector<double> alam_full(nfull, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched == lang::SchedStrategy::EXT) continue;
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t b = blk_state[i * K + r], p = blk_len[i * K + r];
            if (p == 0) continue;
            double rate = 0.0;
            for (std::size_t sidx = 0; sidx < M; ++sidx)
                if (src_arrival(sidx, r) > 0.0)
                    rate += src_arrival(sidx, r) * P(sidx * K + r, i * K + r);
            if (!(rate > 0.0)) continue;
            for (std::size_t a = 0; a < p; ++a) alam_full[b + a] = pie_full[b + a] * rate;
        }
    }

    // ---- drop the disabled placeholders and build the result maps ---------
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < nfull; ++i)
        if (!disabled_state[i]) keep.push_back(i);
    const std::size_t n = keep.size();
    s.nstates = n;
    s.W = Matrix<double>(n, n, 0.0);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) s.W(a, b) = Wfull(keep[a], keep[b]);
    s.alambda.assign(n, 0.0);
    s.x0.assign(n, 0.0);
    s.qa.assign(n, 0);
    s.sa.assign(n, 1.0);
    s.is_source.assign(n, false);
    s.is_inf.assign(n, false);
    s.sqc = Matrix<double>(M * K, n, 0.0);
    s.suc = Matrix<double>(M * K, n, 0.0);
    s.stc = Matrix<double>(M * K, n, 0.0);

    double closed_pop = 0.0;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population)) closed_pop += sn.classes[r].population;

    // Map every surviving state back to its (station, class, phase).
    std::vector<std::size_t> st_of(nfull, 0), cl_of(nfull, 0), ph_of(nfull, 0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t b = blk_state[i * K + r], p = blk_len[i * K + r];
            for (std::size_t a = 0; a < (p == 0 ? 1u : p); ++a) {
                st_of[b + a] = i;
                cl_of[b + a] = r;
                ph_of[b + a] = a;
            }
        }

    for (std::size_t a = 0; a < n; ++a) {
        const std::size_t f = keep[a], i = st_of[f], r = cl_of[f], k = ph_of[f];
        const double c = sn.stations[i].nservers;
        const double servers = std::isfinite(c) ? c : closed_pop;
        s.qa[a] = i;
        s.sa[a] = servers;
        s.alambda[a] = alam_full[f];
        s.is_source[a] = (sn.stations[i].sched == lang::SchedStrategy::EXT);
        s.is_inf[a] = !std::isfinite(c);
        s.sqc(i * K + r, a) = 1.0;
        s.suc(i * K + r, a) = (servers > 0.0) ? 1.0 / servers : 0.0;
        double row = 0.0;
        for (std::size_t cc = 0; cc < sn.service[i][r].D1.cols(); ++cc)
            row += num_traits<T>::to_double(sn.service[i][r].D1(k, cc));
        s.stc(i * K + r, a) = row;
        // The Source carries no mass: its arrivals enter downstream.
        s.x0[a] = s.is_source[a] ? 0.0 : (f < init_sol.size() ? 0.0 : 0.0);
    }

    // The initial state is given in the fluid layout's order (one entry per
    // enabled phase); map it onto the surviving states.
    if (!init_sol.empty()) {
        for (std::size_t a = 0; a < n; ++a) {
            const std::size_t f = keep[a], i = st_of[f], r = cl_of[f], k = ph_of[f];
            if (s.is_source[a] || !L.enabled[i][r]) continue;
            const std::size_t src = L.qidx[i][r] + k;
            if (src < init_sol.size()) s.x0[a] = init_sol[src];
        }
    }

    s.src_arrival = src_arrival;
    double mr = std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) {
            const double v = std::fabs(s.W(a, b));
            if (v > 0.0) mr = std::min(mr, v);
        }
    s.min_rate = std::isfinite(mr) ? mr : 1.0;
    return s;
}

/** The drift dx/dt = W' theta(x) + A_lambda. */
inline std::function<void(double, const double*, double*)> fluid_matrix_drift(
    const FluidMatrixSystem& s) {
    const std::size_t n = s.nstates;
    return [s, n](double, const double* x, double* dx) {
        std::vector<double> theta(n, 0.0);
        detail::fluid_matrix_theta(s, x, theta);
        for (std::size_t i = 0; i < n; ++i) {
            double acc = s.alambda[i];
            for (std::size_t j = 0; j < n; ++j) acc += s.W(j, i) * theta[j];  // W' theta
            dx[i] = acc;
        }
    };
}

/**
 * Is the returned point one of a CONTINUUM of fixed points?
 *
 * A station whose queue exceeds its server count has theta pinned at the server
 * count, so the drift cannot tell one split of the mass between two such
 * stations from another. The test is direct rather than structural: move a
 * little mass from one station to another along a POPULATION-CONSERVING
 * direction and see whether the drift moves at all. Both directions are tried,
 * because the integrator typically stops on the BOUNDARY of the degenerate set,
 * where one of the two does change the drift.
 *
 * The DIRECTIONAL DERIVATIVE is the scale-free quantity to threshold: a live
 * direction moves the drift at the station's own service rate (measured
 * 9.99e-01) and a null one only by the FineTol the share carries (measured
 * 1.00e-08), four orders apart. Source and INF states are excluded: neither can
 * be the pinned coordinate. Returns false whenever the point is not a fixed
 * point at all, so a transient run is never repaired.
 */
inline bool fluid_matrix_degenerate(
    const FluidMatrixSystem& s,
    const std::function<void(double, const double*, double*)>& drift,
    const std::vector<double>& x, std::size_t K) {
    const std::size_t n = s.nstates;
    if (x.size() != n || n == 0 || K == 0) return false;
    double max_abs_x = 1.0;
    for (std::size_t i = 0; i < n; ++i) max_abs_x = std::max(max_abs_x, std::fabs(x[i]));
    std::vector<double> d0(n, 0.0);
    drift(0.0, x.data(), d0.data());
    double max_d0 = 0.0;
    for (std::size_t i = 0; i < n; ++i) max_d0 = std::max(max_d0, std::fabs(d0[i]));
    if (max_d0 > 1e-6 * max_abs_x) return false;
    if (s.sqc.rows() % K != 0) return false;
    const std::size_t M = s.sqc.rows() / K;

    double rate_scale = 0.0;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b)
            rate_scale = std::max(rate_scale, std::fabs(s.W(a, b)));
    rate_scale = std::max(rate_scale, 1e-12);

    const double step = 1e-3 * max_abs_x;
    std::vector<double> xp(n, 0.0), dp(n, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        // PER CLASS, NOT PER STATION: a direction that moves a station's mass
        // across ALL its classes is infeasible where a SelfLoopingClass is
        // pinned at one station, and a well-posed model then reads as
        // degenerate. Source and INF states carry no min() to pin.
        std::vector<std::vector<std::size_t> > groups;
        for (std::size_t i = 0; i < M; ++i) {
            std::vector<std::size_t> members;
            for (std::size_t a = 0; a < n; ++a)
                if (!s.is_source[a] && !s.is_inf[a] && s.sqc(i * K + r, a) > 0.0)
                    members.push_back(a);
            if (!members.empty()) groups.push_back(members);
        }
        if (groups.size() < 2) continue;
        std::vector<double> mass(groups.size(), 0.0);
        for (std::size_t a = 0; a < groups.size(); ++a)
            for (std::size_t idx : groups[a]) mass[a] += x[idx];
        for (std::size_t a = 0; a < groups.size(); ++a) {
            if (mass[a] <= step) continue;
            for (std::size_t b = 0; b < groups.size(); ++b) {
                if (a == b) continue;
                xp = x;
                for (std::size_t idx : groups[a]) xp[idx] -= step * x[idx] / mass[a];
                for (std::size_t idx : groups[b])
                    xp[idx] += (mass[b] > 0.0) ? step * x[idx] / mass[b]
                                               : step / static_cast<double>(groups[b].size());
                drift(0.0, xp.data(), dp.data());
                double max_dd = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                    max_dd = std::max(max_dd, std::fabs(dp[i] - d0[i]));
                if (max_dd / step <= 1e-4 * rate_scale) return true;
            }
        }
    }
    return false;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_MATRIX_H
