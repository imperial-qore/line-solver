/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_ODES_H
#define LINE_SOLVERS_FLUID_FLUID_ODES_H

/**
 * The fluid drift: a port of `solver_fluid_odes.m` and the `ode_jumps_new` /
 * `ode_rate_base` / `ode_rates_closing` triple it drives on the DEFAULT
 * (`closing`) method.
 *
 * WHAT THE STATE VECTOR IS. One entry per (station, class, service phase),
 * laid out station-major then class-major then phase, and holding the mean
 * number of class-r jobs at station i currently in phase k. `q_indices(i,r)`
 * is where that block starts and `Kic(i,r)` is how long it is; a (station,
 * class) pair with no service contributes NO entries at all, which is why the
 * layout has to be computed rather than assumed to be M*K*phases.
 *
 * WHY THE DRIFT FACTORS THE WAY IT DOES. Every event is a single job moving:
 * it either completes at (i,c,ki) and starts at (j,l,kj), or it changes phase
 * within (i,c). So each event's effect on the state is a vector with one -1
 * and one +1, and the drift is
 *
 *     dx/dt = sum over events of jump_e * rate_e(x)
 *
 * The reference stores those jumps as a dense (nstates x nevents) matrix and
 * multiplies. THIS PORT STORES THE TWO INDICES INSTEAD, which is the same
 * arithmetic -- the matrix has exactly two nonzeros per column -- while
 * turning an O(nstates * nevents) product into O(nevents) per evaluation. On
 * a model with a few hundred states that is the difference between a fluid
 * solve dominated by the drift and one dominated by the integrator, and LSODA
 * evaluates the drift thousands of times.
 *
 * The rate of an event factors into a part fixed by the model and a part that
 * depends on the state:
 *
 *     rate_e(x) = rateBase_e * g(x)_{eventIdx_e}
 *
 * `rateBase` folds the phase completion rate, the routing probability and the
 * destination's entry-phase probability into one number computed once;
 * `g(x)` is where the scheduling lives, and is the only thing re-evaluated.
 * That split is the reference's and is what makes the drift cheap.
 *
 * WHAT g(x) IS, PER DISCIPLINE (`ode_rates_closing_factors`). g starts as x
 * itself, which is already right for a delay and for any station whose
 * population is below its server count, and is then corrected:
 *   INF      nothing, unless a load-dependent alpha(n_i) scales the station
 *   EXT      the source keeps unit mass per class, so the first phase absorbs
 *            whatever the other phases do not hold
 *   PS/FCFS  the servers are shared: the block is scaled by psi(n_i)/n_i, with
 *            psi = min(n_i,c)*alpha(n_i) the work the station clears
 *   DPS      the same capacity psi, split by the WEIGHTED share w_j x_j / sum
 *   GPS      the server is split by weight among the BACKLOGGED classes, then
 *            equally among that class's own jobs
 * Any other discipline is left as x, exactly as the reference leaves it -- which
 * means it is integrated as an infinite server, and is why the featset gate
 * refuses the disciplines that have no branch here.
 *
 * ALL THREE NON-LINEAR TERMS ABOVE TAKE A SECOND MOMENT WHEN ONE IS SUPPLIED
 * (`FluidClosure`): min() through `fluid_capacity_closure`, the PS/DPS ratio
 * through `fluid_share_closure`, the GPS indicator through `fluid_gps_share`.
 * With no closure every one of them collapses to its value at the mean, so the
 * first-order methods take exactly the same code path.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_closures.h"
#include "line/util/error.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** Where each (station, class) block sits in the state vector. */
struct FluidLayout {
    std::size_t nstates = 0;                        ///< length of the state vector
    std::vector<std::vector<std::size_t>> qidx;     ///< 0-based first index of (i,r)
    std::vector<std::vector<std::size_t>> kic;      ///< phases held by (i,r); 0 when disabled
    std::vector<std::vector<bool>> enabled;         ///< whether (i,r) is served at all
};

/**
 * One event of the drift.
 *
 * `minus` and `plus` are the state entries the event moves a job out of and
 * into. They can coincide -- a self-routing class completing and restarting in
 * the same phase -- and then the event contributes nothing, which is exactly
 * what the reference's jump column of all zeros contributes.
 */
struct FluidEvent {
    std::size_t minus = 0;
    std::size_t plus = 0;
    std::size_t event_idx = 0;  ///< state entry whose g(x) drives this rate
    double rate_base = 0.0;     ///< the model-fixed part of the rate
};

/**
 * The second moment the drift closes its non-linear terms with, i.e.
 * `options.config.moment_sigma2` and `options.config.moment_cov`.
 *
 * EMPTY IS THE FIRST-ORDER CLOSURE, and not a missing input: every closure
 * evaluates at the mean when it is given no variance, which is what the
 * `closing`, `matrix`, `statedep`, `softmin` and `tbi` methods do. Only
 * `solver_fluid_moments` fills this, and it fills it from the covariance it
 * converged to rather than from the one its own solve produced -- see there for
 * why the two must be the same one.
 *
 * `sigma2[i]` is the population variance of station i and closes its min(); the
 * capacity SHARE is a ratio of coordinates, so closing it needs `cov[i]`, the
 * covariance block over station i's own coordinates, and not just their total.
 */
struct FluidClosure {
    std::vector<double> sigma2;           ///< per station; empty selects first order
    std::vector<Matrix<double>> cov;      ///< per station, 0x0 keeps the plug-in share
    /** `any(sigma2 > 0)`, the reference's GLOBAL gaussian flag. */
    bool gaussian() const {
        for (std::size_t i = 0; i < sigma2.size(); ++i)
            if (sigma2[i] > 0.0) return true;
        return false;
    }
    double sigma2_of(std::size_t i) const { return i < sigma2.size() ? sigma2[i] : 0.0; }
    const Matrix<double>* cov_of(std::size_t i) const {
        if (i >= cov.size() || cov[i].rows() == 0) return nullptr;
        return &cov[i];
    }
};

/** The assembled drift: the layout, the events, and the per-station schedule. */
/**
 * Port of `solver_fluid_ratemult.m`'s output: a per-event MULTIPLIER on a time
 * grid, evaluated by `fluid_interpcols`.
 *
 * WHY A MULTIPLIER AND NOT A SECOND DRIFT. The closing rate is
 * `rate = rate_base .* theta(x)` and `rate_base` is LINEAR in the station-class
 * service or arrival rate, so any time variation of that rate -- an NHPP source
 * intensity, an inter-layer demand trajectory injected by the coupled LN
 * transient -- reduces exactly to a scalar factor per event per instant. That is
 * what lets three independent sources compose by elementwise product on the
 * union of their grids instead of each needing its own drift.
 *
 * EMPTY IS THE AUTONOMOUS DRIFT, not a missing input. `fluid_drift` skips the
 * lookup entirely when there is no multiplier, so a model with no time-varying
 * source integrates exactly the function it integrated before this existed.
 */
struct FluidRateMult {
    std::vector<double> tgrid;   ///< strictly increasing, one column per entry
    Matrix<double> Mmat;         ///< (nevents x ngrid)
    bool empty() const { return tgrid.empty() || Mmat.cols() == 0; }
};

/**
 * Port of `fluid_interpcols.m`: clamped piecewise-linear interpolation of the
 * columns of `B` at a scalar time, into `out`.
 *
 * CLAMPED, not extrapolated: a time outside the grid takes the boundary column
 * (a zero-order hold). A trajectory supplied over [0, T] therefore holds its
 * last value if the integrator steps past T, rather than continuing a linear
 * ramp to a rate the caller never declared.
 */
inline void fluid_interpcols(const std::vector<double>& tg, const Matrix<double>& B, double tt,
                             std::vector<double>& out) {
    const std::size_t nr = B.rows();
    out.assign(nr, 1.0);
    if (tg.empty() || B.cols() == 0) return;
    if (tt <= tg.front()) {
        for (std::size_t i = 0; i < nr; ++i) out[i] = B(i, 0);
        return;
    }
    if (tt >= tg.back()) {
        for (std::size_t i = 0; i < nr; ++i) out[i] = B(i, B.cols() - 1);
        return;
    }
    std::size_t j = 0;
    for (std::size_t k = 0; k < tg.size(); ++k)
        if (tg[k] <= tt) j = k;
    if (j + 1 >= tg.size()) {
        for (std::size_t i = 0; i < nr; ++i) out[i] = B(i, j);
        return;
    }
    const double w = (tt - tg[j]) / (tg[j + 1] - tg[j]);
    for (std::size_t i = 0; i < nr; ++i) out[i] = (1.0 - w) * B(i, j) + w * B(i, j + 1);
}

struct FluidOdeSystem {
    FluidLayout layout;
    std::vector<FluidEvent> events;
    /**
     * How many leading entries of `events` are DEPARTURES (a job completing at
     * one block and starting at another); the rest are phase changes within a
     * block. The passage-time construction needs the distinction, and it
     * cannot be recovered by inspecting the indices -- a self-routing class
     * produces a departure whose two endpoints lie in the same block, exactly
     * like a phase change.
     */
    std::size_t n_departures = 0;
    std::vector<lang::SchedStrategy> sched;      ///< per station
    std::vector<double> nservers;                ///< per station, already finite
    std::vector<std::vector<double>> weight;     ///< per station, per class (DPS, GPS)
    /**
     * `sn.lldscaling(i,:)` per station, EMPTY when the station has none or when
     * every entry is one -- the reference's `if all(lldrow == 1), lldrow = []`,
     * which keeps a station without load dependence on the plain branch rather
     * than on an interpolation that would return 1 anyway.
     */
    std::vector<std::vector<double>> lld;
    /** The second moment the closures read; empty is the first-order drift. */
    FluidClosure closure;
    /** `solver_fluid_ratemult`'s multiplier; empty is the autonomous drift. */
    FluidRateMult ratemult;
};

namespace detail {

/** True when the (station, class) pair has a usable service process. */
template <class T>
bool fluid_service_defined(const qn::NetworkStruct<T>& sn, std::size_t i, std::size_t r) {
    if (sn.disabled[i][r]) return false;
    const lang::Distrib<T>& d = sn.service[i][r];
    return d.D0.rows() > 0 && d.D0.rows() == d.D0.cols();
}

/**
 * mu and phi of a phase-type service process, as `Markovian.getMu`/`getPhi`
 * define them: mu(k) is the total rate out of phase k and phi(k) the share of
 * it that is a completion. The reference's guard for a zero diagonal -- an
 * Immediate distribution -- is reproduced, since dividing by it would produce
 * a NaN that then poisons the whole drift.
 */
template <class T>
void fluid_mu_phi(const lang::Distrib<T>& d, std::vector<double>& mu, std::vector<double>& phi) {
    const std::size_t n = d.D0.rows();
    mu.assign(n, 0.0);
    phi.assign(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const double d0kk = num_traits<T>::to_double(d.D0(k, k));
        double rowD1 = 0.0;
        for (std::size_t j = 0; j < d.D1.cols(); ++j) rowD1 += num_traits<T>::to_double(d.D1(k, j));
        mu[k] = -d0kk;
        phi[k] = (d0kk == 0.0) ? 1.0 : rowD1 / (-d0kk);
    }
}

/** The entry-phase distribution of a service process, `map_pie`. */
template <class T>
std::vector<double> fluid_pie(const lang::Distrib<T>& d) {
    const std::size_t n = d.D0.rows();
    if (n == 0) return std::vector<double>{1.0};
    if (n == 1) return std::vector<double>{1.0};
    mam::Map<double> m;
    m.D0 = Matrix<double>(n, n, 0.0);
    m.D1 = Matrix<double>(n, n, 0.0);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) {
            m.D0(a, b) = num_traits<T>::to_double(d.D0(a, b));
            m.D1(a, b) = num_traits<T>::to_double(d.D1(a, b));
        }
    return mam::map_pie(m);
}

}  // namespace detail

/**
 * Port of the layout half of `solver_fluid_odes.m`.
 *
 * The cumulative index runs even over disabled pairs, which contribute zero
 * phases: the reference records `q_indices` for them too, so that a later
 * lookup of a disabled pair lands on the next block's start rather than out of
 * range. Reproduced deliberately.
 */
template <class T>
FluidLayout fluid_layout(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    FluidLayout L;
    L.qidx.assign(M, std::vector<std::size_t>(K, 0));
    L.kic.assign(M, std::vector<std::size_t>(K, 0));
    L.enabled.assign(M, std::vector<bool>(K, false));
    std::size_t cursor = 0;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            L.qidx[i][r] = cursor;
            if (detail::fluid_service_defined(sn, i, r)) {
                L.enabled[i][r] = true;
                L.kic[i][r] = sn.service[i][r].D0.rows();
            }
            cursor += L.kic[i][r];
        }
    L.nstates = cursor;
    return L;
}

/**
 * Build the drift of `sn`: the port of `ode_jumps_new` and `ode_rate_base`
 * fused into one pass.
 *
 * The two reference functions walk the same nested loops in the same order and
 * their outputs are matched element by element, so building them together is
 * the only way to keep them aligned by construction rather than by comment.
 *
 * `rt` is read in STATION space, by the STOCHASTIC COMPLEMENT of `sn.rt` on the
 * station rows. `sn.rt` is indexed by stateful node, and the reference indexes
 * it with `(i-1)*K+c` for i over stations, which is the same thing only when
 * every stateful node is a station. Mapping the station index through
 * `stateful_of_station` is not enough either: it reads the DIRECT station pair
 * and so drops every path that traverses a non-station stateful node. On a
 * cache model that is the whole flow -- Think -> Cache -> Queue has a zero
 * direct entry -- and the drift then has no outgoing route at all, so the ODE
 * returns its initial condition with every job parked at the reference station.
 * Absorbing those nodes is exact because they hold no jobs.
 */
template <class T>
FluidOdeSystem fluid_ode_system(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    FluidOdeSystem sys;
    sys.layout = fluid_layout(sn);
    const FluidLayout& L = sys.layout;

    // The service processes, lowered once to plain doubles.
    std::vector<std::vector<std::vector<double>>> mu(M, std::vector<std::vector<double>>(K));
    std::vector<std::vector<std::vector<double>>> phi(M, std::vector<std::vector<double>>(K));
    std::vector<std::vector<std::vector<double>>> pie(M, std::vector<std::vector<double>>(K));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) {
                pie[i][r] = std::vector<double>{1.0};
                continue;
            }
            detail::fluid_mu_phi(sn.service[i][r], mu[i][r], phi[i][r]);
            pie[i][r] = detail::fluid_pie(sn.service[i][r]);
        }

    // Station-space routing, rt((i,c) -> (j,l)), through the stochastic
    // complement of the stateful-indexed matrix on the station rows.
    const std::size_t S = sn.nof_stateful();
    const bool have_rt = sn.rt.rows() == S * K;
    std::vector<std::size_t> keep_idx;
    keep_idx.reserve(M * K);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
        for (std::size_t c = 0; c < K; ++c) keep_idx.push_back(isf * K + c);
    }
    Matrix<double> rt_st(M * K, M * K, 0.0);
    if (have_rt) {
        Matrix<double> rt_full(S * K, S * K, 0.0);
        for (std::size_t a = 0; a < S * K; ++a)
            for (std::size_t b = 0; b < S * K; ++b)
                rt_full(a, b) = num_traits<T>::to_double(sn.rt(a, b));
        rt_st = mc::dtmc_stochcomp(rt_full, keep_idx);
    }
    const auto route = [&](std::size_t i, std::size_t c, std::size_t j, std::size_t l) -> double {
        if (!have_rt) return 0.0;
        return rt_st(i * K + c, j * K + l);
    };

    // ---- departures: (i,c,ki) completes and the job starts at (j,l,kj) -----
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            if (!L.enabled[i][c]) continue;
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t l = 0; l < K; ++l) {
                    const double p = route(i, c, j, l);
                    if (!(p > 0.0)) continue;
                    for (std::size_t ki = 0; ki < L.kic[i][c]; ++ki)
                        for (std::size_t kj = 0; kj < L.kic[j][l]; ++kj) {
                            FluidEvent e;
                            e.minus = L.qidx[i][c] + ki;
                            e.plus = L.qidx[j][l] + kj;
                            e.event_idx = L.qidx[i][c] + ki;
                            const double pj = kj < pie[j][l].size() ? pie[j][l][kj] : 0.0;
                            e.rate_base = phi[i][c][ki] * mu[i][c][ki] * p * pj;
                            sys.events.push_back(e);
                        }
                }
        }

    sys.n_departures = sys.events.size();

    // ---- phase changes within (i,c): rate is the off-diagonal of D0 --------
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            if (!L.enabled[i][c]) continue;
            const lang::Distrib<T>& d = sn.service[i][c];
            // EVERY source phase, the last included. Bounding ki at kic-1 is
            // valid only for an acyclic PH and drops the LAST ROW of D0, which
            // a general MAP or an MMPP2 carries: on the 2-phase MAP
            // D0=[-5 1; 2 -4], D1=[3 1; 1 1] the missing row cost 14% of the
            // arrival rate (2.7586 against the exact 3.2). MATLAB
            // ode_jumps_new and the JAR PassageTimeODE both iterate all rows.
            for (std::size_t ki = 0; ki < L.kic[i][c]; ++ki)
                for (std::size_t kp = 0; kp < L.kic[i][c]; ++kp) {
                    if (kp == ki) continue;
                    FluidEvent e;
                    e.minus = L.qidx[i][c] + ki;
                    e.plus = L.qidx[i][c] + kp;
                    e.event_idx = L.qidx[i][c] + ki;
                    e.rate_base = num_traits<T>::to_double(d.D0(ki, kp));
                    sys.events.push_back(e);
                }
        }

    // ---- per-station scheduling data --------------------------------------
    sys.sched.resize(M);
    sys.nservers.resize(M);
    sys.weight.assign(M, std::vector<double>(K, 1.0));
    double closed_pop = 0.0;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population)) closed_pop += sn.classes[r].population;
    for (std::size_t i = 0; i < M; ++i) {
        sys.sched[i] = sn.stations[i].sched;
        const double c = sn.stations[i].nservers;
        // A delay has infinitely many servers; the reference substitutes the
        // closed population, which is the most that can ever be in service.
        sys.nservers[i] = std::isfinite(c) ? c : closed_pop;
        if (sn.stations[i].sched == lang::SchedStrategy::DPS ||
            sn.stations[i].sched == lang::SchedStrategy::GPS)
            for (std::size_t r = 0; r < K && r < sn.stations[i].schedparam.size(); ++r)
                sys.weight[i][r] = num_traits<T>::to_double(sn.stations[i].schedparam[r]);
    }

    // Load dependence, lowered per station and dropped where it is the identity.
    sys.lld.assign(M, std::vector<double>());
    for (std::size_t i = 0; i < M; ++i) {
        const std::vector<T>& row = sn.stations[i].lldscaling;
        bool all_one = true;
        for (std::size_t k = 0; k < row.size(); ++k)
            if (num_traits<T>::to_double(row[k]) != 1.0) all_one = false;
        if (row.empty() || all_one) continue;
        sys.lld[i].resize(row.size());
        for (std::size_t k = 0; k < row.size(); ++k)
            sys.lld[i][k] = num_traits<T>::to_double(row[k]);
    }
    return sys;
}

/**
 * The reference's dense jump matrix D, (nstates x nevents), rebuilt from the
 * two-index event form this port stores instead.
 *
 * `solver_fluid_odes.m` discards D once it has composed the right-hand side and
 * the drift never needs it back, which is why the drift does not keep it. The
 * covariance equation of the linear noise approximation does: its diffusion
 * matrix is D*diag(r(x))*D', which cannot be recovered from F alone, and its
 * reachable subspace is range(D). Built on demand, in the event order
 * `fluid_ode_system` emits, so column e is event e.
 */
inline Matrix<double> fluid_jump_matrix(const FluidOdeSystem& sys) {
    Matrix<double> D(sys.layout.nstates, sys.events.size(), 0.0);
    for (std::size_t e = 0; e < sys.events.size(); ++e) {
        D(sys.events[e].minus, e) -= 1.0;
        D(sys.events[e].plus, e) += 1.0;
    }
    return D;
}

/**
 * Port of `ode_rates_closing_factors`: the state-dependent factor g(x), in place.
 *
 * `g` must already be a copy of `x` on entry, which is the reference's
 * `rates = x` and is what makes INF and every under-loaded station correct
 * with no work.
 *
 * THE CLOSURE VARIANCE IS READ FROM THE SYSTEM, not passed separately, because
 * `fluid_drift` closes over the system alone and LSODA holds that callback for
 * the whole integration; a closure supplied beside it would have to be captured
 * somewhere else and could then disagree with the Jacobian, which reads it from
 * here.
 *
 * WHAT THE `gaussian` FLAG IS. The reference's `any(sigma2 > 0)` is GLOBAL over
 * the stations, not per station, so one station carrying a variance sends every
 * PS/FCFS station down the closure branch -- with its own sigma2(i), which may be
 * zero, in which case `fluid_capacity_closure` returns min(n_i,c) and the branch
 * reproduces the first-order value. Reproduced as the same global flag so that
 * the one place the two differ, the zero floor on a negative E[min], is reached
 * on the same models as in the reference.
 */
inline void fluid_rates_closing_factors(const FluidOdeSystem& sys, const double* x,
                                        std::vector<double>& g) {
    const FluidLayout& L = sys.layout;
    const std::size_t M = L.qidx.size();
    const std::size_t K = M ? L.qidx[0].size() : 0;
    const FluidClosure& cl = sys.closure;
    const bool gaussian = cl.gaussian();

    for (std::size_t i = 0; i < M; ++i) {
        const std::vector<double>& lld = sys.lld[i];
        const double s2 = cl.sigma2_of(i);
        const Matrix<double>* Ci = cl.cov_of(i);
        const std::size_t blo = K ? L.qidx[i][0] : 0;
        const std::size_t bhi = K ? L.qidx[i][K - 1] + L.kic[i][K - 1] : 0;  // one past the end
        switch (sys.sched[i]) {
            case lang::SchedStrategy::INF: {
                // Without load dependence each job is served at its own rate and
                // the share is the identity; alpha(n_i) scales the whole station.
                if (lld.empty()) break;
                double ni = 0.0;
                for (std::size_t p = blo; p < bhi; ++p) ni += x[p];
                if (!(ni > 0.0)) break;
                const ClosureValue h = fluid_capacity_closure(ni, sys.nservers[i], s2, lld, true);
                for (std::size_t p = blo; p < bhi; ++p) g[p] = x[p] / ni * h.h;
                break;
            }
            case lang::SchedStrategy::EXT: {
                // The source holds unit mass per class at all times: phase one
                // carries whatever the later phases do not.
                for (std::size_t k = 0; k < K; ++k) {
                    if (!L.enabled[i][k]) continue;
                    const std::size_t b = L.qidx[i][k], n = L.kic[i][k];
                    double rest = 0.0;
                    for (std::size_t p = 1; p < n; ++p) rest += x[b + p];
                    g[b] = 1.0 - rest;
                }
                break;
            }
            case lang::SchedStrategy::PS:
            case lang::SchedStrategy::FCFS: {
                if (K == 0) break;
                double ni = 0.0;
                for (std::size_t p = blo; p < bhi; ++p) ni += x[p];
                if ((gaussian || !lld.empty()) && ni > 0.0) {
                    const ClosureValue h =
                        fluid_capacity_closure(ni, sys.nservers[i], s2, lld, false);
                    if (Ci == nullptr) {
                        for (std::size_t p = blo; p < bhi; ++p) g[p] = x[p] / ni * h.h;
                    } else {
                        // THE SHARE AND THE CAPACITY ARE CLOSED JOINTLY. What
                        // the station clears is S_j*psi(N), and both factors move
                        // with N, so the product needs Cov(S_j,N)*psi'(n) on top
                        // of the two separate closures; see `fluid_share_closure`.
                        // With unit weights this is the DPS branch below.
                        const std::size_t nb = bhi - blo;
                        std::vector<double> xb(x + blo, x + bhi), wv(nb, 1.0);
                        const ShareValue sh = fluid_share_closure(xb, wv, *Ci, false, true);
                        std::vector<double> rb(nb, 0.0);
                        for (std::size_t p = 0; p < nb; ++p)
                            rb[p] = sh.s[p] * h.h + h.dh * sh.cn[p];
                        fluid_project_rate(rb, xb, lld.empty(), h.h);
                        for (std::size_t p = 0; p < nb; ++p) g[blo + p] = rb[p];
                    }
                } else if (ni > sys.nservers[i]) {  // min = ni is handled by g = x
                    const double s = sys.nservers[i] / ni;
                    for (std::size_t p = blo; p < bhi; ++p) g[p] = x[p] * s;
                }
                break;
            }
            case lang::SchedStrategy::DPS: {
                // DPS is PS with a weighted share: the class-k coordinates get
                // w_k*x/xi of the station capacity psi(xi) instead of x/xi of it.
                //
                // THE DENOMINATOR CARRIES NO ADDITIVE GUARD. It used to seed the
                // sum with mean(w) to keep the ratio finite on an empty station;
                // that term never cancels, so the shares summed to
                // 1 - mean(w)/xi instead of 1 and the utilization was depressed
                // by that factor. The capacity was also taken as the full server
                // count rather than psi(xi), so an underloaded station was served
                // at full rate. Both are the PS/FCFS branch's rules here, guarded
                // by xi > 0 the way that branch guards, and equal weights now
                // reduce DPS to PS identically.
                if (K == 0) break;
                double wsum = 0.0;
                for (std::size_t k = 0; k < K; ++k) wsum += sys.weight[i][k];
                if (wsum <= 0.0) break;
                const std::size_t nb = bhi - blo;
                std::vector<double> wv(nb, 0.0);
                for (std::size_t k = 0; k < K; ++k) {
                    if (!L.enabled[i][k]) continue;
                    const std::size_t b = L.qidx[i][k] - blo;
                    for (std::size_t p = 0; p < L.kic[i][k]; ++p)
                        wv[b + p] = sys.weight[i][k] / wsum;
                }
                double xi = 0.0, wx = 0.0;
                for (std::size_t p = 0; p < nb; ++p) {
                    xi += x[blo + p];
                    wx += wv[p] * x[blo + p];
                }
                if (!(xi > 0.0) || !(wx > 0.0)) break;
                const ClosureValue psi =
                    fluid_capacity_closure(xi, sys.nservers[i], s2, lld, false);
                const std::vector<double> xb(x + blo, x + bhi);
                const ShareValue sh = fluid_share_closure(
                    xb, wv, Ci ? *Ci : Matrix<double>(0, 0, 0.0), false, true);
                std::vector<double> rb(nb, 0.0);
                for (std::size_t p = 0; p < nb; ++p)
                    rb[p] = sh.s[p] * psi.h + psi.dh * sh.cn[p];
                fluid_project_rate(rb, xb, lld.empty(), psi.h);
                for (std::size_t p = 0; p < nb; ++p) g[blo + p] = rb[p];
                break;
            }
            case lang::SchedStrategy::GPS: {
                // GPS splits the server by WEIGHT among the BACKLOGGED classes,
                // then equally among that class's own jobs. The share is a
                // function of the backlog indicator, so `fluid_gps_share` closes
                // it over the 2^K patterns using P(X_k >= 1). No capacity term
                // multiplies it: GPS is single-server and the indicator already
                // carries the idle server, so the shares sum to
                // 1 - P(station empty) by design.
                if (K == 0) break;
                if (sys.nservers[i] > 1.0)
                    throw UnsupportedError(
                        "ode_rates_closing_factors: multi-server GPS stations are not supported, as "
                        "in the reference: the backlog closure splits ONE server by weight");
                std::vector<double> xk(K, 0.0), vk(K, 0.0), wk(K, 0.0);
                for (std::size_t k = 0; k < K; ++k) {
                    wk[k] = sys.weight[i][k];
                    if (!L.enabled[i][k]) continue;
                    const std::size_t b = L.qidx[i][k], n = L.kic[i][k];
                    for (std::size_t p = 0; p < n; ++p) xk[k] += x[b + p];
                    if (Ci == nullptr) continue;
                    double v = 0.0;
                    for (std::size_t p = 0; p < n; ++p)
                        for (std::size_t q = 0; q < n; ++q)
                            v += (*Ci)(b - blo + p, b - blo + q);
                    vk[k] = std::max(0.0, v);
                }
                const ShareValue sk = fluid_gps_share(xk, wk, vk, false);
                double a = 1.0;
                if (!lld.empty()) {
                    double ni = 0.0;
                    for (std::size_t p = blo; p < bhi; ++p) ni += x[p];
                    a = fluid_lld_scaling(lld, ni).h;
                }
                for (std::size_t k = 0; k < K; ++k) {
                    if (!L.enabled[i][k] || !(xk[k] > 0.0)) continue;
                    const std::size_t b = L.qidx[i][k], n = L.kic[i][k];
                    for (std::size_t p = 0; p < n; ++p)
                        g[b + p] = x[b + p] / xk[k] * sk.s[k] * a;
                }
                break;
            }
            default:
                break;  // as the reference leaves it: g = x
        }
    }
}

/** The reference's `ode_rates_closing` name, kept for the first-order callers. */
inline void fluid_rates_closing(const FluidOdeSystem& sys, const double* x, std::vector<double>& g) {
    fluid_rates_closing_factors(sys, x, g);
}

/**
 * The drift dx/dt, ready to hand to the integrator.
 *
 * Returned by value as a closure over a copy of the system, so the caller can
 * let the builder go out of scope; LSODA holds the callback for the whole
 * integration.
 */
inline std::function<void(double, const double*, double*)> fluid_drift(const FluidOdeSystem& sys) {
    const std::size_t n = sys.layout.nstates;
    return [sys, n](double t, const double* x, double* dx) {
        std::vector<double> g(x, x + n);
        fluid_rates_closing(sys, x, g);
        for (std::size_t i = 0; i < n; ++i) dx[i] = 0.0;
        // The autonomous arm is kept SEPARATE rather than multiplied by a vector
        // of ones: the multiplier is absent on every model that has no
        // time-varying source, and a per-step interpolation there would be paid
        // by every fluid solve in the port for nothing.
        if (sys.ratemult.empty()) {
            for (const FluidEvent& e : sys.events) {
                const double r = e.rate_base * g[e.event_idx];
                if (r == 0.0) continue;
                dx[e.minus] -= r;
                dx[e.plus] += r;
            }
            return;
        }
        std::vector<double> mult;
        fluid_interpcols(sys.ratemult.tgrid, sys.ratemult.Mmat, t, mult);
        for (std::size_t k = 0; k < sys.events.size(); ++k) {
            const FluidEvent& e = sys.events[k];
            const double m = k < mult.size() ? mult[k] : 1.0;
            const double r = m * e.rate_base * g[e.event_idx];
            if (r == 0.0) continue;
            dx[e.minus] -= r;
            dx[e.plus] += r;
        }
    };
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_ODES_H
