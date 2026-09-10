/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_MOMENTS_H
#define LINE_SOLVERS_FLUID_FLUID_MOMENTS_H

/**
 * The second-order fluid methods: `fluid_moment_terms.m`, `fluid_lyapunov.m`,
 * `fluid_drift_jacobian.m`, `fluid_refine_meanfield.m` and
 * `solver_fluid_moments.m`, which back `options.method` `minnormal` and
 * `refined`.
 *
 * WHY THE PORT NEEDED THESE AT ALL, given that the first-order methods already
 * answered every model. `minnormal` is what the REFERENCE'S `default` resolves to
 * wherever it applies (`fluid_resolve_default_method.m`), so without it the port
 * answered a different method than the reference under the same name -- and
 * answered it less accurately, since the first-order closure replaces
 * E[min(X,c)] by min(E[X],c) and is worst exactly at rho ~ 1. On the reference's
 * own sweep (Delay(Z=1) -> Queue(PS, c=2), N=6) the exact CTMC queue length is
 * 1.95137, `closing` returns 2.00000 and `minnormal` 1.96063.
 *
 * WHAT THE SECOND MOMENT IS, and why a fluid solver has one. The closing ODEs are
 * a density-dependent Markov population process
 *
 *     dx/dt = F(x) = D r(x),     r_e(x) = rateBase_e g_e(x)
 *
 * whose fluctuation process Z = X - x* obeys, to leading order,
 * dZ = A Z dt + sqrt(D diag(r) D') dW with A = dF/dx. That linear noise
 * approximation has a stationary covariance, the solution of the Lyapunov
 * equation A Sigma + Sigma A' + D diag(r) D' = 0, and THAT is the second moment
 * reported through `getMoments`. `solver_fluid_odes.m` throws D and r away once
 * it has composed F, which is why `fluid_moment_terms` rebuilds them: the
 * diffusion matrix cannot be recovered from F alone.
 *
 * THE TWO METHODS DIFFER IN WHICH FIXED POINT THEY EXPAND ABOUT, and mixing them
 * would count the same term twice. `minnormal` solves mean and covariance
 * self-consistently, so its fixed point already RESUMS the O(1/N) correction --
 * expanding E[F(X)] to second order and setting it to zero reproduces the Gast
 * correction equation exactly. `refined` therefore recomputes the base point with
 * the FIRST-order closure and adds the correction to that.
 */

#include "line/util/line_console.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_closures.h"
#include "line/solvers/fluid/fluid_nonhyperbolic.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"
#include "line/util/svd.h"
#include "line/util/sylvester.h"

namespace line {
namespace fluid {

/**
 * The reference's `MException('LINE:FluidNonHyperbolic')`, as a type.
 *
 * IT HAS TO BE DISTINGUISHABLE FROM EVERY OTHER FAILURE. A non-hyperbolic fluid
 * fixed point cannot be detected before the mean is solved, so
 * `fluid_minnormal_applicable` cannot decline it in advance; the runner catches
 * THIS exception, and only this one, to fall back to a first-order method when
 * the moment closure was RESOLVED from `default` rather than asked for by name.
 * Catching a plain NumericError there would also swallow a singular Jacobian or a
 * failed integration, which are defects and not model properties.
 */
// FluidNonHyperbolicError now lives in fluid_nonhyperbolic.h, so that
// solver_fluid.h can raise it without including this header (which includes it).

/** What `fluid_lyapunov` reports about the fixed point it linearized at. */
struct FluidLyapunovInfo {
    std::size_t rank = 0;
    double max_real_eig = -std::numeric_limits<double>::infinity();
    bool stable = true;
};

namespace detail {

/**
 * MATLAB's `orth`: an orthonormal basis of the column space, from the SVD, over
 * the singular values above `max(size(A))*eps*sigma_1`.
 */
inline Matrix<double> fluid_orth(const Matrix<double>& A) {
    if (A.rows() == 0 || A.cols() == 0) return Matrix<double>(A.rows(), 0, 0.0);
    const SvdFactors f = svd_full(A);
    const double eps = std::numeric_limits<double>::epsilon();
    const double s1 = f.s.empty() ? 0.0 : f.s[0];
    const double tol = static_cast<double>(std::max(A.rows(), A.cols())) * eps * s1;
    std::size_t keep = 0;
    for (std::size_t i = 0; i < f.s.size(); ++i)
        if (f.s[i] > tol) ++keep;
    Matrix<double> V(A.rows(), keep, 0.0);
    for (std::size_t j = 0; j < keep; ++j)
        for (std::size_t i = 0; i < A.rows(); ++i) V(i, j) = f.U(i, j);
    return V;
}

/** Symmetrize in place: (M + M')/2. */
inline void fluid_symmetrize(Matrix<double>& M) {
    for (std::size_t i = 0; i < M.rows(); ++i)
        for (std::size_t j = i + 1; j < M.cols(); ++j) {
            const double v = 0.5 * (M(i, j) + M(j, i));
            M(i, j) = v;
            M(j, i) = v;
        }
}

}  // namespace detail

/**
 * `a + step*(b - a)` for a closure, entry by entry.
 *
 * A 0x0 covariance block is read as the zero matrix and stays 0x0 when both
 * sides are. The matrix half of the damped variance step; the twin of the
 * scalar `sigma2` blend, so the two stay consistent.
 */
inline FluidClosure fluid_blend_closure(const FluidClosure& a, const FluidClosure& b, double step) {
    FluidClosure out;
    out.sigma2.assign(b.sigma2.size(), 0.0);
    for (std::size_t i = 0; i < b.sigma2.size(); ++i) {
        const double ai = i < a.sigma2.size() ? a.sigma2[i] : 0.0;
        out.sigma2[i] = ai + step * (b.sigma2[i] - ai);
    }
    out.cov.assign(b.cov.size(), Matrix<double>(0, 0, 0.0));
    for (std::size_t i = 0; i < b.cov.size(); ++i) {
        const Matrix<double> zero(0, 0, 0.0);
        const Matrix<double>& ai = i < a.cov.size() ? a.cov[i] : zero;
        const Matrix<double>& bi = b.cov[i];
        if (ai.rows() == 0 && bi.rows() == 0) continue;
        if (ai.rows() == 0) {
            Matrix<double> m(bi.rows(), bi.cols(), 0.0);
            for (std::size_t r = 0; r < bi.rows(); ++r)
                for (std::size_t c = 0; c < bi.cols(); ++c) m(r, c) = step * bi(r, c);
            out.cov[i] = m;
        } else if (bi.rows() == 0) {
            Matrix<double> m(ai.rows(), ai.cols(), 0.0);
            for (std::size_t r = 0; r < ai.rows(); ++r)
                for (std::size_t c = 0; c < ai.cols(); ++c) m(r, c) = (1.0 - step) * ai(r, c);
            out.cov[i] = m;
        } else {
            Matrix<double> m(bi.rows(), bi.cols(), 0.0);
            for (std::size_t r = 0; r < bi.rows(); ++r)
                for (std::size_t c = 0; c < bi.cols(); ++c)
                    m(r, c) = ai(r, c) + step * (bi(r, c) - ai(r, c));
            out.cov[i] = m;
        }
    }
    return out;
}

/**
 * Port of `fluid_lyapunov.m`: the stationary covariance of the linear noise
 * approximation.
 *
 * A IS SINGULAR WHENEVER THE MODEL CONSERVES POPULATION -- every closed class
 * contributes a left null vector -- so the equation has no unique solution on the
 * full state space. It has one on the reachable subspace, which is exactly
 * range(D): the state moves only along jump directions, so the fluctuation lives
 * there and nowhere else. A = D diag(rateBase) G and Qdiff = D diag(r) D' both map
 * into range(D) as well, so restricting to an orthonormal basis of it is an EXACT
 * reduction and not an approximation, and the reduced equation is nonsingular
 * whenever the fixed point is stable.
 */
inline Matrix<double> fluid_lyapunov(const Matrix<double>& A, const Matrix<double>& Qdiff,
                                    const Matrix<double>& D, FluidLyapunovInfo& info,
                                    double tol = -1.0) {
    if (tol < 0.0) tol = std::sqrt(std::numeric_limits<double>::epsilon());
    const std::size_t n = A.rows();
    const Matrix<double> V = detail::fluid_orth(D);
    if (V.cols() == 0) {
        info.rank = 0;
        info.max_real_eig = -std::numeric_limits<double>::infinity();
        info.stable = true;
        return Matrix<double>(n, n, 0.0);
    }

    const Matrix<double> Vt = V.transpose();
    Matrix<double> Ar = matmul(matmul(Vt, A), V);
    Matrix<double> Qr = matmul(matmul(Vt, Qdiff), V);
    detail::fluid_symmetrize(Qr);

    const std::vector<std::complex<double>> ev = eig_values(Ar);
    double max_re = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < ev.size(); ++i) max_re = std::max(max_re, ev[i].real());
    info.rank = V.cols();
    info.max_real_eig = max_re;
    info.stable = max_re < -tol;
    if (!info.stable)
        throw FluidNonHyperbolicError(
            "fluid_lyapunov: the fluid fixed point is not exponentially stable on the reachable "
            "subspace (largest Jacobian eigenvalue has real part " +
            std::to_string(max_re) +
            "), so the linear noise approximation has no stationary covariance. This happens at an "
            "unstable model or at a drift kink; use method 'closing' for the mean only");

    // Ar W + W Ar' = -Qr, i.e. the Sylvester equation with B = Ar'.
    Matrix<double> negQr(Qr.rows(), Qr.cols(), 0.0);
    for (std::size_t i = 0; i < Qr.rows(); ++i)
        for (std::size_t j = 0; j < Qr.cols(); ++j) negQr(i, j) = -Qr(i, j);
    Matrix<double> W = sylvester_solve(Ar, Ar.transpose(), negQr);
    detail::fluid_symmetrize(W);
    Matrix<double> Sigma = matmul(matmul(V, W), Vt);
    detail::fluid_symmetrize(Sigma);
    return Sigma;
}

/**
 * Port of `fluid_moment_terms.m`: the event representation of the fluid
 * population process, plus the drift, rate and Jacobian handles the covariance
 * equation needs.
 *
 * The state layout, the events and the rate factors are already
 * `fluid_ode_system`'s; what this adds is the dense jump matrix, the event
 * classification the throughput is read from, the per-station and per-class index
 * blocks, and the projection that makes an OPEN model solvable.
 *
 * OPEN AND MIXED MODELS: THE COVARIANCE LIVES ON THE QUEUE COORDINATES ONLY. The
 * closing representation models a Source as an EXT pseudo-station holding unit
 * mass, so its coordinate is a normalisation constant and not a job count;
 * building D diag(r) D' over it would invent noise for a direction that carries no
 * population. Projecting those coordinates away leaves exactly the right open
 * event set, because the closing form already emits the correct events: with a
 * single-phase source the EXT rate factor is 1 - sum(of nothing) = 1 identically,
 * so an arrival is a CONSTANT-rate event whose jump, once the source row is
 * dropped, is a lone +1 into the destination queue -- the canonical exogenous
 * Poisson arrival with diffusion intensity lambda -- and the return leg LINE
 * routes Sink -> Source becomes a lone -1. The EXT row of the Jacobian is
 * identically zero for a single-phase source, so A restricted to the kept
 * coordinates IS the Jacobian of the projected drift.
 *
 * A MULTI-PHASE SOURCE IS REFUSED: those coordinates track the phase of ONE
 * arrival process, a single Markov chain rather than a population, so their
 * fluctuations are O(1) and no linear noise approximation applies to them at any
 * scale.
 */
struct FluidMomentTerms {
    FluidOdeSystem sys;
    Matrix<double> D;                                  ///< (nstate x nevents)
    std::size_t nstate = 0;
    std::vector<bool> ev_is_departure;                 ///< the leading n_departures events
    std::vector<std::size_t> ev_station, ev_class;     ///< 0-based, from the event's coordinate
    /**
     * `emap(e, o)`: expected firings of the ORIGINAL event o per firing of the
     * reduced event e; the identity when no immediate coordinate was eliminated.
     * The classification above is indexed by ORIGINAL event, so a throughput is
     * read as `r' * (emap * indicator_over_original_events)`.
     */
    Matrix<double> emap;
    /// Projector taking an initial condition onto the surviving coordinates.
    Matrix<double> absorb;
    std::vector<std::vector<std::size_t>> station_block;
    std::vector<std::vector<std::vector<std::size_t>>> class_block;
    std::vector<std::size_t> cov_idx;                  ///< coordinates carrying a real population
    std::vector<double> S;                             ///< servers, INF substituted, lld peak folded
    std::vector<bool> is_ext;
    /// stations whose occupancy cannot reach their server count, where min(n,c) is
    /// the identity and the closure must stay first order
    std::vector<bool> min_exact;
};

/**
 * True when the immediate reduction folded coordinate `s` away, so the reduced
 * drift holds no mass there and no event lands on it.
 *
 * `absorb` is the projector the reduction returns: the identity on a surviving
 * coordinate and the absorption distribution on an eliminated one, so a zero
 * diagonal is exactly the eliminated case. It is empty when nothing was
 * eliminated, where every coordinate survives.
 */
inline bool fluid_coord_eliminated(const FluidMomentTerms& t, std::size_t s) {
    if (t.absorb.rows() == 0 || s >= t.absorb.rows()) return false;
    return t.absorb(s, s) == 0.0;
}

template <class T>
FluidMomentTerms fluid_moment_terms(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    FluidMomentTerms t;
    t.sys = fluid_ode_system(sn);

    // THE MOMENT CLOSURE READS THE SAME REDUCED EVENT SET AS EVERY OTHER ROUTE.
    // It used to refuse the reduction, on the grounds that it needs the
    // untransformed event set; what it actually needs is to be able to say which
    // (station,class) each event is a completion of, and `emap` carries exactly
    // that across the composition -- an event folded through an immediate
    // coordinate keeps a row with weight on every original event it stands for,
    // including the two completions a pass-through realises at once. The
    // diffusion D*diag(r)*D' is then the diffusion of the reduced process, which
    // is the right one: the eliminated coordinate holds O(1/InfRate) mass and
    // contributes noise of the same order.
    const FluidOdeSystem sys0 = t.sys;
    if (fluid_hide_immediate(sn, opt)) {
        const FluidImmediateResult ir = fluid_eliminate_immediate(sn, t.sys);
        if (ir.eliminated) {
            t.sys = ir.sys;
            t.emap = ir.emap;
            t.absorb = ir.absorb;
        }
    }
    const FluidLayout& L = t.sys.layout;
    t.nstate = L.nstates;
    t.D = fluid_jump_matrix(t.sys);

    // A delay serves every job at once, so the reference substitutes the closed
    // population for its server count -- and floors it at one, since a pure open
    // model has no closed population and the utilization divisor would vanish.
    double npop = 0.0;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population)) npop += sn.classes[r].population;
    t.S.assign(M, 0.0);
    t.is_ext.assign(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        const double c = sn.stations[i].nservers;
        t.S[i] = std::isfinite(c) ? c : std::max(npop, 1.0);
        t.is_ext[i] = sn.stations[i].sched == lang::SchedStrategy::EXT;
    }

    // Index blocks, and the projection of the EXT source pool.
    t.station_block.assign(M, std::vector<std::size_t>());
    t.class_block.assign(M, std::vector<std::vector<std::size_t>>(K));
    std::vector<bool> keep(t.nstate, true);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < K; ++r) {
            for (std::size_t k = 0; k < L.kic[i][r]; ++k) {
                t.class_block[i][r].push_back(L.qidx[i][r] + k);
                t.station_block[i].push_back(L.qidx[i][r] + k);
            }
            if (!t.is_ext[i] || L.kic[i][r] == 0) continue;
            if (L.kic[i][r] > 1)
                throw UnsupportedError(
                    "fluid_moment_terms: the moment-closure methods need a Poisson arrival stream, "
                    "but the source of class " +
                    std::to_string(r + 1) + " is a " + std::to_string(L.kic[i][r]) +
                    "-phase process. Those coordinates track the phase of a single arrival process "
                    "rather than a population, so they carry no linear noise approximation. Use an "
                    "exponential inter-arrival time, or method 'matrix'");
            for (std::size_t k = 0; k < L.kic[i][r]; ++k) keep[L.qidx[i][r] + k] = false;
        }
    }
    for (std::size_t a = 0; a < t.nstate; ++a)
        if (keep[a]) t.cov_idx.push_back(a);

    // A STATION THAT CANNOT FILL ITS SERVERS HAS NOTHING TO CLOSE. min(n_i,c_i) is
    // the identity on the whole support whenever the occupancy of station i is
    // bounded above by its server count, and there the Gaussian closure is not an
    // improvement on the first-order one, it is an ERROR: it spreads a normal
    // marginal over n_i > c_i, mass the station can never hold, and returns
    // E[min(n_i,c_i)] < n_i. On a closed model with one job per chain the exact
    // answer is R = D at every queue (a job cannot queue behind itself), which the
    // first-order closure reproduces to machine precision while the closure reads
    // 0.4758 against 0.5 on the queue length. The bound is the total population of
    // every chain that VISITS the station -- a station may declare a service time
    // for every class while the routing never sends most of them there -- and an
    // open chain contributes an infinite population and never qualifies.
    // `solver_fluid_moments` holds the drift variance of these stations at zero,
    // exactly as it does for the delay stations, whose min() is likewise absent.
    t.min_exact.assign(M, false);
    if (!sn.chains.empty()) {
        for (std::size_t i = 0; i < M; ++i) {
            const lang::SchedStrategy s = sn.stations[i].sched;
            if (t.is_ext[i] || s == lang::SchedStrategy::INF ||
                !std::isfinite(sn.stations[i].nservers))
                continue;
            const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
            double bound = 0.0;
            bool infinite = false;
            std::vector<bool> covered(K, false);
            for (std::size_t ch = 0; ch < sn.chains.size(); ++ch) {
                bool here = false;
                double pop = 0.0;
                bool pop_inf = false;
                const bool has_vis = ch < sn.visits.size() && sn.visits[ch].rows() > isf;
                for (std::size_t r = 0; r < K; ++r) {
                    if (!sn.chains[ch][r]) continue;
                    covered[r] = true;
                    const double n = sn.classes[r].population;
                    if (std::isfinite(n)) pop += n; else pop_inf = true;
                    if (has_vis)
                        here = here || num_traits<T>::to_double(sn.visits[ch](isf, r)) > 0.0;
                    else
                        here = here || t.sys.layout.kic[i][r] > 0;
                }
                if (here) {
                    if (pop_inf) infinite = true;
                    bound += pop;
                }
            }
            bool uncovered = false;
            for (std::size_t r = 0; r < K; ++r)
                if (t.sys.layout.kic[i][r] > 0 && !covered[r]) uncovered = true;
            if (uncovered) continue;  // a class outside every chain carries no bound
            t.min_exact[i] = !infinite && bound <= sn.stations[i].nservers +
                                                       lang::GlobalConstants::FineTol;
        }
    }

    // Event classification: `ode_rate_base` emits every service completion first
    // and then every intra-PH phase change, so the leading `n_departures` events
    // are the departures. Summing their rates at (i,c) gives the class-c
    // throughput at station i EXACTLY, because the routing probabilities and the
    // entry-phase vector each sum to one over the destinations enumerated there.
    std::vector<std::size_t> coord_station(t.nstate, 0), coord_class(t.nstate, 0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            for (std::size_t k = 0; k < L.kic[i][r]; ++k) {
                coord_station[L.qidx[i][r] + k] = i;
                coord_class[L.qidx[i][r] + k] = r;
            }
    // Classified on the ORIGINAL events, which is what `emap` maps onto. Without a
    // reduction `emap` is the identity and the two indexings coincide.
    const std::size_t nev0 = sys0.events.size();
    t.ev_is_departure.assign(nev0, false);
    t.ev_station.assign(nev0, 0);
    t.ev_class.assign(nev0, 0);
    for (std::size_t e = 0; e < nev0; ++e) {
        t.ev_is_departure[e] = e < sys0.n_departures;
        t.ev_station[e] = coord_station[sys0.events[e].event_idx];
        t.ev_class[e] = coord_class[sys0.events[e].event_idx];
    }
    if (t.emap.rows() == 0) {
        t.emap = Matrix<double>(t.sys.events.size(), nev0, 0.0);
        for (std::size_t e = 0; e < t.sys.events.size() && e < nev0; ++e) t.emap(e, e) = 1.0;
    }
    return t;
}

/** The rate factors g(x) under a closure: `terms.factorFcn`. */
inline std::vector<double> fluid_moment_factors(const FluidMomentTerms& t,
                                                const std::vector<double>& x,
                                                const FluidClosure& cl) {
    FluidOdeSystem sys = t.sys;
    sys.closure = cl;
    std::vector<double> g = x;
    fluid_rates_closing_factors(sys, x.data(), g);
    return g;
}

/** The event rates r(x) under a closure: `terms.ratesFcn`. */
inline std::vector<double> fluid_moment_rates(const FluidMomentTerms& t,
                                              const std::vector<double>& x,
                                              const FluidClosure& cl) {
    const std::vector<double> g = fluid_moment_factors(t, x, cl);
    std::vector<double> r(t.sys.events.size(), 0.0);
    for (std::size_t e = 0; e < r.size(); ++e)
        r[e] = t.sys.events[e].rate_base * g[t.sys.events[e].event_idx];
    return r;
}

/** The drift F(x) = D r(x) under a closure: `terms.driftFcn`. */
inline std::vector<double> fluid_moment_drift(const FluidMomentTerms& t,
                                              const std::vector<double>& x,
                                              const FluidClosure& cl) {
    const std::vector<double> r = fluid_moment_rates(t, x, cl);
    std::vector<double> f(t.nstate, 0.0);
    for (std::size_t e = 0; e < r.size(); ++e) {
        if (r[e] == 0.0) continue;
        f[t.sys.events[e].minus] -= r[e];
        f[t.sys.events[e].plus] += r[e];
    }
    return f;
}

/**
 * Port of `fluid_drift_jacobian.m`: the analytic Jacobian of the fluid drift.
 *
 * IT MUST MIRROR `ode_rates_closing_factors` BRANCH BY BRANCH. The Jacobian
 * drives the Lyapunov equation and the 1/N refinement, so a branch that exists
 * there and not here linearizes a drift that was never integrated -- and any
 * policy with no case there keeps g = x and so contributes the identity here.
 * With sigma2 = 0 the derivative of the occupancy factor is the indicator of the
 * unsaturated region, the a.e. derivative of the first-order closure; with
 * sigma2 > 0 it is the smooth derivative the closure returns.
 */
inline Matrix<double> fluid_drift_jacobian(const FluidMomentTerms& t, const std::vector<double>& x,
                                           const FluidClosure& cl) {
    const FluidOdeSystem& sys = t.sys;
    const FluidLayout& L = sys.layout;
    const std::size_t M = L.qidx.size();
    const std::size_t K = M ? L.qidx[0].size() : 0;
    const std::size_t n = t.nstate;
    const bool gaussian = cl.gaussian();

    Matrix<double> G = eye<double>(n);  // INF, EXT phases 2.., and every policy without a case
    const auto zero_rows = [&](const std::vector<std::size_t>& rows) {
        for (std::size_t a = 0; a < rows.size(); ++a)
            for (std::size_t j = 0; j < n; ++j) G(rows[a], j) = 0.0;
    };

    for (std::size_t i = 0; i < M; ++i) {
        const std::vector<double>& lld = sys.lld[i];
        const double s2 = cl.sigma2_of(i);
        const Matrix<double>* Ci = cl.cov_of(i);
        const std::vector<std::size_t>& blk = t.station_block[i];
        const std::size_t nb = blk.size();
        if (nb == 0) continue;
        double ni = 0.0;
        for (std::size_t a = 0; a < nb; ++a) ni += x[blk[a]];

        switch (sys.sched[i]) {
            case lang::SchedStrategy::INF: {
                if (lld.empty() || !(ni > 0.0)) break;
                const ClosureValue h = fluid_capacity_closure(ni, t.S[i], s2, lld, true);
                const double f = h.h / ni, fp = (ni * h.dh - h.h) / (ni * ni);
                zero_rows(blk);
                for (std::size_t a = 0; a < nb; ++a)
                    for (std::size_t b = 0; b < nb; ++b)
                        G(blk[a], blk[b]) = (a == b ? f : 0.0) + x[blk[a]] * fp;
                break;
            }
            case lang::SchedStrategy::EXT: {
                for (std::size_t r = 0; r < K; ++r) {
                    if (!L.enabled[i][r]) continue;
                    const std::size_t b = L.qidx[i][r], nn = L.kic[i][r];
                    for (std::size_t j = 0; j < n; ++j) G(b, j) = 0.0;
                    for (std::size_t p = 1; p < nn; ++p) G(b, b + p) = -1.0;
                }
                break;
            }
            case lang::SchedStrategy::PS:
            case lang::SchedStrategy::FCFS: {
                if (!(ni > 0.0)) break;  // g = x on an empty station
                double h = 0.0, dh = 0.0;
                if (gaussian || !lld.empty()) {
                    const ClosureValue cv = fluid_capacity_closure(ni, t.S[i], s2, lld, false);
                    h = cv.h;
                    dh = cv.dh;
                    if (Ci != nullptr) {
                        // g = s(x_blk)*h(ni) + h'(ni)*cn(x_blk), the joint closure
                        // of the share and the capacity; reached only from this
                        // branch, exactly as in the drift. Differentiating it with
                        // C held fixed adds h'*dcn and h''*cn to the product rule.
                        const double d2h = fluid_capacity_closure(ni, t.S[i], s2, lld, false).d2h;
                        std::vector<double> xb(nb, 0.0), wv(nb, 1.0);
                        for (std::size_t a = 0; a < nb; ++a) xb[a] = x[blk[a]];
                        const ShareValue sh = fluid_share_closure(xb, wv, *Ci, true, true);
                        zero_rows(blk);
                        for (std::size_t a = 0; a < nb; ++a)
                            for (std::size_t b = 0; b < nb; ++b)
                                G(blk[a], blk[b]) = sh.ds(a, b) * h + sh.s[a] * dh +
                                                    dh * sh.dcn(a, b) + sh.cn[a] * d2h;
                        break;
                    }
                } else if (ni > t.S[i] - lang::GlobalConstants::FineTol * std::max(1.0, ni)) {
                    // THE SATURATION TEST CARRIES A BAND, and it is a cross-codebase
                    // requirement: a saturated fixed point sits exactly at ni = c, and
                    // each engine's ODE stops on its own residual (MATLAB 1.0004, this
                    // port 1 - 1.8e-13 on the same model). A strict ni > c reads
                    // saturated in one and unsaturated in the other, which flips this
                    // whole station block between a zero row and the identity, and with
                    // it the hyperbolicity verdict `fluid_lyapunov` returns and the
                    // method SolverFluid ends up answering with. See
                    // `fluid_min_closure`, whose degenerate branch carries the band too.
                    h = t.S[i];
                    dh = 0.0;
                } else {
                    break;  // g = x, the identity is already in place
                }
                // g_j = x_j h(ni)/ni -> dg_j/dx_m = delta_jm f + x_j f',
                // f = h/ni, f' = (ni dh - h)/ni^2
                const double f = h / ni, fp = (ni * dh - h) / (ni * ni);
                zero_rows(blk);
                for (std::size_t a = 0; a < nb; ++a)
                    for (std::size_t b = 0; b < nb; ++b)
                        G(blk[a], blk[b]) = (a == b ? f : 0.0) + x[blk[a]] * fp;
                break;
            }
            case lang::SchedStrategy::DPS: {
                double wsum = 0.0;
                for (std::size_t r = 0; r < K; ++r) wsum += sys.weight[i][r];
                if (wsum <= 0.0) break;
                std::vector<double> wv(nb, 0.0), xb(nb, 0.0);
                for (std::size_t r = 0; r < K; ++r) {
                    if (!L.enabled[i][r]) continue;
                    const std::size_t b = L.qidx[i][r] - blk[0];
                    for (std::size_t p = 0; p < L.kic[i][r]; ++p)
                        wv[b + p] = sys.weight[i][r] / wsum;
                }
                double wx = 0.0;
                for (std::size_t a = 0; a < nb; ++a) {
                    xb[a] = x[blk[a]];
                    wx += wv[a] * xb[a];
                }
                if (!(ni > 0.0) || !(wx > 0.0)) break;  // g = x on an empty station
                const ClosureValue psi = fluid_capacity_closure(ni, t.S[i], s2, lld, false);
                const ShareValue sh = fluid_share_closure(
                    xb, wv, Ci ? *Ci : Matrix<double>(0, 0, 0.0), true, true);
                zero_rows(blk);
                for (std::size_t a = 0; a < nb; ++a)
                    for (std::size_t b = 0; b < nb; ++b)
                        G(blk[a], blk[b]) = sh.ds(a, b) * psi.h + sh.s[a] * psi.dh +
                                            psi.dh * sh.dcn(a, b) + sh.cn[a] * psi.d2h;
                break;
            }
            case lang::SchedStrategy::GPS: {
                if (sys.nservers[i] > 1.0)
                    throw UnsupportedError(
                        "fluid_drift_jacobian: multi-server GPS stations are not supported, as in "
                        "the reference");
                std::vector<double> xk(K, 0.0), vk(K, 0.0), wk(K, 0.0);
                for (std::size_t r = 0; r < K; ++r) {
                    wk[r] = sys.weight[i][r];
                    const std::vector<std::size_t>& bk = t.class_block[i][r];
                    for (std::size_t a = 0; a < bk.size(); ++a) xk[r] += x[bk[a]];
                    if (Ci == nullptr) continue;
                    double v = 0.0;
                    for (std::size_t a = 0; a < bk.size(); ++a)
                        for (std::size_t b = 0; b < bk.size(); ++b)
                            v += (*Ci)(bk[a] - blk[0], bk[b] - blk[0]);
                    vk[r] = std::max(0.0, v);
                }
                const ShareValue sk = fluid_gps_share(xk, wk, vk, true);
                ClosureValue a1;
                a1.h = 1.0;
                a1.dh = 0.0;
                if (!lld.empty()) a1 = fluid_lld_scaling(lld, ni);
                zero_rows(blk);
                // g_j = (x_j/x_k) s_k a for coordinate j of class k, so
                //   dg_j/dx_l = [delta_jl/x_k - x_j/x_k^2] s_k a       (l in class k)
                //             + (x_j/x_k) ds_k/dx_m a                 (l in class m)
                //             + (x_j/x_k) s_k da/dxi                  (l in the station)
                for (std::size_t r = 0; r < K; ++r) {
                    const std::vector<std::size_t>& bk = t.class_block[i][r];
                    if (bk.empty() || !(xk[r] > 0.0)) continue;
                    for (std::size_t a = 0; a < bk.size(); ++a) {
                        for (std::size_t b = 0; b < bk.size(); ++b)
                            G(bk[a], bk[b]) += (a == b ? sk.s[r] * a1.h / xk[r] : 0.0) -
                                               (sk.s[r] * a1.h / (xk[r] * xk[r])) * x[bk[a]];
                        for (std::size_t m = 0; m < K; ++m) {
                            const std::vector<std::size_t>& bm = t.class_block[i][m];
                            for (std::size_t b = 0; b < bm.size(); ++b)
                                G(bk[a], bm[b]) += (a1.h * sk.ds(r, m) / xk[r]) * x[bk[a]];
                        }
                        if (a1.dh != 0.0)
                            for (std::size_t b = 0; b < nb; ++b)
                                G(bk[a], blk[b]) += (sk.s[r] * a1.dh / xk[r]) * x[bk[a]];
                    }
                }
                break;
            }
            default:
                break;
        }
    }

    // A = D (rateBase .* G(eventIdx,:)), assembled through the two-index event
    // form: every column of D has one -1 and one +1, so this is the same product
    // without building the (nstate x nevents) intermediate.
    Matrix<double> A(n, n, 0.0);
    for (std::size_t e = 0; e < sys.events.size(); ++e) {
        const FluidEvent& ev = sys.events[e];
        if (ev.rate_base == 0.0) continue;
        for (std::size_t j = 0; j < n; ++j) {
            const double v = ev.rate_base * G(ev.event_idx, j);
            if (v == 0.0) continue;
            A(ev.minus, j) -= v;
            A(ev.plus, j) += v;
        }
    }
    return A;
}

/**
 * Every station whose population sits ON the saturation kink n_i = c_i of the
 * first-order rate factor, in increasing order, empty when none does.
 *
 * With sigma2 = 0 the occupancy factor is min(n_i, c_i), whose derivative is the
 * indicator of the unsaturated region: slope 1 below c_i, slope 0 above, and NO
 * derivative at c_i itself. `fluid_drift_jacobian` resolves the tie onto the
 * saturated side, as its MATLAB, python and JAR twins do, so it silently returns
 * one one-sided value there; which side a fixed point lands on is decided by the
 * integrator's rounding residue rather than by the model. Callers that need a
 * differentiable drift consult this instead of trusting the tie-break.
 *
 * Only the branches that take the indicator derivative can sit on a kink: a
 * positive sigma2 or a load-dependent row makes the closure smooth, and an
 * infinite server never saturates. Twin of
 * `FluidRateFactors.driftKinkStations` in the JAR.
 */
inline std::vector<std::size_t> fluid_kink_stations(const FluidMomentTerms& t,
                                                    const std::vector<double>& x,
                                                    const FluidClosure& cl) {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < cl.sigma2.size(); ++i)
        if (cl.sigma2[i] > 0.0) return out;  // the Gaussian closure has no kink
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon());
    for (std::size_t i = 0; i < t.station_block.size(); ++i) {
        const lang::SchedStrategy sc = t.sys.sched[i];
        if (sc == lang::SchedStrategy::INF || sc == lang::SchedStrategy::EXT) continue;
        if (!t.sys.lld[i].empty()) continue;  // psi is piecewise quadratic and smooth
        const double c = t.S[i];
        if (!std::isfinite(c) || c <= 0.0) continue;
        const std::vector<std::size_t>& blk = t.station_block[i];
        if (blk.empty()) continue;
        double ni = 0.0;
        for (std::size_t a = 0; a < blk.size(); ++a) ni += x[blk[a]];
        if (!(ni > 0.0)) continue;  // g = x on an empty station
        if (std::fabs(ni - c) <= tol * std::max(1.0, c)) out.push_back(i);
    }
    return out;
}

/**
 * A copy of `x` with every station in `kink` moved to `c_i*(1 + rel)`, i.e.
 * strictly onto one side of its kink. The station's coordinates are scaled
 * together, so the phase mix and every other station are untouched.
 */
inline std::vector<double> fluid_nudge_off_kink(const FluidMomentTerms& t,
                                                const std::vector<double>& x,
                                                const std::vector<std::size_t>& kink, double rel) {
    std::vector<double> y = x;
    for (std::size_t k = 0; k < kink.size(); ++k) {
        const std::vector<std::size_t>& blk = t.station_block[kink[k]];
        double ni = 0.0;
        for (std::size_t a = 0; a < blk.size(); ++a) ni += y[blk[a]];
        if (!(ni > 0.0)) continue;
        const double scale = t.S[kink[k]] * (1.0 + rel) / ni;
        for (std::size_t a = 0; a < blk.size(); ++a) y[blk[a]] *= scale;
    }
    return y;
}

/**
 * `local_lyapunov` of `solver_fluid_moments.m`: the covariance on the coordinates
 * that carry a real population, scattered back to full size.
 *
 * For a closed model `cov_idx` is every coordinate and this is the plain solve.
 * For an open or mixed model it drops the EXT source pool; the zeros left on the
 * dropped rows keep the station and class block indexing unchanged downstream.
 */
inline Matrix<double> fluid_moment_lyapunov(const FluidMomentTerms& t,
                                            const Matrix<double>& A,
                                            const std::vector<double>& r,
                                            const Matrix<double>* clampT = nullptr) {
    const std::size_t nc = t.cov_idx.size(), nev = r.size();
    Matrix<double> Dc(nc, nev, 0.0);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t e = 0; e < nev; ++e) Dc(a, e) = t.D(t.cov_idx[a], e);
    // CLAMPT, when given, is the tangent space of the caps that CLAMP: a cap that
    // holds the job upstream or loses it fixes its own combination of the state
    // while it binds, so that combination does not fluctuate. Projecting the jump
    // directions is enough to state the reduced problem, because FLUID_LYAPUNOV
    // restricts everything to range(D) already. See the DAE route's clamp tangent.
    if (clampT && clampT->rows() == nc && clampT->cols() == nc) {
        Matrix<double> Dp(nc, nev, 0.0);
        for (std::size_t a = 0; a < nc; ++a)
            for (std::size_t e = 0; e < nev; ++e) {
                double acc = 0.0;
                for (std::size_t b = 0; b < nc; ++b) acc += (*clampT)(a, b) * Dc(b, e);
                Dp(a, e) = acc;
            }
        Dc = Dp;
    }
    Matrix<double> Qc(nc, nc, 0.0);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t b = 0; b < nc; ++b) {
            double acc = 0.0;
            for (std::size_t e = 0; e < nev; ++e) acc += Dc(a, e) * r[e] * Dc(b, e);
            Qc(a, b) = acc;
        }
    Matrix<double> Ac(nc, nc, 0.0);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t b = 0; b < nc; ++b) Ac(a, b) = A(t.cov_idx[a], t.cov_idx[b]);

    FluidLyapunovInfo info;
    const Matrix<double> Sc = fluid_lyapunov(Ac, Qc, Dc, info);
    Matrix<double> Sigma(t.nstate, t.nstate, 0.0);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t b = 0; b < nc; ++b) Sigma(t.cov_idx[a], t.cov_idx[b]) = Sc(a, b);
    return Sigma;
}

/** What `fluid_refine_meanfield` reports about the correction it computed. */
struct FluidRefineInfo {
    std::size_t rank = 0;
    double stepsize = 0.0;
    double residual = 0.0;
    double condition = 0.0;
};

/**
 * Port of `fluid_refine_meanfield.m`: the O(1/N) refined mean field correction of
 * Gast (POMACS 2017).
 *
 * The correction V solves A V + (1/2) sum_{jk} Sigma_jk d2F/dx_j dx_k = 0. The
 * Hessian contraction is evaluated WITHOUT forming the tensor: with
 * Sigma = sum_m lam_m v_m v_m', the contraction is sum_m lam_m d2F/dv_m^2 and each
 * directional second derivative is one central second difference, so the cost is
 * O(rank(Sigma)) drift evaluations rather than O(n^2). Because Sigma scales with
 * the population, V is the O(1/N) term written directly in job counts and no
 * explicit density rescaling is needed.
 *
 * THE DRIFT MUST BE TWICE DIFFERENTIABLE. The first-order closure is only
 * piecewise linear -- second derivative zero away from the kink and a delta at it
 * -- so a zero variance is REFUSED rather than silently returning a null
 * correction.
 */
inline std::vector<double> fluid_refine_meanfield(const FluidMomentTerms& t,
                                                  const std::vector<double>& x,
                                                  const FluidClosure& cl,
                                                  const Matrix<double>& Sigma,
                                                  FluidRefineInfo& info, double epsrel = 1e-4) {
    if (!cl.gaussian())
        throw InputError(
            "fluid_refine_meanfield: the refined mean field expansion needs a twice-differentiable "
            "drift, but the first-order closure is only piecewise linear. Reach this function "
            "through method 'refined', which converges the Gaussian closure first");

    const std::size_t n = x.size();
    Matrix<double> Sig = Sigma;
    detail::fluid_symmetrize(Sig);
    // Sigma is symmetric positive semidefinite, so its SVD IS its
    // eigendecomposition: the singular values are the eigenvalues and the left
    // singular vectors the eigenvectors. Using it avoids a second, symmetric
    // eigensolver for a matrix that already has one.
    const SvdFactors f = svd_full(Sig);
    const double eps = std::numeric_limits<double>::epsilon();
    const double lmax = f.s.empty() ? 0.0 : f.s[0];
    std::vector<std::size_t> keep;
    for (std::size_t m = 0; m < f.s.size(); ++m)
        if (f.s[m] > lmax * std::sqrt(eps) && f.s[m] > 0.0) keep.push_back(m);

    double xnorm = 0.0;
    for (std::size_t a = 0; a < n; ++a) xnorm += x[a] * x[a];
    xnorm = std::sqrt(xnorm);
    const double scale = std::max(1.0, xnorm);
    const double step = epsrel * scale;

    const std::vector<double> F0 = fluid_moment_drift(t, x, cl);
    std::vector<double> b(n, 0.0);
    for (std::size_t idx = 0; idx < keep.size(); ++idx) {
        const std::size_t m = keep[idx];
        std::vector<double> xp = x, xm = x;
        for (std::size_t a = 0; a < n; ++a) {
            xp[a] += step * f.U(a, m);
            xm[a] -= step * f.U(a, m);
        }
        const std::vector<double> Fp = fluid_moment_drift(t, xp, cl);
        const std::vector<double> Fm = fluid_moment_drift(t, xm, cl);
        for (std::size_t a = 0; a < n; ++a)
            b[a] += f.s[m] * (Fp[a] - 2.0 * F0[a] + Fm[a]) / (step * step);
    }
    for (std::size_t a = 0; a < n; ++a) b[a] *= 0.5;

    // Solve A V = -b on the reachable subspace, where A is invertible.
    const Matrix<double> A = fluid_drift_jacobian(t, x, cl);
    const Matrix<double> Vb = detail::fluid_orth(t.D);
    const Matrix<double> Vbt = Vb.transpose();
    const Matrix<double> Ar = matmul(matmul(Vbt, A), Vb);
    const std::vector<double> sv = svd_values(Ar);
    const double cond = (sv.empty() || sv.back() == 0.0)
                            ? std::numeric_limits<double>::infinity()
                            : sv.front() / sv.back();
    if (!std::isfinite(cond) || cond > 1.0 / std::sqrt(eps))
        throw FluidNonHyperbolicError(
            "fluid_refine_meanfield: the fluid Jacobian is numerically singular on the reachable "
            "subspace (condition number " +
            std::to_string(cond) +
            "), so the refinement equation A V = -b has no meaningful solution. The fixed point "
            "sits at a drift kink or the model is marginally stable; use method 'minnormal', which "
            "resums the same correction without inverting A");

    std::vector<double> rhs(Vb.cols(), 0.0);
    for (std::size_t j = 0; j < Vb.cols(); ++j) {
        double acc = 0.0;
        for (std::size_t a = 0; a < n; ++a) acc += Vb(a, j) * b[a];
        rhs[j] = -acc;
    }
    const Matrix<double> Arinv = inverse(Ar);
    std::vector<double> vr(Vb.cols(), 0.0);
    for (std::size_t j = 0; j < Vb.cols(); ++j) {
        double acc = 0.0;
        for (std::size_t k = 0; k < Vb.cols(); ++k) acc += Arinv(j, k) * rhs[k];
        vr[j] = acc;
    }
    std::vector<double> V(n, 0.0);
    for (std::size_t a = 0; a < n; ++a) {
        double acc = 0.0;
        for (std::size_t j = 0; j < Vb.cols(); ++j) acc += Vb(a, j) * vr[j];
        V[a] = acc;
    }

    // The refinement is the next term of an asymptotic expansion, so it is only
    // meaningful while it stays small against the leading term; a correction the
    // size of the fixed point means the expansion has not kicked in at this
    // population, and returning it would be worse than refusing.
    double vnorm = 0.0, resid = 0.0;
    for (std::size_t a = 0; a < n; ++a) vnorm += V[a] * V[a];
    vnorm = std::sqrt(vnorm);
    for (std::size_t a = 0; a < n; ++a) {
        double acc = b[a];
        for (std::size_t j = 0; j < n; ++j) acc += A(a, j) * V[j];
        resid += acc * acc;
    }
    info.rank = keep.size();
    info.stepsize = step;
    info.residual = std::sqrt(resid);
    info.condition = cond;
    if (vnorm > 0.5 * std::max(xnorm, std::sqrt(eps)))
        throw FluidNonHyperbolicError(
            "fluid_refine_meanfield: the 1/N refinement (norm " + std::to_string(vnorm) +
            ") is not small against the mean-field fixed point (norm " + std::to_string(xnorm) +
            "), so the asymptotic expansion is outside its range of validity at this population. "
            "Use method 'minnormal'");
    return V;
}

/**
 * Port of `solver_fluid_moments.m`: the second-order fluid analysis backing
 * `minnormal` and `refined`.
 *
 * THE OUTER FIXED POINT IS OVER THE VARIANCE, not over the mean. Each sweep
 * solves the mean at the current closure variance -- through the ordinary closing
 * integration, which is why the closure travels on `FluidOptions` -- then solves
 * the Lyapunov equation at that mean and reads a new variance off the covariance
 * blocks. `sigma2` alone is not enough to iterate on: the DPS and PS capacity
 * share is a RATIO of coordinates, so closing it needs the covariance BETWEEN
 * them, and the blocks are carried through the same fixed point and compared in
 * the same convergence test -- a block sum can converge while the off-diagonals
 * the share closure reads are still moving.
 *
 * THE METRICS ARE READ AT THE VARIANCE THE MEAN SOLVE USED, not at the variance
 * that solve produced. Using the latter evaluates the rate functions away from
 * their own fixed point and throughput stops balancing: with the variance held at
 * zero for the mean solve, Tput came back 2.000000 at the delay against 1.949745
 * at the queue on Delay -> Queue(PS,c=2), N=6, a 2.5% gap in a closed cycle where
 * the two must be equal. The two differ only within the outer tolerance once the
 * fixed point has converged.
 *
 * A DELAY'S VARIANCE IS KEPT FOR REPORTING AND EXCLUDED FROM THE DRIFT: there is
 * no min() to close at an infinite server, so letting it in would perturb a term
 * that is exactly linear.
 */
template <class T>
FluidSolution solver_fluid_moments(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "solver_fluid_moments: the fluid solver integrates its drift with LSODA, whose "
            "coefficients assume double precision; rerun with --arith double");

    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::string m = opt.method;
    if (m.size() > 6 && m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    if (!(m == "minnormal" || m == "refined"))
        throw UnsupportedError("solver_fluid_moments: '" + opt.method +
                               "' is not a moment-closure method; only 'minnormal' and 'refined' "
                               "are solved here");
    // `fluid_moment_terms.m:114`. The closure solves a STATIONARY Lyapunov
    // equation, so it needs an autonomous drift; a time-varying rate multiplier
    // leaves no fixed point for a stationary covariance to sit at.
    if (detail::fluid_has_time_varying_rates(opt))
        throw UnsupportedError(
            "solver_fluid_moments: the moment closures require an AUTONOMOUS drift, but "
            "options.config.rate_traj / nhpp_sched / rate_sched make the rates time-varying. Use "
            "method 'closing' or 'matrix'");
    // The EXT projection covers `minnormal` only. `fluid_refine_meanfield` solves
    // its correction on orth(D) over the FULL state and would add a perturbation to
    // the source pool mass, which is a normalisation constant rather than a
    // population; only minnormal was validated open, so refined keeps the
    // closed-model restriction instead of being declared on an untested path.
    if (m == "refined")
        for (std::size_t r = 0; r < K; ++r)
            if (!std::isfinite(sn.classes[r].population))
                throw UnsupportedError(
                    "solver_fluid_moments: the 'refined' method supports closed models only: its "
                    "1/N correction is solved over the full state, including the source pool. Use "
                    "method 'minnormal' for open or mixed models");

    const FluidMomentTerms terms = fluid_moment_terms(sn, opt);

    // The covariance is a dense nstate-by-nstate object and the Lyapunov solve is
    // cubic in it, so refuse rather than silently crawl. The same cap decides
    // whether `default` resolves here at all (`fluid_minnormal_applicable`).
    const std::size_t maxstate = opt.moment_maxstate;
    if (terms.nstate > maxstate)
        throw UnsupportedError(
            "solver_fluid_moments: the moment-closure methods solve a " +
            std::to_string(terms.nstate) + "x" + std::to_string(terms.nstate) +
            " Lyapunov equation, above the limit of " + std::to_string(maxstate) +
            " set by moment_maxstate. Raise that limit or use method 'closing'");

    // A station whose share is a ratio needs its covariance BLOCK, not only the
    // block sum: PS, FCFS, DPS and GPS all read one.
    std::vector<bool> share_sched(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        const lang::SchedStrategy s = terms.sys.sched[i];
        share_sched[i] = s == lang::SchedStrategy::PS || s == lang::SchedStrategy::FCFS ||
                         s == lang::SchedStrategy::DPS || s == lang::SchedStrategy::GPS;
    }

    std::size_t outer_max = 20;
    if (opt.iter_max > 0) outer_max = std::min<std::size_t>(outer_max, std::max<std::size_t>(2, opt.iter_max));
    // THE CLOSURE IS JUDGED FAR TIGHTER THAN CoarseTol, so it must not stop there.
    // Converged only to 1e-3 this alternation is not a fixed point to two machines:
    // on mqn_singleserver_ps the MATLAB twin answered 42.3962 on two hosts and
    // 42.8207 on a third, 1e-2 relative apart, because the transient iterate below
    // fell on opposite sides. 1e-6 is the loosest that reproduces; min(), not
    // assignment, so a caller may still ask tighter.
    //
    // The INNER mean solve is deliberately left alone here, unlike in the MATLAB
    // and Java twins. FluidOptions::iter_tol carries the OPPOSITE sense in this
    // codebase -- 0, the default, runs to iter_max and is the TIGHTEST setting,
    // while a positive value stops early -- so handing it mom_tol would loosen the
    // solve the other ports tighten.
    double mom_tol = 1e-6;
    if (opt.iter_tol > 0.0) mom_tol = std::min(mom_tol, opt.iter_tol);
    const double outer_tol = mom_tol;

    FluidClosure cl;  // the variance the NEXT mean solve will use
    cl.sigma2.assign(M, 0.0);
    cl.cov.assign(M, Matrix<double>(0, 0, 0.0));
    FluidClosure used = cl;  // the variance the LAST mean solve actually used
    // The last closure whose Lyapunov solve SUCCEEDED, and the floor on the step
    // taken toward the next one. See the damping in the loop below.
    FluidClosure okcl;
    okcl.sigma2.assign(M, 0.0);
    okcl.cov.assign(M, Matrix<double>(0, 0, 0.0));
    const double damp_min = 1.0 / 64.0;
    Matrix<double> Sigma(terms.nstate, terms.nstate, 0.0);
    std::vector<double> x;
    std::size_t iters = 0, outer = 0;

    // A DELAY STATION HAS NO min() TO CLOSE, so its variance must never reach the
    // drift -- only the report. This mask used to be applied to `drift_cl` after the
    // loop and nowhere inside it, so every mean solve of the fixed point ran with the
    // delay variance switched on. The rate factor there is mu*n, which the Gaussian
    // correction turns into something that does not vanish with n: the coordinate is
    // driven NEGATIVE, the drift is conservative so another coordinate grows to match,
    // and the trajectory leaves the simplex for good. On CQN_Cox_CS_9 (Delay + PS +
    // PS(c=5), N=6) the first window past sigma2 = 0 moved 8.7e3 of mass and the drift
    // norm reached 5.4e9.
    std::vector<bool> no_drift_var(M, false);
    for (std::size_t i = 0; i < M; ++i)
        no_drift_var[i] = terms.sys.sched[i] == lang::SchedStrategy::INF ||
                          terms.sys.sched[i] == lang::SchedStrategy::EXT;

    line::util::LineConsole::loop("iterating the moment closure (at most %zu passes)", outer_max);
    for (outer = 1; outer <= outer_max; ++outer) {
        line::util::LineConsole::iter(static_cast<long>(outer),
                                      "closure pass %zu: %zu ODE iterations so far", outer, iters);
        // A TRANSIENT ITERATE MUST NOT VETO THE METHOD. The Lyapunov gate asks
        // whether the linear noise approximation has a stationary covariance at the
        // point THIS iterate landed on; a fixed point that fails it is a model the
        // closure cannot answer, but an intermediate iterate that fails it is only a
        // variance step that overshot. On mqn_singleserver_ps iterate 1 was stable at
        // -4.93e-03, iterate 2 declined at +1.07e+01, and the fixed point the fallback
        // then found was stable at -4.92e-03. So a failing iterate RETREATS toward the
        // last closure that succeeded, halving until the LNA is defined again; only a
        // step below damp_min, or a failure at the seed where there is nothing to
        // retreat toward, is the model's own non-hyperbolicity and still throws.
        double step = 1.0;
        FluidClosure trycl;
        for (;;) {
            trycl = fluid_blend_closure(okcl, cl, step);
            used = trycl;
            // A station whose occupancy cannot reach its server count has min(n,c) = n on
            // the whole support, so the closure there must stay first order: see
            // `fluid_moment_terms`, which decides it from the chain populations. Its share
            // closure follows, because mu_r*(n_r/n)*min(n,c) collapses to mu_r*n_r once the
            // min is the identity. The covariance is still solved for these stations and
            // still reported, it just does not enter the drift, exactly as at the delay
            // stations below.
            for (std::size_t i = 0; i < M; ++i) {
                if (!terms.min_exact[i] && !no_drift_var[i]) continue;
                if (i < used.sigma2.size()) used.sigma2[i] = 0.0;
                if (i < used.cov.size()) used.cov[i] = Matrix<double>(0, 0, 0.0);
            }
            FluidOptions mo = opt;
            mo.method = "closing";  // the closure enters through the drift, not the name
            mo.closure = used;
            const FluidSolution mean = detail::fluid_dispatch(sn, mo);
            iters += mean.iters;
            x = mean.xvec;

            const std::vector<double> r = fluid_moment_rates(terms, x, used);

            // A POINT ON A SATURATION KINK HAS NO JACOBIAN. `fluid_drift_jacobian`
            // resolves the tie onto the saturated side, so a verdict read off it would
            // depend on which side the integrator stopped. The VERDICT, not the point,
            // has to be side-independent: both one-sided Jacobians are ordinary
            // matrices, so ASK BOTH and decline only when a side fails. Refusing at
            // every kink instead throws away models the reference solves -- the first
            // outer iterate runs at sigma2 = 0 and a saturated model's first-order
            // fixed point lands on the kink by construction. Later iterates carry a
            // positive sigma2 and are smooth, so this costs two Jacobians on the seed
            // and nothing after it. Twin of
            // `FluidRateFactors.driftKinkStation`/`nudgedOffKink` in the JAR.
            const std::vector<std::size_t> kink = fluid_kink_stations(terms, x, used);
            if (!kink.empty()) {
                const double probe[2] = {-1e-6, 1e-6};
                for (std::size_t side = 0; side < 2; ++side) {
                    const std::vector<double> xs = fluid_nudge_off_kink(terms, x, kink, probe[side]);
                    try {
                        fluid_moment_lyapunov(terms, fluid_drift_jacobian(terms, xs, used),
                                              fluid_moment_rates(terms, xs, used));
                    } catch (const FluidNonHyperbolicError& e) {
                        throw FluidNonHyperbolicError(
                            "solver_fluid_moments: the fluid fixed point sits on the saturation kink "
                            "of station " + std::to_string(kink[0] + 1) + " (population equals its " +
                            std::to_string(terms.S[kink[0]]) +
                            " servers) and the two one-sided drift Jacobians there disagree on "
                            "hyperbolicity, so which of them the linear noise approximation would use "
                            "is decided by the integrator's rounding residue rather than by the model. "
                            "This is the saturated boundary of a continuum of equilibria; use method "
                            "'closing' for the mean only. Underlying: " + std::string(e.what()));
                    }
                }
            }

            const Matrix<double> A = fluid_drift_jacobian(terms, x, used);
            try {
                Sigma = fluid_moment_lyapunov(terms, A, r);
                break;
            } catch (const FluidNonHyperbolicError&) {
                bool at_seed = true;
                for (std::size_t i = 0; i < M && at_seed; ++i)
                    if (cl.sigma2[i] != okcl.sigma2[i]) at_seed = false;
                for (std::size_t i = 0; i < M && at_seed; ++i)
                    if (cl.cov[i].rows() != 0) at_seed = false;
                if (step <= damp_min || at_seed) throw;
                step = step / 2.0;
            }
        }
        okcl = trycl;

        std::vector<double> s2new(M, 0.0);
        std::vector<Matrix<double>> covnew(M, Matrix<double>(0, 0, 0.0));
        for (std::size_t i = 0; i < M; ++i) {
            const std::vector<std::size_t>& blk = terms.station_block[i];
            if (blk.empty()) continue;
            double acc = 0.0;
            for (std::size_t a = 0; a < blk.size(); ++a)
                for (std::size_t b = 0; b < blk.size(); ++b) acc += Sigma(blk[a], blk[b]);
            s2new[i] = std::max(0.0, acc);
            if (!share_sched[i]) continue;
            Matrix<double> B(blk.size(), blk.size(), 0.0);
            for (std::size_t a = 0; a < blk.size(); ++a)
                for (std::size_t b = 0; b < blk.size(); ++b) B(a, b) = Sigma(blk[a], blk[b]);
            covnew[i] = B;
        }

        double l1new = 0.0, l1diff = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            l1new += std::fabs(s2new[i]);
            l1diff += std::fabs(s2new[i] - trycl.sigma2[i]);
        }
        double delta = l1diff / std::max(1.0, l1new);
        // `norm(dc,1)` on a MATRIX is the maximum absolute COLUMN SUM, not the
        // entrywise sum that the same call gives on a vector. Summing every entry
        // instead overstates the residual, so the loop ran past the reference's
        // break and settled on a different closure fixed point: on the 3-station
        // 2-class PS model of the parity corpus that was Tput 0.8232419 against
        // 0.8234276, a 2.3e-4 gap that no tolerance change could close because
        // both sides were converged, just to different points.
        for (std::size_t i = 0; i < M; ++i) {
            if (covnew[i].rows() == 0) continue;
            double dn = 0.0, nn = 0.0;
            for (std::size_t b = 0; b < covnew[i].cols(); ++b) {
                double dcol = 0.0, ncol = 0.0;
                for (std::size_t a = 0; a < covnew[i].rows(); ++a) {
                    const double old = (trycl.cov[i].rows() == covnew[i].rows()) ? trycl.cov[i](a, b) : 0.0;
                    dcol += std::fabs(covnew[i](a, b) - old);
                    ncol += std::fabs(covnew[i](a, b));
                }
                dn = std::max(dn, dcol);
                nn = std::max(nn, ncol);
            }
            delta = std::max(delta, dn / std::max(1.0, nn));
        }
        cl.sigma2 = s2new;
        cl.cov = covnew;
        if (delta < outer_tol) break;
    }
    const std::size_t outer_iters = std::min(outer, outer_max);

    // The drift closure: the variance the mean solve used, with the delays already
    // excluded on the way in by `no_drift_var`, because they have no min() to close.
    //
    // NOT const, and the reason is the `refined` branch below: it re-points this at
    // the MEAN-FIELD closure once it has corrected the base point, and that choice is
    // read after the branch by `gfac`. So the variable carries the drift closure of
    // whichever method ran -- converged for `minnormal`, zero for `refined` -- and
    // making it const compiles only if that distinction is dropped, which would read
    // `refined`'s rate factors at a variance its correction has already resummed.
    FluidClosure drift_cl = used;

    std::vector<double> refinement;
    std::vector<double> r;
    if (m == "refined") {
        // The refinement is a truncated expansion about the MEAN-FIELD fixed
        // point, not about the Gaussian one: adding it to the `minnormal` point
        // would count the same O(1/N) term twice, since the Gaussian closure
        // already resums it. So the base point is recomputed with the first-order
        // closure while the Hessian and the Jacobian are taken from the SMOOTH
        // Gaussian drift -- the hard min being only piecewise linear and, at
        // saturation, kinked exactly at the fixed point.
        FluidOptions mfo = opt;
        mfo.method = "closing";
        mfo.closure = FluidClosure();
        const FluidSolution mf = detail::fluid_dispatch(sn, mfo);
        iters += mf.iters;
        const std::vector<double> xmf = mf.xvec;

        const Matrix<double> A = fluid_drift_jacobian(terms, xmf, drift_cl);
        Sigma = fluid_moment_lyapunov(terms, A, fluid_moment_rates(terms, xmf, drift_cl));
        FluidRefineInfo rinfo;
        // A LINEAR DRIFT NEEDS NO REFINEMENT, and that is not the degenerate
        // call `fluid_refine_meanfield` refuses. When every station is either
        // an infinite server or `min_exact` -- min(n,c) is the identity on the
        // reachable set, the population bound never reaching c -- the drift is
        // exactly affine there, its Hessian vanishes and the O(1/N) correction
        // is identically zero. The mask above then zeroes all of the drift
        // closure, which `gaussian()` reads as "the caller handed me the first
        // order closure" and the refinement rejects. Settle it here, where the
        // reason for the zero is known: a null correction, not an error.
        // Delay + PS(c=2) at N=2 is the smallest case.
        bool drift_is_linear = true;
        for (std::size_t i = 0; i < M && drift_is_linear; ++i)
            if (!terms.min_exact[i] && !no_drift_var[i]) drift_is_linear = false;
        if (drift_is_linear)
            refinement.assign(xmf.size(), 0.0);
        else
            refinement = fluid_refine_meanfield(terms, xmf, drift_cl, Sigma, rinfo);
        x = xmf;
        for (std::size_t a = 0; a < x.size(); ++a) {
            x[a] += refinement[a];
            if (x[a] < 0.0) x[a] = 0.0;
        }
        // The corrected point is a correction OF the mean-field fixed point, so
        // its rates are read with the mean-field (zero) variance.
        drift_cl = FluidClosure();
        r = fluid_moment_rates(terms, x, drift_cl);
        for (std::size_t i = 0; i < M; ++i) {
            const std::vector<std::size_t>& blk = terms.station_block[i];
            if (blk.empty()) continue;
            double acc = 0.0;
            for (std::size_t a = 0; a < blk.size(); ++a)
                for (std::size_t b = 0; b < blk.size(); ++b) acc += Sigma(blk[a], blk[b]);
            cl.sigma2[i] = std::max(0.0, acc);
        }
    } else {
        r = fluid_moment_rates(terms, x, drift_cl);
    }

    // ---- performance measures, read off the event representation -------------
    const std::vector<double> gfac = fluid_moment_factors(terms, x, drift_cl);

    // A load-dependent station clears alpha(n) times the nominal work, so its
    // utilization normalises by the PEAK scaling (T*S/peak, as in the CTMC).
    std::vector<double> Seff = terms.S;
    for (std::size_t i = 0; i < M; ++i) {
        const std::vector<double>& lld = terms.sys.lld[i];
        for (std::size_t k = 0; k < lld.size(); ++k) Seff[i] = std::max(Seff[i], lld[k]);
    }

    FluidSolution out;
    out.method = m;
    out.iters = iters;
    out.xvec = x;
    out.QN = Matrix<double>(M, K, 0.0);
    out.UN = Matrix<double>(M, K, 0.0);
    out.RN = Matrix<double>(M, K, 0.0);
    out.TN = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const std::vector<std::size_t>& blk = terms.class_block[i][k];
            if (blk.empty()) continue;
            double q = 0.0, g = 0.0;
            for (std::size_t a = 0; a < blk.size(); ++a) {
                q += x[blk[a]];
                g += gfac[blk[a]];
            }
            out.QN(i, k) = q;
            out.UN(i, k) = (terms.sys.sched[i] == lang::SchedStrategy::INF) ? q : g / Seff[i];
            // Summed over ORIGINAL events through `emap`: a reduced event folded
            // through an immediate coordinate is a completion at more than one
            // (station,class), and its rate has to reach every one of them.
            double tn = 0.0;
            for (std::size_t e = 0; e < r.size() && e < terms.emap.rows(); ++e) {
                double w = 0.0;
                for (std::size_t o = 0; o < terms.ev_is_departure.size(); ++o)
                    if (terms.ev_is_departure[o] && terms.ev_station[o] == i
                        && terms.ev_class[o] == k)
                        w += terms.emap(e, o);
                if (w != 0.0) tn += r[e] * w;
            }
            out.TN(i, k) = tn;
            // TN is zero only to the integrator's accuracy: a class that never visits leaves
            // a ~1e-20 residue in TN too, and a strict > 0 test then divides residue by residue.
            if (tn > lang::GlobalConstants::Zero) out.RN(i, k) = q / tn;
        }

    // A Source and a Sink report no queue length, utilization or response time,
    // the same rule `solver_fluid` applies to the first-order methods and
    // `NetworkSolver.zeroSourceMetrics` applies in the reference. The closing
    // representation holds UNIT MASS at an EXT source so that `g` can read
    // `1 - rest` (see `fluid_odes.h`), and that coordinate is a normalisation
    // constant rather than a job count: `fluid_moment_terms` already projects it
    // out of the covariance, but `class_block` still spans it, so the metric loop
    // above reads the pool as a queue. Measured on m3 (Exp(0.5) -> Erlang(2) PS)
    // the source row came back QLen 0.2497 and RespT 0.4995 against 0 and 0 in
    // the reference, and the pool also entered `CN` below.
    for (std::size_t i = 0; i < M; ++i) {
        const qn::NodeType nt = sn.stations[i].nodetype;
        if (nt != qn::NodeType::Source && nt != qn::NodeType::Sink) continue;
        for (std::size_t k = 0; k < K; ++k) {
            out.QN(i, k) = 0.0;
            out.UN(i, k) = 0.0;
            out.RN(i, k) = 0.0;
        }
    }

    // ---- the moment report --------------------------------------------------
    FluidMomentReport rep;
    rep.Sigma = Sigma;
    rep.QVar = Matrix<double>(M, K, 0.0);
    rep.QStd = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const std::vector<std::size_t>& blk = terms.class_block[i][k];
            if (blk.empty()) continue;
            double acc = 0.0;
            for (std::size_t a = 0; a < blk.size(); ++a)
                for (std::size_t b = 0; b < blk.size(); ++b) acc += Sigma(blk[a], blk[b]);
            rep.QVar(i, k) = std::max(0.0, acc);
            rep.QStd(i, k) = std::sqrt(rep.QVar(i, k));
        }
    rep.sigma2 = cl.sigma2;
    rep.refinement = refinement;
    rep.outer_iters = outer_iters;
    rep.class_block = terms.class_block;
    out.has_moments = true;
    out.moments = rep;
    out.closure = drift_cl;

    // System throughput and response time, per chain reference station.
    out.XN.assign(K, 0.0);
    out.CN.assign(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t rs = sn.classes[k].refstat;
        if (rs >= 1 && rs <= M) out.XN[k] = out.TN(rs - 1, k);
        double q = 0.0;
        for (std::size_t i = 0; i < M; ++i) q += out.QN(i, k);
        if (out.XN[k] > 0.0) out.CN[k] = q / out.XN[k];
    }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_MOMENTS_H
