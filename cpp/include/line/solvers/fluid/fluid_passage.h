/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_PASSAGE_H
#define LINE_SOLVERS_FLUID_FLUID_PASSAGE_H

/**
 * Response-time distribution by tagged fluid: a port of
 * `solver_fluid_passage_time.m`, which is what `@@SolverFLD/getCdfRespT`
 * delegates to.
 *
 * THE IDEA. The fluid drift gives means, not distributions. To recover one,
 * mark the fluid that is at the station of interest RIGHT NOW and watch it
 * leave: if `z(t)` is how much of that marked fluid is still there at time t
 * and `z(0)` is how much there was, then
 *
 *     P(response time <= t) = 1 - z(t) / z(0)
 *
 * because a marked drop has completed exactly when it is no longer in the
 * block. This is a passage time read off a deterministic trajectory, which is
 * why it needs no state space.
 *
 * HOW THE MARKING IS DONE. The reference adds a whole extra CLASS to the model,
 * a copy of class c that routes as c does everywhere except at the observed
 * station, where it routes into the ORIGINAL classes. That is a rebuild of the
 * routing table for a class that, by construction, can only ever hold mass at
 * one station: the moment it completes there it becomes untagged. This port
 * therefore adds the tagged block at THAT STATION ONLY and mirrors the events
 * sourced there, which is the same system with none of the table surgery.
 *
 * WHAT THE TAGGED FLUID STILL PARTICIPATES IN. Its own service, and the
 * station's occupancy: a processor-sharing station splits its capacity over
 * everything present, marked or not, so the tagged block must be counted in
 * `ni` or the marked fluid would drain as if the station were emptier than it
 * is. That coupling is the reason this cannot be computed from the mean
 * trajectory alone.
 *
 * THE INITIAL STATE, as the reference sets it: the whole (station, class)
 * block is moved into PHASE ONE of the tagged block, and the original block
 * starts empty. Fluid arriving afterwards is untagged and is not measured.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"
#include "line/util/lsoda.h"

namespace line {
namespace fluid {

/** The response-time CDF of one (station, class), sampled on a grid. */
struct FluidPassage {
    std::vector<double> t;
    std::vector<double> cdf;
    double fluid0 = 0.0;  ///< the marked mass at t = 0; zero means nothing to measure
};

/**
 * Response-time CDF at station `ist` for class `cls`, both 1-based.
 *
 * `x_steady` is the converged fluid state the marking starts from -- the
 * reference passes `options.init_sol = odeStateVec`, i.e. the steady state, so
 * the distribution is the stationary one.
 *
 * `closure` is the variance the mean solve closed its drift at, i.e.
 * `FluidSolution::closure`. This is a SECOND solve on that solve's fixed point,
 * so it has to be driven by the same drift: closing `min(n_i,c_i)` at zero
 * variance drains a station the mean solve holds below capacity at full rate,
 * and the distribution then contradicts the mean the same solver reports. A
 * default-constructed closure is the first-order drift, which is what every
 * first-order method wants.
 */
template <class T>
FluidPassage fluid_passage_time(const qn::NetworkStruct<T>& sn, const std::vector<double>& x_steady,
                                std::size_t ist, std::size_t cls, double tol = 1e-4,
                                std::size_t points = 201,
                                const FluidClosure& closure = FluidClosure()) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    if (ist == 0 || ist > M) throw InputError("fluid_passage_time: station out of range");
    if (cls == 0 || cls > K) throw InputError("fluid_passage_time: class out of range");
    const std::size_t i = ist - 1, c = cls - 1;

    FluidOdeSystem sys = fluid_ode_system(sn);
    // Only the per-STATION variance travels: the marked block relabels a
    // station population, which leaves `sigma2` meaningful, whereas the
    // coordinate covariance is indexed by a layout the tagged block extends, so
    // the class share stays the plug-in ratio.
    sys.closure.sigma2 = closure.sigma2;
    const FluidLayout& L = sys.layout;
    if (x_steady.size() != L.nstates)
        throw InputError("fluid_passage_time: the initial state has the wrong length");
    if (sn.stations[i].nodetype == qn::NodeType::Source)
        throw InputError("fluid_passage_time: a Source has no response time");

    // THE DEGENERATE CURVE IS STILL A CURVE. A class that holds no fluid at
    // this station completes instantly, and the reference reports that as F = 1
    // over the SAME grid every other class gets (`RT{i,c,2} = ones(size(fullt))`
    // when `fluid_c == 0`), not as a single point. A one-point curve is not a
    // distribution a consumer can read: `diff(F)' * t(2:end)`, which is how the
    // mean is taken off it, is empty there -- on cdf_respt_closed_threeclasses,
    // whose Class2 and Class3 are declared with population 0, that aborted the
    // example with "the size of the right side is 0-by-1". Two points carry the
    // same law and the same mean (0) and are readable.
    const double flat_end = 100.0 / detail::fluid_slow_rate(sn, sys.layout, tol);
    const std::size_t P = L.kic[i][c];
    FluidPassage out;
    if (P == 0) {  // the class is not served here: nothing to measure
        out.t.push_back(0.0);
        out.t.push_back(flat_end);
        out.cdf.push_back(1.0);
        out.cdf.push_back(1.0);
        return out;
    }

    const std::size_t n = L.nstates;      // untagged states
    const std::size_t tag0 = n;           // the tagged block starts here
    const std::size_t nt = n + P;         // augmented size

    // The marked mass: the whole block, collapsed into phase one.
    double fluid0 = 0.0;
    for (std::size_t k = 0; k < P; ++k) fluid0 += x_steady[L.qidx[i][c] + k];
    out.fluid0 = fluid0;
    if (!(fluid0 > 0.0)) {  // an empty block completes instantly
        out.t.push_back(0.0);
        out.t.push_back(flat_end);
        out.cdf.push_back(1.0);
        out.cdf.push_back(1.0);
        return out;
    }

    std::vector<double> y0(nt, 0.0);
    for (std::size_t s = 0; s < n; ++s) y0[s] = x_steady[s];
    for (std::size_t k = 0; k < P; ++k) y0[L.qidx[i][c] + k] = 0.0;  // block emptied
    y0[tag0] = fluid0;                                               // all into phase one

    // The tagged copies of every event sourced at (i, c): a departure takes
    // marked fluid OUT of the system being measured and delivers it untagged;
    // a phase change keeps it marked.
    struct TagEvent {
        std::size_t minus, plus;
        std::size_t event_idx;
        double rate_base;
        bool leaves;  // true when the fluid stops being marked
    };
    std::vector<TagEvent> tev;
    const std::size_t base = L.qidx[i][c];
    for (std::size_t e = 0; e < sys.events.size(); ++e) {
        const FluidEvent& ev = sys.events[e];
        if (ev.event_idx < base || ev.event_idx >= base + P) continue;  // not sourced here
        const std::size_t k = ev.event_idx - base;
        TagEvent t2;
        t2.event_idx = tag0 + k;
        t2.rate_base = ev.rate_base;
        t2.minus = tag0 + k;
        if (e < sys.n_departures) {
            t2.plus = ev.plus;  // arrives untagged wherever the class routes
            t2.leaves = true;
        } else {
            t2.plus = tag0 + (ev.plus - base);  // stays marked, next phase
            t2.leaves = false;
        }
        tev.push_back(t2);
    }

    // The drift of the augmented system. The only subtlety is that station i's
    // occupancy must include the tagged block, so its service sharing is right.
    const std::size_t Kc = K;
    const LsodaRhs f = [&sys, &L, &tev, i, c, base, tag0, P, n, nt, Kc](double, const double* x,
                                                                       double* dx) {
        std::vector<double> xb(x, x + n);
        // Fold the marked mass back into the block it came from, purely to
        // compute the station's occupancy and therefore its service share.
        double tagsum = 0.0;
        for (std::size_t k = 0; k < P; ++k) tagsum += x[tag0 + k];
        std::vector<double> g(xb);
        for (std::size_t k = 0; k < P; ++k) g[base + k] += x[tag0 + k];
        std::vector<double> gg(g);
        fluid_rates_closing(sys, g.data(), gg);
        // The share the station gives to this block, as a factor.
        double blk = 0.0, gblk = 0.0;
        for (std::size_t k = 0; k < P; ++k) {
            blk += g[base + k];
            gblk += gg[base + k];
        }
        const double share = (blk > 0.0) ? gblk / blk : 1.0;

        for (std::size_t s = 0; s < nt; ++s) dx[s] = 0.0;
        // Untagged events, evaluated on the folded state so the untagged part
        // of the block gets its correct share too.
        for (const FluidEvent& ev : sys.events) {
            double drive = gg[ev.event_idx];
            if (ev.event_idx >= base && ev.event_idx < base + P) {
                // Only the UNTAGGED part of this block drives untagged events.
                const std::size_t k = ev.event_idx - base;
                drive = xb[base + k] * share;
            }
            const double r = ev.rate_base * drive;
            if (r == 0.0) continue;
            dx[ev.minus] -= r;
            dx[ev.plus] += r;
        }
        // Tagged events, driven by the marked mass at the same share.
        for (const TagEvent& te : tev) {
            const double r = te.rate_base * x[te.event_idx] * share;
            if (r == 0.0) continue;
            dx[te.minus] -= r;
            dx[te.plus] += r;
        }
        (void)tagsum;
        (void)Kc;
    };

    // Integrate until the marked fluid is gone. The horizon follows the
    // reference: 100 events of the slowest rate, extended while mass remains.
    const double min_rate = detail::fluid_slow_rate(sn, L, tol);
    const double t_end = 100.0 / min_rate;

    std::vector<double> grid(points);
    for (std::size_t j = 0; j < points; ++j)
        grid[j] = t_end * static_cast<double>(j) / static_cast<double>(points - 1);

    LsodaOptions lopt;
    lopt.rtol = tol;
    lopt.atol = tol;

    // THE GRID IS REFINED WHERE THE CDF JUMPS, and it has to be. The curve is
    // read back by quadrature -- SolverLN's `moment3` takes its first three
    // moments off it -- and a uniform grid over a horizon set by the SLOWEST
    // rate resolves the fast rise near the origin with a handful of points, so
    // every moment comes out biased HIGH: the mass that arrives inside the first
    // interval is charged at that interval's midpoint. On a two-layer LQN the
    // entry service times came out 1 to 7 per cent above the JAR's until this
    // loop was added, with nothing in the output to say why.
    //
    // The rule is the reference's (SolverFluid.passageTime): while some adjacent
    // pair of CDF values differs by more than `kMaxJump`, insert points inside
    // the offending intervals and integrate again. The reference walks its own
    // ODE output grid and refines ONE interval per round; here the grid is
    // supplied to LSODA, so every offending interval is split in the same round
    // and the whole trajectory is re-integrated -- same fixed point, fewer
    // rounds. Both caps are bounds on work, not on accuracy: convergence is the
    // jump test, and hitting a cap leaves a coarser curve rather than a wrong one.
    const double kMaxJump = 0.0005;
    const int kMaxRounds = 5;
    const std::size_t kMaxPoints = 20001;
    const std::size_t kSplit = 20;

    auto integrate_on = [&](const std::vector<double>& g) {
        const LsodaSolution s = fluid_integrate_grid(f, y0, g, lopt);
        FluidPassage p;
        p.fluid0 = fluid0;
        p.t.reserve(s.y.size());
        p.cdf.reserve(s.y.size());
        for (std::size_t j = 0; j < s.y.size(); ++j) {
            double z = 0.0;
            for (std::size_t k = 0; k < P; ++k) z += std::max(0.0, s.y[j][tag0 + k]);
            double v = 1.0 - z / fluid0;
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            p.t.push_back(s.t[j]);
            p.cdf.push_back(v);
        }
        return p;
    };

    FluidPassage cur = integrate_on(grid);
    for (int round = 0; round < kMaxRounds; ++round) {
        if (cur.t.size() >= kMaxPoints) break;
        std::vector<double> next;
        next.reserve(cur.t.size() * 2);
        bool refined = false;
        for (std::size_t j = 0; j + 1 < cur.t.size(); ++j) {
            next.push_back(cur.t[j]);
            if (cur.cdf[j + 1] - cur.cdf[j] > kMaxJump && cur.t[j + 1] > cur.t[j]) {
                refined = true;
                const double dt = (cur.t[j + 1] - cur.t[j]) / static_cast<double>(kSplit);
                for (std::size_t s = 1; s < kSplit; ++s) next.push_back(cur.t[j] + dt * s);
            }
        }
        if (!cur.t.empty()) next.push_back(cur.t.back());
        if (!refined || next.size() > kMaxPoints) break;
        cur = integrate_on(next);
    }

    out.t.swap(cur.t);
    out.cdf.swap(cur.cdf);
    return out;
}

/**
 * Port of `@@SolverFLD/getTranCdfPassT`: the same passage-time distribution
 * started from the model's INITIAL state rather than from its steady state.
 *
 * IT IS THE SAME ALGORITHM AT A DIFFERENT STARTING POINT, which is exactly what
 * the reference is: `getCdfRespT` passes the converged `odeStateVec`,
 * `getTranCdfPassT` passes `solver_fluid_initsol(sn, options)`. The two are kept
 * as separate named entry points because the quantity differs -- one is the
 * stationary response-time law, the other the law seen by a job marked while the
 * system is still where the model says it starts -- and a caller that had to
 * assemble the second by hand would have to know that, which is what a name is
 * for.
 *
 * WHAT IS NOT REPRODUCED, and it is a MODEL-LAYER restriction rather than an
 * algorithmic one: the reference first collapses `sn.state` to the first row of
 * its prior and ERRORS when more than one initial state carries non-zero prior
 * mass, because a passage time from a mixture of starting states is not one
 * distribution. `fluid_default_initsol` is the closed form of the single-state
 * decode (see fluid_closing.h) and there is no prior to collapse here, so the
 * refusal has no input in this port.
 */
template <class T>
FluidPassage fluid_tran_passage_time(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                     std::size_t cls, double tol = 1e-4,
                                     std::size_t points = 201) {
    const FluidLayout L = fluid_layout(sn);
    return fluid_passage_time(sn, detail::fluid_default_initsol(sn, L), ist, cls, tol, points);
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_PASSAGE_H
