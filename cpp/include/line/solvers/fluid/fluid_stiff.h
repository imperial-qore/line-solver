/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_STIFF_H
#define LINE_SOLVERS_FLUID_FLUID_STIFF_H

/**
 * Port of `ode_eliminate_immediate.m`, `eliminate_immediate_matrix.m` and
 * `ode_solve_stiff.m`: the two answers to an IMMEDIATE transition.
 *
 * WHERE THE STIFFNESS COMES FROM. An Immediate distribution fires at
 * `GlobalConstants.Immediate`, 1e8, while the rest of the model runs at rates
 * of order one. The drift then has a mode with time constant 1e-8 alongside
 * modes with time constant 1, and an explicit step controller is pinned to the
 * fastest one for the WHOLE integration, long after that mode has died: the run
 * either crawls or goes unstable. That is stiffness, and it is a property of
 * the equations rather than of the integrator.
 *
 * THE TWO ANSWERS ARE NOT EQUALLY GOOD, AND THE REFERENCE PREFERS THE FIRST.
 * Eliminating the immediate transitions removes the fast mode from the system,
 * so what is left is not stiff at all. Falling back to a stiff integrator keeps
 * the fast mode and pays an implicit solve to stay stable across it. The first
 * is ALGEBRA and the second is NUMERICS, which is why the first is exact:
 *
 *   ELIMINATION IS EXACT. The stochastic complement of a generator over a
 *   retained set is the generator of the process WATCHED ONLY ON THAT SET,
 *   S = Q11 + Q12 (-Q22)^-1 Q21. Its stationary law is the original's,
 *   conditioned on the retained set and renormalized. It is an identity, not an
 *   approximation, and it holds at any rate ratio. The one place a gap appears
 *   is when the answer is compared against the UNREDUCED system: the immediate
 *   states hold stationary mass of order 1/Immediate, about 1e-8, which the
 *   reduced system does not carry and the conditioning divides out.
 *
 *   THE STIFF INTEGRATOR ONLY COPES. It keeps the fast mode and controls the
 *   error on it, so its answer carries the integrator's tolerance and nothing
 *   better.
 *
 * A CAVEAT THE REFERENCE DOES NOT STATE. Reconstructing a generator from the
 * jump/rate representation reads the drift as `dx/dt = x W`, which is what the
 * fluid drift IS wherever the state-dependent factor `g(x)` is x itself: an
 * infinite server, or any station holding fewer jobs than it has servers. At a
 * SATURATED station g is `min(n_i, S_i)/n_i` times x and the reconstructed W is
 * not the drift there, so the complement is exact for the linear part only.
 * That is the reference's behaviour, reproduced rather than corrected, and it
 * is why the elimination is applied to the transitions and not to the metrics.
 *
 * WHAT WAS ALREADY HERE, AND IS NOT REBUILT. `mc::ctmc_stochcomp` is the
 * complement, with the shared LU factorization the reference gets from
 * backslash. `ode_rosenbrock4` in util/ode.h is the stiff integrator: a
 * four-stage L-stable Rosenbrock method, the family MATLAB's ode23s belongs to,
 * with an embedded estimate for step control. And the fluid path's ordinary
 * integrator is LSODA, which already switches itself from Adams to BDF when it
 * detects stiffness; the point of this file is therefore not "a stiff solver
 * exists" but that the stiffness can be removed before anyone integrates.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/mc/dtmc_solve.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/linalg.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/util/error.h"
#include "line/util/lsoda.h"
#include "line/util/matrix.h"
#include "line/util/ode.h"

namespace line {
namespace fluid {

/**
 * The reference's two thresholds, which do NOT agree and are not meant to.
 *
 * `ode_eliminate_immediate` knows the rates it is looking at are `rateBase`
 * entries built from an Immediate distribution, so it asks for a rate within 1%
 * of the constant itself. `eliminate_immediate_matrix` is handed a generator
 * whose entries have already been multiplied by routing probabilities and
 * summed, so it asks only for an order of magnitude below it. A transition
 * routed with probability 0.5 is invisible to the first rule and caught by the
 * second; both are kept as they are, because tightening the first would start
 * eliminating transitions the reference integrates.
 */
inline double fluid_immediate_transition_tol() {
    return lang::GlobalConstants::Immediate * 0.99;
}
inline double fluid_immediate_state_tol() { return lang::GlobalConstants::Immediate / 10.0; }

/** What an elimination attempt produced. */
struct FluidImmediateResult {
    FluidOdeSystem sys;                  ///< reduced, or the input on a fallback
    std::vector<std::size_t> state_map;  ///< reduced position -> original index, 0-based
    bool eliminated = false;             ///< false when the input is returned unchanged
    std::size_t n_immediate = 0;         ///< immediate transitions detected
    /**
     * Why the elimination was abandoned, empty when it was not attempted or
     * when it succeeded. The reference warns here and carries on with a system
     * that is still stiff; a caller that cannot see the warning would integrate
     * it without knowing to switch to the stiff arm, which is the whole
     * decision this file exists to inform.
     */
    std::string fallback;
    /**
     * `emap(e, o)`: expected firings of the ORIGINAL event o per firing of the
     * reduced event e; the identity when nothing was eliminated. A caller maps a
     * per-event quantity with `new = emap * old`, which is what lets a throughput
     * read off a reduced event set stay exact -- an event folded through an
     * immediate coordinate is a completion at more than one (station,class).
     */
    Matrix<double> emap;
    /**
     * Projector for the initial condition: identity on the timed rows, the
     * absorption distribution on the immediate ones. Mass parked on an eliminated
     * coordinate would otherwise be frozen there for the whole integration.
     */
    Matrix<double> absorb;
};

namespace detail {

/**
 * The generator the jump/rate representation encodes, `W` in the reference.
 *
 * The row is the GATING index, `event_idx`, and not the coordinate the event
 * takes mass from: the rate of an event is `rate_base * g(x)_{event_idx}`, so
 * that is the index the rate is proportional to and therefore the only one a
 * generator row can be. The two coincide in every system `fluid_ode_system`
 * builds, and the caller below refuses to complement one where they do not
 * rather than silently complementing a matrix that is not the drift.
 *
 * An event whose endpoints coincide contributes nothing to the drift (it adds
 * and removes the same mass), and it is dropped here for the same reason the
 * reference's row-sum subtraction cancels it.
 */
inline Matrix<double> fluid_generator_from_events(const FluidOdeSystem& sys, std::size_t n) {
    Matrix<double> W(n, n, 0.0);
    for (const FluidEvent& e : sys.events) {
        if (e.minus == e.plus) continue;
        W(e.event_idx, e.plus) += e.rate_base;
    }
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            if (j != i) s += W(i, j);
        W(i, i) = -s;
    }
    return W;
}

/**
 * Port of `generator_to_jumps.m`, with the reduced indices mapped back.
 *
 * Only strictly positive off-diagonal entries become events. A complement can
 * emit an entry of size 1e-300 where two paths nearly cancel, and the reference
 * keeps those too; what it never keeps is the diagonal, which is not a
 * transition but the negative of the row sum the jumps rebuild by themselves.
 */
inline std::vector<FluidEvent> fluid_events_from_generator(
    const Matrix<double>& W, const std::vector<std::size_t>& state_map) {
    std::vector<FluidEvent> out;
    for (std::size_t i = 0; i < W.rows(); ++i)
        for (std::size_t j = 0; j < W.cols(); ++j) {
            if (i == j || !(W(i, j) > 0.0)) continue;
            FluidEvent e;
            e.minus = state_map[i];
            e.plus = state_map[j];
            e.event_idx = state_map[i];
            e.rate_base = W(i, j);
            out.push_back(e);
        }
    return out;
}

}  // namespace detail

/**
 * Stochastic complementation of the INSTANTANEOUS coordinates of a fluid drift,
 * the twin of `ode_eliminate_immediate.m`.
 *
 * A coordinate whose exit rate is `GlobalConstants::Immediate` is not a fast
 * coordinate, it is an INSTANTANEOUS one: the rate is LINE's stand-in for
 * infinity, written by SolverLN for the branch of an activity that takes no time
 * (an entry called with probability y < 1 carries a second PH phase at InfRate
 * entered with probability 1-y). Integrating it is meaningless work no
 * integrator does well.
 *
 * THE REDUCTION IS EXACT. The instantaneous coordinates F are absorbed into the
 * timed ones S by the absorption probabilities of the embedded jump chain
 * restricted to F, so the flow that would enter F is routed straight to where F
 * would have sent it.
 *
 * WHY THIS IS A STRUCTURAL COMPOSITION AND NOT A GENERATOR ROUND TRIP, which is
 * what this function used to be. Every event is a single -1 at `event_idx` and a
 * single +1 at `plus`, so a path through F composes to ONE event, -1 at the
 * original source and +1 at the absorbing coordinate, that keeps the original
 * source's GATING. Rebuilding the events from a reduced generator loses that
 * identity -- and with it `n_departures`, which this function used to report as
 * zero because it was "no longer recoverable". It is recoverable: `emap(e, o)`
 * is the expected number of times the ORIGINAL event o fires per firing of the
 * reduced event e, so a caller maps any per-event quantity with
 * `new = emap * old` and gets an exact rate accounting. That is what lets the
 * moment-closure methods, which used to refuse the reduction outright, read
 * their throughputs off a reduced event set.
 *
 * A COMPOSED EVENT CAN BE A DEPARTURE AT TWO STATIONS AT ONCE: a job that leaves
 * a delay, passes through a queue's immediate phase and returns has completed at
 * both, and both throughputs must count it. `emap` gives it a row with weight on
 * both original events, and the null jump it composes to (-1 and +1 on the same
 * coordinate) correctly contributes nothing to the drift and nothing to the
 * diffusion.
 *
 * `absorb` projects an initial condition onto the surviving coordinates: mass
 * parked on an eliminated one would otherwise be frozen there for the whole
 * integration, because nothing moves it any more.
 */
/**
 * Whether this model's fluid drift is built on the stochastic complement of its
 * INSTANTANEOUS coordinates.
 *
 * Every fluid route that builds its drift from the station/class/phase event set
 * or from the linear generator asks here rather than reading the flag directly,
 * so the answer is the same across matrix, closing, statedep, tbi, minnormal,
 * refined and dae. The flag defaults to TRUE: a coordinate whose exit rate is
 * `GlobalConstants::Immediate` is LINE's stand-in for infinity, and integrating
 * it is meaningless work no integrator does well.
 *
 * THE STOCHASTIC PETRI NET ROUTE IS THE ONE EXCEPTION, and it is not a refusal.
 * It carries immediate firings as ALGEBRAIC unknowns of an index-1 DAE, a
 * stronger treatment than absorbing them, and never builds the event set this
 * reduction acts on, so the answer here is simply false.
 */
template <class T, class Opt>
bool fluid_hide_immediate(const qn::NetworkStruct<T>& sn, const Opt& opt) {
    if (!opt.hide_immediate) return false;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Transition) return false;
    return true;
}

inline FluidImmediateResult fluid_eliminate_immediate(
    const FluidOdeSystem& sys, double imm_tol = fluid_immediate_transition_tol()) {
    const std::size_t n = sys.layout.nstates;
    const std::size_t ne = sys.events.size();
    FluidImmediateResult out;
    out.sys = sys;
    out.state_map.resize(n);
    for (std::size_t i = 0; i < n; ++i) out.state_map[i] = i;
    out.emap = Matrix<double>(ne, ne, 0.0);
    for (std::size_t e = 0; e < ne; ++e) out.emap(e, e) = 1.0;
    out.absorb = Matrix<double>(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) out.absorb(i, i) = 1.0;

    // The immediate coordinates are the SOURCES of the immediate events: it is
    // the coordinate that empties instantaneously, not the event.
    std::vector<bool> is_imm(n, false);
    for (const FluidEvent& e : sys.events)
        if (e.rate_base >= imm_tol) {
            ++out.n_immediate;
            if (e.event_idx < n) is_imm[e.event_idx] = true;
        }
    if (out.n_immediate == 0) return out;  // nothing to eliminate, and no warning to give

    for (const FluidEvent& e : sys.events)
        if (e.minus != e.event_idx) {
            // The rate would be proportional to one coordinate while the mass
            // left another, which no branching chain can express.
            out.fallback =
                "an event draws its rate from a coordinate other than the one it removes mass "
                "from, so the drift is not a generator and cannot be complemented";
            return out;
        }

    // A coordinate with no outflow cannot be complemented away, and one whose
    // outflow is entirely a self-loop would make the fundamental matrix
    // singular. Both are dropped rather than guessed at.
    for (std::size_t f = 0; f < n; ++f) {
        if (!is_imm[f]) continue;
        double tot = 0.0;
        bool leaves = false;
        for (const FluidEvent& e : sys.events)
            if (e.event_idx == f) {
                tot += e.rate_base;
                if (e.plus != f) leaves = true;
            }
        if (!(tot > 0.0) || !leaves) is_imm[f] = false;
    }

    std::vector<std::size_t> Fidx, Sidx;
    for (std::size_t i = 0; i < n; ++i) (is_imm[i] ? Fidx : Sidx).push_back(i);
    if (Fidx.empty() || Sidx.size() <= 1) {
        out.fallback =
            "every fluid coordinate but at most one sources an immediate transition, so the "
            "complement would be a trivial system; integrating the original stiff drift instead";
        return out;
    }
    const std::size_t nF = Fidx.size(), nS = Sidx.size();
    std::vector<std::size_t> posF(n, 0), posS(n, 0);
    for (std::size_t a = 0; a < nF; ++a) posF[Fidx[a]] = a;
    for (std::size_t b = 0; b < nS; ++b) posS[Sidx[b]] = b;

    // Branching of the embedded jump chain out of each immediate coordinate. The
    // probabilities are the rate shares, so a coordinate carrying both an
    // immediate and an ordinary exit gives the ordinary one its (vanishing)
    // share rather than being special-cased.
    Matrix<double> PFF(nF, nF, 0.0), PFS(nF, nS, 0.0), cnt(nF, ne, 0.0);
    for (std::size_t a = 0; a < nF; ++a) {
        const std::size_t f = Fidx[a];
        double tot = 0.0;
        for (const FluidEvent& e : sys.events)
            if (e.event_idx == f) tot += e.rate_base;
        for (std::size_t o = 0; o < ne; ++o) {
            const FluidEvent& e = sys.events[o];
            if (e.event_idx != f) continue;
            const double p = e.rate_base / tot;
            cnt(a, o) += p;
            if (is_imm[e.plus]) PFF(a, posF[e.plus]) += p;
            else PFS(a, posS[e.plus]) += p;
        }
    }

    // Fundamental matrix of the instantaneous chain. A closed cycle of immediate
    // coordinates has no absorption distribution and is left unreduced.
    Matrix<double> ImP(nF, nF, 0.0);
    for (std::size_t a = 0; a < nF; ++a)
        for (std::size_t b = 0; b < nF; ++b)
            ImP(a, b) = (a == b ? 1.0 : 0.0) - PFF(a, b);
    Matrix<double> Nfm;
    try {
        Nfm = inverse(ImP);
    } catch (const std::exception& ex) {
        out.fallback = std::string("the immediate coordinates form a closed cycle, so they have "
                                   "no absorption distribution: ") + ex.what();
        return out;
    }
    const Matrix<double> Aabs = matmul(Nfm, PFS);
    const Matrix<double> expcnt = matmul(Nfm, cnt);
    for (std::size_t a = 0; a < nF; ++a) {
        double rowsum = 0.0;
        for (std::size_t b = 0; b < nS; ++b) {
            if (!std::isfinite(Aabs(a, b))) rowsum = std::numeric_limits<double>::quiet_NaN();
            rowsum += Aabs(a, b);
        }
        if (!(rowsum > 0.5)) {
            out.fallback =
                "stochastic complementation produced no absorption distribution, which is what a "
                "closed set of immediate coordinates looks like once it has been solved";
            return out;
        }
    }

    // Compose the event list. An event sourced in F is dropped: its flow is
    // already carried by whichever event feeds F. Departures stay first, so
    // `n_departures` survives the composition -- a composed event is a departure
    // exactly when the event that fed the immediate coordinate was one.
    std::vector<FluidEvent> kept_dep, kept_other;
    std::vector<std::vector<double>> emap_dep, emap_other;
    for (std::size_t o = 0; o < ne; ++o) {
        const FluidEvent& e = sys.events[o];
        if (is_imm[e.event_idx]) continue;
        const bool is_dep = o < sys.n_departures;
        std::vector<FluidEvent>& bucket = is_dep ? kept_dep : kept_other;
        std::vector<std::vector<double>>& rows = is_dep ? emap_dep : emap_other;
        if (!is_imm[e.plus]) {
            bucket.push_back(e);
            std::vector<double> row(ne, 0.0);
            row[o] = 1.0;
            rows.push_back(row);
            continue;
        }
        // The event feeds an immediate coordinate: one event per absorbing
        // destination, keeping the original source and so the original gating,
        // since the rate of the composed flow IS the rate of the inflow.
        const std::size_t a = posF[e.plus];
        for (std::size_t b = 0; b < nS; ++b) {
            if (!(Aabs(a, b) > 0.0)) continue;
            FluidEvent ce = e;
            ce.plus = Sidx[b];
            ce.rate_base = e.rate_base * Aabs(a, b);
            bucket.push_back(ce);
            // Weighting every absorbing branch by the SAME unconditional expected
            // counts is what makes the rate accounting exact: the branch rates
            // sum back to e.rate_base, so the mapped total is e.rate_base times
            // the counts.
            std::vector<double> row(ne, 0.0);
            for (std::size_t oo = 0; oo < ne; ++oo) row[oo] = expcnt(a, oo);
            row[o] += 1.0;
            rows.push_back(row);
        }
    }

    out.sys.events.clear();
    out.sys.events.insert(out.sys.events.end(), kept_dep.begin(), kept_dep.end());
    out.sys.events.insert(out.sys.events.end(), kept_other.begin(), kept_other.end());
    out.sys.n_departures = kept_dep.size();

    const std::size_t nnew = out.sys.events.size();
    out.emap = Matrix<double>(nnew, ne, 0.0);
    for (std::size_t e = 0; e < emap_dep.size(); ++e)
        for (std::size_t o = 0; o < ne; ++o)
            if (emap_dep[e][o] != 0.0) out.emap(e, o) = emap_dep[e][o];
    for (std::size_t e = 0; e < emap_other.size(); ++e)
        for (std::size_t o = 0; o < ne; ++o)
            if (emap_other[e][o] != 0.0) out.emap(emap_dep.size() + e, o) = emap_other[e][o];

    out.absorb = Matrix<double>(n, n, 0.0);
    for (std::size_t b = 0; b < nS; ++b) out.absorb(Sidx[b], Sidx[b]) = 1.0;
    for (std::size_t a = 0; a < nF; ++a)
        for (std::size_t b = 0; b < nS; ++b)
            if (Aabs(a, b) > 0.0) out.absorb(Fidx[a], Sidx[b]) = Aabs(a, b);

    out.state_map = Sidx;
    out.eliminated = true;
    return out;
}

/**
 * The same, at the reference's own signature, which carries `sn`.
 *
 * It is here for the arithmetic gate rather than for the struct: the fluid
 * drift is integrated by LSODA and by the Rosenbrock method below, both of
 * which are double by construction, so a non-double model has no business
 * reaching either. The complement ITSELF is field arithmetic and would be exact
 * over the rationals, which is why the gate sits on the entry point that knows
 * what the model is made of and not on the complement.
 */
template <class T>
FluidImmediateResult fluid_eliminate_immediate(const qn::NetworkStruct<T>& sn,
                                               const FluidOdeSystem& sys,
                                               double imm_tol = fluid_immediate_transition_tol()) {
    (void)sn;
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "fluid_eliminate_immediate: the fluid solver integrates its drift with LSODA, whose "
            "coefficients assume double precision; rerun with --arith double");
    return fluid_eliminate_immediate(sys, imm_tol);
}

/** What the matrix-level elimination produced. */
template <class T>
struct FluidImmediateMatrix {
    Matrix<T> W;                         ///< reduced, or the input on a fallback
    std::vector<std::size_t> state_map;  ///< reduced position -> original index, 0-based
    bool eliminated = false;
    std::string fallback;
    /**
     * `emap(e, o)`: expected firings of the ORIGINAL event o per firing of the
     * reduced event e; the identity when nothing was eliminated. A caller maps a
     * per-event quantity with `new = emap * old`, which is what lets a throughput
     * read off a reduced event set stay exact -- an event folded through an
     * immediate coordinate is a completion at more than one (station,class).
     */
    Matrix<double> emap;
    /**
     * Projector for the initial condition: identity on the timed rows, the
     * absorption distribution on the immediate ones. Mass parked on an eliminated
     * coordinate would otherwise be frozen there for the whole integration.
     */
    Matrix<double> absorb;
};

/**
 * Port of `eliminate_immediate_matrix.m`: the same elimination on a generator
 * that is already assembled.
 *
 * It detects immediate STATES directly, by the largest rate on their row,
 * rather than immediate transitions. That is the only rule available once the
 * routing probabilities have been folded in and parallel edges summed, and it
 * catches a case the transition rule misses: an immediate transition taken with
 * probability one half carries rate 5e7, which is nowhere near Immediate but is
 * still seven orders above the rest of the model.
 */
template <class T>
FluidImmediateMatrix<T> fluid_eliminate_immediate_matrix(
    const Matrix<T>& W, double imm_tol = fluid_immediate_state_tol()) {
    const std::size_t n = W.rows();
    if (W.cols() != n)
        throw InputError("fluid_eliminate_immediate_matrix: the generator is not square");

    FluidImmediateMatrix<T> out;
    out.W = W;
    out.state_map.resize(n);
    for (std::size_t i = 0; i < n; ++i) out.state_map[i] = i;

    std::vector<std::size_t> timed;
    std::size_t n_imm = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double mx = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            mx = std::max(mx, std::fabs(num_traits<T>::to_double(W(i, j))));
        if (mx >= imm_tol) ++n_imm;
        else timed.push_back(i);
    }
    if (n_imm == 0) return out;
    if (timed.size() <= 1) {
        out.fallback =
            "at most one state of the generator is timed, so the complement would be a trivial "
            "system; the original generator is returned unreduced";
        return out;
    }

    try {
        Matrix<T> S = mc::ctmc_stochcomp(W, timed).S;
        for (std::size_t i = 0; i < S.rows(); ++i)
            for (std::size_t j = 0; j < S.cols(); ++j)
                if (!std::isfinite(num_traits<T>::to_double(S(i, j)))) {
                    out.fallback =
                        "stochastic complementation produced non-finite rates, which is what a "
                        "singular immediate block looks like once it has been solved";
                    return out;
                }
        out.W = S;
        out.state_map = timed;
        out.eliminated = true;
    } catch (const std::exception& ex) {
        out.fallback = std::string("stochastic complementation failed: ") + ex.what();
    }
    return out;
}

/**
 * Controls for the stiff arm.
 *
 * `stiff` is `options.stiff`, the reference's choice between its accurate and
 * its fast stiff solver. THAT CHOICE DOES NOT SURVIVE THE PORT as a choice of
 * method: MATLAB picks between ode15s and ode23s, and this port has exactly one
 * stiff method, so claiming either name would be a fiction. What the two arms
 * differ in that IS expressible here is the NonNegative handling, which the
 * reference clears for the fast arm because ode23s cannot honour it, and that
 * is what the flag selects below. Tolerances come from the caller on both arms,
 * as they do in the reference.
 */
struct FluidStiffOptions {
    /** Refuses any named MATLAB solver: see the gate below. */
    std::string solver = "default";
    bool stiff = true;             ///< options.stiff: keep the nonnegativity projection
    double rtol = 1e-4;            ///< the value solver_fluid.h hands LSODA, FluidOptions::tol
    double atol = 1e-4;
    std::size_t max_steps = 100000;
    bool store_trajectory = false;
    /** Per-accepted-step stop; see OdeOptions::step_stop and solver_fluid.h. */
    std::function<bool(const double&, const std::vector<double>&)> step_stop;
};

/**
 * Port of `ode_solve_stiff.m`.
 *
 * The reference dispatches through `options.odesolvers.*StiffOdeSolver`, which
 * are function handles a caller may point at any MATLAB integrator. There is no
 * such registry here and no way to honour an arbitrary one, so a solver asked
 * for BY NAME is refused rather than served by the method that happens to be
 * present under a name it does not have.
 *
 * THE NONNEGATIVITY IS A PROJECTION, NOT A CONSTRAINT. MATLAB's NonNegative
 * option is enforced by the step controller, which rejects a step that would
 * take a component below zero. Clamping the accepted trajectory afterwards
 * cannot rescue such a step; it prevents a negative mass from feeding back into
 * the drift as a negative rate, which is what `solver_fluid.h` already does
 * after every integration leg. The difference matters on a drift that is only
 * defined for nonnegative states, and is stated rather than hidden.
 */
inline OdeSolution<double> fluid_ode_solve_stiff(
    const std::function<void(double, const double*, double*)>& f, double t0, double t1,
    const std::vector<double>& y0, const FluidStiffOptions& opt = FluidStiffOptions()) {
    if (!(opt.solver == "default" || opt.solver == "rosenbrock4"))
        throw UnsupportedError(
            "fluid_ode_solve_stiff: the '" + opt.solver +
            "' integrator is a MATLAB solver handle this port does not carry; the stiff arm here "
            "is one four-stage L-stable Rosenbrock method (util/ode.h), asked for as 'default'");
    if (y0.empty()) throw InputError("fluid_ode_solve_stiff: the initial state is empty");
    if (!(t1 > t0))
        throw InputError(
            "fluid_ode_solve_stiff: the horizon must be positive; the fluid iteration marches "
            "forwards and a reversed range is a caller error rather than a backwards solve");

    OdeOptions<double> o;
    o.rtol = opt.rtol;
    o.atol = opt.atol;
    o.max_steps = opt.max_steps;
    o.store_trajectory = opt.store_trajectory;
    o.step_stop = opt.step_stop;

    const std::size_t n = y0.size();
    // util/ode.h takes and returns whole vectors, the fluid drift writes into a
    // raw buffer; the adapter is the only thing between them.
    const auto g = [&f, n](const double& t, const std::vector<double>& y) {
        std::vector<double> dy(n, 0.0);
        f(t, y.data(), dy.data());
        return dy;
    };
    OdeSolution<double> sol = ode_rosenbrock4<double>(g, t0, t1, y0, o);
    if (opt.stiff)
        for (std::size_t i = 0; i < sol.y.size(); ++i)
            for (std::size_t j = 0; j < sol.y[i].size(); ++j)
                if (sol.y[i][j] < 0.0) sol.y[i][j] = 0.0;
    return sol;
}

/**
 * The step budget LSODA gets on a fluid leg before the stiff arm takes over.
 *
 * `LsodaOptions::max_steps` defaults to the JAR's raised `mxstep` of 1e7, which
 * is a GIVE-UP threshold: nothing follows it, so it is set high enough that a
 * slow-but-finishing integration is not cut off. Here it is a SWITCH threshold
 * instead, because `fluid_integrate_leg` finishes the leg by another method, and
 * a switch wants to be cheap. A fluid leg that is going to finish takes tens to
 * a few thousand steps (measured: 3 to 121 on the LQN layers of `randomLQN`), so
 * 1e6 keeps a factor of a thousand in hand while costing a fraction of a second
 * to discover the Adams stability wall instead of the ~20 s that 1e7 cost.
 */
inline LsodaOptions fluid_lsoda(const LsodaOptions& lopt) {
    LsodaOptions out = lopt;
    if (out.max_steps > 1000000) out.max_steps = 1000000;
    return out;
}

/**
 * One integration leg, with the reference's retry on a failed solve.
 *
 * WHY A RETRY EXISTS AT ALL. `FluidOptions::stiff` is false in this port where
 * the reference's `options.stiff` is true, on the argument that LSODA switches
 * to BDF on its own stiffness detector and so already covers ode15s. That holds
 * on most drifts and NOT on all of them: an LQN layer whose entry carries an
 * Immediate rate has 1e8 in `sn.rates`, so the drift's Jacobian eigenvalue is
 * -1e8 while the trajectory is otherwise O(1), and LSODA has been measured
 * staying in ADAMS across such a leg -- 1e7 accepted steps at h = 5.6e-9,
 * reaching t = 0.056 of a horizon of 20.9 before `mxstep` stopped it. Adams is
 * stability-limited to h < 2/1e8 there, so no step count would have finished.
 *
 * WHAT THE FAILURE USED TO COST. Every caller took `final_state()` regardless:
 * `solver_fluid`'s closing loop broke out and reported the state it happened to
 * reach as the fixed point, while the matrix arm and the transient getter did
 * not read `success` at all. A layer that never left its initial transient was
 * published as a converged answer, which is the one outcome an integrator
 * failure must not produce.
 *
 * The reference retries too -- `solver_fluid_matrix.m` re-solves once with
 * `hide_immediate` toggled, `solver_fluid_iteration.m` catches the ODE error and
 * re-solves from the default initial state -- so retrying is the reference's
 * shape, not an invention. What is retried here is the INTEGRATOR: the same leg,
 * the same tolerances, run by the L-stable Rosenbrock arm above, which is the
 * family ode23s belongs to and is not stability-limited. A leg LSODA completes
 * is untouched, so no result that already converged moves.
 */
inline std::vector<double> fluid_integrate_leg(
    const std::function<void(double, const double*, double*)>& f, double t0, double t1,
    const std::vector<double>& y0, const LsodaOptions& lopt) {
    const LsodaSolution s = lsoda_integrate(f, y0, std::vector<double>{t0, t1}, fluid_lsoda(lopt));
    if (s.success) return s.final_state();
    FluidStiffOptions sopt;
    sopt.rtol = lopt.rtol;
    sopt.atol = lopt.atol;
    // A Rosenbrock failure is raised, not swallowed: both integrators refusing
    // the same leg is a statement about the model, and answering it with the
    // partial trajectory of either one would be a fabricated fixed point.
    const OdeSolution<double> ss = fluid_ode_solve_stiff(f, t0, t1, y0, sopt);
    return ss.final_state();
}

/**
 * The same retry over a whole output grid, for the callers that ask LSODA for a
 * trajectory rather than an endpoint.
 *
 * `lsoda_integrate` stops at the FIRST failed interval and returns the grid it
 * reached, so a caller that reads `s.y` without reading `s.success` silently
 * publishes a trajectory that stops short of the horizon. On failure the grid is
 * re-walked one leg at a time through `fluid_integrate_leg`, which restarts the
 * integrator at every output point -- a small accuracy cost paid only where the
 * single-call integration had already given up.
 */
inline LsodaSolution fluid_integrate_grid(
    const std::function<void(double, const double*, double*)>& f,
    const std::vector<double>& y0, const std::vector<double>& grid,
    const LsodaOptions& lopt) {
    LsodaSolution s = lsoda_integrate(f, y0, grid, fluid_lsoda(lopt));
    if (s.success) return s;
    LsodaSolution out;
    out.t.assign(1, grid[0]);
    out.y.assign(1, y0);
    for (std::size_t j = 1; j < grid.size(); ++j) {
        out.y.push_back(grid[j] == grid[j - 1]
                            ? out.y.back()
                            : fluid_integrate_leg(f, grid[j - 1], grid[j], out.y.back(), lopt));
        out.t.push_back(grid[j]);
    }
    out.method = "rosenbrock4";  // the arm that finished it, not the one that gave up
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_STIFF_H
