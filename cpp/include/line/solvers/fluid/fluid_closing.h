/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_CLOSING_H
#define LINE_SOLVERS_FLUID_FLUID_CLOSING_H

/**
 * Port of `solver_fluid_initsol.m`, and of the entry point of
 * `solver_fluid_closing.m` that consumes it.
 *
 * WHY THE INITIAL CONDITION IS NOT AN IMPLEMENTATION DETAIL. The fluid answer
 * is a fixed point of the drift, and a drift that is not convex has more than
 * one: a closed model whose queue can either drain or saturate settles wherever
 * it was pushed from. The integration therefore SELECTS a fixed point, and what
 * selects it is y0. The reference does not guess: it decodes the model's own
 * initial state, `sn.state`, into the fluid coordinates, so the ODE starts
 * where the model says the system starts. Reproducing that exactly is the whole
 * point of this file; a start that merely conserves the population is a
 * different question with a possibly different answer.
 *
 * WHERE `sn.state` WENT. This port's NetworkStruct carries no state field, so
 * there is nothing to read. What the reference reads, however, is the state
 * `Network.initDefault` wrote: every closed class at its reference station,
 * one job per class in service at each Source, everything else empty, and every
 * started job in PHASE ONE -- `State.fromMarginalAndStarted` writes
 * `init(1) = si(r)` and never enumerates the phase assignment.
 *
 * THE ONE RULE THAT IS EASY TO MISS. Outside a Source, phase one of a block
 * does NOT get `kir(r,1)`: it gets `nir(r) - sum_{k>=2} kir(r,k)`, the jobs in
 * service in phase one PLUS everyone waiting in the buffer. The fluid state has
 * no waiting room -- a queue is one mass per phase -- so the buffer has to land
 * somewhere, and the reference restarts it in phase one. On the state
 * `initDefault` writes every `kir(r,k>=2)` is zero, so that term returns the
 * station's whole population to phase one and the decode collapses to the
 * closed form `detail::fluid_default_initsol` computes directly. The rule is
 * recorded because it is what makes the two the same vector, not because two
 * of them are kept.
 *
 * WHAT THIS FILE DOES NOT REDO. The drift, the restarting integration and the
 * Q/U/R/T extraction of `solver_fluid_closing.m` are already ported in
 * `fluid_odes.h` and `solver_fluid.h` (`fluid_closing_metrics`,
 * `detail::fluid_dispatch`). `solver_fluid_closing` below is the reference
 * function's SHAPE -- seed the initial condition, gate the method, integrate --
 * over that machinery, so there is one copy of the metric rules and not two.
 * The reference's second output, the expanded `state` cell, is not ported: no
 * caller in `solver_fluid_analyzer.m` reads it.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"

namespace line {
namespace fluid {

namespace detail {

/**
 * The disciplines `solver_fluid_initsol.m` knows how to decode.
 *
 * The reference's switch has two arms and an error arm, and the error arm is
 * not a gap to be filled by falling through to the first arm: a discipline
 * outside this list either holds its buffer in an encoding whose waiting jobs
 * cannot be attributed to a phase (POLLING carries a controller, SRPT a
 * remaining-work order) or is preemptive-resume, where a waiting job's phase is
 * recorded per job and restarting it in phase one throws that away. Both would
 * start the drift somewhere the model never is.
 *
 * GPS IS DECODED HERE AND STILL REFUSED BY MOST METHODS, which is not a
 * contradiction: its state encoding folds into one mass per phase exactly like
 * PS's, so `solver_fluid_initsol.m` lists it, while its capacity SHARE needs the
 * backlog probability that only `minnormal` supplies. The refusal therefore
 * belongs to the featset gate (`fluid_feature_set`), per method, and not here.
 *
 * Checked BEFORE the state row is decoded, which the reference does after; for
 * an accepted discipline the two are the same, and for a refused one this
 * reports the discipline rather than whatever the encoder happens to say first.
 */
inline void fluid_initsol_check_sched(lang::SchedStrategy s, std::size_t ist) {
    switch (s) {
        case lang::SchedStrategy::EXT:
        case lang::SchedStrategy::FCFS:
        case lang::SchedStrategy::SIRO:
        case lang::SchedStrategy::PS:
        case lang::SchedStrategy::INF:
        case lang::SchedStrategy::DPS:
        case lang::SchedStrategy::GPS:
        case lang::SchedStrategy::HOL:
        case lang::SchedStrategy::LCFS:
        case lang::SchedStrategy::LCFSPR:
            return;
        default:
            break;
    }
    throw UnsupportedError(
        std::string("solver_fluid_initsol: station ") + std::to_string(ist) + " is scheduled '" +
        lang::sched_to_text(s) +
        "', whose state encoding this port cannot fold into one fluid mass per service phase; "
        "only ext, fcfs, siro, ps, inf, dps, gps, hol, lcfs and lcfspr are decoded");
}

}  // namespace detail

/**
 * Port of `solver_fluid_initsol.m`: the ODE's initial condition, in the layout
 * the drift indexes.
 *
 * THIS IS THE DISCIPLINE GATE OVER `detail::fluid_default_initsol`, AND NOT A
 * SECOND CONSTRUCTION OF THE SAME VECTOR. There used to be two: this one
 * synthesized the `initDefault` marginal, pushed it through `from_marginal`,
 * took the FIRST row and decoded it back through `to_marginal`, while the
 * analyzer's own default wrote the population straight into the reference
 * station's phase-one entry. They agreed on a one-phase model and disagreed on
 * every other, so a multi-phase TRANSIENT was integrated from one of two
 * different initial conditions depending on which entry point was called. The
 * closed form is the correct one, for two reasons:
 *
 *   THE ROW WAS THE WRONG ROW. `from_marginal` enumerates every phase
 *   assignment of the jobs in service in `multichoose` order, whose first row
 *   is `[0,...,0,n]` -- all of the mass in the LAST phase. The reference never
 *   enumerates: `State.fromMarginalAndStarted` writes `init(1) = si(r)` and so
 *   always starts its jobs in phase ONE, which is also what its
 *   `space(end:-1:1,:)` reordering exists to guarantee. Seeding an Erlang-2
 *   station's whole population in phase two starts the drift half a service
 *   ahead of where the model says the system is.
 *
 *   AND IT WAS PAID FOR IN FACTORIALS. The FCFS branch of `from_marginal`
 *   enumerates the buffer PERMUTATIONS, so decoding the initial state of a
 *   closed station holding N jobs cost O(N!) before the first integration step.
 *
 * What the decode contributes that the closed form does not is the by-name
 * refusal above, which is checked here per station and is the reference's own
 * error arm; it is kept, and is the reason this wrapper exists at all.
 */
template <class T>
std::vector<double> fluid_initsol(const qn::NetworkStruct<T>& sn, const FluidLayout& L) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    for (std::size_t i = 0; i < M; ++i) {
        // A station with no enabled class occupies no fluid coordinates, which
        // is how the reference's `isnan(ist)` skip comes out here.
        bool any = false;
        for (std::size_t r = 0; r < K; ++r) any = any || L.enabled[i][r];
        if (any) detail::fluid_initsol_check_sched(sn.stations[i].sched, i + 1);
    }
    return detail::fluid_default_initsol(sn, L);
}

/** The same, for a caller that has not built the layout itself. */
template <class T>
std::vector<double> fluid_initsol(const qn::NetworkStruct<T>& sn) {
    return fluid_initsol(sn, fluid_layout(sn));
}

/**
 * Port of `solver_fluid_closing.m`: the closing family's entry point.
 *
 * IT RETURNS THE UNCORRECTED TABLE, as the reference function does. The
 * utilization and response time that reach a user go through
 * `fluid_analyzer_correct`, which `solver_fluid_analyzer.m` applies AFTER the
 * method switch and to every branch alike; applying it here as well would
 * either double it or fork it. Callers who want the analyzer's answer call
 * `solver_fluid`, which is that function.
 *
 * The method is gated rather than forwarded because the analyzer's switch sends
 * only these names here; `matrix` and `pnorm` are a different drift and
 * `mfq`, `diffusion` and `rmf` are different solvers entirely, and answering
 * for them under this name would report one method's number as another's.
 */
template <class T>
FluidSolution solver_fluid_closing(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "solver_fluid_closing: the fluid solver integrates its drift with LSODA, whose "
            "coefficients assume double precision; rerun with --arith double");

    std::string m = opt.method;
    if (m.size() > 6 && m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    // `default` at THIS entry point means the closing drift: the analyzer's own
    // default resolves to the matrix method, and a caller who wants that calls
    // `solver_fluid`.
    if (m == "default") m = "closing";
    if (!(m == "closing" || m == "statedep" || m == "softmin" || m == "tbi"))
        throw UnsupportedError("solver_fluid_closing: the '" + opt.method +
                               "' method is not part of the closing family; 'matrix' and 'pnorm' "
                               "are solved by solver_fluid_matrix and 'mfq', 'diffusion' and 'rmf' "
                               "by their own solvers, all reachable through solver_fluid");

    FluidOptions o = opt;
    o.method = m;
    // The reference fills `options.init_sol` in the analyzer, before the switch,
    // and re-fills it after every phase refitting; an empty one here means the
    // caller has not overridden it, not that the drift may start anywhere.
    if (o.init_sol.empty()) o.init_sol = fluid_initsol(sn, fluid_layout(sn));
    return detail::fluid_dispatch(sn, o);
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_CLOSING_H
