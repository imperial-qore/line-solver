/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SSA_DISPATCH_H
#define LINE_SOLVERS_SSA_SSA_DISPATCH_H

/**
 * @file
 * @ingroup line_solvers
 * The SolverSSA entry surface: a port of `@@SolverSSA/runAnalyzer.m`'s method
 * whitelist, of `solver_ssa_analyzer.m`'s eligibility gate (`isNrmEligible` and
 * the per-feature `*NrmOK` predicates) and of `solver_ssa_analyzer_nrm.m`'s own
 * scheduling validation.
 *
 * EVERY METHOD THE REFERENCE OFFERS IS REACHED FROM HERE. `SolverSSA.m` line 55
 * lists `default`, `ssa`, `serial`, `para`/`parallel` and `nrm`; the JAR's
 * `SolverSSA.listValidMethods` and native Python's list the same set. All five
 * resolve to one of the three ported engines: the NRM of `solver_ssa_nrm.h`,
 * the event-driven engine of `solver_ssa_serial.h`, and the replica mean of
 * `solver_ssa_parallel.h`. `para`/`parallel` prefers the NRM whenever the model
 * is eligible and replicates the serial engine otherwise, exactly as the
 * reference does; `ssa` is the reference's alias for `serial`, NOT for the NRM;
 * `default` is the reference's ladder -- the NRM when eligible, the serial
 * engine otherwise -- as `solver_ssa` sets out below.
 * `solver_ssa` returns the metric table alone, so a caller who needs the sample
 * path, the per-replica tables or their standard errors calls the engine's own
 * entry.
 *
 * WITHIN `nrm` THERE ARE TWO ENGINES, and both are reached from here. The
 * reference picks between them on `options.config.state_space_gen`
 * (`solver_ssa_analyzer_nrm.m` lines 52-68): `none` and `default` take the
 * plain engine of `solver_ssa_nrm.m`, which integrates the metrics along the
 * sample path, and any other value takes the tabulating engine of
 * `solver_ssa_nrm_space.m`, which records the distinct states visited and forms
 * the means as `pi * A`. `SsaOptions::state_space_gen` is that switch, and the
 * eligibility gate below runs before it, exactly as the reference's does: it is
 * the NRM's gate, not one engine's. The space engine then applies its OWN
 * further refusals (an open model, a phase-type service, a space above the cap),
 * each by name, in `NrmSpaceEngine::check`.
 *
 * WHY THE GATE IS SO LONG, AND WHAT IT IS FOR. It is `isNrmEligible`
 * (`solver_ssa_analyzer.m` lines 365-386) and it answers ONE question: can the
 * NRM run this model? Under `method='nrm'` a no is a refusal BY NAME, because a
 * caller who spelled out the estimator asked for that estimator; under
 * `default` and `parallel` the same no is a fallback to the serial engine, as
 * the reference's is. The two share a body (`raise`), so they cannot drift.
 *
 * WHAT THE GATE MEANS SINCE THE FALLBACK WAS RESTORED. It is the NRM's reach and
 * NOT SolverSSA's: everything it rejects, the serial engine runs (a Fork, a
 * BLOCKING region rule, global dependence), so a rejection under `default`
 * costs the caller the NRM's SPEED and not the answer. Under `method='nrm'` it
 * is still a refusal by name, because an estimator that was asked for by name
 * must be the one that runs. What the NRM lacks falls in two groups:
 *
 *   NOT REPRESENTABLE  the struct has no field at all, so a model using the
 *                      feature cannot be built (balking, reneging patience,
 *                      SelfLoopingClass -- `JobClassType` is OPEN or CLOSED
 *                      only). WRROBIN / JSQ / SQ / RL routing belongs here too
 *                      by a different route: the enumerators exist, but
 *                      `NetworkStruct::refresh_routing` refuses those
 *                      state-dependent strategies when the struct is built, so
 *                      the gate's own test is reached only by a caller who
 *                      filled `routing` by hand. RROBIN is the exception: the
 *                      refresh now EXPANDS it uniformly, for QNA/MNA, which
 *                      recover the determinism from the split degree. An
 *                      RROBIN model therefore reaches this gate for real, and
 *                      the refusal below is what keeps the NRM from simulating
 *                      it as random routing.
 *   REFUSED IN EVERY CODEBASE'S NRM   Fork/Join, the blocking region rules
 *                      BAS / BBS / RSRD, and global (Whittle) dependence; the
 *                      reference's `isNrmEligible` excludes the same set.
 *                      Retrial orbits (`retrialparam`) and G-network signals
 *                      (`issignal` and friends) are refused downstream, by the
 *                      engines that would have to simulate them.
 *
 * The sub-engines that have LEFT the refusal list, all in `solver_ssa_nrm.h`
 * unless named: the *PRIO shares on 2026-08-15; the stochastic Petri net path
 * on 2026-09-13 (`solver_ssa_nrm_spn.h`, which still refuses the four firing
 * laws `spn_nrm_supported` names); and on 2026-09-24 class- and joint-dependent
 * scaling (`cdfac`), finite capacity regions under DROP and WAITQ
 * (`fcr_outcome` / `fcr_release_cascade`), PAS / OI (`pas_rate` /
 * `pas_depart`), POLLING (`PollCtrl` and one switchover reaction per station)
 * and the Cache (`cache_fire`, drawing each READ from `after_event_cache`). The
 * tabulating engine of `solver_ssa_nrm_space.h` does not carry those five and
 * refuses them itself, by name.
 *
 * Nothing falls through silently: a model with a construct the reaction table
 * does not model would run to completion and return numbers that are simply
 * not the model's, which is why the gate stays even where the serial engine can
 * answer.
 */

#include <string>
#include <vector>

#include "line/api/sn/sn_pn_firing_rates.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/ssa/solver_ssa_nrm.h"
#include "line/solvers/ssa/solver_ssa_nrm_space.h"
#include "line/solvers/ssa/solver_ssa_nrm_spn.h"
#include "line/solvers/ssa/solver_ssa_parallel.h"
#include "line/solvers/ssa/solver_ssa_serial.h"
#include "line/solvers/ssa/ssa_types.h"
#include "line/util/error.h"

namespace line {
namespace ssa {

namespace detail {

/**
 * The scheduling whitelist of `solver_ssa_analyzer_nrm.m` lines 15-26,
 * restricted to what this port's rate laws cover.
 *
 * The three priority-weighted shares PSPRIO / DPSPRIO / GPSPRIO joined the
 * whitelist on 2026-08-15, and PAS / OI and POLLING on 2026-09-24, when their
 * sub-engines (the pass-and-swap list and the polling controller) were ported.
 *
 * `raise` is what makes this both the REFUSAL and the ELIGIBILITY test, as the
 * reference's single `isNrmEligible` is: `raise=true` names the offending
 * station and throws, `raise=false` answers false. One body, so the message a
 * caller of `method='nrm'` reads and the predicate `default` and `parallel`
 * branch on can never disagree.
 */
template <class T>
bool ssa_check_scheduling(const qn::NetworkStruct<T>& sn, bool raise = true) {
    using lang::SchedStrategy;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        switch (s) {
            case SchedStrategy::INF:
            case SchedStrategy::EXT:
            case SchedStrategy::PS:
            case SchedStrategy::LPS:
            case SchedStrategy::DPS:
            case SchedStrategy::GPS:
            case SchedStrategy::FCFS:
            case SchedStrategy::LCFS:
            case SchedStrategy::SIRO:
            case SchedStrategy::HOL:
            case SchedStrategy::SEPT:
            case SchedStrategy::LEPT:
            case SchedStrategy::LCFSPR:
            case SchedStrategy::PSPRIO:
            case SchedStrategy::DPSPRIO:
            case SchedStrategy::GPSPRIO:
            // Pass-and-swap and order-independent: served since 2026-09-24 through
            // the ordered-list sub-engine (`pas_rate` / `pas_depart`).
            case SchedStrategy::PAS:
            case SchedStrategy::OI:
            // Polling: served since 2026-09-24 through the controller sub-engine
            // (`PollCtrl`, one switchover reaction per polling station).
            case SchedStrategy::POLLING:
                continue;
            default:
                if (!raise) return false;
                throw UnsupportedError("SolverSSA(method='nrm'): the scheduling policy '" +
                                       std::string(lang::sched_to_text(s)) + "' at station '" +
                                       sn.stations[i].name +
                                       "' is not supported by the NRM in any codebase");
        }
    }
    return true;
}

/**
 * Node kinds the NRM reaction grid cannot express in this port.
 *
 * A Place and a Transition are NOT among them since 2026-09-13:
 * `solver_ssa_nrm_spn.h` ports the reference's SPN builder and run loop, so a
 * net carrying them is routed to that engine instead of refused. What it can
 * and cannot read is `spn_nrm_supported`, asked separately below, because the
 * reference asks it separately too -- its four conditions are about the FIRING
 * LAWS of the net, not about the presence of a node kind.
 */
template <class T>
bool ssa_check_nodes(const qn::NetworkStruct<T>& sn, bool raise = true) {
    for (const qn::NodeDef& nd : sn.nodes) {
        switch (nd.nodetype) {
            // Cache: served since 2026-09-24 through the cache-read reaction
            // (`NrmEngine::cache_fire`).
            case qn::NodeType::Fork:
            case qn::NodeType::Join:
                if (!raise) return false;
                throw UnsupportedError("SolverSSA(method='nrm'): node '" + nd.name +
                                       "' makes this a fork-join model, which the NRM does not "
                                       "handle in any codebase (isNrmEligible excludes it); the "
                                       "reference falls back to the serial engine, whose own "
                                       "fork handler (sn.fjsync / State.afterFJEvent) is not "
                                       "ported either");
            default:
                break;
        }
    }
    return true;
}

/**
 * Routing strategies this port resolves, plus what the model layer itself
 * cannot represent.
 *
 * RROBIN, WRROBIN, JSQ and SQ are resolved AT FIRING TIME by the NRM engine
 * (`build_state_dependent_dest` / `resolve_state_dependent_dest`) from the
 * declared out-arcs and the live population, which is what the reference does
 * and what makes a round-robin dispatcher a dispatcher rather than a coin.
 * `refresh_routing` still expands them into a probability split, because a
 * matrix solver has nothing else to read; the engine ignores that split for
 * these nodes and walks the arcs itself.
 *
 * SDR is the one that stays refused: its routing is a function the struct does
 * not carry, so there is nothing here to evaluate.
 */
template <class T>
bool ssa_check_routing(const qn::NetworkStruct<T>& sn, bool raise = true) {
    for (const qn::NodeDef& nd : sn.nodes)
        for (std::size_t r = 0; r < nd.routing.size(); ++r) {
            const qn::RoutingStrategy rs = nd.routing[r];
            if (rs == qn::RoutingStrategy::PROB || rs == qn::RoutingStrategy::RAND ||
                rs == qn::RoutingStrategy::DISABLED || rs == qn::RoutingStrategy::RROBIN ||
                rs == qn::RoutingStrategy::WRROBIN || rs == qn::RoutingStrategy::JSQ ||
                rs == qn::RoutingStrategy::SQ)
                continue;
            if (!raise) return false;
            throw UnsupportedError(
                "SolverSSA(method='nrm'): node '" + nd.name + "' routes class '" +
                sn.classes[r].name + "' by '" + std::string(lang::routing_to_text(rs)) +
                "'. Its destination is a function the C++ NetworkStruct does not carry, so "
                "there is nothing to resolve at firing time");
        }
    return true;
}

/**
 * `phaseNrmOK` (solver_ssa_analyzer.m lines 407-...): where a non-exponential
 * service process may sit.
 *
 * The INF / PS family expands exactly, because every job present is in service
 * and the class share splits across the phases in the ratio kir/nir. The
 * non-preemptive buffered family expands through the auxiliary in-service
 * multiset. EXT is deliberately excluded and the reason is worth repeating: a
 * phase-type ARRIVAL process is not a service law, and expanding it fires one
 * arrival per PHASE instead of one per RENEWAL, so an Erlang-2 source doubles
 * lambda. LCFSPR is excluded because preempt-resume would have to remember the
 * preempted job's phase.
 */
template <class T>
bool ssa_check_phases(const qn::NetworkStruct<T>& sn, bool raise = true) {
    using lang::SchedStrategy;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        const bool exact = s == SchedStrategy::INF || s == SchedStrategy::PS ||
                           s == SchedStrategy::LPS || s == SchedStrategy::DPS ||
                           s == SchedStrategy::GPS || s == SchedStrategy::FCFS ||
                           s == SchedStrategy::LCFS || s == SchedStrategy::SIRO ||
                           s == SchedStrategy::HOL || s == SchedStrategy::SEPT ||
                           s == SchedStrategy::LEPT;
        // A MAP/MMPP2 service at FCFS chains its draws through a memory phase (and, with c > 1,
        // the landing phase of each draw, `ctmc_mapc`), which the in-service multiset has no
        // column for: the NRM would restart every draw from the entry law.
        if (s == SchedStrategy::FCFS)
            for (std::size_t r = 0; r < sn.nclasses; ++r) {
                if (sn.disabled[i][r]) continue;
                const lang::ProcessType pt = sn.service[i][r].type;
                if (pt != lang::ProcessType::MAP && pt != lang::ProcessType::MMPP2) continue;
                if (!raise) return false;
                throw UnsupportedError(
                    "SolverSSA(method='nrm'): class '" + sn.classes[r].name + "' has MAP service "
                    "at FCFS station '" + sn.stations[i].name + "', whose draws are chained "
                    "through a memory phase the NRM does not carry; ask for 'default', which "
                    "falls back to the serial engine");
            }
        if (exact) continue;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (sn.disabled[i][r]) continue;
            const lang::ProcessType pt = sn.service[i][r].type;
            if (pt == lang::ProcessType::EXP || pt == lang::ProcessType::DISABLED ||
                pt == lang::ProcessType::IMMEDIATE)
                continue;
            if (!raise) return false;
            throw UnsupportedError(
                "SolverSSA(method='nrm'): class '" + sn.classes[r].name + "' has non-exponential "
                "service at station '" + sn.stations[i].name + "', whose '" +
                std::string(lang::sched_to_text(s)) +
                "' discipline the NRM phase expansion does not cover; ask for 'default', which "
                "falls back to the serial engine as the reference does");
        }
    }
    return true;
}

/**
 * Global (Whittle) dependence, refused because the NRM does not SCALE by it.
 *
 * Class- and joint-dependent scaling are served since 2026-09-24: `cdfac` in
 * `solver_ssa_nrm.h` evaluates eta_i(n) .* beta_{i,r}(n) on the station's
 * per-class population at firing time, as the reference does, and utilization
 * is normalized by the declared peaks. The reference's gate (`ssa_nrm_guards.gd`)
 * excludes only the global handle, and so does this one.
 */
template <class T>
bool ssa_check_cdscaling(const qn::NetworkStruct<T>& sn, bool raise = true) {
    // The NRM's propensities receive the per-station population slice, not the
    // whole population matrix a global handle reads.
    if (sn.gdscaling) {
        if (!raise) return false;
        throw UnsupportedError(
            "SolverSSA(method='nrm'): the model declares a global dependence "
            "(setGlobalDependence). The NRM builds its propensities from the per-station "
            "population slice and never sees the whole population matrix the handle reads, so "
            "the sample path would run at the unscaled rates; ask for method='serial', which "
            "applies it");
    }
    return true;
}

/**
 * Finite capacity regions: served since 2026-09-24, under DROP and WAITQ.
 *
 * `solver_ssa_nrm.h` ports `fcrPrecompute` / `fcrRefusingRegion` /
 * `fcrReleaseCascade`: an admission gate on every arrival into a region and a
 * head-of-line release cascade from a per-region FIFO after every firing. What
 * stays refused is a BLOCKING region rule (BAS, BBS, RSRD), which no codebase's
 * NRM carries.
 */
template <class T>
bool ssa_check_regions(const qn::NetworkStruct<T>& sn, bool raise = true) {
    for (const typename qn::NetworkStruct<T>::Region& rg : sn.regions)
        for (std::size_t r = 0; r < rg.rule.size(); ++r) {
            const qn::DropStrategy d = rg.rule[r];
            if (d != qn::DropStrategy::BAS && d != qn::DropStrategy::BBS &&
                d != qn::DropStrategy::RSRD)
                continue;
            if (!raise) return false;
            throw UnsupportedError(
                "SolverSSA(method='nrm'): finite capacity region '" + rg.name +
                "' declares a blocking rule for class '" + sn.classes[r].name +
                "'. The NRM region gate implements DROP and WAITQ only, as the reference's does");
        }
    return true;
}

/**
 * Synchronous calls (REPLY signals), `ssa_nrm_guards.reply`: refused by the NRM in every codebase.
 *
 * A calling job leaves its caller and keeps the server there until the reply returns. The
 * serial engine carries that hold in the reply block of the local row (`after_event`), while
 * the NRM's reactions see per-class counts only and have no held-server counter, so the caller
 * would free its server at once. `default` and `parallel` therefore take the serial engine,
 * and an explicit `nrm` falls back to it with a warning, as the reference does.
 */
template <class T>
bool ssa_check_reply(const qn::NetworkStruct<T>& sn, bool raise = true) {
    for (std::size_t ind = 0; ind < sn.replyblock.size(); ++ind)
        for (std::size_t k = 0; k < sn.replyblock[ind].size(); ++k) {
            if (!sn.replyblock[ind][k]) continue;
            if (!raise) return false;
            throw UnsupportedError(
                "SolverSSA(method='nrm'): class '" + sn.classes[k].name + "' makes a synchronous "
                "call (REPLY signal) from node '" + sn.nodes[ind].name + "'. The NRM has no "
                "held-server counter for a pending reply; use method='serial'");
        }
    return true;
}

/**
 * `isNrmEligible` (`solver_ssa_analyzer.m` lines 365-386): can the NRM run this
 * model at all?
 *
 * It is the SAME tests the refusal above takes, asked with `raise=false`,
 * so the predicate cannot drift from the message. It is deliberately the C++
 * NRM's reach and not the reference's: this port's NRM covers less, and a
 * `default` that consulted MATLAB's wider predicate would send a model to a
 * refusal the serial engine can actually answer.
 */
template <class T>
bool ssa_nrm_eligible(const qn::NetworkStruct<T>& sn) {
    // Immediate feedback keeps the fed-back job on the server it just used, which
    // the serial engine expresses through `after_event`'s `no_promote` arc. The
    // NRM has no such arc: its reactions see per-class counts, not which job holds
    // the server, so a self-loop there is a class switch with re-queueing. An
    // explicit `method='nrm'` still RUNS the model with that approximation,
    // warning in as many words (`ssa_nrm_supports`, the GATE, deliberately does
    // not test this); this only keeps `default` and `parallel` from PREFERRING an
    // approximation over the exact serial engine.
    if (sn.has_immediate_feedback()) return false;
    return ssa_check_nodes(sn, false) && ssa_check_scheduling(sn, false) &&
           ssa_check_routing(sn, false) && ssa_check_phases(sn, false) &&
           ssa_check_cdscaling(sn, false) && ssa_check_regions(sn, false) &&
           ssa_check_reply(sn, false) && spn_nrm_supported(sn, false);
}

/**
 * The same tests as a SENTENCE, so a report can say why `nrm` is withheld.
 *
 * `ssa_nrm_eligible` above answers the dispatch's question (may I PREFER the
 * NRM?) and nothing answered the gate's (may I OFFER it?), so `ssa.nrm` was
 * reported runnable on every model and an explicit request then raised. The
 * catch is the ADAPTER between those two callers and not error suppression: the
 * checks already carry the exact wording a user should see, and re-deriving it
 * here is how the message and the predicate drift apart.
 */
template <class T>
std::string ssa_nrm_supports(const qn::NetworkStruct<T>& sn) {
    try {
        ssa_check_nodes(sn, true);
        ssa_check_scheduling(sn, true);
        ssa_check_routing(sn, true);
        ssa_check_phases(sn, true);
        ssa_check_cdscaling(sn, true);
        ssa_check_regions(sn, true);
        ssa_check_reply(sn, true);
        spn_nrm_supported(sn, true);
    } catch (const UnsupportedError& e) {
        return e.what();
    }
    return "";
}

/**
 * The tabulating engine's knobs from the caller's.
 *
 * `state_max` keeps its own default because `SsaOptions` carries no cap: a
 * caller who wants a different one calls `solver_ssa_nrm_space_analyzer`
 * directly, which is also where the tabulated path and the propensity table
 * survive rather than being reduced to the metrics.
 */
inline SsaNrmSpaceOptions ssa_space_options(const SsaOptions& o) {
    SsaNrmSpaceOptions s;
    static_cast<SsaOptions&>(s) = o;
    return s;
}

/**
 * The serial and replicated engines' knobs from the caller's.
 *
 * `cutoff`, `state_max`, `nreplicas` and `eventcache` keep their
 * own defaults for the same reason `state_max` does above: `SsaOptions` has no
 * field for them, and `nreplicas = 8` is `SolverOptions('SSA')`'s own default
 * rather than a number chosen here. A caller who wants a different R calls
 * `solver_ssa_parallel` directly, which is also where the R per-replica tables
 * and their standard errors survive rather than being reduced to the mean.
 */
inline SsaSerialOptions ssa_serial_options(const SsaOptions& o) {
    SsaSerialOptions s;
    static_cast<SsaOptions&>(s) = o;
    return s;
}

inline SsaParallelOptions ssa_parallel_options(const SsaOptions& o) {
    SsaParallelOptions p;
    static_cast<SsaOptions&>(p) = o;
    return p;
}

/**
 * One serial run, keeping the cache write-back the metric table does not carry.
 *
 * Written once because `solver_ssa` reaches the serial engine from two arms
 * (the `default` fallback and the named method), and a caller that got its
 * cache shares from one arm and not the other would report the offered split
 * for the same model under a different spelling of the same request.
 */
template <class T>
SsaSolution ssa_serial_avg(const qn::NetworkStruct<T>& sn, const SsaOptions& opt,
                           std::vector<SsaCacheRatio>* cache) {
    SsaSerialSolution<T> s = solver_ssa_serial_analyzer(sn, ssa_serial_options(opt));
    if (cache) *cache = s.cache;
    return s.avg;
}

/**
 * `@@SolverSSA/runAnalyzer.m:227`: a Place counts TOKENS, an engine counts
 * FIRING EVENTS, and the two are the same number only where every arc has
 * multiplicity one and no immediate transition moves the tokens.
 *
 * WHAT THE ENGINES HAND OVER. Every SPN engine accumulates a Place's
 * throughput as the summed propensity of the TIMED modes consuming from it, so
 * an arc of weight 2 counts one departure where two tokens left, and an
 * IMMEDIATE consumer counts none at all -- a firing that takes zero time is not
 * a timed event and no engine, here or in the reference, ever sees it as one.
 * A Place drained only by an immediate transition therefore arrives here with
 * TN = 0, which then makes RespT 0 and leaves the neighbouring Place's ArvR 0
 * once the caller derives it from this table.
 *
 * `sn_pn_avg_rates` is the reference's repair and recovers the missing rates
 * from the flow-balance equations, which is why it runs at the SOLVER level and
 * on the ORIGINAL struct: both engines fold a fork-join tagging back before
 * returning, so the table is already in the caller's own coordinates here, as
 * it is after `solver_tr_fjtag_analyzer`'s lift in the reference. A model with
 * no Place, or one whose firing rates cannot be recovered, is returned
 * untouched.
 *
 * AN IS LEFT TO THE CALLER, exactly as `[TN,~,RN] = sn_pn_avg_rates(...)` leaves
 * it in the reference: `SsaSolution` carries no arrival-rate table, and every
 * consumer derives one with `sn_get_arvr_from_tput` from the TN repaired here.
 */
template <class T>
SsaSolution ssa_pn_token_rates(const qn::NetworkStruct<T>& sn, const SsaSolution& in) {
    bool has_place = false;
    for (std::size_t a = 0; a < sn.nodes.size() && !has_place; ++a)
        if (sn.nodes[a].nodetype == qn::NodeType::Place) has_place = true;
    if (!has_place || in.TN.rows() == 0) return in;

    const std::size_t M = in.TN.rows(), K = in.TN.cols();
    Matrix<T> QN(M, K), TN(M, K), RN(M, K);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            QN(i, k) = num_traits<T>::from_double(in.QN(i, k));
            TN(i, k) = num_traits<T>::from_double(in.TN(i, k));
            RN(i, k) = num_traits<T>::from_double(in.RN(i, k));
        }
    const api::SnPnAvgRates<T> pn = api::sn_pn_avg_rates(sn, QN, TN, Matrix<T>(), RN);

    SsaSolution out = in;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            out.TN(i, k) = num_traits<T>::to_double(pn.TN(i, k));
            out.RN(i, k) = num_traits<T>::to_double(pn.RN(i, k));
        }
    return out;
}

}  // namespace detail

/**
 * Port of `SolverSSA.listValidMethods`.
 *
 * Six names for three engines, because the reference spells the same engine
 * more than one way: 'ssa' and 'serial' are the serial trajectory, 'para' and
 * 'parallel' the replicated one, 'nrm' the next-reaction method, and 'default'
 * is the ladder `solver_ssa` walks -- the NRM when it can run the model and the
 * serial engine when it cannot. Every name here is dispatched by `solver_ssa`
 * below, which refuses anything else by name.
 */
inline std::vector<std::string> list_valid_methods() {
    return {"default", "ssa", "serial", "para", "pana", "parallel", "nrm"};
}

/**
 * `solver_ssa_analyzer_nrm.m`: run the NRM and return the metric table.
 *
 * The reference's post-processing (`QN(isnan(QN)) = 0` and the rest) is inside
 * the engine already: it never produces a NaN, because every division is
 * guarded at the point it is taken.
 */
template <class T>
SsaSolution solver_ssa_nrm_analyzer(const qn::NetworkStruct<T>& sn, const SsaOptions& opt,
                                    std::vector<SsaCacheRatio>* cache = nullptr) {
    // `if constexpr`, not a run-time test: the engine reaches `dist_to_map`,
    // whose APH fit static_asserts on transcendental arithmetic, so a Rational
    // instantiation would fail to COMPILE rather than refuse. The gate has to
    // keep the body from being instantiated at all.
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        (void)cache;
        throw UnsupportedError(
            "solver_ssa_nrm: an SSA sample path is generated from exponential clocks, which are "
            "logarithms of uniform draws; there is no exact value to compute and a wider float "
            "carries no information the Monte Carlo error does not swamp. Rerun with --arith "
            "double");
    } else {
        detail::ssa_check_nodes(sn);
        detail::ssa_check_scheduling(sn);
        detail::ssa_check_routing(sn);
        detail::ssa_check_phases(sn);
        detail::ssa_check_cdscaling(sn);
        detail::ssa_check_regions(sn);
        detail::ssa_check_reply(sn);
        detail::spn_nrm_supported(sn);
        // `solver_ssa_nrm.m:174`: a net carrying a Transition is a STOCHASTIC
        // PETRI NET and goes to the SPN builder and run loop, which reads a
        // marking of tokens rather than a grid of jobs in service. The branch is
        // taken before the space/plain switch below because the tabulating
        // engine tabulates the queueing state vector, which an SPN does not
        // have; the reference branches at the same point and for the same
        // reason.
        for (const qn::NodeDef& nd : sn.nodes)
            if (nd.nodetype == qn::NodeType::Transition) {
                NrmSpnEngine<T> spn(sn, opt);
                return spn.run();
            }
        // The reference's engine switch, taken AFTER the gate above because the
        // gate is the NRM's and not one engine's.
        if (opt.state_space_gen != "none" && opt.state_space_gen != "default")
            return solver_ssa_nrm_space_analyzer(sn, detail::ssa_space_options(opt)).avg;
        NrmEngine<T> eng(sn, opt);
        SsaSolution out = eng.run();
        if (cache) *cache = eng.cache();
        return out;
    }
}

/**
 * `solver_ssa_analyzer.m`: choose the method.
 *
 * The ladder is the reference's, in its order:
 *
 *   `default`   the NRM when it is eligible, ELSE the serial engine, and only a
 *               model neither engine runs is refused -- by the NRM's message,
 *               which names the missing sub-engine.
 *   `nrm`       the NRM alone, reached by naming the estimator: a model it
 *               cannot run is refused rather than silently answered by the
 *               other engine, because a caller who spelled out an estimator
 *               asked for that estimator's variance as well as its mean.
 *   `ssa`       the reference's alias for `serial` (line 128), not for the NRM.
 *   `serial`    one run of the event-driven engine.
 *   `para`      the NRM when eligible (the reference prefers one fast run over
 *   `parallel`  R replicated ones, lines 143-157), else the replica mean of
 *               `nreplicas` independent serial runs.
 *
 * `default`'s FALLBACK IS THE REFERENCE'S (lines 66-78) and was restored on
 * 2026-07-31, when the serial engine gained fork-join and the finite capacity
 * regions. It had been held back while the serial engine covered less than the
 * NRM gate rejects, on the ground that a second failure further downstream is
 * less informative than the NRM's own message. That ground is gone for every
 * construct the serial engine now runs, and where it still holds -- a model
 * NEITHER engine covers -- the NRM's message is what the caller reads, because
 * `serial_can_run` is asked BEFORE the fallback is taken rather than after it
 * has failed.
 *
 * `cache`, when given, receives the serial engine's cache write-back -- the
 * realized hit and miss shares of every Cache node, which are a SOLVER RESULT
 * and the only thing that tells a node table apart from the 1/2-1/2 `link()`
 * offers. It is an out-parameter rather than a field of `SsaSolution` because
 * that struct is the metric table the three engines share, and only two of them
 * have a cache to report, the serial engine and (since 2026-09-24) the NRM; the
 * parallel engine averages replicas that carry none.
 */
template <class T>
SsaSolution solver_ssa_engine(const qn::NetworkStruct<T>& sn, const SsaOptions& opt,
                              std::vector<SsaCacheRatio>* cache) {
    // A fed-back job keeps its server, which the reference encodes in
    // `State.afterEvent`'s immfeed self-loop, and only the SERIAL engine has that
    // arc -- here as in the reference. `solver_ssa_analyzer_nrm.m:29` warns in as
    // many words that the NRM "does not model immediate feedback (immfeed);
    // self-loops are treated as class-switching with re-queueing" and then runs,
    // so the name stays offerable and only the approximation is announced.
    // `ssa_nrm_eligible` excludes immediate feedback, so `default` and `parallel`
    // reach the exact serial engine and only an EXPLICIT `nrm` warns.
    if (sn.has_immediate_feedback()) {
        const std::string& mm = opt.method;
        const bool runs_nrm =
            mm == "nrm" || ((mm == "default" || mm == "para" || mm == "pana" || mm == "parallel") &&
                            detail::ssa_nrm_eligible(sn));
        if (runs_nrm)
            std::cerr << "[LINE] Warning: SolverSSA(method=nrm) does not model immediate feedback "
                         "(immfeed); self-loops are treated as class-switching with re-queueing. "
                         "Use method='serial' for immediate feedback."
                      << std::endl;
    }
    const std::string& m = opt.method;
    if (m == "default") {
        // The reference's ladder: the NRM when it can run the model, the serial
        // engine when it cannot, and the NRM's own refusal -- which names the
        // sub-engine that is missing -- when neither can. Asking
        // `serial_can_run` first is what keeps that message: falling through to
        // the serial engine and letting IT fail would report whichever guard it
        // hit, about a model the caller never asked it to run.
        if (detail::ssa_nrm_eligible(sn)) return solver_ssa_nrm_analyzer(sn, opt, cache);
        if (serial_detail::serial_can_run(sn))
            return detail::ssa_serial_avg(sn, opt, cache);
        return solver_ssa_nrm_analyzer(sn, opt, cache);  // raises, by name
    }
    if (m == "nrm") {
        // `solver_ssa_analyzer.m`'s explicit 'nrm' arm: a synchronous call falls back to the
        // serial engine with a warning rather than refusing, in every codebase.
        if (!detail::ssa_check_reply(sn, false)) {
            std::cerr << "[LINE] Warning: SolverSSA: NRM does not support synchronous calls "
                         "(REPLY signals); falling back to the serial method."
                      << std::endl;
            return detail::ssa_serial_avg(sn, opt, cache);
        }
        return solver_ssa_nrm_analyzer(sn, opt, cache);
    }
    if (m == "ssa" || m == "serial") return detail::ssa_serial_avg(sn, opt, cache);
    if (m == "para" || m == "pana" || m == "parallel") {
        // The reference's own preference, `solver_ssa_analyzer.m` lines 143-157:
        // an NRM-eligible model runs once on the NRM rather than R times here.
        if (detail::ssa_nrm_eligible(sn)) return solver_ssa_nrm_analyzer(sn, opt, cache);
        return solver_ssa_parallel_analyzer(sn, detail::ssa_parallel_options(opt)).avg;
    }
    throw UnsupportedError("SolverSSA: '" + m +
                           "' is not a valid method; the reference offers 'default', 'ssa', "
                           "'serial', 'para', 'parallel' and 'nrm', all of which this port "
                           "implements");
}

/**
 * `@@SolverSSA/runAnalyzer` itself: the engine the method selects, then the
 * result assembly the reference does before it hands the table to a caller.
 *
 * Only one step of that assembly belongs here rather than to the caller. A
 * Place is measured in TOKENS and every engine counts FIRING EVENTS, so
 * `ssa_pn_token_rates` rewrites the Place rows of TN and RN exactly where
 * `runAnalyzer.m:227` does; the arrival rates and residence times each consumer
 * still derives for itself, as before. A model with no Place is untouched, so
 * this costs an ordinary queueing network one scan of the node list.
 */
template <class T>
SsaSolution solver_ssa(const qn::NetworkStruct<T>& sn, const SsaOptions& opt,
                       std::vector<SsaCacheRatio>* cache = nullptr) {
    return detail::ssa_pn_token_rates(sn, solver_ssa_engine(sn, opt, cache));
}

/**
 * `CacheMetrics` from the serial engine's cache write-back.
 *
 * The shares are matched to their Cache node BY INDEX, not by position: the
 * base struct walks every Cache the model declares while the engine reports
 * only the nodes it simulated, so pairing the two off in order would write one
 * cache's measurement onto another on any model holding more than one.
 *
 * `residt` becomes `latency` unchanged, NaN included. The reference warns that
 * retrieval latency is not implemented and reports NaN in every codebase, so
 * carrying the NaN IS parity; dropping the field would report "not computed"
 * about a quantity the engine did state.
 */
template <class T>
solvers::CacheMetrics<T> cache_metrics_of_ssa(const qn::NetworkStruct<T>& sn,
                                              const std::vector<SsaCacheRatio>& cache) {
    solvers::CacheMetrics<T> out = solvers::cache_metrics_of(
        sn, std::vector<T>(), std::vector<T>(), std::vector<T>(), std::vector<T>(), Matrix<T>(),
        Matrix<T>(), std::vector<T>());
    for (std::size_t c = 0; c < out.caches.size(); ++c) {
        for (std::size_t j = 0; j < cache.size(); ++j) {
            if (cache[j].node != out.caches[c].node) continue;
            solvers::CacheNodeMetrics<T>& m = out.caches[c];
            m.hitprob.clear();
            m.missprob.clear();
            m.delayedprob.clear();
            m.latency.clear();
            for (std::size_t r = 0; r < cache[j].hitprob.size(); ++r)
                m.hitprob.push_back(num_traits<T>::from_double(cache[j].hitprob[r]));
            for (std::size_t r = 0; r < cache[j].missprob.size(); ++r)
                m.missprob.push_back(num_traits<T>::from_double(cache[j].missprob[r]));
            for (std::size_t r = 0; r < cache[j].delayedprob.size(); ++r)
                m.delayedprob.push_back(num_traits<T>::from_double(cache[j].delayedprob[r]));
            for (std::size_t r = 0; r < cache[j].residt.size(); ++r)
                m.latency.push_back(num_traits<T>::from_double(cache[j].residt[r]));
            break;
        }
    }
    return out;
}

/**
 * The struct with the cache split the SIMULATION MEASURED, visits rebuilt.
 *
 * A cache splits the read stream into a hit stream and a miss stream, and that
 * split IS routing: the visit ratios of everything downstream depend on it. The
 * base struct carries only what `link()` offered, an even share over the hit and
 * miss classes, because the split is a RESULT and cannot be known before the
 * solve. Anything derived from visits after the solve therefore has to be taken
 * on the rewritten struct, not on the base one -- `sn_get_residt_from_respt`
 * reported ResidT = RespT/2 for both classes of tut06_cache_lru_zipf (0.1 and
 * 0.5 against the reference's 0.16475 and 0.17625), which is the even split
 * showing through, not a residence time.
 *
 * The rewrite is the one `da_cacheqn` performs between passes, with the measured
 * shares in place of the analytical ones: the cache row of `rtnodes` is cleared
 * and the hit and miss mass sent to every connected successor, after which
 * `da_recompute_visits_from_rtnodes` rebuilds `visits` and `nodevisits`. A share
 * the engine left undefined (NaN, a class that does not read this cache) leaves
 * that row alone rather than zeroing a routing the model does have.
 */
template <class T>
qn::NetworkStruct<T> sn_with_ssa_cache_split(const qn::NetworkStruct<T>& base,
                                             const std::vector<SsaCacheRatio>& cache) {
    if (cache.empty()) return base;
    qn::NetworkStruct<T> sn = base;
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;
    if (I == 0 || K == 0 || sn.rtnodes.rows() < I * K) return base;
    const T zero = num_traits<T>::from_int(0);
    // The connectivity is read BEFORE any row is rewritten, as in `da_cacheqn`:
    // a cache's own read self-switch is not a downstream node.
    std::vector<std::vector<bool> > conn(I, std::vector<bool>(I, false));
    for (std::size_t a = 0; a < I; ++a)
        for (std::size_t b = 0; b < I; ++b) {
            if (a == b) continue;
            for (std::size_t r = 0; r < K && !conn[a][b]; ++r)
                for (std::size_t s = 0; s < K && !conn[a][b]; ++s)
                    if (sn.rtnodes(a * K + r, b * K + s) > zero) conn[a][b] = true;
        }
    bool rewrote = false;
    for (std::size_t j = 0; j < cache.size(); ++j) {
        const std::size_t ind = cache[j].node;  // 1-based node index
        if (ind == 0 || ind > I) continue;
        typename std::map<std::size_t, qn::CacheParam<T> >::const_iterator ci =
            sn.nodeparam.find(ind);
        if (ci == sn.nodeparam.end()) continue;
        const std::size_t nd = ind - 1;
        for (std::size_t r = 0; r < K; ++r) {
            if (ci->second.hitclass.size() <= r || ci->second.missclass.size() <= r) continue;
            const std::size_t hc = ci->second.hitclass[r], mc = ci->second.missclass[r];
            if (hc == 0 || mc == 0 || hc > K || mc > K) continue;
            if (cache[j].hitprob.size() <= r || cache[j].missprob.size() <= r) continue;
            const double hp = cache[j].hitprob[r], mp = cache[j].missprob[r];
            if (!(hp == hp) || !(mp == mp)) continue;  // NaN: the engine measured none
            for (std::size_t col = 0; col < I * K; ++col) sn.rtnodes(nd * K + r, col) = zero;
            for (std::size_t jnd = 0; jnd < I; ++jnd) {
                if (!conn[nd][jnd]) continue;
                sn.rtnodes(nd * K + r, jnd * K + (hc - 1)) = num_traits<T>::from_double(hp);
                sn.rtnodes(nd * K + r, jnd * K + (mc - 1)) = num_traits<T>::from_double(mp);
            }
            rewrote = true;
        }
    }
    if (!rewrote) return base;
    sn.da_recompute_visits_from_rtnodes();
    return sn;
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SSA_DISPATCH_H
