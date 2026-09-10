/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_SOLVER_AG_RUNNER_H
#define LINE_SOLVERS_AG_SOLVER_AG_RUNNER_H

/**
 * @file solver_ag_runner.h
 * @brief The gates and the dispatch of the agent-based (RCAT) solver.
 *
 * The twin of solver_mam_runner.h for SolverAG: the valid method list, the
 * structural gate that decides whether the RCAT decomposition can represent the
 * model at all, and the conversion every AG method needs before the fixed point.
 *
 * These gates used to live in SolverMAM, which was answering for two unrelated
 * decompositions at once: RCAT decomposes the MODEL into cooperating agents,
 * every other MAM method decomposes its TRAFFIC. They share the shape of the
 * answer and nothing else.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "line/api/sn/sn_nonmarkov_toph.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ag/ag_types.h"
#include "line/solvers/ag/solver_ag.h"
#include "line/util/error.h"

namespace line {
namespace ag {

namespace runner_detail {

/** Every method SolverAG serves. */
inline std::vector<std::string> list_valid_methods() {
    return {"default", "inap", "inapplus", "inapinf", "exact"};
}

/** An unlisted method is refused rather than silently resolved. */
inline void check_method(const std::string& method) {
    const std::vector<std::string> valid = list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) == valid.end()) {
        std::string names;
        for (std::size_t i = 0; i < valid.size(); ++i) {
            names += (i ? ", " : "") + valid[i];
        }
        throw InputError("SolverAG: unknown method '" + method + "'; valid methods are "
                         + names);
    }
}

/**
 * `SolverAG.supportsModelMethod`: the structural gate of the RCAT decomposition,
 * as the REASON it refuses, empty when the model is admissible.
 *
 * RCAT gives every agent a service-phase and an arrival-phase dimension, so any
 * law with a genuine (D0,D1) Markovian representation is admissible; see
 * ag::rcat_supports_process for what is not. A signal is a trigger with no
 * service, so its service entry is never read; only its Source arrival rate is,
 * and that one must stay exponential because the removal is folded into the
 * agent as a scalar rate.
 *
 * A STRING RATHER THAN A THROW, the shape `ba::method_refusal` and
 * `mam::mam_model_method_refusal` already carry: the AUTO report has to ASK the
 * question without raising, so that a pair it offers is a pair the run accepts.
 * `check_model_method` below is the throwing form the analyzer keeps.
 */
template <class T>
std::string method_refusal(const qn::NetworkStruct<T>& L, const std::string& method) {
    using lang::ProcessType;
    for (std::size_t i = 0; i < L.nstations; ++i) {
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            // has_service_law, not disabled: a Join has rates = Inf and no
            // service law, and would otherwise be refused out of hand.
            if (!L.has_service_law(i, r) || !(L.rates(i, r) > num_traits<T>::from_int(0)))
                continue;
            const bool is_signal = r < L.issignal.size() && L.issignal[r];
            const bool is_source = L.stations[i].nodetype == qn::NodeType::Source;
            if (is_signal && !is_source) continue;
            if (is_signal) {
                if (L.service[i][r].type == ProcessType::EXP) continue;
                return std::string(
                    "SolverAG: the " + method +
                    " method needs an exponential signal arrival process (a removal signal "
                    "is folded into the agent as a scalar rate), but station '" +
                    L.stations[i].name + "' class '" + L.classes[r].name + "' is " +
                    lang::process_to_text(L.service[i][r].type) +
                    ". Use SolverMAM (-s mam, method 'dec.source') for such models");
            }
            if (rcat_supports_process(L.service[i][r].type)) continue;
            return std::string(
                "SolverAG: the " + method +
                " method supports processes with a Markovian (D0,D1) representation only "
                "(RCAT builds a phase dimension per agent out of it), but station '" +
                L.stations[i].name + "' class '" + L.classes[r].name + "' is " +
                lang::process_to_text(L.service[i][r].type) +
                ". Use SolverMAM (-s mam, method 'dec.source') for such models");
        }
    }
    // ROUTING MUST BE STATE INDEPENDENT, and the featset gate above CANNOT check
    // this, which is why the explicit loop stays. build_rcat reads sn.rt once and
    // stores it as RcatAction::prob; for RROBIN, WRROBIN, JSQ, SQ or SDR that entry
    // is only a state-independent PLACEHOLDER (SyncEvent::prob in state_events.h
    // says so, and the generator is meant to re-evaluate the split per state), so
    // solving from it answers a DIFFERENT model at a plausible number.
    //
    // WHY THE FEATSET MISSES IT: `used_lang_features` registers a routing feature
    // only inside the Queue/Delay branch and the Router/Dispatcher one. The Source
    // branch deliberately does not -- "MATLAB registers no discipline and no
    // routing for it" -- so a JSQ declared ON A SOURCE, which is the canonical JSQ
    // model, never reaches the DECLARED-vs-USED comparison at all. That omission is
    // faithful to the reference, so it cannot simply be widened here; it means the
    // gate and this loop cover different ground and both are needed.
    for (std::size_t nd = 0; nd < L.nodes.size(); ++nd) {
        for (std::size_t r = 0; r < L.nodes[nd].routing.size(); ++r) {
            const lang::RoutingStrategy rs = L.nodes[nd].routing[r];
            if (rs == lang::RoutingStrategy::PROB || rs == lang::RoutingStrategy::RAND
                || rs == lang::RoutingStrategy::DISABLED) {
                continue;
            }
            return std::string("SolverAG: the ") + method + " method needs state-independent "
                "routing (RCAT folds the split into a fixed action probability), but node '" +
                L.nodes[nd].name + "' routes class '" + L.classes[r].name + "' by " +
                lang::routing_to_text(rs) +
                ", whose share depends on the network state. Use SolverCTMC, SolverSSA or "
                "SolverLDES, which re-evaluate the split per state";
        }
    }

    // build_rcat never reads sn.nservers, so a c-server station would be driven
    // at rho = lambda/mu instead of lambda/(c mu).
    for (std::size_t i = 0; i < L.nstations; ++i) {
        if (std::isfinite(L.stations[i].nservers) && L.stations[i].nservers > 1.0) {
            return std::string(
                "SolverAG: the " + method +
                " method supports single-server stations only (RCAT does not model "
                "sn.nservers), but station '" + L.stations[i].name + "' has " +
                std::to_string(static_cast<long long>(L.stations[i].nservers)) +
                " servers. Use SolverMAM (-s mam, method 'dec.source') for multiserver "
                "models");
        }
    }

    // A FINITE BUFFER IS NOT SOMETHING RCAT CAN CARRY, and it was being answered
    // rather than refused: nothing under solvers/ag reads st.cap or st.classcap,
    // so a capped station was decomposed as an unbounded one and the table
    // reported the UNCONSTRAINED figures (a closed 2-job tandem with cap 1 on the
    // second queue returned the same numbers with and without the cap, 1.09 jobs
    // in a buffer of 1). Lowering the agent's level bound (nlev = njobs(r)+1) to
    // the buffer would not fix it: the top-level boundary is a self-loop, which
    // LOSES the arrival, whereas a closed job refused at a full buffer must BLOCK
    // the upstream departure, and that coupling is exactly the independence RCAT
    // assumes. see _kb/06-solver-catalog.md
    return qn::binding_capacity_reason(std::string("SolverAG"), L);
}

/** The refusal above as the analyzer takes it: raise when there is one. */
template <class T>
void check_model_method(const qn::NetworkStruct<T>& L, const std::string& method) {
    const std::string why = method_refusal(L, method);
    if (!why.empty()) throw UnsupportedError(why);
}

}  // namespace runner_detail

/**
 * Every method SolverAG serves, as the other families expose theirs.
 *
 * `auto_family_methods` asks each family for its own list rather than keeping a
 * second copy, so the "ag" row of the AUTO report is this list and cannot drift
 * from the one `check_method` refuses against.
 */
inline std::vector<std::string> list_valid_methods() {
    return runner_detail::list_valid_methods();
}

/**
 * The gates and the dispatch of SolverAG's runAnalyzer.
 *
 * The conversion is not method-dependent, unlike SolverMAM's: EVERY AG method
 * builds a CTMC per agent out of (D0,D1), so a preserved Det would reach it with
 * no matrix at all and be read back as its mean rate, and a concentrated matrix
 * exponential is not a generator at all.
 */
template <class T>
AgResult<T> solver_ag_solve(const qn::NetworkStruct<T>& L, const AgOptions& opt) {
    runner_detail::check_method(opt.method);
    runner_detail::check_model_method(L, opt.method);
    // runAnalyzerChecks' universal feature gate. `ag_feature_set` transcribes
    // SolverAG.getFeatureSet name for name and has existed since the port landed,
    // but NOTHING CALLED IT until 2026-08-19: its only two references were its own
    // definition and a comment in mam_feature_set, so every other family here was
    // gated and AG was not.
    //
    // AFTER check_model_method, following SolverCTMC: a solver's own check names
    // WHY a construct cannot be served, and the gate only names the construct. For
    // a matrix exponential that is the difference between "supports processes with
    // a Markovian (D0,D1) representation only (RCAT builds a phase dimension per
    // agent out of it)" and a bare feature name -- the set withholds ME either way,
    // so ordering decides only which message the caller reads. Reversing these two
    // lines is a silent downgrade of every AG refusal that has a specific reason.
    qn::feature_gate("SolverAG", qn::ag_feature_set(opt.method), L);

    // sn_nonmarkov_toph evaluates densities and fits moments, so it exists only
    // for a transcendental arithmetic and STATIC-ASSERTS otherwise. The guard is
    // `if constexpr` rather than a run-time test because the CLI instantiates
    // this template for the exact Rational backend as well, and an unguarded
    // call would fail to COMPILE there rather than refuse at run time. Under an
    // exact arithmetic solver_ag itself raises the refusal, which is where it
    // belongs: it is the fixed point's tolerance and the QBD's square root that
    // cannot be represented, not the conversion.
    qn::NetworkStruct<T> converted;
    const qn::NetworkStruct<T>* Lp = &L;
    if constexpr (num_traits<T>::has_transcendental) {
        if (api::sn_has_nonmarkov(L, /*preserve_det=*/false)) {
            converted = L;
            api::NonmarkovOptions no;
            no.order = opt.nonmkv_order;
            no.preserve_det = false;
            no.phfit = api::PhFit::Ph;
            api::sn_nonmarkov_toph(converted, no);
            Lp = &converted;
        }
    }

    return solver_ag(*Lp, opt);
}

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_SOLVER_AG_RUNNER_H
