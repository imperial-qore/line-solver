/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AUTO_AUTO_METHODS_H
#define LINE_SOLVERS_AUTO_AUTO_METHODS_H

/**
 * `SolverAUTO.listValidMethods`: the method names THIS MODEL can actually run.
 *
 * WHY IT IS A SEPARATE HEADER. `solver_auto.h` is deliberately light -- a
 * struct, a ranking table and the feature sets -- and it is included by the CLI,
 * the facade, the examples and `test_all_headers`. Answering this question needs
 * every family's method registry, i.e. the MVA, NC, CTMC, fluid, MAM, BA, SSA,
 * LDES, JMT and QNS runners, so putting it there would make a heavy include of a
 * cheap one for every caller that only wants the ranking.
 *
 * THE GATE IS THE FEATURE SET, asked of every candidate rather than of the first
 * feasible one. `auto_supports` already asks it for the seven ranked slots, and
 * it is the same question `chooseSolverRanked` asks before delegating: is every
 * feature the model uses declared by the solver that would run this method? A
 * family whose feature set refuses the model contributes nothing, and a method
 * whose own set refuses it is not offered. Both matter, because a per-method set
 * is where the deltas live -- MVA's `rqna` consumes MAP arrivals that the rest of
 * the envelope does not, MAM's LoadDependence holds only for the methods that
 * read `lldscaling`.
 *
 * WHAT A FLAT FEATURE SET CANNOT SAY, and is therefore asked separately, exactly
 * as the reference's `supportsModelMethod` overrides do:
 *   product form      'exact' is not exact without one (SolverMVA's OI/PAS
 *                     stations excepted); folded into `auto_supports`
 *   a chain that fits the CTMC slot is screened for state-space size
 *   a binding buffer  setCapacity/classCap, which MVA, NC, FLD and QNS refuse
 *                     through `check_binding_capacity` and no feature name
 *                     describes
 *   the model SHAPE   `mva::list_valid_methods(L)` and `ba::list_valid_methods(L)`
 *                     are themselves model-aware: the queueing-system closed
 *                     forms are offered on a two-station open model and nowhere
 *                     else, the QRF reduction bounds on a single-class closed
 *                     network of single servers, the three open-network bounds
 *                     on a fully open one. Asking the registry for the model is
 *                     what keeps that rule in one place.
 *
 * A FAMILY WHOSE EVERY METHOD IS REFUSED LOSES ITS BARE TOKEN TOO: `nc` alone
 * delegates to SolverNC, which is exactly the rejection the per-method gate just
 * returned.
 *
 * WHAT IS NOT LISTED, and why each is a token this port would refuse anyway:
 *   ln, env, lqns, uq   take a LayeredNetwork, an Environment or an inner-solver
 *                       factory rather than a Network. The reference's own
 *                       `buildFamilySolver` cannot construct them from here
 *                       either, and its `familyAcceptsModelClass` drops them.
 *   ldes without an engine  `auto_supports` folds `auto_solver_is_available` in,
 *                       and the CLI refuses `--method ldes` by name when no
 *                       engine is found beside the binary. Listing it there
 *                       would be a claim the very next call denies.
 * JMT IS LISTED although it is not a ranked candidate. It has no slot in
 * `AutoSolver` -- the reference removed it as one, LDES subsuming its feature
 * set -- but `-s auto --method jmt` reaches it as an explicit method name, and
 * `listValidMethods` answers about the method names the caller may ask for, not about
 * the ones the ranking would pick.
 *
 * THERE IS NO `list_all_methods` TWIN. In the reference the model-independent
 * list exists so that a method-NAME check can reject an unknown method name with the
 * family's own explanation rather than a flat "unsupported"; `auto_resolve_token`
 * here performs no such check -- it splits the method name and lets the family refuse
 * the submethod by name -- so a second list would have no caller.
 */

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "line/lang/qn/feature_set.h"
#include "line/util/error.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ag/solver_ag_runner.h"
#include "line/solvers/auto/solver_auto.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
// The two shape predicates the report asks for, which the analyzer header above
// does not pull in: cftp and mdd each build no generator and are restricted to a
// shape the state-space guard says nothing about.
#include "line/solvers/ctmc/solver_ctmc_cftp.h"
#include "line/solvers/ctmc/solver_ctmc_mdd_analyzer.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/solvers/wrappers/ldes/ldes_options.h"
#include "line/solvers/wrappers/qns/solver_qns.h"

namespace line {
namespace autosolver {

/**
 * The method families that solve a flat Network, in `SolverAUTO.familyNames`'
 * own order -- which is also the order an unqualified algorithm name is looked
 * up in, so it is not arbitrary.
 *
 * "ag" SITS AFTER "mam", whose RCAT names it took over, and is a family for the
 * METHOD NAME and REPORT tables only: `--method ag.inap` and the model help
 * reach SolverAG through it. It is deliberately NOT one of the ranked
 * `AutoSolver` slots and NOT in the feature-set union the automatic ranking
 * takes, so a G-network is still refused by "default" and has to be asked for
 * by name; see _kb/06-solver-catalog.md, "SolverAG owns the RCAT methods".
 */
inline std::vector<std::string> auto_network_family_names() {
    return {"mva", "nc", "ctmc", "fluid", "mam", "ag", "ba", "ssa", "ldes", "jmt", "qns"};
}

/**
 * The prefixes under which a family advertises a SECOND SPELLING of a method it
 * already declares plainly.
 *
 * THIS IS A DECLARATION, not a derivation, and belongs beside
 * `auto_family_metrics` and `auto_method_class` for the same reason: the
 * knowledge lives in the solver's own dispatch (the mva registry strips a
 * leading "amva." before selecting an algorithm) and nothing exposes it, so a
 * family that gains or loses an alias spelling must be edited into all four
 * copies in the SAME change. An omission does not fail; it puts the same
 * algorithm in the report twice.
 */
inline std::vector<std::string> auto_method_alias_prefixes(const std::string& family) {
    if (family == "mva") return {"amva."};
    return {};
}

/**
 * Is `name` a second spelling of another method this family declares?
 *
 * The remainder has to be declared too, which is what keeps the rule from eating
 * a genuine method that merely starts with the prefix: it is an alias only when
 * the thing it aliases is there beside it.
 */
inline bool auto_is_method_alias(const std::string& family, const std::string& name,
                                 const std::vector<std::string>& declared) {
    for (const std::string& prefix : auto_method_alias_prefixes(family)) {
        if (name.size() > prefix.size() && name.compare(0, prefix.size(), prefix) == 0 &&
            std::find(declared.begin(), declared.end(), name.substr(prefix.size())) !=
                declared.end())
            return true;
    }
    return false;
}

namespace detail {

/** The ranked slot a family occupies, or false when it has none. */
inline bool auto_slot_of_family(const std::string& family, AutoSolver& out) {
    if (family == "mva") { out = AutoSolver::MVA; return true; }
    if (family == "nc") { out = AutoSolver::NC; return true; }
    if (family == "mam") { out = AutoSolver::MAM; return true; }
    if (family == "fluid") { out = AutoSolver::FLUID; return true; }
    if (family == "ssa") { out = AutoSolver::SSA; return true; }
    if (family == "ctmc") { out = AutoSolver::CTMC; return true; }
    if (family == "ldes") { out = AutoSolver::LDES; return true; }
    return false;
}

/**
 * The finite-buffer gate as a predicate.
 *
 * `has_binding_capacity` is the gate's own question asked without raising, which
 * is what this needs: on the solve path a capped station the solver cannot
 * honour has to stop the run rather than be reported unconstrained, but here the
 * same verdict is a yes or a no. The families are the ones whose own runners ask
 * it -- SolverMVA, SolverNC, SolverFLD and SolverQNS -- and no feature name
 * describes a capacity, which is why it cannot ride in the feature set.
 *
 * NC EXEMPTS ITS MEM ALGORITHM, which solves the censored GE/GE/c/0;N queue and
 * therefore DOES honour the buffer. The exemption is on the LITERAL method name
 * that `nc_dispatch` branches on, not on a resolved one: nothing resolves
 * 'default' into 'mem', so exempting a default run would advertise a token that
 * dispatches somewhere the buffer is ignored.
 */
template <class T>
bool capacity_admits(const std::string& family, const qn::NetworkStruct<T>& sn,
                     const std::string& method) {
    if (family != "mva" && family != "nc" && family != "fluid" && family != "qns") return true;
    if (!qn::has_binding_capacity(sn)) return true;
    if (family == "nc" && method == "mem") {
        const nc::MemSupport ms = nc::solver_nc_mem_supports(sn);
        return ms.supported && ms.blocking;
    }
    return false;
}

/**
 * The feature set a family declares for a method, and the label its refusals
 * are reported under; false when the family solves no flat Network.
 *
 * It is the switch `auto_supports` opens on, lifted out so that the bool gate
 * and the reason-returning one read the same table. A second copy of it is how
 * a family gains a feature in one answer and not in the other.
 */
inline bool declared_feature_set(const std::string& family, const std::string& method,
                                 qn::FeatureSet& declared, std::string& label) {
    if (family == "mva") { declared = qn::mva_feature_set(method); label = "SolverMVA"; return true; }
    if (family == "nc") { declared = qn::nc_feature_set(method); label = "SolverNC"; return true; }
    if (family == "mam") { declared = qn::mam_feature_set(method); label = "SolverMAM"; return true; }
    if (family == "fluid") { declared = qn::fluid_feature_set(method); label = "SolverFLD"; return true; }
    if (family == "ssa") { declared = qn::ssa_feature_set(method); label = "SolverSSA"; return true; }
    if (family == "ctmc") { declared = qn::ctmc_feature_set(method); label = "SolverCTMC"; return true; }
    if (family == "ldes") { declared = qn::ldes_feature_set(method); label = "SolverLDES"; return true; }
    if (family == "ag") { declared = qn::ag_feature_set(method); label = "SolverAG"; return true; }
    if (family == "ba") { declared = qn::ba_feature_set(method); label = "SolverBA"; return true; }
    if (family == "jmt") { declared = qn::jmt_feature_set(method); label = "SolverJMT"; return true; }
    if (family == "qns") { declared = qn::qns_feature_set(method); label = "SolverQNS"; return true; }
    return false;
}

}  // namespace detail

/**
 * The methods a family declares ON THIS MODEL, empty when it declares none.
 *
 * The registries are asked rather than copied: `mva::list_valid_methods` and
 * `ba::list_valid_methods` take the struct and narrow themselves, and a second
 * copy of either shape rule here is how the two drift apart.
 */
template <class T>
std::vector<std::string> auto_family_methods(const std::string& family,
                                             const qn::NetworkStruct<T>& sn) {
    if (family == "mva") return mva::list_valid_methods(sn);
    if (family == "nc") return nc::list_valid_methods();
    if (family == "ctmc") return ctmc::list_valid_methods();
    if (family == "fluid") return fluid::fluid_list_valid_methods();
    if (family == "mam") return mam::list_valid_methods();
    if (family == "ag") return ag::list_valid_methods();
    if (family == "ba") return ba::list_valid_methods(sn);
    if (family == "ssa") return ssa::list_valid_methods();
    if (family == "ldes") return ldes::list_valid_methods();
    if (family == "jmt") return jmt::jmt_list_valid_methods();
    if (family == "qns") return qns::list_valid_methods();
    return std::vector<std::string>();
}

/**
 * `Solver.supportsModelMethod` for a family method name: may THIS model run THIS
 * method of THIS family?
 *
 * The seven ranked slots defer to `auto_supports`, which is the same call
 * `chooseSolverRanked` makes and already carries the product-form rule, the
 * CTMC state-space screen and the LDES engine probe. The three families with no
 * slot are gated on their own declared set, since they are reachable by explicit
 * token and by no ranking.
 */
template <class T>
std::string auto_family_refusal(const std::string& family, const qn::NetworkStruct<T>& sn,
                                const std::string& method) {
    if (!detail::capacity_admits(family, sn, method))
        return "SolverAUTO: this method ignores the finite station capacity this model sets "
               "(setCapacity / classCap), so it would report a capped station as unbounded.";

    qn::FeatureSet declared;
    std::string label;
    if (!detail::declared_feature_set(family, method, declared, label))
        return "SolverAUTO: the '" + family + "' method family does not solve a flat Network.";

    const qn::SupportResult r = qn::feature_set_supports(label, declared, qn::used_lang_features(sn));
    if (!r.ok) return r.reason;

    // NC's structural per-method rules, the ones no feature name can carry:
    // whether the method has a route on THIS model at all. `nc_method_refusal`
    // is the single copy of them, asked here and by `solver_nc_solve`, so a pair
    // this report offers is a pair the run accepts. Asked BEFORE `auto_supports`
    // because it names the offending thing, where the fallback below can only
    // say that some structural check refused the model.
    if (family == "nc") {
        const std::string ncr = nc::nc_method_refusal(sn, method);
        if (!ncr.empty()) return ncr;
    }

    // THE STRUCTURAL RULES OF THE TWO WRAPPER-ADJACENT FAMILIES, asked of the
    // families themselves rather than restated here. `jmt` has no ranked slot,
    // so the branch below never reaches it, and `auto_supports` asks the MAM
    // slot for its feature set only -- which is how findSolver came to offer
    // jmt.jmva.<alg> on a multi-server model, jmt.replication with no horizon
    // and mam.retrial on a model that declares no orbit, each of which then
    // raised when it was run. The two calls are the solvers' OWN predicates,
    // the same ones their runners raise, taken with the default options a probe
    // solver carries.
    if (family == "jmt") {
        const std::string why = jmt::jmt_method_refusal(sn, method, jmt::JmtOptions());
        if (!why.empty()) return why;
    }
    if (family == "mam") {
        const std::string why = mam::mam_model_method_refusal(sn, method);
        if (!why.empty()) return why;
    }
    // AG's own structural gate, on the same terms: RCAT needs a Markovian
    // (D0,D1) service law, state-independent routing, one server per station and
    // no binding buffer, and none of the four is a feature name. It has no
    // ranked slot, so nothing below reaches it.
    if (family == "ag") {
        const std::string why = ag::runner_detail::method_refusal(sn, method);
        if (!why.empty()) return why;
    }
    // QNS on the same terms: immediate feedback is outside BOTH conversions,
    // and the multiserver approximations `qnsolver -m` does not offer are a
    // rule about the method rather than about the model.
    if (family == "qns") {
        const std::string why = qns::method_refusal(sn, method);
        if (!why.empty()) return why;
    }

    // THE REMAINING THREE FAMILIES, on the same terms. `ba` has no ranked slot
    // either, so nothing below reaches it; `ctmc` and `fluid` do have one, but
    // `auto_supports` asks each for a rule about the SOLVER (state-space size,
    // product form) and never about the METHOD, which is where a class count
    // and a horizon live. Each call is the family's own predicate, the same one
    // its analyzer raises, so a pair this report offers is a pair the run
    // accepts and the two cannot drift.
    if (family == "ba") {
        const std::string why = ba::method_refusal(sn, method);
        if (!why.empty()) return "SolverBA: " + why;
    }
    if (family == "ctmc") {
        // Only the two shape-restricted methods have one; the rest are gated by
        // the state-space size below, which is a question about the solver.
        const std::string why =
            (method == "cftp" || method == "cftp.approx")
                ? ctmc::solver_ctmc_cftp_supports(sn)
                : (method == "mdd" ? ctmc::solver_ctmc_mdd_supports(sn) : std::string());
        if (!why.empty()) return why;
    }
    // THE FORK-JOIN MODEL CLASS, which is a question about the MODEL and not
    // about the method, so both exact families have to clear it whichever name
    // was asked for: each tag-augments through `fj_tag`, whose first line is
    // `sn_fj_validate`. The feature set cannot state it -- Fork and Join are
    // declared, and the rules are about how they are WIRED -- and the pairing in
    // particular is a DECLARATION on the Join rather than something the routing
    // implies, so a Join that names no fork leaves an unmatched Fork behind.
    if (family == "ctmc" || family == "ssa") {
        const std::string why = qn::sn_fj_supports(sn);
        if (!why.empty()) return why;
    }
    // `nrm` IS THE ONE SSA METHOD WITH A MODEL CLASS OF ITS OWN, and the test
    // already existed: `ssa_nrm_eligible` decides whether the dispatch may
    // PREFER the NRM, while nothing decided whether the name could be OFFERED.
    // So `ssa.nrm` was reported runnable on every model and an explicit request
    // then raised from the analyzer. `ssa_nrm_supports` is the same six checks
    // read as a sentence.
    if (family == "ssa" && method == "nrm") {
        const std::string why = ssa::detail::ssa_nrm_supports(sn);
        if (!why.empty()) return why;
    }
    if (family == "fluid") {
        // The horizon is an OPTION, not a model feature, so it cannot be a
        // feature-set delta. A probe carries the defaults, under which the
        // time-varying methods have no finite end and the refusal says so.
        const std::string why = fluid::fluid_qsys_horizon_supports(method, fluid::FluidOptions());
        if (!why.empty()) return why;
        // Fork AND open is a CONJUNCTION of two declared names, which no feature
        // set can state; `dae` is the one method the MMT fixed point has no
        // route for on an open model.
        const std::string fj = fluid::detail::fluid_forkjoin_supports(sn, method);
        if (!fj.empty()) return fj;
    }

    // WHAT auto_supports STILL REFUSES once the feature set has passed is one of
    // its two residual rules, and each is named rather than reported as a bare
    // no. Asking auto_supports rather than repeating the rules keeps it the one
    // authority; the branches below only turn its verdict into a sentence.
    AutoSolver slot = AutoSolver::MVA;
    if (detail::auto_slot_of_family(family, slot)) {
        if (!auto_solver_is_available(slot))
            return std::string(auto_solver_name(slot)) +
                   ": no engine for this solver is available in this build.";
        if (!auto_supports(slot, sn, method)) {
            if (method == "exact" && (slot == AutoSolver::MVA || slot == AutoSolver::NC))
                return std::string(auto_solver_name(slot)) +
                       ": 'exact' needs a product-form model, and this one has no product-form "
                       "solution.";
            if (slot == AutoSolver::CTMC)
                return "SolverCTMC: the state space of this model is too large to enumerate.";
            return label + ": this solver refuses the model through its own structural check.";
        }
    }
    return "";
}

template <class T>
bool auto_family_supports(const std::string& family, const qn::NetworkStruct<T>& sn,
                          const std::string& method) {
    return auto_family_refusal(family, sn, method).empty();
}

// ---------------------------------------------------------------------------
// findSolver: which solvers and methods can analyze this model
// Port of @SolverAUTO/findSolver.m and its native python and JAR twins.
// ---------------------------------------------------------------------------

/**
 * One row of `auto_find_solver`: a (family, method) pair this model can be
 * asked for, whether it runs, what kind of answer it returns and which
 * measures it can report.
 *
 * The fields are the columns of the MATLAB table, of the native python frame
 * and of the JAR's SolverCandidate, under the same names, so a report can be
 * compared across the four codebases row for row.
 */
struct SolverCandidate {
    std::string solver;        ///< the method family: "mva", "ctmc", "ldes", ...
    std::string method;        ///< the method name to pass, "mva.exact"
    bool runnable = false;     ///< the model passes this method's own support gate
    std::string method_class;  ///< "exact", "approx", "bound" or "simulation"
    std::string metrics;       ///< the measure groups the family answers, comma-joined
    std::string reason;        ///< why a refused pair was refused, "" when runnable
};

/**
 * The measure groups `auto_find_solver` reports on, in report order.
 *
 * A group is a family of accessors that stand or fall together: a solver that
 * returns getCdfRespT returns getCdfPassT and getPerctRespT as well, because
 * all three read the same passage time, so listing the three separately would
 * say nothing extra.
 */
inline std::vector<std::string> auto_metric_groups() {
    return {"avg", "tran", "cdf", "prob", "tranprob", "sample",
            "cache", "loss", "orbit", "moment", "sens"};
}

namespace detail {

inline bool name_in(const std::string& n, std::initializer_list<const char*> names) {
    for (const char* c : names)
        if (n == c) return true;
    return false;
}

inline bool starts_with(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

}  // namespace detail

/**
 * The measure group an accessor belongs to, "" when the name belongs to none.
 *
 * A group name maps to itself, so `auto_find_solver(sn, "cdf")` and
 * `auto_find_solver(sn, "getCdfRespT")` ask the same question.
 *
 * THIS IS NOT `auto_choose_solver`'S TABLE, although both are keyed by accessor
 * name. That one maps an accessor to a RANKING, i.e. which candidate should be
 * preferred; this one maps it to a CAPABILITY question, i.e. which candidates
 * can answer it at all. The two differ wherever a family can serve a measure
 * but is never the one AUTO would pick for it.
 */
inline std::string auto_metric_group_of(const std::string& name) {
    if (name.empty()) return "";
    for (const std::string& g : auto_metric_groups())
        if (g == name) return g;
    if (name == "any" || name == "all") return "";
    if (detail::name_in(name, {"getTranAvg", "getTranAvgVar", "tranAvg"})) return "tran";
    if (detail::name_in(name, {"getCdfRespT", "getCdfPassT", "getPerctRespT", "getTranCdfPassT",
                               "getTranCdfRespT", "getCdfSysRespT"}))
        return "cdf";
    if (detail::name_in(name, {"getTranProb", "getTranProbSys", "getTranProbAggr",
                               "getTranProbSysAggr"}))
        return "tranprob";
    if (detail::name_in(name, {"getProb", "getProbAggr", "getProbSys", "getProbSysAggr",
                               "getProbMarg", "getProbNormConstAggr"}))
        return "prob";
    if (detail::name_in(name, {"sample", "sampleSys", "sampleAggr", "sampleSysAggr"}))
        return "sample";
    if (detail::name_in(name, {"getAvgCacheTable", "getAvgCacheT", "getAvgItemTable",
                               "getAvgItemT", "cacheAvgT", "itemAvgT", "aCaT", "aIT"}))
        return "cache";
    if (detail::name_in(name, {"getAvgLossTable", "getAvgLossT", "getAvgRegionLossTable",
                               "getAvgRegionLossT", "lossAvgT", "regionLossAvgT", "aLT", "aRLT"}))
        return "loss";
    if (detail::name_in(name, {"getAvgOrbitTable", "getAvgOrbitT", "getAvgOrbit", "orbitAvgT",
                               "aOT"}))
        return "orbit";
    if (detail::name_in(name, {"getMomentTable", "getMomentChainTable", "getMomentStationTable",
                               "getMomentT", "getMomentChainT", "getMomentStationT", "momentT",
                               "momentChainT", "momentStationT", "mT", "mCT", "mST"}))
        return "moment";
    if (detail::name_in(name, {"getSensitivityTable", "getSensitivityT", "sensitivityT", "sT",
                               "getSensitivity", "getSensitivityRanking"}))
        return "sens";
    // Everything else in the accessor surface is a mean measure: getAvg, its
    // chain, node and system forms, their handles and their short aliases.
    if (detail::starts_with(name, "getAvg") || detail::starts_with(name, "avg") ||
        detail::name_in(name, {"getAvgSysRespT", "getAvgSysTput", "aT", "aNT", "aCT", "aST",
                               "aNCT"}))
        return "avg";
    return "";
}

/**
 * The measure groups a method family can answer.
 *
 * Every family answers "avg", which is what a solver is for; the rest is the
 * capability declaration this header owns.
 *
 * SOURCES, so that a claim here can be checked rather than trusted: "tran" is
 * the reference's supportsTransientAnalysis, which FLD, CTMC, LDES and JMT
 * override to true and no one else does. "cdf", "prob", "tranprob" and "sample"
 * are the families that carry an implementation of the corresponding accessor
 * rather than inheriting the base refusal. The remaining five groups are
 * computed from a solver's own results, so no per-solver entry point marks
 * them: their lists are the reference chooseSolverHeur's rankings for the same
 * accessors, which is where AUTO already records who can serve them.
 *
 * A family that gains or loses a measure must be edited here in the same
 * change, the way a solver that gains a feature is edited into its feature set:
 * an omission here does not fail, it silently hides the family from a caller
 * asking for that measure.
 */
inline std::vector<std::string> auto_family_metrics(const std::string& family) {
    if (family == "mva") return {"avg", "prob", "cache", "orbit", "moment", "sens"};
    if (family == "nc") return {"avg", "cdf", "prob", "cache", "moment", "sens"};
    if (family == "ctmc")
        return {"avg", "tran", "cdf", "prob", "tranprob", "sample", "cache", "loss", "orbit",
                "moment"};
    if (family == "fluid") return {"avg", "tran", "cdf", "prob", "cache", "sens"};
    if (family == "mam") return {"avg", "cdf"};
    // The RCAT fixed point reports means only; the passage-time law it answers
    // is the base exponential fit, not its own.
    if (family == "ag") return {"avg"};
    // A bound brackets the mean measures and nothing else.
    if (family == "ba") return {"avg"};
    if (family == "ssa") return {"avg", "cdf", "prob", "sample", "loss"};
    if (family == "ldes")
        return {"avg", "tran", "cdf", "prob", "sample", "cache", "loss", "orbit"};
    if (family == "jmt") return {"avg", "tran", "cdf", "prob", "tranprob", "sample"};
    return {"avg"};
}

/**
 * `Solver.isStochasticMethod` for a family method name: does this method return
 * seed-dependent estimates?
 *
 * SSA and LDES are simulators outright; JMT is one except through its
 * analytical JMVA engine, whose own sampling variants are stochastic again; NC
 * has the Monte Carlo, importance-sampling and MCMC routes, and answers for
 * itself through `nc::is_stochastic_method` rather than through a second copy
 * of that token list here.
 */
inline bool auto_is_stochastic_method(const std::string& family, const std::string& method) {
    if (family == "ssa" || family == "ldes") return true;
    if (family == "nc") return nc::is_stochastic_method(method);
    if (family == "jmt") {
        std::string tok;
        std::vector<std::string> toks;
        for (char ch : method) {
            if (ch == '.' || ch == '/') {
                toks.push_back(tok);
                tok.clear();
            } else {
                tok += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
            }
        }
        toks.push_back(tok);
        bool jmva = false;
        for (const std::string& t : toks)
            if (t == "jmva") jmva = true;
        if (!jmva) return true;
        for (const std::string& t : toks)
            if (t == "ls" || t == "mci" || t == "imci" || t == "sampling") return true;
        return false;
    }
    return false;
}

namespace detail {
inline std::string exact_if(bool cond) { return cond ? "exact" : "approx"; }
}  // namespace detail

/**
 * What KIND of answer a method returns: "exact", "approx", "bound" or
 * "simulation".
 *
 * "simulation" is not decided here: `is_stochastic` is the solver's own verdict,
 * which already tokenizes qualified and runtime-resolved names.
 *
 * "exact" IS CLAIMED ONLY WHERE IT IS TRUE OF THIS MODEL, never of the
 * algorithm in the abstract. Exactness of a normalizing constant or of mean
 * value analysis is a property of the product-form model it is computed on, and
 * of the QBD shape for the matrix analytic methods, so both conditions are
 * passed in and a method that needs one reports "approx" without it. The bias
 * is deliberate: an under-claimed "approx" costs a user a better method they
 * could have had, an over-claimed "exact" costs them a wrong number they
 * trusted.
 *
 * A CACHE IS THE THIRD CONDITION, and it was the over-claim the bias above
 * exists to prevent. `has_product_form` answers about the QUEUEING network and
 * knows nothing of a cache: the hit/miss split is a class switch whose
 * probabilities are not routing data but the output of a cache model, so a
 * network holding one reads as product form and "mva.exact" was labelled exact
 * on it. Measured on the tut06 shape with an LRU cache: exact MVA returns QLen
 * 0.2516 at the hit station where the CTMC returns 0.3022 and simulation 0.3023,
 * a 17% error under a label that says there is none. The analytic families are
 * conditioned on it; SolverCTMC is NOT, because its state space carries the
 * cache contents and it is exact there, which is what the two numbers show.
 */
inline std::string auto_method_class(const std::string& family, const std::string& method,
                                     bool is_stochastic, bool is_product_form,
                                     bool is_qbd_shape, bool has_cache = false) {
    // Bounds are what SolverBA is for; every one of its methods returns a
    // bracket rather than an estimate.
    if (family == "ba") return "bound";
    if (is_stochastic) return "simulation";
    if (family == "ctmc") {
        // The generator is solved as written, so every state-space route is
        // exact. "cftp.approx" says in its own name that it is not, and "mdd"
        // is exact on a product-form model and an approximation otherwise.
        if (method == "cftp.approx") return "approx";
        if (method == "mdd") return detail::exact_if(is_product_form);
        return "exact";
    }
    if (family == "nc") {
        // The normalizing-constant routes that evaluate G exactly rather than
        // expanding or estimating it. The asymptotic expansions (pana, le,
        // kt, bk, gm, ...) and the non-product-form "morrison" are
        // approximations by construction and are left out.
        if (detail::name_in(method, {"exact", "divdiff", "ca", "comom", "comomld", "rec", "ms",
                                     "cub", "rgf"}))
            return detail::exact_if(is_product_form && !has_cache);
        return "approx";
    }
    if (family == "mva") {
        // Exact MVA; every "amva.*" arm is an approximation, and so are the
        // open-network QNA transforms.
        if (detail::name_in(method, {"exact", "mva"}))
            return detail::exact_if(is_product_form && !has_cache);
        return "approx";
    }
    if (family == "jmt") {
        // JMVA's exact algorithms. "jsim" and "replication" are simulation and
        // never reach here.
        if (detail::name_in(method, {"jmva.mva", "jmva.recal", "jmva.comom", "jmva.treeconv"}))
            return detail::exact_if(is_product_form && !has_cache);
        return "approx";
    }
    if (family == "mam") {
        // The QBD is solved exactly on the shape it is stated for, one queueing
        // station fed by a Source. Everything named "dec.*" is a decomposition
        // of a larger network into such queues, hence an approximation of it.
        if (detail::name_in(method, {"default", "mna", "ldqbd", "bgchain", "retrial"}))
            return detail::exact_if(is_qbd_shape);
        return "approx";
    }
    if (family == "ag") {
        // Every RCAT arm estimates the reversed rate of each synchronising
        // action and iterates to a fixed point, an approximation by
        // construction; SolverAG's "exact" is a vestigial alias that warns and
        // runs "inap", so nothing here is claimed exact.
        return "approx";
    }
    return "approx";
}

/**
 * `SolverAUTO.findSolver`: which solvers and solver methods can analyze this
 * model, and for the ones that cannot, why not.
 *
 * THE GATE IS NOT A SECOND ONE. It is `auto_family_supports`, the gate
 * `auto_choose_solver` applies before delegating, asked of every candidate
 * instead of of the first feasible one -- which is exactly what
 * `auto_list_valid_methods` already did; that function is now the method column
 * of the runnable rows, so the two cannot disagree. What is new is that the
 * REASON is kept rather than discarded, and that the answer carries the two
 * facts a caller needs in order to choose among the survivors: whether the
 * method is exact on this model, and which measures it can report.
 *
 * `metric` narrows the report to the pairs that answer one measure, named
 * either by its group ("cdf") or by the accessor that returns it
 * ("getCdfRespT"); "" or "any" keeps every pair. `show_all` keeps the refused
 * pairs too; by default only the runnable ones are listed, since a caller
 * asking what it can run has no use for the rows that say it cannot.
 */
template <class T>
std::vector<SolverCandidate> auto_find_solver(const qn::NetworkStruct<T>& sn,
                                              const std::string& metric = std::string(),
                                              bool show_all = false) {
    const std::string group = auto_metric_group_of(metric);
    if (group.empty() && !metric.empty() && metric != "any" && metric != "all") {
        std::string groups;
        for (const std::string& g : auto_metric_groups()) {
            if (!groups.empty()) groups += ", ";
            groups += g;
        }
        // InputError and not a bare runtime_error: this is a CALLER mistake, and
        // the CLI reports the two differently -- an input error prints
        // "line-cli: <message>" and exits 2, anything else prints "unexpected
        // failure" and exits 3, which is what a typo'd measure was getting.
        throw InputError("'" + metric + "' names no measure. Pass a group (" + groups +
                         ") or the accessor that returns it, e.g. 'getCdfRespT'.");
    }

    // The two model properties an exactness claim can rest on, evaluated once:
    // a method whose exactness needs one of them reports "approx" without it.
    const bool is_product_form = sn.has_product_form();
    std::size_t nsources = 0;
    for (const qn::Station<T>& st : sn.stations)
        if (st.nodetype == lang::NodeType::Source) ++nsources;
    const std::vector<double> njobs = sn.njobs();
    bool all_open = !njobs.empty();
    for (double n : njobs)
        if (!std::isinf(n)) all_open = false;
    const bool is_qbd_shape = all_open && (sn.nstations - nsources) == 1;
    bool has_cache = false;
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == qn::NodeType::Cache) has_cache = true;

    std::vector<SolverCandidate> rows;
    for (const std::string& family : auto_network_family_names()) {
        const std::vector<std::string> groups = auto_family_metrics(family);
        if (!group.empty() &&
            std::find(groups.begin(), groups.end(), group) == groups.end())
            continue;
        std::string metric_list;
        for (const std::string& g : groups) {
            if (!metric_list.empty()) metric_list += ",";
            metric_list += g;
        }
        const std::vector<std::string> declared = auto_family_methods(family, sn);
        for (const std::string& name : declared) {
            if (detail::starts_with(name, family + ".")) {
                // A spelling already qualified with its own family. The fluid
                // registry declares both "dae" and "fluid.dae" so that its own
                // gate takes either, and prefixing the family again yields
                // "fluid.fluid.dae": a token that does resolve, but that names
                // the same method twice and would double every fluid row.
                continue;
            }
            if (auto_is_method_alias(family, name, declared)) {
                // The same duplication under a different prefix. The mva
                // registry advertises every AMVA name twice, plain and
                // "amva."-prefixed, and the dispatch strips the prefix, so the
                // two spellings are one algorithm; that alone was 20 of the 49
                // mva rows of a report. The plain spelling is the one kept.
                continue;
            }
            const std::string reason = auto_family_refusal(family, sn, name);
            const bool ok = reason.empty();
            if (!ok && !show_all) continue;
            SolverCandidate row;
            row.solver = family;
            row.method = family + "." + name;
            row.runnable = ok;
            row.method_class = auto_method_class(
                family, name, auto_is_stochastic_method(family, name), is_product_form,
                is_qbd_shape, has_cache);
            row.metrics = metric_list;
            row.reason = ok ? std::string() : reason;
            rows.push_back(row);
        }
    }
    return rows;
}

/** Alias of `auto_find_solver`: the same table, asked for by method. */
template <class T>
std::vector<SolverCandidate> auto_find_method(const qn::NetworkStruct<T>& sn,
                                              const std::string& metric = std::string(),
                                              bool show_all = false) {
    return auto_find_solver(sn, metric, show_all);
}

/** Alias of `auto_find_solver`: what can this model be solved with? */
template <class T>
std::vector<SolverCandidate> auto_help(const qn::NetworkStruct<T>& sn,
                                       const std::string& metric = std::string(),
                                       bool show_all = false) {
    return auto_find_solver(sn, metric, show_all);
}

/**
 * The rows as an aligned text table, the form the CLI and a console caller
 * want. The reason column is last and unpadded, since it is the only one whose
 * width is unbounded.
 */
inline std::string auto_find_solver_table(const std::vector<SolverCandidate>& rows) {
    if (rows.empty()) return "No solver method can analyze this model.\n";
    std::vector<std::array<std::string, 6>> cells;
    cells.push_back({"Solver", "Method", "Runnable", "Class", "Metrics", "Reason"});
    for (const SolverCandidate& r : rows)
        cells.push_back({r.solver, r.method, r.runnable ? "true" : "false", r.method_class,
                         r.metrics, r.reason});
    std::array<std::size_t, 6> width{};
    for (const auto& row : cells)
        for (std::size_t c = 0; c < 5; ++c) width[c] = std::max(width[c], row[c].size());
    std::string out;
    for (const auto& row : cells) {
        std::string line;
        for (std::size_t c = 0; c < 5; ++c) line += row[c] + std::string(width[c] - row[c].size() + 2, ' ');
        line += row[5];
        while (!line.empty() && line.back() == ' ') line.pop_back();
        out += line + "\n";
    }
    return out;
}

/**
 * `SolverAUTO.listValidMethods`: every method name this model can be asked for.
 *
 * The selection intents come first and unconditionally: they name a RANKING and
 * not an algorithm, and finding a family that supports the model is the
 * ranking's own job. `bound` is the one that is not quite a ranking -- it names
 * SolverBA with method 'auto' -- and it is listed unconditionally all the same,
 * because MATLAB and the JAR list it that way and a caller reading 'bound' as
 * "give me bounds" should get SolverBA's own refusal rather than a missing
 * method name. Then each family in `auto_network_family_names` order, bare method name
 * first and its qualified methods after it.
 */
template <class T>
std::vector<std::string> auto_list_valid_methods(const qn::NetworkStruct<T>& sn) {
    // IT IS THE RUNNABLE ROWS OF `auto_find_solver`, projected onto their
    // method name. The narrowing used to be written out a second time here,
    // and a second copy of one gate is how two answers to one question start to
    // differ; `auto_find_solver` owns it now, and this adds only the method names
    // that name no single method: the selection intents and each family's bare
    // token.
    //
    // A family that declares methods and keeps none loses its bare method name too:
    // it delegates to the solver the per-method gate just refused. That falls
    // out of the projection, since such a family contributes no row to name.
    std::vector<std::string> out{"accurate", "auto", "bound", "default",
                                 "exact",    "fast", "heur",  "sim"};
    for (const SolverCandidate& row : auto_find_solver(sn)) {
        out.push_back(row.method);
        out.push_back(row.solver);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

}  // namespace autosolver
}  // namespace line

#endif  // LINE_SOLVERS_AUTO_AUTO_METHODS_H
