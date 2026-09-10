/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The universal language-feature gate: the port of SolverFeatureSet.m,
 * getUsedLangFeatures.m and NetworkSolver.runAnalyzerChecks.
 *
 * WHAT IS ORACLED. Three things, and the third is the point of the exercise:
 *   1. the registry is byte-identical to SolverFeatureSet.fields, in order, so
 *      that a name cannot go missing (SUPPORTS iterates the registry, so an
 *      unregistered capability is invisible and passes as if absent);
 *   2. the used-feature derivation, one construct at a time, against what
 *      getUsedLangFeatures.m emits for the same construct;
 *   3. the REASON strings, because a refusal that does not name the offending
 *      construct leaves the user with no way to find it.
 *
 * The structs here are built field by field rather than through NetworkBuilder:
 * the derivation reads the struct, and a hand-built struct isolates it from the
 * refresh, which is what the parity regressions cover.
 */
#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/lang/qn/network_builder.h"

using line::Matrix;
using line::UnsupportedError;
using line::lang::RemovalPolicy;
using line::lang::SignalType;
using line::qn::CacheParam;
using line::qn::DropStrategy;
using line::qn::Feature;
using line::qn::feature_count;
using line::qn::feature_gate;
using line::qn::feature_name;
using line::qn::feature_of_process;
using line::qn::feature_of_replacement;
using line::qn::feature_of_routing;
using line::qn::feature_of_sched;
using line::qn::feature_of_server_section;
using line::qn::feature_phrase;
using line::qn::FeatureSet;
using line::qn::feature_generalization;
using line::qn::feature_set_supports;
using line::qn::SupportResult;
using line::qn::JobClass;
using line::qn::JobClassType;
using line::qn::NetworkStruct;
using line::qn::NodeType;
using line::qn::ProcessType;
using line::qn::ReplacementStrategy;
using line::qn::RetrialParam;
using line::qn::RoutingStrategy;
using line::qn::SchedStrategy;
using line::qn::Station;
using line::qn::TransitionParam;
using line::qn::used_lang_features;

using Dist = line::lang::Distrib<double>;
using SN = NetworkStruct<double>;

namespace {

/** A closed class, the minimum a station needs to carry a service process. */
JobClass closed_class(const std::string& nm, double njobs) {
    JobClass c;
    c.name = nm;
    c.type = JobClassType::CLOSED;
    c.population = njobs;
    return c;
}

JobClass open_class(const std::string& nm) {
    JobClass c;
    c.name = nm;
    c.type = JobClassType::OPEN;
    c.population = std::numeric_limits<double>::infinity();
    return c;
}

/** A station of a given kind and discipline, with routing set for every class. */
std::size_t add_st(SN& sn, const std::string& nm, NodeType ty, SchedStrategy sched,
                   RoutingStrategy rs = RoutingStrategy::PROB) {
    Station<double> st;
    st.name = nm;
    st.nodetype = ty;
    st.sched = sched;
    const std::size_t ist = sn.add_station(st);
    sn.nodes[sn.station_to_node[ist - 1] - 1].routing.assign(sn.classes.size(), rs);
    return ist;
}

/** One closed class at one FCFS queue with exponential service: the baseline. */
SN baseline() {
    SN sn;
    sn.add_class(closed_class("c1", 2.0));
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.set_service(q, 1, Dist::exp_mean(1.0));
    sn.nclasses = sn.classes.size();
    sn.nstations = sn.stations.size();
    return sn;
}

bool holds(const FeatureSet& u, std::initializer_list<Feature> fs) {
    for (Feature f : fs)
        if (!u.has(f)) return false;
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// The registry
// ---------------------------------------------------------------------------

TEST_CASE("registry size and names match SolverFeatureSet.fields") {
    // 172 = numel(SolverFeatureSet.fields), which is 169, PLUS the three names
    // this port emitted from `used_lang_features` and MATLAB had no counterpart
    // for: HeteroServers, DepartureDiscipline and JoinPartial. JoinPartial is
    // no longer one of them -- it joined SolverFeatureSet.fields, FeatureSet.java
    // and the python FIELDS when the quorum join was gated in all four
    // codebases, so the registries now differ by two names. The order is the
    // reference's, so Region sits between Buffer and Linkage and
    // SchedStrategy_REF is last. The ORDER is not a cross-codebase invariant --
    // MATLAB and the JAR carry the same 169 names in different positions -- so
    // only the names and the two anchors below are asserted.
    // 2026-08-09: 170 -> 171 with RoutingStrategy_SDR, the Krzesinski (1987)
    // state-dependent routing, declared in SolverFeatureSet.m:105.
    // 2026-08-09: 171 -> 172 on the merge, ServerParallelism having joined the
    // registry alongside HeteroServers on master in the same window. BOTH
    // BRANCHES WROTE 171 -- each counting only its own addition -- so the line
    // merged without a conflict and the number had to be recomputed by hand.
    // 2026-08-14: 172 -> 175 with CacheTask, ItemEntry and
    // ActivityPrecedence_POST_CACHE, the three layered cache-queueing names the
    // JAR SolverLDES.getLNFeatureSet DECLARES and no registry carried, so that
    // whole feature set threw on construction. Appended after SchedStrategy_REF
    // so every earlier index is unchanged.
    // 2026-08-14: 175 -> 178 with ForkFanoutVector, ForkFanoutRandom and
    // ForkBranchProbability, the variable-forking-level names added to all four
    // registries at once. Appended at the tail for the same reason, and MARKED
    // in used_lang_features -- a name declared but never marked is a name no
    // solver can refuse.
    // 2026-08-22: the count is UNCHANGED at 178 but the registries no longer
    // differ: HeteroServers and DepartureDiscipline, the last two names this
    // port carried alone, joined SolverFeatureSet.fields, FeatureSet.java and
    // the python FIELDS (appended at THEIR tails, so their indices differ from
    // the positions used here -- the order was never a cross-codebase
    // invariant). All four registries now hold the same 178 names.
    // 2026-09-05: 178 -> 180 with MultiServer and FiniteCapacity, appended at
    // the TAIL of LINE_QN_FEATURE_LIST so every earlier index is unchanged.
    CHECK(feature_count() == 180);
    CHECK(std::string(feature_name(Feature::CacheTask)) == "CacheTask");
    CHECK(std::string(feature_name(Feature::ActivityPrecedence_POST_CACHE)) ==
          "ActivityPrecedence_POST_CACHE");
    CHECK(std::string(feature_name(Feature::ForkFanoutVector)) == "ForkFanoutVector");
    CHECK(std::string(feature_name(Feature::ForkBranchProbability)) == "ForkBranchProbability");
    CHECK(std::string(feature_name(Feature::ClassSwitch)) == "ClassSwitch");
    CHECK(std::string(feature_name(Feature::Region)) == "Region");
    CHECK(std::string(feature_name(Feature::JointDependence)) == "JointDependence");
    CHECK(std::string(feature_name(Feature::SchedStrategy_REF)) == "SchedStrategy_REF");
    CHECK(std::string(feature_name(Feature::SignalType_CATASTROPHE)) == "SignalType_CATASTROPHE");
    CHECK(static_cast<int>(Feature::ClassSwitch) == 0);
    CHECK(static_cast<int>(Feature::SchedStrategy_REF) == 171);
}

TEST_CASE("every registry name is non-empty and unique") {
    std::vector<std::string> names;
    for (std::size_t i = 0; i < feature_count(); ++i)
        names.push_back(feature_name(static_cast<Feature>(i)));
    for (const std::string& n : names) {
        INFO("name = ", n);
        CHECK_FALSE(n.empty());
        CHECK(n != "Unknown");
    }
    std::vector<std::string> sorted = names;
    std::sort(sorted.begin(), sorted.end());
    CHECK(std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end());
}

// ---------------------------------------------------------------------------
// FeatureSet
// ---------------------------------------------------------------------------

TEST_CASE("FeatureSet set, has, unset and registry-ordered listing") {
    FeatureSet s;
    CHECK(s.empty());
    CHECK_FALSE(s.has(Feature::Region));
    s.set({Feature::SchedStrategy_REF, Feature::Region, Feature::ClassSwitch});
    CHECK_FALSE(s.empty());
    CHECK(holds(s, {Feature::Region, Feature::ClassSwitch, Feature::SchedStrategy_REF}));
    // list() is in registry order, not insertion order.
    const std::vector<Feature> l = s.list();
    REQUIRE(l.size() == 3);
    CHECK(l[0] == Feature::ClassSwitch);
    CHECK(l[1] == Feature::Region);
    CHECK(l[2] == Feature::SchedStrategy_REF);
    s.unset(Feature::Region);
    CHECK_FALSE(s.has(Feature::Region));
    CHECK(s.has(Feature::ClassSwitch));
}

TEST_CASE("Feature::COUNT is inert in set, unset and has") {
    FeatureSet s;
    s.set(Feature::COUNT);
    CHECK(s.empty());
    CHECK_FALSE(s.has(Feature::COUNT));
    s.unset(Feature::COUNT);
    CHECK(s.empty());
}

// ---------------------------------------------------------------------------
// supports and its reason strings
// ---------------------------------------------------------------------------

TEST_CASE("supports accepts a declared superset") {
    FeatureSet declared, used;
    declared.set({Feature::Queue, Feature::Exp, Feature::SchedStrategy_FCFS, Feature::Region});
    used.set({Feature::Queue, Feature::Exp});
    const auto r = feature_set_supports("SolverMVA", declared, used);
    CHECK(r.ok);
    CHECK(r.missing.empty());
    CHECK(r.reason.empty());
}

TEST_CASE("a single missing feature is named in the reason") {
    FeatureSet declared, used;
    declared.set({Feature::Queue, Feature::Delay, Feature::Exp});
    used.set({Feature::Queue, Feature::Exp, Feature::Region});
    const auto r = feature_set_supports("SolverFluid", declared, used);
    CHECK_FALSE(r.ok);
    REQUIRE(r.missing.size() == 1);
    CHECK(r.missing[0] == Feature::Region);
    CHECK(r.reason ==
          "SolverFluid: this model uses a Finite Capacity Region, which this solver does not "
          "support");
}

TEST_CASE("several missing features are listed in registry order") {
    FeatureSet declared, used;
    declared.set(Feature::Queue);
    // Deliberately set in reverse registry order to prove the report is sorted.
    used.set({Feature::Queue, Feature::Retrial, Feature::Region, Feature::Fork});
    const auto r = feature_set_supports("SolverNC", declared, used);
    CHECK_FALSE(r.ok);
    REQUIRE(r.missing.size() == 3);
    CHECK(r.missing[0] == Feature::Fork);
    CHECK(r.missing[1] == Feature::Region);
    CHECK(r.missing[2] == Feature::Retrial);
    CHECK(r.reason ==
          "SolverNC: this model uses features which this solver does not support: a Fork node, "
          "a Finite Capacity Region, a retrial orbit");
}

TEST_CASE("the comparison is one-sided: only used-and-undeclared counts") {
    FeatureSet declared, used;
    declared.set({Feature::Region, Feature::Fork, Feature::Join});
    used.set(Feature::Queue);
    const auto r = feature_set_supports("SolverCTMC", declared, used);
    // The three declared-but-unused features are silent; only Queue is reported.
    CHECK_FALSE(r.ok);
    REQUIRE(r.missing.size() == 1);
    CHECK(r.missing[0] == Feature::Queue);
}

TEST_CASE("phrases are derived for the structured families") {
    CHECK(feature_phrase(Feature::SchedStrategy_LCFSPR) ==
          "the LCFSPR scheduling discipline");
    CHECK(feature_phrase(Feature::RoutingStrategy_JSQ) == "the JSQ routing strategy");
    CHECK(feature_phrase(Feature::ReplacementStrategy_LRU) ==
          "the LRU cache replacement policy");
    CHECK(feature_phrase(Feature::SignalType_NEGATIVE) == "NEGATIVE signals");
    CHECK(feature_phrase(Feature::ActivityPrecedence_PRE_AND) ==
          "the PRE_AND activity precedence");
    // Unnamed leaf names fall back to the registry name.
    CHECK(feature_phrase(Feature::Erlang) == "the Erlang feature");
    // The hand-written phrases are the ones a user must be able to act on.
    CHECK(feature_phrase(Feature::Region) == "a Finite Capacity Region");
    CHECK(feature_phrase(Feature::Retrial) == "a retrial orbit");
    CHECK(feature_phrase(Feature::ClassDependence) == "class-dependent service rates");
    CHECK(feature_phrase(Feature::CacheRetrieval) == "a delayed-hit cache retrieval system");
}

// ---------------------------------------------------------------------------
// Enum to feature
// ---------------------------------------------------------------------------

TEST_CASE("Immediate and Disabled are not user-facing distributions") {
    CHECK(feature_of_process(ProcessType::IMMEDIATE) == Feature::COUNT);
    CHECK(feature_of_process(ProcessType::DISABLED) == Feature::COUNT);
    CHECK(feature_of_process(ProcessType::NONE) == Feature::COUNT);
    CHECK(feature_of_process(ProcessType::EXP) == Feature::Exp);
    CHECK(feature_of_process(ProcessType::COX2) == Feature::Cox2);
    CHECK(feature_of_process(ProcessType::MMAP) == Feature::MMAP);
}

TEST_CASE("FIRING and DISABLED routing map to no registry name") {
    CHECK(feature_of_routing(RoutingStrategy::FIRING) == Feature::COUNT);
    CHECK(feature_of_routing(RoutingStrategy::DISABLED) == Feature::COUNT);
    CHECK(feature_of_routing(RoutingStrategy::PROB) == Feature::RoutingStrategy_PROB);
}

TEST_CASE("FORK and NONE scheduling map to no registry name, HOL covers FCFSPRIO") {
    CHECK(feature_of_sched(SchedStrategy::FORK) == Feature::COUNT);
    CHECK(feature_of_sched(SchedStrategy::NONE) == Feature::COUNT);
    CHECK(feature_of_sched(SchedStrategy::HOL) == Feature::SchedStrategy_HOL);
    CHECK(feature_of_sched(SchedStrategy::OI) == Feature::SchedStrategy_OI);
    CHECK(feature_of_sched(SchedStrategy::REF) == Feature::SchedStrategy_REF);
}

TEST_CASE("the server section follows the discipline, as Queue.m assigns it") {
    CHECK(feature_of_server_section(SchedStrategy::PS) == Feature::SharedServer);
    CHECK(feature_of_server_section(SchedStrategy::DPS) == Feature::SharedServer);
    CHECK(feature_of_server_section(SchedStrategy::INF) == Feature::InfiniteServer);
    CHECK(feature_of_server_section(SchedStrategy::FCFS) == Feature::Server);
    CHECK(feature_of_server_section(SchedStrategy::SIRO) == Feature::Server);
    // PreemptiveServer and PollingServer are not registry names.
    CHECK(feature_of_server_section(SchedStrategy::LCFSPR) == Feature::COUNT);
    CHECK(feature_of_server_section(SchedStrategy::POLLING) == Feature::COUNT);
}

TEST_CASE("a specialization falls back to its generalization, one way only") {
    // Cox2 and Trace name a SPECIAL CASE of another registry entry rather than a
    // capability of their own, so the recorder marks the specific name and the
    // gate resolves it against the general one. Without the fallback, marking
    // the specific name would refuse the model at every solver that declares
    // only the general one -- which is every solver accepting it today.
    CHECK(feature_generalization(Feature::Cox2) == Feature::Coxian);
    CHECK(feature_generalization(Feature::Trace) == Feature::Replayer);
    // Everything else stands on its own; nothing else may acquire a fallback
    // silently, since a fallback WIDENS what a declared set accepts.
    for (std::size_t i = 0; i < feature_count(); ++i) {
        const Feature f = static_cast<Feature>(i);
        if (f == Feature::Cox2 || f == Feature::Trace) continue;
        INFO("feature = ", feature_name(f));
        CHECK(feature_generalization(f) == Feature::COUNT);
    }

    FeatureSet used;
    used.set(Feature::Cox2);
    FeatureSet general;
    general.set(Feature::Coxian);
    CHECK(feature_set_supports("s", general, used).ok);

    // One way only: declaring the SPECIAL case does not buy the general one, so
    // a solver restricted to two phases still refuses a five-phase Coxian.
    FeatureSet specific;
    specific.set(Feature::Cox2);
    FeatureSet used_general;
    used_general.set(Feature::Coxian);
    CHECK_FALSE(feature_set_supports("s", specific, used_general).ok);

    // And neither declared is still a refusal, naming the specific feature.
    const SupportResult r = feature_set_supports("s", FeatureSet(), used);
    CHECK_FALSE(r.ok);
    CHECK(r.missing.size() == 1);
    CHECK(r.missing[0] == Feature::Cox2);
}

TEST_CASE("a two-phase Coxian is what the Cox2 entry names") {
    // MATLAB has no Cox2 OBJECT -- Cox2.fitMeanAndSCV returns a Coxian -- so the
    // phase count is the only reading of the entry the four codebases can share,
    // and the one under which it is reachable at all. Dist::coxian already types
    // the two-phase case as COX2, so the derivation needs no special case.
    SN sn;
    sn.add_class(closed_class("c1", 1.0));
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.set_service(q, 1, Dist::coxian(std::vector<double>{2.0, 3.0},
                                      std::vector<double>{0.4, 1.0}));
    const FeatureSet u2 = used_lang_features(sn);
    CHECK(u2.has(Feature::Cox2));
    CHECK_FALSE(u2.has(Feature::Coxian));

    SN sn3;
    sn3.add_class(closed_class("c1", 1.0));
    const std::size_t q3 = add_st(sn3, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn3.set_service(q3, 1, Dist::coxian(std::vector<double>{2.0, 3.0, 4.0},
                                        std::vector<double>{0.3, 0.5, 1.0}));
    const FeatureSet u3 = used_lang_features(sn3);
    CHECK(u3.has(Feature::Coxian));
    CHECK_FALSE(u3.has(Feature::Cox2));

    // The gating a Coxian-declaring solver applies is unchanged by the split:
    // neither model is refused FOR ITS SERVICE PROCESS. Asserted on the missing
    // list rather than on ok, so that the structural features these hand-built
    // structs also carry (routing, sections) do not decide the outcome.
    FeatureSet coxian_only;
    coxian_only.set(Feature::Coxian);
    for (const FeatureSet* u : {&u2, &u3}) {
        const std::vector<Feature> missing = feature_set_supports("s", coxian_only, *u).missing;
        INFO("phases = ", (u == &u2 ? 2 : 3));
        CHECK(std::find(missing.begin(), missing.end(), Feature::Cox2) == missing.end());
        CHECK(std::find(missing.begin(), missing.end(), Feature::Coxian) == missing.end());
    }
}

TEST_CASE("replacement strategies map one to one") {
    CHECK(feature_of_replacement(ReplacementStrategy::RR) == Feature::ReplacementStrategy_RR);
    CHECK(feature_of_replacement(ReplacementStrategy::QLRU) == Feature::ReplacementStrategy_QLRU);
    CHECK(feature_of_replacement(ReplacementStrategy::HLRU) == Feature::ReplacementStrategy_HLRU);
}

// ---------------------------------------------------------------------------
// used_lang_features
// ---------------------------------------------------------------------------

TEST_CASE("baseline closed queueing model") {
    const SN sn = baseline();
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::ClosedClass, Feature::Queue, Feature::Buffer, Feature::Dispatcher,
                    Feature::Server, Feature::SchedStrategy_FCFS, Feature::RoutingStrategy_PROB,
                    Feature::Exp}));
    // Nothing the model does not have.
    CHECK_FALSE(u.has(Feature::OpenClass));
    CHECK_FALSE(u.has(Feature::Region));
    CHECK_FALSE(u.has(Feature::Fork));
    CHECK_FALSE(u.has(Feature::Retrial));
    CHECK_FALSE(u.has(Feature::SharedServer));
    CHECK_FALSE(u.has(Feature::InfiniteServer));
    // SelfLoopingClass is in the registry but getUsedLangFeatures never emits it.
    CHECK_FALSE(u.has(Feature::SelfLoopingClass));
}

TEST_CASE("a Delay carries Delay and InfiniteServer, never Queue") {
    SN sn;
    sn.add_class(closed_class("c1", 1.0));
    const std::size_t d = add_st(sn, "Think", NodeType::Delay, SchedStrategy::INF);
    sn.set_service(d, 1, Dist::erlang(2.0, 3));
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::Delay, Feature::InfiniteServer, Feature::SchedStrategy_INF,
                    Feature::Erlang}));
    CHECK_FALSE(u.has(Feature::Queue));
    CHECK_FALSE(u.has(Feature::Server));
}

TEST_CASE("a PS queue takes SharedServer, not Server") {
    SN sn;
    sn.add_class(closed_class("c1", 1.0));
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::PS);
    sn.set_service(q, 1, Dist::exp_mean(0.5));
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::SharedServer, Feature::SchedStrategy_PS}));
    CHECK_FALSE(u.has(Feature::Server));
}

TEST_CASE("an open model reports Source, Sink and the arrival distribution") {
    SN sn;
    sn.add_class(open_class("o1"));
    const std::size_t src = add_st(sn, "Src", NodeType::Source, SchedStrategy::EXT);
    sn.sourceIdx = src;
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.set_service(src, 1, Dist::hyperexp(0.4, 1.0, 3.0));
    sn.set_service(q, 1, Dist::exp_mean(0.2));
    sn.sinkNode = sn.add_node("Sink", NodeType::Sink, false);
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::OpenClass, Feature::Source, Feature::RandomSource,
                    Feature::ServiceTunnel, Feature::Dispatcher, Feature::Sink,
                    Feature::HyperExp, Feature::Exp}));
    CHECK_FALSE(u.has(Feature::ClosedClass));
    // The Source branch registers no discipline: MATLAB never emits SchedStrategy_EXT.
    CHECK_FALSE(u.has(Feature::SchedStrategy_EXT));
}

TEST_CASE("a disabled (station, class) pair contributes nothing") {
    SN sn;
    sn.add_class(closed_class("c1", 1.0));
    sn.add_class(closed_class("c2", 1.0));
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.set_service(q, 1, Dist::exp_mean(1.0));
    sn.set_service(q, 2, Dist::pareto(3.0, 1.0));
    sn.nodes[sn.station_to_node[q - 1] - 1].routing[1] = RoutingStrategy::JSQ;
    sn.disabled.assign(1, std::vector<bool>{false, true});
    const FeatureSet u = used_lang_features(sn);
    CHECK(u.has(Feature::Exp));
    CHECK_FALSE(u.has(Feature::Pareto));
    CHECK_FALSE(u.has(Feature::RoutingStrategy_JSQ));
    CHECK(u.has(Feature::RoutingStrategy_PROB));
}

TEST_CASE("a finite capacity region is reported") {
    SN sn = baseline();
    CHECK_FALSE(used_lang_features(sn).has(Feature::Region));
    SN::Region reg;
    reg.members = {true};
    reg.cap = {{-1.0, 4.0}};
    sn.regions.push_back(reg);
    CHECK(used_lang_features(sn).has(Feature::Region));
}

TEST_CASE("a retrial orbit is reported from retrialparam and from the drop rule") {
    SN sn = baseline();
    CHECK_FALSE(used_lang_features(sn).has(Feature::Retrial));
    SN rp = sn;
    RetrialParam<double> par;
    par.retrial_proc.assign(1, Dist::exp_mean(1.0));
    par.retrial_rate.assign(1, 2.0);
    par.max_attempts.assign(1, 0);
    rp.retrialparam[1] = par;
    CHECK(used_lang_features(rp).has(Feature::Retrial));

    SN dr = sn;
    dr.droprule.assign(1, std::vector<DropStrategy>{DropStrategy::RETRIAL});
    CHECK(used_lang_features(dr).has(Feature::Retrial));

    SN wq = sn;
    wq.droprule.assign(1, std::vector<DropStrategy>{DropStrategy::WAITQ});
    CHECK_FALSE(used_lang_features(wq).has(Feature::Retrial));
}

TEST_CASE("load dependence and class dependence come from the station handles") {
    SN sn = baseline();
    CHECK_FALSE(used_lang_features(sn).has(Feature::LoadDependence));
    CHECK_FALSE(used_lang_features(sn).has(Feature::ClassDependence));
    sn.stations[0].lldscaling = {1.0, 1.8, 2.4};
    sn.stations[0].cdscaling = [](const std::vector<double>& n) { return n; };
    const FeatureSet u = used_lang_features(sn);
    CHECK(u.has(Feature::LoadDependence));
    CHECK(u.has(Feature::ClassDependence));
    // JointDependence has no NetworkStruct field and can never be emitted.
    CHECK_FALSE(u.has(Feature::JointDependence));
}

TEST_CASE("G-network signal classes report kind, type, batch removal and policy") {
    SN sn;
    sn.add_class(closed_class("c1", 2.0));
    sn.add_class(open_class("sig"));
    const std::size_t q = add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.set_service(q, 1, Dist::exp_mean(1.0));
    sn.issignal = {false, true};
    sn.signaltype = {SignalType::REPLY, SignalType::NEGATIVE};
    sn.signalrempolicy = {RemovalPolicy::RANDOM, RemovalPolicy::LCFS};
    sn.signalremdist = {{}, {0.0, 0.5, 0.5}};
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::OpenSignal, Feature::SignalType_NEGATIVE,
                    Feature::SignalBatchRemoval, Feature::SignalRemovalPolicy}));
    CHECK_FALSE(u.has(Feature::ClosedSignal));
    // The non-signal class must not contribute its signaltype entry.
    CHECK_FALSE(u.has(Feature::SignalType_REPLY));
}

TEST_CASE("a closed signal class is reported as ClosedSignal") {
    SN sn;
    sn.add_class(closed_class("sig", 1.0));
    add_st(sn, "Q", NodeType::Queue, SchedStrategy::FCFS);
    sn.issignal = {true};
    sn.signaltype = {SignalType::CATASTROPHE};
    sn.signalrempolicy = {RemovalPolicy::RANDOM};
    sn.signalremdist = {{}};
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::ClosedSignal, Feature::SignalType_CATASTROPHE}));
    CHECK_FALSE(u.has(Feature::OpenSignal));
    CHECK_FALSE(u.has(Feature::SignalBatchRemoval));
    CHECK_FALSE(u.has(Feature::SignalRemovalPolicy));
}

TEST_CASE("a Cache reports its replacement policy and its retrieval system") {
    SN sn;
    sn.add_class(closed_class("c1", 1.0));
    const std::size_t cn = sn.add_node("C", NodeType::Cache, true);
    CacheParam<double> par;
    par.nitems = 4;
    par.itemcap = {2};
    par.replacestrat = ReplacementStrategy::LRU;
    sn.nodeparam[cn] = par;
    FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::Cache, Feature::CacheClassSwitcher, Feature::Buffer,
                    Feature::Dispatcher, Feature::ReplacementStrategy_LRU}));
    CHECK_FALSE(u.has(Feature::CacheRetrieval));

    sn.nodeparam[cn].retrieval_capacity = 2;
    u = used_lang_features(sn);
    CHECK(u.has(Feature::CacheRetrieval));
}

TEST_CASE("fork and join are reported from the node kinds and from sn.fj") {
    SN sn;
    sn.add_class(open_class("o1"));
    const std::size_t f = sn.add_node("F", NodeType::Fork, false);
    const std::size_t j = sn.add_node("J", NodeType::Join, true);
    FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::Fork, Feature::Forker, Feature::Join, Feature::Joiner,
                    Feature::ServiceTunnel}));

    // The fj pairing alone is enough, which is what survives the transform.
    SN only;
    only.add_class(open_class("o1"));
    only.fj.emplace_back(f, j);
    u = used_lang_features(only);
    CHECK(holds(u, {Feature::Fork, Feature::Forker, Feature::Join, Feature::Joiner}));
}

TEST_CASE("a class switch is reported from the switch matrix alone") {
    // This port never materialises a ClassSwitch node, so csmatrix is the only
    // trace it leaves and the gate must key on it.
    SN sn = baseline();
    CHECK_FALSE(used_lang_features(sn).has(Feature::ClassSwitch));
    sn.csmatrix[1] = Matrix<double>(1, 1, 1.0);
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::ClassSwitch, Feature::StatelessClassSwitcher}));
}

/**
 * A Router emits its routing strategy AND NOTHING ELSE.
 *
 * The previous expectation here was `Router, Buffer, Dispatcher, ServiceTunnel`
 * plus the strategy, which no reference codebase emits: getUsedLangFeatures.m:
 * 77-78, Network.java:3184-3190 and network.py:1983-1987 each record only
 * RoutingStrategy.toFeature(...). Since neither mva_feature_set nor
 * nc_feature_set declares Router -- correctly, SolverMVA.m:231-257 does not
 * either -- the extra names refused every dispatcher model those two solvers
 * solve in MATLAB and python. The assertion was changed because it encoded the
 * defect, not the reference.
 */
TEST_CASE("a Router reports its routing strategy and nothing else") {
    SN sn;
    sn.add_class(open_class("o1"));
    const std::size_t rn = sn.add_node("R", NodeType::Router, false);
    sn.nodes[rn - 1].routing.assign(1, RoutingStrategy::RROBIN);
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::RoutingStrategy_RROBIN}));
    CHECK_FALSE(u.has(Feature::Router));
    CHECK_FALSE(u.has(Feature::Buffer));
    CHECK_FALSE(u.has(Feature::Dispatcher));
    CHECK_FALSE(u.has(Feature::ServiceTunnel));
}

TEST_CASE("a Logger reports its log tunnel") {
    SN sn = baseline();
    sn.add_node("L", NodeType::Logger, true);
    const FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::Logger, Feature::LogTunnel}));
}

TEST_CASE("an SPN reports Place, Transition and only finite inhibitor arcs") {
    SN sn;
    sn.add_class(closed_class("tok", 3.0));
    add_st(sn, "P1", NodeType::Place, SchedStrategy::INF);
    const std::size_t tn = sn.add_node("T1", NodeType::Transition, true);
    TransitionParam<double> par;
    par.nmodes = 1;
    par.modenames = {"m1"};
    // One (node x class) matrix per mode; this net is single class.
    par.enabling.assign(1, Matrix<double>(1, 1, 1.0));
    par.firing.assign(1, Matrix<double>(1, 1, 1.0));
    const double inf = std::numeric_limits<double>::infinity();
    par.inhibiting.assign(1, Matrix<double>(1, 1, inf));
    sn.transparam[tn] = par;

    FeatureSet u = used_lang_features(sn);
    CHECK(holds(u, {Feature::Place, Feature::Storage, Feature::Linkage, Feature::Transition,
                    Feature::Enabling, Feature::Timing, Feature::Firing}));
    // An all-infinite inhibiting row is "never blocks": a plain SPN must not be
    // gated out of the solvers that lack inhibitor arcs.
    CHECK_FALSE(u.has(Feature::Inhibiting));
    // QueueingPlace has no NetworkStruct field and can never be emitted.
    CHECK_FALSE(u.has(Feature::QueueingPlace));

    sn.transparam[tn].inhibiting.assign(1, Matrix<double>(1, 1, 2.0));
    u = used_lang_features(sn);
    CHECK(u.has(Feature::Inhibiting));
}

TEST_CASE("the features NetworkStruct cannot express are never emitted") {
    // Each of these is a real MATLAB feature with no struct field; a model using
    // one reaches a C++ solver unflagged, and the fix is a field, not a proxy.
    SN sn = baseline();
    sn.stations[0].cdscaling = [](const std::vector<double>& n) { return n; };
    sn.stations[0].lldscaling = {1.0, 2.0};
    const FeatureSet u = used_lang_features(sn);
    CHECK_FALSE(u.has(Feature::JointDependence));
    CHECK_FALSE(u.has(Feature::SetupDelayOff));
    CHECK_FALSE(u.has(Feature::Reneging));
    CHECK_FALSE(u.has(Feature::Balking));
    CHECK_FALSE(u.has(Feature::Breakdown));
    CHECK_FALSE(u.has(Feature::QueueingPlace));
    CHECK_FALSE(u.has(Feature::BatchArrival));
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

TEST_CASE("feature_gate passes a model the solver declares") {
    const SN sn = baseline();
    FeatureSet declared;
    declared.set({Feature::ClosedClass, Feature::Queue, Feature::Buffer, Feature::Dispatcher,
                  Feature::Server, Feature::SchedStrategy_FCFS, Feature::RoutingStrategy_PROB,
                  Feature::Exp});
    CHECK_NOTHROW(feature_gate("SolverMVA", declared, sn));
}

TEST_CASE("feature_gate refuses a region by name") {
    SN sn = baseline();
    SN::Region reg;
    reg.members = {true};
    reg.cap = {{-1.0, 4.0}};
    sn.regions.push_back(reg);
    FeatureSet declared;
    declared.set({Feature::ClosedClass, Feature::Queue, Feature::Buffer, Feature::Dispatcher,
                  Feature::Server, Feature::SchedStrategy_FCFS, Feature::RoutingStrategy_PROB,
                  Feature::Exp});
    CHECK_THROWS_AS(feature_gate("SolverFluid", declared, sn), UnsupportedError);
    try {
        feature_gate("SolverFluid", declared, sn);
        FAIL("feature_gate did not refuse a finite capacity region");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()) ==
              "SolverFluid: this model uses a Finite Capacity Region, which this solver does not "
              "support");
    }
}

// ---------------------------------------------------------------------------
// The declared side: the seven per-solver sets
//
// Each case pins the exact declared CARDINALITY, so that a name silently added
// or lost cannot pass, then a few must-have and must-not-have members that
// carry the meaning: the MATLAB names that define the solver's envelope and the
// C++ narrowings whose justifying refusal sites are named in the header.
// ---------------------------------------------------------------------------

TEST_CASE("mva_feature_set is SolverMVA.getFeatureSet") {
    const FeatureSet f = line::qn::mva_feature_set("default");
    // 2026-08-04: 59 -> 60 with Cox2. MATLAB has no Cox2 FEATURE to copy --
    // `Cox2` there is a static factory that returns a `Coxian`, so a model never
    // reports the name -- while this port carries ProcessType::COX2 as its own
    // enumerator and feature_of_process maps it to Feature::Cox2. Withholding it
    // refused a two-phase Coxian written as Cox2 and admitted the identical
    // distribution written as Coxian.
    // 2026-08-14: 60 -> 61 with JoinPartial. The mmt fork-join transform now
    // reads joinRequired and synchronises on E[X_(k)] instead of E[X_(n)], so
    // the quorum is approximated rather than silently ignored.
    // 2026-08-16: 61 -> 66 with the five size-based disciplines. MATLAB's
    // getFeatureSet declares them and mva_dispatch.h branch 3 serves them; the
    // earlier removal rested on a MATLAB state that has since changed.
    // 2026-08-29: 66 -> 67 with SchedStrategy_FCFSPRPRIO, which SolverMVA.m:321
    // and SolverMVA.java:187 both declare and whose preemptive-priority arm the
    // runner offers when a station actually uses it.
    // 2026-09-05: 67 -> 69 with MultiServer and FiniteCapacity. SolverMVA.m
    // declares MultiServer in the base envelope and grants FiniteCapacity to
    // 'default'/'sqd'; before this the C++ set named neither, so the newly
    // recorded features refused every M/M/c and every binding buffer.
    CHECK(f.list().size() == 69);
    CHECK(f.has(Feature::SchedStrategy_FCFSPRPRIO));
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::JoinPartial));
    CHECK(f.has(Feature::ClosedClass));
    CHECK(f.has(Feature::SelfLoopingClass));
    CHECK(f.has(Feature::Cache));
    CHECK(f.has(Feature::CacheRetrieval));
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Joiner));
    // MATLAB's getFeatureSet declares the five size-based disciplines, and
    // mva_dispatch.h branch 3 serves them; the 2026-07-29 assertion that it
    // declared none of them expired when SolverMVA.m listed them.
    CHECK(f.has(Feature::SchedStrategy_SRPT));
    CHECK(f.has(Feature::SchedStrategy_PSJF));
    CHECK(f.has(Feature::SchedStrategy_FB));
    CHECK(f.has(Feature::SchedStrategy_LRPT));
    CHECK(f.has(Feature::SchedStrategy_SETF));
    CHECK(f.has(Feature::SchedStrategy_SJF));
    CHECK(f.has(Feature::SchedStrategy_POLLING));
    CHECK(f.has(Feature::LoadDependence));
    CHECK(f.has(Feature::ClassDependence));
    CHECK(f.has(Feature::MMAP));
    CHECK(f.has(Feature::BMAP));
    // The five constructs MVA neither models nor refuses today.
    CHECK_FALSE(f.has(Feature::Region));
    CHECK_FALSE(f.has(Feature::Retrial));
    CHECK_FALSE(f.has(Feature::Logger));
    CHECK_FALSE(f.has(Feature::Router));
    CHECK_FALSE(f.has(Feature::Place));
    // Not in the MATLAB base set; only the rqna delta adds them.
    CHECK_FALSE(f.has(Feature::MAP));
    CHECK_FALSE(f.has(Feature::MMPP2));
    CHECK_FALSE(f.has(Feature::RAP));
    CHECK_FALSE(f.has(Feature::PH));
    CHECK_FALSE(f.has(Feature::Gamma));
    CHECK_FALSE(f.has(Feature::SchedStrategy_GPS));
    CHECK_FALSE(f.has(Feature::OpenSignal));
}

TEST_CASE("mva_feature_set applies the rqna per-method delta") {
    const FeatureSet f = line::qn::mva_feature_set("rqna");
    // base 67, plus MAP, MMPP2 and RAP (MMAP is already in), minus the two
    // closed classes, the three dependence flags rqna is on the load-dependent
    // refusal list for, SchedStrategy_OI and _PAS (solver_mva_oi is reached
    // under 'default'/'exact' only) and the five fork-join names (a Join carries
    // no service process, so the index-of-dispersion curve rqna reads off every
    // station does not exist for it): 67 + 3 - 12 = 58.
    CHECK(f.list().size() == 58);
    CHECK_FALSE(f.has(Feature::LoadDependence));
    CHECK_FALSE(f.has(Feature::Join));
    CHECK_FALSE(f.has(Feature::SchedStrategy_OI));
    CHECK(f.has(Feature::MAP));
    CHECK(f.has(Feature::MMPP2));
    CHECK(f.has(Feature::MMAP));
    CHECK(f.has(Feature::RAP));
    CHECK(f.has(Feature::OpenClass));
    CHECK_FALSE(f.has(Feature::ClosedClass));
    CHECK_FALSE(f.has(Feature::SelfLoopingClass));
}

TEST_CASE("mva_feature_set applies the qna per-method delta") {
    const FeatureSet f = line::qn::mva_feature_set("qna");
    // base 67 minus the two closed classes, plus RoutingStrategy_RROBIN
    // (SolverMVA.m:149: round-robin dispatching enters as a deterministic
    // traffic split, npfqn_traffic_split_rr, for this method only), minus the
    // three dependence flags qna is on the load-dependent refusal list for,
    // SchedStrategy_OI and _PAS, and the fourteen disciplines outside solver_qna's
    // INF/PS/FCFS station loop: 67 + 1 - 20 = 48. Fork/Join STAY -- that loop has
    // an explicit Join arm, which is what separates qna from rqna here.
    // 2026-09-05: 48 -> 49. qna keeps MultiServer from the base envelope; it is
    // on MATLAB's withdrawal list for rqna, not for qna.
    CHECK(f.list().size() == 49);
    CHECK(f.has(Feature::Join));
    CHECK_FALSE(f.has(Feature::SchedStrategy_SIRO));
    CHECK_FALSE(f.has(Feature::LoadDependence));
    CHECK(f.has(Feature::OpenClass));
    CHECK(f.has(Feature::RoutingStrategy_RROBIN));
    CHECK_FALSE(f.has(Feature::ClosedClass));
    CHECK_FALSE(f.has(Feature::SelfLoopingClass));
    // The rqna-only non-renewal additions stay out under qna.
    CHECK_FALSE(f.has(Feature::MAP));
    CHECK_FALSE(f.has(Feature::MMPP2));
    CHECK_FALSE(f.has(Feature::RAP));
}

TEST_CASE("mva_feature_set gives every other method the base envelope less OI/PAS") {
    // 'default' is the ONLY name that keeps the whole envelope, because it is the
    // only one (with 'exact') that mva_dispatch routes to solver_mva_oi: an
    // order-independent or pass-and-swap station is served there and refused by
    // name everywhere else, so advertising the two disciplines under any other
    // method would name a pair that cannot run. 'exact' drops the class- and
    // joint-dependent scalings on top, its recursion carrying load dependence
    // alone. Nothing else may take a delta -- a new one appearing here is the
    // regression this case exists to catch.
    const FeatureSet base = line::qn::mva_feature_set("default");
    for (const char* m : {"amva", "mm1", "gig1.kingman", "sqd", ""}) {
        CAPTURE(m);
        const std::string name(m);
        FeatureSet want = base;
        want.unset(Feature::SchedStrategy_OI);
        want.unset(Feature::SchedStrategy_PAS);
        // 2026-09-05, with FiniteCapacity and MultiServer now registered names:
        // 'sqd' is the ONE name here that keeps FiniteCapacity, because it names
        // solver_sqd, the only ladder that HONOURS the buffers rather than
        // solving them away (SolverMVA.m:393-395). 'mm1' and 'gig1.kingman' are
        // single-server closed forms and drop MultiServer on top, as they do in
        // MATLAB's withdrawal list. These are the intended deltas; any OTHER one
        // appearing here is still the regression this case exists to catch.
        if (name != "sqd") want.unset(Feature::FiniteCapacity);
        if (name == "mm1" || name == "gig1.kingman") want.unset(Feature::MultiServer);
        CHECK(line::qn::mva_feature_set(m).list() == want.list());
    }
    FeatureSet want_exact = base;
    want_exact.unset(Feature::ClassDependence);
    want_exact.unset(Feature::JointDependence);
    // 'exact' is the one name MATLAB withholds FiniteCapacity from on the
    // single-station M/M/1/K shape, so it drops it here too. It keeps OI/PAS.
    want_exact.unset(Feature::FiniteCapacity);
    CHECK(line::qn::mva_feature_set("exact").list() == want_exact.list());
    CHECK(line::qn::mva_feature_set("default").has(Feature::SchedStrategy_OI));
}

TEST_CASE("nc_feature_set is SolverNC.getFeatureSet") {
    const FeatureSet f = line::qn::nc_feature_set("default");
    // 49 -> 50 when Geometric was declared for the discrete-time route
    // (opt.slotted, solver_nc_dt); the continuous routes take it by its two
    // moments like any other renewal law.
    // 2026-08-04: 50 -> 51 with Cox2, for the reason given on mva_feature_set.
    // 2026-08-09: 51 -> 52 with RoutingStrategy_SDR. `solver_nc` intercepts an
    // sdr model before the convolution and MVA analyzers and evaluates the
    // Krzesinski eq. (16) product form, so the claim is honoured.
    // 2026-08-14: 52 -> 53 with JoinPartial, for the reason given on
    // mva_feature_set -- nc shares the mmt fork-join transform.
    // 2026-09-04: 53 -> 62. SchedStrategy_DPS, which SolverNC.m declared on
    // 2026-08-29 with Morrison's asymptotics as the NC default, plus the eight
    // Petri-net constructs the 'rec' route (solver_nc_spn) reads. Both were
    // ported to the JAR and python in their own commits and this set was the
    // one left behind, so the CHECK_FALSE on DPS below was asserting a
    // divergence rather than the reference.
    // 2026-09-05: 62 -> 64 with MultiServer (base) and FiniteCapacity, which
    // SolverNC.m grants to 'mem'/'default'/'exact'.
    CHECK(f.list().size() == 64);
    CHECK(f.has(Feature::JoinPartial));
    // Morrison's closed think+DPS shape only (nc_is_dps_model); every other DPS
    // model is refused imperatively, as MATLAB's runAnalyzer does.
    CHECK(f.has(Feature::SchedStrategy_DPS));
    // The 'rec' route walks the reachable marking set in a decision diagram.
    CHECK(f.has(Feature::Place));
    CHECK(f.has(Feature::Transition));
    CHECK(f.has(Feature::Firing));
    // A queueing Place stays out: spn_pf decides the product-form class.
    CHECK_FALSE(f.has(Feature::QueueingPlace));
    CHECK(f.has(Feature::Geometric));
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::RoutingStrategy_SDR));
    // Region is declared DELIBERATELY: the loss-network case is exact and the
    // queueing-station case stays as solver_nc_runner.h:315-332.
    CHECK(f.has(Feature::Region));
    CHECK(f.has(Feature::ClosedClass));
    CHECK(f.has(Feature::OpenClass));
    CHECK(f.has(Feature::Cache));
    CHECK(f.has(Feature::CacheRetrieval));
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::SchedStrategy_OI));
    CHECK(f.has(Feature::LoadDependence));
    CHECK(f.has(Feature::ClassDependence));
    // MATLAB's NC refuses these disciplines and distributions; ground truth wins.
    CHECK_FALSE(f.has(Feature::SchedStrategy_GPS));
    CHECK_FALSE(f.has(Feature::SchedStrategy_HOL));
    CHECK_FALSE(f.has(Feature::Gamma));
    CHECK_FALSE(f.has(Feature::Pareto));
    CHECK_FALSE(f.has(Feature::Weibull));
    CHECK_FALSE(f.has(Feature::MMAP));
    CHECK_FALSE(f.has(Feature::BMAP));
    CHECK_FALSE(f.has(Feature::Replayer));
    CHECK_FALSE(f.has(Feature::ReplacementStrategy_LRU));
    CHECK_FALSE(f.has(Feature::Retrial));
    CHECK_FALSE(f.has(Feature::Router));
    CHECK_FALSE(f.has(Feature::Logger));
    CHECK(line::qn::nc_feature_set("exact").list() == f.list());
}

TEST_CASE("ctmc_feature_set is SolverCTMC.getFeatureSet in full") {
    const FeatureSet f = line::qn::ctmc_feature_set("default");
    // 104 MATLAB names, none withheld.
    // 2026-07-29: 99 -> 104 when the five remaining preempt disciplines were
    // declared in SolverCTMC.m:159, the generator having always supported them.
    // 2026-07-30: the five names this set used to withhold -- Fork, Forker, Join,
    // Joiner and ClassDependence -- are now declared, `fj_tag` and `cd_factor`
    // having landed. The four fork-join names are asserted PRESENT below for the
    // same reason they were asserted absent before: a featset name is a claim.
    // 2026-07-31: 104 -> 103. RoutingStrategy_RROBIN was WITHDRAWN when
    // refresh_routing started expanding round robin (uniformly, for QNA/MNA):
    // this port's generator carries no round-robin pointer, so a declared
    // RROBIN would now answer a random-routing model instead of being refused
    // when the struct is built.
    // 2026-08-04: 103 -> 104 with GlobalDependence, phi(n) over the full
    // population matrix, which solver_ctmc.h tabulates per state.
    // 2026-08-04: 104 -> 105 with Cox2, for the reason given on mva_feature_set.
    // 2026-08-09: 105 -> 106 with RoutingStrategy_SDR. The generator reads the
    // eq. (10) probabilities per state through `qn::rt_state`, so the route is
    // exact rather than a uniform stand-in.
    // 2026-08-16: 106 -> 107. RoutingStrategy_RROBIN is BACK, and the 2026-07-31
    // withdrawal above is what changed: the port now carries the dispatch
    // pointer. `refresh_local_vars` allocates it, `from_marginal_node` and
    // `append_local_vars` enumerate it, `after_event_station` and
    // `after_event_router` advance it on every departure, and the generator
    // reads the destination out of the ACTIVE node's post-departure row. The
    // declaration is a claim about a dispatcher this generator now carries, not
    // about a uniform stand-in.
    // 2026-09-05: 107 -> 109 with MultiServer and FiniteCapacity, both in the
    // CTMC base envelope; cftp/cftp.approx/mdd withdraw FiniteCapacity, as in
    // MATLAB, since neither construction carries a buffer.
    CHECK(f.list().size() == 109);
    CHECK(f.has(Feature::GlobalDependence));
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::RoutingStrategy_SDR));
    CHECK(f.has(Feature::RoutingStrategy_RROBIN));
    // WRROBIN is a dispatcher too, and shares the pointer machinery: its slot
    // holds a POSITION in the weighted cycle rather than a destination, because
    // a repeated outlink must advance once per repetition. JSQ and SQ remain
    // declared on the older ground -- the reference spreads them uniformly
    // (getRoutingMatrix.m:117) and refresh_routing does the same -- so those two
    // are still a claim about the reference's approximation.
    CHECK(f.has(Feature::RoutingStrategy_WRROBIN));
    CHECK(f.has(Feature::SchedStrategy_LCFSPI));
    CHECK(f.has(Feature::SchedStrategy_FCFSPR));
    CHECK(f.has(Feature::Region));
    CHECK(f.has(Feature::Retrial));
    CHECK(f.has(Feature::Router));
    CHECK(f.has(Feature::Place));
    CHECK(f.has(Feature::Transition));
    CHECK(f.has(Feature::Inhibiting));
    CHECK(f.has(Feature::OpenSignal));
    CHECK(f.has(Feature::SignalType_REPLY));
    CHECK(f.has(Feature::SchedStrategy_LCFSPRPRIO));
    CHECK(f.has(Feature::RoutingStrategy_JSQ));
    CHECK(f.has(Feature::ReplacementStrategy_QLRU));
    CHECK(f.has(Feature::LoadDependence));
    CHECK(f.has(Feature::JointDependence));
    // Fork-join is solved on the tag-augmented copy built by `fj_tag`.
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Forker));
    CHECK(f.has(Feature::Join));
    CHECK(f.has(Feature::Joiner));
    // Class dependence is applied by `cd_factor` in the event layer.
    CHECK(f.has(Feature::ClassDependence));
    CHECK_FALSE(f.has(Feature::Logger));
    CHECK_FALSE(f.has(Feature::SchedStrategy_EXT));
    CHECK(line::qn::ctmc_feature_set("gpu").list() == f.list());
}

TEST_CASE("ssa_feature_set is the SOLVER's reach, both engines together") {
    // The set is what `solver_ssa` under `default` can answer, not what one
    // engine can: the NRM runs the models it is eligible for and the serial
    // engine runs the rest, so a construct EITHER covers belongs here. AUTO
    // picks SolverSSA on the strength of it, which is why an entry is added only
    // once a test drives the construct end to end.
    const FeatureSet f = line::qn::ssa_feature_set("default");
    CHECK(f.has(Feature::Router));
    CHECK(f.has(Feature::Retrial));
    CHECK(f.has(Feature::Replayer));
    CHECK(f.has(Feature::SchedStrategy_EXT));
    CHECK(f.has(Feature::SchedStrategy_LCFSPR));
    CHECK(f.has(Feature::SchedStrategy_SEPT));
    CHECK(f.has(Feature::SchedStrategy_LPS));
    CHECK(f.has(Feature::OpenSignal));
    CHECK(f.has(Feature::SignalRemovalPolicy));
    CHECK(f.has(Feature::LoadDependence));
    // Reached through the serial engine, each with an end-to-end test behind it
    // (test_ssa_fj_fcr.cpp for the first four groups, test_ssa_serial.cpp for
    // the disciplines).
    CHECK(f.has(Feature::Cache));
    CHECK(f.has(Feature::CacheClassSwitcher));
    CHECK(f.has(Feature::CacheRetrieval));
    CHECK(f.has(Feature::Place));
    CHECK(f.has(Feature::Transition));
    CHECK(f.has(Feature::Storage));
    CHECK(f.has(Feature::Inhibiting));
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Joiner));
    CHECK(f.has(Feature::Region));
    CHECK(f.has(Feature::SchedStrategy_PAS));
    CHECK(f.has(Feature::SchedStrategy_OI));
    CHECK(f.has(Feature::SchedStrategy_POLLING));
    CHECK(f.has(Feature::SchedStrategy_PSPRIO));
    CHECK(f.has(Feature::SchedStrategy_DPSPRIO));
    CHECK(f.has(Feature::SchedStrategy_GPSPRIO));
    CHECK(f.has(Feature::SchedStrategy_LCFSPRPRIO));
    CHECK(f.has(Feature::SchedStrategy_FCFSPRPRIO));
    CHECK(f.has(Feature::ClassDependence));
    // 2026-08-16: THE DISPATCHERS ARE IN. The claim above -- that
    // `refresh_routing` rejects them and that no round-robin pointer is
    // resolved anywhere -- was true of neither by then: the expansion accepts
    // all four, the NRM keeps its own cursor, and the serial engine reads the
    // pointer the state now carries.
    CHECK(f.has(Feature::RoutingStrategy_RROBIN));
    CHECK(f.has(Feature::RoutingStrategy_WRROBIN));
    CHECK(f.has(Feature::RoutingStrategy_JSQ));
    CHECK(f.has(Feature::RoutingStrategy_SQ));
    CHECK_FALSE(f.has(Feature::Logger));
    // The serial engine RUNS these, and the set still withholds them: MATLAB's
    // `SolverSSA.getFeatureSet` does not declare them, and AUTO selects on this
    // set, so declaring them would diverge the two AUTOs by declaration alone.
    CHECK_FALSE(f.has(Feature::SchedStrategy_FCFSPR));
    CHECK_FALSE(f.has(Feature::SchedStrategy_LCFSPI));
    // The set does not depend on the method NAME: `nrm` names an estimator, not
    // a smaller model class, and a caller who spells it out on a model the NRM
    // cannot run is refused by name rather than by a featset that lied about it.
    CHECK(line::qn::ssa_feature_set("nrm").list() == f.list());
}

TEST_CASE("fluid_feature_set is SolverFLD.getFeatureSet, per METHOD") {
    // The declaration is now method-aware, as `@SolverFLD/getMethodFeatureSet` is:
    // the base set carries LoadDependence and GPS and every method subtracts what
    // its own drift cannot evaluate. `default` is the literal name here, which
    // MATLAB's getMethodFeatureSet also treats as "not the closing family" --
    // runAnalyzer resolves it before the gate ever sees it.
    const FeatureSet f = line::qn::fluid_feature_set("default");
    // 48 -> 50 with the time-inhomogeneous pair MAPt/PHt, which the reference
    // restricts in no method either: the first-order drifts integrate the
    // width-weighted nominal pair and `kp` integrates the schedule itself.
    // 50 -> 55 with the fork-join five, which `SolverFLD.m:200-202` declares as of
    // 2026-08-15: the MMT fixed point (`mva::fj_fixed_point`) drives a fluid inner
    // solve and emits only nodes the drift already carries.
    // 2026-09-05: 55 -> 57 with MultiServer and FiniteCapacity in the base
    // envelope; diffusion and mfq withdraw MultiServer, as in MATLAB.
    CHECK(f.list().size() == 57);
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Join));
    CHECK(f.has(Feature::Forker));
    CHECK(f.has(Feature::Joiner));
    CHECK(f.has(Feature::JoinPartial));
    CHECK(f.has(Feature::MAPt));
    CHECK(f.has(Feature::PHt));
    CHECK(f.has(Feature::Cache));
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::NHPP));
    CHECK(f.has(Feature::MAP));
    CHECK(f.has(Feature::Weibull));
    CHECK(f.has(Feature::ReplacementStrategy_SFIFO));
    CHECK(f.has(Feature::Replayer));
    CHECK(f.has(Feature::SchedStrategy_SIRO));
    // Deliberately undeclared upstream (SolverFLD.m:95-101, flow not conserved).
    CHECK_FALSE(f.has(Feature::CacheRetrieval));
    // The silent sites the transcription alone closes.
    CHECK_FALSE(f.has(Feature::Region));
    CHECK_FALSE(f.has(Feature::Router));
    CHECK_FALSE(f.has(Feature::Retrial));
    CHECK_FALSE(f.has(Feature::Place));
    CHECK_FALSE(f.has(Feature::OpenSignal));
    CHECK_FALSE(f.has(Feature::ClassDependence));
    CHECK_FALSE(f.has(Feature::SchedStrategy_POLLING));
    CHECK_FALSE(f.has(Feature::MMAP));
    CHECK_FALSE(f.has(Feature::PH));
    CHECK_FALSE(f.has(Feature::ReplacementStrategy_LRU));
    // Load dependence multiplies the scheduling share, which only the closing
    // family evaluates; HOL needs the class priorities, which only the mfq
    // priority branch reads; GPS needs a backlog PROBABILITY, which only the
    // second-order closure supplies.
    CHECK_FALSE(f.has(Feature::LoadDependence));
    CHECK_FALSE(f.has(Feature::SchedStrategy_HOL));
    CHECK_FALSE(f.has(Feature::SchedStrategy_GPS));
    CHECK(line::qn::fluid_feature_set("closing").has(Feature::LoadDependence));
    CHECK(line::qn::fluid_feature_set("minnormal").has(Feature::LoadDependence));
    CHECK(line::qn::fluid_feature_set("minnormal").has(Feature::SchedStrategy_GPS));
    CHECK(line::qn::fluid_feature_set("mfq").has(Feature::SchedStrategy_HOL));
    CHECK_FALSE(line::qn::fluid_feature_set("refined").has(Feature::SchedStrategy_GPS));
    CHECK_FALSE(line::qn::fluid_feature_set("matrix").has(Feature::LoadDependence));
    // The closing family has no branch for SIRO, LCFS or LCFSPR, and a station
    // with no branch is integrated as an INFINITE SERVER: the metric reader used
    // to catch LCFS/LCFSPR but the moment closure reads the event representation
    // directly and never reached that guard, and SIRO was accepted AS FCFS.
    for (const char* m : {"closing", "statedep", "softmin", "tbi", "minnormal", "refined", "dae"}) {
        const FeatureSet g = line::qn::fluid_feature_set(m);
        CHECK_FALSE(g.has(Feature::SchedStrategy_SIRO));
        CHECK_FALSE(g.has(Feature::SchedStrategy_LCFS));
        CHECK_FALSE(g.has(Feature::SchedStrategy_LCFSPR));
    }
    // `dae` is the ONE method that widens the set: a Finite Capacity Region is a
    // linear inequality on the state, which the DAE form carries as an algebraic
    // equation and no ODE method can carry at all. Every other method must go on
    // refusing it -- an ODE integrated through a cap it cannot see returns the
    // unconstrained answer silently -- so this pair of checks is what keeps the
    // widening from leaking.
    CHECK(line::qn::fluid_feature_set("dae").has(Feature::Region));
    // The Petri-net constructs are `dae`'s alone: it is the only fluid method
    // that states the P-invariants, an immediate transition's firing FLOW and a
    // bounded place as EQUATIONS (fluid_petri.h). Every other method would build
    // its drift from the station/class/phase encoding, where a Place contributes
    // no coordinate, and report zeros for the whole net without a warning.
    CHECK(line::qn::fluid_feature_set("dae").has(Feature::Place));
    CHECK(line::qn::fluid_feature_set("dae").has(Feature::Transition));
    CHECK(line::qn::fluid_feature_set("dae").has(Feature::Firing));
    CHECK_FALSE(line::qn::fluid_feature_set("closing").has(Feature::Place));
    CHECK_FALSE(line::qn::fluid_feature_set("minnormal").has(Feature::Transition));
    CHECK(line::qn::fluid_feature_set("fluid.dae").has(Feature::Region));
    for (const char* m : {"default", "closing", "matrix", "pnorm", "minnormal", "refined",
                          "statedep", "softmin", "tbi", "diffusion", "mfq", "kp", "rmf"}) {
        CAPTURE(m);
        CHECK_FALSE(line::qn::fluid_feature_set(m).has(Feature::Region));
    }
    // ... and it narrows in the other direction, for the reason `minnormal` keeps
    // GPS and it does not: DPS closes on the covariance BETWEEN a station's class
    // coordinates, where the DAE form carries one scalar variance per station.
    CHECK_FALSE(line::qn::fluid_feature_set("dae").has(Feature::SchedStrategy_DPS));
    CHECK_FALSE(line::qn::fluid_feature_set("dae").has(Feature::SchedStrategy_GPS));
    CHECK(line::qn::fluid_feature_set("closing").has(Feature::SchedStrategy_DPS));
    // The moment-closure envelope is inherited, load dependence included.
    CHECK(line::qn::fluid_feature_set("dae").has(Feature::LoadDependence));
    // `fluid.<name>` is the same method under its qualified spelling.
    CHECK(line::qn::fluid_feature_set("fluid.minnormal").list() ==
          line::qn::fluid_feature_set("minnormal").list());
}

TEST_CASE("ba_feature_set carries the Petri-net constructs of the spnlp family") {
    const FeatureSet f = line::qn::ba_feature_set("default");
    // solver_ba_spnlp is indexed by the MARKING rather than by demands, and
    // reads the enabling, inhibiting and firing arcs out of the node
    // parameters, so a Petri net is a model SolverBA can bound. SolverBA.m,
    // SolverBA.java and solver_ba.py all declare these eight; this set was the
    // one left behind, and because auto_methods.h consults it, spnlp was
    // unreachable through SolverAUTO here alone.
    CHECK(f.has(Feature::Place));
    CHECK(f.has(Feature::Transition));
    CHECK(f.has(Feature::Linkage));
    CHECK(f.has(Feature::Enabling));
    CHECK(f.has(Feature::Inhibiting));
    CHECK(f.has(Feature::Timing));
    CHECK(f.has(Feature::Firing));
    CHECK(f.has(Feature::Storage));
    // A place with an embedded queue has local state the relaxation has no
    // variable for, and solver_ba_spnlp refuses it by name.
    CHECK_FALSE(f.has(Feature::QueueingPlace));
    // The queueing side of the set is untouched by the widening.
    CHECK(f.has(Feature::Queue));
    CHECK(f.has(Feature::ClosedClass));
    CHECK(f.has(Feature::SchedStrategy_FCFS));
}

TEST_CASE("mam_feature_set widens SolverMAM.getFeatureSet by LoadDependence for default/ldqbd") {
    const FeatureSet f = line::qn::mam_feature_set("default");
    // the union of the four setTrue calls is 55 names, MINUS FCFSPRPRIO (no MAM
    // analyzer serves it, see the header note), PLUS LoadDependence (widening,
    // 'default'/'ldqbd' only) = 54 - 1 + 1... i.e. 54 names, plus LoadDependence.
    // 55 -> 57 with Geometric and DiscreteUniform, the two lattice laws the
    // slotted path (mam_dispatch.h branch -1, solver_mam_dt) takes; off the
    // lattice the continuous analyzers take them through their moment fits.
    // 2026-08-04: 57 -> 58 with Cox2, for the reason given on mva_feature_set.
    // 2026-08-19: 58 -> 53. The five G-network names (OpenSignal, ClosedSignal,
    // SignalType_NEGATIVE, SignalType_CATASTROPHE, SignalBatchRemoval) are now
    // declared by the RCAT methods ALONE: solver_mam_ag.h is the only analyzer
    // that reads issignal, so on any other method the union promised a model it
    // would have solved with the signals turned into ordinary customers.
    // 2026-09-05: 53 -> 55 with MultiServer and FiniteCapacity in the base
    // envelope; ldqbd and bgchain withdraw FiniteCapacity, as in MATLAB.
    CHECK(f.list().size() == 55);
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::Geometric));
    CHECK(f.has(Feature::DiscreteUniform));
    // solver_mam_ldqbd.h:152 consumes st.lldscaling, so refusing it would be false.
    CHECK(f.has(Feature::LoadDependence));
    CHECK(f.has(Feature::Retrial));
    CHECK(f.has(Feature::SetupDelayOff));
    CHECK(f.has(Feature::BMAP));
    CHECK(f.has(Feature::PH));
    CHECK(f.has(Feature::RAP));
    CHECK(f.has(Feature::DMAP));
    CHECK(f.has(Feature::ME));
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Joiner));
    // The G-network names belong to the RCAT methods; see the signal test below.
    CHECK_FALSE(f.has(Feature::OpenSignal));
    CHECK_FALSE(f.has(Feature::SignalBatchRemoval));
    // No MAM analyzer serves FCFSPRPRIO: solve_fcfs_station dispatches only FCFS
    // and HOL (solver_mam_basic.h), so declaring it would only let a model past
    // the gate to fail deeper, with a less specific message.
    CHECK_FALSE(f.has(Feature::SchedStrategy_FCFSPRPRIO));
    CHECK_FALSE(f.has(Feature::Region));
    CHECK_FALSE(f.has(Feature::Router));
    CHECK_FALSE(f.has(Feature::Logger));
    CHECK_FALSE(f.has(Feature::Cache));
    CHECK_FALSE(f.has(Feature::Replayer));
    CHECK_FALSE(f.has(Feature::ClassDependence));
    CHECK_FALSE(f.has(Feature::SchedStrategy_LCFS));
    CHECK_FALSE(f.has(Feature::SchedStrategy_LCFSPR));
    CHECK_FALSE(f.has(Feature::SchedStrategy_SIRO));
    CHECK_FALSE(f.has(Feature::SchedStrategy_DPS));
    CHECK_FALSE(f.has(Feature::SignalType_REPLY));
    // 'default' and 'ldqbd' both reach solver_mam_ldqbd, so both declare
    // LoadDependence. Since 2026-09-05 the two sets are no longer IDENTICAL:
    // ldqbd withdraws FiniteCapacity, as SolverMAM.m does, so it is exactly one
    // name short of 'default' and that name is the only difference.
    {
        const FeatureSet ld = line::qn::mam_feature_set("ldqbd");
        CHECK(ld.list().size() == f.list().size() - 1);
        CHECK_FALSE(ld.has(Feature::FiniteCapacity));
        CHECK(ld.has(Feature::MultiServer));
        CHECK(ld.has(Feature::LoadDependence));
    }
}

TEST_CASE("mam_feature_set withholds LoadDependence outside default/ldqbd") {
    // dec.source, dec.poisson and mna all reach solver_mam_basic or
    // solver_mna, neither of which reads st.lldscaling; declaring
    // LoadDependence there would let a load-dependent model past the gate to
    // be silently solved at nominal rates.
    for (const char* m : {"dec.source", "dec.poisson", "mna"}) {
        const std::string name(m);
        const FeatureSet f = line::qn::mam_feature_set(m);
        CHECK_FALSE(f.has(Feature::LoadDependence));
        CHECK_FALSE(f.has(Feature::SchedStrategy_FCFSPRPRIO));
        // 52 = the 53 of `default` minus LoadDependence; 'mna' is the one method
        // that resolves the round-robin split and adds it straight back. The
        // RCAT methods are no longer SolverMAM methods at all -- they moved to
        // SolverAG, whose envelope is ag_feature_set.
        // 2026-09-05: 53/52 -> 55/54 with MultiServer and FiniteCapacity, both
        // in the MAM base envelope and neither withdrawn by these three.
        CHECK(f.list().size() == (name == "mna" ? 55u : 54u));
    }
}

TEST_CASE("the G-network names belong to ag_feature_set and to no MAM method") {
    // solver_ag.h's build_rcat is the only code in LINE that reads `issignal`.
    // It moved to SolverAG with the RCAT methods, so the names moved with it:
    // any MAM method declaring them would promise what it cannot deliver, which
    // is an answer with every signal turned into an ordinary customer.
    for (const char* m : {"default", "inap", "inapplus", "inapinf", "exact"}) {
        const FeatureSet f = line::qn::ag_feature_set(m);
        CHECK(f.has(Feature::OpenSignal));
        CHECK(f.has(Feature::ClosedSignal));
        CHECK(f.has(Feature::SignalType_NEGATIVE));
        CHECK(f.has(Feature::SignalType_CATASTROPHE));
        CHECK(f.has(Feature::SignalBatchRemoval));
    }
    for (const char* m : {"default", "dec.source", "dec.poisson", "dec.mmap", "mna", "ldqbd",
                          "bgchain"}) {
        const FeatureSet f = line::qn::mam_feature_set(m);
        CHECK_FALSE(f.has(Feature::OpenSignal));
        CHECK_FALSE(f.has(Feature::ClosedSignal));
        CHECK_FALSE(f.has(Feature::SignalType_NEGATIVE));
        CHECK_FALSE(f.has(Feature::SignalType_CATASTROPHE));
        CHECK_FALSE(f.has(Feature::SignalBatchRemoval));
    }
    // REPLY is a layered construct and is declared by neither solver.
    CHECK_FALSE(line::qn::ag_feature_set("inap").has(Feature::SignalType_REPLY));
    // Every AG method shares one envelope: they differ in how the reversed rate
    // is read off an agent, not in what a model may contain.
    CHECK(line::qn::ag_feature_set("inap").list()
          == line::qn::ag_feature_set("inapinf").list());
}

TEST_CASE("SolverMAM refuses a G-network, and names the RCAT methods that moved") {
    line::qn::NetworkStruct<double> sn;
    sn.issignal = std::vector<bool>{false, true};
    // No MAM algorithm reads `issignal`, so every MAM method refuses outright
    // rather than answering with the signals turned into ordinary customers.
    for (const char* m : {"default", "dec.source", "dec.poisson", "mna", "ldqbd", "bgchain"}) {
        CHECK_THROWS_AS(line::mam::runner_detail::check_model_method(sn, std::string(m)),
                        UnsupportedError);
    }
    // The RCAT names are refused too, but by a DIFFERENT gate and for a
    // different reason: they moved to SolverAG, and a caller carrying an old
    // options.method is told where they went rather than told "unknown method".
    for (const char* m : {"inap", "inapplus", "inapinf", "exact"}) {
        CHECK_THROWS_AS(line::mam::runner_detail::check_model_method(sn, std::string(m)),
                        UnsupportedError);
    }
    // A model with no signal class is untouched by the signal gate.
    line::qn::NetworkStruct<double> plain;
    plain.issignal = std::vector<bool>{false, false};
    for (const char* m : {"default", "dec.source"})
        CHECK_NOTHROW(line::mam::runner_detail::check_model_method(plain, std::string(m)));
    // ... but the moved names are still redirected on it.
    CHECK_THROWS_AS(line::mam::runner_detail::check_model_method(plain, std::string("inap")),
                    UnsupportedError);
}

TEST_CASE("mam_feature_set applies the mna round-robin delta") {
    // SolverMAM.m:60: solver_mna_open resolves round-robin dispatching as a
    // deterministic traffic split (npfqn_traffic_split_rr); no other MAM
    // analyzer does, and the closed branch is refused in check_model_method.
    CHECK(line::qn::mam_feature_set("mna").has(Feature::RoutingStrategy_RROBIN));
    for (const char* m : {"default", "ldqbd", "dec.source", "dec.poisson"})
        CHECK_FALSE(line::qn::mam_feature_set(m).has(Feature::RoutingStrategy_RROBIN));
}

TEST_CASE("env_feature_set is SolverENV.supports's inline set") {
    const FeatureSet f = line::qn::env_feature_set("meanfield");
    // 2026-07-31: 28 -> 27, RoutingStrategy_RROBIN withdrawn with the CTMC
    // stage solver's own declaration (env_dispatch.h:164 routes to it).
    CHECK(f.list().size() == 27);
    CHECK(f.has(Feature::ClassSwitch));
    CHECK(f.has(Feature::Cox2));
    CHECK(f.has(Feature::Coxian));
    CHECK_FALSE(f.has(Feature::RoutingStrategy_RROBIN));
    CHECK(f.has(Feature::SchedStrategy_FCFS));
    CHECK(f.has(Feature::ClosedClass));
    CHECK(f.has(Feature::OpenClass));
    // The narrowest declared set in the codebase: no APH, no Det, no fork-join.
    CHECK_FALSE(f.has(Feature::APH));
    CHECK_FALSE(f.has(Feature::Det));
    CHECK_FALSE(f.has(Feature::MAP));
    CHECK_FALSE(f.has(Feature::Replayer));
    CHECK_FALSE(f.has(Feature::Cache));
    CHECK_FALSE(f.has(Feature::Fork));
    CHECK_FALSE(f.has(Feature::Region));
    CHECK_FALSE(f.has(Feature::Router));
    CHECK_FALSE(f.has(Feature::SelfLoopingClass));
    CHECK_FALSE(f.has(Feature::SchedStrategy_LCFS));
    CHECK_FALSE(f.has(Feature::LoadDependence));
    CHECK(line::qn::env_feature_set("statevec").list() == f.list());
}

TEST_CASE("every declared set is a subset of the registry and none is empty") {
    const std::vector<FeatureSet> all = {
        line::qn::mva_feature_set("default"),   line::qn::nc_feature_set("default"),
        line::qn::ctmc_feature_set("default"),  line::qn::ssa_feature_set("default"),
        line::qn::fluid_feature_set("default"), line::qn::mam_feature_set("default"),
        line::qn::env_feature_set("default")};
    for (const FeatureSet& f : all) {
        CHECK_FALSE(f.empty());
        CHECK(f.list().size() < feature_count());
        // The LQN-only and placeholder names belong to no NetworkStruct solver.
        CHECK_FALSE(f.has(Feature::Host));
        CHECK_FALSE(f.has(Feature::Task));
        CHECK_FALSE(f.has(Feature::Activity));
        CHECK_FALSE(f.has(Feature::SchedStrategy_REF));
        CHECK_FALSE(f.has(Feature::Disabled));
        CHECK_FALSE(f.has(Feature::Immediate));
        CHECK_FALSE(f.has(Feature::Logger));
        CHECK_FALSE(f.has(Feature::LogTunnel));
        CHECK_FALSE(f.has(Feature::QueueingPlace));
        CHECK_FALSE(f.has(Feature::BatchArrival));
        CHECK_FALSE(f.has(Feature::Trace));
        CHECK_FALSE(f.has(Feature::Normal));
        CHECK_FALSE(f.has(Feature::Poisson));
        CHECK_FALSE(f.has(Feature::Zipf));
        CHECK_FALSE(f.has(Feature::SchedStrategy_EDF));
        // SchedStrategy_FCFSPR left this list on 2026-07-29, see below. FCFSPRIO
        // stays: it is HOL's alias, and every set that serves HOL declares it
        // under the HOL name.
        CHECK_FALSE(f.has(Feature::SchedStrategy_FCFSPRIO));
    }
    // Retrial is declared by CTMC, SSA and MAM only, which is the tightening
    // recorded at feature_set.h:580-585.
    CHECK(all[2].has(Feature::Retrial));
    CHECK(all[3].has(Feature::Retrial));
    CHECK(all[5].has(Feature::Retrial));
    CHECK_FALSE(all[0].has(Feature::Retrial));
    CHECK_FALSE(all[1].has(Feature::Retrial));
    CHECK_FALSE(all[4].has(Feature::Retrial));
    CHECK_FALSE(all[6].has(Feature::Retrial));

    // 2026-07-29: the CTMC generator carries a buffer_is_tag_phase_pairs arm for
    // ALL EIGHT members of the preempt family, in both after_event_station_arv
    // and after_event_station_dep, mirroring afterEventStation.m:1043-1330.
    // Declaring only LCFSPR, LCFSPRPRIO and FCFSPRPRIO gated off five chains the
    // generator builds correctly, so the declaration was corrected to match the
    // generator in both codebases (SolverCTMC.m:159 and ctmc_feature_set).
    // FCFSPR therefore left the never-declared list above. What that assertion
    // was worth is kept here and sharpened: the five are CTMC's ALONE, since no
    // other solver has a rate law for them.
    const Feature preempt[] = {
        Feature::SchedStrategy_FCFSPR, Feature::SchedStrategy_LCFSPI,
        Feature::SchedStrategy_FCFSPI, Feature::SchedStrategy_LCFSPIPRIO,
        Feature::SchedStrategy_FCFSPIPRIO};
    for (std::size_t p = 0; p < 5; ++p) {
        CHECK(all[2].has(preempt[p]));
        for (std::size_t i = 0; i < all.size(); ++i)
            if (i != 2) CHECK_FALSE(all[i].has(preempt[p]));
    }
}

TEST_CASE("a declared set gates a real model through feature_gate") {
    const SN sn = baseline();
    CHECK_NOTHROW(feature_gate("SolverMVA", line::qn::mva_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverNC", line::qn::nc_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverCTMC", line::qn::ctmc_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverSSA", line::qn::ssa_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverFluid", line::qn::fluid_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverMAM", line::qn::mam_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverENV", line::qn::env_feature_set("default"), sn));
    // A closed class is exactly what the rqna and qna deltas withdraw.
    CHECK_THROWS_AS(feature_gate("SolverMVA", line::qn::mva_feature_set("rqna"), sn),
                    UnsupportedError);
    CHECK_THROWS_AS(feature_gate("SolverMVA", line::qn::mva_feature_set("qna"), sn),
                    UnsupportedError);
}

TEST_CASE("a finite capacity region is refused by every solver but NC, CTMC and SSA") {
    SN sn = baseline();
    SN::Region reg;
    reg.members = {true};
    reg.cap = {{-1.0, 4.0}};
    sn.regions.push_back(reg);
    CHECK_NOTHROW(feature_gate("SolverNC", line::qn::nc_feature_set("default"), sn));
    CHECK_NOTHROW(feature_gate("SolverCTMC", line::qn::ctmc_feature_set("default"), sn));
    // SSA joined them on 2026-07-31: the serial engine censors the path under
    // DROP and runs the shared token-FIFO relation under WAITQ, so the region is
    // honoured rather than ignored. The reference declares `Region` for SolverSSA
    // too, on the strength of its own NRM.
    CHECK_NOTHROW(feature_gate("SolverSSA", line::qn::ssa_feature_set("default"), sn));
    CHECK_THROWS_AS(feature_gate("SolverMVA", line::qn::mva_feature_set("default"), sn),
                    UnsupportedError);
    CHECK_THROWS_AS(feature_gate("SolverFluid", line::qn::fluid_feature_set("default"), sn),
                    UnsupportedError);
    CHECK_THROWS_AS(feature_gate("SolverMAM", line::qn::mam_feature_set("default"), sn),
                    UnsupportedError);
    CHECK_THROWS_AS(feature_gate("SolverENV", line::qn::env_feature_set("default"), sn),
                    UnsupportedError);
}


// ---------------------------------------------------------------------------
// The gate as WIRED into the six runner sites
//
// These models are built through `qn::Network` and taken as `get_struct()`, so
// they are REFRESHED structs, not the hand-built ones the derivation cases
// above use. That distinction is load bearing: a hand-built struct carries no
// chains, no visits and no routing table, and a runner that indexes them faults
// long before any conclusion about the runner could be drawn. Only a refreshed
// struct tests what these cases claim to test.
// ---------------------------------------------------------------------------

namespace {

using NetD = line::qn::Network<double>;

/** Think -> Q -> Think, one closed class: the smallest model every solver takes. */
NetD closed_dq() {
    NetD m("cqn2");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(2.0));
    m.set_service(q, c, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * `Network.cluster`: Source -> Dispatcher (Router) -> two PS servers -> Sink,
 * with RAND dispatching. The shape MNetwork.m:970 builds and the one the
 * Router emission used to refuse.
 */
NetD cluster_open_2() {
    NetD m("Cluster");
    const std::size_t src = m.add_source("Source");
    const std::size_t disp = m.add_router("Dispatcher");
    const std::size_t s1 = m.add_queue("Station1", SchedStrategy::PS);
    const std::size_t s2 = m.add_queue("Station2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1", 0);
    m.set_arrival(src, c, Dist::exp_mean(1.0 / 0.4));
    m.set_service(s1, c, Dist::exp_mean(1.0));
    m.set_service(s2, c, Dist::exp_mean(1.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c, c, src, disp, 1.0);
    P.set(c, c, disp, s1, 1.0);
    P.set(c, c, disp, s2, 1.0);
    P.set(c, c, s1, snk, 1.0);
    P.set(c, c, s2, snk, 1.0);
    m.link(P);
    m.set_routing(disp, c, RoutingStrategy::RAND);
    return m;
}

/** The same model under a finite capacity region over both stations. */
NetD closed_dq_region() {
    NetD m = closed_dq();
    m.add_region(std::vector<std::size_t>{1, 2}, std::vector<double>{-1.0}, 2.0,
                 std::vector<DropStrategy>{DropStrategy::DROP});
    return m;
}

}  // namespace

/**
 * A dispatcher model is ADMITTED by SolverMVA and SolverNC.
 *
 * The regression for the Router emission: `used_lang_features` used to emit
 * Feature::Router, which neither featset declares, so every
 * Network.cluster* model was refused with "this model uses a Router node,
 * which this solver does not support" while MATLAB and python solved it. The
 * Router carries no jobs and is eliminated by the stochastic complement, so the
 * two servers split the arrivals evenly: lambda 0.4 over two servers of mean
 * service 1.0 gives utilization 0.2 at each.
 */
TEST_CASE("a Router model is admitted by SolverMVA and SolverNC") {
    NetD m = cluster_open_2();
    const SN sn = m.get_struct();

    const FeatureSet u = used_lang_features(sn);
    CHECK_FALSE(u.has(Feature::Router));
    CHECK(u.has(Feature::RoutingStrategy_RAND));

    line::mva::MvaOptions mopt;
    Matrix<double> init;
    const line::mva::AvgResult<double> rm = line::mva::solver_mva_run_analyzer(sn, mopt, init);
    CHECK(rm.UN(1, 0) == doctest::Approx(0.2).epsilon(1e-9));
    CHECK(rm.UN(2, 0) == doctest::Approx(0.2).epsilon(1e-9));
    CHECK(rm.TN(1, 0) == doctest::Approx(0.2).epsilon(1e-9));

    line::nc::NcSolverOptions nopt;
    const line::mva::AvgResult<double> rn = line::nc::solver_nc_run_analyzer(sn, nopt);
    CHECK(rn.UN(1, 0) == doctest::Approx(0.2).epsilon(1e-9));
    CHECK(rn.UN(2, 0) == doctest::Approx(0.2).epsilon(1e-9));
}

TEST_CASE("mva resolve_method reproduces SolverMVA.resolveMethod") {
    const SN closed = closed_dq().get_struct();
    // A closed model never upgrades: rqna is an open-network analyzer.
    CHECK(line::mva::resolve_method(closed, "default") == "default");
    // A named method is never rewritten, whatever the model.
    CHECK(line::mva::resolve_method(closed, "exact") == "exact");
    CHECK(line::mva::resolve_method(closed, "qna") == "qna");
    CHECK(line::mva::resolve_method(closed, "rqna") == "rqna");
}

TEST_CASE("the MVA runner refuses a finite capacity region by name") {
    const SN sn = closed_dq_region().get_struct();
    line::mva::MvaOptions opt;
    Matrix<double> init;
    CHECK_THROWS_AS(line::mva::solver_mva_run_analyzer(sn, opt, init), UnsupportedError);
    try {
        line::mva::solver_mva_run_analyzer(sn, opt, init);
        FAIL("the MVA runner did not refuse a finite capacity region");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("Finite Capacity Region") != std::string::npos);
        CHECK(std::string(e.what()).find("SolverMVA") != std::string::npos);
    }
    // The same model without the region is solved, so the refusal is the
    // region's doing and not the model's.
    CHECK_NOTHROW(line::mva::solver_mva_run_analyzer(closed_dq().get_struct(), opt, init));
}

TEST_CASE("the fluid runner refuses a finite capacity region by name") {
    // A Finite Capacity Region is NO LONGER a blanket fluid refusal. `dae`
    // carries the cap as an algebraic equation beside the drift, so
    // `fluid_feature_set("dae")` declares Feature::Region and
    // `fluid_resolve_method` sends `default` down that route on a region model.
    // The refusal therefore splits in three, and all three are the reference's
    // (verified against MATLAB R2026a on the same model, 2026-08-22).
    const SN sn = closed_dq_region().get_struct();
    line::fluid::FluidOptions opt;

    // 1. A DROP rule is not a constraint on this drift -- it changes the EVENT
    //    SET rather than throttling the admission flow -- so `default` reaches
    //    the dae route and is refused there, naming the region AND the rule.
    //    MATLAB: "Region 1 applies drop to class 1. Only a waiting queue is a
    //    constraint on the fluid drift...". Naming the real blocker is the
    //    point: the old blanket message said only "Finite Capacity Region",
    //    which does not tell the user that it is the DROP rule that has to go.
    CHECK_THROWS_AS(line::fluid::solver_fluid_run_analyzer(sn, opt), UnsupportedError);
    try {
        line::fluid::solver_fluid_run_analyzer(sn, opt);
        FAIL("the fluid runner did not refuse a DROP region");
    } catch (const UnsupportedError& e) {
        const std::string msg = e.what();
        CHECK(msg.find("region 1") != std::string::npos);
        CHECK(msg.find("drop rule other than a waiting queue") != std::string::npos);
    }

    // 2. Every method that is NOT dae still refuses ANY region through the
    //    universal feature gate, by the registry name: an ODE integrated
    //    through a cap it cannot see returns the unconstrained answer silently.
    //    This is the assertion this test has always made, on the methods where
    //    it still holds.
    line::fluid::FluidOptions ode;
    ode.method = "closing";
    try {
        line::fluid::solver_fluid_run_analyzer(sn, ode);
        FAIL("an ODE fluid method did not refuse a finite capacity region");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("Finite Capacity Region") != std::string::npos);
    }

    // 3. A WAITQ region IS expressible, and `default` now solves it rather than
    //    refusing it -- the capability the split exists for. MATLAB reports
    //    "method: default/dae" on this model and returns a table.
    NetD waitq = closed_dq();
    waitq.add_region(std::vector<std::size_t>{1, 2}, std::vector<double>{-1.0}, 2.0,
                     std::vector<DropStrategy>{DropStrategy::WAITQ});
    CHECK_NOTHROW(line::fluid::solver_fluid_run_analyzer(waitq.get_struct(), opt));

    // and the same model without any region is solved, so none of the above is
    // the model's own doing.
    CHECK_NOTHROW(line::fluid::solver_fluid_run_analyzer(closed_dq().get_struct(), opt));
}

TEST_CASE("the MAM runner refuses a finite capacity region by name") {
    const SN sn = closed_dq_region().get_struct();
    line::mam::MamOptions opt;
    try {
        line::mam::solver_mam_solve(sn, opt);
        FAIL("the MAM runner did not refuse a finite capacity region");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("Finite Capacity Region") != std::string::npos);
    }
}

TEST_CASE("NC declares Region, so the gate defers to its imperative split") {
    // SolverNC.m:185-190: NC solves the OPEN single-Delay loss network exactly,
    // so the gate must NOT pre-empt solver_nc_runner.h:315-332. The refusal here
    // must therefore name the queueing-station shape, never the gate's phrase.
    const SN sn = closed_dq_region().get_struct();
    line::nc::NcSolverOptions opt;
    try {
        line::nc::solver_nc_solve(sn, opt);
        FAIL("SolverNC did not refuse a region on queueing stations");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.find("Finite Capacity Region") != std::string::npos);
        // the imperative message, not the gate's
        CHECK(w.find("which this solver does not support") == std::string::npos);
        CHECK(w.find("queueing stations") != std::string::npos);
    }
}

TEST_CASE("CTMC declares Region, so the gate lets a region model through") {
    // CTMC represents the region exactly; the gate must be inert here. Compared
    // against the same model without the region, which must also solve.
    const line::ctmc::CtmcOptions opt;
    CHECK_NOTHROW(line::ctmc::solver_ctmc_analyzer(closed_dq().get_struct(), opt));
    CHECK_NOTHROW(line::ctmc::solver_ctmc_analyzer(closed_dq_region().get_struct(), opt));
    // and through the WAITQ-aware entry, which carries its own gate
    CHECK_NOTHROW(line::ctmc::solver_ctmc_analyzer_any(closed_dq_region().get_struct(), opt));
}

TEST_CASE("feature_gate emits both NetworkSolver.m:184-190 message forms") {
    const SN sn = closed_dq().get_struct();
    FeatureSet declared;  // declares nothing, so any model is refused
    // Form 1: the user named the method, so the solver alone is blamed.
    try {
        feature_gate("SolverMVA", declared, sn, "exact", "exact");
        FAIL("the plain form did not refuse");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.compare(0, 10, "SolverMVA:") == 0);
        CHECK(w.find("method") == std::string::npos);
    }
    // Form 2: resolveMethod changed the method, so the message names it.
    try {
        feature_gate("SolverMVA", declared, sn, "default", "rqna");
        FAIL("the resolved form did not refuse");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).compare(0, 26, "SolverMVA's 'rqna' method:") == 0);
    }
    // Omitting both arguments keeps the plain form, which is what the other
    // five call sites rely on.
    try {
        feature_gate("SolverNC", declared, sn);
        FAIL("the default form did not refuse");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).compare(0, 9, "SolverNC:") == 0);
    }
}

TEST_CASE("the method whitelist still precedes the feature gate at MVA") {
    // 'rqna' is advertised only for a fully open model (list_valid_methods), so
    // a closed model asking for it is refused by the METHOD gate first. This
    // pins the ordering: the feature gate is additive and never pre-empts a
    // more specific refusal that already exists.
    const SN sn = closed_dq().get_struct();
    line::mva::MvaOptions opt;
    opt.method = "rqna";
    Matrix<double> init;
    try {
        line::mva::solver_mva_run_analyzer(sn, opt, init);
        FAIL("SolverMVA did not refuse rqna on a closed model");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.find("method is unsupported") != std::string::npos);
        CHECK(w.find("does not support") == std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// The finite station/class capacity gate (BUG-39)
//
// SolverMVA.m:143-153 lays SolverMVA.supportsFiniteCapacity (SolverMVA.m:158-186)
// on top of the universal feature gate, which cannot see capacity: there is no
// LINE_QN_FEATURE_LIST enumerator for it, because it is a NUMBER on a station
// and not a construct. Without the check a model built with setCapacity is
// solved as if uncapacitated and returns a WRONG NUMBER, not a refusal.
//
// Recorded 2026-07-29. The reference draws three lines these cases pin:
//   1. only a capacity that can BIND is refused (NetworkSolver.m:1189-1194);
//   2. a BAS model is exempt, because solver_mva_analyzer sends it to SQD,
//      which does honour the buffers -- but the exemption predicate is the
//      single-CLASS sn_is_bas_model, not the single-CHAIN dispatch test;
//   3. a single-station M/M/1/K with tail drop is exempt under every method
//      BUT 'exact', because the mg1k.mgs branch that solves it is an
//      approximation away from scv=1.
// ---------------------------------------------------------------------------

namespace {

/** closed_dq with a station capacity on the queue. */
NetD closed_dq_cap(double cap) {
    NetD m("cqn2cap");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(2.0));
    m.set_service(q, c, Dist::exp_rate(3.0));
    m.set_capacity(q, cap);
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** closed_dq with a per-class buffer on the queue. */
NetD closed_dq_classcap(double ccap) {
    NetD m("cqn2ccap");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(2.0));
    m.set_service(q, c, Dist::exp_rate(3.0));
    m.set_class_capacity(q, c, ccap);
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue(K, tail drop) -> Sink, exponential throughout: the M/M/1/K. */
NetD mm1k_loss(double lambda, double mu, double K) {
    NetD m("mm1k");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, K);
    m.set_drop_rule(q, c, DropStrategy::DROP);
    line::qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/**
 * Think -> B1(cap 2, BAS) -> B2(cap 2, BAS) -> Think, `nclasses` closed classes
 * of 2 jobs each. With two classes the return leg SWITCHES class, so the two sit
 * in ONE chain: that is what separates the single-chain dispatch predicate from
 * the single-class gate predicate.
 */
NetD bas_cycle(std::size_t nclasses) {
    NetD m("bascycle");
    const std::size_t d = m.add_delay("Think");
    const std::size_t b1 = m.add_queue("B1", SchedStrategy::FCFS);
    const std::size_t b2 = m.add_queue("B2", SchedStrategy::FCFS);
    m.set_capacity(b1, 2.0);
    m.set_capacity(b2, 2.0);
    std::vector<std::size_t> cls;
    for (std::size_t j = 0; j < nclasses; ++j) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(j + 1), 2.0, d);
        m.set_service(d, c, Dist::exp_rate(2.0));
        m.set_service(b1, c, Dist::exp_rate(1.5));
        m.set_service(b2, c, Dist::exp_rate(2.5));
        m.set_drop_rule(b1, c, DropStrategy::BAS);
        m.set_drop_rule(b2, c, DropStrategy::BAS);
        cls.push_back(c);
    }
    line::qn::RoutingMatrix<double> P;
    for (std::size_t j = 0; j < nclasses; ++j) {
        const std::size_t c = cls[j];
        const std::size_t nxt = cls[(j + 1) % nclasses];
        P.set(c, c, d, b1, 1.0);
        P.set(c, c, b1, b2, 1.0);
        P.set(c, nxt, b2, d, 1.0);  // the class switch that merges the chains
    }
    m.link(P);
    return m;
}

line::mva::AvgResult<double> mva_of(const SN& sn, const std::string& method = "default") {
    line::mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return line::mva::solver_mva_run_analyzer(sn, opt, init);
}

}  // namespace

TEST_CASE("the MVA runner refuses a binding station capacity by name") {
    // cap = 1 < 2 jobs, so a job CAN be turned away: the product-form answer is
    // not the answer to this model. NetworkSolver.m:1189-1194.
    try {
        mva_of(closed_dq_cap(1.0).get_struct());
        FAIL("SolverMVA did not refuse a binding station capacity");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.compare(0, 10, "SolverMVA:") == 0);
        CHECK(w.find("setCapacity=1") != std::string::npos);
        CHECK(w.find("'Q'") != std::string::npos);
        CHECK(w.find("SolverCTMC") != std::string::npos);
    }
}

TEST_CASE("a station capacity that cannot bind is not a refusal") {
    // setCapacity(N) on a station of an N-job closed model is a common idiom and
    // a no-op: no job is ever turned away, so the product-form answer is exact.
    // Pinned against the same model with no capacity at all.
    const line::mva::AvgResult<double> base = mva_of(closed_dq().get_struct());
    for (double cap : {2.0, 5.0}) {
        const line::mva::AvgResult<double> r = mva_of(closed_dq_cap(cap).get_struct());
        CHECK(r.QN(1, 0) == doctest::Approx(base.QN(1, 0)).epsilon(1e-12));
        CHECK(r.RN(1, 0) == doctest::Approx(base.RN(1, 0)).epsilon(1e-12));
        CHECK(r.TN(1, 0) == doctest::Approx(base.TN(1, 0)).epsilon(1e-12));
    }
}

TEST_CASE("the MVA runner refuses a binding per-class capacity by name") {
    // The per-class test is against that class's OWN population, not the total.
    try {
        mva_of(closed_dq_classcap(1.0).get_struct());
        FAIL("SolverMVA did not refuse a binding per-class capacity");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.compare(0, 10, "SolverMVA:") == 0);
        CHECK(w.find("classCap=1 for class 1") != std::string::npos);
        CHECK(w.find("'Q'") != std::string::npos);
    }
    // classCap = the class population cannot bind either.
    CHECK_NOTHROW(mva_of(closed_dq_classcap(2.0).get_struct()));
}

TEST_CASE("an open class makes any finite capacity binding") {
    // njobs is Inf for an open class, so cap < totalJobs holds for every finite
    // cap. This model is NOT the M/M/1/K shape (two queues), so no exemption
    // applies and the refusal must fire.
    NetD m("openpair");
    const std::size_t s = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, Dist::exp_rate(0.5));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    m.set_capacity(q2, 100.0);
    m.set_drop_rule(q2, c, DropStrategy::DROP);
    line::qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    try {
        mva_of(m.get_struct());
        FAIL("SolverMVA did not refuse a finite capacity reachable by an open class");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("setCapacity=100") != std::string::npos);
    }
}

TEST_CASE("the capacity gate exempts a single-class BAS model, not a multiclass one") {
    // SQD honours the buffers, so a BAS model is exempt -- but the exemption
    // predicate is sn_is_bas_model (single CLASS), narrower than the analyzer's
    // dispatch test mva_is_bas_model (single CHAIN). The Smith decomposition
    // models ONE circulating population, so a chain built from two classes is
    // dispatched to SQD yet is NOT exempt from the refusal.
    const SN one = bas_cycle(1).get_struct();
    CHECK(line::mva::sn_is_bas_model(one));
    CHECK(line::mva::mva_is_bas_model(one));
    const line::mva::AvgResult<double> r = mva_of(one);
    CHECK(r.actualmethod == "sqd");

    const SN two = bas_cycle(2).get_struct();
    CHECK_FALSE(line::mva::sn_is_bas_model(two));  // two classes
    CHECK(line::mva::mva_is_bas_model(two));       // but one chain
    try {
        mva_of(two);
        FAIL("SolverMVA did not refuse a multiclass BAS model's finite buffers");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find("setCapacity=2") != std::string::npos);
    }
}

TEST_CASE("the capacity gate exempts M/M/1/K under every method but exact") {
    const SN sn = mm1k_loss(1.0, 2.0, 3.0).get_struct();
    CHECK(line::mva::sn_is_mm1k_loss(sn));

    // rho = 1/2, K = 3. The MacGregor Smith loss probability reduces to the
    // exact M/M/1/K one at scv = 1: Ploss = 1/15, so the carried rate is 14/15
    // and the truncated distribution gives Lsys = 11/15, R = Lsys/Tq = 11/14.
    const line::mva::AvgResult<double> r = mva_of(sn);
    CHECK(r.actualmethod == "mg1k.mgs");
    CHECK(r.QN(1, 0) == doctest::Approx(11.0 / 15.0).epsilon(1e-10));
    CHECK(r.UN(1, 0) == doctest::Approx(7.0 / 15.0).epsilon(1e-10));
    CHECK(r.RN(1, 0) == doctest::Approx(11.0 / 14.0).epsilon(1e-10));
    CHECK(r.TN(1, 0) == doctest::Approx(14.0 / 15.0).epsilon(1e-10));
    CHECK(r.TN(0, 0) == doctest::Approx(1.0).epsilon(1e-10));

    // and it is NOT the unconstrained M/M/1, which would hold rho/(1-rho) = 1.
    CHECK(r.QN(1, 0) < 1.0);

    // 'exact' is not exempt: the branch is an approximation away from scv=1, so
    // the reference makes MVA refuse rather than answer. Since 2026-09-05 it is
    // the FEATURE gate that refuses, not the capacity check: getMethodFeatureSet
    // withholds FiniteCapacity from 'exact', and SolverMVA.supportsModelMethod
    // asks the base featset BEFORE supportsFiniteCapacity (SolverMVA.m:408-411),
    // so the specific setCapacity sentence is unreachable under this name in the
    // reference too. What must hold is that it refuses and says the buffer is why.
    try {
        mva_of(sn, "exact");
        FAIL("SolverMVA solved an M/M/1/K under method='exact'");
    } catch (const UnsupportedError& e) {
        const std::string msg(e.what());
        CHECK_MESSAGE(msg.find("buffer") != std::string::npos, msg);
    }
}

TEST_CASE("the feature gate still precedes the capacity check") {
    // A model that fails BOTH must be refused by the FEATURE, whose message
    // names the construct to remove. supportsModelMethod tests capacity only
    // `if bool` (SolverMVA.m:150-152), so the order is load bearing.
    NetD m = closed_dq_cap(1.0);
    m.add_region(std::vector<std::size_t>{1, 2}, std::vector<double>{-1.0}, 2.0,
                 std::vector<DropStrategy>{DropStrategy::DROP});
    try {
        mva_of(m.get_struct());
        FAIL("SolverMVA did not refuse a region plus a binding capacity");
    } catch (const UnsupportedError& e) {
        const std::string w = e.what();
        CHECK(w.find("Finite Capacity Region") != std::string::npos);
        // The region message names setCapacity as the single-station REMEDY
        // (runAnalyzer.m:24 does too), so the bare word cannot separate the two
        // gates; check_binding_capacity's own wording can.
        CHECK(w.find("finite station capacity (setCapacity=") == std::string::npos);
    }
}

TEST_CASE("check_binding_capacity exempts a Cache model wholesale") {
    // Cache.m sets classCap=1 on the retrieval queues it builds, and the cache
    // analyzers solve those rather than treating them as buffers, so one Cache
    // node exempts the whole model (NetworkSolver.m:1205-1207).
    SN sn = baseline();
    sn.stations[0].classcap.assign(1, 1.0);  // 1 < the 2-job population: binding
    CHECK_THROWS_AS(line::qn::check_binding_capacity("SolverMVA", sn), UnsupportedError);
    sn.add_node("C", NodeType::Cache, true);
    CHECK_NOTHROW(line::qn::check_binding_capacity("SolverMVA", sn));
}

TEST_CASE("check_binding_capacity ignores a Source and a zero class buffer") {
    // A Source holds no jobs, so its capacity can never bind; and classCap = 0
    // marks a class the station never serves, which the reference skips with
    // its `ccap(r) > 0` test rather than reporting as a zero-size buffer.
    SN sn = baseline();
    const std::size_t s = add_st(sn, "Source", NodeType::Source, SchedStrategy::EXT);
    sn.nstations = sn.stations.size();
    sn.stations[s - 1].cap = 0.0;
    CHECK_NOTHROW(line::qn::check_binding_capacity("SolverMVA", sn));
    sn.stations[0].classcap.assign(1, 0.0);
    CHECK_NOTHROW(line::qn::check_binding_capacity("SolverMVA", sn));
}

// ---------------------------------------------------------------------------
// The reference's warning text, carried out through the runner
//
// mva_dispatch already pins that DispatchResult::warning is set by the SJN
// starvation cap and empty otherwise. These two cases pin the LAST hop, which
// the dispatch-level ones cannot see: a caller who uses solver_mva_run_analyzer rather
// than mva_dispatch must get the same text. Dropped there, a capped answer
// reaches the user looking authoritative exactly where MATLAB declines to
// stand behind it.
// ---------------------------------------------------------------------------

namespace {

/** Delay -> SJF queue -> Delay, one closed class: the SJN model. */
NetD sjn_closed(const std::string& name, double njobs) {
    NetD m(name);
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the MVA runner carries the reference warning verbatim") {
    // N = 6 at rho -> 1 drives the SJN station into the starvation regime and
    // the utilization cap binds. The answer is USABLE -- the population law
    // holds exactly below -- so the text is carried, not thrown.
    const line::mva::AvgResult<double> r = mva_of(sjn_closed("sjnwarn", 6.0).get_struct(), "amva");
    CHECK(r.actualmethod == "sjn.amva");
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(6.0).epsilon(1e-12));
    CHECK_FALSE(r.warning.empty());
    CHECK(r.warning.find("starvation regime") != std::string::npos);
    CHECK(r.warning.find("SolverLDES") != std::string::npos);
}

TEST_CASE("the MVA runner leaves the warning empty when the reference is silent") {
    // Set by the CONDITION, not by the branch: an SJN solve that does not hit
    // the cap must be silent, or the field becomes noise and stops informing.
    const line::mva::AvgResult<double> quiet = mva_of(sjn_closed("sjnquiet", 3.0).get_struct());
    CHECK(quiet.actualmethod == "sjn.mva");
    CHECK(quiet.warning.empty());
    // and so must an ordinary product-form solve that never reaches SJN at all
    CHECK(mva_of(closed_dq().get_struct()).warning.empty());
}

// ---------------------------------------------------------------------------
// The last two names this port carried alone (aligned 2026-08-22)
// ---------------------------------------------------------------------------

TEST_CASE("heterogeneous server pools are marked and gated by their own name") {
    // A pooled station is NOT a station of the same total size: the pools carry
    // their own class compatibilities and their own rates. Every solver that
    // reads only `nservers` therefore has to refuse rather than flatten them.
    SN sn = baseline();
    typename line::qn::Station<double>::ServerType fast;
    fast.name = "fast";
    fast.count = 1.0;
    typename line::qn::Station<double>::ServerType slow;
    slow.name = "slow";
    slow.count = 2.0;
    sn.stations[0].server_types.push_back(fast);
    sn.stations[0].server_types.push_back(slow);

    const FeatureSet used = line::qn::used_lang_features(sn);
    CHECK(used.has(Feature::HeteroServers));
    // and it is absent from the same model without the pools
    CHECK_FALSE(line::qn::used_lang_features(baseline()).has(Feature::HeteroServers));

    // JMT serialises the pools (serverTypesNames / serverTypesNumOfServers /
    // serverTypesCompatibilities), so it accepts; the analytical solvers do not.
    CHECK(line::qn::jmt_feature_set("default").has(Feature::HeteroServers));
    CHECK_FALSE(line::qn::mva_feature_set("default").has(Feature::HeteroServers));
    CHECK_FALSE(line::qn::nc_feature_set("default").has(Feature::HeteroServers));
    CHECK_FALSE(line::qn::ctmc_feature_set("default").has(Feature::HeteroServers));
    // 2026-09-04: THE NATIVE LDES ENGINE NOW SIMULATES THE POOLS, so it declares
    // the name like the other three. Each job occupies one server of one
    // compatible pool and is served at that pool's own law -- `free_slot_for`
    // picks the pool by the station's heterogeneous scheduling policy and
    // `buffer_pop_for_slot` takes the first waiting job in service order that
    // the slot's pool accepts. What the engine still refuses BY NAME is the
    // pools under a sharing discipline, which no engine allocates per pool.
    CHECK(line::qn::ldes_feature_set("default").has(Feature::HeteroServers));

    CHECK(feature_set_supports("SolverJMT", line::qn::jmt_feature_set("default"), used).ok);
    const SupportResult r =
        feature_set_supports("SolverMVA", line::qn::mva_feature_set("default"), used);
    CHECK_FALSE(r.ok);
    CHECK(r.reason.find("HeteroServers") != std::string::npos);
}

TEST_CASE("a non-NORMAL depository departure discipline is refused by every solver") {
    // A FIFO depository releases a served token only after the tokens that
    // entered service before it, so which output transitions are enabled
    // depends on the arrival order and not only on the marking. NO engine in
    // ANY codebase implements it, so no feature set declares it and the model
    // is refused rather than served as if it were NORMAL.
    SN sn = baseline();
    sn.stations[0].departure_discipline.assign(sn.classes.size(),
                                               line::lang::DepartureDiscipline::NORMAL);
    CHECK_FALSE(line::qn::used_lang_features(sn).has(Feature::DepartureDiscipline));

    sn.stations[0].departure_discipline[0] = line::lang::DepartureDiscipline::FIFO;
    const FeatureSet used = line::qn::used_lang_features(sn);
    CHECK(used.has(Feature::DepartureDiscipline));

    for (const FeatureSet& f : {line::qn::jmt_feature_set("default"),
                                line::qn::ldes_feature_set("default"),
                                line::qn::ctmc_feature_set("default"),
                                line::qn::ssa_feature_set("default"),
                                line::qn::mva_feature_set("default")}) {
        CHECK_FALSE(f.has(Feature::DepartureDiscipline));
        const SupportResult res = feature_set_supports("s", f, used);
        CHECK_FALSE(res.ok);
        CHECK(res.reason.find("DepartureDiscipline") != std::string::npos);
    }
}
