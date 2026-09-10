/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The per-method support gates of SolverJMT and SolverMAM.
 *
 * WHAT THIS PINS. `model.help()` / `findSolver()` reports one row per (solver,
 * method) pair and asks the per-method gate -- the same one the AUTO ranking
 * consults before delegating. Until the rules below reached it, that gate was
 * much weaker than what the ANALYZERS enforce at run time, so the report
 * offered pairs that then threw or, worse, answered with a table of zeros:
 *
 *   jmt jmva.*            the JSIM envelope was declared for the JMVA engine,
 *                         which reads a station type, a per-chain demand, a
 *                         per-chain visit count, the populations and a
 *                         reference station and NOTHING ELSE. On a cache model
 *                         all eight closed-form jmva methods returned an
 *                         entirely zero table, jmva.mva labelled 'exact'.
 *   jmt jmva.<closed>     open classes, load-dependent and multi-server
 *                         stations are each refused by JMT itself
 *   jmt replication       a finite horizon
 *   mam dec.mmap          open models only, and no INF station
 *   mam retrial           a retrial topology must be PRESENT
 *
 * TWO ROUTES REACH THE REPORT, and both are asserted here. What the feature
 * registry can NAME rides in `qn::jmt_feature_set(method)` and
 * `qn::mam_feature_set(method)`; what it cannot -- a server count, a horizon,
 * a missing orbit -- is a structural predicate the ANALYZER ITSELF calls
 * (`jmt::jmt_method_refusal`, `mam::mam_model_method_refusal`), so the gate and
 * the run are one body of rules rather than two copies.
 *
 * Each case asserts the refusal AND its converse: a model the method IS derived
 * for must keep it. Over-tightening a gate hides a method the user could have
 * run, which is the same defect with the sign flipped.
 *
 * The MATLAB twins are @SolverJMT/supportsModelMethod.m, jmtMethodRefusal.m and
 * @SolverMAM/supportsModelMethod.m; the Python twin is
 * python/tests/test_gate_jmt_mam.py and the JAR twin
 * jar/src/test/java/jline/solvers/wrappers/jmt/JmtMamGateTest.java.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"

using line::qn::CacheParam;
using line::qn::Feature;
using line::qn::feature_set_supports;
using line::qn::FeatureSet;
using line::qn::NetworkStruct;
using line::qn::NodeType;
using line::qn::used_lang_features;
using line::lang::ReplacementStrategy;
using line::lang::SchedStrategy;

namespace {

using NetD = line::qn::Network<double>;
using Dist = line::lang::Distrib<double>;
using SN = NetworkStruct<double>;

const std::vector<std::string> kJmvaClosedOnly = {
    "jmva.amva", "jmva.recal", "jmva.comom", "jmva.chow",
    "jmva.bs",   "jmva.aql",   "jmva.lin",   "jmva.dmlin"};
const std::vector<std::string> kJmvaExact = {"jmva", "jmva.mva"};

/** Source -> FCFS Queue -> Sink, one open class: open, single server. */
NetD mm1() {
    NetD m("mm1");
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Delay -> FCFS Queue, N = 3: the shape every jmva algorithm serves. */
NetD closed_single_server() {
    NetD m("repairmen");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The same with c = 3, which the closed-form jmva algorithms refuse. */
NetD closed_multiserver() {
    NetD m("multiserver");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    m.set_number_of_servers(q, 3.0);
    const std::size_t c = m.add_closed_class("C", 4.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Two open classes at a head-of-line priority queue. */
NetD open_hol() {
    NetD m("prio");
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::HOL);
    const std::size_t k = m.add_sink("K");
    const std::size_t hi = m.add_open_class("Hi", 0);
    const std::size_t lo = m.add_open_class("Lo", 1);
    m.set_arrival(s, hi, Dist::exp_rate(0.4));
    m.set_arrival(s, lo, Dist::exp_rate(0.4));
    m.set_service(q, hi, Dist::exp_rate(2.0));
    m.set_service(q, lo, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(hi, hi, s, q, 1.0);
    P.set(hi, hi, q, k, 1.0);
    P.set(lo, lo, s, q, 1.0);
    P.set(lo, lo, q, k, 1.0);
    m.link(P);
    return m;
}

/**
 * Delay -> FCFS Queue with cap 2 and N = 4: the buffer BINDS.
 *
 * LINE blocks a closed job that finds no room; no JMT drop strategy reproduces
 * that, and the JMVA document has no capacity element at all.
 */
NetD closed_binding_buffer() {
    NetD m("blk");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    m.set_capacity(q, 2.0);
    const std::size_t c = m.add_closed_class("C", 4.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * M/M/1/2: an OPEN binding buffer, which JSIM DOES simulate.
 *
 * JMT's queue section carries the drop rule directly, so a refused arrival is
 * lost in JMT exactly as it is in LINE. The over-tightening guard.
 */
NetD open_binding_buffer() {
    NetD m("mm1k");
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("K");
    m.set_capacity(q, 2.0);
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/**
 * Two closed classes at one PS queue: sn.cap is DERIVED as 2 x N here.
 *
 * The multi-class shape is the one that makes "finite cap" and "binding cap"
 * different questions, so it is the model the buffer rule must leave alone.
 */
NetD cqn_two_class() {
    NetD m("cqn2");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(2.0));
    m.set_service(q, c1, Dist::exp_rate(2.0));
    m.set_service(q, c2, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Source -> Fork -> two FCFS queues -> Join -> Sink, one open class.
 *
 * `solver_mam_decmmap`'s sweep uses `solver_mam_traffic`, the PLAIN traffic
 * step, so a Fork or a Join is refused by the traffic step itself rather than
 * synchronized; `mam_dispatch` sends this shape to `solver_mam_fj` from
 * 'default' and 'dec.source' and never to dec.mmap.
 */
NetD open_forkjoin() {
    NetD m("fj");
    const std::size_t s = m.add_source("S");
    const std::size_t fk = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("J", fk);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, Dist::exp_rate(0.5));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(s, fk, 1.0);
    P.set(fk, q1, 1.0);
    P.set(fk, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, k, 1.0);
    m.link(P);
    return m;
}

/**
 * A closed network carrying an LRU Cache node.
 *
 * The cache is attached to the finished struct rather than routed through the
 * builder, exactly as test_feature_set.cpp does: what is under test here is
 * whether the DECLARED envelope admits the Cache name, and a full cache network
 * would add routing that the assertion does not depend on.
 */
SN cache_struct() {
    SN sn = closed_single_server().get_struct();
    const std::size_t cn = sn.add_node("Cache", NodeType::Cache, true);
    CacheParam<double> par;
    par.nitems = 4;
    par.itemcap = {2};
    par.replacestrat = ReplacementStrategy::LRU;
    sn.nodeparam[cn] = par;
    return sn;
}

/** Does `declared` accept everything `sn` uses? */
bool accepts(const FeatureSet& declared, const SN& sn) {
    return feature_set_supports("gate", declared, used_lang_features(sn)).ok;
}

}  // namespace

// ---------------------------------------------------------------------------
// SolverJMT: two engines, two envelopes
// ---------------------------------------------------------------------------

TEST_CASE("jmt_feature_set: the JMVA envelope drops what write_jmva cannot emit") {
    const FeatureSet jsim = line::qn::jmt_feature_set("jsim");
    const FeatureSet jmva = line::qn::jmt_feature_set("jmva");

    // Constructs the JMVA document has no element for at all.
    const Feature absent[] = {Feature::Cache,          Feature::CacheClassSwitcher,
                              Feature::ReplacementStrategy_LRU, Feature::Fork,
                              Feature::Join,           Feature::Region,
                              Feature::Place,          Feature::Transition,
                              Feature::Reneging,       Feature::Balking,
                              Feature::HeteroServers,  Feature::SetupDelayOff,
                              Feature::ServerParallelism};
    for (Feature f : absent) {
        CHECK_MESSAGE(jsim.has(f), line::qn::feature_name(f));
        CHECK_MESSAGE(!jmva.has(f), line::qn::feature_name(f));
    }
    // The writer emits NO discipline, so only the BCMP station types survive.
    const Feature no_discipline[] = {Feature::SchedStrategy_HOL, Feature::SchedStrategy_DPS,
                                     Feature::SchedStrategy_GPS, Feature::SchedStrategy_POLLING,
                                     Feature::SchedStrategy_SRPT,
                                     Feature::SchedStrategy_FCFSPRPRIO};
    for (Feature f : no_discipline) CHECK_MESSAGE(!jmva.has(f), line::qn::feature_name(f));
    const Feature bcmp[] = {Feature::SchedStrategy_INF, Feature::SchedStrategy_PS,
                            Feature::SchedStrategy_FCFS, Feature::SchedStrategy_LCFSPR};
    for (Feature f : bcmp) CHECK_MESSAGE(jmva.has(f), line::qn::feature_name(f));
    // Mean visit counts are not a join-the-shortest-queue dispatcher.
    CHECK_FALSE(jmva.has(Feature::RoutingStrategy_JSQ));
    CHECK_FALSE(jmva.has(Feature::RoutingStrategy_RROBIN));
    CHECK(jmva.has(Feature::RoutingStrategy_PROB));
    // A mean demand is all JMVA reads, so the laws stay.
    const Feature laws[] = {Feature::Exp, Feature::Erlang, Feature::HyperExp, Feature::Pareto,
                            Feature::Replayer};
    for (Feature f : laws) CHECK_MESSAGE(jmva.has(f), line::qn::feature_name(f));
}

TEST_CASE("jmt_feature_set: the closed-only algorithms drop two more names") {
    const FeatureSet base = line::qn::jmt_feature_set("jmva");
    CHECK(base.has(Feature::OpenClass));
    CHECK(base.has(Feature::LoadDependence));
    for (const std::string& method : kJmvaExact) {
        CHECK_FALSE(line::qn::jmva_is_closed_only(method));
        CHECK_MESSAGE(line::qn::jmt_feature_set(method).has(Feature::OpenClass), method);
        CHECK_MESSAGE(line::qn::jmt_feature_set(method).has(Feature::LoadDependence), method);
    }
    for (const std::string& method : kJmvaClosedOnly) {
        CHECK_MESSAGE(line::qn::jmva_is_closed_only(method), method);
        const FeatureSet f = line::qn::jmt_feature_set(method);
        CHECK_MESSAGE(!f.has(Feature::OpenClass), method);
        CHECK_MESSAGE(!f.has(Feature::LoadDependence), method);
        CHECK_MESSAGE(f.has(Feature::ClosedClass), method);
        CHECK_MESSAGE(f.has(Feature::Queue), method);
    }
    // The simulation methods keep the whole JSIM envelope.
    const std::string sim[] = {"default", "jsim", "replication"};
    for (const std::string& method : sim)
        CHECK_MESSAGE(line::qn::jmt_feature_set(method).has(Feature::Cache), method);
}

TEST_CASE("jmt_feature_set: a cache model reaches the simulator and not the analytical engine") {
    const SN sn = cache_struct();
    CHECK(accepts(line::qn::jmt_feature_set("jsim"), sn));
    for (const std::string& method : kJmvaExact)
        CHECK_MESSAGE(!accepts(line::qn::jmt_feature_set(method), sn), method);
    for (const std::string& method : kJmvaClosedOnly)
        CHECK_MESSAGE(!accepts(line::qn::jmt_feature_set(method), sn), method);
    // and the refusal NAMES the cache rather than reporting a bare no
    const std::string reason =
        line::autosolver::auto_family_refusal("jmt", sn, std::string("jmva"));
    CHECK(reason.find("Cache") != std::string::npos);
}

TEST_CASE("jmt_feature_set: an open model keeps the exact engine and loses the rest") {
    const SN open_sn = mm1().get_struct();
    for (const std::string& method : kJmvaExact)
        CHECK_MESSAGE(accepts(line::qn::jmt_feature_set(method), open_sn), method);
    for (const std::string& method : kJmvaClosedOnly)
        CHECK_MESSAGE(!accepts(line::qn::jmt_feature_set(method), open_sn), method);
    // a priority model reaches the simulator alone
    const SN hol_sn = open_hol().get_struct();
    CHECK(accepts(line::qn::jmt_feature_set("jsim"), hol_sn));
    for (const std::string& method : kJmvaExact)
        CHECK_MESSAGE(!accepts(line::qn::jmt_feature_set(method), hol_sn), method);
}

// ---------------------------------------------------------------------------
// SolverJMT: what the registry CANNOT name
// ---------------------------------------------------------------------------

TEST_CASE("jmt_method_refusal: the server count is structural") {
    const SN single = closed_single_server().get_struct();
    const SN multi = closed_multiserver().get_struct();
    line::jmt::JmtOptions opt;
    for (const std::string& method : kJmvaClosedOnly) {
        opt.method = method;
        CHECK_MESSAGE(line::jmt::jmt_method_refusal(single, method, opt).empty(), method);
        const std::string reason = line::jmt::jmt_method_refusal(multi, method, opt);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK(reason.find("multi-server") != std::string::npos);
    }
    // The exact engine writes an <ldstation> for it and is unaffected.
    for (const std::string& method : kJmvaExact) {
        opt.method = method;
        CHECK_MESSAGE(line::jmt::jmt_method_refusal(multi, method, opt).empty(), method);
    }
    // The report carries the same sentence.
    const std::string offered =
        line::autosolver::auto_family_refusal("jmt", multi, std::string("jmva.amva"));
    CHECK(offered.find("multi-server") != std::string::npos);
    CHECK(line::autosolver::auto_family_refusal("jmt", multi, std::string("jmva")).empty());
}

TEST_CASE("jmt_method_refusal: replication needs a finite horizon") {
    const SN sn = mm1().get_struct();
    line::jmt::JmtOptions unbounded;
    unbounded.method = "replication";
    const std::string reason = line::jmt::jmt_method_refusal(sn, "replication", unbounded);
    CHECK_FALSE(reason.empty());
    CHECK(reason.find("finite timespan") != std::string::npos);

    line::jmt::JmtOptions bounded = unbounded;
    bounded.max_simulated_time = 10.0;
    CHECK(line::jmt::jmt_method_refusal(sn, "replication", bounded).empty());

    // The report probes with the default options, so the row is not offered.
    CHECK_FALSE(
        line::autosolver::auto_family_refusal("jmt", sn, std::string("replication")).empty());
    CHECK(line::autosolver::auto_family_refusal("jmt", sn, std::string("jsim")).empty());
}

TEST_CASE("jmt_method_refusal: a binding buffer is refused by BOTH engines") {
    // NEITHER engine carries it, for opposite reasons: JSIM because no JMT drop
    // strategy reproduces LINE's blocking -- "waiting queue" does not enforce
    // the size at all and "BAS blocking" completes the service first -- and
    // JMVA because its document has no capacity element, so it answered 2.19
    // jobs at a station that can hold 2 (exact: 1.33).
    const SN sn = closed_binding_buffer().get_struct();
    line::jmt::JmtOptions opt;
    const std::string jsim_reason = line::jmt::jmt_method_refusal(sn, "jsim", opt);
    CHECK_FALSE(jsim_reason.empty());
    CHECK(jsim_reason.find("binds for the closed class") != std::string::npos);
    // ONE PREDICATE, TWO CALLERS: the sentence is the JSIM WRITER's own.
    CHECK(jsim_reason == line::io::jmt_buffer_capacity_refusal(sn, false));

    for (const std::string& method : kJmvaExact) {
        const std::string reason = line::jmt::jmt_method_refusal(sn, method, opt);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK(reason.find("no capacity element") != std::string::npos);
    }
    for (const std::string& method : kJmvaClosedOnly)
        CHECK_MESSAGE(!line::jmt::jmt_method_refusal(sn, method, opt).empty(), method);

    // and the report carries it
    CHECK_FALSE(
        line::autosolver::auto_family_refusal("jmt", sn, std::string("jsim")).empty());
    CHECK_FALSE(
        line::autosolver::auto_family_refusal("jmt", sn, std::string("jmva")).empty());
}

TEST_CASE("jmt_method_refusal: an open loss buffer keeps the simulator") {
    // The over-tightening guard: a refused OPEN arrival is LOST, which JMT's
    // queue section expresses directly, so jsim genuinely simulates an M/M/1/K.
    const SN sn = open_binding_buffer().get_struct();
    line::jmt::JmtOptions opt;
    CHECK(line::jmt::jmt_method_refusal(sn, "jsim", opt).empty());
    CHECK(line::jmt::jmt_method_refusal(sn, "default", opt).empty());
    CHECK(line::autosolver::auto_family_refusal("jmt", sn, std::string("jsim")).empty());
    // The analytical engine still cannot: it has no capacity element at all.
    for (const std::string& method : kJmvaExact) {
        const std::string reason = line::jmt::jmt_method_refusal(sn, method, opt);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK(reason.find("no capacity element") != std::string::npos);
    }
}

TEST_CASE("jmt_method_refusal: an uncapped model is untouched by the buffer rule") {
    // refresh_capacity DERIVES a finite sn.cap for a station nobody capped, so
    // the rule tests the cap against the population that can REACH the station;
    // reading "finite" as "binding" would refuse every closed model.
    line::jmt::JmtOptions opt;
    const SN single = closed_single_server().get_struct();
    const SN two = cqn_two_class().get_struct();
    CHECK(line::jmt::jmt_method_refusal(single, "jsim", opt).empty());
    CHECK(line::jmt::jmt_method_refusal(two, "jsim", opt).empty());
    for (const std::string& method : kJmvaExact) {
        CHECK_MESSAGE(line::jmt::jmt_method_refusal(single, method, opt).empty(), method);
        CHECK_MESSAGE(line::jmt::jmt_method_refusal(two, method, opt).empty(), method);
    }
    // A method that names NEITHER engine gets no verdict: write_jmva is also
    // SolverQNS's writer and calls the predicate with QNS's own method names.
    CHECK(line::jmt::jmt_method_refusal(closed_binding_buffer().get_struct(), "conway",
                                        opt).empty());
}

TEST_CASE("jmt gate: the methods a model can genuinely run are still offered") {
    const SN closed = closed_single_server().get_struct();
    for (const std::string& method : kJmvaExact)
        CHECK_MESSAGE(line::autosolver::auto_family_refusal("jmt", closed, method).empty(), method);
    for (const std::string& method : kJmvaClosedOnly)
        CHECK_MESSAGE(line::autosolver::auto_family_refusal("jmt", closed, method).empty(), method);
    CHECK(line::autosolver::auto_family_refusal("jmt", closed, std::string("jsim")).empty());
    // An open product-form network keeps the exact engine and the simulator.
    const SN open_sn = mm1().get_struct();
    for (const std::string& method : kJmvaExact)
        CHECK_MESSAGE(line::autosolver::auto_family_refusal("jmt", open_sn, method).empty(), method);
    CHECK(line::autosolver::auto_family_refusal("jmt", open_sn, std::string("jsim")).empty());
}

// ---------------------------------------------------------------------------
// SolverMAM: dec.mmap and retrial
// ---------------------------------------------------------------------------

TEST_CASE("mam_feature_set: dec.mmap is an open-network method with no INF branch") {
    const FeatureSet base = line::qn::mam_feature_set("dec.source");
    CHECK(base.has(Feature::ClosedClass));
    CHECK(base.has(Feature::SchedStrategy_INF));

    const FeatureSet f = line::qn::mam_feature_set("dec.mmap");
    CHECK_FALSE(f.has(Feature::ClosedClass));
    CHECK_FALSE(f.has(Feature::SelfLoopingClass));
    CHECK_FALSE(f.has(Feature::SchedStrategy_INF));
    // ... and it keeps the open shape it is written for
    CHECK(f.has(Feature::OpenClass));
    CHECK(f.has(Feature::SchedStrategy_FCFS));
    CHECK(f.has(Feature::SchedStrategy_HOL));
    CHECK(f.has(Feature::SchedStrategy_PS));
    CHECK(f.has(Feature::Queue));
    CHECK(f.has(Feature::Source));

    CHECK_FALSE(accepts(f, closed_single_server().get_struct()));
    CHECK(accepts(f, mm1().get_struct()));
    CHECK(accepts(f, open_hol().get_struct()));
}

TEST_CASE("mam_feature_set: dec.mmap has no fork-join branch either") {
    // The sweep has no synchronization, so a Fork is not a smaller model: the
    // dispatch sends this shape to solver_mam_fj from 'default' and
    // 'dec.source', and dec.mmap has no such route.
    const FeatureSet f = line::qn::mam_feature_set("dec.mmap");
    CHECK_FALSE(f.has(Feature::Fork));
    CHECK_FALSE(f.has(Feature::Join));
    CHECK_FALSE(f.has(Feature::Forker));
    CHECK_FALSE(f.has(Feature::Joiner));
    CHECK(line::qn::mam_feature_set("dec.source").has(Feature::Fork));

    const SN fj = open_forkjoin().get_struct();
    CHECK_FALSE(accepts(f, fj));
    const std::string reason =
        line::autosolver::auto_family_refusal("mam", fj, std::string("dec.mmap"));
    CHECK_FALSE(reason.empty());
    // ... and the routes that ARE written for it stay
    CHECK(line::autosolver::auto_family_refusal("mam", fj, std::string("default")).empty());
    CHECK(line::autosolver::auto_family_refusal("mam", fj, std::string("dec.source")).empty());
}

TEST_CASE("mam gate: dec.mmap is not offered on a closed model but stays on an open one") {
    const SN closed = closed_single_server().get_struct();
    const std::string reason =
        line::autosolver::auto_family_refusal("mam", closed, std::string("dec.mmap"));
    CHECK_FALSE(reason.empty());
    // feature_set_supports words its refusals as PHRASES rather than registry
    // names, and both offending declarations have to appear: the closed class
    // mix and the Delay's INF discipline.
    CHECK(reason.find("closed job class") != std::string::npos);
    CHECK(reason.find("INF scheduling discipline") != std::string::npos);
    CHECK(line::autosolver::auto_family_refusal("mam", mm1().get_struct(),
                                                std::string("dec.mmap")).empty());
    CHECK(line::autosolver::auto_family_refusal("mam", open_hol().get_struct(),
                                                std::string("dec.mmap")).empty());
    // The analyzer refuses the same model by name, so the two cannot drift.
    line::mam::MamOptions opt;
    opt.method = "dec.mmap";
    CHECK_THROWS_AS(line::mam::solver_mam_solve(closed, opt), line::UnsupportedError);
}

TEST_CASE("mam_retrial_refusal: an absent orbit is a MUST BE PRESENT rule") {
    // A feature set states which constructs are ACCEPTED, so it can refuse a
    // model for HAVING something and never for LACKING it: this rule can only
    // be structural.
    const SN open_sn = mm1().get_struct();
    const std::string open_reason = line::mam::mam_retrial_refusal(open_sn);
    CHECK_FALSE(open_reason.empty());
    CHECK(open_reason.find("retrial") != std::string::npos);
    // The runner and the report repeat the analyzer's own sentence.
    CHECK(line::mam::mam_model_method_refusal(open_sn, std::string("retrial")) == open_reason);
    CHECK(line::autosolver::auto_family_refusal("mam", open_sn, std::string("retrial")) ==
          open_reason);
    // A closed model is refused for the reason qsys_is_retrial reports.
    const std::string closed_reason = line::mam::mam_retrial_refusal(closed_single_server().get_struct());
    CHECK(closed_reason.find("open queueing model") != std::string::npos);

    line::mam::MamOptions opt;
    opt.method = "retrial";
    CHECK_THROWS_AS(line::mam::solver_mam_solve(open_sn, opt), line::UnsupportedError);
}

TEST_CASE("mam gate: the methods a model can genuinely run are still offered") {
    const SN closed = closed_single_server().get_struct();
    const std::string kept[] = {"default", "dec.source", "dec.poisson", "ldqbd"};
    for (const std::string& method : kept)
        CHECK_MESSAGE(line::autosolver::auto_family_refusal("mam", closed, method).empty(), method);
    const SN open_sn = mm1().get_struct();
    const std::string kept_open[] = {"default", "dec.source", "dec.poisson", "dec.mmap"};
    for (const std::string& method : kept_open)
        CHECK_MESSAGE(line::autosolver::auto_family_refusal("mam", open_sn, method).empty(), method);
}
