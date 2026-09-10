/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The per-method support gate of SolverNC.
 *
 * WHAT THIS PINS. `model.help()` / `findSolver()` reports one row per (solver,
 * method) pair, and it builds the nc rows from `nc::list_valid_methods()`
 * screened by `qn::nc_feature_set(method)` and `autosolver::auto_family_refusal`.
 * Until the rules below reached that route the report was far weaker than what
 * the RUNNER enforces, so it offered pairs that then threw: measured on a plain
 * M/M/1, 15 of the 42 nc.* rows called runnable raised, and on a single-class
 * closed network 10 did.
 *
 * TWO ROUTES REACH THE REPORT, and both are asserted here.
 *
 *   - What the feature registry CAN name rides in `qn::nc_feature_set(method)`:
 *     a closed population for 'is' and the six load-dependent evaluators (drop
 *     OpenClass), no think time for 'divdiff' (drop SchedStrategy_INF). A
 *     feature set says "I ACCEPT this construct", so it can refuse a model for
 *     HAVING one and never for lacking one -- which is why the next group cannot
 *     live there.
 *   - What it CANNOT -- requires a cache, requires state-dependent routing,
 *     requires a loss network, requires exactly two stations, requires normal
 *     usage -- is `nc::nc_method_refusal`, ONE predicate that `solver_nc_solve`
 *     throws on and that `auto_family_refusal` asks, so the gate and the run are
 *     one body of rules rather than two copies that drift.
 *
 * Each case asserts the refusal AND its converse: a model the method IS derived
 * for must keep it. Over-tightening hides an answer the caller could have had,
 * which is the same defect with the sign flipped.
 *
 * The MATLAB twin is matlab/src/solvers/NC/nc_method_refusal.m, the python twin
 * python/tests/test_gate_nc.py and the JAR twin
 * jar/src/test/java/jline/solvers/nc/SolverNCGateTest.java.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/nc/solver_nc_runner.h"

using line::lang::SchedStrategy;
using line::qn::Feature;
using line::qn::feature_set_supports;
using line::qn::NetworkStruct;
using line::qn::used_lang_features;

namespace {

using NetD = line::qn::Network<double>;
using Dist = line::lang::Distrib<double>;
using SN = NetworkStruct<double>;

/** Source -> FCFS Queue -> Sink, one OPEN class. */
NetD mm1() {
    NetD m("ncGateMM1");
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

/** Delay -> FCFS Queue, one closed class, N = 3. Saturated: the queue's offered
 *  load 1.5 exceeds its unit saturation rate, so the model is NOT in normal
 *  usage and PANACEA's asymptotic expansion does not apply to it. */
NetD repairmen() {
    NetD m("ncGateRepairmen");
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

/** The same closed shape with a three-server queue, N = 4: mu(n) = min(n,3) is 3
 *  at saturation against an offered load of 2, so this one IS in normal usage
 *  and 'panald' must survive the gate. */
NetD multiserver() {
    NetD m("ncGateMultiserver");
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

/** The repairmen shape with a rate LATTICE on the queue, N = 4. It matters
 *  because `pfqn_ncld` evaluates 'pana' and 'panald' with the same
 *  `pfqn_panaceald`, so on a load-dependent model the load-INDEPENDENT name
 *  reaches the load-dependent expansion. Saturation rate 1.9 against an offered
 *  load of 2, so this one is NOT in normal usage and both names must be refused;
 *  every other nc method on it runs. */
NetD loaddep() {
    NetD m("ncGateLoadDep");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 4.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    std::vector<double> alpha;
    alpha.push_back(1.0);
    alpha.push_back(1.6);
    alpha.push_back(1.8);
    alpha.push_back(1.9);
    m.set_load_dependence(q, alpha);
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Q1 -> Q2, one closed class, N = 4: TWO queueing stations, which is
 *  what the mmint2/gleint/comomld recursions are not stated for. */
NetD cqn3() {
    NetD m("ncGateCqn3");
    const std::size_t d = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 4.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

/** Source -> Q1 -> Q2 -> Sink: two queueing stations but NO closed population,
 *  so `pfqn_nc` answers with the exact open formulas before its method switch
 *  and the single-station rule must stay inactive. */
NetD open_tandem() {
    NetD m("ncGateOpenTandem");
    const std::size_t s = m.add_source("S");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, Dist::exp_rate(0.5));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    return m;
}

/** Queue1 -> Queue2, one closed class, NO delay station: zero think time, which
 *  is the shape the divided-difference closed form 'divdiff' is derived for
 *  (Casale, SIGMETRICS 2017, Eqs. 15-16). */
NetD cyclic_delay_free() {
    NetD m("ncGateCyclic");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2.0, q1);
    m.set_service(q1, c, Dist::exp_rate(1.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    return m;
}

/** Does the method's declared envelope accept everything `sn` uses? */
bool accepts(const std::string& method, const SN& sn) {
    return feature_set_supports("SolverNC", line::qn::nc_feature_set(method),
                                used_lang_features(sn))
        .ok;
}

/** The reason the REPORT would print, "" when the pair is offered. It is what
 *  `auto_family_refusal` composes: the feature set first, the structural
 *  predicate second. */
std::string gate(const SN& sn, const std::string& method) {
    const line::qn::SupportResult r = feature_set_supports(
        "SolverNC", line::qn::nc_feature_set(method), used_lang_features(sn));
    if (!r.ok) return r.reason;
    return line::nc::nc_method_refusal(sn, method);
}

/** The predicate asked the RUN's question: "" when the reference performs the
 *  method by name, even where the report declines to offer it. */
std::string run_question(const SN& sn, const std::string& method) {
    return line::nc::nc_method_refusal(sn, method, false, /*for_report=*/false);
}

/** Did the RUNNER refuse the same pair? The message when it did, else "". */
std::string runner_refusal(const SN& sn, const std::string& method) {
    line::nc::NcSolverOptions opt;
    opt.method = method;
    try {
        line::nc::solver_nc_run_analyzer(sn, opt);
    } catch (const std::exception& e) {
        return std::string(e.what());
    }
    return "";
}

bool mentions(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

}  // namespace

// ---------------------------------------------------------------------------
// What the registry CANNOT name: nc::nc_method_refusal
// ---------------------------------------------------------------------------

TEST_CASE("nc gate: the cache tokens name the missing Cache node") {
    NetD m = mm1();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"rayint", "spm"};
    for (int i = 0; i < 2; ++i) {
        const std::string reason = gate(sn, tokens[i]);
        CHECK_FALSE(reason.empty());
        CHECK(mentions(reason, "Cache node"));
        CHECK(mentions(reason, "SPM saddle point"));
    }
}

TEST_CASE("nc gate: the loss-network tokens name the missing region") {
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    CHECK(mentions(gate(sn, "ms"), "loss network"));
    CHECK(mentions(gate(sn, "erlangfp"), "loss network"));
}

TEST_CASE("nc gate: rec names both routes it has") {
    NetD m = repairmen();
    const std::string reason = gate(m.get_struct(), "rec");
    CHECK_FALSE(reason.empty());
    CHECK(mentions(reason, "MDD-rec"));
    CHECK(mentions(reason, "Petri net"));
}

TEST_CASE("nc gate: sdr names the routing the model does not declare") {
    NetD m = mm1();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"sdr", "sdr.mva"};
    for (int i = 0; i < 2; ++i) {
        const std::string reason = gate(sn, tokens[i]);
        CHECK_FALSE(reason.empty());
        CHECK(mentions(reason, "state-dependent routing"));
    }
}

TEST_CASE("nc gate: morrison names the shape it is derived for") {
    NetD m = repairmen();
    const std::string reason = gate(m.get_struct(), "morrison");
    CHECK_FALSE(reason.empty());
    CHECK(mentions(reason, "DPS"));
    CHECK(mentions(reason, "two stations"));
}

TEST_CASE("nc gate: panald names normal usage on a saturated closed model") {
    NetD m = repairmen();
    const std::string reason = gate(m.get_struct(), "panald");
    CHECK_FALSE(reason.empty());
    CHECK(mentions(reason, "normal usage"));
}

TEST_CASE("nc gate: pana is refused on a load-dependent model outside normal usage") {
    // `pfqn_ncld`'s case label covers both names: on a rate lattice the
    // load-independent NAME is evaluated by the load-dependent kernel, so it
    // throws the panald refusal. The gate has to know that aliasing.
    NetD m = loaddep();
    const SN& sn = m.get_struct();
    const std::string reason = gate(sn, "pana");
    CHECK_FALSE(reason.empty());
    CHECK(mentions(reason, "normal usage"));
    CHECK(mentions(reason, "evaluates it as 'panald'"));
    CHECK_FALSE(runner_refusal(sn, "pana").empty());
}

TEST_CASE("nc gate: the load-dependent model keeps every other method") {
    // The converse, and the whole point of not gating 'pana' off the lattice.
    NetD m = loaddep();
    const SN& sn = m.get_struct();
    const char* kept[] = {"default", "exact",  "ca",  "clw", "comom", "comomld", "le",
                          "ble",    "mmint2", "gleint", "kt", "bkt", "lekt", "cub", "rd",  "nrl",
                          "nrp",     "nre",    "propfair", "ger", "rgf"};
    for (int i = 0; i < 21; ++i) {
        INFO("method ", kept[i]);
        CHECK(gate(sn, kept[i]).empty());
    }
}

TEST_CASE("nc gate: pana is left alone off the load-dependent route") {
    // Without a rate lattice 'pana' takes its own `pfqn_nc` arm, which warns
    // and returns an empty constant rather than throwing, so the gate must not
    // refuse it there -- that would be over-tightening on a path this change is
    // not about.
    NetD a = repairmen();
    NetD b = multiserver();
    CHECK(gate(a.get_struct(), "pana").empty());
    CHECK(gate(b.get_struct(), "pana").empty());
}

TEST_CASE("nc gate: the single-station recursions name the station count") {
    // `pfqn_nc` states "a model with a delay and a single queueing station" for
    // mmint2/gleint, and `pfqn_comomrm_ld` refuses with "accepts at most a single
    // queueing station". Neither is a feature name: it is a COUNT.
    NetD m = cqn3();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"mmint2", "gleint", "comomld"};
    for (int i = 0; i < 3; ++i) {
        const std::string reason = gate(sn, tokens[i]);
        INFO("method ", tokens[i]);
        CHECK_FALSE(reason.empty());
        CHECK(mentions(reason, "single queueing station"));
        CHECK(mentions(reason, "has 2"));
    }
}

TEST_CASE("nc gate: mmint2 and gleint are gated for the report only") {
    // THE RULING (2026-07-25, register row N1, reaffirmed when this gate was
    // added): the report answers "should this be offered" and the run answers
    // "what does the reference do". `pfqn_nc` answers mmint2/gleint outside their
    // shape with an empty constant and a ZERO TABLE, which `test_nc.cpp`'s "a
    // method outside its own domain returns the reference's zero table" pins, so
    // a caller who names the method keeps that answer while the report stops
    // offering it.
    //
    // comomld is NOT in that bucket: `pfqn_comomrm_ld` refuses natively, so it is
    // refused on both paths.
    NetD m = cqn3();
    const SN& sn = m.get_struct();
    const char* report_only[] = {"mmint2", "gleint"};
    for (int i = 0; i < 2; ++i) {
        INFO("method ", report_only[i]);
        CHECK_FALSE(gate(sn, report_only[i]).empty());        // the report declines
        CHECK(run_question(sn, report_only[i]).empty());      // the run does not
        CHECK(runner_refusal(sn, report_only[i]).empty());    // and it really runs
    }
    CHECK_FALSE(gate(sn, "comomld").empty());
    CHECK_FALSE(run_question(sn, "comomld").empty());
    CHECK_FALSE(runner_refusal(sn, "comomld").empty());
}

TEST_CASE("nc gate: the single-station recursions survive on one queueing station") {
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"mmint2", "gleint", "comomld"};
    for (int i = 0; i < 3; ++i) {
        INFO("method ", tokens[i]);
        CHECK(gate(sn, tokens[i]).empty());
    }
}

TEST_CASE("nc gate: the station-count rule is inactive without a closed population") {
    // An open network never reaches `pfqn_nc`'s method switch, so a two-queue
    // OPEN tandem runs these names correctly and must keep them.
    NetD m = open_tandem();
    const SN& sn = m.get_struct();
    CHECK(gate(sn, "mmint2").empty());
    CHECK(gate(sn, "gleint").empty());
}

TEST_CASE("nc gate: panald survives where the expansion does apply") {
    // The converse: a saturation rate of 3 against an offered load of 2 IS
    // normal usage, so the row must stay and the run must go through.
    NetD m = multiserver();
    const SN& sn = m.get_struct();
    CHECK(gate(sn, "panald").empty());
    CHECK(runner_refusal(sn, "panald").empty());
}

// ---------------------------------------------------------------------------
// What the registry CAN name: qn::nc_feature_set(method)
// ---------------------------------------------------------------------------

TEST_CASE("nc gate: divdiff refuses a think time by feature name") {
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    CHECK_FALSE(accepts("divdiff", sn));
    CHECK(accepts("default", sn));
    CHECK_FALSE(line::qn::nc_feature_set("divdiff").has(Feature::SchedStrategy_INF));
    CHECK(line::qn::nc_feature_set("default").has(Feature::SchedStrategy_INF));
}

TEST_CASE("nc gate: divdiff survives on a delay-free closed model") {
    NetD m = cyclic_delay_free();
    const SN& sn = m.get_struct();
    CHECK(accepts("divdiff", sn));
    CHECK(gate(sn, "divdiff").empty());
    CHECK(runner_refusal(sn, "divdiff").empty());
}

TEST_CASE("nc gate: the closed-population methods refuse an open chain") {
    NetD m = mm1();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"is", "rd", "nrp", "nrl", "nre", "comomld", "panald"};
    for (int i = 0; i < 7; ++i) {
        CHECK_FALSE(line::qn::nc_feature_set(tokens[i]).has(Feature::OpenClass));
        CHECK_FALSE(accepts(tokens[i], sn));
        // The reason names the feature in the registry's own words.
        CHECK(mentions(gate(sn, tokens[i]), "an open job class"));
    }
    CHECK(line::qn::nc_feature_set("default").has(Feature::OpenClass));
}

TEST_CASE("nc gate: the closed-population methods survive on a closed model") {
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    const char* tokens[] = {"is", "rd", "nrp", "nrl", "nre", "comomld"};
    for (int i = 0; i < 6; ++i) {
        CHECK(accepts(tokens[i], sn));
        CHECK(gate(sn, tokens[i]).empty());
    }
}

// ---------------------------------------------------------------------------
// One predicate, two callers
// ---------------------------------------------------------------------------

TEST_CASE("nc gate: every refused pair is refused by the runner too") {
    NetD models[6] = {mm1(), repairmen(), cyclic_delay_free(), loaddep(), cqn3(),
                      open_tandem()};
    const std::vector<std::string> declared = line::nc::list_valid_methods();
    for (int mi = 0; mi < 6; ++mi) {
        const SN& sn = models[mi].get_struct();
        for (std::size_t i = 0; i < declared.size(); ++i) {
            if (gate(sn, declared[i]).empty()) continue;
            // A REPORT-ONLY refusal, and the ruling says so: the reference performs
            // this one by name (a warning and a zero table), which is exactly why
            // the report declines to offer it.
            if (run_question(sn, declared[i]).empty()) continue;
            const std::string msg = runner_refusal(sn, declared[i]);
            INFO("model ", mi, " method ", declared[i]);
            CHECK_FALSE(msg.empty());
        }
    }
}

TEST_CASE("nc gate: a closed product-form network keeps the methods it can run") {
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    const char* kept[] = {"default", "exact", "ca",   "comom", "le",  "ble", "mmint2",
                          "gleint",  "pana", "propfair", "cub", "kt", "bkt", "lekt", "clw"};
    for (int i = 0; i < 15; ++i) {
        INFO("method ", kept[i]);
        CHECK(gate(sn, kept[i]).empty());
    }
}

TEST_CASE("nc gate: the AUTO report route carries the same verdicts") {
    // `auto_family_refusal` is the call `auto_find_solver` makes for every row,
    // so the report and the predicate must agree token for token.
    NetD m = repairmen();
    const SN& sn = m.get_struct();
    CHECK(mentions(line::autosolver::auto_family_refusal("nc", sn, "rec"), "MDD-rec"));
    CHECK(mentions(line::autosolver::auto_family_refusal("nc", sn, "morrison"), "DPS"));
    CHECK(line::autosolver::auto_family_refusal("nc", sn, "comom").empty());
    CHECK(mentions(line::autosolver::auto_family_refusal("nc", sn, "divdiff"),
                   "INF scheduling discipline"));
}
