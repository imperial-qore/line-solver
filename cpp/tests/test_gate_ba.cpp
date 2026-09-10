/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The per-method support gate of SolverBA.
 *
 * WHAT THIS PINS. `model.help()` / `findSolver()` reports one row per (solver,
 * method) pair, and it builds the ba rows from `ba::list_valid_methods(sn)`
 * screened by `qn::ba_feature_set(method)` -- the route
 * `autosolver::auto_family_methods` and `auto_family_refusal` take. Until the
 * rules below reached that route, the report was much weaker than what the
 * ANALYZER enforces at run time, so it offered pairs that then threw: measured
 * on a two-class closed network, 30 of the 36 ba.* rows called runnable threw,
 * and on a closed multiserver model 31 of 38 did.
 *
 * TWO ROUTES REACH THE REPORT, and both are asserted here. What the feature
 * registry can NAME rides in `qn::ba_feature_set(method)`: a delay station for
 * sb/harel/sib/scb/lr, a closed class for the three open families. What it
 * cannot -- a class count and a server count -- is `ba::method_refusal`, ONE
 * predicate that `solver_ba_analyzer` throws on and that `list_valid_methods(sn)`
 * projects, so the gate and the run are one body of rules rather than two
 * copies. The four codebases agree on what a refusal DOES since the 2026-07-25
 * zero-table convention was overturned (user ruling 2026-09-04): every family
 * refuses by name, and `ba_out_of_domain` is gone with it.
 *
 * A SECOND PREDICATE, `ba::method_degenerate`, answers a different question:
 * inside the domain, does the formula still say anything. Only the report reads
 * it, because a vacuous bound is a VALID bound and a caller naming it is
 * entitled to the answer.
 *
 * Each case asserts the refusal AND its converse: a model the bounds ARE derived
 * for must keep every method. Over-tightening hides a bound the user could have
 * had, which is the same defect with the sign flipped.
 *
 * The MATLAB twin is matlab/src/solvers/BA/ba_method_refusal.m with
 * @SolverBA/SolverBA.m; the python twin is python/tests/test_gate_ba.py and the
 * JAR twin jar/src/test/java/jline/solvers/ba/SolverBAGateTest.java.
 */
#include <algorithm>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/feature_set.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ba/solver_ba_runner.h"

using line::UnsupportedError;
using line::lang::SchedStrategy;
using line::qn::Feature;
using line::qn::feature_set_supports;
using line::qn::FeatureSet;
using line::qn::NetworkStruct;
using line::qn::used_lang_features;

namespace {

using NetD = line::qn::Network<double>;
using Dist = line::lang::Distrib<double>;
using SN = NetworkStruct<double>;

/** Queue1 -> Queue2, ONE closed class, no delay, one server each: the shape
 *  every bound family in this solver is derived for. */
NetD cyclic_delay_free() {
    NetD m("baGateCyclic");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 3.0, q1);
    m.set_service(q1, c, Dist::exp_rate(1.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    return m;
}

/** The same shape with TWO closed classes, still delay-free, so the CLASS
 *  premise is the only one this model violates. */
NetD cyclic_two_class() {
    NetD m("baGate2Class");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, q1);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, q1);
    m.set_service(q1, c1, Dist::exp_rate(1.0));
    m.set_service(q1, c2, Dist::exp_rate(2.0));
    m.set_service(q2, c1, Dist::exp_rate(2.0));
    m.set_service(q2, c2, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, q1, 1.0);
    P.set(c2, c2, q1, q2, 1.0);
    P.set(c2, c2, q2, q1, 1.0);
    m.link(P);
    return m;
}

/** Single-class closed with a three-server station: only ssd, ldbcmp and the
 *  auto composite that draws on them survive. */
NetD multiserver() {
    NetD m("baGateMultiserver");
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

/** Single-class closed WITH a delay station: the think-time premise of
 *  harel/sb/scb/sib is what this one violates, and nothing else. */
NetD with_delay() {
    NetD m("baGateDelay");
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

/**
 * Two delay stations and nothing else: no queueing station at all, so the
 * ldbcmp bottleneck the open-network occupancy is built on does not exist. The
 * general form of the shape a Petri net presents, every Place being an INF
 * station -- and the shape that made `method_degenerate` THROW from inside
 * `list_valid_methods` before it was guarded.
 */
NetD all_delay_closed() {
    NetD m("baGateAllDelay");
    const std::size_t d1 = m.add_delay("D1");
    const std::size_t d2 = m.add_delay("D2");
    const std::size_t c = m.add_closed_class("C", 2.0, d1);
    m.set_service(d1, c, Dist::exp_rate(1.0));
    m.set_service(d2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d1, d2, 1.0);
    P.set(d2, d1, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue -> Sink, one open class, with the two laws as given. */
NetD open_mm1(const Dist& arrival, const Dist& service, const std::string& name) {
    NetD m(name);
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, arrival);
    m.set_service(q, c, service);
    line::qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

NetD mm1() { return open_mm1(Dist::exp_rate(1.0), Dist::exp_rate(2.0), "baGateMM1"); }

/** Erlang SERVICE: what all three open families refuse at a queueing station. */
NetD mm1_erlang_service() {
    return open_mm1(Dist::exp_rate(1.0), Dist::erlang(3.0 / 0.5, 3), "baGateErlSvc");
}

/**
 * Erlang SOURCE with exponential service. THE CONVERSE MODEL: 'snc' consumes
 * the arrival law and answers this one, so a gate that dropped the Erlang
 * feature outright would hide a bound the user could have had. 'bpt' and 'bgt'
 * must still be refused -- they read the mean alone and would bound the Poisson
 * system instead.
 */
NetD mm1_erlang_source() {
    return open_mm1(Dist::erlang(3.0 / 1.0, 3), Dist::exp_rate(2.0), "baGateErlSrc");
}

/**
 * Delay + two PS queues, one class, N = 4: exactly the ldbcmp regime boundary
 * N == Qhat, where the bound degenerates to the trivial X >= 0.
 */
NetD ldbcmp_boundary() {
    NetD m("baGateLdbcmp");
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

/** Does the method's declared envelope accept everything `sn` uses? */
bool accepts(const std::string& method, const SN& sn) {
    return feature_set_supports("SolverBA", line::qn::ba_feature_set(method),
                                used_lang_features(sn))
        .ok;
}

/** The methods the REPORT offers: the model-aware list screened by the feature
 *  set, which is what `auto_family_methods` + `auto_family_refusal` compose. */
std::vector<std::string> offered(const SN& sn) {
    std::vector<std::string> out;
    const std::vector<std::string> declared = line::ba::list_valid_methods(sn);
    for (std::size_t i = 0; i < declared.size(); ++i)
        if (accepts(declared[i], sn)) out.push_back(declared[i]);
    return out;
}

bool contains(const std::vector<std::string>& v, const std::string& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
}

bool mentions(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

}  // namespace

// ---------------------------------------------------------------------------
// What the registry CANNOT name: ba::method_refusal
// ---------------------------------------------------------------------------

TEST_CASE("ba gate: a multiclass model refuses the single-class families by name") {
    NetD m = cyclic_two_class();
    const SN& sn = m.get_struct();
    const std::string methods[] = {"aba.upper", "gb.lower",    "pbh.upper", "harel.lower",
                                   "ssd.upper", "ldbcmp.lower", "auto.upper", "default"};
    for (const std::string& method : methods) {
        const std::string reason = line::ba::method_refusal(sn, method);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK_MESSAGE(mentions(reason, "supports single-class closed networks only"), reason);
        CHECK_FALSE_MESSAGE(contains(line::ba::list_valid_methods(sn), method), method);
    }
    // The reason names the method that will RUN, not the alias asked for.
    CHECK(mentions(line::ba::method_refusal(sn, "default"), "'gb.upper'"));
}

TEST_CASE("ba gate: a multiserver model refuses the single-server families by name") {
    NetD m = multiserver();
    const SN& sn = m.get_struct();
    const std::string methods[] = {"aba.upper", "bjb.lower", "gb.upper",
                                   "pbh.upper", "bjbk.lower", "default"};
    for (const std::string& method : methods) {
        const std::string reason = line::ba::method_refusal(sn, method);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK_MESSAGE(mentions(reason, "does not support multi-server stations"), reason);
        // 'ssd' IS the multiserver bound, so the refusal points at it.
        CHECK_MESSAGE(mentions(reason, "use 'ssd'"), reason);
    }
}

TEST_CASE("ba gate: the multiclass chain families refuse a multiserver model, without naming ssd") {
    // mwba/cub/mbjb/looping survive a multiclass model but not a multiserver
    // one, and 'ssd' is single-class, so it is no alternative for them.
    NetD m = multiserver();
    const SN& sn = m.get_struct();
    const std::string methods[] = {"mwba.upper", "cub.upper", "mbjb.lower", "looping.lower"};
    for (const std::string& method : methods) {
        const std::string reason = line::ba::method_refusal(sn, method);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK_MESSAGE(mentions(reason, "does not support multi-server stations"), reason);
        CHECK_FALSE_MESSAGE(mentions(reason, "use 'ssd'"), reason);
    }
}

TEST_CASE("ba gate: an open model refuses the closed families by name") {
    NetD m = mm1();
    const SN& sn = m.get_struct();
    const std::string methods[] = {"aba.upper", "gb.upper", "mwba.upper", "cub.upper"};
    for (const std::string& method : methods) {
        const std::string reason = line::ba::method_refusal(sn, method);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK_MESSAGE(mentions(reason, "closed networks only"), reason);
    }
}

// ---------------------------------------------------------------------------
// What the registry CAN name: qn::ba_feature_set
// ---------------------------------------------------------------------------

TEST_CASE("ba_feature_set: the think-time families drop the delay-station names") {
    const FeatureSet base = line::qn::ba_feature_set("default");
    CHECK(base.has(Feature::SchedStrategy_INF));
    CHECK(base.has(Feature::Delay));
    CHECK(base.has(Feature::ClosedClass));

    const std::string delay_free[] = {"sb.upper", "sb.lower",  "harel.upper", "harel.lower",
                                      "sib.lower", "scb.upper", "lr.upper",   "lr"};
    for (const std::string& method : delay_free) {
        const FeatureSet f = line::qn::ba_feature_set(method);
        CHECK_FALSE_MESSAGE(f.has(Feature::SchedStrategy_INF), method);
        CHECK_FALSE_MESSAGE(f.has(Feature::Delay), method);
        CHECK_FALSE_MESSAGE(f.has(Feature::InfiniteServer), method);
        // ... and each keeps the closed core it does serve.
        CHECK_MESSAGE(f.has(Feature::ClosedClass), method);
        CHECK_MESSAGE(f.has(Feature::Queue), method);
        CHECK_MESSAGE(f.has(Feature::Exp), method);
    }
}

TEST_CASE("ba_feature_set: the three open families drop the closed-model names") {
    const std::string open_only[] = {"bpt.lower", "bgt.upper", "snc.upper"};
    for (const std::string& method : open_only) {
        const FeatureSet f = line::qn::ba_feature_set(method);
        CHECK_FALSE_MESSAGE(f.has(Feature::ClosedClass), method);
        CHECK_FALSE_MESSAGE(f.has(Feature::SchedStrategy_INF), method);
        CHECK_MESSAGE(f.has(Feature::OpenClass), method);
        CHECK_MESSAGE(f.has(Feature::Source), method);
    }
}

TEST_CASE("ba gate: a delay station is what refuses the think-time families on a real model") {
    NetD m = with_delay();
    const SN& sn = m.get_struct();
    const std::string refused[] = {"sb.upper",  "sb.lower",  "harel.upper",
                                   "harel.lower", "sib.lower", "scb.upper", "lr.upper"};
    for (const std::string& method : refused) {
        CHECK_FALSE_MESSAGE(accepts(method, sn), method);
        // ... and the STRUCTURAL half stays silent about it: a delay station is
        // a construct the model HAS, which is what a feature set is for.
        CHECK_MESSAGE(line::ba::method_refusal(sn, method).empty(), method);
    }
    const std::string kept[] = {"gb.upper", "aba.lower", "pbh.upper", "default"};
    for (const std::string& method : kept) CHECK_MESSAGE(accepts(method, sn), method);
}

// ---------------------------------------------------------------------------
// One predicate, several callers
// ---------------------------------------------------------------------------

TEST_CASE("ba gate: the run throws the sentence the gate reports") {
    NetD two = cyclic_two_class();
    NetD ms = multiserver();
    const SN& sn2 = two.get_struct();
    const SN& snm = ms.get_struct();
    struct Case {
        const SN* sn;
        const char* method;
    };
    // EVERY family throws its reason, the five noniterative single-class ones
    // included: they used to answer a model outside their class domain with a
    // zero table, and that convention was overturned on 2026-09-04.
    const Case cases[] = {{&sn2, "aba.upper"}, {&sn2, "harel.upper"},
                          {&snm, "gb.lower"},  {&snm, "cub.upper"}};
    for (const Case& c : cases) {
        const std::string reason = line::ba::method_refusal(*c.sn, c.method);
        CHECK_MESSAGE(!reason.empty(), c.method);
        line::ba::BaOptions opt;
        opt.method = c.method;
        // Same predicate, same string, whichever gate the caller meets first.
        try {
            line::ba::solver_ba_run_analyzer(*c.sn, opt);
            FAIL_CHECK("solver_ba_run_analyzer answered an inapplicable model");
        } catch (const UnsupportedError& e) {
            CHECK_MESSAGE(mentions(std::string(e.what()), reason), e.what());
        }
    }
}

TEST_CASE("ba gate: an inapplicable family throws instead of answering with zeros") {
    // The retired `ba_out_of_domain` used to answer here: aba/bjb/pb/sb/gb and
    // mwba returned zeroed metrics with an EMPTY C and X on a model outside
    // their class domain, reproducing a reference that has since changed. The
    // 2026-07-25 alignment was overturned on 2026-09-04, so the run refuses and
    // the report never names the pair in the first place.
    NetD two = cyclic_two_class();
    const SN& sn = two.get_struct();
    const std::string zeroing[] = {"aba.upper", "aba.lower", "bjb.upper",
                                   "pb.lower",  "sb.upper",  "gb.upper", "default"};
    for (const std::string& method : zeroing) {
        CHECK_FALSE_MESSAGE(contains(line::ba::list_valid_methods(sn), method), method);
        CHECK_FALSE_MESSAGE(contains(offered(sn), method), method);
        line::ba::BaOptions opt;
        opt.method = method;
        CHECK_THROWS_AS(line::ba::solver_ba_run_analyzer(sn, opt), UnsupportedError);
    }
}

TEST_CASE("ba gate: the model-aware list is a projection of the same predicate") {
    NetD models[] = {cyclic_delay_free(), cyclic_two_class(), multiserver(), with_delay(), mm1()};
    for (NetD& m : models) {
        const SN& sn = m.get_struct();
        const std::vector<std::string> declared = line::ba::list_valid_methods(sn);
        for (std::size_t i = 0; i < declared.size(); ++i)
            CHECK_MESSAGE(line::ba::method_refusal(sn, declared[i]).empty(), declared[i]);
    }
}

TEST_CASE("ba gate: the two copies of the alias table agree") {
    // `qn::ba_resolve_method_name` sits beside the feature set because
    // solver_feature_sets.h is included BY the solver and cannot include it
    // back. It is one line, and this is what keeps it honest.
    const std::string names[] = {"default", "auto", "lr", "qr", "gb.upper", "spnlp.upper"};
    for (const std::string& n : names)
        CHECK(line::ba::resolve_method(n) == line::qn::ba_resolve_method_name(n));
}

// ---------------------------------------------------------------------------
// The converse: the bounds still answer the models they are derived for
// ---------------------------------------------------------------------------

TEST_CASE("ba gate: a delay-free single-class closed network keeps every method") {
    // Over-tightening is as bad as the leak. The only names this model may lose
    // are the three OPEN families and the four spnlp ones, both of which were
    // already gated before this change.
    NetD m = cyclic_delay_free();
    const std::vector<std::string> listed = offered(m.get_struct());
    const std::string absent[] = {"bpt.lower",   "bgt.upper",   "snc.upper",     "spnlp.upper",
                                  "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower"};
    const std::vector<std::string> all = line::ba::list_valid_methods();
    for (std::size_t i = 0; i < all.size(); ++i) {
        const bool expected_absent =
            std::find(absent, absent + sizeof(absent) / sizeof(*absent), all[i]) !=
            absent + sizeof(absent) / sizeof(*absent);
        CHECK_MESSAGE(contains(listed, all[i]) == !expected_absent, all[i]);
    }
}

TEST_CASE("ba gate: the multiserver model keeps the bounds stated for it") {
    NetD m = multiserver();
    const std::vector<std::string> listed = offered(m.get_struct());
    const std::string kept[] = {"ssd.upper", "ssd.lower", "ldbcmp.lower", "auto.upper",
                                "auto.lower"};
    for (const std::string& method : kept) CHECK_MESSAGE(contains(listed, method), method);
}

TEST_CASE("ba gate: the two-class model keeps the multiclass bounds") {
    NetD m = cyclic_two_class();
    const std::vector<std::string> listed = offered(m.get_struct());
    const std::string kept[] = {"mwba.upper", "mwba.lower",    "cub.upper",
                                "mbjb.lower", "looping.upper", "looping.lower"};
    for (const std::string& method : kept) CHECK_MESSAGE(contains(listed, method), method);
}

TEST_CASE("ba gate: every offered method actually runs") {
    // The audit condition, in the small: no offered ba.* pair may throw.
    NetD models[] = {cyclic_delay_free(), cyclic_two_class(), multiserver(), with_delay(), mm1()};
    for (NetD& m : models) {
        const SN& sn = m.get_struct();
        const std::vector<std::string> list = offered(sn);
        CHECK_FALSE(list.empty());
        for (std::size_t i = 0; i < list.size(); ++i) {
            line::ba::BaOptions opt;
            opt.method = list[i];
            CHECK_NOTHROW(line::ba::solver_ba_run_analyzer(sn, opt));
        }
    }
}

// ---------------------------------------------------------------------------
// The exponential-service premise of the three open families
// ---------------------------------------------------------------------------

TEST_CASE("ba gate: Erlang service refuses all three open families") {
    NetD m = mm1_erlang_service();
    const SN& sn = m.get_struct();
    const std::string methods[] = {"bpt.lower", "bgt.upper", "snc.upper"};
    for (const std::string& method : methods) {
        const bool by_feature = !accepts(method, sn);
        const bool by_structure = !line::ba::method_refusal(sn, method).empty();
        const bool refused = by_feature || by_structure;
        CHECK_MESSAGE(refused, method);
    }
    CHECK(offered(sn).empty());
}

TEST_CASE("ba_feature_set: bpt and bgt drop every law but Exp") {
    // Registry-expressible: both read the mean alone, so a non-exponential law
    // ANYWHERE -- source included -- is silently bounded as if it were Poisson
    // rather than refused. Measured: swapping the Exp(1) source of an M/M/1 for
    // an Erlang of the same mean leaves bgt.upper at QLen 32.6667 and bpt.lower
    // at 1.0, digit for digit.
    const Feature laws[] = {Feature::APH,       Feature::Coxian,  Feature::Cox2,
                            Feature::Erlang,    Feature::HyperExp, Feature::PH,
                            Feature::Det,       Feature::Lognormal, Feature::Pareto,
                            Feature::Uniform,   Feature::Weibull};
    const std::string methods[] = {"bpt.lower", "bgt.upper"};
    for (const std::string& method : methods) {
        const FeatureSet f = line::qn::ba_feature_set(method);
        CHECK_MESSAGE(f.has(Feature::Exp), method);
        for (Feature law : laws) CHECK_FALSE_MESSAGE(f.has(law), method);
    }
    // snc keeps them: it CONSUMES the arrival law.
    const FeatureSet s = line::qn::ba_feature_set("snc.upper");
    for (Feature law : laws) CHECK(s.has(law));
}

TEST_CASE("ba gate: an Erlang source keeps snc and drops bpt and bgt") {
    // The converse of the delta above, and the reason snc's rule is structural:
    // no feature name can say "Erlang at a Queue but not at a Source", so
    // dropping the law would refuse a model snc answers.
    NetD m = mm1_erlang_source();
    const SN& sn = m.get_struct();
    CHECK(accepts("snc.upper", sn));
    CHECK(line::ba::method_refusal(sn, "snc.upper").empty());
    const std::vector<std::string> only_snc = offered(sn);
    CHECK(only_snc.size() == 1u);
    CHECK(contains(only_snc, "snc.upper"));
    CHECK_FALSE(accepts("bpt.lower", sn));
    CHECK_FALSE(accepts("bgt.upper", sn));
    // ... and it really does run, which is what makes the refusal of the other
    // two a judgement about correctness rather than about coverage.
    line::ba::BaOptions opt;
    opt.method = "snc.upper";
    CHECK_NOTHROW(line::ba::solver_ba_run_analyzer(sn, opt));
}

TEST_CASE("ba gate: the snc service rule is structural and skips the Source") {
    NetD svc = mm1_erlang_service();
    NetD src = mm1_erlang_source();
    CHECK(mentions(line::ba::method_refusal(svc.get_struct(), "snc.upper"),
                   "requires exponential service"));
    CHECK(line::ba::method_refusal(src.get_struct(), "snc.upper").empty());
    // bpt and bgt carry no structural rule at all: theirs is the feature set.
    CHECK(line::ba::method_refusal(svc.get_struct(), "bpt.lower").empty());
    CHECK(line::ba::method_refusal(svc.get_struct(), "bgt.upper").empty());
}

// ---------------------------------------------------------------------------
// A degenerate bound is not offered, but is still answered
// ---------------------------------------------------------------------------

TEST_CASE("ba gate: ldbcmp is withheld at its regime boundary") {
    NetD m = ldbcmp_boundary();
    const SN& sn = m.get_struct();
    const std::string why = line::ba::method_degenerate(sn, "ldbcmp.lower");
    CHECK_MESSAGE(mentions(why, "Qhat=4"), why);
    CHECK_MESSAGE(mentions(why, "N=4"), why);
    CHECK_MESSAGE(mentions(why, "trivial bound X >= 0"), why);
    // The APPLICABILITY predicate stays silent: the model is inside the
    // method's domain, which is exactly why this is a separate question.
    CHECK(line::ba::method_refusal(sn, "ldbcmp.lower").empty());
    CHECK_FALSE(contains(line::ba::list_valid_methods(sn), "ldbcmp.lower"));
    CHECK_FALSE(contains(offered(sn), "ldbcmp.lower"));
}

TEST_CASE("ba gate: the run still answers a degenerate bound when named") {
    // Not a refusal: X >= 0 IS a lower bound, so the run publishes it -- which
    // tests/test_ba.cpp pins at this boundary. What changed is that nothing
    // offers it.
    NetD m = ldbcmp_boundary();
    line::ba::BaOptions opt;
    opt.method = "ldbcmp.lower";
    const line::mva::AvgResult<double> r = line::ba::solver_ba_run_analyzer(m.get_struct(), opt);
    CHECK(r.TN(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("ba gate: ldbcmp survives where its bound says something") {
    NetD models[] = {with_delay(), multiserver()};
    for (NetD& m : models) {
        const SN& sn = m.get_struct();
        CHECK(line::ba::method_degenerate(sn, "ldbcmp.lower").empty());
        CHECK(contains(offered(sn), "ldbcmp.lower"));
    }
}

TEST_CASE("ba gate: a model with no queueing station is answered, not thrown") {
    // A PREDICATE MUST NOT THROW. `method_degenerate` is asked once per name by
    // `list_valid_methods`, which runs it before the Petri sieve has dropped
    // anything, so a model whose every station is an INF server reached
    // `pfqn_ldbcmp` with an empty demand vector and took the whole listing down
    // with it -- which is how `tests/test_spn_lpbnd.cpp` caught it.
    NetD m = all_delay_closed();
    const SN& sn = m.get_struct();
    const std::string why = line::ba::method_degenerate(sn, "ldbcmp.lower");
    CHECK_MESSAGE(mentions(why, "no queueing station"), why);
    const std::vector<std::string> listed = line::ba::list_valid_methods(sn);
    CHECK_FALSE(listed.empty());
    CHECK_FALSE(contains(listed, "ldbcmp.lower"));
}

TEST_CASE("ba gate: no other method has a degenerate regime") {
    NetD models[] = {cyclic_delay_free(), cyclic_two_class(), multiserver(),
                     with_delay(),        mm1(),              ldbcmp_boundary(),
                     all_delay_closed()};
    const std::vector<std::string> all = line::ba::list_valid_methods();
    for (NetD& m : models) {
        const SN& sn = m.get_struct();
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (all[i] == "ldbcmp.lower") continue;
            CHECK_MESSAGE(line::ba::method_degenerate(sn, all[i]).empty(), all[i]);
        }
    }
}
