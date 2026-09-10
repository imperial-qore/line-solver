/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The SolverMVA class surface: the method gate, the two conversions and the
 * metric filter that sit around one dispatch.
 *
 * The numbers are MATLAB's `[Q,U,R,T,A,W] = SolverMVA(model).getAvg()`, which
 * is the right comparison for this layer: the arrival rate and the residence
 * time are computed HERE and not by any analyzer, so an analyzer-level test
 * cannot see them.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

mva::AvgResult<double> run(qn::Network<double>& m, const std::string& method = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

qn::Network<double> closed_model() {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("the runner reproduces every column MATLAB's getAvg returns") {
    SUBCASE("a closed model") {
        qn::Network<double> m = closed_model();
        const mva::AvgResult<double> r = run(m);
        CHECK(r.actualmethod == "exact");
        CHECK(r.QN(0, 0) == doctest::Approx(1.578947368).epsilon(1e-9));
        CHECK(r.QN(1, 0) == doctest::Approx(1.421052632).epsilon(1e-9));
        CHECK(r.UN(1, 0) == doctest::Approx(0.789473684).epsilon(1e-9));
        CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-9));
        // residence time, computed by the runner and not by the analyzer
        CHECK(r.WN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.WN(1, 0) == doctest::Approx(0.9).epsilon(1e-9));
        // arrival rate, from the class-expanded routing sn.rt
        CHECK(r.AN(1, 0) == doctest::Approx(1.578947368).epsilon(1e-9));
        CHECK(r.TN(1, 0) == doctest::Approx(1.578947368).epsilon(1e-9));
    }
    SUBCASE("an open model, where the Source has no arrivals of its own") {
        qn::Network<double> m("mm1");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t o = m.add_open_class("O1");
        m.set_arrival(s, o, D::exp_rate(1.0));
        m.set_service(q, o, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(s, q, 1.0);
        P.set(q, k, 1.0);
        m.link(P);
        const mva::AvgResult<double> r = run(m);
        CHECK(r.actualmethod == "mm1");
        CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(r.RN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.WN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.AN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.TN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        // a Source is not fed by anything
        CHECK(r.AN(0, 0) == doctest::Approx(0.0));
    }
}

TEST_CASE("the method whitelist is the solver's, not the analyzer's") {
    qn::Network<double> m = closed_model();
    const std::vector<std::string> valid = mva::list_valid_methods(m.get_struct());
    auto has = [&](const std::string& s) {
        return std::find(valid.begin(), valid.end(), s) != valid.end();
    };
    CHECK(has("default"));
    CHECK(has("exact"));
    CHECK(has("lin"));
    CHECK(has("schmidt-ext"));
    // marie is listed for a closed model whose classes are routed alike
    CHECK(has("marie"));
    // qna and rqna are open-network analyzers and are not advertised here
    CHECK_FALSE(has("qna"));
    CHECK_FALSE(has("rqna"));
    // aql, qsa and tay are listed for a single-server model, as SolverMVA.m:56
    // does, and each is refused by name once a station gains servers
    CHECK(has("aql"));
    CHECK(has("amva.aql"));
    CHECK(has("tay"));
    CHECK_NOTHROW(run(m, "aql"));
    CHECK_THROWS_AS(run(m, "nosuchmethod"), UnsupportedError);
}

TEST_CASE("a Finite Capacity Region is refused, not silently dropped") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    // A region does not need a binding cap for the gate to fire, since MVA
    // cannot enforce ANY region regardless of its shape.
    m.add_region({q}, {-1.0}, 5.0);
    CHECK_THROWS_AS(run(m), UnsupportedError);
    try {
        run(m);
        FAIL("expected UnsupportedError");
    } catch (const UnsupportedError& e) {
        const std::string what = e.what();
        CHECK(what.find("Finite Capacity Region") != std::string::npos);
    }
}

// Immediate feedback is the one construct the reference APPROXIMATES rather
// than refusing: a fed-back job keeps its server, which no mean-value argument
// can express, so MVA and NC warn and count the visit as an ordinary re-entry
// (@SolverMVA/runAnalyzer.m:26-27). The failure this pins is a silent one --
// the numbers are unchanged and only the warning says they are approximate, so
// a dropped warning leaves the table looking authoritative.
TEST_CASE("immediate feedback is warned about, not refused") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const mva::AvgResult<double> before = run(m);
    CHECK(before.warning.empty());
    CHECK_FALSE(m.get_struct().has_immediate_feedback());

    m.set_immediate_feedback(q, c);
    CHECK(m.get_struct().has_immediate_feedback());

    const mva::AvgResult<double> after = run(m);
    CHECK(after.warning.find("immediate feedback") != std::string::npos);
    // The reference solves; the numbers are the same table it would print
    // without the flag, which is exactly what the warning is warning about.
    CHECK(after.QN(0, 0) == doctest::Approx(before.QN(0, 0)));

    SUBCASE("the class-wide spelling reaches the same matrix") {
        qn::Network<double> n("cqn2");
        const std::size_t d2 = n.add_delay("Delay");
        const std::size_t q2 = n.add_queue("Queue", SchedStrategy::FCFS);
        const std::size_t c2 = n.add_closed_class("C1", 3.0, d2);
        n.set_service(d2, c2, D::exp_rate(1.0));
        n.set_service(q2, c2, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P2;
        P2.set(d2, q2, 1.0);
        P2.set(q2, d2, 1.0);
        n.link(P2);
        n.set_class_immediate_feedback(c2);
        // sn.immfeed is the OR, so a class-level flag marks EVERY station
        CHECK(n.get_struct().immfeed[0][0]);
        CHECK(n.get_struct().immfeed[1][0]);
        CHECK(run(n).warning.find("immediate feedback") != std::string::npos);
    }
}

TEST_CASE("an open single-class two-station model additionally lists the closed forms") {
    qn::Network<double> m("mm1");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, D::exp_rate(1.0));
    m.set_service(q, o, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const std::vector<std::string> valid = mva::list_valid_methods(m.get_struct());
    auto has = [&](const std::string& s2) {
        return std::find(valid.begin(), valid.end(), s2) != valid.end();
    };
    CHECK(has("mm1"));
    CHECK(has("gig1.klb"));
    CHECK(has("qna"));
    CHECK(has("rqna"));
    // and marie, a closed-network method, is withheld
    CHECK_FALSE(has("marie"));
}

TEST_CASE("a fork-join model runs the shared fixed point, not the bare dispatch") {
    // Two branches with different rates, so the join actually synchronises: the
    // Fork becomes a router, the Join a delay carrying E[max] - mean, and the
    // branches the circulating job did not take are carried by the auxiliary
    // open class the transform mints. Numbers are MATLAB's SolverMVA getAvg.
    qn::Network<double> m("fj");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    const mva::AvgResult<double> r = run(m);
    // the transform's auxiliary columns are merged back, so the result is over
    // the ONE original class and the FOUR original stations
    CHECK(r.QN.rows() == 4);
    CHECK(r.QN.cols() == 1);
    // Q sums to 3.078 and not to the population 2: a forked job is present on
    // every branch at once, so the queue lengths of a fork-join model are not a
    // partition of N. The reference reports the same.
    CHECK(r.QN(0, 0) == doctest::Approx(0.992220497).epsilon(1e-8));
    CHECK(r.QN(1, 0) == doctest::Approx(0.882210744).epsilon(1e-8));
    CHECK(r.QN(2, 0) == doctest::Approx(0.465068005).epsilon(1e-8));
    CHECK(r.QN(3, 0) == doctest::Approx(0.738217004).epsilon(1e-8));
    // both branches are loaded: without the auxiliary open stream the branch
    // the circulating job did not take would carry no load at all
    CHECK(r.UN(1, 0) == doctest::Approx(0.496110097).epsilon(1e-8));
    CHECK(r.UN(2, 0) == doctest::Approx(0.330740065).epsilon(1e-8));
    // the Join is an infinite server, so it holds jobs but is never utilized
    CHECK(r.UN(3, 0) == doctest::Approx(0.0).epsilon(1e-8));
    // and its response time IS the synchronisation delay the fixed point set
    CHECK(r.RN(3, 0) == doctest::Approx(0.372002610).epsilon(1e-8));
    CHECK(r.RN(1, 0) == doctest::Approx(0.889127986).epsilon(1e-8));
    CHECK(r.TN(0, 0) == doctest::Approx(0.992220265).epsilon(1e-8));
    // the Join sees both branches arrive, so twice the class throughput
    CHECK(r.AN(3, 0) == doctest::Approx(1.984440390).epsilon(1e-8));
}

}  // namespace
