/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Several forks in one model, and one fork nested inside another.
 *
 * `fj_mmt` used to refuse both by name ("the reference transform assumes one
 * fork per layer and its nesting bookkeeping is not ported"). Every number below
 * is MATLAB `SolverMVA(model).getAvg()` on the same model, which is the only
 * oracle that can settle these: the transform's own bookkeeping -- which fork an
 * auxiliary class belongs to, which fork is OUTER on a class, and which node's
 * visits measure a nested fork's firing rate -- has no closed form to check
 * against, so agreement with the reference is the specification.
 *
 * WHAT EACH MODEL ISOLATES.
 *
 *   twofork  Two forks in SERIES, neither inside the other. Both are outer, so
 *            both write their own synchronisation delay onto the original class,
 *            at their own join. This is the case that needs `fjforkmap`: without
 *            it the two auxiliary blocks are indistinguishable and the second
 *            fork's join is driven from the first fork's branch times.
 *   nestfork One fork INSIDE a branch of another. Only the outer fork writes the
 *            original class's delay; the inner join's delay is charged by
 *            `fj_find_paths` as it walks the enclosing branch, and the inner
 *            fork's firing rate is measured at its PARENT's node. This is the
 *            case that needs `fj_sort_forks`.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

mva::AvgResult<double> run(qn::Network<double>& m) {
    mva::MvaOptions opt;
    opt.method = "default";
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

/**
 * Delay -> Fork1 -{A1, A2}- Join1 -> Fork2 -{B1, B2}- Join2 -> Delay.
 * Station order: Delay, A1, A2, Join1, B1, B2, Join2.
 */
qn::Network<double> twofork() {
    qn::Network<double> m("twofork");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f1 = m.add_fork("Fork1");
    const std::size_t a1 = m.add_queue("A1", SchedStrategy::PS);
    const std::size_t a2 = m.add_queue("A2", SchedStrategy::PS);
    const std::size_t j1 = m.add_join("Join1", f1);
    const std::size_t f2 = m.add_fork("Fork2");
    const std::size_t b1 = m.add_queue("B1", SchedStrategy::PS);
    const std::size_t b2 = m.add_queue("B2", SchedStrategy::PS);
    const std::size_t j2 = m.add_join("Join2", f2);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(a1, c, D::exp_rate(2.0));
    m.set_service(a2, c, D::exp_rate(3.0));
    m.set_service(b1, c, D::exp_rate(4.0));
    m.set_service(b2, c, D::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f1, 1.0);
    P.set(f1, a1, 1.0);
    P.set(f1, a2, 1.0);
    P.set(a1, j1, 1.0);
    P.set(a2, j1, 1.0);
    P.set(j1, f2, 1.0);
    P.set(f2, b1, 1.0);
    P.set(f2, b2, 1.0);
    P.set(b1, j2, 1.0);
    P.set(b2, j2, 1.0);
    P.set(j2, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Delay -> Fork1 -{A1, (Fork2 -{B1,B2}- Join2)}- Join1 -> Delay.
 * Station order: Delay, A1, B1, B2, Join2, Join1.
 */
qn::Network<double> nestfork() {
    qn::Network<double> m("nestfork");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f1 = m.add_fork("Fork1");
    const std::size_t a1 = m.add_queue("A1", SchedStrategy::PS);
    const std::size_t f2 = m.add_fork("Fork2");
    const std::size_t b1 = m.add_queue("B1", SchedStrategy::PS);
    const std::size_t b2 = m.add_queue("B2", SchedStrategy::PS);
    const std::size_t j2 = m.add_join("Join2", f2);
    const std::size_t j1 = m.add_join("Join1", f1);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(a1, c, D::exp_rate(2.0));
    m.set_service(b1, c, D::exp_rate(4.0));
    m.set_service(b2, c, D::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f1, 1.0);
    P.set(f1, a1, 1.0);
    P.set(f1, f2, 1.0);
    P.set(f2, b1, 1.0);
    P.set(f2, b2, 1.0);
    P.set(b1, j2, 1.0);
    P.set(b2, j2, 1.0);
    P.set(a1, j1, 1.0);
    P.set(j2, j1, 1.0);
    P.set(j1, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fj_mmt: two forks each get their own record, join and auxiliary block") {
    qn::Network<double> m = twofork();
    const mva::FjMmt<double> tr = mva::fj_mmt(m.get_struct());

    REQUIRE(tr.forks.size() == 2);
    CHECK(tr.forks[0].joinNode != 0);
    CHECK(tr.forks[1].joinNode != 0);
    CHECK(tr.forks[0].joinStation != tr.forks[1].joinStation);
    // Both joins were converted, whichever fork claims them.
    CHECK(tr.joinStations.size() == 2);
    // Each fork forks the one class into two branches.
    CHECK(tr.forks[0].origfanout[1] == 2);
    CHECK(tr.forks[1].origfanout[1] == 2);

    // One class, one chain, two forks: one auxiliary class per fork, and
    // fjforkmap is what tells them apart -- fjclassmap alone maps both to class 1.
    REQUIRE(tr.auxclasses.size() == 2);
    CHECK(tr.fjclassmap[tr.auxclasses[0]] == 1);
    CHECK(tr.fjclassmap[tr.auxclasses[1]] == 1);
    CHECK(tr.fjforkmap[tr.auxclasses[0]] != tr.fjforkmap[tr.auxclasses[1]]);

    // Neither fork is inside the other, so both are outer on the class and both
    // are their own parent.
    CHECK(tr.forks[0].outer[1]);
    CHECK(tr.forks[1].outer[1]);
    CHECK(tr.forks[0].parent == 0);
    CHECK(tr.forks[1].parent == 1);

    // No Fork or Join survives the transform, which is what stops the inner
    // solver's own fork checks from firing.
    for (const auto& nd : tr.V.nodes) {
        CHECK(nd.nodetype != lang::NodeType::Fork);
        CHECK(nd.nodetype != lang::NodeType::Join);
    }
    CHECK(tr.V.fj.empty());
}

TEST_CASE("fj_mmt: a nested fork is marked inner and reparented onto its enclosure") {
    qn::Network<double> m = nestfork();
    const mva::FjMmt<double> tr = mva::fj_mmt(m.get_struct());

    REQUIRE(tr.forks.size() == 2);
    // Fork1 (node order: it is added first) encloses Fork2.
    const std::size_t outerIdx = 0, innerIdx = 1;
    CHECK(tr.forks[outerIdx].outer[1]);
    CHECK_FALSE(tr.forks[innerIdx].outer[1]);
    // The inner fork's firing rate is measured at its enclosure's node.
    CHECK(tr.forks[innerIdx].parent == outerIdx);
    CHECK(tr.forks[outerIdx].parent == outerIdx);
    // The outer fork splits into A1 and Fork2, so its fan-out is 2, not 3.
    CHECK(tr.forks[outerIdx].origfanout[1] == 2);
    CHECK(tr.forks[innerIdx].origfanout[1] == 2);
}

TEST_CASE("two forks in series agree with MATLAB SolverMVA") {
    qn::Network<double> m = twofork();
    const mva::AvgResult<double> r = run(m);
    REQUIRE(r.QN.rows() == 7);
    REQUIRE(r.QN.cols() == 1);

    // MATLAB SolverMVA(model).getAvg(), station order as built.
    const double Q[7] = {0.8578274838, 0.6892125082, 0.3798724274, 0.5792963439,
                         0.2613484643, 0.1991721546, 0.2344576433};
    const double U[7] = {0.8578272630, 0.4289136084, 0.2859424056, 0.0,
                         0.2144568042, 0.1715654434, 0.0};
    const double R[7] = {1.0000002574, 0.8034397774, 0.4428309337, 0.3376532782,
                         0.3046632925, 0.2321821349, 0.1366578482};
    const double Tp[7] = {0.8578272630, 0.8578272169, 0.8578272169, 0.8578272630,
                          0.8578272169, 0.8578272169, 0.8578272630};
    for (std::size_t i = 0; i < 7; ++i) {
        INFO("station ", i);
        CHECK(r.QN(i, 0) == doctest::Approx(Q[i]).epsilon(1e-7));
        CHECK(r.UN(i, 0) == doctest::Approx(U[i]).epsilon(1e-7));
        CHECK(r.RN(i, 0) == doctest::Approx(R[i]).epsilon(1e-7));
        CHECK(r.TN(i, 0) == doctest::Approx(Tp[i]).epsilon(1e-7));
    }

    // The two synchronisation delays are DIFFERENT, which is the whole point:
    // Join1 balances branches of rate 2 and 3, Join2 branches of rate 4 and 5, so
    // a transform that drove both joins from one fork's branch times would report
    // the same delay twice.
    CHECK(r.RN(3, 0) > r.RN(6, 0));
    // Both branches of both forks carry load; an unmodelled branch would be idle.
    for (std::size_t i : {1u, 2u, 4u, 5u}) CHECK(r.UN(i, 0) > 0.1);
    // Each join sees both of its branches arrive, so twice the class throughput.
    CHECK(r.AN(3, 0) == doctest::Approx(2.0 * r.TN(0, 0)).epsilon(1e-6));
    CHECK(r.AN(6, 0) == doctest::Approx(2.0 * r.TN(0, 0)).epsilon(1e-6));
}

TEST_CASE("a fork nested in a fork agrees with MATLAB SolverMVA") {
    qn::Network<double> m = nestfork();
    const mva::AvgResult<double> r = run(m);
    REQUIRE(r.QN.rows() == 6);

    // MATLAB SolverMVA(model).getAvg(): Delay, A1, B1, B2, Join2, Join1.
    const double Q[6] = {0.9976662120, 0.8908970454, 0.3263707969, 0.2444214671,
                         0.2912789044, 0.7409525311};
    const double U[6] = {0.9976660562, 0.4988330401, 0.2494164668, 0.1995331734, 0.0, 0.0};
    const double R[6] = {1.0000001562, 0.8929811918, 0.3271343721, 0.2449933141,
                         0.1167841459, 0.3713429502};
    const double Tp[6] = {0.9976660562, 0.9976660802, 0.9976658670, 0.9976658670,
                          0.9976660802, 0.9976660562};
    for (std::size_t i = 0; i < 6; ++i) {
        INFO("station ", i);
        CHECK(r.QN(i, 0) == doctest::Approx(Q[i]).epsilon(1e-7));
        CHECK(r.UN(i, 0) == doctest::Approx(U[i]).epsilon(1e-7));
        CHECK(r.RN(i, 0) == doctest::Approx(R[i]).epsilon(1e-7));
        CHECK(r.TN(i, 0) == doctest::Approx(Tp[i]).epsilon(1e-7));
    }

    // The outer join delays more than the inner one: it waits for A1 (mean 0.5)
    // against a branch that is itself a parallel pair, and the inner maximum is
    // already charged into that branch's time.
    CHECK(r.RN(5, 0) > r.RN(4, 0));
    // The inner join is an infinite server and is never utilized, but it DOES
    // hold jobs -- a nested join whose delay was never set would hold none.
    CHECK(r.UN(4, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(r.QN(4, 0) > 0.1);
}

TEST_CASE("a fork with no join and no Sink is refused by the model layer") {
    // THE REFUSAL IS NARROW, and the narrowness is the point: `fj_nojoin` ships
    // in all three reference suites as an OPEN model whose fork branches each
    // end at the Sink, and MATLAB and native Python both solve it. What is
    // refused is a join-less fork whose siblings can neither merge nor depart,
    // as in the CLOSED model below, where the population is not conserved.
    //
    // NOT a parity gap, though the reference's fjFixedPoint does carry a
    // join-less branch. MATLAB accepts the link and then dies inside the solve:
    // `sortForks` calls `nestedForks(f, [])`, whose `startNode == endNode` test
    // cannot hold against an empty join, so on a closed model it recurses until
    // MATLAB reports "Out of memory. The likely cause is an infinite recursion"
    // (measured 2026-07-30 on exactly this model). Refusing by name at the model
    // layer is the better answer, so this pins the refusal rather than a number.
    qn::Network<double> m("nojoin");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);

    m.link(P);
    // The refusal is in the struct refresh, not in link(): the fork/join pairing
    // is a property of `sn.fj`, which link() does not build.
    std::string what;
    bool threw = false;
    try {
        m.get_struct();
    } catch (const line::Error& e) {
        threw = true;
        what = e.what();
    }
    CHECK(threw);
    // The refusal must name the Fork, so the caller knows which node to close.
    CHECK(what.find("Fork") != std::string::npos);
    CHECK(what.find("Join") != std::string::npos);
}

TEST_CASE("a fork with no join whose branches reach the Sink is accepted") {
    // `fj_nojoin` of the reference suites: Source -> Fork -> {Q1, Q2, Q3} ->
    // Sink. The siblings never merge, but they DO depart, so the model is well
    // formed and MATLAB, native Python and this port all build it.
    qn::Network<double> m("fj_nojoin");
    const std::size_t src = m.add_source("Source");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("class1");
    m.set_arrival(src, c, D::exp_rate(0.5));
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, snk, 1.0);
    P.set(q2, snk, 1.0);
    m.link(P);
    CHECK_NOTHROW(m.get_struct());
    CHECK(m.get_struct().fj.empty());
}
