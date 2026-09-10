/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * sn_to_qrf_blocking, the derivation of the QRF BAS blocking tables from an sn.
 *
 * The tables used to be hand-built by the caller and the solver refused a model
 * without them, which made `qrf.bas` unreachable from a plain
 * `SolverBA(model,'qrf.bas')`. Everything the enumeration needs is implied by
 * the model, so the expected values below are not a record of what this port
 * happens to produce: they are what the blocking structure of each model
 * forces, checked against the reference tables for cqn_bas_blocking.
 */

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_to_qrf_blocking.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ba/solver_ba_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** cqn_bas_blocking: Queue1 -BAS-> Queue2(cap 1), N = 2. */
qn::Network<double> model_cqn_bas_blocking() {
    qn::Network<double> m("cqn_bas_blocking");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2.0, q1);
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(0.8));
    m.set_capacity(q2, 1.0);
    m.set_drop_rule(q1, c, qn::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    return m;
}

/** Two feeders into one capped queue: Q1 -> {Q2, Q3(cap 1)}, Q2 -> Q3, Q3 -> Q1. */
qn::Network<double> model_two_feeders(double N, double cap3) {
    qn::Network<double> m("twoFeeders");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", N, q1);
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(0.9));
    m.set_service(q3, c, D::exp_rate(1.2));
    m.set_capacity(q3, cap3);
    m.set_drop_rule(q1, c, qn::DropStrategy::BAS);
    m.set_drop_rule(q2, c, qn::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 0.5);
    P.set(q1, q3, 0.5);
    P.set(q2, q3, 1.0);
    P.set(q3, q1, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("sn_to_qrf_blocking reproduces the cqn_bas_blocking tables") {
    qn::Network<double> m = model_cqn_bas_blocking();
    const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(m.get_struct(), 2);

    REQUIRE(blk.msg.empty());
    CHECK(blk.f == 2);   // Queue2 is the one binding buffer
    CHECK(blk.MR == 2);  // empty configuration, plus Queue1 blocked
    CHECK(blk.ZM == 1);
    REQUIRE(blk.blockers.size() == 1);
    CHECK(blk.blockers[0] == 1);

    REQUIRE(blk.F.size() == 2);
    CHECK(blk.F[0] == 2);  // unbounded: capped by the population
    CHECK(blk.F[1] == 1);

    // Invariant 1: configuration 1 is the empty one, which ZERO4/5/7/8 assume.
    CHECK(blk.ZZ[0] == 0);
    CHECK(blk.BB[0][0] == 0);
    CHECK(blk.BB[0][1] == 0);
    CHECK(blk.MM[0][0] == 0);

    // Configuration 2: Queue1 blocked, and it is the head that takes the slot.
    CHECK(blk.ZZ[1] == 1);
    CHECK(blk.BB[1][0] == 1);
    CHECK(blk.BB[1][1] == 0);
    CHECK(blk.MM[1][0] == 1);

    // MM1 is indexed by the queue that BECOMES blocked, not by f: adding
    // Queue1 to the empty configuration reaches configuration 2.
    CHECK(blk.MM1[0][0] == 2);
    CHECK(blk.MM1[0][1] == 0);
    CHECK(blk.MM1[1][0] == 0);
}

TEST_CASE("sn_to_qrf_blocking enumerates every order with two feeders") {
    qn::Network<double> m = model_two_feeders(4.0, 1.0);
    const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(m.get_struct(), 3);

    REQUIRE(blk.msg.empty());
    CHECK(blk.f == 3);
    CHECK(blk.ZM == 2);
    // 1 empty + 2 singletons + 2 orders of {Q1,Q2}
    CHECK(blk.MR == 5);
    REQUIRE(blk.blockers.size() == 2);

    for (int m_i = 0; m_i < blk.MR; ++m_i) {
        int depth = 0;
        for (std::size_t i = 0; i < blk.BB[m_i].size(); ++i) depth += blk.BB[m_i][i];
        CHECK(blk.ZZ[m_i] == depth);
        // f is never one of the blocked queues
        CHECK(blk.BB[m_i][blk.f - 1] == 0);
    }

    // Invariant 3: the head is preserved along a successor edge, so a
    // configuration below the maximum depth reaches a deeper one for every
    // feeder not already blocked, and that successor keeps MM(m,0).
    for (int m_i = 0; m_i < blk.MR; ++m_i) {
        if (blk.ZZ[m_i] >= blk.ZM) continue;
        for (std::size_t b = 0; b < blk.blockers.size(); ++b) {
            const int j = blk.blockers[b];
            if (blk.BB[m_i][j - 1]) continue;
            const int mp = blk.MM1[m_i][j - 1];
            REQUIRE(mp >= 1);
            CHECK(blk.ZZ[mp - 1] == blk.ZZ[m_i] + 1);
            CHECK(blk.BB[mp - 1][j - 1] == 1);
            if (blk.ZZ[m_i] > 0) CHECK(blk.MM[mp - 1][0] == blk.MM[m_i][0]);
        }
    }
}

TEST_CASE("sn_to_qrf_blocking refuses what qrf.bas cannot express") {
    SUBCASE("more than one binding buffer") {
        qn::Network<double> m = model_two_feeders(4.0, 1.0);
        qn::NetworkStruct<double> sn = m.get_struct();
        sn.stations[1].cap = 1.0;  // a second binding buffer
        sn.classcap[1][0] = 1.0;
        const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(sn, 3);
        CHECK_FALSE(blk.msg.empty());
        CHECK(blk.msg.find("single finite-capacity queue") != std::string::npos);
        CHECK(blk.msg.find("qrf.rsrd") != std::string::npos);
    }

    SUBCASE("an oversized enumeration is refused, never truncated") {
        qn::Network<double> m = model_two_feeders(4.0, 1.0);
        const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(m.get_struct(), 3, 10.0);
        CHECK_FALSE(blk.msg.empty());
        CHECK(blk.msg.find("cannot be truncated") != std::string::npos);
    }
}

TEST_CASE("sn_to_qrf_blocking reports no blocking state when none is reachable") {
    // f holds 3 of the 3 jobs at capacity, so nothing can ever be held behind it.
    qn::Network<double> m = model_two_feeders(3.0, 3.0);
    const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(m.get_struct(), 3);
    REQUIRE(blk.msg.empty());
    CHECK(blk.ZM == 0);
    CHECK(blk.MR == 1);
    CHECK(blk.ZZ[0] == 0);
}

TEST_CASE("a blocked model routes 'default' to the QRF BAS bound") {
    qn::Network<double> m = model_cqn_bas_blocking();
    const qn::NetworkStruct<double> sn = m.get_struct();

    // The routing itself, and its one-line reason when it does not apply.
    const std::pair<std::string, std::string> routed = ba::blocking_default(sn);
    CHECK(routed.first == "qrf.bas");
    CHECK(routed.second.empty());

    // 'default' used to be refused outright on this model: gb.upper is
    // blocking-blind, and every other family was dropped from the list.
    ba::BaOptions opt;
    opt.method = "default";
    const mva::AvgResult<double> r = ba::solver_ba_run_analyzer(sn, opt);
    REQUIRE(r.UN.rows() == 2);
    // The exact utilizations of the 3-state chain. Queue1 SERVES in states
    // (2,0) and (1,1) and is BLOCKED in the third, holding a job it has already
    // finished, so pi0 + pi1 = 0.590164; Queue2 is busy in the last two,
    // pi1 + pi2 = 0.737705.
    //
    // Queue1 used to read a vacuous 1.0 here, because the objective summed p2
    // over every blocking configuration -- P(n >= 1), which counts a blocked
    // server as busy. The objective is the e variables now, restricted by UEFF
    // to the configurations where the station is not blocked.
    CHECK(r.UN(0, 0) == doctest::Approx(0.590164).epsilon(1e-4));
    CHECK(r.UN(1, 0) == doctest::Approx(0.737705).epsilon(1e-4));

    // And the model's own default is offered back by the method list. Every
    // OTHER survivor must still be a blocking bound: 'default' is the one entry
    // whose resolved form ('gb.upper') is blind, which is precisely why the
    // runner rewrites it rather than dispatching it.
    const std::vector<std::string> valid = ba::list_valid_methods(sn);
    CHECK(std::find(valid.begin(), valid.end(), "default") != valid.end());
    for (std::size_t i = 0; i < valid.size(); ++i) {
        if (valid[i] == "default") continue;
        CHECK_FALSE(ba::ignores_blocking(ba::resolve_method(valid[i])));
    }
}

TEST_CASE("a blocked model outside the QRF shape still refuses, and says why") {
    // A multiserver capped queue: qrf.bas models every station as one server.
    qn::Network<double> m = model_two_feeders(4.0, 1.0);
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.stations[2].nservers = 2.0;

    const std::pair<std::string, std::string> routed = ba::blocking_default(sn);
    CHECK(routed.first.empty());
    CHECK(routed.second.find("multiserver") != std::string::npos);

    ba::BaOptions opt;
    opt.method = "default";
    CHECK_THROWS_AS(ba::solver_ba_run_analyzer(sn, opt), UnsupportedError);
}

TEST_CASE("getBounds refuses the routed default as one-sided, not blind") {
    qn::Network<double> m = model_cqn_bas_blocking();
    ba::BaOptions opt;
    opt.method = "default";
    // qrf.bas is solved in the 'max' direction alone, so there is no bracket --
    // and saying "does not support blocking" here would contradict the run.
    CHECK_THROWS_AS(ba::ba_bounds(m.get_struct(), opt), UnsupportedError);
}
