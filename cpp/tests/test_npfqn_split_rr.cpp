/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Deterministic (round-robin) traffic split degrees, and the two analyzers that
 * consume them.
 *
 * The degree matrix is checked against the four topologies
 * `npfqn_traffic_split_rr.m` distinguishes: no dispatcher at all, a station
 * dispatching itself, a station feeding a dedicated round-robin router, and a
 * router shared by two upstream streams, which interleaves them and so admits
 * no deterministic rule.
 *
 * The analyzer rows pin the CONSEQUENCE: a one-in-k dispatch splits a stream
 * into flows that are LESS variable than Bernoulli thinning gives, so the
 * downstream queue is shorter. A model whose degrees are all 1 must reproduce
 * the numbers the port produced before the correction existed.
 */

#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_traffic_split_rr.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;
using qn::RoutingStrategy;

namespace {

using D = Distrib<double>;

/**
 * Source -> Q1 -> {Q2, Q3} -> Sink, with the split at Q1 governed by `rs`.
 * At PROB the two branches carry one half each, which is also what the uniform
 * expansion of RROBIN produces, so the two models differ ONLY in the split
 * degree the traffic equations see.
 */
qn::Network<double> fork_out(RoutingStrategy rs) {
    qn::Network<double> m("rr");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, D::exp_rate(1.0));
    m.set_service(q1, o, D::erlang(4.0, 2));  // scv 1/2, a non-renewal split
    m.set_service(q2, o, D::exp_rate(2.0));
    m.set_service(q3, o, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 0.5);
    P.set(q1, q3, 0.5);
    P.set(q2, snk, 1.0);
    P.set(q3, snk, 1.0);
    m.link(P);
    m.set_routing(q1, o, rs);
    return m;
}

}  // namespace

TEST_CASE("npfqn_traffic_split_rr is all ones without a round-robin dispatcher") {
    qn::Network<double> m = fork_out(RoutingStrategy::PROB);
    const Matrix<double> kRR = npfqn::npfqn_traffic_split_rr(m.get_struct());
    CHECK(kRR.rows() == m.get_struct().nstations);
    CHECK(kRR.cols() == 1u);
    for (std::size_t i = 0; i < kRR.rows(); ++i) CHECK(kRR(i, 0) == doctest::Approx(1.0));
}

TEST_CASE("a station that dispatches round-robin has degree equal to its outgoing links") {
    qn::Network<double> m = fork_out(RoutingStrategy::RROBIN);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const Matrix<double> kRR = npfqn::npfqn_traffic_split_rr(sn);
    // stations are Source, Q1, Q2, Q3; only Q1 dispatches, over two links
    CHECK(kRR(0, 0) == doctest::Approx(1.0));
    CHECK(kRR(1, 0) == doctest::Approx(2.0));
    CHECK(kRR(2, 0) == doctest::Approx(1.0));
    CHECK(kRR(3, 0) == doctest::Approx(1.0));
    // the expansion is uniform: round robin visits each link equally often
    const std::size_t K = sn.nclasses, q1 = 2;  // 1-based node index of Q1
    CHECK(sn.rtnodes((q1 - 1) * K, (3 - 1) * K) == doctest::Approx(0.5));
    CHECK(sn.rtnodes((q1 - 1) * K, (4 - 1) * K) == doctest::Approx(0.5));
}

TEST_CASE("a dedicated round-robin router lends its degree to the station upstream") {
    qn::Network<double> m("rr-router");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t rtr = m.add_router("R");
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t q4 = m.add_queue("Q4", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, D::exp_rate(1.0));
    for (std::size_t q : {q1, q2, q3, q4}) m.set_service(q, o, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, rtr, 1.0);
    P.set(rtr, q2, 1.0 / 3.0);
    P.set(rtr, q3, 1.0 / 3.0);
    P.set(rtr, q4, 1.0 / 3.0);
    P.set(q2, snk, 1.0);
    P.set(q3, snk, 1.0);
    P.set(q4, snk, 1.0);
    m.link(P);
    m.set_routing(rtr, o, RoutingStrategy::RROBIN);
    const Matrix<double> kRR = npfqn::npfqn_traffic_split_rr(m.get_struct());
    // Q1 feeds the router with probability one and is its only feed, so its
    // departures are the ones dispatched one-in-three
    CHECK(kRR(1, 0) == doctest::Approx(3.0));
    CHECK(kRR(0, 0) == doctest::Approx(1.0));
    CHECK(kRR(2, 0) == doctest::Approx(1.0));
}

TEST_CASE("a round-robin router shared by two streams admits no deterministic split") {
    qn::Network<double> m("rr-shared");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t rtr = m.add_router("R");
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t q4 = m.add_queue("Q4", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, D::exp_rate(1.0));
    for (std::size_t q : {q1, q2, q3, q4}) m.set_service(q, o, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 0.5);
    P.set(src, q2, 0.5);
    P.set(q1, rtr, 1.0);
    P.set(q2, rtr, 1.0);
    P.set(rtr, q3, 0.5);
    P.set(rtr, q4, 0.5);
    P.set(q3, snk, 1.0);
    P.set(q4, snk, 1.0);
    m.link(P);
    m.set_routing(rtr, o, RoutingStrategy::RROBIN);
    const Matrix<double> kRR = npfqn::npfqn_traffic_split_rr(m.get_struct());
    // the pointer is advanced by both streams, so neither sees a clean one-in-2
    for (std::size_t i = 0; i < kRR.rows(); ++i) CHECK(kRR(i, 0) == doctest::Approx(1.0));
}

TEST_CASE("QNA reports shorter queues under round-robin than under a probabilistic split") {
    qn::Network<double> mp = fork_out(RoutingStrategy::PROB);
    qn::Network<double> mr = fork_out(RoutingStrategy::RROBIN);
    mva::MvaOptions opt;
    opt.method = "qna";
    Matrix<double> init;
    const mva::AvgResult<double> rp = mva::solver_mva_run_analyzer(mp.get_struct(), opt, init);
    const mva::AvgResult<double> rr = mva::solver_mva_run_analyzer(mr.get_struct(), opt, init);
    CHECK(rp.actualmethod == "qna");
    CHECK(rr.actualmethod == "qna");
    // MATLAB SolverMVA(model,'method','qna').getAvgQLen(), stations in order
    // Source, Q1, Q2, Q3. Q1 is upstream of the split and is unchanged; the two
    // branches drop from 1/3 to 5/16 because a one-in-2 dispatch feeds them a
    // less variable stream than a fair coin does.
    CHECK(rp.QN(1, 0) == doctest::Approx(0.874999998125).epsilon(1e-9));
    CHECK(rp.QN(2, 0) == doctest::Approx(0.333333333194).epsilon(1e-9));
    CHECK(rp.QN(3, 0) == doctest::Approx(0.333333333194).epsilon(1e-9));
    CHECK(rr.QN(1, 0) == doctest::Approx(0.874999998125).epsilon(1e-9));
    CHECK(rr.QN(2, 0) == doctest::Approx(0.312500000104).epsilon(1e-9));
    CHECK(rr.QN(3, 0) == doctest::Approx(0.312500000104).epsilon(1e-9));
    // the dispatcher does not change the rate a branch receives, only the
    // variability of its arrival stream
    CHECK(rr.TN(2, 0) == doctest::Approx(rp.TN(2, 0)).epsilon(1e-9));
}

TEST_CASE("MNA resolves the same split and refuses it on a closed model") {
    qn::Network<double> mp = fork_out(RoutingStrategy::PROB);
    qn::Network<double> mr = fork_out(RoutingStrategy::RROBIN);
    mam::MamOptions opt;
    opt.method = "mna";
    const mva::AvgResult<double> rp = mam::solver_mam_run_analyzer(mp.get_struct(), opt);
    const mva::AvgResult<double> rr = mam::solver_mam_run_analyzer(mr.get_struct(), opt);
    // MATLAB SolverMAM(model,'method','mna').getAvgQLen(). The branch queues sit
    // below their QNA values because the station solve is MMAP/PH/1 rather than
    // a Whitt formula; the round-robin correction moves them the same way.
    CHECK(rp.QN(2, 0) == doctest::Approx(0.329037112293).epsilon(1e-9));
    CHECK(rp.QN(3, 0) == doctest::Approx(0.329037112293).epsilon(1e-9));
    CHECK(rr.QN(2, 0) == doctest::Approx(0.289637243717).epsilon(1e-9));
    CHECK(rr.QN(3, 0) == doctest::Approx(0.289637243717).epsilon(1e-9));

    // closed: solver_mna_closed has no round-robin correction, so the model is
    // refused by name rather than solved as if the dispatcher were random
    qn::Network<double> mc("rr-closed");
    const std::size_t d = mc.add_delay("Think");
    const std::size_t q1 = mc.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = mc.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = mc.add_closed_class("C1", 3.0, d);
    mc.set_service(d, c, D::exp_rate(1.0));
    mc.set_service(q1, c, D::exp_rate(2.0));
    mc.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 0.5);
    P.set(d, q2, 0.5);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    mc.link(P);
    mc.set_routing(d, c, RoutingStrategy::RROBIN);
    CHECK_THROWS_AS(mam::solver_mam_run_analyzer(mc.get_struct(), opt), UnsupportedError);
}
