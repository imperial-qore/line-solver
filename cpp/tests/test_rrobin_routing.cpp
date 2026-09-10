/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ssa/ssa_dispatch.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;
using lang::RoutingStrategy;
using D = Distrib<double>;

namespace {
qn::Network<double> rr_model(RoutingStrategy rs) {
    qn::Network<double> m("rr");
    const std::size_t d = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2, d, 0);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q1, 0.5);
    P.set(c, c, d, q2, 0.5);
    P.set(c, c, q1, d, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);
    m.set_routing(d, c, rs);
    return m;
}
}  // namespace

TEST_CASE("rrobin: the pointer is a state coordinate") {
    qn::Network<double> m = rr_model(RoutingStrategy::RROBIN);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.rr_var_slot(1, 1) == 1u);
    const std::vector<std::size_t> ol = sn.rr_outlinks(1, 1);
    REQUIRE(ol.size() == 2u);
    CHECK(ol[0] == 2u);
    CHECK(ol[1] == 3u);
    std::vector<double> var(1, 2.0);
    sn.rr_advance(1, 1, var);
    CHECK(var[0] == doctest::Approx(3.0));
    sn.rr_advance(1, 1, var);
    CHECK(var[0] == doctest::Approx(2.0));
}

TEST_CASE("rrobin: CTMC matches MATLAB and differs from RAND") {
    ctmc::CtmcOptions opt;
    const mva::AvgResult<double> rr =
        ctmc::solver_ctmc_run_analyzer(rr_model(RoutingStrategy::RROBIN).get_struct(), opt);
    // MATLAB SolverCTMC on the same model
    CHECK(rr.QN(0, 0) == doctest::Approx(1.3933649289099526).epsilon(1e-9));
    CHECK(rr.QN(1, 0) == doctest::Approx(0.36966824644549762).epsilon(1e-9));
    CHECK(rr.QN(2, 0) == doctest::Approx(0.23696682464454977).epsilon(1e-9));
    CHECK(rr.TN(1, 0) == doctest::Approx(0.69668246445497628).epsilon(1e-9));
    CHECK(rr.TN(2, 0) == doctest::Approx(0.69668246445497628).epsilon(1e-9));

    const mva::AvgResult<double> rd =
        ctmc::solver_ctmc_run_analyzer(rr_model(RoutingStrategy::RAND).get_struct(), opt);
    CHECK(rd.QN(0, 0) == doctest::Approx(1.3509933774834435).epsilon(1e-9));
    // the dispatcher is not a coin: the two must not agree
    CHECK(std::abs(rr.QN(0, 0) - rd.QN(0, 0)) > 1e-6);
}

TEST_CASE("wrrobin: the weighted cycle, and the 2:1 split it forces") {
    qn::Network<double> m("wrr");
    const std::size_t d = m.add_delay("D");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2, d, 0);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q1, 0.5);
    P.set(c, c, d, q2, 0.5);
    P.set(c, c, q1, d, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);
    m.set_routing(d, c, RoutingStrategy::WRROBIN);
    std::map<std::size_t, double> w;
    w[q1] = 2.0;
    w[q2] = 1.0;
    m.set_routing_weights(d, c, w);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    // MATLAB's weighted_outlinks on the same model is [2 2 3]
    const std::vector<std::size_t> cy = sn.rr_weighted_outlinks(1, 1);
    REQUIRE(cy.size() == 3u);
    CHECK(cy[0] == 2u);
    CHECK(cy[1] == 2u);
    CHECK(cy[2] == 3u);
    // the slot holds a POSITION, so it advances 1 -> 2 -> 3 -> 1
    std::vector<double> var(1, 1.0);
    sn.rr_advance(1, 1, var);
    CHECK(var[0] == doctest::Approx(2.0));
    sn.rr_advance(1, 1, var);
    CHECK(var[0] == doctest::Approx(3.0));
    sn.rr_advance(1, 1, var);
    CHECK(var[0] == doctest::Approx(1.0));

    ctmc::CtmcOptions opt;
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, opt);
    CHECK(r.QN(0, 0) == doctest::Approx(1.3273237771538349).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(0.52490025121915163).epsilon(1e-9));
    CHECK(r.QN(2, 0) == doctest::Approx(0.14777597162701345).epsilon(1e-9));
    // the weights ARE the split: twice the throughput at Q1
    CHECK(r.TN(1, 0) == doctest::Approx(2.0 * r.TN(2, 0)).epsilon(1e-9));
}

TEST_CASE("rrobin: a Router NODE dispatches the same way a station does") {
    qn::Network<double> m("rrnode");
    const std::size_t d = m.add_delay("D");
    const std::size_t rt = m.add_router("R");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2, d, 0);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, rt, 1.0);
    P.set(c, c, rt, q1, 0.5);
    P.set(c, c, rt, q2, 0.5);
    P.set(c, c, q1, d, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);
    m.set_routing(rt, c, RoutingStrategy::RROBIN);

    ctmc::CtmcOptions opt;
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);
    // MATLAB gives the SAME numbers as the station-dispatch model: a Router is
    // transparent, so where the pointer lives cannot change the answer.
    //
    // 1e-7 and not 1e-9 because a Router holds the job for 1/Immediate rather
    // than for no time at all, and that residue is O(1e-8) here. The reference
    // carries the same finite rate, so the agreement is exactly as close as the
    // two idealizations are.
    CHECK(r.QN(0, 0) == doctest::Approx(1.3933649289099526).epsilon(1e-7));
    CHECK(r.QN(1, 0) == doctest::Approx(0.36966824644549762).epsilon(1e-7));
    CHECK(r.QN(2, 0) == doctest::Approx(0.23696682464454977).epsilon(1e-7));
    CHECK(r.TN(1, 0) == doctest::Approx(r.TN(2, 0)).epsilon(1e-9));
}

TEST_CASE("rrobin: the SSA serial engine reads the same pointer") {
    // The serial engine needs no cursor of its own: the pointer is a coordinate
    // of the state it walks, so the sample path must converge on the exact CTMC
    // answer rather than on the uniform-split one.
    qn::Network<double> mrr = rr_model(RoutingStrategy::RROBIN);
    ctmc::CtmcOptions copt;
    const mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(mrr.get_struct(), copt);

    ssa::SsaOptions sopt;
    sopt.samples = 200000;
    sopt.seed = 23000;
    sopt.method = "serial";
    const ssa::SsaSolution sim = ssa::solver_ssa(mrr.get_struct(), sopt);
    // 3 per cent of the exact value: a sample path of this length on a two-job
    // closed chain, not a claim about the estimator's variance.
    CHECK(sim.QN(1, 0) == doctest::Approx(exact.QN(1, 0)).epsilon(0.03));
    CHECK(sim.QN(2, 0) == doctest::Approx(exact.QN(2, 0)).epsilon(0.03));
    // the dispatcher forces equal throughputs; the uniform split would too, so
    // the discriminating check is the queue lengths above
    CHECK(sim.TN(1, 0) == doctest::Approx(sim.TN(2, 0)).epsilon(0.03));
}

TEST_CASE("ctmc getters: getStartRate and getPreemptRate reach the caller") {
    // The two derived rates were computed by solver_ctmc_avg_from_pi and
    // reachable through nothing: solver_ctmc_run_analyzer returns an mva::AvgResult,
    // which has no slot for them. MATLAB exposes them as plain getters.
    qn::Network<double> m("sp");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 2, d, 0);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ctmc::CtmcOptions opt;
    const Matrix<double> St = ctmc::ctmc_get_start_rate(m.get_struct(), opt);
    const Matrix<double> Pr = ctmc::ctmc_get_preempt_rate(m.get_struct(), opt);
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);
    // MATLAB SolverCTMC on the same model
    CHECK(St(0, 0) == doctest::Approx(1.2).epsilon(1e-12));
    CHECK(St(1, 0) == doctest::Approx(1.2).epsilon(1e-12));
    CHECK(Pr(0, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(Pr(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    // the identity the two exist to be checked against: at a lossless station
    // with no in-service abandonment, every job starts service once per entry
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(St(i, 0) == doctest::Approx(r.TN(i, 0) + Pr(i, 0)).epsilon(1e-9));
}
