/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * getProbAggr / getProbSysAggr (WS-C): the per-class joint at a station and the
 * whole-system joint, fitted to the MVA means over the model's default initial
 * state (populations at their reference stations). Reference numbers are MATLAB
 * SolverMVA getProbAggr / getProbSysAggr.
 */

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_prob.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

TEST_CASE("getProbAggr/SysAggr on a closed product-form network match MATLAB") {
    // Delay(Z=1) + PS Queue(mu=2), N=2. MVA: Q(Delay)=1.2, Q(Q)=0.8.
    // Default state: both jobs at Delay -> nir(Delay)=2, nir(Q)=0.
    // getProbAggr(Delay) = C(2,2)(0.6)^2 = 0.36; getProbAggr(Q) = C(2,0)(0.6)^2 = 0.36.
    qn::Network<double> m("pf");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const mva::AggrResult<double> pa1 = mva::solver_mva_get_prob_aggr(sn, r, 1);
    const mva::AggrResult<double> pa2 = mva::solver_mva_get_prob_aggr(sn, r, 2);
    const mva::AggrResult<double> ps = mva::solver_mva_get_prob_sys_aggr(sn, r);
    CHECK(pa1.P == doctest::Approx(0.36).epsilon(1e-9));
    CHECK(pa2.P == doctest::Approx(0.36).epsilon(1e-9));
    CHECK(ps.P == doctest::Approx(0.36).epsilon(1e-9));
}

TEST_CASE("getProbAggr/SysAggr on an open M/M/1 match MATLAB") {
    // Source(1) -> FCFS Q(mu=2) -> Sink. rho = 0.5. Default state: no jobs.
    // getProbAggr(Q) = (1-rho) = 0.5; getProbSysAggr = 0.5.
    qn::Network<double> m("mm1");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, s, q, 1.0);
    P.set(c, c, q, k, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // station 2 is the Queue (station 1 is the Source)
    const mva::AggrResult<double> pa = mva::solver_mva_get_prob_aggr(sn, r, 2);
    const mva::AggrResult<double> ps = mva::solver_mva_get_prob_sys_aggr(sn, r);
    CHECK(pa.P == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(ps.P == doctest::Approx(0.5).epsilon(1e-9));
}

}  // namespace
