/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The multiclass open polling analyzer (ladder branch 5).
 *
 * Numbers are MATLAB's `SolverMVA(model).getAvg()` on a two-class exhaustive
 * polling system. Utilization and throughput are exact (rho_r = lambda_r/mu_r
 * and T_r = lambda_r regardless of the discipline); the waiting time carries the
 * pseudo-conservation law, so R and Q are what pin the analyzer.
 */

#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::PollingType;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

qn::Network<double> polling_model(PollingType pt) {
    qn::Network<double> m("polling");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::POLLING);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, D::exp_rate(0.3));
    m.set_arrival(src, c2, D::exp_rate(0.2));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    // MATLAB's setSwitchover(Exp(...)) does not propagate into the polling
    // nodeparam on this path: its reference runs with the Immediate default,
    // which is Exp(1e8) (mean 1e-8), so the port uses the same near-zero
    // switchover to compare like with like. Non-zero switchover moments are
    // exercised at the polling_qsys_* API level in test_polling_takagi.
    m.set_switchover(q, c1, D::exp_rate(1e8));
    m.set_switchover(q, c2, D::exp_rate(1e8));
    m.set_polling_type(q, pt);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("exhaustive polling reproduces the reference waiting times") {
    qn::Network<double> m = polling_model(PollingType::EXHAUSTIVE);
    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    CHECK(r.actualmethod == "stationtime");
    // stations: Source (0), Q (1); the Sink is not a station
    CHECK(r.QN(1, 0) == doctest::Approx(0.1856625809).epsilon(1e-8));
    CHECK(r.QN(1, 1) == doctest::Approx(0.0938465647).epsilon(1e-8));
    CHECK(r.RN(1, 0) == doctest::Approx(0.6188752696).epsilon(1e-8));
    CHECK(r.RN(1, 1) == doctest::Approx(0.4692328236).epsilon(1e-8));
    // exact regardless of the discipline
    CHECK(r.UN(1, 0) == doctest::Approx(0.15).epsilon(1e-9));
    CHECK(r.UN(1, 1) == doctest::Approx(0.2 / 3.0).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(0.3).epsilon(1e-9));
    CHECK(r.TN(1, 1) == doctest::Approx(0.2).epsilon(1e-9));
    // Little's law holds per class at the polling station
    CHECK(r.QN(1, 0) == doctest::Approx(r.TN(1, 0) * r.RN(1, 0)).epsilon(1e-9));
    CHECK(r.QN(1, 1) == doctest::Approx(r.TN(1, 1) * r.RN(1, 1)).epsilon(1e-9));
}

TEST_CASE("gated polling reproduces the reference waiting times") {
    qn::Network<double> m = polling_model(PollingType::GATED);
    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    CHECK(r.actualmethod == "stationtime");
    // gated serves only what was present at the visit instant, so its waiting
    // times differ from exhaustive; MATLAB SolverMVA getAvg on the same model
    CHECK(r.QN(1, 0) == doctest::Approx(0.1883299417).epsilon(1e-8));
    CHECK(r.QN(1, 1) == doctest::Approx(0.0898455255).epsilon(1e-8));
    CHECK(r.RN(1, 0) == doctest::Approx(0.6277664724).epsilon(1e-8));
    CHECK(r.RN(1, 1) == doctest::Approx(0.4492276276).epsilon(1e-8));
}

TEST_CASE("K-limited polling is refused above K = 1") {
    qn::Network<double> m("kpoll");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::POLLING);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, D::exp_rate(0.3));
    m.set_arrival(src, c2, D::exp_rate(0.2));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    m.set_switchover(q, c1, D::exp_rate(10.0));
    m.set_switchover(q, c2, D::exp_rate(10.0));
    m.set_polling_type(q, PollingType::KLIMITED, 3);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    Matrix<double> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
}

}  // namespace
