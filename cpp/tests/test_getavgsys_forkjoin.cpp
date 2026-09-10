/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * `getAvgSys` on an OPEN chain through a fork-join section.
 *
 * This case used to be refused by name: the reference fills the join station's
 * response time with `d0`, the mean of the maximum of the parallel branch
 * times, and reaching it needs `ModelAdapter.pathsCS` to enumerate every path
 * from the fork to the join. Both halves are ported in
 * `solvers/solver_chain_tables.h` and this file pins them.
 *
 * ORACLE: live MATLAB, 2026-08-15, `SolverMVA(model).getAvgSys()` on the model
 * of `matlab/examples/basic/forkJoin/fj_basic_open.m`.
 *
 * WHY d0 IS ALSO CHECKED IN CLOSED FORM. The MVA response times feeding the
 * inclusion-exclusion are themselves iterates, so agreeing with MATLAB to 1e-9
 * on CN could in principle mean both codebases made the same mistake in the
 * order statistic. On this model the branch times are the exact M/M/1 values,
 * 1/(mu - lambda), so d0 has the closed form
 *
 *     1/l1 + 1/l2 - 1/(l1 + l2),  l1 = 0.95, l2 = 1.95,
 *
 * and the two agree to 7 digits -- the residue being the MVA iterate, not the
 * order statistic.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/solver_chain_tables.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Source -> Fork -{Queue1, Queue2}- Join -> Sink, one open class. */
qn::Network<double> fj_basic_open() {
    qn::Network<double> m("fj_basic_open");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("class1");
    m.set_arrival(src, c, D::exp_rate(0.05));
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("getAvgSys serves an OPEN chain through a fork-join") {
    qn::Network<double> m = fj_basic_open();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    mva::MvaOptions opt;
    opt.method = "default";
    Matrix<double> init;
    const mva::AvgResult<double> avg = mva::solver_mva_run_analyzer(sn, opt, init);

    const solvers::SysResult<double> sys = solvers::solver_get_avg_sys(sn, avg);
    REQUIRE(sys.CN.size() == 1u);
    // MATLAB SolverMVA(model).getAvgSys(): CN 1.22061419381533, XN 0.05.
    CHECK(sys.CN[0] == doctest::Approx(1.22061419381533).epsilon(1e-9));
    CHECK(sys.XN[0] == doctest::Approx(0.05).epsilon(1e-12));

    // The closed form of the same order statistic on the exact M/M/1 branch
    // times: this is the independent check that d0 is a MAXIMUM and not a sum.
    const double l1 = 1.0 - 0.05, l2 = 2.0 - 0.05;
    const double d0 = 1.0 / l1 + 1.0 / l2 - 1.0 / (l1 + l2);
    CHECK(std::fabs(sys.CN[0] - d0) < 1e-4);
    // A SUM of the branch times would be 1.5654, and the LONGER branch alone
    // 1.0526; both are far outside that band, so the check discriminates.
    CHECK(sys.CN[0] > 1.0 / l1);
    CHECK(sys.CN[0] < 1.0 / l1 + 1.0 / l2);
}

TEST_CASE("exp_max_mean is the mean of the maximum of independent exponentials") {
    namespace cd = solvers::chain_detail;
    // Two branches of mean 1 and 2: 1 + 2 - 1/(1 + 1/2) = 4/3.
    std::vector<double> two;
    two.push_back(1.0);
    two.push_back(2.0);
    CHECK(cd::exp_max_mean(two) == doctest::Approx(1.0 + 2.0 - 1.0 / 1.5).epsilon(1e-14));

    // n identical branches of mean 1: the harmonic number H_n, the classical
    // coupon-collector value, which no term-by-term slip reproduces.
    for (std::size_t n = 1; n <= 6; ++n) {
        std::vector<double> eq(n, 1.0);
        double harmonic = 0.0;
        for (std::size_t k = 1; k <= n; ++k) harmonic += 1.0 / static_cast<double>(k);
        CAPTURE(n);
        CHECK(cd::exp_max_mean(eq) == doctest::Approx(harmonic).epsilon(1e-12));
    }

    // One branch is its own maximum.
    CHECK(cd::exp_max_mean(std::vector<double>(1, 3.5)) == doctest::Approx(3.5).epsilon(1e-14));
    CHECK(cd::exp_max_mean(std::vector<double>()) == 0.0);

    // The maximum is at least the largest mean and at most their sum, always.
    std::vector<double> mixed;
    mixed.push_back(0.5);
    mixed.push_back(4.0);
    mixed.push_back(1.25);
    const double got = cd::exp_max_mean(mixed);
    CHECK(got >= 4.0);
    CHECK(got <= 0.5 + 4.0 + 1.25);

    // A zero-length branch has no exponential rate and is refused, not divided by.
    std::vector<double> bad;
    bad.push_back(1.0);
    bad.push_back(0.0);
    CHECK_THROWS_AS(cd::exp_max_mean(bad), line::NumericError);
    // Past 20 branches the alternating sum has no significant digits left.
    CHECK_THROWS_AS(cd::exp_max_mean(std::vector<double>(21, 1.0)), line::UnsupportedError);
}

TEST_CASE("a CLOSED fork-join chain still takes Little's law, not the path walk") {
    qn::Network<double> m("fj_closed");
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

    mva::MvaOptions opt;
    opt.method = "default";
    Matrix<double> init;
    const mva::AvgResult<double> avg = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const solvers::SysResult<double> sys = solvers::solver_get_avg_sys(m.get_struct(), avg);
    REQUIRE(sys.CN.size() == 1u);
    // N / X, with N = 2: the closed branch never enters the path walk, so this
    // must hold exactly whatever the branch times are.
    CHECK(sys.CN[0] == doctest::Approx(2.0 / sys.XN[0]).epsilon(1e-12));
}

TEST_CASE("a NESTED fork inside an open branch is collapsed before the outer walk") {
    // Source -> Fork1 -{ Q1 , Fork2 -{Q2,Q3}- Join2 }- Join1 -> Sink.
    //
    // This is the branch of pathsCS that mutates RN mid-walk: the inner join
    // takes the inner section's own d0 and the inner branch stations are
    // zeroed, so the outer walk crosses the whole inner section as one station.
    // Getting that wrong does not throw -- it double-counts Q2 and Q3 along the
    // outer path -- so only the reference value catches it.
    qn::Network<double> m("nested");
    const std::size_t src = m.add_source("Source");
    const std::size_t f1 = m.add_fork("Fork1");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t f2 = m.add_fork("Fork2");
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t j2 = m.add_join("Join2", f2);
    const std::size_t j1 = m.add_join("Join1", f1);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("c1");
    m.set_arrival(src, c, D::exp_rate(0.05));
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    m.set_service(q3, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f1, 1.0);
    P.set(f1, q1, 1.0);
    P.set(f1, f2, 1.0);
    P.set(f2, q2, 1.0);
    P.set(f2, q3, 1.0);
    P.set(q2, j2, 1.0);
    P.set(q3, j2, 1.0);
    P.set(j2, j1, 1.0);
    P.set(q1, j1, 1.0);
    P.set(j1, snk, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    opt.method = "default";
    Matrix<double> init;
    const mva::AvgResult<double> avg = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const solvers::SysResult<double> sys = solvers::solver_get_avg_sys(m.get_struct(), avg);
    REQUIRE(sys.CN.size() == 1u);
    // MATLAB SolverMVA(model).getAvgSys() on this model, 2026-08-15.
    CHECK(sys.CN[0] == doctest::Approx(1.29936014710338).epsilon(1e-9));
    CHECK(sys.XN[0] == doctest::Approx(0.05).epsilon(1e-12));
}
