/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression: the dec.source fixed point must conserve the closed population at
 * a multiserver station.
 *
 * `solver_mam_basic` converges a per-chain surrogate arrival rate so that the
 * queue lengths sum to the chain population. Its post-loop second pass then
 * rescaled QN to that population, floored the response time at one full service
 * time (RN = max(S, QN/TN)) and restated QN = RN*TN from an untouched TN --
 * discarding the rescale. With one server the floor is rarely active and the
 * pass is a near no-op; with c > 1 the loop's surrogate response time sits well
 * below S, so the floor roughly doubled RN and both QN and TN came out inflated
 * by the same per-class factor. This model (10 servers, N = 2) returned
 * sum(QN) = 3.52 instead of 2.
 *
 * THE ORACLE IS THE EXACT CTMC, NOT RECORDED OUTPUT. With 10 servers and 2 jobs
 * no job ever queues, so R = S is the exact answer and MAM's own floor makes it
 * reproduce the CTMC outright rather than approximately. The population identity
 * is checked separately, at every server count, so that the two cannot cover for
 * each other: it is a conservation law the decomposition owes at any c.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mam/solver_mam_basic.h"

namespace qn = line::qn;
namespace mam = line::mam;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay + two PS queues, two closed classes of one job each, class-dependent routing. */
qn::Network<double> cqn_multiserver(double nservers) {
    static const double rate[3][2] = {{96.0, 60.0}, {70.0, 84.0}, {96.0, 1.0}};
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t ca = m.add_closed_class("ClassA", 1.0, d);
    const std::size_t cb = m.add_closed_class("ClassB", 1.0, d);
    const std::size_t node[3] = {d, q1, q2};
    const std::size_t cls[2] = {ca, cb};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t r = 0; r < 2; ++r) m.set_service(node[i], cls[r], Dist::exp_rate(rate[i][r]));
        if (i > 0) m.set_number_of_servers(node[i], nservers);
    }
    qn::RoutingMatrix<double> P;
    P.set(ca, ca, d, q1, 0.6);
    P.set(ca, ca, d, q2, 0.4);
    P.set(ca, ca, q1, d, 1.0);
    P.set(ca, ca, q2, d, 1.0);
    P.set(cb, cb, d, q1, 1.0);
    P.set(cb, cb, q1, d, 1.0);
    P.set(cb, cb, q2, d, 1.0);
    m.link(P);
    return m;
}

double column_sum(const line::Matrix<double>& m, std::size_t col) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.rows(); ++i) s += m(i, col);
    return s;
}

}  // namespace

TEST_CASE("mam closed: the chain population survives the second pass at any server count") {
    for (double c : {1.0, 2.0, 10.0}) {
        qn::Network<double> m = cqn_multiserver(c);
        const line::mva::MvaSolution<double> s =
            mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
        // one job per class, and each class is its own chain
        CHECK(column_sum(s.Q, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(column_sum(s.Q, 1) == doctest::Approx(1.0).epsilon(1e-6));
    }
}

TEST_CASE("mam closed: 10 servers and 2 jobs never queue, so MAM reproduces the CTMC") {
    qn::Network<double> m = cqn_multiserver(10.0);
    const line::mva::MvaSolution<double> s =
        mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
    qn::Network<double> mc = cqn_multiserver(10.0);
    const line::mva::AvgResult<double> e = ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());

    for (std::size_t i = 0; i < s.Q.rows(); ++i) {
        for (std::size_t r = 0; r < s.Q.cols(); ++r) {
            CHECK(s.Q(i, r) == doctest::Approx(e.QN(i, r)).epsilon(1e-9));
            CHECK(s.R(i, r) == doctest::Approx(e.RN(i, r)).epsilon(1e-9));
            CHECK(s.Tp(i, r) == doctest::Approx(e.TN(i, r)).epsilon(1e-6));
        }
    }
}
