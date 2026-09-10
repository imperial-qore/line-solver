/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The state-probability half of the SolverMVA surface.
 *
 * The numbers are MATLAB's `getProbMarg` / `getProbNormConstAggr`, and the two
 * open cases are additionally checkable in closed form without either codebase:
 * an M/M/infinity with mean 0.5 is Poisson(0.5), and a single-class M/M/1 at
 * rho = 0.25 is 0.75 * 0.25^n. Both agree, so the reference values are not
 * merely self-consistent.
 */

#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_prob.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

TEST_CASE("getProbMarg fits the Schmidt binomial on a closed class") {
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

    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> avg = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const mva::MargResult<double> pm =
        mva::solver_mva_get_prob_marg(m.get_struct(), avg, 2, 1, std::vector<long>());

    // n = 0..3, since the class population is 3
    REQUIRE(pm.P.size() == 4);
    CHECK(pm.P[0] == doctest::Approx(0.145793847500).epsilon(1e-9));
    CHECK(pm.P[1] == doctest::Approx(0.393643388249).epsilon(1e-9));
    CHECK(pm.P[2] == doctest::Approx(0.354279049424).epsilon(1e-9));
    CHECK(pm.P[3] == doctest::Approx(0.106283714827).epsilon(1e-9));
    CHECK(pm.logP[0] == doctest::Approx(-1.925561658517).epsilon(1e-9));
    CHECK(pm.logP[3] == doctest::Approx(-2.241643205491).epsilon(1e-9));
    // a binomial is a proper law, so it sums to one; MVA's mean is only an
    // approximation of the true marginal, but the FIT is exact
    double s = 0.0;
    for (double p : pm.P) s += p;
    CHECK(s == doctest::Approx(1.0).epsilon(1e-12));

    // and the states can be selected
    const mva::MargResult<double> sel =
        mva::solver_mva_get_prob_marg(m.get_struct(), avg, 2, 1, std::vector<long>{0, 3});
    REQUIRE(sel.P.size() == 2);
    CHECK(sel.P[0] == doctest::Approx(0.145793847500).epsilon(1e-9));
    CHECK(sel.P[1] == doctest::Approx(0.106283714827).epsilon(1e-9));
    CHECK_THROWS_AS(mva::solver_mva_get_prob_marg(m.get_struct(), avg, 2, 1, std::vector<long>{4}),
                    InputError);
    CHECK_THROWS_AS(mva::solver_mva_get_prob_marg(m.get_struct(), avg, 9, 1, std::vector<long>()),
                    InputError);

    // the normalizing constant, which only the exact recursion carries
    CHECK(mva::solver_mva_get_prob_norm_const_aggr(m.get_struct(), opt) ==
          doctest::Approx(-0.233614851182).epsilon(1e-9));
}

TEST_CASE("getProbMarg uses the exact open product form, per discipline") {
    qn::Network<double> m("oqn");
    const std::size_t src = m.add_source("Source");
    const std::size_t dl = m.add_delay("Delay");
    const std::size_t qq = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, D::exp_rate(1.0));
    m.set_service(dl, o, D::exp_rate(2.0));
    m.set_service(qq, o, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, dl, 1.0);
    P.set(dl, qq, 1.0);
    P.set(qq, sk, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    Matrix<double> init;
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mva::AvgResult<double> avg = mva::solver_mva_run_analyzer(sn, opt, init);
    // station 2 is the Delay and station 3 the Queue, matching MATLAB's order
    CHECK(avg.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-8));
    CHECK(avg.UN(2, 0) == doctest::Approx(0.25).epsilon(1e-8));

    const std::vector<long> nn{0, 1, 2, 3, 4};
    // an infinite server carries an independent Poisson per class, here with
    // mean 0.5, so P(n) = e^-0.5 0.5^n / n!
    const mva::MargResult<double> pd = mva::solver_mva_get_prob_marg(sn, avg, 2, 1, nn);
    REQUIRE(pd.P.size() == 5);
    CHECK(pd.P[0] == doctest::Approx(0.606530659995).epsilon(1e-8));
    CHECK(pd.P[1] == doctest::Approx(0.303265329715).epsilon(1e-8));
    CHECK(pd.P[2] == doctest::Approx(0.075816332358).epsilon(1e-8));
    CHECK(pd.P[3] == doctest::Approx(0.012636055381).epsilon(1e-8));
    CHECK(pd.P[4] == doctest::Approx(0.001579506921).epsilon(1e-8));

    // a queueing station carries the multiclass geometric; with one class it
    // collapses to the M/M/1 law (1-rho) rho^n at rho = 0.25
    const mva::MargResult<double> pq = mva::solver_mva_get_prob_marg(sn, avg, 3, 1, nn);
    REQUIRE(pq.P.size() == 5);
    CHECK(pq.P[0] == doctest::Approx(0.75).epsilon(1e-9));
    CHECK(pq.P[1] == doctest::Approx(0.1875).epsilon(1e-9));
    CHECK(pq.P[2] == doctest::Approx(0.046875).epsilon(1e-9));
    CHECK(pq.P[3] == doctest::Approx(0.01171875).epsilon(1e-9));
    CHECK(pq.P[4] == doctest::Approx(0.0029296875).epsilon(1e-9));

    // the Source is not a queue and has no distribution
    const mva::MargResult<double> ps =
        mva::solver_mva_get_prob_marg(sn, avg, 1, 1, std::vector<long>());
    REQUIRE(ps.P.size() == 1);
    CHECK(ps.P[0] == doctest::Approx(1.0));

    // with no states given, each open case picks the reference's own range
    const mva::MargResult<double> auto_q =
        mva::solver_mva_get_prob_marg(sn, avg, 3, 1, std::vector<long>());
    CHECK(auto_q.P.size() > 5);
    double s = 0.0;
    for (double p : auto_q.P) s += p;
    CHECK(s == doctest::Approx(1.0).epsilon(1e-9));
}

}  // namespace
