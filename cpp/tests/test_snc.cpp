/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Stochastic network calculus: the api domain and the SolverBA 'snc.upper' arm.
 *
 * Every reference number here is the MATLAB answer of the same model under
 * `SolverBA(model,'method','snc.upper')`, which is the ground truth for the
 * port, and each is also reproduced by the JAR and by native Python. The
 * codebases differ only in the refinement step of the Chernoff search (MATLAB
 * `fminbnd`, scipy's bounded Brent, and the golden section of `snc_thetaopt`
 * here and in the JAR), which agrees far inside the tolerances used.
 *
 * The exact M/M/1 values are quoted alongside each bound, so that the cases
 * document the tightness rather than only pinning numbers: this is a
 * policy-robust bound, loose on the mean and asymptotically exact in the tail.
 *
 * Tolerances are ABSOLUTE, written as an explicit fabs difference rather than
 * doctest::Approx, whose epsilon is relative.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/snc/snc_bound_backlog.h"
#include "line/api/snc/snc_bound_delay.h"
#include "line/api/snc/snc_conv.h"
#include "line/api/snc/snc_env_cpoisson.h"
#include "line/api/snc/snc_env_map.h"
#include "line/api/snc/snc_env_poisson.h"
#include "line/api/snc/snc_env_tokenbucket.h"
#include "line/api/snc/snc_mean_backlog.h"
#include "line/api/snc/snc_mean_delay.h"
#include "line/api/snc/snc_perc_backlog.h"
#include "line/api/snc/snc_perc_delay.h"
#include "line/api/snc/snc_srv_exp.h"
#include "line/api/snc/snc_srv_rate.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ba/solver_ba_snc.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using Dist = Distrib<double>;

/** Source -> FCFS Queue -> Sink, one open class. */
qn::Network<double> mm1(double lambda, double mu) {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run_snc(qn::Network<double>& m) {
    ba::BaOptions opt;
    opt.method = "snc.upper";
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

}  // namespace

TEST_CASE("snc: the M/M/1 backlog bound decays at the exact rate") {
    // In JOB units the optimal theta tends to log(mu/lambda), which is exactly
    // what makes P{Q>n} <= C*(lambda/mu)^n reproduce the exact decay rate.
    const double lambda = 0.5, mu = 1.0;
    const snc::Envelope arv = snc::snc_env_poisson_fn(lambda);
    const snc::Envelope srv = snc::snc_srv_exp_fn(mu);
    double prevTheta = 0.0;
    const double n[4] = {5.0, 10.0, 20.0, 40.0};
    for (int i = 0; i < 4; ++i) {
        const snc::SncResult r = snc::snc_bound_backlog(arv, srv, n[i]);
        const double exact = std::pow(lambda / mu, n[i] + 1.0);
        CHECK(r.value >= exact);          // it is a bound
        CHECK(r.theta > prevTheta);       // approached from below
        CHECK(r.theta < std::log(mu / lambda) + 1e-12);
        prevTheta = r.theta;
    }
    CHECK(std::fabs(prevTheta - std::log(mu / lambda)) <= 0.05);
}

TEST_CASE("snc: the constant-rate element would UNDERSTATE an Exp server") {
    // The trap snc_srv_exp exists to avoid: a job-counting arrival envelope
    // paired with snc_srv_rate is an M/D/1, whose delay is smaller than the
    // M/M/1 it is meant to bound.
    const snc::Envelope arv = snc::snc_env_poisson_fn(0.6);
    const double md1 = snc::snc_mean_delay(arv, snc::snc_srv_rate_fn(1.0)).value;
    const double mm1b = snc::snc_mean_delay(arv, snc::snc_srv_exp_fn(1.0)).value;
    CHECK(md1 < mm1b);
    CHECK(mm1b > 1.0 / (1.0 - 0.6));  // the M/M/1 arm is above the exact mean
}

TEST_CASE("snc: a one-phase MAP envelope equals the Poisson envelope") {
    Matrix<double> D0(1, 1, -0.6), D1(1, 1, 0.6);
    const double theta[4] = {0.1, 0.5, 1.0, 2.0};
    for (int i = 0; i < 4; ++i) {
        const snc::Env a = snc::snc_env_map(D0, D1, theta[i]);
        const snc::Env b = snc::snc_env_poisson(0.6, theta[i]);
        CHECK(std::fabs(a.sigma - b.sigma) <= 1e-12);
        CHECK(std::fabs(a.rho - b.rho) <= 1e-9);
    }
}

TEST_CASE("snc: an MMPP(2) carries a positive burst term") {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -3.0;
    D0(0, 1) = 1.0;
    D0(1, 0) = 1.0;
    D0(1, 1) = -2.0;
    D1(0, 0) = 2.0;
    D1(1, 1) = 1.0;
    // MATLAB: sigma = 0.525613, rho = 1.591380 at theta = 0.1.
    const snc::Env e = snc::snc_env_map(D0, D1, 0.1);
    CHECK(std::fabs(e.sigma - 0.525613) <= 1e-5);
    CHECK(std::fabs(e.rho - 1.591380) <= 1e-5);
}

TEST_CASE("snc: the token bucket converges to the deterministic delay from above") {
    // b = 5 on C = 1 has deterministic delay b/C = 5; MATLAB reports 5.004605,
    // 5.013816 and 5.027631 at eps = 1e-2, 1e-6 and 1e-12.
    const snc::Envelope tb = snc::snc_env_tokenbucket_fn(5.0, 0.5);
    const snc::Envelope srv = snc::snc_srv_rate_fn(1.0);
    CHECK(std::fabs(snc::snc_perc_delay(tb, srv, 1e-2).value - 5.004605) <= 1e-4);
    CHECK(std::fabs(snc::snc_perc_delay(tb, srv, 1e-6).value - 5.013816) <= 1e-4);
    CHECK(std::fabs(snc::snc_perc_delay(tb, srv, 1e-12).value - 5.027631) <= 1e-4);
}

TEST_CASE("snc: the compound Poisson reads the M/M/1 in units of work") {
    // MATLAB and native Python: d(1e-3) = 22.3947 at lambda = 0.5, mu = 1, the
    // work-unit reading of the M/M/1 (compound Poisson arrivals of Exp work on
    // a constant-rate server), which is a DIFFERENT pairing from the job-unit
    // one and answers a different number for the same queue.
    const snc::Envelope arv = snc::snc_env_cpoisson_fn(0.5, 1.0);
    const snc::Envelope srv = snc::snc_srv_rate_fn(1.0);
    const snc::SncResult d = snc::snc_perc_delay(arv, srv, 1e-3);
    CHECK(std::fabs(d.value - 22.394666) <= 1e-3);
    // and the forward bound at that level returns the same eps
    CHECK(std::fabs(snc::snc_bound_delay(arv, srv, d.value).value - 1e-3) <= 1e-8);
}

TEST_CASE("snc: pay bursts only once beats hop-by-hop") {
    const double rates[3] = {1.5, 1.2, 1.0};
    const snc::Envelope arv = snc::snc_env_poisson_fn(0.6);
    const snc::Envelope endToEnd = [&](double theta) {
        snc::Env acc = snc::snc_srv_exp(rates[0], theta);
        for (int i = 1; i < 3; ++i) acc = snc::snc_conv(acc, snc::snc_srv_exp(rates[i], theta), theta);
        return acc;
    };
    const double concat = snc::snc_perc_delay(arv, endToEnd, 1e-3).value;
    double hopByHop = 0.0;
    for (int i = 0; i < 3; ++i)
        hopByHop += snc::snc_perc_delay(arv, snc::snc_srv_exp_fn(rates[i]), 1e-3 / 3.0).value;
    // MATLAB: 42.5195 concatenated against 65.5128 summed per hop.
    CHECK(std::fabs(concat - 42.5195) <= 1e-2);
    CHECK(std::fabs(hopByHop - 65.5128) <= 1e-2);
    CHECK(concat < hopByHop);
}

TEST_CASE("snc: an unstable composition reports no feasible theta") {
    const snc::Envelope arv = snc::snc_env_poisson_fn(1.5);
    const snc::Envelope srv = snc::snc_srv_exp_fn(1.0);
    CHECK(std::fabs(snc::snc_bound_delay(arv, srv, 10.0).value - 1.0) <= 0.0);
    CHECK(!std::isfinite(snc::snc_perc_delay(arv, srv, 1e-3).value));
    CHECK(!std::isfinite(snc::snc_mean_backlog(arv, srv).value));
}

TEST_CASE("solver_ba_snc: the M/M/1 mean bound matches MATLAB") {
    const double rho[3] = {0.3, 0.6, 0.9};
    const double expectedR[3] = {4.950039, 13.205142, 87.992154};
    for (int i = 0; i < 3; ++i) {
        CAPTURE(rho[i]);
        qn::Network<double> m = mm1(rho[i], 1.0);
        const mva::AvgResult<double> r = run_snc(m);
        CHECK(std::fabs(r.RN(1, 0) - expectedR[i]) <= 1e-4 * expectedR[i]);
        CHECK(r.RN(1, 0) > 1.0 / (1.0 - rho[i]));         // it is a bound
        CHECK(std::fabs(r.QN(1, 0) - rho[i] * r.RN(1, 0)) <= 1e-9);  // Little's law
        CHECK(std::fabs(r.UN(1, 0) - rho[i]) <= 1e-12);   // U is exact
        CHECK(std::fabs(r.TN(1, 0) - rho[i]) <= 1e-12);   // so is T
    }
}

TEST_CASE("solver_ba_snc: the quantiles match MATLAB") {
    qn::Network<double> m = mm1(0.6, 1.0);
    const ba::SncPercentiles p = ba::solver_ba_snc_perc(m.get_struct(), 1e-3);
    // MATLAB: 29.6082 slots and 23.6449 jobs.
    CHECK(std::fabs(p.D(1, 0) - 29.608188) <= 1e-4);
    CHECK(std::fabs(p.B(1, 0) - 23.644943) <= 1e-4);
    // The Source row carries no traffic and stays NaN rather than zero.
    CHECK(std::isnan(p.D(0, 0)));
}

TEST_CASE("solver_ba_snc: a tandem degrades hop by hop") {
    qn::Network<double> m("tandem");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, Dist::exp_rate(0.6));
    m.set_service(q1, c, Dist::exp_rate(1.2));
    m.set_service(q2, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run_snc(m);
    // MATLAB: 7.3397 at Q1 and 20.0114 at Q2, against the exact Jackson 1.6667
    // and 2.5. The second hop is looser because its arrival envelope is the
    // DEPARTURE envelope of the first, which carries the burst the server added.
    CHECK(std::fabs(r.RN(1, 0) - 7.3397) <= 1e-3);
    CHECK(std::fabs(r.RN(2, 0) - 20.0114) <= 1e-3);
}

TEST_CASE("solver_ba_snc: two classes share one server by blind multiplexing") {
    qn::Network<double> m("shared");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Shared", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t a = m.add_open_class("ClassA");
    const std::size_t b = m.add_open_class("ClassB");
    m.set_arrival(src, a, Dist::exp_rate(0.3));
    m.set_arrival(src, b, Dist::exp_rate(0.3));
    m.set_service(q, a, Dist::exp_rate(1.0));
    m.set_service(q, b, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;  // set(r, s, i, j, p): classes first, then nodes
    P.set(a, a, src, q, 1.0);
    P.set(a, a, q, k, 1.0);
    P.set(b, b, src, q, 1.0);
    P.set(b, b, q, k, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run_snc(m);
    // MATLAB: 24.0303 per class, against an exact 2.5 for the aggregate M/M/1.
    CHECK(std::fabs(r.RN(1, 0) - 24.030304) <= 1e-3);
    CHECK(std::fabs(r.RN(1, 1) - 24.030304) <= 1e-3);

    // The backlog quantile of a shared station equals that of the solo station
    // at the aggregate rate: with zero burst terms the bound sees the envelopes
    // only through rhoS - rhoA, and blind multiplexing subtracts the cross rate
    // from rhoS exactly as aggregating the classes would add it to rhoA. The
    // DELAY quantile does separate them, since it scales the level by rhoS.
    const ba::SncPercentiles p = ba::solver_ba_snc_perc(m.get_struct(), 1e-3);
    CHECK(std::fabs(p.B(1, 0) - 23.644943) <= 1e-4);
    CHECK(std::fabs(p.D(1, 0) - 55.901808) <= 1e-3);
}

TEST_CASE("solver_ba_snc: the gates refuse what the calculus cannot carry") {
    // closed model
    qn::Network<double> closed("closed");
    const std::size_t d = closed.add_delay("Think");
    const std::size_t q = closed.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = closed.add_closed_class("C", 3.0, d);
    closed.set_service(d, c, Dist::exp_rate(1.0));
    closed.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> Pc;
    Pc.set(d, q, 1.0);
    Pc.set(q, d, 1.0);
    closed.link(Pc);
    CHECK_THROWS_AS(run_snc(closed), UnsupportedError);

    // probabilistic split downstream of the Source
    qn::Network<double> split("split");
    const std::size_t s2 = split.add_source("Source");
    const std::size_t qa = split.add_queue("QA", SchedStrategy::FCFS);
    const std::size_t qb = split.add_queue("QB", SchedStrategy::FCFS);
    const std::size_t k2 = split.add_sink("Sink");
    const std::size_t c2 = split.add_open_class("C");
    split.set_arrival(s2, c2, Dist::exp_rate(0.4));
    split.set_service(qa, c2, Dist::exp_rate(1.0));
    split.set_service(qb, c2, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> Ps;
    Ps.set(s2, qa, 1.0);
    Ps.set(qa, qb, 0.5);
    Ps.set(qa, k2, 0.5);
    Ps.set(qb, k2, 1.0);
    split.link(Ps);
    CHECK_THROWS_AS(run_snc(split), UnsupportedError);

    // unequal service rates among the classes sharing a station
    qn::Network<double> unequal("unequal");
    const std::size_t s3 = unequal.add_source("Source");
    const std::size_t q3 = unequal.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k3 = unequal.add_sink("Sink");
    const std::size_t a3 = unequal.add_open_class("A");
    const std::size_t b3 = unequal.add_open_class("B");
    unequal.set_arrival(s3, a3, Dist::exp_rate(0.2));
    unequal.set_arrival(s3, b3, Dist::exp_rate(0.2));
    unequal.set_service(q3, a3, Dist::exp_rate(1.0));
    unequal.set_service(q3, b3, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> Pu;
    Pu.set(a3, a3, s3, q3, 1.0);
    Pu.set(a3, a3, q3, k3, 1.0);
    Pu.set(b3, b3, s3, q3, 1.0);
    Pu.set(b3, b3, q3, k3, 1.0);
    unequal.link(Pu);
    CHECK_THROWS_AS(run_snc(unequal), UnsupportedError);
}

TEST_CASE("solver_ba_snc: list_valid_methods narrows to the open families") {
    qn::Network<double> m = mm1(0.6, 1.0);
    const std::vector<std::string> v = ba::list_valid_methods(m.get_struct());
    CHECK(v.size() == 3);
    CHECK(std::find(v.begin(), v.end(), "snc.upper") != v.end());
    CHECK(std::find(v.begin(), v.end(), "bpt.lower") != v.end());
    CHECK(std::find(v.begin(), v.end(), "bgt.upper") != v.end());
}
