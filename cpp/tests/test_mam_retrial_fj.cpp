/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The two impatience/fork-join routes of SolverMAM: `solver_mam_retrial` and
 * `solver_mam_fj`.
 *
 * THE RETRIAL ORACLES ARE CLOSED FORMS, NOT RECORDED OUTPUT. The M/M/1 retrial
 * queue with a linear retrial rate has a known mean orbit length,
 *
 *     E[N_orbit] = rho / (1 - rho) * (rho + lambda / nu),   rho = lambda / mu,
 *
 * (Falin and Templeton, "Retrial Queues", Chapman and Hall 1997), which reduces
 * to the M/M/1 waiting line rho^2 / (1 - rho) as nu -> infinity, i.e. when an
 * orbiting job retries instantly and the orbit IS the waiting line. That limit
 * is the "a retrial queue with instantaneous retrials is the corresponding
 * non-retrial queue" identity, and it is checked separately from the formula so
 * that the two cannot cover for each other. Alongside them sits the work
 * conservation identity: no job is lost when the orbit abandonment rate and the
 * batch rejection probability are zero, so the throughput is the arrival rate
 * and the utilization is lambda * b1 / N for ANY retrial rate, which is what
 * pins the multiserver phase-type case where no closed form is at hand.
 *
 * Each of these was confirmed against MATLAB's own qsys_bmapphnn_retrial before
 * being written down; the numbers quoted in the comments are that run's, and
 * they agree with the closed forms to every digit printed.
 *
 * THE FORK-JOIN TESTS ARE MOSTLY GATE TESTS: a predicate is checked by the
 * models it must accept and reject, and an extraction by whether the descriptor
 * it builds carries the rates the model declared. The last one runs the engine.
 *
 * APPROVED CHANGE (2026-07-31): "the analyzer refuses by naming the engine it
 * wraps" asserted that `solver_mam_fj` throws, pinning the ABSENCE of mainFJ.
 * The engine is now ported (`api/fj/fj_codes.h`), so the row is replaced by the
 * positive test rather than relaxed. The FJ_codes layer itself is validated
 * against MATLAB's own `returnRT1` / `returnRT2` / `mainFJ` on the parameter
 * sets of `main_example.m`; see `_kb/06-solver-catalog.md` for the measured
 * agreement.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_fj.h"
#include "line/solvers/mam/solver_mam_retrial.h"

namespace qn = line::qn;
namespace mam = line::mam;
using line::UnsupportedError;
using line::lang::DropStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/**
 * Source -> bufferless Queue with an orbit -> Sink, one open class.
 *
 * `set_capacity(q, servers)` is what makes the station bufferless: with the
 * capacity equal to the server count there is no waiting line, so a blocked job
 * has nowhere to go but the orbit. Both that and the RETRIAL drop rule are what
 * `qsys_is_retrial` keys on.
 */
qn::Network<double> retrial_queue(double lambda, const Dist& service, double nu, int servers) {
    qn::Network<double> m("retrial");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Orbit", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, service);
    m.set_number_of_servers(q, static_cast<double>(servers));
    m.set_capacity(q, static_cast<double>(servers));
    m.set_drop_rule(q, c, DropStrategy::RETRIAL);
    m.set_retrial(q, c, Dist::exp_rate(nu), nu);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    return m;
}

/** Source -> Fork -> K parallel FCFS queues -> Join -> Sink, one open class. */
qn::Network<double> open_fork_join(double lambda, const std::vector<double>& branch_rates) {
    qn::Network<double> m("openfj");
    const std::size_t src = m.add_source("Src");
    const std::size_t f = m.add_fork("Fork");
    std::vector<std::size_t> qs;
    for (std::size_t k = 0; k < branch_rates.size(); ++k)
        qs.push_back(m.add_queue("Q" + std::to_string(k + 1), SchedStrategy::FCFS));
    const std::size_t j = m.add_join("Join", f);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    for (std::size_t k = 0; k < qs.size(); ++k)
        m.set_service(qs[k], c, Dist::exp_rate(branch_rates[k]));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    for (std::size_t k = 0; k < qs.size(); ++k) {
        P.set(f, qs[k], 1.0);
        P.set(qs[k], j, 1.0);
    }
    P.set(j, snk, 1.0);
    m.link(P);
    return m;
}

/** The Falin-Templeton mean orbit length of the M/M/1 retrial queue. */
double falin_orbit(double lambda, double mu, double nu) {
    const double rho = lambda / mu;
    return rho / (1.0 - rho) * (rho + lambda / nu);
}

}  // namespace

TEST_CASE("mam retrial: the M/M/1 orbit reproduces the Falin-Templeton mean") {
    // MATLAB, qsys_bmapphnn_retrial({-0.5,0.5}, 1, -1, 1, 1, 0, 0, 0, 'MaxLevel', 400):
    // L_orbit = 1.0000000000, N_server = 0.5, Utilization = 0.5, Throughput = 0.5,
    // with 2.8e-24 of mass left at the top level.
    const double lambda = 0.5, mu = 1.0, nu = 1.0;
    qn::Network<double> m = retrial_queue(lambda, Dist::exp_rate(mu), nu, 1);
    const line::mva::MvaSolution<double> s =
        mam::solver_mam_retrial(m.get_struct(), mam::MamOptions());

    const double rho = lambda / mu;
    // the station holds the orbit AND the job in service
    CHECK(s.Q(1, 0) == doctest::Approx(falin_orbit(lambda, mu, nu) + rho).epsilon(1e-9));
    // work conservation: nothing is lost, so every arrival is eventually served
    CHECK(s.U(1, 0) == doctest::Approx(rho).epsilon(1e-9));
    CHECK(s.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-9));
    // Little's law over the whole station, orbit included
    CHECK(s.R(1, 0) == doctest::Approx(s.Q(1, 0) / lambda).epsilon(1e-12));
    CHECK(s.X[0] == doctest::Approx(lambda).epsilon(1e-9));
    CHECK(s.C[0] == doctest::Approx(s.R(1, 0)).epsilon(1e-12));
    // `totiter` for this analyzer is the truncation level that produced the answer
    CHECK(s.iter >= 100);
}

TEST_CASE("mam retrial: the orbit mean tracks the retrial rate across two decades") {
    // The retrial rate is the only parameter separating these three solves, and
    // the closed form says the orbit shrinks as nu grows while the utilization
    // does not move at all. MATLAB reports L_orbit = 0.8, 0.4 and 0.2666933333
    // for nu = 0.5, 2 and 1e4.
    const double lambda = 0.4, mu = 1.0;
    const double rates[3] = {0.5, 2.0, 1e4};
    for (int i = 0; i < 3; ++i) {
        qn::Network<double> m = retrial_queue(lambda, Dist::exp_rate(mu), rates[i], 1);
        const line::mva::MvaSolution<double> s =
            mam::solver_mam_retrial(m.get_struct(), mam::MamOptions());
        CHECK(s.Q(1, 0) ==
              doctest::Approx(falin_orbit(lambda, mu, rates[i]) + lambda / mu).epsilon(1e-8));
        CHECK(s.U(1, 0) == doctest::Approx(lambda / mu).epsilon(1e-9));
        CHECK(s.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-9));
    }
}

TEST_CASE("mam retrial: an instantaneous retrial collapses the orbit onto the M/M/1 queue") {
    // With nu = 1e4 an orbiting job re-attempts far faster than anything else in
    // the model happens, so the orbit is the M/M/1 waiting line and the station
    // holds rho/(1-rho) jobs. The gap is O(lambda/nu) = 4e-5, which is why this
    // is checked to 1e-4 and not to machine precision.
    const double lambda = 0.4, mu = 1.0, rho = lambda / mu;
    qn::Network<double> m = retrial_queue(lambda, Dist::exp_rate(mu), 1e4, 1);
    const line::mva::MvaSolution<double> s =
        mam::solver_mam_retrial(m.get_struct(), mam::MamOptions());
    CHECK(s.Q(1, 0) == doctest::Approx(rho / (1.0 - rho)).epsilon(1e-4));
    CHECK(s.R(1, 0) == doctest::Approx(1.0 / (mu - lambda)).epsilon(1e-4));
}

TEST_CASE("mam retrial: the utilization law holds for phase-type service on two servers") {
    // No closed form is at hand for an Erlang(2) orbit queue on two servers, but
    // work conservation still is: with no abandonment and no batch loss the
    // carried rate is the offered rate and the utilization is lambda*b1/N.
    // Erlang(2) at phase rate 2 has mean 1, so b1 = 1. MATLAB,
    // qsys_bmapphnn_retrial({-0.8,0.8}, [1 0], [-2 2;0 -2], 2, 1, 0, 0, 1):
    // Utilization = 0.4000000000, Throughput = 0.8000000000.
    const double lambda = 0.8;
    qn::Network<double> m = retrial_queue(lambda, Dist::erlang(2.0, 2), 1.0, 2);
    const line::mva::MvaSolution<double> s =
        mam::solver_mam_retrial(m.get_struct(), mam::MamOptions());
    CHECK(s.U(1, 0) == doctest::Approx(0.4).epsilon(1e-9));
    CHECK(s.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-9));
    // the station holds strictly more than the busy servers: the orbit is not empty
    CHECK(s.Q(1, 0) > 0.8);
}

TEST_CASE("mam retrial: a station without an orbit is not claimed") {
    // A plain M/M/1 has an unbounded waiting line and no retrial drop rule, so
    // the detector must decline it rather than answering it as a retrial queue.
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK_FALSE(mam::mam_retrial_detect(sn).ok);
    CHECK_THROWS_AS(mam::solver_mam_retrial(sn, mam::MamOptions()), UnsupportedError);
}

TEST_CASE("mam retrial: a non-phase-type service is refused rather than Erlang-fitted") {
    // The reference replaces a Det service by an Erlang of matching mean and
    // warns; there is no warning channel here, and the substitute is a different
    // model.
    qn::Network<double> m = retrial_queue(0.5, Dist::det(1.0), 1.0, 1);
    CHECK_THROWS_AS(mam::solver_mam_retrial(m.get_struct(), mam::MamOptions()), UnsupportedError);
}

TEST_CASE("mam retrial: a bounded retry count is refused by name") {
    // The generator indexes levels by the orbit population alone, so a job on
    // its last attempt cannot be told apart from a persistent one.
    qn::Network<double> m("retrial_limited");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Orbit", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(1.0));
    m.set_capacity(q, 1.0);
    m.set_drop_rule(q, c, DropStrategy::RETRIAL_WITH_LIMIT);
    m.set_retrial(q, c, Dist::exp_rate(1.0), 1.0, 3);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(mam::mam_retrial_detect(sn).ok);  // the topology IS a retrial topology
    CHECK_THROWS_AS(mam::solver_mam_retrial(sn, mam::MamOptions()), UnsupportedError);
}

TEST_CASE("mam retrial: a multiclass model is declined, as the reference declines it") {
    qn::Network<double> m("retrial2c");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Orbit", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dist::exp_rate(0.3));
    m.set_arrival(src, c2, Dist::exp_rate(0.2));
    m.set_service(q, c1, Dist::exp_rate(1.0));
    m.set_service(q, c2, Dist::exp_rate(1.0));
    m.set_capacity(q, 1.0);
    m.set_drop_rule(q, c1, DropStrategy::RETRIAL);
    m.set_drop_rule(q, c2, DropStrategy::RETRIAL);
    m.set_retrial(q, c1, Dist::exp_rate(1.0), 1.0);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mam::MamRetrialInfo info = mam::mam_retrial_detect(sn);
    CHECK_FALSE(info.ok);
    CHECK(info.why.find("single class") != std::string::npos);
}

TEST_CASE("mam fj: the homogeneity gate accepts an open two-branch fork-join") {
    qn::Network<double> m = open_fork_join(0.5, std::vector<double>{2.0, 2.0});
    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(m.get_struct());
    CHECK(info.ok);
    CHECK(info.K == 2);
    CHECK(info.queueNodes.size() == 2);
    CHECK(info.forkNode != 0);
    CHECK(info.joinNode != 0);
}

TEST_CASE("mam fj: heterogeneous branches are rejected, since one queue no longer stands for K") {
    qn::Network<double> m = open_fork_join(0.5, std::vector<double>{2.0, 3.0});
    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(m.get_struct());
    CHECK_FALSE(info.ok);
    CHECK(info.why.find("heterogeneous") != std::string::npos);
}

TEST_CASE("mam fj: a single-branch fork carries the model's own rates into the descriptor") {
    // K = 1 is the degenerate fork: there is nothing to synchronise, so the
    // descriptor the extraction builds must be the plain Source and Queue of the
    // unforked model, rate for rate.
    const double lambda = 0.5, mu = 2.0;
    qn::Network<double> m = open_fork_join(lambda, std::vector<double>{mu});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(sn);
    REQUIRE(info.ok);
    CHECK(info.K == 1);
    const mam::MamFjParams<double> par = mam::mam_fj_extract_params(sn, info);
    CHECK(par.K == 1);
    REQUIRE(par.arrival.size() == 1);
    CHECK(par.arrival[0].lambda == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(par.service[0].mu == doctest::Approx(mu).epsilon(1e-12));
    // choice 1 is the algorithm's code for an exponential
    CHECK(par.arrival[0].choice == 1);
    CHECK(par.service[0].choice == 1);
}

TEST_CASE("mam fj: an unstable branch is refused where the model can still be named") {
    // Every branch sees the whole arrival stream, so lambda >= mu is instability
    // per branch. mainFJ errors on it; the reference only warns here and lets
    // the engine object one call later.
    qn::Network<double> m = open_fork_join(2.0, std::vector<double>{1.0, 1.0});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(sn);
    REQUIRE(info.ok);
    CHECK_THROWS(mam::mam_fj_extract_params(sn, info));
}

TEST_CASE("mam fj: a closed fork-join is outside the class the approximation is defined on") {
    qn::Network<double> m("closedfj");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(m.get_struct());
    CHECK_FALSE(info.ok);
    CHECK(info.why.find("open") != std::string::npos);
}

TEST_CASE("mam fj: the analyzer reports the reference's own metric tuple") {
    // lambda = 0.5 into two exponential branches of rate 2. Every branch sees
    // the whole stream, so the reference's M/M/1 stand-ins are exact here:
    // U = 0.25, Q = rho/(1-rho) = 1/3, R = 1/(mu-lambda) = 2/3. The Join then
    // carries the synchronisation delay, which is the only number in the tuple
    // that mainFJ produces.
    qn::Network<double> m = open_fork_join(0.5, std::vector<double>{2.0, 2.0});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    mam::MamOptions opt;
    opt.fj_accuracy = 20;  // the default 100 costs 25x the matrices for the same 3 digits
    const line::mva::MvaSolution<double> s = mam::solver_mam_fj(sn, opt);

    const mam::MamFjInfo info = mam::mam_fj_is_homogeneous(sn);
    REQUIRE(info.ok);
    REQUIRE(info.K == 2);
    for (std::size_t k = 0; k < info.K; ++k) {
        const std::size_t st = sn.nodes[info.queueNodes[k] - 1].station;
        CHECK(s.U(st - 1, 0) == doctest::Approx(0.25).epsilon(1e-12));
        CHECK(s.Tp(st - 1, 0) == doctest::Approx(0.5).epsilon(1e-12));
        CHECK(s.Q(st - 1, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
        CHECK(s.R(st - 1, 0) == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
    }
    // The synchronisation delay is positive and Little's law holds on it: with
    // two branches the fork-join mean is strictly above one branch's own.
    const std::size_t jst = sn.nodes[info.joinNode - 1].station;
    REQUIRE(jst != 0);
    const double sync = s.R(jst - 1, 0);
    CHECK(sync > 0.0);
    CHECK(s.Q(jst - 1, 0) == doctest::Approx(0.5 * sync).epsilon(1e-12));
    CHECK(s.Tp(jst - 1, 0) == doctest::Approx(0.5).epsilon(1e-12));
}

TEST_CASE("mam fj: the stored percentiles are the four solver_mam_fj.m keeps") {
    qn::Network<double> m = open_fork_join(0.5, std::vector<double>{2.0, 2.0});
    mam::MamOptions opt;
    opt.fj_accuracy = 20;
    const std::vector<double> want = {0.50, 0.90, 0.99};
    const std::vector<std::vector<double> > p =
        mam::solver_mam_fj_percentiles(m.get_struct(), opt, want);
    REQUIRE(p.size() == 1);
    REQUIRE(p[0].size() == 3);
    // A response-time quantile function is nondecreasing in the level.
    CHECK(p[0][0] > 0.0);
    CHECK(p[0][1] > p[0][0]);
    CHECK(p[0][2] > p[0][1]);

    // A level the reference does not store is INTERPOLATED off the stored four,
    // linearly, as getPerctRespT.m's interp1(..., 'linear', 'extrap') does. 0.70
    // is the midpoint of the [0.50, 0.90] segment.
    const std::vector<std::vector<double> > q =
        mam::solver_mam_fj_percentiles(m.get_struct(), opt, std::vector<double>{0.70});
    REQUIRE(q.size() == 1);
    CHECK(q[0][0] == doctest::Approx(0.5 * (p[0][0] + p[0][1])).epsilon(1e-12));
}
