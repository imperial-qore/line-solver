/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The CTMC generator assembly. The oracle is analytic: a finite-buffer M/M/1
 * has a geometric stationary distribution, so the generator built from the
 * ported state space and event handlers is checked against a closed form
 * rather than against itself.
 */
#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc.h"

namespace qn = line::qn;
using line::lang::EventType;
using line::lang::SchedStrategy;

namespace {

/** Source -> FCFS Queue (capacity K) -> Sink, one open class. */
qn::Network<double> mm1k(double lambda, double mu, int K) {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the generator of an M/M/1/K matches the analytic chain") {
    const double lambda = 0.6, mu = 1.0;
    const int K = 3;
    qn::Network<double> m = mm1k(lambda, mu, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{static_cast<std::size_t>(K)});
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);
    const line::ctmc::CtmcResult<double> r = line::ctmc::solver_ctmc(sn, ss, sy);

    // Every row of an infinitesimal generator sums to zero: that is the
    // defining property, and it is what a mismatched successor index would
    // break first.
    for (std::size_t i = 0; i < r.Q.rows(); ++i) {
        double s = 0;
        for (std::size_t j = 0; j < r.Q.cols(); ++j) s += r.Q(i, j);
        CHECK(s == doctest::Approx(0.0).epsilon(1e-12));
    }
    // Off-diagonals are rates, so non-negative; the diagonal is their negative.
    for (std::size_t i = 0; i < r.Q.rows(); ++i)
        for (std::size_t j = 0; j < r.Q.cols(); ++j)
            if (i != j) CHECK(r.Q(i, j) >= 0.0);

    const std::vector<double> pi = line::mc::ctmc_solve(r.Q);
    REQUIRE(pi.size() == r.Q.rows());

    // Aggregate the stationary law by queue length and compare with the
    // truncated geometric p_n = rho^n (1-rho) / (1-rho^(K+1)).
    std::vector<double> byq(K + 1, 0.0);
    for (std::size_t i = 0; i < pi.size(); ++i) {
        const qn::NetState<double>& st = r.space[i];
        // Stateful nodes are Source then Queue; the Queue block is [buf | srv].
        const std::vector<double>& row = st.local[1];
        double nq = 0;
        for (std::size_t j = 0; j < row.size(); ++j)
            nq += row[j] > 0 ? 1.0 : 0.0;  // one tag per waiting job, one per server
        REQUIRE(static_cast<int>(nq) <= K);
        byq[static_cast<std::size_t>(nq)] += pi[i];
    }
    const double rho = lambda / mu;
    double norm = 0;
    for (int i = 0; i <= K; ++i) norm += std::pow(rho, i);
    for (int i = 0; i <= K; ++i)
        CHECK(byq[i] == doctest::Approx(std::pow(rho, i) / norm).epsilon(1e-9));
}

TEST_CASE("CTMC mean metrics of an M/M/1/K match the analytic values") {
    const double lambda = 0.6, mu = 1.0;
    const int K = 4;
    qn::Network<double> m = mm1k(lambda, mu, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{static_cast<std::size_t>(K)});
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);
    const line::ctmc::CtmcResult<double> r = line::ctmc::solver_ctmc(sn, ss, sy);
    const std::vector<double> pi = line::mc::ctmc_solve(r.Q);
    const line::ctmc::CtmcAvg<double> avg = line::ctmc::solver_ctmc_avg_from_pi(sn, r, pi);

    // Analytic M/M/1/K: p_n = rho^n / sum, so Q = sum n p_n, the carried
    // throughput is mu(1-p_0) and utilization is 1-p_0.
    const double rho = lambda / mu;
    double norm = 0;
    for (int i = 0; i <= K; ++i) norm += std::pow(rho, i);
    double q = 0;
    for (int i = 0; i <= K; ++i) q += i * std::pow(rho, i) / norm;
    const double p0 = 1.0 / norm;

    // Station 1 is the Source, station 2 the Queue.
    CHECK(avg.QN(1, 0) == doctest::Approx(q).epsilon(1e-9));
    CHECK(avg.UN(1, 0) == doctest::Approx(1.0 - p0).epsilon(1e-9));
    CHECK(avg.TN(1, 0) == doctest::Approx(mu * (1.0 - p0)).epsilon(1e-9));
    // Little's law must hold on the ported quantities, not just analytically.
    CHECK(avg.RN(1, 0) == doctest::Approx(q / (mu * (1.0 - p0))).epsilon(1e-9));

    // A Source holds no jobs: its Inf marginal is an encoding sentinel, and
    // reading it as a queue length is what made an earlier CTMC report Q = Inf.
    CHECK(avg.QN(0, 0) == doctest::Approx(0.0));
    CHECK(std::isfinite(avg.RN(0, 0)));
}

TEST_CASE("CTMC reproduces the M/M/1 queue length as the buffer grows") {
    const double lambda = 0.5, mu = 1.0;
    const int K = 12;
    qn::Network<double> m = mm1k(lambda, mu, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{static_cast<std::size_t>(K)});
    const line::ctmc::CtmcResult<double> r =
        line::ctmc::solver_ctmc(sn, ss, qn::refresh_sync(sn));
    const std::vector<double> pi = line::mc::ctmc_solve(r.Q);
    const line::ctmc::CtmcAvg<double> avg = line::ctmc::solver_ctmc_avg_from_pi(sn, r, pi);

    // With a deep buffer the truncated chain approaches M/M/1: Q = rho/(1-rho)
    // = 1 and U = rho. The residual is the tail beyond K, so it shrinks
    // geometrically rather than being an approximation of the method.
    CHECK(avg.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-3));
    CHECK(avg.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-3));
    CHECK(avg.TN(1, 0) == doctest::Approx(0.5).epsilon(1e-3));
}

TEST_CASE("the reachable space of an M/M/1/K is the whole encoded space") {
    const int K = 4;
    qn::Network<double> m = mm1k(0.6, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::NetState<double>> full =
        qn::space_generator(sn, std::vector<std::size_t>{static_cast<std::size_t>(K)});
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);

    // Start empty: the Source's block plus an all-zero queue.
    qn::NetState<double> init = full[0];
    const std::vector<qn::NetState<double>> reach =
        line::ctmc::reachable_space_generator(sn, init, sy);

    // An ordinary FCFS queue can reach every encoded state, so the two agree.
    // They diverge only where the encoding is wider than the dynamics.
    CHECK(reach.size() == full.size());

    // The generator over the reachable space gives the same answers.
    const line::ctmc::CtmcResult<double> r = line::ctmc::solver_ctmc(sn, reach, sy);
    const std::vector<double> pi = line::mc::ctmc_solve(r.Q);
    const line::ctmc::CtmcAvg<double> avg = line::ctmc::solver_ctmc_avg_from_pi(sn, r, pi);
    const double rho = 0.6;
    double norm = 0, q = 0;
    for (int i = 0; i <= K; ++i) norm += std::pow(rho, i);
    for (int i = 0; i <= K; ++i) q += i * std::pow(rho, i) / norm;
    CHECK(avg.QN(1, 0) == doctest::Approx(q).epsilon(1e-9));
}

TEST_CASE("the reachable space of a closed network is smaller than the encoded one") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, line::lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<qn::NetState<double>> full =
        qn::space_generator(sn, std::vector<std::size_t>{0});
    const std::vector<qn::Sync<double>> sy = qn::refresh_sync(sn);
    const std::vector<qn::NetState<double>> reach =
        line::ctmc::reachable_space_generator(sn, full[0], sy);

    // Population is conserved, so every reachable state holds exactly 2 jobs.
    CHECK(!reach.empty());
    CHECK(reach.size() <= full.size());

    const line::ctmc::CtmcResult<double> r = line::ctmc::solver_ctmc(sn, reach, sy);
    for (std::size_t i = 0; i < r.Q.rows(); ++i) {
        double s = 0;
        for (std::size_t j = 0; j < r.Q.cols(); ++j) s += r.Q(i, j);
        CHECK(s == doctest::Approx(0.0).epsilon(1e-12));
    }
    const std::vector<double> pi = line::mc::ctmc_solve(r.Q);
    const line::ctmc::CtmcAvg<double> avg = line::ctmc::solver_ctmc_avg_from_pi(sn, r, pi);
    // Closed network: the two stations hold the whole population between them.
    CHECK(avg.QN(0, 0) + avg.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
    // A closed chain has one throughput, seen identically at both stations.
    CHECK(avg.TN(0, 0) == doctest::Approx(avg.TN(1, 0)).epsilon(1e-9));
}

TEST_CASE("an SPN transition fires atomically across its arcs") {
    // P1 -> T1 -> P2, one token, one mode: the classic marked graph.
    qn::Network<double> m("spn");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(p2, c, line::lang::Distrib<double>::exp_rate(1.0));

    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    // enabling/firing are indexed by NODE: consume one from P1, produce to P2.
    tp.enabling.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(3, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.enabling[0](p1 - 1, 0) = 1.0;
    tp.firing[0](p2 - 1, 0) = 1.0;
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(line::lang::TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(line::lang::Distrib<double>::exp_rate(2.0));
    const std::size_t tr = m.add_transition("T1", tp);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<qn::GlobalSync<double>> gs = qn::refresh_global_sync(sn);
    // One ENABLE and one FIRE for the single mode.
    REQUIRE(gs.size() == 2);
    CHECK(gs[0].active.event == EventType::ENABLE);
    CHECK(gs[1].active.event == EventType::FIRE);
    CHECK(gs[1].active.node == tr);

    // The FIRE carries a PRE on P1 and a POST on P2: the arcs move together.
    std::size_t npre = 0, npost = 0;
    for (std::size_t j = 0; j < gs[1].passive.size(); ++j) {
        if (gs[1].passive[j].event == EventType::PRE) {
            ++npre;
            CHECK(gs[1].passive[j].node == p1);
        }
        if (gs[1].passive[j].event == EventType::POST) {
            ++npost;
            CHECK(gs[1].passive[j].node == p2);
        }
    }
    CHECK(npre == 1);
    CHECK(npost == 1);
    // An ENABLE only READS the marking, so it moves no token at all.
    for (std::size_t j = 0; j < gs[0].passive.size(); ++j)
        CHECK(gs[0].passive[j].event == EventType::LOCAL);
}

TEST_CASE("an SPN firing moves the token and is marked a completion") {
    qn::Network<double> m("spn2");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(p2, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(3, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.enabling[0](p1 - 1, 0) = 1.0;
    tp.firing[0](p2 - 1, 0) = 1.0;
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(line::lang::TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(line::lang::Distrib<double>::exp_rate(2.0));
    m.add_transition("T1", tp);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::GlobalSync<double>> gs = qn::refresh_global_sync(sn);

    // Build a state with the token at P1 and one server of the mode running.
    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{0});
    REQUIRE(!ss.empty());
    qn::NetState<double> st = ss[0];
    // Stateful order is P1, P2, T1. Put the token at P1 and run the server.
    st.local[0] = std::vector<double>{1.0};
    st.local[1] = std::vector<double>{0.0};
    st.local[2] = std::vector<double>{0.0, 1.0, 0.0};  // idle=0, phase1=1, fired=0

    const qn::GlobalOutcome<double> go = qn::after_global_event(sn, st, gs[1]);
    REQUIRE(go.space.size() == 1);
    // The token left P1 and arrived at P2 in ONE transition: there is no state
    // in which it has left one place without reaching the other.
    CHECK(go.space[0].local[0][0] == doctest::Approx(0.0));
    CHECK(go.space[0].local[1][0] == doctest::Approx(1.0));
    // The firing server returned to the idle pool.
    CHECK(go.space[0].local[2][0] == doctest::Approx(1.0));
    CHECK(go.space[0].local[2][1] == doctest::Approx(0.0));
    CHECK(go.rate[0] == doctest::Approx(2.0));
    // The outcome is a COMPLETION. Callers must not re-derive this from the
    // markings: a firing that returns exactly what it consumed leaves every
    // marking invariant yet still completed.
    REQUIRE(go.completion.size() == 1);
    CHECK(go.completion[0]);
}

TEST_CASE("an inhibited SPN mode cannot fire") {
    qn::Network<double> m("spn3");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, line::lang::Distrib<double>::exp_rate(1.0));
    m.set_service(p2, c, line::lang::Distrib<double>::exp_rate(1.0));
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(3, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(3, 1, 0.0));
    tp.enabling[0](p1 - 1, 0) = 1.0;
    tp.firing[0](p2 - 1, 0) = 1.0;
    tp.inhibiting[0](p1 - 1, 0) = 1.0;  // inhibited as soon as P1 holds a token
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(line::lang::TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(line::lang::Distrib<double>::exp_rate(2.0));
    m.add_transition("T1", tp);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::GlobalSync<double>> gs = qn::refresh_global_sync(sn);

    const std::vector<qn::NetState<double>> ss =
        qn::space_generator(sn, std::vector<std::size_t>{0});
    qn::NetState<double> st = ss[0];
    st.local[0] = std::vector<double>{1.0};
    st.local[1] = std::vector<double>{0.0};
    st.local[2] = std::vector<double>{0.0, 1.0, 0.0};

    // The inhibitor threshold is met, so the mode is disabled and no firing
    // transition exists -- the token stays put.
    CHECK(qn::after_global_event(sn, st, gs[1]).empty());
}
