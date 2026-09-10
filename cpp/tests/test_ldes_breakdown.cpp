/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Server breakdowns in the native LDES engine.
 *
 * THE ORACLE IS AN EXACT CTMC BUILT IN THIS FILE, not a reference value from
 * another codebase and not a simulation. A finite-buffer M/M/1/K whose server
 * fails and is repaired on independent exponential clocks is a small,
 * two-dimensional chain (jobs, server status) that can be written down and
 * solved by `ctmc_solve`, so the engine is checked against arithmetic rather
 * than against agreement. That matters here more than usual: a breakdown that
 * was implemented as an EVICTION rather than an interruption, or whose failure
 * clock ran only while the server was busy, would still produce a stable
 * simulation with plausible utilizations.
 *
 *   states (n, s), n = 0..K jobs, s = 1 up / 0 down
 *   (n,s) -> (n+1,s)   lambda        n < K
 *   (n,1) -> (n-1,1)   mu            n >= 1
 *   (n,0) -> (n-1,0)   mu_down       n >= 1, when a degraded rate is declared
 *   (n,1) -> (n,0)     xi            the server fails whether or not it serves
 *   (n,0) -> (n,1)     eta
 *
 * THE FAILURE CLOCK RUNS WHILE THE SERVER IS IDLE, which is the reference's
 * semantics -- `State.afterEventStation`'s FAILURE branch says so in as many
 * words ("an up server fails at the memoryless rate breakdownMu, whether or not
 * it is serving"). Modelling it as operation-dependent instead would change the
 * availability from eta/(xi+eta) to something load-dependent, so the two are
 * separated by the idle-server availability check below.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ldes/ldes_engine.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

/** M/M/1/K with a breaking server, as the LDES engine sees it. */
qn::Network<double> mm1k_breakdown(double lambda, double mu, std::size_t K, double xi, double eta,
                                   double mu_down) {
    qn::Network<double> m("mm1k-breakdown");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, D::exp_rate(lambda));
    m.set_service(q, c, D::exp_rate(mu));
    m.set_capacity(q, static_cast<double>(K));
    std::vector<D> down;
    if (mu_down > 0.0) down.push_back(D::exp_rate(mu_down));
    m.set_breakdown(q, D::exp_rate(xi), D::exp_rate(eta), down);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    return m;
}

/** What the exact chain says about the same model. */
struct Exact {
    double qlen = 0.0;   ///< mean number in system at the queue
    double tput = 0.0;   ///< departure rate
    double avail = 0.0;  ///< P(server up)
    /**
     * P(a job is held in the server), which is what the engine reports as
     * utilization -- NOT the fraction of time the server is working.
     *
     * The two differ by exactly the outages: a job holds the slot across a
     * failure (it is interrupted, not evicted), so the server is occupied and
     * idle at the same time. That is why `U = T E[S]` does NOT hold on a
     * breakdown model, and it is not a defect: on the fixture below the
     * occupancy is 0.653009 and T E[S] is 0.389074.
     */
    double occupancy = 0.0;
};

Exact exact_mm1k_breakdown(double lambda, double mu, std::size_t K, double xi, double eta,
                           double mu_down) {
    // Index (n, s) as 2n + s, so state 0 is (0 jobs, down) and 1 is (0, up).
    const std::size_t n = 2 * (K + 1);
    Matrix<double> Q(n, n, 0.0);
    auto idx = [](std::size_t jobs, std::size_t up) { return 2 * jobs + up; };
    for (std::size_t j = 0; j <= K; ++j)
        for (std::size_t s = 0; s <= 1; ++s) {
            const std::size_t a = idx(j, s);
            if (j < K) Q(a, idx(j + 1, s)) += lambda;
            if (j >= 1) {
                const double rate = (s == 1) ? mu : mu_down;
                if (rate > 0.0) Q(a, idx(j - 1, s)) += rate;
            }
            if (s == 1)
                Q(a, idx(j, 0)) += xi;
            else
                Q(a, idx(j, 1)) += eta;
        }
    for (std::size_t a = 0; a < n; ++a) {
        double off = 0.0;
        for (std::size_t b = 0; b < n; ++b)
            if (a != b) off += Q(a, b);
        Q(a, a) = -off;
    }
    const std::vector<double> pi = mc::ctmc_solve(Q);

    Exact e;
    for (std::size_t j = 0; j <= K; ++j)
        for (std::size_t s = 0; s <= 1; ++s) {
            const double p = pi[idx(j, s)];
            e.qlen += static_cast<double>(j) * p;
            if (s == 1) e.avail += p;
            if (j >= 1) {
                e.tput += p * ((s == 1) ? mu : mu_down);
                e.occupancy += p;
            }
        }
    return e;
}

ldes::LdesOptions opts(std::size_t samples, long seed = 23000) {
    ldes::LdesOptions o;
    o.samples = samples;
    o.seed = seed;
    return o;
}

}  // namespace

TEST_CASE("LDES breakdown: a stopped server matches the exact chain") {
    // The server is down a third of the time and delivers no service there, so
    // the queue is far longer than the same model without breakdowns -- which
    // is what makes this discriminating rather than a tolerance check.
    const double lambda = 0.4, mu = 1.0, xi = 0.2, eta = 0.4;
    const std::size_t K = 8;
    const Exact ex = exact_mm1k_breakdown(lambda, mu, K, xi, eta, 0.0);
    CHECK(ex.avail == doctest::Approx(eta / (xi + eta)).epsilon(1e-12));

    qn::Network<double> m = mm1k_breakdown(lambda, mu, K, xi, eta, 0.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
    CHECK(r.QN(1, 0) == doctest::Approx(ex.qlen).epsilon(0.05));
    CHECK(r.TN(1, 0) == doctest::Approx(ex.tput).epsilon(0.03));
    // THE UTILIZATION IS THE OCCUPANCY, and the gap to T E[S] is the outage.
    // Asserting T E[S] here instead would be asserting that the interrupted job
    // gives its server back, which is the opposite of the semantics.
    CHECK(r.UN(1, 0) == doctest::Approx(ex.occupancy).epsilon(0.03));
    CHECK(ex.occupancy > 1.5 * ex.tput);

    // The same model with a server that never fails is a DIFFERENT answer, and
    // by a wide margin: without this the test would pass on an engine that
    // ignored the breakdown entirely.
    const Exact nofail = exact_mm1k_breakdown(lambda, mu, K, 0.0, 1.0, 0.0);
    CHECK(nofail.qlen < 0.5 * ex.qlen);
}

TEST_CASE("LDES breakdown: a DEGRADED server is served, not stopped") {
    // downServiceRates is a slower server and not a stopped one, so the answer
    // must sit strictly between the no-failure and the full-outage models.
    const double lambda = 0.4, mu = 1.0, xi = 0.2, eta = 0.4, mu_down = 0.3;
    const std::size_t K = 8;
    const Exact ex = exact_mm1k_breakdown(lambda, mu, K, xi, eta, mu_down);
    const Exact stopped = exact_mm1k_breakdown(lambda, mu, K, xi, eta, 0.0);
    const Exact nofail = exact_mm1k_breakdown(lambda, mu, K, 0.0, 1.0, 0.0);
    REQUIRE(nofail.qlen < ex.qlen);
    REQUIRE(ex.qlen < stopped.qlen);

    qn::Network<double> m = mm1k_breakdown(lambda, mu, K, xi, eta, mu_down);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
    CHECK(r.QN(1, 0) == doctest::Approx(ex.qlen).epsilon(0.05));
    CHECK(r.TN(1, 0) == doctest::Approx(ex.tput).epsilon(0.03));
    // Strictly between the two, by more than the estimator's noise.
    CHECK(r.QN(1, 0) > nofail.qlen);
    CHECK(r.QN(1, 0) < stopped.qlen);
}

TEST_CASE("LDES breakdown: the outage INTERRUPTS the job, it does not evict it") {
    // A single job cycling a closed Delay -> Queue loop. With no eviction the
    // throughput is 1 / (E[think] + E[completion time]), where the completion
    // time of an exponential service under time-dependent failures is its own
    // exact chain -- built above. An engine that discarded the interrupted job
    // and redrew its service would report a HIGHER throughput, because the
    // exponential redraw discards the work already done for free.
    qn::Network<double> m("closed-breakdown");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Jobs", 1.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.0));
    m.set_breakdown(q, D::exp_rate(0.5), D::exp_rate(0.5));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    // The exact chain of the whole cycle: (where, server) with where = think or
    // queue. Think -> Queue at rate 1, Queue -> Think at mu when up and never
    // when down, plus the failure and repair clocks.
    Matrix<double> Q(4, 4, 0.0);
    // 0 = (think, up), 1 = (think, down), 2 = (queue, up), 3 = (queue, down)
    Q(0, 2) = 1.0;  Q(0, 1) = 0.5;
    Q(1, 3) = 1.0;  Q(1, 0) = 0.5;
    Q(2, 0) = 1.0;  Q(2, 3) = 0.5;
    Q(3, 2) = 0.5;
    for (std::size_t a = 0; a < 4; ++a) {
        double off = 0.0;
        for (std::size_t b = 0; b < 4; ++b)
            if (a != b) off += Q(a, b);
        Q(a, a) = -off;
    }
    const std::vector<double> pi = mc::ctmc_solve(Q);
    const double xput = pi[2] * 1.0;  // completions happen only from (queue, up)
    CHECK(r.TN(1, 0) == doctest::Approx(xput).epsilon(0.04));
    // Population conservation: the one job is either thinking or queueing.
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE("LDES breakdown: refused under slotted time, and both clocks are required") {
    qn::Network<double> m = mm1k_breakdown(0.4, 1.0, 4, 0.2, 0.4, 0.0);
    ldes::LdesOptions o = opts(1000);
    o.slotted = true;
    o.slot_length = 1.0;
    // The failure and repair epochs are drawn on a continuous clock and cannot
    // land on the lattice, so the run is refused rather than rounded onto it.
    CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), o), UnsupportedError);

    // A server that fails and is never repaired is an absorbing model, not a
    // breakdown, and the builder declines to infer an infinite repair time.
    qn::Network<double> b("bad");
    b.add_source("S");
    const std::size_t q = b.add_queue("Q", lang::SchedStrategy::FCFS);
    // A CLASS IS NEEDED FOR THE DEGRADED-SERVICE CHECK TO RUN AT ALL: the
    // per-class loop is bounded by the class count, so on a classless model the
    // Erlang below is never looked at and the refusal never fires.
    b.add_open_class("C");
    CHECK_THROWS_AS(b.set_breakdown(q, D::exp_rate(0.5), D::disabled_dist()), InputError);
    // A phase-type degraded service would need its own phase block in the joint
    // chain, which no codebase builds, so it is refused rather than approximated.
    std::vector<D> erl(1, D::erlang(2.0, 2));
    CHECK_THROWS_AS(b.set_breakdown(q, D::exp_rate(0.5), D::exp_rate(0.5), erl), UnsupportedError);
}
