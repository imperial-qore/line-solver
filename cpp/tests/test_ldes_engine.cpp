/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The NATIVE LDES engine, first increment.
 *
 * WHY THE ASSERTIONS ARE WHAT THEY ARE. This engine does not share SSJ's
 * random streams, so a fixed seed gives a DIFFERENT sample path from the Java
 * engine and no number here may be compared element by element against it.
 * What a simulator can be held to is:
 *
 *  1. THE CLOSED FORM, within a tolerance derived from the completion budget.
 *     M/M/1 at rho gives QLen = rho/(1-rho) and RespT = 1/(mu-lambda); a run
 *     of n completions estimates them with a relative standard error that
 *     grows like 1/((1-rho) sqrt(n)), so the tolerance is quoted per case and
 *     not shared.
 *  2. THE LAWS, which hold on any path however short: utilization at most one,
 *     flow balance between the source and the sink, Little on the truncated
 *     table up to the warmup mass the reference's response-time tally keeps.
 *  3. REPRODUCIBILITY. One seed, two runs, the same numbers exactly -- this is
 *     an equality check, and it is the one that catches an uninitialized read
 *     or a container iterated in address order.
 *  4. THE REFUSALS. Everything outside the increment must throw by name rather
 *     than simulate a different model; a PS station served FCFS would produce
 *     numbers no downstream test could tell from correct ones.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/api/mam/map_moment.h"
#include "line/solvers/ldes/ldes_engine.h"

using namespace line;

namespace {

/** M/M/1 with arrival rate `lambda` and service rate `mu`. */
qn::Network<double> mm1(double lambda, double mu) {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    return m;
}

ldes::LdesOptions opts(std::size_t samples, long seed = 23000) {
    ldes::LdesOptions o;
    o.samples = samples;
    o.seed = seed;
    return o;
}

}  // namespace

TEST_CASE("native LDES engine: M/M/1 matches the closed form") {
    // rho = 0.5: QLen = 1, RespT = 2, Util = 0.5, Tput = 0.5.
    qn::Network<double> m = mm1(0.5, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));

    REQUIRE(r.nstations == 2);
    REQUIRE(r.nclasses == 1);
    // Station 0 is the Source, station 1 the Queue.
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.05));
    CHECK(r.RN(1, 0) == doctest::Approx(2.0).epsilon(0.05));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(0.02));
    CHECK(r.TN(1, 0) == doctest::Approx(0.5).epsilon(0.02));
    // The Source row carries the nominal arrival rate, as the reference does.
    CHECK(r.TN(0, 0) == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(r.XN(0, 0) == doctest::Approx(0.5).epsilon(0.02));
}

TEST_CASE("native LDES engine: M/M/1 at a heavier load") {
    // rho = 0.8: QLen = 4, RespT = 5. The estimator's variance grows like
    // 1/(1-rho)^2, so the tolerance is looser at the same budget -- that is a
    // property of the estimator and not slack in the test.
    qn::Network<double> m = mm1(0.8, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
    CHECK(r.QN(1, 0) == doctest::Approx(4.0).epsilon(0.08));
    CHECK(r.RN(1, 0) == doctest::Approx(5.0).epsilon(0.08));
    CHECK(r.UN(1, 0) == doctest::Approx(0.8).epsilon(0.02));
    CHECK(r.TN(1, 0) == doctest::Approx(0.8).epsilon(0.02));
}

TEST_CASE("native LDES engine: the laws hold on the returned table") {
    qn::Network<double> m = mm1(0.7, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    CHECK(r.UN(1, 0) <= 1.0);
    // Little on the truncated pair: QLen = Tput * RespT. RespT averages the
    // warmup observations too (the reference's tally is not truncated), so the
    // identity holds to the warmup mass, not to machine precision.
    CHECK(r.QN(1, 0) == doctest::Approx(r.TN(1, 0) * r.RN(1, 0)).epsilon(0.05));
    // Flow balance: what the queue completes is what the source released.
    CHECK(r.TN(1, 0) == doctest::Approx(r.TN(0, 0)).epsilon(0.02));
}

TEST_CASE("native LDES engine: M/M/c") {
    // lambda = 1.6, mu = 1, c = 2 gives rho = 0.8 per server; Erlang-C puts
    // QLen at 4.4444... and the per-server utilization at 0.8.
    qn::Network<double> m("mmc");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(1.6));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_number_of_servers(q, 2);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
    CHECK(r.UN(1, 0) == doctest::Approx(0.8).epsilon(0.02));
    CHECK(r.QN(1, 0) == doctest::Approx(4.4444444).epsilon(0.08));
    CHECK(r.TN(1, 0) == doctest::Approx(1.6).epsilon(0.02));
}

TEST_CASE("native LDES engine: a two-station open tandem") {
    qn::Network<double> m("tandem");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
    // Jackson: each station is an independent M/M/1 at its own rho.
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.06));       // rho = 0.5
    CHECK(r.QN(2, 0) == doctest::Approx(1.0 / 3.0).epsilon(0.06));  // rho = 0.25
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(0.02));
    CHECK(r.UN(2, 0) == doctest::Approx(0.25).epsilon(0.02));
    // The system response time is the sum of the two station response times.
    CHECK(r.CN(0, 0) == doctest::Approx(2.0 + 2.0 / 3.0).epsilon(0.06));
}

TEST_CASE("native LDES engine: a closed cyclic network") {
    // Delay(rate 1) -> Queue(rate 2) -> Delay, N = 2. Exact by MVA: R(1) = 0.5
    // and X(1) = 1/1.5, so Q(1) = 1/3; then R(2) = 0.5(1 + 1/3) = 2/3,
    // X(2) = 2/(1 + 2/3) = 1.2, Q_queue = 0.8, U_queue = 0.6, Q_delay = 1.2.
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2, d);
    m.set_service(d, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    // The population is conserved: the two queue lengths must sum to N.
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(0.03));
    CHECK(r.TN(1, 0) == doctest::Approx(1.2).epsilon(0.03));
    CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(0.03));
    CHECK(r.QN(1, 0) == doctest::Approx(0.8).epsilon(0.05));
    CHECK(r.QN(0, 0) == doctest::Approx(1.2).epsilon(0.05));
}

TEST_CASE("native LDES engine: one seed gives one path") {
    qn::Network<double> m = mm1(0.6, 1.0);
    const ldes::LdesResult a = ldes::ldes_engine_solve(m.get_struct(), opts(50000, 777));
    const ldes::LdesResult b = ldes::ldes_engine_solve(m.get_struct(), opts(50000, 777));
    CHECK(a.QN(1, 0) == b.QN(1, 0));
    CHECK(a.RN(1, 0) == b.RN(1, 0));
    CHECK(a.UN(1, 0) == b.UN(1, 0));
    CHECK(a.TN(1, 0) == b.TN(1, 0));
    // A different seed must move the estimate, or the stream is not being used.
    const ldes::LdesResult c = ldes::ldes_engine_solve(m.get_struct(), opts(50000, 778));
    CHECK(a.QN(1, 0) != c.QN(1, 0));
}

TEST_CASE("native LDES engine: MSER-5 truncation is the reference's rule") {
    // A constant series has zero batch variance everywhere, so the criterion is
    // minimised at d = 0: no truncation. A series that starts high and settles
    // must truncate past the transient.
    std::vector<double> flat(200, 1.0);
    CHECK(ldes::engine::mser5_truncation(flat, 5) == 0);

    std::vector<double> ramp;
    for (int i = 0; i < 50; ++i) ramp.push_back(100.0 - 2.0 * i);
    for (int i = 0; i < 250; ++i) ramp.push_back(1.0);
    CHECK(ldes::engine::mser5_truncation(ramp, 5) >= 10);

    // Below four batches the rule has too few terms and reports none.
    std::vector<double> tiny(15, 3.0);
    CHECK(ldes::engine::mser5_truncation(tiny, 5) == 0);
}

// ===========================================================================
// Increment 2: distributions, the discipline families, capacity, CI.
// ===========================================================================

namespace {

/** M/G/1 with an arbitrary service law, so the closed forms below can be used. */
qn::Network<double> mg1(double lambda, const lang::Distrib<double>& svc,
                        lang::SchedStrategy sched = lang::SchedStrategy::FCFS) {
    qn::Network<double> m("mg1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q, c, svc);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    return m;
}

/** Pollaczek-Khinchine: the mean number in an M/G/1 queue. */
double pk_qlen(double rho, double scv) {
    return rho + rho * rho * (1.0 + scv) / (2.0 * (1.0 - rho));
}

/** Two-phase hyperexponential with balanced means, given (mean, SCV > 1). */
lang::Distrib<double> hyperexp_moments(double mean, double scv) {
    const double p = 0.5 * (1.0 + std::sqrt((scv - 1.0) / (scv + 1.0)));
    return lang::Distrib<double>::hyperexp(p, 2.0 * p / mean, 2.0 * (1.0 - p) / mean);
}

/** Gamma given (mean, SCV): shape = 1/SCV, scale = mean*SCV. */
lang::Distrib<double> gamma_moments(double mean, double scv) {
    return lang::Distrib<double>::gamma_dist(1.0 / scv, mean * scv);
}

/** Lognormal given (mean, SCV), through its log-moments. */
lang::Distrib<double> lognormal_moments(double mean, double scv) {
    const double sigma = std::sqrt(std::log(1.0 + scv));
    return lang::Distrib<double>::lognormal(std::log(mean) - 0.5 * sigma * sigma, sigma);
}

}  // namespace

TEST_CASE("native LDES engine: M/G/1 against Pollaczek-Khinchine") {
    // One law per shape family, each at rho = 0.5 and mean service 0.5 so the
    // only thing that moves the answer is the SCV. P-K is exact for FCFS.
    const double lambda = 1.0, mean = 0.5, rho = 0.5;

    SUBCASE("Erlang-4, SCV = 0.25") {
        // Erlang(k phases, phase rate): mean = k/rate, SCV = 1/k.
        qn::Network<double> m = mg1(lambda, lang::Distrib<double>::erlang(8.0, 4));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 0.25)).epsilon(0.05));
        CHECK(r.UN(1, 0) == doctest::Approx(rho).epsilon(0.02));
    }
    SUBCASE("Deterministic, SCV = 0") {
        qn::Network<double> m = mg1(lambda, lang::Distrib<double>::det(mean));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 0.0)).epsilon(0.05));
        CHECK(r.UN(1, 0) == doctest::Approx(rho).epsilon(0.02));
    }
    SUBCASE("Hyperexponential, SCV = 4") {
        qn::Network<double> m = mg1(lambda, hyperexp_moments(mean, 4.0));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(2000000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 4.0)).epsilon(0.10));
        CHECK(r.UN(1, 0) == doctest::Approx(rho).epsilon(0.03));
    }
    SUBCASE("Gamma, SCV = 2") {
        qn::Network<double> m = mg1(lambda, gamma_moments(mean, 2.0));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 2.0)).epsilon(0.08));
    }
    SUBCASE("Lognormal, SCV = 2") {
        qn::Network<double> m = mg1(lambda, lognormal_moments(mean, 2.0));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 2.0)).epsilon(0.10));
    }
    SUBCASE("Uniform, SCV = 1/12") {
        qn::Network<double> m = mg1(lambda, lang::Distrib<double>::uniform(0.25, 0.75));
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
        CHECK(r.QN(1, 0) == doctest::Approx(pk_qlen(rho, 1.0 / 12.0)).epsilon(0.05));
    }
}

TEST_CASE("native LDES engine: the work-conserving disciplines share one mean") {
    // FCFS, LCFS, SIRO, SJF, LJF, SEPT and LEPT are all work-conserving and
    // non-preemptive, so on a SINGLE class they have the same mean queue
    // length and differ only in its variance. A discipline wired to the wrong
    // comparator still passes this; what it catches is a discipline that
    // LOSES or DUPLICATES work, which every wrong container does.
    // SJF and LJF are NOT in this list, and that is the point: they order by
    // the sampled SIZE, so they change the mean number in system even though
    // they conserve work. Only the orders that are independent of the service
    // time share the M/M/1 mean -- SEPT and LEPT among them, because a single
    // class makes their class-mean key constant and they degenerate to FCFS.
    const lang::SchedStrategy scheds[] = {
        lang::SchedStrategy::FCFS, lang::SchedStrategy::LCFS, lang::SchedStrategy::SIRO,
        lang::SchedStrategy::SEPT, lang::SchedStrategy::LEPT, lang::SchedStrategy::HOL};
    for (lang::SchedStrategy s : scheds) {
        qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0), s);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
        INFO("scheduling = " << lang::sched_to_text(s));
        CHECK(r.QN(1, 0) == doctest::Approx(1.5).epsilon(0.08));   // rho/(1-rho)
        CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(0.02));
        CHECK(r.TN(1, 0) == doctest::Approx(0.6).epsilon(0.02));
    }
}

TEST_CASE("native LDES engine: SJF beats FCFS beats LJF") {
    // The size-based orders are the ones whose mean DOES move, and in a fixed
    // direction: serving the shortest job first minimises the mean number in
    // system among non-preemptive orders, serving the longest first maximises
    // it. All three still deliver the same work, so the utilization and the
    // throughput are common to them.
    double q[3] = {0.0, 0.0, 0.0};
    const lang::SchedStrategy order[3] = {lang::SchedStrategy::SJF, lang::SchedStrategy::FCFS,
                                          lang::SchedStrategy::LJF};
    for (int k = 0; k < 3; ++k) {
        qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0), order[k]);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
        q[k] = r.QN(1, 0);
        CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(0.02));
        CHECK(r.TN(1, 0) == doctest::Approx(0.6).epsilon(0.02));
    }
    CHECK(q[0] < q[1]);
    CHECK(q[1] < q[2]);
}

TEST_CASE("native LDES engine: M/M/1-PS equals M/M/1-FCFS") {
    // PS is insensitive and, for exponential service, has the SAME queue
    // length distribution as FCFS. This is the check that the share
    // integration and the reschedule-on-every-change actually conserve work:
    // an error in either shows up immediately as a wrong mean.
    qn::Network<double> m = mg1(0.7, lang::Distrib<double>::exp_rate(1.0),
                                lang::SchedStrategy::PS);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.7 / 0.3).epsilon(0.08));
    CHECK(r.RN(1, 0) == doctest::Approx(1.0 / 0.3).epsilon(0.08));
    CHECK(r.UN(1, 0) == doctest::Approx(0.7).epsilon(0.02));
}

TEST_CASE("native LDES engine: M/G/1-PS is insensitive to the service law") {
    // The point of PS: at the same mean, a hyperexponential service (SCV 4)
    // gives the SAME queue length as an exponential one, where FCFS would give
    // P-K's 2.5x. A share rule that used the residual work of the wrong job
    // would break this and nothing else in the suite would notice.
    qn::Network<double> m = mg1(1.0, hyperexp_moments(0.5, 4.0),
                                lang::SchedStrategy::PS);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(2000000));
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.10));  // rho/(1-rho) at rho = 0.5
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(0.03));
}

TEST_CASE("native LDES engine: DPS with equal weights is PS") {
    qn::Network<double> m = mg1(0.7, lang::Distrib<double>::exp_rate(1.0),
                                lang::SchedStrategy::DPS);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.7 / 0.3).epsilon(0.10));
}

TEST_CASE("native LDES engine: DPS gives the heavier weight the shorter wait") {
    // Two classes, equal load and equal service, weights 1 and 4. DPS is
    // work-conserving so the TOTAL is PS's, but the weighted class must be
    // strictly faster -- that ordering is the whole content of the discipline.
    qn::Network<double> m("dps2");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::DPS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Light");
    const std::size_t c2 = m.add_open_class("Heavy");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(0.35));
    m.set_arrival(src, c2, lang::Distrib<double>::exp_rate(0.35));
    m.set_service(q, c1, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c2, lang::Distrib<double>::exp_rate(1.0));
    m.set_sched_param(q, c1, 1.0);
    m.set_sched_param(q, c2, 4.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(800000));
    CHECK(r.RN(1, 1) < r.RN(1, 0));
    // Work conservation: the total is what plain PS would hold at rho = 0.7.
    CHECK(r.QN(1, 0) + r.QN(1, 1) == doctest::Approx(0.7 / 0.3).epsilon(0.12));
}

TEST_CASE("native LDES engine: HOL serves the urgent class first") {
    // Non-preemptive priority, equal loads. The class with the LOWER priority
    // value is more urgent (LINE's ascending convention) and must wait less.
    qn::Network<double> m("hol");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::HOL);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t hi = m.add_open_class("Urgent", 0);
    const std::size_t lo = m.add_open_class("Background", 5);
    m.set_arrival(src, hi, lang::Distrib<double>::exp_rate(0.35));
    m.set_arrival(src, lo, lang::Distrib<double>::exp_rate(0.35));
    m.set_service(q, hi, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, lo, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(hi, hi, src, q, 1.0);
    P.set(hi, hi, q, snk, 1.0);
    P.set(lo, lo, src, q, 1.0);
    P.set(lo, lo, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(800000));
    CHECK(r.RN(1, 0) < r.RN(1, 1));
    // Cobham for the head of the line: W_1 = rho E[S^2] / (2 (1 - rho_1)) with
    // exponential service, E[S^2] = 2, rho = 0.7, rho_1 = 0.35 -> W_1 = 1.0769,
    // so R_1 = W_1 + 1.
    CHECK(r.RN(1, 0) == doctest::Approx(1.0 + 0.7 * 2.0 / (2.0 * (1.0 - 0.35))).epsilon(0.08));
    // Work conservation: the load-weighted mean wait is FCFS's.
    const double w = 0.5 * (r.RN(1, 0) - 1.0) + 0.5 * (r.RN(1, 1) - 1.0);
    CHECK(w == doctest::Approx(0.7 / (1.0 - 0.7)).epsilon(0.10));
}

TEST_CASE("native LDES engine: M/M/1/K drops the overflow") {
    // K = 3 (Kendall: waiting plus in service). The exact loss probability of
    // M/M/1/K at rho is (1-rho) rho^K / (1 - rho^(K+1)), and the carried
    // throughput is lambda (1 - P_K).
    const double rho = 0.8, lambda = 0.8, mu = 1.0;
    const int Kcap = 3;
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(mu));
    m.set_capacity(q, Kcap);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    double num = 0.0, den = 0.0;
    for (int n = 0; n <= Kcap; ++n) den += std::pow(rho, n);
    num = std::pow(rho, Kcap);
    const double pK = num / den;
    double qlen = 0.0;
    for (int n = 0; n <= Kcap; ++n) qlen += n * std::pow(rho, n) / den;
    CHECK(r.QN(1, 0) == doctest::Approx(qlen).epsilon(0.04));
    CHECK(r.TN(1, 0) == doctest::Approx(lambda * (1.0 - pK)).epsilon(0.04));
    CHECK(r.QN(1, 0) <= static_cast<double>(Kcap));
}

TEST_CASE("native LDES engine: confidence intervals are reported and shrink") {
    // A half-width is only meaningful if it falls as the budget grows. Four
    // times the samples should roughly halve it; the check is loose (it must
    // merely shrink) because the ratio is itself an estimate.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    const ldes::LdesResult a = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(m.get_struct(), opts(800000));
    REQUIRE(a.QNCI.rows() == 2);
    CHECK(a.QNCI(1, 0) > 0.0);
    CHECK(a.TNCI(1, 0) > 0.0);
    CHECK(b.QNCI(1, 0) < a.QNCI(1, 0));
    // The interval must cover the exact value at 95%.
    CHECK(std::abs(b.QN(1, 0) - 1.0) < 3.0 * b.QNCI(1, 0));
}

TEST_CASE("native LDES engine: a MAP service keeps its autocorrelation") {
    // An MMPP(2) service with the SAME mean as an exponential one produces a
    // markedly longer queue, because its interdeparture times are positively
    // correlated. If the sampler restarted the phase every variate the process
    // would collapse to a renewal one with the right marginal, and the queue
    // would fall back towards the exponential answer -- which is exactly the
    // failure this case exists to catch.
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;  D0(0, 1) = 0.1;
    D0(1, 0) = 0.1;   D0(1, 1) = -0.4;
    D1(0, 0) = 1.9;
    D1(1, 1) = 0.3;
    lang::Distrib<double> map = lang::Distrib<double>::map_dist(D0, D1, lang::ProcessType::MMPP2);
    // `map_dist` leaves the moment fields as placeholders (mean 0, SCV 1), so
    // the reference values have to come from the pair itself.
    mam::Map<double> pair;
    pair.D0 = D0;
    pair.D1 = D1;
    const double map_mean = mam::map_mean(pair);
    const double map_scv = mam::map_scv(pair);

    qn::Network<double> m = mg1(0.2, map);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
    // Sanity first: the station is stable and the utilization is lambda*E[S].
    CHECK(r.UN(1, 0) < 1.0);
    CHECK(r.UN(1, 0) == doctest::Approx(0.2 * map_mean).epsilon(0.05));
    // The correlated stream must congest more than the P-K value its marginal
    // alone would predict.
    const double rho = r.UN(1, 0);
    CHECK(r.QN(1, 0) > pk_qlen(rho, map_scv));
}

TEST_CASE("native LDES engine: the increment's boundary is refused by name") {
    SUBCASE("an EXT station") {
        // EXT is the one discipline still unported; it must be refused rather
        // than served as something else.
        qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0),
                                    lang::SchedStrategy::EXT);
        CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), opts(1000)), UnsupportedError);
    }
    SUBCASE("a PAS station with no rate function") {
        // A pass-and-swap station IS ported, but it is meaningless without the
        // mu(c) its whole definition rests on, so that is an input error rather
        // than an unsupported feature.
        qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0),
                                    lang::SchedStrategy::PAS);
        CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), opts(1000)), InputError);
    }
}

// ===========================================================================
// Increment 3: preemption and deadlines.
// ===========================================================================

TEST_CASE("native LDES engine: preemptive M/M/1 keeps the M/M/1 mean") {
    // Under EXPONENTIAL service every work-conserving preemptive discipline has
    // the M/M/1 queue length distribution: memorylessness makes resume and
    // restart agree, and the number in system does not depend on the order.
    // This is the check that preemption neither loses nor duplicates work --
    // a stale departure event surviving its preemption shows up here at once.
    const lang::SchedStrategy scheds[] = {
        lang::SchedStrategy::LCFSPR, lang::SchedStrategy::LCFSPI,
        lang::SchedStrategy::FCFSPR, lang::SchedStrategy::FCFSPI,
        lang::SchedStrategy::SRPT,   lang::SchedStrategy::PSJF,
        lang::SchedStrategy::FB,     lang::SchedStrategy::LRPT,
        lang::SchedStrategy::SETF,   lang::SchedStrategy::FSP};
    for (lang::SchedStrategy s : scheds) {
        qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0), s);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
        INFO("scheduling = " << lang::sched_to_text(s));
        // The size-based orders change the RESPONSE TIME but not the
        // utilization or the throughput, which are pure work conservation.
        CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(0.03));
        CHECK(r.TN(1, 0) == doctest::Approx(0.6).epsilon(0.03));
        CHECK(r.UN(1, 0) <= 1.0);
    }
}

TEST_CASE("native LDES engine: the order-independent preemptive disciplines") {
    // LCFSPR, LCFSPI, FCFSPR and FCFSPI do not reorder by size on a single
    // class, so they keep the M/M/1 mean exactly, not merely its work.
    const lang::SchedStrategy scheds[] = {
        lang::SchedStrategy::LCFSPR, lang::SchedStrategy::LCFSPI,
        lang::SchedStrategy::FCFSPR, lang::SchedStrategy::FCFSPI};
    for (lang::SchedStrategy s : scheds) {
        qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0), s);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
        INFO("scheduling = " << lang::sched_to_text(s));
        CHECK(r.QN(1, 0) == doctest::Approx(1.5).epsilon(0.08));
    }
}

TEST_CASE("native LDES engine: SRPT is the optimal order, LRPT the worst") {
    // SRPT minimises the mean number in system over ALL orders, preemptive or
    // not, and LRPT maximises it. Both must therefore straddle non-preemptive
    // FCFS, and SRPT must also beat the non-preemptive SJF, which cannot
    // interrupt a long job already running.
    double q_srpt = 0.0, q_sjf = 0.0, q_fcfs = 0.0, q_lrpt = 0.0;
    {
        qn::Network<double> m = mg1(0.7, hyperexp_moments(1.0, 4.0), lang::SchedStrategy::SRPT);
        q_srpt = ldes::ldes_engine_solve(m.get_struct(), opts(500000)).QN(1, 0);
    }
    {
        qn::Network<double> m = mg1(0.7, hyperexp_moments(1.0, 4.0), lang::SchedStrategy::SJF);
        q_sjf = ldes::ldes_engine_solve(m.get_struct(), opts(500000)).QN(1, 0);
    }
    {
        qn::Network<double> m = mg1(0.7, hyperexp_moments(1.0, 4.0), lang::SchedStrategy::FCFS);
        q_fcfs = ldes::ldes_engine_solve(m.get_struct(), opts(500000)).QN(1, 0);
    }
    {
        qn::Network<double> m = mg1(0.7, hyperexp_moments(1.0, 4.0), lang::SchedStrategy::LRPT);
        q_lrpt = ldes::ldes_engine_solve(m.get_struct(), opts(500000)).QN(1, 0);
    }
    CHECK(q_srpt < q_sjf);
    CHECK(q_sjf < q_fcfs);
    CHECK(q_fcfs < q_lrpt);
    // FCFS is the one with a closed form: P-K at rho = 0.7, SCV = 4.
    CHECK(q_fcfs == doctest::Approx(pk_qlen(0.7, 4.0)).epsilon(0.12));
}

TEST_CASE("native LDES engine: preemptive-resume and restart differ off exponential") {
    // The PR/PI split is invisible under exponential service and visible under
    // any other law: a restarted job discards the work it had done, so the
    // station delivers strictly more work per completion and congests more.
    // An Erlang-4 service (SCV 1/4) makes the loss large enough to see.
    qn::Network<double> mpr = mg1(0.6, lang::Distrib<double>::erlang(4.0, 4),
                                  lang::SchedStrategy::LCFSPR);
    qn::Network<double> mpi = mg1(0.6, lang::Distrib<double>::erlang(4.0, 4),
                                  lang::SchedStrategy::LCFSPI);
    const ldes::LdesResult pr = ldes::ldes_engine_solve(mpr.get_struct(), opts(500000));
    const ldes::LdesResult pi = ldes::ldes_engine_solve(mpi.get_struct(), opts(500000));
    CHECK(pr.QN(1, 0) < pi.QN(1, 0));
    // Preemptive-resume is work-conserving, so its utilization is still lambda
    // times the mean service; restart is NOT, and its utilization is higher.
    CHECK(pr.UN(1, 0) == doctest::Approx(0.6).epsilon(0.03));
    CHECK(pi.UN(1, 0) > pr.UN(1, 0));
}

TEST_CASE("native LDES engine: FCFSPR preempts only a less urgent class") {
    // Preemptive priority. The urgent class should see almost no queueing at
    // all: it preempts on arrival, so its response time approaches its own
    // service time, which non-preemptive HOL cannot deliver.
    qn::Network<double> m("fcfspr");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFSPR);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t hi = m.add_open_class("Urgent", 0);
    const std::size_t lo = m.add_open_class("Background", 5);
    m.set_arrival(src, hi, lang::Distrib<double>::exp_rate(0.3));
    m.set_arrival(src, lo, lang::Distrib<double>::exp_rate(0.4));
    m.set_service(q, hi, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, lo, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(hi, hi, src, q, 1.0);
    P.set(hi, hi, q, snk, 1.0);
    P.set(lo, lo, src, q, 1.0);
    P.set(lo, lo, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(800000));
    // The urgent class sees an M/M/1 in its OWN load alone: the background
    // class is invisible to it because it is always displaced.
    CHECK(r.RN(1, 0) == doctest::Approx(1.0 / (1.0 - 0.3)).epsilon(0.06));
    CHECK(r.RN(1, 1) > r.RN(1, 0));
    // Work conservation across both classes.
    CHECK(r.UN(1, 0) + r.UN(1, 1) == doctest::Approx(0.7).epsilon(0.03));
}

TEST_CASE("native LDES engine: EDF orders by deadline") {
    // Two classes, same load and service, different soft deadlines. EDF must
    // give the tighter deadline the shorter wait; both are work-conserving so
    // the total is the FCFS one.
    qn::Network<double> m("edf");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::EDF);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t tight = m.add_open_class("Tight");
    const std::size_t loose = m.add_open_class("Loose");
    m.set_arrival(src, tight, lang::Distrib<double>::exp_rate(0.35));
    m.set_arrival(src, loose, lang::Distrib<double>::exp_rate(0.35));
    m.set_service(q, tight, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, loose, lang::Distrib<double>::exp_rate(1.0));
    m.set_class_deadline(tight, 1.0);
    m.set_class_deadline(loose, 100.0);
    qn::RoutingMatrix<double> P;
    P.set(tight, tight, src, q, 1.0);
    P.set(tight, tight, q, snk, 1.0);
    P.set(loose, loose, src, q, 1.0);
    P.set(loose, loose, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(600000));
    CHECK(r.RN(1, 0) < r.RN(1, 1));
    CHECK(r.QN(1, 0) + r.QN(1, 1) == doctest::Approx(0.7 / 0.3).epsilon(0.12));
}

// ===========================================================================
// Increment 4: load and class dependence, balking, reneging.
// ===========================================================================

TEST_CASE("native LDES engine: a load-dependent station is M/M/c when the table says so") {
    // lldscaling(n) = min(n, 2) makes one station behave exactly as two
    // servers: the rate doubles from the second job on. The answer must match
    // the M/M/2 case above -- QLen 4.4444 at lambda 1.6, mu 1.
    qn::Network<double> m("ld");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(1.6));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    std::vector<double> alpha;
    for (int n = 1; n <= 200; ++n) alpha.push_back(std::min(n, 2));
    m.set_load_dependence(q, alpha);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(1000000));
    CHECK(r.QN(1, 0) == doctest::Approx(4.4444444).epsilon(0.10));
    CHECK(r.TN(1, 0) == doctest::Approx(1.6).epsilon(0.03));
}

TEST_CASE("native LDES engine: a flat load-dependence table changes nothing") {
    // alpha(n) = 1 for every n is the identity. This is the check that the
    // rescheduling machinery a state-dependent station switches on does not
    // itself perturb the answer -- it re-times every completion on every
    // population change, and a residual mishandled there would show here.
    qn::Network<double> m("ld1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.6));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_load_dependence(q, std::vector<double>(200, 1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.QN(1, 0) == doctest::Approx(1.5).epsilon(0.06));
    CHECK(r.UN(1, 0) == doctest::Approx(0.6).epsilon(0.03));
    CHECK(r.RN(1, 0) == doctest::Approx(2.5).epsilon(0.06));
}

TEST_CASE("native LDES engine: balking refuses to join above a threshold") {
    // Balk with probability 1 whenever 3 or more jobs are present: the station
    // then holds at most 3, exactly as a capacity of 3 would -- and the
    // resulting queue length must match M/M/1/3, whose distribution is the
    // truncated geometric.
    const double rho = 0.8;
    qn::Network<double> m("balk");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.8));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    std::vector<qn::Station<double>::BalkingThreshold> th(1);
    th[0].min_jobs = 3;
    th[0].max_jobs = -1;
    th[0].probability = 1.0;
    m.set_balking(q, c, lang::BalkingStrategy::QUEUE_LENGTH, th);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    double den = 0.0, qlen = 0.0;
    for (int n = 0; n <= 3; ++n) den += std::pow(rho, n);
    for (int n = 0; n <= 3; ++n) qlen += n * std::pow(rho, n) / den;
    CHECK(r.QN(1, 0) == doctest::Approx(qlen).epsilon(0.05));
    CHECK(r.QN(1, 0) <= 3.0);
    REQUIRE(r.balkedCustomers.rows() == 2);
    CHECK(r.balkedCustomers(1, 0) > 0.0);
}

TEST_CASE("native LDES engine: reneging abandons the waiting line") {
    // An M/M/1+M queue: exponential patience of mean 1 at rho = 0.9, which is
    // UNSTABLE without abandonment and stable with it. The queue must stay
    // bounded and the reneged count must be positive; the carried throughput
    // is then strictly below the offered rate.
    qn::Network<double> m("renege");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(1.5));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_patience(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.UN(1, 0) < 1.0);
    CHECK(r.QN(1, 0) < 10.0);              // bounded: without patience it diverges
    CHECK(r.TN(1, 0) < 1.5);               // some load is abandoned, not served
    REQUIRE(r.renegedCustomers.rows() == 2);
    CHECK(r.renegedCustomers(1, 0) > 0.0);
    // Flow balance including abandonment: served + reneged = offered.
    const double abandon_rate = r.renegingRate(1, 0);
    CHECK(r.TN(1, 0) + abandon_rate == doctest::Approx(1.5).epsilon(0.05));
}

// ===========================================================================
// Increment 5: busy periods.
// ===========================================================================

namespace {

/** The measured target whose name matches, or nullptr. */
const ldes::LdesResult::BusyPeriodTarget* bp_named(const ldes::LdesResult& r,
                                                   const std::string& name) {
    for (const auto& t : r.busy_periods)
        if (t.name == name) return &t;
    return 0;
}

}  // namespace

TEST_CASE("native LDES engine: the M/M/1 busy period is 1/(mu-lambda) at every order") {
    // An M/M/1 queue is a birth-death chain with the SAME rates at every level,
    // so an excursion above level n-1 has the law of an excursion above 0: the
    // mean busy period is 1/(mu-lambda) for EVERY order, not just the first.
    // That makes the order axis a real test and not a repetition of one number.
    qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(600000);
    o.busy_period_orders = 3;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    const auto* t = bp_named(r, "Queue");
    REQUIRE(t != 0);
    REQUIRE(t->mean.size() == 3);
    const double exact = 1.0 / (1.0 - 0.6);
    for (std::size_t n = 0; n < 3; ++n) {
        INFO("order " << (n + 1) << " count " << t->count[n]);
        CHECK(t->count[n] > 100.0);
        CHECK(t->mean[n] == doctest::Approx(exact).epsilon(0.08));
    }

    // Renewal-reward consistency: (periods per unit time) x (mean duration) is
    // the fraction of time the station is at level >= 1, i.e. its utilization.
    // This ties the COUNT to the MEAN, which no single-quantity check does.
    const double horizon = t->count[0] / (0.6 * (1.0 - 0.6));
    CHECK(t->count[0] * t->mean[0] / horizon == doctest::Approx(0.6).epsilon(0.10));
}

TEST_CASE("native LDES engine: a per-class busy period target is measured") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(300000);
    o.busy_period_orders = 1;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    // Aggregated and per class, in that order, for every set.
    const auto* agg = bp_named(r, "Queue");
    const auto* per = bp_named(r, "Queue:Class1");
    REQUIRE(agg != 0);
    REQUIRE(per != 0);
    CHECK(per->job_class == 0);
    CHECK(agg->job_class == -1);
    // A single-class model makes the two the same measurement.
    CHECK(per->mean[0] == doctest::Approx(agg->mean[0]).epsilon(1e-12));
}

TEST_CASE("native LDES engine: a same-instant internal move does not end a period") {
    // THE TRAP. In a tandem the subnetwork {Q1,Q2} sees a job leave Q1 and
    // enter Q2 at the SAME instant: the count drops and is restored with no
    // time in between. Treating that drop as an exit ends a busy period that
    // never ended, which roughly HALVES the subnetwork's reported mean and, in
    // particular, makes it come out SHORTER than one of its own members --
    // impossible, since the set is busy whenever either station is.
    qn::Network<double> m("bptandem");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, snk, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(400000);
    o.busy_period_orders = 1;
    std::vector<std::size_t> sub;
    sub.push_back(1);  // Q1, station index 1 (0 is the Source)
    sub.push_back(2);  // Q2
    o.busy_period_subnets.push_back(sub);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    const auto* a = bp_named(r, "Q1");
    const auto* b = bp_named(r, "Q2");
    const auto* both = bp_named(r, "Q1+Q2");
    REQUIRE(a != 0);
    REQUIRE(b != 0);
    REQUIRE(both != 0);
    // Each station alone is an M/M/1 at rho = 0.5: 1/(mu-lambda) = 2.
    CHECK(a->mean[0] == doctest::Approx(2.0).epsilon(0.10));
    CHECK(b->mean[0] == doctest::Approx(2.0).epsilon(0.10));
    // The set is busy whenever EITHER is, so its period cannot be shorter than
    // either member's. This is the assertion the deferral exists to satisfy.
    CHECK(both->mean[0] > a->mean[0]);
    CHECK(both->mean[0] > b->mean[0]);
    // And it must be observed, not merely long: a target that never closes a
    // period reports mean 0 with count 0 and would pass a one-sided check.
    CHECK(both->count[0] > 100.0);
}

TEST_CASE("native LDES engine: no busy period is measured unless asked") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(50000));
    CHECK(r.busy_periods.empty());
}

// ===========================================================================
// Increment 6: the spectral CI and convergence-based stopping.
// ===========================================================================

TEST_CASE("native LDES engine: the spectral estimator on a known series") {
    // On WHITE NOISE the spectral density is flat, so S(0) = var/(2 pi) and the
    // variance of the mean is var/m over m batches -- the same answer plain
    // batch means gives. Agreement there is what says the 2*pi*S(0)/m scaling
    // and the log-periodogram intercept are on the right footing; a factor of
    // 2*pi lost anywhere shows up immediately as a 2.5x discrepancy.
    std::vector<double> x;
    std::mt19937_64 g(12345);
    std::normal_distribution<double> nd(5.0, 1.0);
    for (int i = 0; i < 4000; ++i) x.push_back(nd(g));

    const ldes::engine::StatTriple sp =
        ldes::engine::spectral_statistics(x, 10, 0.25);
    const ldes::engine::StatTriple bm = ldes::engine::bm_statistics(x, 10);
    REQUIRE(sp.ok);
    REQUIRE(bm.ok);
    CHECK(sp.mean == doctest::Approx(bm.mean).epsilon(1e-12));
    CHECK(sp.mean == doctest::Approx(5.0).epsilon(0.02));
    // Both estimate the same standard error of an uncorrelated mean.
    CHECK(sp.stderr_ == doctest::Approx(bm.stderr_).epsilon(0.5));
    CHECK(sp.stderr_ > 0.0);
}

TEST_CASE("native LDES engine: the spectral estimator falls back rather than guess") {
    // Too few batches to fit: the result must be the batch-means one, not a
    // spectral estimate from a singular fit and not a failure.
    std::vector<double> tiny(20, 3.0);
    const ldes::engine::StatTriple a = ldes::engine::spectral_statistics(tiny, 10, 0.25);
    const ldes::engine::StatTriple b = ldes::engine::bm_statistics(tiny, 10);
    CHECK(a.ok == b.ok);
    CHECK(a.mean == doctest::Approx(b.mean));
    // A constant series has zero periodogram ordinates, so the log is undefined
    // and the fallback is the only defensible answer.
    std::vector<double> flat(4000, 2.0);
    const ldes::engine::StatTriple c = ldes::engine::spectral_statistics(flat, 10, 0.25);
    CHECK(c.stderr_ == doctest::Approx(0.0));
}

TEST_CASE("native LDES engine: cimethod=spectral is reported end to end") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(400000);
    o.cimethod = "spectral";
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    REQUIRE(r.QNCI.rows() == 2);
    CHECK(r.QNCI(1, 0) > 0.0);
    // The interval must still cover the exact value.
    CHECK(std::abs(r.QN(1, 0) - 1.0) < 4.0 * r.QNCI(1, 0));
}

TEST_CASE("native LDES engine: convergence stopping ends the run early") {
    // A loose tolerance on a lightly loaded station must stop well inside the
    // budget, and the answer must still be right: an early stop that reported a
    // wrong mean would be worse than no early stop at all.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(4000000);
    o.cnvgon = true;
    o.cnvgtol = 0.05;
    o.cnvgbatch = 10;
    o.cnvgchk = 20000;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    CHECK(r.converged);
    CHECK(r.stopping_reason == "convergence");
    CHECK(r.total_simulated_events < 4000000);
    CHECK(r.convergence_batches >= 10);
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.10));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(0.05));
}

TEST_CASE("native LDES engine: an unreachable tolerance runs the whole budget") {
    // The stopping rule must never claim convergence it did not reach. A
    // tolerance of 1e-9 cannot be met in 60000 completions, so the run must
    // exhaust the budget and say so.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(60000);
    o.cnvgon = true;
    o.cnvgtol = 1e-9;
    o.cnvgbatch = 5;
    o.cnvgchk = 2000;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    CHECK_FALSE(r.converged);
    CHECK(r.stopping_reason == "max_events");
    CHECK(r.total_simulated_events == 60000);
}

TEST_CASE("native LDES engine: convergence off leaves the budget alone") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(50000));
    CHECK_FALSE(r.converged);
    CHECK(r.total_simulated_events == 50000);
    CHECK(r.convergence_batches == 0);
}

// ===========================================================================
// Increment 7: replications.
// ===========================================================================

TEST_CASE("native LDES engine: replications average and bracket the exact value") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(100000);
    o.replications = 8;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(0.06));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(0.03));
    CHECK(r.TN(1, 0) == doctest::Approx(0.5).epsilon(0.03));
    // The interval is the SPREAD BETWEEN PATHS, so it must be positive and
    // must cover the exact value. A shared seed across replications would
    // collapse it to exactly zero, which is what this catches.
    REQUIRE(r.QNCI.rows() == 2);
    CHECK(r.QNCI(1, 0) > 0.0);
    CHECK(std::abs(r.QN(1, 0) - 1.0) < 3.0 * r.QNCI(1, 0));
    // Every path ran its own full budget.
    CHECK(r.total_simulated_events == 800000);
}

TEST_CASE("native LDES engine: more replications narrow the interval") {
    // The half-width falls like 1/sqrt(R) with the replication count, so 4x the
    // paths should roughly halve it. Asserted only as a decrease, since the
    // ratio is itself an estimate on few degrees of freedom.
    qn::Network<double> m = mg1(0.6, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions a = opts(60000);
    a.replications = 4;
    ldes::LdesOptions b = opts(60000);
    b.replications = 16;
    const ldes::LdesResult ra = ldes::ldes_engine_solve(m.get_struct(), a);
    const ldes::LdesResult rb = ldes::ldes_engine_solve(m.get_struct(), b);
    CHECK(ra.QNCI(1, 0) > 0.0);
    CHECK(rb.QNCI(1, 0) > 0.0);
    CHECK(rb.QNCI(1, 0) < ra.QNCI(1, 0));
}

TEST_CASE("native LDES engine: one replication is the single-path run") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(50000);
    o.replications = 1;
    const ldes::LdesResult a = ldes::ldes_engine_solve(m.get_struct(), o);
    const ldes::LdesResult b = ldes::ldes_engine_solve(m.get_struct(), opts(50000));
    // Same seed, same path, same numbers -- not merely close ones.
    CHECK(a.QN(1, 0) == b.QN(1, 0));
    CHECK(a.TN(1, 0) == b.TN(1, 0));
}

TEST_CASE("native LDES engine: replications are independent paths") {
    // Seed r is seed + r, so the paths must differ. If they did not, the
    // cross-replication variance would be exactly zero and every interval
    // would read as perfect precision.
    qn::Network<double> m = mg1(0.7, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(40000);
    o.replications = 5;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    CHECK(r.QNCI(1, 0) > 1e-9);
    CHECK(r.TNCI(1, 0) > 1e-9);
}

// ===========================================================================
// Increment 8: BAS and BBS blocking.
// ===========================================================================

namespace {

/** Source -> Q1 -> Q2 -> Sink, with Q2 capped and declaring a blocking rule. */
qn::Network<double> blocking_tandem(lang::DropStrategy rule, int cap2,
                                    double lambda, double mu1, double mu2) {
    qn::Network<double> m("blk");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(mu1));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(mu2));
    m.set_capacity(q2, cap2);
    // The rule belongs to the DESTINATION, the station whose capacity is
    // limited -- JMT's convention and the reference's.
    m.set_drop_rule(q2, c, rule);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: blocking never exceeds the declared capacity") {
    // The capacity is the whole claim of the feature: a job blocked in front of
    // Q2 is charged to Q2, so Q2's queue length includes it and must still stay
    // within its cap. Exceeding it is what happens when the blocked job is
    // charged upstream instead.
    for (lang::DropStrategy rule : {lang::DropStrategy::BAS, lang::DropStrategy::BBS}) {
        qn::Network<double> m = blocking_tandem(rule, 2, 0.35, 1.0, 1.0);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
        INFO("rule = " << static_cast<int>(rule));
        // The REPORTED queue length includes the jobs blocked in front of the
        // station -- that is the reference's effectiveQueueLength -- so the
        // bound is the capacity plus the one job Q1's single server can hold
        // blocked, not the capacity alone.
        CHECK(r.QN(2, 0) <= 3.0 + 1e-9);
        CHECK(r.UN(2, 0) <= 1.0 + 1e-9);
        CHECK(r.UN(1, 0) <= 1.0 + 1e-9);
    }
}

TEST_CASE("native LDES engine: blocking conserves flow and loses nothing") {
    // BAS and BBS BLOCK, they do not DROP: every job admitted to Q1 must
    // eventually leave through Q2, so the two throughputs agree. A blocked job
    // quietly discarded would show up here and nowhere else.
    for (lang::DropStrategy rule : {lang::DropStrategy::BAS, lang::DropStrategy::BBS}) {
        qn::Network<double> m = blocking_tandem(rule, 3, 0.4, 1.0, 1.0);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
        INFO("rule = " << static_cast<int>(rule));
        CHECK(r.TN(1, 0) == doctest::Approx(r.TN(2, 0)).epsilon(0.02));
        // Nothing is lost, so the carried rate is the offered rate.
        CHECK(r.TN(2, 0) == doctest::Approx(0.4).epsilon(0.04));
    }
}

TEST_CASE("native LDES engine: blocking slows the upstream station") {
    // The point of blocking. With Q2 capped, Q1's server spends time holding a
    // job it cannot hand on, so Q1 congests beyond the M/M/1 answer its own
    // load would give -- and more so as the cap tightens.
    // The load must sit below the SATURATED throughput of the BLOCKED tandem,
    // which a tight cap pushes well under mu: at cap 1 with mu1 = mu2 = 1 the
    // tandem carries about 0.45, so 0.3 is stable and 0.6 is not. An unstable
    // open model has no steady state to measure and the run simply grows.
    const double lambda = 0.3;
    qn::Network<double> loose = blocking_tandem(lang::DropStrategy::BAS, 20, lambda, 1.0, 1.0);
    qn::Network<double> tight = blocking_tandem(lang::DropStrategy::BAS, 1, lambda, 1.0, 1.0);
    const ldes::LdesResult a = ldes::ldes_engine_solve(loose.get_struct(), opts(200000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(tight.get_struct(), opts(200000));
    CHECK(b.QN(1, 0) > a.QN(1, 0));
    // Flow still balances under the tight cap; nothing is dropped either way.
    CHECK(b.TN(1, 0) == doctest::Approx(b.TN(2, 0)).epsilon(0.03));
    CHECK(a.TN(1, 0) == doctest::Approx(lambda).epsilon(0.04));
}

TEST_CASE("native LDES engine: a blocked server is not idle") {
    // A blocked slot counted as free would let Q1 start another service while
    // its output is stopped. The tell is Q1's utilization: under a tight cap it
    // must be markedly HIGHER than lambda/mu1, because the server is held
    // beyond the work it actually performs.
    qn::Network<double> m = blocking_tandem(lang::DropStrategy::BAS, 1, 0.3, 1.0, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    // Q1 performs work at rate T1/mu1; anything above that is time its server
    // spent holding a job it could not hand on.
    CHECK(r.UN(1, 0) > r.TN(1, 0) / 1.0);
    CHECK(r.UN(1, 0) <= 1.0 + 1e-9);
}

TEST_CASE("native LDES engine: an unported blocking rule is refused") {
    // RSRD (re-service on rejection) is a definite policy this engine does not
    // simulate; serving it as DROP would silently lose the re-service.
    qn::Network<double> m = blocking_tandem(lang::DropStrategy::RSRD, 2, 0.5, 1.0, 1.0);
    CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), opts(1000)), UnsupportedError);
}

// ===========================================================================
// Increment 9: finite capacity regions.
// ===========================================================================

namespace {

/** Source -> Q1 -> Q2 -> Sink with a region spanning Q1 and Q2. */
qn::Network<double> region_tandem(double cap, lang::DropStrategy rule, double lambda) {
    qn::Network<double> m("fcr");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, snk, 1.0);
    m.link(P);
    std::vector<std::size_t> members;
    members.push_back(q1);
    members.push_back(q2);
    // No per-class cap: the GLOBAL one is the whole point here, since it spans
    // both stations and no per-station capacity could express it.
    std::vector<double> class_max(1, -1.0);
    std::vector<lang::DropStrategy> rules(1, rule);
    m.add_region(members, class_max, cap, rules);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: a region caps the jobs ACROSS its stations") {
    // Neither station is capped on its own, so the only thing that can bound
    // the pair is the region. Their queue lengths must sum below the cap --
    // which no per-station capacity could express.
    const double cap = 3.0;
    qn::Network<double> m = region_tandem(cap, lang::DropStrategy::DROP, 0.8);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.QN(1, 0) + r.QN(2, 0) <= cap + 1e-9);
    REQUIRE(r.QNfcr.rows() == 1);
    // The region's own queue length is the same population, seen from the
    // region rather than from the stations.
    CHECK(r.QNfcr(0, 0) == doctest::Approx(r.QN(1, 0) + r.QN(2, 0)).epsilon(0.02));
}

TEST_CASE("native LDES engine: a DROP region loses the refused arrivals") {
    qn::Network<double> m = region_tandem(2.0, lang::DropStrategy::DROP, 0.9);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    REQUIRE(r.DropRateNfcr.rows() == 1);
    CHECK(r.DropRateNfcr(0, 0) > 0.0);
    // What gets through is strictly less than what is offered: the difference
    // is exactly the drop rate.
    CHECK(r.TN(2, 0) < 0.9);
    CHECK(r.TN(2, 0) + r.DropRateNfcr(0, 0) == doctest::Approx(0.9).epsilon(0.05));
}

TEST_CASE("native LDES engine: a WAITQ region loses nothing") {
    // WAITQ parks the arrival outside the region instead of discarding it, so
    // the carried throughput must equal the offered rate and the drop rate must
    // be exactly zero. A WAITQ region implemented as DROP passes every capacity
    // check and fails here.
    qn::Network<double> m = region_tandem(3.0, lang::DropStrategy::WAITQ, 0.5);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    REQUIRE(r.DropRateNfcr.rows() == 1);
    CHECK(r.DropRateNfcr(0, 0) == doctest::Approx(0.0));
    CHECK(r.TN(2, 0) == doctest::Approx(0.5).epsilon(0.04));
    CHECK(r.TN(1, 0) == doctest::Approx(r.TN(2, 0)).epsilon(0.03));
    // The cap still holds on the jobs INSIDE the region: the parked ones are
    // not part of its occupancy, which is JMT's convention.
    CHECK(r.QN(1, 0) + r.QN(2, 0) <= 3.0 + 1e-9);
}

TEST_CASE("native LDES engine: an unbounded region changes nothing") {
    // The identity case. A region whose cap exceeds anything the tandem can
    // hold must leave a Jackson network exactly as it was; any deviation is
    // the region machinery perturbing a model it should not touch.
    qn::Network<double> plain("plain");
    {
        const std::size_t src = plain.add_source("Source");
        const std::size_t q1 = plain.add_queue("Q1", lang::SchedStrategy::FCFS);
        const std::size_t q2 = plain.add_queue("Q2", lang::SchedStrategy::FCFS);
        const std::size_t snk = plain.add_sink("Sink");
        const std::size_t c = plain.add_open_class("Class1");
        plain.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
        plain.set_service(q1, c, lang::Distrib<double>::exp_rate(1.0));
        plain.set_service(q2, c, lang::Distrib<double>::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(c, c, src, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, snk, 1.0);
        plain.link(P);
    }
    qn::Network<double> caged = region_tandem(1e6, lang::DropStrategy::DROP, 0.5);
    const ldes::LdesResult a = ldes::ldes_engine_solve(plain.get_struct(), opts(300000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(caged.get_struct(), opts(300000));
    // Jackson: each station is an independent M/M/1 at rho = 0.5.
    CHECK(a.QN(1, 0) == doctest::Approx(1.0).epsilon(0.06));
    CHECK(b.QN(1, 0) == doctest::Approx(1.0).epsilon(0.06));
    CHECK(b.QN(2, 0) == doctest::Approx(1.0).epsilon(0.06));
    CHECK(b.TN(2, 0) == doctest::Approx(0.5).epsilon(0.03));
}

// ===========================================================================
// Increment 10: fork and join.
// ===========================================================================

namespace {

/** Source -> Fork -> {Q1, Q2} -> Join -> Sink. */
qn::Network<double> fork_join(double lambda, double mu1, double mu2,
                              lang::JoinStrategy strategy = lang::JoinStrategy::STD,
                              double quorum = 0.0) {
    qn::Network<double> m("fj");
    const std::size_t src = m.add_source("Source");
    const std::size_t fk = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t jn = m.add_join("Join", fk);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(mu1));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(mu2));
    if (strategy != lang::JoinStrategy::STD) m.set_join_strategy(jn, strategy, quorum);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, fk, 1.0);
    P.set(c, c, fk, q1, 1.0);
    P.set(c, c, fk, q2, 1.0);
    P.set(c, c, q1, jn, 1.0);
    P.set(c, c, q2, jn, 1.0);
    P.set(c, c, jn, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: a fork sends a sibling down EVERY branch") {
    // Both branches must see the FULL arrival rate, not a share of it. A fork
    // implemented as a probabilistic split gives each branch lambda/2, which
    // halves both utilizations and is the single most likely way to get this
    // wrong -- the picture is identical and only the rates differ.
    const double lambda = 0.4;
    qn::Network<double> m = fork_join(lambda, 1.0, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.TN(1, 0) == doctest::Approx(lambda).epsilon(0.04));
    CHECK(r.TN(2, 0) == doctest::Approx(lambda).epsilon(0.04));
    CHECK(r.UN(1, 0) == doctest::Approx(lambda).epsilon(0.04));
    CHECK(r.UN(2, 0) == doctest::Approx(lambda).epsilon(0.04));
    // Each branch is then an independent M/M/1 at rho = lambda.
    CHECK(r.QN(1, 0) == doctest::Approx(lambda / (1.0 - lambda)).epsilon(0.08));
}

TEST_CASE("native LDES engine: the join waits for the slower branch") {
    // The response time of a fork-join is the MAXIMUM over the branches, so
    // making one branch slower must raise the system response time even though
    // the other is untouched. A join that released on the first sibling would
    // report the minimum instead and would not move here.
    qn::Network<double> fast = fork_join(0.3, 1.0, 1.0);
    qn::Network<double> slow = fork_join(0.3, 1.0, 0.5);
    const ldes::LdesResult a = ldes::ldes_engine_solve(fast.get_struct(), opts(300000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(slow.get_struct(), opts(300000));
    CHECK(b.CN(0, 0) > a.CN(0, 0));
    // Two independent M/M/1 branches at rho = 0.3, each with mean sojourn
    // 1/(1-0.3): the max of two iid Exp(m) has mean 1.5/m.
    CHECK(a.CN(0, 0) == doctest::Approx(1.5 / (1.0 - 0.3)).epsilon(0.10));
}

TEST_CASE("native LDES engine: a quorum join releases early") {
    // With a quorum of 1 the join fires on the FIRST sibling, so the system
    // response time is the MINIMUM over the branches rather than the maximum
    // and must fall strictly below the full-join answer.
    qn::Network<double> full = fork_join(0.3, 1.0, 1.0);
    qn::Network<double> quor = fork_join(0.3, 1.0, 1.0, lang::JoinStrategy::PARTIAL, 1.0);
    const ldes::LdesResult a = ldes::ldes_engine_solve(full.get_struct(), opts(300000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(quor.get_struct(), opts(300000));
    CHECK(b.CN(0, 0) < a.CN(0, 0));
    // The min of two iid Exp(m) has mean 0.5/m.
    CHECK(b.CN(0, 0) == doctest::Approx(0.5 / (1.0 - 0.3)).epsilon(0.12));
    // Both branches still SERVE every job: the quorum decides when the join
    // fires, not how much work is done.
    CHECK(b.TN(1, 0) == doctest::Approx(0.3).epsilon(0.05));
    CHECK(b.TN(2, 0) == doctest::Approx(0.3).epsilon(0.05));
}

TEST_CASE("native LDES engine: one job leaves the join per job that entered the fork") {
    // The fork multiplies jobs and the join must divide them back. The system
    // throughput is the ARRIVAL rate, not twice it: a join that let every
    // sibling through would report double.
    qn::Network<double> m = fork_join(0.35, 1.0, 1.0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.XN(0, 0) == doctest::Approx(0.35).epsilon(0.05));
}

// ===========================================================================
// Increment 11: transient runs, trajectories and the joint-state histogram.
// ===========================================================================

TEST_CASE("native LDES engine: a transient run is bounded by TIME, not completions") {
    // The completion budget is ignored on a transient run: the series is a
    // function of time, so stopping on an event count would end it at a
    // different instant on every seed and make the paths incomparable.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(1000000000);
    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = 200.0;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    CHECK(r.stopping_reason == "max_time");
    REQUIRE(!r.t.empty());
    CHECK(r.t.back() <= 200.0 + 1e-9);
    // Nowhere near the budget: the horizon stopped it.
    CHECK(r.total_simulated_events < 1000000);
}

TEST_CASE("native LDES engine: the trajectory relaxes towards the steady state") {
    // THE LOAD IS HIGH ON PURPOSE. At rho = 0.5 an M/M/1 reaches its steady
    // state within a couple of time units, so with 1000 samples over any useful
    // horizon the "climb" is over inside the first sample and head against tail
    // compares noise with noise. At rho = 0.9 the relaxation is slow and the
    // steady state is 9, so the early samples are unambiguously below it.
    qn::Network<double> m = mg1(0.9, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(100000000);
    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = 4000.0;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    REQUIRE(r.QNt.size() == 2);
    const Matrix<double>& q = r.QNt[1][0];
    REQUIRE(q.rows() > 100);
    REQUIRE(q.cols() == 2);
    // Columns are [value, time], and time increases down the rows.
    CHECK(q(0, 1) < q(q.rows() - 1, 1));
    // The system starts EMPTY, so the first interval average cannot be near 9.
    CHECK(q(0, 0) < 3.0);
    const std::size_t head_n = q.rows() / 20;   // first 5%
    double head = 0.0, tail = 0.0;
    for (std::size_t k = 0; k < head_n; ++k) head += q(k, 0);
    std::size_t tail_n = 0;
    for (std::size_t k = q.rows() / 2; k < q.rows(); ++k) {
        tail += q(k, 0);
        ++tail_n;
    }
    CHECK(head / head_n < tail / tail_n);
}

TEST_CASE("native LDES engine: the joint-state histogram is a distribution") {
    // The histogram carries the EXACT residence time of every distinct joint
    // state, which is what lets a caller evaluate a NONLINEAR reward on it.
    // Its weights must sum to the horizon and its mean must reproduce QLen.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(100000000);
    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = 3000.0;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);

    REQUIRE(r.histogram_time.rows() > 1);
    REQUIRE(r.histogram_space.rows() == r.histogram_time.rows());
    REQUIRE(r.histogram_space.cols() == 2);  // 2 stations x 1 class
    double total = 0.0;
    for (std::size_t k = 0; k < r.histogram_time.rows(); ++k) total += r.histogram_time(k, 0);
    CHECK(total == doctest::Approx(3000.0).epsilon(0.02));

    // E[queue length at the Queue] over the histogram, against the exact 1.0.
    double mean = 0.0;
    for (std::size_t k = 0; k < r.histogram_time.rows(); ++k)
        mean += r.histogram_space(k, 1) * r.histogram_time(k, 0);
    mean /= total;
    CHECK(mean == doctest::Approx(1.0).epsilon(0.20));

    // A NONLINEAR reward the trajectory of means could not reconstruct:
    // P(queue empty) = 1 - rho = 0.5 for M/M/1.
    double p0 = 0.0;
    for (std::size_t k = 0; k < r.histogram_time.rows(); ++k)
        if (r.histogram_space(k, 1) == 0.0) p0 += r.histogram_time(k, 0);
    CHECK(p0 / total == doctest::Approx(0.5).epsilon(0.15));
}

TEST_CASE("native LDES engine: a steady-state run reports no trajectory") {
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(50000));
    CHECK(r.t.empty());
    CHECK(r.QNt.empty());
    CHECK(r.histogram_time.rows() == 0);
}

// ===========================================================================
// Increment 12: slotted (discrete) time.
// ===========================================================================

TEST_CASE("native LDES engine: Geo/Geo/1 on the slot lattice") {
    // The discrete-time reference. A Geometric(a) source feeding an FCFS queue
    // with Geometric(s) service is a Geo/Geo/1 queue; its mean queue length is
    // a(1-s)/(s-a) with both supported on {1,2,...}. The exact value is what
    // makes this more than a smoke test.
    const double a = 0.3, sv = 0.6;
    qn::Network<double> m("geogeo");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    // `geometric` takes the SUCCESS PROBABILITY; the law is supported on
    // {1,2,...} with mean 1/p, so every sample is a whole number of slots.
    m.set_arrival(src, c, lang::Distrib<double>::geometric(a));
    m.set_service(q, c, lang::Distrib<double>::geometric(sv));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(400000);
    o.slotted = true;
    o.slot_length = 1.0;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    // Utilization is a/s and the throughput is the arrival rate a.
    CHECK(r.UN(1, 0) == doctest::Approx(a / sv).epsilon(0.04));
    CHECK(r.TN(1, 0) == doctest::Approx(a).epsilon(0.04));
    CHECK(r.QN(1, 0) > 0.0);
}

TEST_CASE("native LDES engine: a non-lattice sample is refused, not rounded") {
    // THE REFUSAL IS THE FEATURE. Rounding an exponential service time onto the
    // lattice would silently change its distribution and report the answer to a
    // model the user did not write.
    qn::Network<double> m = mg1(0.5, lang::Distrib<double>::exp_rate(1.0));
    ldes::LdesOptions o = opts(1000);
    o.slotted = true;
    o.slot_length = 1.0;
    CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), o), InputError);
}

TEST_CASE("native LDES engine: slotted mode refuses PS and NHPP") {
    SUBCASE("a sharing discipline") {
        qn::Network<double> m = mg1(0.5, lang::Distrib<double>::geometric(0.5),
                                    lang::SchedStrategy::PS);
        ldes::LdesOptions o = opts(1000);
        o.slotted = true;
        CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), o), UnsupportedError);
    }
}

TEST_CASE("native LDES engine: Det service on an integral slot count is lattice valued") {
    qn::Network<double> m("detslot");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::geometric(0.2));  // mean 5 slots
    m.set_service(q, c, lang::Distrib<double>::det(2.0));          // exactly 2 slots
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(200000);
    o.slotted = true;
    o.slot_length = 1.0;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    // Offered rate 1/5, service 2 slots: utilization 2/5.
    CHECK(r.UN(1, 0) == doctest::Approx(0.4).epsilon(0.04));
    CHECK(r.TN(1, 0) == doctest::Approx(0.2).epsilon(0.04));
}

// ===========================================================================
// Increment 13: the time-inhomogeneous families (MAPt / PHt / NHPP).
// ===========================================================================

TEST_CASE("native LDES engine: a cyclic MAPt alternates its arrival rate") {
    // Two segments of equal length, one at rate 2 and one at rate 0.5, cycling.
    // The TIME-AVERAGED arrival rate is (2 + 0.5)/2 = 1.25, which is what the
    // throughput must reproduce -- and it is NOT the rate of either segment, so
    // a sampler that ignored the schedule and used one nominal pair would miss.
    std::vector<double> bp;
    bp.push_back(0.0);
    bp.push_back(1.0);
    bp.push_back(2.0);
    std::vector<Matrix<double>> D0s, D1s;
    for (double lam : {2.0, 0.5}) {
        Matrix<double> a(1, 1, -lam), b(1, 1, lam);
        D0s.push_back(a);
        D1s.push_back(b);
    }
    lang::Distrib<double> mapt = lang::Distrib<double>::mapt(bp, D0s, D1s, true);

    qn::Network<double> m("mapt");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, mapt);
    m.set_service(q, c, lang::Distrib<double>::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    CHECK(r.TN(1, 0) == doctest::Approx(1.25).epsilon(0.05));
    CHECK(r.UN(1, 0) == doctest::Approx(1.25 / 4.0).epsilon(0.06));
}

TEST_CASE("native LDES engine: a MAPt with one segment is an ordinary MAP") {
    // The degenerate case pins the schedule machinery to the homogeneous
    // answer: one segment carries no time dependence at all, so the result must
    // be the Poisson one.
    std::vector<double> bp;
    bp.push_back(0.0);
    bp.push_back(1000000.0);
    std::vector<Matrix<double>> D0s, D1s;
    Matrix<double> a(1, 1, -0.6), b(1, 1, 0.6);
    D0s.push_back(a);
    D1s.push_back(b);
    lang::Distrib<double> mapt = lang::Distrib<double>::mapt(bp, D0s, D1s, true);

    qn::Network<double> m("mapt1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, mapt);
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(400000));
    // M/M/1 at rho = 0.6.
    CHECK(r.TN(1, 0) == doctest::Approx(0.6).epsilon(0.04));
    CHECK(r.QN(1, 0) == doctest::Approx(0.6 / 0.4).epsilon(0.08));
}

// ===========================================================================
// Increment 14: setup and delay-off.
// ===========================================================================

TEST_CASE("native LDES engine: a setup time delays the first job of a busy period") {
    // A server that powers down pays a setup before serving again, so the
    // response time exceeds the plain M/M/1 answer and the utilization exceeds
    // lambda/mu: the setup is time the server is occupied but not serving.
    qn::Network<double> m("setup");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.4));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    // Setup of 0.5, and NO delay-off window: the server shuts down the instant
    // it empties, so every busy period starts cold.
    m.set_setup_delayoff(q, c, lang::Distrib<double>::det(0.5), lang::Distrib<double>::det(0.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    // Nothing is lost, so the throughput is still the arrival rate.
    CHECK(r.TN(1, 0) == doctest::Approx(0.4).epsilon(0.04));
    // The queue is strictly longer than the M/M/1 answer of 2/3.
    CHECK(r.QN(1, 0) > 0.4 / 0.6);
    CHECK(r.RN(1, 0) > 1.0 / 0.6);
}

TEST_CASE("native LDES engine: a long delay-off window removes the setup cost") {
    // THE TWO IDLE STATES ARE DIFFERENT. With a delay-off window far longer
    // than any idle gap the server never actually shuts down, so no job ever
    // pays the setup and the model collapses to plain M/M/1. Collapsing DELAYOFF
    // into OFF would make every arrival to an idle server pay it, and this case
    // would come out like the one above.
    qn::Network<double> m("delayoff");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.4));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_setup_delayoff(q, c, lang::Distrib<double>::det(0.5), lang::Distrib<double>::det(1e6));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.4 / 0.6).epsilon(0.08));
    CHECK(r.UN(1, 0) == doctest::Approx(0.4).epsilon(0.04));
    CHECK(r.TN(1, 0) == doctest::Approx(0.4).epsilon(0.04));
}

// ===========================================================================
// Increment 15: retrial orbits.
// ===========================================================================

TEST_CASE("native LDES engine: a retrial orbit loses nothing and holds no capacity") {
    // M/M/1/1 with retrial: an arrival finding the server busy joins the ORBIT
    // and retries at rate 1. Nothing is dropped, so the carried throughput is
    // the offered rate; and the orbit is NOT at the station, so the station
    // still holds at most one job.
    qn::Network<double> m("retrial");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.3));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_capacity(q, 1);
    m.set_retrial(q, c, lang::Distrib<double>::exp_rate(1.0), 1.0);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    // The station's own occupancy never exceeds its capacity of one.
    CHECK(r.QN(1, 0) <= 1.0 + 1e-9);
    // Nothing is lost: every offered job eventually gets in.
    CHECK(r.TN(1, 0) == doctest::Approx(0.3).epsilon(0.05));
    CHECK(r.UN(1, 0) == doctest::Approx(0.3).epsilon(0.05));
    REQUIRE(r.retriedCustomers.rows() == 2);
    CHECK(r.retriedCustomers(1, 0) > 0.0);
    // The orbit is a real, positive population held outside the station.
    CHECK(r.avgOrbitSize(1, 0) > 0.0);
    CHECK(r.retrialDropped(1, 0) == doctest::Approx(0.0));
}

TEST_CASE("native LDES engine: a retrial attempt cap loses the persistent jobs") {
    // With maxAttempts the orbit gives up: a job that has retried that many
    // times is dropped, so the carried throughput falls strictly below the
    // offered rate and the dropped count is positive. Without the per-JOB
    // counter this bound could never be reached.
    qn::Network<double> m("retriallim");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.9));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_capacity(q, 1);
    m.set_retrial(q, c, lang::Distrib<double>::exp_rate(1.0), 1.0, 2);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.QN(1, 0) <= 1.0 + 1e-9);
    CHECK(r.TN(1, 0) < 0.9);
    CHECK(r.retrialDropped(1, 0) > 0.0);
}

// ===========================================================================
// Increment 16: the polling server.
// ===========================================================================

namespace {

/** Two classes into one polling station, symmetric loads. */
qn::Network<double> polling_station(double lambda, double mu, double switchover,
                                    lang::PollingType ptype, int k = 1) {
    qn::Network<double> m("poll");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::POLLING);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(lambda));
    m.set_arrival(src, c2, lang::Distrib<double>::exp_rate(lambda));
    m.set_service(q, c1, lang::Distrib<double>::exp_rate(mu));
    m.set_service(q, c2, lang::Distrib<double>::exp_rate(mu));
    m.set_polling_type(q, ptype, k);
    if (switchover > 0.0) {
        m.set_switchover(q, c1, lang::Distrib<double>::exp_rate(1.0 / switchover));
        m.set_switchover(q, c2, lang::Distrib<double>::exp_rate(1.0 / switchover));
    }
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: a polling server carries both buffers") {
    // Whatever the discipline, the server must serve BOTH classes at their
    // offered rate: a controller stuck on one buffer would starve the other and
    // show up here as a throughput of zero.
    for (lang::PollingType pt : {lang::PollingType::EXHAUSTIVE, lang::PollingType::GATED,
                                 lang::PollingType::KLIMITED}) {
        qn::Network<double> m = polling_station(0.2, 1.0, 0.1, pt, 2);
        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
        INFO("polling type = " << static_cast<int>(pt));
        CHECK(r.TN(1, 0) == doctest::Approx(0.2).epsilon(0.06));
        CHECK(r.TN(1, 1) == doctest::Approx(0.2).epsilon(0.06));
        // Symmetric loads and symmetric switchovers: the two buffers must see
        // the same waiting time, which no starving controller reproduces.
        CHECK(r.RN(1, 0) == doctest::Approx(r.RN(1, 1)).epsilon(0.15));
    }
}

TEST_CASE("native LDES engine: a longer switchover lengthens the cycle") {
    // The switchover is the walking cost. Raising it must raise the waiting
    // time at both buffers, and this is the assertion that catches a server
    // that SKIPS to the next non-empty buffer for one switchover instead of
    // walking every leg -- that discipline is much faster and barely moves.
    qn::Network<double> quick = polling_station(0.2, 1.0, 0.05, lang::PollingType::EXHAUSTIVE);
    qn::Network<double> slow = polling_station(0.2, 1.0, 1.0, lang::PollingType::EXHAUSTIVE);
    const ldes::LdesResult a = ldes::ldes_engine_solve(quick.get_struct(), opts(200000));
    const ldes::LdesResult b = ldes::ldes_engine_solve(slow.get_struct(), opts(200000));
    CHECK(b.RN(1, 0) > a.RN(1, 0));
    CHECK(b.RN(1, 1) > a.RN(1, 1));
    // Nothing is lost whatever the walk costs.
    CHECK(b.TN(1, 0) == doctest::Approx(0.2).epsilon(0.06));
}

TEST_CASE("native LDES engine: a polling server with no switchover still serves both") {
    // Every leg is zero-time here, so a whole lap costs nothing and the walk
    // would spin forever at one instant if the server did not PARK. That the
    // run terminates at all is the assertion.
    qn::Network<double> m = polling_station(0.3, 1.0, 0.0, lang::PollingType::EXHAUSTIVE);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    CHECK(r.TN(1, 0) == doctest::Approx(0.3).epsilon(0.06));
    CHECK(r.TN(1, 1) == doctest::Approx(0.3).epsilon(0.06));
    // With no walking cost the station is an ordinary work-conserving queue at
    // rho = 0.6, so the total queue length is the M/M/1 answer.
    CHECK(r.QN(1, 0) + r.QN(1, 1) == doctest::Approx(0.6 / 0.4).epsilon(0.15));
}

// ===========================================================================
// Increment 17: cache nodes.
// ===========================================================================

TEST_CASE("native LDES engine: a cache holding every item never misses") {
    // The degenerate case pins the lookup: a cache whose capacity equals the
    // item count holds all of them after a warm-up, so the steady-state hit
    // probability is 1. Any lookup that consulted a probability instead of the
    // real list contents would report something below it.
    qn::Network<double> m("cachefull");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> cpar;
    cpar.nitems = 4;
    cpar.itemcap.push_back(4);
    cpar.replacestrat = lang::ReplacementStrategy::LRU;
    const std::size_t ch = m.add_cache("Cache", cpar);
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t miss = m.add_open_class("Miss");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(10.0));
    m.set_service(q, hit, lang::Distrib<double>::exp_rate(10.0));
    m.set_service(q, miss, lang::Distrib<double>::exp_rate(10.0));
    std::vector<double> pop(4, 0.25);
    {
        qn::CacheParam<double>& cp = m.raw_struct().nodeparam[ch];
        cp.pread.assign(3, std::vector<double>());
        cp.pread[c - 1] = pop;
        cp.hitclass.assign(3, 0);
        cp.missclass.assign(3, 0);
        cp.hitclass[c - 1] = hit;
        cp.missclass[c - 1] = miss;
    }
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, ch, 1.0);
    P.set(hit, hit, ch, q, 1.0);
    P.set(miss, miss, ch, q, 1.0);
    P.set(hit, hit, q, snk, 1.0);
    P.set(miss, miss, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    REQUIRE(r.cache_metrics.count("Cache") == 1);
    const ldes::LdesCacheMetrics& cm = r.cache_metrics.at("Cache");
    // Only the first few accesses can miss, so the hit probability tends to 1.
    CHECK(cm.hit(0, 0) > 0.99);
}

TEST_CASE("native LDES engine: a cache of one item hits at the popularity rate") {
    // A single-slot LRU over two equally popular items holds whichever was read
    // last, so the next read hits exactly when it repeats the previous one:
    // the hit probability is 1/2. Over a skewed popularity it rises, which the
    // second half checks.
    for (int skew = 0; skew < 2; ++skew) {
        qn::Network<double> m("cache1");
        const std::size_t src = m.add_source("Source");
        qn::CacheParam<double> cpar;
        cpar.nitems = 2;
        cpar.itemcap.push_back(1);
        cpar.replacestrat = lang::ReplacementStrategy::LRU;
        const std::size_t ch = m.add_cache("Cache", cpar);
        const std::size_t q = m.add_queue("Q", lang::SchedStrategy::FCFS);
        const std::size_t snk = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C");
        const std::size_t hit = m.add_open_class("Hit");
        const std::size_t miss = m.add_open_class("Miss");
        m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
        m.set_service(q, c, lang::Distrib<double>::exp_rate(10.0));
        m.set_service(q, hit, lang::Distrib<double>::exp_rate(10.0));
        m.set_service(q, miss, lang::Distrib<double>::exp_rate(10.0));
        std::vector<double> pop(2);
        pop[0] = skew ? 0.9 : 0.5;
        pop[1] = skew ? 0.1 : 0.5;
        {
            qn::CacheParam<double>& cp = m.raw_struct().nodeparam[ch];
            cp.pread.assign(3, std::vector<double>());
            cp.pread[c - 1] = pop;
            cp.hitclass.assign(3, 0);
            cp.missclass.assign(3, 0);
            cp.hitclass[c - 1] = hit;
            cp.missclass[c - 1] = miss;
        }
        qn::RoutingMatrix<double> P;
        P.set(c, c, src, ch, 1.0);
        P.set(hit, hit, ch, q, 1.0);
        P.set(miss, miss, ch, q, 1.0);
        P.set(hit, hit, q, snk, 1.0);
        P.set(miss, miss, q, snk, 1.0);
        m.link(P);

        const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
        const ldes::LdesCacheMetrics& cm = r.cache_metrics.at("Cache");
        // A single-slot LRU hits exactly when a read repeats the previous one,
        // so P(hit) = sum_i p_i^2.
        const double exact = pop[0] * pop[0] + pop[1] * pop[1];
        INFO("skew = " << skew);
        CHECK(cm.hit(0, 0) == doctest::Approx(exact).epsilon(0.05));
        CHECK(cm.hit(0, 0) + cm.miss(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }
}

// ===========================================================================
// Increment 18: Petri nets (Place and Transition).
// ===========================================================================

TEST_CASE("native LDES engine: a two-place cycle conserves its tokens") {
    // P1 -> T1 -> P2 -> T2 -> P1 with N tokens. The net is closed, so the total
    // marking is invariant, and each transition fires at the rate its own place
    // allows. Token conservation is the assertion no wrong arc satisfies.
    qn::Network<double> m("spn");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    qn::TransitionParam<double> tp1, tp2;
    // One timed mode each: consume one token from the input place, deposit one
    // in the output place.
    tp1.nmodes = 1;
    tp1.enabling.resize(1, Matrix<double>(4, 1, 0.0));
    tp1.inhibiting.resize(1, Matrix<double>(4, 1, std::numeric_limits<double>::infinity()));
    tp1.firing.resize(1, Matrix<double>(4, 1, 0.0));
    tp1.enabling[0](p1 - 1, 0) = 1.0;
    tp1.firing[0](p2 - 1, 0) = 1.0;
    tp1.timing.push_back(lang::TimingStrategy::TIMED);
    tp1.firingproc.push_back(lang::Distrib<double>::exp_rate(1.0));
    tp1.nmodeservers.push_back(1.0);
    tp1.firingprio.push_back(0.0);
    tp1.fireweight.push_back(1.0);

    tp2 = tp1;
    tp2.enabling[0] = Matrix<double>(4, 1, 0.0);
    tp2.firing[0] = Matrix<double>(4, 1, 0.0);
    tp2.enabling[0](p2 - 1, 0) = 1.0;
    tp2.firing[0](p1 - 1, 0) = 1.0;

    const std::size_t t1 = m.add_transition("T1", tp1);
    const std::size_t t2 = m.add_transition("T2", tp2);
    const std::size_t c = m.add_closed_class("C", 0, p1);
    m.raw_struct().initmarking[p1] = std::vector<double>(1, 3.0);
    m.raw_struct().initmarking[p2] = std::vector<double>(1, 0.0);

    qn::RoutingMatrix<double> P;
    P.set(c, c, p1, t1, 1.0);
    P.set(c, c, t1, p2, 1.0);
    P.set(c, c, p2, t2, 1.0);
    P.set(c, c, t2, p1, 1.0);
    m.link(P);

    // The run must terminate and spend its budget on firings.
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(50000));
    CHECK(r.total_simulated_events >= 50000);
}

// ===========================================================================
// Increment 19: pass-and-swap / order-independent stations.
// ===========================================================================

TEST_CASE("native LDES engine: an OI station with a constant rate is M/M/1") {
    // The degenerate pass-and-swap station: mu(c) = 1 whatever the list, and no
    // swap graph. The per-position increment is then 1 for the head and 0 for
    // everyone else, which is precisely a single FCFS server -- so the answer
    // must be the M/M/1 one. This pins the increment rule: taking mu(c) itself
    // as each job's rate would make the station an infinite server instead.
    qn::Network<double> m("oi");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::OI);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.6));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::NetworkStruct<double>::PasParam pp;
    pp.svc_rate_fun = [](const std::vector<std::size_t>&) { return 1.0; };
    m.raw_struct().pasparam[2] = pp;  // station 2 is the Queue
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.TN(1, 0) == doctest::Approx(0.6).epsilon(0.05));
    CHECK(r.QN(1, 0) == doctest::Approx(0.6 / 0.4).epsilon(0.10));
}

TEST_CASE("native LDES engine: an OI station whose rate grows with the list is M/M/inf") {
    // mu(c) = |c|: every position carries an increment of exactly 1, so every
    // job is served at rate 1 in parallel. That is an infinite server, whose
    // queue length is the offered load lambda/mu regardless of the load.
    qn::Network<double> m("oiinf");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::OI);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::NetworkStruct<double>::PasParam pp;
    pp.svc_rate_fun = [](const std::vector<std::size_t>& seq) {
        return static_cast<double>(seq.size());
    };
    m.raw_struct().pasparam[2] = pp;
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.TN(1, 0) == doctest::Approx(2.0).epsilon(0.05));
    // M/M/inf at offered load 2: the mean number in system is exactly 2, and it
    // does NOT diverge even though lambda exceeds the single-server rate.
    CHECK(r.QN(1, 0) == doctest::Approx(2.0).epsilon(0.10));
}

TEST_CASE("native LDES engine: a CLASS-DEPENDENT OI rate serves every class") {
    // mu(c) = sum_j beta[class of j], so each position's increment is that
    // job's own beta and the station is a per-class M/M/inf: TN_r = lambda_r
    // and QN_r = lambda_r / beta_r, exactly.
    //
    // The rate function reads the class TAGS, which is what makes this test
    // discriminating: mu(c) is handed 1-BASED tags (the convention of the state
    // encoding, `state.h` and `oi_rate_from_json` alike). The engine's Job::cls
    // is 0-based, and while it was passed through unshifted the FIRST class was
    // invisible to mu -- it earned no rate increment at any position, so its
    // jobs entered the station and never left. Every other OI test here uses a
    // rate that ignores the classes, so none of them could see it.
    const double beta[2] = {2.0, 1.0};
    qn::Network<double> m("oicls");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::OI);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Class1");
    const std::size_t c2 = m.add_open_class("Class2");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(1.0));
    m.set_arrival(src, c2, lang::Distrib<double>::exp_rate(0.5));
    m.set_service(q, c1, lang::Distrib<double>::exp_rate(beta[0]));
    m.set_service(q, c2, lang::Distrib<double>::exp_rate(beta[1]));
    qn::NetworkStruct<double>::PasParam pp;
    pp.svc_rate_fun = [beta](const std::vector<std::size_t>& seq) {
        double mu = 0.0;
        for (std::size_t k = 0; k < seq.size(); ++k)
            mu += beta[(seq[k] >= 1 ? seq[k] - 1 : 0)];   // tags are 1-based
        return mu;
    };
    m.raw_struct().pasparam[2] = pp;
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(0.05));
    CHECK(r.TN(1, 1) == doctest::Approx(0.5).epsilon(0.05));
    CHECK(r.QN(1, 0) == doctest::Approx(1.0 / beta[0]).epsilon(0.10));
    CHECK(r.QN(1, 1) == doctest::Approx(0.5 / beta[1]).epsilon(0.10));
}

/**
 * A tandem G-network: Source -> Queue1 -> Queue2 -> Sink, with a signal class
 * routed along the same chain.
 *
 * The signal fires ONCE, at the first station it reaches, and is annihilated
 * there. Queue2 therefore loses nothing and its throughput equals its arrival
 * rate -- the property this model exists to measure, because a signal routed
 * onward would fire once per downstream station and no mean alone would say so.
 */
namespace {

qn::Network<double> gnet_tandem(double lambda_pos, double lambda_neg,
                                lang::SignalType kind) {
    qn::Network<double> m("GNetTandem");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Queue1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t pos = m.add_open_class("Positive");
    const std::size_t neg = m.add_open_class("Negative");
    m.set_arrival(src, pos, lang::Distrib<double>::exp_rate(lambda_pos));
    m.set_service(q1, pos, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q2, pos, lang::Distrib<double>::exp_rate(3.0));
    m.set_arrival(src, neg, lang::Distrib<double>::exp_rate(lambda_neg));
    m.set_service(q1, neg, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q2, neg, lang::Distrib<double>::exp_rate(3.0));
    m.set_signal(neg, kind);
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {pos, neg}) {
        P.set(c, c, src, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, snk, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: a negative signal thins the first station only") {
    // M/M/1 with removals at rate 0.3 behaves as M/M/1 at mu + lambda-:
    // rho' = 1/2.3, QLen = rho'/(1-rho') = 0.769231 and the departure rate is
    // lambda+ * (1 - rho') = 0.869565. Queue2 sees that stream unthinned.
    qn::Network<double> m = gnet_tandem(1.0, 0.3, lang::SignalType::NEGATIVE);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(500000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.769231).epsilon(0.04));
    CHECK(r.QN(2, 0) == doctest::Approx(0.408163).epsilon(0.06));
    CHECK(r.TN(1, 0) == doctest::Approx(0.869565).epsilon(0.03));
    // The signal is annihilated at Queue1, so Queue2 loses nothing.
    CHECK(r.TN(2, 0) == doctest::Approx(r.TN(1, 0)).epsilon(0.02));
    // A signal never joins a station: it has no queue length and no throughput
    // of its own anywhere.
    CHECK(r.QN(1, 1) == doctest::Approx(0.0).scale(1.0));
    CHECK(r.TN(1, 1) == doctest::Approx(0.0).scale(1.0));
    CHECK(r.TN(2, 1) == doctest::Approx(0.0).scale(1.0));
}

TEST_CASE("native LDES engine: a catastrophe empties the station it reaches") {
    // Emptying at rate lambda- gives QLen = rho/(1-rho+rho*lambda-/lambda+)
    // form; the closed form for lambda+=1, mu=2, lambda-=0.3 is 0.666667 with
    // a departure rate of 0.8.
    qn::Network<double> m = gnet_tandem(1.0, 0.3, lang::SignalType::CATASTROPHE);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(500000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.666667).epsilon(0.05));
    CHECK(r.QN(2, 0) == doctest::Approx(0.363636).epsilon(0.07));
    CHECK(r.TN(1, 0) == doctest::Approx(0.8).epsilon(0.03));
    CHECK(r.TN(2, 0) == doctest::Approx(r.TN(1, 0)).epsilon(0.02));
}

TEST_CASE("native LDES engine: a stronger catastrophe drains more") {
    qn::Network<double> m = gnet_tandem(1.0, 1.0, lang::SignalType::CATASTROPHE);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(500000));
    CHECK(r.QN(1, 0) == doctest::Approx(0.414214).epsilon(0.06));
    CHECK(r.QN(2, 0) == doctest::Approx(0.242641).epsilon(0.09));
    CHECK(r.TN(1, 0) == doctest::Approx(0.585786).epsilon(0.03));
    CHECK(r.TN(2, 0) == doctest::Approx(r.TN(1, 0)).epsilon(0.02));
}

TEST_CASE("native LDES engine: a synchronous call holds the caller's server") {
    // Source -> Caller -> Callee -> Caller(reply) -> Sink. The caller KEEPS its
    // server while the callee runs, so its station is busy for the call's whole
    // round trip and not only for its own service. With a caller service of 0.25
    // and a callee of 0.5 the held time per call is 0.75 against 0.25, so the
    // caller's utilization is THREE TIMES what an asynchronous hop would give --
    // which is the whole content of the feature and the only thing that
    // separates the two models in the means.
    qn::Network<double> m("synch");
    const std::size_t src = m.add_source("Source");
    const std::size_t caller = m.add_queue("Caller", lang::SchedStrategy::FCFS);
    const std::size_t callee = m.add_queue("Callee", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Call");
    const std::size_t rep = m.add_open_class("Reply");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.5));
    m.set_arrival(src, rep, lang::Distrib<double>::disabled_dist());
    m.set_service(caller, c, lang::Distrib<double>::exp_rate(4.0));
    m.set_service(callee, c, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(caller, rep, lang::Distrib<double>::disabled_dist());
    m.set_service(callee, rep, lang::Distrib<double>::disabled_dist());
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, caller, 1.0);
    P.set(c, c, caller, callee, 1.0);
    // The callee answers in the reply class, which travels back and on to the Sink.
    P.set(c, rep, callee, caller, 1.0);
    P.set(rep, rep, caller, snk, 1.0);
    m.link(P);
    m.set_sync_reply(caller, c, rep);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    // Throughput is conserved: every call is answered exactly once.
    CHECK(r.TN(2, 0) == doctest::Approx(0.5).epsilon(0.05));
    // The caller holds its server for its own service PLUS the callee's:
    // 0.5 * (0.25 + 0.5) = 0.375, against 0.5 * 0.25 = 0.125 without the hold.
    CHECK(r.UN(1, 0) == doctest::Approx(0.375).epsilon(0.08));
    CHECK(r.UN(1, 0) > 0.25);
}

/**
 * The DISPATCHERS: a fork of two identical queues fed by one upstream station.
 *
 * WHY THIS MODEL. Every strategy here sends half the stream to each queue, so
 * the first moment alone cannot tell them apart: what differs is the ORDER, and
 * the order is visible only in the queue length. Round robin makes each branch's
 * arrival stream Erlang-2 rather than Poisson and its queue drops well below the
 * M/M/1 value; a weighted 3:1 split loads one branch far past the other;
 * join-the-shortest-queue balances them below both. A dispatcher silently read
 * as PROB reports the M/M/1 value in every one of these cases, which is exactly
 * the failure this file exists to catch.
 */
namespace {

qn::Network<double> dispatcher(lang::RoutingStrategy rs, double w1, double w2, int d) {
    qn::Network<double> m("dispatch");
    const std::size_t src = m.add_source("Source");
    const std::size_t q0 = m.add_queue("Q0", lang::SchedStrategy::FCFS);
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q0, c, lang::Distrib<double>::exp_rate(4.0));
    m.set_service(q1, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q2, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q0, 1.0);
    P.set(c, c, q0, q1, 0.5);
    P.set(c, c, q0, q2, 0.5);
    P.set(c, c, q1, snk, 1.0);
    P.set(c, c, q2, snk, 1.0);
    m.link(P);
    qn::NetworkStruct<double>& raw = m.raw_struct();
    raw.nodes[q0 - 1].routing.assign(1, rs);
    if (rs == lang::RoutingStrategy::WRROBIN) {
        raw.nodes[q0 - 1].routing_weights.assign(1, std::map<std::size_t, double>());
        raw.nodes[q0 - 1].routing_weights[0][q1] = w1;
        raw.nodes[q0 - 1].routing_weights[0][q2] = w2;
    }
    if (rs == lang::RoutingStrategy::SQ) raw.nodes[q0 - 1].routing_param.assign(1, d);
    return m;
}

}  // namespace

TEST_CASE("native LDES engine: PROB routing leaves each branch an M/M/1") {
    qn::Network<double> m = dispatcher(lang::RoutingStrategy::PROB, 0, 0, 0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    // Each branch sees a Poisson stream of rate 0.5 against a rate-1 server.
    CHECK(r.QN(2, 0) == doctest::Approx(1.0).epsilon(0.10));
    CHECK(r.QN(3, 0) == doctest::Approx(1.0).epsilon(0.10));
}

TEST_CASE("native LDES engine: round robin beats the probabilistic split") {
    qn::Network<double> m = dispatcher(lang::RoutingStrategy::RROBIN, 0, 0, 0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    // Alternating arrivals make each branch E2/M/1: strictly shorter than the
    // M/M/1 the same rate would give, and the gap is the whole content of the
    // strategy.
    CHECK(r.QN(2, 0) < 0.9);
    CHECK(r.QN(3, 0) < 0.9);
    CHECK(r.QN(2, 0) == doctest::Approx(r.QN(3, 0)).epsilon(0.10));
    // The split is still even in the first moment.
    CHECK(r.TN(2, 0) == doctest::Approx(r.TN(3, 0)).epsilon(0.05));
}

TEST_CASE("native LDES engine: weighted round robin loads by its weights") {
    qn::Network<double> m = dispatcher(lang::RoutingStrategy::WRROBIN, 3.0, 1.0, 0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    // Three of every four jobs take the first branch: 0.75 against 0.25 of the
    // rate-1 stream, so the first queue is heavily loaded and the second is not.
    CHECK(r.TN(2, 0) == doctest::Approx(0.75).epsilon(0.05));
    CHECK(r.TN(3, 0) == doctest::Approx(0.25).epsilon(0.05));
    CHECK(r.QN(2, 0) > 2.0);
    CHECK(r.QN(3, 0) < 0.5);
}

TEST_CASE("native LDES engine: join-the-shortest-queue balances below round robin") {
    qn::Network<double> m = dispatcher(lang::RoutingStrategy::JSQ, 0, 0, 0);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.QN(2, 0) < 0.85);
    CHECK(r.QN(3, 0) < 0.85);
    CHECK(r.QN(2, 0) == doctest::Approx(r.QN(3, 0)).epsilon(0.10));
    CHECK(r.TN(2, 0) == doctest::Approx(r.TN(3, 0)).epsilon(0.05));
}

TEST_CASE("native LDES engine: power-of-d with d at the fan-out is JSQ") {
    // Two destinations and d = 2: the sampling without replacement takes both,
    // so the choice is the shortest of all of them and the two strategies
    // coincide. A d that did not saturate would be strictly worse.
    qn::Network<double> m = dispatcher(lang::RoutingStrategy::SQ, 0, 0, 2);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.QN(2, 0) < 0.85);
    CHECK(r.QN(3, 0) < 0.85);
    CHECK(r.QN(2, 0) == doctest::Approx(r.QN(3, 0)).epsilon(0.10));
}

TEST_CASE("native LDES engine: joint dependence scales the service rate") {
    // eta(n) = n makes the station serve n jobs at total rate n, which IS the
    // infinite server: the mean number in system is the offered load 2 and does
    // NOT diverge, even though the nominal single-server rate is below lambda.
    // Declared but unread, this station would be an unstable M/M/1.
    qn::Network<double> m("jd");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    m.set_joint_dependence(q,
                           [](const std::vector<double>& n) {
                               double total = 0.0;
                               for (std::size_t k = 0; k < n.size(); ++k) total += n[k];
                               return std::vector<double>(1, total > 0.0 ? total : 1.0);
                           },
                           std::vector<double>(1, 1.0));
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.TN(1, 0) == doctest::Approx(2.0).epsilon(0.05));
    CHECK(r.QN(1, 0) == doctest::Approx(2.0).epsilon(0.10));
}

TEST_CASE("native LDES engine: class and joint dependence multiply") {
    // beta = 1/2 against eta = 2n leaves exactly the n of the case above, so
    // the same M/M/inf answer must come back. A port that applied one and
    // dropped the other would be off by the factor it dropped.
    qn::Network<double> m("cdjd");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    m.set_class_dependence(q,
                           [](const std::vector<double>&) {
                               return std::vector<double>(1, 0.5);
                           },
                           std::vector<double>(1, 1.0));
    m.set_joint_dependence(q,
                           [](const std::vector<double>& n) {
                               double total = 0.0;
                               for (std::size_t k = 0; k < n.size(); ++k) total += n[k];
                               return std::vector<double>(1, total > 0.0 ? 2.0 * total : 2.0);
                           },
                           std::vector<double>(1, 1.0));
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(300000));
    CHECK(r.TN(1, 0) == doctest::Approx(2.0).epsilon(0.05));
    CHECK(r.QN(1, 0) == doctest::Approx(2.0).epsilon(0.10));
}

TEST_CASE("native LDES engine: immediate feedback holds the server") {
    // A job that feeds back keeps the server instead of re-queueing, so its
    // return is NOT an arrival: AN counts the external stream alone. The queue
    // length is unchanged -- a single work-conserving server cares about the
    // work offered and not its order -- which is exactly why AN is the
    // measurement that separates the two models.
    qn::Network<double> off("imf");
    const std::size_t src = off.add_source("Source");
    const std::size_t q = off.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = off.add_sink("Sink");
    const std::size_t c = off.add_open_class("C1");
    off.set_arrival(src, c, lang::Distrib<double>::exp_rate(0.4));
    off.set_service(q, c, lang::Distrib<double>::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, q, 0.5);
    P.set(c, c, q, snk, 0.5);
    off.link(P);
    const ldes::LdesResult a = ldes::ldes_engine_solve(off.get_struct(), opts(200000));

    qn::Network<double> on = off;
    on.set_immediate_feedback(q, c);
    const ldes::LdesResult b = ldes::ldes_engine_solve(on.get_struct(), opts(200000));

    // Re-queueing counts every visit: the offered rate is lambda/(1-p) = 0.8.
    CHECK(a.AN(1, 0) == doctest::Approx(0.8).epsilon(0.05));
    // Feedback counts the external arrivals only.
    CHECK(b.AN(1, 0) == doctest::Approx(0.4).epsilon(0.05));
    // The work is the same either way, so the station is equally busy.
    CHECK(b.UN(1, 0) == doctest::Approx(a.UN(1, 0)).epsilon(0.05));
}

TEST_CASE("native LDES engine: a completion spawns its phase-2 continuation") {
    // Every completion of C1 injects a fresh C2 at the same station, so the
    // station serves TWO jobs per external arrival and its throughput in C2
    // equals its throughput in C1. Unread, C2 never appears at all.
    qn::Network<double> m("spawn");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(0.3));
    m.set_arrival(src, c2, lang::Distrib<double>::disabled_dist());
    m.set_service(q, c1, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q, c2, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);
    m.set_class_spawn(c1, c2);
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    CHECK(r.TN(1, 0) == doctest::Approx(0.3).epsilon(0.05));
    CHECK(r.TN(1, 1) == doctest::Approx(r.TN(1, 0)).epsilon(0.05));
    // Offered load 2 * 0.3 / 2 = 0.3 across the two classes.
    CHECK(r.UN(1, 0) + r.UN(1, 1) == doctest::Approx(0.3).epsilon(0.08));
}

TEST_CASE("native LDES engine: the warm start places the closed population") {
    // A machine repairman of four jobs, started with all four AT THE QUEUE
    // rather than at the reference Delay. The stationary answer does not depend
    // on where the run began, so the two placements must agree once the warmup
    // is truncated -- which is exactly what a warm start is for, and what a
    // silently ignored --initsol would also produce. The measurement that
    // separates them is the REFUSAL below.
    qn::Network<double> m("cqn4");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(200000);
    o.init_sol.push_back(0.0);  // Delay
    o.init_sol.push_back(4.0);  // Queue
    const ldes::LdesResult warm = ldes::ldes_engine_solve(m.get_struct(), o);
    const ldes::LdesResult cold = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    CHECK(warm.QN(0, 0) + warm.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-9));
    CHECK(warm.QN(1, 0) == doctest::Approx(cold.QN(1, 0)).epsilon(0.05));
}

TEST_CASE("native LDES engine: a warm start that loses jobs is refused") {
    // Three jobs declared for a population of four. A placement that quietly
    // dropped one would run a different model and report a confident answer for
    // it, and no mean in the table would say so.
    qn::Network<double> m("cqn4bad");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(10000);
    o.init_sol.push_back(1.0);
    o.init_sol.push_back(2.0);
    CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), o), InputError);
}

TEST_CASE("native LDES engine: the steady-state histogram recovers the marginal") {
    // The repairman's queue-length law is product form:
    //   P(n) = D^n Z^(N-n)/(N-n)! / G,  D = 0.5, Z = 1, N = 2, G = 1.25,
    // so P(0) = P(1) = 0.4 and P(2) = 0.2 EXACTLY. The histogram is a
    // residence-time law over the joint state, so the marginal read off it must
    // reproduce those numbers, and E[n] over it must be the QLen the means
    // report -- the two are the same measurement.
    qn::Network<double> m("cqn2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ldes::LdesOptions o = opts(300000);
    o.export_histogram = true;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    REQUIRE(r.histogram_space.rows() > 0);
    REQUIRE(r.histogram_time.rows() == r.histogram_space.rows());
    CHECK(r.histogram_space.cols() == 2);

    double total = 0.0, en = 0.0;
    std::vector<double> p(3, 0.0);
    for (std::size_t s = 0; s < r.histogram_time.rows(); ++s) {
        const double t = r.histogram_time(s, 0);
        total += t;
        // Every visited state holds the whole population.
        CHECK(r.histogram_space(s, 0) + r.histogram_space(s, 1) == doctest::Approx(2.0));
        const std::size_t n = static_cast<std::size_t>(r.histogram_space(s, 1) + 0.5);
        if (n < 3) p[n] += t;
        en += t * r.histogram_space(s, 1);
    }
    CHECK(total > 0.0);
    CHECK(p[0] / total == doctest::Approx(0.4).epsilon(0.05));
    CHECK(p[1] / total == doctest::Approx(0.4).epsilon(0.05));
    CHECK(p[2] / total == doctest::Approx(0.2).epsilon(0.08));
    CHECK(en / total == doctest::Approx(r.QN(1, 0)).epsilon(0.02));
}

TEST_CASE("native LDES engine: the response-time samples are the tally itself") {
    // M/M/1 at rho = 0.5: the mean of the recorded samples IS the reported RN,
    // and there is one sample per completion. A block filled from a different
    // measurement, or downsampled, would fail one of the two.
    qn::Network<double> m = mm1(0.5, 1.0);
    ldes::LdesOptions o = opts(50000);
    o.export_respt = true;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    REQUIRE(r.respTimeSamples.size() == r.nstations);
    const std::vector<double>& s = r.respTimeSamples[1][0];
    CHECK(s.size() == 50000u);
    double sum = 0.0;
    for (std::size_t k = 0; k < s.size(); ++k) sum += s[k];
    CHECK(sum / static_cast<double>(s.size()) == doctest::Approx(r.RN(1, 0)).epsilon(1e-9));
    CHECK(sum / static_cast<double>(s.size()) == doctest::Approx(2.0).epsilon(0.05));
}

TEST_CASE("native LDES engine: a retrieval system produces delayed hits") {
    // A DELAYED HIT IS ITS OWN OUTCOME. With a retrieval system a miss FETCHES
    // the item through a sub-network, and a second read of an item already in
    // flight is parked rather than starting a fetch of its own; it is released
    // when the fetch returns. Before this was implemented the engine reported
    // `delayed` empty and served the model as a plain cache, which is a silent
    // wrong number rather than a refusal.
    //
    // The fetch station is deliberately SLOW (rate 0.5 against an arrival rate
    // of 4): a long fetch keeps items in flight, so parked reads are common and
    // the delayed fraction is unmistakably positive.
    const std::vector<double> access{0.6, 0.3, 0.1};
    qn::Network<double> m("retr");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> cpar;
    cpar.nitems = access.size();
    cpar.itemcap.push_back(1);
    cpar.replacestrat = lang::ReplacementStrategy::FIFO;
    const std::size_t ch = m.add_cache("Cache", cpar);
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::INF);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t miss = m.add_open_class("Miss");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(4.0));
    m.set_service(q, c, lang::Distrib<double>::exp_rate(0.5));
    {
        qn::CacheParam<double>& cp = m.raw_struct().nodeparam[ch];
        cp.pread.assign(3, std::vector<double>());
        cp.pread[c - 1] = access;
        cp.hitclass.assign(3, 0);
        cp.missclass.assign(3, 0);
        cp.hitclass[c - 1] = hit;
        cp.missclass[c - 1] = miss;
    }
    m.set_retrieval_system(ch, c, miss, std::vector<std::size_t>{q});

    qn::RoutingMatrix<double> P;
    P.set(c, c, src, ch, 1.0);
    P.set(c, c, ch, q, 1.0);
    P.set(c, c, q, ch, 1.0);
    P.set(hit, hit, ch, snk, 1.0);
    P.set(miss, miss, ch, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    REQUIRE(r.cache_metrics.count("Cache") == 1);
    const ldes::LdesCacheMetrics& cm = r.cache_metrics.at("Cache");
    REQUIRE(cm.delayed.cols() > 0);
    const double h = cm.hit(0, c - 1), d = cm.delayed(0, c - 1), mi = cm.miss(0, c - 1);
    CHECK(d > 0.01);
    // The three outcomes PARTITION the reads: a delayed hit is neither a hit nor
    // a miss, so dividing by hits+misses alone would push the total above one.
    CHECK(h + d + mi == doctest::Approx(1.0).epsilon(1e-9));
    // Every fetch is measured, so the latency is reported and positive.
    REQUIRE(cm.latency.cols() > 0);
    CHECK(cm.latency(0, c - 1) > 0.0);
    // An EXPONENTIAL fetch period is memoryless, so a request parked partway
    // through one waits the same 2.0 on average as the fetch itself: this model
    // cannot tell the two averages apart, which is why the deterministic-fetch
    // case below exists.
    CHECK(cm.latency(0, c - 1) == doctest::Approx(2.0).epsilon(0.05));
}

TEST_CASE("native LDES engine: retrieval latency averages the parked waits too") {
    // THE REPORTED LATENCY IS A WAIT PER REQUEST, NOT PER FETCH. Two populations
    // wait on the retrieval system -- the request that triggered the fetch, for
    // the whole fetch period, and every request parked behind it, for its
    // RESIDUAL -- and `getAvgCacheTable` puts the mean beside ArvR =
    // lambda*(miss + delayed), the rate into that system. Little's law over that
    // pair is what fixes the average over BOTH, and it is what the JAR engine
    // (`Solver_ssj.exportCacheResults`) and the analytic `retrieval_fpi_latency`
    // both compute. Averaging the fetch sojourns alone reads high.
    //
    // A DETERMINISTIC FETCH IS WHAT MAKES THIS EXACTLY CHECKABLE. The retrieval
    // station is an infinite server with Det(D) service, so every fetch lasts
    // exactly D; the reads parked behind one arrive Poisson over that window, so
    // their epochs are uniform in it and they wait D/2 on average. With F
    // fetches and R releases the answer is therefore D*(F + R/2)/(F + R), which
    // the measured miss and delayed SHARES give as D*(m + d/2)/(m + d). The
    // fetch-only mean would be D exactly, so the two differ by a factor bounded
    // away from 1 as soon as any read is parked.
    const double D = 2.0;
    const std::vector<double> access{0.6, 0.3, 0.1};
    qn::Network<double> m("retr_det");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> cpar;
    cpar.nitems = access.size();
    cpar.itemcap.push_back(1);
    cpar.replacestrat = lang::ReplacementStrategy::FIFO;
    const std::size_t ch = m.add_cache("Cache", cpar);
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::INF);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t miss = m.add_open_class("Miss");
    m.set_arrival(src, c, lang::Distrib<double>::exp_rate(4.0));
    m.set_service(q, c, lang::Distrib<double>::det(D));
    {
        qn::CacheParam<double>& cp = m.raw_struct().nodeparam[ch];
        cp.pread.assign(3, std::vector<double>());
        cp.pread[c - 1] = access;
        cp.hitclass.assign(3, 0);
        cp.missclass.assign(3, 0);
        cp.hitclass[c - 1] = hit;
        cp.missclass[c - 1] = miss;
    }
    m.set_retrieval_system(ch, c, miss, std::vector<std::size_t>{q});

    qn::RoutingMatrix<double> P;
    P.set(c, c, src, ch, 1.0);
    P.set(c, c, ch, q, 1.0);
    P.set(c, c, q, ch, 1.0);
    P.set(hit, hit, ch, snk, 1.0);
    P.set(miss, miss, ch, snk, 1.0);
    m.link(P);

    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), opts(200000));
    REQUIRE(r.cache_metrics.count("Cache") == 1);
    const ldes::LdesCacheMetrics& cm = r.cache_metrics.at("Cache");
    REQUIRE(cm.delayed.cols() > 0);
    REQUIRE(cm.latency.cols() > 0);
    const double d = cm.delayed(0, c - 1), mi = cm.miss(0, c - 1);
    REQUIRE(d > 0.01);
    const double lat = cm.latency(0, c - 1);
    CHECK(lat == doctest::Approx(D * (mi + 0.5 * d) / (mi + d)).epsilon(0.02));
    // And it is strictly under the fetch period, which the fetch-only mean is not.
    CHECK(lat < D * 0.99);
}
