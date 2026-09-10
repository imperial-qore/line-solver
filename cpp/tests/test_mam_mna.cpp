/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `solver_mna_open` and `solver_mna_closed`: the two analyzers behind
 * SolverMAM's `mna` method, and every refusal they make, asserted by name.
 *
 * THE ORACLES, in decreasing order of independence.
 *
 *  1. CLOSED FORMS the analyzer has no way to know about. A single open FCFS
 *     queue is fed an arrival stream whose fitted SCV is EXACTLY 1 -- the
 *     source's own SCV reaches it through a splitting factor of exactly
 *     1 + 1*(1-1) -- so the station solve is an exact M/G/1 and its mean queue
 *     length must be Pollaczek-Khinchine's. The finite-buffer branch must
 *     likewise reproduce the M/M/1/K stationary law, and the exponential
 *     tandem must reproduce Burke's theorem, because Whitt's departure-SCV
 *     formula returns 1 for an M/M/1 and APHFrom2Moments at cv2 = 1 returns an
 *     order-2 REPRESENTATION of the exponential (its transform collapses to
 *     lambda/(s+lambda)), so the downstream queue is fed genuine Poisson.
 *  2. IDENTITIES the analyzer does not impose: flow balance across a tandem,
 *     the utilization law U = T S / c, and the closed analyzer's population
 *     conservation.
 *  3. MATLAB, for the cases with no closed form. Every such value is labelled
 *     with the exact call that produced it, run against
 *     matlab/src/solvers/MAM/solver_mna_{open,closed}.m at
 *     SolverMAM.defaultOptions with options.method = 'mna'.
 *
 * The analyzers are called DIRECTLY rather than through `mam_dispatch`: the
 * dispatch tries the exact MAP/MAP/1 fast path ahead of any method, so an
 * M/M/1 asked for by method name never reaches MNA at all.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mna.h"

using namespace line;
using Dd = lang::Distrib<double>;
using lang::SchedStrategy;

namespace {

/** Source -> Queue -> Sink, one open class. */
qn::Network<double> mna_sq(const std::string& name, SchedStrategy sched, const Dd& arrival,
                           const Dd& service, double servers = 1.0, double cap = -1.0) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, arrival);
    m.set_service(q, o, service);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    if (cap > 0.0) m.set_capacity(q, cap);
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Closed Delay -> Queue -> Delay, one class. */
qn::Network<double> mna_dq(const std::string& name, double N, const Dd& think, const Dd& service,
                           SchedStrategy sched = SchedStrategy::FCFS, double servers = 1.0) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", sched);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    const std::size_t c = m.add_closed_class("C1", N, d);
    m.set_service(d, c, think);
    m.set_service(q, c, service);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

mva::MvaSolution<double> run_open(qn::Network<double>& m,
                                  const mam::MnaConfig& cfg = mam::MnaConfig()) {
    return mam::solver_mna_open(m.get_struct(), mam::MamOptions(), cfg);
}

mva::MvaSolution<double> run_closed(qn::Network<double>& m) {
    return mam::solver_mna_closed(m.get_struct(), mam::MamOptions());
}

const double kEps = 1e-9;

}  // namespace

// ---------------------------------------------------------------------------
// The open analyzer against closed forms
// ---------------------------------------------------------------------------

TEST_CASE("mna open: a single M/M/1 queue is exact") {
    // lambda = 1, mu = 2: E[N] = rho/(1-rho) = 1, R = 1/(mu-lambda) = 1.
    qn::Network<double> m = mna_sq("mnaA", SchedStrategy::FCFS, Dd::exp_rate(1.0),
                                   Dd::exp_rate(2.0));
    const mva::MvaSolution<double> r = run_open(m);
    CHECK(r.method == "mna");
    CHECK(r.Q(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.U(1, 0) == doctest::Approx(0.5).epsilon(kEps));
    CHECK(r.R(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.Tp(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.Tp(0, 0) == doctest::Approx(1.0).epsilon(kEps));
    // The Source carries no queue of its own.
    CHECK(r.Q(0, 0) == doctest::Approx(0.0));
    CHECK(r.C[0] == doctest::Approx(1.0).epsilon(kEps));
    // THE REFERENCE NEVER ASSIGNS X. See solver_mna.h; both analyzers return
    // the zeros they initialised it with, whatever the model.
    CHECK(r.X[0] == doctest::Approx(0.0));
    // Two sweeps: the first fills a1 and a2 from zero, the second reproduces
    // them exactly because an M/M/1's departure SCV is its arrival SCV.
    CHECK(r.iter == 2);
}

TEST_CASE("mna open: a single M/E2/1 queue is the Pollaczek-Khinchine mean") {
    // lambda = 1, Erlang-2 service of mean 0.5 and SCV 0.5, so rho = 0.5 and
    // L = rho + rho^2 (1 + Cs^2) / (2 (1 - rho)) = 0.5 + 0.375 = 0.875.
    qn::Network<double> m = mna_sq("mnaB", SchedStrategy::FCFS, Dd::exp_rate(1.0),
                                   Dd::erlang_fit(0.5, 0.5));
    const mva::MvaSolution<double> r = run_open(m);
    const double rho = 0.5, cs2 = 0.5;
    const double L = rho + rho * rho * (1.0 + cs2) / (2.0 * (1.0 - rho));
    CHECK(L == doctest::Approx(0.875));
    CHECK(r.Q(1, 0) == doctest::Approx(L).epsilon(kEps));
    CHECK(r.R(1, 0) == doctest::Approx(L).epsilon(kEps));  // Little at lambda = 1
    CHECK(r.U(1, 0) == doctest::Approx(0.5).epsilon(kEps));
}

TEST_CASE("mna open: an exponential tandem recovers Burke's theorem") {
    // Source(1) -> M/M/1(mu=2) -> M/M/1(mu=4) -> Sink. Both queues must report
    // their own M/M/1 answer, 1 and 1/3, because the departure stream of the
    // first is Poisson. MNA reaches that only because Whitt's formula returns
    // d2 = 1 for an M/M/1; the residue is the FineTol the reference adds to
    // every service rate before dividing.
    qn::Network<double> m("mnaC");
    const std::size_t s = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Dd::exp_rate(1.0));
    m.set_service(q1, o, Dd::exp_rate(2.0));
    m.set_service(q2, o, Dd::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    const mva::MvaSolution<double> r = run_open(m);

    CHECK(r.Q(1, 0) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(r.Q(2, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-8));
    // Flow balance and the utilization law, neither of which the sweep imposes.
    CHECK(r.Tp(1, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.Tp(2, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.U(1, 0) == doctest::Approx(r.Tp(1, 0) * 0.5).epsilon(kEps));
    CHECK(r.U(2, 0) == doctest::Approx(r.Tp(2, 0) * 0.25).epsilon(kEps));
    // MATLAB solver_mna_open reports 0.333333333666667 for the second queue,
    // i.e. the exact 1/3 perturbed at the ninth digit by that FineTol.
    CHECK(r.Q(2, 0) == doctest::Approx(0.333333333666667).epsilon(1e-13));
    CHECK(r.C[0] == doctest::Approx(1.33333333366667).epsilon(1e-12));
}

TEST_CASE("mna open: a finite buffer takes the exact M/M/1/K law") {
    // lambda = 0.6, mu = 1, K = 4. The M/M/1/K stationary law is
    // p_n = rho^n / sum_j rho^j, and the branch must report its mean and its
    // carried throughput lambda (1 - p_K).
    const double lambda = 0.6, mu = 1.0;
    const int K = 4;
    qn::Network<double> m = mna_sq("mnaD", SchedStrategy::FCFS, Dd::exp_rate(lambda),
                                   Dd::exp_rate(mu), 1.0, K);
    const mva::MvaSolution<double> r = run_open(m);

    const double rho = lambda / mu;
    double norm = 0.0, meanN = 0.0;
    for (int i = 0; i <= K; ++i) norm += std::pow(rho, i);
    for (int i = 0; i <= K; ++i) meanN += i * std::pow(rho, i) / norm;
    const double pK = std::pow(rho, K) / norm;

    CHECK(r.Q(1, 0) == doctest::Approx(meanN).epsilon(kEps));
    CHECK(r.Tp(1, 0) == doctest::Approx(lambda * (1.0 - pK)).epsilon(kEps));
    CHECK(r.U(1, 0) == doctest::Approx(lambda * (1.0 - pK) / mu).epsilon(kEps));
    // Little's law on the CARRIED stream, which the branch reconstructs from a
    // single class-independent Wq rather than reading off the queue solver.
    CHECK(r.Q(1, 0) == doctest::Approx(r.Tp(1, 0) * r.R(1, 0)).epsilon(kEps));
    // The Source still reports the OFFERED rate: nothing in MNA drops it.
    CHECK(r.Tp(0, 0) == doctest::Approx(lambda).epsilon(kEps));
}

TEST_CASE("mna open: a delay station passes Poisson through untouched") {
    // Source(1) -> Delay(mu=4) -> M/M/1(mu=2) -> Sink. The delay reports
    // U = Q = T S = 0.25 and leaves the downstream queue an exact M/M/1.
    qn::Network<double> m("mnaE");
    const std::size_t s = m.add_source("Src");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, Dd::exp_rate(1.0));
    m.set_service(d, o, Dd::exp_rate(4.0));
    m.set_service(q, o, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, d, 1.0);
    P.set(d, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const mva::MvaSolution<double> r = run_open(m);
    CHECK(r.Q(1, 0) == doctest::Approx(0.25).epsilon(kEps));
    CHECK(r.U(1, 0) == doctest::Approx(0.25).epsilon(kEps));
    CHECK(r.R(1, 0) == doctest::Approx(0.25).epsilon(kEps));
    CHECK(r.Q(2, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r.C[0] == doctest::Approx(1.25).epsilon(kEps));
}

TEST_CASE("mna open: a PS station reports nothing, as the reference's does") {
    // solver_mna_open.m's PS branch assigns to TN/UN/QN/RN, which are fresh
    // undefined variables rather than the T/U/Q/R it returns. The station keeps
    // the throughput the first sweep seeded and reports no queue at all.
    // Measured: MATLAB solver_mna_open returns Q = [0;0], U = [0;0], R = [0;0],
    // T = [1;1] for Source(1) -> PS(mu=2) -> Sink.
    qn::Network<double> m = mna_sq("mnaF", SchedStrategy::PS, Dd::exp_rate(1.0),
                                   Dd::exp_rate(2.0));
    const mva::MvaSolution<double> r = run_open(m);
    CHECK(r.Q(1, 0) == doctest::Approx(0.0));
    CHECK(r.U(1, 0) == doctest::Approx(0.0));
    CHECK(r.R(1, 0) == doctest::Approx(0.0));
    CHECK(r.C[0] == doctest::Approx(0.0));
    // The seeded throughput survives, which is what makes the omission visible
    // rather than a station that was simply never reached.
    CHECK(r.Tp(1, 0) == doctest::Approx(1.0).epsilon(kEps));
}

TEST_CASE("mna open: two classes at one queue against MATLAB") {
    // Source(0.4, 0.6) -> M/M/1(mu = 2, 1) -> Sink. This is NOT the marked
    // Poisson M/G/1: the sweep hands each class its own two-moment fit, and the
    // superposition formula leaves class r with a2 = lambda_r / lambda, so both
    // fitted streams are hypoexponential rather than Poisson. Values are
    // MATLAB solver_mna_open's.
    qn::Network<double> m("mnaG");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(s, c1, Dd::exp_rate(0.4));
    m.set_arrival(s, c2, Dd::exp_rate(0.6));
    m.set_service(q, c1, Dd::exp_rate(2.0));
    m.set_service(q, c2, Dd::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t r : {c1, c2}) {
        P.set(r, r, s, q, 1.0);
        P.set(r, r, q, k, 1.0);
    }
    m.link(P);
    const mva::MvaSolution<double> r = run_open(m);
    CHECK(r.Q(1, 0) == doctest::Approx(1.28585435733447).epsilon(1e-9));
    CHECK(r.Q(1, 1) == doctest::Approx(2.19531355930531).epsilon(1e-9));
    CHECK(r.R(1, 0) == doctest::Approx(3.21463589333618).epsilon(1e-9));
    CHECK(r.R(1, 1) == doctest::Approx(3.65885593217552).epsilon(1e-9));
    // The utilization law holds per class, and the total is the true rho = 0.8.
    CHECK(r.U(1, 0) == doctest::Approx(0.2).epsilon(kEps));
    CHECK(r.U(1, 1) == doctest::Approx(0.6).epsilon(kEps));
    CHECK(r.Tp(1, 0) == doctest::Approx(0.4).epsilon(kEps));
    CHECK(r.Tp(1, 1) == doctest::Approx(0.6).epsilon(kEps));
}

// ---------------------------------------------------------------------------
// The closed analyzer
// ---------------------------------------------------------------------------

TEST_CASE("mna closed: the Delay+Queue bisection against MATLAB") {
    // Z = 1, S = 0.5. Values are MATLAB solver_mna_closed's, which stops the
    // outer bisection after 16 halvings of the bracket [0, min mu] = [0, 2].
    qn::Network<double> m2 = mna_dq("mnaH2", 2.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const mva::MvaSolution<double> r2 = run_closed(m2);
    CHECK(r2.iter == 16);
    CHECK(r2.Q(0, 0) == doctest::Approx(1.12310526541279).epsilon(1e-9));
    CHECK(r2.Q(1, 0) == doctest::Approx(0.876894734587207).epsilon(1e-9));
    CHECK(r2.U(0, 0) == doctest::Approx(1.12310526541279).epsilon(1e-9));
    CHECK(r2.U(1, 0) == doctest::Approx(0.561553955078125).epsilon(1e-9));
    CHECK(r2.R(0, 0) == doctest::Approx(1.0).epsilon(kEps));
    CHECK(r2.R(1, 0) == doctest::Approx(0.780776977539062).epsilon(1e-9));
    CHECK(r2.Tp(0, 0) == doctest::Approx(1.12310791015625).epsilon(1e-9));
    CHECK(r2.Tp(1, 0) == doctest::Approx(1.12310791015625).epsilon(1e-9));
    CHECK(r2.C[0] == doctest::Approx(1.78077697753906).epsilon(1e-9));
    CHECK(r2.X[0] == doctest::Approx(0.0));

    qn::Network<double> m4 = mna_dq("mnaH4", 4.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    const mva::MvaSolution<double> r4 = run_closed(m4);
    CHECK(r4.iter == 16);
    CHECK(r4.Q(0, 0) == doctest::Approx(1.6089165944223).epsilon(1e-9));
    CHECK(r4.Q(1, 0) == doctest::Approx(2.3910834055777).epsilon(1e-9));
    CHECK(r4.U(1, 0) == doctest::Approx(0.804473876953125).epsilon(1e-9));
    CHECK(r4.R(1, 0) == doctest::Approx(1.48614503316516).epsilon(1e-9));
    CHECK(r4.Tp(1, 0) == doctest::Approx(1.60894775390625).epsilon(1e-9));

    // Erlang-2 service, so the station solve is a genuine MMAP/PH/1 rather than
    // a Markovian one. N = 3, mean 0.5, SCV 0.5.
    qn::Network<double> me = mna_dq("mnaHE", 3.0, Dd::exp_rate(1.0), Dd::erlang_fit(0.5, 0.5));
    const mva::MvaSolution<double> re = run_closed(me);
    CHECK(re.iter == 16);
    CHECK(re.Q(0, 0) == doctest::Approx(1.51602766853582).epsilon(1e-9));
    CHECK(re.Q(1, 0) == doctest::Approx(1.48397233146418).epsilon(1e-9));
    CHECK(re.U(1, 0) == doctest::Approx(0.758026123046875).epsilon(1e-9));
    CHECK(re.R(1, 0) == doctest::Approx(0.978855704459141).epsilon(1e-9));

    // Two servers, N = 5, mu = 1.5.
    qn::Network<double> ms =
        mna_dq("mnaHS", 5.0, Dd::exp_rate(1.0), Dd::exp_rate(1.5), SchedStrategy::FCFS, 2.0);
    const mva::MvaSolution<double> rs = run_closed(ms);
    CHECK(rs.iter == 17);
    CHECK(rs.Q(0, 0) == doctest::Approx(1.34586656666781).epsilon(1e-9));
    CHECK(rs.Q(1, 0) == doctest::Approx(3.65413343333219).epsilon(1e-9));
    CHECK(rs.U(1, 0) == doctest::Approx(0.448616027832031).epsilon(1e-9));
    CHECK(rs.R(1, 0) == doctest::Approx(2.71507854034843).epsilon(1e-9));
}

TEST_CASE("mna closed: population conservation and the utilization law") {
    // Neither is imposed by the station solves: conservation comes from the
    // terminal renormalization and the utilization law from T and S alone, so
    // a station solve that drifted would still have to satisfy both.
    for (double N : {1.0, 2.0, 4.0, 8.0}) {
        qn::Network<double> m = mna_dq("mnaI", N, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
        const mva::MvaSolution<double> r = run_closed(m);
        CHECK(r.Q(0, 0) + r.Q(1, 0) == doctest::Approx(N).epsilon(1e-12));
        CHECK(r.U(1, 0) == doctest::Approx(r.Tp(1, 0) * 0.5).epsilon(1e-12));
        // A delay reports its queue length as its utilization.
        CHECK(r.U(0, 0) == doctest::Approx(r.Q(0, 0)).epsilon(1e-12));
        // The bisection never leaves its bracket [0, min mu] = [0, 2], so the
        // queue can never be driven past saturation.
        CHECK(r.Tp(1, 0) <= 2.0 + 1e-12);
        CHECK(r.U(1, 0) <= 1.0 + 1e-12);
    }
}

TEST_CASE("mna closed: the PS branch stalls at the throughput bracket") {
    // With Z = 1 and S = 0.5 the bracket's upper end mu = 2 already drives
    // U = 1 exactly, so the geometric bound (U - U^(N+1))/(1-Uden) is 0/1e-8 =
    // 0, the population target is never met, and the bisection halves an
    // interval whose ends have both become 2. The reference runs to iter_max
    // and reports the delay holding every job. Values are MATLAB's.
    qn::Network<double> m = mna_dq("mnaJ", 3.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0),
                                   SchedStrategy::PS);
    const mva::MvaSolution<double> r = run_closed(m);
    CHECK(r.iter == 100);  // options.iter_max for SolverMAM
    CHECK(r.Q(0, 0) == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(r.Q(1, 0) == doctest::Approx(0.0));
    CHECK(r.U(0, 0) == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(r.U(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.Tp(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.Tp(1, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.C[0] == doctest::Approx(1.0).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

namespace {

/** Runs `f` and asserts it refused with `needle` in the message. */
template <class F>
void refuses(F f, const std::string& needle) {
    try {
        f();
        FAIL("expected a refusal naming: ", needle);
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find(needle) != std::string::npos);
    }
}

}  // namespace

TEST_CASE("mna: each analyzer refuses the other's population regime by name") {
    qn::Network<double> closed = mna_dq("mnaK1", 2.0, Dd::exp_rate(1.0), Dd::exp_rate(2.0));
    refuses([&] { run_open(closed); }, "belongs to solver_mna_closed");

    qn::Network<double> open = mna_sq("mnaK2", SchedStrategy::FCFS, Dd::exp_rate(1.0),
                                      Dd::exp_rate(2.0));
    refuses([&] { run_closed(open); }, "belongs to solver_mna_open");
}

TEST_CASE("mna: the dead etaqa departure-SCV option is refused by name") {
    qn::Network<double> m = mna_sq("mnaL", SchedStrategy::FCFS, Dd::exp_rate(1.0),
                                   Dd::exp_rate(2.0));
    mam::MnaConfig cfg;
    cfg.dep_scv = "etaqa";
    refuses([&] { run_open(m, cfg); }, "qbd_depproc_jointmom");
    cfg.dep_scv = "whitt";
    refuses([&] { run_open(m, cfg); }, "unknown config.dep_scv");
}

TEST_CASE("mna: a Fork-Join model is refused by name") {
    qn::Network<double> m("mnaM");
    const std::size_t d = m.add_delay("D");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dd::exp_rate(1.0));
    m.set_service(q1, c, Dd::exp_rate(2.0));
    m.set_service(q2, c, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    refuses([&] { run_closed(m); }, "Fork nodes are not supported yet");
}

TEST_CASE("mna: exact arithmetic is refused by name") {
    qn::Network<Rational> m("mnaN");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("C1");
    m.set_arrival(s, o, lang::Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, o, lang::Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(s, q, num_traits<Rational>::from_int(1));
    P.set(q, k, num_traits<Rational>::from_int(1));
    m.link(P);
    const qn::NetworkStruct<Rational>& L = m.get_struct();
    refuses([&] { mam::solver_mna_open(L, mam::MamOptions()); }, "--arith double or --arith real");
    refuses([&] { mam::solver_mna_closed(L, mam::MamOptions()); },
            "--arith double or --arith real");
}
