/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `solver_mvac_analyzer` (solver_mvac.h) and `solver_sqd` (solver_mva.h).
 *
 * THEY ARE NEVER COMPARED WITH EACH OTHER, and not only because they take
 * different models. MVAC is EXACT for its class -- it returns what the classic
 * population recursion returns, by a chain recursion that shares no code with
 * it -- so it is checked against exact MVA and against closed forms at the
 * tolerance of an identity. SQD is an APPROXIMATION, the only handler here that
 * represents blocking at all, and there is no exact solver in this port to
 * compare it against on a blocking model; it is checked against the
 * conservation law it must satisfy whatever its error, and its MATLAB numbers
 * are pinned separately in test_mva_special.cpp.
 *
 * Oracles are closed forms or cross-solver identities. Nothing here is a value
 * read back out of the implementation under test.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/mva/solver_mvac.h"

using namespace line;
using lang::Distrib;
using lang::DropStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/**
 * Two identical single-server PS queues in a closed cycle, no think time.
 * Balanced by construction, which is what gives the closed form its symmetry.
 */
template <class T>
qn::Network<T> balanced_pair(double njobs, const T& rate) {
    qn::Network<T> m("balanced");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, q1);
    m.set_service(q1, c, Distrib<T>::exp_rate(rate));
    m.set_service(q2, c, Distrib<T>::exp_rate(rate));
    qn::RoutingMatrix<T> P;
    P.set(c, c, q1, q2, num_traits<T>::from_int(1));
    P.set(c, c, q2, q1, num_traits<T>::from_int(1));
    m.link(P);
    return m;
}

/** Delay -> PS -> PS closed cycle, two classes, one class per chain. */
qn::Network<double> mvac_cycle(double n1, double n2) {
    qn::Network<double> m("mvaccycle");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const double Sd[2] = {1.0, 0.5}, S1[2] = {0.6, 0.3}, S2[2] = {0.25, 0.4};
    const double pop[2] = {n1, n2};
    qn::RoutingMatrix<double> P;
    for (std::size_t j = 0; j < 2; ++j) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(j + 1), pop[j], d);
        m.set_service(d, c, D::exp_rate(1.0 / Sd[j]));
        m.set_service(q1, c, D::exp_rate(1.0 / S1[j]));
        m.set_service(q2, c, D::exp_rate(1.0 / S2[j]));
        P.set(c, c, d, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

/** Delay -> FCFS(cap) -> FCFS(cap) closed cycle with blocking after service. */
qn::Network<double> bas_cycle(double njobs, double cap1, double cap2) {
    qn::Network<double> m("bascycle");
    const std::size_t d = m.add_delay("Think");
    const std::size_t b1 = m.add_queue("B1", SchedStrategy::FCFS);
    const std::size_t b2 = m.add_queue("B2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, D::exp_rate(2.0));
    m.set_service(b1, c, D::exp_rate(1.5));
    m.set_service(b2, c, D::exp_rate(2.5));
    m.set_capacity(b1, cap1);
    m.set_drop_rule(b1, c, DropStrategy::BAS);
    m.set_capacity(b2, cap2);
    m.set_drop_rule(b2, c, DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(d, b1, 1.0);
    P.set(b1, b2, 1.0);
    P.set(b2, d, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> mva_run(qn::Network<double>& m, const std::string& method = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

mva::MvaSolution<double> mvac_direct(qn::Network<double>& m) {
    return mva::solver_mvac_analyzer(m.get_struct(), mva::MvaOptions());
}

/**
 * The wired SQD entry point takes the chain demands the analyzer has already
 * built, rather than rebuilding them, so the caller supplies them here too.
 */
template <class T>
mva::MvaSolution<T> sqd_direct(qn::Network<T>& m) {
    const qn::NetworkStruct<T>& sn = m.get_struct();
    return mva::solver_sqd(sn, mva::sn_get_demands_chain(sn));
}

}  // namespace

// ---------------------------------------------------------------------------
// solver_mvac.h
// ---------------------------------------------------------------------------

TEST_CASE("mva sqd/mvac: MVAC and the population recursion agree on a product form") {
    // Both are exact on this model and share no code: MVAC recurs on the chains
    // with self-looping single-customer replacements, exact MVA on the
    // population vector. Any gap is a defect in one of them, not a tolerance.
    qn::Network<double> m = mvac_cycle(3.0, 2.0);
    const mva::MvaSolution<double> a = mvac_direct(m);
    const mva::AvgResult<double> e = mva_run(m, "exact");
    const double N[2] = {3.0, 2.0};
    for (std::size_t j = 0; j < 2; ++j) {
        double n = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            CHECK(a.Q(i, j) == doctest::Approx(e.QN(i, j)).epsilon(1e-9));
            CHECK(a.U(i, j) == doctest::Approx(e.UN(i, j)).epsilon(1e-9));
            CHECK(a.R(i, j) == doctest::Approx(e.RN(i, j)).epsilon(1e-9));
            CHECK(a.Tp(i, j) == doctest::Approx(e.TN(i, j)).epsilon(1e-9));
            n += a.Q(i, j);
        }
        CHECK(a.X[j] == doctest::Approx(e.XN[j]).epsilon(1e-9));
        // a closed class holds its whole population somewhere
        CHECK(n == doctest::Approx(N[j]).epsilon(1e-9));
    }
    CHECK(std::isnan(a.lG));  // MVAC forms no normalizing constant
}

TEST_CASE("mva sqd/mvac: MVAC reproduces the balanced two-station closed form") {
    // Two identical single servers, no think time. By symmetry Q_i = N/2, and
    // the MVA recursion then gives a cycle time of D(N+1), so X = N/(D(N+1))
    // and U_i = N/(N+1). This also exercises the path with NO delay station,
    // where the api is handed a zero think time rather than a folded one.
    const double Dq = 0.5, N = 4.0;
    qn::Network<double> m = balanced_pair<double>(N, 1.0 / Dq);
    const mva::MvaSolution<double> r = mvac_direct(m);
    const double X = N / (Dq * (N + 1.0));
    CHECK(r.X[0] == doctest::Approx(X).epsilon(1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.Q(i, 0) == doctest::Approx(N / 2.0).epsilon(1e-12));
        CHECK(r.Tp(i, 0) == doctest::Approx(X).epsilon(1e-12));
        CHECK(r.U(i, 0) == doctest::Approx(N / (N + 1.0)).epsilon(1e-12));
        CHECK(r.R(i, 0) == doctest::Approx(Dq * (1.0 + (N - 1.0) / 2.0)).epsilon(1e-12));
    }
}

TEST_CASE("mva sqd/mvac: MVAC is exact-capable arithmetic") {
    // The header claims the recursion forms no normalizing constant and no
    // logarithm, and carries no arithmetic gate; under Rational it must return
    // the closed form as an exact fraction. With D = 1/2 and N = 4,
    // X = N/(D(N+1)) = 8/5 and Q_i = 2, both exactly representable.
    qn::Network<Rational> m =
        balanced_pair<Rational>(4.0, num_traits<Rational>::from_int(2));
    const mva::MvaSolution<Rational> r =
        mva::solver_mvac_analyzer(m.get_struct(), mva::MvaOptions());
    CHECK(r.X[0] == Rational(8) / Rational(5));
    CHECK(r.Q(0, 0) == Rational(2));
    CHECK(r.Q(1, 0) == Rational(2));
}

TEST_CASE("mva sqd/mvac: MVAC refuses every model outside its class") {
    SUBCASE("a multiserver queue, which the SSFR arrival theorem has no term for") {
        qn::Network<double> m("mvacms");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
        m.set_number_of_servers(q, 2.0);
        const std::size_t c = m.add_closed_class("C1", 3.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(mvac_direct(m), UnsupportedError);
    }
    SUBCASE("an open class, which has no population to recur on") {
        qn::Network<double> m("mvacopen");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("O1");
        m.set_arrival(s, c, D::exp_rate(0.5));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(s, q, 1.0);
        P.set(q, k, 1.0);
        m.link(P);
        CHECK_THROWS_AS(mvac_direct(m), UnsupportedError);
    }
    SUBCASE("a discipline outside the product form") {
        qn::Network<double> m("mvacdps");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", SchedStrategy::DPS);
        const std::size_t c = m.add_closed_class("C1", 2.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        m.set_sched_param(q, c, 1.0);
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(mvac_direct(m), UnsupportedError);
    }
    SUBCASE("a model of delays only, which gives the recursion nothing to recur on") {
        qn::Network<double> m("mvacdelays");
        const std::size_t d1 = m.add_delay("Z1");
        const std::size_t d2 = m.add_delay("Z2");
        const std::size_t c = m.add_closed_class("C1", 2.0, d1);
        m.set_service(d1, c, D::exp_rate(1.0));
        m.set_service(d2, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d1, d2, 1.0);
        P.set(d2, d1, 1.0);
        m.link(P);
        CHECK_THROWS_AS(mvac_direct(m), UnsupportedError);
    }
    SUBCASE("a fractional population, which cannot be removed one customer at a time") {
        qn::Network<double> m = mvac_cycle(2.5, 1.0);
        CHECK_THROWS_AS(mvac_direct(m), UnsupportedError);
    }
}

// ---------------------------------------------------------------------------
// solver_sqd (solver_mva.h)
// ---------------------------------------------------------------------------

TEST_CASE("mva sqd/mvac: blocked jobs are still in the network") {
    // The conservation law is the one identity the approximation cannot trade
    // away: a job held at a full downstream buffer occupies its upstream server
    // but has not left the network, so the queue lengths still sum to the
    // population. The tolerance is the decomposition's own, not machine
    // precision, and no value here is compared against an exact solver -- there
    // is none for a blocking model.
    for (double n : {3.0, 5.0}) {
        qn::Network<double> m = bas_cycle(n, 2.0, 2.0);
        const mva::MvaSolution<double> r = sqd_direct(m);
        double tot = 0.0;
        for (std::size_t i = 0; i < 3; ++i) tot += r.Q(i, 0);
        CHECK(tot == doctest::Approx(n).epsilon(1e-6));
        // and the throughput is one circulating flow: every station sees it
        CHECK(r.Tp(1, 0) == doctest::Approx(r.Tp(0, 0)).epsilon(1e-9));
        CHECK(r.Tp(2, 0) == doctest::Approx(r.Tp(0, 0)).epsilon(1e-9));
    }
}

TEST_CASE("mva sqd/mvac: SQD refuses a model it has no single population for") {
    SUBCASE("a multichain model") {
        qn::Network<double> m("bastwochain");
        const std::size_t d = m.add_delay("Think");
        const std::size_t b = m.add_queue("B1", SchedStrategy::FCFS);
        m.set_capacity(b, 2.0);
        qn::RoutingMatrix<double> P;
        for (std::size_t j = 0; j < 2; ++j) {
            const std::size_t c = m.add_closed_class("C" + std::to_string(j + 1), 2.0, d);
            m.set_service(d, c, D::exp_rate(1.0));
            m.set_service(b, c, D::exp_rate(2.0));
            m.set_drop_rule(b, c, DropStrategy::BAS);
            P.set(c, c, d, b, 1.0);
            P.set(c, c, b, d, 1.0);
        }
        m.link(P);
        // The reference warns and returns NaN here; this port refuses by name,
        // so the caller cannot carry a NaN onward as if it were a result.
        CHECK_THROWS_AS(sqd_direct(m), UnsupportedError);
    }
    SUBCASE("an open class") {
        qn::Network<double> m("basopen");
        const std::size_t s = m.add_source("Source");
        const std::size_t b = m.add_queue("B1", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("O1");
        m.set_arrival(s, c, D::exp_rate(0.5));
        m.set_service(b, c, D::exp_rate(2.0));
        m.set_capacity(b, 2.0);
        m.set_drop_rule(b, c, DropStrategy::BAS);
        qn::RoutingMatrix<double> P;
        P.set(s, b, 1.0);
        P.set(b, k, 1.0);
        m.link(P);
        CHECK_THROWS_AS(sqd_direct(m), UnsupportedError);
    }
}

TEST_CASE("mva sqd/mvac: SQD refuses exact arithmetic by name") {
    // The effective-rate calibration evaluates exp, log and real powers, so
    // there is no exact-field version of the answer. Contrast MVAC above, which
    // is exact-capable and runs under the same backend.
    qn::Network<Rational> m =
        balanced_pair<Rational>(2.0, num_traits<Rational>::from_int(2));
    CHECK_THROWS_AS(sqd_direct(m), UnsupportedError);
}
