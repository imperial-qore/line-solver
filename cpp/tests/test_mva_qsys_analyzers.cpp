/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The two open M/G/1 analyzers the dispatch carries -- HOL priority
 * (mva_dispatch.h, from solver_mva_qsys_prio_analyzer.m) and the size-based
 * family (mva_dispatch.h, from solver_mva_qsys_sizebased_analyzer.m) -- driven
 * directly, against oracles that do NOT come from the implementations.
 *
 * WHY THIS FILE IS SEPARATE FROM test_mva_dispatch.cpp. That file pins the
 * ROUTING -- which analyzer a model reaches -- and carries MATLAB reference
 * numbers. For the size-based branch its own comment records that its expected
 * values were read back out of the C++ analyzer, because SolverMVA's featset
 * does not list SRPT/PSJF/FB/LRPT/SETF and the branch is unreachable from the
 * MATLAB solver. A read-back number pins a regression but cannot detect that
 * the port was wrong on the day it landed. Asserted here instead are laws the
 * answer must satisfy whatever the closed form computes.
 *
 * 2026-07-29: mva_feature_set was narrowed to MATLAB's then-59 declared names,
 * so the five size-based disciplines became undeclared in the C++ gate too and
 * the branch was unreachable from SolverMVA in either codebase.
 *
 * 2026-08-16: reversed. MATLAB declares all five, and the JAR and native Python
 * gained the analyzer, so the gate admits them again and the cases at the FOOT
 * of this file drive the whole route end to end against MATLAB numbers. The
 * cases below still call mva::solver_mva_qsys_sizebased_analyzer directly and
 * assert laws rather than numbers, which is the check a reference value cannot
 * make. The laws are:
 *
 *   - a single priority class collapses Cobham to Pollaczek-Khinchine, hence
 *     to the already-verified qsys_mm1 for exponential service;
 *   - Kleinrock's conservation law: with class-independent service HOL is work
 *     conserving and blind to the split, so the AGGREGATE queue length is the
 *     M/M/1 value at the pooled arrival rate and depends on neither the load
 *     split nor the priority order (Kleinrock, Queueing Systems II, Sec. 3.4);
 *   - Little's law per class at the queue;
 *   - Schrage's optimality of SRPT and pessimality of LRPT among work-
 *     conserving M/G/1 disciplines, which bracket the other three.
 *
 * THE GLUE, NOT THE FORMULAS, IS UNDER TEST. The closed forms already have
 * MATLAB-pinned unit tests in test_qsys_mg1_disciplines.cpp. What only these
 * can catch is the sn-to-api layer: the priority permutation and its inverse
 * on writeback, the per-class chain lookup for the visit ratio, and the
 * assembly of Q/U/R/T/C/X. The permutation is silent when wrong -- every
 * number is a legal Cobham value, merely on the wrong class -- so the ordering
 * case below declares the classes in reverse of their priority order.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/mva_dispatch.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/**
 * Source -> Queue -> Sink with one open class per (arrival, service) pair, each
 * in its own chain. `prios` is the declared class priority, lowest served first.
 */
qn::Network<double> open_multiclass(const std::string& name, SchedStrategy sched,
                                    const std::vector<double>& lambda,
                                    const std::vector<double>& mu,
                                    const std::vector<int>& prios) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t k = m.add_sink("Sink");
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 0; r < lambda.size(); ++r) {
        const std::size_t c = m.add_open_class("C" + std::to_string(r), prios[r]);
        m.set_arrival(s, c, D::exp_rate(lambda[r]));
        m.set_service(q, c, D::exp_rate(mu[r]));
        P.set(c, c, s, q, 1.0);
        P.set(c, c, q, k, 1.0);
    }
    m.link(P);
    return m;
}

/** The dispatch reports the algorithm alongside the metrics, so this IS its result. */
template <class T>
using Ran = mva::DispatchResult<T>;

template <class T>
Ran<T> prio_of(qn::Network<T>& m) {
    return mva::solver_mva_qsys_prio_analyzer(m.get_struct(), mva::MvaOptions());
}

template <class T>
Ran<T> sizebased_of(qn::Network<T>& m) {
    return mva::solver_mva_qsys_sizebased_analyzer(m.get_struct(), mva::MvaOptions());
}

TEST_CASE("qsys-analyzers: one HOL class collapses to the verified M/M/1 path") {
    // The K=1 reduction is what fixes the whole assembly -- R, Q, U, X and both
    // throughputs -- against an independent path in one shot.
    const double lambda = 0.6, mu = 2.0;
    qn::Network<double> m = open_multiclass("hol1", SchedStrategy::HOL, {lambda}, {mu}, {0});
    const Ran<double> r = prio_of(m);
    const double W = qsys::qsys_mm1(lambda, mu).W;

    CHECK(r.actualmethod == "mg1.prio");
    CHECK(r.sol.R(1, 0) == doctest::Approx(W).epsilon(1e-12));
    CHECK(r.sol.C[0] == doctest::Approx(W).epsilon(1e-12));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(lambda * W).epsilon(1e-12));
    CHECK(r.sol.U(1, 0) == doctest::Approx(lambda / mu).epsilon(1e-12));
    CHECK(r.sol.X[0] == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.Tp(0, 0) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-12));
}

TEST_CASE("qsys-analyzers: HOL priority obeys Kleinrock's conservation law") {
    // Each case holds mu = 2 and the pooled rate at 0.8, so all five must
    // return the same total 0.4/0.6 = 2/3 while the per-class values differ
    // widely. The last two share a load vector and differ only in priority
    // order, which isolates the permutation from everything else.
    struct Case {
        std::vector<double> lambda;
        std::vector<int> prio;
    };
    const Case cases[] = {
        {{0.4, 0.4}, {0, 1}},
        {{0.1, 0.7}, {0, 1}},
        {{0.7, 0.1}, {0, 1}},
        {{0.2, 0.3, 0.3}, {0, 1, 2}},
        {{0.2, 0.3, 0.3}, {2, 0, 1}},
    };
    const double mu = 2.0, pooled = 0.8;
    const double Qtot_ref = qsys::qsys_mm1(pooled, mu).W * pooled;

    for (const Case& c : cases) {
        const std::vector<double> mus(c.lambda.size(), mu);
        qn::Network<double> m =
            open_multiclass("holcons", SchedStrategy::HOL, c.lambda, mus, c.prio);
        const Ran<double> r = prio_of(m);
        double Qtot = 0.0, Utot = 0.0;
        for (std::size_t k = 0; k < c.lambda.size(); ++k) {
            Qtot += r.sol.Q(1, k);
            Utot += r.sol.U(1, k);
            // Q and R are written from separate expressions, so per-class
            // Little's law is a real constraint on the assembly.
            CHECK(r.sol.Q(1, k) ==
                  doctest::Approx(r.sol.Tp(1, k) * r.sol.R(1, k)).epsilon(1e-12));
            CHECK(r.sol.Tp(0, k) == doctest::Approx(c.lambda[k]).epsilon(1e-12));
        }
        CHECK(Qtot == doctest::Approx(Qtot_ref).epsilon(1e-9));
        CHECK(Utot == doctest::Approx(pooled / mu).epsilon(1e-12));
    }
}

TEST_CASE("qsys-analyzers: HOL results follow the declared priority, not the class order") {
    // Classes declared in REVERSE priority order, so a port that dropped the
    // inverse permutation returns the same two numbers with columns swapped
    // and still passes every conservation law above.
    //
    // Cobham by hand for lambda = 0.4 each, mu = 2, exponential (cs = 1):
    //   rho_i = 0.2, B_0 = (1/2) sum lambda_i (1+1)/mu_i^2 = 0.2
    //   W_high = 0.2/(1   * 0.8) + 0.5 = 3/4
    //   W_low  = 0.2/(0.8 * 0.6) + 0.5 = 11/12
    qn::Network<double> m =
        open_multiclass("holperm", SchedStrategy::HOL, {0.4, 0.4}, {2.0, 2.0}, {5, 1});
    const Ran<double> r = prio_of(m);

    CHECK(r.sol.R(1, 0) == doctest::Approx(11.0 / 12.0).epsilon(1e-12));  // prio 5, served last
    CHECK(r.sol.R(1, 1) == doctest::Approx(3.0 / 4.0).epsilon(1e-12));    // prio 1, served first
    CHECK(r.sol.Q(1, 0) == doctest::Approx(11.0 / 30.0).epsilon(1e-12));
    CHECK(r.sol.Q(1, 1) == doctest::Approx(3.0 / 10.0).epsilon(1e-12));
}

TEST_CASE("qsys-analyzers: HOL priority stays exact under Rational arithmetic") {
    // Unlike its qsys neighbours the Cobham formula is rational in its inputs,
    // and for exponential service sqrt(scv) = 1 is exact, so this analyzer is
    // deliberately NOT gated on transcendental arithmetic and must return the
    // exact rationals derived by hand above, not a rounded double lifted into
    // Rational.
    qn::Network<Rational> m("holexact");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::HOL);
    const std::size_t k = m.add_sink("Sink");
    const Rational two = num_traits<Rational>::from_int(2);
    const Rational lam = Rational(2, 5);
    const std::size_t hi = m.add_open_class("Hi", 0);
    const std::size_t lo = m.add_open_class("Lo", 1);
    qn::RoutingMatrix<Rational> P;
    for (std::size_t c : {hi, lo}) {
        m.set_arrival(s, c, Distrib<Rational>::exp_rate(lam));
        m.set_service(q, c, Distrib<Rational>::exp_rate(two));
        P.set(c, c, s, q, num_traits<Rational>::from_int(1));
        P.set(c, c, q, k, num_traits<Rational>::from_int(1));
    }
    m.link(P);
    const Ran<Rational> r = prio_of(m);

    CHECK(r.actualmethod == "mg1.prio");
    CHECK(r.sol.R(1, 0) == Rational(3) / Rational(4));
    CHECK(r.sol.R(1, 1) == Rational(11) / Rational(12));
    CHECK(r.sol.Q(1, 0) == Rational(3) / Rational(10));
    CHECK(r.sol.Q(1, 1) == Rational(11) / Rational(30));
    // The conservation law again, now with no rounding to hide behind.
    CHECK(r.sol.Q(1, 0) + r.sol.Q(1, 1) == Rational(2) / Rational(3));
}

TEST_CASE("qsys-analyzers: the size-based disciplines bracket between SRPT and LRPT") {
    // Two classes with DIFFERENT mean sizes, which is what makes a size-based
    // discipline differ from FCFS at all, and which exercises the per-class
    // chain lookup for the visit ratio: the reference read the visit of the
    // chain whose number equalled the Source's STATION index, so every class
    // outside chain 1 saw a visit of 0 and the analyzer rejected its own
    // model. A throw here, not a wrong number, is that regression.
    //
    // Schrage (1968): SRPT minimizes the mean response time over all work-
    // conserving M/G/1 disciplines and LRPT maximizes it. With the arrival
    // rates fixed, Little's law turns that into a bracket on the aggregate
    // queue length, so the five are comparable without knowing any value.
    const std::vector<double> lambda = {0.3, 0.3}, mu = {2.0, 1.0};
    const struct Case {
        SchedStrategy sched;
        const char* actual;
    } cases[] = {
        {SchedStrategy::SRPT, "mg1.srpt"}, {SchedStrategy::PSJF, "mg1.psjf"},
        {SchedStrategy::FB, "mg1.fb"},     {SchedStrategy::SETF, "mg1.setf"},
        {SchedStrategy::LRPT, "mg1.lrpt"},
    };
    double Qtot[5] = {0, 0, 0, 0, 0};

    for (std::size_t i = 0; i < 5; ++i) {
        qn::Network<double> m =
            open_multiclass("sb", cases[i].sched, lambda, mu, std::vector<int>(2, 0));
        const Ran<double> r = sizebased_of(m);
        CHECK(r.actualmethod == cases[i].actual);
        for (std::size_t k = 0; k < 2; ++k) {
            CHECK(r.sol.Tp(0, k) == doctest::Approx(lambda[k]).epsilon(1e-12));
            CHECK(r.sol.Tp(1, k) == doctest::Approx(lambda[k]).epsilon(1e-12));
            CHECK(r.sol.X[k] == doctest::Approx(lambda[k]).epsilon(1e-12));
            CHECK(r.sol.U(1, k) == doctest::Approx(lambda[k] / mu[k]).epsilon(1e-12));
            // With a unit visit ratio R is W, so Little's law must hold; the
            // feedback case, where the reference's own assembly would break
            // it, is never dispatched to this branch.
            CHECK(r.sol.Q(1, k) == doctest::Approx(r.sol.X[k] * r.sol.R(1, k)).epsilon(1e-12));
            CHECK(r.sol.C[k] == doctest::Approx(r.sol.R(1, k)).epsilon(1e-12));
            Qtot[i] += r.sol.Q(1, k);
        }
        // Work conserving, so none may report a queue shorter than the
        // utilization it also reports.
        CHECK(Qtot[i] > lambda[0] / mu[0] + lambda[1] / mu[1]);
    }
    for (std::size_t i = 1; i < 5; ++i) {
        CHECK(Qtot[0] <= Qtot[i]);  // SRPT is optimal
        CHECK(Qtot[4] >= Qtot[i]);  // LRPT is pessimal
    }
}

TEST_CASE("qsys-analyzers: the size-based analyzer refuses what it cannot cover") {
    SUBCASE("a discipline that is not size based is named in the refusal") {
        qn::Network<double> m =
            open_multiclass("sbfcfs", SchedStrategy::FCFS, {0.3, 0.3}, {2.0, 1.0}, {0, 0});
        CHECK_THROWS_AS(sizebased_of(m), UnsupportedError);
    }
    SUBCASE("a non-positive service rate is an input error, not a silent zero") {
        qn::Network<double> m =
            open_multiclass("sbzero", SchedStrategy::SRPT, {0.3, 0.3}, {2.0, 1.0}, {0, 0});
        qn::NetworkStruct<double> L = m.get_struct();
        L.rates(1, 1) = 0.0;
        CHECK_THROWS_AS(mva::solver_mva_qsys_sizebased_analyzer(L, mva::MvaOptions()), InputError);
    }
    SUBCASE("exact arithmetic, which the quadrature-based closed forms cannot serve") {
        qn::Network<Rational> m("sbexact");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::SRPT);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t a = m.add_open_class("A");
        const std::size_t b = m.add_open_class("B");
        qn::RoutingMatrix<Rational> P;
        for (std::size_t c : {a, b}) {
            m.set_arrival(s, c, Distrib<Rational>::exp_rate(Rational(3, 10)));
            m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
            P.set(c, c, s, q, num_traits<Rational>::from_int(1));
            P.set(c, c, q, k, num_traits<Rational>::from_int(1));
        }
        m.link(P);
        CHECK_THROWS_AS(sizebased_of(m), UnsupportedError);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// The size-based family THROUGH THE GATE (2026-08-16)
// ---------------------------------------------------------------------------
//
// The cases above call solver_mva_qsys_sizebased_analyzer directly and assert
// laws rather than numbers, because when they were written the branch was
// unreachable from SolverMVA in both codebases. That has changed: MATLAB's
// getFeatureSet declares the five, mva_feature_set declares them again, and the
// JAR and native Python gained the same analyzer on 2026-08-16. So the whole
// route -- feature_gate, resolve_method, mva_dispatch branch 3, the closed form
// and the writeback -- can now be driven end to end and pinned to MATLAB.
//
// The fixture is the one the JAR (MvaSizeBasedTest) and Python
// (test_mva_sizebased.py) use, so the three ports are compared against ONE set
// of reference numbers: C1 Exp(1) at rate 0.3, C2 Erlang(mean 2, SCV 0.5) at
// rate 0.2, hence rho = 0.7 with the classes differing in both mean size and
// variability, which is what the size-based orderings act on.

namespace {

qn::Network<double> sizebased_fixture(SchedStrategy sched) {
    qn::Network<double> m("sb");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1", 0);
    const std::size_t c2 = m.add_open_class("C2", 0);
    m.set_arrival(s, c1, D::exp_rate(0.3));
    m.set_arrival(s, c2, D::exp_rate(0.2));
    m.set_service(q, c1, D::exp_rate(1.0));
    m.set_service(q, c2, D::erlang(1.0, 2));  // mean 2, SCV 0.5
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, s, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, s, q, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);
    return m;
}

/** The full solver route, gate included. */
mva::AvgResult<double> sizebased_run(qn::Network<double>& m) {
    const qn::NetworkStruct<double>& L = m.get_struct();
    qn::feature_gate("SolverMVA", qn::mva_feature_set("default"), L);
    return mva::solver_mva_run_analyzer(L, mva::MvaOptions(), Matrix<double>());
}

void check_sizebased(SchedStrategy sched, double q1, double q2, double r1, double r2) {
    qn::Network<double> m = sizebased_fixture(sched);
    const mva::AvgResult<double> r = sizebased_run(m);
    CHECK(r.QN(1, 0) == doctest::Approx(q1).epsilon(1e-9));
    CHECK(r.QN(1, 1) == doctest::Approx(q2).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(r1).epsilon(1e-9));
    CHECK(r.RN(1, 1) == doctest::Approx(r2).epsilon(1e-9));
    // discipline-independent: rho_k = lambda_k/mu_k, and the flow is lossless
    CHECK(r.UN(1, 0) == doctest::Approx(0.3).epsilon(1e-12));
    CHECK(r.UN(1, 1) == doctest::Approx(0.4).epsilon(1e-12));
    CHECK(r.TN(1, 0) == doctest::Approx(0.3).epsilon(1e-12));
    CHECK(r.TN(1, 1) == doctest::Approx(0.2).epsilon(1e-12));
    // Little's law at the station
    CHECK(r.QN(1, 0) == doctest::Approx(r.TN(1, 0) * r.RN(1, 0)).epsilon(1e-9));
    CHECK(r.QN(1, 1) == doctest::Approx(r.TN(1, 1) * r.RN(1, 1)).epsilon(1e-9));
}

}  // namespace

TEST_CASE("qsys-sizebased: SRPT through the gate matches MATLAB") {
    check_sizebased(SchedStrategy::SRPT, 0.4985567765823074, 0.8212546141301067,
                    1.6618559219410247, 4.1062730706505333);
}

TEST_CASE("qsys-sizebased: PSJF through the gate matches MATLAB") {
    check_sizebased(SchedStrategy::PSJF, 0.61224489795918369, 3.333333333333333,
                    2.0408163265306123, 16.666666666666664);
}

TEST_CASE("qsys-sizebased: FB through the gate matches MATLAB") {
    check_sizebased(SchedStrategy::FB, 0.63587346812077949, 2.1712143660569709,
                    2.1195782270692649, 10.856071830284854);
}

TEST_CASE("qsys-sizebased: LRPT through the gate matches MATLAB") {
    check_sizebased(SchedStrategy::LRPT, 2.1333333333333333, 0.66666666666666674,
                    7.1111111111111107, 3.3333333333333335);
}

TEST_CASE("qsys-sizebased: SETF through the gate matches MATLAB") {
    check_sizebased(SchedStrategy::SETF, 1.2256856111717471, 2.8758520284392772,
                    4.0856187039058236, 14.379260142196385);
}

TEST_CASE("qsys-sizebased: the gate no longer refuses the five disciplines") {
    // The regression this guards is the featset, not the arithmetic: the branch
    // and the closed forms were present and reachable code held behind a closed
    // gate for two and a half weeks.
    for (SchedStrategy s : {SchedStrategy::SRPT, SchedStrategy::PSJF, SchedStrategy::FB,
                            SchedStrategy::LRPT, SchedStrategy::SETF}) {
        qn::Network<double> m = sizebased_fixture(s);
        CHECK_NOTHROW(
            qn::feature_gate("SolverMVA", qn::mva_feature_set("default"), m.get_struct()));
    }
}
