/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The RCAT analyzers of SolverMAM: INAP, INAP+ and the matrix-geometric
 * INAPINF (solver_ag.h), and the refusal of the optimization-based autocat
 * (solver_ag_autocat.h).
 *
 * THESE METHODS ARE APPROXIMATIONS and the oracles here are chosen so that
 * nothing pretends otherwise. RCAT is EXACT exactly when the reversed rate of
 * every synchronizing action is state-independent, which is what a product-form
 * model gives it. So the models below are deliberately the ones where the
 * product form is real -- an open tandem of exponential queues, i.e. a Jackson
 * network -- because there and only there is there a closed form to check
 * against at a tolerance that means anything. On a non-product-form model these
 * analyzers return marginals that are not the model's, and no test here claims
 * otherwise.
 *
 * WHAT MAKES THE TANDEM AN ORACLE RATHER THAN A REGRESSION. For
 * Source -> Q1 -> Q2 -> Sink the first component sees no passive action at all,
 * so its marginal is the truncated geometric in rho1 = lambda/mu1 whatever the
 * reversed rate iterate is. INAP then estimates the action rate as the mean of
 * A(i,j) pi(i)/pi(j) over the subdiagonal, every term of which is exactly
 * mu1 * rho1 = lambda; INAP+ estimates it as mu1 (1 - pi(0)), which is the
 * throughput. Both therefore land on lambda, and the second component becomes
 * the truncated geometric in lambda/mu2. Those two truncated geometrics are
 * written out independently in `trunc_geom` below from the arrival and service
 * rates alone, so the expected values never come from the analyzer.
 *
 * THE PHASE-TYPE CASES HAVE THEIR OWN ORACLE, and it is not a product form.
 * An isolated M/PH/1 carries no synchronizing action, so whatever the
 * reversed-rate iterate does its marginal is the M/G/1 one and its mean is the
 * Pollaczek-Khinchine value rho + rho^2 (1+scv) / (2 (1-rho)), written out in
 * `pk` below from the arrival rate and the SCV alone. That is what pins the
 * QBD construction -- the level/phase blocks of qbd_mapmap1 -- rather than the
 * fixed point, and it is what the scalar birth-death chain could not deliver:
 * before the phase dimension existed this analyzer returned the M/M/1 answer
 * for EVERY scv.
 *
 * THE TRUNCATION IS QUANTIFIED, NOT IGNORED. INAP and INAP+ cut every open
 * component at maxStates = 100 states, so their marginals differ from the
 * infinite M/M/1 by O(rho^100), which at the rho used here is below 1e-30 and
 * far under any tolerance asserted. INAPINF removes the truncation outright by
 * solving the scalar QBD root, and is checked against the exact untruncated
 * M/M/1 moments instead.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/ag/solver_ag.h"
#include "line/solvers/ag/solver_ag_autocat.h"
#include "line/solvers/ag/solver_ag_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Source -> Queue -> Sink, one open class: no synchronizing action exists. */
qn::Network<double> ag_mm1(double lambda, double mu) {
    qn::Network<double> m("agmm1");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, D::exp_rate(lambda));
    m.set_service(q, c, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> Q1 -> Q2 -> Sink: one action, and the RCAT product form is exact. */
qn::Network<double> ag_tandem(double lambda, double mu1, double mu2) {
    qn::Network<double> m("agtandem");
    const std::size_t s = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, D::exp_rate(lambda));
    m.set_service(q1, c, D::exp_rate(mu1));
    m.set_service(q2, c, D::exp_rate(mu2));
    qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    return m;
}

/**
 * Moments of the geometric law truncated to n states, pi_j = r^j / sum_k r^k.
 * Written from r and n alone, which is what makes it an oracle for the
 * truncated components INAP and INAP+ build.
 */
struct TruncGeom {
    double U, Q;
};
TruncGeom trunc_geom(double r, std::size_t n) {
    double Z = 0.0, q = 0.0, p = 1.0;
    for (std::size_t j = 0; j < n; ++j) {
        Z += p;
        q += static_cast<double>(j) * p;
        p *= r;
    }
    return {1.0 - 1.0 / Z, q / Z};
}

ag::AgResult<double> run(qn::Network<double>& m, const std::string& method) {
    ag::AgOptions opt;
    opt.method = method;
    return ag::solver_ag(m.get_struct(), opt);
}

TEST_CASE("mam-ag: an isolated open queue has no action and is the truncated M/M/1") {
    // Source and Sink are not RCAT components, so the only component is Q1 and
    // there is nothing to synchronize: the analyzer solves it once, with no
    // fixed point, and must report iter = 0.
    const double lambda = 0.5, mu = 2.0;
    qn::Network<double> m = ag_mm1(lambda, mu);
    const ag::AgResult<double> r = run(m, "inap");
    const TruncGeom g = trunc_geom(lambda / mu, 100);

    CHECK(r.actualmethod == "inap");
    CHECK(r.sol.iter == 0);
    // This branch goes through the general ctmc_solve rather than the
    // birth-death recursion, and the 100-state truncated geometric spans some
    // 60 orders of magnitude, so the tail entries are noise. The moments are
    // not: they are carried by the first few states. Hence 1e-10 here and
    // 1e-12 on the branches that use the exact recursion.
    CHECK(r.sol.U(1, 0) == doctest::Approx(g.U).epsilon(1e-10));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(g.Q).epsilon(1e-10));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(mu * g.U).epsilon(1e-10));
    // R is assembled from Little's law, so this pins the assembly, not the law.
    CHECK(r.sol.R(1, 0) == doctest::Approx(g.Q / (mu * g.U)).epsilon(1e-10));
    // An open class takes its system throughput from the Source rate directly.
    CHECK(r.sol.X[0] == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.C[0] == doctest::Approx(r.sol.R(1, 0)).epsilon(1e-12));
}

TEST_CASE("mam-ag: INAP recovers the exact reversed rate on a Jackson tandem") {
    // The reversed rate of the single action must come out at lambda, which
    // makes Q2 an M/M/1 fed at the source rate. Both marginals are then the
    // truncated geometrics below, and rate conservation across the two
    // stations holds to the truncation error, O(rho^100).
    const double lambda = 0.4, mu1 = 2.0, mu2 = 1.0;
    qn::Network<double> m = ag_tandem(lambda, mu1, mu2);
    const ag::AgResult<double> r = run(m, "inap");
    const TruncGeom g1 = trunc_geom(lambda / mu1, 100);
    const TruncGeom g2 = trunc_geom(lambda / mu2, 100);

    CHECK(r.sol.U(1, 0) == doctest::Approx(g1.U).epsilon(1e-12));
    CHECK(r.sol.U(2, 0) == doctest::Approx(g2.U).epsilon(1e-12));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(g1.Q).epsilon(1e-12));
    CHECK(r.sol.Q(2, 0) == doctest::Approx(g2.Q).epsilon(1e-12));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.Tp(2, 0) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.X[0] == doctest::Approx(lambda).epsilon(1e-12));
}

TEST_CASE("mam-ag: INAP+ conserves rate through the tandem") {
    // INAP+ estimates the action rate as mu1 (1 - pi1(0)), i.e. the throughput
    // of the upstream component, so rate conservation is what this estimator
    // asserts by construction and the downstream station must agree with it.
    const double lambda = 0.4, mu1 = 2.0, mu2 = 1.0;
    qn::Network<double> m = ag_tandem(lambda, mu1, mu2);
    const ag::AgResult<double> r = run(m, "inapplus");

    CHECK(r.actualmethod == "inapplus");
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.Tp(2, 0) == doctest::Approx(r.sol.Tp(1, 0)).epsilon(1e-9));
    CHECK(r.sol.U(1, 0) == doctest::Approx(trunc_geom(lambda / mu1, 100).U).epsilon(1e-12));
}

TEST_CASE("mam-ag: INAPINF drops the truncation and returns the exact M/M/1 marginals") {
    // Each open component is solved on its infinite state space by the scalar
    // QBD root, which for a birth-death component with birth f and death b and
    // no catastrophe drain is exactly f/b. So the marginals are the UNtruncated
    // geometrics and the moments are the textbook M/M/1 ones, with no 100-state
    // cut anywhere in them.
    const double lambda = 0.4, mu1 = 2.0, mu2 = 1.0;
    const double rho1 = lambda / mu1, rho2 = lambda / mu2;
    qn::Network<double> m = ag_tandem(lambda, mu1, mu2);
    const ag::AgResult<double> r = run(m, "inapinf");

    CHECK(r.actualmethod == "inapinf");
    CHECK(r.sol.U(1, 0) == doctest::Approx(rho1).epsilon(1e-12));
    CHECK(r.sol.U(2, 0) == doctest::Approx(rho2).epsilon(1e-12));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(rho1 / (1.0 - rho1)).epsilon(1e-12));
    CHECK(r.sol.Q(2, 0) == doctest::Approx(rho2 / (1.0 - rho2)).epsilon(1e-12));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.Tp(2, 0) == doctest::Approx(lambda).epsilon(1e-12));
    // Remark 2: the residual is zero exactly when the product form is real, and
    // on a Jackson tandem it is. This is the assertion that the method found a
    // genuine product form rather than merely a converged average.
    CHECK(r.rcat_residual < 1e-10);
}

TEST_CASE("mam-ag: the method aliases resolve as the reference does") {
    qn::Network<double> m = ag_tandem(0.4, 2.0, 1.0);
    // 'default' is INAP, and 'exact' is too: autocat left the tree, so the
    // reference warns and falls back rather than erroring.
    CHECK(run(m, "default").actualmethod == "inap");
    CHECK(run(m, "exact").actualmethod == "inap");
    // Aliasing, not an oracle: the same code path must give the same numbers.
    CHECK(run(m, "exact").sol.U(2, 0) == doctest::Approx(run(m, "inap").sol.U(2, 0)));
}

/** Source -> Q1 -> Q2 -> Sink with the given arrival and service laws. */
qn::Network<double> ag_tandem_dist(const D& arr, const D& s1, const D& s2) {
    qn::Network<double> m("agtd");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, arr);
    m.set_service(q1, c, s1);
    m.set_service(q2, c, s2);
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    return m;
}

/** Source -> Q1 -> Sink with the given arrival and service laws. */
qn::Network<double> ag_single(const D& arr, const D& svc) {
    qn::Network<double> m("agsd");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, arr);
    m.set_service(q, c, svc);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** M/G/1 mean number in system, the Pollaczek-Khinchine formula. */
double pk(double rho, double scv) {
    return rho + rho * rho * (1.0 + scv) / (2.0 * (1.0 - rho));
}

TEST_CASE("mam-ag: an M/PH/1 component is the exact Pollaczek-Khinchine mean") {
    // Erlang(k) of mean 1 has scv = 1/k. The component is a QBD of 100 levels
    // and k phases, and its first moment must be the M/G/1 one to the
    // truncation error, which at these rho is far below the tolerance.
    for (double rho : {0.5, 0.8})
        for (std::size_t k : {std::size_t(2), std::size_t(4)}) {
            const double scv = 1.0 / static_cast<double>(k);
            qn::Network<double> m =
                ag_single(D::exp_rate(rho), D::erlang_fit(1.0, scv));
            for (const std::string& meth : {std::string("inap"), std::string("inapinf")}) {
                const ag::AgResult<double> r = run(m, meth);
                CHECK(r.sol.Q(1, 0) == doctest::Approx(pk(rho, scv)).epsilon(1e-9));
                CHECK(r.sol.U(1, 0) == doctest::Approx(rho).epsilon(1e-9));
                // The throughput is the SERVICE completion rate averaged over
                // the marginal, not mu times P(busy): with phases those differ.
                CHECK(r.sol.Tp(1, 0) == doctest::Approx(rho).epsilon(1e-9));
            }
        }
}

TEST_CASE("mam-ag: one phase per level is the birth-death chain, unchanged") {
    // An Erlang(1) IS an exponential, so the phase machinery must reproduce the
    // scalar path entry for entry rather than merely approximate it.
    const double rho = 0.5;
    for (const std::string& meth :
         {std::string("inap"), std::string("inapplus"), std::string("inapinf")}) {
        qn::Network<double> me = ag_single(D::exp_rate(rho), D::exp_rate(1.0));
        qn::Network<double> m1 = ag_single(D::exp_rate(rho), D::erlang(1.0, 1));
        CHECK(run(m1, meth).sol.Q(1, 0) == doctest::Approx(run(me, meth).sol.Q(1, 0)).epsilon(1e-12));
    }
}

TEST_CASE("mam-ag: a phase-type arrival is not answered as a Poisson one") {
    // Er2/M/1 with mean interarrival 2 and mean service 1. The M/M/1 answer at
    // the same rho is exactly 1, so a Poisson reading of the Source is visible
    // here; the exact value is the one SolverCTMC and the MATLAB and Python
    // twins agree on.
    qn::Network<double> m = ag_single(D::erlang_fit(2.0, 0.5), D::exp_rate(1.0));
    for (const std::string& meth : {std::string("inap"), std::string("inapinf")})
        CHECK(run(m, meth).sol.Q(1, 0) == doctest::Approx(0.8090169943749).epsilon(1e-9));

    qn::Network<double> m2 =
        ag_single(lang::hyperexp_fit_mean_scv<double>(2.0, 4.0), D::exp_rate(1.0));
    for (const std::string& meth : {std::string("inap"), std::string("inapinf")})
        CHECK(run(m2, meth).sol.Q(1, 0) == doctest::Approx(1.13396527562).epsilon(1e-9));
}

TEST_CASE("mam-ag: a tandem with phase-type service conserves flow") {
    // The upstream component is an isolated M/PH/1 whatever the fixed point
    // does, so the reversed rate must come out at lambda and both stations must
    // carry it. This is the assertion the mean-of-ratios estimator fails on a
    // phase-expanded component -- it returned 1.27 against the exact 0.5 --
    // which is why such a component takes the rate-conserving one.
    const double lambda = 0.5;
    qn::Network<double> m =
        ag_tandem_dist(D::exp_rate(lambda), D::erlang_fit(1.0, 0.5), D::exp_rate(2.0));
    for (const std::string& meth :
         {std::string("inap"), std::string("inapplus"), std::string("inapinf")}) {
        const ag::AgResult<double> r = run(m, meth);
        CHECK(r.sol.Q(1, 0) == doctest::Approx(pk(lambda, 0.5)).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(lambda).epsilon(1e-9));
        CHECK(r.sol.Tp(2, 0) == doctest::Approx(lambda).epsilon(1e-9));
    }
}

TEST_CASE("mam-ag: INAPINF drops the truncation on a phase-expanded component") {
    // A HyperExp of scv 4 at rho 0.8 has a tail slow enough that 100 levels lose
    // more than a percent of the first moment. INAPINF solves the component by
    // Neuts' R instead, so it lands on the exact P-K mean.
    const double lambda = 0.8;
    qn::Network<double> m = ag_tandem_dist(D::exp_rate(lambda),
                                           lang::hyperexp_fit_mean_scv<double>(1.0, 4.0),
                                           D::exp_rate(4.0));
    const double exact = pk(lambda, 4.0);
    const double q_inap = run(m, "inap").sol.Q(1, 0);
    const double q_inf = run(m, "inapinf").sol.Q(1, 0);
    CHECK(q_inf == doctest::Approx(exact).epsilon(1e-9));
    CHECK(std::fabs(q_inap - exact) > std::fabs(q_inf - exact));
}

TEST_CASE("mam-ag: 'para' aliases 'parallel', and 'threads' is gone by name") {
    // The alias pair is SolverSSA's ({'para','parallel'} in
    // solver_ssa_analyzer.m, the same two in solver_ssa_parallel.h), so one
    // spelling convention covers both solvers rather than each inventing its own.
    CHECK(ag::exec_is_parallel(ag::exec_parallel()));
    CHECK(ag::exec_is_parallel(ag::exec_para()));
    CHECK_FALSE(ag::exec_is_parallel("threads"));
    CHECK_FALSE(ag::exec_is_parallel(ag::exec_serial()));

    // 'threads' was this backend's name until 2026-08-19. It must be refused BY
    // NAME: falling through to the generic message tells a caller with an old
    // script that a backend which still exists does not.
    qn::Network<double> m = ag_tandem(0.4, 2.0, 1.0);
    ag::AgOptions opt;
    opt.method = "inap";
    opt.exec = "threads";
    try {
        ag::solver_ag_solve(m.get_struct(), opt);
        FAIL("expected 'threads' to be refused");
    } catch (const InputError& e) {
        const std::string what = e.what();
        CHECK(what.find("renamed to 'parallel'") != std::string::npos);
    }
    // And the alias actually runs, giving the serial answer -- this backend is
    // asserted bit-identical to serial, being the same code in the same process.
    opt.exec = ag::exec_para();
    const ag::AgResult<double> par = ag::solver_ag_solve(m.get_struct(), opt);
    opt.exec = ag::exec_serial();
    const ag::AgResult<double> ser = ag::solver_ag_solve(m.get_struct(), opt);
    // Bit-for-bit, not approximately: this backend is the same code in the same
    // process, so `==` on the doubles is the assertion the header's claim earns.
    // doctest's Approx cannot express it -- even .epsilon(0) keeps a relative
    // tolerance -- so compare the values directly.
    CHECK(par.sol.Q(2, 0) == ser.sol.Q(2, 0));
    CHECK(par.sol.U(2, 0) == ser.sol.U(2, 0));
    CHECK(par.sol.Tp(2, 0) == ser.sol.Tp(2, 0));
}

TEST_CASE("mam-ag: what the RCAT analyzers refuse by name") {
    SUBCASE("a method that is not an RCAT variant") {
        qn::Network<double> m = ag_tandem(0.4, 2.0, 1.0);
        CHECK_THROWS_AS(run(m, "dec.source"), UnsupportedError);
    }
    SUBCASE("exact arithmetic, which a tolerance-stopped fixed point cannot serve") {
        qn::Network<Rational> m("agex");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(s, c, Distrib<Rational>::exp_rate(Rational(1, 2)));
        m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
        qn::RoutingMatrix<Rational> P;
        P.set(s, q, num_traits<Rational>::from_int(1));
        P.set(q, k, num_traits<Rational>::from_int(1));
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        CHECK_THROWS_AS(ag::solver_ag(m.get_struct(), opt), UnsupportedError);
    }
    SUBCASE("state-dependent routing, whose split is not a fixed action probability") {
        // build_rcat reads sn.rt ONCE and stores it as RcatAction::prob, but for
        // JSQ that entry is only a state-independent placeholder: the real split
        // depends on the queue lengths. Solving from the placeholder answers a
        // different model at a plausible number, so it is refused.
        //
        // The reference refuses it by FEATSET -- SolverAG.getFeatureSet declares
        // only RoutingStrategy_PROB and RoutingStrategy_RAND, and the JAR and
        // Python twins mirror that -- while this port has no ag featset, so until
        // 2026-08-19 C++ was the one codebase that accepted these and answered.
        qn::Network<double> m("agjsq");
        const std::size_t src = m.add_source("Source");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(src, c, D::exp_rate(1.0));
        m.set_service(q1, c, D::exp_rate(2.0));
        m.set_service(q2, c, D::exp_rate(2.0));
        m.set_routing(src, c, lang::RoutingStrategy::JSQ);
        qn::RoutingMatrix<double> P;
        P.set(src, q1, 0.5);
        P.set(src, q2, 0.5);
        P.set(q1, k, 1.0);
        P.set(q2, k, 1.0);
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        try {
            ag::solver_ag_solve(m.get_struct(), opt);
            FAIL("expected a refusal for state-dependent routing");
        } catch (const UnsupportedError& e) {
            const std::string what = e.what();
            // The EXPLICIT loop's wording, not the featset gate's: routing on a
            // SOURCE never reaches used_lang_features (its Source branch omits
            // routing, for MATLAB parity), so the declared-vs-used comparison
            // cannot see this and the dedicated check is what refuses it.
            CHECK(what.find("state-independent routing") != std::string::npos);
            CHECK(what.find("jsq") != std::string::npos);
        }
    }
    SUBCASE("a discipline the declared feature set withholds") {
        // The universal gate itself, which AG gained on 2026-08-19: ag_feature_set
        // had transcribed SolverAG.getFeatureSet since the port landed but nothing
        // ever called it, so AG was the one family here running ungated.
        //
        // SIRO is the probe because it is a SCHEDULING discipline -- the declared
        // set lists SchedStrategy_INF, _PS and _FCFS only -- and unlike a process
        // or a routing strategy no other check in this runner looks at it. So a
        // refusal here can only have come from feature_gate.
        qn::Network<double> m("agsiro");
        const std::size_t src = m.add_source("Source");
        const std::size_t q = m.add_queue("Q1", SchedStrategy::SIRO);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(src, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(src, q, 1.0);
        P.set(q, k, 1.0);
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        try {
            ag::solver_ag_solve(m.get_struct(), opt);
            FAIL("expected the feature gate to refuse an undeclared discipline");
        } catch (const UnsupportedError& e) {
            const std::string what = e.what();
            CHECK(what.find("SolverAG") != std::string::npos);
            CHECK(what.find("SIRO") != std::string::npos);
        }
    }
    SUBCASE("a binding finite buffer, which the decomposition cannot carry") {
        // Nothing under solvers/ag reads st.cap or st.classcap, so before this
        // gate a capped station was decomposed as an UNBOUNDED one and answered:
        // this very model returned the same figures with and without the cap
        // (QLen 1.09 in a buffer of 1, against the exact 0.5556).
        //
        // Lowering the agent's level bound would not have fixed it. nlev is the
        // class population by construction, and capping it at the buffer makes the
        // top-level boundary -- a SELF-LOOP -- lose the arrival, whereas a closed
        // job refused at a full buffer must BLOCK the upstream departure. That
        // coupling is exactly the independence RCAT assumes.
        qn::Network<double> m("agcap");
        const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("Class1", 2.0, q1);
        m.set_service(q1, c, D::exp_rate(1.0));
        m.set_service(q2, c, D::exp_rate(0.8));
        m.set_capacity(q2, 1.0);
        qn::RoutingMatrix<double> P;
        P.set(q1, q2, 1.0);
        P.set(q2, q1, 1.0);
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        try {
            ag::solver_ag_solve(m.get_struct(), opt);
            FAIL("expected SolverAG to refuse a binding finite capacity");
        } catch (const UnsupportedError& e) {
            const std::string what = e.what();
            CHECK(what.find("SolverAG") != std::string::npos);
            CHECK(what.find("capacity") != std::string::npos);
        }
    }
    SUBCASE("a capacity that cannot BIND, which stays solvable") {
        // setCapacity(N) on an N-job closed model is a no-op the gate must not
        // mistake for a buffer: only a capacity BELOW the reachable population is
        // refused.
        qn::Network<double> m("agnocap");
        const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("Class1", 2.0, q1);
        m.set_service(q1, c, D::exp_rate(1.0));
        m.set_service(q2, c, D::exp_rate(0.8));
        m.set_capacity(q2, 2.0);
        qn::RoutingMatrix<double> P;
        P.set(q1, q2, 1.0);
        P.set(q2, q1, 1.0);
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        CHECK_NOTHROW(ag::solver_ag_solve(m.get_struct(), opt));
    }
    SUBCASE("a matrix-exponential, which is not a generator at all") {
        // Moved here from tests/test_mam.cpp on 2026-08-19 with the RCAT methods.
        // A(0,1)=2 makes row 0 of A sum to +1, so A is no generator: a CTMC
        // assembled from it has a signed stationary vector, and answering would
        // report numbers that are not a distribution. Left unguarded the analyzer
        // does answer -- Q=97.5, U=1.0 on this model -- so the gate is the whole
        // protection.
        //
        // Through `solver_ag_solve`, NOT the `run` helper above: the process gate
        // lives in solver_ag_runner (check_model_method, the twin of
        // SolverAG.supportsModelMethod), and `run` calls the analyzer directly and
        // so bypasses every gate. A refusal case has to enter where the gate is.
        Matrix<double> A(2, 2, 0.0);
        A(0, 0) = -1.0;
        A(0, 1) = 2.0;
        A(1, 1) = -3.0;
        qn::Network<double> m("agme");
        const std::size_t src = m.add_source("Source");
        const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(src, c, D::exp_rate(1.0));
        m.set_service(q, c, D::me({1.0, 0.0}, A));
        qn::RoutingMatrix<double> P;
        P.set(src, q, 1.0);
        P.set(q, k, 1.0);
        m.link(P);
        ag::AgOptions opt;
        opt.method = "inap";
        try {
            ag::solver_ag_solve(m.get_struct(), opt);
            FAIL("expected a refusal for a matrix-exponential service process");
        } catch (const UnsupportedError& e) {
            CHECK(std::string(e.what()).find("Markovian") != std::string::npos);
        }
    }
    SUBCASE("autocat, which is dead upstream and needs an LP/NLP stack") {
        qn::Network<double> m = ag_tandem(0.4, 2.0, 1.0);
        ag::AgOptions opt;
        opt.method = "exact";
        CHECK_THROWS_AS(ag::solver_ag_autocat(m.get_struct(), opt), UnsupportedError);
    }
}

}  // namespace
