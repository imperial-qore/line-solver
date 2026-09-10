/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The two MVA analyzers that are not mean-value analysis, both in solver_mva.h:
 * `solver_mva_sum` (SUM/ESUM) and `solver_mva_lcfsqn` (the closed LCFS +
 * LCFS-PR pair).
 *
 * THE TWO HAVE DIFFERENT ERROR CLASSES AND ARE NEVER COMPARED WITH EACH OTHER.
 * LCFS-QN is EXACT for its model class, so it is checked against the NC
 * analyzer of the SAME model class, which reaches the answer through permanents
 * and shares no code with the recursion, at the tolerance of an identity. SUM is
 * an APPROXIMATION with no error bound, so it is checked only where it is
 * provably exact (a population that cannot queue) or against identities it must
 * satisfy whatever its error: the population constraint it roots, Little's law
 * and the utilization law.
 *
 * The MATLAB numbers for both are pinned separately, in test_mva_special.cpp
 * and test_mva_dispatch.cpp, through `mva::solver_mva_run_analyzer`. What these add is
 * the laws and the cross-solver identities that a recorded value cannot express,
 * plus the gate each analyzer sits behind.
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
#include "line/solvers/nc/solver_nc_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** LCFS -> LCFS-PR closed tandem; `mu[i][r]` is the service RATE at station i. */
template <class T>
qn::Network<T> lcfs_pair(const std::vector<std::vector<double>>& mu,
                         const std::vector<double>& pop) {
    qn::Network<T> m("lcfspair");
    const std::size_t q1 = m.add_queue("LCFS", SchedStrategy::LCFS);
    const std::size_t q2 = m.add_queue("LCFSPR", SchedStrategy::LCFSPR);
    qn::RoutingMatrix<T> P;
    for (std::size_t r = 0; r < pop.size(); ++r) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(r + 1), pop[r], q1);
        m.set_service(q1, c, Distrib<T>::exp_rate(num_traits<T>::from_double(mu[0][r])));
        m.set_service(q2, c, Distrib<T>::exp_rate(num_traits<T>::from_double(mu[1][r])));
        P.set(c, c, q1, q2, num_traits<T>::from_int(1));
        P.set(c, c, q2, q1, num_traits<T>::from_int(1));
    }
    m.link(P);
    return m;
}

/** Delay -> FCFS(2 servers) -> PS closed cycle, one class per chain. */
qn::Network<double> sum_cycle(const std::vector<double>& pop, double scv1) {
    qn::Network<double> m("sumcycle");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    m.set_number_of_servers(q1, 2.0);
    const double Sd[2] = {1.0, 0.5}, S1[2] = {0.8, 0.4}, S2[2] = {1.0 / 3.0, 0.25};
    qn::RoutingMatrix<double> P;
    for (std::size_t j = 0; j < pop.size(); ++j) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(j + 1), pop[j], d);
        m.set_service(d, c, D::exp_rate(1.0 / Sd[j]));
        m.set_service(q1, c, D::erlang_fit(S1[j], scv1));
        m.set_service(q2, c, D::exp_rate(1.0 / S2[j]));
        P.set(c, c, d, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

mva::AvgResult<double> mva_run(qn::Network<double>& m, const std::string& method = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

mva::AvgResult<double> nc_run(qn::Network<double>& m) {
    nc::NcSolverOptions opt;
    return nc::solver_nc_run_analyzer(m.get_struct(), opt);
}

/**
 * The LCFS analyzer entered directly, below the gate `solver_mva` applies. The
 * station pair is located here rather than assumed, so a helper that built the
 * model differently cannot silently pass the wrong indices.
 */
template <class T>
mva::MvaSolution<T> lcfsqn_direct(qn::Network<T>& m) {
    const qn::NetworkStruct<T>& sn = m.get_struct();
    std::size_t lcfs = 0, lcfspr = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched == SchedStrategy::LCFS) lcfs = i + 1;
        if (sn.stations[i].sched == SchedStrategy::LCFSPR) lcfspr = i + 1;
    }
    REQUIRE(lcfs != 0);
    REQUIRE(lcfspr != 0);
    return mva::solver_mva_lcfsqn(sn, lcfs, lcfspr);
}

/** `solver_mva`'s LCFS gate, which is where a model outside the class is refused. */
mva::MvaSolution<double> lcfsqn_gated(qn::Network<double>& m) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    mva::MvaOptions opt;
    opt.method = "exact";
    return mva::solver_mva(sn, mva::sn_get_demands_chain(sn), opt);
}

mva::MvaSolution<double> sum_direct(qn::Network<double>& m, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    return mva::solver_mva_sum(m.get_struct(), opt);
}

}  // namespace

// ---------------------------------------------------------------------------
// the closed LCFS + LCFS-PR pair (solver_mva.h)
// ---------------------------------------------------------------------------

TEST_CASE("mva sum/lcfsqn: the LCFS recursion and the LCFS permanents agree") {
    // Three classes, one job each: the recursion walks the population lattice,
    // the NC route sums permanents over the boundary position. Both are exact
    // for this model class, so any gap is a port defect in one of them.
    qn::Network<double> m =
        lcfs_pair<double>({{2.0, 0.5, 1.0 / 3.0}, {3.0, 1.5, 0.25}}, {1.0, 1.0, 1.0});
    const mva::MvaSolution<double> a = lcfsqn_direct(m);
    const mva::AvgResult<double> b = nc_run(m);
    CHECK(b.actualmethod == "default/lcfsqn.ca");
    for (std::size_t r = 0; r < 3; ++r) {
        CHECK(a.Q(0, r) == doctest::Approx(b.QN(0, r)).epsilon(1e-9));
        CHECK(a.Q(1, r) == doctest::Approx(b.QN(1, r)).epsilon(1e-9));
        CHECK(a.U(0, r) == doctest::Approx(b.UN(0, r)).epsilon(1e-9));
        CHECK(a.U(1, r) == doctest::Approx(b.UN(1, r)).epsilon(1e-9));
        CHECK(a.Tp(0, r) == doctest::Approx(b.TN(0, r)).epsilon(1e-9));
        CHECK(a.R(1, r) == doctest::Approx(b.RN(1, r)).epsilon(1e-9));
        // every job of a class is at one of the two stations
        CHECK(a.Q(0, r) + a.Q(1, r) == doctest::Approx(1.0).epsilon(1e-9));
    }
    CHECK(std::isnan(a.lG));  // this route computes no normalizing constant
}

TEST_CASE("mva sum/lcfsqn: the two LCFS routes agree when a class has two jobs") {
    // The NC route expands a class of N_r > 1 into exchangeable copies against
    // G prod_r N_r!; the MVA recursion never expands anything. The two agreeing
    // here is what distinguishes a correct expansion from a coincidence at one
    // job per class.
    qn::Network<double> m = lcfs_pair<double>({{2.0, 0.5}, {3.0, 1.5}}, {2.0, 1.0});
    const mva::MvaSolution<double> a = lcfsqn_direct(m);
    const mva::AvgResult<double> b = nc_run(m);
    CHECK(a.Q(0, 0) == doctest::Approx(b.QN(0, 0)).epsilon(1e-9));
    CHECK(a.Q(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-9));
    CHECK(a.Q(0, 1) == doctest::Approx(b.QN(0, 1)).epsilon(1e-9));
    CHECK(a.Q(1, 1) == doctest::Approx(b.QN(1, 1)).epsilon(1e-9));
    CHECK(a.Tp(0, 0) == doctest::Approx(b.TN(0, 0)).epsilon(1e-9));
    CHECK(a.Tp(0, 1) == doctest::Approx(b.TN(0, 1)).epsilon(1e-9));
    CHECK(a.Q(0, 0) + a.Q(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(a.Q(0, 1) + a.Q(1, 1) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("mva sum/lcfsqn: a class of zero population does not disturb the others") {
    // REGRESSION. `solver_mva_lcfsqn` left alpha and beta at zero for a class
    // with no jobs, copying solver_mva_lcfsqn.m:31-34, which is safe in MATLAB
    // because pfqn_lcfsqn_mva.m validates nothing. The C++ api rejects a
    // non-positive alpha for EVERY class regardless of its population
    // (pfqn_lcfsqn_mva.h:97-99), so a closed LCFS model carrying an empty class
    // threw InputError where MATLAB answers. It must solve, and the empty class
    // must be inert: the two-class model has to return exactly what the model
    // without that class returns.
    qn::Network<double> two = lcfs_pair<double>({{2.0, 3.0}, {3.0, 4.0}}, {2.0, 0.0});
    qn::Network<double> one = lcfs_pair<double>({{2.0}, {3.0}}, {2.0});
    mva::MvaSolution<double> a;
    REQUIRE_NOTHROW(a = lcfsqn_direct(two));
    const mva::MvaSolution<double> b = lcfsqn_direct(one);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(a.Q(i, 0) == doctest::Approx(b.Q(i, 0)).epsilon(1e-12));
        CHECK(a.U(i, 0) == doctest::Approx(b.U(i, 0)).epsilon(1e-12));
        CHECK(a.R(i, 0) == doctest::Approx(b.R(i, 0)).epsilon(1e-12));
        CHECK(a.Tp(i, 0) == doctest::Approx(b.Tp(i, 0)).epsilon(1e-12));
        // the empty class carries no jobs, no work and no flow anywhere
        CHECK(a.Q(i, 1) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(a.U(i, 1) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(a.Tp(i, 1) == doctest::Approx(0.0).epsilon(1e-12));
    }
    CHECK(a.X[0] == doctest::Approx(b.X[0]).epsilon(1e-12));
    CHECK(a.X[1] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("mva sum/lcfsqn: one circulating job gives the cycle-time closed form") {
    // With a single job there is no queueing anywhere and the discipline stops
    // mattering: the job alternates between the two stations, so the cycle time
    // is alpha + beta and X = 1/(alpha+beta), Q_i = X S_i, U_i = Q_i, R_i = S_i.
    const double mu1 = 2.0, mu2 = 5.0;
    const double alpha = 1.0 / mu1, beta = 1.0 / mu2;
    qn::Network<double> m = lcfs_pair<double>({{mu1}, {mu2}}, {1.0});
    const mva::MvaSolution<double> r = lcfsqn_direct(m);
    const double X = 1.0 / (alpha + beta);
    CHECK(r.X[0] == doctest::Approx(X).epsilon(1e-12));
    CHECK(r.Tp(0, 0) == doctest::Approx(X).epsilon(1e-12));
    CHECK(r.Tp(1, 0) == doctest::Approx(X).epsilon(1e-12));
    CHECK(r.Q(0, 0) == doctest::Approx(X * alpha).epsilon(1e-12));
    CHECK(r.Q(1, 0) == doctest::Approx(X * beta).epsilon(1e-12));
    CHECK(r.U(0, 0) == doctest::Approx(X * alpha).epsilon(1e-12));
    CHECK(r.U(1, 0) == doctest::Approx(X * beta).epsilon(1e-12));
    CHECK(r.R(0, 0) == doctest::Approx(alpha).epsilon(1e-12));
    CHECK(r.R(1, 0) == doctest::Approx(beta).epsilon(1e-12));
    CHECK(r.C[0] == doctest::Approx(alpha + beta).epsilon(1e-12));
}

TEST_CASE("mva sum/lcfsqn: the LCFS recursion is exact-capable arithmetic") {
    // The header claims the recursion is exact and carries no arithmetic gate,
    // so it must run under Rational and land on the double answer. The NC
    // sibling cannot be compared here: it reports lG = log(G) and is refused by
    // name without transcendentals.
    qn::Network<double> md = lcfs_pair<double>({{2.0, 4.0}, {5.0, 8.0}}, {1.0, 1.0});
    qn::Network<Rational> mr = lcfs_pair<Rational>({{2.0, 4.0}, {5.0, 8.0}}, {1.0, 1.0});
    const mva::MvaSolution<double> a = lcfsqn_direct(md);
    const mva::MvaSolution<Rational> b = lcfsqn_direct(mr);
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(num_traits<Rational>::to_double(b.Q(0, r)) ==
              doctest::Approx(a.Q(0, r)).epsilon(1e-12));
        CHECK(num_traits<Rational>::to_double(b.Q(1, r)) ==
              doctest::Approx(a.Q(1, r)).epsilon(1e-12));
        CHECK(num_traits<Rational>::to_double(b.Tp(0, r)) ==
              doctest::Approx(a.Tp(0, r)).epsilon(1e-12));
    }
}

TEST_CASE("mva sum/lcfsqn: the LCFS gate refuses every model outside its class") {
    SUBCASE("no LCFS station at all is not a refusal, it is another model class") {
        qn::Network<double> m = sum_cycle({2.0}, 1.0);
        CHECK_NOTHROW(lcfsqn_gated(m));  // solved by the ordinary product-form path
    }
    SUBCASE("an LCFS station without its LCFS-PR partner") {
        qn::Network<double> m("lcfsalone");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", SchedStrategy::LCFS);
        const std::size_t c = m.add_closed_class("C1", 2.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(lcfsqn_gated(m), UnsupportedError);
    }
    SUBCASE("two LCFS stations, which the two-station closed form has no term for") {
        qn::Network<double> m("lcfstwo");
        const std::size_t q1 = m.add_queue("L1", SchedStrategy::LCFS);
        const std::size_t q2 = m.add_queue("L2", SchedStrategy::LCFS);
        const std::size_t q3 = m.add_queue("P", SchedStrategy::LCFSPR);
        const std::size_t c = m.add_closed_class("C1", 2.0, q1);
        m.set_service(q1, c, D::exp_rate(1.0));
        m.set_service(q2, c, D::exp_rate(2.0));
        m.set_service(q3, c, D::exp_rate(3.0));
        qn::RoutingMatrix<double> P;
        P.set(q1, q2, 1.0);
        P.set(q2, q3, 1.0);
        P.set(q3, q1, 1.0);
        m.link(P);
        CHECK_THROWS_AS(lcfsqn_gated(m), UnsupportedError);
    }
    // A fractional population is NOT refused here: the gate rounds it (llround)
    // rather than rejecting it, which is what the reference does. Asserting a
    // refusal would pin behaviour neither codebase has.
    SUBCASE("an open class, which the closed form has no external world for") {
        qn::Network<double> m("lcfsopen");
        const std::size_t s = m.add_source("Source");
        const std::size_t q1 = m.add_queue("LCFS", SchedStrategy::LCFS);
        const std::size_t q2 = m.add_queue("LCFSPR", SchedStrategy::LCFSPR);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("O1");
        m.set_arrival(s, c, D::exp_rate(0.5));
        m.set_service(q1, c, D::exp_rate(2.0));
        m.set_service(q2, c, D::exp_rate(3.0));
        qn::RoutingMatrix<double> P;
        P.set(s, q1, 1.0);
        P.set(q1, q2, 1.0);
        P.set(q2, k, 1.0);
        m.link(P);
        CHECK_THROWS_AS(lcfsqn_gated(m), UnsupportedError);
    }
}

// ---------------------------------------------------------------------------
// the summation method (solver_mva.h)
// ---------------------------------------------------------------------------

TEST_CASE("mva sum/lcfsqn: the summation method is exact at a population of one") {
    // sum_closed's node function short-circuits at K <= m: a station that can
    // never hold more than its server count queues nothing, so K_i = U_i and the
    // population constraint collapses to X (Z + S) = 1. The ESUM SCV correction
    // is switched off with it, which is why the Erlang service below must not
    // move the answer at all.
    SUBCASE("an Erlang FCFS station, where the correction is switched off") {
        qn::Network<double> m("sumone");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("C1", 1.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::erlang_fit(0.5, 0.5));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        const mva::MvaSolution<double> r = sum_direct(m, "sum");
        const double X = 1.0 / (1.0 + 0.5);
        CHECK(r.X[0] == doctest::Approx(X).epsilon(1e-6));
        CHECK(r.Q(0, 0) == doctest::Approx(X * 1.0).epsilon(1e-6));
        CHECK(r.Q(1, 0) == doctest::Approx(X * 0.5).epsilon(1e-6));
        CHECK(r.U(1, 0) == doctest::Approx(X * 0.5).epsilon(1e-6));
        CHECK(r.Q(0, 0) + r.Q(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
    }
    SUBCASE("an exponential station, where exact MVA is available to agree") {
        qn::Network<double> m("sumoneexp");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("C1", 1.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        const mva::MvaSolution<double> s = sum_direct(m, "sum");
        const mva::AvgResult<double> e = mva_run(m, "exact");
        CHECK(s.Q(0, 0) == doctest::Approx(e.QN(0, 0)).epsilon(1e-6));
        CHECK(s.Q(1, 0) == doctest::Approx(e.QN(1, 0)).epsilon(1e-6));
        CHECK(s.X[0] == doctest::Approx(e.XN[0]).epsilon(1e-6));
    }
}

TEST_CASE("mva sum/lcfsqn: the summation method satisfies the laws it does not assume") {
    // A two-class closed model over a delay, a two-server FCFS station (the
    // Erlang-C branch of the node function) and a PS station. SUM's error here
    // is unbounded, so nothing is compared against another analyzer; what must
    // hold whatever the error is the population constraint the method roots, and
    // the two laws the deaggregation reconstructs the class-level metrics with.
    // Those are properties of the sn-to-api glue, which is what breaks if a row
    // or a chain is misindexed.
    qn::Network<double> m = sum_cycle({4.0, 2.0}, 0.5);
    const mva::MvaSolution<double> r = sum_direct(m, "esum");
    const double N[2] = {4.0, 2.0};
    const double Sd[2] = {1.0, 0.5}, S1[2] = {0.8, 0.4}, S2[2] = {1.0 / 3.0, 0.25};
    for (std::size_t j = 0; j < 2; ++j) {
        double n = 0.0;
        for (std::size_t i = 0; i < 3; ++i) n += r.Q(i, j);
        // the constraint sum_i K_i = N is what the bisection roots, so it holds
        // to the method's own tolerance and not to machine precision
        CHECK(n == doctest::Approx(N[j]).epsilon(1e-5));
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(r.Q(i, j) == doctest::Approx(r.R(i, j) * r.Tp(i, j)).epsilon(1e-9));
        // U = T S / m at the finite-server stations; the delay reports T S,
        // since an infinite server has no per-server share to divide by
        CHECK(r.U(1, j) == doctest::Approx(r.Tp(1, j) * S1[j] / 2.0).epsilon(1e-9));
        CHECK(r.U(2, j) == doctest::Approx(r.Tp(2, j) * S2[j]).epsilon(1e-9));
        CHECK(r.U(0, j) == doctest::Approx(r.Tp(0, j) * Sd[j]).epsilon(1e-9));
    }
}

TEST_CASE("mva sum/lcfsqn: the closing method approaches the open product form") {
    // One open class through two exponential single servers, where the product
    // form gives rho/(1-rho) = 1/3 and 1/2. The closing method replaces the
    // external world by a station of population Kclosed = 5000, so it approaches
    // those from below; the gap is the truncation, and the tolerance says so.
    qn::Network<double> m("sumopen");
    const std::size_t s = m.add_source("Source");
    const std::size_t qa = m.add_queue("Qa", SchedStrategy::FCFS);
    const std::size_t qb = m.add_queue("Qb", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, D::exp_rate(0.5));
    m.set_service(qa, o, D::exp_rate(2.0));
    m.set_service(qb, o, D::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(s, qa, 1.0);
    P.set(qa, qb, 1.0);
    P.set(qb, k, 1.0);
    m.link(P);
    const mva::MvaSolution<double> r = sum_direct(m, "sum");
    CHECK(r.Q(1, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-3));
    CHECK(r.Q(2, 0) == doctest::Approx(0.5).epsilon(1e-3));
    CHECK(r.Q(1, 0) < 1.0 / 3.0);
    CHECK(r.Q(2, 0) < 0.5);
}

TEST_CASE("mva sum/lcfsqn: the summation method refuses a discipline it cannot model") {
    // SUM has a node function for the insensitive disciplines and for FCFS/SIRO
    // and for nothing else. Handing a DPS station one of those functions would
    // return a number for a different station.
    qn::Network<double> m("sumdps");
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
    CHECK_THROWS_AS(sum_direct(m, "sum"), UnsupportedError);
    CHECK_THROWS_AS(sum_direct(m, "esum"), UnsupportedError);
}

TEST_CASE("mva sum/lcfsqn: the summation method refuses exact arithmetic by name") {
    // The bisection stops on a tolerance and the Erlang-C branch is evaluated on
    // the way, so there is no exact-field version of the answer; the analyzer
    // says so rather than returning the iterate its stopping rule happened to
    // select. Contrast the LCFS recursion above, which is exact and runs.
    qn::Network<Rational> m("sumrat");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(c, c, d, q, num_traits<Rational>::from_int(1));
    P.set(c, c, q, d, num_traits<Rational>::from_int(1));
    m.link(P);
    mva::MvaOptions opt;
    opt.method = "sum";
    CHECK_THROWS_AS(mva::solver_mva_sum(m.get_struct(), opt), UnsupportedError);
}
