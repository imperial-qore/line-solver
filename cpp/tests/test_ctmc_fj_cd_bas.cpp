/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The three constructs `ctmc_check_support` used to refuse: fork-join, class and
 * joint dependence, and true-BAS blocking.
 *
 * The oracles are chosen so that none of them is "whatever the code prints".
 * Fork-join is checked against the identity that makes the tag transform correct
 * -- a fork whose branches are the SAME station is a network in which the fork is
 * pure bookkeeping, so its throughput must equal the un-forked model's scaled by
 * the branch count -- and against a two-branch closed model whose join-free
 * counterpart SolverMVA solves exactly. Class dependence is checked against the
 * exponential-server identity beta(n) = min(n, c), which turns a single-server PS
 * station into an exact c-server one, so SolverMVA on the same model with c
 * servers must agree. BAS is checked against the conservation laws it must obey:
 * the population is conserved, and the blocked job is reported at its destination.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/fj_tag.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> Fork -> {two Queues} -> Join -> back to the Delay, one closed class. */
qn::Network<double> fj_closed(double zrate, double mu1, double mu2, double njobs) {
    qn::Network<double> m("fj");
    const std::size_t d = m.add_delay("Think");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("J", f);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(zrate));
    m.set_service(q1, c, Dist::exp_rate(mu1));
    m.set_service(q2, c, Dist::exp_rate(mu2));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Delay -> PS Queue -> Delay, four closed jobs, with BOTH handles on the queue:
 * beta(n) = min(n,2) declaring peak 2, and eta(n) = 3 flat declaring peak 3. Every
 * constant is an integer or a dyadic rational so the model is representable
 * without rounding in any arithmetic.
 */
template <class T>
qn::Network<T> cdjd_model() {
    typedef line::num_traits<T> nt;
    typedef line::lang::Distrib<T> D;
    qn::Network<T> m("cdjd");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, k, D::exp_rate(nt::from_int(1)));
    m.set_service(q, k, D::exp_rate(nt::from_double(1.5)));
    m.set_class_dependence(q, [](const std::vector<T>& n) {
        return std::vector<T>(1, std::min(n[0], line::num_traits<T>::from_int(2)));
    }, std::vector<T>(1, nt::from_int(2)));
    m.set_joint_dependence(q, [](const std::vector<T>&) {
        return std::vector<T>(1, line::num_traits<T>::from_int(3));
    }, std::vector<T>(1, nt::from_int(3)));
    qn::RoutingMatrix<T> P;
    P.set(d, q, nt::from_int(1));
    P.set(q, d, nt::from_int(1));
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fj_tag mints one auxiliary class per branch and tag") {
    qn::Network<double> m = fj_closed(1.0, 2.0, 3.0, 2);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const qn::FjTagged<double> t = qn::fj_tag(sn);

    // 2 branches x 2 tags = 4 auxiliary classes on top of the original one.
    CHECK(t.korig == 1);
    CHECK(t.V.classes.size() == 5);
    CHECK(t.fjsync.size() == 2);  // one firing entry per tag
    for (std::size_t a = 1; a < 5; ++a) {
        CHECK(t.fjclassmap[a] == 1);
        CHECK(t.V.classes[a].population == 0.0);
    }
    CHECK(t.fjclassmap[0] == 0);
    // Both branches appear in every firing entry, and the tags are distinct.
    CHECK(t.fjsync[0].branchheads.size() == 2);
    CHECK(t.fjsync[0].tag == 1);
    CHECK(t.fjsync[1].tag == 2);
    CHECK(t.fjsync[0].auxclasses != t.fjsync[1].auxclasses);
    // The Fork became stateful, which is what gives it a state row to hold the
    // parent job in.
    CHECK(t.V.stateful_index(t.fjsync[0].fork) != 0);
}

TEST_CASE("a degenerate one-branch fork is refused by name") {
    qn::Network<double> m("fj1");
    const std::size_t d = m.add_delay("Think");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("J", f);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(q1, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    CHECK_THROWS(qn::fj_tag(m.get_struct()));
}

TEST_CASE("one job in a fork-join cycle matches the closed form exactly") {
    // WHY N = 1 IS THE ORACLE. With a single parent there is no queueing anywhere,
    // so every quantity is an expectation of independent exponentials:
    //   cycle time    Z + E[max(S1, S2)],  E[max] = 1/mu1 + 1/mu2 - 1/(mu1+mu2)
    //   Delay QLen    Z / cycle                       (Little on the think time)
    //   branch+join   2*E[max] / cycle
    // The last one is the test that matters: each sibling EXISTS from the firing to
    // the join, i.e. for E[max] whether it is in service or waiting at the join, so
    // the sibling-time per cycle is exactly 2*E[max]. A join that released a parent
    // early would shorten it; one that never released would make it unbounded.
    const double mu1 = 2.0, mu2 = 3.0, Z = 1.0;
    qn::Network<double> m = fj_closed(1.0 / Z, mu1, mu2, 1);
    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);

    const double emax = 1.0 / mu1 + 1.0 / mu2 - 1.0 / (mu1 + mu2);
    const double cycle = Z + emax;
    CHECK(r.CN[0] == doctest::Approx(cycle).epsilon(1e-6));
    CHECK(r.XN[0] == doctest::Approx(1.0 / cycle).epsilon(1e-6));
    CHECK(r.QN(0, 0) == doctest::Approx(Z / cycle).epsilon(1e-6));
    double sib = 0;
    for (std::size_t i = 1; i < r.QN.rows(); ++i) sib += r.QN(i, 0);
    CHECK(sib == doctest::Approx(2.0 * emax / cycle).epsilon(1e-6));
    // The branch services themselves are 1/mu each, per cycle.
    CHECK(r.QN(1, 0) == doctest::Approx(1.0 / mu1 / cycle).epsilon(1e-6));
    CHECK(r.QN(2, 0) == doctest::Approx(1.0 / mu2 / cycle).epsilon(1e-6));
}

TEST_CASE("the parent throughput is the same all the way round the cycle") {
    const double N = 2;
    qn::Network<double> m = fj_closed(1.0, 2.0, 3.0, N);
    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);

    CHECK(r.TN(0, 0) > 0.0);
    // A closed cycle is traversed once per parent completion, so the Delay, both
    // branches and the Join all carry the SAME parent throughput. The branch
    // columns are the FOLDED sibling throughputs, which is why this is a real
    // check on the foldback and not a tautology.
    for (std::size_t i = 1; i < r.TN.rows(); ++i)
        CHECK(r.TN(i, 0) == doctest::Approx(r.TN(0, 0)).epsilon(1e-6));

    // Little's law on the PARENT class closes the loop: N parents circulate, so
    // N = X * (cycle time), and the cycle time is the Delay plus the sibling
    // residence divided by the branch count -- a forked cycle traverses both
    // branches CONCURRENTLY, so the summed queue lengths exceed N and the total is
    // NOT a population. That is why no conservation law is asserted here.
    CHECK(r.CN[0] == doctest::Approx(N / r.XN[0]).epsilon(1e-9));
}

TEST_CASE("class dependence beta(n)=min(n,c) reproduces a c-server station") {
    // A single-server PS station whose rate scales with its own population up to
    // c is EXACTLY a c-server station for exponential service, so SolverMVA on
    // the multiserver model is an independent oracle for the scaled CTMC.
    const double c = 2.0, N = 4;
    qn::Network<double> cd("cd");
    const std::size_t d = cd.add_delay("Think");
    const std::size_t q = cd.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = cd.add_closed_class("C1", N, d);
    cd.set_service(d, k, Dist::exp_rate(1.0));
    cd.set_service(q, k, Dist::exp_rate(1.5));
    cd.set_class_dependence(q, [c](const std::vector<double>& n) {
        return std::vector<double>(1, std::min(n[0], c));
    }, std::vector<double>(1, c));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    cd.link(P);

    qn::Network<double> ms("ms");
    const std::size_t d2 = ms.add_delay("Think");
    const std::size_t q2 = ms.add_queue("Q", SchedStrategy::PS);
    const std::size_t k2 = ms.add_closed_class("C1", N, d2);
    ms.set_service(d2, k2, Dist::exp_rate(1.0));
    ms.set_service(q2, k2, Dist::exp_rate(1.5));
    ms.set_number_of_servers(q2, c);
    qn::RoutingMatrix<double> P2;
    P2.set(d2, q2, 1.0);
    P2.set(q2, d2, 1.0);
    ms.link(P2);

    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> a = ctmc::solver_ctmc_run_analyzer(cd.get_struct(), opt);
    const line::mva::AvgResult<double> b = ctmc::solver_ctmc_run_analyzer(ms.get_struct(), opt);
    CHECK(a.QN(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-9));
    CHECK(a.TN(1, 0) == doctest::Approx(b.TN(1, 0)).epsilon(1e-9));
    // Utilization at a class-dependent station is T/mu/peak, which is the same
    // T*S/c the multiserver station reports.
    CHECK(a.UN(1, 0) == doctest::Approx(b.UN(1, 0)).epsilon(1e-9));
}

TEST_CASE("class dependence without a declared peak is a model defect") {
    qn::Network<double> cd("cdnopeak");
    const std::size_t d = cd.add_delay("Think");
    const std::size_t q = cd.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = cd.add_closed_class("C1", 2, d);
    cd.set_service(d, k, Dist::exp_rate(1.0));
    cd.set_service(q, k, Dist::exp_rate(1.5));
    cd.set_class_dependence(q, [](const std::vector<double>& n) {
        return std::vector<double>(1, std::min(n[0], 2.0));
    });
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    cd.link(P);
    ctmc::CtmcOptions opt;
    CHECK_THROWS(ctmc::solver_ctmc_run_analyzer(cd.get_struct(), opt));
}

TEST_CASE("joint dependence scales the generator and is reported on its own peak") {
    // eta(n) = min(n_1, 2) reads the class-1 marginal only, so on a SINGLE-class
    // model it is numerically the same scaling as the class-dependent case above
    // and must give the same chain -- what differs is only which declared peak
    // the utilization column is normalized by, and here they coincide too.
    const double N = 4, c = 2.0;
    qn::Network<double> jd("jd");
    const std::size_t d = jd.add_delay("Think");
    const std::size_t q = jd.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = jd.add_closed_class("C1", N, d);
    jd.set_service(d, k, Dist::exp_rate(1.0));
    jd.set_service(q, k, Dist::exp_rate(1.5));
    jd.set_joint_dependence(q, [c](const std::vector<double>& n) {
        return std::vector<double>(1, std::min(n[0], c));
    }, std::vector<double>(1, c));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    jd.link(P);

    qn::Network<double> ms("ms2");
    const std::size_t d2 = ms.add_delay("Think");
    const std::size_t q2 = ms.add_queue("Q", SchedStrategy::PS);
    const std::size_t k2 = ms.add_closed_class("C1", N, d2);
    ms.set_service(d2, k2, Dist::exp_rate(1.0));
    ms.set_service(q2, k2, Dist::exp_rate(1.5));
    ms.set_number_of_servers(q2, c);
    qn::RoutingMatrix<double> P2;
    P2.set(d2, q2, 1.0);
    P2.set(q2, d2, 1.0);
    ms.link(P2);

    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> a = ctmc::solver_ctmc_run_analyzer(jd.get_struct(), opt);
    const line::mva::AvgResult<double> b = ctmc::solver_ctmc_run_analyzer(ms.get_struct(), opt);
    CHECK(a.QN(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-9));
    CHECK(a.TN(1, 0) == doctest::Approx(b.TN(1, 0)).epsilon(1e-9));
}

TEST_CASE("the class and joint declared peaks fold multiplicatively, not by overwrite") {
    // Both handles at ONE station: beta(n) = min(n,2) with peak 2, eta(n) = 3
    // flat with peak 3. The two scale the same nominal rate, so the generator
    // carries their product and the peak attainable rate is mu*2*3 = 9.
    const double mu = 1.5;
    qn::Network<double> m = cdjd_model<double>();

    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);
    // The reported utilization must be T/(mu*cdpeak*jdpeak); an overwrite would
    // drop the cd factor and report exactly twice this.
    CHECK(r.UN(1, 0) == doctest::Approx(r.TN(1, 0) / (mu * 2.0 * 3.0)).epsilon(1e-9));
    CHECK(r.UN(1, 0) < 1.0);
    // The chain itself carries the product: a station serving at 3*min(n,2)*mu is
    // exactly a 2-server station whose per-server rate is 3*mu.
    qn::Network<double> ms("cdjdms");
    const std::size_t d2 = ms.add_delay("Think");
    const std::size_t q2 = ms.add_queue("Q", SchedStrategy::PS);
    const std::size_t k2 = ms.add_closed_class("C1", 4.0, d2);
    ms.set_service(d2, k2, Dist::exp_rate(1.0));
    ms.set_service(q2, k2, Dist::exp_rate(3.0 * mu));
    ms.set_number_of_servers(q2, 2.0);
    qn::RoutingMatrix<double> P2;
    P2.set(d2, q2, 1.0);
    P2.set(q2, d2, 1.0);
    ms.link(P2);
    const line::mva::AvgResult<double> b = ctmc::solver_ctmc_run_analyzer(ms.get_struct(), opt);
    CHECK(r.QN(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(b.TN(1, 0)).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(b.UN(1, 0)).epsilon(1e-9));
    // Python-native SolverCTMC on the same model, which fixes the absolute scale
    // rather than only the internal consistency the checks above pin.
    CHECK(r.QN(1, 0) == doctest::Approx(0.756428426807).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(3.24357157319).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.360396841466).epsilon(1e-9));
}

TEST_CASE("the cd/jd peak product survives an exact instantiation") {
    // The golden above is a reference NUMBER, so on its own it cannot separate a
    // correct answer from one a rounding path happens to reproduce. Rerunning the
    // same model over a field with no rounding does: the declared peaks are
    // integers and the rates are rationals, so `Rational` computes the fixed
    // point exactly and any drift would be arithmetic, not modelling. This case
    // lives here rather than in test_cross_arith.cpp because what it pins is the
    // cd/jd contract, not the linear-algebra path underneath.
    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<line::Rational> e =
        ctmc::solver_ctmc_run_analyzer(cdjd_model<line::Rational>().get_struct(), opt);
    const line::mva::AvgResult<line::Real50> h =
        ctmc::solver_ctmc_run_analyzer(cdjd_model<line::Real50>().get_struct(), opt);
    typedef line::num_traits<line::Rational> ent;
    typedef line::num_traits<line::Real50> hnt;
    CHECK(ent::to_double(e.QN(1, 0)) == doctest::Approx(0.756428426807).epsilon(1e-11));
    CHECK(ent::to_double(e.TN(1, 0)) == doctest::Approx(3.24357157319).epsilon(1e-11));
    CHECK(ent::to_double(e.UN(1, 0)) == doctest::Approx(0.360396841466).epsilon(1e-11));
    CHECK(hnt::to_double(h.QN(1, 0)) == doctest::Approx(ent::to_double(e.QN(1, 0))).epsilon(1e-14));
    CHECK(hnt::to_double(h.UN(1, 0)) == doctest::Approx(ent::to_double(e.UN(1, 0))).epsilon(1e-14));
    // Util = T/(mu*cdpeak*jdpeak) holds with NO tolerance at all in exact arithmetic
    CHECK(e.UN(1, 0) == e.TN(1, 0) / (ent::from_double(1.5) * ent::from_int(6)));
}

TEST_CASE("true BAS reserves a blocked marker and conserves the population") {
    // Two FCFS queues in a closed cycle; the SECOND has a finite capacity and
    // declares BAS, so the first must hold a completed job when the second is
    // full. Under the destination declaration form the marker still belongs to
    // the upstream station, which is what `refresh_bas_blocking` must find.
    const double N = 3;
    qn::Network<double> m("bas");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", N, q1);
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(1.0));
    m.set_capacity(q2, 2);
    m.set_drop_rule(q2, c, line::lang::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.isbasblocking.size() >= 2);
    CHECK(sn.isbasblocking[q1 - 1]);        // the upstream station holds the job
    CHECK(!sn.isbasblocking[q2 - 1]);       // the full destination does not
    CHECK(sn.isbasdestination[1][0]);       // refusing at Q2 blocks Q1
    CHECK(sn.nvars[q1 - 1][2 * sn.nclasses] == 1);

    ctmc::CtmcOptions opt;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, opt);
    double tot = 0;
    for (std::size_t i = 0; i < r.QN.rows(); ++i) tot += r.QN(i, 0);
    CHECK(tot == doctest::Approx(N).epsilon(1e-6));
    // The blocked job is reported AT ITS DESTINATION, so Q2's mean queue length
    // legitimately exceeds its declared buffer of 2: the held job is physically at
    // Q1's server and counted at Q2 by the reference's convention. What it may
    // never exceed is the population.
    CHECK(r.QN(1, 0) > 2.0);
    CHECK(r.QN(1, 0) < N);
    CHECK(r.TN(0, 0) == doctest::Approx(r.TN(1, 0)).epsilon(1e-6));
}

TEST_CASE("an unbounded destination reserves no marker") {
    // BAS is declared, but nothing downstream can ever fill, so the marker would
    // widen every state for a transition that can never fire.
    qn::Network<double> m("basinf");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2, q1);
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(1.0));
    m.set_drop_rule(q1, c, line::lang::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    for (std::size_t i = 0; i < sn.isbasblocking.size(); ++i) CHECK(!sn.isbasblocking[i]);
}
