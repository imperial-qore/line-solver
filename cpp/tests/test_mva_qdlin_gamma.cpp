/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `amva.qdlin` takes the CLASS-AGGREGATE Linearizer correction, not the
 * per-class one, and must therefore not be an alias of `amva.lin`.
 *
 * `solver_amvald.m` allocates the (K,M,K) per-class gamma array for `qdlin` but
 * writes `gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) - sum_r Q(k,r)/Nt` into it with
 * two subscripts, which MATLAB linear-indexes to `(s,k,1)`; slices `2..K` stay
 * zero while every reader indexes gamma per class. MATLAB, the JAR and native
 * Python all do this. THIS PORT DID NOT until 2026-09-04: it filled gamma per
 * class for `lin` and `qdlin` alike, with no other branch separating the two,
 * so C++ `qdlin` was a bit-exact alias of C++ `lin` and disagreed with the
 * other three codebases on every multichain model. On the two-chain model
 * below it returned `X = 0.838443047588` for C1 -- the `lin` figure -- where
 * the other three return `0.829094887846`.
 *
 * Every reference value here is native-Python `SolverMVA(model, 'qdlin')` on
 * the identical model, and the port now reproduces it to twelve decimals. The
 * tolerance is the fixed point's own `iter_tol`, 1e-6: the last digits of an
 * iterate stopped on a tolerance are an artifact of the stopping test, so a
 * reference is only meaningful alongside the tolerance that produced it.
 *
 * The K = 1 case is the other half of the property and is asserted too. With
 * one chain the aggregate correction IS the per-class one, so `lin` and `qdlin`
 * must stay identical there; a fix that made them differ everywhere would be as
 * wrong as the alias was.
 */
#include <cmath>
#include <cstddef>
#include <string>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

// Station rows, in the order the builder creates them.
constexpr std::size_t kThink = 0, kQ1 = 1, kQ2 = 2, kQ3 = 3;
constexpr std::size_t kC1 = 0, kC2 = 1;

const double kTol = 1e-6;  // the AMVA fixed point's iter_tol

/** Think -> Q1 -> Q2 -> Q3 -> Think, two closed chains. */
qn::Network<double> two_chain() {
    qn::Network<double> m("qdlin2");
    const std::size_t dly = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 4.0, dly);
    const std::size_t c2 = m.add_closed_class("C2", 3.0, dly);
    m.set_service(dly, c1, D::exp_rate(1.0 / 2.0));
    m.set_service(dly, c2, D::exp_rate(1.0 / 1.0));
    m.set_service(q1, c1, D::exp_rate(1.0 / 0.4));
    m.set_service(q2, c1, D::exp_rate(1.0 / 0.6));
    m.set_service(q3, c1, D::exp_rate(1.0 / 0.2));
    m.set_service(q1, c2, D::exp_rate(1.0 / 0.9));
    m.set_service(q2, c2, D::exp_rate(1.0 / 0.1));
    m.set_service(q3, c2, D::exp_rate(1.0 / 0.5));
    qn::RoutingMatrix<double> P;
    const std::size_t cls[2] = {c1, c2};
    for (std::size_t k = 0; k < 2; ++k) {
        P.set(cls[k], cls[k], dly, q1, 1.0);
        P.set(cls[k], cls[k], q1, q2, 1.0);
        P.set(cls[k], cls[k], q2, q3, 1.0);
        P.set(cls[k], cls[k], q3, dly, 1.0);
    }
    m.link(P);
    return m;
}

/** Think -> Q1 -> Q2 -> Think, one closed chain, Q1 a two-server station.
 *
 * The second server is what keeps `lin` on solver_amvald: a single-server
 * product-form model sends the lin family to pfqn_linearizermx instead, and the
 * two kernels are not the same computation, so the comparison would say nothing
 * about the gamma fill.
 */
qn::Network<double> one_chain_multiserver() {
    qn::Network<double> m("qdlin1");
    const std::size_t dly = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, dly);
    m.set_number_of_servers(q1, 2.0);
    m.set_service(dly, c1, D::exp_rate(1.0 / 1.0));
    m.set_service(q1, c1, D::exp_rate(1.0 / 0.5));
    m.set_service(q2, c1, D::exp_rate(1.0 / 0.3));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, dly, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, dly, 1.0);
    m.link(P);
    return m;
}

mva::MvaSolution<double> solve(qn::Network<double>& m, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_analyzer(m.get_struct(), opt, init);
}

}  // namespace

TEST_CASE("amva.qdlin reproduces the class-aggregate correction of the reference") {
    qn::Network<double> m = two_chain();
    const mva::MvaSolution<double> r = solve(m, "amva.qdlin");
    CHECK(r.method == "qdlin");

    CHECK(r.X[kC1] == doctest::Approx(0.829094887846).epsilon(kTol));
    CHECK(r.X[kC2] == doctest::Approx(0.636814002325).epsilon(kTol));

    CHECK(r.Q(kQ1, kC1) == doctest::Approx(1.202347197630).epsilon(kTol));
    CHECK(r.Q(kQ2, kC1) == doctest::Approx(0.845540644436).epsilon(kTol));
    CHECK(r.Q(kQ3, kC1) == doctest::Approx(0.293922388308).epsilon(kTol));
    CHECK(r.Q(kQ1, kC2) == doctest::Approx(1.713825783471).epsilon(kTol));
    CHECK(r.Q(kQ2, kC2) == doctest::Approx(0.132121612353).epsilon(kTol));
    CHECK(r.Q(kQ3, kC2) == doctest::Approx(0.517238825234).epsilon(kTol));

    CHECK(r.U(kQ1, kC1) == doctest::Approx(0.331637955139).epsilon(kTol));
    CHECK(r.U(kQ2, kC2) == doctest::Approx(0.063681400233).epsilon(kTol));

    // The delay is charged the think time at the chain's throughput.
    CHECK(r.Q(kThink, kC1) == doctest::Approx(1.658189779988).epsilon(kTol));
    CHECK(r.R(kThink, kC2) == doctest::Approx(1.0).epsilon(kTol));

    // Little's law over the whole network, per chain.
    const double n[2] = {4.0, 3.0};
    for (std::size_t c = 0; c < 2; ++c) {
        double q = 0.0;
        for (std::size_t i = 0; i < 4; ++i) q += r.Q(i, c);
        CHECK(q == doctest::Approx(n[c]).epsilon(1e-6));
    }
}

TEST_CASE("amva.qdlin is NOT an alias of amva.lin on a multichain model") {
    qn::Network<double> mq = two_chain();
    qn::Network<double> ml = two_chain();
    const mva::MvaSolution<double> q = solve(mq, "amva.qdlin");
    const mva::MvaSolution<double> l = solve(ml, "amva.lin");
    CHECK(l.method == "lin");

    // The per-class fill this port used to apply to both is exactly the lin
    // figure, so an alias would land here.
    CHECK(l.X[kC1] == doctest::Approx(0.838455752174).epsilon(kTol));
    CHECK(l.Q(kQ1, kC1) == doctest::Approx(1.159794545847).epsilon(kTol));

    CHECK(std::fabs(q.X[kC1] - l.X[kC1]) > 1e-3);
    CHECK(std::fabs(q.Q(kQ1, kC1) - l.Q(kQ1, kC1)) > 1e-2);
}

TEST_CASE("with one chain the two corrections coincide and the methods agree") {
    qn::Network<double> mq = one_chain_multiserver();
    qn::Network<double> ml = one_chain_multiserver();
    const mva::MvaSolution<double> q = solve(mq, "amva.qdlin");
    const mva::MvaSolution<double> l = solve(ml, "amva.lin");

    CHECK(q.Q(kQ1, kC1) == doctest::Approx(0.784408538911).epsilon(kTol));
    CHECK(q.Q(kQ2, kC1) == doctest::Approx(0.646725716767).epsilon(kTol));
    // Identical, not merely close: with K = 1 the aggregate correction and the
    // per-class one are the same number, so the two methods run the same
    // arithmetic on the same iterates.
    for (std::size_t i = 0; i < 3; ++i) CHECK(q.Q(i, kC1) == doctest::Approx(l.Q(i, kC1)).epsilon(1e-14));
    CHECK(q.X[kC1] == doctest::Approx(l.X[kC1]).epsilon(1e-14));
}
