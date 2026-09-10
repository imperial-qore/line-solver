/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The one MIXED model solver_amva keeps in its product-form branch.
 *
 * `cond3` in `solver_amva.m:146` is not "no open classes": a STRICT product-form
 * network solved by `lin` stays with the product-form kernels, where
 * `pfqn_linearizermx` takes the open classes in closed form (`X_r = lambda_r`,
 * `U^o(i) = sum_r lambda_r L(i,r)`) and inflates the closed demands by
 * `1/(1-U^o(i))`. Every other open model goes to `solver_amvald`, the only path
 * carrying the open-class corrections.
 *
 * The port omitted that clause until 2026-09-04 and sent every open model to
 * `solver_amvald`, so an explicit `amva.lin` ran a DIFFERENT algorithm here than
 * in MATLAB and the JAR -- on this model X of the closed class came back
 * 0.83146863611949 after 373 iterations instead of the 0.84186985221022 below
 * after 63. The `pfqn_linearizermx` call site also passed a zero arrival-rate
 * vector, which would have solved the closed subnetwork with no open
 * interference at all had it ever been reached.
 *
 * The model is the network form of the mixed case of
 * `cpp/tests/test_pfqn_linearizer.cpp`: two PS queues, class O1 open at
 * lambda = 0.4 with demands 0.5 and 0.4, class C1 closed at N = 3 with demands
 * 0.3 and 0.6 and think time 2. So the demand matrix reaching the kernel is that
 * test's `demands()` exactly, and every value below is the MATLAB reference it
 * already asserts, reached through the SOLVER rather than through a direct
 * kernel call. The tolerance is pinned to the kernel test's 1e-8 for the same
 * reason it is there: the last digits of a fixed point stopped on a tolerance
 * are an artifact of the stopping test, so the reference is only meaningful
 * alongside the tolerance that produced it.
 */
#include <cmath>
#include <string>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

// Station rows, in the order the builder creates them; a Sink is not a station.
constexpr std::size_t kSrc = 0, kQ1 = 1, kQ2 = 2, kThink = 3;
// Chain columns: one class per chain, the open class declared first.
constexpr std::size_t kOpen = 0, kClosed = 1;

const double kTol = 1e-8;  // the tolerance the Linearizer fixed point stops on

/** Src -> Q1 -> Q2 -> Sink for O1, Think -> Q1 -> Q2 -> Think for C1. */
qn::Network<double> mixed_model() {
    qn::Network<double> m("mixlin");
    const std::size_t src = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t dly = m.add_delay("Think");
    const std::size_t snk = m.add_sink("Snk");
    const std::size_t oc = m.add_open_class("O1");
    const std::size_t cc = m.add_closed_class("C1", 3.0, dly);
    m.set_arrival(src, oc, D::exp_rate(0.4));
    m.set_service(q1, oc, D::exp_rate(1.0 / 0.5));
    m.set_service(q2, oc, D::exp_rate(1.0 / 0.4));
    m.set_service(q1, cc, D::exp_rate(1.0 / 0.3));
    m.set_service(q2, cc, D::exp_rate(1.0 / 0.6));
    m.set_service(dly, cc, D::exp_rate(1.0 / 2.0));
    qn::RoutingMatrix<double> P;
    P.set(oc, oc, src, q1, 1.0);
    P.set(oc, oc, q1, q2, 1.0);
    P.set(oc, oc, q2, snk, 1.0);
    P.set(cc, cc, dly, q1, 1.0);
    P.set(cc, cc, q1, q2, 1.0);
    P.set(cc, cc, q2, dly, 1.0);
    m.link(P);
    return m;
}

mva::MvaSolution<double> solve_lin(qn::Network<double>& m) {
    mva::MvaOptions opt;
    opt.method = "amva.lin";
    opt.tol = kTol;
    opt.iter_max = 1000;
    Matrix<double> init;
    return mva::solver_mva_analyzer(m.get_struct(), opt, init);
}

}  // namespace

TEST_CASE("solver_amva keeps a mixed product-form model on the 'lin' kernel") {
    qn::Network<double> m = mixed_model();
    const mva::MvaSolution<double> r = solve_lin(m);
    CHECK(r.method == "lin");

    // An open chain's throughput is its arrival rate, and the kernel returns it
    // unchanged; the closed chain is the Linearizer fixed point.
    CHECK(r.X[kOpen] == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(r.X[kClosed] == doctest::Approx(0.84186985221022).epsilon(kTol));

    CHECK(r.Q(kQ1, kOpen) == doctest::Approx(0.34914456045114).epsilon(kTol));
    CHECK(r.Q(kQ1, kClosed) == doctest::Approx(0.39657824180458).epsilon(kTol));
    CHECK(r.Q(kQ2, kOpen) == doctest::Approx(0.36565372452857).epsilon(kTol));
    CHECK(r.Q(kQ2, kClosed) == doctest::Approx(0.91968205377498).epsilon(kTol));

    // The utilization law holds on the ORIGINAL demands, open classes included:
    // 0.4*0.5 and 0.4*0.4 are the arrival rate times the demand, undiscounted.
    CHECK(r.U(kQ1, kOpen) == doctest::Approx(0.2).epsilon(1e-12));
    CHECK(r.U(kQ2, kOpen) == doctest::Approx(0.16).epsilon(1e-12));
    CHECK(r.U(kQ1, kClosed) == doctest::Approx(0.25256095566307).epsilon(kTol));
    CHECK(r.U(kQ2, kClosed) == doctest::Approx(0.50512191132613).epsilon(kTol));

    CHECK(r.R(kQ1, kOpen) == doctest::Approx(0.87286140112786).epsilon(kTol));
    CHECK(r.R(kQ1, kClosed) == doctest::Approx(0.47106834953575).epsilon(kTol));
    CHECK(r.R(kQ2, kOpen) == doctest::Approx(0.91413431132142).epsilon(kTol));
    CHECK(r.R(kQ2, kClosed) == doctest::Approx(1.0924278276036).epsilon(kTol));
}

TEST_CASE("solver_amva reconstructs the delay and Source rows of a mixed model") {
    qn::Network<double> m = mixed_model();
    const mva::MvaSolution<double> r = solve_lin(m);

    // The delay is not a row of the kernel's demand matrix: solver_amva charges it
    // the ORIGINAL think time at the chain's throughput, Q = X Z, and reports the
    // same figure as its utilization.
    CHECK(r.Q(kThink, kClosed) == doctest::Approx(1.6837397044204).epsilon(kTol));
    CHECK(r.U(kThink, kClosed) == doctest::Approx(1.6837397044204).epsilon(kTol));
    CHECK(r.R(kThink, kClosed) == doctest::Approx(2.0).epsilon(kTol));
    CHECK(r.Q(kThink, kOpen) == doctest::Approx(0.0));

    // A Source is a station but not a queueing one: it carries the open chain's
    // throughput and nothing else, so charging interarrival time as residence
    // would inflate the response time by the whole open cycle.
    CHECK(r.Tp(kSrc, kOpen) == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(r.Q(kSrc, kOpen) == doctest::Approx(0.0));
    CHECK(r.R(kSrc, kOpen) == doctest::Approx(0.0));

    // The closed population is conserved across the two queues and the delay.
    const double q = r.Q(kQ1, kClosed) + r.Q(kQ2, kClosed) + r.Q(kThink, kClosed);
    CHECK(q == doctest::Approx(3.0).epsilon(1e-9));

    // Cycle time is N/X less the think time actually spent at the delay. On an open
    // chain the reference divides Inf by the arrival rate and reports Inf, which is
    // what MATLAB's `C = N./X - sum(Z0,1)` gives; it is not a residence-time sum.
    CHECK(r.C[kClosed] == doctest::Approx(1.5634961771394).epsilon(kTol));
    CHECK(std::isinf(r.C[kOpen]));
}
