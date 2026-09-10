/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Elimination of the VANISHING states, BUG-98.
 *
 * This port did not eliminate them AT ALL. `hide_immediate` had no counterpart,
 * `ctmc_find_vanishing_states` was not ported, and `solver_ctmc_ratecomplement`
 * was ported but had ZERO callers and said so in its own comment. A Router or
 * Fork holding a job, a Join whose sibling set is complete, and an immediate SPN
 * mode therefore all kept their GlobalConstants::Immediate row, so `-a states`
 * and `-a gen` returned a state space the other three codebases had already
 * complemented away.
 *
 * THE MEASURED SYMPTOM, and the repro this file uses: `fj_tiny_closed` --
 * Delay + Fork + two FCFS queues + Join, one closed class with N=1 -- printed
 * SIX states where native python's generator is 4x4, the two extras carrying
 * 6.12245e-09 of the mass each against the four tangible 0.612245, 0.122449,
 * 0.081633 and 0.183673. `-a avg` agreed to the digits printed for exactly that
 * reason, which is why the gap stayed invisible unless the chain itself was
 * asked for.
 *
 * The mean metrics are pinned alongside the state count, because eliminating the
 * states is only half the change: an action that fires ONLY from vanishing
 * states -- the fork firing and the join rendezvous both do -- has no tangible
 * row to be read off, so `arv_rates` and `dep_rates` have to be RATE
 * COMPLEMENTED back onto the tangible states or the throughput would silently
 * drop to zero. The closed form for N=1 is
 * X = 1/(Z + E[max(S1,S2)]) = 1/(1 + 1/2 + 1/3 - 1/5) = 0.612244898.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** Delay -> Fork -> {Queue1, Queue2} -> Join -> Delay, one closed job. */
qn::Network<double> fj_tiny_closed() {
    qn::Network<double> m("fj_tiny_closed");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t f = m.add_fork("Fork1");
    const std::size_t j = m.add_join("Join1", f);
    const std::size_t c = m.add_closed_class("Class1", 1.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Router -> Queue -> Delay, one closed job. */
qn::Network<double> router_closed() {
    qn::Network<double> m("router_closed");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t r = m.add_router("Router1");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, r, 1.0);
    P.set(c, c, r, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the fork-join chain keeps only its tangible states") {
    qn::Network<double> m = fj_tiny_closed();
    const ctmc::CtmcSolution<double> s =
        ctmc::solver_ctmc_analyzer(m.get_struct(), ctmc::CtmcOptions());

    // Four tangible states, not the six the enumeration produces: the
    // fork-occupied and the join-firable rows are vanishing and are
    // complemented out.
    CHECK(s.chain.space.size() == 4);
    CHECK(s.chain.Q.rows() == 4);
    CHECK(s.pi.size() == 4);

    // No 1e-8-sojourn row survives: every remaining state carries real mass.
    double total = 0.0;
    for (std::size_t i = 0; i < s.pi.size(); ++i) {
        CHECK(s.pi[i] > 1e-4);
        total += s.pi[i];
    }
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));

    // The generator no longer carries a GlobalConstants::Immediate rate.
    for (std::size_t i = 0; i < s.chain.Q.rows(); ++i)
        for (std::size_t j = 0; j < s.chain.Q.cols(); ++j)
            CHECK(std::fabs(s.chain.Q(i, j)) < 1e6);
}

TEST_CASE("the fork-join mean metrics survive the elimination") {
    // The rate complement is what this pins. A fork firing and a join rendezvous
    // fire ONLY from vanishing states, so restricting arv_rates/dep_rates to the
    // tangible rows alone would zero the throughput outright.
    qn::Network<double> m = fj_tiny_closed();
    const line::mva::AvgResult<double> r =
        ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());

    const double X = 0.612244898;   // 1/(Z + E[max(S1,S2)]), N = 1
    CHECK(r.QN(0, 0) == doctest::Approx(0.612244898).epsilon(1e-9));  // Delay
    CHECK(r.QN(1, 0) == doctest::Approx(0.306122449).epsilon(1e-9));  // Queue1
    CHECK(r.QN(2, 0) == doctest::Approx(0.204081633).epsilon(1e-9));  // Queue2
    CHECK(r.TN(0, 0) == doctest::Approx(X).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(X).epsilon(1e-9));
    CHECK(r.TN(2, 0) == doctest::Approx(X).epsilon(1e-9));
}

TEST_CASE("a router pass-through leaves no vanishing state behind") {
    // A Router is stateful but performs no service, so a job at one is in
    // transit: every state that holds one is vanishing. The chain must be the
    // plain Delay+Queue chain of three levels, N = 2.
    qn::Network<double> m = router_closed();
    const ctmc::CtmcSolution<double> s =
        ctmc::solver_ctmc_analyzer(m.get_struct(), ctmc::CtmcOptions());
    CHECK(s.chain.space.size() == 3);
    for (std::size_t i = 0; i < s.pi.size(); ++i) CHECK(s.pi[i] > 1e-4);

    // And the metrics are the exact closed Delay+Queue ones: Z = 1, S = 0.5,
    // N = 2 gives X = 1.2, Q_queue = 0.8, R_queue = 2/3.
    const line::mva::AvgResult<double> r =
        ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());
    CHECK(r.TN(0, 0) == doctest::Approx(1.2).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(0.8).epsilon(1e-9));
}
