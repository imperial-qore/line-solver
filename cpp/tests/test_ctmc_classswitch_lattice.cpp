/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The population lattice is CAPACITY-BOUND, PERF-3.
 *
 * `capc` already zeroes a (node, class) pair the class never visits, and
 * `space_generator` uses it to drop the composed rows that breach it -- but it
 * dropped them ONE AT A TIME, after `space_closed_multi_cs` had spread every
 * class over every slot with no bound. On a class-switching chain
 * Source -> Q1(A) -> Q2(B) -> Q3(C), where each queue serves exactly one class
 * and the others are Disabled, that is (cutoff+1)^(3*3) candidates for
 * (cutoff+1)^3 reachable states -- and the rejection is not free, since each
 * candidate costs a `from_marginal`.
 *
 * The bound is now applied at the BRANCH: a slot of capacity 0 takes only the
 * zero item, so the enumeration is the reachable lattice rather than the full
 * one. The JAR has done this since `spaceClosedMultiCSBounded`; this is the
 * port catching up.
 *
 * WHAT MUST NOT CHANGE is the state space itself. The bound removes only rows
 * the composition was going to reject anyway, so the count and the metrics are
 * identical -- which is what the assertions below pin, alongside the count that
 * shows the lattice was the thing that shrank.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** Source -> Q1(A) -> Q2(B) -> Q3(C) -> Sink, one class switch per hop. */
qn::Network<double> switching_chain() {
    qn::Network<double> m("cs_chain");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t A = m.add_open_class("A");
    const std::size_t B = m.add_open_class("B");
    const std::size_t C = m.add_open_class("C");
    m.set_arrival(src, A, D::exp_rate(0.3));
    m.set_service(q1, A, D::exp_rate(0.7));
    m.set_service(q2, B, D::exp_rate(0.6));
    m.set_service(q3, C, D::exp_rate(0.5));
    qn::RoutingMatrix<double> P;
    P.set(A, A, src, q1, 1.0);
    P.set(A, B, q1, q2, 1.0);
    P.set(B, C, q2, q3, 1.0);
    P.set(C, C, q3, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the class-switching lattice enumerates only the reachable slots") {
    qn::Network<double> m = switching_chain();
    const qn::NetworkStruct<double> sn = m.get_struct();

    // Each queue holds one class at a time, so a state is (nA@Q1, nB@Q2, nC@Q3)
    // and the space is exactly (cutoff+1)^3. The Source contributes one row.
    for (int cutoff = 2; cutoff <= 4; ++cutoff) {
        CAPTURE(cutoff);
        const std::vector<std::size_t> cut(3, static_cast<std::size_t>(cutoff));
        const std::vector<qn::NetState<double>> space =
            qn::space_generator(sn, cut, static_cast<std::size_t>(1) << 22);
        const std::size_t want = static_cast<std::size_t>((cutoff + 1) * (cutoff + 1) * (cutoff + 1));
        CHECK(space.size() == want);
    }
}

TEST_CASE("the bounded lattice leaves the metrics alone") {
    // The bound removes only rows the composition rejected anyway, so every
    // reported number is unchanged. Flow balance is the sharpest statement
    // available on an open chain: one class per hop, so all three queues carry
    // the SAME throughput, and the Source carries it too.
    qn::Network<double> m = switching_chain();
    ctmc::CtmcOptions opt;
    opt.cutoff = 3;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);

    const double X = r.TN(0, 0);
    CHECK(X > 0.0);
    CHECK(r.TN(1, 0) == doctest::Approx(X).epsilon(1e-9));
    CHECK(r.TN(2, 1) == doctest::Approx(X).epsilon(1e-9));
    CHECK(r.TN(3, 2) == doctest::Approx(X).epsilon(1e-9));

    // A class only ever sits at its own queue.
    CHECK(r.QN(1, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.QN(1, 2) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.QN(2, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.QN(3, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // Utilization is X*E[S] at each single-server queue.
    CHECK(r.UN(1, 0) == doctest::Approx(X / 0.7).epsilon(1e-9));
    CHECK(r.UN(2, 1) == doctest::Approx(X / 0.6).epsilon(1e-9));
    CHECK(r.UN(3, 2) == doctest::Approx(X / 0.5).epsilon(1e-9));
}

TEST_CASE("the bound is what shrinks the lattice, not the composition") {
    // The DISCRIMINATING assertion: the composed state space was always
    // (cutoff+1)^3, because the composition rejected the impossible rows. What
    // the bound changes is how many candidates were generated to be rejected.
    // Three slots (Q1, Q2, Q3) and three classes, each class admissible at one
    // slot only.
    const std::size_t Mp = 3;
    std::vector<std::vector<bool>> chains(1, std::vector<bool>(3, true));  // one chain, A/B/C switch
    std::vector<std::vector<long>> caps(3, std::vector<long>(Mp, 0));
    caps[0][0] = 3;   // class A only at Q1
    caps[1][1] = 3;   // class B only at Q2
    caps[2][2] = 3;   // class C only at Q3

    const std::vector<std::size_t> n(3, 1);   // one job of each class in the chain
    const std::vector<std::vector<double>> unbounded =
        qn::space_closed_multi_cs<double>(Mp, n, chains);
    const std::vector<std::vector<double>> bounded =
        qn::space_closed_multi_cs<double>(Mp, n, chains, caps);

    CHECK(bounded.size() < unbounded.size());
    // Every bounded row must be one the unbounded enumeration also produced:
    // the bound removes candidates, it does not invent them.
    for (std::size_t i = 0; i < bounded.size(); ++i) {
        bool found = false;
        for (std::size_t j = 0; j < unbounded.size() && !found; ++j)
            if (bounded[i] == unbounded[j]) found = true;
        CHECK(found);
    }
    // And every bounded row respects the caps it was given.
    for (std::size_t i = 0; i < bounded.size(); ++i)
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t slot = 0; slot < Mp; ++slot)
                CHECK(bounded[i][r * Mp + slot] <= static_cast<double>(caps[r][slot]));
}
