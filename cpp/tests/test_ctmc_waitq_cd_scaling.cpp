/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression: class/joint-dependent service rate scaling inside the WAITQ
 * (finite capacity region) generator.
 *
 * The default generator (`solver_ctmc.h`) and the WAITQ generator
 * (`solver_ctmc_waitq.h`) both dispatch every active-node transition through
 * the SAME `qn::after_event` -> `after_event_station` ->
 * `after_event_station_dep` / `after_event_station_phase` pipeline
 * (`state_events.h`), and `cd_factor` is folded into the rate INSIDE that
 * pipeline (state_events.h:831-832, :1298), not re-applied by either caller.
 * Grepping the WAITQ file for the literal text "cd_factor" therefore finds
 * nothing even though the scaling is exercised on every transition the WAITQ
 * walk emits, because it never needed a call site of its own: it is baked
 * into the `EventOutcome::rate` the shared dispatcher returns. MATLAB is the
 * same shape -- `solver_ctmc_fcr_waitq.m` calls `State.afterEventHashed`,
 * which calls the identical `State.afterEvent` the default generator uses, and
 * grepping `solver_ctmc_fcr_waitq.m` for cdscaling/jdscaling likewise finds
 * nothing.
 *
 * The oracle: wrap a class-dependent, joint-dependent PS queue in a WAITQ
 * region whose caps are all the -1 (unbounded) sentinel. An unbounded region
 * refuses nothing and parks nothing, so the augmented chain must reproduce
 * the region-less chain's means EXACTLY -- and it can only do that if the
 * dependence-scaled rates it walks on agree with the default generator's,
 * state for state. A WAITQ path that walked on UNSCALED rates would still
 * produce a normalized, plausible-looking solution, but a DIFFERENT one; this
 * test would catch exactly that regression.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::DropStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/**
 * Think -> PS Queue -> Think, one closed class, four jobs. The queue carries
 * BOTH a class-dependence handle beta(n) = min(n, 2) (peak 2) and a
 * joint-dependence handle eta(n) = 3 flat (peak 3), so every DEP/PHASE rate
 * `after_event_station_dep` emits at the queue is scaled by both factors.
 */
qn::Network<double> cdjd_cqn(double njobs) {
    qn::Network<double> m("cdjd_cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(1.5));
    m.set_class_dependence(q, [](const std::vector<double>& n) {
        return std::vector<double>(1, std::min(n[0], 2.0));
    }, std::vector<double>(1, 2.0));
    m.set_joint_dependence(q, [](const std::vector<double>&) {
        return std::vector<double>(1, 3.0);
    }, std::vector<double>(1, 3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("ctmc waitq: an unbounded WAITQ region preserves class/joint-dependent rates") {
    const double N = 4.0;
    qn::Network<double> plain = cdjd_cqn(N);
    qn::Network<double> regd = cdjd_cqn(N);
    // Every node is a member and every cap is -1 (unbounded), so nothing is
    // ever refused or parked -- the region is inert on the STATE SPACE, which
    // is exactly what isolates rate scaling as the only remaining way the two
    // solutions could differ.
    regd.add_region(std::vector<std::size_t>{1, 2}, std::vector<double>{-1.0}, -1.0,
                    std::vector<DropStrategy>{DropStrategy::WAITQ});

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(plain.get_struct(), opt);
    const qn::NetworkStruct<double>& sn1 = regd.get_struct();
    REQUIRE(ctmc::ctmc_has_waitq_region(sn1));
    const ctmc::WaitqSolution<double> d1 = ctmc::solver_ctmc_waitq_analyzer(sn1, opt);

    for (std::size_t i = 0; i < sn1.nstations; ++i) {
        CHECK(d1.sol.avg.QN(i, 0) == doctest::Approx(d0.avg.QN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.UN(i, 0) == doctest::Approx(d0.avg.UN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.TN(i, 0) == doctest::Approx(d0.avg.TN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.RN(i, 0) == doctest::Approx(d0.avg.RN(i, 0)).epsilon(1e-9));
    }
    CHECK(d1.sol.avg.XN[0] == doctest::Approx(d0.avg.XN[0]).epsilon(1e-9));

    // No token was ever parked, confirming the region truly stayed inert
    // rather than the two throughputs agreeing on a WAITQ chain that quietly
    // diverged elsewhere.
    double parked = 0;
    for (std::size_t r = 0; r < d1.parked.size(); ++r) parked += d1.parked[r];
    CHECK(parked == doctest::Approx(0.0).epsilon(1e-12));
}
