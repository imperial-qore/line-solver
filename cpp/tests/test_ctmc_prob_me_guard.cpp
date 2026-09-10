/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression: the ME (matrix-exponential) gate on the four steady-state
 * probability functions.
 *
 * `assert_phase_type_states` (solver_ctmc_prob.h) refuses a query whose answer
 * would be a per-state or per-marginal probability read off a SIGNED
 * stationary vector -- which is what a matrix-exponential service process
 * produces. The reference (`@SolverCTMC/getProb*.m`) calls the guard from
 * INSIDE `getProb`, `getProbAggr`, `getProbSys` and `getProbSysAggr`, i.e. the
 * wrappers around `solver_ctmc_marg`, `_margaggr`, `_joint` and `_jointaggr`.
 * This port had the guard defined and correctly wired into every
 * `getTranProb*` accessor and `sampleSys`, but never into these four
 * steady-state functions themselves, so a caller reaching them directly (or
 * through a future `getProb*` wrapper) on an ME model got back `d.pi[s]` --
 * a signed number -- formatted as a plain probability with no complaint.
 */
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/util/error.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> FCFS Queue -> Delay, one closed class. */
qn::Network<double> cqn(double njobs) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * The same topology with a matrix-exponential service at the queue: D0 carries
 * a NEGATIVE off-diagonal, so the stationary vector is a signed measure.
 * Patched in place because no builder entry point produces an ME.
 */
qn::NetworkStruct<double>& make_me(qn::Network<double>& m) {
    qn::NetworkStruct<double>& sn = const_cast<qn::NetworkStruct<double>&>(m.get_struct());
    line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = -0.5;  // the signature of an ME rather than a PH
    D0(1, 1) = -3.0;
    D1(0, 1) = 2.5;
    D1(1, 0) = 3.0;
    sn.service[1][0].D0 = D0;
    sn.service[1][0].D1 = D1;
    sn.service[1][0].disabled = false;
    return sn;
}

}  // namespace

TEST_CASE("ctmc prob: an ME service model is refused by all four steady-state queries") {
    qn::Network<double> m = cqn(2.0);
    qn::NetworkStruct<double>& sn = make_me(m);
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    REQUIRE(!d.chain.space.empty());

    const line::qn::NetState<double>& s0 = d.chain.space[0];
    CHECK_THROWS_AS(ctmc::solver_ctmc_joint(sn, d, s0), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::solver_ctmc_jointaggr(sn, d, s0), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::solver_ctmc_marg(sn, d, s0), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::solver_ctmc_margaggr(sn, d, s0), line::UnsupportedError);
}

TEST_CASE("ctmc prob: an ordinary phase-type model passes the same four queries") {
    // The guard must not over-refuse: a plain exponential service is a
    // (trivial) phase-type process and every query below has to run to
    // completion and return a genuine probability.
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    REQUIRE(!d.chain.space.empty());

    const line::qn::NetState<double>& s0 = d.chain.space[0];
    double joint = 0, jointaggr = 0;
    std::vector<double> marg, margaggr;
    CHECK_NOTHROW(joint = ctmc::solver_ctmc_joint(sn, d, s0));
    CHECK_NOTHROW(jointaggr = ctmc::solver_ctmc_jointaggr(sn, d, s0));
    CHECK_NOTHROW(marg = ctmc::solver_ctmc_marg(sn, d, s0));
    CHECK_NOTHROW(margaggr = ctmc::solver_ctmc_margaggr(sn, d, s0));
    CHECK(joint >= 0.0);
    CHECK(jointaggr >= 0.0);
    CHECK(!marg.empty());
    CHECK(!margaggr.empty());
}
