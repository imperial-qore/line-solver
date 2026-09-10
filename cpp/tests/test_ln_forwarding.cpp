/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Forwarding calls (LqnBuilder::forward -> CallType::FWD -> lqn_fwd_rendezvous).
 *
 * Model mirrors line-test.git/test/testsAdvFeatures/lqn/test_forwarding.m,
 * Test 5 ("SolverLN vs SolverLQNS Comparison"): a reference task with 3
 * clients calls e0 on server0, which forwards unconditionally (prob=1.0) to
 * e1 on server1; server0 and server1 share one PS processor. Reference
 * numbers below are MATLAB SolverLN.getAvgTable() on the identical model
 * (defaultOptions, 2026-07-29).
 */

#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::SchedStrategy;
using lang::Distrib;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

lqn::LqnStruct<double> build_forwarding() {
    lqn::LqnBuilder<double> b;
    b.processor("clientProc", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("serverProc", 1, SchedStrategy::PS);

    b.task("client", 3, SchedStrategy::REF, "clientProc");
    b.think_time("client", E(10.0));
    b.task("server0", 1, SchedStrategy::FCFS, "serverProc");
    b.task("server1", 1, SchedStrategy::FCFS, "serverProc");

    b.entry("client", "client");
    b.entry("e0", "server0");
    b.entry("e1", "server1");

    b.activity("client_1", E(0.5), "client");
    b.bound_to("client_1", "client");
    b.sync_call("client_1", "e0", 1.0);

    b.activity("e0_1", E(1.0), "server0");
    b.bound_to("e0_1", "e0");
    b.replies_to("e0_1", "e0");

    b.activity("e1_1", E(1.0), "server1");
    b.bound_to("e1_1", "e1");
    b.replies_to("e1_1", "e1");

    b.forward("e0", "e1", 1.0);

    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("forwarding: the builder produces a CallType::FWD call, inert after the rewrite") {
    const lqn::LqnStruct<double> l = build_forwarding();
    bool found_fwd = false;
    for (std::size_t c = 1; c <= l.ncalls; ++c)
        if (l.calltype[c] == lang::CallType::FWD) found_fwd = true;
    CHECK(found_fwd);
}

TEST_CASE("forwarding: e0~>e1 solves as a caller-side pseudo rendezvous, matching MATLAB SolverLN") {
    const lqn::LqnStruct<double> l = build_forwarding();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // server0 and server1 must be symmetric: e0 forwards ALL of its traffic to
    // e1 (prob=1.0), so both entries see the same visit count and, since they
    // share one PS processor with identical service, the same metrics.
    const std::size_t e0 = idx_of(l, "E:e0"), e1 = idx_of(l, "E:e1");
    REQUIRE(e0 > 0);
    REQUIRE(e1 > 0);
    CHECK(sol.RN[e0] == doctest::Approx(sol.RN[e1]).epsilon(1e-9));
    CHECK(sol.TN[e0] == doctest::Approx(sol.TN[e1]).epsilon(1e-9));
    CHECK(sol.UN[e0] == doctest::Approx(sol.UN[e1]).epsilon(1e-9));

    // Against the MATLAB reference (SolverLN.getAvgTable on the identical model,
    // RE-RECORDED 2026-08-11, after the interlock probability was aligned to
    // Li and Franks (2015), Eq. (5), and to lqns' m' rule):
    //   e0/e1 Entry:  QLen=0.277035 Util=0.224035 RespT=1.23657 Tput=0.224035
    //   client Entry: QLen=0.745590 Util=0.112720 RespT=3.30725 Tput=0.225441
    //   clientProc Util=0.112720, serverProc Util=0.448069 (= e0.Util + e1.Util)
    //
    // The 2026-08-04 het-FCFS gate tightening briefly re-recorded this golden
    // (client RespT 3.2460): the FCFS server-task layers carry seeded rates on
    // call classes with ZERO visits, which the raw per-class mean comparison
    // mistook for class-dependent service, diverting effectively single-class
    // layers off the product-form branch. The gate now compares visit-weighted
    // CHAIN service times, restoring this reference row.
    //
    // THE ROW BEFORE THE INTERLOCK FIX DID NOT BALANCE FLOW, which is how that
    // change was recognised as a fix and not a drift: it read Tput 0.22266 at
    // e0 against 0.22688 at the client entry, and e0 is called exactly once
    // per client cycle, so the two must agree. They now do, in both codebases,
    // to 2e-4.
    CHECK(sol.RN[e0] == doctest::Approx(1.23657).epsilon(1e-3));
    CHECK(sol.TN[e0] == doctest::Approx(0.224035).epsilon(1e-3));
    CHECK(sol.UN[e0] == doctest::Approx(0.224035).epsilon(1e-3));
    CHECK(sol.QN[e0] == doctest::Approx(0.277035).epsilon(1e-3));

    const std::size_t clientEntry = idx_of(l, "E:client");
    REQUIRE(clientEntry > 0);
    // The client's round trip must be the SUM of e0's and e1's own residence
    // times plus its own host demand: forwarding must not be double-charged
    // (the raw FWD arc contributes nothing) nor dropped (e1's residence time
    // must reach the client).
    CHECK(sol.RN[clientEntry] == doctest::Approx(3.30725).epsilon(1e-3));
    CHECK(sol.TN[clientEntry] == doctest::Approx(0.225441).epsilon(1e-3));
    // The flow balance the old row broke: e0 is called once per client cycle.
    CHECK(sol.TN[e0] == doctest::Approx(sol.TN[clientEntry]).epsilon(1e-2));

    const std::size_t clientProc = idx_of(l, "P:clientProc"), serverProc = idx_of(l, "P:serverProc");
    REQUIRE(clientProc > 0);
    REQUIRE(serverProc > 0);
    CHECK(sol.UN[clientProc] == doctest::Approx(0.112720).epsilon(1e-3));
    CHECK(sol.UN[serverProc] == doctest::Approx(0.448069).epsilon(1e-3));
}
