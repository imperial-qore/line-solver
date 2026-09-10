/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * A fork on the same task that receives an entry-level open arrival.
 *
 * The Server task carries two entries with an AND fork/join each: SE is called
 * by the closed Client, OE takes an exogenous Poisson stream. Both live on one
 * task, so ONE layer holds a fork and a Source at once -- which the fork-join
 * transform cannot serve, because fj_mmt mints its own Source/Sink pair and
 * would detach the open stream already routed through the layer's own Source.
 * All four codebases refuse it by the same name (MATLAB buildLayersRecursive.m,
 * JAR SolverLN.java, python solver_ln.py, and build_fork_views() here).
 *
 * The refusal is asserted, not the numbers, so it cannot regress into a silent
 * wrong answer: python used to return a Server task throughput of 5e7 on this
 * model and the JAR a NullPointerException. lqns 6.2.28 DOES solve it (valid),
 * so this is a LINE feature gap, not an ill-formed model: Client throughput
 * 0.413391, Server task 0.513391, OE throughput 0.1 with open-wait 1.22917;
 * lqsim (T=5e5, seed 1234) gives 0.4154, 0.52242, 1.10546. Note that the open
 * entry takes NO reply activity: an open-arrival entry is send-no-reply, and
 * lqns rejects both a reply on it and a rendezvous sharing it.
 *
 * The flat fork-join with the same mix of an open and a closed class IS solved:
 * see test_fj_mixed_openclosed.cpp.
 */

#include <limits>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

/**
 * @param open_arrival  add the second entry (OE) with the Poisson stream; with
 *                      false the model is the same fork/join without any open
 *                      traffic, which is the control the solver does accept
 */
lqn::LqnStruct<double> build_fork_open(bool open_arrival) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", std::numeric_limits<double>::infinity(), SchedStrategy::INF);

    b.task("Client", 1, SchedStrategy::REF, "P1");
    b.think_time("Client", E(1.0));
    b.task("Server", 1, SchedStrategy::FCFS, "P2");
    b.think_time("Server", D::immediate());

    b.entry("CE", "Client");
    b.entry("SE", "Server");

    b.activity("CA", E(0.5), "Client");
    b.bound_to("CA", "CE");
    b.sync_call("CA", "SE", 1.0);

    b.activity("RA1", E(0.2), "Server");
    b.bound_to("RA1", "SE");
    b.activity("RA2", E(0.3), "Server");
    b.activity("RA3", E(0.4), "Server");
    b.activity("RA4", E(0.1), "Server");
    b.replies_to("RA4", "SE");
    b.and_fork("RA1", {"RA2", "RA3"});
    b.and_join({"RA2", "RA3"}, "RA4");

    if (open_arrival) {
        b.entry("OE", "Server");
        b.open_arrival("OE", E(10.0));
        b.activity("OA1", E(0.2), "Server");
        b.bound_to("OA1", "OE");
        b.activity("OA2", E(0.3), "Server");
        b.activity("OA3", E(0.4), "Server");
        b.activity("OA4", E(0.1), "Server");
        b.and_fork("OA1", {"OA2", "OA3"});
        b.and_join({"OA2", "OA3"}, "OA4");
    }
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("a layer with a fork and an entry open arrival is refused by name") {
    const lqn::LqnStruct<double> l = build_fork_open(true);
    // the model itself is well formed: two entries on Server, one of them open
    CHECK(l.nentries == 3);
    CHECK(l.nacts == 9);
    CHECK(l.has_arrival[idx_of(l, "E:OE")]);

    ln::LnOptions opt;
    // The refusal is the ROUTING encoding's: it is fj_mmt's own Source that
    // collides with the layer's. 'srvn.ph' composes the fork into the entry law
    // and never mints one, so it SOLVES this model, and the 'srvn' alias now
    // resolves to it here -- which is why the encoding is named rather than
    // taken. MATLAB refuses and solves it the same way round.
    opt.method = "srvn.cs";
    CHECK_THROWS_AS(ln::SolverLN<double>(l, opt), UnsupportedError);

    // the alias reaches the phase-type encoding, which does serve it: MATLAB
    // with method='srvn.ph' gives Client 0.4, Server task 0.554545, OE 0.110909
    ln::LnOptions ph;
    ln::SolverLN<double> s(l, ph);
    CHECK(s.state_lnmethod() == "srvn.ph");
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);
    CHECK(sol.TN[idx_of(l, "R:Client")] == doctest::Approx(0.4).epsilon(1e-3));
    CHECK(sol.TN[idx_of(l, "T:Server")] == doctest::Approx(0.554545).epsilon(1e-3));
    CHECK(sol.TN[idx_of(l, "E:OE")] == doctest::Approx(0.110909).epsilon(1e-3));
}

TEST_CASE("the same fork/join without the open stream is accepted") {
    const lqn::LqnStruct<double> l = build_fork_open(false);
    ln::LnOptions opt;
    // a Fork node is a routing-encoding object, as above
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    // Client and Server layers, plus the two host layers collapse to the hosts
    // that carry tasks: the fork lives in the Server layer and is transformed
    CHECK(s.nlayers() >= 2);
    bool anyfork = false;
    for (const qn::Layer<double>& L : s.layers())
        if (L.has_fork()) anyfork = true;
    CHECK(anyfork);
}
