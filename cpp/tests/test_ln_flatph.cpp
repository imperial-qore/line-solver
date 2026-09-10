/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Method `flat.ph`: the squashed layering with the composed phase-type encoding.
 *
 * ONE submodel holds a station for every processor and every called task, and a
 * caller task is one closed class that visits each server it uses once per
 * invocation, carrying there the composed law of the demand it places on that
 * server. It is the same composition `srvn.ph` performs; what changes is that
 * the servers contend inside one network instead of meeting each other through
 * surrogate delays, so the client delay keeps only the think times.
 *
 * The reference values are MATLAB's, and the two models are the twins the JAR
 * and native-Python tests use, so the four codebases are pinned to one set of
 * numbers. The tolerance is the band the codebases' AMVA fixed points already
 * differ by, not machine precision: measured 3e-6 on the two-task model and
 * 2.3e-4 on the three-deep chain.
 */

#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** T1 -> T2, a single call level with a three-thread callee. */
lqn::LqnStruct<double> build_twotasks() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 3, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(0.5), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 2.0);
    b.activity("A2", D::exp_mean(1.0 / 3.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    return b.build();
}

/** T1 -> T2 -> T3, so the middle task is at once a server and a caller. */
lqn::LqnStruct<double> build_serialchain() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 4, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.2), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", D::exp_mean(0.2), "T2");
    b.bound_to("A2", "E2");
    b.sync_call("A2", "E3", 1.0);
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.2), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

ln::LnSolution<double> solve(const lqn::LqnStruct<double>& l, const std::string& method) {
    ln::LnOptions opt;
    opt.method = method;
    return ln::SolverLN<double>(l, opt).get_ensemble_avg();
}

}  // namespace

TEST_CASE("flat.ph: one layer holds every server, with one class per caller") {
    const lqn::LqnStruct<double> l = build_serialchain();
    ln::LnOptions opt;
    opt.method = "flat.ph";
    ln::SolverLN<double> s(l, opt);
    s.get_ensemble_avg();

    // the label reports what was BUILT, and the alias 'flat' does not reach it
    CHECK(s.state_lnmethod() == "flat.ph");
    // ONE submodel, not one per server
    CHECK(s.nlayers() == 1);
    const qn::Layer<double>& L = s.layers()[0];
    CHECK(L.flat);
    // three processors and two called tasks, all stations of that one model
    CHECK(L.host_stations.size() == 3);
    CHECK(L.task_stations.size() == 2);
    // one closed class per caller task, not one per entry, activity and call:
    // that collapse is the whole point of the composed encoding
    CHECK(L.classes.size() == 3);
    // every served element resolves to a station of its own
    CHECK(L.server_idx_of[idx_of(l, "P:P1")] != L.server_idx_of[idx_of(l, "P:P2")]);
    CHECK(L.server_idx_of[idx_of(l, "T:T2")] != L.server_idx_of[idx_of(l, "T:T3")]);
}

TEST_CASE("flat.ph is the MATLAB answer") {
    {
        const lqn::LqnStruct<double> l = build_twotasks();
        const ln::LnSolution<double> sol = solve(l, "flat.ph");
        CHECK(sol.TN[idx_of(l, "R:T1")] == doctest::Approx(1.391663).epsilon(1e-3));
        CHECK(sol.TN[idx_of(l, "T:T2")] == doctest::Approx(2.783326).epsilon(1e-3));
        CHECK(sol.UN[idx_of(l, "P:P1")] == doctest::Approx(0.695830).epsilon(1e-3));
        CHECK(sol.UN[idx_of(l, "P:P2")] == doctest::Approx(0.927773).epsilon(1e-3));
    }
    {
        const lqn::LqnStruct<double> l = build_serialchain();
        const ln::LnSolution<double> sol = solve(l, "flat.ph");
        CHECK(sol.TN[idx_of(l, "R:T1")] == doctest::Approx(2.126442).epsilon(1e-3));
        CHECK(sol.TN[idx_of(l, "T:T2")] == doctest::Approx(2.126442).epsilon(1e-3));
        CHECK(sol.TN[idx_of(l, "T:T3")] == doctest::Approx(2.126442).epsilon(1e-3));
    }
}

TEST_CASE("the flat alias stays on the routing encoding") {
    // The composed law folds every call of an invocation into ONE visit, so the
    // dispatch order of a routed call group has nowhere to be expressed.
    // Squashing does not recover it, which is why 'flat' does not probe 'flat.ph'.
    const lqn::LqnStruct<double> l = build_serialchain();
    ln::LnOptions opt;
    opt.method = "flat";
    ln::SolverLN<double> s(l, opt);
    s.get_ensemble_avg();
    CHECK(s.state_lnmethod() == "flat.cs");
}

TEST_CASE("flat.ph refuses what a single submodel cannot carry") {
    // a replicated task: the copies need a station each, and only a submodel of
    // its own can hold both readings
    lqn::LqnStruct<double> l = build_serialchain();
    l.repl[idx_of(l, "T:T2")] = 2.0;
    ln::LnOptions opt;
    opt.method = "flat.ph";
    CHECK_THROWS_AS(ln::SolverLN<double>(l, opt).get_ensemble_avg(), UnsupportedError);
}
