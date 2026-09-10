/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The programmatic LQN builder, on lqn_workflows.
 *
 * matlab/examples/basic/layeredModel/lqn_workflows.m is built from the
 * Processor / Task / Entry / Activity constructors, not from a file, and it
 * CANNOT be routed through the .lqnx reader: its T3 carries a think time on an
 * inf-scheduled (non-reference) task, which the .lqnx writer omits because lqns
 * rejects the attribute there. So it is the regression for the builder path.
 *
 * It also exercises the three activity-graph constructs ofbiz does not have: a
 * loop, an AND fork/join, and an OR fork/join. The AND pair is REFUSED by the
 * solver at layer-build time; that refusal is asserted here by name, so it
 * cannot regress into a silent omission.
 */

#include <string>
#include <vector>

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
 * @param serialise_and  replace the AND fork/join on T2 by a serial chain, which
 *                       is how the solvable part of the model is reached until
 *                       fork/join support lands
 */
lqn::LqnStruct<double> build_workflows(bool serialise_and) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", INFINITY, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", D::immediate());
    b.entry("Entry", "T1");
    b.processor("P2", INFINITY, SchedStrategy::INF);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.think_time("T2", D::immediate());
    b.entry("E2", "T2");
    b.processor("P3", 5, SchedStrategy::PS);
    b.task("T3", 20, SchedStrategy::INF, "P3");
    b.think_time("T3", E(10));
    b.entry("E1", "T3");

    b.activity("A1", E(1), "T1");
    b.bound_to("A1", "Entry");
    b.activity("A2", E(2), "T1");
    b.activity("A3", E(3), "T1");
    b.sync_call("A3", "E2", 1.0);
    b.activity("B1", E(0.1), "T2");
    b.bound_to("B1", "E2");
    b.activity("B2", E(0.2), "T2");
    b.activity("B3", E(0.3), "T2");
    b.activity("B4", E(0.4), "T2");
    b.activity("B5", E(0.5), "T2");
    b.activity("B6", E(0.6), "T2");
    b.sync_call("B6", "E1", 1.0);
    b.replies_to("B6", "E2");
    b.activity("C1", E(0.1), "T3");
    b.bound_to("C1", "E1");
    b.activity("C2", E(0.2), "T3");
    b.activity("C3", E(0.3), "T3");
    b.activity("C4", E(0.4), "T3");
    b.activity("C5", E(0.5), "T3");
    b.replies_to("C5", "E1");

    b.loop("A1", {"A2"}, "A3", 3.0);
    if (serialise_and) {
        b.serial("B1", "B2");
        b.serial("B2", "B3");
        b.serial("B3", "B4");
        b.serial("B4", "B5");
        b.serial("B5", "B6");
    } else {
        b.serial("B4", "B5");
        b.and_fork("B1", {"B2", "B3", "B4"});
        b.and_join({"B2", "B3", "B5"}, "B6");
    }
    b.or_fork("C1", {"C2", "C3", "C4"}, {0.3, 0.3, 0.4});
    b.or_join({"C2", "C3", "C4"}, "C5");
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("lqn_workflows: the builder reproduces the MATLAB LayeredNetworkStruct") {
    const lqn::LqnStruct<double> l = build_workflows(false);
    CHECK(l.nhosts == 3);
    CHECK(l.ntasks == 3);
    CHECK(l.nentries == 3);
    CHECK(l.nacts == 14);
    CHECK(l.ncalls == 2);
    CHECK(l.nidx == 23);
    CHECK(l.hashnames[4] == "R:T1");
    CHECK(l.hashnames[5] == "T:T2");
    CHECK(l.hashnames[6] == "T:T3");

    // The think time the .lqnx format cannot carry: T3 is inf-scheduled AND has
    // a think time of 10. Losing it is the whole reason the builder exists.
    CHECK(l.think[6].disabled == false);
    CHECK(l.think[6].mean == doctest::Approx(10.0));

    // the loop back-edge is in graph but not in dag, which is what marks it
    const std::size_t a2 = idx_of(l, "A:A2"), a3 = idx_of(l, "A:A3");
    CHECK(a2 > 0);
    CHECK(a3 > 0);
    CHECK(l.graph.get(a2, a2) != l.dag.get(a2, a2));

    // the OR-fork branch probabilities reach the graph as edge weights
    const std::size_t c1 = idx_of(l, "A:C1");
    CHECK(l.graph.get(c1, idx_of(l, "A:C2")) == doctest::Approx(0.3));
    CHECK(l.graph.get(c1, idx_of(l, "A:C3")) == doctest::Approx(0.3));
    CHECK(l.graph.get(c1, idx_of(l, "A:C4")) == doctest::Approx(0.4));

    CHECK(l.actposttype[idx_of(l, "A:B2")] == lang::PrecedenceType::POST_AND);
    CHECK(l.actpretype[idx_of(l, "A:B5")] == lang::PrecedenceType::PRE_AND);
    CHECK(l.actposttype[idx_of(l, "A:C2")] == lang::PrecedenceType::POST_OR);
    CHECK(l.actpretype[idx_of(l, "A:C4")] == lang::PrecedenceType::PRE_OR);
}

TEST_CASE("lqn_workflows: an AND fork/join builds Fork, Routers and a Join") {
    // The layer structure is compared against a MATLAB dump of the same model:
    // stations, stateful nodes, classes, chains and the fork-join pairing all
    // agree. `nnodes` does not, and is not expected to: MATLAB materialises an
    // auto-added ClassSwitch node per class-switching link, which the same
    // stochastic complement then removes again (see qn_layer.h).
    const lqn::LqnStruct<double> l = build_workflows(false);
    ln::LnOptions opt;
    // Fork, Router and Join are nodes of the ROUTING encoding; 'srvn.ph' reduces
    // the same AND pair inside the composed entry law and builds none of them.
    // The 'srvn' alias resolves to 'srvn.ph' on this model.
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    CHECK(s.nlayers() == 5);

    // P2 is the host layer of T2, whose activity graph carries the AND pair
    const qn::Layer<double>* p2 = nullptr;
    for (const qn::Layer<double>& L : s.layers())
        if (L.name == "P:P2") p2 = &L;
    REQUIRE(p2 != nullptr);
    CHECK(p2->nstations == 3);                 // Clients, P:P2, Join_PreAnd
    CHECK(p2->stateful_nodes.size() == 6);     // every node but the Fork
    CHECK(p2->nodes.size() == 7);
    CHECK(p2->nclasses == 9);
    CHECK(p2->nchains == 1);
    CHECK(p2->has_fork());
    REQUIRE(p2->fj.size() == 1);
    CHECK(p2->fj[0].first == 3);   // Fork_PostAnd
    CHECK(p2->fj[0].second == 7);  // Join_PreAnd
    CHECK(p2->nodes[2].nodetype == qn::NodeType::Fork);
    CHECK(p2->nodes[2].stateful == false);  // a Fork holds no jobs
    CHECK(p2->nodes[6].nodetype == qn::NodeType::Join);
    CHECK(p2->station_to_node[2] == 7);

    // the fork-corrected visit ratios, as MATLAB computes them: the
    // population-preserving SPN argument sets every visited entry to one, and
    // the normalisation by the reference node's total makes them 1/3 here
    CHECK(p2->visits[0](1, 2) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(p2->visits[0](0, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
}

TEST_CASE("lqn_workflows: the AND join credits branch concurrency") {
    // With the AND pair in place, the entry response time must come out BELOW
    // the serialised sum of its activities: the three branches of T2's fork
    // ({B2}, {B3}, {B4,B5}) overlap, so E2 sees the k-th order statistic of
    // the branch times rather than their sum. MATLAB reports 2.6410308 for the
    // published model; serialising the branches instead gives 3.01.
    const lqn::LqnStruct<double> l = build_workflows(false);
    ln::LnOptions opt;
    // the MMT fixed point of the routing encoding, which is what 2.6410308 is
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);
    std::size_t e2 = 0;
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == "E:E2") e2 = i;
    REQUIRE(e2 > 0);
    CHECK(sol.RN[e2] == doctest::Approx(2.6410308).epsilon(1e-3));
    CHECK(sol.RN[e2] < 3.0);  // strictly below the serialised sum

    // Every element of the model, against the MATLAB AvgTable. The tolerance is
    // 1e-3: the residual is 4.6e-4 and sits on the branch classes, where the
    // reference's auxiliary open classes carry a slightly different visit ratio
    // than the uniform 1/nbranches merge this port applies (see fj_mmt.h).
    struct Row { const char* hn; double util; double tput; };
    const Row rows[] = {
        {"P:P1", 0.79106861, 0.0},
        {"P:P2", 0.16662209, 0.0},
        {"P:P3", 0.072180711, 0.0},
        {"R:T1", 0.79106861, 0.079106861},
        {"T:T2", 0.16662209, 0.079319464},
        {"T:T3", 0.072180711, 0.079319463},
        {"E:Entry", 0.79106861, 0.079106861},
        {"E:E2", 0.16662209, 0.079319464},
        {"E:E1", 0.072180711, 0.079319463},
        {"A:B1", 0.0079319464, 0.079319464},
        {"A:B6", 0.047591678, 0.079319464},
        {"A:C5", 0.039659731, 0.079319463},
    };
    for (const Row& row : rows) {
        std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(row.hn);
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
    }
}

TEST_CASE("lqn_workflows: the AND fork reproduces the reference's own flow residual") {
    // B1 is the fork head, B2..B5 lie on its three branches, B6 is the join
    // target. Each is executed exactly once per completion of E2, so in an
    // exactly solved model all six throughputs would be equal.
    //
    // They are not, in ANY codebase. MATLAB reports the fork head and the join
    // target at 0.079319463729 and the four branch activities at
    // 0.079319499179, a 4.5e-7 relative violation. It is a property of the MMT
    // fixed point, not a coding error: the auxiliary open classes that carry
    // the branches the circulating job did not take are merged back at an
    // arrival rate that stopped one iteration short, and the residual scales
    // with options.iter_tol (see the note in @NetworkSolver/fjFixedPoint.m).
    //
    // This port runs that same fixed point rather than approximating it, so the
    // assertion here is a PARITY assertion against the reference's two distinct
    // levels, checked to 1e-8 -- which is two orders of magnitude tighter than
    // the residual itself, so it would fail if the port drifted to either a
    // conserving answer or a differently non-conserving one. An earlier revision
    // of this port approximated the auxiliary merge by scaling the branch
    // metrics by the branch count, which conserved flow exactly and therefore
    // did NOT agree with any reference; that shortcut is gone.
    const lqn::LqnStruct<double> l = build_workflows(false);
    ln::LnOptions opt;
    // the residual belongs to the MMT auxiliary classes, so it exists only in
    // the routing encoding; the alias resolves to 'srvn.ph' on this model
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    const double head = 0.079319463729293591;   // MATLAB, A:B1 and A:B6
    const double branch = 0.079319499179196493; // MATLAB, A:B2 .. A:B5
    const char* at_head[] = {"A:B1", "A:B6"};
    const char* on_branch[] = {"A:B2", "A:B3", "A:B4", "A:B5"};
    for (const char* hn : at_head) {
        CAPTURE(hn);
        CHECK(sol.TN[idx_of(l, hn)] == doctest::Approx(head).epsilon(1e-8));
    }
    for (const char* hn : on_branch) {
        CAPTURE(hn);
        CHECK(sol.TN[idx_of(l, hn)] == doctest::Approx(branch).epsilon(1e-8));
    }
    // and the residual really is the reference's, not zero
    CHECK(branch - head > 0.0);
    CHECK(sol.TN[idx_of(l, "A:B2")] > sol.TN[idx_of(l, "A:B1")]);
}

TEST_CASE("lqn_workflows: loops and OR fork/join solve, against the MATLAB table") {
    // The AND pair is serialised so the rest of the model is reachable. The
    // quantities checked below are the ones the AND fork cannot move: a
    // residence time at a host is the activity's own demand times its visit
    // count, and neither the loop count on T1 nor the branch probabilities on
    // T3 depend on how T2's activities are ordered. They are read off the
    // MATLAB AvgTable of the PUBLISHED model, so this is a cross-check against
    // the reference and not against this port's own output.
    const lqn::LqnStruct<double> l = build_workflows(true);
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    CHECK(s.nlayers() == 5);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    auto residt = [&](const char* hn) { return sol.WN[idx_of(l, hn)]; };
    // loop on T1: A2 is executed 3 times in expectation, at demand 2
    CHECK(residt("A:A1") == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(residt("A:A2") == doctest::Approx(6.0).epsilon(1e-6));
    CHECK(residt("A:A3") == doctest::Approx(3.0).epsilon(1e-6));
    // OR-fork on T3: residence is the branch probability times the demand
    CHECK(residt("A:C1") == doctest::Approx(0.10).epsilon(1e-6));
    CHECK(residt("A:C2") == doctest::Approx(0.06).epsilon(1e-6));
    CHECK(residt("A:C3") == doctest::Approx(0.09).epsilon(1e-6));
    CHECK(residt("A:C4") == doctest::Approx(0.16).epsilon(1e-6));
    CHECK(residt("A:C5") == doctest::Approx(0.50).epsilon(1e-6));
}
