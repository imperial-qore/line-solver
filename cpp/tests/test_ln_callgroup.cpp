/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Routed call groups: one dispatch with n destinations.
 *
 * `sync_call_round_robin` / `sync_call_jsq` issue ONE call per invocation whose
 * destination cycles over, or is chosen among, several target entries. The
 * members are ordinary sync calls of mean `total/n`, so the aggregate call rate
 * is the same as the probabilistic twin's; what the group removes is the
 * VARIANCE of the branching, and that is the only thing to test for -- a
 * comparison of mean call rates would pass on a model with no dispatcher at
 * all.
 *
 * The construct is representable only under the squashed layering, because
 * under `srvn` each target lives in a submodel of its own and is replaced, in
 * the caller's submodel, by a surrogate delay: no node ever has arcs to more
 * than one of them. It further needs a layer solver that resolves the strategy
 * from the STATE, which in this port is `ssa` alone: `refresh_routing` expands
 * RROBIN and JSQ into a uniform probability split so that the matrix solvers
 * have a matrix to read, and returning that split under a round-robin label
 * would report a deterministic policy as a coin. Both refusals are asserted
 * here by name so they cannot regress into a silent answer.
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/solvers/ssa/ssa_dispatch.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** 0 = three independent calls of mean 1/3, 1 = round robin, 2 = JSQ. */
lqn::LqnBuilder<double> group_builder(int mode) {
    lqn::LqnBuilder<double> b;
    b.processor("PC", 1, SchedStrategy::INF);
    b.processor("PS", 1, SchedStrategy::PS);
    b.task("TC", 10, SchedStrategy::REF, "PC");
    b.think_time("TC", D::exp_rate(1.0 / 5));
    b.task("TS1", 5, SchedStrategy::FCFS, "PS");
    b.task("TS2", 5, SchedStrategy::FCFS, "PS");
    b.task("TS3", 5, SchedStrategy::FCFS, "PS");
    b.entry("EC", "TC");
    b.entry("ES1", "TS1");
    b.entry("ES2", "TS2");
    b.entry("ES3", "TS3");
    b.activity("AC", D::exp_rate(2), "TC");
    b.bound_to("AC", "EC");
    if (mode == 0) {
        b.sync_call("AC", "ES1", 1.0 / 3);
        b.sync_call("AC", "ES2", 1.0 / 3);
        b.sync_call("AC", "ES3", 1.0 / 3);
    } else if (mode == 1) {
        b.sync_call_round_robin("AC", {"ES1", "ES2", "ES3"}, 1.0);
    } else {
        b.sync_call_jsq("AC", {"ES1", "ES2", "ES3"}, 1.0);
    }
    b.activity("AS1", D::exp_rate(1), "TS1");
    b.bound_to("AS1", "ES1");
    b.replies_to("AS1", "ES1");
    b.activity("AS2", D::exp_rate(1), "TS2");
    b.bound_to("AS2", "ES2");
    b.replies_to("AS2", "ES2");
    b.activity("AS3", D::exp_rate(1), "TS3");
    b.bound_to("AS3", "ES3");
    b.replies_to("AS3", "ES3");
    return b;
}

lqn::LqnStruct<double> build_group(int mode) { return group_builder(mode).build(); }

std::string temp_path(const char* stem) {
    const char* base = std::getenv("TMPDIR");
    return (base && *base ? std::string(base) : std::string("/tmp")) + "/" + stem;
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

/** Router -> three identical FCFS queues, the dispatcher on its own. */
qn::Network<double> build_router(RoutingStrategy rs) {
    qn::Network<double> m("rr");
    const std::size_t s = m.add_source("Source");
    const std::size_t r = m.add_router("R");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, D::exp_rate(1.2));
    m.set_service(q1, c, D::exp_rate(0.6));
    m.set_service(q2, c, D::exp_rate(0.6));
    m.set_service(q3, c, D::exp_rate(0.6));
    qn::RoutingMatrix<double> P;
    P.set(c, c, s, r, 1.0);
    P.set(c, c, r, q1, 1.0 / 3);
    P.set(c, c, r, q2, 1.0 / 3);
    P.set(c, c, r, q3, 1.0 / 3);
    P.set(c, c, q1, k, 1.0);
    P.set(c, c, q2, k, 1.0);
    P.set(c, c, q3, k, 1.0);
    m.link(P);
    if (rs != RoutingStrategy::PROB) m.set_routing(r, c, rs);
    return m;
}

double total_qlen(const ssa::SsaSolution& s) {
    double t = 0.0;
    for (std::size_t i = 1; i <= 3; ++i) t += s.QN(i, 0);
    return t;
}

}  // namespace

TEST_CASE("state-dependent routing: a dispatcher is not a coin") {
    ssa::SsaOptions o;
    o.samples = 300000;
    o.seed = 23000;
    qn::Network<double> mp = build_router(RoutingStrategy::PROB);
    qn::Network<double> mr = build_router(RoutingStrategy::RROBIN);
    qn::Network<double> mj = build_router(RoutingStrategy::JSQ);
    const double qp = total_qlen(ssa::solver_ssa(mp.get_struct(), o));
    const double qr = total_qlen(ssa::solver_ssa(mr.get_struct(), o));
    const double qj = total_qlen(ssa::solver_ssa(mj.get_struct(), o));

    // Round robin removes the variance of the branching and JSQ additionally
    // steers away from the longest queue, so the three are strictly ordered.
    // A uniform expansion of the strategy would make all three equal, which is
    // exactly the defect this guards.
    CHECK(qj < qr);
    CHECK(qr < qp);
    // MATLAB SolverSSA on the same model reads 6.0407 / 4.2628 / 3.1834. The
    // two engines share no random stream, so these are compared loosely, in
    // expectation, and not as a reproduction of a path.
    CHECK(qp == doctest::Approx(6.04).epsilon(0.10));
    CHECK(qr == doctest::Approx(4.26).epsilon(0.10));
    CHECK(qj == doctest::Approx(3.18).epsilon(0.10));

    // JSQ over identical servers BALANCES them: the tie among equal queues is
    // split uniformly, and taking the first candidate instead skews the load
    // onto it (measured 1.72 / 1.47 / 1.23 before the split was uniform).
    const ssa::SsaSolution sj = ssa::solver_ssa(mj.get_struct(), o);
    const double lo = std::min(std::min(sj.QN(1, 0), sj.QN(2, 0)), sj.QN(3, 0));
    const double hi = std::max(std::max(sj.QN(1, 0), sj.QN(2, 0)), sj.QN(3, 0));
    CHECK(hi - lo < 0.10);
}

TEST_CASE("call groups: the builder records one dispatch, not n calls") {
    const lqn::LqnStruct<double> l = build_group(1);
    REQUIRE(l.callgroups.size() == 1);
    CHECK(l.callgroups[0].strategy == RoutingStrategy::RROBIN);
    CHECK(l.callgroups[0].targets.size() == 3);
    CHECK(l.callgroups[0].caller == idx_of(l, "A:AC"));
    // the members are ordinary sync calls of mean 1/n, so the aggregate call
    // rate matches the probabilistic twin's exactly
    CHECK(l.ncalls == 3);
    double total = 0.0;
    for (std::size_t c = 1; c <= l.ncalls; ++c) total += l.callproc_mean[c];
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));

    const lqn::LqnStruct<double> lj = build_group(2);
    CHECK(lj.callgroups[0].strategy == RoutingStrategy::JSQ);
    // and a group of one is not a dispatch decision
    lqn::LqnBuilder<double> b;
    b.processor("P", 1, SchedStrategy::PS);
    b.task("T", 1, SchedStrategy::REF, "P");
    b.entry("E", "T");
    b.activity("A", D::exp_rate(1), "T");
    b.bound_to("A", "E");
    CHECK_THROWS_AS(b.sync_call_round_robin("A", {"E"}, 1.0), InputError);
}

TEST_CASE("call groups: the layer carries a Router the strategy owns") {
    const lqn::LqnStruct<double> l = build_group(1);
    ln::LnOptions o;
    o.method = "flat";
    o.layer_solver = "ssa";
    o.layer_ssa.samples = 20000;
    o.iter_max = 2;
    ln::SolverLN<double> s(l, o);
    const qn::Layer<double>& L = s.layers()[0];

    // one Router, and the strategy sits on it for the dispatch class alone
    std::size_t router = 0, ndecl = 0;
    for (std::size_t i = 0; i < L.nodes.size(); ++i)
        if (L.nodes[i].nodetype == NodeType::Router) router = i + 1;
    REQUIRE(router != 0);
    std::size_t dispatch = 0;
    for (std::size_t r = 0; r < L.nodes[router - 1].routing.size(); ++r)
        if (L.nodes[router - 1].routing[r] == RoutingStrategy::RROBIN) {
            ++ndecl;
            dispatch = r + 1;
        }
    CHECK(ndecl == 1);

    // its arcs are exactly the three targets, and the hop does NOT switch class
    // -- a state-dependent routing function is evaluated at zero off the class
    // diagonal, so the switch goes on the return arc instead
    std::size_t arcs = 0;
    for (std::size_t j = 1; j <= L.nodes.size(); ++j)
        if (L.get_route(dispatch, dispatch, router, j) > 0.0) ++arcs;
    CHECK(arcs == 3);
}

TEST_CASE("call groups: the .lqnx dialect round-trips the dispatch") {
    // The member calls are ordinary synch-calls in the document, so the group
    // element carries the grouping ALONE: a reader that re-issued the calls
    // from it would double the call rate, and one that ignored it would report
    // a dispatcher as three independent coins.
    const lqn::LqnBuilder<double> b = group_builder(2);
    const std::string path = temp_path("line_callgroup_roundtrip.lqnx");
    lqn::write_lqnx(b.model(), path, "callgroup");

    const lqn::LqnStruct<double> before = build_group(2);
    const lqn::LqnStruct<double> after = lqn::read_lqnx<double>(path);
    REQUIRE(after.callgroups.size() == before.callgroups.size());
    REQUIRE(after.callgroups.size() == 1);
    CHECK(after.callgroups[0].strategy == before.callgroups[0].strategy);
    CHECK(after.callgroups[0].caller == before.callgroups[0].caller);
    REQUIRE(after.callgroups[0].targets.size() == before.callgroups[0].targets.size());
    for (std::size_t i = 0; i < before.callgroups[0].targets.size(); ++i)
        CHECK(after.names[after.callgroups[0].targets[i]] ==
              before.names[before.callgroups[0].targets[i]]);
    // the aggregate call rate is untouched by the trip
    CHECK(after.ncalls == before.ncalls);
    double total = 0.0;
    for (std::size_t c = 1; c <= after.ncalls; ++c) total += after.callproc_mean[c];
    CHECK(total == doctest::Approx(1.0).epsilon(1e-9));

    // a document with no group reads back with none, and RROBIN survives too
    lqn::write_lqnx(group_builder(0).model(), path, "plain");
    CHECK(lqn::read_lqnx<double>(path).callgroups.empty());
    lqn::write_lqnx(group_builder(1).model(), path, "rrobin");
    const lqn::LqnStruct<double> rr = lqn::read_lqnx<double>(path);
    REQUIRE(rr.callgroups.size() == 1);
    CHECK(rr.callgroups[0].strategy == RoutingStrategy::RROBIN);
    std::remove(path.c_str());
}

TEST_CASE("call groups: refused by name where they cannot be dispatched") {
    const lqn::LqnStruct<double> l = build_group(1);
    // under srvn the targets never share a submodel
    {
        ln::LnOptions o;
        o.method = "srvn.cs";
        o.layer_solver = "ssa";
        CHECK_THROWS_AS(ln::SolverLN<double>(l, o), UnsupportedError);
    }
    // and a matrix layer solver would answer with the uniform split
    {
        ln::LnOptions o;
        o.method = "flat";
        o.layer_solver = "mva";
        CHECK_THROWS_AS(ln::SolverLN<double>(l, o), UnsupportedError);
    }
    // a model with no group is unaffected by either gate
    {
        ln::LnOptions o;
        o.method = "srvn.cs";
        CHECK_NOTHROW(ln::SolverLN<double>(build_group(0), o));
    }
}
