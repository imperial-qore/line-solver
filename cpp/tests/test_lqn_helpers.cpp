/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Unit tests of line/solvers/ln/lqn_helpers.h: the forwarding rewrite and the
 * LQNS overtaking chain.
 *
 * WHY THE ORACLES ARE WHAT THEY ARE. Neither routine has a published table of
 * numbers to check against at this level, so every expected value below is
 * either an identity the model must satisfy whatever the implementation, or a
 * hand substitution into the LQNS rate equations quoted in the source. Nothing
 * is a number this implementation produced.
 *
 * For the forwarding rewrite the identities are structural: a chain of length
 * one with probability one must leave the model indistinguishable from a direct
 * rendezvous to the forwarded-to entry, a chain of length two must carry the
 * PRODUCT of its probabilities, a target the caller already calls directly must
 * gain visits rather than a second arc, and a chain that loops back to the
 * caller's own task must produce nothing at all.
 *
 * For the overtaking chain they are limits and invariances: a server with no
 * work in the tested phase cannot be found in it, a client that never calls the
 * server cannot overtake at it, and every rate in the chain is a ratio of times,
 * so the answer cannot move when all times are rescaled together.
 *
 * The LqnBuilder has no forwarding API (the .lqnx reader is the only producer of
 * FWD calls), so the tests append the FWD calls to the built struct exactly as
 * lqn_reader.h does.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/lqn_helpers.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** Absolute element index of a hashname such as "E:E2" or "A:A1". */
std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

/** Append a forwarding call, mirroring the FWD branch of lqn_reader.h. */
void add_fwd(lqn::LqnStruct<double>& l, std::size_t src_e, std::size_t dst_e, double prob) {
    ++l.ncalls;
    l.callpair_src.push_back(src_e);
    l.callpair_dst.push_back(dst_e);
    l.calltype.push_back(CallType::FWD);
    l.callproc_mean.push_back(prob);
    l.callnames.push_back(l.names[src_e] + "~>" + l.names[dst_e]);
    l.callhashnames.push_back(l.hashnames[src_e] + "~>" + l.hashnames[dst_e]);
    l.taskgraph.set(l.parent[src_e], l.parent[dst_e], 1.0);
    l.graph.set(src_e, dst_e, 1.0);
}

/**
 * A chain of `ntasks` single-entry tasks, each on its own processor, with T1 the
 * reference task and its activity A1 making one synchronous call to E2.
 * Forwarding arcs are added by the tests on top of this.
 */
lqn::LqnStruct<double> chain_model(std::size_t ntasks, double base_mean) {
    lqn::LqnBuilder<double> b;
    for (std::size_t i = 1; i <= ntasks; ++i) {
        const std::string s = std::to_string(i);
        b.processor("P" + s, 1, SchedStrategy::PS);
        b.task("T" + s, i == 1 ? 3 : 1, i == 1 ? SchedStrategy::REF : SchedStrategy::FCFS,
               "P" + s);
    }
    for (std::size_t i = 1; i <= ntasks; ++i)
        b.entry("E" + std::to_string(i), "T" + std::to_string(i));
    for (std::size_t i = 1; i <= ntasks; ++i) {
        const std::string s = std::to_string(i);
        b.activity("A" + s, D::exp_mean(0.1 * double(i)), "T" + s);
        b.bound_to("A" + s, "E" + s);
        b.replies_to("A" + s, "E" + s);
    }
    b.sync_call("A1", "E2", base_mean);
    return b.build();
}

/** Mean of the SYNC call from `src` to `dst`, or -1 when there is none. */
double sync_mean(const lqn::LqnStruct<double>& l, std::size_t src, std::size_t dst) {
    for (std::size_t c = 1; c <= l.ncalls; ++c)
        if (l.calltype[c] == CallType::SYNC && l.callpair_src[c] == src &&
            l.callpair_dst[c] == dst)
            return l.callproc_mean[c];
    return -1.0;
}

}  // namespace

TEST_CASE("lqnfwd: a model without forwarding is left untouched") {
    lqn::LqnStruct<double> l = chain_model(3, 3.0);
    const std::size_t ncalls0 = l.ncalls;
    const std::vector<double> means0 = l.callproc_mean;

    ln::lqn_fwd_rendezvous(l);

    CHECK(l.ncalls == ncalls0);
    CHECK(l.callproc_mean == means0);
}

TEST_CASE("lqnfwd: a forwarding chain of length one is an ordinary rendezvous") {
    lqn::LqnStruct<double> l = chain_model(3, 3.0);
    const std::size_t a1 = idx_of(l, "A:A1"), e3 = idx_of(l, "E:E3");
    // "R:T1", not "T:T1": chain_model declares task 1 SchedStrategy::REF, and a
    // reference task is hashed with the R prefix (getStruct.m:129 against :132
    // for every other task). idx_of returns 0 on a miss, so the wrong prefix
    // does not fail the lookup -- it silently asks about element 0.
    const std::size_t t1 = idx_of(l, "R:T1"), t3 = idx_of(l, "T:T3");
    add_fwd(l, idx_of(l, "E:E2"), e3, 1.0);
    const std::size_t ncalls0 = l.ncalls;

    ln::lqn_fwd_rendezvous(l);

    // E2 forwards everything it is sent to E3, so A1 blocks on E3 once per call
    // it makes to E2: the same arc a direct rendezvous would have produced.
    REQUIRE(l.ncalls == ncalls0 + 1);
    CHECK(sync_mean(l, a1, e3) == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(l.calltype[l.ncalls] == CallType::SYNC);
    CHECK(l.callsof[a1].back() == l.ncalls);

    // and the caller must now be visible as a client of T3 everywhere the layer
    // builder looks, or T3's layer is built without it
    CHECK(l.iscaller.get(t1, t3));
    CHECK(l.iscaller.get(a1, e3));
    CHECK(l.issynccaller.get(t1, t3));
    CHECK(l.issynccaller.get(a1, e3));
    CHECK(l.graph.get(a1, e3) != 0.0);
    CHECK(l.taskgraph.get(t1, t3) != 0.0);
}

TEST_CASE("lqnfwd: a chain of length two carries the product of its probabilities") {
    lqn::LqnStruct<double> l = chain_model(4, 3.0);
    const std::size_t a1 = idx_of(l, "A:A1");
    const std::size_t e2 = idx_of(l, "E:E2"), e3 = idx_of(l, "E:E3"), e4 = idx_of(l, "E:E4");
    add_fwd(l, e2, e3, 0.5);
    add_fwd(l, e3, e4, 0.4);
    const std::size_t ncalls0 = l.ncalls;

    ln::lqn_fwd_rendezvous(l);

    // half the calls reach E3, and 40% of those go on to E4
    REQUIRE(l.ncalls == ncalls0 + 2);
    CHECK(sync_mean(l, a1, e3) == doctest::Approx(3.0 * 0.5).epsilon(1e-12));
    CHECK(sync_mean(l, a1, e4) == doctest::Approx(3.0 * 0.5 * 0.4).epsilon(1e-12));
    CHECK(sync_mean(l, a1, e2) == doctest::Approx(3.0).epsilon(1e-12));
}

TEST_CASE("lqnfwd: forwarded work joins an arc the caller already has") {
    lqn::LqnBuilder<double> b;
    for (std::size_t i = 1; i <= 3; ++i) {
        const std::string s = std::to_string(i);
        b.processor("P" + s, 1, SchedStrategy::PS);
        b.task("T" + s, i == 1 ? 3 : 1, i == 1 ? SchedStrategy::REF : SchedStrategy::FCFS,
               "P" + s);
        b.entry("E" + s, "T" + s);
    }
    for (std::size_t i = 1; i <= 3; ++i) {
        const std::string s = std::to_string(i);
        b.activity("A" + s, D::exp_mean(0.1), "T" + s);
        b.bound_to("A" + s, "E" + s);
        b.replies_to("A" + s, "E" + s);
    }
    b.sync_call("A1", "E2", 3.0);
    b.sync_call("A1", "E3", 2.0);
    lqn::LqnStruct<double> l = b.build();
    const std::size_t a1 = idx_of(l, "A:A1"), e3 = idx_of(l, "E:E3");
    add_fwd(l, idx_of(l, "E:E2"), e3, 0.5);
    const std::size_t ncalls0 = l.ncalls;

    ln::lqn_fwd_rendezvous(l);

    // A second arc between the same pair would be a second class in T3's layer
    // competing with itself; the visits belong on the arc that is already there.
    CHECK(l.ncalls == ncalls0);
    CHECK(sync_mean(l, a1, e3) == doctest::Approx(2.0 + 3.0 * 0.5).epsilon(1e-12));
}

TEST_CASE("lqnfwd: a chain that returns to the caller's own task yields no arc") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 3, SchedStrategy::REF, "P1");
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E1b", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.1), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.activity("A1b", D::exp_mean(0.1), "T1");
    b.bound_to("A1b", "E1b");
    b.replies_to("A1b", "E1b");
    b.activity("A2", D::exp_mean(0.1), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.1), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    b.sync_call("A1", "E2", 3.0);
    lqn::LqnStruct<double> l = b.build();
    const std::size_t a1 = idx_of(l, "A:A1"), e1b = idx_of(l, "E:E1b");
    add_fwd(l, idx_of(l, "E:E2"), idx_of(l, "E:E3"), 1.0);
    add_fwd(l, idx_of(l, "E:E3"), e1b, 1.0);
    const std::size_t ncalls0 = l.ncalls;

    ln::lqn_fwd_rendezvous(l);

    // T1 would otherwise be entered as its own client, which is not a queueing
    // relation the layer for T1 can express.
    CHECK(l.ncalls == ncalls0 + 1);  // only the arc to E3
    CHECK(sync_mean(l, a1, e1b) == -1.0);
}

// ---------------------------------------------------------------------------
// lqn_overtake_markov
// ---------------------------------------------------------------------------

TEST_CASE("lqnovt: a server with no residence in the tested phase is never overtaken") {
    // rows p=0..2: [nSlices service y_ij y_ik t_k]
    Matrix<double> cp({{1.0, 0.5, 0.0, 0.0, 0.0},
                       {4.0, 1.2, 2.0, 1.0, 0.7},
                       {2.0, 0.8, 1.0, 0.0, 0.0}});
    const std::vector<double> y_aj{3.0, 2.0, 1.0};

    // xj is the residence of the server phase being tested: with none, there is
    // no interval in which an arrival can find the server there.
    CHECK(ln::lqn_overtake_markov(cp, 0.8, 0.0, y_aj) == doctest::Approx(0.0).epsilon(1e-14));
}

TEST_CASE("lqnovt: a client that never calls the server cannot overtake at it") {
    Matrix<double> cp({{1.0, 0.5, 0.0, 0.0, 0.0},
                       {3.0, 1.2, 0.0, 2.0, 0.7},
                       {2.0, 0.8, 0.0, 1.0, 0.4}});
    const std::vector<double> y_aj{0.0, 0.0, 0.0};

    CHECK(ln::lqn_overtake_markov(cp, 0.8, 0.6, y_aj) == doctest::Approx(0.0).epsilon(1e-14));
}

TEST_CASE("lqnovt: the probability is invariant to the unit of time") {
    Matrix<double> cp({{1.0, 0.5, 0.0, 0.0, 0.0},
                       {4.0, 1.2, 2.0, 1.0, 0.7},
                       {2.0, 0.8, 1.0, 0.0, 0.0}});
    const std::vector<double> y_aj{3.0, 2.0, 1.0};
    const double prVisit = 0.8, xj = 0.6;
    const double base = ln::lqn_overtake_markov(cp, prVisit, xj, y_aj);

    // Every rate in the chain is xj/(xj+something) or a ratio of call counts,
    // so measuring all times in units seven times smaller cannot move a
    // probability. Anything that breaks this has an absolute time somewhere.
    const double k = 7.0;
    Matrix<double> cpk = cp;
    for (std::size_t p = 0; p < cpk.rows(); ++p) {
        cpk(p, 1) *= k;  // host residence of the phase
        cpk(p, 4) *= k;  // mean time at the other servers
    }

    CHECK(ln::lqn_overtake_markov(cpk, prVisit, xj * k, y_aj) ==
          doctest::Approx(base).epsilon(1e-12));
    CHECK(base > 0.0);
}

TEST_CASE("lqnovt: the single-phase client matches the hand-solved chain") {
    // One client phase (plus the think slice), one call to the server task and
    // none to any other, so setRates collapses to
    //   b0 = xj/(xj+Z),  q1 = xj/(xj+S/(1+y)),
    //   b1 = q1*prVisit/(y+1),  c1 = q1*y/(y+1),  d0 = d1 = 0,
    // and the two-state chain leaves prOt = c1/(1 - b1*b0), the y_aj ratio being
    // one because the single phase issues every call.
    const double Z = 3.0, S = 4.0, y = 1.0, xj = 2.0, prVisit = 0.5;
    Matrix<double> cp({{1.0, Z, 0.0, 0.0, 0.0}, {1.0 + y, S, y, 0.0, 0.0}});
    const std::vector<double> y_aj{y, y};

    const double b0 = xj / (xj + Z);
    const double q1 = xj / (xj + S / (1.0 + y));
    const double b1 = q1 * prVisit / (y + 1.0);
    const double c1 = q1 * y / (y + 1.0);
    const double expected = c1 / (1.0 - b1 * b0);

    CHECK(ln::lqn_overtake_markov(cp, prVisit, xj, y_aj) ==
          doctest::Approx(expected).epsilon(1e-12));
}
