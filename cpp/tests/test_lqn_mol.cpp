/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Tests of line/api/lqn/lqn_mol.h, the Method-of-Layers solve of an entry-only
 * LQN on the SRVN decomposition, with pfqn_qdamva as the layer solver.
 *
 * WHERE THE ORACLES COME FROM. Three kinds, and no number here was ever read
 * back out of this implementation.
 *
 *  1. EXACT LAWS. Identities that hold at ANY fixed point whatever the AMVA
 *     error -- the utilization law at an entry, flow balance across a
 *     synchronous call, a host utilization being the sum of the entry
 *     utilizations it carries and lying in [0,1], a task's mean busy threads
 *     never exceeding its thread count, and throughput never exceeding the
 *     bottleneck rate. These are true assertions, not tolerances.
 *
 *  2. THE MATLAB REFERENCE. matlab/src/api/lqn/lqn_mol.m is ground truth for
 *     this port (CLAUDE.md: MATLAB is the reference implementation). The
 *     tables below were produced by running it under R2026a on 2026-09-04 on
 *     exactly the models built here, printed at %.14g. The port reproduces all
 *     eight to 14 significant digits AND lands on the same iteration count, so
 *     the tolerance is 1e-12 relative rather than an approximation band: this
 *     asserts that the two implementations walk the same fixed point, not
 *     merely that they end up near each other.
 *
 *  3. A CLOSED FORM. At one customer nothing can queue, so the cycle time is
 *     the sum of the declared demands and the throughput is its reciprocal --
 *     computed by hand, independent of both implementations.
 *
 * NOT tested here: that lqn_mol is wired into any solver, because it is not.
 * It is an API-layer entry point, as in MATLAB.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/lqn/lqn_mol.h"
#include "line/lang/lqn/lqn_builder.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** Absolute element index of a hashname such as "E:E2" or "P:P1". */
std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

// ---------------------------------------------------------------------------
// the eight fixtures, mirroring the MATLAB run that produced the tables below
// ---------------------------------------------------------------------------

/** Two-tier client-server: T1 (REF, N threads, think Z) calls E2 on T2. */
lqn::LqnStruct<double> two_tier(double N, double Z, double d1, double d2, double y, double thr,
                                double c1, double c2) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", c1, SchedStrategy::PS);
    b.processor("P2", c2, SchedStrategy::PS);
    b.task("T1", N, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(Z));
    b.task("T2", thr, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(d1), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", y);
    b.activity("A2", D::exp_mean(d2), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    return b.build();
}

lqn::LqnStruct<double> m_molfix() { return two_tier(10, 5.0, 1.0, 0.8, 2.5, 2, 1, 1); }
lqn::LqnStruct<double> m_multisrv() { return two_tier(20, 1.0, 1.0, 0.8, 2.0, 4, 2, 3); }
lqn::LqnStruct<double> m_sat() { return two_tier(50, 5.0, 1.0, 0.8, 2.5, 5, 1, 1); }

/** matlab/examples/basic/layeredModel/lqn_basic.m, three tiers on multicore hosts. */
lqn::LqnStruct<double> m_basic() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 2, SchedStrategy::PS);
    b.processor("P2", 3, SchedStrategy::PS);
    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_rate(1.0 / 2.0));
    b.task("T2", 50, SchedStrategy::FCFS, "P1");
    b.think_time("T2", D::exp_rate(1.0 / 3.0));
    b.task("T3", 25, SchedStrategy::FCFS, "P2");
    b.think_time("T3", D::exp_rate(1.0 / 4.0));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", D::exp_rate(10.0), "T1");
    b.bound_to("AS1", "E1");
    b.replies_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1);
    b.activity("AS2", D::exp_rate(20.0), "T2");
    b.bound_to("AS2", "E2");
    b.replies_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5);
    b.activity("AS3", D::exp_rate(50.0), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b.build();
}

/** T1 -> T2 -> T3, so a callee is itself a caller. */
lqn::LqnStruct<double> m_chain3() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 20, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 3, SchedStrategy::FCFS, "P2");
    b.task("T3", 2, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.5), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 1);
    b.activity("A2", D::exp_mean(0.4), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.sync_call("A2", "E3", 2);
    b.activity("A3", D::exp_mean(0.3), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** One activity calling two different tasks, so a layer has two callees. */
lqn::LqnStruct<double> m_fanout() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 15, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(3.0));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.task("T3", 2, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.5), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 1);
    b.sync_call("A1", "E3", 2);
    b.activity("A2", D::exp_mean(0.6), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.4), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** Two entries on ONE task, which is what exercises the entry-share update. */
lqn::LqnStruct<double> m_multientry() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 12, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T2");
    b.activity("A1", D::exp_mean(0.5), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 1);
    b.sync_call("A1", "E3", 3);
    b.activity("A2", D::exp_mean(0.6), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.2), "T2");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** An INF task never queues for a thread, so the software layer is skipped. */
lqn::LqnStruct<double> m_inftask() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 10, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(4.0));
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(1.0), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 2);
    b.activity("A2", D::exp_mean(0.5), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    return b.build();
}

// ---------------------------------------------------------------------------
// the MATLAB reference tables
// ---------------------------------------------------------------------------

/** One reported row of `lqn_mol.m`, NaN written as a nan sentinel. */
struct Row {
    const char* hashname;
    double QN, UN, RN, TN;
};

struct Fixture {
    const char* name;
    lqn::LqnStruct<double> (*build)();
    std::size_t iter;  ///< iteration count the reference converged in
    std::vector<Row> rows;
};

const double NA = std::numeric_limits<double>::quiet_NaN();

/**
 * Produced by `lqn_mol.m` under MATLAB R2026a on 2026-09-04, printed at %.14g.
 * Hosts carry a utilization only; tasks carry no response time; entries carry
 * all four. An activity repeats its entry's row, which the mask test asserts.
 */
std::vector<Fixture> fixtures() {
    std::vector<Fixture> f;
    f.push_back({"molfix", m_molfix, 59, {
        {"P:P1", NA, 0.48793059982728, NA, NA},
        {"P:P2", NA, 0.97586119965455, NA, NA},
        {"R:T1", 7.5603470008636, 0.48793059982728, NA, 0.48793059982728},
        {"T:T2", 1.905723440048, 0.97586119965455, NA, 1.2198264995682},
        {"E:E1", 7.5603470008636, 0.48793059982728, 15.494717903612, 0.48793059982728},
        {"E:E2", 1.905723440048, 0.97586119965455, 1.5622905722434, 1.2198264995682}}});
    f.push_back({"lqn_basic", m_basic, 27, {
        {"P:P1", NA, 0.96744908657238, NA, NA},
        {"P:P2", NA, 0.42997737180995, NA, NA},
        {"R:T1", 24.201357691403, 0.64496605771492, NA, 12.899321154298},
        {"T:T2", 8.9270569129465, 0.32248302885746, NA, 12.899321154298},
        {"T:T3", 1.2899322001865, 0.42997737180995, NA, 64.496605771492},
        {"E:E1", 24.201357691403, 0.64496605771492, 1.8761729708031, 12.899321154298},
        {"E:E2", 8.9270569129465, 0.32248302885746, 0.69205633429568, 12.899321154298},
        {"E:E3", 1.2899322001865, 0.42997737180995, 0.020000001314127, 64.496605771492}}});
    f.push_back({"chain3", m_chain3, 154, {
        {"P:P1", NA, 0.76138077093716, NA, NA},
        {"P:P2", NA, 0.60910461674973, NA, NA},
        {"P:P3", NA, 0.9136569251246, NA, NA},
        {"R:T1", 16.954476916251, 0.76138077093716, NA, 1.5227615418743},
        {"T:T2", 2.940002951681, 0.60910461674973, NA, 1.5227615418743},
        {"T:T3", 1.6820774372842, 0.9136569251246, NA, 3.0455230837487},
        {"E:E1", 16.954476916251, 0.76138077093716, 11.134032775337, 1.5227615418743},
        {"E:E2", 2.940002951681, 0.60910461674973, 1.9307047563483, 1.5227615418743},
        {"E:E3", 1.6820774372842, 0.9136569251246, 0.55231150479865, 3.0455230837487}}});
    f.push_back({"fanout", m_fanout, 57, {
        {"P:P1", NA, 0.61014685674094, NA, NA},
        {"P:P2", NA, 0.73217622808912, NA, NA},
        {"P:P3", NA, 0.9762349707855, NA, NA},
        {"R:T1", 11.339118859554, 0.61014685674094, NA, 1.2202937134819},
        {"T:T2", 1.1550139345947, 0.73217622808912, NA, 1.2202937134819},
        {"T:T3", 1.9071499368363, 0.9762349707855, NA, 2.4405874269637},
        {"E:E1", 11.339118859554, 0.61014685674094, 9.2921226539801, 1.2202937134819},
        {"E:E2", 1.1550139345947, 0.73217622808912, 0.94650486340632, 1.2202937134819},
        {"E:E3", 1.9071499368363, 0.9762349707855, 0.78143069810408, 2.4405874269637}}});
    f.push_back({"multientry", m_multientry, 59, {
        {"P:P1", NA, 0.41255997925524, NA, NA},
        {"P:P2", NA, 0.99014395021258, NA, NA},
        {"R:T1", 10.349760082979, 0.41255997925524, NA, 0.82511995851048},
        {"T:T2", 1.9609664377822, 0.99014395021258, NA, 3.3004798340419},
        {"E:E1", 10.349760082979, 0.41255997925524, 12.543339881952, 0.82511995851048},
        {"E:E2", 0.9804832188911, 0.49507197510629, 1.1882917250736, 0.82511995851048},
        {"E:E3", 0.9804832188911, 0.49507197510629, 0.39609724169121, 2.4753598755315}}});
    f.push_back({"inftask", m_inftask, 52, {
        {"P:P1", NA, 0.83333383203266, NA, NA},
        {"P:P2", NA, 0.83333383203266, NA, NA},
        {"R:T1", 6.6666646718693, 0.83333383203266, NA, 0.83333383203266},
        {"T:T2", 3.3333422089391, 0.83333383203266, NA, 1.6666676640653},
        {"E:E1", 6.6666646718693, 0.83333383203266, 7.9999928187339, 0.83333383203266},
        {"E:E2", 3.3333422089391, 0.83333383203266, 2.0000041284826, 1.6666676640653}}});
    f.push_back({"multisrv", m_multisrv, 102, {
        {"P:P1", NA, 0.92679568568711, NA, NA},
        {"P:P2", NA, 0.98858206473291, NA, NA},
        {"R:T1", 18.146408628626, 0.92679568568711, NA, 1.8535913713742},
        {"T:T2", 3.8233345324695, 0.98858206473291, NA, 3.7071827427484},
        {"E:E1", 18.146408628626, 0.92679568568711, 9.7898646426976, 1.8535913713742},
        {"E:E2", 3.8233345324695, 0.98858206473291, 1.0313315522274, 3.7071827427484}}});
    f.push_back({"satN50", m_sat, 91, {
        {"P:P1", NA, 0.49985106774989, NA, NA},
        {"P:P2", NA, 0.99970213549978, NA, NA},
        {"R:T1", 47.500744661251, 0.49985106774989, NA, 0.49985106774989},
        {"T:T2", 4.9925402345185, 0.99970213549978, NA, 1.2496276693747},
        {"E:E1", 47.500744661251, 0.49985106774989, 95.029795324992, 0.49985106774989},
        {"E:E2", 4.9925402345185, 0.99970213549978, 3.9952222224853, 1.2496276693747}}});
    return f;
}

/** Compare against a reference cell, treating NaN as a value that must match. */
void check_cell(const std::string& where, double got, double want) {
    if (std::isnan(want)) {
        INFO(where << ": expected NaN, got " << got);
        CHECK(std::isnan(got));
    } else {
        INFO(where << ": expected " << want << ", got " << got);
        CHECK(got == doctest::Approx(want).epsilon(1e-12));
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// agreement with the MATLAB reference
// ---------------------------------------------------------------------------

TEST_CASE("lqnmol: every reported cell matches lqn_mol.m to 1e-12") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        for (std::size_t j = 0; j < fx[k].rows.size(); ++j) {
            const Row& row = fx[k].rows[j];
            const std::size_t i = idx_of(lsn, row.hashname);
            INFO("model " << fx[k].name << ", element " << row.hashname);
            REQUIRE(i != 0);
            const std::string w = std::string(fx[k].name) + "/" + row.hashname;
            check_cell(w + "/QN", r.QN[i], row.QN);
            check_cell(w + "/UN", r.UN[i], row.UN);
            check_cell(w + "/RN", r.RN[i], row.RN);
            check_cell(w + "/TN", r.TN[i], row.TN);
        }
    }
}

TEST_CASE("lqnmol: the fixed point is reached in the same number of sweeps as the reference") {
    // Two implementations of one iteration can agree on the answer while
    // disagreeing on the path to it. Matching the sweep count as well says the
    // relaxation, the residual and the layer solver all behave identically.
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(fx[k].build());
        INFO("model " << fx[k].name);
        CHECK(r.info.iter == fx[k].iter);
        CHECK(r.info.resid < 1e-6);   // converged
        CHECK(r.info.iter < 200);     // not merely exhausted
    }
}

TEST_CASE("lqnmol: an activity repeats the row of the entry it is bound to") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        for (std::size_t e = lsn.eshift + 1; e <= lsn.eshift + lsn.nentries; ++e) {
            const std::size_t a = lsn.actsof[e][0];
            INFO("model " << fx[k].name << ", entry " << lsn.hashnames[e]);
            CHECK(r.QN[a] == doctest::Approx(r.QN[e]).epsilon(1e-14));
            CHECK(r.UN[a] == doctest::Approx(r.UN[e]).epsilon(1e-14));
            CHECK(r.RN[a] == doctest::Approx(r.RN[e]).epsilon(1e-14));
            CHECK(r.TN[a] == doctest::Approx(r.TN[e]).epsilon(1e-14));
        }
    }
}

TEST_CASE("lqnmol: the reported columns carry the LQNS NaN mask") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        INFO("model " << fx[k].name);
        for (std::size_t h = 1; h <= lsn.nhosts; ++h) {
            // A processor has a utilization and nothing else, as in LQNS.
            CHECK(std::isnan(r.QN[h]));
            CHECK_FALSE(std::isnan(r.UN[h]));
            CHECK(std::isnan(r.RN[h]));
            CHECK(std::isnan(r.TN[h]));
        }
        for (std::size_t t = lsn.tshift + 1; t <= lsn.tshift + lsn.ntasks; ++t) {
            CHECK(std::isnan(r.RN[t]));  // a task has no response time
            CHECK_FALSE(std::isnan(r.TN[t]));
        }
        for (std::size_t e = lsn.eshift + 1; e <= lsn.eshift + lsn.nentries; ++e) {
            CHECK_FALSE(std::isnan(r.QN[e]));
            CHECK_FALSE(std::isnan(r.UN[e]));
            CHECK_FALSE(std::isnan(r.RN[e]));
            CHECK_FALSE(std::isnan(r.TN[e]));
        }
    }
}

// ---------------------------------------------------------------------------
// laws that hold at any fixed point
// ---------------------------------------------------------------------------

TEST_CASE("lqnmol: the utilization law holds at every entry") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        for (std::size_t e = lsn.eshift + 1; e <= lsn.eshift + lsn.nentries; ++e) {
            INFO("model " << fx[k].name << ", entry " << lsn.hashnames[e]);
            CHECK(r.QN[e] == doctest::Approx(r.TN[e] * r.RN[e]).epsilon(1e-12));
        }
    }
}

TEST_CASE("lqnmol: flow balance holds across every synchronous call") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        // The calling ACTIVITY's entry issues the call; on an entry-only model
        // that is the unique entry the activity is bound to.
        std::vector<std::size_t> entry_of_act(lsn.nidx + 1, 0);
        for (std::size_t e = lsn.eshift + 1; e <= lsn.eshift + lsn.nentries; ++e)
            entry_of_act[lsn.actsof[e][0]] = e;
        for (std::size_t c = 1; c <= lsn.ncalls; ++c) {
            const std::size_t src = entry_of_act[lsn.callpair_src[c]];
            const std::size_t dst = lsn.callpair_dst[c];
            INFO("model " << fx[k].name << ", call " << lsn.callhashnames[c]);
            CHECK(r.TN[dst] ==
                  doctest::Approx(lsn.callproc_mean[c] * r.TN[src]).epsilon(1e-9));
        }
    }
}

TEST_CASE("lqnmol: a host utilization is the sum of its entries' and lies in [0,1]") {
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        for (std::size_t h = 1; h <= lsn.nhosts; ++h) {
            double u = 0.0;
            for (std::size_t j = 0; j < lsn.tasksof[h].size(); ++j)
                for (std::size_t m = 0; m < lsn.entriesof[lsn.tasksof[h][j]].size(); ++m)
                    u += r.UN[lsn.entriesof[lsn.tasksof[h][j]][m]];
            INFO("model " << fx[k].name << ", host " << lsn.hashnames[h]);
            CHECK(r.UN[h] == doctest::Approx(u).epsilon(1e-12));
            // LINE scales a station utilization into [0,1] whatever its
            // multiplicity, so a multicore host is normalized, not summed.
            CHECK(r.UN[h] <= 1.0 + 1e-9);
            CHECK(r.UN[h] >= 0.0);
        }
    }
}

TEST_CASE("lqnmol: a task never holds more busy threads than it has") {
    // QN at a task is sum(X*S) over its entries, i.e. mean busy threads by
    // Little's law. It cannot exceed the thread count whatever the AMVA says.
    const std::vector<Fixture> fx = fixtures();
    for (std::size_t k = 0; k < fx.size(); ++k) {
        const lqn::LqnStruct<double> lsn = fx[k].build();
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        for (std::size_t t = lsn.tshift + 1; t <= lsn.tshift + lsn.ntasks; ++t) {
            if (lsn.isref[t] || lsn.sched[t] == SchedStrategy::INF) continue;
            INFO("model " << fx[k].name << ", task " << lsn.hashnames[t]);
            CHECK(r.QN[t] <= lsn.mult[t] + 1e-9);
        }
    }
}

// ---------------------------------------------------------------------------
// closed form and asymptotics
// ---------------------------------------------------------------------------

TEST_CASE("lqnmol: one customer cannot queue, so the cycle is the declared demand") {
    // N=1: the single thread of T1 is the only job in the model, so no station
    // can ever hold two. Cycle = Z + d1 + y*d2 = 5 + 1 + 2.5*0.8 = 8 exactly.
    const lqn::LqnStruct<double> lsn = two_tier(1, 5.0, 1.0, 0.8, 2.5, 1, 1, 1);
    const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
    CHECK(r.RN[idx_of(lsn, "E:E2")] == doctest::Approx(0.8).epsilon(1e-12));
    CHECK(r.RN[idx_of(lsn, "E:E1")] == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(r.TN[idx_of(lsn, "R:T1")] == doctest::Approx(1.0 / 8.0).epsilon(1e-12));
    CHECK(r.TN[idx_of(lsn, "E:E2")] == doctest::Approx(2.5 / 8.0).epsilon(1e-12));
    CHECK(r.UN[idx_of(lsn, "P:P2")] == doctest::Approx(2.5 * 0.8 / 8.0).epsilon(1e-12));
}

TEST_CASE("lqnmol: throughput rises with population towards the bottleneck rate") {
    // P2 receives 2.5 calls of 0.8s per reference cycle, so it needs 2.0s of
    // service per cycle against P1's 1.0s and T2's 2 threads. The processor
    // binds first, capping X(T1) at 1/2.0 = 0.5 for any population.
    const double bound = 0.5;
    const double N[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000};
    double prev = 0.0;
    for (std::size_t k = 0; k < sizeof(N) / sizeof(N[0]); ++k) {
        const lqn::LqnStruct<double> lsn = two_tier(N[k], 5.0, 1.0, 0.8, 2.5, 2, 1, 1);
        const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
        const double X = r.TN[idx_of(lsn, "R:T1")];
        INFO("N = " << N[k] << ", X = " << X);
        CHECK(r.info.iter < 200);
        CHECK(X <= bound + 1e-9);   // never exceeds the bottleneck rate
        CHECK(X > prev);            // strictly increasing in the population
        prev = X;
    }
    // and it actually gets there rather than saturating early or low
    CHECK(prev == doctest::Approx(bound).epsilon(1e-5));
}

TEST_CASE("lqnmol: a single-thread server saturates at its thread, not its processor") {
    // T2 holds one thread of 0.8s service on a 4-core host, so the thread is
    // the constraint: 1/0.8 = 1.25 calls/s, and X(T1) = 1.25/2.5 = 0.5.
    const lqn::LqnStruct<double> lsn = two_tier(500, 5.0, 1.0, 0.8, 2.5, 1, 1, 4);
    const lqn::LqnMolResult<double> r = lqn::lqn_mol(lsn);
    CHECK(r.TN[idx_of(lsn, "R:T1")] == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(r.TN[idx_of(lsn, "R:T1")] <= 0.5 + 1e-9);
    // the one thread is the thing that is full, and it is never over-full
    CHECK(r.QN[idx_of(lsn, "T:T2")] <= 1.0 + 1e-9);
    CHECK(r.QN[idx_of(lsn, "T:T2")] == doctest::Approx(1.0).epsilon(1e-4));
}

// ---------------------------------------------------------------------------
// refusals: every out-of-scope feature is named, not silently approximated
// ---------------------------------------------------------------------------

TEST_CASE("lqnmol: an asynchronous call is refused") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(1.0), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.async_call("A1", "E2", 1);
    b.activity("A2", D::exp_mean(0.5), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    CHECK_THROWS_AS(lqn::lqn_mol(b.build()), UnsupportedError);
}

TEST_CASE("lqnmol: an open arrival is refused") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.open_arrival("E2", D::exp_rate(0.2));
    b.activity("A1", D::exp_mean(1.0), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 1);
    b.activity("A2", D::exp_mean(0.5), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    CHECK_THROWS_AS(lqn::lqn_mol(b.build()), UnsupportedError);
}

TEST_CASE("lqnmol: a scheduling discipline the decomposition cannot represent is refused") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::SIRO);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(1.0), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    b.sync_call("A1", "E2", 1);
    b.activity("A2", D::exp_mean(0.5), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    CHECK_THROWS_AS(lqn::lqn_mol(b.build()), UnsupportedError);
}

TEST_CASE("lqnmol: an activity graph is refused") {
    // lqn_serial.m: T1 runs AS1 then AS2 in series, so E1 binds a graph.
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 10, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(100.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("AS1", D::exp_mean(1.6), "T1");
    b.bound_to("AS1", "E1");
    b.activity("AS2", D::immediate(), "T1");
    b.sync_call("AS2", "E2", 1);
    b.replies_to("AS2", "E1");
    b.activity("AS3", D::exp_mean(5.0), "T2");
    b.bound_to("AS3", "E2");
    b.activity("AS4", D::exp_mean(1.0), "T2");
    b.replies_to("AS4", "E2");
    b.serial("AS1", "AS2");
    b.serial("AS3", "AS4");
    CHECK_THROWS_AS(lqn::lqn_mol(b.build()), UnsupportedError);
}
