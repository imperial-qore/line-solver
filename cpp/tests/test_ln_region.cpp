/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * An admission constraint on a task: n(E2) + n(E3) <= cap.
 *
 * Mirrors line-test.git/test/testsAdvFeatures/fcr/test_fcr_lincon_ln.m. There
 * is no external oracle for this shape -- JMT's FCR XML cannot express a
 * general A n <= b and lqns has no admission-constraint concept -- so the
 * assertions are MATLAB SolverLN plus the structural and conservation
 * properties the reference test checks.
 *
 * The constrained layer CANNOT go through MVA, which refuses a Region by name,
 * so SolverLN must route it to CTMC on its own (SolverLN.m:1047-1089).
 *
 * Reference numbers are MATLAB SolverLN.getAvgTable() (defaultOptions,
 * 2026-07-29), where layer T:T2 reports nregions=1 and solver=SolverCTMC.
 */

#include <limits>
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

lqn::LqnStruct<double> build_region(double cap) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);

    b.task("T1", 2, SchedStrategy::REF, "P1");
    b.think_time("T1", E(1.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");

    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T2");

    b.activity("A1", E(0.4), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 1.0);

    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");

    b.activity("A3", E(1.0), "T2");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");

    // n(E2) + n(E3) <= cap, a passive resource of T2 as a whole
    b.add_constraint("T2", {"E2", "E3"}, {1.0, 1.0}, cap);
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("region: the constraint reaches the struct in ELEMENT columns") {
    const lqn::LqnStruct<double> l = build_region(1.0);
    const std::size_t t2 = idx_of(l, "T:T2");
    REQUIRE(t2 > 0);
    REQUIRE(t2 < l.lincon_A.size());
    // one row, one column per ENTRY of T2 -- not per class
    CHECK(l.lincon_A[t2].rows() == 1);
    CHECK(l.lincon_A[t2].cols() == l.entriesof[t2].size());
    CHECK(l.lincon_A[t2](0, 0) == doctest::Approx(1.0));
    CHECK(l.lincon_A[t2](0, 1) == doctest::Approx(1.0));
    REQUIRE(l.lincon_b[t2].size() == 1);
    CHECK(l.lincon_b[t2][0] == doctest::Approx(1.0));
    // a host with no constraint stays empty
    CHECK(l.lincon_A[idx_of(l, "P:P1")].rows() == 0);
}

TEST_CASE("region: exactly one layer carries it, and it is not solved by MVA") {
    const lqn::LqnStruct<double> l = build_region(1.0);
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);

    std::size_t nregionlayers = 0;
    for (const qn::Layer<double>& L : s.layers()) {
        CAPTURE(L.name);
        if (!L.regions.empty()) {
            ++nregionlayers;
            CHECK(L.name == "T:T2");
            REQUIRE(L.regions.size() == 1);
            // the region spans the server station, and its columns are CLASSES
            CHECK(L.regions[0].lincon_A.cols() == L.nclasses);
            CHECK(L.regions[0].members[L.serverIdx - 1]);
            // an entry column expands to the CALL classes that target it, so
            // the row must be non-empty after the translation
            double rowsum = 0.0;
            for (std::size_t k = 0; k < L.regions[0].lincon_A.cols(); ++k)
                rowsum += L.regions[0].lincon_A(0, k);
            CHECK(rowsum == doctest::Approx(2.0));
        }
    }
    CHECK(nregionlayers == 1);

    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);
    // MVA refuses a Region by name, so reaching an answer at all proves the
    // CTMC dispatch fired.
    CHECK(sol.TN[idx_of(l, "E:E2")] > 0.0);
}

TEST_CASE("region: the constrained model matches MATLAB SolverLN") {
    const lqn::LqnStruct<double> l = build_region(1.0);
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // MATLAB AvgTable: P1 Util 0.17481, P2 Util 0.87404,
    //   T1 QLen 1.563 Tput 0.43702; T2 QLen 0.87404 Tput 0.87404
    //   E1 RespT 3.5765; E2/E3 RespT 1 Tput 0.43702; A1 RespT 2.4
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:P1", 0.17481, 0.0, 0.0},
        {"P:P2", 0.87404, 0.0, 0.0},
        {"R:T1", 0.17481, 0.43702, 0.0},
        {"T:T2", 0.87404, 0.87404, 0.0},
        {"E:E1", 0.17481, 0.43702, 3.5765},
        {"E:E2", 0.43702, 0.43702, 1.0},
        {"E:E3", 0.43702, 0.43702, 1.0},
        {"A:A1", 0.17481, 0.43702, 2.4},
        {"A:A2", 0.43702, 0.43702, 1.0},
        {"A:A3", 0.43702, 0.43702, 1.0},
    };
    for (const Row& row : rows) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(std::string(row.hn));
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
        if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
    }

    // The constraint must actually bind: n(E2) + n(E3) <= 1.
    CHECK(sol.QN[idx_of(l, "E:E2")] + sol.QN[idx_of(l, "E:E3")] <= 1.0 + 1e-6);
    // Flow balance across the layer boundary: A1 calls each entry once.
    CHECK(sol.TN[idx_of(l, "E:E2")] ==
          doctest::Approx(sol.TN[idx_of(l, "A:A1")]).epsilon(1e-6));
}

TEST_CASE("region: throughput is monotone in the cap") {
    double prev = -1.0;
    for (double cap = 1.0; cap <= 2.0; cap += 1.0) {
        const lqn::LqnStruct<double> l = build_region(cap);
        ln::LnOptions opt;
        ln::SolverLN<double> s(l, opt);
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        const double x = sol.TN[idx_of(l, "E:E2")];
        CAPTURE(cap);
        CHECK(x > 0.0);
        CHECK(x >= prev - 1e-9);  // relaxing the constraint cannot slow the model
        prev = x;
    }
}
