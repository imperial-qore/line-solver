/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * A synchronous call made a non-unit mean number of times per visit.
 *
 * `A2` calls entry `E3` five times for every one of its own executions. In the
 * layer that holds E3's server the caller therefore has to enter the call class
 * five times, and since a self-loop on that class would re-enter its service as
 * well, the reference manufactures the visit count with a second class,
 * `<call>.Aux`, that carries no service and exists only as a routing waypoint.
 * See route_sync_call in solver_ln.h for the three-way construction and why the
 * direction of the inequality decides which class carries the time.
 *
 * The observable that pins it down is flow conservation: E3 must complete
 * exactly five times for every completion of E2. A wrong number of visits is
 * not a small numerical error, it is a throughput ratio that is not five.
 *
 * The model is lqn_basic with the multiplicities dropped to one, so the layers
 * are single-server and the multiserver AMVA approximation does not enter; the
 * reference values are MATLAB SolverLN(SolverMVA) on exactly this model.
 */

#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** lqn_basic, single-server, with the five-fold call from AS2 to E3. */
lqn::LqnStruct<double> build_multicall(double calls_to_e3) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P1");
    b.task("T3", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", D::exp_mean(0.1), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", D::exp_mean(0.05), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", calls_to_e3);
    b.replies_to("AS2", "E2");
    b.activity("AS3", D::exp_mean(0.02), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("multi-call: the callee is entered calls-mean times per caller completion") {
    for (double cm : {1.0, 2.0, 5.0}) {
        CAPTURE(cm);
        const lqn::LqnStruct<double> l = build_multicall(cm);
        ln::SolverLN<double> s(l, ln::LnOptions());
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        const std::size_t e2 = idx_of(l, "E:E2"), e3 = idx_of(l, "E:E3");
        REQUIRE(e2 > 0);
        REQUIRE(e3 > 0);
        REQUIRE(sol.TN[e2] > 0.0);
        // the defining property, independent of any reference table
        CHECK(sol.TN[e3] / sol.TN[e2] == doctest::Approx(cm).epsilon(1e-5));
    }
}

TEST_CASE("multi-call: the whole AvgTable against MATLAB SolverLN(SolverMVA)") {
    const lqn::LqnStruct<double> l = build_multicall(5.0);
    ln::LnOptions opt;
    // the encoding this table was recorded under; the 'srvn' alias now resolves
    // to 'srvn.ph' here, which reaches its own fixed point (T1 QLen 40, E1 RespT
    // 8, T3 Tput 25 in MATLAB) rather than the routing one below
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    struct Row {
        const char* hn;
        double qlen, util, respt, residt, tput;
    };
    // MATLAB getAvgTable, this model, SolverLN(@(m) SolverMVA(m)). RE-RECORDED
    // 2026-08-11, after the interlock probability was aligned to Li and Franks
    // (2015), Eq. (5), and to lqns' m' rule.
    const Row rows[] = {
        {"P:P1", 0.0, 0.757549508549238, 0.0, 0.0, 0.0},
        {"P:P2", 0.0, 0.499049051376959, 0.0, 0.0, 0.0},
        {"R:T1", 39.8, 0.508024858475075, 0.0, 0.289670675354027, 5.08024858475075},
        {"T:T2", 0.997395328303713, 0.249524650074163, 0.0, 0.0998590776517012,
         4.99049300148326},
        {"T:T3", 0.499049051376959, 0.499049051376959, 0.0, 0.02, 24.9524525688479},
        {"E:E1", 39.8, 0.508024858475075, 7.84012220254295, 0.0, 5.08024858475075},
        {"E:E2", 0.997395328303713, 0.249524650074163, 0.2, 0.0, 4.99049300148326},
        {"E:E3", 0.499049051376959, 0.499049051376959, 0.02, 0.0, 24.9524525688479},
        {"A:AS1", 40.6, 0.508024858475075, 8.0, 0.289670675354027, 5.08024858475075},
        {"A:AS2", 1.0, 0.249524650074163, 0.200362704739934, 0.0998590776517012,
         4.99049300148326},
        {"A:AS3", 0.499049051376959, 0.499049051376959, 0.02, 0.02, 24.9524525688479},
    };
    // CoarseTol, not the arithmetic's precision, for two reasons that both floor
    // the achievable agreement at 1e-3. The reference numbers are the table
    // getAvgTable PRINTS, and it snaps a value to one decimal place whenever it
    // is already within CoarseTol of one -- 39.9232 is reported as 39.9, 0.999909
    // as 1 -- so several rows carry no more than three digits of information. On
    // top of that MATLAB and this port stop at slightly different iterates of the
    // same Picard sequence. The sharp assertion in this file is the throughput
    // ratio above, which is exact by construction and needs no reference table.
    for (const Row& r : rows) {
        const std::size_t i = idx_of(l, r.hn);
        CAPTURE(r.hn);
        REQUIRE(i > 0);
        if (r.qlen > 1e-7) CHECK(sol.QN[i] == doctest::Approx(r.qlen).epsilon(1e-3));
        if (r.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(r.util).epsilon(1e-3));
        if (r.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(r.respt).epsilon(1e-3));
        if (r.residt > 1e-7) CHECK(sol.WN[i] == doctest::Approx(r.residt).epsilon(1e-3));
        if (r.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(r.tput).epsilon(1e-3));
    }
}
