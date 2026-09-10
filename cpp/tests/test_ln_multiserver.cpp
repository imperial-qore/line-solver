/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * A layered model whose stations have more than one server.
 *
 * This is lqn_basic: two PS processors of multiplicity 2 and 3, three tasks of
 * multiplicity 50, 50 and 25, think times on the non-reference tasks, and a
 * five-fold call. Every station in it is a multiserver, which is what the
 * earlier models -- lqn_ofbiz, lqn_workflows, lqn_serial -- did not exercise:
 * their layers are single-server or infinite-server throughout, so the whole
 * multiserver branch of solver_amvald went unmeasured until this model.
 *
 * Two things it pins down, both invisible on a single-server layer:
 *
 *  1. The soft-minimum multiserver term takes `1 + interpTotArvlQlen + mean(g)`,
 *     where g is the Linearizer's own correction. gamma is zero at a
 *     single-server station, so dropping the term costs nothing there and
 *     costs 40% here.
 *  2. solver_amvald warm-starts its queue lengths from options.init_sol. A
 *     saturated layer stops on the iteration cap rather than on its tolerance,
 *     so the iterate it returns depends on where it started; ignoring the warm
 *     start agrees at the first outer iteration and diverges from the second.
 *
 * Reference: MATLAB SolverLN(@(m) SolverMVA(m)) on this model.
 */

#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

lqn::LqnStruct<double> build_basic() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 2, SchedStrategy::PS);
    b.processor("P2", 3, SchedStrategy::PS);
    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 50, SchedStrategy::FCFS, "P1");
    b.think_time("T2", D::exp_mean(3.0));
    b.task("T3", 25, SchedStrategy::FCFS, "P2");
    b.think_time("T3", D::exp_mean(4.0));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", D::exp_mean(0.1), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", D::exp_mean(0.05), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5.0);
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

TEST_CASE("lqn_basic: multiserver layers against MATLAB SolverLN(SolverMVA)") {
    const lqn::LqnStruct<double> l = build_basic();
    ln::LnOptions opt;
    // the encoding this table was recorded under; the 'srvn' alias now resolves
    // to 'srvn.ph' here, and the case below asserts that ensemble separately
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    struct Row {
        const char* hn;
        double qlen, util, respt, residt, tput;
    };
    const Row rows[] = {
        // Re-recorded 2026-08-14 against MATLAB LN(@SolverMVA,'method','srvn.cs'),
        // after 58eb739f1 took the arriving chain out of the AMVA arrival queue.
        // The 2026-08-11 recording that preceded it was the pre-58eb739f1 fixed
        // point; only the third digit moves, since that fix changes the rate the
        // multiserver correction is charged at, not the regime of this model.
        {"P:P1", 0.0, 0.997132269606374, 0.0, 0.0, 0.0},
        {"P:P2", 0.0, 0.44274892067581, 0.0, 0.0, 0.0},
        {"R:T1", 23.2658534957177, 0.66507056543551, 0.0, 1.07741372970453, 13.3},
        {"T:T2", 8.66752554811615, 0.332061704170864, 0.0, 0.55228201783017,
         13.2824681668346},
        {"T:T3", 1.32824651225174, 0.44274892067581, 0.0, 0.0199999962390167, 66.4},
        {"E:E1", 23.2658534957177, 0.66507056543551, 1.74912668706684, 0.0, 13.3},
        {"E:E2", 8.66752554811615, 0.332061704170864, 0.652553835570891, 0.0,
         13.2824681668346},
        {"E:E3", 1.32824651225174, 0.44274892067581, 0.0199999962390167, 0.0, 66.4},
        {"A:AS1", 23.4, 0.66507056543551, 1.7592870351184,
         1.07741372970453, 13.3},
        {"A:AS2", 8.73491505045899, 0.332061704170864, 0.65762740333679,
         0.55228201783017, 13.2824681668346},
        {"A:AS3", 1.32824651225169, 0.44274892067581, 0.019999996239016,
         0.0199999962390167, 66.4},
    };
    // 1e-3, for the same two reasons as test_ln_multicall.cpp: the reference
    // numbers are the printed, one-decimal-snapped table (8.7 here is one such),
    // and the two codebases stop at slightly different iterates.
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

    // The five-fold call makes E3's throughput five times E2's, because nothing
    // in this model saturates the callee. It used to settle at 0.468 instead,
    // when T3's declared think time of 4 was charged per request and capped it
    // at 25/4.02 = 6.22 completions; that reading is gone from every codebase
    // and no oracle ever supported it (lqsim 66.5, LDES 66.955, lqns 75.6).
    const std::size_t e2 = idx_of(l, "E:E2"), e3 = idx_of(l, "E:E3");
    CHECK(sol.TN[e3] / sol.TN[e2] == doctest::Approx(5.0).epsilon(1e-3));
}

TEST_CASE("lqn_basic: the default encoding matches MATLAB srvn.ph") {
    const lqn::LqnStruct<double> l = build_basic();
    ln::SolverLN<double> s(l, ln::LnOptions());  // the alias, resolving to srvn.ph
    CHECK(s.state_lnmethod() == "srvn.ph");
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    // MATLAB AvgTable with method='srvn.ph', re-recorded 2026-08-14 after
    // 58eb739f1 (the 2026-08-11 recording was its predecessor). The composed
    // entry law puts one class per caller task in a layer, so the multiserver
    // corrections act on a different class mix and every row moves in the third
    // digit against the routing table above -- a different fixed point of the
    // same model, not a looser one.
    struct Row { const char* hn; double qlen, util, respt, residt, tput; };
    const Row rows[] = {
        {"P:P1", 0.0, 0.997106824896032, 0.0, 0.0, 0.0},
        {"P:P2", 0.0, 0.443158596173747, 0.0, 0.0, 0.0},
        {"R:T1", 23.4, 0.664737881458665, 0.0, 1.08476552346506, 13.3},
        {"T:T2", 8.7258653331067, 0.332368943437367, 0.0, 0.556066458158175, 13.3},
        {"T:T3", 1.32947578852125, 0.443158596173747, 0.0, 0.02, 66.5},
        {"E:E1", 23.4, 0.664737881458665, 1.76088083938466, 0.0, 13.3},
        {"E:E2", 8.7258653331067, 0.332368943437367, 0.656338799502716, 0.0, 13.3},
        {"E:E3", 1.32947578852125, 0.443158596173747, 0.02, 0.0, 66.5},
    };
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
