/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Unit tests of line/solvers/ln/lqn_analyzers.h: the reduced overtaking CTMC,
 * the input mapping onto the LQNS overtaking chain, and the Robbins-Monro /
 * Polyak-Ruppert controller for noisy layer solvers.
 *
 * WHY THE ORACLES ARE WHAT THEY ARE. None of the three routines is a solver, so
 * none of them has an AvgTable to check; every expected value below is an
 * identity that holds whatever the implementation, or a formula quoted from the
 * MATLAB source with the file and lines named. Nothing is a number read back out
 * of this implementation.
 *
 * The three-state chain is pinned by its own global balance equations, which is
 * the definition of a stationary law and says nothing about how it is computed.
 * The mapping layer is pinned by the column contract its reference documents
 * (nSlices, service, y_ij, y_ik, t_k, per client phase): the test assembles that
 * matrix by hand from the model's declared parameters and requires the mapping
 * to feed the chain the same one. The controller is pinned by the two exact
 * properties it exists to have, that its running average IS the arithmetic mean
 * and that it refuses to stop until the drift has been under tolerance for the
 * configured number of consecutive iterations.
 *
 * NOT tested here: that any of this is wired into SolverLN, because it is not.
 * solver_ln.h refuses phase-2 activities by name, so no solved layered model in
 * this suite can reach the overtaking path at all.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/lqn_analyzers.h"

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

/** Index of the call from `src` to `dst`, whatever its type. */
std::size_t call_of(const lqn::LqnStruct<double>& l, std::size_t src, std::size_t dst) {
    for (std::size_t c = 1; c <= l.ncalls; ++c)
        if (l.callpair_src[c] == src && l.callpair_dst[c] == dst) return c;
    return 0;
}

/**
 * T1 (reference, think time Z) calls E2 on T2 twice and E3 on T3 once. T2 is the
 * server whose entry is tested, T3 is the "other task" whose delay lands in t_k.
 */
lqn::LqnStruct<double> three_task_client(double Z, double y_to_e2, double y_to_e3,
                                         bool async_to_e2 = false) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 4, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(Z));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.1), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    if (async_to_e2)
        b.async_call("A1", "E2", y_to_e2);
    else
        b.sync_call("A1", "E2", y_to_e2);
    b.sync_call("A1", "E3", y_to_e3);
    b.activity("A2", D::exp_mean(0.2), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.3), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** One layer's results, all metrics 1x1 so the running average is inspectable. */
ln::LayerResult<double> flat_result(double q) {
    ln::LayerResult<double> r;
    r.QN = Matrix<double>(1, 1, q);
    r.UN = Matrix<double>(1, 1, 0.0);
    r.RN = Matrix<double>(1, 1, 0.0);
    r.TN = Matrix<double>(1, 1, 0.0);
    r.WN = Matrix<double>(1, 1, 0.0);
    return r;
}

}  // namespace

// ---------------------------------------------------------------------------
// overtake_prob
// ---------------------------------------------------------------------------

TEST_CASE("lqnanl: the three-state overtaking chain satisfies its balance equations") {
    const double S1 = 0.4, S2 = 0.6, lambda = 0.8;
    const ln::OvertakeCtmcState<double> pi = ln::lqn_overtake_ctmc(S1, S2, lambda);

    // The defining equations of the stationary law of the cycle
    // idle -> phase 1 -> phase 2 -> idle: every cut carries the same flow, and
    // the three probabilities sum to one. Both are independent of how the law
    // was obtained.
    CHECK(pi.idle + pi.phase1 + pi.phase2 == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(lambda * pi.idle == doctest::Approx(pi.phase1 / S1).epsilon(1e-14));
    CHECK(pi.phase1 / S1 == doctest::Approx(pi.phase2 / S2).epsilon(1e-14));
    CHECK(ln::lqn_overtake_prob(S1, S2, lambda, 1.0) ==
          doctest::Approx(pi.phase2).epsilon(1e-14));
}

TEST_CASE("lqnanl: a phase that cannot be observed carries no overtaking") {
    // Each guard removes one of the three intervals the chain is built from, so
    // the answer is zero by absence of the event, not by a small number.
    CHECK(ln::lqn_overtake_prob(0.4, 0.0, 0.8, 1.0) == 0.0);
    CHECK(ln::lqn_overtake_prob(0.0, 0.6, 0.8, 1.0) == 0.0);
    CHECK(ln::lqn_overtake_prob(0.4, 0.6, 0.0, 1.0) == 0.0);
}

TEST_CASE("lqnanl: a saturated single server splits its busy time in the phase ratio") {
    const double S1 = 0.4, S2 = 0.6;

    // As the arrival rate grows the idle state vanishes and the server alternates
    // phase 1 and phase 2 forever, so by renewal reward the time fraction in
    // phase 2 tends to S2/(S1+S2). The approach is monotone in lambda.
    const double lo = ln::lqn_overtake_prob(S1, S2, 0.1, 1.0);
    const double mid = ln::lqn_overtake_prob(S1, S2, 1.0, 1.0);
    const double hi = ln::lqn_overtake_prob(S1, S2, 10.0, 1.0);
    CHECK(lo < mid);
    CHECK(mid < hi);
    CHECK(ln::lqn_overtake_prob(S1, S2, 1e9, 1.0) ==
          doctest::Approx(S2 / (S1 + S2)).epsilon(1e-6));
}

TEST_CASE("lqnanl: the multiserver branch is the phase fraction scaled by the load") {
    // overtake_prob.m lines 77 to 91: rho = lambda*(S1+S2)/c, and the answer is
    // (S2/(S1+S2))*rho below saturation and S2/(S1+S2) at or above it. This
    // pins the branch selection; the approximation itself is the reference's.
    const double S1 = 0.4, S2 = 0.6, frac = S2 / (S1 + S2);

    // rho = 1*1.0/2 = 0.5
    CHECK(ln::lqn_overtake_prob(S1, S2, 1.0, 2.0) == doctest::Approx(frac * 0.5).epsilon(1e-12));
    // rho = 4*1.0/2 = 2, saturated
    CHECK(ln::lqn_overtake_prob(S1, S2, 4.0, 2.0) == doctest::Approx(frac).epsilon(1e-12));
    // the two branches must agree where they meet, at rho = 2*1.0/2 = 1
    CHECK(ln::lqn_overtake_prob(S1, S2, 2.0, 2.0) == doctest::Approx(frac).epsilon(1e-12));
}

TEST_CASE("lqnanl: an infinite-multiplicity server is never overtaken") {
    const double inf = std::numeric_limits<double>::infinity();
    // Every arrival gets an idle server, so the load per server is zero and the
    // reference's own expression evaluates to zero rather than to a limit.
    CHECK(ln::lqn_overtake_prob(0.4, 0.6, 5.0, inf) == 0.0);
}

// ---------------------------------------------------------------------------
// overtake_prob_markov
// ---------------------------------------------------------------------------

TEST_CASE("lqnanl: the overtaking mapping needs a phase-2 residence to work with") {
    lqn::LqnStruct<double> l = three_task_client(3.0, 2.0, 1.0);
    std::vector<double> servt(l.nidx + 1, 0.5), tput(l.nidx + 1, 1.0),
        callresidt(l.ncalls + 1, 0.7);

    CHECK(ln::lqn_overtake_prob_markov(l, servt, callresidt, tput, idx_of(l, "E:E2"), 0.0) ==
          0.0);
}

TEST_CASE("lqnanl: an entry nobody calls cannot be overtaken") {
    lqn::LqnStruct<double> l = three_task_client(3.0, 2.0, 1.0);
    std::vector<double> servt(l.nidx + 1, 0.5), tput(l.nidx + 1, 1.0),
        callresidt(l.ncalls + 1, 0.7);

    // E1 is the reference task's own entry, so no call has it as a destination.
    CHECK(ln::lqn_overtake_prob_markov(l, servt, callresidt, tput, idx_of(l, "E:E1"), 0.6) ==
          0.0);
}

TEST_CASE("lqnanl: an asynchronous caller contributes no overtaking") {
    lqn::LqnStruct<double> l = three_task_client(3.0, 2.0, 1.0, /*async_to_e2=*/true);
    std::vector<double> servt(l.nidx + 1, 0.5), tput(l.nidx + 1, 1.0),
        callresidt(l.ncalls + 1, 0.7);

    // Nobody blocks on a send-no-reply, so the caller has no second outstanding
    // request that could pass a first one. The reference reaches the same answer
    // by a different route: its caller selection is unfiltered, so E1 IS taken as
    // a client entry, and the per-phase scan then finds no synchronous calls to
    // the server task and drops the contribution.
    CHECK(ln::lqn_overtake_prob_markov(l, servt, callresidt, tput, idx_of(l, "E:E2"), 0.6) ==
          0.0);
}

TEST_CASE("lqnanl: the mapping builds the client phases its reference documents") {
    const double Z = 3.0, y2 = 2.0, y3 = 1.0, xj = 0.6;
    lqn::LqnStruct<double> l = three_task_client(Z, y2, y3);
    const std::size_t a1 = idx_of(l, "A:A1"), e1 = idx_of(l, "E:E1"), e2 = idx_of(l, "E:E2");
    // "R:T1": three_task_client declares T1 SchedStrategy::REF, and a reference
    // task hashes with the R prefix (getStruct.m:129). With "T:T1" idx_of misses
    // and returns 0, so `tput[t1] = 5.0` below wrote element 0 and the oracle's
    // prVisit was 3/5 while the function, reading the real task index, saw a
    // throughput nobody had set.
    const std::size_t e3 = idx_of(l, "E:E3"), t1 = idx_of(l, "R:T1");

    std::vector<double> servt(l.nidx + 1, 0.0), tput(l.nidx + 1, 0.0),
        callresidt(l.ncalls + 1, 0.0);
    servt[a1] = 1.2;
    tput[t1] = 5.0;
    tput[e1] = 3.0;
    callresidt[call_of(l, a1, e3)] = 0.9;

    // Columns are [nSlices service y_ij y_ik t_k] per client phase, row 0 the
    // think slice, as overtake_prob_markov.m documents at lines 61 to 96. A1 is
    // the client's only activity and is phase 1: it is cut into one slice per
    // rendezvous plus one, two of its calls reach the server task T2 and one
    // reaches T3, whose mean delay is that call's residence time.
    Matrix<double> expect({{1.0, Z, 0.0, 0.0, 0.0},
                           {1.0 + y2 + y3, servt[a1], y2, y3, callresidt[call_of(l, a1, e3)]}});
    const std::vector<double> y_aj{y2, y2};
    const double prVisit = tput[e1] / tput[t1];
    const double oracle = ln::lqn_overtake_markov(expect, prVisit, xj, y_aj);

    CHECK(ln::lqn_overtake_prob_markov(l, servt, callresidt, tput, e2, xj) ==
          doctest::Approx(oracle).epsilon(1e-12));
    CHECK(oracle > 0.0);
}

TEST_CASE("lqnanl: the mapping separates a two-phase client into its phases") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 4, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(2.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", D::exp_mean(0.1), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 2.0);
    b.activity("A1b", D::exp_mean(0.1), "T1");
    b.sync_call("A1b", "E2", 1.0);
    b.serial("A1", "A1b");
    b.replies_to("A1b", "E1");
    b.activity("A2", D::exp_mean(0.2), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    lqn::LqnStruct<double> l = b.build();

    const std::size_t a1 = idx_of(l, "A:A1"), a1b = idx_of(l, "A:A1b");
    const std::size_t e2 = idx_of(l, "E:E2");
    // The builder has no phase setter because solver_ln.h refuses phase-2
    // activities outright; the whole point of the overtaking model is the case
    // it refuses, so the phase is set on the struct here.
    l.actphase[a1b - l.ashift] = 2;

    std::vector<double> servt(l.nidx + 1, 0.0), tput(l.nidx + 1, 0.0),
        callresidt(l.ncalls + 1, 0.0);
    servt[a1] = 1.2;
    servt[a1b] = 0.4;
    const double xj = 0.5;

    // Each phase gets its own row: its own host residence, its own slice count
    // and its own share of the calls, and y_aj records the split that the chain
    // conditions on.
    Matrix<double> expect({{1.0, 2.0, 0.0, 0.0, 0.0},
                           {1.0 + 2.0, servt[a1], 2.0, 0.0, 0.0},
                           {1.0 + 1.0, servt[a1b], 1.0, 0.0, 0.0}});
    const std::vector<double> y_aj{3.0, 2.0, 1.0};
    // both throughputs left at zero, so the reference falls back to prVisit = 1
    const double raw = ln::lqn_overtake_markov(expect, 1.0, xj, y_aj);
    // lqn_overtake_markov is the KERNEL and lqn_overtake_prob_markov is the
    // PROBABILITY, which is the kernel clamped to [0,1]: overtake_prob_markov.m
    // ends with max(0, min(1, prOt)) and the port reproduces it at
    // lqn_analyzers.h:272-273. Native Python clamps the same quantity at both
    // return sites of _overtake_prob, from a different algorithm. By PASTA this
    // IS a probability (the steady-state chance of being in phase 2), so the
    // raw 1.126005 this model produces has no admissible reading; comparing the
    // two functions unclamped asserts an equality that holds only where the
    // kernel happens to stay below 1.
    CHECK(raw > 1.0);  // this model saturates, which is what makes the clamp load bearing
    const double oracle = std::min(1.0, raw);

    CHECK(ln::lqn_overtake_prob_markov(l, servt, callresidt, tput, e2, xj) ==
          doctest::Approx(oracle).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// convergedStoch
// ---------------------------------------------------------------------------

TEST_CASE("lqnanl: the Polyak running average is the arithmetic mean of the samples") {
    ln::LnStochConfig cfg;
    cfg.burnin = 0;
    cfg.conseq = 1000;  // never stop, so every sample is folded in
    ln::LnStochController<double> ctl(cfg);
    const std::vector<double> jobs{10.0};

    double sum = 0.0;
    for (int it = 1; it <= 4; ++it) {
        const double q = double(it);
        sum += q;
        std::vector<ln::LayerResult<double>> latest{flat_result(q)};
        const std::vector<double> servt{q, 2.0 * q}, residt{3.0 * q, 4.0 * q};
        ctl.update(it, latest, jobs, servt, residt);

        // m_k = m_{k-1} + (x_k - m_{k-1})/k is the mean of x_1..x_k exactly, and
        // that identity is the only reason the running form is allowed at all.
        const double mean = sum / double(it);
        CHECK(ctl.averaging_count() == it);
        CHECK(ctl.averaged_results()[0].QN(0, 0) == doctest::Approx(mean).epsilon(1e-14));
        CHECK(ctl.averaged_servt()[1] == doctest::Approx(2.0 * mean).epsilon(1e-14));
        CHECK(ctl.averaged_residt()[0] == doctest::Approx(3.0 * mean).epsilon(1e-14));
    }
}

TEST_CASE("lqnanl: the Robbins-Monro step decays on the documented schedule") {
    ln::LnStochConfig cfg;
    cfg.burnin = 5;
    cfg.a0 = 1.0;
    cfg.alpha = 0.6;
    cfg.conseq = 1000;
    cfg.relax_burnin = 0.5;
    ln::LnStochController<double> ctl(cfg);
    const std::vector<double> jobs{10.0};
    const std::vector<double> none;

    double prev = 1.0;
    for (int it = 1; it <= 10; ++it) {
        std::vector<ln::LayerResult<double>> latest{flat_result(1.0)};
        ctl.update(it, latest, jobs, none, none);
        const double w = ctl.relax_omega();
        if (it < cfg.burnin) {
            // untouched during burn-in: whatever init chose stays in force
            CHECK(w == doctest::Approx(0.5).epsilon(1e-14));
        } else {
            // convergedStoch.m line 49: omega = min(1, a0/(it-burnin+1)^alpha)
            const double k = double(it - cfg.burnin + 1);
            CHECK(w == doctest::Approx(std::min(1.0, cfg.a0 / std::pow(k, cfg.alpha)))
                           .epsilon(1e-14));
            CHECK(w <= 1.0);
            CHECK(w <= prev);
            prev = w;
        }
    }
}

TEST_CASE("lqnanl: burn-in never reports convergence") {
    ln::LnStochConfig cfg;
    cfg.burnin = 5;
    cfg.conseq = 1;
    cfg.iter_tol = 1e9;  // any drift would pass the test, if one were computed
    ln::LnStochController<double> ctl(cfg);
    const std::vector<double> jobs{10.0};
    const std::vector<double> none;

    for (int it = 1; it <= 5; ++it) {
        std::vector<ln::LayerResult<double>> latest{flat_result(1.0)};
        CHECK_FALSE(ctl.update(it, latest, jobs, none, none));
        // no average exists yet, so the iteration cannot be said to have drifted
        // by a small amount; the reference records an infinite error instead
        CHECK(std::isinf(ctl.iteration_error()[std::size_t(it)]));
        CHECK(ctl.averaging_count() == 0);
        CHECK(ctl.averaging_start() == -1);
    }
}

TEST_CASE("lqnanl: stopping needs the drift under tolerance consecutively") {
    ln::LnStochConfig cfg;
    cfg.burnin = 2;
    cfg.conseq = 3;
    cfg.iter_tol = 1e-3;
    ln::LnStochController<double> ctl(cfg);
    const std::vector<double> jobs{10.0};
    const std::vector<double> none;

    // Noise-free layers: the average stops moving immediately, so the drift is
    // zero from the first averaged iteration onwards and the only thing left
    // holding the test back is the consecutive-iteration requirement.
    for (int it = 1; it <= cfg.burnin + cfg.conseq; ++it) {
        std::vector<ln::LayerResult<double>> latest{flat_result(1.0)};
        CHECK_FALSE(ctl.update(int(it), latest, jobs, none, none));
    }
    std::vector<ln::LayerResult<double>> latest{flat_result(1.0)};
    CHECK(ctl.update(int(cfg.burnin + cfg.conseq + 1), latest, jobs, none, none));
    CHECK(ctl.averaging_start() == cfg.burnin + 1);
}

TEST_CASE("lqnanl: a drift at the last iteration blocks the stop") {
    ln::LnStochConfig cfg;
    cfg.burnin = 2;
    cfg.conseq = 3;
    cfg.iter_tol = 1e-3;
    ln::LnStochController<double> ctl(cfg);
    const std::vector<double> jobs{10.0};
    const std::vector<double> none;

    for (int it = 1; it <= cfg.burnin + cfg.conseq; ++it) {
        std::vector<ln::LayerResult<double>> latest{flat_result(1.0)};
        ctl.update(int(it), latest, jobs, none, none);
    }
    // Same iteration that stopped the run above, but this layer moved: the
    // window must include the current iteration, or a run stops on the strength
    // of iterations that predate the disturbance.
    std::vector<ln::LayerResult<double>> latest{flat_result(100.0)};
    CHECK_FALSE(ctl.update(int(cfg.burnin + cfg.conseq + 1), latest, jobs, none, none));
}
