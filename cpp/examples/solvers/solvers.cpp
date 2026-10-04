/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/solvers/`: four solver-comparison walkthroughs.
 *
 * Each one solves ONE model with several engines and prints the theoretical
 * answer beside them, so the value of the example is the commentary, not the
 * table. That commentary is kept section for section; the reference's Greek
 * letters and check glyphs become their ASCII spellings, which is the only
 * change.
 *
 * A PLACE CARRIES AN INERT SERVICE HERE. `Place` in MATLAB and Python needs no
 * service process, but a `qn::NetworkStruct` station with no service marks the
 * class disabled at it; `state_events.h:152` then exempts a Place from that
 * rule, so the declared distribution is never read. It is declared for the same
 * reason `tests/test_spn_rec.cpp` declares one -- to keep the class visiting the
 * place -- and it does not change a firing rate.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"

namespace line {
namespace examples {

namespace {

const char* kGroup = "solvers";

void rule(char c, int n) {
    std::string s(static_cast<std::size_t>(n), c);
    std::printf("%s\n", s.c_str());
}

/** The station index of a named station, which is how the reference selects rows. */
std::size_t station_named(const Sn& sn, const std::string& name) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == name) return i;
    throw InputError("solvers: the model has no station named '" + name + "'");
}

/** The AvgResult of `-s ctmc -a avg`, at the analyzer the CLI's avg arm uses. */
mva::AvgResult<double> ctmc_avg(const Sn& sn, const ctmc::CtmcOptions& opt,
                                std::size_t* states = nullptr,
                                std::size_t* cutoff = nullptr) {
    const ctmc::CtmcAnySolution<double> a = ctmc::solver_ctmc_analyzer_any(sn, opt);
    if (states) *states = a.sol.chain.space.size();
    if (cutoff) {
        std::size_t c = 0;
        for (std::size_t i = 0; i < a.sol.cutoff.size(); ++i)
            c = a.sol.cutoff[i] > c ? a.sol.cutoff[i] : c;
        *cutoff = c;
    }
    return ctmc::solver_ctmc_avg_table(sn, a.sol, opt.method);
}

/**
 * One mode of a Petri-net transition, over `nnodes` nodes.
 *
 * `enabling` is what the mode consumes on firing (the PRE arcs) and `firing`
 * what it produces (the POST arcs); a NEGATIVE firing outcome, which the
 * reference writes at the place a token leaves, is inert here because the
 * consumption is already the enabling arc's. It is transcribed all the same so
 * the two files read alike.
 */
struct SpnMode {
    qn::TransitionParam<double> tp;

    SpnMode(const std::string& name, std::size_t nnodes, const D& firing_proc,
            std::size_t nclasses = 1) {
        tp.nmodes = 1;
        tp.modenames.push_back(name);
        // (node x class): a coloured net moves each class along its own arc.
        tp.enabling.assign(1, Matrix<double>(nnodes, nclasses, 0.0));
        tp.inhibiting.assign(
            1, Matrix<double>(nnodes, nclasses, std::numeric_limits<double>::infinity()));
        tp.firing.assign(1, Matrix<double>(nnodes, nclasses, 0.0));
        tp.nmodeservers.push_back(1.0);
        tp.firingprio.push_back(0.0);
        tp.fireweight.push_back(1.0);
        tp.firingphases.push_back(1);
        tp.timing.push_back(lang::TimingStrategy::TIMED);
        tp.firingproc.push_back(firing_proc);
    }

    SpnMode& enable(std::size_t node, double count, std::size_t cls = 1) {
        tp.enabling[0](node - 1, cls - 1) = count;
        return *this;
    }
    SpnMode& fire(std::size_t node, double count, std::size_t cls = 1) {
        tp.firing[0](node - 1, cls - 1) = count;
        return *this;
    }
};

// ---------------------------------------------------------------------------
// closed_qn_nc_vs_mva
// ---------------------------------------------------------------------------

/**
 * The closed product-form QN both exact solvers are run on.
 *
 * Queue1 mu1 = 1.0, Queue2 mu2 = 0.8, Delay mean 2.0, population N = 3.
 */
Net closed_qn() {
    Net m("closed_qn");
    Queue q1(m, "queue1", SchedStrategy::FCFS);
    Queue q2(m, "queue2", SchedStrategy::FCFS);
    Delay delay(m, "delay");
    ClosedClass jobs(m, "jobs", 3.0, q1);
    m.set_service(q1, jobs, D::exp_mean(1.0));
    q2.set_service(jobs, D::exp_mean(1.25));  // mu = 0.8
    delay.set_service(jobs, D::exp_mean(2.0));
    Routing P;
    cyclic(P, jobs, {q1, q2, delay});
    m.link(P);
    return m;
}

/**
 * Closed Product-Form QN: NC vs MVA Solver Comparison.
 *
 * Demonstrates product-form network analysis using two exact solvers:
 * - NC (Normalizing Constant): Direct exact calculation
 * - MVA (Mean Value Analysis): Iterative convergence to exact solution
 */
void closed_qn_nc_vs_mva() {
    rule('=', 70);
    note("CLOSED PRODUCT-FORM QN: NC vs MVA SOLVER COMPARISON");
    rule('=', 70);

    // Create and solve with NC
    std::printf("\n");
    rule('-', 70);
    note("NC Solver (Normalizing Constant - EXACT)");
    rule('-', 70);
    Net model_nc = closed_qn();
    const nc::NcSolverOptions ncopt;
    const mva::AvgResult<double> result_nc = nc::solver_nc_run_analyzer(model_nc.get_struct(), ncopt);
    std::printf("\nResults:\n");
    print_avg(model_nc.get_struct(), result_nc);

    // Create and solve with MVA
    std::printf("\n");
    rule('-', 70);
    note("MVA Solver (Mean Value Analysis - ITERATIVE)");
    rule('-', 70);
    Net model_mva = closed_qn();
    const mva::MvaOptions mvaopt;
    Matrix<double> init;
    const mva::AvgResult<double> result_mva =
        mva::solver_mva_run_analyzer(model_mva.get_struct(), mvaopt, init);
    std::printf("\nResults:\n");
    print_avg(model_mva.get_struct(), result_mva);

    // Comparison summary
    std::printf("\n");
    rule('=', 70);
    note("VERIFICATION");
    rule('=', 70);

    const Sn& sn = model_nc.get_struct();
    const std::size_t q1 = station_named(sn, "queue1"), q2 = station_named(sn, "queue2");
    const std::size_t dl = station_named(sn, "delay");

    std::printf("\nRESULTS COMPARISON:\n");
    std::printf("  Looking at the tables above:\n");
    std::printf("  Queue1: NC QLen=%.5f, MVA QLen=%.5f\n", result_nc.QN(q1, 0),
                result_mva.QN(q1, 0));
    std::printf("  Queue2: NC QLen=%.5f, MVA QLen=%.5f\n", result_nc.QN(q2, 0),
                result_mva.QN(q2, 0));
    std::printf("  Delay:  NC QLen=%.5f, MVA QLen=%.5f\n", result_nc.QN(dl, 0),
                result_mva.QN(dl, 0));
    std::printf("  Utilizations: %.5f, %.5f, %.5f\n", result_nc.UN(q1, 0), result_nc.UN(q2, 0),
                result_nc.UN(dl, 0));
    std::printf("  Throughputs:  %.5f\n", result_nc.TN(q1, 0));

    std::printf("\nPOPULATION CONSERVATION:\n");
    std::printf("  Total QLen = %.5f + %.5f + %.5f = %.3f\n", result_nc.QN(q1, 0),
                result_nc.QN(q2, 0), result_nc.QN(dl, 0),
                result_nc.QN(q1, 0) + result_nc.QN(q2, 0) + result_nc.QN(dl, 0));

    std::printf("\nFLOW CONSERVATION:\n");
    std::printf("  All stations have throughput = %.5f\n", result_nc.TN(q1, 0));

    std::printf("\nPRODUCT-FORM NETWORK PROPERTIES:\n");
    std::printf("  Jackson's theorem conditions satisfied\n");
    std::printf("  FCFS disciplines at queue stations\n");
    std::printf("  Exponential service times\n");
    std::printf("  Single job class with fixed population (N=3)\n");
    std::printf("  NC and MVA both converge to exact solution\n");

    std::printf("\nSOLVER COMPARISON:\n");
    std::printf("  NC (Normalizing Constant): EXACT solver\n");
    std::printf("    - Direct calculation of steady-state probabilities\n");
    std::printf("    - Computes normalizing constant Z\n");
    std::printf("    - Exact results immediately\n");
    std::printf("  \n");
    std::printf("  MVA (Mean Value Analysis): ITERATIVE solver\n");
    std::printf("    - Iterates using mean value equations\n");
    std::printf("    - Converges to exact solution for product-form\n");
    std::printf("    - Efficient for closed networks\n");

    std::printf("\n");
    rule('=', 70);
}
LINE_EXAMPLE(kGroup, closed_qn_nc_vs_mva);

// ---------------------------------------------------------------------------
// closed_qn_nc_mva_ctmc_spn
// ---------------------------------------------------------------------------

/** The same closed QN as a Stochastic Petri Net (SPN). */
Net closed_qn_spn() {
    Net m("closed_qn_spn");

    // Places represent station buffers
    Place p1(m, "p1");  // Queue1 buffer
    Place p2(m, "p2");  // Queue2 buffer
    Place p3(m, "p3");  // Delay buffer

    // Closed class with 3 jobs starting at p1
    ClosedClass jobs(m, "jobs", 3.0, p1);
    p1.set_service(jobs, Exp(1.0));
    p2.set_service(jobs, Exp(1.0));
    p3.set_service(jobs, Exp(1.0));

    // Transitions represent service completions; the arc vectors are sized
    // against the FULL node count, so the places and the transitions together.
    const std::size_t nnodes = 6;
    // Transition t1: Service at Queue1 (mu1 = 1.0)
    const std::size_t t1 =
        m.add_transition("t1", SpnMode("service1", nnodes, D::exp_mean(1.0))
                                   .enable(p1, 1.0)
                                   .fire(p1, -1.0)
                                   .fire(p2, 1.0)
                                   .tp);
    // Transition t2: Service at Queue2 (mu2 = 0.8)
    const std::size_t t2 =
        m.add_transition("t2", SpnMode("service2", nnodes, D::exp_mean(1.25))
                                   .enable(p2, 1.0)
                                   .fire(p2, -1.0)
                                   .fire(p3, 1.0)
                                   .tp);
    // Transition t3: Delay service (mean = 2.0)
    const std::size_t t3 =
        m.add_transition("t3", SpnMode("service3", nnodes, D::exp_mean(2.0))
                                   .enable(p3, 1.0)
                                   .fire(p3, -1.0)
                                   .fire(p1, 1.0)
                                   .tp);

    // Routing (Places to Transitions to Places)
    Routing P;
    cyclic(P, jobs, {p1, t1, p2, t2, p3, t3});
    m.link(P);

    // Set initial state: all 3 jobs at p1
    p1.set_initial_marking(std::vector<double>{3.0});
    p2.set_initial_marking(std::vector<double>{0.0});
    p3.set_initial_marking(std::vector<double>{0.0});
    return m;
}

void print_results_header(const std::string& solver_name, const std::string& model_type) {
    std::printf("\n");
    rule('-', 70);
    note(model_type.empty() ? solver_name + " Solver" : solver_name + " Solver (" + model_type + ")");
    rule('-', 70);
}

/**
 * Closed Product-Form QN: NC vs MVA vs CTMC, over a traditional QN and an SPN.
 *
 * Demonstrates that the CTMC solver handles Stochastic Petri Net
 * representations of closed product-form networks, matching the exact
 * solutions from NC and MVA.
 */
void closed_qn_nc_mva_ctmc_spn() {
    rule('=', 70);
    note("CLOSED PRODUCT-FORM QN: NC vs MVA vs CTMC SOLVER COMPARISON");
    note("Supporting Both Traditional QN and SPN Representations");
    rule('=', 70);

    Net model_traditional = closed_qn();
    Net model_spn = closed_qn_spn();

    // NC Solver (EXACT) - Traditional QN only
    print_results_header("NC", "Traditional QN - Normalizing Constant (EXACT)");
    const mva::AvgResult<double> result_nc =
        nc::solver_nc_run_analyzer(model_traditional.get_struct(), nc::NcSolverOptions());
    std::printf("\nResults:\n");
    print_avg(model_traditional.get_struct(), result_nc);

    // MVA Solver (EXACT) - Traditional QN only
    print_results_header("MVA", "Traditional QN - Mean Value Analysis (EXACT)");
    Net model_traditional2 = closed_qn();
    Matrix<double> init;
    const mva::AvgResult<double> result_mva =
        mva::solver_mva_run_analyzer(model_traditional2.get_struct(), mva::MvaOptions(), init);
    std::printf("\nResults:\n");
    print_avg(model_traditional2.get_struct(), result_mva);

    // CTMC Solver - Traditional QN
    print_results_header("CTMC", "Traditional QN (State-Space Approximation)");
    Net model_traditional3 = closed_qn();
    std::size_t states_trad = 0;
    const mva::AvgResult<double> result_ctmc_trad =
        ctmc_avg(model_traditional3.get_struct(), ctmc::CtmcOptions(), &states_trad);
    std::printf("\nResults (%zu states):\n", states_trad);
    print_avg(model_traditional3.get_struct(), result_ctmc_trad);

    // CTMC Solver - SPN Model
    print_results_header("CTMC", "SPN Representation (State-Space Approximation)");
    std::size_t states_spn = 0;
    const mva::AvgResult<double> result_ctmc_spn =
        ctmc_avg(model_spn.get_struct(), ctmc::CtmcOptions(), &states_spn);
    std::printf("\nResults (%zu states):\n", states_spn);
    print_avg(model_spn.get_struct(), result_ctmc_spn);

    // Verification and Comparison
    std::printf("\n");
    rule('=', 70);
    note("VERIFICATION AND COMPARISON");
    rule('=', 70);

    // Expected values from NC/MVA (ground truth)
    const char* stations[3] = {"queue1", "queue2", "delay"};
    const char* spn_places[3] = {"p1", "p2", "p3"};
    const char* metrics[4] = {"QLen", "Util", "RespT", "Tput"};
    const double expected[3][4] = {{0.80954, 0.53644, 1.51313, 0.53644},
                                   {1.11759, 0.67055, 2.08844, 0.53644},
                                   {1.07287, 1.07287, 2.0, 0.53644}};

    std::printf("\nRESULTS COMPARISON (tolerance 1%% for exact/CTMC):\n");
    std::printf("\n  Station      Metric    NC         MVA        CTMC(Trad) CTMC(SPN)  Expected\n");
    std::printf("  ");
    rule('-', 88);

    const Sn& sn_trad = model_traditional.get_struct();
    const Sn& sn_spn = model_spn.get_struct();
    bool has_errors = false;

    for (int s = 0; s < 3; ++s) {
        const std::size_t it = station_named(sn_trad, stations[s]);
        // The SPN's places stand for the traditional model's stations, in order.
        const std::size_t ip = station_named(sn_spn, spn_places[s]);
        for (int mi = 0; mi < 4; ++mi) {
            auto pick = [&](const mva::AvgResult<double>& r, std::size_t i) {
                switch (mi) {
                    case 0: return r.QN(i, 0);
                    case 1: return r.UN(i, 0);
                    case 2: return r.RN(i, 0);
                    default: return r.TN(i, 0);
                }
            };
            const double v_nc = pick(result_nc, it), v_mva = pick(result_mva, it);
            const double v_trad = pick(result_ctmc_trad, it), v_spn = pick(result_ctmc_spn, ip);
            const double exp_val = expected[s][mi];
            std::printf("  %-12s %-8s %-10.5g %-10.5g %-10.5g %-10.5g %-10.5g\n", stations[s],
                        metrics[mi], v_nc, v_mva, v_trad, v_spn, exp_val);
            const double pair[2] = {v_trad, v_spn};
            const char* who[2] = {"CTMC_Traditional", "CTMC_SPN"};
            for (int w = 0; w < 2; ++w) {
                if (exp_val == 0.0) continue;
                const double rel_err = std::fabs(pair[w] - exp_val) / std::fabs(exp_val);
                if (rel_err > 0.01) {
                    std::printf("    WARNING: %s %s/%s deviates %.1f%%\n", who[w], stations[s],
                                metrics[mi], rel_err * 100.0);
                    has_errors = true;
                }
            }
        }
    }

    std::printf("\nPOPULATION CONSERVATION:\n");
    {
        const double a = result_nc.QN(station_named(sn_trad, "queue1"), 0);
        const double b = result_nc.QN(station_named(sn_trad, "queue2"), 0);
        const double c = result_nc.QN(station_named(sn_trad, "delay"), 0);
        std::printf("  NC:              %.5f + %.5f + %.5f = %.5f, about 3.0\n", a, b, c,
                    a + b + c);
    }

    std::printf("\nFLOW CONSERVATION (Throughput):\n");
    {
        const char* names[4] = {"NC", "MVA", "CTMC_Traditional", "CTMC_SPN"};
        const double tputs[4] = {result_nc.TN(station_named(sn_trad, "queue1"), 0),
                                 result_mva.TN(station_named(sn_trad, "queue1"), 0),
                                 result_ctmc_trad.TN(station_named(sn_trad, "queue1"), 0),
                                 result_ctmc_spn.TN(station_named(sn_spn, "p1"), 0)};
        double lo = tputs[0], hi = tputs[0];
        for (int i = 0; i < 4; ++i) {
            std::printf("  %-20s: %.5g\n", names[i], tputs[i]);
            lo = tputs[i] < lo ? tputs[i] : lo;
            hi = tputs[i] > hi ? tputs[i] : hi;
        }
        if (lo > 0.0 && (hi - lo) / lo < 0.01) std::printf("  All solvers agree on throughput\n");
        else std::printf("  WARNING: Throughput varies by %.1f%%\n", (hi - lo) / lo * 100.0);
    }

    std::printf("\nSOLVER CHARACTERISTICS:\n");
    std::printf("  NC (Normalizing Constant):\n");
    std::printf("    - Exact solver for product-form networks\n");
    std::printf("    - Direct calculation via normalizing constant Z\n");
    std::printf("    - Traditional QN representation\n\n");
    std::printf("  MVA (Mean Value Analysis):\n");
    std::printf("    - Exact solver for product-form networks\n");
    std::printf("    - Iterative mean value equations\n");
    std::printf("    - Converges to exact solution\n");
    std::printf("    - Traditional QN representation\n\n");
    std::printf("  CTMC (Continuous-Time Markov Chain):\n");
    std::printf("    - State-space solver (approximate via state space truncation/cutoff)\n");
    std::printf("    - Supports both traditional QN and SPN representations\n");
    std::printf("    - Works on arbitrary non-product-form networks\n");
    std::printf("    - For product-form networks with sufficient cutoff, results match exact\n\n");
    std::printf("  SPN (Stochastic Petri Net) support:\n");
    std::printf("    - Places represent station buffers\n");
    std::printf("    - Transitions represent service completions\n");
    std::printf("    - Enabling conditions define job availability\n");
    std::printf("    - Firing outcomes specify job routing\n");
    std::printf("    - CTMC enumerates SPN states and firings\n");

    std::printf("\nSUMMARY:\n");
    if (has_errors)
        std::printf("  Some solvers show deviations >1%% (may indicate state space truncation)\n");
    else
        std::printf("  All solvers agree within 1%% tolerance\n");
    std::printf("  CTMC correctly handles both traditional QN and SPN representations\n");
    std::printf("  SPN-based model produces identical results to traditional QN\n");

    std::printf("\n");
    rule('=', 70);
    std::printf("\n");
}
LINE_EXAMPLE(kGroup, closed_qn_nc_mva_ctmc_spn);

// ---------------------------------------------------------------------------
// ctmc_spn_vs_nc_closed_qn
// ---------------------------------------------------------------------------

/**
 * The same closed product-form network solved twice: as a queueing network by
 * NC, and as a Petri net by CTMC, station by station.
 *
 * This is the two-solver form of `closed_qn_nc_mva_ctmc_spn`: it drops MVA and
 * the traditional-QN CTMC pass and reports the per-station relative deviation
 * of the SPN's CTMC answer from NC's exact one, so the claim being checked is
 * that the modelling STYLE does not move the numbers.
 */
void ctmc_spn_vs_nc_closed_qn() {
    rule('=', 72);
    note("Closed QN: CTMC (SPN) vs NC Solver Comparison");
    rule('=', 72);

    std::printf("\n");
    rule('=', 72);
    note("Approach 1: Traditional Closed QN (Solved with NC Solver)");
    rule('=', 72);

    Net model_qn = closed_qn();
    const Sn& sn_qn = model_qn.get_struct();
    std::printf("\nSolving with NC (Normalizing Constant) solver...\n");
    const mva::AvgResult<double> result_nc = nc::solver_nc_run_analyzer(sn_qn, nc::NcSolverOptions());
    std::printf("NC Results:\n");
    print_avg(sn_qn, result_nc);

    std::printf("\n");
    rule('=', 72);
    note("Approach 2: Same Network as SPN (Solved with CTMC)");
    rule('=', 72);

    Net model_spn = closed_qn_spn();
    const Sn& sn_spn = model_spn.get_struct();
    std::printf("\nSolving with CTMC solver...\n");
    ctmc::CtmcOptions copt;
    copt.method = "exact";
    const mva::AvgResult<double> result_ctmc = ctmc_avg(sn_spn, copt);
    std::printf("CTMC Results:\n");
    print_avg(sn_spn, result_ctmc);

    std::printf("\n");
    rule('=', 72);
    note("COMPARISON: NC Solver (Traditional) vs CTMC (SPN)");
    rule('=', 72);
    std::printf("\nDetailed Metrics Comparison:\n");
    rule('-', 72);

    const char* stations_nc[3] = {"queue1", "queue2", "delay"};
    const char* stations_spn[3] = {"p1", "p2", "p3"};
    for (int i = 0; i < 3; ++i) {
        const std::size_t it = station_named(sn_qn, stations_nc[i]);
        const std::size_t ip = station_named(sn_spn, stations_spn[i]);
        const double nc_qlen = result_nc.QN(it, 0), ctmc_qlen = result_ctmc.QN(ip, 0);
        const double nc_util = result_nc.UN(it, 0), ctmc_util = result_ctmc.UN(ip, 0);
        const double qlen_diff =
            std::fabs(nc_qlen - ctmc_qlen) / std::max(std::fabs(nc_qlen), 1e-6) * 100.0;
        const double util_diff =
            std::fabs(nc_util - ctmc_util) / std::max(std::fabs(nc_util), 1e-6) * 100.0;
        std::printf("\n%s / %s:\n", stations_nc[i], stations_spn[i]);
        std::printf("  %s QLen:  NC=%.6f, CTMC=%.6f (diff=%.3f%%)\n",
                    qlen_diff >= 1.0 ? "X" : "OK", nc_qlen, ctmc_qlen, qlen_diff);
        std::printf("  %s Util:  NC=%.6f, CTMC=%.6f (diff=%.3f%%)\n",
                    util_diff >= 1.0 ? "X" : "OK", nc_util, ctmc_util, util_diff);
    }

    std::printf("\n");
    rule('-', 72);
    note("System-Level Metrics:");
    rule('-', 72);

    double nc_total = 0.0, ctmc_total = 0.0;
    for (std::size_t i = 0; i < sn_qn.nstations; ++i) nc_total += result_nc.QN(i, 0);
    for (std::size_t i = 0; i < sn_spn.nstations; ++i) ctmc_total += result_ctmc.QN(i, 0);
    std::printf("\nTotal QLen (should be 3.0):\n");
    std::printf("  NC:   %.6f\n", nc_total);
    std::printf("  CTMC: %.6f\n", ctmc_total);

    const double nc_tput = result_nc.TN(station_named(sn_qn, "queue1"), 0);
    const double ctmc_tput = result_ctmc.TN(station_named(sn_spn, "p1"), 0);
    const double tput_diff =
        std::fabs(nc_tput - ctmc_tput) / std::max(std::fabs(nc_tput), 1e-6) * 100.0;
    std::printf("\n%s System Throughput:\n", tput_diff >= 1.0 ? "X" : "OK");
    std::printf("  NC:   %.6f\n", nc_tput);
    std::printf("  CTMC: %.6f\n", ctmc_tput);
    std::printf("  Difference: %.3f%%\n", tput_diff);

    std::printf("\n");
    rule('=', 72);
    note("Summary");
    rule('=', 72);
    std::printf("\nOK Both approaches (NC solver on traditional QN, CTMC on SPN)\n");
    std::printf("  give identical results!\n");
    std::printf("\nOK Product-form property verified:\n");
    std::printf("  - Total population conserved (=3 jobs)\n");
    std::printf("  - System throughput matches\n");
    std::printf("  - Individual station metrics match\n");
    std::printf("\nOK CTMC successfully models closed QN as SPN\n");
    std::printf("OK Different modeling approaches yield same answers\n");
    rule('=', 72);
}
LINE_EXAMPLE(kGroup, ctmc_spn_vs_nc_closed_qn);

// ---------------------------------------------------------------------------
// ctmc_spn_mm1
// ---------------------------------------------------------------------------

/**
 * M/M/1 Queue Modeled as Stochastic Petri Net (SPN) with CTMC Solver.
 *
 * - Source: generates arrivals (lambda = 0.8)
 * - Place P: customers waiting in queue
 * - Transition T_serve: service start (mu = 1.0)
 * - Place S: customer in service
 * - Transition T_complete: service completion
 * - Sink: customers depart after service completion
 */
void ctmc_spn_mm1() {
    rule('=', 70);
    note("M/M/1 Queue as Stochastic Petri Net (SPN)");
    rule('=', 70);

    Net m("spn_mm1");

    // SPN Components: a Source for open arrivals, a Sink for departures
    Source source(m, "source");
    Sink sink(m, "sink");

    // SPN places represent queue states
    Place p_queue(m, "queue");      // Customers waiting
    Place p_service(m, "service");  // Customer in service

    // Job class
    OpenClass jobclass(m, "jobs");

    // Define arrival process (lambda = 0.8, mean = 1.25)
    source.set_arrival(jobclass, D::exp_mean(1.0 / 0.8));
    p_queue.set_service(jobclass, Exp(1.0));
    p_service.set_service(jobclass, Exp(1.0));

    // "begin_service": 1 customer in queue, service empty; queue -> service.
    // "complete_service": 1 customer in service; remove it, and route to the sink.
    const std::size_t nnodes = 6;
    const std::size_t t_begin =
        m.add_transition("begin_service", SpnMode("begin", nnodes, D::exp_mean(1.0))
                                              .enable(p_queue, 1.0)
                                              .enable(p_service, 0.0)
                                              .fire(p_queue, -1.0)
                                              .fire(p_service, 1.0)
                                              .tp);
    const std::size_t t_finish =
        m.add_transition("complete_service", SpnMode("finish", nnodes, D::exp_mean(1.0))
                                                 .enable(p_service, 1.0)
                                                 .fire(p_service, -1.0)
                                                 .tp);

    // Routing: source -> queue -> begin_service -> service -> complete -> sink
    Routing P;
    serial(P, jobclass, {source, p_queue, t_begin, p_service, t_finish, sink});
    m.link(P);

    // Solve with CTMC
    std::printf("\nSolving with CTMC solver...\n");
    std::printf("(Note: open network with CTMC uses state space truncation)\n");
    std::size_t states = 0, cutoff = 0;
    const mva::AvgResult<double> avg =
        ctmc_avg(m.get_struct(), ctmc::CtmcOptions(), &states, &cutoff);
    std::printf("\nCTMC Results for SPN M/M/1 (%zu states, cutoff %zu):\n", states, cutoff);
    print_avg(m.get_struct(), avg);

    // Extract key metrics
    const Sn& sn = m.get_struct();
    const std::size_t iq = station_named(sn, "queue"), is = station_named(sn, "service");
    const double qlen_queue = avg.QN(iq, 0);
    const double qlen_service = avg.QN(is, 0);
    const double util_service = avg.UN(is, 0);
    const double tput = avg.TN(iq, 0);

    std::printf("\n");
    rule('=', 70);
    note("Key Metrics:");
    rule('=', 70);
    std::printf("  Queue (waiting):       QLen=%.6f\n", qlen_queue);
    std::printf("  Service (in progress): QLen=%.6f, Util=%.6f\n", qlen_service, util_service);
    std::printf("  System Throughput:     %.6f\n", tput);

    // Theoretical M/M/1 for comparison
    const double lambda_rate = 0.8, mu = 1.0;
    const double rho = lambda_rate / mu;
    const double L = rho / (1 - rho);       // Average number in system
    const double Lq = L - rho;              // Average number waiting
    const double W = 1 / (mu * (1 - rho));  // Average time in system
    const double Wq = W - 1 / mu;           // Average waiting time

    std::printf("\n");
    rule('=', 70);
    note("Theoretical M/M/1 (lambda=0.8, mu=1.0, rho=0.8):");
    rule('=', 70);
    std::printf("  Utilization (rho):      %.6f\n", rho);
    std::printf("  Queue Length (Lq):      %.6f\n", Lq);
    std::printf("  System Length (L):      %.6f\n", L);
    std::printf("  Response Time (W):      %.6f\n", W);
    std::printf("  Waiting Time (Wq):      %.6f\n", Wq);

    std::printf("\n");
    rule('=', 70);
    note("Summary");
    rule('=', 70);
    std::printf("CTMC analyzes M/M/1 modeled as SPN\n");
    std::printf("SPN feature support in CTMC is working\n");
    std::printf("  Differences from theory are due to CTMC state space truncation\n");
    std::printf("  (open networks with infinite arrivals truncated at cutoff=%zu)\n", cutoff);
    rule('=', 70);
}
LINE_EXAMPLE(kGroup, ctmc_spn_mm1);

// ---------------------------------------------------------------------------
// ctmc_tandem_mm1
// ---------------------------------------------------------------------------

/** Compare a computed value with its theoretical counterpart. */
double percent_diff(double computed, double theoretical) {
    if (theoretical == 0.0) return computed == 0.0 ? 0.0 : 100.0;
    return std::fabs(computed - theoretical) / theoretical * 100.0;
}

/**
 * Tandem M/M/1 Queues (Series of Two M/M/1 Stations).
 *
 * - Station 1: lambda1 = 0.6, mu1 = 1.0, rho1 = 0.6
 * - Station 2: lambda2 = 0.6, mu2 = 1.2, rho2 = 0.5
 * The system is product-form: pi(n1,n2) = pi1(n1) pi2(n2).
 */
void ctmc_tandem_mm1() {
    rule('=', 70);
    note("Tandem M/M/1 Queues (Series of Two Stations)");
    rule('=', 70);

    Net m("tandem_mm1");
    Source source(m, "source");
    Queue queue1(m, "queue1", SchedStrategy::FCFS);
    Queue queue2(m, "queue2", SchedStrategy::FCFS);
    Sink sink(m, "sink");
    OpenClass jobclass(m, "jobs");

    // External arrivals lambda = 0.6; mu1 = 1.0 and mu2 = 1.2 downstream.
    source.set_arrival(jobclass, D::exp_mean(1 / 0.6));
    queue1.set_service(jobclass, D::exp_mean(1.0));
    queue2.set_service(jobclass, D::exp_mean(1 / 1.2));

    Routing P;
    serial(P, jobclass, {source, queue1, queue2, sink});
    m.link(P);

    std::printf("\nSolving tandem M/M/1 system with CTMC...\n");
    std::printf("Parameters:\n");
    std::printf("  Station 1: lambda=0.6, mu=1.0, rho=0.6\n");
    std::printf("  Station 2: lambda=0.6, mu=1.2, rho=0.5\n\n");

    std::size_t states = 0, cutoff = 0;
    const mva::AvgResult<double> avg =
        ctmc_avg(m.get_struct(), ctmc::CtmcOptions(), &states, &cutoff);
    std::printf("CTMC Results for Tandem M/M/1 (%zu states, cutoff %zu):\n", states, cutoff);
    print_avg(m.get_struct(), avg);

    const Sn& sn = m.get_struct();
    const std::size_t i1 = station_named(sn, "queue1"), i2 = station_named(sn, "queue2");
    const double q1_qlen = avg.QN(i1, 0), q1_util = avg.UN(i1, 0);
    const double q1_respt = avg.RN(i1, 0), q1_tput = avg.TN(i1, 0);
    const double q2_qlen = avg.QN(i2, 0), q2_util = avg.UN(i2, 0);
    const double q2_respt = avg.RN(i2, 0), q2_tput = avg.TN(i2, 0);

    std::printf("\n");
    rule('=', 70);
    note("CTMC Results Summary:");
    rule('=', 70);
    std::printf("\nQueue 1:\n");
    std::printf("  Queue Length:   %.6f\n", q1_qlen);
    std::printf("  Utilization:    %.6f\n", q1_util);
    std::printf("  Response Time:  %.6f\n", q1_respt);
    std::printf("  Throughput:     %.6f\n", q1_tput);
    std::printf("\nQueue 2:\n");
    std::printf("  Queue Length:   %.6f\n", q2_qlen);
    std::printf("  Utilization:    %.6f\n", q2_util);
    std::printf("  Response Time:  %.6f\n", q2_respt);
    std::printf("  Throughput:     %.6f\n", q2_tput);

    // Theoretical M/M/1 results (product-form)
    const double lambda_rate = 0.6, mu1 = 1.0, mu2 = 1.2;
    const double rho1 = lambda_rate / mu1, rho2 = lambda_rate / mu2;
    const double L1 = rho1 / (1 - rho1), L2 = rho2 / (1 - rho2);
    const double W1 = 1 / (mu1 * (1 - rho1)), W2 = 1 / (mu2 * (1 - rho2));
    const double L_system = L1 + L2, W_system = W1 + W2;

    std::printf("\n");
    rule('=', 70);
    note("Theoretical Results (Product-Form M/M/1-M/M/1):");
    rule('=', 70);
    std::printf("\nQueue 1 (lambda=0.6, mu=1.0, rho=0.6):\n");
    std::printf("  Queue Length (L1):      %.6f\n", L1);
    std::printf("  Response Time (W1):     %.6f\n", W1);
    std::printf("  Utilization (rho1):     %.6f\n", rho1);
    std::printf("\nQueue 2 (lambda=0.6, mu=1.2, rho=0.5):\n");
    std::printf("  Queue Length (L2):      %.6f\n", L2);
    std::printf("  Response Time (W2):     %.6f\n", W2);
    std::printf("  Utilization (rho2):     %.6f\n", rho2);
    std::printf("\nSystem Totals:\n");
    std::printf("  Total Queue Length:     %.6f\n", L_system);
    std::printf("  Total Response Time:    %.6f\n", W_system);

    // Comparison and Validation
    std::printf("\n");
    rule('=', 70);
    note("COMPARISON: CTMC vs Theory");
    rule('=', 70);
    auto row = [](const char* label, double computed, double theoretical) {
        const double d = percent_diff(computed, theoretical);
        std::printf("%s %s: CTMC=%.6f, Theory=%.6f (diff=%.2f%%)\n", d <= 5.0 ? "OK  " : "FAIL",
                    label, computed, theoretical, d);
    };
    std::printf("\nQueue 1:\n");
    row("Queue Length", q1_qlen, L1);
    row("Response Time", q1_respt, W1);
    row("Utilization", q1_util, rho1);
    std::printf("\nQueue 2:\n");
    row("Queue Length", q2_qlen, L2);
    row("Response Time", q2_respt, W2);
    row("Utilization", q2_util, rho2);
    std::printf("\nSystem:\n");
    row("Total QLen", q1_qlen + q2_qlen, L_system);

    std::printf("\n");
    rule('=', 70);
    note("Summary");
    rule('=', 70);
    std::printf("\nCTMC analyzes the tandem M/M/1-M/M/1 system\n");
    std::printf("Product-form property verified:\n");
    std::printf("  - Each station behaves independently as M/M/1\n");
    std::printf("  - System decomposition holds\n");
    std::printf("Results match theoretical predictions up to the state-space cutoff\n");
    rule('=', 70);
}
LINE_EXAMPLE(kGroup, ctmc_tandem_mm1);

}  // namespace

}  // namespace examples
}  // namespace line
