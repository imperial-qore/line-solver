/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/workflowModels/`: seven activity workflows reduced to a
 * phase-type distribution.
 *
 * Each builds the reference's precedence graph on `line/lang/workflow/workflow.h`
 * -- `add_activity`, `add_precedence` over
 * `Serial / AndFork / AndJoin / OrFork / OrJoin / Loop`, then `to_ph()` -- and
 * prints the same mean, SCV and phase count the reference `__main__` prints.
 *
 * TWO THINGS DO NOT CROSS, and both are named at their site rather than
 * silently substituted:
 *
 *   - `wf_complex` ends by drawing 10000 samples from the composed law and
 *     printing their moments. Those are a property of numpy's RNG STREAM at
 *     `np.random.seed(1241)`, not of the model, so no other implementation can
 *     reproduce them digit for digit; the analytic moments of the same law are
 *     printed instead and the sampling call is refused by name.
 *
 *   - An LQN fork-join is NOT a substitute for `AndFork`/`AndJoin` here: its
 *     branches contend for a host, so it yields a response time under
 *     contention rather than the order statistic of independent activity times
 *     that `to_ph` returns.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/workflow/workflow.h"

namespace line {
namespace examples {

namespace {

namespace wfl = line::workflow;
using Dist = line::lang::Distrib<double>;
using Wf = wfl::Workflow<double>;

/** E[max(X, Y)] for independent exponentials, the identity the scripts print. */
double exp_max_mean(double mean_x, double mean_y) {
    return mean_x + mean_y - 1.0 / (1.0 / mean_x + 1.0 / mean_y);
}

/**
 * Mean and SCV of a phase-type law, by -alpha S^-1 e and 2 alpha S^-2 e.
 *
 * Gauss with partial pivoting on S rather than an explicit inverse: the
 * composed generators here are block triangular and badly scaled once an
 * immediate phase of rate 1e8 is in them.
 */
void ph_moments(const Dist& d, double& mean, double& scv, std::size_t& phases) {
    const std::size_t n = d.D0.rows();
    phases = n;

    std::vector<std::vector<double>> M(n, std::vector<double>(n + 1, 0.0));
    std::vector<double> x(n, 0.0), y(n, 0.0);

    auto solve = [&](const std::vector<double>& rhs, std::vector<double>& out) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) M[i][j] = d.D0(i, j);
            M[i][n] = rhs[i];
        }
        for (std::size_t c = 0; c < n; ++c) {
            std::size_t piv = c;
            for (std::size_t r = c + 1; r < n; ++r)
                if (std::abs(M[r][c]) > std::abs(M[piv][c])) piv = r;
            std::swap(M[c], M[piv]);
            for (std::size_t r = 0; r < n; ++r) {
                if (r == c) continue;
                const double f = M[r][c] / M[c][c];
                for (std::size_t j = c; j <= n; ++j) M[r][j] -= f * M[c][j];
            }
        }
        out.assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) out[i] = M[i][n] / M[i][i];
    };

    solve(std::vector<double>(n, -1.0), x);  // x = -S^-1 e
    std::vector<double> minus_x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) minus_x[i] = -x[i];
    solve(minus_x, y);  // y = S^-2 e

    double m1 = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        m1 += d.params[i] * x[i];
        m2 += 2.0 * d.params[i] * y[i];
    }
    mean = m1;
    scv = (m2 - m1 * m1) / (m1 * m1);
}

/**
 * The phase-type moments an example prints, printed AND recorded.
 *
 * The seven `wf_*` goldens hold nothing but these three numbers, keyed
 * ('PH','mean'), ('PH','SCV') and ('PH','phases') under the shape key `WF` --
 * a workflow is not a solver, so no solver name owns them. The Python twin
 * recovers the same rows by hooking `Workflow.getMean` / `getSCV` / `toPH`; a
 * C++ example computes them itself, so it declares them itself. `scv` is
 * recorded only where the reference prints it: `wf_branch` and its siblings
 * print the mean and the order alone, and their goldens hold exactly that.
 */
void print_ph_moments(double mean, double scv, std::size_t phases, bool with_scv) {
    std::printf("Computed PH mean: %.4f\n", mean);
    derived("WF", "PH", "mean", mean);
    if (with_scv) {
        std::printf("Computed PH SCV: %.4f\n", scv);
        derived("WF", "PH", "SCV", scv);
    }
    std::printf("Number of phases: %zu\n", phases);
    derived("WF", "PH", "phases", static_cast<double>(phases));
}

}  // namespace

/** A -> B -> C, three exponential activities in series. */
void wf_serial() {
    Wf wf("SerialWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(1.5));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.add_precedence(Wf::Serial("B", "C"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    std::printf("Serial Workflow: A -> B -> C\n");
    std::printf("Activity means: A=1.0, B=2.0, C=1.5\n");
    std::printf("Expected total mean: %.2f\n", 1.0 + 2.0 + 1.5);
    print_ph_moments(mean, scv, phases, true);
}

/** A -> [B || C] -> D, an AND fork/join over two exponential branches. */
void wf_parallel() {
    Wf wf("ParallelWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(3.0));
    wf.add_activity("D", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::AndFork("A", {"B", "C"}));
    wf.add_precedence(Wf::AndJoin({"B", "C"}, "D"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    const double expected_max = exp_max_mean(2.0, 3.0);
    std::printf("Parallel Workflow: A -> [B || C] -> D\n");
    std::printf("Activity means: A=1.0, B=2.0, C=3.0, D=0.5\n");
    std::printf("Expected max(B,C) mean: %.4f\n", expected_max);
    std::printf("Expected total mean: %.4f\n", 1.0 + expected_max + 0.5);
    print_ph_moments(mean, scv, phases, false);
}

/** A -> [B at 60% | C at 40%] -> D, an OR fork/join. */
void wf_branch() {
    Wf wf("BranchingWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(5.0));
    wf.add_activity("D", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::OrFork("A", {"B", "C"}, {0.6, 0.4}));
    wf.add_precedence(Wf::OrJoin({"B", "C"}, "D"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    const double expected_branch = 0.6 * 2.0 + 0.4 * 5.0;
    std::printf("Branching Workflow: A -> [B(60%%) | C(40%%)] -> D\n");
    std::printf("Activity means: A=1.0, B=2.0, C=5.0, D=0.5\n");
    std::printf("Expected branch mean: 0.6*2.0 + 0.4*5.0 = %.2f\n", expected_branch);
    std::printf("Expected total mean: %.2f\n", 1.0 + expected_branch + 0.5);
    print_ph_moments(mean, scv, phases, false);
}

/**
 * A -> [B x 3] -> C, a loop body repeated a GEOMETRIC number of times of mean 3.
 *
 * Not a 3-fold convolution: POST_LOOP takes the back edge with probability
 * 1-1/3, so the composed law keeps the ORDER of the body and its mean is still
 * 3 times the body mean.
 */
void wf_loop() {
    Wf wf("LoopWorkflow");
    wf.add_activity("A", Dist::exp_mean(1.0));
    wf.add_activity("B", Dist::exp_mean(2.0));
    wf.add_activity("C", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::Loop("A", {"B", "C"}, 3.0));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    std::printf("Loop Workflow: A -> [B x 3] -> C\n");
    std::printf("Activity means: A=1.0, B=2.0, C=0.5\n");
    std::printf("Expected total mean: 1.0 + 3*2.0 + 0.5 = %.2f\n", 1.0 + 3.0 * 2.0 + 0.5);
    print_ph_moments(mean, scv, phases, false);
}

/** A -> B -> C with Erlang activities, i.e. an SCV below one. */
void wf_erlang() {
    Wf wf("ErlangWorkflow");
    wf.add_activity("A", lang::erlang_fit_mean_order<double>(2.0, 4));
    wf.add_activity("B", lang::erlang_fit_mean_order<double>(3.0, 2));
    wf.add_activity("C", Dist::exp_mean(1.0));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.add_precedence(Wf::Serial("B", "C"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    std::printf("Erlang Workflow: A -> B -> C\n");
    std::printf("Activity A: Erlang(mean=2.0, phases=4), SCV=0.25\n");
    std::printf("Activity B: Erlang(mean=3.0, phases=2), SCV=0.50\n");
    std::printf("Activity C: Exp(mean=1.0), SCV=1.00\n");
    std::printf("Expected total mean: %.2f\n", 2.0 + 3.0 + 1.0);
    print_ph_moments(mean, scv, phases, true);
}

/** A -> B -> C with APH activities spanning three variability regimes. */
void wf_aph() {
    Wf wf("APHWorkflow");
    wf.add_activity("A", lang::aph_fit_mean_scv<double>(1.0, 0.5));
    wf.add_activity("B", lang::aph_fit_mean_scv<double>(2.0, 2.0));
    wf.add_activity("C", lang::aph_fit_mean_scv<double>(1.5, 1.0));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.add_precedence(Wf::Serial("B", "C"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    std::printf("APH Workflow: A -> B -> C\n");
    std::printf("Activity A: APH(mean=1.0, SCV=0.5)\n");
    std::printf("Activity B: APH(mean=2.0, SCV=2.0)\n");
    std::printf("Activity C: APH(mean=1.5, SCV=1.0)\n");
    std::printf("Expected total mean: %.2f\n", 1.0 + 2.0 + 1.5);
    print_ph_moments(mean, scv, phases, true);
}

/** A -> B -> [C || D] -> F -> G, serial and AND blocks combined. */
void wf_complex() {
    Wf wf("ComplexWorkflow");
    wf.add_activity("A", Dist::exp_mean(0.5));
    wf.add_activity("B", Dist::exp_mean(1.0));
    wf.add_activity("C", Dist::exp_mean(2.0));
    wf.add_activity("D", Dist::exp_mean(1.5));
    wf.add_activity("F", Dist::exp_mean(1.0));
    wf.add_activity("G", Dist::exp_mean(0.5));
    wf.add_precedence(Wf::Serial("A", "B"));
    wf.add_precedence(Wf::AndFork("B", {"C", "D"}));
    wf.add_precedence(Wf::AndJoin({"C", "D"}, "F"));
    wf.add_precedence(Wf::Serial("F", "G"));

    double mean = 0.0, scv = 0.0;
    std::size_t phases = 0;
    ph_moments(wf.to_ph(), mean, scv, phases);

    const double expected_max = exp_max_mean(2.0, 1.5);
    std::printf("Complex Workflow: A -> B -> [C || D] -> F -> G\n");
    std::printf("Activity means: A=0.5, B=1.0, C=2.0, D=1.5, F=1.0, G=0.5\n");
    std::printf("Expected max(C,D) mean: %.4f\n", expected_max);
    std::printf("Expected total mean: %.4f\n", 0.5 + 1.0 + expected_max + 1.0 + 0.5);
    print_ph_moments(mean, scv, phases, false);

    // The reference then draws 10000 samples at np.random.seed(1241) and prints
    // their moments; those are a property of that RNG stream, so the ANALYTIC
    // moments of the same law are printed here instead.
    std::printf("\nAnalytic moments of the composed law (the reference samples 10000 draws):\n");
    std::printf("Law mean: %.4f\n", mean);
    std::printf("Law std: %.4f\n", std::sqrt(scv) * mean);
    // TODO(cpp): np.random.seed(1241); samples = wf.sample(10000)
    na("Workflow.sample", "the reference prints numpy RNG-stream statistics at seed 1241, which "
                          "no other implementation reproduces digit for digit; the analytic "
                          "moments of the same law are printed above");
}

LINE_EXAMPLE("basic/workflowModels", wf_serial);
LINE_EXAMPLE("basic/workflowModels", wf_parallel);
LINE_EXAMPLE("basic/workflowModels", wf_branch);
LINE_EXAMPLE("basic/workflowModels", wf_loop);
LINE_EXAMPLE("basic/workflowModels", wf_erlang);
LINE_EXAMPLE("basic/workflowModels", wf_aph);
LINE_EXAMPLE("basic/workflowModels", wf_complex);

}  // namespace examples
}  // namespace line
