/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/initState/`: the initial state a transient starts
 * from, and the priors that select it.
 *
 * THE INITIAL STATE IS THE DECLARED ONE. A station's (statespace, stateprior)
 * pair, which `Network::set_state_prior` writes, is what CTMC, SSA (serial),
 * FLD and the JMT writer read in place of the default marking; with nothing
 * declared, every closed class starts at its reference station in phase one.
 * `initFromMarginalAndStarted` is spelled below as exactly that declaration.
 * A fluid transient takes its prior as an ODE initial condition instead
 * (`fluid_initsol_from_marginal`, `fluid_initsol_from_state_prior`).
 *
 * The transient facades are `fluid_tran_avg` and `ctmc_tran_avg`; the LDES
 * warm start is `ldes_warm_start.h`, the port of `warmStartPlacement`.
 */

#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/qn/state.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/wrappers/ldes/ldes_warm_start.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"

namespace line {
namespace examples {

namespace {

/**
 * `Network.initFromMarginalAndStarted(n, s)`: at each station, the FIRST state
 * with `n[i]` jobs of which `s[i]` are in service, declared as a one-row state
 * space with prior 1 -- the (statespace, stateprior) pair every initial-state
 * reader of this port takes over the default marking.
 */
void init_from_marginal_and_started(Net& m, const std::vector<std::vector<std::size_t> >& n,
                                    const std::vector<std::vector<std::size_t> >& s) {
    const Sn sn = m.get_struct();
    for (std::size_t i = 1; i <= sn.nstations; ++i) {
        std::vector<std::size_t> ph(sn.nclasses, 1);
        for (std::size_t r = 0; r < sn.nclasses; ++r) ph[r] = sn.phases_of(i, r + 1);
        const std::vector<std::vector<double> > rows =
            qn::from_marginal_and_started(sn, i, n[i - 1], s[i - 1], ph);
        if (rows.empty())
            throw std::invalid_argument("initFromMarginalAndStarted: station '" +
                                        sn.stations[i - 1].name + "' admits no such state");
        Matrix<double> row(1, rows[0].size());
        for (std::size_t c = 0; c < rows[0].size(); ++c) row(0, c) = rows[0][c];
        m.set_state_prior(sn.node_of_station(i), row, std::vector<double>(1, 1.0));
    }
}

/**
 * The last value and length of one (station, class) transient curve.
 *
 * `prior` NAMES THE LINE THE SHARED PARSER READS. The reference prints
 * `SteadyStateQLen[FLD/Prior1]: 4.300000` and `compare_parity.parse_cdf_...`
 * keys the golden on exactly that spelling, so a twin that printed only its own
 * `at t=end` caption reported the same number under a name nothing matches --
 * "solver FLD missing from output", a row that FAILED having compared no cell.
 * Empty leaves the extra line out.
 */
void report_curve(const TranAvg& tr, std::size_t station, std::size_t cls, std::size_t nclasses,
                  const std::string& label, const std::string& prior = std::string(),
                  const std::string& solver = "FLD") {
    const std::size_t col = station * nclasses + cls;
    std::printf("%s: %zu time points\n", label.c_str(), tr.t.size());
    if (tr.t.empty()) return;
    kv(label + " at t=0", tr.QNt.front()[col]);
    kv(label + " at t=end", tr.QNt.back()[col]);
    if (prior.empty()) return;
    std::printf("SteadyStateQLen[%s/%s]: %.6f\n", solver.c_str(), prior.c_str(),
                tr.QNt.back()[col]);
    // THE LAST POINT OF A TRANSIENT CURVE, which is one element of what
    // `getTranAvg` returns and not a result any getter reports on its own --
    // so the example, which is what selects it, declares the key. The
    // golden keys it ('Prior<k>', 'QLen') under the solver whose curve it
    // came from.
    derived(solver, prior, "QLen", tr.QNt.back()[col]);
}

/**
 * The two-node class-switching routing both `init_state_fcfs_nonexp` and
 * `init_state_ps` declare, with node 1 the delay and node 2 the queue.
 */
void link_classswitch(Net& m, std::size_t c1, std::size_t c2, std::size_t d, std::size_t q) {
    Routing P;
    P.set(c1, c1, d, d, 0.3);
    P.set(c1, c1, d, q, 0.1);
    P.set(c1, c1, q, d, 0.2);
    P.set(c1, c2, d, d, 0.6);
    P.set(c1, c2, q, d, 0.8);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c1, q, d, 1.0);
    m.link(P);
}

/** `build_model(n)` of `ldes_warmstart`: the closed near-balanced tandem. */
Net warmstart_tandem(double n) {
    Net m("ldesWarmStart");
    Delay think(m, "Think");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    ClosedClass jobs(m, "Jobs", n, think);
    think.set_service(jobs, Exp(1.0));
    q1.set_service(jobs, Exp(1.0));
    q2.set_service(jobs, Exp(0.98));  // near-balanced bottleneck
    Routing P;
    P.set(jobs, jobs, think, q1, 1.0);
    P.set(jobs, jobs, q1, q2, 1.0);
    P.set(jobs, jobs, q2, think, 1.0);
    m.link(P);
    return m;
}

}  // namespace

void init_state_fcfs_exp() {
    note("=== Initial State - FCFS Queue with Exponential Service ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue1", SchedStrategy::FCFS);
    ClosedClass k(m, "Class1", 5, q);
    m.set_service(d, k, Exp(1.0));
    m.set_service(q, k, Exp(0.7));
    link_serial(m, {d, q});

    SolverOpts o;
    o.samples = 10000;
    o.stiff = true;
    o.timespan_end = 40.0;

    note("--- Prior 1: Default initialization ---");
    note("Initial state is the default one: all 5 Class1 jobs at their reference station Queue1.");
    section("CTMC");
    report_curve(ctmc_tran_avg(m, o), 1, 0, 1, "QLen[Queue1,Class1]");
    section("FLD");
    report_curve(fluid_tran_avg(m, o), 1, 0, 1, "QLen[Queue1,Class1]", "Prior1");

    note("\n--- Prior 2: Prior on state with 3 jobs in station 2 ---");
    // `initFromMarginal([[2], [3]])`: 2 jobs at the Delay, 3 at Queue1, all in
    // phase one, which is the first moment the ODE actually needs.
    SolverOpts o2 = o;
    o2.init_sol = fluid_initsol_from_marginal(m, {{2.0}, {3.0}});
    section("FLD");
    report_curve(fluid_tran_avg(m, o2), 1, 0, 1, "QLen[Queue1,Class1]", "Prior2");

    note("\nNote: Different initial state priors lead to different transient behaviors.");
    note("      Default initialization starts with all jobs at their reference station.");
}

LINE_EXAMPLE("advanced/initState", init_state_fcfs_exp);

void init_state_fcfs_nonexp() {
    note("=== Initial State - FCFS Queue with Non-Exponential Service ===");

    Net m("model");
    Delay d(m, "Delay");
    Queue q(m, "Queue1", SchedStrategy::FCFS);
    m.set_number_of_servers(q, 3.0);
    ClosedClass k1(m, "Class1", 3, q);
    ClosedClass k2(m, "Class2", 2, q);
    m.set_service(d, k1, Exp(1.0));
    m.set_service(d, k2, Exp(1.0));
    m.set_service(q, k1, Exp(1.2));
    m.set_service(q, k2, erlang_fit(1.0, 0.5));
    link_classswitch(m, k1, k2, d, q);

    SolverOpts o;
    o.samples = 10000;
    o.stiff = true;
    o.timespan_end = 5.0;

    note("--- Prior 1: Default initialization ---");
    note("Initial state is the default one: every closed class at its reference station Queue1.");
    section("CTMC");
    report_curve(ctmc_tran_avg(m, o), 1, 0, 2, "QLen[Queue1,Class1]");
    section("FLD");
    report_curve(fluid_tran_avg(m, o), 1, 0, 2, "QLen[Queue1,Class1]", "Prior1");

    note("\n--- Prior 2: First state with marginal [0,0; 4,1] ---");
    SolverOpts o2 = o;
    o2.init_sol = fluid_initsol_from_marginal(m, {{0.0, 0.0}, {4.0, 1.0}});
    section("FLD");
    report_curve(fluid_tran_avg(m, o2), 1, 0, 2, "QLen[Queue1,Class1]", "Prior2");

    note("\n--- Prior 3: Uniform prior over states with marginal [0,0; 4,1] ---");
    SolverOpts o3 = o;
    o3.init_sol = fluid_initsol_from_state_prior(m, {{0.0, 0.0}, {4.0, 1.0}}, 2);
    section("FLD");
    report_curve(fluid_tran_avg(m, o3), 1, 0, 2, "QLen[Queue1,Class1]", "Prior3");

    note("\nNote: This example shows three types of initial state specification:");
    note("  1. Default: All jobs at reference stations");
    note("  2. Marginal: First state matching the marginal distribution");
    note("  3. Uniform: Uniform distribution over all states with given marginal");
}

LINE_EXAMPLE("advanced/initState", init_state_fcfs_nonexp);

void init_state_ps() {
    note("=== Initial State - PS Queue with Class Switching ===");
    note("This example shows solver execution on a 2-class 2-node class-switching model");
    note("with specified initial state.");

    Net m("model");
    Delay d(m, "InfiniteServer");
    Queue q(m, "Queue1", SchedStrategy::PS);
    m.set_number_of_servers(q, 2.0);
    ClosedClass k1(m, "Class1", 3, d);
    ClosedClass k2(m, "Class2", 1, d);
    m.set_service(d, k1, Exp(3.0));
    m.set_service(d, k2, Exp(0.5));
    m.set_service(q, k1, Exp(0.1));
    m.set_service(q, k2, Exp(1.0));
    link_classswitch(m, k1, k2, d, q);

    note("\nInitial state specification the reference sets:");
    note("  Marginal: [[2,1], [1,0]] - Job distribution across stations and classes");
    note("  Started:  [[0,0], [1,0]] - Number of jobs in service at each station");
    init_from_marginal_and_started(m, {{2, 1}, {1, 0}}, {{0, 0}, {1, 0}});

    SolverOpts o;
    o.seed = 23000;
    o.samples = 100000;

    section("CTMC");
    print_avg(solve_avg("CTMC", m, o));

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));

    section("SSA");
    print_avg(solve_avg("SSA", m, o));

    section("FLD");
    print_avg(solve_avg("FLD", m, o));

    section("MVA");
    print_avg(solve_avg("MVA", m, o));

    section("NC");
    print_avg(solve_avg("NC", m, o));

    note("\nNote: Initial state affects transient behavior but not steady-state metrics.");
    note("      MVA and NC compute steady-state only, so initial state has no effect.");
}

LINE_EXAMPLE("advanced/initState", init_state_ps);

void ldes_warmstart() {
    const double N = 100;
    note("=== LDES warm start from an auxiliary solver ===");

    Net exact_model = warmstart_tandem(N);
    SolverOpts o;
    o.method = "exact";
    const AvgTable exact = solve_avg("MVA", exact_model, o);
    print_vector("Exact mean queue lengths", exact.column("QLen"));

    // SolverLDES(m, SolverMVA(m, method='exact')): the warm placement is computed
    // ONCE from the auxiliary solve, and the derived init_sol reused across runs.
    const auto t0 = std::chrono::steady_clock::now();
    Net m = warmstart_tandem(N);
    mva::MvaOptions mo;
    mo.method = "exact";
    Matrix<double> none;
    const mva::AvgResult<double> aux = mva::solver_mva_run_analyzer(m.get_struct(), mo, none);
    ldes::LdesOptions warm;
    ldes::init_from_placement(warm, ldes::warm_start_placement_from_qlen(m.get_struct(), aux.QN));
    const double init_time =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    char label[64];
    std::snprintf(label, sizeof label, "Warm placement from SolverMVA (%.3fs)", init_time);
    print_vector(label, warm.init_sol);

    // SolverLDES(ms, SolverCTMC(ms, 'exact')): the mode of the exact stationary
    // law, feasible on a small instance.
    Net ms = warmstart_tandem(20);
    ctmc::CtmcOptions co;
    co.method = "exact";
    ldes::LdesOptions warm20;
    ldes::init_from_placement(warm20, ldes::warm_start_placement_from_ctmc(ms.get_struct(), co));
    print_vector("CTMC distribution-mode placement (N=20)", warm20.init_sol);

    const std::vector<double> ex = exact.column("QLen");
    double exsum = 0.0;
    for (std::size_t i = 0; i < ex.size(); ++i) exsum += ex[i];
    // run_set: the average L1 relative error and the total wall clock over the seeds.
    const long seeds[3] = {23000, 23001, 23002};
    auto run_set = [&](std::size_t samples, bool warm_start, double& rtime) {
        double err = 0.0;
        rtime = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            Net mk = warmstart_tandem(N);
            const auto s0 = std::chrono::steady_clock::now();
            ldes::LdesOptions lo = warm_start ? warm : ldes::LdesOptions();
            lo.samples = samples;
            lo.seed = seeds[k];
            const ldes::LdesResult r = ldes::solver_ldes(mk.get_struct(), lo);
            rtime += std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
            double e = 0.0;
            for (std::size_t i = 0; i < ex.size(); ++i) e += std::fabs(r.QN(i, 0) - ex[i]);
            err += e / exsum;
        }
        return err / 3.0;
    };

    std::printf("\n%10s  %16s  %16s\n", "samples", "COLD err/time", "WARM err/time");
    const std::size_t budgets[3] = {10000, 50000, 200000};
    std::size_t cold_at = 0, warm_at = 0;
    for (std::size_t k = 0; k < 3; ++k) {
        double ct = 0.0, wt = 0.0;
        const double ce = run_set(budgets[k], false, ct);
        const double we = run_set(budgets[k], true, wt);
        if (cold_at == 0 && ce < 0.10) cold_at = budgets[k];
        if (warm_at == 0 && we < 0.10) warm_at = budgets[k];
        std::printf("%10zu  %7.2f%% %6.2fs  %7.2f%% %6.2fs\n", budgets[k], 100 * ce, ct, 100 * we,
                    wt);
        if (cold_at != 0 && warm_at != 0) break;
    }
    auto at = [](std::size_t b) { return b ? std::to_string(b) : std::string("None"); };
    std::printf("\nSamples to reach 10%% error: COLD=%s, WARM=%s\n", at(cold_at).c_str(),
                at(warm_at).c_str());
}

LINE_EXAMPLE("advanced/initState", ldes_warmstart);

}  // namespace examples
}  // namespace line
