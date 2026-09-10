/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/initState/`: the initial state a transient starts
 * from, and the priors that select it.
 *
 * THIS PORT HAS NO `sn.state`. `fluid_closing.h` records the design decision in
 * full: the NetworkStruct carries no state field, and every solver seeds itself
 * from the state `initDefault` WOULD have written -- every closed class at its
 * reference station, in phase one. `initFromMarginal`,
 * `initFromMarginalAndStarted` and `setStatePrior` on a station therefore have
 * no counterpart, so the blocks that change the prior are refused by name
 * rather than answered from the default state under the prior's caption.
 *
 * The C++ transient facade is `fluid_tran_avg` alone. `SolverCTMC.getTranAvg`
 * exists in the library (`ctmc::solver_ctmc_transient_analyzer`) but has no
 * entry point in `example_util.h`, so the CTMC transient blocks name that gap
 * instead of quietly reporting the fluid curve.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"

namespace line {
namespace examples {

namespace {

/** The refusal the changed-prior blocks share, named after the absent method. */
void na_prior(const std::string& method, const std::string& what) {
    na(method, "this port's NetworkStruct carries no `state` field, so " + what +
                   " has no counterpart; every solve here starts from the default initial state "
                   "(each closed class at its reference station, in phase one)");
}

/** The CTMC transient gap, named where the reference calls getTranAvg. */
void na_ctmc_tran() {
    na("CTMC",
       "getTranAvg has no entry point in example_util.h; the library carries "
       "ctmc::solver_ctmc_transient_analyzer but the example facade exposes only the fluid "
       "transient");
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
                  const std::string& label, const std::string& prior = std::string()) {
    const std::size_t col = station * nclasses + cls;
    std::printf("%s: %zu time points\n", label.c_str(), tr.t.size());
    if (tr.t.empty()) return;
    kv(label + " at t=0", tr.QNt.front()[col]);
    kv(label + " at t=end", tr.QNt.back()[col]);
    if (!prior.empty())
        std::printf("SteadyStateQLen[FLD/%s]: %.6f\n", prior.c_str(), tr.QNt.back()[col]);
        // THE LAST POINT OF A TRANSIENT CURVE, which is one element of what
        // `getTranAvg` returns and not a result any getter reports on its own --
        // so the example, which is what selects it, declares the key. The
        // golden keys it ('Prior<k>', 'QLen') under the solver whose curve it
        // came from.
        derived("FLD", prior, "QLen", tr.QNt.back()[col]);
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
    // TODO(cpp): model.initDefault(); QNt_ctmc, _, _ = CTMC(model, options).getTranAvg(Qt, Ut, Tt)
    na_ctmc_tran();
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
    // TODO(cpp): model.initDefault(); QNt_ctmc_1, _, _ = CTMC(model, options).getTranAvg(Qt, Ut, Tt)
    na_ctmc_tran();
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
    // TODO(cpp): model.initFromMarginalAndStarted([[2, 1], [1, 0]], [[0, 0], [1, 0]])
    na_prior("Network.initFromMarginalAndStarted",
             "initFromMarginalAndStarted([[2,1], [1,0]], [[0,0], [1,0]])");
    note("The blocks below therefore run from the default initial state. That leaves the CTMC,");
    note("MVA and NC answers unchanged (they are steady-state on an irreducible chain) and makes");
    note("the SSA and FLD blocks start from a different point of the same model.");

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

    // TODO(cpp): proto = SolverLDES(m, SolverMVA(m, method='exact')); init_sol = proto.options.init_sol
    // TODO(cpp): proto_ctmc = SolverLDES(ms, SolverCTMC(ms, 'exact', force=True))
    // TODO(cpp): run_set(N, exact, b, seeds, None) and run_set(N, exact, b, seeds, init_sol)
    // The LDES engine runs here (`ldes_avg`); what this facade has no entry
    // point for is the WARM START -- `SolverLDES(model, SolverMVA)` writing
    // options.init_sol, which the CLI reaches as `--ldes-initsol`.
    na("LDES",
       "the warm placement is SolverLDES(model, SolverMVA) writing options.init_sol, and this "
       "example facade exposes no init_sol argument; the engine itself runs (see the other LDES "
       "blocks) and `line-cli -s ldes --ldes-initsol` takes the placement");
}

LINE_EXAMPLE("advanced/initState", ldes_warmstart);

}  // namespace examples
}  // namespace line
