/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `advanced/randomEnv`: a queueing network modulated by a random environment.
 *
 * Each stage of the environment holds a COMPLETE network, and what couples the
 * stages is that the jobs present at a switch are carried into the next stage.
 * `env::solver_env` is the analyzer selection of `SolverENV.init`, and the
 * examples here drive its mean-field arm with fluid stages, which is the one
 * coupling this port carries.
 *
 * THREE THINGS THE REFERENCE ASKS FOR THAT ARE REFUSED BY NAME rather than
 * approximated, because each would answer a different question under the same
 * label:
 *
 *  - `options.method = 'default'` is the reference's spelling of the DEFAULT
 *    coupling, which is the mean-field one; `EnvOptions::method` names the
 *    coupling directly, so `meanfield` is the same request and not a
 *    substitution. It is the only renaming below.
 *  - `options.timespan = [0, Inf]`. The mean-field exit metrics are a
 *    Riemann-Stieltjes sum of the stage trajectory against the holding-time
 *    CDF, so an infinite horizon has no grid to sum over and `SolverEnv::init`
 *    refuses it. Where the reference sets a finite stage horizon that value is
 *    used; where it sets Inf the port's default horizon is used and printed on
 *    the banner, so the quadrature the number came from is stated.
 *  - `ENV(env, @(m) CTMC(m, 'exact', ...))` -- a TRANSIENT CTMC stage solve.
 *    The mean-field coupling in this port takes fluid stages only, and the
 *    state-vector coupling that does take CTMC stages is a DIFFERENT coupling,
 *    so running it under a `method='default'` example would report one
 *    analyzer's numbers under another's name.
 */

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/qn/environment.h"
#include "line/solvers/env/env_dispatch.h"
#include "line/io/map2renv.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

using Env = env::Environment<double>;

/** `envModel.getStageTable()`: the stage index, its name and its category. */
void stage_table(const Env& e) {
    std::printf("%-8s %-16s %-14s %-16s\n", "Stage", "Name", "Type", "Model");
    for (std::size_t s = 0; s < e.nstages(); ++s)
        // A LAYERED stage carries no NetworkStruct and so has no model name of
        // its own; it is named by what it is, as the reference's stage table
        // names a LayeredNetwork by its class rather than by a station list.
        std::printf("%-8zu %-16s %-14s %-16s\n", s + 1, e.stage(s).name.c_str(),
                    e.stage(s).type.c_str(),
                    e.is_lqn(s) ? "LayeredNetwork" : e.stage(s).model.name.c_str());
}

/**
 * The environment-blended AvgTable.
 *
 * RespT and ArvR are NaN and ResidT is QLen/Tput, which is exactly what
 * `@SolverENV/getEnsembleAvg` returns: ENV blends per-stage metrics over the
 * environment process and computes no response time or arrival rate at all.
 * Printing QLen/Tput under RespT would invent a number the reference declines
 * to give.
 */
void print_env_avg(const Sn& sn, const env::EnvAnalyzerSolution<double>& r) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    avg_rows(sn, [&](std::size_t i, std::size_t c, int k) {
        switch (k) {
            case 0: return r.QN(i, c);
            case 1: return r.UN(i, c);
            case 2: return nan;
            case 3: return r.QN(i, c) / r.TN(i, c);
            case 4: return nan;
            default: return r.TN(i, c);
        }
    });
}

/** The banner the CLI's `-s env` arm prints: the quadrature is part of the answer. */
void env_banner(const Env& e, const env::EnvOptions& o,
                const env::EnvAnalyzerSolution<double>& r) {
    std::printf("SolverENV method=%s stages=%zu horizon=%.6g points=%zu iters=%d%s\n",
                r.method.c_str(), e.nstages(), o.timespan_end, o.tran_points, r.iterations,
                r.converged ? "" : " (NOT CONVERGED)");
}

/** `renv_genqn.m`: a Delay and a PS Queue in a cycle, one closed class of N jobs. */
Net renv_genqn_model(double rate_delay, double rate_queue, double N) {
    Net m("qn1");
    Delay d(m, "Queue1");
    Queue q(m, "Queue2", SchedStrategy::PS);
    ClosedClass c(m, "Class1", N, d);
    d.set_service(c, Exp(rate_delay));
    q.set_service(c, Exp(rate_queue));
    Routing P;
    cyclic(P, c, {d, q});
    m.link(P);
    return m;
}

/**
 * The renv_basic base network with the server running at `svc`.
 *
 * The model NAME stays `BaseModel` in both stages: the reference builds one
 * network and hands each stage a `copy()` of it with the service overwritten,
 * and a copy keeps the name it was copied from.
 */
Net renv_basic_stage(double svc) {
    Net m("BaseModel");
    Delay d(m, "ThinkTime");
    Queue q(m, "Fast/Slow Server", SchedStrategy::FCFS);
    ClosedClass c(m, "Jobs", 5, d);
    d.set_service(c, Exp(1.0));
    q.set_service(c, Exp(svc));
    Routing P;
    cyclic(P, c, {d, q});
    m.link(P);
    return m;
}

/** The renv_node_breakdown base model: an open M/M/1 whose server can fail. */
Net renv_breakdown_base() {
    Net m("ServerWithFailures");
    Source src(m, "Arrivals");
    Queue q(m, "Server", SchedStrategy::FCFS);
    Sink snk(m, "Departures");
    OpenClass c(m, "Jobs");
    src.set_arrival(c, Exp(0.8));
    q.set_service(c, Exp(2.0));
    q.set_number_of_servers(1.0);
    Routing P;
    serial(P, c, {src, q, snk});
    m.link(P);
    return m;
}

/** `Environment.addNodeFailureRepair` on the fixed two-stage slot pair. */
Env renv_breakdown_env(const std::string& nm, const Net& base) {
    Env e(nm, 2);
    Net b = base;
    e.add_node_failure_repair(0, 1, b.get_struct(), "Server", Exp(0.1), Exp(1.0),
                             Exp(0.5));
    return e;
}

}  // namespace

// ---------------------------------------------------------------------------
// renv_basic
// ---------------------------------------------------------------------------

/**
 * A server alternating between a Fast and a Slow mode, solved by ENV over fluid
 * stages and compared against the STEADY state of each stage taken alone.
 */
void renv_basic() {
    Net fast = renv_basic_stage(4.0);
    Net slow = renv_basic_stage(1.0);

    Env e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, Exp(0.5));
    e.add_transition(1, 0, Exp(1.0));

    note("Environment stages:");
    stage_table(e);

    env::EnvOptions o;
    o.method = "meanfield";
    o.timespan_end = 100.0;
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);

    note("\n--- Environment-Averaged Results ---");
    env_banner(e, o, r);
    // The reference prints these two blocks with no solver banner, so the
    // declaration is made to the recorder alone. An ENSEMBLE IS NAMED BY ITS
    // MEMBER: `ENV(FLD)` is a mean-field coupling over fluid stages, which is
    // what the golden's `FLD` row holds and what the comparator reconciles the
    // two spellings of.
    attribute("ENV(FLD)");
    print_env_avg(e.stage(0).model, r);

    note("\n--- Individual Stage Analysis (MVA) ---");
    for (std::size_t s = 0; s < e.nstages(); ++s) {
        std::printf("\nStage %zu:\n", s);
        const Sn& sn = e.stage(s).model;
        mva::MvaOptions mo;
        Matrix<double> init;
        // One table per stage under one key; the golden holds the FIRST, which
        // is why `renv_basic` is listed in `corpus.json`'s multi-model set.
        attribute("MVA");
        print_avg(sn, mva::solver_mva_run_analyzer(sn, mo, init));
    }
}

// ---------------------------------------------------------------------------
// renv_twostages_repairmen
// ---------------------------------------------------------------------------

/** Two stages, UP and DOWN, over a closed Delay/PS pair holding one job. */
void renv_twostages_repairmen() {
    const std::size_t E = 2;
    const double N = 1.0;
    const double rate[2][2] = {{2.0, 1.0}, {1.0, 2.0}};
    const char* env_name[2] = {"Stage1", "Stage2"};
    const char* env_type[2] = {"UP", "DOWN"};

    Env e("MyEnv", E);
    std::vector<Net> sub;
    for (std::size_t s = 0; s < E; ++s) sub.push_back(renv_genqn_model(rate[0][s], rate[1][s], N));
    for (std::size_t s = 0; s < E; ++s)
        e.set_stage(s, env_name[s], env_type[s], sub[s].get_struct());

    const double env_rates[2][2] = {{0.0, 1.0}, {0.5, 0.5}};
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b)
            if (env_rates[a][b] > 0.0) e.add_transition(a, b, Exp(env_rates[a][b]));

    note("Stage Table:");
    stage_table(e);

    env::EnvOptions o;
    o.method = "meanfield";
    o.iter_max = 100;
    o.iter_tol = 0.01;
    // The reference's ENV timespan is [0, Inf]; the horizon that reaches a
    // stage is the fluid solver's own, which it sets to 1e3.
    o.timespan_end = 1e3;
    // The reference catches the solve and PRINTS the failure rather than
    // stopping, so a stage the coupling cannot take costs this block and not
    // every block after it.
    try {
        const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
        section("ENV");
        env_banner(e, o, r);
        const Sn& sn = e.stage(0).model;
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                std::printf("QN[%s,%s] = %.6g  UN = %.6g  TN = %.6g\n",
                            sn.stations[i].name.c_str(), sn.classes[c].name.c_str(), r.QN(i, c),
                            r.UN(i, c), r.TN(i, c));
            }
        note("\nAverage Table:");
        print_env_avg(sn, r);
    } catch (const std::exception& ex) {
        std::printf("Error during solving: %s\n", ex.what());
        note("Note: Environment solver features may not be fully implemented");
    }

    // NOT a gap in this port: the reference leaves its CTMC alternative
    // commented out, as the MATLAB script it follows does, so there is no call
    // to record and nothing was refused.
    note("\n## Alternative: CTMC Solver (Commented)");
    note("The MATLAB version also shows an alternative using CTMC, which is commented out in "
         "the original.");
}

// ---------------------------------------------------------------------------
// renv_threestages_repairmen
// ---------------------------------------------------------------------------

/**
 * Three stages on a circulant transition graph with Erlang holding times.
 *
 * The reference's point is `envSolver.getGenerator()`, the infinitesimal
 * generator of the joint (environment, network) chain, which it obtains from
 * CTMC stage solvers. Neither piece exists on this path, so both are refused by
 * name and the model itself is still built and reported.
 */
void renv_threestages_repairmen() {
    const std::size_t E = 3;
    const double N = 2.0;
    const char* env_name[3] = {"Stage1", "Stage2", "Stage3"};
    const char* env_type[3] = {"UP", "DOWN", "FAST"};

    // rate = ones(M,E); rate(M,:) = 1:E; rate(1,:) = E:-1:1, with M = 2.
    const double rate_delay[3] = {3.0, 2.0, 1.0};
    const double rate_queue[3] = {1.0, 2.0, 3.0};

    note("Rate matrix:");
    for (std::size_t s = 0; s < E; ++s)
        std::printf("  stage %zu: Queue1 %.6g  Queue2 %.6g\n", s + 1, rate_delay[s],
                    rate_queue[s]);

    Env e("MyEnv", E);
    std::vector<Net> sub;
    for (std::size_t s = 0; s < E; ++s) sub.push_back(renv_genqn_model(rate_delay[s], rate_queue[s], N));
    for (std::size_t s = 0; s < E; ++s)
        e.set_stage(s, env_name[s], env_type[s], sub[s].get_struct());

    // circul(3) is the right circular shift squared: 1 -> 2 -> 3 -> 1 at rate 1.
    const double env_rates[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}};
    note("Environment transition rates (circulant matrix):");
    for (std::size_t a = 0; a < E; ++a)
        std::printf("  %.6g %.6g %.6g\n", env_rates[a][0], env_rates[a][1], env_rates[a][2]);
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b)
            if (env_rates[a][b] > 0.0)
                e.add_transition(a, b, lang::erlang_fit_mean_order<double>(1.0 / env_rates[a][b],
                                                                          a + b + 2));

    note("The metasolver considers an environment with 3 stages and a queueing network with 2 "
         "stations.");
    note("This example illustrates the computation of the infinitesimal generator of the system.");
    stage_table(e);

    // TODO(cpp): solvers = [CTMC(envSubModel[e]) for e in range(E)]
    // TODO(cpp): envSolver = ENV(envModel, solvers)
    // TODO(cpp): infGen, stageInfGen = envSolver.generator(); print(infGen)
    na("SolverENV.getGenerator",
       "the joint (environment, network) infinitesimal generator and the per-stage generators "
       "are read off enumerated CTMC stage solves; this port carries no getGenerator on the ENV "
       "path");
    // `ENV(envModel, [CTMC(sub) for sub in stages])`: the reference's stage
    // solver is the exact chain, under the DEFAULT (mean-field) coupling.
    env::EnvOptions o;
    o.method = "meanfield";
    o.stage_solver = "ctmc";
    o.iter_max = 100;
    o.iter_tol = 0.01;
    o.timespan_end = 1e3;
    // No catch here: a failed solve is this example's result, and swallowing it
    // into a printed line would report `ok` with no table. Letting it throw is
    // what `line-examples` reports as FAILED and the harvest records as an error.
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
    section("ENV(CTMC)", "ENV");
    env_banner(e, o, r);
    print_env_avg(e.stage(0).model, r);
}

// ---------------------------------------------------------------------------
// renv_fourstages_repairmen
// ---------------------------------------------------------------------------

/** Four stages with APH holding times, over the same closed Delay/PS pair. */
void renv_fourstages_repairmen() {
    const std::size_t E = 4;
    const double N = 30.0;
    const char* env_name[4] = {"Stage1", "Stage2", "Stage3", "Stage4"};
    const char* env_type[4] = {"UP", "DOWN", "FAST", "SLOW"};

    // rate = ones(M,E); rate(M,:) = 1:E; rate(1,:) = E:-1:1, with M = 3: the
    // generator reads rows 1 and 2, so the queue runs at 1 in every stage.
    const double rate_delay[4] = {4.0, 3.0, 2.0, 1.0};
    const double rate_queue[4] = {1.0, 1.0, 1.0, 1.0};

    note("Rate matrix:");
    for (std::size_t s = 0; s < E; ++s)
        std::printf("  stage %zu: Queue1 %.6g  Queue2 %.6g\n", s + 1, rate_delay[s],
                    rate_queue[s]);

    Env e("MyEnv", E);
    std::vector<Net> sub;
    for (std::size_t s = 0; s < E; ++s) sub.push_back(renv_genqn_model(rate_delay[s], rate_queue[s], N));
    for (std::size_t s = 0; s < E; ++s)
        e.set_stage(s, env_name[s], env_type[s], sub[s].get_struct());

    const double env_rates[4][4] = {{0.0, 0.5, 0.0, 0.0},
                                    {0.0, 0.0, 0.5, 0.5},
                                    {0.5, 0.0, 0.0, 0.5},
                                    {0.5, 0.5, 0.0, 0.0}};
    note("Environment transition rates:");
    for (std::size_t a = 0; a < E; ++a)
        std::printf("  %.6g %.6g %.6g %.6g\n", env_rates[a][0], env_rates[a][1], env_rates[a][2],
                    env_rates[a][3]);
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b)
            if (env_rates[a][b] > 0.0)
                e.add_transition(a, b,
                                 lang::aph_fit_mean_scv<double>(1.0 / env_rates[a][b], 0.5));

    note("The metasolver considers an environment with 4 stages and a queueing network with 3 "
         "stations.");
    note("Every time the stage changes, the queueing network will modify the service rates of "
         "the stations.");
    stage_table(e);

    env::EnvOptions o;
    o.method = "meanfield";
    o.iter_max = 100;
    o.iter_tol = 0.05;
    // The reference sets both the ENV and the fluid horizon to Inf, which the
    // Riemann-Stieltjes exit quadrature has no grid for; the banner states the
    // finite horizon the numbers below were summed over.
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);

    note("AvgTable =");
    env_banner(e, o, r);
    attribute("ENV(FLD)");
    print_env_avg(e.stage(0).model, r);
}

// ---------------------------------------------------------------------------
// renv_node_breakdown
// ---------------------------------------------------------------------------

/**
 * The node breakdown / repair convenience API, in the reference's five blocks.
 *
 * `addNodeBreakdown` grows the reference's stage graph as breakdowns are
 * declared; this port sizes the environment at construction, so the UP and DOWN
 * slots are named at the call instead. Everything else -- the degraded copy of
 * the model, the two arcs, the pair of reset policies -- is the reference's.
 */
void renv_node_breakdown_blocks() {
    const Net base = renv_breakdown_base();

    note("======================================================================");
    note("Example 1: Using addNodeFailureRepair with node object");
    note("======================================================================");
    Env e1 = renv_breakdown_env("ServerEnv1", base);
    e1.init();
    note("\nStage table for env1:");
    stage_table(e1);

    note("\n======================================================================");
    note("Example 2: Using separate breakdown and repair calls with node object");
    note("======================================================================");
    Env e2("ServerEnv2", 2);
    {
        Net b = base;
        e2.add_node_breakdown(0, 1, b.get_struct(), "Server", Exp(0.1), Exp(0.5));
    }
    e2.add_node_repair("Server", Exp(1.0));
    e2.init();
    note("\nStage table for env2:");
    stage_table(e2);

    note("\n======================================================================");
    note("Example 3: With custom reset policies using node object");
    note("======================================================================");
    Env e3("ServerEnv3", 2);
    {
        Net b = base;
        e3.add_node_failure_repair(0, 1, b.get_struct(), "Server", Exp(0.1),
                                   Exp(1.0), Exp(0.5), "clear", "keep");
    }
    e3.init();
    note("\nStage table for env3:");
    stage_table(e3);

    note("\n======================================================================");
    note("Example 4: Modifying reset policies after creation using node object");
    note("======================================================================");
    Env e4 = renv_breakdown_env("ServerEnv4", base);
    // The breakdown arc flushes the buffer, the repair arc carries it over;
    // `keep` is the EMPTY reset, which is how this port spells the identity.
    e4.set_reset(0, 1, env::env_reset_policy("clear"));
    e4.set_reset(1, 0, env::env_reset_policy("keep"));
    e4.init();
    note("\nStage table for env4:");
    stage_table(e4);

    note("\n======================================================================");
    note("Example 5: Solving environment model with ENV solver");
    note("======================================================================");
    Env e5 = renv_breakdown_env("ServerEnv5", base);
    e5.init();

    env::EnvOptions o;
    o.method = "meanfield";
    o.iter_tol = 0.01;
    o.iter_max = 100;
    o.timespan_end = 1000.0;
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e5, o);

    note("\nAverage Performance Metrics:");
    attribute("ENV(FLD)");
    env_banner(e5, o, r);
    print_env_avg(e5.stage(0).model, r);

    note("\nInterpretation:");
    note("- The system alternates between UP (operational) and DOWN (failed) states");
    note("- UP state: Server processes jobs at rate 2.0");
    note("- DOWN state: Server processes jobs at reduced rate 0.5");
    note("- Breakdown occurs at rate 0.1 (mean time to failure = 10 time units)");
    note("- Repair occurs at rate 1.0 (mean time to repair = 1 time unit)");
    note("- Results show averaged performance across both states");
    note("");

    note("======================================================================");
    note("All examples completed successfully!");
    note("======================================================================");
}

/**
 * The reference runs the five blocks inside ONE try/except and prints the
 * failure; without it the first block that throws would silently drop the four
 * below it and the closing banner.
 */
void renv_node_breakdown() {
    try {
        renv_node_breakdown_blocks();
    } catch (const std::exception& ex) {
        std::printf("Error running examples: %s\n", ex.what());
    }
}

// ---------------------------------------------------------------------------
// renv_lqn_twostages
// ---------------------------------------------------------------------------

namespace {

/**
 * The LQN both stages are built from: a client task calling a database task.
 *
 * `db_mean` is the whole difference between the stages -- the DB activity's host
 * demand -- so the two stages are structurally identical and their layer blocks
 * line up entry for entry, which is what the coupling blends across.
 */
lqn::LqnStruct<double> renv_lqn_model(double db_mean) {
    lqn::LqnBuilder<double> b;
    b.processor("ClientProcessor", 1, SchedStrategy::PS);
    b.processor("DBProcessor", 1, SchedStrategy::PS);
    b.task("ClientTask", 5, SchedStrategy::REF, "ClientProcessor");
    b.think_time("ClientTask", D::exp_mean(5.0));
    b.task("DBTask", std::numeric_limits<double>::infinity(), SchedStrategy::INF, "DBProcessor");
    b.entry("ClientEntry", "ClientTask");
    b.entry("DBEntry", "DBTask");
    b.activity("ClientActivity", D::exp_mean(1.0), "ClientTask");
    b.bound_to("ClientActivity", "ClientEntry");
    b.sync_call("ClientActivity", "DBEntry", 2.5);
    b.activity("DBActivity", D::exp_mean(db_mean), "DBTask");
    b.bound_to("DBActivity", "DBEntry");
    b.replies_to("DBActivity", "DBEntry");
    return b.build();
}

/** The two-stage UP/DOWN environment over LQN stages, at the given switch rates. */
Env renv_lqn_env(const std::string& nm, double a, double b) {
    Env e(nm, 2);
    e.set_lqn_stage(0, "UP", "operational", renv_lqn_model(0.8));
    e.set_lqn_stage(1, "DOWN", "degraded", renv_lqn_model(3.0));
    e.add_transition(0, 1, Exp(a));
    e.add_transition(1, 0, Exp(b));
    return e;
}

/** The EnvOptions the reference's `lnFactory` amounts to: fluid layers over [0, T]. */
env::EnvOptions renv_lqn_options(double T, int iter_max, double iter_tol) {
    env::EnvOptions o;
    o.method = "meanfield";
    o.timespan_end = T;
    o.iter_max = iter_max;
    o.iter_tol = iter_tol;
    o.lqn.layer_solver = "fluid";
    return o;
}

/**
 * `aggregateStationClassNames`: the labels of the block-diagonal aggregate.
 *
 * A layered stage has no stations of its own, so the aggregate row and column
 * names are the LAYER's, prefixed by the layer name to disambiguate the client
 * Delay that every layer carries.
 */
void lqn_aggregate_names(const ln::SolverLN<double>& s, std::vector<std::string>& rows,
                         std::vector<std::string>& cols) {
    const ln::LnLayerBlocks b = s.layer_blocks();
    rows.assign(b.M, std::string());
    cols.assign(b.K, std::string());
    const std::vector<qn::Layer<double>>& L = s.layers();
    for (std::size_t l = 0; l < L.size(); ++l) {
        for (std::size_t i = 0; i < b.msz[l]; ++i)
            rows[b.roff[l] + i] = L[l].name + "." + L[l].stations[i].name;
        for (std::size_t r = 0; r < b.ksz[l]; ++r)
            cols[b.coff[l] + r] = L[l].name + "." + L[l].classes[r].name;
    }
}

/**
 * The environment-blended AvgTable of a LAYERED environment.
 *
 * `avg_rows` cannot serve it: it walks the stations of a `NetworkStruct`, which
 * a layered stage does not have. The columns are the flat table's and carry the
 * same NaNs -- ENV computes no response time and no arrival rate -- and the
 * off-diagonal (station, class) pairs, which pair one layer's station with
 * another layer's class and stand for nothing, are dropped exactly as the
 * reference drops an unvisited pair.
 */
void print_env_avg_lqn(const ln::SolverLN<double>& s,
                       const env::EnvAnalyzerSolution<double>& r) {
    std::vector<std::string> rows, cols;
    lqn_aggregate_names(s, rows, cols);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    avg_header();
    namespace parity = line::examples::parity;
    if (parity::enabled()) parity::begin_table("avg", "Station", "JobClass");
    for (std::size_t i = 0; i < rows.size() && i < r.QN.rows(); ++i)
        for (std::size_t c = 0; c < cols.size() && c < r.QN.cols(); ++c) {
            const double q = r.QN(i, c), u = r.UN(i, c), t = r.TN(i, c);
            if (q == 0.0 && u == 0.0 && t == 0.0) continue;
            const double w = q / t;
            std::printf("%-16s %-14s %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",
                        rows[i].c_str(), cols[c].c_str(), q, u, nan, w, nan, t);
            if (!parity::enabled()) continue;
            std::vector<parity::Cell> cells;
            cells.push_back(parity::Cell{"QLen", q});
            cells.push_back(parity::Cell{"Util", u});
            cells.push_back(parity::Cell{"RespT", nan});
            cells.push_back(parity::Cell{"ResidT", w});
            cells.push_back(parity::Cell{"ArvR", nan});
            cells.push_back(parity::Cell{"Tput", t});
            parity::add_row(rows[i], cols[c], cells);
        }
}

/** The finite total of a blended metric, the reference's `sum(TN(isfinite(TN)))`. */
double finite_sum(const Matrix<double>& A) {
    double s = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            if (std::isfinite(A(i, j))) s += A(i, j);
    return s;
}

/** `envTputSum`: solve one LQN-in-ENV and return its finite aggregate throughput. */
double renv_lqn_tput(Env& e, double T) {
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, renv_lqn_options(T, 10, 0.03));
    return finite_sum(r.TN);
}

/**
 * `stageAggregate`: one LQN stage alone, in the same block-diagonal layout.
 *
 * The last point of the layered transient over the same horizon, which is what
 * the coupling would carry into the next stage had there been one -- so the two
 * are read on the same terms.
 */
void renv_lqn_stage_alone(const lqn::LqnStruct<double>& model, double T, double& sumQ,
                          double& sumT) {
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = T;
    ln::SolverLN<double> s(model, lo);
    const ln::LnTranSolution tr = s.get_tran_avg();
    sumQ = 0.0;
    sumT = 0.0;
    for (std::size_t l = 0; l < tr.layers.size(); ++l) {
        const ln::LnTranLayer& L = tr.layers[l];
        if (L.t.empty()) continue;
        const std::size_t last = L.t.size() - 1;
        for (std::size_t i = 0; i < L.QN.size(); ++i)
            for (std::size_t r = 0; r < L.QN[i].size(); ++r) {
                if (last < L.QN[i][r].size()) sumQ += L.QN[i][r][last];
                if (last < L.TN[i][r].size()) sumT += L.TN[i][r][last];
            }
    }
}

}  // namespace

/**
 * A LayeredNetwork operating in a two-stage random environment.
 *
 * WHAT MAKES THIS DIFFERENT FROM EVERY OTHER EXAMPLE HERE is only the stage
 * model: the coupling is the same mean-field fixed point, carrying the same
 * marginal mean queue lengths across a switch. A layered stage has no stations
 * of its own, so the (station, class) view the metrics are blended over is the
 * BLOCK-DIAGONAL UNION of the layers SolverLN builds -- `layerBlocks` in the
 * reference -- and the handoff splits an aggregate marginal into per-layer
 * blocks on the way in and reassembles the layer transients on the way out.
 *
 * THE ORACLE IS THE BRACKET, and it is the reference's own: a slower database
 * slows the whole system, so the environment-averaged total throughput must lie
 * between the two single-stage answers, and must RISE with the stationary
 * probability of the fast stage. That is a property no quadrature choice can
 * fake, which is why the reference asserts it rather than a number: a handoff
 * that is silently inert restarts every stage from the default state, reports
 * warm-up averages, and lands OUTSIDE the bracket.
 *
 * `options.method = 'default'` is again the reference's spelling of this
 * coupling; `meanfield` is the same request under the name this port gives it.
 */
void renv_lqn_twostages() {
    const double T = 50.0;

    Env e = renv_lqn_env("DBReliability", 0.2, 1.0);
    note("Environment stages:");
    stage_table(e);

    const env::EnvOptions o = renv_lqn_options(T, 20, 0.02);
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);

    note("\n--- Environment-Averaged Results ---");
    env_banner(e, o, r);
    std::printf("ENV over LQN ran: aggregate size = %zu x %zu\n", r.QN.rows(), r.QN.cols());
    // The labels come from a solver built on the SAME stage model, so they are
    // the labels of the very blocks the aggregate was assembled from.
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = T;
    ln::SolverLN<double> labels(renv_lqn_model(0.8), lo);
    attribute("ENV(LN(FLD))");
    print_env_avg_lqn(labels, r);

    // (1) Finite, and the closed population is conserved: the aggregate sums
    // each layer's own five client jobs, and the environment moves none of them.
    double upQ = 0.0, upT = 0.0, dnQ = 0.0, dnT = 0.0;
    renv_lqn_stage_alone(renv_lqn_model(0.8), T, upQ, upT);
    renv_lqn_stage_alone(renv_lqn_model(3.0), T, dnQ, dnT);
    const double envQ = finite_sum(r.QN), envT = finite_sum(r.TN);
    note("\n--- Single-Stage References (same block-diagonal layout) ---");
    std::printf("Total throughput  UP=%.4f  DOWN=%.4f  ENV=%.4f\n", upT, dnT, envT);
    std::printf("Aggregate Q (sum) UP=%.4f  DOWN=%.4f  ENV=%.4f\n", upQ, dnQ, envQ);
    derived_here("Aggregate", "SumQ", envQ);
    derived_here("Aggregate", "SumT", envT, "Tput");

    const double lo_x = std::min(upT, dnT), hi_x = std::max(upT, dnT);
    const double tol = 1e-2 * std::max(1.0, hi_x);
    const bool conserved = std::fabs(envQ - upQ) < 1e-2 && std::fabs(envQ - dnQ) < 1e-2;
    const bool bracketed = envT >= lo_x - tol && envT <= hi_x + tol;

    // (2) Monotone in P(UP) = b/(a+b): more time in the fast stage, more work
    // done. Both settings must stay inside the same bracket.
    Env e_low = renv_lqn_env("R", 1.0, 0.2);   // P(UP) = 0.167
    Env e_high = renv_lqn_env("R", 0.2, 1.0);  // P(UP) = 0.833
    const double x_low = renv_lqn_tput(e_low, T);
    const double x_high = renv_lqn_tput(e_high, T);
    std::printf("Monotonicity   P(UP)=0.167 -> %.4f   P(UP)=0.833 -> %.4f\n", x_low, x_high);
    derived_here("Monotone", "Xlow", x_low, "Tput");
    derived_here("Monotone", "Xhigh", x_high, "Tput");
    const bool monotone = x_high > x_low + 1e-3 && x_low >= lo_x - tol && x_high <= hi_x + tol;

    // (3) A three-stage environment, which exercises the E > 2 coupling: the
    // entry marginal of MID mixes the exits of both its neighbours by probOrig.
    Env e3("R3", 3);
    e3.set_lqn_stage(0, "UP", "operational", renv_lqn_model(0.8));
    e3.set_lqn_stage(1, "MID", "degraded", renv_lqn_model(1.6));
    e3.set_lqn_stage(2, "DOWN", "failed", renv_lqn_model(3.0));
    e3.add_transition(0, 1, Exp(0.3));
    e3.add_transition(1, 2, Exp(0.3));
    e3.add_transition(2, 1, Exp(0.6));
    e3.add_transition(1, 0, Exp(0.6));
    const double x3 = renv_lqn_tput(e3, T);
    std::printf("Three-stage ENV throughput = %.4f (bracket [%.4f, %.4f])\n", x3, lo_x, hi_x);
    derived_here("ThreeStage", "X3", x3, "Tput");
    const bool bracketed3 = x3 >= lo_x - tol && x3 <= hi_x + tol;

    note(conserved && bracketed && monotone && bracketed3
             ? "PASS: ENV-over-LQN meanfield ran; throughput bracketed, monotone in P(UP), "
               "3-stage bracketed."
             : "FAIL: ENV-over-LQN violated one of conservation, bracketing or monotonicity.");
}

// ---------------------------------------------------------------------------
// renv_genqn
// ---------------------------------------------------------------------------

/**
 * The generator model of the repairmen family, solved on its own.
 *
 * The reference is a bare factory in MATLAB; its Python twin's `__main__` is the
 * example body, and it solves the (1.0, 0.5, N=4) instance with MVA and prints
 * one table. `attribute` rather than `section`, because the reference prints no
 * solver banner -- Python labels the record from the solver class instead.
 */
void renv_genqn() {
    Net m = renv_genqn_model(1.0, 0.5, 4.0);
    const Sn& sn = m.get_struct();
    mva::MvaOptions mo;
    Matrix<double> init;
    attribute("MVA");
    print_avg(sn, mva::solver_mva_run_analyzer(sn, mo, init));
}

// ---------------------------------------------------------------------------
// renv_map_fallback / example_mapqn2renv
// ---------------------------------------------------------------------------

namespace {

/** `MMPP2(lambda0, lambda1, sigma0, sigma1)` as the (D0, D1) the reader builds. */
D mmpp2(double l0, double l1, double s0, double s1) {
    Matrix<double> D0(2, 2), D1(2, 2);
    D1(0, 0) = l0;
    D1(1, 1) = l1;
    D0(0, 0) = -(l0 + s0);
    D0(0, 1) = s0;
    D0(1, 0) = s1;
    D0(1, 1) = -(l1 + s1);
    return D::map_dist(D0, D1, lang::ProcessType::MMPP2);
}

/** Think -> FCFS queue with MMPP(2) service, one closed class of five. */
Net mmpp_closed(const std::string& name) {
    Net m(name);
    Delay think(m, "Think");
    Queue q(m, "Q1", SchedStrategy::FCFS);
    ClosedClass c(m, "C1", 5.0, think);
    think.set_service(c, Exp(1.0));
    q.set_service(c, mmpp2(1.0, 10.0, 0.2, 0.3));
    Routing P;
    cyclic(P, c, {think, q});
    m.link(P);
    return m;
}

}  // namespace

/**
 * The random-environment image of a closed network with MMPP(2) service.
 *
 * `map2renv` turns each modulating phase into a stage whose network runs at that
 * phase's completion rate, with the modulating generator's off-diagonal as the
 * stage arcs. The image is exact: the environment and the MMPP describe the same
 * process, which is what makes it a legitimate fallback for a solver that cannot
 * take the MMPP directly.
 */
void renv_map_fallback() {
    Net m = mmpp_closed("mmppClosed");
    io::Map2RenvInfo<double> info;
    Env e = io::map2renv(m.get_struct(), &info);

    std::printf("Random-environment image: %zu stages, phase orders %s, MMPP image: %d\n",
                info.nstages, info.orders.empty() ? "[]" : "[2]",
                info.is_mmpp ? 1 : 0);
    stage_table(e);

    section("CTMC");
    ctmc::CtmcOptions copt;
    print_avg(m.get_struct(), ctmc::solver_ctmc_run_analyzer(m.get_struct(), copt));

    // The two environment LIMITS the reference reaches through
    // `options.config.map_env_method`. Both are available here, but note what
    // differs: MATLAB's mapEnvApprox solves each stage with the SAME solver
    // class it was called on (MVA here), while SolverEnvLimit in this port
    // solves stages with the fluid engine and refuses any other stage solver.
    // For a five-job exponential closed stage the fluid steady state is not the
    // exact MVA one, so these numbers are the port's, not the reference's.
    for (const char* limit : {"dec", "avg"}) {
        env::EnvOptions o;
        o.method = limit;
        o.iter_max = 100;
        o.iter_tol = 0.01;
        section(std::string("ENV(") + limit + ", fluid stages)", "ENV");
        const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
        env_banner(e, o, r);
        print_env_avg(e.stage(0).model, r);
    }

    // TODO(cpp): AvgTableMVA = MVA(model).getAvgTable()
    // TODO(cpp): options.config.map_env = 'off' -> MVA refuses the MMPP service
    na("SolverMVA map_env fallback",
       "the reference's MVA INTERCEPTS an unsupported MAP/MMPP service and re-solves it "
       "through the random-environment image automatically; this port has no map_env knob "
       "and mva_feature_set('default') declares neither MAP nor MMPP2, so solver_mva_run_analyzer "
       "refuses the model outright rather than falling back. The image itself is above, and "
       "for this model the reference's own 'auto' selection resolves to the 'avg' limit");
}

/**
 * The same conversion reached through `MAPQN2RENV`, the reference's named entry
 * point, with the stage networks then solved one at a time.
 *
 * `mapqn2renv` is a pure delegate to `map2renv` in both codebases; it is kept
 * because the reference names it, and because the per-stage solve below is the
 * part that shows what the image is FOR -- each stage is an ordinary closed
 * network that any solver can take.
 */
void example_mapqn2renv() {
    note("=== MMPP2 Closed QN to Random Environment Transformation ===");
    Net m = mmpp_closed("MMPP_ClosedQN");
    kv("Population", 5.0);
    note("Topology: Delay -> MMPP2 Queue -> Delay (cyclic), think time Exp(1.0)");

    io::Map2RenvInfo<double> info;
    Env e = io::mapqn2renv(m.get_struct(), &info);
    stage_table(e);

    env::EnvOptions o;
    o.iter_max = 100;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;

    section("ENV(FLD)", "ENV");
    const env::EnvAnalyzerSolution<double> rf = env::solver_env(e, o);
    env_banner(e, o, rf);
    print_env_avg(e.stage(0).model, rf);

    section("ENV(CTMC)", "ENV");
    env::EnvOptions oc = o;
    oc.stage_solver = "ctmc";
    oc.stage_cutoff = 100.0;
    const env::EnvAnalyzerSolution<double> rc = env::solver_env(e, oc);
    env_banner(e, oc, rc);
    print_env_avg(e.stage(0).model, rc);

    // Each stage is an ordinary closed network once the modulation is carried by
    // the environment rather than by the service process, which is the point of
    // the transformation: MVA cannot take the MMPP, but it can take these.
    section("MVA (per stage)");
    mva::MvaOptions mo;
    Matrix<double> init;
    for (std::size_t s = 0; s < e.nstages(); ++s) {
        note("  stage " + e.stage(s).name);
        print_avg(e.stage(s).model, mva::solver_mva_run_analyzer(e.stage(s).model, mo, init));
    }

    // TODO(cpp): solver_jmt = JMT(model); print(solver_jmt.getAvgTable())
    na("JMT", "the reference opens with the ORIGINAL MMPP2 model solved by JMT as the "
              "reference point the two environment solves are read against; this port's "
              "example corpus has no JMT wrapper for a closed MMPP service model");
}

LINE_EXAMPLE("advanced/randomEnv", renv_basic);
LINE_EXAMPLE("advanced/randomEnv", renv_genqn);
LINE_EXAMPLE("advanced/randomEnv", renv_map_fallback);
LINE_EXAMPLE("advanced/randomEnv", example_mapqn2renv);
LINE_EXAMPLE("advanced/randomEnv", renv_twostages_repairmen);
LINE_EXAMPLE("advanced/randomEnv", renv_threestages_repairmen);
LINE_EXAMPLE("advanced/randomEnv", renv_fourstages_repairmen);
LINE_EXAMPLE("advanced/randomEnv", renv_node_breakdown);
LINE_EXAMPLE("advanced/randomEnv", renv_lqn_twostages);

}  // namespace examples
}  // namespace line
