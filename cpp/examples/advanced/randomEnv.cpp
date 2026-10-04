/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `advanced/randomEnv`: a queueing network modulated by a random environment.
 *
 * Each stage of the environment holds a COMPLETE network, and what couples the
 * stages is that the jobs present at a switch are carried into the next stage.
 * `env::solver_env` is the analyzer selection of `SolverENV.init`; most examples
 * here drive its MEAN-FIELD arm, which collapses the carried state to marginal
 * mean queue lengths, and `renv_container_terminal` drives the STATE-VECTOR arm,
 * which carries the whole joint law across a switch. Either arm takes fluid or
 * CTMC stages (`EnvOptions::stage_solver`); only a CACHE stage is fluid-only,
 * and it is refused by name under any other stage solver.
 *
 * TWO THINGS THE REFERENCE ASKS FOR THAT ARE REFUSED BY NAME rather than
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
 *
 * ONE UTILIZATION CONVENTION worth stating, since `renv_container_terminal`
 * reads a conditional per-stage law: `solver_ctmc_avg_from_pi` reports `max` of
 * the arrival- and departure-based utilization, exactly as
 * `solver_ctmc_avg_from_pi.m` does, and the two estimates coincide only in
 * equilibrium. A stage conditional law is not one, so a Util column read that
 * way sits slightly above `T*E[S]/c` while QLen and Tput stay exact.
 */

#include <algorithm>
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
#include "line/solvers/env/env_generator.h"
#include "line/io/map2renv.h"
#include "line/solvers/ctmc/ctmc_stationary.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/map_env.h"
#include "line/solvers/map_env_stages.h"
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
 * CTMC stage solvers; `env::env_get_generator` (`env_generator.h`) assembles it
 * block by block as the reference does, and the ENV solve follows.
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

    // `ENV(envModel, @(model) CTMC(model, 'exact', soptions)).getGenerator()`.
    {
        const env::EnvGenerator<double> g = env::env_get_generator(e, ctmc::CtmcOptions());
        section("ENV.getGenerator", "ENV");
        print_matrix("infGen", g.Q);
        for (std::size_t s = 0; s < E; ++s)
            print_matrix("stageInfGen{" + std::to_string(s + 1) + "}", g.stage_Q[s]);
    }
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
    // `options.config.map_env_method`, here with the port's DEFAULT stage
    // backend, which is still the fluid one. For a five-job exponential closed
    // stage the fluid steady state is not the exact MVA one, so these two
    // tables and the `map_env` run below differ by the stage solver alone.
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

    // `MVA(model).getAvgTable()` on the model AS WRITTEN: the MMPP2 service is
    // not in `mva_feature_set('default')`, so the solve is intercepted and the
    // random-environment image above is solved instead, EVERY STAGE WITH MVA.
    // That is `mapEnvApprox` and it is what makes the numbers the reference's
    // rather than the fluid backend's.
    section("MVA through the map_env fallback", "MVA");
    mva::MvaOptions mopt;
    Matrix<double> minit;
    const solvers::MapEnvConfig mecfg;  // mode="default", method="auto"
    const mva::AvgResult<double> rme = solvers::run_avg<double>(
        m.get_struct(), "SolverMVA",
        qn::mva_feature_set(mva::resolve_method(m.get_struct(), mopt.method), m.get_struct()),
        mecfg,
        [&mopt, &minit](const qn::NetworkStruct<double>& stage) {
            return mva::solver_mva_run_analyzer(stage, mopt, minit);
        },
        solvers::mva_stage_fn(mopt), mopt.method);
    std::printf("actualmethod=%s\n", rme.actualmethod.c_str());
    print_avg(m.get_struct(), rme);

    // `options.config.map_env = 'off'` restores the plain rejection, which is
    // the point of the knob: a caller who wants to KNOW the model is out of
    // reach must be able to switch the substitution off.
    solvers::MapEnvConfig offcfg;
    offcfg.mode = "off";
    try {
        Matrix<double> oinit;
        solvers::run_avg<double>(
            m.get_struct(), "SolverMVA",
            qn::mva_feature_set(mva::resolve_method(m.get_struct(), mopt.method), m.get_struct()),
            offcfg,
            [&mopt, &oinit](const qn::NetworkStruct<double>& stage) {
                return mva::solver_mva_run_analyzer(stage, mopt, oinit);
            },
            solvers::mva_stage_fn(mopt), mopt.method);
        note("map_env='off': the model was accepted, which it should not have been");
    } catch (const std::exception& ex) {
        std::printf("map_env='off' -> %s\n", ex.what());
    }
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

    // The reference point the environment solves are read against: the ORIGINAL
    // MMPP2 model, simulated.
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m));
}

// ---------------------------------------------------------------------------
// renv_container_terminal
// ---------------------------------------------------------------------------

namespace {

/** One hour-stage network: N carriers cycling yard-Delay -> quay-crane Queue. */
Net terminal_model(double stage_rate, double crane_rate, double n_cranes, double N) {
    Net m("Terminal");
    Delay yard(m, "Yard");
    Queue cranes(m, "QuayCranes", SchedStrategy::FCFS);
    cranes.set_number_of_servers(n_cranes);
    ClosedClass containers(m, "Containers", N, yard);
    yard.set_service(containers, Exp(stage_rate));
    cranes.set_service(containers, Exp(crane_rate));
    Routing P;
    cyclic(P, containers, {yard, cranes});
    m.link(P);
    return m;
}

/**
 * The same hour-stage with the terminal internals collapsed into one closed
 * load-dependent FES serving at mu(n) when n containers are inside it.
 */
Net terminal_fes_model(double stage_rate, const std::vector<double>& fes_rate) {
    const double N = double(fes_rate.size());
    Net m("TerminalFES");
    Delay yard(m, "Yard");
    Queue fesq(m, "TerminalFES", SchedStrategy::PS);
    ClosedClass containers(m, "Containers", N, yard);
    yard.set_service(containers, Exp(stage_rate));
    fesq.set_service(containers, Exp(fes_rate[0]));  // the base rate mu(1)
    std::vector<double> alpha(fes_rate.size());
    for (std::size_t i = 0; i < fes_rate.size(); ++i) alpha[i] = fes_rate[i] / fes_rate[0];
    fesq.set_load_dependence(alpha);  // alpha(n) = mu(n)/mu(1)
    Routing P;
    cyclic(P, containers, {yard, fesq});
    m.link(P);
    return m;
}

/**
 * Norton flow-equivalent rates: the throughput of the isolated internal
 * subnetwork (quay cranes -> stacking cranes, both processor sharing) at
 * n = 1..N containers, with the rest of the terminal short-circuited.
 */
std::vector<double> fes_rate_curve(double mu_quay, double mu_stack, std::size_t N) {
    std::vector<double> mu(N, 0.0);
    for (std::size_t n = 1; n <= N; ++n) {
        Net sub("TerminalInternals");
        Delay ref(sub, "ShortCircuit");
        Queue quay(sub, "QuayCranes", SchedStrategy::PS);
        Queue stack(sub, "StackingCranes", SchedStrategy::PS);
        ClosedClass cls(sub, "Containers", double(n), ref);
        ref.set_service(cls, Exp(1e6));  // the short circuit, ~instantaneous
        quay.set_service(cls, Exp(mu_quay));
        stack.set_service(cls, Exp(mu_stack));
        Routing P;
        cyclic(P, cls, {ref, quay, stack});
        sub.link(P);
        mva::MvaOptions o;
        o.method = "exact";
        Matrix<double> init;
        mu[n - 1] = mva::solver_mva_run_analyzer(sub.get_struct(), o, init).TN(0, 0);
    }
    return mu;
}

/** The three metrics the exact joint chain is read for. */
struct JointAvg {
    Matrix<double> QN, UN, TN;
};

/**
 * Day-averaged metrics from the exact joint (hour x network-state) CTMC.
 *
 * THE JOINT GENERATOR IS BUILT HERE because this port exposes no `getGenerator`
 * on the ENV entry: the reference asks its SolverENV for the flattened
 * generator and solves that, and the flattening is the loop below -- each
 * stage's own generator on the diagonal block, the environment's rate on the
 * IDENTITY of the off-diagonal one, since a stage switch carries the network
 * state across unchanged. That identity is also why every stage must enumerate
 * the same space, which holds here (the 24 stages differ in their rates alone)
 * and is asserted rather than assumed.
 *
 * The per-stage metrics come from `solver_ctmc_avg_from_pi` on the CONDITIONAL
 * law of each block, which is the reference's own reduction. One consequence is
 * worth naming: that function reports `max` of the arrival- and
 * departure-based utilization estimates, and the two coincide only in
 * equilibrium. A stage conditional law is not one, so the Util column sits
 * slightly above `T*E[S]/c` here while QLen and Tput are exact.
 */
JointAvg exact_joint_metrics(std::vector<Net>& stage, const std::vector<double>& env_rate) {
    const std::size_t E = stage.size();
    std::vector<ctmc::CtmcSolution<double> > sol;
    ctmc::CtmcOptions co;
    co.method = "exact";
    for (std::size_t e = 0; e < E; ++e)
        sol.push_back(ctmc::solver_ctmc_analyzer(stage[e].get_struct(), co));

    const std::size_t ns = sol[0].chain.Q.rows();
    for (std::size_t e = 1; e < E; ++e)
        if (sol[e].chain.Q.rows() != ns)
            throw InputError(
                "renv_container_terminal: the stages do not enumerate one state space, so the "
                "environment arcs are not the identity on it");

    Matrix<double> Q(E * ns, E * ns, 0.0);
    for (std::size_t e = 0; e < E; ++e) {
        for (std::size_t i = 0; i < ns; ++i)
            for (std::size_t j = 0; j < ns; ++j)
                if (i != j) Q(e * ns + i, e * ns + j) = sol[e].chain.Q(i, j);
        const std::size_t h = (e + 1) % E;  // the cyclic daily arc
        for (std::size_t i = 0; i < ns; ++i) Q(e * ns + i, h * ns + i) += env_rate[e];
    }
    for (std::size_t i = 0; i < E * ns; ++i) {  // re-close the rows, as ctmc_makeinfgen does
        Q(i, i) = 0.0;
        double row = 0.0;
        for (std::size_t j = 0; j < E * ns; ++j) row += Q(i, j);
        Q(i, i) = -row;
    }

    const std::vector<double> pi = ctmc::ctmc_stationary(Q).pi;
    const Sn& sn0 = stage[0].get_struct();
    JointAvg out;
    out.QN = Matrix<double>(sn0.nstations, sn0.nclasses, 0.0);
    out.UN = Matrix<double>(sn0.nstations, sn0.nclasses, 0.0);
    out.TN = Matrix<double>(sn0.nstations, sn0.nclasses, 0.0);
    for (std::size_t e = 0; e < E; ++e) {
        std::vector<double> blk(pi.begin() + e * ns, pi.begin() + (e + 1) * ns);
        double prob = 0.0;
        for (std::size_t i = 0; i < ns; ++i) prob += blk[i];
        if (prob > 0.0)
            for (std::size_t i = 0; i < ns; ++i) blk[i] /= prob;
        const ctmc::CtmcAvg<double> a =
            ctmc::solver_ctmc_avg_from_pi(stage[e].get_struct(), sol[e].chain, blk);
        for (std::size_t i = 0; i < out.QN.rows(); ++i)
            for (std::size_t c = 0; c < out.QN.cols(); ++c) {
                out.QN(i, c) += prob * a.QN(i, c);
                out.UN(i, c) += prob * a.UN(i, c);
                out.TN(i, c) += prob * a.TN(i, c);
            }
    }
    return out;
}

/** One station's metric, totalled over the classes it carries. */
double station_sum(const Matrix<double>& A, std::size_t st) {
    double s = 0.0;
    for (std::size_t c = 0; c < A.cols(); ++c) s += A(st, c);
    return s;
}

/** One row of the analyzer comparison the reference prints. */
void compare_row(const char* label, const Matrix<double>& Q, const Matrix<double>& U,
                 const Matrix<double>& T, std::size_t st) {
    std::printf("%-10s %12.5f %12.5f %12.5f\n", label, station_sum(Q, st), station_sum(U, st),
                station_sum(T, st));
}

/** The two error lines that read the couplings against the exact chain. */
void compare_errors(const Matrix<double>& Qsv, const Matrix<double>& Qmf, const Matrix<double>& Qex,
                    const Matrix<double>& Usv, const Matrix<double>& Umf, const Matrix<double>& Uex,
                    std::size_t st) {
    std::printf("\nQLen error vs exact:  statevec = %.3e , meanfield = %.3e\n",
                std::fabs(station_sum(Qsv, st) - station_sum(Qex, st)),
                std::fabs(station_sum(Qmf, st) - station_sum(Qex, st)));
    std::printf("Util error vs exact:  statevec = %.3e , meanfield = %.3e\n",
                std::fabs(station_sum(Usv, st) - station_sum(Uex, st)),
                std::fabs(station_sum(Umf, st) - station_sum(Uex, st)));
}

/** The 24 hourly stages wired into the fixed daily cycle 1 -> 2 -> ... -> 24 -> 1. */
Env daily_cycle(const std::string& nm, std::vector<Net>& stage,
                const std::vector<double>& env_rate) {
    const std::size_t E = stage.size();
    Env e(nm, E);
    for (std::size_t h = 0; h < E; ++h) {
        char hour[16];
        std::snprintf(hour, sizeof(hour), "Hour%02zu", h);
        e.set_stage(h, hour, "operational", stage[h].get_struct());
    }
    for (std::size_t h = 0; h < E; ++h) e.add_transition(h, (h + 1) % E, Exp(env_rate[h]));
    return e;
}

}  // namespace

/**
 * A daily-cycle container terminal, solved by the STATE-VECTOR coupling and
 * read against the exact joint chain.
 *
 * The Rotterdam terminal of Dhingra et al.: container handling demand varies
 * over the 24 hours of a day, each hour is one environment stage with its own
 * demand intensity and random (exponential) duration, and the environment
 * visits the stages in the fixed daily cycle. The semi-open SOQN of the
 * original study is rendered as the equivalent finite-token CLOSED network that
 * an environment stage must be: N straddle carriers cycle between a yard
 * staging Delay and a multi-server quay-crane Queue, and the hourly demand
 * modulates the yard staging rate, so the cranes congest during the peaks.
 *
 * WHAT THE THREE ROWS MEASURE. `statevec` carries the whole joint distribution
 * across a stage switch, `meanfield` collapses it to marginal mean queue
 * lengths, and `exact` is the stationary law of the full (hour x network-state)
 * chain built above. The state-vector blend reproduces the exact answer to nine
 * digits here, which is the property the example exists to show; the mean-field
 * collapse does not.
 *
 * The second section repeats it with the internal handling (quay cranes ->
 * stacking cranes) collapsed by Norton's theorem into one closed
 * load-dependent FES, so `sn.lldscaling` reaches the state-vector analyzer's
 * CTMC backend. The third analyses the OPEN counterpart, where the hourly
 * demand is an MMPP(24) into a multi-server queue with an unbounded buffer:
 * that is outside the closed-CTMC state-vector path and is what MAM solves
 * exactly as a QBD.
 */
void renv_container_terminal() {
    const double daily_rate_hr[24] = {6,   30,  40,  62,  76,  79,  119, 164, 152, 130, 79, 70,
                                      57,  57,  113, 130, 162, 202, 148, 118, 92,  62,  36, 8};
    const double daily_duration_hr[24] = {0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62,
                                          0.36, 1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57,
                                          0.34, 0.22, 0.47, 1.33, 1.56, 0.93, 0.71, 1.11};
    const std::size_t E = 24;        // hourly environment stages
    const double N = 6.0;            // straddle carriers circulating (the closed tokens)
    const double n_cranes = 2.0;     // quay cranes, a multi-server FCFS station
    const double crane_rate = 16.0;  // moves per hour served by one crane

    std::vector<Net> stage;
    std::vector<double> env_rate(E);
    for (std::size_t h = 0; h < E; ++h) {
        // The per-token yard completion rate of hour h: the closed-network
        // image of that hour's external demand intensity.
        stage.push_back(terminal_model(daily_rate_hr[h] / N, crane_rate, n_cranes, N));
        env_rate[h] = 1.0 / daily_duration_hr[h];
    }
    Env e = daily_cycle("RotterdamDailyCycle", stage, env_rate);

    std::printf("Rotterdam container terminal: %zu hourly stages, %g carriers, %g cranes.\n", E, N,
                n_cranes);

    // The reference's stage solver is `CTMC(m, 'exact', 'timespan', [0,T])` with
    // T = 12 hours, which covers the hour-duration CDF past 99.9% for every
    // stage, so the quadrature below is not truncating the blend.
    env::EnvOptions base;
    base.iter_max = 100;
    base.iter_tol = 1e-5;
    base.stage_solver = "ctmc";
    base.timespan_end = 12.0;

    env::EnvOptions osv = base;
    osv.method = "statevec";
    const env::EnvAnalyzerSolution<double> sv = env::solver_env(e, osv);

    env::EnvOptions omf = base;
    omf.method = "meanfield";
    const env::EnvAnalyzerSolution<double> mf = env::solver_env(e, omf);

    const JointAvg ex = exact_joint_metrics(stage, env_rate);

    const std::size_t crane = 1;  // the QuayCranes station
    std::printf("\n=== Day-averaged quay-crane metrics (container terminal) ===\n");
    std::printf("%-10s %12s %12s %12s\n", "analyzer", "QLen", "Util", "Tput");
    compare_row("exact", ex.QN, ex.UN, ex.TN, crane);
    compare_row("statevec", sv.QN, sv.UN, sv.TN, crane);
    compare_row("meanfield", mf.QN, mf.UN, mf.TN, crane);
    compare_errors(sv.QN, mf.QN, ex.QN, sv.UN, mf.UN, ex.UN, crane);

    std::printf("\nDay-averaged crane-queue table (state-vector analyzer):\n");
    section("ENV(statevec, CTMC stages)", "ENV");
    env_banner(e, osv, sv);
    print_env_avg(e.stage(0).model, sv);

    std::printf("\n=== Closed FES variant (Norton-aggregated terminal internals) ===\n");
    const std::vector<double> fes_rate = fes_rate_curve(16.0, 20.0, std::size_t(N));
    std::printf("Norton FES rate curve mu(n) = [");
    for (std::size_t i = 0; i < fes_rate.size(); ++i)
        std::printf("%s%.5f", i ? " " : "", fes_rate[i]);
    std::printf("]\n");

    std::vector<Net> stage_fes;
    for (std::size_t h = 0; h < E; ++h)
        stage_fes.push_back(terminal_fes_model(daily_rate_hr[h] / N, fes_rate));
    Env efes = daily_cycle("RotterdamDailyCycleFES", stage_fes, env_rate);

    const env::EnvAnalyzerSolution<double> svf = env::solver_env(efes, osv);
    const env::EnvAnalyzerSolution<double> mff = env::solver_env(efes, omf);
    const JointAvg exf = exact_joint_metrics(stage_fes, env_rate);

    const std::size_t fes = 1;  // the load-dependent FES station
    std::printf("%-10s %12s %12s %12s\n", "analyzer", "QLen", "Util", "Tput");
    compare_row("exact", exf.QN, exf.UN, exf.TN, fes);
    compare_row("statevec", svf.QN, svf.UN, svf.TN, fes);
    compare_row("meanfield", mff.QN, mff.UN, mff.TN, fes);
    compare_errors(svf.QN, mff.QN, exf.QN, svf.UN, mff.UN, exf.UN, fes);

    // The open counterpart: the 24-hour demand as a Markov-modulated Poisson
    // stream rather than a finite carrier pool. The ring of hours IS the
    // modulating generator and the hourly demand its per-phase arrival rate, so
    // the model is an MMPP(24)/M/c with an unbounded buffer -- outside the
    // closed-CTMC state-vector path, and exactly the QBD that MAM takes.
    std::printf("\n=== Open system via MAM (MMPP(24)/M/c, infinite buffer) ===\n");
    Matrix<double> D0(E, E, 0.0), D1(E, E, 0.0);
    for (std::size_t h = 0; h < E; ++h) {
        D0(h, h) = -env_rate[h] - daily_rate_hr[h];
        D0(h, (h + 1) % E) = env_rate[h];
        D1(h, h) = daily_rate_hr[h];
    }
    const double c_open = 8.0;  // an open stream is not throttled by a token pool
    Net open_model("OpenTerminal");
    Source ships(open_model, "Ships");
    Queue quay_open(open_model, "QuayCranes", SchedStrategy::FCFS);
    Sink done(open_model, "Departures");
    OpenClass oc(open_model, "Containers");
    ships.set_arrival(oc, MAP(D0, D1));
    quay_open.set_service(oc, Exp(crane_rate));
    quay_open.set_number_of_servers(c_open);
    Routing Po;
    serial(Po, oc, {ships, quay_open, done});
    open_model.link(Po);

    double num = 0.0, den = 0.0;
    for (std::size_t h = 0; h < E; ++h) {
        num += daily_duration_hr[h] * daily_rate_hr[h];
        den += daily_duration_hr[h];
    }
    const double mean_lambda = num / den;
    std::printf("MMPP mean arrival = %.2f/hr, %g cranes @ %.0f/hr, rho = %.3f\n", mean_lambda,
                c_open, crane_rate, mean_lambda / (c_open * crane_rate));
    section("MAM");
    print_avg(open_model.get_struct(),
              mam::solver_mam_run_analyzer(open_model.get_struct(), mam::MamOptions()));
}

// ---------------------------------------------------------------------------
// renv_rotterdam_blending
// ---------------------------------------------------------------------------

namespace {

/** `PoissonSOQNSpec`: the geometry of the Rotterdam semi-open network. */
struct SoqnParams {
    std::size_t num_entry_servers, num_stacks, num_exit_servers;
    double entry_service, travel_to_stack, stack_service, travel_to_exit, exit_service;
    SoqnParams()
        : num_entry_servers(6),
          num_stacks(29),
          num_exit_servers(6),
          entry_service(6.0),
          travel_to_stack(5.6),
          stack_service(6.0),
          travel_to_exit(5.6),
          exit_service(6.0) {}
};

/**
 * The flow-equivalent throughput curve of one SOQN subnetwork at l = 1..N
 * tokens, by exact MVA on the underlying closed network.
 *
 * `upstream` is S1, entry gates -> travel -> 29 stacks -> back; the downstream
 * S2 is travel-to-exit <-> exit gates. `mu[0]` is left at zero so the curve is
 * indexed by occupancy directly, which is how the QBD blocks read it.
 */
std::vector<double> fes_curve(std::size_t N, const SoqnParams& p, bool upstream) {
    std::vector<double> mu(N + 1, 0.0);
    for (std::size_t l = 1; l <= N; ++l) {
        Net m(upstream ? "S1" : "S2");
        if (upstream) {
            Queue eg(m, "EntryGates", SchedStrategy::FCFS);
            eg.set_number_of_servers(double(p.num_entry_servers));
            Delay tv(m, "TravelToStack");
            std::vector<std::size_t> st;
            for (std::size_t i = 0; i < p.num_stacks; ++i) {
                Queue q(m, "Stack" + std::to_string(i + 1), SchedStrategy::FCFS);
                q.set_number_of_servers(1.0);
                st.push_back(q);
            }
            ClosedClass cls(m, "Trucks", double(l), eg);
            eg.set_service(cls, Exp(1.0 / p.entry_service));
            tv.set_service(cls, Exp(1.0 / p.travel_to_stack));
            for (std::size_t i = 0; i < st.size(); ++i)
                m.set_service(st[i], cls, Exp(1.0 / p.stack_service));
            Routing P;
            P.set(cls, cls, eg, tv, 1.0);
            for (std::size_t i = 0; i < st.size(); ++i) {
                P.set(cls, cls, tv, st[i], 1.0 / double(p.num_stacks));
                P.set(cls, cls, st[i], eg, 1.0);
            }
            m.link(P);
        } else {
            Delay tv(m, "TravelToExit");
            Queue xg(m, "ExitGates", SchedStrategy::FCFS);
            xg.set_number_of_servers(double(p.num_exit_servers));
            ClosedClass cls(m, "Trucks", double(l), tv);
            tv.set_service(cls, Exp(1.0 / p.travel_to_exit));
            xg.set_service(cls, Exp(1.0 / p.exit_service));
            Routing P;
            P.set(cls, cls, tv, xg, 1.0);
            P.set(cls, cls, xg, tv, 1.0);
            m.link(P);
        }
        mva::MvaOptions o;
        o.method = "exact";
        Matrix<double> init;
        // The gate station's throughput IS the subnetwork's completion rate at
        // this occupancy, since every token passes it exactly once per cycle.
        mu[l] = mva::solver_mva_run_analyzer(m.get_struct(), o, init).TN(upstream ? 0 : 1, 0);
    }
    return mu;
}

/**
 * A banded LU factorization with NO PIVOTING, in column-major band storage.
 *
 * The unpivoted form is what makes the blend affordable: every one of the 24
 * resolvents is factored ONCE and then applied on every sweep of the fixed
 * point, so the iteration costs band solves rather than 24 fresh solves per
 * sweep. Skipping the pivot search is safe here and not a shortcut: the matrix
 * is `(sI - Q)^T` with `s > 0` and `Q` a generator, so `sI - Q` is strictly
 * row diagonally dominant and its transpose is strictly column diagonally
 * dominant, which is exactly the condition under which Gaussian elimination
 * without pivoting is stable.
 */
struct Banded {
    std::size_t n, kl, ku, ld;
    std::vector<double> a;
    Banded(std::size_t n_, std::size_t kl_, std::size_t ku_)
        : n(n_), kl(kl_), ku(ku_), ld(kl_ + ku_ + 1), a((kl_ + ku_ + 1) * n_, 0.0) {}
    double& at(std::size_t i, std::size_t j) { return a[j * ld + (i + ku - j)]; }
    double at(std::size_t i, std::size_t j) const { return a[j * ld + (i + ku - j)]; }
    void factor() {
        for (std::size_t j = 0; j + 1 < n; ++j) {
            const double d = at(j, j);
            const std::size_t imax = std::min(n - 1, j + kl);
            const std::size_t kmax = std::min(n - 1, j + ku);
            for (std::size_t i = j + 1; i <= imax; ++i) {
                const double l = at(i, j) / d;
                at(i, j) = l;
                if (l == 0.0) continue;
                for (std::size_t k = j + 1; k <= kmax; ++k) at(i, k) -= l * at(j, k);
            }
        }
    }
    void solve(std::vector<double>& b) const {
        for (std::size_t j = 0; j + 1 < n; ++j) {
            const std::size_t imax = std::min(n - 1, j + kl);
            for (std::size_t i = j + 1; i <= imax; ++i) b[i] -= at(i, j) * b[j];
        }
        for (std::size_t jj = n; jj-- > 0;) {
            const std::size_t kmax = std::min(n - 1, jj + ku);
            for (std::size_t k = jj + 1; k <= kmax; ++k) b[jj] -= at(jj, k) * b[k];
            b[jj] /= at(jj, jj);
        }
    }
};

/**
 * `(sI - Q)^T` of the level-dependent QBD, banded and ready to factor.
 *
 * The state is `(n, k)`: level `n` counts the jobs upstream of S2 (those inside
 * S1 plus the external backlog) and phase `k` the S2 occupancy, so an S1
 * completion moves `(n,k) -> (n-1,k+1)`, an S2 completion `(n,k) -> (n,k-1)`
 * and an arrival `(n,k) -> (n+1,k)`. S1 serves at `mu1(min(n, N-k))`, since a
 * token can only be inside S1 if a token of the pool is free to hold it, and
 * the top level `Mtr` is the truncation, where arrivals are dropped.
 *
 * Everything is written TRANSPOSED, because the resolvent is applied to a ROW
 * vector `pi` and the solve wants a column system.
 */
Banded soqn_resolvent_matrix(std::size_t N, double lam, const std::vector<double>& mu1,
                             const std::vector<double>& mu2, std::size_t tail_factor, double s) {
    const std::size_t M = N + 1, Mtr = N + tail_factor * N, dim = (Mtr + 1) * M;
    Banded B(dim, M, M - 1);
    for (std::size_t n = 0; n <= Mtr; ++n) {
        const double lam_eff = (n == Mtr) ? 0.0 : lam;
        for (std::size_t k = 0; k <= N; ++k) {
            const std::size_t row = n * M + k;
            const double m1 = mu1[std::min(n, N - k)], m2 = mu2[k];
            B.at(row, row) = s + lam_eff + m1 + m2;
            if (k >= 1) B.at(n * M + (k - 1), row) = -m2;
            if (n < Mtr) B.at((n + 1) * M + k, row) = -lam;
            if (n >= 1 && k + 1 <= N) B.at((n - 1) * M + (k + 1), row) = -m1;
        }
    }
    return B;
}

/** The two Table-13 measures the blend is read for. */
struct BlendResult {
    double W, Qex;
};

/**
 * The cyclic resolvent blend over the 24 hourly environments.
 *
 * One visit to stage h applies `s*pi*(sI-Q_h)^{-1}` with `s = 1/E[duration]`,
 * which is the LAW OF THE STATE AT A MEMORYLESS EXIT: the resolvent is the
 * time-average over an exponential horizon, so no quadrature grid is needed and
 * the exit vector is exact rather than sampled. Entry vectors are chained
 * around the cycle to an L1 fixed point, and the day average weights each stage
 * by its time fraction.
 */
BlendResult blend_soqn(std::size_t N, const std::vector<double>& lam,
                       const std::vector<double>& dur_min, const std::vector<double>& frac_w,
                       const std::vector<double>& mu1, const std::vector<double>& mu2,
                       std::size_t tail_factor) {
    const std::size_t K = lam.size(), M = N + 1, Mtr = N + tail_factor * N, dim = (Mtr + 1) * M;
    std::vector<Banded> A;  // factored ONCE; the sweeps below only back-solve
    for (std::size_t h = 0; h < K; ++h) {
        Banded B = soqn_resolvent_matrix(N, lam[h], mu1, mu2, tail_factor, 1.0 / dur_min[h]);
        B.factor();
        A.push_back(B);
    }
    std::vector<std::vector<double> > pi_enter(K, std::vector<double>(dim, 0.0));
    for (std::size_t h = 0; h < K; ++h) pi_enter[h][0] = 1.0;  // start every stage empty

    for (int it = 0; it < 200; ++it) {
        const std::vector<std::vector<double> > prev = pi_enter;
        for (std::size_t h = 0; h < K; ++h) {
            const double s = 1.0 / dur_min[h];
            std::vector<double> rhs(dim);
            for (std::size_t i = 0; i < dim; ++i) rhs[i] = s * pi_enter[h][i];
            A[h].solve(rhs);
            pi_enter[(h + 1) % K] = rhs;  // the exit law of h is the entry law of h+1
        }
        double l1 = 0.0;
        for (std::size_t h = 0; h < K; ++h) {
            double d = 0.0;
            for (std::size_t i = 0; i < dim; ++i) d += std::fabs(pi_enter[h][i] - prev[h][i]);
            l1 = std::max(l1, d);
        }
        if (l1 < 1e-10) break;
    }

    std::vector<double> pi_avg(dim, 0.0);
    for (std::size_t h = 0; h < K; ++h) {
        const double s = 1.0 / dur_min[h];
        std::vector<double> rhs(dim);
        for (std::size_t i = 0; i < dim; ++i) rhs[i] = s * pi_enter[h][i];
        A[h].solve(rhs);
        for (std::size_t i = 0; i < dim; ++i) pi_avg[i] += frac_w[h] * rhs[i];
    }
    double tot = 0.0;
    for (std::size_t i = 0; i < dim; ++i) {
        // The truncation leaves the blend a hair off a probability vector; the
        // clamp and the renormalization are the reference's own closing step.
        if (pi_avg[i] < 0.0) pi_avg[i] = 0.0;
        tot += pi_avg[i];
    }
    for (std::size_t i = 0; i < dim; ++i) pi_avg[i] /= tot;

    double QLex = 0.0, QL1 = 0.0, QL2 = 0.0, tput = 0.0;
    for (std::size_t n = 0; n <= Mtr; ++n)
        for (std::size_t k = 0; k <= N; ++k) {
            const double p = pi_avg[n * M + k];
            if (p == 0.0) continue;
            // A token pool of N holds n + k jobs at most; the excess waits
            // OUTSIDE the network, which is the external queue the study reports.
            QLex += ((n + k > N) ? double(n + k - N) : 0.0) * p;
            QL1 += double(std::min(n, N - k)) * p;
            QL2 += double(k) * p;
            if (k >= 1) tput += mu2[k] * p;
        }
    BlendResult r;
    r.W = (QLex + QL1 + QL2) / tput;  // Little's law over the whole sojourn
    r.Qex = QLex;
    return r;
}

}  // namespace

/**
 * The 24-environment exponential blending accuracy result of the Rotterdam
 * container-terminal study, against its published simulation.
 *
 * The terminal is a SEMI-OPEN network: trucks arrive as a Poisson stream at the
 * hour's rate, a pool of N tokens admits them, and a truck that finds no token
 * waits outside. Inside, two flow-equivalent servers in tandem stand for the
 * upstream half (entry gates, travel, 29 stacks) and the downstream half
 * (travel to exit, exit gates), each calibrated by exact MVA on the closed
 * subnetwork it replaces. That gives a level-dependent QBD whose blocks are the
 * same ones the ENV state-vector analyzer assembles internally, which is why
 * this example builds them directly rather than through `solver_env`: the
 * quantities it reports (external wait and external queue) live OUTSIDE the
 * token pool and so outside any stage network's own metric table.
 *
 * The comparison is against Table 13 of the study's discrete-event simulation.
 */
void renv_rotterdam_blending() {
    const double rate_hr[24] = {6,   30,  40,  62,  76,  79,  119, 164, 152, 130, 79, 70,
                                57,  57,  113, 130, 162, 202, 148, 118, 92,  62,  36, 8};
    const double dur_hr[24] = {0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62,
                               0.36, 1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57,
                               0.34, 0.22, 0.47, 1.33, 1.56, 0.93, 0.71, 1.11};
    const std::size_t K = 24;
    std::vector<double> lam(K), dur_min(K), frac_w(K);
    double tot_dur = 0.0;
    for (std::size_t h = 0; h < K; ++h) {
        lam[h] = rate_hr[h] * 0.5 / 60.0;  // arrivals per MINUTE
        dur_min[h] = dur_hr[h] * 60.0;
        tot_dur += dur_min[h];
    }
    for (std::size_t h = 0; h < K; ++h) frac_w[h] = dur_min[h] / tot_dur;

    // The study's simulated ground truth, N = 24..34 (Table 13).
    const double sim_W[11] = {154.0661, 118.3936, 99.962,  88.2379, 80.9272, 75.0887,
                              70.745,   67.4326,  64.9719, 62.7913, 61.2534};
    const double sim_Qex[11] = {89.8508, 63.3970, 49.3650, 40.8783, 35.3546, 30.6543,
                                27.4111, 24.7218, 22.3622, 20.3893, 18.9252};
    const std::size_t tail_factor = 15;  // M_trunc = N*(1 + tailFactor)
    const SoqnParams p;

    std::printf("Rotterdam SOQN exponential blending vs Dhingra simulation (tailFactor=%zu)\n",
                tail_factor);
    std::printf("%4s | %10s %10s %7s | %10s %10s %7s\n", "N", "W_blend", "W_sim", "err%", "Qex_bl",
                "Qex_sim", "err%");
    const std::size_t Nlist[3] = {24, 28, 34};
    for (int i = 0; i < 3; ++i) {
        const std::size_t N = Nlist[i];
        const std::vector<double> mu1 = fes_curve(N, p, true);
        const std::vector<double> mu2 = fes_curve(N, p, false);
        const BlendResult r = blend_soqn(N, lam, dur_min, frac_w, mu1, mu2, tail_factor);
        const std::size_t j = N - 24;
        std::printf("%4zu | %10.4f %10.4f %6.2f | %10.4f %10.4f %6.2f\n", N, r.W, sim_W[j],
                    100.0 * std::fabs(r.W - sim_W[j]) / sim_W[j], r.Qex, sim_Qex[j],
                    100.0 * std::fabs(r.Qex - sim_Qex[j]) / sim_Qex[j]);
    }
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
LINE_EXAMPLE("advanced/randomEnv", renv_container_terminal);
LINE_EXAMPLE("advanced/randomEnv", renv_rotterdam_blending);

}  // namespace examples
}  // namespace line
