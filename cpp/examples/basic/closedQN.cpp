/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/closedQN`: the closed-network half of the basic
 * gallery.
 *
 * Every entry is its reference script's model and its `__main__` block, in the
 * same order and with the same solvers and methods. What the port cannot answer
 * it refuses BY NAME: JMT, LDES and QNS have no C++ engine, and neither
 * `getAvgChainTable` nor `State.fromMarginalAndRunning` have a C++
 * counterpart. Two of these scripts
 * (`cqn_mmpp2_service`, `cqn_twoclass_erl`) run ONLY wrapper solvers, so they
 * build their model and refuse both calls rather than answering with a native
 * solver the reference never asked for -- an analytical number carrying a
 * simulator's name is exactly what the parity harness exists to prevent.
 *
 * NONE OF THESE SCRIPTS DRAWS A RANDOM PARAMETER, so nothing here is a
 * materialized sample; the seeds that appear are solver seeds.
 */

#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "example_util.h"
#include "example_node_table.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/qn/state.h"
#include "line/solvers/auto/solver_auto.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/wrappers/qns/solver_qns.h"

namespace line {
namespace examples {

namespace {

/** `MMPP2(lambda0, lambda1, sigma0, sigma1)` as the (D0, D1) the reader builds. */
D mmpp2(double l0, double l1, double s0, double s1) {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D1(0, 0) = l0;
    D1(1, 1) = l1;
    D0(0, 0) = -(l0 + s0);
    D0(0, 1) = s0;
    D0(1, 0) = s1;
    D0(1, 1) = -(l1 + s1);
    return D::map_dist(D0, D1, lang::ProcessType::MMPP2);
}

/** The one MVA call the scripts make over and over, with its own method. */
void run_mva(const Sn& sn, const std::string& method = "default",
             const std::string& multiserver = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    opt.multiserver = multiserver;
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, opt, init));
}

void run_nc(const Sn& sn, const std::string& method = "default") {
    nc::NcSolverOptions opt;
    opt.method = method;
    print_avg(sn, nc::solver_nc_run_analyzer(sn, opt));
}

void run_ctmc(const Sn& sn) { print_avg(sn, ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions())); }

void run_ssa(const Sn& sn, std::size_t samples, unsigned long seed) {
    ssa::SsaOptions opt;
    opt.seed = seed;
    if (samples) opt.samples = samples;
    print_avg_sim(sn, ssa::solver_ssa(sn, opt));
}

void run_fluid(const Sn& sn) {
    print_avg_sim(sn, fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions()));
}

void run_mam(const Sn& sn) { print_avg(sn, mam::solver_mam_run_analyzer(sn, mam::MamOptions())); }

/**
 * `LINE(model)`, i.e. SolverAUTO: the engine `chooseSolver` ranks first on this
 * model, then that engine's own table.
 *
 * The choice is the port's, so the banner names the engine that actually
 * answered -- an AUTO row that did not say which solver produced it would be
 * unattributable, and the ranking differs from MATLAB's wherever a slot has no
 * C++ engine.
 */
void run_auto(const Sn& sn, unsigned long seed) {
    const autosolver::AutoChoice c = autosolver::auto_choose_avg_solver_ex(sn);
    const std::string name = autosolver::auto_solver_name(c.solver);
    const std::string method = c.method.empty() ? "default" : c.method;
    section("AUTO -> " + name);
    if (c.solver == autosolver::AutoSolver::MVA) {
        run_mva(sn, method);
    } else if (c.solver == autosolver::AutoSolver::NC) {
        run_nc(sn, method);
    } else if (c.solver == autosolver::AutoSolver::CTMC) {
        run_ctmc(sn);
    } else if (c.solver == autosolver::AutoSolver::MAM) {
        run_mam(sn);
    } else if (c.solver == autosolver::AutoSolver::FLUID) {
        run_fluid(sn);
    } else if (c.solver == autosolver::AutoSolver::SSA) {
        run_ssa(sn, 0, seed);
    } else {
        na("AUTO -> " + name, "the ranking chose an engine this port does not build");
    }
}

/**
 * One solver block of a script whose reference loop CATCHES a solver failure.
 *
 * `cqn_repairmen_multi`, `cqn_scheduling_dps` and `cqn_threeclass_hyperl` each
 * wrap every solver in try/except and print the error rather than stopping, so a
 * solver that refuses the model is part of what they report. Letting the
 * exception escape would end the example at the first refusal and hide every row
 * below it. `banner` empty means the callee prints its own; `err_max` is the
 * truncation the reference applies to the message.
 */
void guarded(const std::string& banner, const std::function<void()>& fn,
             const std::string& err_prefix, std::size_t err_max = 0,
             const std::string& golden_key = std::string()) {
    if (!banner.empty()) {
        if (golden_key.empty()) section(banner);
        else section(banner, golden_key);
    }
    try {
        fn();
    } catch (const std::exception& e) {
        std::string msg = e.what();
        if (err_max && msg.size() > err_max) msg = msg.substr(0, err_max);
        std::printf("%s%s\n", err_prefix.c_str(), msg.c_str());
    }
}

/** Print the rows `State.fromMarginal` produced, as the reference displays them. */
void print_states(const std::string& label, const std::vector<std::vector<double>>& rows) {
    std::printf("%s: %zu state(s)\n", label.c_str(), rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        std::printf("  [");
        for (std::size_t j = 0; j < rows[i].size(); ++j)
            std::printf("%s%g", j ? " " : "", rows[i][j]);
        std::printf("]\n");
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// cqn_bas_blocking
// ---------------------------------------------------------------------------

/**
 * Blocking-After-Service: Queue2 holds one job, and a Queue1 completion that
 * finds it full keeps its server until room appears.
 *
 * `sqd` (Smith queue decomposition) is the handler for it; MVA's `default`
 * routes a BAS model there too, and the reference names it explicitly.
 */
void cqn_bas_blocking() {
    Net m("cqn_bas_blocking");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    ClosedClass c1(m, "Class1", 2.0, q1, 0);

    q1.set_service(c1, Exp(1.0));
    q2.set_service(c1, Exp(0.8));
    q2.set_capacity(1.0);
    q1.set_drop_rule(c1, lang::DropStrategy::BAS);

    Routing P;
    cyclic(P, c1, {q1, q2});
    m.link(P);
    const Sn& sn = m.get_struct();

    // The method is presentation: the golden keys this table plain `MVA`.
    section("MVA (method=sqd)", "MVA");
    run_mva(sn, "sqd");
}

// ---------------------------------------------------------------------------
// cqn_bcmp_theorem
// ---------------------------------------------------------------------------

namespace {

/**
 * The one model of `cqn_bcmp_theorem`, at the discipline under test.
 *
 * The point of the example is that PS, FCFS and LCFS-PR give the SAME numbers
 * here -- the BCMP theorem's product-form disciplines -- so the three models
 * differ in exactly one argument and in nothing else.
 */
Net bcmp_model(const std::string& name, SchedStrategy sched) {
    Net m(name);
    Delay delay(m, "Delay");
    Queue queue(m, "Queue1", sched);
    queue.set_number_of_servers(1.0);

    ClosedClass c1(m, "Class1", 2.0, delay, 0);
    ClosedClass c2(m, "Class2", 2.0, delay, 0);

    delay.set_service(c1, Erlang(3.0, 2));
    delay.set_service(c2, HyperExp(0.5, 3.0, 10.0));
    queue.set_service(c1, Exp(1.0));
    queue.set_service(c2, Exp(1.0));

    Routing P;
    cyclic(P, c1, {delay, queue});
    cyclic(P, c2, {delay, queue});
    m.link(P);
    return m;
}

}  // namespace

/** The three product-form disciplines, each solved exactly by CTMC. */
void cqn_bcmp_theorem() {
    const char* names[3] = {"PS scheduling model", "FCFS scheduling model",
                            "LCFS-PR scheduling model"};
    const SchedStrategy scheds[3] = {SchedStrategy::PS, SchedStrategy::FCFS,
                                     SchedStrategy::LCFSPR};
    for (std::size_t k = 0; k < 3; ++k) {
        Net m = bcmp_model(names[k], scheds[k]);
        const Sn& sn = m.get_struct();
        note(std::string("\nMODEL: ") + names[k]);
        // The reference prints the model name and the table, with no solver
        // banner, so the declaration is made to the recorder alone. All three
        // models answer identically -- that IS the BCMP theorem this example
        // demonstrates -- so which of the three the golden holds does not
        // change a digit.
        attribute("CTMC");
        run_ctmc(sn);
    }
}

// ---------------------------------------------------------------------------
// cqn_lcfs_multiclass
// ---------------------------------------------------------------------------

/**
 * Two queues, one LCFS and one LCFS-PR, with one job in each of three classes.
 *
 * Reference: G. Casale, "A family of multiclass LCFS queueing networks with
 * order-dependent product-form solutions", QUESTA 2026. MVA and CTMC must agree
 * here, which is what the example is for.
 */
void cqn_lcfs_multiclass() {
    // mu(i,r), the service RATE of class r at station i: 1./[[1,3,5],[2,4,6]].
    const double mu[2][3] = {{1.0 / 1.0, 1.0 / 3.0, 1.0 / 5.0},
                             {1.0 / 2.0, 1.0 / 4.0, 1.0 / 6.0}};

    Net m("LCFS Multiclass Model");
    Queue q1(m, "Queue1", SchedStrategy::LCFS);
    Queue q2(m, "Queue2", SchedStrategy::LCFSPR);

    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < 3; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1), 1.0, q1, 0));
    for (std::size_t r = 0; r < 3; ++r) {
        q1.set_service(cls[r], Exp(mu[0][r]));
        q2.set_service(cls[r], Exp(mu[1][r]));
    }

    Routing P;
    for (std::size_t r = 0; r < 3; ++r) cyclic(P, cls[r], {q1, q2});
    m.link(P);
    const Sn& sn = m.get_struct();

    note("Solving LCFS+LCFS-PR network with MVA...");
    section("MVA");
    run_mva(sn, "exact");

    note("\nSolving with CTMC for validation...");
    section("CTMC");
    run_ctmc(sn);

    note("\n=== Comparison ===");
    note("MVA and CTMC results should match for this product-form network.");
}

// ---------------------------------------------------------------------------
// cqn_mmpp2_service
// ---------------------------------------------------------------------------

/**
 * A closed network whose Class2 service at Queue1 is an MMPP(2), so the service
 * process is autocorrelated and no product form applies.
 *
 * Class1 routes probabilistically out of the Delay, Class2 uniformly at random
 * over the SAME connections -- including the Delay's self-loop, which is a link
 * the reference adds and which the uniform split therefore counts.
 */
void cqn_mmpp2_service() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);

    ClosedClass c1(m, "Class1", 1.0, delay, 0);
    ClosedClass c2(m, "Class2", 1.0, delay, 0);

    delay.set_service(c1, Erlang(3.0, 2));
    delay.set_service(c2, HyperExp(0.5, 3.0, 10.0));
    q1.set_service(c1, HyperExp(0.1, 1.0, 10.0));
    q1.set_service(c2, mmpp2(1.0, 2.0, 3.0, 4.0));
    q2.set_service(c1, HyperExp(0.1, 1.0, 10.0));
    q2.set_service(c2, Erlang(1.0, 2));

    Routing P;
    P.set(c1, c1, delay, q1, 0.3);
    P.set(c1, c1, delay, q2, 0.7);
    P.set(c1, c1, q1, delay, 1.0);
    P.set(c1, c1, q2, delay, 1.0);
    // RAND needs the CONNECTIONS only; the refresh replaces them by the uniform
    // split, so the Delay's self-link is one of the three destinations.
    m.set_routing(delay, c2, RoutingStrategy::RAND);
    m.set_routing(q1, c2, RoutingStrategy::RAND);
    m.set_routing(q2, c2, RoutingStrategy::RAND);
    P.set(c2, c2, delay, delay, 1.0);
    P.set(c2, c2, delay, q1, 1.0);
    P.set(c2, c2, delay, q2, 1.0);
    P.set(c2, c2, q1, delay, 1.0);
    P.set(c2, c2, q2, delay, 1.0);
    m.link(P);
    m.get_struct();

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000))
    na("LDES", "SolverLDES(model, seed=23000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_multichain_cs
// ---------------------------------------------------------------------------

/**
 * Five closed classes in three chains, two of which switch at the WebServer.
 *
 * The reference prints the CHAIN table, and so does this twin. It used to refuse
 * it -- "this port has no chain getter" -- and print the per-class table
 * instead, which was already stale when it was written: `solver_get_avg_chain`
 * is what `line-cli -a chain` has been calling, and the two tables share no row
 * key, so all 24 of this golden's cells went unmeasured.
 */
void cqn_multichain_cs() {
    Net m("cqn_multichain");
    Delay think(m, "ThinkingTime");
    Queue web(m, "WebServer", SchedStrategy::FCFS);
    Queue app(m, "AppServer", SchedStrategy::PS);
    Queue data(m, "DataServer", SchedStrategy::FCFS);

    const double pop[5] = {2.0, 1.0, 3.0, 2.0, 2.0};
    const char* names[5] = {"HighPriorityA", "HighPriorityB", "RegularA", "RegularB",
                            "Background"};
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < 5; ++r) cls.push_back(m.add_closed_class(names[r], pop[r], think));

    const std::size_t st[4] = {think, web, app, data};
    const double svc[4][5] = {{1.5, 1.6, 2.0, 2.1, 3.0},
                              {0.4, 0.42, 0.5, 0.52, 0.8},
                              {0.7, 0.75, 0.9, 0.95, 1.2},
                              {0.5, 0.52, 0.6, 0.62, 1.0}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t r = 0; r < 5; ++r) m.set_service(st[i], cls[r], D::exp_mean(svc[i][r]));

    Routing P;
    for (std::size_t r = 0; r < 5; ++r) {
        P.set(cls[r], cls[r], think, web, 1.0);
        P.set(cls[r], cls[r], app, data, 1.0);
        P.set(cls[r], cls[r], data, think, 1.0);
    }
    // The class switch happens on the WebServer -> AppServer edge alone.
    P.set(cls[0], cls[0], web, app, 0.7);
    P.set(cls[0], cls[1], web, app, 0.3);
    P.set(cls[1], cls[1], web, app, 0.75);
    P.set(cls[1], cls[0], web, app, 0.25);
    P.set(cls[2], cls[2], web, app, 0.65);
    P.set(cls[2], cls[3], web, app, 0.35);
    P.set(cls[3], cls[3], web, app, 0.7);
    P.set(cls[3], cls[2], web, app, 0.3);
    P.set(cls[4], cls[4], web, app, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("MVA");
    mva::MvaOptions mopt;
    const Matrix<double> minit;
    print_avg_chain(sn, mva::solver_mva_run_analyzer(sn, mopt, minit));

    section("SSA");
    ssa::SsaOptions sopt;
    sopt.samples = 200000;
    sopt.seed = 23000;
    const ssa::SsaSolution sr = ssa::solver_ssa(sn, sopt);
    print_avg_chain(sn, solvers::avg_result_from_sim<double>(sn, sr.QN, sr.UN, sr.RN, sr.TN, sr.CN,
                                                            sr.XN, sr.method));
}

// ---------------------------------------------------------------------------
// cqn_multiserver
// ---------------------------------------------------------------------------

/** Four closed classes over two three-server FCFS queues, with class switching. */
void cqn_multiserver() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    q1.set_number_of_servers(3.0);
    q2.set_number_of_servers(3.0);

    const double pop[4] = {2.0, 2.0, 2.0, 1.0};
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < 4; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1), pop[r], delay, 0));

    delay.set_service(cls[0], Exp(1.0));
    delay.set_service(cls[1], Exp(1.0));
    delay.set_service(cls[2], Exp(10.0));
    delay.set_service(cls[3], Exp(1.0));

    q1.set_service(cls[0], Exp(1.0));
    q1.set_service(cls[1], Erlang(1.0, 2));
    q1.set_service(cls[2], Exp(10.0));
    q1.set_service(cls[3], Exp(1.0));

    q2.set_service(cls[0], Disabled());
    q2.set_service(cls[1], Disabled());
    q2.set_service(cls[2], Erlang(1.0, 2));
    q2.set_service(cls[3], Exp(1.0));

    Routing P;
    P.set(cls[0], cls[0], delay, q1, 0.5);
    P.set(cls[0], cls[0], q1, delay, 1.0);
    P.set(cls[0], cls[0], q2, delay, 1.0);
    P.set(cls[0], cls[1], delay, q1, 0.5);

    P.set(cls[1], cls[1], q1, delay, 1.0);
    P.set(cls[1], cls[1], q2, delay, 1.0);
    P.set(cls[1], cls[0], delay, q1, 1.0);

    P.set(cls[2], cls[2], delay, q1, 0.25);
    P.set(cls[2], cls[2], delay, q2, 0.25);
    P.set(cls[2], cls[2], q1, delay, 1.0);
    P.set(cls[2], cls[2], q2, delay, 1.0);
    P.set(cls[2], cls[3], delay, q1, 0.5);

    P.set(cls[3], cls[3], q1, delay, 1.0);
    P.set(cls[3], cls[3], q2, delay, 1.0);
    P.set(cls[3], cls[2], delay, q1, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    const std::size_t ist = m.station_index(q1);
    std::vector<std::size_t> n = {2, 1, 1, 1}, phases(4, 1);
    for (std::size_t r = 0; r < 4; ++r) phases[r] = sn.phases_of(ist, r + 1);

    note("N/A: State.fromMarginalAndRunning has no C++ counterpart; State.fromMarginal and "
         "State.fromMarginalAndStarted are ported");
    print_states("State (Marginal)", qn::from_marginal(sn, ist, n, phases));
    // The started counts the reference names for this station: two class-1 jobs
    // present with one in service, the other three classes fully in service.
    const std::vector<std::size_t> started = {1, 1, 1, 1};
    print_states("State (MarginalAndStarted)",
                 qn::from_marginal_and_started(sn, ist, n, started, phases));

    section("MVA");
    run_mva(sn);
}

// ---------------------------------------------------------------------------
// cqn_multiserver_nc
// ---------------------------------------------------------------------------

/**
 * A Delay and one three-server FCFS queue, five jobs.
 *
 * With `exact` the multiserver station is rewritten as the load-dependent rate
 * lattice min(n, c) and solved by comomld, which is exact; CTMC is the ground
 * truth it is checked against.
 */
void cqn_multiserver_nc() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    q1.set_number_of_servers(3.0);

    ClosedClass c1(m, "Class1", 5.0, delay, 0);
    delay.set_service(c1, D::exp_mean(1.0));
    q1.set_service(c1, D::exp_mean(0.8));

    Routing P;
    cyclic(P, c1, {delay, q1});
    m.link(P);
    const Sn& sn = m.get_struct();

    section("CTMC");
    run_ctmc(sn);
    section("MVA");
    run_mva(sn);
    section("NC");
    run_nc(sn, "exact");
}

// ---------------------------------------------------------------------------
// cqn_oneline
// ---------------------------------------------------------------------------

/**
 * `Network.cyclicPsInf(N, D, Z)`, written out.
 *
 * The factory is not ported, so its model is built here: one Delay per row of Z
 * and one PS queue per row of D, named `Delay<i>` and `Queue<i>` as the factory
 * names them, cyclically routed, with the first station the reference station of
 * every class and the service rate 1/D(i,r).
 */
void cqn_oneline() {
    const double Dm[2][2] = {{10.0, 5.0}, {5.0, 9.0}};
    const double N[2] = {1.0, 2.0};
    const double Z[2][2] = {{91.0, 92.0}, {93.0, 94.0}};

    Net m("Model");
    std::vector<std::size_t> st;
    st.push_back(m.add_delay("Delay1"));
    st.push_back(m.add_delay("Delay2"));
    st.push_back(m.add_queue("Queue1", SchedStrategy::PS));
    st.push_back(m.add_queue("Queue2", SchedStrategy::PS));

    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < 2; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1), N[r], st[0]));
    for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t i = 0; i < 2; ++i) {
            m.set_service(st[i], cls[r], Exp(1.0 / Z[i][r]));
            m.set_service(st[i + 2], cls[r], Exp(1.0 / Dm[i][r]));
        }

    Routing P;
    for (std::size_t r = 0; r < 2; ++r) cyclic(P, cls[r], st);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("MVA");
    run_mva(sn, "exact");
}

// ---------------------------------------------------------------------------
// cqn_repairmen
// ---------------------------------------------------------------------------

/**
 * The repairman problem: ten machines cycling between operation (the Delay) and
 * repair (the queue), with a 0.7 chance of going straight back to operation.
 */
void cqn_repairmen() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    ClosedClass c1(m, "Class1", 10.0, delay, 0);

    delay.set_service(c1, D::exp_mean(1.0));
    q1.set_service(c1, D::exp_mean(1.5));

    Routing P;
    P.set(c1, c1, delay, delay, 0.7);
    P.set(c1, c1, delay, q1, 0.3);
    P.set(c1, c1, q1, delay, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("CTMC");
    run_ctmc(sn);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("SSA");
    run_ssa(sn, 5000, 23000);
    section("FLD");
    run_fluid(sn);
    section("MVA");
    run_mva(sn);
    section("NC");
    run_nc(sn, "exact");
    section("MAM");
    run_mam(sn);
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000, samples=5000))
    na("LDES", "SolverLDES(model, seed=23000, samples=5000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_repairmen_multi
// ---------------------------------------------------------------------------

/**
 * The repairman problem with two classes and three repairmen.
 *
 * The two MVA rows differ only in `options.config.multiserver`: `softmin` is the
 * smoothed min(n, c) correction and `seidmann` the flow-equivalent one, and the
 * gap between them is the whole point of the comparison.
 */
void cqn_repairmen_multi() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    q1.set_number_of_servers(3.0);

    ClosedClass c1(m, "Class1", 4.0, delay, 0);
    ClosedClass c2(m, "Class2", 2.0, delay, 0);

    delay.set_service(c1, Exp(1.0));
    delay.set_service(c2, Exp(1.0));
    q1.set_service(c1, Exp(1.0));
    q1.set_service(c2, Exp(10.0));

    Routing P;
    cyclic(P, c1, {delay, q1});
    cyclic(P, c2, {delay, q1});
    m.link(P);
    const Sn& sn = m.get_struct();

    // The reference catches each solver's failure and prints it, so a refusal
    // does not cost the rows below it.
    guarded("CTMC", [&]() { run_ctmc(sn); }, "  Error: ");

    // The reference guards these four with `if QNS.isAvailable()`, and so does
    // this: LINE ships no qnsolver binary, so on a host without one the rows are
    // reported absent rather than left to fail.
    if (qns::is_available()) {
        const char* qns_methods[] = {"conway", "reiser", "rolia", "zhou"};
        for (const char* qm : qns_methods) {
            qns::QnsOptions o;
            o.method = qm;
            guarded(std::string("QNS (") + qm + ")",
                    [&]() { print_avg(sn, qns::solver_qns_run_analyzer(sn, o)); }, "  Error: ");
        }
    } else {
        na("QNS", "SolverQNS (conway, reiser, rolia, zhou) shells out to qnsolver, which is not "
                  "installed on this host; LINE ships no copy of it");
    }

    // THE GOLDEN'S `MVA` ROW IS THIS ONE. The banner says which of the two
    // multiserver rules ran, because the reference prints both and they differ
    // by 39% on the Delay queue length; the golden holds softmin, so softmin is
    // what carries the plain key and seidmann keeps its qualified name.
    guarded("MVA (multiserver=softmin)", [&]() { run_mva(sn, "default", "softmin"); },
            "  Error: ", 0, "MVA");
    guarded("MVA (multiserver=seidmann)", [&]() { run_mva(sn, "default", "seidmann"); },
            "  Error: ");
    guarded("NC", [&]() { run_nc(sn); }, "  Error: ");
}

// ---------------------------------------------------------------------------
// cqn_scheduling_dps
// ---------------------------------------------------------------------------

/**
 * PS against DPS on the same two-class population.
 *
 * The weights are set at both queues and matter at ONE of them: PS ignores them
 * by definition, DPS shares the server in proportion to them, so Class2's weight
 * of 5 against Class1's 1 is what separates the two rows.
 */
void cqn_scheduling_dps() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::DPS);

    ClosedClass c1(m, "Class1", 2.0, delay, 0);
    ClosedClass c2(m, "Class2", 1.0, delay, 0);

    note("Network created with mixed scheduling:");
    note("  Delay: Delay node");
    note("  Queue1: PS (Processor Sharing)");
    note("  Queue2: DPS (Discriminatory Processor Sharing)");
    note("\nClass populations: Class1=2, Class2=1");

    delay.set_service(c1, D::exp_mean(1.0 / 3.0));
    delay.set_service(c2, D::exp_mean(1.0 / 0.5));

    q1.set_service(c1, D::exp_mean(1.0 / 0.1));
    q1.set_sched_param(c1, 5.0);
    q1.set_service(c2, D::exp_mean(1.0 / 1.0));
    q1.set_sched_param(c2, 1.0);

    q2.set_service(c1, D::exp_mean(1.0 / 0.1));
    q2.set_sched_param(c1, 1.0);
    q2.set_service(c2, D::exp_mean(1.0 / 1.0));
    q2.set_sched_param(c2, 5.0);

    note("Service parameters configured:");
    note("  Delay: Class1=Exp(3), Class2=Exp(0.5)");
    note("  Queue1 (PS): Class1=Exp(0.1), Class2=Exp(1) - weights ignored");
    note("  Queue2 (DPS): Class1=Exp(0.1,w=1), Class2=Exp(1,w=5) - Class2 has priority");

    Routing P;
    P.set(c1, c1, delay, q1, 0.3);
    P.set(c1, c1, delay, q2, 0.7);
    P.set(c1, c1, q1, delay, 1.0);
    P.set(c1, c1, q2, delay, 1.0);
    P.set(c2, c2, delay, q1, 0.7);
    P.set(c2, c2, delay, q2, 0.3);
    P.set(c2, c2, q1, delay, 1.0);
    P.set(c2, c2, q2, delay, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    note("Routing configured:");
    note("  Class1: 30% Queue1, 70% Queue2 (more traffic to DPS queue)");
    note("  Class2: 70% Queue1, 30% Queue2 (more traffic to PS queue)");
    note("  All queues return 100% to Delay");
    note("=== Multi-Solver Analysis ===");
    note("Configured 5 solvers for comparison");

    // The reference catches a solver that refuses the model and prints this
    // line instead; it catches RuntimeError alone, and std::exception is the
    // nearest C++ equivalent.
    const std::string refused = "Solver not supported for this model: ";
    guarded("CTMC", [&]() { run_ctmc(sn); }, refused);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
    guarded("FLD", [&]() { run_fluid(sn); }, refused);
    guarded("MVA", [&]() { run_mva(sn); }, refused);
    // TODO(cpp): solver_list.append(LDES(model, verbose=True, samples=10000, seed=23000))
    na("LDES", "SolverLDES(model, samples=10000, seed=23000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_threeclass_hyperl
// ---------------------------------------------------------------------------

/**
 * Three classes over a Delay and a two-server PS queue, one of them EMPTY.
 *
 * Class2 carries population zero and exists only as the destination of the class
 * switch out of Class1, so it is a live class with no jobs of its own.
 */
void cqn_threeclass_hyperl() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    q1.set_number_of_servers(2.0);

    ClosedClass c1(m, "Class1", 2.0, delay, 0);
    ClosedClass c2(m, "Class2", 0.0, delay, 0);
    ClosedClass c3(m, "Class3", 1.0, delay, 0);

    delay.set_service(c1, Erlang(3.0, 2));
    delay.set_service(c2, HyperExp(0.5, 3.0, 10.0));
    delay.set_service(c3, Exp(1.0));
    q1.set_service(c1, HyperExp(0.1, 1.0, 10.0));
    q1.set_service(c2, Exp(2.0));
    q1.set_service(c3, Exp(3.0));

    Routing P;
    P.set(c1, c1, delay, delay, 0.3);
    P.set(c1, c1, delay, q1, 0.1);
    P.set(c1, c1, q1, delay, 0.2);
    P.set(c1, c2, delay, delay, 0.6);
    P.set(c1, c2, q1, delay, 0.8);
    P.set(c2, c1, q1, delay, 1.0);
    P.set(c2, c2, delay, q1, 1.0);
    cyclic(P, c3, {delay, q1});
    m.link(P);
    const Sn& sn = m.get_struct();

    // The reference catches each solver's failure, truncates the message to 100
    // characters and moves on, so a refusal does not cost the rows below it.
    guarded("CTMC", [&]() { run_ctmc(sn); }, "Error with CTMC: ", 100);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 5000)));
    guarded("SSA", [&]() { run_ssa(sn, 5000, 23000); }, "Error with SSA: ", 100);
    guarded("FLD", [&]() { run_fluid(sn); }, "Error with FLD: ", 100);
    guarded("MVA", [&]() { run_mva(sn); }, "Error with MVA: ", 100);
    guarded("NC", [&]() { run_nc(sn, "exact"); }, "Error with NC: ", 100);
    guarded("MAM", [&]() { run_mam(sn); }, "Error with MAM: ", 100);
    guarded("", [&]() { run_auto(sn, 23000); }, "Error with LINE: ", 100);
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000, samples=5000, verbose=False))
    na("LDES", "SolverLDES(model, seed=23000, samples=5000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_twoclass_erl
// ---------------------------------------------------------------------------

/**
 * Two classes switched into each other by a ClassSwitch node, dispatched from
 * the delay by round robin (Class1) and weighted round robin (Class2).
 *
 * Only the two simulators are asked for, so both are refused by name. The delay
 * is a Queue with INF scheduling rather than a Delay node, as the reference
 * builds it; `add_queue` leaves the server count at one there, so it is set to
 * infinity on the struct -- `set_number_of_servers` is a deliberate no-op on an
 * INF station and cannot raise it.
 */
void cqn_twoclass_erl() {
    Net m("model");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    Queue delay(m, "Delay", SchedStrategy::INF);
    m.raw_struct().stations[m.station_index(delay) - 1].nservers =
        std::numeric_limits<double>::infinity();

    ClosedClass c1(m, "Class1", 15.0, delay, 0);
    ClosedClass c2(m, "Class2", 5.0, delay, 0);

    Matrix<double> C(2, 2, 0.0);
    C(0, 1) = 1.0;
    C(1, 0) = 1.0;
    const std::size_t cs = m.add_class_switch("CS", C);

    q1.set_service(c1, D::exp_mean(1.5));
    q1.set_service(c2, lang::erlang_fit_mean_order<double>(1.5, 2));
    q2.set_service(c1, lang::erlang_fit_mean_order<double>(1.5, 2));
    q2.set_service(c2, D::exp_mean(1.5));
    delay.set_service(c1, D::exp_mean(1.0));
    delay.set_service(c2, D::exp_mean(1.0));

    Routing P;
    for (std::size_t r = 0; r < 2; ++r) {
        const std::size_t cr = r == 0 ? c1 : c2;
        P.set(cr, cr, q1, cs, 1.0);
        P.set(cr, cr, q2, cs, 1.0);
        P.set(cr, cr, cs, delay, 1.0);
        P.set(cr, cr, delay, q1, 1.0);
        P.set(cr, cr, delay, q2, 1.0);
    }
    m.link(P);

    m.set_routing(cs, c1, RoutingStrategy::RAND);
    m.set_routing(q1, c1, RoutingStrategy::RAND);
    m.set_routing(q2, c1, RoutingStrategy::RAND);
    m.set_routing(delay, c1, RoutingStrategy::RROBIN);

    m.set_routing(cs, c2, RoutingStrategy::RAND);
    m.set_routing(q1, c2, RoutingStrategy::RAND);
    m.set_routing(q2, c2, RoutingStrategy::RAND);
    m.set_routing(delay, c2, RoutingStrategy::WRROBIN);
    std::map<std::size_t, double> w;
    w[q1] = 1.0;
    w[q2] = 2.0;
    m.set_routing_weights(delay, c2, w);
    m.get_struct();

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000, verbose=True))
    na("LDES", "SolverLDES(model, seed=23000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_twoclass_hyperl
// ---------------------------------------------------------------------------

/**
 * Two classes over a Delay and a PS queue, with a class switch out of Class1.
 *
 * The solver list and its options are the ones the JAR test scenario fixes, so
 * the four codebases compare row by row on the same run lengths.
 */
void cqn_twoclass_hyperl() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);

    ClosedClass c1(m, "Class1", 2.0, delay, 0);
    ClosedClass c2(m, "Class2", 2.0, delay, 0);

    delay.set_service(c1, Erlang(3.0, 2));
    delay.set_service(c2, HyperExp(0.5, 3.0, 10.0));
    q1.set_service(c1, HyperExp(0.1, 1.0, 10.0));
    q1.set_service(c2, Exp(1.0));

    Routing P;
    P.set(c1, c1, delay, delay, 0.3);
    P.set(c1, c1, delay, q1, 0.1);
    P.set(c1, c1, q1, delay, 0.2);
    P.set(c1, c2, delay, delay, 0.6);
    P.set(c1, c2, q1, delay, 0.8);
    P.set(c2, c1, q1, delay, 1.0);
    P.set(c2, c2, delay, q1, 1.0);
    m.link(P);
    const Sn& sn = m.get_struct();

    section("CTMC");
    run_ctmc(sn);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 5000)));
    section("SSA");
    run_ssa(sn, 5000, 23000);
    section("FLD");
    run_fluid(sn);
    section("MVA");
    run_mva(sn, "exact");
    section("NC");
    run_nc(sn, "exact");
    section("MAM");
    run_mam(sn);
    // TODO(cpp): solver = np.append(solver, LDES(model, seed=23000, samples=5000))
    na("LDES", "SolverLDES(model, seed=23000, samples=5000) has no C++ engine in this port");
}

// ---------------------------------------------------------------------------
// cqn_twoqueues_multi
// ---------------------------------------------------------------------------

/** A Delay and two FCFS queues in series, ten jobs in each of two classes. */
void cqn_twoqueues_multi() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);

    ClosedClass c1(m, "Class1", 10.0, delay, 0);
    ClosedClass c2(m, "Class2", 10.0, delay, 0);

    const std::size_t st[3] = {delay, q1, q2};
    const double svc[3] = {1.0, 1.5, 3.0};
    for (std::size_t i = 0; i < 3; ++i) {
        m.set_service(st[i], c1, D::exp_mean(svc[i]));
        m.set_service(st[i], c2, D::exp_mean(svc[i]));
    }

    Routing P;
    cyclic(P, c1, {delay, q1, q2});
    cyclic(P, c2, {delay, q1, q2});
    m.link(P);
    const Sn& sn = m.get_struct();

    section("MVA");
    run_mva(sn);
}

LINE_EXAMPLE("basic/closedQN", cqn_bas_blocking);
LINE_EXAMPLE("basic/closedQN", cqn_bcmp_theorem);
LINE_EXAMPLE("basic/closedQN", cqn_lcfs_multiclass);
LINE_EXAMPLE("basic/closedQN", cqn_mmpp2_service);
LINE_EXAMPLE("basic/closedQN", cqn_multichain_cs);
LINE_EXAMPLE("basic/closedQN", cqn_multiserver);
LINE_EXAMPLE("basic/closedQN", cqn_multiserver_nc);
LINE_EXAMPLE("basic/closedQN", cqn_oneline);
LINE_EXAMPLE("basic/closedQN", cqn_repairmen);
LINE_EXAMPLE("basic/closedQN", cqn_repairmen_multi);
LINE_EXAMPLE("basic/closedQN", cqn_scheduling_dps);
LINE_EXAMPLE("basic/closedQN", cqn_threeclass_hyperl);
LINE_EXAMPLE("basic/closedQN", cqn_twoclass_erl);
LINE_EXAMPLE("basic/closedQN", cqn_twoclass_hyperl);
LINE_EXAMPLE("basic/closedQN", cqn_twoqueues_multi);

}  // namespace examples
}  // namespace line
