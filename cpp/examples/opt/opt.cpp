/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/opt/opt_*.py`: the nine optimisation tutorials.
 *
 * Each tutorial builds and solves the same optimization problem as its MATLAB
 * and native-Python counterpart through the public C++ LineOpt surface.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/opt/bisection_solver.h"
#include "line/opt/decomposition.h"
#include "line/opt/layered_variables.h"
#include "line/opt/line_opt_solver.h"
#include "line/opt/pareto.h"
#include "line/opt/sensitivity.h"

namespace line {
namespace examples {

namespace {

using Lqn = lqn::LqnBuilder<double>;

void variable(const std::string& text) { std::printf("  variable  : %s\n", text.c_str()); }
void objective(const std::string& text) { std::printf("  objective : %s\n", text.c_str()); }
void constraint(const std::string& text) { std::printf("  constraint: %s\n", text.c_str()); }

/** The model the tutorial optimises: its stations, servers and service rates. */
void describe(Net& m) {
    const Sn& sn = m.get_struct();
    std::printf("MODEL: %s\n", sn.name.c_str());
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        std::printf("  station %-12s servers=%g", sn.stations[i].name.c_str(),
                    sn.stations[i].nservers);
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (sn.disabled[i][r]) continue;
            std::printf("  %s: rate=%g", sn.classes[r].name.c_str(), sn.rates(i, r));
        }
        std::printf("\n");
    }
    for (std::size_t r = 0; r < sn.nclasses; ++r)
        std::printf("  class   %-12s population=%g\n", sn.classes[r].name.c_str(),
                    sn.classes[r].population);
}

/** Source -> Queue -> Sink with one open class, the M/M/c of five tutorials. */
Net mmc(const std::string& model_name, double arrival_rate, double service_rate) {
    Net m(model_name);
    Source source(m, "Arrivals");
    Queue queue(m, "Server", SchedStrategy::FCFS);
    Sink sink(m, "Departures");
    OpenClass jobs(m, "Jobs");
    source.set_arrival(jobs, Exp(arrival_rate));
    queue.set_service(jobs, Exp(service_rate));
    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Tutorial 1: server sizing with differential evolution
// ---------------------------------------------------------------------------

void opt_server_sizing() {
    Net m = mmc("MMc", 3.0, 1.0);
    describe(m);
    variable("ServerAllocation(Server) in [1, 10]");
    objective("MinimizeCost, 10.0 per server of Server");
    constraint("Utilization(Server) <= 0.5");
    note("  solver    : differential evolution, seed 42");
    note("  reference answer: 6 servers, cost 60.0");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::ServerAllocation>("Server", 1, 10))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}}))
        .add_constraint(std::make_shared<opt::UtilizationConstraint>("Server", 0.5));
    opt::LineOptSolverOptions options;
    options.set_seed(42).set_max_iterations(60);
    const opt::OptimizationResult result = opt::LineOptSolver(problem, options).solve();
    std::printf("  result    : servers=%g cost=%g feasible=%s evaluations=%zu\n",
                opt::scalar_value(result.variable_values.at("Server_servers")),
                result.objective_value, result.feasible ? "true" : "false",
                result.model_evaluations);
}

// ---------------------------------------------------------------------------
// Tutorial 2: the same sizing, solved exactly by bisection
// ---------------------------------------------------------------------------

void opt_bisection_sizing() {
    Net m = mmc("MMc", 3.0, 1.0);
    describe(m);
    variable("ServerAllocation(Server) in [1, 10]");
    objective("MinimizeCost, 10.0 per server of Server");
    constraint("Utilization(Server) <= 0.5");
    note("  solver    : BisectionSolver, monotone feasibility in the server count");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::ServerAllocation>("Server", 1, 10))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}}))
        .add_constraint(std::make_shared<opt::UtilizationConstraint>("Server", 0.5));
    const opt::OptimizationResult result = opt::BisectionSolver(problem).solve();
    std::printf("  result    : servers=%g cost=%g feasible=%s evaluations=%zu\n",
                opt::scalar_value(result.variable_values.at("Server_servers")),
                result.objective_value, result.feasible ? "true" : "false",
                result.model_evaluations);
}

// ---------------------------------------------------------------------------
// Tutorial 3: continuous service-rate tuning
// ---------------------------------------------------------------------------

void opt_service_rate() {
    Net m = mmc("MM1", 3.0, 4.0);
    describe(m);
    variable("ServiceRate(Server, Jobs) in [3.5, 8.0]");
    objective("MinimizeCost, 20.0 per unit of rate at Server");
    constraint("ResponseTime(Server, Jobs) <= 0.5");
    note("  solver    : differential evolution, seed 42, 60 iterations");
    note("  reference answer: rate 5.0, cost 100.0");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::ServiceRate>(
               "Server", "Jobs", 3.5, 8.0))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>(),
            std::map<std::string, double>{{"Server", 20.0}}))
        .add_constraint(std::make_shared<opt::ResponseTimeConstraint>(
            "Server", "Jobs", 0.5));
    opt::LineOptSolverOptions options;
    options.set_seed(42).set_max_iterations(60);
    const opt::OptimizationResult result = opt::LineOptSolver(problem, options).solve();
    std::printf("  result    : rate=%g cost=%g feasible=%s\n",
                opt::scalar_value(result.variable_values.at("Server_Jobs_rate")),
                result.objective_value, result.feasible ? "true" : "false");
}

// ---------------------------------------------------------------------------
// Tutorial 4: load balancing across heterogeneous servers
// ---------------------------------------------------------------------------

void opt_load_balancing() {
    Net m("LoadBalance");
    Source source(m, "S");
    Queue fast(m, "Fast", SchedStrategy::PS);
    Queue slow(m, "Slow", SchedStrategy::PS);
    Sink sink(m, "K");
    OpenClass jobs(m, "Jobs");
    source.set_arrival(jobs, Exp(2.0));
    fast.set_service(jobs, Exp(4.0));
    slow.set_service(jobs, Exp(2.5));
    Routing P;
    // Two outgoing links and no probabilities set: the reference's default RAND.
    P.set(jobs, jobs, source, fast, 0.5);
    P.set(jobs, jobs, source, slow, 0.5);
    P.set(jobs, jobs, fast, sink, 1.0);
    P.set(jobs, jobs, slow, sink, 1.0);
    m.link(P);

    describe(m);
    variable("RoutingProbabilities(Jobs, from S, to [Fast, Slow])");
    objective("MinimizeSystemResponseTime(Jobs)");
    note("  solver    : differential evolution, seed 42, 60 iterations");
    note("  reference answer: 0.743 to Fast, system response time 0.4264");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::RoutingProbabilities>(
               "Jobs", "S", std::vector<std::string>{"Fast", "Slow"}))
        .set_objective(std::make_shared<opt::MinimizeSystemResponseTime>("Jobs"));
    opt::LineOptSolverOptions options;
    options.set_seed(42).set_max_iterations(60);
    const opt::OptimizationResult result = opt::LineOptSolver(problem, options).solve();
    const opt::Value& probabilities = result.variable_values.at("Jobs_routing_from_S");
    std::printf("  result    : fast=%g slow=%g response=%g feasible=%s\n",
                probabilities[0], probabilities[1], result.objective_value,
                result.feasible ? "true" : "false");
}

// ---------------------------------------------------------------------------
// Tutorial 5: population sizing in a closed network
// ---------------------------------------------------------------------------

void opt_population_sizing() {
    Net m("Interactive");
    Delay think(m, "Think");
    Queue server(m, "AppServer", SchedStrategy::PS);
    ClosedClass users(m, "Users", 1, think);
    think.set_service(users, Exp(1.0));
    server.set_service(users, Exp(2.0));
    Routing P;
    cyclic(P, users, {think, server});
    m.link(P);

    describe(m);
    variable("JobPopulation(Users) in [1, 50]");
    objective("MinimizeCost (pure feasibility sizing)");
    constraint("ResponseTime(AppServer, Users) <= 2.0");
    note("  solver    : BisectionSolver, direction max_feasible");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::JobPopulation>("Users", 1, 50))
        .set_objective(std::make_shared<opt::MinimizeCost>())
        .add_constraint(std::make_shared<opt::ResponseTimeConstraint>(
            "AppServer", "Users", 2.0));
    const opt::OptimizationResult result =
        opt::BisectionSolver(problem, "max_feasible").solve();
    std::printf("  result    : users=%g feasible=%s evaluations=%zu\n",
                opt::scalar_value(result.variable_values.at("Users_population")),
                result.feasible ? "true" : "false", result.model_evaluations);
}

// ---------------------------------------------------------------------------
// Tutorial 6: robust sizing under a peak-load scenario
// ---------------------------------------------------------------------------

void opt_robust_sizing() {
    Net base = mmc("MMc", 3.0, 1.0);
    Net peak = mmc("MMc", 4.5, 1.0);
    note("base scenario (average load)");
    describe(base);
    note("peak scenario");
    describe(peak);
    variable("ServerAllocation(Server) in [1, 12]");
    objective("MinimizeCost, 10.0 per server of Server");
    constraint("Utilization(Server) <= 0.5 on EVERY scenario");
    note("  solver    : BisectionSolver over the base model plus one scenario");
    note("  reference answer: 9 servers, against 6 for the average load alone");
    opt::OptimizationProblem problem(base);
    problem.add_variable(std::make_shared<opt::ServerAllocation>("Server", 1, 12))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}}))
        .add_constraint(std::make_shared<opt::UtilizationConstraint>("Server", 0.5))
        .add_scenario(peak);
    const opt::OptimizationResult result = opt::BisectionSolver(problem).solve();
    std::printf("  result    : servers=%g cost=%g feasible=%s\n",
                opt::scalar_value(result.variable_values.at("Server_servers")),
                result.objective_value, result.feasible ? "true" : "false");
}

// ---------------------------------------------------------------------------
// Tutorial 7: the cost-performance Pareto frontier
// ---------------------------------------------------------------------------

void opt_pareto_frontier() {
    Net m = mmc("MMc", 3.0, 1.0);
    describe(m);
    variable("ServerAllocation(Server) in [1, 12]");
    objective("MinimizeCost, 10.0 per server of Server");
    constraint("Utilization(Server) <= eps, eps in {0.3, 0.4, 0.5, 0.6, 0.75, 0.9}");
    note("  solver    : ParetoSweep, epsilon-constraint method, bisection per point");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::ServerAllocation>("Server", 1, 12))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}}));
    opt::ParetoSweep sweep(problem,
        [](double epsilon) {
            return std::make_shared<opt::UtilizationConstraint>("Server", epsilon);
        }, {0.3, 0.4, 0.5, 0.6, 0.75, 0.9}, "bisection");
    sweep.solve();
    for (const opt::ParetoPoint& point : sweep.frontier())
        std::printf("  frontier  : epsilon=%g cost=%g servers=%g\n", point.epsilon,
                    point.objective_value,
                    opt::scalar_value(point.result.variable_values.at("Server_servers")));
}

// ---------------------------------------------------------------------------
// Tutorial 8: block coordinate descent over two variable types
// ---------------------------------------------------------------------------

void opt_decomposition() {
    Net m = mmc("MMc", 3.0, 1.0);
    describe(m);
    variable("ServerAllocation(Server) in [1, 10]");
    variable("ServiceRate(Server, Jobs) in [1.0, 4.0]");
    objective("MinimizeCost, 10.0 per server and 20.0 per unit of rate at Server");
    constraint("ResponseTime(Server, Jobs) <= 0.5");
    note("  solver    : DecompositionWorkflow, auto_decompose, 4 cycles, tolerance 1e-3");
    opt::OptimizationProblem problem(m);
    problem.add_variable(std::make_shared<opt::ServerAllocation>("Server", 1, 10))
        .add_variable(std::make_shared<opt::ServiceRate>("Server", "Jobs", 1.0, 4.0))
        .set_objective(std::make_shared<opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}},
            std::map<std::string, double>{{"Server", 20.0}}))
        .add_constraint(std::make_shared<opt::ResponseTimeConstraint>(
            "Server", "Jobs", 0.5));
    opt::LineOptSolverOptions options;
    options.set_seed(42).set_max_iterations(30).set_popsize(10);
    opt::DecompositionWorkflow workflow(problem);
    workflow.auto_decompose().set_solver_options(options);
    const opt::WorkflowResult result = workflow.solve_sequential(4, 1e-3);
    std::printf("  result    : converged=%s cycles=%zu objective=%g\n",
                result.converged ? "true" : "false", result.cycles_completed,
                result.final_objective);
}

// ---------------------------------------------------------------------------
// Tutorial 9: LQN host-demand tuning with layer freezing
// ---------------------------------------------------------------------------

namespace {

/** The reference's `build()`: the three-task client/server LQN it tunes. */
Lqn lqn_hostdemand_model() {
    Lqn m;
    m.processor("P1", 2, SchedStrategy::PS);
    m.processor("P2", 3, SchedStrategy::PS);
    // The reference's think times are Exp(1/2), Exp(1/3) and Exp(1/4) BY RATE.
    m.task("T1", 50, SchedStrategy::REF, "P1");
    m.think_time("T1", Exp(1.0 / 2.0));
    m.task("T2", 50, SchedStrategy::FCFS, "P1");
    m.think_time("T2", Exp(1.0 / 3.0));
    m.task("T3", 25, SchedStrategy::FCFS, "P2");
    m.think_time("T3", Exp(1.0 / 4.0));
    m.entry("E1", "T1");
    m.entry("E2", "T2");
    m.entry("E3", "T3");
    m.activity("AS1", Exp(10.0), "T1");
    m.bound_to("AS1", "E1");
    m.sync_call("AS1", "E2", 1.0);
    m.activity("AS2", Exp(20.0), "T2");
    m.bound_to("AS2", "E2");
    m.sync_call("AS2", "E3", 5.0);
    m.replies_to("AS2", "E2");
    m.activity("AS3", Exp(50.0), "T3");
    m.bound_to("AS3", "E3");
    m.replies_to("AS3", "E3");
    return m;
}

}  // namespace

void opt_lqn_hostdemand() {
    const Lqn m = lqn_hostdemand_model();
    const lqn::LqnStruct<double> ln = m.build();
    std::printf("MODEL: LQN-Basic, %zu processors, %zu tasks, %zu entries, %zu activities\n",
                ln.nhosts, ln.ntasks, ln.nentries, ln.nacts);
    variable("HostDemand(AS1) in [0.02, 0.2]");
    variable("HostDemand(AS2) in [0.01, 0.1]");
    objective("MinimizeSystemResponseTime(T1)");
    constraint("Utilization(P1) <= 0.95");
    note("  solver    : gradient, lqn_gradient in {fd, partial_sens, partial_plus_fd}, seed 1");
    note("  then      : frozen_layers=[P2], and the layered workflow with auto_freeze");
    const auto make_problem = [&m]() {
        opt::OptimizationProblem problem(m.model());
        problem.add_variable(std::make_shared<opt::HostDemand>("AS1", 0.02, 0.2))
            .add_variable(std::make_shared<opt::HostDemand>("AS2", 0.01, 0.1))
            .set_objective(std::make_shared<opt::MinimizeSystemResponseTime>("T1"))
            .add_constraint(std::make_shared<opt::UtilizationConstraint>("P1", 0.95));
        return problem;
    };
    for (const std::string mode : {"fd", "partial_sens", "partial_plus_fd"}) {
        opt::OptimizationProblem problem = make_problem();
        opt::LineOptSolverOptions options;
        options.set_seed(1).set_optimizer("gradient").set_lqn_gradient(mode)
            .set_max_iterations(8).set_gradient_restarts(1).set_time_limit(180.0);
        const opt::OptimizationResult result = opt::LineOptSolver(problem, options).solve();
        std::printf("  %-15s objective=%g feasible=%s evaluations=%zu\n", mode.c_str(),
                    result.objective_value, result.feasible ? "true" : "false",
                    result.model_evaluations);
    }

    opt::OptimizationProblem frozen_problem = make_problem();
    opt::LineOptSolverOptions frozen_options;
    frozen_options.set_seed(1).set_optimizer("gradient").set_lqn_gradient("fd")
        .set_frozen_layers({"P2"}).set_max_iterations(8).set_gradient_restarts(1);
    const opt::OptimizationResult frozen =
        opt::LineOptSolver(frozen_problem, frozen_options).solve();
    std::printf("  frozen P2      objective=%g feasible=%s\n", frozen.objective_value,
                frozen.feasible ? "true" : "false");

    opt::OptimizationProblem layered_problem = make_problem();
    opt::DecompositionWorkflow workflow(layered_problem);
    opt::LineOptSolverOptions workflow_options;
    workflow_options.set_seed(1).set_optimizer("gradient").set_lqn_gradient("fd")
        .set_max_iterations(8).set_gradient_restarts(1);
    workflow.set_solver_options(workflow_options);
    const opt::WorkflowResult layered = workflow.solve_layered(6, 1e-3, true, 1e-2);
    std::printf("  layered        objective=%g converged=%s cycles=%zu\n",
                layered.final_objective, layered.converged ? "true" : "false",
                layered.cycles_completed);
}

void opt_sensitivity_report() {
    Net model("Tandem");
    Source source(model, "Arrivals");
    Queue frontend(model, "Frontend", SchedStrategy::FCFS);
    Queue backend(model, "Backend", SchedStrategy::FCFS);
    Sink sink(model, "Departures");
    OpenClass jobs(model, "Jobs");
    source.set_arrival(jobs, Exp(1.0));
    frontend.set_service(jobs, Exp(2.0));
    backend.set_service(jobs, Exp(1.5));
    Routing routing;
    serial(routing, {source, frontend, backend, sink});
    model.link(routing);

    opt::LineEvaluator evaluator(model, {});
    const opt::EvaluationResult result = evaluator.evaluate({});
    std::printf("Solver used : %s\n", result.solver_used.c_str());
    std::printf("Feasible    : %s\n", result.feasible ? "true" : "false");
    std::printf("Frontend Q  : %.4f\n", result.queue_length("Frontend"));
    std::printf("Backend  Q  : %.4f\n", result.queue_length("Backend"));
    const auto sensitivity = opt::compute_model_sensitivities(model);
    if (!sensitivity) {
        note("No analytic sensitivities available for this model.");
        return;
    }
    for (const std::string kind : {"QLen", "RespT", "Util", "Tput"}) {
        const opt::SensitivityData::MetricMap* metrics = sensitivity->for_kind(kind);
        if (!metrics) continue;
        std::printf("\n=== d(%s)/d(rate) ===\n", kind.c_str());
        for (const auto& metric : *metrics)
            for (const auto& parameter : metric.second)
                std::printf("  %-24s %-32s % .6f\n", metric.first.c_str(),
                            parameter.first.c_str(), parameter.second);
    }
}

LINE_EXAMPLE("opt", opt_server_sizing);
LINE_EXAMPLE("opt", opt_bisection_sizing);
LINE_EXAMPLE("opt", opt_service_rate);
LINE_EXAMPLE("opt", opt_load_balancing);
LINE_EXAMPLE("opt", opt_population_sizing);
LINE_EXAMPLE("opt", opt_robust_sizing);
LINE_EXAMPLE("opt", opt_pareto_frontier);
LINE_EXAMPLE("opt", opt_decomposition);
// ---------------------------------------------------------------------------
// opt_lqn_host_demand
// ---------------------------------------------------------------------------

namespace {

/** A three-tier bookstore: client thinks, app server renders, database queries. */
Lqn bookstore_model() {
    Lqn b;
    b.processor("ClientCPU", 1, SchedStrategy::INF);
    b.task("Client", 4, SchedStrategy::REF, "ClientCPU");
    b.think_time("Client", Exp(1.0));
    b.entry("Browse", "Client");
    b.processor("AppCPU", 1, SchedStrategy::PS);
    b.task("AppServer", 1, SchedStrategy::FCFS, "AppCPU");
    b.entry("Render", "AppServer");
    b.processor("DbCPU", 1, SchedStrategy::PS);
    b.task("Database", 1, SchedStrategy::FCFS, "DbCPU");
    b.entry("Query", "Database");

    b.activity("BrowseAct", Exp(2.0), "Client");
    b.bound_to("BrowseAct", "Browse");
    b.sync_call("BrowseAct", "Render", 1.0);
    b.activity("RenderAct", Exp(4.0), "AppServer");
    b.bound_to("RenderAct", "Render");
    b.sync_call("RenderAct", "Query", 1.0);
    b.replies_to("RenderAct", "Render");
    b.activity("QueryAct", Exp(5.0), "Database");
    b.bound_to("QueryAct", "Query");
    b.replies_to("QueryAct", "Query");
    return b;
}

}  // namespace

/**
 * Tune two host demands against a SPENDING budget rather than a utilization cap.
 *
 * This is a different question from `opt_lqn_hostdemand`, which caps a processor
 * and asks for the fastest system: here the two demands are bought from one
 * budget, so the optimizer trades them against each other. Faster is always
 * better for response time, so the optimum sits at the lower bound of both and
 * the budget is slack -- the constraint binds only once the bounds are widened,
 * which is what makes it a budget and not a cap.
 */
void opt_lqn_host_demand() {
    const Lqn m = bookstore_model();
    const lqn::LqnStruct<double> ln = m.build();
    std::printf("MODEL: BookstoreLQN, %zu processors, %zu tasks, %zu entries, %zu activities\n",
                ln.nhosts, ln.ntasks, ln.nentries, ln.nacts);
    variable("HostDemand(RenderAct) in [0.05, 0.50]");
    variable("HostDemand(QueryAct) in [0.05, 0.50]");
    objective("MinimizeSystemResponseTime(Client)");
    constraint("Budget: 1.0*RenderAct + 1.0*QueryAct <= 0.45");
    note("  solver    : gradient, lqn_gradient partial_plus_fd, fd_refresh 3, seed 7");

    opt::OptimizationProblem problem(m.model());
    std::map<std::string, double> cost;
    cost["RenderAct_hostdemand"] = 1.0;
    cost["QueryAct_hostdemand"] = 1.0;
    problem.add_variable(std::make_shared<opt::HostDemand>("RenderAct", 0.05, 0.50))
        .add_variable(std::make_shared<opt::HostDemand>("QueryAct", 0.05, 0.50))
        .set_objective(std::make_shared<opt::MinimizeSystemResponseTime>("Client"))
        .add_constraint(std::make_shared<opt::BudgetConstraint>(0.45, cost));

    opt::LineOptSolverOptions options;
    options.set_seed(7)
        .set_optimizer("gradient")
        .set_lqn_gradient("partial_plus_fd")
        .set_fd_refresh(3)
        .set_max_iterations(8);
    const opt::OptimizationResult result = opt::LineOptSolver(problem, options).solve();

    std::printf("Render demand   : %.4f\n",
                opt::scalar_value(result.variable_values.at("RenderAct_hostdemand")));
    std::printf("Query demand    : %.4f\n",
                opt::scalar_value(result.variable_values.at("QueryAct_hostdemand")));
    std::printf("System RespT    : %.4f\n", result.objective_value);
    std::printf("Feasible        : %d\n", result.feasible ? 1 : 0);
    std::printf("Total violation : %.4f\n", result.total_violation());
    std::printf("LINE evaluations: %zu\n", result.model_evaluations);
}

LINE_EXAMPLE("opt", opt_lqn_hostdemand);
LINE_EXAMPLE("opt", opt_lqn_host_demand);
LINE_EXAMPLE("opt", opt_sensitivity_report);

}  // namespace examples
}  // namespace line
