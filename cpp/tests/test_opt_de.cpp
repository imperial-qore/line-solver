/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#include "doctest.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/lang/lqn/lqn_builder.h"
#include "line/opt/de/differential_evolution.h"
#include "line/opt/de/numpy_random_state.h"
#include "line/opt/bisection_solver.h"
#include "line/opt/decomposition.h"
#include "line/opt/line_opt_solver.h"
#include "line/opt/layered_variables.h"
#include "line/opt/pareto.h"
#include "line/opt/sensitivity.h"

namespace {

using Json = nlohmann::json;

Json rng_golden() {
    std::ifstream in(std::string(LINE_MP_REPO_ROOT) +
                     "/jar/src/test/resources/opt/opt_rng_golden.json");
    REQUIRE_MESSAGE(bool(in), "jar/src/test/resources/opt/opt_rng_golden.json is missing");
    Json out;
    in >> out;
    return out;
}

const std::array<std::uint64_t, 4> seeds = {{0, 1, 42, 12345}};

}  // namespace

TEST_CASE("line-opt MT19937 reproduces NumPy's seeded state") {
    const Json golden = rng_golden().at("seed_state");
    for (std::uint64_t seed : seeds) {
        line::opt::de::NumpyRandomState rng(seed);
        const auto& key = rng.generator().state_key();
        const Json& expected = golden.at(std::to_string(seed));
        for (std::size_t i = 0; i < 8; ++i)
            CHECK(static_cast<std::uint64_t>(key[i]) == expected.at("first8").at(i).get<std::uint64_t>());
        for (std::size_t i = 0; i < 4; ++i)
            CHECK(static_cast<std::uint64_t>(key[620 + i]) ==
                  expected.at("last4").at(i).get<std::uint64_t>());
    }
}

TEST_CASE("line-opt RandomState reproduces NumPy draws and word consumption") {
    const Json golden = rng_golden();

    SUBCASE("random_sample") {
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("random_sample").at(std::to_string(seed));
            for (const Json& value : expected)
                CHECK(rng.random_sample() == value.get<double>());
        }
    }

    SUBCASE("randint") {
        const std::array<std::uint64_t, 12> high = {{2, 3, 5, 10, 16, 17, 100, 1000, 7, 13, 64, 65}};
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("randint").at(std::to_string(seed));
            for (std::size_t i = 0; i < high.size(); ++i)
                CHECK(rng.randint(0, high[i]) == expected.at(i).get<std::uint64_t>());
        }
    }

    SUBCASE("shuffle") {
        const std::array<std::size_t, 7> sizes = {{2, 3, 4, 5, 7, 10, 17}};
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("shuffle").at(std::to_string(seed));
            for (std::size_t n : sizes) {
                std::vector<std::size_t> values(n);
                for (std::size_t i = 0; i < n; ++i) values[i] = i;
                rng.shuffle(values);
                for (std::size_t i = 0; i < n; ++i)
                    CHECK(values[i] == expected.at(std::to_string(n)).at(i).get<std::size_t>());
            }
        }
    }

    SUBCASE("permutation") {
        const std::array<std::size_t, 4> sizes = {{3, 5, 7, 10}};
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("permutation").at(std::to_string(seed));
            for (std::size_t n : sizes) {
                const std::vector<std::size_t> values = rng.permutation(n);
                for (std::size_t i = 0; i < n; ++i)
                    CHECK(values[i] == expected.at(std::to_string(n)).at(i).get<std::size_t>());
            }
        }
    }

    SUBCASE("uniform dither") {
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("uniform_dither").at(std::to_string(seed));
            for (const Json& value : expected)
                CHECK(rng.uniform(0.5, 1.0) == value.get<double>());
        }
    }

    SUBCASE("uniform grid uses C-order fill") {
        for (std::uint64_t seed : seeds) {
            line::opt::de::NumpyRandomState rng(seed);
            const Json& expected = golden.at("uniform_grid").at(std::to_string(seed));
            const std::size_t total = expected.at("shape").at(0).get<std::size_t>() *
                                      expected.at("shape").at(1).get<std::size_t>();
            const std::vector<double> values = rng.uniform(0.0, 1.0, total);
            for (std::size_t i = 0; i < total; ++i)
                CHECK(values[i] == expected.at("flat").at(i).get<double>());
        }
    }
}

TEST_CASE("line-opt Differential Evolution reproduces the MATLAB SciPy trajectory") {
    const Json golden = rng_golden().at("de_trajectory");
    const auto sphere = [](const std::vector<double>& x) {
        double value = 0.0;
        for (double xi : x) value += (xi - 0.3) * (xi - 0.3);
        return value;
    };

    for (std::uint64_t seed : {UINT64_C(42), UINT64_C(7)}) {
        line::opt::de::DifferentialEvolution solver(
            sphere, {-1.0, -1.0, 0.0}, {1.0, 1.0, 2.0}, "best1bin", 5, 8,
            0.5, 1.0, 0.7, 0.0, seed);
        const line::opt::de::DifferentialEvolutionResult result = solver.solve();
        const Json& expected = golden.at(std::to_string(seed));
        REQUIRE(result.x.size() == 3);
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(result.x[j] == doctest::Approx(expected.at("x").at(j).get<double>()).epsilon(1e-12));
        CHECK(result.fun == doctest::Approx(expected.at("fun").get<double>()).epsilon(1e-12));
        REQUIRE(result.best_per_generation.size() == expected.at("best_per_gen").size());
        for (std::size_t i = 0; i < result.best_per_generation.size(); ++i)
            for (std::size_t j = 0; j < result.best_per_generation[i].size(); ++j)
                CHECK(result.best_per_generation[i][j] ==
                      doctest::Approx(expected.at("best_per_gen").at(i).at(j).get<double>())
                          .epsilon(1e-12));
    }
}

TEST_CASE("line-opt flat variables and objective arithmetic match MATLAB") {
    line::qn::Network<double> model("mm1");
    auto src=model.add_source("Src"),q=model.add_queue("Server",line::lang::SchedStrategy::FCFS),sink=model.add_sink("Snk");
    auto jobs=model.add_open_class("Jobs");model.set_arrival(src,jobs,line::lang::Distrib<double>::exp_rate(1));model.set_service(q,jobs,line::lang::Distrib<double>::exp_rate(2));line::qn::RoutingMatrix<double>P;P.set(jobs,jobs,src,q,1);P.set(jobs,jobs,q,sink,1);model.link(P);
    line::opt::ServerAllocation servers("Server",1,8);CHECK(line::opt::scalar_value(servers.decode({0}))==1);CHECK(line::opt::scalar_value(servers.decode({1}))==8);servers.apply(model,{5});CHECK(model.get_struct().stations[1].nservers==5);
    line::opt::ServiceRate rate("Server","Jobs",.5,5);rate.apply(model,{3});CHECK(model.get_struct().rates(1,0)==doctest::Approx(3));
    line::opt::EvaluationResult er;er.set_response_time("Server","Jobs",.8);er.set_throughput("Server","Jobs",5);er.utilizations["Server"]=.75;line::opt::VariableValues vv{{"Server_servers",{3}}};
    CHECK(line::opt::ResponseTimeConstraint("Server","Jobs",.5).evaluate(er,vv)==doctest::Approx(.3));CHECK(line::opt::UtilizationConstraint("Server",.9).evaluate(er,vv)==0);CHECK(line::opt::ThroughputConstraint("Server","",10).evaluate(er,vv)==5);
    line::opt::MinimizeCost cost({{"Server",100}});cost.constraints.push_back(std::make_shared<line::opt::ResponseTimeConstraint>("Server","Jobs",.5));CHECK(cost.evaluate(er,vv)==300);CHECK(cost.evaluate_with_penalty(er,vv)==doctest::Approx(300300));
}

TEST_CASE("line-opt bisection sizes the MATLAB M/M/c fixture") {
    line::qn::Network<double> model("MMc");auto src=model.add_source("Arrivals"),q=model.add_queue("Server",line::lang::SchedStrategy::FCFS),sink=model.add_sink("Departures");auto jobs=model.add_open_class("Jobs");model.set_arrival(src,jobs,line::lang::Distrib<double>::exp_rate(3));model.set_service(q,jobs,line::lang::Distrib<double>::exp_rate(1));line::qn::RoutingMatrix<double>P;P.set(jobs,jobs,src,q,1);P.set(jobs,jobs,q,sink,1);model.link(P);
    line::opt::OptimizationProblem problem(model);problem.add_variable(std::make_shared<line::opt::ServerAllocation>("Server",1,10)).set_objective(std::make_shared<line::opt::MinimizeCost>(std::map<std::string,double>{{"Server",10}})).add_constraint(std::make_shared<line::opt::UtilizationConstraint>("Server",.5));
    auto result=line::opt::BisectionSolver(problem).solve();CHECK(line::opt::scalar_value(result.variable_values.at("Server_servers"))==6);CHECK(result.objective_value==doctest::Approx(60));CHECK(result.feasible);CHECK(result.model_evaluations==4);

    line::opt::OptimizationProblem continuous(model);
    continuous.add_variable(std::make_shared<line::opt::ServiceRate>("Server", "Jobs", 0.5, 5.0))
        .set_objective(std::make_shared<line::opt::MinimizeCost>());
    CHECK_THROWS_AS(line::opt::BisectionSolver{continuous}, line::InputError);
}

TEST_CASE("line-opt main solver sizes the MATLAB M/M/c fixture") {
    line::qn::Network<double> model("MMc");
    auto src = model.add_source("Arrivals");
    auto queue = model.add_queue("Server", line::lang::SchedStrategy::FCFS);
    auto sink = model.add_sink("Departures");
    auto jobs = model.add_open_class("Jobs");
    model.set_arrival(src, jobs, line::lang::Distrib<double>::exp_rate(3));
    model.set_service(queue, jobs, line::lang::Distrib<double>::exp_rate(1));
    line::qn::RoutingMatrix<double> routing;
    routing.set(jobs, jobs, src, queue, 1);
    routing.set(jobs, jobs, queue, sink, 1);
    model.link(routing);

    line::opt::OptimizationProblem problem(model);
    problem.add_variable(std::make_shared<line::opt::ServerAllocation>("Server", 1, 10))
        .set_objective(std::make_shared<line::opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10}}))
        .add_constraint(std::make_shared<line::opt::UtilizationConstraint>("Server", 0.5));
    line::opt::LineOptSolverOptions options;
    options.set_seed(42).set_max_iterations(60);
    const auto result = line::opt::LineOptSolver(problem, options).solve();
    CHECK(line::opt::scalar_value(result.variable_values.at("Server_servers")) == 6);
    CHECK(result.objective_value == doctest::Approx(60));
    CHECK(result.feasible);
}

TEST_CASE("line-opt projected gradient reaches the MATLAB M/M/1 rate optimum") {
    line::qn::Network<double> model("MM1");
    auto src = model.add_source("Arrivals");
    auto queue = model.add_queue("Server", line::lang::SchedStrategy::FCFS);
    auto sink = model.add_sink("Departures");
    auto jobs = model.add_open_class("Jobs");
    model.set_arrival(src, jobs, line::lang::Distrib<double>::exp_rate(3));
    model.set_service(queue, jobs, line::lang::Distrib<double>::exp_rate(4));
    line::qn::RoutingMatrix<double> routing;
    routing.set(jobs, jobs, src, queue, 1);
    routing.set(jobs, jobs, queue, sink, 1);
    model.link(routing);

    line::opt::OptimizationProblem problem(model);
    problem.add_variable(
               std::make_shared<line::opt::ServiceRate>("Server", "Jobs", 3.5, 8.0))
        .set_objective(std::make_shared<line::opt::MinimizeCost>(
            std::map<std::string, double>(),
            std::map<std::string, double>{{"Server", 20}}))
        .add_constraint(
            std::make_shared<line::opt::ResponseTimeConstraint>("Server", "Jobs", 0.5));
    line::opt::LineOptSolverOptions options;
    options.set_seed(42).set_optimizer("gradient").set_max_iterations(60)
        .set_gradient_restarts(2);
    const auto result = line::opt::LineOptSolver(problem, options).solve();
    CHECK(line::opt::scalar_value(result.variable_values.at("Server_Jobs_rate")) ==
          doctest::Approx(5.0).epsilon(0.02));
    CHECK(result.feasible);
}

TEST_CASE("line-opt Pareto sweep reproduces the MATLAB cost frontier") {
    line::qn::Network<double> model("MMc");
    auto src = model.add_source("Arrivals");
    auto queue = model.add_queue("Server", line::lang::SchedStrategy::FCFS);
    auto sink = model.add_sink("Departures");
    auto jobs = model.add_open_class("Jobs");
    model.set_arrival(src, jobs, line::lang::Distrib<double>::exp_rate(3));
    model.set_service(queue, jobs, line::lang::Distrib<double>::exp_rate(1));
    line::qn::RoutingMatrix<double> routing;
    routing.set(jobs, jobs, src, queue, 1);
    routing.set(jobs, jobs, queue, sink, 1);
    model.link(routing);

    line::opt::OptimizationProblem problem(model);
    problem.add_variable(std::make_shared<line::opt::ServerAllocation>("Server", 1, 10))
        .set_objective(std::make_shared<line::opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10}}));
    line::opt::ParetoSweep sweep(
        problem,
        [](double epsilon) {
            return std::make_shared<line::opt::UtilizationConstraint>("Server", epsilon);
        },
        {0.3, 0.4, 0.5, 0.6, 0.75, 0.9}, "bisection");
    sweep.solve();
    const auto frontier = sweep.frontier();
    const std::array<std::array<double, 2>, 5> expected =
        {{{0.3, 100}, {0.4, 80}, {0.5, 60}, {0.6, 50}, {0.75, 40}}};
    REQUIRE(frontier.size() == expected.size());
    for (std::size_t i = 0; i < frontier.size(); ++i) {
        CHECK(frontier[i].epsilon == expected[i][0]);
        CHECK(frontier[i].objective_value == expected[i][1]);
    }
}

TEST_CASE("line-opt class mapping selects one parallel station for the mapped class") {
    line::qn::Network<double> model("mapping");
    const auto source = model.add_source("Source");
    const auto first = model.add_queue("First");
    const auto second = model.add_queue("Second");
    const auto sink = model.add_sink("Sink");
    const auto mapped = model.add_open_class("Mapped");
    const auto other = model.add_open_class("Other");
    model.set_arrival(source, mapped, line::lang::Distrib<double>::exp_rate(1));
    model.set_arrival(source, other, line::lang::Distrib<double>::exp_rate(1));
    for (const auto jobclass : {mapped, other}) {
        model.set_service(first, jobclass, line::lang::Distrib<double>::exp_rate(3));
        model.set_service(second, jobclass, line::lang::Distrib<double>::exp_rate(3));
    }
    line::qn::RoutingMatrix<double> routing;
    for (const auto jobclass : {mapped, other}) {
        routing.set(jobclass, jobclass, source, first, 0.5);
        routing.set(jobclass, jobclass, source, second, 0.5);
        routing.set(jobclass, jobclass, first, sink, 1.0);
        routing.set(jobclass, jobclass, second, sink, 1.0);
    }
    model.link(routing);

    line::opt::ClassServiceMapping variable("Mapped", {"First", "Second"});
    CHECK(line::opt::scalar_value(variable.decode({0.0})) == 0.0);
    CHECK(line::opt::scalar_value(variable.decode({0.5})) == 1.0);
    CHECK(line::opt::scalar_value(variable.decode({1.0})) == 1.0);
    variable.apply(model, {1.0});
    const auto& structure = model.get_struct();
    CHECK(structure.get_route(mapped, mapped, source, first) == 0.0);
    CHECK(structure.get_route(mapped, mapped, source, second) == 1.0);
    CHECK(structure.get_route(other, other, source, first) == doctest::Approx(0.5));
    CHECK(structure.get_route(other, other, source, second) == doctest::Approx(0.5));
}

TEST_CASE("line-opt decomposition groups variables in MATLAB order") {
    line::qn::Network<double> model("decomposition");
    const auto source = model.add_source("Source");
    const auto queue = model.add_queue("Server");
    const auto sink = model.add_sink("Sink");
    const auto jobs = model.add_open_class("Jobs");
    model.set_arrival(source, jobs, line::lang::Distrib<double>::exp_rate(1));
    model.set_service(queue, jobs, line::lang::Distrib<double>::exp_rate(3));
    line::qn::RoutingMatrix<double> routing;
    routing.set(jobs, jobs, source, queue, 1.0);
    routing.set(jobs, jobs, queue, sink, 1.0);
    model.link(routing);

    const auto rate = std::make_shared<line::opt::ServiceRate>(
        "Server", "Jobs", 2.0, 5.0);
    const auto servers = std::make_shared<line::opt::ServerAllocation>(
        "Server", 1, 4);
    line::opt::OptimizationProblem problem(model);
    problem.add_variable(rate).add_variable(servers).set_objective(
        std::make_shared<line::opt::MinimizeCost>(
            std::map<std::string, double>{{"Server", 10.0}}));
    line::opt::DecompositionWorkflow workflow(problem);
    workflow.auto_decompose();
    REQUIRE(workflow.subproblems().size() == 2);
    CHECK(workflow.subproblems()[0].variable_type == "server_allocation");
    CHECK(workflow.subproblems()[1].variable_type == "service_rate");

    workflow.set_dependency("service_rate", "server_allocation");
    const auto ordered = workflow.execution_order();
    REQUIRE(ordered.size() == 2);
    CHECK(ordered[0].name == "service_rate");
    CHECK(ordered[1].name == "server_allocation");
}

TEST_CASE("line-opt open sensitivities match the M/M/1 closed form") {
    line::qn::Network<double> model("sensitivity");
    const auto source = model.add_source("Source");
    const auto queue = model.add_queue("Server");
    const auto sink = model.add_sink("Sink");
    const auto jobs = model.add_open_class("Jobs");
    model.set_arrival(source, jobs, line::lang::Distrib<double>::exp_rate(3.0));
    model.set_service(queue, jobs, line::lang::Distrib<double>::exp_rate(4.0));
    line::qn::RoutingMatrix<double> routing;
    routing.set(jobs, jobs, source, queue, 1.0);
    routing.set(jobs, jobs, queue, sink, 1.0);
    model.link(routing);

    const auto sensitivity = line::opt::compute_model_sensitivities(model);
    REQUIRE(sensitivity.has_value());
    line::opt::LineEvaluator evaluator(model, {});
    const auto evaluation = evaluator.evaluate({});
    REQUIRE(evaluation.sensitivities != nullptr);
    const std::string metric = "Server||Jobs";
    const std::string parameter = "rate||Server||Jobs";
    REQUIRE(sensitivity->for_kind("RespT") != nullptr);
    CHECK(sensitivity->for_kind("RespT")->at(metric).at(parameter) ==
          doctest::Approx(-1.0));
    CHECK(sensitivity->for_kind("QLen")->at(metric).at(parameter) ==
          doctest::Approx(-3.0));
    CHECK(sensitivity->for_kind("Tput")->at(metric).at(parameter) == 0.0);
    CHECK(sensitivity->for_kind("Util")->at("Server").at(parameter) ==
          doctest::Approx(-0.1875));
}

TEST_CASE("line-opt layered variables mutate a copied LQN model") {
    line::lqn::LqnBuilder<double> builder;
    builder.processor("ClientCpu", INFINITY, line::lang::SchedStrategy::INF);
    builder.processor("ServerCpu", 1, line::lang::SchedStrategy::PS);
    builder.task("Client", 3, line::lang::SchedStrategy::REF, "ClientCpu");
    builder.think_time("Client", line::lang::Distrib<double>::exp_mean(1.0));
    builder.task("Server", 1, line::lang::SchedStrategy::FCFS, "ServerCpu");
    builder.entry("ClientEntry", "Client");
    builder.entry("ServerEntry", "Server");
    builder.activity("ClientAct", line::lang::Distrib<double>::exp_mean(0.1), "Client");
    builder.bound_to("ClientAct", "ClientEntry");
    builder.sync_call("ClientAct", "ServerEntry", 1.0);
    builder.activity("ServerAct", line::lang::Distrib<double>::exp_mean(0.5), "Server");
    builder.bound_to("ServerAct", "ServerEntry");
    builder.replies_to("ServerAct", "ServerEntry");
    line::lqn::LqnModel<double> model = builder.model();

    line::opt::HostDemand demand("ServerAct", 0.2, 0.8);
    line::opt::TaskMultiplicity multiplicity("Server", 1, 5);
    line::opt::TaskReplication replication("Server", 1, 4);
    line::opt::ProcessorMultiplicity processors("ServerCpu", 1, 8);
    line::opt::TaskThinkTime task_think("Client", 0.5, 2.0);
    line::opt::ActivityThinkTime activity_think("ServerAct", 0.0, 0.4);
    demand.apply(model, demand.decode({0.5}));
    multiplicity.apply(model, {4.0});
    replication.apply(model, {3.0});
    processors.apply(model, {6.0});
    task_think.apply(model, {1.5});
    activity_think.apply(model, {0.25});
    CHECK(model.acts[1].hostdem.mean == doctest::Approx(0.5));
    CHECK(model.tasks[1].mult == 4.0);
    CHECK(model.tasks[1].repl == 3.0);
    CHECK(model.procs[1].mult == 6.0);
    CHECK(model.tasks[0].thinktime.mean == doctest::Approx(1.5));
    CHECK(model.acts[1].thinktime.mean == doctest::Approx(0.25));
    CHECK(demand.layers(model) == std::vector<std::string>{"ServerCpu"});
    CHECK(multiplicity.layers(model) ==
          std::vector<std::string>{"Server", "ServerCpu"});

    line::opt::LineEvaluator evaluator(
        model, {std::make_shared<line::opt::HostDemand>(demand)});
    const auto result = evaluator.evaluate({{"ServerAct_hostdemand", {0.4}}});
    CHECK(result.feasible);
    CHECK(result.solver_used == "SolverLN");
    CHECK(std::isfinite(result.system_response_time("Client")));
    const auto layered_sensitivity = evaluator.evaluate_layered_sensitivities(
        {{"ServerAct_hostdemand", {0.4}}});
    REQUIRE(layered_sensitivity.has_value());
    CHECK(layered_sensitivity->find("ServerCpu||Server") != layered_sensitivity->end());
}

TEST_CASE("line-opt freezes layered variables at their model values") {
    line::lqn::LqnBuilder<double> builder;
    builder.processor("Cpu", 1, line::lang::SchedStrategy::PS);
    builder.task("Client", 2, line::lang::SchedStrategy::REF, "Cpu");
    builder.think_time("Client", line::lang::Distrib<double>::exp_mean(1.0));
    builder.entry("Entry", "Client");
    builder.activity("Activity", line::lang::Distrib<double>::exp_mean(0.25), "Client");
    builder.bound_to("Activity", "Entry");
    line::opt::OptimizationProblem problem(builder.model());
    problem.add_variable(std::make_shared<line::opt::HostDemand>(
               "Activity", 0.1, 0.5))
        .set_objective(std::make_shared<line::opt::MinimizeSystemResponseTime>("Client"));
    line::opt::LineOptSolverOptions options;
    options.set_frozen_layers({"Cpu"});
    const auto result = line::opt::LineOptSolver(problem, options).solve();
    CHECK(result.terminated_by == "empty");
    CHECK(result.variable_values.empty());
}
