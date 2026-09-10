/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_DECOMPOSITION_H
#define LINE_OPT_DECOMPOSITION_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "line/opt/line_opt_solver.h"

namespace line {
namespace opt {

struct SubProblem {
    std::string name;
    std::string variable_type;
    std::vector<VariablePtr> variables;
    VariableValues fixed_values;

    std::vector<std::string> variable_names() const {
        std::vector<std::string> names;
        names.reserve(variables.size());
        for (const VariablePtr& variable : variables) names.push_back(variable->name());
        return names;
    }
};

struct SubProblemResult {
    std::string name;
    OptimizationResult result;
    VariableValues variables_fixed;
};

struct WorkflowResult {
    double final_objective = std::numeric_limits<double>::infinity();
    std::map<std::string, SubProblemResult> subproblem_results;
    std::size_t cycles_completed = 0;
    bool converged = false;
    double total_solve_time = 0.0;
    std::vector<double> objective_history;
    VariableValues final_variable_values;
    std::vector<std::string> frozen_layers;
    std::size_t model_evaluations = 0;

    const SubProblemResult* subproblem_result(const std::string& name) const {
        const auto found = subproblem_results.find(name);
        return found == subproblem_results.end() ? nullptr : &found->second;
    }
    const Value* final_variable_value(const std::string& name) const {
        const auto found = final_variable_values.find(name);
        return found == final_variable_values.end() ? nullptr : &found->second;
    }
};

class DecompositionWorkflow {
public:
    explicit DecompositionWorkflow(const OptimizationProblem& problem) : problem_(problem) {}

    static const std::vector<std::string>& default_order() {
        static const std::vector<std::string> order = {
            "server_allocation", "station_replicas", "service_rate", "job_population",
            "class_priority", "routing", "class_mapping", "processor_multiplicity",
            "task_multiplicity", "task_replication", "host_demand", "think_time"};
        return order;
    }

    DecompositionWorkflow& set_solver_options(LineOptSolverOptions options) {
        solver_options_ = std::move(options);
        return *this;
    }

    DecompositionWorkflow& auto_decompose() {
        std::map<std::string, std::vector<VariablePtr>> by_type;
        std::vector<std::string> discovered;
        for (const VariablePtr& variable : problem_.variables()) {
            const std::string type = variable->type();
            if (by_type.find(type) == by_type.end()) discovered.push_back(type);
            by_type[type].push_back(variable);
        }
        subproblems_.clear();
        std::set<std::string> used;
        for (const std::string& type : default_order()) {
            const auto found = by_type.find(type);
            if (found == by_type.end()) continue;
            subproblems_.push_back({type, type, found->second, {}});
            used.insert(type);
        }
        for (const std::string& type : discovered)
            if (used.find(type) == used.end())
                subproblems_.push_back({type, type, by_type.at(type), {}});
        return *this;
    }

    DecompositionWorkflow& set_dependency(const std::string& from_problem,
                                          const std::string& to_problem) {
        dependencies_[to_problem].push_back(from_problem);
        return *this;
    }

    DecompositionWorkflow& add_subproblem(std::string name,
                                          std::vector<VariablePtr> variables,
                                          std::vector<std::string> after = {}) {
        const std::string type = variables.empty() ? "custom" : variables.front()->type();
        const std::string dependency_target = name;
        subproblems_.push_back({std::move(name), type, std::move(variables), {}});
        for (const std::string& dependency : after)
            set_dependency(dependency, dependency_target);
        return *this;
    }

    const std::vector<SubProblem>& subproblems() const { return subproblems_; }

    std::vector<SubProblem> execution_order() const {
        if (dependencies_.empty()) return subproblems_;
        std::map<std::string, std::size_t> indegree;
        std::map<std::string, std::vector<std::string>> adjacent;
        std::map<std::string, SubProblem> by_name;
        for (const SubProblem& subproblem : subproblems_) {
            indegree[subproblem.name] = 0;
            adjacent[subproblem.name] = {};
            by_name[subproblem.name] = subproblem;
        }
        for (const auto& target : dependencies_) {
            if (indegree.find(target.first) == indegree.end()) continue;
            for (const std::string& source : target.second) {
                if (adjacent.find(source) == adjacent.end()) continue;
                adjacent[source].push_back(target.first);
                ++indegree[target.first];
            }
        }
        std::vector<std::string> queue;
        for (const SubProblem& subproblem : subproblems_)
            if (indegree[subproblem.name] == 0) queue.push_back(subproblem.name);
        std::vector<SubProblem> ordered;
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const std::string current = queue[head];
            ordered.push_back(by_name.at(current));
            for (const std::string& successor : adjacent[current])
                if (--indegree[successor] == 0) queue.push_back(successor);
        }
        return ordered.size() == subproblems_.size() ? ordered : subproblems_;
    }

    WorkflowResult solve_sequential(std::size_t max_cycles = 10, double tolerance = 0.01) const {
        const auto start = std::chrono::steady_clock::now();
        WorkflowResult out;
        if (subproblems_.empty()) {
            out.converged = true;
            return out;
        }
        const std::vector<SubProblem> ordered = execution_order();
        VariableValues fixed_values;
        double previous_objective = std::numeric_limits<double>::infinity();
        for (std::size_t cycle = 1; cycle <= max_cycles; ++cycle) {
            for (const SubProblem& subproblem : ordered) {
                OptimizationProblem partial = create_partial_problem(subproblem, fixed_values);
                OptimizationResult solved = LineOptSolver(partial, solver_options_).solve();
                SubProblemResult subresult{subproblem.name, solved, fixed_values};
                out.model_evaluations += solved.model_evaluations;
                out.subproblem_results[subproblem.name] = std::move(subresult);
                for (const auto& value : solved.variable_values)
                    fixed_values[value.first] = value.second;
            }
            const double current_objective = evaluate_full_objective(fixed_values);
            ++out.model_evaluations;
            out.objective_history.push_back(current_objective);
            if (std::abs(current_objective - previous_objective) < tolerance) {
                out.converged = true;
                break;
            }
            previous_objective = current_objective;
            out.cycles_completed = cycle;
        }
        if (!out.objective_history.empty()) out.final_objective = out.objective_history.back();
        out.final_variable_values = std::move(fixed_values);
        out.total_solve_time = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        return out;
    }

    WorkflowResult solve_hierarchical() const { return solve_sequential(1, 0.0); }

    WorkflowResult solve_layered(std::size_t max_cycles = 10, double tolerance = 0.01,
                                 bool auto_freeze = true, double freeze_tolerance = 1e-3,
                                 std::vector<std::string> frozen_layers = {}) const {
        if (!problem_.is_layered()) return solve_sequential(max_cycles, tolerance);
        const auto start = std::chrono::steady_clock::now();
        WorkflowResult out;
        const lqn::LqnModel<double>& model = problem_.layered_model();
        std::map<std::string, std::vector<VariablePtr>> groups;
        std::vector<std::string> group_order;
        for (const VariablePtr& variable : problem_.variables()) {
            const std::vector<std::string> owned = variable->layers(model);
            const std::string layer = owned.empty() ? "_nolayer" : owned.front();
            if (groups.find(layer) == groups.end()) group_order.push_back(layer);
            groups[layer].push_back(variable);
        }

        VariableValues fixed_values;
        std::map<std::string, std::vector<double>> previous_signatures;
        bool have_previous_signatures = false;
        double previous_objective = std::numeric_limits<double>::infinity();
        LineEvaluator full_evaluator(problem_.model_variant(), problem_.variables(),
                                     problem_.fixed_variables());
        for (std::size_t cycle = 1; cycle <= max_cycles; ++cycle) {
            for (const std::string& layer : group_order) {
                if (std::find(frozen_layers.begin(), frozen_layers.end(), layer) !=
                    frozen_layers.end())
                    continue;
                const SubProblem subproblem{layer, groups[layer].front()->type(),
                                            groups[layer], {}};
                OptimizationProblem partial = create_partial_problem(subproblem, fixed_values);
                OptimizationResult solved = LineOptSolver(partial, solver_options_).solve();
                out.model_evaluations += solved.model_evaluations;
                out.subproblem_results[layer] = {layer, solved, fixed_values};
                for (const auto& value : solved.variable_values)
                    fixed_values[value.first] = value.second;
            }

            const EvaluationResult evaluation = full_evaluator.evaluate(fixed_values);
            ++out.model_evaluations;
            double current_objective = std::numeric_limits<double>::infinity();
            std::map<std::string, std::vector<double>> signatures;
            if (evaluation.feasible) {
                VariableValues all_values;
                for (const auto& fixed : problem_.fixed_variables())
                    all_values[fixed.first->name()] = fixed.second;
                for (const auto& value : fixed_values) all_values[value.first] = value.second;
                current_objective = problem_.objective()->evaluate_with_penalty(
                    evaluation, all_values, solver_options_.penalty_weight);
                for (const ConstraintPtr& constraint : problem_.constraints())
                    current_objective += constraint->evaluate(evaluation, all_values) *
                                         solver_options_.penalty_weight;
                signatures = layer_signatures(evaluation, group_order);
            }

            if (auto_freeze && have_previous_signatures) {
                std::map<std::string, double> moved;
                bool active_moved = false;
                for (const std::string& layer : group_order) {
                    const auto before = previous_signatures.find(layer);
                    const auto after = signatures.find(layer);
                    moved[layer] = signature_delta(
                        before == previous_signatures.end() ? std::vector<double>()
                                                            : before->second,
                        after == signatures.end() ? std::vector<double>() : after->second);
                    if (std::find(frozen_layers.begin(), frozen_layers.end(), layer) ==
                            frozen_layers.end() &&
                        moved[layer] > freeze_tolerance)
                        active_moved = true;
                }
                for (const std::string& layer : group_order) {
                    const auto frozen = std::find(frozen_layers.begin(), frozen_layers.end(), layer);
                    if (frozen != frozen_layers.end()) {
                        if (active_moved) frozen_layers.erase(frozen);
                    } else if (moved[layer] < freeze_tolerance) {
                        frozen_layers.push_back(layer);
                    }
                }
            }

            previous_signatures = std::move(signatures);
            have_previous_signatures = true;
            out.objective_history.push_back(current_objective);
            if (std::abs(current_objective - previous_objective) < tolerance) {
                out.converged = true;
                break;
            }
            previous_objective = current_objective;
            out.cycles_completed = cycle;
            if (frozen_layers.size() >= group_order.size()) {
                out.converged = true;
                break;
            }
        }
        if (!out.objective_history.empty()) out.final_objective = out.objective_history.back();
        out.final_variable_values = std::move(fixed_values);
        std::sort(frozen_layers.begin(), frozen_layers.end());
        frozen_layers.erase(std::unique(frozen_layers.begin(), frozen_layers.end()),
                            frozen_layers.end());
        out.frozen_layers = std::move(frozen_layers);
        out.total_solve_time = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        return out;
    }

private:
    static std::map<std::string, std::vector<double>> layer_signatures(
        const EvaluationResult& result, const std::vector<std::string>& layers) {
        std::map<std::string, std::vector<double>> out;
        for (const std::string& layer : layers) {
            const double utilization = result.utilization(layer);
            const double queue_length = result.queue_length(layer);
            const double throughput = result.throughput(layer);
            const double response_time = result.response_time(layer);
            const double probes[] = {utilization, queue_length, throughput};
            bool present = false;
            for (double probe : probes)
                if (probe != 0.0 && !std::isinf(probe)) present = true;
            if (present)
                out[layer] = {utilization, queue_length, throughput, response_time};
        }
        return out;
    }

    static double signature_delta(const std::vector<double>& before,
                                  const std::vector<double>& after) {
        if (before.empty() || after.empty()) return std::numeric_limits<double>::infinity();
        double delta = 0.0;
        const std::size_t size = std::min(before.size(), after.size());
        for (std::size_t i = 0; i < size; ++i) {
            if (!std::isfinite(before[i]) || !std::isfinite(after[i])) continue;
            delta = std::max(delta, std::abs(after[i] - before[i]) /
                                    (std::abs(before[i]) + 1e-12));
        }
        return delta;
    }

    OptimizationProblem create_partial_problem(const SubProblem& subproblem,
                                               const VariableValues& fixed_values) const {
        OptimizationProblem partial(problem_.model_variant());
        for (const VariablePtr& variable : subproblem.variables) partial.add_variable(variable);
        std::set<std::string> subproblem_names;
        for (const VariablePtr& variable : subproblem.variables)
            subproblem_names.insert(variable->name());
        std::vector<std::pair<VariablePtr, Value>> fixed = problem_.fixed_variables();
        for (const VariablePtr& variable : problem_.variables()) {
            if (subproblem_names.find(variable->name()) != subproblem_names.end()) continue;
            const auto found = fixed_values.find(variable->name());
            if (found != fixed_values.end()) fixed.push_back({variable, found->second});
        }
        partial.set_fixed_variables(std::move(fixed));
        partial.set_objective(problem_.objective());
        for (const ConstraintPtr& constraint : problem_.constraints())
            partial.add_constraint(constraint);
        for (const Scenario& scenario : problem_.scenarios())
            partial.add_scenario(scenario.model, scenario.weight);
        return partial;
    }

    double evaluate_full_objective(const VariableValues& values) const {
        LineEvaluator evaluator(problem_.model_variant(), problem_.variables(),
                                problem_.fixed_variables());
        const EvaluationResult result = evaluator.evaluate(values);
        if (!result.feasible) return std::numeric_limits<double>::infinity();
        VariableValues all_values;
        for (const auto& fixed : problem_.fixed_variables())
            all_values[fixed.first->name()] = fixed.second;
        for (const auto& value : values) all_values[value.first] = value.second;
        double objective = problem_.objective()->evaluate_with_penalty(
            result, all_values, solver_options_.penalty_weight);
        for (const ConstraintPtr& constraint : problem_.constraints())
            objective += constraint->evaluate(result, all_values) * solver_options_.penalty_weight;
        return objective;
    }

    const OptimizationProblem& problem_;
    std::vector<SubProblem> subproblems_;
    std::map<std::string, std::vector<std::string>> dependencies_;
    LineOptSolverOptions solver_options_;
};

}  // namespace opt
}  // namespace line

#endif
