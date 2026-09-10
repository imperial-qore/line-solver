/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_LINE_OPT_SOLVER_H
#define LINE_OPT_LINE_OPT_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "line/opt/de/differential_evolution.h"
#include "line/opt/problem.h"
#include "line/util/error.h"

namespace line {
namespace opt {

struct LineOptSolverOptions {
    std::string strategy = "best1bin";
    std::size_t popsize = 15;
    double mutation_low = 0.5;
    double mutation_high = 1.0;
    double recombination = 0.7;
    double tolerance = 0.01;
    std::size_t max_iterations = 100;
    double time_limit = 300.0;
    std::optional<std::uint64_t> seed;
    bool verbose = false;
    double penalty_weight = 1e6;
    std::string scenario_aggregation = "worst";
    std::string optimizer = "evolution";
    double fd_step = 1e-6;
    double fd_step_layered = 1e-4;
    std::size_t gradient_restarts = 4;
    std::string lqn_gradient = "fd";
    std::size_t fd_refresh = 5;
    std::vector<std::string> frozen_layers;

    LineOptSolverOptions& set_strategy(std::string value) {
        strategy = std::move(value);
        return *this;
    }
    LineOptSolverOptions& set_popsize(std::size_t value) {
        popsize = value;
        return *this;
    }
    LineOptSolverOptions& set_mutation(double low, double high) {
        mutation_low = low;
        mutation_high = high;
        return *this;
    }
    LineOptSolverOptions& set_recombination(double value) {
        recombination = value;
        return *this;
    }
    LineOptSolverOptions& set_tolerance(double value) {
        tolerance = value;
        return *this;
    }
    LineOptSolverOptions& set_max_iterations(std::size_t value) {
        max_iterations = value;
        return *this;
    }
    LineOptSolverOptions& set_time_limit(double value) {
        time_limit = value;
        return *this;
    }
    LineOptSolverOptions& set_seed(std::uint64_t value) {
        seed = value;
        return *this;
    }
    LineOptSolverOptions& clear_seed() {
        seed.reset();
        return *this;
    }
    LineOptSolverOptions& set_verbose(bool value) {
        verbose = value;
        return *this;
    }
    LineOptSolverOptions& set_penalty_weight(double value) {
        penalty_weight = value;
        return *this;
    }
    LineOptSolverOptions& set_scenario_aggregation(std::string value) {
        scenario_aggregation = std::move(value);
        return *this;
    }
    LineOptSolverOptions& set_optimizer(std::string value) {
        optimizer = std::move(value);
        return *this;
    }
    LineOptSolverOptions& set_fd_step(double value) {
        fd_step = value;
        return *this;
    }
    LineOptSolverOptions& set_fd_step_layered(double value) {
        fd_step_layered = value;
        return *this;
    }
    LineOptSolverOptions& set_gradient_restarts(std::size_t value) {
        gradient_restarts = value;
        return *this;
    }
    LineOptSolverOptions& set_lqn_gradient(std::string value) {
        if (value != "fd" && value != "partial_sens" && value != "partial_plus_fd")
            throw InputError("lqn_gradient must be 'fd', 'partial_sens', or 'partial_plus_fd'");
        lqn_gradient = std::move(value);
        return *this;
    }
    LineOptSolverOptions& set_fd_refresh(std::size_t value) {
        fd_refresh = value;
        return *this;
    }
    LineOptSolverOptions& set_frozen_layers(std::vector<std::string> value) {
        frozen_layers = std::move(value);
        return *this;
    }
};

class LineOptSolver {
public:
    LineOptSolver(const OptimizationProblem& problem, LineOptSolverOptions options = {},
                  LineEvaluator::SolveFunction solve = {})
        : problem_(problem), options_(std::move(options)), solve_(std::move(solve)) {
        const std::vector<std::string> errors = problem_.validate();
        if (!errors.empty()) throw InputError("LineOptSolver: invalid optimization problem");
        free_variables_ = problem_.variables();
        fixed_variables_ = problem_.fixed_variables();
        freeze_layers();
        for (const auto& item : fixed_variables_)
            fixed_values_[item.first->name()] = item.second;

        evaluators_.emplace_back(problem_.model_variant(), free_variables_, fixed_variables_, solve_);
        scenario_weights_.push_back(1.0);
        for (const Scenario& scenario : problem_.scenarios()) {
            evaluators_.emplace_back(scenario.model, free_variables_, fixed_variables_, solve_);
            scenario_weights_.push_back(scenario.weight);
        }
        caches_.resize(evaluators_.size());
    }

    OptimizationResult solve() {
        start_ = Clock::now();
        iterations_ = 0;
        best_value_ = std::numeric_limits<double>::infinity();
        best_x_.clear();
        convergence_history_.clear();
        gradient_calls_ = 0;
        lqn_sensitivity_cache_.clear();
        resolved_seed_ = options_.seed.value_or(static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now().time_since_epoch()).count()) & UINT64_C(0x7fffffff));
        for (auto& cache : caches_) cache.clear();

        const std::vector<std::pair<double, double>> bounds = evaluators_[0].bounds();
        if (bounds.empty()) return empty_result();
        if (should_use_gradient()) return solve_gradient(bounds.size());
        return solve_evolution(bounds.size());
    }

private:
    using Clock = std::chrono::steady_clock;
    struct TimeLimit {};

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }

    VariableValues merge_values(const VariableValues& values) const {
        VariableValues out = fixed_values_;
        for (const auto& item : values) out[item.first] = item.second;
        return out;
    }

    const EvaluationResult& evaluate(std::size_t scenario, const VariableValues& values) {
        auto found = caches_[scenario].find(values);
        if (found != caches_[scenario].end()) return found->second;
        EvaluationResult result = evaluators_[scenario].evaluate(values);
        return caches_[scenario].emplace(values, std::move(result)).first->second;
    }

    double aggregate_scenarios(const std::vector<double>& values) const {
        if (values.size() == 1) return values[0];
        if (options_.scenario_aggregation == "mean") {
            double weighted = 0.0, weights = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i) {
                weighted += scenario_weights_[i] * values[i];
                weights += scenario_weights_[i];
            }
            return weighted / weights;
        }
        if (options_.scenario_aggregation != "worst")
            throw InputError("LineOptSolver: unknown scenario aggregation '" +
                             options_.scenario_aggregation + "'");
        return *std::max_element(values.begin(), values.end());
    }

    double objective_function(const std::vector<double>& x) {
        if (elapsed() >= options_.time_limit) throw TimeLimit();
        const VariableValues values = evaluators_[0].decode(x);
        const VariableValues all_values = merge_values(values);
        std::vector<double> scenario_values(evaluators_.size(), 0.0);
        for (std::size_t i = 0; i < evaluators_.size(); ++i) {
            const EvaluationResult& result = evaluate(i, values);
            if (!result.feasible) return std::numeric_limits<double>::infinity();
            double value = problem_.objective()->evaluate_with_penalty(
                result, all_values, options_.penalty_weight);
            for (const ConstraintPtr& constraint : problem_.constraints())
                value += constraint->evaluate(result, all_values) * options_.penalty_weight;
            scenario_values[i] = value;
        }
        const double total = aggregate_scenarios(scenario_values);
        if (total < best_value_) {
            best_value_ = total;
            best_x_ = x;
        }
        return total;
    }

    OptimizationResult solve_evolution(std::size_t dimension) {
        std::vector<double> low(dimension, 0.0), high(dimension, 1.0);
        de::DifferentialEvolution engine(
            [this](const std::vector<double>& x) { return objective_function(x); }, low, high,
            options_.strategy, options_.popsize, options_.max_iterations,
            options_.mutation_low, options_.mutation_high, options_.recombination,
            options_.tolerance, resolved_seed_);
        engine.set_callback([this](const std::vector<double>&, std::size_t nit) {
            iterations_ = nit;
            convergence_history_.push_back(best_value_);
            return elapsed() >= options_.time_limit;
        });

        bool timed_out = false;
        std::vector<double> x;
        double value = std::numeric_limits<double>::infinity();
        try {
            const de::DifferentialEvolutionResult result = engine.solve();
            if (!best_x_.empty() && best_value_ <= result.fun) {
                x = best_x_;
                value = best_value_;
            } else {
                x = result.x;
                value = result.fun;
            }
        } catch (const TimeLimit&) {
            timed_out = true;
            if (best_x_.empty()) return empty_result();
            x = best_x_;
            value = best_value_;
        }
        OptimizationResult result = build_result(x, value);
        if (timed_out) result.terminated_by = "time_limit";
        return result;
    }

    bool all_continuous() const {
        if (free_variables_.empty()) return false;
        for (const VariablePtr& variable : free_variables_) {
            const std::string type = variable->type();
            if (type != "service_rate" && type != "routing" && type != "host_demand" &&
                type != "think_time")
                return false;
        }
        return true;
    }

    bool should_use_gradient() const {
        if (options_.optimizer == "gradient") return true;
        if (options_.optimizer == "evolution") return false;
        if (options_.optimizer == "auto") return all_continuous();
        throw InputError("LineOptSolver: unknown optimizer '" + options_.optimizer + "'");
    }

    std::vector<double> finite_difference_gradient(const std::vector<double>& x) {
        const double h = evaluators_[0].is_layered() ? options_.fd_step_layered
                                                      : options_.fd_step;
        std::vector<double> gradient(x.size(), 0.0);
        std::optional<double> centre;
        for (std::size_t i = 0; i < x.size(); ++i) {
            std::vector<double> plus = x, minus = x;
            plus[i] = std::min(1.0, x[i] + h);
            minus[i] = std::max(0.0, x[i] - h);
            const double fp = objective_function(plus);
            const double fm = objective_function(minus);
            if (std::isfinite(fp) && std::isfinite(fm) && plus[i] > minus[i]) {
                gradient[i] = (fp - fm) / (plus[i] - minus[i]);
                continue;
            }
            if (!centre) centre = objective_function(x);
            if (std::isfinite(fp) && std::isfinite(*centre) && plus[i] > x[i])
                gradient[i] = (fp - *centre) / (plus[i] - x[i]);
            else if (std::isfinite(fm) && std::isfinite(*centre) && x[i] > minus[i])
                gradient[i] = (*centre - fm) / (x[i] - minus[i]);
        }
        return gradient;
    }

    std::vector<double> projected_gradient_descent(std::vector<double> x) {
        for (double& value : x) value = std::clamp(value, 0.0, 1.0);
        double f = objective_function(x);
        for (std::size_t iter = 0; iter < options_.max_iterations; ++iter) {
            const std::vector<double> gradient = objective_gradient(x);
            double norm2 = 0.0;
            for (double value : gradient) norm2 += value * value;
            if (std::sqrt(norm2) < 1e-9) break;
            double step = 1.0;
            bool improved = false;
            for (std::size_t search = 0; search < 30; ++search) {
                std::vector<double> next(x.size());
                for (std::size_t i = 0; i < x.size(); ++i)
                    next[i] = std::clamp(x[i] - step * gradient[i], 0.0, 1.0);
                const double fn = objective_function(next);
                if (std::isfinite(fn) && fn < f - 1e-12) {
                    x = std::move(next);
                    f = fn;
                    improved = true;
                    break;
                }
                step *= 0.5;
            }
            if (!improved) break;
        }
        return x;
    }

    OptimizationResult solve_gradient(std::size_t dimension) {
        const std::size_t starts = std::max<std::size_t>(1, options_.gradient_restarts);
        de::NumpyRandomState random(resolved_seed_);
        std::vector<double> best;
        double best_fun = std::numeric_limits<double>::infinity();
        bool timed_out = false;
        for (std::size_t restart = 0; restart < starts; ++restart) {
            if (elapsed() >= options_.time_limit) {
                timed_out = true;
                break;
            }
            std::vector<double> initial =
                restart == 0 ? std::vector<double>(dimension, 0.5)
                             : random.uniform(0.0, 1.0, dimension);
            try {
                std::vector<double> candidate = projected_gradient_descent(initial);
                const double value = objective_function(candidate);
                if (std::isfinite(value) && value < best_fun) {
                    best = std::move(candidate);
                    best_fun = value;
                }
            } catch (const TimeLimit&) {
                timed_out = true;
                break;
            }
        }
        if (best.empty() || (!best_x_.empty() && best_value_ < best_fun)) {
            if (!best_x_.empty()) {
                best = best_x_;
                best_fun = best_value_;
            } else {
                return empty_result();
            }
        }
        OptimizationResult result = build_result(best, best_fun);
        if (timed_out) result.terminated_by = "time_limit";
        return result;
    }

    void freeze_layers() {
        if (!problem_.is_layered() || options_.frozen_layers.empty()) return;
        const lqn::LqnModel<double>& model = problem_.layered_model();
        std::set<std::string> frozen(options_.frozen_layers.begin(), options_.frozen_layers.end());
        std::set<std::string> already;
        for (const auto& fixed : fixed_variables_) already.insert(fixed.first->name());
        std::vector<VariablePtr> kept;
        for (const VariablePtr& variable : free_variables_) {
            const std::vector<std::string> variable_layers = variable->layers(model);
            bool intersects = false;
            for (const std::string& layer : variable_layers)
                if (frozen.find(layer) != frozen.end()) intersects = true;
            if (!intersects) {
                kept.push_back(variable);
                continue;
            }
            if (already.find(variable->name()) != already.end()) continue;
            const std::optional<Value> value = variable->current_value(model);
            if (value) {
                fixed_variables_.push_back({variable, *value});
                already.insert(variable->name());
            }
        }
        free_variables_ = std::move(kept);
    }

    std::vector<double> objective_gradient(const std::vector<double>& x) {
        if (!evaluators_[0].is_layered()) return finite_difference_gradient(x);
        const std::string& mode = options_.lqn_gradient;
        if (mode == "partial_sens" || mode == "partial_plus_fd") {
            ++gradient_calls_;
            const std::size_t refresh = std::max<std::size_t>(1, options_.fd_refresh);
            if (!(mode == "partial_plus_fd" && gradient_calls_ % refresh == 0)) {
                const std::optional<std::vector<double>> gradient =
                    lqn_analytic_gradient(x);
                if (gradient) return *gradient;
            }
        }
        return finite_difference_gradient(x);
    }

    double scalar_objective(const EvaluationResult& result,
                            const VariableValues& values) const {
        double out = problem_.objective()->evaluate_with_penalty(
            result, values, options_.penalty_weight);
        for (const ConstraintPtr& constraint : problem_.constraints())
            out += constraint->evaluate(result, values) * options_.penalty_weight;
        return out;
    }

    double scalar_metric_derivative(const EvaluationResult& base,
                                    const VariableValues& values,
                                    const std::string& kind,
                                    const std::string& key) const {
        const double h = options_.fd_step;
        EvaluationResult plus = base, minus = base;
        std::map<std::string, double>* plus_map = nullptr;
        std::map<std::string, double>* minus_map = nullptr;
        if (kind == "RespT") {
            plus_map = &plus.response_times;
            minus_map = &minus.response_times;
        } else if (kind == "QLen") {
            plus_map = &plus.queue_lengths;
            minus_map = &minus.queue_lengths;
        } else if (kind == "Tput") {
            plus_map = &plus.throughputs;
            minus_map = &minus.throughputs;
        } else if (kind == "Util") {
            plus_map = &plus.utilizations;
            minus_map = &minus.utilizations;
        } else {
            return 0.0;
        }
        const auto found = plus_map->find(key);
        if (found == plus_map->end()) return 0.0;
        (*plus_map)[key] = found->second + h;
        (*minus_map)[key] = found->second - h;
        return (scalar_objective(plus, values) - scalar_objective(minus, values)) /
               (2.0 * h);
    }

    std::optional<std::vector<double>> lqn_analytic_gradient(
        const std::vector<double>& x) {
        if (evaluators_.size() != 1 || free_variables_.empty()) return std::nullopt;
        for (const VariablePtr& variable : free_variables_)
            if (!variable->supports_sensitivity()) return std::nullopt;
        const VariableValues values = evaluators_[0].decode(x);
        const VariableValues all_values = merge_values(values);
        const EvaluationResult& result = evaluate(0, values);
        if (!result.feasible) return std::nullopt;

        auto cached = lqn_sensitivity_cache_.find(values);
        if (cached == lqn_sensitivity_cache_.end())
            cached = lqn_sensitivity_cache_.emplace(
                values, evaluators_[0].evaluate_layered_sensitivities(values)).first;
        if (!cached->second) return std::nullopt;
        const LineEvaluator::LayerSensitivityMap& sensitivities = *cached->second;
        const lqn::LqnModel<double>& model = problem_.layered_model();
        std::vector<double> gradient(x.size(), 0.0);
        std::size_t offset = 0;
        for (const VariablePtr& variable : free_variables_) {
            const auto row = sensitivities.find(variable->sensitivity_key(model));
            if (row != sensitivities.end()) {
                const std::map<std::string, std::string> targets =
                    variable->sensitivity_metric_targets(model);
                double scalar_rate = 0.0;
                const struct {
                    const char* kind;
                    double LayerSensitivity::*field;
                } metrics[] = {{"Tput", &LayerSensitivity::throughput},
                               {"RespT", &LayerSensitivity::response_time},
                               {"QLen", &LayerSensitivity::queue_length},
                               {"Util", &LayerSensitivity::utilization}};
                for (const auto& metric : metrics) {
                    const auto target = targets.find(metric.kind);
                    if (target == targets.end()) continue;
                    scalar_rate += scalar_metric_derivative(
                        result, all_values, metric.kind, target->second) *
                        row->second.*(metric.field);
                }
                const auto value = all_values.find(variable->name());
                if (value != all_values.end())
                    gradient[offset] = scalar_rate * variable->rate_jacobian(value->second) *
                                       variable->decode_jacobian(x[offset]);
            }
            offset += variable->dimension();
        }
        return gradient;
    }

    std::vector<ConstraintPtr> all_constraints() const {
        std::vector<ConstraintPtr> constraints = problem_.objective()->constraints;
        constraints.insert(constraints.end(), problem_.constraints().begin(),
                           problem_.constraints().end());
        return constraints;
    }

    OptimizationResult build_result(const std::vector<double>& x, double objective_value) {
        OptimizationResult out;
        out.objective_value = objective_value;
        out.variable_values = evaluators_[0].decode(x);
        out.iterations = iterations_;
        out.solve_time = elapsed();
        out.convergence_history = convergence_history_;
        for (const LineEvaluator& evaluator : evaluators_)
            out.model_evaluations += evaluator.evaluation_count();

        const VariableValues all_values = merge_values(out.variable_values);
        out.feasible = true;
        for (std::size_t i = 0; i < evaluators_.size(); ++i) {
            const EvaluationResult& result = evaluate(i, out.variable_values);
            if (!result.feasible) {
                out.feasible = false;
                continue;
            }
            for (const ConstraintPtr& constraint : all_constraints()) {
                const double violation = constraint->evaluate(result, all_values);
                if (violation > 0.0) {
                    out.feasible = false;
                    out.constraint_violations[constraint->name()] = std::max(
                        out.constraint_violations[constraint->name()], violation);
                }
            }
        }
        if (elapsed() >= options_.time_limit)
            out.terminated_by = "time_limit";
        else if (iterations_ >= options_.max_iterations)
            out.terminated_by = "iterations";
        else
            out.terminated_by = "convergence";
        return out;
    }

    OptimizationResult empty_result() const {
        OptimizationResult out;
        out.objective_value = 0.0;
        out.feasible = true;
        out.terminated_by = "empty";
        return out;
    }

    const OptimizationProblem& problem_;
    LineOptSolverOptions options_;
    LineEvaluator::SolveFunction solve_;
    std::vector<VariablePtr> free_variables_;
    std::vector<std::pair<VariablePtr, Value>> fixed_variables_;
    VariableValues fixed_values_;
    std::vector<LineEvaluator> evaluators_;
    std::vector<double> scenario_weights_;
    std::vector<std::map<VariableValues, EvaluationResult>> caches_;
    Clock::time_point start_;
    std::size_t iterations_ = 0;
    double best_value_ = std::numeric_limits<double>::infinity();
    std::vector<double> best_x_;
    std::vector<double> convergence_history_;
    std::size_t gradient_calls_ = 0;
    std::uint64_t resolved_seed_ = 0;
    std::map<VariableValues, std::optional<LineEvaluator::LayerSensitivityMap>>
        lqn_sensitivity_cache_;
};

}  // namespace opt
}  // namespace line

#endif
