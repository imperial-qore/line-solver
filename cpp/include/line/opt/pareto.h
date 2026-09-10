/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_PARETO_H
#define LINE_OPT_PARETO_H

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "line/opt/bisection_solver.h"
#include "line/opt/line_opt_solver.h"

namespace line {
namespace opt {

struct ParetoPoint {
    double epsilon = 0.0;
    double objective_value = std::numeric_limits<double>::infinity();
    bool feasible = false;
    OptimizationResult result;
};

class ParetoSweep {
public:
    using ConstraintFactory = std::function<ConstraintPtr(double)>;

    ParetoSweep(const OptimizationProblem& problem, ConstraintFactory constraint_factory,
                std::vector<double> epsilons, std::string solver = "de",
                LineEvaluator::SolveFunction solve = {})
        : problem_(problem),
          constraint_factory_(std::move(constraint_factory)),
          epsilons_(std::move(epsilons)),
          solver_(std::move(solver)),
          solve_(std::move(solve)) {
        if (!constraint_factory_)
            throw InputError("ParetoSweep: the constraint factory is empty");
        if (solver_ != "de" && solver_ != "bisection")
            throw InputError("ParetoSweep: solver must be 'de' or 'bisection'");
    }

    const std::vector<ParetoPoint>& solve(const LineOptSolverOptions& options = {}) {
        points_.clear();
        for (double epsilon : epsilons_) {
            OptimizationProblem clone = clone_problem(epsilon);
            OptimizationResult result =
                solver_ == "bisection" ? BisectionSolver(clone, "min_feasible", solve_).solve()
                                       : LineOptSolver(clone, options, solve_).solve();
            points_.push_back(
                {epsilon, result.objective_value, result.feasible, std::move(result)});
        }
        return points_;
    }

    const std::vector<ParetoPoint>& points() const { return points_; }

    std::vector<ParetoPoint> frontier() const {
        std::vector<ParetoPoint> out;
        for (const ParetoPoint& point : points_) {
            if (!point.feasible) continue;
            bool dominated = false;
            for (const ParetoPoint& other : points_) {
                if (!other.feasible) continue;
                if (other.objective_value <= point.objective_value &&
                    other.epsilon <= point.epsilon &&
                    (other.objective_value < point.objective_value ||
                     other.epsilon < point.epsilon)) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) out.push_back(point);
        }
        std::sort(out.begin(), out.end(),
                  [](const ParetoPoint& a, const ParetoPoint& b) {
                      return a.epsilon < b.epsilon;
                  });
        return out;
    }

private:
    OptimizationProblem clone_problem(double epsilon) const {
        OptimizationProblem clone(problem_.model_variant());
        for (const VariablePtr& variable : problem_.variables()) clone.add_variable(variable);
        clone.set_objective(problem_.objective());
        for (const ConstraintPtr& constraint : problem_.constraints())
            clone.add_constraint(constraint);
        clone.set_fixed_variables(problem_.fixed_variables());
        for (const Scenario& scenario : problem_.scenarios())
            clone.add_scenario(scenario.model, scenario.weight);
        clone.add_constraint(constraint_factory_(epsilon));
        return clone;
    }

    const OptimizationProblem& problem_;
    ConstraintFactory constraint_factory_;
    std::vector<double> epsilons_;
    std::string solver_;
    LineEvaluator::SolveFunction solve_;
    std::vector<ParetoPoint> points_;
};

}  // namespace opt
}  // namespace line

#endif
