#ifndef LINE_OPT_PROBLEM_H
#define LINE_OPT_PROBLEM_H

#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "line/opt/evaluator.h"
#include "line/opt/objectives.h"

namespace line {
namespace opt {

using OptimizationModel = LineEvaluator::Model;

struct Scenario {
    OptimizationModel model;
    double weight = 1.0;
};

class OptimizationProblem {
public:
    explicit OptimizationProblem(qn::Network<double> model) : model_(std::move(model)) {}
    explicit OptimizationProblem(lqn::LqnModel<double> model) : model_(std::move(model)) {}
    explicit OptimizationProblem(OptimizationModel model) : model_(std::move(model)) {}

    OptimizationProblem& add_variable(VariablePtr variable) {
        variables_.push_back(std::move(variable));
        return *this;
    }
    OptimizationProblem& set_objective(ObjectivePtr objective) {
        objective_ = std::move(objective);
        return *this;
    }
    OptimizationProblem& add_constraint(ConstraintPtr constraint) {
        constraints_.push_back(std::move(constraint));
        return *this;
    }
    OptimizationProblem& set_fixed_variables(
        std::vector<std::pair<VariablePtr, Value>> variables) {
        fixed_ = std::move(variables);
        return *this;
    }
    OptimizationProblem& add_scenario(qn::Network<double> model, double weight = 1.0) {
        scenarios_.push_back({std::move(model), weight});
        return *this;
    }
    OptimizationProblem& add_scenario(lqn::LqnModel<double> model, double weight = 1.0) {
        scenarios_.push_back({std::move(model), weight});
        return *this;
    }
    OptimizationProblem& add_scenario(OptimizationModel model, double weight = 1.0) {
        scenarios_.push_back({std::move(model), weight});
        return *this;
    }

    bool is_layered() const {
        return std::holds_alternative<lqn::LqnModel<double>>(model_);
    }
    std::vector<std::string> validate() const {
        std::vector<std::string> errors;
        if (variables_.empty()) errors.push_back("No decision variables defined");
        if (!objective_) errors.push_back("Objective function is not set");
        for (const VariablePtr& variable : variables_) {
            const std::string type = variable->type();
            const bool layered_variable =
                type == "host_demand" || type == "think_time" ||
                type == "task_multiplicity" || type == "task_replication" ||
                type == "processor_multiplicity";
            if (is_layered() && !layered_variable)
                errors.push_back("Variable '" + variable->name() + "' (" + type +
                                 ") is a flat-network variable but the model is a LayeredNetwork");
            else if (!is_layered() && layered_variable)
                errors.push_back("Variable '" + variable->name() + "' is a LayeredNetwork "
                                 "variable but the model is a flat Network");
        }
        for (const Scenario& scenario : scenarios_)
            if (scenario.model.index() != model_.index())
                errors.push_back("Scenario model kind differs from the base model");
        return errors;
    }

    const qn::Network<double>& model() const {
        return std::get<qn::Network<double>>(model_);
    }
    const lqn::LqnModel<double>& layered_model() const {
        return std::get<lqn::LqnModel<double>>(model_);
    }
    const OptimizationModel& model_variant() const { return model_; }
    const std::vector<VariablePtr>& variables() const { return variables_; }
    const ObjectivePtr& objective() const { return objective_; }
    const std::vector<ConstraintPtr>& constraints() const { return constraints_; }
    const std::vector<std::pair<VariablePtr, Value>>& fixed_variables() const { return fixed_; }
    const std::vector<Scenario>& scenarios() const { return scenarios_; }

private:
    OptimizationModel model_;
    std::vector<VariablePtr> variables_;
    ObjectivePtr objective_;
    std::vector<ConstraintPtr> constraints_;
    std::vector<std::pair<VariablePtr, Value>> fixed_;
    std::vector<Scenario> scenarios_;
};

}  // namespace opt
}  // namespace line

#endif
