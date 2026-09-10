#ifndef LINE_OPT_EVALUATOR_H
#define LINE_OPT_EVALUATOR_H

#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "line/opt/results.h"
#include "line/opt/sensitivity.h"
#include "line/opt/variables.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/solver_chain_tables.h"

namespace line {
namespace opt {

struct LayerSensitivity {
    double throughput = 0.0;
    double response_time = 0.0;
    double queue_length = 0.0;
    double utilization = 0.0;
};

class LineEvaluator {
public:
    using Model = std::variant<qn::Network<double>, lqn::LqnModel<double>>;
    using SolveFunction = std::function<mva::AvgResult<double>(
        const qn::NetworkStruct<double>&)>;
    using LayerSensitivityMap = std::map<std::string, LayerSensitivity>;

    LineEvaluator(qn::Network<double> model, std::vector<VariablePtr> variables,
                  std::vector<std::pair<VariablePtr, Value>> fixed = {},
                  SolveFunction solve = {})
        : LineEvaluator(Model(std::move(model)), std::move(variables), std::move(fixed),
                        std::move(solve)) {}

    LineEvaluator(lqn::LqnModel<double> model, std::vector<VariablePtr> variables,
                  std::vector<std::pair<VariablePtr, Value>> fixed = {})
        : LineEvaluator(Model(std::move(model)), std::move(variables), std::move(fixed), {}) {}

    LineEvaluator(Model model, std::vector<VariablePtr> variables,
                  std::vector<std::pair<VariablePtr, Value>> fixed = {},
                  SolveFunction solve = {})
        : base_(std::move(model)), variables_(std::move(variables)),
          fixed_(std::move(fixed)), solve_(std::move(solve)) {
        if (!solve_)
            solve_ = [](const qn::NetworkStruct<double>& structure) {
                return mva::solver_mva_run_analyzer(structure, mva::MvaOptions(), Matrix<double>());
            };
        std::size_t dimension = 0;
        for (const VariablePtr& variable : variables_) {
            offsets_.push_back(dimension);
            dimension += variable->dimension();
        }
        dimension_ = dimension;
    }

    bool is_layered() const {
        return std::holds_alternative<lqn::LqnModel<double>>(base_);
    }
    std::size_t evaluation_count() const { return count_; }
    std::vector<std::pair<double, double>> bounds() const {
        return std::vector<std::pair<double, double>>(dimension_, {0.0, 1.0});
    }
    VariableValues decode(const std::vector<double>& x) const {
        VariableValues out;
        for (std::size_t i = 0; i < variables_.size(); ++i) {
            const auto begin = x.begin() + offsets_[i];
            out[variables_[i]->name()] = variables_[i]->decode(
                std::vector<double>(begin, begin + variables_[i]->dimension()));
        }
        return out;
    }

    EvaluationResult evaluate(const VariableValues& values) {
        ++count_;
        const auto start = std::chrono::steady_clock::now();
        EvaluationResult out;
        try {
            if (is_layered())
                evaluate_layered(values, out);
            else
                evaluate_flat(values, out);
            out.feasible = true;
        } catch (const std::exception&) {
            out.feasible = false;
        }
        out.solve_time = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        return out;
    }

    std::optional<LayerSensitivityMap> evaluate_layered_sensitivities(
        const VariableValues& values) const {
        if (!is_layered()) return std::nullopt;
        try {
            lqn::LqnModel<double> model = std::get<lqn::LqnModel<double>>(base_);
            apply_layered(model, values);
            lqn::LqnStruct<double> structure = lqn::lqn_finalize(model);
            ln::SolverLN<double> solver(structure, ln::LnOptions());
            const ln::LnSensTable<double> table = solver.get_sensitivity_table(
                sens::SensOptions());
            LayerSensitivityMap out;
            for (const auto& row : table.rows)
                out[metric_key(layered_name(row.station), layered_name(row.jobclass))] = {
                    row.dTput, row.dRespT, row.dQLen, row.dUtil};
            return out;
        } catch (const std::exception&) {
            return std::nullopt;
        }
    }

private:
    static std::string layered_name(const std::string& value) {
        if (value.size() > 2 && value[1] == ':' &&
            std::string("PTERA").find(value[0]) != std::string::npos)
            return value.substr(2);
        return value;
    }

    void evaluate_flat(const VariableValues& values, EvaluationResult& out) const {
        qn::Network<double> model = std::get<qn::Network<double>>(base_);
        for (const auto& fixed : fixed_) fixed.first->apply(model, fixed.second);
        for (const VariablePtr& variable : variables_) {
            const auto found = values.find(variable->name());
            if (found != values.end()) variable->apply(model, found->second);
        }
        const qn::NetworkStruct<double>& structure = model.get_struct();
        const mva::AvgResult<double> average = solve_(structure);
        out.solver_used = average.actualmethod;
        for (std::size_t i = 0; i < structure.nstations; ++i) {
            double utilization = 0.0;
            for (std::size_t r = 0; r < structure.nclasses; ++r) {
                const std::string& station = structure.stations[i].name;
                const std::string& jobclass = structure.classes[r].name;
                if (!average.RN.empty()) out.set_response_time(station, jobclass, average.RN(i, r));
                if (!average.TN.empty()) out.set_throughput(station, jobclass, average.TN(i, r));
                if (!average.QN.empty()) out.set_queue_length(station, jobclass, average.QN(i, r));
                if (!average.UN.empty()) utilization += average.UN(i, r);
            }
            out.utilizations[structure.stations[i].name] = utilization;
        }
        const auto system = solvers::solver_get_avg_sys(structure, average);
        for (std::size_t c = 0; c < structure.nchains; ++c) {
            double response = c < system.CN.size()
                                  ? system.CN[c]
                                  : std::numeric_limits<double>::quiet_NaN();
            const double throughput = c < system.XN.size() ? system.XN[c] : 0.0;
            bool open = false;
            double jobs = 0.0;
            for (const std::size_t member : structure.inchain[c]) {
                open = open || std::isinf(structure.classes[member - 1].population);
                for (std::size_t i = 0; i < structure.nstations; ++i)
                    jobs += out.queue_length(structure.stations[i].name,
                                             structure.classes[member - 1].name);
            }
            if (open && throughput > 0.0) response = jobs / throughput;
            for (const std::size_t member : structure.inchain[c]) {
                out.system_response_times[structure.classes[member - 1].name] = response;
                out.system_throughputs[structure.classes[member - 1].name] = throughput;
            }
        }
        const std::optional<SensitivityData> sensitivity =
            compute_model_sensitivities(model);
        if (sensitivity)
            out.sensitivities = std::make_shared<SensitivityData>(*sensitivity);
    }

    void apply_layered(lqn::LqnModel<double>& model, const VariableValues& values) const {
        for (const auto& fixed : fixed_) fixed.first->apply(model, fixed.second);
        for (const VariablePtr& variable : variables_) {
            const auto found = values.find(variable->name());
            if (found != values.end()) variable->apply(model, found->second);
        }
    }

    void evaluate_layered(const VariableValues& values, EvaluationResult& out) const {
        lqn::LqnModel<double> model = std::get<lqn::LqnModel<double>>(base_);
        apply_layered(model, values);
        const lqn::LqnStruct<double> structure = lqn::lqn_finalize(model);
        ln::SolverLN<double> solver(structure, ln::LnOptions());
        const ln::LnSolution<double> average = solver.get_ensemble_avg();
        out.solver_used = "SolverLN";
        for (std::size_t i = 1; i <= structure.nidx; ++i) {
            const std::string& name = structure.names[i];
            if (i < average.RN.size() && i < average.defined_R.size() &&
                average.defined_R[i] && std::isfinite(average.RN[i]))
                out.set_response_time(name, name, average.RN[i]);
            if (i < average.TN.size() && i < average.defined_T.size() &&
                average.defined_T[i] && std::isfinite(average.TN[i]))
                out.set_throughput(name, name, average.TN[i]);
            if (i < average.QN.size() && i < average.defined_Q.size() &&
                average.defined_Q[i] && std::isfinite(average.QN[i]))
                out.set_queue_length(name, name, average.QN[i]);
            if (i < average.UN.size() && i < average.defined_U.size() &&
                average.defined_U[i] && std::isfinite(average.UN[i]))
                out.utilizations[name] = average.UN[i];
        }
        for (std::size_t task = 1; task <= structure.ntasks; ++task) {
            const std::size_t index = structure.tshift + task;
            if (index >= structure.isref.size() || !structure.isref[index]) continue;
            const std::string& name = structure.names[index];
            const double throughput = out.throughput(name, name);
            if (throughput > 0.0) out.system_throughputs[name] = throughput;
            double response = 0.0;
            bool have_response = false;
            for (const std::size_t entry : structure.entriesof[index]) {
                const double value = out.response_time(structure.names[entry],
                                                       structure.names[entry]);
                if (std::isfinite(value)) {
                    response += value;
                    have_response = true;
                }
            }
            if (have_response) out.system_response_times[name] = response;
        }
    }

    Model base_;
    std::vector<VariablePtr> variables_;
    std::vector<std::pair<VariablePtr, Value>> fixed_;
    SolveFunction solve_;
    std::vector<std::size_t> offsets_;
    std::size_t dimension_ = 0;
    std::size_t count_ = 0;
};

}  // namespace opt
}  // namespace line

#endif
