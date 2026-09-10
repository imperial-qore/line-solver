/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_LAYERED_VARIABLES_H
#define LINE_OPT_LAYERED_VARIABLES_H

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "line/opt/variables.h"

namespace line {
namespace opt {

namespace layered_detail {

inline std::size_t processor(const lqn::LqnModel<double>& model, const std::string& name) {
    for (std::size_t i = 0; i < model.procs.size(); ++i)
        if (model.procs[i].name == name) return i;
    return model.procs.size();
}
inline std::size_t task(const lqn::LqnModel<double>& model, const std::string& name) {
    for (std::size_t i = 0; i < model.tasks.size(); ++i)
        if (model.tasks[i].name == name) return i;
    return model.tasks.size();
}
inline std::size_t activity(const lqn::LqnModel<double>& model, const std::string& name) {
    for (std::size_t i = 0; i < model.acts.size(); ++i)
        if (model.acts[i].name == name) return i;
    return model.acts.size();
}
inline std::string task_processor(const lqn::LqnModel<double>& model,
                                  const std::string& task_name) {
    const std::size_t index = task(model, task_name);
    if (index >= model.tasks.size() || model.tasks[index].proc_slot >= model.procs.size())
        return {};
    return model.procs[model.tasks[index].proc_slot].name;
}
inline std::string activity_task(const lqn::LqnModel<double>& model,
                                 const std::string& activity_name) {
    const std::size_t index = activity(model, activity_name);
    if (index >= model.acts.size() || model.acts[index].task_slot >= model.tasks.size())
        return {};
    return model.tasks[model.acts[index].task_slot].name;
}
inline std::string activity_processor(const lqn::LqnModel<double>& model,
                                      const std::string& activity_name) {
    return task_processor(model, activity_task(model, activity_name));
}
inline Value integer_decode(const std::vector<double>& x, int low, int high) {
    return {static_cast<double>(std::clamp<int>(
        static_cast<int>(std::llround(low + x.at(0) * (high - low))), low, high))};
}
inline Value continuous_decode(const std::vector<double>& x, double low, double high) {
    return {low + x.at(0) * (high - low)};
}

}  // namespace layered_detail

class ProcessorMultiplicity final : public DecisionVariable {
public:
    ProcessorMultiplicity(std::string processor, int low, int high, std::string name = "")
        : DecisionVariable(name.empty() ? processor + "_multiplicity" : std::move(name)),
          processor_(std::move(processor)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::integer_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::processor(model, processor_);
        if (i < model.procs.size()) model.procs[i].mult = std::round(scalar_value(value));
    }
    std::string type() const override { return "processor_multiplicity"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>&) const override {
        return {processor_};
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::processor(model, processor_);
        if (i >= model.procs.size() || !std::isfinite(model.procs[i].mult)) return std::nullopt;
        return Value{std::round(model.procs[i].mult)};
    }
private:
    std::string processor_;
    int low_, high_;
};

class TaskMultiplicity final : public DecisionVariable {
public:
    TaskMultiplicity(std::string task, int low, int high, std::string name = "")
        : DecisionVariable(name.empty() ? task + "_multiplicity" : std::move(name)),
          task_(std::move(task)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::integer_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i < model.tasks.size()) model.tasks[i].mult = std::round(scalar_value(value));
    }
    std::string type() const override { return "task_multiplicity"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>& model) const override {
        std::vector<std::string> out{task_};
        const std::string processor = layered_detail::task_processor(model, task_);
        if (!processor.empty()) out.push_back(processor);
        return out;
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i >= model.tasks.size() || !std::isfinite(model.tasks[i].mult)) return std::nullopt;
        return Value{std::round(model.tasks[i].mult)};
    }
private:
    std::string task_;
    int low_, high_;
};

class TaskReplication final : public DecisionVariable {
public:
    TaskReplication(std::string task, int low, int high, std::string name = "")
        : DecisionVariable(name.empty() ? task + "_replication" : std::move(name)),
          task_(std::move(task)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::integer_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i < model.tasks.size()) model.tasks[i].repl = std::round(scalar_value(value));
    }
    std::string type() const override { return "task_replication"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>& model) const override {
        std::vector<std::string> out{task_};
        const std::string processor = layered_detail::task_processor(model, task_);
        if (!processor.empty()) out.push_back(processor);
        return out;
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i >= model.tasks.size() || !std::isfinite(model.tasks[i].repl)) return std::nullopt;
        return Value{std::round(model.tasks[i].repl)};
    }
private:
    std::string task_;
    int low_, high_;
};

class HostDemand final : public DecisionVariable {
public:
    HostDemand(std::string activity, double low, double high, std::string name = "")
        : DecisionVariable(name.empty() ? activity + "_hostdemand" : std::move(name)),
          activity_(std::move(activity)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::continuous_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::activity(model, activity_);
        if (i < model.acts.size())
            model.acts[i].hostdem = lang::Distrib<double>::exp_mean(scalar_value(value));
    }
    std::string type() const override { return "host_demand"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>& model) const override {
        const std::string processor = layered_detail::activity_processor(model, activity_);
        return processor.empty() ? std::vector<std::string>()
                                 : std::vector<std::string>{processor};
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::activity(model, activity_);
        if (i >= model.acts.size() || model.acts[i].hostdem.disabled) return std::nullopt;
        return Value{model.acts[i].hostdem.mean};
    }
    bool supports_sensitivity() const override { return true; }
    std::string sensitivity_key(const lqn::LqnModel<double>& model) const override {
        const std::string processor = layered_detail::activity_processor(model, activity_);
        const std::string task = layered_detail::activity_task(model, activity_);
        return processor.empty() || task.empty() ? std::string()
                                                  : processor + "||" + task;
    }
    std::map<std::string,std::string> sensitivity_metric_targets(
        const lqn::LqnModel<double>& model) const override {
        const std::string processor = layered_detail::activity_processor(model, activity_);
        return {{"Util", processor}, {"Tput", metric_key(activity_, activity_)},
                {"QLen", metric_key(activity_, activity_)},
                {"RespT", metric_key(activity_, activity_)}};
    }
    double rate_jacobian(const Value& value) const override {
        const double demand = scalar_value(value);
        return demand <= 0.0 ? 0.0 : -1.0 / (demand * demand);
    }
    double decode_jacobian(double) const override { return high_ - low_; }
private:
    std::string activity_;
    double low_, high_;
};

class TaskThinkTime final : public DecisionVariable {
public:
    TaskThinkTime(std::string task, double low, double high, std::string name = "")
        : DecisionVariable(name.empty() ? task + "_thinktime" : std::move(name)),
          task_(std::move(task)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::continuous_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i < model.tasks.size())
            model.tasks[i].thinktime = lang::Distrib<double>::exp_mean(scalar_value(value));
    }
    std::string type() const override { return "think_time"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>& model) const override {
        std::vector<std::string> out{task_};
        const std::string processor = layered_detail::task_processor(model, task_);
        if (!processor.empty()) out.push_back(processor);
        return out;
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::task(model, task_);
        if (i >= model.tasks.size() || model.tasks[i].thinktime.disabled) return std::nullopt;
        return Value{model.tasks[i].thinktime.mean};
    }
private:
    std::string task_;
    double low_, high_;
};

class ActivityThinkTime final : public DecisionVariable {
public:
    ActivityThinkTime(std::string activity, double low, double high, std::string name = "")
        : DecisionVariable(name.empty() ? activity + "_thinktime" : std::move(name)),
          activity_(std::move(activity)), low_(low), high_(high) {}
    Value decode(const std::vector<double>& x) const override {
        return layered_detail::continuous_decode(x, low_, high_);
    }
    void apply(lqn::LqnModel<double>& model, const Value& value) const override {
        const std::size_t i = layered_detail::activity(model, activity_);
        if (i < model.acts.size())
            model.acts[i].thinktime = lang::Distrib<double>::exp_mean(scalar_value(value));
    }
    std::string type() const override { return "think_time"; }
    std::vector<std::string> layers(const lqn::LqnModel<double>& model) const override {
        std::vector<std::string> out;
        const std::string task = layered_detail::activity_task(model, activity_);
        const std::string processor = layered_detail::activity_processor(model, activity_);
        if (!task.empty()) out.push_back(task);
        if (!processor.empty()) out.push_back(processor);
        return out;
    }
    std::optional<Value> current_value(const lqn::LqnModel<double>& model) const override {
        const std::size_t i = layered_detail::activity(model, activity_);
        if (i >= model.acts.size() || model.acts[i].thinktime.disabled) return std::nullopt;
        return Value{model.acts[i].thinktime.mean};
    }
private:
    std::string activity_;
    double low_, high_;
};

}  // namespace opt
}  // namespace line

#endif
