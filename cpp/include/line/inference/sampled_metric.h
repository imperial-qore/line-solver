/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_INFERENCE_SAMPLED_METRIC_H
#define LINE_INFERENCE_SAMPLED_METRIC_H

/**
 * Observed metric data supplied to parameter estimators.
 *
 * Port of MATLAB `SampledMetric`, with the JAR's explicit format and condition
 * types. Nodes and classes use the C++ model builder's 1-based indices; class
 * zero is the aggregate sentinel corresponding to an empty MATLAB class.
 */

#include <cstddef>
#include <optional>
#include <vector>

#include "line/lang/lang_types.h"

namespace line {
namespace infer {

enum class SampledFormat { TIMESERIES, TRACE };

struct ConditionEvent {
    std::size_t node = 0;
    std::size_t jobclass = 0;
    lang::EventType event = lang::EventType::INIT;

    ConditionEvent() = default;
    ConditionEvent(std::size_t node_index, std::size_t class_index, lang::EventType event_type)
        : node(node_index), jobclass(class_index), event(event_type) {}

    bool operator==(const ConditionEvent& other) const {
        return node == other.node && jobclass == other.jobclass && event == other.event;
    }
    bool operator!=(const ConditionEvent& other) const { return !(*this == other); }
};

class SampledMetric {
  public:
    lang::MetricType type;
    std::vector<double> t;
    std::vector<double> data;
    std::size_t node;
    std::size_t jobclass;
    std::optional<ConditionEvent> cond;
    SampledFormat format;

    SampledMetric(lang::MetricType metric_type, const std::vector<double>& times,
                  const std::vector<double>& observations, std::size_t node_index,
                  std::size_t class_index = 0)
        : type(metric_type),
          t(times),
          data(observations),
          node(node_index),
          jobclass(class_index),
          format(SampledFormat::TIMESERIES) {}

    void set_conditional(const ConditionEvent& event) { cond = event; }
    void set_trace() { format = SampledFormat::TRACE; }

    bool is_aggregate() const { return jobclass == 0; }
    bool is_conditional() const { return cond.has_value(); }
    bool is_trace() const { return format == SampledFormat::TRACE; }

    SampledMetric copy() const { return *this; }
};

}  // namespace infer
}  // namespace line

#endif  // LINE_INFERENCE_SAMPLED_METRIC_H
