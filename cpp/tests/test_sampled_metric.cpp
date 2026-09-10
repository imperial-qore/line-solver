/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/inference/sampled_metric.h"

using line::infer::ConditionEvent;
using line::infer::SampledFormat;
using line::infer::SampledMetric;
using line::lang::EventType;
using line::lang::MetricType;

TEST_CASE("SampledMetric class and aggregate forms match the reference") {
    SampledMetric per_class(MetricType::ArvR, std::vector<double>{0.0, 1.0},
                            std::vector<double>{2.0, 3.0}, 2, 1);
    CHECK(per_class.type == MetricType::ArvR);
    CHECK((per_class.t == std::vector<double>{0.0, 1.0}));
    CHECK((per_class.data == std::vector<double>{2.0, 3.0}));
    CHECK(per_class.node == 2);
    CHECK(per_class.jobclass == 1);
    CHECK(per_class.format == SampledFormat::TIMESERIES);
    CHECK_FALSE(per_class.is_aggregate());
    CHECK_FALSE(per_class.is_conditional());
    CHECK_FALSE(per_class.is_trace());

    SampledMetric aggregate(MetricType::Util, std::vector<double>{0.0},
                            std::vector<double>{0.75}, 2);
    CHECK(aggregate.is_aggregate());
    CHECK(aggregate.jobclass == 0);
}

TEST_CASE("SampledMetric condition and trace flags are independent") {
    SampledMetric sample(MetricType::QLen, std::vector<double>{0.0},
                         std::vector<double>{4.0}, 3);
    const ConditionEvent arrival(3, 2, EventType::ARV);
    sample.set_conditional(arrival);
    REQUIRE(sample.cond.has_value());
    CHECK(*sample.cond == arrival);
    CHECK(sample.is_conditional());
    CHECK_FALSE(sample.is_trace());

    sample.set_trace();
    CHECK(sample.is_conditional());
    CHECK(sample.is_trace());
    CHECK(sample.format == SampledFormat::TRACE);
}

TEST_CASE("SampledMetric copy owns its sampled arrays") {
    SampledMetric original(MetricType::RespT, std::vector<double>{1.0, 2.0},
                           std::vector<double>{0.2, 0.4}, 4, 1);
    original.set_conditional(ConditionEvent(4, 1, EventType::DEP));
    original.set_trace();
    SampledMetric clone = original.copy();

    original.t[0] = 99.0;
    original.data[0] = 88.0;
    original.cond = ConditionEvent(7, 3, EventType::ARV);

    CHECK((clone.t == std::vector<double>{1.0, 2.0}));
    CHECK((clone.data == std::vector<double>{0.2, 0.4}));
    REQUIRE(clone.cond.has_value());
    CHECK((*clone.cond == ConditionEvent(4, 1, EventType::DEP)));
    CHECK(clone.is_trace());
}

TEST_CASE("MetricType preserves MATLAB numbering and text") {
    CHECK(static_cast<int>(MetricType::ResidT) == 0);
    CHECK(static_cast<int>(MetricType::RespT) == 1);
    CHECK(static_cast<int>(MetricType::QLen) == 3);
    CHECK(static_cast<int>(MetricType::Tput) == 13);
    CHECK(static_cast<int>(MetricType::ArvR) == 14);
    CHECK(static_cast<int>(MetricType::Util) == 16);
    CHECK(static_cast<int>(MetricType::SysTard) == 22);
    CHECK(std::string(line::lang::metric_to_text(MetricType::QLen)) == "Number of Customers");
    CHECK(std::string(line::lang::metric_to_text(MetricType::ArvR)) == "Arrival Rate");
    CHECK(std::string(line::lang::metric_to_text(MetricType::SysTard)) == "System Tardiness");
}
