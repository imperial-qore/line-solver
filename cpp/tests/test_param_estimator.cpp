/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/inference/param_estimator.h"
#include "line/solvers/mva/solver_mva_runner.h"

using line::Matrix;
using line::infer::ConditionEvent;
using line::infer::EstimatorOptions;
using line::infer::ParamEstimator;
using line::infer::SampledMetric;
using line::lang::Distrib;
using line::lang::EventType;
using line::lang::MetricType;
using line::lang::SchedStrategy;
using line::qn::Network;

namespace {

struct EstimatorModel {
    Network<double> model{"param-estimator"};
    std::size_t queue = model.add_queue("Queue", SchedStrategy::PS);
    std::size_t c1 = model.add_closed_class("Class1", 2.0, queue);
    std::size_t c2 = model.add_closed_class("Class2", 1.0, queue);

    EstimatorModel() {
        model.set_number_of_servers(queue, 2.0);
        model.set_service(queue, c1, Distrib<double>::exp_mean(1.0));
        model.set_service(queue, c2, Distrib<double>::hyperexp(0.5, 0.5, 0.5));
        line::qn::RoutingMatrix<double> routing;
        routing.set(c1, c1, queue, queue, 1.0);
        routing.set(c2, c2, queue, queue, 1.0);
        model.link(routing);
    }
};

struct ClosedEstimatorModel {
    Network<double> model{"closed-param-estimator"};
    std::size_t delay = model.add_delay("Delay");
    std::size_t queue = model.add_queue("Queue", SchedStrategy::PS);
    std::size_t c1 = model.add_closed_class("Class1", 5.0, delay);

    explicit ClosedEstimatorModel(double demand = 0.4) {
        model.set_service(delay, c1, Distrib<double>::exp_mean(1.0));
        model.set_service(queue, c1, Distrib<double>::exp_mean(demand));
        line::qn::RoutingMatrix<double> routing;
        routing.set(c1, c1, delay, queue, 1.0);
        routing.set(c1, c1, queue, delay, 1.0);
        model.link(routing);
    }
};

SampledMetric sample(MetricType type, const std::vector<double>& t,
                     const std::vector<double>& data, std::size_t node, std::size_t cls = 0) {
    return SampledMetric(type, t, data, node, cls);
}

}  // namespace

TEST_CASE("ParamEstimator stores and retrieves class and aggregate samples") {
    EstimatorModel fixture;
    ParamEstimator estimator(fixture.model);
    estimator.add_samples(sample(MetricType::ArvR, {0.0}, {2.0}, fixture.queue, fixture.c1));
    estimator.add_samples(sample(MetricType::Util, {0.0}, {0.4}, fixture.queue));
    SampledMetric q1 = sample(MetricType::QLen, {0.0}, {3.0}, fixture.queue, fixture.c1);
    q1.set_conditional(ConditionEvent(fixture.queue, fixture.c1, EventType::ARV));
    estimator.add_samples(q1);
    estimator.add_samples(sample(MetricType::QLen, {0.0}, {4.0}, fixture.queue, fixture.c1));

    REQUIRE(estimator.get_arvr(fixture.queue, fixture.c1) != nullptr);
    CHECK(estimator.get_arvr(fixture.queue, fixture.c1)->data[0] == doctest::Approx(2.0));
    CHECK(estimator.get_util(fixture.queue, fixture.c1) == nullptr);
    REQUIRE(estimator.get_aggr_util(fixture.queue) != nullptr);
    CHECK(estimator.get_aggr_util(fixture.queue)->is_aggregate());

    const std::vector<SampledMetric*> all = estimator.get_qlen(fixture.queue, fixture.c1);
    CHECK(all.size() == 2);
    const ConditionEvent arrival(fixture.queue, fixture.c1, EventType::ARV);
    const std::vector<SampledMetric*> conditional =
        estimator.get_qlen(fixture.queue, fixture.c1, &arrival);
    REQUIRE(conditional.size() == 1);
    CHECK(conditional[0]->data[0] == doctest::Approx(3.0));
}

TEST_CASE("ParamEstimator interpolation is MATLAB not-a-knot spline with extrapolation") {
    EstimatorModel fixture;
    ParamEstimator estimator(fixture.model);
    estimator.add_samples(sample(MetricType::ArvR, {1.0, 2.0, 3.0, 4.0},
                                 {1.0, 8.0, 27.0, 64.0}, fixture.queue, fixture.c1));
    estimator.add_samples(
        sample(MetricType::Util, {0.0, 5.0}, {10.0, 20.0}, fixture.queue));
    estimator.interpolate();

    const SampledMetric* cubic = estimator.get_arvr(fixture.queue, fixture.c1);
    REQUIRE(cubic != nullptr);
    CHECK((cubic->t == std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}));
    const std::vector<double> cubic_expected{0.0, 1.0, 8.0, 27.0, 64.0, 125.0};
    REQUIRE(cubic->data.size() == cubic_expected.size());
    for (std::size_t i = 0; i < cubic_expected.size(); ++i)
        CHECK(cubic->data[i] == doctest::Approx(cubic_expected[i]).epsilon(1e-12));
    const SampledMetric* linear = estimator.get_aggr_util(fixture.queue);
    REQUIRE(linear != nullptr);
    const std::vector<double> linear_expected{10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
    REQUIRE(linear->data.size() == linear_expected.size());
    for (std::size_t i = 0; i < linear_expected.size(); ++i)
        CHECK(linear->data[i] == doctest::Approx(linear_expected[i]).epsilon(1e-12));
}

TEST_CASE("ParamEstimator NNLS passive solve handles collinear arrival traces") {
    Matrix<double> arrivals(3, 2, 0.0);
    for (std::size_t i = 0; i < arrivals.rows(); ++i) {
        arrivals(i, 0) = static_cast<double>(i + 1);
        arrivals(i, 1) = static_cast<double>(i + 1);
    }
    std::vector<bool> passive{true, true};
    const std::vector<double> demand =
        line::infer::detail::least_squares_passive(arrivals, {1.0, 2.0, 3.0}, passive);

    CHECK(passive[0]);
    CHECK_FALSE(passive[1]);
    CHECK(demand[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(demand[1] == doctest::Approx(0.0));
}

TEST_CASE("ParamEstimator auto method follows MATLAB precedence") {
    EstimatorModel fixture;
    ParamEstimator estimator(fixture.model);
    SampledMetric arrivals =
        sample(MetricType::ArvR, {0.0}, {1.0}, fixture.queue, fixture.c1);
    arrivals.set_trace();
    SampledMetric response =
        sample(MetricType::RespT, {0.0}, {1.0}, fixture.queue, fixture.c1);
    response.set_trace();
    estimator.add_samples(arrivals);
    estimator.add_samples(response);
    estimator.add_samples(sample(MetricType::Util, {0.0}, {0.5}, fixture.queue));
    estimator.add_samples(sample(MetricType::QLen, {0.0}, {1.0}, fixture.queue));

    CHECK(estimator.auto_method() == "erps");
    CHECK(estimator.options.method == "erps");
    CHECK(estimator.has_estimator("erps"));
}

TEST_CASE("ParamEstimator UBR estimates demands and preserves service families") {
    EstimatorModel fixture;
    ParamEstimator estimator(fixture.model);
    const std::vector<double> t{0.0, 1.0, 2.0};
    estimator.add_samples(sample(MetricType::ArvR, t, {1.0, 2.0, 3.0}, fixture.queue,
                                 fixture.c1));
    estimator.add_samples(sample(MetricType::Util, t, {0.2, 0.4, 0.6}, fixture.queue,
                                 fixture.c1));
    estimator.add_samples(sample(MetricType::ArvR, t, {1.0, 0.5, 1.0}, fixture.queue,
                                 fixture.c2));
    estimator.add_samples(
        sample(MetricType::Util, t, {0.325, 0.4625, 0.725}, fixture.queue));

    const Matrix<double> result = estimator.estimate_at({fixture.queue});
    REQUIRE(result.rows() == 1);
    REQUIRE(result.cols() == 2);
    CHECK(result(0, 0) == doctest::Approx(0.4).epsilon(1e-12));
    CHECK(result(0, 1) == doctest::Approx(0.25).epsilon(1e-12));

    const auto& sn = fixture.model.raw_struct();
    const std::size_t station = sn.nodes[fixture.queue - 1].station;
    CHECK(sn.service[station - 1][fixture.c1 - 1].type == line::lang::ProcessType::EXP);
    CHECK(sn.service[station - 1][fixture.c2 - 1].type == line::lang::ProcessType::HYPEREXP);
    CHECK(sn.service[station - 1][fixture.c1 - 1].mean == doctest::Approx(0.4));
    CHECK(sn.service[station - 1][fixture.c2 - 1].mean == doctest::Approx(0.25));

    line::mva::MvaOptions options;
    Matrix<double> init;
    const auto solved = line::mva::solver_mva_run_analyzer(fixture.model.get_struct(), options, init);
    CHECK(solved.QN(0, 0) + solved.QN(0, 1) == doctest::Approx(3.0).epsilon(1e-10));
    CHECK(solved.UN(0, 0) == doctest::Approx(solved.TN(0, 0) * 0.4 / 2.0).epsilon(1e-10));
    CHECK(solved.UN(0, 1) == doctest::Approx(solved.TN(0, 1) * 0.25 / 2.0).epsilon(1e-10));
}

TEST_CASE("ParamEstimator model update rejects a family without MATLAB setMean") {
    EstimatorModel fixture;
    fixture.model.set_service(fixture.queue, fixture.c2, Distrib<double>::det(2.0));
    EstimatorOptions options;
    options.method = "fixed";
    ParamEstimator estimator(fixture.model, options);
    estimator.register_estimator("fixed", [](ParamEstimator&, const std::vector<std::size_t>&) {
        Matrix<double> values(1, 2, 1.0);
        return values;
    });
    CHECK_THROWS_WITH_AS(estimator.estimate_at({fixture.queue}),
                         "ParamEstimator.estimate_at: Det has no MATLAB setMean contract",
                         line::UnsupportedError);
}

TEST_CASE("ParamEstimator typed dispatch updates every requested station") {
    EstimatorModel fixture;
    const std::size_t queue2 = fixture.model.add_queue("Queue2", SchedStrategy::FCFS);
    fixture.model.set_service(queue2, fixture.c1, Distrib<double>::exp_mean(3.0));
    fixture.model.set_service(queue2, fixture.c2, Distrib<double>::exp_mean(4.0));
    EstimatorOptions options;
    options.method = "fixed";
    ParamEstimator estimator(fixture.model, options);
    estimator.register_estimator("fixed", [](ParamEstimator&, const std::vector<std::size_t>&) {
        Matrix<double> values(2, 2, 0.0);
        values(0, 0) = 0.5;
        values(0, 1) = 0.75;
        values(1, 0) = 1.25;
        values(1, 1) = 1.5;
        return values;
    });

    const Matrix<double> result = estimator.estimate_at({fixture.queue, queue2});
    CHECK(result(1, 1) == doctest::Approx(1.5));
    const auto& sn = fixture.model.raw_struct();
    const std::size_t s1 = sn.nodes[fixture.queue - 1].station;
    const std::size_t s2 = sn.nodes[queue2 - 1].station;
    CHECK(sn.service[s1 - 1][fixture.c1 - 1].mean == doctest::Approx(0.5));
    CHECK(sn.service[s1 - 1][fixture.c2 - 1].mean == doctest::Approx(0.75));
    CHECK(sn.service[s2 - 1][fixture.c1 - 1].mean == doctest::Approx(1.25));
    CHECK(sn.service[s2 - 1][fixture.c2 - 1].mean == doctest::Approx(1.5));
}

TEST_CASE("ParamEstimator exposes MATLAB defaults and metric descriptions") {
    const EstimatorOptions options = ParamEstimator::default_options();
    CHECK(options.method == "ubr");
    CHECK(options.variant == "default");
    CHECK(options.iter_max == 1000);
    CHECK(options.tol == doctest::Approx(1e-3));
    CHECK(options.open_population == 100);
    CHECK(ParamEstimator::get_required_metrics("ubr") ==
          "ArvR (per-class) + Util (per-class or aggregate)");
    CHECK(ParamEstimator::get_required_metrics("missing") == "Unknown method: missing");
}

TEST_CASE("ParamEstimator installs every MATLAB estimator adapter") {
    EstimatorModel fixture;
    ParamEstimator estimator(fixture.model);
    const std::vector<std::string> methods{"ubr",   "ubo",  "erps",  "ekf",  "mcmc",
                                           "mle",   "mlps", "fmlps", "qmle", "gibbs"};
    for (const std::string& method : methods) CHECK(estimator.has_estimator(method));
}

TEST_CASE("ParamEstimator QMLE marshals queue means and closed population") {
    ClosedEstimatorModel fixture;
    EstimatorOptions options;
    options.method = "qmle";
    ParamEstimator estimator(fixture.model, options);
    estimator.add_samples(sample(MetricType::QLen, {0.0, 1.0, 2.0}, {2.0, 2.0, 2.0},
                                 fixture.queue, fixture.c1));
    const Matrix<double> result = estimator.estimate_at({fixture.queue});
    CHECK(result(0, 0) == doctest::Approx(2.0 / 3.0 / 2.6).epsilon(1e-12));
}

TEST_CASE("ParamEstimator UBO solves the reference quadratic as weighted NNLS") {
    EstimatorModel fixture;
    EstimatorOptions options;
    options.method = "ubo";
    ParamEstimator estimator(fixture.model, options);
    const std::vector<double> t{0.0, 1.0, 2.0};
    estimator.add_samples(sample(MetricType::ArvR, t, {0.5, 1.0, 0.5}, fixture.queue,
                                 fixture.c1));
    estimator.add_samples(sample(MetricType::ArvR, t, {0.25, 0.5, 0.25}, fixture.queue,
                                 fixture.c2));
    estimator.add_samples(sample(MetricType::RespT, t, {0.5333333333333333, 0.8,
                                                         0.5333333333333333},
                                 fixture.queue, fixture.c1));
    estimator.add_samples(sample(MetricType::RespT, t, {0.2666666666666667, 0.4,
                                                         0.2666666666666667},
                                 fixture.queue, fixture.c2));
    estimator.add_samples(sample(MetricType::Util, t, {0.125, 0.25, 0.125}, fixture.queue));
    const Matrix<double> result = estimator.estimate_at({fixture.queue});
    CHECK(result(0, 0) == doctest::Approx(0.4).epsilon(1e-10));
    CHECK(result(0, 1) == doctest::Approx(0.2).epsilon(1e-10));
}

TEST_CASE("ParamEstimator ERPS uses arrival-conditional aggregate queue length") {
    EstimatorModel fixture;
    EstimatorOptions options;
    options.method = "erps";
    ParamEstimator estimator(fixture.model, options);
    const std::vector<double> t{0.0, 1.0};
    estimator.add_samples(sample(MetricType::RespT, t, {0.2, 0.4}, fixture.queue,
                                 fixture.c1));
    estimator.add_samples(sample(MetricType::RespT, t, {0.8, 0.4}, fixture.queue,
                                 fixture.c2));
    SampledMetric q1 = sample(MetricType::QLen, t, {1.0, 2.0}, fixture.queue);
    q1.set_conditional(ConditionEvent(fixture.queue, fixture.c1, EventType::ARV));
    estimator.add_samples(q1);
    SampledMetric q2 = sample(MetricType::QLen, t, {2.0, 1.0}, fixture.queue);
    q2.set_conditional(ConditionEvent(fixture.queue, fixture.c2, EventType::ARV));
    estimator.add_samples(q2);
    const Matrix<double> result = estimator.estimate_at({fixture.queue});
    CHECK(result(0, 0) == doctest::Approx(0.3).epsilon(1e-12));
    CHECK(result(0, 1) == doctest::Approx(0.6).epsilon(1e-12));
}

TEST_CASE("ParamEstimator MCMC is seeded and preserves the node-class shape") {
    ClosedEstimatorModel fixture;
    EstimatorOptions options;
    options.method = "mcmc";
    options.random_seed = 17;
    ParamEstimator first(fixture.model, options);
    first.add_samples(sample(MetricType::QLen, {0.0, 1.0, 2.0, 3.0},
                               {1.0, 2.0, 1.5, 2.5}, fixture.queue));
    const Matrix<double> a = first.estimate_at({fixture.queue});

    ClosedEstimatorModel fixture2;
    ParamEstimator second(fixture2.model, options);
    second.add_samples(sample(MetricType::QLen, {0.0, 1.0, 2.0, 3.0},
                                {1.0, 2.0, 1.5, 2.5}, fixture2.queue));
    const Matrix<double> b = second.estimate_at({fixture2.queue});
    CHECK(a.rows() == 1);
    CHECK(a.cols() == 1);
    CHECK(a(0, 0) > 0.0);
    CHECK(a(0, 0) == doctest::Approx(b(0, 0)).epsilon(1e-15));
}

TEST_CASE("ParamEstimator MLE fits MVA response time and busy servers") {
    ClosedEstimatorModel truth(0.4);
    line::mva::MvaOptions mva_options;
    Matrix<double> init;
    const auto measured = line::mva::solver_mva_run_analyzer(truth.model.get_struct(), mva_options, init);

    ClosedEstimatorModel target(1.0);
    EstimatorOptions options;
    options.method = "mle";
    options.iter_max = 300;
    options.x0 = {0.2};
    ParamEstimator estimator(target.model, options);
    const std::vector<double> t{0.0, 1.0, 2.0};
    estimator.add_samples(sample(MetricType::ArvR, t, {1.0, 1.0, 1.0}, target.queue,
                                 target.c1));
    estimator.add_samples(sample(MetricType::RespT, t,
                                 {measured.RN(1, 0), measured.RN(1, 0), measured.RN(1, 0)},
                                 target.queue, target.c1));
    estimator.add_samples(sample(MetricType::Util, t,
                                 {measured.UN(1, 0), measured.UN(1, 0), measured.UN(1, 0)},
                                 target.queue));
    const Matrix<double> result = estimator.estimate_at({target.queue});
    CHECK(result(0, 0) == doctest::Approx(0.4).epsilon(1e-5));
}

TEST_CASE("ParamEstimator EKF consumes MVA measurements sequentially") {
    ClosedEstimatorModel truth(0.4);
    line::mva::MvaOptions mva_options;
    Matrix<double> init;
    const auto measured = line::mva::solver_mva_run_analyzer(truth.model.get_struct(), mva_options, init);

    ClosedEstimatorModel target(1.0);
    EstimatorOptions options;
    options.method = "ekf";
    options.iter_max = 8;
    options.x0 = {0.2};
    ParamEstimator estimator(target.model, options);
    const std::vector<double> t{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    estimator.add_samples(sample(MetricType::ArvR, t, std::vector<double>(8, 1.0),
                                 target.queue, target.c1));
    estimator.add_samples(sample(MetricType::RespT, t,
                                 std::vector<double>(8, measured.RN(1, 0)), target.queue,
                                 target.c1));
    estimator.add_samples(sample(MetricType::Util, t,
                                 std::vector<double>(8, measured.UN(1, 0)), target.queue));
    const Matrix<double> result = estimator.estimate_at({target.queue});
    CHECK(result(0, 0) == doctest::Approx(0.4).epsilon(0.03));
}

TEST_CASE("ParamEstimator trace adapters execute Gibbs MLPS and FMLPS leaves") {
    const std::vector<double> clock{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const std::vector<double> arrivals{0.1, 1.1, 2.1, 3.1, 4.1, 5.1};
    const std::vector<double> response(6, 0.3);

    for (const std::string& method : {std::string("mlps"), std::string("fmlps"),
                                      std::string("gibbs")}) {
        ClosedEstimatorModel fixture;
        EstimatorOptions options;
        options.method = method;
        options.random_seed = 3;
        options.tol = 0.05;
        options.gibbs.data_needed = 0;
        options.gibbs.likelihood_sample = 8;
        options.gibbs.nsamples = 4;
        options.gibbs.block = 2;
        ParamEstimator estimator(fixture.model, options);
        SampledMetric arv = sample(MetricType::ArvR, clock, arrivals, fixture.queue,
                                   fixture.c1);
        arv.set_trace();
        SampledMetric resp = sample(MetricType::RespT, clock, response, fixture.queue,
                                    fixture.c1);
        resp.set_trace();
        estimator.add_samples(arv);
        estimator.add_samples(resp);
        if (method == "gibbs")
            estimator.add_samples(sample(MetricType::Tput, clock,
                                         std::vector<double>(clock.size(), 5.0 / 1.4),
                                         fixture.queue, fixture.c1));
        const Matrix<double> result = estimator.estimate_at({fixture.queue});
        INFO(method);
        CHECK(result.rows() == 1);
        CHECK(result.cols() == 1);
        CHECK(result(0, 0) > 0.0);
        CHECK(std::isfinite(result(0, 0)));
    }
}
