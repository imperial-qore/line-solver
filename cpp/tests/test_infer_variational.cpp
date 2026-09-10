/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * infer_variational, the Perez-Casale (AAP 53(3), 2021) variational estimator.
 *
 * The estimator carries no random-number stream: the expectations over the
 * other transitions are taken on a Halton lattice mapped through the inverse
 * marginal c.d.f. The fixture below must therefore reproduce the MATLAB,
 * Java and Python implementations digit for digit, and the golden values are
 * the MATLAB reference's.
 *
 * The second case is the one place where the method is exact rather than
 * approximate: with a single transition there is nothing to decouple, so the
 * marginal must reproduce the transient of the underlying birth process.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_variational.h"
#include "line/inference/param_estimator.h"

using line::Matrix;
using line::infer::infer_variational;
using line::infer::VariationalOptions;
using line::infer::VariationalResult;
using line::infer::VariationalSpec;

namespace {

VariationalSpec<double> closed_loop_spec() {
    const double N = 10.0, lam = 0.5;
    const double obsQ[10] = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
    VariationalSpec<double> spec;
    spec.arcs = {{{1, 2, 1}}, {{2, 1, 1}}};
    spec.x0 = Matrix<double>(2, 1, 0.0);
    spec.x0(0, 0) = N;
    spec.sched = {0, 1};
    spec.nservers = {1.0, 1.0};
    spec.routeprob = {1.0, 1.0};
    spec.arcparam = {0, 1};
    spec.arcrate = {lam, 0.0};
    spec.alpha0 = {2.0};
    spec.beta0 = {1.0};
    spec.obsTimes.resize(10);
    spec.obsData = Matrix<double>(10, 2, 0.0);
    for (std::size_t k = 0; k < 10; ++k) {
        spec.obsTimes[k] = static_cast<double>(k + 1);
        spec.obsData(k, 0) = N - obsQ[k];
        spec.obsData(k, 1) = obsQ[k];
    }
    spec.obsRange = {N, N};
    spec.epsilon = 0.1;
    spec.capacity = {N, N};
    return spec;
}

double count_mean(const Matrix<double>& Y, std::size_t g) {
    double m = 0.0;
    for (std::size_t y = 0; y < Y.cols(); ++y) m += Y(g, y) * static_cast<double>(y);
    return m;
}

}  // namespace

TEST_CASE("infer_variational reproduces the MATLAB closed-loop fixture") {
    VariationalOptions<double> opt;
    opt.ngrid = 51;
    opt.nsamples = 32;
    opt.ymax = 60;
    opt.iter_max = 5;
    opt.tol = 0.0;
    opt.delta = 1e-3;
    const VariationalResult<double> out = infer_variational(closed_loop_spec(), opt);

    CHECK(out.alpha[0] == doctest::Approx(20.8738918926).epsilon(0).scale(1).epsilon(1e-9));
    CHECK(out.beta[0] == doctest::Approx(10.4687500000).epsilon(1e-9));
    CHECK(out.rates[0] == doctest::Approx(1.9939240017).epsilon(1e-9));
    CHECK(out.iter == 5);

    const double expected_bound[5] = {-89.392514468, -164.0463074709, -163.7108480833,
                                      -174.6500856726, -163.2336393829};
    REQUIRE(out.bound.size() == 5);
    for (std::size_t i = 0; i < 5; ++i) {
        CHECK(out.bound[i] == doctest::Approx(expected_bound[i]).epsilon(1e-8));
    }

    const std::size_t last = out.Y[0].rows() - 1;
    CHECK(count_mean(out.Y[0], last) == doctest::Approx(22.7916783314).epsilon(1e-9));
    CHECK(count_mean(out.Y[1], last) == doctest::Approx(18.8738918926).epsilon(1e-9));

    const std::size_t gl = out.qlen.rows() - 1;
    CHECK(out.qlen(gl, 0) == doctest::Approx(6.0822135613).epsilon(1e-9));
    CHECK(out.qlen(gl, 1) == doctest::Approx(3.9177864387).epsilon(1e-9));
    // the two stations hold the whole closed population at every epoch
    for (std::size_t g = 0; g < out.qlen.rows(); ++g) {
        CHECK(out.qlen(g, 0) + out.qlen(g, 1) == doctest::Approx(10.0).epsilon(1e-12));
    }
}

TEST_CASE("infer_variational is exact when there is a single transition") {
    const double N = 50.0, lam = 0.1, tmax = 20.0;
    VariationalSpec<double> spec;
    spec.arcs = {{{1, 0, 1}}};
    spec.x0 = Matrix<double>(1, 1, N);
    spec.sched = {0};
    spec.nservers = {1.0};
    spec.routeprob = {1.0};
    spec.arcparam = {0};
    spec.arcrate = {lam};
    spec.obsTimes = {tmax};
    spec.obsData = Matrix<double>(1, 1, VariationalSpec<double>::unobserved());
    spec.obsRange = {N};
    spec.epsilon = 0.2;
    spec.capacity = {N};

    VariationalOptions<double> opt;
    opt.ngrid = 201;
    opt.nsamples = 50;
    opt.iter_max = 4;
    opt.tmax = tmax;
    const VariationalResult<double> out = infer_variational(spec, opt);

    const double exact = N * (1.0 - std::exp(-lam * tmax));
    CHECK(count_mean(out.Y[0], out.Y[0].rows() - 1) == doctest::Approx(exact).epsilon(1e-3));
    CHECK(out.tailmass < 1e-6);
}

TEST_CASE("infer_variational reaches the api fixture through ParamEstimator") {
    // The same fixture driven through the estimator, so that the translation
    // from the network to the transition set is pinned to the api golden:
    // same transitions in the same order, same priors, same observations.
    const double N = 10.0;
    const double obsQ[10] = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
    line::qn::Network<double> model{"vi"};
    const std::size_t think = model.add_delay("Think");
    const std::size_t queue = model.add_queue("Q", line::lang::SchedStrategy::FCFS);
    const std::size_t cl = model.add_closed_class("C", N, think);
    model.set_service(think, cl, line::lang::Distrib<double>::exp_mean(2.0));
    model.set_service(queue, cl, line::lang::Distrib<double>::exp_mean(0.5));
    line::qn::RoutingMatrix<double> routing;
    routing.set(cl, cl, think, queue, 1.0);
    routing.set(cl, cl, queue, think, 1.0);
    model.link(routing);

    std::vector<double> t(10), qd(10), dd(10);
    for (std::size_t k = 0; k < 10; ++k) {
        t[k] = static_cast<double>(k + 1);
        qd[k] = obsQ[k];
        dd[k] = N - obsQ[k];
    }
    line::infer::ParamEstimator est(model);
    est.add_samples(line::infer::SampledMetric(line::lang::MetricType::QLen, t, dd, think, cl));
    est.add_samples(line::infer::SampledMetric(line::lang::MetricType::QLen, t, qd, queue, cl));
    est.options.method = "vi";
    est.options.epsilon = 0.1;
    est.options.prior_shape = 2.0;  // with the model rate 2.0 this gives Gamma(2,1)
    est.options.variational.ngrid = 51;
    est.options.variational.nsamples = 32;
    est.options.variational.ymax = 60;
    est.options.variational.iter_max = 5;
    est.options.variational.tol = 0.0;
    est.options.variational.delta = 1e-3;

    const Matrix<double> out = est.estimate_at({queue});
    CHECK(est.options.posterior_alpha[0] == doctest::Approx(20.8738918926).epsilon(1e-9));
    CHECK(est.options.posterior_beta[0] == doctest::Approx(10.4687500000).epsilon(1e-9));
    CHECK(out(0, 0) == doctest::Approx(10.4687500000 / 20.8738918926).epsilon(1e-9));
}
