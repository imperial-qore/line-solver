/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Liu-Whitt time-varying many-server fluid queue and the network of them.
 *
 * The checks are against models that are exactly known: with constant data the
 * trajectory must settle on the stationary G/GI/s+GI fluid model of Whitt
 * (2006); with unlimited staffing it must reproduce the exact Mt/M/infinity
 * mean; and the network must satisfy its own traffic equations, which for a
 * tandem and for a feedback loop have closed-form fixed points.
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_ggisgi_fluid.h"
#include "line/api/qsys/qsys_gtmtst_fluid.h"
#include "line/api/qsys/qsys_mtginf.h"

using line::qsys::TvFluidOptions;
using F = std::function<double(const double&)>;

namespace {
const double THETA = 0.5;
F patience() {
    return [](const double& x) { return std::exp(-THETA * x); };
}
F patiencePdf() {
    return [](const double& x) { return THETA * std::exp(-THETA * x); };
}
F constFn(double c) {
    return [c](const double&) { return c; };
}
}  // namespace

TEST_CASE("with constant data the fluid queue settles on the stationary fluid model") {
    TvFluidOptions<double> o;
    o.dt = 0.01;
    o.pdf = patiencePdf();
    const auto r = line::qsys::qsys_gtmtst_fluid<double>(constFn(110.0), constFn(100.0),
                                                         constFn(1.0), patience(), 60.0, o);
    const auto st = line::qsys::qsys_ggisgi_fluid<double>(110.0, 1.0, 100u, patience());
    const std::size_t n = r.times.size() - 1;
    CHECK(r.regime[n] == 1);
    CHECK(r.w[n] == doctest::Approx(st.offeredWait).epsilon(1e-6));
    CHECK(r.Q[n] == doctest::Approx(st.meanQueueLength).epsilon(1e-6));
    CHECK(r.B[n] == doctest::Approx(st.meanNumberInService).epsilon(1e-9));
    // In equilibrium everything that cannot be served abandons.
    CHECK(r.alpha[n] == doctest::Approx(110.0 - 100.0).epsilon(1e-6));
    CHECK(r.utilization[n] == doctest::Approx(1.0));
}

TEST_CASE("with unlimited staffing the fluid queue is the exact Mt/M/infinity mean") {
    F lam = [](const double& t) { return 10.0 + 5.0 * std::sin(t); };
    TvFluidOptions<double> o;
    o.dt = 0.002;
    const auto r = line::qsys::qsys_gtmtst_fluid<double>(lam, constFn(1e6), constFn(2.0),
                                                         patience(), 20.0, o);
    std::vector<double> ts = {1.0, 5.0, 10.0, 20.0};
    F gc = [](const double& x) { return std::exp(-2.0 * x); };
    const auto exact = line::qsys::qsys_mtginf<double>(lam, gc, 0.5, ts, 0.0);
    for (std::size_t k = 0; k < ts.size(); ++k) {
        const std::size_t i = static_cast<std::size_t>(std::llround(ts[k] / 0.002));
        CHECK(r.B[i] == doctest::Approx(exact.meanNumber[k]).epsilon(1e-5));
        CHECK(r.regime[i] == 0);          // never overloaded with that much capacity
        CHECK(r.Q[i] == doctest::Approx(0.0));
    }
}

TEST_CASE("the fluid queue switches regime when the servers fill") {
    // A rate that starts below capacity and rises above it must produce exactly
    // one switch into overload, and the queue must be empty before it.
    F lam = [](const double& t) { return 50.0 + 20.0 * t; };
    TvFluidOptions<double> o;
    o.dt = 0.005;
    const auto r = line::qsys::qsys_gtmtst_fluid<double>(lam, constFn(100.0), constFn(1.0),
                                                         patience(), 8.0, o);
    int switches = 0;
    for (std::size_t i = 1; i < r.regime.size(); ++i)
        if (r.regime[i] != r.regime[i - 1]) ++switches;
    CHECK(switches == 1);
    CHECK(r.regime.front() == 0);
    CHECK(r.regime.back() == 1);
    CHECK(r.Q.front() == doctest::Approx(0.0));
    CHECK(r.Q.back() > 0.0);
    CHECK(r.w.back() > 0.0);
}

TEST_CASE("the tandem network routes exactly the service completion flow") {
    std::vector<F> lams = {constFn(110.0), constFn(0.0)};
    std::vector<F> ss = {constFn(100.0), constFn(80.0)};
    std::vector<F> mus = {constFn(1.0), constFn(1.0)};
    std::vector<F> fcs = {patience(), patience()};
    std::vector<std::vector<double>> P = {{0.0, 1.0}, {0.0, 0.0}};
    const auto net = line::qsys::npfqn_gtmtst_fluid<double>(lams, ss, mus, fcs, P, 60.0, 0.02, {},
                                                            {}, 1e-9, 100);
    const std::size_t n = net.times.size() - 1;
    // Queue 1 is overloaded, so it completes exactly s*mu and that is what
    // reaches queue 2.
    CHECK(net.arrivalRates[1][n] == doctest::Approx(net.queues[0].sigma[n]).epsilon(1e-12));
    CHECK(net.arrivalRates[1][n] == doctest::Approx(100.0).epsilon(1e-9));
    // Queue 2 must then agree with the stationary fluid model fed at that rate.
    const auto ref = line::qsys::qsys_ggisgi_fluid<double>(100.0, 1.0, 80u, patience());
    CHECK(net.queues[1].Q[n] == doctest::Approx(ref.meanQueueLength).epsilon(1e-6));
    CHECK(net.queues[1].w[n] == doctest::Approx(ref.offeredWait).epsilon(1e-6));
}

TEST_CASE("a feedback loop converges to its closed-form arrival rate") {
    // 30% of the completions return, so an underloaded queue settles at
    // lambda = 50/(1-0.3).
    std::vector<F> lams = {constFn(50.0)};
    std::vector<F> ss = {constFn(100.0)};
    std::vector<F> mus = {constFn(1.0)};
    std::vector<F> fcs = {patience()};
    std::vector<std::vector<double>> P = {{0.3}};
    const auto net = line::qsys::npfqn_gtmtst_fluid<double>(lams, ss, mus, fcs, P, 40.0, 0.02, {},
                                                            {}, 1e-10, 100);
    const std::size_t n = net.times.size() - 1;
    CHECK(net.arrivalRates[0][n] == doctest::Approx(50.0 / 0.7).epsilon(1e-8));
    CHECK(net.queues[0].B[n] == doctest::Approx(50.0 / 0.7).epsilon(1e-6));
    CHECK(net.queues[0].regime[n] == 0);
    CHECK(net.iterations > 1);
}

TEST_CASE("the routing matrix must be substochastic") {
    std::vector<F> lams = {constFn(1.0)};
    std::vector<F> ss = {constFn(10.0)};
    std::vector<F> mus = {constFn(1.0)};
    std::vector<F> fcs = {patience()};
    std::vector<std::vector<double>> bad = {{1.5}};
    CHECK_THROWS(line::qsys::npfqn_gtmtst_fluid<double>(lams, ss, mus, fcs, bad, 1.0, 0.1, {}, {},
                                                        1e-6, 10));
}
