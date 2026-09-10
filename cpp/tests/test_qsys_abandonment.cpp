/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Whitt family of multiserver models: abandonment (Erlang A and the
 * M/GI/s/r+GI engineering solution), the G/GI/s+GI fluid limit, the Halfin-Whitt
 * QED regime and the Mt/G/infinity queue.
 *
 * The checks are IDENTITIES rather than stored numbers wherever one exists,
 * because each of these is exact somewhere: Erlang A satisfies Little's law and
 * rate conservation to machine precision, the fluid model has a closed form
 * under exponential patience and is the s -> Inf limit of Erlang A, alpha(beta)
 * is the s -> Inf limit of Erlang C, and the Mt/G/Inf mean has a closed form for
 * a sinusoidal rate with exponential service.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_gig1_bnds_extremal.h"
#include "line/api/qsys/qsys_ggingi_tga.h"
#include "line/api/qsys/qsys_ggisgi_fluid.h"
#include "line/api/qsys/qsys_ggnm_diffusion.h"
#include "line/api/qsys/qsys_mgisrgi_whitt.h"
#include "line/api/qsys/qsys_mmk_qed.h"
#include "line/api/qsys/qsys_mtginf.h"

using line::qsys::MgisrgiOptions;
using line::qsys::Patience;
using line::qsys::QedCriterion;

namespace {
constexpr double LAM = 102.0, MU = 1.0, THETA = 0.5;
constexpr unsigned S = 100;
constexpr double R = 200.0;
}  // namespace

TEST_CASE("qsys_erlanga satisfies Little's law and rate conservation exactly") {
    const auto res = line::qsys::qsys_erlanga<double>(LAM, MU, THETA, S, R);
    // E[Q] = lambda_eff E[W]: the queue length comes from the birth-death chain
    // and the wait from the customer-experience recursions, so their agreement
    // is a real check that eqs. (7.10)-(7.29) were transcribed correctly.
    const double lamEff = LAM * (1.0 - res.probLoss);
    CHECK(res.meanQueueLength == doctest::Approx(lamEff * res.meanWait).epsilon(1e-12));
    // Every abandonment is a waiting customer whose exponential patience fired.
    CHECK(res.abandonRate == doctest::Approx(THETA * res.meanQueueLength).epsilon(1e-12));
    // Everything served left a busy server.
    CHECK(res.throughput ==
          doctest::Approx(MU * res.utilization * static_cast<double>(S)).epsilon(1e-12));
    CHECK(res.probServed + res.probAbandon == doctest::Approx(1.0));
    CHECK(res.exponentialPatience);
}

TEST_CASE("the three patience forms agree when the patience is exponential") {
    const auto viaRate = line::qsys::qsys_erlanga<double>(LAM, MU, THETA, S, R);
    const auto viaHazard = line::qsys::qsys_mgisrgi_whitt<double>(
        LAM, MU, S, R, Patience<double>::hazard([](const double&) { return THETA; }));
    // The ccdf form integrates the hazard, so it agrees only once the missing
    // factor lambda of the printed eq. (3.6) is restored; that is the point of
    // the divergence noted in the header.
    const auto viaCcdf = line::qsys::qsys_mgisrgi_whitt<double>(
        LAM, MU, S, R,
        Patience<double>::ccdf([](const double& t) { return std::exp(-THETA * t); }));
    CHECK(viaHazard.probAbandon == doctest::Approx(viaRate.probAbandon).epsilon(1e-12));
    CHECK(viaCcdf.probAbandon == doctest::Approx(viaRate.probAbandon).epsilon(1e-12));
    CHECK(viaCcdf.meanWaitServed == doctest::Approx(viaRate.meanWaitServed).epsilon(1e-12));
}

TEST_CASE("a patience law with zero hazard at the origin abandons far less") {
    // Erlang-2 patience of mean 2 has h(0) = 0, so nobody abandons immediately;
    // exponential patience of the same mean abandons from the first instant.
    // Whitt's point is that the mean is not what matters, the hazard at 0 is.
    const double l2 = 1.0;
    const auto e2 = line::qsys::qsys_mgisrgi_whitt<double>(
        LAM, MU, S, R,
        Patience<double>::hazard([l2](const double& t) { return l2 * (l2 * t) / (1.0 + l2 * t); }));
    const auto ex = line::qsys::qsys_erlanga<double>(LAM, MU, 0.5, S, R);
    CHECK(e2.probAbandon < ex.probAbandon);
    CHECK(e2.meanQueueLength > ex.meanQueueLength);
}

TEST_CASE("the waiting-time cdfs are proper and consistent with the split") {
    MgisrgiOptions opts;
    opts.wPoints = {0.01, 0.1, 0.5, 5.0};
    const auto res = line::qsys::qsys_erlanga<double>(LAM, MU, THETA, S, R, opts);
    for (std::size_t i = 0; i < opts.wPoints.size(); ++i) {
        CHECK(res.cdfWaitServed[i] >= 0.0);
        CHECK(res.cdfWaitServed[i] <= 1.0 + 1e-9);
        CHECK(res.cdfWaitAbandon[i] >= 0.0);
        CHECK(res.cdfWaitAbandon[i] <= 1.0 + 1e-9);
        if (i > 0) CHECK(res.cdfWait[i] >= res.cdfWait[i - 1] - 1e-12);
    }
    // Far out in the tail everyone has left, one way or the other.
    CHECK(res.cdfWait.back() == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("qsys_ggisgi_fluid matches its closed form under exponential patience") {
    // F^c(w) = 1/rho gives w = log(rho)/theta and Q = lambda(1-1/rho)/theta.
    const double rho = 1.1, mu = 1.0, theta = 0.5;
    const unsigned s = 1000;
    const double lam = rho * s * mu;
    const auto f = line::qsys::qsys_ggisgi_fluid<double>(
        lam, mu, s, [theta](const double& t) { return std::exp(-theta * t); });
    CHECK(f.regime == "overloaded");
    CHECK(f.offeredWait == doctest::Approx(std::log(rho) / theta).epsilon(1e-9));
    CHECK(f.meanQueueLength == doctest::Approx(lam * (1.0 - 1.0 / rho) / theta).epsilon(1e-8));
    CHECK(f.probAbandon == doctest::Approx(1.0 - 1.0 / rho).epsilon(1e-12));
    CHECK(f.utilization == doctest::Approx(1.0));
}

TEST_CASE("qsys_ggisgi_fluid is the many-server limit of Erlang A") {
    const double rho = 1.1, mu = 1.0, theta = 0.5;
    double prevErr = 1e9;
    for (unsigned s : {100u, 1000u, 10000u}) {
        const double lam = rho * s * mu;
        const auto f = line::qsys::qsys_ggisgi_fluid<double>(
            lam, mu, s, [theta](const double& t) { return std::exp(-theta * t); });
        const auto a = line::qsys::qsys_erlanga<double>(lam, mu, theta, s);
        const double err = std::fabs(f.meanQueueLength - a.meanQueueLength) / a.meanQueueLength;
        CHECK(err < prevErr);          // the gap closes as the system grows
        prevErr = err;
    }
    CHECK(prevErr < 1e-4);
}

TEST_CASE("the fluid abandonment probability ignores the patience law beyond rho") {
    // Corollary 3.1(i): P(ab) = 1 - 1/rho whatever F is, while the wait does
    // depend on F beyond its mean.
    const double lam = 1100.0, mu = 1.0;
    const unsigned s = 1000;
    const auto ex = line::qsys::qsys_ggisgi_fluid<double>(
        lam, mu, s, [](const double& t) { return std::exp(-t / 2.0); });
    const auto e2 = line::qsys::qsys_ggisgi_fluid<double>(
        lam, mu, s, [](const double& t) { return (1.0 + t) * std::exp(-t); });
    CHECK(ex.probAbandon == doctest::Approx(e2.probAbandon).epsilon(1e-12));
    CHECK(e2.offeredWait > ex.offeredWait);
}

TEST_CASE("the underloaded fluid model is the infinite-server one") {
    const auto u = line::qsys::qsys_ggisgi_fluid<double>(
        80.0, 1.0, 100u, [](const double& t) { return std::exp(-t / 2.0); });
    CHECK(u.regime == "underloaded");
    CHECK(u.meanQueueLength == doctest::Approx(0.0));
    CHECK(u.meanNumberInService == doctest::Approx(80.0));
    CHECK(u.utilization == doctest::Approx(0.8));
    CHECK(u.probAbandon == doctest::Approx(0.0));
}

TEST_CASE("alpha(beta) is the many-server limit of Erlang C") {
    for (double beta : {0.1, 0.5, 1.0, 2.0}) {
        const double target = line::qsys::qsys_mmk_qed_alpha<double>(beta);
        double prevErr = 1e9;
        for (unsigned s : {100u, 1000u, 10000u}) {
            const double lam = (1.0 - beta / std::sqrt(static_cast<double>(s))) * s;
            const double c = line::qsys::qsys_mmk_qed_erlangc<double>(s, lam, 1.0);
            const double err = std::fabs(c - target);
            CHECK(err < prevErr);
            prevErr = err;
        }
    }
}

TEST_CASE("alpha(beta) decreases and stays finite where the normal density underflows") {
    CHECK(line::qsys::qsys_mmk_qed_alpha<double>(0.0) == doctest::Approx(1.0));
    double prev = 1.0;
    for (double beta : {0.5, 1.0, 2.0, 5.0, 20.0, 40.0}) {
        const double a = line::qsys::qsys_mmk_qed_alpha<double>(beta);
        CHECK(a <= prev);
        CHECK(a >= 0.0);
        CHECK(std::isfinite(a));
        prev = a;
    }
}

TEST_CASE("square-root staffing lands where the exact Erlang C rule does") {
    for (double target : {0.5, 0.2, 0.05}) {
        const auto approx = line::qsys::qsys_mmk_qed_staffing<double>(1000.0, 1.0, target);
        const auto exact = line::qsys::qsys_mmk_qed_staffing<double>(
            1000.0, 1.0, target, QedCriterion::Delay, 0.0, 0.0, true);
        CHECK(approx.numServers == exact.numServers);
        CHECK(line::qsys::qsys_mmk_qed_erlangc<double>(exact.numServers, 1000.0, 1.0) <= target);
    }
}

TEST_CASE("the service-level criterion meets its deadline") {
    const auto sl = line::qsys::qsys_mmk_qed_staffing<double>(1000.0, 1.0, 0.0,
                                                              QedCriterion::ServiceLevel, 0.02, 0.8);
    CHECK(sl.serviceLevel >= 0.8);
    CHECK(sl.numServers > 1000u);
}

TEST_CASE("qsys_mtginf matches the sinusoidal closed form") {
    // lambda(t) = a + b sin(g t), S ~ Exp(mu):
    //   m(t) = a/mu + b (mu sin(gt) - g cos(gt))/(mu^2+g^2)
    const double a = 10.0, b = 5.0, g = 1.0, mu = 2.0;
    std::function<double(const double&)> lam = [&](const double& t) {
        return a + b * std::sin(g * t);
    };
    std::function<double(const double&)> gc = [&](const double& x) { return std::exp(-mu * x); };
    std::function<double(const double&)> pdf = [&](const double& x) {
        return mu * std::exp(-mu * x);
    };
    std::vector<double> ts;
    for (int i = 0; i < 5; ++i) ts.push_back(i * 2.0 * M_PI / 4.0);
    const auto r = line::qsys::qsys_mtginf<double>(lam, gc, 1.0 / mu, ts,
                                                   -std::numeric_limits<double>::infinity(),
                                                   2.0 / (mu * mu), pdf);
    for (std::size_t i = 0; i < ts.size(); ++i) {
        const double closed =
            a / mu + b * (mu * std::sin(g * ts[i]) - g * std::cos(g * ts[i])) / (mu * mu + g * g);
        CHECK(r.meanNumber[i] == doctest::Approx(closed).epsilon(1e-8));
        CHECK(r.varNumber[i] == doctest::Approx(r.meanNumber[i]));   // Poisson
        // delta(t) = E[lambda(t-S)] has its own closed form here.
        const double dclosed =
            a + b * (mu * mu * std::sin(g * ts[i]) - mu * g * std::cos(g * ts[i])) / (mu * mu + g * g);
        CHECK(r.departureRate[i] == doctest::Approx(dclosed).epsilon(1e-8));
    }
    CHECK(r.hasLag);
    CHECK(r.meanLag == doctest::Approx(1.0 / mu));      // E[S_e] = E[S] for exponential service
}

TEST_CASE("the mtginf departure rate agrees whether it is integrated or differenced") {
    const double a = 10.0, b = 5.0, g = 1.0, mu = 2.0;
    std::function<double(const double&)> lam = [&](const double& t) {
        return a + b * std::sin(g * t);
    };
    std::function<double(const double&)> gc = [&](const double& x) { return std::exp(-mu * x); };
    std::function<double(const double&)> pdf = [&](const double& x) {
        return mu * std::exp(-mu * x);
    };
    const std::vector<double> ts = {0.0, 1.0, 2.0, 3.0};
    const auto withPdf = line::qsys::qsys_mtginf<double>(
        lam, gc, 1.0 / mu, ts, -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN(), pdf);
    const auto byBalance = line::qsys::qsys_mtginf<double>(lam, gc, 1.0 / mu, ts);
    for (std::size_t i = 0; i < ts.size(); ++i)
        CHECK(withPdf.departureRate[i] == doctest::Approx(byBalance.departureRate[i]).epsilon(1e-6));
}

TEST_CASE("an mtginf queue started empty fills from zero") {
    const double a = 10.0, b = 5.0, mu = 2.0;
    std::function<double(const double&)> lam = [&](const double& t) {
        return a + b * std::sin(t);
    };
    std::function<double(const double&)> gc = [&](const double& x) { return std::exp(-mu * x); };
    const std::vector<double> ts = {0.0, 0.1, 0.5, 1.0};
    const auto r = line::qsys::qsys_mtginf<double>(lam, gc, 1.0 / mu, ts, 0.0);
    CHECK(r.meanNumber[0] == doctest::Approx(0.0));
    for (std::size_t i = 1; i < ts.size(); ++i) CHECK(r.meanNumber[i] > r.meanNumber[i - 1]);
}

TEST_CASE("the extremal GI/GI/1 bounds reproduce Table 1 of Chen and Whitt (2020)") {
    // ca^2 = cs^2 = 4 with E[U] = 1, i.e. lambda = 1 and mu = 1/rho.
    struct Row {
        double rho, lb, hta, ub, ubClosed, delta, daley, kingman;
    };
    const Row rows[] = {{0.10, 0.000, 0.044, 0.422, 0.422, 0.000, 0.444, 2.244},
                        {0.20, 0.000, 0.200, 0.904, 0.906, 0.007, 1.000, 2.600},
                        {0.50, 0.750, 2.000, 3.470, 3.510, 0.203, 4.000, 5.000},
                        {0.80, 6.000, 12.800, 14.917, 15.017, 0.629, 16.000, 16.400},
                        {0.90, 15.750, 32.400, 34.721, 34.843, 0.807, 36.000, 36.200}};
    for (const Row& r : rows) {
        const auto b = line::qsys::qsys_gig1_bnds_extremal<double>(1.0, 1.0 / r.rho, 2.0, 2.0,
                                                                   4000, 4000);
        CHECK(b.lowerBound == doctest::Approx(r.lb).epsilon(1e-3));
        CHECK(b.heavyTraffic == doctest::Approx(r.hta).epsilon(1e-3));
        CHECK(b.upperBound == doctest::Approx(r.ub).epsilon(1e-3));
        CHECK(b.upperBoundClosed == doctest::Approx(r.ubClosed).epsilon(1e-3));
        CHECK(b.delta == doctest::Approx(r.delta).epsilon(1e-2));
        CHECK(b.upperBoundDaley == doctest::Approx(r.daley).epsilon(1e-3));
        CHECK(b.upperBoundKingman == doctest::Approx(r.kingman).epsilon(1e-3));
    }
}

TEST_CASE("the extremal bounds are ordered and scale with the time unit") {
    const auto b = line::qsys::qsys_gig1_bnds_extremal<double>(0.7, 1.0, 1.5, 1.5, 2000, 2000);
    CHECK(b.lowerBound <= b.heavyTraffic);
    CHECK(b.heavyTraffic <= b.upperBound);
    CHECK(b.upperBound <= b.upperBoundClosed);
    CHECK(b.upperBoundClosed <= b.upperBoundDaley);
    CHECK(b.upperBoundDaley <= b.upperBoundKingman);
    CHECK(b.relativeWidth > 0.0);
    CHECK(b.relativeWidth < 1.0);
    // Doubling both rates halves every waiting time.
    const auto slow = line::qsys::qsys_gig1_bnds_extremal<double>(1.0, 2.0, 2.0, 2.0, 1, 1, true);
    const auto fast = line::qsys::qsys_gig1_bnds_extremal<double>(2.0, 4.0, 2.0, 2.0, 1, 1, true);
    CHECK(slow.upperBoundClosed == doctest::Approx(2.0 * fast.upperBoundClosed).epsilon(1e-12));
}

TEST_CASE("the M/M/1 mean wait lies inside the two-moment interval") {
    // At ca^2 = cs^2 = 1 the exact M/M/1 wait rho^2/((1-rho)lambda) is one of the
    // values two moments admit, so it must sit between the bounds. It is NOT the
    // upper end: Kingman is 1.7 against the exact 0.9 at rho = 0.6, which is the
    // whole reason the tight bound was worth computing.
    const double rho = 0.6, lambda = 1.0;
    const auto b = line::qsys::qsys_gig1_bnds_extremal<double>(lambda, lambda / rho, 1.0, 1.0,
                                                               3000, 3000);
    const double mm1 = rho * rho / ((1.0 - rho) * lambda);
    CHECK(b.lowerBound <= mm1);
    CHECK(mm1 <= b.upperBound);
    CHECK(b.upperBound < b.upperBoundKingman);
}

TEST_CASE("the G/GI/n/m diffusion reduces to Halfin-Whitt with an unbounded room") {
    for (double beta : {0.5, 1.0, 2.0}) {
        const unsigned s = 1000;
        const double lam = (1.0 - beta / std::sqrt(static_cast<double>(s))) * s;
        const auto d = line::qsys::qsys_ggnm_diffusion<double>(
            lam, 1.0, s, std::numeric_limits<double>::infinity(), 1.0, 1.0);
        CHECK(d.probDelay == doctest::Approx(line::qsys::qsys_mmk_qed_alpha<double>(beta))
                                 .epsilon(1e-12));
        CHECK(d.peakedness == doctest::Approx(1.0));
        CHECK(d.probBlock == doctest::Approx(0.0));
    }
}

TEST_CASE("the diffusion peakedness ranks the service laws as the theory says") {
    // omega_G is 1 for deterministic service, 1/2 for exponential, and falls as
    // service gets more variable; z = 1 + (ca^2-1)omega then orders the delay
    // probabilities when ca^2 > 1.
    const double lam = 950.0;
    std::function<double(const double&)> det = [](const double& x) { return x < 1.0 ? 1.0 : 0.0; };
    std::function<double(const double&)> erl2 = [](const double& x) {
        return (1.0 + 2.0 * x) * std::exp(-2.0 * x);
    };
    const auto d = line::qsys::qsys_ggnm_diffusion<double>(
        lam, 1.0, 1000u, std::numeric_limits<double>::infinity(), 2.0, 0.0, det);
    const auto e = line::qsys::qsys_ggnm_diffusion<double>(
        lam, 1.0, 1000u, std::numeric_limits<double>::infinity(), 2.0, 1.0 / std::sqrt(2.0), erl2);
    const auto m = line::qsys::qsys_ggnm_diffusion<double>(
        lam, 1.0, 1000u, std::numeric_limits<double>::infinity(), 2.0, 1.0);
    CHECK(d.peakednessWeight == doctest::Approx(1.0).epsilon(1e-3));
    CHECK(m.peakednessWeight == doctest::Approx(0.5));
    CHECK(e.peakednessWeight == doctest::Approx(0.625).epsilon(1e-6));
    CHECK(d.probDelay > e.probDelay);
    CHECK(e.probDelay > m.probDelay);
}

TEST_CASE("the heavily-loaded Gaussian approximation beats the fluid model") {
    // Against the exact Erlang A the truncated Gaussian must land between the
    // fluid answer and the truth, and closer to the truth.
    const double th = 0.5, mu = 1.0, rho = 1.05;
    std::function<double(const double&)> fc = [th](const double& x) { return std::exp(-th * x); };
    std::function<double(const double&)> fp = [th](const double& x) {
        return th * std::exp(-th * x);
    };
    for (unsigned n : {50u, 100u, 400u}) {
        const double lam = rho * n * mu;
        const auto t = line::qsys::qsys_ggingi_tga<double>(lam, mu, n, 1.0, 1.0, fc, fp);
        const auto e = line::qsys::qsys_erlanga<double>(lam, mu, th, n);
        const auto f = line::qsys::qsys_ggisgi_fluid<double>(lam, mu, n, fc);
        const double errTga = std::fabs(t.meanQueueLength - e.meanQueueLength);
        const double errFluid = std::fabs(f.meanQueueLength - e.meanQueueLength);
        CHECK(errTga < errFluid);
        CHECK(t.meanQueueLength > f.meanQueueLength);        // the fluid model undershoots
        CHECK(t.regime == "overloaded");
    }
}

TEST_CASE("the heavily-loaded approximation converges to the exact Erlang A") {
    const double th = 0.5, mu = 1.0, rho = 1.2;
    std::function<double(const double&)> fc = [th](const double& x) { return std::exp(-th * x); };
    std::function<double(const double&)> fp = [th](const double& x) {
        return th * std::exp(-th * x);
    };
    const unsigned n = 1000;
    const double lam = rho * n * mu;
    const auto t = line::qsys::qsys_ggingi_tga<double>(lam, mu, n, 1.0, 1.0, fc, fp);
    const auto e = line::qsys::qsys_erlanga<double>(lam, mu, th, n);
    CHECK(t.meanQueueLength == doctest::Approx(e.meanQueueLength).epsilon(1e-4));
    CHECK(t.probAbandon == doctest::Approx(e.probAbandon).epsilon(2e-3));
    CHECK(t.meanWait == doctest::Approx(e.meanWaitServed).epsilon(2e-3));
}

TEST_CASE("an underloaded heavily-loaded model reports no queue") {
    std::function<double(const double&)> fc = [](const double& x) { return std::exp(-0.5 * x); };
    const auto u = line::qsys::qsys_ggingi_tga<double>(80.0, 1.0, 100u, 1.5, 1.0, fc);
    CHECK(u.regime == "underloaded");
    CHECK(u.meanQueueLength == doctest::Approx(0.0));
    CHECK(u.meanNumber == doctest::Approx(80.0));
    CHECK(u.sigmaX > 0.0);
}
