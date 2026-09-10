/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Cohen (1979) Sect. 9 workload distribution of a load-dependent PS station
 * with blocking (line/api/qsys/qsys_ldps_workload.h).
 *
 * Oracles, in the order the task prescribes.
 *  (a) Closed forms. The state law is the truncated load-dependent
 *      birth-death law p_h propto rho^h/prod_k alpha(k), which for the
 *      instance used is the exact rational vector (125, 75, 30, 10, 3)/243 and
 *      is asserted as such. For an EXPONENTIAL required service time the
 *      equilibrium residual is exponential again by memorylessness, so the
 *      workload given h in system is Erlang(h, beta) and the whole CDF is the
 *      Erlang mixture sum_h p_h Erlang_h(t) -- a closed form for the very
 *      quantity the grid computes.
 *  (b) Invariants: F is non-decreasing, F(0) equals the atom p_0 exactly, F
 *      approaches 1 at the top of the support, and p is a probability vector.
 *      Plus the convergence rate: the reference states second order for an
 *      absolutely continuous B and first order when B has an atom, and both
 *      regimes are measured here rather than assumed (4.0x per grid doubling
 *      for Exp, about 2x for Det).
 *  (c) MATLAB, digit for digit at 1e-15 relative on both the automatic and the
 *      user-supplied grid.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_ldps_workload.h"

using line::Real50;
using line::qsys::LdpsWorkloadResult;
using line::qsys::WorkloadServiceLaw;
using line::qsys::qsys_ldps_workload;

namespace {

/** Exponential of mean 1/rate. */
WorkloadServiceLaw<double> exp_law(double rate) {
    WorkloadServiceLaw<double> B;
    B.mean = 1.0 / rate;
    B.scv = 1.0;
    B.cdf = [rate](const double& x) { return 1.0 - std::exp(-rate * x); };
    return B;
}

/** Constant service time d. */
WorkloadServiceLaw<double> det_law(double d) {
    WorkloadServiceLaw<double> B;
    B.mean = d;
    B.scv = 0.0;
    B.cdf = [d](const double& x) { return x < d ? 0.0 : 1.0; };
    return B;
}

/** Erlang(h, scale) CDF, h >= 1: 1 - sum_{k<h} e^{-x/s}(x/s)^k/k!. */
double erlang_cdf(double x, unsigned h, double s) {
    const double z = x / s;
    double term = std::exp(-z), acc = term;
    for (unsigned k = 1; k < h; ++k) {
        term *= z / static_cast<double>(k);
        acc += term;
    }
    return 1.0 - acc;
}

/** The Erlang-mixture oracle for exponential required service. */
double erlang_mixture(double x, const std::vector<double>& p, double beta) {
    double f = p[0];
    for (std::size_t h = 1; h < p.size(); ++h) f += p[h] * erlang_cdf(x, static_cast<unsigned>(h), beta);
    return f;
}

const std::vector<double> kAlpha = {1.0, 1.5, 1.8, 2.0};
const double kLambda = 1.2;
const std::size_t kN = 4;

}  // namespace

// ---------------------------------------------------------------------------
// (a) closed forms
// ---------------------------------------------------------------------------

TEST_CASE("the state law is the exact truncated load-dependent birth-death law") {
    // rho = lambda beta = 0.6, alpha = (1, 1.5, 1.8, 2). The unnormalized
    // weights are 1, 0.6, 0.24, 0.08, 0.024 summing to 1.944 = 243/125, so
    // p = (125, 75, 30, 10, 3)/243.
    const auto r = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN);
    const double den = 243.0;
    const std::vector<double> exact = {125.0 / den, 75.0 / den, 30.0 / den, 10.0 / den,
                                       3.0 / den};
    REQUIRE(r.p.size() == kN + 1);
    double mass = 0.0;
    for (std::size_t h = 0; h <= kN; ++h) {
        CHECK(r.p[h] == doctest::Approx(exact[h]).epsilon(1e-15));
        mass += r.p[h];
    }
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-15));
    // The atom at the origin: F(0) = p_0 exactly, no quadrature involved.
    CHECK(r.F[0] == doctest::Approx(exact[0]).epsilon(1e-15));
    CHECK(r.t[0] == 0.0);
}

TEST_CASE("exponential service reproduces the Erlang mixture") {
    const std::vector<double> tq = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    const auto r = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN, tq,
                                      static_cast<std::size_t>(2001));
    const double beta = 0.5;
    for (std::size_t j = 0; j < tq.size(); ++j) {
        INFO("t = ", tq[j]);
        CHECK(r.F[j] == doctest::Approx(erlang_mixture(tq[j], r.p, beta)).epsilon(2e-6));
    }
    // Exact at the origin, where the mixture is the atom alone.
    CHECK(r.F[0] == doctest::Approx(r.p[0]).epsilon(1e-15));
}

// ---------------------------------------------------------------------------
// (b) invariants and the two convergence regimes
// ---------------------------------------------------------------------------

TEST_CASE("the workload CDF is a CDF") {
    for (int which = 0; which < 2; ++which) {
        const WorkloadServiceLaw<double> B = which == 0 ? exp_law(2.0) : det_law(0.5);
        INFO("service law index ", which);
        const auto r = qsys_ldps_workload(kLambda, B, kAlpha, kN);
        REQUIRE(r.F.size() == r.t.size());
        for (std::size_t i = 1; i < r.F.size(); ++i) {
            CHECK(r.t[i] > r.t[i - 1]);
            CHECK(r.F[i] >= r.F[i - 1] - 1e-15);
        }
        CHECK(r.F.front() == doctest::Approx(r.p[0]).epsilon(1e-15));
        CHECK(r.F.back() == doctest::Approx(1.0).epsilon(3e-3));
    }
}

TEST_CASE("the grid error is second order for a continuous B and first for an atom") {
    const std::vector<double> tq = {0.25, 0.5, 1.0, 2.0};
    SUBCASE("exponential: 4x per grid doubling") {
        std::vector<double> err;
        for (std::size_t ng : {501u, 1001u, 2001u, 4001u}) {
            const auto r = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN, tq, ng);
            double e = 0.0;
            for (std::size_t j = 0; j < tq.size(); ++j)
                e = std::max(e, std::fabs(r.F[j] - erlang_mixture(tq[j], r.p, 0.5)));
            err.push_back(e);
        }
        for (std::size_t k = 1; k < err.size(); ++k) {
            INFO("ratio = ", err[k - 1] / err[k]);
            CHECK(err[k - 1] / err[k] > 3.7);
            CHECK(err[k - 1] / err[k] < 4.3);
        }
    }
    SUBCASE("deterministic: the mass defect halves, not quarters") {
        // Four residuals of a constant service, each uniform on [0, 0.5], so
        // the workload has bounded support 2 and the CDF plateaus there. The
        // distance of the plateau from 1 is the discretization defect.
        std::vector<double> defect;
        for (std::size_t ng : {501u, 1001u, 2001u, 4001u}) {
            const auto r = qsys_ldps_workload(kLambda, det_law(0.5), kAlpha, kN,
                                              std::vector<double>{2.5}, ng);
            defect.push_back(std::fabs(1.0 - r.F[0]));
        }
        for (std::size_t k = 1; k < defect.size(); ++k) {
            INFO("ratio = ", defect[k - 1] / defect[k]);
            CHECK(defect[k - 1] / defect[k] > 1.8);
            CHECK(defect[k - 1] / defect[k] < 2.4);
        }
    }
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("qsys_ldps_workload agrees with MATLAB digit for digit") {
    SUBCASE("exponential service, automatic grid") {
        const auto r = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN);
        CHECK(r.t.size() == 2001u);
        CHECK(r.t.back() == doctest::Approx(10.0).epsilon(1e-15));
        CHECK(r.F.back() == doctest::Approx(1.0000015818111168).epsilon(1e-14));
        CHECK(r.F.front() == doctest::Approx(0.51440329218106995).epsilon(1e-15));
    }
    SUBCASE("exponential service, user grid") {
        const std::vector<double> tq = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
        const auto r = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN, tq,
                                          static_cast<std::size_t>(2001));
        const std::vector<double> ref = {0.51440329218106995, 0.64759358664087197,
                                         0.74566404416178456, 0.86967785548238419,
                                         0.96789206403881989, 0.99843555599153122,
                                         0.99999896116835618};
        for (std::size_t j = 0; j < tq.size(); ++j)
            CHECK(r.F[j] == doctest::Approx(ref[j]).epsilon(1e-14));
    }
    SUBCASE("deterministic service") {
        const auto r = qsys_ldps_workload(kLambda, det_law(0.5), kAlpha, kN);
        CHECK(r.t.back() == doctest::Approx(5.0).epsilon(1e-15));
        CHECK(r.F.back() == doctest::Approx(0.99817991216563784).epsilon(1e-14));
        const std::vector<double> tq = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
        const auto q = qsys_ldps_workload(kLambda, det_law(0.5), kAlpha, kN, tq,
                                          static_cast<std::size_t>(2001));
        const std::vector<double> ref = {0.51440329218106995, 0.68504715843739916,
                                         0.89090853793606584, 0.98446466636273244,
                                         0.99708872269958837, 0.99708872269958837,
                                         0.99708872269958837};
        for (std::size_t j = 0; j < tq.size(); ++j)
            CHECK(q.F[j] == doctest::Approx(ref[j]).epsilon(1e-14));
    }
}

TEST_CASE("qsys_ldps_workload rejects an instance outside the assumed model") {
    const WorkloadServiceLaw<double> B = exp_law(2.0);
    CHECK_THROWS_AS(qsys_ldps_workload(-1.0, B, kAlpha, kN), line::InputError);
    CHECK_THROWS_AS(qsys_ldps_workload(kLambda, B, kAlpha, static_cast<std::size_t>(0)),
                    line::InputError);
    // alpha too short for N: every busy level needs a rate.
    CHECK_THROWS_AS(qsys_ldps_workload(kLambda, B, std::vector<double>{1.0, 1.5}, kN),
                    line::InputError);
    // A non-positive rate scaling would stall the stage.
    CHECK_THROWS_AS(
        qsys_ldps_workload(kLambda, B, std::vector<double>{1.0, 0.0, 1.8, 2.0}, kN),
        line::InputError);
    WorkloadServiceLaw<double> bad = B;
    bad.mean = 0.0;
    CHECK_THROWS_AS(qsys_ldps_workload(kLambda, bad, kAlpha, kN), line::InputError);
}

TEST_CASE("qsys_ldps_workload at Real50 reproduces its own double result") {
    WorkloadServiceLaw<Real50> B;
    B.mean = Real50(1) / Real50(2);
    B.scv = Real50(1);
    B.cdf = [](const Real50& x) { return Real50(1) - exp(Real50(-2) * x); };
    std::vector<Real50> al;
    for (double v : kAlpha) al.push_back(Real50(v));
    const std::vector<Real50> tq = {Real50(1), Real50(2)};
    const auto rq = qsys_ldps_workload(Real50(kLambda), B, al, kN, tq,
                                       static_cast<std::size_t>(501));
    const auto rd = qsys_ldps_workload(kLambda, exp_law(2.0), kAlpha, kN,
                                       std::vector<double>{1.0, 2.0},
                                       static_cast<std::size_t>(501));
    for (std::size_t j = 0; j < 2; ++j)
        CHECK(static_cast<double>(rq.F[j]) == doctest::Approx(rd.F[j]).epsilon(1e-13));
}
