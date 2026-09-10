/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Exhaustive and gated polling (Takagi 1988, eqs. 15 and 20) and the M/M/c/c
 * retrial fixed point.
 *
 * Oracles, in the order the task prescribes.
 *  (a) Closed forms the models collapse to. A one-queue polling system is the
 *      M/G/1 queue with multiple vacations, W = lambda b2/(2(1-rho)) +
 *      E[R^2]/(2 E[R]) under exhaustive service and that plus rho R/(1-rho)
 *      under gated service. Both are asserted as EXACT RATIONALS, digit for
 *      digit (45/112 and 53/112 on the instance used), which is the sharpest
 *      form of the check and is only available because the port stays in the
 *      field. The retrial fixed point collapses to r = 1/6, B = 1/4 on
 *      lambda = 1/2, mu = 2, c = 1, where the fixed-point equation is a
 *      quadratic that can be solved by hand.
 *  (b) The invariant. Boxma and Groenendijk's pseudo-conservation law couples
 *      the two disciplines through one identity and is asserted on every
 *      instance, at 1e-12 relative in double.
 *  (c) MATLAB, at 1e-13 relative.
 *
 * No simulator is used: LINE's LDES polling server parks at the last visited
 * queue rather than roving, so it sits a few percent BELOW Takagi by
 * construction and cannot adjudicate an absolute value.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/polling/polling_qsys_exhaustive.h"
#include "line/api/qsys/qsys_mmcc_retrial_fp.h"

using line::Rational;
using line::Real50;
using line::polling::PollingMoments;
using line::polling::polling_qsys_exhaustive;
using line::polling::polling_qsys_gated;

namespace {

/** Build a moment set from double literals. */
template <class T>
PollingMoments<T> mk(const std::vector<double>& lambda, const std::vector<double>& b,
                     const std::vector<double>& b2, const std::vector<double>& r,
                     const std::vector<double>& delta2) {
    PollingMoments<T> m;
    for (double v : lambda) m.lambda.push_back(line::num_traits<T>::from_double(v));
    for (double v : b) m.b.push_back(line::num_traits<T>::from_double(v));
    for (double v : b2) m.b2.push_back(line::num_traits<T>::from_double(v));
    for (double v : r) m.r.push_back(line::num_traits<T>::from_double(v));
    for (double v : delta2) m.delta2.push_back(line::num_traits<T>::from_double(v));
    return m;
}

/**
 * Right-hand side of the pseudo-conservation law of Boxma and Groenendijk
 * (1987):
 *   sum_i rho_i W_i = rho/(2(1-rho)) sum_i lambda_i b2_i
 *                   + rho (delta2tot + R^2)/(2R)
 *                   + R (rho^2 - sum_i rho_i^2)/(2(1-rho))
 *                   + sum_i E[M_i],
 * with E[M_i] the mean work at queue i when it is polled: zero under
 * exhaustive service, rho_i^2 R/(1-rho) under gated service.
 */
double pcl_rhs(const PollingMoments<double>& m, bool gated) {
    const std::size_t n = m.size();
    double rho = 0.0, R = 0.0, D2 = 0.0, sumLb2 = 0.0, sumRho2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double ri = m.lambda[i] * m.b[i];
        rho += ri;
        sumRho2 += ri * ri;
        R += m.r[i];
        D2 += m.delta2[i];
        sumLb2 += m.lambda[i] * m.b2[i];
    }
    double v = rho / (2.0 * (1.0 - rho)) * sumLb2 + rho * (D2 + R * R) / (2.0 * R) +
               R * (rho * rho - sumRho2) / (2.0 * (1.0 - rho));
    if (gated) v += sumRho2 * R / (1.0 - rho);
    return v;
}

double pcl_lhs(const PollingMoments<double>& m, const std::vector<double>& W) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.size(); ++i) s += m.lambda[i] * m.b[i] * W[i];
    return s;
}

/** Symmetric two-queue system, exponential everything (MATLAB reference case). */
PollingMoments<double> case2q() {
    return mk<double>({0.6, 0.2}, {1.0, 1.0}, {2.0, 2.0}, {1.0, 1.0}, {1.0, 1.0});
}

/** Three queues, Erlang and exponential switchovers, rho = 0.395. */
PollingMoments<double> case3q() {
    return mk<double>({0.3, 0.5, 0.4}, {0.4, 0.15, 0.5}, {0.32, 0.03375, 0.5},
                      {0.2, 0.15, 0.1}, {0.02, 0.0225, 0.0033333333333333305});
}

/** Four queues, rho = 0.4, all five moment vectors asymmetric. */
PollingMoments<double> case4q() {
    return mk<double>({0.4, 0.25, 1.0 / 3.0, 0.2}, {0.5, 0.2, 0.3, 0.25},
                      {0.5, 0.06, 0.18, 0.08333333333333334}, {0.15, 0.05, 0.25, 0.1},
                      {0.01125, 0.000625, 0.0625, 0.005});
}

}  // namespace

// ---------------------------------------------------------------------------
// (a) closed-form collapses, exact in Rational
// ---------------------------------------------------------------------------

TEST_CASE("one-queue polling is the M/G/1 queue with multiple vacations, exactly") {
    // lambda = 1/2, b = 1/4, b2 = 3/32, R ~ Erlang-2 with mean 1/2 so
    // delta2 = 1/8. rho = 1/8.
    //   W_exhaustive = lambda b2/(2(1-rho)) + E[R^2]/(2 E[R])
    //                = (3/64)/(7/4) + (3/8) = 3/112 + 3/8 = 45/112
    //   W_gated      = W_exhaustive + rho R/(1-rho) = 45/112 + 1/14 = 53/112
    PollingMoments<Rational> m;
    m.lambda.push_back(Rational(1, 2));
    m.b.push_back(Rational(1, 4));
    m.b2.push_back(Rational(3, 32));
    m.r.push_back(Rational(1, 2));
    m.delta2.push_back(Rational(1, 8));

    const std::vector<Rational> We = polling_qsys_exhaustive(m);
    const std::vector<Rational> Wg = polling_qsys_gated(m);
    CHECK(We.size() == 1u);
    CHECK(We[0] == Rational(45, 112));
    CHECK(Wg[0] == Rational(53, 112));
    // Gated exceeds exhaustive by exactly rho R/(1 - rho) = 1/14.
    CHECK(Rational(Wg[0] - We[0]) == Rational(1, 14));
    // MATLAB returns 0.4017857142857143 and 0.4732142857142857 for the same
    // instance, which are these rationals rounded.
    CHECK(static_cast<double>(We[0]) == doctest::Approx(0.4017857142857143).epsilon(1e-15));
    CHECK(static_cast<double>(Wg[0]) == doctest::Approx(0.4732142857142857).epsilon(1e-15));
}

TEST_CASE("one-queue polling stays exact for an exponential switchover too") {
    // Exponential switchover of mean 1 has delta2 = 1, so E[R^2]/(2R) = 1 and
    // W_exhaustive = lambda b2/(2(1-rho)) + 1. With lambda = 1/2, b = 1,
    // b2 = 2, rho = 1/2: W = 1/2/(2*1/2)*2 ... = 1 + 1 = 2 exactly.
    PollingMoments<Rational> m;
    m.lambda.push_back(Rational(1, 2));
    m.b.push_back(Rational(1));
    m.b2.push_back(Rational(2));
    m.r.push_back(Rational(1));
    m.delta2.push_back(Rational(1));
    CHECK(polling_qsys_exhaustive(m)[0] == Rational(2));
    // Gated adds rho R/(1-rho) = 1.
    CHECK(polling_qsys_gated(m)[0] == Rational(3));
}

// ---------------------------------------------------------------------------
// (b) the pseudo-conservation law, on every instance
// ---------------------------------------------------------------------------

TEST_CASE("exhaustive and gated both satisfy the pseudo-conservation law") {
    for (const PollingMoments<double>& m : {case2q(), case3q(), case4q()}) {
        const std::vector<double> We = polling_qsys_exhaustive(m);
        const std::vector<double> Wg = polling_qsys_gated(m);
        INFO("n = ", m.size());
        CHECK(pcl_lhs(m, We) == doctest::Approx(pcl_rhs(m, false)).epsilon(1e-12));
        CHECK(pcl_lhs(m, Wg) == doctest::Approx(pcl_rhs(m, true)).epsilon(1e-12));
        // Gated carries strictly more work at a polling instant, so its
        // load-weighted mean waiting time is strictly larger.
        CHECK(pcl_lhs(m, Wg) > pcl_lhs(m, We));
        for (std::size_t i = 0; i < m.size(); ++i) {
            CHECK(We[i] > 0.0);
            CHECK(Wg[i] > 0.0);
        }
    }
}

TEST_CASE("the pseudo-conservation law holds exactly in rational arithmetic") {
    // Two queues, everything a small fraction, so the whole computation is a
    // finite sequence of exact field operations and the identity holds with a
    // zero residual rather than to a tolerance.
    PollingMoments<Rational> m;
    m.lambda = {Rational(3, 5), Rational(1, 5)};
    m.b = {Rational(1), Rational(1)};
    m.b2 = {Rational(2), Rational(2)};
    m.r = {Rational(1), Rational(1)};
    m.delta2 = {Rational(1), Rational(1)};
    const std::vector<Rational> We = polling_qsys_exhaustive(m);
    const std::vector<Rational> Wg = polling_qsys_gated(m);

    const Rational rho(4, 5), R(2), D2(2);
    Rational sumLb2(0), sumRho2(0), lhsE(0), lhsG(0);
    for (std::size_t i = 0; i < 2; ++i) {
        const Rational ri = m.lambda[i] * m.b[i];
        sumLb2 += m.lambda[i] * m.b2[i];
        sumRho2 += ri * ri;
        lhsE += ri * We[i];
        lhsG += ri * Wg[i];
    }
    const Rational one(1), two(2);
    const Rational base = rho / (two * (one - rho)) * sumLb2 + rho * (D2 + R * R) / (two * R) +
                          R * (rho * rho - sumRho2) / (two * (one - rho));
    CHECK(lhsE == base);
    CHECK(lhsG == Rational(base + sumRho2 * R / (one - rho)));
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("exhaustive and gated agree with MATLAB") {
    SUBCASE("two queues, lambda = (0.6, 0.2), exponential service and switchover") {
        const PollingMoments<double> m = case2q();
        const std::vector<double> We = polling_qsys_exhaustive(m);
        const std::vector<double> Wg = polling_qsys_gated(m);
        CHECK(We[0] == doctest::Approx(5.4999999999999991).epsilon(1e-13));
        CHECK(We[1] == doctest::Approx(11.499999999999998).epsilon(1e-13));
        CHECK(Wg[0] == doctest::Approx(12.77005347593583).epsilon(1e-13));
        CHECK(Wg[1] == doctest::Approx(9.6898395721925148).epsilon(1e-13));
    }
    SUBCASE("three queues, rho = 0.395") {
        const PollingMoments<double> m = case3q();
        const std::vector<double> We = polling_qsys_exhaustive(m);
        const std::vector<double> Wg = polling_qsys_gated(m);
        CHECK(We[0] == doctest::Approx(0.64823434675384695).epsilon(1e-13));
        CHECK(We[1] == doctest::Approx(0.70224422154401744).epsilon(1e-13));
        CHECK(We[2] == doctest::Approx(0.58186831009305007).epsilon(1e-13));
        CHECK(Wg[0] == doctest::Approx(0.71737485594634465).epsilon(1e-13));
        CHECK(Wg[1] == doctest::Approx(0.68269503048730473).epsilon(1e-13));
        CHECK(Wg[2] == doctest::Approx(0.77094842229819904).epsilon(1e-13));
    }
    SUBCASE("four queues, rho = 0.4, fully asymmetric") {
        const PollingMoments<double> m = case4q();
        const std::vector<double> We = polling_qsys_exhaustive(m);
        const std::vector<double> Wg = polling_qsys_gated(m);
        CHECK(We[0] == doctest::Approx(0.64675916438288406).epsilon(1e-13));
        CHECK(We[1] == doctest::Approx(0.78036785087517624).epsilon(1e-13));
        CHECK(We[2] == doctest::Approx(0.76336362918724476).epsilon(1e-13));
        CHECK(We[3] == doctest::Approx(0.7900854049359689).epsilon(1e-13));
        CHECK(Wg[0] == doctest::Approx(0.88270902621576974).epsilon(1e-13));
        CHECK(Wg[1] == doctest::Approx(0.76637505320399957).epsilon(1e-13));
        CHECK(Wg[2] == doctest::Approx(0.81601889769691738).epsilon(1e-13));
        CHECK(Wg[3] == doctest::Approx(0.76330155158959112).epsilon(1e-13));
    }
}

TEST_CASE("polling rejects an unstable or degenerate system") {
    // rho >= 1: the reference silently returns negative waiting times, e.g.
    // W = [-6.86, 12.44, -5.59] at rho = 2.483. The port refuses instead.
    PollingMoments<double> m = mk<double>({0.3, 0.5}, {4.0, 3.0}, {32.0, 18.0}, {0.2, 0.15},
                                          {0.02, 0.0225});
    CHECK_THROWS_AS(polling_qsys_exhaustive(m), line::NumericError);
    CHECK_THROWS_AS(polling_qsys_gated(m), line::NumericError);
    // Zero total switchover time: the station-time method divides by R.
    PollingMoments<double> z = mk<double>({0.3}, {1.0}, {2.0}, {0.0}, {0.0});
    CHECK_THROWS_AS(polling_qsys_exhaustive(z), line::InputError);
    CHECK_THROWS_AS(polling_qsys_gated(z), line::InputError);
}

TEST_CASE("polling at Real50 reproduces the double result") {
    const PollingMoments<double> md = case3q();
    PollingMoments<Real50> m;
    for (std::size_t i = 0; i < md.size(); ++i) {
        m.lambda.push_back(Real50(md.lambda[i]));
        m.b.push_back(Real50(md.b[i]));
        m.b2.push_back(Real50(md.b2[i]));
        m.r.push_back(Real50(md.r[i]));
        m.delta2.push_back(Real50(md.delta2[i]));
    }
    const std::vector<Real50> We = polling_qsys_exhaustive(m);
    const std::vector<double> Wed = polling_qsys_exhaustive(md);
    for (std::size_t i = 0; i < md.size(); ++i)
        CHECK(static_cast<double>(We[i]) == doctest::Approx(Wed[i]).epsilon(1e-13));
}

// ---------------------------------------------------------------------------
// M/M/c/c retrial fixed point
// ---------------------------------------------------------------------------

TEST_CASE("the retrial fixed point has a closed form at c = 1") {
    // lambda = 1/2, mu = 2, c = 1. With y = lambda + r, B = y/(2+y) and
    // r = y B gives y^2 = (y - 1/2)(2 + y), i.e. y = 2/3, r = 1/6, B = 1/4.
    const auto r = line::qsys::qsys_mmcc_retrial_fp(0.5, 2.0, 1u);
    CHECK(r.converged);
    CHECK(r.retrialRate == doctest::Approx(1.0 / 6.0).epsilon(1e-9));
    CHECK(r.blockingProbability == doctest::Approx(0.25).epsilon(1e-9));
    // MATLAB: B = 0.24999999998127839, r = 0.16666666660010085, niter = 26.
    CHECK(r.retrialRate == doctest::Approx(0.16666666660010085).epsilon(1e-9));
    CHECK(r.blockingProbability == doctest::Approx(0.24999999998127839).epsilon(1e-9));
}

TEST_CASE("the retrial fixed point agrees with MATLAB and satisfies its own equation") {
    const auto r = line::qsys::qsys_mmcc_retrial_fp(2.0, 1.0, 3u);
    CHECK(r.converged);
    CHECK(r.blockingProbability == doctest::Approx(0.3620285958654505).epsilon(1e-9));
    CHECK(r.retrialRate == doctest::Approx(1.13493674947779).epsilon(1e-9));
    // Residual of r = (lambda + r) B((lambda + r)/mu, c) at the returned r.
    const double B = line::qsys::erlang_b_recursive((2.0 + r.retrialRate) / 1.0, 3u);
    CHECK((2.0 + r.retrialRate) * B == doctest::Approx(r.retrialRate).epsilon(1e-9));
}

TEST_CASE("the retrial fixed point reports non-convergence rather than a wrong number") {
    // lambda = 5 with c = 4 servers of rate 1: the orbit is unstable and r
    // diverges. MATLAB exhausts its 10000 iterations and silently returns
    // r = 10029.175621327076; the port returns the same iterate but flags it.
    const auto r = line::qsys::qsys_mmcc_retrial_fp(5.0, 1.0, 4u);
    CHECK_FALSE(r.converged);
    CHECK(r.iterations == 10000u);
    CHECK(r.retrialRate > 1e3);
    CHECK(r.blockingProbability > 0.99);
}

TEST_CASE("erlang_b_recursive matches its own closed form") {
    // B(a, 1) = a/(1 + a), exactly, in any field.
    CHECK(line::qsys::erlang_b_recursive(Rational(2), 1u) == Rational(2, 3));
    // B(a, 2) = a^2/(2 + 2a + a^2).
    CHECK(line::qsys::erlang_b_recursive(Rational(3), 2u) == Rational(9, 17));
    // B(1/2, 3) = (a^3/6)/(1 + a + a^2/2 + a^3/6) = (1/48)/(79/48) = 1/79.
    CHECK(line::qsys::erlang_b_recursive(Rational(1, 2), 3u) == Rational(1, 79));
}
