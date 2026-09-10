/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Elementary closed-form queueing systems: M/M/1, M/M/k, M/G/1, G/M/1 and
 * M/G/infinity. Oracles are the textbook formulas themselves, evaluated
 * independently of the code under test, plus the exact identities the port
 * must satisfy (M/M/k at k = 1 is M/M/1; M/G/1 at cs = 1 is M/M/1).
 */
#include <cmath>

#include "doctest.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mginf.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_gm1.h"

using line::Rational;
using line::Real50;


namespace {
constexpr double REL_TOL = 1e-9;
}  // namespace

TEST_CASE("qsys_mm1 lambda=1/2 mu=1 gives W=2 exactly in rational arithmetic") {
    auto r = line::qsys::qsys_mm1(Rational(1, 2), Rational(1));
    CHECK(r.W == Rational(2));
    CHECK(r.rhohat == Rational(1, 2));

    auto rd = line::qsys::qsys_mm1(0.5, 1.0);
    CHECK(rd.W == doctest::Approx(2.0).epsilon(REL_TOL));

    auto rr = line::qsys::qsys_mm1(Real50(1) / Real50(2), Real50(1));
    CHECK(static_cast<double>(rr.W) == doctest::Approx(2.0).epsilon(REL_TOL));
}

TEST_CASE("qsys_mm1 W = 1/(mu-lambda), exact and double agree") {
    // lambda = 3/7, mu = 5/4: W = 1/(5/4-3/7) = 28/23.
    auto rq = line::qsys::qsys_mm1(Rational(3, 7), Rational(5, 4));
    CHECK(rq.W == Rational(28, 23));
    auto rd = line::qsys::qsys_mm1(3.0 / 7.0, 5.0 / 4.0);
    CHECK(rd.W == doctest::Approx(28.0 / 23.0).epsilon(REL_TOL));
}

TEST_CASE("qsys_mm1 W increases with rho") {
    double prev = 0.0;
    for (double lam = 0.1; lam < 0.95; lam += 0.1) {
        double W = line::qsys::qsys_mm1(lam, 1.0).W;
        CHECK(W > prev);
        prev = W;
    }
}

TEST_CASE("qsys_mm1 rejects the pole at rho == 1") {
    CHECK_THROWS_AS(line::qsys::qsys_mm1(Rational(1), Rational(1)), line::InputError);
}

TEST_CASE("qsys_mmk with k=1 equals qsys_mm1 exactly") {
    const Rational lam(3, 7), mu(5, 4);
    auto a = line::qsys::qsys_mmk(lam, mu, 1u);
    auto b = line::qsys::qsys_mm1(lam, mu);
    CHECK(a.W == b.W);
    CHECK(a.rhohat == b.rhohat);

    // Same identity at a second operating point.
    const Rational lam2(1, 3), mu2(1);
    CHECK(line::qsys::qsys_mmk(lam2, mu2, 1u).W == line::qsys::qsys_mm1(lam2, mu2).W);
}

TEST_CASE("qsys_mmk matches an independent Erlang-C evaluation") {
    // lambda = 2, mu = 1, k = 3: rho = 2/3.
    const double lam = 2.0, mu = 1.0;
    const int k = 3;
    const double rho = lam / mu / k;
    double S = 0.0, f = 1.0;
    for (int j = 0; j < k; ++j) {
        if (j > 0) f *= j;
        S += std::pow(k * rho, j) / f;
    }
    double kfact = 1.0;
    for (int j = 2; j <= k; ++j) kfact *= j;
    const double C = 1.0 / (1.0 + (1.0 - rho) * kfact / std::pow(k * rho, k) * S);
    const double Wexp = (rho / (1.0 - rho) * C + k * rho) / lam;

    CHECK(line::qsys::qsys_mmk(lam, mu, 3u).W == doctest::Approx(Wexp).epsilon(REL_TOL));
    CHECK(static_cast<double>(line::qsys::qsys_mmk(Rational(2), Rational(1), 3u).W) ==
          doctest::Approx(Wexp).epsilon(REL_TOL));
    CHECK(static_cast<double>(line::qsys::qsys_mmk(Real50(2), Real50(1), 3u).W) ==
          doctest::Approx(Wexp).epsilon(REL_TOL));
}

TEST_CASE("qsys_mmk W increases with rho at fixed k") {
    double prev = 0.0;
    for (double lam = 0.5; lam < 2.9; lam += 0.4) {
        double W = line::qsys::qsys_mmk(lam, 1.0, 3u).W;
        CHECK(W > prev);
        prev = W;
    }
}

TEST_CASE("qsys_mmk rejects k=0") {
    CHECK_THROWS_AS(line::qsys::qsys_mmk(0.5, 1.0, 0u), line::InputError);
}

TEST_CASE("qsys_mg1 with cs=1 reduces to M/M/1 exactly") {
    const Rational lam(1, 2), mu(1);
    auto r = line::qsys::qsys_mg1(lam, mu, Rational(1));
    CHECK(r.W == line::qsys::qsys_mm1(lam, mu).W);
    // rhohat = Q/(1+Q) with Q = rho/(1-rho) = 1, so rhohat = 1/2 = rho.
    CHECK(r.rhohat == Rational(1, 2));
}

TEST_CASE("qsys_mg1 Pollaczek-Khinchine, exact vs double") {
    // lambda = 1/2, mu = 1, cs = 2: Q = 1/2 + 1/4/1 + (1/4)(4)/1/1 = 1/2+1/4+1 ... exact below.
    const Rational lam(1, 2), mu(1), cs(2);
    auto rq = line::qsys::qsys_mg1(lam, mu, cs);
    // Q = rho + rho^2/(2(1-rho)) + lambda^2 cs^2/mu^2/(2(1-rho))
    //   = 1/2 + (1/4)/1 + (1/4*4)/1 = 1/2 + 1/4 + 1 = 7/4;  W = Q/lambda = 7/2.
    CHECK(rq.W == Rational(7, 2));
    CHECK(rq.rhohat == Rational(7, 11));
    auto rd = line::qsys::qsys_mg1(0.5, 1.0, 2.0);
    CHECK(rd.W == doctest::Approx(3.5).epsilon(REL_TOL));
    CHECK(rd.rhohat == doctest::Approx(7.0 / 11.0).epsilon(REL_TOL));
}

TEST_CASE("qsys_gm1 is 1/((1-sigma) mu)") {
    CHECK(line::qsys::qsys_gm1(Rational(1, 2), Rational(1)) == Rational(2));
    CHECK(line::qsys::qsys_gm1(0.25, 2.0) == doctest::Approx(1.0 / 0.75 / 2.0).epsilon(REL_TOL));
    CHECK(static_cast<double>(line::qsys::qsys_gm1(Real50(1) / Real50(4), Real50(2))) ==
          doctest::Approx(1.0 / 0.75 / 2.0).epsilon(REL_TOL));
    CHECK_THROWS_AS(line::qsys::qsys_gm1(Rational(1), Rational(1)), line::InputError);
}

TEST_CASE("qsys_gm1 W increases with sigma") {
    double prev = 0.0;
    for (double s = 0.1; s < 0.95; s += 0.1) {
        double W = line::qsys::qsys_gm1(s, 1.0);
        CHECK(W > prev);
        prev = W;
    }
}

TEST_CASE("qsys_mginf is Poisson(rho) and independent of the service distribution") {
    const double lam = 3.0, mu = 2.0, rho = 1.5;
    auto r = line::qsys::qsys_mginf(lam, mu);
    CHECK(r.L == doctest::Approx(rho).epsilon(REL_TOL));
    CHECK(r.Lq == 0.0);
    CHECK(r.W == doctest::Approx(0.5).epsilon(REL_TOL));
    CHECK(r.Wq == 0.0);
    CHECK(r.p0 == doctest::Approx(std::exp(-rho)).epsilon(REL_TOL));
    CHECK_FALSE(r.has_pk);

    auto r3 = line::qsys::qsys_mginf(lam, mu, 3u);
    CHECK(r3.has_pk);
    CHECK(r3.pk == doctest::Approx(std::exp(-rho) * rho * rho * rho / 6.0).epsilon(REL_TOL));

    // High precision agrees with double.
    auto rr = line::qsys::qsys_mginf(Real50(3), Real50(2), 3u);
    CHECK(static_cast<double>(rr.pk) == doctest::Approx(r3.pk).epsilon(REL_TOL));
    CHECK(static_cast<double>(rr.p0) == doctest::Approx(r.p0).epsilon(REL_TOL));
}
