/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Finite-capacity and batch queueing systems: Geo^X/Geo/1, M/M/1/K, M/G/1/K
 * (exact and MacGregor Smith), M/M/c/K and M^X/M/1.
 *
 * Oracles, in decreasing strength:
 *   - exact closed forms in rational arithmetic, where the algorithm is a
 *     finite field computation;
 *   - the collapse identities: Geo^X/Geo/1 at beta = 1 is Geo/Geo/1, M/G/1/K
 *     at an exponential density is M/M/1/K, M/G/1/K-MGS at scv = 1 is
 *     M/M/1/K, M/M/c/K at c = 1 and K -> inf is M/M/1;
 *   - loss probabilities lie in [0,1] and decay to zero as K grows;
 *   - MATLAB reference values, obtained by running
 *       matlab -singleCompThread -batch "addpath(genpath('matlab/src')); ..."
 *     against matlab/src/api/qsys, at the tolerance each method guarantees.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/dqsys/dqsys_geogeo1.h"
#include "line/api/dqsys/dqsys_geoxgeo1.h"
#include "line/api/qsys/qsys_mg1k_loss.h"
#include "line/api/qsys/qsys_mg1k_loss_mgs.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mm1k_loss.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/api/qsys/qsys_mxm1.h"

using line::Rational;
using line::Real50;
using line::dqsys::GeoConvention;

TEST_CASE("Geo^X/Geo/1 with a unit batch is Geo/Geo/1, exactly") {
    const Rational a(1, 4), s(1, 2), one(1);
    for (GeoConvention c : {GeoConvention::LAS_DA, GeoConvention::EAS}) {
        auto bx = line::dqsys::dqsys_geoxgeo1(a, one, s, c);
        auto g = line::dqsys::dqsys_geogeo1(a, s, c);
        CHECK(bx.meanQueueLength == g.meanQueueLength);
        CHECK(bx.meanSojournTime == g.meanSojournTime);
        CHECK(bx.meanWaitingTime == g.meanWaitingTime);
        CHECK(bx.meanServiceTime == g.meanServiceTime);
        CHECK(bx.meanWaitingQueue == g.meanWaitingQueue);
        CHECK(bx.utilization == g.utilization);
    }
}

TEST_CASE("Geo^X/Geo/1 obeys Little's law and the epoch shift, exactly") {
    // a = 1/10, beta = 1/2 (mean batch 2), s = 9/10: lambda = 1/5.
    const Rational a(1, 10), beta(1, 2), s(9, 10);
    auto las = line::dqsys::dqsys_geoxgeo1(a, beta, s, GeoConvention::LAS_DA);
    auto eas = line::dqsys::dqsys_geoxgeo1(a, beta, s, GeoConvention::EAS);
    CHECK(las.arrivalRate == Rational(1, 5));
    CHECK(las.utilization == Rational(2, 9));
    CHECK(las.meanQueueLength == las.arrivalRate * las.meanSojournTime);
    CHECK(las.meanWaitingQueue == las.arrivalRate * las.meanWaitingTime);
    CHECK(las.meanSojournTime == las.meanWaitingTime + las.meanServiceTime);
    // The EAS epoch drops exactly one slot of sojourn and lambda jobs.
    CHECK(eas.meanQueueLength == las.meanQueueLength - las.arrivalRate);
    CHECK(eas.meanSojournTime == las.meanSojournTime - Rational(1));
    // The waiting time is epoch independent.
    CHECK(eas.meanWaitingTime == las.meanWaitingTime);
}

TEST_CASE("Geo^X/Geo/1 matches MATLAB dqsys_geoxgeo1(0.1, 0.5, 0.9)") {
    // Rational closed form, so the only error is MATLAB's double rounding:
    // agreement is asserted at 1e-14 relative.
    auto r = line::dqsys::dqsys_geoxgeo1(Rational(1, 10), Rational(1, 2), Rational(9, 10));
    CHECK(static_cast<double>(r.arrivalRate) == doctest::Approx(0.20000000000000001).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanQueueLength) ==
          doctest::Approx(0.51428571428571435).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanWaitingTime) ==
          doctest::Approx(1.4603174603174605).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanSojournTime) ==
          doctest::Approx(2.5714285714285716).epsilon(1e-14));
    auto e = line::dqsys::dqsys_geoxgeo1(Rational(1, 10), Rational(1, 2), Rational(9, 10),
                                       GeoConvention::EAS);
    CHECK(static_cast<double>(e.meanQueueLength) ==
          doctest::Approx(0.31428571428571433).epsilon(1e-14));
    CHECK(static_cast<double>(e.meanSojournTime) ==
          doctest::Approx(1.5714285714285716).epsilon(1e-14));
}

TEST_CASE("Geo^X/Geo/1 pgf is one at z = 1 and rejects bad moments") {
    auto r = line::dqsys::dqsys_geoxgeo1(Rational(1, 10), Rational(1, 2), Rational(9, 10));
    CHECK(line::dqsys::dqsys_geoxgeo1_pgf(r, Rational(1), Rational(1)) == Rational(1));
    CHECK_THROWS_AS(line::dqsys::dqsys_geoxgeo1_pgf(r, Rational(0), Rational(1)), line::InputError);
    // E[X(X-1)] below E[X]^2-E[X] describes no random variable.
    CHECK_THROWS_AS(line::dqsys::dqsys_geoxgeo1_moments(Rational(1, 10), Rational(2), Rational(0),
                                                      Rational(9, 10)),
                    line::InputError);
    // A mean batch below one is not a batch.
    CHECK_THROWS_AS(line::dqsys::dqsys_geoxgeo1_moments(Rational(1, 10), Rational(1, 2), Rational(0),
                                                      Rational(9, 10)),
                    line::InputError);
    // Unstable load.
    CHECK_THROWS_AS(
        line::dqsys::dqsys_geoxgeo1(Rational(1, 2), Rational(1, 2), Rational(9, 10)),
        line::InputError);
}

TEST_CASE("M/M/1/K loss is exact and matches MATLAB") {
    // MATLAB: qsys_mm1k_loss(2,3,5) -> 0.048120300751879688, rho = 2/3.
    auto r = line::qsys::qsys_mm1k_loss(Rational(2), Rational(3), 5u);
    CHECK(r.utilization == Rational(2, 3));
    // Exact value: (1-rho) rho^K / (1-rho^(K+1)) with rho = 2/3, K = 5.
    const Rational rho(2, 3);
    const Rational expected = (Rational(1) - rho) * line::num_pow_int(rho, 5u) /
                              (Rational(1) - line::num_pow_int(rho, 6u));
    CHECK(r.lossProbability == expected);
    CHECK(static_cast<double>(r.lossProbability) ==
          doctest::Approx(0.048120300751879688).epsilon(1e-14));
}

TEST_CASE("M/M/1/K loss lies in [0,1] and vanishes as K grows") {
    const Rational lambda(2), mu(3);
    Rational prev(1);
    for (unsigned K = 1; K <= 40; ++K) {
        auto r = line::qsys::qsys_mm1k_loss(lambda, mu, K);
        CHECK(r.lossProbability >= Rational(0));
        CHECK(r.lossProbability <= Rational(1));
        CHECK(r.lossProbability < prev);  // strictly decreasing in K
        prev = r.lossProbability;
    }
    CHECK(prev < Rational(1, 1000000));
    CHECK_THROWS_AS(line::qsys::qsys_mm1k_loss(Rational(1), Rational(1), 4u), line::InputError);
}

TEST_CASE("M/G/1/K loss reproduces M/M/1/K at an exponential density") {
    // The embedded chain must collapse onto the closed form. MATLAB's own run
    // agrees with the closed form to 1.1e-11 absolute, which is the accuracy
    // its integral() tolerances (RelTol 1e-6) buy; the port is asserted at the
    // same 1e-9 absolute.
    const double mu = 3.0;
    auto r = line::qsys::qsys_mg1k_loss<double>(2.0, [&](double t) { return mu * std::exp(-mu * t); },
                                                5u);
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    CHECK(r.lossProbability == doctest::Approx(0.048120300751879688).epsilon(1e-9));
    // MATLAB qsys_mg1k_loss(2,@(t)3*exp(-3*t),5) -> 0.048120300751869016.
    CHECK(std::abs(r.lossProbability - 0.048120300751869016) < 1e-9);
}

TEST_CASE("M/G/1/K loss matches MATLAB on an Erlang-2 service") {
    // Erlang-2 with rate 6: mean 1/3, scv 1/2, so rho = 2/3 as above but with
    // half the service variability, and the loss must be strictly smaller.
    auto r = line::qsys::qsys_mg1k_loss<double>(
        2.0, [](double t) { return 36.0 * t * std::exp(-6.0 * t); }, 5u);
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    // MATLAB: 0.029920129229096482, at its RelTol 1e-6 quadrature.
    CHECK(r.lossProbability == doctest::Approx(0.029920129229096482).epsilon(1e-7));
    CHECK(r.lossProbability < 0.048120300751879688);
    CHECK(r.lossProbability > 0.0);
}

TEST_CASE("M/G/1/K loss decays as capacity grows") {
    const double mu = 3.0;
    auto f = [&](double t) { return mu * std::exp(-mu * t); };
    double prev = 1.0;
    for (unsigned K : {2u, 3u, 5u, 8u, 12u}) {
        auto r = line::qsys::qsys_mg1k_loss<double>(2.0, f, K);
        CHECK(r.lossProbability >= 0.0);
        CHECK(r.lossProbability <= 1.0);
        CHECK(r.lossProbability < prev);
        prev = r.lossProbability;
    }
    CHECK(prev < 0.01);
    CHECK_THROWS_AS(line::qsys::qsys_mg1k_loss<double>(2.0, f, 1u), line::InputError);
}

TEST_CASE("MacGregor Smith M/G/1/K reduces to M/M/1/K at scv one") {
    // At scv = 1 the exponents are exactly K and K+1, so the approximation is
    // the exact M/M/1/K formula; assert that at round-off (1e-12 relative).
    auto approx = line::qsys::qsys_mg1k_loss_mgs(2.0, 3.0, 1.0, 5u);
    auto exact = line::qsys::qsys_mm1k_loss(2.0, 3.0, 5u);
    CHECK(approx.lossProbability == doctest::Approx(exact.lossProbability).epsilon(1e-12));
    // MATLAB qsys_mg1k_loss_mgs(2,3,1,5) -> 0.048120300751879688.
    CHECK(approx.lossProbability == doctest::Approx(0.048120300751879688).epsilon(1e-12));
    // MATLAB qsys_mg1k_loss_mgs(2,3,2,5) -> 0.081727019440218485.
    auto v2 = line::qsys::qsys_mg1k_loss_mgs(2.0, 3.0, 2.0, 5u);
    CHECK(v2.lossProbability == doctest::Approx(0.081727019440218485).epsilon(1e-13));
    CHECK(v2.lossProbability > approx.lossProbability);  // more variability, more loss
}

TEST_CASE("MacGregor Smith loss stays a probability and decays with K") {
    double prev = 1.0;
    for (unsigned K : {2u, 4u, 8u, 16u, 32u}) {
        auto r = line::qsys::qsys_mg1k_loss_mgs(2.0, 3.0, 2.0, K);
        CHECK(r.lossProbability >= 0.0);
        CHECK(r.lossProbability <= 1.0);
        CHECK(r.lossProbability < prev);
        prev = r.lossProbability;
    }
    CHECK(prev < 1e-3);
}

TEST_CASE("M/M/c/K is exact and matches MATLAB qsys_mmck(1,1,2,4)") {
    auto r = line::qsys::qsys_mmck(Rational(1), Rational(1), 2u, 4u);
    // The distribution is a probability vector, exactly.
    Rational mass(0);
    for (const Rational& v : r.queueLengthDist) {
        CHECK(v >= Rational(0));
        mass += v;
    }
    CHECK(mass == Rational(1));
    // Little's law on the carried load, exactly.
    CHECK(r.meanQueueLength == r.throughput * r.meanSojournTime);
    CHECK(r.meanQueueLengthQ == r.throughput * r.meanWaitingTime);
    // MATLAB reference values, exact rational arithmetic against MATLAB
    // doubles: 1e-14 relative.
    CHECK(static_cast<double>(r.meanQueueLength) ==
          doctest::Approx(1.1304347826086956).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanQueueLengthQ) ==
          doctest::Approx(0.17391304347826086).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanWaitingTime) ==
          doctest::Approx(0.1818181818181818).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanSojournTime) ==
          doctest::Approx(1.1818181818181817).epsilon(1e-14));
    CHECK(static_cast<double>(r.utilization) == doctest::Approx(0.47826086956521741).epsilon(1e-14));
    CHECK(static_cast<double>(r.throughput) == doctest::Approx(0.95652173913043481).epsilon(1e-14));
    CHECK(static_cast<double>(r.lossProbability) ==
          doctest::Approx(0.043478260869565216).epsilon(1e-14));
}

TEST_CASE("M/M/c/K blocking goes to zero and W to the M/M/1 value as K grows") {
    const Rational lambda(1, 2), mu(1);
    Rational prev(1);
    for (unsigned K = 1; K <= 30; ++K) {
        auto r = line::qsys::qsys_mmck(lambda, mu, 1u, K);
        CHECK(r.lossProbability >= Rational(0));
        CHECK(r.lossProbability <= Rational(1));
        CHECK(r.lossProbability < prev);
        prev = r.lossProbability;
    }
    auto big = line::qsys::qsys_mmck(lambda, mu, 1u, 60u);
    auto mm1 = line::qsys::qsys_mm1(lambda, mu);
    CHECK(static_cast<double>(big.meanSojournTime) ==
          doctest::Approx(static_cast<double>(mm1.W)).epsilon(1e-12));
    CHECK_THROWS_AS(line::qsys::qsys_mmck(Rational(1), Rational(1), 3u, 2u), line::InputError);
}

TEST_CASE("M/M/c/K at c = K is the Erlang-B loss system") {
    // With no waiting room, p_K is the Erlang-B blocking of a = lambda/mu.
    const Rational a(3), one(1);
    const unsigned c = 4;
    auto r = line::qsys::qsys_mmck(a, one, c, c);
    Rational num = line::num_pow_int(a, c) / line::num_factorial<Rational>(c);
    Rational den(0);
    for (unsigned n = 0; n <= c; ++n)
        den += line::num_pow_int(a, n) / line::num_factorial<Rational>(n);
    CHECK(r.lossProbability == num / den);
    CHECK(r.meanQueueLengthQ == Rational(0));  // nobody ever waits
}

TEST_CASE("M^X/M/1 is exact and matches MATLAB") {
    // MATLAB qsys_mxm1(0.3, 2, 2, 6) -> W 1.4285714285714286, Wq 0.9285714285714286,
    // U 0.3, Q 0.8571428571428571.
    auto r = line::qsys::qsys_mxm1(Rational(3, 10), Rational(2), Rational(2), Rational(6));
    CHECK(r.U == Rational(3, 10));
    CHECK(r.W == r.Wq + Rational(1, 2));
    CHECK(r.Q == Rational(3, 10) * Rational(2) * r.W);  // Little's law, exactly
    CHECK(static_cast<double>(r.W) == doctest::Approx(1.4285714285714286).epsilon(1e-14));
    CHECK(static_cast<double>(r.Wq) == doctest::Approx(0.9285714285714286).epsilon(1e-14));
    CHECK(static_cast<double>(r.Q) == doctest::Approx(0.8571428571428571).epsilon(1e-14));
}

TEST_CASE("M^X/M/1 with a unit batch is M/M/1, exactly") {
    const Rational lambda(1, 2), mu(1), one(1);
    auto bx = line::qsys::qsys_mxm1(lambda, mu, one, one);  // X == 1 a.s.
    auto mm1 = line::qsys::qsys_mm1(lambda, mu);
    CHECK(bx.W == mm1.W);
    CHECK(bx.U == mm1.rhohat);
}

TEST_CASE("M^X/M/1 accepts the three MATLAB input forms consistently") {
    // Batch on {1,2,3} uniform: E[X] = 2, E[X^2] = 14/3, Var = 2/3.
    const Rational lambda(1, 10), mu(1);
    std::vector<Rational> sizes = {Rational(1), Rational(2), Rational(3)};
    std::vector<Rational> pmf = {Rational(1), Rational(1), Rational(1)};  // unnormalized
    auto a = line::qsys::qsys_mxm1_pmf(lambda, mu, sizes, pmf);
    auto b = line::qsys::qsys_mxm1(lambda, mu, Rational(2), Rational(14, 3));
    auto c = line::qsys::qsys_mxm1_variance(lambda, mu, Rational(2), Rational(2, 3));
    CHECK(a.W == b.W);
    CHECK(a.W == c.W);
    CHECK_THROWS_AS(line::qsys::qsys_mxm1(Rational(1), Rational(1), Rational(2), Rational(6)),
                    line::InputError);  // rho = 2 >= 1
}

TEST_CASE("the finite-capacity family agrees between double and Real50") {
    auto d = line::qsys::qsys_mmck(1.0, 1.0, 2u, 4u);
    auto h = line::qsys::qsys_mmck(Real50(1), Real50(1), 2u, 4u);
    CHECK(static_cast<double>(h.meanQueueLength) == doctest::Approx(d.meanQueueLength).epsilon(1e-14));
    auto dl = line::qsys::qsys_mg1k_loss_mgs(2.0, 3.0, 2.0, 5u);
    auto hl = line::qsys::qsys_mg1k_loss_mgs(Real50(2), Real50(3), Real50(2), 5u);
    CHECK(static_cast<double>(hl.lossProbability) ==
          doctest::Approx(dl.lossProbability).epsilon(1e-14));
}
