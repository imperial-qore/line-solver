/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Polling waiting times. Reference values come from the MATLAB implementations
 * (polling_qsys_1limited.m, polling_qsys_decrementing.m) evaluated on the same
 * moment inputs; the structural properties (symmetry, monotonicity in load,
 * the stability pole) are asserted directly.
 */
#include <vector>

#include "doctest.h"
#include "line/api/polling/polling_qsys_1limited.h"

using line::Rational;
using line::polling::PollingMoments;
using line::polling::polling_qsys_1limited;
using line::polling::polling_qsys_decrementing;

namespace {

/** Symmetric 3-queue system: Poisson arrivals 1/10, exponential service 1. */
template <class T>
PollingMoments<T> symmetric3() {
    PollingMoments<T> m;
    for (int i = 0; i < 3; ++i) {
        m.lambda.push_back(line::num_traits<T>::from_rational(1, 10));
        m.b.push_back(line::num_traits<T>::from_int(1));
        m.b2.push_back(line::num_traits<T>::from_int(2));            // exponential: 2 b^2
        m.r.push_back(line::num_traits<T>::from_rational(1, 5));     // switchover mean
        m.delta2.push_back(line::num_traits<T>::from_rational(1, 25));  // exponential switchover
    }
    return m;
}

}  // namespace

TEST_CASE("1-limited polling is symmetric on a symmetric system, exactly") {
    const std::vector<Rational> W = polling_qsys_1limited(symmetric3<Rational>());
    CHECK(W.size() == 3);
    CHECK(W[0] == W[1]);
    CHECK(W[1] == W[2]);
    CHECK(W[0] > Rational(0));
}

TEST_CASE("1-limited waiting time grows with the load") {
    PollingMoments<Rational> light = symmetric3<Rational>();
    PollingMoments<Rational> heavy = symmetric3<Rational>();
    for (std::size_t i = 0; i < heavy.size(); ++i) heavy.lambda[i] = Rational(1, 5);
    const Rational wl = polling_qsys_1limited(light)[0];
    const Rational wh = polling_qsys_1limited(heavy)[0];
    CHECK(wh > wl);
}

TEST_CASE("1-limited rejects an unstable queue at the pole") {
    // 1 - rho - lambda_i R = 0 exactly: three queues with b = 1 and r = 1/5
    // give rho = 3 lambda and R = 3/5, so the denominator vanishes at
    // lambda (3 + 3/5) = 1, i.e. lambda = 5/18: rho = 5/6, lambda R = 1/6.
    PollingMoments<Rational> m = symmetric3<Rational>();
    for (std::size_t i = 0; i < m.size(); ++i) m.lambda[i] = Rational(5, 18);
    CHECK_THROWS_AS(polling_qsys_1limited(m), line::NumericError);
}

TEST_CASE("decrementing polling matches the closed form and is symmetric") {
    const PollingMoments<Rational> m = symmetric3<Rational>();
    const std::vector<Rational> W = polling_qsys_decrementing(m);
    CHECK(W[0] == W[1]);
    CHECK(W[1] == W[2]);

    // Hand evaluation: N=3, lam=1/10, b=1, b2=2, r=1/5, d2=1/25.
    // rho = 3/10; denom = 2(1 - 3/10 - (1/10)(1/5)(3 - 3/10)) = 2(7/10 - 27/500)
    //      = 2 * 323/500 = 323/250.
    // W = (1/25)/(2/5) + (3*(1/10)*2*(1 - 1/50) + (1/5 + (1/10)(1/25))(3 - 3/10)) / (323/250)
    const Rational rho(3, 10);
    const Rational denom = Rational(2) * (Rational(1) - rho - Rational(1, 10) * Rational(1, 5) *
                                                                 (Rational(3) - rho));
    const Rational expected =
        Rational(1, 25) / (Rational(2) * Rational(1, 5)) +
        (Rational(3) * Rational(1, 10) * Rational(2) * (Rational(1) - Rational(1, 10) * Rational(1, 5)) +
         (Rational(1, 5) + Rational(1, 10) * Rational(1, 25)) * (Rational(3) - rho)) /
            denom;
    CHECK(W[0] == expected);
}

TEST_CASE("decrementing polling rejects an asymmetric system") {
    PollingMoments<Rational> m = symmetric3<Rational>();
    m.lambda[1] = Rational(1, 9);
    CHECK_THROWS_AS(polling_qsys_decrementing(m), line::InputError);
}

TEST_CASE("polling formulas agree between double and exact arithmetic") {
    const PollingMoments<Rational> mq = symmetric3<Rational>();
    const PollingMoments<double> md = symmetric3<double>();
    CHECK(static_cast<double>(polling_qsys_1limited(mq)[0]) ==
          doctest::Approx(polling_qsys_1limited(md)[0]).epsilon(1e-12));
    CHECK(static_cast<double>(polling_qsys_decrementing(mq)[0]) ==
          doctest::Approx(polling_qsys_decrementing(md)[0]).epsilon(1e-12));
}
