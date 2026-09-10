/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Geo/Geo/1. Oracles: the pmf must be a probability distribution summing to
 * one, Little's law must hold on the slot time scale, and the two conventions
 * must give the documented distinct answers. All exact in the rational
 * instantiation.
 */
#include <vector>

#include "doctest.h"
#include "line/api/dqsys/dqsys_geogeo1.h"

using line::Rational;
using line::dqsys::GeoConvention;
using line::dqsys::dqsys_geogeo1;
using line::dqsys::dqsys_geogeo1_pmf;

TEST_CASE("Geo/Geo/1 LAS_DA closed forms are exact") {
    // a = 1/4, s = 1/2: rho = 1/2, r = a(1-s)/(s(1-a)) = (1/4)(1/2)/((1/2)(3/4)) = 1/3.
    const Rational a(1, 4), s(1, 2);
    auto r = dqsys_geogeo1(a, s);
    CHECK(r.utilization == Rational(1, 2));
    CHECK(r.ratio == Rational(1, 3));
    CHECK(r.emptyProb == Rational(1, 2));
    CHECK(r.meanQueueLength == a * (Rational(1) - a) / (s - a));  // 3/4
    CHECK(r.meanQueueLength == Rational(3, 4));
    CHECK(r.meanSojournTime == Rational(3));
    CHECK(r.meanServiceTime == Rational(2));
    CHECK(r.meanWaitingTime == Rational(1));
    CHECK(r.meanWaitingQueue == Rational(1, 4));
}

TEST_CASE("the two conventions genuinely differ") {
    const Rational a(1, 4), s(1, 2);
    auto las = dqsys_geogeo1(a, s, GeoConvention::LAS_DA);
    auto eas = dqsys_geogeo1(a, s, GeoConvention::EAS);
    CHECK(las.emptyProb != eas.emptyProb);
    CHECK(las.meanQueueLength != eas.meanQueueLength);
    CHECK(las.meanServiceTime != eas.meanServiceTime);
    // EAS: empty probability is 1 - r, queue length a(1-s)/(s-a) = 1/2.
    CHECK(eas.emptyProb == Rational(2, 3));
    CHECK(eas.meanQueueLength == Rational(1, 2));
    // The waiting time does not depend on the convention.
    CHECK(las.meanWaitingTime == eas.meanWaitingTime);
}

TEST_CASE("the pmf is a distribution and reproduces the mean") {
    const Rational a(1, 4), s(1, 2);
    for (GeoConvention c : {GeoConvention::LAS_DA, GeoConvention::EAS}) {
        auto r = dqsys_geogeo1(a, s, c);
        // Geometric tail: sum over 0..K plus the exact remainder is 1.
        Rational mass(0), mean(0);
        const int K = 60;
        for (int n = 0; n <= K; ++n) {
            const Rational p = dqsys_geogeo1_pmf(r, n);
            CHECK(p >= Rational(0));
            mass += p;
            mean += Rational(n) * p;
        }
        // With ratio 1/3 the truncated tail beyond K = 60 is below 3^-60.
        CHECK(mass <= Rational(1));
        CHECK(Rational(1) - mass < Rational(1, 1000000000000000000LL));
        CHECK(static_cast<double>(mean) ==
              doctest::Approx(static_cast<double>(r.meanQueueLength)).epsilon(1e-12));
    }
}

TEST_CASE("Little's law holds on the slot time scale") {
    const Rational a(1, 5), s(2, 5);
    auto las = dqsys_geogeo1(a, s);
    // L = lambda W with lambda = a.
    CHECK(las.meanQueueLength == a * las.meanSojournTime);
    CHECK(las.meanWaitingQueue == a * las.meanWaitingTime);
    // Sojourn time is waiting plus service.
    CHECK(las.meanSojournTime == las.meanWaitingTime + las.meanServiceTime);
}

TEST_CASE("Geo/Geo/1 rejects unstable or out-of-range parameters") {
    CHECK_THROWS_AS(dqsys_geogeo1(Rational(1, 2), Rational(1, 4)), line::InputError);  // a >= s
    CHECK_THROWS_AS(dqsys_geogeo1(Rational(0), Rational(1, 2)), line::InputError);
    CHECK_THROWS_AS(dqsys_geogeo1(Rational(1, 4), Rational(3, 2)), line::InputError);
    CHECK_THROWS_AS(dqsys_geogeo1_pmf(dqsys_geogeo1(Rational(1, 4), Rational(1, 2)), -1),
                    line::InputError);
}

TEST_CASE("Geo/Geo/1 agrees between double and exact arithmetic") {
    auto q = dqsys_geogeo1(Rational(1, 4), Rational(1, 2));
    auto d = dqsys_geogeo1(0.25, 0.5);
    CHECK(static_cast<double>(q.meanQueueLength) ==
          doctest::Approx(d.meanQueueLength).epsilon(1e-12));
    CHECK(static_cast<double>(q.meanWaitingTime) ==
          doctest::Approx(d.meanWaitingTime).epsilon(1e-12));
}
