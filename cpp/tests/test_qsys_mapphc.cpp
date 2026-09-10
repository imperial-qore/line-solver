/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The exact MAP/PH/c queue.
 *
 * THE ORACLES ARE INDEPENDENT OF THE IMPLEMENTATION. With Poisson arrivals and
 * exponential service the model collapses to M/M/c, whose queue length, delay
 * probability, mean wait, second moment and waiting-time CCDF all have closed
 * forms (Erlang-C). With a correlated MMPP2 arrival and exponential service it
 * must reproduce qsys_mapmc, which solves a DIFFERENT chain: Q-MAM's
 * level-dependent QBD, whose phase is the arrival phase alone. The MATLAB
 * reference reports 1.5218751032 for that instance.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mapmc.h"
#include "line/api/qsys/qsys_mapphc.h"

namespace qsys = line::qsys;
namespace mam = line::mam;
using line::Matrix;

namespace {

mam::Map<double> poisson_map(double lam) {
    mam::Map<double> m;
    m.D0 = Matrix<double>(1, 1, -lam);
    m.D1 = Matrix<double>(1, 1, lam);
    return m;
}

mam::Map<double> mmpp2_map() {
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -(1.8 + 0.15);
    m.D0(0, 1) = 0.15;
    m.D0(1, 0) = 0.25;
    m.D0(1, 1) = -(0.4 + 0.25);
    m.D1(0, 0) = 1.8;
    m.D1(1, 1) = 0.4;
    return m;
}

}  // namespace

TEST_CASE("qsys_mapphc: exponential service collapses onto Erlang-C") {
    const double lam = 1.2, mu = 1.0;
    const unsigned c = 3;
    std::vector<double> alpha(1, 1.0);
    Matrix<double> S(1, 1, -mu);
    std::vector<double> pts;
    pts.push_back(0.1);
    pts.push_back(0.5);
    pts.push_back(1.0);
    const qsys::MapPhcResult<double> r = qsys::qsys_mapphc(poisson_map(lam), alpha, S, c, 500, 3, pts);

    const double a = lam / mu, rho = a / c;
    double sum = 0.0, fact = 1.0;
    for (unsigned n = 0; n < c; ++n) {
        if (n > 0) fact *= n;
        sum += std::pow(a, static_cast<double>(n)) / fact;
    }
    const double factC = fact * c;
    const double p0 = 1.0 / (sum + std::pow(a, static_cast<double>(c)) / (factC * (1 - rho)));
    const double cErl = std::pow(a, static_cast<double>(c)) / (factC * (1 - rho)) * p0;
    const double lq = cErl * rho / (1 - rho);

    CHECK(r.meanQueueLength == doctest::Approx(lq + a).epsilon(1e-9));
    CHECK(r.meanWaitingTime == doctest::Approx(lq / lam).epsilon(1e-9));
    CHECK(r.probWait == doctest::Approx(cErl).epsilon(1e-9));
    CHECK(r.waitingTimeMoments[1] ==
          doctest::Approx(2 * cErl / std::pow(c * mu - lam, 2.0)).epsilon(1e-9));
    for (std::size_t i = 0; i < pts.size(); ++i) {
        CHECK(r.waitingTimeCCDF[i] ==
              doctest::Approx(cErl * std::exp(-(c * mu - lam) * pts[i])).epsilon(1e-9));
    }
    CHECK(r.phaseCount == 1u);
}

TEST_CASE("qsys_mapphc: a correlated arrival reproduces the independent MAP/M/c chain") {
    const unsigned c = 3;
    std::vector<double> alpha(1, 1.0);
    Matrix<double> S(1, 1, -1.0);
    const qsys::MapPhcResult<double> r = qsys::qsys_mapphc(mmpp2_map(), alpha, S, c);
    const qsys::MapMcResult<double> ref = qsys::qsys_mapmc(mmpp2_map(), 1.0, c);
    CHECK(r.meanQueueLength == doctest::Approx(ref.meanQueueLength).epsilon(1e-6));
    CHECK(r.meanWaitingTime == doctest::Approx(ref.meanWaitingTime).epsilon(1e-6));
    CHECK(r.meanQueueLength == doctest::Approx(1.5218751032).epsilon(1e-6));
}

TEST_CASE("qsys_mapphc: the phase count is the multiset count, not the power count") {
    // Erlang-2 service, ms = 2: the ordered space at c = 3 is 8, the multiset
    // space is binomial(4,3) = 4.
    std::vector<double> alpha(2, 0.0);
    alpha[0] = 1.0;
    Matrix<double> S(2, 2, 0.0);
    S(0, 0) = -2.0;
    S(0, 1) = 2.0;
    S(1, 1) = -2.0;
    const qsys::MapPhcResult<double> r = qsys::qsys_mapphc(poisson_map(0.9), alpha, S, 3);
    CHECK(r.phaseCount == 4u);
    // Little's law on the station ties the two means together
    CHECK(r.meanSojournTime == doctest::Approx(r.meanQueueLength / 0.9).epsilon(1e-8));
}
