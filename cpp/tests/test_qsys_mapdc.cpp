/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP/D/c by Crommelin's embedded lattice chain (line/api/qsys/qsys_mapdc.h).
 *
 * Oracles, in the order the task prescribes.
 *  (a) Closed form. Poisson arrivals with c = 1 give the M/D/1 queue, whose
 *      Pollaczek-Khinchin mean queue length rho + rho^2/(2(1-rho)) and mean
 *      waiting time rho s/(2(1-rho)) are exact. On rho = 2/3 that is L = 4/3
 *      and Wq = 1/3 -- the value against which the reference's superseded
 *      one-point quadrature returned 0.4444444443, and which this port must
 *      reproduce.
 *  (b) Invariants of the construction, asserted on every instance:
 *      sum_n P(N = n) = 1, and sum_n min(n, c) P(N = n) = lambda s, i.e. the
 *      mean departures per interval equal the mean arrivals. The second is the
 *      c-server form of Little's law applied to the servers and is what pins
 *      the multiserver grouping: an off-by-one in the sub-level bookkeeping
 *      moves it immediately.
 *  (c) c = 1 against qsys_mapd1, which reaches the same numbers through the
 *      DEPARTURE-epoch chain, a different construction with a different
 *      truncation; then MATLAB's Q-MAM at its own accuracy.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mapd1.h"
#include "line/api/qsys/qsys_mapdc.h"

using line::Matrix;
using line::Real50;
using line::mam::Map;

namespace {

/** Assemble a MAP from row-major (D0, D1) literals. */
template <class T>
Map<T> mkmap(const std::vector<std::vector<double>>& d0,
             const std::vector<std::vector<double>>& d1) {
    const std::size_t n = d0.size();
    Map<T> m;
    m.D0 = Matrix<T>(n, n);
    m.D1 = Matrix<T>(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            m.D0(i, j) = line::num_traits<T>::from_double(d0[i][j]);
            m.D1(i, j) = line::num_traits<T>::from_double(d1[i][j]);
        }
    return m;
}

/** Poisson arrivals of rate lambda as a one-phase MAP. */
template <class T>
Map<T> poisson(double lambda) {
    return mkmap<T>({{-lambda}}, {{lambda}});
}

/** Correlated MMPP2, lambda = 7/6. */
template <class T>
Map<T> corr_map() {
    return mkmap<T>({{-2.5, 0.2}, {0.1, -0.7}}, {{2.3, 0.0}, {0.0, 0.6}});
}

/** Erlang-2 renewal arrivals, lambda = 2. */
template <class T>
Map<T> erlang2_map() {
    return mkmap<T>({{-4.0, 4.0}, {0.0, -4.0}}, {{0.0, 0.0}, {4.0, 0.0}});
}

double dist_mass(const std::vector<double>& p) {
    double s = 0.0;
    for (double v : p) s += v;
    return s;
}

/** sum_n min(n, c) P(N = n), the mean number of busy servers. */
double busy_servers(const std::vector<double>& p, unsigned c) {
    double s = 0.0;
    for (std::size_t n = 0; n < p.size(); ++n)
        s += static_cast<double>(std::min<std::size_t>(n, c)) * p[n];
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// (a) the M/D/1 closed form, which is also the reference-defect regression
// ---------------------------------------------------------------------------

TEST_CASE("MAP/D/c at c = 1 with Poisson arrivals reproduces Pollaczek-Khinchin") {
    // lambda = 2, s = 1/3, rho = 2/3. L = rho + rho^2/(2(1-rho)) = 4/3,
    // Wq = rho s/(2(1-rho)) = 1/3.
    const auto r = line::qsys::qsys_mapdc(poisson<double>(2.0), 1.0 / 3.0, 1u);
    const double rho = 2.0 / 3.0;
    CHECK(r.utilization == doctest::Approx(rho).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(4.0 / 3.0).epsilon(1e-11));
    CHECK(r.meanWaitingTime == doctest::Approx(1.0 / 3.0).epsilon(1e-11));
    CHECK(r.meanSojournTime == doctest::Approx(1.0 / 3.0 + 1.0 / 3.0).epsilon(1e-11));
    // The waiting time the reference used to return, from a one-point
    // left-rectangle quadrature, was 0.4444444443: +33%. It is now Little's
    // law on both sides and this assertion is what keeps it that way.
    CHECK(std::fabs(r.meanWaitingTime - 0.4444444443040794) > 0.1);
    // MATLAB after the fix: L = 1.33333333034558, Wq = 0.333333331839459, at
    // its own Q-MAM truncation accuracy of about 2e-9.
    CHECK(r.meanQueueLength == doctest::Approx(1.33333333034558).epsilon(1e-8));
    CHECK(r.meanWaitingTime == doctest::Approx(0.333333331839459).epsilon(1e-8));
    // Deterministic service is the least variable at a given mean, so M/D/1
    // must queue strictly less than M/M/1 at the same load.
    CHECK(r.meanQueueLength < rho / (1.0 - rho));
    CHECK(r.queueLengthDist[0] == doctest::Approx(1.0 - rho).epsilon(1e-11));
}

// ---------------------------------------------------------------------------
// (b) invariants
// ---------------------------------------------------------------------------

TEST_CASE("MAP/D/c satisfies its two construction identities on every instance") {
    struct Case {
        Map<double> arrival;
        double s;
        unsigned c;
        double lambda;
    };
    const std::vector<Case> cases = {{poisson<double>(2.0), 1.0 / 3.0, 1u, 2.0},
                                     {poisson<double>(2.0), 0.6, 2u, 2.0},
                                     {poisson<double>(2.0), 1.2, 3u, 2.0},
                                     {corr_map<double>(), 0.8, 2u, 7.0 / 6.0},
                                     {erlang2_map<double>(), 0.9, 2u, 2.0}};
    for (const Case& cs : cases) {
        INFO("c = ", cs.c, ", s = ", cs.s);
        const auto r = line::qsys::qsys_mapdc(cs.arrival, cs.s, cs.c, 400u, 4096u, 20000u, 1e-14);
        CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-10));
        // Mean busy servers = mean arrivals per interval = lambda s.
        CHECK(busy_servers(r.queueLengthDist, cs.c) ==
              doctest::Approx(cs.lambda * cs.s).epsilon(1e-9));
        // Utilization is that divided by c.
        CHECK(r.utilization == doctest::Approx(cs.lambda * cs.s / cs.c).epsilon(1e-14));
        // Little's law is how meanSojournTime is defined, so it must close.
        CHECK(r.meanQueueLength == doctest::Approx(cs.lambda * r.meanSojournTime).epsilon(1e-13));
        CHECK(r.meanWaitingTime > 0.0);
        for (double p : r.queueLengthDist) CHECK(p >= -1e-15);
    }
}

// ---------------------------------------------------------------------------
// (c) cross-checks: qsys_mapd1, then MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("MAP/D/c at c = 1 agrees with qsys_mapd1's departure-epoch chain") {
    // Two independent constructions of the same queue: the lattice chain here
    // and the M/G/1-type chain embedded at departures in qsys_mapd1. They
    // share only the counting matrices A_k, so agreement to 1e-10 exercises
    // everything downstream of those.
    SUBCASE("Poisson, s = 1/3") {
        const auto a = line::qsys::qsys_mapdc(poisson<double>(2.0), 1.0 / 3.0, 1u);
        const auto b = line::qsys::qsys_mapd1(poisson<double>(2.0), 1.0 / 3.0);
        CHECK(a.meanQueueLength == doctest::Approx(b.meanQueueLength).epsilon(1e-10));
        CHECK(a.meanWaitingTime == doctest::Approx(b.meanWaitingTime).epsilon(1e-10));
    }
    SUBCASE("correlated MMPP2, s = 0.4") {
        const auto a = line::qsys::qsys_mapdc(corr_map<double>(), 0.4, 1u);
        const auto b = line::qsys::qsys_mapd1(corr_map<double>(), 0.4);
        CHECK(a.meanQueueLength == doctest::Approx(b.meanQueueLength).epsilon(1e-10));
        CHECK(a.meanWaitingTime == doctest::Approx(b.meanWaitingTime).epsilon(1e-10));
        // and therefore with the MATLAB Q-MAM value quoted in qsys_mapd1.h.
        CHECK(a.meanQueueLength == doctest::Approx(1.118300079597477).epsilon(1e-8));
    }
    SUBCASE("Erlang-2 arrivals, s = 0.3") {
        const auto a = line::qsys::qsys_mapdc(erlang2_map<double>(), 0.3, 1u);
        const auto b = line::qsys::qsys_mapd1(erlang2_map<double>(), 0.3);
        CHECK(a.meanQueueLength == doctest::Approx(b.meanQueueLength).epsilon(1e-10));
        CHECK(a.meanQueueLength == doctest::Approx(0.7758216464528508).epsilon(1e-8));
    }
}

TEST_CASE("MAP/D/c agrees with the MATLAB reference for c > 1") {
    // Q-MAM truncates its queue-length vector at maxNumComp and loses about
    // 8e-11 of mass on these instances, so its numbers get a 1e-8 tolerance,
    // as they do for qsys_mapd1.
    SUBCASE("M/D/2, lambda = 2, s = 0.6, rho = 0.6") {
        const auto r = line::qsys::qsys_mapdc(poisson<double>(2.0), 0.6, 2u);
        CHECK(r.meanQueueLength == doctest::Approx(1.55164329134673).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.175821645673367).epsilon(1e-7));
        CHECK(r.meanSojournTime == doctest::Approx(0.775821645673367).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(0.6).epsilon(1e-14));
        CHECK(r.queueLengthDist[0] == doctest::Approx(0.238685365188794).epsilon(1e-8));
        CHECK(r.queueLengthDist[1] == doctest::Approx(0.322629269626472).epsilon(1e-8));
        CHECK(r.queueLengthDist[2] == doctest::Approx(0.231148685357196).epsilon(1e-8));
        // Deterministic service must wait strictly less than exponential
        // service at the same load. Erlang C at a = 1.2, c = 2 is 0.45, so
        // Wq(M/M/2) = 0.45/(2/0.6 - 2) = 0.3375.
        CHECK(r.meanWaitingTime < 0.3375);
        // and strictly more than half of it: M/D/c waiting is not exactly
        // half the M/M/c value away from heavy traffic.
        CHECK(r.meanWaitingTime > 0.5 * 0.3375 * 0.9);
    }
    SUBCASE("M/D/3, lambda = 2, s = 1.2, rho = 0.8") {
        const auto r = line::qsys::qsys_mapdc(poisson<double>(2.0), 1.2, 3u);
        CHECK(r.meanQueueLength == doctest::Approx(3.72936053685491).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.664680268427454).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(0.8).epsilon(1e-14));
        CHECK(r.queueLengthDist[0] == doctest::Approx(0.0498414304563823).epsilon(1e-8));
        CHECK(r.queueLengthDist[3] == doctest::Approx(0.181953124332401).epsilon(1e-8));
    }
    SUBCASE("MMPP2/D/2, s = 0.8") {
        const auto r = line::qsys::qsys_mapdc(corr_map<double>(), 0.8, 2u);
        CHECK(r.meanQueueLength == doctest::Approx(1.50708585214034).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.491787873263148).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(7.0 / 6.0 * 0.8 / 2.0).epsilon(1e-14));
        CHECK(r.queueLengthDist[0] == doctest::Approx(0.40302492018309).epsilon(1e-8));
    }
    SUBCASE("Erlang-2/D/2, s = 0.9, rho = 0.9") {
        const auto r = line::qsys::qsys_mapdc(erlang2_map<double>(), 0.9, 2u);
        CHECK(r.meanQueueLength == doctest::Approx(3.59993858928319).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.899969294641594).epsilon(1e-8));
        CHECK(r.queueLengthDist[0] == doctest::Approx(0.0276440900540898).epsilon(1e-7));
    }
}

TEST_CASE("MAP/D/c rejects an unstable or degenerate instance") {
    CHECK_THROWS_AS(line::qsys::qsys_mapdc(poisson<double>(2.0), 1.1, 2u), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mapdc(poisson<double>(2.0), -0.1, 2u), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mapdc(poisson<double>(2.0), 0.3, 0u), line::InputError);
}

TEST_CASE("MAP/D/c at Real50 reproduces its own double result") {
    const auto rd = line::qsys::qsys_mapdc(poisson<double>(2.0), 0.6, 2u);
    const auto rq = line::qsys::qsys_mapdc(poisson<Real50>(2.0), Real50(0.6), 2u);
    CHECK(static_cast<double>(rq.meanQueueLength) ==
          doctest::Approx(rd.meanQueueLength).epsilon(1e-12));
    CHECK(static_cast<double>(rq.meanWaitingTime) ==
          doctest::Approx(rd.meanWaitingTime).epsilon(1e-11));
}
