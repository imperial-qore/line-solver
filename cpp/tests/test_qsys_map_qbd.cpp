/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP-driven queues solved on the port's own QBD and M/G/1-type machinery:
 * MAP/M/1, MAP/M/c, MAP/D/1, MAP/PH/1, MAP/MAP/1 and PH/PH/1.
 *
 * The MATLAB references wrap Q-MAM (Q_CT_MAP_M_C, Q_CT_MAP_D_C) and BUTools
 * (MMAPPH1FCFS); the port does not transcribe either, so what is asserted here
 * is that the QUANTITIES agree, not that the algorithms do. Three kinds of
 * oracle appear below.
 *
 *  - Textbook collapses. A Poisson arrival MAP with exponential or
 *    deterministic service must reproduce M/M/1, Erlang-C and the
 *    Pollaczek-Khinchine M/D/1 formulas. These are exact and are asserted at
 *    1e-12 relative, a margin over the measured 1e-13..1e-15.
 *  - MATLAB values, recorded to full precision by running
 *      matlab -batch "lineStart; qsys_...(...)"
 *    on the inputs named in each test. Two tolerances are used, and which one
 *    applies is stated at every call site:
 *      * 1e-12 relative where the reference is BUTools-based (qsys_phph1,
 *        qsys_mapph1, qsys_mapmap1 with a renewal service). Measured agreement
 *        there is below 1e-15 on every metric.
 *      * 1e-8 relative where the reference is Q-MAM-based (qsys_mapm1,
 *        qsys_mapmc, and the meanQueueLength of qsys_mapd1). Q-MAM truncates
 *        its level vector at maxNumComp = 500, which costs it 1e-9..1e-10 of
 *        relative accuracy: on the M/M/1 and M/M/c instances below the port
 *        matches the textbook formula to 1e-13 while MATLAB is 2e-9 away from
 *        it, so the 1e-8 tolerance is bounding the REFERENCE's error and must
 *        not be tightened to chase MATLAB.
 *  - Internal identities that hold for every instance: Little's law
 *    L = lambda W, W = Wq + E[S], utilization = lambda E[S]/c, the
 *    queue-length distribution summing to 1, and for MAP/D/1 the additional
 *    exact identity p_0 = 1 - rho.
 *
 * qsys_mapmap1's correlated-service case is deliberately NOT asserted against
 * qsys_mapmap1.m: the reference replaces the service MAP by its embedded PH and
 * so drops the service correlation. It is asserted against LINE's own MATLAB
 * qbd_mapmap1 instead, which solves the genuine MAP/MAP/1 queue.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mapd1.h"
#include "line/api/qsys/qsys_mapm1.h"
#include "line/api/qsys/qsys_mapmap1.h"
#include "line/api/qsys/qsys_mapmc.h"
#include "line/api/qsys/qsys_mapph1.h"
#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_phph1.h"

using line::Matrix;
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

/** Assemble a square matrix from a row-major literal. */
template <class T>
Matrix<T> mkmat(const std::vector<std::vector<double>>& a) {
    const std::size_t n = a.size();
    Matrix<T> M(n, a[0].size());
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j)
            M(i, j) = line::num_traits<T>::from_double(a[i][j]);
    return M;
}

/** Poisson arrival MAP of rate lambda. */
template <class T>
Map<T> poisson(double lambda) {
    return mkmap<T>({{-lambda}}, {{lambda}});
}

/** The correlated MMPP2 used throughout, lambda = 7/6. */
template <class T>
Map<T> corr_map() {
    return mkmap<T>({{-2.5, 0.2}, {0.1, -0.7}}, {{2.3, 0.0}, {0.0, 0.6}});
}

/** Erlang-2 renewal arrival MAP with phase rate 4, lambda = 2. */
template <class T>
Map<T> erlang2_map() {
    return mkmap<T>({{-4.0, 4.0}, {0.0, -4.0}}, {{0.0, 0.0}, {4.0, 0.0}});
}

/** Total mass of a queue-length distribution, as a double. */
template <class T>
double dist_mass(const std::vector<T>& p) {
    double t = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) t += line::num_traits<T>::to_double(p[i]);
    return t;
}

}  // namespace

// ---------------------------------------------------------------------------
// MAP/MAP/1, MAP/PH/1, PH/PH/1
// ---------------------------------------------------------------------------

TEST_CASE("MAP/MAP/1 with Poisson arrivals and exponential service is M/M/1") {
    // lambda = 2, mu = 3, rho = 2/3, L = 2, W = 1, Wq = 2/3. Exact oracle, so
    // the tolerance is the textbook one, 1e-12.
    const auto r = line::qsys::qsys_mapmap1(poisson<double>(2.0), poisson<double>(3.0));
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
    // MATLAB qsys_mapmap1([-2],[2],[-3],[3]) reports 1.999999999999999 and
    // 0.9999999999999989; BUTools-based, so 1e-12 relative.
    CHECK(r.meanQueueLength == doctest::Approx(1.999999999999999).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(0.9999999999999989).epsilon(1e-12));
    // Little's law and W = Wq + E[S] with E[S] = 1/3.
    CHECK(r.meanQueueLength == doctest::Approx(2.0 * r.meanSojournTime).epsilon(1e-13));
    CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0 / 3.0).epsilon(1e-13));
}

TEST_CASE("MAP/MAP/1 with a renewal service reproduces the BUTools reference") {
    // Correlated MMPP2 arrivals, Erlang-2 renewal service of mean 1/3. The
    // reference's PH reduction of a renewal service is exact, so the two model
    // the same system; BUTools-based, 1e-12 relative.
    const Map<double> svc = mkmap<double>({{-6.0, 6.0}, {0.0, -6.0}}, {{0.0, 0.0}, {6.0, 0.0}});
    const auto r = line::qsys::qsys_mapmap1(corr_map<double>(), svc);
    CHECK(r.meanQueueLength == doctest::Approx(0.8364637529649164).epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(0.3836355977794533).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(0.7169689311127866).epsilon(1e-12));
    CHECK(r.utilization == doctest::Approx(0.3888888888888888).epsilon(1e-14));
    // Little's law with lambda = 7/6 and E[S] = 1/3.
    CHECK(r.meanQueueLength == doctest::Approx((7.0 / 6.0) * r.meanSojournTime).epsilon(1e-13));
    CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0 / 3.0).epsilon(1e-13));
    CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-8));
}

TEST_CASE("MAP/MAP/1 honours the service correlation the reference discards") {
    // Correlated arrivals AND a correlated service MAP (lambda_s = 7/3, so
    // rho = 1/2). qsys_mapmap1.m replaces the service MAP by its embedded PH,
    // dropping the correlation, and reports 2.0305355320573; that value is NOT
    // asserted here because it answers a different question. The oracle is
    // LINE's own MATLAB qbd_mapmap1 on the same pair of MAPs, which solves the
    // genuine MAP/MAP/1 queue and returns QN = 2.547554867631855,
    // UN = 0.5000000000000002, RN = 2.183618457970162. Tolerance 1e-5
    // relative: MATLAB's qbd_mapmap1 sums the level distribution only until the
    // accumulated mass reaches 1 - 1e-10 and drops the tail times its level
    // index, which costs it 1.05e-6 here, while this port sums the geometric
    // tail in closed form.
    const Map<double> svc = mkmap<double>({{-5.0, 0.4}, {0.2, -1.4}}, {{4.6, 0.0}, {0.0, 1.2}});
    const auto r = line::qsys::qsys_mapmap1(corr_map<double>(), svc, 4000);
    CHECK(r.meanQueueLength == doctest::Approx(2.547554867631855).epsilon(1e-5));
    CHECK(r.meanSojournTime == doctest::Approx(2.183618457970162).epsilon(1e-5));
    CHECK(r.utilization == doctest::Approx(0.5).epsilon(1e-13));
    // The discarded correlation is worth 20% of the mean queue length, so the
    // port must sit well away from the reference's PH-reduced answer.
    CHECK(r.meanQueueLength > 1.15 * 2.0305355320573);
    CHECK(r.meanQueueLength == doctest::Approx((7.0 / 6.0) * r.meanSojournTime).epsilon(1e-13));
    CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("MAP/PH/1 matches the BUTools reference") {
    // BUTools-based reference throughout, 1e-12 relative.
    SUBCASE("Poisson arrivals and exponential service collapse to M/M/1") {
        const auto r = line::qsys::qsys_mapph1(poisson<double>(2.0), std::vector<double>{1.0},
                                               mkmat<double>({{-3.0}}));
        CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
    }
    SUBCASE("correlated arrivals, Erlang-2 service") {
        // MATLAB qsys_mapph1(D0c, D1c, [1 0], [-6 6; 0 -6]).
        const auto r = line::qsys::qsys_mapph1(corr_map<double>(), std::vector<double>{1.0, 0.0},
                                               mkmat<double>({{-6.0, 6.0}, {0.0, -6.0}}));
        CHECK(r.meanQueueLength == doctest::Approx(0.8364637529649164).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(0.3836355977794533).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(0.7169689311127866).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(0.3888888888888888).epsilon(1e-14));
    }
    SUBCASE("Erlang-2 arrivals, hyperexponential service") {
        // MATLAB qsys_mapph1(D0e2, D1e2, [0.6 0.4], diag(-8, -1.6)). E[S] =
        // 0.6/8 + 0.4/1.6 = 0.325, lambda = 2, so rho = 0.65.
        const auto r = line::qsys::qsys_mapph1(erlang2_map<double>(),
                                               std::vector<double>{0.6, 0.4},
                                               mkmat<double>({{-8.0, 0.0}, {0.0, -1.6}}));
        CHECK(r.meanQueueLength == doctest::Approx(2.133647429543165).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(0.7418237147715814).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(1.066823714771581).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(0.65).epsilon(1e-14));
        CHECK(r.meanQueueLength == doctest::Approx(2.0 * r.meanSojournTime).epsilon(1e-13));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 0.325).epsilon(1e-13));
    }
}

TEST_CASE("PH/PH/1 matches the BUTools reference") {
    // BUTools-based reference throughout, 1e-12 relative.
    SUBCASE("exponential over exponential is M/M/1") {
        const auto r = line::qsys::qsys_phph1(std::vector<double>{1.0}, mkmat<double>({{-2.0}}),
                                              std::vector<double>{1.0}, mkmat<double>({{-3.0}}));
        CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
    }
    SUBCASE("Erlang-2 over Erlang-2") {
        // MATLAB qsys_phph1([1 0], [-4 4; 0 -4], [1 0], [-6 6; 0 -6]).
        const auto r = line::qsys::qsys_phph1(std::vector<double>{1.0, 0.0},
                                              mkmat<double>({{-4.0, 4.0}, {0.0, -4.0}}),
                                              std::vector<double>{1.0, 0.0},
                                              mkmat<double>({{-6.0, 6.0}, {0.0, -6.0}}));
        CHECK(r.meanQueueLength == doctest::Approx(1.250000000000001).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(0.2916666666666675).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(0.6250000000000008).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
        // Both processes are less variable than Poisson, so E2/E2/1 must queue
        // less than the M/M/1 at the same rates.
        CHECK(r.meanQueueLength < 2.0);
    }
    SUBCASE("hyperexponential arrivals over exponential service") {
        // MATLAB qsys_phph1([0.7 0.3], diag(-5,-1), [1], [-10/3]). Mean
        // interarrival 0.44, lambda = 25/11, E[S] = 0.3, rho = 15/22.
        const auto r = line::qsys::qsys_phph1(std::vector<double>{0.7, 0.3},
                                              mkmat<double>({{-5.0, 0.0}, {0.0, -1.0}}),
                                              std::vector<double>{1.0},
                                              mkmat<double>({{-1.0 / 0.3}}));
        CHECK(r.meanQueueLength == doctest::Approx(3.273624198148769).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(1.140394647185455).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(1.440394647185455).epsilon(1e-12));
        CHECK(r.utilization == doctest::Approx(15.0 / 22.0).epsilon(1e-13));
        CHECK(r.meanQueueLength ==
              doctest::Approx((25.0 / 11.0) * r.meanSojournTime).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 0.3).epsilon(1e-13));
        // Hyperexponential arrivals are burstier than Poisson at the same rate.
        CHECK(r.meanQueueLength > (15.0 / 22.0) / (1.0 - 15.0 / 22.0));
    }
}

// ---------------------------------------------------------------------------
// MAP/M/1 and MAP/M/c
// ---------------------------------------------------------------------------

TEST_CASE("MAP/M/1 with Poisson arrivals reproduces M/M/1 exactly") {
    // lambda = 2, mu = 3. The exact oracle is the textbook formula, asserted at
    // 1e-12. MATLAB's Q-MAM value 1.999999994584355 is 2.7e-9 away from it, so
    // it is asserted only at the 1e-8 Q-MAM tolerance.
    const auto r = line::qsys::qsys_mapm1(poisson<double>(2.0), 3.0);
    CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(1.999999994584355).epsilon(1e-8));
    CHECK(r.meanWaitingTime == doctest::Approx(0.6666666665442794).epsilon(1e-8));
    // The level distribution is the M/M/1 geometric (1-rho) rho^n.
    CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(r.queueLengthDist[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(r.queueLengthDist[3] == doctest::Approx((1.0 / 3.0) * std::pow(2.0 / 3.0, 3.0))
                                      .epsilon(1e-12));
}

TEST_CASE("MAP/M/1 with non-Poisson arrivals matches MATLAB and MAP/MAP/1") {
    // Q-MAM-based reference, so 1e-8 relative; the port is cross-checked
    // against qsys_mapmap1 with an exponential service MAP, an independent
    // route (Kronecker QBD plus the closed-form factorial moment), at 1e-12.
    SUBCASE("correlated MMPP2 arrivals, mu = 2") {
        const auto r = line::qsys::qsys_mapm1(corr_map<double>(), 2.0);
        CHECK(r.meanQueueLength == doctest::Approx(2.795385119831413).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(1.896044396731732).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(2.396044396731732).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(7.0 / 12.0).epsilon(1e-14));
        const auto q = line::qsys::qsys_mapmap1(corr_map<double>(), poisson<double>(2.0), 4000);
        CHECK(r.meanQueueLength == doctest::Approx(q.meanQueueLength).epsilon(1e-12));
        CHECK(r.meanWaitingTime == doctest::Approx(q.meanWaitingTime).epsilon(1e-12));
        CHECK(r.meanQueueLength == doctest::Approx((7.0 / 6.0) * r.meanSojournTime).epsilon(1e-13));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 0.5).epsilon(1e-13));
    }
    SUBCASE("Erlang-2 arrivals, mu = 3") {
        const auto r = line::qsys::qsys_mapm1(erlang2_map<double>(), 3.0);
        CHECK(r.meanQueueLength == doctest::Approx(1.568729300352212).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.4510313188166638).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(0.7843646521499972).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
        const auto q = line::qsys::qsys_mapmap1(erlang2_map<double>(), poisson<double>(3.0), 4000);
        CHECK(r.meanQueueLength == doctest::Approx(q.meanQueueLength).epsilon(1e-12));
        // Erlang-2 arrivals are less bursty than Poisson at lambda = 2.
        CHECK(r.meanQueueLength < 2.0);
    }
}

TEST_CASE("MAP/M/c with Poisson arrivals reproduces Erlang-C exactly") {
    // lambda = 1.2, mu = 1, c = 2. Erlang-C: P0 = 1/4, Lq = 0.675, L = 1.875,
    // Wq = 0.5625, W = 1.5625. Exact oracle at 1e-12; MATLAB's Q-MAM value
    // 1.874999996296936 is 2.0e-9 away and is asserted only at 1e-8.
    const auto r = line::qsys::qsys_mapmc(poisson<double>(1.2), 1.0, 2u);
    CHECK(r.utilization == doctest::Approx(0.6).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(1.875).epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(0.5625).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(1.5625).epsilon(1e-12));
    CHECK(r.queueLengthDist[0] == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-13));
    CHECK(r.meanQueueLength == doctest::Approx(1.874999996296936).epsilon(1e-8));
    CHECK(r.meanWaitingTime == doctest::Approx(0.5624999998489539).epsilon(1e-8));
    // The port must agree with the port's own M/M/c closed form.
    CHECK(r.meanSojournTime == doctest::Approx(line::qsys::qsys_mmk(1.2, 1.0, 2u).W).epsilon(1e-12));
}

TEST_CASE("MAP/M/c with non-Poisson arrivals matches MATLAB") {
    // Q-MAM-based reference, 1e-8 relative.
    SUBCASE("correlated MMPP2 arrivals, mu = 1, c = 2") {
        const auto r = line::qsys::qsys_mapmc(corr_map<double>(), 1.0, 2u);
        CHECK(r.meanQueueLength == doctest::Approx(3.121850390835039).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(1.675871772757968).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(2.675871772757968).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(7.0 / 12.0).epsilon(1e-14));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0).epsilon(1e-13));
        // Little's law on the queued jobs: L = Lq + lambda E[S] with c = 2.
        CHECK(r.meanQueueLength ==
              doctest::Approx((7.0 / 6.0) * r.meanWaitingTime + (7.0 / 6.0)).epsilon(1e-12));
        // Correlated arrivals must queue more than Poisson at the same rate.
        CHECK(r.meanSojournTime > line::qsys::qsys_mmk(7.0 / 6.0, 1.0, 2u).W);
    }
    SUBCASE("Erlang-2 arrivals, mu = 1.5, c = 2") {
        const auto r = line::qsys::qsys_mapmc(erlang2_map<double>(), 1.5, 2u);
        CHECK(r.meanQueueLength == doctest::Approx(2.022056024561071).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(0.3443613471500779).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(1.011028013816745).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-14));
        CHECK(r.meanQueueLength == doctest::Approx(2.0 * r.meanSojournTime).epsilon(1e-12));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0 / 1.5).epsilon(1e-13));
    }
    SUBCASE("three servers, correlated arrivals, mu = 0.5") {
        const auto r = line::qsys::qsys_mapmc(corr_map<double>(), 0.5, 3u, 4000);
        CHECK(r.meanQueueLength == doctest::Approx(9.441571472754669).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(6.092775552492106).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(8.092775552492107).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(7.0 / 9.0).epsilon(1e-14));
        CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("MAP/M/c rejects an unstable or degenerate instance") {
    CHECK_THROWS_AS(line::qsys::qsys_mapmc(poisson<double>(3.0), 1.0, 2u), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mapmc(poisson<double>(1.0), -1.0, 1u), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mapm1(poisson<double>(4.0), 3.0), line::InputError);
}

// ---------------------------------------------------------------------------
// MAP/D/1
// ---------------------------------------------------------------------------

TEST_CASE("MAP/D/1 with Poisson arrivals reproduces the M/D/1 formulas") {
    // lambda = 2, s = 1/3, rho = 2/3. Pollaczek-Khinchine:
    // L = rho + rho^2/(2(1-rho)) = 4/3, Wq = rho s/(2(1-rho)) = 1/3. Exact
    // oracle, asserted at 1e-12; the port measures 2e-13 against it. MATLAB's
    // Q-MAM meanQueueLength 1.333333330345585 is 2.2e-9 away from the textbook
    // value, so it gets the 1e-8 Q-MAM tolerance.
    const auto r = line::qsys::qsys_mapd1(poisson<double>(2.0), 1.0 / 3.0);
    const double rho = 2.0 / 3.0;
    CHECK(r.utilization == doctest::Approx(rho).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(rho + rho * rho / (2.0 * (1.0 - rho)))
                                   .epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(rho * (1.0 / 3.0) / (2.0 * (1.0 - rho)))
                                   .epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0 / 3.0).epsilon(1e-13));
    CHECK(r.meanQueueLength == doctest::Approx(1.333333330345585).epsilon(1e-8));
    // Deterministic service is the least variable service at a given mean, so
    // M/D/1 must queue strictly less than the M/M/1 at the same rho.
    CHECK(r.meanQueueLength < rho / (1.0 - rho));
    // Identities of the construction: the distribution is a probability vector
    // and the empty-system probability is exactly 1 - rho.
    CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.queueLengthDist[0] == doctest::Approx(1.0 - rho).epsilon(1e-13));
    // MATLAB's meanWaitingTime is NOT asserted: qsys_mapdc integrates the
    // waiting-time CDF with a left-rectangle rule of step s/numSteps and
    // numSteps defaulting to 1, returning 0.4444444443040794 for the exact
    // 1/3. See the defect note in qsys_mapd1.h.
}

TEST_CASE("MAP/D/1 with non-Poisson arrivals matches the MATLAB queue length") {
    // Q-MAM-based reference for meanQueueLength, 1e-8 relative. Its
    // meanWaitingTime is defective (see qsys_mapd1.h) and is not asserted; the
    // port's own Wq is checked against Little's law instead.
    SUBCASE("correlated MMPP2 arrivals, s = 0.4") {
        const auto r = line::qsys::qsys_mapd1(corr_map<double>(), 0.4);
        const double lambda = 7.0 / 6.0, rho = lambda * 0.4;
        CHECK(r.meanQueueLength == doctest::Approx(1.118300079597477).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(rho).epsilon(1e-14));
        CHECK(r.meanWaitingTime ==
              doctest::Approx((r.meanQueueLength - rho) / lambda).epsilon(1e-14));
        CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 0.4).epsilon(1e-14));
        CHECK(r.meanQueueLength == doctest::Approx(lambda * r.meanSojournTime).epsilon(1e-13));
        CHECK(dist_mass(r.queueLengthDist) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(r.queueLengthDist[0] == doctest::Approx(1.0 - rho).epsilon(1e-13));
    }
    SUBCASE("Erlang-2 arrivals, s = 0.3") {
        const auto r = line::qsys::qsys_mapd1(erlang2_map<double>(), 0.3);
        const double lambda = 2.0, rho = 0.6;
        CHECK(r.meanQueueLength == doctest::Approx(0.7758216464528508).epsilon(1e-8));
        CHECK(r.utilization == doctest::Approx(rho).epsilon(1e-14));
        CHECK(r.meanQueueLength == doctest::Approx(lambda * r.meanSojournTime).epsilon(1e-13));
        CHECK(r.queueLengthDist[0] == doctest::Approx(1.0 - rho).epsilon(1e-13));
        // Erlang-2/D/1 is less variable on both sides than M/D/1 at rho = 0.6,
        // whose mean queue length is 0.6 + 0.36/0.8 = 1.05.
        CHECK(r.meanQueueLength < 1.05);
    }
}

TEST_CASE("MAP/D/1 rejects an unstable or degenerate instance") {
    CHECK_THROWS_AS(line::qsys::qsys_mapd1(poisson<double>(2.0), 0.6), line::InputError);  // rho>1
    CHECK_THROWS_AS(line::qsys::qsys_mapd1(poisson<double>(2.0), -0.1), line::InputError);
}

// ---------------------------------------------------------------------------
// Real50 instantiation
// ---------------------------------------------------------------------------

TEST_CASE("the MAP queue family instantiates and agrees at Real50") {
    // Extended precision must reproduce the double results to within double's
    // own accuracy, which shows that the double values are not
    // precision-limited. Measured agreement is 2e-14; asserted at 1e-11.
    typedef line::Real50 R;
    const auto rd = line::qsys::qsys_mapm1(corr_map<double>(), 2.0);
    const auto rr = line::qsys::qsys_mapm1(corr_map<R>(), R(2));
    CHECK(line::num_traits<R>::to_double(rr.meanQueueLength) ==
          doctest::Approx(rd.meanQueueLength).epsilon(1e-11));
    CHECK(line::num_traits<R>::to_double(rr.meanWaitingTime) ==
          doctest::Approx(rd.meanWaitingTime).epsilon(1e-11));

    const auto dd = line::qsys::qsys_mapd1(corr_map<double>(), 0.4);
    const auto dr = line::qsys::qsys_mapd1(corr_map<R>(), R("0.4"));
    CHECK(line::num_traits<R>::to_double(dr.meanQueueLength) ==
          doctest::Approx(dd.meanQueueLength).epsilon(1e-11));
    CHECK(line::num_traits<R>::to_double(dr.utilization) ==
          doctest::Approx(7.0 / 6.0 * 0.4).epsilon(1e-13));

    const auto md = line::qsys::qsys_mapmap1(corr_map<double>(), poisson<double>(3.0));
    const auto mr = line::qsys::qsys_mapmap1(corr_map<R>(), poisson<R>(3.0));
    CHECK(line::num_traits<R>::to_double(mr.meanQueueLength) ==
          doctest::Approx(md.meanQueueLength).epsilon(1e-11));

    const auto pd = line::qsys::qsys_phph1(std::vector<double>{1.0, 0.0},
                                           mkmat<double>({{-4.0, 4.0}, {0.0, -4.0}}),
                                           std::vector<double>{1.0, 0.0},
                                           mkmat<double>({{-6.0, 6.0}, {0.0, -6.0}}));
    const auto pr = line::qsys::qsys_phph1(std::vector<R>{R(1), R(0)},
                                           mkmat<R>({{-4.0, 4.0}, {0.0, -4.0}}),
                                           std::vector<R>{R(1), R(0)},
                                           mkmat<R>({{-6.0, 6.0}, {0.0, -6.0}}));
    CHECK(line::num_traits<R>::to_double(pr.meanQueueLength) ==
          doctest::Approx(pd.meanQueueLength).epsilon(1e-11));
    CHECK(line::num_traits<R>::to_double(pr.meanQueueLength) ==
          doctest::Approx(1.25).epsilon(1e-11));
}
