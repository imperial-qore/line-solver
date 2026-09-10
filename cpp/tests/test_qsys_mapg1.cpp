/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The MAP/G/1 queue by PH moment matching. Oracles, in order of strength:
 *   1. A CLOSED FORM that does not depend on the fit at all. For POISSON
 *      arrivals the Pollaczek-Khinchine mean wait lambda m2 / (2(1 - rho))
 *      depends only on the first two service moments, so any fit that
 *      preserves them must reproduce it exactly. This tests the whole chain
 *      (branch selection, PH construction, QBD solution) without depending on
 *      which PH the fit chose.
 *   2. The moments of the constructed PH: the Erlang and hyperexponential
 *      branches are explicit constructions whose mean and cv2 are known.
 *   3. MATLAB qsys_mapg1, which routes through BUTools MMAPPH1FCFS, for the
 *      MAP arrival cases where no closed form exists.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapg1.h"

using line::Matrix;
using line::Real50;
using line::mam::Map;
using line::qsys::MapG1Result;
using line::qsys::MapG1ServiceFit;
using line::qsys::qsys_mapg1;
using line::qsys::qsys_mapg1_service_fit;

namespace {

/** Poisson arrivals of rate lambda, as a MAP. */
template <class T>
Map<T> poisson(const T& lambda) {
    Map<T> m;
    m.D0 = Matrix<T>(1, 1, T(-lambda));
    m.D1 = Matrix<T>(1, 1, lambda);
    return m;
}

/** The MMPP2 arrival of the MATLAB reference runs, lambda = 7/6. */
template <class T>
Map<T> mmpp2() {
    using nt = line::num_traits<T>;
    Map<T> m;
    m.D0 = Matrix<T>{{nt::from_double(-2.5), nt::from_double(0.2)},
                     {nt::from_double(0.1), nt::from_double(-0.7)}};
    m.D1 = Matrix<T>{{nt::from_double(2.3), nt::from_double(0)},
                     {nt::from_double(0), nt::from_double(0.6)}};
    return m;
}

/** Pollaczek-Khinchine mean waiting time of an M/G/1 queue. */
double pk_wait(double lambda, double m1, double m2) {
    return lambda * m2 / (2.0 * (1.0 - lambda * m1));
}

const double M1 = 1.0 / 3.0;
constexpr double TOL = 1e-12;

}  // namespace

TEST_CASE("qsys_mapg1 collapses to M/M/1 on one and on two moments") {
    // Exp(3) service: cv2 = 1 exactly, so both the 1-moment and the 2-moment
    // branch must produce the same exponential.
    const std::vector<double> one{M1};
    const std::vector<double> two{M1, 2.0 * M1 * M1};
    for (const std::vector<double>& mom : {one, two}) {
        const MapG1Result<double> r = qsys_mapg1(poisson(2.0), mom);
        CHECK(r.servicePhases == 1);
        CHECK(r.fitKind == MapG1ServiceFit::Exponential);
        // M/M/1 with rho = 2/3: L = 2, W = 1, Wq = 2/3. MATLAB agrees to 1e-15.
        CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(TOL));
        CHECK(r.meanSojournTime == doctest::Approx(1.0).epsilon(TOL));
        CHECK(r.meanWaitingTime == doctest::Approx(2.0 / 3.0).epsilon(TOL));
        CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(TOL));
    }
}

TEST_CASE("qsys_mapg1 reproduces Pollaczek-Khinchine on every service branch") {
    // The mean wait of an M/G/1 queue depends only on lambda, m1 and m2, so a
    // fit that preserves them must land on the closed form whatever PH it
    // picked. Each entry names the branch it exercises.
    struct Case {
        const char* name;
        double m1;
        double m2;
        std::size_t phases;
        MapG1ServiceFit kind;
    };
    const std::vector<Case> cases{
        {"cv2 = 1, exponential", M1, 2.0 * M1 * M1, 1, MapG1ServiceFit::Exponential},
        {"cv2 = 1/4, Erlang-4", 0.25, 1.25 * 0.0625, 4, MapG1ServiceFit::Erlang},
        {"cv2 = 4, hyperexponential", 0.25, 5.0 * 0.0625, 2, MapG1ServiceFit::Hyperexponential},
    };
    for (const Case& c : cases) {
        INFO(c.name);
        const MapG1Result<double> r = qsys_mapg1(poisson(2.0), {c.m1, c.m2});
        CHECK(r.servicePhases == c.phases);
        CHECK(r.fitKind == c.kind);
        CHECK(r.meanWaitingTime == doctest::Approx(pk_wait(2.0, c.m1, c.m2)).epsilon(1e-11));
        // and the fitted PH really does carry the requested two moments
        CHECK(line::mam::map_moment(r.serviceFit, 1) == doctest::Approx(c.m1).epsilon(TOL));
        CHECK(line::mam::map_moment(r.serviceFit, 2) == doctest::Approx(c.m2).epsilon(TOL));
    }
    // MATLAB qsys_mapg1 on the same three: 0.666666666666666, 0.15625,
    // 0.625000000000002.
    CHECK(qsys_mapg1(poisson(2.0), {0.25, 1.25 * 0.0625}).meanQueueLength ==
          doctest::Approx(0.812499999999996).epsilon(1e-11));
    CHECK(qsys_mapg1(poisson(2.0), {0.25, 5.0 * 0.0625}).meanQueueLength ==
          doctest::Approx(1.75).epsilon(1e-11));
}

TEST_CASE("qsys_mapg1 three-moment branch matches MATLAB on a MAP arrival") {
    // Interior moment set (cv2 = 0.6, m3 = 3.6 m1^3): the port's aph_fit and
    // MATLAB's agree here on the order (2) and on the sub-generator
    // ([-10.8541019662499 10.8541019662499; 0 -4.14589803375029]), so the queue
    // results agree too.
    const std::vector<double> mom{M1, 1.6 * M1 * M1, 3.6 * M1 * M1 * M1};
    const MapG1Result<double> r = qsys_mapg1(mmpp2<double>(), mom);
    CHECK(r.fitKind == MapG1ServiceFit::Acyclic);
    CHECK(r.servicePhases == 2);
    CHECK(r.serviceFit.D0(0, 0) == doctest::Approx(-10.8541019662499).epsilon(1e-12));
    CHECK(r.serviceFit.D0(1, 1) == doctest::Approx(-4.14589803375029).epsilon(1e-12));
    // MATLAB qsys_mapg1(A0, A1, mom)
    CHECK(r.meanQueueLength == doctest::Approx(0.858359337155326).epsilon(1e-10));
    CHECK(r.meanWaitingTime == doctest::Approx(0.402403241371232).epsilon(1e-10));
    CHECK(r.meanSojournTime == doctest::Approx(0.735736574704565).epsilon(1e-10));
    CHECK(r.utilization == doctest::Approx(0.388888888888889).epsilon(1e-12));
    // the queue length distribution, entry by entry
    const double qld[4] = {0.611111111111111, 0.19840954130397, 0.0814806347403443,
                           0.0432634808557709};
    REQUIRE(r.queueLengthDist.size() >= 4);
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(r.queueLengthDist[i] == doctest::Approx(qld[i]).epsilon(1e-10));

    // cv2 > 1 interior set, same agreement
    const std::vector<double> mom2{M1, 2.5 * M1 * M1, 12.0 * M1 * M1 * M1};
    const MapG1Result<double> r2 = qsys_mapg1(mmpp2<double>(), mom2);
    CHECK(r2.meanQueueLength == doctest::Approx(1.03782328887253).epsilon(1e-10));
    CHECK(r2.meanWaitingTime == doctest::Approx(0.556229485700265).epsilon(1e-10));

    // two-moment branch under a MAP arrival, where P-K does not apply
    const MapG1Result<double> r3 = qsys_mapg1(mmpp2<double>(), {0.25, 1.25 * 0.0625});
    CHECK(r3.meanQueueLength == doctest::Approx(0.448488770069416).epsilon(1e-10));
    CHECK(r3.meanSojournTime == doctest::Approx(0.384418945773785).epsilon(1e-10));
}

TEST_CASE("qsys_mapg1 recovers the exact law when the moments are on the APH(2) boundary") {
    // (m1, m2, m3) = (1/3, 1/6, 1/9) are the moments of Erlang(2) with mean
    // 1/3, and that third moment sits exactly ON the APH(2) lower bound. The
    // port's aph_fit returns the Erlang(2) itself, order 2 with both rates 6;
    // MATLAB returns an order-3 APH with the same three moments (see the
    // divergence reported for mam::aph_fit). Both are valid three-moment
    // matches, so this asserts the property rather than the representation.
    const std::vector<double> mom{M1, 1.5 * M1 * M1, 3.0 * M1 * M1 * M1};
    MapG1ServiceFit kind = MapG1ServiceFit::Exponential;
    const Map<double> fit = qsys_mapg1_service_fit(mom, kind);
    CHECK(kind == MapG1ServiceFit::Acyclic);
    CHECK(fit.order() == 2);
    CHECK(fit.D0(0, 0) == doctest::Approx(-6.0).epsilon(1e-12));
    CHECK(fit.D0(1, 1) == doctest::Approx(-6.0).epsilon(1e-12));
    for (unsigned k = 1; k <= 3; ++k) {
        INFO("moment ", k);
        CHECK(line::mam::map_moment(fit, k) == doctest::Approx(mom[k - 1]).epsilon(1e-12));
    }
    // Under POISSON arrivals the representation cannot matter: the mean wait
    // is P-K, and MATLAB's order-3 fit gives the same 0.5 / 1.66666666666667.
    const MapG1Result<double> r = qsys_mapg1(poisson(2.0), mom);
    CHECK(r.meanWaitingTime == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(r.meanQueueLength == doctest::Approx(1.66666666666667).epsilon(1e-11));
}

TEST_CASE("qsys_mapg1 reproduces the Erlang-100 fit of the cv2 <= 0 branch") {
    // REFERENCE DEFECT, qsys_mapg1.m line 128: k = max(1, round(1/max(cv2,
    // 0.01))) pins the denominator at 0.01 for every cv2 <= 0, so a
    // deterministic service silently becomes an Erlang-100 with cv2 = 0.01.
    // MATLAB reproduction: qsys_mapg1(-2, 2, [1/3, (1/3)^2]) returns
    // meanQueueLength 1.33999999999909, not the M/D/1 value.
    const MapG1Result<double> r = qsys_mapg1(poisson(2.0), {M1, M1 * M1});
    CHECK(r.servicePhases == 100);
    CHECK(r.fitKind == MapG1ServiceFit::Erlang);
    // The fitted law has cv2 = 1/100, not 0.
    const double m1 = line::mam::map_moment(r.serviceFit, 1);
    const double m2 = line::mam::map_moment(r.serviceFit, 2);
    CHECK(m1 == doctest::Approx(M1).epsilon(1e-12));
    CHECK(m2 / (m1 * m1) - 1.0 == doctest::Approx(0.01).epsilon(1e-10));
    // P-K on the law that is actually solved, not on the requested moments
    CHECK(r.meanWaitingTime == doctest::Approx(pk_wait(2.0, m1, m2)).epsilon(1e-9));
    // MATLAB 1.33999999999909 through BUTools; the port sums the QBD tail in
    // closed form and gives 1.33999999999998. Same quantity, 7e-13 apart.
    CHECK(r.meanQueueLength == doctest::Approx(1.34).epsilon(1e-9));
    // and the exact M/D/1 value, which neither computes, is 4/3
    CHECK(std::abs(r.meanQueueLength - 4.0 / 3.0) > 1e-4);
}

TEST_CASE("qsys_mapg1 rejects malformed input and runs at high precision") {
    CHECK_THROWS_AS(qsys_mapg1(poisson(2.0), std::vector<double>{}), line::InputError);
    CHECK_THROWS_AS(qsys_mapg1(poisson(2.0), std::vector<double>{0.0}), line::InputError);

    // Dyadic moments, so cv2 is exactly 1 in every arithmetic: m1 = 1/2,
    // m2 = 1/2 is Exp(2), and with lambda = 1 the queue is M/M/1 at rho = 1/2.
    const MapG1Result<Real50> r = qsys_mapg1(
        poisson(Real50(1)), std::vector<Real50>{Real50(0.5), Real50(0.5)});
    CHECK(r.servicePhases == 1);
    CHECK(static_cast<double>(r.meanQueueLength) == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(static_cast<double>(r.meanWaitingTime) == doctest::Approx(0.5).epsilon(1e-14));
}

TEST_CASE("qsys_mapg1 cv2 == 1 is a float equality and is precision dependent") {
    // REFERENCE DEFECT, qsys_mapg1.m line 137: the exponential branch is
    // selected by `cv2 == 1`, an exact comparison on a computed quantity.
    // With the moments of Exp(3) written as (1/3, 2/9) the quotient rounds to
    // exactly 1 in IEEE double, so MATLAB takes the exponential branch and
    // returns the M/M/1 answer; at 50 digits the same moments give a cv2 just
    // below 1, the cv2 < 1 branch is taken, k = max(2, round(1/cv2)) = 2, and
    // the service becomes an ERLANG-2. The queue answer moves from L = 2 to
    // L = 5/3, a 17 per cent step from a rounding decision. Reproduced rather
    // than patched: a tolerance here would be a modelling choice the reference
    // does not make.
    const std::vector<double> d{M1, 2.0 * M1 * M1};
    CHECK(qsys_mapg1(poisson(2.0), d).servicePhases == 1);  // double: exponential
    const std::vector<Real50> r50{Real50(1) / Real50(3), Real50(2) / Real50(9)};
    const MapG1Result<Real50> r = qsys_mapg1(poisson(Real50(2)), r50);
    CHECK(r.servicePhases == 2);  // Real50: Erlang-2
    CHECK(r.fitKind == MapG1ServiceFit::Erlang);
    CHECK(static_cast<double>(r.meanQueueLength) == doctest::Approx(5.0 / 3.0).epsilon(1e-12));
}
