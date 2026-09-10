/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * qsys_mdc_crommelin: M/D/c by the embedded chain at service-time multiples.
 * Ported to C++ on 2026-08-01; the JAR and native Python already had it.
 *
 * ORACLES.
 *  (a) M/D/1 is closed form. Pollaczek-Khinchine with SCV 0 gives
 *      Lq = rho^2 / (2(1-rho)), L = rho + Lq, and those are exact, not fitted.
 *  (b) qsys_mapdc with a Poisson arrival is an INDEPENDENT M/D/c algorithm
 *      (Q-MAM's Q_CT_MAP_D_C), so agreeing with it at c > 1 checks the multi-
 *      server branch against something that shares no code with this one.
 *  (c) Little's law inside the returned struct, which no truncation can fix if
 *      the stationary vector is wrong.
 */
#include <cmath>
#include <cstddef>

#include "doctest.h"
#include "line/api/qsys/qsys_mapdc.h"
#include "line/api/qsys/qsys_mdc_crommelin.h"
#include "line/api/mam/map_moment.h"

namespace qsys = line::qsys;
namespace mam = line::mam;

TEST_CASE("M/D/1 reproduces the Pollaczek-Khinchine closed form") {
    const double rhos[] = {0.25, 0.5, 0.8};
    for (std::size_t i = 0; i < 3; ++i) {
        const double rho = rhos[i], s = 1.0, lambda = rho / s;
        const qsys::MDcCrommelinResult<double> r = qsys::qsys_mdc_crommelin(lambda, s, 1u);
        const double Lq = rho * rho / (2.0 * (1.0 - rho));
        CHECK(r.utilization == doctest::Approx(rho).epsilon(1e-12));
        CHECK(r.meanWaitingQueue == doctest::Approx(Lq).epsilon(1e-8));
        CHECK(r.meanQueueLength == doctest::Approx(rho + Lq).epsilon(1e-8));
        CHECK(r.meanWaitingTime == doctest::Approx(Lq / lambda).epsilon(1e-8));
        CHECK(r.meanSojournTime == doctest::Approx(Lq / lambda + s).epsilon(1e-8));
    }
}

TEST_CASE("Little's law holds inside the returned struct") {
    const qsys::MDcCrommelinResult<double> r = qsys::qsys_mdc_crommelin(1.2, 0.5, 2u);
    // L = lambda W, with W the sojourn time; the queue-length mean and the
    // waiting mean are computed from the same pi but by different sums.
    CHECK(r.meanQueueLength == doctest::Approx(1.2 * r.meanSojournTime).epsilon(1e-8));
    CHECK(r.meanWaitingQueue == doctest::Approx(1.2 * r.meanWaitingTime).epsilon(1e-10));
}

TEST_CASE("M/D/c agrees with the independent Q-MAM MAP/D/c route") {
    // Poisson(lambda) as a one-phase MAP, so qsys_mapdc solves the same system
    // by an algorithm that shares no code with the embedded chain.
    const double lambda = 1.5, s = 1.0;
    const unsigned c = 2;
    mam::Map<double> arv;
    arv.D0 = line::Matrix<double>(1, 1, -lambda);
    arv.D1 = line::Matrix<double>(1, 1, lambda);

    const qsys::MDcCrommelinResult<double> a = qsys::qsys_mdc_crommelin(lambda, s, c);
    const qsys::MapDcResult<double> b = qsys::qsys_mapdc(arv, s, c);
    CHECK(a.meanQueueLength == doctest::Approx(b.meanQueueLength).epsilon(1e-4));
    CHECK(a.meanSojournTime == doctest::Approx(b.meanSojournTime).epsilon(1e-4));
}

TEST_CASE("an explicit truncation is honoured and refusals are by name") {
    const qsys::MDcCrommelinResult<double> t = qsys::qsys_mdc_crommelin(0.5, 1.0, 1u, 400);
    CHECK(t.meanQueueLength == doctest::Approx(0.75).epsilon(1e-8));

    CHECK_THROWS_AS(qsys::qsys_mdc_crommelin(0.0, 1.0, 1u), line::InputError);
    CHECK_THROWS_AS(qsys::qsys_mdc_crommelin(1.0, 0.0, 1u), line::InputError);
    CHECK_THROWS_AS(qsys::qsys_mdc_crommelin(1.0, 1.0, 0u), line::InputError);
    CHECK_THROWS_AS(qsys::qsys_mdc_crommelin(2.0, 1.0, 1u), line::InputError);  // rho >= 1
}
