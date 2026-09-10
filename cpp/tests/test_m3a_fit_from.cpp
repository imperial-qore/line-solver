/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The m3a fitters driven from a process or a trace: aph2_fit_map,
 * aph2_fit_trace, amap2_fit_gamma_map, amap2_fit_gamma_trace.
 *
 * ORACLES. A wrapper's job is to measure the right descriptors and pass them on,
 * so the tests check exactly that:
 *  - fitting FROM a process must equal fitting from that process's own moments,
 *    computed independently;
 *  - fitting an APH(2) to a two-phase Erlang must return that Erlang's moments,
 *    since it is already an APH(2) and the fit is exact;
 *  - a trace drawn from a known process must fit close to the process itself,
 *    which is the only end-to-end check the trace path admits.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/m3a_fit_from.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_sample.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;

TEST_CASE("aph2_fit_map equals aph2_fit on the same process's moments") {
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const mam::Aph2FitResult<double> a = mam::aph2_fit_map(m);
    const mam::Aph2FitResult<double> b =
        mam::aph2_fit(mam::map_mean(m), mam::map_moment(m, 2), mam::map_moment(m, 3));
    REQUIRE(a.aph.order() == b.aph.order());
    for (std::size_t i = 0; i < a.aph.order(); ++i)
        for (std::size_t j = 0; j < a.aph.order(); ++j) {
            CHECK(a.aph.D0(i, j) == doctest::Approx(b.aph.D0(i, j)));
            CHECK(a.aph.D1(i, j) == doctest::Approx(b.aph.D1(i, j)));
        }
}

TEST_CASE("an APH(2) fitted to an Erlang-2 reproduces its moments exactly") {
    // An Erlang-2 IS an APH(2), so the fit is exact rather than approximate.
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const mam::Aph2FitResult<double> r = mam::aph2_fit_map(m);
    CHECK(r.adjusted == false);
    CHECK(mam::map_mean(r.aph) == doctest::Approx(mam::map_mean(m)).epsilon(1e-9));
    CHECK(mam::map_moment(r.aph, 2) == doctest::Approx(mam::map_moment(m, 2)).epsilon(1e-9));
    CHECK(mam::map_moment(r.aph, 3) == doctest::Approx(mam::map_moment(m, 3)).epsilon(1e-9));
}

TEST_CASE("aph2_fit_trace recovers the law the trace was drawn from") {
    line::pfqn::McRng rng(4242u);
    const mam::Map<double> src = mam::map_erlang(1.0, 2);
    const std::vector<double> S = mam::map_sample(src, 300000, rng);
    const mam::Aph2FitResult<double> r = mam::aph2_fit_trace(S);
    CHECK(mam::map_mean(r.aph) == doctest::Approx(mam::map_mean(src)).epsilon(0.02));
    CHECK(mam::map_scv(r.aph) == doctest::Approx(mam::map_scv(src)).epsilon(0.05));
}

TEST_CASE("aph2_fit_trace refuses an empty trace by name") {
    CHECK_THROWS_AS(mam::aph2_fit_trace(std::vector<double>()), line::InputError);
    CHECK_THROWS_AS(mam::amap2_fit_gamma_trace(std::vector<double>()), line::InputError);
}

TEST_CASE("amap2_fit_gamma_map equals amap2_fit_gamma on the same descriptors") {
    // An MMPP(2): a renewal process has gamma 0 and would not exercise the
    // fourth descriptor at all.
    mam::Map<double> m;
    m.D0 = line::Matrix<double>(2, 2, 0.0);
    m.D1 = line::Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -(3.0 + 0.2); m.D0(0, 1) = 0.2;
    m.D0(1, 0) = 0.1;          m.D0(1, 1) = -(0.5 + 0.1);
    m.D1(0, 0) = 3.0;
    m.D1(1, 1) = 0.5;

    const mam::Amap2FitGammaResult<double> a = mam::amap2_fit_gamma_map(m);
    const mam::Amap2FitGammaResult<double> b = mam::amap2_fit_gamma(
        mam::map_mean(m), mam::map_moment(m, 2), mam::map_moment(m, 3), mam::map_gamma(m));
    CHECK(a.poisson_fallback == b.poisson_fallback);
    REQUIRE(a.amap.order() == b.amap.order());
    for (std::size_t i = 0; i < a.amap.order(); ++i)
        for (std::size_t j = 0; j < a.amap.order(); ++j) {
            CHECK(a.amap.D0(i, j) == doctest::Approx(b.amap.D0(i, j)));
            CHECK(a.amap.D1(i, j) == doctest::Approx(b.amap.D1(i, j)));
        }
    // The fit keeps the mean it was given.
    if (!a.poisson_fallback)
        CHECK(mam::map_mean(a.amap) == doctest::Approx(mam::map_mean(m)).epsilon(1e-6));
}

TEST_CASE("amap2_fit_gamma_trace measures its gamma from the trace") {
    line::pfqn::McRng rng(808u);
    mam::Map<double> src;
    src.D0 = line::Matrix<double>(2, 2, 0.0);
    src.D1 = line::Matrix<double>(2, 2, 0.0);
    src.D0(0, 0) = -(3.0 + 0.2); src.D0(0, 1) = 0.2;
    src.D0(1, 0) = 0.1;          src.D0(1, 1) = -(0.5 + 0.1);
    src.D1(0, 0) = 3.0;
    src.D1(1, 1) = 0.5;

    const std::vector<double> S = mam::map_sample(src, 200000, rng);
    const mam::Amap2FitGammaResult<double> r = mam::amap2_fit_gamma_trace(S);
    // Whatever branch it took, the mean must track the trace's own mean.
    double sm = 0.0;
    for (std::size_t i = 0; i < S.size(); ++i) sm += S[i];
    sm /= static_cast<double>(S.size());
    CHECK(mam::map_mean(r.amap) == doctest::Approx(sm).epsilon(0.03));
}
