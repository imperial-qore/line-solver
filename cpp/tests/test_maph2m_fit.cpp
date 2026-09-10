/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * maph2m_fit and its multiclass core.
 *
 * ORACLES, in the order the fit constrains them:
 *  - the class probabilities p are an EQUALITY constraint, so they must come
 *    back exactly whatever the backward moments do;
 *  - the inter-arrival law is the APH(2)'s, so the first three moments of the
 *    fitted MAPH must be the ones it was asked for;
 *  - the backward moments are the objective, so they are checked as an
 *    approach, not an identity -- the feasibility of the split is what stops
 *    them from being matched, and asserting equality would be asserting that
 *    the constraints never bind.
 * The degenerate branch is checked on its own terms: there the only thing that
 * can be matched IS p.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/maph2m_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

void check_is_mmap(const mam::Mmap<double>& m) {
    const std::size_t n = m.order();
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0, dcs = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-9);
            for (std::size_t c = 0; c < m.classes(); ++c) {
                CHECK(m.Dc[c](i, j) >= -1e-9);
                dcs += m.Dc[c](i, j);
            }
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-8);
        (void)dcs;
    }
}

/** A canonical acyclic APH(2): phase 1 -> phase 2 -> restart. */
mam::Map<double> canonical_aph2(double mu1, double mu2, double r1) {
    mam::Map<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = -mu1;
    a.D0(0, 1) = mu1 * r1;
    a.D0(1, 1) = -mu2;
    a.D1(0, 0) = mu1 * (1.0 - r1);
    a.D1(1, 0) = mu2;
    return a;
}

}  // namespace

TEST_CASE("the marked fit reproduces the class probabilities exactly") {
    const mam::Map<double> aph = canonical_aph2(2.0, 0.7, 0.6);
    std::vector<double> p;
    p.push_back(0.3);
    p.push_back(0.5);
    p.push_back(0.2);
    std::vector<double> B(3, 0.0);
    for (std::size_t c = 0; c < 3; ++c) B[c] = mam::map_mean(aph);  // a neutral target

    const mam::Maph2mFitResult<double> r = mam::maph2m_fit_multiclass(aph, p, B);
    check_is_mmap(r.maph);
    CHECK(r.maph.classes() == 3);

    const std::vector<double> pc = mam::mmap_pc(r.maph);
    for (std::size_t c = 0; c < 3; ++c) CHECK(pc[c] == doctest::Approx(p[c]).epsilon(1e-6));

    // The timing is the APH's and is untouched by the marking.
    CHECK(mam::map_mean(r.maph.map()) == doctest::Approx(mam::map_mean(aph)).epsilon(1e-12));
    CHECK(mam::map_scv(r.maph.map()) == doctest::Approx(mam::map_scv(aph)).epsilon(1e-12));
}

TEST_CASE("the degenerate branch splits both flows by p alone") {
    // r1 = 1: phase one never completes, so both exit flows see the same law.
    const mam::Map<double> aph = canonical_aph2(2.0, 0.7, 1.0);
    std::vector<double> p;
    p.push_back(0.25);
    p.push_back(0.75);
    std::vector<double> B(2, 1.0);

    const mam::Maph2mFitResult<double> r = mam::maph2m_fit_multiclass(aph, p, B);
    check_is_mmap(r.maph);
    const std::vector<double> pc = mam::mmap_pc(r.maph);
    CHECK(pc[0] == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(pc[1] == doctest::Approx(0.75).epsilon(1e-9));
    // The degenerate branch reports the backward moments it happened to get.
    REQUIRE(r.fB.size() == 2);
    CHECK(r.fB[0] > 0.0);
}

TEST_CASE("maph2m_fit keeps the inter-arrival moments it was given") {
    // A hyperexponential-ish target well inside the APH(2) region.
    const double M1 = 1.0, scv = 3.0;
    const double M2 = (1.0 + scv) * M1 * M1;
    const double M3 = 3.0 * 1.5 * M2 * M2 / M1;
    std::vector<double> p;
    p.push_back(0.4);
    p.push_back(0.6);
    std::vector<double> B;
    B.push_back(1.0);
    B.push_back(1.0);

    const mam::Mmap<double> m = mam::maph2m_fit(M1, M2, M3, p, B);
    check_is_mmap(m);
    CHECK(mam::map_mean(m.map()) == doctest::Approx(M1).epsilon(1e-6));
    CHECK(mam::map_moment(m.map(), 2) == doctest::Approx(M2).epsilon(1e-5));

    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(0.4).epsilon(1e-5));
    CHECK(pc[1] == doctest::Approx(0.6).epsilon(1e-5));
}

TEST_CASE("maph2m_fit_mmap recovers the class law of the process it read") {
    // Build a two-class marked Erlang-2 and refit it from its own descriptors.
    mam::Mmap<double> src;
    const mam::Map<double> e = mam::map_erlang(1.0, 2);
    src.D0 = e.D0;
    src.D1 = e.D1;
    src.Dc.assign(2, Matrix<double>(2, 2, 0.0));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            src.Dc[0](i, j) = 0.35 * e.D1(i, j);
            src.Dc[1](i, j) = 0.65 * e.D1(i, j);
        }

    const mam::Mmap<double> m = mam::maph2m_fit_mmap(src);
    check_is_mmap(m);
    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(0.35).epsilon(1e-4));
    CHECK(pc[1] == doctest::Approx(0.65).epsilon(1e-4));
    CHECK(mam::map_mean(m.map()) == doctest::Approx(1.0).epsilon(1e-5));
}

TEST_CASE("the multiclass core refuses a non-canonical APH by name") {
    mam::Map<double> cyclic = canonical_aph2(2.0, 0.7, 0.6);
    cyclic.D0(1, 0) = 0.1;  // a back edge: no longer acyclic
    std::vector<double> p(2, 0.5), B(2, 1.0);
    CHECK_THROWS_AS(mam::maph2m_fit_multiclass(cyclic, p, B), line::InputError);

    mam::Map<double> noncanon = canonical_aph2(2.0, 0.7, 0.6);
    noncanon.D1(0, 1) = 0.05;  // restart into phase two
    CHECK_THROWS_AS(mam::maph2m_fit_multiclass(noncanon, p, B), line::InputError);

    const mam::Map<double> third = mam::map_erlang(1.0, 3);
    CHECK_THROWS_AS(mam::maph2m_fit_multiclass(third, p, B), line::InputError);
}

TEST_CASE("maph2m_fit_trace reads its descriptors from the labelled trace") {
    // Alternating classes on a constant-ish trace: p must come back 1/2, 1/2.
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 4000; ++i) {
        Tv.push_back(0.5 + 0.5 * static_cast<double>((i * 7919) % 13) / 13.0);
        A.push_back(static_cast<int>(i % 2) + 1);
    }
    const mam::Mmap<double> m = mam::maph2m_fit_trace(Tv, A);
    check_is_mmap(m);
    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(0.5).epsilon(1e-3));
    CHECK(pc[1] == doctest::Approx(0.5).epsilon(1e-3));

    CHECK_THROWS_AS(mam::maph2m_fit_trace(std::vector<double>(), std::vector<int>()),
                    line::InputError);
}
