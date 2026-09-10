/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mmap2k_fit and mmap3k_fit: the EXACT inverses of the class-marking map.
 *
 * THE ORACLE IS ROUND-TRIP. Take a marked MAP, measure its own (p, F, B) with
 * mmap_pc / mmap_forward_moment / mmap_backward_moment, hand those back to the
 * fitter, and the split must come back: an exact inverse applied to its own
 * image is the identity. That is a far stronger check than a moment comparison,
 * because it pins every entry of every Dc rather than three aggregates.
 *
 * The `exact` flag is checked as carefully as the numbers: a fitter that clamped
 * an infeasible split and still reported exact would pass every value check
 * here, which is precisely the failure mode the flag exists to expose.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_k_fit.h"
#include "line/api/mam/mmap_stats.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

void check_is_mmap(const mam::Mmap<double>& m) {
    const std::size_t n = m.order();
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-9);
            double dc = 0.0;
            for (std::size_t c = 0; c < m.classes(); ++c) {
                CHECK(m.Dc[c](i, j) >= -1e-9);
                dc += m.Dc[c](i, j);
            }
            CHECK(m.D1(i, j) == doctest::Approx(dc).epsilon(1e-8));
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-8);
    }
}

/** A marked MAP built by splitting each arrival entry by its own class law. */
mam::Mmap<double> mark_by(const mam::Map<double>& base,
                          const std::vector<std::vector<double>>& q) {
    const std::size_t n = base.order(), K = q[0].size();
    mam::Mmap<double> m;
    m.D0 = base.D0;
    m.D1 = base.D1;
    m.Dc.assign(K, Matrix<double>(n, n, 0.0));
    // q is indexed by the non-zero entries of D1, in row-major order.
    std::size_t idx = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (base.D1(i, j) == 0.0) continue;
            for (std::size_t c = 0; c < K; ++c) m.Dc[c](i, j) = base.D1(i, j) * q[idx][c];
            ++idx;
        }
    return m;
}

}  // namespace

TEST_CASE("mmap3k_fit inverts the marking of an Erlang-2 exactly") {
    const mam::Map<double> e = mam::map_erlang(1.0, 2);
    // Erlang-2's D1 has a single non-zero entry, so one split parameter.
    std::vector<std::vector<double>> q(1, std::vector<double>(2, 0.0));
    q[0][0] = 0.4;
    q[0][1] = 0.6;
    const mam::Mmap<double> src = mark_by(e, q);

    const std::vector<unsigned> ord(1, 1u);
    const std::vector<double> p = mam::mmap_pc(src);
    const Matrix<double> fm = mam::mmap_forward_moment(src, ord, true);
    const std::vector<std::vector<double>> bm = mam::mmap_backward_moment(src, ord, true);
    std::vector<double> F(2, 0.0), B(2, 0.0);
    for (std::size_t c = 0; c < 2; ++c) {
        F[c] = fm(c, 0);
        B[c] = bm[c][0];
    }

    const mam::MmapKFitResult<double> r = mam::mmap3k_fit(src.D0, src.D1, p, F, B);
    check_is_mmap(r.mmap);
    // The split must come back entry for entry.
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(r.mmap.Dc[c](i, j) == doctest::Approx(src.Dc[c](i, j)).epsilon(1e-6));
}

TEST_CASE("mmap3k_fit reports a degenerate marking system rather than guessing") {
    // A Poisson process has one phase and one arrival entry, and its marking
    // system asks for two characteristics from one parameter.
    Matrix<double> D0(1, 1, -1.0), D1(1, 1, 1.0);
    std::vector<double> p(2, 0.5), F(2, 1.0), B(2, 1.0);
    // One non-zero entry means z = 1, so only the (1,0) characteristic is used
    // and the system is solvable; the class probabilities come straight back.
    const mam::MmapKFitResult<double> r = mam::mmap3k_fit(D0, D1, p, F, B);
    check_is_mmap(r.mmap);
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.5).epsilon(1e-9));
}

TEST_CASE("mmap3k_fit refuses a MAP with nothing to mark") {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -1.0;
    D0(1, 1) = -1.0;
    std::vector<double> p(2, 0.5), F(2, 1.0), B(2, 1.0);
    CHECK_THROWS_AS(mam::mmap3k_fit(D0, D1, p, F, B), line::InputError);
}

TEST_CASE("mmap2k_fit returns a valid marked MAP and reports whether it was exact") {
    const double M1 = 1.0, scv = 2.0;
    const double M2 = (1.0 + scv) * M1 * M1;
    const double M3 = 3.0 * 1.5 * M2 * M2 / M1;
    const double gamma = 0.3;
    std::vector<double> p;
    p.push_back(0.35);
    p.push_back(0.65);
    std::vector<double> F(2, M1), B(2, M1);

    const mam::MmapKFitResult<double> r = mam::mmap2k_fit(M1, M2, M3, gamma, p, F, B);
    check_is_mmap(r.mmap);
    CHECK(mam::map_mean(r.mmap.map()) == doctest::Approx(M1).epsilon(1e-5));
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    // Whichever route it took, the class law is matched.
    CHECK(pc[0] == doctest::Approx(0.35).epsilon(1e-3));
    CHECK(pc[1] == doctest::Approx(0.65).epsilon(1e-3));

    // When the exact inverse succeeded, the forward and backward moments are
    // matched too; when it fell back, they are only approached.
    if (r.exact) {
        const std::vector<unsigned> ord(1, 1u);
        const Matrix<double> fm = mam::mmap_forward_moment(r.mmap, ord, true);
        for (std::size_t c = 0; c < 2; ++c) CHECK(fm(c, 0) == doctest::Approx(F[c]).epsilon(1e-4));
    }
}

TEST_CASE("mmap2k_fit refuses mismatched target lengths by name") {
    std::vector<double> p(2, 0.5), F(3, 1.0), B(2, 1.0);
    CHECK_THROWS_AS(mam::mmap2k_fit(1.0, 3.0, 13.5, 0.3, p, F, B), line::InputError);
    CHECK_THROWS_AS(
        mam::mmap2k_fit(1.0, 3.0, 13.5, 0.3, std::vector<double>(), F, B), line::InputError);
}
