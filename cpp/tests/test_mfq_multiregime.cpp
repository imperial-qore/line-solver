/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The multi-regime feedback fluid queue (line/api/mam/mfq_multiregime.h) and
 * the real Schur infrastructure it needed (line/util/eig.h: schur_decomposition
 * and schur_reorder).
 *
 * Oracles, in the order the task prescribes.
 *  (a) A closed form the model collapses to. One regime with the threshold
 *      pushed far out is the ordinary homogeneous fluid queue on [0, inf), and
 *      on Q = [-2 2; 1 -1] with drifts (+1, -1) that queue is solvable by hand:
 *      the Riccati equation reduces to Psi^2 - 3 Psi + 2 = 0 with minimal root
 *      1, giving K = -1, ini = 1/3 and clo = [1 1], so the density is exactly
 *      pi(x) = (1/3) e^{-x} [1 1] with an atom of 1/3 at level zero. The port
 *      reproduces every one of those to 1e-12 through a completely different
 *      route (censoring, ordered Schur form, two Sylvester solves and a 6 x 6
 *      boundary system), which is the strongest check available here.
 *  (b) Invariants. Total mass one; cdf non-decreasing; cdfm - cdf equal to the
 *      atom at the point and zero elsewhere; pdfd the derivative of pdf; and
 *      the increment of cdf across an interval equal to the integral of pdf
 *      over it. For the Schur routines: A = Z T Z^T, Z orthogonal, T
 *      quasi-triangular, and after reordering the diagonal blocks in
 *      descending key order with 2 x 2 blocks left intact.
 *  (c) MATLAB, digit for digit, on a two-regime instance with different
 *      generators and different drift vectors per regime.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mfq_multiregime.h"
#include "line/api/mam/mfq_solve.h"
#include "line/util/eig.h"
#include "line/util/linalg.h"

using line::Matrix;
using line::RealSchur;
using line::mam::mfq_multiregime;

namespace {

Matrix<double> mat(const std::vector<std::vector<double>>& a) {
    Matrix<double> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) m(i, j) = a[i][j];
    return m;
}

/** max |Z T Z^T - A|. */
double schur_residual(const Matrix<double>& A, const RealSchur& s) {
    const std::size_t n = A.rows();
    Matrix<double> Zt(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Zt(i, j) = s.Z(j, i);
    const Matrix<double> R = line::matmul(line::matmul(s.Z, s.T), Zt);
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) m = std::max(m, std::fabs(R(i, j) - A(i, j)));
    return m;
}

/** max |Z^T Z - I|. */
double orthogonality_defect(const RealSchur& s) {
    const std::size_t n = s.Z.rows();
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            double d = 0.0;
            for (std::size_t k = 0; k < n; ++k) d += s.Z(k, i) * s.Z(k, j);
            m = std::max(m, std::fabs(d - (i == j ? 1.0 : 0.0)));
        }
    return m;
}

/** The sign-class key the reference's ordschur call uses. */
std::vector<double> sign_key(const Matrix<double>& T) {
    const std::size_t n = T.rows();
    std::vector<double> k(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double d = T(i, i);
        k[i] = (std::fabs(d) < 1e-10 ? 5.0 : 0.0) + (d < 0.0 ? 2.0 : 0.0) + (d > 0.0 ? 1.0 : 0.0);
    }
    for (std::size_t i = 0; i + 1 < n; ++i)
        if (T(i + 1, i) != 0.0) k[i + 1] = k[i];
    return k;
}

/** The two-regime instance checked against MATLAB. */
struct TwoRegime {
    std::vector<Matrix<double>> Q;
    std::vector<std::vector<double>> R;
    std::vector<double> T;
};

TwoRegime two_regime() {
    TwoRegime t;
    t.Q = {mat({{-2.0, 2.0}, {1.0, -1.0}}), mat({{-3.0, 3.0}, {2.0, -2.0}})};
    t.R = {{1.0, -1.0}, {0.5, -2.0}};
    t.T = {1.0, 3.0};
    return t;
}

}  // namespace

// ---------------------------------------------------------------------------
// the Schur infrastructure added to util/eig.h
// ---------------------------------------------------------------------------

TEST_CASE("schur_decomposition factors A = Z T Z^T with Z orthogonal") {
    // Deterministic pseudo-random matrices, several orders, including sizes
    // that reliably produce complex conjugate pairs and hence 2 x 2 blocks.
    unsigned seed = 12345u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return static_cast<double>((seed >> 16) % 2001u) / 250.0 - 4.0;
    };
    for (std::size_t n : {2u, 3u, 5u, 8u}) {
        Matrix<double> A(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = next();
        const RealSchur s = line::schur_decomposition(A);
        INFO("n = ", n);
        CHECK(schur_residual(A, s) < 1e-12);
        CHECK(orthogonality_defect(s) < 1e-13);
        // Quasi-triangular: nothing below the first subdiagonal, and no two
        // consecutive non-zero subdiagonal entries.
        for (std::size_t i = 2; i < n; ++i)
            for (std::size_t j = 0; j + 2 <= i; ++j) CHECK(std::fabs(s.T(i, j)) < 1e-13);
        for (std::size_t i = 0; i + 2 < n; ++i)
            {
                const bool adjacent_blocks =
                    (s.T(i + 1, i) != 0.0) && (s.T(i + 2, i + 1) != 0.0);
                CHECK_FALSE(adjacent_blocks);
            }
    }
}

TEST_CASE("schur_reorder sorts the diagonal blocks by descending key") {
    unsigned seed = 999u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return static_cast<double>((seed >> 16) % 2001u) / 250.0 - 4.0;
    };
    for (std::size_t n : {3u, 5u, 7u}) {
        Matrix<double> A(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = next();
        const RealSchur s = line::schur_decomposition(A);
        const RealSchur r = line::schur_reorder(s, sign_key(s.T));
        INFO("n = ", n);
        // The factorization survives the reordering.
        CHECK(schur_residual(A, r) < 1e-11);
        CHECK(orthogonality_defect(r) < 1e-12);
        // The keys now come out non-increasing.
        const std::vector<double> k = sign_key(r.T);
        for (std::size_t i = 1; i < n; ++i) CHECK(k[i] <= k[i - 1]);
        // The spectrum is unchanged as a multiset of the diagonal.
        std::vector<double> before, after;
        for (std::size_t i = 0; i < n; ++i) {
            before.push_back(s.T(i, i));
            after.push_back(r.T(i, i));
        }
        std::sort(before.begin(), before.end());
        std::sort(after.begin(), after.end());
        for (std::size_t i = 0; i < n; ++i)
            CHECK(after[i] == doctest::Approx(before[i]).epsilon(1e-10));
    }
}

TEST_CASE("schur_reorder is the identity when the key is already sorted") {
    const Matrix<double> A = mat({{-2.0, 2.0}, {1.0, -1.0}});
    const RealSchur s = line::schur_decomposition(A);
    const std::vector<double> k = sign_key(s.T);
    std::vector<double> sorted = k;
    std::sort(sorted.rbegin(), sorted.rend());
    if (k == sorted) {
        const RealSchur r = line::schur_reorder(s, k);
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) {
                CHECK(r.T(i, j) == doctest::Approx(s.T(i, j)).epsilon(1e-14));
                CHECK(r.Z(i, j) == doctest::Approx(s.Z(i, j)).epsilon(1e-14));
            }
    }
}

// ---------------------------------------------------------------------------
// (a) the closed form the model collapses to
// ---------------------------------------------------------------------------

TEST_CASE("one regime with a far threshold is the homogeneous fluid queue") {
    // Q = [-2 2; 1 -1], drifts (+1, -1). By hand: Psi = 1, K = -1, ini = 1/3,
    // clo = [1 1], so pi(x) = (1/3) e^{-x} [1 1] and the atom at zero is
    // (0, 1/3). The port reaches it through censoring, an ordered Schur form,
    // two Sylvester solves and a 6 x 6 boundary system.
    const std::vector<Matrix<double>> Q = {mat({{-2.0, 2.0}, {1.0, -1.0}})};
    const std::vector<std::vector<double>> R = {{1.0, -1.0}};
    const std::vector<double> T = {40.0};
    const std::vector<double> pts = {0.5, 1.0, 2.0, 3.0};
    const auto r = mfq_multiregime(Q, R, {}, {}, T, pts, std::vector<double>{0.0});
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double exact = std::exp(-pts[i]) / 3.0;
        INFO("x = ", pts[i]);
        CHECK(r.pdf[i][0] == doctest::Approx(exact).epsilon(1e-12));
        CHECK(r.pdf[i][1] == doctest::Approx(exact).epsilon(1e-12));
        // The density decays as e^{-x}, so its derivative is its negative.
        CHECK(r.pdfd[i][0] == doctest::Approx(-exact).epsilon(1e-11));
    }
    // The atom at level zero sits entirely in the down-drift state.
    CHECK(r.cdfm[0][0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.cdfm[0][1] == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    // and the up-drift state carries no mass at zero, so Cdf and Cdfm differ
    // there by exactly the atom.
    CHECK(r.cdf[0][1] == doctest::Approx(0.0).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// (b) invariants
// ---------------------------------------------------------------------------

TEST_CASE("the multi-regime law is a probability law") {
    const TwoRegime t = two_regime();
    // Cdfm at the top threshold is the whole mass.
    const auto top = mfq_multiregime(t.Q, t.R, {}, {}, t.T, {}, std::vector<double>{3.0});
    double total = 0.0;
    for (double v : top.cdfm[0]) total += v;
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));

    // Non-decreasing in each state.
    const std::vector<double> pts = {0.0, 0.4, 1.0, 1.6, 2.4, 3.0};
    const auto r = mfq_multiregime(t.Q, t.R, {}, {}, t.T, {}, pts);
    for (std::size_t i = 1; i < pts.size(); ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(r.cdf[i][j] >= r.cdf[i - 1][j] - 1e-14);
    // Cdfm - Cdf is the atom at the point: non-negative everywhere, and zero
    // away from a threshold.
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            const double atom = r.cdfm[i][j] - r.cdf[i][j];
            CHECK(atom >= -1e-14);
            const bool at_threshold = (pts[i] == 0.0 || pts[i] == 1.0 || pts[i] == 3.0);
            if (!at_threshold) CHECK(atom == doctest::Approx(0.0).epsilon(1e-12));
        }
}

TEST_CASE("pdfd is the derivative of pdf, and pdf the derivative of cdf") {
    const TwoRegime t = two_regime();
    const double h = 1e-6;
    for (double p : {0.4, 0.8, 1.7, 2.6}) {
        INFO("p = ", p);
        const auto c = mfq_multiregime(t.Q, t.R, {}, {}, t.T, {p - h, p, p + h},
                                       std::vector<double>{p - h, p + h});
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK((c.pdf[2][j] - c.pdf[0][j]) / (2 * h) ==
                  doctest::Approx(c.pdfd[1][j]).epsilon(1e-5));
            CHECK((c.cdf[1][j] - c.cdf[0][j]) / (2 * h) ==
                  doctest::Approx(c.pdf[1][j]).epsilon(1e-5));
        }
    }
}

TEST_CASE("the cdf increment over an interval is the integral of the pdf") {
    // Inside one regime, so no atom and no kink is crossed.
    const TwoRegime t = two_regime();
    const double a = 1.2, b = 2.8;
    const std::size_t M = 4001;
    std::vector<double> grid(M);
    for (std::size_t i = 0; i < M; ++i)
        grid[i] = a + (b - a) * static_cast<double>(i) / (M - 1);
    const auto r = mfq_multiregime(t.Q, t.R, {}, {}, t.T, grid, std::vector<double>{a, b});
    for (std::size_t j = 0; j < 2; ++j) {
        double quad = 0.0;
        for (std::size_t i = 1; i < M; ++i)
            quad += 0.5 * (r.pdf[i - 1][j] + r.pdf[i][j]) * (grid[i] - grid[i - 1]);
        INFO("state ", j);
        CHECK(r.cdf[1][j] - r.cdf[0][j] == doctest::Approx(quad).epsilon(1e-8));
    }
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("mfq_multiregime agrees with MATLAB on two regimes") {
    const TwoRegime t = two_regime();
    const std::vector<double> pdfpts = {0.2, 0.9, 1.0, 1.5, 2.5, 3.0};
    const std::vector<double> cdfpts = {0.0, 0.5, 1.0, 2.0, 3.0};
    const auto r = mfq_multiregime(t.Q, t.R, {}, {}, t.T, pdfpts, cdfpts);

    const std::vector<std::vector<double>> refPdf = {
        {0.334424512764062, 0.334424512764062},
        {0.166070298266236, 0.166070298266236},
        {0.300533239791365, 0.0751333099478413},
        {0.0246692705747101, 0.00616731764367756},
        {0.000166220237638526, 4.15550594096709e-05},
        {1.36441879778515e-05, 3.4110469945022e-06}};
    const std::vector<std::vector<double>> refPdfd = {
        {-0.334424512764062, -0.334424512764062},
        {-0.166070298266236, -0.166070298266236},
        {-1.50266619895682, -0.375666549739206},
        {-0.12334635287355, -0.0308365882183876},
        {-0.000831101188192475, -0.000207775297048119},
        {-6.82209398890999e-05, -1.7055234972275e-05}};
    const std::vector<std::vector<double>> refCdf = {
        {0.0, 0.0},
        {0.160719249788173, 0.56918627207457},
        {0.258200402390714, 0.66666742467711},
        {0.317902054940751, 0.681592837814619},
        {0.318304321511391, 0.681693404457279}};
    const std::vector<std::vector<double>> refCdfm = {
        {0.0, 0.408467022286396},
        {0.160719249788173, 0.56918627207457},
        {0.258200402390714, 0.66666742467711},
        {0.317902054940751, 0.681592837814619},
        {0.318306595542721, 0.681693404457279}};

    for (std::size_t i = 0; i < pdfpts.size(); ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            INFO("pdf point ", i, " state ", j);
            CHECK(r.pdf[i][j] == doctest::Approx(refPdf[i][j]).epsilon(1e-11));
            CHECK(r.pdfd[i][j] == doctest::Approx(refPdfd[i][j]).epsilon(1e-11));
        }
    for (std::size_t i = 0; i < cdfpts.size(); ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            INFO("cdf point ", i, " state ", j);
            // MATLAB returns cdfm(0,0) as -6.26e-17, a rounding of the exact
            // zero the port produces; compared as an absolute quantity.
            CHECK(std::fabs(r.cdf[i][j] - refCdf[i][j]) < 1e-11);
            CHECK(std::fabs(r.cdfm[i][j] - refCdfm[i][j]) < 1e-11);
        }
}

TEST_CASE("mfq_multiregime rejects a malformed instance") {
    const TwoRegime t = two_regime();
    // One threshold per regime.
    CHECK_THROWS_AS(mfq_multiregime(t.Q, t.R, {}, {}, std::vector<double>{1.0}, {},
                                    std::vector<double>{1.0}),
                    line::InputError);
    // Negative evaluation point.
    CHECK_THROWS_AS(
        mfq_multiregime(t.Q, t.R, {}, {}, t.T, std::vector<double>{-1.0}, {}),
        line::InputError);
    // A regime in which every state has zero drift has no dynamics at all.
    std::vector<std::vector<double>> Rbad = t.R;
    Rbad[0] = {0.0, 0.0};
    CHECK_THROWS_AS(mfq_multiregime(t.Q, Rbad, {}, {}, t.T, {}, std::vector<double>{1.0}),
                    line::InputError);
}
