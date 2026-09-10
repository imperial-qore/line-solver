/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Level-dependent (multi-regime) fluid queues, first and second order
 * (line/api/mam/mfq_ld_solve.h). This closes the mfq family: the blocks it
 * returns are the ones mfq_ld_mean and mfq_ld_distr already consume.
 *
 * Oracles, in the order the task prescribes.
 *  (a) Two closed forms. One regime with a far threshold is the homogeneous
 *      fluid queue, hand-solvable on Q = [-2 2; 1 -1] with drifts (+1,-1):
 *      Psi = 1, K = -1, an atom of 1/3 at zero and density (1/3) e^{-x}. And
 *      the QBD identity itself, which is the thing that went wrong first and is
 *      therefore asserted directly rather than inferred: the R this port feeds
 *      to the exponent must satisfy Fm + R Lm + R^2 Bm = 0 on the actual triple
 *      the algorithm builds.
 *  (b) Invariants: total mass one; masses non-negative; and a CROSS-CHECK
 *      against mfq_multiregime, which solves the same first-order model by a
 *      completely different route (ordered Schur form and Sylvester solves
 *      rather than a QBD per direction). Two independently ported functions
 *      meeting on one law is the strongest evidence available here that neither
 *      is self-consistently wrong.
 *  (c) MATLAB, on a first-order two-regime instance, a genuinely SECOND-order
 *      one with non-zero variance in both regimes, and a single-regime one.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mfq_ld_distr.h"
#include "line/api/mam/mfq_ld_mean.h"
#include "line/api/mam/mfq_ld_solve.h"
#include "line/api/mam/mfq_multiregime.h"

using line::Matrix;
using line::mam::FluidBoundary;
using line::mam::FluidDistrKind;
using line::mam::LevelDependentFluidBlocks;
using line::mam::mfq_ld_distr;
using line::mam::mfq_ld_mean;
using line::mam::mfq_ld_solve;

namespace {

Matrix<double> mat(const std::vector<std::vector<double>>& a) {
    Matrix<double> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) m(i, j) = a[i][j];
    return m;
}

Matrix<double> dg(const std::vector<double>& v) {
    Matrix<double> m(v.size(), v.size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i) m(i, i) = v[i];
    return m;
}

/** Total mass of the law, read off Cdfm at the top threshold. */
double total_mass(const LevelDependentFluidBlocks<double>& b) {
    const auto r = mfq_ld_distr(b, FluidDistrKind::Cdfm, std::vector<double>{b.Thr.back()});
    double s = 0.0;
    for (double v : r[0]) s += v;
    return s;
}

const std::vector<std::vector<double>> kQ1 = {{-2.0, 2.0}, {1.0, -1.0}};
const std::vector<std::vector<double>> kQ2 = {{-3.0, 3.0}, {2.0, -2.0}};

std::vector<Matrix<double>> two_Q() { return {mat(kQ1), mat(kQ2)}; }
std::vector<Matrix<double>> two_R() { return {dg({1.0, -1.0}), dg({0.5, -2.0})}; }
std::vector<Matrix<double>> zero_S() { return {dg({0.0, 0.0}), dg({0.0, 0.0})}; }

}  // namespace

// ---------------------------------------------------------------------------
// (a) the QBD identity, pinned rather than inferred
// ---------------------------------------------------------------------------

TEST_CASE("the QBD R satisfies Fm + R Lm + R^2 Bm = 0 on the triple the algorithm builds") {
    // This is the regression for the mistake that cost the first attempt. The
    // reference calls QBD_CR(Bm, Lm, Fm), and QBD_CR's error message ("the
    // matrix A0+A1+A2 has to be (sub)stochastic") invites the reading that it
    // wants the DISCRETE form R = A0 + R A1 + R^2 A2 with A0 = Bm, which would
    // put a -I shift on L. The triples this caller builds sum to a GENERATOR,
    // and that reading produces a completely different, wrong K. The triple
    // below is the exact forward one for regime 1 of the first-order instance:
    // Q = [-2 2; 1 -1], R = diag(1, -1), S = 0, giving c = 2.
    const Matrix<double> Bm = mat({{0.0, 0.0}, {0.0, 1.0}});
    const Matrix<double> Lm = mat({{-1.0, 0.0}, {0.5, -1.5}});
    const Matrix<double> Fm = mat({{0.0, 1.0}, {0.0, 0.0}});
    const Matrix<double> R =
        line::mam::ld_solve_detail::qbd_cr_R(Bm, Lm, Fm, 1e-14);
    // Residual of the equation the port relies on.
    Matrix<double> res = Fm;
    const Matrix<double> RL = line::matmul(R, Lm);
    const Matrix<double> RRB = line::matmul(line::matmul(R, R), Bm);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) res(i, j) += RL(i, j) + RRB(i, j);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(res(i, j)) < 1e-11);
    // MATLAB's QBD_CR returns R = [0.5 1; 0 0] here, from which the reference
    // reads K = (R(1,1) - 1) c = -1 and Psi = R(1,2) = 1.
    CHECK(R(0, 0) == doctest::Approx(0.5).epsilon(1e-11));
    CHECK(R(0, 1) == doctest::Approx(1.0).epsilon(1e-11));
    // The discrete reading would have needed R to satisfy this instead; it does
    // not, by a wide margin, which is what makes the two readings separable.
    Matrix<double> alt = Bm;
    const Matrix<double> RRF = line::matmul(line::matmul(R, R), Fm);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) alt(i, j) += RL(i, j) + RRF(i, j) - R(i, j);
    double worst = 0.0;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) worst = std::max(worst, std::fabs(alt(i, j)));
    CHECK(worst > 1.0);
}

TEST_CASE("one regime with a far threshold is the homogeneous fluid queue") {
    // Hand-solved: Psi^2 - 3 Psi + 2 = 0 with minimal root 1, so K = -1, the
    // atom at zero is (0, 1/3) and the density is (1/3) e^{-x} [1 1].
    const std::vector<Matrix<double>> Q = {mat(kQ1)};
    const std::vector<Matrix<double>> R = {dg({1.0, -1.0})};
    const std::vector<Matrix<double>> S = {dg({0.0, 0.0})};
    const auto b = mfq_ld_solve(Q, R, S, std::vector<double>{40.0});
    REQUIRE(b.KF.size() == 1u);
    CHECK(b.KF[0].rows() == 1u);
    CHECK(b.KF[0](0, 0) == doctest::Approx(-1.0).epsilon(1e-10));
    CHECK(b.masses[0][0] == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(b.masses[0][1] == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    const auto pdf = mfq_ld_distr(b, FluidDistrKind::Pdf, std::vector<double>{0.5, 1.0, 2.0});
    const std::vector<double> xs = {0.5, 1.0, 2.0};
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const double exact = std::exp(-xs[i]) / 3.0;
        INFO("x = ", xs[i]);
        CHECK(pdf[i][0] == doctest::Approx(exact).epsilon(1e-9));
        CHECK(pdf[i][1] == doctest::Approx(exact).epsilon(1e-9));
    }
    CHECK(total_mass(b) == doctest::Approx(1.0).epsilon(1e-10));
}

// ---------------------------------------------------------------------------
// (b) invariants, and the cross-check against mfq_multiregime
// ---------------------------------------------------------------------------

TEST_CASE("the level-dependent law is a probability law") {
    struct Case {
        std::vector<Matrix<double>> S;
        const char* name;
    };
    const std::vector<Case> cases = {{zero_S(), "first order"},
                                     {{dg({0.5, 0.5}), dg({0.3, 0.8})}, "second order"}};
    for (const Case& c : cases) {
        INFO(c.name);
        const auto b = mfq_ld_solve(two_Q(), two_R(), c.S, std::vector<double>{1.0, 3.0});
        CHECK(total_mass(b) == doctest::Approx(1.0).epsilon(1e-11));
        for (const std::vector<double>& m : b.masses)
            for (double v : m) CHECK(v > -1e-10);
        // The Cdfm is non-decreasing across the whole support.
        const std::vector<double> pts = {0.0, 0.5, 1.0, 1.8, 2.5, 3.0};
        const auto cdf = mfq_ld_distr(b, FluidDistrKind::Cdfm, pts);
        for (std::size_t i = 1; i < pts.size(); ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(cdf[i][j] >= cdf[i - 1][j] - 1e-10);
        // and the mean sits inside the support.
        const double mu = mfq_ld_mean(b);
        CHECK(mu > 0.0);
        CHECK(mu < 3.0);
    }
}

TEST_CASE("the first-order solve agrees with mfq_multiregime on the same model") {
    // Two independently ported functions on one physical model. mfq_ld_solve
    // builds a QBD per direction per regime; mfq_multiregime builds an ordered
    // real Schur form and two Sylvester solves. They share only the model.
    const auto b = mfq_ld_solve(two_Q(), two_R(), zero_S(), std::vector<double>{1.0, 3.0});
    const auto ld = mfq_ld_distr(b, FluidDistrKind::Cdfm, std::vector<double>{3.0});

    const std::vector<std::vector<double>> Rvec = {{1.0, -1.0}, {0.5, -2.0}};
    const auto mr = line::mam::mfq_multiregime(two_Q(), Rvec, {}, {},
                                               std::vector<double>{1.0, 3.0}, {},
                                               std::vector<double>{3.0});
    for (std::size_t j = 0; j < 2; ++j) {
        INFO("state ", j);
        CHECK(ld[0][j] == doctest::Approx(mr.cdfm[0][j]).epsilon(1e-10));
    }
    // and both are the same probability law, not merely equal to each other.
    double a = 0.0, c = 0.0;
    for (std::size_t j = 0; j < 2; ++j) {
        a += ld[0][j];
        c += mr.cdfm[0][j];
    }
    CHECK(a == doctest::Approx(1.0).epsilon(1e-11));
    CHECK(c == doctest::Approx(1.0).epsilon(1e-11));
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("mfq_ld_solve agrees with MATLAB") {
    // Tolerances are 1e-10 rather than 1e-13 because the port reaches R by the
    // linearly convergent qbd_R where MATLAB uses cyclic reduction; the blocks
    // land about 2e-14 apart and that propagates. Tightening would need a
    // quadratically convergent QBD solver, not a smaller tolerance.
    SUBCASE("first order, two regimes") {
        const auto b = mfq_ld_solve(two_Q(), two_R(), zero_S(), std::vector<double>{1.0, 3.0});
        CHECK(b.masses[0][1] == doctest::Approx(0.40846702228639625).epsilon(1e-10));
        CHECK(std::fabs(b.masses[1][0]) < 1e-9);
        CHECK(std::fabs(b.masses[1][1]) < 1e-9);
        CHECK(b.masses[2][0] == doctest::Approx(2.274031329625606e-06).epsilon(1e-8));
        CHECK(b.KF[0](0, 0) == doctest::Approx(-1.0000000000000002).epsilon(1e-10));
        CHECK(b.KF[1](0, 0) == doctest::Approx(-5.0).epsilon(1e-10));
        CHECK(b.cloF[0](0, 1) == doctest::Approx(1.0).epsilon(1e-10));
        CHECK(b.cloF[1](0, 1) == doctest::Approx(0.25).epsilon(1e-10));
        CHECK(b.cloB[0](0, 0) == doctest::Approx(0.5).epsilon(1e-10));
        CHECK(b.cloB[1](0, 0) == doctest::Approx(2.0 / 3.0).epsilon(1e-10));
        CHECK(mfq_ld_mean(b) == doctest::Approx(0.30602344367107848).epsilon(1e-11));
        const auto r = mfq_ld_distr(b, FluidDistrKind::Cdfm, std::vector<double>{3.0});
        CHECK(r[0][0] == doctest::Approx(0.31830659554272089).epsilon(1e-10));
        CHECK(r[0][1] == doctest::Approx(0.68169340445727922).epsilon(1e-10));
    }
    SUBCASE("second order, non-zero variance in both regimes") {
        const std::vector<Matrix<double>> S = {dg({0.5, 0.5}), dg({0.3, 0.8})};
        const auto b = mfq_ld_solve(two_Q(), two_R(), S, std::vector<double>{1.0, 3.0});
        // Every regime is now fully second order, so both exponents are 2 x 2.
        REQUIRE(b.KF[0].rows() == 2u);
        REQUIRE(b.KB[0].rows() == 2u);
        CHECK(b.iniF[0][0] == doctest::Approx(0.14961505685939519).epsilon(1e-9));
        CHECK(b.iniF[0][1] == doctest::Approx(1.8365746676209236).epsilon(1e-10));
        CHECK(b.iniB[0][0] == doctest::Approx(0.014637918675819673).epsilon(1e-9));
        CHECK(b.iniF[1][0] == doctest::Approx(0.23853530737254458).epsilon(1e-10));
        CHECK(b.iniF[1][1] == doctest::Approx(0.18698742343641714).epsilon(1e-10));
        CHECK(b.KF[0](0, 0) == doctest::Approx(-1.1383690972587435).epsilon(1e-10));
        CHECK(b.KF[0](0, 1) == doctest::Approx(5.1383690972587432).epsilon(1e-10));
        CHECK(b.KB[0](1, 0) == doctest::Approx(2.5691845486293676).epsilon(1e-10));
        CHECK(b.KF[1](1, 1) == doctest::Approx(-5.4316271137897445).epsilon(1e-10));
        // A second-order queue puts no mass at an interior threshold.
        for (const std::vector<double>& m : b.masses)
            for (double v : m) CHECK(std::fabs(v) < 1e-9);
        CHECK(mfq_ld_mean(b) == doctest::Approx(0.59089112156266232).epsilon(1e-10));
        const auto r = mfq_ld_distr(b, FluidDistrKind::Cdfm, std::vector<double>{3.0});
        CHECK(r[0][0] == doctest::Approx(0.32175270467503359).epsilon(1e-10));
        CHECK(r[0][1] == doctest::Approx(0.67824729532496641).epsilon(1e-10));
    }
    SUBCASE("single regime") {
        const std::vector<Matrix<double>> Q = {mat(kQ1)};
        const std::vector<Matrix<double>> R = {dg({1.0, -1.0})};
        const std::vector<Matrix<double>> S = {dg({0.0, 0.0})};
        const auto b = mfq_ld_solve(Q, R, S, std::vector<double>{2.0});
        CHECK(b.masses[0][1] == doctest::Approx(0.3575262944985847).epsilon(1e-10));
        CHECK(b.masses[1][0] == doctest::Approx(0.024192961165251266).epsilon(1e-9));
        CHECK(b.iniF[0][0] == doctest::Approx(0.3575262944985847).epsilon(1e-10));
        CHECK(b.KF[0](0, 0) == doctest::Approx(-1.0000000000000002).epsilon(1e-10));
        CHECK(mfq_ld_mean(b) == doctest::Approx(0.47312297734465641).epsilon(1e-10));
        const auto r = mfq_ld_distr(b, FluidDistrKind::Cdfm, std::vector<double>{2.0});
        CHECK(r[0][0] == doctest::Approx(1.0 / 3.0).epsilon(1e-10));
        CHECK(r[0][1] == doctest::Approx(2.0 / 3.0).epsilon(1e-10));
    }
}

TEST_CASE("mfq_ld_solve takes boundary flags per background state, not per regime") {
    // The reference's default is zeros(K,N) and crashes for K >= 2 on its own
    // arguments; the documented contract, which this port implements, is one
    // flag per BACKGROUND STATE. A wrongly sized vector is rejected by name.
    const std::vector<double> T = {1.0, 3.0};
    CHECK_NOTHROW(mfq_ld_solve(two_Q(), two_R(), zero_S(), T));
    const std::vector<FluidBoundary> ok(2, FluidBoundary::Reflective);
    CHECK_NOTHROW(mfq_ld_solve(two_Q(), two_R(), zero_S(), T, ok, ok,
                               std::vector<Matrix<double>>(), 1e-14));
    const std::vector<FluidBoundary> wrong(4, FluidBoundary::Reflective);
    CHECK_THROWS_AS(mfq_ld_solve(two_Q(), two_R(), zero_S(), T, wrong, wrong,
                                 std::vector<Matrix<double>>(), 1e-14),
                    line::InputError);
    // and the defaults really are reflective, i.e. the same as passing them.
    const auto a = mfq_ld_solve(two_Q(), two_R(), zero_S(), T);
    const auto c = mfq_ld_solve(two_Q(), two_R(), zero_S(), T, ok, ok,
                                std::vector<Matrix<double>>(), 1e-14);
    CHECK(mfq_ld_mean(a) == doctest::Approx(mfq_ld_mean(c)).epsilon(1e-13));
}

TEST_CASE("mfq_ld_solve rejects a malformed instance") {
    CHECK_THROWS_AS(mfq_ld_solve(two_Q(), two_R(), zero_S(), std::vector<double>{1.0}),
                    line::InputError);
    CHECK_THROWS_AS(
        mfq_ld_solve(two_Q(), two_R(), zero_S(), std::vector<double>()), line::InputError);
    // A regime with neither drift nor variance anywhere has no dynamics.
    const std::vector<Matrix<double>> Rbad = {dg({0.0, 0.0}), dg({0.5, -2.0})};
    CHECK_THROWS_AS(
        mfq_ld_solve(two_Q(), Rbad, zero_S(), std::vector<double>{1.0, 3.0}),
        line::InputError);
}
