/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Theorem 1 of Casale, Harrison and Ong (Perform. Eval. 2021), carried over the
 * divided-difference form of Casale (SIGMETRICS 2017), Cor. 3.2.
 *
 * The closed form is EXACT, so pfqn_gld -- the load-dependent convolution, a
 * structurally different algorithm -- is the oracle throughout, and
 * pfqn_explicit is the oracle on the fixed-rate degeneration mu = 1, which the
 * load-dependent route must reproduce.
 *
 * The one case that is not an agreement check is the near-tie: two induced
 * demands that are equal in exact arithmetic can land two ulps apart in doubles,
 * so the eps-relative redundancy test misses the tie and Eq. (15) divides by the
 * ulp. That is not a defect to paper over -- it is what `tol` controls -- but the
 * loss report has to say so, and a caller holding a budget has to be refused.
 */
#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_explicit.h"
#include "line/api/pfqn/pfqn_explicit_ld.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::pfqn::ExplicitResult;
using line::pfqn::pfqn_explicit;
using line::pfqn::pfqn_explicit_ld;
using line::pfqn::NcldMethod;
using line::pfqn::pfqn_gld;
using line::pfqn::pfqn_ncld;

namespace {

Matrix<double> mat(const std::vector<std::vector<double> >& v) {
    Matrix<double> m(v.size(), v[0].size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v[0].size(); ++j) m(i, j) = v[i][j];
    return m;
}

/** mu(n) = min(n,c) at each of M centers. */
Matrix<double> ms_rates(std::size_t M, std::size_t Nt, int c) {
    Matrix<double> mu(M, Nt, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t n = 1; n <= Nt; ++n)
            mu(i, n - 1) = static_cast<double>(std::min<int>(static_cast<int>(n), c));
    return mu;
}

Matrix<double> ones(std::size_t M, std::size_t Nt) { return Matrix<double>(M, Nt, 1.0); }

void agrees_with_gld(const Matrix<double>& L, const std::vector<int>& N,
                     const Matrix<double>& mu) {
    const double lref = pfqn_gld(L, N, mu).lG;
    const ExplicitResult<double> got = pfqn_explicit_ld(L, N, mu);
    REQUIRE(got.valid);
    REQUIRE(std::isfinite(got.lG));
    CHECK(got.lG == doctest::Approx(lref).epsilon(1e-9));
}

}  // namespace

TEST_CASE("pfqn_explicit_ld: multiserver multiclass matches pfqn_gld") {
    const Matrix<double> L = mat({{1.2, 0.7}, {0.4, 1.9}});
    const std::vector<int> N = {3, 2};
    Matrix<double> mu(2, 5, 0.0);
    for (int n = 1; n <= 5; ++n) {
        mu(0, n - 1) = std::min(n, 2);
        mu(1, n - 1) = std::min(n, 3);
    }
    agrees_with_gld(L, N, mu);
}

TEST_CASE("pfqn_explicit_ld: fixed-rate multiclass matches pfqn_gld") {
    const Matrix<double> L = mat({{1.0, 0.7, 0.3}, {0.5, 1.3, 0.9}, {0.9, 0.4, 1.7}});
    const std::vector<int> N = {2, 1, 2};
    agrees_with_gld(L, N, ones(3, 5));
}

TEST_CASE("pfqn_explicit_ld: arbitrary limited load dependence matches pfqn_gld") {
    // rates that neither increase nor follow a multi-server shape, settling at s=3
    const Matrix<double> L = mat({{1.1}, {0.6}, {2.3}});
    const std::vector<int> N = {6};
    const Matrix<double> mu = mat({{0.5, 1.4, 2.2, 2.2, 2.2, 2.2},
                                   {1.0, 0.8, 1.7, 1.7, 1.7, 1.7},
                                   {2.0, 1.1, 0.9, 0.9, 0.9, 0.9}});
    agrees_with_gld(L, N, mu);
}

TEST_CASE("pfqn_explicit_ld: rates that never settle are still exact") {
    // an infinite server never satisfies alpha(n)=alpha(s) for a finite s, but
    // s = |N| is admissible because no larger population occurs
    const Matrix<double> L = mat({{0.9, 1.4}, {1.1, 0.5}});
    const std::vector<int> N = {2, 2};
    Matrix<double> mu(2, 4, 0.0);
    for (int n = 1; n <= 4; ++n) {
        mu(0, n - 1) = n;
        mu(1, n - 1) = std::min(n, 2);
    }
    agrees_with_gld(L, N, mu);
}

TEST_CASE("pfqn_explicit_ld: a station a class never visits") {
    // The oracle is the rational convolution over the state space, G = 0.755325,
    // computed offline: it is independent of pfqn_gld, which in the MATLAB, JAR
    // and python ports used to be WRONG on exactly this shape (fixed 2026-09-03,
    // see _kb/07-cross-language-parity.md). THIS port was always right, and the
    // second CHECK keeps it that way.
    const Matrix<double> L = mat({{0.0, 1.3}, {0.9, 0.4}});
    const std::vector<int> N = {2, 2};
    const Matrix<double> mu = ms_rates(2, 4, 2);
    const double lexact = std::log(0.755325);
    const ExplicitResult<double> got = pfqn_explicit_ld(L, N, mu);
    REQUIRE(got.valid);
    CHECK(got.lG == doctest::Approx(lexact).epsilon(1e-9));
    CHECK(pfqn_gld(L, N, mu).lG == doctest::Approx(lexact).epsilon(1e-9));
}

TEST_CASE("pfqn_explicit_ld: repeated scaled demands take Eq. (16)") {
    const Matrix<double> L = mat({{1.3, 0.8}, {1.3, 0.8}, {1.3, 0.8}});
    const std::vector<int> N = {2, 2};
    const Matrix<double> mu = ms_rates(3, 4, 2);
    CHECK(pfqn_explicit_ld(L, N, mu).method == "repeated");
    agrees_with_gld(L, N, mu);
}

TEST_CASE("pfqn_explicit_ld: the fixed-rate degeneration reproduces pfqn_explicit") {
    const Matrix<double> L = mat({{1.0, 0.7}, {0.5, 1.3}, {0.9, 0.4}});
    const std::vector<int> N = {3, 2};
    const ExplicitResult<double> a = pfqn_explicit(L, N);
    const ExplicitResult<double> b = pfqn_explicit_ld(L, N, ones(3, 5));
    CHECK(a.method == b.method);
    CHECK(b.lG == doctest::Approx(a.lG).epsilon(1e-13));
}

TEST_CASE("pfqn_explicit_ld: the single-class route skips the outer sum") {
    // h_theta(N) is homogeneous of degree N in theta, so the divided difference
    // is the identity at R=1 and the constant must still match the convolution
    const Matrix<double> L = mat({{1.4}, {0.9}, {0.35}});
    const std::vector<int> N = {7};
    agrees_with_gld(L, N, ms_rates(3, 7, 2));
}

TEST_CASE("pfqn_explicit_ld: an ulp-wide tie is missed at eps and reported") {
    // both scaled demands are 1.95 at t=[2 3], but land two ulps apart in doubles
    const Matrix<double> L = mat({{0.0, 1.3}, {0.9, 0.7}});
    const std::vector<int> N = {2, 3};
    const Matrix<double> mu = ms_rates(2, 5, 2);
    // the rational convolution over the state space, G = 1.14240375, pinned rather
    // than read off pfqn_gld so this case does not depend on that recursion
    const double lref = std::log(1.14240375);

    const ExplicitResult<double> eps = pfqn_explicit_ld(L, N, mu);
    CHECK(eps.method == "distinct");
    CHECK(eps.lossDigits > 15.0);
    CHECK(std::fabs(eps.lG - lref) > 1.0);

    const ExplicitResult<double> merged = pfqn_explicit_ld(L, N, mu, 1e-12);
    CHECK(merged.method == "repeated");
    CHECK(merged.lG == doctest::Approx(lref).epsilon(1e-9));

    // a caller holding a budget is refused rather than handed the number
    const ExplicitResult<double> budgeted = pfqn_explicit_ld(
        L, N, mu, std::numeric_limits<double>::epsilon(), "auto", 8.0);
    CHECK_FALSE(budgeted.valid);
}

TEST_CASE("pfqn_explicit_ld: inadmissible arguments are refused") {
    const Matrix<double> L = mat({{1.0, 0.5}, {0.5, 1.0}});
    const std::vector<int> N = {2, 2};
    CHECK_THROWS_AS(pfqn_explicit_ld(L, N, ones(2, 4), std::numeric_limits<double>::epsilon(),
                                     "bogus"),
                    line::InputError);
    // a rate lattice shorter than the population cannot answer at |N|
    CHECK_THROWS_AS(pfqn_explicit_ld(L, N, ones(2, 3)), line::InputError);
    // and a rate of zero would divide by zero in phi_k
    Matrix<double> bad = ones(2, 4);
    bad(0, 2) = 0.0;
    CHECK_THROWS_AS(pfqn_explicit_ld(L, N, bad), line::InputError);
}

TEST_CASE("pfqn_ncld: the 'divdiff' token routes through pfqn_explicit_ld") {
    const Matrix<double> L = mat({{1.2, 0.7}, {0.4, 1.9}});
    const std::vector<int> N = {3, 2};
    Matrix<double> mu(2, 5, 0.0);
    for (int n = 1; n <= 5; ++n) {
        mu(0, n - 1) = std::min(n, 2);
        mu(1, n - 1) = std::min(n, 3);
    }
    const Matrix<double> Z(1, 2, 0.0);
    const double lref = pfqn_ncld(L, N, Z, mu, NcldMethod::Exact, 0.0).lG;
    const auto got = pfqn_ncld(L, N, Z, mu, NcldMethod::Divdiff, 0.0);
    CHECK(got.method == "divdiff.ld/distinct");
    CHECK(got.lG == doctest::Approx(lref).epsilon(1e-9));

    // a think time would have to enter g_sigma, whose closed form covers queues only
    Matrix<double> Zt(1, 2, 0.5);
    CHECK_THROWS_AS(pfqn_ncld(L, N, Zt, mu, NcldMethod::Divdiff, 0.0), line::UnsupportedError);
}

TEST_CASE("pfqn_gld: a zero demand on a load-dependent model") {
    // Cross-codebase regression. The MATLAB, JAR and python pfqn_gld dropped a
    // class with jobs and no demand at the ONLY station from the single-station
    // sum, instead of zeroing the constant, and their recursion reaches that base
    // case with the full population every time it peels a station. This port
    // convolves station by station and never had the defect; the constants below
    // are the state-space sums the other three now agree with.
    const std::vector<int> N = {2, 2};
    const Matrix<double> mu = ms_rates(2, 4, 2);
    CHECK(pfqn_gld(mat({{0.0, 1.3}, {0.9, 0.4}}), N, mu).lG ==
          doctest::Approx(std::log(0.755325)).epsilon(1e-9));
    CHECK(pfqn_gld(mat({{0.0, 0.0}, {0.9, 0.4}}), N, mu).lG ==
          doctest::Approx(-2.3309845675200272).epsilon(1e-9));
    CHECK(pfqn_gld(mat({{0.0, 1.3}, {0.9, 0.7}}), std::vector<int>{2, 3},
                   ms_rates(2, 5, 2))
              .lG == doctest::Approx(std::log(1.14240375)).epsilon(1e-9));
    // and the fixed-rate model, which the defect never reached
    CHECK(pfqn_gld(mat({{0.0, 1.3}, {0.9, 0.4}}), N, Matrix<double>(2, 4, 1.0)).lG ==
          doctest::Approx(1.2267416163786375).epsilon(1e-9));
}
