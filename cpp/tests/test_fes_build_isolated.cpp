/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * fes_build_isolated: demands and visit ratios of an isolated subnetwork.
 *
 * ORACLES.
 *  (a) MATLAB, on a three-station two-class subset with hand-built routing
 *      blocks, to 15 digits.
 *  (b) The defining fixed point, exactly at Rational: the visit vector of each
 *      class must satisfy v P = v and sum to one with identically zero
 *      residual. That is stronger than agreement with the reference, which
 *      solves the same system in double.
 *  (c) The MVA normalization, L(i,k) = v(i,k) / (rate(i,k) v(0,k)), checked
 *      entry by entry against the visits the same call returned.
 *
 * The reference is driven through an `sn` stub carrying only the five fields
 * it reads; the port takes the extracted rates directly, so the comparison is
 * of the algorithm and not of the extraction.
 */
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/fes/fes_build_isolated.h"

using line::Matrix;
using line::Rational;
using namespace line::fes;

namespace {

/** The class-1 and class-2 routing blocks of the MATLAB fixture. */
template <class T>
Matrix<T> fixture_S() {
    const std::size_t M = 3, K = 2;
    const Matrix<T> P1{{line::num_traits<T>::from_int(0), line::num_traits<T>::from_rational(3, 5),
                        line::num_traits<T>::from_rational(2, 5)},
                       {line::num_traits<T>::from_rational(1, 2), line::num_traits<T>::from_int(0),
                        line::num_traits<T>::from_rational(1, 2)},
                       {line::num_traits<T>::from_rational(3, 10),
                        line::num_traits<T>::from_rational(7, 10), line::num_traits<T>::from_int(0)}};
    const Matrix<T> P2{{line::num_traits<T>::from_rational(1, 10),
                        line::num_traits<T>::from_rational(1, 2),
                        line::num_traits<T>::from_rational(2, 5)},
                       {line::num_traits<T>::from_rational(1, 5),
                        line::num_traits<T>::from_rational(3, 10),
                        line::num_traits<T>::from_rational(1, 2)},
                       {line::num_traits<T>::from_rational(3, 5),
                        line::num_traits<T>::from_rational(1, 10),
                        line::num_traits<T>::from_rational(3, 10)}};
    Matrix<T> S(M * K, M * K, line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            S(i * K + 0, j * K + 0) = P1(i, j);
            S(i * K + 1, j * K + 1) = P2(i, j);
        }
    return S;
}

/** rates = [2 4; 1 1/2; 5 5/2], the MATLAB fixture. */
template <class T>
Matrix<T> fixture_rates() {
    Matrix<T> r(3, 2);
    r(0, 0) = line::num_traits<T>::from_int(2);
    r(0, 1) = line::num_traits<T>::from_int(4);
    r(1, 0) = line::num_traits<T>::from_int(1);
    r(1, 1) = line::num_traits<T>::from_rational(1, 2);
    r(2, 0) = line::num_traits<T>::from_int(5);
    r(2, 1) = line::num_traits<T>::from_rational(5, 2);
    return r;
}

}  // namespace

TEST_CASE("fes_build_isolated reproduces the MATLAB demands and visits") {
    // Oracle (a). MATLAB fes_build_isolated(sn, [1 2 3], S) on the fixture.
    const FesIsolated<double> r = fes_build_isolated(fixture_rates<double>(), fixture_S<double>());

    const double refV[3][2] = {{0.291479820627803, 0.323529411764706},
                               {0.394618834080718, 0.286764705882353},
                               {0.313901345291480, 0.389705882352941}};
    const double refL[3][2] = {{0.5, 0.25},
                               {1.35384615384615, 1.77272727272727},
                               {0.215384615384615, 0.481818181818182}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
            INFO("i=", i, " k=", k);
            CHECK(r.visits(i, k) == doctest::Approx(refV[i][k]).epsilon(1e-12));
            CHECK(r.L(i, k) == doctest::Approx(refL[i][k]).epsilon(1e-12));
        }
}

TEST_CASE("the visit ratios are the exact stationary vector of each class") {
    // Oracle (b): v P = v and sum(v) = 1, with identically zero residual at
    // Rational. The reference can only satisfy this to its double solve.
    const Matrix<Rational> S = fixture_S<Rational>();
    const FesIsolated<Rational> r = fes_build_isolated(fixture_rates<Rational>(), S);
    const std::size_t M = 3, K = 2;

    for (std::size_t k = 0; k < K; ++k) {
        Matrix<Rational> P(M, M);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) P(i, j) = S(i * K + k, j * K + k);

        Rational total(0);
        for (std::size_t i = 0; i < M; ++i) total += r.visits(i, k);
        CHECK(total == Rational(1));  // exact

        for (std::size_t j = 0; j < M; ++j) {
            Rational vp(0);
            for (std::size_t i = 0; i < M; ++i) vp += r.visits(i, k) * P(i, j);
            CHECK(vp == r.visits(j, k));  // exact fixed point, zero residual
        }
    }

    // Oracle (c): the MVA normalization, entry by entry.
    const Matrix<Rational> rates = fixture_rates<Rational>();
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i)
            CHECK(r.L(i, k) == r.visits(i, k) / rates(i, k) / r.visits(0, k));  // exact
}

TEST_CASE("a disabled service yields a zero demand") {
    // MATLAB with sn.rates(3,2) = 0 returns L(3,2) = 0 and leaves the rest of
    // the column unchanged; the visits are untouched because they come from
    // the routing alone.
    Matrix<double> rates = fixture_rates<double>();
    rates(2, 1) = 0.0;
    const FesIsolated<double> r = fes_build_isolated(rates, fixture_S<double>());
    CHECK(r.L(2, 1) == doctest::Approx(0.0));
    CHECK(r.L(0, 1) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(r.L(1, 1) == doctest::Approx(1.77272727272727).epsilon(1e-12));
    // and the visits are unchanged by the rate
    CHECK(r.visits(2, 1) == doctest::Approx(0.389705882352941).epsilon(1e-12));

    SUBCASE("a non-finite rate is disabled too") {
        Matrix<double> bad = fixture_rates<double>();
        bad(1, 0) = std::numeric_limits<double>::infinity();
        const FesIsolated<double> q = fes_build_isolated(bad, fixture_S<double>());
        CHECK(q.L(1, 0) == doctest::Approx(0.0));
    }
}

TEST_CASE("a class with no routing falls back to the uniform visit vector") {
    // This is the branch where the port and the reference reach the same
    // answer by different routes. Zeroing the class-2 block makes every row
    // sum fall below FineTol, so each row is repaired to a self-loop and
    // P_2 becomes the identity; then A = I - I + e e'/M = e e'/M has rank one.
    // MATLAB detects that with rank(A) ~= M and substitutes the uniform
    // distribution; the port's LU hits an exactly zero pivot and takes the
    // same fallback. MATLAB returns [1/3 1/3 1/3] and so does the port.
    Matrix<Rational> S = fixture_S<Rational>();
    const std::size_t M = 3, K = 2;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) S(i * K + 1, j * K + 1) = Rational(0);

    const FesIsolated<Rational> r = fes_build_isolated(fixture_rates<Rational>(), S);
    for (std::size_t i = 0; i < M; ++i) CHECK(r.visits(i, 1) == Rational(1, 3));  // exact
    // class 1 is untouched
    CHECK(line::num_traits<Rational>::to_double(r.visits(0, 0)) ==
          doctest::Approx(0.291479820627803).epsilon(1e-12));
}

TEST_CASE("a row that does not sum to one is renormalized") {
    // The reference divides any row whose sum differs from one by more than
    // FineTol. Scaling the whole class-1 block by 2 must therefore leave the
    // visit ratios unchanged.
    Matrix<Rational> S = fixture_S<Rational>();
    const Matrix<Rational> base = S;
    const std::size_t M = 3, K = 2;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) S(i * K + 0, j * K + 0) *= Rational(2);

    const FesIsolated<Rational> scaled = fes_build_isolated(fixture_rates<Rational>(), S);
    const FesIsolated<Rational> plain = fes_build_isolated(fixture_rates<Rational>(), base);
    for (std::size_t i = 0; i < M; ++i) CHECK(scaled.visits(i, 0) == plain.visits(i, 0));  // exact
}

TEST_CASE("fes_build_isolated rejects malformed input") {
    const Matrix<double> S = fixture_S<double>();
    CHECK_THROWS_AS(fes_build_isolated(Matrix<double>(0, 0), S), line::InputError);
    CHECK_THROWS_AS(fes_build_isolated(Matrix<double>(3, 0), S), line::InputError);
    // the stochastic complement must cover M_sub * K
    CHECK_THROWS_AS(fes_build_isolated(fixture_rates<double>(), Matrix<double>(4, 4, 0.0)),
                    line::InputError);
}

TEST_CASE("fes_build_isolated instantiates at Real50") {
    using R50 = line::Real50;
    const FesIsolated<R50> r = fes_build_isolated(fixture_rates<R50>(), fixture_S<R50>());
    R50 tot = line::num_traits<R50>::from_int(0);
    for (std::size_t i = 0; i < 3; ++i) tot += r.visits(i, 0);
    CHECK(line::num_traits<R50>::to_double(tot) == doctest::Approx(1.0).epsilon(1e-30));
    CHECK(line::num_traits<R50>::to_double(r.L(0, 0)) == doctest::Approx(0.5).epsilon(1e-30));
}
