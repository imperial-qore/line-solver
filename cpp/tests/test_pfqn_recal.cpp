/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Tests for the RECAL port. Three kinds of check:
 *   - agreement with the convolution algorithm on the same model, in all three
 *     arithmetics, and bit-for-bit rational equality in the exact one;
 *   - agreement with a direct enumeration of the product-form state space,
 *     itself carried out in exact arithmetic and independent of both
 *     algorithms;
 *   - the edge cases of the MATLAB contract (empty population, single class,
 *     no think time) and the station-multiplicity argument.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_recal.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_recal;

namespace {

constexpr double EXACT_TOL = 1e-10;

/**
 * Product-form balance function of one queueing station holding cnt[r] jobs of
 * class r: (sum_r cnt_r)! / prod_r cnt_r! * prod_r L(j,r)^cnt_r.
 */
Rational station_term(const Matrix<Rational>& L, std::size_t j, const std::vector<int>& cnt) {
    unsigned tot = 0;
    for (int c : cnt) tot += static_cast<unsigned>(c);
    Rational f = line::num_factorial<Rational>(tot);
    for (std::size_t r = 0; r < cnt.size(); ++r) {
        f /= line::num_factorial<Rational>(static_cast<unsigned>(cnt[r]));
        f *= line::num_pow_int(L(j, r), static_cast<unsigned>(cnt[r]));
    }
    return f;
}

/** Sum of prod_j station_term over every allocation of rem to stations j..M-1. */
Rational station_sum(const Matrix<Rational>& L, std::size_t j, const std::vector<int>& rem) {
    if (j + 1 == L.rows()) return station_term(L, j, rem);
    const std::size_t R = rem.size();
    Rational s(0);
    std::vector<int> c(R, 0);
    bool more = true;
    while (more) {
        std::vector<int> rest(R);
        for (std::size_t r = 0; r < R; ++r) rest[r] = rem[r] - c[r];
        s += station_term(L, j, c) * station_sum(L, j + 1, rest);
        more = line::next_pop(c, rem);
    }
    return s;
}

/**
 * G(N) by direct enumeration of the state space: every split of the population
 * between the delay and the queueing stations, and every allocation of the
 * remainder over the stations. Exact, and independent of RECAL and of
 * convolution alike.
 */
Rational enumerate_g(const Matrix<Rational>& L, const std::vector<int>& N,
                     const std::vector<Rational>& Z) {
    const std::size_t R = N.size();
    Rational g(0);
    std::vector<int> d(R, 0);
    bool more = true;
    while (more) {
        Rational zterm(1);
        for (std::size_t r = 0; r < R; ++r) {
            zterm *= line::num_pow_int(Z[r], static_cast<unsigned>(d[r]));
            zterm /= line::num_factorial<Rational>(static_cast<unsigned>(d[r]));
        }
        if (zterm != Rational(0)) {
            std::vector<int> rem(R);
            for (std::size_t r = 0; r < R; ++r) rem[r] = N[r] - d[r];
            g += zterm * station_sum(L, 0, rem);
        }
        more = line::next_pop(d, N);
    }
    return g;
}

/** Three stations, two classes, both classes thinking. */
template <class T>
void model_with_think(Matrix<T>& L, Matrix<T>& Z) {
    L = Matrix<T>(3, 2);
    L(0, 0) = line::num_traits<T>::from_rational(1, 2);
    L(0, 1) = line::num_traits<T>::from_rational(3, 10);
    L(1, 0) = line::num_traits<T>::from_rational(2, 5);
    L(1, 1) = line::num_traits<T>::from_rational(3, 5);
    L(2, 0) = line::num_traits<T>::from_rational(1, 4);
    L(2, 1) = line::num_traits<T>::from_rational(7, 10);
    Z = Matrix<T>(1, 2);
    Z(0, 0) = line::num_traits<T>::from_rational(3, 10);
    Z(0, 1) = line::num_traits<T>::from_rational(1, 5);
}

/** The same demands with no delay, exercising the Mz = M code path. */
template <class T>
void model_no_think(Matrix<T>& L) {
    Matrix<T> Z;
    model_with_think(L, Z);
}

}  // namespace

TEST_CASE("pfqn_recal agrees with pfqn_ca with think time, all arithmetics") {
    const std::vector<int> N{2, 3};

    Matrix<double> Ld, Zd;
    model_with_think(Ld, Zd);
    const auto rd = pfqn_recal(Ld, N, Zd);
    const auto cd = pfqn_ca(Ld, N, Zd);
    CHECK(rd.lG == doctest::Approx(cd.lG).epsilon(EXACT_TOL));
    CHECK(rd.G == doctest::Approx(cd.G).epsilon(EXACT_TOL));

    Matrix<Real50> Lr, Zr;
    model_with_think(Lr, Zr);
    const auto rr = pfqn_recal(Lr, N, Zr);
    const auto cr = pfqn_ca(Lr, N, Zr);
    CHECK(rr.lG == doctest::Approx(cr.lG).epsilon(EXACT_TOL));

    // Exact arithmetic: the two algorithms must return the identical rational,
    // not merely the same value to rounding.
    Matrix<Rational> Lq, Zq;
    model_with_think(Lq, Zq);
    const auto rq = pfqn_recal(Lq, N, Zq);
    const auto cq = pfqn_ca(Lq, N, Zq);
    CHECK(rq.G == cq.G);
    CHECK(rq.lG == doctest::Approx(cq.lG).epsilon(1e-14));
    CHECK(rd.lG == doctest::Approx(rq.lG).epsilon(EXACT_TOL));
    CHECK(static_cast<double>(rr.lG) == doctest::Approx(rq.lG).epsilon(1e-14));
}

TEST_CASE("pfqn_recal agrees with pfqn_ca without think time, all arithmetics") {
    const std::vector<int> N{3, 2};

    Matrix<double> Ld;
    model_no_think(Ld);
    CHECK(pfqn_recal(Ld, N).lG == doctest::Approx(pfqn_ca(Ld, N).lG).epsilon(EXACT_TOL));

    Matrix<Real50> Lr;
    model_no_think(Lr);
    CHECK(pfqn_recal(Lr, N).lG == doctest::Approx(pfqn_ca(Lr, N).lG).epsilon(EXACT_TOL));

    Matrix<Rational> Lq;
    model_no_think(Lq);
    CHECK(pfqn_recal(Lq, N).G == pfqn_ca(Lq, N).G);
}

TEST_CASE("pfqn_recal matches a direct enumeration of the state space") {
    const std::vector<int> N{2, 3};

    Matrix<Rational> Lq, Zq;
    model_with_think(Lq, Zq);
    std::vector<Rational> Z{Zq(0, 0), Zq(0, 1)};
    const Rational expected = enumerate_g(Lq, N, Z);
    CHECK(pfqn_recal(Lq, N, Zq).G == expected);

    // Same demands with no delay: the enumeration collapses to the station sum.
    std::vector<Rational> Z0{Rational(0), Rational(0)};
    const Rational expected0 = enumerate_g(Lq, N, Z0);
    CHECK(pfqn_recal(Lq, N).G == expected0);

    // A single-class instance, where the enumeration is a plain composition sum.
    Matrix<Rational> L1(2, 1);
    L1(0, 0) = Rational(1, 2);
    L1(1, 0) = Rational(1, 3);
    const std::vector<int> N1{5};
    const Rational e1 = enumerate_g(L1, N1, std::vector<Rational>{Rational(0)});
    CHECK(pfqn_recal(L1, N1).G == e1);
}

TEST_CASE("pfqn_recal station multiplicity replicates a station exactly") {
    // Two identical stations plus a third, unit multiplicities.
    Matrix<Rational> Lexp(3, 2);
    Lexp(0, 0) = Rational(1, 2);
    Lexp(0, 1) = Rational(3, 10);
    Lexp(1, 0) = Rational(1, 2);
    Lexp(1, 1) = Rational(3, 10);
    Lexp(2, 0) = Rational(2, 5);
    Lexp(2, 1) = Rational(3, 5);
    const std::vector<int> N{2, 2};

    // The same network written with two stations and a multiplicity of two.
    Matrix<Rational> Lcon(2, 2);
    Lcon(0, 0) = Rational(1, 2);
    Lcon(0, 1) = Rational(3, 10);
    Lcon(1, 0) = Rational(2, 5);
    Lcon(1, 1) = Rational(3, 5);
    const std::vector<int> m0{2, 1};

    const Rational expanded = pfqn_ca(Lexp, N).G;  // convolution never consolidates
    CHECK(pfqn_recal(Lexp, N).G == expanded);
    CHECK(pfqn_recal(Lcon, N, Matrix<Rational>(), m0).G == expanded);
    CHECK(pfqn_recal(Lexp, N).G == enumerate_g(Lexp, N, std::vector<Rational>{Rational(0), Rational(0)}));
}

TEST_CASE("pfqn_recal edge cases follow the MATLAB contract") {
    Matrix<double> L(2, 1, 1.0);

    // Zero population: G = 1, lG = 0.
    CHECK(pfqn_recal(L, std::vector<int>{0}).G == 1.0);
    CHECK(pfqn_recal(L, std::vector<int>{0}).lG == 0.0);

    // Zero population with think time present is still 1.
    Matrix<double> Z(1, 1);
    Z(0, 0) = 4.0;
    CHECK(pfqn_recal(L, std::vector<int>{0}, Z).G == 1.0);

    // Single class, no think time: two unit-demand stations give G = N + 1.
    CHECK(pfqn_recal(L, std::vector<int>{4}).G == doctest::Approx(5.0));

    // Single station, single class, no think time: G = L^N.
    Matrix<Rational> Ls(1, 1);
    Ls(0, 0) = Rational(1, 2);
    CHECK(pfqn_recal(Ls, std::vector<int>{3}).G == Rational(1, 8));

    // Delay-only network (M = 0): G = Z^N / N!.
    Matrix<double> noL;
    Matrix<double> Zd(1, 1);
    Zd(0, 0) = 2.0;
    CHECK(pfqn_recal(noL, std::vector<int>{3}, Zd).G == doctest::Approx(8.0 / 6.0));

    // A class with zero population contributes nothing.
    Matrix<Rational> L2(2, 2);
    L2(0, 0) = Rational(1, 2);
    L2(0, 1) = Rational(1, 3);
    L2(1, 0) = Rational(1, 5);
    L2(1, 1) = Rational(1, 7);
    Matrix<Rational> L2a(2, 1);
    L2a(0, 0) = Rational(1, 2);
    L2a(1, 0) = Rational(1, 5);
    CHECK(pfqn_recal(L2, std::vector<int>{3, 0}).G == pfqn_recal(L2a, std::vector<int>{3}).G);
}

TEST_CASE("pfqn_recal rejects malformed input") {
    Matrix<double> L(2, 2, 1.0);
    CHECK_THROWS_AS(pfqn_recal(L, std::vector<int>{1}), line::InputError);
    CHECK_THROWS_AS(pfqn_recal(L, std::vector<int>{-1, 1}), line::InputError);
    CHECK_THROWS_AS(pfqn_recal(L, std::vector<int>{1, 1}, Matrix<double>(), std::vector<int>{1}),
                    line::InputError);
    CHECK_THROWS_AS(pfqn_recal(L, std::vector<int>{1, 1}, Matrix<double>(), std::vector<int>{0, 1}),
                    line::InputError);
}

TEST_CASE("pfqn_recal scaling keeps lG finite where the unscaled recursion overflows") {
    Matrix<double> L(2, 1);
    L(0, 0) = 1e30;
    L(1, 0) = 1e30;
    const auto r = pfqn_recal(L, std::vector<int>{40});
    CHECK(std::isfinite(r.lG));
    CHECK(r.lG == doctest::Approx(pfqn_ca(L, std::vector<int>{40}).lG).epsilon(1e-12));

    // The exact path needs no scaling and must agree with the scaled double.
    const line::BigInt e30("1000000000000000000000000000000");
    Matrix<Rational> Lq(2, 1);
    Lq(0, 0) = Rational(e30);
    Lq(1, 0) = Lq(0, 0);
    const auto rq = pfqn_recal(Lq, std::vector<int>{40});
    CHECK(rq.G == pfqn_ca(Lq, std::vector<int>{40}).G);
    CHECK(rq.lG == doctest::Approx(r.lG).epsilon(1e-12));
}
