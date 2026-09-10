/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Exact MVA. Oracles, in order of strength:
 *   1. Little's law and the utilization law, which any correct solution obeys.
 *   2. pfqn_ca, an independent algorithm for the same normalizing constant.
 *   3. Closed-form single-class results.
 * Every case runs in all three arithmetics; the exact instantiation must give
 * the same answer as double to rounding.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mva.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_mva;

namespace {

constexpr double TOL = 1e-9;

template <class T>
Matrix<T> demands_2x2() {
    Matrix<T> L(2, 2);
    L(0, 0) = line::num_traits<T>::from_rational(1, 2);
    L(0, 1) = line::num_traits<T>::from_rational(3, 10);
    L(1, 0) = line::num_traits<T>::from_rational(2, 5);
    L(1, 1) = line::num_traits<T>::from_rational(3, 5);
    return L;
}

template <class T>
Matrix<T> think_2() {
    Matrix<T> Z(1, 2);
    Z(0, 0) = line::num_traits<T>::from_rational(3, 10);
    Z(0, 1) = line::num_traits<T>::from_rational(1, 5);
    return Z;
}

}  // namespace

TEST_CASE("pfqn_mva obeys Little's law and the utilization law") {
    const std::vector<int> N{3, 2};
    auto r = pfqn_mva(demands_2x2<double>(), N, think_2<double>());

    for (std::size_t s = 0; s < 2; ++s) {
        // sum_i Q(i,s) + X(s) Z(s) = N(s)
        double q = 0.0;
        for (std::size_t i = 0; i < 2; ++i) q += r.QN(i, s);
        CHECK(q + r.XN[s] * think_2<double>()(0, s) == doctest::Approx(N[s]).epsilon(TOL));
        // U(i,s) = X(s) L(i,s)
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(r.UN(i, s) == doctest::Approx(r.XN[s] * demands_2x2<double>()(i, s)).epsilon(TOL));
    }
}

TEST_CASE("pfqn_mva normalizing constant matches pfqn_ca") {
    const std::vector<int> N{3, 2};

    auto mva_d = pfqn_mva(demands_2x2<double>(), N, think_2<double>());
    auto ca_d = pfqn_ca(demands_2x2<double>(), N, think_2<double>());
    CHECK(mva_d.lG == doctest::Approx(ca_d.lG).epsilon(TOL));

    auto mva_q = pfqn_mva(demands_2x2<Rational>(), N, think_2<Rational>());
    auto ca_q = pfqn_ca(demands_2x2<Rational>(), N, think_2<Rational>());
    // Both are exact rationals, so they must be equal as values, not merely close.
    CHECK(mva_q.G == ca_q.G);
    CHECK(mva_q.lG == doctest::Approx(ca_d.lG).epsilon(TOL));

    auto mva_r = pfqn_mva(demands_2x2<Real50>(), N, think_2<Real50>());
    CHECK(mva_r.lG == doctest::Approx(ca_d.lG).epsilon(TOL));
}

TEST_CASE("pfqn_mva exact equals double to rounding on every metric") {
    const std::vector<int> N{4, 3};
    auto d = pfqn_mva(demands_2x2<double>(), N, think_2<double>());
    auto q = pfqn_mva(demands_2x2<Rational>(), N, think_2<Rational>());
    for (std::size_t s = 0; s < 2; ++s) {
        CHECK(static_cast<double>(q.XN[s]) == doctest::Approx(d.XN[s]).epsilon(TOL));
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(static_cast<double>(q.QN(i, s)) == doctest::Approx(d.QN(i, s)).epsilon(TOL));
            CHECK(static_cast<double>(q.CN(i, s)) == doctest::Approx(d.CN(i, s)).epsilon(TOL));
            CHECK(static_cast<double>(q.UN(i, s)) == doctest::Approx(d.UN(i, s)).epsilon(TOL));
        }
    }
}

TEST_CASE("pfqn_mva single class single station matches the closed form") {
    // One station, demand D, no delay: X = 1/D and Q = N for any N.
    Matrix<double> L(1, 1);
    L(0, 0) = 0.25;
    const std::vector<int> N{5};
    auto r = pfqn_mva(L, N);
    CHECK(r.XN[0] == doctest::Approx(4.0).epsilon(TOL));
    CHECK(r.QN(0, 0) == doctest::Approx(5.0).epsilon(TOL));
    CHECK(r.UN(0, 0) == doctest::Approx(1.0).epsilon(TOL));

    // Exact: G(N) = D^N, so lG = N log D.
    Matrix<Rational> Lq(1, 1);
    Lq(0, 0) = Rational(1, 4);
    auto rq = pfqn_mva(Lq, N);
    CHECK(rq.G == Rational(1, 1024));
    CHECK(rq.lG == doctest::Approx(5.0 * std::log(0.25)).epsilon(TOL));
}

TEST_CASE("pfqn_mva multiplicity enters the residence time") {
    // mi(i) scales the "one more job" term: C = L (mi + Q).
    Matrix<double> L(1, 1);
    L(0, 0) = 0.5;
    const std::vector<int> N{2};
    auto r1 = pfqn_mva(L, N, Matrix<double>(), std::vector<int>{1});
    auto r2 = pfqn_mva(L, N, Matrix<double>(), std::vector<int>{2});
    CHECK(r2.CN(0, 0) > r1.CN(0, 0));
    CHECK(r1.XN[0] == doctest::Approx(2.0).epsilon(TOL));
}

TEST_CASE("pfqn_mva empty and degenerate populations follow the MATLAB contract") {
    Matrix<double> L(2, 1, 1.0);
    auto z = pfqn_mva(L, std::vector<int>{0});
    CHECK(z.XN[0] == 0.0);
    CHECK(z.lG == 0.0);

    CHECK_THROWS_AS(pfqn_mva(L, std::vector<int>{-1}), line::InputError);
}
