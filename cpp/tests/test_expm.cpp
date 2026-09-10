/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Matrix exponential. The oracles are the cases whose exponential is known in
 * closed form: exp(0) = I, a diagonal matrix (elementwise exp), a nilpotent
 * matrix (a finite polynomial, so the result must be exact to rounding), the
 * 2x2 rotation generator (cos/sin), a Jordan block (t e^t off-diagonal), and
 * the two-state generator. The group law exp(A)exp(-A) = I and the identity
 * exp(A+B) = exp(A)exp(B) for commuting A,B are checked as well, together with
 * agreement between the double and Real50 instantiations.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/util/expm.h"

using line::Matrix;
using line::Real50;
using line::expm;

namespace {

constexpr double TOL = 1e-13;

double maxdiff(const Matrix<double>& A, const Matrix<double>& B) {
    double m = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const double d = std::fabs(A(i, j) - B(i, j));
            if (d > m) m = d;
        }
    return m;
}

}  // namespace

TEST_CASE("expm of the zero matrix is the identity, exactly") {
    Matrix<double> Z(4, 4, 0.0);
    Matrix<double> E = expm(Z);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) CHECK(E(i, j) == (i == j ? 1.0 : 0.0));
}

TEST_CASE("expm of a diagonal matrix exponentiates the diagonal") {
    Matrix<double> D(3, 3, 0.0);
    D(0, 0) = -1.0;
    D(1, 1) = 0.5;
    D(2, 2) = -7.25;
    Matrix<double> E = expm(D);
    CHECK(E(0, 0) == doctest::Approx(std::exp(-1.0)).epsilon(TOL));
    CHECK(E(1, 1) == doctest::Approx(std::exp(0.5)).epsilon(TOL));
    CHECK(E(2, 2) == doctest::Approx(std::exp(-7.25)).epsilon(TOL));
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            if (i != j) CHECK(std::fabs(E(i, j)) < 1e-15);
}

TEST_CASE("expm of a nilpotent matrix is the finite polynomial, to rounding") {
    // N strictly upper triangular 4x4 with N^4 = 0:
    // exp(N) = I + N + N^2/2 + N^3/6.
    Matrix<double> N(4, 4, 0.0);
    N(0, 1) = 1.0;
    N(1, 2) = 1.0;
    N(2, 3) = 1.0;
    Matrix<double> E = expm(N);
    Matrix<double> ref(4, 4, 0.0);
    for (std::size_t i = 0; i < 4; ++i) ref(i, i) = 1.0;
    ref(0, 1) = ref(1, 2) = ref(2, 3) = 1.0;
    ref(0, 2) = ref(1, 3) = 0.5;
    ref(0, 3) = 1.0 / 6.0;
    CHECK(maxdiff(E, ref) < 1e-15);

    // Scaled by 3, still nilpotent: exp(3N) = I + 3N + 9N^2/2 + 27N^3/6.
    Matrix<double> N3 = N;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) N3(i, j) *= 3.0;
    Matrix<double> E3 = expm(N3);
    Matrix<double> ref3(4, 4, 0.0);
    for (std::size_t i = 0; i < 4; ++i) ref3(i, i) = 1.0;
    ref3(0, 1) = ref3(1, 2) = ref3(2, 3) = 3.0;
    ref3(0, 2) = ref3(1, 3) = 4.5;
    ref3(0, 3) = 4.5;
    CHECK(maxdiff(E3, ref3) < 1e-14);
}

TEST_CASE("expm of the rotation generator is the rotation matrix") {
    // A = [[0,-th],[th,0]] -> exp(A) = [[cos th, -sin th],[sin th, cos th]].
    const double th = 0.7;
    Matrix<double> A{{0.0, -th}, {th, 0.0}};
    Matrix<double> E = expm(A);
    CHECK(E(0, 0) == doctest::Approx(std::cos(th)).epsilon(TOL));
    CHECK(E(0, 1) == doctest::Approx(-std::sin(th)).epsilon(TOL));
    CHECK(E(1, 0) == doctest::Approx(std::sin(th)).epsilon(TOL));
    CHECK(E(1, 1) == doctest::Approx(std::cos(th)).epsilon(TOL));
    // A large angle exercises the scaling-and-squaring path.
    const double big = 37.0;
    Matrix<double> B{{0.0, -big}, {big, 0.0}};
    Matrix<double> EB = expm(B);
    CHECK(EB(0, 0) == doctest::Approx(std::cos(big)).epsilon(1e-11));
    CHECK(EB(1, 0) == doctest::Approx(std::sin(big)).epsilon(1e-11));
}

TEST_CASE("expm of a Jordan block matches the closed form") {
    // J = [[a,1],[0,a]] -> exp(J) = e^a [[1,1],[0,1]].
    const double a = -0.75;
    Matrix<double> J{{a, 1.0}, {0.0, a}};
    Matrix<double> E = expm(J);
    CHECK(E(0, 0) == doctest::Approx(std::exp(a)).epsilon(TOL));
    CHECK(E(0, 1) == doctest::Approx(std::exp(a)).epsilon(TOL));
    CHECK(E(1, 0) == doctest::Approx(0.0).epsilon(1e-15));
    CHECK(E(1, 1) == doctest::Approx(std::exp(a)).epsilon(TOL));
}

TEST_CASE("expm of a two-state generator matches the closed form") {
    // Q = [[-a,a],[b,-b]], exp(Qt) = (1/(a+b)) [[b,a],[b,a]]
    //                               + (e^{-(a+b)t}/(a+b)) [[a,-a],[-b,b]].
    const double a = 2.0, b = 3.0, t = 1.3;
    Matrix<double> Q{{-a, a}, {b, -b}};
    Matrix<double> E = expm(Q, t);
    const double s = a + b;
    const double d = std::exp(-s * t);
    CHECK(E(0, 0) == doctest::Approx((b + a * d) / s).epsilon(TOL));
    CHECK(E(0, 1) == doctest::Approx((a - a * d) / s).epsilon(TOL));
    CHECK(E(1, 0) == doctest::Approx((b - b * d) / s).epsilon(TOL));
    CHECK(E(1, 1) == doctest::Approx((a + b * d) / s).epsilon(TOL));
    // Rows of exp(Qt) sum to one for any generator Q.
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(E(i, 0) + E(i, 1) == doctest::Approx(1.0).epsilon(1e-14));
}

TEST_CASE("expm satisfies the group law and the commuting sum rule") {
    Matrix<double> A{{-1.5, 0.4, 0.2}, {1.0, -2.0, 1.0}, {0.3, 0.7, -1.0}};
    Matrix<double> mA = A;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) mA(i, j) = -A(i, j);
    Matrix<double> P = line::matmul(expm(A), expm(mA));
    Matrix<double> I = line::eye<double>(3);
    CHECK(maxdiff(P, I) < 1e-13);

    // exp((s+t)A) = exp(sA) exp(tA), commuting by construction.
    const double s = 0.9, t = 2.4;
    Matrix<double> lhs = expm(A, s + t);
    Matrix<double> rhs = line::matmul(expm(A, s), expm(A, t));
    CHECK(maxdiff(lhs, rhs) < 1e-13);
}

TEST_CASE("expm at Real50 is more accurate than at double") {
    // exp of a diagonal entry is the scalar exponential to the working
    // precision: at 50 digits the error must be far below the double epsilon.
    Matrix<Real50> D(2, 2, Real50(0));
    D(0, 0) = Real50(-1);
    D(1, 1) = Real50(-1);
    Matrix<Real50> E = expm(D);
    const Real50 err = abs(E(0, 0) - exp(Real50(-1)));
    CHECK(static_cast<double>(err) < 1e-40);

    // A nontrivial matrix agrees with the double instantiation to 1e-13.
    Matrix<double> Ad{{-3.0, 3.0}, {1.0, -1.0}};
    Matrix<Real50> Ar(2, 2, Real50(0));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) Ar(i, j) = Real50(Ad(i, j));
    Matrix<double> Ed = expm(Ad);
    Matrix<Real50> Er = expm(Ar);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(std::fabs(Ed(i, j) - static_cast<double>(Er(i, j))) < 1e-13);
}

TEST_CASE("expm rejects malformed input") {
    Matrix<double> R(2, 3, 0.0);
    CHECK_THROWS_AS(expm(R), line::InputError);
    Matrix<double> Z(0, 0);
    CHECK_THROWS_AS(expm(Z), line::InputError);
}
