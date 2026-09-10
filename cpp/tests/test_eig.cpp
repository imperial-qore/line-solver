/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * LAPACK-backed eigenvalues and singular values. Oracles are matrices whose
 * spectrum is known in closed form; the subdominant modulus is checked because
 * that is the quantity the NCD machinery consumes.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/util/eig.h"

using line::Matrix;

TEST_CASE("eigenvalues of matrices with a known spectrum") {
    Matrix<double> D(3, 3, 0.0);
    D(0, 0) = 2.0; D(1, 1) = -5.0; D(2, 2) = 0.5;
    std::vector<double> mods;
    for (const std::complex<double>& z : line::eig_values(D)) mods.push_back(std::abs(z));
    std::sort(mods.begin(), mods.end());
    CHECK(mods[0] == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(mods[1] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(mods[2] == doctest::Approx(5.0).epsilon(1e-12));
    CHECK(line::spectral_radius(D) == doctest::Approx(5.0).epsilon(1e-12));
    CHECK(line::subdominant_modulus(D) == doctest::Approx(2.0).epsilon(1e-12));

    // Rotation generator [[0,-1],[1,0]]: eigenvalues +-i, both of modulus 1.
    Matrix<double> R(2, 2, 0.0);
    R(0, 1) = -1.0;
    R(1, 0) = 1.0;
    std::vector<std::complex<double>> e = line::eig_values(R);
    CHECK(std::abs(e[0]) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(std::abs(e[1]) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(std::fabs(e[0].real()) < 1e-12);
}

TEST_CASE("a stochastic matrix has unit spectral radius and a smaller subdominant") {
    Matrix<double> P(3, 3);
    P(0, 0) = 0.5; P(0, 1) = 0.25; P(0, 2) = 0.25;
    P(1, 0) = 1.0 / 3; P(1, 1) = 1.0 / 3; P(1, 2) = 1.0 / 3;
    P(2, 0) = 1.0 / 6; P(2, 1) = 0.5; P(2, 2) = 1.0 / 3;
    CHECK(line::spectral_radius(P) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(line::subdominant_modulus(P) < 1.0);
    CHECK(line::subdominant_modulus(P) > 0.0);
}

TEST_CASE("singular values and rank on a known matrix") {
    // [[3,0],[0,2],[0,0]] has singular values 3 and 2 and rank 2.
    Matrix<double> A(3, 2, 0.0);
    A(0, 0) = 3.0;
    A(1, 1) = 2.0;
    std::vector<double> s = line::svd_values(A);
    REQUIRE(s.size() == 2);
    CHECK(s[0] == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(s[1] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(line::matrix_rank(A) == 2);

    // A rank-deficient matrix: second row is twice the first.
    Matrix<double> B(2, 2);
    B(0, 0) = 1.0; B(0, 1) = 2.0;
    B(1, 0) = 2.0; B(1, 1) = 4.0;
    CHECK(line::matrix_rank(B) == 1);
    CHECK(line::svd_values(B)[1] < 1e-14 * line::svd_values(B)[0]);
    // Frobenius identity: sum of squared singular values equals sum of squares.
    const std::vector<double> sb = line::svd_values(B);
    CHECK(sb[0] * sb[0] + sb[1] * sb[1] == doctest::Approx(1.0 + 4.0 + 4.0 + 16.0).epsilon(1e-12));
}
