/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * api/perm: the permanent of a matrix, by four algorithms.
 *
 * A PERMANENT HAS NO CHEAP INDEPENDENT CHECK, so the verification is that an
 * O(n!) enumeration of the definition and three O(2^n) formulas return the same
 * number on random matrices. That is what these tests do, and it is stronger
 * than any single stored value: a bug that survives it would have to be present
 * in the definition itself.
 *
 * Two hand-computable cases anchor the rest: the all-ones matrix has permanent
 * n!, and a triangular matrix has the product of its diagonal -- the same as its
 * determinant, since only one permutation contributes.
 */
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "doctest.h"
#include "line/api/perm/permanent.h"

namespace perm = line::perm;
using line::Matrix;

TEST_CASE("the all-ones matrix has permanent n!") {
    double fact = 1.0;
    for (std::size_t n = 1; n <= 7; ++n) {
        fact *= static_cast<double>(n);
        const Matrix<double> m(n, n, 1.0);
        CHECK(perm::permanent(m, perm::PermMethod::Naive) == doctest::Approx(fact));
        CHECK(perm::permanent(m, perm::PermMethod::Ryser) == doctest::Approx(fact));
        CHECK(perm::permanent(m, perm::PermMethod::RyserGray) == doctest::Approx(fact));
        CHECK(perm::permanent(m, perm::PermMethod::Multiplicity) == doctest::Approx(fact));
    }
}

TEST_CASE("a triangular matrix has the product of its diagonal") {
    // Only the identity permutation contributes, so the permanent equals the
    // determinant here -- the one family where the two agree.
    const std::size_t n = 5;
    Matrix<double> m(n, n, 0.0);
    double prod = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        m(i, i) = 1.0 + 0.5 * static_cast<double>(i);
        prod *= m(i, i);
        for (std::size_t j = i + 1; j < n; ++j) m(i, j) = 0.3;  // strictly upper
    }
    CHECK(perm::permanent(m, perm::PermMethod::Naive) == doctest::Approx(prod));
    CHECK(perm::permanent(m, perm::PermMethod::Ryser) == doctest::Approx(prod).epsilon(1e-9));
    CHECK(perm::permanent(m, perm::PermMethod::RyserGray) == doctest::Approx(prod).epsilon(1e-9));
}

TEST_CASE("the four algorithms agree on random matrices") {
    std::mt19937_64 g(20260801u);
    std::uniform_real_distribution<double> u(0.1, 2.0);
    for (std::size_t n = 1; n <= 7; ++n) {
        for (int rep = 0; rep < 3; ++rep) {
            Matrix<double> m(n, n, 0.0);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) m(i, j) = u(g);
            const double a = perm::permanent(m, perm::PermMethod::Naive);
            CHECK(perm::permanent(m, perm::PermMethod::Ryser) ==
                  doctest::Approx(a).epsilon(1e-9));
            CHECK(perm::permanent(m, perm::PermMethod::RyserGray) ==
                  doctest::Approx(a).epsilon(1e-9));
            CHECK(perm::permanent(m, perm::PermMethod::Multiplicity) ==
                  doctest::Approx(a).epsilon(1e-9));
        }
    }
}

TEST_CASE("repeated columns are where the default method earns its place") {
    // This is the shape the library actually forms: a class of N jobs
    // contributes N identical columns, so the distinct-column count stays small
    // while n grows.
    const std::size_t n = 8;
    Matrix<double> m(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            m(i, j) = (j < 4) ? (0.5 + 0.1 * static_cast<double>(i))
                              : (1.0 + 0.2 * static_cast<double>(i));

    const double byGray = perm::permanent(m, perm::PermMethod::RyserGray);
    const double byMult = perm::permanent(m, perm::PermMethod::Multiplicity);
    CHECK(byMult == doctest::Approx(byGray).epsilon(1e-9));
    CHECK(byMult > 0.0);

    // Two distinct columns of multiplicity four each: the box product is
    // 5 x 5 = 25 terms against 2^8 = 256 subsets.
    Matrix<double> uniq;
    std::vector<std::size_t> mult;
    perm::permdetail::unique_columns(m, &uniq, &mult);
    REQUIRE(mult.size() == 2u);
    CHECK(mult[0] == 4u);
    CHECK(mult[1] == 4u);
}

TEST_CASE("an empty matrix has permanent one, the empty product") {
    const Matrix<double> e(0, 0, 0.0);
    CHECK(perm::permanent(e, perm::PermMethod::Naive) == doctest::Approx(1.0));
    CHECK(perm::permanent(e, perm::PermMethod::Ryser) == doctest::Approx(1.0));
    CHECK(perm::permanent(e, perm::PermMethod::RyserGray) == doctest::Approx(1.0));
    CHECK(perm::permanent(e, perm::PermMethod::Multiplicity) == doctest::Approx(1.0));
}

TEST_CASE("a zero row makes the permanent zero") {
    Matrix<double> m(4, 4, 1.0);
    for (std::size_t j = 0; j < 4; ++j) m(2, j) = 0.0;
    // Every permutation picks something from row 2, so every term vanishes.
    CHECK(perm::permanent(m, perm::PermMethod::Naive) == doctest::Approx(0.0));
    CHECK(std::fabs(perm::permanent(m, perm::PermMethod::RyserGray)) < 1e-9);
    CHECK(std::fabs(perm::permanent(m, perm::PermMethod::Multiplicity)) < 1e-9);
}

TEST_CASE("snap_to_lattice snaps near-equal columns onto one lattice") {
    // The multiplicity method keys on EXACT column equality, so two demands
    // differing in the last bits are two distinct columns and the saving is
    // lost. Snapping recovers it; it is a deliberate perturbation, which is why
    // it is a separate call rather than something `permanent` does silently.
    Matrix<double> m(3, 4, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        m(i, 0) = 0.5;
        m(i, 1) = 0.5 + 1e-7;  // a hair away from column 0
        m(i, 2) = 1.25;
        m(i, 3) = 1.25 - 1e-7;
    }
    Matrix<double> u0;
    std::vector<std::size_t> k0;
    perm::permdetail::unique_columns(m, &u0, &k0);
    CHECK(k0.size() == 4u);  // four "distinct" columns before snapping

    const Matrix<double> s = perm::snap_to_lattice(m, 0.001);
    Matrix<double> u1;
    std::vector<std::size_t> k1;
    perm::permdetail::unique_columns(s, &u1, &k1);
    CHECK(k1.size() == 2u);  // two after
    CHECK(k1[0] == 2u);
    CHECK(k1[1] == 2u);

    CHECK_THROWS_AS(perm::snap_to_lattice(m, 0.0), line::InputError);
}

TEST_CASE("the refusals are by name") {
    const Matrix<double> nonsq(2, 3, 1.0);
    CHECK_THROWS_AS(perm::permanent(nonsq, perm::PermMethod::Naive), line::InputError);
    CHECK_THROWS_AS(perm::permanent(nonsq, perm::PermMethod::Ryser), line::InputError);
    CHECK_THROWS_AS(perm::permanent(nonsq, perm::PermMethod::RyserGray), line::InputError);
    CHECK_THROWS_AS(perm::permanent(nonsq, perm::PermMethod::Multiplicity), line::InputError);

    // The enumerations refuse a size they cannot finish rather than hanging.
    const Matrix<double> huge(13, 13, 1.0);
    CHECK_THROWS_AS(perm::permanent(huge, perm::PermMethod::Naive), line::InputError);
    const Matrix<double> vast(31, 31, 1.0);
    CHECK_THROWS_AS(perm::permanent(vast, perm::PermMethod::Ryser), line::InputError);
    CHECK_THROWS_AS(perm::permanent(vast, perm::PermMethod::RyserGray), line::InputError);
}
