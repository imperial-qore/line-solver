/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * perm_heur and perm_bethe: approximate permanents.
 *
 * Every golden was measured in native Python on the same inputs. The two
 * approximations carry DIFFERENT guarantees and the tests keep them apart:
 * `perm_heur` has no error bound in either direction, while `perm_bethe` is a
 * LOWER bound for a nonnegative matrix -- so only the second is checked against
 * the exact value as a bound, and the first only against its own reference.
 *
 * THE ALL-ONES FAMILY IS WHERE THE HEURISTIC IS EXACT: the matrix is already
 * doubly stochastic after scaling, so the mean-field and capacity estimates
 * coincide with n!. That is a useful anchor precisely because it is the one
 * case where a heuristic can be held to an exact number.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/perm/perm_approx.h"
#include "line/api/perm/permanent.h"
#include "line/api/perm/perm_sampling.h"

namespace perm = line::perm;
using line::Matrix;

TEST_CASE("the heuristic is exact on the all-ones matrix") {
    double fact = 1.0;
    for (std::size_t n = 2; n <= 6; ++n) {
        fact *= static_cast<double>(n);
        const Matrix<double> m(n, n, 1.0);
        // n! for n = 2..6: 2, 6, 24, 120, 720. Native Python agrees.
        CHECK(perm::perm_heur(m) == doctest::Approx(fact).epsilon(1e-6));
    }
}

TEST_CASE("the Bethe estimate matches native Python and stays a lower bound") {
    // Native Python on the all-ones family.
    const double golden[5] = {1.0, 2.37037, 8.10915, 36.0288, 196.549};
    double fact = 1.0;
    for (std::size_t n = 2; n <= 6; ++n) {
        fact *= static_cast<double>(n);
        const Matrix<double> m(n, n, 1.0);
        const double b = perm::perm_bethe(m);
        CHECK(b == doctest::Approx(golden[n - 2]).epsilon(1e-5));
        // A LOWER bound: this is the guarantee the heuristic does not carry.
        CHECK(b <= fact + 1e-9);
        CHECK(b > 0.0);
    }
}

TEST_CASE("a dense asymmetric matrix, where the two diverge") {
    Matrix<double> m(4, 4, 0.0);
    const double v[4][4] = {{1, 2, 3, 4}, {2, 1, 4, 3}, {3, 4, 1, 2}, {4, 3, 2, 1}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) m(i, j) = v[i][j];

    const double exact = perm::permanent(m, perm::PermMethod::Naive);
    CHECK(exact == doctest::Approx(1092.0));
    // Native Python, same matrix.
    CHECK(perm::perm_heur(m) == doctest::Approx(937.5).epsilon(1e-5));
    CHECK(perm::perm_bethe(m) == doctest::Approx(325.717).epsilon(1e-4));
    // The Bethe value is a bound; the heuristic happens to be below here but
    // is NOT guaranteed to be, which is why only the first is asserted as one.
    CHECK(perm::perm_bethe(m) < exact);
}

TEST_CASE("an exact zero is REFUSED by both approximations") {
    // This test used to pin the floored behaviour, and the numbers it pinned
    // are why that behaviour was removed. On this matrix (exact permanent 3):
    //
    //   perm_heur  returned 2.4644822097  -- the van der Waerden bound of the
    //              Sinkhorn limit of the FLOORED matrix, not of this one;
    //   perm_bethe returned 2748880111.1  -- nine orders of magnitude above an
    //              exact 3, with the lower-bound property gone.
    //
    // Both came from nudging a zero to a small eps before scaling. That
    // substitution is not invertible: every permutation takes one entry per
    // row, so a floored matrix has permanent n! eps times the permanent of the
    // rest where the truth may be 0, and n! outruns eps by n = 18. The old
    // comment here already said what a caller needs to know -- that these are
    // meaningful only on a STRICTLY POSITIVE matrix -- so the contract is now
    // enforced instead of documented.
    Matrix<double> m(3, 3, 1.0);
    m(0, 0) = 0.0;
    m(2, 1) = 0.0;
    const double exact = perm::permanent(m, perm::PermMethod::Naive);
    CHECK(exact == doctest::Approx(3.0));

    CHECK_THROWS_AS(perm::perm_heur(m), line::InputError);
    CHECK_THROWS_AS(perm::perm_bethe(m), line::InputError);

    // The exact engine is correct on zeros and remains the documented route.
    CHECK(perm::permanent(m, perm::PermMethod::Ryser) == doctest::Approx(3.0));
}

TEST_CASE("Sinkhorn scaling really does approach double stochasticity") {
    Matrix<double> m(4, 4, 0.0);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            m(i, j) = 1.0 + static_cast<double>(i) + 2.0 * static_cast<double>(j);

    Matrix<double> B;
    std::vector<double> r, c;
    perm::sinkhorn_scaling(m, &B, &r, &c);
    for (std::size_t i = 0; i < 4; ++i) {
        double rs = 0.0, cs = 0.0;
        for (std::size_t j = 0; j < 4; ++j) {
            rs += B(i, j);
            cs += B(j, i);
        }
        CHECK(rs == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(cs == doctest::Approx(1.0).epsilon(1e-6));
    }
}

TEST_CASE("a negative entry is refused, with its position named") {
    Matrix<double> m(3, 3, 1.0);
    m(1, 2) = -0.5;
    CHECK_THROWS_AS(perm::perm_heur(m), line::InputError);
    CHECK_THROWS_AS(perm::perm_bethe(m), line::InputError);

    const Matrix<double> nonsq(2, 3, 1.0);
    CHECK_THROWS_AS(perm::perm_heur(nonsq), line::InputError);
    CHECK_THROWS_AS(perm::perm_bethe(nonsq), line::InputError);

    const Matrix<double> empty(0, 0, 0.0);
    CHECK_THROWS_AS(perm::perm_heur(empty), line::InputError);
    CHECK_THROWS_AS(perm::perm_bethe(empty), line::InputError);
}

// ---------------------------------------------------------------------------
// Zeros in the demand matrix: the approximations require full support.
// An eps floor is not invertible -- perm(max(A,eps)) = n! eps perm(rest)
// against a true permanent of 0, and n! outruns eps by n = 18.
// ---------------------------------------------------------------------------

TEST_CASE("the approximations refuse a matrix with a structural zero") {
    const Matrix<double> zero_row({{1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, {4.0, 5.0, 6.0}});
    CHECK(perm::permanent(zero_row, perm::PermMethod::Ryser) == doctest::Approx(0.0));
    CHECK_THROWS_AS(perm::perm_bethe(zero_row), line::InputError);
    CHECK_THROWS_AS(perm::perm_heur(zero_row), line::InputError);
    CHECK_THROWS_AS(perm::AdaPartSampler{zero_row}, line::InputError);
    CHECK_THROWS_AS(perm::HuberLawSampler{zero_row}, line::InputError);
}

TEST_CASE("the approximations terminate on a matrix with a positive permanent") {
    // [[1,2],[0,3]] has permanent 3 > 0 but no total support. The samplers used
    // to spin here: a zero row sum makes the guarded Sinkhorn normalization a
    // no-op, the margin error never falls, and the loop had no cap.
    const Matrix<double> a({{1.0, 2.0}, {0.0, 3.0}});
    CHECK(perm::permanent(a, perm::PermMethod::Ryser) == doctest::Approx(3.0));
    CHECK_THROWS_AS(perm::perm_bethe(a), line::InputError);
    CHECK_THROWS_AS(perm::perm_heur(a), line::InputError);
    CHECK_THROWS_AS(perm::AdaPartSampler{a}, line::InputError);
    CHECK_THROWS_AS(perm::HuberLawSampler{a}, line::InputError);
}

TEST_CASE("the maximum weight assignment is not the row-by-row greedy") {
    // The greedy takes the larger entry of row 0 and leaves row 1 with the
    // smaller one. That weight is alpha3, which sets the Huber-Law flooring
    // level alpha1, so the optimum is what the method's guarantee needs.
    const Matrix<double> w({{std::log(1.0), std::log(2.0)},
                            {std::log(1e-9), std::log(3.0)}});
    const std::vector<std::size_t> a = perm::samplingdetail::max_weight_assignment(w);
    CHECK(a[0] == 0u);
    CHECK(a[1] == 1u);
    const double optimal = w(0, a[0]) + w(1, a[1]);
    const double greedy = w(0, 1) + w(1, 0);
    CHECK(optimal > greedy);
}

TEST_CASE("the saddle point is exact on a single column") {
    // h == 1 leaves no direction after the homogeneity is quotiented out, so the
    // Laplace factor is empty and the expansion returns the exact n! prod_k a_k0.
    for (std::size_t n = 1; n <= 6; ++n) {
        Matrix<double> col(n, 1, 0.0);
        double product = 1.0, fact = 1.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double v = 0.1 + 0.13 * static_cast<double>(i + 1);
            col(i, 0) = v;
            product *= v;
        }
        for (std::size_t i = 2; i <= n; ++i) fact *= static_cast<double>(i);
        CHECK(perm::perm_spm(col, std::vector<std::size_t>(1, n)) ==
              doctest::Approx(fact * product).epsilon(1e-12));
    }
}

TEST_CASE("the saddle point matches its closed form on the all-ones matrix") {
    // J_n scales to itself (xi = 1), so phi = n log n, the reduced Laplacian
    // I - J/n has determinant 1/n, and the estimate is (2 pi)^(-(n-1)/2)
    // n^(n+1/2). MATLAB, the JAR and native Python are held to the same anchor.
    for (std::size_t n = 2; n <= 10; ++n) {
        const Matrix<double> m(n, n, 1.0);
        const double nd = static_cast<double>(n);
        const double expected = std::pow(2.0 * M_PI, -0.5 * (nd - 1.0)) * std::pow(nd, nd + 0.5);
        CHECK(perm::perm_spm(m) == doctest::Approx(expected).epsilon(1e-9));
    }
}

TEST_CASE("the saddle point lands nearer than the capacity it corrects") {
    // log_capacity is the Gurvits capacity, an upper bound of the permanent. The
    // Gaussian factor is what turns that e^n-scale bound into a usable estimate.
    unsigned seed = 4711u;
    auto rnd = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return static_cast<double>((seed >> 16) & 0x7fff) / 32768.0;
    };
    for (std::size_t n = 3; n <= 7; ++n) {
        Matrix<double> m(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) m(i, j) = rnd() + 0.05;
        const double reference = perm::permanent(m);
        const perm::PermSpmResult r = perm::perm_spm_expand(m, std::vector<std::size_t>());
        const double capacity = std::exp(r.log_capacity);
        CHECK(capacity >= reference * (1.0 - 1e-9));
        CHECK(std::fabs(r.value - reference) < std::fabs(capacity - reference));
        CHECK(r.log_value == doctest::Approx(std::log(r.value)).epsilon(1e-12));
        // measured bias at unit multiplicities, near (e/sqrt(2 pi))^n
        const double bias = std::pow(std::exp(1.0) / std::sqrt(2.0 * M_PI), static_cast<double>(n));
        CHECK(r.value >= reference);
        CHECK(r.value <= reference * bias * 1.1);
    }
}

TEST_CASE("the saddle point sharpens as the column multiplicities grow") {
    // With h fixed this is a genuine asymptotic expansion in min(m). The matrix
    // whose 2k rows all equal (a, b) with m = (k, k) has permanent (2k)! (ab)^k,
    // and the expansion returns exactly 4^k / (C(2k,k) sqrt(pi k)) times it: an
    // analytic ratio, free of the matrix, decreasing to 1 like 1 + 1/(8k).
    const double a = 0.7, b = 1.3;
    double previous = std::numeric_limits<double>::infinity();
    for (std::size_t k = 1; k <= 8; ++k) {
        Matrix<double> rows(2 * k, 2, 0.0);
        for (std::size_t i = 0; i < 2 * k; ++i) {
            rows(i, 0) = a;
            rows(i, 1) = b;
        }
        double exact = 1.0;
        for (std::size_t i = 2; i <= 2 * k; ++i) exact *= static_cast<double>(i);
        exact *= std::pow(a * b, static_cast<double>(k));
        const std::vector<std::size_t> mult(2, k);
        const double ratio = perm::perm_spm(rows, mult) / exact;
        double binomial = 1.0;
        for (std::size_t i = 1; i <= k; ++i)
            binomial = binomial * static_cast<double>(k + i) / static_cast<double>(i);
        const double kd = static_cast<double>(k);
        const double closed = std::pow(4.0, kd) / (binomial * std::sqrt(M_PI * kd));
        CHECK(ratio == doctest::Approx(closed).epsilon(1e-9));
        CHECK(std::fabs(ratio - (1.0 + 1.0 / (8.0 * kd))) < 0.01 / kd);
        CHECK(ratio < previous);
        previous = ratio;
    }
}

TEST_CASE("the saddle point refuses what it cannot expand") {
    Matrix<double> zero_row(3, 3, 0.0);
    const double v[3][3] = {{1, 2, 3}, {0, 0, 0}, {4, 5, 6}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) zero_row(i, j) = v[i][j];
    CHECK_THROWS(perm::perm_spm(zero_row));

    Matrix<double> negative(2, 2, 1.0);
    negative(0, 1) = -1.0;
    CHECK_THROWS(perm::perm_spm(negative));

    const Matrix<double> wide(4, 2, 1.0);
    CHECK_THROWS(perm::perm_spm(wide));                                 // not square, no m
    CHECK_THROWS(perm::perm_spm(wide, std::vector<std::size_t>(2, 1))); // sum(m) != rows
    CHECK(perm::perm_spm(wide, std::vector<std::size_t>(2, 2)) > 0.0);  // sum(m) == rows
}
