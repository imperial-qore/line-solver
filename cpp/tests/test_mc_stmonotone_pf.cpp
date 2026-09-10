/**
 * The three api/mc entries that had no C++ twin until 2026-08-01:
 * ctmc_stmonotone, ctmc_testpf_kolmogorov and ctmc_pseudostochcomp.
 *
 * The Kolmogorov case also pins the reference defect the port fixed: the JAR
 * took the reverse product on the time-reversed generator, whose stationary
 * ratios telescope to 1, so the test could never return false. A pure
 * three-cycle is the smallest witness.
 */

#include "doctest.h"

#include <cmath>
#include <vector>

#include "line/api/mc/ctmc_pseudostochcomp.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_stmonotone.h"
#include "line/api/mc/ctmc_testpf_kolmogorov.h"
#include "line/util/matrix.h"

using namespace line;

namespace {

/** Birth-death generator on n states with the given up and down rates. */
Matrix<double> birth_death(std::size_t n, double up, double down) {
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = up;
        Q(i + 1, i) = down;
    }
    return mc::ctmc_makeinfgen(Q);
}

}  // namespace

TEST_CASE("dtmc_stmonotone dominates P and is monotone in the row index") {
    // A non-monotone stochastic matrix: row 1's tail sits below row 0's.
    Matrix<double> P(3, 3, 0.0);
    P(0, 0) = 0.1; P(0, 1) = 0.2; P(0, 2) = 0.7;
    P(1, 0) = 0.6; P(1, 1) = 0.3; P(1, 2) = 0.1;
    P(2, 0) = 0.2; P(2, 1) = 0.2; P(2, 2) = 0.6;

    const Matrix<double> Q = mc::dtmc_stmonotone(P);

    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK(Q(i, j) >= -1e-12);
            s += Q(i, j);
        }
        CHECK(s == doctest::Approx(1.0));
    }

    // Every row tail of Q dominates the same tail of P (st-order upper bound)...
    for (std::size_t l = 0; l < 3; ++l)
        for (std::size_t i = 0; i < 3; ++i) {
            double tq = 0.0, tp = 0.0;
            for (std::size_t k = l; k < 3; ++k) { tq += Q(i, k); tp += P(i, k); }
            CHECK(tq >= tp - 1e-12);
        }

    // ... and the tails are nondecreasing in i, which is st-monotonicity.
    for (std::size_t l = 0; l < 3; ++l)
        for (std::size_t i = 1; i < 3; ++i) {
            double a = 0.0, b = 0.0;
            for (std::size_t k = l; k < 3; ++k) { a += Q(i - 1, k); b += Q(i, k); }
            CHECK(b >= a - 1e-12);
        }
}

TEST_CASE("dtmc_stmonotone is the identity on an already monotone chain") {
    Matrix<double> P(3, 3, 0.0);
    P(0, 0) = 0.5; P(0, 1) = 0.3; P(0, 2) = 0.2;
    P(1, 0) = 0.3; P(1, 1) = 0.3; P(1, 2) = 0.4;
    P(2, 0) = 0.1; P(2, 1) = 0.2; P(2, 2) = 0.7;

    const Matrix<double> Q = mc::dtmc_stmonotone(P);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(Q(i, j) == doctest::Approx(P(i, j)));
}

TEST_CASE("ctmc_stmonotone returns a generator bounding the uniformized chain") {
    const Matrix<double> Q = birth_death(4, 2.0, 1.0);
    const Matrix<double> Qub = mc::ctmc_stmonotone(Q);

    for (std::size_t i = 0; i < 4; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 4; ++j) s += Qub(i, j);
        CHECK(std::fabs(s) < 1e-12);
        CHECK(Qub(i, i) <= 1e-12);
    }
    // A birth-death chain is already st-monotone, so the bound is its own
    // uniformization mapped back: same generator up to the rounding.
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(Qub(i, j) == doctest::Approx(Q(i, j) / mc::ctmc_maxabs(Q)).epsilon(1e-9));
}

TEST_CASE("ctmc_testpf_kolmogorov accepts a reversible birth-death chain") {
    CHECK(mc::ctmc_testpf_kolmogorov(birth_death(4, 2.0, 1.0)) == true);
}

TEST_CASE("ctmc_testpf_kolmogorov rejects a one-directional three-cycle") {
    // 0 -> 1 -> 2 -> 0 with no reverse edges: the reverse product is zero while
    // the forward product is one, so the criterion must fail. This is the case
    // the pre-fix reference passed.
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 1) = 1.0;
    Q(1, 2) = 1.0;
    Q(2, 0) = 1.0;
    CHECK(mc::ctmc_testpf_kolmogorov(mc::ctmc_makeinfgen(Q)) == false);
}

TEST_CASE("ctmc_testpf_kolmogorov rejects an imbalanced triangle") {
    // Both directions present, but the two cycle products differ (2 vs 1).
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 1) = 2.0; Q(1, 0) = 1.0;
    Q(1, 2) = 1.0; Q(2, 1) = 1.0;
    Q(2, 0) = 1.0; Q(0, 2) = 1.0;
    CHECK(mc::ctmc_testpf_kolmogorov(mc::ctmc_makeinfgen(Q)) == false);
}

TEST_CASE("ctmc_pseudostochcomp splits the blocks and closes the rows") {
    const Matrix<double> Q = birth_death(4, 2.0, 1.0);
    std::vector<std::size_t> keep;
    keep.push_back(0);
    keep.push_back(1);
    const mc::PseudoStochCompResult<double> r = mc::ctmc_pseudostochcomp(Q, keep);

    CHECK(r.Q11.rows() == 2);
    CHECK(r.Q22.rows() == 2);
    CHECK(r.Q11(0, 1) == doctest::Approx(2.0));
    CHECK(r.Q12(1, 0) == doctest::Approx(2.0));
    CHECK(r.Q21(0, 1) == doctest::Approx(1.0));

    // The rank-one return term restores what Q12 leaks, so S closes its rows.
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += r.S(i, j);
        CHECK(std::fabs(s) < 1e-9);
    }

    // Only state 1 leaks into the complement, so only its row is corrected.
    CHECK(std::fabs(r.Tm(0, 0)) < 1e-12);
    CHECK(std::fabs(r.Tm(0, 1)) < 1e-12);
    CHECK(r.Tm(1, 0) + r.Tm(1, 1) == doctest::Approx(2.0));
}

TEST_CASE("ctmc_pseudostochcomp defaults to the first half of the states") {
    const Matrix<double> Q = birth_death(5, 2.0, 1.0);
    const mc::PseudoStochCompResult<double> r = mc::ctmc_pseudostochcomp(Q);
    CHECK(r.Q11.rows() == 3);  // ceil(5/2)
    CHECK(r.Q22.rows() == 2);
}
