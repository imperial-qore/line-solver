/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Boundaries of the SPM cache normalizing constant.
 *
 * cache_spm returns Ehat(m) = E(m) * prod_l m_l!, the same constant cache_erec
 * evaluates exactly, so cache_erec is the oracle throughout. Two boundaries
 * used to be mishandled in every codebase:
 *
 *   m_l == 0     list l has xi(l)=0, a boundary of the Laplace integral rather
 *                than a direction of it. Left in, the -(1/2) sum_l log xi(l)
 *                prefactor gained ~+17 per empty list: n=12 items at
 *                gamma=(.8,.6,.4) and m=(0,3,3) returned lZ=24.814 against an
 *                exact 9.127, a factor of 6.5e8.
 *   n == sum(m)  every item is cached, the multipliers diverge and
 *                cache_xi_iter cannot converge. It used to be run anyway --
 *                MATLAB and the JAR hung, this port raised NumericError from
 *                the sweep cap and lost the exact Z with it.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_prob_spm.h"
#include "line/api/cache/cache_spm.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::cache::cache_erec;
using line::cache::cache_spm;

namespace {

/** Uniform popularity: 12 items, three lists. */
Matrix<double> uniform_gamma() {
    Matrix<double> g(12, 3, 0.0);
    for (std::size_t i = 0; i < 12; ++i) {
        g(i, 0) = 0.8;
        g(i, 1) = 0.6;
        g(i, 2) = 0.4;
    }
    return g;
}

/** Non-uniform popularity, monotone across the list index as SPM assumes. */
Matrix<double> skewed_gamma() {
    Matrix<double> g(12, 3, 0.0);
    for (std::size_t i = 0; i < 12; ++i) {
        const double gi = 0.3 + 2.7 * ((i + 1.0) / 12.0);
        g(i, 0) = gi;
        g(i, 1) = 0.6 * gi;
        g(i, 2) = 0.3 * gi;
    }
    return g;
}

/** gamma restricted to the given columns, in order. */
Matrix<double> cols(const Matrix<double>& g, const std::vector<std::size_t>& keep) {
    Matrix<double> out(g.rows(), keep.size(), 0.0);
    for (std::size_t i = 0; i < g.rows(); ++i)
        for (std::size_t a = 0; a < keep.size(); ++a) out(i, a) = g(i, keep[a]);
    return out;
}

}  // namespace

TEST_CASE("a zero-capacity list leaves the SPM expansion") {
    const std::vector<Matrix<double>> profiles{uniform_gamma(), skewed_gamma()};
    const std::vector<std::vector<int>> caps{{0, 3, 3}, {3, 0, 3}, {0, 0, 6}, {4, 4, 0}};
    const std::vector<std::vector<std::size_t>> keeps{{1, 2}, {0, 2}, {2}, {0, 1}};

    for (const Matrix<double>& g : profiles) {
        for (std::size_t c = 0; c < caps.size(); ++c) {
            std::vector<int> mk;
            for (std::size_t l : keeps[c]) mk.push_back(caps[c][l]);

            const auto with_empty = cache_spm(g, caps[c]);
            const auto without = cache_spm(cols(g, keeps[c]), mk);
            // dropping an empty list is exact, so the two agree bit for bit
            CHECK(with_empty.lZ == without.lZ);
            // and SPM's own error, not a 1e8 blow-up, separates it from exact
            CHECK(std::abs(with_empty.lZ - std::log(cache_erec(g, caps[c]))) < 0.30);
            // the dropped lists carry xi = 0, the root of their capacity equation
            for (std::size_t l = 0; l < 3; ++l) CHECK((with_empty.xi[l] == 0.0) == (caps[c][l] == 0));
        }
    }
}

TEST_CASE("an empty cache is the unit constant") {
    const auto r = cache_spm(uniform_gamma(), std::vector<int>{0, 0, 0});
    CHECK(r.Z == 1.0);
    CHECK(r.lZ == 0.0);
    for (std::size_t l = 0; l < 3; ++l) CHECK(r.xi[l] == 0.0);
}

TEST_CASE("a full cache returns the exact constant and terminates") {
    const std::vector<Matrix<double>> profiles{uniform_gamma(), skewed_gamma()};
    const std::vector<std::vector<int>> caps{{4, 4, 4}, {12, 0, 0}, {6, 5, 1}};

    for (const Matrix<double>& g : profiles) {
        for (const std::vector<int>& m : caps) {
            const auto r = cache_spm(g, m);
            CHECK(r.Z == doctest::Approx(cache_erec(g, m)).epsilon(1e-12));
            CHECK(r.lZ == doctest::Approx(std::log(cache_erec(g, m))).epsilon(1e-12));
            // the multipliers' limit, not a number
            for (std::size_t l = 0; l < 3; ++l) CHECK(std::isinf(r.xi[l]));
        }
    }
}

TEST_CASE("the interior saddle is unchanged by the boundary handling") {
    const Matrix<double> g = uniform_gamma();
    const std::vector<std::vector<int>> ms{{3, 3, 3}, {4, 3, 2}, {2, 2, 2}, {6, 4, 1}};
    const std::vector<double> want{13.34844416820355, 14.04836686545269, 10.23839884632498,
                                   15.878606855246987};
    for (std::size_t c = 0; c < ms.size(); ++c)
        CHECK(cache_spm(g, ms[c]).lZ == doctest::Approx(want[c]).epsilon(1e-12));
}

TEST_CASE("cache_prob_spm stays a distribution when a list holds one item") {
    // cache_prob_spm evaluates E at oner(m,l), so m_0 = 1 drives cache_spm to 0
    const std::vector<int> m{1, 3, 2};
    for (const Matrix<double>& g : {uniform_gamma(), skewed_gamma()}) {
        const Matrix<double> prob = line::cache::cache_prob_spm(g, m);
        for (std::size_t i = 0; i < prob.rows(); ++i) {
            double rowsum = 0.0;
            for (std::size_t j = 0; j < prob.cols(); ++j) {
                CHECK(prob(i, j) >= 0.0);
                CHECK(prob(i, j) <= 1.0);
                rowsum += prob(i, j);
            }
            CHECK(rowsum == doctest::Approx(1.0).epsilon(1e-9));
        }
        // sum_i P(item i in list l) = m_l holds exactly for E, closely for SPM
        for (std::size_t l = 0; l < 3; ++l) {
            double occ = 0.0;
            for (std::size_t i = 0; i < prob.rows(); ++i) occ += prob(i, l + 1);
            CHECK(occ == doctest::Approx(static_cast<double>(m[l])).epsilon(0.20));
        }
    }
}
