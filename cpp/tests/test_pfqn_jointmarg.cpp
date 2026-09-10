/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_jointmarg: the joint law of the per-station TOTAL queue lengths, and in
 * particular its permanent ENGINES.
 *
 * The engines carry different guarantees and the tests keep them apart. "exact"
 * is the reference the others are measured against. "spm" is the saddle-point
 * expansion and is the only engine that does NOT expand the demand matrix to
 * order sum(N): it takes the row-replicated matrix with the class populations
 * as column multiplicities, which is the regime the expansion is asymptotic in.
 * That is what these tests pin -- not a golden number, but the LAW the error
 * obeys: it must overestimate, must fall like 1/min(N), and must be nearly the
 * same factor at every state, so that renormalizing a sweep removes most of it.
 * The same three assertions are made in the native Python and JAR suites.
 */
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_jointmarg.h"

namespace pfqn = line::pfqn;
using line::Matrix;

namespace {

/** The 3-station 2-class demand matrix the saddle-point engine is calibrated on. */
Matrix<double> spm_demands() {
    const double v[3][2] = {{0.286, 0.437}, {1.001, 0.782}, {0.294, 0.633}};
    Matrix<double> L(3, 2, 0.0);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) L(i, j) = v[i][j];
    return L;
}

}  // namespace

TEST_CASE("the spm engine error falls as the populations grow") {
    const Matrix<double> L = spm_demands();
    double previous = std::numeric_limits<double>::infinity();
    for (int k = 1; k <= 4; ++k) {
        const std::vector<int> N{k, k};
        double worst = 0.0, mass_exact = 0.0, mass_spm = 0.0;
        for (int i = 0; i <= 2 * k; ++i)
            for (int j = 0; i + j <= 2 * k; ++j) {
                const std::vector<int> n{i, j, 2 * k - i - j};
                const double ex = pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "exact");
                const double sp = pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "spm");
                mass_exact += ex;
                mass_spm += sp;
                if (ex > 1e-12) {
                    // the expansion overestimates the permanent, hence the probability
                    CHECK(sp >= ex * (1.0 - 1e-9));
                    worst = std::max(worst, std::fabs(sp - ex) / ex);
                }
            }
        // G is passed as 1 here, so the two masses are comparable but not 1
        CHECK(mass_spm > mass_exact);
        CHECK(worst < 0.2 / static_cast<double>(k));   // tracks 1/(8 min N), here R = 2
        CHECK(worst < previous);
        previous = worst;
    }
}

TEST_CASE("the spm engine bias cancels under renormalization") {
    const Matrix<double> L = spm_demands();
    double previous = std::numeric_limits<double>::infinity();
    for (int k = 1; k <= 3; ++k) {
        const std::vector<int> N{k, k};
        std::vector<double> exact, spm;
        for (int i = 0; i <= 2 * k; ++i)
            for (int j = 0; i + j <= 2 * k; ++j) {
                const std::vector<int> n{i, j, 2 * k - i - j};
                exact.push_back(pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "exact"));
                spm.push_back(pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "spm"));
            }
        double mass_exact = 0.0, mass_spm = 0.0, raw = 0.0, tvd = 0.0;
        for (std::size_t t = 0; t < spm.size(); ++t) {
            mass_exact += exact[t];
            mass_spm += spm[t];
        }
        for (std::size_t t = 0; t < spm.size(); ++t) {
            if (exact[t] > 1e-12) raw = std::max(raw, std::fabs(spm[t] - exact[t]) / exact[t]);
            tvd += std::fabs(spm[t] / mass_spm - exact[t] / mass_exact);
        }
        tvd *= 0.5;
        CHECK(tvd < 0.1 * raw);
        CHECK(tvd < previous);
        previous = tvd;
    }
}

TEST_CASE("the spm engine refuses a structural zero and an unknown name") {
    Matrix<double> L = spm_demands();
    L(0, 0) = 0.0;
    const std::vector<int> N{2, 1}, n{1, 1, 1};
    // A zero demand that reaches the replicated matrix has no full support.
    CHECK_THROWS(pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "spm"));
    // The exact engine is unaffected on the same state.
    CHECK(pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "exact") > 0.0);
    CHECK_THROWS(pfqn::pfqn_jointmarg<double>(n, spm_demands(), N, {}, 1.0, "nope"));
}

TEST_CASE("a class holding no jobs does not reach the spm full-support check") {
    // Its column is dropped rather than passed with multiplicity zero, so a zero
    // demand there is irrelevant to the permanent and must not be refused.
    Matrix<double> L = spm_demands();
    L(1, 1) = 0.0;
    const std::vector<int> N{3, 0}, n{1, 1, 1};
    const double ex = pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "exact");
    const double sp = pfqn::pfqn_jointmarg<double>(n, L, N, {}, 1.0, "spm");
    CHECK(sp == doctest::Approx(ex).epsilon(1e-12));
}
