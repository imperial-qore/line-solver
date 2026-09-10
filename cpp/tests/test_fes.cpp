/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Flow-equivalent server aggregation. Oracles, in order of strength:
 *   1. Norton's theorem (Chandy, Herzog and Woo 1975): for a single-class
 *      product-form network, replacing a subnetwork by a load-dependent
 *      station whose rates are the isolated subnetwork throughputs leaves the
 *      normalizing constant UNCHANGED. In exact arithmetic that is an equality
 *      of rationals, not an approximation, so the table produced here is
 *      checked against pfqn_ca on the unaggregated model term for term.
 *   2. Closed forms of the throughput table on one- and two-station
 *      subnetworks, computed by hand.
 *   3. The definition of beta, beta_r(n) = X_r(n) |n| / n_r.
 * Every case runs at double and at Rational and the two must agree.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/fes/fes_beta_handle.h"
#include "line/api/fes/fes_compute_throughputs.h"
#include "line/api/fes/ljd_linearize.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_mva.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::fes::fes_beta_handle;
using line::fes::fes_compute_throughputs;
using line::fes::ljd_linearize;

namespace {

constexpr double TOL = 1e-9;

}  // namespace

// ---------------------------------------------------------------------------
// ljd_linearize
// ---------------------------------------------------------------------------

TEST_CASE("ljd_linearize is the 1-based mixed-radix index with clamping") {
    const std::vector<int> cutoffs{2, 3};
    CHECK(ljd_linearize({0, 0}, cutoffs) == 1);
    CHECK(ljd_linearize({1, 0}, cutoffs) == 2);
    CHECK(ljd_linearize({2, 0}, cutoffs) == 3);
    CHECK(ljd_linearize({0, 1}, cutoffs) == 4);   // 1 + 0 + 1*(2+1)
    CHECK(ljd_linearize({2, 3}, cutoffs) == 12);  // 1 + 2 + 3*3
    // beyond the cutoff the index saturates
    CHECK(ljd_linearize({5, 9}, cutoffs) == ljd_linearize({2, 3}, cutoffs));
    CHECK_THROWS_AS(ljd_linearize({1}, cutoffs), line::InputError);
}

// ---------------------------------------------------------------------------
// fes_compute_throughputs
// ---------------------------------------------------------------------------

TEST_CASE("fes_compute_throughputs matches the single-station closed form") {
    // One queue with demand 1/2 and no delay: X(n) = 2 for every n >= 1.
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(1, 2);
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1}, std::vector<bool>{false},
                                             std::vector<int>{3});
    REQUIRE(tbl.size() == 1);
    REQUIRE(tbl[0].size() == 4);
    CHECK(tbl[0][0] == Rational(0));
    for (std::size_t n = 1; n <= 3; ++n) CHECK(tbl[0][n] == Rational(2));
}

TEST_CASE("fes_compute_throughputs honours the server count of a multiserver station") {
    // REGRESSION (2026-08-25). `mi` was handed straight to pfqn_mva, whose mi is
    // NOT a server count: it enters only as the additive term of
    // C(i,s)=L(i,s)*(mi[i]+Qarv), so mi=c INFLATED the residence time by c --
    // making the station SLOWER where c servers should make it faster. One
    // 3-server station of unit demand holding one job answered X=1/3 against the
    // true 1. The kernel goes through pfqn_mvams now. See
    // _kb/03-api-layer.md ("pfqn_mva's mi is NOT a server count").
    const int c = 3;
    Matrix<Rational> L(1, 2);
    L(0, 0) = Rational(1);
    L(0, 1) = Rational(1);
    const std::vector<int> cutoffs{3, 3};
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{c}, std::vector<bool>{false},
                                             cutoffs);
    // Equal exponential rates at a c-server station: the total completion rate is
    // min(|n|,c) and class r takes the share n_r/|n|.
    for (int n1 = 0; n1 <= 3; ++n1)
        for (int n2 = 0; n2 <= 3; ++n2) {
            const int tot = n1 + n2;
            if (tot == 0) continue;
            const std::size_t idx = ljd_linearize(std::vector<int>{n1, n2}, cutoffs) - 1;
            const Rational tref(std::min(tot, c), tot);
            CHECK(tbl[0][idx] == tref * Rational(n1));
            CHECK(tbl[1][idx] == tref * Rational(n2));
        }
    // The signature of the defect was insensitivity to c, so pin that the
    // single-server table DIFFERS.
    const auto single = fes_compute_throughputs(L, std::vector<int>{1}, std::vector<bool>{false},
                                                cutoffs);
    const std::size_t two_two = ljd_linearize(std::vector<int>{2, 2}, cutoffs) - 1;
    CHECK(tbl[0][two_two] != single[0][two_two]);
}

TEST_CASE("fes_compute_throughputs matches a two-station closed form") {
    // Demands 1/3 and 1/5. G(n) = sum_{a+b=n} (1/3)^a (1/5)^b and
    // X(n) = G(n-1)/G(n), which is what MVA returns.
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(1, 3);
    L(1, 0) = Rational(1, 5);
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1, 1},
                                             std::vector<bool>{false, false}, std::vector<int>{3});
    // X(1) = 1/(1/3 + 1/5) = 15/8
    CHECK(tbl[0][1] == Rational(15, 8));
    // G(1) = 8/15, G(2) = 1/9 + 1/15 + 1/25 = 49/225, X(2) = G(1)/G(2)
    CHECK(tbl[0][2] == Rational(8, 15) / Rational(49, 225));
}

TEST_CASE("fes_compute_throughputs handles a pure-delay subnetwork") {
    // Only delay stations: X_k = n_k / Z_k with Z from the delay demands.
    Matrix<Rational> L(1, 2);
    L(0, 0) = Rational(1, 2);
    L(0, 1) = Rational(1, 4);
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1}, std::vector<bool>{true},
                                             std::vector<int>{2, 1});
    // index of n = (2,1) is 1 + 2 + 1*3 = 6
    const std::size_t idx = ljd_linearize({2, 1}, {2, 1});
    CHECK(tbl[0][idx - 1] == Rational(4));  // 2 / (1/2)
    CHECK(tbl[1][idx - 1] == Rational(4));  // 1 / (1/4)
    // a class with no jobs contributes no throughput
    const std::size_t idx0 = ljd_linearize({2, 0}, {2, 1});
    CHECK(tbl[1][idx0 - 1] == Rational(0));
}

TEST_CASE("fes_compute_throughputs reproduces Norton's theorem exactly") {
    // Full single-class network of three queues with demands 1/2, 1/3, 1/5.
    const int N = 4;
    Matrix<Rational> Lfull(3, 1);
    Lfull(0, 0) = Rational(1, 2);
    Lfull(1, 0) = Rational(1, 3);
    Lfull(2, 0) = Rational(1, 5);
    const auto full = line::pfqn::pfqn_ca(Lfull, std::vector<int>{N});

    // Aggregate stations 2 and 3 into a flow-equivalent server whose
    // load-dependent rates are the isolated subnetwork throughputs X_sub(n).
    Matrix<Rational> Lsub(2, 1);
    Lsub(0, 0) = Rational(1, 3);
    Lsub(1, 0) = Rational(1, 5);
    const auto tbl = fes_compute_throughputs(Lsub, std::vector<int>{1, 1},
                                             std::vector<bool>{false, false}, std::vector<int>{N});

    // Composite: station 1 unchanged, plus a rate-1-demand load-dependent
    // station serving at mu(k) = X_sub(k).
    Matrix<Rational> Lc(2, 1);
    Lc(0, 0) = Rational(1, 2);
    Lc(1, 0) = Rational(1);
    Matrix<Rational> mu(2, static_cast<std::size_t>(N), Rational(1));
    for (int k = 1; k <= N; ++k) mu(1, static_cast<std::size_t>(k - 1)) = tbl[0][static_cast<std::size_t>(k)];
    const auto composite = line::pfqn::pfqn_gld(Lc, std::vector<int>{N}, mu);

    // Norton: the two normalizing constants are the same rational number.
    CHECK(composite.G == full.G);

    // and therefore so is the throughput X(N) = G(N-1)/G(N)
    const auto full1 = line::pfqn::pfqn_ca(Lfull, std::vector<int>{N - 1});
    Matrix<Rational> mu1(2, static_cast<std::size_t>(N), Rational(1));
    for (int k = 1; k <= N; ++k) mu1(1, static_cast<std::size_t>(k - 1)) = tbl[0][static_cast<std::size_t>(k)];
    const auto composite1 = line::pfqn::pfqn_gld(Lc, std::vector<int>{N - 1}, mu1);
    CHECK(full1.G / full.G == composite1.G / composite.G);
}

TEST_CASE("fes_compute_throughputs multiclass table equals pfqn_mva state by state") {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5;
    L(0, 1) = 0.3;
    L(1, 0) = 0.4;
    L(1, 1) = 0.6;
    const std::vector<int> cutoffs{2, 2};
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1, 1},
                                             std::vector<bool>{false, false}, cutoffs);
    for (int n0 = 0; n0 <= 2; ++n0)
        for (int n1 = 0; n1 <= 2; ++n1) {
            const std::vector<int> n{n0, n1};
            const std::size_t idx = ljd_linearize(n, cutoffs);
            if (n0 + n1 == 0) {
                CHECK(tbl[0][idx - 1] == 0.0);
                CHECK(tbl[1][idx - 1] == 0.0);
                continue;
            }
            const auto r = line::pfqn::pfqn_mva(L, n);
            for (std::size_t k = 0; k < 2; ++k)
                CHECK(tbl[k][idx - 1] == doctest::Approx(n[k] > 0 ? r.XN[k] : 0.0).epsilon(TOL));
        }
}

TEST_CASE("fes_compute_throughputs exact equals double to rounding") {
    Matrix<double> Ld(2, 1);
    Ld(0, 0) = 1.0 / 3.0;
    Ld(1, 0) = 0.2;
    Matrix<Rational> Lq(2, 1);
    Lq(0, 0) = Rational(1, 3);
    Lq(1, 0) = Rational(1, 5);
    const auto d = fes_compute_throughputs(Ld, std::vector<int>{1, 1},
                                           std::vector<bool>{false, false}, std::vector<int>{5});
    const auto q = fes_compute_throughputs(Lq, std::vector<int>{1, 1},
                                           std::vector<bool>{false, false}, std::vector<int>{5});
    for (std::size_t n = 0; n < d[0].size(); ++n)
        CHECK(static_cast<double>(q[0][n]) == doctest::Approx(d[0][n]).epsilon(TOL));
}

TEST_CASE("fes_compute_throughputs rejects inconsistent inputs") {
    Matrix<double> L(2, 1, 1.0);
    CHECK_THROWS_AS(fes_compute_throughputs(L, std::vector<int>{1, 1}, std::vector<bool>{false},
                                            std::vector<int>{2}),
                    line::InputError);
    CHECK_THROWS_AS(fes_compute_throughputs(L, std::vector<int>{1, 1},
                                            std::vector<bool>{false, false}, std::vector<int>{2, 2}),
                    line::InputError);
}

// ---------------------------------------------------------------------------
// fes_beta_handle
// ---------------------------------------------------------------------------

TEST_CASE("fes_beta_handle returns X_r(n) |n| / n_r") {
    // Single class, table X(n) = 2 for n >= 1: beta(n) = 2 n / n = 2.
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(1, 2);
    const std::vector<int> cutoffs{3};
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1}, std::vector<bool>{false}, cutoffs);
    const auto beta = fes_beta_handle(tbl, cutoffs);
    for (int n = 1; n <= 3; ++n) {
        const std::vector<Rational> v = beta({n});
        REQUIRE(v.size() == 1);
        CHECK(v[0] == Rational(2));
    }
    // n_r = 0 is never consulted by the recurrence and comes back as 1
    CHECK(beta({0})[0] == Rational(1));
    // beyond the cutoff the scaling saturates on the clamped table entry, but
    // the |n| / n_r factor still uses the unclamped population
    CHECK(beta({6})[0] == Rational(2));
}

TEST_CASE("fes_beta_handle multiclass matches the definition entry by entry") {
    Matrix<Rational> L(2, 2);
    L(0, 0) = Rational(1, 2);
    L(0, 1) = Rational(3, 10);
    L(1, 0) = Rational(2, 5);
    L(1, 1) = Rational(3, 5);
    const std::vector<int> cutoffs{2, 2};
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1, 1},
                                             std::vector<bool>{false, false}, cutoffs);
    const auto beta = fes_beta_handle(tbl, cutoffs);
    for (int n0 = 0; n0 <= 2; ++n0)
        for (int n1 = 0; n1 <= 2; ++n1) {
            const std::vector<int> n{n0, n1};
            const std::size_t idx = ljd_linearize(n, cutoffs);
            const std::vector<Rational> v = beta(n);
            for (std::size_t r = 0; r < 2; ++r) {
                if (n[r] == 0) {
                    CHECK(v[r] == Rational(1));
                } else {
                    CHECK(v[r] == tbl[r][idx - 1] * Rational(n0 + n1) / Rational(n[r]));
                }
            }
        }
}

TEST_CASE("fes_beta_handle pads and truncates the population vector") {
    Matrix<Rational> L(1, 2);
    L(0, 0) = Rational(1, 2);
    L(0, 1) = Rational(1, 2);
    const std::vector<int> cutoffs{2, 2};
    const auto tbl = fes_compute_throughputs(L, std::vector<int>{1}, std::vector<bool>{false}, cutoffs);
    const auto beta = fes_beta_handle(tbl, cutoffs);
    // a short vector is zero padded, so beta({1}) == beta({1,0})
    CHECK(beta({1})[0] == beta({1, 0})[0]);
    // a long vector is truncated
    CHECK(beta({1, 1, 7})[0] == beta({1, 1})[0]);
}

TEST_CASE("fes_beta_handle exact equals double to rounding") {
    Matrix<double> Ld(1, 1, 1.0 / 3.0);
    Matrix<Rational> Lq(1, 1, Rational(1, 3));
    const std::vector<int> cutoffs{4};
    const auto bd = fes_beta_handle(fes_compute_throughputs(Ld, std::vector<int>{1},
                                                            std::vector<bool>{false}, cutoffs),
                                    cutoffs);
    const auto bq = fes_beta_handle(fes_compute_throughputs(Lq, std::vector<int>{1},
                                                            std::vector<bool>{false}, cutoffs),
                                    cutoffs);
    for (int n = 1; n <= 4; ++n)
        CHECK(static_cast<double>(bq({n})[0]) == doctest::Approx(bd({n})[0]).epsilon(TOL));
}
