/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Achievable-region LP relaxation of a multiclass open Markovian network.
 *
 * Reference values come from the MATLAB reference
 * (matlab/src/api/npfqn/npfqn_bnd_bpt.m) and agree with the JAR and python
 * ports to the digits printed here. M/M/1 is exact, so the first block is a
 * closed-form check rather than a pinned golden.
 *
 * The tolerances below are ABSOLUTE, so they are written as an explicit
 * fabs difference rather than doctest::Approx, whose epsilon is relative.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_bnd_bpt.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using namespace line;

namespace {

double solve(const std::vector<double>& l0, const std::vector<double>& mu, const Matrix<double>& P,
             const std::vector<std::size_t>& st, const std::vector<double>& c) {
    return npfqn::npfqn_bnd_bpt(l0, mu, P, st, c).zlb;
}

}  // namespace

TEST_CASE("npfqn_bnd_bpt is exact on M/M/1") {
    const double rho[4] = {0.3, 0.5, 0.7, 0.9};
    for (int i = 0; i < 4; ++i) {
        CAPTURE(rho[i]);
        Matrix<double> P(1, 1, 0.0);
        const double z = solve({rho[i]}, {1.0}, P, {0}, {1.0});
        CHECK(std::fabs(z - 1.0 / (1.0 - rho[i])) <= 1e-9);
    }
}

TEST_CASE("npfqn_bnd_bpt with two classes at one station") {
    Matrix<double> P(2, 2, 0.0);
    const std::vector<double> l0 = {0.4, 0.4}, mu = {2.0, 1.0};
    const std::vector<std::size_t> st = {0, 0};
    CHECK(std::fabs(solve(l0, mu, P, st, {1.0, 1.0}) - 3.4375) <= 1e-9);
    CHECK(std::fabs(solve(l0, mu, P, st, {1.0, 0.0}) - 0.625) <= 1e-9);
    CHECK(std::fabs(solve(l0, mu, P, st, {0.0, 1.0}) - 1.0 / 0.6) <= 1e-9);
    // the bound must sit below the cmu-optimal (nonpreemptive priority to
    // class 1) sum of mean sojourn times, 3.6875
    CHECK(solve(l0, mu, P, st, {1.0, 1.0}) <= 3.6875 + 1e-12);
}

TEST_CASE("npfqn_bnd_bpt on a tandem: station 1 is tight, station 2 is not") {
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    const std::vector<double> l0 = {0.5, 0.0}, mu = {1.0, 1.0};
    const std::vector<std::size_t> st = {0, 1};
    CHECK(std::fabs(solve(l0, mu, P, st, {1.0, 1.0}) - 3.0) <= 1e-9);
    CHECK(std::fabs(solve(l0, mu, P, st, {1.0, 0.0}) - 2.0) <= 1e-9);  // exact M/M/1
    CHECK(std::fabs(solve(l0, mu, P, st, {0.0, 1.0}) - 1.0) <= 1e-9);  // service time only
}

TEST_CASE("npfqn_bnd_bpt on the Lu-Kumar re-entrant line, 2 stations and 4 classes") {
    Matrix<double> P(4, 4, 0.0);
    P(0, 1) = 1.0;
    P(1, 2) = 1.0;
    P(2, 3) = 1.0;
    const std::vector<double> l0 = {1.0, 0.0, 0.0, 0.0};
    const std::vector<double> mu = {1 / 0.3, 1 / 0.6, 1 / 0.3, 1 / 0.1};
    const std::vector<std::size_t> st = {0, 1, 1, 0};
    CHECK(std::fabs(solve(l0, mu, P, st, {1.0, 1.0, 1.0, 1.0}) - 7.3) <= 1e-9);
    const double xref[4] = {3.0 / 7.0, 0.6, 0.3, 0.1};
    for (int r = 0; r < 4; ++r) {
        CAPTURE(r);
        std::vector<double> e(4, 0.0);
        e[static_cast<std::size_t>(r)] = 1.0;
        CHECK(std::fabs(solve(l0, mu, P, st, e) - xref[r]) <= 1e-9);
    }
}

TEST_CASE("npfqn_bnd_bpt on M/M/1 in exact arithmetic") {
    typedef Rational R;
    Matrix<R> P(1, 1, R(0));
    const std::vector<R> l0(1, R(7) / R(10)), mu(1, R(1)), c(1, R(1));
    const std::vector<std::size_t> st(1, 0);
    const npfqn::BndBpt<R> r = npfqn::npfqn_bnd_bpt(l0, mu, P, st, c);
    CHECK(r.zlb == R(10) / R(3));
}

TEST_CASE("npfqn_bnd_bpt refuses a saturated station rather than answering") {
    Matrix<double> P(1, 1, 0.0);
    CHECK_THROWS_AS(solve({1.2}, {1.0}, P, {0}, {1.0}), std::exception);
}
