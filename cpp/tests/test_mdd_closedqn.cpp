/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `mdd_closedqn`, the MDD-stored exact closed-network solve.
 *
 * THE ORACLE IS THE GORDON-NEWELL CLOSED FORM, not another codebase: on a
 * 2-station cyclic network with mu = {1, 2} and N = 3 the stationary law is
 * pi(n1) proportional to 2^n1, so QLen = {34/15, 11/15}, X = 14/15 at both
 * stations and U = {14/15, 7/15} are exact rationals. The multiserver/IS case
 * checks the population invariant and the cyclic equal-throughput identity,
 * which a wrong min(n_i, c_i) busy-server term breaks.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mdd/mdd_closedqn.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::mdd::MddClosedQnResult;
using line::mdd::mdd_closedqn;

TEST_CASE("mdd_closedqn matches the Gordon-Newell closed form on a 2-station cycle") {
    const std::vector<double> mu{1.0, 2.0};
    Matrix<double> P(2, 2);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;
    const std::vector<double> servers{1.0, 1.0};
    const MddClosedQnResult<double> r = mdd_closedqn(mu, P, servers, 3);

    CHECK(r.states.size() == 4);  // compositions of 3 over 2 stations
    CHECK(r.QLen[0] == doctest::Approx(34.0 / 15.0).epsilon(1e-12));
    CHECK(r.QLen[1] == doctest::Approx(11.0 / 15.0).epsilon(1e-12));
    CHECK(r.X[0] == doctest::Approx(14.0 / 15.0).epsilon(1e-12));
    CHECK(r.X[1] == doctest::Approx(14.0 / 15.0).epsilon(1e-12));
    CHECK(r.U[0] == doctest::Approx(14.0 / 15.0).epsilon(1e-12));
    CHECK(r.U[1] == doctest::Approx(7.0 / 15.0).epsilon(1e-12));
    CHECK(r.stats.num_states == 4);
}

TEST_CASE("mdd_closedqn multiserver and IS stations keep the closed-network identities") {
    const std::vector<double> mu{1.0, 2.0, 3.0};
    Matrix<double> P(3, 3);
    P(0, 1) = 1.0;
    P(1, 2) = 1.0;
    P(2, 0) = 1.0;
    const std::vector<double> servers{2.0, std::numeric_limits<double>::infinity(), 1.0};
    const int N = 4;
    const MddClosedQnResult<double> r = mdd_closedqn(mu, P, servers, N);

    CHECK(r.states.size() == 15);  // C(N+M-1, M-1)
    double totq = 0.0;
    for (std::size_t i = 0; i < 3; ++i) totq += r.QLen[i];
    CHECK(totq == doctest::Approx(static_cast<double>(N)).epsilon(1e-10));
    // cyclic routing: every station sees the same throughput
    CHECK(r.X[1] == doctest::Approx(r.X[0]).epsilon(1e-10));
    CHECK(r.X[2] == doctest::Approx(r.X[0]).epsilon(1e-10));
    // IS utilisation is a mean count, not a ratio
    CHECK(r.U[1] == doctest::Approx(r.QLen[1]).epsilon(1e-12));

    // reusing the diagram must reproduce the solve with no reachability phase
    const MddClosedQnResult<double> r2 = mdd_closedqn(mu, P, servers, N, &r.mdd);
    CHECK(r2.time_reach == 0.0);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(r2.QLen[i] == doctest::Approx(r.QLen[i]).epsilon(1e-14));
}
