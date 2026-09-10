/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Birman-Kogan asymptotics. The multiprogramming model of Table 3 of Birman and
 * Kogan (Stochastic Models 8(3):543-563, 1992) is fully specified in print --
 * J = 2 device groups of 10 and 50 stations, service times .015 and .030, five
 * jobs per chain, branching .50/.50 and .6/.4 alternating, cpu rate mu_k = M*mu0
 * with M = 50 -- and the paper tabulates four algorithms on it. The rows below
 * pin the algorithms rather than this port: the saddle point column is fixed by
 * Algorithm 1's bottleneck classification, and the two load concealment columns by the
 * fixed point of Algorithm 2.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_bk.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_kt.h"
#include "line/api/pfqn/pfqn_mva.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

const double TH[2] = {1 / 0.015, 1 / 0.030};
const int MG[2] = {10, 50};
const int MPAR = 50;

double branch(int k, int j) { return (k % 2 == 0) ? 0.5 : (j == 0 ? 0.6 : 0.4); }

/// The Table 3 network: K dedicated cpus followed by the two device groups.
Matrix<double> multiprogramming(int K, double mu0) {
    const std::size_t M = static_cast<std::size_t>(K) + MG[0] + MG[1];
    Matrix<double> L(M, K, 0.0);
    for (int k = 0; k < K; ++k) L(k, k) = 1.0 / (MPAR * mu0);
    std::size_t row = K;
    for (int j = 0; j < 2; ++j)
        for (int c = 0; c < MG[j]; ++c) {
            for (int k = 0; k < K; ++k) L(row, k) = branch(k, j) / (TH[j] * MG[j]);
            ++row;
        }
    return L;
}

}  // namespace

TEST_CASE("pfqn_bk reproduces the saddle point column of Birman-Kogan Table 3") {
    // U_k = x_k^0 / mu_0k by Corollary 1, capped at one for a saturated chain.
    // x^0 does not depend on mu_0k, so the mu_0 = 4 row fixes every other row of
    // the column by a pure rescaling.
    const int Ks[4] = {2, 3, 4, 5};
    const double want[4][2] = {{90.9, 95.0}, {83.2, 86.0}, {75.6, 76.9}, {69.7, 70.2}};
    for (int i = 0; i < 4; ++i) {
        const int K = Ks[i];
        Matrix<double> L = multiprogramming(K, 4.0);
        std::vector<double> N(K, 5.0), Z(K, 0.0);
        BkResult<double> sp = pfqn_bk(L, N, Z);
        for (int k = 0; k < 2; ++k) {
            const double U = std::min(sp.X[k] * L(k, k), 1.0) * 100;
            CHECK(U == doctest::Approx(want[i][k]).epsilon(0).epsilon(1e-3));
        }
        CHECK(sp.B.empty());  // no chain saturates at mu_0 = 4
    }
    // At mu_0 = 2 every chain's dedicated station is a bottleneck, so Algorithm 1
    // pins all of them and Corollary 1 returns a utilization of one.
    Matrix<double> Lsat = multiprogramming(2, 2.0);
    std::vector<double> Nsat(2, 5.0), Zsat(2, 0.0);
    BkResult<double> sat = pfqn_bk(Lsat, Nsat, Zsat);
    CHECK(sat.A.empty());
    CHECK(sat.B.size() == 2);
    for (int k = 0; k < 2; ++k) CHECK(sat.X[k] * Lsat(k, k) >= 1.0);
}

TEST_CASE("pfqn_bklc reproduces the two load concealment columns of Birman-Kogan Table 3") {
    struct Row {
        int K;
        double mu0;
        double mva[2];
        double ue[2];
    };
    const Row rows[4] = {{2, 4.0, {70.2, 72.3}, {69.8, 72.0}},
                         {3, 4.0, {67.2, 69.1}, {66.8, 68.7}},
                         {2, 2.0, {93.1, 94.0}, {93.4, 94.4}},
                         {5, 1.5, {95.9, 96.4}, {96.3, 96.8}}};
    for (int i = 0; i < 4; ++i) {
        Matrix<double> L = multiprogramming(rows[i].K, rows[i].mu0);
        std::vector<double> N(rows[i].K, 5.0), Z(rows[i].K, 0.0);
        BkLcResult<double> tm = pfqn_bklc(L, N, Z, "mva", 1e-12, 2000);
        BkLcResult<double> tu = pfqn_bklc(L, N, Z, "ue", 1e-12, 2000);
        for (int k = 0; k < 2; ++k) {
            CHECK(tm.X[k] * L(k, k) * 100 == doctest::Approx(rows[i].mva[k]).epsilon(0).epsilon(2e-3));
            CHECK(tu.X[k] * L(k, k) * 100 == doctest::Approx(rows[i].ue[k]).epsilon(0).epsilon(3e-3));
        }
    }
}

TEST_CASE("pfqn_bk equals the plain saddle point where the paper's structure is absent") {
    // With a think time in every chain there are no dedicated stations to keep
    // out of the exponent, so Proposition 1 IS the multidimensional saddle point
    // that pfqn_kt already computes.
    Matrix<double> L(3, 2, 0.0);
    const double v[3][2] = {{1.0, 0.6}, {0.8, 1.2}, {0.5, 0.9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) L(i, j) = v[i][j];
    std::vector<double> N = {10.0, 8.0}, Z = {2.0, 1.0};
    CHECK(pfqn_bk(L, N, Z).lG == doctest::Approx(pfqn_kt(L, N, Z).lG).epsilon(1e-12));
}

TEST_CASE("pfqn_bkue is exact on one slow station against a group of replicas") {
    // The regime the uniform expansion is built for: a single dominant pole
    // against M >> 1 identical stations.
    Matrix<double> L(13, 1, 0.1);
    L(0, 0) = 0.9;
    std::vector<double> Lv(13, 0.1);
    Lv[0] = 0.9;
    std::vector<int> Ni(1, 20);
    std::vector<double> N(1, 20.0), Z(1, 0.0);
    Matrix<double> Zm(1, 1, 0.0);
    const double exact = pfqn_ca(L, Ni, Zm).lG;
    CHECK(pfqn_bkue(Lv, 20.0, 0.0).lG == doctest::Approx(exact).epsilon(0).epsilon(1e-4));
    CHECK(pfqn_bk(L, N, Z).lG == doctest::Approx(exact).epsilon(0).epsilon(1e-4));
    // and it degenerates to the plain saddle point with no group to stand against
    std::vector<double> Ld = {0.5, 0.4, 0.3};
    Matrix<double> Lm(3, 1, 0.0);
    for (int i = 0; i < 3; ++i) Lm(i, 0) = Ld[i];
    std::vector<double> N60(1, 60.0), Z5(1, 5.0);
    CHECK(pfqn_bkue(Ld, 60.0, 5.0).lG == doctest::Approx(pfqn_kt(Lm, N60, Z5).lG).epsilon(1e-9));
}
