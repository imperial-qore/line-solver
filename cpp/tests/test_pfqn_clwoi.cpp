/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_clwoi: the order-independent member of the CLW inversion family.
 *
 * Expected values come from MATLAB (lineStart; then pfqn_clwoi), quoted to 17
 * significant digits. Like the rest of the family the routine is inherently
 * inexact (it inverts a generating function on a contour), so the MATLAB
 * comparison is at 1e-8 rather than bit-level: MATLAB sums the innermost chain
 * vectorized while this port accumulates it sequentially, which moves the last
 * digits without moving the answer. The independent check is against the exact
 * convolution pfqn_ncoi, at the accuracy the method actually has.
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_clwoi.h"
#include "line/api/pfqn/pfqn_ncoi.h"

using line::Matrix;
using line::num_traits;
using line::pfqn::OiRate;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_clw;
using line::pfqn::pfqn_clwoi;
using line::pfqn::pfqn_ncoi;

namespace {

double relerr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

int busy(const std::vector<int>& n) {
    int c = 0;
    for (std::size_t r = 0; r < n.size(); ++r)
        if (n[r] > 0) c++;
    return c;
}

Matrix<double> mat(std::size_t rows, std::size_t cols, const double* v) {
    Matrix<double> m(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) m(i, j) = v[i * cols + j];
    return m;
}

}  // namespace

TEST_CASE("pfqn_clwoi inverts one OI station") {
    // MATLAB pfqn_clwoi([1 2], [5 4], {@(n) sum(n>0)}).
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 2.0;
    std::vector<int> N(2);
    N[0] = 5;
    N[1] = 4;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return static_cast<double>(busy(n)); });

    const line::pfqn::ClwResult<double> r = pfqn_clwoi(Z, N, mu);
    CHECK(relerr(r.G, 19.016666704936238) < 1e-8);
    CHECK(relerr(r.lG, 2.9453157896523976) < 1e-8);
    // the exact convolution, at the accuracy the inversion has
    CHECK(relerr(r.G, pfqn_ncoi(Z, N, mu).G) < 1e-7);
}

TEST_CASE("pfqn_clwoi inverts two OI stations with visits") {
    // MATLAB pfqn_clwoi([0.5 1.5], [4 3], {@(n) 1+0.5*sum(n>0),
    //        @(n) 2*any(n>0)+0.25*(n(2)>0)}, [1 0.5; 0.8 1.2]).
    std::vector<double> Z(2);
    Z[0] = 0.5;
    Z[1] = 1.5;
    std::vector<int> N(2);
    N[0] = 4;
    N[1] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return 1.0 + 0.5 * busy(n); });
    mu.push_back([](const std::vector<int>& n) {
        return 2.0 * (busy(n) > 0 ? 1.0 : 0.0) + 0.25 * (n[1] > 0 ? 1.0 : 0.0);
    });
    const double vv[4] = {1.0, 0.5, 0.8, 1.2};
    const Matrix<double> vis = mat(2, 2, vv);

    const line::pfqn::ClwResult<double> r = pfqn_clwoi(Z, N, mu, vis);
    CHECK(relerr(r.G, 6.0565705741983091) < 1e-8);
    CHECK(relerr(r.lG, 1.8011437280442113) < 1e-8);
    CHECK(relerr(r.G, pfqn_ncoi(Z, N, mu, vis).G) < 1e-7);
}

TEST_CASE("pfqn_clwoi drops a zero-population chain exactly") {
    // MATLAB pfqn_clwoi([1 2 1], [4 0 3], {@(n) 1+(n(1)>0)+2*(n(3)>0)}, [1 0.5 0.8]).
    // The contour count 2 l_j K_j vanishes for K_j = 0, so the chain must be
    // removed by restricting the generating function to z_2 = 0, not inverted.
    std::vector<double> Z(3);
    Z[0] = 1.0;
    Z[1] = 2.0;
    Z[2] = 1.0;
    std::vector<int> N(3);
    N[0] = 4;
    N[1] = 0;
    N[2] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) {
        return 1.0 + (n[0] > 0 ? 1.0 : 0.0) + 2.0 * (n[2] > 0 ? 1.0 : 0.0);
    });
    const double vv[3] = {1.0, 0.5, 0.8};
    const Matrix<double> vis = mat(1, 3, vv);

    const line::pfqn::ClwResult<double> r = pfqn_clwoi(Z, N, mu, vis);
    CHECK(std::isfinite(r.G));
    CHECK(relerr(r.G, 0.2226763117835889) < 1e-8);
    CHECK(relerr(r.G, pfqn_ncoi(Z, N, mu, vis).G) < 1e-7);
}

TEST_CASE("pfqn_clwoi on load-independent stations is the ordinary constant") {
    // Unit rates make F_i(z) = 1/(1 - sum_r v_{i,r} z_r), the factor pfqn_clw
    // inverts, so the OI form must land on the same constant as pfqn_ca and
    // pfqn_clw on the demands carried by the visit matrix.
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 0.5;
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>&) { return 1.0; });
    mu.push_back([](const std::vector<int>&) { return 1.0; });
    const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
    const Matrix<double> vis = mat(2, 2, Lv);

    const line::pfqn::ClwResult<double> r = pfqn_clwoi(Z, N, mu, vis);
    CHECK(relerr(r.G, 0.11608708335106621) < 1e-8);
    const double Zd[2] = {1.0, 0.5};
    const double ca = pfqn_ca(mat(2, 2, Lv), N, mat(1, 2, Zd)).G;
    CHECK(relerr(r.G, ca) < 1e-9);
    CHECK(relerr(r.G, pfqn_clw(mat(2, 2, Lv), N, Z).G) < 1e-8);
}

TEST_CASE("pfqn_clwoi rejects a rate that is not support-only") {
    // mu(n) = 1 + sum(n) is a balanced-fairness rate, not an OI one: its
    // transform is not the finite rational function the inversion assumes, so
    // the routine must refuse rather than return a wrong constant.
    std::vector<double> Z(2, 1.0);
    std::vector<int> N(2, 3);
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return 1.0 + n[0] + n[1]; });
    CHECK_THROWS_AS(pfqn_clwoi(Z, N, mu), line::InputError);
}

TEST_CASE("pfqn_clwoi rejects a rate that varies at a single interior state") {
    // mu(n) = 1 except at n = 3: the support is the same everywhere, so the rate
    // agrees at the indicator n = 1 and at the fullest n = 4 and mid n = 2 states
    // that the old sampled check probed. Only a scan of the whole lattice sees it.
    std::vector<double> Z(1, 1.0);
    std::vector<int> N(1, 4);
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return n[0] == 3 ? 2.0 : 1.0; });
    CHECK_THROWS_AS(pfqn_clwoi(Z, N, mu), line::InputError);
}

TEST_CASE("pfqn_clwoi returns the trivial constants on empty populations") {
    std::vector<double> Z(2, 1.0);
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return static_cast<double>(busy(n)); });
    std::vector<int> N(2, 0);
    CHECK(pfqn_clwoi(Z, N, mu).G == 1.0);
    N[0] = -1;
    CHECK(pfqn_clwoi(Z, N, mu).G == 0.0);
}
