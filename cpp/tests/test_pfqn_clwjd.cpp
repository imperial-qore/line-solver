/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_clwjd: the limited joint-dependent member of the CLW inversion family,
 * plus the pfqn_ncjd / pfqn_mvajd aliases.
 *
 * Expected values come from MATLAB (lineStart; then pfqn_clwjd), quoted to 17
 * significant digits, and are checked at 1e-8: like the rest of the family the
 * routine inverts a generating function on a contour, and MATLAB sums the
 * innermost chain vectorized where this port accumulates it sequentially. The
 * independent check is against the exact convolution pfqn_ncjd.
 *
 * The structural checks are the two reductions: cutoff 1 must reproduce
 * pfqn_clwoi on a support-only rate, and a rate saturating at the server count
 * must reproduce pfqn_clw_lld on the load-dependent encoding of the same model.
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_clwjd.h"
#include "line/api/pfqn/pfqn_clwoi.h"
#include "line/api/pfqn/pfqn_mvajd.h"
#include "line/api/pfqn/pfqn_ncjd.h"

using line::Matrix;
using line::pfqn::OiRate;
using line::pfqn::pfqn_clw_lld;
using line::pfqn::pfqn_clwjd;
using line::pfqn::pfqn_clwoi;
using line::pfqn::pfqn_mvaoi;
using line::pfqn::pfqn_mvajd;
using line::pfqn::pfqn_ncjd;

namespace {

double relerr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

int busy(const std::vector<int>& n) {
    int c = 0;
    for (std::size_t r = 0; r < n.size(); ++r)
        if (n[r] > 0) c++;
    return c;
}

double clipsum(const std::vector<int>& n, int c) {
    double s = 0;
    for (std::size_t r = 0; r < n.size(); ++r) s += std::min(n[r], c);
    return s;
}

Matrix<double> mat(std::size_t rows, std::size_t cols, const double* v) {
    Matrix<double> m(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) m(i, j) = v[i * cols + j];
    return m;
}

Matrix<int> imat(std::size_t rows, std::size_t cols, const int* v) {
    Matrix<int> m(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) m(i, j) = v[i * cols + j];
    return m;
}

}  // namespace

TEST_CASE("pfqn_clwjd at cutoff 1 is pfqn_clwoi") {
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 2.0;
    std::vector<int> N(2);
    N[0] = 5;
    N[1] = 4;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return static_cast<double>(busy(n)); });
    const int lv[2] = {1, 1};

    const line::pfqn::ClwResult<double> jd =
        pfqn_clwjd(Z, N, mu, Matrix<double>(), imat(1, 2, lv));
    const line::pfqn::ClwResult<double> oi = pfqn_clwoi(Z, N, mu);
    // A support-only rate clipped at 1 IS the support indicator, so the two
    // recursions traverse the same regions in the same order: bit equality, not
    // mere closeness, is the right assertion.
    CHECK(jd.G == oi.G);
    CHECK(relerr(jd.G, pfqn_ncjd(Z, N, mu).G) < 1e-7);
}

TEST_CASE("pfqn_clwjd inverts a rate saturating at two jobs per class") {
    // MATLAB pfqn_clwjd([1 2], [5 4], {@(n) 1+sum(min(n,2))}, [], [2 2]).
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 2.0;
    std::vector<int> N(2);
    N[0] = 5;
    N[1] = 4;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return 1.0 + clipsum(n, 2); });
    const int lv[2] = {2, 2};

    const line::pfqn::ClwResult<double> r =
        pfqn_clwjd(Z, N, mu, Matrix<double>(), imat(1, 2, lv));
    CHECK(relerr(r.G, 0.28304853784201872) < 1e-8);
    CHECK(relerr(r.G, pfqn_ncjd(Z, N, mu).G) < 1e-7);

    SUBCASE("the all-N cutoff is exact too, just slower") {
        const line::pfqn::ClwResult<double> full = pfqn_clwjd(Z, N, mu);
        CHECK(relerr(full.G, pfqn_ncjd(Z, N, mu).G) < 1e-7);
    }
}

TEST_CASE("pfqn_clwjd reproduces pfqn_clw_lld on a multiserver station") {
    // min(sum n, c) is a function of the clipped vector once every cutoff is c,
    // so the joint-dependent encoding must land on the load-dependent constant.
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 0.5;
    std::vector<int> N(2);
    N[0] = 4;
    N[1] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) {
        return static_cast<double>(std::min(n[0] + n[1], 3));
    });
    const double Dv[2] = {0.1, 0.2};
    const int lv[2] = {3, 3};

    const line::pfqn::ClwResult<double> r = pfqn_clwjd(Z, N, mu, mat(1, 2, Dv), imat(1, 2, lv));
    const double muv[7] = {1, 2, 3, 3, 3, 3, 3};
    const line::pfqn::ClwResult<double> lld =
        pfqn_clw_lld(mat(1, 2, Dv), N, Z, mat(1, 7, muv));
    CHECK(relerr(r.G, lld.G) < 1e-8);
    CHECK(relerr(r.G, pfqn_ncjd(Z, N, mu, mat(1, 2, Dv)).G) < 1e-7);
}

TEST_CASE("pfqn_clwjd with two stations and asymmetric cutoffs") {
    // MATLAB pfqn_clwjd([0.5 1.5], [4 3], {...}, [1 0.5; 0.8 1.2], [2 1; 1 1]).
    std::vector<double> Z(2);
    Z[0] = 0.5;
    Z[1] = 1.5;
    std::vector<int> N(2);
    N[0] = 4;
    N[1] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) {
        return 1.0 + 0.5 * std::min(n[0], 2) + 0.25 * std::min(n[1], 1);
    });
    mu.push_back([](const std::vector<int>& n) { return 2.0 + std::min(n[0], 1); });
    const double vv[4] = {1.0, 0.5, 0.8, 1.2};
    const int lv[4] = {2, 1, 1, 1};

    const line::pfqn::ClwResult<double> r = pfqn_clwjd(Z, N, mu, mat(2, 2, vv), imat(2, 2, lv));
    CHECK(relerr(r.G, 3.0279157447846402) < 1e-8);
    CHECK(relerr(r.G, pfqn_ncjd(Z, N, mu, mat(2, 2, vv)).G) < 1e-7);
}

TEST_CASE("pfqn_clwjd rejects a rate that does not saturate at the cutoff") {
    std::vector<double> Z(2, 1.0);
    std::vector<int> N(2, 3);
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return 1.0 + clipsum(n, 3); });
    const int lv[2] = {1, 1};
    CHECK_THROWS_AS(pfqn_clwjd(Z, N, mu, Matrix<double>(), imat(1, 2, lv)), line::InputError);
}

TEST_CASE("pfqn_ncjd and pfqn_mvajd are the joint-dependent names of the OI routines") {
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 2.0;
    std::vector<int> N(2);
    N[0] = 4;
    N[1] = 3;
    std::vector<OiRate<double>> mu;
    mu.push_back([](const std::vector<int>& n) { return 1.0 + clipsum(n, 2); });

    CHECK(pfqn_ncjd(Z, N, mu).G == line::pfqn::pfqn_ncoi(Z, N, mu).G);
    const line::pfqn::MvaoiResult<double> a = pfqn_mvaoi(Z, N, mu, Matrix<double>());
    const line::pfqn::MvaoiResult<double> b = pfqn_mvajd(Z, N, mu, Matrix<double>());
    for (std::size_t r = 0; r < a.X.size(); ++r) CHECK(a.X[r] == b.X[r]);
    // and the mean-value throughput is the normalizing-constant ratio
    std::vector<int> Nm = N;
    Nm[0] -= 1;
    CHECK(relerr(a.X[0], pfqn_ncjd(Z, Nm, mu).G / pfqn_ncjd(Z, N, mu).G) < 1e-9);
}
