/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_pnt / map_pntiter: the counting probabilities of a MAP by uniformization.
 *
 * THE ORACLES ARE IDENTITIES, not another implementation, because the reference
 * fails them (see the header note): any counting law must satisfy
 *   P_0(t)          = exp(D0 t)
 *   sum_n P_n(t)    = exp((D0 + D1) t)
 *   sum_n n pi P_n(t) e = lambda t   with pi the phase-stationary law map_prob
 * and on a Poisson process P_n(t) is the Poisson pmf. MATLAB's map_pntiter
 * misses the first by 1.9e-1, the second by 5.0e-1 and returns 0 for the third
 * on an Erlang-2 MAP of mean 1 at t = 1.3.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_pnt.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

/** exp(A t) by scaling and squaring around a truncated Taylor series. */
Matrix<double> expmt(const Matrix<double>& A, double t) {
    const std::size_t n = A.rows();
    double nrm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double r = 0.0;
        for (std::size_t j = 0; j < n; ++j) r += std::fabs(A(i, j) * t);
        nrm = std::max(nrm, r);
    }
    int s = 0;
    while (nrm > 0.5) { nrm /= 2.0; ++s; }
    const double h = t / std::pow(2.0, s);

    Matrix<double> R(n, n, 0.0), term(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) { R(i, i) = 1.0; term(i, i) = 1.0; }
    for (int k = 1; k <= 40; ++k) {
        Matrix<double> nx(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                double v = 0.0;
                for (std::size_t q = 0; q < n; ++q) v += term(i, q) * A(q, j) * h;
                nx(i, j) = v / static_cast<double>(k);
            }
        term = nx;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) R(i, j) += term(i, j);
    }
    for (int k = 0; k < s; ++k) {
        Matrix<double> sq(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                double v = 0.0;
                for (std::size_t q = 0; q < n; ++q) v += R(i, q) * R(q, j);
                sq(i, j) = v;
            }
        R = sq;
    }
    return R;
}

double maxdiff(const Matrix<double>& A, const Matrix<double>& B) {
    double d = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) d = std::max(d, std::fabs(A(i, j) - B(i, j)));
    return d;
}

}  // namespace

TEST_CASE("P_n(t) of a Poisson process is the Poisson pmf") {
    const double lam = 1.5, t = 2.0;
    mam::Map<double> m;
    m.D0 = Matrix<double>(1, 1, -lam);
    m.D1 = Matrix<double>(1, 1, lam);
    const std::vector<Matrix<double>> P = mam::map_pnt(m, 6, t);
    double logf = 0.0;
    for (std::size_t n = 0; n <= 6; ++n) {
        if (n > 0) logf += std::log(static_cast<double>(n));
        const double want = std::exp(-lam * t + static_cast<double>(n) * std::log(lam * t) - logf);
        CHECK(P[n](0, 0) == doctest::Approx(want).epsilon(1e-10));
    }
}

TEST_CASE("P_0(t) is exp(D0 t) on a genuine two-phase MAP") {
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const double t = 1.3;
    const std::vector<Matrix<double>> P = mam::map_pnt(m, 0, t);
    CHECK(maxdiff(P[0], expmt(m.D0, t)) < 1e-9);
}

TEST_CASE("the counting probabilities sum to exp((D0+D1) t)") {
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const double t = 1.3;
    const std::size_t na = 30;
    const std::vector<Matrix<double>> P = mam::map_pnt(m, na, t);
    Matrix<double> S(m.order(), m.order(), 0.0);
    for (std::size_t n = 0; n <= na; ++n)
        for (std::size_t i = 0; i < m.order(); ++i)
            for (std::size_t j = 0; j < m.order(); ++j) S(i, j) += P[n](i, j);
    CHECK(maxdiff(S, expmt(mam::map_infgen(m), t)) < 1e-9);
}

TEST_CASE("the mean count from the stationary phase is lambda t") {
    // The start vector must be the PHASE-process stationary law map_prob, which
    // solves pi (D0 + D1) = 0. map_pie is the embedded law at arrival epochs and
    // is a different vector; starting from it gives a different, correct, mean.
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const double t = 1.3;
    const std::size_t na = 30;
    const std::vector<Matrix<double>> P = mam::map_pnt(m, na, t);
    const std::vector<double> pie = mam::map_prob(m);
    double mc = 0.0;
    for (std::size_t n = 0; n <= na; ++n) {
        double s = 0.0;
        for (std::size_t i = 0; i < m.order(); ++i)
            for (std::size_t j = 0; j < m.order(); ++j) s += pie[i] * P[n](i, j);
        mc += static_cast<double>(n) * s;
    }
    CHECK(mc == doctest::Approx(mam::map_lambda(m) * t).epsilon(1e-8));
}

TEST_CASE("the identities hold on a correlated MMPP(2) as well") {
    // Two phases with different arrival rates and slow modulation: the counting
    // law is genuinely phase dependent, which is what the reference gets wrong.
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -(3.0 + 0.2); m.D0(0, 1) = 0.2;
    m.D0(1, 0) = 0.1;          m.D0(1, 1) = -(0.5 + 0.1);
    m.D1(0, 0) = 3.0;
    m.D1(1, 1) = 0.5;

    const double t = 2.0;
    const std::size_t na = 40;
    const std::vector<Matrix<double>> P = mam::map_pnt(m, na, t);
    CHECK(maxdiff(P[0], expmt(m.D0, t)) < 1e-9);

    Matrix<double> S(2, 2, 0.0);
    for (std::size_t n = 0; n <= na; ++n)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) S(i, j) += P[n](i, j);
    CHECK(maxdiff(S, expmt(mam::map_infgen(m), t)) < 1e-9);

    const std::vector<double> pie = mam::map_prob(m);
    double mc = 0.0;
    for (std::size_t n = 0; n <= na; ++n) {
        double s = 0.0;
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) s += pie[i] * P[n](i, j);
        mc += static_cast<double>(n) * s;
    }
    CHECK(mc == doctest::Approx(mam::map_lambda(m) * t).epsilon(1e-7));
}

TEST_CASE("map_pntiter returns the highest count and honours an explicit M") {
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const double t = 1.3;
    const Matrix<double> a = mam::map_pntiter(m, 3, t);
    const Matrix<double> b = mam::map_pntiter(m, 3, t, 6);
    const Matrix<double> c = mam::map_pntiter(m, 3, t, 0);  // no squaring at all
    CHECK(maxdiff(a, b) < 1e-9);
    CHECK(maxdiff(a, c) < 1e-9);
}
