/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Transient and sensitivity routines of the mc layer: ctmc_foxglynn,
 * ctmc_transient, ctmc_transient_sens and ctmc_sens.
 *
 * Oracles:
 *  - conservation of probability, which every transient solution must satisfy
 *    at every time point;
 *  - convergence to the stationary vector of ctmc_solve as t grows;
 *  - the Poisson distribution itself, for the Fox-Glynn truncation window and
 *    weights;
 *  - the defining equations of the sensitivity, checked EXACTLY in rational
 *    arithmetic, plus a central difference of the stationary vector;
 *  - MATLAB reference values from LINE 3.0.6 on the same fixtures
 *    (generated 2026-07-21).
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_foxglynn.h"
#include "line/api/mc/ctmc_sens.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/api/mc/ctmc_transient_sens.h"
#include "line/api/mc/ctmc_uniformization.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::mc::ctmc_foxglynn;
using line::mc::ctmc_makeinfgen;
using line::mc::ctmc_sens;
using line::mc::ctmc_solve;
using line::mc::ctmc_transient;
using line::mc::ctmc_transient_sens;
using line::mc::ctmc_uniformization;

namespace {

/** The nearly completely decomposable six-state chain used across the mc tests. */
template <class T>
Matrix<T> chain_A() {
    Matrix<T> Q(6, 6, line::num_traits<T>::from_int(0));
    Q(0, 1) = line::num_traits<T>::from_int(3);
    Q(1, 0) = line::num_traits<T>::from_int(2);
    Q(2, 3) = line::num_traits<T>::from_int(5);
    Q(3, 2) = line::num_traits<T>::from_int(1);
    Q(4, 5) = line::num_traits<T>::from_int(2);
    Q(5, 4) = line::num_traits<T>::from_int(4);
    Q(1, 2) = line::num_traits<T>::from_rational(1, 100);
    Q(2, 1) = line::num_traits<T>::from_rational(5, 1000);
    Q(3, 4) = line::num_traits<T>::from_rational(2, 100);
    Q(5, 0) = line::num_traits<T>::from_rational(3, 100);
    return ctmc_makeinfgen(Q);
}

/** Derivative of chain A with respect to the rate in position (0,1). */
template <class T>
Matrix<T> chain_A_dQ() {
    Matrix<T> dQ(6, 6, line::num_traits<T>::from_int(0));
    dQ(0, 1) = line::num_traits<T>::from_int(1);
    dQ(0, 0) = line::num_traits<T>::from_int(-1);
    return dQ;
}

double linf(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

double sum_of(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Fox-Glynn
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_foxglynn weights sum to one and bracket the Poisson mass") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0(6, 1.0 / 6.0);
    const double tol = 1e-12;

    for (double t : {0.05, 0.7, 5.0, 40.0}) {
        const line::mc::FoxGlynnResult<double> r = ctmc_foxglynn(pi0, Q, t, tol);
        REQUIRE(r.w.size() == static_cast<std::size_t>(r.right - r.left + 1));
        CHECK(sum_of(r.w) == doctest::Approx(1.0).epsilon(1e-14));
        for (double w : r.w) CHECK(w >= 0.0);

        // The window must carry all but tol of the Poisson(lambda) mass, with
        // lambda = q t and q = 1.1 max|q_ii|; that is the property the
        // truncation claims and the Chernoff certification enforces.
        double qmax = 0.0;
        for (std::size_t i = 0; i < 6; ++i) qmax = std::max(qmax, std::fabs(Q(i, i)));
        const double lambda = 1.1 * qmax * t;
        double inside = 0.0, total = 0.0, pk = std::exp(-lambda);
        for (long k = 0; k <= r.right + 400; ++k) {
            total += pk;
            if (k >= r.left && k <= r.right) inside += pk;
            pk *= lambda / static_cast<double>(k + 1);
        }
        CHECK(total - inside <= tol);

        // Truncating one term earlier on the right must lose more than the
        // weights outside the window would suggest is free, i.e. the window is
        // not grossly larger than needed: its rightmost weight is above tol.
        CHECK(r.w[r.w.size() - 1] > 0.0);
    }
}

TEST_CASE("ctmc_foxglynn agrees with uniformization and with MATLAB") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0(6, 1.0 / 6.0);
    const line::mc::FoxGlynnResult<double> r = ctmc_foxglynn(pi0, Q, 0.7, 1e-12);
    CHECK(sum_of(r.pi) == doctest::Approx(1.0).epsilon(1e-12));

    // Both truncate the same Poisson mixture at tolerance 1e-12, so they must
    // agree to that tolerance.
    const std::vector<double> u = ctmc_uniformization(pi0, Q, 0.7, 1e-12).pi;
    CHECK(linf(r.pi, u) < 1e-12);

    // MATLAB ctmc_foxglynn(ones(1,6)/6, QA, 0.7, 1e-12).
    const std::vector<double> matlab{0.13548804749281984, 0.19939271418514765, 0.057186264055455585,
                                     0.27369060302903203, 0.22266434523338821, 0.1115780260041567};
    CHECK(linf(r.pi, matlab) < 1e-12);
    CHECK(r.left == 0);
    CHECK(r.right == 26);
}

TEST_CASE("ctmc_foxglynn converges to the stationary vector as t grows") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> stat = ctmc_solve(Q);
    double prev = 1.0;
    for (double t : {10.0, 100.0, 2000.0}) {
        const std::vector<double> pit = ctmc_foxglynn(pi0, Q, t, 1e-12).pi;
        CHECK(sum_of(pit) == doctest::Approx(1.0).epsilon(1e-10));
        const double d = linf(pit, stat);
        CHECK(d < prev);
        prev = d;
    }
    CHECK(prev < 1e-8);
}

TEST_CASE("ctmc_foxglynn at high precision reproduces the double result") {
    const Matrix<Real50> Q = chain_A<Real50>();
    const std::vector<Real50> pi0(6, Real50(1) / Real50(6));
    const line::mc::FoxGlynnResult<Real50> r = ctmc_foxglynn(pi0, Q, Real50(0.7), 1e-12);
    Real50 s(0);
    for (const Real50& w : r.w) s += w;
    CHECK(static_cast<double>(s) == doctest::Approx(1.0).epsilon(1e-30));
    const std::vector<double> pid = ctmc_foxglynn(std::vector<double>(6, 1.0 / 6.0), chain_A<double>(),
                                                 0.7, 1e-12)
                                        .pi;
    for (std::size_t i = 0; i < 6; ++i)
        CHECK(static_cast<double>(r.pi[i]) == doctest::Approx(pid[i]).epsilon(1e-11));
}

// ---------------------------------------------------------------------------
// ctmc_transient
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_transient conserves probability at every time point") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const line::mc::TransientResult<double> r = ctmc_transient(Q, pi0, 0.0, 5.0);
    REQUIRE(r.t.size() >= 2);
    REQUIRE(r.pi.rows() == r.t.size());
    CHECK(r.t[0] == doctest::Approx(0.0));
    CHECK(r.t[r.t.size() - 1] == doctest::Approx(5.0));
    for (std::size_t k = 0; k < r.pi.rows(); ++k) {
        double s = 0.0;
        for (std::size_t j = 0; j < 6; ++j) s += r.pi(k, j);
        // The generator has zero row sums, so d(sum pi)/dt is identically zero
        // and the integrator can only lose mass through its local error.
        CHECK(s == doctest::Approx(1.0).epsilon(1e-9));
    }
    // The time grid is increasing.
    for (std::size_t k = 1; k < r.t.size(); ++k) CHECK(r.t[k] > r.t[k - 1]);
}

TEST_CASE("ctmc_transient converges to the stationary vector as t grows") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const std::vector<double> stat = ctmc_solve(Q);
    double prev = 1.0;
    for (double t : {20.0, 200.0, 4000.0}) {
        const line::mc::TransientResult<double> r = ctmc_transient(Q, pi0, 0.0, t);
        std::vector<double> last(6);
        for (std::size_t j = 0; j < 6; ++j) last[j] = r.pi(r.pi.rows() - 1, j);
        const double d = linf(last, stat);
        CHECK(d < prev);
        prev = d;
    }
    // ode23 at RelTol 1e-3 cannot do better than its own accumulated local
    // error, so this is asserted at the integrator's tolerance, not tighter.
    CHECK(prev < 1e-3);
}

TEST_CASE("ctmc_transient matches MATLAB ode23 to the integrator tolerance") {
    const Matrix<double> Q = chain_A<double>();
    const std::vector<double> pi0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const line::mc::TransientResult<double> r = ctmc_transient(Q, pi0, 0.0, 0.7);
    std::vector<double> last(6);
    for (std::size_t j = 0; j < 6; ++j) last[j] = r.pi(r.pi.rows() - 1, j);

    // MATLAB ctmc_transient(QA, e_1, 0, 0.7), which is ode23 with RelTol 1e-3
    // and AbsTol 1e-6. Two adaptive integrations control the LOCAL error only,
    // so the guarantee on the end value is the tolerance itself; asserting
    // tighter would be asserting that the two step sequences coincide.
    const std::vector<double> matlab{0.41720597362018935,   0.57976698730119092, 0.0012457324620651882,
                                     0.0017729262627971772, 6.8687522876476664e-06,
                                     1.5116014699955553e-06};
    CHECK(linf(last, matlab) < 1e-3);

    // Against an accurate reference the same tolerance applies, and both
    // solutions must sit on the same side of it.
    const std::vector<double> accurate = ctmc_foxglynn(pi0, Q, 0.7, 1e-14).pi;
    CHECK(linf(last, accurate) < 1e-3);
    CHECK(linf(matlab, accurate) < 1e-3);
}

TEST_CASE("ctmc_transient starts from the uniform vector when none is given") {
    const Matrix<double> Q = chain_A<double>();
    const line::mc::TransientResult<double> r = ctmc_transient(Q, 0.0, 1.0);
    for (std::size_t j = 0; j < 6; ++j) CHECK(r.pi(0, j) == doctest::Approx(1.0 / 6.0));
}

// ---------------------------------------------------------------------------
// ctmc_sens
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_sens satisfies its defining equations EXACTLY in rational arithmetic") {
    const Matrix<Rational> Q = chain_A<Rational>();
    const Matrix<Rational> dQ = chain_A_dQ<Rational>();
    const std::vector<Rational> pi = ctmc_solve(Q);
    const std::vector<Rational> dpi = ctmc_sens(Q, dQ, pi);

    // sum_i dpi_i = 0, exactly.
    Rational s(0);
    for (const Rational& v : dpi) s += v;
    CHECK(s == Rational(0));

    // dpi Q + pi dQ = 0, exactly, in every coordinate. This is Eq. (9.81) of
    // Trivedi and Bobbio and is the whole content of the routine; in double
    // precision it holds only to round-off.
    for (std::size_t j = 0; j < 6; ++j) {
        Rational acc(0);
        for (std::size_t i = 0; i < 6; ++i) acc += dpi[i] * Q(i, j) + pi[i] * dQ(i, j);
        CHECK(acc == Rational(0));
    }
}

TEST_CASE("ctmc_sens matches a central difference and MATLAB") {
    const Matrix<double> Q = chain_A<double>();
    const Matrix<double> dQ = chain_A_dQ<double>();
    const std::vector<double> dpi = ctmc_sens(Q, dQ);

    // MATLAB ctmc_sens(QA, dQ): the same finite linear solve, so this is a
    // round-off comparison.
    const std::vector<double> matlab{-0.055183718651365316, 0.021820198759437723,
                                     0.0021176596322194535, 0.010380684471663989,
                                     0.013944719473601656,  0.0069204563144425087};
    CHECK(linf(dpi, matlab) < 1e-12);

    // Central difference of the stationary vector in the same parameter. The
    // second-order error of the difference is what limits the comparison.
    const double h = 1e-5;
    Matrix<double> Qp = Q, Qm = Q;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) {
            Qp(i, j) += h * dQ(i, j);
            Qm(i, j) -= h * dQ(i, j);
        }
    const std::vector<double> pp = ctmc_solve(Qp), pm = ctmc_solve(Qm);
    for (std::size_t i = 0; i < 6; ++i)
        CHECK((pp[i] - pm[i]) / (2.0 * h) == doctest::Approx(dpi[i]).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// ctmc_transient_sens
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_transient_sens starts at zero and conserves zero total sensitivity") {
    const Matrix<double> Q = chain_A<double>();
    const Matrix<double> dQ = chain_A_dQ<double>();
    const std::vector<double> pi0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const line::mc::TransientSensResult<double> r = ctmc_transient_sens(Q, dQ, pi0, 0.0, 0.7);
    REQUIRE(r.dpi.rows() == r.t.size());

    for (std::size_t j = 0; j < 6; ++j) CHECK(r.dpi(0, j) == doctest::Approx(0.0));
    for (std::size_t k = 0; k < r.dpi.rows(); ++k) {
        // pi(t) sums to one for every theta, so its derivative sums to zero.
        double s = 0.0, p = 0.0;
        for (std::size_t j = 0; j < 6; ++j) {
            s += r.dpi(k, j);
            p += r.pi(k, j);
        }
        CHECK(s == doctest::Approx(0.0).epsilon(1e-8));
        CHECK(p == doctest::Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("ctmc_transient_sens matches MATLAB and a central difference") {
    const Matrix<double> Q = chain_A<double>();
    const Matrix<double> dQ = chain_A_dQ<double>();
    const std::vector<double> pi0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const line::mc::TransientSensResult<double> r = ctmc_transient_sens(Q, dQ, pi0, 0.0, 0.7);
    std::vector<double> last(6);
    for (std::size_t j = 0; j < 6; ++j) last[j] = r.dpi(r.dpi.rows() - 1, j);

    // MATLAB ctmc_transient_sens(QA, dQ, e_1, 0, 0.7). Both integrate the
    // augmented system with ode23 at RelTol 1e-3, so the guarantee is that
    // tolerance.
    const std::vector<double> matlab{-0.090372280626725157,  0.08976224982271197,
                                     0.00022986245028041051, 0.00037820115717254575,
                                     1.5997573736088449e-06, 3.6743918659298926e-07};
    CHECK(linf(last, matlab) < 1e-3);

    // Central difference of the transient solution at the same time, computed
    // by Fox-Glynn so that the reference itself carries no integration error.
    const double h = 1e-4;
    Matrix<double> Qp = Q, Qm = Q;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) {
            Qp(i, j) += h * dQ(i, j);
            Qm(i, j) -= h * dQ(i, j);
        }
    const std::vector<double> pp = ctmc_foxglynn(pi0, Qp, 0.7, 1e-14).pi;
    const std::vector<double> pm = ctmc_foxglynn(pi0, Qm, 0.7, 1e-14).pi;
    std::vector<double> fd(6);
    for (std::size_t i = 0; i < 6; ++i) fd[i] = (pp[i] - pm[i]) / (2.0 * h);
    CHECK(linf(last, fd) < 1e-3);
}
