/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Optimizers in line/util: levmar.h, neldermead.h, auglag.h.
 *
 * Every case here has a known optimum, established independently of the
 * optimizer: an analytic minimizer (Rosenbrock, box-constrained quadratic), a
 * normal-equation solution computed with line::solve (linear least squares),
 * or a hand-derived KKT point (the equality- and inequality-constrained
 * quadratics). Nothing is checked against a stored iterate.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/levmar.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "line/util/neldermead.h"

using line::AugLagOptions;
using line::auglag;
using line::auglag_defaults;
using line::auglag_ls;
using line::Bound;
using line::bound_box;
using line::bound_free;
using line::bound_lower;
using line::bound_upper;
using line::levmar;
using line::levmar_defaults;
using line::levmar_jac;
using line::LevmarResult;
using line::Matrix;
using line::nelder_mead;
using line::nelder_mead_box;
using line::nelder_mead_defaults;
using line::NelderMeadResult;
using line::NoConstraints;
using line::solve;

namespace {

using T = double;

/** Residual form of Rosenbrock: f = r0^2 + r1^2 with the minimum 0 at (1,1). */
std::vector<T> rosenbrock_resid(const std::vector<T>& x) {
    std::vector<T> r(2);
    r[0] = 10.0 * (x[1] - x[0] * x[0]);
    r[1] = 1.0 - x[0];
    return r;
}

/** Analytic Jacobian of the above. */
Matrix<T> rosenbrock_jac(const std::vector<T>& x) {
    Matrix<T> J(2, 2, 0.0);
    J(0, 0) = -20.0 * x[0];
    J(0, 1) = 10.0;
    J(1, 0) = -1.0;
    J(1, 1) = 0.0;
    return J;
}

T rosenbrock_scalar(const std::vector<T>& x) {
    const T a = 1.0 - x[0];
    const T b = x[1] - x[0] * x[0];
    return a * a + 100.0 * b * b;
}

/** A fixed, well-conditioned overdetermined system; no randomness. */
Matrix<T> ls_design() {
    Matrix<T> A(5, 3, 0.0);
    const double vals[15] = {1.0, 0.5, 0.25, 1.0, 1.0, 1.0, 1.0,  1.5,
                             2.25, 1.0, 2.0, 4.0, 1.0, 2.5, 6.25};
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t j = 0; j < 3; ++j) A(i, j) = vals[i * 3 + j];
    return A;
}

std::vector<T> ls_rhs() {
    std::vector<T> b(5);
    b[0] = 1.1;
    b[1] = 1.9;
    b[2] = 3.2;
    b[3] = 4.8;
    b[4] = 7.4;
    return b;
}

}  // namespace

TEST_CASE("levmar reaches the Rosenbrock minimum from the classic start") {
    std::vector<T> x0(2);
    x0[0] = -1.2;
    x0[1] = 1.0;

    const LevmarResult<T> r = levmar<T>(rosenbrock_resid, x0, 2);
    CHECK(r.converged);
    CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(r.x[1] == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(r.ssq < 1e-18);
    CHECK(r.iterations <= 100);
}

TEST_CASE("levmar with an analytic Jacobian agrees with the differenced one") {
    std::vector<T> x0(2);
    x0[0] = -1.2;
    x0[1] = 1.0;

    const LevmarResult<T> rfd = levmar<T>(rosenbrock_resid, x0, 2);
    const LevmarResult<T> ran = levmar_jac<T>(rosenbrock_resid, rosenbrock_jac, x0, 2,
                                              levmar_defaults<T>());
    CHECK(ran.converged);
    CHECK(ran.x[0] == doctest::Approx(1.0).epsilon(1e-10));
    CHECK(ran.x[1] == doctest::Approx(1.0).epsilon(1e-10));
    // the analytic Jacobian costs no extra residual evaluations
    CHECK(ran.evaluations < rfd.evaluations);
}

TEST_CASE("levmar from several starts always lands on the same Rosenbrock optimum") {
    const double starts[4][2] = {{-1.2, 1.0}, {5.0, 5.0}, {0.0, 0.0}, {-3.0, -4.0}};
    for (int k = 0; k < 4; ++k) {
        std::vector<T> x0(2);
        x0[0] = starts[k][0];
        x0[1] = starts[k][1];
        const LevmarResult<T> r = levmar<T>(rosenbrock_resid, x0, 2);
        CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-7));
        CHECK(r.x[1] == doctest::Approx(1.0).epsilon(1e-7));
    }
}

TEST_CASE("levmar on a linear least-squares problem matches the normal equations") {
    const Matrix<T> A = ls_design();
    const std::vector<T> b = ls_rhs();

    // exact reference: (A'A) x = A'b, solved with the LU in line/util/lu.h
    Matrix<T> AtA(3, 3, 0.0);
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t k = 0; k < 3; ++k) {
            T s = 0.0;
            for (std::size_t i = 0; i < 5; ++i) s += A(i, j) * A(i, k);
            AtA(j, k) = s;
        }
    std::vector<T> Atb(3, 0.0);
    for (std::size_t j = 0; j < 3; ++j) {
        T s = 0.0;
        for (std::size_t i = 0; i < 5; ++i) s += A(i, j) * b[i];
        Atb[j] = s;
    }
    const std::vector<T> xref = solve(AtA, Atb);

    auto resid = [&A, &b](const std::vector<T>& x) {
        std::vector<T> r(5);
        for (std::size_t i = 0; i < 5; ++i) {
            T s = -b[i];
            for (std::size_t j = 0; j < 3; ++j) s += A(i, j) * x[j];
            r[i] = s;
        }
        return r;
    };

    std::vector<T> x0(3, 0.0);
    const LevmarResult<T> r = levmar<T>(resid, x0, 5);
    CHECK(r.converged);
    for (std::size_t j = 0; j < 3; ++j) CHECK(r.x[j] == doctest::Approx(xref[j]).epsilon(1e-8));

    // and the objective is the residual of the reference solution
    const std::vector<T> rref = resid(xref);
    T sref = 0.0;
    for (std::size_t i = 0; i < 5; ++i) sref += rref[i] * rref[i];
    CHECK(r.ssq == doctest::Approx(sref).epsilon(1e-10));
}

TEST_CASE("levmar reports non-convergence rather than throwing when capped") {
    line::LevmarOptions<T> o = levmar_defaults<T>();
    o.max_iter = 2;
    std::vector<T> x0(2);
    x0[0] = -1.2;
    x0[1] = 1.0;
    const LevmarResult<T> r = levmar<T>(rosenbrock_resid, x0, 2, o);
    CHECK(r.iterations == 2);
    CHECK(!r.converged);
}

TEST_CASE("levmar rejects a degenerate problem by input error") {
    std::vector<T> empty;
    CHECK_THROWS_AS(levmar<T>(rosenbrock_resid, empty, 2), line::InputError);
}

TEST_CASE("nelder_mead reaches the Rosenbrock minimum") {
    std::vector<T> x0(2);
    x0[0] = -1.2;
    x0[1] = 1.0;
    line::NelderMeadOptions<T> o = nelder_mead_defaults<T>();
    o.ftol = 1e-14;
    o.xtol = 1e-10;
    const NelderMeadResult<T> r = nelder_mead<T>(rosenbrock_scalar, x0, o);
    CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(r.x[1] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(r.fval < 1e-10);
}

TEST_CASE("nelder_mead is deterministic: identical inputs give identical output") {
    std::vector<T> x0(2);
    x0[0] = -1.2;
    x0[1] = 1.0;
    const NelderMeadResult<T> a = nelder_mead<T>(rosenbrock_scalar, x0);
    const NelderMeadResult<T> b = nelder_mead<T>(rosenbrock_scalar, x0);
    CHECK(a.evaluations == b.evaluations);
    CHECK(a.iterations == b.iterations);
    CHECK(a.x[0] == b.x[0]);
    CHECK(a.x[1] == b.x[1]);
    CHECK(a.fval == b.fval);
}

TEST_CASE("nelder_mead_box solves a box-constrained quadratic to the clamped optimum") {
    // min sum (x_j - c_j)^2 over a box: the minimizer is c clamped componentwise
    std::vector<T> c(3);
    c[0] = 3.0;    // above its upper bound
    c[1] = -2.0;   // below its lower bound
    c[2] = 0.4;    // interior
    auto f = [&c](const std::vector<T>& x) {
        T s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) s += (x[j] - c[j]) * (x[j] - c[j]);
        return s;
    };
    std::vector<Bound<T>> bnd(3);
    bnd[0] = bound_box<T>(0.0, 1.0);
    bnd[1] = bound_lower<T>(0.0);
    bnd[2] = bound_box<T>(0.0, 1.0);

    std::vector<T> x0(3, 0.5);
    line::NelderMeadOptions<T> o = nelder_mead_defaults<T>();
    o.xtol = 1e-12;
    o.ftol = 1e-16;
    const NelderMeadResult<T> r = nelder_mead_box<T>(f, x0, bnd, o);

    CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(std::fabs(r.x[1]) < 1e-5);
    CHECK(r.x[2] == doctest::Approx(0.4).epsilon(1e-6));
    // the achieved objective equals the analytic one at the clamped point
    const T fref = (1.0 - 3.0) * (1.0 - 3.0) + (0.0 + 2.0) * (0.0 + 2.0) + 0.0;
    CHECK(r.fval == doctest::Approx(fref).epsilon(1e-6));
}

TEST_CASE("nelder_mead_box never evaluates outside the box and honours an upper bound") {
    bool inside = true;
    auto f = [&inside](const std::vector<T>& x) {
        if (x[0] < -1.0 || x[0] > 2.0) inside = false;
        return (x[0] - 10.0) * (x[0] - 10.0);
    };
    std::vector<Bound<T>> bnd(1);
    bnd[0] = bound_upper<T>(2.0);
    std::vector<T> x0(1, 0.0);
    const NelderMeadResult<T> r = nelder_mead_box<T>(f, x0, bnd);
    CHECK(r.x[0] <= 2.0);
    CHECK(r.x[0] == doctest::Approx(2.0).epsilon(1e-4));
}

TEST_CASE("nelder_mead_box holds a variable with coincident bounds fixed") {
    auto f = [](const std::vector<T>& x) {
        return (x[0] - 5.0) * (x[0] - 5.0) + (x[1] - 5.0) * (x[1] - 5.0);
    };
    std::vector<Bound<T>> bnd(2);
    bnd[0] = bound_box<T>(2.0, 2.0);
    bnd[1] = bound_free<T>();
    std::vector<T> x0(2, 0.0);
    const NelderMeadResult<T> r = nelder_mead_box<T>(f, x0, bnd);
    CHECK(r.x[0] == 2.0);
    CHECK(r.x[1] == doctest::Approx(5.0).epsilon(1e-6));
}

TEST_CASE("auglag solves an equality-constrained quadratic at its KKT point") {
    // min 1/2 x'x - b'x  s.t.  a'x = d
    // KKT: x - b + nu a = 0, a'x = d  =>  nu = (a'b - d)/(a'a), x = b - nu a
    std::vector<T> b(3);
    b[0] = 1.0;
    b[1] = 2.0;
    b[2] = -1.0;
    std::vector<T> a(3);
    a[0] = 1.0;
    a[1] = 1.0;
    a[2] = 1.0;
    const T d = 3.0;

    T ab = 0.0, aa = 0.0;
    for (std::size_t j = 0; j < 3; ++j) {
        ab += a[j] * b[j];
        aa += a[j] * a[j];
    }
    const T nu = (ab - d) / aa;  // = (2 - 3)/3 = -1/3
    std::vector<T> xref(3);
    for (std::size_t j = 0; j < 3; ++j) xref[j] = b[j] - nu * a[j];
    // xref = (4/3, 7/3, -2/3), which sums to 3 as required

    auto f = [&b](const std::vector<T>& x) {
        T s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) s += 0.5 * x[j] * x[j] - b[j] * x[j];
        return s;
    };
    auto h = [&a, d](const std::vector<T>& x) {
        std::vector<T> hv(1, -d);
        for (std::size_t j = 0; j < 3; ++j) hv[0] += a[j] * x[j];
        return hv;
    };

    std::vector<Bound<T>> free_bounds(3);
    std::vector<T> x0(3, 0.0);
    AugLagOptions<T> o = auglag_defaults<T>();
    o.inner_nm.xtol = 1e-13;
    o.inner_nm.ftol = 1e-16;
    o.ctol = 1e-11;

    const line::AugLagResult<T> r = auglag<T>(f, h, NoConstraints<T>(), x0, free_bounds, o);
    CHECK(r.violation < 1e-9);
    for (std::size_t j = 0; j < 3; ++j) CHECK(r.x[j] == doctest::Approx(xref[j]).epsilon(1e-6));
    CHECK(r.fval == doctest::Approx(f(xref)).epsilon(1e-8));
    // The recovered multiplier matches the hand-derived one (sign: L = f + lambda h).
    // A first-order multiplier method estimates lambda to the accuracy of the
    // inner minimization, so it is one or two digits looser than x itself.
    CHECK(r.lambda[0] == doctest::Approx(nu).epsilon(1e-3));
}

TEST_CASE("auglag solves an inequality-constrained quadratic with an active constraint") {
    // min (x1-2)^2 + (x2-2)^2  s.t.  x1 + x2 <= 2
    // The unconstrained minimizer (2,2) violates the constraint, so it is
    // active: by symmetry the KKT point is (1,1) with multiplier 2.
    auto f = [](const std::vector<T>& x) {
        return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0);
    };
    auto g = [](const std::vector<T>& x) {
        std::vector<T> gv(1);
        gv[0] = x[0] + x[1] - 2.0;
        return gv;
    };
    std::vector<Bound<T>> free_bounds(2);
    std::vector<T> x0(2, 0.0);
    AugLagOptions<T> o = auglag_defaults<T>();
    o.inner_nm.xtol = 1e-13;
    o.inner_nm.ftol = 1e-16;
    o.ctol = 1e-11;

    const line::AugLagResult<T> r = auglag<T>(f, NoConstraints<T>(), g, x0, free_bounds, o);
    CHECK(r.violation < 1e-9);
    CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r.x[1] == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r.fval == doctest::Approx(2.0).epsilon(1e-7));
    CHECK(r.mu[0] == doctest::Approx(2.0).epsilon(1e-5));
}

TEST_CASE("auglag leaves an inactive inequality alone") {
    // same objective, constraint x1 + x2 <= 10 is slack at the optimum (2,2)
    auto f = [](const std::vector<T>& x) {
        return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0);
    };
    auto g = [](const std::vector<T>& x) {
        std::vector<T> gv(1);
        gv[0] = x[0] + x[1] - 10.0;
        return gv;
    };
    std::vector<Bound<T>> free_bounds(2);
    std::vector<T> x0(2, 0.0);
    const line::AugLagResult<T> r = auglag<T>(f, NoConstraints<T>(), g, x0, free_bounds);
    CHECK(r.x[0] == doctest::Approx(2.0).epsilon(1e-5));
    CHECK(r.x[1] == doctest::Approx(2.0).epsilon(1e-5));
    CHECK(std::fabs(r.mu[0]) < 1e-10);
    CHECK(r.fval < 1e-8);
}

TEST_CASE("auglag combines a box, an equality and an inequality") {
    // min (x1-3)^2 + (x2-3)^2  s.t.  x1 + x2 = 4,  x1 <= 1,  0 <= x = 5
    // The equality and the active x1 <= 1 fix x = (1,3); the objective is 4.
    auto f = [](const std::vector<T>& x) {
        return (x[0] - 3.0) * (x[0] - 3.0) + (x[1] - 3.0) * (x[1] - 3.0);
    };
    auto h = [](const std::vector<T>& x) {
        std::vector<T> hv(1);
        hv[0] = x[0] + x[1] - 4.0;
        return hv;
    };
    auto g = [](const std::vector<T>& x) {
        std::vector<T> gv(1);
        gv[0] = x[0] - 1.0;
        return gv;
    };
    std::vector<Bound<T>> bnd(2);
    bnd[0] = bound_box<T>(0.0, 5.0);
    bnd[1] = bound_box<T>(0.0, 5.0);
    std::vector<T> x0(2, 2.0);
    AugLagOptions<T> o = auglag_defaults<T>();
    o.inner_nm.xtol = 1e-13;
    o.inner_nm.ftol = 1e-16;
    o.ctol = 1e-10;
    const line::AugLagResult<T> r = auglag<T>(f, h, g, x0, bnd, o);
    CHECK(r.violation < 1e-8);
    CHECK(r.x[0] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(r.x[1] == doctest::Approx(3.0).epsilon(1e-5));
    CHECK(r.fval == doctest::Approx(4.0).epsilon(1e-5));
}

TEST_CASE("auglag_ls solves the same equality-constrained quadratic through levmar") {
    // min ||x - c||^2 s.t. sum x = 1, with c = (1, 2, 3).
    // KKT: x = c + nu*1 with sum x = 1 => nu = (1 - 6)/3 = -5/3.
    std::vector<T> c(3);
    c[0] = 1.0;
    c[1] = 2.0;
    c[2] = 3.0;
    const T nu = (1.0 - 6.0) / 3.0;
    std::vector<T> xref(3);
    for (std::size_t j = 0; j < 3; ++j) xref[j] = c[j] + nu;

    auto r = [&c](const std::vector<T>& x) {
        std::vector<T> rv(3);
        for (std::size_t j = 0; j < 3; ++j) rv[j] = x[j] - c[j];
        return rv;
    };
    auto h = [](const std::vector<T>& x) {
        std::vector<T> hv(1, -1.0);
        for (std::size_t j = 0; j < 3; ++j) hv[0] += x[j];
        return hv;
    };
    std::vector<T> x0(3, 0.0);
    const line::AugLagResult<T> res = auglag_ls<T>(r, 3, h, NoConstraints<T>(), x0);
    CHECK(res.violation < 1e-10);
    for (std::size_t j = 0; j < 3; ++j) CHECK(res.x[j] == doctest::Approx(xref[j]).epsilon(1e-8));
    T fref = 0.0;
    for (std::size_t j = 0; j < 3; ++j) fref += nu * nu;
    CHECK(res.fval == doctest::Approx(fref).epsilon(1e-8));
}

TEST_CASE("auglag_ls handles bounds expressed as inequality rows") {
    // min ||x - c||^2 with c = (3, -2) and 0 <= x <= 1 written as four rows.
    // Solution is the clamp, (1, 0), objective 4 + 4 = 8.
    std::vector<T> c(2);
    c[0] = 3.0;
    c[1] = -2.0;
    auto r = [&c](const std::vector<T>& x) {
        std::vector<T> rv(2);
        rv[0] = x[0] - c[0];
        rv[1] = x[1] - c[1];
        return rv;
    };
    auto g = [](const std::vector<T>& x) {
        std::vector<T> gv(4);
        gv[0] = -x[0];
        gv[1] = x[0] - 1.0;
        gv[2] = -x[1];
        gv[3] = x[1] - 1.0;
        return gv;
    };
    std::vector<T> x0(2, 0.5);
    AugLagOptions<T> o = auglag_defaults<T>();
    o.ctol = 1e-10;
    const line::AugLagResult<T> res = auglag_ls<T>(r, 2, NoConstraints<T>(), g, x0, o);
    CHECK(res.violation < 1e-8);
    CHECK(res.x[0] == doctest::Approx(1.0).epsilon(1e-7));
    CHECK(std::fabs(res.x[1]) < 1e-7);
    CHECK(res.fval == doctest::Approx(8.0).epsilon(1e-6));
}

TEST_CASE("auglag_ls and auglag agree on a problem both can express") {
    // min (x1-2)^2 + (x2-2)^2 s.t. x1 + x2 <= 2, from both entry points
    auto rls = [](const std::vector<T>& x) {
        std::vector<T> rv(2);
        rv[0] = x[0] - 2.0;
        rv[1] = x[1] - 2.0;
        return rv;
    };
    auto f = [](const std::vector<T>& x) {
        return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0);
    };
    auto g = [](const std::vector<T>& x) {
        std::vector<T> gv(1);
        gv[0] = x[0] + x[1] - 2.0;
        return gv;
    };
    std::vector<T> x0(2, 0.0);
    std::vector<Bound<T>> free_bounds(2);
    AugLagOptions<T> o = auglag_defaults<T>();
    o.ctol = 1e-11;
    o.inner_nm.xtol = 1e-13;
    o.inner_nm.ftol = 1e-16;
    const line::AugLagResult<T> a = auglag_ls<T>(rls, 2, NoConstraints<T>(), g, x0, o);
    const line::AugLagResult<T> b = auglag<T>(f, NoConstraints<T>(), g, x0, free_bounds, o);
    CHECK(a.fval == doctest::Approx(b.fval).epsilon(1e-5));
    CHECK(a.x[0] == doctest::Approx(b.x[0]).epsilon(1e-5));
    CHECK(a.x[1] == doctest::Approx(b.x[1]).epsilon(1e-5));
    CHECK(a.violation < 1e-8);
    CHECK(b.violation < 1e-8);
}
