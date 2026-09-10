/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The adaptive stiff integrator of line/util/ode.h.
 *
 * The oracles are problems whose solution is known in closed form -- scalar
 * exponential decay, a linear system whose flow is a matrix exponential that
 * expm.h computes independently, and a non-autonomous problem with an
 * elementary solution -- plus two properties that no closed form is needed for:
 * the Robertson stiff problem conserves mass exactly, and a problem with a
 * 1e6:1 spread of time constants must cost a number of steps that an explicit
 * method could not come close to. The coefficient set is checked against the
 * conditions it was derived from, so a mistyped digit cannot survive.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/util/expm.h"
#include "line/util/matrix.h"
#include "line/util/ode.h"

using line::Matrix;
using line::OdeOptions;
using line::Real50;

namespace {

/** R(z) of the method's stability function, evaluated by running one step on
 *  y' = z y with h = 1 and y0 = 1. */
double stability_function(double z) {
    const auto f = [z](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{z * y[0]};
    };
    const auto J = [z](const double& t, const std::vector<double>& y) {
        (void)t;
        (void)y;
        Matrix<double> M(1, 1, z);
        return M;
    };
    OdeOptions<double> opt;
    opt.h_init = 1.0;
    opt.h_min = 1.0;
    opt.h_max = 1.0;
    opt.rtol = 1e100;  // accept the single step whatever the estimate says
    opt.atol = 1e100;
    opt.store_trajectory = false;
    const std::vector<double> y0{1.0};
    return line::ode_rosenbrock4(f, J, 0.0, 1.0, y0, opt).final_state()[0];
}

}  // namespace

TEST_CASE("ode: the coefficient set satisfies the conditions it was derived from") {
    const line::ode_detail::Ros4<double> C;

    // The abscissae are the row sums of a by construction; what has to hold
    // beyond that is that they are finite and ordered sensibly, and that the
    // first stage is explicit in y.
    CHECK(C.alpha[0] == 0.0);
    for (int i = 1; i < 4; ++i) {
        CHECK(C.alpha[i] > 0.0);
        CHECK(C.alpha[i] < 1.5);
    }
    CHECK(C.gamma > 0.0);

    // The linear order conditions. On y' = J y the method IS the implicit
    // Runge-Kutta with matrix B = a + gamma (constant diagonal gamma), so
    // b^T B^(k-1) 1 = 1/k! for k = 1..4 is necessary and sufficient for order
    // four on linear problems, and these four are checked directly.
    Matrix<double> B(4, 4, 0.0);
    for (int i = 0; i < 4; ++i) {
        B(i, i) = C.gamma;
        for (int j = 0; j < i; ++j) B(i, j) = C.a[i][j] + C.gam[i][j];
    }
    std::vector<double> v(4, 1.0);
    const double fact[4] = {1.0, 2.0, 6.0, 24.0};
    for (int k = 0; k < 4; ++k) {
        double s = 0.0;
        for (int i = 0; i < 4; ++i) s += C.b[i] * v[i];
        CHECK(s == doctest::Approx(1.0 / fact[k]).epsilon(1e-12));
        std::vector<double> w(4, 0.0);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) w[i] += B(i, j) * v[j];
        v = w;
    }

    // The embedded weights are a genuine lower-order formula: consistent
    // (weights sum to one) and different from b.
    double sb = 0.0, diff = 0.0;
    for (int i = 0; i < 4; ++i) {
        sb += C.bhat[i];
        diff += std::fabs(C.b[i] - C.bhat[i]);
    }
    CHECK(sb == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(diff > 1e-3);
}

TEST_CASE("ode: the method is A-stable and L-stable") {
    // |R(z)| <= 1 along the negative real axis, and R -> 0 at infinity. An
    // A-stable but not L-stable method would leave the fastest modes ringing
    // at |R| = 1 instead of damping them, which defeats the purpose here.
    for (double z : {-0.01, -0.1, -1.0, -3.0, -10.0, -100.0, -1e4, -1e8, -1e12})
        CHECK(std::fabs(stability_function(z)) <= 1.0 + 1e-12);
    CHECK(std::fabs(stability_function(-1e10)) < 1e-6);
    CHECK(std::fabs(stability_function(-1e14)) < 1e-10);
}

TEST_CASE("ode: exponential decay against the closed form") {
    const double lambda = 3.5;
    const auto f = [lambda](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{-lambda * y[0]};
    };
    OdeOptions<double> opt;
    opt.rtol = 1e-10;
    opt.atol = 1e-14;
    const std::vector<double> y0{2.0};
    const line::OdeSolution<double> s = line::ode_rosenbrock4(f, 0.0, 2.0, y0, opt);
    CHECK(s.final_time() == doctest::Approx(2.0).epsilon(1e-14));
    CHECK(s.final_state()[0] == doctest::Approx(2.0 * std::exp(-lambda * 2.0)).epsilon(1e-9));
    // Every stored point must lie on the solution, not just the endpoint.
    for (std::size_t i = 0; i < s.t.size(); ++i)
        CHECK(s.y[i][0] == doctest::Approx(2.0 * std::exp(-lambda * s.t[i])).epsilon(1e-8));
}

TEST_CASE("ode: linear system against expm") {
    // A stiff-ish non-symmetric 3x3 with eigenvalues of very different size;
    // the flow is exp(A t) y0 and expm.h computes it independently.
    Matrix<double> A(3, 3, 0.0);
    A(0, 0) = -100.0; A(0, 1) = 1.0;   A(0, 2) = 0.0;
    A(1, 0) = 2.0;    A(1, 1) = -3.0;  A(1, 2) = 0.5;
    A(2, 0) = 0.0;    A(2, 1) = 1.0;   A(2, 2) = -0.25;

    const auto f = [&A](const double& t, const std::vector<double>& y) {
        (void)t;
        std::vector<double> dy(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) dy[i] += A(i, j) * y[j];
        return dy;
    };
    const auto J = [&A](const double& t, const std::vector<double>& y) {
        (void)t;
        (void)y;
        return A;
    };

    const std::vector<double> y0{1.0, -0.5, 2.0};
    for (double tf : {0.05, 0.5, 4.0}) {
        Matrix<double> At(3, 3, 0.0);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) At(i, j) = A(i, j) * tf;
        const Matrix<double> E = line::expm(At);
        std::vector<double> ref(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) ref[i] += E(i, j) * y0[j];

        OdeOptions<double> opt;
        opt.rtol = 1e-11;
        opt.atol = 1e-14;
        opt.store_trajectory = false;
        const std::vector<double> y = line::ode_rosenbrock4(f, J, 0.0, tf, y0, opt).final_state();
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(y[i] == doctest::Approx(ref[i]).epsilon(1e-8));
    }
}

TEST_CASE("ode: order four on a nonlinear non-autonomous problem") {
    // y1' = y1 cos t + 0.3 y2^2 + 0.1 t, y2' = -0.5 y2 + 0.2 y1^3 has no
    // closed form, so the reference is the method's own answer at a step size
    // far below the ones being measured. Halving the step must cut the error
    // by about 2^4; the non-autonomous terms make this a test of the h^2 f_t
    // stage term as well, which an autonomous problem would not exercise.
    const auto f = [](const double& t, const std::vector<double>& y) {
        return std::vector<double>{y[0] * std::cos(t) + 0.3 * y[1] * y[1] + 0.1 * t,
                                   -0.5 * y[1] + 0.2 * y[0] * y[0] * y[0]};
    };
    const auto J = [](const double& t, const std::vector<double>& y) {
        Matrix<double> M(2, 2, 0.0);
        M(0, 0) = std::cos(t);
        M(0, 1) = 0.6 * y[1];
        M(1, 0) = 0.6 * y[0] * y[0];
        M(1, 1) = -0.5;
        return M;
    };
    const std::vector<double> y0{1.0, 0.7};
    const double tf = 1.5;

    // The reference is a fixed-step run at h/16 of the finest measured step,
    // so its own error is 16^4 = 65536 times smaller and does not pollute the
    // rates. A tolerance-driven reference would be the wrong tool here: the
    // controller is calibrated on the ORDER-TWO embedded estimate, so a
    // tolerance of 1e-13 would ask for a step of about 1e-4 and a very large
    // number of them, which measures the controller rather than the order.
    OdeOptions<double> ref_opt;
    ref_opt.rtol = 1e100;
    ref_opt.atol = 1e100;
    ref_opt.h_init = tf / 2560;
    ref_opt.h_min = tf / 2560 / 1e6;
    ref_opt.h_max = tf / 2560;
    ref_opt.store_trajectory = false;
    const std::vector<double> ref = line::ode_rosenbrock4(f, J, 0.0, tf, y0, ref_opt).final_state();

    // Fixed-step runs: pin h_init = h_min = h_max and make the tolerance
    // unreachable-free so no step is rejected.
    std::vector<double> errs;
    for (int nst : {20, 40, 80, 160}) {
        OdeOptions<double> opt;
        opt.rtol = 1e100;
        opt.atol = 1e100;
        opt.h_init = tf / nst;
        opt.h_min = tf / nst / 1e6;
        opt.h_max = tf / nst;
        opt.store_trajectory = false;
        const std::vector<double> y = line::ode_rosenbrock4(f, J, 0.0, tf, y0, opt).final_state();
        errs.push_back(std::max(std::fabs(y[0] - ref[0]), std::fabs(y[1] - ref[1])));
    }
    for (std::size_t i = 0; i + 1 < errs.size(); ++i) {
        const double rate = std::log2(errs[i] / errs[i + 1]);
        CHECK(rate > 3.6);
        CHECK(rate < 4.4);
    }
}

TEST_CASE("ode: Robertson conserves mass to 1e-10 over a long horizon") {
    // The Robertson problem, the standard stiff test:
    //   y1' = -0.04 y1 + 1e4 y2 y3
    //   y2' =  0.04 y1 - 1e4 y2 y3 - 3e7 y2^2
    //   y3' =  3e7 y2^2
    // The three rates span nine orders of magnitude and y1+y2+y3 is invariant.
    // Conservation is not imposed anywhere in the integrator, so it is a real
    // check on the step: an unstable or badly damped method loses it long
    // before it loses accuracy.
    const auto f = [](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{-0.04 * y[0] + 1.0e4 * y[1] * y[2],
                                   0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1],
                                   3.0e7 * y[1] * y[1]};
    };
    const auto J = [](const double& t, const std::vector<double>& y) {
        (void)t;
        Matrix<double> M(3, 3, 0.0);
        M(0, 0) = -0.04;
        M(0, 1) = 1.0e4 * y[2];
        M(0, 2) = 1.0e4 * y[1];
        M(1, 0) = 0.04;
        M(1, 1) = -1.0e4 * y[2] - 6.0e7 * y[1];
        M(1, 2) = -1.0e4 * y[1];
        M(2, 1) = 6.0e7 * y[1];
        return M;
    };

    const std::vector<double> y0{1.0, 0.0, 0.0};
    OdeOptions<double> opt;
    opt.rtol = 1e-8;
    opt.atol = 1e-12;
    opt.h_init = 1e-8;
    const line::OdeSolution<double> s = line::ode_rosenbrock4(f, J, 0.0, 1.0e6, y0, opt);

    // Mass conservation along the WHOLE trajectory, not merely at the end.
    double worst = 0.0;
    for (const std::vector<double>& y : s.y)
        worst = std::max(worst, std::fabs(y[0] + y[1] + y[2] - 1.0));
    CHECK(worst < 1e-10);

    // The solution itself, against MATLAB. The oracle is
    //   ode15s(f, [0 1e6], [1;0;0], odeset('RelTol',1e-12,'AbsTol',1e-16,
    //                                      'Jacobian',J))
    // which returns
    //   y1 = 2.031483925146904e-03
    //   y2 = 8.142277784052899e-09
    //   y3 = 9.979685079325762e-01
    // (scipy's Radau at rtol 1e-12 agrees with those to 1e-10 relative, so the
    // oracle is not an artifact of one code). This run asks for rtol 1e-8, so
    // 1e-6 relative is the honest tolerance to hold it to.
    CHECK(s.final_state()[0] == doctest::Approx(2.031483925146904e-03).epsilon(1e-6));
    // `.scale(0.0)` because Approx compares |lhs-rhs| < eps*(scale + max|.|)
    // and scale defaults to 1: on an expected value of 8e-09 the relative
    // epsilon silently became an ABSOLUTE 1e-6, so the row passed for anything
    // in [-1e-6, 1e-6] -- including zero, five orders of magnitude out.
    CHECK(s.final_state()[1] == doctest::Approx(8.142277784052899e-09).epsilon(1e-6).scale(0.0));
    CHECK(s.final_state()[2] == doctest::Approx(9.979685079325762e-01).epsilon(1e-9));
    // Non-negativity is a property of the exact solution and a common casualty
    // of an over-long step.
    for (const std::vector<double>& y : s.y)
        for (double v : y) CHECK(v > -1e-12);
}

TEST_CASE("ode: the stiff path costs what an explicit method could not") {
    // y' = A y with eigenvalues -1e6 and -1: after the fast mode dies the
    // solution is smooth, but an explicit method stays pinned to the fast
    // scale for the whole integration. Forward Euler needs h < 2/1e6 for
    // stability alone, so 1e6 steps over [0,1]; the check is that the stiff
    // integrator does the same job in fewer than a thousand.
    Matrix<double> A(2, 2, 0.0);
    A(0, 0) = -1.0e6;
    A(0, 1) = 1.0e6 - 1.0;
    A(1, 1) = -1.0;
    const auto f = [&A](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{A(0, 0) * y[0] + A(0, 1) * y[1], A(1, 1) * y[1]};
    };
    const auto J = [&A](const double& t, const std::vector<double>& y) {
        (void)t;
        (void)y;
        return A;
    };
    const std::vector<double> y0{1.0, 1.0};
    OdeOptions<double> opt;
    opt.rtol = 1e-8;
    opt.atol = 1e-12;
    opt.store_trajectory = false;
    const line::OdeSolution<double> s = line::ode_rosenbrock4(f, J, 0.0, 1.0, y0, opt);

    // y2 = exp(-t); y1 = exp(-t) exactly for this A and this initial state
    // (the fast mode is not excited), so the answer is known.
    CHECK(s.final_state()[1] == doctest::Approx(std::exp(-1.0)).epsilon(1e-9));
    CHECK(s.final_state()[0] == doctest::Approx(std::exp(-1.0)).epsilon(1e-7));
    CHECK(s.steps < 1000);
    const std::size_t explicit_steps = static_cast<std::size_t>(1.0 / (2.0 / 1.0e6));
    CHECK(explicit_steps > 100 * s.steps);

    // And the explicit alternative really does blow up at the step the stiff
    // method is taking: one forward Euler step of size 1/steps on this A
    // amplifies the fast mode by |1 + h*(-1e6)| >> 1.
    const double h_stiff = 1.0 / static_cast<double>(s.steps);
    CHECK(std::fabs(1.0 + h_stiff * A(0, 0)) > 100.0);
}

TEST_CASE("ode: numeric and analytic Jacobians agree") {
    // The numeric Jacobian is a central difference, so it must agree with the
    // analytic one to about eps^(2/3), and the two integrations must agree to
    // the requested tolerance rather than to that of the difference.
    const auto f = [](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{-2.0 * y[0] + y[0] * y[1], -5.0 * y[1] + 0.5 * y[0] * y[0]};
    };
    const auto J = [](const double& t, const std::vector<double>& y) {
        (void)t;
        Matrix<double> M(2, 2, 0.0);
        M(0, 0) = -2.0 + y[1];
        M(0, 1) = y[0];
        M(1, 0) = y[0];
        M(1, 1) = -5.0;
        return M;
    };
    const std::vector<double> y{0.7, -0.3};
    const Matrix<double> Jn = line::ode_numeric_jacobian<double>(f, 0.0, y, f(0.0, y));
    const Matrix<double> Ja = J(0.0, y);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(Jn(i, j) == doctest::Approx(Ja(i, j)).epsilon(1e-9));

    const std::vector<double> y0{1.0, 1.0};
    OdeOptions<double> opt;
    opt.rtol = 1e-10;
    opt.atol = 1e-13;
    opt.store_trajectory = false;
    const std::vector<double> ya = line::ode_rosenbrock4(f, J, 0.0, 3.0, y0, opt).final_state();
    const std::vector<double> yn = line::ode_rosenbrock4(f, 0.0, 3.0, y0, opt).final_state();
    for (std::size_t i = 0; i < 2; ++i) CHECK(yn[i] == doctest::Approx(ya[i]).epsilon(1e-7));
}

TEST_CASE("ode: step control rejects and recovers") {
    // A problem with a sharp transient: the controller must reject at least
    // one step (the initial guess is deliberately far too large) and still
    // land on the right answer.
    const auto f = [](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{-50.0 * (y[0] - 1.0)};
    };
    OdeOptions<double> opt;
    opt.rtol = 1e-10;
    opt.atol = 1e-14;
    opt.h_init = 0.5;
    const std::vector<double> y0{0.0};
    const line::OdeSolution<double> s = line::ode_rosenbrock4(f, 0.0, 1.0, y0, opt);
    CHECK(s.final_state()[0] == doctest::Approx(1.0 - std::exp(-50.0)).epsilon(1e-9));
    CHECK(s.steps > 1);
}

TEST_CASE("ode: input validation") {
    const auto f = [](const double& t, const std::vector<double>& y) {
        (void)t;
        return std::vector<double>{-y[0]};
    };
    OdeOptions<double> opt;
    const std::vector<double> y0{1.0};
    CHECK_THROWS_AS(line::ode_rosenbrock4(f, 1.0, 0.0, y0, opt), line::InputError);
    CHECK_THROWS_AS(line::ode_rosenbrock4(f, 0.0, 1.0, std::vector<double>(), opt),
                    line::InputError);
    OdeOptions<double> bad = opt;
    bad.rtol = 0.0;
    CHECK_THROWS_AS(line::ode_rosenbrock4(f, 0.0, 1.0, y0, bad), line::InputError);
    OdeOptions<double> tiny = opt;
    tiny.max_steps = 2;
    tiny.h_init = 1e-6;
    tiny.h_max = 1e-6;
    CHECK_THROWS_AS(line::ode_rosenbrock4(f, 0.0, 1.0, y0, tiny), line::NumericError);
}

TEST_CASE("ode: Real50 instantiation") {
    const Real50 lambda("2.5");
    const auto f = [&lambda](const Real50& t, const std::vector<Real50>& y) {
        (void)t;
        return std::vector<Real50>{-lambda * y[0]};
    };
    const auto J = [&lambda](const Real50& t, const std::vector<Real50>& y) {
        (void)t;
        (void)y;
        Matrix<Real50> M(1, 1, Real50(-lambda));
        return M;
    };
    const std::vector<Real50> y0{Real50(1)};
    const Real50 ref = exp(Real50(-lambda));

    SUBCASE("adaptive") {
        OdeOptions<Real50> opt;
        opt.rtol = Real50("1e-14");
        opt.atol = Real50("1e-20");
        opt.store_trajectory = false;
        const std::vector<Real50> y =
            line::ode_rosenbrock4(f, J, Real50(0), Real50(1), y0, opt).final_state();
        CHECK(static_cast<double>(abs(Real50(y[0] - ref))) < 1e-12);
    }

    SUBCASE("fixed step, past what double can represent") {
        // 20000 steps of h = 5e-5: the truncation error of an order-four
        // method is then far below the double unit roundoff, so the answer can
        // only be that accurate if the COEFFICIENTS are accurate to more than
        // double as well. This is the check that the 32-digit strings are
        // doing real work; with coefficients rounded to double the error would
        // stall at about 1e-16.
        OdeOptions<Real50> opt;
        opt.rtol = Real50("1e100");
        opt.atol = Real50("1e100");
        const Real50 h = Real50(1) / Real50(20000);
        opt.h_init = h;
        opt.h_min = h / Real50(1000);
        opt.h_max = h;
        opt.store_trajectory = false;
        const std::vector<Real50> y =
            line::ode_rosenbrock4(f, J, Real50(0), Real50(1), y0, opt).final_state();
        CHECK(static_cast<double>(abs(Real50(y[0] - ref))) < 1e-17);
    }
}
