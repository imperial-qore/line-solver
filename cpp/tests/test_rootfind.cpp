/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Scalar root finding. The oracles are functions whose roots are known in
 * closed form: polynomials with rational roots, a root sitting exactly on a
 * bracket endpoint, a double root (no sign change, so only Newton can see it),
 * and a transcendental equation whose root is checked by residual. Determinism
 * is asserted by running the same solve twice and requiring bit-identical
 * answers.
 */
#include <cmath>

#include "doctest.h"
#include "line/util/rootfind.h"

using line::Rational;
using line::Real50;
using line::RootResult;
using line::bracket_expand;
using line::root_bisect;
using line::root_brent;
using line::root_newton;
using line::root_newton_safe;

TEST_CASE("bisection finds a simple polynomial root") {
    // x^2 - 2 on [0,2]: root sqrt(2).
    auto f = [](double x) { return x * x - 2.0; };
    RootResult<double> r = root_bisect<double>(f, 0.0, 2.0, 1e-14);
    CHECK(r.converged);
    CHECK(r.root == doctest::Approx(std::sqrt(2.0)).epsilon(1e-13));
    CHECK(std::fabs(r.value) < 1e-13);
}

TEST_CASE("bisection returns a root sitting on a bracket endpoint") {
    auto f = [](double x) { return x * x - 4.0; };
    RootResult<double> lo = root_bisect<double>(f, 2.0, 5.0, 1e-12);
    CHECK(lo.converged);
    CHECK(lo.root == 2.0);
    CHECK(lo.value == 0.0);
    CHECK(lo.iterations == 0u);

    RootResult<double> hi = root_bisect<double>(f, -5.0, -2.0, 1e-12);
    CHECK(hi.converged);
    CHECK(hi.root == -2.0);
    CHECK(hi.value == 0.0);
}

TEST_CASE("bisection rejects a bracket without a sign change") {
    auto f = [](double x) { return x * x + 1.0; };
    CHECK_THROWS_AS(root_bisect<double>(f, -1.0, 1.0, 1e-12), line::InputError);
    // A double root is the interesting case: f >= 0 everywhere, no bracket.
    auto g = [](double x) { return (x - 1.0) * (x - 1.0); };
    CHECK_THROWS_AS(root_bisect<double>(g, 0.0, 3.0, 1e-12), line::InputError);
}

TEST_CASE("Brent finds the same roots as bisection, in fewer iterations") {
    auto f = [](double x) { return x * x * x - 2.0 * x - 5.0; };  // classic, root 2.0945514815...
    RootResult<double> b = root_bisect<double>(f, 2.0, 3.0, 1e-15);
    RootResult<double> r = root_brent<double>(f, 2.0, 3.0, 1e-15);
    CHECK(std::fabs(r.value) < 1e-12);
    CHECK(r.root == doctest::Approx(2.0945514815423265).epsilon(1e-12));
    CHECK(b.root == doctest::Approx(2.0945514815423265).epsilon(1e-12));
    CHECK(r.iterations < b.iterations);
}

TEST_CASE("Brent handles an endpoint root and a rejected bracket") {
    auto f = [](double x) { return x - 3.0; };
    RootResult<double> r = root_brent<double>(f, 3.0, 9.0, 1e-12);
    CHECK(r.root == 3.0);
    CHECK(r.value == 0.0);
    auto g = [](double x) { return std::exp(x) + 1.0; };
    CHECK_THROWS_AS(root_brent<double>(g, -1.0, 1.0, 1e-12), line::InputError);
}

TEST_CASE("Newton converges on a double root, where no bracket exists") {
    // f(x) = (x-1)^2, f'(x) = 2(x-1): Newton halves the error each step.
    auto f = [](double x) { return (x - 1.0) * (x - 1.0); };
    auto df = [](double x) { return 2.0 * (x - 1.0); };
    RootResult<double> r = root_newton<double>(f, df, 2.0, 1e-12);
    CHECK(r.converged);
    CHECK(r.root == doctest::Approx(1.0).epsilon(1e-11));
    CHECK(r.value < 1e-22);

    // A triple root, also invisible to bracketing on one side.
    auto f3 = [](double x) { return (x - 2.0) * (x - 2.0) * (x - 2.0); };
    auto df3 = [](double x) { return 3.0 * (x - 2.0) * (x - 2.0); };
    RootResult<double> r3 = root_newton<double>(f3, df3, 5.0, 1e-12);
    CHECK(r3.root == doctest::Approx(2.0).epsilon(1e-10));
}

TEST_CASE("Newton reports a vanishing derivative rather than dividing by zero") {
    auto f = [](double x) { return x * x + 1.0; };
    auto df = [](double x) { return 2.0 * x; };
    CHECK_THROWS_AS(root_newton<double>(f, df, 0.0, 1e-12), line::NumericError);
}

TEST_CASE("safeguarded Newton keeps the bracket when Newton would escape") {
    // f(x) = atan(x) - 0.5: plain Newton from a large start diverges, the
    // bracketed variant cannot leave [-10,10].
    auto f = [](double x) { return std::atan(x) - 0.5; };
    auto df = [](double x) { return 1.0 / (1.0 + x * x); };
    RootResult<double> r = root_newton_safe<double>(f, df, -10.0, 10.0, 1e-14);
    CHECK(r.converged);
    CHECK(r.root == doctest::Approx(std::tan(0.5)).epsilon(1e-12));
    CHECK(r.root > -10.0);
    CHECK(r.root < 10.0);

    RootResult<double> e = root_newton_safe<double>(f, df, std::tan(0.5), 10.0, 1e-14);
    CHECK(e.root == doctest::Approx(std::tan(0.5)).epsilon(1e-15));
    auto g = [](double x) { return x * x + 3.0; };
    auto dg = [](double x) { return 2.0 * x; };
    CHECK_THROWS_AS(root_newton_safe<double>(g, dg, -1.0, 1.0, 1e-12), line::InputError);
}

TEST_CASE("bracket_expand grows the upper end until the sign changes") {
    auto f = [](double x) { return x - 137.0; };
    double b = 1.0;
    bracket_expand<double>(f, 0.0, b);
    CHECK(b >= 137.0);
    RootResult<double> r = root_bisect<double>(f, 0.0, b, 1e-12);
    CHECK(r.root == doctest::Approx(137.0).epsilon(1e-13));

    auto g = [](double x) { return 1.0 + 0.0 * x; };  // never changes sign
    double b2 = 1.0;
    CHECK_THROWS_AS(bracket_expand<double>(g, 0.0, b2, 20), line::NumericError);
}

TEST_CASE("root finding is deterministic and backend independent") {
    auto f = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    RootResult<double> a = root_brent<double>(f, 2.0, 3.0, 1e-15);
    RootResult<double> b = root_brent<double>(f, 2.0, 3.0, 1e-15);
    CHECK(a.root == b.root);  // bit-identical, no global state
    CHECK(a.iterations == b.iterations);

    // Same problem at 50 digits: the double answer is a prefix of it.
    auto fr = [](Real50 x) { return x * x * x - Real50(2) * x - Real50(5); };
    RootResult<Real50> r = root_brent<Real50>(fr, Real50(2), Real50(3), Real50("1e-40"));
    CHECK(std::fabs(static_cast<double>(r.root) - a.root) < 1e-12);
    CHECK(static_cast<double>(abs(r.value)) < 1e-38);
}

TEST_CASE("bisection instantiates at exact rational arithmetic") {
    // 4 x^2 - 1 has the exact root 1/2, which bisection on [0,1] hits after a
    // single midpoint evaluation.
    // The `-> Rational` is belt and braces. It was load bearing until
    // `Rational` was pinned to et_off (line/num/number.h): under Boost < 1.72
    // a deduced return type here is an expression node over the Rational(4)
    // and Rational(1) temporaries, which die with the return statement, and
    // this case segfaulted every 20.04 worker on 2026-08-26.
    auto f = [](Rational x) -> Rational { return Rational(4) * x * x - Rational(1); };
    RootResult<Rational> r = root_bisect<Rational>(f, Rational(0), Rational(1), Rational(1, 1000));
    CHECK(r.converged);
    CHECK(r.root == Rational(1, 2));
    CHECK(r.value == Rational(0));
}
