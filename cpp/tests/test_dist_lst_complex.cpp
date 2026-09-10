/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Laplace-Stieltjes transform off the real axis.
 *
 * A transform is evaluated at a COMPLEX argument by everything that inverts it
 * or locates its roots, so `dist_lst` carries a complex overload beside the real
 * one, the C++ counterpart of the JAR widening sn.lst to
 * SerializableFunction<Complex, Complex>.
 *
 * THE ORACLES ARE THE CLOSED FORMS THEMSELVES, continued to the complex plane:
 * mu/(mu+s) for the exponential, (r/(r+s))^k for the Erlang, exp(-sd) for the
 * deterministic law, and (e^-sa - e^-sb)/(s(b-a)) for the uniform. Each is also
 * checked to agree with the REAL overload on the real axis, which is what makes
 * the two a single function rather than two implementations that happen to
 * share a name.
 */
#include <complex>

#include "doctest.h"
#include "line/lang/distribution.h"

namespace lang = line::lang;
using Dist = line::lang::Distrib<double>;
using cdbl = std::complex<double>;

TEST_CASE("dist_lst complex: exponential matches its closed form") {
    const double mu = 1.4;
    const Dist d = Dist::exp_rate(mu);
    const cdbl s(0.5, 1.0);
    const cdbl ref = cdbl(mu, 0.0) / (s + mu);
    const cdbl got = lang::dist_lst(d, s);
    CHECK(got.real() == doctest::Approx(ref.real()).epsilon(1e-12));
    CHECK(got.imag() == doctest::Approx(ref.imag()).epsilon(1e-12));
    // the two overloads agree on the real axis
    CHECK(lang::dist_lst(d, cdbl(0.5, 0.0)).real() ==
          doctest::Approx(lang::dist_lst(d, 0.5)).epsilon(1e-12));
}

TEST_CASE("dist_lst complex: Erlang matches its closed form") {
    const Dist d = Dist::erlang(2.0, 2);
    const cdbl s(0.3, -0.7);
    const cdbl one = cdbl(2.0, 0.0) / (s + 2.0);
    const cdbl ref = one * one;
    const cdbl got = lang::dist_lst(d, s);
    CHECK(got.real() == doctest::Approx(ref.real()).epsilon(1e-12));
    CHECK(got.imag() == doctest::Approx(ref.imag()).epsilon(1e-12));
}

TEST_CASE("dist_lst complex: deterministic, where a CDF sum would see nothing") {
    const double v = 0.75;
    const Dist d = Dist::det(v);
    const cdbl s(0.4, 2.0);
    const cdbl ref = std::exp(-s * v);
    const cdbl got = lang::dist_lst(d, s);
    CHECK(got.real() == doctest::Approx(ref.real()).epsilon(1e-12));
    CHECK(got.imag() == doctest::Approx(ref.imag()).epsilon(1e-12));
}

TEST_CASE("dist_lst complex: uniform matches its closed form") {
    const double a = 0.1, b = 0.9;
    const Dist d = Dist::uniform(a, b);
    const cdbl s(0.6, 0.8);
    const cdbl ref = (std::exp(-s * a) - std::exp(-s * b)) / (s * (b - a));
    const cdbl got = lang::dist_lst(d, s);
    CHECK(got.real() == doctest::Approx(ref.real()).epsilon(1e-12));
    CHECK(got.imag() == doctest::Approx(ref.imag()).epsilon(1e-12));
    CHECK(lang::dist_lst(d, cdbl(0.0, 0.0)).real() == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("dist_lst complex: a law with no complex closed form still transforms") {
    // Gamma has one here; the point of the check is that the complex and real
    // overloads answer the same on the real axis for a non-phase-type family.
    const Dist d = Dist::gamma_dist(2.0, 0.5);
    CHECK(lang::dist_lst(d, cdbl(0.5, 0.0)).real() ==
          doctest::Approx(lang::dist_lst(d, 0.5)).epsilon(1e-9));
}
