/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Moment transforms. The oracles are the algebraic identities the transforms
 * are defined by: each pair is a mutual inverse, the tables satisfy their
 * defining recurrences and known values, and the round trips must be the
 * IDENTITY in exact arithmetic. That last point is the whole argument for the
 * exact backend here, since these alternating sums lose all significance in
 * double at modest order.
 */
#include <vector>

#include "doctest.h"
#include "line/api/moment/moment_binomial_from_factorial.h"
#include "line/api/moment/moment_binomial_from_negbinomial.h"
#include "line/api/moment/moment_binotrans.h"
#include "line/api/moment/moment_binotransinv.h"
#include "line/api/moment/moment_central_from_raw.h"
#include "line/api/moment/moment_factorial_from_binomial.h"
#include "line/api/moment/moment_factorial_from_raw.h"
#include "line/api/moment/moment_factorial_from_upfactorial.h"
#include "line/api/moment/moment_lah.h"
#include "line/api/moment/moment_negbinomial_from_binomial.h"
#include "line/api/moment/moment_negbinomial_from_upfactorial.h"
#include "line/api/moment/moment_raw_from_central.h"
#include "line/api/moment/moment_raw_from_factorial.h"
#include "line/api/moment/moment_raw_from_upfactorial.h"
#include "line/api/moment/moment_stirling1.h"
#include "line/api/moment/moment_stirling2.h"
#include "line/api/moment/moment_stirlingcycle.h"
#include "line/api/moment/moment_upfactorial_from_factorial.h"
#include "line/api/moment/moment_upfactorial_from_negbinomial.h"
#include "line/api/moment/moment_upfactorial_from_raw.h"

using line::Rational;
using namespace line::moment;

namespace {

/** Raw moments of the geometric-like sequence m_i = i! / 2^i, a valid moment set. */
std::vector<Rational> sample_raw(int n) {
    std::vector<Rational> m(n + 1);
    for (int i = 0; i <= n; ++i) m[i] = Rational(1) / line::num_pow_int(Rational(2), i);
    m[0] = Rational(1);
    return m;
}

}  // namespace

TEST_CASE("Stirling and Lah tables match their known values") {
    // Stirling second kind: S(4,2) = 7, S(5,3) = 25.
    auto S = moment_stirling2<Rational>(5);
    CHECK(S(4, 2) == Rational(7));
    CHECK(S(5, 3) == Rational(25));
    // Unsigned Stirling first kind (cycle): sigma(4,2) = 11, sigma(5,2) = 50.
    auto sig = moment_stirlingcycle<Rational>(5);
    CHECK(sig(4, 2) == Rational(11));
    CHECK(sig(5, 2) == Rational(50));
    // Signed first kind alternates: s(4,2) = +11, s(5,2) = -50.
    auto s1 = moment_stirling1<Rational>(5);
    CHECK(s1(4, 2) == Rational(11));
    CHECK(s1(5, 2) == Rational(-50));
    // Lah: L(n,k) = C(n-1,k-1) n!/k!, so L(4,2) = 36 and L(5,3) = 120.
    // Values cross-checked against MATLAB moment_lah(5).
    auto L = moment_lah<Rational>(5);
    CHECK(L(4, 2) == Rational(36));
    CHECK(L(5, 3) == Rational(120));
    CHECK(L(5, 5) == Rational(1));
    // Row sums of the second-kind table are the Bell numbers 1,1,2,5,15,52.
    const long bell[] = {1, 1, 2, 5, 15, 52};
    for (int i = 0; i <= 5; ++i) {
        Rational r(0);
        for (int j = 0; j <= i; ++j) r += S(i, j);
        CHECK(r == Rational(bell[i]));
    }
}

TEST_CASE("raw and central moments are exact mutual inverses") {
    const std::vector<Rational> m = sample_raw(6);
    const std::vector<Rational> mc = moment_central_from_raw(m);
    const std::vector<Rational> back = moment_raw_from_central(mc, m[1]);
    for (std::size_t i = 0; i < m.size(); ++i) CHECK(back[i] == m[i]);
    // The first central moment vanishes identically.
    CHECK(mc[1] == Rational(0));
    // The second central moment is the variance.
    CHECK(mc[2] == m[2] - m[1] * m[1]);
}

TEST_CASE("raw and factorial moments are exact mutual inverses") {
    const std::vector<Rational> m = sample_raw(6);
    const std::vector<Rational> f = moment_factorial_from_raw(m);
    const std::vector<Rational> back = moment_raw_from_factorial(f);
    for (std::size_t i = 0; i < m.size(); ++i) CHECK(back[i] == m[i]);
    // f_1 = m_1 always, since the first falling factorial is the identity.
    CHECK(f[1] == m[1]);
}

TEST_CASE("raw and rising-factorial moments are exact mutual inverses") {
    const std::vector<Rational> m = sample_raw(6);
    const std::vector<Rational> fp = moment_upfactorial_from_raw(m);
    const std::vector<Rational> back = moment_raw_from_upfactorial(fp);
    for (std::size_t i = 0; i < m.size(); ++i) CHECK(back[i] == m[i]);
}

TEST_CASE("falling and rising factorial moments are exact mutual inverses") {
    std::vector<Rational> f(6);
    f[0] = Rational(1);
    for (std::size_t i = 1; i < f.size(); ++i) f[i] = Rational(static_cast<long>(i), 3);
    const std::vector<Rational> fp = moment_upfactorial_from_factorial(f);
    const std::vector<Rational> back = moment_factorial_from_upfactorial(fp);
    for (std::size_t i = 0; i < f.size(); ++i) CHECK(back[i] == f[i]);
}

TEST_CASE("binomial and factorial moment scalings are exact mutual inverses") {
    std::vector<Rational> f(6);
    for (std::size_t i = 0; i < f.size(); ++i) f[i] = Rational(static_cast<long>(i) + 1, 5);
    const std::vector<Rational> b = moment_binomial_from_factorial(f);
    const std::vector<Rational> back = moment_factorial_from_binomial(b);
    for (std::size_t i = 0; i < f.size(); ++i) CHECK(back[i] == f[i]);
    CHECK(b[3] == f[3] / Rational(6));  // 3! = 6

    const std::vector<Rational> bm = moment_negbinomial_from_upfactorial(f);
    const std::vector<Rational> backfp = moment_upfactorial_from_negbinomial(bm);
    for (std::size_t i = 0; i < f.size(); ++i) CHECK(backfp[i] == f[i]);
}

TEST_CASE("binomial and negative-binomial moments are exact mutual inverses") {
    std::vector<Rational> bm(6);
    bm[0] = Rational(1);
    for (std::size_t i = 1; i < bm.size(); ++i) bm[i] = Rational(1, static_cast<long>(i) + 1);
    const std::vector<Rational> b = moment_binomial_from_negbinomial(bm);
    const std::vector<Rational> back = moment_negbinomial_from_binomial(b);
    for (std::size_t i = 0; i < bm.size(); ++i) CHECK(back[i] == bm[i]);
}

TEST_CASE("binomial transform is an exact involution") {
    std::vector<Rational> x(7);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = Rational(static_cast<long>(i) * i + 1, 7);
    const std::vector<Rational> y = moment_binotrans(x);
    const std::vector<Rational> back = moment_binotransinv(y);
    for (std::size_t i = 0; i < x.size(); ++i) CHECK(back[i] == x[i]);
}

TEST_CASE("exact arithmetic recovers central moments that double cannot") {
    // A shifted exponential, X = c + Y with c = 1e6 and Y unit-mean
    // exponential. Its raw moments m_i = sum_k C(i,k) c^(i-k) k! are dominated
    // by the shift, while the central moments are those of Y and are O(1): the
    // conversion is a difference of huge nearly-equal terms, which is where
    // double fails outright. The reference values and the size of the double
    // error were confirmed against MATLAB moment_central_from_raw.
    const int n = 6;
    const Rational c(1000000);

    std::vector<Rational> m(n + 1);
    for (int i = 0; i <= n; ++i) {
        Rational s(0);
        for (int k = 0; k <= i; ++k)
            s += line::num_nck<Rational>(i, k) * line::num_pow_int(c, static_cast<unsigned>(i - k)) *
                 line::num_factorial<Rational>(static_cast<unsigned>(k));
        m[i] = s;
    }

    // Central moments of the exponential: mc_i = sum_k (-1)^(i-k) C(i,k) k!.
    std::vector<Rational> expected(n + 1);
    for (int i = 0; i <= n; ++i) {
        Rational s(0);
        for (int k = 0; k <= i; ++k) {
            const Rational term = line::num_nck<Rational>(i, k) *
                                  line::num_factorial<Rational>(static_cast<unsigned>(k));
            s += ((i - k) % 2 == 0) ? term : -term;
        }
        expected[i] = s;
    }

    const std::vector<Rational> mc = moment_central_from_raw(m);
    for (int i = 0; i <= n; ++i) CHECK(mc[i] == expected[i]);  // exact, to the last digit

    // Same conversion in double: the third central moment comes out with the
    // wrong sign and the sixth is off by 18 orders of magnitude.
    std::vector<double> md(n + 1);
    for (int i = 0; i <= n; ++i) md[i] = static_cast<double>(m[i]);
    const std::vector<double> mcd = moment_central_from_raw(md);
    const double rel3 = std::fabs(mcd[3] - static_cast<double>(expected[3])) /
                        std::fabs(static_cast<double>(expected[3]));
    const double rel6 = std::fabs(mcd[6] - static_cast<double>(expected[6])) /
                        std::fabs(static_cast<double>(expected[6]));
    CHECK(rel3 > 1.0);
    CHECK(rel6 > 1e6);
}
