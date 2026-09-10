/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Throughput and queue-length bounds. The oracle for a bound is the exact
 * solution it brackets: pfqn_mva gives the true throughput, and every lower
 * bound must sit below it while every upper bound sits above. Checking that
 * relation in EXACT arithmetic is the point of the exercise, since a bound
 * violated only by rounding cannot be told apart from a real violation.
 */
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_qzgblow.h"
#include "line/api/pfqn/pfqn_qzgbup.h"
#include "line/api/pfqn/pfqn_ssd.h"
#include "line/api/pfqn/pfqn_xzabalow.h"
#include "line/api/pfqn/pfqn_xzabaup.h"
#include "line/api/pfqn/pfqn_xzgsblow.h"
#include "line/api/pfqn/pfqn_xzgsbup.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

namespace {

/** Single-class model: demands 1/2, 1/3, 1/5 with think time 1. */
template <class T>
std::vector<T> demands() {
    return {line::num_traits<T>::from_rational(1, 2), line::num_traits<T>::from_rational(1, 3),
            line::num_traits<T>::from_rational(1, 5)};
}

/** The same model as a (M x 1) matrix for pfqn_mva. */
template <class T>
Matrix<T> demand_matrix() {
    const std::vector<T> d = demands<T>();
    Matrix<T> L(d.size(), 1);
    for (std::size_t i = 0; i < d.size(); ++i) L(i, 0) = d[i];
    return L;
}

}  // namespace

TEST_CASE("ABA bounds bracket the exact throughput, exactly") {
    const std::vector<Rational> L = demands<Rational>();
    const Rational Z(1);
    for (int n = 1; n <= 6; ++n) {
        const Rational N(n);
        const Rational lo = pfqn_xzabalow(L, N, Z);
        const Rational hi = pfqn_xzabaup(L, N, Z);

        Matrix<Rational> Zm(1, 1);
        Zm(0, 0) = Z;
        auto mva = pfqn_mva(demand_matrix<Rational>(), std::vector<int>{n}, Zm);
        const Rational X = mva.XN[0];

        INFO("N = ", n);
        CHECK(lo <= X);   // exact comparison of rationals
        CHECK(X <= hi);
        CHECK(lo <= hi);
    }
}

TEST_CASE("ABA bounds match their closed forms") {
    const std::vector<Rational> L = demands<Rational>();
    const Rational Z(1), N(4);
    // sum(L) = 1/2 + 1/3 + 1/5 = 31/30; lower = N/(Z + N sum L) = 4/(1 + 4*31/30).
    CHECK(pfqn_xzabalow(L, N, Z) == Rational(4) / (Rational(1) + Rational(4) * Rational(31, 30)));
    // upper = min(1/Lmax, N/(sum L + Z)) = min(2, 4/(61/30)) = 120/61.
    CHECK(pfqn_xzabaup(L, N, Z) == Rational(120, 61));
}

TEST_CASE("geometric queue bounds bracket the exact queue length") {
    const std::vector<Rational> L = demands<Rational>();
    const Rational Z(1);
    const int n = 5;
    Matrix<Rational> Zm(1, 1);
    Zm(0, 0) = Z;
    auto mva = pfqn_mva(demand_matrix<Rational>(), std::vector<int>{n}, Zm);

    for (std::size_t i = 0; i < L.size(); ++i) {
        const Rational qlo = pfqn_qzgblow(L, Rational(n), Z, i);
        const Rational qhi = pfqn_qzgbup(L, Rational(n), Z, i);
        INFO("station ", i);
        CHECK(qlo <= qhi);
        CHECK(qlo <= mva.QN(i, 0));
        CHECK(mva.QN(i, 0) <= qhi + Rational(1));  // the upper bound is loose but must dominate
    }
}

TEST_CASE("geometric square-root bounds bracket the exact throughput") {
    const std::vector<double> L{0.5, 1.0 / 3.0, 0.2};
    const double Z = 1.0;
    for (int n = 2; n <= 6; ++n) {
        const double lo = pfqn_xzgsblow(L, static_cast<double>(n), Z);
        const double hi = pfqn_xzgsbup(L, static_cast<double>(n), Z);

        Matrix<double> Zm(1, 1);
        Zm(0, 0) = Z;
        auto mva = pfqn_mva(demand_matrix<double>(), std::vector<int>{n}, Zm);
        INFO("N = ", n);
        CHECK(lo <= mva.XN[0] * (1.0 + 1e-12));
        CHECK(mva.XN[0] <= hi * (1.0 + 1e-12));
        // The GSB pair must be at least as tight as the ABA pair.
        CHECK(lo >= pfqn_xzabalow(L, static_cast<double>(n), Z) * (1.0 - 1e-12));
    }
}

// No comma in this name on purpose: doctest treats commas in -tc as filter
// separators, so a comma here would make the case impossible to select.
TEST_CASE("pfqn_ssd brackets the exact throughput with a think time") {
    const std::vector<Rational> L = demands<Rational>();
    const Rational Z(1), N(5);

    // 2026-07-29: this case used to assert that Z>0 was REFUSED. That refusal
    // was correct about the bare +Z form, which is not a bound because the BJB
    // optimistic step needs sum_k Q_k(N-1) = N-1 and that fails once Z X(N-1)
    // jobs sit at the terminal: on this very model it gave Xhi = 1.4657980456,
    // BELOW the exact 1.6300060713. The reference has since scaled the queueing
    // term by the Lazowska Table 5.2 terminal-workload factor instead of
    // disaggregating the delay, which restores the bound rather than dodging
    // it. Observed MATLAB run today on L = [1/2, 1/3, 1/5], N = 5, Z = 1:
    // Xlo = 1.3480689823, Xhi = 1.8291463983, both bracketing 1.6300060713.
    auto z1 = pfqn_ssd(L, N, Z);
    CHECK(static_cast<double>(z1.Xlo) == doctest::Approx(1.3480689823).epsilon(1e-9));
    CHECK(static_cast<double>(z1.Xhi) == doctest::Approx(1.8291463983).epsilon(1e-9));
    CHECK(z1.Xlo <= z1.Xhi);

    Matrix<Rational> Z1m(1, 1);
    Z1m(0, 0) = Z;
    auto exact1 = pfqn_mva(demand_matrix<Rational>(), std::vector<int>{5}, Z1m);
    CHECK(z1.Xlo <= exact1.XN[0]);
    CHECK(exact1.XN[0] <= z1.Xhi);

    // A multiserver model takes the same corrected form.
    std::vector<Rational> two(L.size(), Rational(2));
    auto z1ms = pfqn_ssd(L, N, Z, two);
    CHECK(z1ms.Xlo <= z1ms.Xhi);

    // Z=0 must be unchanged: both correction factors collapse to 1.
    const Rational Z0(0);
    auto b = pfqn_ssd(L, N, Z0);
    CHECK(b.Xlo <= b.Xhi);
    Matrix<Rational> Zm(1, 1);
    Zm(0, 0) = Z0;
    auto mva = pfqn_mva(demand_matrix<Rational>(), std::vector<int>{5}, Zm);
    CHECK(b.Xlo <= mva.XN[0]);
    CHECK(mva.XN[0] <= b.Xhi);
}

TEST_CASE("bounds agree between double and exact arithmetic") {
    const std::vector<Rational> Lq = demands<Rational>();
    const std::vector<double> Ld{0.5, 1.0 / 3.0, 0.2};
    const Rational Zq(1), Nq(4);
    CHECK(static_cast<double>(pfqn_xzabalow(Lq, Nq, Zq)) ==
          doctest::Approx(pfqn_xzabalow(Ld, 4.0, 1.0)).epsilon(1e-12));
    CHECK(static_cast<double>(pfqn_qzgblow(Lq, Nq, Zq, 0)) ==
          doctest::Approx(pfqn_qzgblow(Ld, 4.0, 1.0, 0)).epsilon(1e-12));
    // pfqn_ssd took Z=0 only until 2026-07-29, when the Lazowska Table 5.2
    // correction gave it a think time. Z=0 is kept because it is the Theorem 5
    // form that the correction must reduce to, and Z=1 is added because that is
    // the branch the correction actually introduced: it is +, * and / only, so
    // it must agree across arithmetics just as tightly.
    CHECK(static_cast<double>(pfqn_ssd(Lq, Nq, Rational(0)).Xlo) ==
          doctest::Approx(pfqn_ssd(Ld, 4.0, 0.0).Xlo).epsilon(1e-12));
    CHECK(static_cast<double>(pfqn_ssd(Lq, Nq, Zq).Xlo) ==
          doctest::Approx(pfqn_ssd(Ld, 4.0, 1.0).Xlo).epsilon(1e-12));
    CHECK(static_cast<double>(pfqn_ssd(Lq, Nq, Zq).Xhi) ==
          doctest::Approx(pfqn_ssd(Ld, 4.0, 1.0).Xhi).epsilon(1e-12));
}
