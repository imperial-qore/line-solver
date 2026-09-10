/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The second family of pfqn bounds: the convolutional and performance bound
 * hierarchies, the successively improving bounds, the load-dependent BCMP
 * bound, the multiclass composite bound, the robust box bounds and the
 * square-root non-iterative approximation.
 *
 * The oracle for a bound is the exact solution it brackets, so every bound is
 * checked against pfqn_mva on several models. Where a bound does NOT bracket,
 * the violation is pinned with the numbers that produce it rather than dropped,
 * following the pattern pfqn_ssd already sets in test_pfqn_bounds.cpp.
 *
 * Most of these are finite rational computations, so the bracket is checked in
 * EXACT arithmetic wherever the bound admits it: a bound violated only by
 * rounding cannot be told apart from a real violation.
 *
 * Models:
 *   B  single class, L = [1/2, 1/3, 1/5], N = 5, Z = 1 (and Z = 0)
 *   A  M = 3, R = 2, L = [1 .5; .7 1.2; .3 .9], N = [2 3], Z = [1 .5]
 *   C  repairman, one station, L = [.6 .4], N = [3 2], Z = [1 2]
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_cbh.h"
#include "line/api/pfqn/pfqn_ldbcmp.h"
#include "line/api/pfqn/pfqn_mcub.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mwrbb.h"
#include "line/api/pfqn/pfqn_pbh.h"
#include "line/api/pfqn/pfqn_scb.h"
#include "line/api/pfqn/pfqn_sib.h"
#include "line/api/pfqn/pfqn_sqni.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

namespace {

template <class T>
std::vector<T> demandsB() {
    return {line::num_traits<T>::from_rational(1, 2), line::num_traits<T>::from_rational(1, 3),
            line::num_traits<T>::from_rational(1, 5)};
}

template <class T>
Matrix<T> matrixB() {
    const std::vector<T> d = demandsB<T>();
    Matrix<T> L(3, 1);
    for (std::size_t i = 0; i < 3; ++i) L(i, 0) = d[i];
    return L;
}

/** Exact throughput of model B at population n with think time z. */
template <class T>
T exactB(int n, const T& z) {
    Matrix<T> Zm(1, 1);
    Zm(0, 0) = z;
    return pfqn_mva(matrixB<T>(), std::vector<int>{n}, Zm).XN[0];
}

Matrix<double> modelA() { return Matrix<double>{{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}}; }

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_cbh
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_cbh is exact at level M-1 and brackets elsewhere, exactly") {
    const std::vector<Rational> L = demandsB<Rational>();
    const Rational Z(1);
    // MATLAB pfqn_cbh([1/2;1/3;1/5],5,1,2) and level 3 both return
    // 1.630006071256988 for BOTH sides, which is exactly pfqn_mva's value:
    // once c = max(1, M-level) collapses to one, the filled column is the
    // single-server column, which is itself exact, and every remaining station
    // is convolved exactly. The hierarchy is exact from level M-1 upwards.
    const Rational X = exactB<Rational>(5, Z);
    for (int level = 2; level <= 3; ++level) {
        const CbhBounds<Rational> b = pfqn_cbh(L, 5, Z, level);
        INFO("level ", level);
        CHECK(b.Xlo == X);  // exact equality of rationals
        CHECK(b.Xhi == X);
    }
    CHECK(static_cast<double>(X) == doctest::Approx(1.630006071256988).epsilon(1e-12));

    // Level 1 fills two of the three stations from the balanced estimate and is
    // a genuine bracket. MATLAB at Z = 0: pfqn_cbh(L,5,0,1) -> [1.728357955409774,
    // 1.934090524217816] with pfqn_mva 1.855757234273755.
    const Rational Z0(0);
    const CbhBounds<Rational> b1 = pfqn_cbh(L, 5, Z0, 1);
    const Rational X0 = exactB<Rational>(5, Z0);
    CHECK(b1.Xlo <= X0);
    CHECK(X0 <= b1.Xhi);
    CHECK(static_cast<double>(b1.Xlo) == doctest::Approx(1.728357955409774).epsilon(1e-12));
    CHECK(static_cast<double>(b1.Xhi) == doctest::Approx(1.934090524217816).epsilon(1e-12));
    CHECK(static_cast<double>(X0) == doctest::Approx(1.855757234273755).epsilon(1e-12));
}

TEST_CASE("pfqn_cbh brackets across populations, in exact arithmetic") {
    const std::vector<Rational> L = demandsB<Rational>();
    for (int n = 1; n <= 8; ++n) {
        for (int zi = 0; zi <= 1; ++zi) {
            const Rational Z(zi);
            const CbhBounds<Rational> b = pfqn_cbh(L, n, Z, 1);
            const Rational X = exactB<Rational>(n, Z);
            INFO("N = ", n, " Z = ", zi);
            CHECK(b.Xlo <= X);
            CHECK(X <= b.Xhi);
        }
    }
}

// ---------------------------------------------------------------------------
// pfqn_pbh, pfqn_pbk, pfqn_bjbk
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_pbh matches MATLAB and brackets the exact throughput") {
    const std::vector<Rational> L = demandsB<Rational>();
    const Rational Z(1);
    const PbhBounds<Rational> b = pfqn_pbh(L, 5, Z, 2);
    // MATLAB pfqn_pbh([1/2;1/3;1/5],5,1,2) ->
    //   Xlo 1.536694147857819, Xhi 1.678336718843350,
    //   Qlo [1.674215473776748 0.9148395603294232 0.4522777735472179]
    CHECK(static_cast<double>(b.Xlo) == doctest::Approx(1.536694147857819).epsilon(1e-12));
    CHECK(static_cast<double>(b.Xhi) == doctest::Approx(1.678336718843350).epsilon(1e-12));
    CHECK(static_cast<double>(b.Qlo[0]) == doctest::Approx(1.674215473776748).epsilon(1e-12));
    CHECK(static_cast<double>(b.Qlo[1]) == doctest::Approx(0.9148395603294232).epsilon(1e-12));
    CHECK(static_cast<double>(b.Qlo[2]) == doctest::Approx(0.4522777735472179).epsilon(1e-12));
    const Rational X = exactB<Rational>(5, Z);
    CHECK(b.Xlo <= X);
    CHECK(X <= b.Xhi);
}

TEST_CASE("the PBH hierarchy tightens with the level and stays a bracket") {
    const std::vector<Rational> L = demandsB<Rational>();
    for (int zi = 0; zi <= 1; ++zi) {
        const Rational Z(zi);
        const Rational X = exactB<Rational>(6, Z);
        Rational prevLo(0), prevHi(0);
        for (int level = 1; level <= 6; ++level) {
            const PbhBounds<Rational> b = pfqn_pbh(L, 6, Z, level);
            INFO("Z = ", zi, " level = ", level);
            CHECK(b.Xlo <= X);
            CHECK(X <= b.Xhi);
            if (level > 1) {
                CHECK(b.Xlo >= prevLo);  // nested, monotone in the level
                CHECK(b.Xhi <= prevHi);
            }
            prevLo = b.Xlo;
            prevHi = b.Xhi;
        }
        // At level N the hierarchy is exact MVA.
        CHECK(prevLo == X);
        CHECK(prevHi == X);
    }
}

TEST_CASE("pfqn_pbk and pfqn_bjbk are the same function as pfqn_pbh") {
    // REFERENCE REDUNDANCY, pinned rather than hidden: pfqn_pbk.m and
    // pfqn_bjbk.m are one-line forwarders to pfqn_pbh with identical
    // arguments, so PB(k) and BJB(k) return the same numbers for every input.
    const std::vector<Rational> L = demandsB<Rational>();
    const Rational Z0(0);
    for (int k = 1; k <= 4; ++k) {
        const PbhBounds<Rational> a = pfqn_pbk(L, 5, Z0, k);
        const PbhBounds<Rational> b = pfqn_bjbk(L, 5, Z0, k);
        const PbhBounds<Rational> c = pfqn_pbh(L, 5, Z0, k);
        INFO("k = ", k);
        CHECK(a.Xlo == c.Xlo);
        CHECK(a.Xhi == c.Xhi);
        CHECK(b.Xlo == c.Xlo);
        CHECK(b.Xhi == c.Xhi);
    }
    // MATLAB pfqn_bjbk(L,5,0,2) -> [1.753846153846154 1.933471933471933]
    const PbhBounds<Rational> b2 = pfqn_bjbk(L, 5, Z0, 2);
    CHECK(static_cast<double>(b2.Xlo) == doctest::Approx(1.753846153846154).epsilon(1e-12));
    CHECK(static_cast<double>(b2.Xhi) == doctest::Approx(1.933471933471933).epsilon(1e-12));
    // MATLAB pfqn_pbk(L,5,0,3) -> [1.810619172575512 1.878446268839603]
    const PbhBounds<Rational> b3 = pfqn_pbk(L, 5, Z0, 3);
    CHECK(static_cast<double>(b3.Xlo) == doctest::Approx(1.810619172575512).epsilon(1e-12));
    CHECK(static_cast<double>(b3.Xhi) == doctest::Approx(1.878446268839603).epsilon(1e-12));
    // MATLAB pfqn_pbh(L,5,0,1) -> [1.648351648351648 2]
    const PbhBounds<Rational> b1 = pfqn_pbh(L, 5, Z0, 1);
    CHECK(static_cast<double>(b1.Xlo) == doctest::Approx(1.648351648351648).epsilon(1e-12));
    CHECK(b1.Xhi == Rational(2));  // 1/Lmax, exactly
}

// ---------------------------------------------------------------------------
// pfqn_sib
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_sib matches MATLAB and brackets the exact throughput at Z = 0") {
    const std::vector<double> L = demandsB<double>();
    const double X = static_cast<double>(exactB<Rational>(5, Rational(0)));
    CHECK(X == doctest::Approx(1.855757234273755).epsilon(1e-12));
    // MATLAB pfqn_sib([1/2;1/3;1/5],5,0,1) -> [1.760695427409565
    // 1.933471933471933], which the port reproduces exactly: at level 1 the
    // Theorem-3.5 correction term is empty and the two agree.
    const SibBounds<double> b1 = pfqn_sib(L, 5, 0.0, 1);
    CHECK(b1.Xlo == doctest::Approx(1.760695427409565).epsilon(1e-12));
    CHECK(b1.Xhi == doctest::Approx(1.933471933471933).epsilon(1e-12));
    for (int lv = 1; lv <= 4; ++lv) {
        const SibBounds<double> b = pfqn_sib(L, 5, 0.0, lv);
        INFO("level ", lv);
        CHECK(b.Xlo <= X * (1.0 + 1e-12));
        CHECK(X <= b.Xhi * (1.0 + 1e-12));
        CHECK(b.Wlo <= b.Whi);
    }
    // DELIBERATE DIVERGENCE FROM MATLAB ABOVE LEVEL 1, pinned here with both
    // numbers. MATLAB's nested phi_u1 overwrites the parent's `eta` (see the
    // note in pfqn_sib.h), which inflates the Theorem-3.5 upper bound on phi by
    // the ratio of the two etas, 3/4 over 1/2. The inflated value then loses to
    // the Section-2 baseline in the `min`, so MATLAB returns the baseline at
    // every level >= 2: [1.740695097407274 1.933471933471933] for levels 2, 3
    // and 4 alike -- a LOOSER bracket than its own level 1, which is not what a
    // hierarchy is supposed to do. The port computes the intended value, so its
    // lower bound tightens monotonically instead. If pfqn_sib.m is ever fixed
    // (declare the helpers as subfunctions, or rename the local eta), these
    // become equalities and this comment should go.
    const SibBounds<double> b2 = pfqn_sib(L, 5, 0.0, 2);
    const SibBounds<double> b3 = pfqn_sib(L, 5, 0.0, 3);
    const SibBounds<double> b4 = pfqn_sib(L, 5, 0.0, 4);
    CHECK(b2.Xlo == doctest::Approx(1.818193451777910).epsilon(1e-12));
    CHECK(b3.Xlo == doctest::Approx(1.836935306200480).epsilon(1e-12));
    CHECK(b4.Xlo == doctest::Approx(1.839965957477760).epsilon(1e-12));
    CHECK(b2.Xlo > b1.Xlo);   // monotone, unlike the reference
    CHECK(b3.Xlo > b2.Xlo);
    CHECK(b4.Xlo > b3.Xlo);
    CHECK(b4.Xlo < X);        // and still a lower bound
    // MATLAB returns its baseline 1.740695097407274 at all three levels, which
    // is below the port's level-1 value; both are valid bounds, the port's are
    // tighter.
    CHECK(1.740695097407274 < b1.Xlo);
}

TEST_CASE("pfqn_sib brackets across populations and rejects a delay") {
    const std::vector<double> L = demandsB<double>();
    for (int n = 2; n <= 10; ++n) {
        const double X = static_cast<double>(exactB<Rational>(n, Rational(0)));
        const SibBounds<double> b = pfqn_sib(L, n, 0.0, 3);
        INFO("N = ", n);
        CHECK(b.Xlo <= X * (1.0 + 1e-12));
        CHECK(X <= b.Xhi * (1.0 + 1e-12));
    }
    // The reference refuses a delay rather than return an invalid bracket.
    CHECK_THROWS_AS(pfqn_sib(L, 5, 1.0, 3), line::InputError);
    CHECK_THROWS_AS(pfqn_sib(L, 1, 0.0, 3), line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_ldbcmp
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_ldbcmp matches MATLAB; its bound is valid but very loose here") {
    const std::vector<Rational> L = demandsB<Rational>();
    const LdBcmpBound<Rational> b = pfqn_ldbcmp(L, Rational(5), Rational(1));
    // MATLAB pfqn_ldbcmp([1/2;1/3;1/5],5,1) ->
    //   Xlo 0.5003678303552243, Rhi 9.992648800883879, Qhat 4.666666666666666
    CHECK(b.applicable);
    CHECK(static_cast<double>(b.Xlo) == doctest::Approx(0.5003678303552243).epsilon(1e-9));
    CHECK(static_cast<double>(b.Rhi) == doctest::Approx(9.992648800883879).epsilon(1e-9));
    CHECK(b.Qhat == Rational(14, 3));
    // It is a LOWER bound on the throughput, and it holds -- but at N = 5 the
    // asymptotic closed-open equivalence is far from its regime, so the bound
    // sits at 31% of the exact 1.630006071256988.
    const Rational X = exactB<Rational>(5, Rational(1));
    CHECK(b.Xlo <= X);
    CHECK(static_cast<double>(b.Xlo / X) < 0.35);
}

TEST_CASE("pfqn_ldbcmp is asymptotically exact and reports inapplicability") {
    const std::vector<Rational> L = demandsB<Rational>();
    // As N grows the bound approaches the exact throughput from below. The
    // fixed point raises (Dmax X)^N, whose exact rational form has a
    // denominator that squares on every iteration, so the sweep runs in
    // double: the exactness of the bound is established at N = 5 above, and
    // what is being checked here is asymptotic behaviour, not arithmetic.
    const std::vector<double> Ld = demandsB<double>();
    double prevRatio = 0.0;
    for (int n = 10; n <= 40; n += 10) {
        const LdBcmpBound<double> b = pfqn_ldbcmp(Ld, static_cast<double>(n), 1.0);
        const double X = static_cast<double>(exactB<Rational>(n, Rational(1)));
        INFO("N = ", n);
        REQUIRE(b.applicable);
        CHECK(b.Xlo <= X);
        const double ratio = b.Xlo / X;
        CHECK(ratio > prevRatio);  // monotone approach to the exact value
        prevRatio = ratio;
    }
    CHECK(prevRatio > 0.85);
    // N < Qhat is outside the regime; MATLAB returns NaN, the port flags.
    const LdBcmpBound<Rational> bad = pfqn_ldbcmp(L, Rational(1), Rational(1));
    CHECK_FALSE(bad.applicable);
}

// ---------------------------------------------------------------------------
// pfqn_mcub
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_mcub matches MATLAB but its upper bound is NOT a bound here") {
    Matrix<double> L = modelA();
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 0.5};
    const McubBounds<double> b = pfqn_mcub(L, N, Z);
    // MATLAB pfqn_mcub([1 .5;.7 1.2;.3 .9],[2 3],[1 .5]) ->
    //   Xub [0.7775768535262206 0.6666666666666667]
    //   Xlb [0.2857142857142857 0.3797468354430379]
    CHECK(b.Xub[0] == doctest::Approx(0.7775768535262206).epsilon(1e-12));
    CHECK(b.Xub[1] == doctest::Approx(0.6666666666666667).epsilon(1e-12));
    CHECK(b.Xlb[0] == doctest::Approx(0.2857142857142857).epsilon(1e-12));
    CHECK(b.Xlb[1] == doctest::Approx(0.3797468354430379).epsilon(1e-12));

    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 0.5;
    const auto mva = pfqn_mva(L, std::vector<int>{2, 3}, Zm);
    for (std::size_t r = 0; r < 2; ++r) {
        INFO("class ", r);
        CHECK(b.Xlb[r] <= mva.XN[r] * (1.0 + 1e-12));  // the BJB lower bound holds
        CHECK(mva.XN[r] <= b.Xub[r] * (1.0 + 1e-12));  // and so does the composite upper
    }

    // The bracket is exact in rational arithmetic, since eqs. (10) and (13)-(16)
    // are finite rational expressions.
    Matrix<Rational> Lq{{Rational(1), Rational(1, 2)},
                        {Rational(7, 10), Rational(6, 5)},
                        {Rational(3, 10), Rational(9, 10)}};
    Matrix<Rational> Zq(1, 2);
    Zq(0, 0) = Rational(1);
    Zq(0, 1) = Rational(1, 2);
    const McubBounds<Rational> bq =
        pfqn_mcub(Lq, std::vector<Rational>{Rational(2), Rational(3)},
                  std::vector<Rational>{Rational(1), Rational(1, 2)});
    const auto mvaq = pfqn_mva(Lq, std::vector<int>{2, 3}, Zq);
    for (std::size_t r = 0; r < 2; ++r) {
        INFO("class ", r);
        CHECK(bq.Xlb[r] <= mvaq.XN[r]);
        CHECK(mvaq.XN[r] <= bq.Xub[r]);
    }
    CHECK(bq.Xub[1] == Rational(2, 3));  // 1/Lmax for class 2, exactly
}

TEST_CASE("pfqn_mcub tightens relative to the per-class asymptotic bound") {
    // The composite upper bound is meant to beat 1/max_k L_kr at moderate load.
    Matrix<Rational> Lq{{Rational(1), Rational(1, 2)},
                        {Rational(7, 10), Rational(6, 5)},
                        {Rational(3, 10), Rational(9, 10)}};
    const McubBounds<Rational> b =
        pfqn_mcub(Lq, std::vector<Rational>{Rational(4), Rational(4)},
                  std::vector<Rational>{Rational(0), Rational(0)});
    // Class 1 bottleneck demand is 1, so the asymptotic bound is 1.
    CHECK(b.Xub[0] < Rational(1));
    // Class 2 bottleneck demand is 6/5, so the asymptotic bound is 5/6.
    CHECK(b.Xub[1] < Rational(5, 6));
    Matrix<Rational> Zq;
    const auto mvaq = pfqn_mva(Lq, std::vector<int>{4, 4}, Zq);
    for (std::size_t r = 0; r < 2; ++r) {
        INFO("class ", r);
        CHECK(b.Xlb[r] <= mvaq.XN[r]);
        CHECK(mvaq.XN[r] <= b.Xub[r]);
    }
}

// ---------------------------------------------------------------------------
// pfqn_mwrbb
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_mwrbb matches MATLAB and brackets the exact product-form solution") {
    Matrix<double> V{{1.0, 1.0}, {1.0, 1.0}};
    Matrix<double> S{{0.4, 0.2}, {0.1, 0.5}};
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 1.0};
    const MwrbbBounds<double> b = pfqn_mwrbb(V, S, N, Z);
    // MATLAB pfqn_mwrbb(V,S,[2 3],[1 1]) ->
    //   Xlo [0.4878048780487804 0.7317073170731706]
    //   Xup [1.3333333333333330 1.7647058823529410]
    CHECK(b.Xlo[0] == doctest::Approx(0.4878048780487804).epsilon(1e-9));
    CHECK(b.Xlo[1] == doctest::Approx(0.7317073170731706).epsilon(1e-9));
    CHECK(b.Xup[0] == doctest::Approx(1.3333333333333330).epsilon(1e-9));
    CHECK(b.Xup[1] == doctest::Approx(1.7647058823529410).epsilon(1e-9));

    // The bounds are distribution-insensitive, so they must contain the
    // exponential product-form solution as a special case.
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 1.0;
    const auto mva = pfqn_mva(S, std::vector<int>{2, 3}, Zm);
    for (std::size_t c = 0; c < 2; ++c) {
        INFO("class ", c);
        CHECK(b.Xlo[c] <= mva.XN[c] * (1.0 + 1e-9));
        CHECK(mva.XN[c] <= b.Xup[c] * (1.0 + 1e-9));
    }

    // Processor sharing (Lemma 2) tightens the lower bound relative to FIFO.
    const std::vector<MwrbbSched> ps(2, MwrbbSched::Ps);
    const MwrbbBounds<double> bps = pfqn_mwrbb(V, S, N, Z, ps, std::vector<int>());
    // MATLAB with sched = [1;1] -> Xlo [0.5714285714285714 0.7428571428571364]
    CHECK(bps.Xlo[0] == doctest::Approx(0.5714285714285714).epsilon(1e-9));
    CHECK(bps.Xlo[1] == doctest::Approx(0.7428571428571364).epsilon(1e-8));
    CHECK(bps.Xlo[0] >= b.Xlo[0]);
    CHECK(bps.Xlo[1] >= b.Xlo[1]);
    for (std::size_t c = 0; c < 2; ++c) {
        INFO("PS class ", c);
        CHECK(bps.Xlo[c] <= mva.XN[c] * (1.0 + 1e-9));
    }
}

TEST_CASE("pfqn_mwrbb brackets under every discipline code") {
    Matrix<double> V{{1.0, 1.0}, {1.0, 1.0}};
    Matrix<double> S{{0.4, 0.2}, {0.1, 0.5}};
    const std::vector<double> N{2.0, 3.0}, Z{1.0, 1.0};
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 1.0;
    const auto mva = pfqn_mva(S, std::vector<int>{2, 3}, Zm);
    const MwrbbSched codes[5] = {MwrbbSched::Fifo, MwrbbSched::Ps, MwrbbSched::PrioNonPreemptive,
                                 MwrbbSched::PrioPreemptive, MwrbbSched::Aba};
    for (int k = 0; k < 5; ++k) {
        const std::vector<MwrbbSched> sched(2, codes[k]);
        const std::vector<int> prio{0, 1};
        const MwrbbBounds<double> b = pfqn_mwrbb(V, S, N, Z, sched, prio);
        for (std::size_t c = 0; c < 2; ++c) {
            INFO("discipline ", k, " class ", c);
            CHECK(b.Xlo[c] <= b.Xup[c]);
            // The exponential product-form solution is a FIFO/PS instance; the
            // priority codes describe a different station, so only the
            // discipline-independent upper bound must dominate it there.
            CHECK(mva.XN[c] <= b.Xup[c] * (1.0 + 1e-9));
            if (codes[k] == MwrbbSched::Fifo || codes[k] == MwrbbSched::Ps ||
                codes[k] == MwrbbSched::Aba)
                CHECK(b.Xlo[c] <= mva.XN[c] * (1.0 + 1e-9));
        }
    }
}

// ---------------------------------------------------------------------------
// pfqn_sqni
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_sqni matches MATLAB on the repairman model C") {
    const std::vector<double> N{3.0, 2.0}, L{0.6, 0.4}, Z{1.0, 2.0};
    const SqniResult<double> r = pfqn_sqni(N, L, Z);
    // MATLAB pfqn_sqni([3 2],[.6 .4],[1 2]) ->
    //   X = [1.163105284445113 0.6448191940253879]
    //   Q = [1.836894715554887 0.7103616119492242]
    CHECK(r.X[0] == doctest::Approx(1.163105284445113).epsilon(1e-12));
    CHECK(r.X[1] == doctest::Approx(0.6448191940253879).epsilon(1e-12));
    CHECK(r.Q[0] == doctest::Approx(1.836894715554887).epsilon(1e-12));
    CHECK(r.Q[1] == doctest::Approx(0.7103616119492242).epsilon(1e-12));

    // Accuracy against exact MVA on the same model. The approximation is not a
    // bound: class 1 is over-predicted and class 2 under-predicted here, and
    // the station is saturated (U > 1), which is what makes this a hard case.
    Matrix<double> Lm(1, 2);
    Lm(0, 0) = 0.6;
    Lm(0, 1) = 0.4;
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 1.0;
    Zm(0, 1) = 2.0;
    const auto mva = pfqn_mva(Lm, std::vector<int>{3, 2}, Zm);
    const double e0 = std::fabs(r.X[0] - mva.XN[0]) / mva.XN[0];
    const double e1 = std::fabs(r.X[1] - mva.XN[1]) / mva.XN[1];
    INFO("relative errors ", e0, " ", e1);
    CHECK(e0 < 0.75);
    CHECK(e1 < 0.75);
    // Little's law must hold identically for every class, by construction:
    // Q_r = N_r - X_r Z_r.
    for (std::size_t c = 0; c < 2; ++c) {
        INFO("class ", c);
        CHECK(r.Q[c] + r.X[c] * Z[c] == doctest::Approx(N[c]).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_sqni handles the unit population and zero-think-time classes") {
    // N = 1 total: the closed form is exact, X_r = N_r/(Z_r + L_r).
    const SqniResult<double> r =
        pfqn_sqni(std::vector<double>{1.0, 0.0}, std::vector<double>{0.6, 0.4},
                  std::vector<double>{1.0, 2.0});
    CHECK(r.X[0] == doctest::Approx(1.0 / 1.6).epsilon(1e-12));
    CHECK(r.X[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.Q[0] == doctest::Approx(0.6 / 1.6).epsilon(1e-12));

    // A class with Z = 0 is filled in from the total queue length after the
    // main loop, so its queue length is its whole population.
    const SqniResult<double> r2 =
        pfqn_sqni(std::vector<double>{3.0, 2.0}, std::vector<double>{0.6, 0.4},
                  std::vector<double>{1.0, 0.0});
    CHECK(r2.Q[1] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r2.X[1] > 0.0);
}

// ---- SCB: Dowdy et al. (1992) single-class bounds of multiclass networks ----
//
// The oracle is the paper itself: Table I (p.200) for the aggregation gap,
// Section 4.7 for the class-count bound, and the Section-2 example plus the
// eight Table-II rows for the bracket. A row is bracketed by solving the
// MULTICLASS model exactly and checking that the exactly-weighted single-class
// aggregate's [Xlo, Xhi] contains it -- the bracketed object is the multiclass
// system, not the single-class model, which is what makes scb unlike every
// other bound in this file.

TEST_CASE("pfqn_scbgap reproduces Table I and the 50% limit") {
    // Table I, rows N = 1..5, columns K = 1..5, in percent.
    const int pub[5][5] = {{0, 0, 0, 0, 0},    {0, 33, 33, 33, 33}, {0, 25, 40, 40, 40},
                           {0, 20, 33, 43, 43}, {0, 17, 29, 38, 44}};
    for (long N = 1; N <= 5; ++N)
        for (long K = 1; K <= 5; ++K) {
            INFO("N=", N, " K=", K);
            CHECK(static_cast<int>(std::lround(100 * pfqn_scbgap<double>(N, K))) == pub[N - 1][K - 1]);
        }
    // The bound is below 50% everywhere and approaches it as N = K -> inf.
    CHECK(pfqn_scbgap<double>(1000, 1000) < 0.5);
    CHECK(pfqn_scbgap<double>(1000, 1000) > 0.499);
    // Exact arithmetic: the gap is a ratio of integers, so it is representable.
    CHECK(pfqn_scbgap<Rational>(4, 3) == Rational(1, 3));
    // Expression (4) at r < N, and the Theorem-5 form agreeing with (3) at r=N=K.
    CHECK(pfqn_scbgap<double>(8, 5, 2, false) == doctest::Approx(1.0 / 3.0));
    CHECK(pfqn_scbgap<double>(5, 5, 5, true) == doctest::Approx(pfqn_scbgap<double>(5, 5, 5, false)));
    // The undominated form is defined for r <= K, is strictly tighter below it
    // and coincides at r = K; past K it stops being a bound, so it is refused.
    for (long r = 2; r < 5; ++r) {
        INFO("r=", r);
        CHECK(pfqn_scbgap<double>(8, 5, r, true) < pfqn_scbgap<double>(8, 5, r, false));
    }
    CHECK(pfqn_scbgap<double>(8, 5, 5, true) == doctest::Approx(pfqn_scbgap<double>(8, 5, 5, false)));
    CHECK_THROWS_AS(pfqn_scbgap<double>(8, 5, 6, true), line::InputError);
    // Merging a single class cannot change anything.
    CHECK(pfqn_scbgap<double>(6, 4, 1, false) == 0.0);
}

TEST_CASE("pfqn_usumbound and pfqn_minclasses reproduce the Section 4.7 example") {
    // K = 2 devices, N = 3 customers: one class admits at most 2N/(N+1) = 1.5.
    CHECK(pfqn_usumbound<double>(1, 2, 3) == doctest::Approx(1.5));
    CHECK(pfqn_usumbound<double>(2, 2, 3) == doctest::Approx(2.0));
    CHECK(pfqn_usumbound<Rational>(1, 2, 3) == Rational(3, 2));
    // A measured 1.6 refutes the single-class assumption; 1.4 does not.
    CHECK(pfqn_minclasses<double>(1.6, 2, 3) == 2);
    CHECK(pfqn_minclasses<double>(1.4, 2, 3) == 1);
    // Above min(N,K) no class structure can explain the measurement.
    CHECK(pfqn_minclasses<double>(2.5, 2, 3) == -1);
    // Nondecreasing in R, which is what makes the inversion well posed.
    for (long R = 1; R < 5; ++R)
        CHECK(pfqn_usumbound<double>(R, 4, 5) <= pfqn_usumbound<double>(R + 1, 4, 5));
}

TEST_CASE("pfqn_scb brackets the multiclass throughput of the Section 2 example") {
    // Single-class demands (0.114, 0.040, 0.062) at N = 4; the paper's
    // multiclass counterpart runs at X_R = 8.761 and the aggregate at 8.152.
    const std::vector<double> L{0.114, 0.040, 0.062};
    const ScbBounds<double> b = pfqn_scb(L, 4);
    CHECK(b.Xlo == doctest::Approx(8.151894490678735).epsilon(1e-12));
    CHECK(b.Xlo <= 8.7615);
    CHECK(b.Xhi >= 8.7615);
    // Corollary 1 makes the utilization ratio uniform across devices.
    for (std::size_t k = 0; k < 3; ++k) {
        INFO("device ", k);
        CHECK(b.Ulo[k] == doctest::Approx(b.Xlo * L[k]).epsilon(1e-12));
        CHECK(b.Uhi[k] / b.Ulo[k] == doctest::Approx(b.Xhi / b.Xlo).epsilon(1e-12));
    }
    // The single-server cap binds here, so the busiest device sits exactly at 1.
    CHECK(b.Uhi[0] == doctest::Approx(1.0).epsilon(1e-12));

    // The lower side IS the exact single-class solution, not an approximation.
    Matrix<double> Lm(3, 1);
    for (std::size_t k = 0; k < 3; ++k) Lm(k, 0) = L[k];
    const auto mva = pfqn_mva(Lm, std::vector<int>{4});
    CHECK(b.Xlo == doctest::Approx(mva.XN[0]).epsilon(1e-12));

    // N = 1: no aggregation error is possible, so the bracket collapses.
    const ScbBounds<double> b1 = pfqn_scb(L, 1);
    CHECK(b1.Xhi == doctest::Approx(b1.Xlo).epsilon(1e-12));

    // Exact arithmetic: the recursion and the scaling stay in the field.
    const std::vector<Rational> Lr{Rational(1, 2), Rational(1, 3), Rational(1, 5)};
    const ScbBounds<Rational> br = pfqn_scb(Lr, 5);
    CHECK(br.Xlo <= br.Xhi);
    CHECK(br.Uhi[0] <= Rational(1));
}
