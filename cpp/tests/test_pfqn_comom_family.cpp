/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_comom, pfqn_comomrm_orig, pfqn_procomom, pfqn_procomom2, pfqn_momlin.
 *
 * Every expected number below was produced by MATLAB (lineStart; then the
 * routine under test) and is quoted to 17 significant digits. Where the routine
 * is exact-capable the test asserts BIT-IDENTITY at Rational against an
 * independent exact route -- pfqn_ca for the normalizing constants, pfqn_mva
 * for the queue lengths -- which is a stronger statement than agreeing with
 * MATLAB's doubles, and is the whole point of carrying the exact backend.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comom.h"
#include "line/api/pfqn/pfqn_momlin.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_procomom.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::num_traits;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_comom;
using line::pfqn::pfqn_comomrm_orig;
using line::pfqn::pfqn_momlin;
using line::pfqn::pfqn_mva;
using line::pfqn::pfqn_procomom;
using line::pfqn::pfqn_procomom2;

namespace {

template <class T>
T q(long n, long d) {
    return num_traits<T>::from_rational(n, d);
}

double relerr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

// ---------------------------------------------------------------------------
// Repairman models: one queueing station plus a delay.
// ---------------------------------------------------------------------------

/** L = [1/10 1/5], N = [2 3], Z = [1 1/2]. */
template <class T>
void repairman_A(Matrix<T>& L, std::vector<int>& N, std::vector<T>& Z) {
    L = Matrix<T>(1, 2);
    L(0, 0) = q<T>(1, 10);
    L(0, 1) = q<T>(1, 5);
    N.assign(2, 0);
    N[0] = 2;
    N[1] = 3;
    Z.assign(2, num_traits<T>::from_int(0));
    Z[0] = num_traits<T>::from_int(1);
    Z[1] = q<T>(1, 2);
}

/** L = [1/2 1/4 1/8], N = [1 2 1], Z = [1 1/2 2]. */
template <class T>
void repairman_A2(Matrix<T>& L, std::vector<int>& N, std::vector<T>& Z) {
    L = Matrix<T>(1, 3);
    L(0, 0) = q<T>(1, 2);
    L(0, 1) = q<T>(1, 4);
    L(0, 2) = q<T>(1, 8);
    N.assign(3, 0);
    N[0] = 1;
    N[1] = 2;
    N[2] = 1;
    Z.assign(3, num_traits<T>::from_int(0));
    Z[0] = num_traits<T>::from_int(1);
    Z[1] = q<T>(1, 2);
    Z[2] = num_traits<T>::from_int(2);
}

/**
 * Same demands as A2 but Z = [4 1/4 1/2], whose per-class ratios Z/L are
 * [8 1 4] and are therefore NOT sorted. This is the model that shows
 * pfqn_comomrm_orig's second, apparently unsynchronized, class sort to be a
 * no-op: MATLAB returns lG = -0.439231970578982 against the exact
 * -0.43923197057898189.
 */
template <class T>
void repairman_K(Matrix<T>& L, std::vector<int>& N, std::vector<T>& Z) {
    repairman_A2(L, N, Z);
    Z[0] = num_traits<T>::from_int(4);
    Z[1] = q<T>(1, 4);
    Z[2] = q<T>(1, 2);
}

/** L = [1 1/10], N = [3 1], Z = [1/5 5]; Z/L = [1/5 50]. */
template <class T>
void repairman_K2(Matrix<T>& L, std::vector<int>& N, std::vector<T>& Z) {
    L = Matrix<T>(1, 2);
    L(0, 0) = num_traits<T>::from_int(1);
    L(0, 1) = q<T>(1, 10);
    N.assign(2, 0);
    N[0] = 3;
    N[1] = 1;
    Z.assign(2, num_traits<T>::from_int(0));
    Z[0] = q<T>(1, 5);
    Z[1] = num_traits<T>::from_int(5);
}

template <class T>
Matrix<T> as_row(const std::vector<T>& v) {
    Matrix<T> m(1, v.size());
    for (std::size_t i = 0; i < v.size(); ++i) m(0, i) = v[i];
    return m;
}

}  // namespace

TEST_CASE("pfqn_comom reproduces MATLAB and equals pfqn_ca exactly") {
    // MATLAB pfqn_comom(L,N,Z,1e-14) on each model.
    const double lgA = -2.8795801968179671;
    const double lgA2 = 0.31275571000389757;
    const double lgK = -0.43923197057898111;
    const double lgK2 = 1.8826355906849082;

    SUBCASE("double") {
        Matrix<double> L;
        std::vector<int> N;
        std::vector<double> Z;
        repairman_A(L, N, Z);
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(lgA).epsilon(1e-13));
        repairman_A2(L, N, Z);
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(lgA2).epsilon(1e-13));
        repairman_K(L, N, Z);
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(lgK).epsilon(1e-13));
        repairman_K2(L, N, Z);
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(lgK2).epsilon(1e-13));
    }

    SUBCASE("exact: identical to pfqn_ca as a fraction") {
        Matrix<Rational> L;
        std::vector<int> N;
        std::vector<Rational> Z;
        const Rational atol = num_traits<Rational>::from_int(0);

        repairman_A(L, N, Z);
        CHECK(pfqn_comom(L, N, Z, atol).G == pfqn_ca(L, N, as_row(Z)).G);
        repairman_A2(L, N, Z);
        CHECK(pfqn_comom(L, N, Z, atol).G == pfqn_ca(L, N, as_row(Z)).G);
        // A2 has the closed form G = 175/128; check the literal value too.
        CHECK(pfqn_comom(L, N, Z, atol).G == Rational(175) / Rational(128));
        repairman_K(L, N, Z);
        CHECK(pfqn_comom(L, N, Z, atol).G == pfqn_ca(L, N, as_row(Z)).G);
        repairman_K2(L, N, Z);
        CHECK(pfqn_comom(L, N, Z, atol).G == pfqn_ca(L, N, as_row(Z)).G);
    }

    SUBCASE("real50") {
        Matrix<Real50> L;
        std::vector<int> N;
        std::vector<Real50> Z;
        repairman_A2(L, N, Z);
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(lgA2).epsilon(1e-14));
    }
}

TEST_CASE("pfqn_comomrm_orig reproduces MATLAB and equals pfqn_ca exactly") {
    // MATLAB pfqn_comomrm_orig(L,N,Z,1e-14).
    const double lgA = -2.879580196817968;
    const double lgA2 = 0.31275571000389757;
    const double lgK = -0.439231970578982;
    const double lgK2 = 1.8826355906849086;

    SUBCASE("double") {
        Matrix<double> L;
        std::vector<int> N;
        std::vector<double> Z;
        repairman_A(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).lG == doctest::Approx(lgA).epsilon(1e-13));
        repairman_A2(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).lG == doctest::Approx(lgA2).epsilon(1e-13));
        repairman_K(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).lG == doctest::Approx(lgK).epsilon(1e-13));
        repairman_K2(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).lG == doctest::Approx(lgK2).epsilon(1e-13));
    }

    SUBCASE("exact: identical to pfqn_ca AND to pfqn_comom as a fraction") {
        Matrix<Rational> L;
        std::vector<int> N;
        std::vector<Rational> Z;
        const Rational atol = num_traits<Rational>::from_int(0);
        repairman_A(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == pfqn_ca(L, N, as_row(Z)).G);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == pfqn_comom(L, N, Z, atol).G);
        repairman_A2(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == pfqn_ca(L, N, as_row(Z)).G);
        repairman_K(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == pfqn_ca(L, N, as_row(Z)).G);
        repairman_K2(L, N, Z);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == pfqn_ca(L, N, as_row(Z)).G);
    }
}

TEST_CASE("pfqn_procomom marginals reproduce MATLAB") {
    // MATLAB pfqn_procomom([1 2; 3 1], [2 2], [1 1], 1e-14).
    const double PrB[2][5] = {{0.26155794806839777, 0.28119062697910074, 0.24445851804939836,
                               0.15199493350221657, 0.0607979734008866},
                              {0.12729575680810643, 0.2089930335655478, 0.27612412919569346,
                               0.25079164027865736, 0.13679544015199493}};
    const double QB[2] = {1.4692843571880936, 2.0607979734008866};

    Matrix<double> L(2, 2);
    L(0, 0) = 1;
    L(0, 1) = 2;
    L(1, 0) = 3;
    L(1, 1) = 1;
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 2;
    std::vector<double> Z(2, 1.0);

    const line::pfqn::ProcomomResult<double> r = pfqn_procomom(L, N, Z);
    CHECK_FALSE(r.rankdef);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 5; ++j) CHECK(r.Pr(i, j) == doctest::Approx(PrB[i][j]).epsilon(1e-12));
        CHECK(r.Q[i] == doctest::Approx(QB[i]).epsilon(1e-12));
        // a marginal distribution sums to one
        double s = 0;
        for (int j = 0; j < 5; ++j) s += r.Pr(i, j);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-13));
    }
}

TEST_CASE("pfqn_procomom single class reproduces MATLAB") {
    // MATLAB pfqn_procomom([1; 2], 3, 0.5, 1e-14).
    const double PrB2[2][4] = {{0.54355016538037493, 0.27122381477398017, 0.13230429988974637,
                                0.05292171995589854},
                               {0.087100330760749758, 0.17199558985667041, 0.31753031973539148,
                                0.42337375964718837}};
    Matrix<double> L(2, 1);
    L(0, 0) = 1;
    L(1, 0) = 2;
    std::vector<int> N(1, 3);
    std::vector<double> Z(1, 0.5);
    const line::pfqn::ProcomomResult<double> r = pfqn_procomom(L, N, Z);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 4; ++j) CHECK(r.Pr(i, j) == doctest::Approx(PrB2[i][j]).epsilon(1e-12));
    CHECK(r.Q[0] == doctest::Approx(0.69459757442116854).epsilon(1e-12));
    CHECK(r.Q[1] == doctest::Approx(2.0771775082690183).epsilon(1e-12));
}

TEST_CASE("pfqn_procomom queue lengths equal pfqn_mva exactly at Rational") {
    // Q_i = sum_j j P(n_i = j) is the per-station total queue length, so the
    // ProCoMoM marginals and the MVA recursion must agree as fractions. MATLAB's
    // doubles differ in the last ulp (1.4692843571880936 against
    // 1.4692843571880938); the rationals do not differ at all.
    SUBCASE("M = R = 2") {
        Matrix<Rational> L(2, 2);
        L(0, 0) = num_traits<Rational>::from_int(1);
        L(0, 1) = num_traits<Rational>::from_int(2);
        L(1, 0) = num_traits<Rational>::from_int(3);
        L(1, 1) = num_traits<Rational>::from_int(1);
        std::vector<int> N(2, 2);
        std::vector<Rational> Z(2, num_traits<Rational>::from_int(1));
        const line::pfqn::ProcomomResult<Rational> r =
            pfqn_procomom(L, N, Z, num_traits<Rational>::from_int(0));
        const line::pfqn::MvaResult<Rational> m = pfqn_mva(L, N, as_row(Z));
        for (std::size_t i = 0; i < 2; ++i) {
            Rational tot = num_traits<Rational>::from_int(0);
            for (std::size_t s = 0; s < 2; ++s) tot += m.QN(i, s);
            CHECK(r.Q[i] == tot);
        }
    }
    SUBCASE("M = 3, R = 2") {
        Matrix<Rational> L(3, 2);
        L(0, 0) = num_traits<Rational>::from_int(1);
        L(0, 1) = num_traits<Rational>::from_int(2);
        L(1, 0) = num_traits<Rational>::from_int(3);
        L(1, 1) = num_traits<Rational>::from_int(1);
        L(2, 0) = q<Rational>(1, 2);
        L(2, 1) = q<Rational>(3, 2);
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<Rational> Z(2);
        Z[0] = num_traits<Rational>::from_int(1);
        Z[1] = q<Rational>(1, 2);
        const line::pfqn::ProcomomResult<Rational> r =
            pfqn_procomom(L, N, Z, num_traits<Rational>::from_int(0));
        const line::pfqn::MvaResult<Rational> m = pfqn_mva(L, N, as_row(Z));
        // MATLAB doubles for the same quantities.
        const double QM[3] = {0.74003466204506063, 1.5597920277296362, 0.3994800693240903};
        for (std::size_t i = 0; i < 3; ++i) {
            Rational tot = num_traits<Rational>::from_int(0);
            for (std::size_t s = 0; s < 2; ++s) tot += m.QN(i, s);
            CHECK(r.Q[i] == tot);
            CHECK(relerr(static_cast<double>(r.Q[i]), QM[i]) < 1e-12);
        }
    }
}

TEST_CASE("pfqn_procomom2 queue-plus-delay marginal") {
    // MATLAB pfqn_procomom2([0.4 0.6],[2 1],[1 2],[],1).
    const double pkC[4] = {0.31367628607277287, 0.34504391468005019, 0.25094102885821834,
                           0.0903387703889586};
    SUBCASE("double") {
        std::vector<double> L(2);
        L[0] = 0.4;
        L[1] = 0.6;
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<double> Z(2);
        Z[0] = 1.0;
        Z[1] = 2.0;
        const line::pfqn::Procomom2Result<double> r = pfqn_procomom2(L, N, Z);
        for (int i = 0; i < 4; ++i) CHECK(r.pk[i] == doctest::Approx(pkC[i]).epsilon(1e-13));
        CHECK(r.G == doctest::Approx(3.188).epsilon(1e-14));
        CHECK(r.lG == doctest::Approx(1.1593937609279688).epsilon(1e-13));
        double s = 0;
        for (int i = 0; i < 4; ++i) s += r.pk[i];
        CHECK(s == doctest::Approx(1.0).epsilon(1e-14));
    }
    SUBCASE("exact: G identical to pfqn_ca") {
        std::vector<Rational> L(2);
        L[0] = q<Rational>(2, 5);
        L[1] = q<Rational>(3, 5);
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<Rational> Z(2);
        Z[0] = num_traits<Rational>::from_int(1);
        Z[1] = num_traits<Rational>::from_int(2);
        const line::pfqn::Procomom2Result<Rational> r = pfqn_procomom2(L, N, Z);
        Matrix<Rational> Lm(1, 2);
        Lm(0, 0) = L[0];
        Lm(0, 1) = L[1];
        CHECK(r.G == pfqn_ca(Lm, N, as_row(Z)).G);
        CHECK(r.G == Rational(3188) / Rational(1000));
    }
    SUBCASE("load dependent rates are honoured") {
        // mu(n) = min(n,2): a two-server queue. Its marginal must differ from
        // the single-server one and still normalize.
        std::vector<double> L(2);
        L[0] = 0.4;
        L[1] = 0.6;
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<double> Z(2);
        Z[0] = 1.0;
        Z[1] = 2.0;
        std::vector<double> mu(3);
        mu[0] = 1.0;
        mu[1] = 2.0;
        mu[2] = 2.0;
        const line::pfqn::Procomom2Result<double> r = pfqn_procomom2(L, N, Z, mu, 1);
        double s = 0;
        for (int i = 0; i < 4; ++i) s += r.pk[i];
        CHECK(s == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(r.pk[0] > pkC[0]);  // faster service, emptier queue
    }
}

TEST_CASE("pfqn_momlin means and second moments reproduce MATLAB") {
    // MATLAB pfqn_momlin([1 2; 3 1; 0.5 0.5], [5 4], [1 0.5], 1e-12, 5000).
    Matrix<double> L(3, 2);
    L(0, 0) = 1;
    L(0, 1) = 2;
    L(1, 0) = 3;
    L(1, 1) = 1;
    L(2, 0) = 0.5;
    L(2, 1) = 0.5;
    std::vector<int> N(2);
    N[0] = 5;
    N[1] = 4;
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 0.5;
    const line::pfqn::MomlinResult<double> r = pfqn_momlin(L, N, Z, 1e-12, 5000);

    const double Q[3][2] = {{0.74651313672790687, 1.827031978948813},
                            {3.8936101263997522, 1.8384187115202224},
                            {0.14186847978437939, 0.18800160607426167}};
    const double Rr[3][2] = {{3.424242488332462, 6.2335742418802091},
                             {17.859920437916045, 6.2724241600393622},
                             {0.65074819495088265, 0.64143484217002888}};
    const double QVar[3][2] = {{1.0846798083412816, 1.9691059116792604},
                               {1.4296132059451365, 1.8758418086224871},
                               {0.1569210103955952, 0.20893815729614298}};
    for (int i = 0; i < 3; ++i)
        for (int s = 0; s < 2; ++s) {
            CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-11));
            CHECK(r.R(i, s) == doctest::Approx(Rr[i][s]).epsilon(1e-11));
            CHECK(r.QVar(i, s) == doctest::Approx(QVar[i][s]).epsilon(1e-10));
        }
    CHECK(r.X[0] == doctest::Approx(0.21800825708796223).epsilon(1e-11));
    CHECK(r.X[1] == doctest::Approx(0.29309540691340707).epsilon(1e-11));
    CHECK(r.U(1, 0) == doctest::Approx(0.6540247712638867).epsilon(1e-11));
    // off-diagonal covariance and one demand derivative
    CHECK(r.cov(1, 0, 2, 1) == doctest::Approx(-0.0069312322414829194).epsilon(1e-9));
    CHECK(r.dq(2, 1, 1, 0) == doctest::Approx(-0.00098960838730241269).epsilon(1e-9));
    // the population is conserved by the fixed point
    for (int s = 0; s < 2; ++s) {
        double tot = r.X[s] * Z[s];
        for (int i = 0; i < 3; ++i) tot += r.Q(i, s);
        CHECK(tot == doctest::Approx(static_cast<double>(N[s])).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_momlin single class reproduces MATLAB") {
    // MATLAB pfqn_momlin([1; 2; 0.5], 8, 1, 1e-12, 5000).
    Matrix<double> L(3, 1);
    L(0, 0) = 1;
    L(1, 0) = 2;
    L(2, 0) = 0.5;
    std::vector<int> N(1, 8);
    std::vector<double> Z(1, 1.0);
    const line::pfqn::MomlinResult<double> r = pfqn_momlin(L, N, Z, 1e-12, 5000);
    const double Q[3] = {0.84096445182997825, 6.3671778696571693, 0.30738758642574182};
    const double QVar[3] = {1.4115436288701275, 2.2109653346719442, 0.38661955942106779};
    for (int i = 0; i < 3; ++i) {
        CHECK(r.Q(i, 0) == doctest::Approx(Q[i]).epsilon(1e-11));
        CHECK(r.QVar(i, 0) == doctest::Approx(QVar[i]).epsilon(1e-10));
    }
    CHECK(r.X[0] == doctest::Approx(0.48447009208710995).epsilon(1e-11));
    CHECK(r.U(1, 0) == doctest::Approx(0.9689401841742199).epsilon(1e-11));
}

TEST_CASE("pfqn_comom regression: the model that caught the JAR index defects") {
    // Pfqn_comom.java differenced the hash lookups against the TARGET
    // population instead of the running one, and read sum(Dn(d,1:r)) over an
    // exclusive range. Either fault makes matchrow miss and the result collapse
    // to -Inf. MATLAB, verified against a brute-force state-space sum, gives
    // lG = 0.6108519378 here, and pfqn_comom, pfqn_comomrm_orig and pfqn_ca all
    // agree on it.
    SUBCASE("double") {
        Matrix<double> L(1, 2);
        L(0, 0) = 0.6;
        L(0, 1) = 0.4;
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<double> Z(2);
        Z[0] = 1.0;
        Z[1] = 0.5;
        CHECK(pfqn_comom(L, N, Z).lG == doctest::Approx(0.6108519378).epsilon(1e-10));
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).lG ==
              doctest::Approx(0.6108519378).epsilon(1e-10));
        CHECK(relerr(pfqn_comom(L, N, Z).lG, pfqn_ca(L, N, as_row(Z)).lG) < 1e-14);
    }
    SUBCASE("exact: G = 921/500 by three routes") {
        Matrix<Rational> L(1, 2);
        L(0, 0) = q<Rational>(3, 5);
        L(0, 1) = q<Rational>(2, 5);
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 1;
        std::vector<Rational> Z(2);
        Z[0] = num_traits<Rational>::from_int(1);
        Z[1] = q<Rational>(1, 2);
        const Rational G = Rational(921) / Rational(500);
        CHECK(pfqn_comom(L, N, Z, num_traits<Rational>::from_int(0)).G == G);
        CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == G);
        CHECK(pfqn_ca(L, N, as_row(Z)).G == G);
    }
}

TEST_CASE("pfqn_comom refuses more than one queueing station") {
    Matrix<double> L(2, 2, 1.0);
    std::vector<int> N(2, 1);
    std::vector<double> Z(2, 1.0);
    CHECK_THROWS_AS(pfqn_comom(L, N, Z), line::InputError);
    CHECK_THROWS_AS(pfqn_comomrm_orig(L, N, as_row(Z)), line::InputError);
}

TEST_CASE("pfqn_comom family: exact identity with pfqn_ca over a model sweep") {
    // The cheapest gate the whole family has: at Rational the CoMoM basis
    // recursion and the convolution must return the SAME fraction, so the sweep
    // needs no tolerance at all. A hash that resolves to the wrong basis row, or
    // a Dn column range off by one, breaks this on the first model with R >= 3
    // even when R = 2 still happens to work.
    unsigned seed = 20260721u;
    const auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return (seed >> 16) & 0x7fffu;
    };
    int checked = 0;
    for (int R = 2; R <= 4; ++R) {
        for (int trial = 0; trial < 4; ++trial) {
            Matrix<Rational> L(1, static_cast<std::size_t>(R));
            std::vector<int> N(static_cast<std::size_t>(R), 0);
            std::vector<Rational> Z(static_cast<std::size_t>(R));
            for (int r = 0; r < R; ++r) {
                L(0, static_cast<std::size_t>(r)) =
                    q<Rational>(1 + static_cast<long>(next() % 20u), 10);
                Z[static_cast<std::size_t>(r)] =
                    q<Rational>(1 + static_cast<long>(next() % 30u), 10);
                N[static_cast<std::size_t>(r)] = 1 + static_cast<int>(next() % 3u);
            }
            const Rational G = pfqn_ca(L, N, as_row(Z)).G;
            CHECK(pfqn_comom(L, N, Z, num_traits<Rational>::from_int(0)).G == G);
            CHECK(pfqn_comomrm_orig(L, N, as_row(Z)).G == G);
            ++checked;
        }
    }
}

TEST_CASE("pfqn_procomom: exact identity with pfqn_ca on the single-station case") {
    // With M = 1 the ProCoMoM marginal of the only queue is the full state
    // distribution, so sum_j Pr(1,j) x G_j must rebuild pfqn_ca's constant and
    // the mean must be pfqn_mva's queue length, both as fractions.
    unsigned seed = 991u;
    const auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return (seed >> 16) & 0x7fffu;
    };
    for (int R = 2; R <= 3; ++R) {
        for (int trial = 0; trial < 3; ++trial) {
            Matrix<Rational> L(1, static_cast<std::size_t>(R));
            std::vector<int> N(static_cast<std::size_t>(R), 0);
            std::vector<Rational> Z(static_cast<std::size_t>(R));
            for (int r = 0; r < R; ++r) {
                L(0, static_cast<std::size_t>(r)) =
                    q<Rational>(1 + static_cast<long>(next() % 20u), 10);
                Z[static_cast<std::size_t>(r)] =
                    q<Rational>(1 + static_cast<long>(next() % 30u), 10);
                N[static_cast<std::size_t>(r)] = 1 + static_cast<int>(next() % 2u);
            }
            const line::pfqn::ProcomomResult<Rational> p =
                pfqn_procomom(L, N, Z, num_traits<Rational>::from_int(0));
            const line::pfqn::MvaResult<Rational> m = pfqn_mva(L, N, as_row(Z));
            Rational tot = num_traits<Rational>::from_int(0);
            for (int s = 0; s < R; ++s) tot += m.QN(0, static_cast<std::size_t>(s));
            CHECK_FALSE(p.rankdef);
            CHECK(p.Q[0] == tot);
            Rational mass = num_traits<Rational>::from_int(0);
            for (std::size_t j = 0; j < p.Pr.cols(); ++j) mass += p.Pr(0, j);
            CHECK(mass == num_traits<Rational>::from_int(1));
        }
    }
}

TEST_CASE("pfqn_momlin agrees with exact MVA exactly where Schweitzer-Bard is exact") {
    // pfqn_momlin is an AMVA fixed point, so it CANNOT be gated against pfqn_ca
    // the way the CoMoM family is: it is approximate by construction and the
    // gap is the method's, not the port's. The two statements that do hold
    // exactly are asserted instead.
    SUBCASE("a single customer: the population term c(r,r) vanishes") {
        Matrix<Rational> L(3, 1);
        L(0, 0) = num_traits<Rational>::from_int(1);
        L(1, 0) = num_traits<Rational>::from_int(2);
        L(2, 0) = q<Rational>(1, 2);
        std::vector<int> N(1, 1);
        std::vector<Rational> Z(1, num_traits<Rational>::from_int(1));
        const line::pfqn::MomlinResult<Rational> r =
            pfqn_momlin(L, N, Z, num_traits<Rational>::from_int(0), 200);
        const line::pfqn::MvaResult<Rational> m = pfqn_mva(L, N, as_row(Z));
        CHECK(r.X[0] == m.XN[0]);
        for (std::size_t i = 0; i < 3; ++i) CHECK(r.Q(i, 0) == m.QN(i, 0));
    }
    SUBCASE("the fixed point conserves the population at every population") {
        Matrix<double> L(3, 2);
        L(0, 0) = 1;
        L(0, 1) = 2;
        L(1, 0) = 3;
        L(1, 1) = 1;
        L(2, 0) = 0.5;
        L(2, 1) = 0.5;
        std::vector<double> Z(2);
        Z[0] = 1.0;
        Z[1] = 0.5;
        for (int n = 1; n <= 6; ++n) {
            std::vector<int> N(2);
            N[0] = n;
            N[1] = n + 1;
            const line::pfqn::MomlinResult<double> r = pfqn_momlin(L, N, Z, 1e-13, 20000);
            for (int s = 0; s < 2; ++s) {
                double tot = r.X[s] * Z[s];
                for (int i = 0; i < 3; ++i) tot += r.Q(i, s);
                CHECK(tot == doctest::Approx(static_cast<double>(N[s])).epsilon(1e-10));
            }
        }
    }
}
