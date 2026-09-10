/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Exactness tests for the normalizing-constant family: convolution with class
 * dependence, the CoMoM repairman variants, the load-dependent single-class
 * recursion and the effective-capacity terms of the mixed solver.
 *
 * The governing assertion is that every EXACT method returns the IDENTICAL
 * rational G as pfqn_ca on the same model. That is asserted with operator==
 * on line::Rational, not with a tolerance: two exact algorithms that disagree
 * in the last bit disagree, and a tolerance would hide it.
 */
#include <algorithm>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_comomrm_ld.h"
#include "line/api/pfqn/pfqn_comomrm_ms.h"
#include "line/api/pfqn/pfqn_conv.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_gldsingle.h"
#include "line/api/pfqn/pfqn_lldsingle.h"
#include "line/api/pfqn/pfqn_ldmx_ec.h"
#include "line/api/pfqn/pfqn_mu_ms.h"
#include "line/api/pfqn/pfqn_nc_sanitize.h"

using line::Matrix;
using line::Rational;

namespace {

Matrix<Rational> rmat(std::initializer_list<std::initializer_list<int>> rows, int den) {
    std::vector<std::vector<int>> v;
    for (const auto& r : rows) v.emplace_back(r);
    Matrix<Rational> m(v.size(), v.empty() ? 0 : v[0].size());
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v[i].size(); ++j) m(i, j) = Rational(v[i][j], den);
    return m;
}

}  // namespace

TEST_CASE("pfqn_conv with no class dependence reproduces pfqn_ca exactly") {
    // Demands and think times as exact tenths, so the two algorithms are
    // compared on identical rational inputs.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}, {2, 2}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);

    const std::vector<std::vector<int>> pops = {{1, 0}, {0, 3}, {2, 2}, {4, 3}, {3, 5}};
    for (const std::vector<int>& N : pops) {
        const Rational gca = line::pfqn::pfqn_ca(L, N, Z).G;
        const Rational gcv = line::pfqn::pfqn_conv(L, N, Z).G;
        CHECK(gcv == gca);
    }
}

TEST_CASE("pfqn_conv with a unit class-dependent scaling reproduces pfqn_ca exactly") {
    // beta = 1 is the documented "no correction" value, so the class-dependent
    // convolution path must return the load-independent constant. This
    // exercises the O(P^2) direct convolution branch, not the Buzen shortcut.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {3, 2};

    std::vector<line::pfqn::CdScaling<Rational>> cd(2);
    cd[0] = [](const std::vector<Rational>&) {
        return std::vector<Rational>{Rational(1)};
    };
    // Station 1 stays load independent, so the two branches are mixed.
    const Rational gca = line::pfqn::pfqn_ca(L, N, Z).G;
    CHECK(line::pfqn::pfqn_conv(L, N, Z, cd).G == gca);
}

TEST_CASE("pfqn_conv class dependence halves the demand when beta = 2") {
    // A constant beta = 2 divides the effective demand of that station by two,
    // so the constant equals the one of the model with that row halved.
    const Matrix<Rational> L = rmat({{6, 4}, {3, 7}}, 10);
    const Matrix<Rational> Lhalf = rmat({{3, 2}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 3};

    std::vector<line::pfqn::CdScaling<Rational>> cd(2);
    cd[0] = [](const std::vector<Rational>&) {
        return std::vector<Rational>{Rational(2)};
    };
    CHECK(line::pfqn::pfqn_conv(L, N, Z, cd).G == line::pfqn::pfqn_ca(Lhalf, N, Z).G);
}

TEST_CASE("pfqn_comomrm matches pfqn_ca exactly on repairman models") {
    // One queueing station plus a delay, which is the whole domain of CoMoM-RM.
    const Matrix<Rational> L = rmat({{7, 4, 5}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9, 2}}, 10);

    const std::vector<std::vector<int>> pops = {
        {2, 0, 0}, {1, 1, 0}, {3, 2, 0}, {2, 2, 2}, {4, 1, 3}, {5, 4, 2}};
    for (const std::vector<int>& N : pops) {
        const Rational gca = line::pfqn::pfqn_ca(L, N, Z).G;
        CHECK(line::pfqn::pfqn_comomrm(L, N, Z).G == gca);
    }
}

TEST_CASE("pfqn_comomrm handles zero-think-time classes") {
    // The zero-think-time head is seeded from a ratio of factorials rather than
    // from the recursion, so it needs its own coverage.
    const Matrix<Rational> L = rmat({{7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{0, 5}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {3, 1}, {1, 4}}) {
        CHECK(line::pfqn::pfqn_comomrm(L, N, Z).G == line::pfqn::pfqn_ca(L, N, Z).G);
    }
}

TEST_CASE("pfqn_comomrm with multiplicity m equals m replicated stations") {
    // The multiplicity argument stands for m identical single-server stations,
    // so it must agree with the convolution over an explicitly replicated L.
    const Matrix<Rational> L1 = rmat({{7, 4}}, 10);
    const Matrix<Rational> L3 = rmat({{7, 4}, {7, 4}, {7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 1}, {3, 2}}) {
        CHECK(line::pfqn::pfqn_comomrm(L1, N, Z, 3).G == line::pfqn::pfqn_ca(L3, N, Z).G);
    }
}

TEST_CASE("pfqn_comomrm_ms matches pfqn_ca on the single-server case") {
    // S = 1, m = 1 degenerates the multiserver lattice to all ones.
    const Matrix<Rational> L = rmat({{7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {4, 1}, {3, 3}}) {
        const line::pfqn::ComomRmResult<Rational> r = line::pfqn::pfqn_comomrm_ms(L, N, Z, 1, 1);
        CHECK(r.G == line::pfqn::pfqn_ca(L, N, Z).G);
        Rational psum(0);
        for (const Rational& p : r.prob) psum += p;
        CHECK(psum == Rational(1));
    }
}

TEST_CASE("pfqn_comomrm_ms matches pfqn_gld on a multiserver station") {
    // A c-server station is the load-dependent station with rate min(k, c), so
    // the CoMoM-MS constant must equal the generalized convolution's.
    const Matrix<Rational> L = rmat({{7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9}}, 10);
    const int S = 3;
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {3, 3}, {4, 2}}) {
        int Nt = 0;
        for (int v : N) Nt += v;
        // pfqn_gld takes no think time: the delay is a row with rate lattice k.
        Matrix<Rational> Lg(2, 2), mug(2, static_cast<std::size_t>(Nt));
        for (std::size_t r = 0; r < 2; ++r) {
            Lg(0, r) = L(0, r);
            Lg(1, r) = Z(0, r);
        }
        for (int k = 1; k <= Nt; ++k) {
            mug(0, static_cast<std::size_t>(k - 1)) = Rational(k < S ? k : S);
            mug(1, static_cast<std::size_t>(k - 1)) = Rational(k);
        }
        CHECK(line::pfqn::pfqn_comomrm_ms(L, N, Z, 1, S).G == line::pfqn::pfqn_gld(Lg, N, mug).G);
    }
}

TEST_CASE("pfqn_comomrm_ld matches pfqn_gld on an arbitrary rate lattice") {
    const Matrix<Rational> L = rmat({{7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {3, 2}, {2, 4}}) {
        int Nt = 0;
        for (int v : N) Nt += v;
        // A non-monotone, non-multiserver lattice, to avoid accidentally
        // testing only the min(k,c) shape.
        Matrix<Rational> mu(1, static_cast<std::size_t>(Nt));
        for (int k = 1; k <= Nt; ++k)
            mu(0, static_cast<std::size_t>(k - 1)) = Rational(3 + (k % 3), 2);

        Matrix<Rational> Lg(2, 2), mug(2, static_cast<std::size_t>(Nt));
        for (std::size_t r = 0; r < 2; ++r) {
            Lg(0, r) = L(0, r);
            Lg(1, r) = Z(0, r);
        }
        for (int k = 1; k <= Nt; ++k) {
            mug(0, static_cast<std::size_t>(k - 1)) = mu(0, static_cast<std::size_t>(k - 1));
            mug(1, static_cast<std::size_t>(k - 1)) = Rational(k);
        }
        CHECK(line::pfqn::pfqn_comomrm_ld(L, N, Z, mu).G == line::pfqn::pfqn_gld(Lg, N, mug).G);
    }
}

TEST_CASE("pfqn_gldsingle matches pfqn_gld and pfqn_ca in one class") {
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(7, 10);
    L(1, 0) = Rational(4, 10);
    for (int N : {1, 2, 5, 8}) {
        Matrix<Rational> mu(2, static_cast<std::size_t>(N));
        for (int k = 1; k <= N; ++k) {
            mu(0, static_cast<std::size_t>(k - 1)) = Rational(1);
            mu(1, static_cast<std::size_t>(k - 1)) = Rational(1);
        }
        const std::vector<int> Nv = {N};
        CHECK(line::pfqn::pfqn_gldsingle(L, N, mu).G == line::pfqn::pfqn_ca(L, Nv).G);
        CHECK(line::pfqn::pfqn_gldsingle(L, N, mu).G == line::pfqn::pfqn_gld(L, Nv, mu).G);
    }
}

TEST_CASE("pfqn_gldsingle handles a delay row expressed as the rate lattice k") {
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(7, 10);
    L(1, 0) = Rational(3, 10);  // delay demand
    for (int N : {1, 3, 6}) {
        Matrix<Rational> mu(2, static_cast<std::size_t>(N));
        for (int k = 1; k <= N; ++k) {
            mu(0, static_cast<std::size_t>(k - 1)) = Rational(1);
            mu(1, static_cast<std::size_t>(k - 1)) = Rational(k);
        }
        Matrix<Rational> Lq(1, 1), Zq(1, 1);
        Lq(0, 0) = Rational(7, 10);
        Zq(0, 0) = Rational(3, 10);
        const std::vector<int> Nv = {N};
        CHECK(line::pfqn::pfqn_gldsingle(L, N, mu).G == line::pfqn::pfqn_ca(Lq, Nv, Zq).G);
    }
}

TEST_CASE("pfqn_lldsingle reproduces pfqn_gldsingle exactly on every rate shape") {
    // Rational arithmetic, so "==" is exact equality of the constant and not a
    // tolerance: the capped sweep must perform a SUBSET of the full sweep's
    // operations, never a rearrangement of them.
    Matrix<Rational> L(3, 1);
    L(0, 0) = Rational(7, 10);
    L(1, 0) = Rational(4, 10);
    L(2, 0) = Rational(9, 10);
    for (int N : {1, 2, 5, 9, 14}) {
        const std::size_t Nu = static_cast<std::size_t>(N);
        Matrix<Rational> mu(3, Nu);
        for (int k = 1; k <= N; ++k) {
            const std::size_t j = static_cast<std::size_t>(k - 1);
            mu(0, j) = Rational(1);                    // load independent, s = 1
            mu(1, j) = Rational(std::min(k, 3));       // multiserver, s = 3
            mu(2, j) = Rational(k);                    // never settles, s = N
        }
        CHECK(line::pfqn::pfqn_lldsingle(L, N, mu).G == line::pfqn::pfqn_gldsingle(L, N, mu).G);
    }
}

TEST_CASE("pfqn_lldsingle agrees with pfqn_gldsingle on a rate row that settles late") {
    // The threshold sits strictly inside the row, which is the case the cap is
    // for: below it every offset is distinct, above it they all collapse.
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(5, 10);
    L(1, 0) = Rational(8, 10);
    for (int N : {4, 7, 11}) {
        const std::size_t Nu = static_cast<std::size_t>(N);
        Matrix<Rational> mu(2, Nu);
        for (int k = 1; k <= N; ++k) {
            const std::size_t j = static_cast<std::size_t>(k - 1);
            mu(0, j) = Rational(std::min(k, 4));
            // 1, 2, 5/2, 5/2, ... : settles at 3 without a multiserver shape
            mu(1, j) = k == 1 ? Rational(1) : (k == 2 ? Rational(2) : Rational(5, 2));
        }
        CHECK(line::pfqn::pfqn_lldsingle(L, N, mu).G == line::pfqn::pfqn_gldsingle(L, N, mu).G);
    }
}

TEST_CASE("pfqn_mu_ms returns min(k,c) for a single station") {
    // With m = 1 the aggregate is the station itself, whose rate is min(k, c).
    const std::vector<Rational> mu = line::pfqn::pfqn_mu_ms<Rational>(7, 1, 3);
    REQUIRE(mu.size() == 7u);
    for (int k = 1; k <= 7; ++k) CHECK(mu[static_cast<std::size_t>(k - 1)] == Rational(k < 3 ? k : 3));
}

TEST_CASE("pfqn_mu_ms aggregate of m single-server stations reproduces the exact constant") {
    // The aggregate rate lattice of m identical stations must, when fed to
    // pfqn_gldsingle as one load-dependent station, give the same constant as
    // the convolution over the m stations kept apart.
    const int m = 3;
    for (int N : {1, 2, 4, 6}) {
        const std::vector<Rational> mu = line::pfqn::pfqn_mu_ms<Rational>(N, m, 1);
        Matrix<Rational> Lagg(1, 1), Lrep(static_cast<std::size_t>(m), 1);
        Lagg(0, 0) = Rational(7, 10);
        for (int i = 0; i < m; ++i) Lrep(static_cast<std::size_t>(i), 0) = Rational(7, 10);
        Matrix<Rational> muagg(1, static_cast<std::size_t>(N));
        for (int k = 1; k <= N; ++k)
            muagg(0, static_cast<std::size_t>(k - 1)) = mu[static_cast<std::size_t>(k - 1)];
        const std::vector<int> Nv = {N};
        CHECK(line::pfqn::pfqn_gldsingle(Lagg, N, muagg).G == line::pfqn::pfqn_ca(Lrep, Nv).G);
    }
}

TEST_CASE("pfqn_nc_sanitize factors the model without changing the constant") {
    // Sanitizing must be a change of variables: G_original = Gremaind * G_sanitized.
    Matrix<Rational> L(2, 4);
    Matrix<Rational> Z(1, 4);
    // class 0: normal; class 1: empty (N = 0); class 2: normal; class 3: normal
    const int lv[2][4] = {{7, 3, 0, 5}, {4, 6, 0, 2}};
    const int zv[4] = {3, 5, 9, 0};
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 4; ++r) L(i, r) = Rational(lv[i][r], 10);
    for (std::size_t r = 0; r < 4; ++r) Z(0, r) = Rational(zv[r], 10);
    // class 2 has zero demand everywhere, so it factors out as Z^N/N!.
    const std::vector<int> N = {2, 0, 3, 2};

    const line::pfqn::NcSanitizeResult<Rational> san = line::pfqn::pfqn_nc_sanitize(L, N, Z);
    const Rational gfull = line::pfqn::pfqn_ca(L, N, Z).G;
    const Rational gsan = line::pfqn::pfqn_ca(san.L, san.N, san.Z).G;
    CHECK(san.Gremaind * gsan == gfull);
    // The zero-think-time class must come first in the sanitized order.
    CHECK(san.classIndex.front() == 3u);
}

TEST_CASE("pfqn_ldmx_ec reduces to 1/(1-rho) at a load-independent station") {
    // With mu identically 1 the saturation level is b = 1, so the head sums are
    // empty and E(n) = 1/(1-Lo)^{n+1}, EC(n) = 1/(1-Lo).
    Matrix<Rational> D(2, 1), mu(2, 4);
    D(0, 0) = Rational(1, 4);
    D(1, 0) = Rational(1, 5);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 4; ++k) mu(i, k) = Rational(1);
    const std::vector<Rational> lambda = {Rational(1)};

    const line::pfqn::LdmxEcResult<Rational> ec = line::pfqn::pfqn_ldmx_ec(lambda, D, mu);
    for (std::size_t i = 0; i < 2; ++i) {
        const Rational rho = ec.Lo[i];
        const Rational inv = Rational(1) / (Rational(1) - rho);
        for (std::size_t n = 0; n <= 4; ++n) CHECK(ec.E(i, n) == line::num_pow_int(inv, static_cast<unsigned>(n + 1)));
        for (std::size_t n = 0; n < 4; ++n) CHECK(ec.EC(i, n) == inv);
    }
}

TEST_CASE("pfqn_ldmx_ec E-function matches its defining series") {
    // E_i(n) = sum_{n0>=0} C(n+n0,n0) Lo^{n0} prod_{j=n+1}^{n+n0} C(j), summed
    // to convergence in exact arithmetic against a truncation whose tail is
    // bounded well below the assertion.  The check is done in double because
    // the series is infinite; the ALGORITHM stays exact.
    Matrix<Rational> D(1, 1), mu(1, 5);
    D(0, 0) = Rational(1, 5);
    const int rates[5] = {1, 2, 3, 3, 3};  // saturates at b = 3
    for (std::size_t k = 0; k < 5; ++k) mu(0, k) = Rational(rates[k]);
    const std::vector<Rational> lambda = {Rational(1)};

    const line::pfqn::LdmxEcResult<Rational> ec = line::pfqn::pfqn_ldmx_ec(lambda, D, mu);
    const double Lo = static_cast<double>(ec.Lo[0]);
    for (int n = 0; n <= 5; ++n) {
        double series = 0.0;
        for (int n0 = 0; n0 <= 4000; ++n0) {
            double term = line::nck(n + n0, n0);
            for (int j = 1; j <= n0; ++j) {
                const int idx = n + j;
                term *= 1.0 / static_cast<double>(idx <= 5 ? rates[idx - 1] : rates[4]);
            }
            series += term * std::pow(Lo, n0);
        }
        CHECK(static_cast<double>(ec.E(0, static_cast<std::size_t>(n))) ==
              doctest::Approx(series).epsilon(1e-9));
    }
}

#include "line/api/pfqn/pfqn_nc.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/api/pfqn/pfqn_ncldmx.h"

TEST_CASE("pfqn_nc dispatches to methods that all return the identical exact G") {
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 3};
    const Rational gref = line::pfqn::pfqn_ca(L, N, Z).G;
    for (line::pfqn::NcMethod m :
         {line::pfqn::NcMethod::Ca, line::pfqn::NcMethod::Exact, line::pfqn::NcMethod::Mva}) {
        const auto r = line::pfqn::pfqn_nc(L, N, Z, m);
        INFO("method ", r.method);
        CHECK(r.G == gref);
    }
}

TEST_CASE("pfqn_nc refuses the unported methods by name instead of substituting") {
    // A multi-station model under 'default' routes to cub/le in the reference,
    // neither of which is in this tree. Answering with the convolution instead
    // would be indistinguishable from a correct answer until it is not.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 3};
    CHECK_THROWS_AS(line::pfqn::pfqn_nc(L, N, Z, line::pfqn::NcMethod::Default),
                    line::UnsupportedError);
    CHECK_THROWS_AS(line::pfqn::pfqn_nc_refuse("le"), line::UnsupportedError);
}

TEST_CASE("pfqn_nc default and comom solve the repairman model exactly") {
    const Matrix<Rational> L = rmat({{7, 4}}, 10);
    const Matrix<Rational> Z = rmat({{3, 9}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {3, 2}, {4, 3}}) {
        const Rational gref = line::pfqn::pfqn_ca(L, N, Z).G;
        CHECK(line::pfqn::pfqn_nc(L, N, Z, line::pfqn::NcMethod::Default).G == gref);
        CHECK(line::pfqn::pfqn_nc(L, N, Z, line::pfqn::NcMethod::Comom).G == gref);
    }
}

TEST_CASE("pfqn_nc reduces the degenerate models in closed form") {
    // Single station, no delay: the multinomial closed form.
    const Matrix<Rational> L1 = rmat({{5, 4}}, 10);
    const std::vector<int> N = {2, 3};
    CHECK(line::pfqn::pfqn_nc(L1, N, Matrix<Rational>(), line::pfqn::NcMethod::Ca).G ==
          line::pfqn::pfqn_ca(L1, N).G);
    // Three identical replicas, no delay: the multiset closed form.
    const Matrix<Rational> L3 = rmat({{5, 4}, {5, 4}, {5, 4}}, 10);
    CHECK(line::pfqn::pfqn_nc(L3, N, Matrix<Rational>(), line::pfqn::NcMethod::Ca).G ==
          line::pfqn::pfqn_ca(L3, N).G);
    // Delay only.
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    CHECK(line::pfqn::pfqn_nc(Matrix<Rational>(), N, Z, line::pfqn::NcMethod::Ca).G ==
          line::pfqn::pfqn_ca(Matrix<Rational>(), N, Z).G);
}

TEST_CASE("pfqn_nc with MVA satisfies Little's law and the utilization law") {
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {3, 2};
    const auto r = line::pfqn::pfqn_nc(L, N, Z, line::pfqn::NcMethod::Mva);
    REQUIRE(r.X.size() == 2u);
    for (std::size_t c = 0; c < 2; ++c) {
        // Population conservation: sum_i Q_ic + X_c Z_c = N_c, exactly.
        Rational q(0);
        for (std::size_t i = 0; i < 2; ++i) q += r.Q(i, c);
        CHECK(q + r.X[c] * Z(0, c) == Rational(N[c]));
        // Utilization law: U_ic = X_c L_ic, checked through Q at the bottleneck
        // by the residence-time identity Q = X * L * (1 + queue ahead) >= X * L.
        for (std::size_t i = 0; i < 2; ++i) CHECK(r.Q(i, c) >= r.X[c] * L(i, c));
    }
}

TEST_CASE("pfqn_ncld with unit rates reproduces the load-independent constant") {
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{2, 2}, {3, 2}, {1, 4}}) {
        int Nt = 0;
        for (int v : N) Nt += v;
        Matrix<Rational> mu(2, static_cast<std::size_t>(Nt));
        for (std::size_t i = 0; i < 2; ++i)
            for (int k = 0; k < Nt; ++k) mu(i, static_cast<std::size_t>(k)) = Rational(1);
        const auto r = line::pfqn::pfqn_ncld(L, N, Z, mu);
        INFO("method ", r.method);
        CHECK(r.G == line::pfqn::pfqn_ca(L, N, Z).G);
    }
}

TEST_CASE("pfqn_ncld single-class path agrees with the multiclass one") {
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(7, 10);
    L(1, 0) = Rational(4, 10);
    Matrix<Rational> Z(1, 1);
    Z(0, 0) = Rational(3, 10);
    for (int n : {2, 4, 7}) {
        Matrix<Rational> mu(2, static_cast<std::size_t>(n));
        for (int k = 1; k <= n; ++k) {
            mu(0, static_cast<std::size_t>(k - 1)) = Rational(k < 2 ? k : 2);
            mu(1, static_cast<std::size_t>(k - 1)) = Rational(1);
        }
        const std::vector<int> N = {n};
        const auto r = line::pfqn::pfqn_ncld(L, N, Z, mu);
        // The delay is folded into a row with rate lattice k, so the same model
        // written for pfqn_gld directly must give the same constant.
        Matrix<Rational> Lg(3, 1), mug(3, static_cast<std::size_t>(n));
        Lg(0, 0) = L(0, 0);
        Lg(1, 0) = L(1, 0);
        Lg(2, 0) = Z(0, 0);
        for (int k = 1; k <= n; ++k) {
            mug(0, static_cast<std::size_t>(k - 1)) = mu(0, static_cast<std::size_t>(k - 1));
            mug(1, static_cast<std::size_t>(k - 1)) = mu(1, static_cast<std::size_t>(k - 1));
            mug(2, static_cast<std::size_t>(k - 1)) = Rational(k);
        }
        CHECK(r.G == line::pfqn::pfqn_gld(Lg, N, mug).G);
    }
}

TEST_CASE("pfqn_ncldmx open prefactor is prod 1/(1-rho) with unit rates") {
    // Load-independent limit: E_i(0) = 1/(1-rho_i), and with no closed jobs the
    // closed-conditional constant is 1.
    Matrix<Rational> D(2, 2);
    D(0, 0) = Rational(1, 4);
    D(0, 1) = Rational(1, 2);
    D(1, 0) = Rational(1, 5);
    D(1, 1) = Rational(3, 10);
    const std::vector<int> N = {-1, 0};  // class 0 open, class 1 empty closed
    const std::vector<Rational> lambda = {Rational(1), Rational(0)};
    Matrix<Rational> mu(2, 2);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 2; ++k) mu(i, k) = Rational(1);

    const auto r = line::pfqn::pfqn_ncldmx(lambda, D, N, Matrix<Rational>(), mu);
    Rational expect(1);
    for (std::size_t i = 0; i < 2; ++i) expect /= (Rational(1) - D(i, 0));
    CHECK(r.Gopen == expect);
    CHECK(r.G == Rational(1));
}

TEST_CASE("pfqn_ncldmx with no open class is the purely closed constant") {
    // With lambda = 0 the effective capacity collapses to 1/mu, so the mixed
    // solver must return exactly what pfqn_ncld returns on the same model.
    Matrix<Rational> D(2, 2);
    D(0, 0) = Rational(5, 10);
    D(0, 1) = Rational(4, 10);
    D(1, 0) = Rational(3, 10);
    D(1, 1) = Rational(7, 10);
    const std::vector<int> N = {2, 2};
    const std::vector<Rational> lambda = {Rational(0), Rational(0)};
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(3, 10);
    Z(0, 1) = Rational(6, 10);
    Matrix<Rational> mu(2, 4);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 4; ++k) mu(i, k) = Rational(1);

    const auto r = line::pfqn::pfqn_ncldmx(lambda, D, N, Z, mu);
    CHECK(r.Gopen == Rational(1));
    CHECK(r.G == line::pfqn::pfqn_ca(D, N, Z).G);
}

TEST_CASE("pfqn_ncldmx mean measures reproduce pfqn_mvaldmx") {
    // pfqn_mvaldmx is exact on these models, so the normalizing-constant route
    // must return its throughputs and queue lengths, on single-server rates,
    // multiserver rates and a genuinely nonlinear rate row alike.
    struct Case {
        const char* name;
        std::vector<double> lambda;
        std::vector<std::vector<double> > D;
        std::vector<int> N;
        std::vector<double> Z;
        std::vector<std::vector<double> > mu;
    };
    const std::vector<Case> cases = {
        {"single server",
         {0.0, 0.2},
         {{1.0, 0.8}, {0.6, 0.4}},
         {5, -1},
         {0.0, 0.0},
         {{1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1}}},
        {"multiserver rates",
         {0.0, 0.2},
         {{1.0, 0.8}, {0.6, 0.4}},
         {5, -1},
         {0.0, 0.0},
         {{1, 2, 3, 3, 3, 3}, {1, 2, 2, 2, 2, 2}}},
        {"nonlinear rates",
         {0.0, 0.1},
         {{1.0, 0.9}, {0.7, 0.3}},
         {5, -1},
         {0.0, 0.0},
         {{1, 1.7, 2.2, 2.5, 2.6, 2.6}, {1, 1.4, 1.4, 1.4, 1.4, 1.4}}},
        {"think time",
         {0.0, 0.2},
         {{1.0, 0.8}, {0.6, 0.4}},
         {6, -1},
         {2.0, 0.0},
         {{1, 2, 2, 2, 2, 2, 2}, {1, 1, 1, 1, 1, 1, 1}}},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const std::size_t M = c.D.size();
        const std::size_t R = c.D[0].size();
        Matrix<double> D(M, R, 0.0), mu(M, c.mu[0].size(), 1.0), Z(1, R, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < R; ++r) D(i, r) = c.D[i][r];
            for (std::size_t k = 0; k < c.mu[i].size(); ++k) mu(i, k) = c.mu[i][k];
        }
        for (std::size_t r = 0; r < R; ++r) Z(0, r) = c.Z[r];
        const auto nc = line::pfqn::pfqn_ncldmx(c.lambda, D, c.N, Z, mu);
        const auto mx = line::pfqn::pfqn_mvaldmx(c.lambda, D, c.N, Z, mu);
        CHECK(nc.lG == doctest::Approx(mx.lG).epsilon(1e-9));
        for (std::size_t r = 0; r < R; ++r)
            CHECK(nc.XN[r] == doctest::Approx(mx.XN[r]).epsilon(1e-9));
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                CHECK(nc.QN(i, r) == doctest::Approx(mx.QN(i, r)).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_ncldmx single-server open queue length is the classical formula") {
    // b_i = 1 leaves no marginal correction, so the open queue length collapses
    // to lambda_r D_ir (1 + Q_i^closed) / (1 - rho_i).
    Matrix<double> D(2, 2, 0.0);
    D(0, 0) = 1.0; D(0, 1) = 0.8;
    D(1, 0) = 0.6; D(1, 1) = 0.4;
    const std::vector<int> N = {5, -1};
    const std::vector<double> lambda = {0.0, 0.2};
    Matrix<double> mu(2, 6, 1.0), Z(1, 2, 0.0);
    const auto nc = line::pfqn::pfqn_ncldmx(lambda, D, N, Z, mu);
    for (std::size_t i = 0; i < 2; ++i) {
        const double rho = lambda[1] * D(i, 1);
        const double expect = lambda[1] * D(i, 1) * (1.0 + nc.QN(i, 0)) / (1.0 - rho);
        CHECK(nc.QN(i, 1) == doctest::Approx(expect).epsilon(1e-10));
    }
}

#include "line/api/pfqn/pfqn_joint.h"
#include "line/api/pfqn/pfqn_lcfsqn_ca.h"
#include "line/api/pfqn/pfqn_lcfsqn_mva.h"
#include "line/api/pfqn/pfqn_lcfsqn_nc.h"
#include "line/api/pfqn/pfqn_pas_nc.h"
#include "line/api/pfqn/pfqn_ncoi.h"
#include "line/api/pfqn/pfqn_perm.h"

TEST_CASE("pfqn_perm reproduces the permanent of a small expanded matrix") {
    // With all multiplicities 1 the formula is the ordinary permanent, which is
    // enumerated here directly over the permutations.
    Matrix<Rational> A(3, 3);
    const int a[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) A(i, j) = Rational(a[i][j]);
    Rational expect(0);
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& p : perms) {
        Rational t(1);
        for (std::size_t i = 0; i < 3; ++i) t *= A(i, p[i]);
        expect += t;
    }
    CHECK(line::pfqn::pfqn_perm(A, std::vector<int>{1, 1, 1}) == expect);
}

TEST_CASE("pfqn_lcfsqn_nc closed form equals the pfqn_lcfsqn_ca recursion") {
    // The two share no code: one is a sum of permanents, the other a lattice
    // recursion. Agreement as an exact rational is the oracle for the two
    // reference defects the closed form had to be corrected for.
    const std::vector<Rational> alpha = {Rational(1, 2), Rational(1, 3), Rational(2, 5)};
    const std::vector<Rational> beta = {Rational(1, 5), Rational(3, 7), Rational(1, 4)};
    const std::vector<std::vector<int>> pops = {
        {1, 1, 1}, {2, 1, 1}, {2, 2, 1}, {3, 2, 1}, {1, 2, 3}, {2, 2, 2}, {4, 1, 1}};
    for (const std::vector<int>& N : pops) {
        INFO("N = ", N[0], ",", N[1], ",", N[2]);
        CHECK(line::pfqn::pfqn_lcfsqn_nc(alpha, beta, N) ==
              line::pfqn::pfqn_lcfsqn_ca(alpha, beta, N).G);
    }
}

TEST_CASE("pfqn_lcfsqn_mva base case matches the two-job closed form") {
    // One job of one class: it is at the LCFS station with probability
    // alpha/(alpha+beta), and the throughput is 1/(alpha+beta).
    const std::vector<Rational> alpha = {Rational(1, 2)};
    const std::vector<Rational> beta = {Rational(1, 5)};
    const auto r = line::pfqn::pfqn_lcfsqn_mva(alpha, beta, std::vector<int>{1});
    const Rational den = alpha[0] + beta[0];
    CHECK(r.Q(0, 0) == alpha[0] / den);
    CHECK(r.Q(1, 0) == beta[0] / den);
    CHECK(r.T_[0] == Rational(1) / den);
}

TEST_CASE("pfqn_lcfsqn_mva conserves the population and closes the utilizations") {
    const std::vector<Rational> alpha = {Rational(1, 2), Rational(1, 3)};
    const std::vector<Rational> beta = {Rational(1, 5), Rational(3, 7)};
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{1, 1}, {2, 1}, {2, 2}, {3, 2}}) {
        const auto r = line::pfqn::pfqn_lcfsqn_mva(alpha, beta, N);
        INFO("N = ", N[0], ",", N[1]);
        for (std::size_t c = 0; c < 2; ++c) {
            // Every class-c job is at one of the two stations: Q sums to N_c.
            CHECK(r.Q(0, c) + r.Q(1, c) == Rational(N[c]));
            // Utilization law U = X * S, by construction of U.
            CHECK(r.U(0, c) == r.T_[c] * alpha[c]);
            CHECK(r.U(1, c) == r.T_[c] * beta[c]);
        }
        // Back probabilities: each is a probability, and each station's total
        // cannot exceed one since at most one job is at the back of a queue.
        // (They do NOT sum to one across the two stations: a station can be
        // empty, and the two totals are normalized by different denominators
        // W_k, one per class.)
        Rational b0(0), b1(0);
        for (std::size_t c = 0; c < 2; ++c) {
            CHECK(r.B(0, c) > Rational(0));
            CHECK(r.B(1, c) > Rational(0));
            b0 += r.B(0, c);
            b1 += r.B(1, c);
            // Q(s,c) = B(s,c) + sum_r B(s,r) Q_{n-e_r}(s,c), a sum of
            // nonnegative terms, so the queue length dominates the back
            // probability at every station.
            CHECK(r.Q(0, c) >= r.B(0, c));
            CHECK(r.Q(1, c) >= r.B(1, c));
        }
        CHECK(b0 <= Rational(1));
        CHECK(b1 <= Rational(1));
    }
}

TEST_CASE("pfqn_ncoi with a constant unit rate is the load-independent constant") {
    // An OI station whose rate is 1 for every occupancy has balance function
    // Phi(n) = |n|!/prod n_r!, which is the single-server station with unit
    // demands, so the constant must equal pfqn_ca on that model.
    const std::vector<Rational> Z = {Rational(3, 10), Rational(6, 10)};
    std::vector<line::pfqn::OiRate<Rational>> mu(1);
    mu[0] = [](const std::vector<int>&) { return Rational(1); };
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{1, 1}, {2, 1}, {2, 2}, {3, 2}}) {
        Matrix<Rational> L(1, 2);
        L(0, 0) = Rational(1);
        L(0, 1) = Rational(1);
        Matrix<Rational> Zm(1, 2);
        Zm(0, 0) = Z[0];
        Zm(0, 1) = Z[1];
        CHECK(line::pfqn::pfqn_ncoi(Z, N, mu).G == line::pfqn::pfqn_ca(L, N, Zm).G);
    }
}

TEST_CASE("pfqn_ncoi with two constant-rate stations convolves them") {
    // Two OI stations at constant rates c1 and c2 are two single-server
    // stations with demands 1/c1 and 1/c2.
    const std::vector<Rational> Z = {Rational(3, 10), Rational(6, 10)};
    const Rational c1(2), c2(5, 2);
    std::vector<line::pfqn::OiRate<Rational>> mu(2);
    mu[0] = [c1](const std::vector<int>&) { return c1; };
    mu[1] = [c2](const std::vector<int>&) { return c2; };
    Matrix<Rational> L(2, 2), Zm(1, 2);
    for (std::size_t r = 0; r < 2; ++r) {
        L(0, r) = Rational(1) / c1;
        L(1, r) = Rational(1) / c2;
        Zm(0, r) = Z[r];
    }
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{1, 1}, {2, 1}, {2, 2}}) {
        CHECK(line::pfqn::pfqn_ncoi(Z, N, mu).G == line::pfqn::pfqn_ca(L, N, Zm).G);
    }
}

TEST_CASE("pfqn_ncoi on the count lattice equals the microstate pfqn_pas_nc") {
    // The macrostate convolution is legitimate only because an OI rank rate is
    // permutation-invariant: Phi(n) is the sum of the ordered-prefix weights of
    // every ordering of the multiset n. Exact rationals, so the two routines
    // must agree bit for bit, not merely to a tolerance.
    const std::vector<Rational> Z = {Rational(7, 10), Rational(13, 10)};
    std::vector<line::pfqn::OiRate<Rational>> mu(2);
    mu[0] = [](const std::vector<int>& n) {
        Rational s(0);
        if (n[0] > 0) s += Rational(4, 5);
        if (n[1] > 0) s += Rational(6, 5);
        return s;
    };
    mu[1] = [](const std::vector<int>& n) {
        Rational s(0);
        if (n[0] > 0) s += Rational(3, 2);
        if (n[1] > 0) s += Rational(3, 5);
        return s;
    };
    for (const std::vector<int>& N : std::vector<std::vector<int>>{{1, 1}, {2, 1}, {2, 2}, {3, 2}}) {
        CHECK(line::pfqn::pfqn_ncoi(Z, N, mu).G == line::pfqn::pfqn_pas_nc(Z, N, mu).G);
    }
}

TEST_CASE("pfqn_pas_nc restricts the constant to the communicating class") {
    // Four single-job classes on a two-station P&S cycle whose swap graph forces
    // a TOTAL placement order: the communicating class then holds exactly one
    // ordering per split, so G_C = sum_k (1/mu1)^k (1/mu2)^(4-k). The downstream
    // station traverses its chain in the opposite direction, hence the
    // transpose. MATLAB reference (SolverCTMC state space): 3.166240678.
    const std::size_t R = 4;
    line::pfqn::PlacementOrder P(R, std::vector<int>(R, 0));
    line::pfqn::PlacementOrder Pt(R, std::vector<int>(R, 0));
    for (std::size_t i = 0; i < R; ++i) {
        for (std::size_t j = i + 1; j < R; ++j) {
            P[i][j] = 1;
            Pt[j][i] = 1;
        }
    }
    const Rational mu1(1), mu2(13, 10);
    std::vector<line::pfqn::OiRate<Rational>> mu(2);
    mu[0] = [mu1](const std::vector<int>&) { return mu1; };
    mu[1] = [mu2](const std::vector<int>&) { return mu2; };
    const std::vector<Rational> Z(R, Rational(0));
    const std::vector<int> N(R, 1);

    Rational expect(0);
    for (std::size_t k = 0; k <= R; ++k) {
        Rational term(1);
        for (std::size_t a = 0; a < k; ++a) term /= mu1;
        for (std::size_t a = k; a < R; ++a) term /= mu2;
        expect += term;
    }
    std::vector<line::pfqn::PlacementOrder> prec;
    prec.push_back(P);
    prec.push_back(Pt);
    CHECK(line::pfqn::pfqn_pas_nc(Z, N, mu, prec).G == expect);

    // Without the order every ordering counts, so the plain OI constant is larger.
    CHECK(line::pfqn::pfqn_ncoi(Z, N, mu).G > expect);
}

TEST_CASE("pfqn_joint per-class probabilities sum to one") {
    // The strongest self-contained oracle for a joint distribution: summing the
    // weight over every reachable occupancy must give exactly the normalizing
    // constant, i.e. the probabilities must sum to exactly 1 as rationals.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 2};
    const Rational G = line::pfqn::pfqn_ca(L, N, Z).G;

    Rational tot(0);
    for (int a = 0; a <= N[0]; ++a)
        for (int b = 0; b <= N[1]; ++b)
            for (int c = 0; a + c <= N[0]; ++c)
                for (int d = 0; b + d <= N[1]; ++d) {
                    Matrix<int> n(2, 2);
                    n(0, 0) = a;
                    n(0, 1) = b;
                    n(1, 0) = c;
                    n(1, 1) = d;
                    tot += line::pfqn::pfqn_joint(n, L, N, Z, G);
                }
    CHECK(tot == Rational(1));
}

TEST_CASE("pfqn_joint total-queue-length probabilities sum to one") {
    // Same oracle for the permanent-based form, which is a genuinely different
    // computation: it sums the product-form weight over every per-class split.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 2};
    const Rational G = line::pfqn::pfqn_ca(L, N, Z).G;
    const int Nt = N[0] + N[1];

    Rational tot(0);
    for (int m0 = 0; m0 <= Nt; ++m0)
        for (int m1 = 0; m0 + m1 <= Nt; ++m1)
            tot += line::pfqn::pfqn_joint_total(std::vector<int>{m0, m1}, L, N, Z, G);
    CHECK(tot == Rational(1));
}

TEST_CASE("pfqn_joint total form is the marginalized per-class form") {
    // The two forms must agree state by state, not only in total.
    const Matrix<Rational> L = rmat({{5, 4}, {3, 7}}, 10);
    const Matrix<Rational> Z = rmat({{3, 6}}, 10);
    const std::vector<int> N = {2, 2};
    const Rational G = line::pfqn::pfqn_ca(L, N, Z).G;

    for (int m0 = 0; m0 <= 4; ++m0)
        for (int m1 = 0; m0 + m1 <= 4; ++m1) {
            Rational marg(0);
            for (int a = 0; a <= N[0] && a <= m0; ++a)
                for (int c = 0; c <= N[0] - a; ++c) {
                    const int b = m0 - a, d = m1 - c;
                    if (b < 0 || d < 0 || b > N[1] || d > N[1] - b) continue;
                    Matrix<int> n(2, 2);
                    n(0, 0) = a;
                    n(0, 1) = b;
                    n(1, 0) = c;
                    n(1, 1) = d;
                    marg += line::pfqn::pfqn_joint(n, L, N, Z, G);
                }
            INFO("m = ", m0, ",", m1);
            CHECK(line::pfqn::pfqn_joint_total(std::vector<int>{m0, m1}, L, N, Z, G) == marg);
        }
}
