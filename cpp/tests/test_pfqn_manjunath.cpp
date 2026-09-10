/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Manjunath-Sikdar transform for product-form queueing networks, against
 * oracles that do not come out of the implementation.
 *
 * TWO ORACLES. With no extra rows the transform computes the ordinary
 * closed-network normalizing constant, for which `pfqn_ca`'s convolution
 * recursion is an independent exact algorithm sharing no code; the two are
 * algebraic identities for the same sum, so agreement to 1e-13 is the correct
 * expectation and not a tolerance chosen to pass. With extra rows `pfqn_ca` has
 * nothing to say, and the oracle becomes `bcmp_enum` below, which sums the BCMP
 * product form over the enumerated state space and applies each row by direct
 * comparison -- the very enumeration the transform exists to avoid, so a
 * coefficient-domain defect cannot hide behind a shared traversal.
 *
 * THE EXACT PATH IS ITS OWN ORACLE. Under Arith::Exact the whole computation is
 * rational, so G is not merely close to the double answer, it is the answer the
 * double path is rounding. Comparing the two therefore bounds the accumulated
 * floating-point error rather than restating it.
 *
 * Ported from matlab/src/api/pfqn/test_pfqn_manjunath.m; the enumeration oracle
 * is its `bcmp_enum`.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_manjunath.h"
#include "line/num/number.h"
#include "line/util/error.h"

using line::Matrix;
using line::Rational;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_manjunath;
using line::pfqn::PfqnManjunathResult;

namespace {

Matrix<double> mk(const std::vector<std::vector<double>>& v) {
    if (v.empty()) return Matrix<double>();
    Matrix<double> m(v.size(), v[0].size());
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v[0].size(); ++j) m(i, j) = v[i][j];
    return m;
}

double factorial_d(int n) {
    double f = 1.0;
    for (int k = 2; k <= n; ++k) f *= k;
    return f;
}

/** Every way of splitting n indistinguishable jobs over k stations. */
void compositions(int n, int k, std::vector<std::vector<int>>& out) {
    out.clear();
    if (k <= 0) return;
    if (k == 1) {
        out.push_back(std::vector<int>(1, n));
        return;
    }
    for (int a = 0; a <= n; ++a) {
        std::vector<std::vector<int>> sub;
        compositions(n - a, k - 1, sub);
        for (std::size_t s = 0; s < sub.size(); ++s) {
            std::vector<int> row(1, a);
            row.insert(row.end(), sub[s].begin(), sub[s].end());
            out.push_back(row);
        }
    }
}

/**
 * Direct sum of the BCMP product form over the enumerated closed state space,
 * keeping the states that satisfy every extra row.
 */
double bcmp_enum(const Matrix<double>& L, const std::vector<int>& N, const Matrix<double>& Z,
                 const std::vector<std::vector<double>>& A, const std::vector<double>& b,
                 const std::string& sense) {
    const int M = static_cast<int>(L.empty() ? 0 : L.rows());
    const int Mz = static_cast<int>(Z.empty() ? 0 : Z.rows());
    const int S = M + Mz;
    const int R = static_cast<int>(N.size());
    std::vector<std::vector<std::vector<int>>> alloc(R);
    for (int r = 0; r < R; ++r) compositions(N[r], S, alloc[r]);

    std::vector<std::size_t> idx(R, 0);
    double g = 0.0;
    while (true) {
        // n(i,r), read column by column when the rows are applied.
        std::vector<int> n(S * R, 0);
        for (int r = 0; r < R; ++r)
            for (int i = 0; i < S; ++i) n[i + S * r] = alloc[r][idx[r]][i];

        bool ok = true;
        for (std::size_t j = 0; j < b.size() && ok; ++j) {
            double v = 0.0;
            for (int c = 0; c < S * R; ++c) v += A[j][c] * n[c];
            if (sense[j] == 'E')
                ok = std::fabs(v - b[j]) < 1e-9;
            else if (sense[j] == 'L')
                ok = v <= b[j] + 1e-9;
            else
                ok = v > b[j] + 1e-9;
        }
        if (ok) {
            double t = 1.0;
            for (int i = 0; i < M; ++i) {
                int ni = 0;
                for (int r = 0; r < R; ++r) ni += n[i + S * r];
                t *= factorial_d(ni);
                for (int r = 0; r < R; ++r)
                    t *= std::pow(L(i, r), n[i + S * r]) / factorial_d(n[i + S * r]);
            }
            for (int k = 0; k < Mz; ++k)
                for (int r = 0; r < R; ++r)
                    t *= std::pow(Z(k, r), n[M + k + S * r]) / factorial_d(n[M + k + S * r]);
            g += t;
        }

        int d = R - 1;
        for (; d >= 0; --d) {
            if (++idx[d] < alloc[d].size()) break;
            idx[d] = 0;
        }
        if (d < 0) break;
    }
    return g;
}

// The shared constrained fixture: three queueing stations and one delay, two
// classes, so the occupancy is 4-by-2 and every extra row has 8 coefficients.
const std::size_t S_C = 4;

std::vector<double> row_queue1() {
    std::vector<double> a(S_C * 2, 0.0);
    a[0] = 1;
    a[0 + S_C] = 1;
    return a;
}
std::vector<double> row_budget() {
    std::vector<double> a(S_C * 2, 0.0);
    a[0] = 2;
    a[1] = 1;
    a[0 + S_C] = 1;
    a[1 + S_C] = 3;
    return a;
}
std::vector<double> row_delay_class1() {
    std::vector<double> a(S_C * 2, 0.0);
    a[3] = 1;
    return a;
}

}  // namespace

TEST_CASE("pfqn_manjunath: with no extra rows it IS the closed-network constant") {
    struct C {
        std::vector<std::vector<double>> L;
        std::vector<int> N;
        std::vector<std::vector<double>> Z;
    };
    const std::vector<C> cs = {
        {{{1}, {2}}, {4}, {}},
        {{{1, 2}, {3, 1}}, {2, 3}, {}},
        {{{1, 2}, {3, 1}}, {2, 3}, {{0.5, 1.5}}},
        {{{0.4, 0.2}, {0.9, 0.7}, {0.1, 1.1}}, {3, 2}, {{1, 2}}},
        {{}, {2, 1}, {{1, 3}}},
        {{{1, 2}, {3, 1}}, {0, 0}, {{1, 1}}},
        {{{5, 1}, {1, 6}}, {6, 5}, {{2, 3}}},
    };
    for (std::size_t k = 0; k < cs.size(); ++k) {
        const Matrix<double> L = mk(cs[k].L), Z = mk(cs[k].Z);
        const double lG_ca = pfqn_ca<double>(L, cs[k].N, Z).lG;
        const double lG_mj = pfqn_manjunath<double>(L, cs[k].N, Z).lG;
        INFO("case ", k);
        CHECK(std::fabs(lG_ca - lG_mj) / std::max(1.0, std::fabs(lG_ca)) < 1e-13);
    }
}

TEST_CASE("pfqn_manjunath: extra rows match brute-force enumeration") {
    const Matrix<double> L = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Z = mk({{1, 2}});
    const std::vector<int> N = {3, 2};
    const std::vector<double> a1 = row_queue1(), a2 = row_budget(), a3 = row_delay_class1();

    struct Tst {
        std::vector<std::vector<double>> A;
        std::vector<long> b;
        std::string s;
    };
    const std::vector<Tst> ts = {
        {{a1}, {2}, "L"},        {{a1}, {2}, "E"},         {{a1}, {1}, "G"},
        {{a2}, {6}, "L"},        {{a1, a2}, {2, 6}, "LL"}, {{a1, a2}, {2, 6}, "EL"},
        {{a1, a2}, {1, 6}, "GL"},{{a1, a2}, {1, 5}, "GG"}, {{a1, a3}, {2, 1}, "LL"},
    };
    for (std::size_t k = 0; k < ts.size(); ++k) {
        std::vector<double> bd(ts[k].b.begin(), ts[k].b.end());
        const double want = bcmp_enum(L, N, Z, ts[k].A, bd, ts[k].s);
        const PfqnManjunathResult<double> got =
            pfqn_manjunath<double>(L, N, Z, mk(ts[k].A), ts[k].b, ts[k].s);
        INFO("row set ", ts[k].s, " case ", k);
        CHECK(std::fabs(want - got.G) / std::fabs(want) < 1e-12);
        CHECK(std::fabs(std::exp(got.lG) - got.G) / got.G < 1e-12);
    }
}

TEST_CASE("pfqn_manjunath: the exact path is the answer the double path rounds") {
    // Arith::Exact carries the whole series as rationals, so its G is not an
    // approximation of the double answer -- the double answer approximates it.
    Matrix<Rational> L(3, 2);
    const double lv[3][2] = {{1, 2}, {3, 1}, {0.5, 0.5}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) L(i, j) = Rational(lv[i][j]);
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(1);
    Z(0, 1) = Rational(2);
    const std::vector<int> N = {3, 2};

    const std::vector<double> a1 = row_queue1();
    Matrix<Rational> Ar(1, S_C * 2);
    for (std::size_t c = 0; c < S_C * 2; ++c) Ar(0, c) = Rational(a1[c]);

    const Matrix<double> Ld = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Zd = mk({{1, 2}});
    const Matrix<double> Ad = mk({a1});

    for (long cap = 0; cap <= 5; ++cap) {
        const PfqnManjunathResult<Rational> ex =
            pfqn_manjunath<Rational>(L, N, Z, Ar, std::vector<long>{cap}, std::string("L"));
        const PfqnManjunathResult<double> db =
            pfqn_manjunath<double>(Ld, N, Zd, Ad, std::vector<long>{cap}, std::string("L"));
        const double exd = static_cast<double>(ex.G);
        INFO("cap ", cap);
        CHECK(std::fabs(exd - db.G) / std::max(1e-300, std::fabs(exd)) < 1e-13);
        CHECK(ex.peak_states == db.peak_states);
    }
}

TEST_CASE("pfqn_manjunath: the population row is redundant, its complement empty") {
    const Matrix<double> L = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Z = mk({{1, 2}});
    const std::vector<int> N = {3, 2};
    const long total = N[0] + N[1];
    const double Gref = pfqn_manjunath<double>(L, N, Z).G;

    std::vector<double> ones(S_C * 2, 1.0);
    const Matrix<double> A = mk({ones});
    CHECK(std::fabs(pfqn_manjunath<double>(L, N, Z, A, std::vector<long>{total},
                                           std::string("E")).G - Gref) / Gref < 1e-13);
    CHECK(std::fabs(pfqn_manjunath<double>(L, N, Z, A, std::vector<long>{total},
                                           std::string("L")).G - Gref) / Gref < 1e-13);
    CHECK(pfqn_manjunath<double>(L, N, Z, A, std::vector<long>{total}, std::string("G")).G
          == 0.0);
}

TEST_CASE("pfqn_manjunath: trivial rows are decided rather than carried") {
    const Matrix<double> L = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Z = mk({{1, 2}});
    const std::vector<int> N = {3, 2};
    const double Gref = pfqn_manjunath<double>(L, N, Z).G;

    const Matrix<double> zero = mk({std::vector<double>(S_C * 2, 0.0)});
    CHECK(std::fabs(pfqn_manjunath<double>(L, N, Z, zero, std::vector<long>{0},
                                           std::string("E")).G - Gref) / Gref < 1e-13);
    CHECK(pfqn_manjunath<double>(L, N, Z, zero, std::vector<long>{3}, std::string("E")).G == 0.0);
    CHECK(pfqn_manjunath<double>(L, N, Z, zero, std::vector<long>{-1}, std::string("L")).G == 0.0);
    const Matrix<double> A = mk({row_queue1()});
    CHECK(std::fabs(pfqn_manjunath<double>(L, N, Z, A, std::vector<long>{-1},
                                           std::string("G")).G - Gref) / Gref < 1e-13);
}

TEST_CASE("pfqn_manjunath: refusals name the reason") {
    const Matrix<double> L = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Z = mk({{1, 2}});
    const std::vector<int> N = {3, 2};

    std::vector<double> frac(S_C * 2, 0.0);
    frac[0] = 0.5;
    CHECK_THROWS_AS(pfqn_manjunath<double>(L, N, Z, mk({frac}), std::vector<long>{1},
                                           std::string("L")),
                    line::InputError);
    std::vector<double> neg(S_C * 2, 0.0);
    neg[0] = -1.0;
    CHECK_THROWS_AS(pfqn_manjunath<double>(L, N, Z, mk({neg}), std::vector<long>{1},
                                           std::string("L")),
                    line::InputError);
    CHECK_THROWS_AS(pfqn_manjunath<double>(L, N, Z, mk({row_queue1()}), std::vector<long>{1},
                                           std::string("X")),
                    line::InputError);
    // A cap far below what one job needs is a legitimate empty set, not an error.
    CHECK(pfqn_manjunath<double>(L, {-1, 2}, Z).G == 0.0);
}

TEST_CASE("pfqn_manjunath: the peak is the class lattice when unconstrained") {
    const Matrix<double> L = mk({{1, 2}, {3, 1}, {0.5, 0.5}});
    const Matrix<double> Z = mk({{1, 2}});
    const std::vector<int> N = {3, 2};
    // No extra rows leaves only the class axes live, so the realised cost must be
    // exactly the lattice pfqn_ca walks.
    CHECK(pfqn_manjunath<double>(L, N, Z).peak_states ==
          static_cast<std::size_t>((N[0] + 1) * (N[1] + 1)));
}

// ---------------------------------------------------------------------------
// The per-class decomposition (PfqnManjunathOptions::stats)
// ---------------------------------------------------------------------------
// Reference instance: PS queue (demands 1, 2) + delay (think 2, 4), N = [4 4]
// both starting at the delay, with 2*n1 + 3*n2 <= 10 on the queue occupancy.
// The expectations are the stationary law of an INDEPENDENTLY built exact CTMC
// under HOLD truncation (a refused admission is a deleted transition), which
// agrees with the truncated product form to 8.2e-17, so these are not the
// routine restating itself.
TEST_CASE("pfqn_manjunath: the per-class decomposition matches the exact HOLD chain") {
    const Matrix<double> L = mk({{1.0, 2.0}});
    const Matrix<double> Z = mk({{2.0, 4.0}});
    const std::vector<int> N = {4, 4};
    std::vector<double> a(4, 0.0);
    a[0] = 2;
    a[2] = 3;
    line::pfqn::PfqnManjunathOptions o;
    o.stats = true;
    const PfqnManjunathResult<double> r =
        pfqn_manjunath<double>(L, N, Z, mk({a}), std::vector<long>{10}, std::string("L"), o);

    REQUIRE(r.has_stats);
    const double refQ[2] = {1.88612099644128, 1.50177935943061};
    const double refX[2] = {0.544483985765125, 0.224199288256228};
    const double refU[2] = {0.544483985765125, 0.448398576512456};
    const double refT[2] = {1.08896797153025, 0.896797153024911};
    const double refB[2] = {1.02491103202847, 1.60142348754448};
    const double refD[2] = {2.11387900355872, 2.49822064056939};
    for (int r2 = 0; r2 < 2; ++r2) {
        INFO("class ", r2 + 1);
        CHECK(std::fabs(r.Q[r2] - refQ[r2]) / refQ[r2] < 1e-12);
        CHECK(std::fabs(r.X[r2] - refX[r2]) / refX[r2] < 1e-12);
        CHECK(std::fabs(r.U[r2] - refU[r2]) / refU[r2] < 1e-12);
        CHECK(std::fabs(r.think[r2] - refT[r2]) / refT[r2] < 1e-12);
        CHECK(std::fabs(r.blocked[r2] - refB[r2]) / refB[r2] < 1e-12);
        CHECK(std::fabs(r.delay[r2] - refD[r2]) / refD[r2] < 1e-12);
        // A blocked job never leaves the delay, so nothing escapes the
        // accounting: queue + thinking + held is exactly N.
        CHECK(std::fabs(r.Q[r2] + r.think[r2] + r.blocked[r2] - N[r2]) < 1e-12);
        CHECK(std::fabs(r.delay[r2] - r.think[r2] - r.blocked[r2]) < 1e-12);
    }
}

TEST_CASE("pfqn_manjunath: unconstrained throughput is the classical G(N-e_r)/G(N)") {
    const Matrix<double> L = mk({{1.0, 2.0}});
    const Matrix<double> Z = mk({{2.0, 4.0}});
    const std::vector<int> N = {4, 4};
    line::pfqn::PfqnManjunathOptions o;
    o.stats = true;
    const PfqnManjunathResult<double> r = pfqn_manjunath<double>(L, N, Z, o);
    const double lG = pfqn_ca<double>(L, N, Z).lG;
    for (int c = 0; c < 2; ++c) {
        std::vector<int> Nr(N);
        Nr[c] -= 1;
        const double lGr = pfqn_ca<double>(L, Nr, Z).lG;
        INFO("class ", c + 1);
        CHECK(std::fabs(r.X[c] - std::exp(lGr - lG)) < 1e-13);
        // Nothing can be held when nothing is constrained.
        CHECK(std::fabs(r.blocked[c]) < 1e-12);
    }
}

TEST_CASE("pfqn_manjunath: the decomposition refuses every other configuration") {
    const Matrix<double> L = mk({{1.0, 2.0}});
    const Matrix<double> Z = mk({{2.0, 4.0}});
    const std::vector<int> N = {4, 4};
    line::pfqn::PfqnManjunathOptions o;
    o.stats = true;

    CHECK_THROWS_AS(pfqn_manjunath<double>(L, N, mk({{2.0, 4.0}, {1.0, 1.0}}),
                                           Matrix<double>(1, 6), std::vector<long>{10},
                                           std::string("L"), o),
                    line::InputError);
    // Two queueing stations: the delay->q1->q2->delay cycle makes the chain
    // irreversible, so Kelly truncation no longer holds.
    CHECK_THROWS_AS(pfqn_manjunath<double>(mk({{1.0, 2.0}, {1.0, 2.0}}), N, Z,
                                           Matrix<double>(1, 6), std::vector<long>{10},
                                           std::string("L"), o),
                    line::UnsupportedError);
    std::vector<double> ad(4, 0.0);
    ad[1] = 1;  // column 1 is (delay, class 1)
    CHECK_THROWS_AS(pfqn_manjunath<double>(L, N, Z, mk({ad}), std::vector<long>{10},
                                           std::string("L"), o),
                    line::InputError);
    // G is still produced when the decomposition is not asked for.
    line::pfqn::PfqnManjunathOptions plain;
    CHECK(pfqn_manjunath<double>(mk({{1.0, 2.0}, {1.0, 2.0}}), N, Z, plain).G > 0.0);
}
