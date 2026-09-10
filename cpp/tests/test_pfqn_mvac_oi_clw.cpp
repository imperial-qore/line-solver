/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_mvac, pfqn_mvacld, pfqn_mvaoi, pfqn_mvaoi_marg, pfqn_clw, pfqn_clw_lld.
 *
 * Expected values come from MATLAB (lineStart; then the routine under test),
 * quoted to 17 significant digits. The two exact-capable families additionally
 * assert BIT-IDENTITY at Rational against pfqn_mva, which reaches the same
 * numbers by the population recursion rather than by the chain recursion. The
 * CLW pair is inherently inexact (it inverts a generating function on a
 * contour) and is checked against MATLAB tightly and against the exact pfqn_ca
 * at the accuracy the method actually has.
 */
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mvac.h"
#include "line/api/pfqn/pfqn_mvacld.h"
#include "line/api/pfqn/pfqn_mvaoi.h"
#include "line/api/pfqn/pfqn_mvaoi_marg.h"
#include "line/api/pfqn/pfqn_ncoi.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::num_traits;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_clw;
using line::pfqn::pfqn_clw_lld;
using line::pfqn::pfqn_mva;
using line::pfqn::pfqn_mvac;
using line::pfqn::pfqn_mvacld;
using line::pfqn::pfqn_mvaoi;
using line::pfqn::pfqn_mvaoi_marg;

namespace {

double relerr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

template <class T>
Matrix<T> mat(std::size_t rows, std::size_t cols, const double* v) {
    Matrix<T> m(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j) m(i, j) = num_traits<T>::from_double(v[i * cols + j]);
    return m;
}

}  // namespace

TEST_CASE("pfqn_mvac reproduces MATLAB and pfqn_mva") {
    // MATLAB pfqn_mvac([1 2; 3 1], [2 3], [1 0.5]).
    const double Lv[4] = {1, 2, 3, 1};
    const double Zv[2] = {1, 0.5};
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 3;
    const double X[2] = {0.17316008294345051, 0.34404859902676904};
    const double Q[2][2] = {{0.56593274911690239, 1.7599433150288113},
                            {1.2609071679396473, 1.068032385457804}};
    const double U[2][2] = {{0.17316008294345051, 0.68809719805353808},
                            {0.51948024883035149, 0.34404859902676904}};
    const double C[2][2] = {{3.2682633289204479, 5.1153916045793206},
                            {7.2817426886508603, 3.1043067417772123}};

    SUBCASE("double") {
        const line::pfqn::MvacResult<double> r =
            pfqn_mvac(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zv));
        for (int s = 0; s < 2; ++s) {
            CHECK(r.X[s] == doctest::Approx(X[s]).epsilon(1e-13));
            for (int i = 0; i < 2; ++i) {
                CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-13));
                CHECK(r.U(i, s) == doctest::Approx(U[i][s]).epsilon(1e-13));
                CHECK(r.C(i, s) == doctest::Approx(C[i][s]).epsilon(1e-13));
            }
        }
    }
    SUBCASE("exact: identical to pfqn_mva as fractions") {
        const Matrix<Rational> L = mat<Rational>(2, 2, Lv);
        const Matrix<Rational> Z = mat<Rational>(1, 2, Zv);
        const line::pfqn::MvacResult<Rational> r = pfqn_mvac(L, N, Z);
        const line::pfqn::MvaResult<Rational> m = pfqn_mva(L, N, Z);
        for (std::size_t s = 0; s < 2; ++s) {
            CHECK(r.X[s] == m.XN[s]);
            for (std::size_t i = 0; i < 2; ++i) {
                CHECK(r.Q(i, s) == m.QN(i, s));
                CHECK(r.U(i, s) == m.UN(i, s));
            }
        }
    }
}

TEST_CASE("pfqn_mvac without a delay, three classes") {
    // MATLAB pfqn_mvac([1 2 0.5; 3 1 1.5; 0.5 0.5 2], [1 2 1], zeros(1,3)).
    const double Lv[9] = {1, 2, 0.5, 3, 1, 1.5, 0.5, 0.5, 2};
    const double Zv[3] = {0, 0, 0};
    std::vector<int> N(3);
    N[0] = 1;
    N[1] = 2;
    N[2] = 1;
    const double X[3] = {0.11003757380568976, 0.28502415458937197, 0.13929146537842191};
    const double Q[3][3] = {{0.25684380032206122, 1.0445517981749866, 0.17136339237788514},
                            {0.64975845410628019, 0.72302737520128824, 0.47222222222222221},
                            {0.093397745571658627, 0.23242082662372518, 0.35641438539989267}};
    const line::pfqn::MvacResult<double> r =
        pfqn_mvac(mat<double>(3, 3, Lv), N, mat<double>(1, 3, Zv));
    for (int s = 0; s < 3; ++s) {
        CHECK(r.X[s] == doctest::Approx(X[s]).epsilon(1e-13));
        for (int i = 0; i < 3; ++i) CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-13));
    }
    // Exact route: with no IS center every chain is resolved by a re-execution
    // of part 1, so this instance exercises the label-interchange path.
    const line::pfqn::MvacResult<Rational> re =
        pfqn_mvac(mat<Rational>(3, 3, Lv), N, mat<Rational>(1, 3, Zv));
    const line::pfqn::MvaResult<Rational> m =
        pfqn_mva(mat<Rational>(3, 3, Lv), N, mat<Rational>(1, 3, Zv));
    for (std::size_t s = 0; s < 3; ++s) {
        CHECK(re.X[s] == m.XN[s]);
        for (std::size_t i = 0; i < 3; ++i) CHECK(re.Q(i, s) == m.QN(i, s));
    }
}

TEST_CASE("pfqn_mvac reports empty classes as pfqn_mva does") {
    const double Lv[4] = {1, 2, 3, 1};
    const double Zv[2] = {1, 0.5};
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 0;
    const line::pfqn::MvacResult<double> r =
        pfqn_mvac(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zv));
    CHECK(r.X[1] == 0.0);
    CHECK(r.Q(0, 1) == 0.0);
    CHECK(r.C(0, 1) == doctest::Approx(2.0));  // the bare demand, as in pfqn_mva
    CHECK(r.C(1, 1) == doctest::Approx(1.0));
}

TEST_CASE("pfqn_mvacld reproduces MATLAB and pfqn_mvac at unit rates") {
    const double Lv[4] = {1, 2, 3, 1};
    const double Zv[2] = {1, 0.5};
    const double muv[8] = {1, 2, 2, 2, 1, 1, 1, 1};
    std::vector<int> N(2, 2);

    // MATLAB pfqn_mvacld([1 2; 3 1], [2 2], [1 0.5], [1 2 2 2; 1 1 1 1]).
    const double X[2] = {0.2014806980433633, 0.35536753040719199};
    const double Q[2][2] = {{0.23744050766790059, 0.76044420941300894},
                            {1.5610787942887363, 1.0618720253833951}};
    const double U[2] = {0.65388683236382861, 0.95980962453728191};
    const double C[2] = {8.9265091863517068, 5.1279761904761898};
    const double pij[2][5] = {{0.34611316763617139, 0.39555790586991013, 0.1853516657852988,
                               0.06028556319407722, 0.012691697514542572},
                              {0.040190375462718135, 0.11528291909042834, 0.25436277102062399,
                               0.36171337916446328, 0.22845055526176633}};

    SUBCASE("load dependent") {
        const line::pfqn::MvacldResult<double> r =
            pfqn_mvacld(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zv), mat<double>(2, 4, muv));
        for (int s = 0; s < 2; ++s) {
            CHECK(r.X[s] == doctest::Approx(X[s]).epsilon(1e-13));
            CHECK(r.C[s] == doctest::Approx(C[s]).epsilon(1e-13));
            for (int i = 0; i < 2; ++i) CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-13));
        }
        for (int i = 0; i < 2; ++i) {
            CHECK(r.U[i] == doctest::Approx(U[i]).epsilon(1e-13));
            double s = 0;
            for (int n = 0; n < 5; ++n) {
                CHECK(r.pij(i, n) == doctest::Approx(pij[i][n]).epsilon(1e-12));
                s += r.pij(i, n);
            }
            CHECK(s == doctest::Approx(1.0).epsilon(1e-13));  // (25) is self-normalizing
        }
    }

    SUBCASE("unit rates collapse to pfqn_mvac, exactly") {
        const double ones[8] = {1, 1, 1, 1, 1, 1, 1, 1};
        const line::pfqn::MvacldResult<Rational> r = pfqn_mvacld(
            mat<Rational>(2, 2, Lv), N, mat<Rational>(1, 2, Zv), mat<Rational>(2, 4, ones));
        const line::pfqn::MvacResult<Rational> c =
            pfqn_mvac(mat<Rational>(2, 2, Lv), N, mat<Rational>(1, 2, Zv));
        for (std::size_t s = 0; s < 2; ++s) {
            CHECK(r.X[s] == c.X[s]);
            for (std::size_t i = 0; i < 2; ++i) CHECK(r.Q(i, s) == c.Q(i, s));
        }
        // MATLAB pfqn_mvacld with mu = ones(2,4).
        CHECK(relerr(static_cast<double>(r.X[0]), 0.1966200254406687) < 1e-13);
        CHECK(relerr(static_cast<double>(r.X[1]), 0.28275486098491726) < 1e-13);
        CHECK(relerr(static_cast<double>(r.Q(0, 0)), 0.48628021079411232) < 1e-13);
        CHECK(relerr(static_cast<double>(r.Q(1, 1)), 0.78720697801199346) < 1e-13);
    }
}

TEST_CASE("pfqn_mvacld with two delay centers") {
    // MATLAB pfqn_mvacld([1 2; 3 1], [2 2], [0.5 0.5; 0.25 0.25], [1 2 2 2; 1 1 1 1]).
    const double Lv[4] = {1, 2, 3, 1};
    const double Zv[4] = {0.5, 0.5, 0.25, 0.25};
    const double muv[8] = {1, 2, 2, 2, 1, 1, 1, 1};
    std::vector<int> N(2, 2);
    const line::pfqn::MvacldResult<double> r =
        pfqn_mvacld(mat<double>(2, 2, Lv), N, mat<double>(2, 2, Zv), mat<double>(2, 4, muv));
    const double X[2] = {0.20833117351760216, 0.33860594195212784};
    const double Q[2][2] = {{0.24303708636205956, 0.72365656138727974},
                            {1.6007145334997388, 1.0223889821486245}};
    const double U[2] = {0.64024115505616197, 0.96359946250493433};
    const double C[2] = {8.8500995253406831, 5.1565708902496468};
    for (int s = 0; s < 2; ++s) {
        CHECK(r.X[s] == doctest::Approx(X[s]).epsilon(1e-13));
        CHECK(r.C[s] == doctest::Approx(C[s]).epsilon(1e-13));
        for (int i = 0; i < 2; ++i) CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-13));
    }
    for (int i = 0; i < 2; ++i) CHECK(r.U[i] == doctest::Approx(U[i]).epsilon(1e-13));
}

TEST_CASE("pfqn_mvaoi one OI station reproduces MATLAB") {
    // MATLAB pfqn_mvaoi([1.0 0.5], [2 1], {@(n) min(sum(n),2)}, [0.5 0.25]).
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 0.5;
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 1;
    std::vector<std::function<double(const std::vector<int>&)>> mu(1);
    mu[0] = [](const std::vector<int>& n) {
        int s = 0;
        for (int x : n) s += x;
        return static_cast<double>(s < 2 ? s : 2);
    };
    const double Dv[2] = {0.5, 0.25};
    const line::pfqn::MvaoiResult<double> r =
        pfqn_mvaoi(Z, N, mu, mat<double>(1, 2, Dv), true);

    CHECK(r.X[0] == doctest::Approx(0.71287128712871295).epsilon(1e-12));
    CHECK(r.X[1] == doctest::Approx(0.51485148514851486).epsilon(1e-12));
    CHECK(r.Qoi(0, 0) == doctest::Approx(0.79207920792079212).epsilon(1e-12));
    CHECK(r.Qoi(0, 1) == doctest::Approx(0.5544554455445545).epsilon(1e-12));
    CHECK(r.Qli(0, 0) == doctest::Approx(0.4950495049504951).epsilon(1e-12));
    CHECK(r.Qli(0, 1) == doctest::Approx(0.18811881188118812).epsilon(1e-12));
    CHECK(r.Qdelay[0] == doctest::Approx(0.71287128712871295).epsilon(1e-12));
    CHECK(r.Qdelay[1] == doctest::Approx(0.25742574257425743).epsilon(1e-12));
    CHECK(r.Soi(0, 0) == doctest::Approx(0.71287128712871295).epsilon(1e-12));
    CHECK(r.Soi(0, 1) == doctest::Approx(0.51485148514851486).epsilon(1e-12));
    // the population is conserved across delay, LI queue and OI station
    for (int s = 0; s < 2; ++s)
        CHECK(r.Qdelay[s] + r.Qli(0, s) + r.Qoi(0, s) ==
              doctest::Approx(static_cast<double>(N[s])).epsilon(1e-12));
}

TEST_CASE("pfqn_mvaoi two OI stations, and exactly at Rational") {
    // MATLAB pfqn_mvaoi([1.0 0.5],[2 1],{@(n) min(sum(n),2), @(n) sum(n)*0.5+1},[0.5 0.25]).
    const double X[2] = {0.56620465939250952, 0.37157180772633441};
    const double Qoi[2][2] = {{0.60159245060454136, 0.38926570333235033},
                              {0.46947803007962263, 0.29961663226186969}};
    const double Qli[2] = {0.36272485992332643, 0.1253317605426128};
    const double Soi[2][2] = {{0.56620465939250963, 0.37157180772633436},
                              {0.46947803007962247, 0.29961663226186963}};
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 1;
    const double Dv[2] = {0.5, 0.25};

    SUBCASE("double") {
        std::vector<double> Z(2);
        Z[0] = 1.0;
        Z[1] = 0.5;
        std::vector<std::function<double(const std::vector<int>&)>> mu(2);
        mu[0] = [](const std::vector<int>& n) {
            int s = 0;
            for (int x : n) s += x;
            return static_cast<double>(s < 2 ? s : 2);
        };
        mu[1] = [](const std::vector<int>& n) {
            int s = 0;
            for (int x : n) s += x;
            return 0.5 * s + 1.0;
        };
        const line::pfqn::MvaoiResult<double> r =
            pfqn_mvaoi(Z, N, mu, mat<double>(1, 2, Dv), true);
        for (int s = 0; s < 2; ++s) {
            CHECK(r.X[s] == doctest::Approx(X[s]).epsilon(1e-12));
            CHECK(r.Qli(0, s) == doctest::Approx(Qli[s]).epsilon(1e-12));
            for (int i = 0; i < 2; ++i) {
                CHECK(r.Qoi(i, s) == doctest::Approx(Qoi[i][s]).epsilon(1e-12));
                CHECK(r.Soi(i, s) == doctest::Approx(Soi[i][s]).epsilon(1e-12));
            }
        }
    }

    SUBCASE("exact") {
        // Every rate is rational, so the whole recursion stays in the field and
        // the result is a fraction that rounds to MATLAB's double.
        std::vector<Rational> Z(2);
        Z[0] = num_traits<Rational>::from_int(1);
        Z[1] = num_traits<Rational>::from_rational(1, 2);
        std::vector<std::function<Rational(const std::vector<int>&)>> mu(2);
        mu[0] = [](const std::vector<int>& n) {
            int s = 0;
            for (int x : n) s += x;
            return num_traits<Rational>::from_int(s < 2 ? s : 2);
        };
        // The trailing return type is LOAD-BEARING at Rational: cpp_rational has
        // expression templates ON, so `a + b` is an expression node holding
        // REFERENCES to the two temporaries. A deduced return type makes that node
        // the lambda's return type and std::function converts it to Rational after
        // the temporaries have died -- a use-after-free that Boost 1.74 tolerates
        // and 1.71 turns into a SIGSEGV. `-> Rational` converts inside the return
        // statement, while the operands are still alive.
        mu[1] = [](const std::vector<int>& n) -> Rational {
            int s = 0;
            for (int x : n) s += x;
            return num_traits<Rational>::from_rational(s, 2) + num_traits<Rational>::from_int(1);
        };
        const line::pfqn::MvaoiResult<Rational> r =
            pfqn_mvaoi(Z, N, mu, mat<Rational>(1, 2, Dv), true);
        for (std::size_t s = 0; s < 2; ++s) {
            CHECK(relerr(static_cast<double>(r.X[s]), X[s]) < 1e-13);
            for (std::size_t i = 0; i < 2; ++i)
                CHECK(relerr(static_cast<double>(r.Qoi(i, s)), Qoi[i][s]) < 1e-13);
        }
        // Population conservation holds EXACTLY, with no tolerance at all.
        for (std::size_t s = 0; s < 2; ++s) {
            Rational tot = r.Qdelay[s] + r.Qli(0, s);
            for (std::size_t i = 0; i < 2; ++i) tot += r.Qoi(i, s);
            CHECK(tot == num_traits<Rational>::from_int(N[s]));
        }
    }
}

TEST_CASE("pfqn_mvaoi_marg reaches the same numbers by the marginal route") {
    // MATLAB pfqn_mvaoi_marg([1.0 0.5; 0.5 0.25; 1 1], [2 1], [true false false],
    //                        {[], [], @(micro) min(numel(micro),2)}).
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 1;
    const double Dv[6] = {1.0, 0.5, 0.5, 0.25, 1, 1};
    std::vector<bool> isDelay(3, false);
    isDelay[0] = true;
    std::vector<std::function<double(const std::vector<int>&)>> mu(3);
    mu[2] = [](const std::vector<int>& micro) {
        return static_cast<double>(micro.size() < 2 ? micro.size() : 2);
    };
    const line::pfqn::MvaoiMargResult<double> r =
        pfqn_mvaoi_marg(mat<double>(3, 2, Dv), N, isDelay, mu);
    CHECK(r.X[0] == doctest::Approx(0.71287128712869696).epsilon(1e-11));
    CHECK(r.X[1] == doctest::Approx(0.51485148514852563).epsilon(1e-11));
    const double Q[3][2] = {{0.71287128712869696, 0.25742574257426282},
                            {0.495049504950484, 0.18811881188119209},
                            {0.79207920792078301, 0.55445544554455739}};
    for (int i = 0; i < 3; ++i)
        for (int s = 0; s < 2; ++s) CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-11));
    // Agreement with the mean-value form, which shares no code with this one.
    std::vector<double> Z(2);
    Z[0] = 1.0;
    Z[1] = 0.5;
    std::vector<std::function<double(const std::vector<int>&)>> muc(1);
    muc[0] = [](const std::vector<int>& n) {
        int s = 0;
        for (int x : n) s += x;
        return static_cast<double>(s < 2 ? s : 2);
    };
    const double Dli[2] = {0.5, 0.25};
    const line::pfqn::MvaoiResult<double> c =
        pfqn_mvaoi(Z, N, muc, mat<double>(1, 2, Dli), false);
    for (int s = 0; s < 2; ++s) {
        CHECK(relerr(r.X[s], c.X[s]) < 1e-12);
        CHECK(relerr(r.Q(2, s), c.Qoi(0, s)) < 1e-12);
        CHECK(relerr(r.Q(1, s), c.Qli(0, s)) < 1e-11);
    }
}

TEST_CASE("pfqn_mvaoi_marg with two OI stations") {
    // MATLAB pfqn_mvaoi_marg([1.0 0.5; 0.5 0.25; 0 0; 0 0], [2 1],
    //   [true false false false], {[],[],@(mi) min(numel(mi),2), @(mi) numel(mi)*0.5+1}).
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 1;
    const double Dv[8] = {1.0, 0.5, 0.5, 0.25, 0, 0, 0, 0};
    std::vector<bool> isDelay(4, false);
    isDelay[0] = true;
    std::vector<std::function<double(const std::vector<int>&)>> mu(4);
    mu[2] = [](const std::vector<int>& mi) {
        return static_cast<double>(mi.size() < 2 ? mi.size() : 2);
    };
    mu[3] = [](const std::vector<int>& mi) { return 0.5 * static_cast<double>(mi.size()) + 1.0; };
    const line::pfqn::MvaoiMargResult<double> r =
        pfqn_mvaoi_marg(mat<double>(4, 2, Dv), N, isDelay, mu);
    CHECK(r.X[0] == doctest::Approx(0.56620465939252973).epsilon(1e-11));
    CHECK(r.X[1] == doctest::Approx(0.37157180772635634).epsilon(1e-11));
    const double Q[4][2] = {{0.56620465939252973, 0.18578590386317817},
                            {0.36272485992334474, 0.12533176054261896},
                            {0.60159245060454791, 0.3892657033323858},
                            {0.4694780300796374, 0.29961663226189628}};
    for (int i = 0; i < 4; ++i)
        for (int s = 0; s < 2; ++s) CHECK(r.Q(i, s) == doctest::Approx(Q[i][s]).epsilon(1e-11));
}

/*
 * TOLERANCES FOR THE CLW PAIR. The inversion is an ALTERNATING sum of 2 l_j K_j
 * contour evaluations per chain, so its double path is roundoff dominated and
 * agreeing with MATLAB's doubles beyond about 1e-10 is not achievable by any
 * reordering of the same arithmetic. Measured, on the p = 2 model:
 *
 *   port (double)  0.11608708328781989
 *   MATLAB         0.11608708328228583   (4.8e-11 apart)
 *   port (Real50)  0.11608708333417295
 *   pfqn_ca exact  0.11608708333333334   (7.2e-12 apart)
 *
 * The high-precision instantiation lands on the exact constant to 7.2e-12,
 * which is the 10^-gamma_1 = 1e-11 aliasing floor the method is configured for.
 * That is the assertion that actually pins the algorithm: it says the port
 * computes the right integral, and that both double runs are merely rounding
 * it. The double checks below are therefore held at the level the arithmetic
 * supports and no tighter, and the exact-oracle checks carry the weight.
 */
TEST_CASE("pfqn_clw reproduces MATLAB and approximates pfqn_ca") {
    SUBCASE("p = 2") {
        // MATLAB pfqn_clw([0.1 0.2; 0.3 0.05], [2 3], [1.0 0.5]).
        const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
        const double Zv[2] = {1.0, 0.5};
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 3;
        std::vector<double> Z(Zv, Zv + 2);
        const line::pfqn::ClwResult<double> r = pfqn_clw(mat<double>(2, 2, Lv), N, Z);
        CHECK(relerr(r.G, 0.11608708328228583) < 1e-10);
        CHECK(relerr(r.lG, -2.1534146515728034) < 1e-10);
        // and the exact constant, to the aliasing the method admits
        const double ca = pfqn_ca(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zv)).G;
        CHECK(relerr(r.G, ca) < 1e-9);
    }
    SUBCASE("queue multiplicities") {
        // MATLAB pfqn_clw([0.1 0.2; 0.3 0.05], [2 3], [1.0 0.5], [2;1]).
        const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
        const double Zv[2] = {1.0, 0.5};
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 3;
        std::vector<double> Z(Zv, Zv + 2);
        std::vector<long> m(2);
        m[0] = 2;
        m[1] = 1;
        const line::pfqn::ClwResult<double> r = pfqn_clw(mat<double>(2, 2, Lv), N, Z, m);
        CHECK(relerr(r.G, 0.28466458343747597) < 1e-9);
        CHECK(relerr(r.lG, -1.2564436920237849) < 1e-9);
    }
    SUBCASE("p = 1") {
        // MATLAB pfqn_clw([0.1; 0.3; 0.05], 6, 1.5).
        const double Lv[3] = {0.1, 0.3, 0.05};
        std::vector<int> N(1, 6);
        std::vector<double> Z(1, 1.5);
        const double Zv[1] = {1.5};
        const line::pfqn::ClwResult<double> r = pfqn_clw(mat<double>(3, 1, Lv), N, Z);
        CHECK(relerr(r.G, 0.12740260937696071) < 1e-13);
        CHECK(relerr(r.lG, -2.0604030542880398) < 1e-13);
        const double ca = pfqn_ca(mat<double>(3, 1, Lv), N, mat<double>(1, 1, Zv)).G;
        CHECK(relerr(r.G, ca) < 1e-10);
    }
    SUBCASE("p = 3") {
        // MATLAB pfqn_clw([0.1 0.2 0.15; 0.3 0.05 0.1], [2 2 1], [1.0 0.5 0.2]).
        const double Lv[6] = {0.1, 0.2, 0.15, 0.3, 0.05, 0.1};
        const double Zv[3] = {1.0, 0.5, 0.2};
        std::vector<int> N(3);
        N[0] = 2;
        N[1] = 2;
        N[2] = 1;
        std::vector<double> Z(Zv, Zv + 3);
        const line::pfqn::ClwResult<double> r = pfqn_clw(mat<double>(2, 3, Lv), N, Z);
        // p = 3 doubles the nesting depth and the cancellation with it: the
        // reference and the port sit 4.1e-7 apart and 6.8e-7 / 2.6e-7 from the
        // exact constant respectively, the port being the closer of the two.
        CHECK(relerr(r.G, 0.23645858843155063) < 1e-6);
        CHECK(relerr(r.lG, -1.441982188151596) < 1e-6);
        const double ca = pfqn_ca(mat<double>(2, 3, Lv), N, mat<double>(1, 3, Zv)).G;
        CHECK(relerr(r.G, ca) < 1e-6);
        CHECK(relerr(r.G, ca) < relerr(0.23645858843155063, ca));
    }
    SUBCASE("a zero-population chain, where the reference returns NaN") {
        // REFERENCE DEFECT: MATLAB pfqn_clw divides by 2 l_j K_j with K_j = 0.
        // The port drops the chain, as pfqn_clw_lld already does, and returns
        // the p = 1 constant.
        const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
        const double Zv[2] = {1.0, 0.5};
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 0;
        std::vector<double> Z(Zv, Zv + 2);
        const line::pfqn::ClwResult<double> r = pfqn_clw(mat<double>(2, 2, Lv), N, Z);
        CHECK(std::isfinite(r.G));
        const double ca = pfqn_ca(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zv)).G;
        CHECK(relerr(r.G, ca) < 1e-9);
    }
    SUBCASE("real50 tracks the double result") {
        const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
        std::vector<int> N(2);
        N[0] = 2;
        N[1] = 3;
        std::vector<Real50> Z(2);
        Z[0] = num_traits<Real50>::from_double(1.0);
        Z[1] = num_traits<Real50>::from_double(0.5);
        const line::pfqn::ClwResult<Real50> r = pfqn_clw(mat<Real50>(2, 2, Lv), N, Z);
        // THE strong check: at 50 digits the roundoff disappears and what is
        // left is the exact constant to within the configured aliasing.
        const double Zd[2] = {1.0, 0.5};
        const double ca = pfqn_ca(mat<double>(2, 2, Lv), N, mat<double>(1, 2, Zd)).G;
        CHECK(relerr(static_cast<double>(r.G), ca) < 1e-10);
        CHECK(relerr(static_cast<double>(r.G), ca) < relerr(0.11608708328228583, ca));
    }
}

TEST_CASE("pfqn_clw_lld reproduces MATLAB") {
    const double Lv[4] = {0.1, 0.2, 0.3, 0.05};
    const double Zv[2] = {1.0, 0.5};
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 3;
    std::vector<double> Z(Zv, Zv + 2);

    SUBCASE("two-server first queue") {
        // MATLAB pfqn_clw_lld(L, N, Z, [1 2 2 2 2; 1 1 1 1 1]).
        const double muv[10] = {1, 2, 2, 2, 2, 1, 1, 1, 1, 1};
        const line::pfqn::ClwResult<double> r =
            pfqn_clw_lld(mat<double>(2, 2, Lv), N, Z, mat<double>(2, 5, muv));
        CHECK(relerr(r.G, 0.08338854169556939) < 1e-10);
        CHECK(relerr(r.lG, -2.4842443687976452) < 1e-10);
    }
    SUBCASE("three-server first queue") {
        // MATLAB pfqn_clw_lld(L, N, Z, [1 2 3 3 3; 1 1 1 1 1]).
        const double muv[10] = {1, 2, 3, 3, 3, 1, 1, 1, 1, 1};
        const line::pfqn::ClwResult<double> r =
            pfqn_clw_lld(mat<double>(2, 2, Lv), N, Z, mat<double>(2, 5, muv));
        CHECK(relerr(r.G, 0.081573217597482123) < 1e-10);
        CHECK(relerr(r.lG, -2.5062542866037738) < 1e-10);
    }
    SUBCASE("unit rates collapse to pfqn_clw") {
        const double ones[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        const line::pfqn::ClwResult<double> r =
            pfqn_clw_lld(mat<double>(2, 2, Lv), N, Z, mat<double>(2, 5, ones));
        const line::pfqn::ClwResult<double> c = pfqn_clw(mat<double>(2, 2, Lv), N, Z);
        CHECK(relerr(r.G, 0.11608708328228583) < 1e-10);
        // F_i(x) = 1/(1-x) is the load-independent case, so the LLD path must
        // reproduce pfqn_clw bit for bit, not merely closely.
        CHECK(r.G == c.G);
    }
}

TEST_CASE("pfqn_mvaoi weights the demand base case by the per-station visits") {
    // MATLAB pfqn_mvaoi([1.0 0.5], [2 1], {@(n) min(sum(n),2), @(n) 0.5*sum(n)+1},
    //                   [0.5 0.25], [1 0.6; 0.8 1.2]).
    // Visits enter only the N_r = 1 base case theta_{i,r} = v_{i,r}/mu_i(...);
    // the N_r >= 2 ratio case cancels them, which is why a wrong wiring shows up
    // on small populations rather than large ones.
    const double Zv[2] = {1.0, 0.5};
    std::vector<double> Z(Zv, Zv + 2);
    std::vector<int> N(2);
    N[0] = 2;
    N[1] = 1;
    std::vector<std::function<double(const std::vector<int>&)>> mu;
    mu.push_back([](const std::vector<int>& n) {
        return static_cast<double>(std::min(n[0] + n[1], 2));
    });
    mu.push_back([](const std::vector<int>& n) { return 0.5 * (n[0] + n[1]) + 1.0; });
    const double Dv[2] = {0.5, 0.25};
    const double vv[4] = {1.0, 0.6, 0.8, 1.2};

    const line::pfqn::MvaoiResult<double> r =
        pfqn_mvaoi(Z, N, mu, mat<double>(1, 2, Dv), mat<double>(2, 2, vv), true);
    CHECK(relerr(r.X[0], 0.59268901682165487) < 1e-12);
    CHECK(relerr(r.X[1], 0.41442152342388788) < 1e-12);
    CHECK(relerr(r.Qoi(0, 0), 0.61860836158061938) < 1e-12);
    CHECK(relerr(r.Qoi(0, 1), 0.26161258643381502) < 1e-12);
    CHECK(relerr(r.Qoi(1, 0), 0.40236038832938303) < 1e-12);
    CHECK(relerr(r.Qoi(1, 1), 0.38941223565868249) < 1e-12);
    CHECK(relerr(r.Qli(0, 0), 0.38634223326834305) < 1e-12);
    CHECK(relerr(r.Qli(0, 1), 0.14176441619555860) < 1e-12);

    SUBCASE("an empty visit matrix is unit visits") {
        const line::pfqn::MvaoiResult<double> a = pfqn_mvaoi(Z, N, mu, mat<double>(1, 2, Dv));
        const double ones[4] = {1, 1, 1, 1};
        const line::pfqn::MvaoiResult<double> b =
            pfqn_mvaoi(Z, N, mu, mat<double>(1, 2, Dv), mat<double>(2, 2, ones));
        for (std::size_t s = 0; s < a.X.size(); ++s) CHECK(a.X[s] == b.X[s]);
    }

    SUBCASE("without LI queues the throughput is the visit-weighted G ratio") {
        std::vector<std::function<double(const std::vector<int>&)>> mu1;
        mu1.push_back(mu[0]);
        const double v1[2] = {1.0, 0.6};
        const Matrix<double> vis = mat<double>(1, 2, v1);
        const line::pfqn::MvaoiResult<double> s =
            pfqn_mvaoi(Z, N, mu1, Matrix<double>(), vis, false);
        std::vector<line::pfqn::OiRate<double>> murate;
        murate.push_back(mu1[0]);
        const double G = line::pfqn::pfqn_ncoi(Z, N, murate, vis).G;
        for (std::size_t rr = 0; rr < 2; ++rr) {
            std::vector<int> Nm = N;
            Nm[rr] -= 1;
            CHECK(relerr(s.X[rr], line::pfqn::pfqn_ncoi(Z, Nm, murate, vis).G / G) < 1e-12);
        }
    }
}

/*
 * The two speed-ups of CLW Sections 2.4 and 3, as applied by pfqn_clw: Euler
 * summation of the inner sums, and dimension reduction by decomposition. The
 * expected values are the paper's own Tables I and III, which pin the algorithm
 * rather than the port; they agree across MATLAB, JAR, Python native and C++.
 * Mirrors line-test.git/test/testsAPI/test_pfqn_clw_acceleration.m.
 */
TEST_CASE("pfqn_clw applies the CLW accelerations") {
    const double LOG10 = std::log(10.0);

    // Example 8.1: p = 1, q' = 10 distinct queues of multiplicity 5
    Matrix<double> L81(10, 1);
    for (int i = 0; i < 10; ++i) L81(static_cast<std::size_t>(i), 0) = 0.1 * (i + 1);
    const std::vector<long> m5(10, 5);
    const std::vector<double> Z81(1, 5.0);

    SUBCASE("Euler summation reproduces Table I") {
        // without it the last row alone would cost 4e7 contour points
        const int K[8] = {2, 20, 200, 2000, 20000, 200000, 2000000, 20000000};
        const double mant[8] = {5.377500, 1.906584, 1.381312, 1.284918,
                                1.541538, 1.569301, 1.572100, 1.572380};
        const int expo[8] = {2, 13, 26, 31, 35, 39, 43, 47};
        for (int t = 0; t < 8; ++t) {
            const std::vector<int> N(1, K[t]);
            const line::pfqn::ClwResult<double> r = pfqn_clw(L81, N, Z81, m5);
            CHECK(std::abs(r.lG / LOG10 - (std::log10(mant[t]) + expo[t])) < 1e-6);
        }
    }

    SUBCASE("Euler summation does not move the answer") {
        // the order m is doubled until |E(m,n) - E(m,n+1)| settles, so the
        // acceleration must agree with the exact sum past its n+m = 31 threshold
        line::pfqn::ClwOptions off;
        off.euler = false;
        const int K[3] = {40, 100, 250};
        for (int t = 0; t < 3; ++t) {
            const std::vector<int> N(1, K[t]);
            const line::pfqn::ClwResult<double> a = pfqn_clw(L81, N, Z81, m5);
            const line::pfqn::ClwResult<double> b =
                pfqn_clw(L81, N, Z81, m5, off);
            CHECK(relerr(a.lG, b.lG) < 1e-9);
        }
    }

    // Examples 8.3 and 8.4: chain 1 is a hub visiting all ten queues, chain g+1
    // has queue g to itself, so the interdependence graph is a star and D = {1}
    // leaves ten single-variable components (the paper's Figure 1)
    Matrix<double> Ls(10, 11);
    std::vector<double> Zs(11, 0.0);
    Zs[0] = 50.0;
    for (int g = 1; g <= 10; ++g) {
        Ls(static_cast<std::size_t>(g - 1), 0) = 1 + 0.1 * g;
        Ls(static_cast<std::size_t>(g - 1), static_cast<std::size_t>(g)) = 0.1 * g;
        Zs[static_cast<std::size_t>(g)] = 5.0 * (g + 1) - 10.0;
    }

    SUBCASE("dimension reduction reaches Table III") {
        // eleven chains: unreachable without the reduction. Mantissa and
        // exponent are kept apart, since 1.937826e683 is not a double.
        const int K1[4] = {2, 20, 200, 2000};
        const double mant[4] = {1.235628, 7.503087, 5.970503, 1.937826};
        const int expo[4] = {25, 45, 129, 683};
        for (int t = 0; t < 4; ++t) {
            std::vector<int> N(11, 2);
            N[0] = K1[t];
            const line::pfqn::ClwResult<double> r = pfqn_clw(Ls, N, Zs, m5);
            CHECK(std::abs(r.lG / LOG10 - (std::log10(mant[t]) + expo[t])) < 1e-5);
        }
    }

    SUBCASE("Table III rows 5-8 under the paper's own scale tuning") {
        // page 956: beta_j in [0.8, 1.2] on the largest examples
        const int K1[4] = {2, 20, 200, 2000};
        const double b1[4] = {0.8, 0.8, 0.8, 0.95};
        const double mant[4] = {3.004462, 1.677866, 8.032122, 1.617153};
        const int expo[4] = {107, 133, 260, 926};
        for (int t = 0; t < 4; ++t) {
            std::vector<int> N(11, 0);
            N[0] = K1[t];
            for (int g = 1; g <= 10; ++g) N[static_cast<std::size_t>(g)] = 5 * g;
            line::pfqn::ClwOptions opt;
            opt.beta.assign(11, 1.0);
            opt.beta[0] = b1[t];
            const line::pfqn::ClwResult<double> r = pfqn_clw(Ls, N, Zs, m5, opt);
            CHECK(std::abs(r.lG / LOG10 - (std::log10(mant[t]) + expo[t])) < 1e-3);
        }
    }

    SUBCASE("the reduced and the full inversion agree with exact convolution") {
        // three leaves keep the undecomposed inversion affordable
        Matrix<double> L(3, 4);
        std::vector<double> Z(4, 0.0);
        Z[0] = 5.0;
        for (int g = 1; g <= 3; ++g) {
            L(static_cast<std::size_t>(g - 1), 0) = 1 + 0.1 * g;
            L(static_cast<std::size_t>(g - 1), static_cast<std::size_t>(g)) = 0.1 * g;
            Z[static_cast<std::size_t>(g)] = 5.0 * g - 5.0;
        }
        std::vector<int> N(4, 3);
        N[0] = 4;
        line::pfqn::ClwOptions off;
        off.dimred = false;
        const double Zr[4] = {5.0, 0.0, 5.0, 10.0};
        const double lca = pfqn_ca(L, N, mat<double>(1, 4, Zr)).lG;
        CHECK(relerr(pfqn_clw(L, N, Z).lG, lca) < 1e-7);
        CHECK(relerr(pfqn_clw(L, N, Z, std::vector<long>(), off).lG, lca) < 1e-6);
    }

    SUBCASE("dimension reduction is inert on a coupled model") {
        // every chain visits every queue, so no subset D reduces the dimension
        // and the classical path must be taken to the last bit
        const double Lv[6] = {0.1, 0.2, 0.15, 0.3, 0.05, 0.1};
        const Matrix<double> L = mat<double>(2, 3, Lv);
        std::vector<int> N(3);
        N[0] = 2; N[1] = 2; N[2] = 1;
        std::vector<double> Z(3);
        Z[0] = 1.0; Z[1] = 0.5; Z[2] = 0.2;
        line::pfqn::ClwOptions off;
        off.dimred = false;
        const double a = pfqn_clw(L, N, Z).lG;
        const double b = pfqn_clw(L, N, Z, std::vector<long>(), off).lG;
        CHECK(a == b);
        CHECK(relerr(a, -1.441982188151596) < 1e-6);
    }

    SUBCASE("a disconnected model needs no committed variable") {
        // the graph is already disconnected, so the reduction is exact with D
        // empty and the constant is the product of the two single-chain ones
        const double Lv[8] = {0.4, 0.0, 0.2, 0.0, 0.0, 0.5, 0.0, 0.3};
        std::vector<int> N(2);
        N[0] = 3; N[1] = 4;
        std::vector<double> Z(2);
        Z[0] = 1.0; Z[1] = 2.0;
        const double both = pfqn_clw(mat<double>(4, 2, Lv), N, Z).lG;
        const double L1v[2] = {0.4, 0.2}, L2v[2] = {0.5, 0.3};
        const double one =
            pfqn_clw(mat<double>(2, 1, L1v), std::vector<int>(1, 3), std::vector<double>(1, 1.0)).lG;
        const double two =
            pfqn_clw(mat<double>(2, 1, L2v), std::vector<int>(1, 4), std::vector<double>(1, 2.0)).lG;
        CHECK(relerr(both, one + two) < 1e-9);
    }
}
