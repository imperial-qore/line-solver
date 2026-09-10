/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The remaining QBD and MAP family: qbd_rap / qbd_raprap1, qbd_bmapbmap1,
 * qbd_setupdelayoff, the ETAQA departure processes, the level-dependent QBD,
 * mmap_compress, mmdp_isfeasible, the MAP derivative pair and the MAP/M/1-PS
 * sojourn family.
 *
 * ORACLES, in the order the task prefers them.
 *
 * (a) EXACTLY SOLVABLE COLLAPSES. A QBD whose blocks are 1 x 1 is an M/M/1
 *     queue and every quantity has a closed form: R = rho, QN = rho/(1-rho),
 *     UN = rho. A RAP that is a PH, a BMAP whose batch size is one, and a
 *     setup delay of unbounded rate all collapse the same way. The departure
 *     process of a stable M/M/1 is Poisson by Burke's theorem, which pins the
 *     joint moments of qbd_depproc_jointmom to 1/lambda, 2/lambda^2 and
 *     1/lambda^2 without any reference implementation being involved.
 *
 * (b) INVARIANTS. Every R is checked against its own defining equation, every
 *     G against A0 G^2 + A1 G + A2 = 0, and every stationary vector against
 *     the balance equations of the generator it came from. For the
 *     level-dependent QBD the block-tridiagonal generator is assembled
 *     explicitly and pi Q = 0 is asserted as an EXACT rational identity, which
 *     is a strictly stronger statement than agreement with MATLAB to 1e-10.
 *
 * (c) MATLAB. Where neither of the above pins the value -- the ETAQA
 *     descriptors and the MAP/M/1-PS sojourn curves -- the numbers below were
 *     produced by the reference and are quoted with their provenance.
 */
#include <cmath>
#include <utility>
#include <vector>

#include "doctest.h"
#include "line/api/mam/ldqbd.h"
#include "line/api/mam/map_joint_derivative.h"
#include "line/api/mam/map_m1ps.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmdp_isfeasible.h"
#include "line/api/mam/qbd_bmapbmap1.h"
#include "line/api/mam/qbd_depproc.h"
#include "line/api/mam/qbd_rap.h"
#include "line/api/mam/qbd_setupdelayoff.h"

using line::Matrix;
using line::Rational;
using line::mam::Map;
using line::mam::Mmap;
using namespace line::mam;

namespace {

/** 1 x 1 matrix holding v. */
template <class T>
Matrix<T> m1(const T& v) {
    Matrix<T> A(1, 1);
    A(0, 0) = v;
    return A;
}

/** ||A||_inf. */
template <class T>
T ninf(const Matrix<T>& A) {
    T w = line::num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i) {
        T s = line::num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < A.cols(); ++j) s += line::num_abs(T(A(i, j)));
        if (s > w) w = s;
    }
    return w;
}

/** Residual ||A0 + R A1 + R^2 A2||_inf of the rate equation. */
template <class T>
double r_residual(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2,
                  const Matrix<T>& R) {
    Matrix<T> res = A0;
    const Matrix<T> t1 = line::matmul(R, A1);
    const Matrix<T> t2 = line::matmul(line::matmul(R, R), A2);
    for (std::size_t i = 0; i < res.rows(); ++i)
        for (std::size_t j = 0; j < res.cols(); ++j) res(i, j) += t1(i, j) + t2(i, j);
    return line::num_traits<T>::to_double(ninf(res));
}

/** Residual ||A0 G^2 + A1 G + A2||_inf. */
template <class T>
double g_residual(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2,
                  const Matrix<T>& G) {
    Matrix<T> res = A2;
    const Matrix<T> t1 = line::matmul(A0, line::matmul(G, G));
    const Matrix<T> t2 = line::matmul(A1, G);
    for (std::size_t i = 0; i < res.rows(); ++i)
        for (std::size_t j = 0; j < res.cols(); ++j) res(i, j) += t1(i, j) + t2(i, j);
    return line::num_traits<T>::to_double(ninf(res));
}

/** The rank-one RAP QBD of Bean and Nielsen (2010), their running example. */
struct BeanNielsen {
    Matrix<double> A0, A1, A2, B0, B1;
};

BeanNielsen bean_nielsen(double g) {
    BeanNielsen b;
    const Matrix<double> A1{{-1.0, 0.0, 0.0}, {-2.0 / 3.0, -1.0, 1.0}, {2.0 / 3.0, -1.0, -1.0}};
    const Matrix<double> Da{{14.0 / 5.0, -9.0 / 10.0, -9.0 / 10.0},
                            {26.0 / 15.0, -8.0 / 15.0, -8.0 / 15.0},
                            {58.0 / 15.0, -19.0 / 15.0, -19.0 / 15.0}};
    Matrix<double> Ds(3, 3);
    const double u[3] = {1.0, 2.0 / 3.0, 4.0 / 3.0};
    const double v[3] = {3.0, -1.0, -1.0};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) Ds(i, j) = u[i] * v[j];
    b.A1 = A1;
    b.A0 = Matrix<double>(3, 3);
    b.A2 = Matrix<double>(3, 3);
    b.B1 = Matrix<double>(3, 3);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            b.A0(i, j) = g * Da(i, j);
            b.A2(i, j) = (1.0 - g) * Ds(i, j);
            b.B1(i, j) = g * A1(i, j);
        }
    b.B0 = b.A0;
    return b;
}

}  // namespace

// ---------------------------------------------------------------------------
// map_ccdf_derivative and map_jointpdf_derivative
// ---------------------------------------------------------------------------

TEST_CASE("MAP derivatives at zero collapse to the exponential closed form") {
    // For Exp(lambda), pie = 1, D0 = -lambda, D1 = lambda, so the i-th CCDF
    // derivative is (-lambda)^i and the joint density derivative of orders
    // (i_1..i_k) is prod_j (-lambda)^{i_j} lambda^k.
    const Rational lam(3, 2);
    const Map<Rational> m = map_exponential(lam);
    CHECK(map_ccdf_derivative(m, 0u) == Rational(1));
    CHECK(map_ccdf_derivative(m, 1u) == -lam);
    CHECK(map_ccdf_derivative(m, 2u) == lam * lam);
    CHECK(map_ccdf_derivative(m, 3u) == -lam * lam * lam);

    std::vector<unsigned> iset;
    CHECK(map_jointpdf_derivative(m, iset) == Rational(1));
    iset.push_back(0u);
    CHECK(map_jointpdf_derivative(m, iset) == lam);
    iset.push_back(1u);  // (0, 1)
    CHECK(map_jointpdf_derivative(m, iset) == lam * (-lam) * lam);
    std::vector<unsigned> two;
    two.push_back(2u);
    two.push_back(2u);
    CHECK(map_jointpdf_derivative(m, two) == lam * lam * lam * lam * lam * lam);
}

TEST_CASE("MAP CCDF derivatives satisfy the moment identity") {
    // Integrating the CCDF gives the mean, and more generally
    // E[T^k] = k! pie (-D0)^-k e; the derivative pie D0^i e is its inverse
    // counterpart. The identity that ties them is pie D0^0 e = 1 and, for an
    // order-2 MAP, pie D0 e = -pie D0 (-D0)^0 e = -1/E[T] only when the MAP is
    // a renewal exponential. What holds in general is that nu_1 is the
    // negative of the initial hazard rate, pie (-D0) e.
    Matrix<Rational> D0(2, 2), D1(2, 2);
    D0(0, 0) = Rational(-4);
    D0(0, 1) = Rational(1);
    D0(1, 0) = Rational(0);
    D0(1, 1) = Rational(-2);
    D1(0, 0) = Rational(2);
    D1(0, 1) = Rational(1);
    D1(1, 0) = Rational(1);
    D1(1, 1) = Rational(1);
    const Map<Rational> m{D0, D1};
    const std::vector<Rational> pie = map_pie(m);
    Rational hz(0);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) hz += pie[i] * D0(i, j);
    CHECK(map_ccdf_derivative(m, 1u) == hz);
    // The joint derivative of orders (0) is pie D1 e, which for a MAP is the
    // rate at which the embedded chain leaves pie, not 1; computing it
    // independently pins the product order of the implementation.
    std::vector<unsigned> zero_one(1, 0u);
    Rational direct(0);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) direct += pie[i] * D1(i, j);
    CHECK(map_jointpdf_derivative(m, zero_one) == direct);
    // Orders (1) instead insert D0 once: pie D0 D1 e.
    std::vector<unsigned> one_one(1, 1u);
    Rational direct1(0);
    {
        const std::vector<Rational> t = line::vecmul(line::vecmul(pie, D0), D1);
        for (const Rational& v : t) direct1 += v;
    }
    CHECK(map_jointpdf_derivative(m, one_one) == direct1);
}

// ---------------------------------------------------------------------------
// mmdp_isfeasible
// ---------------------------------------------------------------------------

TEST_CASE("mmdp_isfeasible accepts a generator with a diagonal rate matrix") {
    Matrix<Rational> Q(2, 2), R(2, 2, Rational(0));
    Q(0, 0) = Rational(-1);
    Q(0, 1) = Rational(1);
    Q(1, 0) = Rational(2);
    Q(1, 1) = Rational(-2);
    R(0, 0) = Rational(1, 2);
    R(1, 1) = Rational(3);
    CHECK(mmdp_isfeasible(Q, R, Rational(0)));

    SUBCASE("a nonzero row sum is rejected") {
        Matrix<Rational> Qb = Q;
        Qb(0, 1) = Rational(2);
        CHECK_FALSE(mmdp_isfeasible(Qb, R, Rational(0)));
    }
    SUBCASE("a negative off-diagonal is rejected") {
        Matrix<Rational> Qb = Q;
        Qb(0, 0) = Rational(1);
        Qb(0, 1) = Rational(-1);
        CHECK_FALSE(mmdp_isfeasible(Qb, R, Rational(0)));
    }
    SUBCASE("an off-diagonal entry in R is rejected") {
        Matrix<Rational> Rb = R;
        Rb(0, 1) = Rational(1, 1000);
        CHECK_FALSE(mmdp_isfeasible(Q, Rb, Rational(0)));
        // and accepted once the tolerance covers it, which is the whole point
        // of exposing the tolerance rather than hard-wiring 1e-10
        CHECK(mmdp_isfeasible(Q, Rb, Rational(1, 100)));
    }
    SUBCASE("a negative deterministic time is rejected") {
        Matrix<Rational> Rb = R;
        Rb(1, 1) = Rational(-1);
        CHECK_FALSE(mmdp_isfeasible(Q, Rb, Rational(0)));
    }
    SUBCASE("a non-square Q is rejected") {
        Matrix<Rational> Qb(2, 3, Rational(0));
        CHECK_FALSE(mmdp_isfeasible(Qb, R, Rational(0)));
    }
}

// ---------------------------------------------------------------------------
// map_compute_R and the MAP/M/1-PS family
// ---------------------------------------------------------------------------

TEST_CASE("map_compute_R collapses to rho for Poisson arrivals") {
    // C = -lambda, D = lambda is Poisson, and the rate equation
    // D + R(C - mu) + mu R^2 = 0 becomes mu R^2 - (lambda + mu) R + lambda = 0
    // whose root in [0,1) is rho = lambda/mu.
    const double lambda = 0.8, mu = 1.0;
    const Matrix<double> C = m1(-lambda), D = m1(lambda);
    const Matrix<double> R1 = map_compute_R(C, D, mu);
    const Matrix<double> R2 = map_compute_R_quadratic(C, D, mu);
    CHECK(R1(0, 0) == doctest::Approx(lambda / mu).epsilon(1e-9));
    CHECK(R2(0, 0) == doctest::Approx(lambda / mu).epsilon(1e-14));
    // The closed-form scalar branch of map_compute_R_quadratic is exact.
    CHECK(std::fabs(R2(0, 0) - 0.8) < 1e-15);
    CHECK(map_compute_R_residual(C, D, mu, R1) < 1e-9);
    CHECK(map_compute_R_residual(C, D, mu, R2) < 1e-15);
}

TEST_CASE("the two map_compute_R splittings reach the same fixed point") {
    // A genuine 2-phase MAP: the two iterations use different splittings and
    // different iteration caps, so agreement is evidence about the fixed point
    // rather than about the code being shared.
    Matrix<double> C(2, 2), D(2, 2);
    C(0, 0) = -3.0;
    C(0, 1) = 0.5;
    C(1, 0) = 0.2;
    C(1, 1) = -1.2;
    D(0, 0) = 2.0;
    D(0, 1) = 0.5;
    D(1, 0) = 0.6;
    D(1, 1) = 0.4;
    const double mu = 6.0;
    const Matrix<double> R1 = map_compute_R(C, D, mu);
    const Matrix<double> R2 = map_compute_R_quadratic(C, D, mu);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(R1(i, j) - R2(i, j)) < 1e-9);
    // MATLAB map_compute_R(C, D, 6.0) on this instance, to all 15 digits.
    CHECK(R1(0, 0) == doctest::Approx(0.288177805028462).epsilon(1e-14));
    CHECK(R1(0, 1) == doctest::Approx(0.128488861610385).epsilon(1e-14));
    CHECK(R1(1, 0) == doctest::Approx(0.0903110565855864).epsilon(1e-14));
    CHECK(R1(1, 1) == doctest::Approx(0.0763556100713033).epsilon(1e-14));
    CHECK(map_compute_R_residual(C, D, mu, R1) < 1e-9);
    CHECK(map_compute_R_residual(C, D, mu, R2) < 1e-9);
}

TEST_CASE("map_m1ps_h_recursive satisfies its own recursion exactly") {
    // The recursion is finite and rational, so it instantiates at Rational and
    // the residual is identically zero, not merely small.
    Matrix<Rational> C(2, 2), D(2, 2);
    C(0, 0) = Rational(-3);
    C(0, 1) = Rational(1, 2);
    C(1, 0) = Rational(1, 5);
    C(1, 1) = Rational(-6, 5);
    D(0, 0) = Rational(2);
    D(0, 1) = Rational(1, 2);
    D(1, 0) = Rational(3, 5);
    D(1, 1) = Rational(2, 5);
    const Rational mu(6);
    const std::size_t N = 5, K = 4;
    const std::vector<std::vector<std::vector<Rational>>> h =
        map_m1ps_h_recursive(C, D, mu, N, K);

    Rational theta(0);
    for (std::size_t i = 0; i < 2; ++i)
        if (line::num_abs(Rational(C(i, i))) > theta) theta = line::num_abs(Rational(C(i, i)));
    Matrix<Rational> tIC = C;
    for (std::size_t i = 0; i < 2; ++i) tIC(i, i) += theta;

    for (std::size_t n = 0; n <= N; ++n) {
        CHECK(h[n][0][0] == Rational(1));
        CHECK(h[n][0][1] == Rational(1));
    }
    for (std::size_t k = 0; k + 1 <= K; ++k)
        for (std::size_t n = 0; n <= N; ++n) {
            std::vector<Rational> rhs = line::mulvec(tIC, h[n][k]);
            if (n > 0) {
                const Rational c = Rational(static_cast<long>(n)) * mu /
                                   Rational(static_cast<long>(n + 1));
                for (std::size_t i = 0; i < 2; ++i) rhs[i] += c * h[n - 1][k][i];
            }
            if (n < N) {
                const std::vector<Rational> t = line::mulvec(D, h[n + 1][k]);
                for (std::size_t i = 0; i < 2; ++i) rhs[i] += t[i];
            }
            for (std::size_t i = 0; i < 2; ++i)
                CHECK(h[n][k + 1][i] == rhs[i] / (theta + mu));  // exact, no tolerance
        }
}

TEST_CASE("the M/M/1-PS sojourn CCDF integrates to the known mean") {
    // For M/M/1-PS the mean sojourn time is 1/(mu - lambda), independent of the
    // discipline (PS is work conserving and service is exponential). The CCDF
    // must start at 1, decrease, and integrate to that mean.
    const double lambda = 0.8, mu = 1.0;
    const Matrix<double> C = m1(-lambda), D = m1(lambda);
    std::vector<double> x;
    const double h = 0.05;
    const std::size_t np = 2401;  // 0 .. 120
    for (std::size_t k = 0; k < np; ++k) x.push_back(h * static_cast<double>(k));

    const MapM1psResult<double> r = map_m1ps_sojourn(C, D, mu, x);
    REQUIRE(r.w_bar.size() == np);
    CHECK(r.w_bar[0] == doctest::Approx(1.0).epsilon(1e-9));
    for (std::size_t k = 1; k < np; ++k) CHECK(r.w_bar[k] <= r.w_bar[k - 1] + 1e-12);

    double mean = 0.0;
    for (std::size_t k = 1; k < np; ++k) mean += 0.5 * h * (r.w_bar[k] + r.w_bar[k - 1]);
    CHECK(mean == doctest::Approx(1.0 / (mu - lambda)).epsilon(2e-3));

    // The other entry point is a different truncation rule and a different
    // splitting for R, so it is an independent computation of the same curve.
    const MapM1psResult<double> r2 = map_m1ps_cdfrespt(C, D, mu, x);
    for (std::size_t k = 0; k < np; ++k) CHECK(std::fabs(r.w_bar[k] - r2.w_bar[k]) < 1e-8);
}

TEST_CASE("map_m1ps_cdfrespt agrees with MATLAB and map_m1ps_sojourn does not") {
    // MATLAB map_m1ps_cdfrespt(-0.8, 0.8, 1, [0 .5 1 2 5 10]), all 15 digits.
    // The port reproduces it exactly.
    const Matrix<double> C = m1(-0.8), D = m1(0.8);
    std::vector<double> x;
    x.push_back(0.0);
    x.push_back(0.5);
    x.push_back(1.0);
    x.push_back(2.0);
    x.push_back(5.0);
    x.push_back(10.0);
    const double ref[6] = {0.999999999998497, 0.831840811326327, 0.710219465572534,
                           0.543535585559104, 0.294793542137958, 0.13817076928751};
    const MapM1psResult<double> r2 = map_m1ps_cdfrespt(C, D, 1.0, x);
    for (std::size_t k = 0; k < 6; ++k)
        CHECK(r2.w_bar[k] == doctest::Approx(ref[k]).epsilon(1e-12));

    // REFERENCE DEFECT. map_m1ps_sojourn.m OVERWRITES the rate matrix R with
    // the Poisson truncation upper index inside its per-point loop
    //     L = max(0, floor(mean_val - 10*sqrt(mean_val)));
    //     R = ceil(mean_val + 10*sqrt(mean_val));
    // and then keeps using R as the rate matrix in
    // weight = pi_0 * (R^n) * D. On this very instance MATLAB returns
    //   0.200000000379427  3.98171991111919e+148  1.51131679024013e+135
    //   7.45737240016127e+152  4.27105806216854e+178  2.93689832612765e+200
    // for a complementary distribution function, which must lie in [0, 1].
    // At x = 0 the branch mean_val = 0 sets R = 0, so only the n = 0 term
    // survives and the answer is pi_0 D e / lambda = 0.2 instead of 1; at
    // every later point R is an integer around 11 and the level weights grow
    // like 11^n. The port keeps the two names apart, so its map_m1ps_sojourn
    // returns the same curve as map_m1ps_cdfrespt to 1e-9, which is the only
    // reading under which the two documented-identical functions agree.
    const MapM1psResult<double> r1 = map_m1ps_sojourn(C, D, 1.0, x);
    for (std::size_t k = 0; k < 6; ++k) {
        CHECK(r1.w_bar[k] == doctest::Approx(ref[k]).epsilon(1e-9));
        CHECK(r1.w_bar[k] <= 1.0 + 1e-9);
        CHECK(r1.w_bar[k] >= 0.0);
    }
}

// ---------------------------------------------------------------------------
// qbd_rap and qbd_raprap1
// ---------------------------------------------------------------------------

TEST_CASE("qbd_rap on scalar blocks is the M/M/1 queue") {
    const double lambda = 0.6, mu = 1.0;
    const double rho = lambda / mu;
    const QbdRapResult<double> r =
        qbd_rap(m1(lambda), m1(-(lambda + mu)), m1(mu), m1(lambda), m1(-lambda), std::size_t(12));
    // G is the MINIMAL NONNEGATIVE solution, which for scalar blocks with
    // lambda < mu is the root 1 and not the root lambda/mu: the port returns
    // 0.99999999999997224442, nowhere near 0.6. It does not reach 1 to within
    // 1e-14, and neither does MATLAB, which returns 0.99999999999997202238 on
    // the same input with the IDENTICAL residual 1.121325e-14 -- both stop as
    // soon as the residual falls under qbd_rap's resTol = 1e-10 * blockScale,
    // one Newton step short of convergence. R agrees with MATLAB to all 20
    // printed digits (0.59999999999998987477) and QN to 1.4999999999999365.
    CHECK(std::fabs(r.G(0, 0) - 1.0) < 3e-14);
    CHECK(std::fabs(r.G(0, 0) - rho) > 0.39);  // the other root is not returned
    CHECK(r.R(0, 0) == doctest::Approx(0.59999999999998987477).epsilon(1e-16));
    CHECK(r.QN == doctest::Approx(1.4999999999999364952).epsilon(1e-15));
    CHECK(r.spr == doctest::Approx(rho).epsilon(1e-12));
    CHECK(r.QN == doctest::Approx(rho / (1.0 - rho)).epsilon(1e-12));
    CHECK(r.pi0[0] == doctest::Approx(1.0 - rho).epsilon(1e-12));
    for (std::size_t n = 0; n <= 12; ++n)
        CHECK(r.levelProb[n] == doctest::Approx((1.0 - rho) * std::pow(rho, static_cast<double>(n)))
                                    .epsilon(1e-12));
    // Cancellation noise, only the magnitude means anything: port 1.09912e-14, MATLAB 1.121325e-14
    CHECK(g_residual(m1(lambda), m1(-(lambda + mu)), m1(mu), r.G) < 1e-12);
    CHECK(r_residual(m1(lambda), m1(-(lambda + mu)), m1(mu), r.R) < 1e-12);
}

TEST_CASE("qbd_rap reproduces the Bean and Nielsen rank-one example") {
    const BeanNielsen b = bean_nielsen(0.25);
    const QbdRapResult<double> r = qbd_rap(b.A0, b.A1, b.A2, b.B0, b.B1, std::size_t(8));

    // The defining equations, which hold whatever the blocks are.
    CHECK(g_residual(b.A0, b.A1, b.A2, r.G) < 1e-12);
    CHECK(r_residual(b.A0, b.A1, b.A2, r.R) < 1e-12);
    // G is a stochastic-analogue matrix: G e = e.
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) s += r.G(i, j);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
    // The level distribution is a probability distribution with a geometric
    // tail of ratio Sp(R), and the closed-form mean agrees with the series.
    CHECK(r.spr < 1.0);
    double mass = 0.0, mean = 0.0;
    for (std::size_t n = 0; n < r.levelProb.size(); ++n) {
        mass += r.levelProb[n];
        mean += static_cast<double>(n) * r.levelProb[n];
    }
    CHECK(mass > 0.0);
    CHECK(mass <= 1.0 + 1e-12);
    CHECK(mean <= r.QN + 1e-12);
    // pi_0 (B1 + R A2) = 0, the boundary equation of Theorem 7.
    Matrix<double> V = b.B1;
    const Matrix<double> RA2 = line::matmul(r.R, b.A2);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) V(i, j) += RA2(i, j);
    const std::vector<double> bal = line::vecmul(r.pi0, V);
    for (std::size_t j = 0; j < 3; ++j) CHECK(std::fabs(bal[j]) < 1e-10);

    // MATLAB qbd_rap(g*Da, A1, (1-g)*Ds, g*Da, g*A1, 8) with g = 0.25, all 15
    // digits. Table 1 of Bean and Nielsen (2010).
    const double ref[9] = {0.673567977915804,    0.217506326201979,    0.0726190476190478,
                           0.0242044002249317,   0.00806816589133422,  0.00268938808906102,
                           0.000896462705376736, 0.000298820901641861, 9.96069672164603e-05};
    for (std::size_t n = 0; n < 9; ++n)
        CHECK(r.levelProb[n] == doctest::Approx(ref[n]).epsilon(1e-11));
    CHECK(r.QN == doctest::Approx(0.489817737501274).epsilon(1e-12));
    CHECK(r.spr == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    // pi_0 is NOT a probability vector: a RAP QBD has signed level vectors,
    // and only the level MARGINALS pi_n e are probabilities.
    CHECK(r.pi0[0] == doctest::Approx(1.17874396135266).epsilon(1e-12));
    CHECK(std::fabs(r.pi0[1]) < 1e-12);
    CHECK(r.pi0[2] == doctest::Approx(-0.505175983436854).epsilon(1e-12));
    // The rank-one closed form: every row of G is v/(v e) with v = [3 -1 -1].
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(r.G(i, 0) == doctest::Approx(3.0).epsilon(1e-12));
        CHECK(r.G(i, 1) == doctest::Approx(-1.0).epsilon(1e-12));
        CHECK(r.G(i, 2) == doctest::Approx(-1.0).epsilon(1e-12));
    }
    CHECK(r.R(0, 0) == doctest::Approx(0.744444444444445).epsilon(1e-12));
    CHECK(r.R(2, 2) == doctest::Approx(-0.427777777777778).epsilon(1e-12));
}

TEST_CASE("qbd_rap rejects a non-conservative or unstable model") {
    Matrix<double> A0 = m1(1.0), A1 = m1(-2.0), A2 = m1(1.0);
    // (A0+A1+A2) e = 0 holds, but Sp(R) = 1 exactly: null recurrent.
    CHECK_THROWS_AS(qbd_rap(A0, A1, A2, A0, m1(-1.0), std::size_t(3)), line::NumericError);
    // A non-conservative repeating block.
    CHECK_THROWS_AS(qbd_rap(m1(1.0), m1(-2.0), m1(0.5), m1(1.0), m1(-1.0), std::size_t(3)),
                    line::InputError);
}

TEST_CASE("qbd_raprap1 on two exponentials is the M/M/1 queue") {
    const double lambda = 0.5, mu = 1.25;
    const double rho = lambda / mu;
    const QbdRapRap1Result<double> r =
        qbd_raprap1(map_exponential(lambda), map_exponential(mu));
    CHECK(r.XN == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.UN == doctest::Approx(rho).epsilon(1e-9));
    CHECK(r.eta == doctest::Approx(rho).epsilon(1e-9));
    CHECK(r.QN == doctest::Approx(rho / (1.0 - rho)).epsilon(1e-6));
    CHECK(r.core.QN == doctest::Approx(rho / (1.0 - rho)).epsilon(1e-12));
    // MATLAB qbd_raprap1({[-0.5],[0.5]}, {[-1.25],[1.25]}): the truncated
    // series stops at the same level, so the value agrees to all 15 digits.
    CHECK(r.QN == doctest::Approx(0.666666665465703).epsilon(1e-14));
    // The truncated series can only UNDERSTATE the closed form.
    CHECK(r.QN <= r.core.QN + 1e-12);
    CHECK(r_residual(r.F, r.L, r.B, r.R) < 1e-12);
}

TEST_CASE("qbd_raprap1 rescales the service process to a target utilization") {
    const double lambda = 2.0, target = 0.4;
    const QbdRapRap1Result<double> r =
        qbd_raprap1(map_exponential(lambda), map_exponential(1.0), target);
    CHECK(r.UN == doctest::Approx(target).epsilon(1e-8));
    CHECK(r.XN == doctest::Approx(lambda).epsilon(1e-12));
    // MATLAB qbd_raprap1({[-2],[2]}, {[-1],[1]}, 0.4).
    CHECK(r.QN == doctest::Approx(0.666666665465703).epsilon(1e-14));
}

// ---------------------------------------------------------------------------
// qbd_bmapbmap1
// ---------------------------------------------------------------------------

TEST_CASE("qbd_bmapbmap1 with unit batches reproduces the MAP/MAP/1 blocks") {
    const Map<Rational> a = map_exponential(Rational(3, 5));
    const Map<Rational> s = map_exponential(Rational(1));
    std::vector<Rational> p1(1, Rational(1));
    const QbdBmapBmap1Blocks<Rational> b = qbd_bmapbmap1(a, p1, s);
    const QbdMapMap1Blocks<Rational> ref = qbd_mapmap1_blocks(a, s);
    REQUIRE(b.A1.size() == 1u);
    CHECK(b.A1[0](0, 0) == ref.F(0, 0));
    CHECK(b.A0(0, 0) == ref.L(0, 0));
    CHECK(b.Am1(0, 0) == ref.B(0, 0));
    CHECK(b.A0bar(0, 0) == ref.Lbar(0, 0));
    // Reference note: B0 is a Kronecker SUM with the identity, so it differs
    // from A0bar by an unconditional unit rate on every diagonal entry.
    CHECK(b.B0(0, 0) == ref.Lbar(0, 0) + Rational(1));
    CHECK(b.B1[0](0, 0) == b.A1[0](0, 0));
}

TEST_CASE("qbd_bmapbmap1 splits D1 by the batch-size law") {
    const Map<Rational> a = map_exponential(Rational(2));
    const Map<Rational> s = map_exponential(Rational(3));
    std::vector<Rational> p;
    p.push_back(Rational(1, 4));
    p.push_back(Rational(1, 2));
    p.push_back(Rational(1, 4));
    const QbdBmapBmap1Blocks<Rational> b = qbd_bmapbmap1(a, p, s);
    REQUIRE(b.A1.size() == 3u);
    Rational total(0);
    for (std::size_t k = 0; k < 3; ++k) total += b.A1[k](0, 0);
    // The up-blocks partition the arrival matrix exactly.
    CHECK(total == a.D1(0, 0));
    CHECK(b.A1[1](0, 0) == a.D1(0, 0) * Rational(1, 2));
}

// ---------------------------------------------------------------------------
// qbd_setupdelayoff
// ---------------------------------------------------------------------------

TEST_CASE("the Coxian phase reproduces its target rate and SCV") {
    // Every branch of the fit, checked through the moments of the PH it
    // defines rather than through its parameters.
    const double rates[2] = {2.0, 0.5};
    const double scvs[4] = {0.25, 0.75, 1.0, 3.0};
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 4; ++b) {
            const Matrix<double> T0 = coxian_phase_subgen(rates[a], scvs[b]);
            const std::size_t n = T0.rows();
            // Build the renewal MAP with entry vector [1 0 ...].
            Matrix<double> D1(n, n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < n; ++j) s += T0(i, j);
                D1(i, 0) = -s;
            }
            const Map<double> ph{T0, D1};
            CHECK(map_mean(ph) == doctest::Approx(1.0 / rates[a]).epsilon(1e-10));
            // The Erlang branch can only realize 1/n, so its SCV is the
            // attainable one below the target, never above.
            const double got = map_scv(ph);
            if (scvs[b] > 0.5)
                CHECK(got == doctest::Approx(scvs[b]).epsilon(1e-8));
            else
                CHECK(got <= scvs[b] + 1e-12);
        }
}

TEST_CASE("qbd_setupdelayoff collapses to M/M/1 when setup and delay-off vanish") {
    // With an instantaneous setup and an instantaneous delay-off the server is
    // always available on an arrival, so the queue is a plain M/M/1.
    const double lambda = 0.6, mu = 1.0;
    const double QN = qbd_setupdelayoff(lambda, mu, 1e7, 1.0, 1e7, 1.0);
    CHECK(QN == doctest::Approx(lambda / mu / (1.0 - lambda / mu)).epsilon(1e-5));
    CHECK(QN == doctest::Approx(1.50000005883156).epsilon(1e-13));  // MATLAB
}

TEST_CASE("qbd_setupdelayoff grows monotonically as the setup slows down") {
    const double lambda = 0.6, mu = 1.0;
    const double base = lambda / mu / (1.0 - lambda / mu);
    double prev = qbd_setupdelayoff(lambda, mu, 1e7, 1.0, 1.0, 1.0);
    CHECK(prev >= base - 1e-6);
    const double alphas[3] = {5.0, 1.0, 0.25};
    for (std::size_t k = 0; k < 3; ++k) {
        const double q = qbd_setupdelayoff(lambda, mu, alphas[k], 1.0, 1.0, 1.0);
        CHECK(q >= prev - 1e-9);
        prev = q;
    }
    // A non-exponential setup is admitted on every branch of the Coxian fit.
    CHECK(qbd_setupdelayoff(lambda, mu, 1.0, 0.3, 1.0, 2.0) > base);
    CHECK(qbd_setupdelayoff(lambda, mu, 1.0, 0.75, 1.0, 0.75) > base);

    // MATLAB qbd_setupdelayoff(0.6, 1, alpharate, alphascv, betarate, betascv),
    // all 15 digits: the Erlang, exponential and hyperexponential branches of
    // the phase fit and the exponential delay-off.
    CHECK(qbd_setupdelayoff(lambda, mu, 1.0, 1.0, 1.0, 1.0) ==
          doctest::Approx(1.93636363196485).epsilon(1e-13));
    CHECK(qbd_setupdelayoff(lambda, mu, 5.0, 1.0, 1.0, 1.0) ==
          doctest::Approx(1.57813953175353).epsilon(1e-13));
    CHECK(qbd_setupdelayoff(lambda, mu, 0.25, 1.0, 1.0, 1.0) ==
          doctest::Approx(3.53999999408965).epsilon(1e-13));
    CHECK(qbd_setupdelayoff(lambda, mu, 1.0, 0.3, 1.0, 2.0) ==
          doctest::Approx(1.89187499603791).epsilon(1e-13));
    CHECK(qbd_setupdelayoff(lambda, mu, 1.0, 0.75, 1.0, 0.75) ==
          doctest::Approx(1.90757237888086).epsilon(1e-13));
}

// ---------------------------------------------------------------------------
// ETAQA departure processes
// ---------------------------------------------------------------------------

TEST_CASE("qbd_depproc_jointmom reproduces Burke's theorem for M/M/1") {
    // The departure process of a stationary M/M/1 queue is Poisson of the same
    // rate, so inter-departure times are i.i.d. Exp(lambda):
    //   E[X] = 1/lambda, E[X^2] = 2/lambda^2, E[X_0 X_1] = 1/lambda^2.
    // This is an oracle no reference implementation is involved in.
    const double lambda = 0.6, mu = 1.0;
    std::vector<std::pair<unsigned, unsigned>> iset;
    iset.push_back(std::make_pair(1u, 0u));
    iset.push_back(std::make_pair(2u, 0u));
    iset.push_back(std::make_pair(1u, 1u));
    const std::vector<double> jm =
        qbd_depproc_jointmom(map_exponential(lambda), map_exponential(mu), iset);
    REQUIRE(jm.size() == 3u);
    CHECK(jm[0] == doctest::Approx(1.0 / lambda).epsilon(1e-9));
    CHECK(jm[1] == doctest::Approx(2.0 / (lambda * lambda)).epsilon(1e-9));
    CHECK(jm[2] == doctest::Approx(1.0 / (lambda * lambda)).epsilon(1e-9));
    // SCV 1 and lag-1 autocorrelation 0, the Poisson signature.
    const double scv = (jm[1] - jm[0] * jm[0]) / (jm[0] * jm[0]);
    const double acf = (jm[2] - jm[0] * jm[0]) / (jm[1] - jm[0] * jm[0]);
    CHECK(scv == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(std::fabs(acf) < 1e-8);

    // REFERENCE DEFECT, and the reason this test has no MATLAB column.
    // qbd_depproc_jointmom.m does
    //     pi = QBD_pi(B, L0, R);  v0 = pi(1,:);
    // but QBD_pi.m closes with pi = reshape(pi', 1, []), so it returns a FLAT
    // ROW vector of length (levels * m), never a (levels x m) matrix. pi(1,:)
    // is therefore the WHOLE level distribution, not the level-0 phase vector,
    // and z comes out of length 3*levels*m instead of 3*m. The very next
    // product, z * factorial(i) * (-M0)^(-i-1), then fails with "Incorrect
    // dimensions for matrix multiplication" at line 70. It is not specific to
    // scalar blocks: it was reproduced both with the M/M/1 instance above
    // (z of length 138 against a 3 x 3 M0) and with a 2-phase arrival MAP.
    // The function is unusable in MATLAB for any input for which the level
    // series runs past level 0.
}

TEST_CASE("the ETAQA departure descriptor is NOT a conservative generator") {
    // Reference defect 1, reproduced and measured. See the header of
    // qbd_depproc.h: the aggregate row keeps F G twice and Lhat is written one
    // block too early, so (D0 + D1) e picks up exactly F e = lambda in the two
    // affected block rows and is zero everywhere else.
    const double lambda = 0.6, mu = 1.0;
    const std::size_t n = 4;
    const Map<double> Dp = qbd_depproc_etaqa(map_exponential(lambda), map_exponential(mu), n);
    REQUIRE(Dp.D0.rows() == n + 1);
    CHECK(qbd_depproc_residual(Dp) == doctest::Approx(lambda).epsilon(1e-9));

    std::vector<double> rowsum(n + 1, 0.0);
    for (std::size_t i = 0; i <= n; ++i)
        for (std::size_t j = 0; j <= n; ++j) rowsum[i] += Dp.D0(i, j) + Dp.D1(i, j);
    for (std::size_t i = 0; i + 2 < n + 1; ++i) CHECK(std::fabs(rowsum[i]) < 1e-12);
    CHECK(rowsum[n - 1] == doctest::Approx(lambda).epsilon(1e-9));
    CHECK(rowsum[n] == doctest::Approx(lambda).epsilon(1e-9));

    // The block structure that IS as documented: the boundary row is [L0, F],
    // every explicit level has L on the diagonal and B one block down.
    CHECK(Dp.D0(0, 0) == doctest::Approx(-lambda));
    CHECK(Dp.D0(0, 1) == doctest::Approx(lambda));
    CHECK(Dp.D0(1, 1) == doctest::Approx(-(lambda + mu)));
    CHECK(Dp.D1(1, 0) == doctest::Approx(mu));
    CHECK(Dp.D1(n, n - 1) == doctest::Approx(mu + lambda));  // Bbar = B + F G, G = 1
    CHECK(Dp.D1(n, n) == doctest::Approx(lambda));           // Bhat = F G
}

TEST_CASE("the PS ETAQA descriptor splits the completions by the sharing factor") {
    const double lambda = 0.6, mu = 1.0;
    const std::size_t n = 4;
    const Map<double> Dp = qbd_depproc_etaqa_ps(map_exponential(lambda), map_exponential(mu), n);
    REQUIRE(Dp.D0.rows() == n + 1);
    // At level j a completion is a departure with probability 1/j.
    for (std::size_t j = 1; j + 1 <= n; ++j) {
        CHECK(Dp.D1(j, j - 1) == doctest::Approx(mu / static_cast<double>(j)));
        CHECK(Dp.D0(j, j - 1) == doctest::Approx(mu * (1.0 - 1.0 / static_cast<double>(j))));
        CHECK(Dp.D0(j, j - 1) + Dp.D1(j, j - 1) == doctest::Approx(mu));
    }
    // Reference defect 2: the aggregate down-block is NOT split, it appears in
    // full in both D0 and D1, so the sum is 2 Bbar there rather than Bbar.
    CHECK(Dp.D0(n, n - 1) == doctest::Approx(mu + lambda));
    CHECK(Dp.D1(n, n - 1) == doctest::Approx(mu + lambda));
    // MATLAB qbd_depproc_etaqa_ps: the aggregate local block is
    // L + Bhat (n-1)/n = -1.6 + 0.6*3/4 = -1.15, the aggregate departure block
    // is Bhat/n = 0.15, and the row sums are 0, 0, 0, 0.6, 2.2.
    CHECK(Dp.D0(n, n) == doctest::Approx(-1.15).epsilon(1e-12));
    CHECK(Dp.D1(n, n) == doctest::Approx(0.15).epsilon(1e-12));
    CHECK(qbd_depproc_residual(Dp) == doctest::Approx(2.2).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// Level-dependent QBD
// ---------------------------------------------------------------------------

namespace {

/** Assemble the block-tridiagonal generator of an LD-QBD, levels 0..N. */
template <class T>
Matrix<T> ldqbd_generator(const std::vector<Matrix<T>>& q0, const std::vector<Matrix<T>>& q1,
                          const std::vector<Matrix<T>>& q2) {
    const std::size_t N = q1.size() - 1;
    std::vector<std::size_t> off(N + 2, 0);
    for (std::size_t n = 0; n <= N; ++n) off[n + 1] = off[n] + q1[n].rows();
    const std::size_t dim = off[N + 1];
    Matrix<T> Q(dim, dim, line::num_traits<T>::from_int(0));
    for (std::size_t n = 0; n <= N; ++n) {
        for (std::size_t i = 0; i < q1[n].rows(); ++i)
            for (std::size_t j = 0; j < q1[n].cols(); ++j) Q(off[n] + i, off[n] + j) = q1[n](i, j);
        if (n < N)
            for (std::size_t i = 0; i < q0[n].rows(); ++i)
                for (std::size_t j = 0; j < q0[n].cols(); ++j)
                    Q(off[n] + i, off[n + 1] + j) = q0[n](i, j);
        if (n > 0)
            for (std::size_t i = 0; i < q2[n].rows(); ++i)
                for (std::size_t j = 0; j < q2[n].cols(); ++j)
                    Q(off[n] + i, off[n - 1] + j) = q2[n](i, j);
    }
    return Q;
}

/** Flatten the per-level stationary vectors into one vector over all states. */
template <class T>
std::vector<T> flatten(const LdqbdPi<T>& p) {
    std::vector<T> v;
    for (std::size_t n = 0; n < p.pi_level.size(); ++n)
        for (const T& x : p.pi_level[n]) v.push_back(x);
    return v;
}

}  // namespace

TEST_CASE("ldqbd solves a finite M/M/1/K exactly in rational arithmetic") {
    // Scalar levels: Q0 = lambda, Q2 = mu, Q1 = -(lambda + mu) inside,
    // -lambda at level 0 and -mu at level K. Every R^(n) is rho and the
    // stationary law is the truncated geometric, both EXACTLY.
    const Rational lambda(3, 5), mu(1);
    const std::size_t K = 6;
    std::vector<Matrix<Rational>> q0, q1, q2;
    for (std::size_t n = 0; n <= K; ++n) {
        Rational loc(0);
        if (n == 0)
            loc = -lambda;
        else if (n == K)
            loc = -mu;
        else
            loc = -(lambda + mu);
        q1.push_back(m1(loc));
        q0.push_back(m1(lambda));
        q2.push_back(m1(mu));
    }
    const LdqbdResult<Rational> r = ldqbd(q0, q1, q2);
    const Rational rho = lambda / mu;
    for (std::size_t n = 1; n <= K; ++n) {
        REQUIRE(r.R[n].rows() == 1u);
        CHECK(r.R[n](0, 0) == rho);  // exact
    }
    Rational norm(0);
    for (std::size_t n = 0; n <= K; ++n) norm += line::num_pow_int(rho, static_cast<unsigned>(n));
    for (std::size_t n = 0; n <= K; ++n)
        CHECK(r.pi.pi[n] == line::num_pow_int(rho, static_cast<unsigned>(n)) / norm);  // exact

    // MATLAB ldqbd on the same chain: R^(n) = 0.6 for every n and
    // pi = [0.411519924148648 0.246911954489189 0.148147172693513
    //       0.0888883036161079 0.0533329821696647 0.0319997893017988
    //       0.0191998735810793], which is the rational answer above rounded.
    const double refpi[7] = {0.411519924148648,  0.246911954489189, 0.148147172693513,
                             0.0888883036161079, 0.0533329821696647, 0.0319997893017988,
                             0.0191998735810793};
    for (std::size_t n = 0; n <= K; ++n)
        CHECK(line::num_traits<Rational>::to_double(r.pi.pi[n]) ==
              doctest::Approx(refpi[n]).epsilon(1e-12));

    // pi Q = 0 on the assembled generator, exactly.
    const Matrix<Rational> Q = ldqbd_generator(q0, q1, q2);
    const std::vector<Rational> bal = line::vecmul(flatten(r.pi), Q);
    for (std::size_t j = 0; j < bal.size(); ++j) CHECK(bal[j] == Rational(0));
}

TEST_CASE("ldqbd solves a level-dependent M/M/c/K exactly") {
    // Service rate min(n, c) mu, so the down-blocks genuinely depend on the
    // level and R is no longer constant. The oracle is the Erlang-B / Erlang-C
    // birth-death product form, plus pi Q = 0.
    const Rational lambda(4), mu(1);
    const std::size_t c = 3, K = 8;
    std::vector<Matrix<Rational>> q0, q1, q2;
    for (std::size_t n = 0; n <= K; ++n) {
        const Rational down =
            Rational(static_cast<long>(n < c ? n : c)) * mu;
        const Rational up = (n == K) ? Rational(0) : lambda;
        q1.push_back(m1(Rational(-(up + down))));
        q0.push_back(m1(lambda));
        q2.push_back(m1(down));
    }
    const LdqbdResult<Rational> r = ldqbd(q0, q1, q2);

    std::vector<Rational> ref(K + 1, Rational(1));
    for (std::size_t n = 1; n <= K; ++n)
        ref[n] = ref[n - 1] * lambda / (Rational(static_cast<long>(n < c ? n : c)) * mu);
    Rational tot(0);
    for (std::size_t n = 0; n <= K; ++n) tot += ref[n];
    for (std::size_t n = 0; n <= K; ++n) CHECK(r.pi.pi[n] == ref[n] / tot);  // exact

    const Matrix<Rational> Q = ldqbd_generator(q0, q1, q2);
    const std::vector<Rational> bal = line::vecmul(flatten(r.pi), Q);
    for (std::size_t j = 0; j < bal.size(); ++j) CHECK(bal[j] == Rational(0));
}

TEST_CASE("ldqbd handles matrix levels of heterogeneous order") {
    // Level 0 has two phases, the higher levels three, so R^(1) is 2 x 3 and
    // R^(n) is 3 x 3 above it. The oracle is the balance equation of the
    // assembled generator, checked exactly.
    const std::size_t N = 4;
    std::vector<Matrix<Rational>> q0, q1, q2;
    Matrix<Rational> up0(2, 3, Rational(0)), loc0(2, 2, Rational(0));
    up0(0, 0) = Rational(1);
    up0(1, 1) = Rational(2);
    loc0(0, 0) = Rational(-2);
    loc0(0, 1) = Rational(1);
    loc0(1, 0) = Rational(1, 2);
    loc0(1, 1) = Rational(-5, 2);
    q0.push_back(up0);
    q1.push_back(loc0);
    q2.push_back(Matrix<Rational>(2, 2, Rational(0)));  // unused level-0 down-block

    Matrix<Rational> upn(3, 3, Rational(0)), locn(3, 3, Rational(0)), dn(3, 3, Rational(0));
    for (std::size_t i = 0; i < 3; ++i) {
        upn(i, i) = Rational(1);
        dn(i, i) = Rational(3);
    }
    dn(0, 1) = Rational(1, 2);
    locn(0, 0) = Rational(-9, 2);
    locn(0, 1) = Rational(0);
    locn(1, 1) = Rational(-4);
    locn(1, 2) = Rational(0);
    locn(2, 2) = Rational(-4);
    Matrix<Rational> down1(3, 2, Rational(0));
    down1(0, 0) = Rational(3);
    down1(1, 1) = Rational(3);
    down1(2, 0) = Rational(3);

    for (std::size_t n = 1; n <= N; ++n) {
        q1.push_back(locn);
        q2.push_back(n == 1 ? down1 : dn);
        if (n < N) q0.push_back(upn);
    }
    // The top level cannot go up, so its local block must absorb the up-rate.
    Matrix<Rational> locTop = locn;
    for (std::size_t i = 0; i < 3; ++i) locTop(i, i) += Rational(1);
    q1[N] = locTop;
    // Close every row of the generator so it is conservative.
    Matrix<Rational> Q = ldqbd_generator(q0, q1, q2);
    for (std::size_t i = 0; i < Q.rows(); ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < Q.cols(); ++j) s += Q(i, j);
        Q(i, i) -= s;
    }
    // Push the diagonal repair back into the level blocks.
    std::vector<std::size_t> off(N + 2, 0);
    for (std::size_t n = 0; n <= N; ++n) off[n + 1] = off[n] + q1[n].rows();
    for (std::size_t n = 0; n <= N; ++n)
        for (std::size_t i = 0; i < q1[n].rows(); ++i)
            for (std::size_t j = 0; j < q1[n].cols(); ++j) q1[n](i, j) = Q(off[n] + i, off[n] + j);

    const LdqbdResult<Rational> r = ldqbd(q0, q1, q2);
    REQUIRE(r.R[1].rows() == 2u);
    REQUIRE(r.R[1].cols() == 3u);
    REQUIRE(r.R[2].rows() == 3u);

    Rational mass(0);
    for (std::size_t n = 0; n <= N; ++n) mass += r.pi.pi[n];
    CHECK(mass == Rational(1));  // exact normalization
    const std::vector<Rational> bal = line::vecmul(flatten(r.pi), Q);
    for (std::size_t j = 0; j < bal.size(); ++j) CHECK(bal[j] == Rational(0));  // exact balance
}

TEST_CASE("ldqbd_R refuses an exactly singular level rather than pseudo-inverting") {
    std::vector<Matrix<Rational>> q0, q1, q2;
    q0.push_back(m1(Rational(1)));
    q0.push_back(m1(Rational(1)));
    q1.push_back(m1(Rational(-1)));
    q1.push_back(m1(Rational(-2)));
    q1.push_back(m1(Rational(0)));  // singular top level
    q2.push_back(m1(Rational(0)));
    q2.push_back(m1(Rational(1)));
    q2.push_back(m1(Rational(1)));
    CHECK_THROWS_AS(ldqbd_R(q0, q1, q2), line::NumericError);
    // In double the reference's pinv fallback applies and returns zero there.
    std::vector<Matrix<double>> d0, d1, d2;
    for (std::size_t k = 0; k < q0.size(); ++k)
        d0.push_back(m1(line::num_traits<Rational>::to_double(q0[k](0, 0))));
    for (std::size_t k = 0; k < q1.size(); ++k)
        d1.push_back(m1(line::num_traits<Rational>::to_double(q1[k](0, 0))));
    for (std::size_t k = 0; k < q2.size(); ++k)
        d2.push_back(m1(line::num_traits<Rational>::to_double(q2[k](0, 0))));
    const std::vector<Matrix<double>> Rd = ldqbd_R(d0, d1, d2);
    CHECK(Rd[2](0, 0) == doctest::Approx(0.0));
}

// ---------------------------------------------------------------------------
// mmap_compress
// ---------------------------------------------------------------------------

namespace {

/** A two-class marked MAP with genuine autocorrelation. */
template <class T>
Mmap<T> two_class_mmap() {
    Mmap<T> m;
    m.D0 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D0(0, 0) = line::num_traits<T>::from_rational(-5, 1);
    m.D0(0, 1) = line::num_traits<T>::from_rational(1, 2);
    m.D0(1, 0) = line::num_traits<T>::from_rational(1, 5);
    m.D0(1, 1) = line::num_traits<T>::from_rational(-3, 2);
    Matrix<T> Dc1(2, 2), Dc2(2, 2);
    Dc1(0, 0) = line::num_traits<T>::from_rational(3, 1);
    Dc1(0, 1) = line::num_traits<T>::from_rational(1, 2);
    Dc1(1, 0) = line::num_traits<T>::from_rational(1, 10);
    Dc1(1, 1) = line::num_traits<T>::from_rational(1, 5);
    Dc2(0, 0) = line::num_traits<T>::from_rational(1, 1);
    Dc2(0, 1) = line::num_traits<T>::from_rational(0, 1);
    Dc2(1, 0) = line::num_traits<T>::from_rational(1, 2);
    Dc2(1, 1) = line::num_traits<T>::from_rational(1, 2);
    m.Dc.push_back(Dc1);
    m.Dc.push_back(Dc2);
    m.D1 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) m.D1(i, j) += m.Dc[c](i, j);
    return m;
}

}  // namespace

TEST_CASE("mmap_backward_moment reproduces the aggregate moments exactly") {
    // The M3A mixture law: M_k = sum_c B(c,k) p_c. Exact at Rational, which is
    // the point: a discrepancy in the 16th digit would be indistinguishable
    // from a wrong normalization otherwise.
    const Mmap<Rational> m = two_class_mmap<Rational>();
    REQUIRE(mmap_isfeasible(m));
    const std::vector<Rational> p = mmap_pc(m);
    std::vector<unsigned> orders;
    orders.push_back(1u);
    orders.push_back(2u);
    orders.push_back(3u);
    const std::vector<std::vector<Rational>> B = mmap_backward_moment(m, orders, true);
    for (std::size_t h = 0; h < 3; ++h) {
        Rational acc(0);
        for (std::size_t c = 0; c < 2; ++c) acc += B[c][h] * p[c];
        CHECK(acc == map_moment(m.map(), orders[h]));  // exact
    }
    // MATLAB mmap_backward_moment(MMAP, [1 2 3], 1), all 15 digits.
    const double refB[2][3] = {{0.284219703574542, 0.229034614387709, 0.373088751231497},
                               {0.51051051051051, 0.62657251846441, 1.23717581716121}};
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t h = 0; h < 3; ++h)
            CHECK(line::num_traits<Rational>::to_double(B[c][h]) ==
                  doctest::Approx(refB[c][h]).epsilon(1e-12));
    // The unnormalized form sums to the same moment without the weights.
    const std::vector<std::vector<Rational>> Bu = mmap_backward_moment(m, orders, false);
    for (std::size_t h = 0; h < 3; ++h) {
        Rational acc(0);
        for (std::size_t c = 0; c < 2; ++c) acc += Bu[c][h];
        CHECK(acc == map_moment(m.map(), orders[h]));  // exact
    }
}

TEST_CASE("mmap_mixture is a feasible MMAP with the prescribed class weights") {
    std::vector<Rational> alpha;
    alpha.push_back(Rational(1, 3));
    alpha.push_back(Rational(2, 3));
    std::vector<Map<Rational>> comps;
    comps.push_back(map_exponential(Rational(2)));
    comps.push_back(map_exponential(Rational(5)));
    const Mmap<Rational> mix = mmap_mixture(alpha, comps);
    CHECK(mmap_isfeasible(mix));  // exact partition of D1
    const std::vector<Rational> p = mmap_pc(mix);
    CHECK(p[0] == alpha[0]);  // exact
    CHECK(p[1] == alpha[1]);  // exact
    // The mixture of exponentials has the mixture mean, exactly.
    CHECK(map_mean(mix.map()) == alpha[0] / Rational(2) + alpha[1] / Rational(5));
}

TEST_CASE("mmap_compress preserves the class probabilities and the mean") {
    const Mmap<double> m = two_class_mmap<double>();
    const std::vector<double> p = mmap_pc(m);
    const double M1 = map_mean(m.map());
    const Mmap<double> c = mmap_compress(m);
    CHECK(mmap_isfeasible_tol(c, 1e-10));
    const std::vector<double> pc = mmap_pc(c);
    REQUIRE(pc.size() == 2u);
    CHECK(pc[0] == doctest::Approx(p[0]).epsilon(1e-9));
    CHECK(pc[1] == doctest::Approx(p[1]).epsilon(1e-9));
    // M1 is always exact: aph2_adjust never alters the first moment.
    CHECK(map_mean(c.map()) == doctest::Approx(M1).epsilon(1e-9));
    // MATLAB mmap_compress on this MMAP: order 4, mean 0.36734693877551,
    // scv 1.77944611277946, class probabilities 0.63265306122449 and
    // 0.36734693877551, per-class rates 1.72222222222222 and 1. The compressed
    // descriptor agrees entrywise, e.g. D0(0,0) = -5.02834192437788 and
    // D1(0,2) = 1.61514154048716.
    CHECK(c.order() == 4u);
    CHECK(map_mean(c.map()) == doctest::Approx(0.36734693877551).epsilon(1e-13));
    CHECK(map_scv(c.map()) == doctest::Approx(1.77944611277946).epsilon(1e-12));
    CHECK(c.D0(0, 0) == doctest::Approx(-5.02834192437788).epsilon(1e-12));
    CHECK(c.D0(0, 1) == doctest::Approx(0.631567730829488).epsilon(1e-12));
    CHECK(c.D0(2, 3) == doctest::Approx(2.30611970215571).epsilon(1e-12));
    CHECK(c.D1(0, 0) == doctest::Approx(2.78163265306123).epsilon(1e-12));
    CHECK(c.D1(0, 2) == doctest::Approx(1.61514154048716).epsilon(1e-12));
    // The per-class rates follow, lambda_c = p_c / M1.
    const std::vector<double> lam = mmap_count_lambda(c);
    const std::vector<double> lam0 = mmap_count_lambda(m);
    CHECK(lam[0] == doctest::Approx(lam0[0]).epsilon(1e-8));
    CHECK(lam[1] == doctest::Approx(lam0[1]).epsilon(1e-8));
    // And the autocorrelation is gone by construction: the result is renewal.
    // MATLAB map_acf(mmap_compress(MMAP), 1) = -1.24782989116879e-16 on this
    // instance, against 0.0739911109309197 for the uncompressed MMAP.
    std::vector<unsigned> lag1(1, 1u);
    CHECK(std::fabs(map_acf(c.map(), lag1)[0]) < 1e-9);
    CHECK(map_acf(m.map(), lag1)[0] == doctest::Approx(0.0739911109309197).epsilon(1e-12));
}

TEST_CASE("mmap_compress refuses the methods that are not ported") {
    const Mmap<double> m = two_class_mmap<double>();
    CHECK_THROWS_AS(mmap_compress(m, MmapCompressMethod::Mamap2), line::UnsupportedError);
    CHECK_THROWS_AS(mmap_compress(m, MmapCompressMethod::M3ppExactDelta), line::UnsupportedError);
}

// ---------------------------------------------------------------------------
// High-precision instantiation
// ---------------------------------------------------------------------------

TEST_CASE("the family instantiates at Real50") {
    // A double-only test cannot catch the Boost expression-template trap
    // (num_traits<expression<...>> is an incomplete type), and it cannot show
    // that the algorithms carry more digits when the arithmetic does. Both are
    // checked here on the M/M/1 collapse, where every answer is a rational
    // number whose decimal expansion is known: rho = 3/5, QN = 3/2.
    using R50 = line::Real50;
    const R50 lambda = line::num_traits<R50>::from_int(3) / line::num_traits<R50>::from_int(5);
    const R50 mu = line::num_traits<R50>::from_int(1);

    const QbdRapResult<R50> r = qbd_rap(m1(lambda), m1(R50(-(lambda + mu))), m1(mu), m1(lambda),
                                        m1(R50(-lambda)), std::size_t(6));
    // The arithmetic carries 50 digits, but qbd_rap does NOT reach them: its
    // Newton iteration stops at the reference's hard-wired residual target
    // 1e-10 * blockScale, so the last accepted step leaves about 1e-14 in R
    // and 6e-14 in QN whatever the working precision. Measured: 1 - G is
    // 2.773380e-14 at Real50 against 2.775558e-14 at double, i.e. 36 extra
    // digits of arithmetic buy nothing, which is the proof that the limit is
    // the stopping rule and not the number type. The tolerance is reproduced
    // rather than scaled to the arithmetic, because it is what the reference
    // and the JAR and Python ports agree on; ldqbd below is the contrast,
    // being a finite recursion with no tolerance at all.
    CHECK(line::num_traits<R50>::to_double(line::num_abs(R50(r.R(0, 0) - lambda))) < 1e-13);
    CHECK(line::num_traits<R50>::to_double(
              line::num_abs(R50(r.QN - line::num_traits<R50>::from_int(3) /
                                          line::num_traits<R50>::from_int(2)))) < 1e-12);

    // The setup-delay QBD, whose Coxian fit takes a square root.
    const R50 qn = qbd_setupdelayoff(lambda, mu, line::num_traits<R50>::from_int(10000000),
                                     line::num_traits<R50>::from_int(1),
                                     line::num_traits<R50>::from_int(10000000),
                                     line::num_traits<R50>::from_int(1));
    CHECK(line::num_traits<R50>::to_double(qn) == doctest::Approx(1.5).epsilon(1e-6));

    // The ETAQA descriptor and the departure joint moments.
    const Map<R50> Dp = qbd_depproc_etaqa(map_exponential(lambda), map_exponential(mu),
                                          std::size_t(4));
    CHECK(line::num_traits<R50>::to_double(qbd_depproc_residual(Dp)) ==
          doctest::Approx(0.6).epsilon(1e-12));
    std::vector<std::pair<unsigned, unsigned>> iset;
    iset.push_back(std::make_pair(1u, 0u));
    const std::vector<R50> jm =
        qbd_depproc_jointmom(map_exponential(lambda), map_exponential(mu), iset);
    CHECK(line::num_traits<R50>::to_double(jm[0]) == doctest::Approx(1.0 / 0.6).epsilon(1e-12));

    // The level-dependent recursion, which is exact in any arithmetic.
    std::vector<Matrix<R50>> q0, q1, q2;
    const std::size_t K = 4;
    for (std::size_t n = 0; n <= K; ++n) {
        R50 loc = line::num_traits<R50>::from_int(0);
        if (n == 0)
            loc = -lambda;
        else if (n == K)
            loc = -mu;
        else
            loc = -(lambda + mu);
        q1.push_back(m1(loc));
        q0.push_back(m1(lambda));
        q2.push_back(m1(mu));
    }
    const LdqbdResult<R50> ld = ldqbd(q0, q1, q2);
    // No tolerance anywhere in the recursion: this is accurate to the working
    // precision, 45 digits beyond what qbd_rap can offer above.
    CHECK(line::num_traits<R50>::to_double(line::num_abs(R50(ld.R[K](0, 0) - lambda))) < 1e-45);

    // The MAP/M/1-PS coefficients and the compression, which needs aph2_fit.
    const std::vector<std::vector<std::vector<R50>>> h =
        map_m1ps_h_recursive(m1(R50(-lambda)), m1(lambda), mu, std::size_t(3), std::size_t(3));
    CHECK(line::num_traits<R50>::to_double(h[0][0][0]) == doctest::Approx(1.0));
    const Mmap<R50> comp = mmap_compress(two_class_mmap<R50>());
    CHECK(comp.order() == 4u);
    CHECK(line::num_traits<R50>::to_double(map_mean(comp.map())) ==
          doctest::Approx(0.36734693877551).epsilon(1e-12));

    // And the predicate, which needs no transcendental function at all.
    Matrix<R50> Q(2, 2), Rd(2, 2, line::num_traits<R50>::from_int(0));
    Q(0, 0) = line::num_traits<R50>::from_int(-1);
    Q(0, 1) = line::num_traits<R50>::from_int(1);
    Q(1, 0) = line::num_traits<R50>::from_int(2);
    Q(1, 1) = line::num_traits<R50>::from_int(-2);
    Rd(0, 0) = line::num_traits<R50>::from_int(1);
    CHECK(mmdp_isfeasible(Q, Rd, R50(line::num_traits<R50>::from_int(0))));
}
