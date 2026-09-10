/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * QBD family: the rate matrix R, the fundamental matrix G, the boundary solve
 * and the MAP/MAP/1 queue, plus the MAP transformations they consume.
 *
 * Three kinds of oracle are used, in decreasing order of strength.
 *
 *  1. Closed forms that must hold exactly. An M/M/1 queue written as a QBD has
 *     R = rho and pi_k = (1-rho) rho^k; supplying that R as an exact rational
 *     makes the boundary solve, the level distribution and the queue-length
 *     moments exact rational identities, with no tolerance at all.
 *  2. The defining equation itself. The residual ||F + R L + R^2 B||_inf is
 *     evaluated rather than convergence being taken on trust, and it is
 *     evaluated at the tolerance the method actually guarantees: cyclic
 *     reduction is driven to 1e-14 on its own internal check and reaches a
 *     residual of order 1e-15, whereas successive substitution stops when two
 *     iterates differ by 1e-12 in the 1-norm and, contracting at the caudal
 *     rate eta = 0.7155, is left about tol/(1-eta) ~ 3.5e-12 from the true R.
 *     Those two are therefore asserted at 1e-13 and 1e-10 respectively.
 *  3. MATLAB reference values, obtained by running
 *       matlab -singleCompThread -batch "addpath(genpath('matlab/src'));
 *         addpath(genpath('matlab/lib')); ..."
 *     on the arrival MAP and service MAP defined below. MATLAB carries these in
 *     double arithmetic, so the comparison is made at 1e-12 relative for the
 *     quadratically convergent quantities and at 1e-9 for the linearly
 *     convergent qbd_R; asserting tighter would be asserting against MATLAB's
 *     own rounding.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/mam/qbd_r.h"

using line::Matrix;
using line::Rational;
using line::mam::Map;
using namespace line::mam;

namespace {

/** Two-state MMPP: phase rates q01, q10, arrival rates l0, l1. */
template <class T>
Map<T> mmpp2(const T& l0, const T& l1, const T& q01, const T& q10) {
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D1 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D0(0, 0) = -(l0 + q01);
    m.D0(0, 1) = q01;
    m.D0(1, 0) = q10;
    m.D0(1, 1) = -(l1 + q10);
    m.D1(0, 0) = l0;
    m.D1(1, 1) = l1;
    return m;
}

/**
 * The reference instance: arrival MMPP with rates 4 and 1 and switching rates
 * 1/2 and 1/3 (lambda = 2.2, SCV 1.7406), service Erlang-2 with mean 1/4
 * (rate 4), hence rho = 0.55.
 */
Map<double> ref_arrival() {
    return mmpp2<double>(4.0, 1.0, 0.5, 1.0 / 3.0);
}
Map<double> ref_service() { return map_erlang<double>(0.25, 2); }

/** M/M/1 as a QBD: F = lambda, L = -(lambda+mu), B = mu, Lbar = -lambda. */
template <class T>
void mm1_blocks(const T& lambda, const T& mu, Matrix<T>& B, Matrix<T>& L, Matrix<T>& F,
                Matrix<T>& Lbar) {
    B = Matrix<T>(1, 1, mu);
    L = Matrix<T>(1, 1, T(-(lambda + mu)));
    F = Matrix<T>(1, 1, lambda);
    Lbar = Matrix<T>(1, 1, T(-lambda));
}

}  // namespace

// ---------------------------------------------------------------------------
// M/M/1 as a QBD
// ---------------------------------------------------------------------------

TEST_CASE("M/M/1 as a QBD has R = rho and a geometric level distribution, exactly") {
    // lambda = 3, mu = 5, rho = 3/5. R is the 1x1 matrix rho, which is the
    // minimal non-negative root of lambda - R(lambda+mu) + R^2 mu = 0; the
    // other root is 1. Supplying it exactly makes every downstream quantity an
    // exact rational identity: no tolerance appears anywhere in this case.
    const Rational lambda(3), mu(5), rho(3, 5);
    Matrix<Rational> B, L, F, Lbar;
    mm1_blocks<Rational>(lambda, mu, B, L, F, Lbar);
    Matrix<Rational> R(1, 1, rho);

    // The defining equation holds identically, not approximately.
    CHECK(qbd_R_residual(B, L, F, R) == Rational(0));

    const Matrix<Rational> pq = qbd_pi(B, Lbar, R, static_cast<std::size_t>(12), Rational(0));
    REQUIRE(pq.rows() == 12);
    // pi_k = (1 - rho) rho^k
    Rational pk = Rational(1) - rho;
    for (std::size_t k = 0; k < pq.rows(); ++k) {
        CHECK(pq(k, 0) == pk);
        pk *= rho;
    }

    // Mean queue length rho/(1-rho) = 3/2, and E[N(N-1)] = 2 rho^2/(1-rho)^2.
    const std::vector<Rational> pi0(1, pq(0, 0));
    CHECK(qbd_qlen_factmoment(pi0, R, 1) == Rational(3, 2));
    CHECK(qbd_qlen_factmoment(pi0, R, 2) == Rational(2) * rho * rho /
                                                ((Rational(1) - rho) * (Rational(1) - rho)));
    // E[N^2] = rho(1+rho)/(1-rho)^2 = (3/5)(8/5)/(4/25) = 6
    CHECK(qbd_qlen_moment(pi0, R, 2) == Rational(6));
    CHECK(qbd_qlen_moment(pi0, R, 1) == Rational(3, 2));
}

TEST_CASE("the three R iterations recover rho on M/M/1") {
    const double lambda = 3.0, mu = 5.0, rho = 0.6;
    Matrix<double> B, L, F, Lbar;
    mm1_blocks<double>(lambda, mu, B, L, F, Lbar);

    // Cyclic reduction is driven to 1e-14 internally and is quadratically
    // convergent, so the fixed point is reached to machine precision.
    const QbdFundMat<double> fm = qbd_fundmat(B, L, F);
    CHECK(fm.R(0, 0) == doctest::Approx(rho).epsilon(1e-14));
    CHECK(qbd_R_residual(B, L, F, fm.R) < 1e-14);
    // For a birth-death chain G is the probability of ever going down, = 1.
    CHECK(fm.G(0, 0) == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(qbd_G_residual(B, L, F, fm.G) < 1e-14);

    // Logarithmic reduction, same guarantee.
    const Matrix<double> Rlr = qbd_R_logred(B, L, F);
    CHECK(Rlr(0, 0) == doctest::Approx(rho).epsilon(1e-13));
    CHECK(qbd_R_residual(B, L, F, Rlr) < 1e-13);

    // Successive substitution stops at 1e-12 on ||R_k - R_{k+1}||_1 and
    // contracts at rate rho = 0.6, so it is within tol/(1-rho) = 2.5e-12.
    const Matrix<double> Rss = qbd_R(B, L, F);
    CHECK(Rss(0, 0) == doctest::Approx(rho).epsilon(1e-11));
    CHECK(qbd_R_residual(B, L, F, Rss) < 1e-11);

    CHECK(qbd_caudal(fm.R) == doctest::Approx(rho).epsilon(1e-13));
}

TEST_CASE("MAP/MAP/1 with Poisson arrivals and exponential service is M/M/1") {
    // lambda = 0.6, mu = 1: QN = rho/(1-rho) = 1.5, UN = 0.6, RN = 2.5.
    const Map<double> a = map_exponential<double>(0.6);
    const Map<double> s = map_exponential<double>(1.0);
    const QbdMapMap1Result<double> r = qbd_mapmap1(a, s);
    CHECK(r.XN == doctest::Approx(0.6).epsilon(1e-14));
    CHECK(r.UN == doctest::Approx(0.6).epsilon(1e-12));
    // Closed form, so it is exact up to double rounding, unlike the truncated
    // level sum that MATLAB reports (1.499999997037549 for the same instance).
    CHECK(r.QN == doctest::Approx(1.5).epsilon(1e-12));
    CHECK(r.RN == doctest::Approx(2.5).epsilon(1e-12));
    CHECK(r.eta == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(qbd_R_residual(r.B, r.L, r.F, r.R) < 1e-14);

    // The truncated sum reproduces MATLAB's under-count, which is what tells
    // us the two are computing the same thing and only the tail differs.
    const double trunc = qbd_mapmap1_qlen_truncated(r);
    CHECK(trunc < 1.5);
    CHECK(1.5 - trunc < 1e-7);
}

// ---------------------------------------------------------------------------
// The reference instance, against MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("MAP descriptors of the reference arrival MAP match MATLAB") {
    // MATLAB: map_lambda 2.2, map_scv 1.740571428571429,
    //         map_skew 2.999445333567394, map_kurt 16.20096536840022,
    //         map_idc 3.356363636363637, map_moment(.,3) 1.231168831168831.
    // These are rational in the entries of (D0, D1) except map_skew, so the
    // exact and double paths must agree to double rounding; the comparison is
    // at 1e-13 relative, which is the accuracy MATLAB's own double arithmetic
    // has on a 2x2 solve.
    const Map<double> a = ref_arrival();
    CHECK(map_lambda(a) == doctest::Approx(2.2).epsilon(1e-13));
    CHECK(map_scv(a) == doctest::Approx(1.740571428571429).epsilon(1e-13));
    CHECK(map_skew(a) == doctest::Approx(2.999445333567394).epsilon(1e-13));
    CHECK(map_kurt(a) == doctest::Approx(16.20096536840022).epsilon(1e-13));
    CHECK(map_idc(a) == doctest::Approx(3.356363636363637).epsilon(1e-13));
    CHECK(map_moment(a, 3) == doctest::Approx(1.231168831168831).epsilon(1e-13));
    // map_joint(MAPa, [1 1], [1 1]) = 0.2590723562152134
    CHECK(map_joint(a, std::vector<unsigned>{1, 1}, std::vector<unsigned>{1, 1}) ==
          doctest::Approx(0.2590723562152134).epsilon(1e-13));
}

TEST_CASE("map_sum, map_renewal and map_scale match MATLAB") {
    const Map<double> a = ref_arrival();
    // map_sum(MAPa,3): mean 1.363636363636364, SCV 0.7317305344995141
    const Map<double> s3 = map_sum(a, 3);
    CHECK(map_mean(s3) == doctest::Approx(1.363636363636364).epsilon(1e-13));
    CHECK(map_scv(s3) == doctest::Approx(0.7317305344995141).epsilon(1e-13));
    // map_renewal preserves the marginal (SCV 1.740571428571429) and kills the
    // correlation: the lag-1 ACF drops to roundoff, the renewal value 0.
    const Map<double> rn = map_renewal(a);
    CHECK(map_scv(rn) == doctest::Approx(1.740571428571429).epsilon(1e-13));
    CHECK(std::fabs(map_acf(rn, std::vector<unsigned>{1})[0]) < 1e-12);
    // map_scale to mean 2 leaves the SCV invariant.
    const Map<double> sc = map_scale(a, 2.0);
    CHECK(map_mean(sc) == doctest::Approx(2.0).epsilon(1e-13));
    CHECK(map_scv(sc) == doctest::Approx(1.74057142857143).epsilon(1e-12));
}

TEST_CASE("cyclic and logarithmic reduction reproduce the MATLAB R and G") {
    const Map<double> a = ref_arrival();
    const Map<double> s = ref_service();
    const QbdMapMap1Blocks<double> blk = qbd_mapmap1_blocks(a, s);

    // MATLAB qbd_fundmat(B,L,F,'GR'), printed at 17 significant digits.
    const double Rref[4][4] = {
        {0.53816741149145442, 0.34846799677389373, 0.15500808729875568, 0.15153200322610633},
        {0.14550477417185245, 0.41537959948375219, 0.072762575634554599, 0.084620400516247807},
        {0.0091664872003227871, 0.0089609270859155026, 0.13481886045689551, 0.1160390729140845},
        {0.0028190341204719555, 0.0050040732179565438, 0.014682493336261749, 0.11999592678204346}};
    const double Gref[4][4] = {
        {0.69693599354778735, 0.0, 0.3030640064522126, 0.0},
        {0.83075919896750439, 0.0, 0.16924080103249559, 0.0},
        {0.071687416687324021, 0.0, 0.92831258331267585, 0.0},
        {0.040032585743652337, 0.0, 0.95996741425634768, 0.0}};

    const QbdFundMat<double> fm = qbd_fundmat(blk.B, blk.L, blk.F);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            // Absolute tolerance: several reference entries are exactly zero,
            // where a relative comparison is meaningless.
            CHECK(fm.R(i, j) == doctest::Approx(Rref[i][j]).epsilon(1e-12).scale(1.0));
            CHECK(fm.G(i, j) == doctest::Approx(Gref[i][j]).epsilon(1e-12).scale(1.0));
        }
    // MATLAB reports residuals 1.39e-15 for R and 1.78e-15 for G on this
    // instance; cyclic reduction guarantees nothing better than its 1e-14
    // internal check, so 1e-13 is the honest assertion.
    CHECK(qbd_R_residual(blk.B, blk.L, blk.F, fm.R) < 1e-13);
    CHECK(qbd_G_residual(blk.B, blk.L, blk.F, fm.G) < 1e-13);

    // Logarithmic reduction: MATLAB residual 7.77e-16, same guarantee.
    const Matrix<double> Rlr = qbd_R_logred(blk.B, blk.L, blk.F);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(Rlr(i, j) == doctest::Approx(Rref[i][j]).epsilon(1e-12).scale(1.0));
    CHECK(qbd_R_residual(blk.B, blk.L, blk.F, Rlr) < 1e-13);

    // Successive substitution: MATLAB residual 6.35e-12 on this instance, and
    // it stops on a 1e-12 iterate difference contracting at eta = 0.7155, so
    // the honest bound is tol/(1-eta) ~ 3.5e-12. Assert 1e-10, not tighter.
    const Matrix<double> Rss = qbd_R(blk.B, blk.L, blk.F);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(Rss(i, j) == doctest::Approx(Rref[i][j]).epsilon(1e-9).scale(1e-10));
    CHECK(qbd_R_residual(blk.B, blk.L, blk.F, Rss) < 1e-10);

    // Caudal characteristic: MATLAB QBD_Caudal 0.7155290134371834, max |eig(R)|
    // 0.7155290134371777. The bisection in QBD_Caudal stops at a 1e-15 bracket,
    // so the two MATLAB routes already differ at 6e-15; 1e-12 is the tolerance
    // that means anything here.
    CHECK(qbd_caudal(fm.R) == doctest::Approx(0.7155290134371777).epsilon(1e-12));
}

TEST_CASE("qbd_mapmap1 reproduces the MATLAB MAP/MAP/1 solution") {
    const Map<double> a = ref_arrival();
    const Map<double> s = ref_service();
    const QbdMapMap1Result<double> r = qbd_mapmap1(a, s);

    // MATLAB: XN 2.2, UN 0.55, eta 0.7155290134371834,
    //         pi0 = [0.086081961314618455 -1.94e-17 0.36391803868538142 6.93e-17]
    CHECK(r.XN == doctest::Approx(2.2).epsilon(1e-13));
    CHECK(r.UN == doctest::Approx(0.55).epsilon(1e-12));
    CHECK(r.eta == doctest::Approx(0.7155290134371834).epsilon(1e-12));
    CHECK(r.pi0[0] == doctest::Approx(0.086081961314618455).epsilon(1e-11));
    CHECK(r.pi0[2] == doctest::Approx(0.36391803868538142).epsilon(1e-11));
    // The odd phases carry no boundary mass: the service Erlang cannot be in
    // its second phase while the system is empty.
    CHECK(line::num_abs(r.pi0[1]) < 1e-14);
    CHECK(line::num_abs(r.pi0[3]) < 1e-14);

    // MATLAB QN 1.742822351190845 is the truncated level sum, stopped once the
    // accumulated mass reaches 1 - 1e-10 (68 levels). Our truncated sum must
    // reproduce it closely; the closed form is strictly larger, by the tail
    // mass times its level index, which at eta = 0.7155 and 1e-10 residual
    // mass is of order 1e-8.
    const double trunc = qbd_mapmap1_qlen_truncated(r);
    CHECK(trunc == doctest::Approx(1.742822351190845).epsilon(1e-9));
    CHECK(r.QN > trunc);
    CHECK(r.QN - 1.742822351190845 < 1e-6);
    CHECK(r.RN == doctest::Approx(r.QN / 2.2).epsilon(1e-14));

    // Little's law and the utilization identity, both structural.
    CHECK(r.UN == doctest::Approx(map_lambda(a) / map_lambda(s)).epsilon(1e-12));
    CHECK(qbd_R_residual(r.B, r.L, r.F, r.R) < 1e-13);
    // U = L + R B is a generator of the taboo process: its row sums vanish
    // only after adding F, but L + R B + F(...) is not needed here; what must
    // hold is that U has non-positive diagonal and non-negative off-diagonal.
    for (std::size_t i = 0; i < r.U.rows(); ++i) {
        CHECK(r.U(i, i) <= 1e-14);
        for (std::size_t j = 0; j < r.U.cols(); ++j)
            if (i != j) CHECK(r.U(i, j) >= -1e-14);
    }
}

TEST_CASE("qbd_rg agrees with qbd_mapmap1 and honours the utilization argument") {
    const Map<double> a = ref_arrival();
    const Map<double> s = ref_service();
    const QbdRg<double> rg = qbd_rg(a, s);
    const QbdMapMap1Result<double> r = qbd_mapmap1(a, s);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(rg.R(i, j) == doctest::Approx(r.R(i, j)).epsilon(1e-14).scale(1.0));

    // With util = 0.8 the service MAP is rescaled to mean util/lambda_a, so the
    // realized utilization is exactly 0.8.
    const QbdMapMap1Result<double> r8 =
        qbd_mapmap1(a, s, 0.8, static_cast<std::size_t>(20000));
    CHECK(map_lambda(a) / map_lambda(r8.service) == doctest::Approx(0.8).epsilon(1e-12));
    CHECK(r8.UN == doctest::Approx(0.8).epsilon(1e-8));
    CHECK(r8.QN > r.QN);  // a busier queue is a longer queue
}

// ---------------------------------------------------------------------------
// MAP transformations, exact
// ---------------------------------------------------------------------------

TEST_CASE("map_erlang and map_kurt are exact in rational arithmetic") {
    // Erlang-2 with mean 1: SCV 1/2, third moment 3/2 * ... , kurtosis 3 + 6/k.
    const Map<Rational> e2 = map_erlang<Rational>(Rational(1), 2);
    CHECK(map_mean(e2) == Rational(1));
    CHECK(map_scv(e2) == Rational(1, 2));
    CHECK(map_kurt(e2) == Rational(6));  // 3 + 6/k with k = 2
    CHECK(map_isfeasible(e2));

    // A k-fold self-convolution of an exponential is an Erlang-k.
    const Map<Rational> ex = map_exponential_mean<Rational>(Rational(1, 3));
    const Map<Rational> s3 = map_sum(ex, 3);
    CHECK(map_mean(s3) == Rational(1));
    CHECK(map_scv(s3) == Rational(1, 3));
    CHECK(map_isfeasible(s3));
}

TEST_CASE("map_normalize, map_scale and map_renewal preserve the MAP identities exactly") {
    const Map<Rational> a =
        mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    CHECK(map_lambda(a) == Rational(11, 5));  // 2.2
    CHECK(map_isfeasible(a));

    const Map<Rational> sc = map_scale(a, Rational(2));
    CHECK(map_mean(sc) == Rational(2));
    CHECK(map_scv(sc) == map_scv(a));  // scaling time leaves the SCV invariant
    CHECK(map_isfeasible(sc));

    const Map<Rational> rn = map_renewal(a);
    CHECK(map_mean(rn) == map_mean(a));
    CHECK(map_scv(rn) == map_scv(a));
    // A renewal MAP has no autocorrelation: with map_acf.m's normalization the
    // coefficient is exactly zero at every lag, in exact arithmetic.
    const std::vector<Rational> acf = map_acf(rn, std::vector<unsigned>{1, 2, 5});
    for (const Rational& v : acf) CHECK(v == Rational(0));
}

TEST_CASE("map2ph / ph2map round-trips a renewal MAP exactly") {
    const Map<Rational> e2 = map_erlang<Rational>(Rational(1), 2);
    const PhType<Rational> ph = map2ph(e2);
    const Map<Rational> back = ph2map(ph);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK(back.D0(i, j) == e2.D0(i, j));
            CHECK(back.D1(i, j) == e2.D1(i, j));
        }
}

TEST_CASE("map_stochcomp on the full phase set is the identity, exactly") {
    const Map<Rational> a =
        mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    const Map<Rational> full = map_stochcomp(a, std::vector<std::size_t>{0, 1});
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK(full.D0(i, j) == a.D0(i, j));
            CHECK(full.D1(i, j) == a.D1(i, j));
        }
    // Censoring an MMPP down to one phase leaves a Poisson process whose rate
    // is the censored one; the result is still a feasible MAP.
    const Map<Rational> one = map_stochcomp(a, std::vector<std::size_t>{0});
    CHECK(one.order() == 1);
    CHECK(map_isfeasible(one));
}

TEST_CASE("map_sumind and map_mixture build feasible MAPs with the right mean") {
    std::vector<Map<Rational>> ms;
    ms.push_back(map_exponential_mean<Rational>(Rational(1, 2)));
    ms.push_back(map_erlang<Rational>(Rational(1), 2));
    // The sum of independent inter-arrival times has the sum of the means.
    const Map<Rational> si = map_sumind(ms);
    CHECK(map_mean(si) == Rational(3, 2));
    CHECK(map_isfeasible(si));

    // A mixture with weights (1/4, 3/4) has the mixed mean.
    std::vector<Rational> alpha{Rational(1, 4), Rational(3, 4)};
    const Map<Rational> mx = map_mixture(alpha, ms);
    CHECK(map_mean(mx) == Rational(1, 4) * Rational(1, 2) + Rational(3, 4) * Rational(1));
    CHECK(map_isfeasible(mx));
}

TEST_CASE("map_hyperexp matches its target moments") {
    // Two-phase hyperexponential, mean 1 and SCV 4. The fit is exact in the
    // first two moments by construction, so the only error is the square root;
    // 1e-12 relative is what double arithmetic supports here.
    const Map<double> h = map_hyperexp<double>(1.0, 4.0);
    CHECK(map_mean(h) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(map_scv(h) == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(map_isfeasible(h, 1e-12));
}

TEST_CASE("the boundary solve is exact given R, at high precision too") {
    // Same M/M/1 QBD carried at 100 decimal digits: pi_k = (1-rho) rho^k must
    // hold to the working precision, which is what separates an error in the
    // boundary solve from an error in R.
    using R100 = line::Real100;
    const R100 lambda(3), mu(5);
    const R100 rho = R100(3) / R100(5);
    Matrix<R100> B, L, F, Lbar;
    mm1_blocks<R100>(lambda, mu, B, L, F, Lbar);
    const QbdFundMat<R100> fm = qbd_fundmat(B, L, F, 200u, R100("1e-80"));
    CHECK(line::num_abs(R100(fm.R(0, 0) - rho)) < R100("1e-70"));
    CHECK(qbd_R_residual(B, L, F, fm.R) < R100("1e-70"));

    const Matrix<R100> pq = qbd_pi(B, Lbar, fm.R, static_cast<std::size_t>(10), R100(0));
    R100 pk = R100(1) - rho;
    for (std::size_t k = 0; k < pq.rows(); ++k) {
        CHECK(line::num_abs(R100(pq(k, 0) - pk)) < R100("1e-68"));
        pk *= rho;
    }
    const std::vector<R100> pi0(1, pq(0, 0));
    CHECK(line::num_abs(R100(qbd_qlen_factmoment(pi0, fm.R, 1) - R100(3) / R100(2))) <
          R100("1e-66"));
}
