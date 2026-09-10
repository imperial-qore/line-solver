/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The pfqn sensitivity family. Every reference value below was produced by
 * running the MATLAB routine of the same name on the same model, at the
 * tolerance that routine guarantees: the analytic kernels are asserted at
 * 1e-12 relative because they are exact up to floating-point rounding, the
 * Linearizer at the accuracy band its own reference reports.
 *
 * Beyond the reference values the file asserts the two properties that give
 * the derivatives their meaning: agreement with a central finite difference of
 * the solver being differentiated, and, where the derivative is a rational
 * function of the inputs, exactness at line::Rational.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_sens.h"
#include "line/api/pfqn/pfqn_sens_ldmx_ec.h"
#include "line/api/pfqn/pfqn_sens_linearizer.h"
#include "line/api/pfqn/pfqn_sens_mom.h"
#include "line/api/pfqn/pfqn_sens_mva.h"
#include "line/api/pfqn/pfqn_sens_mvaldmx.h"
#include "line/api/pfqn/pfqn_sens_respt.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

namespace {

Matrix<double> demands2x2() {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5; L(0, 1) = 0.3;
    L(1, 0) = 0.4; L(1, 1) = 0.6;
    return L;
}

const std::vector<int> N21{2, 1};
const std::vector<double> Z12{1.0, 2.0};

/** Mean queue lengths and throughputs from exact MVA, for finite differences. */
MvaResult<double> mva(const Matrix<double>& L, const std::vector<int>& N,
                      const std::vector<double>& Z) {
    Matrix<double> Zm(1, Z.size());
    for (std::size_t r = 0; r < Z.size(); ++r) Zm(0, r) = Z[r];
    return pfqn_mva(L, N, Zm, std::vector<int>());
}

}  // namespace

TEST_CASE("pfqn_sens matches the MATLAB reference on a two-station two-class model") {
    // MATLAB: sens = pfqn_sens([0.5 0.3;0.4 0.6],[2 1],[1 2]).
    const auto s = pfqn_sens(demands2x2(), N21, Z12);
    REQUIRE(s.params.size() == 6);
    CHECK(s.XN[0] == doctest::Approx(0.87524106215695).epsilon(1e-12));
    CHECK(s.XN[1] == doctest::Approx(0.298175344904317).epsilon(1e-12));
    CHECK(s.QN(0, 0) == doctest::Approx(0.609701824655096).epsilon(1e-12));
    CHECK(s.QN(1, 0) == doctest::Approx(0.515057113187954).epsilon(1e-12));
    CHECK(s.QN(0, 1) == doctest::Approx(0.142857142857143).epsilon(1e-12));
    CHECK(s.QN(1, 1) == doctest::Approx(0.260792167334223).epsilon(1e-12));
    CHECK(s.CN(0, 0) == doctest::Approx(0.696610169491525).epsilon(1e-12));
    CHECK(s.CN(1, 1) == doctest::Approx(0.874626865671642).epsilon(1e-12));

    // parameter order: L(0,0) L(0,1) L(1,0) L(1,1) Z(0) Z(1)
    CHECK(s.params[0].type == 'L');
    CHECK(s.params[0].station == 0);
    CHECK(s.params[0].cls == 0);
    CHECK(s.params[4].type == 'Z');
    CHECK(s.params[4].station == -1);

    CHECK(s.dX(0, 0) == doctest::Approx(-0.592565128457221).epsilon(1e-12));
    CHECK(s.dX(1, 0) == doctest::Approx(-0.00756584113833031).epsilon(1e-10));
    CHECK(s.dX(0, 2) == doctest::Approx(-0.607787037460662).epsilon(1e-12));
    CHECK(s.dX(1, 3) == doctest::Approx(-0.129602990738711).epsilon(1e-12));
    CHECK(s.dX(0, 4) == doctest::Approx(-0.335843682944074).epsilon(1e-12));
    CHECK(s.dX(1, 5) == doctest::Approx(-0.0889085363088083).epsilon(1e-12));

    CHECK(s.dQ[0](0, 0) == doctest::Approx(0.995141818919317).epsilon(1e-12));
    CHECK(s.dQ[0](1, 0) == doctest::Approx(-0.402576690462095).epsilon(1e-12));
    CHECK(s.dQ[0](0, 1) == doctest::Approx(0.0839214190349037).epsilon(1e-12));
    CHECK(s.dQ[2](1, 0) == doctest::Approx(1.11100790053828).epsilon(1e-12));
    CHECK(s.dQ[4](0, 0) == doctest::Approx(-0.296282564228611).epsilon(1e-12));
    CHECK(s.dQ[5](1, 1) == doctest::Approx(-0.0777617944432264).epsilon(1e-12));

    CHECK(s.QVar(0, 0) == doctest::Approx(0.497570909459658).epsilon(1e-12));
    CHECK(s.QVar(1, 1) == doctest::Approx(0.192779612791342).epsilon(1e-12));
    CHECK(s.QTotVar[0] == doctest::Approx(0.703941308086399).epsilon(1e-12));
    CHECK(s.QTotVar[1] == doctest::Approx(0.753049734886369).epsilon(1e-12));
    // MATLAB s.QCov(1,1,2,2), i.e. Cov[n(0,0),n(1,1)] in 0-based indices
    CHECK(s.QCov(0 * 2 + 0, 1 * 2 + 1) == doctest::Approx(-0.0343948683791215).epsilon(1e-10));
    CHECK(std::fabs(s.QCovAsym) < 1e-14);
}

TEST_CASE("pfqn_sens agrees with a central finite difference of exact MVA") {
    // The derivative is analytic, the comparison is not: a central difference
    // of step eps carries an O(eps^2) truncation error plus an O(1e-16/eps)
    // cancellation error, minimized near eps ~ 1e-5 for these magnitudes. The
    // band below is what that trade-off allows, not what the kernel achieves.
    const Matrix<double> L = demands2x2();
    const auto s = pfqn_sens(L, N21, Z12);
    const double eps = 1e-5;
    for (std::size_t p = 0; p < s.params.size(); ++p) {
        Matrix<double> Lp = L, Lm = L;
        std::vector<double> Zp = Z12, Zm = Z12;
        if (s.params[p].type == 'L') {
            const std::size_t i = static_cast<std::size_t>(s.params[p].station);
            Lp(i, s.params[p].cls) += eps;
            Lm(i, s.params[p].cls) -= eps;
        } else {
            Zp[s.params[p].cls] += eps;
            Zm[s.params[p].cls] -= eps;
        }
        const auto hi = mva(Lp, N21, Zp);
        const auto lo = mva(Lm, N21, Zm);
        for (std::size_t r = 0; r < 2; ++r) {
            const double fd = (hi.XN[r] - lo.XN[r]) / (2 * eps);
            CHECK(s.dX(r, p) == doctest::Approx(fd).epsilon(1e-7).scale(1e-3));
            for (std::size_t i = 0; i < 2; ++i) {
                const double fdq = (hi.QN(i, r) - lo.QN(i, r)) / (2 * eps);
                CHECK(s.dQ[p](i, r) == doctest::Approx(fdq).epsilon(1e-7).scale(1e-3));
            }
        }
    }
}

TEST_CASE("pfqn_sens repairman model matches MATLAB and the general kernel") {
    // MATLAB: pfqn_sens([0.6 0.4],[2 1],[1 2]) dispatches to the CoMoM kernel.
    Matrix<double> L(1, 2);
    L(0, 0) = 0.6; L(0, 1) = 0.4;
    const auto s = pfqn_sens(L, N21, Z12);
    REQUIRE(s.params.size() == 4);
    CHECK(s.XN[0] == doctest::Approx(1.01190476190476).epsilon(1e-12));
    CHECK(s.XN[1] == doctest::Approx(0.362103174603175).epsilon(1e-12));
    CHECK(s.QN(0, 0) == doctest::Approx(0.988095238095239).epsilon(1e-12));
    CHECK(s.QN(0, 1) == doctest::Approx(0.275793650793651).epsilon(1e-12));
    CHECK(s.dX(0, 0) == doctest::Approx(-0.971986016628878).epsilon(1e-11));
    CHECK(s.dX(1, 1) == doctest::Approx(-0.249664391219451).epsilon(1e-11));
    CHECK(s.dX(0, 2) == doctest::Approx(-0.428713151927441).epsilon(1e-11));
    CHECK(s.dX(1, 3) == doctest::Approx(-0.131118709057697).epsilon(1e-11));
    CHECK(s.dQ[0](0, 0) == doctest::Approx(0.971986016628873).epsilon(1e-11));
    CHECK(s.dQ[1](0, 1) == doctest::Approx(0.499328782438901).epsilon(1e-11));
    CHECK(s.dQ[2](0, 0) == doctest::Approx(-0.583191609977327).epsilon(1e-11));
    CHECK(s.dQ[3](0, 1) == doctest::Approx(-0.0998657564877804).epsilon(1e-11));
    CHECK(s.QVar(0, 0) == doctest::Approx(0.583191609977324).epsilon(1e-11));

    // The two kernels are two exact evaluations of the same derivative, so at
    // double they may differ only by rounding.
    const auto d = pfqn_sens_dmva(L, N21, Z12, std::vector<int>());
    for (std::size_t p = 0; p < 4; ++p) {
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(s.dX(r, p) == doctest::Approx(d.dX(r, p)).epsilon(1e-10).scale(1e-6));
            CHECK(s.dQ[p](0, r) == doctest::Approx(d.dQ[p](0, r)).epsilon(1e-10).scale(1e-6));
        }
    }
}

TEST_CASE("the two pfqn_sens kernels agree exactly in rational arithmetic") {
    // Both are exact differentiations of the same finite recursion, so over the
    // rationals they are not merely close: they are the same number.
    Matrix<Rational> L(1, 2);
    L(0, 0) = Rational(3, 5);
    L(0, 1) = Rational(2, 5);
    const std::vector<Rational> Z{Rational(1), Rational(2)};
    const auto c = pfqn_sens_comom(L, N21, Z);
    const auto d = pfqn_sens_dmva(L, N21, Z, std::vector<int>());
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(c.XN[r] == d.XN[r]);
        CHECK(c.QN(0, r) == d.QN(0, r));
        for (std::size_t p = 0; p < 4; ++p) {
            CHECK(c.dX(r, p) == d.dX(r, p));
            CHECK(c.dQ[p](0, r) == d.dQ[p](0, r));
            CHECK(c.dU[p](0, r) == d.dU[p](0, r));
        }
    }
}

TEST_CASE("pfqn_sens reproduces the exactly differentiated single-job model") {
    // One station, one class, N = 1: X = 1/(Z+L), Q = L/(Z+L), so
    // dQ/dL = Z/(Z+L)^2, dX/dL = dX/dZ = -1/(Z+L)^2, dQ/dZ = -L/(Z+L)^2.
    // With L and Z rational the derivatives are rational and must match exactly.
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(3, 7);
    const std::vector<Rational> Z{Rational(5, 4)};
    const std::vector<int> N{1};
    const auto s = pfqn_sens_dmva(L, N, Z, std::vector<int>());
    const Rational den = Z[0] + L(0, 0);
    const Rational den2 = den * den;
    CHECK(s.XN[0] == Rational(1) / den);
    CHECK(s.QN(0, 0) == L(0, 0) / den);
    CHECK(s.dQ[0](0, 0) == Z[0] / den2);
    CHECK(s.dX(0, 0) == -Rational(1) / den2);
    CHECK(s.dQ[1](0, 0) == -L(0, 0) / den2);
    CHECK(s.dX(0, 1) == -Rational(1) / den2);
    // dU/dL = X + L dX/dL, an independent route to the same rational
    CHECK(s.dU[0](0, 0) == s.XN[0] + L(0, 0) * s.dX(0, 0));
}

TEST_CASE("pfqn_sens_mva matches MATLAB and closes the covariance identity exactly") {
    // MATLAB: mom = pfqn_sens_mva([0.5 0.3;0.4 0.6],[2 1],[1 2]).
    const auto m = pfqn_sens_mva(demands2x2(), N21, Z12);
    CHECK(m.QCov[0](0, 0) == doctest::Approx(0.497570909459658).epsilon(1e-12));
    CHECK(m.QCov[1](0, 0) == doctest::Approx(0.444403160215313).epsilon(1e-12));
    CHECK(m.QCov[0](1, 0) == doctest::Approx(0.0419607095174518).epsilon(1e-11));
    CHECK(m.QCov[1](1, 0) == doctest::Approx(0.0579334809398571).epsilon(1e-11));
    CHECK(m.QCov[0](1, 1) == doctest::Approx(0.122448979591837).epsilon(1e-12));
    CHECK(m.QCov[1](1, 1) == doctest::Approx(0.192779612791342).epsilon(1e-12));
    CHECK(m.QTotVar[0] == doctest::Approx(0.703941308086399).epsilon(1e-12));
    CHECK(m.QTotVar[1] == doctest::Approx(0.753049734886369).epsilon(1e-12));
    // the covariance matrix must be symmetric and its diagonal non-negative
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(m.QCov[i](0, 1) == doctest::Approx(m.QCov[i](1, 0)).epsilon(1e-14));
        CHECK(m.QVar(i, 0) >= 0.0);
        CHECK(m.QVar(i, 1) >= 0.0);
    }

    // Cov[n(i,r),n(i,s)] = L(i,s) dQ(i,r)/dL(i,s) is the identity both routes
    // rest on; over the rationals the two exact expressions must coincide.
    Matrix<Rational> Lq(2, 2);
    Lq(0, 0) = Rational(1, 2); Lq(0, 1) = Rational(3, 10);
    Lq(1, 0) = Rational(2, 5); Lq(1, 1) = Rational(3, 5);
    const std::vector<Rational> Zq{Rational(1), Rational(2)};
    const auto mq = pfqn_sens_mva(Lq, N21, Zq, std::vector<int>());
    const auto sq = pfqn_sens_dmva(Lq, N21, Zq, std::vector<int>());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r)
            for (std::size_t s = 0; s < 2; ++s)
                CHECK(mq.QCov[i](r, s) == Lq(i, s) * sq.dQ[i * 2 + s](i, r));
    CHECK(mq.QCovAsym == Rational(0));
}

TEST_CASE("pfqn_sens_mom matches the MATLAB reference for totals and per class") {
    // MATLAB: pfqn_sens_mom([0.5 0.3;0.4 0.6],[2 1],[1 2]) with the default
    // grouping, i.e. the per-station totals.
    const auto mm = pfqn_sens_mom(demands2x2(), N21, Z12);
    CHECK(mm.m(0, 0) == doctest::Approx(0.752558967512239).epsilon(1e-12));
    CHECK(mm.m(1, 0) == doctest::Approx(0.775849280522178).epsilon(1e-12));
    CHECK(mm.Var(0, 0) == doctest::Approx(0.703941308086399).epsilon(1e-12));
    CHECK(mm.Var(1, 0) == doctest::Approx(0.753049734886369).epsilon(1e-12));
    CHECK(mm.M2(0, 0) == doctest::Approx(1.27028630766949).epsilon(1e-12));
    CHECK(mm.M2(1, 0) == doctest::Approx(1.35499184097315).epsilon(1e-12));
    CHECK(mm.M3(0, 0) == doctest::Approx(2.50600801068091).epsilon(1e-12));
    CHECK(mm.M3(1, 0) == doctest::Approx(2.76961875092716).epsilon(1e-12));
    CHECK(mm.Cov(0, 1) == doctest::Approx(-0.296674588353213).epsilon(1e-11));
    CHECK(mm.Cov(1, 0) == doctest::Approx(-0.296674588353213).epsilon(1e-11));

    // The variance of the total must reproduce pfqn_sens_mva's QTotVar, which
    // is a completely different recursion.
    const auto mv = pfqn_sens_mva(demands2x2(), N21, Z12);
    CHECK(mm.Var(0, 0) == doctest::Approx(mv.QTotVar[0]).epsilon(1e-12));
    CHECK(mm.Var(1, 0) == doctest::Approx(mv.QTotVar[1]).epsilon(1e-12));

    // MATLAB: pfqn_sens_mom(L,N,Z,[],[1 2]), one group per class.
    const std::vector<int> perClass{0, 1};
    const auto mc = pfqn_sens_mom(demands2x2(), N21, Z12, std::vector<int>(), perClass);
    CHECK(mc.m(0, 0) == doctest::Approx(0.609701824655096).epsilon(1e-12));
    CHECK(mc.m(1, 0) == doctest::Approx(0.515057113187954).epsilon(1e-12));
    CHECK(mc.m(0, 1) == doctest::Approx(0.142857142857143).epsilon(1e-12));
    CHECK(mc.m(1, 1) == doctest::Approx(0.260792167334223).epsilon(1e-12));
    CHECK(mc.Var(0, 0) == doctest::Approx(0.497570909459658).epsilon(1e-12));
    CHECK(mc.Var(1, 1) == doctest::Approx(0.192779612791342).epsilon(1e-12));
    CHECK(mc.M3(0, 0) == doctest::Approx(1.38851802403204).epsilon(1e-12));
    CHECK(mc.M3(1, 1) == doctest::Approx(0.260792167334224).epsilon(1e-11));
    // the per-class second moments must reproduce pfqn_sens_mva's variances
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r)
            CHECK(mc.Var(i, r) == doctest::Approx(mv.QVar(i, r)).epsilon(1e-11));
}

TEST_CASE("pfqn_sens_mom agrees with a finite difference of the mean queue length") {
    // Var[Q_i] = x_i dm_i/dx_i, i.e. the derivative of the total queue length
    // at station i with respect to a rescaling of station i's demand column.
    const Matrix<double> L = demands2x2();
    const auto mm = pfqn_sens_mom(L, N21, Z12);
    const double eps = 1e-5;
    for (std::size_t h = 0; h < 2; ++h) {
        Matrix<double> Lp = L, Lm = L;
        for (std::size_t r = 0; r < 2; ++r) {
            Lp(h, r) *= (1.0 + eps);
            Lm(h, r) *= (1.0 - eps);
        }
        const auto hi = mva(Lp, N21, Z12);
        const auto lo = mva(Lm, N21, Z12);
        for (std::size_t i = 0; i < 2; ++i) {
            double sh = 0, sl = 0;
            for (std::size_t r = 0; r < 2; ++r) {
                sh += hi.QN(i, r);
                sl += lo.QN(i, r);
            }
            const double fd = (sh - sl) / (2 * eps);
            CHECK(mm.dm(i * 1 + 0, h * 1 + 0) == doctest::Approx(fd).epsilon(1e-7).scale(1e-3));
        }
    }
}

TEST_CASE("pfqn_sens_mom moments are exact rationals") {
    Matrix<Rational> L(2, 2);
    L(0, 0) = Rational(1, 2); L(0, 1) = Rational(3, 10);
    L(1, 0) = Rational(2, 5); L(1, 1) = Rational(3, 5);
    const std::vector<Rational> Z{Rational(1), Rational(2)};
    const auto mm = pfqn_sens_mom(L, N21, Z, std::vector<int>(), std::vector<int>());
    const auto mv = pfqn_sens_mva(L, N21, Z, std::vector<int>());
    // Var of the total from the second-order pass equals the sum of the exact
    // per-station covariances, as rationals.
    for (std::size_t i = 0; i < 2; ++i) CHECK(mm.Var(i, 0) == mv.QTotVar[i]);
    // E[Q^2] = Var + m^2 exactly
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(mm.M2(i, 0) == mm.Var(i, 0) + mm.m(i, 0) * mm.m(i, 0));
    CHECK(mm.CovAsym == Rational(0));
}

TEST_CASE("pfqn_sens_linearizer matches MATLAB and stays within its accuracy band") {
    // MATLAB: pfqn_sens_linearizer([0.5 0.3;0.4 0.6],[4 3],[1 2]) gives
    // m = [1.98180603934677 2.5807421101326], Var = [2.73327650402531
    // 3.31376366497691], M3 = [26.8124222923514 44.4981244002154] in 76 CORE
    // iterations. Asserted at 1e-6 relative: this is a fixed point stopped on
    // the reference's own tolerance 1/(4000+16 sum(n)), so the trailing digits
    // belong to the stopping test, not to the algorithm.
    const std::vector<int> N{4, 3};
    const auto sl = pfqn_sens_linearizer(demands2x2(), N, Z12);
    CHECK(sl.m[0] == doctest::Approx(1.98180603934677).epsilon(1e-6));
    CHECK(sl.m[1] == doctest::Approx(2.5807421101326).epsilon(1e-6));
    CHECK(sl.Var[0] == doctest::Approx(2.73327650402531).epsilon(1e-6));
    CHECK(sl.Var[1] == doctest::Approx(3.31376366497691).epsilon(1e-6));
    CHECK(sl.M2[0] == doctest::Approx(6.66083168161663).epsilon(1e-6));
    CHECK(sl.M3[0] == doctest::Approx(26.8124222923514).epsilon(1e-6));
    CHECK(sl.M3[1] == doctest::Approx(44.4981244002154).epsilon(1e-6));
    CHECK(sl.iter == 76u);

    // Against the exact pfqn_sens_mom on the same model. The reference reports
    // errors below 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3] over its
    // benchmark set; this model sits well inside those bands, and the bands are
    // what is asserted -- the Linearizer is an approximation and is expected to
    // disagree with the exact answer.
    const auto ex = pfqn_sens_mom(demands2x2(), N, Z12);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(std::fabs(sl.m[i] - ex.m(i, 0)) / ex.m(i, 0) < 0.021);
        CHECK(std::fabs(sl.M2[i] - ex.M2(i, 0)) / ex.M2(i, 0) < 0.041);
        CHECK(std::fabs(sl.M3[i] - ex.M3(i, 0)) / ex.M3(i, 0) < 0.062);
    }
    // Cov is symmetric by construction; CovAsym measures how far the fixed
    // point was from enforcing it, and must stay small on a benign model.
    CHECK(sl.CovAsym < 0.05);
    CHECK(sl.Cov(0, 1) == doctest::Approx(sl.Cov(1, 0)).epsilon(1e-14));
}

TEST_CASE("pfqn_sens_respt matches the MATLAB reference on a two-server center") {
    // MATLAB: pfqn_sens_respt([0.5;0.4],[1 1;2 1],[2 1],[1 2],[1;2],3).
    const std::vector<double> S{0.5, 0.4};
    Matrix<double> V(2, 2);
    V(0, 0) = 1.0; V(0, 1) = 1.0;
    V(1, 0) = 2.0; V(1, 1) = 1.0;
    const std::vector<int> b{1, 2};
    const auto r = pfqn_sens_respt(S, V, N21, Z12, b, 3);
    CHECK(r.W(0, 0) == doctest::Approx(0.723988439306358).epsilon(1e-12));
    CHECK(r.W(1, 0) == doctest::Approx(0.409248554913295).epsilon(1e-12));
    CHECK(r.W(0, 1) == doctest::Approx(0.752707581227437).epsilon(1e-12));
    CHECK(r.W(1, 1) == doctest::Approx(0.423104693140794).epsilon(1e-12));
    CHECK(r.WM[1](0, 0) == doctest::Approx(0.984104046242775).epsilon(1e-12));
    CHECK(r.WM[1](1, 1) == doctest::Approx(0.347725631768953).epsilon(1e-12));
    CHECK(r.WM[2](0, 0) == doctest::Approx(1.92052023121387).epsilon(1e-12));
    CHECK(r.WM[2](1, 1) == doctest::Approx(0.422815884476534).epsilon(1e-12));
    CHECK(r.m[0] == doctest::Approx(0.806524951688076).epsilon(1e-12));
    CHECK(r.m[1] == doctest::Approx(0.777083096510174).epsilon(1e-12));
    CHECK(r.Var[0] == doctest::Approx(0.752836815706844).epsilon(1e-12));
    CHECK(r.Var[1] == doctest::Approx(0.595186990139417).epsilon(1e-12));
    CHECK(r.p(0, 0) == doctest::Approx(0.449244060475162).epsilon(1e-12));
    CHECK(r.p(1, 0) == doctest::Approx(0.412072297374105).epsilon(1e-12));
    CHECK(r.p(1, 1) == doctest::Approx(0.420597931112879).epsilon(1e-12));
    CHECK(r.XN[0] == doctest::Approx(0.786631806297601).epsilon(1e-12));
    CHECK(r.XN[1] == doctest::Approx(0.314880072752075).epsilon(1e-12));

    // W(i,l) = w_i(l)/V(i,l) is the independent check of the t = 1 case of (4.5)
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t l = 0; l < 2; ++l)
            CHECK(r.W(i, l) == doctest::Approx(r.Wresid(i, l) / V(i, l)).epsilon(1e-11));
    // Var[W] >= 0 and the raw moments increase
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t l = 0; l < 2; ++l) {
            CHECK(r.WVar(i, l) >= -1e-14);
            CHECK(r.WM[1](i, l) >= r.WM[0](i, l) * r.WM[0](i, l) - 1e-14);
        }
}

TEST_CASE("pfqn_sens_respt single-server moments are exact rationals") {
    // With b = 1 the coefficients a_{t,0}(0) vanish and the marginal correction
    // disappears (Remark 4.2), so the moments are rational functions of S and V.
    const std::vector<Rational> S{Rational(1, 2), Rational(2, 5)};
    Matrix<Rational> V(2, 2);
    V(0, 0) = Rational(1); V(0, 1) = Rational(1);
    V(1, 0) = Rational(2); V(1, 1) = Rational(1);
    const std::vector<Rational> Z{Rational(1), Rational(2)};
    const auto r = pfqn_sens_respt(S, V, N21, Z, std::vector<int>(), 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t l = 0; l < 2; ++l)
            CHECK(r.W(i, l) * V(i, l) == r.Wresid(i, l));
    // the marginal p_i(0) is not used at b = 1 but must still be a probability
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.p(i, 0) >= Rational(0));
        CHECK(r.p(i, 0) <= Rational(1));
    }
}

TEST_CASE("pfqn_sens_ldmx_ec matches MATLAB and differentiates the load exactly") {
    // MATLAB: [EC,E,Ep,Lo,dEC,dE,dEp] = pfqn_sens_ldmx_ec([0.1 0],
    //         [0.5 0.3;0.4 0.6],[1 2 2;1 1 1]).
    const std::vector<double> lambda{0.1, 0.0};
    Matrix<double> D = demands2x2();
    Matrix<double> mu(2, 3);
    mu(0, 0) = 1.0; mu(0, 1) = 2.0; mu(0, 2) = 2.0;
    mu(1, 0) = 1.0; mu(1, 1) = 1.0; mu(1, 2) = 1.0;
    const auto e = pfqn_sens_ldmx_ec(lambda, D, mu);
    CHECK(e.Lo[0] == doctest::Approx(0.05).epsilon(1e-14));
    CHECK(e.Lo[1] == doctest::Approx(0.04).epsilon(1e-14));
    CHECK(e.EC(0, 0) == doctest::Approx(1.00062539086929).epsilon(1e-12));
    CHECK(e.EC(1, 0) == doctest::Approx(1.04166666666667).epsilon(1e-12));
    CHECK(e.EC(0, 1) == doctest::Approx(0.512820512820513).epsilon(1e-12));
    CHECK(e.EC(0, 2) == doctest::Approx(0.512820512820513).epsilon(1e-12));
    CHECK(e.E(0, 0) == doctest::Approx(1.05128205128205).epsilon(1e-12));
    CHECK(e.E(0, 1) == doctest::Approx(1.05193951347798).epsilon(1e-12));
    CHECK(e.E(1, 3) == doctest::Approx(1.17737569926698).epsilon(1e-12));
    CHECK(e.Eprime(0, 0) == doctest::Approx(1.02564102564103).epsilon(1e-12));
    CHECK(e.Eprime(1, 1) == doctest::Approx(1.08506944444444).epsilon(1e-12));
    CHECK(e.dEC(0, 0) == doctest::Approx(0.0250312793213085).epsilon(1e-11));
    CHECK(e.dEC(1, 0) == doctest::Approx(1.08506944444445).epsilon(1e-12));
    CHECK(e.dE(0, 0) == doctest::Approx(1.05193951347797).epsilon(1e-12));
    CHECK(e.dE(1, 3) == doctest::Approx(4.90573208027906).epsilon(1e-12));
    CHECK(e.dEprime(0, 0) == doctest::Approx(0.525969756738987).epsilon(1e-12));
    CHECK(e.dEprime(1, 3) == doctest::Approx(4.90573208027906).epsilon(1e-12));

    // dE/dLo against a central difference in the open load, obtained by moving
    // the open class's demand: Lo(i) = lambda(0) D(i,0).
    const double h = 1e-6;
    for (std::size_t i = 0; i < 2; ++i) {
        Matrix<double> Dp = D, Dm = D;
        Dp(i, 0) += h / lambda[0];
        Dm(i, 0) -= h / lambda[0];
        const auto hi = pfqn_sens_ldmx_ec(lambda, Dp, mu);
        const auto lo = pfqn_sens_ldmx_ec(lambda, Dm, mu);
        for (std::size_t n = 0; n < 4; ++n) {
            const double fd = (hi.E(i, n) - lo.E(i, n)) / (2 * h);
            CHECK(e.dE(i, n) == doctest::Approx(fd).epsilon(1e-6).scale(1e-6));
        }
        for (std::size_t n = 0; n < 3; ++n) {
            const double fd = (hi.EC(i, n) - lo.EC(i, n)) / (2 * h);
            CHECK(e.dEC(i, n) == doctest::Approx(fd).epsilon(1e-6).scale(1e-6));
        }
    }
}

TEST_CASE("pfqn_sens_mvaldmx matches the MATLAB reference on a mixed model") {
    // MATLAB: pfqn_sens_mvaldmx([0.2 0],[0.5 0.3;0.4 0.6],[Inf 2],[0 1],
    //         [1 2 2;1 1 1],[1;1]). Class 1 is open, class 2 closed.
    const std::vector<double> lambda{0.2, 0.0};
    Matrix<double> D = demands2x2();
    const std::vector<int> N{-1, 2};  // -1 marks the open class (MATLAB Inf)
    const std::vector<double> Z{0.0, 1.0};
    Matrix<double> mu(2, 3);
    mu(0, 0) = 1.0; mu(0, 1) = 2.0; mu(0, 2) = 2.0;
    mu(1, 0) = 1.0; mu(1, 1) = 1.0; mu(1, 2) = 1.0;
    const auto q = pfqn_sens_mvaldmx(lambda, D, N, Z, mu);
    CHECK(q.XN[0] == doctest::Approx(0.2).epsilon(1e-14));
    CHECK(q.XN[1] == doctest::Approx(0.920372508176005).epsilon(1e-12));
    CHECK(q.QN(0, 0) == doctest::Approx(0.102714496911566).epsilon(1e-12));
    CHECK(q.QN(1, 0) == doctest::Approx(0.156581941889852).epsilon(1e-12));
    CHECK(q.QN(0, 1) == doctest::Approx(0.278935160090696).epsilon(1e-12));
    CHECK(q.QN(1, 1) == doctest::Approx(0.800692331733299).epsilon(1e-12));
    CHECK(q.UN(0, 0) == doctest::Approx(0.1).epsilon(1e-12));
    CHECK(q.UN(1, 1) == doctest::Approx(0.552223504905603).epsilon(1e-12));
    CHECK(q.CN(0, 0) == doctest::Approx(0.513572484557832).epsilon(1e-12));
    CHECK(q.CN(1, 1) == doctest::Approx(0.869965502685551).epsilon(1e-12));
    CHECK(q.QVar(0, 0) == doctest::Approx(0.104645615644411).epsilon(1e-11));
    CHECK(q.QVar(1, 0) == doctest::Approx(0.174435818915237).epsilon(1e-11));
    CHECK(q.QVar(0, 1) == doctest::Approx(0.245889699454148).epsilon(1e-11));
    CHECK(q.QVar(1, 1) == doctest::Approx(0.560482904873818).epsilon(1e-11));
    CHECK(q.QTotVar[0] == doctest::Approx(0.356668649301303).epsilon(1e-11));
    CHECK(q.QTotVar[1] == doctest::Approx(0.832394011593197).epsilon(1e-11));
    CHECK(std::fabs(q.QCovAsym) < 1e-14);

    // The utilization law holds for the open class, whose throughput is its
    // arrival rate, whatever the load dependence at the station.
    CHECK(q.QN(0, 0) == doctest::Approx(q.XN[0] * q.CN(0, 0)).epsilon(1e-12));
    CHECK(q.QN(1, 0) == doctest::Approx(q.XN[0] * q.CN(1, 0)).epsilon(1e-12));
    CHECK(q.QN(0, 1) == doctest::Approx(q.XN[1] * q.CN(0, 1)).epsilon(1e-12));
    // the closed class's queue lengths sum to its population
    CHECK(q.QN(0, 1) + q.QN(1, 1) + q.XN[1] * Z[1] == doctest::Approx(2.0).epsilon(1e-10));
}
