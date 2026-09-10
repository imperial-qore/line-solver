/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Markovian fluid queues: the fundamental matrices, the general stationary
 * solution, and the two sojourn-time representations built on them
 * (line/api/mam/mfq_solve.h, mfq_sojourn.h, mfq_fluflu_sojourn.h).
 *
 * Oracles, in the order the task prescribes.
 *  (a) A hand-solvable closed form. On Q = [-2 2; 1 -1], Rin = diag(3,1),
 *      Rout = diag(2,2) the normalized blocks are scalars, the Riccati
 *      equation reduces to Psi^2 - 3 Psi + 2 = 0 whose minimal non-negative
 *      root is Psi = 1, and the whole solution follows by hand: K = -1,
 *      clo = [1 1], mass0 = [0, 1/3], ini = 1/3, mean level 2/3, fluid arrival
 *      rate 5/3 and mean sojourn time exactly 2/5. Every one of those is
 *      asserted, so the test pins the algebra and not just a number.
 *  (b) Invariants. The Riccati residual Fpm + Fpp Psi + Psi Fmm + Psi Fmp Psi
 *      must vanish; the stationary law must have total mass one,
 *      mass0 e + ini (-K)^-1 clo e = 1; and Little's law for fluid,
 *      E[sojourn] = E[level]/lambda, must hold. The ME and PH representations
 *      must describe the SAME distribution, which is a real check because they
 *      are built through disjoint code paths (a triangular similarity against
 *      a diagonal one).
 *  (c) MATLAB, digit for digit, on a two-state instance, a three-state
 *      instance with a ZERO-drift state (which exercises the censoring of the
 *      zero-drift block) and the same instance with an explicit level-zero
 *      boundary generator Q0.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mfq_fluflu_sojourn.h"
#include "line/api/mam/mfq_ld_distr.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mam/mfq_ld_mean.h"
#include "line/api/mam/mfq_sojourn.h"
#include "line/util/expm.h"

using line::Matrix;
using line::Real50;
using line::mam::MeRepresentation;
using line::mam::mfq_fluflu_sojourn;
using line::mam::mfq_fundamental;
using line::mam::mfq_general_solve;
using line::mam::mfq_sojourn;

namespace {

template <class T>
Matrix<T> mat(const std::vector<std::vector<double>>& a) {
    Matrix<T> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j)
            m(i, j) = line::num_traits<T>::from_double(a[i][j]);
    return m;
}

template <class T>
Matrix<T> diag(const std::vector<double>& v) {
    Matrix<T> m(v.size(), v.size(), line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < v.size(); ++i) m(i, i) = line::num_traits<T>::from_double(v[i]);
    return m;
}

/** k-th raw moment of the (possibly defective) ME law: (-1)^k k! alpha A^-k e. */
double me_moment(const MeRepresentation<double>& r, unsigned k) {
    Matrix<double> P = line::eye<double>(r.A.rows());
    const Matrix<double> Ai = line::inverse(r.A);
    for (unsigned i = 0; i < k; ++i) P = line::matmul(P, Ai);
    const std::vector<double> v = line::mulvec(P, line::ones<double>(r.A.rows()));
    double s = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) s += r.alpha[i] * v[i];
    double f = 1.0;
    for (unsigned i = 2; i <= k; ++i) f *= static_cast<double>(i);
    return (k % 2 == 0 ? 1.0 : -1.0) * f * s;
}

double me_cdf(const MeRepresentation<double>& r, double t) {
    const Matrix<double> E = line::expm(r.A, t);
    const std::vector<double> v = line::mulvec(E, line::ones<double>(r.A.rows()));
    double s = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) s += r.alpha[i] * v[i];
    return 1.0 - s;
}

const std::vector<std::vector<double>> kQ2 = {{-2.0, 2.0}, {1.0, -1.0}};
const std::vector<std::vector<double>> kQ3 = {{-3.0, 2.0, 1.0}, {1.0, -2.0, 1.0},
                                              {2.0, 1.0, -3.0}};

}  // namespace

// ---------------------------------------------------------------------------
// (a) the hand-solvable instance
// ---------------------------------------------------------------------------

TEST_CASE("the two-state fluid queue is solvable by hand and the port matches") {
    // Normalized blocks are scalars: Fpp = -2, Fpm = 2, Fmp = 1, Fmm = -1, so
    // Psi solves Psi^2 - 3 Psi + 2 = 0 and the minimal non-negative root is 1.
    const Matrix<double> Fpp = mat<double>({{-2.0}});
    const Matrix<double> Fpm = mat<double>({{2.0}});
    const Matrix<double> Fmp = mat<double>({{1.0}});
    const Matrix<double> Fmm = mat<double>({{-1.0}});
    const auto ff = mfq_fundamental(Fpp, Fpm, Fmp, Fmm, 1e-14, 150u,
                                    line::mam::RiccatiMethod::ADDA);
    CHECK(ff.converged);
    CHECK(ff.Psi(0, 0) == doctest::Approx(1.0).epsilon(1e-13));
    CHECK(ff.K(0, 0) == doctest::Approx(-1.0).epsilon(1e-13));
    CHECK(ff.U(0, 0) == doctest::Approx(0.0).epsilon(1e-13));

    // and the full stationary solution, all of it exact rationals.
    const Matrix<double> Q = mat<double>(kQ2);
    const Matrix<double> Rnet = diag<double>({1.0, -1.0});
    const auto s = mfq_general_solve(Q, Rnet);
    CHECK(s.mass0[0] == doctest::Approx(0.0).epsilon(1e-13));
    CHECK(s.mass0[1] == doctest::Approx(1.0 / 3.0).epsilon(1e-13));
    CHECK(s.ini[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-13));
    CHECK(s.K(0, 0) == doctest::Approx(-1.0).epsilon(1e-13));
    CHECK(s.clo(0, 0) == doctest::Approx(1.0).epsilon(1e-13));
    CHECK(s.clo(0, 1) == doctest::Approx(1.0).epsilon(1e-13));

    // The sojourn time then has mean E[level]/lambda = (2/3)/(5/3) = 2/5.
    const auto me = mfq_sojourn(Q, diag<double>({3.0, 1.0}), diag<double>({2.0, 2.0}));
    CHECK(me_moment(me, 1) == doctest::Approx(0.4).epsilon(1e-12));
    // alpha sums to 1 - P(the level is empty on arrival): the sojourn law is
    // defective, the missing 0.2 being the atom at zero.
    double asum = 0.0;
    for (double v : me.alpha) asum += v;
    CHECK(asum == doctest::Approx(0.8).epsilon(1e-12));
    CHECK(me_cdf(me, 0.0) == doctest::Approx(1.0 - asum).epsilon(1e-12));
}

TEST_CASE("SDA reaches the same Psi as ADDA") {
    // The same doubling iteration with a common shift: a different sequence of
    // iterates converging to the same minimal solution.
    const Matrix<double> Fpp = mat<double>({{-2.5, 0.5}, {0.4, -1.6}});
    const Matrix<double> Fpm = mat<double>({{1.6, 0.4}, {0.7, 0.5}});
    const Matrix<double> Fmp = mat<double>({{0.9, 0.3}, {0.2, 1.1}});
    const Matrix<double> Fmm = mat<double>({{-1.9, 0.7}, {0.6, -1.9}});
    const auto a = mfq_fundamental(Fpp, Fpm, Fmp, Fmm, 1e-14, 150u,
                                   line::mam::RiccatiMethod::ADDA);
    const auto d = mfq_fundamental(Fpp, Fpm, Fmp, Fmm, 1e-14, 150u,
                                   line::mam::RiccatiMethod::SDA);
    REQUIRE(a.converged);
    REQUIRE(d.converged);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(a.Psi(i, j) == doctest::Approx(d.Psi(i, j)).epsilon(1e-11));
    // The Riccati residual, the definition of Psi, must vanish for both.
    for (const auto* r : {&a, &d}) {
        const Matrix<double> res = line::mam::mfq_detail::add(
            line::mam::mfq_detail::add(Fpm, line::matmul(Fpp, r->Psi)),
            line::mam::mfq_detail::add(line::matmul(r->Psi, Fmm),
                                       line::matmul(line::matmul(r->Psi, Fmp), r->Psi)));
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(res(i, j)) < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// (b) invariants
// ---------------------------------------------------------------------------

TEST_CASE("the stationary fluid law has total mass one and satisfies Little's law") {
    struct Case {
        std::vector<std::vector<double>> Q;
        std::vector<double> rin, rout;
    };
    const std::vector<Case> cases = {{kQ2, {3.0, 1.0}, {2.0, 2.0}},
                                     {kQ3, {2.0, 1.0, 2.0}, {1.0, 3.0, 2.0}},
                                     {kQ3, {4.0, 1.0, 1.0}, {2.0, 3.0, 2.0}}};
    for (const Case& c : cases) {
        const Matrix<double> Q = mat<double>(c.Q);
        const Matrix<double> Rin = diag<double>(c.rin);
        const Matrix<double> Rout = diag<double>(c.rout);
        const std::size_t N = c.rin.size();
        Matrix<double> Rnet = Rin;
        for (std::size_t i = 0; i < N; ++i) Rnet(i, i) -= Rout(i, i);
        const auto s = mfq_general_solve(Q, Rnet);
        const Matrix<double> negKinv = line::inverse(line::mam::mfq_detail::scale(s.K, -1.0));

        // Total mass: the atom at zero plus the density integrated over the level.
        double mass = 0.0;
        for (double v : s.mass0) mass += v;
        const std::vector<double> iniKi = line::vecmul(s.ini, negKinv);
        const std::vector<double> tail = line::vecmul(iniKi, s.clo);
        for (double v : tail) mass += v;
        CHECK(mass == doctest::Approx(1.0).epsilon(1e-11));

        // Mean level: ini (-K)^-2 clo e. Fluid arrival rate: mass0 Rin summed
        // plus the density weighted by Rin.
        double level = 0.0;
        {
            const std::vector<double> t = line::vecmul(line::vecmul(iniKi, negKinv), s.clo);
            for (double v : t) level += v;
        }
        double lambda = 0.0;
        {
            const std::vector<double> a = line::vecmul(s.mass0, Rin);
            const std::vector<double> b = line::vecmul(tail, Rin);
            for (std::size_t j = 0; j < N; ++j) lambda += a[j] + b[j];
        }
        const auto me = mfq_sojourn(Q, Rin, Rout);
        INFO("N = ", N);
        CHECK(me_moment(me, 1) == doctest::Approx(level / lambda).epsilon(1e-10));
    }
}

TEST_CASE("the ME and PH representations describe the same sojourn law") {
    const Matrix<double> Q = mat<double>(kQ3);
    const Matrix<double> Rin = diag<double>({2.0, 1.0, 2.0});
    const Matrix<double> Rout = diag<double>({1.0, 3.0, 2.0});
    const auto me = mfq_sojourn(Q, Rin, Rout);
    const auto ph = mfq_sojourn(Q, Rin, Rout, Matrix<double>(0, 0), true, 1e-14);
    for (unsigned k = 1; k <= 3; ++k) {
        INFO("moment ", k);
        CHECK(me_moment(me, k) == doctest::Approx(me_moment(ph, k)).epsilon(1e-10));
    }
    for (double t : {0.1, 0.5, 1.0, 2.5}) {
        INFO("t = ", t);
        CHECK(me_cdf(me, t) == doctest::Approx(me_cdf(ph, t)).epsilon(1e-10));
    }
    // A CDF: non-decreasing, starting at the atom, approaching one.
    CHECK(me_cdf(me, 0.0) < me_cdf(me, 0.5));
    CHECK(me_cdf(me, 0.5) < me_cdf(me, 2.5));
    CHECK(me_cdf(me, 40.0) == doctest::Approx(1.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("mfq_sojourn agrees with MATLAB") {
    SUBCASE("two states, no zero-drift state") {
        const Matrix<double> Q = mat<double>(kQ2);
        const auto me = mfq_sojourn(Q, diag<double>({3.0, 1.0}), diag<double>({2.0, 2.0}));
        CHECK(me.A.rows() == 2u);
        double asum = 0.0;
        for (double v : me.alpha) asum += v;
        CHECK(asum == doctest::Approx(0.8).epsilon(1e-13));
        CHECK(me_moment(me, 1) == doctest::Approx(0.4).epsilon(1e-13));
        CHECK(me_moment(me, 2) == doctest::Approx(0.399999999999999).epsilon(1e-12));
        CHECK(me_moment(me, 3) == doctest::Approx(0.599999999999999).epsilon(1e-12));
        CHECK(me_cdf(me, 0.5) == doctest::Approx(0.705696447062846).epsilon(1e-12));
    }
    SUBCASE("three states, one of them zero drift") {
        // Rin(3) = Rout(3) = 2, so state 3 is censored out before the Riccati
        // solve. mass0, ini, K and clo are all exact rationals here:
        // (0, 9/32, 3/32), 15/32, -3/2 and (1, 1/2, 1/2).
        const Matrix<double> Q = mat<double>(kQ3);
        const Matrix<double> Rin = diag<double>({2.0, 1.0, 2.0});
        const Matrix<double> Rout = diag<double>({1.0, 3.0, 2.0});
        Matrix<double> Rnet = Rin;
        for (std::size_t i = 0; i < 3; ++i) Rnet(i, i) -= Rout(i, i);
        const auto s = mfq_general_solve(Q, Rnet);
        CHECK(s.mass0[0] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(s.mass0[1] == doctest::Approx(9.0 / 32.0).epsilon(1e-12));
        CHECK(s.mass0[2] == doctest::Approx(3.0 / 32.0).epsilon(1e-12));
        CHECK(s.ini[0] == doctest::Approx(15.0 / 32.0).epsilon(1e-12));
        CHECK(s.K(0, 0) == doctest::Approx(-1.5).epsilon(1e-12));
        CHECK(s.clo(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(s.clo(0, 1) == doctest::Approx(0.5).epsilon(1e-12));
        CHECK(s.clo(0, 2) == doctest::Approx(0.5).epsilon(1e-12));

        const auto me = mfq_sojourn(Q, Rin, Rout);
        double asum = 0.0;
        for (double v : me.alpha) asum += v;
        CHECK(asum == doctest::Approx(0.699999999999999).epsilon(1e-12));
        CHECK(me_moment(me, 1) == doctest::Approx(0.266666666666666).epsilon(1e-12));
        CHECK(me_moment(me, 2) == doctest::Approx(0.195981087470448).epsilon(1e-12));
        CHECK(me_cdf(me, 1.0) == doctest::Approx(0.952444761106016).epsilon(1e-12));
    }
    SUBCASE("an explicit level-zero boundary generator") {
        const Matrix<double> Q = mat<double>(kQ3);
        const Matrix<double> Q0 = mat<double>({{-1.0, 1.0, 0.0}, {2.0, -3.0, 1.0},
                                               {1.0, 1.0, -2.0}});
        const auto me = mfq_sojourn(Q, diag<double>({2.0, 1.0, 2.0}),
                                    diag<double>({1.0, 3.0, 2.0}), Q0, false, 1e-14);
        double asum = 0.0;
        for (double v : me.alpha) asum += v;
        CHECK(asum == doctest::Approx(0.744680851063829).epsilon(1e-12));
        CHECK(me_moment(me, 1) == doctest::Approx(0.283687943262411).epsilon(1e-12));
        CHECK(me_moment(me, 2) == doctest::Approx(0.208490518585583).epsilon(1e-12));
        CHECK(me_cdf(me, 1.0) == doctest::Approx(0.949409320325548).epsilon(1e-12));
    }
}

TEST_CASE("mfq_fluflu_sojourn agrees with MATLAB in both service conventions") {
    const Matrix<double> Qin = mat<double>(kQ2);
    const Matrix<double> Rin = diag<double>({3.0, 1.0});
    const Matrix<double> Qout = mat<double>({{-1.0, 1.0}, {2.0, -2.0}});
    const Matrix<double> Rout = diag<double>({4.0, 1.0});
    SUBCASE("service continues at a zero server level") {
        const auto me = mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, false);
        const auto ph = mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, false, true, 1e-14);
        CHECK(me_moment(me, 1) == doctest::Approx(0.147480963363268).epsilon(1e-12));
        CHECK(me_moment(me, 2) == doctest::Approx(0.108275686460907).epsilon(1e-12));
        CHECK(me_cdf(me, 0.5) == doctest::Approx(0.897097984978141).epsilon(1e-12));
        CHECK(me_moment(ph, 1) == doctest::Approx(me_moment(me, 1)).epsilon(1e-12));
        CHECK(me_cdf(ph, 0.5) == doctest::Approx(me_cdf(me, 0.5)).epsilon(1e-12));
    }
    SUBCASE("service stops at a zero server level") {
        const auto me = mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, true);
        CHECK(me_moment(me, 1) == doctest::Approx(0.103036518918824).epsilon(1e-12));
        CHECK(me_moment(me, 2) == doctest::Approx(0.0756460329662891).epsilon(1e-12));
        CHECK(me_cdf(me, 0.5) == doctest::Approx(0.928108244102875).epsilon(1e-12));
        // Stopping service when the server's own fluid runs out cannot make a
        // drop leave later on average in this instance: the weighting removes
        // the mass that arrives while the server is idle.
        const auto go = mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, false);
        CHECK(me_moment(me, 1) < me_moment(go, 1));
    }
}

TEST_CASE("only an all-up-drift model is rejected, and an all-down one is solved") {
    const Matrix<double> Q = mat<double>(kQ2);
    // Every state drifting UP: the level grows without bound, so there is no
    // stationary law. The MATLAB reference does not check this and dies inside
    // its boundary solve with "Incorrect dimensions for matrix multiplication";
    // the port refuses by name instead.
    CHECK_THROWS_AS(mfq_general_solve(Q, diag<double>({1.0, 2.0})), line::InputError);

    // Every state drifting DOWN is a perfectly good model, not an error: the
    // level is identically zero and the answer is the degenerate law, a point
    // mass at zero equal to the stationary distribution of Q and no density.
    // A fluid priority queue reaches this whenever the states feeding a class
    // never together exceed the service rate, so rejecting it would have made
    // mfq_prio_queue unusable. MATLAB returns mass0 = (1/3, 2/3) here, an
    // empty ini and a 0 x 0 K.
    const auto s = mfq_general_solve(Q, diag<double>({-1.0, -2.0}));
    CHECK(s.ini.empty());
    CHECK(s.K.rows() == 0u);
    CHECK(s.clo.rows() == 0u);
    REQUIRE(s.mass0.size() == 2u);
    CHECK(s.mass0[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-13));
    CHECK(s.mass0[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-13));
    // and that mass0 really is the stationary law of Q, not a coincidence.
    const std::vector<double> pi = line::mc::ctmc_solve(Q);
    for (std::size_t j = 0; j < 2; ++j)
        CHECK(s.mass0[j] == doctest::Approx(pi[j]).epsilon(1e-13));
}

TEST_CASE("mfq_sojourn at Real50 reproduces its own double result") {
    const Matrix<double> Qd = mat<double>(kQ2);
    const auto md = mfq_sojourn(Qd, diag<double>({3.0, 1.0}), diag<double>({2.0, 2.0}));
    const auto mq = mfq_sojourn(mat<Real50>(kQ2), diag<Real50>({3.0, 1.0}),
                                diag<Real50>({2.0, 2.0}));
    REQUIRE(mq.A.rows() == md.A.rows());
    for (std::size_t i = 0; i < md.A.rows(); ++i) {
        CHECK(static_cast<double>(mq.alpha[i]) == doctest::Approx(md.alpha[i]).epsilon(1e-11));
        for (std::size_t j = 0; j < md.A.rows(); ++j)
            CHECK(static_cast<double>(mq.A(i, j)) == doctest::Approx(md.A(i, j)).epsilon(1e-11));
    }
}

// ---------------------------------------------------------------------------
// mfq_ld_mean
// ---------------------------------------------------------------------------

TEST_CASE("mfq_ld_mean integrates a singular regime exactly") {
    // KF = 0 makes the two regime integrals J0 = L I and J1 = (L^2/2) I, which
    // no inversion-based formula could reach: it is exactly the case the
    // nilpotent augmentation exists for. With L = 2, masses (0.2, 0.1) at
    // level 0 and (0.1, 0.1) at level 2, iniF = (0.3, 0.2), cloF = I and a
    // zero backward part, E[X] = 2 (0.1 + 0.1) + (0.3 + 0.2) 2 = 1.4 exactly.
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {2.0};
    b.masses = {{0.2, 0.1}, {0.1, 0.1}};
    b.iniF = {{0.3, 0.2}};
    b.KF = {Matrix<double>(2, 2, 0.0)};
    b.cloF = {line::eye<double>(2)};
    b.iniB = {{0.0, 0.0}};
    b.KB = {line::mam::mfq_detail::scale(line::eye<double>(2), -1.0)};
    b.cloB = {Matrix<double>(2, 2, 0.0)};
    CHECK(line::mam::mfq_ld_mean(b) == doctest::Approx(1.4).epsilon(1e-13));
}

TEST_CASE("mfq_ld_mean recovers the closed-form mean level of a homogeneous queue") {
    // One regime stretched over [0, 200] with the blocks of the hand-solvable
    // two-state queue and no backward part: the formula must converge to the
    // homogeneous closed form ini (-K)^-2 clo e = 2/3.
    const Matrix<double> Q = mat<double>(kQ2);
    const auto s = mfq_general_solve(Q, diag<double>({1.0, -1.0}));
    const std::size_t Np = s.ini.size();
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {200.0};
    b.masses = {s.mass0, std::vector<double>(s.mass0.size(), 0.0)};
    b.iniF = {s.ini};
    b.KF = {s.K};
    b.cloF = {s.clo};
    b.iniB = {std::vector<double>(Np, 0.0)};
    b.KB = {line::mam::mfq_detail::scale(line::eye<double>(Np), -1.0)};
    b.cloB = {Matrix<double>(Np, s.clo.cols(), 0.0)};
    const double closed = 2.0 / 3.0;
    CHECK(line::mam::mfq_ld_mean(b) == doctest::Approx(closed).epsilon(1e-12));
    // MATLAB on the same blocks: 0.666666666666629.
    CHECK(line::mam::mfq_ld_mean(b) == doctest::Approx(0.666666666666629).epsilon(1e-12));
}

TEST_CASE("mfq_ld_mean agrees with MATLAB on a two-regime instance") {
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {1.5, 4.0};
    b.masses = {{0.1, 0.05}, {0.02, 0.03}, {0.04, 0.01}};
    b.iniF = {{0.3, 0.2}, {0.15, 0.25}};
    b.KF = {mat<double>({{-1.5, 0.5}, {0.25, -2.0}}), mat<double>({{-0.8, 0.0}, {0.3, -1.1}})};
    b.cloF = {mat<double>({{1.0, 0.5}, {0.25, 1.0}}), mat<double>({{0.6, 0.4}, {0.2, 0.9}})};
    b.iniB = {{0.12, 0.08}, {0.05, 0.11}};
    b.KB = {mat<double>({{-2.2, 0.3}, {0.1, -1.7}}), mat<double>({{-1.3, 0.2}, {0.4, -0.9}})};
    b.cloB = {mat<double>({{0.9, 0.1}, {0.3, 0.7}}), mat<double>({{1.1, 0.2}, {0.05, 0.8}})};
    CHECK(line::mam::mfq_ld_mean(b) == doctest::Approx(2.27771820045627).epsilon(1e-12));
}

TEST_CASE("mfq_ld_mean rejects inconsistent block lists") {
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {1.0};
    b.masses = {{0.5, 0.5}};  // needs K+1 = 2
    b.iniF = {{0.1, 0.1}};
    b.KF = {line::eye<double>(2)};
    b.cloF = {line::eye<double>(2)};
    b.iniB = {{0.0, 0.0}};
    b.KB = {line::eye<double>(2)};
    b.cloB = {line::eye<double>(2)};
    CHECK_THROWS_AS(line::mam::mfq_ld_mean(b), line::InputError);
    b.Thr.clear();
    CHECK_THROWS_AS(line::mam::mfq_ld_mean(b), line::InputError);
}

// ---------------------------------------------------------------------------
// mfq_ld_distr
// ---------------------------------------------------------------------------

namespace {

/** The two-regime block set used against MATLAB, KF singular by construction. */
line::mam::LevelDependentFluidBlocks<double> ld_blocks() {
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {1.5, 4.0};
    b.masses = {{0.10, 0.05}, {0.02, 0.03}, {0.04, 0.01}};
    b.iniF = {{0.30, 0.20}, {0.15, 0.25}};
    // Zero row sums: these are the directions the deflation has to handle.
    b.KF = {mat<double>({{-1.5, 1.5}, {0.25, -0.25}}),
            mat<double>({{-0.8, 0.8}, {0.3, -0.3}})};
    b.cloF = {mat<double>({{1.0, 0.5}, {0.25, 1.0}}), mat<double>({{0.6, 0.4}, {0.2, 0.9}})};
    b.iniB = {{0.12, 0.08}, {0.05, 0.11}};
    b.KB = {mat<double>({{-2.2, 0.3}, {0.1, -1.7}}), mat<double>({{-1.3, 0.2}, {0.4, -0.9}})};
    b.cloB = {mat<double>({{0.9, 0.1}, {0.3, 0.7}}), mat<double>({{1.1, 0.2}, {0.05, 0.8}})};
    return b;
}

}  // namespace

TEST_CASE("mfq_ld_distr integrates a singular exponent against its closed form") {
    // One regime, no backward part, KF = [[-a, a], [b, -b]] with a = 1.5,
    // b = 0.25 and s = a + b. Then exp(KF u) = P + (I - P) e^{-s u} with
    // P = e pi, pi = (b, a)/s, so
    //   int_0^L exp(KF u) du = P L + (I - P)(1 - e^{-s L})/s
    // in closed form. KF is singular, so the port must take the deflation
    // route to get there, and this is the check that it does.
    const double a = 1.5, bb = 0.25, s = a + bb, L = 0.9;
    line::mam::LevelDependentFluidBlocks<double> b;
    b.Thr = {3.0};
    b.masses = {{0.2, 0.1}, {0.0, 0.0}};
    b.iniF = {{0.3, 0.2}};
    b.KF = {mat<double>({{-a, a}, {bb, -bb}})};
    b.cloF = {line::eye<double>(2)};
    b.iniB = {{0.0, 0.0}};
    b.KB = {mat<double>({{-2.2, 0.3}, {0.1, -1.7}})};
    b.cloB = {Matrix<double>(2, 2, 0.0)};

    // Closed-form integral.
    const double pi0 = bb / s, pi1 = a / s;
    Matrix<double> P(2, 2, 0.0);
    P(0, 0) = pi0; P(0, 1) = pi1; P(1, 0) = pi0; P(1, 1) = pi1;
    Matrix<double> J(2, 2, 0.0);
    const double decay = (1.0 - std::exp(-s * L)) / s;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            J(i, j) = P(i, j) * L + ((i == j ? 1.0 : 0.0) - P(i, j)) * decay;
    const std::vector<double> want = line::vecmul(b.iniF[0], J);

    const auto got = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdf,
                                             std::vector<double>{L});
    // Cdf at L excludes no atom here (L is not a threshold), and the mass at
    // level 0 is included because L > 0.
    for (std::size_t j = 0; j < 2; ++j)
        CHECK(got[0][j] == doctest::Approx(want[j] + b.masses[0][j]).epsilon(1e-11));
}

TEST_CASE("Cdfm minus Cdf is exactly the atom at the point, and zero elsewhere") {
    const auto b = ld_blocks();
    // 1.5 and 4.0 are thresholds and carry atoms; 0 carries the level-zero
    // atom; 0.5, 2.0 and 5.0 carry none.
    const std::vector<double> pts = {0.0, 0.5, 1.5, 2.0, 4.0, 5.0};
    const auto lo = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdf, pts);
    const auto hi = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdfm, pts);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        INFO("p = ", pts[i]);
        std::vector<double> atom(2, 0.0);
        if (pts[i] == 0.0) atom = b.masses[0];
        else if (pts[i] == 1.5) atom = b.masses[1];
        else if (pts[i] == 4.0) atom = b.masses[2];
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(hi[i][j] - lo[i][j] == doctest::Approx(atom[j]).epsilon(1e-11));
    }
}

TEST_CASE("Pdf is the derivative of Cdf and Pdfd the derivative of Pdf") {
    const auto b = ld_blocks();
    // Central differences well inside a regime, away from the atoms.
    const double h = 1e-5;
    for (double p : {0.7, 2.6, 3.4}) {
        INFO("p = ", p);
        const auto up = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdf,
                                                std::vector<double>{p + h});
        const auto dn = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdf,
                                                std::vector<double>{p - h});
        const auto pdf = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Pdf,
                                                 std::vector<double>{p});
        const auto pu = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Pdf,
                                                std::vector<double>{p + h});
        const auto pd = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Pdf,
                                                std::vector<double>{p - h});
        const auto pdfd = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Pdfd,
                                                  std::vector<double>{p});
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK((up[0][j] - dn[0][j]) / (2 * h) == doctest::Approx(pdf[0][j]).epsilon(1e-7));
            CHECK((pu[0][j] - pd[0][j]) / (2 * h) == doctest::Approx(pdfd[0][j]).epsilon(1e-6));
        }
    }
}

TEST_CASE("mfq_ld_distr and mfq_ld_mean agree on the same blocks") {
    // Two ported functions, two disjoint code paths: mfq_ld_mean integrates
    // x pi(x) in closed form through a nilpotent augmentation, while this
    // quadrature reads pi(x) off mfq_ld_distr on a fine grid and adds the
    // atoms by hand. They must land on the same E[X].
    //
    // The quadrature is run PER REGIME, not over one grid spanning both. The
    // density is discontinuous at a threshold, so a panel straddling 1.5
    // makes the trapezoidal rule first order there and caps the agreement at
    // about 1e-5 however fine the grid is. Splitting at the kink restores
    // second order and the two agree to 1e-9. That is a property of the test's
    // quadrature, not of either ported function.
    const auto b = ld_blocks();
    const std::vector<double> edge = {0.0, 1.5, 4.0};
    const double inset = 1e-9;
    const std::size_t M = 8001;
    double quad = 0.0;
    for (std::size_t r = 0; r + 1 < edge.size(); ++r) {
        const double lo = edge[r] + inset, hi = edge[r + 1] - inset;
        std::vector<double> grid(M);
        for (std::size_t i = 0; i < M; ++i)
            grid[i] = lo + (hi - lo) * static_cast<double>(i) / (M - 1);
        const auto pdf = line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Pdf, grid);
        for (std::size_t i = 1; i < M; ++i) {
            double a0 = 0.0, a1 = 0.0;
            for (std::size_t j = 0; j < 2; ++j) {
                a0 += grid[i - 1] * pdf[i - 1][j];
                a1 += grid[i] * pdf[i][j];
            }
            quad += 0.5 * (a0 + a1) * (grid[i] - grid[i - 1]);
        }
    }
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t j = 0; j < 2; ++j) quad += edge[k] * b.masses[k][j];
    CHECK(quad == doctest::Approx(line::mam::mfq_ld_mean(b)).epsilon(1e-9));
}

TEST_CASE("mfq_ld_distr agrees with MATLAB on all four functionals") {
    const auto b = ld_blocks();
    const std::vector<double> pts = {0.0, 0.5, 1.5, 2.0, 4.0, 5.0};
    const std::vector<std::vector<double>> refPdf = {
        {0.357583809516312, 0.357118359647807},  {0.269255532860157, 0.432115561514628},
        {0.152452774972538, 0.299614970984423},  {0.15080604934843, 0.314992417545488},
        {0.185182455910655, 0.402146930111681},  {0.197231408821993, 0.517976005359482}};
    const std::vector<std::vector<double>> refPdfd = {
        {-0.286050849971734, 0.211160900017862},
        {-0.0888929770349855, 0.107123030737703},
        {-0.00910085774826709, 0.0335553235855513},
        {0.00197517083130355, 0.0291481139143934},
        {0.0263992984982793, 0.0768383768771509},
        {-0.0218252485555463, 0.16753730491858}};
    const std::vector<std::vector<double>> refCdf = {
        {0.0, 0.0},                             {0.252655428003363, 0.249447555265573},
        {0.523079295668364, 0.727092911423565}, {0.618664075346444, 0.91083584592206},
        {0.945960416686505, 1.61259122029769},  {1.14090574936178, 2.06518509741431}};
    const std::vector<std::vector<double>> refCdfm = {
        {0.1, 0.05},                            {0.252655428003363, 0.249447555265573},
        {0.543079295668364, 0.757092911423565}, {0.618664075346444, 0.91083584592206},
        {0.985960416686505, 1.62259122029769},  {1.14090574936178, 2.06518509741431}};
    struct Item {
        line::mam::FluidDistrKind kind;
        const std::vector<std::vector<double>>* ref;
    };
    const std::vector<Item> items = {{line::mam::FluidDistrKind::Pdf, &refPdf},
                                     {line::mam::FluidDistrKind::Pdfd, &refPdfd},
                                     {line::mam::FluidDistrKind::Cdf, &refCdf},
                                     {line::mam::FluidDistrKind::Cdfm, &refCdfm}};
    for (const Item& it : items) {
        const auto got = line::mam::mfq_ld_distr(b, it.kind, pts);
        REQUIRE(got.size() == pts.size());
        for (std::size_t i = 0; i < pts.size(); ++i)
            for (std::size_t j = 0; j < 2; ++j) {
                INFO("point index ", i, " state ", j);
                CHECK(got[i][j] == doctest::Approx((*it.ref)[i][j]).epsilon(1e-11));
            }
    }
}

TEST_CASE("mfq_ld_distr rejects a negative evaluation point") {
    const auto b = ld_blocks();
    CHECK_THROWS_AS(
        line::mam::mfq_ld_distr(b, line::mam::FluidDistrKind::Cdf, std::vector<double>{-0.1}),
        line::InputError);
}
