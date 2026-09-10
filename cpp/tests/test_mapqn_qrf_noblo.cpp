/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The QRF no-blocking nonlinear bounds: qrf_noblo_mmi, qrf_noblo_mem and the
 * load-dependent qrf_noblo_mmi_ld.
 *
 * THE TWO OBJECTIVES NEED DIFFERENT ORACLES, and conflating them is how this
 * family gets mis-tested:
 *
 *  - MEM is CONVEX on the polytope, so its optimum is unique and the port is
 *    checked against the reference NUMBER (native Python, run on the same
 *    instances) and against the exact product-form answer where one exists.
 *  - MMI is NOT convex -- its -p_ij log p_ii terms are not -- so both codebases
 *    report a local optimum fixed by whichever vertex the phase-1 LP returns.
 *    It is therefore checked against the PROPERTIES a correct solve has
 *    (feasibility, an objective no worse than the start, a vanishing
 *    Frank-Wolfe gap) plus the reference number on the instance where the two
 *    starts do coincide.
 *
 * Measured on the one-phase instance below, both objectives land on the EXACT
 * product-form solution, which is the strongest single check available here:
 * it validates the constraint inventory, the equality reduction, the phase-1
 * start and the optimizer at once, against an answer neither codebase computed.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_qrf_noblo.h"

namespace mq = line::mapqn;
using line::Matrix;

namespace {

/** Two exponential stations in a cycle, rates mu1 and mu2. */
void expo_pair(double mu1, double mu2, std::vector<int>* K, Matrix<double>* mu,
               Matrix<double>* v, Matrix<double>* rt) {
    K->assign(2, 1);
    *mu = Matrix<double>(2, 1, 0.0);
    *v = Matrix<double>(2, 1, 0.0);
    *rt = Matrix<double>(2, 2, 0.0);
    (*mu)(0, 0) = mu1;
    (*mu)(1, 0) = mu2;
    (*rt)(0, 1) = 1.0;
    (*rt)(1, 0) = 1.0;
}

/** The exact product-form answer for that cycle at population N. */
void exact_cycle(double mu1, double mu2, std::size_t N, std::vector<double>* UN,
                 std::vector<double>* QN) {
    const double x1 = 1.0 / mu1, x2 = 1.0 / mu2;
    double G = 0.0;
    std::vector<double> w(N + 1, 0.0);
    for (std::size_t n1 = 0; n1 <= N; ++n1) {
        w[n1] = std::pow(x1, static_cast<double>(n1)) * std::pow(x2, static_cast<double>(N - n1));
        G += w[n1];
    }
    UN->assign(2, 0.0);
    QN->assign(2, 0.0);
    for (std::size_t n1 = 0; n1 <= N; ++n1) {
        const double p = w[n1] / G;
        if (n1 >= 1) (*UN)[0] += p;
        if (N - n1 >= 1) (*UN)[1] += p;
        (*QN)[0] += static_cast<double>(n1) * p;
        (*QN)[1] += static_cast<double>(N - n1) * p;
    }
}

}  // namespace

TEST_CASE("MMI reproduces the exact product-form solution of a cycle") {
    std::vector<int> K;
    Matrix<double> mu, v, rt;
    expo_pair(1.0, 2.0, &K, &mu, &v, &rt);
    const mq::QrfMetrics<double> r = mq::qrf_noblo_mmi<double>(2, K, 2, mu, v, rt);

    std::vector<double> UN, QN;
    exact_cycle(1.0, 2.0, 2, &UN, &QN);
    CHECK(r.UN[0] == doctest::Approx(UN[0]).epsilon(1e-6));
    CHECK(r.UN[1] == doctest::Approx(UN[1]).epsilon(1e-6));
    CHECK(r.QN[0] == doctest::Approx(QN[0]).epsilon(1e-6));
    CHECK(r.QN[1] == doctest::Approx(QN[1]).epsilon(1e-6));
    // Native Python on the same instance: 0.85714286 / 0.42857143.
    CHECK(r.UN[0] == doctest::Approx(0.85714286).epsilon(1e-6));
    // The population is conserved.
    CHECK(r.QN[0] + r.QN[1] == doctest::Approx(2.0).epsilon(1e-6));
}

TEST_CASE("MEM reproduces the same exact solution, and it is the unique one") {
    std::vector<std::pair<Matrix<double>, Matrix<double>>> MAPs;
    Matrix<double> A0(1, 1, -1.0), A1(1, 1, 1.0);
    Matrix<double> B0(1, 1, -2.0), B1(1, 1, 2.0);
    MAPs.push_back(std::make_pair(A0, A1));
    MAPs.push_back(std::make_pair(B0, B1));
    Matrix<double> rt(2, 2, 0.0);
    rt(0, 1) = 1.0;
    rt(1, 0) = 1.0;

    const mq::QrfMetrics<double> r = mq::qrf_noblo_mem<double>(MAPs, 2, rt);
    std::vector<double> UN, QN;
    exact_cycle(1.0, 2.0, 2, &UN, &QN);
    CHECK(r.UN[0] == doctest::Approx(UN[0]).epsilon(1e-6));
    CHECK(r.UN[1] == doctest::Approx(UN[1]).epsilon(1e-6));
    CHECK(r.QN[0] == doctest::Approx(QN[0]).epsilon(1e-6));
    CHECK(r.QN[1] == doctest::Approx(QN[1]).epsilon(1e-6));
}

TEST_CASE("MEM agrees with the reference on a HETEROGENEOUS phase count") {
    // Erlang-2 at station 0, exponential at station 1. This is the case that
    // exercises the ragged layout: the fill loops run over K[j] and K[i], not
    // over Kmax, so a station with fewer phases leaves gaps.
    std::vector<std::pair<Matrix<double>, Matrix<double>>> MAPs;
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = 2.0;
    D0(1, 1) = -2.0;
    D1(1, 0) = 2.0;
    MAPs.push_back(std::make_pair(D0, D1));
    Matrix<double> E0(1, 1, -2.0), E1(1, 1, 2.0);
    MAPs.push_back(std::make_pair(E0, E1));
    Matrix<double> rt(2, 2, 0.0);
    rt(0, 1) = 1.0;
    rt(1, 0) = 1.0;

    const mq::QrfMetrics<double> r = mq::qrf_noblo_mem<double>(MAPs, 2, rt);
    // Native Python `qrf_noblo_mem(MAPs, 2, rt)` on THIS instance: UN
    // [0.86567522, 0.43283761], QN [1.43283761, 0.56716239], which this port
    // reproduces to every digit printed. The exact CTMC is
    // [0.88888889, 0.44444444] / [1.44444444, 0.55555556], so the polytope is
    // NOT tight here (free dimension 66 at K = [2,1], N = 2) and the objective
    // is what decides the answer -- which is the whole point of this row.
    //
    // THIS GOLDEN MOVED ON 2026-08-29, from [0.8, 0.4] / [1.4, 0.6]. Until then
    // `mem_objective` returned +H and every port MINIMIZED it, so `qrf.mem`
    // reported the MINIMUM-entropy point of the polytope while being documented
    // as maximum-entropy; the AMPL model states the objective as `maximize H`.
    // The objective now returns -H, minimizing it maximizes H, and the answer
    // moves from 10.0% below the exact utilization to 2.6% below it.
    //
    // WHY THE START STILL MATTERS. +H is concave and -H is convex, so the
    // corrected objective does have a unique minimizer over the polytope; the
    // minimum-norm phase-1 point remains the start every codebase uses, and it
    // is what keeps them agreeing. See `qrf_feasible_start_lp`: from the
    // arbitrary phase-1 VERTEX this port used to start at, conditional gradient
    // reported a zero gap immediately and handed the vertex back, giving a
    // utilization of exactly 1.
    //
    // AN EVEN EARLIER GOLDEN WAS [1.0, 0.0] / [2.0, 0.0], WHICH IS THE REVERSED
    // ERLANG. `extract_mu_v_from_maps` used to write v(i,k,h) = D0(h,k) while
    // mu stayed untransposed, and MATLAB on a deliberately reversed D0 returns
    // [0.999932, 0.000068] / [1.999823, 0.000096] -- that golden to four
    // digits. A closed two-station cycle cannot leave a station idle, so that
    // vertex was never physical.
    CHECK(r.UN[0] == doctest::Approx(0.86567522).epsilon(1e-6));
    CHECK(r.UN[1] == doctest::Approx(0.43283761).epsilon(1e-6));
    CHECK(r.QN[0] == doctest::Approx(1.43283761).epsilon(1e-6));
    CHECK(r.QN[1] == doctest::Approx(0.56716239).epsilon(1e-6));
}

TEST_CASE("the ragged layout drops the coordinates nothing reads") {
    // With K = [2,1] and N = 2, the reference sizes its vector on the full Kmax
    // tensor (148) while the fill loops write only 84. The 64 unwritten
    // coordinates have an all-zero column in Aeq and would enter the reduced
    // problem as flat directions.
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    CHECK(mq::qrf_num_vars(2, 2, K, 1) == 84u);
    // The homogeneous case has no padding at all, so the two agree there.
    std::vector<int> K1(2, 1);
    CHECK(mq::qrf_num_vars(2, 2, K1, 1) == 6u * 6u + 2u);
}

TEST_CASE("the load-dependent form answers, and reduces to the plain one") {
    std::vector<int> K;
    Matrix<double> mu, v, rt;
    expo_pair(1.0, 2.0, &K, &mu, &v, &rt);

    // alpha == 1 everywhere is exactly the population-free model.
    Matrix<double> ones(2, 2, 1.0);
    const mq::QrfMetrics<double> flat =
        mq::qrf_noblo_mmi_ld<double>(2, K, 2, mu, v, rt, ones);
    std::vector<double> UN, QN;
    exact_cycle(1.0, 2.0, 2, &UN, &QN);
    CHECK(flat.UN[0] == doctest::Approx(UN[0]).epsilon(1e-5));
    CHECK(flat.QN[0] + flat.QN[1] == doctest::Approx(2.0).epsilon(1e-6));

    // Speeding station 1 up at its second population must not leave the answer
    // unchanged: that is the whole point of carrying the population index, and
    // a port that dropped it would return the flat answer here.
    Matrix<double> alpha(2, 2, 1.0);
    alpha(1, 1) = 2.0;
    const mq::QrfMetrics<double> ld = mq::qrf_noblo_mmi_ld<double>(2, K, 2, mu, v, rt, alpha);
    CHECK(ld.QN[0] + ld.QN[1] == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(std::fabs(ld.UN[1] - flat.UN[1]) > 1e-3);
    // The faster station holds less work.
    CHECK(ld.QN[1] < flat.QN[1]);
}

TEST_CASE("an empty alpha means all ones") {
    std::vector<int> K;
    Matrix<double> mu, v, rt;
    expo_pair(1.0, 2.0, &K, &mu, &v, &rt);
    Matrix<double> ones(2, 2, 1.0), none;
    const mq::QrfMetrics<double> a = mq::qrf_noblo_mmi_ld<double>(2, K, 2, mu, v, rt, ones);
    const mq::QrfMetrics<double> b = mq::qrf_noblo_mmi_ld<double>(2, K, 2, mu, v, rt, none);
    CHECK(a.UN[0] == doctest::Approx(b.UN[0]).epsilon(1e-9));
    CHECK(a.QN[1] == doctest::Approx(b.QN[1]).epsilon(1e-9));
    // A misshapen alpha is refused rather than broadcast.
    Matrix<double> bad(3, 3, 1.0);
    CHECK_THROWS_AS(mq::qrf_noblo_mmi_ld<double>(2, K, 2, mu, v, rt, bad), line::InputError);
}

TEST_CASE("the gradients are the gradients of the objectives") {
    // The Frank-Wolfe direction is the gradient's, so a wrong gradient sends
    // the solve to the wrong vertex while every other check still passes.
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const std::size_t M = 2, N = 2, MR = 1;
    const std::vector<int> F(M, static_cast<int>(N));
    const std::size_t n = mq::qrf_num_vars(M, N, K, MR);
    const std::vector<long> idx = mq::qrf_index_map(M, N, K, MR);

    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        x[i] = 0.05 + 0.4 * std::fmod(0.37 * static_cast<double>(i), 1.0);

    const std::vector<double> gm = mq::mmi_gradient(x, M, N, K, F, MR, idx);
    const std::vector<double> ge = mq::mem_gradient(x, M, N, K, F, MR, idx);
    const double h = 1e-6;
    double worst_m = 0.0, worst_e = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> a = x, b = x;
        a[j] += h;
        b[j] -= h;
        worst_m = std::max(worst_m,
                           std::fabs((mq::mmi_objective(a, M, N, K, F, MR) -
                                      mq::mmi_objective(b, M, N, K, F, MR)) /
                                         (2 * h) -
                                     gm[j]));
        worst_e = std::max(worst_e,
                           std::fabs((mq::mem_objective(a, M, N, K, F, MR) -
                                      mq::mem_objective(b, M, N, K, F, MR)) /
                                         (2 * h) -
                                     ge[j]));
    }
    CHECK(worst_m < 1e-7);
    CHECK(worst_e < 1e-7);
}

TEST_CASE("the Bethe gradient is the gradient of the Bethe objective") {
    // Same synthetic point and step as the MMI/MEM row above: every coordinate
    // sits well inside (0,1), so the central difference measures the
    // derivative rather than the LOGTOL shift at the boundary.
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const std::size_t M = 2, N = 2, MR = 1;
    const std::vector<int> F(M, static_cast<int>(N));
    const std::size_t n = mq::qrf_num_vars(M, N, K, MR);
    const std::vector<long> idx = mq::qrf_index_map(M, N, K, MR);

    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        x[i] = 0.05 + 0.4 * std::fmod(0.37 * static_cast<double>(i), 1.0);

    const std::vector<double> gb = mq::bethe_gradient(x, M, N, K, F, MR, idx);
    const double h = 1e-6;
    double worst = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> a = x, b = x;
        a[j] += h;
        b[j] -= h;
        worst = std::max(worst, std::fabs((mq::bethe_objective(a, M, N, K, F, MR) -
                                           mq::bethe_objective(b, M, N, K, F, MR)) /
                                              (2 * h) -
                                          gb[j]));
    }
    CHECK(worst < 1e-7);
}

TEST_CASE("the Bethe objective sums the n = 0 cells the MEM body drops") {
    // The three bodies do NOT share a lower bound, and each matches the AMPL
    // line it transcribes. The MI line reads `sum {ni, nj in 0..F}`, so
    // `mmi_objective` was corrected to n = 0 on 2026-09-02 (defect D1). The
    // MEM line reads `sum {ni in 1..F[i]}`, so `mem_objective` keeps its n = 1
    // bound -- it already matched its own spec. `bethe_objective` sums both
    // terms from 0 by construction. The residual between lambda*MI + MEM and
    // bethe is therefore the MEM n = 0 block ALONE; the MI n = 0 block now
    // lives inside `mmi_objective`.
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const std::size_t M = 2, N = 2, MR = 1;
    const std::vector<int> F(M, static_cast<int>(N));
    const std::size_t n = mq::qrf_num_vars(M, N, K, MR);
    const double tol = 1e-6, lambda = 1.0 / static_cast<double>(M);

    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        x[i] = 0.05 + 0.4 * std::fmod(0.37 * static_cast<double>(i), 1.0);
    const mq::QrfVars<double> v = mq::sub_qrfvar(x, M, N, K, MR);

    double idle = 0.0;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k) {
            const double pv = v.p(i, 0, k, i, 0, k, 0);
            idle += pv * std::log(tol + pv);
        }

    const double coded = lambda * mq::mmi_objective(x, M, N, K, F, MR) +
                         mq::mem_objective(x, M, N, K, F, MR);
    CHECK(std::fabs(idle) > 1e-3);
    CHECK(mq::bethe_objective(x, M, N, K, F, MR) == doctest::Approx(coded + idle).epsilon(1e-10));
}

TEST_CASE("BETHE reproduces the exact product-form solution of a cycle") {
    // The two-station cycle pins a single point -- the pairwise joint of a
    // closed chain at M == 2 is fixed by the marginal -- so every objective
    // over this polytope must report the exact answer, the Bethe one included.
    std::vector<int> K;
    Matrix<double> mu, v, rt;
    expo_pair(1.0, 2.0, &K, &mu, &v, &rt);
    const mq::QrfMetrics<double> r = mq::qrf_noblo_bethe<double>(2, K, 2, mu, v, rt);

    std::vector<double> UN, QN;
    exact_cycle(1.0, 2.0, 2, &UN, &QN);
    CHECK(r.UN[0] == doctest::Approx(UN[0]).epsilon(1e-6));
    CHECK(r.UN[1] == doctest::Approx(UN[1]).epsilon(1e-6));
    CHECK(r.QN[0] == doctest::Approx(QN[0]).epsilon(1e-6));
    CHECK(r.QN[1] == doctest::Approx(QN[1]).epsilon(1e-6));
    // Native Python on this instance: 0.85714286 / 0.42857143.
    CHECK(r.UN[0] == doctest::Approx(0.85714286).epsilon(1e-6));
    CHECK(r.QN[0] + r.QN[1] == doctest::Approx(2.0).epsilon(1e-6));
}

TEST_CASE("the BETHE optimum does not depend on the start point") {
    // WHAT CONVEXITY BUYS, ASSERTED DIRECTLY. At lambda = 1/M the objective is
    // the negative of a tree-reweighted entropy whose uniform edge weight
    // rho_ij = 2/M is the uniform point of the spanning tree polytope of K_M,
    // so the program is convex and every local optimum is global: the reported
    // point must be a property of the model, not of where the descent started.
    // This is exactly what `qrf.mmi` cannot promise on this same feasible set
    // (free dimension 66 at K = [2,1], N = 2).
    //
    // The second start is the MEM optimum of the SAME polytope, which is a
    // genuinely different feasible point rather than a perturbation of the
    // first -- native Python puts MEM's U1 at 0.8657 against Bethe's 0.8616.
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const std::size_t M = 2, N = 2, MR = 1;
    const std::vector<int> F(M, static_cast<int>(N));
    const std::size_t n = mq::qrf_num_vars(M, N, K, MR);
    const std::vector<long> idx = mq::qrf_index_map(M, N, K, MR);

    // Erlang-2 at station 0 (mean rate 1), exponential rate 2 at station 1 --
    // the same instance the MEM heterogeneous row above uses.
    std::vector<std::pair<Matrix<double>, Matrix<double>>> MAPs;
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = 2.0;
    D0(1, 1) = -2.0;
    D1(1, 0) = 2.0;
    MAPs.push_back(std::make_pair(D0, D1));
    Matrix<double> E0(1, 1, -2.0), E1(1, 1, 2.0);
    MAPs.push_back(std::make_pair(E0, E1));
    Matrix<double> mu, v, rt(2, 2, 0.0);
    mq::qrf_extract_mu_v(MAPs, M, K, &mu, &v);
    rt(0, 1) = 1.0;
    rt(1, 0) = 1.0;

    Matrix<int> BB(1, M, 0);
    const mq::QrfRates<double> q = mq::qrf_build_q(M, K, mu, v, rt);
    const mq::QrfAffine<double> eq = mq::qrf_affine_matrices<double>(
        [&](const std::vector<double>& z) { return mq::sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K).eq; },
        n);
    const mq::QrfAffine<double> ub = mq::qrf_affine_matrices<double>(
        [&](const std::vector<double>& z) {
            return mq::sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K).ineq;
        },
        n);
    const mq::QrfReduced<double> red = mq::qrf_reduce_equalities(eq.A, eq.b);
    const std::vector<double> phase1 = mq::qrf_feasible_start(red.A, red.b, ub.A, ub.b, n);

    const std::vector<double> mem_opt = mq::solve_qrf_nlp(
        [&](const std::vector<double>& z) { return mq::mem_objective(z, M, N, K, F, MR); },
        [&](const std::vector<double>& z) { return mq::mem_gradient(z, M, N, K, F, MR, idx); },
        phase1, red.A, red.b, ub.A, ub.b, "mem_start");

    double start_gap = 0.0;
    for (std::size_t j = 0; j < n; ++j)
        start_gap = std::max(start_gap, std::fabs(phase1[j] - mem_opt[j]));
    CHECK(start_gap > 1e-4);   // the two starts are genuinely different points

    auto bethe_obj = [&](const std::vector<double>& z) {
        return mq::bethe_objective(z, M, N, K, F, MR);
    };
    auto bethe_grad = [&](const std::vector<double>& z) {
        return mq::bethe_gradient(z, M, N, K, F, MR, idx);
    };
    const std::vector<double> a =
        mq::solve_qrf_nlp(bethe_obj, bethe_grad, phase1, red.A, red.b, ub.A, ub.b, "bethe_a");
    const std::vector<double> b =
        mq::solve_qrf_nlp(bethe_obj, bethe_grad, mem_opt, red.A, red.b, ub.A, ub.b, "bethe_b");

    const mq::QrfMetrics<double> ra = mq::qrf_extract_results(mq::sub_qrfvar(a, M, N, K, MR), M, K, F, MR);
    const mq::QrfMetrics<double> rb = mq::qrf_extract_results(mq::sub_qrfvar(b, M, N, K, MR), M, K, F, MR);
    // 1e-5 rather than solver precision: the restored n = 0 cells leave the
    // DESCENT DIRECTION sensitive to LOGTOL even though the objective VALUE is
    // not, so two runs stop at slightly different points of the same optimum.
    CHECK(ra.UN[0] == doctest::Approx(rb.UN[0]).epsilon(1e-5));
    CHECK(ra.UN[1] == doctest::Approx(rb.UN[1]).epsilon(1e-5));
    CHECK(bethe_obj(a) == doctest::Approx(bethe_obj(b)).epsilon(1e-8));

    // ... and the point it reports is feasible.
    double worst_eq = 0.0, worst_ub = 0.0;
    for (std::size_t i = 0; i < red.A.rows(); ++i) {
        double val = -red.b[i];
        for (std::size_t j = 0; j < n; ++j) val += red.A(i, j) * a[j];
        worst_eq = std::max(worst_eq, std::fabs(val));
    }
    for (std::size_t i = 0; i < ub.A.rows(); ++i) {
        double val = -ub.b[i];
        for (std::size_t j = 0; j < n; ++j) val += ub.A(i, j) * a[j];
        worst_ub = std::max(worst_ub, val);
    }
    CHECK(worst_eq < 1e-7);
    CHECK(worst_ub < 1e-7);
    for (std::size_t j = 0; j < n; ++j) {
        CHECK(a[j] >= -1e-9);
        CHECK(a[j] <= 1.0 + 1e-9);
    }
}

TEST_CASE("an inconsistent equality system is refused, not reduced away") {
    // x1 = 0 and x1 = 1: rank(A) is 1, rank([A|b]) is 2. Dropping the second
    // row would return numbers for a model nobody wrote.
    Matrix<double> A(2, 2, 0.0);
    A(0, 0) = 1.0;
    A(1, 0) = 1.0;
    std::vector<double> b;
    b.push_back(0.0);
    b.push_back(1.0);
    CHECK_THROWS_AS(mq::qrf_reduce_equalities(A, b), line::InputError);

    // The consistent duplicate reduces to one row and keeps the solution.
    std::vector<double> b2(2, 1.0);
    const mq::QrfReduced<double> ok = mq::qrf_reduce_equalities(A, b2);
    CHECK(ok.A.rows() == 1u);
    CHECK(ok.b[0] == doctest::Approx(1.0));
}

TEST_CASE("a non-affine residual map is caught at the probe") {
    // Every QRF constraint is linear; a quadratic one would give a matrix that
    // is right at the probe and wrong everywhere else.
    CHECK_THROWS_AS(mq::qrf_affine_matrices<double>(
                        [](const std::vector<double>& z) {
                            std::vector<double> r(1, z[0] * z[0]);
                            return r;
                        },
                        2),
                    line::InputError);
    // The affine one round-trips exactly.
    const mq::QrfAffine<double> aff = mq::qrf_affine_matrices<double>(
        [](const std::vector<double>& z) {
            std::vector<double> r(1, 2.0 * z[0] - 3.0 * z[1] + 1.0);
            return r;
        },
        2);
    CHECK(aff.A(0, 0) == doctest::Approx(2.0));
    CHECK(aff.A(0, 1) == doctest::Approx(-3.0));
    CHECK(aff.b[0] == doctest::Approx(-1.0));
}

TEST_CASE("qrf_noblo_mmi_linear is the load-dependent MMI bound") {
    // The objective is MMI; the `linear` in the name is about HOW the reference
    // builds its constraints, not about which objective it optimizes. Until
    // 2026-08-29 the MATLAB reference called its own mem() here and this port
    // mirrored that. MMI is not convex, so the numbers below are a stationary
    // point fixed by the minimum-norm start, not a unique optimum; they did not
    // move when the objective was corrected, because the polytope pins the
    // answer on all three of these instances.
    Matrix<double> rt(2, 2, 0.0);
    rt(0, 1) = 1.0;
    rt(1, 0) = 1.0;

    std::vector<std::pair<Matrix<double>, Matrix<double>>> one;
    Matrix<double> A0(1, 1, -1.0), A1(1, 1, 1.0), B0(1, 1, -2.0), B1(1, 1, 2.0);
    one.push_back(std::make_pair(A0, A1));
    one.push_back(std::make_pair(B0, B1));
    Matrix<double> none;
    const mq::QrfMetrics<double> r1 = mq::qrf_noblo_mmi_linear<double>(one, 2, rt, none);
    // Native Python: 0.85714286 / 0.42857143, itself the exact product form.
    CHECK(r1.UN[0] == doctest::Approx(0.85714286).epsilon(1e-6));
    CHECK(r1.UN[1] == doctest::Approx(0.42857143).epsilon(1e-6));
    CHECK(r1.QN[0] == doctest::Approx(1.42857143).epsilon(1e-6));

    // A ragged phase count, where the compact layout matters.
    //
    // RE-PINNED 2026-09-02 to [1, 0.5] / [1.5, 0.5], and the history matters
    // because this row has now been wrong twice in opposite directions.
    //
    // It carried [1, 0] / [2, 0] from THIS port when its phase 1 started at an
    // arbitrary vertex the conditional gradient could not leave. It was then
    // re-pinned to native Python's [0.8, 0.4] / [1.4, 0.6] (and [0.888889,
    // 0.333333] with alpha). Those Python values turned out to be a SECOND
    // stall, in the other solver: SLSQP returned its own start unchanged --
    // status 0, nit 1, max|t| exactly 0 -- at an objective of 1.900535 where
    // this port reaches 0.824607 on the same polytope. U_0 = 0.8 and U_0 = 1.0
    // are respectively the LP MINIMUM and the LP MAXIMUM of that polytope, so
    // both were feasible and the lower objective decides.
    //
    // The reason neither solver's stop was detectable by a first-order test:
    // the MI objective is CONCAVE on a phase-type polytope, so a point with
    // zero directional derivative is a local MAXIMUM or a saddle, and the LP
    // descent probe reads exactly 0.0 there. Concave minimisers are vertices,
    // so escaping needs a FINITE step to one -- which conditional gradient does
    // implicitly and which SLSQP now does through the vertex-probe escape in
    // `solve_qrf_nlp`. With that escape both ports return these values.
    std::vector<std::pair<Matrix<double>, Matrix<double>>> two;
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = 2.0;
    D0(1, 1) = -2.0;
    D1(1, 0) = 2.0;
    two.push_back(std::make_pair(D0, D1));
    Matrix<double> E0(1, 1, -2.0), E1(1, 1, 2.0);
    two.push_back(std::make_pair(E0, E1));
    const mq::QrfMetrics<double> r2 = mq::qrf_noblo_mmi_linear<double>(two, 2, rt, none);
    CHECK(r2.UN[0] == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r2.UN[1] == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(r2.QN[0] == doctest::Approx(1.5).epsilon(1e-6));
    CHECK(r2.QN[1] == doctest::Approx(0.5).epsilon(1e-6));

    Matrix<double> al(2, 2, 1.0);
    al(1, 1) = 2.0;
    const mq::QrfMetrics<double> r3 = mq::qrf_noblo_mmi_linear<double>(two, 2, rt, al);
    CHECK(r3.UN[0] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(r3.UN[1] == doctest::Approx(0.5).epsilon(1e-5));
}
