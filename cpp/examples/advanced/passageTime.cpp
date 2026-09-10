/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * First passage times in a Markov chain, and the exact cycle time along an
 * overtake-free path of a closed tree-like product-form network.
 *
 * Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
 * in Large Markov Chains", 2002.
 *
 * Twin of matlab/examples/advanced/passageTime/passage_firstpassage.m and
 * python/examples/advanced/passageTime/passage_firstpassage.py.
 */

#include <cmath>
#include <cstdio>
#include <vector>

#include "examples_common.h"
#include "line/api/mc/ctmc_passage.h"
#include "line/api/pfqn/pfqn_cyclet_ofree.h"

namespace line {
namespace examples {

namespace {

double interp(const std::vector<double>& x, const std::vector<double>& y, double xq) {
    for (std::size_t i = 1; i < x.size(); ++i)
        if (x[i] >= xq) {
            const double w = (xq - x[i - 1]) / (x[i] - x[i - 1]);
            return y[i - 1] + w * (y[i] - y[i - 1]);
        }
    return y.empty() ? 0.0 : y.back();
}

}  // namespace

void passage_firstpassage() {
    using line::Matrix;
    namespace mc = line::mc;
    namespace pfqn = line::pfqn;

    // 1. The time for an M/M/1/K queue to fill from empty. This is a first
    // passage into a STATE SET, which no response-time getter can express: it
    // is the chain reaching a marking, not a job finishing service.
    const std::size_t K = 6, n = K + 1;
    const double lambda = 1.0, mu = 1.5;
    Matrix<double> Q(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i + 1 < n) Q(i, i + 1) = lambda;
        if (i > 0) Q(i, i - 1) = mu;
    }
    for (std::size_t i = 0; i < n; ++i) {
        double r = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            if (j != i) r += Q(i, j);
        Q(i, i) = -r;
    }

    std::vector<double> pi0(n, 0.0);
    pi0[0] = 1.0;                                   // start empty
    const std::vector<std::size_t> target{n - 1};   // the full buffer

    const mc::PassageMoments<double> pm = mc::ctmc_passage_moments(Q, pi0, target, 3);
    const double var = pm.m[1] - pm.m[0] * pm.m[0];
    std::printf("Time to fill an M/M/1/%zu from empty (lambda=%g, mu=%g)\n", K, lambda, mu);
    std::printf("  mean            = %.6f\n", pm.m[0]);
    std::printf("  variance        = %.6f\n", var);
    std::printf("  coeff. of var.  = %.6f\n", std::sqrt(var) / pm.m[0]);
    std::printf("  from state K-1  = %.6f (one arrival away)\n", pm.mall(n - 2, 0));

    std::vector<double> tset(400);
    for (std::size_t i = 0; i < tset.size(); ++i)
        tset[i] = 4.0 * pm.m[0] * double(i) / double(tset.size() - 1);
    const mc::PassageCurve<double> c = mc::ctmc_passage_time(Q, pi0, target, tset);
    std::printf("  P(fill <= mean) = %.6f\n", interp(tset, c.F, pm.m[0]));

    const std::vector<double> h = mc::ctmc_hitting_time(Q, target);
    std::printf("  hitting times   =");
    for (double v : h) std::printf(" %.3f", v);
    std::printf("\n");

    // 2. The same passage with a non-exponential sojourn (semi-Markov). The
    // embedded chain is unchanged; only the holding-time law moves. Nothing in
    // a generator can express this, which is why the semi-Markov route exists.
    Matrix<double> P(n, n);
    std::vector<double> rate(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        rate[i] = -Q(i, i);
        if (rate[i] > 0.0) {
            for (std::size_t j = 0; j < n; ++j) P(i, j) = (j == i) ? 0.0 : Q(i, j) / rate[i];
        } else {
            P(i, i) = 1.0;
        }
    }
    Matrix<double> hmom(n, 3);
    for (std::size_t i = 0; i < n; ++i) {
        const double d = 1.0 / rate[i];
        hmom(i, 0) = d;
        hmom(i, 1) = d * d;
        hmom(i, 2) = d * d * d;
    }
    const mc::PassageMoments<double> mD = mc::smp_passage_moments(P, hmom, pi0, target, 3);
    const double varD = mD.m[1] - mD.m[0] * mD.m[0];
    std::printf("\nSame embedded chain, DETERMINISTIC sojourns of equal mean:\n");
    std::printf("  mean            = %.6f (unchanged, as it must be)\n", mD.m[0]);
    std::printf("  coeff. of var.  = %.6f (was %.6f)\n", std::sqrt(varD) / mD.m[0],
                std::sqrt(var) / pm.m[0]);

    // 3. The cycle time of the tree-like network of Fig. 6 of the paper.
    const std::vector<double> mu6{3.0, 5.0, 4.0, 6.0, 2.0, 1.0};
    const double p12 = 0.2, p13 = 0.5, p14 = 0.3;
    const std::vector<double> v{1.0, p12, p13, p14, p12, p14};
    const std::vector<std::vector<std::size_t>> paths{{0, 2}, {0, 1, 4}, {0, 3, 5}};
    std::vector<double> tt;
    for (double t = 0.0; t <= 40.0 + 1e-9; t += 0.25) tt.push_back(t);
    const pfqn::CycletResult r =
        pfqn::pfqn_cyclet_ofree(v, mu6, 18, paths, tt, "auto", 3, {p13, p12, p14});
    std::printf("\nTree network of Fig. 6, N = 18 customers\n");
    std::printf("  moments  = %.5f  %.4f  %.3f\n", r.mom[0], r.mom[1], r.mom[2]);
    std::printf("  paper    = 6.12717  53.3067  612.887\n");
    std::printf("  routes   =");
    for (const pfqn::CycletPathInfo& i : r.info) std::printf(" %s", i.method.c_str());
    std::printf("\n  P(cycle <= mean) = %.6f\n", interp(tt, r.F, r.mom[0]));
}
LINE_EXAMPLE("advanced/passageTime", passage_firstpassage);

}  // namespace examples
}  // namespace line
