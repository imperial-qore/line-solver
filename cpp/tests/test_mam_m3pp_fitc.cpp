/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The M3PP counting-process fitters: the covariance-matching M3PP(2, 2), the
 * superposition and interleaving composites, and the trace and theoretical
 * entry points.
 *
 * WHAT IS ASSERTED, AND WHY IT IS EXACT. These are fitters, so the fitted
 * process generally does NOT equal the process the characteristics came from,
 * and nothing here compares generators against a reference implementation's.
 * What is exact, and is what these cases pin, is the algebra the fit is built
 * on:
 *
 *  1. THE COVARIANCE SPLIT IS AN EXACT INVERSION when the request is feasible.
 *     Given the underlying MMPP(2) of a genuine M3PP(2, 2), the reference
 *     inverts a QUADRATIC to recover the per-phase marking probabilities from
 *     the count covariance. Handed that process's own covariance the inversion
 *     must return its own (q1, q2), not an approximation of them, and the
 *     covariance recomputed from the result must equal the request. Both are
 *     asserted to 1e-10.
 *  2. A MARKING IS A PARTITION. sum_c Dc = D1 holds entry by entry for every
 *     composite here, because each class matrix is a diagonal scaling of D1 by
 *     probabilities that sum to one. Asserted exactly on the interleaved
 *     process, where the class matrices are assembled level by level.
 *  3. THE PER-CLASS RATES ARE MATCHED EXACTLY, by construction, in every
 *     fitter of this family: the rate enters the closed forms as a linear
 *     constraint that is satisfied identically, never as a least-squares term.
 *     This is the one property that survives every method and every composite,
 *     so it is checked on all of them.
 *  4. THE ORDERS ARE STRUCTURAL. Superposing L two-phase processes gives the
 *     PRODUCT chain, order 2^L; interleaving lumps them onto a birth-death
 *     chain of order L + 1. Those are properties of the construction, not of
 *     the data.
 *  5. mmpp2_fitc_theoretical ROUND TRIPS an MMPP(2): fed the exact counting
 *     characteristics of a two-phase MMPP, the fit reproduces it up to the
 *     phase relabelling, since the characteristics determine the four
 *     parameters.
 *
 * The Poisson regime is deliberately absent: a Poisson trace has IDC = 1 at
 * every scale, the fitters then return a REDUCIBLE MMPP(2), and MATLAB's
 * ctmc_solve refuses it with the same message the native Python one does. That
 * regime is out of the family's domain by construction, not a case to pin.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mam/m3pp22_fitc_cov.h"
#include "line/api/mam/m3pp2m_fitc_trace.h"
#include "line/api/mam/m3pp2m_interleave.h"
#include "line/api/mam/m3pp_superpos_fitc.h"
#include "line/api/mam/map_rand.h"
#include "line/api/npfqn/npfqn_traffic_merge.h"

using namespace line;
using namespace line::mam;

namespace {

/** A genuine M3PP(2, 2): an MMPP(2) with a per-phase Bernoulli marking. */
Mmap<double> genuine_m3pp22(double l1, double l2, double r1, double r2, double q1, double q2) {
    Mmap<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 1) = r1;
    m.D0(1, 0) = r2;
    m.D1(0, 0) = l1;
    m.D1(1, 1) = l2;
    m.D0(0, 0) = -(l1 + r1);
    m.D0(1, 1) = -(l2 + r2);
    Matrix<double> A(2, 2, 0.0), B(2, 2, 0.0);
    A(0, 0) = q1 * l1;
    A(1, 1) = q2 * l2;
    B(0, 0) = (1.0 - q1) * l1;
    B(1, 1) = (1.0 - q2) * l2;
    m.Dc.push_back(A);
    m.Dc.push_back(B);
    return m;
}

/** Worst entrywise deviation of sum_c Dc from D1. */
double marking_defect(const Mmap<double>& m) {
    double worst = 0.0;
    for (std::size_t i = 0; i < m.order(); ++i)
        for (std::size_t j = 0; j < m.order(); ++j) {
            double s = 0.0;
            for (std::size_t c = 0; c < m.classes(); ++c) s += m.Dc[c](i, j);
            const double d = s - m.D1(i, j);
            worst = std::max(worst, std::fabs(d));
        }
    return worst;
}

}  // namespace

TEST_CASE("m3pp22_fitc_approx_cov_multiclass inverts its own covariance exactly") {
    const double l1 = 2.0, l2 = 0.5, r1 = 0.3, r2 = 0.7, q1 = 0.7, q2 = 0.3, t3 = 1.0;
    const Mmap<double> truth = genuine_m3pp22(l1, l2, r1, r2, q1, q2);
    const std::vector<double> ai = mmap_count_mean(truth, 1.0);
    const double st3 = mmap_count_mcov(truth, t3)(0, 1);

    const M3pp22FitcCovResult<double> f =
        m3pp22_fitc_approx_cov_multiclass(truth.map(), ai, st3, t3);

    CHECK_FALSE(f.clamped);
    CHECK(f.mmap.Dc[0](0, 0) / l1 == doctest::Approx(q1).epsilon(1e-10));
    CHECK(f.mmap.Dc[0](1, 1) / l2 == doctest::Approx(q2).epsilon(1e-10));
    CHECK(mmap_count_mcov(f.mmap, t3)(0, 1) == doctest::Approx(st3).epsilon(1e-10));
    CHECK(marking_defect(f.mmap) < 1e-12);

    const std::vector<double> got = mmap_count_mean(f.mmap, 1.0);
    CHECK(got[0] == doctest::Approx(ai[0]).epsilon(1e-12));
    CHECK(got[1] == doctest::Approx(ai[1]).epsilon(1e-12));
}

TEST_CASE("m3pp22_fitc_approx_cov_multiclass short-circuits the degenerate cases") {
    Map<double> poisson;
    poisson.D0 = Matrix<double>(1, 1, -3.0);
    poisson.D1 = Matrix<double>(1, 1, 3.0);
    std::vector<double> ai;
    ai.push_back(1.0);
    ai.push_back(2.0);
    const M3pp22FitcCovResult<double> f =
        m3pp22_fitc_approx_cov_multiclass(poisson, ai, 0.0, 1.0);
    CHECK(f.degenerate);
    CHECK(f.mmap.Dc[0](0, 0) == doctest::Approx(1.0));
    CHECK(f.mmap.Dc[1](0, 0) == doctest::Approx(2.0));

    std::vector<double> one;
    one.push_back(3.0);
    const Map<double> mmpp = genuine_m3pp22(2.0, 0.5, 0.3, 0.7, 1.0, 1.0).map();
    const M3pp22FitcCovResult<double> g = m3pp22_fitc_approx_cov_multiclass(
        mmpp, std::vector<double>(1, map_lambda(mmpp)), 0.0, 1.0);
    CHECK(g.degenerate);
    CHECK(g.mmap.classes() == 1u);
}

TEST_CASE("m3pp_superpos_fitc gives the product chain and matches every class rate") {
    std::vector<double> av, btv, binfv, m3tv;
    av.push_back(1.0);
    av.push_back(2.0);
    btv.push_back(2.0);
    btv.push_back(3.0);
    binfv.push_back(5.0);
    binfv.push_back(6.0);
    m3tv.push_back(0.5);
    m3tv.push_back(0.8);

    const M3ppSuperposResult<double> s = m3pp_superpos_fitc(av, btv, binfv, m3tv, 1.0, 100.0);
    CHECK(s.parts.size() == 2u);
    CHECK(s.mmap.classes() == 2u);
    CHECK(s.mmap.order() == 4u);  // 2 x 2, the product chain

    const std::vector<double> got = mmap_count_mean(s.mmap, 1.0);
    CHECK(got[0] == doctest::Approx(av[0]).epsilon(1e-10));
    CHECK(got[1] == doctest::Approx(av[1]).epsilon(1e-10));
    CHECK(marking_defect(s.mmap) < 1e-12);
}

TEST_CASE("m3pp_superpos_fitc_theoretical reproduces the per-class rates of its target") {
    const Mmap<double> truth = genuine_m3pp22(2.0, 0.5, 0.3, 0.7, 0.7, 0.3);
    const std::vector<double> want = mmap_count_mean(truth, 1.0);

    const M3ppSuperposResult<double> s = m3pp_superpos_fitc_theoretical(truth, 1.0, 100.0);
    const std::vector<double> got = mmap_count_mean(s.mmap, 1.0);
    CHECK(s.mmap.order() == 4u);
    CHECK(got[0] == doctest::Approx(want[0]).epsilon(1e-10));
    CHECK(got[1] == doctest::Approx(want[1]).epsilon(1e-10));
}

TEST_CASE("m3pp2m_interleave lumps L two-phase processes onto L + 1 levels") {
    std::vector<Mmap<double>> parts;
    parts.push_back(genuine_m3pp22(2.0, 0.5, 0.9, 0.4, 0.7, 0.3));
    parts.push_back(genuine_m3pp22(3.0, 0.8, 1.4, 0.9, 0.4, 0.6));

    const Mmap<double> s = m3pp2m_interleave(parts);
    CHECK(s.order() == 3u);     // L + 1
    CHECK(s.classes() == 4u);   // the two class lists, concatenated
    CHECK(marking_defect(s) < 1e-12);

    // the generator rows sum to zero, so it is a proper MAP
    for (std::size_t i = 0; i < s.order(); ++i) {
        double row = 0.0;
        for (std::size_t j = 0; j < s.order(); ++j) row += s.D0(i, j) + s.D1(i, j);
        CHECK(std::fabs(row) < 1e-12);
    }
}

TEST_CASE("m3pp22_interleave_fitc matches every per-class rate exactly") {
    Matrix<double> av(2, 2, 0.0);
    av(0, 0) = 0.6;
    av(0, 1) = 0.4;
    av(1, 0) = 1.2;
    av(1, 1) = 0.8;
    std::vector<double> btv, binfv, stv;
    btv.push_back(2.0);
    btv.push_back(2.5);
    binfv.push_back(6.0);
    binfv.push_back(7.0);
    stv.push_back(1.0);
    stv.push_back(1.0);

    const M3pp22InterleaveResult<double> r = m3pp22_interleave_fitc(av, btv, binfv, stv, 1.0);
    CHECK(r.parts.size() == 2u);
    CHECK(r.mmap.order() == 3u);
    CHECK(r.mmap.classes() == 4u);
    CHECK(marking_defect(r.mmap) < 1e-12);

    const std::vector<double> got = mmap_count_mean(r.mmap, 1.0);
    CHECK(got[0] == doctest::Approx(0.6).epsilon(1e-10));
    CHECK(got[1] == doctest::Approx(0.4).epsilon(1e-10));
    CHECK(got[2] == doctest::Approx(1.2).epsilon(1e-10));
    CHECK(got[3] == doctest::Approx(0.8).epsilon(1e-10));
}

TEST_CASE("m3pp22_interleave_fitc refuses an infeasible IDC pair by name") {
    Matrix<double> av(1, 2, 0.0);
    av(0, 0) = 0.5;
    av(0, 1) = 0.5;
    std::vector<double> btv, binfv, stv;
    btv.push_back(0.5);  // IDC(t) below one: sub-Poisson, outside the family
    binfv.push_back(3.0);
    stv.push_back(0.0);
    CHECK_THROWS_AS(m3pp22_interleave_fitc(av, btv, binfv, stv, 1.0), InputError);
}

TEST_CASE("m3pp2m_fitc_theoretical matches the per-class rates on every method") {
    const Mmap<double> truth = genuine_m3pp22(2.0, 0.5, 0.3, 0.7, 0.7, 0.3);
    const std::vector<double> want = mmap_count_mean(truth, 1.0);

    const char* methods[] = {"exact_delta", "approx_delta", "approx_cov", "approx_ag"};
    for (int k = 0; k < 4; ++k) {
        const Mmap<double> f = m3pp2m_fitc_theoretical(truth, std::string(methods[k]), 1.0, 1e4);
        CAPTURE(methods[k]);
        CHECK(f.order() == 2u);
        CHECK(f.classes() == 2u);
        const std::vector<double> got = mmap_count_mean(f, 1.0);
        CHECK(got[0] == doctest::Approx(want[0]).epsilon(1e-8));
        CHECK(got[1] == doctest::Approx(want[1]).epsilon(1e-8));
    }
    CHECK_THROWS_AS(m3pp2m_fitc_theoretical(truth, std::string("nonesuch"), 1.0, 1e4), InputError);
}

TEST_CASE("mmpp2_fitc_theoretical round trips an MMPP(2)") {
    const Map<double> truth = genuine_m3pp22(2.0, 0.5, 0.3, 0.7, 1.0, 1.0).map();
    const Mmpp2FitcResult<double> f = mmpp2_fitc_theoretical(truth);
    CHECK_FALSE(f.degenerate);
    CHECK(f.map.order() == 2u);
    // the same four parameters, up to which phase is listed first
    const double a = map_lambda(f.map), want = map_lambda(truth);
    CHECK(a == doctest::Approx(want).epsilon(1e-8));
    const double hi = std::max(f.map.D1(0, 0), f.map.D1(1, 1));
    const double lo = std::min(f.map.D1(0, 0), f.map.D1(1, 1));
    CHECK(hi == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(lo == doctest::Approx(0.5).epsilon(1e-6));
}

TEST_CASE("m3pp2m_fitc_trace reads the rate and the class split off the trace") {
    // an MMPP(2)-modulated two-class trace: fast phase favours class 1
    std::vector<double> Tv;
    std::vector<int> A;
    unsigned long seed = 20260801UL;
    auto unif = [&seed]() {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        return double((seed >> 11) & ((1ULL << 53) - 1)) / double(1ULL << 53);
    };
    const double lam[2] = {4.0, 0.4}, sw[2] = {0.05, 0.05};
    int ph = 0;
    double total = 0.0;
    long n1 = 0;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        const double te = -std::log(1.0 - unif()) / lam[ph];
        const double ts = -std::log(1.0 - unif()) / sw[ph];
        double dt = te;
        if (ts < te) {
            ph = 1 - ph;
            dt = ts;
        }
        Tv.push_back(dt);
        total += dt;
        const int c = unif() < (ph == 0 ? 0.3 : 0.7) ? 1 : 2;
        A.push_back(c);
        if (c == 1) ++n1;
    }
    const double a = double(n) / total;
    const double p1 = double(n1) / double(n);

    const M3pp2mFitcTraceResult<double> f = m3pp2m_fitc_trace(Tv, A, std::string("approx_ag"));
    CHECK(f.a == doctest::Approx(a).epsilon(1e-12));
    CHECK(f.bt1 > 1.0);  // the trace is over-dispersed, which is what the family needs
    const std::vector<double> got = mmap_count_mean(f.mmap, 1.0);
    CHECK(got[0] == doctest::Approx(a * p1).epsilon(1e-6));
    CHECK(got[1] == doctest::Approx(a * (1.0 - p1)).epsilon(1e-6));

    const M3ppSuperposResult<double> s = m3pp_superpos_fitc_trace(Tv, A);
    const std::vector<double> sr = mmap_count_mean(s.mmap, 1.0);
    CHECK(sr[0] == doctest::Approx(a * p1).epsilon(1e-6));
    CHECK(sr[1] == doctest::Approx(a * (1.0 - p1)).epsilon(1e-6));

    CHECK_THROWS_AS(m3pp2m_fitc_trace(Tv, A, std::string("nonesuch")), InputError);
}

TEST_CASE("m3pp_rand marks a random MMPP without disturbing it") {
    std::mt19937 gen(7);
    const Mmap<double> m = m3pp_rand<double>(3, 4, gen);
    CHECK(m.order() == 3u);
    CHECK(m.classes() == 4u);
    CHECK(marking_defect(m) < 1e-12);
}

TEST_CASE("npfqn_traffic_merge serves the interpos merge") {
    std::vector<Mmap<double>> flows;
    flows.push_back(genuine_m3pp22(2.0, 0.5, 0.3, 0.7, 0.7, 0.3));
    flows.push_back(genuine_m3pp22(3.0, 0.8, 0.4, 0.6, 0.4, 0.6));
    npfqn::MergeConfig cfg;
    cfg.merge = npfqn::Merge::Interpos;
    cfg.compress = npfqn::Compress::None;

    const Mmap<double> s = npfqn::npfqn_traffic_merge(flows, cfg);
    CHECK(s.order() == 3u);    // two flows lumped onto three levels
    CHECK(s.classes() == 4u);  // the two class lists, concatenated
    CHECK(marking_defect(s) < 1e-10);
}
