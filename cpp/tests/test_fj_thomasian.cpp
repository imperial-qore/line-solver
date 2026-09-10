/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The fork-join formulas of A. Thomasian, "Analysis of Fork/Join and Related
 * Queueing Systems", ACM Computing Surveys 47(2), Article 17, 2014.
 *
 * The oracles are independent of the port: the survey's own worked numbers (the
 * independent-server utilization 0.9286, the team-service capacities 1.6 and
 * 4/3, the two task-system makespans 26.77 and 10.1), and exact identities the
 * formulas must reproduce (unit fork degrees turn fj_qgb into the ordinary
 * geometric bound and fj_amva into exact mean value analysis; the moment
 * recurrence must agree with inclusion-exclusion; a Coxian with q = 0 is an
 * exponential, whose maximum is H_K; a unit batch on one server is an M/M/1).
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/fj/fj_amva.h"
#include "line/api/fj/fj_char_max_blom.h"
#include "line/api/fj/fj_char_max_discrete.h"
#include "line/api/fj/fj_cox_fit.h"
#include "line/api/fj/fj_dag_makespan.h"
#include "line/api/fj/fj_delay_opt.h"
#include "line/api/fj/fj_dispersion.h"
#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_ism_green.h"
#include "line/api/fj/fj_lst_max_het.h"
#include "line/api/fj/fj_qgb.h"
#include "line/api/fj/fj_respt_bulk.h"
#include "line/api/fj/fj_respt_closed.h"
#include "line/api/fj/fj_respt_nosplit.h"
#include "line/api/fj/fj_serialization.h"
#include "line/api/fj/fj_tsm_capacity.h"
#include "line/api/fj/fj_xmax_coxian.h"
#include "line/api/fj/fj_xmax_het.h"
#include "line/api/fj/fj_xmax_hz.h"
#include "line/api/fj/fj_xmax_hz_het.h"
#include "line/api/fj/fj_xmax_moments_het.h"

using namespace line;
using namespace line::fj;

TEST_CASE("fj_qgb with unit fork degrees is the ordinary geometric bound") {
    std::vector<double> D;
    D.push_back(1.0);
    D.push_back(2.0);
    D.push_back(3.0);
    const std::vector<unsigned> P(3, 1u);
    const FJQgbResult<double> g = fj_qgb<double>(D, P, 7u, 1.5);
    const double denom = 1.5 + 6.0 + 3.0 * 7.0;
    for (int i = 0; i < 3; ++i) {
        const double y = D[i] * 7.0 / denom;
        CHECK(g.Q[i] == doctest::Approx(y / (1 - y) - std::pow(y, 8.0) / (1 - y)).epsilon(1e-12));
    }
}

TEST_CASE("fj_amva with unit fork degrees is exact mean value analysis") {
    std::vector<double> D;
    D.push_back(1.0);
    D.push_back(2.0);
    D.push_back(3.0);
    const std::vector<unsigned> P(3, 1u);
    const FJAmvaResult<double> a = fj_amva<double>(D, P, 7u, 1.5);
    double Q[3] = {0, 0, 0}, X = 0;
    for (unsigned m = 1; m <= 7; ++m) {
        double R[3], tot = 0;
        for (int i = 0; i < 3; ++i) {
            R[i] = D[i] * (1 + Q[i]);
            tot += R[i];
        }
        X = m / (1.5 + tot);
        for (int i = 0; i < 3; ++i) Q[i] = X * R[i];
    }
    CHECK(a.X == doctest::Approx(X).epsilon(1e-12));
    for (int i = 0; i < 3; ++i) CHECK(a.Q[i] == doctest::Approx(Q[i]).epsilon(1e-12));
}

TEST_CASE("fj_respt_closed is tight at two branches") {
    const FJResptClosedResult<double> r = fj_respt_closed<double>(2u, 0.4, 5u);
    CHECK(r.R == doctest::Approx(0.4 * (1.5 + 4)).epsilon(1e-12));
    CHECK(r.exact);
}

TEST_CASE("fj_xmax_het matches the textbook answers") {
    std::vector<double> two;
    two.push_back(1.0);
    two.push_back(2.0);
    CHECK(fj_xmax_het<double>(two) == doctest::Approx(1 + 0.5 - 1.0 / 3.0).epsilon(1e-12));
    CHECK(fj_xmax_het<double>(std::vector<double>(4, 2.0)) ==
          doctest::Approx(fj_harmonic<double>(4) / 2.0).epsilon(1e-12));
}

TEST_CASE("the moment recurrence agrees with inclusion-exclusion") {
    std::vector<double> lam;
    lam.push_back(1.0);
    lam.push_back(2.0);
    lam.push_back(3.0);
    lam.push_back(5.0);
    const std::vector<double> m = fj_xmax_moments_het<double>(lam, 3u);
    for (unsigned n = 1; n <= 3; ++n)
        CHECK(m[n - 1] == doctest::Approx(fj_xmax_het<double>(lam, n)).epsilon(1e-10));
}

TEST_CASE("the transform is one at the origin and its slope is the mean") {
    std::vector<double> lam;
    lam.push_back(1.0);
    lam.push_back(2.0);
    lam.push_back(3.0);
    lam.push_back(5.0);
    CHECK(fj_lst_max_het<double>(lam, 0.0) == doctest::Approx(1.0).epsilon(1e-12));
    const double h = 1e-5;
    CHECK(-(fj_lst_max_het<double>(lam, h) - 1.0) / h ==
          doctest::Approx(fj_xmax_het<double>(lam)).epsilon(1e-4));
}

TEST_CASE("the Harrison-Zertal closed form is exact for the exponential") {
    const double m1 = 0.7;
    CHECK(fj_xmax_hz<double>(m1, 2 * m1 * m1, 6u).Xmax ==
          doctest::Approx(fj_harmonic<double>(6) * m1).epsilon(1e-12));
}

TEST_CASE("the Harrison-Zertal recurrence is exact for i.i.d. exponentials") {
    const unsigned K = 4;
    const double lam = 1.3;
    std::vector<double> m1(K, 1.0 / lam), m2(K, 2.0 / (lam * lam));
    std::vector<std::function<double(const double&)> > cdf;
    for (unsigned i = 0; i < K; ++i)
        cdf.push_back([lam](const double& t) { return t > 0 ? 1 - std::exp(-lam * t) : 0.0; });
    CHECK(fj_xmax_hz_het<double>(m1, m2, cdf) ==
          doctest::Approx(fj_harmonic<double>(K) / lam).epsilon(1e-6));
}

TEST_CASE("the characteristic maximum bounds the lattice maximum") {
    const FJCharMaxDiscreteResult<double> g =
        fj_char_max_discrete<double>(8u, FJDiscreteDist::Geometric, 0.6);
    CHECK(g.MK >= g.exact - 1e-9);
    const FJCharMaxDiscreteResult<double> p =
        fj_char_max_discrete<double>(8u, FJDiscreteDist::Poisson, 4.0);
    CHECK(p.MK >= p.exact - 1e-9);
}

TEST_CASE("the Blom position sits inside the Kruskal-Weiss bracket") {
    const FJCharMaxBlomResult<double> b = fj_char_max_blom<double>(20u);
    CHECK(b.bracket_available);
    CHECK(b.mK > b.lo);
    CHECK(b.mK < b.hi);
}

TEST_CASE("the Coxian fit reproduces its targets and degenerates to the exponential") {
    const FJCoxFitResult<double> f = fj_cox_fit<double>(2.0, 1.5);
    const FJXmaxCoxianResult<double> c1 = fj_xmax_coxian<double>(1u, f.mu1, f.mu2, f.q);
    CHECK(c1.m1 == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(c1.c2 == doctest::Approx(1.5).epsilon(1e-12));
    CHECK(c1.Xmax == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(fj_xmax_coxian<double>(5u, 1.0, 1.0, 0.0).Xmax ==
          doctest::Approx(fj_harmonic<double>(5)).epsilon(1e-10));
}

TEST_CASE("the dispersion of two exponential branches") {
    const double mu = 1.7;
    const std::vector<unsigned> shape(2, 1u);
    const std::vector<double> rate(2, mu);
    const FJDispersionResult<double> d = fj_dispersion<double>(shape, rate);
    CHECK(d.Emax == doctest::Approx(1.5 / mu).epsilon(1e-7));
    CHECK(d.Emin == doctest::Approx(0.5 / mu).epsilon(1e-7));
    CHECK(d.Edisp == doctest::Approx(1.0 / mu).epsilon(1e-7));
}

TEST_CASE("delaying never increases the dispersion") {
    std::vector<unsigned> shape;
    shape.push_back(1);
    shape.push_back(3);
    shape.push_back(2);
    std::vector<double> rate;
    rate.push_back(1.0);
    rate.push_back(2.0);
    rate.push_back(0.8);
    const double d0 = fj_dispersion<double>(shape, rate).Edisp;
    const FJDelayOptResult<double> o = fj_delay_opt<double>(shape, rate);
    CHECK(o.Edisp <= d0 + 1e-9);
    double dmin = o.d[0];
    for (std::size_t i = 1; i < o.d.size(); ++i)
        if (o.d[i] < dmin) dmin = o.d[i];
    CHECK(dmin == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("no splitting collapses to M/M/1 at one task") {
    CHECK(fj_respt_nosplit<double>(1u, 0.5, 1.4).R == doctest::Approx(1.0 / 0.9).epsilon(1e-12));
}

TEST_CASE("bulk arrivals collapse to M/M/1 at a unit batch on one server") {
    const FJResptBulkResult<double> b = fj_respt_bulk<double>(1u, 0.5, 1.4, 1u);
    const double rho = 0.5 / 1.4;
    CHECK(b.Q == doctest::Approx(rho / (1 - rho)).epsilon(1e-7));
    CHECK(b.Rreq == doctest::Approx(1.0 / 0.9).epsilon(1e-7));
}

TEST_CASE("the independent server model reproduces the survey's utilization") {
    std::vector<double> c;
    c.push_back(0.1);
    c.push_back(0.2);
    c.push_back(0.3);
    c.push_back(0.4);
    const FJIsmGreenResult<double> g = fj_ism_green<double>(1.0, 1.4, 4u, c);
    // The survey quotes rho = 0.9286 for this instance
    CHECK(g.rho == doctest::Approx(0.9285714285714286).epsilon(1e-9));
    CHECK(g.W > 0);
    CHECK(g.pq > 0);
    CHECK(g.pq < 1);
    CHECK(g.pd > 0);
    CHECK(g.pd <= 1);
}

TEST_CASE("the team service capacity matches the survey's worked example") {
    const std::vector<double> f(4, 0.25), x(4, 1.0);
    std::vector<unsigned> r;
    r.push_back(1);
    r.push_back(2);
    r.push_back(3);
    r.push_back(4);
    const FJTsmCapacityResult<double> t = fj_tsm_capacity<double>(4u, f, r, x);
    // Lambda_max = s / sum f r x = 4 / 2.5 = 1.6, attained because every state
    // carrying probability is full capacity
    CHECK(t.Lmax == doctest::Approx(1.6).epsilon(1e-12));
    CHECK(t.Llp == doctest::Approx(1.6).epsilon(1e-7));
    std::vector<double> fs;
    fs.push_back(1e-6);
    fs.push_back(0.5 - 1e-6);
    fs.push_back(0.5 - 1e-6);
    fs.push_back(1e-6);
    // The skewed frequency vector of the survey reaches only 4/3
    CHECK(fj_tsm_capacity<double>(4u, fs, r, x).Llp ==
          doctest::Approx(4.0 / 3.0).epsilon(1e-4));
    std::vector<double> f2(2, 0.5), x2;
    x2.push_back(0.5);
    x2.push_back(1.0 / 3.0);
    std::vector<unsigned> r2;
    r2.push_back(1);
    r2.push_back(2);
    const FJTsmCapacityResult<double> t3 = fj_tsm_capacity<double>(2u, f2, r2, x2);
    CHECK(t3.fcfs_available);
    CHECK(t3.Lfcfs ==
          doctest::Approx(2 * 2.0 * 3.0 / (0.25 * 3 + 2 * 0.25 * 2 + 2 * 0.25 * 5)).epsilon(1e-12));
}

TEST_CASE("the serialization blocking probability") {
    std::vector<double> Rs;
    Rs.push_back(0.2);
    Rs.push_back(0.3);
    const FJSerializationResult<double> s = fj_serialization<double>(Rs, 1.0, 5u);
    CHECK(s.P[0] == doctest::Approx(1 - std::pow(1 - 0.2 / 1.5, 4.0)).epsilon(1e-12));
}

TEST_CASE("the task-system makespan reproduces both worked examples") {
    Matrix<double> pred(2, 2, 0.0);
    Matrix<double> rate(2, 2, 0.0);
    rate(0, 0) = 1.0 / 10;
    rate(0, 1) = 1.0 / 15;
    rate(1, 0) = 1.0 / 20;
    rate(1, 1) = 1.0 / 30;
    const FJDagMakespanResult<double> d = fj_dag_makespan<double>(pred, rate);
    // The coupled two-task example of the survey gives 26.77 by hand, 80/3 exactly
    CHECK(d.C == doctest::Approx(80.0 / 3).epsilon(1e-12));
    CHECK(d.Cend[0] == doctest::Approx(40.0 / 3).epsilon(1e-12));
    CHECK(d.Cend[1] == doctest::Approx(70.0 / 3).epsilon(1e-12));
    Matrix<double> rate2(2, 2, 0.0);
    rate2(0, 0) = 1.0;
    rate2(0, 1) = 1.0;
    rate2(1, 0) = 0.1;
    rate2(1, 1) = 0.1;
    // The state-truncation example of the survey gives 10.1
    CHECK(fj_dag_makespan<double>(pred, rate2).C ==
          doctest::Approx(1 / 1.1 + (1 / 1.1) * 10 + (0.1 / 1.1) * 1).epsilon(1e-12));
}

TEST_CASE("the task-system makespan of a chain is the sum of the means") {
    Matrix<double> pc(3, 3, 0.0);
    pc(0, 1) = 1;
    pc(1, 2) = 1;
    Matrix<double> rc(3, 3, 0.0);
    for (int k = 0; k < 3; ++k) {
        rc(0, k) = 1.0;
        rc(1, k) = 2.0;
        rc(2, k) = 4.0;
    }
    const FJDagMakespanResult<double> d = fj_dag_makespan<double>(pc, rc);
    CHECK(d.C == doctest::Approx(1 + 0.5 + 0.25).epsilon(1e-12));
    CHECK(d.I[2] == doctest::Approx(1.5).epsilon(1e-12));
}
