/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mamap2m_fit_fb_multiclass, mamap2m_fit_gamma_fb and mamap2m_fit_trace.
 *
 * ORACLES, again chosen to be independent of the construction:
 *  - the class probabilities are an equality constraint of every branch and
 *    must come back exactly;
 *  - the inter-arrival law belongs to the AMAP and must survive the marking;
 *  - the forward and backward moments are the objective, so they are checked as
 *    an approach and, where the branch identifies only one of them, only that
 *    one is checked at all;
 *  - the result must be a marked MAP: Dc non-negative, sum_c Dc = D1, rows of
 *    D0 + D1 summing to zero.
 *
 * THE BRANCH SELECTION IS PART OF THE CONTRACT, so each degenerate branch is
 * driven deliberately rather than left to chance: a fit that silently took the
 * Poisson branch would satisfy every moment check above and still be wrong.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mamap2m_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

void check_is_mmap(const mam::Mmap<double>& m) {
    const std::size_t n = m.order();
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-9);
            double dc = 0.0;
            for (std::size_t c = 0; c < m.classes(); ++c) {
                CHECK(m.Dc[c](i, j) >= -1e-9);
                dc += m.Dc[c](i, j);
            }
            CHECK(m.D1(i, j) == doctest::Approx(dc).epsilon(1e-9));
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-8);
    }
}

/** First canonical acyclic form: D1(1,2) = 0, a positive decay. */
mam::Map<double> form1(double mu1, double mu2, double r1, double r2) {
    mam::Map<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = -mu1;
    a.D0(0, 1) = mu1 * r1;
    a.D0(1, 1) = -mu2;
    a.D1(0, 0) = mu1 * (1.0 - r1);
    a.D1(1, 0) = mu2 * (1.0 - r2);
    a.D1(1, 1) = mu2 * r2;
    return a;
}

/** Second canonical acyclic form: D1(1,1) = 0, a negative decay. */
mam::Map<double> form2(double mu1, double mu2, double r1, double r2) {
    mam::Map<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = -mu1;
    a.D0(0, 1) = mu1 * r1;
    a.D0(1, 1) = -mu2;
    a.D1(0, 1) = mu1 * (1.0 - r1);
    a.D1(1, 0) = mu2 * (1.0 - r2);
    a.D1(1, 1) = mu2 * r2;
    return a;
}

}  // namespace

TEST_CASE("the general branch reproduces the class law and the timing") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p;
    p.push_back(0.3);
    p.push_back(0.7);
    const double mean = mam::map_mean(a);
    std::vector<double> F(2, mean), B(2, mean);

    const mam::Mamap2mFitResult<double> r = mam::mamap2m_fit_fb_multiclass(a, p, F, B);
    check_is_mmap(r.mmap);
    CHECK(r.mmap.order() == 2);
    CHECK(r.mmap.classes() == 2);

    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.3).epsilon(1e-5));
    CHECK(pc[1] == doctest::Approx(0.7).epsilon(1e-5));

    // The marking cannot move the inter-arrival law.
    CHECK(mam::map_mean(r.mmap.map()) == doctest::Approx(mean).epsilon(1e-12));
    CHECK(mam::map_scv(r.mmap.map()) == doctest::Approx(mam::map_scv(a)).epsilon(1e-12));
    REQUIRE(r.fF.size() == 2);
    REQUIRE(r.fB.size() == 2);
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(r.fF[c] > 0.0);
        CHECK(r.fB[c] > 0.0);
    }
}

TEST_CASE("the second canonical form is accepted and marked in its own flows") {
    const mam::Map<double> a = form2(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p;
    p.push_back(0.45);
    p.push_back(0.55);
    const double mean = mam::map_mean(a);
    std::vector<double> F(2, mean), B(2, mean);

    const mam::Mamap2mFitResult<double> r = mam::mamap2m_fit_fb_multiclass(a, p, F, B);
    check_is_mmap(r.mmap);
    // Form 2 has no (1,1) arrival flow, so no class may claim one.
    for (std::size_t c = 0; c < 2; ++c) CHECK(r.mmap.Dc[c](0, 0) == doctest::Approx(0.0));
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.45).epsilon(1e-5));
}

TEST_CASE("the canonical phase-type branch delegates to the MAPH fitter") {
    // r2 = 0 in form 1: the MAP is an APH(2), so the backward moments are what
    // gets fitted and the (2,2) flow must be empty in every class.
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.0);
    std::vector<double> p;
    p.push_back(0.25);
    p.push_back(0.75);
    const double mean = mam::map_mean(a);
    std::vector<double> F(2, mean), B(2, mean);

    const mam::Mamap2mFitResult<double> r = mam::mamap2m_fit_fb_multiclass(a, p, F, B);
    check_is_mmap(r.mmap);
    for (std::size_t c = 0; c < 2; ++c) CHECK(r.mmap.Dc[c](1, 1) == doctest::Approx(0.0));
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.25).epsilon(1e-4));
    CHECK(pc[1] == doctest::Approx(0.75).epsilon(1e-4));
}

TEST_CASE("the degenerate phase-type branch splits every flow by p") {
    // Form 2 with r2 = 0 and r1 = 1.
    const mam::Map<double> a = form2(2.0, 0.6, 1.0, 0.0);
    std::vector<double> p;
    p.push_back(0.2);
    p.push_back(0.8);
    std::vector<double> F(2, 1.0), B(2, 1.0);

    const mam::Mamap2mFitResult<double> r = mam::mamap2m_fit_fb_multiclass(a, p, F, B);
    check_is_mmap(r.mmap);
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.2).epsilon(1e-6));
    CHECK(pc[1] == doctest::Approx(0.8).epsilon(1e-6));
}

TEST_CASE("a collapsed AMAP gives a marked Poisson process") {
    // h1 - h2 + h2 r1 = 0 trips the Poisson guard while leaving the chain
    // irreducible: h1 = 1, h2 = 2, r1 = 1/2 gives 1 - 2 + 1 = 0. The other two
    // triggers (r1 = 0, r2 = 1) also fire but leave a state absorbing, and
    // map_mean refuses such a generator before the guard is ever reached.
    const mam::Map<double> a = form1(1.0, 0.5, 0.5, 0.4);
    std::vector<double> p;
    p.push_back(0.4);
    p.push_back(0.6);
    std::vector<double> F(2, 1.0), B(2, 1.0);

    const mam::Mamap2mFitResult<double> r = mam::mamap2m_fit_fb_multiclass(a, p, F, B);
    CHECK(r.mmap.order() == 1);
    check_is_mmap(r.mmap);
    CHECK(mam::map_mean(r.mmap.map()) == doctest::Approx(mam::map_mean(a)).epsilon(1e-12));
    CHECK(mam::map_scv(r.mmap.map()) == doctest::Approx(1.0).epsilon(1e-12));
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.4).epsilon(1e-12));
}

TEST_CASE("a non-canonical or cyclic underlying MAP is refused by name") {
    std::vector<double> p(2, 0.5), F(2, 1.0), B(2, 1.0);
    mam::Map<double> cyc = form1(2.0, 0.6, 0.5, 0.4);
    cyc.D0(1, 0) = 0.1;
    CHECK_THROWS_AS(mam::mamap2m_fit_fb_multiclass(cyc, p, F, B), line::InputError);

    mam::Map<double> both = form1(2.0, 0.6, 0.5, 0.4);
    both.D1(0, 1) = 0.05;  // neither D1(1,1) nor D1(1,2) is zero
    CHECK_THROWS_AS(mam::mamap2m_fit_fb_multiclass(both, p, F, B), line::InputError);

    const mam::Map<double> third = mam::map_erlang(1.0, 3);
    CHECK_THROWS_AS(mam::mamap2m_fit_fb_multiclass(third, p, F, B), line::InputError);
}

TEST_CASE("mamap2m_fit_gamma_fb keeps the inter-arrival moments it was given") {
    const double M1 = 1.0, scv = 2.0;
    const double M2 = (1.0 + scv) * M1 * M1;
    const double M3 = 3.0 * 1.5 * M2 * M2 / M1;
    const double gamma = 0.3;
    std::vector<double> p;
    p.push_back(0.35);
    p.push_back(0.65);
    std::vector<double> F(2, M1), B(2, M1);

    const mam::Mmap<double> m = mam::mamap2m_fit_gamma_fb(M1, M2, M3, gamma, p, F, B);
    check_is_mmap(m);
    CHECK(mam::map_mean(m.map()) == doctest::Approx(M1).epsilon(1e-5));
    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(0.35).epsilon(1e-4));
    CHECK(pc[1] == doctest::Approx(0.65).epsilon(1e-4));
}

TEST_CASE("mamap2m_fit_trace reads its descriptors from the labelled trace") {
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 3000; ++i) {
        Tv.push_back(0.4 + 0.6 * static_cast<double>((i * 7919) % 17) / 17.0);
        A.push_back(static_cast<int>(i % 3 == 0) + 1);
    }
    const mam::Mmap<double> m = mam::mamap2m_fit_trace(Tv, A);
    check_is_mmap(m);
    double sm = 0.0;
    for (std::size_t i = 0; i < Tv.size(); ++i) sm += Tv[i];
    sm /= static_cast<double>(Tv.size());
    CHECK(mam::map_mean(m.map()) == doctest::Approx(sm).epsilon(0.05));

    CHECK_THROWS_AS(mam::mamap2m_fit_trace(std::vector<double>(), std::vector<int>()),
                    line::InputError);
}

TEST_CASE("mamap2m_fit_gamma_fb_trace fits through the (F, B) pair alone") {
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 3000; ++i) {
        Tv.push_back(0.4 + 0.6 * static_cast<double>((i * 7919) % 17) / 17.0);
        A.push_back(static_cast<int>(i % 3 == 0) + 1);
    }
    const mam::Mmap<double> m = mam::mamap2m_fit_gamma_fb_trace(Tv, A);
    check_is_mmap(m);
    double sm = 0.0;
    for (std::size_t i = 0; i < Tv.size(); ++i) sm += Tv[i];
    sm /= static_cast<double>(Tv.size());
    CHECK(mam::map_mean(m.map()) == doctest::Approx(sm).epsilon(0.05));
    // the class probabilities are matched exactly by every descriptor pair;
    // class 1 labels the i % 3 != 0 arrivals, two thirds of the trace
    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-3));

    // sigma-preferring weights exercise the (F,S)/(B,S) branches of mamap2m_fit
    std::vector<double> w;
    w.push_back(1.0);
    w.push_back(1.0);
    w.push_back(5.0);
    const mam::Mmap<double> ms = mam::mamap2m_fit_trace(Tv, A, w);
    check_is_mmap(ms);
    const std::vector<double> pcs = mam::mmap_pc(ms);
    CHECK(pcs[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-3));
}
