/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mamap22_fit_bs_multiclass: the backward-plus-sigma fit of a MAMAP(2,2).
 *
 * ORACLE: ROUND-TRIP THROUGH THE EXACT BRANCH. The closed-form inverse is
 * exact whenever its answer lands in the unit box, so measuring a marked MAP's
 * own (p, B, S) and handing them back must return that marking and set the
 * `exact` flag. That pins the coefficient tables, the inverse and the assembly
 * in one check.
 *
 * The degenerate branches are driven deliberately, as in the F+B tests, because
 * a fit that silently took a degenerate branch would still satisfy the class
 * probabilities and look right.
 *
 * The unported repair is asserted as a NAMED refusal, not as a wrong number: a
 * caller that hits it must be told which solver is missing.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mamap22_fit_bs.h"
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

/** Mark a canonical AMAP(2) with the three per-flow class-1 shares. */
mam::Mmap<double> mark(const mam::Map<double>& a, int form, double q1, double q2, double q3) {
    mam::Mmap<double> m;
    m.D0 = a.D0;
    m.D1 = a.D1;
    m.Dc.assign(2, Matrix<double>(2, 2, 0.0));
    if (form == 1) {
        m.Dc[0](0, 0) = a.D1(0, 0) * q1;
        m.Dc[1](0, 0) = a.D1(0, 0) * (1.0 - q1);
    } else {
        m.Dc[0](0, 1) = a.D1(0, 1) * q1;
        m.Dc[1](0, 1) = a.D1(0, 1) * (1.0 - q1);
    }
    m.Dc[0](1, 0) = a.D1(1, 0) * q2;
    m.Dc[1](1, 0) = a.D1(1, 0) * (1.0 - q2);
    m.Dc[0](1, 1) = a.D1(1, 1) * q3;
    m.Dc[1](1, 1) = a.D1(1, 1) * (1.0 - q3);
    return m;
}

}  // namespace

TEST_CASE("the closed form inverts a marking it was given, in form 1") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    const mam::Mmap<double> src = mark(a, 1, 0.3, 0.55, 0.42);

    const std::vector<unsigned> ord(1, 1u);
    const std::vector<double> p = mam::mmap_pc(src);
    const std::vector<std::vector<double>> bm = mam::mmap_backward_moment(src, ord, true);
    std::vector<double> B(2, 0.0);
    for (std::size_t c = 0; c < 2; ++c) B[c] = bm[c][0];
    const Matrix<double> S = mam::mmap_sigma(src);

    const mam::Mamap22FitResult<double> r = mam::mamap22_fit_bs_multiclass(a, p, B, S);
    check_is_mmap(r.mmap);
    CHECK(r.exact == true);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(r.mmap.Dc[c](i, j) == doctest::Approx(src.Dc[c](i, j)).epsilon(1e-6));
    CHECK(r.fB[0] == doctest::Approx(B[0]).epsilon(1e-6));
    CHECK(r.fS(0, 0) == doctest::Approx(S(0, 0)).epsilon(1e-6));
}

TEST_CASE("the closed form inverts a marking it was given, in form 2") {
    const mam::Map<double> a = form2(2.0, 0.6, 0.5, 0.4);
    const mam::Mmap<double> src = mark(a, 2, 0.4, 0.35, 0.6);

    const std::vector<unsigned> ord(1, 1u);
    const std::vector<double> p = mam::mmap_pc(src);
    const std::vector<std::vector<double>> bm = mam::mmap_backward_moment(src, ord, true);
    std::vector<double> B(2, 0.0);
    for (std::size_t c = 0; c < 2; ++c) B[c] = bm[c][0];
    const Matrix<double> S = mam::mmap_sigma(src);

    const mam::Mamap22FitResult<double> r = mam::mamap22_fit_bs_multiclass(a, p, B, S);
    check_is_mmap(r.mmap);
    CHECK(r.exact == true);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(r.mmap.Dc[c](i, j) == doctest::Approx(src.Dc[c](i, j)).epsilon(1e-6));
    // Form 2 has no (1,1) arrival flow to mark.
    for (std::size_t c = 0; c < 2; ++c) CHECK(r.mmap.Dc[c](0, 0) == doctest::Approx(0.0));
}

TEST_CASE("the canonical phase-type branch delegates to the MAPH fitter") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.0);  // r2 = 0
    std::vector<double> p;
    p.push_back(0.3);
    p.push_back(0.7);
    std::vector<double> B(2, mam::map_mean(a));
    Matrix<double> S(2, 2, 0.25);

    const mam::Mamap22FitResult<double> r = mam::mamap22_fit_bs_multiclass(a, p, B, S);
    check_is_mmap(r.mmap);
    for (std::size_t c = 0; c < 2; ++c) CHECK(r.mmap.Dc[c](1, 1) == doctest::Approx(0.0));
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.3).epsilon(1e-4));
}

TEST_CASE("the degenerate phase-type branch splits every flow by p") {
    const mam::Map<double> a = form2(2.0, 0.6, 1.0, 0.0);  // r1 = 1, r2 = 0
    std::vector<double> p;
    p.push_back(0.25);
    p.push_back(0.75);
    std::vector<double> B(2, 1.0);
    Matrix<double> S(2, 2, 0.0625);

    const mam::Mamap22FitResult<double> r = mam::mamap22_fit_bs_multiclass(a, p, B, S);
    check_is_mmap(r.mmap);
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.25).epsilon(1e-6));
    CHECK(pc[1] == doctest::Approx(0.75).epsilon(1e-6));
}

TEST_CASE("an infeasible closed form is refused by name, and adjust=false clamps") {
    // Targets far outside what the marking can reach.
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p;
    p.push_back(0.5);
    p.push_back(0.5);
    std::vector<double> B(2, 50.0);   // a backward moment the form cannot carry
    Matrix<double> S(2, 2, 0.9);

    CHECK_THROWS_AS(mam::mamap22_fit_bs_multiclass(a, p, B, S), line::UnsupportedError);

    // With the repair declined, the clamped closed form comes back as a valid
    // marked MAP; it simply does not match the targets.
    const mam::Mamap22FitResult<double> r =
        mam::mamap22_fit_bs_multiclass(a, p, B, S, std::vector<double>(), std::vector<double>(),
                                       false);
    check_is_mmap(r.mmap);
    CHECK(r.exact == false);
}

TEST_CASE("the two-class restriction and the canonical form are enforced by name") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p3(3, 1.0 / 3.0), B3(3, 1.0);
    Matrix<double> S(2, 2, 0.25);
    CHECK_THROWS_AS(mam::mamap22_fit_bs_multiclass(a, p3, B3, S), line::InputError);

    std::vector<double> p(2, 0.5), B(2, 1.0);
    mam::Map<double> both = form1(2.0, 0.6, 0.5, 0.4);
    both.D1(0, 1) = 0.05;
    CHECK_THROWS_AS(mam::mamap22_fit_bs_multiclass(both, p, B, S), line::InputError);

    const mam::Map<double> third = mam::map_erlang(1.0, 3);
    CHECK_THROWS_AS(mam::mamap22_fit_bs_multiclass(third, p, B, S), line::InputError);
}

TEST_CASE("mamap22_fit_gamma_bs round-trips the descriptors of a marked AMAP(2)") {
    // the source marking is exactly representable, so the gamma_bs sweep over
    // the AMAP(2) forms must land on a fit reproducing P exactly and (B, S)
    // to fitting accuracy
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    const mam::Mmap<double> src = mark(a, 1, 0.3, 0.55, 0.42);

    const mam::Mmap<double> fit = mam::mamap22_fit_gamma_bs_mmap(src);
    check_is_mmap(fit);
    const std::vector<double> p = mam::mmap_pc(src);
    const std::vector<double> pf = mam::mmap_pc(fit);
    CHECK(pf[0] == doctest::Approx(p[0]).epsilon(1e-6));
    CHECK(pf[1] == doctest::Approx(p[1]).epsilon(1e-6));
    const std::vector<unsigned> ord(1, 1u);
    const std::vector<std::vector<double>> bs = mam::mmap_backward_moment(src, ord, true);
    const std::vector<std::vector<double>> bf = mam::mmap_backward_moment(fit, ord, true);
    CHECK(bf[0][0] == doctest::Approx(bs[0][0]).epsilon(1e-4));
    const Matrix<double> Ss = mam::mmap_sigma(src);
    const Matrix<double> Sf = mam::mmap_sigma(fit);
    CHECK(Sf(0, 0) == doctest::Approx(Ss(0, 0)).epsilon(1e-4));
    // the underlying inter-arrival moments are matched by the AMAP(2) fit
    CHECK(mam::map_moment(fit.map(), 1u) ==
          doctest::Approx(mam::map_moment(src.map(), 1u)).epsilon(1e-6));
}

TEST_CASE("mamap22_fit_gamma_bs_trace matches the class probabilities of the trace") {
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 3000; ++i) {
        Tv.push_back(0.4 + 0.6 * static_cast<double>((i * 7919) % 17) / 17.0);
        A.push_back(static_cast<int>(i % 3 == 0) + 1);
    }
    const mam::Mmap<double> m = mam::mamap22_fit_gamma_bs_trace(Tv, A);
    check_is_mmap(m);
    const std::vector<double> pc = mam::mmap_pc(m);
    CHECK(pc[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-3));
    double sm = 0.0;
    for (std::size_t i = 0; i < Tv.size(); ++i) sm += Tv[i];
    sm /= static_cast<double>(Tv.size());
    CHECK(mam::map_mean(m.map()) == doctest::Approx(sm).epsilon(0.05));

    CHECK_THROWS_AS(mam::mamap22_fit_gamma_bs_trace(std::vector<double>(), std::vector<int>()),
                    line::InputError);
}
