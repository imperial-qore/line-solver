/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mamap22_fit_fs_multiclass: the forward-plus-sigma fit of a MAMAP(2,2).
 *
 * Same oracle as the backward twin -- ROUND-TRIP through the exact branch --
 * plus the two things that are specific to this fitter:
 *  - the canonical phase-type branch substitutes the ordinary mean for the
 *    backward targets it has no forward equivalent of, and must SAY SO; the
 *    test asserts the diagnostic is present, because a caller who asked for a
 *    forward fit and silently got a mean-matched backward one has no other way
 *    to find out;
 *  - the two unported YALMIP paths are asserted as named refusals.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mamap22_fit_fs.h"
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
    const Matrix<double> fm = mam::mmap_forward_moment(src, ord, true);
    std::vector<double> F(2, 0.0);
    for (std::size_t c = 0; c < 2; ++c) F[c] = fm(c, 0);
    const Matrix<double> S = mam::mmap_sigma(src);

    const mam::Mamap22FsFitResult<double> r = mam::mamap22_fit_fs_multiclass(a, p, F, S);
    check_is_mmap(r.mmap);
    CHECK(r.exact == true);
    CHECK(r.warning.empty());
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(r.mmap.Dc[c](i, j) == doctest::Approx(src.Dc[c](i, j)).epsilon(1e-6));
    CHECK(r.fF[0] == doctest::Approx(F[0]).epsilon(1e-6));
    CHECK(r.fS(0, 0) == doctest::Approx(S(0, 0)).epsilon(1e-6));
}

TEST_CASE("the closed form inverts a marking it was given, in form 2") {
    const mam::Map<double> a = form2(2.0, 0.6, 0.5, 0.4);
    const mam::Mmap<double> src = mark(a, 2, 0.4, 0.35, 0.6);

    const std::vector<unsigned> ord(1, 1u);
    const std::vector<double> p = mam::mmap_pc(src);
    const Matrix<double> fm = mam::mmap_forward_moment(src, ord, true);
    std::vector<double> F(2, 0.0);
    for (std::size_t c = 0; c < 2; ++c) F[c] = fm(c, 0);
    const Matrix<double> S = mam::mmap_sigma(src);

    const mam::Mamap22FsFitResult<double> r = mam::mamap22_fit_fs_multiclass(a, p, F, S);
    check_is_mmap(r.mmap);
    CHECK(r.exact == true);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                CHECK(r.mmap.Dc[c](i, j) == doctest::Approx(src.Dc[c](i, j)).epsilon(1e-6));
}

TEST_CASE("the canonical phase-type branch says it substituted its targets") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.0);  // r2 = 0
    std::vector<double> p;
    p.push_back(0.3);
    p.push_back(0.7);
    std::vector<double> F(2, mam::map_mean(a));
    Matrix<double> S(2, 2, 0.25);

    const mam::Mamap22FsFitResult<double> r = mam::mamap22_fit_fs_multiclass(a, p, F, S);
    check_is_mmap(r.mmap);
    // The reference warns here; the port carries the warning out.
    CHECK(!r.warning.empty());
    CHECK(r.warning.find("ordinary mean") != std::string::npos);
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.3).epsilon(1e-4));
}

TEST_CASE("the degenerate phase-type branch splits every flow by p") {
    const mam::Map<double> a = form2(2.0, 0.6, 1.0, 0.0);
    std::vector<double> p;
    p.push_back(0.2);
    p.push_back(0.8);
    std::vector<double> F(2, 1.0);
    Matrix<double> S(2, 2, 0.04);

    const mam::Mamap22FsFitResult<double> r = mam::mamap22_fit_fs_multiclass(a, p, F, S);
    check_is_mmap(r.mmap);
    const std::vector<double> pc = mam::mmap_pc(r.mmap);
    CHECK(pc[0] == doctest::Approx(0.2).epsilon(1e-6));
}

TEST_CASE("the unported repairs are refused by name") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p(2, 0.5), F(2, 50.0);
    Matrix<double> S(2, 2, 0.9);
    CHECK_THROWS_AS(mam::mamap22_fit_fs_multiclass(a, p, F, S), line::UnsupportedError);

    // Declining the repair returns the clamped closed form as a valid MAP.
    const mam::Mamap22FsFitResult<double> r = mam::mamap22_fit_fs_multiclass(
        a, p, F, S, std::vector<double>(), std::vector<double>(), false);
    check_is_mmap(r.mmap);
    CHECK(r.exact == false);

    // The sigma-weighted arm of the gamma < 0 degeneracy IS ported: its
    // feasible set is an interval and its objective is (x - S11)^2, so the
    // projection is the exact global optimum of the reference's YALMIP program.
    const mam::Map<double> b = form2(2.0, 0.6, 0.5, 0.0);
    std::vector<double> fw;
    fw.push_back(1.0);
    fw.push_back(2.0);  // sigma weighted above the forward moment
    const mam::Mamap22FsFitResult<double> s = mam::mamap22_fit_fs_multiclass(
        b, p, F, S, std::vector<double>(), fw);
    check_is_mmap(s.mmap);
    // S11 = 0.9 is above the p1^2 = 0.25 ceiling, so it is clamped onto it and
    // the marking degenerates to q = p, which reproduces S11 = p1^2 exactly.
    CHECK(s.fS(0, 0) == doctest::Approx(0.25).epsilon(1e-9));
}

TEST_CASE("the two-class restriction and the canonical form are enforced by name") {
    const mam::Map<double> a = form1(2.0, 0.6, 0.5, 0.4);
    std::vector<double> p3(3, 1.0 / 3.0), F3(3, 1.0);
    Matrix<double> S(2, 2, 0.25);
    CHECK_THROWS_AS(mam::mamap22_fit_fs_multiclass(a, p3, F3, S), line::InputError);

    std::vector<double> p(2, 0.5), F(2, 1.0);
    const mam::Map<double> third = mam::map_erlang(1.0, 3);
    CHECK_THROWS_AS(mam::mamap22_fit_fs_multiclass(third, p, F, S), line::InputError);
}
