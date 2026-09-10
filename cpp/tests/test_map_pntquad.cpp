/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_pntquad and mmap_count_moment.
 *
 * map_pntquad reaches the same P_n(t) as map_pnt by a different route --
 * quadrature on the forward equations against uniformization -- so the two
 * agreeing to the integrator's tolerance is evidence neither is wrong. It is
 * the strongest check available for either, because no third implementation
 * exists to compare against.
 *
 * mmap_count_moment is checked against the identity that ties it to the
 * unmarked process: the per-class counting MEANS must sum to the mean of the
 * aggregate counting process, since every arrival carries exactly one class.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_pnt.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

double maxdiff(const Matrix<double>& A, const Matrix<double>& B) {
    double d = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) d = std::max(d, std::fabs(A(i, j) - B(i, j)));
    return d;
}

}  // namespace

TEST_CASE("map_pntquad agrees with the uniformization route on an Erlang MAP") {
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const double t = 1.3;
    const std::size_t na = 4;
    const std::vector<Matrix<double>> a = mam::map_pnt(m, na, t);
    const std::vector<Matrix<double>> b = mam::map_pntquad(m, na, t);
    REQUIRE(a.size() == b.size());
    for (std::size_t n = 0; n <= na; ++n) CHECK(maxdiff(a[n], b[n]) < 1e-7);
}

TEST_CASE("map_pntquad agrees with the uniformization route on an MMPP(2)") {
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -(3.0 + 0.2); m.D0(0, 1) = 0.2;
    m.D0(1, 0) = 0.1;          m.D0(1, 1) = -(0.5 + 0.1);
    m.D1(0, 0) = 3.0;
    m.D1(1, 1) = 0.5;

    const double t = 2.0;
    const std::size_t na = 6;
    const std::vector<Matrix<double>> a = mam::map_pnt(m, na, t);
    const std::vector<Matrix<double>> b = mam::map_pntquad(m, na, t);
    for (std::size_t n = 0; n <= na; ++n) CHECK(maxdiff(a[n], b[n]) < 1e-6);
}

TEST_CASE("map_pntquad reproduces the Poisson pmf") {
    const double lam = 1.5, t = 2.0;
    mam::Map<double> m;
    m.D0 = Matrix<double>(1, 1, -lam);
    m.D1 = Matrix<double>(1, 1, lam);
    const std::vector<Matrix<double>> P = mam::map_pntquad(m, 5, t);
    double logf = 0.0;
    for (std::size_t n = 0; n <= 5; ++n) {
        if (n > 0) logf += std::log(static_cast<double>(n));
        const double want = std::exp(-lam * t + static_cast<double>(n) * std::log(lam * t) - logf);
        CHECK(P[n](0, 0) == doctest::Approx(want).epsilon(1e-6));
    }
}

TEST_CASE("the per-class counting means sum to the aggregate counting mean") {
    // Two classes splitting an Erlang-2's arrivals 30/70.
    mam::Mmap<double> m;
    const mam::Map<double> e = mam::map_erlang(1.0, 2);
    m.D0 = e.D0;
    m.D1 = e.D1;
    m.Dc.assign(2, Matrix<double>(2, 2, 0.0));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            m.Dc[0](i, j) = 0.3 * e.D1(i, j);
            m.Dc[1](i, j) = 0.7 * e.D1(i, j);
        }

    const double t = 3.0;
    const std::vector<unsigned> orders(1, 1u);
    const Matrix<double> pc = mam::mmap_count_moment(m, t, orders);
    REQUIRE(pc.rows() == 1);
    REQUIRE(pc.cols() == 2);
    const std::vector<double> agg = mam::map_count_moment(e, t, orders);
    CHECK(pc(0, 0) + pc(0, 1) == doctest::Approx(agg[0]).epsilon(1e-8));

    // Each class's own mean count is its share of the aggregate rate over t.
    CHECK(pc(0, 0) == doctest::Approx(0.3 * agg[0]).epsilon(1e-6));
    CHECK(pc(0, 1) == doctest::Approx(0.7 * agg[0]).epsilon(1e-6));
}

TEST_CASE("mmap_count_moment refuses an unmarked process by name") {
    mam::Mmap<double> m;
    m.D0 = Matrix<double>(1, 1, -1.0);
    m.D1 = Matrix<double>(1, 1, 1.0);
    CHECK_THROWS_AS(mam::mmap_count_moment(m, 1.0, std::vector<unsigned>(1, 1u)),
                    line::InputError);
}
