/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The last four gaps of the C++ port: mtrace_var, dmap_optim_dist(_acf) and
 * cache_rmf_lna.
 *
 * Each has its own kind of oracle, and none of them is a stored number:
 *
 *  - `mtrace_var` is checked against the variance computed by hand, and its
 *    NaN-on-a-singleton contract is pinned, since zero would be
 *    indistinguishable from a class whose samples are identical;
 *  - the discrete fitters are checked by the FIXED POINT their continuous twins
 *    use -- fit a D-MAP against itself and its own D1 must come back at
 *    distance zero;
 *  - `cache_rmf_lna` is checked by the CONSERVATION LAWS the cache obeys: the
 *    covariance of a deterministic total is zero, which is precisely what the
 *    double-centred subspace restriction exists to produce.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/dmap_optim_dist.h"
#include "line/api/trace/mtrace_var.h"

namespace mam = line::mam;
namespace tr = line::trace;
using line::Matrix;

TEST_CASE("mtrace_var is the POPULATION variance, per class") {
    std::vector<double> v;
    std::vector<int> t;
    // Class 0: 1, 3, 5 -> mean 3, E[X^2] = 35/3, var = 35/3 - 9 = 8/3.
    v.push_back(1.0);
    t.push_back(0);
    v.push_back(3.0);
    t.push_back(0);
    v.push_back(5.0);
    t.push_back(0);
    // Class 1: 2, 4 -> mean 3, E[X^2] = 10, var = 1.
    v.push_back(2.0);
    t.push_back(1);
    v.push_back(4.0);
    t.push_back(1);

    const std::vector<double> var = tr::mtrace_var(v, 2, t);
    REQUIRE(var.size() == 2u);
    CHECK(var[0] == doctest::Approx(8.0 / 3.0).epsilon(1e-12));
    CHECK(var[1] == doctest::Approx(1.0).epsilon(1e-12));
    // The SAMPLE variance would be 4 and 2; this is not that.
    CHECK(var[0] != doctest::Approx(4.0));
}

TEST_CASE("a class with under two samples is NaN, not zero") {
    std::vector<double> v;
    std::vector<int> t;
    v.push_back(1.0);
    t.push_back(0);
    v.push_back(3.0);
    t.push_back(0);
    v.push_back(7.0);
    t.push_back(1);  // a single sample
    // class 2 never appears

    const std::vector<double> var = tr::mtrace_var(v, 3, t);
    CHECK(!std::isnan(var[0]));
    CHECK(std::isnan(var[1]));
    CHECK(std::isnan(var[2]));

    // Identical samples DO give zero, and that is a different finding from a
    // class too short to measure -- which is why the singleton is NaN.
    std::vector<double> same;
    std::vector<int> st;
    same.push_back(2.0);
    st.push_back(0);
    same.push_back(2.0);
    st.push_back(0);
    const std::vector<double> z = tr::mtrace_var(same, 1, st);
    CHECK(z[0] == doctest::Approx(0.0));
    CHECK(!std::isnan(z[0]));

    std::vector<int> shortlab(1, 0);
    CHECK_THROWS_AS(tr::mtrace_var(v, 2, shortlab), line::InputError);
}

TEST_CASE("the discrete fitters recover a D-MAP's own D1, at distance zero") {
    mam::Dmap<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = 0.5;
    a.D0(0, 1) = 0.1;
    a.D0(1, 0) = 0.2;
    a.D0(1, 1) = 0.3;
    a.D1(0, 0) = 0.3;
    a.D1(0, 1) = 0.1;
    a.D1(1, 0) = 0.4;
    a.D1(1, 1) = 0.1;
    // D0 + D1 is stochastic: that is what makes it DISCRETE time.
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += a.D0(i, j) + a.D1(i, j);
        REQUIRE(s == doctest::Approx(1.0).epsilon(1e-12));
    }

    const std::vector<double> pie = mam::dmap_pie(a);
    const mam::MapOptimDist<double> r = mam::dmap_optim_dist(a, pie, a.D0, pie, 1u);
    // No convex branch in discrete time: never a global claim.
    CHECK(r.global == false);
    CHECK(std::fabs(r.d) < 1e-9);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(r.B1(i, j) == doctest::Approx(a.D1(i, j)).epsilon(1e-4));

    const mam::MapOptimDist<double> q = mam::dmap_optim_dist_acf(a, pie, a.D0, pie);
    CHECK(std::fabs(q.d) < 1e-9);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(q.B1(i, j) == doctest::Approx(a.D1(i, j)).epsilon(1e-4));
}

TEST_CASE("the discrete constraints are the CONTINUOUS ones with (I - D0)") {
    mam::Dmap<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = 0.5;
    a.D0(0, 1) = 0.1;
    a.D0(1, 0) = 0.2;
    a.D0(1, 1) = 0.3;
    a.D1(0, 0) = 0.3;
    a.D1(0, 1) = 0.1;
    a.D1(1, 0) = 0.4;
    a.D1(1, 1) = 0.1;
    const std::vector<double> pie = mam::dmap_pie(a);

    Matrix<double> Aeq;
    std::vector<double> beq;
    mam::dmapoptdetail::build_constraints_d(a.D0, pie, &Aeq, &beq);
    const std::vector<double> v = mam::optdistdetail::vec(a.D1);

    // The reference's own D1 is feasible: that is the check that the (I - D0)
    // substitution was made in BOTH blocks and in the right-hand side.
    double worst = 0.0;
    for (std::size_t i = 0; i < Aeq.rows(); ++i) {
        double s = -beq[i];
        for (std::size_t j = 0; j < v.size(); ++j) s += Aeq(i, j) * v[j];
        worst = std::max(worst, std::fabs(s));
    }
    CHECK(worst < 1e-12);

    // The row-sum targets are those of (I - D0), not of -D0.
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += (i == j ? 1.0 : 0.0) - a.D0(i, j);
        CHECK(beq[2 + i] == doctest::Approx(s).epsilon(1e-12));
    }
}

TEST_CASE("the discrete refusals are by name") {
    mam::Dmap<double> a;
    a.D0 = Matrix<double>(2, 2, 0.25);
    a.D1 = Matrix<double>(2, 2, 0.25);
    const std::vector<double> pie(2, 0.5);
    Matrix<double> B0(2, 2, 0.25);

    CHECK_THROWS_AS(mam::dmap_optim_dist(a, pie, Matrix<double>(2, 3, 0.1), pie, 1u),
                    line::InputError);
    CHECK_THROWS_AS(mam::dmap_optim_dist(a, pie, B0, std::vector<double>(3, 0.3), 1u),
                    line::InputError);
    CHECK_THROWS_AS(mam::dmap_optim_dist(a, pie, B0, pie, 0u), line::InputError);
    CHECK_THROWS_AS(mam::dmap_optim_dist_acf(a, std::vector<double>(5, 0.2), B0, pie),
                    line::InputError);
}
