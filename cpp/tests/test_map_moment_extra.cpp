/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_factorial_moment, map_joint_moment, mmap_infgen, qbd_blocks_mapmap1.
 *
 * The goldens are native Python's on the same MAP. Two structural checks carry
 * more than the numbers:
 *
 *  - the first factorial moment IS the mean interarrival time, which
 *    `map_mean` computes by an unrelated route, so the two agreeing pins the
 *    embedded law and the resolvent together;
 *  - `joint(1,1)` differs from `mean^2` exactly when the MAP is correlated, and
 *    EQUALS IT FOR A POISSON PROCESS, whose interarrivals are independent. That
 *    identity is what exposed a defect in the reference: native Python puts
 *    `D1` between the two resolvents where the embedded KERNEL
 *    `P = (-D0)^-1 D1` belongs, so its joint moment is short by one factor of
 *    the mean and is dimensionally wrong. At rate 1.5 it returns 0.6667, the
 *    mean, where the answer is 0.4444. This port uses the kernel and the test
 *    pins the identity.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_moment_extra.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

mam::Map<double> corr_map() {
    mam::Map<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = -2.0;
    a.D0(0, 1) = 0.6;
    a.D0(1, 1) = -1.0;
    a.D1(0, 0) = 1.0;
    a.D1(0, 1) = 0.4;
    a.D1(1, 0) = 0.7;
    a.D1(1, 1) = 0.3;
    return a;
}

}  // namespace

TEST_CASE("the factorial moments match native Python") {
    const mam::Map<double> a = corr_map();
    CHECK(mam::map_factorial_moment(a, 1) == doctest::Approx(0.8585858586).epsilon(1e-9));
    CHECK(mam::map_factorial_moment(a, 2) == doctest::Approx(1.575757576).epsilon(1e-9));
    CHECK(mam::map_factorial_moment(a, 3) == doctest::Approx(4.515151515).epsilon(1e-9));

    // The first factorial moment IS the mean, reached by an unrelated route.
    CHECK(mam::map_factorial_moment(a, 1) == doctest::Approx(mam::map_mean(a)).epsilon(1e-9));
    // The zeroth is the total mass, one.
    CHECK(mam::map_factorial_moment(a, 0) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("a Poisson process pins the joint moment: it must be mean squared") {
    // Independent interarrivals, so E[X_n X_n+1] = E[X]^2 exactly. This is the
    // identity that exposed the reference's missing resolvent: native Python
    // returns 0.6667 here, the MEAN, against the correct 0.4444.
    mam::Map<double> p;
    p.D0 = Matrix<double>(1, 1, -1.5);
    p.D1 = Matrix<double>(1, 1, 1.5);
    const double m = mam::map_mean(p);
    CHECK(m == doctest::Approx(1.0 / 1.5).epsilon(1e-12));
    CHECK(mam::map_joint_moment(p, 1, 1) == doctest::Approx(m * m).epsilon(1e-9));
    CHECK(mam::map_joint_moment(p, 1, 1) == doctest::Approx(0.4444444444).epsilon(1e-8));
    // And it is dimensionally a time SQUARED: doubling the rate quarters it.
    mam::Map<double> q;
    q.D0 = Matrix<double>(1, 1, -3.0);
    q.D1 = Matrix<double>(1, 1, 3.0);
    CHECK(mam::map_joint_moment(q, 1, 1) ==
          doctest::Approx(mam::map_joint_moment(p, 1, 1) / 4.0).epsilon(1e-8));
}

TEST_CASE("a correlated MAP's joint moment is NOT the product of marginals") {
    const mam::Map<double> a = corr_map();
    // The kernel form. Native Python returns 0.8585858586 and 1.475353535 for
    // these, which are its own values short by one factor of the mean.
    CHECK(mam::map_joint_moment(a, 1, 1) == doctest::Approx(0.7376767677).epsilon(1e-8));
    CHECK(mam::map_joint_moment(a, 2, 1) == doctest::Approx(1.354444444).epsilon(1e-8));

    const double mc = mam::map_mean(a);
    CHECK(std::fabs(mam::map_joint_moment(a, 1, 1) - mc * mc) > 1e-6);
}

TEST_CASE("the MMAP generator closes") {
    const mam::Map<double> a = corr_map();
    std::vector<Matrix<double>> Dk;
    Dk.push_back(a.D1);
    const Matrix<double> Q = mam::mmap_infgen(a.D0, Dk);
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += Q(i, j);
        CHECK(std::fabs(s) < 1e-12);
    }

    // Splitting D1 in two marked halves leaves the generator unchanged: the
    // marking is a labelling of arrivals, not a change of the phase process.
    std::vector<Matrix<double>> half;
    Matrix<double> h1 = a.D1, h2 = a.D1;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            h1(i, j) *= 0.3;
            h2(i, j) *= 0.7;
        }
    half.push_back(h1);
    half.push_back(h2);
    const Matrix<double> Q2 = mam::mmap_infgen(a.D0, half);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(Q2(i, j) == doctest::Approx(Q(i, j)).epsilon(1e-12));
}

TEST_CASE("the MAP/MAP/1 QBD blocks form a generator") {
    const mam::Map<double> a = corr_map();
    Matrix<double> B, L, F;
    mam::qbd_blocks_mapmap1(a.D0, a.D1, a.D0, a.D1, &B, &L, &F);
    REQUIRE(B.rows() == 4u);

    // B + L + F is the generator of a level-independent QBD, so its rows sum
    // to zero -- which is the only check that catches a Kronecker factor
    // ordered the wrong way in ONE of the three blocks.
    for (std::size_t i = 0; i < 4; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 4; ++j) s += B(i, j) + L(i, j) + F(i, j);
        CHECK(std::fabs(s) < 1e-12);
    }
    // The forward block raises the level and carries only arrivals, so its
    // total mass is the arrival D1's, replicated over the service phases.
    double fsum = 0.0, d1sum = 0.0;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) fsum += F(i, j);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) d1sum += a.D1(i, j);
    CHECK(fsum == doctest::Approx(2.0 * d1sum).epsilon(1e-12));
}

TEST_CASE("the refusals are by name") {
    const mam::Map<double> a = corr_map();
    mam::Map<double> bad = a;
    bad.D0 = Matrix<double>(2, 3, 0.0);
    CHECK_THROWS_AS(mam::map_factorial_moment(bad, 1), line::InputError);
    CHECK_THROWS_AS(mam::map_joint_moment(bad, 1, 1), line::InputError);

    std::vector<Matrix<double>> wrong;
    wrong.push_back(Matrix<double>(3, 3, 0.1));
    CHECK_THROWS_AS(mam::mmap_infgen(a.D0, wrong), line::InputError);

    Matrix<double> B, L, F;
    CHECK_THROWS_AS(
        mam::qbd_blocks_mapmap1(Matrix<double>(2, 3, 0.0), a.D1, a.D0, a.D1, &B, &L, &F),
        line::InputError);
}
