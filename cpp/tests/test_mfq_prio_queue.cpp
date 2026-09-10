/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The fluid priority queue (line/api/mam/mfq_prio_queue.h).
 *
 * Oracles, in the order the task prescribes.
 *  (a) A collapse checked against an INDEPENDENT ported function. With a single
 *      class the priority queue is an ordinary fluid queue served at the
 *      constant rate d, and the sojourn time of a drop in that queue is exactly
 *      what mfq_sojourn returns for Rin = diag(R) and Rout = d I. The two reach
 *      it through disjoint constructions -- the priority path goes through the
 *      workload solve and the canonical similarity, mfq_sojourn through the
 *      Kronecker product and the TransformToOnes similarity -- and they agree to
 *      all fifteen digits. Nothing in either path resembles the other, which is
 *      what makes the check worth something.
 *  (b) Invariants: the CDFs are non-decreasing, bounded by one and approach it;
 *      the moment sequences satisfy E[X^2] >= E[X]^2; a HIGHER priority class
 *      has strictly smaller sojourn time and fluid level than a lower one; and
 *      the erlangization converges as its order grows, which is the accuracy
 *      knob of the distribution path and is measured rather than assumed.
 *      Separately, the degenerate model in which no class can ever exceed the
 *      service rate must return zero moments and unit distributions rather than
 *      throwing -- that path exists only because mfq_general_solve treats an
 *      empty up-drift set as an answer.
 *  (c) MATLAB, on a fully non-degenerate two-class instance, for all four
 *      measures and both classes.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mfq_prio_queue.h"
#include "line/api/mam/mfq_sojourn.h"

using line::Matrix;
using line::mam::FluidPrioOptions;
using line::mam::FluidPrioResult;
using line::mam::MeRepresentation;
using line::mam::mfq_prio_queue;
using line::mam::mfq_sojourn;

namespace {

Matrix<double> mat(const std::vector<std::vector<double>>& a) {
    Matrix<double> m(a.size(), a[0].size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) m(i, j) = a[i][j];
    return m;
}

/** k-th raw moment of an ME law: (-1)^k k! alpha A^-k e. */
double me_moment(const MeRepresentation<double>& r, unsigned k) {
    Matrix<double> P = line::eye<double>(r.A.rows());
    const Matrix<double> Ai = line::inverse(r.A);
    for (unsigned i = 0; i < k; ++i) P = line::matmul(P, Ai);
    const std::vector<double> v = line::mulvec(P, line::ones<double>(r.A.rows()));
    double s = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) s += r.alpha[i] * v[i];
    double f = 1.0;
    for (unsigned i = 2; i <= k; ++i) f *= static_cast<double>(i);
    return (k % 2 == 0 ? 1.0 : -1.0) * f * s;
}

const std::vector<std::vector<double>> kQ3 = {{-3.0, 2.0, 1.0}, {1.0, -2.0, 1.0},
                                              {2.0, 1.0, -3.0}};
/** Row 2 is the HIGHEST priority. Both sub-workloads are stable and non-empty. */
const std::vector<std::vector<double>> kR2 = {{0.5, 1.0, 0.5}, {2.5, 0.5, 0.2}};
const double kD = 2.0;

}  // namespace

// ---------------------------------------------------------------------------
// (a) the single-class collapse, against an independent ported function
// ---------------------------------------------------------------------------

TEST_CASE("a one-class priority queue is the fluid queue mfq_sojourn solves") {
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat({{2.5, 0.5, 0.2}});
    FluidPrioOptions o;
    o.stMoms = 3;
    const FluidPrioResult r = mfq_prio_queue(Q, R, kD, o);
    REQUIRE(r.stMoms.size() == 1u);

    Matrix<double> Rin(3, 3, 0.0), Rout(3, 3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        Rin(i, i) = R(0, i);
        Rout(i, i) = kD;
    }
    const MeRepresentation<double> me = mfq_sojourn(Q, Rin, Rout);
    for (unsigned k = 1; k <= 3; ++k) {
        INFO("moment ", k);
        CHECK(r.stMoms[0][k - 1] == doctest::Approx(me_moment(me, k)).epsilon(1e-12));
    }
    // and the values themselves, so a joint regression in both is still caught.
    CHECK(r.stMoms[0][0] == doctest::Approx(0.0755393362587057).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// (b) invariants
// ---------------------------------------------------------------------------

TEST_CASE("the priority measures are distributions and moments of one") {
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    FluidPrioOptions o;
    o.stMoms = 3;
    o.flMoms = 3;
    o.stDistr = {0.05, 0.2, 0.5, 1.0, 2.0, 8.0};
    o.flDistr = {0.05, 0.2, 0.5, 1.0, 2.0, 8.0};
    const FluidPrioResult r = mfq_prio_queue(Q, R, kD, o);
    REQUIRE(r.stDistr.size() == 2u);

    for (std::size_t k = 0; k < 2; ++k) {
        INFO("class ", k + 1);
        for (const std::vector<double>* cdf : {&r.stDistr[k], &r.flDistr[k]}) {
            for (std::size_t i = 0; i < cdf->size(); ++i) {
                CHECK((*cdf)[i] >= -1e-12);
                CHECK((*cdf)[i] <= 1.0 + 1e-9);
                if (i > 0) CHECK((*cdf)[i] >= (*cdf)[i - 1] - 1e-9);
            }
            // far out in the tail the mass is essentially all collected
            CHECK(cdf->back() > 0.95);
        }
        // Raw moments of a non-negative variable: E[X^2] >= E[X]^2.
        CHECK(r.stMoms[k][1] >= r.stMoms[k][0] * r.stMoms[k][0]);
        CHECK(r.flMoms[k][1] >= r.flMoms[k][0] * r.flMoms[k][0]);
        for (double v : r.stMoms[k]) CHECK(v > 0.0);
        for (double v : r.flMoms[k]) CHECK(v > 0.0);
    }
}

TEST_CASE("the higher priority class waits strictly less") {
    // Row 2 is the HIGHEST priority. If a port silently reversed the class
    // order this assertion is what fails, and it is the failure that matters:
    // reversing it produces numbers that look entirely plausible.
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    FluidPrioOptions o;
    o.stMoms = 1;
    o.flMoms = 1;
    const FluidPrioResult r = mfq_prio_queue(Q, R, kD, o);
    CHECK(r.stMoms[1][0] < r.stMoms[0][0]);
    CHECK(r.flMoms[1][0] < r.flMoms[0][0]);
    // and by a wide margin on this instance, not a rounding.
    CHECK(r.stMoms[0][0] / r.stMoms[1][0] > 5.0);
}

TEST_CASE("the erlangization converges as its order grows") {
    // The distribution path replaces the deterministic horizon by an Erlang of
    // order L, an O(1/L) approximation. Successive refinements must settle,
    // and the increments must shrink; this measures that rather than assuming
    // the default 200 is enough.
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    std::vector<double> vals;
    for (std::size_t L : {50u, 100u, 200u, 400u}) {
        FluidPrioOptions o;
        o.stDistr = {1.0};
        o.erlMaxOrder = L;
        o.classes = {1};  // the lower priority class, which uses erlangization
        vals.push_back(mfq_prio_queue(Q, R, kD, o).stDistr[0][0]);
    }
    for (std::size_t i = 1; i < vals.size(); ++i) {
        INFO("L index ", i, " value ", vals[i]);
        CHECK(vals[i] > 0.0);
        CHECK(vals[i] < 1.0);
    }
    // Successive increments shrink, and roughly halve as L doubles.
    const double d1 = std::fabs(vals[1] - vals[0]);
    const double d2 = std::fabs(vals[2] - vals[1]);
    const double d3 = std::fabs(vals[3] - vals[2]);
    CHECK(d2 < d1);
    CHECK(d3 < d2);
    CHECK(d1 / d2 > 1.5);
    CHECK(d2 / d3 > 1.5);
}

TEST_CASE("a queue that can never build is answered, not refused") {
    // Every class's aggregate input stays below the service rate, so no class
    // ever queues: all moments are zero and all distributions are one. This
    // reaches mfq_general_solve with an empty up-drift set, which returns the
    // degenerate law rather than throwing -- without that, the priority queue
    // would be unusable on a large class of legitimate inputs. MATLAB returns
    // the same, with its own float noise (3.1e-33, -2.3e-17, -1.1e-17 on the
    // first class), so these are compared ABSOLUTELY: a relative test against
    // an exact zero fails on one ulp.
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat({{1.5, 0.5, 0.0}, {0.5, 1.0, 0.5}});
    FluidPrioOptions o;
    o.flMoms = 3;
    o.stMoms = 3;
    o.stDistr = {0.2, 1.0};
    o.flDistr = {0.2, 1.0};
    FluidPrioResult r;
    REQUIRE_NOTHROW(r = mfq_prio_queue(Q, R, 3.0, o));
    for (std::size_t k = 0; k < 2; ++k) {
        INFO("class ", k + 1);
        for (double v : r.flMoms[k]) CHECK(std::fabs(v) < 1e-12);
        for (double v : r.stMoms[k]) CHECK(std::fabs(v) < 1e-12);
        for (double v : r.stDistr[k]) CHECK(std::fabs(v - 1.0) < 1e-12);
        for (double v : r.flDistr[k]) CHECK(std::fabs(v - 1.0) < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// (c) MATLAB
// ---------------------------------------------------------------------------

TEST_CASE("mfq_prio_queue agrees with MATLAB on a two-class instance") {
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    FluidPrioOptions o;
    o.flMoms = 3;
    o.stMoms = 3;
    o.stDistr = {0.2, 0.5, 1.0, 2.0};
    o.flDistr = {0.2, 0.5, 1.0, 2.0};
    const FluidPrioResult r = mfq_prio_queue(Q, R, kD, o);

    const std::vector<std::vector<double>> refFlMoms = {
        {0.50495597577952056, 0.75921280238646116, 1.7687521158620767},
        {0.079316303071640989, 0.03065071530447773, 0.017766833153347861}};
    const std::vector<std::vector<double>> refStMoms = {
        {0.70254744456281126, 1.5887870179580448, 5.6744524159564236},
        {0.075539336258705708, 0.014595578716417968, 0.0042301983698447296}};
    const std::vector<std::vector<double>> refStDistr = {
        {0.4502218005703616, 0.60508016979309764, 0.75631334846768916, 0.89829814597708757},
        {0.90135398689677682, 0.99557955149430732, 0.99997500934342454, 0.99999999920126936}};
    const std::vector<std::vector<double>> refFlDistr = {
        {0.4758706181121593, 0.66254828851360825, 0.82826941608775573, 0.95230465215356541},
        {0.85419360421455603, 0.96913473757005864, 0.99767926453451128, 0.99998687990529778}};

    for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t i = 0; i < 3; ++i) {
            INFO("class ", k + 1, " moment ", i + 1);
            CHECK(r.flMoms[k][i] == doctest::Approx(refFlMoms[k][i]).epsilon(1e-11));
            CHECK(r.stMoms[k][i] == doctest::Approx(refStMoms[k][i]).epsilon(1e-11));
        }
        for (std::size_t i = 0; i < 4; ++i) {
            INFO("class ", k + 1, " point ", i);
            CHECK(r.stDistr[k][i] == doctest::Approx(refStDistr[k][i]).epsilon(1e-11));
            CHECK(r.flDistr[k][i] == doctest::Approx(refFlDistr[k][i]).epsilon(1e-11));
        }
    }
}

TEST_CASE("mfq_prio_queue honours the classes option") {
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    FluidPrioOptions all, one;
    all.stMoms = 2;
    one.stMoms = 2;
    one.classes = {2};
    const FluidPrioResult ra = mfq_prio_queue(Q, R, kD, all);
    const FluidPrioResult ro = mfq_prio_queue(Q, R, kD, one);
    REQUIRE(ra.stMoms.size() == 2u);
    REQUIRE(ro.stMoms.size() == 1u);
    CHECK(ro.classes[0] == 2u);
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(ro.stMoms[0][i] == doctest::Approx(ra.stMoms[1][i]).epsilon(1e-13));
}

TEST_CASE("mfq_prio_queue rejects a malformed instance") {
    const Matrix<double> Q = mat(kQ3);
    const Matrix<double> R = mat(kR2);
    FluidPrioOptions o;
    o.stMoms = 1;
    CHECK_THROWS_AS(mfq_prio_queue(Q, R, 0.0, o), line::InputError);
    CHECK_THROWS_AS(mfq_prio_queue(Q, mat({{1.0, -1.0, 1.0}}), kD, o), line::InputError);
    CHECK_THROWS_AS(mfq_prio_queue(Q, mat({{1.0, 2.0}}), kD, o), line::InputError);
    FluidPrioOptions bad = o;
    bad.classes = {3};
    CHECK_THROWS_AS(mfq_prio_queue(Q, R, kD, bad), line::InputError);
}
