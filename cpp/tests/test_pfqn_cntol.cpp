/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Chandy-Neuse population-scaled termination test (pfqn_cntol).
 *
 * The cutoff 1/(4000 + 16*sum(N)) and the metric max_{i,r}|dQ(i,r)|/N_r it is
 * compared against are published in K. M. Chandy, D. Neuse, "Linearizer: A
 * Heuristic Algorithm for Queuing Network Models of Computing Systems",
 * Commun. ACM 25(2):126-134, 1982, p.129 and appendix; LQNS runs the same
 * expression from SchweitzerCommon. It is opt-in here: the sentinel is a NaN
 * tolerance and the default tolerances are untouched, so the figures below
 * double as a guard that the sentinel path is the only thing it changes.
 *
 * Every reference value is MATLAB's, printed by pfqn_linearizer / pfqn_bs with
 * tol = 'cn' on the same model at %.12f. They are asserted at 1e-9, which is
 * looser than the observed agreement but is all a fixed point stopped at a
 * 2.4e-4 cutoff entitles one to claim.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_cntol.h"
#include "line/api/pfqn/pfqn_linearizer.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

Matrix<double> demands() {
    Matrix<double> L(3, 2);
    L(0, 0) = 1.0;  L(0, 1) = 2.0;
    L(1, 0) = 3.0;  L(1, 1) = 1.0;
    L(2, 0) = 0.5;  L(2, 1) = 0.5;
    return L;
}

const std::vector<int> kN{8, 1};
const std::vector<double> kNd{8.0, 1.0};

std::vector<SchedStrategy> allPS() {
    return std::vector<SchedStrategy>(3, SchedStrategy::PS);
}

}  // namespace

TEST_CASE("pfqn_cntol: the cutoff is the published expression") {
    CHECK(pfqn_cntol_total(10.0) == doctest::Approx(1.0 / (4000.0 + 160.0)).epsilon(1e-15));
    CHECK(pfqn_cntol(kN) == doctest::Approx(1.0 / (4000.0 + 16.0 * 9.0)).epsilon(1e-15));
    CHECK(pfqn_cntol(kNd) == doctest::Approx(1.0 / (4000.0 + 16.0 * 9.0)).epsilon(1e-15));
    // Below 0.00025 even at a population of one, as the paper's appendix states.
    CHECK(pfqn_cntol_total(1.0) < 0.00025);
    // Decreasing in the population, which is the point of the scaling.
    CHECK(pfqn_cntol_total(100.0) < pfqn_cntol_total(10.0));
}

TEST_CASE("pfqn_cntol: NaN is the sentinel and nothing else is") {
    CHECK(is_cntol(std::numeric_limits<double>::quiet_NaN()));
    CHECK_FALSE(is_cntol(1e-8));
    CHECK_FALSE(is_cntol(0.0));
}

TEST_CASE("pfqn_linearizer under the Chandy-Neuse test reproduces MATLAB") {
    const Matrix<double> L = demands();
    const Matrix<double> Z;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const LinearizerResult<double> cn = pfqn_linearizer(L, kN, Z, allPS(), nan, 1000, Matrix<double>());
    CHECK(cn.X[0] == doctest::Approx(0.304718171165).epsilon(1e-9));
    CHECK(cn.X[1] == doctest::Approx(0.084038879640).epsilon(1e-9));
    CHECK(cn.Q(0, 0) == doctest::Approx(0.555262542605).epsilon(1e-9));
    CHECK(cn.Q(1, 0) == doctest::Approx(7.255702544745).epsilon(1e-9));
    CHECK(cn.Q(2, 0) == doctest::Approx(0.189034912650).epsilon(1e-9));
    CHECK(cn.Q(0, 1) == doctest::Approx(0.251901005767).epsilon(1e-9));
    CHECK(cn.Q(1, 1) == doctest::Approx(0.697720815306).epsilon(1e-9));
    CHECK(cn.Q(2, 1) == doctest::Approx(0.050378178927).epsilon(1e-9));

    // The looser published cutoff must stop earlier than the 1e-8 default, and the two
    // answers must differ: an identical result would mean the sentinel did nothing.
    const LinearizerResult<double> def =
        pfqn_linearizer(L, kN, Z, allPS(), 1e-8, 1000, Matrix<double>());
    CHECK(cn.totiter < def.totiter);
    CHECK(std::fabs(cn.X[0] - def.X[0]) > 1e-12);
    CHECK(cn.X[0] == doctest::Approx(def.X[0]).epsilon(1e-4));
}

TEST_CASE("pfqn_bs under the Chandy-Neuse test reproduces MATLAB") {
    const Matrix<double> L = demands();
    const std::vector<double> Z(2, 0.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const AmvaResult<double> cn = pfqn_bs(L, kNd, Z, std::vector<AmvaSched>(), nan, 1000);
    CHECK(cn.XN[0] == doctest::Approx(0.301066244185).epsilon(1e-9));
    CHECK(cn.XN[1] == doctest::Approx(0.083878140075).epsilon(1e-9));

    const AmvaResult<double> def = pfqn_bs(L, kNd, Z, std::vector<AmvaSched>(), 1e-6, 1000);
    CHECK(cn.iterations < def.iterations);
    CHECK(cn.XN[0] == doctest::Approx(def.XN[0]).epsilon(1e-4));
}

TEST_CASE("pfqn_cntol: an empty class does not divide by zero") {
    // An empty class has N_r = 0, so the published metric would divide by zero
    // there; it must be skipped rather than guarded downstream.
    const Matrix<double> L = demands();
    const std::vector<double> Nz{8.0, 0.0};
    const std::vector<double> Z(2, 0.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const AmvaResult<double> r = pfqn_bs(L, Nz, Z, std::vector<AmvaSched>(), nan, 1000);
    for (std::size_t cls = 0; cls < 2; ++cls) CHECK(std::isfinite(r.XN[cls]));
    for (std::size_t i = 0; i < r.QN.rows(); ++i)
        for (std::size_t cls = 0; cls < r.QN.cols(); ++cls) CHECK(std::isfinite(r.QN(i, cls)));
    CHECK(r.XN[1] == doctest::Approx(0.0).epsilon(1e-15));
}
