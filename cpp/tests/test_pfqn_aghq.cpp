/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Adaptive Gauss-Hermite quadrature of the simplex factor of the McKenna-Mitra
 * integral. The reference values come from MATLAB pfqn_aghq on the same models;
 * the structural identities (q = 1 reproduces pfqn_le, M = 1 is exact, an
 * all-zero Z is the Z = 0 branch) are asserted directly, because those are what
 * a port is most likely to get wrong.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_aghq.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_le.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

Matrix<double> demands3x2() {
    Matrix<double> L(3, 2);
    L(0, 0) = 1.0; L(0, 1) = 0.5;
    L(1, 0) = 0.7; L(1, 1) = 1.2;
    L(2, 0) = 0.3; L(2, 1) = 0.9;
    return L;
}

std::vector<double> pop23() {
    std::vector<double> N;
    N.push_back(2.0);
    N.push_back(3.0);
    return N;
}

/** pfqn_ca takes an integer population and a Matrix think time. */
std::vector<int> ipop23() {
    std::vector<int> N;
    N.push_back(2);
    N.push_back(3);
    return N;
}

}  // namespace

TEST_CASE("pfqn_aghq at one node is the logistic expansion") {
    // A single node sits at the mode with weight sqrt(2 pi), so the rule collapses
    // to Cas17 eq. (34) up to the tolerance of the shared fixed point.
    const Matrix<double> L = demands3x2();
    const std::vector<double> N = pop23();
    const std::vector<double> Z0;
    CHECK(pfqn_aghq(L, N, Z0, static_cast<std::size_t>(1)).lG ==
          doctest::Approx(pfqn_le(L, N, Z0).lG).epsilon(1e-9));
}

TEST_CASE("pfqn_aghq converges in q towards the exact constant") {
    const Matrix<double> L = demands3x2();
    const std::vector<double> N = pop23();
    const std::vector<double> Z0;
    const double exact = pfqn_ca(L, ipop23(), Matrix<double>()).lG;
    const double e2 = std::fabs(pfqn_aghq(L, N, Z0, static_cast<std::size_t>(2)).lG - exact);
    const double e7 = std::fabs(pfqn_aghq(L, N, Z0, static_cast<std::size_t>(7)).lG - exact);
    CHECK(e7 <= e2);
    CHECK(e7 < 0.02);
}

TEST_CASE("the rule is exact at a single station") {
    // M = 1 leaves a point mass on the simplex, so there is nothing left to close.
    Matrix<double> L(1, 2);
    L(0, 0) = 2.0;
    L(0, 1) = 3.0;
    std::vector<double> N;
    N.push_back(4.0);
    N.push_back(1.0);
    std::vector<int> Ni;
    Ni.push_back(4);
    Ni.push_back(1);
    const double exact = pfqn_ca(L, Ni, Matrix<double>()).lG;
    CHECK(pfqn_aghq(L, N).lG == doctest::Approx(exact).epsilon(1e-10));
}

TEST_CASE("an all-zero think time is the Z = 0 branch") {
    // pfqn_nc always passes sum(Z,1), so this is the shape a delay-free model
    // arrives in; the count of Laplaced directions must not disagree with it.
    const Matrix<double> L = demands3x2();
    const std::vector<double> N = pop23();
    const std::vector<double> Zzero(2, 0.0);
    CHECK(pfqn_aghq(L, N, Zzero).lG ==
          doctest::Approx(pfqn_aghq(L, N, std::vector<double>()).lG).epsilon(1e-12));
    CHECK(pfqn_le(L, N, Zzero).lG ==
          doctest::Approx(pfqn_le(L, N, std::vector<double>()).lG).epsilon(1e-12));
}

TEST_CASE("the rule instantiates at extended precision") {
    // The registry declares {Double, Real}; Real50 must compile and agree with the
    // double instantiation, which the cases above pin to the MATLAB reference.
    typedef line::Real50 R;
    Matrix<R> L(3, 2);
    const Matrix<double> Ld = demands3x2();
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t r = 0; r < 2; ++r) L(i, r) = R(Ld(i, r));
    std::vector<R> N;
    N.push_back(R(2));
    N.push_back(R(3));
    const R lG = pfqn_aghq(L, N, std::vector<R>()).lG;
    CHECK(static_cast<double>(lG) ==
          doctest::Approx(pfqn_aghq(demands3x2(), pop23(), std::vector<double>()).lG)
              .epsilon(1e-9));
}
