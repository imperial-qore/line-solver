/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * AQL approximate MVA. The reference values come from MATLAB pfqn_aql on the
 * same model; the operational laws and the comparison against exact MVA and
 * against Bard-Schweitzer are asserted directly.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_aql.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_mva.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

Matrix<double> demands() {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5; L(0, 1) = 0.3;
    L(1, 0) = 0.4; L(1, 1) = 0.6;
    return L;
}

}  // namespace

TEST_CASE("AQL matches the MATLAB reference values") {
    // MATLAB: pfqn_aql([0.5 0.3;0.4 0.6],[4 3],[1 2]) gives
    // X = [1.141921860520, 0.648077557030], Q(1,1) = 1.437775148488,
    // Q(2,2) = 1.145884145121, in 28 iterations. The port takes the same 28
    // iterations and agrees to 9.3e-9 relative at worst. The assertions sit at
    // the tolerance the method actually guarantees (TOL = 1e-7 on the relative
    // change of Q), not tighter: the last digits of a fixed point stopped on a
    // tolerance are an artifact of the stopping test, not a property of the
    // algorithm.
    const std::vector<double> N{4.0, 3.0};
    const std::vector<double> Z{1.0, 2.0};
    auto a = pfqn_aql(demands(), N, Z);
    CHECK(a.converged);
    CHECK(a.XN[0] == doctest::Approx(1.141921860520).epsilon(1e-7));
    CHECK(a.XN[1] == doctest::Approx(0.648077557030).epsilon(1e-7));
    CHECK(a.QN(0, 0) == doctest::Approx(1.437775148488).epsilon(1e-7));
    CHECK(a.QN(1, 1) == doctest::Approx(1.145884145121).epsilon(1e-7));
}

TEST_CASE("AQL obeys the utilization law and stays under the capacity bound") {
    const std::vector<double> N{4.0, 3.0};
    const std::vector<double> Z{1.0, 2.0};
    auto a = pfqn_aql(demands(), N, Z);
    for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(a.UN(i, r) == doctest::Approx(a.XN[r] * demands()(i, r)).epsilon(1e-9));
        CHECK(a.XN[r] > 0.0);
        CHECK(a.XN[r] <= 1.0 / 0.4 + 1e-9);  // capacity of the slowest station
    }
}

TEST_CASE("AQL is closer to exact MVA than Bard-Schweitzer on this model") {
    const std::vector<double> N{4.0, 3.0};
    const std::vector<double> Z{0.0, 0.0};
    auto exact = pfqn_mva(demands(), std::vector<int>{4, 3});
    auto aql = pfqn_aql(demands(), N, Z);
    auto bs = pfqn_bs(demands(), N, Z);
    for (std::size_t r = 0; r < 2; ++r) {
        const double ea = std::fabs(aql.XN[r] - exact.XN[r]) / exact.XN[r];
        const double eb = std::fabs(bs.XN[r] - exact.XN[r]) / exact.XN[r];
        INFO("class ", r, ": AQL ", ea, " vs BS ", eb);
        CHECK(ea <= eb + 1e-12);
        CHECK(ea < 0.05);
    }
}

TEST_CASE("AQL degenerates gracefully on a single job") {
    Matrix<double> Ls(2, 1);
    Ls(0, 0) = 0.5;
    Ls(1, 0) = 0.4;
    auto a = pfqn_aql(Ls, std::vector<double>{1.0}, std::vector<double>{0.0});
    auto e = pfqn_mva(Ls, std::vector<int>{1});
    CHECK(a.XN[0] == doctest::Approx(e.XN[0]).epsilon(1e-6));
}
