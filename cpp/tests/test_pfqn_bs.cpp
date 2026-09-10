/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Bard-Schweitzer AMVA. Oracles: the exact pfqn_mva solution (which BS must
 * reproduce at N = 1 and approach as the model becomes balanced), the
 * operational laws that any consistent solution obeys, and MATLAB reference
 * values.
 */
#include <vector>

#include "doctest.h"
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

TEST_CASE("Bard-Schweitzer is exact for one class with one job, not for several") {
    // Single class, N = 1: the (N-1)/N factor kills the only queueing term, so
    // BS coincides with exact MVA. MATLAB pfqn_bs and pfqn_mva both give
    // 1.111111111111 on this model.
    Matrix<double> Ls(2, 1);
    Ls(0, 0) = 0.5;
    Ls(1, 0) = 0.4;
    auto bs1 = pfqn_bs(Ls, std::vector<double>{1.0}, std::vector<double>{0.0});
    auto ex1 = pfqn_mva(Ls, std::vector<int>{1});
    CHECK(bs1.XN[0] == doctest::Approx(ex1.XN[0]).epsilon(1e-9));
    CHECK(bs1.XN[0] == doctest::Approx(1.111111111111).epsilon(1e-11));

    // Multiclass N = [1,1] is NOT exact: the cross-class terms still use the
    // queue lengths at the full population. MATLAB pfqn_bs returns
    // [0.749999990625, 0.740740740741] where exact MVA gives 0.75 for both,
    // i.e. a 1.2 percent error on the second class. The port reproduces the
    // MATLAB values to 12 digits.
    const std::vector<double> N{1.0, 1.0};
    const std::vector<double> Z{0.0, 0.0};
    auto bs = pfqn_bs(demands(), N, Z);
    auto exact = pfqn_mva(demands(), std::vector<int>{1, 1});
    CHECK(bs.converged);
    CHECK(bs.XN[0] == doctest::Approx(0.749999990625).epsilon(1e-11));
    CHECK(bs.XN[1] == doctest::Approx(0.740740740741).epsilon(1e-11));
    CHECK(exact.XN[0] == doctest::Approx(0.75).epsilon(1e-12));
    CHECK(std::fabs(bs.XN[1] - exact.XN[1]) / exact.XN[1] < 0.02);
}

TEST_CASE("Bard-Schweitzer obeys the operational laws") {
    const std::vector<double> N{4.0, 3.0};
    const std::vector<double> Z{1.0, 2.0};
    auto bs = pfqn_bs(demands(), N, Z);
    CHECK(bs.converged);

    for (std::size_t r = 0; r < 2; ++r) {
        // Little's law over the whole system, including the delay.
        double q = 0.0;
        for (std::size_t i = 0; i < 2; ++i) q += bs.QN(i, r);
        CHECK(q + bs.XN[r] * Z[r] == doctest::Approx(N[r]).epsilon(1e-6));
        // Utilization law.
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(bs.UN(i, r) == doctest::Approx(bs.XN[r] * demands()(i, r)).epsilon(1e-9));
        // Residence time consistency, R = Q/X.
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(bs.RN(i, r) == doctest::Approx(bs.QN(i, r) / bs.XN[r]).epsilon(1e-9));
    }
}

TEST_CASE("Bard-Schweitzer approaches the exact solution and stays bounded") {
    const std::vector<double> N{6.0, 4.0};
    const std::vector<double> Z{0.0, 0.0};
    auto bs = pfqn_bs(demands(), N, Z);
    auto exact = pfqn_mva(demands(), std::vector<int>{6, 4});
    for (std::size_t r = 0; r < 2; ++r) {
        // BS is known to be accurate to a few percent on models of this size.
        const double rel = std::fabs(bs.XN[r] - exact.XN[r]) / exact.XN[r];
        INFO("class ", r, " relative error ", rel);
        CHECK(rel < 0.10);
        // Never above the capacity bound.
        CHECK(bs.XN[r] <= 1.0 / 0.4 + 1e-9);
    }
}

TEST_CASE("FCFS and PS variants differ on heterogeneous demands") {
    const std::vector<double> N{4.0, 4.0};
    const std::vector<double> Z{0.0, 0.0};
    const std::vector<AmvaSched> fcfs{AmvaSched::FCFS, AmvaSched::FCFS};
    auto ps = pfqn_bs(demands(), N, Z);
    auto fc = pfqn_bs(demands(), N, Z, fcfs);
    CHECK(ps.converged);
    CHECK(fc.converged);
    CHECK(ps.XN[0] != doctest::Approx(fc.XN[0]).epsilon(1e-9));
}

TEST_CASE("an empty class contributes nothing and leaves no NaN") {
    const std::vector<double> N{3.0, 0.0};
    const std::vector<double> Z{0.0, 0.0};
    auto bs = pfqn_bs(demands(), N, Z);
    CHECK(bs.XN[1] == 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(bs.QN(i, 1) == 0.0);
        CHECK(bs.RN(i, 1) == 0.0);  // 0/0 would be NaN; the contract says 0
    }
    CHECK(bs.XN[0] > 0.0);
}
