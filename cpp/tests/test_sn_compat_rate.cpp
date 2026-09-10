/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The activated-server rate of a compatibility structure, and the scaling
 * SolverLN carries onto the layer station.
 *
 * Pool layout, the 3x2 of Dorsman and Gardner Fig. 1 reduced to two operands:
 * s1 serves op0 only, s2 serves both, s3 serves op1 only.
 */
#include <algorithm>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_compat_rate.h"

namespace api = line::api;
using line::Matrix;

namespace {

Matrix<double> compat() {
    Matrix<double> c(3, 2, 0.0);
    c(0, 0) = 1;
    c(1, 0) = 1;
    c(1, 1) = 1;
    c(2, 1) = 1;
    return c;
}

const std::vector<double> COUNTS = {1, 1, 1};
const std::vector<double> RATES = {1, 1, 1};

}  // namespace

TEST_CASE("sn_compat_rate: activated-server rate at the integer states") {
    CHECK(api::sn_compat_rate(compat(), COUNTS, RATES, {2, 0}) == doctest::Approx(2.0));
    CHECK(api::sn_compat_rate(compat(), COUNTS, RATES, {2, 1}) == doctest::Approx(3.0));
    CHECK(api::sn_compat_rate(compat(), COUNTS, RATES, {0, 1}) == doctest::Approx(2.0));
    CHECK(api::sn_compat_rate(compat(), COUNTS, RATES, {0, 0}) == doctest::Approx(0.0));
    CHECK(api::sn_compat_peak(COUNTS, RATES) == doctest::Approx(3.0));
}

TEST_CASE("sn_compat_rate: min(1,.) equals the hard indicator on the lattice") {
    // The invariant that keeps the order-independent law intact. The fractional
    // relaxation exists only so a mean-value solver can see the structure; CTMC
    // and the simulators evaluate at integer states alone and must not be able
    // to tell the two laws apart.
    const Matrix<double> c = compat();
    for (int n0 = 0; n0 < 4; ++n0) {
        for (int n1 = 0; n1 < 4; ++n1) {
            const std::vector<double> n = {double(n0), double(n1)};
            double hard = 0.0;
            for (std::size_t t = 0; t < 3; ++t)
                for (std::size_t j = 0; j < 2; ++j)
                    if (n[j] > 0 && c(t, j) != 0) {
                        hard += COUNTS[t] * RATES[t];
                        break;
                    }
            CHECK(api::sn_compat_rate(c, COUNTS, RATES, n) == doctest::Approx(hard));
        }
    }
}

TEST_CASE("sn_compat_rate: a pool scales below one compatible job") {
    // only op0 present at half a job: pools s1 and s2 reach it, each at 0.5
    CHECK(api::sn_compat_rate(compat(), COUNTS, RATES, {0.5, 0.0}) == doctest::Approx(1.0));
    // a pool sees the COMBINED load of the operands it can reach
    Matrix<double> both(1, 2, 1.0);
    const std::vector<double> one = {1};
    CHECK(api::sn_compat_rate(both, one, one, {0.4, 0.7}) == doctest::Approx(1.0));
    CHECK(api::sn_compat_rate(both, one, one, {0.4, 0.2}) == doctest::Approx(0.6));
}

TEST_CASE("sn_compat_scaling: a fully-compatible pool is neutral everywhere") {
    // What makes the lowering safe: the compatibility graph is the only thing
    // eta expresses. The statement of it is the EXACTNESS IDENTITY
    // min(N,S) * (peak/S) * eta(n) == mu(n) -- the solver applies min(N,S)
    // servers at the average server rate and eta corrects the product to the
    // activated-server rate, so a pool whose graph takes nothing away is one
    // where eta carries the whole remainder and nothing else.
    //
    // eta ITSELF IS ONE ONLY AT SATURATION, and above one below it: a pool of S
    // servers facing one compatible job clears S, not 1, because every one of
    // them works on it and the first to finish cancels the rest. That speed-up
    // is what the identity has to carry, and it is why the denominator damps by
    // min(1, N/S) rather than by min(1, N) -- the latter cancelled it, and
    // LDES, which simulates mu(n) directly, disagreed by that factor.
    Matrix<double> both(1, 2, 1.0);
    const std::vector<double> three = {3}, one = {1};
    const std::vector<std::vector<double>> states = {
        {2, 1}, {0.61, 0.90}, {0.2, 0.3}, {0, 1}, {5, 5}};
    const double peak = api::sn_compat_peak(three, one);
    for (std::size_t k = 0; k < states.size(); ++k) {
        const double total = states[k][0] + states[k][1];
        const double eta = api::sn_compat_scaling(both, three, one, states[k]);
        CHECK(std::min(total, three[0]) * (peak / three[0]) * eta ==
              doctest::Approx(api::sn_compat_rate(both, three, one, states[k])));
        if (total >= three[0])
            CHECK(eta == doctest::Approx(1.0));  // saturated: nothing left to carry
        else
            CHECK(eta > 1.0);                    // the redundancy speed-up
    }
}

TEST_CASE("sn_compat_scaling: a partial graph scales below neutral") {
    // the same three servers, with every operand reachable from every one
    Matrix<double> full(1, 2, 1.0);
    const std::vector<double> three = {3}, one = {1};
    const double whole = api::sn_compat_scaling(full, three, one, {0.61, 0.90});
    const double eta = api::sn_compat_scaling(compat(), COUNTS, RATES, {0.61, 0.90});
    CHECK(eta < whole);
    // mu = 0.61 + min(1, 1.51) + 0.90 = 2.51, over peak * min(1, 1.51/3)
    CHECK(eta == doctest::Approx(2.51 / 1.51));
    CHECK(whole == doctest::Approx(3.0 / 1.51));
}
