/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * wf_pattern_updater: the phase-type algebra behind the workflow collapse.
 * Ported from native Python on 2026-08-01, which stays the reference for
 * api/wf. The JAR's four convolutions returned the first branch's law until
 * they were ported on the same date; the same oracles now pin all three.
 *
 * ORACLES. Every convolution has a closed-form moment its result must match,
 * and none of them depends on the representation the formula happens to build:
 *  - sequence: the means add;
 *  - parallel: the max of two independent exponentials has mean
 *    1/l1 + 1/l2 - 1/(l1+l2);
 *  - loop: a geometric number of repetitions has mean m/(1-p);
 *  - branch: the means mix with the branch probabilities.
 * The mean of an (alpha, T) law is alpha (-T)^-1 e, computed here by a solve,
 * so the test never reuses the assembly it is checking.
 */
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "doctest.h"
#include "line/api/wf/wf_pattern_updater.h"
#include "line/util/lu.h"

namespace wf = line::wf;
using line::Matrix;

namespace {

/** An exponential of rate lambda as an (alpha, T) pair. */
wf::ServiceParameters<double> expo(double lambda) {
    wf::ServiceParameters<double> p;
    p.alpha.assign(1, 1.0);
    p.T_ = Matrix<double>(1, 1, -lambda);
    return p;
}

/** Erlang-k of total mean k/rate. */
wf::ServiceParameters<double> erlang(double rate, std::size_t k) {
    wf::ServiceParameters<double> p;
    p.alpha.assign(k, 0.0);
    p.alpha[0] = 1.0;
    p.T_ = Matrix<double>(k, k, 0.0);
    for (std::size_t i = 0; i < k; ++i) {
        p.T_(i, i) = -rate;
        if (i + 1 < k) p.T_(i, i + 1) = rate;
    }
    return p;
}

/** alpha (-T)^-1 e, the mean of the law, by a solve rather than by assembly. */
double ph_mean(const wf::ServiceParameters<double>& p) {
    const std::size_t n = p.alpha.size();
    Matrix<double> A(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = -p.T_(i, j);
    std::vector<double> e(n, 1.0);
    const std::vector<double> x = line::solve(A, e);
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i) m += p.alpha[i] * x[i];
    return m;
}

}  // namespace

TEST_CASE("convolve_sequence adds the means") {
    std::vector<wf::ServiceParameters<double>> ps;
    ps.push_back(expo(2.0));      // mean 0.5
    ps.push_back(erlang(4.0, 3)); // mean 0.75
    ps.push_back(expo(5.0));      // mean 0.2
    const wf::ServiceParameters<double> c = wf::convolve_sequence(ps);
    CHECK(c.alpha.size() == 5);
    CHECK(ph_mean(c) == doctest::Approx(0.5 + 0.75 + 0.2).epsilon(1e-10));

    // Degenerate inputs take the reference's fallbacks.
    CHECK(wf::convolve_sequence(std::vector<wf::ServiceParameters<double>>()).alpha.size() == 1);
    std::vector<wf::ServiceParameters<double>> one(1, expo(3.0));
    CHECK(ph_mean(wf::convolve_sequence(one)) == doctest::Approx(1.0 / 3.0));
}

TEST_CASE("convolve_parallel gives the maximum of two exponentials") {
    // E[max(X1,X2)] = 1/l1 + 1/l2 - 1/(l1+l2) for independent exponentials.
    const double l1 = 2.0, l2 = 3.0;
    std::vector<wf::ServiceParameters<double>> ps;
    ps.push_back(expo(l1));
    ps.push_back(expo(l2));
    const wf::ServiceParameters<double> c = wf::convolve_parallel(ps);
    CHECK(c.alpha.size() == 1 * 1 + 1 + 1);
    const double want = 1.0 / l1 + 1.0 / l2 - 1.0 / (l1 + l2);
    CHECK(ph_mean(c) == doctest::Approx(want).epsilon(1e-10));
}

TEST_CASE("convolve_parallel dominates each branch and is symmetric in the mean") {
    std::vector<wf::ServiceParameters<double>> ab, ba;
    ab.push_back(erlang(4.0, 2));  // mean 0.5
    ab.push_back(expo(1.0));       // mean 1.0
    ba.push_back(expo(1.0));
    ba.push_back(erlang(4.0, 2));
    const double m1 = ph_mean(wf::convolve_parallel(ab));
    const double m2 = ph_mean(wf::convolve_parallel(ba));
    CHECK(m1 == doctest::Approx(m2).epsilon(1e-10));
    CHECK(m1 > 1.0);  // the max exceeds the slower branch's own mean
}

TEST_CASE("convolve_loop inflates the mean by the geometric factor") {
    const double p = 0.4;
    const wf::ServiceParameters<double> c = wf::convolve_loop(erlang(4.0, 2), p);
    CHECK(ph_mean(c) == doctest::Approx(0.5 / (1.0 - p)).epsilon(1e-10));
    // Outside (0,1) the law is returned unchanged, as in the reference.
    CHECK(ph_mean(wf::convolve_loop(erlang(4.0, 2), 0.0)) == doctest::Approx(0.5));
    CHECK(ph_mean(wf::convolve_loop(erlang(4.0, 2), 1.0)) == doctest::Approx(0.5));
}

TEST_CASE("convolve_branches mixes the means by the branch probabilities") {
    std::vector<wf::ServiceParameters<double>> ps;
    ps.push_back(expo(2.0));       // mean 0.5
    ps.push_back(erlang(4.0, 4));  // mean 1.0
    std::vector<double> probs;
    probs.push_back(0.25);
    probs.push_back(0.75);
    const wf::ServiceParameters<double> c = wf::convolve_branches(ps, probs);
    CHECK(c.alpha.size() == 5);
    CHECK(ph_mean(c) == doctest::Approx(0.25 * 0.5 + 0.75 * 1.0).epsilon(1e-10));

    // Unnormalized probabilities are renormalized, so the answer is the same.
    std::vector<double> raw;
    raw.push_back(1.0);
    raw.push_back(3.0);
    CHECK(ph_mean(wf::convolve_branches(ps, raw)) == doctest::Approx(ph_mean(c)).epsilon(1e-12));

    // A zero total falls back to a uniform choice.
    std::vector<double> zeros(2, 0.0);
    CHECK(ph_mean(wf::convolve_branches(ps, zeros)) == doctest::Approx(0.75).epsilon(1e-12));
}

TEST_CASE("update_patterns collapses a chain and convolves its service laws") {
    // 1 -> 2 -> 3 -> 4 with 2,3,4 service nodes: the sequence 2-3-4 collapses
    // onto node 2, whose law becomes the convolution of the three.
    Matrix<double> L(3, 3, 0.0);
    L(0, 0) = 1; L(0, 1) = 2; L(0, 2) = 1.0;
    L(1, 0) = 2; L(1, 1) = 3; L(1, 2) = 1.0;
    L(2, 0) = 3; L(2, 1) = 4; L(2, 2) = 1.0;

    std::vector<int> svc;
    svc.push_back(2);
    svc.push_back(3);
    svc.push_back(4);
    std::map<int, wf::ServiceParameters<double>> params;
    params[2] = expo(2.0);  // 0.5
    params[3] = expo(4.0);  // 0.25
    params[4] = expo(5.0);  // 0.2

    const wf::UpdatedWorkflow<double> w = wf::update_patterns(
        L, svc, std::vector<int>(), std::vector<int>(), std::vector<int>(), params);

    CHECK(w.linkMatrix.rows() < L.rows());          // the chain's edges are gone
    REQUIRE(w.serviceParameters.count(2) == 1);
    CHECK(w.serviceParameters.count(3) == 0);       // absorbed
    CHECK(w.serviceParameters.count(4) == 0);
    CHECK(ph_mean(w.serviceParameters.at(2)) == doctest::Approx(0.95).epsilon(1e-10));

    const wf::UpdateStats<double> s = wf::get_update_stats(L, w);
    CHECK(s.originalLinks == 3);
    CHECK(s.linksReduced == static_cast<long>(3 - w.linkMatrix.rows()));
    CHECK(s.reductionRatio == doctest::Approx(static_cast<double>(s.linksReduced) / 3.0));
}

TEST_CASE("update_patterns is the identity when nothing is detected") {
    Matrix<double> L(1, 3, 0.0);
    L(0, 0) = 1; L(0, 1) = 2; L(0, 2) = 1.0;
    std::map<int, wf::ServiceParameters<double>> params;
    params[2] = expo(1.0);
    const wf::UpdatedWorkflow<double> w = wf::update_patterns(
        L, std::vector<int>(1, 2), std::vector<int>(), std::vector<int>(), std::vector<int>(),
        params);
    CHECK(w.linkMatrix.rows() == 1);
    CHECK(w.serviceParameters.size() == 1);
    CHECK(ph_mean(w.serviceParameters.at(2)) == doctest::Approx(1.0));
}
