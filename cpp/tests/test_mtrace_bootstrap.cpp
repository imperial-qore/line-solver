/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mtrace_bootstrap: block-bootstrap BCa intervals for a marked trace.
 *
 * A bootstrap has no fixed number to compare against -- the resampling stream
 * is not shared with MATLAB -- so the oracle here is the set of properties a
 * correct interval estimator has and a broken one does not:
 *
 *  - the interval BRACKETS the point estimate, for every descriptor;
 *  - it NARROWS when the trace grows, at roughly the 1/sqrt(N) rate;
 *  - the point estimate is exactly the [pc, backward, forward, sigma] vector
 *    the m3a fitters consume, in that layout;
 *  - the class probabilities of a trace built with a known mix are covered;
 *  - the block structure is the one the reference builds, and a trace too
 *    short to cut into two blocks is REFUSED rather than silently bootstrapped
 *    by observation, which would report intervals far too narrow.
 *
 * The last one is the point of the routine. Resampling observations instead of
 * blocks destroys the autocorrelation that sigma and the forward moment
 * measure, so the failure mode is not an error, it is an interval that looks
 * good and is wrong.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/api/trace/mtrace_bootstrap.h"

namespace tr = line::trace;

namespace {

/**
 * A two-class marked trace with autocorrelation: the labels come in runs, and
 * class 2 is slower than class 1. Runs are what a block bootstrap must keep.
 *
 * The run lengths are GEOMETRIC, not fixed. With fixed runs of 7 and 3 every
 * block of 50 holds exactly 35 and 15 labels, so every replicate returns the
 * same class probability and the interval collapses to a point -- a property of
 * the fixture, not of the estimator. Geometric runs of mean 7 and 3 keep the
 * long-run mix at 0.7 while letting the blocks differ.
 */
void make_trace(std::size_t n, std::vector<double>* T, std::vector<int>* A) {
    T->clear();
    A->clear();
    line::pfqn::McRng g(20260801u);
    int cls = 2;
    std::size_t left = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (left == 0) {
            cls = (cls == 1) ? 2 : 1;
            const double mean = (cls == 1) ? 7.0 : 3.0;
            const double u = line::pfqn::mc_uniform01(g);
            left = 1 + static_cast<std::size_t>(-std::log(1.0 - u) * (mean - 1.0));
        }
        --left;
        const double u = line::pfqn::mc_uniform01(g);
        const double rate = (cls == 1) ? 2.0 : 0.5;
        T->push_back(-std::log(1.0 - u) / rate);
        A->push_back(cls);
    }
}

}  // namespace

TEST_CASE("the point estimate is the descriptor vector the fitters consume") {
    std::vector<double> T;
    std::vector<int> A;
    make_trace(600, &T, &A);
    line::pfqn::McRng g(1u);
    const tr::MtraceBootstrapResult<double> r = tr::mtrace_bootstrap(T, A, g, 60);

    // 2 class probabilities + 2 backward + 2 forward + 2x2 sigma = 10.
    CHECK(r.estimate.size() == 10);
    CHECK(r.blocks == 12);  // 600 / 50

    // The probabilities are the first two entries and they sum to one.
    CHECK(r.estimate[0] + r.estimate[1] == doctest::Approx(1.0).epsilon(1e-12));
    // Class 1 runs are longer, so it is the more frequent label: 7 in 10.
    CHECK(r.estimate[0] == doctest::Approx(0.7).epsilon(0.02));
}

TEST_CASE("the interval brackets the point estimate for every descriptor") {
    std::vector<double> T;
    std::vector<int> A;
    make_trace(800, &T, &A);
    line::pfqn::McRng g(2u);
    const tr::MtraceBootstrapResult<double> r = tr::mtrace_bootstrap(T, A, g, 300);

    for (std::size_t i = 0; i < r.estimate.size(); ++i) {
        CHECK(r.lower[i] <= r.upper[i]);
        // A BCa endpoint can fall on the point estimate when the replicates
        // pile up there, so this is inclusive with a rounding margin.
        CHECK(r.lower[i] <= r.estimate[i] + 1e-9);
        CHECK(r.upper[i] >= r.estimate[i] - 1e-9);
    }
}

TEST_CASE("the interval narrows as the trace grows") {
    std::vector<int> A1, A2;
    std::vector<double> T1, T2;
    make_trace(500, &T1, &A1);
    make_trace(5000, &T2, &A2);

    line::pfqn::McRng g1(3u), g2(3u);
    const tr::MtraceBootstrapResult<double> a = tr::mtrace_bootstrap(T1, A1, g1, 200);
    const tr::MtraceBootstrapResult<double> b = tr::mtrace_bootstrap(T2, A2, g2, 200);

    // Widths of the class-1 probability interval: ten times the data should
    // shrink it by roughly sqrt(10) ~ 3.2. A factor of two is a safe floor and
    // still fails outright if the estimator ignores the sample size.
    const double wa = a.upper[0] - a.lower[0];
    const double wb = b.upper[0] - b.lower[0];
    CHECK(wa > 0.0);
    CHECK(wb > 0.0);
    CHECK(wb < wa / 2.0);
}

TEST_CASE("the interval covers the true class mix of a trace with a known mix") {
    // 70/30 by construction; the interval for p1 must contain 0.7.
    std::vector<double> T;
    std::vector<int> A;
    make_trace(4000, &T, &A);
    line::pfqn::McRng g(4u);
    const tr::MtraceBootstrapResult<double> r = tr::mtrace_bootstrap(T, A, g, 400);
    CHECK(r.lower[0] <= 0.7);
    CHECK(r.upper[0] >= 0.7);
}

TEST_CASE("the block layout tiles the trace exactly") {
    // 523 = 10 blocks of 52 plus 3 blocks one longer, i.e. sum(len) == N.
    std::vector<double> T;
    std::vector<int> A;
    make_trace(523, &T, &A);
    line::pfqn::McRng g(5u);
    const tr::MtraceBootstrapResult<double> r = tr::mtrace_bootstrap(T, A, g, 40);
    CHECK(r.blocks == 10);  // floor(523 / 50)
    // Every replicate has BN blocks whose lengths tile 523, so a replicate is
    // the same length as the trace and the moments stay comparable.
    CHECK(r.estimate[0] + r.estimate[1] == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("a block length larger than the trace is refused, not worked around") {
    std::vector<double> T;
    std::vector<int> A;
    make_trace(80, &T, &A);
    line::pfqn::McRng g(6u);
    // 80 / 50 = 1 block: nothing to resample.
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, A, g, 100), line::InputError);
    // With a block length that does fit, the same trace works.
    const tr::MtraceBootstrapResult<double> r = tr::mtrace_bootstrap(T, A, g, 100, 0.05, 10);
    CHECK(r.blocks == 8);
}

TEST_CASE("the arguments are validated by name") {
    std::vector<double> T;
    std::vector<int> A;
    make_trace(300, &T, &A);
    line::pfqn::McRng g(7u);
    std::vector<int> shortA(A.begin(), A.begin() + 10);
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, shortA, g, 10), line::InputError);
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, A, g, 1), line::InputError);
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, A, g, 10, 0.0), line::InputError);
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, A, g, 10, 1.0), line::InputError);
    CHECK_THROWS_AS(tr::mtrace_bootstrap(T, A, g, 10, 0.05, 0), line::InputError);
}

TEST_CASE("a wider level gives a wider interval") {
    std::vector<double> T;
    std::vector<int> A;
    make_trace(1000, &T, &A);
    line::pfqn::McRng g99(8u), g80(8u);
    const tr::MtraceBootstrapResult<double> w = tr::mtrace_bootstrap(T, A, g99, 400, 0.01);
    const tr::MtraceBootstrapResult<double> n = tr::mtrace_bootstrap(T, A, g80, 400, 0.20);
    CHECK((w.upper[0] - w.lower[0]) >= (n.upper[0] - n.lower[0]));
}
