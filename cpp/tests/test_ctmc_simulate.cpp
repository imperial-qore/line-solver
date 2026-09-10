/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * ctmc_simulate: the sample-path simulator of a CTMC given its generator.
 *
 * ORACLES. A simulator cannot be checked path-for-path against MATLAB (a
 * different generator and a different bits-to-deviates mapping), so every
 * check here is distributional or structural:
 *
 *  (a) the ergodic theorem. The time-average occupancy of a long path must
 *      converge to the exact stationary law, which ctmc_solve gives in closed
 *      form. This is the check that ties the holding times and the jump chain
 *      together: a wrong holding-time mean or a wrong jump matrix moves it.
 *  (b) the embedded jump chain. The empirical transition frequencies must
 *      converge to Q(i,j)/sum_{k!=i} Q(i,k), checked per row.
 *  (c) the holding-time law. The mean sojourn in state i must converge to
 *      -1/Q(i,i), and its second moment to 2/Q(i,i)^2 (exponential).
 *  (d) the initial distribution. The first state of many independent
 *      one-step runs must follow pi0. This is the assertion the reference
 *      fails, and it is stated against pi0 itself rather than against
 *      MATLAB's realized law.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_simulate.h"
#include "line/api/mc/ctmc_solve.h"

using line::Matrix;
using namespace line::mc;

namespace {

/** Two-state generator with stationary law [2/3, 1/3]. */
Matrix<double> two_state() {
    Matrix<double> Q(2, 2);
    Q(0, 0) = -1.0;
    Q(0, 1) = 1.0;
    Q(1, 0) = 2.0;
    Q(1, 1) = -2.0;
    return Q;
}

/** An irreducible three-state generator with distinct rates. */
Matrix<double> three_state() {
    Matrix<double> Q(3, 3);
    Q(0, 0) = -3.0;
    Q(0, 1) = 2.0;
    Q(0, 2) = 1.0;
    Q(1, 0) = 1.0;
    Q(1, 1) = -5.0;
    Q(1, 2) = 4.0;
    Q(2, 0) = 3.0;
    Q(2, 1) = 1.0;
    Q(2, 2) = -4.0;
    return Q;
}

}  // namespace

TEST_CASE("the simulated time-average occupancy converges to the stationary law") {
    // Oracle (a): the ergodic theorem against the exact ctmc_solve answer.
    const Matrix<double> Q = three_state();
    const std::vector<double> exact = ctmc_solve(Q);

    line::pfqn::McRng rng(20260721u);
    std::vector<double> pi0(3, 0.0);
    pi0[0] = 1.0;
    const CtmcPath<double> path = ctmc_simulate(Q, pi0, std::size_t(400000), rng);
    REQUIRE(path.states.size() == 400000u);
    REQUIRE(path.sojourn.size() == 400000u);

    std::vector<double> occ(3, 0.0);
    double total = 0.0;
    for (std::size_t k = 0; k < path.states.size(); ++k) {
        occ[path.states[k]] += path.sojourn[k];
        total += path.sojourn[k];
    }
    CHECK(total > 0.0);
    for (std::size_t i = 0; i < 3; ++i) CHECK(occ[i] / total == doctest::Approx(exact[i]).epsilon(0.01));

    // Oracle (c): the holding times are exponential of mean -1/Q(i,i), so the
    // mean and the second moment both pin down the same rate.
    std::vector<double> s1(3, 0.0), s2(3, 0.0);
    std::vector<std::size_t> cnt(3, 0);
    for (std::size_t k = 0; k < path.states.size(); ++k) {
        const std::size_t i = path.states[k];
        s1[i] += path.sojourn[k];
        s2[i] += path.sojourn[k] * path.sojourn[k];
        ++cnt[i];
    }
    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE(cnt[i] > 1000u);
        const double mean = -1.0 / Q(i, i);
        CHECK(s1[i] / static_cast<double>(cnt[i]) == doctest::Approx(mean).epsilon(0.02));
        CHECK(s2[i] / static_cast<double>(cnt[i]) ==
              doctest::Approx(2.0 * mean * mean).epsilon(0.05));
    }
}

TEST_CASE("the embedded jump chain matches Q(i,j) normalized off the diagonal") {
    // Oracle (b).
    const Matrix<double> Q = three_state();
    line::pfqn::McRng rng(4242u);
    const CtmcPath<double> path = ctmc_simulate(Q, std::vector<double>(), std::size_t(300000), rng);

    Matrix<double> cnt(3, 3, 0.0);
    std::vector<double> row(3, 0.0);
    for (std::size_t k = 0; k + 1 < path.states.size(); ++k) {
        cnt(path.states[k], path.states[k + 1]) += 1.0;
        row[path.states[k]] += 1.0;
    }
    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE(row[i] > 1000.0);
        double off = 0.0;
        for (std::size_t j = 0; j < 3; ++j)
            if (j != i) off += Q(i, j);
        for (std::size_t j = 0; j < 3; ++j) {
            const double want = (j == i) ? 0.0 : Q(i, j) / off;
            CHECK(cnt(i, j) / row[i] == doctest::Approx(want).epsilon(0.03));
        }
        // A CTMC never jumps to itself in the embedded chain.
        CHECK(cnt(i, i) == 0.0);
    }
}

TEST_CASE("the initial state is drawn from pi0, which the reference does not do") {
    // Oracle (d), and the regression for the reference defect. Each run is one
    // step, so the recorded state is exactly the initial draw.
    const Matrix<double> Q = three_state();
    const std::size_t runs = 200000;

    SUBCASE("uniform pi0") {
        std::vector<double> pi0(3, 1.0 / 3.0);
        line::pfqn::McRng rng(11u);
        std::vector<double> cnt(3, 0.0);
        for (std::size_t t = 0; t < runs; ++t) {
            const CtmcPath<double> p = ctmc_simulate(Q, pi0, std::size_t(1), rng);
            cnt[p.states[0]] += 1.0;
        }
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(cnt[i] / static_cast<double>(runs) == doctest::Approx(1.0 / 3.0).epsilon(0.02));

        // What MATLAB realizes instead: [1/2, 1/3, 1/6]. Measured over 400000
        // draws in the reference: [0.4997 0.3335 0.1668].
        const std::vector<double> ref = ctmc_simulate_reference_initial_law(pi0);
        CHECK(ref[0] == doctest::Approx(0.5));
        CHECK(ref[1] == doctest::Approx(1.0 / 3.0));
        CHECK(ref[2] == doctest::Approx(1.0 / 6.0));
        // The last state loses half its mass, which is the signature of the bug.
        CHECK(ref[2] == doctest::Approx(pi0[2] / 2.0));
        // and the law is still a law
        CHECK(ref[0] + ref[1] + ref[2] == doctest::Approx(1.0));
    }

    SUBCASE("a state of probability zero is never selected") {
        std::vector<double> pi0(3, 0.0);
        pi0[0] = 0.9;
        pi0[2] = 0.1;
        line::pfqn::McRng rng(12u);
        std::vector<double> cnt(3, 0.0);
        for (std::size_t t = 0; t < runs; ++t) {
            const CtmcPath<double> p = ctmc_simulate(Q, pi0, std::size_t(1), rng);
            cnt[p.states[0]] += 1.0;
        }
        CHECK(cnt[1] == 0.0);
        CHECK(cnt[0] / static_cast<double>(runs) == doctest::Approx(0.9).epsilon(0.01));
        CHECK(cnt[2] / static_cast<double>(runs) == doctest::Approx(0.1).epsilon(0.03));

        // MATLAB on the same pi0 gives [0.9497 0.0000 0.0503] over 400000
        // draws, and [0.9488 0.0000 0.0512] through the entry point itself.
        // The impossible state stays impossible only because it ties with
        // state 0 on the cumulative value (c = [0.9, 0.9, 1]) and min() takes
        // the first minimizer; state 0 therefore absorbs the whole window up
        // to the midpoint 0.95 and the LAST state keeps only half its mass.
        // Note this is NOT the (p_k + p_{k+1})/2 rule that holds when the
        // cumulative values are distinct -- a naive reading predicts 0.05 for
        // state 1 here, and the reference gives 0.
        const std::vector<double> ref = ctmc_simulate_reference_initial_law(pi0);
        CHECK(ref[0] == doctest::Approx(0.95));
        CHECK(ref[1] == doctest::Approx(0.0));
        CHECK(ref[2] == doctest::Approx(0.05));
    }

    SUBCASE("pi0 need not be normalized") {
        std::vector<double> pi0(2, 0.0);
        pi0[0] = 3.0;  // 3:1, i.e. [0.75, 0.25]
        pi0[1] = 1.0;
        line::pfqn::McRng rng(13u);
        std::vector<double> cnt(2, 0.0);
        for (std::size_t t = 0; t < 100000; ++t) {
            const CtmcPath<double> p = ctmc_simulate(two_state(), pi0, std::size_t(1), rng);
            cnt[p.states[0]] += 1.0;
        }
        CHECK(cnt[0] / 100000.0 == doctest::Approx(0.75).epsilon(0.02));
    }
}

TEST_CASE("an empty pi0 draws a random initial distribution") {
    // The reference's `r = rand(n,1); pi0 = r/sum(r)` branch. Over many runs
    // every state must be reachable as a starting point.
    const Matrix<double> Q = three_state();
    line::pfqn::McRng rng(77u);
    std::vector<double> cnt(3, 0.0);
    for (std::size_t t = 0; t < 20000; ++t) {
        const CtmcPath<double> p = ctmc_simulate(Q, std::vector<double>(), std::size_t(1), rng);
        cnt[p.states[0]] += 1.0;
    }
    for (std::size_t i = 0; i < 3; ++i) CHECK(cnt[i] > 1000.0);
}

TEST_CASE("ctmc_simulate rejects what it cannot simulate") {
    Matrix<double> Q = two_state();
    line::pfqn::McRng rng(1u);

    SUBCASE("a non-square generator") {
        Matrix<double> bad(2, 3, 0.0);
        CHECK_THROWS_AS(ctmc_simulate(bad, std::vector<double>(), std::size_t(1), rng),
                        line::InputError);
    }
    SUBCASE("pi0 of the wrong length") {
        CHECK_THROWS_AS(ctmc_simulate(Q, std::vector<double>(3, 1.0), std::size_t(1), rng),
                        line::InputError);
    }
    SUBCASE("pi0 with no mass") {
        CHECK_THROWS_AS(ctmc_simulate(Q, std::vector<double>(2, 0.0), std::size_t(1), rng),
                        line::InputError);
    }
    SUBCASE("an absorbing state") {
        // The reference silently jumps to state 1 with an infinite holding
        // time here; the port refuses and names the state.
        Matrix<double> abs2(2, 2, 0.0);
        abs2(0, 0) = -1.0;
        abs2(0, 1) = 1.0;
        CHECK_THROWS_AS(ctmc_simulate(abs2, std::vector<double>(), std::size_t(1), rng),
                        line::NumericError);
    }
    SUBCASE("zero steps is legal and returns an empty path") {
        const CtmcPath<double> p = ctmc_simulate(Q, std::vector<double>(), std::size_t(0), rng);
        CHECK(p.states.empty());
        CHECK(p.sojourn.empty());
    }
}

TEST_CASE("ctmc_simulate instantiates at Real50") {
    // The generator is double-backed (mt19937_64) but the holding times are
    // computed in the working arithmetic, so the path is representable at any
    // precision that has a logarithm.
    using R50 = line::Real50;
    Matrix<R50> Q(2, 2);
    Q(0, 0) = line::num_traits<R50>::from_int(-1);
    Q(0, 1) = line::num_traits<R50>::from_int(1);
    Q(1, 0) = line::num_traits<R50>::from_int(2);
    Q(1, 1) = line::num_traits<R50>::from_int(-2);

    line::pfqn::McRng rng(5u);
    const CtmcPath<R50> p = ctmc_simulate(Q, std::vector<R50>(), std::size_t(50000), rng);
    REQUIRE(p.states.size() == 50000u);
    R50 occ0 = line::num_traits<R50>::from_int(0);
    R50 tot = line::num_traits<R50>::from_int(0);
    for (std::size_t k = 0; k < p.states.size(); ++k) {
        if (p.states[k] == 0) occ0 += p.sojourn[k];
        tot += p.sojourn[k];
    }
    CHECK(line::num_traits<R50>::to_double(R50(occ0 / tot)) ==
          doctest::Approx(2.0 / 3.0).epsilon(0.02));
    for (std::size_t k = 0; k < p.sojourn.size(); ++k)
        CHECK(line::num_traits<R50>::to_double(p.sojourn[k]) > 0.0);
}

TEST_CASE("the reference initial-state law is characterized exactly") {
    // ctmc_simulate_reference_initial_law is the claim "this is what MATLAB
    // actually samples". It is checked here against the frequencies measured
    // in MATLAB over 400000 draws of
    //     [~,st] = min(abs(rand-cumsum(pi0)));
    // on seven pi0 vectors, covering distinct cumulative values, repeated
    // cumulative values (a zero-probability entry), an interior run of zeros,
    // and a degenerate point mass. Agreement across all of them is what makes
    // the defect report a characterization rather than an anecdote.
    struct Case {
        std::vector<double> pi0;
        std::vector<double> measured;
    };
    std::vector<Case> cases;
    cases.push_back({{0.5, 0.5}, {0.7498, 0.2502}});
    cases.push_back({{1.0 / 3, 1.0 / 3, 1.0 / 3}, {0.4997, 0.3335, 0.1668}});
    cases.push_back({{0.9, 0.0, 0.1}, {0.9497, 0.0000, 0.0503}});
    cases.push_back({{0.1, 0.2, 0.3, 0.4}, {0.2003, 0.2509, 0.3502, 0.1986}});
    cases.push_back({{0.4, 0.0, 0.0, 0.6}, {0.6998, 0.0, 0.0, 0.3002}});
    cases.push_back({{0.25, 0.25, 0.25, 0.25}, {0.3753, 0.2498, 0.2489, 0.1260}});
    cases.push_back({{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});

    for (std::size_t ci = 0; ci < cases.size(); ++ci) {
        INFO("case ", ci);
        const std::vector<double> got = ctmc_simulate_reference_initial_law(cases[ci].pi0);
        REQUIRE(got.size() == cases[ci].pi0.size());
        double sum = 0.0;
        for (std::size_t i = 0; i < got.size(); ++i) {
            // 400000 draws, so the sampling error is about 1e-3.
            CHECK(got[i] == doctest::Approx(cases[ci].measured[i]).epsilon(0.01).scale(1.0));
            sum += got[i];
            CHECK(got[i] >= 0.0);
        }
        CHECK(sum == doctest::Approx(1.0));  // still a probability law
        // A zero-probability state must never be selected, which the exact
        // tie rule delivers and a midpoint reading would not.
        for (std::size_t i = 0; i < got.size(); ++i)
            if (cases[ci].pi0[i] == 0.0) CHECK(got[i] == doctest::Approx(0.0));
    }
}
