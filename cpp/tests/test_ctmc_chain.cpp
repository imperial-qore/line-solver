/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverCTMC applied to a USER-SUPPLIED chain (`solver_ctmc_chain`,
 * `solver_ctmc_chain_transient`) and the single stationary entry point
 * (`ctmc_stationary`).
 *
 * THE ORACLE IS CLOSED FORM, not a second solver: a two-state chain has an
 * algebraic stationary vector and an algebraic transient, so a wrong
 * uniformization constant or a transposed balance equation shows up as a wrong
 * number rather than as a plausible one. The reducible cases are the point of
 * `ctmc_stationary` and are checked against absorption probabilities, which is
 * the one thing a weak-component split followed by renormalization gets wrong:
 * it spreads mass over states the declared start can never be absorbed in.
 */
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/ctmc_stationary.h"
#include "line/solvers/ctmc/solver_ctmc_chain.h"

using line::Matrix;
using line::ctmc::MarkovChainModel;

namespace {

/** Two-state CTMC: rate `a` out of state 0, rate `b` out of state 1. */
Matrix<double> two_state_gen(double a, double b) {
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -a;
    Q(0, 1) = a;
    Q(1, 0) = b;
    Q(1, 1) = -b;
    return Q;
}

}  // namespace

TEST_CASE("solver_ctmc_chain solves a MarkovProcess against the closed form") {
    const double a = 1.5, b = 0.5;
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(two_state_gen(a, b));
    const line::ctmc::CtmcChainSolution<double> s = line::ctmc::solver_ctmc_chain(m);

    REQUIRE(s.pi.size() == 2);
    CHECK(s.pi[0] == doctest::Approx(b / (a + b)).epsilon(1e-10));
    CHECK(s.pi[1] == doctest::Approx(a / (a + b)).epsilon(1e-10));
    // For a CTMC the returned infgen IS the generator, untouched.
    CHECK(s.infgen(0, 1) == doctest::Approx(a));
    // No state space was attached, so the state indices stand in, 1-based.
    REQUIRE(s.state_space.rows() == 2);
    CHECK(s.state_space(0, 0) == doctest::Approx(1.0));
    CHECK(s.state_space(1, 0) == doctest::Approx(2.0));
    CHECK(s.runtime >= 0.0);
}

TEST_CASE("solver_ctmc_chain solves a MarkovChain and returns P-I as the generator") {
    Matrix<double> P(2, 2, 0.0);
    P(0, 0) = 0.5;
    P(0, 1) = 0.5;
    P(1, 0) = 0.25;
    P(1, 1) = 0.75;
    const MarkovChainModel<double> m = MarkovChainModel<double>::chain(P);
    const line::ctmc::CtmcChainSolution<double> s = line::ctmc::solver_ctmc_chain(m);

    REQUIRE(s.pi.size() == 2);
    CHECK(s.pi[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-10));
    CHECK(s.pi[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-10));
    // P - I carries the SAME stationary vector, which is why one name serves
    // both chain kinds; the off-diagonal is untouched and the diagonal shifted.
    CHECK(s.infgen(0, 0) == doctest::Approx(-0.5));
    CHECK(s.infgen(0, 1) == doctest::Approx(0.5));
    CHECK(s.infgen(1, 1) == doctest::Approx(-0.25));
}

TEST_CASE("an attached state space is returned instead of the indices") {
    Matrix<double> space(2, 2, 0.0);
    space(0, 0) = 3.0;
    space(0, 1) = 0.0;
    space(1, 0) = 2.0;
    space(1, 1) = 1.0;
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(two_state_gen(1.0, 1.0), space);
    const line::ctmc::CtmcChainSolution<double> s = line::ctmc::solver_ctmc_chain(m);
    REQUIRE(s.state_space.rows() == 2);
    REQUIRE(s.state_space.cols() == 2);
    CHECK(s.state_space(0, 0) == doctest::Approx(3.0));
    CHECK(s.state_space(1, 1) == doctest::Approx(1.0));
}

TEST_CASE("solver_ctmc_chain_transient integrates the CTMC forward equations") {
    const double a = 2.0, b = 1.0;
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(two_state_gen(a, b));
    std::vector<double> pi0(2, 0.0);
    pi0[0] = 1.0;

    const line::ctmc::CtmcChainTransientSolution<double> r =
        line::ctmc::solver_ctmc_chain_transient(m, pi0, 0.0, 2.0);

    REQUIRE(r.t.size() >= 2);
    REQUIRE(r.pi_t.rows() == r.t.size());
    CHECK(r.t.front() == doctest::Approx(0.0));
    // Row 0 is the INITIAL CONDITION, not the first integrated point.
    CHECK(r.pi_t(0, 0) == doctest::Approx(1.0).epsilon(1e-10));

    // p0(t) = b/(a+b) + a/(a+b) exp(-(a+b) t) for the two-state chain.
    //
    // The tolerance is the INTEGRATOR'S, not a fudge: neither the reference nor
    // the port passes options to `ctmc_transient`, so this runs ode23 at its
    // default rtol = 1e-3 and a tighter assertion would be testing the step
    // controller rather than the forward equations. Mass conservation is
    // checked separately and IS sharp, because it is a structural invariant of
    // the integrator rather than an accuracy claim.
    for (std::size_t k = 0; k < r.t.size(); ++k) {
        const double t = r.t[k];
        const double exact = b / (a + b) + a / (a + b) * std::exp(-(a + b) * t);
        CHECK(r.pi_t(k, 0) == doctest::Approx(exact).epsilon(2e-3));
        CHECK(r.pi_t(k, 0) + r.pi_t(k, 1) == doctest::Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("solver_ctmc_chain_transient advances a DTMC one step per unit of time") {
    Matrix<double> P(2, 2, 0.0);
    P(0, 0) = 0.5;
    P(0, 1) = 0.5;
    P(1, 0) = 0.25;
    P(1, 1) = 0.75;
    const MarkovChainModel<double> m = MarkovChainModel<double>::chain(P);
    std::vector<double> pi0(2, 0.0);
    pi0[0] = 1.0;

    const line::ctmc::CtmcChainTransientSolution<double> r =
        line::ctmc::solver_ctmc_chain_transient(m, pi0, 0.0, 3.0);

    // The timespan holds the integer steps 0..3, so four rows and not a grid.
    REQUIRE(r.t.size() == 4);
    for (std::size_t k = 0; k < 4; ++k) CHECK(r.t[k] == doctest::Approx(double(k)));
    CHECK(r.pi_t(0, 0) == doctest::Approx(1.0));
    CHECK(r.pi_t(1, 0) == doctest::Approx(0.5));
    CHECK(r.pi_t(1, 1) == doctest::Approx(0.5));
    // step 2: [0.5,0.5] P = [0.375, 0.625]
    CHECK(r.pi_t(2, 0) == doctest::Approx(0.375));
    CHECK(r.pi_t(2, 1) == doctest::Approx(0.625));
    // step 3: [0.375,0.625] P = [0.34375, 0.65625]
    CHECK(r.pi_t(3, 0) == doctest::Approx(0.34375));

    // A non-integer lower end starts at the first step INSIDE the span.
    const line::ctmc::CtmcChainTransientSolution<double> r2 =
        line::ctmc::solver_ctmc_chain_transient(m, pi0, 1.5, 3.0);
    REQUIRE(r2.t.size() == 2);
    CHECK(r2.t[0] == doctest::Approx(2.0));
    CHECK(r2.pi_t(0, 0) == doctest::Approx(0.375));
}

TEST_CASE("an empty initial distribution is the uniform one") {
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(two_state_gen(1.0, 1.0));
    const line::ctmc::CtmcChainTransientSolution<double> r =
        line::ctmc::solver_ctmc_chain_transient(m, std::vector<double>(), 0.0, 1.0);
    CHECK(r.pi_t(0, 0) == doctest::Approx(0.5));
    CHECK(r.pi_t(0, 1) == doctest::Approx(0.5));
}

TEST_CASE("the transient refuses what it cannot answer") {
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(two_state_gen(1.0, 1.0));
    std::vector<double> bad(2, 0.3);  // sums to 0.6
    CHECK_THROWS_AS(line::ctmc::solver_ctmc_chain_transient(m, bad, 0.0, 1.0), line::InputError);

    std::vector<double> pi0(2, 0.5);
    // An infinite horizon is not a transient question.
    CHECK_THROWS_AS(
        line::ctmc::solver_ctmc_chain_transient(m, pi0, 0.0, std::numeric_limits<double>::infinity()),
        line::InputError);

    const MarkovChainModel<double> d = MarkovChainModel<double>::chain(
        Matrix<double>{{0.5, 0.5}, {0.25, 0.75}});
    // [1.2, 1.8] contains no integer step, so there is nothing to report.
    CHECK_THROWS_AS(line::ctmc::solver_ctmc_chain_transient(d, pi0, 1.2, 1.8), line::InputError);
}

TEST_CASE("ctmc_stationary reproduces the closed form on an irreducible chain") {
    const double a = 3.0, b = 1.0;
    const line::ctmc::CtmcStationaryResult<double> s = line::ctmc::ctmc_stationary(two_state_gen(a, b));
    REQUIRE(s.pi.size() == 2);
    CHECK(s.pi[0] == doctest::Approx(b / (a + b)).epsilon(1e-10));
    CHECK(s.pi[1] == doctest::Approx(a / (a + b)).epsilon(1e-10));
    // One BSCC is the degenerate case, and it is NOT a mixture to warn about.
    CHECK(s.nbscc == 1);
    CHECK(s.warning.empty());
    CHECK_FALSE(s.seeded);
}

TEST_CASE("ctmc_stationary weights the BSCCs by absorption from the declared start") {
    // 0 --2--> 1 (absorbing), 0 --1--> 2 (absorbing). From state 0 the chain is
    // absorbed in 1 with probability 2/3 and in 2 with probability 1/3.
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 1) = 2.0;
    Q(0, 2) = 1.0;
    Q(0, 0) = -3.0;

    const line::ctmc::CtmcStationaryResult<double> s = line::ctmc::ctmc_stationary(Q, 0);
    REQUIRE(s.pi.size() == 3);
    CHECK(s.seeded);
    CHECK(s.nbscc == 2);
    // A seeded solve is the model's answer, so no warning is raised.
    CHECK(s.warning.empty());
    // State 0 is TRANSIENT and carries no stationary mass, whatever it was seeded with.
    CHECK(s.pi[0] == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(s.pi[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
    CHECK(s.pi[2] == doctest::Approx(1.0 / 3.0).epsilon(1e-6));
}

TEST_CASE("an unseeded reducible mixture is answered but flagged") {
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 1) = 2.0;
    Q(0, 2) = 1.0;
    Q(0, 0) = -3.0;

    const line::ctmc::CtmcStationaryResult<double> s = line::ctmc::ctmc_stationary(Q);
    CHECK_FALSE(s.seeded);
    CHECK(s.nbscc == 2);
    // The number is still produced -- the caller may have no state to declare --
    // but it is the invented start's answer and must say so.
    CHECK_FALSE(s.warning.empty());
    CHECK(s.warning.find("closed communicating classes") != std::string::npos);
    double total = 0.0;
    for (std::size_t i = 0; i < s.pi.size(); ++i) total += s.pi[i];
    CHECK(total == doctest::Approx(1.0).epsilon(1e-8));
}

TEST_CASE("ctmc_stationary refuses a malformed seed") {
    CHECK_THROWS_AS(line::ctmc::ctmc_stationary(two_state_gen(1.0, 1.0), 7), line::InputError);
    Matrix<double> notsquare(2, 3, 0.0);
    CHECK_THROWS_AS(line::ctmc::ctmc_stationary(notsquare), line::InputError);
}
