/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * The MarkovProcess / MarkovChain object surface
 * (`cpp/include/line/lang/processes/markov_chain.h`).
 *
 * ORACLE: live MATLAB, run 2026-08-15 through the matlab MCP server on the
 * fixture below, with `lineStart` from `matlab/`. Every literal here was
 * printed by `%.15g` from the reference method of the same name, so a failure
 * is a divergence from LINE 3.0.7 and not from a hand computation. Where the
 * reference has a closed form the test asserts that too -- the fixture's
 * stationary vector is (12,4,5)/21 exactly, which is what makes the Cramer
 * determinants readable as integers.
 *
 * The fixture is a 3-state CTMC that is deliberately NOT reversible and has one
 * zero off-diagonal entry (state 3 cannot jump to state 2), so `toEmbedded` has
 * a row the uniformized `toDTMC` does not, and the two conversions cannot be
 * confused by accident:
 *
 *     Q = [-2  1  1;  1 -3  2;  4  0 -4]
 *
 * The aggregation and stochastic-complement arms are checked as DISPATCH: this
 * header chooses among `api/mc` primitives that `test_mc_aggregation.cpp`
 * already validates against MATLAB, so re-asserting their numerics here would
 * pin the same fact twice and drift apart. What is asserted is that the surface
 * reaches the primitive the reference's `switch` reaches, which is the part
 * this header is responsible for.
 */
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/lang/processes/markov_chain.h"
#include "line/util/lu.h"

using line::Matrix;
using line::lang::processes::MarkovChainModel;

namespace proc = line::lang::processes;

namespace {

/** The fixture generator, as MATLAB was given it. */
Matrix<double> fixture_q() {
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 0) = -2.0; Q(0, 1) = 1.0; Q(0, 2) = 1.0;
    Q(1, 0) = 1.0;  Q(1, 1) = -3.0; Q(1, 2) = 2.0;
    Q(2, 0) = 4.0;  Q(2, 1) = 0.0;  Q(2, 2) = -4.0;
    return Q;
}

MarkovChainModel<double> fixture_process() {
    Matrix<double> space(3, 1, 0.0);
    space(0, 0) = 1.0; space(1, 0) = 2.0; space(2, 0) = 3.0;
    return MarkovChainModel<double>::process(fixture_q(), space);
}

void check_vec(const std::vector<double>& got, const std::vector<double>& want, double tol) {
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
        CAPTURE(i);
        CHECK(std::fabs(got[i] - want[i]) <= tol);
    }
}

}  // namespace

TEST_CASE("chain_solve reproduces MarkovProcess.solve") {
    const std::vector<double> pi = proc::chain_solve(fixture_process());
    check_vec(pi, {0.571428571428571, 0.19047619047619, 0.238095238095238}, 1e-12);
    // The exact law is (12,4,5)/21, which is what makes the Cramer test below
    // read as integers; assert it independently of the printed decimals.
    CHECK(std::fabs(pi[0] - 12.0 / 21.0) <= 1e-14);
    CHECK(std::fabs(pi[1] - 4.0 / 21.0) <= 1e-14);
    CHECK(std::fabs(pi[2] - 5.0 / 21.0) <= 1e-14);
}

TEST_CASE("toEmbedded is the jump chain and does NOT preserve the stationary law") {
    const MarkovChainModel<double> emb = proc::to_embedded(fixture_process());
    REQUIRE(emb.discrete);
    // MATLAB MarkovProcess.toEmbedded on the fixture.
    const double want[3][3] = {{0.0, 0.5, 0.5},
                               {1.0 / 3.0, 0.0, 2.0 / 3.0},
                               {1.0, 0.0, 0.0}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            CAPTURE(i); CAPTURE(j);
            CHECK(std::fabs(emb.mat(i, j) - want[i][j]) <= 1e-14);
        }
    // THE POINT OF THE SEPARATE CONVERSION: the jump chain's stationary law is
    // (3,1.5,2.5)/7, not the CTMC's (12,4,5)/21. A port that answered with the
    // uniformized chain here would still be stochastic and still look right.
    check_vec(proc::chain_solve(emb),
              {0.428571428571429, 0.214285714285714, 0.357142857142857}, 1e-12);
    CHECK(std::fabs(proc::chain_solve(emb)[0] - proc::chain_solve(fixture_process())[0]) > 0.1);
    // The state space rides along, as setStateSpace does in the reference.
    REQUIRE(emb.state_space.rows() == 3);
}

TEST_CASE("toEmbedded leaves an absorbing state absorbing") {
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -1.0; Q(0, 1) = 1.0;  // state 1 absorbs
    const MarkovChainModel<double> emb = proc::to_embedded(MarkovChainModel<double>::process(Q));
    CHECK(std::fabs(emb.mat(0, 1) - 1.0) <= 1e-14);
    CHECK(std::fabs(emb.mat(1, 1) - 1.0) <= 1e-14);
    CHECK(std::fabs(emb.mat(1, 0)) <= 1e-14);
}

TEST_CASE("toDTMC uniformizes and DOES preserve the stationary law") {
    const MarkovChainModel<double> dt = proc::to_markov_chain(fixture_process(), 10.0);
    // MATLAB mp.toMarkovChain(10).
    const double want[3][3] = {{0.8, 0.1, 0.1}, {0.1, 0.7, 0.2}, {0.4, 0.0, 0.6}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            CAPTURE(i); CAPTURE(j);
            CHECK(std::fabs(dt.mat(i, j) - want[i][j]) <= 1e-14);
        }
    check_vec(proc::chain_solve(dt), proc::chain_solve(fixture_process()), 1e-10);

    // The deterministic default rate is strictly above max|Q|, which is what
    // makes the uniformization valid; the reference reaches the same guarantee
    // with a random draw it cannot reproduce.
    const double q = proc::default_uniformization_rate(fixture_q());
    CHECK(q > 4.0);
    check_vec(proc::chain_solve(proc::to_markov_chain(fixture_process())),
              proc::chain_solve(fixture_process()), 1e-10);
}

TEST_CASE("toMarkovProcess and toTimeReversed round-trip the two kinds") {
    const MarkovChainModel<double> emb = proc::to_embedded(fixture_process());
    const MarkovChainModel<double> back = proc::to_markov_process(emb);
    REQUIRE_FALSE(back.discrete);
    for (std::size_t i = 0; i < 3; ++i) CHECK(back.mat(i, i) < 0.0);

    // MATLAB mp.toTimeReversed().getGenerator().
    const MarkovChainModel<double> rev = proc::to_time_reversed(fixture_process());
    const double want[3][3] = {{-2.0, 1.0 / 3.0, 5.0 / 3.0},
                               {3.0, -3.0, 0.0},
                               {2.4, 1.6, -4.0}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            CAPTURE(i); CAPTURE(j);
            CHECK(std::fabs(rev.mat(i, j) - want[i][j]) <= 1e-12);
        }
    // A time reversal carries the same stationary law, by construction.
    check_vec(proc::chain_solve(rev), proc::chain_solve(fixture_process()), 1e-12);
}

TEST_CASE("getProbState is Cramer's rule and returns both determinants") {
    const MarkovChainModel<double> m = fixture_process();
    const double num[3] = {12.0, 4.0, 5.0};
    for (std::size_t s = 0; s < 3; ++s) {
        CAPTURE(s);
        const proc::ProbStateResult<double> r = proc::get_prob_state(m, s);
        CHECK(std::fabs(r.den - 21.0) <= 1e-11);
        CHECK(std::fabs(r.num - num[s]) <= 1e-11);
        CHECK(std::fabs(r.pi_i - num[s] / 21.0) <= 1e-13);
    }
    // Addressed by the state row rather than by the index.
    std::vector<double> st(1, 2.0);
    CHECK(std::fabs(proc::get_prob_state(m, st).pi_i - 4.0 / 21.0) <= 1e-13);
}

TEST_CASE("lu_det returns zero on a singular matrix rather than throwing") {
    Matrix<double> A(2, 2, 1.0);  // rank one
    CHECK(std::fabs(line::lu_det(A)) <= 1e-15);
    Matrix<double> B(2, 2, 0.0);
    B(0, 0) = 2.0; B(0, 1) = 1.0; B(1, 0) = 1.0; B(1, 1) = 3.0;
    CHECK(std::fabs(line::lu_det(B) - 5.0) <= 1e-14);
    // The sign is the parity of the swaps: this one needs exactly one.
    Matrix<double> C(2, 2, 0.0);
    C(0, 1) = 1.0; C(1, 0) = 1.0;
    CHECK(std::fabs(line::lu_det(C) + 1.0) <= 1e-14);
}

TEST_CASE("solveRelative pins the reference state at one") {
    const MarkovChainModel<double> m = fixture_process();
    check_vec(proc::solve_relative(m, 0), {1.0, 1.0 / 3.0, 5.0 / 12.0}, 1e-12);
    check_vec(proc::solve_relative(m, 1), {3.0, 1.0, 1.25}, 1e-12);
}

TEST_CASE("transient, timeAverage and their DTMC counterparts") {
    const MarkovChainModel<double> m = fixture_process();
    const std::vector<double> pi0 = {1.0, 0.0, 0.0};
    const std::vector<double> at_half = {0.608782910966201, 0.179548497347113,
                                         0.211668591686327};

    const proc::TransientAtResult<double> u = proc::chain_transient_at(m, pi0, 0.5);
    check_vec(u.pi, at_half, 1e-12);
    CHECK(u.kmax == 19u);

    // Fox-Glynn answers the same question by different weights: same law, its
    // own truncation point. Agreement to 1e-12 is the real assertion.
    const proc::TransientAtResult<double> f = proc::chain_transient_at(m, pi0, 0.5, "foxglynn");
    check_vec(f.pi, at_half, 1e-12);
    CHECK(f.kmax == 20u);

    const proc::TimeAverageOut<double> ta = proc::time_average(m, pi0, 0.5);
    check_vec(ta.pi_time_avg,
              {0.734504453251639, 0.125135819519092, 0.140359727229249}, 1e-12);
    check_vec(ta.pi_exit, at_half, 1e-12);

    // MarkovChain.transient counts STEPS and returns the whole trajectory.
    const MarkovChainModel<double> emb = proc::to_embedded(m);
    const Matrix<double> traj = proc::chain_transient_steps(emb, pi0, 3);
    REQUIRE(traj.rows() == 4);
    const double want[4][3] = {{1.0, 0.0, 0.0},
                               {0.0, 0.5, 0.5},
                               {2.0 / 3.0, 0.0, 1.0 / 3.0},
                               {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    for (std::size_t k = 0; k < 4; ++k)
        for (std::size_t j = 0; j < 3; ++j) {
            CAPTURE(k); CAPTURE(j);
            CHECK(std::fabs(traj(k, j) - want[k][j]) <= 1e-13);
        }

    // MarkovChain.transientUnif reads the same matrix on a CONTINUOUS clock.
    const proc::TransientAtResult<double> tu = proc::chain_transient_unif(emb, pi0, 0.5);
    check_vec(tu.pi, {0.66206490974292, 0.156146165018171, 0.181788925237946}, 1e-12);
    CHECK(tu.kmax == 11u);
}

TEST_CASE("sens differentiates the stationary law") {
    Matrix<double> dQ(3, 3, 0.0);
    dQ(0, 0) = -1.0; dQ(0, 1) = 1.0;
    check_vec(proc::chain_sens(fixture_process(), dQ),
              {-0.163265306122449, 0.136054421768707, 0.0272108843537415}, 1e-12);
}

TEST_CASE("hittingTime measures time for a CTMC and steps for a DTMC") {
    const MarkovChainModel<double> m = fixture_process();
    const std::vector<std::size_t> target(1, 2);
    check_vec(proc::hitting_time(m, target), {0.8, 0.6, 0.0}, 1e-12);
    // The jump chain counts steps, so the same target reads differently. This
    // is the pair that catches a port dispatching hitting_time on the wrong kind.
    check_vec(proc::hitting_time(proc::to_embedded(m), target), {1.8, 1.6, 0.0}, 1e-12);
}

TEST_CASE("stochComp complements a CTMC by its generator and a DTMC by its matrix") {
    const MarkovChainModel<double> m = fixture_process();
    const std::vector<std::size_t> I = {0, 1};
    const proc::StochCompOut<double> c = proc::stoch_comp_full(m, I);
    // MATLAB mp.stochCompFull([1 2]): S and the return-path term T.
    CHECK(std::fabs(c.S(0, 0) + 1.0) <= 1e-12);
    CHECK(std::fabs(c.S(0, 1) - 1.0) <= 1e-12);
    CHECK(std::fabs(c.S(1, 0) - 3.0) <= 1e-12);
    CHECK(std::fabs(c.S(1, 1) + 3.0) <= 1e-12);
    CHECK(std::fabs(c.T12(0, 0) - 1.0) <= 1e-12);
    CHECK(std::fabs(c.T12(1, 0) - 2.0) <= 1e-12);
    // S = A11 + T12, the identity the block form exists to expose.
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            CAPTURE(i); CAPTURE(j);
            CHECK(std::fabs(c.S(i, j) - (c.A11(i, j) + c.T12(i, j))) <= 1e-12);
        }

    // MATLAB E.stochComp([1 2]) on the embedded DTMC.
    const Matrix<double> Sd = proc::stoch_comp(proc::to_embedded(m), I);
    CHECK(std::fabs(Sd(0, 0) - 0.5) <= 1e-12);
    CHECK(std::fabs(Sd(0, 1) - 0.5) <= 1e-12);
    CHECK(std::fabs(Sd(1, 0) - 1.0) <= 1e-12);
    CHECK(std::fabs(Sd(1, 1)) <= 1e-12);
}

TEST_CASE("aggregate reaches the primitive the reference's switch reaches") {
    // A nearly completely decomposable chain: two blocks, coupling of order 0.01.
    Matrix<double> Q(4, 4, 0.0);
    Q(0, 1) = 1.0;  Q(0, 2) = 0.01;
    Q(1, 0) = 2.0;  Q(1, 3) = 0.02;
    Q(2, 0) = 0.03; Q(2, 3) = 3.0;
    Q(3, 1) = 0.04; Q(3, 2) = 4.0;
    const MarkovChainModel<double> m = MarkovChainModel<double>::process(Q);
    const std::vector<std::vector<std::size_t>> MS = {{0, 1}, {2, 3}};

    const proc::AggregateResult<double> c = proc::aggregate(m, MS);
    check_vec(c.p, line::mc::ctmc_courtois(m.mat, MS).p, 1e-14);
    CHECK(std::fabs(c.eps - line::mc::ctmc_courtois(m.mat, MS).eps) <= 1e-14);

    const double ten = 10.0;
    const proc::AggregateResult<double> k = proc::aggregate(m, MS, "kms", &ten);
    check_vec(k.p, line::mc::ctmc_kms(m.mat, MS, 10).p, 1e-14);

    const proc::AggregateResult<double> t = proc::aggregate(m, MS, "takahashi", &ten);
    check_vec(t.p, line::mc::ctmc_takahashi(m.mat, MS, 10).p, 1e-14);

    // EVERY ARM MUST LAND ON THE EXACT LAW, which is what says the dispatch
    // reached an aggregation method at all: a surface that returned the uniform
    // vector, or the wrong block's estimate, would still be a distribution and
    // would still have the right length. This partition is decomposable enough
    // that all three are at machine precision, so no arm is allowed to be
    // merely close -- do not read the tolerance as a claim about Courtois's
    // accuracy in general, which is governed by eps.
    const std::vector<double> exact = line::mc::ctmc_solve(m.mat);
    for (std::size_t i = 0; i < 4; ++i) {
        CAPTURE(i);
        CHECK(std::fabs(c.p[i] - exact[i]) <= 1e-9);
        CHECK(std::fabs(k.p[i] - exact[i]) <= 1e-9);
        CHECK(std::fabs(t.p[i] - exact[i]) <= 1e-9);
    }
    CHECK(c.eps < c.epsMAX);

    // 'multi' needs the second-level partition, so the scalar overload refuses
    // it by name rather than silently defaulting.
    CHECK_THROWS_AS(proc::aggregate(m, MS, "multi"), line::InputError);
    CHECK_THROWS_AS(proc::aggregate(m, MS, "nosuch"), line::InputError);
    const std::vector<std::vector<std::size_t>> MSS = {{0, 1}};
    CHECK(proc::aggregate_multi(m, MS, MSS).p.size() == 4u);
}

TEST_CASE("the kind guards refuse a method of the other class") {
    const MarkovChainModel<double> ctmc = fixture_process();
    const MarkovChainModel<double> dtmc = proc::to_embedded(ctmc);
    CHECK_THROWS_AS(proc::to_embedded(dtmc), line::InputError);
    CHECK_THROWS_AS(proc::to_markov_chain(dtmc), line::InputError);
    CHECK_THROWS_AS(proc::to_markov_process(ctmc), line::InputError);
    CHECK_THROWS_AS(proc::chain_transient_unif(ctmc, std::vector<double>(), 1.0),
                    line::InputError);
    CHECK_THROWS_AS(proc::chain_transient_steps(ctmc, std::vector<double>(), 1), line::InputError);
    CHECK_THROWS_AS(proc::time_average(dtmc, std::vector<double>(), 1.0), line::InputError);
    CHECK_THROWS_AS(proc::solve_relative(dtmc), line::InputError);
    CHECK_THROWS_AS(proc::aggregate(dtmc, {{0, 1}}), line::InputError);
}

TEST_CASE("isFeasible accepts both kinds and rejects a broken matrix") {
    CHECK(proc::is_feasible(fixture_process()));
    CHECK(proc::is_feasible(proc::to_embedded(fixture_process())));
    Matrix<double> bad(2, 2, 0.0);
    bad(0, 0) = 1.0; bad(0, 1) = -3.0;  // a negative off-diagonal rate
    MarkovChainModel<double> m;
    m.mat = bad;
    m.discrete = false;
    CHECK_FALSE(proc::is_feasible(m));
}

TEST_CASE("fromSampleSysAggr estimates a DTMC from an observed trajectory") {
    // Trajectory a,b,a,b,c over one state variable: a->b twice, b->a once,
    // b->c once, so row a is (0,1,0) and row b is (1/2,0,1/2).
    Matrix<double> obs(5, 1, 0.0);
    obs(0, 0) = 7.0; obs(1, 0) = 8.0; obs(2, 0) = 7.0; obs(3, 0) = 8.0; obs(4, 0) = 9.0;
    const MarkovChainModel<double> est = proc::from_sample_sys_aggr(obs);
    REQUIRE(est.discrete);
    REQUIRE(est.mat.rows() == 3);
    // First-appearance order: 7, 8, 9.
    CHECK(std::fabs(est.state_space(0, 0) - 7.0) <= 1e-15);
    CHECK(std::fabs(est.state_space(2, 0) - 9.0) <= 1e-15);
    CHECK(std::fabs(est.mat(0, 1) - 1.0) <= 1e-14);
    CHECK(std::fabs(est.mat(1, 0) - 0.5) <= 1e-14);
    CHECK(std::fabs(est.mat(1, 2) - 0.5) <= 1e-14);
    // Every row sums to one, since dtmc_makestochastic normalized the counts
    // and the unvisited-out state 9 was given a law rather than a zero row.
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) s += est.mat(i, j);
        CAPTURE(i);
        CHECK(std::fabs(s - 1.0) <= 1e-13);
    }
    CHECK_THROWS_AS(proc::from_sample_sys_aggr(Matrix<double>(1, 1, 0.0)), line::InputError);
}

TEST_CASE("sample draws a path of each kind") {
    std::mt19937_64 gen(20260815u);
    const MarkovChainModel<double> m = fixture_process();
    const proc::ChainPath<double> p = proc::chain_sample(m, std::vector<double>(), 50, gen);
    CHECK(p.states.size() == p.sojourn.size());
    REQUIRE(p.states.size() > 1);
    for (std::size_t i = 0; i < p.states.size(); ++i) {
        CAPTURE(i);
        CHECK(p.states[i] < 3u);
        CHECK(p.sojourn[i] > 0.0);
        // State 3 has no jump to state 2, so that pair must never occur.
        if (i + 1 < p.states.size()) {
            const bool forbidden_pair = (p.states[i] == 2u) && (p.states[i + 1] == 1u);
            CHECK_FALSE(forbidden_pair);
        }
    }

    const proc::ChainPath<double> d =
        proc::chain_sample(proc::to_embedded(m), std::vector<double>(), 20, gen);
    CHECK(d.sojourn.empty());
    REQUIRE(d.states.size() > 1);
    for (std::size_t i = 0; i < d.states.size(); ++i) CHECK(d.states[i] < 3u);
}

TEST_CASE("rand builds a valid chain of each kind") {
    std::mt19937_64 gen(11u);
    const MarkovChainModel<double> p = proc::rand_process<double>(5, gen);
    CHECK_FALSE(p.discrete);
    CHECK(proc::is_feasible(p));
    const MarkovChainModel<double> c = proc::rand_chain<double>(5, gen);
    CHECK(c.discrete);
    CHECK(proc::is_feasible(c));
}
