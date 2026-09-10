/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Aggregation-disaggregation and reducibility in the mc layer:
 * ctmc_courtois, ctmc_kms, ctmc_takahashi, ctmc_multi, ctmc_gmres,
 * ctmc_gmres_multi, ctmc_bicgstab, ctmc_bicgstab_multi, ctmc_solve_reducible, dtmc_solve_reducible,
 * ctmc_solve_reducible_blkdecomp, stronglyconncomp, ctmc_randomization,
 * dtmc_makestochastic and ctmc_rand.
 *
 * Oracles:
 *  - the exact stationary vector from ctmc_solve, which every method must
 *    reproduce to the accuracy it claims (Courtois to the degree of coupling,
 *    the iterative methods to convergence, the reducible solvers exactly);
 *  - closed-form absorption probabilities on a reducible chain, computed by
 *    hand from the embedded jump chain;
 *  - MATLAB reference values obtained by running LINE 3.0.6 on the same
 *    fixtures (matlab/src/api/mc, generated 2026-07-21).
 *
 * The GMRES fixture has a ZERO leading pivot, so its elimination cannot avoid
 * a row swap; a diagonally dominant fixture would exercise no pivoting at all.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_bicgstab.h"
#include "line/api/mc/ctmc_gmres.h"
#include "line/api/mc/ctmc_gmres_multi.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_multi.h"
#include "line/api/mc/ctmc_rand.h"
#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/ctmc_solve_reducible_blkdecomp.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/util/lu.h"

using line::Matrix;
using line::Rational;
using line::mc::ctmc_courtois;
using line::mc::ctmc_gmres;
using line::mc::ctmc_gmres_multi;
using line::mc::ctmc_kms;
using line::mc::ctmc_makeinfgen;
using line::mc::ctmc_multi;
using line::mc::ctmc_rand;
using line::mc::ctmc_randomization;
using line::mc::ctmc_solve;
using line::mc::ctmc_solve_reducible;
using line::mc::ctmc_solve_reducible_blkdecomp;
using line::mc::ctmc_takahashi;
using line::mc::dtmc_makestochastic;
using line::mc::dtmc_solve_reducible;
using line::mc::LcgUniform;
using line::mc::stronglyconncomp;

namespace {

using Blocks = std::vector<std::vector<std::size_t>>;

/**
 * Nearly completely decomposable chain A: three macro-states of two states,
 * intra-block rates of order one and coupling of order 0.01.
 */
template <class T>
Matrix<T> ncd_chain_A() {
    Matrix<T> Q(6, 6, line::num_traits<T>::from_int(0));
    Q(0, 1) = line::num_traits<T>::from_int(3);
    Q(1, 0) = line::num_traits<T>::from_int(2);
    Q(2, 3) = line::num_traits<T>::from_int(5);
    Q(3, 2) = line::num_traits<T>::from_int(1);
    Q(4, 5) = line::num_traits<T>::from_int(2);
    Q(5, 4) = line::num_traits<T>::from_int(4);
    Q(1, 2) = line::num_traits<T>::from_rational(1, 100);
    Q(2, 1) = line::num_traits<T>::from_rational(5, 1000);
    Q(3, 4) = line::num_traits<T>::from_rational(2, 100);
    Q(5, 0) = line::num_traits<T>::from_rational(3, 100);
    return ctmc_makeinfgen(Q);
}

/**
 * Chain C: the same skeleton, but with two couplings landing on the SAME state,
 * so the largest column sum of the coupling matrix (MATLAB's eps) and its
 * largest row sum (the JAR's) differ and can be told apart.
 */
template <class T>
Matrix<T> ncd_chain_C() {
    Matrix<T> Q(6, 6, line::num_traits<T>::from_int(0));
    Q(0, 1) = line::num_traits<T>::from_int(3);
    Q(1, 0) = line::num_traits<T>::from_int(2);
    Q(2, 3) = line::num_traits<T>::from_int(5);
    Q(3, 2) = line::num_traits<T>::from_int(1);
    Q(4, 5) = line::num_traits<T>::from_int(2);
    Q(5, 4) = line::num_traits<T>::from_int(4);
    Q(1, 2) = line::num_traits<T>::from_rational(1, 100);
    Q(3, 4) = line::num_traits<T>::from_rational(2, 100);
    Q(4, 0) = line::num_traits<T>::from_rational(4, 100);
    Q(5, 0) = line::num_traits<T>::from_rational(3, 100);
    return ctmc_makeinfgen(Q);
}

/**
 * Reducible chain B: transient class {0,1}, recurrent class {2,3} and the
 * absorbing state 4.
 */
template <class T>
Matrix<T> reducible_chain_B() {
    Matrix<T> Q(5, 5, line::num_traits<T>::from_int(0));
    Q(0, 1) = line::num_traits<T>::from_int(1);
    Q(0, 4) = line::num_traits<T>::from_rational(1, 4);
    Q(1, 0) = line::num_traits<T>::from_int(1);
    Q(1, 2) = line::num_traits<T>::from_rational(1, 2);
    Q(2, 3) = line::num_traits<T>::from_int(2);
    Q(3, 2) = line::num_traits<T>::from_int(1);
    return ctmc_makeinfgen(Q);
}

Blocks contiguous_pairs() { return Blocks{{0, 1}, {2, 3}, {4, 5}}; }

double linf(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

double sum_of(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Building blocks
// ---------------------------------------------------------------------------

TEST_CASE("dtmc_makestochastic normalizes rows and closes dead ones") {
    Matrix<double> P{{1.0, 3.0}, {0.0, 0.0}};
    Matrix<double> S = dtmc_makestochastic(P);
    CHECK(S(0, 0) == doctest::Approx(0.25));
    CHECK(S(0, 1) == doctest::Approx(0.75));
    CHECK(S(1, 0) == doctest::Approx(0.0));
    CHECK(S(1, 1) == doctest::Approx(1.0));

    // Exact: a rational row normalizes to rationals with no rounding at all.
    Matrix<Rational> Pq{{Rational(1), Rational(3)}, {Rational(0), Rational(0)}};
    Matrix<Rational> Sq = dtmc_makestochastic(Pq);
    CHECK(Sq(0, 0) == Rational(1, 4));
    CHECK(Sq(0, 1) == Rational(3, 4));
    CHECK(Sq(1, 1) == Rational(1));
}

TEST_CASE("ctmc_randomization produces a stochastic matrix with the same stationary vector") {
    const Matrix<double> Q = ncd_chain_A<double>();
    const line::mc::RandomizationResult<double> r = ctmc_randomization(Q);
    for (std::size_t i = 0; i < 6; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 6; ++j) {
            CHECK(r.P(i, j) >= 0.0);
            s += r.P(i, j);
        }
        CHECK(s == doctest::Approx(1.0).epsilon(1e-14));
    }
    // The uniformized chain has the same stationary distribution as Q.
    const std::vector<double> piQ = ctmc_solve(Q);
    const std::vector<double> piP = line::mc::dtmc_solve(r.P);
    CHECK(linf(piQ, piP) < 1e-12);

    // Exact: the rate is (21/20) max|Q|, so the rational path stays rational.
    const Matrix<Rational> Qq = ncd_chain_A<Rational>();
    const line::mc::RandomizationResult<Rational> rq = ctmc_randomization(Qq);
    CHECK(rq.q == Rational(21, 20) * Rational(5005, 1000));
    CHECK(line::mc::dtmc_solve(rq.P) == ctmc_solve(Qq));
}

TEST_CASE("ctmc_rand is reproducible for a given generator and is a valid generator") {
    LcgUniform g1(12345), g2(12345), g3(999);
    const Matrix<double> A = ctmc_rand<double>(5, g1);
    const Matrix<double> B = ctmc_rand<double>(5, g2);
    const Matrix<double> C = ctmc_rand<double>(5, g3);
    bool identical = true, differs = false;
    for (std::size_t i = 0; i < 5; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 5; ++j) {
            if (A(i, j) != B(i, j)) identical = false;
            if (A(i, j) != C(i, j)) differs = true;
            s += A(i, j);
            if (i != j) CHECK(A(i, j) >= 0.0);
        }
        CHECK(s == doctest::Approx(0.0).epsilon(1e-14));
        CHECK(A(i, i) <= 0.0);
    }
    CHECK(identical);
    CHECK(differs);

    // The variates are dyadic, so the exact instantiation is exact.
    LcgUniform g4(12345);
    const Matrix<Rational> R = ctmc_rand<Rational>(5, g4);
    for (std::size_t i = 0; i < 5; ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < 5; ++j) s += R(i, j);
        CHECK(s == Rational(0));
    }
}

TEST_CASE("stronglyconncomp finds the components and their recurrence") {
    const Matrix<double> Q = reducible_chain_B<double>();
    Matrix<double> Adj = Q;
    for (std::size_t i = 0; i < 5; ++i) Adj(i, i) = 0.0;
    const line::mc::SccResult s = stronglyconncomp(Adj);
    REQUIRE(s.numSCC() == 3);
    // MATLAB stronglyconncomp on this chain returns scc = [1 1 2 2 3] with
    // recurrent = [false true true]; components are numbered by decreasing size.
    CHECK(s.scc[0] == 1);
    CHECK(s.scc[1] == 1);
    CHECK(s.scc[2] == 2);
    CHECK(s.scc[3] == 2);
    CHECK(s.scc[4] == 3);
    CHECK(s.recurrent[0] == false);
    CHECK(s.recurrent[1] == true);
    CHECK(s.recurrent[2] == true);
}

// ---------------------------------------------------------------------------
// Reducible solvers: exact, and correct on a genuinely reducible chain
// ---------------------------------------------------------------------------

TEST_CASE("reducible solvers reproduce ctmc_solve EXACTLY on an irreducible chain") {
    const Matrix<Rational> Q = ncd_chain_A<Rational>();
    const std::vector<Rational> pi = ctmc_solve(Q);

    // A single strongly connected component: both routines must return the very
    // same rationals, not merely agree to a tolerance.
    const std::vector<Rational> pr = ctmc_solve_reducible(Q).pi;
    const std::vector<Rational> pb = ctmc_solve_reducible_blkdecomp(Q).pi;
    REQUIRE(pr.size() == pi.size());
    REQUIRE(pb.size() == pi.size());
    for (std::size_t i = 0; i < pi.size(); ++i) {
        CHECK(pr[i] == pi[i]);
        CHECK(pb[i] == pi[i]);
    }
    Rational s(0);
    for (const Rational& v : pi) s += v;
    CHECK(s == Rational(1));
}

TEST_CASE("ctmc_solve_reducible_blkdecomp matches the hand-computed absorption") {
    // Jump chain: from 0 absorb into state 4 with probability 1/5, else go to
    // 1; from 1 enter the recurrent class {2,3} with probability 1/3, else go
    // to 0. Hence h0 = 4/7 and h1 = 5/7 for reaching {2,3}; starting uniformly
    // in the transient class gives 9/14 there and 5/14 at state 4. Within
    // {2,3} the stationary vector is (1/3, 2/3).
    const Matrix<Rational> Q = reducible_chain_B<Rational>();
    const line::mc::BlkDecompResult<Rational> r = ctmc_solve_reducible_blkdecomp(Q);
    CHECK(r.pi[0] == Rational(0));
    CHECK(r.pi[1] == Rational(0));
    CHECK(r.pi[2] == Rational(9, 14) * Rational(1, 3));
    CHECK(r.pi[3] == Rational(9, 14) * Rational(2, 3));
    CHECK(r.pi[4] == Rational(5, 14));
    Rational s(0);
    for (const Rational& v : r.pi) s += v;
    CHECK(s == Rational(1));

    // MATLAB ctmc_solve_reducible_blkdecomp on the same generator, to the
    // double round-off of the reference (it solves the same finite system).
    const line::mc::BlkDecompResult<double> rd = ctmc_solve_reducible_blkdecomp(reducible_chain_B<double>());
    const std::vector<double> matlab{0.0, 0.0, 0.2142857142857143, 0.4285714285714286, 0.35714285714285721};
    CHECK(linf(rd.pi, matlab) < 1e-14);
}

TEST_CASE("dtmc_solve_reducible lumps the reducible chain as MATLAB does") {
    // The uniformized chain at rate 2.1, the matrix MATLAB's ctmc_solve_reducible
    // hands to dtmc_solve_reducible.
    const Matrix<Rational> Q = reducible_chain_B<Rational>();
    const Matrix<Rational> P = ctmc_randomization(Q, Rational(21, 10)).P;
    const line::mc::ReducibleResult<Rational> r = dtmc_solve_reducible(P);

    // The lumped chain absorbs into {2,3} with probability 2/3 and into {4}
    // with probability 1/3, so pi = (0, 0, 2/9, 4/9, 1/3). This is EXACT here,
    // where MATLAB reaches it by a spectral decomposition with a power-method
    // fallback capped at 1000 iterations.
    CHECK(r.pi[0] == Rational(0));
    CHECK(r.pi[1] == Rational(0));
    CHECK(r.pi[2] == Rational(2, 9));
    CHECK(r.pi[3] == Rational(4, 9));
    CHECK(r.pi[4] == Rational(1, 3));

    // Rows of pis: starting in the transient component, in {2,3}, in {4}.
    CHECK(r.pis(1, 2) == Rational(1, 3));
    CHECK(r.pis(1, 3) == Rational(2, 3));
    CHECK(r.pis(2, 4) == Rational(1));

    // MATLAB dtmc_solve_reducible on the same stochastic matrix.
    const line::mc::ReducibleResult<double> rd =
        dtmc_solve_reducible(ctmc_randomization(reducible_chain_B<double>(), 2.1).P);
    const std::vector<double> matlab{0.0, 0.0, 0.22222222222222227, 0.44444444444444453,
                                     0.33333333333333337};
    CHECK(linf(rd.pi, matlab) < 1e-14);
}

/**
 * The single-transient-component branch reads a row of `pis` BY COMPONENT INDEX,
 * so every row has to be filled -- including the rows of components that carry
 * no starting mass. Filling only the rows with positive weight left that branch
 * reading a row of zeros and returning an all-zero "distribution".
 *
 * The fixture is the routing chain the H-T fork-join transform produces for the
 * second fork of a two-fork model (fj_serialfjs_closed): three isolated states
 * with no incoming mass, one transient state and one three-cycle. The reference
 * value is `dtmc_solve_reducible(P, [], struct('tol', GlobalConstants.FineTol))`
 * in MATLAB on the same matrix.
 */
TEST_CASE("dtmc_solve_reducible fills every component row, not only the weighted ones") {
    Matrix<Rational> P(7, 7, Rational(0));
    P(3, 5) = Rational(1);  // the three-cycle 3 -> 5 -> 6 -> 3
    P(5, 6) = Rational(1);
    P(6, 3) = Rational(1);
    P(4, 5) = Rational(1);  // the one transient state, absorbed into the cycle
    // states 0, 1 and 2 are isolated: no incoming mass, so weight zero
    const line::mc::ReducibleResult<Rational> r = dtmc_solve_reducible(P);
    CHECK(r.pi[0] == Rational(0));
    CHECK(r.pi[1] == Rational(0));
    CHECK(r.pi[2] == Rational(0));
    CHECK(r.pi[3] == Rational(1, 3));
    CHECK(r.pi[4] == Rational(0));
    CHECK(r.pi[5] == Rational(1, 3));
    CHECK(r.pi[6] == Rational(1, 3));
}

TEST_CASE("ctmc_solve_reducible and the block decomposition disagree, as the references do") {
    // MATLAB's two reducible solvers give different answers on this chain:
    // ctmc_solve_reducible lumps transition MASS between components and loses
    // the sojourn structure inside the transient class, while the block
    // decomposition solves for the expected sojourn and is exact. The port
    // reproduces both, so the divergence is visible rather than papered over.
    const Matrix<Rational> Q = reducible_chain_B<Rational>();
    const std::vector<Rational> viaUniformization = ctmc_solve_reducible(Q).pi;
    const std::vector<Rational> viaBlocks = ctmc_solve_reducible_blkdecomp(Q).pi;
    CHECK(viaUniformization[4] == Rational(1, 3));   // MATLAB 0.33333333333333337
    CHECK(viaBlocks[4] == Rational(5, 14));          // MATLAB 0.35714285714285721
    CHECK(viaUniformization[4] != viaBlocks[4]);
}

// ---------------------------------------------------------------------------
// Aggregation-disaggregation
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_courtois approximates the stationary vector to the degree of coupling") {
    const Matrix<double> Q = ncd_chain_A<double>();
    const std::vector<double> exact = ctmc_solve(Q);
    const line::mc::CourtoisResult<double> r = ctmc_courtois(Q, contiguous_pairs());

    CHECK(sum_of(r.p) == doctest::Approx(1.0).epsilon(1e-12));
    // Courtois's bound: the error is O(eps), the degree of coupling, and eps
    // must be well below epsMAX for the decomposition to be meaningful.
    CHECK(r.eps < r.epsMAX);
    CHECK(linf(r.p, exact) < 10.0 * r.eps);

    // MATLAB ctmc_courtois(QA, MS): the construction is a finite sequence of
    // linear solves, so agreement is at double round-off, not at eps.
    const std::vector<double> matlab{0.20895522388059701, 0.31343283582089548, 0.029850746268656719,
                                     0.14925373134328357, 0.19900497512437809, 0.099502487562189046};
    CHECK(linf(r.p, matlab) < 1e-12);
    CHECK(r.eps == doctest::Approx(0.0057085771371485662).epsilon(1e-12));
    CHECK(r.epsMAX == doctest::Approx(0.42914228628514339).epsilon(1e-10));
    CHECK(r.q == doctest::Approx(5.2552500000000002).epsilon(1e-14));
}

TEST_CASE("the NCD index is the ROW sum, and the old column sum is still reported") {
    // The degree of coupling ||B||_inf is the largest ROW sum, which is what
    // Ctmc_courtois.java always computed and what MATLAB computes since
    // 2026-08-15 (`max(sum(B,2))`; it was `max(sum(B))`, the largest COLUMN
    // sum, before that). On chain C two couplings land on state 0, which
    // separates the two definitions -- and it is what makes this fixture worth
    // keeping. Live MATLAB on this generator, 2026-08-15: eps 0.007619047619047619,
    // column max 0.013333333333333334.
    const line::mc::CourtoisResult<double> r = ctmc_courtois(ncd_chain_C<double>(), contiguous_pairs());
    CHECK(r.eps == doctest::Approx(0.007619047619047619).epsilon(1e-12));        // MATLAB and JAR
    CHECK(r.epsRowMax == doctest::Approx(0.007619047619047619).epsilon(1e-12));
    CHECK(r.epsColMax == doctest::Approx(0.013333333333333334).epsilon(1e-12));  // pre-2026-08-15
    // The column sum UNDERSTATES the coupling here, which is why it was the
    // wrong one to report: a partition looks 1.75x more decomposable than it is.
    CHECK(r.epsColMax > r.eps);
    CHECK(r.epsMAX == doctest::Approx(0.4285714285714286).epsilon(1e-10));

    const std::vector<double> matlab{0.26252983293556093, 0.3937947494033413, 0.039379474940334135,
                                     0.19689737470167068, 0.07159904534606204, 0.03579952267303102};
    CHECK(linf(r.p, matlab) < 1e-12);
}

TEST_CASE("ctmc_kms converges to the exact stationary vector") {
    const Matrix<double> Q = ncd_chain_A<double>();
    const std::vector<double> exact = ctmc_solve(Q);
    const line::mc::KmsResult<double> r0 = ctmc_kms(Q, contiguous_pairs(), 0);
    const line::mc::KmsResult<double> r3 = ctmc_kms(Q, contiguous_pairs(), 3);

    // Zero sweeps is exactly the Courtois starting point.
    CHECK(linf(r0.p, ctmc_courtois(Q, contiguous_pairs()).p) < 1e-14);
    // Three sweeps drive the iterate to the exact vector: the fixed point of
    // the aggregation-disaggregation map is the stationary distribution.
    CHECK(sum_of(r3.p) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(linf(r3.p, exact) < 1e-12);

    // MATLAB ctmc_kms(QA, MS, 3).
    const std::vector<double> matlab{0.2093991671624047, 0.31261154074955577, 0.030339083878643813,
                                     0.14872099940511671, 0.19978187586753635, 0.099147332936742591};
    CHECK(linf(r3.p, matlab) < 1e-12);
}

TEST_CASE("ctmc_takahashi converges to the exact stationary vector") {
    const Matrix<double> Q = ncd_chain_A<double>();
    const std::vector<double> exact = ctmc_solve(Q);
    const line::mc::TakahashiResult<double> r = ctmc_takahashi(Q, contiguous_pairs(), 3);
    CHECK(sum_of(r.p) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(linf(r.p, exact) < 1e-11);

    // MATLAB ctmc_takahashi(QA, MS, 3); the reference's own iterate carries
    // about 1e-12 of accumulated round-off against ctmc_solve.
    const std::vector<double> matlab{0.20939916716240181, 0.3126115407495515, 0.030339083878644257,
                                     0.1487209994051189, 0.19978187586753954, 0.099147332936744187};
    CHECK(linf(r.p, matlab) < 1e-11);
}

TEST_CASE("aggregation methods address macro-states by index, not by position") {
    // Chain A with the states interleaved so that no macro-state is a
    // contiguous range. The answer must be the answer for the natural ordering,
    // permuted; a routine that walks contiguous block offsets instead of the
    // index sets silently solves a different partition of the same sizes.
    const Matrix<double> Q = ncd_chain_A<double>();
    const std::size_t sig[6] = {0, 2, 4, 1, 3, 5};  // new position -> old state
    Matrix<double> QP(6, 6);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) QP(i, j) = Q(sig[i], sig[j]);
    const Blocks MSP{{0, 3}, {1, 4}, {2, 5}};

    const std::vector<double> pNat = ctmc_courtois(Q, contiguous_pairs()).p;
    const std::vector<double> pPerm = ctmc_courtois(QP, MSP).p;
    for (std::size_t i = 0; i < 6; ++i) CHECK(pPerm[i] == doctest::Approx(pNat[sig[i]]).epsilon(1e-12));

    // MATLAB ctmc_courtois on the permuted chain with the interleaved blocks.
    const std::vector<double> matlabCourtois{0.20895522388059701, 0.029850746268656719,
                                             0.19900497512437809, 0.31343283582089548,
                                             0.14925373134328357, 0.099502487562189046};
    CHECK(linf(pPerm, matlabCourtois) < 1e-12);

    // MATLAB ctmc_takahashi on the same permuted chain, two sweeps.
    const std::vector<double> matlabTaka{0.20939916716240176, 0.03033908387864425,
                                         0.19978187586753951, 0.31261154074955139,
                                         0.14872099940511888, 0.099147332936744159};
    CHECK(linf(ctmc_takahashi(QP, MSP, 2).p, matlabTaka) < 1e-11);
}

TEST_CASE("ctmc_multi solves the coarse chain by a second Courtois step") {
    const Matrix<double> Q = ncd_chain_A<double>();
    const Blocks MSS{{0, 1}, {2}};
    const line::mc::MultiResult<double> r = ctmc_multi(Q, contiguous_pairs(), MSS);
    CHECK(sum_of(r.p) == doctest::Approx(1.0).epsilon(1e-12));

    // MATLAB ctmc_multi(QA, MS, MSS). The two-level estimate is far from the
    // stationary vector here because the COARSE chain is not itself nearly
    // completely decomposable; that is a property of the method on this
    // partition, and the port reproduces the reference value.
    const std::vector<double> matlab{0.019801980198019806, 0.029702970297029702,
                                     0.05940594059405941,  0.29702970297029702,
                                     0.39603960396039595,  0.19801980198019797};
    CHECK(linf(r.p, matlab) < 1e-12);
    // The plain Courtois estimate it also returns is the accurate one.
    CHECK(linf(r.pcourt, ctmc_solve(Q)) < 0.01);
}

// ---------------------------------------------------------------------------
// GMRES
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_gmres solves a system whose pivot sequence is non-trivial") {
    // A(0,0) = 0, so the elimination must swap rows; the answer is checked
    // against the direct solve and against MATLAB's ctmc_gmres.
    Matrix<double> A{{0.0, 2.0, 1.0}, {1.0, 0.0, 3.0}, {4.0, 5.0, 0.0}};
    const std::vector<double> b{1.0, 2.0, 3.0};
    const line::mc::GmresResult<double> g = ctmc_gmres(A, b);
    CHECK(g.flag == 0);
    const std::vector<double> direct = line::solve(A, b);
    CHECK(linf(g.x, direct) < 1e-10);

    const std::vector<double> matlab{0.44827586206896541, 0.24137931034482762, 0.51724137931034475};
    CHECK(linf(g.x, matlab) < 1e-10);

    // The residual it reports is the true one, so it can be checked directly.
    double res = 0.0, bn = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) s += A(i, j) * g.x[j];
        res += (b[i] - s) * (b[i] - s);
        bn += b[i] * b[i];
    }
    CHECK(std::sqrt(res / bn) < 1e-10);
}

TEST_CASE("ctmc_gmres solves a stiff birth-death balance system") {
    // The normalization row is O(1) while the rates span three decades, which
    // is the case row equilibration exists for.
    const std::size_t n = 40;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = 2.0;
        Q(i + 1, i) = 3.0 + 0.05 * static_cast<double>(i);
    }
    const Matrix<double> G = ctmc_makeinfgen(Q);
    Matrix<double> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = (i == n - 1) ? 1.0 : G(j, i);
    std::vector<double> b(n, 0.0);
    b[n - 1] = 1.0;

    const line::mc::GmresResult<double> g = ctmc_gmres(A, b);
    CHECK(g.flag == 0);
    CHECK(linf(g.x, ctmc_solve(G)) < 1e-9);
}

TEST_CASE("ctmc_gmres_multi solves every column against one factorization") {
    Matrix<double> A{{0.0, 2.0, 1.0}, {1.0, 0.0, 3.0}, {4.0, 5.0, 0.0}};
    Matrix<double> B{{1.0, 0.0}, {2.0, 1.0}, {3.0, 1.0}};
    // The tolerance here is 1e-10 rather than the default 1e-12 because this
    // port measures the TRUE residual (right preconditioning) where MATLAB
    // measures the preconditioned one (left preconditioning): the second column
    // of this system settles at a true relative residual of about 2e-12, which
    // is the round-off floor of a matrix with this conditioning, so 1e-12 is
    // simply unattainable in double precision and MATLAB's flag = 0 on the same
    // call reflects a different, preconditioner-dependent quantity.
    const line::mc::GmresMultiResult<double> r = ctmc_gmres_multi(A, B, 1e-10);
    REQUIRE(r.flag == 0);
    REQUIRE(r.X.rows() == 3);
    REQUIRE(r.X.cols() == 2);

    // MATLAB ctmc_gmres_multi(AG, BG), column-major.
    const double matlab[3][2] = {{0.44827586206896541, 0.37931034482758619},
                                 {0.24137931034482762, -0.10344827586206884},
                                 {0.51724137931034475, 0.20689655172413801}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(r.X(i, j) == doctest::Approx(matlab[i][j]).epsilon(1e-9));

    // Every column must also satisfy its own equation.
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 3; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < 3; ++j) s += A(i, j) * r.X(j, c);
            CHECK(s == doctest::Approx(B(i, c)).epsilon(1e-9));
        }
}

TEST_CASE("ctmc_bicgstab reproduces the direct solve on the non-trivial pivot fixture") {
    // The same system as the GMRES fixture: A(0,0) = 0, so the elimination must
    // swap rows. Both methods must return the same solution, which is the point
    // of sharing the equilibration, the reordering and the factorization.
    Matrix<double> A{{0.0, 2.0, 1.0}, {1.0, 0.0, 3.0}, {4.0, 5.0, 0.0}};
    const std::vector<double> b{1.0, 2.0, 3.0};
    const line::mc::BicgstabResult<double> s = line::mc::ctmc_bicgstab(A, b);
    CHECK(s.flag == 0);
    CHECK(linf(s.x, line::solve(A, b)) < 1e-10);

    const std::vector<double> matlab{0.44827586206896541, 0.24137931034482762, 0.51724137931034475};
    CHECK(linf(s.x, matlab) < 1e-10);

    // The residual it reports is the true one, so it can be checked directly.
    double res = 0.0, bn = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        double t = 0.0;
        for (std::size_t j = 0; j < 3; ++j) t += A(i, j) * s.x[j];
        res += (b[i] - t) * (b[i] - t);
        bn += b[i] * b[i];
    }
    CHECK(std::sqrt(res / bn) < 1e-10);

    // Two matrix-vector products per complete iteration, so the reported count
    // is even and comparable with the one GMRES reports.
    CHECK(s.iter % 2 == 0);
}

TEST_CASE("ctmc_bicgstab solves a stiff birth-death balance system") {
    // The normalization row is O(1) while the rates span three decades, which
    // is the case row equilibration exists for.
    const std::size_t n = 40;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = 2.0;
        Q(i + 1, i) = 3.0 + 0.05 * static_cast<double>(i);
    }
    const Matrix<double> G = ctmc_makeinfgen(Q);
    Matrix<double> A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = (i == n - 1) ? 1.0 : G(j, i);
    std::vector<double> b(n, 0.0);
    b[n - 1] = 1.0;

    const line::mc::BicgstabResult<double> s = line::mc::ctmc_bicgstab(A, b);
    CHECK(s.flag == 0);
    CHECK(linf(s.x, ctmc_solve(G)) < 1e-9);
}

TEST_CASE("ctmc_bicgstab matches the closed-form M/M/1/K stationary vector") {
    // An ANALYTICAL oracle rather than a recorded baseline: the stationary
    // vector of M/M/1/K is the truncated geometric, so a consistent error in
    // both Krylov kernels cannot pass this test the way it would pass a
    // cross-check against one another.
    const std::size_t K = 400;
    const double lambda = 0.7, mu = 1.0;
    const std::size_t n = K + 1;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = lambda;
        Q(i + 1, i) = mu;
    }
    const Matrix<double> G = ctmc_makeinfgen(Q);
    // Qnnz(:,end) = 1 then solve Qnnz' x = e_n, which is what ctmc_solve does.
    Matrix<double> A(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = (i == n - 1) ? 1.0 : G(j, i);
    std::vector<double> b(n, 0.0);
    b[n - 1] = 1.0;

    const double rho = lambda / mu;
    std::vector<double> exact(n);
    double geom = 0.0, power = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        exact[i] = power;
        geom += power;
        power *= rho;
    }
    for (std::size_t i = 0; i < n; ++i) exact[i] /= geom;

    const line::mc::BicgstabResult<double> s = line::mc::ctmc_bicgstab(A, b);
    CHECK(s.flag == 0);
    CHECK(linf(s.x, exact) < 1e-12);

    const line::mc::GmresResult<double> g = ctmc_gmres(A, b);
    CHECK(g.flag == 0);
    CHECK(linf(g.x, exact) < 1e-12);
}

TEST_CASE("ctmc_bicgstab_multi solves every column against one factorization") {
    Matrix<double> A{{0.0, 2.0, 1.0}, {1.0, 0.0, 3.0}, {4.0, 5.0, 0.0}};
    Matrix<double> B{{1.0, 0.0}, {2.0, 1.0}, {3.0, 1.0}};
    // Tolerance 1e-10 for the reason given on the GMRES multi-column test: this
    // port measures the TRUE residual, whose round-off floor on this system is
    // about 2e-12, so the default 1e-12 is not attainable in double precision.
    const line::mc::BicgstabMultiResult<double> r = line::mc::ctmc_bicgstab_multi(A, B, 1e-10);
    REQUIRE(r.flag == 0);
    REQUIRE(r.X.rows() == 3);
    REQUIRE(r.X.cols() == 2);

    // MATLAB ctmc_gmres_multi(AG, BG) on the same system, column-major: the
    // solution does not depend on which Krylov method reaches it.
    const double matlab[3][2] = {{0.44827586206896541, 0.37931034482758619},
                                 {0.24137931034482762, -0.10344827586206884},
                                 {0.51724137931034475, 0.20689655172413801}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(r.X(i, j) == doctest::Approx(matlab[i][j]).epsilon(1e-9));

    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t i = 0; i < 3; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < 3; ++j) s += A(i, j) * r.X(j, c);
            CHECK(s == doctest::Approx(B(i, c)).epsilon(1e-9));
        }
}

TEST_CASE("ctmc_solve takes the Krylov path above GMRES_MIN_STATES and stays exact") {
    // The gate is a SIZE, so the only way to exercise it is a generator past it.
    // M/M/1/K is the fixture because its stationary vector is the truncated
    // geometric in closed form: the answer is checked against arithmetic, not
    // against the LU this path replaces, so a fault shared by both would not pass.
    const std::size_t K = line::mc::GMRES_MIN_STATES + 2000;
    const double lambda = 0.7, mu = 1.0;
    const std::size_t n = K + 1;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = lambda;
        Q(i + 1, i) = mu;
    }
    const Matrix<double> G = ctmc_makeinfgen(Q);
    const std::vector<double> pi = ctmc_solve(G);
    REQUIRE(pi.size() == n);

    const double rho = lambda / mu;
    std::vector<double> exact(n);
    double total = 0.0, power = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        exact[i] = power;
        total += power;
        power *= rho;
    }
    for (std::size_t i = 0; i < n; ++i) exact[i] /= total;
    CHECK(linf(pi, exact) < 1e-11);
}
