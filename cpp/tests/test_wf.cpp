/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Workflow pattern detectors (api/wf).
 *
 * Oracles: the structural measures are read off graphs drawn by hand here, so
 * every count, chain and pattern is known independently of the code, and the
 * counting identities (total nodes = sum of chain lengths, mean = total /
 * count) are asserted exactly in the rational instantiation.
 *
 * Reference values were produced by running
 * jar/src/main/java/jline/api/wf/*.java against common/jline.jar on the same
 * hand-drawn graphs. The detectors are exact combinatorics, so the port must
 * agree with the JAR to the last bit; the entropy-based diversity metrics are
 * compared at 1e-15 relative, the rounding level of a sum of p log p, since
 * the method carries no tolerance of its own.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/wf/wf_branch_detector.h"
#include "line/api/wf/wf_loop_detector.h"
#include "line/api/wf/wf_parallel_detector.h"
#include "line/api/wf/wf_sequence_detector.h"

using line::Matrix;
using line::Rational;

namespace {

std::vector<int> ids(int a, int b = -1, int c = -1) {
    std::vector<int> v;
    v.push_back(a);
    if (b >= 0) v.push_back(b);
    if (c >= 0) v.push_back(c);
    return v;
}

/** 1 -> 2 -> 3 -> 4 -> 9, services {2,3,4}: one chain of three. */
template <class T>
Matrix<T> sequence_graph() {
    Matrix<T> m(4, 3, line::num_traits<T>::from_int(0));
    const int e[4][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 9}};
    for (std::size_t i = 0; i < 4; ++i) {
        m(i, 0) = line::num_traits<T>::from_int(e[i][0]);
        m(i, 1) = line::num_traits<T>::from_int(e[i][1]);
        m(i, 2) = line::num_traits<T>::from_int(1);
    }
    return m;
}

/** 1 -> 6(fork) -> {2,3,5} -> 7(join) -> 9: one parallel pattern of three. */
template <class T>
Matrix<T> parallel_graph() {
    const int e[8][2] = {{1, 6}, {6, 2}, {6, 3}, {6, 5}, {2, 7}, {3, 7}, {5, 7}, {7, 9}};
    Matrix<T> m(8, 3, line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < 8; ++i) {
        m(i, 0) = line::num_traits<T>::from_int(e[i][0]);
        m(i, 1) = line::num_traits<T>::from_int(e[i][1]);
        m(i, 2) = line::num_traits<T>::from_int(1);
    }
    return m;
}

/** 1 -> 2 -> 8(router), 8 -> 2 with p = 3/10, 8 -> 4 with 7/10, 4 -> 9. */
template <class T>
Matrix<T> loop_graph() {
    const int e[5][2] = {{1, 2}, {2, 8}, {8, 2}, {8, 4}, {4, 9}};
    Matrix<T> m(5, 3, line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < 5; ++i) {
        m(i, 0) = line::num_traits<T>::from_int(e[i][0]);
        m(i, 1) = line::num_traits<T>::from_int(e[i][1]);
        m(i, 2) = line::num_traits<T>::from_int(1);
    }
    m(2, 2) = line::num_traits<T>::from_rational(3, 10);
    m(3, 2) = line::num_traits<T>::from_rational(7, 10);
    return m;
}

/** 1 -> 8, 8 -> 3 with p = 1/4, 8 -> 4 with 3/4, both into 7(join) -> 9. */
template <class T>
Matrix<T> branch_graph() {
    const int e[6][2] = {{1, 8}, {8, 3}, {8, 4}, {3, 7}, {4, 7}, {7, 9}};
    Matrix<T> m(6, 3, line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < 6; ++i) {
        m(i, 0) = line::num_traits<T>::from_int(e[i][0]);
        m(i, 1) = line::num_traits<T>::from_int(e[i][1]);
        m(i, 2) = line::num_traits<T>::from_int(1);
    }
    m(1, 2) = line::num_traits<T>::from_rational(1, 4);
    m(2, 2) = line::num_traits<T>::from_rational(3, 4);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// sequences
// ---------------------------------------------------------------------------

TEST_CASE("detect_sequences finds the hand-drawn chain") {
    const Matrix<Rational> g = sequence_graph<Rational>();
    const std::vector<std::vector<int>> chains =
        line::wf::detect_sequences(g, ids(2, 3, 4));
    // JAR reference: sequences = [[2, 3, 4]]
    REQUIRE(chains.size() == 1u);
    REQUIRE(chains[0].size() == 3u);
    CHECK(chains[0][0] == 2);
    CHECK(chains[0][1] == 3);
    CHECK(chains[0][2] == 4);
}

TEST_CASE("sequence statistics obey the counting identities exactly") {
    const Matrix<Rational> g = sequence_graph<Rational>();
    const std::vector<std::vector<int>> chains = line::wf::detect_sequences(g, ids(2, 3, 4));
    const line::wf::SequenceStats<Rational> s = line::wf::get_sequence_stats<Rational>(chains);
    // JAR reference: {avgLength=3.0, maxLength=3, minLength=3,
    //                 numSequences=1, totalNodes=3}
    CHECK(s.numSequences == 1u);
    CHECK(s.totalNodes == 3u);
    CHECK(s.avgLength == Rational(3));
    CHECK(s.maxLength == 3u);
    CHECK(s.minLength == 3u);
    // STRUCTURAL LAW, exact: total = sum of the lengths, mean = total / count.
    std::size_t total = 0;
    for (std::size_t i = 0; i < chains.size(); ++i) total += chains[i].size();
    CHECK(s.totalNodes == total);
    CHECK(s.avgLength * Rational(static_cast<long>(s.numSequences)) ==
          Rational(static_cast<long>(s.totalNodes)));
    CHECK(s.minLength <= s.maxLength);
}

TEST_CASE("detect_sequences separates two disjoint chains") {
    // 2 -> 3 and 5 -> 4: two endpoints each, so two chains of two.
    Matrix<Rational> m(2, 3, Rational(0));
    m(0, 0) = Rational(2);
    m(0, 1) = Rational(3);
    m(0, 2) = Rational(1);
    m(1, 0) = Rational(5);
    m(1, 1) = Rational(4);
    m(1, 2) = Rational(1);
    std::vector<int> svc = ids(2, 3, 4);
    svc.push_back(5);
    const std::vector<std::vector<int>> chains = line::wf::detect_sequences(m, svc);
    REQUIRE(chains.size() == 2u);
    CHECK(chains[0][0] == 2);
    CHECK(chains[0][1] == 3);
    CHECK(chains[1][0] == 5);
    CHECK(chains[1][1] == 4);
}

TEST_CASE("a workflow with no service-to-service edge has no sequence") {
    const Matrix<Rational> g = parallel_graph<Rational>();
    CHECK(line::wf::detect_sequences(g, ids(2, 3, 5)).empty());
    const line::wf::SequenceStats<Rational> s =
        line::wf::get_sequence_stats<Rational>(std::vector<std::vector<int>>());
    CHECK(s.numSequences == 0u);
    CHECK(s.totalNodes == 0u);
    CHECK(s.avgLength == Rational(0));
}

TEST_CASE("validate_sequence follows the edges, not the JAR's broken Pair") {
    // The JAR returns false here (jline.util.Pair has no equals/hashCode), so
    // this is the one deliberate disagreement; see the header. The port must
    // accept a chain whose consecutive pairs really are edges, and reject one
    // whose pairs are not.
    const Matrix<Rational> g = sequence_graph<Rational>();
    CHECK(line::wf::validate_sequence(ids(2, 3, 4), g));
    CHECK(line::wf::validate_sequence(ids(1, 2, 3), g));
    CHECK_FALSE(line::wf::validate_sequence(ids(2, 4, 3), g));
    CHECK_FALSE(line::wf::validate_sequence(ids(4, 3, 2), g));
    CHECK_FALSE(line::wf::validate_sequence(ids(2), g));
}

// ---------------------------------------------------------------------------
// parallel
// ---------------------------------------------------------------------------

TEST_CASE("detect_parallel finds the three branches of the hand-drawn fork") {
    const Matrix<Rational> g = parallel_graph<Rational>();
    const std::vector<std::vector<int>> p =
        line::wf::detect_parallel(g, ids(2, 3, 5), ids(6), ids(7));
    // JAR reference: parallel = [[2, 3, 5]]
    REQUIRE(p.size() == 1u);
    REQUIRE(p[0].size() == 3u);
    CHECK(p[0][0] == 2);
    CHECK(p[0][1] == 3);
    CHECK(p[0][2] == 5);
    CHECK(line::wf::validate_parallel_pattern(p[0], g, ids(6), ids(7)));
}

TEST_CASE("parallel statistics obey the counting identities exactly") {
    const Matrix<Rational> g = parallel_graph<Rational>();
    const std::vector<std::vector<int>> p =
        line::wf::detect_parallel(g, ids(2, 3, 5), ids(6), ids(7));
    const line::wf::ParallelStats<Rational> s = line::wf::get_parallel_stats<Rational>(p);
    // JAR reference: {avgParallelism=3.0, maxParallelism=3, numPatterns=1,
    //                 totalParallelNodes=3}
    CHECK(s.numPatterns == 1u);
    CHECK(s.totalParallelNodes == 3u);
    CHECK(s.avgParallelism == Rational(3));
    CHECK(s.maxParallelism == 3u);
    CHECK(s.avgParallelism * Rational(static_cast<long>(s.numPatterns)) ==
          Rational(static_cast<long>(s.totalParallelNodes)));
}

TEST_CASE("a single-branch fork-join is not a parallel pattern") {
    // 1 -> 6 -> 2 -> 7 -> 9: only one path, so no pattern and no validation.
    const int e[4][2] = {{1, 6}, {6, 2}, {2, 7}, {7, 9}};
    Matrix<Rational> m(4, 3, Rational(1));
    for (std::size_t i = 0; i < 4; ++i) {
        m(i, 0) = Rational(e[i][0]);
        m(i, 1) = Rational(e[i][1]);
    }
    CHECK(line::wf::detect_parallel(m, ids(2), ids(6), ids(7)).empty());
    CHECK_FALSE(line::wf::validate_parallel_pattern(ids(2), m, ids(6), ids(7)));
    const line::wf::ParallelStats<Rational> s =
        line::wf::get_parallel_stats<Rational>(std::vector<std::vector<int>>());
    CHECK(s.numPatterns == 0u);
    CHECK(s.avgParallelism == Rational(0));
}

// ---------------------------------------------------------------------------
// loops
// ---------------------------------------------------------------------------

TEST_CASE("detect_loops finds the rework loop and only it") {
    const Matrix<Rational> g = loop_graph<Rational>();
    const std::vector<int> loops = line::wf::detect_loops(g, ids(2, 4), ids(8));
    // JAR reference: loops = [2]
    REQUIRE(loops.size() == 1u);
    CHECK(loops[0] == 2);
    CHECK(line::wf::validate_loop_pattern(2, g, ids(8)));
    CHECK_FALSE(line::wf::validate_loop_pattern(4, g, ids(8)));
    // JAR reference: loopProb2 = 0.3, loopProb4 = 0.0
    CHECK(line::wf::get_loop_probability(2, g, ids(8)) == Rational(3, 10));
    CHECK(line::wf::get_loop_probability(4, g, ids(8)) == Rational(0));
}

TEST_CASE("the expected number of loop iterations is exact") {
    // JAR reference: getExpectedLoopIterations(0.3) = 1.4285714285714286 = 10/7
    const line::wf::ExpectedIterations<Rational> e =
        line::wf::get_expected_loop_iterations(Rational(3, 10));
    CHECK_FALSE(e.infinite);
    CHECK(e.value == Rational(10, 7));
    // GEOMETRIC LAW, exact: (1 - p) E = 1.
    CHECK((Rational(1) - Rational(3, 10)) * e.value == Rational(1));
    // p = 0 means the activity runs once.
    CHECK(line::wf::get_expected_loop_iterations(Rational(0)).value == Rational(1));
    // p = 1 diverges; the exact field has no infinity, hence the flag.
    CHECK(line::wf::get_expected_loop_iterations(Rational(1)).infinite);
    CHECK(line::wf::get_expected_loop_iterations(Rational(3, 2)).infinite);
}

TEST_CASE("loop statistics match the JAR, exactly") {
    const Matrix<Rational> g = loop_graph<Rational>();
    const std::vector<int> loops = line::wf::detect_loops(g, ids(2, 4), ids(8));
    const line::wf::LoopStats<Rational> s = line::wf::get_loop_stats(loops, g, ids(8));
    // JAR reference: {avgExpectedIterations=1.4285714285714286,
    //                 avgLoopProbability=0.3, maxExpectedIterations=1.4285714285714286,
    //                 maxLoopProbability=0.3, minLoopProbability=0.3, numLoops=1}
    CHECK(s.numLoops == 1u);
    CHECK(s.avgLoopProbability == Rational(3, 10));
    CHECK(s.maxLoopProbability == Rational(3, 10));
    CHECK(s.minLoopProbability == Rational(3, 10));
    CHECK(s.avgExpectedIterations == Rational(10, 7));
    CHECK(s.maxExpectedIterations == Rational(10, 7));
    // A single loop makes mean, max and min coincide; that identity is the
    // cheapest check that the three are not silently swapped.
    CHECK(s.avgLoopProbability == s.maxLoopProbability);
    CHECK(s.avgExpectedIterations == s.maxExpectedIterations);
}

TEST_CASE("detect_loops finds a service inside a strongly connected component") {
    // 1 -> 2 -> 8(router), 8 -> 2, 8 -> 7(join), 7 -> 3, 3 -> 7, 7 -> 9:
    // node 2 sits on the simple loop and node 3 only in the SCC {3,7}.
    const int e[7][2] = {{1, 2}, {2, 8}, {8, 2}, {8, 7}, {7, 3}, {3, 7}, {7, 9}};
    Matrix<Rational> m(7, 3, Rational(1));
    for (std::size_t i = 0; i < 7; ++i) {
        m(i, 0) = Rational(e[i][0]);
        m(i, 1) = Rational(e[i][1]);
    }
    m(2, 2) = Rational(3, 10);
    m(3, 2) = Rational(7, 10);
    // JAR reference: loopsScc = [2, 3]
    const std::vector<int> with = line::wf::detect_loops(m, ids(2, 3), ids(8), ids(7));
    REQUIRE(with.size() == 2u);
    CHECK(with[0] == 2);
    CHECK(with[1] == 3);
    // Without join nodes the SCC search is skipped, as in the two-argument
    // Java overload, so only the simple loop is reported.
    const std::vector<int> without = line::wf::detect_loops(m, ids(2, 3), ids(8));
    REQUIRE(without.size() == 1u);
    CHECK(without[0] == 2);
}

// ---------------------------------------------------------------------------
// branches
// ---------------------------------------------------------------------------

TEST_CASE("detect_branches finds the hand-drawn choice and its join") {
    const Matrix<Rational> g = branch_graph<Rational>();
    const std::vector<line::wf::BranchPattern<Rational>> b =
        line::wf::detect_branches(g, ids(3, 4), ids(7));
    // JAR reference: nodes=[3, 4] probs=[0.25, 0.75] fork=8 join=7
    REQUIRE(b.size() == 1u);
    REQUIRE(b[0].branchNodes.size() == 2u);
    CHECK(b[0].branchNodes[0] == 3);
    CHECK(b[0].branchNodes[1] == 4);
    CHECK(b[0].probabilities[0] == Rational(1, 4));
    CHECK(b[0].probabilities[1] == Rational(3, 4));
    CHECK(b[0].forkNode == 8);
    CHECK(b[0].hasJoinNode);
    CHECK(b[0].joinNode == 7);
    CHECK(line::wf::validate_branch_pattern(b[0], g));

    // CONSERVATION: a detected branch is a probability distribution, exactly.
    Rational total(0);
    for (std::size_t k = 0; k < b[0].probabilities.size(); ++k) total += b[0].probabilities[k];
    CHECK(total == Rational(1));

    // JAR reference: most = 4 0.75, least = 3 0.25
    const line::wf::BranchAlternative<Rational> most = line::wf::find_most_probable_branch(b[0]);
    const line::wf::BranchAlternative<Rational> least = line::wf::find_least_probable_branch(b[0]);
    CHECK(most.valid);
    CHECK(most.node == 4);
    CHECK(most.probability == Rational(3, 4));
    CHECK(least.valid);
    CHECK(least.node == 3);
    CHECK(least.probability == Rational(1, 4));
}

TEST_CASE("a branch whose probabilities do not sum to one is rejected") {
    Matrix<Rational> g = branch_graph<Rational>();
    g(2, 2) = Rational(1, 4);  // 1/4 + 1/4 = 1/2, outside the 1e-2 slack
    CHECK(line::wf::detect_branches(g, ids(3, 4), ids(7)).empty());

    // The slack itself: 0.25 + 0.746 is 4e-3 short and still admitted.
    Matrix<double> gd = branch_graph<double>();
    gd(2, 2) = 0.746;
    CHECK(line::wf::detect_branches(gd, ids(3, 4), ids(7)).size() == 1u);
    gd(2, 2) = 0.73;  // 2e-2 short, rejected
    CHECK(line::wf::detect_branches(gd, ids(3, 4), ids(7)).empty());
}

TEST_CASE("branch diversity matches the JAR") {
    // JAR reference on the same pattern:
    //   balance=1.3333333333333333, entropy=0.5623351446188083,
    //   gini=0.5, normalizedEntropy=0.8112781244591328
    const Matrix<double> g = branch_graph<double>();
    const std::vector<line::wf::BranchPattern<double>> b =
        line::wf::detect_branches(g, ids(3, 4), ids(7));
    REQUIRE(b.size() == 1u);
    const line::wf::BranchDiversity<double> d = line::wf::calculate_branch_diversity(b[0]);
    // Sums of p log p in double: agreement at the rounding level, the method
    // has no tolerance of its own.
    CHECK(d.entropy == doctest::Approx(0.5623351446188083).epsilon(1e-15));
    CHECK(d.normalizedEntropy == doctest::Approx(0.8112781244591328).epsilon(1e-15));
    CHECK(d.gini == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(d.balance == doctest::Approx(1.3333333333333333).epsilon(1e-15));
    // Entropy is bounded by log(n), so the normalized value is in [0,1].
    CHECK(d.normalizedEntropy <= 1.0);
    CHECK(d.entropy <= std::log(2.0));

    // JAR reference (equiprobable three-way branch): balance=3.0,
    // entropy=1.0986122886681096, gini=0.0, normalizedEntropy=0.9999999999999998
    const int e[8][2] = {{1, 8}, {8, 2}, {8, 3}, {8, 4}, {2, 7}, {3, 7}, {4, 7}, {7, 9}};
    Matrix<double> m(8, 3, 1.0);
    for (std::size_t i = 0; i < 8; ++i) {
        m(i, 0) = e[i][0];
        m(i, 1) = e[i][1];
    }
    m(1, 2) = 1.0 / 3;
    m(2, 2) = 1.0 / 3;
    m(3, 2) = 1.0 / 3;
    const std::vector<line::wf::BranchPattern<double>> b3 =
        line::wf::detect_branches(m, ids(2, 3, 4), ids(7));
    REQUIRE(b3.size() == 1u);
    const line::wf::BranchDiversity<double> d3 = line::wf::calculate_branch_diversity(b3[0]);
    CHECK(d3.entropy == doctest::Approx(1.0986122886681096).epsilon(1e-15));
    CHECK(d3.normalizedEntropy == doctest::Approx(0.9999999999999998).epsilon(1e-15));
    CHECK(d3.gini == doctest::Approx(0.0).epsilon(1e-15));
    CHECK(d3.balance == doctest::Approx(3.0).epsilon(1e-15));
}

TEST_CASE("branch statistics match the JAR") {
    // JAR reference: {avgBalance=1.3333333333333333, avgBranches=2.0,
    //                 avgEntropy=0.5623351446188083, maxBranches=2,
    //                 minBranches=2, numPatterns=1, totalBranchNodes=2}
    const Matrix<double> g = branch_graph<double>();
    const std::vector<line::wf::BranchPattern<double>> b =
        line::wf::detect_branches(g, ids(3, 4), ids(7));
    const line::wf::BranchStats<double> s = line::wf::get_branch_stats(b);
    CHECK(s.numPatterns == 1u);
    CHECK(s.totalBranchNodes == 2u);
    CHECK(s.avgBranches == doctest::Approx(2.0).epsilon(1e-15));
    CHECK(s.maxBranches == 2u);
    CHECK(s.minBranches == 2u);
    CHECK(s.avgEntropy == doctest::Approx(0.5623351446188083).epsilon(1e-15));
    CHECK(s.avgBalance == doctest::Approx(1.3333333333333333).epsilon(1e-15));

    const line::wf::BranchStats<double> empty =
        line::wf::get_branch_stats(std::vector<line::wf::BranchPattern<double>>());
    CHECK(empty.numPatterns == 0u);
    CHECK(empty.avgEntropy == 0.0);
}

TEST_CASE("a one-alternative pattern is refused by branch diversity") {
    line::wf::BranchPattern<double> p;
    p.branchNodes.push_back(3);
    p.probabilities.push_back(1.0);
    p.forkNode = 8;
    // The JAR divides by (n - 1) and returns NaN; the port refuses instead.
    CHECK_THROWS_AS(line::wf::calculate_branch_diversity(p), line::InputError);
}

TEST_CASE("the detectors reject a link matrix without a probability column") {
    Matrix<Rational> m(2, 2, Rational(1));
    CHECK_THROWS_AS(line::wf::detect_sequences(m, ids(2, 3)), line::InputError);
    CHECK_THROWS_AS(line::wf::detect_parallel(m, ids(2), ids(6), ids(7)), line::InputError);
    CHECK_THROWS_AS(line::wf::detect_loops(m, ids(2), ids(8)), line::InputError);
    CHECK_THROWS_AS(line::wf::detect_branches(m, ids(3, 4), ids(7)), line::InputError);
}

TEST_CASE("the detectors instantiate at high precision") {
    // Same graphs at 50 decimal digits: the combinatorics are identical and
    // the loop algebra is now accurate far past double, so (1 - p) E = 1 holds
    // to the working precision rather than to 1e-16.
    using line::Real50;
    const Matrix<Real50> g = loop_graph<Real50>();
    const std::vector<int> loops = line::wf::detect_loops(g, ids(2, 4), ids(8));
    REQUIRE(loops.size() == 1u);
    CHECK(loops[0] == 2);
    const Real50 p = line::wf::get_loop_probability(2, g, ids(8));
    const line::wf::ExpectedIterations<Real50> e = line::wf::get_expected_loop_iterations(p);
    CHECK_FALSE(e.infinite);
    const Real50 residual = Real50((Real50(1) - p) * e.value - Real50(1));
    CHECK(static_cast<double>(line::num_abs(residual)) < 1e-45);

    const Matrix<Real50> s = sequence_graph<Real50>();
    const std::vector<std::vector<int>> chains = line::wf::detect_sequences(s, ids(2, 3, 4));
    REQUIRE(chains.size() == 1u);
    CHECK(chains[0][2] == 4);
    CHECK(line::wf::validate_sequence(ids(2, 3, 4), s));

    // The entropy metrics are gated on transcendental arithmetic, which Real50
    // has: the equiprobable two-way branch must give exactly log 2 there.
    const Matrix<Real50> b = branch_graph<Real50>();
    Matrix<Real50> half(b);
    half(1, 2) = Real50(1) / Real50(2);
    half(2, 2) = Real50(1) / Real50(2);
    const std::vector<line::wf::BranchPattern<Real50>> bp =
        line::wf::detect_branches(half, ids(3, 4), ids(7));
    REQUIRE(bp.size() == 1u);
    const line::wf::BranchDiversity<Real50> d = line::wf::calculate_branch_diversity(bp[0]);
    const Real50 log2 = line::wf::detail::num_log(Real50(2));
    CHECK(static_cast<double>(line::num_abs(Real50(d.entropy - log2))) < 1e-45);
    CHECK(static_cast<double>(line::num_abs(Real50(d.normalizedEntropy - Real50(1)))) < 1e-45);
}
