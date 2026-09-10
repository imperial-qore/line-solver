/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_BRANCH_DETECTOR_H
#define LINE_API_WF_WF_BRANCH_DETECTOR_H

/**
 * Branch (probabilistic choice) pattern detection in a workflow network.
 *
 * Templated port of jar/src/main/java/jline/api/wf/Wf_branch_detector.java (no
 * MATLAB counterpart, so the JAR is the reference). A branch point is a node
 * with more than one outgoing edge of which at least two lead to a service
 * node; the branch is accepted when the probabilities of those service edges
 * sum to one within 1e-2, and it is annotated with the first join node
 * reachable from every branch alternative.
 *
 * The six public methods of the Java class are kept 1:1: detect_branches,
 * validate_branch_pattern, calculate_branch_diversity, get_branch_stats,
 * find_most_probable_branch, find_least_probable_branch.
 *
 * ARITHMETIC. detect_branches, validate_branch_pattern and the two extremal
 * queries only add and compare probabilities, so they are finite field
 * computations and instantiate at exact arithmetic; the 1e-2 slack on the
 * branch probability sum is a structural admission threshold inherited from
 * the reference, not a rounding allowance, so it is kept in every
 * instantiation. calculate_branch_diversity and get_branch_stats compute the
 * Shannon entropy of the branch probabilities and therefore need log; both
 * carry the transcendental gate.
 *
 * REFERENCE DEFECT (JAR): calculateBranchDiversity divides the Gini sum by
 * (n-1), so a one-alternative pattern yields a division by zero (NaN in Java).
 * The port rejects n < 2 with an InputError instead of returning a NaN;
 * detect_branches never emits such a pattern, so no detected pattern is
 * affected.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "line/api/wf/wf_link_matrix.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** Mirrors the Java BranchPattern. */
template <class T>
struct BranchPattern {
    std::vector<int> branchNodes;
    std::vector<T> probabilities;
    int forkNode = -1;
    bool hasJoinNode = false;  ///< the Java Integer may be null
    int joinNode = -1;
};

/** Mirrors the Java calculateBranchDiversity map. */
template <class T>
struct BranchDiversity {
    T entropy = num_traits<T>::from_int(0);
    T normalizedEntropy = num_traits<T>::from_int(0);
    T gini = num_traits<T>::from_int(0);
    T balance = num_traits<T>::from_int(0);
};

/** Mirrors the Java getBranchStats map. */
template <class T>
struct BranchStats {
    std::size_t numPatterns = 0;
    std::size_t totalBranchNodes = 0;
    T avgBranches = num_traits<T>::from_int(0);
    std::size_t maxBranches = 0;
    std::size_t minBranches = 0;
    T avgEntropy = num_traits<T>::from_int(0);
    T avgBalance = num_traits<T>::from_int(0);
};

/** One alternative of a branch: the node and its probability. */
template <class T>
struct BranchAlternative {
    bool valid = false;
    int node = -1;
    T probability = num_traits<T>::from_int(0);
};

namespace detail {

/** Nodes reachable from startNode; the search does not continue past stopSet. */
template <class T>
std::set<int> wf_reachable_stop(int startNode,
                                const std::map<int, std::vector<std::pair<int, T>>>& adj,
                                const std::set<int>& stopSet) {
    std::set<int> reachable;
    std::set<int> visited;
    std::vector<int> queue;
    std::size_t head = 0;
    queue.push_back(startNode);
    while (head < queue.size()) {
        const int current = queue[head++];
        if (!visited.insert(current).second) continue;
        typename std::map<int, std::vector<std::pair<int, T>>>::const_iterator it =
            adj.find(current);
        if (it == adj.end()) continue;
        for (std::size_t k = 0; k < it->second.size(); ++k) {
            const int nb = it->second[k].first;
            reachable.insert(nb);
            if (!stopSet.count(nb)) queue.push_back(nb);
        }
    }
    return reachable;
}

/** log(v), ADL-visible so a non-double T can supply its own overload. */
template <class T>
inline T num_log(const T& v) {
    using std::log;
    return T(log(v));
}

}  // namespace detail

/**
 * @param linkMatrix   (nedges x 3) edge list
 * @param serviceNodes ids of the service nodes
 * @param joinNodes    ids of the join nodes
 */
template <class T>
std::vector<BranchPattern<T>> detect_branches(const Matrix<T>& linkMatrix,
                                              const std::vector<int>& serviceNodes,
                                              const std::vector<int>& joinNodes) {
    detail::wf_check(linkMatrix);
    const std::set<int> serviceSet(serviceNodes.begin(), serviceNodes.end());
    const std::set<int> joinSet(joinNodes.begin(), joinNodes.end());
    const std::map<int, std::vector<std::pair<int, T>>> adj = detail::wf_adjacency_prob(linkMatrix);
    const T one = num_traits<T>::from_int(1);
    const T slack = num_traits<T>::from_rational(1, 100);

    std::vector<BranchPattern<T>> patterns;
    for (typename std::map<int, std::vector<std::pair<int, T>>>::const_iterator it = adj.begin();
         it != adj.end(); ++it) {
        if (it->second.size() <= 1) continue;
        std::vector<std::pair<int, T>> targets;
        for (std::size_t k = 0; k < it->second.size(); ++k)
            if (serviceSet.count(it->second[k].first)) targets.push_back(it->second[k]);
        if (targets.size() < 2) continue;

        T total = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < targets.size(); ++k) total += targets[k].second;
        if (num_abs(T(total - one)) > slack) continue;

        BranchPattern<T> p;
        p.forkNode = it->first;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            p.branchNodes.push_back(targets[k].first);
            p.probabilities.push_back(targets[k].second);
        }

        // Common join point: intersect the reachable sets of the alternatives,
        // prefer a join node, else the smallest common successor.
        std::set<int> common = detail::wf_reachable_stop(p.branchNodes[0], adj, joinSet);
        for (std::size_t k = 1; k < p.branchNodes.size(); ++k) {
            const std::set<int> r = detail::wf_reachable_stop(p.branchNodes[k], adj, joinSet);
            std::set<int> inter;
            std::set_intersection(common.begin(), common.end(), r.begin(), r.end(),
                                  std::inserter(inter, inter.begin()));
            common.swap(inter);
        }
        std::set<int> joinPoints;
        std::set_intersection(common.begin(), common.end(), joinSet.begin(), joinSet.end(),
                              std::inserter(joinPoints, joinPoints.begin()));
        if (!joinPoints.empty()) {
            p.hasJoinNode = true;
            p.joinNode = *joinPoints.begin();
        } else if (!common.empty()) {
            p.hasJoinNode = true;
            p.joinNode = *common.begin();
        }

        patterns.push_back(p);
    }
    return patterns;
}

/** Probabilities sum to one within 1e-2 and every alternative is a fork successor. */
template <class T>
bool validate_branch_pattern(const BranchPattern<T>& pattern, const Matrix<T>& linkMatrix) {
    detail::wf_check(linkMatrix);
    const T one = num_traits<T>::from_int(1);
    T total = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < pattern.probabilities.size(); ++k) total += pattern.probabilities[k];
    if (num_abs(T(total - one)) > num_traits<T>::from_rational(1, 100)) return false;
    if (pattern.forkNode < 0) return false;

    const std::map<int, std::vector<std::pair<int, T>>> adj = detail::wf_adjacency_prob(linkMatrix);
    typename std::map<int, std::vector<std::pair<int, T>>>::const_iterator it =
        adj.find(pattern.forkNode);
    if (it == adj.end()) return false;
    std::set<int> forkTargets;
    for (std::size_t k = 0; k < it->second.size(); ++k) forkTargets.insert(it->second[k].first);
    for (std::size_t k = 0; k < pattern.branchNodes.size(); ++k)
        if (!forkTargets.count(pattern.branchNodes[k])) return false;
    return true;
}

/**
 * Shannon entropy of the branch probabilities, the same entropy normalized by
 * log(n), the Gini coefficient of the probability vector, and the reciprocal
 * of the largest probability.
 */
template <class T>
BranchDiversity<T> calculate_branch_diversity(const BranchPattern<T>& pattern) {
    static_assert(num_traits<T>::has_transcendental,
                  "calculate_branch_diversity requires transcendental arithmetic: the entropy of "
                  "the branch probabilities is a sum of p log p");
    const std::vector<T>& probs = pattern.probabilities;
    if (probs.size() < 2)
        throw InputError(
            "calculate_branch_diversity: needs at least two alternatives, the Gini coefficient "
            "divides by (n - 1)");
    const std::size_t n = probs.size();
    const T zero = num_traits<T>::from_int(0);

    BranchDiversity<T> d;
    for (std::size_t k = 0; k < n; ++k)
        if (probs[k] > zero) d.entropy -= T(probs[k] * detail::num_log(probs[k]));
    d.normalizedEntropy =
        d.entropy / detail::num_log(num_traits<T>::from_int(static_cast<long>(n)));

    std::vector<T> sorted(probs);
    std::sort(sorted.begin(), sorted.end());
    T sumProbs = zero;
    for (std::size_t k = 0; k < n; ++k) sumProbs += probs[k];
    T gini = zero;
    for (std::size_t k = 0; k < n; ++k)
        gini += num_traits<T>::from_int(static_cast<long>(2 * (k + 1)) -
                                        static_cast<long>(n) - 1) *
                sorted[k];
    const T giniDen = T(num_traits<T>::from_int(static_cast<long>(n) - 1) * sumProbs);
    gini /= giniDen;
    d.gini = num_abs(gini);

    T maxProb = zero;
    for (std::size_t k = 0; k < n; ++k)
        if (probs[k] > maxProb) maxProb = probs[k];
    d.balance = num_traits<T>::from_int(1) / (maxProb == zero ? num_traits<T>::from_int(1) : maxProb);
    return d;
}

/** Count, total, mean/max/min alternatives, and the mean entropy and balance. */
template <class T>
BranchStats<T> get_branch_stats(const std::vector<BranchPattern<T>>& patterns) {
    static_assert(num_traits<T>::has_transcendental,
                  "get_branch_stats requires transcendental arithmetic: it averages the entropy "
                  "returned by calculate_branch_diversity");
    BranchStats<T> stats;
    stats.numPatterns = patterns.size();
    std::size_t total = 0;
    for (std::size_t i = 0; i < patterns.size(); ++i) total += patterns[i].branchNodes.size();
    stats.totalBranchNodes = total;
    if (patterns.empty()) return stats;

    stats.maxBranches = patterns[0].branchNodes.size();
    stats.minBranches = patterns[0].branchNodes.size();
    for (std::size_t i = 1; i < patterns.size(); ++i) {
        const std::size_t sz = patterns[i].branchNodes.size();
        if (sz > stats.maxBranches) stats.maxBranches = sz;
        if (sz < stats.minBranches) stats.minBranches = sz;
    }
    const T np = num_traits<T>::from_int(static_cast<long>(patterns.size()));
    stats.avgBranches = num_traits<T>::from_int(static_cast<long>(total)) / np;

    T sumE = num_traits<T>::from_int(0);
    T sumB = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < patterns.size(); ++i) {
        const BranchDiversity<T> d = calculate_branch_diversity(patterns[i]);
        sumE += d.entropy;
        sumB += d.balance;
    }
    stats.avgEntropy = sumE / np;
    stats.avgBalance = sumB / np;
    return stats;
}

/** The alternative with the largest probability. */
template <class T>
BranchAlternative<T> find_most_probable_branch(const BranchPattern<T>& pattern) {
    BranchAlternative<T> r;
    if (pattern.branchNodes.empty() || pattern.probabilities.empty()) return r;
    std::size_t best = 0;
    for (std::size_t k = 1; k < pattern.probabilities.size(); ++k)
        if (pattern.probabilities[k] > pattern.probabilities[best]) best = k;
    r.valid = true;
    r.node = pattern.branchNodes[best];
    r.probability = pattern.probabilities[best];
    return r;
}

/** The alternative with the smallest probability. */
template <class T>
BranchAlternative<T> find_least_probable_branch(const BranchPattern<T>& pattern) {
    BranchAlternative<T> r;
    if (pattern.branchNodes.empty() || pattern.probabilities.empty()) return r;
    std::size_t best = 0;
    for (std::size_t k = 1; k < pattern.probabilities.size(); ++k)
        if (pattern.probabilities[k] < pattern.probabilities[best]) best = k;
    r.valid = true;
    r.node = pattern.branchNodes[best];
    r.probability = pattern.probabilities[best];
    return r;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_BRANCH_DETECTOR_H
