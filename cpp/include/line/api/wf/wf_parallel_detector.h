/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_PARALLEL_DETECTOR_H
#define LINE_API_WF_WF_PARALLEL_DETECTOR_H

/**
 * Parallel (fork-join) pattern detection in a workflow network.
 *
 * Templated port of jar/src/main/java/jline/api/wf/Wf_parallel_detector.java
 * (no MATLAB counterpart, so the JAR is the reference). A fork f and a join j
 * form a pair when a breadth-first path count from f, forbidden to pass
 * through any other fork or join, reaches j along more than one path; the
 * parallel branches of that pair are then the service nodes that are both
 * reachable from f without passing through j and able to reach j without
 * passing through f. A pair contributes a pattern only when it has at least
 * two such service nodes.
 *
 * detect_parallel, validate_parallel_pattern and get_parallel_stats are the
 * three public methods of the Java class, kept 1:1.
 *
 * Pure graph traversal; the only arithmetic is the mean branch count, one
 * division. Finite field computation, instantiates at exact arithmetic, no
 * transcendental gate.
 */

#include <cstddef>
#include <deque>
#include <map>
#include <set>
#include <vector>

#include "line/api/wf/wf_link_matrix.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace wf {

/** Mirrors the Java getParallelStats map. */
template <class T>
struct ParallelStats {
    std::size_t numPatterns = 0;
    std::size_t totalParallelNodes = 0;
    T avgParallelism = num_traits<T>::from_int(0);
    std::size_t maxParallelism = 0;
};

namespace detail {

/**
 * More than one fork-to-join path avoiding every other fork and join. The path
 * count is accumulated in breadth-first order, as in the reference; it is a
 * lower bound on the true path count, not the exact one, but the test is only
 * "more than one".
 */
inline bool wf_valid_fork_join_pair(int fork, int join, const std::map<int, std::vector<int>>& adj,
                                    const std::set<int>& forkSet, const std::set<int>& joinSet) {
    std::deque<int> queue;
    std::set<int> visited;
    std::map<int, long> pathCount;
    queue.push_back(fork);
    pathCount[fork] = 1;

    while (!queue.empty()) {
        const int current = queue.front();
        queue.pop_front();
        if (!visited.insert(current).second) continue;
        std::map<int, std::vector<int>>::const_iterator it = adj.find(current);
        if (it == adj.end()) continue;
        for (std::size_t k = 0; k < it->second.size(); ++k) {
            const int nb = it->second[k];
            if (nb == join) {
                pathCount[join] += pathCount[current];
            } else if (!visited.count(nb) && !forkSet.count(nb) && !joinSet.count(nb)) {
                queue.push_back(nb);
                pathCount[nb] += pathCount[current];
            }
        }
    }
    std::map<int, long>::const_iterator jt = pathCount.find(join);
    return jt != pathCount.end() && jt->second > 1;
}

/** Nodes reachable from startNode without entering endNode. */
template <class T>
std::set<int> wf_reachable(const Matrix<T>& linkMatrix, int startNode, int endNode) {
    std::set<int> reachable;
    std::set<int> visited;
    std::deque<int> queue;
    queue.push_back(startNode);
    while (!queue.empty()) {
        const int current = queue.front();
        queue.pop_front();
        if (visited.count(current) || current == endNode) continue;
        visited.insert(current);
        for (std::size_t i = 0; i < linkMatrix.rows(); ++i) {
            const int s = wf_id(linkMatrix, i, 0);
            const int e = wf_id(linkMatrix, i, 1);
            if (s == current && e != endNode) {
                reachable.insert(e);
                queue.push_back(e);
            }
        }
    }
    return reachable;
}

/** Nodes that reach targetNode without passing through startNode. */
template <class T>
std::set<int> wf_can_reach(const Matrix<T>& linkMatrix, int targetNode, int startNode) {
    const std::map<int, std::vector<int>> radj = wf_reverse_adjacency(linkMatrix);
    std::set<int> canReach;
    std::set<int> visited;
    std::deque<int> queue;
    queue.push_back(targetNode);
    while (!queue.empty()) {
        const int current = queue.front();
        queue.pop_front();
        if (visited.count(current) || current == startNode) continue;
        visited.insert(current);
        std::map<int, std::vector<int>>::const_iterator it = radj.find(current);
        if (it == radj.end()) continue;
        for (std::size_t k = 0; k < it->second.size(); ++k) {
            const int pred = it->second[k];
            if (pred != startNode) {
                canReach.insert(pred);
                queue.push_back(pred);
            }
        }
    }
    return canReach;
}

}  // namespace detail

/**
 * @param linkMatrix   (nedges x 3) edge list
 * @param serviceNodes ids of the service nodes
 * @param forkNodes    ids of the fork nodes
 * @param joinNodes    ids of the join nodes
 * @return one list of parallel service nodes per detected fork-join pair
 */
template <class T>
std::vector<std::vector<int>> detect_parallel(const Matrix<T>& linkMatrix,
                                              const std::vector<int>& serviceNodes,
                                              const std::vector<int>& forkNodes,
                                              const std::vector<int>& joinNodes) {
    detail::wf_check(linkMatrix);
    const std::set<int> forkSet(forkNodes.begin(), forkNodes.end());
    const std::set<int> joinSet(joinNodes.begin(), joinNodes.end());
    const std::set<int> serviceSet(serviceNodes.begin(), serviceNodes.end());
    const std::map<int, std::vector<int>> adj = detail::wf_adjacency(linkMatrix);

    std::vector<std::vector<int>> patterns;
    for (std::size_t a = 0; a < forkNodes.size(); ++a) {
        for (std::size_t b = 0; b < joinNodes.size(); ++b) {
            const int fork = forkNodes[a];
            const int join = joinNodes[b];
            if (!detail::wf_valid_fork_join_pair(fork, join, adj, forkSet, joinSet)) continue;

            const std::set<int> fromFork = detail::wf_reachable(linkMatrix, fork, join);
            const std::set<int> toJoin = detail::wf_can_reach(linkMatrix, join, fork);
            std::vector<int> parallelServices;
            for (std::set<int>::const_iterator it = fromFork.begin(); it != fromFork.end(); ++it)
                if (toJoin.count(*it) && serviceSet.count(*it)) parallelServices.push_back(*it);
            if (parallelServices.size() > 1) patterns.push_back(parallelServices);
        }
    }
    return patterns;
}

/**
 * A pattern is valid when its nodes have exactly one common fork predecessor
 * and exactly one common join successor.
 */
template <class T>
bool validate_parallel_pattern(const std::vector<int>& pattern, const Matrix<T>& linkMatrix,
                               const std::vector<int>& forkNodes,
                               const std::vector<int>& joinNodes) {
    detail::wf_check(linkMatrix);
    if (pattern.size() < 2) return false;
    const std::set<int> forkSet(forkNodes.begin(), forkNodes.end());
    const std::set<int> joinSet(joinNodes.begin(), joinNodes.end());

    std::set<int> sources, targets;
    for (std::size_t p = 0; p < pattern.size(); ++p) {
        for (std::size_t i = 0; i < linkMatrix.rows(); ++i) {
            const int s = detail::wf_id(linkMatrix, i, 0);
            const int e = detail::wf_id(linkMatrix, i, 1);
            if (e == pattern[p] && forkSet.count(s)) sources.insert(s);
            if (s == pattern[p] && joinSet.count(e)) targets.insert(e);
        }
    }
    return sources.size() == 1 && targets.size() == 1;
}

/** Count, total, mean and maximum degree of parallelism. */
template <class T>
ParallelStats<T> get_parallel_stats(const std::vector<std::vector<int>>& patterns) {
    ParallelStats<T> stats;
    stats.numPatterns = patterns.size();
    std::size_t total = 0;
    std::size_t mx = 0;
    for (std::size_t i = 0; i < patterns.size(); ++i) {
        total += patterns[i].size();
        if (patterns[i].size() > mx) mx = patterns[i].size();
    }
    stats.totalParallelNodes = total;
    stats.maxParallelism = mx;
    if (!patterns.empty())
        stats.avgParallelism = num_traits<T>::from_int(static_cast<long>(total)) /
                               num_traits<T>::from_int(static_cast<long>(patterns.size()));
    return stats;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_PARALLEL_DETECTOR_H
