/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_WF_WF_LOOP_DETECTOR_H
#define LINE_API_WF_WF_LOOP_DETECTOR_H

/**
 * Loop pattern detection in a workflow network.
 *
 * Templated port of jar/src/main/java/jline/api/wf/Wf_loop_detector.java (no
 * MATLAB counterpart, so the JAR is the reference). Two mechanisms:
 *
 *  - a SIMPLE loop is a service node with an edge to a router that has an edge
 *    back to it, the two-hop rework loop that a workflow model builds for a
 *    "repeat the activity with probability p" construct;
 *  - a COMPLEX loop is a service node inside a strongly connected component of
 *    more than one node that also contains a router or a join, found by
 *    Tarjan's algorithm. Only looked for when join nodes are supplied, exactly
 *    as in the reference.
 *
 * The five public methods of the Java class are kept 1:1: detect_loops,
 * get_loop_probability, validate_loop_pattern, get_expected_loop_iterations,
 * get_loop_stats.
 *
 * Traversal plus, in the statistics, the geometric mean number of iterations
 * 1/(1-p) and averages: sums, one division each. Finite field computation, so
 * this instantiates at exact arithmetic and the expected iteration count of a
 * rational loop probability is an exact rational. No transcendental gate.
 *
 * The reference returns GlobalConstants.Inf when the loop probability reaches
 * one. An exact rational field has no infinity, so get_expected_loop_iterations
 * returns ExpectedIterations with an explicit `infinite` flag instead of a
 * sentinel value; get_loop_stats drops the infinite entries, which is what the
 * reference's Double.isFinite filter does.
 */

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

/** 1/(1-p), with the p >= 1 divergence reported rather than encoded. */
template <class T>
struct ExpectedIterations {
    bool infinite = false;
    T value = num_traits<T>::from_int(0);
};

/** Mirrors the Java getLoopStats map. */
template <class T>
struct LoopStats {
    std::size_t numLoops = 0;
    T avgLoopProbability = num_traits<T>::from_int(0);
    T maxLoopProbability = num_traits<T>::from_int(0);
    T minLoopProbability = num_traits<T>::from_int(0);
    T avgExpectedIterations = num_traits<T>::from_int(0);
    T maxExpectedIterations = num_traits<T>::from_int(0);
};

namespace detail {

/** service -> router -> service, the two-hop rework loop. */
template <class T>
bool wf_in_simple_loop(int serviceNode, const std::map<int, std::vector<std::pair<int, T>>>& adj,
                       const std::set<int>& routerSet) {
    typename std::map<int, std::vector<std::pair<int, T>>>::const_iterator it = adj.find(serviceNode);
    if (it == adj.end()) return false;
    for (std::size_t a = 0; a < it->second.size(); ++a) {
        const int router = it->second[a].first;
        if (!routerSet.count(router)) continue;
        typename std::map<int, std::vector<std::pair<int, T>>>::const_iterator jt = adj.find(router);
        if (jt == adj.end()) continue;
        for (std::size_t b = 0; b < jt->second.size(); ++b)
            if (jt->second[b].first == serviceNode) return true;
    }
    return false;
}

struct TarjanState {
    std::map<int, int> index;
    std::map<int, int> lowlink;
    std::set<int> onStack;
    std::vector<int> stack;
    std::vector<std::vector<int>> sccs;
    int counter = 0;
};

inline void wf_strong_connect(int node, const std::map<int, std::set<int>>& graph, TarjanState& st) {
    st.index[node] = st.counter;
    st.lowlink[node] = st.counter;
    st.counter++;
    st.stack.push_back(node);
    st.onStack.insert(node);

    std::map<int, std::set<int>>::const_iterator it = graph.find(node);
    if (it != graph.end()) {
        for (std::set<int>::const_iterator nb = it->second.begin(); nb != it->second.end(); ++nb) {
            if (!st.index.count(*nb)) {
                wf_strong_connect(*nb, graph, st);
                if (st.lowlink[*nb] < st.lowlink[node]) st.lowlink[node] = st.lowlink[*nb];
            } else if (st.onStack.count(*nb)) {
                if (st.index[*nb] < st.lowlink[node]) st.lowlink[node] = st.index[*nb];
            }
        }
    }

    if (st.lowlink[node] == st.index[node]) {
        std::vector<int> scc;
        int w;
        do {
            w = st.stack.back();
            st.stack.pop_back();
            st.onStack.erase(w);
            scc.push_back(w);
        } while (w != node);
        st.sccs.push_back(scc);
    }
}

}  // namespace detail

/**
 * @param linkMatrix   (nedges x 3) edge list
 * @param serviceNodes ids of the service nodes
 * @param routerNodes  ids of the router nodes
 * @param joinNodes    ids of the join nodes; empty disables the SCC search,
 *                     matching the two-argument Java overload
 * @return the service nodes that sit on a loop, in detection order, distinct
 */
template <class T>
std::vector<int> detect_loops(const Matrix<T>& linkMatrix, const std::vector<int>& serviceNodes,
                              const std::vector<int>& routerNodes,
                              const std::vector<int>& joinNodes = std::vector<int>()) {
    detail::wf_check(linkMatrix);
    const std::set<int> routerSet(routerNodes.begin(), routerNodes.end());
    const std::map<int, std::vector<std::pair<int, T>>> adj = detail::wf_adjacency_prob(linkMatrix);

    std::vector<int> loopNodes;
    for (std::size_t i = 0; i < serviceNodes.size(); ++i)
        if (detail::wf_in_simple_loop(serviceNodes[i], adj, routerSet))
            loopNodes.push_back(serviceNodes[i]);

    if (!joinNodes.empty()) {
        const std::set<int> serviceSet(serviceNodes.begin(), serviceNodes.end());
        const std::set<int> joinSet(joinNodes.begin(), joinNodes.end());

        std::map<int, std::set<int>> graph;
        for (std::size_t i = 0; i < linkMatrix.rows(); ++i)
            graph[detail::wf_id(linkMatrix, i, 0)].insert(detail::wf_id(linkMatrix, i, 1));

        detail::TarjanState st;
        for (std::map<int, std::set<int>>::const_iterator it = graph.begin(); it != graph.end(); ++it)
            if (!st.index.count(it->first)) detail::wf_strong_connect(it->first, graph, st);

        for (std::size_t s = 0; s < st.sccs.size(); ++s) {
            if (st.sccs[s].size() <= 1) continue;
            std::vector<int> inScc;
            bool hasRouterOrJoin = false;
            for (std::size_t k = 0; k < st.sccs[s].size(); ++k) {
                const int n = st.sccs[s][k];
                if (serviceSet.count(n)) inScc.push_back(n);
                if (routerSet.count(n) || joinSet.count(n)) hasRouterOrJoin = true;
            }
            if (!inScc.empty() && hasRouterOrJoin)
                loopNodes.insert(loopNodes.end(), inScc.begin(), inScc.end());
        }
    }

    std::vector<int> distinct;
    std::set<int> seen;
    for (std::size_t i = 0; i < loopNodes.size(); ++i)
        if (seen.insert(loopNodes[i]).second) distinct.push_back(loopNodes[i]);
    return distinct;
}

/**
 * Probability on the router-to-service edge that closes the loop, 0 when the
 * node is not on a simple loop.
 */
template <class T>
T get_loop_probability(int serviceNode, const Matrix<T>& linkMatrix,
                       const std::vector<int>& routerNodes) {
    detail::wf_check(linkMatrix);
    const std::set<int> routerSet(routerNodes.begin(), routerNodes.end());
    for (std::size_t i = 0; i < linkMatrix.rows(); ++i) {
        const int start = detail::wf_id(linkMatrix, i, 0);
        const int end = detail::wf_id(linkMatrix, i, 1);
        if (start != serviceNode || !routerSet.count(end)) continue;
        for (std::size_t j = 0; j < linkMatrix.rows(); ++j)
            if (detail::wf_id(linkMatrix, j, 0) == end &&
                detail::wf_id(linkMatrix, j, 1) == serviceNode)
                return linkMatrix(j, 2);
    }
    return num_traits<T>::from_int(0);
}

/** True when the node still has the service -> router -> service structure. */
template <class T>
bool validate_loop_pattern(int loopNode, const Matrix<T>& linkMatrix,
                           const std::vector<int>& routerNodes) {
    detail::wf_check(linkMatrix);
    const std::set<int> routerSet(routerNodes.begin(), routerNodes.end());
    return detail::wf_in_simple_loop(loopNode, detail::wf_adjacency_prob(linkMatrix), routerSet);
}

/** Mean number of visits of a geometric loop, 1/(1-p). */
template <class T>
ExpectedIterations<T> get_expected_loop_iterations(const T& loopProbability) {
    const T one = num_traits<T>::from_int(1);
    ExpectedIterations<T> r;
    if (loopProbability >= one) {
        r.infinite = true;
        return r;
    }
    r.value = one / (one - loopProbability);
    return r;
}

/** Count and moments of the loop probabilities and iteration counts. */
template <class T>
LoopStats<T> get_loop_stats(const std::vector<int>& loopNodes, const Matrix<T>& linkMatrix,
                            const std::vector<int>& routerNodes) {
    LoopStats<T> stats;
    stats.numLoops = loopNodes.size();

    std::vector<T> probs;
    for (std::size_t i = 0; i < loopNodes.size(); ++i)
        probs.push_back(get_loop_probability(loopNodes[i], linkMatrix, routerNodes));

    if (!probs.empty()) {
        T sum = num_traits<T>::from_int(0);
        T mx = probs[0];
        T mn = probs[0];
        for (std::size_t i = 0; i < probs.size(); ++i) {
            sum += probs[i];
            if (probs[i] > mx) mx = probs[i];
            if (probs[i] < mn) mn = probs[i];
        }
        stats.avgLoopProbability = sum / num_traits<T>::from_int(static_cast<long>(probs.size()));
        stats.maxLoopProbability = mx;
        stats.minLoopProbability = mn;
    }

    std::vector<T> iters;
    for (std::size_t i = 0; i < probs.size(); ++i) {
        const ExpectedIterations<T> e = get_expected_loop_iterations(probs[i]);
        if (!e.infinite) iters.push_back(e.value);
    }
    if (!iters.empty()) {
        T sum = num_traits<T>::from_int(0);
        T mx = iters[0];
        for (std::size_t i = 0; i < iters.size(); ++i) {
            sum += iters[i];
            if (iters[i] > mx) mx = iters[i];
        }
        stats.avgExpectedIterations = sum / num_traits<T>::from_int(static_cast<long>(iters.size()));
        stats.maxExpectedIterations = mx;
    }
    return stats;
}

}  // namespace wf
}  // namespace line

#endif  // LINE_API_WF_WF_LOOP_DETECTOR_H
