/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_GAMMA_H
#define LINE_API_CACHE_GAMMA_H

/**
 * Access factors of a multi-list cache whose lists form a general access GRAPH.
 *
 * Templated port of jar/src/main/java/jline/api/cache/Cache_gamma.java. MATLAB
 * has no counterpart: it carries only cache_gamma_lp, the specialization to a
 * tree ("linear path"), which recovers the path by walking the unique parent
 * relation and rejects a node with two parents.
 *
 * Here the structure is only required to be reachable. The path is the
 * BREADTH-FIRST shortest path in the access graph of item i, so a node with
 * several parents is admissible and the first shortest path found in node order
 * is the one taken. Along that path,
 *
 *   gamma(i,j) = (sum_v lambda(v,i,0)) prod_{edges (a,b)} sum_v lambda(v,i,a) R{v,i}(a,b)
 *
 * THREE DIVERGENCES FROM cache_gamma_lp, all faithful to the JAR and all of
 * them changing the number, so a caller must not treat the two as substitutes:
 *   - THE DESTINATION IS NODE j, NOT NODE j+1. cache_gamma_lp walks to node l+1
 *     for column l, node 0 being the miss list; this walks to node j. Column 0
 *     therefore has the trivial one-node path and carries NO edge factor at all,
 *     where the tree version carries the miss-to-first-list edge.
 *   - the leading factor is the aggregate miss-node request rate sum_v
 *     lambda(v,i,0), whereas the tree version starts the product at one;
 *   - each edge factor reads lambda(v,i,a) at the SOURCE node a alone, whereas
 *     the tree version sums lambda(v,i,t) over every t <= a.
 * The first of these is hard to read as anything but an off-by-one against the
 * shared meaning of gamma. It is reproduced rather than corrected because this
 * routine exists only in the JAR -- MATLAB has no cache_gamma to arbitrate, and
 * no solver in any codebase calls it (only a JUnit import does), so "fixing" it
 * would leave the C++ disagreeing with the sole reference that defines it. Use
 * cache_gamma_lp for the access factors a cache solver consumes.
 *
 * An unreachable node gives gamma(i,j) = 0, which unlike the tree version is a
 * REACHABLE branch here: the BFS genuinely returns no path.
 *
 * THE GRAPH IS READ FROM USER 0 ONLY. The JAR takes R.get(0).get(i) for the
 * adjacency and then sums the per-user rates along that one path, so a model
 * whose users route an item differently is analysed on the first user's graph.
 * Reproduced rather than corrected, since changing it would silently move the
 * answer for every such model.
 *
 * Arithmetic: EXACT-CAPABLE. Sums and products only; the BFS is pure integer
 * bookkeeping and tests adjacency against zero, which is exact in any T.
 */

#include <cstddef>
#include <deque>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_gamma, mirroring the JAR's Ret.cacheGamma. */
template <class T>
struct CacheGammaGraphResult {
    Matrix<T> gamma;  ///< (n x h) access factors
    std::size_t u;    ///< number of user streams
    std::size_t n;    ///< number of items
    std::size_t h;    ///< number of lists
};

namespace detail {

/**
 * Breadth-first shortest path from source to destination over the nonzero
 * entries of an adjacency matrix; empty when the destination is unreachable.
 */
template <class T>
std::vector<std::size_t> cache_bfs_path(const Matrix<T>& adjacency, std::size_t source,
                                        std::size_t destination) {
    const std::size_t n = adjacency.rows();
    if (source >= n || destination >= n) return std::vector<std::size_t>();
    const T zero = num_traits<T>::from_int(0);
    std::vector<bool> visited(n, false);
    std::vector<int> parent(n, -1);
    std::deque<std::size_t> queue;
    queue.push_back(source);
    visited[source] = true;
    while (!queue.empty()) {
        const std::size_t current = queue.front();
        queue.pop_front();
        if (current == destination) {
            std::vector<std::size_t> path;
            int node = static_cast<int>(destination);
            while (node != -1) {
                path.insert(path.begin(), static_cast<std::size_t>(node));
                node = parent[node];
            }
            return path;
        }
        for (std::size_t next = 0; next < n; ++next)
            if (!visited[next] && adjacency(current, next) > zero) {
                visited[next] = true;
                parent[next] = static_cast<int>(current);
                queue.push_back(next);
            }
    }
    return std::vector<std::size_t>();
}

}  // namespace detail

/**
 * @param lambda (u) matrices of size (n x (h+1)): lambda[v](i,t) is the rate at
 *               which user v requests item i while it sits at node t
 * @param R      (u x n) routing matrices of size ((h+1) x (h+1))
 */
template <class T>
CacheGammaGraphResult<T> cache_gamma(const std::vector<Matrix<T>>& lambda,
                                     const std::vector<std::vector<Matrix<T>>>& R) {
    if (lambda.empty()) throw InputError("cache_gamma: no user streams");
    const std::size_t u = lambda.size();
    const std::size_t n = lambda[0].rows();
    if (lambda[0].cols() == 0) throw InputError("cache_gamma: empty lambda");
    const std::size_t h = lambda[0].cols() - 1;
    if (R.size() != u) throw InputError("cache_gamma: R and lambda disagree on the user count");
    for (std::size_t v = 0; v < u; ++v)
        if (R[v].size() != n)
            throw InputError("cache_gamma: R and lambda disagree on the item count");

    const T zero = num_traits<T>::from_int(0);
    CacheGammaGraphResult<T> res;
    res.u = u;
    res.n = n;
    res.h = h;
    res.gamma = Matrix<T>(n, h, zero);

    for (std::size_t i = 0; i < n; ++i) {
        const Matrix<T>& graph = R[0][i];
        for (std::size_t j = 0; j < h; ++j) {
            const std::vector<std::size_t> Pj = detail::cache_bfs_path(graph, 0, j);
            if (Pj.empty()) continue;  // unreachable list: no access factor
            T g = zero;
            for (std::size_t v = 0; v < u; ++v) g += lambda[v](i, 0);
            for (std::size_t li = 1; li < Pj.size(); ++li) {
                const std::size_t a = Pj[li - 1];
                const std::size_t b = Pj[li];
                T y = zero;
                for (std::size_t v = 0; v < u; ++v) y += lambda[v](i, a) * R[v][i](a, b);
                g *= y;
            }
            res.gamma(i, j) = g;
        }
    }
    return res;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_GAMMA_H
