/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_GAMMA_LP_H
#define LINE_API_CACHE_GAMMA_LP_H

/**
 * Access factors of a tree-structured multi-list cache.
 *
 * Templated port of matlab/src/api/cache/cache_gamma_lp.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_gamma_lp.java.
 *
 * The lists of the cache form a tree rooted at node 0 ("not cached"), list l
 * being node l+1. The access factor gamma(i,l) of item i at list l is the
 * product, along the unique path from the root to node l+1, of the aggregate
 * request flow crossing each edge:
 *
 *   gamma(i,l) = prod_{edges (a,b) of the path} sum_v sum_{t<=a} lambda(v,i,t) R{v,i}(a,b)
 *
 * The path is recovered by walking up from node l+1 through the parent
 * relation, the parent of a node being the unique earlier node with a nonzero
 * routing probability into it; more than one parent means the structure is not
 * a tree and is an error, as in both reference implementations.
 *
 * Only sums and products, so this instantiates at exact arithmetic and the
 * access factors of a rational model are exact rationals.
 *
 * REFERENCE DEFECT (both codebases): the `isempty(Pij)` branch that sets
 * gamma(i,l) = 0 is unreachable, because Pij is seeded with the node itself
 * and so is never empty. A list disconnected from the root therefore does not
 * yield a zero access factor; it yields the empty product 1. Reproduced here
 * so that the three codebases agree, and flagged rather than fixed.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_gamma_lp, mirroring [gamma,u,n,h]. */
template <class T>
struct CacheGammaResult {
    Matrix<T> gamma;  ///< (n x h) access factors
    std::size_t u;    ///< number of user streams
    std::size_t n;    ///< number of items
    std::size_t h;    ///< number of lists
    /** Parent list of each list, 0-based, -1 for lists rooted in the miss list. */
    std::vector<int> parent;
};

namespace detail {

/** Unique node a < j with R(a,j) != 0, or -1 when j is a root. */
template <class T>
int cache_parent(const Matrix<T>& R, std::size_t j) {
    const T zero = num_traits<T>::from_int(0);
    int parent = -1;
    for (std::size_t i = 0; i < j; ++i) {
        if (R(i, j) != zero) {
            if (parent >= 0)
                throw InputError(
                    "cache_gamma_lp: a cache list has more than one parent, but the structure must "
                    "be a tree");
            parent = static_cast<int>(i);
        }
    }
    return parent;
}

}  // namespace detail

/**
 * @param lambda (u) matrices of size (n x (h+1)): lambda[v](i,t) is the rate at
 *               which user v requests item i while it sits at node t
 * @param R      (u x n) routing matrices of size ((h+1) x (h+1))
 */
template <class T>
CacheGammaResult<T> cache_gamma_lp(const std::vector<Matrix<T>>& lambda,
                                   const std::vector<std::vector<Matrix<T>>>& R) {
    if (lambda.empty()) throw InputError("cache_gamma_lp: no user streams");
    const std::size_t u = lambda.size();
    const std::size_t n = lambda[0].rows();
    if (lambda[0].cols() == 0) throw InputError("cache_gamma_lp: empty lambda");
    const std::size_t h = lambda[0].cols() - 1;
    if (R.size() != u) throw InputError("cache_gamma_lp: R and lambda disagree on the user count");
    for (std::size_t v = 0; v < u; ++v)
        if (R[v].size() != n)
            throw InputError("cache_gamma_lp: R and lambda disagree on the item count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> gamma(n, h, zero);

    for (std::size_t i = 0; i < n; ++i) {
        // Rvi: routing of item i aggregated over users, used only for the tree
        // structure (the reference sums the matrices before taking parents).
        Matrix<T> Rvi(R[0][i].rows(), R[0][i].cols(), zero);
        for (std::size_t v = 0; v < u; ++v)
            for (std::size_t a = 0; a < Rvi.rows(); ++a)
                for (std::size_t b = 0; b < Rvi.cols(); ++b) Rvi(a, b) += R[v][i](a, b);

        for (std::size_t l = 0; l < h; ++l) {
            std::vector<std::size_t> path;
            path.push_back(l + 1);
            int pr = detail::cache_parent(Rvi, l + 1);
            while (pr >= 0) {
                path.insert(path.begin(), static_cast<std::size_t>(pr));
                pr = detail::cache_parent(Rvi, static_cast<std::size_t>(pr));
            }

            T g = one;
            for (std::size_t li = 1; li < path.size(); ++li) {
                const std::size_t a = path[li - 1];
                const std::size_t b = path[li];
                T y = zero;
                for (std::size_t v = 0; v < u; ++v)
                    for (std::size_t t = 0; t <= a; ++t) y += lambda[v](i, t) * R[v][i](t, b);
                g *= y;
            }
            gamma(i, l) = g;
        }
    }

    // Tree structure of the lists, read off item 0's routing matrix aggregated
    // over users -- the same matrix the gamma loop walks.
    Matrix<T> Rtot(R[0][0].rows(), R[0][0].cols(), zero);
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t a = 0; a < Rtot.rows(); ++a)
            for (std::size_t b = 0; b < Rtot.cols(); ++b) Rtot(a, b) += R[v][0](a, b);
    std::vector<int> parent(h, -1);
    for (std::size_t l = 0; l < h; ++l) {
        const int pr = detail::cache_parent(Rtot, l + 1);
        parent[l] = (pr < 0) ? -1 : (pr - 1);  // list indices, -1 = miss list
    }

    CacheGammaResult<T> r;
    r.gamma = gamma;
    r.u = u;
    r.n = n;
    r.h = h;
    r.parent = parent;
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_GAMMA_LP_H
