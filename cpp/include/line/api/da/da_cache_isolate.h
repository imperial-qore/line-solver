/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DA_DA_CACHE_ISOLATE_H
#define LINE_API_DA_DA_CACHE_ISOLATE_H

/**
 * Isolated-cache input construction for the decomposition methods.
 *
 * Templated port of matlab/src/api/da/da_cache_isolate.m. The decomposition
 * driver replaces a cache node embedded in a queueing network by the same
 * cache in isolation, driven by the current per-class arrival rates. This
 * routine builds that isolated model: it spreads each class rate lambda(v)
 * over the items through the class read distribution pread{v}, attaches the
 * access-cost (routing) matrices Rcost, and returns the access factors gamma
 * that every cache algorithm of the family consumes.
 *
 * The MATLAB source reads the cache parameters off sn.nodeparam of the cache
 * node; the port takes only the four fields it actually uses, in CacheParam,
 * so nothing of the NetworkStruct layer is needed here.
 *
 * A note on lambda_cache: the reference fills every list position l = 1..h+1
 * of item k with the same lambda(v) pread{v}(k). That is deliberate - the
 * request rate for an item does not depend on which list currently holds it -
 * and is reproduced verbatim, because cache_gamma_lp reads lambda(v,i,t) as
 * "rate at which v requests item i while it sits at node t".
 *
 * Only products and sums (here and in cache_gamma_lp), so this is a finite
 * field computation: instantiated at exact arithmetic the access factors of a
 * rational cache model are exact rationals. No transcendental gate.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_gamma_lp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace da {

/**
 * The fields of sn.nodeparam{cache} that da_cache_isolate reads.
 *
 * pread[v] is the read distribution of class v over the n items; an EMPTY
 * vector encodes MATLAB's NaN placeholder, i.e. "class v does not read this
 * cache", and leaves that class's rates at zero.
 *
 * accost[v][k] is the ((h+1) x (h+1)) list-to-list routing matrix of class v
 * on item k. Leaving accost empty selects the reference default, the linear
 * cache in which an item moves from list l to list l+1 on a hit and stays in
 * the last list once it gets there.
 */
template <class T>
struct CacheParam {
    std::vector<int> itemcap;                    ///< (h) list capacities
    std::size_t nitems = 0;                      ///< n
    std::vector<std::vector<T>> pread;           ///< (u) x (n), empty row = NaN
    std::vector<std::vector<Matrix<T>>> accost;  ///< (u) x (n) of (h+1)x(h+1), or empty
};

/** Return value of da_cache_isolate, mirroring [gamma,lambda_cache,Rcost]. */
template <class T>
struct CacheIsolateResult {
    Matrix<T> gamma;                            ///< (n x h) access factors
    std::vector<Matrix<T>> lambda_cache;        ///< (u) matrices of size n x (h+1)
    std::vector<std::vector<Matrix<T>>> Rcost;  ///< (u x n) of (h+1)x(h+1)
};

/**
 * @param ch     cache node parameters
 * @param lambda (u) per-class arrival rates at the cache
 */
template <class T>
CacheIsolateResult<T> da_cache_isolate(const CacheParam<T>& ch, const std::vector<T>& lambda) {
    const std::size_t h = ch.itemcap.size();
    const std::size_t n = ch.nitems;
    const std::size_t u = lambda.size();
    if (h == 0) throw InputError("da_cache_isolate: the cache has no lists");
    if (n == 0) throw InputError("da_cache_isolate: the cache has no items");
    if (u == 0) throw InputError("da_cache_isolate: no arrival rates given");
    if (ch.pread.size() != u)
        throw InputError("da_cache_isolate: pread and lambda disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    CacheIsolateResult<T> r;
    r.lambda_cache.assign(u, Matrix<T>(n, h + 1, zero));
    for (std::size_t v = 0; v < u; ++v) {
        if (ch.pread[v].empty()) continue;  // MATLAB's isnan(pread{v}) branch
        if (ch.pread[v].size() != n)
            throw InputError("da_cache_isolate: pread has the wrong number of items");
        for (std::size_t k = 0; k < n; ++k) {
            const T rate = lambda[v] * ch.pread[v][k];
            for (std::size_t l = 0; l <= h; ++l) r.lambda_cache[v](k, l) = rate;
        }
    }

    if (!ch.accost.empty()) {
        if (ch.accost.size() != u)
            throw InputError("da_cache_isolate: accost and lambda disagree on the class count");
        for (std::size_t v = 0; v < u; ++v) {
            if (ch.accost[v].size() != n)
                throw InputError("da_cache_isolate: accost has the wrong number of items");
            for (std::size_t k = 0; k < n; ++k)
                if (ch.accost[v][k].rows() != h + 1 || ch.accost[v][k].cols() != h + 1)
                    throw InputError("da_cache_isolate: an accost matrix is not (h+1) x (h+1)");
        }
        r.Rcost = ch.accost;
    } else {
        // Default linear cache routing: items flow from list l to list l+1,
        // and the last list is absorbing.
        Matrix<T> Rmat(h + 1, h + 1, zero);
        for (std::size_t l = 0; l < h; ++l) Rmat(l, l + 1) = one;
        Rmat(h, h) = one;
        r.Rcost.assign(u, std::vector<Matrix<T>>(n, Rmat));
    }

    r.gamma = cache::cache_gamma_lp(r.lambda_cache, r.Rcost).gamma;
    return r;
}

}  // namespace da
}  // namespace line

#endif  // LINE_API_DA_DA_CACHE_ISOLATE_H
