/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CACHE_METRICS_H
#define LINE_SOLVERS_CACHE_METRICS_H

/**
 * What a solver observed about the Cache nodes of a model.
 *
 * THE AvgTable CANNOT CARRY THIS AND IS NOT MEANT TO. A cache's answer is a hit
 * probability per (node, read class), a per-list breakdown of it, and a per-item
 * occupancy -- three different index spaces, none of them (station, class). The
 * reference keeps them in `getAvgCacheTable` and `getAvgItemTable` for exactly
 * that reason, and this struct is what those two tables are built from.
 *
 * EVERY FIELD IS OPTIONAL AND ABSENT MEANS NOT COMPUTED, never zero. A hit
 * probability of 0 is a cache that never hits; an empty vector is a solver that
 * did not measure one, and the tables print NaN there. The distinction matters
 * because the four cache branches of SolverNC compute different subsets: the
 * integrated caching-queueing network gives hit and miss but no per-item law,
 * the non-reentrant Source-Cache-Sink model gives the per-item law and the
 * per-list breakdown, and only the retrieval branches give a delayed-hit
 * fraction at all.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/util/matrix.h"

namespace line {
namespace solvers {

/** One Cache node's measured behaviour. */
template <class T>
struct CacheNodeMetrics {
    std::size_t node = 0;              ///< 1-based node index of the Cache
    /**
     * The Cache node's NAME, which is what a cross-language payload must key on.
     *
     * `node` is an index into THIS struct's node order, and that order is not
     * the model.json declaration order: on retrieval_simple the JSON declares
     * Source, Cache, Queue, Sink while the struct holds Source, Queue, Sink,
     * Cache, so the Cache is 2 to MATLAB and 4 here. A host restoring results
     * by index therefore wrote onto the Sink, found no Cache, and silently kept
     * the previous solver's numbers -- see CPPLINE.restoreCacheResults.
     */
    std::string name;
    std::vector<double> itemcap;       ///< (h) capacity of each list
    std::vector<double> itemsize;      ///< (n) storage cost per item, EMPTY without setItemSizes
    std::size_t nitems = 0;
    std::vector<T> hitprob;            ///< (K) TRUE hit fraction, EMPTY = not computed
    std::vector<T> missprob;           ///< (K)
    std::vector<T> delayedprob;        ///< (K) delayed-hit fraction, EMPTY off a retrieval system
    std::vector<T> latency;            ///< (K) expected retrieval latency, EMPTY = not computed
    Matrix<T> hitproblist;             ///< (K x h) per-list hit fraction, EMPTY = not computed
    Matrix<T> itemprob;                ///< (n x h+1), column 0 = miss; EMPTY = not computed
    std::vector<T> listcost;           ///< (h) mean storage cost held by each list
    /**
     * (n) mean secondary requests waiting on the in-flight fetch of each item,
     * and the same including the request that triggered the fetch. EMPTY off a
     * retrieval system and off every solver that does not form the per-item law
     * of block B; only the exact chain does, which is why the reference computes
     * them in SolverCTMC alone. They are the DelayedHitQLen columns of
     * `getAvgItemTable`.
     */
    std::vector<T> delayedhitqlen, delayedhitqlenfull;
};

/** Every Cache node of the model, in node order; empty on a model with none. */
template <class T>
struct CacheMetrics {
    std::vector<CacheNodeMetrics<T>> caches;
    bool empty() const { return caches.empty(); }
};

/**
 * Assemble `CacheMetrics` from what a cache analyzer returned.
 *
 * ONE PLACE, EVERY BRANCH OF EVERY SOLVER. Each cache analyzer computes a
 * different subset -- the retrieval ones give a delayed-hit fraction, the
 * non-reentrant one gives the per-item law, the integrated one gives neither --
 * and each leaves the rest empty. Assembling the struct here rather than in each
 * branch keeps "absent means not computed" a single rule, which is what lets
 * `getAvgCacheTable` print NaN in exactly the right places.
 *
 * IT LIVES HERE, BESIDE THE STRUCT IT BUILDS, rather than in a solver's runner:
 * SolverNC and SolverMVA both have cache branches and must assemble the answer
 * the same way, and a runner including another runner to borrow the helper is
 * how the two would drift apart.
 *
 * THE MODEL'S CACHES ARE READ FROM THE STRUCT, not from the analyzer: the caps,
 * the item count and the item sizes are model parameters and are reported even
 * where the solve measured nothing, so a caller can see the cache it described.
 */
template <class T>
CacheMetrics<T> cache_metrics_of(const qn::NetworkStruct<T>& sn, const std::vector<T>& hitprob,
                                 const std::vector<T>& missprob,
                                 const std::vector<T>& delayedprob, const std::vector<T>& latency,
                                 const Matrix<T>& hitproblist, const Matrix<T>& itemprob,
                                 const std::vector<T>& listcost) {
    CacheMetrics<T> out;
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        if (sn.nodes[ind - 1].nodetype != qn::NodeType::Cache) continue;
        const typename std::map<std::size_t, qn::CacheParam<T>>::const_iterator it =
            sn.nodeparam.find(ind);
        if (it == sn.nodeparam.end()) continue;
        CacheNodeMetrics<T> m;
        m.node = ind;
        m.name = sn.nodes[ind - 1].name;
        m.nitems = it->second.nitems;
        for (std::size_t l = 0; l < it->second.itemcap.size(); ++l)
            m.itemcap.push_back(static_cast<double>(it->second.itemcap[l]));
        for (std::size_t i = 0; i < it->second.itemsize.size(); ++i)
            m.itemsize.push_back(static_cast<double>(it->second.itemsize[i]));
        m.hitprob = hitprob;
        m.missprob = missprob;
        m.delayedprob = delayedprob;
        m.latency = latency;
        m.hitproblist = hitproblist;
        m.itemprob = itemprob;
        m.listcost = listcost;
        // The scalar hit and miss fractions of the non-reentrant branch, derived
        // from the per-list breakdown it does report. NOT invented: a read
        // either finds the item in some list or misses, so the row sum IS the
        // hit probability and 1 minus it IS the miss probability.
        if (m.hitprob.empty() && hitproblist.rows() > 0) {
            const T one = num_traits<T>::from_int(1);
            for (std::size_t r = 0; r < hitproblist.rows(); ++r) {
                T acc = num_traits<T>::from_int(0);
                bool any = false;
                for (std::size_t l = 0; l < hitproblist.cols(); ++l) {
                    const double v = num_traits<T>::to_double(hitproblist(r, l));
                    if (std::isnan(v)) continue;
                    acc = T(acc + hitproblist(r, l));
                    any = true;
                }
                m.hitprob.push_back(any ? acc : num_traits<T>::from_double(
                                                    std::numeric_limits<double>::quiet_NaN()));
                m.missprob.push_back(any ? T(one - acc)
                                         : num_traits<T>::from_double(
                                               std::numeric_limits<double>::quiet_NaN()));
            }
        }
        out.caches.push_back(m);
    }
    return out;
}

/**
 * The same, for the integrated caching-queueing branch, whose hit and miss
 * probabilities are (ncaches x nclasses) rather than one vector per model.
 */
template <class T>
CacheMetrics<T> cache_metrics_of_matrix(const qn::NetworkStruct<T>& sn, const Matrix<T>& hitprob,
                                        const Matrix<T>& missprob) {
    CacheMetrics<T> out =
        cache_metrics_of(sn, std::vector<T>(), std::vector<T>(), std::vector<T>(), std::vector<T>(),
                         Matrix<T>(), Matrix<T>(), std::vector<T>());
    for (std::size_t c = 0; c < out.caches.size(); ++c) {
        if (c >= hitprob.rows()) break;
        out.caches[c].hitprob.clear();
        out.caches[c].missprob.clear();
        for (std::size_t r = 0; r < hitprob.cols(); ++r) {
            out.caches[c].hitprob.push_back(hitprob(c, r));
            out.caches[c].missprob.push_back(missprob(c, r));
        }
    }
    return out;
}

}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_CACHE_METRICS_H
