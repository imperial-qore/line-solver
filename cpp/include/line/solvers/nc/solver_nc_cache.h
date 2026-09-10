/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_CACHE_H
#define LINE_SOLVERS_NC_SOLVER_NC_CACHE_H

/**
 * Port of `solver_nc_cache_analyzer.m`: the NON-REENTRANT cache, a model that
 * is exactly a Source, a Cache and a Sink.
 *
 * WHAT THE ANALYZER ACTUALLY COMPUTES. There is no queueing here at all. Each
 * class reads item k with probability `pread(v,k)`, the cache holds `itemcap`
 * items per list, and the question is only which items are resident. The answer
 * is a per-item occupancy `pij` -- column 0 the miss probability, column 1+j the
 * probability that the item sits on list j -- from which the miss RATE follows
 * by weighting with the read rates. The throughput table is then just the
 * source rate split between each class's hit and miss classes.
 *
 * THREE ALGORITHMS, AND ONLY ONE OF THEM IS EXACT.
 *
 *   exact     `cache_prob_erec`, the exact recursion. REFUSED for any
 *             replacement policy outside the exchangeable (product-form)
 *             family: RR and FIFO have a product form, LRU / h-LRU / q-LRU /
 *             CLIMB do not, and the recursion would silently return the
 *             exchangeable answer for them.
 *   sampling  `cache_miss_is` / `cache_prob_is`, importance sampling.
 *   default   the SPM saddle point, which `spm` and `rayint` also name. With
 *             per-item storage costs it is `cache_spm_size`, the size-tilted
 *             expansion, reported as `spm.size`; without them it is
 *             `cache_miss_spm` / `cache_prob_spm`, reported as `spm`.
 *
 * WHY THE PER-LIST BREAKDOWN IS NaN OUTSIDE THE EXACT BRANCH, which is a
 * deliberate refusal to report a number rather than an omission. Only the exact
 * recursion produces a miss column and per-list columns from ONE consistent
 * solution, so that the per-list rows sum to the aggregate hit. The approximate
 * algorithms derive the miss and the per-list columns from different expansions
 * and their breakdown does not form a distribution for more than one list.
 *
 * THE PER-ITEM TABLE IS ALWAYS TAKEN FROM THE EXACT RECURSION, even in the
 * approximate branches: the cache is product-form, so the exact per-item table
 * is available regardless of how the aggregate miss rate was obtained. It is
 * skipped above 10 items, where the recursion stops being tractable.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/cache/cache_gamma_lp.h"
#include <sstream>

#include "line/api/cache/cache_cost.h"
#include "line/api/cache/cache_miss_is.h"
#include "line/api/cache/cache_miss_spm.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_prob_is.h"
#include "line/api/cache/cache_prob_spm.h"
#include "line/api/cache/cache_spm_size.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** What the cache analyzer returns beyond the usual metric table. */
template <class T>
struct NcCacheSolution {
    NcSolution<T> sol;
    Matrix<T> pij;          ///< (n x h+1) per-item occupancy, column 0 = miss
    Matrix<T> itemprob;     ///< (n x h+1) the same from the EXACT recursion, NaN above 10 items
    Matrix<T> hitproblist;  ///< (u x h) access-weighted per-list hit probability, NaN if not exact
    std::vector<T> missrate;  ///< (u) per-class miss rate
    /**
     * (u) per-class miss and hit PROBABILITIES, `missrate` divided by the read
     * class's arrival rate and its complement.
     *
     * Carried beside the rate because `getAvgCacheTable` reports a probability
     * and the two differ by a factor no consumer downstream can recover: the
     * arrival rate is the Source's and is gone by the time the table is built.
     * NaN for a class with no read stream, which is not a hit ratio of zero.
     */
    std::vector<T> missprob, hitprob;
    std::vector<T> listcost;  ///< (h) mean storage cost held by each list, EMPTY without item sizes
    /**
     * Storage cost caps that block a promotion path, so that the exact
     * recursion normalizes over MORE states than the cache can reach. The port
     * has no warning channel, so this is a FLAG, like `SjnResult::capped`:
     * a non-empty vector means cross-check the answer with SolverLDES.
     */
    std::vector<cache::CacheBlockedPair> costcap_blocked;
    /** True when the requested method had no cost-capped counterpart and was switched. */
    bool costcap_method_switched = false;
};

/**
 * Port of `solver_nc_cache_analyzer.m`.
 *
 * @param sn  the refreshed struct; must be a Source-Cache-Sink model
 * @param opt solver controls; `method` selects exact / sampling / spm
 */
template <class T>
NcCacheSolution<T> solver_nc_cache_analyzer(const qn::NetworkStruct<T>& sn,
                                            const NcSolverOptions& opt) {
    NcCacheSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_cache_analyzer: the cache occupancy is a normalizing constant formed in "
            "logarithms and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const double dnan = std::numeric_limits<double>::quiet_NaN();
        const std::size_t K = sn.nclasses;

        std::size_t cacheNode = 0, sourceStation = 0;
        for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
            if (sn.nodes[i].nodetype == qn::NodeType::Cache) cacheNode = i + 1;
            if (sn.nodes[i].nodetype == qn::NodeType::Source)
                sourceStation = sn.nodes[i].station;
        }
        if (cacheNode == 0) throw UnsupportedError("solver_nc_cache_analyzer: no Cache node");
        if (sourceStation == 0) throw UnsupportedError("solver_nc_cache_analyzer: no Source node");
        const auto itp = sn.nodeparam.find(cacheNode);
        if (itp == sn.nodeparam.end())
            throw UnsupportedError("solver_nc_cache_analyzer: the Cache node carries no parameters");
        const qn::CacheParam<T>& ch = itp->second;

        std::vector<T> sourceRate(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (!sn.disabled[sourceStation - 1][r]) sourceRate[r] = sn.rates(sourceStation - 1, r);

        const std::vector<int>& m = ch.itemcap;
        const std::size_t n = ch.nitems;
        // `n < m + 2` in the reference, where `m` is the CAPACITY VECTOR, not
        // the list count -- and MATLAB's `if` on a vector requires every entry,
        // so the gate fires only when the item count is short for EVERY list.
        // Comparing against the number of lists instead lets a cache with more
        // capacity than items through, where the recursion has no headroom.
        {
            bool shortForAll = !m.empty();
            for (int cap : m)
                if (static_cast<long>(n) >= static_cast<long>(cap) + 2) shortForAll = false;
            if (shortForAll)
                throw UnsupportedError(
                    "solver_nc_cache_analyzer: NC requires the number of items to exceed the "
                    "cache capacity at least by 2; this cache holds " + std::to_string(n) +
                    " items with capacity " + std::to_string(m[0]));
        }
        const std::size_t h = m.size();
        const std::size_t u = K;

        // lambda[v](k,l): the rate at which class v requests item k while it
        // sits at node l. The reference fills every list column with the same
        // rate -- the read rate does not depend on where the item currently is.
        std::vector<Matrix<T>> lambda(u, Matrix<T>(n, h + 1, zero));
        for (std::size_t v = 0; v < u; ++v)
            if (!ch.pread[v].empty())
                for (std::size_t k = 0; k < n && k < ch.pread[v].size(); ++k)
                    for (std::size_t l = 0; l <= h; ++l)
                        lambda[v](k, l) = T(sourceRate[v] * ch.pread[v][k]);

        // The access-cost matrices. Absent, the reference installs the DEFAULT
        // LINEAR routing: an item moves from list l to list l+1 on a hit, and
        // the last list is absorbing.
        std::vector<std::vector<Matrix<T>>> R = ch.accost;
        if (R.empty()) {
            R.assign(u, std::vector<Matrix<T>>(n, Matrix<T>(h + 1, h + 1, zero)));
            for (std::size_t v = 0; v < u; ++v)
                for (std::size_t k = 0; k < n; ++k) {
                    Matrix<T> Rm(h + 1, h + 1, zero);
                    for (std::size_t l = 0; l < h; ++l) Rm(l, l + 1) = one;
                    Rm(h, h) = one;
                    R[v][k] = Rm;
                }
        }

        const cache::CacheGammaResult<T> gr = cache::cache_gamma_lp(lambda, R);

        // lambda(:,:,1) as the (u x n) matrix the miss routines want.
        Matrix<T> lam1(u, n, zero);
        for (std::size_t v = 0; v < u; ++v)
            for (std::size_t k = 0; k < n; ++k) lam1(v, k) = lambda[v](k, 0);

        // Per-item storage costs and per-list cost caps (ton21cache Sec. IX).
        const std::vector<int>& sigma = ch.itemsize;
        const std::vector<int>& costcap = ch.costcap;
        if (!costcap.empty()) {
            if (sigma.empty())
                throw InputError(
                    "solver_nc_cache_analyzer: storage cost caps require per-item sizes");
            if (sigma.size() != n)
                throw InputError(
                    "solver_nc_cache_analyzer: the item size vector must have one entry per item");
            if (costcap.size() != h)
                throw InputError(
                    "solver_nc_cache_analyzer: the cost cap vector must have one entry per cache "
                    "list");
            out.costcap_blocked =
                cache::cache_cost_pathcheck(gr.gamma, sigma, costcap, gr.parent);
            if (!out.costcap_blocked.empty()) {
                // Carried on the solution's warning channel, not just as a flag:
                // the reference warns here, and a flag no caller reads says nothing
                std::ostringstream w;
                w << "storage cost caps block the promotion path of item "
                  << (out.costcap_blocked[0].item + 1) << " into list "
                  << (out.costcap_blocked[0].list + 1) << " at list "
                  << (out.costcap_blocked[0].blocking_list + 1) << " (and "
                  << (out.costcap_blocked.size() - 1)
                  << " further pairs). The exact recursion normalizes over all size-feasible "
                     "states, which is then a strict superset of the states the cache can reach; "
                     "cross-check with SolverLDES.";
                out.sol.warning = w.str();
            }
        }

        std::string cacheMethod = opt.method;
        // 'rayint' is an alias of 'spm': on a cache both name the SPM saddle point,
        // and the method name stays live for solver_nc_retrieval_analyzer's delayed-hit
        // expansion.
        if (cacheMethod == "rayint" || cacheMethod == "spm") cacheMethod = "default";
        // The SPM family serves its size-tilted form (cache_spm_size) once the items
        // carry storage costs. The saddle escapes to infinity at sum(m) = n, so the
        // size-free saddle point takes over there rather than the exact recursion,
        // which would refuse every replacement policy outside RR/FIFO.
        long msum = 0;
        for (std::size_t j = 0; j < h; ++j) msum += m[j];
        const bool useSpmSize =
            cacheMethod == "default" && !sigma.empty() && msum < static_cast<long>(n);
        if (!costcap.empty() && !useSpmSize && cacheMethod != "exact" &&
            cacheMethod != "sampling") {
            // The size-free SPM and the mean-field methods have no cost-capped
            // counterpart: the k - sigma_i e_j argument couples item sizes into the
            // recursion graph, which only cache_spm_size and the exact path carry.
            double lattice = static_cast<double>(n);
            for (std::size_t j = 0; j < h; ++j)
                lattice *= static_cast<double>(m[j] + 1) * static_cast<double>(costcap[j] + 1);
            cacheMethod = (lattice <= 1e6) ? "exact" : "sampling";
            out.costcap_method_switched = true;
            std::ostringstream w;
            w << "method '" << opt.method << "' does not support storage cost caps; using '"
              << cacheMethod << "' instead.";
            out.sol.warning = out.sol.warning.empty() ? w.str() : out.sol.warning + " " + w.str();
        }

        std::vector<T> missRate(u, zero);
        bool exact = false;
        std::string method;
        if (cacheMethod == "exact") {
            // cache_prob_erec is exact only for the exchangeable family.
            if (ch.replacestrat != lang::ReplacementStrategy::RR &&
                ch.replacestrat != lang::ReplacementStrategy::FIFO)
                throw UnsupportedError(
                    "solver_nc_cache_analyzer: NC does not support the exact solution of this "
                    "cache replacement policy -- only RR and FIFO are exchangeable, and a "
                    "recency-based policy (LRU, h-LRU, q-LRU, CLIMB) would silently receive the "
                    "exchangeable answer. Use the default (approximate) method or SolverCTMC");
            out.pij = cache::cache_prob_erec(gr.gamma, m, sigma, costcap);
            for (std::size_t v = 0; v < u; ++v) {
                T acc = zero;
                for (std::size_t k = 0; k < n; ++k) acc += T(lambda[v](k, 0) * out.pij(k, 0));
                missRate[v] = acc;
            }
            exact = true;
            method = "exact";
        } else if (useSpmSize) {
            // Size-tilted SPM: a 2h Newton solve whose cost does not grow with the
            // (m,k) lattice the exact recursion walks. O(1/n), so it wants room
            // between the occupancies and n.
            std::vector<int> raycap = costcap;
            if (raycap.empty()) {
                // Sizes but no caps: cap each list at the dearest load it can hold,
                // which is exactly slack, so the cost coordinate leaves the saddle
                // and the expansion degenerates to the size-free one.
                std::vector<int> srt = sigma;
                std::sort(srt.begin(), srt.end(), std::greater<int>());
                raycap.assign(h, 0);
                for (std::size_t j = 0; j < h; ++j) {
                    int top = 0;
                    for (int a = 0; a < m[j]; ++a) top += srt[static_cast<std::size_t>(a)];
                    raycap[j] = top;
                }
            }
            const cache::CacheSpmSizeResult<T> ray =
                cache::cache_spm_size<T>(gr.gamma, m, sigma, raycap);
            out.pij = ray.pij;
            for (std::size_t v = 0; v < u; ++v) {
                T acc = zero;
                for (std::size_t k = 0; k < n; ++k) acc += T(lambda[v](k, 0) * out.pij(k, 0));
                missRate[v] = acc;
            }
            method = "spm.size";
        } else if (cacheMethod == "sampling") {
            const cache::CacheMissIsResult<T> mi =
                cache::cache_miss_is(gr.gamma, m, lam1, opt.samples, opt.seed, sigma, costcap);
            missRate = mi.MU;
            out.pij = cache::cache_prob_is(gr.gamma, m, opt.samples, opt.seed, sigma, costcap);
            method = "sampling";
        } else {
            // Size-free SPM, the default/spm/rayint branch with no item sizes.
            const cache::CacheMissSpmResult<T> ms = cache::cache_miss_spm(gr.gamma, m, lam1);
            missRate = ms.MU;
            out.pij = cache::cache_prob_spm(gr.gamma, m);
            method = "spm";
        }
        out.missrate = missRate;
        // The probabilities the cache table reports, formed where the arrival
        // rate is still in scope.
        out.missprob.assign(u, num_traits<T>::from_double(dnan));
        out.hitprob.assign(u, num_traits<T>::from_double(dnan));
        for (std::size_t v = 0; v < u; ++v) {
            if (v >= sourceRate.size()) break;
            const double lam = num_traits<T>::to_double(sourceRate[v]);
            if (!(lam > 0.0)) continue;
            out.missprob[v] = T(missRate[v] / sourceRate[v]);
            out.hitprob[v] = T(num_traits<T>::from_int(1) - out.missprob[v]);
        }

        // The metric table. There is no queueing, so only the throughput row of
        // the Source and the per-class hit/miss split carry anything.
        const std::size_t M = sn.nstations;
        out.sol.sol.Q = Matrix<T>(M, K, zero);
        out.sol.sol.U = Matrix<T>(M, K, zero);
        out.sol.sol.R = Matrix<T>(M, K, zero);
        out.sol.sol.Tp = Matrix<T>(M, K, zero);
        out.sol.sol.C.assign(K, zero);
        out.sol.sol.X.assign(K, zero);
        for (std::size_t r = 0; r < K; ++r) out.sol.sol.Tp(sourceStation - 1, r) = sourceRate[r];
        for (std::size_t r = 0; r < K; ++r) {
            if (r >= ch.hitclass.size() || r >= ch.missclass.size()) continue;
            const std::size_t hc = ch.hitclass[r], mc = ch.missclass[r];
            if (hc == 0 || mc == 0) continue;
            out.sol.sol.X[mc - 1] = T(out.sol.sol.X[mc - 1] + missRate[r]);
            out.sol.sol.X[hc - 1] = T(out.sol.sol.X[hc - 1] + (sourceRate[r] - missRate[r]));
        }
        out.sol.sol.lG = 0.0;
        out.sol.sol.iter = 1;
        out.sol.sol.method = method;
        out.sol.actualmethod = method;

        // Per-list hit probability, access-weighted. Reported ONLY from the
        // exact branch, where the per-list columns and the miss column come
        // from one solution; see the header.
        out.hitproblist = Matrix<T>(u, h, num_traits<T>::from_double(dnan));
        if (exact)
            for (std::size_t v = 0; v < u; ++v) {
                if (ch.pread[v].empty()) continue;
                for (std::size_t l = 0; l < h; ++l) {
                    T acc = zero;
                    for (std::size_t k = 0; k < n && k < ch.pread[v].size(); ++k)
                        acc += T(ch.pread[v][k] * out.pij(k, l + 1));
                    out.hitproblist(v, l) = acc;
                }
            }

        // Per-item occupancy. Always the EXACT recursion, since the cache is
        // product-form; skipped above 10 items where it stops being tractable.
        if (n > 10) {
            out.itemprob = Matrix<T>(n, h + 1, num_traits<T>::from_double(dnan));
        } else if (exact) {
            out.itemprob = out.pij;
        } else {
            out.itemprob = cache::cache_prob_erec(gr.gamma, m, sigma, costcap);
        }

        // Mean storage cost held by each list, K_j = sum_i sigma_i pi_ij. Kept on
        // the inner NcSolution too, since the runner returns only that and the
        // cost would otherwise never leave this function.
        if (!sigma.empty() && out.pij.rows() == n && out.pij.cols() == h + 1) {
            out.listcost = cache::cache_cost(gr.gamma, m, sigma, costcap, out.pij);
            out.sol.listcost = out.listcost;
        }
        return out;
    }
}

/** True when the model is exactly a Source, a Cache and a Sink. */
template <class T>
bool nc_is_noreentrant_cache(const qn::NetworkStruct<T>& sn) {
    if (sn.nodes.size() != 3) return false;
    int src = 0, ca = 0, snk = 0;
    for (const qn::NodeDef& nd : sn.nodes) {
        if (nd.nodetype == qn::NodeType::Source) ++src;
        else if (nd.nodetype == qn::NodeType::Cache) ++ca;
        else if (nd.nodetype == qn::NodeType::Sink) ++snk;
    }
    if (!(src == 1 && ca == 1 && snk == 1)) return false;
    for (const qn::JobClass& c : sn.classes)
        if (std::isfinite(c.population) && c.population > 0.0) return false;
    return true;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_CACHE_H
