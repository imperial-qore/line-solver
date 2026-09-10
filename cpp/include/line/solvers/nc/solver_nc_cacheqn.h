/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_H
#define LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_H

/**
 * Port of `solver_nc_cacheqn_analyzer.m`: the INTEGRATED caching-queueing
 * network, where a Cache sits inside a queueing network rather than between a
 * Source and a Sink.
 *
 * WHY THIS NEEDS A FIXED POINT AND THE NON-REENTRANT CACHE DOES NOT. There, the
 * read rate is the Source's and is known up front. Here the flow reaching the
 * cache depends on the queueing network's throughputs, which depend in turn on
 * how the cache splits that flow between its hit and miss classes. So the two
 * are alternated: solve the caches in ISOLATION at the current arrival rates,
 * relabel the routing with the resulting hit/miss split, solve the queueing
 * network, and repeat. That alternation is `da_cacheqn`, shared with SolverMVA.
 *
 * WHAT THIS FILE SUPPLIES, and it is only two things: the isolated-cache MISS
 * ALGORITHM and the NETWORK SOLVER. Both differ from SolverMVA's:
 *
 *   miss     exact -> `cache_prob_erec`, otherwise `cache_miss_spm`
 *            (SolverMVA uses `cache_mva` / `cache_miss_fpi`)
 *   network  `solver_nc_analyzer`, or `solver_ncld_analyzer` under scaling
 *            (SolverMVA uses its own analyzers)
 *
 * THE DRIVER TAKES A HANDLE AGAIN. `da_cacheqn.m` has always taken `missfun` as
 * an argument; the C++ driver had collapsed it to a boolean, which hardcoded
 * SolverMVA's pair into a function both solvers share. The handle is restored,
 * defaulting to the MVA pair so that solver's call site is unchanged.
 */

#include <functional>
#include <vector>

#include "line/api/cache/cache_miss_spm.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/da/da_cacheqn.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_dispatch.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

/** What the integrated analyzer returns: the metrics plus the converged split. */
template <class T>
struct NcCacheqnSolution {
    NcSolution<T> sol;
    Matrix<T> hitprob;   ///< (ncaches x nclasses) converged hit probability
    Matrix<T> missprob;  ///< (ncaches x nclasses)
    std::vector<Matrix<T> > itemprob;  ///< per cache, (n x h+1); EMPTY = not computed
    /**
     * The struct whose cache self-switch carries the CONVERGED split rather
     * than the offered one, from which the runner derives ArvR and ResidT.
     */
    qn::NetworkStruct<T> refreshed;
};

/**
 * Port of `solver_nc_cacheqn_analyzer.m`.
 *
 * @param sn  the refreshed struct, carrying at least one Cache node
 * @param opt solver controls; `method == "exact"` selects the exact recursion
 */
template <class T>
NcCacheqnSolution<T> solver_nc_cacheqn_analyzer(const qn::NetworkStruct<T>& sn,
                                                const NcSolverOptions& opt) {
    NcCacheqnSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_cacheqn_analyzer: the integrated caching-queueing decomposition alternates "
            "two tolerance-stopped solves and needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);

        // The reference's admissibility gate, per cache. Note it compares
        // against the SUM of the list capacities here, where the non-reentrant
        // analyzer compares against the capacity vector entrywise.
        for (const auto& kv : sn.nodeparam) {
            const qn::CacheParam<T>& ch = kv.second;
            long capsum = 0;
            for (int c : ch.itemcap) capsum += c;
            if (static_cast<long>(ch.nitems) < capsum + 2)
                throw UnsupportedError(
                    "solver_nc_cacheqn_analyzer: NC requires the number of items to exceed the "
                    "cache capacity at least by 2; this cache holds " +
                    std::to_string(ch.nitems) + " items with total capacity " +
                    std::to_string(capsum));
        }

        const bool exact = (opt.method == "exact");
        // EXACT IS EXACT ONLY FOR THE EXCHANGEABLE FAMILY, and this analyzer used
        // to take `--method exact` on any policy: `cache_prob_erec` would then
        // return the RR/FIFO answer for an LRU cache SILENTLY, which is the one
        // outcome `solver_nc_cache_analyzer` refuses by name for the isolated
        // model (see its guard, same argument and near-identical wording). The
        // two analyzers differ in whether the cache is embedded in a queueing
        // network, not in which replacement policies the recursion is valid for.
        if (exact) {
            for (const auto& kv : sn.nodeparam) {
                const qn::CacheParam<T>& ch = kv.second;
                if (ch.replacestrat != lang::ReplacementStrategy::RR &&
                    ch.replacestrat != lang::ReplacementStrategy::FIFO)
                    throw UnsupportedError(
                        "solver_nc_cacheqn_analyzer: NC does not support the exact solution of "
                        "this cache replacement policy -- only RR and FIFO are exchangeable, and "
                        "a recency-based policy (LRU, h-LRU, q-LRU, CLIMB) would silently receive "
                        "the exchangeable answer. Use the default (approximate) method or "
                        "SolverCTMC");
            }
        }
        const NcSolverOptions netopt = opt;

        // The network solver: the ordinary NC analyzer, or the load-dependent
        // one when the model carries scaling. Matches the reference's netsolve.
        std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [netopt](const qn::NetworkStruct<T>& snit) -> mva::MvaSolution<T> {
            return nc_dispatch(snit, netopt).sol;
        };

        // The isolated-cache miss algorithm, which is what distinguishes this
        // analyzer from SolverMVA's.
        std::function<std::vector<T>(const Matrix<T>&, const std::vector<int>&,
                                     const std::vector<Matrix<T> >&, const qn::CacheParam<T>&)>
            missfun = [exact](const Matrix<T>& gamma, const std::vector<int>& m,
                              const std::vector<Matrix<T> >& lambda_cache,
                              const qn::CacheParam<T>& ch) -> std::vector<T> {
            (void)ch;  // the erec/spm pair reads the cache through `gamma` alone
            const std::size_t u = lambda_cache.size();
            std::vector<T> missrate(u, num_traits<T>::from_int(0));
            if (u == 0) return missrate;
            const std::size_t n = lambda_cache[0].rows();
            if (exact) {
                const Matrix<T> pij = cache::cache_prob_erec(gamma, m);
                for (std::size_t v = 0; v < u; ++v) {
                    T acc = num_traits<T>::from_int(0);
                    for (std::size_t k = 0; k < n && k < pij.rows(); ++k)
                        acc = T(acc + lambda_cache[v](k, 0) * pij(k, 0));
                    missrate[v] = acc;
                }
            } else {
                Matrix<T> lam_un(u, n, num_traits<T>::from_int(0));
                for (std::size_t v = 0; v < u; ++v)
                    for (std::size_t k = 0; k < n; ++k) lam_un(v, k) = lambda_cache[v](k, 0);
                const cache::CacheMissSpmResult<T> mr = cache::cache_miss_spm(gamma, m, lam_un);
                for (std::size_t v = 0; v < mr.MU.size() && v < u; ++v) missrate[v] = mr.MU[v];
            }
            return missrate;
        };

        mva::MvaOptions mopt;
        mopt.method = opt.method;
        mopt.tol = opt.tol;
        mopt.iter_tol = opt.iter_tol;
        mopt.iter_max = opt.iter_max;

        const da::CacheqnResult<T> r = da::da_cacheqn<T>(sn, exact, mopt, netfun, missfun);
        out.sol.sol = r.res;
        out.sol.sol.iter = r.iter;
        out.sol.actualmethod = exact ? "exact" : "spm";
        out.sol.sol.method = out.sol.actualmethod;
        out.hitprob = r.hitprob;
        out.missprob = r.missprob;
        // Per-item occupancy from the CONVERGED access factors, as SolverMVA reports it.
        out.itemprob = da::da_cacheqn_itemprob(r.info);
        // The runner reads ArvR and ResidT off a struct whose cache self-switch
        // carries the CONVERGED split, not the offered one.
        out.refreshed = sn;
        out.refreshed.refresh_cacheqn_actual_visits(r.hitprob, r.missprob);
        (void)zero;
        (void)one;
        return out;
    }
}

/** True when the model has a Cache node and is not the Source-Cache-Sink shape. */
template <class T>
bool nc_is_cacheqn(const qn::NetworkStruct<T>& sn) {
    bool hasCache = false;
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == qn::NodeType::Cache) hasCache = true;
    return hasCache;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_H
