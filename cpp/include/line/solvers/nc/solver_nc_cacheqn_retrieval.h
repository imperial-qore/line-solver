/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_RETRIEVAL_H
#define LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_RETRIEVAL_H

/**
 * Port of `solver_nc_cacheqn_retrieval_analyzer.m`: a CLOSED integrated
 * cache-queueing model whose Cache carries a delayed-hit retrieval system.
 *
 * THIS FILE IS GLUE. All the work is in `da_cacheqn_retrieval`, which alternates
 * the isolated-cache solve with the network solve; this supplies only the
 * NETWORK SOLVER and unpacks the result. The load-dependent path is not
 * optional here: the driver installs a coupon-collector `lldscaling` on the
 * fetch station, so the network solve must be the LOAD-DEPENDENT analyzer
 * whenever scaling is present -- which, after the driver has run, it always is.
 * `nc_dispatch` selects it on exactly that test.
 *
 * THE THREE-WAY SPLIT COLLAPSES ON THIS PATH, deliberately. The reference
 * reports `hitprob` as P(item cached) and folds the delayed-hit fraction into
 * `missprob`, returning `delayedprob = 0`. Only the OPEN analyzer
 * (`solver_nc_retrieval.h`) separates true hits from delayed hits. That is why
 * `hitproblist` and `latency` are NaN here and `itemprob` is empty: they are
 * quantities this path does not compute, and reporting a number for them would
 * be inventing one.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/da/da_cacheqn_retrieval.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_dispatch.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

/** What the closed delayed-hit analyzer returns. */
template <class T>
struct NcCacheqnRetrievalSolution {
    NcSolution<T> sol;
    std::vector<T> hitprob;      ///< (K) P(item cached)
    std::vector<T> missprob;     ///< (K) the delayed fraction is folded in here
    std::vector<T> delayedprob;  ///< (K) zero on this path, by the reference's convention
    std::vector<T> latency;      ///< (K) NaN: not computed on this path
    Matrix<T> hitproblist;       ///< (K x h) NaN: not computed on this path
};

/**
 * Port of `solver_nc_cacheqn_retrieval_analyzer.m`.
 *
 * @param sn  the refreshed struct; one Cache with a retrieval system, closed
 * @param opt solver controls
 */
template <class T>
NcCacheqnRetrievalSolution<T> solver_nc_cacheqn_retrieval_analyzer(
    const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    NcCacheqnRetrievalSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_cacheqn_retrieval_analyzer: the closed delayed-hit decomposition alternates "
            "two tolerance-stopped solves and needs transcendental arithmetic");
    } else {
        const T nanT = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
        const std::size_t K = sn.nclasses;
        const NcSolverOptions netopt = opt;

        std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [netopt](const qn::NetworkStruct<T>& snit) -> mva::MvaSolution<T> {
            return nc_dispatch(snit, netopt).sol;
        };

        mva::MvaOptions mopt;
        mopt.method = opt.method;
        mopt.tol = opt.tol;
        mopt.iter_tol = opt.iter_tol;
        mopt.iter_max = opt.iter_max;

        const da::CacheqnRetrievalResult<T> r =
            da::da_cacheqn_retrieval<T>(sn, netfun, mopt);

        out.sol.sol = r.res;
        out.sol.sol.iter = static_cast<int>(r.iter);
        // The reference reports 'fpi' here: the isolated-cache miss is
        // cache_miss_fpi on this path in BOTH solvers, so the name records the
        // algorithm that decided the split rather than the network method.
        out.sol.sol.method = "fpi";
        out.sol.actualmethod = "fpi";
        out.hitprob = r.hitprob;
        out.missprob = r.missprob;
        out.delayedprob = r.delayedprob;
        out.latency.assign(K, nanT);
        std::size_t h = 0;
        for (const auto& kv : sn.nodeparam)
            if (kv.second.itemcap.size() > h) h = kv.second.itemcap.size();
        out.hitproblist = Matrix<T>(K, h, nanT);
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_CACHEQN_RETRIEVAL_H
