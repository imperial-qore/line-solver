/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_H

/**
 * Integrated caching-queueing analyzer, a port of
 * matlab/src/solvers/MVA/solver_mva_cacheqn_analyzer.m.
 *
 * Delegates to the decomposition-aggregation driver `da_cacheqn`, supplying the
 * isolated-cache miss algorithm (exact `cache_mva` or the FPI approximation) and
 * a network solver. The surrounding queueing network is solved by
 * `solver_mva_analyzer` (or the load-dependent analyzer when the model carries
 * scaling), matching the reference's `netsolve`.
 *
 * ARITHMETIC: transcendental. The fixed-point driver stops on a tolerance and
 * the FPI miss path evaluates the cache access factors, so this refuses under
 * Rational by name.
 */

#include <limits>
#include <vector>

#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_ttl_lrua.h"
#include "line/api/da/da_cacheqn.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/solver_mva.h"

namespace line {
namespace mva {

/**
 * The cache half of the reference's return list: the (ncaches x nclasses) hit
 * and miss split, plus the per-item occupancy the reference writes straight onto
 * each Cache node with `setResultItemProb`.
 */
template <class T>
struct MvaCacheqnCacheOutputs {
    Matrix<T> hitprob;                    ///< (ncaches x nclasses)
    Matrix<T> missprob;                   ///< (ncaches x nclasses)
    std::vector<Matrix<T> > itemprob;     ///< per cache, (n x h+1); EMPTY = not computed
};

template <class T>
MvaSolution<T> solver_mva_cacheqn_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt,
                                           qn::NetworkStruct<T>* refreshed_out = nullptr,
                                           MvaCacheqnCacheOutputs<T>* cache_out = nullptr) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_cacheqn_analyzer: the integrated caching-queueing decomposition needs "
            "transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
        const MvaOptions netopt = opt;
        std::function<MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [netopt](const qn::NetworkStruct<T>& snit) -> MvaSolution<T> {
            bool has_scaling = false;
            for (const auto& st : snit.stations)
                if (!st.lldscaling.empty() || st.cdscaling) has_scaling = true;
            MvaOptions o = netopt;
            if (has_scaling) return solver_mvald_analyzer(snit, o, Matrix<T>()).sol;
            return solver_mva_analyzer(snit, o, Matrix<T>());
        };
        const bool exact = (opt.method == "exact");
        const da::CacheqnResult<T> r = da::da_cacheqn<T>(L, exact, opt, netfun);
        // Hand the runner a struct whose cache self-switch carries the ACTUAL
        // (converged) hit/miss split rather than the offered 1/2-1/2, so ArvR and
        // ResidT match MATLAB's setResultHitProb -> getStruct round trip. Unlike
        // da_cacheqn's over-routed inner struct, the self-switch is normalized
        // (no pass-through inflation), so it is correct whether the cache fans out
        // to one downstream node (cache_replc_routing) or several (gallery).
        if (refreshed_out) {
            *refreshed_out = L;
            refreshed_out->refresh_cacheqn_actual_visits(r.hitprob, r.missprob);
        }
        if (cache_out) {
            cache_out->hitprob = r.hitprob;
            cache_out->missprob = r.missprob;
            // Per-item occupancy from the CONVERGED access factors; SolverNC reads
            // the same law off its own fixed point.
            cache_out->itemprob = da::da_cacheqn_itemprob(r.info);
        }
        return r.res;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_H
