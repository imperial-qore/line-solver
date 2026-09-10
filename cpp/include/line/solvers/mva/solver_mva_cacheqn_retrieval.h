/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_RETRIEVAL_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_RETRIEVAL_H

/**
 * Port of `solver_mva_cacheqn_retrieval_analyzer.m`: a CLOSED integrated
 * cache-queueing model whose Cache carries a delayed-hit retrieval system.
 *
 * THIS FILE IS GLUE, and the twin of `solver_nc_cacheqn_retrieval.h`. All the
 * work is in `da_cacheqn_retrieval`, which alternates the isolated-cache solve
 * with the network solve; this supplies only the NETWORK SOLVER and unpacks the
 * result. The two analyzers differ in that one line and in nothing else, which
 * is why the driver takes `netfun` as its only handle: the isolated-cache miss
 * algorithm is `cache_miss_fpi` in both, and both report `method = 'fpi'`.
 *
 * THE LOAD-DEPENDENT BRANCH IS NOT OPTIONAL. The driver installs a
 * coupon-collector `lldscaling` on the fetch station before the first sweep, so
 * `netsolve`'s scaling test is true on every call and the network solve is
 * always `solver_mvald_analyzer`. The `solver_mva_analyzer` arm is kept because
 * the reference keeps it -- it is what runs if a future driver stops installing
 * the scaling -- not because a model reaches it today.
 *
 * THE THREE-WAY SPLIT COLLAPSES ON THIS PATH, deliberately. The reference
 * reports `hitprob` as P(item cached) and folds the delayed-hit fraction into
 * `missprob`, returning `delayedprob = 0`. Only the OPEN analyzer
 * (`solver_mva_retrieval.h`) separates true hits from delayed hits. That is why
 * `hitproblist` and `latency` are NaN here and there is no `itemprob`: they are
 * quantities this path does not compute, and reporting a number for them would
 * be inventing one.
 *
 * EXPERIMENTAL, in the reference's own words (`da_cacheqn_retrieval.m`,
 * LIMITATIONS): the coalescing throughput benefit is captured in DIRECTION and
 * understated in magnitude, and no closed retrieval example ships in the suite.
 * The port reproduces the method rather than repairing it, so a closed model
 * should be validated against LDES before its absolute numbers are trusted.
 *
 * ARITHMETIC: transcendental. The fixed point stops on a tolerance and the
 * coupon-collector rate is a real power, so this refuses under Rational by name.
 */

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/da/da_cacheqn_retrieval.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mva {

/** What the closed delayed-hit analyzer returns. */
template <class T>
struct MvaCacheqnRetrievalSolution {
    MvaSolution<T> sol;
    std::vector<T> hitprob;      ///< (K) P(item cached), read class only
    std::vector<T> missprob;     ///< (K) the delayed fraction is folded in here
    std::vector<T> delayedprob;  ///< (K) zero on this path, by the reference's convention
    std::vector<T> latency;      ///< (K) NaN: not computed on this path
    Matrix<T> hitproblist;       ///< (K x h) NaN: not computed on this path
};

/**
 * Port of `solver_mva_cacheqn_retrieval_analyzer.m`.
 *
 * @param L   the refreshed struct; one Cache with a retrieval system, closed
 * @param opt solver controls
 */
template <class T>
MvaCacheqnRetrievalSolution<T> solver_mva_cacheqn_retrieval_analyzer(const qn::NetworkStruct<T>& L,
                                                                    const MvaOptions& opt) {
    MvaCacheqnRetrievalSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        throw UnsupportedError(
            "solver_mva_cacheqn_retrieval_analyzer: the closed delayed-hit decomposition alternates "
            "two tolerance-stopped solves and needs transcendental arithmetic; rerun with --arith "
            "double or --arith real");
    } else {
        const T nanT = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
        const std::size_t K = L.nclasses;
        const MvaOptions netopt = opt;

        // `netsolve` of the reference: the load-dependent analyzer whenever the
        // struct carries scaling, which after the driver's first mutation it
        // always does. The call is resolved by ADL at instantiation, as in
        // solver_mva_cacheqn.h, so that this header does not have to include
        // mva_dispatch.h and close a cycle with it.
        std::function<MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [netopt](const qn::NetworkStruct<T>& snit) -> MvaSolution<T> {
            bool has_scaling = false;
            for (const auto& st : snit.stations)
                if (!st.lldscaling.empty() || st.cdscaling) has_scaling = true;
            MvaOptions o = netopt;
            if (has_scaling) return solver_mvald_analyzer(snit, o, Matrix<T>()).sol;
            return solver_mva_analyzer(snit, o, Matrix<T>());
        };

        const da::CacheqnRetrievalResult<T> r = da::da_cacheqn_retrieval<T>(L, netfun, netopt);

        out.sol = r.res;
        out.sol.iter = static_cast<int>(r.iter);
        // The reference reports 'fpi' here: the isolated-cache miss is
        // cache_miss_fpi on this path in BOTH solvers, so the name records the
        // algorithm that decided the split rather than the network method.
        out.sol.method = "fpi";
        out.hitprob = r.hitprob;
        out.missprob = r.missprob;
        out.delayedprob = r.delayedprob;
        out.latency.assign(K, nanT);
        std::size_t h = 0;
        for (const auto& kv : L.nodeparam)
            if (kv.second.itemcap.size() > h) h = kv.second.itemcap.size();
        out.hitproblist = Matrix<T>(K, h, nanT);
        return out;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_CACHEQN_RETRIEVAL_H
