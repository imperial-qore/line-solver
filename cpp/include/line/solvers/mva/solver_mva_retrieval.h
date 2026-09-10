/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_RETRIEVAL_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_RETRIEVAL_H

/**
 * Delayed-hit (retrieval-system) cache analyzer, a port of
 * matlab/src/solvers/MVA/solver_mva_retrieval_analyzer.m.
 *
 * A miss triggers a per-item retrieval class that fetches the item through the
 * retrieval queues and returns to the cache; a read arriving while a fetch is in
 * flight is a DELAYED hit. The fixed-point algorithms decide the hit / miss /
 * delayed-hit ratios (`retrieval_fpi`) and the expected latency and per-item
 * queueing (`retrieval_fpi_latency`); this glue reads the inputs
 * (`cache_retrieval_inputs`) and assembles the per-station mean queue length,
 * utilization, response time and throughput of the retrieval sub-network.
 *
 * ARITHMETIC: transcendental (the FPIs iterate to a tolerance and fit
 * distributions), so it refuses under Rational by name.
 */

#include <limits>
#include <vector>

#include "line/api/retrieval/cache_retrieval_inputs.h"
#include "line/api/retrieval/retrieval_fpi.h"
#include "line/api/retrieval/retrieval_fpi_latency.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/linalg.h"

namespace line {
namespace mva {

/**
 * The cache half of the reference's return list: `hitprob`, `missprob`,
 * `delayedprob`, `latency`, `hitproblist` and `itemprob` of
 * `solver_mva_retrieval_analyzer.m`. Optional, because the analyzer is also
 * called for the station tables alone (cpp/examples/basic/cacheModel.cpp).
 */
template <class T>
struct MvaRetrievalCacheOutputs {
    std::vector<T> hitprob, missprob, delayedprob, latency;
    Matrix<T> hitproblist;  ///< (K x h)
    Matrix<T> itemprob;     ///< (n x h+1), column 0 = miss
};

template <class T>
MvaSolution<T> solver_mva_retrieval_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt,
                                             MvaRetrievalCacheOutputs<T>* cache_out = nullptr) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mva_retrieval_analyzer: the delayed-hit retrieval fixed point needs "
            "transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = L.nstations, K = L.nclasses;
        const retrieval::RetrievalInputs<T> in = retrieval::cache_retrieval_inputs(L);
        const std::size_t n = in.lambda.size(), S = in.queue_nodes.size();
        const std::size_t jobin = in.read_class;

        const retrieval::RetrievalFpiResult<T> fp =
            retrieval::retrieval_fpi(in.m, in.lambda, in.eta, in.gamma);
        const retrieval::RetrievalFpiLatencyResult<T> lat =
            retrieval::retrieval_fpi_latency(in.m, in.lambda, in.gamma, in.station, in.R);

        MvaSolution<T> s;
        s.Q = Matrix<T>(M, K, zero);
        s.U = Matrix<T>(M, K, zero);
        s.R = Matrix<T>(M, K, zero);
        s.Tp = Matrix<T>(M, K, zero);
        s.X.assign(K, zero);
        s.method = "fpi";
        s.iter = 1;

        // per-item weights and class aggregates (single read class, IRM)
        T lamsum = zero;
        for (std::size_t i = 0; i < n; ++i) lamsum = T(lamsum + in.lambda[i]);
        std::vector<T> w(n, zero);
        for (std::size_t i = 0; i < n; ++i) w[i] = lamsum > zero ? T(in.lambda[i] / lamsum) : zero;
        T hitAgg = zero, missAgg = zero, delayedAgg = zero;
        for (std::size_t i = 0; i < n; ++i) {
            T pih = zero;  // sum over lists of the hit ratio
            for (std::size_t l = 0; l < fp.phit.rows(); ++l) pih = T(pih + fp.phit(l, i));
            missAgg = T(missAgg + w[i] * fp.pmiss[i]);
            hitAgg = T(hitAgg + w[i] * pih);
            delayedAgg = T(delayedAgg + w[i] * lat.phi[i]);
        }

        // The cache answer, in the reference's per-class layout: NaN marks a
        // class that does not read this cache, which is what getAvgCacheTable
        // tests to decide the class has no row at all.
        if (cache_out) {
            const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
            cache_out->hitprob.assign(K, nan);
            cache_out->missprob.assign(K, nan);
            cache_out->delayedprob.assign(K, nan);
            cache_out->latency.assign(K, nan);
            if (jobin >= 1 && jobin <= K) {
                cache_out->hitprob[jobin - 1] = hitAgg;
                cache_out->missprob[jobin - 1] = missAgg;
                cache_out->delayedprob[jobin - 1] = delayedAgg;
                cache_out->latency[jobin - 1] = lat.Z;
            }
            // per-list hit fractions for the read class, access-weighted over
            // items: `(phit * w).'` of the reference. Rows sum to hitAgg.
            const std::size_t h = fp.phit.rows();
            cache_out->hitproblist = Matrix<T>(K, h, nan);
            if (jobin >= 1 && jobin <= K)
                for (std::size_t l = 0; l < h; ++l) {
                    T acc = zero;
                    for (std::size_t i = 0; i < n; ++i) acc = T(acc + fp.phit(l, i) * w[i]);
                    cache_out->hitproblist(jobin - 1, l) = acc;
                }
            // per-item occupancy `[pi0(:), phit.']`: column 0 is the miss.
            cache_out->itemprob = Matrix<T>(n, h + 1, zero);
            for (std::size_t i = 0; i < n; ++i) {
                cache_out->itemprob(i, 0) = fp.pmiss[i];
                for (std::size_t l = 0; l < h; ++l) cache_out->itemprob(i, l + 1) = fp.phit(l, i);
            }
        }

        // source throughput and the hit/miss class throughputs
        const std::size_t src_st = L.sourceIdx;
        const T sourceRate =
            L.disabled[src_st - 1][jobin - 1] ? zero : L.rates(src_st - 1, jobin - 1);
        s.Tp(src_st - 1, jobin - 1) = sourceRate;
        const qn::CacheParam<T>& ch = [&]() -> const qn::CacheParam<T>& {
            for (const auto& kv : L.nodeparam)
                if (L.nodes[kv.first - 1].nodetype == qn::NodeType::Cache) return kv.second;
            throw InputError("solver_mva_retrieval_analyzer: no cache");
        }();
        const std::size_t hc = (jobin - 1) < ch.hitclass.size() ? ch.hitclass[jobin - 1] : 0;
        const std::size_t mc = (jobin - 1) < ch.missclass.size() ? ch.missclass[jobin - 1] : 0;
        if (hc > 0) s.X[hc - 1] = T(sourceRate * (hitAgg + delayedAgg));
        if (mc > 0) s.X[mc - 1] = T(sourceRate * missAgg);

        // per-retrieval-station occupancy (delayed hits), throughput, response
        std::vector<std::size_t> psIdx;
        for (std::size_t si = 0; si < S; ++si)
            if (in.station[si].type != retrieval::RetrievalStationType::IS) psIdx.push_back(si);
        for (std::size_t si = 0; si < S; ++si) {
            const std::size_t st = L.nodes[in.queue_nodes[si] - 1].station;
            std::size_t prow = 0;  // IS aggregate row
            if (in.station[si].type != retrieval::RetrievalStationType::IS) {
                for (std::size_t p = 0; p < psIdx.size(); ++p)
                    if (psIdx[p] == si) prow = 1 + p;
            }
            T phi_s = zero;
            if (prow < fp.pdh.rows())
                for (std::size_t i = 0; i < n; ++i) phi_s = T(phi_s + fp.pdh(prow, i));
            // throughput: miss traffic times the per-item station visits
            T tput_s = zero;
            for (std::size_t i = 0; i < n; ++i) {
                Matrix<T> ImP(S, S, zero);
                std::vector<T> a(S, zero);
                for (std::size_t x = 0; x < S; ++x) {
                    a[x] = in.R[i](0, x + 1);
                    for (std::size_t y = 0; y < S; ++y)
                        ImP(x, y) = T((x == y ? num_traits<T>::from_int(1) : zero) -
                                      in.R[i](x + 1, y + 1));
                }
                const Matrix<T> F = inverse(ImP);
                T vis = zero;
                for (std::size_t y = 0; y < S; ++y) vis = T(vis + a[y] * F(y, si));
                tput_s = T(tput_s + sourceRate * w[i] * fp.pmiss[i] * vis);
            }
            s.Q(st - 1, jobin - 1) = phi_s;
            s.U(st - 1, jobin - 1) = phi_s;
            s.Tp(st - 1, jobin - 1) = tput_s;
            if (tput_s > zero) s.R(st - 1, jobin - 1) = T(phi_s / tput_s);
        }
        (void)opt;
        return s;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_RETRIEVAL_H
