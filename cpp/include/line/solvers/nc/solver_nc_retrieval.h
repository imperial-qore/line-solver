/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_RETRIEVAL_H
#define LINE_SOLVERS_NC_SOLVER_NC_RETRIEVAL_H

/**
 * Port of `solver_nc_retrieval_analyzer.m`: the OPEN delayed-hit
 * (retrieval-system) cache.
 *
 * WHAT A DELAYED HIT IS, since the three-way split is the whole point. A miss
 * does not simply fail: it starts a FETCH that circulates the item through the
 * retrieval stations and back into the cache. A read arriving while that fetch
 * is still in flight is neither a hit (the item is not resident) nor an
 * ordinary miss (no second fetch is started) -- it is a DELAYED HIT, and it
 * waits for the fetch already running. So
 *
 *     true hit + delayed hit + miss = 1
 *
 * and only the MISS fraction starts new work. That is why the throughput split
 * below sends hits AND delayed hits to the hit class, and misses alone to the
 * miss class: the delayed hit is served by a fetch the miss already paid for.
 *
 * EXACT, NOT A FIXED POINT. SolverMVA solves this shape with `retrieval_fpi`,
 * iterating to a tolerance. SolverNC uses the product-form recurrences instead:
 * `retrieval_nc` for the normalizing constant and `retrieval_metrics` for the
 * three ratios, both exact. The reference reports `method = 'exact'` on that
 * basis and leaves LATENCY to SolverMVA, returning NaN -- the latency needs
 * `retrieval_fpi_latency`, which is a different algorithm and not this one's.
 *
 * THE RETRIEVAL STATIONS ARE READ OFF phi, NOT SOLVED. `pdh(s,i)` is the mean
 * number of copies of item i being fetched at station s, so summing over items
 * gives the station occupancy directly; the throughput comes from the per-item
 * fetch rate times the visit ratios of the fetch routing, and the response time
 * from Little's law. No queueing solve is involved.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/retrieval/cache_retrieval_inputs.h"
#include "line/api/retrieval/retrieval_metrics.h"
#include "line/api/retrieval/retrieval_nc.h"
#include "line/api/retrieval/retrieval_rayint.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** What the delayed-hit analyzer returns beyond the metric table. */
template <class T>
struct NcRetrievalSolution {
    NcSolution<T> sol;
    std::vector<T> hitprob;      ///< (K) TRUE hit fraction, NaN off the read class
    std::vector<T> missprob;     ///< (K)
    std::vector<T> delayedprob;  ///< (K) delayed-hit fraction
    std::vector<T> latency;      ///< (K) NaN: the latency belongs to SolverMVA
    Matrix<T> hitproblist;       ///< (K x h) per-list hit fraction
    Matrix<T> itemprob;          ///< (n x h+1) column 0 = miss, 1.. = per list
};

/**
 * Port of `solver_nc_retrieval_analyzer.m`.
 *
 * @param sn  the refreshed struct; the Cache must carry a retrieval system
 * @param opt solver controls
 */
template <class T>
NcRetrievalSolution<T> solver_nc_retrieval_analyzer(const qn::NetworkStruct<T>& sn,
                                                    const NcSolverOptions& opt) {
    NcRetrievalSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_retrieval_analyzer: the delayed-hit recurrences form a normalizing "
            "constant in logarithms and need transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const T nanT = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
        const std::size_t K = sn.nclasses, M = sn.nstations;

        const retrieval::RetrievalInputs<T> in = retrieval::cache_retrieval_inputs(sn);
        const std::size_t n = in.lambda.size();
        const std::size_t h = in.m.size();
        const std::size_t S = in.queue_nodes.size();
        std::size_t r = 0;  // the PS-like stations, which carry their own eta column
        for (const retrieval::RetrievalStationPH<T>& st : in.station)
            if (st.type != retrieval::RetrievalStationType::IS) ++r;

        // The ray expansion needs the delayed-hit constant to factorize. It does,
        // EXACTLY, when every fetch station is infinite-server: dividing the
        // retrieval_nc recurrence by prod_k D_k with D_k = 1 + lambda_k eta_{0,k}
        // collapses it onto cache_erec with theta_{k,j} = gamma_{k,j}/D_k, so the
        // delayed-hit cache IS a plain cache with fetch-inflated access factors.
        // A queueing (PS) fetch station breaks this: the (v_s+1) multiplicity ties
        // E(0,m) to the whole moment tower E(1_s,m), E(2_s,m), ..., and replacing it
        // by the retrieval_fpi mean field overestimates E by 13%/140%/830% at
        // n=6/8/10 (measured), growing with n. Refuse rather than return a
        // confident wrong number.
        bool useray = (opt.method == "rayint" || opt.method == "ray");
        if (useray) {
            std::string reason;
            bool has_ps = false;
            for (std::size_t i = 0; i < n && !has_ps; ++i)
                for (std::size_t sc = 1; sc <= r; ++sc)
                    if (num_traits<T>::to_double(in.eta(i, sc)) != 0.0) { has_ps = true; break; }
            long msum_chk = 0;
            for (std::size_t j = 0; j < h; ++j) msum_chk += in.m[j];
            if (has_ps)
                reason = "the retrieval system has a queueing (non infinite-server) fetch station";
            else if (msum_chk >= static_cast<long>(n))
                reason = "the cache is full (sum(m) >= n), where the saddle point escapes to infinity";
            if (!reason.empty()) {
                const std::string w = "SolverNC: method 'rayint' does not apply because " + reason +
                                      "; falling back to the exact recurrences.";
                out.sol.warning = out.sol.warning.empty() ? w : out.sol.warning + " " + w;
                useray = false;
            }
        }

        retrieval::RetrievalMetricsResult<T> mt;
        if (useray) {
            // --- ray (WKB) approximation, infinite-server fetch ---
            std::vector<T> D(n);
            Matrix<T> theta(n, h);
            for (std::size_t i = 0; i < n; ++i) {
                D[i] = one + in.lambda[i] * in.eta(i, 0);
                for (std::size_t j = 0; j < h; ++j) theta(i, j) = in.gamma(i, j) / D[i];
            }
            const retrieval::RetrievalRayintResult<T> ray =
                retrieval::retrieval_rayint(theta, in.m);
            double logD = 0.0;
            for (std::size_t i = 0; i < n; ++i) logD += std::log(num_traits<T>::to_double(D[i]));
            out.sol.sol.lG = logD + num_traits<T>::to_double(ray.log_e);

            // Same saddle as the constant, so the ratios are consistent with lG:
            // pi_{i,j} = theta_{i,j} xi_j / (1 + sum_l theta_{i,l} xi_l), and the
            // out-of-cache mass 1 - sum_j pi_{i,j} splits between a true miss (weight 1)
            // and an outstanding fetch (weight lambda_i eta_{0,i}) in proportion 1:D_i-1.
            mt.pmiss.assign(n, zero);
            mt.phit = Matrix<T>(h, n, zero);
            mt.pdh = Matrix<T>(1, n, zero);
            for (std::size_t i = 0; i < n; ++i) {
                T den = one;
                for (std::size_t j = 0; j < h; ++j) den += theta(i, j) * ray.xi[j];
                T pihit = zero;
                for (std::size_t j = 0; j < h; ++j) {
                    mt.phit(j, i) = theta(i, j) * ray.xi[j] / den;
                    pihit += mt.phit(j, i);
                }
                mt.pmiss[i] = (one - pihit) / D[i];
                mt.pdh(0, i) = in.lambda[i] * in.eta(i, 0) * mt.pmiss[i];
            }
        } else {
            // Exact normalizing constant E(m) = retrieval_nc(0, m, ...).
            const T E = retrieval::retrieval_nc(std::vector<int>(r, 0), in.m, in.lambda, in.eta,
                                                in.gamma);
            out.sol.sol.lG = std::log(num_traits<T>::to_double(E));

            mt = retrieval::retrieval_metrics(in.m, in.lambda, in.eta, in.gamma);
        }

        // Per-item aggregates: hit summed over lists, delayed summed over the
        // fetch stations.
        std::vector<T> pih(n, zero), phid(n, zero);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < mt.phit.rows(); ++j) pih[i] += mt.phit(j, i);
            for (std::size_t s = 0; s < mt.pdh.rows(); ++s) phid[i] += mt.pdh(s, i);
        }

        // The access-weighted item mixture.
        T lamtot = zero;
        for (const T& v : in.lambda) lamtot += v;
        std::vector<T> w(n, zero);
        for (std::size_t i = 0; i < n; ++i) w[i] = lamtot > zero ? T(in.lambda[i] / lamtot) : zero;

        T hitAgg = zero, missAgg = zero, delayedAgg = zero;
        for (std::size_t i = 0; i < n; ++i) {
            hitAgg += T(w[i] * pih[i]);
            missAgg += T(w[i] * mt.pmiss[i]);
            delayedAgg += T(w[i] * phid[i]);
        }

        const std::size_t rc = in.read_class;  // 1-based
        std::size_t sourceStation = 0;
        for (const qn::NodeDef& nd : sn.nodes)
            if (nd.nodetype == qn::NodeType::Source) sourceStation = nd.station;
        if (sourceStation == 0)
            throw UnsupportedError(
                "solver_nc_retrieval_analyzer: the OPEN delayed-hit analyzer needs a Source node; "
                "a closed integrated model is solved by solver_nc_cacheqn_retrieval_analyzer");
        std::vector<T> sourceRate(K, zero);
        for (std::size_t k = 0; k < K; ++k)
            if (!sn.disabled[sourceStation - 1][k]) sourceRate[k] = sn.rates(sourceStation - 1, k);

        out.sol.sol.Q = Matrix<T>(M, K, zero);
        out.sol.sol.U = Matrix<T>(M, K, zero);
        out.sol.sol.R = Matrix<T>(M, K, zero);
        out.sol.sol.Tp = Matrix<T>(M, K, zero);
        out.sol.sol.X.assign(K, zero);
        out.sol.sol.C.assign(K, zero);
        for (std::size_t k = 0; k < K; ++k) out.sol.sol.Tp(sourceStation - 1, k) = sourceRate[k];

        out.hitprob.assign(K, nanT);
        out.missprob.assign(K, nanT);
        out.delayedprob.assign(K, nanT);
        out.latency.assign(K, nanT);
        out.hitprob[rc - 1] = hitAgg;
        out.missprob[rc - 1] = missAgg;
        out.delayedprob[rc - 1] = delayedAgg;

        out.hitproblist = Matrix<T>(K, h, nanT);
        for (std::size_t j = 0; j < h && j < mt.phit.rows(); ++j) {
            T acc = zero;
            for (std::size_t i = 0; i < n; ++i) acc += T(mt.phit(j, i) * w[i]);
            out.hitproblist(rc - 1, j) = acc;
        }

        out.itemprob = Matrix<T>(n, h + 1, zero);
        for (std::size_t i = 0; i < n; ++i) {
            out.itemprob(i, 0) = mt.pmiss[i];
            for (std::size_t j = 0; j < h && j < mt.phit.rows(); ++j)
                out.itemprob(i, j + 1) = mt.phit(j, i);
        }

        // A DELAYED HIT IS SERVED BY A FETCH THE MISS ALREADY STARTED, so it
        // leaves through the HIT class; only the miss fraction starts new work.
        std::size_t cacheNode = 0;
        for (std::size_t i = 0; i < sn.nodes.size(); ++i)
            if (sn.nodes[i].nodetype == qn::NodeType::Cache) cacheNode = i + 1;
        if (cacheNode == 0) throw UnsupportedError("solver_nc_retrieval_analyzer: no Cache node");
        const qn::CacheParam<T>& ch = sn.nodeparam.at(cacheNode);
        if (rc - 1 < ch.hitclass.size() && ch.hitclass[rc - 1] > 0)
            out.sol.sol.X[ch.hitclass[rc - 1] - 1] =
                T(sourceRate[rc - 1] * (hitAgg + delayedAgg));
        if (rc - 1 < ch.missclass.size() && ch.missclass[rc - 1] > 0)
            out.sol.sol.X[ch.missclass[rc - 1] - 1] = T(sourceRate[rc - 1] * missAgg);

        // The retrieval stations, read off phi rather than solved. `pdh` row 0
        // is the aggregate IS row and rows 1..r are the PS-like stations in
        // order, which is why the row index is a running count and not `s`.
        std::vector<std::size_t> psRow(S, 0);
        {
            std::size_t seen = 0;
            for (std::size_t s = 0; s < S; ++s)
                if (in.station[s].type != retrieval::RetrievalStationType::IS)
                    psRow[s] = ++seen;  // 1-based row in pdh
                else
                    psRow[s] = 0;       // the aggregate IS row
        }
        for (std::size_t s = 0; s < S; ++s) {
            const std::size_t nd = in.queue_nodes[s];
            const std::size_t ist = sn.nodes[nd - 1].station;
            if (ist == 0) continue;
            T phi_s = zero;
            if (psRow[s] < mt.pdh.rows())
                for (std::size_t i = 0; i < n; ++i) phi_s += mt.pdh(psRow[s], i);
            // Fetch throughput: the per-item fetch rate times the visit ratios
            // of the fetch routing, v = a (I - P)^{-1}.
            T tput_s = zero;
            for (std::size_t i = 0; i < n; ++i) {
                const Matrix<T>& Ri = in.R[i];
                std::vector<T> a(S, zero);
                Matrix<T> Pm(S, S, zero);
                for (std::size_t x = 0; x < S; ++x) {
                    a[x] = Ri(0, x + 1);
                    for (std::size_t y = 0; y < S; ++y) Pm(x, y) = Ri(x + 1, y + 1);
                }
                Matrix<T> A(S, S, zero);
                for (std::size_t x = 0; x < S; ++x) {
                    for (std::size_t y = 0; y < S; ++y) A(y, x) = T(-Pm(x, y));
                    A(x, x) = T(A(x, x) + one);
                }
                const std::vector<T> vis = solve(A, a);
                tput_s += T(sourceRate[rc - 1] * w[i] * mt.pmiss[i] * vis[s]);
            }
            out.sol.sol.Q(ist - 1, rc - 1) = phi_s;
            out.sol.sol.U(ist - 1, rc - 1) = phi_s;
            out.sol.sol.Tp(ist - 1, rc - 1) = tput_s;
            if (tput_s > zero) out.sol.sol.R(ist - 1, rc - 1) = T(phi_s / tput_s);
        }

        out.sol.sol.iter = 1;
        out.sol.sol.method = useray ? "rayint" : "exact";
        out.sol.actualmethod = useray ? "rayint" : "exact";
        return out;
    }
}

/** True when the model's Cache carries a delayed-hit retrieval system. */
template <class T>
bool nc_has_retrieval(const qn::NetworkStruct<T>& sn) {
    for (const auto& kv : sn.nodeparam)
        if (kv.second.retrieval_capacity > 0) return true;
    return false;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_RETRIEVAL_H
