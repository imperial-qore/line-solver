/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_CACHE_RETRIEVAL_INPUTS_H
#define LINE_API_RETRIEVAL_CACHE_RETRIEVAL_INPUTS_H

/**
 * Extract the delayed-hit retrieval-algorithm inputs from a NetworkStruct, a
 * port of matlab/src/api/retrieval/cache_retrieval_inputs.m.
 *
 * Given a Cache equipped with a retrieval system (set_retrieval_system), rebuilds
 * the inputs the retrieval_* algorithms read: the list capacities m, the per-item
 * arrival rates lambda (readRate * pread), the access factors gamma
 * (cache_gamma_lp), the fetching demands eta, the per-station phase-type service
 * (alpha, T -> RetrievalStationPH), and the per-item routing R over the cache and
 * the retrieval stations. Single read class (IRM); IS / PS / SIRO / FCFS / LCFSPR
 * stations only, with SIRO/FCFS requiring exponential, class-independent service.
 *
 * ARITHMETIC: transcendental (the downstream FPI does), so this refuses under
 * Rational at the analyzer's guard; here it is field arithmetic plus PH means.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/cache/cache_gamma_lp.h"
#include "line/api/retrieval/retrieval_fpi_latency.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** The [m, lambda, gamma, eta, alpha/T (station), R] the retrieval algorithms read. */
template <class T>
struct RetrievalInputs {
    std::vector<int> m;                            ///< (h) list capacities
    std::vector<T> lambda;                         ///< (n) per-item arrival rates
    Matrix<T> gamma;                               ///< (n x h) access factors
    Matrix<T> eta;                                 ///< (n x (r+1)) fetching demands
    std::vector<RetrievalStationPH<T> > station;   ///< (S) per-station PH service + type
    std::vector<Matrix<T> > R;                     ///< (n) routing over cache+stations, (S+1)x(S+1)
    std::size_t read_class = 0;                    ///< 1-based read (job-in) class
    std::vector<std::size_t> queue_nodes;          ///< 1-based retrieval station nodes
};

/**
 * @param sn             the model struct with a retrieval-system cache
 * @param lambda_override the closed-sublayer read rate; <0 selects the Source
 *                        throughput (open model)
 */
template <class T>
RetrievalInputs<T> cache_retrieval_inputs(const qn::NetworkStruct<T>& sn,
                                          double lambda_override = -1.0) {
    const T zero = num_traits<T>::from_int(0);
    // -- locate the cache and its retrieval system --------------------------
    std::size_t ci = 0;
    for (const auto& kv : sn.nodeparam)
        if (sn.nodes[kv.first - 1].nodetype == qn::NodeType::Cache) {
            if (ci != 0) throw InputError("cache_retrieval_inputs: more than one Cache node");
            ci = kv.first;
        }
    if (ci == 0) throw InputError("cache_retrieval_inputs: no Cache node");
    const qn::CacheParam<T>& ch = sn.nodeparam.at(ci);
    if (ch.retrieval_capacity <= 0)
        throw InputError("cache_retrieval_inputs: the Cache has no retrieval system");
    if (ch.retrieval_queues.size() != 1)
        throw UnsupportedError("cache_retrieval_inputs: a single read class is supported");

    RetrievalInputs<T> out;
    out.m = ch.itemcap;
    const std::size_t n = ch.nitems, h = ch.itemcap.size(), K = sn.nclasses;

    const std::size_t rd0 = ch.retrieval_queues.begin()->first;  // 0-based read class
    const std::size_t jobin = rd0 + 1;                           // 1-based
    out.read_class = jobin;
    out.queue_nodes = ch.retrieval_queues.begin()->second;       // 1-based nodes
    const std::vector<std::size_t>& qn = out.queue_nodes;
    const std::size_t S = qn.size();
    if (S == 0) throw InputError("cache_retrieval_inputs: the retrieval system has no stations");

    // -- per-item arrival rates lambda(i) = readRate * pread(i) -------------
    T readRate;
    if (lambda_override >= 0.0) {
        readRate = num_traits<T>::from_double(lambda_override);
    } else {
        const std::size_t src_st = sn.stations[sn.sourceIdx - 1].nodetype == qn::NodeType::Source
                                       ? sn.sourceIdx
                                       : 0;
        if (src_st == 0) throw InputError("cache_retrieval_inputs: open model needs a Source");
        readRate = sn.disabled[src_st - 1][jobin - 1] ? zero : sn.rates(src_st - 1, jobin - 1);
    }
    const std::vector<T>& pread = ch.pread[jobin - 1];
    out.lambda.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i) out.lambda[i] = T(readRate * pread[i]);

    // -- gamma via cache_gamma_lp (one user stream) -------------------------
    Matrix<T> lam3d(n, h + 1, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t l = 0; l <= h; ++l) lam3d(i, l) = out.lambda[i];
    std::vector<std::vector<Matrix<T> > > Rcost(1, std::vector<Matrix<T> >());
    if (!ch.accost.empty() && !ch.accost[jobin - 1].empty()) {
        Rcost[0] = ch.accost[jobin - 1];
    } else {
        // default linear cache: item flows from list l to l+1, absorbs at h.
        Rcost[0].resize(n, Matrix<T>(h + 1, h + 1, zero));
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t l = 0; l < h; ++l) Rcost[0][k](l, l + 1) = num_traits<T>::from_int(1);
            Rcost[0][k](h, h) = num_traits<T>::from_int(1);
        }
    }
    out.gamma = cache::cache_gamma_lp(std::vector<Matrix<T> >(1, lam3d), Rcost).gamma;

    // -- station types ------------------------------------------------------
    std::vector<RetrievalStationType> stype(S);
    for (std::size_t s = 0; s < S; ++s) {
        const std::size_t st = sn.nodes[qn[s] - 1].station;
        switch (sn.stations[st - 1].sched) {
            case qn::SchedStrategy::INF: stype[s] = RetrievalStationType::IS; break;
            case qn::SchedStrategy::PS: stype[s] = RetrievalStationType::PS; break;
            case qn::SchedStrategy::SIRO: stype[s] = RetrievalStationType::SIRO; break;
            case qn::SchedStrategy::FCFS: stype[s] = RetrievalStationType::FCFS; break;
            case qn::SchedStrategy::LCFSPR: stype[s] = RetrievalStationType::LCFSPR; break;
            default:
                throw UnsupportedError(
                    "cache_retrieval_inputs: a retrieval station uses an unsupported scheduling "
                    "policy (only IS/PS/SIRO/FCFS/LCFSPR)");
        }
    }

    // -- per-item PH service (alpha, T) per station, routing R --------------
    out.station.resize(S);
    for (std::size_t s = 0; s < S; ++s) out.station[s].type = stype[s];
    std::vector<std::size_t> fsz(S, 0);
    for (std::size_t s = 0; s < S; ++s) {
        const std::size_t st = sn.nodes[qn[s] - 1].station;
        const std::size_t rc0 = ch.retrieval_classes[0][jobin - 1];  // item-1 retrieval class
        const lang::Distrib<T> d0 = sn.service[st - 1][rc0 - 1];
        fsz[s] = lang::dist_to_map(d0).D0.rows();
        if ((stype[s] == RetrievalStationType::SIRO || stype[s] == RetrievalStationType::FCFS) &&
            fsz[s] > 1)
            throw UnsupportedError(
                "cache_retrieval_inputs: SIRO/FCFS retrieval stations need exponential service");
        out.station[s].alpha = Matrix<T>(n, fsz[s], zero);
        out.station[s].sub.assign(n, Matrix<T>(fsz[s], fsz[s], zero));
    }

    out.R.assign(n, Matrix<T>(S + 1, S + 1, zero));
    auto lin = [&](std::size_t node, std::size_t cls) { return (node - 1) * K + (cls - 1); };
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t rcls = ch.retrieval_classes[i][jobin - 1];
        for (std::size_t s = 0; s < S; ++s) {
            const std::size_t st = sn.nodes[qn[s] - 1].station;
            const lang::Distrib<T> d = sn.service[st - 1][rcls - 1];
            const std::vector<T> pie = lang::dist_pie(d);
            const Matrix<T> D0 = lang::dist_to_map(d).D0;
            for (std::size_t a = 0; a < fsz[s]; ++a) {
                out.station[s].alpha(i, a) = pie[a];
                for (std::size_t b = 0; b < fsz[s]; ++b) out.station[s].sub[i](a, b) = D0(a, b);
            }
        }
        for (std::size_t s = 0; s < S; ++s) {
            out.R[i](0, s + 1) = sn.rtnodes(lin(ci, rcls), lin(qn[s], rcls));
            out.R[i](s + 1, 0) = sn.rtnodes(lin(qn[s], rcls), lin(ci, rcls));
            for (std::size_t sp = 0; sp < S; ++sp)
                out.R[i](s + 1, sp + 1) = sn.rtnodes(lin(qn[s], rcls), lin(qn[sp], rcls));
        }
    }

    // -- eta: IS column 0, PS-family columns 1..r ---------------------------
    std::vector<std::size_t> isIdx, psIdx;
    for (std::size_t s = 0; s < S; ++s)
        (stype[s] == RetrievalStationType::IS ? isIdx : psIdx).push_back(s);
    const std::size_t r = psIdx.size();
    out.eta = Matrix<T>(n, r + 1, zero);
    for (std::size_t i = 0; i < n; ++i) {
        // visits = a (I - Pmat)^-1 over the S stations
        Matrix<T> ImP(S, S, zero);
        std::vector<T> a(S, zero);
        for (std::size_t s = 0; s < S; ++s) {
            a[s] = out.R[i](0, s + 1);
            for (std::size_t sp = 0; sp < S; ++sp)
                ImP(s, sp) = T((s == sp ? num_traits<T>::from_int(1) : zero) - out.R[i](s + 1, sp + 1));
        }
        const Matrix<T> F = inverse(ImP);
        std::vector<T> visits(S, zero);
        for (std::size_t s = 0; s < S; ++s) {
            T acc = zero;
            for (std::size_t sp = 0; sp < S; ++sp) acc = T(acc + a[sp] * F(sp, s));
            visits[s] = acc;
        }
        std::vector<T> tau(S, zero);
        for (std::size_t s = 0; s < S; ++s) {
            Matrix<T> arow(1, fsz[s], zero);
            for (std::size_t a = 0; a < fsz[s]; ++a) arow(0, a) = out.station[s].alpha(i, a);
            tau[s] = detail::ph_mean(arow, out.station[s].sub[i]);
        }
        T is_sum = zero;
        for (std::size_t s : isIdx) is_sum = T(is_sum + visits[s] * tau[s]);
        out.eta(i, 0) = is_sum;
        for (std::size_t p = 0; p < r; ++p) out.eta(i, 1 + p) = T(visits[psIdx[p]] * tau[psIdx[p]]);
    }
    return out;
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_CACHE_RETRIEVAL_INPUTS_H
