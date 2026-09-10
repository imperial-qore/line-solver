/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_CACHE_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_CACHE_H

/**
 * The non-reentrant cache analyzer: a Source-Cache-Sink model.
 *
 * Templated port of `matlab/src/solvers/MVA/solver_mva_cache_analyzer.m`.
 *
 * A read of class r reaches the cache at the Source rate weighted by the class
 * read distribution, and completes as a HIT (switching to `hitclass(r)`) or a
 * MISS (switching to `missclass(r)`). The analyzer computes the per-item
 * occupancy of the cache and from it the miss rate of each read class, then
 * splits the Source throughput between the hit and the miss class.
 *
 * The item occupancy comes from the cache algorithm family, keyed on the
 * replacement policy:
 *
 *   RR, FIFO   `exact`  -> cache_mva   (the exact product-form recursion)
 *              default  -> cache_prob_fpi (the fixed-point approximation)
 *   LRU        -> cache_ttl_lrua (the LRU-A characteristic-time approximation)
 *   HLRU       -> cache_ttl_hlru (the h-LRU / LRU(m) characteristic time)
 *
 * The RR/FIFO inputs -- the access factors `gamma` -- are built by
 * `da_cache_isolate`, which spreads each class rate over the items through its
 * read distribution and attaches the list-to-list access-cost matrices;
 * `qn::CacheParam` carries `pread` and `accost` in exactly the shape it consumes.
 *
 * THE TTL FAMILY IS SINGLE-STREAM. `cache_ttl_lrua` and `cache_ttl_hlru` read
 * only the FIRST user stream in the reference (`lambda(1,i,j)`, `R{1,i}`), so
 * the port passes the first reading class's per-item rates and access graphs,
 * which is what those two consume. A model with more than one reading class is
 * therefore an LRU/h-LRU approximation over its first read law, matching the
 * reference exactly rather than silently aggregating.
 *
 * WHAT IS REFUSED BY NAME. The exact solution of any policy other than RR/FIFO,
 * SFIFO, CLIMB and QLRU, and the marked-MAP LRU sub-branch, follow the
 * reference's own `line_error`s. The marked-MAP branch needs an MMAP Source
 * with per-item marks, which the builder cannot yet construct.
 *
 * Arithmetic: RR/FIFO exact (cache_mva) is field arithmetic and runs at any T;
 * the FPI approximation evaluates a fixed point on a tolerance and needs a
 * transcendental T, so it refuses under Rational.
 */

#include <cmath>
#include <vector>

#include "line/api/cache/cache_mva.h"
#include "line/api/cache/cache_prob_fpi.h"
#include "line/api/cache/cache_ttl_hlru.h"
#include "line/api/cache/cache_ttl_lrua.h"
#include "line/api/da/da_cache_isolate.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mva {

/** What the cache analyzer reports beyond the [Q,U,R,T] block. */
template <class T>
struct CacheResult {
    MvaSolution<T> sol;
    std::vector<T> hitprob;   ///< per class, NaN where the class does not read
    std::vector<T> missprob;
    /**
     * (R x h) per-list hit fractions, access-weighted over items. NaN unless the
     * EXACT recursion ran: only it produces a genuine per-list occupancy, which
     * is why the reference gates this on `pijlist` rather than on `pij`.
     */
    Matrix<T> hitproblist;
    /** (n x h+1) per-item occupancy, column 0 = miss. EMPTY = not computed. */
    Matrix<T> itemprob;
    std::string actualmethod;
};

/** Port of `solver_mva_cache_analyzer.m` for a Source-Cache-Sink model. */
template <class T>
CacheResult<T> solver_mva_cache_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::NodeType;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, R = L.nclasses;

    // the Cache node and its parameters
    std::size_t cnode = 0;
    for (std::size_t i = 0; i < L.nof_nodes(); ++i)
        if (L.nodes[i].nodetype == NodeType::Cache) {
            cnode = i + 1;
            break;
        }
    if (cnode == 0) throw InputError("solver_mva_cache_analyzer: the model has no Cache node");
    const auto it = L.nodeparam.find(cnode);
    if (it == L.nodeparam.end())
        throw InputError("solver_mva_cache_analyzer: the Cache node has no parameters");
    const qn::CacheParam<T>& ch = it->second;

    const std::size_t nan_flag = 0;
    (void)nan_flag;
    const std::size_t n = ch.nitems, h = ch.itemcap.size();

    // the Source, and its per-class read rate
    std::size_t src = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == NodeType::Source) src = i + 1;
    if (src == 0) throw InputError("solver_mva_cache_analyzer: a non-reentrant cache needs a Source");
    std::vector<T> sourceRate(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        const double v = num_traits<T>::to_double(L.rates(src - 1, r));
        if (std::isfinite(v)) sourceRate[r] = L.rates(src - 1, r);
    }

    // isolate the cache: gamma over the items, plus the per-class per-item
    // arrival rates and access graphs the TTL family needs
    da::CacheParam<T> dch;
    dch.nitems = n;
    dch.itemcap = ch.itemcap;
    dch.pread = ch.pread;
    dch.accost = ch.accost;
    const da::CacheIsolateResult<T> iso = da::da_cache_isolate(dch, sourceRate);

    // per-item occupancy pij: column 0 is the miss column
    Matrix<T> pij;
    // The genuine per-list occupancy (n x h), set ONLY by the exact branch. The
    // reference gates both per-list and per-item reporting on this being set.
    Matrix<T> pijlist;
    std::string method;
    if (ch.replacestrat == lang::ReplacementStrategy::RR ||
        ch.replacestrat == lang::ReplacementStrategy::FIFO) {
        if (opt.method == "exact") {
            const cache::CacheMvaResult<T> cm = cache::cache_mva(iso.gamma, ch.itemcap);
            // prepend the miss column, |1 - sum_l pij|, as the reference does
            pij = Matrix<T>(n, h + 1, zero);
            for (std::size_t k = 0; k < n; ++k) {
                T s = zero;
                for (std::size_t l = 0; l < h; ++l) {
                    pij(k, l + 1) = cm.pij(k, l);
                    s += cm.pij(k, l);
                }
                pij(k, 0) = num_abs(T(num_traits<T>::from_int(1) - s));
            }
            pijlist = cm.pij;
            method = "exact";
        } else {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_mva_cache_analyzer: the RR/FIFO fixed-point approximation stops on a "
                    "tolerance and needs a transcendental T; use method 'exact' under Rational");
            } else {
                pij = cache::cache_prob_fpi(iso.gamma, ch.itemcap);
                method = "fpi";
            }
        }
    } else if (ch.replacestrat == lang::ReplacementStrategy::LRU ||
               ch.replacestrat == lang::ReplacementStrategy::HLRU) {
        // The TTL characteristic-time approximations, over the FIRST reading
        // class's stream (see the header note). A marked-MAP Source would take
        // the reference's cache_ttl_lrum_map branch, which is not reachable
        // here, so a plain rate stream is the only case.
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "solver_mva_cache_analyzer: the LRU / h-LRU characteristic-time approximation "
                "solves a fixed point on a tolerance and needs a transcendental T");
        } else {
            std::size_t rd = ch.pread.size();
            for (std::size_t v = 0; v < ch.pread.size(); ++v)
                if (!ch.pread[v].empty()) {
                    rd = v;
                    break;
                }
            if (rd == ch.pread.size())
                throw InputError("solver_mva_cache_analyzer: no class reads the cache");
            if (ch.replacestrat == lang::ReplacementStrategy::HLRU) {
                // cache_ttl_hlru wants the (users x n) request-rate matrix and
                // sums over the user rows: MATLAB passes the full lambda and
                // does lam = sum_v lambda(v,:,min(2,h+1)). Build one row per
                // class (empty-pread rows are zero and contribute nothing).
                const std::size_t u = ch.pread.size();
                Matrix<T> lam_h(u, n, zero);
                for (std::size_t v = 0; v < u; ++v)
                    for (std::size_t k = 0; k < n; ++k)
                        lam_h(v, k) = iso.lambda_cache[v](k, 0);
                pij = cache::cache_ttl_hlru(lam_h, ch.itemcap);
                method = "ttl";
            } else {
                // cache_ttl_lrua reads user 1 literally (lambda(1,...), R{1,...});
                // MATLAB passes the full array but only ever indexes the first
                // user row, so the LRU stream is class index 0, not `rd`.
                std::vector<T> mT(ch.itemcap.size());
                for (std::size_t l = 0; l < ch.itemcap.size(); ++l)
                    mT[l] = num_traits<T>::from_int(ch.itemcap[l]);
                pij = cache::cache_ttl_lrua(
                    iso.lambda_cache[0], iso.Rcost[0], mT,
                    num_traits<T>::from_double(lang::GlobalConstants::FineTol));
                method = "ttl";
            }
        }
    } else {
        throw UnsupportedError(
            "solver_mva_cache_analyzer: replacement policy " +
            std::to_string(static_cast<int>(ch.replacestrat)) +
            " has no MVA cache analyzer; RR, FIFO, LRU and h-LRU are supported");
    }

    // per-class miss rate: the class request rate on each item times the item's
    // miss probability, summed over items
    std::vector<T> missRate(R, zero);
    for (std::size_t v = 0; v < R; ++v) {
        if (v >= ch.pread.size() || ch.pread[v].empty()) continue;
        for (std::size_t k = 0; k < n; ++k)
            missRate[v] += T(sourceRate[v] * ch.pread[v][k] * pij(k, 0));
    }

    // assemble outputs. A Cache is not a queue, so Q/U/R at it are zero; the
    // Source carries its own throughput, and the hit and miss classes carry the
    // split of it.
    MvaSolution<T> out;
    out.Q = Matrix<T>(M, R, zero);
    out.U = Matrix<T>(M, R, zero);
    out.R = Matrix<T>(M, R, zero);
    out.Tp = Matrix<T>(M, R, zero);
    out.C.assign(R, zero);
    out.X.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r) out.Tp(src - 1, r) = sourceRate[r];
    for (std::size_t r = 0; r < R; ++r) {
        if (r >= ch.hitclass.size()) continue;
        const std::size_t hc = ch.hitclass[r], mc = ch.missclass[r];
        if (hc == 0 || mc == 0) continue;
        out.X[mc - 1] = T(out.X[mc - 1] + missRate[r]);
        out.X[hc - 1] = T(out.X[hc - 1] + T(sourceRate[r] - missRate[r]));
    }

    CacheResult<T> res;
    res.sol = out;
    res.sol.method = method;
    res.sol.lG = std::numeric_limits<double>::quiet_NaN();
    res.actualmethod = method;
    const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
    res.hitprob.assign(R, nan);
    res.missprob.assign(R, nan);
    for (std::size_t r = 0; r < R; ++r) {
        if (r >= ch.pread.size() || ch.pread[r].empty()) continue;
        if (!(sourceRate[r] > zero)) continue;
        res.missprob[r] = T(missRate[r] / sourceRate[r]);
        res.hitprob[r] = T(num_traits<T>::from_int(1) - res.missprob[r]);
    }

    // per-list hit probabilities, access-weighted, ONLY where the exact
    // algorithm produced a genuine per-list occupancy
    res.hitproblist = Matrix<T>(R, h, nan);
    if (pijlist.rows() > 0) {
        for (std::size_t v = 0; v < R; ++v) {
            if (v >= ch.pread.size() || ch.pread[v].empty()) continue;
            for (std::size_t l = 0; l < h; ++l) {
                T acc = zero;
                for (std::size_t k = 0; k < n; ++k) acc = T(acc + ch.pread[v][k] * pijlist(k, l));
                res.hitproblist(v, l) = acc;
            }
        }
    }

    // per-item occupancy. RR/FIFO under the FIXED POINT has no per-item law, so
    // the reference re-runs the exact recursion for it -- and refuses past 10
    // items, where that recursion is not affordable, reporting NaN rather than
    // the fixed point's marginals, which are not a per-item distribution.
    if (pijlist.rows() > 0) {
        res.itemprob = pij;
    } else if (pij.cols() == h + 1) {
        if (ch.replacestrat == lang::ReplacementStrategy::RR ||
            ch.replacestrat == lang::ReplacementStrategy::FIFO) {
            if (n > 10) {
                res.itemprob = Matrix<T>(n, h + 1, nan);
            } else {
                const cache::CacheMvaResult<T> cm = cache::cache_mva(iso.gamma, ch.itemcap);
                res.itemprob = Matrix<T>(n, h + 1, zero);
                for (std::size_t k = 0; k < n; ++k) {
                    T s = zero;
                    for (std::size_t l = 0; l < h; ++l) {
                        res.itemprob(k, l + 1) = cm.pij(k, l);
                        s = T(s + cm.pij(k, l));
                    }
                    res.itemprob(k, 0) = num_abs(T(num_traits<T>::from_int(1) - s));
                }
            }
        } else {
            res.itemprob = pij;
        }
    }
    return res;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_CACHE_H
