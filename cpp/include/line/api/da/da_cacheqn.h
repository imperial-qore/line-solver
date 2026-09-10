/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DA_DA_CACHEQN_H
#define LINE_API_DA_DA_CACHEQN_H

/**
 * Decomposition-aggregation driver for integrated cache-queueing models, a port
 * of matlab/src/api/da/da_cacheqn.m.
 *
 * Alternates between (i) solving each cache in isolation given the current
 * per-class arrival rates and (ii) solving the surrounding queueing network with
 * every cache replaced by a class switch that routes according to the current
 * hit/miss probabilities, until the cache arrival rates reach a fixed point
 * (da_fpi under the 1-norm). The cache node is RELABELLED to a ClassSwitch but
 * stays stateful, so it remains a retained node of the stochastic complement --
 * exactly as the reference does. The per-sweep routing rewrite goes straight
 * into `rtnodes`, after which `da_recompute_visits_from_rtnodes` rebuilds `rt`
 * and the visits with the chains held fixed.
 *
 * ARITHMETIC: the miss solve (cache_mva / cache_miss_fpi) and the network solve
 * decide the achievable T; the driver itself is field arithmetic. The tolerance
 * loop is inexact by construction whatever the backend.
 */

#include <cstddef>
#include <limits>
#include <functional>
#include <vector>

#include "line/api/cache/cache_mva.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_ttl_lrua.h"
#include "line/api/cache/cache_miss_fpi.h"
#include "line/api/da/da_cache_isolate.h"
#include "line/api/da/da_fpi.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/matrix.h"

namespace line {
namespace da {

/** The converged isolated-cache inputs, kept for per-item occupancy reporting. */
template <class T>
struct CacheqnInfo {
    std::vector<Matrix<T> > gamma;                     ///< per cache: (n x h)
    std::vector<std::vector<int> > itemcap;            ///< per cache
    std::vector<std::vector<Matrix<T> > > lambda_cache;///< per cache: (u) of n x (h+1)
    std::vector<std::vector<std::vector<Matrix<T> > > > Rcost;
    std::vector<lang::ReplacementStrategy> strat;
};

template <class T>
struct CacheqnResult {
    mva::MvaSolution<T> res;
    Matrix<T> hitprob;   ///< (ncaches x nclasses)
    Matrix<T> missprob;  ///< (ncaches x nclasses)
    int iter = 0;
    CacheqnInfo<T> info;
};

/**
 * Per-item occupancy of every cache from the CONVERGED access factors.
 *
 * This is the EMBEDDED (per-request) law: the stationary law of the cache-content
 * chain seen at request instants, which coincides with the time-stationary one
 * only under PASTA. SolverCTMC reports the time-weighted counterpart instead.
 *
 * RR/FIFO past 10 items is NaN, not the fixed point's marginals: the exact
 * recursion is what makes these a distribution, and an approximation of it is
 * not one. SolverMVA and SolverNC read the same law off their own fixed points,
 * so it is derived here rather than in either analyzer.
 *
 * @return one (n x h+1) matrix per cache, column 0 the miss probability; an
 *         EMPTY matrix where the cache carries no access factors
 */
template <class T>
std::vector<Matrix<T> > da_cacheqn_itemprob(const CacheqnInfo<T>& info) {
    std::vector<Matrix<T> > out(info.gamma.size(), Matrix<T>());
    for (std::size_t ci = 0; ci < info.gamma.size(); ++ci) {
        if (info.gamma[ci].rows() == 0) continue;
        const std::size_t ni = info.gamma[ci].rows();
        const std::size_t hi = info.itemcap[ci].size();
        if (info.strat[ci] == lang::ReplacementStrategy::LRU) {
            std::vector<T> mT(hi);
            for (std::size_t l = 0; l < hi; ++l)
                mT[l] = num_traits<T>::from_int(info.itemcap[ci][l]);
            out[ci] = cache::cache_ttl_lrua(
                info.lambda_cache[ci][0], info.Rcost[ci][0], mT,
                num_traits<T>::from_double(lang::GlobalConstants::FineTol));
        } else if (ni > 10) {
            out[ci] = Matrix<T>(ni, hi + 1,
                                num_traits<T>::from_double(
                                    std::numeric_limits<double>::quiet_NaN()));
        } else {
            out[ci] = cache::cache_prob_erec(info.gamma[ci], info.itemcap[ci]);
        }
    }
    return out;
}

/**
 * @param sn      the model struct, taken by value and mutated (caches relabelled)
 * @param exact   true selects the exact isolated-cache miss, false the
 *                approximation; consulted only when `missfun` is empty
 * @param opt     iteration controls (iter_max, tol)
 * @param netfun  solves the surrounding queueing network on the mutated struct;
 *                its MvaSolution.X is read as the per-class system throughput
 * @param missfun THE ISOLATED-CACHE MISS ALGORITHM, as the reference's `missfun`
 *                handle. Empty selects the MVA pair (`cache_mva` /
 *                `cache_miss_fpi`), which is what SolverMVA passes; SolverNC
 *                supplies `cache_prob_erec` / `cache_miss_spm` instead. The
 *                reference has always taken a handle here -- an earlier version
 *                of this port collapsed it to `exact`, which hardcoded one
 *                solver's pair into a driver both solvers share.
 *
 *                THE FOURTH ARGUMENT IS THE CACHE ITSELF, as the reference's
 *                `missfun(gamma, m, lambda_cache, ch)` passes it. Without it a
 *                handle cannot see the replacement strategy or the access graph
 *                of the cache it is being asked about, so a model with two
 *                caches under different policies would be solved by whichever
 *                one the handle was written for. SolverFLD reads exactly those
 *                two fields; SolverMVA and SolverNC ignore the argument.
 */
template <class T>
CacheqnResult<T> da_cacheqn(qn::NetworkStruct<T> sn, bool exact, const mva::MvaOptions& opt,
                            const std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)>&
                                netfun,
                            const std::function<std::vector<T>(const Matrix<T>&,
                                                               const std::vector<int>&,
                                                               const std::vector<Matrix<T> >&,
                                                               const qn::CacheParam<T>&)>&
                                missfun = nullptr) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;

    std::vector<std::size_t> caches;  // 0-based node indices
    for (std::size_t nd = 0; nd < I; ++nd)
        if (sn.nodes[nd].nodetype == qn::NodeType::Cache) caches.push_back(nd);
    const std::size_t ncaches = caches.size();

    // connmatrix from the ORIGINAL node routing, before the cache rows are
    // rewritten: node ind connects to jnd if any class routes ind -> jnd.
    std::vector<std::vector<bool> > conn(I, std::vector<bool>(I, false));
    for (std::size_t a = 0; a < I; ++a)
        for (std::size_t b = 0; b < I; ++b) {
            if (a == b) continue;  // the cache's read self-switch is not a downstream
            for (std::size_t r = 0; r < K && !conn[a][b]; ++r)
                for (std::size_t s = 0; s < K && !conn[a][b]; ++s)
                    if (sn.rtnodes(a * K + r, b * K + s) > zero) conn[a][b] = true;
        }

    CacheqnResult<T> outr;
    outr.hitprob = Matrix<T>(ncaches, K, zero);
    outr.missprob = Matrix<T>(ncaches, K, zero);
    outr.info.gamma.resize(ncaches);
    outr.info.itemcap.resize(ncaches);
    outr.info.lambda_cache.resize(ncaches);
    outr.info.Rcost.resize(ncaches);
    outr.info.strat.resize(ncaches);

    // seed cache arrival rates and relabel the caches as class switches
    std::vector<T> lambda0(K, zero);
    for (std::size_t ci = 0; ci < ncaches; ++ci) {
        const qn::CacheParam<T>& ch = sn.nodeparam.at(caches[ci] + 1);
        for (std::size_t r = 0; r < ch.hitclass.size(); ++r)
            if (ch.hitclass[r] != 0)
                lambda0[r] = num_traits<T>::from_double(0.5);  // deterministic seed (no rand)
        sn.nodes[caches[ci]].nodetype = qn::NodeType::ClassSwitch;
    }

    da::FpiOptions fopt;
    fopt.iter_max = opt.iter_max;
    fopt.iter_tol = opt.tol;
    // The reference stops on the 1-norm of the increment; this da_fpi uses the
    // max-norm. Both converge to the same fixed point (the point is what the
    // metrics read), so only the iteration count can differ, which is not
    // compared.

    mva::MvaSolution<T> lastres;
    auto sweep = [&](const std::vector<T>& x,
                     std::size_t) -> std::pair<std::vector<T>, std::vector<T> > {
        std::vector<T> lambda = x;
        for (std::size_t ci = 0; ci < ncaches; ++ci) {
            const std::size_t ind = caches[ci];
            const qn::CacheParam<T>& ch = sn.nodeparam.at(ind + 1);

            da::CacheParam<T> dch;
            dch.itemcap = ch.itemcap;
            dch.nitems = ch.nitems;
            dch.pread = ch.pread;
            dch.accost = ch.accost;
            const da::CacheIsolateResult<T> iso = da::da_cache_isolate(dch, lambda);
            outr.info.gamma[ci] = iso.gamma;
            outr.info.itemcap[ci] = ch.itemcap;
            outr.info.lambda_cache[ci] = iso.lambda_cache;
            outr.info.Rcost[ci] = iso.Rcost;
            outr.info.strat[ci] = ch.replacestrat;

            // per-class miss rate of the isolated cache
            std::vector<T> missrate(K, zero);
            if (missfun) {
                const std::vector<T> mr = missfun(iso.gamma, ch.itemcap, iso.lambda_cache, ch);
                for (std::size_t v = 0; v < mr.size() && v < K; ++v) missrate[v] = mr[v];
            } else if (exact) {
                const cache::CacheMvaResult<T> cm = cache::cache_mva(iso.gamma, ch.itemcap);
                const std::size_t n = cm.pij.rows(), h = cm.pij.cols();
                std::vector<T> pmiss(n, zero);
                for (std::size_t k = 0; k < n; ++k) {
                    T s = zero;
                    for (std::size_t l = 0; l < h; ++l) s = T(s + cm.pij(k, l));
                    pmiss[k] = T(one - s);
                }
                for (std::size_t v = 0; v < iso.lambda_cache.size() && v < K; ++v) {
                    T acc = zero;
                    for (std::size_t k = 0; k < n; ++k)
                        acc = T(acc + iso.lambda_cache[v](k, 0) * pmiss[k]);
                    missrate[v] = acc;
                }
            } else {
                // cache_miss_fpi wants lambda as the (u x n) first-list slice.
                const std::size_t u = iso.lambda_cache.size();
                const std::size_t n = ch.nitems;
                Matrix<T> lam_un(u, n, zero);
                for (std::size_t v = 0; v < u; ++v)
                    for (std::size_t k = 0; k < n; ++k) lam_un(v, k) = iso.lambda_cache[v](k, 0);
                const cache::CacheMissResult<T> mr =
                    cache::cache_miss_fpi(iso.gamma, ch.itemcap, lam_un);
                for (std::size_t v = 0; v < mr.MU.size() && v < K; ++v) missrate[v] = mr.MU[v];
            }

            for (std::size_t r = 0; r < K; ++r) {
                if (lambda[r] > zero) {
                    outr.missprob(ci, r) = T(missrate[r] / lambda[r]);
                    outr.hitprob(ci, r) = T(one - outr.missprob(ci, r));
                } else {
                    outr.missprob(ci, r) = zero;
                    outr.hitprob(ci, r) = zero;
                }
            }

            // rewrite the cache-as-classswitch routing into rtnodes: send the
            // hit/miss mass to EVERY connected node, exactly as the reference. The
            // row over-sums when there is more than one downstream node; the
            // reference does not normalise it, and the surplus is absorbed by the
            // pass-through classes (which each queue routes onward to its own
            // successors) at zero residence.
            for (std::size_t r = 0; r < ch.hitclass.size(); ++r) {
                if (ch.hitclass[r] == 0) continue;  // class r does not read this cache
                const std::size_t hc = ch.hitclass[r] - 1, mc = ch.missclass[r] - 1;
                for (std::size_t col = 0; col < I * K; ++col) sn.rtnodes(ind * K + r, col) = zero;
                for (std::size_t jnd = 0; jnd < I; ++jnd) {
                    if (!conn[ind][jnd]) continue;
                    sn.rtnodes(ind * K + r, jnd * K + hc) = outr.hitprob(ci, r);
                    sn.rtnodes(ind * K + r, jnd * K + mc) = outr.missprob(ci, r);
                }
            }
        }
        // The cache switch forces classes onto nodes that never carried them in
        // the base model, where they have no onward routing. The reference's
        // refresh routes every class at a node to that node's physical
        // successors (same class); reproduce it for the classes left unrouted, so
        // a pass-through class exits instead of dead-ending (which would zero its
        // throughput and unbalance the served classes).
        for (std::size_t i = 0; i < I; ++i) {
            if (sn.nodes[i].nodetype == qn::NodeType::Cache ||
                sn.nodes[i].nodetype == qn::NodeType::ClassSwitch)
                continue;  // caches (relabelled) carry their explicit switch row
            std::size_t nsucc = 0;
            for (std::size_t j = 0; j < I; ++j)
                if (conn[i][j]) ++nsucc;
            if (nsucc == 0) continue;  // a Sink has no successors
            const T share = T(one / num_traits<T>::from_int(long(nsucc)));
            for (std::size_t s = 0; s < K; ++s) {
                T rowsum = zero;
                for (std::size_t col = 0; col < I * K; ++col) rowsum = T(rowsum + sn.rtnodes(i * K + s, col));
                if (rowsum > zero) continue;  // already routed
                for (std::size_t j = 0; j < I; ++j)
                    if (conn[i][j]) sn.rtnodes(i * K + s, j * K + s) = share;
            }
        }
        sn.da_recompute_visits_from_rtnodes();

        lastres = netfun(sn);

        // node visits summed over chains
        Matrix<T> nv(I, K, zero);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            for (std::size_t a = 0; a < I; ++a)
                for (std::size_t r = 0; r < K; ++r) nv(a, r) = T(nv(a, r) + sn.nodevisits[c](a, r));

        // update cache arrival rates from the network throughputs
        for (std::size_t ci = 0; ci < ncaches; ++ci) {
            const std::size_t ind = caches[ci];
            const qn::CacheParam<T>& ch = sn.nodeparam.at(ind + 1);
            for (std::size_t r = 0; r < ch.hitclass.size(); ++r) {
                if (ch.hitclass[r] == 0) continue;
                // chain of class r
                std::size_t c = sn.nchains;
                for (std::size_t cc = 0; cc < sn.nchains; ++cc)
                    if (sn.chains[cc][r]) { c = cc; break; }
                if (c == sn.nchains) continue;
                T xsum = zero;
                for (std::size_t k : sn.inchain[c]) xsum = T(xsum + lastres.X[k - 1]);
                const std::size_t rstat = sn.classes[r].refstat;      // 1-based station
                const std::size_t refnode = sn.station_to_node[rstat - 1] - 1;  // 0-based
                const std::size_t refcls =
                    (sn.refclass[c] > 0) ? sn.refclass[c] - 1 : r;  // 0-based
                const T denom = nv(refnode, refcls);
                if (denom > zero) lambda[r] = T(xsum * nv(ind, r) / denom);
            }
        }
        return std::make_pair(lambda, x);
    };

    const da::FpiResult<T> fp = da::da_fpi<T>(sweep, lambda0, fopt);
    outr.iter = static_cast<int>(fp.iterations);
    outr.res = lastres;
    return outr;
}

}  // namespace da
}  // namespace line

#endif  // LINE_API_DA_DA_CACHEQN_H
