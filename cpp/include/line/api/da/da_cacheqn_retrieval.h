/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DA_DA_CACHEQN_RETRIEVAL_H
#define LINE_API_DA_DA_CACHEQN_RETRIEVAL_H

/**
 * Decomposition-aggregation driver for a CLOSED integrated cache-queueing model
 * whose Cache carries a delayed-hit retrieval system.
 *
 * Templated port of `matlab/src/api/da/da_cacheqn_retrieval.m`.
 *
 * THE FIXED POINT. The read rate reaching the cache depends on the network's
 * throughput, which depends on how the cache splits that read rate between hits
 * and fetches. One sweep therefore: solves the ISOLATED cache at the current
 * read rate, rewrites the routing so the cache behaves as a class switch with
 * that split, solves the network, and reads the new read rate back off the
 * converged visits. `da_fpi` drives it under the 1-norm.
 *
 * WHAT MAKES THIS DIFFERENT FROM `da_cacheqn`, and it is not the fixed point.
 * A miss here does not leave the cache: it becomes a per-item RETRIEVAL CLASS
 * that must be FETCHED through a backend station and returned. Two consequences
 * are built into the driver rather than into the caller:
 *
 *  - THE FETCH STATION IS GIVEN A LOAD-DEPENDENT RATE, and the shape of it is
 *    the coupon-collector correction `alpha(k) = k / (n_eff (1 - (1-1/n_eff)^k))`
 *    with `n_eff = nitems - total capacity`. With k fetches in flight, the
 *    number of DISTINCT items among them is below k, because two misses may
 *    want the same item and one fetch serves both. alpha is the ratio, so the
 *    station is faster than k independent fetches would be. That is what makes
 *    the model load dependent, and it is why the caller's network solver must
 *    handle `lldscaling`.
 *  - THE CACHE BECOMES A ClassSwitch for the network solve, its routing rewritten
 *    each sweep: the read class goes to the hit class with probability `hp`, and
 *    to item i's retrieval class at the fetch node with probability
 *    `pread(i) pi0(i)`; the retrieval class returns to the miss class.
 *
 * THE MISS ALGORITHM IS FIXED HERE, deliberately, unlike `da_cacheqn` whose
 * `missfun` is a handle. BOTH callers -- SolverMVA and SolverNC -- use
 * `cache_miss_fpi`, and both report `method = 'fpi'` on this path. Nothing
 * varies between them except the NETWORK solver, so `netfun` is the only handle
 * and adding a `missfun` would be speculative generality. Checked against both
 * call sites before fixing the type.
 *
 * DELAYED HITS FOLD INTO MISS on this path: `delayedprob` is returned as zero
 * and the reported hit is P(item cached). The open analyzer
 * (`solver_nc_retrieval_analyzer`) is the one that separates the three.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "line/api/cache/cache_miss_fpi.h"
#include "line/api/da/da_cache_isolate.h"
#include "line/api/da/da_fpi.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace da {

/** What the delayed-hit decomposition returns. */
template <class T>
struct CacheqnRetrievalResult {
    mva::MvaSolution<T> res;
    std::vector<T> hitprob;      ///< (K) P(item cached), read class only
    std::vector<T> missprob;     ///< (K)
    std::vector<T> delayedprob;  ///< (K) zero on this path; see the header
    std::size_t iter = 0;
    qn::NetworkStruct<T> sn;     ///< the mutated struct, cache relabelled
};

/**
 * Port of `da_cacheqn_retrieval`.
 *
 * @param sn     the model struct, taken by value and MUTATED (the cache becomes
 *               a ClassSwitch and the fetch station gains `lldscaling`)
 * @param netfun solves the surrounding queueing network on the mutated struct
 * @param opt    iteration controls
 */
template <class T>
CacheqnRetrievalResult<T> da_cacheqn_retrieval(
    qn::NetworkStruct<T> sn,
    const std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)>& netfun,
    const mva::MvaOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "da_cacheqn_retrieval requires transcendental arithmetic: it alternates two "
                  "tolerance-stopped solves");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;

    std::vector<std::size_t> caches;
    for (std::size_t nd = 0; nd < I; ++nd)
        if (sn.nodes[nd].nodetype == qn::NodeType::Cache) caches.push_back(nd);
    if (caches.size() != 1)
        throw UnsupportedError("da_cacheqn_retrieval: requires exactly one Cache node");
    const std::size_t ci = caches[0];
    const qn::CacheParam<T> ch = sn.nodeparam.at(ci + 1);

    if (ch.retrieval_queues.empty())
        throw UnsupportedError("da_cacheqn_retrieval: the Cache carries no retrieval system");
    const std::size_t readClass = ch.retrieval_queues.begin()->first + 1;  // 1-based
    const std::vector<std::size_t>& queueNodes = ch.retrieval_queues.begin()->second;
    if (queueNodes.size() != 1)
        throw UnsupportedError(
            "da_cacheqn_retrieval: currently supports a single-station (single-backend) retrieval "
            "system; this cache fetches through " + std::to_string(queueNodes.size()));
    const std::size_t fetchNode = queueNodes[0];           // 1-based
    const std::size_t fetchStation = sn.nodes[fetchNode - 1].station;
    if (fetchStation == 0)
        throw UnsupportedError("da_cacheqn_retrieval: the retrieval node is not a station");

    const std::size_t nitems = ch.nitems;
    long totcap = 0;
    for (int c : ch.itemcap) totcap += c;
    const long n_eff_l = std::max<long>(1, static_cast<long>(nitems) - totcap);
    const double n_eff = static_cast<double>(n_eff_l);

    double Npop_d = 0.0;
    for (const qn::JobClass& c : sn.classes)
        if (std::isfinite(c.population)) Npop_d += c.population;
    const std::size_t Npop = static_cast<std::size_t>(std::llround(Npop_d));
    if (Npop == 0)
        throw UnsupportedError(
            "da_cacheqn_retrieval: the closed delayed-hit decomposition needs a positive closed "
            "population");

    // THE COUPON-COLLECTOR CORRECTION. With k fetches in flight the number of
    // DISTINCT items among them is n_eff (1 - (1-1/n_eff)^k), below k, because
    // two misses may want the same item and one fetch serves both. alpha is the
    // ratio, so the fetch station serves k in-flight misses faster than k
    // independent fetches would.
    std::vector<T> alpha(std::max<std::size_t>(1, Npop), one);
    for (std::size_t k = 1; k <= Npop; ++k) {
        const double d = n_eff * (1.0 - std::pow(1.0 - 1.0 / n_eff, static_cast<double>(k)));
        alpha[k - 1] = num_traits<T>::from_double(static_cast<double>(k) / d);
    }
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        std::vector<T>& lld = sn.stations[i].lldscaling;
        if (lld.size() < Npop) lld.resize(Npop, one);
    }
    sn.stations[fetchStation - 1].lldscaling = alpha;

    // The cache is solved in isolation; for the NETWORK it is a class switch.
    sn.nodes[ci].nodetype = qn::NodeType::ClassSwitch;

    std::vector<T> pread = ch.pread.at(readClass - 1);
    T psum = zero;
    for (const T& v : pread) psum += v;
    if (!(psum > zero))
        throw UnsupportedError("da_cacheqn_retrieval: the read class has no item popularity");
    for (T& v : pread) v = T(v / psum);

    CacheqnRetrievalResult<T> outr;
    outr.hitprob.assign(K, zero);
    outr.missprob.assign(K, zero);
    outr.delayedprob.assign(K, zero);

    da::FpiOptions fpopt;
    fpopt.iter_max = static_cast<std::size_t>(opt.iter_max);
    fpopt.iter_tol = opt.iter_tol;

    std::vector<T> x0(K, zero);
    x0[readClass - 1] = one;  // seed

    const auto sweep = [&](const std::vector<T>& x,
                           std::size_t) -> std::pair<std::vector<T>, std::vector<T> > {
        std::vector<T> lambda = x;

        // 1. the isolated cache at the current read rate. `da_cache_isolate`
        //    takes its own lightweight parameter struct, not the model's.
        da::CacheParam<T> dch;
        dch.itemcap = ch.itemcap;
        dch.nitems = ch.nitems;
        // setRetrievalSystem MINTS one retrieval class per item, so the class
        // count grows past the length of `pread`. An EMPTY row is the
        // documented marker for "this class does not read the cache", which is
        // exactly what a retrieval class is, so the tail is padded rather than
        // the rate vector truncated.
        dch.pread = ch.pread;
        dch.pread.resize(lambda.size());
        dch.accost = ch.accost;
        const da::CacheIsolateResult<T> iso = da::da_cache_isolate(dch, lambda);
        const std::size_t u = iso.lambda_cache.size();
        Matrix<T> lam_un(u, nitems, zero);
        for (std::size_t v = 0; v < u; ++v)
            for (std::size_t k = 0; k < nitems; ++k) lam_un(v, k) = iso.lambda_cache[v](k, 0);
        const cache::CacheMissResult<T> mf =
            cache::cache_miss_fpi(iso.gamma, ch.itemcap, lam_un);
        std::vector<T> pi0(nitems, zero);
        for (std::size_t k = 0; k < nitems && k < mf.pi0.size(); ++k) pi0[k] = mf.pi0[k];
        T nonhit = zero;
        for (std::size_t k = 0; k < nitems; ++k) nonhit += T(pread[k] * pi0[k]);
        const T hp = T(one - nonhit);

        // 2. rewrite the cache-as-classswitch routing: hit straight to the hit
        //    class, miss into item i's retrieval class at the fetch node, and
        //    the retrieval class back out as the miss class.
        const std::size_t r = readClass - 1;
        for (std::size_t col = 0; col < I * K; ++col) sn.rtnodes(ci * K + r, col) = zero;
        const std::size_t hc = ch.hitclass.at(r), mc = ch.missclass.at(r);
        if (hc == 0 || mc == 0)
            throw UnsupportedError(
                "da_cacheqn_retrieval: the read class has no hit or miss class");
        sn.rtnodes(ci * K + r, ci * K + (hc - 1)) = hp;
        for (std::size_t i = 0; i < nitems; ++i) {
            const std::size_t rcls =
                (i < ch.retrieval_classes.size() && r < ch.retrieval_classes[i].size())
                    ? ch.retrieval_classes[i][r]
                    : 0;
            if (rcls == 0) continue;
            sn.rtnodes(ci * K + r, (fetchNode - 1) * K + (rcls - 1)) = T(pread[i] * pi0[i]);
            for (std::size_t col = 0; col < I * K; ++col)
                sn.rtnodes((fetchNode - 1) * K + (rcls - 1), col) = zero;
            sn.rtnodes((fetchNode - 1) * K + (rcls - 1), ci * K + (mc - 1)) = one;
            // drop the unused Cache -> Retrieval edge
            for (std::size_t col = 0; col < I * K; ++col)
                sn.rtnodes(ci * K + (rcls - 1), col) = zero;
        }
        sn.da_recompute_visits_from_rtnodes();

        // 3. the network solve, and the read rate read back off the visits
        outr.res = netfun(sn);

        Matrix<T> nv(I, K, zero);
        for (std::size_t c = 0; c < sn.nchains; ++c)
            for (std::size_t a = 0; a < I; ++a)
                for (std::size_t k = 0; k < K; ++k)
                    nv(a, k) = T(nv(a, k) + sn.nodevisits[c](a, k));
        std::size_t chain = sn.nchains;
        for (std::size_t c = 0; c < sn.nchains; ++c)
            if (sn.chains[c][r]) chain = c;
        T denom = zero;
        if (chain < sn.nchains) {
            const std::size_t refnode = sn.node_of_station(sn.classes[r].refstat);
            denom = sn.refclass[chain] > 0 ? nv(refnode - 1, sn.refclass[chain] - 1)
                                           : nv(refnode - 1, r);
        }
        T Xr = zero;
        if (chain < sn.nchains)
            for (std::size_t k = 0; k < K; ++k)
                if (sn.chains[chain][k] && k < outr.res.X.size()) Xr += outr.res.X[k];
        if (denom > zero) lambda[r] = T(Xr * nv(ci, r) / denom);

        outr.hitprob[r] = hp;
        outr.missprob[r] = nonhit;
        outr.delayedprob[r] = zero;
        return std::make_pair(lambda, x);
    };

    const da::FpiResult<T> fr = da::da_fpi<T>(sweep, x0, fpopt);
    outr.iter = fr.iterations;
    outr.sn = sn;
    return outr;
}

}  // namespace da
}  // namespace line

#endif  // LINE_API_DA_DA_CACHEQN_RETRIEVAL_H
