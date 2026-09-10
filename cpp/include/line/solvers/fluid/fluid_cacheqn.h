/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_CACHEQN_H
#define LINE_SOLVERS_FLUID_FLUID_CACHEQN_H

/**
 * The INTEGRATED caching-queueing network under the fluid solver: ports of
 * `solver_fld_cacheqn_analyzer.m` (steady state) and `solver_fld_cacheqn_tran.m`
 * (transient cache trajectory).
 *
 * WHY THIS IS A FIXED POINT AND NOT TWO SOLVES IN SEQUENCE. A cache embedded in
 * a queueing network splits the flow it receives into a hit stream and a miss
 * stream, and that split IS the routing of the surrounding network. The flow
 * reaching the cache is in turn a network throughput. So neither side can be
 * evaluated first: the two are alternated until the cache arrival rates stop
 * moving. Between passes the routing has changed, so the visits must be rebuilt
 * from the rewritten `rtnodes` -- the reference re-enters `sn_refresh_visits`,
 * and this port's `da_cacheqn` calls `da_recompute_visits_from_rtnodes` in the
 * same place. Reusing the previous pass's visits would leave every downstream
 * metric computed against a routing the model no longer has.
 *
 * WHAT THIS FILE SUPPLIES, and it is only two things, exactly as the MVA and NC
 * siblings do: the isolated-cache MISS ALGORITHM and the NETWORK SOLVER.
 *
 *   miss     `cache_miss_rmf`, the refined (1/N-accurate) mean field, or the
 *            position-resolved `cache_miss_fifo_rmf` / `cache_miss_sfifo_rmf`
 *            where the policy or the access graph calls for them
 *            (SolverMVA uses `cache_mva` / `cache_miss_fpi`, SolverNC
 *             `cache_prob_erec` / `cache_miss_spm`)
 *   network  the fluid `matrix` method, as the reference's `netsolve` calls
 *            `solver_fluid_matrix`
 *
 * WHICH REPLACEMENT POLICIES HAVE A FLUID MODEL AT ALL. Only those with a
 * drift. RANDOM(m) has one; FIFO(m) shares its STEADY STATE with RANDOM(m) on
 * the linear access graph (Gast15 Thm 1, pi_FIFO(m) = pi_RAND(m)) but not its
 * transient; strict FIFO(m) has its own position-resolved drift. LRU, HLRU,
 * CLIMB and QLRU have none: the characteristic-time (FPI) approximation those
 * are solved with is not a fluid method, and substituting it here would return
 * a number produced by a different model under the fluid solver's name. They
 * are refused, as the reference refuses them.
 *
 * FLUID IS AN APPROXIMATION. The drift is the mean-field limit of the queueing
 * network and the cache miss rates carry a 1/N correction; neither is exact at
 * finite population, so nothing here reproduces MVA or NC to integrator
 * tolerance. What IS exact is the structure: hit and miss probability sum to
 * one, and a cache that holds every item never misses.
 *
 * DOUBLE ONLY, and gated with `if constexpr` rather than the runtime check
 * `solver_fluid.h` uses: the fluid result type is built on `Matrix<double>`
 * while `da_cacheqn` is templated on T, so the two only meet when T is double.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/cache/cache_miss_pos_rmf.h"
#include "line/api/cache/cache_miss_rmf.h"
#include "line/api/da/da_cacheqn.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"
#include "line/solvers/fluid/fluid_moments.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** What the steady analyzer returns: the fluid metrics plus the converged split. */
template <class T>
struct FluidCacheqnSolution {
    FluidSolution sol;
    Matrix<T> hitprob;   ///< (ncaches x nclasses), cache order as the node scan
    Matrix<T> missprob;  ///< (ncaches x nclasses)
    int iter = 0;        ///< decomposition sweeps, not integrator steps
    /**
     * The struct whose cache self-switch carries the CONVERGED split rather than
     * the offered one, from which the runner derives ArvR and ResidT. Same
     * contract as the MVA and NC siblings.
     */
    qn::NetworkStruct<T> refreshed;
};

/** One cache's transient, the per-cache slice of the reference's 3-D outputs. */
template <class T>
struct FluidCacheqnTranCache {
    std::size_t node = 0;    ///< 0-based Cache node index
    std::vector<T> t;        ///< time grid of THIS cache
    Matrix<T> hitprob_t;     ///< (nclasses x nt)
    Matrix<T> missprob_t;    ///< (nclasses x nt)
    std::vector<T> arate;    ///< (nclasses) converged arrival rate at the cache
    Matrix<T> xocc;          ///< (nitems*(h+1) x nt) DDPP occupancy trajectory
};

namespace detail {

/** 0-based Cache node indices, in the order `find(sn.nodetype == Cache)` gives. */
template <class T>
std::vector<std::size_t> fluid_cacheqn_nodes(const qn::NetworkStruct<T>& sn) {
    std::vector<std::size_t> caches;
    for (std::size_t nd = 0; nd < sn.nodes.size(); ++nd)
        if (sn.nodes[nd].nodetype == qn::NodeType::Cache) caches.push_back(nd);
    return caches;
}

/**
 * `accost_is_linear` of the reference: true when every per-(user,item) access
 * graph is the linear chain (miss -> list 1, hit in list a -> list a+1,
 * self-loop on the top list), which is the graph `da_cache_isolate` installs
 * when `accost` is empty.
 */
template <class T>
bool fluid_cacheqn_accost_is_linear(const std::vector<std::vector<Matrix<T> > >& accost,
                                    std::size_t h) {
    if (accost.empty()) return true;
    const T one = num_traits<T>::from_int(1);
    const double eps = 1e-9;
    Matrix<T> lin(h + 1, h + 1, num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < h; ++a) lin(a, a + 1) = one;
    lin(h, h) = one;
    for (std::size_t v = 0; v < accost.size(); ++v)
        for (std::size_t k = 0; k < accost[v].size(); ++k) {
            const Matrix<T>& g = accost[v][k];
            if (g.rows() == 0) continue;
            if (g.rows() != h + 1 || g.cols() != h + 1) return false;
            for (std::size_t a = 0; a <= h; ++a)
                for (std::size_t b = 0; b <= h; ++b)
                    if (std::fabs(num_traits<T>::to_double(g(a, b)) -
                                  num_traits<T>::to_double(lin(a, b))) > eps)
                        return false;
        }
    return true;
}

/**
 * `miss_isolated` of the reference: the drift that answers THIS cache.
 *
 * RANDOM(m) always takes `cache_miss_rmf`, which now honours a declared access
 * graph. FIFO(m) is served from `cache_miss_rmf` on the LINEAR chain, where
 * Gast15 Thm 1 makes the two steady states identical, and from its own
 * position-resolved drift otherwise -- the equality is proved for the chain
 * only. Strict FIFO(m) always takes its own drift; it is a different policy,
 * not a spelling of FIFO(m).
 */
template <class T>
std::vector<T> fluid_cacheqn_miss(const std::vector<int>& m, const Matrix<T>& lam,
                                  const qn::CacheParam<T>& ch) {
    const bool linear = fluid_cacheqn_accost_is_linear(ch.accost, m.size());
    switch (ch.replacestrat) {
        case lang::ReplacementStrategy::RR:
            return cache::cache_miss_rmf(std::vector<T>(), m, lam,
                                         T(num_traits<T>::from_int(10000)), ch.accost)
                .MU;
        case lang::ReplacementStrategy::FIFO:
            if (linear) return cache::cache_miss_rmf(std::vector<T>(), m, lam).MU;
            return cache::cache_miss_fifo_rmf(std::vector<T>(), m, lam, ch.accost).MU;
        case lang::ReplacementStrategy::SFIFO:
            return cache::cache_miss_sfifo_rmf(std::vector<T>(), m, lam, ch.accost).MU;
        default:
            throw UnsupportedError(
                "solver_fld_cacheqn: replacement strategy " +
                std::to_string(static_cast<int>(ch.replacestrat)) +
                " has no drift-based fluid model");
    }
}

/**
 * The admissibility gate, run ONCE before the decomposition starts.
 *
 * The reference selects the drift inside `miss_isolated`, per cache and per
 * sweep, and so does `fluid_cacheqn_miss` above; this gate exists to refuse the
 * policies that have no drift at all BEFORE a sweep runs, which is the honest
 * behaviour -- a policy without a drift cannot become supported halfway through
 * a fixed point.
 */
template <class T>
void fluid_cacheqn_gate(const qn::NetworkStruct<T>& sn, const std::vector<std::size_t>& caches,
                        bool transient) {
    const std::string who =
        transient ? std::string("solver_fld_cacheqn_tran") : std::string("solver_fld_cacheqn_analyzer");
    for (std::size_t ci = 0; ci < caches.size(); ++ci) {
        const std::size_t key = caches[ci] + 1;  // nodeparam is keyed 1-based
        if (sn.nodeparam.count(key) == 0)
            throw InputError(who + ": cache node " + std::to_string(key) +
                             " carries no cache parameters");
        const qn::CacheParam<T>& ch = sn.nodeparam.at(key);
        const std::string at = " (cache node " + std::to_string(key) + ")";
        switch (ch.replacestrat) {
            case lang::ReplacementStrategy::RR:
            case lang::ReplacementStrategy::FIFO:
            case lang::ReplacementStrategy::SFIFO:
                break;
            default:
                // Verbatim the reference's refusal: the FPI characteristic-time
                // approximation those policies use is not a fluid method.
                throw UnsupportedError(
                    who + ": SolverFLD supports only RANDOM(m)/FIFO(m) (refined mean field) and "
                          "strict FIFO(m) (position-resolved mean field) cache replacement; "
                          "strategy " +
                    std::to_string(static_cast<int>(ch.replacestrat)) +
                    " has no drift-based fluid model. Use SolverNC/SolverMVA or SolverLDES for "
                    "this cache" + at);
        }
    }
}

/**
 * The isolated-cache miss algorithm shared by both entry points.
 *
 * `lambda_cache[v]` is (nitems x (h+1)) and repeats the same rate down every
 * list position; `cache_miss_rmf` wants the (users x items) slice, which is
 * column 0. This is the same reshaping the MVA and NC siblings do for their own
 * miss routines.
 */
template <class T>
Matrix<T> fluid_cacheqn_lambda_slice(const std::vector<Matrix<T> >& lambda_cache) {
    const std::size_t u = lambda_cache.size();
    const std::size_t n = (u == 0) ? 0 : lambda_cache[0].rows();
    Matrix<T> lam(u, n, num_traits<T>::from_int(0));
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t k = 0; k < n; ++k) lam(v, k) = lambda_cache[v](k, 0);
    return lam;
}

/**
 * The reference's default initial occupancy: the first m(1) items in list 1, the
 * next m(2) in list 2, and every remaining item outside the cache. Duplicated
 * from `cache_miss_rmf`'s steady path because its transient entry point takes
 * the state as an argument and has no default of its own.
 */
template <class T>
std::vector<T> fluid_cacheqn_default_x0(const std::vector<int>& m, std::size_t nitems) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t h = m.size();
    std::vector<T> x0(nitems * (h + 1), zero);
    std::size_t obj = 0;
    for (std::size_t k = 1; k <= h; ++k)
        for (int j = 0; j < m[k - 1]; ++j) {
            ++obj;
            if (obj <= nitems) x0[cache::cache_miss_rmf_index(obj - 1, k, nitems)] = one;
        }
    for (std::size_t i = obj; i < nitems; ++i) x0[cache::cache_miss_rmf_index(i, 0, nitems)] = one;
    return x0;
}

}  // namespace detail

/**
 * Port of `solver_fld_cacheqn_analyzer.m`.
 *
 * Returns the branch result BEFORE `solver_fluid`'s analyzer-level utilization
 * and response-time correction, exactly as the reference's branch returns what
 * `solver_fluid_matrix` gave it and lets `solver_fluid_analyzer` correct once
 * after the method switch. Applying the correction here as well would apply it
 * twice, and it is not idempotent (it rescales U by a share computed from the
 * U it is given).
 */
template <class T>
FluidCacheqnSolution<T> solver_fld_cacheqn_analyzer(const qn::NetworkStruct<T>& sn,
                                                    const FluidOptions& opt) {
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_fld_cacheqn_analyzer: the fluid network solve integrates its drift with LSODA "
            "and the cache drift with an adaptive Rosenbrock, both of which assume double "
            "precision; rerun with --arith double");
    } else {
        const std::vector<std::size_t> caches = detail::fluid_cacheqn_nodes(sn);
        if (caches.empty())
            throw InputError(
                "solver_fld_cacheqn_analyzer: the model has no Cache node; use solver_fluid");
        detail::fluid_cacheqn_gate(sn, caches, false);

        // The network solver, the reference's `netsolve`. `init_sol` is left
        // empty so each pass restarts from the default initial state, which is
        // what `solver_fluid_initsol` returns for the default `sn.state`; the
        // drift forgets it anyway at the fixed point.
        //
        // OPT.METHOD SELECTS THE QUEUEING LAYER, not the cache one. "rmf" keeps
        // the first-order matrix method here, which is the historical
        // behaviour; "minnormal" puts the moment closure in its place, so a
        // cache model reaches the same E[min(X,c)] treatment as any other model
        // and returns a covariance. The cache layer is the refined mean field
        // either way: it has no first-order alternative.
        const bool use_moments = (opt.method == "minnormal");
        FluidSolution last;
        bool solved = false;
        std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [&last, &solved, &opt, use_moments](const qn::NetworkStruct<T>& snit) -> mva::MvaSolution<T> {
            FluidOptions fo = opt;
            fo.method = use_moments ? "minnormal" : "matrix";
            fo.init_sol.clear();
            // fluid_dispatch is the FIRST-ORDER dispatcher; the closure lives in
            // solver_fluid_moments, which fluid_runner reaches directly
            last = use_moments ? solver_fluid_moments(snit, fo) : detail::fluid_dispatch(snit, fo);
            solved = true;
            mva::MvaSolution<T> ms;
            ms.Q = last.QN;
            ms.U = last.UN;
            ms.R = last.RN;
            ms.Tp = last.TN;
            ms.C = last.CN;
            ms.X = last.XN;  // TN at the reference station, as netsolve computes XN
            ms.method = last.method;
            ms.iter = static_cast<int>(last.iters);
            return ms;
        };

        std::function<std::vector<T>(const Matrix<T>&, const std::vector<int>&,
                                     const std::vector<Matrix<T> >&, const qn::CacheParam<T>&)>
            missfun = [](const Matrix<T>& gamma, const std::vector<int>& m,
                         const std::vector<Matrix<T> >& lambda_cache,
                         const qn::CacheParam<T>& ch) -> std::vector<T> {
            (void)gamma;  // every drift here is built from the request rates alone
            const std::size_t u = lambda_cache.size();
            std::vector<T> missrate(u, num_traits<T>::from_int(0));
            if (u == 0) return missrate;
            const std::vector<T> mu =
                detail::fluid_cacheqn_miss(m, detail::fluid_cacheqn_lambda_slice(lambda_cache), ch);
            for (std::size_t v = 0; v < mu.size() && v < u; ++v) missrate[v] = mu[v];
            return missrate;
        };

        // `da_cacheqn` reads the fixed-point tolerance off MvaOptions::tol,
        // whereas MATLAB's da_fpi reads options.iter_tol; mapping iter_tol onto
        // it is what keeps the two stopping on the same increment.
        mva::MvaOptions mopt;
        mopt.tol = opt.iter_tol;
        mopt.iter_max = static_cast<int>(opt.iter_max);

        const da::CacheqnResult<T> r = da::da_cacheqn<T>(sn, false, mopt, netfun, missfun);
        if (!solved)
            throw NumericError(
                "solver_fld_cacheqn_analyzer: the decomposition ran no sweep, so there is no fluid "
                "solution; iter_max must be at least one");

        FluidCacheqnSolution<T> out;
        out.sol = last;
        out.sol.method = use_moments ? "minnormal" : "rmf";
        out.sol.iters = static_cast<std::size_t>(r.iter);
        out.hitprob = r.hitprob;
        out.missprob = r.missprob;
        out.iter = r.iter;

        // REFERENCE DEFECT, reproduced deliberately. The analyzer overwrites the
        // system response time with njobs(k)/XN(k) instead of the sum of the
        // queue lengths over the throughput. For a closed class that is Little's
        // law; for an OPEN class njobs is infinite, so the reference reports an
        // infinite response time, and substituting sum(Q)/X here would report a
        // number this analyzer does not produce.
        const std::size_t M = sn.nstations, K = sn.nclasses;
        out.sol.CN.assign(K, 0.0);
        for (std::size_t k = 0; k < K; ++k) {
            const std::size_t rs = sn.classes[k].refstat;
            if (rs >= 1 && rs <= M && out.sol.XN[k] > 0.0)
                out.sol.CN[k] = sn.classes[k].population / out.sol.XN[k];
        }

        out.refreshed = sn;
        out.refreshed.refresh_cacheqn_actual_visits(r.hitprob, r.missprob);
        return out;
    }
}

/**
 * Port of `solver_fld_cacheqn_tran.m`: the transient counterpart of the analyzer
 * above.
 *
 * The steady solver drives the cache drift to its fixed point; this integrates
 * the SAME drift over [t0,t1] from a given occupancy. The decomposition still
 * runs first and in full, because the drift is parameterised by the converged
 * per-class arrival rates -- a transient computed at the offered rates would be
 * the trajectory of a cache the network does not feed.
 *
 * `x0cell[c]`, when present and non-empty, seeds cache c with a flat DDPP state
 * of length nitems*(h+1); otherwise the reference default is used.
 *
 * WHY THE TIME GRID IS PER CACHE. The reference builds one `tcache` from the
 * FIRST cache and writes every other cache's trajectory into a slice of that
 * width. The grids come from an adaptive integrator on different drifts, so they
 * agree in length only by accident; keeping each cache's own grid is the same
 * information without the coincidence.
 */
template <class T>
std::vector<FluidCacheqnTranCache<T> > solver_fld_cacheqn_tran(
    const qn::NetworkStruct<T>& sn, const FluidOptions& opt, double t0, double t1,
    const std::vector<std::vector<T> >& x0cell = std::vector<std::vector<T> >()) {
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        (void)t0;
        (void)t1;
        (void)x0cell;
        throw UnsupportedError(
            "solver_fld_cacheqn_tran: the cache drift is integrated by an adaptive Rosenbrock and "
            "the network by LSODA, both double precision; rerun with --arith double");
    } else {
        const std::vector<std::size_t> caches = detail::fluid_cacheqn_nodes(sn);
        if (caches.empty())
            throw InputError("solver_fld_cacheqn_tran: the model has no Cache node");
        detail::fluid_cacheqn_gate(sn, caches, true);
        if (!std::isfinite(t0)) t0 = 0.0;  // options.timespan(1) = -Inf, as the reference
        if (!(t1 > t0))
            throw InputError("solver_fld_cacheqn_tran: the timespan must have positive width");

        // Converge the cache arrival rates, and keep the isolated-cache inputs
        // the last sweep built -- gamma, the capacities and the per-item rates
        // ARE the drift's parameters.
        FluidSolution last;
        std::function<mva::MvaSolution<T>(const qn::NetworkStruct<T>&)> netfun =
            [&last, &opt](const qn::NetworkStruct<T>& snit) -> mva::MvaSolution<T> {
            FluidOptions fo = opt;
            fo.method = "matrix";
            fo.init_sol.clear();
            last = detail::fluid_dispatch(snit, fo);
            mva::MvaSolution<T> ms;
            ms.Q = last.QN;
            ms.U = last.UN;
            ms.R = last.RN;
            ms.Tp = last.TN;
            ms.C = last.CN;
            ms.X = last.XN;
            return ms;
        };
        std::function<std::vector<T>(const Matrix<T>&, const std::vector<int>&,
                                     const std::vector<Matrix<T> >&, const qn::CacheParam<T>&)>
            missfun = [](const Matrix<T>& gamma, const std::vector<int>& m,
                         const std::vector<Matrix<T> >& lambda_cache,
                         const qn::CacheParam<T>& ch) -> std::vector<T> {
            (void)gamma;
            const std::size_t u = lambda_cache.size();
            std::vector<T> missrate(u, num_traits<T>::from_int(0));
            if (u == 0) return missrate;
            const std::vector<T> mu =
                detail::fluid_cacheqn_miss(m, detail::fluid_cacheqn_lambda_slice(lambda_cache), ch);
            for (std::size_t v = 0; v < mu.size() && v < u; ++v) missrate[v] = mu[v];
            return missrate;
        };
        mva::MvaOptions mopt;
        mopt.tol = opt.iter_tol;
        mopt.iter_max = static_cast<int>(opt.iter_max);
        const da::CacheqnResult<T> r = da::da_cacheqn<T>(sn, false, mopt, netfun, missfun);

        const std::size_t K = sn.nclasses;
        std::vector<FluidCacheqnTranCache<T> > out(caches.size());
        for (std::size_t ci = 0; ci < caches.size(); ++ci) {
            const std::vector<int>& m = r.info.itemcap[ci];
            const std::vector<Matrix<T> >& lam_c = r.info.lambda_cache[ci];
            const Matrix<T> lam = detail::fluid_cacheqn_lambda_slice(lam_c);
            const std::size_t nitems = lam.cols();
            const qn::CacheParam<T>& ch = sn.nodeparam.at(caches[ci] + 1);
            const bool positional = ch.replacestrat == lang::ReplacementStrategy::FIFO ||
                                    ch.replacestrat == lang::ReplacementStrategy::SFIFO;
            // The two families do not share a state space: RANDOM(m) tracks
            // nitems*(h+1) per-list occupancies, the position-resolved policies
            // nitems*sum(m) per-SLOT ones. A seed is validated against the one
            // its own policy uses, and an absent seed is left to the routine's
            // own default (popularity-ordered, or cold on a declared graph)
            // rather than replaced by the RANDOM(m) layout.
            std::vector<T> x0;
            if (ci < x0cell.size() && !x0cell[ci].empty())
                x0 = x0cell[ci];
            else if (!positional)
                x0 = detail::fluid_cacheqn_default_x0<T>(m, nitems);
            std::size_t want = nitems * (m.size() + 1);
            if (positional) {
                want = 0;
                for (std::size_t l = 0; l < m.size(); ++l)
                    want += static_cast<std::size_t>(m[l]);
                want *= nitems;
            }
            if (!x0.empty() && x0.size() != want)
                throw InputError(
                    "solver_fld_cacheqn_tran: the seed occupancy of cache " + std::to_string(ci + 1) +
                    " is not " + std::to_string(want) + " long");

            cache::CacheMissRmfResult<T> tr;
            switch (ch.replacestrat) {
                case lang::ReplacementStrategy::FIFO:
                    tr = cache::cache_miss_fifo_rmf_transient(std::vector<T>(), m, lam, T(t0),
                                                              T(t1), x0, ch.accost);
                    break;
                case lang::ReplacementStrategy::SFIFO:
                    tr = cache::cache_miss_sfifo_rmf_transient(std::vector<T>(), m, lam, T(t0),
                                                               T(t1), x0, ch.accost);
                    break;
                default:
                    tr = cache::cache_miss_rmf_transient(m, lam, T(t0), T(t1), x0);
                    break;
            }

            FluidCacheqnTranCache<T>& oc = out[ci];
            oc.node = caches[ci];
            oc.t = tr.tout;
            oc.xocc = tr.xtraj;
            const std::size_t nt = tr.tout.size();
            oc.hitprob_t = Matrix<T>(K, nt, num_traits<T>::from_int(0));
            oc.missprob_t = Matrix<T>(K, nt, num_traits<T>::from_int(0));
            oc.arate.assign(K, num_traits<T>::from_int(0));
            for (std::size_t v = 0; v < lam.rows() && v < K; ++v) {
                double rowrate = 0.0;
                for (std::size_t i = 0; i < nitems; ++i) rowrate += lam(v, i);
                oc.arate[v] = rowrate;
                // A class that does not read this cache has neither probability:
                // the reference leaves both rows at zero rather than at 1 and 0.
                if (!(rowrate > 0.0)) continue;
                for (std::size_t j = 0; j < nt; ++j) {
                    double mp = tr.MU_t(v, j) / rowrate;
                    if (mp < 0.0) mp = 0.0;
                    if (mp > 1.0) mp = 1.0;
                    oc.missprob_t(v, j) = mp;
                    oc.hitprob_t(v, j) = 1.0 - mp;
                }
            }
        }
        return out;
    }
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_CACHEQN_H
