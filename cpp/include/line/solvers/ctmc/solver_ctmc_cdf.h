/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `@@SolverCTMC/getCdfRespT.m` and `@@SolverCTMC/getCdfSysRespT.m`: the
 * exact distribution of the response time, not just its mean.
 *
 * THE CONSTRUCTION. Tag one job of a chain (`tag_chain`), build the generator of
 * the tagged model WITH its event filtration, and then read the passage of the
 * tagged job off the filtration as a MAP:
 *
 *   A = map_normalize({Q - A1, A1})   A1 = the tagged job ARRIVING at station i
 *   D = map_normalize({Q - D1, D1})   D1 = the tagged job DEPARTING station i
 *
 * `map_pie(A)` is the state of the WHOLE NETWORK as the tagged job sees it on
 * arrival -- the arrival theorem made exact, with no product form assumed -- and
 * from that state the passage ends at the first D1 event. So D0 = D.D0 is the
 * sub-generator of "the tagged job is still in station i", and
 *
 *   F(t) = 1 - pie_arv exp(D0 t) e
 *
 * is the response-time CDF. `getCdfSysRespT` is the same object with the split
 * taken at the tagged job's ARRIVAL at its own reference station, so the passage
 * is one full cycle of the network rather than one visit to one station.
 *
 * WHY THE FILTRATION IS INDISPENSABLE. Q has already summed every
 * synchronization's contribution into one entry, and A1 is one synchronization's
 * share of it; no post-processing of Q can separate them. That is what
 * `CtmcOptions::keep_filtration` is for, and it is forced on below.
 *
 * WHAT IS DELIBERATELY NOT REPRODUCED. The reference recomputes `expm(D0*t)`
 * from scratch at each of the 100001 grid points (or calls `expmv` where the
 * MATLAB release has it). On a uniform grid that is pure waste: exp(D0 k dt) is
 * exp(D0 dt)^k exactly, so ONE matrix exponential and one vector-matrix product
 * per point give the same curve. The recurrence is stable because exp(D0 dt) is
 * substochastic -- its powers contract -- so the rounding of a step is damped by
 * every step after it rather than amplified.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_CDF_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_CDF_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state_events.h"
#include "line/lang/qn/tag_chain.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/**
 * A CDF as the reference returns it: the value at each point of the time grid.
 *
 * `t` and `F` have the same length, and the curve is TRUNCATED at the first
 * point where F is within tolerance of 1 -- carrying the flat tail out to the
 * end of the grid would multiply the size of the result by a factor that
 * depends only on how far the grid overshoots.
 */
template <class T>
struct CdfCurve {
    std::vector<T> t;
    std::vector<T> F;

    bool empty() const { return t.empty(); }
};

namespace cdf_detail {

/** How finely the reference samples each grid; both are `length(0:T/n:T)`. */
constexpr std::size_t respt_intervals = 100000;
constexpr std::size_t sysrespt_intervals = 10000;

/** The tagged model, its synchronization list, and its solved chain. */
template <class T>
struct TaggedCtmc {
    qn::TaggedChain<T> tag;
    std::vector<qn::Sync<T>> sync;
    CtmcSolution<T> sol;
};

/**
 * Port of `s = inchain(N(inchain)>0)`: which class of the chain to tag.
 *
 * The reference indexes with a LOGICAL mask, so `s` is a vector whenever two
 * classes of the chain both hold jobs, and `getClassByIndex(s)` then errors out.
 * The first such class is taken here instead. It is not an approximation: the
 * tagged twin exists for every class of the chain and the tagged job switches
 * between the twins exactly as an untagged one switches between the originals,
 * so which class it starts in changes only its initial position and not the
 * stationary passage the CDF describes.
 */
template <class T>
std::size_t class_to_tag(const NetworkStruct<T>& sn, std::size_t chain) {
    const std::vector<std::size_t>& ic = sn.inchain[chain - 1];
    for (std::size_t a = 0; a < ic.size(); ++a)
        if (sn.classes[ic[a] - 1].population > 0.0) return ic[a];
    throw InputError("SolverCTMC: chain " + std::to_string(chain) +
                     " holds no jobs, so it has no job to tag");
}

/** Tag one chain, then build and solve the generator of the tagged model. */
template <class T>
TaggedCtmc<T> tagged_ctmc(const NetworkStruct<T>& sn, std::size_t chain, const CtmcOptions& opt) {
    TaggedCtmc<T> out;
    out.tag = qn::tag_chain(sn, chain, class_to_tag(sn, chain));

    CtmcOptions o = opt;
    o.keep_filtration = true;  // the whole method is a split of Q on one event

    // The analyzer is used rather than `solver_ctmc` directly because the tagged
    // model is REDUCIBLE far more often than the original: the twins carry a
    // single job between them, so every lattice state that puts two tagged jobs
    // in the network is enumerated and unreachable. The analyzer keeps the
    // component of the initial state and restricts the filtration with Q, which
    // a raw generator would not do, and map_pie on a reducible generator has no
    // unique answer at all.
    out.sol = solver_ctmc_analyzer(out.tag.V, o);

    // `refresh_sync` is deterministic, so recomputing it here reproduces exactly
    // the list the analyzer filtered on; `filt[a]` belongs to `sync[a]`.
    out.sync = qn::refresh_sync(out.tag.V);
    if (out.sol.chain.filt.size() != out.sync.size())
        throw NumericError(
            "SolverCTMC: the event filtration does not match the synchronization list; the "
            "generator was built without keep_filtration");
    return out;
}

/**
 * The tagged classes whose passage through the reference station counts as a
 * completion, which is the reference's `taggedModel.classes{r}.completes` gate.
 *
 * A CHAIN'S CLASSES ARE ALL SUMMED, not just the one that was tagged. The single
 * tagged job switches class as it circulates, so "the tagged job arrives at
 * station i" is the union of the per-class arrival events over the whole tagged
 * block. The reference reaches the same set by a loop variable that shadows its
 * own outer index, which is why its per-class results within a chain are
 * identical.
 */
template <class T>
std::vector<std::size_t> completing_tagged(const qn::TaggedChain<T>& tag) {
    std::vector<std::size_t> out;
    for (std::size_t a = 0; a < tag.tagged.size(); ++a)
        if (tag.V.classes[tag.tagged[a] - 1].completes) out.push_back(tag.tagged[a]);
    if (out.empty())
        throw UnsupportedError(
            "SolverCTMC: no class of the tagged chain has JobClass::completes set, so the response "
            "time has no completion event to end at. `completes` now defaults to true here as it "
            "does in the other three codebases, so reaching this means a caller cleared it; set "
            "classes of the chain before asking for a response-time CDF");
    return out;
}

/**
 * Sum the filtration matrices of every synchronization matching one of the
 * (class, node) selectors on the chosen half.
 *
 * THE REFERENCE COMPARES A STATION INDEX AGAINST A NODE INDEX here: `ist` runs
 * over stations and `tsn.refstat(r)` is a station, while `ev{v}.passive{1}.node`
 * is a node. The two coincide only on a model whose stations are its first
 * nodes, and silently select the wrong node otherwise, so the selectors are
 * built from `node_of_station` on the way in.
 */
template <class T>
Matrix<T> filtration_sum(const CtmcSolution<T>& sol, const std::vector<qn::Sync<T>>& sync,
                         bool on_passive, EventType ev,
                         const std::vector<std::pair<std::size_t, std::size_t>>& sel) {
    const std::size_t n = sol.chain.Q.rows();
    Matrix<T> S(n, n, num_traits<T>::from_int(0));
    for (std::size_t v = 0; v < sync.size(); ++v) {
        const qn::SyncEvent<T>& half = on_passive ? sync[v].passive : sync[v].active;
        if (half.event != ev) continue;
        bool hit = false;
        for (std::size_t k = 0; k < sel.size(); ++k)
            if (half.cls == sel[k].first && half.node == sel[k].second) hit = true;
        if (!hit) continue;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) S(i, j) += sol.chain.filt[v](i, j);
    }
    return S;
}

/** Whether any rate at all was selected; an empty split has no MAP. */
template <class T>
bool any_rate(const Matrix<T>& S) {
    for (std::size_t i = 0; i < S.rows(); ++i)
        for (std::size_t j = 0; j < S.cols(); ++j)
            if (num_traits<T>::to_double(S(i, j)) != 0) return true;
    return false;
}

/**
 * Port of the reference's grid: `T = 100/min(nonzero rate)`, in `intervals`
 * steps, which is "long enough for 100 events at the slowest rate in the model".
 *
 * The diagonal is INCLUDED in the minimum, as `Q(Q~=0)` includes it: it is minus
 * the total exit rate of a state, and a state that leaves slowly is exactly the
 * one that sets the horizon.
 */
template <class T>
T grid_step(const Matrix<T>& Q, std::size_t intervals) {
    double lo = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < Q.rows(); ++i)
        for (std::size_t j = 0; j < Q.cols(); ++j) {
            const double a = std::fabs(num_traits<T>::to_double(Q(i, j)));
            if (a > GlobalConstants::FineTol && a < lo) lo = a;
        }
    if (!std::isfinite(lo))
        throw NumericError("SolverCTMC: the generator of the tagged model carries no rate above "
                           "the tolerance, so no time horizon can be chosen");
    const double horizon = 100.0 / lo;
    return num_traits<T>::from_double(horizon / static_cast<double>(intervals));
}

/**
 * F(t) = 1 - pie exp(D0 t) e on the uniform grid t = 0, dt, 2dt, ...
 *
 * One exponential, then a vector-matrix product per point: see the file header
 * for why this is the same curve the reference computes point by point.
 */
template <class T>
CdfCurve<T> absorption_cdf(const Matrix<T>& D0, const std::vector<T>& pie, const T& dt,
                           std::size_t intervals, double tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "the response-time CDF needs a matrix exponential and is not available in exact "
                  "rational arithmetic");
    const T one = num_traits<T>::from_int(1);
    const Matrix<T> E = expm(D0, dt);
    std::vector<T> v = pie;
    CdfCurve<T> out;
    for (std::size_t k = 0; k <= intervals; ++k) {
        T mass = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < v.size(); ++i) mass += v[i];
        const T F = T(one - mass);
        out.t.push_back(T(dt * num_traits<T>::from_int(static_cast<long>(k))));
        out.F.push_back(F);
        // The reference keeps the point that crossed the threshold and drops
        // everything after it, so the curve ends ON the crossing.
        if (num_traits<T>::to_double(F) > 1.0 - tol) break;
        v = vecmul(v, E);
    }
    return out;
}

}  // namespace cdf_detail

/**
 * Port of `@@SolverCTMC/getCdfRespT.m`: the per-(station, class) response-time
 * CDF, indexed `[ist-1][r-1]`.
 *
 * Entries are EMPTY where the curve does not exist: a station the chain never
 * visits has no arrival event to condition on, and a class outside every chain
 * that holds jobs is never tagged. Every class of a chain gets the SAME curve,
 * which is what the reference produces -- the tagged job's passage through a
 * station is one quantity per chain, since the job carries its class with it.
 */
template <class T>
std::vector<std::vector<CdfCurve<T>>> solver_ctmc_cdf_respt(const NetworkStruct<T>& sn,
                                                            const CtmcOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "getCdfRespT needs a matrix exponential and is not available in exact rational "
                  "arithmetic");
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const std::vector<double> N = sn.njobs();
    for (std::size_t k = 0; k < K; ++k)
        if (!std::isfinite(N[k]))
            throw UnsupportedError(
                "SolverCTMC: getCdfRespT is presently supported only for closed models; class '" +
                sn.classes[k].name + "' is open");

    std::vector<std::vector<CdfCurve<T>>> RD(M, std::vector<CdfCurve<T>>(K));

    for (std::size_t c = 1; c <= sn.nchains; ++c) {
        const cdf_detail::TaggedCtmc<T> tc = cdf_detail::tagged_ctmc(sn, c, opt);
        const std::vector<std::size_t> tagged = cdf_detail::completing_tagged(tc.tag);
        const Matrix<T>& Q = tc.sol.chain.Q;
        const T dt = cdf_detail::grid_step(Q, cdf_detail::respt_intervals);

        for (std::size_t ist = 1; ist <= M; ++ist) {
            const std::size_t ind = tc.tag.V.node_of_station(ist);
            std::vector<std::pair<std::size_t, std::size_t>> sel;
            for (std::size_t a = 0; a < tagged.size(); ++a) sel.push_back(std::make_pair(tagged[a], ind));

            const Matrix<T> A1 =
                cdf_detail::filtration_sum(tc.sol, tc.sync, true, EventType::ARV, sel);
            // No arrival of the tagged job here: the chain does not visit this
            // station, and conditioning on an event of rate zero is not a
            // degenerate CDF but no CDF at all.
            if (!cdf_detail::any_rate(A1)) continue;
            const Matrix<T> D1 =
                cdf_detail::filtration_sum(tc.sol, tc.sync, false, EventType::DEP, sel);
            if (!cdf_detail::any_rate(D1)) continue;

            mam::Map<T> A;
            A.D0 = Q;
            A.D1 = A1;
            for (std::size_t i = 0; i < A.D0.rows(); ++i)
                for (std::size_t j = 0; j < A.D0.cols(); ++j) A.D0(i, j) -= A1(i, j);
            const std::vector<T> pie_arv = mam::map_pie(mam::map_normalize(A));

            mam::Map<T> D;
            D.D0 = Q;
            D.D1 = D1;
            for (std::size_t i = 0; i < D.D0.rows(); ++i)
                for (std::size_t j = 0; j < D.D0.cols(); ++j) D.D0(i, j) -= D1(i, j);
            const mam::Map<T> Dn = mam::map_normalize(D);

            const CdfCurve<T> curve = cdf_detail::absorption_cdf(
                Dn.D0, pie_arv, dt, cdf_detail::respt_intervals, GlobalConstants::CoarseTol);
            for (std::size_t a = 0; a < sn.inchain[c - 1].size(); ++a)
                RD[ist - 1][sn.inchain[c - 1][a] - 1] = curve;
        }
    }
    return RD;
}

/**
 * Port of `@@SolverCTMC/getCdfSysRespT.m`: the per-chain SYSTEM response-time
 * CDF, indexed by chain.
 *
 * The split is taken at the tagged job's arrival at its OWN reference station,
 * for both halves of the MAP: the state seen on arrival there, and the
 * absorption at the next arrival there. The passage is therefore one full cycle
 * of the network -- MATLAB's system response time -- and not a single visit.
 *
 * The reference has no closed-model guard here, unlike getCdfRespT; `tag_chain`
 * supplies one, since tagging moves a job out of a finite population.
 */
template <class T>
std::vector<CdfCurve<T>> solver_ctmc_cdf_sys_respt(const NetworkStruct<T>& sn,
                                                   const CtmcOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "getCdfSysRespT needs a matrix exponential and is not available in exact "
                  "rational arithmetic");
    std::vector<CdfCurve<T>> RD(sn.nchains);

    for (std::size_t c = 1; c <= sn.nchains; ++c) {
        const cdf_detail::TaggedCtmc<T> tc = cdf_detail::tagged_ctmc(sn, c, opt);
        const std::vector<std::size_t> tagged = cdf_detail::completing_tagged(tc.tag);
        const Matrix<T>& Q = tc.sol.chain.Q;

        // The reference station is read PER CLASS: the twins of one chain need
        // not share one, and it is each class's own reference station whose
        // passage counts as that class's completion.
        std::vector<std::pair<std::size_t, std::size_t>> sel;
        for (std::size_t a = 0; a < tagged.size(); ++a) {
            const std::size_t refstat = tc.tag.V.classes[tagged[a] - 1].refstat;
            sel.push_back(std::make_pair(tagged[a], tc.tag.V.node_of_station(refstat)));
        }

        const Matrix<T> D1 = cdf_detail::filtration_sum(tc.sol, tc.sync, true, EventType::ARV, sel);
        if (!cdf_detail::any_rate(D1))
            throw NumericError(
                "SolverCTMC: the tagged job never arrives at its reference station, so the system "
                "response time of chain " + std::to_string(c) + " has no cycle to measure");

        mam::Map<T> D;
        D.D0 = Q;
        D.D1 = D1;
        for (std::size_t i = 0; i < D.D0.rows(); ++i)
            for (std::size_t j = 0; j < D.D0.cols(); ++j) D.D0(i, j) -= D1(i, j);
        const mam::Map<T> Dn = mam::map_normalize(D);
        const std::vector<T> pie_arv = mam::map_pie(Dn);

        const T dt = cdf_detail::grid_step(Q, cdf_detail::sysrespt_intervals);
        RD[c - 1] = cdf_detail::absorption_cdf(Dn.D0, pie_arv, dt, cdf_detail::sysrespt_intervals,
                                               GlobalConstants::FineTol);
    }
    return RD;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_CDF_H
