/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_LQN_BOXBOUNDS_H
#define LINE_API_LQN_LQN_BOXBOUNDS_H

/**
 * Majumdar-Woodside robust box bounds on the throughput of a layered network.
 *
 * Port of `matlab/src/api/lqn/lqn_boxbounds.m`. The layered model is collapsed
 * onto its PROCESSOR-CONTENTION model -- the stations are the hosts and the
 * classes are the reference-task call chains -- and `pfqn_mwrbb` is evaluated
 * on it. The per-chain demand at a processor is the total host demand executed
 * there during one cycle of the reference task, obtained by walking the
 * entry/activity/call graph and scaling by the mean number of synchronous
 * calls.
 *
 * WHAT THIS IS AND IS NOT. It generalizes the classical LQN Type-1 throughput
 * bound X <= mult/(Z + D_total) -- the no-contention bound `lqns -b` reports --
 * by adding the processor-utilization upper bound and the Majumdar-Woodside
 * lower bound. It does NOT model a software bottleneck (a task with finitely
 * many threads is not a station here) and it assumes activities execute
 * sequentially, so an OR-branch probability and a loop count are ignored rather
 * than averaged. A bound that ignores a resource is still a valid bound on the
 * model it does describe, which is why this is worth having; it is not a
 * substitute for solving the layers.
 *
 * NO ITERATION HAPPENS HERE. This is what `method = 'mwba.upper'` /
 * `'mwba.lower'` reports INSTEAD of a fixed point, so a caller asking for it
 * pays one graph walk and one bound evaluation, and gets NaN for every metric
 * the bound does not define (queue lengths, response times, residence times).
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mwrbb.h"
#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lqn {

/** What lqn_boxbounds returns, in the layout of the reference's `out` struct. */
template <class T>
struct LqnBoxBounds {
    std::vector<std::size_t> refidx;  ///< (R) absolute indices of the reference tasks
    std::vector<T> Xlo, Xup;          ///< (R) per-chain throughput bounds
    /// (nidx+1) bounds propagated to every element; `defined` marks the entries
    /// the reference leaves as NaN, i.e. an element no reference chain visits.
    std::vector<T> TN_lo, TN_up, UN_lo, UN_up;
    std::vector<bool> defined_T, defined_U;
    Matrix<T> D;  ///< (nhosts x R) per-chain demand at each processor
};

namespace detail {

/** SchedStrategy -> the Majumdar-Woodside discipline code, discCode of the reference. */
inline pfqn::MwrbbSched lqn_disc_code(lang::SchedStrategy s) {
    using lang::SchedStrategy;
    switch (s) {
        case SchedStrategy::FCFS:
            return pfqn::MwrbbSched::Fifo;
        case SchedStrategy::PS:
        case SchedStrategy::DPS:
        case SchedStrategy::GPS:
        case SchedStrategy::PSPRIO:
        case SchedStrategy::DPSPRIO:
        case SchedStrategy::GPSPRIO:
            return pfqn::MwrbbSched::Ps;
        case SchedStrategy::HOL:
            return pfqn::MwrbbSched::PrioNonPreemptive;
        case SchedStrategy::FCFSPRPRIO:
        case SchedStrategy::LCFSPRPRIO:
            return pfqn::MwrbbSched::PrioPreemptive;
        default:
            // ABA full contention: discipline-independent, hence the safe answer
            // for a discipline the bound has no shape for.
            return pfqn::MwrbbSched::Aba;
    }
}

template <class T>
void lqn_box_visit_entry(const LqnStruct<T>& lqn, std::size_t eidx, const T& mult, std::size_t nH,
                         std::vector<T>& d, std::vector<T>& vis);

/** One activity: charge its host demand, then follow its synchronous calls. */
template <class T>
void lqn_box_visit_activity(const LqnStruct<T>& lqn, std::size_t aidx, const T& mult,
                            std::size_t nH, std::vector<T>& d, std::vector<T>& vis) {
    vis[aidx] = T(vis[aidx] + mult);
    const std::size_t tidx = lqn.parent[aidx];
    const std::size_t hidx = lqn.parent[tidx];  // hosts are their own shift, hshift = 0
    const T hd = lqn.hostdem[aidx].disabled ? num_traits<T>::from_int(0) : lqn.hostdem[aidx].mean;
    if (hidx >= 1 && hidx <= nH) d[hidx] = T(d[hidx] + mult * hd);
    for (std::size_t cidx : lqn.callsof[aidx])
        if (lqn.calltype[cidx] == lang::CallType::SYNC)
            lqn_box_visit_entry(lqn, lqn.callpair_dst[cidx],
                                T(mult * lqn.callproc_mean[cidx]), nH, d, vis);
}

/** One entry: itself, its task, and every activity of its own task. */
template <class T>
void lqn_box_visit_entry(const LqnStruct<T>& lqn, std::size_t eidx, const T& mult, std::size_t nH,
                         std::vector<T>& d, std::vector<T>& vis) {
    vis[eidx] = T(vis[eidx] + mult);
    const std::size_t tidx = lqn.parent[eidx];
    if (tidx >= 1 && tidx < vis.size()) vis[tidx] = T(vis[tidx] + mult);
    for (std::size_t aidx : lqn.actsof[eidx])
        if (lqn.parent[aidx] == lqn.parent[eidx])
            lqn_box_visit_activity(lqn, aidx, mult, nH, d, vis);
}

}  // namespace detail

/**
 * Evaluate the box bounds of `lqn`.
 *
 * A reference task of infinite multiplicity is treated as one customer, as the
 * reference does: the bound is over a closed chain and an infinite population
 * has no lower bound to report.
 */
template <class T>
LqnBoxBounds<T> lqn_boxbounds(const LqnStruct<T>& lqn) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t nidx = lqn.nidx, nH = lqn.nhosts;

    LqnBoxBounds<T> out;
    for (std::size_t t = 1; t <= lqn.ntasks; ++t) {
        const std::size_t tidx = lqn.tshift + t;
        if (lqn.isref[tidx]) out.refidx.push_back(tidx);
    }
    const std::size_t R = out.refidx.size();
    if (R == 0)
        throw InputError("lqn_boxbounds: the model declares no reference task, so it has no chain");

    Matrix<T> D(nH + 1, R, zero);   // row 0 unused, matching the 1-based host index
    Matrix<T> Vis(nidx + 1, R, zero);
    std::vector<T> Nref(R, one), Zref(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        const std::size_t tidx = out.refidx[r];
        const double mult = lqn.mult[tidx];
        Nref[r] = std::isfinite(mult) ? num_traits<T>::from_double(mult) : one;
        Zref[r] = lqn.think[tidx].disabled ? zero : lqn.think[tidx].mean;
        std::vector<T> d(nH + 1, zero), vis(nidx + 1, zero);
        for (std::size_t eidx : lqn.entriesof[tidx])
            detail::lqn_box_visit_entry(lqn, eidx, one, nH, d, vis);
        for (std::size_t h = 0; h <= nH; ++h) D(h, r) = d[h];
        for (std::size_t i = 0; i <= nidx; ++i) Vis(i, r) = vis[i];
    }

    // The bound's own station table: one row per host, visits as an indicator,
    // service the whole per-cycle demand. That is the reference's V = D > 0,
    // S = D, i.e. the demand is charged once per cycle rather than per visit.
    Matrix<T> V(nH, R, zero), S(nH, R, zero);
    for (std::size_t h = 1; h <= nH; ++h)
        for (std::size_t r = 0; r < R; ++r) {
            S(h - 1, r) = D(h, r);
            V(h - 1, r) = D(h, r) > zero ? one : zero;
        }
    std::vector<pfqn::MwrbbSched> sched(nH);
    for (std::size_t h = 1; h <= nH; ++h) sched[h - 1] = detail::lqn_disc_code(lqn.sched[h]);
    const std::vector<int> prio(R, 0);  // equal class priority, as the reference sets

    const pfqn::MwrbbBounds<T> b = pfqn::pfqn_mwrbb(V, S, Nref, Zref, sched, prio);
    out.Xlo = b.Xlo;
    out.Xup = b.Xup;

    out.TN_lo.assign(nidx + 1, zero);
    out.TN_up.assign(nidx + 1, zero);
    out.UN_lo.assign(nidx + 1, zero);
    out.UN_up.assign(nidx + 1, zero);
    out.defined_T.assign(nidx + 1, false);
    out.defined_U.assign(nidx + 1, false);
    for (std::size_t i = 1; i <= nidx; ++i) {
        bool visited = false;
        for (std::size_t r = 0; r < R; ++r)
            if (Vis(i, r) > zero) visited = true;
        if (!visited) continue;
        out.defined_T[i] = true;
        for (std::size_t r = 0; r < R; ++r) {
            out.TN_lo[i] = T(out.TN_lo[i] + b.Xlo[r] * Vis(i, r));
            out.TN_up[i] = T(out.TN_up[i] + b.Xup[r] * Vis(i, r));
        }
    }
    for (std::size_t h = 1; h <= nH; ++h) {
        out.defined_U[h] = true;
        for (std::size_t r = 0; r < R; ++r) {
            out.UN_lo[h] = T(out.UN_lo[h] + b.Xlo[r] * D(h, r));
            out.UN_up[h] = T(out.UN_up[h] + b.Xup[r] * D(h, r));
        }
    }
    out.D = Matrix<T>(nH, R, zero);
    for (std::size_t h = 1; h <= nH; ++h)
        for (std::size_t r = 0; r < R; ++r) out.D(h - 1, r) = D(h, r);
    return out;
}

}  // namespace lqn
}  // namespace line

#endif  // LINE_API_LQN_LQN_BOXBOUNDS_H
