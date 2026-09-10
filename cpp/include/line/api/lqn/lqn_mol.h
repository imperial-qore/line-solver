/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_MOL_H
#define LINE_API_LQN_MOL_H

/**
 * Method of Layers on the SRVN decomposition of a layered queueing network
 * whose entries carry no activity graph.
 *
 * Port of matlab/src/api/lqn/lqn_mol.m. A compact, self-contained
 * reimplementation of the layered fixed point SolverLN runs, restricted to LQNs
 * in which every entry binds exactly one activity and there are no activity
 * precedences. It decomposes the model the way `lqns --srvn-layering` does --
 * one submodel per processor and one per called task -- and sweeps them in the
 * two phases of Rolia-Sevcik's Method of Layers: all software (task) submodels,
 * then all hardware (processor) ones.
 *
 * Every submodel is a closed multiclass queueing network with ONE station and
 * one class per client task, so it is solved by `pfqn_qdamva` rather than by
 * building a Network. The surrogate client delay of SolverLN collapses into the
 * think-time vector Z of that call.
 *
 * WHERE THIS DIFFERS FROM SolverLN's `srvn.cs`: a submodel here carries one
 * class per client TASK with visit-weighted demands, where `srvn.cs` carries
 * one class per activity and encodes the call multiplicities as routing. On an
 * entry-only model the two agree on the structure and differ only in the
 * aggregation, so the throughputs and processor utilizations track closely
 * while entry response times spread more.
 *
 * SCOPE. Entry-only models. Activity graphs (fork/join, OR-branches, loops,
 * second phases, forwarding), asynchronous calls, caches, setup tasks,
 * admission constraints, replication and open arrivals are REFUSED, not
 * approximated, and named when they are.
 *
 * Arithmetic: TRANSCENDENTAL-GATED, inherited from `pfqn_qdamva`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_qdamva.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lqn {

/** Everything `lqn_mol` reports beyond the four measure vectors. */
template <class T>
struct LqnMolInfo {
    std::size_t iter = 0;
    double resid = 0.0;
    std::vector<T> servt;      ///< (nidx+1) entry response time seen by a caller
    std::vector<T> residt;     ///< (nidx+1) entry processor residence
    std::vector<T> callservt;  ///< (ncalls+1) blocking time per call
    std::vector<T> thinkt;     ///< (nidx+1) task surrogate idle time
    std::vector<T> share;      ///< (nidx+1) entry share of its task's invocations
    std::vector<std::size_t> hostLayers;
    std::vector<std::size_t> taskLayers;
};

/**
 * The four (nidx+1) vectors in the column convention SolverLN and LQNS report,
 * so they line up with `LN(model).getAvgTable` cell for cell.
 *
 * | index | QN (QLen) | UN (Util) | RN (RespT) | TN (Tput) |
 * |-------|-----------|-----------|------------|-----------|
 * | host  | NaN       | processor utilization | NaN | NaN |
 * | task  | sum of entry T*S | sum of entry proc util | NaN | cycle rate |
 * | entry | T*S       | processor utilization | response time | throughput |
 * | act   | as its entry | as its entry | as its entry | as its entry |
 */
template <class T>
struct LqnMolResult {
    std::vector<T> QN, UN, RN, TN;
    LqnMolInfo<T> info;
};

/** Tuning of the outer fixed point. */
struct LqnMolOptions {
    std::size_t iter_max = 200;
    double iter_tol = 1e-6;
    double relax_factor = 0.5;
};

namespace mol_detail {

/** GlobalConstants.FineTol, the reference's own "effectively zero". */
inline double fine_tol() { return 1e-8; }

/** CallType.toText, for the refusal messages. */
inline const char* call_kind(CallType c) {
    switch (c) {
        case CallType::SYNC: return "synchronous";
        case CallType::ASYNC: return "asynchronous";
        case CallType::FWD: return "forwarding";
        default: return "none";
    }
}

/**
 * Refuse every feature this decomposition does not represent, NAMING the
 * element, rather than returning a number that quietly ignores it.
 */
template <class T>
void lqn_mol_assert(const LqnStruct<T>& lsn) {
    for (std::size_t eidx = lsn.eshift + 1; eidx <= lsn.eshift + lsn.nentries; ++eidx) {
        const std::size_t n = (eidx < lsn.actsof.size()) ? lsn.actsof[eidx].size() : 0;
        if (n != 1)
            throw UnsupportedError("lqn_mol: entry " + lsn.hashnames[eidx] + " binds " +
                                   std::to_string(n) +
                                   " activities. lqn_mol solves entry-only models; use SolverLN "
                                   "for an activity graph");
    }
    for (std::size_t a = 1; a <= lsn.nacts; ++a) {
        const std::size_t aidx = lsn.ashift + a;
        // A successor INSIDE the activity band is a precedence; an edge to an
        // entry or a task is the ordinary binding every entry-only model has.
        const std::vector<std::size_t> succ = lsn.graph.succ(aidx);
        for (std::size_t k = 0; k < succ.size(); ++k)
            if (succ[k] > lsn.ashift && succ[k] <= lsn.ashift + lsn.nacts)
                throw UnsupportedError("lqn_mol: activity " + lsn.hashnames[aidx] +
                                       " has an activity precedence. lqn_mol solves entry-only "
                                       "models; use SolverLN for an activity graph");
        if (a < lsn.actphase.size() && lsn.actphase[a] != 1)
            throw UnsupportedError("lqn_mol: activity " + lsn.hashnames[aidx] + " is in phase " +
                                   std::to_string(lsn.actphase[a]) +
                                   ". lqn_mol supports phase 1 only");
    }
    for (std::size_t c = 1; c <= lsn.ncalls; ++c)
        if (lsn.calltype[c] != CallType::SYNC)
            throw UnsupportedError("lqn_mol: call " + lsn.callhashnames[c] + " is " +
                                   std::string(mol_detail::call_kind(lsn.calltype[c])) +
                                   ". lqn_mol supports synchronous calls only");
    for (std::size_t idx = 1; idx <= lsn.tshift + lsn.ntasks; ++idx) {
        if (idx < lsn.iscache.size() && lsn.iscache[idx])
            throw UnsupportedError("lqn_mol: " + lsn.hashnames[idx] +
                                   " is a cache task, which lqn_mol does not model");
        if (idx < lsn.hassetup.size() && lsn.hassetup[idx])
            throw UnsupportedError("lqn_mol: " + lsn.hashnames[idx] +
                                   " has a setup time, which lqn_mol does not model");
        if (idx < lsn.repl.size() && lsn.repl[idx] != 1.0)
            throw UnsupportedError("lqn_mol: " + lsn.hashnames[idx] + " is replicated " +
                                   std::to_string(static_cast<long long>(lsn.repl[idx])) +
                                   " times, which lqn_mol does not model");
        if (idx < lsn.lincon_A.size() && lsn.lincon_A[idx].rows() > 0)
            throw UnsupportedError("lqn_mol: " + lsn.hashnames[idx] +
                                   " carries an admission constraint, which lqn_mol does not "
                                   "model");
        const SchedStrategy s = lsn.sched[idx];
        const bool ok = (s == SchedStrategy::PS || s == SchedStrategy::FCFS ||
                         s == SchedStrategy::INF ||
                         (idx > lsn.nhosts && s == SchedStrategy::REF));
        if (!ok)
            throw UnsupportedError("lqn_mol: " + lsn.hashnames[idx] + " is scheduled " +
                                   std::string(sched_to_text(s)) +
                                   ", which lqn_mol does not model");
    }
    if (!lsn.callgroups.empty())
        throw UnsupportedError(
            "lqn_mol: this model uses routed call groups, which lqn_mol does not model");
    for (std::size_t eidx = lsn.eshift + 1; eidx <= lsn.eshift + lsn.nentries; ++eidx)
        if (eidx < lsn.has_arrival.size() && lsn.has_arrival[eidx])
            throw UnsupportedError("lqn_mol: entry " + lsn.hashnames[eidx] +
                                   " has an open arrival. lqn_mol solves closed models only");
}

/**
 * The queue-dependent rate multiplier row of a c-server station over a
 * population of sum(N).
 *
 * `pfqn_lldfun` SKIPS a constant row, so a single server must come back as a
 * row of ones and not as a scalar 1, or the multiserver term is never applied.
 */
template <class T>
Matrix<T> mol_mu(const std::vector<T>& N, double c) {
    double tot = 0.0;
    for (std::size_t k = 0; k < N.size(); ++k) tot += num_traits<T>::to_double(N[k]);
    const std::size_t smax = std::max<std::size_t>(
        2, static_cast<std::size_t>(std::ceil(std::max(0.0, tot))));
    Matrix<T> mu(1, smax, num_traits<T>::from_int(1));
    if (std::isfinite(c) && c > 1.0)
        for (std::size_t n = 1; n <= smax; ++n)
            mu(0, n - 1) = num_traits<T>::from_double(std::min(static_cast<double>(n), c));
    return mu;
}

}  // namespace mol_detail

/**
 * @param lsn     LayeredNetworkStruct of an entry-only, closed, synchronous LQN
 * @param options iteration cap, tolerance and under-relaxation factor
 */
template <class T>
LqnMolResult<T> lqn_mol(const LqnStruct<T>& lsn, const LqnMolOptions& options = LqnMolOptions()) {
    static_assert(num_traits<T>::has_transcendental, "lqn_mol requires transcendental arithmetic");
    mol_detail::lqn_mol_assert(lsn);

    const T zero = num_traits<T>::from_int(0);
    const std::size_t nidx = lsn.nidx, ncalls = lsn.ncalls;
    const std::size_t e0 = lsn.eshift + 1, e1 = lsn.eshift + lsn.nentries;
    const std::size_t t0 = lsn.tshift + 1, t1 = lsn.tshift + lsn.ntasks;
    const double om = options.relax_factor;

    // ---- static per-entry data: bound activity, host demand, owning task ----
    std::vector<std::size_t> actof(nidx + 1, 0), taskof(nidx + 1, 0), hostof(nidx + 1, 0);
    std::vector<T> dem(nidx + 1, zero);
    for (std::size_t tidx = t0; tidx <= t1; ++tidx) hostof[tidx] = lsn.parent[tidx];
    for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
        actof[eidx] = lsn.actsof[eidx][0];
        const T d = lsn.hostdem[actof[eidx]].mean;
        dem[eidx] = std::isnan(num_traits<T>::to_double(d)) ? zero : d;
        taskof[eidx] = lsn.parent[eidx];
    }

    // ---- static per-call data ----------------------------------------------
    std::vector<std::size_t> callsrc(ncalls + 1, 0), calldst(ncalls + 1, 0);
    std::vector<T> cally(ncalls + 1, zero);
    std::vector<std::size_t> entryOfAct(nidx + 1, 0);
    for (std::size_t eidx = e0; eidx <= e1; ++eidx) entryOfAct[actof[eidx]] = eidx;
    std::vector<std::vector<std::size_t>> callsFrom(nidx + 1), callsTo(nidx + 1);
    for (std::size_t c = 1; c <= ncalls; ++c) {
        callsrc[c] = entryOfAct[lsn.callpair_src[c]];
        calldst[c] = lsn.callpair_dst[c];
        const T y = lsn.callproc_mean[c];
        cally[c] = std::isnan(num_traits<T>::to_double(y)) ? zero : y;
        callsFrom[callsrc[c]].push_back(c);
        callsTo[calldst[c]].push_back(c);
    }

    // ---- populations, from maxmult (mult is wrong for INF tasks) ------------
    std::vector<double> npop(nidx + 1, 1.0);
    for (std::size_t idx = 1; idx <= t1; ++idx) {
        double m = (idx < lsn.maxmult.size()) ? lsn.maxmult[idx] : 1.0;
        if (!std::isfinite(m) || m < 1.0) m = 1.0;
        npop[idx] = m;
    }

    // ---- layer sets --------------------------------------------------------
    // One hardware layer per populated host, one software layer per called
    // non-reference task, as buildLayers.m draws them.
    std::vector<std::size_t> hostLayers, taskLayers;
    for (std::size_t hidx = 1; hidx <= lsn.nhosts; ++hidx)
        if (!lsn.tasksof[hidx].empty()) hostLayers.push_back(hidx);
    std::vector<bool> isCalled(nidx + 1, false);
    for (std::size_t c = 1; c <= ncalls; ++c) isCalled[taskof[calldst[c]]] = true;
    for (std::size_t tidx = t0; tidx <= t1; ++tidx)
        if (!lsn.isref[tidx] && isCalled[tidx]) taskLayers.push_back(tidx);

    // ---- fixed-point state -------------------------------------------------
    std::vector<T> residt = dem, servt(nidx + 1, zero), thinkt(nidx + 1, zero),
                   share(nidx + 1, zero), Xtask(nidx + 1, zero), Xentry(nidx + 1, zero),
                   busyth(nidx + 1, zero), zref(nidx + 1, zero);
    std::vector<T> callservt(ncalls + 1, zero);
    for (std::size_t tidx = t0; tidx <= t1; ++tidx) {
        // lqn_ref_thinktime: a reference task's declared think time, and zero
        // for every other task and for a negative or non-finite one.
        T z = zero;
        if (lsn.isref[tidx] && tidx < lsn.think.size() && !lsn.think[tidx].disabled) {
            z = lsn.think[tidx].mean;
            const double zd = num_traits<T>::to_double(z);
            if (!std::isfinite(zd) || zd < 0.0) z = zero;
        }
        zref[tidx] = z;
        thinkt[tidx] = z;
        const std::vector<std::size_t>& es = lsn.entriesof[tidx];
        if (!es.empty())
            for (std::size_t k = 0; k < es.size(); ++k)
                share[es[k]] = num_traits<T>::from_double(1.0 / static_cast<double>(es.size()));
    }
    // Seed servt bottom-up over the call graph so a callee is priced before its
    // caller; a cycle just leaves the residual demand seeded at 0.
    for (std::size_t eidx = e0; eidx <= e1; ++eidx) servt[eidx] = dem[eidx];
    for (std::size_t pass = 0; pass < std::max<std::size_t>(1, lsn.nentries); ++pass)
        for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
            T s = residt[eidx];
            for (std::size_t k = 0; k < callsFrom[eidx].size(); ++k) {
                const std::size_t c = callsFrom[eidx][k];
                s += T(cally[c] * servt[calldst[c]]);
            }
            servt[eidx] = s;
        }
    for (std::size_t c = 1; c <= ncalls; ++c) callservt[c] = servt[calldst[c]];

    // ---- submodel solvers --------------------------------------------------

    /* Time a thread of `tidx` spends away from `excl` in one cycle: its think
       time, its own processor residence, and its blocking at every callee other
       than `excl`. */
    auto cycle_outside = [&](std::size_t tidx, std::size_t excl) -> T {
        T z = thinkt[tidx];
        for (std::size_t j = 0; j < lsn.entriesof[tidx].size(); ++j) {
            const std::size_t eidx = lsn.entriesof[tidx][j];
            T w = T(share[eidx] * residt[eidx]);
            for (std::size_t k = 0; k < callsFrom[eidx].size(); ++k) {
                const std::size_t c = callsFrom[eidx][k];
                if (taskof[calldst[c]] != excl) w += T(share[eidx] * cally[c] * callservt[c]);
            }
            z += w;
        }
        return z;
    };

    /* LINE scales a station utilization into [0,1] whatever its multiplicity,
       and reports busy SERVERS at an infinite server. */
    auto host_servers = [&](std::size_t hidx) -> double {
        return (lsn.sched[hidx] == SchedStrategy::INF) ? 1.0 : npop[hidx];
    };

    /* Reference tasks set the pace; every other rate follows from the call
       rates, so the entries are visited in call-graph order until stable. */
    auto throughputs = [&]() {
        for (std::size_t t = t0; t <= t1; ++t) {
            if (lsn.isref[t]) {
                T cyc = thinkt[t];
                for (std::size_t j = 0; j < lsn.entriesof[t].size(); ++j) {
                    const std::size_t eidx = lsn.entriesof[t][j];
                    cyc += T(share[eidx] * servt[eidx]);
                }
                Xtask[t] = (num_traits<T>::to_double(cyc) > mol_detail::fine_tol())
                               ? T(num_traits<T>::from_double(npop[t]) / cyc)
                               : zero;
                for (std::size_t j = 0; j < lsn.entriesof[t].size(); ++j) {
                    const std::size_t eidx = lsn.entriesof[t][j];
                    Xentry[eidx] = T(Xtask[t] * share[eidx]);
                }
            } else {
                Xtask[t] = zero;
                for (std::size_t j = 0; j < lsn.entriesof[t].size(); ++j)
                    Xentry[lsn.entriesof[t][j]] = zero;
            }
        }
        for (std::size_t pass = 0; pass < std::max<std::size_t>(1, lsn.ntasks); ++pass) {
            for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
                if (lsn.isref[taskof[eidx]]) continue;
                T x = zero;
                for (std::size_t k = 0; k < callsTo[eidx].size(); ++k) {
                    const std::size_t c = callsTo[eidx][k];
                    x += T(cally[c] * Xentry[callsrc[c]]);
                }
                Xentry[eidx] = x;
            }
            for (std::size_t t = t0; t <= t1; ++t) {
                if (lsn.isref[t]) continue;
                const std::vector<std::size_t>& es = lsn.entriesof[t];
                T sum = zero;
                for (std::size_t j = 0; j < es.size(); ++j) sum += Xentry[es[j]];
                Xtask[t] = sum;
                if (num_traits<T>::to_double(sum) > mol_detail::fine_tol())
                    for (std::size_t j = 0; j < es.size(); ++j)
                        share[es[j]] = T(Xentry[es[j]] / sum);
            }
        }
        // Mean busy threads, by Little's law over the entries the task serves.
        // This is the occupancy the think-time closure needs, and it is exact
        // given the throughputs -- unlike the layer AMVA's own U, which is
        // X*L*g with g a reciprocal rate multiplier and not a server count.
        for (std::size_t t = t0; t <= t1; ++t) {
            T u = zero;
            for (std::size_t j = 0; j < lsn.entriesof[t].size(); ++j) {
                const std::size_t eidx = lsn.entriesof[t][j];
                u += T(Xentry[eidx] * servt[eidx]);
            }
            busyth[t] = u;
        }
    };

    std::size_t iter = 0;
    double resid = std::numeric_limits<double>::infinity();
    while (iter < options.iter_max) {
        ++iter;
        const std::vector<T> servt_prev = servt, thinkt_prev = thinkt;

        // ---- phase 1: software layers (thread contention at each called task)
        for (std::size_t li = 0; li < taskLayers.size(); ++li) {
            const std::size_t tidx = taskLayers[li];
            // The task is the station, its caller tasks the classes.
            std::set<std::size_t> callerset;
            for (std::size_t j = 0; j < lsn.entriesof[tidx].size(); ++j)
                for (std::size_t k = 0; k < callsTo[lsn.entriesof[tidx][j]].size(); ++k)
                    callerset.insert(taskof[callsrc[callsTo[lsn.entriesof[tidx][j]][k]]]);
            const std::vector<std::size_t> callers(callerset.begin(), callerset.end());
            const std::size_t Kc = callers.size();
            if (Kc == 0) continue;
            std::vector<T> gcl(Kc, num_traits<T>::from_int(1));
            if (lsn.sched[tidx] != SchedStrategy::INF) {
                // An infinite-thread task never queues for a thread.
                Matrix<T> L(1, Kc, zero);
                std::vector<T> N(Kc, zero), Z(Kc, zero);
                for (std::size_t k = 0; k < Kc; ++k) {
                    const std::size_t ctask = callers[k];
                    N[k] = num_traits<T>::from_double(npop[ctask]);
                    T d = zero;
                    for (std::size_t j = 0; j < lsn.entriesof[ctask].size(); ++j) {
                        const std::size_t eidx = lsn.entriesof[ctask][j];
                        for (std::size_t m = 0; m < callsFrom[eidx].size(); ++m) {
                            const std::size_t c = callsFrom[eidx][m];
                            if (taskof[calldst[c]] == tidx)
                                d += T(share[eidx] * cally[c] * servt[calldst[c]]);
                        }
                    }
                    L(0, k) = d;
                    Z[k] = cycle_outside(ctask, tidx);
                }
                const pfqn::QdAmvaResult<T> r =
                    pfqn::pfqn_qdamva(L, N, Z, mol_detail::mol_mu(N, npop[tidx]), Matrix<T>());
                for (std::size_t k = 0; k < Kc; ++k)
                    if (num_traits<T>::to_double(L(0, k)) > mol_detail::fine_tol())
                        gcl[k] = T(r.R(0, k) / L(0, k));
            }
            for (std::size_t k = 0; k < Kc; ++k) {
                const std::size_t ctask = callers[k];
                for (std::size_t j = 0; j < lsn.entriesof[ctask].size(); ++j) {
                    const std::size_t eidx = lsn.entriesof[ctask][j];
                    for (std::size_t m = 0; m < callsFrom[eidx].size(); ++m) {
                        const std::size_t c = callsFrom[eidx][m];
                        if (taskof[calldst[c]] != tidx) continue;
                        const T newv = T(gcl[k] * servt[calldst[c]]);
                        callservt[c] = T(num_traits<T>::from_double(om) * newv +
                                         num_traits<T>::from_double(1.0 - om) * callservt[c]);
                    }
                }
            }
        }

        // ---- phase 2: hardware layers (processor contention at each host) ---
        for (std::size_t li = 0; li < hostLayers.size(); ++li) {
            const std::size_t hidx = hostLayers[li];
            // The processor is the station, its tasks the classes.
            const std::vector<std::size_t>& tsks = lsn.tasksof[hidx];
            const std::size_t Kt = tsks.size();
            std::vector<T> f(Kt, num_traits<T>::from_int(1));
            if (lsn.sched[hidx] != SchedStrategy::INF) {
                // A delay processor never queues.
                Matrix<T> L(1, Kt, zero);
                std::vector<T> N(Kt, zero), Z(Kt, zero);
                for (std::size_t k = 0; k < Kt; ++k) {
                    const std::size_t tidx = tsks[k];
                    N[k] = num_traits<T>::from_double(npop[tidx]);
                    T d = zero, z = thinkt[tidx];
                    for (std::size_t j = 0; j < lsn.entriesof[tidx].size(); ++j) {
                        const std::size_t eidx = lsn.entriesof[tidx][j];
                        d += T(share[eidx] * dem[eidx]);
                        for (std::size_t m = 0; m < callsFrom[eidx].size(); ++m) {
                            const std::size_t c = callsFrom[eidx][m];
                            z += T(share[eidx] * cally[c] * callservt[c]);
                        }
                    }
                    L(0, k) = d;
                    Z[k] = z;
                }
                const pfqn::QdAmvaResult<T> r =
                    pfqn::pfqn_qdamva(L, N, Z, mol_detail::mol_mu(N, npop[hidx]), Matrix<T>());
                for (std::size_t k = 0; k < Kt; ++k)
                    if (num_traits<T>::to_double(L(0, k)) > mol_detail::fine_tol())
                        f[k] = T(r.R(0, k) / L(0, k));
            }
            for (std::size_t k = 0; k < Kt; ++k)
                for (std::size_t j = 0; j < lsn.entriesof[tsks[k]].size(); ++j) {
                    const std::size_t eidx = lsn.entriesof[tsks[k]][j];
                    residt[eidx] = T(f[k] * dem[eidx]);
                }
        }

        // ---- recompose entry service times ---------------------------------
        for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
            T s = residt[eidx];
            for (std::size_t k = 0; k < callsFrom[eidx].size(); ++k) {
                const std::size_t c = callsFrom[eidx][k];
                s += T(cally[c] * callservt[c]);
            }
            servt[eidx] = T(num_traits<T>::from_double(om) * s +
                            num_traits<T>::from_double(1.0 - om) * servt[eidx]);
        }

        // ---- throughputs, entry shares, think-time closure ------------------
        throughputs();
        for (std::size_t tidx = t0; tidx <= t1; ++tidx) {
            if (lsn.isref[tidx]) {
                thinkt[tidx] = zref[tidx];
                continue;
            }
            if (num_traits<T>::to_double(Xtask[tidx]) <= mol_detail::fine_tol()) continue;
            // Idle time of a thread per cycle. updateThinkTimes splits this into
            // an INF arm (njobs - util) and a finite arm (njobs*abs(1-util))
            // only because LINE reports busy SERVERS at an infinite server and a
            // busy FRACTION at a finite one; carrying the count in both cases
            // makes the two arms the same expression.
            const double v = std::fabs(npop[tidx] - num_traits<T>::to_double(busyth[tidx])) /
                                 num_traits<T>::to_double(Xtask[tidx]) -
                             num_traits<T>::to_double(zref[tidx]);
            const T newz = num_traits<T>::from_double(std::max(0.0, v));
            thinkt[tidx] = T(num_traits<T>::from_double(om) * newz +
                             num_traits<T>::from_double(1.0 - om) * thinkt[tidx]);
        }

        // Both halves of the state must settle: servt alone can sit still for
        // an iteration while the think times are still moving.
        resid = 0.0;
        for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
            const double a = num_traits<T>::to_double(servt[eidx]);
            const double b = num_traits<T>::to_double(servt_prev[eidx]);
            resid = std::max(resid, std::fabs(a - b) / std::max(1.0, std::fabs(a)));
        }
        for (std::size_t tidx = t0; tidx <= t1; ++tidx) {
            const double a = num_traits<T>::to_double(thinkt[tidx]);
            const double b = num_traits<T>::to_double(thinkt_prev[tidx]);
            resid = std::max(resid, std::fabs(a - b) / std::max(1.0, std::fabs(a)));
        }
        if (resid < options.iter_tol) break;
    }
    throughputs();

    // ---- assemble the reported vectors -------------------------------------
    const double nan = std::numeric_limits<double>::quiet_NaN();
    LqnMolResult<T> out;
    out.QN.assign(nidx + 1, num_traits<T>::from_double(nan));
    out.UN = out.RN = out.TN = out.QN;
    for (std::size_t eidx = e0; eidx <= e1; ++eidx) {
        const std::size_t hidx = hostof[taskof[eidx]];
        const T procutil =
            T(Xentry[eidx] * dem[eidx] / num_traits<T>::from_double(host_servers(hidx)));
        out.QN[eidx] = T(Xentry[eidx] * servt[eidx]);
        out.UN[eidx] = procutil;
        out.RN[eidx] = servt[eidx];
        out.TN[eidx] = Xentry[eidx];
        const std::size_t aidx = actof[eidx];
        out.QN[aidx] = out.QN[eidx];
        out.UN[aidx] = out.UN[eidx];
        out.RN[aidx] = out.RN[eidx];
        out.TN[aidx] = out.TN[eidx];
    }
    for (std::size_t tidx = t0; tidx <= t1; ++tidx) {
        T q = zero, u = zero;
        for (std::size_t j = 0; j < lsn.entriesof[tidx].size(); ++j) {
            q += out.QN[lsn.entriesof[tidx][j]];
            u += out.UN[lsn.entriesof[tidx][j]];
        }
        out.QN[tidx] = q;
        out.UN[tidx] = u;
        out.RN[tidx] = num_traits<T>::from_double(nan);
        out.TN[tidx] = Xtask[tidx];
    }
    for (std::size_t hidx = 1; hidx <= lsn.nhosts; ++hidx) {
        T u = zero;
        for (std::size_t j = 0; j < lsn.tasksof[hidx].size(); ++j)
            u += out.UN[lsn.tasksof[hidx][j]];
        out.QN[hidx] = num_traits<T>::from_double(nan);
        out.UN[hidx] = u;
        out.RN[hidx] = num_traits<T>::from_double(nan);
        // No throughput is defined at a processor, as in LQNS.
        out.TN[hidx] = num_traits<T>::from_double(nan);
    }

    out.info.iter = iter;
    out.info.resid = resid;
    out.info.servt = servt;
    out.info.residt = residt;
    out.info.callservt = callservt;
    out.info.thinkt = thinkt;
    out.info.share = share;
    out.info.hostLayers = hostLayers;
    out.info.taskLayers = taskLayers;
    return out;
}

}  // namespace lqn
}  // namespace line

#endif  // LINE_API_LQN_MOL_H
