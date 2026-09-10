/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_MMAP_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_MMAP_H

/**
 * Port of `solver_mam_basic_mmap.m`, `solver_mam_basic_mmap_inner.m` and
 * `solver_mam_basic_mmap_closed.m`: the MMAP fork-join decomposition, reached as
 * the `dec.source.mmap` method and as branch 2b of the dispatch (every open
 * fork-join model that is not in the homogeneous class `solver_mam_fj` serves).
 *
 * THE METHOD. It is a PARAMETRIC DECOMPOSITION, not a per-station isolation:
 * unlike `dec.source`, which hands every station a rescaled copy of the chain's
 * source process, this one carries a per-node DEPARTURE process table and
 * recomputes the arrival stream at every node from the traffic equations each
 * sweep (`solver_mam_traffic_mmap`, the fork-join aware traffic step). The
 * departure process of an FCFS or PS station is the ETAQA truncation of its own
 * QBD (`qbd_depproc_etaqa`, `qbd_depproc_etaqa_ps`), so the correlation a queue
 * introduces travels downstream instead of being discarded. The fixed point is
 * driven on the station queue lengths by `da_fpi` with a RELATIVE increment
 * norm, and the reference starts testing convergence only from the third sweep
 * (`config.da_miniter = 3`).
 *
 * WHAT MAKES IT THE FORK-JOIN ANALYZER. Two things, and neither is in
 * `dec.source`:
 *  - the traffic step SYNCHRONIZES the flows arriving at a join along one sync
 *    group with `mmap_max` rather than superposing them, so the join's output
 *    process is the slowest branch's, blocking included;
 *  - the join's own metrics are derived AFTER the fixed point from the branch
 *    response times, as the expected maximum of independent exponentials with
 *    rates 1/R_b minus their mean. That difference is the synchronization delay,
 *    and `QN = (sum of branch throughputs) * delay` is Little's law on it.
 *
 * THE CLOSED WRAPPER has no source to fix the arrival rates, so it wraps the
 * inner analyzer in a per-class BISECTION on a surrogate arrival rate, driven
 * against the population, exactly as `solver_mna_closed` does. Three details of
 * the reference are reproduced rather than tidied: the bracket's upper bound is
 * the SLOWEST rate over the finite-server stations (the fallback to the
 * infinite-server ones, and then to 1, is the reference's own); the loop breaks
 * when every bracket has collapsed below the precision floor, undoing its own
 * iteration count as it does so; and a diverged inner call is caught, treated as
 * an overload, and the last successful metrics are restored if the FINAL trial
 * is the one that diverged.
 *
 * SELF-LOOPING CLASSES. `sn.isslc` guards the surrogate-rate zeroing, the queue
 * clamp and the final throughput pin in the reference's closed wrapper, and the
 * PS denominator and the FCFS saturation test in the inner analyzer. The C++
 * `JobClassType` is OPEN or CLOSED only, so no model this port can build enters
 * them, and they are not transcribed -- the same decision `solver_mna.h` and
 * `solver_mam_ag.h` record.
 *
 * THE REFERENCE'S OWN QUIRKS, REPRODUCED. `XN` is initialised to zeros and never
 * assigned by either the inner analyzer or the closed wrapper, so the per-class
 * throughput column is zero however the model is solved; the station
 * throughputs in `TN` are the real ones. The inner analyzer's `try/catch` around
 * the ETAQA departure process, and the closed wrapper's around the whole inner
 * call, are the reference's control flow and not defensive additions: the first
 * falls back to the scaled service process, the second to a bisection step
 * downwards.
 *
 * ARITHMETIC. Double (or real) only, for the reasons `solver_mam_basic.h` lists:
 * the fixed point stops on a tolerance, MMAP[K]/PH[K]/1 runs the ADDA doubling
 * iteration, and both ETAQA departure processes static_assert on transcendental
 * arithmetic.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/mam/qbd_depproc.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mam/solver_mam_bmap.h"  // mam_detect_mmck
#include "line/solvers/mam/solver_mam_traffic.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The `options.config` fields the MMAP decomposition reads on top of
 * `MamOptions`.
 *
 * Neither is a SolverMAM option: `etaqa_trunc` is defaulted by
 * `solver_mam_analyzer.m` before the dispatch and `fj_sync_q_len` by
 * `solver_mam_traffic_mmap.m` itself, so a caller cannot reach either through
 * the solver's option surface. Kept out of `MamOptions` for that reason, as
 * `MnaConfig` is.
 */
struct MmapDecConfig {
    /** `config.etaqa_trunc`, the ETAQA level truncation of the departure process. */
    std::size_t etaqa_trunc = 8;
    /** `config.fj_sync_q_len`, the synchronization queue at a join. */
    std::size_t fj_sync_q_len = 2;
};

namespace basic_mmap_detail {

using lang::GlobalConstants;
using lang::SchedStrategy;

/**
 * The `PH`, `pie` and `D0` tables both MMAP analyzers build before the loop.
 *
 * `PH[i][r]` is the service process of class r at station i AFTER the
 * reference's per-discipline rescaling (divided by the server count at FCFS,
 * HOL and PS, left alone at INF and at a Source, whose entry is the ARRIVAL
 * process the reference never overwrites), `svc[i][r]` is the same law as the
 * (sigma, S) pair `MMAPPH1FCFS` consumes, and `known[i][r]` says whether the
 * station serves the class at all.
 *
 * The NaN guard below is applied at EVERY station, where the reference applies
 * it only inside its four rescaling branches. The difference is invisible: an
 * entry it touches is one whose class the station never serves, whose
 * throughput is therefore zero, so the reference's NaN and this port's mean of
 * 1e8 both leave the product at zero -- the first through the terminal NaN
 * sweep, the second directly.
 *
 * THE REFERENCE'S NaN GUARD IS INCONSISTENT, AND IS REPRODUCED AS WRITTEN. When
 * `D0` comes back NaN -- the class is not served here -- it sets
 * `D0 = -GlobalConstants.Immediate` (an IMMEDIATE service, rate 1e8) but
 * `PH = map_exponential(GlobalConstants.Immediate)`, and `map_exponential` takes
 * a MEAN, so that is a service of mean 1e8, the exact opposite. The two are read
 * in different places -- `D0`/`pie` by the queue solver, `PH` only through
 * `map_mean` in the utilization and surrogate-delay lines -- and every metric
 * that reads `PH` multiplies it by a throughput that is identically zero on such
 * a class, so the contradiction never reaches a reported number. Written the
 * reference's way rather than tidied, because tidying it would change what
 * MMAPPH1FCFS is handed.
 */
template <class T>
struct MmapPhTable {
    std::vector<std::vector<Map<T>>> PH;
    std::vector<std::vector<PhService<T>>> svc;
    std::vector<std::vector<bool>> known;
};

template <class T>
MmapPhTable<T> mmap_ph_table(const qn::NetworkStruct<T>& L) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
    MmapPhTable<T> t;
    t.PH.assign(M, std::vector<Map<T>>(K));
    t.svc.assign(M, std::vector<PhService<T>>(K));
    t.known.assign(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        // EVERY STATION IS FILLED, not only the ones with a branch. The
        // reference's `PH = sn.proc` starts as the whole table and its station
        // loop only OVERWRITES the entries it rescales, so a Source keeps its
        // own ARRIVAL process there -- and both analyzers read it: `dec.mmap`'s
        // trailing surrogate-delay block runs at every station, and the MMAP
        // variant seeds DEP from PH at every station node.
        // FCFS, HOL and PS are the disciplines the reference divides by the
        // server count; INF is not (an infinite server has no queue to speed
        // up), and neither is a Source.
        const bool divide = (sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL ||
                             sc == SchedStrategy::PS);
        const double ns = L.stations[i].nservers;
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.has_service_law(i, r) || !(L.rates(i, r) > zero)) {
                t.PH[i][r] = map_exponential_mean(imm);  // mean 1e8; see the struct note
                t.svc[i][r].sigma.assign(1, one);
                t.svc[i][r].S = Matrix<T>(1, 1, T(-imm));
                continue;
            }
            Map<T> ph = lang::dist_to_map(L.service[i][r]);
            if (divide && std::isfinite(ns) && ns > 0.0)
                ph = map_scale(ph, T(map_mean(ph) / num_traits<T>::from_double(ns)));
            t.PH[i][r] = ph;
            t.svc[i][r].sigma = map_pie(ph);
            t.svc[i][r].S = ph.D0;
            t.known[i][r] = true;
        }
    }
    return t;
}

/**
 * The expected maximum of independent exponentials, by inclusion-exclusion.
 *
 * The reference writes it as
 *   sum_{p=0}^{n-1} (-1)^p sum(1 ./ sum(nchoosek(lambda, p+1), 2))
 * i.e. the alternating sum over subsets of every size, which is the classical
 * E[max] = sum_{S nonempty} (-1)^(|S|+1) / sum_{i in S} lambda_i. Enumerated
 * here over bitmasks, which is the same set of subsets in the same signs; the
 * cost is 2^n either way, and the reference's own nchoosek is what bounds n in
 * practice. A join with more than 20 synchronized branches is refused by name
 * rather than run into a 10^6-term sum.
 */
template <class T>
T exp_max_mean(const std::vector<T>& rate) {
    const std::size_t n = rate.size();
    if (n > 20)
        throw UnsupportedError(
            "solver_mam_basic_mmap: the join synchronizes " + std::to_string(n) +
            " branches, and the reference's expected-maximum formula is an alternating sum over "
            "all 2^n subsets of them; that is not evaluable at this width");
    const T zero = num_traits<T>::from_int(0);
    T acc = zero;
    for (std::size_t mask = 1; mask < (std::size_t(1) << n); ++mask) {
        T s = zero;
        std::size_t bits = 0;
        for (std::size_t i = 0; i < n; ++i)
            if (mask & (std::size_t(1) << i)) {
                s += rate[i];
                ++bits;
            }
        if (!(s > zero)) continue;
        const T term = T(num_traits<T>::from_int(1) / s);
        if (bits % 2 == 1) acc += term;
        else acc -= term;
    }
    return acc;
}

/** `mmap_compress(ARV{ind}, config)` under the arithmetic gate the tree uses. */
template <class T>
Mmap<T> compress_arrival(const Mmap<T>& m) {
    if constexpr (num_traits<T>::has_transcendental) {
        return mmap_compress(m, MmapCompressMethod::MixtureOrder1);
    } else {
        (void)m;
        throw UnsupportedError(
            "solver_mam_basic_mmap: the arrival superposition passed config.space_max and must be "
            "compressed, which fits an APH(2) and needs transcendental arithmetic");
    }
}

}  // namespace basic_mmap_detail

/**
 * Port of `solver_mam_basic_mmap_inner.m`.
 *
 * @param L       the refreshed struct
 * @param opt     the SolverMAM options; `tol` is the fixed point's iter_tol
 * @param cfg     the two `options.config` fields the analyzer defaults itself
 * @param lambda  per-CLASS surrogate arrival rate, the reference's `lambda`
 * @param totiter out: the sweeps the fixed point took
 */
template <class T>
mva::MvaSolution<T> solver_mam_basic_mmap_inner(const qn::NetworkStruct<T>& L,
                                                const MamOptions& opt, const MmapDecConfig& cfg,
                                                const std::vector<T>& lambda,
                                                std::size_t* totiter) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L; (void)opt; (void)cfg; (void)lambda; (void)totiter;
        throw UnsupportedError(
            "solver_mam_basic_mmap_inner: the departure-process fixed point stops on a tolerance, "
            "MMAP[K]/PH[K]/1 runs the ADDA doubling iteration and the ETAQA departure process "
            "needs transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
    using namespace basic_mmap_detail;
    using basic_detail::station_visits;
    using basic_detail::truncate_renorm;
    using basic_detail::zero_nans;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t I = L.nof_nodes(), M = L.nstations, K = L.nclasses;
    const T ftol = num_traits<T>::from_double(GlobalConstants::FineTol);

    if (lambda.size() != K)
        throw InputError("solver_mam_basic_mmap_inner: lambda is not indexed over the classes");

    const Matrix<T> V = station_visits(L);
    Matrix<T> S(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            if (!L.disabled[i][r] && L.rates(i, r) > zero) S(i, r) = T(one / L.rates(i, r));

    const MmapPhTable<T> ph = mmap_ph_table(L);
    const FjSyncMap fj = sn_build_fj_sync_map(L);
    TrafficConfig tcfg = traffic_config(opt);
    tcfg.fj_sync_q_len = cfg.fj_sync_q_len;

    Matrix<T> QN(M, K, zero), UN(M, K, zero), RN(M, K, zero), TN(M, K, zero);
    // The Source's throughput is its declared rate and is set once, before the
    // loop, exactly as the reference's pre-loop station switch does.
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].sched == SchedStrategy::EXT)
            for (std::size_t r = 0; r < K; ++r)
                TN(i, r) = L.disabled[i][r] ? zero : L.rates(i, r);

    DepTable<T> DEP(I, std::vector<Map<T>>(K));

    auto sweep = [&](const std::vector<T>&,
                     std::size_t itnum) -> std::pair<std::vector<T>, std::vector<T>> {
        // ---- the departure table, NODE-indexed --------------------------
        if (itnum == 1) {
            for (std::size_t ind = 0; ind < I; ++ind) {
                const qn::NodeType ty = L.nodes[ind].nodetype;
                const bool isfj = (ty == qn::NodeType::Fork || ty == qn::NodeType::Join);
                if (L.nodes[ind].station != 0 && !isfj) {
                    const std::size_t ist = L.nodes[ind].station - 1;
                    for (std::size_t r = 0; r < K; ++r) {
                        if (V(ist, r) > zero && lambda[r] > zero)
                            DEP[ind][r] =
                                map_scale(ph.PH[ist][r], T(one / T(lambda[r] * V(ist, r))));
                        else
                            DEP[ind][r] = ph.PH[ist][r];
                    }
                } else {
                    for (std::size_t r = 0; r < K; ++r)
                        DEP[ind][r] = map_exponential_mean(
                            lambda[r] > zero
                                ? T(one / lambda[r])
                                : T(one / num_traits<T>::from_double(GlobalConstants::Immediate)));
                }
            }
        }

        std::vector<Mmap<T>> ARV = solver_mam_traffic_mmap(L, DEP, tcfg, fj);

        std::vector<T> xref(M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) xref[i * K + r] = QN(i, r);

        // ---- one isolated-station solve per station ----------------------
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t ind = L.node_of_station(i + 1) - 1;
            const qn::NodeType ty = L.nodes[ind].nodetype;
            const SchedStrategy sc = L.stations[i].sched;
            const double ns = L.stations[i].nservers;

            if (ty == qn::NodeType::Join) {
                // Zeroed here and rebuilt from the branch response times after
                // the fixed point; the throughput is the surrogate rate.
                for (std::size_t r = 0; r < K; ++r) {
                    TN(i, r) = lambda[r];
                    UN(i, r) = zero;
                    QN(i, r) = zero;
                    RN(i, r) = zero;
                }
                continue;
            }
            if (ty != qn::NodeType::Queue) {
                if (sc == SchedStrategy::INF) {
                    if (ARV[ind].order() > 0) {
                        const std::vector<T> lam = mmap_lambda(ARV[ind]);
                        for (std::size_t r = 0; r < K; ++r) TN(i, r) = lam[r];
                    }
                    for (std::size_t r = 0; r < K; ++r)
                        if (TN(i, r) > zero) {
                            UN(i, r) = T(S(i, r) * TN(i, r));
                            QN(i, r) = T(TN(i, r) * S(i, r));
                            RN(i, r) = S(i, r);
                        }
                }
                // EXT: the throughput was pinned before the loop.
                continue;
            }
            if (ARV[ind].order() == 0) continue;
            if (ARV[ind].order() > tcfg.space_max) ARV[ind] = compress_arrival(ARV[ind]);

            bool finiteCapUsed = false;
            std::vector<PhService<T>> sl;
            for (std::size_t r = 0; r < K; ++r) sl.push_back(ph.svc[i][r]);

            if (sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL) {
                if (std::isfinite(L.stations[i].cap)) {
                    const std::size_t capK =
                        static_cast<std::size_t>(std::llround(L.stations[i].cap));
                    T meanQ = zero, lossProb = zero;
                    const MmckDetection<T> det = mam_detect_mmck(L, i + 1, ARV[ind]);
                    if (det.isMmck) {
                        const std::vector<T> lam = mmap_lambda(ARV[ind]);
                        T lamTot = zero;
                        for (const T& v : lam)
                            if (!std::isnan(num_traits<T>::to_double(v))) lamTot += v;
                        const qsys::MmckResult<T> ex = qsys::qsys_mmck(
                            lamTot, det.muRate, static_cast<unsigned>(std::llround(ns)),
                            static_cast<unsigned>(capK));
                        meanQ = ex.meanQueueLength;
                        lossProb = ex.lossProbability;
                    } else {
                        const basic_detail::TruncRenorm<T> tr = truncate_renorm(ARV[ind], sl, capK);
                        meanQ = tr.meanQ;
                        lossProb = tr.lossProb;
                    }
                    const std::vector<T> lam = mmap_lambda(ARV[ind]);
                    std::vector<T> eff(K, zero), Sact(K, zero);
                    T sumTN = zero;
                    for (std::size_t r = 0; r < K; ++r) {
                        const T inflow =
                            std::isnan(num_traits<T>::to_double(lam[r])) ? zero : lam[r];
                        eff[r] = T(inflow * T(one - lossProb));
                        // The PH was divided by the server count, so the actual
                        // per-class service mean multiplies it back.
                        Sact[r] = T(map_mean(ph.PH[i][r]) * num_traits<T>::from_double(ns));
                        sumTN += eff[r];
                    }
                    T Wq = zero;
                    if (sumTN > zero) {
                        T sw = zero;
                        for (std::size_t r = 0; r < K; ++r) {
                            const T c = T(eff[r] * Sact[r]);
                            if (!std::isnan(num_traits<T>::to_double(c))) sw += c;
                        }
                        const T w = T(T(meanQ / sumTN) - T(sw / sumTN));
                        Wq = (w > zero) ? w : zero;
                    }
                    for (std::size_t r = 0; r < K; ++r) {
                        TN(i, r) = eff[r];
                        UN(i, r) = T(TN(i, r) * map_mean(ph.PH[i][r]));
                        if (TN(i, r) > zero) {
                            RN(i, r) = T(Wq + Sact[r]);
                            QN(i, r) = T(TN(i, r) * RN(i, r));
                        } else {
                            RN(i, r) = zero;
                            QN(i, r) = zero;
                        }
                    }
                    finiteCapUsed = true;
                } else {
                    const std::vector<T> lam = mmap_lambda(ARV[ind]);
                    T rho = zero;
                    for (std::size_t r = 0; r < K; ++r) {
                        const T u = T(lam[r] * map_mean(ph.PH[i][r]));
                        if (!std::isnan(num_traits<T>::to_double(u))) rho += u;
                    }
                    if (rho < T(one - ftol)) {
                        // A correlated single-class service is answered by the
                        // exact MAP/MAP/1 QBD, which carries the service phase
                        // across departures; MMAPPH1FCFS would discard it.
                        const bool corr =
                            (K == 1) && (ns == 1.0) && ph.known[i][0] &&
                            std::fabs(num_traits<T>::to_double(
                                map_acf(ph.PH[i][0], std::vector<unsigned>{1})[0])) >
                                GlobalConstants::CoarseTol;
                        if (corr) {
                            Map<T> arv;
                            arv.D0 = ARV[ind].D0;
                            arv.D1 = ARV[ind].Dc[0];
                            QN(i, 0) = qbd_mapmap1(arv, ph.PH[i][0]).QN;
                        } else {
                            const std::vector<T> m = mmapph1fcfs_ncmean(ARV[ind], sl);
                            for (std::size_t r = 0; r < K; ++r)
                                QN(i, r) = m[ARV[ind].classes() == 1 ? 0 : r];
                        }
                    } else {
                        // Bounded rather than left to diverge: an overloaded
                        // station holds its class's own population, or 1/FineTol
                        // when the class is open.
                        for (std::size_t r = 0; r < K; ++r)
                            QN(i, r) = std::isfinite(L.classes[r].population)
                                           ? num_traits<T>::from_double(L.classes[r].population)
                                           : T(one / ftol);
                    }
                    for (std::size_t r = 0; r < K; ++r) TN(i, r) = lam[r];
                }
            } else if (sc == SchedStrategy::PS) {
                const std::vector<T> lam = mmap_lambda(ARV[ind]);
                for (std::size_t r = 0; r < K; ++r) {
                    TN(i, r) = lam[r];
                    // S, NOT the server-count-scaled PH mean: the reference
                    // writes 1./sn.rates here and overwrites UN with the PH mean
                    // in the surrogate-delay block below, so this value only ever
                    // reaches the sharing denominator.
                    UN(i, r) = T(TN(i, r) * S(i, r));
                }
                T usum = zero;
                for (std::size_t r = 0; r < K; ++r) usum += UN(i, r);
                const T uden = (usum < T(one - ftol)) ? usum : T(one - ftol);
                for (std::size_t r = 0; r < K; ++r) QN(i, r) = T(UN(i, r) / T(one - uden));
            }

            if (!finiteCapUsed) {
                for (std::size_t r = 0; r < K; ++r) {
                    UN(i, r) = T(TN(i, r) * map_mean(ph.PH[i][r]));
                    // The jobs at the surrogate delay server the c-fold service
                    // speedup removed.
                    if (std::isfinite(ns))
                        QN(i, r) = T(QN(i, r) + TN(i, r) *
                                                    T(map_mean(ph.PH[i][r]) *
                                                      num_traits<T>::from_double(ns)) *
                                                    num_traits<T>::from_double((ns - 1.0) / ns));
                    RN(i, r) = T(QN(i, r) / TN(i, r));
                }
            }
        }

        // ---- the departure processes for the next sweep -------------------
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t ind = L.node_of_station(i + 1) - 1;
            const qn::NodeType ty = L.nodes[ind].nodetype;
            const SchedStrategy sc = L.stations[i].sched;
            if (ty == qn::NodeType::Join) {
                for (std::size_t r = 0; r < K; ++r)
                    if (TN(i, r) > zero) DEP[ind][r] = map_exponential_mean(T(one / TN(i, r)));
                continue;
            }
            if (ty != qn::NodeType::Queue || ARV[ind].order() == 0) continue;
            const bool fcfs = (sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL);
            if (!fcfs && sc != SchedStrategy::PS) continue;

            T rho = zero;
            for (std::size_t r = 0; r < K; ++r) rho += UN(i, r);
            for (std::size_t r = 0; r < K; ++r) {
                const bool scalable = (V(i, r) > zero && lambda[r] > zero);
                // The PS branch of the reference does nothing at all when the
                // class does not flow through the station; the FCFS branch still
                // recomputes the departure process and only skips the rescaling.
                if (!fcfs && !scalable) continue;
                const Mmap<T> A = mmap_hide_but(ARV[ind], r);
                const Map<T>& Srv = ph.PH[i][r];
                const std::size_t etaqa_sz =
                    (cfg.etaqa_trunc + 1) * A.order() * Srv.D0.rows();
                Map<T> dep = Srv;
                if (etaqa_sz <= tcfg.space_max && rho < T(one - ftol)) {
                    // The reference's own try/catch: an ETAQA truncation that
                    // fails to build falls back to the scaled service process.
                    try {
                        const Map<T> Am{A.D0, A.D1};
                        dep = map_normalize(fcfs ? qbd_depproc_etaqa(Am, Srv, cfg.etaqa_trunc)
                                                 : qbd_depproc_etaqa_ps(Am, Srv, cfg.etaqa_trunc));
                    } catch (const Error&) {
                        dep = Srv;
                    }
                }
                if (scalable) dep = map_scale(dep, T(one / T(lambda[r] * V(i, r))));
                DEP[ind][r] = dep;
            }
        }

        std::vector<T> xnew(M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) xnew[i * K + r] = QN(i, r);
        return std::make_pair(xnew, xref);
    };

    da::FpiOptions fo;
    fo.iter_max = static_cast<std::size_t>(opt.iter_max);
    fo.iter_tol = opt.tol;
    // `config.da_miniter = 3`: the legacy loop tested convergence only from the
    // third sweep, and `config.da_norm` is the relative difference offset by
    // FineTol so a zero reference entry does not divide by zero.
    fo.miniter = 3;
    fo.relative_norm = true;
    fo.relative_eps = GlobalConstants::FineTol;
    const da::FpiResult<T> fr = da::da_fpi<T>(sweep, std::vector<T>(M * K, zero), fo);
    if (totiter) *totiter = fr.iterations;

    // ---- the join, from the branch response times -------------------------
    for (std::size_t j = 0; j < I; ++j) {
        if (L.nodes[j].nodetype != qn::NodeType::Join) continue;
        if (L.nodes[j].station == 0) continue;
        const std::size_t jst = L.nodes[j].station - 1;
        std::map<std::size_t, std::vector<std::size_t>> groups;  // ordered, as unique() is
        for (std::size_t b = 0; b < I; ++b) {
            const std::size_t gid = fj.node_sync[j][b];
            if (gid > 0) groups[gid].push_back(b);
        }
        for (std::size_t r = 0; r < K; ++r) {
            if (!(TN(jst, r) > zero)) continue;
            T syncDelay = zero, joinArrivalRate = zero;
            for (const std::pair<const std::size_t, std::vector<std::size_t>>& g : groups) {
                std::vector<T> branchRt, branchTput;
                for (std::size_t b : g.second) {
                    if (L.nodes[b].station == 0) continue;
                    const std::size_t bst = L.nodes[b].station - 1;
                    if (!(RN(bst, r) > zero)) continue;
                    branchRt.push_back(RN(bst, r));
                    branchTput.push_back(TN(bst, r));
                }
                if (branchRt.size() < 2) continue;
                std::vector<T> rate(branchRt.size(), zero);
                T meanRt = zero;
                for (std::size_t b = 0; b < branchRt.size(); ++b) {
                    rate[b] = T(one / branchRt[b]);
                    meanRt += branchRt[b];
                }
                meanRt /= num_traits<T>::from_int(static_cast<long>(branchRt.size()));
                const T excess = T(exp_max_mean(rate) - meanRt);
                if (excess > zero) syncDelay += excess;
                for (const T& t : branchTput) joinArrivalRate += t;
            }
            RN(jst, r) = syncDelay;
            QN(jst, r) = T(joinArrivalRate * syncDelay);
            UN(jst, r) = zero;
        }
    }

    mva::MvaSolution<T> out;
    out.Q = QN;
    out.U = UN;
    out.R = RN;
    out.Tp = TN;
    out.C.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t i = 0; i < M; ++i) out.C[r] += RN(i, r);
    // X is left at zero: the reference never assigns it. See the file header.
    out.X.assign(K, zero);
    zero_nans(out.Q);
    zero_nans(out.U);
    zero_nans(out.R);
    zero_nans(out.Tp);
    for (std::size_t r = 0; r < K; ++r)
        if (std::isnan(num_traits<T>::to_double(out.C[r]))) out.C[r] = zero;
    out.method = "dec.source.mmap";
    out.iter = static_cast<int>(fr.iterations);
    out.lG = 0.0;
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mam_basic_mmap_closed.m`: the per-class bisection on the
 * surrogate arrival rate that makes the inner analyzer's queue lengths match the
 * closed population.
 */
template <class T>
mva::MvaSolution<T> solver_mam_basic_mmap_closed(const qn::NetworkStruct<T>& L,
                                                 const MamOptions& opt, const MmapDecConfig& cfg,
                                                 std::size_t* totiter) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L; (void)opt; (void)cfg; (void)totiter;
        throw UnsupportedError(
            "solver_mam_basic_mmap_closed: the throughput bisection wraps the inner MMAP "
            "decomposition, which needs transcendental arithmetic; rerun with --arith double or "
            "--arith real");
    } else {
    using namespace basic_mmap_detail;
    using basic_detail::zero_nans;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;

    // REFERENCE INDEXING, CHECKED RATHER THAN ASSUMED, exactly as
    // `solver_mna_closed` checks the same line: the terminal renormalization
    // scales chain c's queue lengths to `sn.njobs(c)`, and `sn.njobs` is
    // CLASS-indexed. That is only correct when each chain holds exactly one
    // class; otherwise a chain would be scaled to another class's population,
    // which is a wrong number with no symptom.
    if (C != K)
        throw UnsupportedError(
            "solver_mam_basic_mmap_closed: the reference renormalizes chain c's queue lengths with "
            "the class-indexed sn.njobs(c), which is only correct when each chain holds exactly "
            "one class; this model has " +
            std::to_string(C) + " chains over " + std::to_string(K) + " classes");

    // ---- the bisection bracket -------------------------------------------
    // The upper bound is the SLOWEST rate the class meets at a finite-server
    // station; with none, the reference falls back to the FASTEST rate over the
    // infinite-server ones, and to 1 when the class is served nowhere.
    std::vector<T> lambda_lb(K, zero), lambda_ub(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        bool any = false;
        double best = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            if (!std::isfinite(L.stations[i].nservers)) continue;
            if (L.disabled[i][r] || !(L.rates(i, r) > zero)) continue;
            const double v = num_traits<T>::to_double(L.rates(i, r));
            if (!any || v < best) { best = v; any = true; }
        }
        if (!any) {
            for (std::size_t i = 0; i < M; ++i) {
                if (std::isfinite(L.stations[i].nservers)) continue;
                if (L.disabled[i][r] || !(L.rates(i, r) > zero)) continue;
                const double v = num_traits<T>::to_double(L.rates(i, r));
                if (!any || v > best) { best = v; any = true; }
            }
        }
        lambda_ub[r] = any ? num_traits<T>::from_double(best) : one;
    }
    std::vector<T> lambda = lambda_ub;

    std::vector<T> QNc(K, zero);
    for (std::size_t r = 0; r < K; ++r)
        // An open class contributes 0: only the closed populations gate the loop.
        QNc[r] = std::isfinite(L.classes[r].population)
                     ? num_traits<T>::from_double(L.classes[r].population)
                     : zero;
    std::vector<T> QN_chain(K, zero);

    MamOptions inner = opt;
    inner.iter_max = std::max(20, (opt.iter_max + 9) / 10);
    // The reference caps the MMAP phase truncation at 16 here: the inner call is
    // made once per bisection step, and the compression above that width is
    // O(dim^3) for no accuracy the bisection can use.
    if (inner.space_max > 16) inner.space_max = 16;

    mva::MvaSolution<T> sol, last;
    bool have_good = false, algorithm_ok = false;
    const double bisect_tol = std::max(opt.tol, 1e-3);
    int it_out = 0;

    for (;;) {
        double gap = 0.0;
        for (std::size_t r = 0; r < K; ++r)
            gap = std::max(gap, std::fabs(num_traits<T>::to_double(T(QN_chain[r] - QNc[r]))));
        if (!(gap > bisect_tol) || it_out >= opt.iter_max) break;
        ++it_out;
        if (it_out > 1) {
            bool bracket_collapsed = true;
            for (std::size_t r = 0; r < K; ++r) {
                if (!std::isfinite(L.classes[r].population) || !(QNc[r] > zero)) continue;
                if (QN_chain[r] < QNc[r]) lambda_lb[r] = lambda[r];
                else lambda_ub[r] = lambda[r];
                lambda[r] = T(T(lambda_lb[r] + lambda_ub[r]) / num_traits<T>::from_int(2));
                // Bisection can still refine class r only while its bracket is
                // wider than the precision floor below which lambda cannot move
                // any reported metric.
                const double width =
                    num_traits<T>::to_double(T(lambda_ub[r] - lambda_lb[r]));
                if (width > GlobalConstants::FineTol *
                                std::max(1.0, std::fabs(num_traits<T>::to_double(lambda_ub[r]))))
                    bracket_collapsed = false;
            }
            if (bracket_collapsed) {
                --it_out;
                break;
            }
        }

        try {
            std::size_t inner_iter = 0;
            sol = solver_mam_basic_mmap_inner(L, inner, cfg, lambda, &inner_iter);
            algorithm_ok = true;
        } catch (const Error&) {
            // The reference's own catch: the inner algorithm diverged (typically
            // MMAPPH1FCFS under saturation). Every chain is treated as
            // overloaded, so the next bisection step drops lambda.
            algorithm_ok = false;
        }

        if (algorithm_ok) {
            for (std::size_t r = 0; r < K; ++r) {
                QN_chain[r] = zero;
                for (std::size_t i = 0; i < M; ++i) QN_chain[r] += sol.Q(i, r);
                if (!std::isfinite(num_traits<T>::to_double(QN_chain[r])))
                    QN_chain[r] = T(one / num_traits<T>::from_double(GlobalConstants::FineTol));
            }
            last = sol;
            have_good = true;
        } else {
            for (std::size_t r = 0; r < K; ++r)
                QN_chain[r] = T(one / num_traits<T>::from_double(GlobalConstants::FineTol));
        }
    }

    // If the LAST trial diverged, fall back to the most recent successful one.
    if (!algorithm_ok && have_good) sol = last;
    if (sol.Q.rows() != M) {
        // The loop never ran a trial: the initial gap was already inside the
        // tolerance, which for a model of zero population is the honest answer.
        sol.Q = Matrix<T>(M, K, zero);
        sol.U = Matrix<T>(M, K, zero);
        sol.R = Matrix<T>(M, K, zero);
        sol.Tp = Matrix<T>(M, K, zero);
        sol.C.assign(K, zero);
        sol.X.assign(K, zero);
    }

    // ---- population redistribution within each chain ----------------------
    for (std::size_t c = 0; c < C; ++c) {
        if (!std::isfinite(L.classes[c].population)) continue;
        T sumQ = zero;
        for (std::size_t k : L.inchain[c])
            for (std::size_t i = 0; i < M; ++i) sumQ += sol.Q(i, k - 1);
        if (!(sumQ > zero)) continue;
        const T Nc = num_traits<T>::from_double(L.classes[c].population);
        for (std::size_t k : L.inchain[c])
            for (std::size_t i = 0; i < M; ++i) sol.Q(i, k - 1) = T(Nc * sol.Q(i, k - 1) / sumQ);
    }

    // An infinite server's utilization IS its queue length.
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].sched == SchedStrategy::INF)
            for (std::size_t r = 0; r < K; ++r) sol.U(i, r) = sol.Q(i, r);

    sol.C.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t i = 0; i < M; ++i) sol.C[r] += sol.R(i, r);
    zero_nans(sol.Q);
    zero_nans(sol.U);
    zero_nans(sol.R);
    zero_nans(sol.Tp);
    for (std::size_t r = 0; r < K; ++r)
        if (std::isnan(num_traits<T>::to_double(sol.C[r]))) sol.C[r] = zero;
    sol.iter = it_out;
    if (totiter) *totiter = static_cast<std::size_t>(it_out);
    return sol;
    }  // if constexpr has_transcendental
}

/**
 * Port of `solver_mam_basic_mmap.m`, the top-level dispatcher of the MMAP
 * fork-join decomposition: an open model goes straight to the inner algorithm
 * with the arrival rates its sources declare, a closed one through the
 * bisection wrapper.
 */
template <class T>
mva::MvaSolution<T> solver_mam_basic_mmap(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    const T zero = num_traits<T>::from_int(0);
    MmapDecConfig cfg;
    if (!L.is_open_model()) {
        std::size_t iter = 0;
        return solver_mam_basic_mmap_closed(L, opt, cfg, &iter);
    }
    // The chain's arrival rate is the total rate its classes are released at by
    // the reference station of its first class; every class of the chain carries
    // that same total, as the reference assigns it.
    std::vector<T> lambda(L.nclasses, zero);
    for (std::size_t c = 0; c < L.nchains; ++c) {
        if (L.inchain[c].empty()) continue;
        const std::size_t rs = L.classes[L.inchain[c][0] - 1].refstat - 1;
        T tot = zero;
        for (std::size_t k : L.inchain[c]) {
            if (L.disabled[rs][k - 1]) continue;
            const double v = num_traits<T>::to_double(L.rates(rs, k - 1));
            if (!std::isfinite(v)) continue;
            tot += L.rates(rs, k - 1);
        }
        for (std::size_t k : L.inchain[c]) lambda[k - 1] = tot;
    }
    std::size_t iter = 0;
    return solver_mam_basic_mmap_inner(L, opt, cfg, lambda, &iter);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_MMAP_H
