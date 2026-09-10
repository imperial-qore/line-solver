/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_DECMMAP_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_DECMMAP_H

/**
 * Port of `solver_mam.m`, the `dec.mmap` method: the per-class departure-process
 * decomposition.
 *
 * IT IS THE MMAP DECOMPOSITION WITHOUT THE FORK-JOIN MACHINERY. The sweep is the
 * same one `solver_mam_basic_mmap_inner` runs -- rebuild every arrival stream
 * from the departure table through the traffic equations, solve each station in
 * isolation, and replace each departure process by the ETAQA truncation of that
 * station's own QBD -- with three differences, all of them the reference's:
 *   - the traffic step is `solver_mam_traffic`, the plain one, so DEP is
 *     STATION-indexed and a Fork or a Join in the model is refused by the
 *     traffic step itself rather than synchronized;
 *   - there is no Join post-processing, because there is no Join;
 *   - the increment norm is `max(|xn-xr|./xr)`, WITHOUT the FineTol offset the
 *     fork-join variant adds, so the first sweep's norm is NaN by construction
 *     (Q starts at zero) and `config.da_miniter = 3` is what carries the loop
 *     past it.
 *
 * WHAT THE REFERENCE REFUSES, AND HOW THIS PORT REFUSES IT. Two model classes
 * get NO ANSWER from `solver_mam.m`, and neither refusal is an error there:
 *   - a station whose discipline is not EXT, FCFS, HOL, FCFSPRPRIO or PS makes
 *     it return `[]` for all six metrics with `method = ''`, after a warning.
 *     Note that INF is NOT in that list: a model with a Delay is refused;
 *   - a model that is not purely open returns the ZERO matrices it initialised,
 *     again after a warning.
 * Returning empty or all-zero metrics under a method name is exactly the
 * silently-wrong answer this port refuses to produce, so both are thrown BY NAME
 * here, with the reference's own reason in the message.
 *
 * THE TRAILING SURROGATE-DELAY BLOCK IS OUTSIDE THE NODE SWITCH in the
 * reference, so it runs at EVERY station and not only at the queues -- including
 * the Source, whose `PH` is its own ARRIVAL process (the EXT branch never
 * assigns one) and whose utilization therefore comes out as lambda * (1/lambda)
 * = 1. Reproduced as written: it is what the reference reports, the surrogate
 * term is zero at the Source because it has one server, and `getAvg`'s metric
 * filter is what decides whether a Source utilization is shown.
 *
 * ARITHMETIC. Double (or real) only, for the reasons `solver_mam_basic.h` lists.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/mam/qbd_depproc.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mam/solver_mam_basic_mmap.h"
#include "line/solvers/mam/solver_mam_bmap.h"  // mam_detect_mmck
#include "line/solvers/mam/solver_mam_traffic.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace decmmap_detail {

using lang::GlobalConstants;
using lang::SchedStrategy;

/** `solver_mam.m`'s opening station loop: the disciplines it has a branch for. */
template <class T>
void check_disciplines(const qn::NetworkStruct<T>& L) {
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        if (sc == SchedStrategy::EXT || sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL ||
            sc == SchedStrategy::FCFSPRPRIO || sc == SchedStrategy::PS)
            continue;
        throw UnsupportedError(
            std::string("SolverMAM: the dec.mmap method has no branch for ") +
            lang::sched_to_text(sc) + " scheduling at station '" + L.stations[i].name +
            "'. solver_mam.m warns and returns EMPTY metrics with a cleared method name for such a "
            "model, which is no answer at all; note that INF is not in its list either, so a model "
            "carrying a Delay belongs to 'dec.source' or 'dec.source.mmap'");
    }
}

}  // namespace decmmap_detail

/**
 * Port of `solver_mam.m`.
 *
 * @param L   the refreshed struct; open classes only
 * @param opt the SolverMAM options; `tol` is the fixed point's iter_tol
 */
template <class T>
mva::MvaSolution<T> solver_mam_decmmap(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L; (void)opt;
        throw UnsupportedError(
            "solver_mam_decmmap: the departure-process fixed point stops on a tolerance, "
            "MMAP[K]/PH[K]/1 runs the ADDA doubling iteration and the ETAQA departure process "
            "needs transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
    using namespace decmmap_detail;
    using namespace basic_mmap_detail;
    using basic_detail::station_visits;
    using basic_detail::truncate_renorm;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;
    const T ftol = num_traits<T>::from_double(GlobalConstants::FineTol);
    MmapDecConfig cfg;

    check_disciplines(L);
    if (!L.is_open_model())
        throw UnsupportedError(
            "SolverMAM: the dec.mmap method solves open models only. solver_mam.m warns and "
            "returns the ZERO matrices it initialised for anything else, which reports an empty "
            "network rather than refusing; use 'dec.source', 'dec.source.mmap' or 'mna'");

    // The chain arrival rate, as solver_mam.m computes it: the total rate the
    // chain's classes are released at by the reference station of its first
    // class, carried by every class of the chain.
    std::vector<T> lambda(K, zero);
    for (std::size_t c = 0; c < C; ++c) {
        if (L.inchain[c].empty()) continue;
        const std::size_t rs = L.classes[L.inchain[c][0] - 1].refstat - 1;
        T tot = zero;
        for (std::size_t k : L.inchain[c]) {
            if (L.disabled[rs][k - 1]) continue;
            if (!std::isfinite(num_traits<T>::to_double(L.rates(rs, k - 1)))) continue;
            tot += L.rates(rs, k - 1);
        }
        for (std::size_t k : L.inchain[c]) lambda[k - 1] = tot;
    }

    const Matrix<T> V = station_visits(L);
    const MmapPhTable<T> ph = mmap_ph_table(L);
    const TrafficConfig tcfg = traffic_config(opt);

    Matrix<T> QN(M, K, zero), UN(M, K, zero), RN(M, K, zero), TN(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].sched == SchedStrategy::EXT)
            for (std::size_t r = 0; r < K; ++r)
                TN(i, r) = L.disabled[i][r] ? zero : L.rates(i, r);

    // STATION-indexed here, node-indexed in the fork-join variant: the plain
    // traffic step reads DEP[station][class].
    DepTable<T> DEP(M, std::vector<Map<T>>(K));

    auto sweep = [&](const std::vector<T>&,
                     std::size_t itnum) -> std::pair<std::vector<T>, std::vector<T>> {
        if (itnum == 1)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r)
                    DEP[i][r] = (lambda[r] > zero && V(i, r) > zero)
                                    ? map_scale(ph.PH[i][r], T(one / T(lambda[r] * V(i, r))))
                                    // The reference divides unconditionally; a
                                    // zero denominator would scale to an
                                    // infinite mean, which the traffic step
                                    // reads as a silent stream. Left as the
                                    // unscaled process instead, which is what
                                    // that stream carries.
                                    : ph.PH[i][r];

        std::vector<Mmap<T>> ARV = solver_mam_traffic(L, DEP, tcfg);

        std::vector<T> xref(M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) xref[i * K + r] = QN(i, r);

        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t ind = L.node_of_station(i + 1) - 1;
            const SchedStrategy sc = L.stations[i].sched;
            const double ns = L.stations[i].nservers;
            bool finiteCapUsed = false;

            if (L.nodes[ind].nodetype == qn::NodeType::Queue && ARV[ind].order() > 0) {
                if (ARV[ind].order() > tcfg.space_max) ARV[ind] = compress_arrival(ARV[ind]);
                const std::vector<T> lam = mmap_lambda(ARV[ind]);
                for (std::size_t r = 0; r < K; ++r) TN(i, r) = lam[r];

                std::vector<PhService<T>> sl;
                for (std::size_t r = 0; r < K; ++r) sl.push_back(ph.svc[i][r]);

                if (sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL ||
                    sc == SchedStrategy::FCFSPRPRIO) {
                    if (std::isfinite(L.stations[i].cap)) {
                        const std::size_t capK =
                            static_cast<std::size_t>(std::llround(L.stations[i].cap));
                        T meanQ = zero, lossProb = zero;
                        const MmckDetection<T> det = mam_detect_mmck(L, i + 1, ARV[ind]);
                        if (det.isMmck) {
                            T lamTot = zero;
                            for (std::size_t r = 0; r < K; ++r)
                                if (!std::isnan(num_traits<T>::to_double(TN(i, r))))
                                    lamTot += TN(i, r);
                            const qsys::MmckResult<T> ex = qsys::qsys_mmck(
                                lamTot, det.muRate, static_cast<unsigned>(std::llround(ns)),
                                static_cast<unsigned>(capK));
                            meanQ = ex.meanQueueLength;
                            lossProb = ex.lossProbability;
                        } else {
                            const basic_detail::TruncRenorm<T> tr =
                                truncate_renorm(ARV[ind], sl, capK);
                            meanQ = tr.meanQ;
                            lossProb = tr.lossProb;
                        }
                        std::vector<T> eff(K, zero), Sact(K, zero);
                        T sumTN = zero;
                        for (std::size_t r = 0; r < K; ++r) {
                            const T inflow =
                                std::isnan(num_traits<T>::to_double(TN(i, r))) ? zero : TN(i, r);
                            eff[r] = T(inflow * T(one - lossProb));
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
                        // No saturation test and no MAP/MAP/1 fast path here:
                        // solver_mam.m goes straight to MMAP[K]/PH[K]/1, and its
                        // fork-join sibling is the one that added both.
                        const std::vector<T> m = mmapph1fcfs_ncmean(ARV[ind], sl);
                        for (std::size_t r = 0; r < K; ++r)
                            QN(i, r) = m[ARV[ind].classes() == 1 ? 0 : r];
                    }
                } else if (sc == SchedStrategy::PS) {
                    for (std::size_t r = 0; r < K; ++r)
                        UN(i, r) = T(TN(i, r) * map_mean(ph.PH[i][r]));
                    T usum = zero;
                    for (std::size_t r = 0; r < K; ++r) usum += UN(i, r);
                    const T uden = (usum < T(one - ftol)) ? usum : T(one - ftol);
                    for (std::size_t r = 0; r < K; ++r) QN(i, r) = T(UN(i, r) / T(one - uden));
                }
            }

            // OUTSIDE the node switch in the reference; see the file header.
            if (!finiteCapUsed) {
                for (std::size_t r = 0; r < K; ++r) {
                    UN(i, r) = T(TN(i, r) * map_mean(ph.PH[i][r]));
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
            if (L.nodes[ind].nodetype != qn::NodeType::Queue || ARV[ind].order() == 0) continue;
            const SchedStrategy sc = L.stations[i].sched;
            const bool fcfs = (sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL ||
                               sc == SchedStrategy::FCFSPRPRIO);
            if (!fcfs && sc != SchedStrategy::PS) continue;
            T rho = zero;
            for (std::size_t r = 0; r < K; ++r) rho += UN(i, r);
            for (std::size_t r = 0; r < K; ++r) {
                const Mmap<T> A = mmap_hide_but(ARV[ind], r);
                const Map<T>& Srv = ph.PH[i][r];
                const std::size_t etaqa_sz = (cfg.etaqa_trunc + 1) * A.order() * Srv.D0.rows();
                Map<T> dep = Srv;
                if (etaqa_sz <= tcfg.space_max && rho < T(one - ftol)) {
                    // The reference's own try/catch, as in the fork-join variant.
                    try {
                        const Map<T> Am{A.D0, A.D1};
                        dep = map_normalize(fcfs ? qbd_depproc_etaqa(Am, Srv, cfg.etaqa_trunc)
                                                 : qbd_depproc_etaqa_ps(Am, Srv, cfg.etaqa_trunc));
                    } catch (const Error&) {
                        dep = Srv;
                    }
                }
                // Unconditional in the reference, which divides by zero when the
                // class does not flow here; left unscaled instead, as at the
                // initialisation above.
                if (lambda[r] > zero && V(i, r) > zero)
                    dep = map_scale(dep, T(one / T(lambda[r] * V(i, r))));
                DEP[i][r] = dep;
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
    fo.miniter = 3;             // config.da_miniter
    fo.relative_norm = true;    // config.da_norm, WITHOUT the FineTol offset
    const da::FpiResult<T> fr = da::da_fpi<T>(sweep, std::vector<T>(M * K, zero), fo);

    mva::MvaSolution<T> out;
    out.Q = QN;
    out.U = UN;
    out.R = RN;
    out.Tp = TN;
    // CN and XN are the zeros solver_mam.m initialises and never assigns.
    out.C.assign(K, zero);
    out.X.assign(K, zero);
    out.method = "dec.mmap";
    out.iter = static_cast<int>(fr.iterations);
    out.lG = 0.0;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_DECMMAP_H
