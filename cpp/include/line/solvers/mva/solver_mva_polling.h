/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_POLLING_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_POLLING_H

/**
 * The multiclass open polling analyzer (ladder branch 5).
 *
 * Templated port of `matlab/src/solvers/MVA/solver_mva_polling_analyzer.m`.
 *
 * A polling system is a single server that visits each class's buffer in turn,
 * switching over between them. Each CLASS is a queue of the polling model, and
 * the mean waiting time of every class follows from the pseudo-conservation
 * laws of its discipline: EXHAUSTIVE (serve until the buffer empties), GATED
 * (serve exactly what was present at the visit), 1-limited (one job per visit)
 * and DECREMENTING. The `polling_qsys_*` routines carry those laws; this
 * analyzer builds their two-moment inputs from the arrival, service and
 * switchover distributions and maps the resulting waiting times back onto the
 * per-class [Q,U,R,T] block.
 *
 * SUPPORTED. EXHAUSTIVE and GATED at any number of servers-of-one; K-limited
 * only at K = 1 (the `polling_qsys_1limited` law); DECREMENTING. The reference
 * refuses K > 1 and any other discipline, and so does this port. `exact` is
 * available only for EXHAUSTIVE / GATED with Poisson arrivals at a single
 * server, where the two-moment law is exact.
 *
 * Arithmetic: TRANSCENDENTAL. The pseudo-conservation laws are field
 * arithmetic, but the input SCVs come from distribution moments that are only
 * defined at a transcendental T.
 */

#include <cmath>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/polling/polling_qsys_1limited.h"
#include "line/api/polling/polling_qsys_exhaustive.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mva {

/** Port of `solver_mva_polling_analyzer.m`. */
template <class T>
MvaSolution<T> solver_mva_polling_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::NodeType;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        throw UnsupportedError(
            "solver_mva_polling_analyzer: the two-moment polling inputs come from distribution "
            "SCVs, which need a transcendental T");
    } else {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;

    std::size_t src = 0, q = 0;
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].nodetype == NodeType::Source) src = i + 1;
        if (L.stations[i].nodetype == NodeType::Queue) q = i + 1;
    }
    if (src == 0 || q == 0)
        throw InputError("solver_mva_polling_analyzer: a polling system needs a Source and a Queue");
    const qn::Station<T>& qst = L.stations[q - 1];
    if (qst.polling_type.empty())
        throw InputError("solver_mva_polling_analyzer: the Queue has no polling parameters");
    const std::size_t qstateful = L.stateful_of_station(q);
    // the queue's visit ratio per class, summed over chains: each open class
    // sits in its own chain, so reading only chain 0 would zero every class but
    // the first
    std::vector<T> Vq(K, zero);
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t r = 0; r < K; ++r) Vq[r] += L.visits[c](qstateful - 1, r);

    // per-class arrival rate at the queue, service moments and switchover
    // moments -- the PollingMoments the qsys routines consume
    polling::PollingMoments<T> pm;
    pm.lambda.assign(K, zero);
    pm.b.assign(K, zero);
    pm.b2.assign(K, zero);
    pm.r.assign(K, zero);
    pm.delta2.assign(K, zero);
    std::vector<T> mu(K, zero), ca2(K, one);
    for (std::size_t r = 0; r < K; ++r) {
        pm.lambda[r] = T(L.rates(src - 1, r) * Vq[r]);
        mu[r] = L.rates(q - 1, r);
        const T b = mu[r] > zero ? T(one / mu[r]) : zero;
        pm.b[r] = b;
        // E[S^2] OFF THE PROCESS REPRESENTATION, not off the nominal SCV. The
        // reference hands the qsys routines `sn.proc{queue_ist}`, and
        // `refreshProcessRepresentations` has already replaced a non-Markovian
        // service by `convertToMAP` -- an Erlang matched to the SCV, 20 phases
        // when the law is deterministic. Reading `E[S]^2 (1 + scv)` instead gave
        // a Det service scv = 0 exactly, i.e. b2 = E[S]^2, where the reference
        // carries the Erlang-20 value E[S]^2 (1 + 1/20); on
        // polling_exhaustive_det that is a 5th-significant-digit gap in W.
        if (!L.service[q - 1][r].disabled && mu[r] > zero) {
            pm.b2[r] = mam::map_moment(lang::dist_to_map(L.service[q - 1][r]), 2);
        } else {
            pm.b2[r] = T(b * b * T(one + L.scv(q - 1, r)));
        }
        ca2[r] = L.scv(src - 1, r);
        // AN ABSENT SWITCHOVER IS `Immediate()`, NOT ZERO. Both references
        // default the cell that way -- `solver_mva_polling_analyzer.m` fills it
        // with `Immediate()` and the Python twin appends one -- and the
        // distinction is the whole answer: the station-time formula divides by
        // the total walk time R, so leaving an unset switchover at zero made R
        // zero and `polling_qsys_exhaustive` refuse a model the reference
        // solves. `polling_exhaustive_exp` sets no switchover at all and is
        // goldened for finite waiting times.
        const bool has_so = r < qst.switchover.size() && !qst.switchover[r].disabled;
        // The reference feeds the qsys routines sn.proc, i.e. the PROCESS
        // representation, and Immediate.getProcess() is the exponential of rate
        // GlobalConstants.Immediate rather than a point mass at zero. Reading
        // `so.mean` instead made an immediate switchover a zero total walk time,
        // which the qsys routines refuse for the same reason.
        const bool immediate = !has_so || qst.switchover[r].is_immediate();
        const T so_mean = immediate ? T(num_traits<T>::from_int(1) /
                                        num_traits<T>::from_double(GlobalConstants::Immediate))
                                    : qst.switchover[r].mean;
        const T so_scv = immediate ? one : qst.switchover[r].scv;
        pm.r[r] = so_mean;
        pm.delta2[r] = T(so_scv * so_mean * so_mean);
    }

    // the discipline, assumed identical across buffers as the reference does
    const lang::PollingType pt = qst.polling_type[0];
    std::string method = opt.method;
    if (method == "exact") {
        bool all_poisson = true;
        for (std::size_t r = 0; r < K; ++r)
            if (std::fabs(num_traits<T>::to_double(ca2[r]) - 1.0) >= 1e-9) all_poisson = false;
        const bool single_server =
            !std::isfinite(qst.nservers) ? false : std::llround(qst.nservers) == 1;
        if (!(all_poisson && single_server &&
              (pt == lang::PollingType::EXHAUSTIVE || pt == lang::PollingType::GATED)))
            throw UnsupportedError(
                "solver_mva_polling_analyzer: the exact method covers only EXHAUSTIVE / GATED "
                "polling with Poisson arrivals at a single server");
        method = "stationtime";
    }

    std::vector<T> W;
    switch (pt) {
        case lang::PollingType::EXHAUSTIVE:
            W = polling::polling_qsys_exhaustive(pm);
            break;
        case lang::PollingType::GATED:
            W = polling::polling_qsys_gated(pm);
            break;
        case lang::PollingType::KLIMITED:
            if (qst.polling_par != 1)
                throw UnsupportedError(
                    "solver_mva_polling_analyzer: K-limited polling is only available at K = 1");
            W = polling::polling_qsys_1limited(pm);
            break;
        case lang::PollingType::DECREMENTING:
            W = polling::polling_qsys_decrementing(pm);
            break;
        default:
            throw UnsupportedError("solver_mva_polling_analyzer: unsupported polling type");
    }

    // R = W + 1/mu, then scale by the queue's visit ratio, as the reference does
    const int k = std::isfinite(qst.nservers) ? static_cast<int>(std::llround(qst.nservers)) : 1;
    MvaSolution<T> out;
    out.Q = Matrix<T>(M, K, zero);
    out.U = Matrix<T>(M, K, zero);
    out.R = Matrix<T>(M, K, zero);
    out.Tp = Matrix<T>(M, K, zero);
    out.C.assign(K, zero);
    out.X.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        const T Rr = T(T(W[r] + (mu[r] > zero ? T(one / mu[r]) : zero)) * Vq[r]);
        out.R(q - 1, r) = Rr;
        out.C[r] = Rr;
        out.Tp(src - 1, r) = pm.lambda[r];
        out.Tp(q - 1, r) = pm.lambda[r];
        out.X[r] = pm.lambda[r];
        out.U(q - 1, r) = mu[r] > zero
                              ? T(pm.lambda[r] / mu[r] / num_traits<T>::from_int(k))
                              : zero;
        out.Q(q - 1, r) = T(pm.lambda[r] * Rr);
    }
    out.method = "stationtime";
    out.lG = 0.0;
    (void)method;
    return out;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_POLLING_H
