/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_FES_COMPUTE_METRICS_H
#define LINE_API_FES_FES_COMPUTE_METRICS_H

/**
 * Per-station metrics of the ISOLATED subnetwork at every population state, the
 * companion of `fes_compute_throughputs`.
 *
 * `fes_compute_throughputs` returns the aggregate throughput X(n) that becomes
 * the flow-equivalent server's rate. That is all the REDUCED model needs, but it
 * is not enough to report the COLLAPSED stations' own metrics: those are
 * recovered by conditioning on the FES population,
 *
 *   E[Q_i] = sum_n P(N_fes = n) * Q_i(n),
 *
 * the Chandy-Herzog-Woo hierarchical decomposition, which is EXACT when the
 * subnetwork is product-form. This header supplies the Q_i(n) and U_i(n) that
 * sum is taken over, on the same lattice and in the same linearized order as the
 * throughput table.
 *
 * Port of `matlab/src/api/fes/fes_compute_metrics.m`.
 */

#include <cstddef>
#include <vector>

#include "line/api/fes/ljd_linearize.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** The per-population tables the conditional sum is taken over. */
template <class T>
struct FesConditionalMetrics {
    /** Indexed by the 0-based linearized population state; each (M_sub x K). */
    std::vector<Matrix<T>> QN;
    /** Indexed by the 0-based linearized population state; each (M_sub x K). */
    std::vector<Matrix<T>> UN;
};

/**
 * @param L       (M_sub x K) service demands of the isolated subnetwork
 * @param mi      (M_sub) servers per station; a delay station's entry is ignored
 * @param isDelay (M_sub) true where the station is a pure delay
 * @param cutoffs (K) per-class population cutoffs
 */
template <class T>
FesConditionalMetrics<T> fes_compute_metrics(const Matrix<T>& L, const std::vector<int>& mi,
                                             const std::vector<bool>& isDelay,
                                             const std::vector<int>& cutoffs) {
    const std::size_t M_sub = L.rows();
    const std::size_t K = L.cols();
    if (isDelay.size() != M_sub)
        throw InputError("fes_compute_metrics: isDelay has the wrong length");
    if (cutoffs.size() != K)
        throw InputError("fes_compute_metrics: cutoffs and demands disagree on the class count");

    const T zero = num_traits<T>::from_int(0);

    std::vector<std::size_t> queueIdx, delayIdx;
    for (std::size_t i = 0; i < M_sub; ++i) (isDelay[i] ? delayIdx : queueIdx).push_back(i);
    const std::size_t M_queue = queueIdx.size();

    Matrix<T> L_queue(M_queue, K, zero);
    std::vector<int> mi_queue;
    for (std::size_t a = 0; a < M_queue; ++a) {
        for (std::size_t k = 0; k < K; ++k) L_queue(a, k) = L(queueIdx[a], k);
        mi_queue.push_back(mi.empty() ? 1 : mi[queueIdx[a]]);
    }
    Matrix<T> Z(1, K, zero);
    for (std::size_t d : delayIdx)
        for (std::size_t k = 0; k < K; ++k) Z(0, k) += L(d, k);

    std::size_t tableSize = 1;
    for (int c : cutoffs) {
        if (c < 0) throw InputError("fes_compute_metrics: negative cutoff");
        tableSize *= static_cast<std::size_t>(c + 1);
    }

    FesConditionalMetrics<T> out;
    out.QN.assign(tableSize, Matrix<T>(M_sub, K, zero));
    out.UN.assign(tableSize, Matrix<T>(M_sub, K, zero));

    for (std::size_t idx = 0; idx < tableSize; ++idx) {
        const std::vector<int> nvec = ljd_delinearize(idx + 1, cutoffs);
        int totalPop = 0;
        for (int v : nvec) totalPop += v;
        if (totalPop == 0) continue;  // an empty subnetwork holds nothing

        Matrix<T> Q(M_sub, K, zero), U(M_sub, K, zero);
        std::vector<T> XN(K, zero);
        if (M_queue > 0) {
            // `mi` is the additive C=L*(mi+Qarv) term of pfqn_mva, not a server
            // count: multiservers go through pfqn_mvams. Its UN is per STATION on
            // the closed multiserver branch and per station-class elsewhere, so it
            // is not read here -- utilization is recomputed analytically as
            // U=X*L/S, the [0,1] convention LINE uses at every queueing station
            // whatever its multiplicity.
            // See _kb/03-api-layer.md (pfqn_mva: mi is not S).
            const pfqn::MvaResult<T> r =
                pfqn::pfqn_mvams(std::vector<T>(K, zero), L_queue, nvec, Z,
                                 std::vector<int>(), mi_queue);
            XN = r.XN;
            for (std::size_t a = 0; a < M_queue; ++a) {
                const int srv = mi_queue.empty() || mi_queue[a] < 1 ? 1 : mi_queue[a];
                for (std::size_t k = 0; k < K; ++k) {
                    if (a < r.QN.rows() && k < r.QN.cols()) Q(queueIdx[a], k) = r.QN(a, k);
                    U(queueIdx[a], k) = T(XN[k] * L_queue(a, k) / num_traits<T>::from_int(srv));
                }
            }
        } else {
            for (std::size_t k = 0; k < K; ++k)
                if (nvec[k] > 0 && Z(0, k) > zero)
                    XN[k] = num_traits<T>::from_int(nvec[k]) / Z(0, k);
        }
        // A delay holds X_k * Z_i(k) jobs by Little's law on a station with no
        // queueing, and its INF "utilization" is that same population.
        for (std::size_t d : delayIdx)
            for (std::size_t k = 0; k < K; ++k) {
                Q(d, k) = T(XN[k] * L(d, k));
                U(d, k) = Q(d, k);
            }
        out.QN[idx] = Q;
        out.UN[idx] = U;
    }
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_FES_COMPUTE_METRICS_H
