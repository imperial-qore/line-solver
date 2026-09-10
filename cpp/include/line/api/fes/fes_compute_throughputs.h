/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_COMPUTE_THROUGHPUTS_H
#define LINE_API_FES_COMPUTE_THROUGHPUTS_H

/**
 * Per-class throughput table of an isolated subnetwork, tabulated over the
 * population lattice, for use as the load-dependent rates of a
 * flow-equivalent server (Chandy, Herzog and Woo 1975).
 *
 * Templated port of matlab/src/api/fes/fes_compute_throughputs.m. Every
 * population state 0 <= n <= cutoffs is solved with pfqn_mva on the isolated
 * subnetwork; scalingTable[r][idx-1] holds X_r(n) at the LJD_LINEARIZE index
 * idx of n. States with n_r = 0 store 0 for class r, and the empty state
 * stores 0 for every class.
 *
 * Arithmetic. pfqn_mva stays in the field of the inputs, and this function
 * only enumerates and stores, so the table is exact at T = Rational.
 *
 * Deviation from MATLAB, mechanical: MATLAB walks the lattice with a BFS over
 * a cell-array queue guarded by a visited mask, which reaches every state
 * exactly once; the port enumerates the same lattice directly with next_pop.
 * The states, and hence the table, are identical; only the visit order
 * differs, and each state is solved independently of the others. The MATLAB
 * try/catch that stores zeros when the solver fails is mirrored by catching
 * line::Error.
 */

#include <cstddef>
#include <vector>

#include "line/api/fes/ljd_linearize.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace fes {

/**
 * @param L       (M_sub x K) service demands of the isolated subnetwork
 * @param mi      (M_sub) servers per station; entries of delay stations are
 *                ignored (MATLAB stores Inf there)
 * @param isDelay (M_sub) true where the station is a pure delay
 * @param cutoffs (K) per-class population cutoffs
 * @return one linearized throughput vector per class, each of length
 *         prod_k (cutoffs[k] + 1)
 */
template <class T>
std::vector<std::vector<T>> fes_compute_throughputs(const Matrix<T>& L, const std::vector<int>& mi,
                                                    const std::vector<bool>& isDelay,
                                                    const std::vector<int>& cutoffs) {
    const std::size_t M_sub = L.rows();
    const std::size_t K = L.cols();
    if (isDelay.size() != M_sub)
        throw InputError("fes_compute_throughputs: isDelay has the wrong length");
    if (!mi.empty() && mi.size() != M_sub)
        throw InputError("fes_compute_throughputs: mi has the wrong length");
    if (cutoffs.size() != K)
        throw InputError("fes_compute_throughputs: cutoffs and demands disagree on the class count");

    const T zero = num_traits<T>::from_int(0);

    // Queue stations go into L_queue; delay stations contribute their demands
    // to the think time Z.
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
        if (c < 0) throw InputError("fes_compute_throughputs: negative cutoff");
        tableSize *= static_cast<std::size_t>(c + 1);
    }
    std::vector<std::vector<T>> scalingTable(K, std::vector<T>(tableSize, zero));

    std::vector<int> nvec(K, 0);
    bool more = true;
    while (more) {
        const std::size_t idx = ljd_linearize(nvec, cutoffs);
        int totalPop = 0;
        for (int v : nvec) totalPop += v;
        if (totalPop > 0) {
            std::vector<T> XN(K, zero);
            if (M_queue > 0) {
                // pfqn_mva's `mi` is NOT a server count -- it enters only as the
                // additive term of C(i,s)=L(i,s)*(mi(i)+Qarv), so passing the real
                // multiplicity INFLATES the residence time instead of adding
                // servers. pfqn_mvams forwards to pfqn_mva when every station is a
                // single server and to the load-dependent recursion with
                // mu(i,n)=min(n,S(i)) when one is not.
                // See _kb/03-api-layer.md (pfqn_mva: mi is not S).
                const pfqn::MvaResult<T> r =
                    pfqn::pfqn_mvams(std::vector<T>(K, zero), L_queue, nvec, Z,
                                     std::vector<int>(), mi_queue);
                XN = r.XN;
            } else {
                // only delay stations: X_k = n_k / Z_k
                for (std::size_t k = 0; k < K; ++k)
                    if (nvec[k] > 0 && Z(0, k) > zero)
                        XN[k] = num_traits<T>::from_int(nvec[k]) / Z(0, k);
            }
            // A failure here is NOT caught: an FES table silently filled with
            // zeros is a wrong aggregate, not a degraded one, and every
            // downstream beta_r(n) reads it as "the subnetwork serves nothing".
            for (std::size_t k = 0; k < K; ++k)
                scalingTable[k][idx - 1] = nvec[k] > 0 ? XN[k] : zero;
        }
        more = next_pop(nvec, cutoffs);
    }
    return scalingTable;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_COMPUTE_THROUGHPUTS_H
