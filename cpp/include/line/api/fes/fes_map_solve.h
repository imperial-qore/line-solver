/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_SOLVE_H
#define LINE_API_FES_MAP_SOLVE_H

/**
 * Closed model left by a MAP flow-equivalent server, and the per-station metrics
 * behind it.
 *
 * Templated port of matlab/src/api/fes/fes_map_solve.m and
 * fes_map_deaggregate.m, mirrored by the JAR and native Python.
 *
 * Closes the aggregation of Section 5.2.1 of Casale, Mi, Cherkasova and Smirni,
 * IEEE Trans. Soft. Eng. 37(5), 2011. Once a subnetwork has been replaced by the
 * load-dependent MAP of `fes_map_aggregate`, the model left is a delay holding the
 * think times and one station, which is a finite level-dependent quasi birth-death
 * process: level k is the number of jobs held by the flow-equivalent server and
 * N-k jobs are thinking. The chain is the same block bidiagonal pair used to
 * measure the inter-departure times, now read as a generator rather than as a MAP,
 * so the delay is a station whose process is scaled by the number of jobs it holds
 * and the marked transitions are the arrivals into the flow-equivalent server. The
 * think time may itself be a MAP, which is how Section 5.3.1 models a bounded
 * flash crowd.
 *
 * The de-aggregation conditions on the population of the aggregate,
 * E[Y_i] = sum_k pk(k) Y_i(k), which is the decomposition step of Chandy, Herzog
 * and Woo, IBM J. Res. Dev. 19(1), 1975: exact for a product-form subnetwork, an
 * approximation when the burstiness the flow-equivalent server carries also
 * matters inside it. The aggregate metrics do not rely on it.
 *
 * ARITHMETIC: field operations only, exact at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/fes/fes_map_interdeparture.h"
#include "line/api/fes/fes_map_levels.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** System metrics of the delay plus flow-equivalent server model. */
template <class T>
struct FesMapSolveResult {
    T X;                ///< system throughput
    T R;                ///< response time of the aggregated subnetwork, N/X - E[Z]
    T Q;                ///< mean jobs held by the flow-equivalent server
    std::vector<T> pk;  ///< law of the jobs held by the flow-equivalent server
};

/** Per-station metrics an aggregate stands for. */
template <class T>
struct FesMapDeaggregateResult {
    std::vector<T> Q;
    std::vector<T> U;
    std::vector<T> X;
    std::vector<T> R;
};

/**
 * @param fes       flow-equivalent server, one MAP per level
 * @param think_map think time process (Z0,Z1)
 * @param n         number of jobs in the closed model
 */
template <class T>
FesMapSolveResult<T> fes_map_solve(const std::vector<mam::Map<T>>& fes,
                                   const mam::Map<T>& think_map, std::size_t n) {
    if (n < 1) throw InputError("fes_map_solve: the population must be at least 1");
    const std::vector<mam::Map<T>> fesLev = fes_map_levels(fes, n);
    const std::vector<mam::Map<T>> think =
        fes_map_levels(think_map, n, std::numeric_limits<double>::infinity());

    const mam::Map<T> T01 = fes_map_interdeparture(think, fesLev, n);
    const Matrix<T> Q = mam::map_infgen(T01);
    const std::size_t dim = Q.rows();
    const T zero = num_traits<T>::from_int(0);

    const std::vector<T> phi = mc::ctmc_solve(Q);
    const std::vector<T> phiT1 = vecmul(phi, T01.D1);

    FesMapSolveResult<T> out;
    out.X = zero;
    for (const T& v : phiT1) out.X += v;

    const std::size_t blk = dim / (n + 1);
    out.pk.assign(n + 1, zero);
    out.Q = zero;
    for (std::size_t k = 0; k <= n; ++k) {
        T s = zero;
        for (std::size_t j = 0; j < blk; ++j) s += phi[k * blk + j];
        out.pk[k] = s;
        out.Q += num_traits<T>::from_int(static_cast<long>(k)) * s;
    }
    out.R = num_traits<T>::from_int(static_cast<long>(n)) / out.X - mam::map_moment(think_map, 1);
    return out;
}

/**
 * @param pk       law of the jobs held by the aggregate, index k = P(k jobs)
 * @param L        service demands of the isolated subnetwork
 * @param mi       servers per station, ignored where isDelay
 * @param isDelay  true where the station is a pure delay
 */
template <class T>
FesMapDeaggregateResult<T> fes_map_deaggregate(const std::vector<T>& pk, const std::vector<T>& L,
                                               const std::vector<int>& mi,
                                               const std::vector<bool>& isDelay) {
    const std::size_t M = L.size();
    const std::size_t n = pk.size() - 1;
    const T zero = num_traits<T>::from_int(0);

    FesMapDeaggregateResult<T> out;
    out.Q.assign(M, zero);
    out.U.assign(M, zero);
    out.X.assign(M, zero);
    out.R.assign(M, zero);

    std::vector<std::size_t> queueIdx;
    T Z = zero;
    for (std::size_t i = 0; i < M; ++i) {
        if (isDelay[i]) {
            Z += L[i];
        } else {
            queueIdx.push_back(i);
        }
    }
    Matrix<T> Lq(queueIdx.size(), 1, zero);
    std::vector<int> miq(queueIdx.size(), 1);
    for (std::size_t j = 0; j < queueIdx.size(); ++j) {
        Lq(j, 0) = L[queueIdx[j]];
        miq[j] = mi[queueIdx[j]];
    }
    Matrix<T> Zm(1, 1, Z);

    for (std::size_t k = 1; k <= n; ++k) {
        if (pk[k] == zero) continue;
        const std::vector<int> N(1, static_cast<int>(k));
        const pfqn::MvaResult<T> r = pfqn::pfqn_mva(Lq, N, Zm, miq);
        const T Xk = r.XN[0];
        for (std::size_t j = 0; j < queueIdx.size(); ++j) {
            out.Q[queueIdx[j]] += pk[k] * r.QN(j, 0);
            out.U[queueIdx[j]] += pk[k] * r.UN(j, 0);
            out.X[queueIdx[j]] += pk[k] * Xk;
        }
        for (std::size_t i = 0; i < M; ++i) {
            if (isDelay[i]) {
                out.Q[i] += pk[k] * Xk * L[i];
                out.U[i] += pk[k] * Xk * L[i];
                out.X[i] += pk[k] * Xk;
            }
        }
    }

    for (std::size_t i = 0; i < M; ++i)
        if (out.X[i] != zero) out.R[i] = out.Q[i] / out.X[i];
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_SOLVE_H
