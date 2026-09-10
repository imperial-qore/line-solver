/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MMAPG1K_H
#define LINE_API_QSYS_QSYS_MMAPG1K_H

/**
 * Exact per-class throughput and loss ratio of an MMAP[K]/G/1/K tail-drop
 * queue. Port of matlab/src/api/qsys/qsys_mmapg1k.m.
 *
 * Two classes of equal arrival rate but different interarrival variability or
 * autocorrelation receive different loss ratios. An aggregate-only finite
 * buffer analysis cannot express that: it returns one blocking probability p
 * and sets T_k = lambda_k (1 - p), making the loss ratio identical across
 * classes BY CONSTRUCTION. What breaks the tie here is the phase resolution of
 * the full-buffer probability. A class-k arrival leaves phase i at rate
 * (D1c_k e)_i, so
 *
 *     lambda_k = pi D1c_k e,    L_k = (pKvec D1c_k e)/lambda_k,
 *     T_k      = lambda_k (1 - L_k),
 *
 * with pi the stationary phase law of the aggregate MAP and pKvec the joint law
 * of (level = K, phase) that qsys_mapg1k returns. No independence between
 * classes is assumed and no PASTA argument is used.
 *
 * This is exact whenever the joint MMAP is available, which it is inside a
 * solver that propagates MMAPs between stations. Use qsys_mapg1k_perflow for
 * the setting where the flows are given as separate MAPs and the joint process
 * would cost prod_n M_n phases.
 *
 * The model assumes a single server and a service law that is iid and
 * independent of class: per-class service would make the departure rate depend
 * on which class holds the server, which this state space does not represent.
 *
 * ARITHMETIC. Inherits the transcendental gate from qsys_mapg1k; the per-class
 * layer on top is finite exact linear algebra.
 *
 * Reference: Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied
 * System Innovation 2026, 9, 112.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapg1k.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mmapg1k, mirroring the MATLAB result struct. */
template <class T>
struct MmapG1kResult {
    std::vector<T> throughput;  ///< per-class throughput
    std::vector<T> lossRatio;   ///< per-class loss ratio, in [0,1]
    std::vector<T> lambda;      ///< per-class arrival rate
    T lambdaAggregate;
    T throughputAggregate;
    T lossAggregate;            ///< aggregate loss ratio of the whole stream
    T p0;
    T pK;
    std::vector<T> pKvec;
    std::vector<T> plevel;
    T meanQueueLength;
    T meanServiceTime;
    T utilization;
    T rho;
};

/**
 * MMAP[K]/G/1/K with tail drop.
 *
 * @param D0  hidden transition matrix of the arrival MMAP (M x M)
 * @param D1c per-class arrival matrices; D0 + sum_k D1c[k] must be an
 *            irreducible generator
 * @param svc service law, shared by all classes
 * @param K   buffer size in packets
 * @param tol convergence tolerance
 * @param nmaxCap cap on the level truncation
 */
template <class T>
MmapG1kResult<T> qsys_mmapg1k(const Matrix<T>& D0, const std::vector<Matrix<T>>& D1c,
                              const ServiceLaw<T>& svc, std::size_t K, const T& tol,
                              std::size_t nmaxCap) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mmapg1k requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (D1c.empty()) throw InputError("qsys_mmapg1k: at least one class is required");
    const std::size_t M = D0.rows();
    if (D0.cols() != M) throw InputError("qsys_mmapg1k: D0 must be square");
    const std::size_t R = D1c.size();
    mam::Map<T> agg;
    agg.D0 = D0;
    agg.D1 = Matrix<T>(M, M, zero);
    for (std::size_t k = 0; k < R; ++k) {
        if (D1c[k].rows() != M || D1c[k].cols() != M)
            throw InputError("qsys_mmapg1k: every per-class D1 must match the order of D0");
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) agg.D1(i, j) += D1c[k](i, j);
    }

    const MapG1kResult<T> r = qsys_mapg1k(agg, svc, K, tol, nmaxCap);
    const std::vector<T> e = ones<T>(M);
    const std::vector<T> pit = mam::map_prob(agg);

    MmapG1kResult<T> out;
    out.throughput.assign(R, zero);
    out.lossRatio.assign(R, zero);
    out.lambda.assign(R, zero);
    out.lambdaAggregate = zero;
    out.throughputAggregate = zero;
    for (std::size_t k = 0; k < R; ++k) {
        const std::vector<T> col = mulvec(D1c[k], e);
        T lam = zero, blocked = zero;
        for (std::size_t i = 0; i < M; ++i) {
            lam += pit[i] * col[i];
            blocked += r.pKvec[i] * col[i];
        }
        out.lambda[k] = lam;
        out.lossRatio[k] = (lam > zero) ? T(blocked / lam) : zero;
        out.throughput[k] = lam * (one - out.lossRatio[k]);
        out.lambdaAggregate += lam;
        out.throughputAggregate += out.throughput[k];
    }
    out.lossAggregate = r.lossProbability;
    out.p0 = r.p0;
    out.pK = r.pK;
    out.pKvec = r.pKvec;
    out.plevel = r.plevel;
    out.meanQueueLength = r.meanQueueLength;
    out.meanServiceTime = r.meanServiceTime;
    out.utilization = r.utilization;
    out.rho = r.rho;
    return out;
}

/** qsys_mmapg1k with the qsys_mapg1k defaults tol = 1e-12, nmax = 200000. */
template <class T>
MmapG1kResult<T> qsys_mmapg1k(const Matrix<T>& D0, const std::vector<Matrix<T>>& D1c,
                              const ServiceLaw<T>& svc, std::size_t K) {
    return qsys_mmapg1k(D0, D1c, svc, K, T(num_traits<T>::from_double(1e-12)),
                        static_cast<std::size_t>(200000));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MMAPG1K_H
