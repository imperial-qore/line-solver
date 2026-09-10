/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPG1K_PERFLOW_H
#define LINE_API_QSYS_QSYS_MAPG1K_PERFLOW_H

/**
 * Per-flow throughput and loss ratio of a tail-drop FIFO buffer fed by N flows
 * of mutually different statistical character. Port of
 * matlab/src/api/qsys/qsys_mapg1k_perflow.m, Theorem 1 of [1].
 *
 * Flow n has its own MAP, so two flows may share an arrival rate and still
 * differ in the shape and the autocorrelation of their interarrival times.
 *
 * METHOD, AND WHY IT IS AN APPROXIMATION. The exact model would track the
 * modulating state of every flow jointly with the buffer, prod_n M_n (K+1)
 * states, which [1] notes is already out of reach at N = 10, M_n = 3, K = 10.
 * Instead ONE model is solved per flow: flow n is kept exactly and the other
 * N-1 flows are replaced by a single Poisson stream of rate lambda - lambda_n,
 * which the Palm-Khinchin limiting theorem on superposed point processes
 * justifies as N grows. The substitution is applied N times, once per flow, so
 * no flow is ever the Poissonized one when its own throughput is computed.
 * Superposing MAP_n with the Poisson background gives ([1], eq. 5)
 *
 *     D0 = D0n - lambdaBar I,     D1 = D1n + lambdaBar I,
 *
 * and the flow throughput is read off the aggregate ([1], eq. 20)
 *
 *     T_n = (1 - p0)/S + pK lambdaBar - lambdaBar,
 *
 * the aggregate departure rate less the background throughput
 * lambdaBar (1 - pK), the background loss ratio being pK by PASTA BECAUSE the
 * background is Poisson. PASTA is used for the background only, never for the
 * flow being measured. Cost is O(N (K M)^3) against O(M^(3N) K^3): linear
 * rather than exponential in the flow count.
 *
 * ACCURACY. [1] reports errors against simulation of the exact model below
 * about 8% for N >= 9 with K >= 20, falling to 2.1% at K = 50, 0.5% at K = 100
 * and 1.2% at N = 900. Errors are largest when flows are few, highly variable
 * and the buffer is small. This is a property of the method, not of the port:
 * the port reproduces the reference's numbers, and the two agree with the
 * published Table 2 of [1] to the precision at which it is printed.
 *
 * ARITHMETIC. Inherits the transcendental gate from qsys_mapg1k.
 *
 * Reference: [1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied
 * System Innovation 2026, 9, 112.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapg1k.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mapg1k_perflow, mirroring the MATLAB result struct. */
template <class T>
struct MapG1kPerflowResult {
    std::vector<T> throughput;  ///< per-flow throughput
    std::vector<T> lossRatio;   ///< per-flow loss ratio
    std::vector<T> lambda;      ///< per-flow arrival rate
    T lambdaAggregate;
    T throughputAggregate;
    T lossAggregate;            ///< sum_n L_n lambda_n / lambda
    std::vector<T> p0;          ///< empty-buffer probability of the n-th model
    std::vector<T> pK;          ///< full-buffer probability of the n-th model
    T meanServiceTime;
    T rho;                      ///< offered load lambdaAggregate * S
};

/**
 * Per-flow analysis of a MAP-fed tail-drop buffer.
 *
 * @param flows the N arrival MAPs, whose modulating orders may differ
 * @param svc   service law, shared by all flows
 * @param K     buffer size in packets, the one in transmission included
 * @param tol convergence tolerance
 * @param nmaxCap cap on the level truncation
 */
template <class T>
MapG1kPerflowResult<T> qsys_mapg1k_perflow(const std::vector<mam::Map<T>>& flows,
                                           const ServiceLaw<T>& svc, std::size_t K, const T& tol,
                                           std::size_t nmaxCap) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapg1k_perflow requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t N = flows.size();
    if (N == 0) throw InputError("qsys_mapg1k_perflow: at least one flow is required");

    std::vector<T> lam(N);
    T lamTot = zero;
    for (std::size_t n = 0; n < N; ++n) {
        lam[n] = mam::map_lambda(flows[n]);
        if (lam[n] <= zero)
            throw InputError("qsys_mapg1k_perflow: every flow must have a positive arrival rate");
        lamTot += lam[n];
    }

    MapG1kPerflowResult<T> out;
    out.throughput.assign(N, zero);
    out.lossRatio.assign(N, zero);
    out.lambda = lam;
    out.lambdaAggregate = lamTot;
    out.p0.assign(N, zero);
    out.pK.assign(N, zero);
    out.throughputAggregate = zero;
    T Smean = zero;
    T lossWeighted = zero;
    for (std::size_t n = 0; n < N; ++n) {
        const T lamBar = lamTot - lam[n];
        const std::size_t Mn = flows[n].D0.rows();
        mam::Map<T> sup;
        sup.D0 = flows[n].D0;
        sup.D1 = flows[n].D1;
        for (std::size_t i = 0; i < Mn; ++i) {
            sup.D0(i, i) -= lamBar;
            sup.D1(i, i) += lamBar;
        }
        const MapG1kResult<T> r = qsys_mapg1k(sup, svc, K, tol, nmaxCap);
        Smean = r.meanServiceTime;
        out.p0[n] = r.p0;
        out.pK[n] = r.pK;
        out.throughput[n] = (one - r.p0) / Smean + r.pK * lamBar - lamBar;
        out.lossRatio[n] = one - out.throughput[n] / lam[n];
        out.throughputAggregate += out.throughput[n];
        lossWeighted += out.lossRatio[n] * lam[n];
    }
    out.lossAggregate = lossWeighted / lamTot;
    out.meanServiceTime = Smean;
    out.rho = lamTot * Smean;
    return out;
}

/** qsys_mapg1k_perflow with the qsys_mapg1k defaults tol = 1e-12, nmax = 200000. */
template <class T>
MapG1kPerflowResult<T> qsys_mapg1k_perflow(const std::vector<mam::Map<T>>& flows,
                                           const ServiceLaw<T>& svc, std::size_t K) {
    return qsys_mapg1k_perflow(flows, svc, K, T(num_traits<T>::from_double(1e-12)),
                               static_cast<std::size_t>(200000));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPG1K_PERFLOW_H
