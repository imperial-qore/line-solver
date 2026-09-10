/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_OI_INSVC_H
#define LINE_API_PFQN_PFQN_OI_INSVC_H

/**
 * Conditional mean number of IN-SERVICE jobs per class at an order-independent
 * station.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_oi_insvc.m. This is the quantity
 * behind LINE's utilization convention at OI stations, U_r = E[sir_r]/c, where
 * sir_r counts class-r JOBS receiving a strictly positive rate (a job served
 * concurrently by several servers counts once). Conditioning on the tail
 * element of the ordering gives the balanced-fairness recursion for the OI
 * balance function and its sir-weighted companion,
 *
 *   Phi(0) = 1,   Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r)
 *   Xi_r(0) = 0,
 *   Xi_r(n) = (1/mu(n)) [ sum_s Xi_r(n - e_s)
 *                         + 1{n_r>0} 1{mu(n) > mu(n - e_r)} Phi(n - e_r) ],
 *
 * and E[sir_r | n] = Xi_r(n)/Phi(n).
 *
 * UNREACHABLE STATES. A composition no server can serve has mu(n) <= 0; the
 * reference leaves Phi and Xi at zero there and the port does the same, rather
 * than dividing by zero or substituting a rate. g is then zero at that state,
 * which is the correct reading: the state carries no weight.
 *
 * ARITHMETIC. Additions and divisions only, so the routine is EXACT in
 * rational arithmetic and is deliberately left ungated. The strict comparison
 * mu(n) > mu(n - e_r) that decides whether the tail job is in service is an
 * exact comparison there, which matters: in floating point two rates that are
 * equal in the model can differ in the last bit and flip that indicator.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_oi_insvc, mirroring [g, Xi, Phi]. */
template <class T>
struct OiInsvcResult {
    Matrix<T> g;                      ///< (prod(N+1) x R) E[sir_r | n]
    Matrix<T> Xi;                     ///< (prod(N+1) x R) sir-weighted balance
    std::vector<T> Phi;               ///< (prod(N+1)) OI balance function
    std::vector<std::size_t> stride;  ///< column-major strides, for indexing
};

/**
 * @param oirate OI total service rate mu(n) for a per-class count vector
 * @param N      (R) closed population vector
 */
template <class T>
OiInsvcResult<T> pfqn_oi_insvc(const std::function<T(const std::vector<int>&)>& oirate,
                               const std::vector<int>& N) {
    if (!oirate) throw InputError("pfqn_oi_insvc: oirate must be callable");
    const std::size_t R = N.size();
    if (R == 0) throw InputError("pfqn_oi_insvc: empty population vector");
    std::vector<std::size_t> shp(R), stride(R, 1);
    std::size_t total = 1;
    for (std::size_t d = 0; d < R; ++d) {
        if (N[d] < 0) throw InputError("pfqn_oi_insvc: negative population");
        shp[d] = static_cast<std::size_t>(N[d]) + 1;
    }
    for (std::size_t d = 1; d < R; ++d) stride[d] = stride[d - 1] * shp[d - 1];
    for (std::size_t d = 0; d < R; ++d) total *= shp[d];

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::vector<int>> subs(total, std::vector<int>(R, 0));
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t li = i;
        for (std::size_t d = 0; d < R; ++d) {
            subs[i][d] = static_cast<int>(li % shp[d]);
            li /= shp[d];
        }
    }
    std::vector<T> muv(total, zero);
    for (std::size_t i = 0; i < total; ++i) {
        int tot = 0;
        for (int v : subs[i]) tot += v;
        if (tot > 0) muv[i] = oirate(subs[i]);
    }

    std::vector<T> Phi(total, zero);
    Matrix<T> Xi(total, R, zero);
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<int>& n = subs[i];
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            Phi[i] = one;
            continue;
        }
        const T mun = muv[i];
        if (mun <= zero) continue;  // unreachable state, zero weight
        T sPhi = zero;
        std::vector<T> sXi(R, zero);
        for (std::size_t s = 0; s < R; ++s) {
            if (n[s] <= 0) continue;
            const std::size_t j = i - stride[s];
            sPhi += Phi[j];
            for (std::size_t r = 0; r < R; ++r) sXi[r] += Xi(j, r);
        }
        Phi[i] = T(sPhi / mun);
        for (std::size_t r = 0; r < R; ++r) {
            T acc = sXi[r];
            if (n[r] > 0) {
                const std::size_t j = i - stride[r];
                if (mun > muv[j]) acc += Phi[j];  // the tail class-r job is in service
            }
            Xi(i, r) = T(acc / mun);
        }
    }

    OiInsvcResult<T> res;
    res.g = Matrix<T>(total, R, zero);
    for (std::size_t i = 0; i < total; ++i)
        if (Phi[i] > zero)
            for (std::size_t r = 0; r < R; ++r) res.g(i, r) = T(Xi(i, r) / Phi[i]);
    res.Xi = Xi;
    res.Phi = Phi;
    res.stride = stride;
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_OI_INSVC_H
