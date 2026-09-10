/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_COMOMRM_LD_H
#define LINE_API_PFQN_COMOMRM_LD_H

/**
 * CoMoM for the repairman model with an arbitrary LOAD-DEPENDENT rate lattice
 * at the single queueing station.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_comomrm_ld.m. The recursion is
 * the bidiagonal transfer-matrix sweep of pfqn_comomrm_ms, from which this
 * routine differs only in where the rate lattice comes from and in the
 * preprocessing that folds pure delay rows out of the demand matrix.
 *
 * Delay detection. When no think time is supplied, a station whose rate row is
 * exactly the identity lattice mu(i,k) = k is an infinite server, and the
 * reference moves its demand row into Z and drops it from L and mu. MATLAB
 * tests this with a 2-norm against a tolerance; this port tests exact equality,
 * for the same reason pfqn_unique merges on exact equality: a tolerant test
 * silently reclassifies a station and changes the model, which has no meaning
 * in the rational field. A caller wanting the tolerant behaviour should round
 * its rate lattice before calling, where the rounding is visible.
 *
 * Arithmetic: EXACT-CAPABLE, unconditionally. See pfqn_comomrm_ms for why the
 * reference's per-step renormalization is dropped rather than reproduced.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comomrm_ms.h"
#include "line/api/pfqn/pfqn_nc_sanitize.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L  (M x R) demands
 * @param N  (R) populations
 * @param Z  (K x R) think times; empty or zero triggers delay detection on mu
 * @param mu (M x >=Nt) load-dependent rate lattice
 */
template <class T>
ComomRmResult<T> pfqn_comomrm_ld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                                 const Matrix<T>& mu) {
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    int Nt = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_comomrm_ld: negative population");
        Nt += v;
    }

    std::size_t M = L.empty() ? 0 : L.rows();
    if (M > 0 && L.cols() != R)
        throw InputError("pfqn_comomrm_ld: L and N disagree on the class count");

    // Aggregate the think-time rows.
    std::vector<T> Zsum(R, zero);
    T Ztot = zero;
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_comomrm_ld: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }
    for (std::size_t r = 0; r < R; ++r) Ztot += Zsum[r];

    Matrix<T> Lq = L;
    Matrix<T> muq = mu;
    if (Ztot == zero && M > 0) {
        if (mu.rows() != M)
            throw InputError("pfqn_comomrm_ld: mu and L disagree on the station count");
        // A row whose rate lattice is 1, 2, 3, ... is an infinite server.
        std::vector<std::size_t> queues, delays;
        for (std::size_t i = 0; i < M; ++i) {
            bool isDelay = Nt > 0;
            for (int k = 1; k <= Nt; ++k)
                if (mu(i, static_cast<std::size_t>(k - 1)) != num_traits<T>::from_int(k)) {
                    isDelay = false;
                    break;
                }
            (isDelay ? delays : queues).push_back(i);
        }
        if (!delays.empty()) {
            for (std::size_t d = 0; d < delays.size(); ++d)
                for (std::size_t r = 0; r < R; ++r) Zsum[r] += L(delays[d], r);
            Lq = Matrix<T>(queues.size(), R);
            muq = Matrix<T>(queues.size(), mu.cols());
            for (std::size_t q = 0; q < queues.size(); ++q) {
                for (std::size_t r = 0; r < R; ++r) Lq(q, r) = L(queues[q], r);
                for (std::size_t k = 0; k < mu.cols(); ++k) muq(q, k) = mu(queues[q], k);
            }
            M = queues.size();
        }
    }

    Matrix<T> Zmat(1, R);
    for (std::size_t r = 0; r < R; ++r) Zmat(0, r) = Zsum[r];

    T Lsum = zero;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lsum += Lq(i, r);
    if (M == 0 || Lsum == zero) {
        // Only delays remain: the model is the trivial multinomial one.
        const NcResult<T> ca = pfqn_ca(Matrix<T>(), N, Zmat);
        ComomRmResult<T> res;
        res.G = ca.G;
        res.lG = ca.lG;
        res.prob.assign(static_cast<std::size_t>(Nt) + 1, zero);
        res.prob.back() = one;
        return res;
    }

    const NcSanitizeResult<T> san = pfqn_nc_sanitize(Lq, N, Zmat);
    if (san.L.rows() > 1)
        throw InputError("pfqn_comomrm_ld: the solver accepts at most a single queueing station");

    int Ntk = 0;
    for (int v : san.N) Ntk += v;
    if (Ntk == 0 || san.L.empty()) {
        ComomRmResult<T> res;
        res.G = san.Gremaind;
        res.lG = san.lGremaind;
        res.prob.assign(static_cast<std::size_t>(Nt) + 1, zero);
        res.prob[0] = one;
        return res;
    }

    const std::size_t Rk = san.N.size();
    std::vector<T> Lv(Rk, zero), Zv(Rk, zero);
    for (std::size_t r = 0; r < Rk; ++r) {
        Lv[r] = san.L(0, r);
        for (std::size_t k = 0; k < san.Z.rows(); ++k) Zv[r] += san.Z(k, r);
    }
    std::vector<T> muv(static_cast<std::size_t>(Ntk), one);
    for (int k = 0; k < Ntk; ++k) muv[static_cast<std::size_t>(k)] = muq(0, static_cast<std::size_t>(k));

    return detail::comomrm_finish(detail::comomrm_bidiag(Lv, san.N, Zv, muv), san.Gremaind);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_COMOMRM_LD_H
