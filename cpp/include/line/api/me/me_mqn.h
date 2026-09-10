/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_MQN_H
#define LINE_API_ME_ME_MQN_H

/**
 * Maximum-entropy algorithm for mixed open/closed multiclass networks.
 *
 * Templated port of matlab/src/api/me/me_mqn.m, cross-checked against
 * jar/src/main/java/jline/api/nc/Me_mqn.java. Kouvatsos (1994) notes that the
 * closed two-stage treatment carries over to mixed networks (Section 2.3) but
 * gives no algorithm; both LINE implementations compose the open (3.2) and
 * closed (3.3) algorithms by product-form-style conditioning:
 *
 *   1. the open classes are solved by the open GE-type fixed point on the
 *      station set, ignoring the closed classes;
 *   2. the closed classes are solved on servers whose capacity is reduced by
 *      the open-class utilization, mu_c(i,r) = mu(i,r) (1 - rho_o(i));
 *   3. the open mean queue lengths are inflated by the closed occupancy,
 *      L_o(i,r) <- L_o(i,r) (1 + Lc(i)), at single-server stations.
 *
 * Steps 2-3 are exact in the BCMP product-form limit, where they reduce to
 * the classical mixed MVA treatment, and are GE-type approximations
 * otherwise. Only single-server and infinite-server stations are supported.
 *
 * ARITHMETIC: composes two tolerance-stopped fixed points.
 *   static_assert(num_traits<T>::has_transcendental)
 */

#include <cstddef>
#include <vector>

#include "line/api/me/me_cqn.h"
#include "line/api/me/me_oqn.h"
#include "line/api/me/me_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace me {

namespace detail {

/** Selects the given columns of an (M x R) matrix. */
template <class T>
Matrix<T> select_cols(const Matrix<T>& A, const std::vector<std::size_t>& cols) {
    Matrix<T> B(A.rows(), cols.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t k = 0; k < cols.size(); ++k) B(i, k) = A(i, cols[k]);
    return B;
}

}  // namespace detail

/**
 * @param M            number of stations
 * @param R            number of classes
 * @param open_classes per-class flag, nonzero for an open class
 * @param lambda0      external arrival rates (M x R), zero for closed classes
 * @param Ca0          external arrival scvs (M x R)
 * @param N            populations (R); the entries of open classes are unused
 * @param mu           service rates (M x R)
 * @param Cs           service scvs (M x R)
 * @param P            routing, R matrices (M x M)
 * @param c            servers per station, 0 for an infinite-server station
 * @param refstat      reference station per class, -1 for the default
 * @param insens       insensitive discipline flags per station
 * @param opt          tolerance and iteration budget
 */
template <class T>
MeResult<T> me_mqn(std::size_t M, std::size_t R, const std::vector<char>& open_classes,
                   const Matrix<T>& lambda0, const Matrix<T>& Ca0, const std::vector<long>& N,
                   const Matrix<T>& mu, const Matrix<T>& Cs, const std::vector<Matrix<T>>& P,
                   const std::vector<long>& c, const std::vector<long>& refstat,
                   const std::vector<char>& insens, const MeOptions& opt = MeOptions()) {
    static_assert(num_traits<T>::has_transcendental, "me_mqn requires transcendental arithmetic");
    detail::check_dims(M, R, mu, Cs, P, c, insens, "me_mqn");
    if (open_classes.size() != R) throw InputError("me_mqn: one open flag per class");
    if (N.size() != R) throw InputError("me_mqn: one population per class");
    if (refstat.size() != R) throw InputError("me_mqn: one reference station per class");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::size_t> oc, cc;
    for (std::size_t r = 0; r < R; ++r) (open_classes[r] ? oc : cc).push_back(r);

    MeResult<T> out;
    out.L = Matrix<T>(M, R, zero);
    out.W = Matrix<T>(M, R, zero);
    out.Ca = Matrix<T>(M, R, one);
    out.Cd = Matrix<T>(M, R, one);
    out.lambda = Matrix<T>(M, R, zero);
    out.rho = Matrix<T>(M, R, zero);
    out.X.assign(R, zero);
    out.converged = true;

    // Step 1: open classes
    std::vector<T> rho_o(M, zero);
    if (!oc.empty()) {
        std::vector<Matrix<T>> Po;
        for (std::size_t k = 0; k < oc.size(); ++k) Po.push_back(P[oc[k]]);
        const MeResult<T> ro =
            me_oqn(M, oc.size(), detail::select_cols(lambda0, oc), detail::select_cols(Ca0, oc),
                   detail::select_cols(mu, oc), detail::select_cols(Cs, oc), Po, c, insens, opt);
        out.iter += ro.iter;
        out.converged = out.converged && ro.converged;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < oc.size(); ++k) {
                out.L(i, oc[k]) = ro.L(i, k);
                out.Ca(i, oc[k]) = ro.Ca(i, k);
                out.Cd(i, oc[k]) = ro.Cd(i, k);
                out.lambda(i, oc[k]) = ro.lambda(i, k);
                out.rho(i, oc[k]) = ro.rho(i, k);
            }
        for (std::size_t i = 0; i < M; ++i) {
            if (detail::is_is(c, i)) continue;
            for (std::size_t k = 0; k < oc.size(); ++k) rho_o[i] += ro.rho(i, k);
        }
        for (std::size_t k = 0; k < oc.size(); ++k) out.X[oc[k]] = ro.X[k];
    }

    // Step 2: closed classes on servers with reduced capacity
    if (!cc.empty()) {
        Matrix<T> mu_c = detail::select_cols(mu, cc);
        for (std::size_t i = 0; i < M; ++i) {
            if (detail::is_is(c, i)) continue;
            T fac = one - rho_o[i];
            if (fac < zero) fac = zero;
            for (std::size_t k = 0; k < cc.size(); ++k) mu_c(i, k) *= fac;
        }
        std::vector<Matrix<T>> Pc;
        std::vector<long> Nc, refc;
        for (std::size_t k = 0; k < cc.size(); ++k) {
            Pc.push_back(P[cc[k]]);
            Nc.push_back(N[cc[k]]);
            refc.push_back(refstat[cc[k]]);
        }
        const MeResult<T> rc = me_cqn(M, cc.size(), Nc, mu_c, detail::select_cols(Cs, cc), Pc, c,
                                      refc, insens, opt);
        out.iter += rc.iter;
        out.converged = out.converged && rc.converged;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < cc.size(); ++k) {
                out.L(i, cc[k]) = rc.L(i, k);
                out.W(i, cc[k]) = rc.W(i, k);
                out.Ca(i, cc[k]) = rc.Ca(i, k);
                out.Cd(i, cc[k]) = rc.Cd(i, k);
                out.lambda(i, cc[k]) = rc.lambda(i, k);
                if (detail::is_is(c, i)) {
                    out.rho(i, cc[k]) = rc.rho(i, k);
                } else {
                    T fac = one - rho_o[i];
                    if (fac < zero) fac = zero;
                    out.rho(i, cc[k]) = rc.rho(i, k) * fac;
                }
            }
        for (std::size_t k = 0; k < cc.size(); ++k) out.X[cc[k]] = rc.X[k];
    }

    // Step 3: inflate the open queue lengths by the closed occupancy
    if (!oc.empty() && !cc.empty()) {
        for (std::size_t i = 0; i < M; ++i) {
            if (detail::is_is(c, i)) continue;
            T Lc_i = zero;
            for (std::size_t k = 0; k < cc.size(); ++k) Lc_i += out.L(i, cc[k]);
            for (std::size_t k = 0; k < oc.size(); ++k) out.L(i, oc[k]) *= (one + Lc_i);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < oc.size(); ++k)
            if (out.lambda(i, oc[k]) > zero)
                out.W(i, oc[k]) = out.L(i, oc[k]) / out.lambda(i, oc[k]);
    return out;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_MQN_H
