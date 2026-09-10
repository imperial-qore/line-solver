/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_JOINT_H
#define LINE_API_PFQN_JOINT_H

/**
 * Joint queue-length probability of a closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_joint.m, whose single entry point
 * dispatches on the SHAPE of its first argument; the two behaviours are split
 * into two named functions here, since C++ can tell an (M) vector from an
 * (M x R) matrix at the call site and a shape-dispatching overload would only
 * hide the distinction.
 *
 * PER-CLASS form, pfqn_joint. For a per-station, per-class occupancy n(i,r),
 * the unnormalized weight is the product of the per-station multinomial terms
 * and the delay term,
 *
 *   F(n) = prod_r Z_r^{n0_r} / n0_r!  *  prod_i [ (sum_r n_ir)! / prod_r n_ir!
 *                                                 * prod_r L(i,r)^{n_ir} ],
 *
 * with n0 = N - sum_i n(i,:) the jobs left at the delay, and the probability is
 * F(n)/G(N).
 *
 * TOTAL form, pfqn_joint_total. When only the per-station TOTALS m(i) are
 * given, the per-class split is unknown and the weight is a PERMANENT: it sums
 * the product-form weight over every assignment of the N_r class-r jobs to the
 * m(i) slots of station i. The reference builds that permanent by expanding
 * every station row m(i) times and every class column N_r times and calling
 * perm(); this port calls pfqn_perm directly on the (sum m) x R matrix with
 * column multiplicities N, which is the same quantity without materializing
 * the expansion. The delay contributes its own 1/n0! as in the reference.
 *
 * Arithmetic: EXACT-CAPABLE, both forms. The reference works in log space
 * throughout (log/gammaln/exp) purely to keep the factorials in range; every
 * quantity involved is a ratio of products of the inputs and of factorials, so
 * it is formed directly here and the probability comes out as an exact
 * rational. Note the cost of the total form: pfqn_perm is prod_r (N_r + 1)
 * evaluations of a length-(sum m) product, so it is a small-model routine.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_jointmarg.h"
#include "line/api/pfqn/pfqn_perm.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * Joint probability of a PER-CLASS occupancy matrix.
 *
 * @param n (M x R) per-station, per-class occupancy
 * @param L (M x R) service demands
 * @param N (R) populations
 * @param Z (K x R) think times, summed over rows; may be empty
 * @param G the normalizing constant G(N)
 */
template <class T>
T pfqn_joint(const Matrix<int>& n, const Matrix<T>& L, const std::vector<int>& N,
             const Matrix<T>& Z, const T& G) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_joint: L and N disagree on the class count");
    if (n.rows() != M || n.cols() != R)
        throw InputError("pfqn_joint: the occupancy matrix has the wrong shape");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (G == zero) throw NumericError("pfqn_joint: the normalizing constant is zero");

    std::vector<T> Zsum(R, zero);
    T Ztot = zero;
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_joint: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }
    for (std::size_t r = 0; r < R; ++r) Ztot += Zsum[r];

    // Jobs left at the delay.
    std::vector<int> n0(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        int used = 0;
        for (std::size_t i = 0; i < M; ++i) used += n(i, r);
        n0[r] = N[r] - used;
        if (n0[r] < 0) throw InputError("pfqn_joint: the occupancy exceeds the population");
    }

    T F = one;
    if (Ztot > zero) {
        for (std::size_t r = 0; r < R; ++r) {
            if (n0[r] == 0) continue;
            F *= num_pow_int(Zsum[r], static_cast<unsigned>(n0[r])) /
                 num_factorial<T>(static_cast<unsigned>(n0[r]));
        }
    } else {
        for (std::size_t r = 0; r < R; ++r)
            if (n0[r] != 0) return zero;  // no delay to hold them
    }
    for (std::size_t i = 0; i < M; ++i) {
        int tot = 0;
        for (std::size_t r = 0; r < R; ++r) tot += n(i, r);
        T term = num_factorial<T>(static_cast<unsigned>(tot));
        for (std::size_t r = 0; r < R; ++r) {
            if (n(i, r) == 0) continue;
            term *= num_pow_int(L(i, r), static_cast<unsigned>(n(i, r))) /
                    num_factorial<T>(static_cast<unsigned>(n(i, r)));
        }
        F *= term;
    }
    return F / G;
}

/**
 * Joint probability of the per-station TOTAL queue lengths.
 *
 * @param m (M) per-station total occupancy; the delay takes the remainder
 * @param L (M x R) service demands
 * @param N (R) populations
 * @param Z (K x R) think times, summed over rows; may be empty
 * @param G the normalizing constant G(N)
 */
template <class T>
T pfqn_joint_total(const std::vector<int>& m, const Matrix<T>& L, const std::vector<int>& N,
                   const Matrix<T>& Z, const T& G) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_joint_total: L and N disagree on the class count");
    if (m.size() != M) throw InputError("pfqn_joint_total: the occupancy vector has the wrong length");
    const T zero = num_traits<T>::from_int(0);
    if (G == zero) throw NumericError("pfqn_joint_total: the normalizing constant is zero");

    std::vector<T> Zsum(R, zero);
    T Ztot = zero;
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_joint_total: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }
    for (std::size_t r = 0; r < R; ++r) Ztot += Zsum[r];

    long Ntot = 0, mtot = 0;
    for (int v : N) Ntot += v;
    for (int v : m) {
        if (v < 0) throw InputError("pfqn_joint_total: negative occupancy");
        mtot += v;
    }
    const long n0 = Ntot - mtot;
    if (n0 < 0) throw InputError("pfqn_joint_total: the occupancy exceeds the population");
    if (n0 > 0 && !(Ztot > zero)) return zero;  // no delay to hold the remainder

    // The aggregated think time is ONE extra infinite-server row holding the
    // remainder, which is exactly the shape pfqn_jointmarg takes. Delegating
    // keeps the identity in one place; a model with SEVERAL delays needs the
    // general form, because each of them carries its own 1/n_j!.
    Matrix<T> Lall(M + 1, R);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lall(i, r) = L(i, r);
    for (std::size_t r = 0; r < R; ++r) Lall(M, r) = Zsum[r];
    std::vector<int> nall(m);
    nall.push_back(static_cast<int>(n0));
    std::vector<std::size_t> infset(1, M);
    return pfqn_jointmarg(nall, Lall, N, infset, G);
}

/** Overload computing G with pfqn_ca first, matching the reference's default. */
template <class T>
T pfqn_joint(const Matrix<int>& n, const Matrix<T>& L, const std::vector<int>& N,
             const Matrix<T>& Z) {
    return pfqn_joint(n, L, N, Z, pfqn_ca(L, N, Z).G);
}

template <class T>
T pfqn_joint_total(const std::vector<int>& m, const Matrix<T>& L, const std::vector<int>& N,
                   const Matrix<T>& Z) {
    return pfqn_joint_total(m, L, N, Z, pfqn_ca(L, N, Z).G);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_JOINT_H
