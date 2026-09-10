/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_CME_H
#define LINE_API_MAM_CME_H

/**
 * Concentrated matrix exponentials, and the two-moment fit built on them.
 *
 * Templated port of matlab/src/lang/processes/CME.m and dist_fit_me.m (mirrored
 * by jline.lang.processes.CME / MEFit and the native Python fit_me_mean_scv).
 *
 * A CME of order 2n+1 is the matrix-exponential law whose squared coefficient of
 * variation is numerically minimal for that order, from the tables of Horvath,
 * Horvath and Telek. Its SCV decays as O(1/n^2), so it reaches far below the
 * Erlang bound 1/order that binds any phase-type of the same order: at 101
 * phases a CME reaches 3.9e-4 where Erlang-101 stops at 9.9e-3. The unit-mean
 * density with n harmonics is
 *   f(x) = mu1 e^{-mu1 x} ( c + sum_k a_k cos(k w mu1 x) + b_k sin(k w mu1 x) ),
 * which is alpha exp(A x) (-A e) for the block-diagonal
 *   A = blkdiag( -mu1, mu1 [-1, -k w; k w, -1], k = 1..n ).
 *
 * The coefficients come from the SAME vendored `iltcme` table the CME inverse
 * Laplace transform reads (api/mam/iltcme_table.h), so the two share one source
 * of truth; a TU using this header must link `src/api/mam/iltcme_table.cpp`.
 *
 * IT IS NOT A PHASE-TYPE. The off-diagonal entries of A are not rates -- the
 * rotation blocks carry a negative one -- so a CTMC assembled from a CME does
 * not describe the model. `sn_is_phasetype` is the test every consumer applies,
 * and `sn_nonmarkov_toph` tags the result ME rather than PH on the strength of
 * it.
 *
 * ARITHMETIC: transcendental. The two-moment fit takes a square root and the
 * table itself is a double table, so this does not instantiate under Rational.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/iltcme_table.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The unit-mean (alpha, A) form of a CME, with the SCV it attains. */
template <class T>
struct CmeRepresentation {
    std::vector<T> alpha;  ///< entry law, sums to one
    Matrix<T> A;           ///< (2n+1) square, generator-shaped but not a generator
    double scv;            ///< the tabulated cv2 of this order
};

/** Every phase count 2n+1 the vendored table realizes, ascending. */
inline std::vector<std::size_t> cme_supported_orders() {
    std::vector<std::size_t> orders;
    for (std::size_t i = 0; i < iltcme::kTableSize; ++i)
        orders.push_back(static_cast<std::size_t>(2 * iltcme::kTable[i].n + 1));
    std::sort(orders.begin(), orders.end());
    orders.erase(std::unique(orders.begin(), orders.end()), orders.end());
    return orders;
}

/**
 * The most concentrated table entry realizing `order` phases.
 *
 * Several entries share an n (the table's 'full' and 'approx' optimizations);
 * the smallest cv2 wins, which is the selection rule matlab_ilt also applies.
 */
inline const iltcme::CmeEntry& cme_table_entry(std::size_t order) {
    if (order < 3 || order % 2 == 0)
        throw InputError("cme_table_entry: the order must be an odd integer 2n+1 with n >= 1");
    if (iltcme::kTableSize == 0) throw InputError("cme_table_entry: the CME table is empty");
    const int n = static_cast<int>((order - 1) / 2);
    const iltcme::CmeEntry* best = 0;
    for (std::size_t i = 0; i < iltcme::kTableSize; ++i)
        if (iltcme::kTable[i].n == n && (best == 0 || iltcme::kTable[i].cv2 < best->cv2))
            best = &iltcme::kTable[i];
    if (best == 0)
        throw InputError("cme_table_entry: no tabulated CME of order " + std::to_string(order));
    return *best;
}

/** The minimal SCV a CME of this order attains. */
inline double cme_min_scv(std::size_t order) { return cme_table_entry(order).cv2; }

/**
 * @param order an odd phase count 2n+1 present in the table
 * @return the unit-mean representation of that order
 */
template <class T>
CmeRepresentation<T> cme_representation(std::size_t order) {
    const iltcme::CmeEntry& e = cme_table_entry(order);
    const std::size_t n = static_cast<std::size_t>(e.n), sz = 2 * n + 1;
    const T zero = num_traits<T>::from_int(0);

    CmeRepresentation<T> r;
    r.scv = e.cv2;
    r.A = Matrix<T>(sz, sz, zero);
    r.alpha.assign(sz, zero);

    const double mu1 = e.mu1, w = e.omega;
    r.A(0, 0) = num_traits<T>::from_double(-mu1);
    r.alpha[0] = num_traits<T>::from_double(e.c);
    for (std::size_t k = 1; k <= n; ++k) {
        const std::size_t i = 2 * k - 1;  // 0-based: MATLAB's 2k
        const double wk = static_cast<double>(k) * w;
        r.A(i, i) = num_traits<T>::from_double(-mu1);
        r.A(i, i + 1) = num_traits<T>::from_double(-wk * mu1);
        r.A(i + 1, i) = num_traits<T>::from_double(wk * mu1);
        r.A(i + 1, i + 1) = num_traits<T>::from_double(-mu1);
        const double d = 2.0 * (1.0 + wk * wk);
        r.alpha[i] = num_traits<T>::from_double(((1.0 + wk) * e.a[k - 1] - (1.0 - wk) * e.b[k - 1]) / d);
        r.alpha[i + 1] = num_traits<T>::from_double(((1.0 - wk) * e.a[k - 1] + (1.0 + wk) * e.b[k - 1]) / d);
    }
    T s = zero;
    for (std::size_t i = 0; i < sz; ++i) s += r.alpha[i];
    if (!(num_traits<T>::to_double(s) > 0.0))
        throw NumericError("cme_representation: the entry law does not normalize");
    for (std::size_t i = 0; i < sz; ++i) r.alpha[i] = T(r.alpha[i] / s);
    return r;
}

/** Assemble the renewal (D0, D1) of a matrix-exponential law (alpha, A). */
template <class T>
Map<T> me_to_map(const std::vector<T>& alpha, const Matrix<T>& A) {
    const std::size_t n = A.rows();
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = A;
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < n; ++j) rowsum += A(i, j);
        for (std::size_t j = 0; j < n; ++j) m.D1(i, j) = T(-rowsum * alpha[j]);
    }
    return m;
}

/**
 * Two-moment matrix-exponential fit for 0 < scv < 1, a port of dist_fit_me.m.
 *
 * The fit is the convolution X = c Y + Z of a scaled unit-mean CME Y with an
 * independent exponential Z. Matching c + d = mean and c^2 sY + d^2 = scv mean^2
 * gives
 *   c = mean (1 - sqrt(1 - (1+sY)(1-scv))) / (1 + sY),   d = mean - c,
 * so every target in [sY/(1+sY), 1] is hit EXACTLY in 2n+2 phases. The
 * exponential tail is what lets the convolution reach up to SCV 1; the
 * concentrated part is what lets it reach far below the Erlang bound.
 *
 * BUDGET-LIMITED IS NOT AN ERROR. When maxPhases cannot buy an order whose reach
 * covers the target, the most concentrated affordable member is returned and the
 * caller gets the closest achievable SCV, rather than a silent Erlang.
 *
 * @param mean      target mean, positive and finite
 * @param scv       target SCV, strictly inside (0,1)
 * @param maxPhases cap on the phase count, 0 for no cap
 */
template <class T>
Map<T> dist_fit_me(double mean, double scv, std::size_t maxPhases = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "dist_fit_me takes a square root of the moment discriminant");
    if (!std::isfinite(mean) || mean <= 0.0)
        throw InputError("dist_fit_me: the mean must be a positive finite number");
    if (!std::isfinite(scv) || scv <= 0.0 || scv >= 1.0)
        throw InputError(
            "dist_fit_me: requires 0 < scv < 1; use a hyperexponential for scv >= 1 and a CME "
            "for scv = 0");

    const std::vector<std::size_t> orders = cme_supported_orders();
    std::size_t bestOrder = 0;
    for (std::size_t i = 0; i < orders.size(); ++i) {
        const std::size_t order = orders[i];
        if (maxPhases > 0 && order + 1 > maxPhases) continue;
        const double sY = cme_min_scv(order);
        bestOrder = order;  // budget-limited: keep the most concentrated one that fits
        if (sY / (1.0 + sY) <= scv) break;
    }
    if (bestOrder == 0)
        throw InputError("dist_fit_me: no CME order fits a budget of " + std::to_string(maxPhases) +
                         " phases; the smallest is 3 phases plus one exponential");

    const CmeRepresentation<T> rep = cme_representation<T>(bestOrder);
    const double sY = rep.scv;
    const double reach = sY / (1.0 + sY);
    const double c = scv < reach ? mean / (1.0 + sY)
                                 : mean * (1.0 - std::sqrt(1.0 - (1.0 + sY) * (1.0 - scv))) /
                                       (1.0 + sY);
    const double d = mean - c;
    const std::size_t n = rep.alpha.size();
    const T zero = num_traits<T>::from_int(0);

    if (d <= mean * 1e-12) {  // the whole mass is in the concentrated part
        Matrix<T> A = rep.A;
        const T inv = num_traits<T>::from_double(1.0 / mean);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = T(A(i, j) * inv);
        return me_to_map(rep.alpha, A);
    }
    if (c <= mean * 1e-12) {  // degenerates to the exponential tail alone
        Matrix<T> A(1, 1, num_traits<T>::from_double(-1.0 / mean));
        std::vector<T> alpha(1, num_traits<T>::from_int(1));
        return me_to_map(alpha, A);
    }

    // Convolution: the exit flow of the CME block feeds the exponential phase.
    Matrix<T> A(n + 1, n + 1, zero);
    std::vector<T> alpha(n + 1, zero);
    const T invc = num_traits<T>::from_double(1.0 / c);
    for (std::size_t i = 0; i < n; ++i) {
        alpha[i] = rep.alpha[i];
        T rowsum = zero;
        for (std::size_t j = 0; j < n; ++j) {
            A(i, j) = T(rep.A(i, j) * invc);
            rowsum += A(i, j);
        }
        A(i, n) = -rowsum;
    }
    A(n, n) = num_traits<T>::from_double(-1.0 / d);
    return me_to_map(alpha, A);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_CME_H
