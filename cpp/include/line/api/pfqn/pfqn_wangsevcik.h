/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_WANGSEVCIK_H
#define LINE_API_PFQN_PFQN_WANGSEVCIK_H

/**
 * Two approximate MVA schemes: Queue-Line and Fraction-Line.
 *
 * Port of `pfqn_qli` and `pfqn_fli` from python/line_solver/api/pfqn/mva.py.
 * PYTHON-ONLY: neither MATLAB nor the JAR carries them as standalone entry
 * points, so native Python is the reference.
 *
 * Reference: W. Wang and K. C. Sevcik, "Performance Models for Multiprogrammed
 * Systems", IBM Research Report RC 5925, 1976.
 *
 * `pfqn_qdlin` USED TO LIVE HERE AND NO LONGER DOES. Its Wang-Sevcik arm scaled
 * the job's own-class contribution by (N_r - 1)/N_r and left the other classes
 * alone, which is Bard-Schweitzer written out, so the function reproduced
 * `pfqn_bs` to iteration tolerance and was neither a Linearizer nor
 * queue-dependent. It is now `line/api/pfqn/pfqn_qdlin.h`, the array-level twin
 * of what SolverMVA computes for method='qdlin'.
 *
 * BOTH ANSWER ONE QUESTION: what queue does an ARRIVING class-r job see?
 * Exact MVA says it sees the queue at population N - e_r, which is why the
 * exact recursion has to walk the whole population lattice. Every approximation
 * here replaces that by a correction applied to the queue at the FULL
 * population, and the two differ only in the correction:
 *
 *  - QLI subtracts a `1/(N_r - 1)` weighted discrepancy between the own-class
 *    queue and its DEMAND-PROPORTIONAL share, so a station that holds more of
 *    the class than its demand warrants is discounted more.
 *  - FLI uses the same proportional share but combines it differently, with a
 *    `2/N_r` coefficient and the share ADDED rather than subtracted.
 *
 * THE `Q_seen` FLOOR AT ZERO IS LOAD-BEARING, not defensive. Both Wang-Sevcik
 * corrections are differences of estimates and can go negative at a lightly
 * loaded station; a negative queue would make the residence time SHORTER than
 * the service demand, which is impossible, and the iteration then diverges away
 * from the fixed point rather than toward it. The reference clamps and so does
 * this.
 *
 * THE FALLBACK ARM IS NOT THE SAME FORMULA. When the denominator vanishes, or
 * when the population is too small for the correction's own divisor (`N_r > 1`
 * for QLI), the reference falls back to `Q_total - Q_own`, i.e. the other
 * classes only. That is the Bard-Schweitzer arrival estimate, and it is a
 * DIFFERENT approximation, so a caller comparing two runs across that boundary
 * is comparing two schemes.
 *
 * ARITHMETIC: field.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** What an approximate MVA sweep reports. */
template <class T>
struct WsResult {
    Matrix<T> Q;             ///< (M x R) mean queue lengths
    Matrix<T> U;             ///< (M x R) utilizations
    Matrix<T> R;             ///< (M x R) residence times
    std::vector<T> X;        ///< (R) class throughputs
    std::vector<T> C;        ///< (R) cycle times
    std::size_t iterations = 0;
};

/** Which arrival-queue correction the sweep applies. */
enum class WsScheme { Qli = 0, Fli };

namespace wsdetail {

/** The shared initial guess: each class spread over the stations by demand. */
template <class T>
Matrix<T> proportional_start(const Matrix<T>& L, const std::vector<T>& N) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> Q(M, R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        T s = zero;
        for (std::size_t k = 0; k < M; ++k) s += L(k, r);
        if (s == zero) s = one;  // a class with no demand anywhere
        for (std::size_t k = 0; k < M; ++k) Q(k, r) = L(k, r) / s * N[r];
    }
    return Q;
}

}  // namespace wsdetail

/**
 * One approximate MVA sweep, by the chosen arrival-queue correction.
 *
 * @param L   (M x R) service demands
 * @param N   (R) class populations
 * @param Z   (R) think times; empty means none
 * @param tol convergence tolerance on the queue lengths
 */
template <class T>
WsResult<T> pfqn_wangsevcik(const Matrix<T>& L, const std::vector<T>& N,
                            const std::vector<T>& Z, WsScheme scheme, double tol = 1e-6,
                            std::size_t max_iter = 1000) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M == 0 || R == 0) throw InputError("pfqn_wangsevcik: the demand matrix is empty");
    if (N.size() != R) throw InputError("pfqn_wangsevcik: N has the wrong length");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_wangsevcik: Z has the wrong length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    WsResult<T> out;
    out.Q = Matrix<T>(M, R, zero);
    out.U = Matrix<T>(M, R, zero);
    out.R = Matrix<T>(M, R, zero);
    out.X.assign(R, zero);
    out.C.assign(R, zero);

    T Ntot = zero;
    for (std::size_t r = 0; r < R; ++r) Ntot += N[r];
    if (!(num_traits<T>::to_double(Ntot) > 0.0)) return out;  // an empty network is all zeros

    out.Q = wsdetail::proportional_start(L, N);
    Matrix<T> Qprev(M, R, zero);
    for (std::size_t k = 0; k < M; ++k)
        for (std::size_t r = 0; r < R; ++r)
            Qprev(k, r) = out.Q(k, r) * num_traits<T>::from_int(10);

    while (out.iterations < max_iter) {
        double gap = 0.0;
        for (std::size_t k = 0; k < M; ++k)
            for (std::size_t r = 0; r < R; ++r)
                gap = std::max(gap, std::fabs(num_traits<T>::to_double(out.Q(k, r) - Qprev(k, r))));
        if (gap <= tol) break;
        ++out.iterations;
        Qprev = out.Q;

        // The per-station totals of the PREVIOUS iterate, which every
        // correction reads; computing them inside the station loop would make
        // the sweep depend on the update order.
        std::vector<T> Qtot(M, zero);
        for (std::size_t k = 0; k < M; ++k)
            for (std::size_t s = 0; s < R; ++s) Qtot[k] += Qprev(k, s);

        for (std::size_t r = 0; r < R; ++r) {
            if (!(num_traits<T>::to_double(N[r]) > 0.0)) continue;

            // The demand-proportional share both Wang-Sevcik schemes use.
            T qliden = zero;
            for (std::size_t m = 0; m < M; ++m)
                if (num_traits<T>::to_double(L(m, r)) > 0.0)
                    qliden += L(m, r) * (one + Qtot[m] - Qprev(m, r));

            for (std::size_t k = 0; k < M; ++k) {
                T qseen = zero;
                {
                    const T qlinum = L(k, r) * (one + Qtot[k] - Qprev(k, r));
                    const bool usable = num_traits<T>::to_double(qliden) > 0.0 &&
                                        (scheme == WsScheme::Qli
                                             ? num_traits<T>::to_double(N[r]) > 1.0
                                             : num_traits<T>::to_double(N[r]) > 0.0);
                    if (!usable) {
                        // Bard-Schweitzer: the other classes only. A DIFFERENT
                        // approximation, not a degenerate case of the same one.
                        qseen = T(Qtot[k] - Qprev(k, r));
                    } else if (scheme == WsScheme::Qli) {
                        qseen = T(Qtot[k] - (one / (N[r] - one)) * (Qprev(k, r) - qlinum / qliden));
                    } else {
                        qseen = T(Qtot[k] - (num_traits<T>::from_int(2) / N[r]) * Qprev(k, r) +
                                  qlinum / qliden);
                    }
                }
                // A negative arrival queue would make the residence time
                // shorter than the demand, which is impossible, and the
                // iteration then walks away from the fixed point.
                if (num_traits<T>::to_double(qseen) < 0.0) qseen = zero;
                out.R(k, r) = L(k, r) * (one + qseen);
            }

            T Rtot = zero;
            for (std::size_t k = 0; k < M; ++k) Rtot += out.R(k, r);
            const T zr = Z.empty() ? zero : Z[r];
            out.X[r] = (num_traits<T>::to_double(zr + Rtot) > 0.0) ? T(N[r] / (zr + Rtot)) : zero;
            for (std::size_t k = 0; k < M; ++k) {
                out.Q(k, r) = out.X[r] * out.R(k, r);
                out.U(k, r) = out.X[r] * L(k, r);
            }
            out.C[r] = Rtot;
        }
    }
    return out;
}

/** Wang-Sevcik Queue-Line. */
template <class T>
WsResult<T> pfqn_qli(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                     double tol = 1e-6, std::size_t max_iter = 1000) {
    return pfqn_wangsevcik(L, N, Z, WsScheme::Qli, tol, max_iter);
}

/** Wang-Sevcik Fraction-Line. */
template <class T>
WsResult<T> pfqn_fli(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                     double tol = 1e-6, std::size_t max_iter = 1000) {
    return pfqn_wangsevcik(L, N, Z, WsScheme::Fli, tol, max_iter);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_WANGSEVCIK_H
