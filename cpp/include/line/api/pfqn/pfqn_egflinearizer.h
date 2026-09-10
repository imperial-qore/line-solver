/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_EGFLINEARIZER_H
#define LINE_API_PFQN_EGFLINEARIZER_H

/**
 * Extended generalized fixed-point Linearizer (De Souza e Silva and Muntz's
 * generalization of Chandy and Neuse's Linearizer, with a per-class scaling
 * exponent alpha_r).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_egflinearizer.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/mva/Pfqn_egflinearizer.java. This
 * is the single implementation behind pfqn_linearizer (alpha == 1) and
 * pfqn_gflinearizer (alpha uniform), which are thin wrappers.
 *
 * The algorithm carries the queue lengths at the full population and at each
 * of the R reduced populations N - e_s, and a correction
 *
 *   Delta(i,r,s) = Q(i,r | N - e_s)/(N - e_s)_r^alpha_r - Q(i,r | N)/N_r^alpha_r
 *
 * held fixed while an inner MVA fixed point (Core) is iterated, then refreshed
 * from the new queue lengths. Three refresh rounds are performed, as in Chandy
 * and Neuse's original; the npasses argument exists so that pfqn_scat can ask
 * for the single round that defines SCAT.
 *
 * Arithmetic: RUNTIME-GATED on the exponent, not compile-time gated.
 *
 * Two separate things could put this algorithm outside an exact field, and they
 * deserve separate answers. The first is N_r^alpha_r, a real power, which is
 * genuinely not a field operation for a general alpha; detail::num_pow_real
 * therefore accepts it exactly when alpha is a non-negative integer -- which
 * covers pfqn_linearizer, where alpha is pinned to 1, and the saturated
 * Gompertz exponent 2 that pfqn_linearizermx produces for any population past
 * about thirteen -- and REFUSES by name otherwise.
 *
 * The second is the stopping rule: the inner Core loop halts on
 * enorm(Q_{k+1} - Q_k) < tol, so what comes back is the iterate the stopping
 * rule selected, not the solution of a finite rational problem. That is a real
 * caveat but it is not a reason to deny the exact backend: an exact run returns
 * that iterate WITHOUT rounding error, which is precisely the quantity one
 * wants when asking how much of a double run's residual is arithmetic and how
 * much is the fixed point itself. Callers comparing the two backends must
 * compare like for like -- same tol, same maxiter, same warm start -- because
 * the two runs stop at different iterates otherwise.
 *
 * Scheduling. ForwardMVA in the reference uses the PS residence-time formula
 * W = D (1 + sum_s Q_1(.,s)) for EVERY discipline; the `type` argument is
 * accepted and carried but does not enter the recursion. The MATLAB source
 * documents why: the FCFS correction needs per-visit service times S = D/V,
 * and only chain-level demands D = V S are available here. That behaviour is
 * reproduced exactly rather than "improved", so that the port agrees with the
 * reference on FCFS models.
 *
 * One correction relative to the reference. An empty class (N_r == 0) makes
 * the MATLAB Update_Delta step evaluate Q/N_r^alpha_r = 0/0 and returns an
 * all-NaN solution. Empty classes are treated here as absent (zero queue
 * length, zero throughput, zero Delta), which is the same convention the
 * reference already applies inside pfqn_bs and inside its own Estimate step,
 * and which the MATLAB comment there calls "required, not cosmetic".
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_cntol.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of the Linearizer family, mirroring [Q,U,W,C,X,totiter]. */
template <class T>
struct LinearizerResult {
    Matrix<T> Q;        ///< (M x R) mean queue length
    Matrix<T> U;        ///< (M x R) utilization
    Matrix<T> W;        ///< (M x R) per-station residence time
    std::vector<T> C;   ///< (R) cycle time, N_r/X_r - Z_r
    std::vector<T> X;   ///< (R) per-class throughput
    int totiter;        ///< total inner iterations across all Core calls
};

namespace detail {

/**
 * base^e.
 *
 * With transcendental arithmetic this is std::pow. Without it -- the exact
 * rational backend -- a real power is not an operation of the field, but the
 * SPECIAL CASE that matters here is: the Gompertz exponent alpha is derived
 * from an integer population and, for any population past about thirteen,
 * saturates at exactly 2 in double, and pfqn_linearizer pins it at exactly 1.
 * Both are integers, and an integer power IS a field operation. So the exact
 * backend computes it by repeated multiplication when the exponent is a
 * non-negative integer and REFUSES otherwise, rather than silently rounding
 * the exponent, which would make the "exact" answer exact about the wrong
 * problem.
 */
template <class T>
T num_pow_real(const T& base, const T& e) {
    if constexpr (num_traits<T>::has_transcendental) {
        using std::pow;
        const T r = pow(base, e);
        return r;
    } else {
        const double ed = num_traits<T>::to_double(e);
        const double er = std::floor(ed + 0.5);
        if (er >= 0.0 && er < 1e6 && e == num_traits<T>::from_double(er))
            return num_pow_int(base, static_cast<unsigned>(er));
        throw UnsupportedError(
            "pfqn_egflinearizer: exact arithmetic cannot evaluate a non-integer real power; the "
            "scaling exponent is " + num_traits<T>::to_string(e) +
            ", which is not an integer, so this population needs the double or real backend");
    }
}

/**
 * Estimate step: the arrival-instant queue lengths at N_1 - e_s implied by the
 * queue lengths at N_1 and the frozen Delta. Returns Q1[s](i,r) for s in 1..R.
 */
template <class T>
std::vector<Matrix<T>> egflin_estimate(std::size_t M, std::size_t R, const std::vector<int>& N_1,
                                       const Matrix<T>& Q, const std::vector<Matrix<T>>& Delta,
                                       const std::vector<T>& alpha) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<Matrix<T>> Q1(R + 1, Matrix<T>(M, R, zero));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 1; s <= R; ++s) {
                const std::vector<int> Ns = oner(N_1, s);
                // zero-population guard rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                if (N_1[r] <= 0 || Ns[r] <= 0) {
                    Q1[s](i, r) = zero;
                } else {
                    const T na = num_pow_real(num_traits<T>::from_int(Ns[r]), alpha[r]);
                    const T n1a = num_pow_real(num_traits<T>::from_int(N_1[r]), alpha[r]);
                    Q1[s](i, r) = na * (Q(i, r) / n1a + Delta[r](i, s - 1));
                }
            }
    return Q1;
}

/** Forward MVA step: the PS residence-time formula, as in the reference. */
template <class T>
void egflin_forward_mva(const Matrix<T>& L, std::size_t M, std::size_t R,
                        const std::vector<int>& N_1, const std::vector<T>& Z,
                        const std::vector<Matrix<T>>& Q1, Matrix<T>& Q, Matrix<T>& W,
                        std::vector<T>& X) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            T acc = one;
            for (std::size_t s = 0; s < R; ++s) acc += Q1[r + 1](i, s);
            W(i, r) = L(i, r) * acc;
        }
    for (std::size_t r = 0; r < R; ++r) {
        T den = Z[r];
        for (std::size_t i = 0; i < M; ++i) den += W(i, r);
        if (N_1[r] <= 0) {
            X[r] = zero;
        } else {
            if (den == zero) throw NumericError("pfqn_egflinearizer: zero total residence time");
            X[r] = num_traits<T>::from_int(N_1[r]) / den;
        }
        for (std::size_t i = 0; i < M; ++i) Q(i, r) = X[r] * W(i, r);
    }
}

/** Core: iterate Estimate / ForwardMVA to the tolerance with Delta frozen. */
template <class T>
int egflin_core(const Matrix<T>& L, std::size_t M, std::size_t R, const std::vector<int>& N_1,
                const std::vector<T>& Z, Matrix<T>& Q, const std::vector<Matrix<T>>& Delta,
                const std::vector<T>& alpha, double tol, int maxiter, Matrix<T>& W,
                std::vector<T>& X, bool cntest = false) {
    int iter = 0;
    if (cntest) {
        // Chandy and Neuse (1982), p.129 and appendix: the cutoff is a function of the
        // population Core is running at, so it is recomputed here rather than once for the
        // whole Linearizer.
        tol = pfqn_cntol(N_1);
    }
    while (true) {
        const Matrix<T> Qlast = Q;
        const std::vector<Matrix<T>> Q1 = egflin_estimate(M, R, N_1, Q, Delta, alpha);
        egflin_forward_mva(L, M, R, N_1, Z, Q1, Q, W, X);
        double e;
        if (cntest) {
            // max_{i,r} |dQ(i,r)| / N_r over the non-empty classes; an empty class would
            // divide by zero and it carries no jobs to converge.
            e = 0.0;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    if (N_1[r] <= 0) continue;
                    const double d = std::fabs(num_traits<T>::to_double(T(Q(i, r) - Qlast(i, r)))) /
                                     static_cast<double>(N_1[r]);
                    if (d > e) e = d;
                }
        } else {
            e = enorm_diff(Q, Qlast);
        }
        const bool done = e < tol || iter > maxiter;
        ++iter;
        if (done) break;
    }
    return iter;
}

}  // namespace detail

/**
 * @param L       (M x R) service demands
 * @param N       (R) population per class
 * @param Z       (K x R) think times, summed over rows; may be empty
 * @param type    (M) scheduling discipline; accepted for interface parity,
 *                but the reference recursion is discipline-independent
 * @param tol     convergence tolerance on the Frobenius norm of dQ; NaN selects the
 *                published Linearizer termination test of Chandy and Neuse,
 *                Commun. ACM 25(2), 1982, p.129, under which each Core call stops
 *                when max_{i,r}|dQ(i,r)|/N_r falls below pfqn_cntol evaluated at
 *                the population Core is running at
 * @param maxiter total inner-iteration budget
 * @param alpha   (R) per-class scaling exponent
 * @param QN0     (M x R) warm start for the Bard-Schweitzer initialization
 * @param npasses number of Delta refresh rounds; 3 is the Chandy-Neuse fixed
 *                rule, pfqn_scat passes 1
 */
template <class T>
LinearizerResult<T> pfqn_egflinearizer(const Matrix<T>& L, const std::vector<int>& N,
                                       const Matrix<T>& Z,
                                       const std::vector<SchedStrategy>& type, double tol,
                                       int maxiter, const std::vector<T>& alpha,
                                       const Matrix<T>& QN0, int npasses = 3) {
    // runtime gating rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)

    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError(
            "pfqn_egflinearizer: demand matrix and population vector disagree on the class count");
    if (alpha.size() != R)
        throw InputError("pfqn_egflinearizer: alpha has the wrong class count");
    if (!type.empty() && type.size() != M)
        throw InputError("pfqn_egflinearizer: scheduling vector has the wrong station count");
    // The Chandy-Neuse cutoff is population-dependent, so it is recomputed inside each
    // egflin_core call rather than once here; tol stays NaN so that the pfqn_bs warm-start
    // inherits the same test. NaN is not caught by the positivity check below (every
    // comparison against NaN is false), which is exactly what the sentinel needs.
    const bool cntest = is_cntol(tol);
    if (!cntest && tol <= 0)
        throw InputError("pfqn_egflinearizer: tolerance must be positive");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_egflinearizer: negative population");

    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> Zs = sum_rows(Z, R);

    LinearizerResult<T> res;
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.W = Matrix<T>(M, R, zero);
    res.C.assign(R, zero);
    res.X.assign(R, zero);
    res.totiter = 0;

    // Delay-only model: every class is served entirely at the delay.
    bool anyDemand = false;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) != zero) anyDemand = true;
    if (M == 0 || !anyDemand) {
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == 0) continue;
            if (Zs[r] == zero)
                throw NumericError(
                    "pfqn_egflinearizer: a class has neither demand nor think time");
            res.X[r] = num_traits<T>::from_int(N[r]) / Zs[r];
            for (std::size_t i = 0; i < M; ++i) res.U(i, r) = res.X[r] * L(i, r);
        }
        return res;
    }

    // Initialize every population slice from Bard-Schweitzer.
    std::vector<Matrix<T>> Q(R + 1, Matrix<T>(M, R, zero));
    for (std::size_t s = 0; s <= R; ++s) {
        const std::vector<int> N_1 = oner(N, s);
        bool feasible = true;
        for (int v : N_1)
            if (v < 0) feasible = false;
        if (!feasible) continue;  // r == s with N_r == 0: the slice is unused
        std::vector<T> Nt(R, zero);
        for (std::size_t r = 0; r < R; ++r) Nt[r] = num_traits<T>::from_int(N_1[r]);
        // QN0 Bard-Schweitzer seed rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        const AmvaResult<T> b = pfqn_bs(L, Nt, Zs);
        Q[s] = b.QN;
    }

    // Delta[r](i,s) with s zero-based over 1..R of the MATLAB third index.
    std::vector<Matrix<T>> Delta(R, Matrix<T>(M, R, zero));

    Matrix<T> W(M, R, zero);
    std::vector<T> X(R, zero);
    for (int I = 0; I < npasses; ++I) {
        for (std::size_t s = 0; s <= R; ++s) {
            const std::vector<int> N_1 = oner(N, s);
            bool feasible = true;
            for (int v : N_1)
                if (v < 0) feasible = false;
            if (!feasible) continue;
            res.totiter += detail::egflin_core(L, M, R, N_1, Zs, Q[s], Delta, alpha, tol,
                                               maxiter - res.totiter, W, X, cntest);
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) {
                    // Absent class: no jobs at any population, so no correction.
                    for (std::size_t s = 0; s < R; ++s) Delta[r](i, s) = zero;
                    continue;
                }
                if (N[r] == 1) {
                    // At N - e_r only class r itself vanishes.
                    Q[r + 1](i, r) = zero;
                }
                const T nra = detail::num_pow_real(num_traits<T>::from_int(N[r]), alpha[r]);
                for (std::size_t s = 1; s <= R; ++s) {
                    const std::vector<int> Ns = oner(N, s);
                    if (Ns[r] > 0) {
                        const T nsa =
                            detail::num_pow_real(num_traits<T>::from_int(Ns[r]), alpha[r]);
                        Delta[r](i, s - 1) = Q[s](i, r) / nsa - Q[0](i, r) / nra;
                    } else {
                        // Chandy-Neuse 0/0 convention rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                        Delta[r](i, s - 1) = -Q[0](i, r) / nra;
                    }
                }
            }
    }

    res.totiter += detail::egflin_core(L, M, R, N, Zs, Q[0], Delta, alpha, tol,
                                       maxiter - res.totiter, W, X, cntest);
    res.Q = Q[0];
    res.W = W;
    res.X = X;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.U(i, r) = X[r] * L(i, r);
    for (std::size_t r = 0; r < R; ++r)
        res.C[r] = N[r] == 0 ? zero : num_traits<T>::from_int(N[r]) / X[r] - Zs[r];
    return res;
}

/** MATLAB defaults: tol = 1e-8, maxiter = 1000, no warm start. */
template <class T>
LinearizerResult<T> pfqn_egflinearizer(const Matrix<T>& L, const std::vector<int>& N,
                                       const Matrix<T>& Z, const std::vector<T>& alpha) {
    return pfqn_egflinearizer(L, N, Z, std::vector<SchedStrategy>(), 1e-8, 1000, alpha,
                              Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_EGFLINEARIZER_H
