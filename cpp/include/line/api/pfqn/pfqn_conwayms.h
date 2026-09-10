/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CONWAYMS_H
#define LINE_API_PFQN_CONWAYMS_H

/**
 * Conway's multiserver Linearizer for chain-dependent FCFS queues (Conway
 * 1989, "Fast Approximate Solution of Queueing Networks with Multi-Server
 * Chain-Dependent FCFS Queues").
 *
 * Templated port of matlab/src/api/pfqn/pfqn_conwayms.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_conwayms.java.
 *
 * The distinguishing feature over pfqn_linearizerms is the pair of conditional
 * service rates at a c-server station, obtained by averaging over the
 * compositions n of the c busy servers among the R chains,
 *
 *   A_i(n)    = multinomial(n) prod_c F_r(i,c)^{n_c},   F_r(i,c) = T_1(c|r) L(i,c) / sum_c' ...
 *   XR(i,r)   = sum_{n in B_r} A_i(n) / (sum_c n_c / L(i,c))  /  sum_{n in B_r} A_i(n)
 *   XE(i,r,c) = the same restricted to n_c >= 1
 *
 * with B_r = { n : sum(n) = c, n <= N - e_r }, which enter the residence time
 * as W = L + PB XR + sum_c XE (Q_1 - L T_1).
 *
 * Arithmetic: TRANSCENDENTAL-GATED, on the fixed-point tolerance alone. The
 * inner Core loop stops on norm(Q_{k+1} - Q_k) < tol, so the returned value
 * depends on the stopping rule. The reference forms A_i(n) as
 * exp(multinomialln(n) + n log F), but that is a convenience: the value is a
 * finite rational in F, and this port computes it as multinomial(n) times an
 * integer power product, which is both exact and free of the 0 * (-Inf) = NaN
 * the reference produces whenever some F_r(i,c) vanishes at n_c == 0.
 *
 * Convergence norm and the FCFS selection follow pfqn_linearizerms: Frobenius
 * rather than spectral (dominating, same fixed point), and the FCFS arm of the
 * single-server residence time is taken only when every station is FCFS, which
 * is what MATLAB's `if type == SchedStrategy.FCFS` on a vector means. The JAR
 * inverts this test; MATLAB is the reference.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/api/pfqn/pfqn_linearizerms.h"  // shares the marginal-probability initialization
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Estimate step: frozen-Delta queue lengths and the auxiliary throughputs. */
template <class T>
void conway_estimate(std::size_t M, std::size_t R, const std::vector<int>& N_1,
                     const Matrix<T>& Q, const std::vector<Matrix<T>>& Delta, const Matrix<T>& W,
                     std::vector<Matrix<T>>& Q1, Matrix<T>& T_1) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 1; s <= R; ++s) {
                const std::vector<int> Ns = oner(N_1, s);
                if (N_1[r] <= 0 || Ns[r] <= 0) {
                    Q1[s](i, r) = zero;
                } else {
                    Q1[s](i, r) = num_traits<T>::from_int(Ns[r]) *
                                  (Q(i, r) / num_traits<T>::from_int(N_1[r]) + Delta[r](i, s - 1));
                }
            }
    // T_1 is Little's law over the queueing part of the cycle, sum_i Q1 / sum_i W,
    // and not the ratio at the FIRST station with a positive residence time: the
    // per-station estimates disagree, so picking one made the answer depend on the
    // station order. The demand matrix carries no order, so a model symmetric under
    // permuting classes and stations together must return equal class throughputs,
    // and with the single-station pick it did not.
    for (std::size_t s = 0; s < R; ++s)
        for (std::size_t r = 0; r < R; ++r) T_1(s, r) = zero;
    for (std::size_t r = 1; r <= R; ++r) {
        const std::vector<int> Nr = oner(N_1, r);
        for (std::size_t s = 0; s < R; ++s) {
            if (N_1[s] <= 0 || Nr[s] <= 0) continue;
            T num = zero;
            T den = zero;
            for (std::size_t i = 0; i < M; ++i) {
                if (W(i, s) > zero) {
                    // Delta is indexed [queued class](station, removed class), as the
                    // Q1 loop above uses it: here class s queues and class r-1 is removed
                    num += num_traits<T>::from_int(Nr[s]) *
                           (Q(i, s) / num_traits<T>::from_int(N_1[s]) + Delta[s](i, r - 1));
                    den += W(i, s);
                }
            }
            // a reduced-population throughput cannot be negative; a negative one
            // makes log(F) complex in the XR/XE sums of the forward step
            if (den > zero) {
                const T t1 = num / den;
                T_1(s, r - 1) = (t1 < zero) ? zero : t1;
            }
        }
    }
}

/**
 * The conditional service rates XR(i,r) and XE(i,r,c) at the multiserver
 * stations, by enumeration of the compositions of the server count.
 */
template <class T>
void conway_conditional_rates(const Matrix<T>& L, std::size_t M, std::size_t R,
                              const std::vector<int>& N_1, const std::vector<int>& nservers,
                              const Matrix<T>& T_1, std::vector<T>& XRflat,
                              std::vector<T>& XEflat) {
    const T zero = num_traits<T>::from_int(0);
    XRflat.assign(M * R, zero);
    XEflat.assign(M * R * R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        // F_r(i,c), the chain-composition weights at station i.
        Matrix<T> F(M, R, zero);
        for (std::size_t i = 0; i < M; ++i) {
            T den = zero;
            for (std::size_t c = 0; c < R; ++c) den += L(i, c) * T_1(c, r);
            if (den == zero) continue;  // station unreachable at this population
            for (std::size_t c = 0; c < R; ++c) F(i, c) = T_1(c, r) * L(i, c) / den;
        }
        const std::vector<int> Nr = oner(N_1, r + 1);
        for (std::size_t i = 0; i < M; ++i) {
            if (nservers[i] <= 1) continue;
            std::vector<int> n(R, 0);
            first_composition(n, nservers[i]);
            T Csum = zero, XRacc = zero;
            std::vector<T> Cx(R, zero), XEacc(R, zero);
            bool more = true;
            while (more) {
                bool inB = true;
                for (std::size_t c = 0; c < R; ++c)
                    if (n[c] > Nr[c]) inB = false;
                if (inB) {
                    // A_i(n) = multinomial(n) prod_{c: n_c > 0} F(i,c)^{n_c}.
                    T Ai = num_multinomial<T>(n);
                    for (std::size_t c = 0; c < R; ++c)
                        if (n[c] > 0) Ai *= num_pow_int(F(i, c), static_cast<unsigned>(n[c]));
                    // sum_c n_c mu(i,c) with mu = 1/L, skipping the absent chains
                    // so that a zero demand never contributes 0 * Inf.
                    T rate = zero;
                    for (std::size_t c = 0; c < R; ++c) {
                        if (n[c] == 0) continue;
                        if (L(i, c) == zero)
                            throw NumericError(
                                "pfqn_conwayms: a chain occupies a server at a station where it "
                                "has zero demand, so its service rate is unbounded");
                        rate += num_traits<T>::from_int(n[c]) / L(i, c);
                    }
                    if (rate != zero) {
                        Csum += Ai;
                        XRacc += Ai / rate;
                        for (std::size_t c = 0; c < R; ++c)
                            if (n[c] >= 1) {
                                Cx[c] += Ai;
                                XEacc[c] += Ai / rate;
                            }
                    }
                }
                more = next_composition(n);
            }
            if (Csum != zero) XRflat[i * R + r] = XRacc / Csum;
            for (std::size_t c = 0; c < R; ++c)
                if (Cx[c] != zero) XEflat[(i * R + r) * R + c] = XEacc[c] / Cx[c];
        }
    }
}

template <class T>
void conway_forward_mva(const Matrix<T>& L, std::size_t M, std::size_t R,
                        const std::vector<int>& N_1, const std::vector<T>& Z,
                        const std::vector<int>& nservers, bool allFCFS,
                        const std::vector<Matrix<T>>& Q1, const Matrix<T>& P_1,
                        const std::vector<T>& PB_1, const Matrix<T>& T_1, Matrix<T>& Q,
                        Matrix<T>& W, std::vector<T>& X, Matrix<T>& P, std::vector<T>& PB) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> XR, XE;
    conway_conditional_rates(L, M, R, N_1, nservers, T_1, XR, XE);

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            if (nservers[i] == 1) {
                W(i, r) = L(i, r);
                for (std::size_t c = 0; c < R; ++c)
                    W(i, r) += (allFCFS ? L(i, c) : L(i, r)) * Q1[r + 1](i, c);
            } else {
                W(i, r) = L(i, r) + PB_1[i] * XR[i * R + r];
                for (std::size_t c = 0; c < R; ++c)
                    W(i, r) += XE[(i * R + r) * R + c] * (Q1[r + 1](i, c) - L(i, c) * T_1(c, r));
            }
        }
    for (std::size_t r = 0; r < R; ++r) {
        T den = Z[r];
        for (std::size_t i = 0; i < M; ++i) den += W(i, r);
        if (N_1[r] <= 0) {
            X[r] = zero;
        } else {
            if (den == zero) throw NumericError("pfqn_conwayms: zero total residence time");
            X[r] = num_traits<T>::from_int(N_1[r]) / den;
        }
        for (std::size_t i = 0; i < M; ++i) Q(i, r) = X[r] * W(i, r);
    }
    // Queue-length marginals. The relations
    //   p_j = A*p_{j-1}/j,  pB = A*(pB + p_{ms-1})/ms,  p_0 = 1 - pB - sum_j p_j
    // with A = sum_s X_s*L_is the mean number of busy servers are solved in closed
    // form rather than iterated. As a Jacobi iteration they amplify by A per sweep,
    // and since the convergence test watches Q and W but not P the routine returned
    // marginals whose mass had run to 334 behind the p_0 = max(0,1-...) floor.
    // conway_estimate hands the same marginals to every reduced population, so the
    // population corrections that pfqn_linearizerms carries here are all zero.
    for (std::size_t i = 0; i < M; ++i) {
        if (nservers[i] <= 1) continue;
        const std::size_t ms = static_cast<std::size_t>(nservers[i]);
        T A = zero;
        for (std::size_t s = 0; s < R; ++s) A += L(i, s) * X[s];
        for (std::size_t k = 0; k < P.cols(); ++k) P(i, k) = zero;
        if (!(A < num_traits<T>::from_int(nservers[i]))) {
            // Saturated: the closed form is singular and its limit is the degenerate
            // marginal, every server busy with probability one. N = m with Z = 0
            // reaches it exactly, so this is a legal input.
            PB[i] = one;
            continue;
        }
        std::vector<T> alpha(ms, zero);
        alpha[0] = one;
        T sumAlpha = zero;
        for (std::size_t j = 1; j < ms; ++j) {
            alpha[j] = A * alpha[j - 1] / num_traits<T>::from_int(static_cast<int>(j));
            sumAlpha += alpha[j];
        }
        const T alphaB = A * alpha[ms - 1] / (num_traits<T>::from_int(nservers[i]) - A);
        const T p0 = one / (one + sumAlpha + alphaB);
        P(i, 0) = p0;
        for (std::size_t j = 1; j < ms; ++j) P(i, j) = alpha[j] * p0;
        PB[i] = alphaB * p0;
    }
}

template <class T>
int conway_core(const Matrix<T>& L, std::size_t M, std::size_t R, const std::vector<int>& N_1,
                const std::vector<T>& Z, const std::vector<int>& nservers, bool allFCFS,
                Matrix<T>& Q, Matrix<T>& P, std::vector<T>& PB,
                const std::vector<Matrix<T>>& Delta, double tol, int maxiter, Matrix<T>& W,
                std::vector<T>& X) {
    const T zero = num_traits<T>::from_int(0);
    W = L;  // the reference seeds the residence times with the demands
    Matrix<T> Wlast;
    bool haveWlast = false;
    int iter = 1;
    while (true) {
        const Matrix<T> Qlast = Q;
        const Matrix<T> P_1 = P;
        const std::vector<T> PB_1 = PB;
        std::vector<Matrix<T>> Q1(R + 1, Matrix<T>(M, R, zero));
        Matrix<T> T_1(R, R, zero);
        conway_estimate(M, R, N_1, Q, Delta, W, Q1, T_1);
        conway_forward_mva(L, M, R, N_1, Z, nservers, allFCFS, Q1, P_1, PB_1, T_1, Q, W, X, P, PB);
        // W must enter the test: Q alone is satisfied on the FIRST sweep whenever Q
        // cannot move (M=1 seeds Q at its own fixed point), and the residence times
        // returned then are still the seed W=L, so T_1=Q/W is unbounded and the
        // throughput exceeds the station's own service capacity.
        double e = std::numeric_limits<double>::infinity();
        if (haveWlast) e = std::max(enorm_diff(Q, Qlast), enorm_diff(W, Wlast));
        Wlast = W;
        haveWlast = true;
        const bool done = e < tol || iter > maxiter;
        ++iter;
        if (done) break;
    }
    return iter;
}

}  // namespace detail

/**
 * @param L        (M x R) service demands
 * @param N        (R) population per class
 * @param Z        (K x R) think times, summed over rows; may be empty
 * @param nservers (M) number of servers per station, at least one
 * @param type     (M) scheduling discipline; empty means all-FCFS, the MATLAB
 *                 default for this routine
 * @param tol      convergence tolerance
 * @param maxiter  total inner-iteration budget
 * @param QN0      (M x R) warm start; empty for the default N/M
 */
template <class T>
LinearizerResult<T> pfqn_conwayms(const Matrix<T>& L, const std::vector<int>& N,
                                  const Matrix<T>& Z, const std::vector<int>& nservers,
                                  const std::vector<SchedStrategy>& type, double tol, int maxiter,
                                  const Matrix<T>& QN0) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_conwayms requires transcendental arithmetic");

    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError(
            "pfqn_conwayms: demand matrix and population vector disagree on the class count");
    if (nservers.size() != M)
        throw InputError("pfqn_conwayms: server-count vector has the wrong station count");
    for (int c : nservers)
        if (c < 1) throw InputError("pfqn_conwayms: server count below one");
    if (!type.empty() && type.size() != M)
        throw InputError("pfqn_conwayms: scheduling vector has the wrong station count");
    if (tol <= 0) throw InputError("pfqn_conwayms: tolerance must be positive");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_conwayms: negative population");

    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> Zs = sum_rows(Z, R);

    LinearizerResult<T> res;
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.W = Matrix<T>(M, R, zero);
    res.C.assign(R, zero);
    res.X.assign(R, zero);
    res.totiter = 0;
    if (M == 0) return res;

    // default scheduling rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    bool allFCFS = true;
    for (std::size_t i = 0; i < type.size(); ++i)
        if (type[i] != SchedStrategy::FCFS) allFCFS = false;

    std::size_t cmax = 1;
    for (int c : nservers) cmax = static_cast<std::size_t>(c) > cmax ? c : cmax;

    // Initial queue lengths: N_1(r)/M, or the supplied warm start.
    std::vector<Matrix<T>> Q(R + 1, Matrix<T>(M, R, zero));
    const T mT = num_traits<T>::from_int(static_cast<long>(M));
    for (std::size_t s = 0; s <= R; ++s) {
        const std::vector<int> N_1 = oner(N, s);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                Q[s](i, r) = QN0.empty() ? num_traits<T>::from_int(N_1[r]) / mT : QN0(i, r);
    }
    if (!QN0.empty() && (QN0.rows() != M || QN0.cols() != R))
        throw InputError("pfqn_conwayms: initial queue lengths have the wrong shape");

    std::vector<Matrix<T>> P(R + 1, Matrix<T>(M, cmax, zero));
    std::vector<std::vector<T>> PB(R + 1, std::vector<T>(M, zero));
    detail::linms_init_marginals(M, R, nservers, Q, N, P, PB);

    std::vector<Matrix<T>> Delta(R, Matrix<T>(M, R, zero));
    Matrix<T> W(M, R, zero);
    std::vector<T> X(R, zero);

    for (int I = 0; I < 2; ++I) {
        for (std::size_t s = 0; s <= R; ++s) {
            const std::vector<int> N_1 = oner(N, s);
            bool feasible = true;
            for (int v : N_1)
                if (v < 0) feasible = false;
            if (!feasible) continue;
            res.totiter += detail::conway_core(L, M, R, N_1, Zs, nservers, allFCFS, Q[s], P[s],
                                               PB[s], Delta, tol, maxiter - res.totiter, W, X);
        }
        // The reference refreshes Delta only for classes with N_s > 2.
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) continue;
                const T nrT = num_traits<T>::from_int(N[r]);
                for (std::size_t s = 1; s <= R; ++s) {
                    if (N[s - 1] <= 2) continue;
                    const std::vector<int> Ns = oner(N, s);
                    if (Ns[r] > 0)
                        Delta[r](i, s - 1) =
                            Q[s](i, r) / num_traits<T>::from_int(Ns[r]) - Q[0](i, r) / nrT;
                    else
                        Delta[r](i, s - 1) = -Q[0](i, r) / nrT;
                }
            }
    }

    // The reference passes the full maxiter, not the remaining budget, to the
    // final Core call; that is reproduced here.
    res.totiter +=
        detail::conway_core(L, M, R, N, Zs, nservers, allFCFS, Q[0], P[0], PB[0], Delta, tol,
                            maxiter, W, X);
    res.Q = Q[0];
    res.W = W;
    res.X = X;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            res.U(i, r) = nservers[i] == 1
                              ? T(X[r] * L(i, r))
                              : T(X[r] * L(i, r) / num_traits<T>::from_int(nservers[i]));
    for (std::size_t r = 0; r < R; ++r)
        res.C[r] = N[r] == 0 ? zero : num_traits<T>::from_int(N[r]) / X[r] - Zs[r];
    return res;
}

/** MATLAB defaults: all stations FCFS, tol = 1e-8, maxiter = 1000. */
template <class T>
LinearizerResult<T> pfqn_conwayms(const Matrix<T>& L, const std::vector<int>& N,
                                  const Matrix<T>& Z, const std::vector<int>& nservers) {
    return pfqn_conwayms(L, N, Z, nservers, std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CONWAYMS_H
