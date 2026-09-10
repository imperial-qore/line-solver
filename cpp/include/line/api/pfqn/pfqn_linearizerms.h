/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LINEARIZERMS_H
#define LINE_API_PFQN_LINEARIZERMS_H

/**
 * Multiserver Linearizer (Krzesinski's Linearizer as described in Conway 1989,
 * with De Souza e Silva and Muntz's presentation of the marginal-probability
 * recursions).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_linearizerms.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/mva/Pfqn_linearizerms.java.
 *
 * Beyond the single-server Linearizer this carries, per station, the marginal
 * probabilities P(i,j) of finding j busy servers and the blocking probability
 * PB(i), estimated at the reduced populations by freezing them, and adds the
 * multiserver waiting term sum_{j<c-1} (c-1-j) P(i,j) to the residence time.
 *
 * Arithmetic: TRANSCENDENTAL-GATED. The inner Core loop stops on
 * norm(Q_{k+1} - Q_k) < tol, so the returned value depends on the stopping
 * rule and is not the solution of a finite rational problem. No transcendental
 * function is called.
 *
 * Convergence norm. MATLAB tests norm(Q - Qlast), the SPECTRAL norm of the
 * difference; this port tests the Frobenius norm, which dominates it. The
 * fixed point is identical and the test is if anything stricter, so no
 * solution accepted here would be rejected by the reference; the alternative
 * would be a singular value decomposition inside the inner loop of an
 * approximation, at every arithmetic, for no change in the answer.
 *
 * MATLAB-vs-JAR disagreement, resolved in MATLAB's favour. Both references
 * select the FCFS residence-time formula from the WHOLE type vector rather
 * than from the current station: MATLAB writes `if type == SchedStrategy.FCFS`
 * on an (M x 1) vector, which MATLAB evaluates as all(type == FCFS). The JAR
 * open-codes that test as `flag = true unless ANY station is FCFS` and then
 * takes the FCFS arm when `flag` holds -- exactly the opposite selection. This
 * port follows MATLAB: the FCFS arm is taken when every station is FCFS. Note
 * that the per-station form (type[i] == FCFS) would be the defensible reading
 * of the algorithm, but neither reference implements it and changing it here
 * would silently disagree with both.
 *
 * One correction relative to MATLAB. The Update_Delta and Estimate steps
 * divide by (N - e_s)_r and by N_r without guarding either against zero, so a
 * class with N_r == 1 produces Q/0 = Inf and then 0*Inf = NaN across the whole
 * solution, and an empty class produces 0/0. The guards of the single-server
 * pfqn_egflinearizer.m -- the Chandy and Neuse (1982) eq. (10) 0/0 convention
 * for (N - e_s)_r == 0, and an absent class contributing nothing -- are
 * applied here as well, since the two files implement the same correction and
 * only one of them carries the guard.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * ForwardMVA of the multiserver Linearizer. P and PB are read from the frozen
 * copies P_1, PB_1 (in the reference these are (M x c x 1+R) and (M x 1+R)
 * arrays whose entries are constant in the third index, so a plain copy of the
 * current P, PB is an exact stand-in) and overwritten with the new estimates.
 */
template <class T>
void linms_forward_mva(const Matrix<T>& L, std::size_t M, std::size_t R,
                       const std::vector<int>& N_1, const std::vector<T>& Z,
                       const std::vector<int>& nservers, bool allFCFS,
                       const std::vector<Matrix<T>>& Q1, const std::vector<Matrix<T>>& P1,
                       const std::vector<std::vector<T>>& PB1, Matrix<T>& Q, Matrix<T>& W,
                       std::vector<T>& X, Matrix<T>& P, std::vector<T>& PB) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < M; ++i) {
        const T c = num_traits<T>::from_int(nservers[i]);
        for (std::size_t r = 0; r < R; ++r) {
            W(i, r) = L(i, r) / c;
            // Zero demand: this class does not visit the station.
            if (L(i, r) == zero) continue;
            for (std::size_t s = 0; s < R; ++s)
                W(i, r) += (allFCFS ? L(i, s) : L(i, r)) / c * Q1[r + 1](i, s);
            // Partially-idle-server correction. It compensates the 1/m scaling of
            // the ARRIVING JOB'S OWN service, so it carries L(i,r)/m and no sum
            // over the other classes: at N = e_r no queueing is possible, and the
            // bracket must collapse to W = L(i,r) exactly. Summing L(i,s) over
            // every class instead (and omitting the /m) inflated it by a factor
            // that grew with both R and m -- 1.5x the demand at m=2, 3x at m=5.
            if (nservers[i] > 1)
                for (int j = 0; j <= nservers[i] - 2; ++j) {
                    const T wgt = num_traits<T>::from_int(nservers[i] - 1 - j);
                    W(i, r) += L(i, r) / c * wgt * P1[r + 1](i, static_cast<std::size_t>(j));
                }
        }
    }
    for (std::size_t r = 0; r < R; ++r) {
        T den = Z[r];
        for (std::size_t i = 0; i < M; ++i) den += W(i, r);
        if (N_1[r] <= 0) {
            X[r] = zero;
        } else {
            if (den == zero) throw NumericError("pfqn_linearizerms: zero total residence time");
            X[r] = num_traits<T>::from_int(N_1[r]) / den;
        }
        for (std::size_t i = 0; i < M; ++i) Q(i, r) = X[r] * W(i, r);
    }
    // Queue-length marginals. The relations
    //   p_j = (A p_{j-1} + d_{j-1})/j,  pB = (A (pB + p_{m-1}) + dB)/m,
    //   p_0 = 1 - pB - sum_j p_j
    // with A = sum_s X_s L_is the mean number of busy servers, are solved in
    // CLOSED FORM rather than iterated: as a Jacobi sweep they amplify by A per
    // pass and diverge once A approaches m, which is what made the corrected
    // residence times blow up (marginals reaching -13, throughput 64x too large)
    // on near-saturated models. Same fixed point wherever the iteration
    // converged, and unconditionally stable for any A < m.
    for (std::size_t i = 0; i < M; ++i) {
        const int ms = nservers[i];
        if (ms <= 1) continue;
        const std::size_t msz = static_cast<std::size_t>(ms);
        T A = zero, dB = zero;
        std::vector<T> d(msz, zero);
        for (std::size_t s = 0; s < R; ++s) {
            const T a_s = L(i, s) * X[s];
            A += a_s;
            for (std::size_t j = 0; j < msz; ++j)
                d[j] += a_s * (P1[s + 1](i, j) - P1[0](i, j));
            dB += a_s * (PB1[s + 1][i] - PB1[0][i]);
        }
        const T msT = num_traits<T>::from_int(ms);
        if (!(A < msT))
            throw NumericError(
                "pfqn_linearizerms: the station offers as many busy servers as it has, so the "
                "model is saturated and its queue-length marginals do not exist");
        std::vector<T> alpha(msz, zero), beta(msz, zero);
        alpha[0] = one;
        for (std::size_t j = 1; j < msz; ++j) {
            const T jT = num_traits<T>::from_int(static_cast<int>(j));
            alpha[j] = A * alpha[j - 1] / jT;
            beta[j] = (A * beta[j - 1] + d[j - 1]) / jT;
        }
        const T alphaB = A * alpha[msz - 1] / (msT - A);
        const T betaB = (A * beta[msz - 1] + dB + d[msz - 1]) / (msT - A);
        T num = one - betaB, den = one + alphaB;
        for (std::size_t j = 1; j < msz; ++j) {
            num -= beta[j];
            den += alpha[j];
        }
        for (std::size_t k = 0; k < P.cols(); ++k) P(i, k) = zero;
        P(i, 0) = num / den;
        for (std::size_t j = 1; j < msz; ++j) P(i, j) = alpha[j] * P(i, 0) + beta[j];
        PB[i] = alphaB * P(i, 0) + betaB;
    }
}

/**
 * Estimate step of the multiserver Linearizer.
 *
 * It returns the queue lengths AND the marginals at each reduced population.
 * The marginals used to be a plain copy of the current P, PB at population N,
 * on the reading that the reference's (M x c x 1+R) arrays are constant in the
 * third index. They are not, once DeltaP/DeltaPB exist: mixing queue lengths
 * reduced to N - e_s with marginals still at N breaks the identity
 * Q + sum_{j<=m-2} (m-1-j) p_j >= m-1 that guarantees W >= D, which is exactly
 * how the residence time came out below the mean service time.
 */
template <class T>
void linms_estimate(std::size_t M, std::size_t R, const std::vector<int>& N_1, const Matrix<T>& Q,
                    const Matrix<T>& P, const std::vector<T>& PB,
                    const std::vector<Matrix<T>>& Delta, const std::vector<Matrix<T>>& DeltaP,
                    const Matrix<T>& DeltaPB, const std::vector<int>& nservers,
                    std::vector<Matrix<T>>& Q1, std::vector<Matrix<T>>& P1,
                    std::vector<std::vector<T>>& PB1) {
    const T zero = num_traits<T>::from_int(0);
    Q1.assign(R + 1, Matrix<T>(M, R, zero));
    P1.assign(R + 1, Matrix<T>(M, P.cols(), zero));
    PB1.assign(R + 1, std::vector<T>(M, zero));
    for (std::size_t i = 0; i < M; ++i) {
        if (nservers[i] > 1) {
            for (std::size_t j = 0; j < P.cols(); ++j) {
                P1[0](i, j) = P(i, j);
                for (std::size_t s = 1; s <= R; ++s)
                    P1[s](i, j) = P(i, j) + DeltaP[s - 1](i, j);
            }
            PB1[0][i] = PB[i];
            for (std::size_t s = 1; s <= R; ++s) PB1[s][i] = PB[i] + DeltaPB(i, s - 1);
        }
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
    }
}

template <class T>
int linms_core(const Matrix<T>& L, std::size_t M, std::size_t R, const std::vector<int>& N_1,
               const std::vector<T>& Z, const std::vector<int>& nservers, bool allFCFS,
               Matrix<T>& Q, Matrix<T>& P, std::vector<T>& PB,
               const std::vector<Matrix<T>>& Delta, const std::vector<Matrix<T>>& DeltaP,
               const Matrix<T>& DeltaPB, double tol, int maxiter, Matrix<T>& W,
               std::vector<T>& X) {
    int iter = 0;
    std::vector<Matrix<T>> Q1, P1;
    std::vector<std::vector<T>> PB1;
    while (true) {
        ++iter;
        const Matrix<T> Qlast = Q;
        linms_estimate(M, R, N_1, Q, P, PB, Delta, DeltaP, DeltaPB, nservers, Q1, P1, PB1);
        linms_forward_mva(L, M, R, N_1, Z, nservers, allFCFS, Q1, P1, PB1, Q, W, X, P, PB);
        const double e = enorm_diff(Q, Qlast);
        if (e < tol || iter > maxiter) break;
    }
    return iter;
}

/** Initial marginal probabilities from the aggregate queue lengths. */
template <class T>
void linms_init_marginals(std::size_t M, std::size_t R, const std::vector<int>& nservers,
                          const std::vector<Matrix<T>>& Q, const std::vector<int>& N,
                          std::vector<Matrix<T>>& P, std::vector<std::vector<T>>& PB) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t s = 0; s <= R; ++s) {
        const std::vector<int> N_1 = oner(N, s);
        long pop = 0;
        for (int v : N_1) pop += v;
        if (pop <= 0) {
            // Empty network: the station is idle with probability one. Skipping
            // left p_0 at zero, i.e. a marginal that is not a distribution.
            for (std::size_t i = 0; i < M; ++i) {
                if (nservers[i] <= 1) continue;
                for (std::size_t k = 0; k < P[s].cols(); ++k) P[s](i, k) = zero;
                P[s](i, 0) = one;
                PB[s][i] = zero;
            }
            continue;
        }
        const T popT = num_traits<T>::from_int(pop);
        const T pop1T = num_traits<T>::from_int(pop + 1);
        for (std::size_t i = 0; i < M; ++i) {
            if (nservers[i] <= 1) continue;
            T qsum = zero;
            for (std::size_t r = 0; r < R; ++r) qsum += Q[s](i, r);
            const T two = num_traits<T>::from_int(2);
            for (int j = 1; j <= nservers[i] - 1; ++j)
                P[s](i, static_cast<std::size_t>(j)) = two * qsum / (popT * pop1T);
            // Fewer jobs than servers: they cannot all be busy, so pB is 0 rather
            // than a division by the vanishing slack. MATLAB guards it as
            // `pop > nservers(i)-1`; throwing here refused models the reference
            // solves, and every population below the server count is one of them.
            if (pop > nservers[i] - 1) {
                const T slack = num_traits<T>::from_int(pop + 1 - nservers[i]);
                PB[s][i] = two * qsum / slack / (popT * pop1T);
            } else {
                PB[s][i] = zero;
            }
            T p0 = one - PB[s][i];
            for (int j = 1; j <= nservers[i] - 1; ++j) p0 -= P[s](i, static_cast<std::size_t>(j));
            P[s](i, 0) = p0;
        }
    }
}

}  // namespace detail

/**
 * @param L        (M x R) service demands
 * @param N        (R) population per class
 * @param Z        (K x R) think times, summed over rows; may be empty
 * @param nservers (M) number of servers per station, at least one
 * @param type     (M) scheduling discipline; the FCFS arm is taken only if
 *                 every station is FCFS, as in MATLAB
 * @param tol      convergence tolerance
 * @param maxiter  total inner-iteration budget
 * @param QN0      (M x R) warm start for the Bard-Schweitzer initialization
 */
template <class T>
LinearizerResult<T> pfqn_linearizerms(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z, const std::vector<int>& nservers,
                                      const std::vector<SchedStrategy>& type, double tol,
                                      int maxiter, const Matrix<T>& QN0) {
    // field-arithmetic rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)

    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError(
            "pfqn_linearizerms: demand matrix and population vector disagree on the class count");
    if (nservers.size() != M)
        throw InputError("pfqn_linearizerms: server-count vector has the wrong station count");
    for (int c : nservers)
        if (c < 1) throw InputError("pfqn_linearizerms: server count below one");
    if (!type.empty() && type.size() != M)
        throw InputError("pfqn_linearizerms: scheduling vector has the wrong station count");
    if (tol <= 0) throw InputError("pfqn_linearizerms: tolerance must be positive");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_linearizerms: negative population");

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

    bool allFCFS = !type.empty();
    for (std::size_t i = 0; i < type.size(); ++i)
        if (type[i] != SchedStrategy::FCFS) allFCFS = false;

    std::size_t cmax = 1;
    for (int c : nservers) cmax = static_cast<std::size_t>(c) > cmax ? c : cmax;

    Matrix<T> Zm(1, R, zero);
    for (std::size_t r = 0; r < R; ++r) Zm(0, r) = Zs[r];

    std::vector<Matrix<T>> Q(R + 1, Matrix<T>(M, R, zero));
    for (std::size_t s = 0; s <= R; ++s) {
        const std::vector<int> N_1 = oner(N, s);
        bool feasible = true;
        for (int v : N_1)
            if (v < 0) feasible = false;
        if (!feasible) continue;
        std::vector<T> Nt(R, zero);
        for (std::size_t r = 0; r < R; ++r) Nt[r] = num_traits<T>::from_int(N_1[r]);
        // pfqn_bs three-argument seed rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        const AmvaResult<T> b = pfqn_bs(L, Nt, Zs);
        Q[s] = b.QN;
    }

    std::vector<Matrix<T>> P(R + 1, Matrix<T>(M, cmax, zero));
    std::vector<std::vector<T>> PB(R + 1, std::vector<T>(M, zero));
    detail::linms_init_marginals(M, R, nservers, Q, N, P, PB);

    std::vector<Matrix<T>> Delta(R, Matrix<T>(M, R, zero));
    // Linearizer corrections for the MARGINALS. Without them the marginals stay
    // at population N while the queue lengths are reduced to N - e_s, which
    // breaks Q + sum_j (m-1-j) p_j >= m-1 and lets W fall below the mean service
    // time. Probabilities do not scale with the population, so the analogue of
    // Delta is a plain difference rather than a per-job rate.
    std::vector<Matrix<T>> DeltaP(R, Matrix<T>(M, cmax, zero));
    Matrix<T> DeltaPB(M, R, zero);
    Matrix<T> W(M, R, zero);
    std::vector<T> X(R, zero);

    for (int I = 0; I < 2; ++I) {
        for (std::size_t s = 0; s <= R; ++s) {
            const std::vector<int> N_1 = oner(N, s);
            bool feasible = true;
            for (int v : N_1)
                if (v < 0) feasible = false;
            if (!feasible) continue;
            res.totiter +=
                detail::linms_core(L, M, R, N_1, Zs, nservers, allFCFS, Q[s], P[s], PB[s], Delta,
                                   DeltaP, DeltaPB, tol, maxiter - res.totiter, W, X);
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) {
                    for (std::size_t s = 0; s < R; ++s) Delta[r](i, s) = zero;
                    continue;
                }
                const T nrT = num_traits<T>::from_int(N[r]);
                for (std::size_t s = 1; s <= R; ++s) {
                    const std::vector<int> Ns = oner(N, s);
                    if (Ns[r] > 0) {
                        Delta[r](i, s - 1) =
                            Q[s](i, r) / num_traits<T>::from_int(Ns[r]) - Q[0](i, r) / nrT;
                    } else {
                        Delta[r](i, s - 1) = -Q[0](i, r) / nrT;
                    }
                }
            }
        // Update_DeltaP: the marginals' counterpart of Update_Delta.
        for (std::size_t i = 0; i < M; ++i) {
            if (nservers[i] <= 1) continue;
            for (std::size_t s = 1; s <= R; ++s) {
                for (std::size_t j = 0; j < static_cast<std::size_t>(nservers[i]); ++j)
                    DeltaP[s - 1](i, j) = P[s](i, j) - P[0](i, j);
                DeltaPB(i, s - 1) = PB[s][i] - PB[0][i];
            }
        }
    }

    res.totiter += detail::linms_core(L, M, R, N, Zs, nservers, allFCFS, Q[0], P[0], PB[0], Delta,
                                      DeltaP, DeltaPB, tol, maxiter - res.totiter, W, X);
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

/** MATLAB defaults: all stations PS, tol = 1e-8, maxiter = 1000, no warm start. */
template <class T>
LinearizerResult<T> pfqn_linearizerms(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z, const std::vector<int>& nservers) {
    return pfqn_linearizerms(L, N, Z, nservers, std::vector<SchedStrategy>(), 1e-8, 1000,
                             Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LINEARIZERMS_H
