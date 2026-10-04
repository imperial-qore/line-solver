/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LOOPING_H
#define LINE_API_PFQN_LOOPING_H

/**
 * @file
 * @ingroup api_pfqn
 * Eager Looping bounds for closed multiclass product-form networks.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_looping.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_looping.java. D. L. Eager,
 * "Bounding Algorithms for Queueing Network Models of Computer Systems",
 * Ph.D. thesis, Tech. Rept. CSRG-156, University of Toronto, 1984. Looping
 * supplies the initial pessimistic and optimistic estimates that the
 * multiple-class performance bound hierarchy starts from, so it carries a pair
 * of bounds rather than a single fixed point.
 *
 * A HEAP H_j is the class-j congestion that the current queue-length LOWER
 * BOUNDS have not yet accounted for; it is charged back at the pessimistic
 * inflation factor V_c = max_k D_ck or the optimistic one L_c = min_k D_ck,
 * which are the largest and smallest delays one customer can inflict. The
 * whole bracket rests on Q_jk(N - 1_c) being a lower bound, which is why the
 * queue lengths are seeded from Little's law at the station,
 *
 *   Q_jk(N - 1_c) = X_j(N - 1_c) R_jk(N - 1_c) >= n_j D_jk / (Z_j + U_j),
 *
 * with R_jk >= D_jk and U_j any UPPER bound on R_j(N - 1_c): the level-0 PBH
 * bound B_j, or the pessimistic R_j(N) of the current iterate, whichever is
 * smaller, since response time is nondecreasing in the population. Each
 * iterate lowers R^(pess), which raises the seed, which lowers R^(pess) again,
 * so the refinement is monotone and every iterate is a bound. It previously
 * refined the queue lengths through the convolution identity of Zahorjan
 * (1980), Q_jk(N - 1_c) = [X_j(N - 1_c)/X_j(N)] Q_jk(N), applied with one
 * class-level ratio at every station; that is not a per-station under-estimate,
 * so the queue lengths stopped being lower bounds, the heaps clamped to zero
 * and R^(opt) collapsed onto R^(pess). The level-0 multiple-class PBH bounds on
 * the mean response time are
 *
 *   J_j(n) = sum_k D_jk,   B_j(n) = sum_k D_jk + (sum(n) - 1) max_k D_jk,
 *
 * i.e. an arriving customer queues behind nobody, respectively behind every
 * other customer in the network at its own worst centre.
 *
 * Arithmetic: sums, products, divisions, maxima and minima only, so each
 * iterate is EXACT in rational arithmetic; the stopping rule selects which
 * iterate is returned.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_looping: the throughput bracket plus its queue lengths. */
template <class T>
struct LoopingBounds {
    std::vector<T> Xlo;  ///< (R) pessimistic (lower) throughput bound
    std::vector<T> Xup;  ///< (R) optimistic (upper) throughput bound
    Matrix<T> Q;         ///< (M x R) queue lengths on the pessimistic side
    Matrix<T> R;         ///< (M x R) residence times
    std::size_t iterations = 0;
    bool converged = false;
};

/**
 * @brief Eager Looping bounds for closed multiclass product-form networks.
 *
 * @param L (M x R) demands, @param N (R) populations,
 * @param Z (R) think times (empty for none)
 * @param tol convergence tolerance on the queue lengths
 * @param maxiter iteration cap
 */
template <class T>
LoopingBounds<T> pfqn_looping(const Matrix<T>& L, const std::vector<T>& N,
                              const std::vector<T>& Z, double tol = 1e-6,
                              std::size_t maxiter = 1000) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_looping: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_looping: Z has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    LoopingBounds<T> r;
    r.Xlo.assign(R, zero);
    r.Xup.assign(R, zero);
    r.Q = Matrix<T>(M, R, zero);
    r.R = Matrix<T>(M, R, zero);
    if (M == 0) return r;

    std::vector<T> Dtot(R, zero), Vpess(R, zero), Lopt(R, zero);
    T Ntot = zero;
    for (std::size_t c = 0; c < R; ++c) {
        T mx = L(0, c), mn = L(0, c);
        for (std::size_t k = 0; k < M; ++k) {
            Dtot[c] += L(k, c);
            if (L(k, c) > mx) mx = L(k, c);
            if (L(k, c) < mn) mn = L(k, c);
        }
        Vpess[c] = mx;
        Lopt[c] = mn;
        Ntot += N[c];
    }
    // level-0 multiple-class PBH bounds on R_j, at population N - 1_c
    std::vector<T> Jbnd(R, zero), Bm(R, zero);
    for (std::size_t c = 0; c < R; ++c) {
        Jbnd[c] = Dtot[c];
        const T span = (Ntot > two) ? T(Ntot - two) : zero;
        Bm[c] = T(Dtot[c] + span * Vpess[c]);
    }

    // Population of class j at the stations when one class-c job is removed, at
    // its optimistic and pessimistic extremes. These cap what the queue-length
    // lower bounds may account for; the remainder is the heap.
    Matrix<T> popt(R, R, zero), ppess(R, R, zero);   // (j, c)
    for (std::size_t c = 0; c < R; ++c)
        for (std::size_t j = 0; j < R; ++j) {
            const T nj = (c == j) ? T(N[j] - one) : N[j];
            if (nj <= zero) continue;
            const T zj = Z.empty() ? zero : Z[j];
            const T dj = T(zj + Jbnd[j]);
            if (dj > zero) popt(j, c) = T(T(Jbnd[j] / dj) * nj);
            const T db = T(zj + Bm[j]);
            if (db > zero) ppess(j, c) = T(T(Bm[j] / db) * nj);
        }

    // Zero is the one queue-length lower bound available before any iterate,
    // and the heaps then carry the whole population. The even split N/M this
    // once used is not a bound: it exceeds the true queue at every
    // below-average station, so the first residence is already not a bound.
    // Qm[c][j][k] = Q_jk(N - 1_c), a LOWER bound
    std::vector<std::vector<std::vector<T> > > Qm(
        R, std::vector<std::vector<T> >(R, std::vector<T>(M, zero)));
    Matrix<T> Hopt(R, R, zero), Hpess(R, R, zero);   // (j, c)
    for (std::size_t c = 0; c < R; ++c)
        for (std::size_t j = 0; j < R; ++j) {
            Hopt(j, c) = popt(j, c);
            Hpess(j, c) = ppess(j, c);
        }
    std::vector<T> Ub(Bm);                           // upper bound on R_j(N - 1_c)

    std::vector<T> Rc(R, zero), Rpess(R, zero), Ropt(R, zero);
    for (std::size_t it = 1; it <= maxiter; ++it) {
        r.iterations = it;
        const Matrix<T> Qprev = r.Q;

        for (std::size_t c = 0; c < R; ++c) {
            if (N[c] == zero) {
                for (std::size_t k = 0; k < M; ++k) r.R(k, c) = zero;
                r.Xlo[c] = zero;
                Rc[c] = zero;
                Rpess[c] = zero;
                continue;
            }
            T rsum = zero;
            for (std::size_t k = 0; k < M; ++k) {
                T qk = zero;
                for (std::size_t j = 0; j < R; ++j) qk += Qm[c][j][k];
                r.R(k, c) = L(k, c) * T(one + qk);
                rsum += r.R(k, c);
            }
            Rc[c] = rsum;
            T hp = zero;
            for (std::size_t j = 0; j < R; ++j) hp += Hpess(j, c);
            const T zc = Z.empty() ? zero : Z[c];
            Rpess[c] = T(Rc[c] + Vpess[c] * hp);
            const T den = T(zc + Rpess[c]);
            if (den == zero) throw NumericError("pfqn_looping: zero pessimistic cycle time");
            r.Xlo[c] = N[c] / den;
        }
        for (std::size_t c = 0; c < R; ++c) {
            if (N[c] == zero) {
                Ropt[c] = zero;
                continue;
            }
            const T zc = Z.empty() ? zero : Z[c];
            bool anySat = false;
            T sat = zero;
            for (std::size_t k = 0; k < M; ++k) {
                T used = zero;
                for (std::size_t j = 0; j < R; ++j)
                    if (j != c) used += r.Xlo[j] * L(k, j);
                const T den = T(one - used);
                if (den > zero) {
                    const T v = T(L(k, c) * N[c] / den - zc);
                    if (!anySat || v > sat) {
                        sat = v;
                        anySat = true;
                    }
                }
            }
            T ho = zero;
            for (std::size_t j = 0; j < R; ++j) ho += Hopt(j, c);
            T best = T(Rc[c] + Lopt[c] * ho);
            if (anySat && sat > best) best = sat;
            if (Dtot[c] > best) best = Dtot[c];
            // an optimistic bound can never exceed the pessimistic one
            Ropt[c] = (best < Rpess[c]) ? best : Rpess[c];
        }
        for (std::size_t c = 0; c < R; ++c)
            for (std::size_t k = 0; k < M; ++k) r.Q(k, c) = r.Xlo[c] * r.R(k, c);

        // R_j(N - 1_c) <= R_j(N) <= R_j^(pess): the pessimistic iterate tightens
        // the upper bound the seed divides by, and never loosens it.
        for (std::size_t j = 0; j < R; ++j)
            if (N[j] > zero && Rpess[j] > zero && Rpess[j] < Ub[j]) Ub[j] = Rpess[j];
        for (std::size_t c = 0; c < R; ++c)
            for (std::size_t j = 0; j < R; ++j) {
                const T nj = (c == j) ? T(N[j] - one) : N[j];
                const T zj = Z.empty() ? zero : Z[j];
                const T du = T(zj + Ub[j]);
                T qsum = zero;
                if (N[j] <= zero || nj <= zero || du <= zero) {
                    for (std::size_t k = 0; k < M; ++k) Qm[c][j][k] = zero;
                } else {
                    const T f = T(nj / du);
                    for (std::size_t k = 0; k < M; ++k) {
                        Qm[c][j][k] = f * L(k, j);
                        qsum += Qm[c][j][k];
                    }
                }
                // The clamp cannot fire: Ub >= J = Dtot, so qsum <= popt by
                // construction. It states the invariant the heaps rest on.
                const T ho = T(popt(j, c) - qsum);
                const T hp = T(ppess(j, c) - qsum);
                Hopt(j, c) = (ho > zero) ? ho : zero;
                Hpess(j, c) = (hp > zero) ? hp : zero;
            }

        bool anyNonEmpty = false;
        double maxdiff = 0.0;
        for (std::size_t c = 0; c < R; ++c) {
            if (N[c] <= zero) continue;
            anyNonEmpty = true;
            for (std::size_t k = 0; k < M; ++k) {
                const double d =
                    std::fabs(num_traits<T>::to_double(T(r.Q(k, c) - Qprev(k, c))));
                if (d > maxdiff) maxdiff = d;
            }
        }
        if (!anyNonEmpty || (it > 1 && maxdiff < tol)) {
            r.converged = true;
            break;
        }
    }

    for (std::size_t c = 0; c < R; ++c) {
        if (N[c] <= zero) continue;
        const T zc = Z.empty() ? zero : Z[c];
        const T den = T(zc + Ropt[c]);
        if (den == zero) throw NumericError("pfqn_looping: zero optimistic cycle time");
        r.Xup[c] = N[c] / den;
    }
    return r;
}

template <class T>
LoopingBounds<T> pfqn_looping(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_looping(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LOOPING_H
