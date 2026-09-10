/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_MWRBB_H
#define LINE_API_PFQN_PFQN_MWRBB_H

/**
 * Majumdar-Woodside robust box bounds on the per-class throughput of a closed
 * multiclass network with mixed scheduling disciplines (Perf. Eval. 32 (1998)
 * 101-136).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mwrbb.m together with its three
 * local functions mwrbb_denom, mwrbb_station_wrest and mwrbb_residence. The
 * upper bound intersects the no-contention bound (eq. 2) with the
 * utilization bound (eq. 3); the lower bound is the throughput guarantee of
 * Theorem 2 (eq. 15), whose per-visit queueing delay depends on the discipline
 * at the station: FIFO (Theorem 1 / Lemma 1), processor sharing (Lemma 2),
 * preemptive priority (Lemma 3) and non-preemptive priority (Lemmas 4-5). The
 * coupled inequalities are resolved by interval narrowing.
 *
 * The bounds are distribution-insensitive (NBUE service only) and routing
 * insensitive: mean visits, mean demands, populations, think times,
 * disciplines and priorities are the whole input.
 *
 * ARITHMETIC. Only sums, products, minima and divisions appear, so the bounds
 * are EXACT in rational arithmetic and are deliberately left ungated. The
 * fixed point is a monotone narrowing whose stopping rule (1e-13 absolute on
 * both bound vectors) is a double constant converted into T, so a higher
 * precision instantiation stops at the same place, not further -- the
 * remaining slack there is the bound's, not the iteration's.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Station discipline codes, matching the MATLAB `sched` argument. */
enum class MwrbbSched { Fifo = 0, Ps = 1, PrioNonPreemptive = 2, PrioPreemptive = 3, Aba = 4 };

/** Return value of pfqn_mwrbb, mirroring [Xlo, Xup, Wlo]. */
template <class T>
struct MwrbbBounds {
    std::vector<T> Xlo;
    std::vector<T> Xup;
    Matrix<T> Wlo;
};

namespace detail {

/**
 * Per-visit residence at station k for class c EXCLUDING the isolated
 * higher-priority 1/f_c term (MATLAB mwrbb_station_wrest).
 */
template <class T>
T mwrbb_station_wrest(std::size_t k, std::size_t c, const Matrix<T>& V, const Matrix<T>& S,
                      const std::vector<T>& N, const std::vector<T>& fup, const T& fc,
                      const std::vector<MwrbbSched>& sched, const std::vector<int>& prio) {
    const std::size_t C = V.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T Vkc = V(k, c), Skc = S(k, c);
    const T scale = T(fc * Vkc);
    const MwrbbSched d = sched[k];

    if (d == MwrbbSched::Fifo) {
        T s = zero;
        for (std::size_t m = 0; m < C; ++m) {
            T pcm = one;
            if (scale != zero) {
                pcm = T(T(fup[m] * V(k, m)) / scale);
                if (pcm > one) pcm = one;
            }
            s += N[m] * S(k, m) * pcm;
        }
        return s;  // the m = c term is the tagged job's own service
    }
    if (d == MwrbbSched::Ps) {
        T dp = zero;
        for (std::size_t m = 0; m < C; ++m) {
            T Ncont = N[m];
            if (m == c) Ncont = T(N[c] - one);
            T term = Skc;
            if (scale != zero) {
                const T alt = T(T(fup[m] * V(k, m) * S(k, m)) / scale);
                if (alt < term) term = alt;
            }
            dp += Ncont * term;
        }
        return T(Skc + dp);
    }
    if (d == MwrbbSched::Aba) {
        T s = zero;
        for (std::size_t m = 0; m < C; ++m) s += N[m] * S(k, m);
        return s;
    }

    // Preemptive (3) or non-preemptive (2) priority.
    T dp = zero;
    for (std::size_t m = 0; m < C; ++m) {
        if (prio[m] != prio[c]) continue;  // higher priority goes through Bh
        T Ncont = N[m];
        if (m == c) Ncont = T(N[c] - one);
        T pcm = one;
        if (scale != zero) {
            pcm = T(T(fup[m] * V(k, m)) / scale);
            if (pcm > one) pcm = one;
        }
        dp += Ncont * S(k, m) * pcm;
    }
    if (d == MwrbbSched::PrioNonPreemptive) {
        // Water-filling over the lower-priority classes, longest demand first.
        std::vector<std::size_t> lower;
        for (std::size_t m = 0; m < C; ++m)
            if (prio[m] > prio[c]) lower.push_back(m);
        std::sort(lower.begin(), lower.end(),
                  [&](std::size_t a, std::size_t b) { return S(k, b) < S(k, a); });
        T budget = one;
        for (std::size_t idx = 0; idx < lower.size(); ++idx) {
            const std::size_t l = lower[idx];
            T al = zero;
            if (N[l] > zero) {
                al = T(budget / N[l]);
                if (scale != zero) {
                    const T capr = T(T(fup[l] * V(k, l)) / scale);
                    if (capr < al) al = capr;
                }
            }
            if (al < zero) al = zero;
            dp += N[l] * al * S(k, l);
            budget = T(budget - N[l] * al);
            if (budget < zero) budget = zero;
        }
    }
    return T(Skc + dp);
}

}  // namespace detail

/**
 * @param V     (K x C) mean visits
 * @param S     (K x C) mean demand per visit
 * @param N     (C) population
 * @param Z     (C) think time, empty for zero
 * @param sched (K) per-station discipline, empty for all FIFO
 * @param prio  (C) class priority, lower value = higher priority; empty for equal
 */
template <class T>
MwrbbBounds<T> pfqn_mwrbb(const Matrix<T>& V, const Matrix<T>& S, const std::vector<T>& N,
                          const std::vector<T>& Z, const std::vector<MwrbbSched>& sched,
                          const std::vector<int>& prio) {
    const std::size_t K = V.rows(), C = V.cols();
    if (S.rows() != K || S.cols() != C) throw InputError("pfqn_mwrbb: V and S have different shapes");
    if (N.size() != C) throw InputError("pfqn_mwrbb: N has the wrong length");
    if (!Z.empty() && Z.size() != C) throw InputError("pfqn_mwrbb: Z has the wrong length");
    if (!sched.empty() && sched.size() != K) throw InputError("pfqn_mwrbb: sched has the wrong length");
    if (!prio.empty() && prio.size() != C) throw InputError("pfqn_mwrbb: prio has the wrong length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::vector<T> Zv = Z.empty() ? std::vector<T>(C, zero) : Z;
    const std::vector<MwrbbSched> sc =
        sched.empty() ? std::vector<MwrbbSched>(K, MwrbbSched::Fifo) : sched;
    const std::vector<int> pr = prio.empty() ? std::vector<int>(C, 0) : prio;

    // No-contention upper bound on the cycle rate f_c = X_c/N_c (eqs. 1-2).
    std::vector<T> fup(C, zero), flo(C, zero);
    for (std::size_t c = 0; c < C; ++c) {
        T s = Zv[c];
        for (std::size_t k = 0; k < K; ++k) s += V(k, c) * S(k, c);
        if (s == zero) throw InputError("pfqn_mwrbb: a class has no demand and no think time");
        fup[c] = T(one / s);
    }

    const T tol = num_traits<T>::from_double(1e-13);
    for (int it = 0; it < 20000; ++it) {
        const std::vector<T> fup_old = fup, flo_old = flo;

        // Utilization-based narrowing of the upper bounds (eq. 3).
        for (std::size_t c = 0; c < C; ++c) {
            T cap = fup[c];
            for (std::size_t k = 0; k < K; ++k) {
                T other = zero;
                for (std::size_t m = 0; m < C; ++m)
                    if (m != c) other += N[m] * V(k, m) * S(k, m) * flo[m];
                const T denomk = T(N[c] * V(k, c) * S(k, c));
                if (denomk > zero) {
                    const T v = T(T(one - other) / denomk);
                    if (v < cap) cap = v;
                }
            }
            if (cap < zero) cap = zero;
            if (cap < fup[c]) fup[c] = cap;
        }

        // lower-bound narrowing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        for (std::size_t c = 0; c < C; ++c) {
            T DEN = Zv[c], Bh = zero;
            const T fc = flo[c];
            for (std::size_t k = 0; k < K; ++k) {
                if (V(k, c) == zero) continue;
                if (sc[k] == MwrbbSched::PrioPreemptive || sc[k] == MwrbbSched::PrioNonPreemptive) {
                    for (std::size_t m = 0; m < C; ++m)
                        if (pr[m] < pr[c]) Bh += N[m] * fup[m] * V(k, m) * S(k, m);
                }
                DEN += V(k, c) * detail::mwrbb_station_wrest(k, c, V, S, N, fup, fc, sc, pr);
            }
            if (DEN == zero) throw NumericError("pfqn_mwrbb: zero cycle time in the lower bound");
            T val = T(T(one - Bh) / DEN);
            if (val < zero) val = zero;
            if (val > flo[c]) flo[c] = val;
        }

        T du = zero, dl = zero;
        for (std::size_t c = 0; c < C; ++c) {
            const T a = num_abs(T(fup[c] - fup_old[c]));
            const T b = num_abs(T(flo[c] - flo_old[c]));
            if (a > du) du = a;
            if (b > dl) dl = b;
        }
        if (du < tol && dl < tol) break;
    }

    MwrbbBounds<T> r;
    r.Xlo.resize(C);
    r.Xup.resize(C);
    for (std::size_t c = 0; c < C; ++c) {
        r.Xlo[c] = T(N[c] * flo[c]);
        r.Xup[c] = T(N[c] * fup[c]);
    }
    r.Wlo = Matrix<T>(K, C, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t k = 0; k < K; ++k) {
            if (V(k, c) == zero) continue;
            T W = detail::mwrbb_station_wrest(k, c, V, S, N, fup, flo[c], sc, pr);
            const T scale = T(flo[c] * V(k, c));
            if ((sc[k] == MwrbbSched::PrioNonPreemptive || sc[k] == MwrbbSched::PrioPreemptive) &&
                scale > zero) {
                for (std::size_t m = 0; m < C; ++m)
                    if (pr[m] < pr[c]) W += T(N[m] * fup[m] * V(k, m) * S(k, m) / scale);
            }
            r.Wlo(k, c) = W;
        }
    return r;
}

template <class T>
MwrbbBounds<T> pfqn_mwrbb(const Matrix<T>& V, const Matrix<T>& S, const std::vector<T>& N,
                          const std::vector<T>& Z) {
    return pfqn_mwrbb(V, S, N, Z, std::vector<MwrbbSched>(), std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_MWRBB_H
