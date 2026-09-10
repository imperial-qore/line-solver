/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_MODULATE_H
#define LINE_API_MAM_MMAP_MODULATE_H

/**
 * Modulate a family of marked MAPs by an environment chain.
 *
 * Templated port of matlab/lib/m3a/m3a/mmap/mmap_modulate.m and
 * mmap_mixture_order2.m.
 *
 * `mmap_modulate(P, HT, MMAPs)` builds the marked arrival process of a system
 * that sits in environment j for a phase-type holding time HT[j], arriving
 * according to MMAPs[j] while it does, and then jumps to environment i with
 * probability P(j,i). The state is the pair (holding-time phase, arrival phase),
 * so environment j contributes a block of order nh(j) nm(j) and the result has
 * order sum_j nh(j) nm(j):
 *
 *   diagonal block j: krons(HT0[j], A0[j]), the two clocks running together;
 *   off-diagonal (j,i): P(j,i) (HT1[j] (x) I) 1 (pie(HT[i]) (x) pie(MMAP[i])),
 *     the holding time expiring and the new environment being entered in its
 *     own initial phase.
 *
 * THE ENVIRONMENT SWITCH EMITS NO ARRIVAL. The off-diagonal blocks appear in D0
 * and in NO Dc, which is what makes the jump a hidden transition. The reference
 * writes the per-class off-diagonal blocks explicitly as `0*P(j,i)*...`, i.e. a
 * zero of the right shape; that is reproduced by simply leaving them zero.
 *
 * A PLAIN MAP IS AUTO-CONVERTED to a one-class MMAP, as the reference does, so a
 * caller may pass either. All the components must then agree on the class count.
 *
 * `mmap_mixture_order2` is the second-order companion: given m^2 two-phase
 * components PHs(i,j) and a second-order transition tensor P2, it builds the
 * marked MAP whose state records the PREVIOUS and the current class, which is
 * what lets a mixture reproduce a lag-1 class correlation that an order-1
 * mixture cannot.
 *
 * ARITHMETIC: field. Kronecker products and block assembly only.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/trace/mtrace_cross_moment.h"
#include "line/api/trace/mtrace_sigma2.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_stats.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param P     (J x J) environment transition probabilities
 * @param HT    the J phase-type holding times, as (D0, D1) pairs
 * @param comps the J marked arrival processes, one per environment
 */
template <class T>
Mmap<T> mmap_modulate(const Matrix<T>& P, const std::vector<Map<T>>& HT,
                      const std::vector<Mmap<T>>& comps) {
    const std::size_t J = HT.size();
    if (comps.size() != J)
        throw InputError("mmap_modulate: the holding-time and MMAP lists must have equal length");
    if (J == 0) throw InputError("mmap_modulate: no environments given");
    if (P.rows() != J || P.cols() != J)
        throw InputError("mmap_modulate: P must be square of the environment count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // A plain MAP arrives here as an Mmap with no per-class blocks; the
    // reference promotes it to a one-class MMAP by copying D1.
    std::vector<Mmap<T>> M = comps;
    for (std::size_t j = 0; j < J; ++j)
        if (M[j].Dc.empty()) M[j].Dc.push_back(M[j].D1);
    const std::size_t K = M[0].classes();
    for (std::size_t j = 1; j < J; ++j)
        if (M[j].classes() != K)
            throw InputError("mmap_modulate: the MMAPs must have the same number of types");

    std::vector<std::size_t> nh(J), nm(J), nq(J), off(J + 1, 0);
    for (std::size_t j = 0; j < J; ++j) {
        nh[j] = HT[j].D0.rows();
        nm[j] = M[j].order();
        nq[j] = nh[j] * nm[j];
        off[j + 1] = off[j] + nq[j];
    }
    const std::size_t N = off[J];

    Mmap<T> out;
    out.D0 = Matrix<T>(N, N, zero);
    out.D1 = Matrix<T>(N, N, zero);
    out.Dc.assign(K, Matrix<T>(N, N, zero));

    // The entry law of each environment: pie(HT) (x) pie(MMAP).
    std::vector<std::vector<T>> entry(J);
    for (std::size_t i = 0; i < J; ++i) {
        const std::vector<T> a = map_pie(HT[i]);
        const std::vector<T> b = map_pie(M[i].map());
        entry[i].assign(nq[i], zero);
        for (std::size_t p = 0; p < nh[i]; ++p)
            for (std::size_t q = 0; q < nm[i]; ++q) entry[i][p * nm[i] + q] = T(a[p] * b[q]);
    }

    for (std::size_t j = 0; j < J; ++j) {
        // Diagonal block: the two clocks run together, krons(HT0, A0).
        const Matrix<T> diag = krons(HT[j].D0, M[j].D0);
        for (std::size_t r = 0; r < nq[j]; ++r)
            for (std::size_t c = 0; c < nq[j]; ++c) out.D0(off[j] + r, off[j] + c) = diag(r, c);
        // Per-class blocks: an arrival of class k leaves the holding-time phase
        // untouched, I (x) Ak.
        for (std::size_t k = 0; k < K; ++k) {
            const Matrix<T> blk = kron(eye<T>(nh[j]), M[j].Dc[k]);
            for (std::size_t r = 0; r < nq[j]; ++r)
                for (std::size_t c = 0; c < nq[j]; ++c)
                    out.Dc[k](off[j] + r, off[j] + c) = blk(r, c);
        }

        // Off-diagonal: the holding time expires and environment i is entered.
        // (HT1 (x) I) 1 is the exit rate out of each (phase, phase) pair.
        const Matrix<T> exit = kron(HT[j].D1, eye<T>(nm[j]));
        std::vector<T> exitRow(nq[j], zero);
        for (std::size_t r = 0; r < nq[j]; ++r) {
            T s = zero;
            for (std::size_t c = 0; c < nq[j]; ++c) s += exit(r, c);
            exitRow[r] = s;
        }
        for (std::size_t i = 0; i < J; ++i) {
            if (i == j) continue;
            for (std::size_t r = 0; r < nq[j]; ++r)
                for (std::size_t c = 0; c < nq[i]; ++c)
                    out.D0(off[j] + r, off[i] + c) = T(P(j, i) * exitRow[r] * entry[i][c]);
        }
    }
    (void)one;
    return mmap_normalize(out);
}

/**
 * Second-order mixture: the state records the previous and the current class.
 *
 * @param PHs (m x m) two-phase components, PHs[i][j] being the sojourn in class
 *            j reached from class i
 * @param P2  (m x m) second-order class transition probabilities
 */
template <class T>
Mmap<T> mmap_mixture_order2(const std::vector<std::vector<Map<T>>>& PHs, const Matrix<T>& P2) {
    const std::size_t m = PHs.size();
    if (m == 0) throw InputError("mmap_mixture_order2: no components given");
    for (std::size_t i = 0; i < m; ++i)
        if (PHs[i].size() != m)
            throw InputError("mmap_mixture_order2: the component table must be square");
    if (P2.rows() != m || P2.cols() != m)
        throw InputError("mmap_mixture_order2: P2 must be square of the component count");
    // The reference indexes in blocks of two, so every component is order two.
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            if (PHs[i][j].D0.rows() != 2)
                throw InputError("mmap_mixture_order2: every component must be of order two");

    const T zero = num_traits<T>::from_int(0);
    const std::size_t N = 2 * m * m;
    Mmap<T> out;
    out.D0 = Matrix<T>(N, N, zero);
    out.D1 = Matrix<T>(N, N, zero);
    out.Dc.assign(m, Matrix<T>(N, N, zero));

    // Diagonal: the sojourn generator of each (previous, current) pair.
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            const std::size_t k = i * m + j, base = 2 * k;
            for (std::size_t r = 0; r < 2; ++r)
                for (std::size_t c = 0; c < 2; ++c) out.D0(base + r, base + c) = PHs[i][j].D0(r, c);
        }

    // Off-diagonal: a completion in pair (i1,j1) marks class j1 and moves to a
    // pair whose PREVIOUS class is j1, with probability P2(i1, i2).
    for (std::size_t i1 = 0; i1 < m; ++i1)
        for (std::size_t j1 = 0; j1 < m; ++j1) {
            const std::size_t k1 = i1 * m + j1;
            for (std::size_t i2 = 0; i2 < m; ++i2)
                for (std::size_t j2 = 0; j2 < m; ++j2) {
                    if (j1 != i2) continue;
                    const std::size_t k2 = i2 * m + j2;
                    const std::vector<T> pie = map_pie(PHs[i2][j2]);
                    for (std::size_t r = 0; r < 2; ++r) {
                        T ex = zero;
                        for (std::size_t c = 0; c < 2; ++c) ex += -PHs[i1][j1].D0(r, c);
                        for (std::size_t c = 0; c < 2; ++c) {
                            const T v = T(P2(i1, i2) * ex * pie[c]);
                            out.Dc[j1](2 * k1 + r, 2 * k2 + c) = v;
                        }
                    }
                }
        }
    return mmap_normalize(out);
}


/**
 * Second-order mixture FITTED from cross moments and a triple sigma, the
 * reference's `mmap_mixture_fit`.
 *
 * Port of matlab/lib/m3a/m3a/mmap/mmap_mixture_fit.m and
 * mmap_mixture_fit_trace.m. Each ordered class pair (i,j) gets its own APH(2)
 * fitted to the cross moments M1(i,j), M2(i,j), M3(i,j) -- the sojourn in class
 * j when it followed class i -- and those are assembled exactly as
 * `mmap_mixture_order2` assembles them, except that the transition weight is
 * the CONDITIONAL second-order probability
 *
 *   p = P2(i1, i2, j2) / sum_h P2(i1, i2, h),
 *
 * so the chain over (previous, current) pairs is stochastic by construction.
 * That normalization is the whole difference from `mmap_mixture_order2`, which
 * takes an already-conditioned two-index weight.
 *
 * @param P2 the triple sigma, (C x C*C) with entry (i, j*C + h) = P2(i,j,h)
 * @param M1,M2,M3 the (C x C) cross moments
 */
template <class T>
Mmap<T> mmap_mixture_fit(const Matrix<T>& P2, const Matrix<T>& M1, const Matrix<T>& M2,
                         const Matrix<T>& M3) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmap_mixture_fit fits an APH(2) per class pair");
    const std::size_t m = M1.rows();
    if (m == 0) throw InputError("mmap_mixture_fit: no classes given");
    if (M1.cols() != m || M2.rows() != m || M2.cols() != m || M3.rows() != m || M3.cols() != m)
        throw InputError("mmap_mixture_fit: the cross-moment tables must be square and equal");
    if (P2.rows() != m || P2.cols() != m * m)
        throw InputError("mmap_mixture_fit: P2 must be (C x C*C), the flattened triple sigma");
    const T zero = num_traits<T>::from_int(0);

    // A pair that was never observed has no cross moment to fit, and fitting a
    // zero first moment would silently produce a degenerate component. The
    // reference does not guard it; refusing by name says which pair is missing.
    std::vector<std::vector<Map<T>>> PHs(m, std::vector<Map<T>>(m));
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            if (!(num_traits<T>::to_double(M1(i, j)) > 0.0))
                throw InputError(
                    "mmap_mixture_fit: the cross moment of class pair (" + std::to_string(i + 1) +
                    ", " + std::to_string(j + 1) +
                    ") is not positive, so that pair was never observed and its component cannot "
                    "be fitted. Supply a trace in which every ordered class pair occurs, or fit "
                    "fewer classes");
            PHs[i][j] = aph2_fit(M1(i, j), M2(i, j), M3(i, j)).aph;
        }

    const std::size_t N = 2 * m * m;
    Mmap<T> out;
    out.D0 = Matrix<T>(N, N, zero);
    out.D1 = Matrix<T>(N, N, zero);
    out.Dc.assign(m, Matrix<T>(N, N, zero));

    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            const std::size_t base = 2 * (i * m + j);
            for (std::size_t r = 0; r < 2; ++r)
                for (std::size_t c = 0; c < 2; ++c) out.D0(base + r, base + c) = PHs[i][j].D0(r, c);
        }

    for (std::size_t i1 = 0; i1 < m; ++i1)
        for (std::size_t j1 = 0; j1 < m; ++j1) {
            const std::size_t k1 = i1 * m + j1;
            // The row of the conditional law: P2(i1, j1, .) normalized.
            T rowtot = zero;
            for (std::size_t h = 0; h < m; ++h) rowtot += P2(i1, j1 * m + h);
            for (std::size_t i2 = 0; i2 < m; ++i2)
                for (std::size_t j2 = 0; j2 < m; ++j2) {
                    if (j1 != i2) continue;
                    if (!(num_traits<T>::to_double(rowtot) > 0.0)) continue;
                    const T pr = T(P2(i1, i2 * m + j2) / rowtot);
                    const std::size_t k2 = i2 * m + j2;
                    const std::vector<T> pie = map_pie(PHs[i2][j2]);
                    for (std::size_t r = 0; r < 2; ++r) {
                        T ex = zero;
                        for (std::size_t c = 0; c < 2; ++c) ex += -PHs[i1][j1].D0(r, c);
                        for (std::size_t c = 0; c < 2; ++c)
                            out.Dc[j1](2 * k1 + r, 2 * k2 + c) = T(pr * ex * pie[c]);
                    }
                }
        }
    return mmap_normalize(out);
}

/**
 * `mmap_mixture_fit` driven from a marked trace: the triple sigma and the cross
 * moments are measured on the trace itself.
 *
 * @param Tv the inter-arrival times
 * @param A  the class of each arrival
 */
template <class T>
Mmap<T> mmap_mixture_fit_trace(const std::vector<T>& Tv, const std::vector<int>& A) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError("mmap_mixture_fit_trace: the trace and its labels must agree in length");
    const Matrix<T> P2 = trace::mtrace_sigma2<T>(A);
    const Matrix<T> M1 = trace::mtrace_cross_moment(Tv, A, 1u).mc;
    const Matrix<T> M2 = trace::mtrace_cross_moment(Tv, A, 2u).mc;
    const Matrix<T> M3 = trace::mtrace_cross_moment(Tv, A, 3u).mc;
    return mmap_mixture_fit(P2, M1, M2, M3);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_MODULATE_H
