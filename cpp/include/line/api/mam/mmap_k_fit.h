/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_K_FIT_H
#define LINE_API_MAM_MMAP_K_FIT_H

/**
 * EXACT inverses of the class-marking map: solve for the split directly instead
 * of optimizing it.
 *
 * Templated port of matlab/lib/m3a/m3a/mamap2m/mmap2k_fit.m and mmap3k_fit.m.
 *
 * `mamap2m_fit_fb_multiclass` treats the split as a quadratic program because
 * the per-class targets are generally unreachable. When the number of free
 * marking parameters EQUALS the number of characteristics being matched, the
 * system is square and can simply be solved; these two functions do that, and
 * report whether the solution landed inside the unit box.
 *
 *  - `mmap2k_fit` inverts an AMAP(2)'s three flows against (p, F, B) per class
 *    in closed form, one class at a time, and falls back to
 *    `mamap2m_fit_gamma_fb` when no AMAP(2) form gives a feasible split.
 *  - `mmap3k_fit` does it for a MAP of any order by building the marking system
 *    M q = y explicitly: row (a,b) of M is pie (-D0)^-a Dc (-D0)^-b e, which is
 *    the characteristic that row matches, and the target vector y carries
 *    p, p F, p B and p B2 for the orders (1,0), (1,1), (2,0) and (3,0).
 *
 * EXACT IS REPORTED, NOT ASSUMED. Both return a flag saying whether the solved
 * split was feasible; when it was not, the entries are clamped into [0,1] and
 * the flag is false, so a caller can tell an exact fit from a projected one.
 * Silently clamping and calling the result exact is the failure mode these
 * functions exist to avoid.
 *
 * THE DEGENERATE LOCUS IS REFUSED, not worked around. `mmap2k_fit` skips an
 * AMAP form whose denominators vanish; `mmap3k_fit` refuses outright when the
 * marking matrix is singular, because there the characteristics do not
 * determine the split and any answer would be arbitrary.
 *
 * ARITHMETIC: transcendental, through amap2_fit_gamma and the linear solves.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/mamap2m_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** A marked MAP and whether the marking system was solved exactly. */
template <class T>
struct MmapKFitResult {
    Mmap<T> mmap;
    bool exact = false;
};

namespace kfitdetail {

/** The closed-form marking inverse of one class, for the two canonical forms. */
template <class T>
bool marking_inverse(int form, const T& h1, const T& h2, const T& r1, const T& r2, const T& p,
                     const T& Fc, const T& Bc, double degentol, T* q1, T* q2, T* q3) {
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    auto small = [&](const T& v) { return std::fabs(num_traits<T>::to_double(v)) < degentol; };

    if (form == 1) {
        const T d1 = T((r2 - one) * (r1 - one) * (h1 * r2 - h2));
        const T d2 = T((r2 - one) * r1);
        const T d3 = T(r1 * r2 * (h1 + h2 * r1 - h2));
        const T e1 = T(h1 + h2 * r1 - h2), e2 = T(h1 * r2 - h2);
        if (small(d1) || small(d2) || small(d3) || small(e1) || small(e2)) return false;
        const T W = T(r1 * r2 - r2 + one);
        *q1 = T(p * W * ((h1 * r2 - h1 - h2) + Bc) / d1);
        *q2 = T(p * W *
                ((h1 * h1 * (r2 - one) + h1 * h2 * r1 * (r2 - one) - h2 * h2 * r1) / (e1 * e2) -
                 Fc / e1 + Bc / e2) /
                d2);
        *q3 = T(p * W * ((h1 + h2 * r1) - Fc) / d3);
    } else {
        const T U = T(h1 * r1 * r2 - h1 * r1 + h1 - h2);
        const T d1 = T((r2 - one) * (r1 - one) * U);
        const T e1 = T(h1 + h2 * r1 - h2);
        const T d2 = T((r2 - one) * e1);
        if (small(d1) || small(d2) || small(r2) || small(U) || small(e1)) return false;
        const T V = T(r1 * r2 - r1 - r2 + two);
        *q1 = T(p * V * ((h1 * r1 * r2 - h1 * r1 - h2) + Bc) / d1);
        *q2 = T(p * V * (h2 - Fc) / d2);
        *q3 = T(p * V *
                ((h1 * h1 + h1 * h2 * r1 * r2 - h2 * h2) / (e1 * U) - Fc / e1 - Bc / U) / r2);
    }
    return true;
}

/** The characteristic orders (a, b) the marking system matches, by MAP order. */
inline std::vector<std::pair<unsigned, unsigned>> marking_orders(std::size_t n, std::size_t z) {
    std::vector<std::pair<unsigned, unsigned>> o;
    if (n == 2) {
        o.push_back(std::make_pair(1u, 0u));
        o.push_back(std::make_pair(1u, 1u));
        o.push_back(std::make_pair(2u, 0u));
    } else if (n == 3) {
        o.push_back(std::make_pair(1u, 0u));
        o.push_back(std::make_pair(1u, 1u));
        o.push_back(std::make_pair(2u, 0u));
        o.push_back(std::make_pair(3u, 0u));
    } else {
        o.push_back(std::make_pair(1u, 0u));
        o.push_back(std::make_pair(1u, 1u));
        unsigned a = 2;
        while (o.size() < n + 1) o.push_back(std::make_pair(a++, 0u));
    }
    if (o.size() > z) o.resize(z);
    return o;
}

}  // namespace kfitdetail

/**
 * Exact marking of an AMAP(2) against per-class (p, F, B).
 *
 * @param M1,M2,M3 the first three moments of the inter-arrival time
 * @param GAMMA    the autocorrelation decay rate
 * @param P,F,B    per-class probabilities, forward and backward moments
 */
template <class T>
MmapKFitResult<T> mmap2k_fit(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                             const std::vector<T>& P, const std::vector<T>& F,
                             const std::vector<T>& B) {
    static_assert(num_traits<T>::has_transcendental, "mmap2k_fit inverts a moment system");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = P.size();
    if (K == 0) throw InputError("mmap2k_fit: no classes given");
    if (F.size() != K || B.size() != K)
        throw InputError("mmap2k_fit: P, F and B must have the same length");
    const double degentol = 1e-8, feastol = 1e-8;

    const Amap2FitGammaResult<T> a = amap2_fit_gamma(M1, M2, M3, GAMMA);
    bool haveBest = false;
    double bestViol = 0.0;
    Map<T> bestMap;
    int bestForm = 1;
    std::vector<std::vector<T>> bestQ;

    for (std::size_t j = 0; j < a.amaps.size(); ++j) {
        const Map<T>& mp = a.amaps[j];
        if (mp.order() != 2) continue;
        if (!(num_abs(T(mp.D0(1, 0))) <= zero)) continue;
        int form;
        if (mp.D1(0, 1) == zero) form = 1;
        else if (mp.D1(0, 0) == zero) form = 2;
        else continue;

        const T h1 = T(-one / mp.D0(0, 0)), h2 = T(-one / mp.D0(1, 1));
        const T r1 = T(mp.D0(0, 1) * h1), r2 = T(mp.D1(1, 1) * h2);

        std::vector<std::vector<T>> q(3, std::vector<T>(K, zero));
        bool ok = true;
        for (std::size_t c = 0; c < K && ok; ++c)
            ok = kfitdetail::marking_inverse(form, h1, h2, r1, r2, P[c], F[c], B[c], degentol,
                                             &q[0][c], &q[1][c], &q[2][c]);
        if (!ok) continue;

        // Feasibility: inside the unit box, and each flow's split summing to one.
        double viol = 0.0;
        for (std::size_t jj = 0; jj < 3; ++jj) {
            T s = zero;
            for (std::size_t c = 0; c < K; ++c) {
                const double v = num_traits<T>::to_double(q[jj][c]);
                viol = std::max(viol, std::max(0.0, -v));
                viol = std::max(viol, std::max(0.0, v - 1.0));
                s += q[jj][c];
            }
            viol = std::max(viol, std::fabs(num_traits<T>::to_double(s) - 1.0));
        }
        if (!haveBest || viol < bestViol) {
            haveBest = true;
            bestViol = viol;
            bestMap = mp;
            bestForm = form;
            bestQ = q;
        }
    }

    MmapKFitResult<T> out;
    if (haveBest && bestViol <= feastol) {
        out.mmap.D0 = bestMap.D0;
        out.mmap.D1 = bestMap.D1;
        out.mmap.Dc.assign(K, Matrix<T>(2, 2, zero));
        for (std::size_t c = 0; c < K; ++c) {
            for (std::size_t jj = 0; jj < 3; ++jj) {
                if (bestQ[jj][c] < zero) bestQ[jj][c] = zero;
                if (bestQ[jj][c] > one) bestQ[jj][c] = one;
            }
            if (bestForm == 1) {
                out.mmap.Dc[c](0, 0) = T(bestMap.D1(0, 0) * bestQ[0][c]);
                out.mmap.Dc[c](1, 0) = T(bestMap.D1(1, 0) * bestQ[1][c]);
                out.mmap.Dc[c](1, 1) = T(bestMap.D1(1, 1) * bestQ[2][c]);
            } else {
                out.mmap.Dc[c](0, 1) = T(bestMap.D1(0, 1) * bestQ[0][c]);
                out.mmap.Dc[c](1, 0) = T(bestMap.D1(1, 0) * bestQ[1][c]);
                out.mmap.Dc[c](1, 1) = T(bestMap.D1(1, 1) * bestQ[2][c]);
            }
        }
        out.exact = true;
        return out;
    }
    out.mmap = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);
    out.exact = false;
    return out;
}

/**
 * Exact marking of an arbitrary MAP by solving the marking system directly.
 *
 * @param D0,D1 the underlying MAP
 * @param P,F,B per-class probabilities, forward and backward moments
 * @param B2    per-class second-order backward moments; required once the MAP
 *              has more than three arrival entries to mark
 */
template <class T>
MmapKFitResult<T> mmap3k_fit(const Matrix<T>& D0, const Matrix<T>& D1, const std::vector<T>& P,
                             const std::vector<T>& F, const std::vector<T>& B,
                             const std::vector<T>& B2 = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental, "mmap3k_fit inverts a moment system");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = P.size(), n = D0.rows();
    if (K == 0) throw InputError("mmap3k_fit: no classes given");
    if (F.size() != K || B.size() != K)
        throw InputError("mmap3k_fit: P, F and B must have the same length");
    if (n == 0 || D1.rows() != n) throw InputError("mmap3k_fit: D0 and D1 disagree");
    const double feastol = 1e-8;

    // The entries of D1 that can carry a mark.
    std::vector<std::pair<std::size_t, std::size_t>> nz;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (!(D1(i, j) == zero)) nz.push_back(std::make_pair(i, j));
    const std::size_t z = nz.size();
    if (z == 0) throw InputError("mmap3k_fit: the MAP has no arrival transitions to mark");

    const std::vector<std::pair<unsigned, unsigned>> orders = kfitdetail::marking_orders(n, z);
    std::vector<T> b2 = B2;
    if (b2.empty()) {
        for (std::size_t i = 0; i < orders.size(); ++i)
            if (orders[i].first == 3)
                throw InputError(
                    "mmap3k_fit: this MAP's marking system includes the third-order backward "
                    "characteristic, so the second-order backward moments B2 are required");
        b2.assign(K, zero);
    }

    Matrix<T> negD0(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negD0(i, j) = -D0(i, j);
    const Matrix<T> A = inverse(negD0);

    // The stationary law at arrival epochs: pie (A D1) = pie, sum pie = 1.
    const Matrix<T> Pemb = matmul(A, D1);
    Matrix<T> Tm(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Tm(i, j) = T(Pemb(j, i) - (i == j ? one : zero));
    for (std::size_t j = 0; j < n; ++j) Tm(n - 1, j) = one;
    std::vector<T> rhs(n, zero);
    rhs[n - 1] = one;
    const std::vector<T> pie = solve(Tm, rhs);

    // M(ii, jj) = pie A^a Dc_jj A^b e for the characteristic (a, b) of row ii.
    std::vector<Matrix<T>> Apow(1, eye<T>(n));
    std::size_t maxpow = 0;
    for (std::size_t i = 0; i < orders.size(); ++i)
        maxpow = std::max<std::size_t>(maxpow, std::max(orders[i].first, orders[i].second));
    for (std::size_t k = 1; k <= maxpow; ++k) Apow.push_back(matmul(Apow[k - 1], A));

    Matrix<T> M(z, z, zero);
    for (std::size_t jj = 0; jj < z; ++jj) {
        Matrix<T> Dc(n, n, zero);
        Dc(nz[jj].first, nz[jj].second) = D1(nz[jj].first, nz[jj].second);
        for (std::size_t ii = 0; ii < orders.size(); ++ii) {
            const Matrix<T> L = matmul(Apow[orders[ii].first], Dc);
            const Matrix<T> R = matmul(L, Apow[orders[ii].second]);
            T acc = zero;
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b) acc += pie[a] * R(a, b);
            M(ii, jj) = acc;
        }
    }

    // The degenerate locus shows up as a singular marking matrix; `solve`
    // reports it, and it is rethrown by name rather than left as a bare
    // linear-algebra failure the caller cannot act on.
    std::vector<std::vector<T>> q(z, std::vector<T>(K, zero));
    for (std::size_t c = 0; c < K; ++c) {
        std::vector<T> y(z, zero);
        for (std::size_t ii = 0; ii < orders.size(); ++ii) {
            const unsigned a = orders[ii].first, b = orders[ii].second;
            if (a == 1 && b == 0) y[ii] = P[c];
            else if (a == 1 && b == 1) y[ii] = T(P[c] * F[c]);
            else if (a == 2 && b == 0) y[ii] = T(P[c] * B[c]);
            else if (a == 3 && b == 0) y[ii] = T(P[c] * b2[c]);
            else
                throw InputError(
                    "mmap3k_fit: no target is supplied for one of the marking characteristics");
        }
        std::vector<T> x;
        try {
            x = solve(M, y);
        } catch (const NumericError&) {
            throw NumericError(
                "mmap3k_fit: the underlying MAP is on the degenerate locus of the marking system, "
                "where the characteristics do not determine the split");
        }
        for (std::size_t jj = 0; jj < z; ++jj) q[jj][c] = x[jj];
    }

    double viol = 0.0;
    for (std::size_t jj = 0; jj < z; ++jj) {
        T s = zero;
        for (std::size_t c = 0; c < K; ++c) {
            const double v = num_traits<T>::to_double(q[jj][c]);
            viol = std::max(viol, std::max(0.0, -v));
            viol = std::max(viol, std::max(0.0, v - 1.0));
            s += q[jj][c];
        }
        viol = std::max(viol, std::fabs(num_traits<T>::to_double(s) - 1.0));
    }

    MmapKFitResult<T> out;
    out.exact = viol <= feastol;
    out.mmap.D0 = D0;
    out.mmap.D1 = D1;
    out.mmap.Dc.assign(K, Matrix<T>(n, n, zero));
    for (std::size_t c = 0; c < K; ++c)
        for (std::size_t jj = 0; jj < z; ++jj) {
            T v = q[jj][c];
            if (v < zero) v = zero;
            if (v > one) v = one;
            out.mmap.Dc[c](nz[jj].first, nz[jj].second) =
                T(D1(nz[jj].first, nz[jj].second) * v);
        }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_K_FIT_H
