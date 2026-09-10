/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_COMPRESS_H
#define LINE_API_MAM_MMAP_COMPRESS_H

/**
 * Compression of a marked MAP into a smaller representation, and the two M3A
 * primitives it is built from: the class-conditional backward moments and the
 * probabilistic mixture of MAPs.
 *
 * Templated port of matlab/src/api/mam/mmap_compress.m (the 'default' /
 * 'mixture' / 'mixture.order1' method), matlab/lib/m3a/m3a/mmap/
 * mmap_backward_moment.m and mmap_mixture.m.
 *
 * ORDER-1 MIXTURE. One component per class, recombined by mmap_mixture with
 * the class probabilities p_c as mixing weights. Component c must carry the
 * law of the inter-arrival time CONDITIONED ON THE ARRIVAL THAT ENDS IT BEING
 * OF CLASS c, because mmap_mixture marks the arrival LEAVING component c with
 * class c, so the class of an arrival and the interval preceding it are both
 * governed by the component active during that interval. That conditional law
 * is the class-c BACKWARD moment set B(c, 1:3),
 *
 *     B(c,k) = k! pie (-D0)^-(k+1) D1^(c) e / p_c,
 *
 * i.e. E[T^k | class of the ending arrival = c]. It is NOT the forward moment
 * and it is NOT the class-c marginal MAP of mmap_maps, whose mean is
 * 1/lambda_c, the time between successive class-c arrivals; mixing those with
 * weights lambda_c/Lambda inflates the mean to K/Lambda.
 *
 * PRESERVED exactly: the aggregate moments 1..3, through the mixture law
 * M_k = sum_c B(k,c) p_c (M1 always, since aph2_adjust never alters M1; M2 and
 * M3 whenever the triple is APH(2)-feasible); the class probabilities p_c and
 * hence the per-class rates lambda_c = p_c / M1; the marking consistency
 * D1 = sum_c D1^(c); and MAP feasibility, since mmap_normalize closes the
 * result.
 *
 * LOST by construction: every autocorrelation. mmap_mixture re-enters each
 * component at its map_pie on every arrival, so the intervals are i.i.d. and
 * the result is a RENEWAL process: acf -> 0, IDC -> the SCV-determined renewal
 * value, and the class sequence becomes i.i.d. Retaining the class-transition
 * matrix sigma is what the order-2 method buys with its K^2 components.
 *
 * ARITHMETIC. mmap_backward_moment and mmap_mixture are finite rational
 * expressions in the descriptor entries and instantiate at Rational, which is
 * what makes the preservation claims above CHECKABLE rather than merely
 * plausible: at exact arithmetic the class probabilities of the compressed
 * MMAP equal those of the original digit for digit, so any drift the tests see
 * is a real modelling loss and not rounding. mmap_compress itself is gated on
 * num_traits<T>::has_transcendental because aph2_fit is.
 *
 * NOT PORTED. The other five methods of mmap_compress.m are dispatch to
 * separate fitting families that are outside this change:
 * 'mixture.order2' (mmap_mixture_fit_mmap), 'mamap2' and 'mamap2.fb'
 * (mamap2m_fit_mmap, mamap2m_fit_gamma_fb_mmap) and the four 'm3pp.*'
 * variants (m3pp2m_fitc_theoretical, which needs the derivest numerical
 * differentiation package). MmapCompressMethod names them so a caller gets
 * UnsupportedError identifying the missing family rather than silently
 * receiving the order-1 mixture.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Class-conditional backward moments of an MMAP (mmap_backward_moment.m).
 *
 * @param orders the moment orders to compute
 * @param normalized true for B(c,k) with M_k = sum_c B(c,k) p_c, i.e. divided
 *        by the class probability p_c (the MATLAB default); false for the
 *        unnormalized form with M_k = sum_c B(c,k)
 * @param m the marked MAP whose backward moments are taken
 * @return B[c][h], the moment of order orders[h] for class c
 */
template <class T>
std::vector<std::vector<T>> mmap_backward_moment(const Mmap<T>& m,
                                                 const std::vector<unsigned>& orders,
                                                 bool normalized) {
    const std::size_t n = m.order();
    const std::size_t C = m.classes();
    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> pie = map_pie(m.map());
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> Minv = inverse(negD0);

    std::vector<std::vector<T>> B(C, std::vector<T>(orders.size(), zero));
    for (std::size_t c = 0; c < C; ++c) {
        T pa = num_traits<T>::from_int(1);
        if (normalized) {
            const std::vector<T> t = vecmul(vecmul(pie, Minv), m.Dc[c]);
            pa = zero;
            for (const T& v : t) pa += v;
            if (pa == zero)
                throw NumericError(
                    "mmap_backward_moment: class with zero arrival probability cannot be "
                    "normalized");
        }
        for (std::size_t h = 0; h < orders.size(); ++h) {
            const unsigned k = orders[h];
            const std::vector<T> t = vecmul(vecmul(pie, matpow(Minv, k + 1)), m.Dc[c]);
            T s = zero;
            for (const T& v : t) s += v;
            B[c][h] = num_factorial<T>(k) / pa * s;
        }
    }
    return B;
}

/** mmap_backward_moment with the MATLAB default, normalized. */
template <class T>
std::vector<std::vector<T>> mmap_backward_moment(const Mmap<T>& m,
                                                 const std::vector<unsigned>& orders) {
    return mmap_backward_moment(m, orders, true);
}

/**
 * Probabilistic mixture of MAPs (mmap_mixture.m).
 *
 * The phase space is the disjoint union of the component phase spaces, the
 * hidden generator is block diagonal, and on completion of an interval in
 * component i the process jumps into component j with probability alpha_j,
 * entering it at its own map_pie. The arrival that LEAVES component i is
 * marked with class i, so the resulting MMAP has one class per component:
 *
 *     D0 = blkdiag(D0^1, ..., D0^I)
 *     D1 block (i,j) = alpha_j (D1^i e) pie^j
 *     D1^(c) block (i,j) = D1 block (i,j) if i == c, else 0.
 *
 * The result is a renewal process by construction; see the header note.
 */
template <class T>
Mmap<T> mmap_mixture(const std::vector<T>& alpha, const std::vector<Map<T>>& maps) {
    const std::size_t I = maps.size();
    if (I == 0) throw InputError("mmap_mixture: no components");
    if (alpha.size() != I) throw InputError("mmap_mixture: one weight per component is required");
    const T zero = num_traits<T>::from_int(0);

    std::vector<std::size_t> sz(I), off(I);
    std::size_t total = 0;
    for (std::size_t i = 0; i < I; ++i) {
        sz[i] = maps[i].order();
        off[i] = total;
        total += sz[i];
    }

    std::vector<std::vector<T>> pies(I);
    for (std::size_t j = 0; j < I; ++j) pies[j] = map_pie(maps[j]);

    Mmap<T> out;
    out.D0 = Matrix<T>(total, total, zero);
    out.D1 = Matrix<T>(total, total, zero);
    out.Dc.assign(I, Matrix<T>(total, total, zero));

    for (std::size_t i = 0; i < I; ++i) {
        for (std::size_t a = 0; a < sz[i]; ++a)
            for (std::size_t b = 0; b < sz[i]; ++b) out.D0(off[i] + a, off[i] + b) = maps[i].D0(a, b);
        // The completion rate out of each phase of component i.
        std::vector<T> t(sz[i], zero);
        for (std::size_t a = 0; a < sz[i]; ++a)
            for (std::size_t b = 0; b < sz[i]; ++b) t[a] += maps[i].D1(a, b);
        for (std::size_t j = 0; j < I; ++j)
            for (std::size_t a = 0; a < sz[i]; ++a)
                for (std::size_t b = 0; b < sz[j]; ++b) {
                    const T v = alpha[j] * t[a] * pies[j][b];
                    out.D1(off[i] + a, off[j] + b) = v;
                    out.Dc[i](off[i] + a, off[j] + b) = v;
                }
    }
    return mmap_normalize(out);
}

/** The compression methods of mmap_compress.m. */
enum class MmapCompressMethod {
    MixtureOrder1,  ///< 'default', 'mixture', 'mixture.order1'
    MixtureOrder2,  ///< 'mixture.order2', not ported
    Mamap2,         ///< 'mamap2', not ported
    Mamap2Fb,       ///< 'mamap2.fb', not ported
    M3ppApproxCov,  ///< 'm3pp.approx_cov', not ported
    M3ppApproxAg,   ///< 'm3pp.approx_ag', not ported
    M3ppExactDelta, ///< 'm3pp.exact_delta', not ported
    M3ppApproxDelta ///< 'm3pp.approx_delta', not ported
};

/**
 * Compress an MMAP (mmap_compress.m).
 *
 * A class that never arrives (p_c <= 1e-14, the reference's
 * GlobalConstants.Zero) gets Exp(1) as its component, exactly as in the
 * reference: it carries zero mixture weight, so any proper MAP leaves the
 * result unchanged and the 0/0 normalization of its backward moments is
 * avoided.
 */
template <class T>
Mmap<T> mmap_compress(const Mmap<T>& in, MmapCompressMethod method) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmap_compress requires transcendental arithmetic");
    if (method != MmapCompressMethod::MixtureOrder1)
        throw UnsupportedError(
            "mmap_compress: only the order-1 mixture is ported; the 'mixture.order2', 'mamap2', "
            "'mamap2.fb' and 'm3pp.*' methods dispatch to fitting families that are not part of "
            "this port (mmap_mixture_fit_mmap, mamap2m_fit_mmap, mamap2m_fit_gamma_fb_mmap, "
            "m3pp2m_fitc_theoretical)");

    const std::size_t K = in.classes();
    if (K == 0) throw InputError("mmap_compress: the MMAP has no classes");
    const std::vector<T> p = mmap_pc(in);
    std::vector<unsigned> orders;
    orders.push_back(1u);
    orders.push_back(2u);
    orders.push_back(3u);

    const T zeroTol = num_traits<T>::from_double(1e-14);
    // never-arriving-class 0/0 row: see _kb/03-api-layer.md (cpp port notes: mam)
    std::vector<Map<T>> comps;
    comps.reserve(K);
    for (std::size_t k = 0; k < K; ++k) {
        if (p[k] <= zeroTol) {
            comps.push_back(map_exponential(T(num_traits<T>::from_int(1))));
            continue;
        }
        Mmap<T> one = in;
        one.Dc.assign(1, in.Dc[k]);  // same D0 and D1, hence the same pie
        const std::vector<std::vector<T>> B = mmap_backward_moment(one, orders, true);
        comps.push_back(aph2_fit(B[0][0], B[0][1], B[0][2]).aph);
    }
    return mmap_normalize(mmap_mixture(p, comps));
}

/** mmap_compress with the default method, the order-1 mixture. */
template <class T>
Mmap<T> mmap_compress(const Mmap<T>& in) {
    return mmap_compress(in, MmapCompressMethod::MixtureOrder1);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_COMPRESS_H
