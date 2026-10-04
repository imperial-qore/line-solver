/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_COMPRESS_H
#define LINE_API_MAM_MMAP_COMPRESS_H

/**
 * @file
 * @ingroup api_mam
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
 * THE OTHER METHODS dispatch to their fitting families, all ported in their
 * own headers: 'mixture.order2' to mmap_mixture_fit_mmap (mmap_modulate.h),
 * 'mamap2' and 'mamap2.fb' to mamap2m_fit_mmap and mamap2m_fit_gamma_fb_mmap
 * (mamap2m_fit.h), and the four 'm3pp.*' variants to m3pp2m_fitc_theoretical
 * (m3pp2m_interleave.h) at the reference's scales t = 1, tinf = 1e6. Those
 * families carry their own acceptance contracts where the reference solves a
 * program (quadprog, YALMIP): read the header of each before comparing digits.
 * Every method ends in mmap_normalize, as mmap_compress.m does.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/m3pp2m_interleave.h"
#include "line/api/mam/mamap2m_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_backward_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_modulate.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

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
    MixtureOrder2,  ///< 'mixture.order2'
    Mamap2,         ///< 'mamap2'
    Mamap2Fb,       ///< 'mamap2.fb'
    M3ppApproxCov,  ///< 'm3pp.approx_cov'
    M3ppApproxAg,   ///< 'm3pp.approx_ag'
    M3ppExactDelta, ///< 'm3pp.exact_delta'
    M3ppApproxDelta ///< 'm3pp.approx_delta'
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
    const T t1 = num_traits<T>::from_int(1), tinf = num_traits<T>::from_double(1e6);
    switch (method) {
        case MmapCompressMethod::MixtureOrder1:
            break;
        case MmapCompressMethod::MixtureOrder2:
            return mmap_normalize(mmap_mixture_fit_mmap(in));
        case MmapCompressMethod::Mamap2:
            return mmap_normalize(mamap2m_fit_mmap(mmap_normalize(in)));
        case MmapCompressMethod::Mamap2Fb:
            return mmap_normalize(mamap2m_fit_gamma_fb_mmap(mmap_normalize(in)));
        case MmapCompressMethod::M3ppApproxCov:
            return mmap_normalize(m3pp2m_fitc_theoretical(in, std::string("approx_cov"), t1, tinf));
        case MmapCompressMethod::M3ppApproxAg:
            return mmap_normalize(m3pp2m_fitc_theoretical(in, std::string("approx_ag"), t1, tinf));
        case MmapCompressMethod::M3ppExactDelta:
            return mmap_normalize(m3pp2m_fitc_theoretical(in, std::string("exact_delta"), t1, tinf));
        case MmapCompressMethod::M3ppApproxDelta:
            return mmap_normalize(m3pp2m_fitc_theoretical(in, std::string("approx_delta"), t1, tinf));
    }

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
