/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_ASSEMBLE_H
#define LINE_API_MAM_MMAP_ASSEMBLE_H

/**
 * The MMAP assembly primitives `solver_mam_basic.m` builds its per-station
 * arrival stream from: `mmap_exponential`, the probabilistic `mmap_mark`, the
 * per-class `mmap_scale`, and `mmap_super_safe`.
 *
 * They are separate from `mmap_lambda.h` because each of them is a DIFFERENT
 * function from the same-named one already there: `mmap_lambda.h`'s
 * `mmap_mark` splits a MAP by per-PHASE weights and its `mmap_scale` takes a
 * single target mean, which are the M3A signatures; the MAM solver calls the
 * kpctoolbox ones, which split an MMAP by per-CLASS probabilities and target
 * one mean per class. Naming them apart is deliberate: two functions with one
 * name and two meanings is exactly the failure mode the parity notes record for
 * `map_normalize`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Order-n MMAP with the given per-class arrival rates (mmap_exponential.m).
 *
 * The per-class matrix is `flip(eye(n)) * lambda_c`, so at n = 1 this is the
 * ordinary marked Poisson stream and at n > 1 it is the n-phase cycle the
 * reference uses as a neutral element of the superposition.
 */
template <class T>
Mmap<T> mmap_exponential_vec(const std::vector<T>& lambda, std::size_t n = 1) {
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) throw InputError("mmap_exponential_vec: order must be positive");
    Mmap<T> m;
    m.D0 = Matrix<T>(n, n, zero);
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t c = 0; c < lambda.size(); ++c) {
        Matrix<T> Dc(n, n, zero);
        for (std::size_t i = 0; i < n; ++i) Dc(i, n - 1 - i) = lambda[c];
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) m.D1(i, j) += Dc(i, j);
        m.Dc.push_back(Dc);
    }
    return mmap_normalize(m);
}

/**
 * Re-mark an MMAP by a (K x R) probability matrix (mmap_mark.m): a type-k
 * arrival is reported as class r with probability prob(k,r).
 */
template <class T>
Mmap<T> mmap_mark_probs(const Mmap<T>& in, const Matrix<T>& prob) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = prob.rows(), R = prob.cols();
    if (K > in.classes())
        throw InputError("mmap_mark_probs: the probability matrix has more input types than the "
                         "MMAP has classes");
    Mmap<T> m;
    m.D0 = in.D0;
    m.D1 = in.D1;
    const std::size_t n = in.order();
    for (std::size_t r = 0; r < R; ++r) {
        Matrix<T> Dr(n, n, zero);
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) Dr(i, j) += in.Dc[k](i, j) * prob(k, r);
        m.Dc.push_back(Dr);
    }
    return m;
}

/**
 * Retarget the per-class MEAN inter-arrival times (mmap_scale.m, vector form).
 *
 * Each class matrix is rescaled by (1/M_c)/lambda_c, then the MMAP is
 * renormalized. The reference calls this "heuristic because it also affects the
 * other classes"; the refinement loop that follows it in MATLAB is dead code
 * behind an unconditional `return`, so the heuristic IS the function and is
 * what is reproduced here. A class with zero rate is zeroed rather than divided.
 */
template <class T>
Mmap<T> mmap_scale_perclass(const Mmap<T>& in, const std::vector<T>& M) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t C = in.classes();
    if (M.size() != C) throw InputError("mmap_scale_perclass: one target mean per class is needed");
    const std::vector<T> l = mmap_count_lambda(in);
    Mmap<T> s;
    s.D0 = in.D0;
    s.D1 = Matrix<T>(in.order(), in.order(), zero);
    for (std::size_t c = 0; c < C; ++c) {
        Matrix<T> Dc = in.Dc[c];
        const T f = (l[c] > zero) ? T(T(num_traits<T>::from_int(1) / M[c]) / l[c]) : zero;
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) {
                Dc(i, j) *= f;
                s.D1(i, j) += Dc(i, j);
            }
        s.Dc.push_back(Dc);
    }
    return mmap_normalize(s);
}

namespace mmap_super_detail {

/** 1-norm of D1, used to detect a component that carries no arrivals at all. */
template <class T>
double arrival_norm(const Mmap<T>& m) {
    double best = 0.0;
    for (std::size_t j = 0; j < m.D1.cols(); ++j) {
        double col = 0.0;
        for (std::size_t i = 0; i < m.D1.rows(); ++i)
            col += std::fabs(num_traits<T>::to_double(m.D1(i, j)));
        if (col > best) best = col;
    }
    return best;
}

/** The order-1 marked Poisson stream carrying this MMAP's per-class rates. */
template <class T>
Mmap<T> to_poisson(const Mmap<T>& m) {
    return mmap_exponential_vec(mmap_count_lambda(m), 1);
}

}  // namespace mmap_super_detail

/**
 * Order-bounded superposition of several MMAPs (mmap_super_safe.m).
 *
 * Components are superposed low-SCV first, and the product order is held at or
 * below `maxorder` by replacing a component with its marked Poisson equivalent.
 * The marks are permuted back into INPUT order afterwards, because `mmap_super`
 * concatenates them in fold order while every caller reads mark k as its own
 * k-th class; without that the SCV sort renames the classes.
 *
 * ONE REFERENCE BRANCH IS REFUSED BY NAME rather than substituted. When the
 * order budget still allows an order-2 component, MATLAB compresses with
 * `mamap2m_fit_gamma_fb_mmap`, an acyclic MAP(2) fit that `mmap_compress.h`
 * records as not ported. Substituting the Poisson fallback there would silently
 * discard the component's variability, so this refuses instead. With the
 * solver's default `space_max = 128` the branch needs an arrival stream of
 * order above 128 (or a product above it with room for a 2-phase factor) to be
 * reachable at all.
 */
template <class T>
Mmap<T> mmap_super_safe(const std::vector<Mmap<T>>& in, std::size_t maxorder) {
    if (maxorder == 0) throw InputError("mmap_super_safe: maxorder must be positive");
    std::vector<Mmap<T>> parts;
    for (const Mmap<T>& m : in) {
        if (m.order() == 0) continue;
        // A component with an all-zero D1 has zero rate and an absorbing phase
        // generator, so map_scv would fail on it; canonicalize to the equivalent
        // order-1 null, whose superposition is the identity.
        if (m.order() > 1 && mmap_super_detail::arrival_norm(m) < 1e-13)
            parts.push_back(mmap_exponential_vec(
                std::vector<T>(m.classes(), num_traits<T>::from_int(0)), 1));
        else
            parts.push_back(m);
    }
    if (parts.empty()) throw InputError("mmap_super_safe: no components to superpose");

    std::vector<std::size_t> order(parts.size());
    std::iota(order.begin(), order.end(), 0u);
    // MATLAB sorts by map_scv, and a ZERO-RATE component (the neutral element
    // the MAM analyzer superposes to reshape a marking) has an infinite mean, so
    // its SCV is NaN there and MATLAB's ascending sort puts NaN LAST. Calling
    // map_scv on it here would divide by a zero rate, so the case is answered
    // with +infinity, which sorts last for the same reason.
    std::vector<double> scv(parts.size());
    for (std::size_t i = 0; i < parts.size(); ++i)
        scv[i] = mmap_super_detail::arrival_norm(parts[i]) > 0.0
                     ? num_traits<T>::to_double(map_scv(parts[i].map()))
                     : std::numeric_limits<double>::infinity();
    std::stable_sort(order.begin(), order.end(),
                     [&scv](std::size_t a, std::size_t b) { return scv[a] < scv[b]; });

    // Mark provenance: a zero-rate component sorts last (SCV +inf above), so a
    // chain that never visits the station used to push its marks ahead of one
    // that does, renaming both chains' classes.
    std::vector<std::size_t> markbase(parts.size() + 1, 0);
    for (std::size_t i = 0; i < parts.size(); ++i)
        markbase[i + 1] = markbase[i] + parts[i].classes();
    std::vector<std::size_t> outorder;

    bool first = true;
    Mmap<T> sup;
    for (std::size_t idx : order) {
        for (std::size_t j = 0; j < parts[idx].classes(); ++j)
            outorder.push_back(markbase[idx] + j);
        Mmap<T> cur = parts[idx];
        if (cur.order() > maxorder) {
            if (maxorder >= 2)
                throw UnsupportedError(
                    "mmap_super_safe: a component of order " + std::to_string(cur.order()) +
                    " exceeds the order budget and the reference compresses it with "
                    "mamap2m_fit_gamma_fb_mmap, which is not ported to C++");
            cur = mmap_super_detail::to_poisson(cur);
        }
        if (first) {
            sup = (maxorder == 1) ? mmap_super_detail::to_poisson(cur) : cur;
            first = false;
            continue;
        }
        if (sup.order() * cur.order() > maxorder) {
            if (sup.order() * 2 <= maxorder)
                throw UnsupportedError(
                    "mmap_super_safe: the superposition exceeds the order budget and the "
                    "reference compresses the next component with mamap2m_fit_gamma_fb_mmap, "
                    "which is not ported to C++");
            sup = mmap_super(sup, mmap_super_detail::to_poisson(cur));
        } else {
            sup = mmap_super(sup, cur);
        }
    }
    // Restore the caller's mark order.
    if (sup.Dc.size() == outorder.size() &&
        !std::is_sorted(outorder.begin(), outorder.end())) {
        std::vector<Matrix<T>> reordered(outorder.size());
        for (std::size_t j = 0; j < outorder.size(); ++j) reordered[outorder[j]] = sup.Dc[j];
        sup.Dc.swap(reordered);
    }
    return sup;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_ASSEMBLE_H
