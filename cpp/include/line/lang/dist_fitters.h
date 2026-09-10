/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_DIST_FITTERS_H
#define LINE_LANG_DIST_FITTERS_H

/**
 * The moment fitters the reference distributions carry as STATIC FACTORIES:
 * `Erlang.fitMeanAndOrder`, `HyperExp.fitMeanAndSCV`, `Coxian.fitMeanAndSCV`,
 * `Cox2.fitCentral`, `APH.fitMeanAndSCV`, `APH.fitCentral`,
 * `Gamma.fitMeanAndSCV`, `Pareto.fitMeanAndSCV`.
 *
 * They live BESIDE `Distrib` rather than inside it because they are the only
 * part of the distribution layer that needs `line::mam`: fitting an APH is
 * `aph_fit`, fitting a two-phase hyperexponential is `map_hyperexp`, and
 * `lang/lang_types.h` sits UNDER `api/` in the include order. A model script
 * that only constructs distributions by parameter never pays for this header.
 *
 * WHAT IS A PORT AND WHAT IS A SUBSTITUTION. Every branch below is the
 * reference's arithmetic transcribed, with one documented exception:
 * `APH.fitCentral`/`fitMeanAndSCV` call BUTools' `APHFrom3Moments`, which this
 * port does not transcribe, so they go through
 * `mam::aph_fit` -- the same Bobbio-Horvath-Telek canonical APH, and the same
 * substitution `api/qsys/qsys_mapg1.h:37-45` already documents.
 *
 * EVERY FITTER IS GATED ON TRANSCENDENTAL ARITHMETIC. Each one takes a square
 * root of a moment discriminant, which has no exact rational counterpart; the
 * gate is `mam::map_hyperexp`'s and is stated here rather than discovered as a
 * link error in the exact instantiation.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/aph_fit.h"
#include "line/api/mam/aph_fit_moments.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace lang {

namespace fitdetail {

/** sqrt and pow found by ADL, so a Boost.Multiprecision T resolves its own. */
template <class T>
T fsqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

template <class T>
T fpow(const T& b, const T& e) {
    using std::pow;
    return pow(b, e);
}

}  // namespace fitdetail

/**
 * A PH distribution from a fitted (D0, D1) pair.
 *
 * The initial vector is recovered from D1 rather than carried alongside it:
 * D1 = (-D0 e) alpha by construction, so row i of D1 is alpha scaled by phase
 * i's exit rate and ANY row with a positive exit rate recovers it.
 *
 * IT CANNOT BE ROW 0. The canonical APH `aph_fit` returns is BIDIAGONAL: phase
 * 1 moves to phase 2 and never completes, so its exit rate is exactly zero and
 * its D1 row is all zeros. Reading alpha off row 0 threw on every APH fit.
 */
template <class T>
Distrib<T> ph_from_map(const mam::Map<T>& m, bool acyclic) {
    const std::size_t n = m.D0.rows();
    for (std::size_t i = 0; i < n; ++i) {
        T out = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) out += m.D1(i, j);
        if (!(out > num_traits<T>::from_int(0))) continue;
        std::vector<T> alpha(n);
        for (std::size_t j = 0; j < n; ++j) alpha[j] = T(m.D1(i, j) / out);
        return Distrib<T>::phase_type(alpha, m.D0, acyclic);
    }
    throw NumericError("ph_from_map: no phase completes, so the fitted PH has no alpha");
}

// ---------------------------------------------------------------------------
// Erlang
// ---------------------------------------------------------------------------

/** `Erlang.fitMeanAndOrder(MEAN, k)`: k phases, each of rate k / MEAN. */
template <class T>
Distrib<T> erlang_fit_mean_order(const T& mean, std::size_t k) {
    if (k == 0) throw InputError("Erlang.fitMeanAndOrder: the order must be positive");
    if (!(mean > num_traits<T>::from_int(0)))
        throw InputError("Erlang.fitMeanAndOrder: the mean must be positive");
    return Distrib<T>::erlang(T(num_traits<T>::from_int(static_cast<int>(k)) / mean), k);
}

// ---------------------------------------------------------------------------
// HyperExp
// ---------------------------------------------------------------------------

/**
 * `HyperExp.fitMeanAndSCV(MEAN, SCV)`, which is `map_hyperexp` at p = 0.99
 * read back as (p, mu1, mu2).
 */
template <class T>
Distrib<T> hyperexp_fit_mean_scv(const T& mean, const T& scv) {
    static_assert(num_traits<T>::has_transcendental,
                  "HyperExp.fitMeanAndSCV needs a square root of the moment discriminant");
    const mam::Map<T> m = mam::map_hyperexp(mean, scv, num_traits<T>::from_double(0.99));
    // map_hyperexp returns D0 = diag(-mu1, -mu2) and D1(i, j) = mu_i * p_j, so
    // the branch probability is read off row 0 and the rates off the diagonal.
    const T mu1 = T(-m.D0(0, 0)), mu2 = T(-m.D0(1, 1));
    const T p = T(m.D1(0, 0) / mu1);
    return Distrib<T>::hyperexp(p, mu1, mu2);
}

/**
 * `HyperExp.fitMeanAndSCVBalanced(MEAN, SCV)`: the balanced-means branch,
 * p / mu1 = (1 - p) / mu2. Both roots are tried in the reference's order.
 */
template <class T>
Distrib<T> hyperexp_fit_mean_scv_balanced(const T& mean, const T& scv) {
    static_assert(num_traits<T>::has_transcendental,
                  "HyperExp.fitMeanAndSCVBalanced needs a square root");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T root = fitdetail::fsqrt<T>(T((scv - one) / (scv + one)));
    T p = T(one / two - root / two);
    T mu1 = T(-(two * (root / two - one / two)) / mean);
    if (!(mu1 > num_traits<T>::from_int(0)) || !(p > num_traits<T>::from_int(0)) || !(p < one)) {
        p = T(root / two + one / two);
        mu1 = T((two * (root / two + one / two)) / mean);
    }
    const T mu2 = T((one - p) / p * mu1);
    return Distrib<T>::hyperexp(p, mu1, mu2);
}

// ---------------------------------------------------------------------------
// Coxian
// ---------------------------------------------------------------------------

/**
 * `Coxian.fitMeanAndSCV(MEAN, SCV)`, branch for branch.
 *
 * SCV below 1/2 is matched by an ERLANG of order ceil(1/SCV) expressed as a
 * Coxian with every phi zero but the last -- the reference's own choice, and
 * the reason this cannot be routed through `Distrib::erlang`: the order it
 * picks matches the mean exactly and the SCV only from below.
 */
template <class T>
Distrib<T> coxian_fit_mean_scv(const T& mean, const T& scv) {
    static_assert(num_traits<T>::has_transcendental, "Coxian.fitMeanAndSCV needs a square root");
    const double c2 = num_traits<T>::to_double(scv);
    const double tol = GlobalConstants::CoarseTol;
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    std::vector<T> mu, phi;
    if (c2 >= 1.0 - tol && c2 <= 1.0 + tol) {
        mu.push_back(T(one / mean));
        phi.push_back(one);
    } else if (c2 > 0.5 + tol && c2 < 1.0 - tol) {
        const T r = fitdetail::fsqrt<T>(T(one + two * (scv - one)));
        mu.push_back(T(two / mean / (one + r)));
        mu.push_back(T(two / mean / (one - r)));
        phi.push_back(num_traits<T>::from_int(0));
        phi.push_back(one);
    } else if (c2 <= 0.5 + tol) {
        const std::size_t n = static_cast<std::size_t>(std::ceil(1.0 / c2));
        const T lambda = T(num_traits<T>::from_int(static_cast<int>(n)) / mean);
        mu.assign(n, lambda);
        phi.assign(n, num_traits<T>::from_int(0));
        phi[n - 1] = one;
    } else {
        // SCV > 1: the two-phase hyperexponential written as a Coxian.
        mu.push_back(T(two / mean));
        mu.push_back(T(mu[0] / (two * scv)));
        phi.push_back(T(one - mu[1] / mu[0]));
        phi.push_back(one);
    }
    return Distrib<T>::coxian(mu, phi);
}

/**
 * `Cox2.fitCentral(MEAN, VAR, SKEW)`: the two-phase Coxian matching three
 * central moments exactly when the moment set admits one.
 *
 * Both roots of the moment condition are tried in the reference's order, and
 * the fallback when neither is feasible is the reference's: `fitMeanAndSCV`
 * above SCV = 1/2 and the exponential of that mean below it.
 */
template <class T>
Distrib<T> cox2_fit_central(const T& mean, const T& var, const T& skew) {
    static_assert(num_traits<T>::has_transcendental, "Cox2.fitCentral needs a square root");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T scv = T(var / (mean * mean));
    const T e1 = mean;
    const T e2 = T((one + scv) * e1 * e1);
    const T e3 = T(-(num_traits<T>::from_int(2) * e1 * e1 * e1 -
                     num_traits<T>::from_int(3) * e1 * e2 -
                     skew * fitdetail::fpow<T>(T(e2 - e1 * e1), num_traits<T>::from_double(1.5))));

    const T disc = T(num_traits<T>::from_int(24) * e1 * e1 * e1 * e3 -
                     num_traits<T>::from_int(27) * e1 * e1 * e2 * e2 -
                     num_traits<T>::from_int(18) * e1 * e2 * e3 +
                     num_traits<T>::from_int(18) * e2 * e2 * e2 + e3 * e3);
    const T den = T(num_traits<T>::from_int(-3) * e2 * e2 + num_traits<T>::from_int(2) * e1 * e3);
    if (num_traits<T>::to_double(disc) >= 0.0 && den != zero) {
        const T s = fitdetail::fsqrt<T>(disc);
        const T a = T(num_traits<T>::from_int(3) * e1 * e2);
        T mu1[2], mu2[2];
        mu1[0] = T(num_traits<T>::from_int(2) * (e3 - a) / den + (a - e3 + s) / den);
        mu1[1] = T(num_traits<T>::from_int(2) * (e3 - a) / den - (e3 - a + s) / den);
        mu2[0] = T(-(a - e3 + s) / den);
        mu2[1] = T((e3 - a + s) / den);
        for (int k = 0; k < 2; ++k) {
            const T phi = T(one - mu2[k] * e1 + mu2[k] / mu1[k]);
            if (num_traits<T>::to_double(phi) >= 0.0 && num_traits<T>::to_double(phi) <= 1.0 &&
                num_traits<T>::to_double(mu1[k]) >= 0.0 && num_traits<T>::to_double(mu2[k]) >= 0.0)
                return Distrib<T>::cox2(mu1[k], mu2[k], phi);
        }
    }
    if (num_traits<T>::to_double(scv) >= 0.5) return coxian_fit_mean_scv(mean, scv);
    return Distrib<T>::exp_mean(mean);
}

/** `Coxian.fitCentral`, which the reference forwards to `Cox2.fitCentral`. */
template <class T>
Distrib<T> coxian_fit_central(const T& mean, const T& var, const T& skew) {
    return cox2_fit_central(mean, var, skew);
}

// ---------------------------------------------------------------------------
// APH
// ---------------------------------------------------------------------------

/** `APH.fitMeanAndSCV(MEAN, SCV)`, through `mam::aph_fit_mean_scv`. */
template <class T>
Distrib<T> aph_fit_mean_scv(const T& mean, const T& scv) {
    return ph_from_map(mam::aph_fit_mean_scv(mean, scv), true);
}

/**
 * `APH.fitCentral(MEAN, VAR, SKEW)`: the three central moments converted to
 * raw ones and matched by a canonical APH.
 */
template <class T>
Distrib<T> aph_fit_central(const T& mean, const T& var, const T& skew) {
    static_assert(num_traits<T>::has_transcendental, "APH.fitCentral needs a square root");
    const T one = num_traits<T>::from_int(1);
    const T scv = T(var / (mean * mean));
    const T e1 = mean;
    const T e2 = T((one + scv) * e1 * e1);
    const T e3 = T(-(num_traits<T>::from_int(2) * e1 * e1 * e1 -
                     num_traits<T>::from_int(3) * e1 * e2 -
                     skew * fitdetail::fpow<T>(T(e2 - e1 * e1), num_traits<T>::from_double(1.5))));
    return ph_from_map(mam::aph_fit(e1, e2, e3).aph, true);
}

// ---------------------------------------------------------------------------
// Gamma and Pareto
// ---------------------------------------------------------------------------

/** `Gamma.fitMeanAndSCV(MEAN, SCV)`: shape 1/SCV, scale MEAN * SCV. */
template <class T>
Distrib<T> gamma_fit_mean_scv(const T& mean, const T& scv) {
    const T shape = T(num_traits<T>::from_int(1) / scv);
    return Distrib<T>::gamma_dist(shape, T(mean / shape));
}

/**
 * `Pareto.fitMeanAndSCV(MEAN, SCV)`: alpha = 1 + sqrt(1 + 1/SCV) and
 * k = MEAN (alpha - 1) / alpha.
 */
template <class T>
Distrib<T> pareto_fit_mean_scv(const T& mean, const T& scv) {
    static_assert(num_traits<T>::has_transcendental, "Pareto.fitMeanAndSCV needs a square root");
    const T one = num_traits<T>::from_int(1);
    const T shape = T(one + fitdetail::fsqrt<T>(T(one + one / scv)));
    return Distrib<T>::pareto(shape, T(mean * (shape - one) / shape));
}

}  // namespace lang
}  // namespace line

#endif  // LINE_LANG_DIST_FITTERS_H
