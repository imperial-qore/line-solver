/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_TYPES_H
#define LINE_API_FJ_TYPES_H

/**
 * Shared return types and arithmetic helpers for the templated fork-join port.
 *
 * The MATLAB fork-join family in matlab/src/api/fj/ returns a scalar in most
 * cases but a pair in a few (fj_bounds, fj_char_max, fj_xmax_normal,
 * fj_xmax_pareto, fj_xmax_approx, fj_order_stat, fj_quorum_moments) and a
 * struct in one (fj_gk_bound). Those aggregates live here, together with the
 * handful of ADL wrappers and the binomial coefficient the family needs; each
 * ported function lives in its own header named after the MATLAB file, as
 * required by the port convention.
 *
 * The binomial coefficient is built multiplicatively so that every partial
 * product is an integer: in exact arithmetic C(n,k) is therefore exact for any
 * n, with no overflow and no rounding, which is what makes the alternating
 * sums in fj_respt_vm and fj_xmax_erlang trustworthy. Those sums cancel
 * catastrophically in double past K ~ 20; the exact instantiation is the only
 * way to see how much of the double answer is left.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** Distribution families for which a G(K) standardized-maximum factor exists. */
enum class FJDistType { Exp, Uniform, Evd, Bound };

/** Bracketing methods for the normal-maximum approximation. */
enum class FJNormalMethod { Johnson, Arnold, Corrected };

/** [Rmax, Rmin] of fj_bounds: pessimistic and optimistic response-time bounds. */
template <class T>
struct FJBoundsResult {
    T Rmax;
    T Rmin;
};

/** [MK, mK] of fj_char_max: characteristic maximum and its threshold. */
template <class T>
struct FJCharMaxResult {
    T MK;
    T mK;
};

/** All four G(K) factors of fj_gk_bound in its 'all' mode. */
template <class T>
struct FJGKBoundResult {
    unsigned K;
    T exponential;
    T uniform;
    T evd;
    T upper_bound;
};

/** [Xmax, GK] of fj_xmax_approx. */
template <class T>
struct FJXmaxApproxResult {
    T Xmax;
    T GK;
};

/** [Xmax, Vmax] of fj_xmax_normal. */
template <class T>
struct FJXmaxNormalResult {
    T Xmax;
    T Vmax;
};

/** [Xmax, MK] of fj_xmax_pareto. */
template <class T>
struct FJXmaxParetoResult {
    T Xmax;
    T MK;
};

/**
 * [F_Yk, E_Yk] of fj_order_stat.
 *
 * The CDF is a polynomial in the base CDF value and so is available in every
 * arithmetic mode; the expected value is a quadrature and is only produced
 * when T carries transcendental functions. mean_available says which.
 */
template <class T>
struct FJOrderStatResult {
    T F_Yk;
    T E_Yk;
    bool mean_available;
};

/** [m, v] of fj_quorum_moments: mean and variance of the k-of-n join time. */
template <class T>
struct FJQuorumMomentsResult {
    T m;
    T v;
};

namespace detail {

/** exp(v), resolved by ADL so double, cpp_bin_float and mpfr all work. */
template <class T>
inline T num_exp(const T& v) {
    using std::exp;
    return exp(v);
}

/** log(v), resolved by ADL. */
template <class T>
inline T num_log(const T& v) {
    using std::log;
    return log(v);
}

/** sqrt(v), resolved by ADL. */
template <class T>
inline T num_sqrt(const T& v) {
    using std::sqrt;
    return sqrt(v);
}

/** base^exponent for a real-valued exponent, resolved by ADL. */
template <class T>
inline T num_pow(const T& base, const T& exponent) {
    using std::pow;
    return pow(base, exponent);
}

/** pi to the working precision of T, from the double literal MATLAB uses. */
template <class T>
inline T num_pi() {
    return num_traits<T>::from_double(3.14159265358979323846264338327950288);
}

/**
 * Binomial coefficient C(n,k) as a value of T, built multiplicatively.
 *
 * Every partial product r * (n-k+i) is divisible by i, so the running value is
 * an integer throughout and the result is exact in any field.
 */
template <class T>
inline T fj_binom(unsigned n, unsigned k) {
    if (k > n) return num_traits<T>::from_int(0);
    const unsigned kk = (k > n - k) ? n - k : k;
    T r = num_traits<T>::from_int(1);
    for (unsigned i = 1; i <= kk; ++i) {
        r *= num_traits<T>::from_int(static_cast<long>(n - kk + i));
        r /= num_traits<T>::from_int(static_cast<long>(i));
    }
    return r;
}

/**
 * Composite Simpson rule on [a,b] with n_intervals panels (n_intervals even).
 *
 * MATLAB uses the adaptive quad in `integral`; the JAR replaces it with a fixed
 * 10001-point composite Simpson everywhere it ports one of these functions
 * (FJ_xmax.fj_xmax_erlang, fj_xmax_pareto, FJ_rmax.fj_rmax_erlang). This port
 * follows the JAR convention so that the three non-MATLAB implementations agree
 * with each other, and keeps the same default panel count.
 */
template <class T, class F>
inline T simpson(const F& f, const T& a, const T& b, unsigned n_intervals = 10000) {
    if (n_intervals % 2 != 0 || n_intervals == 0)
        throw InputError("simpson: the panel count must be a positive even number");
    const T h = (b - a) / num_traits<T>::from_int(static_cast<long>(n_intervals));
    T acc = f(a) + f(b);
    const T two = num_traits<T>::from_int(2), four = num_traits<T>::from_int(4);
    for (unsigned i = 1; i < n_intervals; ++i) {
        const T x = a + h * num_traits<T>::from_int(static_cast<long>(i));
        acc += (i % 2 == 1 ? four : two) * f(x);
    }
    return acc * h / num_traits<T>::from_int(3);
}

/**
 * Bisection on a bracketed sign change, the replacement for MATLAB's fzero.
 *
 * MATLAB's fzero is Brent's method; the JAR replaces it with bisection
 * wherever it ports one of these functions (FJ_char_max, Aoi_fcfs_dm1), and
 * this port follows suit. With max_iter = 200 the bracket is narrowed below
 * any representable tolerance for double and for the 50-digit real type
 * alike, so the two agree with MATLAB to full precision.
 */
template <class T, class F>
inline T bisect(const F& f, T lo, T hi, const char* fn, unsigned max_iter = 200) {
    const T zero = num_traits<T>::from_int(0);
    T flo = f(lo), fhi = f(hi);
    if ((flo > zero && fhi > zero) || (flo < zero && fhi < zero))
        throw NumericError(std::string(fn) + ": the root is not bracketed by the initial interval");
    for (unsigned it = 0; it < max_iter; ++it) {
        const T mid = (lo + hi) / num_traits<T>::from_int(2);
        if (mid == lo || mid == hi) break;
        const T fm = f(mid);
        if (fm == zero) return mid;
        if ((fm > zero) == (flo > zero)) { lo = mid; flo = fm; }
        else { hi = mid; fhi = fm; }
    }
    return (lo + hi) / num_traits<T>::from_int(2);
}

/**
 * The standard normal quantile, MATLAB's sqrt(2)*erfinv(2u-1).
 *
 * C++ has no erfinv, so the inverse is taken by bisection on the complementary
 * error function, which std does provide. This is called once per query in the
 * fork-join family, so the ~60 iterations are free and the result is tight to
 * the last representable bit of double.
 */
inline double normal_quantile(double u) {
    if (!(u > 0.0) || !(u < 1.0))
        throw InputError("normal_quantile: the probability must lie in (0,1)");
    double lo = -40.0, hi = 40.0;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi) break;
        // Phi(mid) = erfc(-mid/sqrt(2))/2
        const double phi = 0.5 * std::erfc(-mid / std::sqrt(2.0));
        if (phi < u) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/** Guard shared by every function that forms a harmonic number. */
inline void require_positive_K(unsigned K, const char* fn) {
    if (K < 1) throw InputError(std::string(fn) + ": K must be a positive integer");
}

}  // namespace detail

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_TYPES_H
