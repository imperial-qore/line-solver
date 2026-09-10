/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_ASYMPT_COMMON_H
#define LINE_API_PFQN_ASYMPT_COMMON_H

/**
 * Shared scalar machinery for the integration / asymptotic members of the
 * pfqn family (pfqn_le, pfqn_lap, pfqn_ls, pfqn_cub, pfqn_kt, pfqn_panacea,
 * the McKenna-Mitra quadratures, pfqn_nrl / pfqn_nrp).
 *
 * Nothing here mirrors a MATLAB file of its own: MATLAB gets log-gamma from
 * gammaln, determinants from det, quadrature nodes from a .mat table and
 * complex arithmetic from the language. Those four facilities have to be
 * supplied explicitly in a templated port, and putting them in one header
 * keeps each ported algorithm a 1:1 image of its MATLAB source.
 *
 * ARITHMETIC. Every function here is inherently inexact -- log-gamma, the
 * Gauss node polynomials and the complex exponential all leave the rationals
 * -- so each is gated on num_traits<T>::has_transcendental. The determinant
 * is the exception: it is a finite sequence of field operations and is
 * therefore left exact and ungated, so a bound or a Hessian determinant can
 * still be evaluated in rational arithmetic.
 *
 * PRECISION CEILING of num_lgamma. For an integer argument the value is the
 * exact sum log 1 + ... + log n, accumulated in T, so a Real<D> instantiation
 * gains the full D digits. For a non-integer argument the Lanczos g=7
 * coefficients are double constants, which caps the relative accuracy near
 * 1e-15 whatever T is. Only pfqn_propfair evaluates log-gamma off the
 * integers, and its own optimizer tolerance is far above that ceiling.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {
namespace detail {

// ---------------------------------------------------------------------------
// Determinant (exact in a field, hence ungated)
// ---------------------------------------------------------------------------

/** Determinant by Gaussian elimination with partial pivoting; 0 if singular. */
template <class T>
T pfqn_det(Matrix<T> A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("pfqn_det: matrix is not square");
    const T zero = num_traits<T>::from_int(0);
    T d = num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        T amax = num_abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const T a = num_abs(A(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == zero) return zero;
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
            d = -d;
        }
        d *= A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            const T f = A(i, k) / A(k, k);
            for (std::size_t j = k; j < n; ++j) A(i, j) -= f * A(k, j);
        }
    }
    return d;
}

/**
 * Logarithm of |det A| accumulated over the pivots of the same elimination.
 *
 * det(A) of an R x R Hessian leaves double range well before its logarithm does
 * (it overflowed at R = 64 in pfqn_kt, turning lG into -inf), so callers that
 * only need log det must never form the determinant first.
 */
template <class T>
T pfqn_logdet(Matrix<T> A) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_logdet requires transcendental arithmetic (logarithm)");
    using std::log;
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("pfqn_logdet: matrix is not square");
    const T zero = num_traits<T>::from_int(0);
    T acc = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        T amax = num_abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const T a = num_abs(A(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == zero) throw InputError("pfqn_logdet: matrix is singular");
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
        }
        acc += log(num_abs(A(k, k)));
        for (std::size_t i = k + 1; i < n; ++i) {
            const T f = A(i, k) / A(k, k);
            for (std::size_t j = k; j < n; ++j) A(i, j) -= f * A(k, j);
        }
    }
    return acc;
}

// ---------------------------------------------------------------------------
// log-gamma / factln
// ---------------------------------------------------------------------------

/** log(n!) accumulated in T; exact up to the rounding of each log(k). */
template <class T>
T num_logfact_int(long n) {
    static_assert(num_traits<T>::has_transcendental,
                  "num_logfact_int requires transcendental arithmetic (logarithm)");
    if (n < 0) throw InputError("num_logfact_int: negative argument");
    using std::log;
    T s = num_traits<T>::from_int(0);
    for (long k = 2; k <= n; ++k) s += log(num_traits<T>::from_int(k));
    return s;
}

/** log Gamma(x) for x > 0, Lanczos g = 7 (see the precision note above). */
template <class T>
T num_lgamma(const T& x) {
    static_assert(num_traits<T>::has_transcendental,
                  "num_lgamma requires transcendental arithmetic");
    using std::log;
    using std::sin;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    if (x <= zero) throw InputError("num_lgamma: argument must be positive");
    // Integer arguments take the exact route, which is what every pfqn caller
    // but pfqn_propfair needs and which carries the full precision of T.
    const double xd = num_traits<T>::to_double(x);
    const double rn = std::floor(xd + 0.5);
    if (rn >= 1.0 && rn <= 1e6 && x == num_traits<T>::from_int(static_cast<long>(rn)))
        return num_logfact_int<T>(static_cast<long>(rn) - 1);

    static const double g[9] = {0.99999999999980993,   676.5203681218851,   -1259.1392167224028,
                                771.32342877765313,    -176.61502916214059, 12.507343278686905,
                                -0.13857109526572012,  9.9843695780195716e-6,
                                1.5056327351493116e-7};
    const T one = num_traits<T>::from_int(1);
    const T z = T(x - one);
    T a = num_traits<T>::from_double(g[0]);
    for (int i = 1; i < 9; ++i)
        a += num_traits<T>::from_double(g[i]) / T(z + num_traits<T>::from_int(i));
    const T t = T(z + num_traits<T>::from_double(7.5));
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    return T(num_traits<T>::from_rational(1, 2) * log(twopi) + T(z + num_traits<T>::from_rational(1, 2)) * log(t) -
             t + log(a));
}

/** MATLAB factln(n) = gammaln(1+n). */
template <class T>
T num_factln(const T& n) {
    return num_lgamma<T>(T(n + num_traits<T>::from_int(1)));
}

/** log(sum_i exp(v_i)), shifted by the maximum so no term overflows. */
template <class T>
T logsumexp(const std::vector<T>& v) {
    static_assert(num_traits<T>::has_transcendental,
                  "logsumexp requires transcendental arithmetic");
    using std::exp;
    using std::log;
    if (v.empty()) throw InputError("logsumexp: empty argument");
    T m = v[0];
    for (const T& x : v)
        if (x > m) m = x;
    // An all -inf input has no finite logarithm; return the maximum unchanged.
    if (!(num_traits<T>::to_double(m) > -std::numeric_limits<double>::infinity())) return m;
    T s = num_traits<T>::from_int(0);
    for (const T& x : v) s += exp(T(x - m));
    return T(m + log(s));
}

// ---------------------------------------------------------------------------
// Gauss quadrature nodes, generated rather than tabulated
// ---------------------------------------------------------------------------

/**
 * The first `count` nodes and weights of the n-point Gauss-Legendre rule on
 * [a,b], by Newton iteration on the Legendre polynomial with the standard
 * w = 2/((1-x^2) P'^2). `count` = 0 returns the whole rule.
 *
 * MATLAB's pfqn_mmint2_gausslegendre loads a table generated once in Julia by
 * the Golub-Welsch tridiagonal eigenvalue method -- and then uses only its
 * LEADING ENTRIES: the table is a 20000-point rule on [0, 1e6]
 * (matlab/src/api/pfqn/gausslegendre-nodes.txt, first node 0.00361, last
 * 999999.98), and the routine takes nodes 1..n with n = max(300, ...). Taking
 * a prefix of a Gauss rule is not itself a Gauss rule; what makes it work is
 * that the McKenna-Mitra integrand carries e^{-u}, so the first 300 nodes of
 * the 20000-point rule already span [0.0036, 557.8] and everything beyond is
 * below 1e-240. Regenerating a genuine 300-point rule on [0, 1e6] instead
 * would put the FIRST node at u = 13.7 and miss the mass entirely -- it
 * returns log G = -3.44 where the answer is 1.63 -- so the prefix is not an
 * implementation detail of the reference but part of its definition.
 *
 * The table cannot be carried across arithmetics without pinning every
 * instantiation to the precision it was generated at, so the rule is
 * regenerated in T; only the requested prefix is computed, since each Newton
 * iteration is independent of the other nodes.
 */
template <class T>
void gauss_legendre(std::size_t n, const T& a, const T& b, std::vector<T>& x, std::vector<T>& w,
                    std::size_t count = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "gauss_legendre requires transcendental arithmetic (cos, sqrt)");
    using std::cos;
    x.assign(n, num_traits<T>::from_int(0));
    w.assign(n, num_traits<T>::from_int(0));
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T half = num_traits<T>::from_rational(1, 2);
    const T pi = num_traits<T>::from_double(3.14159265358979323846264338328);
    const T xm = T(half * T(b + a));
    const T xl = T(half * T(b - a));
    std::size_t m = (n + 1) / 2;
    if (count > 0 && count < m) m = count;  // only the requested prefix
    for (std::size_t i = 0; i < m; ++i) {
        // Tricomi's asymptotic start, then Newton to the precision of T.
        T z = cos(T(pi * T(num_traits<T>::from_int(static_cast<long>(i) + 1) - num_traits<T>::from_rational(1, 4)) /
                    T(num_traits<T>::from_int(static_cast<long>(n)) + half)));
        T pp = one;
        for (int it = 0; it < 200; ++it) {
            T p1 = one, p2 = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < n; ++j) {
                const T p3 = p2;
                p2 = p1;
                const T jj = num_traits<T>::from_int(static_cast<long>(j) + 1);
                p1 = T(T(T(two * jj - one) * z * p2 - T(jj - one) * p3) / jj);
            }
            pp = T(num_traits<T>::from_int(static_cast<long>(n)) * T(z * p1 - p2) / T(z * z - one));
            const T dz = T(p1 / pp);
            z -= dz;
            if (num_traits<T>::to_double(num_abs(T(dz))) < 1e-40) break;
        }
        x[i] = T(xm - xl * z);
        x[n - 1 - i] = T(xm + xl * z);
        const T wi = T(two * xl / T(T(one - z * z) * pp * pp));
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
}

/**
 * n-point Gauss-Laguerre rule for weight exp(-x) on [0,inf), by Newton
 * iteration on the Laguerre polynomial (Numerical Recipes gaulag with
 * alpha = 0). The returned weights are the classical ones, i.e. they already
 * carry the exp(-x) factor, exactly like MATLAB's tabulated pair.
 */
template <class T>
void gauss_laguerre(std::size_t n, std::vector<T>& x, std::vector<T>& w) {
    static_assert(num_traits<T>::has_transcendental,
                  "gauss_laguerre requires transcendental arithmetic");
    using std::exp;
    using std::log;
    x.assign(n, num_traits<T>::from_int(0));
    w.assign(n, num_traits<T>::from_int(0));
    const T one = num_traits<T>::from_int(1);
    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    T z = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0) {
            z = num_traits<T>::from_double(3.0) / T(one + num_traits<T>::from_double(2.4) * nT);
        } else if (i == 1) {
            z += num_traits<T>::from_double(15.0) / T(one + num_traits<T>::from_double(2.5) * nT);
        } else {
            const T ai = num_traits<T>::from_int(static_cast<long>(i) - 1);
            z += T(one + num_traits<T>::from_double(2.55) * ai) /
                 T(num_traits<T>::from_double(1.9) * ai) * T(z - x[i - 2]);
        }
        T pp = one, p2 = num_traits<T>::from_int(0);
        for (int it = 0; it < 300; ++it) {
            T p1 = one;
            p2 = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < n; ++j) {
                const T p3 = p2;
                p2 = p1;
                const T jj = num_traits<T>::from_int(static_cast<long>(j) + 1);
                p1 = T(T(T(num_traits<T>::from_int(2 * static_cast<long>(j) + 1) - z) * p2 -
                         T(jj - one) * p3) /
                       jj);
            }
            // After the loop p1 = L_n(z) and p2 = L_{n-1}(z).
            pp = T(T(nT * p1 - nT * p2) / z);
            const T dz = T(p1 / pp);
            z -= dz;
            if (num_traits<T>::to_double(num_abs(T(dz))) < 1e-40) break;
        }
        x[i] = z;
        // Numerical Recipes weight-form rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        w[i] = T(-one / T(nT * pp * p2));
    }
}

// ---------------------------------------------------------------------------
// The reference's Inf in a field that may not have one
// ---------------------------------------------------------------------------

/**
 * MATLAB writes Inf for a rate or a balance value that is undefined (a zero
 * denominator, a NaN, a magnitude beyond 1e15). An exact rational field has no
 * infinity, so there the marker is 0 instead, which is unambiguous: a real
 * load-dependent service rate is never zero, and a station whose rate is zero
 * cannot serve, which is precisely what Inf encodes on the reciprocal side.
 * Callers must test with is_inf_marker rather than with isfinite.
 */
template <class T>
T num_inf_marker() {
    if constexpr (num_traits<T>::is_exact) {
        return num_traits<T>::from_int(0);
    } else {
        return num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    }
}

/** Companion test for num_inf_marker. */
template <class T>
bool is_inf_marker(const T& v) {
    if constexpr (num_traits<T>::is_exact) {
        return v == num_traits<T>::from_int(0);
    } else {
        return !std::isfinite(num_traits<T>::to_double(v));
    }
}

// ---------------------------------------------------------------------------
// Complex arithmetic over T, for the Norlund-Rice inversion integrands
// ---------------------------------------------------------------------------

/** Minimal complex over T. std::complex<T> is only defined for the built-ins. */
template <class T>
struct Cx {
    T re, im;
    Cx() : re(num_traits<T>::from_int(0)), im(num_traits<T>::from_int(0)) {}
    Cx(const T& r, const T& i) : re(r), im(i) {}
    explicit Cx(const T& r) : re(r), im(num_traits<T>::from_int(0)) {}
};

template <class T>
Cx<T> cx_add(const Cx<T>& a, const Cx<T>& b) {
    return Cx<T>(T(a.re + b.re), T(a.im + b.im));
}

template <class T>
Cx<T> cx_mul(const Cx<T>& a, const Cx<T>& b) {
    return Cx<T>(T(a.re * b.re - a.im * b.im), T(a.re * b.im + a.im * b.re));
}

template <class T>
Cx<T> cx_scale(const Cx<T>& a, const T& s) {
    return Cx<T>(T(a.re * s), T(a.im * s));
}

template <class T>
Cx<T> cx_div(const Cx<T>& a, const Cx<T>& b) {
    const T d = T(b.re * b.re + b.im * b.im);
    if (d == num_traits<T>::from_int(0)) throw NumericError("cx_div: division by zero");
    return Cx<T>(T(T(a.re * b.re + a.im * b.im) / d), T(T(a.im * b.re - a.re * b.im) / d));
}

/** exp(i theta). */
template <class T>
Cx<T> cx_expi(const T& theta) {
    static_assert(num_traits<T>::has_transcendental,
                  "cx_expi requires transcendental arithmetic");
    using std::cos;
    using std::sin;
    return Cx<T>(cos(theta), sin(theta));
}

}  // namespace detail
}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_ASYMPT_COMMON_H
