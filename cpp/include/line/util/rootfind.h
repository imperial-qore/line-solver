/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_ROOTFIND_H
#define LINE_UTIL_ROOTFIND_H

/**
 * Deterministic scalar root finding.
 *
 * This is the replacement for MATLAB's fsolve in the scalar cases the API
 * layer needs (the characteristic-time equations of the TTL cache
 * approximations, matlab/src/api/cache/cache_ttl_lrua.m and
 * cache_t_lrum_map.m, whose per-list residual is monotone in its own time).
 * fsolve is a trust-region method seeded from a caller-supplied -- in
 * cache_ttl_lrua.m, a *randomly generated* -- initial point, and its answer
 * therefore depends on rng state and on Optimization Toolbox availability.
 * Everything here is bracket based and deterministic: the same bracket and
 * tolerance always produce the same digits, with no global state, no toolbox,
 * and no fallback path that silently changes method.
 *
 * Three methods are provided:
 *   root_bisect      - bisection, one bit per iteration, cannot fail once the
 *                      bracket has a sign change
 *   root_brent       - Brent's method (inverse quadratic interpolation, secant
 *                      and bisection), superlinear but never worse than
 *                      bisection because every step is kept inside the bracket
 *   root_newton      - plain Newton from a starting point, for roots of even
 *                      multiplicity, which no bracketing method can see
 *   root_newton_safe - Newton safeguarded by a bracket: the Newton step is
 *                      taken only when it lands inside the current bracket and
 *                      reduces it, otherwise the step is a bisection
 *
 * ARITHMETIC: comparisons, addition, multiplication and division only, so
 * these instantiate at every backend including exact rationals. Note that a
 * root is still only located to the caller's tolerance: exact arithmetic makes
 * the iterates exact, not the answer, since an algebraic root need not be
 * rational at all.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {

/** Outcome of a scalar solve. */
template <class T>
struct RootResult {
    T root;             ///< best estimate of the root
    T value;            ///< f(root)
    T bracket_width;    ///< final |b - a|, zero for the unbracketed Newton
    unsigned iterations;///< iterations actually performed
    bool converged;     ///< tolerance was met before the iteration cap
};

/**
 * Bisection on a bracket with a sign change.
 *
 * @param f       callable T -> T
 * @param a,b     bracket endpoints, in either order
 * @param tol     absolute width of the final bracket
 * @param maxiter iteration cap
 * @throws InputError if f(a) and f(b) have the same sign and neither is a root
 */
template <class T, class F>
RootResult<T> root_bisect(F f, const T& a, const T& b, const T& tol, unsigned maxiter = 200) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    T lo = a < b ? a : b;
    T hi = a < b ? b : a;
    T flo = f(lo);
    T fhi = f(hi);

    RootResult<T> r;
    r.iterations = 0;
    r.bracket_width = T(hi - lo);
    if (flo == zero) {
        r.root = lo;
        r.value = flo;
        r.converged = true;
        return r;
    }
    if (fhi == zero) {
        r.root = hi;
        r.value = fhi;
        r.converged = true;
        return r;
    }
    if ((flo > zero) == (fhi > zero))
        throw InputError("root_bisect: the bracket endpoints do not straddle a root");

    for (unsigned it = 0; it < maxiter; ++it) {
        r.iterations = it + 1;
        const T mid = (lo + hi) / two;
        const T fm = f(mid);
        if (fm == zero) {
            r.root = mid;
            r.value = fm;
            r.bracket_width = zero;
            r.converged = true;
            return r;
        }
        if ((fm > zero) == (flo > zero)) {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
            fhi = fm;
        }
        if (T(hi - lo) <= tol) break;
    }
    const T mid = (lo + hi) / two;
    r.root = mid;
    r.value = f(mid);
    r.bracket_width = T(hi - lo);
    r.converged = T(hi - lo) <= tol;
    return r;
}

/**
 * Brent's method on a bracket with a sign change. Falls back to bisection
 * whenever the interpolated step is not a strict improvement, so the bracket
 * is never lost.
 */
template <class T, class F>
RootResult<T> root_brent(F f, const T& a0, const T& b0, const T& tol, unsigned maxiter = 200) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T half = num_traits<T>::from_rational(1, 2);

    T a = a0;
    T b = b0;
    T fa = f(a);
    T fb = f(b);

    RootResult<T> r;
    r.iterations = 0;
    r.bracket_width = num_abs(T(b - a));
    if (fa == zero) {
        r.root = a;
        r.value = fa;
        r.converged = true;
        return r;
    }
    if (fb == zero) {
        r.root = b;
        r.value = fb;
        r.converged = true;
        return r;
    }
    if ((fa > zero) == (fb > zero))
        throw InputError("root_brent: the bracket endpoints do not straddle a root");

    if (num_abs(T(fa)) < num_abs(T(fb))) {
        T t = a;
        a = b;
        b = t;
        t = fa;
        fa = fb;
        fb = t;
    }
    T c = a;
    T fc = fa;
    T d = T(b - a);
    T e = d;
    bool used_bisect = true;

    for (unsigned it = 0; it < maxiter; ++it) {
        r.iterations = it + 1;
        if ((fb > zero) == (fc > zero)) {
            c = a;
            fc = fa;
            d = T(b - a);
            e = d;
        }
        if (num_abs(T(fc)) < num_abs(T(fb))) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const T m = half * T(c - b);
        if (num_abs(T(m)) <= tol || fb == zero) break;

        bool interpolate = false;
        T s = zero, p = zero, q = zero;
        if (num_abs(T(e)) >= tol && num_abs(T(fa)) > num_abs(T(fb))) {
            s = fb / fa;
            if (a == c) {  // secant
                p = two * m * s;
                q = num_traits<T>::from_int(1) - s;
            } else {  // inverse quadratic
                const T qq = fa / fc;
                const T rr = fb / fc;
                p = s * (two * m * qq * (qq - rr) - T(b - a) * (rr - num_traits<T>::from_int(1)));
                q = (qq - num_traits<T>::from_int(1)) * (rr - num_traits<T>::from_int(1)) *
                    (s - num_traits<T>::from_int(1));
            }
            if (p > zero)
                q = -q;
            else
                p = -p;
            const T lim1 = three * m * q - num_abs(T(tol * q));
            const T lim2 = num_abs(T(e * q));
            interpolate = two * p < (lim1 < lim2 ? lim1 : lim2);
        }
        e = interpolate ? d : m;
        d = interpolate ? p / q : m;
        used_bisect = !interpolate;
        (void)used_bisect;

        a = b;
        fa = fb;
        if (num_abs(T(d)) > tol)
            b += d;
        else
            b += (m > zero ? tol : T(-tol));
        fb = f(b);
    }
    r.root = b;
    r.value = fb;
    r.bracket_width = num_abs(T(c - b));
    r.converged = num_abs(T(half * T(c - b))) <= tol || fb == zero;
    return r;
}

/**
 * Plain Newton from a starting point. The only method here that can find a
 * root of even multiplicity, where f does not change sign and no bracket
 * exists; convergence is then linear rather than quadratic.
 *
 * @param f  callable T -> T
 * @param df callable T -> T, the derivative
 * @param x0 starting point of the iteration
 * @param tol convergence tolerance on the Newton step
 * @param maxiter iteration cap (default 200)
 * @throws NumericError if the derivative vanishes at an iterate
 */
template <class T, class F, class DF>
RootResult<T> root_newton(F f, DF df, const T& x0, const T& tol, unsigned maxiter = 200) {
    const T zero = num_traits<T>::from_int(0);
    T x = x0;
    RootResult<T> r;
    r.iterations = 0;
    r.bracket_width = zero;
    r.converged = false;
    for (unsigned it = 0; it < maxiter; ++it) {
        r.iterations = it + 1;
        const T fx = f(x);
        if (fx == zero) {
            r.root = x;
            r.value = fx;
            r.converged = true;
            return r;
        }
        const T dfx = df(x);
        if (dfx == zero) throw NumericError("root_newton: the derivative vanished at an iterate");
        const T step = fx / dfx;
        x -= step;
        if (num_abs(T(step)) <= tol) {
            r.root = x;
            r.value = f(x);
            r.converged = true;
            return r;
        }
    }
    r.root = x;
    r.value = f(x);
    return r;
}

/**
 * Newton safeguarded by a bracket with a sign change: the Newton step is used
 * only when it stays inside the bracket and at least halves it, otherwise the
 * step is a bisection. Never diverges and never leaves the bracket.
 */
template <class T, class F, class DF>
RootResult<T> root_newton_safe(F f, DF df, const T& a0, const T& b0, const T& tol,
                               unsigned maxiter = 200) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    T lo = a0 < b0 ? a0 : b0;
    T hi = a0 < b0 ? b0 : a0;
    T flo = f(lo);
    T fhi = f(hi);

    RootResult<T> r;
    r.iterations = 0;
    r.bracket_width = T(hi - lo);
    if (flo == zero) {
        r.root = lo;
        r.value = flo;
        r.converged = true;
        return r;
    }
    if (fhi == zero) {
        r.root = hi;
        r.value = fhi;
        r.converged = true;
        return r;
    }
    if ((flo > zero) == (fhi > zero))
        throw InputError("root_newton_safe: the bracket endpoints do not straddle a root");

    T x = (lo + hi) / two;
    for (unsigned it = 0; it < maxiter; ++it) {
        r.iterations = it + 1;
        const T fx = f(x);
        if (fx == zero) {
            r.root = x;
            r.value = fx;
            r.bracket_width = zero;
            r.converged = true;
            return r;
        }
        if ((fx > zero) == (flo > zero)) {
            lo = x;
            flo = fx;
        } else {
            hi = x;
            fhi = fx;
        }
        const T width = T(hi - lo);
        if (width <= tol) break;

        const T dfx = df(x);
        bool take_newton = false;
        T xn = x;
        if (!(dfx == zero)) {
            xn = x - fx / dfx;
            take_newton = xn > lo && xn < hi;
        }
        x = take_newton ? xn : (lo + hi) / two;
    }
    r.root = x;
    r.value = f(x);
    r.bracket_width = T(hi - lo);
    r.converged = T(hi - lo) <= tol;
    return r;
}

/**
 * Expand a bracket to the right until f changes sign, doubling the upper end.
 * Used by the TTL cache fixed points, where the residual is monotone in the
 * characteristic time but no upper bound is known a priori.
 *
 * @throws NumericError if no sign change is found before the cap
 */
template <class T, class F>
void bracket_expand(F f, const T& a, T& b, unsigned maxdoubling = 200) {
    const T zero = num_traits<T>::from_int(0);
    const T two = num_traits<T>::from_int(2);
    const T fa = f(a);
    if (fa == zero) {
        b = a;
        return;
    }
    for (unsigned it = 0; it < maxdoubling; ++it) {
        const T fb = f(b);
        if (fb == zero || (fb > zero) != (fa > zero)) return;
        b *= two;
    }
    throw NumericError("bracket_expand: no sign change found while expanding the bracket");
}

}  // namespace line

#endif  // LINE_UTIL_ROOTFIND_H
