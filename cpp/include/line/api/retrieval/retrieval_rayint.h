/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_RAYINT_H
#define LINE_API_RETRIEVAL_RETRIEVAL_RAYINT_H

/**
 * Ray (WKB) asymptotic expansion of the list-based cache normalizing constant.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_rayint.m, cross-checked
 * against jar/src/main/java/jline/api/retrieval/Retrieval_rayint.java.
 *
 * Approximates the constant cache_erec computes exactly, in the SAME
 * normalization, so the two are interchangeable:
 *
 *   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),   E(0,0)=1
 *
 * Writing E = prod_j m_j! * Et, the relaxation Et ~ H exp(phi/eps) with
 * n = y/eps and m_j = x_j/eps gives the eikonal
 * e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j}, whose rays carry the constants
 * xi_j = e^{-phi_j}. With S(v) = 1 + sum_j gamma_j(v) xi_j,
 *
 *   x_j    = int_0^y gamma_j(v) xi_j / S(v) dv          (the saddle conditions)
 *   phi    = int_0^y log S(v) dv - sum_j x_j log xi_j
 *   H      = (2 pi)^{-h/2} sqrt(S(y)/S(0)) / sqrt(prod_j xi_j * det A)
 *   A_{ik} = d x_i / d xi_k
 *
 * and E ~ prod_j m_j! * eps^{h/2} H exp(phi/eps).
 *
 * DISCRETE (gamma an n x h matrix). The ray integrals are the sums they
 * discretize and the expansion collapses to the Laplace form
 *
 *   Et ~ (2 pi)^{-h/2} exp(sum_k log D_k - sum_j m_j log xi_j) / sqrt(det Sigma)
 *
 * with D_k = 1 + sum_j gamma_{k,j} xi_j, sum_k gamma_{k,j} xi_j / D_k = m_j and
 * Sigma = A * diag(xi) the Hessian in log xi. This is the more accurate of the
 * two forms; the sqrt(S(y)/S(0)) factor is exactly the Euler-Maclaurin term
 * relating sum_k to int dv and is already accounted for.
 *
 * CONTINUUM (gamma a callable profile on v in [0,1]). Composite Simpson
 * quadrature on the profile itself, the form written in the note. Costs roughly
 * a factor two in accuracy but does not need the n rows.
 *
 * ACCURACY. The relative error is O(1/n) at fixed occupancy but is governed by
 * the smallest occupancy rather than by n, tracking
 * 0.14 * (1/min_j m_j + 1/(n - sum_j m_j)), so a per cent needs every m_j and
 * n - sum_j m_j above about 15 and a part in a thousand needs them above about
 * 150. Returned in the result's relerr_est. Lists with m_j = 0 contribute
 * nothing and are dropped before the saddle is solved.
 *
 * ARITHMETIC: the expansion is a Laplace approximation built out of logs, exps
 * and a square root, so it is meaningless at exact arithmetic and is gated on
 * has_transcendental. It is also an APPROXIMATION whatever the arithmetic --
 * widening the type sharpens the saddle solve, never the O(1/n) model error.
 * Use cache_erec or retrieval_nc when the exact constant is wanted.
 *
 * This is the no-fetch (q=0) case, i.e. the same quantity as cache_erec. The
 * delayed-hit extension carrying the fetch coordinates is NOT implemented: its
 * eikonal is known but its amplitude has not been derived.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** Outcome of the expansion. */
template <class T>
struct RetrievalRayintResult {
    /** Normalizing constant, same normalization as cache_erec (may overflow; use log_e). */
    T e;
    /** Natural logarithm of e, safe for large n. */
    T log_e;
    /** Saddle point xi_j, one entry per list (0 for a list of zero capacity). */
    std::vector<T> xi;
    /** log det of the Hessian in log xi. */
    T logdet_sigma;
    /** Estimated relative error, 0.14*(1/min_j m_j + 1/(n - sum_j m_j)). */
    double relerr_est = 0.0;
    /** "saddle", "rayint" or "boundary". */
    std::string method;
    /** Newton iterations used. */
    std::size_t iterations = 0;
};

namespace detail {

/** Hessian in theta = log xi; equals A*diag(xi) with A_{ik} = d x_i / d xi_k. */
template <class T>
Matrix<T> rayint_hessian(const Matrix<T>& g, const std::vector<T>& w, const std::vector<T>& xi) {
    const std::size_t K = g.rows();
    const std::size_t h = xi.size();
    const T one = num_traits<T>::from_int(1);
    Matrix<T> hess(h, h, num_traits<T>::from_int(0));
    std::vector<T> a(h);
    for (std::size_t k = 0; k < K; ++k) {
        T s = one;
        for (std::size_t j = 0; j < h; ++j) {
            a[j] = g(k, j) * xi[j];
            s += a[j];
        }
        for (std::size_t j = 0; j < h; ++j) {
            const T aj = a[j] / s;
            hess(j, j) += w[k] * aj;
            for (std::size_t l = 0; l < h; ++l) hess(j, l) -= w[k] * aj * (a[l] / s);
        }
    }
    return hess;
}

/** Gaussian elimination with partial pivoting; h is the list count, so tiny. */
template <class T>
std::vector<T> rayint_solve(const Matrix<T>& a, const std::vector<T>& b) {
    using std::abs;
    const std::size_t h = b.size();
    Matrix<T> m(h, h + 1);
    for (std::size_t i = 0; i < h; ++i) {
        for (std::size_t j = 0; j < h; ++j) m(i, j) = a(i, j);
        m(i, h) = b[i];
    }
    for (std::size_t c = 0; c < h; ++c) {
        std::size_t piv = c;
        for (std::size_t i = c + 1; i < h; ++i)
            if (abs(m(i, c)) > abs(m(piv, c))) piv = i;
        if (num_traits<T>::to_double(abs(m(piv, c))) == 0.0)
            throw NumericError("retrieval_rayint: the saddle-point Hessian is singular; the ray map "
                               "is degenerate here");
        if (piv != c)
            for (std::size_t j = 0; j <= h; ++j) {
                const T tmp = m(c, j);
                m(c, j) = m(piv, j);
                m(piv, j) = tmp;
            }
        for (std::size_t i = c + 1; i < h; ++i) {
            const T f = m(i, c) / m(c, c);
            for (std::size_t j = c; j <= h; ++j) m(i, j) -= f * m(c, j);
        }
    }
    std::vector<T> x(h);
    for (std::size_t ii = h; ii-- > 0;) {
        T s = m(ii, h);
        for (std::size_t j = ii + 1; j < h; ++j) s -= m(ii, j) * x[j];
        x[ii] = s / m(ii, ii);
    }
    return x;
}

/** log det via Cholesky; the Hessian of a strictly convex objective is positive definite. */
template <class T>
T rayint_logdet(const Matrix<T>& a) {
    using std::log;
    using std::sqrt;
    const std::size_t h = a.rows();
    const T half = num_traits<T>::from_rational(1, 2);
    const T two = num_traits<T>::from_int(2);
    Matrix<T> l(h, h, num_traits<T>::from_int(0));
    T ld = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < h; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T s = half * (a(i, j) + a(j, i));
            for (std::size_t k = 0; k < j; ++k) s -= l(i, k) * l(j, k);
            if (i == j) {
                if (num_traits<T>::to_double(s) <= 0.0)
                    throw NumericError("retrieval_rayint: the saddle-point Hessian is not positive "
                                       "definite; the ray map is singular here");
                l(i, j) = sqrt(s);
                ld += two * log(l(i, j));
            } else {
                l(i, j) = s / l(j, j);
            }
        }
    }
    return ld;
}

/**
 * Newton on theta = log xi for sum_k w_k g_{k,j} xi_j / (1 + sum_l g_{k,l} xi_l) = tgt_j.
 * The objective sum_k w_k log(1 + sum_l g_{k,l} e^{theta_l}) - tgt.theta is strictly
 * convex, so the root is unique and damped Newton converges globally.
 */
template <class T>
std::vector<T> rayint_saddle(const Matrix<T>& g, const std::vector<T>& w,
                             const std::vector<T>& tgt, std::size_t& iterations) {
    using std::abs;
    using std::exp;
    using std::log;
    const std::size_t K = g.rows();
    const std::size_t h = tgt.size();
    const T one = num_traits<T>::from_int(1);

    T wsum = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < K; ++k) wsum += w[k];
    T tsum = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < h; ++j) tsum += tgt[j];
    T slack = one - tsum / wsum;
    const T tiny = num_traits<T>::from_double(1e-9);
    if (num_traits<T>::to_double(slack) < 1e-9) slack = tiny;

    std::vector<T> th(h);
    for (std::size_t j = 0; j < h; ++j) {
        T gb = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < K; ++k) gb += w[k] * g(k, j);
        th[j] = log(tgt[j] / (gb * slack));
    }

    double tmax = 1.0;
    for (std::size_t j = 0; j < h; ++j)
        tmax = std::max(tmax, std::abs(num_traits<T>::to_double(tgt[j])));

    std::vector<T> xi(h);
    std::size_t it = 0;
    for (it = 1; it <= 200; ++it) {
        for (std::size_t j = 0; j < h; ++j) xi[j] = exp(th[j]);
        std::vector<T> grad(h, num_traits<T>::from_int(0));
        for (std::size_t k = 0; k < K; ++k) {
            T s = one;
            for (std::size_t j = 0; j < h; ++j) s += g(k, j) * xi[j];
            for (std::size_t j = 0; j < h; ++j) grad[j] += w[k] * g(k, j) * xi[j] / s;
        }
        double gmax = 0.0;
        for (std::size_t j = 0; j < h; ++j) {
            grad[j] -= tgt[j];
            gmax = std::max(gmax, std::abs(num_traits<T>::to_double(grad[j])));
        }
        if (gmax <= 1e-12 * tmax) break;
        const Matrix<T> hess = rayint_hessian(g, w, xi);
        std::vector<T> d = rayint_solve(hess, grad);
        double dmax = 0.0;
        for (std::size_t j = 0; j < h; ++j) {
            d[j] = -d[j];
            dmax = std::max(dmax, std::abs(num_traits<T>::to_double(d[j])));
        }
        T step = one;
        double stepd = 1.0;
        while (stepd * dmax > 2.0) {
            stepd /= 2.0;
            step = step / num_traits<T>::from_int(2);
        }
        for (std::size_t j = 0; j < h; ++j) th[j] += step * d[j];
    }
    iterations = it;
    for (std::size_t j = 0; j < h; ++j) xi[j] = exp(th[j]);
    return xi;
}

/** log(x!) for a non-negative integer x, by summing logs; the counts here are small. */
template <class T>
T rayint_logfact(int x) {
    using std::log;
    T s = num_traits<T>::from_int(0);
    for (int i = 2; i <= x; ++i) s += log(num_traits<T>::from_int(i));
    return s;
}

/** Shared core: `gmat` non-null selects the discrete form, otherwise the quadrature one. */
template <class T>
RetrievalRayintResult<T> rayint_run(const Matrix<T>* gmat,
                                    const std::function<Matrix<T>(const std::vector<T>&)>* gfun,
                                    int n, const std::vector<int>& m, std::size_t nquad) {
    using std::log;
    using std::sqrt;
    const std::size_t h = m.size();
    if (h == 0) throw InputError("retrieval_rayint: the capacity vector must not be empty");

    long msum = 0;
    for (std::size_t j = 0; j < h; ++j) {
        if (m[j] < 0) throw InputError("retrieval_rayint: list capacities must be non-negative");
        msum += m[j];
    }
    if (gmat != nullptr && gmat->cols() != h)
        throw InputError("retrieval_rayint: gamma and m disagree on the number of lists");

    RetrievalRayintResult<T> out;
    out.xi.assign(h, num_traits<T>::from_int(0));
    out.logdet_sigma = num_traits<T>::from_int(0);

    if (msum > n) {
        out.e = num_traits<T>::from_int(0);
        out.log_e = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
        out.method = "boundary";
        return out;
    }
    if (msum == 0) {
        out.e = num_traits<T>::from_int(1);
        out.log_e = num_traits<T>::from_int(0);
        out.method = "boundary";
        return out;
    }
    if (msum == n)
        throw InputError("retrieval_rayint: the expansion requires sum(m) < n; at sum(m) = n the "
                         "saddle point escapes to infinity, use cache_erec for a full cache");

    // lists of zero capacity contribute nothing and would make the saddle singular
    std::vector<std::size_t> keep;
    for (std::size_t j = 0; j < h; ++j)
        if (m[j] > 0) keep.push_back(j);
    const std::size_t hk = keep.size();

    Matrix<T> g;
    std::vector<T> w;
    std::vector<T> tgt(hk);
    const bool discrete = (gmat != nullptr);
    if (discrete) {
        const std::size_t K = gmat->rows();
        g = Matrix<T>(K, hk);
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t a = 0; a < hk; ++a) g(k, a) = (*gmat)(k, keep[a]);
        w.assign(K, num_traits<T>::from_int(1));
        for (std::size_t a = 0; a < hk; ++a) tgt[a] = num_traits<T>::from_int(m[keep[a]]);
    } else {
        const std::size_t K = nquad;
        std::vector<T> v(K);
        for (std::size_t k = 0; k < K; ++k)
            v[k] = num_traits<T>::from_rational(static_cast<long>(k), static_cast<long>(K - 1));
        const Matrix<T> gfull = (*gfun)(v);
        if (gfull.rows() != K)
            throw InputError("retrieval_rayint: the profile must return one row per evaluation point");
        if (gfull.cols() != h)
            throw InputError("retrieval_rayint: the profile must return one column per cache list");
        g = Matrix<T>(K, hk);
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t a = 0; a < hk; ++a) g(k, a) = gfull(k, keep[a]);
        w.assign(K, num_traits<T>::from_int(0));
        for (std::size_t k = 0; k < K; ++k) {
            const long c = (k == 0 || k == K - 1) ? 1 : ((k % 2 == 1) ? 4 : 2);
            w[k] = num_traits<T>::from_rational(c, static_cast<long>(3 * (K - 1)));
        }
        for (std::size_t a = 0; a < hk; ++a)
            tgt[a] = num_traits<T>::from_rational(m[keep[a]], n);
    }

    std::size_t iters = 0;
    const std::vector<T> xi = detail::rayint_saddle(g, w, tgt, iters);

    const std::size_t K = g.rows();
    const T one = num_traits<T>::from_int(1);
    std::vector<T> s(K);
    for (std::size_t k = 0; k < K; ++k) {
        T acc = one;
        for (std::size_t a = 0; a < hk; ++a) acc += g(k, a) * xi[a];
        s[k] = acc;
    }
    const Matrix<T> hess = detail::rayint_hessian(g, w, xi);
    const T logdet = detail::rayint_logdet(hess);

    const T half = num_traits<T>::from_rational(1, 2);
    const T log2pi = num_traits<T>::from_double(std::log(2.0 * M_PI));
    const T hkT = num_traits<T>::from_int(static_cast<long>(hk));

    T log_et;
    if (discrete) {
        T phi = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < K; ++k) phi += log(s[k]);
        for (std::size_t a = 0; a < hk; ++a) phi -= tgt[a] * log(xi[a]);
        log_et = -half * hkT * log2pi + phi - half * logdet;
        out.method = "saddle";
    } else {
        T phi = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < K; ++k) phi += w[k] * log(s[k]);
        for (std::size_t a = 0; a < hk; ++a) phi -= tgt[a] * log(xi[a]);
        const T nT = num_traits<T>::from_int(n);
        log_et = -half * hkT * log(nT) - half * hkT * log2pi + nT * phi - half * logdet +
                 half * log(s[K - 1] / s[0]);
        out.method = "rayint";
    }

    T logfact = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < h; ++j) logfact += detail::rayint_logfact<T>(m[j]);
    out.log_e = log_et + logfact;
    {
        using std::exp;
        out.e = exp(out.log_e);
    }
    for (std::size_t a = 0; a < hk; ++a) out.xi[keep[a]] = xi[a];
    out.logdet_sigma = logdet;
    out.iterations = iters;
    int mmin = m[keep[0]];
    for (std::size_t a = 1; a < hk; ++a) mmin = std::min(mmin, m[keep[a]]);
    out.relerr_est = 0.14 * (1.0 / static_cast<double>(mmin) +
                             1.0 / static_cast<double>(n - msum));
    return out;
}

}  // namespace detail

/**
 * Discrete (saddle) form.
 *
 * @param gamma access factors gamma(k,j), n x h
 * @param m     cache list capacities, length h
 */
template <class T>
RetrievalRayintResult<T> retrieval_rayint(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "retrieval_rayint is a Laplace approximation built out of logs, exps and a "
                  "square root: it is meaningless at exact arithmetic, and widening the type "
                  "sharpens the saddle solve but never the O(1/n) model error. Use cache_erec "
                  "or retrieval_nc for the exact constant");
    if (gamma.rows() == 0) throw InputError("retrieval_rayint: gamma must be a non-empty n x h matrix");
    return detail::rayint_run<T>(&gamma, nullptr, static_cast<int>(gamma.rows()), m, 0);
}

/**
 * Continuum (ray-integral) form.
 *
 * @param gfun  access-factor profile on v in [0,1], returning v.size() x h
 * @param m     cache list capacities, length h
 * @param n     number of items
 * @param nquad composite Simpson nodes (forced odd, at least 5)
 */
template <class T>
RetrievalRayintResult<T> retrieval_rayint(const std::function<Matrix<T>(const std::vector<T>&)>& gfun,
                                          const std::vector<int>& m, int n,
                                          std::size_t nquad = 4097) {
    static_assert(num_traits<T>::has_transcendental,
                  "retrieval_rayint is a Laplace approximation built out of logs, exps and a "
                  "square root: it is meaningless at exact arithmetic. Use cache_erec or "
                  "retrieval_nc for the exact constant");
    const std::size_t nq = std::max<std::size_t>(5, 2 * (nquad / 2) + 1);
    return detail::rayint_run<T>(nullptr, &gfun, n, m, nq);
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_RAYINT_H
