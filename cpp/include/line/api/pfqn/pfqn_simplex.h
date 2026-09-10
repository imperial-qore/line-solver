/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_SIMPLEX_H
#define LINE_API_PFQN_PFQN_SIMPLEX_H

/**
 * Shared machinery for closures of the simplex factor of the McKenna-Mitra integral,
 * used by pfqn_aghq.h.
 *
 * At Z = 0 the integrand is homogeneous of degree sum(N), so y = v*x separates and the
 * radius integrates exactly to gamma(N+M), leaving an integral over the unit simplex in
 * which ALL of the error of the logistic expansion lives. With Z > 0 that factorisation
 * is gone: the radius cannot be marginalised and is integrated numerically here rather
 * than closed, leaving the same M-1 simplex directions to a closure.
 *
 * ARITHMETIC. Everything here is a Laplace-type approximation or a quadrature of a
 * transcendental integrand, so each entry point is gated on
 * num_traits<T>::has_transcendental exactly as pfqn_le is.
 *
 * EIGENSOLVER. util/eig.h is LAPACK and double-only, while this family is templated on
 * T, so a cyclic Jacobi eigensolver is carried here instead. It is used for the
 * Golub-Welsch construction of the quadrature nodes and for the principal-axis frame of
 * the adaptive Gauss-Hermite rule. Jacobi is chosen over a tridiagonal QL because the
 * matrices are small, it needs no shift strategy, and it is symmetric-exact by
 * construction.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_le.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {
namespace simplex {

/**
 * Eigenvalues and eigenvectors of a symmetric matrix by cyclic Jacobi, ascending.
 * V.col(k) is the unit eigenvector of d[k].
 */
template <class T>
void sym_eig(Matrix<T> A, std::vector<T>& d, Matrix<T>& V) {
    using std::abs;
    using std::sqrt;
    const std::size_t n = A.rows();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    V = Matrix<T>(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) V(i, j) = (i == j) ? one : zero;
    for (int sweep = 0; sweep < 100; ++sweep) {
        T off = zero;
        for (std::size_t p = 0; p + 1 < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) off += A(p, q) * A(p, q);
        if (num_traits<T>::to_double(off) <= 1e-30) break;
        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (num_traits<T>::to_double(abs(A(p, q))) == 0.0) continue;
                T theta = T((A(q, q) - A(p, p)) / (A(p, q) + A(p, q)));
                T t = T(one / T(abs(theta) + sqrt(T(theta * theta + one))));
                if (num_traits<T>::to_double(theta) < 0.0) t = T(zero - t);
                T c = T(one / sqrt(T(t * t + one)));
                T s = T(t * c);
                for (std::size_t k = 0; k < n; ++k) {
                    T akp = A(k, p), akq = A(k, q);
                    A(k, p) = T(c * akp - s * akq);
                    A(k, q) = T(s * akp + c * akq);
                }
                for (std::size_t k = 0; k < n; ++k) {
                    T apk = A(p, k), aqk = A(q, k);
                    A(p, k) = T(c * apk - s * aqk);
                    A(q, k) = T(s * apk + c * aqk);
                }
                for (std::size_t k = 0; k < n; ++k) {
                    T vkp = V(k, p), vkq = V(k, q);
                    V(k, p) = T(c * vkp - s * vkq);
                    V(k, q) = T(s * vkp + c * vkq);
                }
            }
        }
    }
    d.assign(n, zero);
    for (std::size_t i = 0; i < n; ++i) d[i] = A(i, i);
    for (std::size_t i = 0; i + 1 < n; ++i) {  // selection sort, ascending
        std::size_t m = i;
        for (std::size_t j = i + 1; j < n; ++j)
            if (num_traits<T>::to_double(d[j]) < num_traits<T>::to_double(d[m])) m = j;
        if (m != i) {
            T tmp = d[i];
            d[i] = d[m];
            d[m] = tmp;
            for (std::size_t k = 0; k < n; ++k) {
                T v = V(k, i);
                V(k, i) = V(k, m);
                V(k, m) = v;
            }
        }
    }
}

/** Nodes and weights of the Gauss rule with the given Jacobi off-diagonal. */
template <class T>
void golub_welsch(const std::vector<T>& off, const T& mu0, std::vector<T>& x,
                  std::vector<T>& w) {
    const std::size_t n = off.size() + 1;
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> J(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) J(i, j) = zero;
    for (std::size_t k = 0; k + 1 < n; ++k) {
        J(k, k + 1) = off[k];
        J(k + 1, k) = off[k];
    }
    Matrix<T> V(n, n);
    sym_eig(J, x, V);
    w.assign(n, zero);
    for (std::size_t k = 0; k < n; ++k) w[k] = T(mu0 * V(0, k) * V(0, k));
}

/** N-point Gauss-Legendre rule on [-1,1]. */
template <class T>
void gauss_legendre(std::size_t n, std::vector<T>& x, std::vector<T>& w) {
    using std::sqrt;
    std::vector<T> off(n - 1);
    for (std::size_t k = 1; k < n; ++k) {
        T kk = num_traits<T>::from_int(static_cast<long>(k));
        off[k - 1] = T(kk / sqrt(T(num_traits<T>::from_int(4) * kk * kk -
                                   num_traits<T>::from_int(1))));
    }
    golub_welsch(off, num_traits<T>::from_int(2), x, w);
}

/** Q-point Gauss-Hermite rule of the probabilists' weight exp(-z^2/2). */
template <class T>
void gauss_hermite(std::size_t q, std::vector<T>& z, std::vector<T>& w) {
    using std::sqrt;
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    if (q == 1) {
        z.assign(1, num_traits<T>::from_int(0));
        w.assign(1, sqrt(twopi));
        return;
    }
    std::vector<T> off(q - 1);
    for (std::size_t k = 1; k < q; ++k)
        off[k - 1] = sqrt(num_traits<T>::from_int(static_cast<long>(k)));
    golub_welsch(off, T(sqrt(twopi)), z, w);
}

/** Log-integrand of the radial integral in t = log v, Jacobian included. */
template <class T>
T radial_logf(const T& t, const std::vector<T>& c, const std::vector<T>& N,
              const std::vector<T>& Z, std::size_t M) {
    using std::exp;
    using std::log;
    const T tiny = num_traits<T>::from_double(2.2250738585072014e-308);
    T v = exp(t);
    T f = T(num_traits<T>::from_int(0) - v +
            num_traits<T>::from_int(static_cast<long>(M)) * t);
    for (std::size_t r = 0; r < c.size(); ++r) {
        T d = T(Z[r] + v * c[r]);
        if (num_traits<T>::to_double(d) < 2.2250738585072014e-308) d = tiny;
        f += N[r] * log(d);
    }
    return f;
}

/** log J(c) and the moments of the tilted law of the radius. */
template <class T>
struct Radial {
    T lJ;
    std::vector<T> G;
    T vbar;
    Matrix<T> Lam;
};

/**
 * log J(c) = log int_0^inf exp(-v) v^(M-1) prod_r (Z_r + v c_r)^N_r dv, plus the moments
 * of the tilted law of v that the simplex derivatives need: G_r = E[T_r], vbar = E[v]
 * and Lam = cov(T) - diag(E[T^2]/N) = grad^2_c log J, with T_r(v) = N_r v/(Z_r + v c_r).
 * Quadrature runs in t = log v, where the integrand is bounded at both ends, over two
 * Gauss-Legendre panels meeting at the mode, each widened until the log-integrand has
 * fallen 60 nats so the discarded tails are below 1e-26 in relative terms.
 */
template <class T>
Radial<T> radial(const std::vector<T>& c, const std::vector<T>& N, const std::vector<T>& Z,
                 std::size_t M) {
    static_assert(num_traits<T>::has_transcendental,
                  "radial requires transcendental arithmetic (quadrature of an integral)");
    using std::abs;
    using std::exp;
    using std::log;
    using std::sqrt;
    const std::size_t R = c.size();
    const T zero = num_traits<T>::from_int(0);
    const T tiny = num_traits<T>::from_double(2.2250738585072014e-308);
    const T Md = num_traits<T>::from_int(static_cast<long>(M));

    static std::vector<T> vg, wg;
    if (vg.empty()) gauss_legendre<T>(64, vg, wg);

    T Ntot = zero;
    for (std::size_t r = 0; r < R; ++r) Ntot += N[r];
    T t = log(T(Ntot + Md));
    for (int it = 0; it < 200; ++it) {
        T v = exp(t);
        T f1 = T(Md - v), f2 = T(zero - v);
        for (std::size_t r = 0; r < R; ++r) {
            T d = T(Z[r] + v * c[r]);
            if (num_traits<T>::to_double(d) < 2.2250738585072014e-308) d = tiny;
            f1 += N[r] * T(v * c[r]) / d;
            f2 += N[r] * T(v * c[r]) * Z[r] / T(d * d);
        }
        if (num_traits<T>::to_double(f2) > -1e-300) break;
        double stepd = num_traits<T>::to_double(T(zero - f1) / f2);
        if (stepd > 2.0) stepd = 2.0;
        if (stepd < -2.0) stepd = -2.0;
        t += num_traits<T>::from_double(stepd);
        if (stepd < 1e-13 && stepd > -1e-13) break;
    }
    T v = exp(t);
    T f2 = T(zero - v);
    for (std::size_t r = 0; r < R; ++r) {
        T d = T(Z[r] + v * c[r]);
        if (num_traits<T>::to_double(d) < 2.2250738585072014e-308) d = tiny;
        f2 += N[r] * T(v * c[r]) * Z[r] / T(d * d);
    }
    T sig = num_traits<T>::from_int(1);
    if (num_traits<T>::to_double(f2) < -1e-300) sig = T(num_traits<T>::from_int(1) /
                                                        sqrt(T(zero - f2)));
    T fm = radial_logf(t, c, N, Z, M);
    T lim = T(t + num_traits<T>::from_double(745.0));
    T a = T(num_traits<T>::from_double(12.0) * sig);
    if (num_traits<T>::to_double(a) > num_traits<T>::to_double(lim)) a = lim;
    for (int k = 0; k < 60; ++k) {
        if (num_traits<T>::to_double(T(t - a)) <= -745.0) break;
        if (num_traits<T>::to_double(radial_logf(T(t - a), c, N, Z, M)) <
            num_traits<T>::to_double(fm) - 60.0)
            break;
        a = T(num_traits<T>::from_double(1.6) * a);
        if (num_traits<T>::to_double(a) > num_traits<T>::to_double(lim)) a = lim;
    }
    T b = T(num_traits<T>::from_double(12.0) * sig);
    for (int k = 0; k < 60; ++k) {
        if (num_traits<T>::to_double(radial_logf(T(t + b), c, N, Z, M)) <
            num_traits<T>::to_double(fm) - 60.0)
            break;
        b = T(num_traits<T>::from_double(1.6) * b);
    }
    const T half = num_traits<T>::from_double(0.5);
    const std::size_t nq = 2 * vg.size();
    std::vector<T> tt(nq), W(nq), vv(nq), fv(nq);
    for (std::size_t k = 0; k < vg.size(); ++k) {
        tt[k] = T(half * a * vg[k] + T(t - half * a));
        W[k] = T(half * a * wg[k]);
        tt[vg.size() + k] = T(half * b * vg[k] + T(t + half * b));
        W[vg.size() + k] = T(half * b * wg[k]);
    }
    Matrix<T> D(nq, R);
    double mx = -1e308;
    for (std::size_t k = 0; k < nq; ++k) {
        vv[k] = exp(tt[k]);
        T f = T(T(zero - vv[k]) + Md * tt[k]);
        for (std::size_t r = 0; r < R; ++r) {
            T d = T(Z[r] + vv[k] * c[r]);
            if (num_traits<T>::to_double(d) < 2.2250738585072014e-308) d = tiny;
            D(k, r) = d;
            f += N[r] * log(d);
        }
        fv[k] = f;
        if (num_traits<T>::to_double(f) > mx) mx = num_traits<T>::to_double(f);
    }
    T mxT = num_traits<T>::from_double(mx);
    std::vector<T> e(nq);
    T se = zero;
    for (std::size_t k = 0; k < nq; ++k) {
        e[k] = T(W[k] * exp(T(fv[k] - mxT)));
        se += e[k];
    }
    Radial<T> out;
    out.lJ = T(mxT + log(se));
    out.G.assign(R, zero);
    out.vbar = zero;
    Matrix<T> Tm(nq, R);
    std::vector<T> p(nq);
    for (std::size_t k = 0; k < nq; ++k) {
        p[k] = T(e[k] / se);
        out.vbar += p[k] * vv[k];
        for (std::size_t r = 0; r < R; ++r) {
            Tm(k, r) = T(N[r] * vv[k] / D(k, r));
            out.G[r] += p[k] * Tm(k, r);
        }
    }
    Matrix<T> et2(R, R);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t s = 0; s < R; ++s) et2(r, s) = zero;
    for (std::size_t k = 0; k < nq; ++k)
        for (std::size_t r = 0; r < R; ++r) {
            T pt = T(p[k] * Tm(k, r));
            for (std::size_t s = 0; s < R; ++s) et2(r, s) += pt * Tm(k, s);
        }
    out.Lam = Matrix<T>(R, R);
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t s = 0; s < R; ++s)
            out.Lam(r, s) = T(et2(r, s) - out.G[r] * out.G[s]);
        if (num_traits<T>::to_double(N[r]) > 0.0)
            out.Lam(r, r) = T(out.Lam(r, r) - et2(r, r) / N[r]);
    }
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t s = r + 1; s < R; ++s) {
            T m = T(half * T(out.Lam(r, s) + out.Lam(s, r)));
            out.Lam(r, s) = m;
            out.Lam(s, r) = m;
        }
    return out;
}

/** Mode, curvature and log-integrand at the mode of the simplex factor. */
template <class T>
struct Mode {
    std::vector<T> x;
    Matrix<T> A;
    T ld;
    T h0;
};

/**
 * Mode and curvature of h(w) = log J(L'x(w)) + sum_i log x_i with J the exact radial
 * integral. The fixed point x = (1 + x.*(L*G))/vbar is the Z > 0 analogue of
 * pfqn_le_fpi: integrating by parts gives sum_i x_i (L*G)_i = vbar - M, so the update is
 * normalised by construction, and at Z = 0 it reduces to pfqn_le_fpi. The term in the
 * second derivative of x(w) drops at the mode against sum_i x_i == 1.
 */
template <class T>
Mode<T> simplex_mode(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    using std::abs;
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> x;
    T vstart = zero;
    pfqn_le_fpiZ(L, N, Z, x, vstart);  // logistic-expansion mode as a warm start
    std::vector<T> x1(M, num_traits<T>::from_double(1e300));
    for (int it = 0; it < 10000; ++it) {
        double diff = 0.0;
        for (std::size_t i = 0; i < M; ++i)
            diff += num_traits<T>::to_double(abs(T(x[i] - x1[i])));
        if (diff <= 1e-11) break;
        x1 = x;
        std::vector<T> c(R, zero);
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t i = 0; i < M; ++i) c[r] += x1[i] * L(i, r);
        Radial<T> rad = radial(c, N, Z, M);
        T s = zero;
        for (std::size_t i = 0; i < M; ++i) {
            T lg = zero;
            for (std::size_t r = 0; r < R; ++r) lg += L(i, r) * rad.G[r];
            x[i] = T(T(one + x1[i] * lg) / rad.vbar);
            s += x[i];
        }
        for (std::size_t i = 0; i < M; ++i) x[i] = T(x[i] / s);
    }
    std::vector<T> c(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) c[r] += x[i] * L(i, r);
    Radial<T> rad = radial(c, N, Z, M);
    // P = L*Lam*L' - diag(1/x^2); A = -Jm'*P*Jm with Jm = (diag(x)-x*x')(:,1:M-1).
    Matrix<T> P(M, M);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            T acc = zero;
            for (std::size_t r = 0; r < R; ++r) {
                T lr = zero;
                for (std::size_t s = 0; s < R; ++s) lr += rad.Lam(r, s) * L(j, s);
                acc += L(i, r) * lr;
            }
            P(i, j) = acc;
        }
        P(i, i) = T(P(i, i) - one / T(x[i] * x[i]));
    }
    const std::size_t d = M - 1;
    Matrix<T> Jm(M, d);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t aI = 0; aI < d; ++aI)
            Jm(i, aI) = T((i == aI ? x[i] : zero) - x[i] * x[aI]);
    Mode<T> out;
    out.A = Matrix<T>(d, d);
    for (std::size_t aI = 0; aI < d; ++aI)
        for (std::size_t bI = 0; bI < d; ++bI) {
            T acc = zero;
            for (std::size_t i = 0; i < M; ++i) {
                T pj = zero;
                for (std::size_t j = 0; j < M; ++j) pj += P(i, j) * Jm(j, bI);
                acc += Jm(i, aI) * pj;
            }
            out.A(aI, bI) = T(zero - acc);
        }
    const T half = num_traits<T>::from_double(0.5);
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = i + 1; j < d; ++j) {
            T m = T(half * T(out.A(i, j) + out.A(j, i)));
            out.A(i, j) = m;
            out.A(j, i) = m;
        }
    out.x = x;
    out.ld = (d == 0) ? zero : detail::pfqn_logdet(out.A);
    T sum_lx = zero;
    for (std::size_t i = 0; i < M; ++i) sum_lx += log(x[i]);
    out.h0 = T(rad.lJ + sum_lx);
    return out;
}

/** softmax of [w; 0], the logistic parametrisation of the simplex with gauge w_M = 0. */
template <class T>
std::vector<T> softmax_gauge(const std::vector<T>& w) {
    using std::exp;
    const std::size_t M = w.size() + 1;
    std::vector<T> a(M, num_traits<T>::from_int(0));
    double mx = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        a[i] = w[i];
        if (num_traits<T>::to_double(a[i]) > mx) mx = num_traits<T>::to_double(a[i]);
    }
    T mxT = num_traits<T>::from_double(mx);
    std::vector<T> x(M);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < M; ++i) {
        x[i] = exp(T(a[i] - mxT));
        s += x[i];
    }
    for (std::size_t i = 0; i < M; ++i) x[i] = T(x[i] / s);
    return x;
}

/**
 * Log of the tensor Gauss-Hermite sum, accumulated with a running maximum; the
 * det(A)^(-1/2) of the rule is applied by the caller. A tensor rule is NOT invariant to
 * the choice of A^(-1/2): the principal-axis frame is used, as in the reference results.
 */
template <class T, class F>
T aghq_rule(const F& h, const std::vector<T>& w0, const T& h0, const Matrix<T>& A,
            std::size_t q, std::size_t d) {
    using std::exp;
    using std::log;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    if (d == 0) return zero;
    double nodesd = std::pow(static_cast<double>(q), static_cast<double>(d));
    if (nodesd > 1e7)
        throw InputError("pfqn_aghq: the tensor rule needs more than 1e7 nodes; reduce q or use pfqn_le");
    const std::size_t nodes = static_cast<std::size_t>(nodesd);
    std::vector<T> lam;
    Matrix<T> V(d, d);
    sym_eig(A, lam, V);
    Matrix<T> B(d, d);
    for (std::size_t j = 0; j < d; ++j) {
        if (num_traits<T>::to_double(lam[j]) <= 0.0)
            throw InputError("pfqn_aghq: the curvature at the mode is not positive definite");
        T sc = T(num_traits<T>::from_int(1) / sqrt(lam[j]));
        for (std::size_t i = 0; i < d; ++i) B(i, j) = T(V(i, j) * sc);
    }
    std::vector<T> z, wt;
    gauss_hermite<T>(q, z, wt);
    std::vector<T> lwt(q);
    for (std::size_t k = 0; k < q; ++k) lwt[k] = log(wt[k]);
    std::vector<std::size_t> idx(d, 0);
    double lmax = -1e308;
    T s = zero;
    std::vector<T> zz(d), w(d);
    bool first = true;
    for (std::size_t k = 0; k < nodes; ++k) {
        T lw = zero, zsq = zero;
        for (std::size_t j = 0; j < d; ++j) {
            zz[j] = z[idx[j]];
            lw += lwt[idx[j]];
            zsq += zz[j] * zz[j];
        }
        for (std::size_t i = 0; i < d; ++i) {
            T acc = w0[i];
            for (std::size_t j = 0; j < d; ++j) acc += B(i, j) * zz[j];
            w[i] = acc;
        }
        T lt = T(lw + h(w) - h0 + num_traits<T>::from_double(0.5) * zsq);
        double ltd = num_traits<T>::to_double(lt);
        if (first || ltd > lmax) {
            T lmaxT = num_traits<T>::from_double(lmax);
            s = first ? num_traits<T>::from_int(1)
                      : T(s * exp(T(lmaxT - lt)) + num_traits<T>::from_int(1));
            lmax = ltd;
            first = false;
        } else {
            s += exp(T(lt - num_traits<T>::from_double(lmax)));
        }
        for (std::size_t j = d; j-- > 0;) {
            if (++idx[j] < q) break;
            idx[j] = 0;
        }
    }
    return T(num_traits<T>::from_double(lmax) + log(s));
}

}  // namespace simplex
}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_SIMPLEX_H
