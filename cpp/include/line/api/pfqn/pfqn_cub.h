/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_CUB_H
#define LINE_API_PFQN_PFQN_CUB_H

/**
 * Normalizing constant by Grundmann-Moeller cubature over the simplex.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_cub.m together with its two local
 * functions simplexquad and grnmol (Grundmann and Moller, SIAM J. Numer. Anal.
 * 15 (1978) 282-290). With Z = 0 the constant is a single integral of
 * prod_r (u' L(:,r))^{N_r} over the (M-1)-simplex; the degree-(2s+1) rule is
 * EXACT once s >= ceil((sum N - 1)/2), because the integrand is then a
 * polynomial of degree sum(N) <= 2s+1. With Z > 0 an outer integral over the
 * McKenna-Mitra scale variable v is added on a uniform grid, which is where
 * the method stops being exact.
 *
 * The `grnmol` rule is exported, because it is the only working
 * Grundmann-Moeller implementation in the reference tree; see the note on
 * pfqn_grnmol in the report accompanying this port.
 *
 * ARITHMETIC. The rule itself is a weighted sum of integrand values at
 * rational barycentric points, so at Z = 0 and full degree the whole
 * computation would be exact in a field -- were it not for the
 * exp(gammaln(1+sum N+M-1) - sum gammaln(1+N)) prefactor, which MATLAB forms in
 * logarithms. Since that prefactor is a ratio of factorials it could be formed
 * exactly, but the reference does not, and reproducing the reference is the
 * contract; the routine is therefore gated on
 * num_traits<T>::has_transcendental, and the exactness claim above is about
 * the cubature, not about the returned scalar.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_cub, mirroring [Gn, lGn]. */
template <class T>
struct CubResult {
    T G;
    T lG;
};

/**
 * Grundmann-Moeller rule of degrees 1, 3, ..., 2s+1 over the n-simplex with
 * vertices the columns of the identity (MATLAB's grnmol on V = eye(n,n+1)).
 *
 * @param f   integrand, evaluated on the n free barycentric coordinates
 * @param n   simplex dimension
 * @param s   maximum rule order
 * @param tol relative stopping tolerance between consecutive degrees
 * @return    the successive estimates, the last of which is the answer
 */
template <class T>
std::vector<T> grnmol(const std::function<T(const std::vector<T>&)>& f, std::size_t n, int s,
                      const T& tol) {
    // exactness-in-any-field rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    if (n == 0) throw InputError("grnmol: zero-dimensional simplex");
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> Q, Qv;
    const T Vol = T(num_traits<T>::from_int(1) / num_factorial<T>(static_cast<unsigned>(n)));
    int d = 0;
    while (true) {
        const long m = static_cast<long>(n) + 2 * d + 1;
        std::vector<long> al(n, 1);
        long alz = 2 * d + 1;
        T Qs = zero;
        while (true) {
            // barycentric evaluation point rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
            std::vector<T> x(n);
            x[0] = T(num_traits<T>::from_int(alz) / num_traits<T>::from_int(m));
            for (std::size_t j = 1; j < n; ++j)
                x[j] = T(num_traits<T>::from_int(al[j - 1]) / num_traits<T>::from_int(m));
            Qs += f(x);
            for (std::size_t j = 0; j < n; ++j) {
                alz -= 2;
                if (alz > 0) {
                    al[j] += 2;
                    break;
                }
                alz += al[j] + 1;
                al[j] = 1;
            }
            if (alz == 2 * d + 1) break;
        }
        ++d;
        Qv.push_back(T(Vol * Qs));
        T q = zero;
        T p = num_traits<T>::from_int(2);
        for (long k = static_cast<long>(n) + 1; k <= m; ++k) p /= num_traits<T>::from_int(2 * k);
        for (int i = 1; i <= d; ++i) {
            q += num_pow_int(num_traits<T>::from_int(m + 2 - 2 * i), static_cast<unsigned>(2 * d - 1)) *
                 p * Qv[static_cast<std::size_t>(d - i)];
            p = T(-p * num_traits<T>::from_int(m + 1 - i) / num_traits<T>::from_int(i));
        }
        Q.push_back(q);
        // MATLAB's test is abs(Q(d)-Q(d-1)) < tol*Q(d-1), with Q(d-1) SIGNED,
        // so a negative previous estimate never stops the loop. Kept as is.
        if (d > s || (d > 1 && num_abs(T(Q[static_cast<std::size_t>(d) - 1] -
                                         Q[static_cast<std::size_t>(d) - 2])) <
                                   T(tol * Q[static_cast<std::size_t>(d) - 2])))
            break;
    }
    return Q;
}

/**
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param Z     (R) think times, empty or all zero for the exact branch
 * @param order rule degree; the default ceil((sum N - 1)/2) makes the Z = 0
 *              branch exact
 * @param atol  absolute tolerance, also the zero test on sum(Z)
 */
template <class T>
CubResult<T> pfqn_cub(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                      int order, const T& atol) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_cub requires transcendental arithmetic (the factorial prefactor is formed "
                  "in logarithms, and the Z > 0 branch is a quadrature)");
    using std::exp;
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_cub: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    CubResult<T> res;

    long Nt = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_cub: negative population");
        Nt += v;
    }
    if (M == 0 || N.empty() || Nt == 0) {
        res.G = one;
        res.lG = zero;
        return res;
    }
    if (M == 1) throw InputError("pfqn_cub: the simplex is degenerate for a single station");

    T Zsum = zero;
    for (const T& v : Z) Zsum += v;

    if (Z.empty() || Zsum < atol) {
        // Integrand prod_r (u' L(:,r))^{N_r} on the (M-1)-simplex.
        const std::function<T(const std::vector<T>&)> f = [&](const std::vector<T>& x) {
            T last = one;
            for (const T& v : x) last -= v;
            T prod = one;
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) continue;
                T uL = zero;
                for (std::size_t i = 0; i + 1 < M; ++i) uL += x[i] * L(i, r);
                uL += last * L(M - 1, r);
                prod *= num_pow_int(uL, static_cast<unsigned>(N[r]));
            }
            return prod;
        };
        const std::vector<T> Q = grnmol<T>(f, M - 1, order, atol);
        T coeff = detail::num_factln<T>(num_traits<T>::from_int(Nt + static_cast<long>(M) - 1));
        for (std::size_t r = 0; r < R; ++r)
            coeff -= detail::num_factln<T>(num_traits<T>::from_int(N[r]));
        res.G = T(Q.back() * exp(coeff));
        res.lG = log(res.G);
        return res;
    }

    // Z > 0: outer McKenna-Mitra integral on a uniform grid of 1e4 steps.
    const long steps = 10000;
    const T vmax = num_traits<T>::from_int(10 * Nt);
    const T dv = T(vmax / num_traits<T>::from_int(steps));
    T Gn = zero;
    for (long k = 0; k <= steps; ++k) {
        const T v = T(dv * num_traits<T>::from_int(k));
        Matrix<T> Lv(M, R);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) Lv(i, r) = T(L(i, r) * v + Z[r]);
        const std::function<T(const std::vector<T>&)> f = [&](const std::vector<T>& x) {
            T last = one;
            for (const T& val : x) last -= val;
            T s = zero;
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) continue;
                T uL = zero;
                for (std::size_t i = 0; i + 1 < M; ++i) uL += x[i] * Lv(i, r);
                uL += last * Lv(M - 1, r);
                s += num_traits<T>::from_int(N[r]) * log(uL);
            }
            return T(exp(s));
        };
        const std::vector<T> Q = grnmol<T>(f, M - 1, order, atol);
        const T dG =
            T(exp(T(-v)) * num_pow_int(v, static_cast<unsigned>(M - 1)) * Q.back() * dv);
        Gn += dG;
        if (k > 0 && Gn > zero && T(dG / Gn) < atol) break;
    }
    T coeff = zero;
    for (std::size_t r = 0; r < R; ++r) coeff -= detail::num_factln<T>(num_traits<T>::from_int(N[r]));
    res.G = T(Gn * exp(coeff));
    res.lG = log(res.G);
    return res;
}

template <class T>
CubResult<T> pfqn_cub(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    long Nt = 0;
    for (int v : N) Nt += v;
    const int order = static_cast<int>((Nt - 1 + 1) / 2);  // ceil((Nt-1)/2)
    return pfqn_cub(L, N, Z, order, num_traits<T>::from_double(1e-8));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_CUB_H
