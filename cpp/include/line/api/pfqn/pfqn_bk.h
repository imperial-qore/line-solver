/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BK_H
#define LINE_API_PFQN_PFQN_BK_H

/**
 * Birman-Kogan asymptotic evaluation of closed networks with many stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_bk.m, pfqn_bkue.m and
 * pfqn_bklc.m. Birman and Kogan (Communications in Statistics. Stochastic
 * Models 8(3):543-563, 1992) evaluate the multichain partition function by the
 * saddle point method applied to the Cauchy inversion of its generating
 * function. Three algorithms live here:
 *
 *   pfqn_bk      Propositions 1 and 3 with Algorithm 1. Stations that serve a
 *                single chain and appear only once (the paper's dedicated
 *                single servers) stay OUTSIDE the exponent as O(1) algebraic
 *                factors, so their poles may be crossed by the saddle point;
 *                Algorithm 1 detects those chains and pins their coordinate on
 *                the pole, where the residue rather than the saddle carries the
 *                mass. The rest are the paper's large groups of identical
 *                stations and are exponentiated.
 *   pfqn_bkue    The van der Waerden uniform expansion of Section 4, which
 *                keeps one dominant pole and the saddle in a single erfc
 *                formula and so stays accurate on both sides of the crossing.
 *   pfqn_bklc    Algorithm 2, the load concealment reduction of a
 *                multichain network to single chain problems.
 *
 * ARITHMETIC. Logarithms, an error function and a Newton iteration, so gated on
 * num_traits<T>::has_transcendental.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_bk, mirroring [G, lG, X, U, A, B]. */
template <class T>
struct BkResult {
    T G;
    T lG;
    std::vector<T> X;               ///< the saddle point coordinates
    Matrix<T> U;                    ///< (M x R) utilizations
    std::vector<std::size_t> A;     ///< chains whose dedicated station is not saturated
    std::vector<std::size_t> B;     ///< chains whose dedicated station is a bottleneck
};

/** Return value of pfqn_bklc. */
template <class T>
struct BkLcResult {
    std::vector<T> X;
    Matrix<T> Q;
    Matrix<T> U;
    int it;
};

namespace detail {

/** Number of stations sharing each demand row, up to relative rounding. */
template <class T>
std::vector<std::size_t> bk_multiplicity(const std::vector<std::vector<T>>& L) {
    const std::size_t M = L.size();
    std::vector<std::size_t> mult(M, 1);
    if (M == 0) return mult;
    const std::size_t R = L[0].size();
    const double tol = 1e-8;  // GlobalConstants.FineTol
    for (std::size_t i = 0; i < M; ++i) {
        if (mult[i] > 1) continue;
        for (std::size_t j = i + 1; j < M; ++j) {
            double scale = 1.0, diff = 0.0;
            for (std::size_t r = 0; r < R; ++r) {
                const double a = num_traits<T>::to_double(L[i][r]);
                const double b = num_traits<T>::to_double(L[j][r]);
                scale = std::max(scale, std::max(std::abs(a), std::abs(b)));
                diff = std::max(diff, std::abs(a - b));
            }
            if (diff <= tol * scale) {
                ++mult[i];
                ++mult[j];
            }
        }
    }
    return mult;
}

/** Exponent of the integrand, groups only (eq. 11 and 23 in unscaled variables). */
template <class T>
T bk_psi(const std::vector<T>& z, const std::vector<std::vector<T>>& Lg, const std::vector<T>& N,
         const std::vector<T>& Z) {
    using std::log;
    const T one = num_traits<T>::from_int(1);
    T f = num_traits<T>::from_int(0);
    for (std::size_t r = 0; r < z.size(); ++r) f += Z[r] * z[r] - N[r] * log(z[r]);
    for (std::size_t i = 0; i < Lg.size(); ++i) {
        T u = num_traits<T>::from_int(0);
        for (std::size_t r = 0; r < z.size(); ++r) u += Lg[i][r] * z[r];
        f -= log(one - u);
    }
    return f;
}

template <class T>
std::vector<T> bk_grad(const std::vector<T>& z, const std::vector<std::vector<T>>& Lg,
                       const std::vector<T>& N, const std::vector<T>& Z) {
    const T one = num_traits<T>::from_int(1);
    const std::size_t R = z.size();
    std::vector<T> g(R);
    for (std::size_t r = 0; r < R; ++r) g[r] = Z[r] - N[r] / z[r];
    for (std::size_t i = 0; i < Lg.size(); ++i) {
        T u = num_traits<T>::from_int(0);
        for (std::size_t r = 0; r < R; ++r) u += Lg[i][r] * z[r];
        const T d = one / (one - u);
        for (std::size_t r = 0; r < R; ++r) g[r] += d * Lg[i][r];
    }
    return g;
}

template <class T>
Matrix<T> bk_hessian(const std::vector<T>& z, const std::vector<std::vector<T>>& Lg,
                     const std::vector<T>& N) {
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    const std::size_t R = z.size();
    Matrix<T> H(R, R, zero);
    for (std::size_t r = 0; r < R; ++r) H(r, r) = N[r] / (z[r] * z[r]);
    for (std::size_t i = 0; i < Lg.size(); ++i) {
        T u = zero;
        for (std::size_t r = 0; r < R; ++r) u += Lg[i][r] * z[r];
        const T d = one / (one - u);
        const T d2 = d * d;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) H(r, s) += d2 * Lg[i][r] * Lg[i][s];
    }
    return H;
}

}  // namespace detail

/**
 * Birman-Kogan saddle point normalizing constant with bottleneck detection.
 *
 * @param L (M x R) service demands
 * @param N (R) population
 * @param Z (R) think times, may be empty
 */
template <class T>
BkResult<T> pfqn_bk(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_bk requires transcendental arithmetic (saddle point expansion of log G)");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T fineTol = num_traits<T>::from_double(1e-8);

    BkResult<T> res;
    const std::size_t M = L.rows(), R = L.cols();
    res.G = one;
    res.lG = zero;
    res.X.assign(R, zero);
    res.U = Matrix<T>(M, R, zero);
    T Ntot = zero;
    for (std::size_t r = 0; r < N.size(); ++r) Ntot += N[r];
    if (L.empty() || N.empty() || !(Ntot > zero)) {
        for (std::size_t r = 0; r < R; ++r) res.A.push_back(r);
        return res;
    }
    if (N.size() != R) throw InputError("pfqn_bk: L and N disagree on the class count");
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);
    if (Zv.size() != R) throw InputError("pfqn_bk: L and Z disagree on the class count");

    // An empty class contributes a factor of 1 and has no saddle coordinate.
    std::size_t nkeep = 0;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > zero) ++nkeep;
    if (nkeep > 0 && nkeep < R) {
        Matrix<T> Lk(M, nkeep, zero);
        std::vector<T> Nk, Zk;
        std::vector<std::size_t> map;
        std::size_t c = 0;
        for (std::size_t r = 0; r < R; ++r) {
            if (!(N[r] > zero)) continue;
            for (std::size_t i = 0; i < M; ++i) Lk(i, c) = L(i, r);
            Nk.push_back(N[r]);
            Zk.push_back(Zv[r]);
            map.push_back(r);
            ++c;
        }
        BkResult<T> red = pfqn_bk(Lk, Nk, Zk);
        res.G = red.G;
        res.lG = red.lG;
        for (std::size_t k = 0; k < map.size(); ++k) {
            res.X[map[k]] = red.X[k];
            for (std::size_t i = 0; i < M; ++i) res.U(i, map[k]) = red.U(i, k);
        }
        for (std::size_t k = 0; k < red.A.size(); ++k) res.A.push_back(map[red.A[k]]);
        for (std::size_t k = 0; k < red.B.size(); ++k) res.B.push_back(map[red.B[k]]);
        return res;
    }

    // stations with no demand at all do not enter the generating function
    std::vector<std::vector<T>> Lq;
    for (std::size_t i = 0; i < M; ++i) {
        T s = zero;
        for (std::size_t r = 0; r < R; ++r) s += L(i, r);
        if (!(s > zero)) continue;
        std::vector<T> row(R);
        for (std::size_t r = 0; r < R; ++r) row[r] = L(i, r);
        Lq.push_back(row);
    }

    // Dedicated station of each chain: single chain, no identical twin, no think
    // time, and only when the model holds a group of identical stations, since it
    // is against M_j >> 1 replicas that a lone station is an O(1) factor.
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> mu(R, inf);
    std::vector<long> poleRow(R, -1);
    std::vector<bool> isPole(Lq.size(), false);
    std::vector<std::size_t> mult = detail::bk_multiplicity(Lq);
    bool hasGroup = false;
    for (std::size_t i = 0; i < mult.size(); ++i)
        if (mult[i] > 1) hasGroup = true;
    if (hasGroup) {
        for (std::size_t i = 0; i < Lq.size(); ++i) {
            if (mult[i] > 1) continue;
            std::size_t nz = 0, cnt = 0;
            for (std::size_t r = 0; r < R; ++r)
                if (Lq[i][r] > zero) {
                    nz = r;
                    ++cnt;
                }
            if (cnt != 1) continue;
            if (Zv[nz] > fineTol) continue;
            const double cand = 1.0 / num_traits<T>::to_double(Lq[i][nz]);
            if (cand < mu[nz]) {
                mu[nz] = cand;
                poleRow[nz] = static_cast<long>(i);
            }
        }
        for (std::size_t r = 0; r < R; ++r)
            if (poleRow[r] >= 0) isPole[static_cast<std::size_t>(poleRow[r])] = true;
    }
    std::vector<std::vector<T>> Lg;
    for (std::size_t i = 0; i < Lq.size(); ++i)
        if (!isPole[i]) Lg.push_back(Lq[i]);

    // Algorithm 1 as an active set method on the strictly convex psi.
    std::vector<T> muT(R);
    for (std::size_t r = 0; r < R; ++r)
        muT[r] = std::isinf(mu[r]) ? num_traits<T>::from_double(0.0)
                                   : num_traits<T>::from_double(mu[r]);
    std::vector<bool> onBound(R, false);
    std::vector<T> z(R);
    for (std::size_t r = 0; r < R; ++r) {
        T den = Zv[r];
        for (std::size_t i = 0; i < Lg.size(); ++i) den += Lg[i][r];
        if (!(den > zero)) den = fineTol;
        z[r] = N[r] / den;
        if (!std::isinf(mu[r])) {
            const T cap = num_traits<T>::from_double(0.99 * mu[r]);
            if (z[r] > cap) z[r] = cap;
        }
    }
    for (int it = 0; it < 200 && !Lg.empty(); ++it) {
        T umax = zero;
        for (std::size_t i = 0; i < Lg.size(); ++i) {
            T u = zero;
            for (std::size_t r = 0; r < R; ++r) u += Lg[i][r] * z[r];
            if (u > umax) umax = u;
        }
        if (num_traits<T>::to_double(umax) < 0.9) break;
        for (std::size_t r = 0; r < R; ++r) z[r] = z[r] * num_traits<T>::from_double(0.7);
    }
    for (std::size_t outer = 0; outer <= R; ++outer) {
        for (std::size_t r = 0; r < R; ++r)
            if (onBound[r]) z[r] = muT[r];
        std::vector<std::size_t> freeIdx;
        for (std::size_t r = 0; r < R; ++r)
            if (!onBound[r]) freeIdx.push_back(r);
        if (freeIdx.empty()) break;
        const std::size_t nf = freeIdx.size();
        for (int it = 0; it < 500; ++it) {
            std::vector<T> g = detail::bk_grad(z, Lg, N, Zv);
            double gn = 0.0;
            for (std::size_t i = 0; i < nf; ++i) {
                const double v = num_traits<T>::to_double(g[freeIdx[i]]);
                gn += v * v;
            }
            gn = std::sqrt(gn);
            if (gn <= 1e-12 * std::max(1.0, num_traits<T>::to_double(Ntot))) break;
            Matrix<T> H = detail::bk_hessian(z, Lg, N);
            Matrix<T> Hf(nf, nf, zero);
            std::vector<T> rhs(nf);
            for (std::size_t i = 0; i < nf; ++i) {
                for (std::size_t j = 0; j < nf; ++j) Hf(i, j) = H(freeIdx[i], freeIdx[j]);
                rhs[i] = zero - g[freeIdx[i]];
            }
            std::vector<T> dz;
            try {
                dz = solve(Hf, rhs);
            } catch (const std::exception&) {
                break;
            }
            double alpha = 1.0;
            bool ok = false;
            std::vector<T> zt(R);
            while (alpha >= 1e-14) {
                zt = z;
                for (std::size_t i = 0; i < nf; ++i)
                    zt[freeIdx[i]] = z[freeIdx[i]] + num_traits<T>::from_double(alpha) * dz[i];
                ok = true;
                for (std::size_t i = 0; i < nf && ok; ++i) {
                    const std::size_t r = freeIdx[i];
                    if (!(zt[r] > zero)) ok = false;
                    if (ok && !std::isinf(mu[r]) && zt[r] > muT[r]) ok = false;
                }
                for (std::size_t i = 0; i < Lg.size() && ok; ++i) {
                    T u = zero;
                    for (std::size_t r = 0; r < R; ++r) u += Lg[i][r] * zt[r];
                    if (!(u < one)) ok = false;
                }
                if (ok) break;
                alpha /= 2;
            }
            if (!ok) break;
            z = zt;
        }
        // a chain whose descent direction still pushes past its pole is in B
        std::vector<T> g = detail::bk_grad(z, Lg, N, Zv);
        bool any = false;
        for (std::size_t i = 0; i < nf; ++i) {
            const std::size_t r = freeIdx[i];
            if (std::isinf(mu[r])) continue;
            if (num_traits<T>::to_double(z[r]) >= mu[r] * (1 - 1e-9) &&
                num_traits<T>::to_double(g[r]) < 0) {
                onBound[r] = true;
                any = true;
            }
        }
        if (!any) break;
    }
    for (std::size_t r = 0; r < R; ++r)
        if (onBound[r]) z[r] = muT[r];

    res.X = z;
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t i = 0; i < M; ++i) {
            T u = L(i, r) * z[r];
            if (onBound[r] && u > one) u = one;
            res.U(i, r) = u;
        }
        if (onBound[r])
            res.B.push_back(r);
        else
            res.A.push_back(r);
    }

    const T psi0 = detail::bk_psi(z, Lg, N, Zv);
    T lG;
    if (res.A.empty()) {  // eq. (25): the residues carry everything
        lG = psi0;
    } else {
        Matrix<T> H = detail::bk_hessian(z, Lg, N);
        const std::size_t nf = res.A.size();
        Matrix<T> Haa(nf, nf, zero);
        for (std::size_t i = 0; i < nf; ++i)
            for (std::size_t j = 0; j < nf; ++j) Haa(i, j) = H(res.A[i], res.A[j]);
        Matrix<T> LU = Haa;
        lu_factor(LU);
        T logdet = zero;
        for (std::size_t i = 0; i < nf; ++i) logdet += log(num_abs(LU(i, i)));
        lG = psi0 - num_traits<T>::from_double(0.5 * static_cast<double>(nf) * std::log(2 * M_PI)) -
             num_traits<T>::from_double(0.5) * logdet;
        for (std::size_t i = 0; i < nf; ++i) {
            const std::size_t r = res.A[i];
            lG -= log(z[r]);
            if (!std::isinf(mu[r])) lG -= log(one - z[r] / muT[r]);
        }
    }
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

/**
 * Scaled complementary error function exp(x^2)*erfc(x) for x >= 0. The direct
 * product overflows past x ~ 26, where the asymptotic series is already exact to
 * double precision.
 */
inline double bk_erfcx(double x) {
    if (x < 25.0) return std::exp(x * x) * std::erfc(x);
    const double y = 1.0 / (2.0 * x * x);
    double term = 1.0, sum = 1.0;
    for (int k = 1; k <= 12; ++k) {
        term *= -(2 * k - 1) * y;
        sum += term;
    }
    return sum / (x * std::sqrt(M_PI));
}

namespace detail {

template <class T>
T bk_h1(const T& z, const std::vector<T>& D, const T& N, const T& Z) {
    using std::log;
    const T one = num_traits<T>::from_int(1);
    T f = Z * z - N * log(z);
    for (std::size_t i = 0; i < D.size(); ++i) f -= log(one - D[i] * z);
    return f;
}

template <class T>
T bk_h1d1(const T& z, const std::vector<T>& D, const T& N, const T& Z) {
    const T one = num_traits<T>::from_int(1);
    T g = Z - N / z;
    for (std::size_t i = 0; i < D.size(); ++i) g += D[i] / (one - D[i] * z);
    return g;
}

template <class T>
T bk_h1d2(const T& z, const std::vector<T>& D, const T& N) {
    const T one = num_traits<T>::from_int(1);
    T h = N / (z * z);
    for (std::size_t i = 0; i < D.size(); ++i) {
        const T d = one - D[i] * z;
        h += D[i] * D[i] / (d * d);
    }
    return h;
}

template <class T>
T bk_h1d3(const T& z, const std::vector<T>& D, const T& N) {
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    T h = num_traits<T>::from_int(-2) * N / (z * z * z);
    for (std::size_t i = 0; i < D.size(); ++i) {
        const T d = one - D[i] * z;
        h += two * D[i] * D[i] * D[i] / (d * d * d);
    }
    return h;
}

template <class T>
T bk_saddle1(const std::vector<T>& D, const T& N, const T& Z) {
    if (D.empty()) return N / Z;
    double dmax = 0.0;
    for (std::size_t i = 0; i < D.size(); ++i)
        dmax = std::max(dmax, num_traits<T>::to_double(D[i]));
    const double hi = 1.0 / dmax;
    T z = num_traits<T>::from_double(0.5 * hi);
    for (int it = 0; it < 200; ++it) {
        const T g = bk_h1d1(z, D, N, Z);
        if (std::abs(num_traits<T>::to_double(g)) <=
            1e-14 * std::max(1.0, num_traits<T>::to_double(N)))
            break;
        const T dz = (num_traits<T>::from_int(0) - g) / bk_h1d2(z, D, N);
        double alpha = 1.0;
        while (true) {
            const double zt = num_traits<T>::to_double(z) + alpha * num_traits<T>::to_double(dz);
            if (zt > 0 && zt < hi) break;
            alpha /= 2;
            if (alpha < 1e-14) break;
        }
        if (alpha < 1e-14) break;
        z += num_traits<T>::from_double(alpha) * dz;
    }
    return z;
}

}  // namespace detail

/**
 * Birman-Kogan uniform (van der Waerden) expansion for a single chain.
 *
 * @param L (M) service demands, single class
 * @param N population
 * @param Z think time
 */
template <class T>
BkResult<T> pfqn_bkue(const std::vector<T>& L, const T& N, const T& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_bkue requires transcendental arithmetic (uniform expansion of log G)");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    BkResult<T> res;
    res.G = one;
    res.lG = zero;
    if (!(N > zero)) return res;
    std::vector<T> Lv;
    for (std::size_t i = 0; i < L.size(); ++i)
        if (L[i] > zero) Lv.push_back(L[i]);
    if (Lv.empty()) {
        res.lG = N * log(Z) - num_traits<T>::from_double(std::lgamma(num_traits<T>::to_double(N) + 1.0));
        res.G = exp(res.lG);
        return res;
    }
    std::size_t ipole = 0;
    for (std::size_t i = 1; i < Lv.size(); ++i)
        if (Lv[i] > Lv[ipole]) ipole = i;
    const double dmax = num_traits<T>::to_double(Lv[ipole]);
    const double tolL = 1e-8 * std::max(1.0, dmax);
    std::size_t ties = 0;
    std::vector<double> sorted;
    for (std::size_t i = 0; i < Lv.size(); ++i) {
        const double v = num_traits<T>::to_double(Lv[i]);
        if (std::abs(v - dmax) <= tolL) ++ties;
        sorted.push_back(v);
    }
    std::sort(sorted.begin(), sorted.end());
    bool hasGroup = false;
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i] - sorted[i - 1] <= tolL) hasGroup = true;
    const bool hasPole = (ties == 1) && hasGroup;
    std::vector<T> D;
    T zp = zero;
    if (hasPole) {
        for (std::size_t i = 0; i < Lv.size(); ++i)
            if (i != ipole) D.push_back(Lv[i]);
        zp = one / Lv[ipole];
    } else {
        D = Lv;
    }
    const T z0 = detail::bk_saddle1(D, N, Z);
    const T h2 = detail::bk_h1d2(z0, D, N);
    const T h3 = detail::bk_h1d3(z0, D, N);
    if (!hasPole) {
        // No pole to keep out of the exponent: the expansion degenerates to the
        // plain saddle point, and the third derivative term goes with the pole it
        // corrects.
        res.lG = detail::bk_h1(z0, D, N, Z) - log(z0) -
                 num_traits<T>::from_double(0.5) * log(num_traits<T>::from_double(2 * M_PI) * h2);
        res.G = exp(res.lG);
        return res;
    }
    const T t2 = (one / z0 + h3 / (num_traits<T>::from_int(6) * h2)) /
                 num_traits<T>::from_double(
                     std::sqrt(2 * M_PI * num_traits<T>::to_double(h2)));
    const double b2 = std::max(0.0, num_traits<T>::to_double(detail::bk_h1(zp, D, N, Z)) -
                                        num_traits<T>::to_double(detail::bk_h1(z0, D, N, Z)));
    if (num_traits<T>::to_double(zp) >= num_traits<T>::to_double(z0)) {  // saddle before the pole
        res.lG = detail::bk_h1(z0, D, N, Z) +
                 log(num_traits<T>::from_double(0.5 * bk_erfcx(std::sqrt(b2))) + t2);
    } else {  // the pole has been crossed and its residue leads
        res.lG = detail::bk_h1(zp, D, N, Z) +
                 log(num_traits<T>::from_double(1.0 - 0.5 * std::erfc(std::sqrt(b2))) +
                     t2 * num_traits<T>::from_double(std::exp(-b2)));
    }
    res.G = exp(res.lG);
    return res;
}

/**
 * Birman-Kogan load concealment algorithm (Algorithm 2).
 *
 * @param L (M x R) service demands
 * @param N (R) population
 * @param Z (R) think times, may be empty
 * @param method single chain solver, "mva" (default) or "ue"
 * @param tol convergence tolerance on the throughputs
 * @param maxiter maximum number of sweeps
 */
template <class T>
BkLcResult<T> pfqn_bklc(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                            const std::string& method = "mva", double tol = 1e-10,
                            int maxiter = 1000) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_bklc requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.rows(), R = L.cols();
    BkLcResult<T> res;
    res.X.assign(R, zero);
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.it = 0;
    T Ntot = zero;
    for (std::size_t r = 0; r < N.size(); ++r) Ntot += N[r];
    if (L.empty() || !(Ntot > zero)) return res;
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);
    const bool ue = (method == "ue");
    if (tol <= 0) tol = 1e-10;
    if (maxiter <= 0) maxiter = 1000;

    // Step 1: the saddle point utilizations of Corollary 1 seed the iteration
    std::vector<T> X(R, zero);
    try {
        BkResult<T> seed = pfqn_bk(L, N, Zv);
        for (std::size_t r = 0; r < R; ++r) {
            const double v = num_traits<T>::to_double(seed.X[r]);
            X[r] = (std::isfinite(v) && v >= 0) ? seed.X[r] : zero;
        }
    } catch (const std::exception&) {
    }
    for (std::size_t r = 0; r < R; ++r) {
        T cap = zero, sum = zero;
        for (std::size_t i = 0; i < M; ++i) {
            if (L(i, r) > cap) cap = L(i, r);
            sum += L(i, r);
        }
        if (!(X[r] > zero) && N[r] > zero) X[r] = N[r] / (Zv[r] + sum);
        if (cap > zero && X[r] > one / cap) X[r] = one / cap;
    }

    Matrix<T> Q(M, R, zero);
    for (res.it = 1; res.it <= maxiter; ++res.it) {
        std::vector<T> Xold = X;
        for (std::size_t l = 0; l < R; ++l) {
            if (!(N[l] > zero)) {
                X[l] = zero;
                for (std::size_t i = 0; i < M; ++i) Q(i, l) = zero;
                continue;
            }
            // Step 2a: residual capacity left to chain l at every station
            std::vector<T> D(M);
            for (std::size_t i = 0; i < M; ++i) {
                T busy = zero;
                for (std::size_t k = 0; k < R; ++k)
                    if (k != l) busy += L(i, k) * X[k];
                T A = one - busy;
                if (num_traits<T>::to_double(A) < 1e-8) A = num_traits<T>::from_double(1e-8);
                D[i] = L(i, l) / A;
            }
            // Step 2b: solve the single chain network with the thinned rates
            if (ue) {
                std::vector<T> Qi(M, zero);
                T lgPrev = zero;
                T Xl = zero;
                const int Nl = static_cast<int>(std::llround(num_traits<T>::to_double(N[l])));
                for (int nn = 1; nn <= Nl; ++nn) {
                    const T lgn = pfqn_bkue(D, num_traits<T>::from_int(nn), Zv[l]).lG;
                    Xl = num_traits<T>::from_double(
                        std::exp(num_traits<T>::to_double(lgPrev) - num_traits<T>::to_double(lgn)));
                    for (std::size_t i = 0; i < M; ++i) Qi[i] = D[i] * Xl * (one + Qi[i]);
                    lgPrev = lgn;
                }
                X[l] = Xl;
                for (std::size_t i = 0; i < M; ++i) Q(i, l) = Qi[i];
            } else {
                Matrix<T> Dm(M, 1, zero);
                for (std::size_t i = 0; i < M; ++i) Dm(i, 0) = D[i];
                std::vector<int> Nl(1, static_cast<int>(std::llround(num_traits<T>::to_double(N[l]))));
                Matrix<T> Zl(1, 1, Zv[l]);
                MvaResult<T> mva = pfqn_mva(Dm, Nl, Zl);
                X[l] = mva.XN[0];
                for (std::size_t i = 0; i < M; ++i) Q(i, l) = mva.QN(i, 0);
            }
        }
        double diff = 0.0, xmax = 1.0;
        for (std::size_t r = 0; r < R; ++r) {
            diff = std::max(diff, std::abs(num_traits<T>::to_double(X[r]) -
                                           num_traits<T>::to_double(Xold[r])));
            xmax = std::max(xmax, std::abs(num_traits<T>::to_double(X[r])));
        }
        if (diff <= tol * xmax) break;
    }
    if (res.it > maxiter) res.it = maxiter;
    res.X = X;
    res.Q = Q;
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) res.U(i, r) = L(i, r) * X[r];
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BK_H
