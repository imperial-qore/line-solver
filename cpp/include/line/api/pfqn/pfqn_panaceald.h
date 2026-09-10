/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_PANACEALD_H
#define LINE_API_PFQN_PFQN_PANACEALD_H

/**
 * PANACEA normal-usage asymptotic expansion for LOAD-DEPENDENT closed networks
 * (Mitra and McKenna, JACM 33(3):568-592, 1986).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_panaceald.m. The expansion
 * coefficients A_n are linear combinations of partition functions of a
 * PSEUDONETWORK whose load dependence is the phi(n) transform of the original
 * rate lattice {mu_i(n)}; the transform itself is carried by logpsi below, the
 * partition functions by a direct load-dependent convolution over a lattice
 * that never exceeds 2*(terms-1) jobs per class.
 *
 * WHERE THE INFINITE SERVERS GO. A type-3 (infinite-server) row is ABSENT from
 * the pseudonetwork and enters only through the expansion parameter rho_j0.
 * solver_ncld encodes such a row as mu(i,n) = n, so the rows whose rate lattice
 * is exactly 1,2,...,Nt are detected here and folded into the think time, which
 * is why the caller may pass the delay either in Z or as such a row of L.
 *
 * WHEN IT DOES NOT APPLY, and this is the part worth knowing. The expansion is
 * an asymptotic series in the population and converges only in NORMAL USAGE,
 * i.e. alpha_i = 1 - lambda_i / mu_i(Ntot) > 0 at every queueing centre. MATLAB
 * returns NaN outside it and its caller turns that into an error. This port
 * reports it through `normalUsage` plus a `reason`, as pfqn_panacea does, so a
 * value that must not be used cannot be mistaken for one that may -- and so the
 * caller can say WHICH of the four conditions declined.
 *
 * ACCURACY REGIME (measured, 2026-07-24). Relative error on lG falls with the
 * population: 6.7e-4 at N = [10 10], 8.6e-8 at N = [400 400] at a fixed load
 * margin alpha_min ~ 0.70, while the runtime stays flat because the
 * pseudonetworks hold at most 4 jobs. At small N and heavy load it is dominated
 * by pfqn_clw_lld, which is exact there and cheap. Do not calibrate it only on
 * models small enough to have an exact reference.
 *
 * ARITHMETIC. A truncated asymptotic series reported as a logarithm, so gated
 * on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_panaceald, mirroring [Gn, lGn] plus why it declined. */
template <class T>
struct PanaceaLdResult {
    T G;
    T lG;
    bool normalUsage;    ///< false wherever the reference returns NaN
    const char* reason;  ///< which condition declined; nullptr when it applies
};

namespace detail {

/** MATLAB xlogy(e,x): e*log(x), with the convention 0*log(0) = 0. */
template <class T>
T pald_xlogy(const T& e, const T& x) {
    using std::log;
    if (e == num_traits<T>::from_int(0)) return num_traits<T>::from_int(0);
    return T(e * log(x));
}

/**
 * log psi(n) = log sum_{s>=n} [s!/(s-n)!] lambda^{s-n} / prod_{k<=s} mu(k),
 * the mu-free part of the phi(n) transform of eq. (3.7)-(3.8a).
 *
 * The series is split into the exact head s <= K and a geometric tail summed in
 * closed form through the Vandermonde identity; every term is positive, so the
 * whole sum is taken by logsumexp with no cancellation.
 *
 * @param lPirow log prod_{k=1}^{s} mu(k) for s = 0..K, i.e. K+1 entries
 */
template <class T>
T pald_logpsi(std::size_t n, const T& lambda, const std::vector<T>& lPirow, const T& muK,
              const T& alpha, std::size_t K) {
    using std::log;
    std::vector<T> t;
    t.reserve((K >= n ? K - n + 1 : 0) + n + 1);
    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    for (std::size_t s = n; s <= K; ++s) {
        const T sT = num_traits<T>::from_int(static_cast<long>(s));
        t.push_back(T(num_factln<T>(sT) - num_factln<T>(T(sT - nT)) +
                      pald_xlogy(T(sT - nT), lambda) - lPirow[s]));
    }
    const std::size_t Tm = n > K + 1 ? n : K + 1;
    const T TmT = num_traits<T>::from_int(static_cast<long>(Tm));
    for (std::size_t i = 0; i <= n; ++i) {
        const T iT = num_traits<T>::from_int(static_cast<long>(i));
        t.push_back(T(num_factln<T>(nT) + num_factln<T>(TmT) - num_factln<T>(T(nT - iT)) -
                      num_factln<T>(T(TmT - nT + iT)) + pald_xlogy(T(TmT + iT - nT), lambda) +
                      T(num_traits<T>::from_int(static_cast<long>(K)) - TmT - iT) * log(muK) -
                      T(iT + num_traits<T>::from_int(1)) * log(alpha) - lPirow[K]));
    }
    return logsumexp(t);
}

/** Mixed-radix index (0-based) to population vector, MATLAB idx2vec. */
inline void pald_idx2vec(std::size_t idx, const std::vector<int>& sizes, std::vector<int>& v) {
    std::size_t t = idx;
    for (std::size_t r = 0; r < sizes.size(); ++r) {
        v[r] = static_cast<int>(t % static_cast<std::size_t>(sizes[r]));
        t /= static_cast<std::size_t>(sizes[r]);
    }
}

/** Population vector to mixed-radix index (0-based), MATLAB vec2idx. */
inline std::size_t pald_vec2idx(const std::vector<int>& v, const std::vector<int>& sizes) {
    std::size_t idx = 0, mult = 1;
    for (std::size_t r = 0; r < sizes.size(); ++r) {
        idx += mult * static_cast<std::size_t>(v[r]);
        mult *= static_cast<std::size_t>(sizes[r]);
    }
    return idx;
}

/**
 * Partition function of the pseudonetwork at population k, normalized to
 * G(0) = 1 (MATLAB pseudonet).
 *
 * Populations never exceed 2*(terms-1) = 4, so the load-dependent convolution
 * runs directly over the mixed-radix lattice rather than through pfqn_gld.
 *
 * @param gam  (Mq x R) pseudonetwork demands
 * @param k    (R) population of the auxiliary lattice
 * @param mups (Mq x >= sum k) pseudonetwork rate lattice
 */
template <class T>
T pald_pseudonet(const Matrix<T>& gam, const std::vector<int>& k, const Matrix<T>& mups) {
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::size_t> nz;
    for (std::size_t r = 0; r < k.size(); ++r)
        if (k[r] > 0) nz.push_back(r);
    const std::size_t Mq = gam.rows(), Rp = nz.size();
    if (Rp == 0) return one;

    std::vector<int> sizes(Rp), kk(Rp);
    std::size_t npop = 1;
    for (std::size_t a = 0; a < Rp; ++a) {
        kk[a] = k[nz[a]];
        sizes[a] = kk[a] + 1;
        npop *= static_cast<std::size_t>(sizes[a]);
    }

    // Per-station balance terms over the whole lattice.
    Matrix<T> sterm(Mq, npop, zero);
    std::vector<int> m(Rp, 0);
    for (std::size_t i = 0; i < Mq; ++i)
        for (std::size_t jdx = 0; jdx < npop; ++jdx) {
            pald_idx2vec(jdx, sizes, m);
            int sm = 0;
            for (std::size_t r = 0; r < Rp; ++r) sm += m[r];
            T v = num_factln<T>(num_traits<T>::from_int(sm));
            bool finite = true;
            for (std::size_t r = 0; r < Rp && finite; ++r) {
                if (m[r] == 0) continue;
                const T g = gam(i, nz[r]);
                if (!(g > zero)) {
                    finite = false;
                    break;
                }
                v = T(v + num_traits<T>::from_int(m[r]) * log(g) -
                      num_factln<T>(num_traits<T>::from_int(m[r])));
            }
            if (!finite) continue;  // exp(-inf) = 0, the reference's break
            if (sm > 0)
                for (int s = 0; s < sm; ++s) v = T(v - log(mups(i, static_cast<std::size_t>(s))));
            sterm(i, jdx) = exp(v);
        }

    // Load-dependent convolution, one station per pass.
    std::vector<T> g(npop, zero), gnext(npop, zero);
    g[0] = one;
    std::vector<int> n(Rp, 0), diff(Rp, 0);
    for (std::size_t i = 0; i < Mq; ++i) {
        for (std::size_t idx = 0; idx < npop; ++idx) {
            pald_idx2vec(idx, sizes, n);
            T acc = zero;
            for (std::size_t jdx = 0; jdx < npop; ++jdx) {
                pald_idx2vec(jdx, sizes, m);
                bool fits = true;
                for (std::size_t r = 0; r < Rp; ++r) {
                    if (m[r] > n[r]) {
                        fits = false;
                        break;
                    }
                    diff[r] = n[r] - m[r];
                }
                if (!fits) continue;
                acc += g[pald_vec2idx(diff, sizes)] * sterm(i, jdx);
            }
            gnext[idx] = acc;
        }
        g = gnext;
    }
    return g[npop - 1];
}

}  // namespace detail

/**
 * @param L     (M x R) demands; an infinite-server row is recognized by its
 *              rate lattice and folded into the think time
 * @param N     (R) population
 * @param Z     (R) think times, already summed over the delay rows; empty for
 *              none
 * @param mu    (M x >= sum N) load-dependent rates; empty means all ones, and
 *              a short lattice is extended with its last column
 * @param terms 1, 2 or 3 terms of the normal-usage series, as in pfqn_panacea
 */
template <class T>
PanaceaLdResult<T> pfqn_panaceald(const Matrix<T>& L, const std::vector<int>& N,
                                  const std::vector<T>& Z, const Matrix<T>& mu, int terms) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_panaceald requires transcendental arithmetic (asymptotic expansion of "
                  "log G through the phi(n) transform)");
    using std::exp;
    using std::log;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_panaceald: L and N disagree on the class count");
    if (terms < 1 || terms > 3)
        throw InputError(
            "pfqn_panaceald: the terms parameter must be 1, 2 or 3 (higher-order coefficients are "
            "not implemented)");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);
    if (Zv.size() != R) throw InputError("pfqn_panaceald: Z has the wrong length");

    PanaceaLdResult<T> res;
    res.G = zero;
    res.lG = zero;
    res.normalUsage = true;
    res.reason = nullptr;
    const auto decline = [&](const char* why) {
        res.normalUsage = false;
        res.reason = why;
        res.G = zero;
        res.lG = zero;
        return res;
    };

    long Ntl = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) throw InputError("pfqn_panaceald: negative population");
        Ntl += N[r];
    }
    if (Ntl == 0) {
        res.G = one;
        res.lG = zero;
        return res;
    }
    const std::size_t Nt = static_cast<std::size_t>(Ntl);

    // The rate lattice, extended with its last column as the reference does.
    Matrix<T> mux(M, Nt, one);
    if (!mu.empty()) {
        if (mu.rows() != M)
            throw InputError("pfqn_panaceald: mu and L disagree on the station count");
        if (mu.cols() == 0) throw InputError("pfqn_panaceald: mu has no rate column");
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < Nt; ++k)
                mux(i, k) = mu(i, k < mu.cols() ? k : mu.cols() - 1);
    }

    // Type-3 rows: rate lattice exactly 1, 2, ..., Nt.
    const double fineTol = 1e-8;  // GlobalConstants.FineTol
    std::vector<bool> isIS(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        bool all = true;
        for (std::size_t k = 0; k < Nt && all; ++k)
            if (std::fabs(num_traits<T>::to_double(mux(i, k)) - static_cast<double>(k + 1)) >=
                fineTol)
                all = false;
        isIS[i] = all;
    }

    std::vector<T> Ztot = Zv;
    std::vector<std::size_t> qrows;
    for (std::size_t i = 0; i < M; ++i) {
        if (isIS[i]) {
            for (std::size_t r = 0; r < R; ++r) Ztot[r] += L(i, r);
        } else {
            qrows.push_back(i);
        }
    }
    const std::size_t Mq = qrows.size();

    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0 && !(Ztot[r] > zero))
            // No infinite server on the route of a populated class: the
            // expansion parameter rho_j0 is undefined and PANACEA does not
            // apply at all, normal usage or not.
            return decline(
                "a populated class visits no infinite server, so the expansion parameter rho_j0 "
                "is undefined");

    const auto delay_only = [&]() {
        T lG = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T nT = num_traits<T>::from_int(N[r]);
            lG -= detail::num_factln<T>(nT);
            lG += detail::pald_xlogy(nT, Ztot[r]);
        }
        res.lG = lG;
        res.G = exp(lG);
        return res;
    };
    if (Mq == 0) return delay_only();

    Matrix<T> Lq(Mq, R), muq(Mq, Nt);
    for (std::size_t a = 0; a < Mq; ++a) {
        for (std::size_t r = 0; r < R; ++r) Lq(a, r) = L(qrows[a], r);
        for (std::size_t k = 0; k < Nt; ++k) muq(a, k) = mux(qrows[a], k);
    }
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t k = 0; k < Nt; ++k) {
            const double v = num_traits<T>::to_double(muq(a, k));
            if (!(v > 0.0) || !std::isfinite(v))
                return decline("a load-dependent rate is not positive and finite");
        }

    // Offered load per queueing centre, sum_j K_j e_ji / rho_j0.
    Matrix<T> r(Mq, R, zero);
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t j = 0; j < R; ++j)
            if (Ztot[j] > zero) r(a, j) = T(Lq(a, j) / Ztot[j]);
    std::vector<T> lambda(Mq, zero), muK(Mq, one), alpha(Mq, one);
    for (std::size_t a = 0; a < Mq; ++a) {
        T s = zero;
        for (std::size_t j = 0; j < R; ++j) s += r(a, j) * num_traits<T>::from_int(N[j]);
        lambda[a] = s;
        muK[a] = muq(a, Nt - 1);
        alpha[a] = T(one - lambda[a] / muK[a]);
        if (!(alpha[a] > zero))
            return decline(
                "the model is not in normal usage (1 - lambda_i/mu_i(Ntot) <= 0 at some queueing "
                "centre), so the {phi(n)} series diverges");
    }

    // log prod_{k=1}^{s} mu_i(k), s = 0..Nt.
    std::vector<std::vector<T> > lPi(Mq, std::vector<T>(Nt + 1, zero));
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t s = 1; s <= Nt; ++s) lPi[a][s] = T(lPi[a][s - 1] + log(muq(a, s - 1)));

    const std::size_t nmax = static_cast<std::size_t>(2 * (terms - 1));
    Matrix<T> lpsi(Mq, nmax + 1, zero);
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t n = 0; n <= nmax; ++n)
            lpsi(a, n) = detail::pald_logpsi(n, lambda[a], lPi[a], muK[a], alpha[a], Nt);

    // Load dependence of the pseudonetwork centres:
    // psi_i(n) = psi_i(0) n! / prod_{k=1}^{n} mups_i(k).
    Matrix<T> mups(Mq, nmax > 0 ? nmax : 1, one);
    for (std::size_t a = 0; a < Mq; ++a)
        for (std::size_t n = 1; n <= nmax; ++n)
            mups(a, n - 1) = exp(T(log(num_traits<T>::from_int(static_cast<long>(n))) +
                                   lpsi(a, n - 1) - lpsi(a, n)));

    // Expansion coefficients (5.4). The large parameter N cancels identically
    // between beta_j = K_j/N, Gamma = N*r and the 1/N^n scaling, so the demands
    // are taken as r and beta as N.
    std::vector<T> A(3, zero);
    A[0] = one;
    if (terms >= 2) {
        for (std::size_t j = 0; j < R; ++j) {
            std::vector<int> m(R, 0);
            m[j] = 2;
            A[1] -= num_traits<T>::from_int(N[j]) * detail::pald_pseudonet(r, m, mups);
        }
    }
    if (terms >= 3) {
        for (std::size_t j = 0; j < R; ++j) {
            std::vector<int> m(R, 0);
            m[j] = 3;
            A[2] += num_traits<T>::from_int(2) * num_traits<T>::from_int(N[j]) *
                    detail::pald_pseudonet(r, m, mups);
            m.assign(R, 0);
            m[j] = 4;
            A[2] += num_traits<T>::from_int(3) *
                    T(num_traits<T>::from_int(N[j]) * num_traits<T>::from_int(N[j])) *
                    detail::pald_pseudonet(r, m, mups);
            for (std::size_t k = 0; k < R; ++k) {
                if (k == j) continue;
                m.assign(R, 0);
                m[j] = 2;
                m[k] = 2;
                A[2] += num_traits<T>::from_rational(1, 2) * num_traits<T>::from_int(N[j]) *
                        num_traits<T>::from_int(N[k]) * detail::pald_pseudonet(r, m, mups);
            }
        }
    }
    T I = zero;
    for (int t = 0; t < terms; ++t) I += A[t];
    if (!(I > zero))
        return decline("the truncated normal-usage series is not positive, so its logarithm is "
                       "undefined");

    T lG = zero;
    for (std::size_t j = 0; j < R; ++j) {
        const T nT = num_traits<T>::from_int(N[j]);
        lG -= detail::num_factln<T>(nT);
        lG += detail::pald_xlogy(nT, Ztot[j]);
    }
    for (std::size_t a = 0; a < Mq; ++a) lG += lpsi(a, 0);
    lG += log(I);
    if (!std::isfinite(num_traits<T>::to_double(lG)))
        return decline("the expansion evaluated to a non-finite logarithm");
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

/** Overload at the reference's default of three terms. */
template <class T>
PanaceaLdResult<T> pfqn_panaceald(const Matrix<T>& L, const std::vector<int>& N,
                                  const std::vector<T>& Z, const Matrix<T>& mu) {
    return pfqn_panaceald(L, N, Z, mu, 3);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_PANACEALD_H
