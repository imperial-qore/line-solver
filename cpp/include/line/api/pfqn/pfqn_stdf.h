/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_STDF_H
#define LINE_API_PFQN_STDF_H

/**
 * Sojourn-time distribution at multiserver FCFS stations of a closed
 * product-form network (J. McKenna, JACM 1987).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_stdf.m. There is no JAR
 * counterpart, so MATLAB is the only reference.
 *
 * Method. A class-r job arriving at station k sees, by the arrival theorem, the
 * network at population N - e_r. Conditional on finding n jobs already there,
 * its sojourn time is the sum of its own service and of the residual work of
 * the queue ahead, whose distribution at a station with S(k) servers is the
 * convolution
 *
 *   h_k(t | n) = Exp(rate_k)                                    n <  S(k)
 *   h_k(t | n) = Exp(rate_k) + Erlang_{n-S(k)+1}(S(k) rate_k)   n >= S(k),
 *
 * built here as MAPs and evaluated with map_cdf. Writing G_krt for the
 * transform of the sojourn CDF, the paper sums h_k(t | |nvec|) F_k(nvec)
 * G_k(N - e_r - nvec) over the whole population lattice. The reference keeps
 * that form as dead commented-out code and executes instead the equivalent
 * RECURSIVE form for load-dependent models, which is what this port implements:
 * the aggregate constant is re-evaluated on a rate lattice tilted by the ratio
 * of successive CDF levels,
 *
 *   gamma_k(t, n) = mu_k(n) h_k(t | n-1) / h_k(t | n),
 *
 * shifted by pfqn_mushift so that station k is one job ahead, giving
 *
 *   H_krt = h_k(t | 0) G_{-k}(N - e_r)
 *         + sum_s L(k,s) h_k(t | 0) / gamma_k(t,1) Y_ks(t),
 *   F(t)  = min(1, H_krt / G(N - e_r)),
 *
 * with Y_ks(t) the constant of the full model at population N - e_r - e_s on
 * the tilted lattice. The single-station case (M == 1) is dispatched to
 * pfqn_comomrm_ld, everything else to pfqn_mvald, exactly as the reference
 * dispatches it.
 *
 * Guards reproduced verbatim from the reference:
 *   - a time point equal to 0 is replaced by GlobalConstants.FineTol, because
 *     h_k(t | n) is not well defined there. FineTol is 1e-8, set by
 *     matlab/lineStart.m (FINE_TOL); it is NOT 1e-12;
 *   - a not-a-number H_krt becomes FineTol;
 *   - the result is clamped from above by min(1, .). pfqn_stdf_heur does NOT
 *     apply that clamp and can therefore report values above one; the asymmetry
 *     is real and is reproduced in both ports;
 *   - an FCFS station whose per-class service rates differ by more than FineTol
 *     is an invalid model and is rejected.
 *
 * Numerical stability. The reference warns when pfqn_mvald reports an unstable
 * marginal; this port returns that flag on the result instead of writing to a
 * log, since the port has no logging channel.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION, gated on has_transcendental. map_cdf is
 * a matrix exponential, the normalizing constants are combined in the log
 * domain, and the zero-time and not-a-number guards are floating-point
 * substitutions with no meaning in an exact field. Note also that pfqn_mvald
 * and pfqn_comomrm_ld report their logarithms as double, so the log-domain part
 * of this computation carries double precision whatever T is; T still governs
 * the CDF evaluation, the tilted rate lattice and the constants themselves.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_transform.h"
#include "line/api/pfqn/pfqn_comomrm_ld.h"
#include "line/api/pfqn/pfqn_mushift.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** GlobalConstants.FineTol, as set by matlab/lineStart.m. */
static const double kStdfFineTol = 1e-8;

/**
 * Result of pfqn_stdf / pfqn_stdf_heur, mirroring the MATLAB cell array RD.
 * RD[k][r] is a (T x 2) matrix whose first column is the sojourn-time CDF and
 * whose second column is the (guarded) time set, empty where the reference
 * leaves the cell unset.
 */
template <class T>
struct StdfResult {
    std::vector<std::vector<Matrix<T>>> RD;
    std::vector<T> tset;      ///< the time set after the zero guard
    bool isNumStable = true;  ///< false once an aggregate solve reported instability
};

namespace detail {

/**
 * Drop one row from a matrix (MATLAB's L(setdiff(1:M,k),:)). Dropping the only
 * row yields a 0 x cols matrix, not a 0 x 0 one: pfqn_mvald validates its rate
 * lattice by column count even with no stations, so the column count has to
 * survive.
 */
template <class T>
Matrix<T> stdf_drop_row(const Matrix<T>& A, std::size_t k) {
    Matrix<T> B(A.rows() - 1, A.cols());
    std::size_t o = 0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i == k) continue;
        for (std::size_t j = 0; j < A.cols(); ++j) B(o, j) = A(i, j);
        ++o;
    }
    return B;
}

/** The multiserver rate lattice mu(k,n) = min(S(k), n), n = 1, ..., sum(N). */
template <class T>
Matrix<T> stdf_mu(const std::vector<int>& S, int Nt) {
    Matrix<T> mu(S.size(), static_cast<std::size_t>(Nt));
    for (std::size_t k = 0; k < S.size(); ++k)
        for (int n = 1; n <= Nt; ++n)
            mu(k, static_cast<std::size_t>(n - 1)) =
                num_traits<T>::from_int(S[k] < n ? S[k] : n);
    return mu;
}

/** The time set with every zero replaced by FineTol. */
template <class T>
std::vector<T> stdf_guard_tset(const std::vector<T>& tset) {
    const T zero = num_traits<T>::from_int(0);
    const T fine = num_traits<T>::from_double(kStdfFineTol);
    std::vector<T> out = tset;
    for (std::size_t t = 0; t < out.size(); ++t) {
        if (out[t] < zero) throw InputError("pfqn_stdf: negative evaluation time");
        if (out[t] == zero) out[t] = fine;
    }
    return out;
}

/**
 * The log of the normalizing constant of a load-dependent closed model,
 * dispatched as the reference dispatches it: pfqn_comomrm_ld when there is at
 * most one queueing station, pfqn_mvald otherwise. Reports the stability flag
 * that only pfqn_mvald can lower.
 */
template <class T>
double stdf_lg(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
               const Matrix<T>& mu, bool singleStation, bool* stable) {
    if (singleStation) return pfqn_comomrm_ld(L, N, Z, mu).lG;
    const MvaLdResult<T> r = pfqn_mvald(L, N, Z, mu);
    if (!r.isNumStable && stable != nullptr) *stable = false;
    return r.lG;
}

/**
 * CDF of a hypoexponential with the given stage rates, evaluated WITHOUT
 * forming 1 - survival.
 *
 * WHY THIS EXISTS. map_cdf returns F(t) = 1 - pie exp(D0 t) e. For s in
 * [1/2, 2] the subtraction is exact by Sterbenz, so all of F's error is the
 * error already in s, which is of order p u whatever the size of F. At the
 * small probe times the sojourn-time recursion visits, F itself is far below
 * p u and the returned value is pure round-off: the deepest level CDFs come
 * back as +-2.2e-16 where the true values are 1e-25 or smaller, and the tilted
 * rate lattice, which divides by them, is then meaningless. Measured against a
 * 50-digit reference at t = 1e-6 the complement form is 7.6e-05 relative wrong
 * in the port and 2.1e-05 in MATLAB, at a level eight orders above eps and
 * nowhere near any noise floor.
 *
 * THE METHOD. Uniformize the pure-birth stage process at Lambda = max r_i and
 * carry the absorbing state explicitly: a_n, the probability that all p stages
 * are done within n uniformized steps, is nondecreasing in [0, 1] and is
 * accumulated from nonnegative flow, and
 *
 *   F(t) = sum_{n >= p} Poisson(n; Lambda t) a_n,
 *
 * a sum of nonnegative terms with no cancellation anywhere. Every stage
 * transition probability r_i / Lambda lies in [0, 1] and the self-loop
 * 1 - r_i / Lambda is exact for the largest rate and benign otherwise. The
 * series is truncated on a bound for the Poisson tail, itself a positive
 * quantity, so the truncation error is controlled rather than assumed.
 *
 * This is the same construction the JAR reaches for by a different route
 * (Foxglynn uniformization in Map_cdf, taken when D0 has no negative
 * off-diagonal entry). The C++ map_cdf header records that path as existing
 * "only for speed"; that is measurably wrong, and it is also the accurate
 * path at small t.
 *
 * DELIBERATE DIVERGENCE FROM MATLAB, authorized: where this disagrees with the
 * reference in this regime, the reference is the one that is wrong.
 */
template <class T>
T stdf_hypoexp_cdf(const std::vector<T>& r, const T& t) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t p = r.size();
    T lambda = r[0];
    for (std::size_t i = 1; i < p; ++i)
        if (lambda < r[i]) lambda = r[i];
    const T x = T(lambda * t);

    std::vector<T> step(p);
    for (std::size_t i = 0; i < p; ++i) step[i] = T(r[i] / lambda);

    std::vector<T> v(p, zero);
    v[0] = one;
    T a = zero;
    using std::exp;  // ADL takes the T-native exponential, not a double one
    T pois = exp(T(zero - x));
    T F = zero;
    const T tol = num_traits<T>::from_double(1e-20);
    // p stages cannot complete in fewer than p uniformized steps, so a_n is
    // zero until then and the sum only starts contributing at n = p.
    for (std::size_t n = 0; n <= 100000; ++n) {
        if (n > 0) {
            const T out = T(v[p - 1] * step[p - 1]);
            for (std::size_t i = p - 1; i >= 1; --i)
                v[i] = T(v[i] * T(one - step[i])) + T(v[i - 1] * step[i - 1]);
            v[0] = T(v[0] * T(one - step[0]));
            a += out;
            pois = T(pois * x) / num_traits<T>::from_int(static_cast<long>(n));
        }
        F += T(pois * a);
        if (n >= p) {
            const T nx = num_traits<T>::from_int(static_cast<long>(n) + 2);
            if (x < nx) {
                // the Poisson tail beyond n, bounded by its geometric majorant
                const T next = T(pois * x) / num_traits<T>::from_int(static_cast<long>(n) + 1);
                const T tail = T(next / T(one - T(x / nx)));
                if (tail <= T(F * tol) || tail == zero) return F;
            }
        }
    }
    throw NumericError("pfqn_stdf: the level CDF series did not converge");
}

/**
 * The level CDF h_k(t | n): Exp(rate) below the server count, and
 * Exp(rate) + Erlang_{n-S+1}(S rate) at or above it.
 *
 * THE TWO FORMS ARE NOT A FALLBACK PAIR: each is the accurate one in its own
 * regime, and the switch is placed where they exchange that role. The series
 * sums about Lambda t nonnegative terms, so its relative error grows with the
 * term count while the complement's shrinks as F leaves the round-off floor;
 * they cross where F stops being small, that is at Lambda t = p, the mean of
 * the stage process. Below it the complement returns round-off and the series
 * is exact to the last digit; above it the series is a few ulp low where the
 * complement is exact -- measured, at Lambda t = 44.8, as F = 1 - 2.2e-16
 * against a true value that rounds to exactly one.
 */
template <class T>
T stdf_level_cdf(const T& rate, int S, int n, const T& t) {
    std::vector<T> r;
    r.push_back(rate);
    if (n >= S) {
        const T sr = T(num_traits<T>::from_int(S) * rate);
        for (int j = 0; j < n - S + 1; ++j) r.push_back(sr);
    }
    T lambda = r[0];
    for (std::size_t i = 1; i < r.size(); ++i)
        if (lambda < r[i]) lambda = r[i];
    if (T(lambda * t) <= num_traits<T>::from_int(static_cast<long>(r.size())))
        return stdf_hypoexp_cdf(r, t);
    const T one = num_traits<T>::from_int(1);
    mam::Map<T> h;
    if (n < S) {
        h = mam::map_exponential_mean(T(one / rate));
    } else {
        std::vector<mam::Map<T>> parts;
        parts.push_back(mam::map_exponential_mean(T(one / rate)));
        parts.push_back(mam::map_erlang(
            T(num_traits<T>::from_int(n - S + 1) / (num_traits<T>::from_int(S) * rate)),
            static_cast<unsigned>(n - S + 1)));
        h = mam::map_sumind(parts);
    }
    std::vector<T> one_t;
    one_t.push_back(t);
    return mam::map_cdf(h, one_t)[0];
}

/** The tilted rate lattice gamma(t, .) and its mushifted, truncated form. */
template <class T>
void stdf_gamma(const Matrix<T>& mu, const Matrix<T>& hkc, std::size_t t, std::size_t k, int sumNr,
                bool truncate, Matrix<T>& gammat, Matrix<T>& gammak) {
    gammat = mu;
    for (int m = 1; m <= sumNr; ++m)
        gammat(k, static_cast<std::size_t>(m - 1)) =
            mu(k, static_cast<std::size_t>(m - 1)) *
            hkc(t, static_cast<std::size_t>(m - 1)) / hkc(t, static_cast<std::size_t>(m));
    gammak = pfqn_mushift(gammat, k);
    if (!truncate) return;
    // column-truncation width rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const std::size_t keep = sumNr >= 1 ? static_cast<std::size_t>(sumNr - 1) : 0;
    Matrix<T> g(gammak.rows(), keep);
    for (std::size_t i = 0; i < gammak.rows(); ++i)
        for (std::size_t j = 0; j < keep; ++j) g(i, j) = gammak(i, j);
    gammak = g;
}

}  // namespace detail

/**
 * Sojourn-time distribution at the listed FCFS stations.
 *
 * @param L         (M x R) service demands
 * @param N         (R) closed population vector
 * @param Z         (K x R) think times, summed over rows; may be empty
 * @param S         (M) server counts
 * @param fcfsNodes 0-based indices of the FCFS stations to analyze
 * @param rates     (M x R) service rates; the rates of an analyzed FCFS station
 *                  must agree across classes to within FineTol
 * @param tset      evaluation times; a zero entry is replaced by FineTol
 */
template <class T>
StdfResult<T> pfqn_stdf(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                        const std::vector<int>& S, const std::vector<std::size_t>& fcfsNodes,
                        const Matrix<T>& rates, const std::vector<T>& tset) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_stdf requires transcendental arithmetic: the level CDFs are matrix "
                  "exponentials and the normalizing constants are combined in the log domain");

    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    if (R != N.size()) throw InputError("pfqn_stdf: L and N disagree on the class count");
    if (S.size() != M) throw InputError("pfqn_stdf: S has the wrong length");
    if (rates.rows() != M || rates.cols() != R)
        throw InputError("pfqn_stdf: rates has the wrong shape");
    for (std::size_t k = 0; k < M; ++k)
        if (S[k] < 1) throw InputError("pfqn_stdf: the server count must be at least one");

    int Nt = 0;
    for (int n : N) {
        if (n < 0) throw InputError("pfqn_stdf: negative population");
        Nt += n;
    }
    if (Nt < 1) throw InputError("pfqn_stdf: the population must be at least one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T fine = num_traits<T>::from_double(kStdfFineTol);
    const std::size_t nt = tset.size();

    const Matrix<T> mu = detail::stdf_mu<T>(S, Nt);
    const std::vector<T> tv = detail::stdf_guard_tset(tset);

    StdfResult<T> res;
    res.tset = tv;
    res.RD.assign(M, std::vector<Matrix<T>>(R));

    const bool singleStation = (M == 1);

    for (std::size_t ki = 0; ki < fcfsNodes.size(); ++ki) {
        const std::size_t k = fcfsNodes[ki];
        if (k >= M) throw InputError("pfqn_stdf: FCFS station index out of range");

        T lo = rates(k, 0), hi = rates(k, 0);
        for (std::size_t r = 1; r < R; ++r) {
            if (rates(k, r) < lo) lo = rates(k, r);
            if (rates(k, r) > hi) hi = rates(k, r);
        }
        if (T(hi - lo) > fine)
            throw InputError(
                "pfqn_stdf: the FCFS station has distinct per-class service rates, the model is "
                "invalid");
        if (!(rates(k, 0) > zero))
            throw InputError("pfqn_stdf: the FCFS service rate must be strictly positive");

        // Level CDFs h_k(t | n), n = 0, ..., sum(N), computed without forming
        // 1 - survival: the tilted lattice divides by these and the complement
        // returns round-off at the probe times the recursion visits.
        Matrix<T> hkc(nt, static_cast<std::size_t>(Nt) + 1);
        for (int n = 0; n <= Nt; ++n)
            for (std::size_t t = 0; t < nt; ++t)
                hkc(t, static_cast<std::size_t>(n)) =
                    detail::stdf_level_cdf(rates(k, 0), S[k], n, tv[t]);

        const Matrix<T> Lk = detail::stdf_drop_row(L, k);
        const Matrix<T> muk = detail::stdf_drop_row(mu, k);

        for (std::size_t r = 0; r < R; ++r) {
            if (!(L(k, r) > fine)) continue;
            std::vector<int> Nr = N;
            Nr[r] -= 1;
            int sumNr = 0;
            for (int n : Nr) sumNr += n;

            const double lGr = detail::stdf_lg(L, Nr, Z, mu, singleStation, &res.isNumStable);
            const double lGk = detail::stdf_lg(Lk, Nr, Z, muk, singleStation, &res.isNumStable);

            Matrix<T> RD(nt, 2);
            for (std::size_t t = 0; t < nt; ++t) RD(t, 1) = tv[t];

            Matrix<T> gammat, gammak;
            for (std::size_t t = 0; t < nt; ++t) {
                detail::stdf_gamma(mu, hkc, t, k, sumNr, true, gammat, gammak);
                T H = hkc(t, 0) * num_traits<T>::from_double(std::exp(lGk));
                for (std::size_t s = 0; s < R; ++s) {
                    if (Nr[s] <= 0) continue;
                    std::vector<int> Nrs = Nr;
                    Nrs[s] -= 1;
                    const double lY =
                        detail::stdf_lg(L, Nrs, Z, gammak, singleStation, &res.isNumStable);
                    H += L(k, s) * hkc(t, 0) / gammat(k, 0) *
                         num_traits<T>::from_double(std::exp(lY));
                }
                if (!(H == H)) H = fine;  // the reference's isnan guard
                const double lH = num_traits<T>::log_as_double(H);
                const T v = num_traits<T>::from_double(std::exp(lH - lGr));
                RD(t, 0) = v > one ? one : v;  // min(1, .), which the heuristic omits
            }
            res.RD[k][r] = RD;
        }
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_STDF_H
