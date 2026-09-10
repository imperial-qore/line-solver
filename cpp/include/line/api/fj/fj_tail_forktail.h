/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_TAIL_FORKTAIL_H
#define LINE_API_FJ_FJ_TAIL_FORKTAIL_H

/**
 * ForkTail black-box tail-latency approximation for fork-join requests.
 *
 * Templated port of matlab/src/api/fj/fj_tail_forktail.m. No JAR counterpart.
 * Approximates the p-th percentile of the response time of a request that
 * forks into K parallel tasks and joins on the last of them, from the mean
 * and variance of the per-branch task response times ALONE. Each branch is a
 * black box whose response time is fitted by a generalized exponential law
 *
 *   F_T(x) = (1 - exp(-x/beta))^alpha,
 *   E[T] = beta (psi(alpha+1) - psi(1)),
 *   V[T] = beta^2 (psi'(1) - psi'(alpha+1)),
 *
 * and the request response time is the maximum over the branches, taken as
 * the PRODUCT of the branch CDFs. That product is exact only for independent
 * branches, which is the approximation the method rests on.
 *
 * Three routes, exactly as in MATLAB:
 *  - homogeneous (one branch mean, scalar K): the product of K identical CDFs
 *    raises the shape to K alpha and inverts in closed form,
 *    x_p = -beta log(1 - p^(1/(K alpha))).
 *  - random fanout (one branch mean, K a vector with probabilities P): the
 *    request law is the mixture sum_i P_i (1-exp(-x/beta))^(K_i alpha),
 *    bracketed between the closed forms at min(K) and max(K) and inverted
 *    numerically.
 *  - heterogeneous (a vector of branch means): solve
 *    sum_i alpha_i log1p(-exp(-x/beta_i)) = log p. F_X is bounded above by
 *    the CDF of any single branch, so the request percentile is at least the
 *    largest branch percentile, which is where the bracket search starts.
 *
 * This is a HIGH-LOAD result, from the central limit theorem for G/G/m queues
 * in heavy traffic: the reference reports errors within 20% and 15% at 80%
 * and 90% utilization and makes no claim at low load, where the tail is
 * dominated by the service law and the branch dependence is strongest.
 *
 * Reference: M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A
 * Black-Box Fork-Join Tail Latency Prediction Model for User-Facing
 * Datacenter Workloads", ACM HPDC 2018, pp. 206-217.
 *
 * ARITHMETIC: the generalized-exponential fit inverts a ratio of digamma and
 * trigamma values, so the exact instantiation is refused.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/rootfind.h"

namespace line {
namespace fj {

namespace detail {

/** Recurrence threshold: the Stirling tail is below 1e-17 from here on. */
inline int digamma_shift() { return 20; }

/**
 * psi(x) for x > 0, by upward recurrence to x >= 20 then the asymptotic series.
 * MATLAB's psi(0,x); the recurrence psi(x) = psi(x+1) - 1/x is exact, so the
 * only error is the truncation of the Stirling tail. The threshold is 20 and
 * not the more usual 6 or 8 because the first omitted term is B12/(12 x^12),
 * which is still 3e-13 at x = 8 and would cap the fit at six correct digits.
 */
template <class T>
T digamma(const T& x0) {
    const T one = num_traits<T>::from_int(1);
    const T lim = num_traits<T>::from_int(digamma_shift());
    T x = x0, acc = num_traits<T>::from_int(0);
    while (x < lim) {
        acc -= one / x;
        x += one;
    }
    const T inv = one / x, inv2 = inv * inv;
    using std::log;
    T s = log(x) - inv / num_traits<T>::from_int(2);
    T p = inv2;
    s -= p / num_traits<T>::from_int(12);
    p *= inv2;
    s += p / num_traits<T>::from_int(120);
    p *= inv2;
    s -= p / num_traits<T>::from_int(252);
    p *= inv2;
    s += p / num_traits<T>::from_int(240);
    p *= inv2;
    s -= p / num_traits<T>::from_int(132);
    p *= inv2;
    s += p * num_traits<T>::from_int(691) / num_traits<T>::from_int(32760);
    return acc + s;
}

/** psi'(x) for x > 0, the same construction on psi'(x) = psi'(x+1) + 1/x^2. */
template <class T>
T trigamma(const T& x0) {
    const T one = num_traits<T>::from_int(1);
    const T lim = num_traits<T>::from_int(digamma_shift());
    T x = x0, acc = num_traits<T>::from_int(0);
    while (x < lim) {
        acc += one / (x * x);
        x += one;
    }
    const T inv = one / x, inv2 = inv * inv;
    T s = inv + inv2 / num_traits<T>::from_int(2);
    T p = inv2 * inv;
    s += p / num_traits<T>::from_int(6);
    p *= inv2;
    s -= p / num_traits<T>::from_int(30);
    p *= inv2;
    s += p / num_traits<T>::from_int(42);
    p *= inv2;
    s -= p / num_traits<T>::from_int(30);
    p *= inv2;
    s += p * num_traits<T>::from_int(5) / num_traits<T>::from_int(66);
    return acc + s;
}

/**
 * Match a generalized exponential on a mean and a variance.
 *
 * The squared coefficient of variation depends on the SHAPE alone and
 * decreases monotonically in it, so the shape is recovered by a scalar root
 * find on a logarithmic scale and the scale then follows in closed form.
 * SCV = 1 is the exponential case alpha = 1 and is kept exact.
 */
template <class T>
void ge_fit(const T& ET, const T& VT, T& alpha, T& beta) {
    const T one = num_traits<T>::from_int(1);
    const T scv = VT / (ET * ET);
    const T fine = num_traits<T>::from_double(1e-12);
    if (num_abs(T(scv - one)) < fine) {
        alpha = one;
    } else {
        const T tri1 = trigamma(one);
        const T di1 = digamma(one);
        auto residual = [&](const T& la) {
            using std::exp;
            const T a = exp(la);
            const T d = digamma(T(a + one)) - di1;
            return T((tri1 - trigamma(T(a + one))) / (d * d) - scv);
        };
        T lo = num_traits<T>::from_int(-30), hi = num_traits<T>::from_int(30);
        const T lomin = num_traits<T>::from_int(-700), himax = num_traits<T>::from_int(700);
        const T step = num_traits<T>::from_int(30);
        const T zero = num_traits<T>::from_int(0);
        while (residual(lo) < zero && lo > lomin) lo -= step;  // smaller shape, larger SCV
        while (residual(hi) > zero && hi < himax) hi += step;  // larger shape, smaller SCV
        const RootResult<T> rr =
            root_brent<T>(residual, lo, hi, num_traits<T>::from_double(1e-14), 500);
        using std::exp;
        alpha = exp(rr.root);
    }
    beta = ET / (digamma(T(alpha + one)) - digamma(one));
}

}  // namespace detail

/** Mirrors MATLAB's [xp, alpha, beta] return list. */
template <class T>
struct ForkTailResult {
    T xp;                  ///< predicted p-th percentile of the request response time
    std::vector<T> alpha;  ///< fitted shape parameters, one per branch
    std::vector<T> beta;   ///< fitted scale parameters, one per branch
};

/**
 * @param ET per-branch mean task response times; one entry means homogeneous
 * @param VT per-branch variances, same length as ET
 * @param K  fanout(s); a vector of distinct fanouts needs P, ignored when ET
 *           has more than one entry
 * @param P  fanout probabilities, required when K has more than one entry
 * @param p  percentile, a fraction in (0,1) or a percentage in (0,100)
 */
template <class T>
ForkTailResult<T> fj_tail_forktail(const std::vector<T>& ET, const std::vector<T>& VT,
                                   const std::vector<T>& K = std::vector<T>(),
                                   const T& p_in = num_traits<T>::from_int(99),
                                   const std::vector<T>& P = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_tail_forktail requires transcendental arithmetic: the generalized "
                  "exponential fit inverts a ratio of digamma and trigamma values");
    using std::exp;
    using std::log;
    using std::log1p;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T hundred = num_traits<T>::from_int(100);

    std::vector<T> Kv = K;
    if (Kv.empty()) Kv.assign(1, one);
    T p = p_in;
    if (p > one) p = p / hundred;
    if (p <= zero || p >= one)
        throw InputError("fj_tail_forktail: the percentile must lie strictly between 0 and 1");
    if (VT.size() != ET.size())
        throw InputError("fj_tail_forktail: ET and VT must have the same number of entries");
    if (ET.empty()) throw InputError("fj_tail_forktail: ET must be nonempty");
    for (std::size_t i = 0; i < ET.size(); ++i)
        if (ET[i] <= zero || VT[i] <= zero)
            throw InputError(
                "fj_tail_forktail: the task response time mean and variance must be positive");

    const std::size_t nbranch = ET.size();
    ForkTailResult<T> r;
    r.alpha.assign(nbranch, zero);
    r.beta.assign(nbranch, zero);
    for (std::size_t i = 0; i < nbranch; ++i) detail::ge_fit(ET[i], VT[i], r.alpha[i], r.beta[i]);

    if (nbranch == 1 && Kv.size() > 1) {
        // random fanout: mix the homogeneous request laws and invert numerically
        if (P.size() != Kv.size())
            throw InputError(
                "fj_tail_forktail: a vector of fanouts K needs a probability vector P of the "
                "same length");
        T sp = zero;
        for (std::size_t i = 0; i < P.size(); ++i) {
            if (P[i] < zero)
                throw InputError("fj_tail_forktail: the fanout probabilities P must be "
                                 "non-negative and sum to one");
            sp += P[i];
        }
        if (num_abs(T(sp - one)) > num_traits<T>::from_double(1e-6))
            throw InputError(
                "fj_tail_forktail: the fanout probabilities P must be non-negative and sum to one");
        const T alpha = r.alpha[0], beta = r.beta[0];
        T kmin = Kv[0], kmax = Kv[0];
        for (std::size_t i = 1; i < Kv.size(); ++i) {
            if (Kv[i] < kmin) kmin = Kv[i];
            if (Kv[i] > kmax) kmax = Kv[i];
        }
        auto mixres = [&](const T& x) {
            T s = zero;
            for (std::size_t i = 0; i < Kv.size(); ++i)
                s += P[i] * exp(Kv[i] * alpha * log1p(-exp(-x / beta)));
            return T(s - p);
        };
        const T xlo = -beta * log1p(-exp(log(p) / (kmin * alpha)));
        const T xhi = -beta * log1p(-exp(log(p) / (kmax * alpha)));
        const T lo = xlo < xhi ? xlo : xhi;
        const T hi = xlo < xhi ? xhi : xlo;
        if (lo == hi) {
            r.xp = lo;
            return r;
        }
        // G(x)^(K alpha) decreases in K, so the mixture obeys mixres(lo) <= 0
        // <= mixres(hi) exactly -- with EQUALITY when P puts all its mass on
        // kmin or on kmax. There the root sits ON an endpoint, the residual
        // there is roundoff of either sign rather than the strict straddle
        // root_brent demands, and the endpoint is already the answer.
        const T flo = mixres(lo);
        const T fhi = mixres(hi);
        if (!(flo < zero)) {
            r.xp = lo;
            return r;
        }
        if (!(fhi > zero)) {
            r.xp = hi;
            return r;
        }
        const RootResult<T> rr = root_brent<T>(
            mixres, lo, hi, num_traits<T>::from_double(1e-14) * (one + hi), 500);
        r.xp = rr.root;
        return r;
    }

    if (nbranch == 1) {
        // homogeneous: K identical CDFs raise the shape to K alpha
        r.xp = -r.beta[0] * log1p(-exp(log(p) / (Kv[0] * r.alpha[0])));
        return r;
    }

    // heterogeneous: solve prod_i (1-exp(-x/beta_i))^alpha_i = p on the log scale
    const T logp = log(p);
    auto residual = [&](const T& x) {
        T s = zero;
        for (std::size_t i = 0; i < nbranch; ++i)
            s += r.alpha[i] * log1p(-exp(-x / r.beta[i]));
        return T(s - logp);
    };
    T xlo = zero;
    for (std::size_t i = 0; i < nbranch; ++i) {
        const T cand = -r.beta[i] * log1p(-exp(logp / r.alpha[i]));
        if (cand > xlo) xlo = cand;
    }
    T xhi = xlo;
    while (residual(xhi) < zero) {
        xhi = two * xhi;
        if (!std::isfinite(num_traits<T>::to_double(xhi)))
            throw InputError("fj_tail_forktail: could not bracket the ForkTail percentile");
    }
    const T tiny = num_traits<T>::from_double(2.220446049250313e-16);
    if (residual(xlo) > zero) {
        xlo = xlo / two;
        while (residual(xlo) > zero && xlo > tiny) xlo = xlo / two;
    }
    const RootResult<T> rr =
        root_brent<T>(residual, xlo, xhi, num_traits<T>::from_double(1e-14) * (one + xhi), 500);
    r.xp = rr.root;
    return r;
}

/** Homogeneous convenience overload with a scalar mean, variance and fanout. */
template <class T>
ForkTailResult<T> fj_tail_forktail(const T& ET, const T& VT, const T& K,
                                   const T& p = num_traits<T>::from_int(99)) {
    return fj_tail_forktail(std::vector<T>(1, ET), std::vector<T>(1, VT), std::vector<T>(1, K), p);
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_TAIL_FORKTAIL_H
