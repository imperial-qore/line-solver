/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_TAIL_ORDSTAT_H
#define LINE_API_FJ_FJ_TAIL_ORDSTAT_H

/**
 * Tail latency of a k-of-n (QUORUM) fork-join request.
 *
 * Templated port of matlab/src/api/fj/fj_tail_ordstat.m, mirrored by
 * jline.api.fj.FJ_tail_ordstat. The request forks into N parallel tasks and
 * joins on the KREQ-th of them: KREQ = N is the ordinary AND-join and
 * reproduces `fj_tail_forktail` EXACTLY, KREQ = 1 is the first completion.
 *
 * Each branch is the same black box ForkTail uses -- its task response time is
 * fitted by a generalized exponential law F_i(x) = (1-exp(-x/beta_i))^alpha_i
 * matched on the branch mean and variance by `detail::ge_fit` -- and the
 * request completes once KREQ of the N branches have, so its law is the
 * KREQ-th ORDER STATISTIC of independent, not identically distributed branch
 * times,
 *
 *   F_X(x) = P(at least KREQ of the N branches are done by x),
 *
 * the upper tail of a Poisson-binomial with success probabilities F_i(x).
 * It is evaluated by the convolution recurrence over the branches, whose terms
 * are all non-negative so it adds no cancellation, and inverted by bisection.
 *
 * THE BRACKET IS STRUCTURAL, not searched: X_(KREQ) lies between the FIRST
 * completion and the LAST, so the smallest single-branch percentile bounds it
 * below and the AND-join percentile -- which is what `fj_tail_forktail`
 * returns -- bounds it above.
 *
 * At KREQ = N the recurrence reduces TERM BY TERM to prod_i F_i(x), so a full
 * join evaluates exactly as it did before this header existed; the code returns
 * the ForkTail root directly there rather than re-deriving it.
 *
 * BRANCH INDEPENDENCE is assumed, as in ForkTail: the branches of one request
 * are positively correlated through their shared arrival instant, so the true
 * quorum percentile is somewhat larger than this one. The same heavy-traffic
 * caveat applies, see fj_tail_forktail.h.
 *
 * References: M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A
 * Black-Box Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
 * Workloads", ACM HPDC 2018, for the branch law; A. Thomasian, "Analysis of
 * Fork/Join and Related Queueing Systems", ACM Computing Surveys 47(2),
 * Article 17, 2014, Sec. 3, for the quorum.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/fj/fj_tail_forktail.h"

namespace line {
namespace fj {

namespace detail {

/**
 * The generalized-exponential CDF (1-exp(-x/beta))^alpha, clamped into [0,1]:
 * a rounding excursion above one would make the recurrence below emit a
 * negative complement.
 */
template <class T>
T ge_cdf(const T& x, const T& alpha, const T& beta) {
    using std::exp;
    using std::log1p;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (!(x > zero)) return zero;
    const T u = exp(alpha * log1p(-exp(-x / beta)));
    if (!std::isfinite(num_traits<T>::to_double(u))) return zero;
    if (u < zero) return zero;
    if (u > one) return one;
    return u;
}

/**
 * P(at least `kreq` successes) for independent Bernoulli trials of success
 * probabilities `u`, by the convolution recurrence over the trials.
 */
template <class T>
T poissbin_upper(const std::vector<T>& u, std::size_t kreq) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t n = u.size();
    std::vector<T> pmf(n + 1, zero);
    pmf[0] = one;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j >= 1; --j)
            pmf[j] = T(pmf[j] * (one - u[i]) + pmf[j - 1] * u[i]);
        pmf[0] = T(pmf[0] * (one - u[i]));
    }
    T acc = zero;
    for (std::size_t j = kreq; j <= n; ++j) acc += pmf[j];
    return acc;
}

}  // namespace detail

/**
 * @param ET   per-branch mean task response times; one entry means homogeneous
 * @param VT   per-branch variances, same length as ET
 * @param K    fanout, used only when ET has a single entry
 * @param p_in percentile, a fraction in (0,1) or a percentage in (0,100)
 * @param kreq the join fires on the kreq-th branch; 0 means every branch
 */
template <class T>
ForkTailResult<T> fj_tail_ordstat(const std::vector<T>& ET, const std::vector<T>& VT,
                                  std::size_t K = 1,
                                  const T& p_in = num_traits<T>::from_int(99),
                                  std::size_t kreq = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_tail_ordstat requires transcendental arithmetic: the generalized "
                  "exponential fit inverts a ratio of digamma and trigamma values");
    using std::exp;
    using std::log;
    using std::log1p;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T hundred = num_traits<T>::from_int(100);

    if (ET.size() != VT.size())
        throw InputError("fj_tail_ordstat: ET and VT must have the same number of entries");
    if (ET.empty()) throw InputError("fj_tail_ordstat: ET must be nonempty");

    T p = p_in;
    if (p > one) p = T(p / hundred);
    if (!(p > zero) || !(p < one))
        throw InputError("fj_tail_ordstat: the percentile must lie strictly between 0 and 1");
    for (std::size_t i = 0; i < ET.size(); ++i)
        if (!(ET[i] > zero) || !(VT[i] > zero))
            throw InputError(
                "fj_tail_ordstat: the task response time mean and variance must be positive");

    const std::size_t nbranch = ET.size();
    const std::size_t nsib = (nbranch == 1) ? (K < 1 ? 1 : K) : nbranch;
    if (kreq == 0) kreq = nsib;
    if (kreq > nsib)
        throw InputError("fj_tail_ordstat: the quorum must not exceed the branch count");

    // A homogeneous request over nsib branches is the same order statistic as a
    // heterogeneous one whose branches all carry the same moments, and expanding
    // it keeps ONE recurrence instead of a second closed form to keep in step.
    std::vector<T> et = ET, vt = VT;
    if (nbranch == 1 && nsib > 1) {
        et.assign(nsib, ET[0]);
        vt.assign(nsib, VT[0]);
    }

    ForkTailResult<T> r;
    r.alpha.assign(et.size(), zero);
    r.beta.assign(et.size(), zero);
    for (std::size_t i = 0; i < et.size(); ++i)
        detail::ge_fit(et[i], vt[i], r.alpha[i], r.beta[i]);

    // The AND-join percentile bounds every quorum above, and the SMALLEST
    // single-branch percentile bounds it below.
    std::vector<T> kv(1, num_traits<T>::from_int(static_cast<long>(nsib)));
    const ForkTailResult<T> andjoin = fj_tail_forktail<T>(ET, VT, kv, p);
    if (kreq == nsib) return andjoin;
    const T xhi = andjoin.xp;

    T xlo = T(-1);
    for (std::size_t i = 0; i < r.alpha.size(); ++i) {
        const T cand = -r.beta[i] * log1p(-exp(log(p) / r.alpha[i]));
        if (!(xlo > zero) || cand < xlo) xlo = cand;
    }

    auto residual = [&](const T& x) {
        std::vector<T> u(r.alpha.size(), zero);
        for (std::size_t i = 0; i < r.alpha.size(); ++i)
            u[i] = detail::ge_cdf(x, r.alpha[i], r.beta[i]);
        return T(detail::poissbin_upper(u, kreq) - p);
    };

    // A single branch may already meet the percentile; the quorum is met earlier.
    const T tiny = num_traits<T>::from_double(2.220446049250313e-16);
    while (residual(xlo) > zero && xlo > tiny) xlo = T(xlo / two);

    const RootResult<T> rr =
        root_brent<T>(residual, xlo, xhi, num_traits<T>::from_double(1e-14) * (one + xhi), 500);
    r.xp = rr.root;
    return r;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_TAIL_ORDSTAT_H
