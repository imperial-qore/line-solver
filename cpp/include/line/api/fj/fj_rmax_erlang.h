/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RMAX_ERLANG_H
#define LINE_API_FJ_RMAX_ERLANG_H

/**
 * Expected maximum of K M/E_k/1 branch response times.
 *
 * Templated port of matlab/src/api/fj/fj_rmax_erlang.m.
 *
 * The branch mean response time comes from Pollaczek-Khinchine with the
 * Erlang-k SCV 1/k,
 *
 *   R = (k/mu) [1 + rho (1 + 1/k) / (2 (1 - rho))],   rho = lambda k / mu
 *
 * At K = 2 the maximum has the closed form (Thomasian 2014, Eq. 34)
 *
 *   R_2^max = 2 R - sum_{m,n < k} C(m+n,m) mu_R^{m+n} / (2 mu_R)^{m+n+1}
 *
 * with mu_R = k / R the rate that matches the response-time mean. For any
 * other K the MATLAB file moment-matches the response time to an Erlang and
 * integrates 1 - F(t)^K numerically.
 *
 * MIXED ARITHMETIC. The K = 2 branch is rational and exact in any field; the
 * general-K branch needs exp and a quadrature and throws UnsupportedError at
 * exact arithmetic.
 *
 * REFERENCE DEFECT: FJ_rmax.fj_rmax_erlang in
 * jar/src/main/java/jline/api/fj/FJ_rmax.java drops the mu_R^{m+n} numerator
 * from the K = 2 correction, computing C(m+n,m)/(2 mu_R)^{m+n+1} instead. The
 * two agree only at mu_R = 1. MATLAB is ground truth and is what this port
 * follows. Note also that the MATLAB correction simplifies to
 * sum C(m+n,m) / (2^{m+n+1} mu_R), which is how it is evaluated here.
 */

#include "line/api/fj/fj_types.h"
#include "line/api/fj/fj_xmax_erlang.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K      number of branches, K >= 1
 * @param k      Erlang stages of the branch service time, k >= 1
 * @param lambda arrival rate
 * @param mu     per-stage service rate (branch mean service is k/mu)
 * @return       expected maximum of the K branch response times
 */
template <class T>
T fj_rmax_erlang(unsigned K, unsigned k, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_rmax_erlang");
    if (k < 1) throw InputError("fj_rmax_erlang: the stage count k must be a positive integer");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T kt = num_traits<T>::from_int(static_cast<long>(k));

    const T mean_service = kt / mu;
    const T rho = lambda * mean_service;
    if (rho >= one) throw NumericError("fj_rmax_erlang: unstable system, rho >= 1");

    const T cv2 = one / kt;
    const T R_single = mean_service * (one + rho * (one + cv2) / (two * (one - rho)));

    if (K == 2) {
        const T mu_resp = kt / R_single;
        T correction = num_traits<T>::from_int(0);
        for (unsigned m = 0; m < k; ++m)
            for (unsigned n = 0; n < k; ++n)
                correction += detail::fj_binom<T>(m + n, m) * num_pow_int(mu_resp, m + n) /
                              num_pow_int(T(two * mu_resp), m + n + 1);
        return two * R_single - correction;
    }

    if constexpr (num_traits<T>::has_transcendental) {
        // MATLAB's cv2_response: (1/k + rho)/(1 + rho), floored at 1/20.
        const T cv2r_raw = (cv2 + rho) / (one + rho);
        const T floor20 = num_traits<T>::from_rational(1, 20);
        const T cv2r = cv2r_raw > floor20 ? cv2r_raw : floor20;
        const double inv = num_traits<T>::to_double(T(one / cv2r));
        unsigned k_resp = static_cast<unsigned>(std::ceil(inv));
        if (k_resp < 1) k_resp = 1;
        const T mu_resp = num_traits<T>::from_int(static_cast<long>(k_resp)) / R_single;
        const T upper = R_single * num_traits<T>::from_int(20);
        return detail::simpson<T>(
            [&](const T& t) { return T(one - num_pow_int(detail::erlang_cdf(t, k_resp, mu_resp), K)); },
            num_traits<T>::from_int(0), upper);
    } else {
        throw UnsupportedError(
            "fj_rmax_erlang: only K = 2 has a closed form; any other branch count needs a "
            "quadrature of the fitted Erlang response-time CDF and therefore transcendental "
            "arithmetic");
    }
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RMAX_ERLANG_H
