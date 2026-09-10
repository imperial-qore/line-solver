/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GG1_H
#define LINE_API_QSYS_QSYS_GG1_H

/**
 * G/G/1 dispatcher: exact where a two-moment description determines the
 * answer, Allen-Cunneen otherwise.
 *
 * Templated port of matlab/src/api/qsys/qsys_gg1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gg1.java.
 *
 *   ca2 = cs2 = 1        -> qsys_mm1                       (exact)
 *   ca2 = 1              -> qsys_mg1 with cs = sqrt(cs2)   (exact)
 *   cs2 = 1              -> qsys_gm1 at the G/M/1 root     (exact)
 *   otherwise            -> qsys_gig1_approx_allencunneen  (approximation)
 *
 * In the G/M/1 branch the interarrival law is fitted from (lambda, ca2) by a
 * two-moment renewal process -- a balanced-means H2 for ca2 > 1, a Tijms
 * mixture of Erlang-(j-1)/Erlang-j for ca2 < 1, deterministic below 1e-6 --
 * and sigma is the root in (0,1) of sigma = A*(mu(1-sigma)) with A* the
 * interarrival LST. The map T(x) = A*(mu(1-x)) is increasing and the queue
 * root is its smallest fixed point, so the iterates from sigma_0 = rho
 * converge monotonically; MATLAB runs at most 100000 of them and stops at an
 * absolute step below 1e-13, which is reproduced here.
 *
 * ARITHMETIC. The fixed point is driven to a tolerance and the deterministic
 * branch evaluates exp, so the whole function is gated on transcendental
 * arithmetic. The three exact branches are individually available at Rational
 * through qsys_mm1, qsys_mg1 and qsys_gm1; only the dispatcher, whose tolerance
 * test |ca2-1| < 1e-8 is itself inexact, is gated.
 *
 * MATLAB-vs-JAR. The JAR returns L, Lq, W, Wq, p0 (and, in a second overload,
 * a geometric pk) whereas MATLAB returns [W, rhohat]. The W they compute is
 * the same; the port follows the MATLAB return list, since rhohat is what the
 * rest of the qsys family consumes.
 */

#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gm1.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

namespace detail {

/**
 * Root in (0,1) of sigma = A*(mu(1-sigma)) for the two-moment fit of the
 * interarrival LST A*. MATLAB's local qsys_gm1_sigma.
 */
template <class T>
T qsys_gm1_sigma(const T& lambda, const T& mu, const T& ca2) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gm1_sigma requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T tiny = T(num_traits<T>::from_double(1e-6));
    const T step_tol = T(num_traits<T>::from_double(1e-13));

    unsigned jj = 0;
    T p = num_traits<T>::from_int(0), nu = num_traits<T>::from_int(0);
    T p1 = num_traits<T>::from_int(0), l1 = num_traits<T>::from_int(0),
      l2 = num_traits<T>::from_int(0);
    if (ca2 >= one) {
        // hyperexponential H2 with balanced means
        p1 = (one + num_sqrt(T((ca2 - one) / (ca2 + one)))) / two;
        l1 = two * p1 * lambda;
        l2 = two * (one - p1) * lambda;
    } else if (ca2 >= tiny) {
        // mixed-Erlang phase-count rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
        const double inv = 1.0 / num_traits<T>::to_double(ca2);
        jj = static_cast<unsigned>(std::ceil(inv));
        if (jj < 1) jj = 1;
        const T jt = num_traits<T>::from_int(static_cast<long>(jj));
        p = (jt * ca2 - num_sqrt(T(jt * (one + ca2) - jt * jt * ca2))) / (one + ca2);
        nu = (jt - p) * lambda;
    }

    T sigma = lambda / mu;
    for (unsigned it = 0; it < 100000u; ++it) {
        const T s = mu * (one - sigma);
        T signew;
        if (ca2 < tiny) {
            signew = num_exp(T(-s / lambda));  // deterministic interarrival times
        } else if (ca2 < one) {
            signew = p * num_pow_int(T(nu / (s + nu)), jj - 1) +
                     (one - p) * num_pow_int(T(nu / (s + nu)), jj);
        } else {
            signew = p1 * l1 / (s + l1) + (one - p1) * l2 / (s + l2);
        }
        const bool done = num_abs(T(signew - sigma)) < step_tol;
        sigma = signew;
        if (done) break;
    }
    return sigma;
}

}  // namespace detail

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param ca2    squared coefficient of variation of the interarrival time
 * @param cs2    squared coefficient of variation of the service time
 */
template <class T>
QsysResult<T> qsys_gg1(const T& lambda, const T& mu, const T& ca2, const T& cs2) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gg1 requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T tol = T(num_traits<T>::from_double(1e-8));
    const bool ca_markov = num_abs(T(ca2 - one)) < tol;
    const bool cs_markov = num_abs(T(cs2 - one)) < tol;

    if (ca_markov && cs_markov) return qsys_mm1(lambda, mu);
    if (ca_markov) return qsys_mg1(lambda, mu, T(detail::num_sqrt(cs2)));
    if (cs_markov) {
        const T sigma = detail::qsys_gm1_sigma(lambda, mu, ca2);
        const T W = qsys_gm1(sigma, mu);
        return {W, detail::rhohat_from_W(W, lambda)};
    }
    return qsys_gig1_approx_allencunneen(lambda, mu, T(detail::num_sqrt(ca2)),
                                         T(detail::num_sqrt(cs2)));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GG1_H
