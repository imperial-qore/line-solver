/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIGK_RQT_H
#define LINE_API_QSYS_QSYS_GIGK_RQT_H

/**
 * Robust Queueing Theory (RQT) worst-case system time of a G/G/k FCFS queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gigk_rqt.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gigk_rqt.java.
 *
 * The arrival and service processes are not described by distributions but by
 * the polyhedral uncertainty sets
 *
 *   U^a = { T : (sum_{i=k+1}^n T_i - (n-k)/lambda)/(n-k)^(1/alpha_a) >= -Gamma_a }
 *   U^s = { X : (sum_{i=k}^n X_i - (n-k+1)/mu)/(n-k+1)^(1/alpha_s) <= Gamma_s }
 *
 * whose shape follows the (generalized) central limit theorem: alpha = 2 is the
 * finite-variance regime, alpha in (1,2) the heavy-tailed one. Performance
 * analysis is then a worst-case optimization rather than an expectation.
 *
 * W is the closed-form bound of Theorem 3 (Theorem 8 when the two tails differ,
 * with ab = min(alpha_a,alpha_s)),
 *
 *   W <= (ab-1)/ab^(ab/(ab-1)) lambda^(1/(ab-1))
 *        (Gamma_a + Gamma_s/k^(1/ab))^(ab/(ab-1)) / (1-rho)^(1/(ab-1)) + k/lambda,
 *
 * which for k = 1 reduces to Theorem 2 and, at ab = 2, to the Kingman-like form
 * (lambda/4)(Gamma_a+Gamma_s)^2/(1-rho) + 1/lambda. Sworst is the exact worst
 * case over the uncertainty sets, eq. (45): the supremum over the integer
 * x = nu-j+1 >= 1 of
 *
 *   x/mu + Gamma_s x^(1/alpha_s) - k(x-1)/lambda + Gamma_a (k(x-1))^(1/alpha_a).
 *
 * The arrival deviation ADDS to the worst case, since the adversary shortens the
 * interarrival times; the sign printed in eq. (12) is easily misread as a
 * subtraction of the whole arrival bracket, and reading it that way puts Sworst
 * an order of magnitude below W.
 *
 * W is a SYSTEM time (waiting plus service), and its additive term is k/lambda
 * rather than the mean service time 1/mu.
 *
 * ARITHMETIC. Real exponents make this transcendental. At rho >= 1 MATLAB
 * returns Inf; the port raises instead, as the rest of the qsys port does.
 *
 * Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
 * Operations Research 63(3), 676-700.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct GigkRqtResult {
    T W;       ///< closed-form bound on the system time (Theorem 3 / Theorem 8)
    T rhohat;  ///< modified utilization, so that M/M/1 relations still hold
    T Sworst;  ///< exact worst-case system time over the uncertainty sets
};

/**
 * @param lambda  arrival rate
 * @param mu      service rate of each server
 * @param Gamma_a variability parameter of the arrival uncertainty set
 * @param Gamma_s variability parameter of the service uncertainty set
 * @param k       number of servers
 * @param alpha_a arrival tail coefficient in (1,2]
 * @param alpha_s service tail coefficient in (1,2]
 */
template <class T>
GigkRqtResult<T> qsys_gigk_rqt(const T& lambda, const T& mu, const T& Gamma_a, const T& Gamma_s,
                               std::size_t k, const T& alpha_a, const T& alpha_s) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_gigk_rqt requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T kk = num_traits<T>::from_int(static_cast<int>(k));
    if (alpha_a <= one || alpha_a > T(one + one) || alpha_s <= one || alpha_s > T(one + one))
        throw InputError("qsys_gigk_rqt: tail coefficients must lie in (1,2]");

    GigkRqtResult<T> r;
    const T rho = lambda / (kk * mu);
    if (lambda <= zero) {
        r.W = one / mu;
        r.rhohat = zero;
        r.Sworst = one / mu;
        return r;
    }
    if (rho >= one)
        throw InputError("qsys_gigk_rqt: rho must be strictly less than 1 for a finite system time");

    // Theorem 8 collapses to Theorem 3 when the two tails agree.
    const T ab = alpha_a < alpha_s ? alpha_a : alpha_s;
    const T beta = Gamma_a + Gamma_s / detail::num_pow(kk, T(one / ab));
    if (beta <= zero) {
        // a nonpositive effective variability leaves only the deterministic term
        r.W = kk / lambda;
    } else {
        const T e = ab / (ab - one);
        r.W = (ab - one) / detail::num_pow(ab, e) * detail::num_pow(lambda, T(one / (ab - one))) *
                  detail::num_pow(beta, e) / detail::num_pow(T(one - rho), T(one / (ab - one))) +
              kk / lambda;
    }
    r.rhohat = detail::rhohat_from_W(r.W, lambda);

    // Exact worst case, eq. (45), over the integer lattice x >= 1.
    auto obj = [&](const T& x) -> T {
        T y = T(x - one);
        if (y < zero) y = zero;
        return x / mu + Gamma_s * detail::num_pow(x, T(one / alpha_s)) - kk * y / lambda +
               Gamma_a * detail::num_pow(T(kk * y), T(one / alpha_a));
    };
    // the continuous maximizer of the bounding problem, eq. (16), sizes the scan
    double xstar = 1.0;
    if (beta > zero) {
        const T xs = detail::num_pow(T(lambda * beta / (ab * (one - rho))), T(ab / (ab - one)));
        xstar = num_traits<T>::to_double(xs);
    }
    double xhi = std::max(4.0, std::ceil(4.0 * xstar));
    if (!(xhi > 0.0) || !std::isfinite(xhi)) xhi = 4.0;
    const std::size_t NS = 400;
    std::vector<double> xs;
    xs.reserve(NS);
    for (std::size_t i = 0; i < NS; ++i) {
        const double e = std::log10(xhi) * static_cast<double>(i) / static_cast<double>(NS - 1);
        const double x = std::max(1.0, std::floor(std::pow(10.0, e) + 0.5));
        if (xs.empty() || x != xs.back()) xs.push_back(x);  // unique(round(logspace(...)))
    }
    std::size_t imax = 0;
    T best = obj(T(num_traits<T>::from_double(xs[0])));
    for (std::size_t i = 1; i < xs.size(); ++i) {
        const T v = obj(T(num_traits<T>::from_double(xs[i])));
        if (v > best) {
            best = v;
            imax = i;
        }
    }
    // Refine on the continuous relaxation, then round back onto the lattice.
    T lo = num_traits<T>::from_double(xs[imax > 0 ? imax - 1 : 0]);
    T hi = num_traits<T>::from_double(xs[imax + 1 < xs.size() ? imax + 1 : xs.size() - 1]);
    if (hi > lo) {
        const T tolx = num_traits<T>::from_double(1e-8);
        const T invphi = num_traits<T>::from_double(0.6180339887498949);
        T c = hi - (hi - lo) * invphi;
        T d = lo + (hi - lo) * invphi;
        T fc = obj(c), fd = obj(d);
        for (unsigned it = 0; it < 500u && hi - lo > tolx; ++it) {
            if (fc > fd) {
                hi = d;
                d = c;
                fd = fc;
                c = hi - (hi - lo) * invphi;
                fc = obj(c);
            } else {
                lo = c;
                c = d;
                fc = fd;
                d = lo + (hi - lo) * invphi;
                fd = obj(d);
            }
        }
        const double xc = num_traits<T>::to_double(T((lo + hi) / (one + one)));
        const double cand[2] = {std::floor(xc), std::ceil(xc)};
        for (int i = 0; i < 2; ++i) {
            if (cand[i] >= 1.0) {
                const T v = obj(T(num_traits<T>::from_double(cand[i])));
                if (v > best) best = v;
            }
        }
    }
    r.Sworst = best;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIGK_RQT_H
