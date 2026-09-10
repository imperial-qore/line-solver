/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_TANDEM_UB_CIUCU_H
#define LINE_API_QSYS_QSYS_TANDEM_UB_CIUCU_H

/**
 * Tail bounds for a GI/Hn/1 -> ./Hn/1 tandem of two FCFS single servers.
 *
 * Templated port of matlab/src/api/qsys/qsys_tandem_ub_ciucu.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_tandem_ub_ciucu.java. Both
 * stations serve the same hyperexponential law Y, Z ~ sum_i p_i Exp(mu_i), a
 * single phase giving exponential service, and the arrivals are renewal with a
 * light-tailed interarrival time supplied through its Laplace-Stieltjes
 * transform E[e^{-s X}].
 *
 * With theta the positive root of E[e^{theta (Y-X)}] = 1 and
 * alpha = E[X e^{-theta X}], the test function
 *
 *   gamma(u,v) = 1{0<=u<=v} [1 - A e^{-theta u} - (B + C u + D v) e^{-theta v}]
 *
 * satisfies the integral inequality of Theorem 1(b) of the reference once the
 * five sufficient conditions of its Lemma 4 fix
 *
 *   A = 1,  C = theta sum_i p_i/(mu_i-theta) / sum_i p_i mu_i/(mu_i-theta)^2,
 *   D = (-C E[U e^{theta V}]/E[V e^{theta V}]) v 0,  U = Y-X, V = Z-X,
 *   B = C (1/mu_1 - alpha E[e^{theta Z}])    if D = 0,
 *     = (C+D)/(mu_1-theta) - theta/mu_1      if D > 0,
 *
 * with mu_1 the smallest service rate. Corollary 2 then turns gamma into
 *
 *   P(S > x) <= sum_i p_i { e^{-mu_i x}
 *                 + mu_i/(mu_i-theta) (A+B) (e^{-theta x} - e^{-mu_i x})
 *                 + mu_i/(mu_i-theta)^2 (C+D) (((mu_i-theta)x-1) e^{-theta x}
 *                                              + e^{-mu_i x}) }
 *
 * and the corresponding closed form for W when the service is exponential.
 * E[V e^{theta V}] is positive at any stable load, so D is always well defined:
 * h(s) = E[e^{s(Z-X)}] is convex with h(0) = h(theta) = 1, hence h'(theta) > 0.
 *
 * The two exponentials mix a polynomial of degree one in x, which is what lets
 * the bound follow the concave bend of the tail on a linear-log scale where a
 * purely exponential bound cannot. In the M/M/1 -> ./M/1 case the five
 * inequalities hold as equalities, so gamma is the exact joint distribution and
 * both bounds are exact, P(S > x) = (1 + theta x) e^{-theta x}. Away from it the
 * bound stays sharp: against an exact CTMC reference for the Erlang(2)/M/1 ->
 * ./M/1 tandem it is within 2% at P(S>x) = 1e-2 and within 0.6% at 5e-10, with
 * the correct asymptotic slope theta^2/(mu(1-alpha mu)). Accuracy degrades with
 * service variability, to about a factor of two at CV(Y) = 2.
 *
 * Reference: F. Ciucu, S. Mehri, "On the Distribution of Sojourn Times in Tandem
 * Queues", Proc. ACM Meas. Anal. Comput. Syst. 9(2), Article 27, 2025 (ACM
 * SIGMETRICS 2025). Registered in .citations() as 'tandemub'.
 *
 * ARITHMETIC: transcendental. theta is a root of a transcendental equation and
 * the bound is a mix of exponentials, so this instantiates only at backends
 * carrying exp; the bracketing solve itself is field arithmetic and
 * deterministic, so the digits do not depend on a starting point.
 */

#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/rootfind.h"

namespace line {
namespace qsys {

/** Mirrors the struct MATLAB returns from qsys_tandem_ub_ciucu. */
template <class T>
struct TandemUbResult {
    std::vector<T> S;   ///< upper bound on P(S > x), one per threshold, capped at one
    std::vector<T> W;   ///< upper bound on P(W > x), NaN unless the service is exponential
    T theta;            ///< tail decay rate, positive root of E[e^{theta (Y-X)}] = 1
    T alpha;            ///< E[X e^{-theta X}]
    T A;                ///< coefficient A of gamma, fixed by Lemma 4
    T B;                ///< coefficient B of gamma, fixed by Lemma 4
    T C;                ///< coefficient C of gamma, fixed by Lemma 4
    T D;                ///< coefficient D of gamma, fixed by Lemma 4
    std::string analyzer;
};

/**
 * @param x    thresholds at which the tails are bounded, nonnegative
 * @param lst  interarrival transform, s -> E[e^{-s X}] for s >= 0
 * @param p    service phase probabilities, nonnegative and summing to one
 * @param mu   service phase rates, positive
 * @param dlst s -> E[X e^{-s X}], minus the derivative of lst; empty to obtain
 *             it by a Richardson-extrapolated central difference, which costs
 *             four extra transform evaluations and loses roughly four digits
 */
template <class T>
TandemUbResult<T> qsys_tandem_ub_ciucu(const std::vector<T>& x,
                                       const std::function<T(const T&)>& lst,
                                       const std::vector<T>& p, const std::vector<T>& mu,
                                       const std::function<T(const T&)>& dlst =
                                           std::function<T(const T&)>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_tandem_ub_ciucu requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    if (x.empty()) throw InputError("qsys_tandem_ub_ciucu: x must hold at least one threshold");
    if (!lst) throw InputError("qsys_tandem_ub_ciucu: lst must be supplied");
    if (p.empty() || p.size() != mu.size())
        throw InputError("qsys_tandem_ub_ciucu: p and mu must have the same number of phases");
    T psum = zero;
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (p[i] < zero)
            throw InputError("qsys_tandem_ub_ciucu: the phase probabilities p must be nonnegative");
        if (mu[i] <= zero)
            throw InputError("qsys_tandem_ub_ciucu: the service rates mu must be positive");
        psum = psum + p[i];
    }
    const T tolp = num_traits<T>::from_rational(1, 10000000000L);
    if (psum - one > tolp || one - psum > tolp)
        throw InputError("qsys_tandem_ub_ciucu: the phase probabilities p must sum to one");
    for (std::size_t k = 0; k < x.size(); ++k)
        if (x[k] < zero) throw InputError("qsys_tandem_ub_ciucu: the thresholds x must be nonnegative");

    T mu1 = mu[0];
    for (std::size_t i = 1; i < mu.size(); ++i)
        if (mu[i] < mu1) mu1 = mu[i];

    // E[e^{t Y}] of the hyperexponential service law, for t below every rate.
    const auto mgfY = [&p, &mu](const T& t) {
        T acc = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < p.size(); ++i) acc = acc + p[i] * mu[i] / (mu[i] - t);
        return acc;
    };
    const auto residual = [&mgfY, &lst, &one](const T& t) { return T(mgfY(t) * lst(t) - one); };

    // Stability: E[X] > E[Y] is what makes E[e^{t(Y-X)}] - 1 cross zero on (0,mu1).
    const T shrink = one - num_traits<T>::from_rational(1, 1000000000000L);
    const T hi = mu1 * shrink;
    if (residual(hi) <= zero)
        throw InputError("qsys_tandem_ub_ciucu: no positive root of E[e^{theta(Y-X)}]=1 below "
                         "min(mu): the tandem is unstable or the service is not the lighter tail");
    const T ten = num_traits<T>::from_int(10);
    T lo = mu1 * num_traits<T>::from_rational(1, 1000000000000L);
    const T lomin = mu1 * num_traits<T>::from_rational(1, 1000000000000L) /
                    num_traits<T>::from_int(10000);
    while (residual(lo) >= zero && lo > lomin) lo = lo / ten;  // walk below the root at zero
    if (residual(lo) >= zero)
        // E[e^{t(Y-X)}]-1 is convex and vanishes at t=0, so it stays positive on the
        // whole of (0,mu1) exactly when its slope E[Y]-E[X] there is nonnegative.
        throw InputError("qsys_tandem_ub_ciucu: the tandem is unstable, E[X] <= E[Y]: "
                         "theta = 0 is the only root of E[e^{theta(Y-X)}]=1");
    const T tol = num_traits<T>::from_rational(1, 1000000000L) /
                  num_traits<T>::from_int(100000);
    const RootResult<T> rr = root_brent<T, decltype(residual)>(residual, lo, hi, tol);
    const T theta = rr.root;

    T alpha;
    if (dlst) {
        alpha = dlst(theta);
    } else {
        // Richardson-extrapolated central difference of -lst at theta.
        T h = num_traits<T>::from_rational(1, 1000) * (one + theta);
        if (h > theta) h = theta / two;
        const T d1 = (lst(T(theta - h)) - lst(T(theta + h))) / (two * h);
        const T d2 = (lst(T(theta - h / two)) - lst(T(theta + h / two))) / h;
        alpha = (num_traits<T>::from_int(4) * d2 - d1) / num_traits<T>::from_int(3);
    }

    const T EexpZ = mgfY(theta);  // E[e^{theta Z}]
    T EZexp = zero;               // E[Z e^{theta Z}]
    T EY = zero;
    T sumPOverMuMinusTheta = zero;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const T dm = mu[i] - theta;
        EZexp = EZexp + p[i] * mu[i] / (dm * dm);
        EY = EY + p[i] / mu[i];
        sumPOverMuMinusTheta = sumPOverMuMinusTheta + p[i] / dm;
    }
    const T A = one;
    const T C = theta * sumPOverMuMinusTheta / EZexp;
    const T EUeV = EY - alpha * EexpZ;              // E[U e^{theta V}]
    const T EVeV = EZexp / EexpZ - alpha * EexpZ;   // E[V e^{theta V}] > 0
    T D = -C * EUeV / EVeV;
    if (!(D > zero)) D = zero;
    const T B = D > zero ? T((C + D) / (mu1 - theta) - theta / mu1)
                         : T(C * (one / mu1 - alpha * EexpZ));

    TandemUbResult<T> r;
    r.theta = theta;
    r.alpha = alpha;
    r.A = A;
    r.B = B;
    r.C = C;
    r.D = D;
    r.analyzer = "qsys_tandem_ub_ciucu";
    r.S.reserve(x.size());
    r.W.reserve(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) {
        T acc = zero;
        const T et = detail::num_exp(T(-theta * x[k]));
        for (std::size_t i = 0; i < p.size(); ++i) {
            const T m = mu[i];
            const T dm = m - theta;
            const T em = detail::num_exp(T(-m * x[k]));
            acc = acc + p[i] * (em + m / dm * (A + B) * (et - em) +
                                m / (dm * dm) * (C + D) * ((dm * x[k] - one) * et + em));
        }
        r.S.push_back(acc < one ? acc : one);
    }
    if (p.size() == 1) {
        const T beta = lst(mu1);  // E[e^{-mu X}]
        for (std::size_t k = 0; k < x.size(); ++k) {
            const T et = detail::num_exp(T(-theta * x[k]));
            T val;
            if (D == zero) {
                val = (one - two * theta * theta / (mu1 * (mu1 + theta)) +
                       theta * (mu1 - theta) / (mu1 + theta) * x[k]) * et +
                      beta * (theta * mu1 * alpha / (two * (mu1 - theta)) - theta / (two * mu1)) *
                          detail::num_exp(T(-mu1 * x[k]));
            } else {
                val = (one - two * theta / mu1 +
                       two * theta * theta * (two - alpha * mu1) /
                           ((mu1 + theta) * (mu1 + theta) * (one - alpha * mu1)) +
                       theta * theta * (mu1 - theta) /
                           (mu1 * (mu1 + theta) * (one - alpha * mu1)) * x[k]) * et;
            }
            r.W.push_back(val < one ? val : one);
        }
    } else {
        // The W form is Exp-service only; MATLAB returns NaN there and so does this.
        for (std::size_t k = 0; k < x.size(); ++k)
            r.W.push_back(num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN()));
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_TANDEM_UB_CIUCU_H
