/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_QUEST_HEURISTIC_CI_H
#define LINE_API_SIM_SIM_QUEST_HEURISTIC_CI_H

/**
 * Fallback interval used when a QUEST stage test fails.
 *
 * Port of matlab/src/api/sim/sim_quest_heuristic_ci.m. Three intervals are formed
 * and the smallest interval containing all of them is returned, which is the
 * article's prescription:
 *
 *   Two symmetric intervals of half-width
 *     h = max(t_{1-alpha/2,K} sqrt(Ap/nstar), t_{1-alpha/2,K-1} sqrt(Np/nstar)),
 *   one about the full-sample quantile and one about the average of the batched
 *   quantile estimators. Taking the WIDER of the two variance components is
 *   deliberately conservative: a stage test has just failed, so neither
 *   component can be trusted.
 *
 *   Willink's asymmetric interval, which corrects the batched quantile
 *   estimators for skewness through the cube-root transform
 *   G(zeta) = ([1+6 gamma(zeta-gamma)]^(1/3)-1)/(2 gamma) with
 *   gamma = skewness/(6 sqrt(K)), evaluated at both t-quantiles so the two arms
 *   differ.
 *
 * useAutocorr additionally scales the asymmetric arms by
 * max(sqrt((1+phi1)/(1-phi1)), 1), with phi1 the lag-1 autocorrelation of the
 * batch quantiles. It is true for sim_fquest, where they come from one sample
 * path and can stay correlated, and false for sim_firquest, where they come
 * from independent replications and the article drops the correction.
 *
 * Np ARRIVES AS NaN when the caller had a single batch per path. The reference
 * relies on MATLAB's max() omitting NaN, so the interval then rests on the Ap
 * component alone; the max here omits NaN explicitly for the same reason, since
 * a C++ comparison against NaN is false in both directions and the naive form
 * would propagate the NaN into the delivered endpoints.
 *
 * Reference: R. Willink, "A Confidence Interval and Test for the Mean of an
 * Asymmetric Distribution", Commun. Statist. Theory Methods 34, 2005;
 * A. Lolos et al., Proc. Winter Simulation Conference, 2023, step 10, and
 * Proc. Winter Simulation Conference, 2025, equations 8 to 10.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/sim/sim_dist.h"
#include "line/api/sim/sim_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

/** A confidence interval, asymmetric about the point estimate in general. */
template <class T>
struct QuestInterval {
    T lower;
    T upper;
};

namespace detail {

/** MATLAB max() of two scalars, i.e. the larger one with NaN omitted. */
template <class T>
inline T max_omitnan(const T& a, const T& b) {
    if (num_isnan(a)) return b;
    if (num_isnan(b)) return a;
    return a < b ? b : a;
}

/** Willink's skewness-adjustment transform, the identity for tiny skewness. */
template <class T>
inline T willink_g(const T& zeta, const T& gamma) {
    const T lim = num_traits<T>::from_double(0.001);
    if (num_abs(gamma) <= lim) return zeta;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T six = num_traits<T>::from_int(6);
    const T third = num_traits<T>::from_rational(1, 3);
    const T arg = T(one + six * gamma * T(zeta - gamma));
    // the cube root is taken on the reals, the argument may turn negative for a
    // strongly skewed and small batch sample
    const T root = arg >= zero ? num_pow(arg, third) : T(-num_pow(T(-arg), third));
    return T(T(root - one) / T(num_traits<T>::from_int(2) * gamma));
}

}  // namespace detail

/**
 * @param bqe         the K batched quantile estimators, pooled over replications
 *                    for sim_firquest
 * @param centre      the full-sample empirical quantile
 * @param Ap          STS area variance-parameter estimator
 * @param Np          NBQ variance-parameter estimator, may be NaN
 * @param nstar       number of observations Ap and Np were computed from
 * @param alpha       nominal non-coverage
 * @param useAutocorr apply the lag-1 correction to the asymmetric arms
 */
template <class T>
QuestInterval<T> sim_quest_heuristic_ci(const std::vector<T>& bqe, const T& centre, const T& Ap,
                                        const T& Np, std::size_t nstar, double alpha,
                                        bool useAutocorr) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_quest_heuristic_ci: t quantiles and square roots, so exact arithmetic is "
                  "refused");
    const std::size_t K = bqe.size();
    if (K < 3)
        throw InputError("sim_quest_heuristic_ci: the heuristic interval needs at least 3 batch "
                         "quantiles");
    if (nstar < 1)
        throw InputError("sim_quest_heuristic_ci: nstar must be positive");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("sim_quest_heuristic_ci: alpha must be a real scalar in (0,1)");

    const double Kd = static_cast<double>(K);
    const T nst = num_traits<T>::from_int(static_cast<long>(nstar));
    const T tK = num_traits<T>::from_double(sim_tinv(1.0 - alpha / 2.0, Kd));
    const T tKm1 = num_traits<T>::from_double(sim_tinv(1.0 - alpha / 2.0, Kd - 1.0));

    const T half = detail::max_omitnan(T(tK * detail::num_sqrt(T(Ap / nst))),
                                       T(tKm1 * detail::num_sqrt(T(Np / nst))));

    T sum = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i) sum += bqe[i];
    const T bqeBar = T(sum / num_traits<T>::from_int(static_cast<long>(K)));

    const T Km1 = num_traits<T>::from_int(static_cast<long>(K - 1));
    T s2 = num_traits<T>::from_int(0), s2t = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i) {
        const T d = T(bqe[i] - bqeBar);
        const T dt = T(bqe[i] - centre);
        s2 += T(d * d);
        s2t += T(dt * dt);
    }
    const T S2 = T(s2 / Km1);
    const T S2tilde = T(s2t / Km1);

    QuestInterval<T> ci;
    ci.lower = std::min(T(centre - half), T(bqeBar - half));
    ci.upper = std::max(T(centre + half), T(bqeBar + half));
    if (!(S2 > num_traits<T>::from_int(0))) return ci;

    const T sdev = detail::num_sqrt(S2);
    T cube = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i) {
        const T z = T(T(bqe[i] - bqeBar) / sdev);
        cube += T(z * z * z);
    }
    const T skew = T(num_traits<T>::from_double(Kd / ((Kd - 1.0) * (Kd - 2.0))) * cube);
    const T gamma = T(skew / num_traits<T>::from_double(6.0 * std::sqrt(Kd)));

    T varphi = num_traits<T>::from_int(1);
    if (useAutocorr) {
        T lag = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i + 1 < K; ++i)
            lag += T(T(bqe[i] - bqeBar) * T(bqe[i + 1] - bqeBar));
        const T phi1 = T(lag / T(Km1 * S2));
        const T one = num_traits<T>::from_int(1);
        if (num_abs(phi1) < one) {
            const T v = detail::num_sqrt(T(T(one + phi1) / T(one - phi1)));
            varphi = v < one ? one : v;
        }
    }

    const T tq = num_traits<T>::from_double(sim_tinv(1.0 - alpha / 2.0, Kd - 1.0));
    const T scale = T(varphi * detail::num_sqrt(T(S2tilde / num_traits<T>::from_int(
                                                              static_cast<long>(K)))));
    const T G1 = T(detail::willink_g(tq, gamma) * scale);
    const T G2 = T(detail::willink_g(T(-tq), gamma) * scale);

    const T armA = T(centre - G1);
    const T armB = T(centre - G2);
    ci.lower = std::min(ci.lower, std::min(armA, armB));
    ci.upper = std::max(ci.upper, std::max(armA, armB));
    return ci;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_QUEST_HEURISTIC_CI_H
