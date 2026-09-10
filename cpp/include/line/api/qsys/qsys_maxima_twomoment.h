/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MAXIMA_TWOMOMENT_H
#define LINE_API_QSYS_MAXIMA_TWOMOMENT_H

/**
 * Two-moment approximation for the maximum of n iid non-negative variables.
 *
 * Templated port of matlab/src/api/qsys/qsys_maxima_twomoment.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_maxima_twomoment.java.
 *
 * THE SHAPE OF THE ANSWER. For a law with an exponential-like tail the maximum
 * of n samples grows like c~^2 (log n + ...): doubling n ADDS a constant, it
 * does not scale the answer. The two moments buy the SLOPE and an offset:
 *
 *   x_n(q) = c~^2 [log(n eta) - log log(1/q)]                       (1.9)
 *   E[M_n] = c~^2 [log(n eta) + gamma]
 *   cs2 >= 1:  c~^2 = cs2,       eta = (cs2+1)/(2 cs2^2)            (1.11)-(1.12)
 *   cs2 <  1:  c~^2 = sqrt(cs2), eta = exp((1-sqrt(cs2))/sqrt(cs2))
 *
 * WHEN NOT TO USE IT: n must pass n* ~ cs2/q (4.22), because with a highly
 * variable law only about n p of the samples can contend for the maximum.
 * Measured against exact maxima the closed form is within a few percent for
 * n >= 100 at cs2 = 4 and 16, and useless at n = 10 for cs2 = 16 -- exactly what
 * n* predicts.
 *
 * AND WHEN TWO MOMENTS ARE NOT ENOUGH: below cs2 = 1 the maximum is genuinely
 * family-dependent. An Erlang and a shifted exponential with the same two
 * moments have maxima differing by tens of percent, diverging as n grows,
 * because their tails decay at different rates.
 *
 * ARITHMETIC. Logarithms throughout: transcendental only.
 *
 * Reference: C. Crow, D. Goldberg, W. Whitt (2007). Two-moment approximations
 * for maxima. Operations Research 55(3), 532-548.
 */

#include <cmath>
#include <cstddef>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

/** Two-moment description of a maximum. */
template <class T>
struct MaximaResult {
    T value;             ///< the closed-form mean or quantile
    T slope;             ///< mean * c~^2, the coefficient of log n
    T eta;               ///< the offset inside the logarithm
    T threshold;         ///< n*, below which the form should not be used
    bool reliable;       ///< whether n >= n*
    std::string family;  ///< the representative law used
    T exactFittedValue;  ///< the maximum computed exactly from that representative
    bool hasFitted = false;
};

/**
 * @param n           the number of samples
 * @param mean        the mean of the underlying law
 * @param cs2         its squared coefficient of variation
 * @param q           a quantile level in (0,1); non-positive returns the mean
 * @param exactFitted also compute the maximum exactly from the fitted law
 */
template <class T>
MaximaResult<T> qsys_maxima_twomoment(std::size_t n, const T& mean, const T& cs2,
                                      const T& q = num_traits<T>::from_int(0),
                                      bool exactFitted = true) {
    static_assert(num_traits<T>::has_transcendental, "qsys_maxima_twomoment needs logarithms");
    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (n < 1) throw InputError("qsys_maxima_twomoment: at least one sample is required");
    if (mean <= zero) throw InputError("qsys_maxima_twomoment: the mean must be positive");
    if (cs2 <= zero) throw InputError("qsys_maxima_twomoment: the SCV must be positive");
    if (q < zero || q >= one)
        throw InputError("qsys_maxima_twomoment: the quantile level must lie in [0,1)");

    const T EULER = num_traits<T>::from_double(0.5772156649015329);
    MaximaResult<T> r;
    T ct, eta;
    if (cs2 >= one) {
        ct = cs2;
        eta = (cs2 + one) / (two * cs2 * cs2);
        r.family = "H2";
    } else {
        ct = sqrt(cs2);
        eta = exp((one - sqrt(cs2)) / sqrt(cs2));
        r.family = "shifted exponential";
    }
    const T nT = num_traits<T>::from_int(static_cast<long>(n));
    const bool wantMean = (q <= zero);
    const T inner = wantMean ? T(log(nT * eta) + EULER) : T(log(nT * eta) - log(log(one / q)));
    r.value = mean * ct * inner;
    r.slope = mean * ct;
    r.eta = eta;
    r.threshold = cs2 / (wantMean ? num_traits<T>::from_rational(1, 2) : q);     // eq. (4.22)
    r.reliable = nT >= r.threshold;

    if (exactFitted) {
        // Fit the representative law and compute the maximum exactly from F^n.
        T d = zero, m = mean, p1 = zero, l1 = zero, l2 = zero, hi;
        if (cs2 >= one) {
            p1 = (one + sqrt((cs2 - one) / (cs2 + one))) / two;
            l1 = two * p1 / mean;
            l2 = two * (one - p1) / mean;
            hi = num_traits<T>::from_int(40) * mean * (cs2 > one ? cs2 : one);
        } else {
            d = mean * (one - sqrt(cs2));
            m = mean * sqrt(cs2);
            hi = d + num_traits<T>::from_int(40) * m;
        }
        auto ccdf = [&](const T& t) {
            if (cs2 >= one) return T(p1 * exp(-l1 * t) + (one - p1) * exp(-l2 * t));
            return t <= d ? one : T(exp(-(t - d) / m));
        };
        const std::size_t gn = 200000;
        const T h = hi / num_traits<T>::from_int(static_cast<long>(gn));
        r.hasFitted = true;
        if (wantMean) {
            T acc = zero;
            for (std::size_t i = 0; i <= gn; ++i) {
                const T t = num_traits<T>::from_int(static_cast<long>(i)) * h;
                const T v = one - pow(one - ccdf(t), num_traits<T>::from_int(static_cast<long>(n)));
                acc += (i == 0 || i == gn) ? T(v / two) : v;
            }
            r.exactFittedValue = acc * h;
        } else {
            r.exactFittedValue = hi;
            for (std::size_t i = 0; i <= gn; ++i) {
                const T t = num_traits<T>::from_int(static_cast<long>(i)) * h;
                if (pow(one - ccdf(t), num_traits<T>::from_int(static_cast<long>(n))) >= q) {
                    r.exactFittedValue = t;
                    break;
                }
            }
        }
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MAXIMA_TWOMOMENT_H
