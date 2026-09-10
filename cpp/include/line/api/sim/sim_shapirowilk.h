/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_SHAPIROWILK_H
#define LINE_API_SIM_SIM_SHAPIROWILK_H

/**
 * Shapiro-Wilk test for univariate normality.
 *
 * Port of matlab/src/api/sim/sim_shapirowilk.m, i.e. Royston's AS R94 algorithm,
 * valid for 3 <= n <= 5000. The statistic is
 *   W = (sum_i a_i x_(i))^2 / sum_i (x_i - xbar)^2,
 * where x_(i) are the order statistics and a is the antisymmetric weight vector
 * obtained by correcting the normalized expected normal order statistics
 * m_i = Phi^{-1}((i-3/8)/(n+1/4)) in their two extreme components. Small W means
 * departure from normality, so the test is one-sided in W and the p-value is an
 * upper normal tail after Royston's normalizing transform, which has three
 * branches: n = 3 exact, 4 <= n <= 11, and n >= 12.
 *
 * The weights are computed in double throughout. They are a function of n
 * alone -- Phi^{-1} at fixed plotting positions plus two polynomial corrections
 * with published five-digit coefficients -- so refining them past double would
 * refine a constant that is only known to five digits anyway; the sample itself
 * enters in T.
 *
 * Reference: J. P. Royston, "Approximating the Shapiro-Wilk W-test for
 * Non-normality", Statistics and Computing 2, 1992; J. P. Royston, "Remark
 * AS R94", Applied Statistics 44(4), 1995. W and the p-value agree with
 * scipy.stats.shapiro to 5e-10 and 1.5e-7 respectively over n up to 2000.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include "line/api/sim/sim_dist.h"
#include "line/api/sim/sim_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

/** Outcome of the Shapiro-Wilk normality test. */
template <class T>
struct ShapiroWilkResult {
    T W;                   ///< The Shapiro-Wilk statistic
    double pvalue = 1.0;   ///< p-value, small means normality is rejected
    double zscore = 0.0;   ///< Normalized statistic, NaN when n = 3
    bool reject = false;   ///< True when pvalue < alpha
    std::size_t nobs = 0;  ///< Number of observations n
};

namespace detail {

/** Royston AS R94 antisymmetric weight vector, a(n+1-i) = -a(i), zero based. */
inline std::vector<double> shapirowilk_weights(std::size_t n) {
    std::vector<double> a(n, 0.0);
    if (n == 3) {
        a[0] = -std::sqrt(0.5);
        a[1] = 0.0;
        a[2] = std::sqrt(0.5);
        return a;
    }

    static const double c1[6] = {0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056};
    static const double c2[6] = {0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633};

    std::vector<double> m(n, 0.0);
    double mm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        m[i] = sim_norminv((static_cast<double>(i + 1) - 0.375) /
                           (static_cast<double>(n) + 0.25));
        mm += m[i] * m[i];
    }
    const double u = 1.0 / std::sqrt(static_cast<double>(n));

    // MATLAB's polyval(fliplr(c), u), i.e. the ascending-power evaluation
    double p1 = 0.0, p2 = 0.0, uk = 1.0;
    for (int k = 0; k < 6; ++k) {
        p1 += c1[k] * uk;
        p2 += c2[k] * uk;
        uk *= u;
    }

    a = m;
    const double an = m[n - 1] / std::sqrt(mm) + p1;
    if (n > 5) {
        const double anm1 = m[n - 2] / std::sqrt(mm) + p2;
        const double phi = (mm - 2.0 * m[n - 1] * m[n - 1] - 2.0 * m[n - 2] * m[n - 2]) /
                           (1.0 - 2.0 * an * an - 2.0 * anm1 * anm1);
        for (std::size_t i = 2; i + 2 < n; ++i) a[i] = m[i] / std::sqrt(phi);
        a[n - 1] = an;
        a[n - 2] = anm1;
        a[0] = -an;
        a[1] = -anm1;
    } else {
        const double phi = (mm - 2.0 * m[n - 1] * m[n - 1]) / (1.0 - 2.0 * an * an);
        for (std::size_t i = 1; i + 1 < n; ++i) a[i] = m[i] / std::sqrt(phi);
        a[n - 1] = an;
        a[0] = -an;
    }
    return a;
}

}  // namespace detail

/**
 * @param x     the sample, 3 to 5000 finite observations, order irrelevant
 * @param alpha significance level in (0,1), 0.05 by default
 */
template <class T>
ShapiroWilkResult<T> sim_shapirowilk(const std::vector<T>& x, double alpha = 0.05) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_shapirowilk: the weights and the p-value are transcendental, so exact "
                  "arithmetic is refused");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("sim_shapirowilk: alpha must be a real scalar in (0,1)");

    const std::size_t n = x.size();
    if (n < 3)
        throw InputError("sim_shapirowilk: at least 3 observations are required");
    if (n > 5000)
        throw InputError("sim_shapirowilk: the AS R94 approximation is valid up to n = 5000");
    for (std::size_t i = 0; i < n; ++i)
        if (!detail::num_isfinite(x[i]))
            throw InputError("sim_shapirowilk: the sample must be finite");

    std::vector<T> s(x);
    std::sort(s.begin(), s.end());

    T sum = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) sum += s[i];
    const T mean = sum / num_traits<T>::from_int(static_cast<long>(n));
    T ssd = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        const T d = T(s[i] - mean);
        ssd += T(d * d);
    }
    if (!(ssd > num_traits<T>::from_int(0)))
        throw InputError("sim_shapirowilk: the sample is constant, W is undefined");

    const std::vector<double> a = detail::shapirowilk_weights(n);
    T ax = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) ax += T(num_traits<T>::from_double(a[i]) * s[i]);

    ShapiroWilkResult<T> r;
    r.nobs = n;
    r.W = T(T(ax * ax) / ssd);
    const T one = num_traits<T>::from_int(1);
    if (r.W > one) r.W = one;

    const double Wd = num_traits<T>::to_double(r.W);
    if (n == 3) {
        // exact null distribution, W is supported on [3/4, 1]
        const double pi = boost::math::constants::pi<double>();
        double p = 6.0 / pi * (std::asin(std::sqrt(Wd)) - std::asin(std::sqrt(0.75)));
        r.pvalue = std::min(std::max(p, 0.0), 1.0);
        r.zscore = std::numeric_limits<double>::quiet_NaN();
    } else {
        const double nd = static_cast<double>(n);
        double w, mu, sigma;
        if (n <= 11) {
            const double g = -2.273 + 0.459 * nd;
            w = -std::log(g - std::log(1.0 - Wd));
            mu = 0.5440 - 0.39978 * nd + 0.025054 * nd * nd - 0.0006714 * nd * nd * nd;
            sigma = std::exp(1.3822 - 0.77857 * nd + 0.062767 * nd * nd -
                             0.0020322 * nd * nd * nd);
        } else {
            const double ln = std::log(nd);
            w = std::log(1.0 - Wd);
            mu = -1.5861 - 0.31082 * ln - 0.083751 * ln * ln + 0.0038915 * ln * ln * ln;
            sigma = std::exp(-0.4803 - 0.082676 * ln + 0.0030302 * ln * ln);
        }
        r.zscore = (w - mu) / sigma;
        r.pvalue = 1.0 - sim_normcdf(r.zscore);
    }
    r.reject = r.pvalue < alpha;
    return r;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_SHAPIROWILK_H
