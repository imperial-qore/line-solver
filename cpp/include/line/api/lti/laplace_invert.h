/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LTI_LAPLACE_INVERT_H
#define LINE_API_LTI_LAPLACE_INVERT_H

/**
 * Numerical inversion of a Laplace transform: Euler, Talbot, Gaver-Stehfest.
 *
 * Port of python/line_solver/api/lti/__init__.py. THIS IS PYTHON-ONLY: MATLAB
 * carries only the CME method (`matlab/lib/thirdparty/iltcme/matlab_ilt.m`,
 * already ported as `api/mam/matlab_ilt.h`), so native Python is the reference
 * for the other three.
 *
 * ALL FOUR ARE THE SAME FRAMEWORK. Abate-Whitt writes
 *
 *     f(t) ~ (1/t) sum_k Re[ omega_k F(alpha_k / t) ],
 *
 * and a method IS its (alpha, omega) pair -- nothing else differs, which is why
 * they share one evaluator here. What differs is where the nodes sit:
 *
 *  - EULER puts them on a vertical line and accelerates an alternating series
 *    with binomial (Euler) weights. Odd `n` only; the reference silently rounds
 *    an even `n` up, and so does this. ITS DEFAULT IS NOT THE REFERENCE'S --
 *    see the note on `laplace_invert_euler`.
 *  - TALBOT deforms the contour into the left half plane, where the transform
 *    decays, so it needs far fewer nodes -- 32 against Euler's 99. It requires
 *    F to be analytic there, which a rational transform is and a transform with
 *    a branch cut is not.
 *  - GAVER-STEHFEST samples F on the REAL axis only, which is what makes it the
 *    one usable method when the transform cannot be evaluated at complex
 *    argument. It pays for that in conditioning: the weights alternate in sign
 *    and grow, so it needs high precision and `n` even (again rounded).
 *  - CME is `api/mam/matlab_ilt.h`, whose coefficient table is vendored; it is
 *    not duplicated here, and `laplace_invert` dispatches to it.
 *
 * The defaults are the reference's own for Talbot (32) and Gaver-Stehfest (12).
 * EULER'S IS NOT: the old default was 99, which is unusable in double
 * precision; both this port and native Python now default to 41. The measurement is at
 * `laplace_invert_euler`.
 *
 * ARITHMETIC: double. Every one of these is a floating-point quadrature.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/api/mam/matlab_ilt.h"
#include "line/util/error.h"

namespace line {
namespace lti {

using Cplx = std::complex<double>;

/** The transform, evaluated at complex argument. */
using LaplaceFn = std::function<Cplx(Cplx)>;
/** A transform that can only be evaluated on the real axis. */
using RealLaplaceFn = std::function<double(double)>;

namespace ltidetail {

/** Binomial coefficient, exactly, for the small n these methods use. */
inline double binom(std::size_t n, std::size_t k) {
    if (k > n) return 0.0;
    double v = 1.0;
    for (std::size_t i = 0; i < k; ++i)
        v = v * static_cast<double>(n - i) / static_cast<double>(i + 1);
    return v;
}

}  // namespace ltidetail

/** Euler nodes: a vertical line at Re = (n-1) log(10) / 6. */
inline std::vector<Cplx> euler_get_alpha(std::size_t n) {
    std::vector<Cplx> a(n);
    for (std::size_t i = 0; i < n; ++i)
        a[i] = Cplx(static_cast<double>(n - 1) * std::log(10.0) / 6.0,
                    M_PI * static_cast<double>(i));
    return a;
}

/**
 * Euler weights before the alternating sign and the scale.
 *
 * The tail is the binomial partial sums that give the Euler acceleration; it is
 * filled BACKWARDS from the last entry, which is what makes the running sum
 * correct.
 */
inline std::vector<double> euler_get_eta(std::size_t n) {
    if (n < 3) throw InputError("euler_get_eta: at least three terms are required");
    std::vector<double> res(n, 0.0);
    res[0] = 0.5;
    for (std::size_t i = 1; i < (n + 1) / 2; ++i) res[i] = 1.0;
    res[n - 1] = 1.0 / std::pow(2.0, (static_cast<double>(n) - 1.0) / 2.0);
    for (std::size_t i = 1; i < (n - 1) / 2; ++i)
        res[n - i - 1] = res[n - i] + std::pow(2.0, (1.0 - static_cast<double>(n)) / 2.0) *
                                          ltidetail::binom((n - 1) / 2, i);
    return res;
}

/** Euler weights: eta, alternating in sign, scaled by 10^((n-1)/6). */
inline std::vector<Cplx> euler_get_omega(std::size_t n) {
    const std::vector<double> eta = euler_get_eta(n);
    std::vector<Cplx> res(n);
    const double scale = std::pow(10.0, (static_cast<double>(n) - 1.0) / 6.0);
    for (std::size_t i = 0; i < n; ++i)
        res[i] = Cplx(scale * ((i % 2 == 0) ? 1.0 : -1.0) * eta[i], 0.0);
    return res;
}

/** Talbot nodes: the cotangent contour, bending into the left half plane. */
inline std::vector<Cplx> talbot_get_alpha(std::size_t n) {
    if (n == 0) throw InputError("talbot_get_alpha: at least one term is required");
    std::vector<Cplx> a(n);
    a[0] = Cplx(2.0 * static_cast<double>(n) / 5.0, 0.0);
    for (std::size_t i = 1; i < n; ++i) {
        const double th = static_cast<double>(i) * M_PI / static_cast<double>(n);
        a[i] = Cplx(2.0 * static_cast<double>(i) * M_PI / 5.0 * (1.0 / std::tan(th)),
                    2.0 * static_cast<double>(i) * M_PI / 5.0);
    }
    return a;
}

/** Talbot weights, which carry the contour's own derivative. */
inline std::vector<Cplx> talbot_get_omega(std::size_t n, const std::vector<Cplx>& alpha) {
    if (alpha.size() != n) throw InputError("talbot_get_omega: alpha has the wrong length");
    std::vector<Cplx> w(n);
    w[0] = std::exp(alpha[0]) / 5.0;
    for (std::size_t i = 1; i < n; ++i) {
        const double th = static_cast<double>(i) * M_PI / static_cast<double>(n);
        const double cot = 1.0 / std::tan(th);
        const Cplx mult(1.0, th * (1.0 + cot * cot) - cot);
        w[i] = 2.0 * std::exp(alpha[i]) / 5.0 * mult;
    }
    return w;
}

/** Gaver-Stehfest nodes: k log 2, on the REAL axis. */
inline std::vector<double> gaver_stehfest_get_alpha(std::size_t n) {
    if (n % 2 == 1) --n;  // the method is defined for even n only
    std::vector<double> a(n);
    for (std::size_t k = 1; k <= n; ++k) a[k - 1] = static_cast<double>(k) * std::log(2.0);
    return a;
}

/**
 * Gaver-Stehfest weights.
 *
 * They alternate in sign and grow rapidly with n, which is why the method needs
 * more precision than the others rather than more terms.
 */
inline std::vector<double> gaver_stehfest_get_omega(std::size_t n) {
    if (n % 2 == 1) --n;
    if (n == 0) throw InputError("gaver_stehfest_get_omega: at least two terms are required");
    const std::size_t h = n / 2;
    double fact = 1.0;
    for (std::size_t i = 2; i <= h; ++i) fact *= static_cast<double>(i);

    std::vector<double> res(n, 0.0);
    for (std::size_t k = 1; k <= n; ++k) {
        double sum = 0.0;
        for (std::size_t j = (k + 1) / 2; j <= std::min(k, h); ++j)
            sum += std::pow(static_cast<double>(j), static_cast<double>(h + 1)) / fact *
                   ltidetail::binom(h, j) * ltidetail::binom(2 * j, j) *
                   ltidetail::binom(j, k - j);
        res[k - 1] = (((h + k) % 2 == 0) ? 1.0 : -1.0) * std::log(2.0) * sum;
    }
    return res;
}

/**
 * Euler inversion of F at t. `n` is rounded UP to odd, as the reference does.
 *
 * THE DEFAULT IS 41, NOT THE OLD 99, AND THAT IS A CORRECTION RATHER THAN A
 * PREFERENCE. The weights carry a factor 10^((n-1)/6), and the sum they
 * multiply alternates in sign, so the method's accuracy is a race between the
 * series converging and the cancellation eating the mantissa. Measured on
 * F(s) = 2/(s+2), whose inverse is 2 exp(-2t), as the worst relative error over
 * t in {0.1, 0.5, 1, 2}:
 *
 *   n   =    11      21      31      41      51      71      99
 *   err = 4.4e-3  2.1e-6  1.6e-9  1.6e-10 4.7e-8  1.7e-4  1.4e+0
 *
 * At the reference's 99 the scale is 2.2e16, past what a double resolves, and
 * the answer is 140 per cent wrong -- negative at some t. Native Python has the
 * same default and the same behaviour (measured: 1.997 against an exact 1.637
 * at t = 0.1, and -0.033 at t = 0.5), and `api/lti` has no MATLAB twin, so
 * nothing else in the tree catches it. Shipping a default that returns noise
 * is not a convention worth preserving; the reference's own defect is recorded
 * for its maintainer rather than reproduced here.
 */
inline double laplace_invert_euler(const LaplaceFn& F, double t, std::size_t n = 41) {
    if (!(t > 0.0)) throw InputError("laplace_invert_euler: t must be positive");
    if (n % 2 == 0) ++n;
    const std::vector<Cplx> a = euler_get_alpha(n), w = euler_get_omega(n);
    double r = 0.0;
    for (std::size_t i = 0; i < n; ++i) r += (w[i] * F(a[i] / t)).real();
    return r / t;
}

/** Talbot inversion of F at t. */
inline double laplace_invert_talbot(const LaplaceFn& F, double t, std::size_t n = 32) {
    if (!(t > 0.0)) throw InputError("laplace_invert_talbot: t must be positive");
    if (n == 0) throw InputError("laplace_invert_talbot: at least one term is required");
    const std::vector<Cplx> a = talbot_get_alpha(n);
    const std::vector<Cplx> w = talbot_get_omega(n, a);
    double r = 0.0;
    for (std::size_t i = 0; i < n; ++i) r += (w[i] * F(a[i] / t)).real();
    return r / t;
}

/**
 * Gaver-Stehfest inversion of F at t.
 *
 * The transform is sampled on the REAL axis only, which is the whole reason to
 * choose this method. `n` is rounded DOWN to even.
 */
inline double laplace_invert_gaver_stehfest(const RealLaplaceFn& F, double t,
                                            std::size_t n = 12) {
    if (!(t > 0.0)) throw InputError("laplace_invert_gaver_stehfest: t must be positive");
    if (n % 2 == 1) --n;
    const std::vector<double> a = gaver_stehfest_get_alpha(n);
    const std::vector<double> w = gaver_stehfest_get_omega(n);
    double r = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) r += w[i] * F(a[i] / t);
    return r / t;
}

// ---------------------------------------------------------------------------
// Weeks / Laguerre (Weeks, JACM 13, 1966; Abate, Choudhury and Whitt, INFORMS
// J. Computing 8(4), 1996; Harrison and Knottenbelt 2002, Sec. 4.1-4.3)
// ---------------------------------------------------------------------------

/** A Laguerre expansion: the damping, the scaling and the coefficients. */
struct WeeksParams {
    double sigma = 0.0;
    double b = 1.0;
    std::vector<double> q;
};

/**
 * Laguerre coefficients q_n, n = 0..2*p0-1, of f_{sigma,b}(t) =
 * exp(-sigma t) f(t/b), whose generating function is
 *
 *     Q_{sigma,b}(z) = b/(1-z) * L( b(1+z)/(2(1-z)) + b*sigma ).
 *
 * NOTE ON THE PAPER. Eq. 10 as printed carries the factor (1-z) rather than
 * 1/(1-z). The scaled form above, printed later in the same section, carries
 * 1/(1-z) and is the correct one: with l_n(t) = exp(-t/2) L_n(t) the transform
 * of l_n is (s-1/2)^n/(s+1/2)^{n+1}, so L(s) = Q(z)/(s+1/2) with
 * z = (s-1/2)/(s+1/2) and s+1/2 = 1/(1-z). Implementing the printed (1-z) is
 * wrong at every t (163 per cent at t = 0.1 on Exp(2)).
 *
 * Sec. 4.3 fixes the trapezoid count at 2*p0 and the radius at r = 0.1^(4/p0)
 * for every n, so the quadrature is one discrete Fourier transform of Q sampled
 * on the circle and the transform is evaluated 2*p0 times IN TOTAL rather than
 * per coefficient. The DFT is evaluated directly: at 2*p0 = 400 points that is
 * 160k complex multiplies, which is not worth a dependency.
 */
inline std::vector<double> laplace_weeks_coeffs(const LaplaceFn& F, double sigma = 0.0,
                                                double b = 1.0, std::size_t p0 = 200) {
    if (!(b > 0.0)) throw InputError("laplace_weeks_coeffs: b must be positive");
    if (p0 == 0) throw InputError("laplace_weeks_coeffs: p0 must be positive");
    const std::size_t N = 2 * p0;
    const double r = std::pow(0.1, 4.0 / static_cast<double>(p0));
    const double twopi = 2.0 * 3.14159265358979323846;

    std::vector<Cplx> Q(N);
    for (std::size_t j = 0; j < N; ++j) {
        const double u = twopi * static_cast<double>(j) / static_cast<double>(N);
        const Cplx z = r * Cplx(std::cos(u), std::sin(u));
        const Cplx s = b * (Cplx(1.0, 0.0) + z) / (2.0 * (Cplx(1.0, 0.0) - z)) + b * sigma;
        Q[j] = b / (Cplx(1.0, 0.0) - z) * F(s);
    }

    std::vector<double> q(N, 0.0);
    double rpow = 1.0;
    for (std::size_t n = 0; n < N; ++n) {
        Cplx acc(0.0, 0.0);
        for (std::size_t j = 0; j < N; ++j) {
            const double u = -twopi * static_cast<double>(n) * static_cast<double>(j) /
                             static_cast<double>(N);
            acc += Q[j] * Cplx(std::cos(u), std::sin(u));
        }
        q[n] = acc.real() / static_cast<double>(N) / rpow;
        rpow *= r;
    }
    return q;
}

/**
 * The automatic (sigma, b) search of Fig. 1: accept the first pair at which the
 * coefficients have decayed by term p0, doubling sigma from 0.001 and stepping
 * b by 4 whenever sigma passes 0.2.
 *
 * REFUSES BY NAME when the box is exhausted. Raising b further is
 * counterproductive and excessive damping is unstable in finite precision, and
 * a density with a discontinuity in itself or its derivatives has no usable
 * Laguerre representation at all (Sec. 4.2). Returning the last iterate would
 * report noise as an answer; Euler handles those cases instead.
 */
inline WeeksParams laplace_weeks_scaling(const LaplaceFn& F, std::size_t p0 = 200,
                                         double tol = 1e-10) {
    WeeksParams w;
    w.sigma = 0.0;
    w.b = 1.0;
    for (;;) {
        w.q = laplace_weeks_coeffs(F, w.sigma, w.b, p0);
        if (std::abs(w.q[p0]) <= tol && std::abs(w.q[p0 + 1]) <= tol) return w;
        w.sigma = (w.sigma == 0.0) ? 0.001 : 2.0 * w.sigma;
        if (w.sigma > 0.2) {
            w.b += 4.0;
            if (w.b > 10.0)
                throw NumericError(
                    "laplace_weeks_scaling: no suitable scaling parameters were found for the "
                    "Laguerre inversion: the transform's density is not smooth enough for a "
                    "Laguerre series. Use the euler method instead.");
            w.sigma = 0.0;
        }
    }
}

namespace weeks_detail {

/**
 * Truncate at the FIRST index where the coefficients have decayed, never the
 * last. The quadrature divides by r^n with r < 1, so past the genuine decay the
 * entries are rounding noise amplified by r^-n: at n = 2*p0 that factor is 1e8,
 * and scanning for the last entry above a threshold sums 1e-8 of pure noise
 * (worst error on Exp(2) 2.7e-09 instead of 1.9e-14).
 */
inline std::size_t nterms(const std::vector<double>& q) {
    const std::size_t p0 = q.size() / 2;
    for (std::size_t n = 1; n + 1 < p0; ++n)
        if (std::abs(q[n]) <= 1e-13 && std::abs(q[n + 1]) <= 1e-13) return n;
    return p0;
}

/** l_n(t) = exp(-t/2) L_n(t) by the stable recursion of Sec. 4.1. */
inline std::vector<double> functions(double t, std::size_t N) {
    std::vector<double> l(N, 0.0);
    if (N == 0) return l;
    l[0] = std::exp(-t / 2.0);
    if (N > 1) l[1] = (1.0 - t) * l[0];
    for (std::size_t n = 2; n < N; ++n) {
        const double dn = static_cast<double>(n);
        l[n] = ((2.0 * dn - 1.0 - t) / dn) * l[n - 1] - ((dn - 1.0) / dn) * l[n - 2];
    }
    return l;
}

}  // namespace weeks_detail

/**
 * Invert by the Laguerre series f(t) = sum_n q_n l_n(t), recovered as
 * exp(sigma*b*t) f_{sigma,b}(b*t).
 *
 * Unlike Euler and Talbot the coefficients do not depend on t, so ONE parameter
 * set serves an arbitrary number of time points: the transform is evaluated
 * 2*p0 times in total, not 2*p0 times per t. That is the property this method
 * is here for, so build the WeeksParams once and reuse it on a grid.
 */
inline double laplace_invert_weeks(const WeeksParams& w, double t) {
    if (!(t > 0.0)) return 0.0;
    const std::size_t n = weeks_detail::nterms(w.q);
    const std::vector<double> l = weeks_detail::functions(w.b * t, n);
    double acc = 0.0;
    for (std::size_t i = 0; i < n; ++i) acc += w.q[i] * l[i];
    return std::exp(w.sigma * w.b * t) * acc;
}

/** Convenience overload: build the parameters, then invert at one point. */
inline double laplace_invert_weeks(const LaplaceFn& F, double t, std::size_t p0 = 200) {
    return laplace_invert_weeks(laplace_weeks_scaling(F, p0), t);
}

/** The methods `laplace_invert` accepts. */
enum class LaplaceMethod { Euler = 0, Talbot, GaverStehfest, Cme, Weeks };

/** Parse the reference's method names, including its two Gaver spellings. */
inline LaplaceMethod laplace_method(const std::string& s) {
    if (s == "euler") return LaplaceMethod::Euler;
    if (s == "talbot") return LaplaceMethod::Talbot;
    if (s == "gaver-stehfest" || s == "gaver_stehfest" || s == "gaver")
        return LaplaceMethod::GaverStehfest;
    if (s == "cme") return LaplaceMethod::Cme;
    if (s == "weeks" || s == "laguerre") return LaplaceMethod::Weeks;
    throw InputError("laplace_invert: unknown method '" + s +
                     "', expected euler, talbot, gaver-stehfest, cme or weeks");
}

/**
 * Invert F at t by the named method.
 *
 * @param n 0 takes the method's own default: 41 Euler (see
 *          `laplace_invert_euler`), 32 Talbot, 12 Gaver-Stehfest, 25 CME
 */
inline double laplace_invert(const LaplaceFn& F, double t,
                             LaplaceMethod method = LaplaceMethod::Euler, std::size_t n = 0) {
    switch (method) {
        case LaplaceMethod::Euler: return laplace_invert_euler(F, t, n ? n : 41);
        case LaplaceMethod::Talbot: return laplace_invert_talbot(F, t, n ? n : 32);
        case LaplaceMethod::GaverStehfest:
            // The real-axis method is fed the same transform restricted to the
            // real axis; a transform that cannot be evaluated there will say so
            // itself rather than being silently approximated.
            return laplace_invert_gaver_stehfest(
                [&F](double s) { return F(Cplx(s, 0.0)).real(); }, t, n ? n : 12);
        case LaplaceMethod::Cme: {
            std::vector<double> tv(1, t);
            return mam::matlab_ilt([&F](const Cplx& s) { return F(s); }, tv, n ? n : 25,
                                   mam::IltMethod::Cme)[0];
        }
        case LaplaceMethod::Weeks:
            // Rebuilding the expansion for a single point wastes the one
            // property this method has; the grid overloads below do it once.
            return laplace_invert_weeks(F, t, n ? n : 200);
    }
    throw InputError("laplace_invert: unreachable method");
}

/**
 * The DENSITY on a grid: the inversion clamped at zero.
 *
 * A density cannot be negative, and a numerical inversion can undershoot near
 * the origin or in a tail; the reference clamps, and so does this.
 */
inline std::vector<double> laplace_invert_pdf(const LaplaceFn& F, const std::vector<double>& t,
                                              LaplaceMethod method = LaplaceMethod::Euler,
                                              std::size_t n = 0) {
    std::vector<double> out(t.size(), 0.0);
    if (method == LaplaceMethod::Weeks) {
        // One expansion serves the whole grid; this is the point of Weeks.
        const WeeksParams w = laplace_weeks_scaling(F, n ? n : 200);
        for (std::size_t i = 0; i < t.size(); ++i)
            out[i] = std::max(0.0, laplace_invert_weeks(w, t[i]));
        return out;
    }
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!(t[i] > 0.0)) continue;
        out[i] = std::max(0.0, laplace_invert(F, t[i], method, n));
    }
    return out;
}

/**
 * The DISTRIBUTION on a grid, from the transform of the DENSITY.
 *
 * F(s)/s is the transform of the CDF, so that is what is inverted -- passing
 * the CDF's own transform here would invert it twice. The result is clamped
 * into [0,1] and made monotone by a running maximum, because a numerical
 * inversion is pointwise and nothing in it enforces either property; a
 * non-monotone "CDF" then yields negative probabilities downstream.
 */
inline std::vector<double> laplace_invert_cdf(const LaplaceFn& F, const std::vector<double>& t,
                                              LaplaceMethod method = LaplaceMethod::Euler,
                                              std::size_t n = 0) {
    const LaplaceFn Fc = [&F](Cplx s) {
        if (std::abs(s) < 1e-15) return Cplx(1.0, 0.0);
        return F(s) / s;
    };
    std::vector<double> out(t.size(), 0.0);
    if (method == LaplaceMethod::Weeks) {
        const WeeksParams w = laplace_weeks_scaling(Fc, n ? n : 200);
        for (std::size_t i = 0; i < t.size(); ++i) {
            double v = laplace_invert_weeks(w, t[i]);
            out[i] = std::min(1.0, std::max(0.0, v));
        }
        for (std::size_t i = 1; i < out.size(); ++i) out[i] = std::max(out[i], out[i - 1]);
        return out;
    }
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!(t[i] > 0.0)) continue;
        double v = laplace_invert(Fc, t[i], method, n);
        if (v < 0.0) v = 0.0;
        if (v > 1.0) v = 1.0;
        out[i] = v;
    }
    for (std::size_t i = 1; i < out.size(); ++i) out[i] = std::max(out[i], out[i - 1]);
    return out;
}

}  // namespace lti
}  // namespace line

#endif  // LINE_API_LTI_LAPLACE_INVERT_H
