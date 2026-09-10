/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_DISTRIBUTION_H
#define LINE_LANG_DISTRIBUTION_H

/**
 * What `refreshProcessRepresentations` and `refreshLST` compute FROM a
 * distribution: the (D0,D1) pair that reaches `sn.proc`, the arrival-phase
 * vector `sn.pie`, and the Laplace-Stieltjes transform `sn.lst`.
 *
 * These are free functions rather than members of `Distrib` because they need
 * the api layer -- the stationary vector of a MAP is a linear solve
 * (`api/mam/map_moment.h`), the Erlang approximation of a non-Markovian
 * distribution is `map_erlang`, and the Replayer's is an APH fit -- and the
 * model layer would otherwise depend on the api layer wholesale.
 *
 * WHERE THE REFERENCE IS BUG-FOR-BUG REPRODUCED, deliberately. The Weibull and
 * Lognormal transforms in MATLAB are 1000-point RIGHT-ENDPOINT Riemann sums
 * over a truncated interval, not converged quadrature: they are biased low by
 * the tail they drop and by the O(dx) rule. Their values enter M/G/1 waiting
 * times, so replacing them with an accurate integral would move numbers this
 * port is supposed to match. The Pareto transform, by contrast, IS converged in
 * the reference (adaptive Gauss-Kronrod at RelTol 1e-12 over the substitution
 * u = k/x), and is reproduced as such with the ported `num_integral`.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/api/mam/aph_fit.h"
#include "line/api/mam/libqbd_taylor.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/lang/lang_types.h"
#include "line/api/infer/infer_nhpp_ks.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/matrix.h"

namespace line {
namespace lang {

/**
 * The number of Erlang phases `convertToMAP` picks for a non-Markovian
 * distribution: 20 when the SCV is below CoarseTol (a Det, or near one), and
 * otherwise ceil(1/SCV) capped at 100.
 *
 * The cap is what makes the approximation one-sided: a Pareto of SCV 64 gets a
 * single phase (an exponential), so the approximation matches the mean and NOT
 * the SCV whenever the SCV exceeds 1. That is the reference's behaviour and the
 * reason `sn.scv` is read from the distribution rather than from `sn.proc`.
 */
inline unsigned convert_to_map_phases(double scv) {
    if (scv < GlobalConstants::CoarseTol) return 20;
    const double n = std::ceil(1.0 / scv);
    const double capped = n < 1.0 ? 1.0 : (n > 100.0 ? 100.0 : n);
    return static_cast<unsigned>(capped);
}

/**
 * The (D0,D1) pair that reaches `sn.proc`.
 *
 * A type that carries its own representation returns it unchanged. Det,
 * Uniform, Pareto, Gamma, Weibull and Lognormal are replaced by the Erlang
 * approximation of `convertToMAP`; a Replayer is fitted by `aph_fit`, which is
 * what MATLAB's `Replayer.fitAPH` does before taking getProcess.
 */
/**
 * The refusal every lowering of a `Prior` shares.
 *
 * A Prior is a set of models, not one law, so there is no (D0,D1), no transform
 * and no moment of it that a solver could integrate: substituting any single
 * alternative would answer for a model the caller did not describe, and
 * collapsing the set to its mixture would answer for a model nobody described.
 * SolverUQ is the one consumer, and it replaces the Prior before the design
 * point is solved.
 */
inline void reject_prior(const char* who) {
    throw UnsupportedError(std::string(who) +
                           ": the distribution is a Prior, which is a weighted set of alternative "
                           "MODELS rather than one law; solve the model with SolverUQ, which "
                           "replaces each Prior by one alternative per design point");
}

/**
 * The first two moments of a DISCRETE-time MAP, from its own law.
 *
 * With alpha the arrival-epoch stationary vector -- `dmap_pie`, the stationary
 * vector of (I - D0)^-1 D1 -- the interarrival count has
 * P(N = k) = alpha D0^(k-1) D1 e, so
 *
 *   E[N]      = alpha (I - D0)^-1 e            (MATLAB `DMAP.getMean`)
 *   E[N(N-1)] = 2 alpha D0 (I - D0)^-2 e
 *
 * THE SCV IS NOT MATLAB'S INHERITED ONE. `DMAP` declares no getSCV, so it falls
 * through to `Markovian.getSCV` = map_scv({D0,D1}), a CONTINUOUS-time formula
 * that reads D0 + D1 as a generator; for a DMAP that matrix is stochastic, so
 * the stationary solve behind it is singular and the number it returns is not
 * the SCV of anything. Reproducing it would propagate an undefined value into
 * every AMVA path, so the discrete second moment is computed here and the
 * reference defect is recorded in BUGS.md.
 */
template <class T>
void dmap_refresh_moments(Distrib<T>& d) {
    const std::size_t n = d.D0.rows();
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    Matrix<T> ImD0(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) ImD0(i, j) = T((i == j ? one : T(0)) - d.D0(i, j));
    // alpha: the stationary vector of P = (I - D0)^-1 D1, formed column by
    // column through one factorization rather than by inverting ImD0.
    Matrix<T> LU = ImD0;
    const std::vector<std::size_t> piv = lu_factor(LU);
    Matrix<T> P(n, n);
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<T> col(n);
        for (std::size_t i = 0; i < n; ++i) col[i] = d.D1(i, j);
        lu_solve(LU, piv, col);
        for (std::size_t i = 0; i < n; ++i) P(i, j) = col[i];
    }
    const std::vector<T> alpha = mc::dtmc_solve(P);
    std::vector<T> y(n, one);
    lu_solve(LU, piv, y);  // (I - D0)^-1 e
    std::vector<T> z = y;
    lu_solve(LU, piv, z);  // (I - D0)^-2 e
    T m1 = num_traits<T>::from_int(0), fac2 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        m1 += T(alpha[i] * y[i]);
        T dz = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) dz += T(d.D0(i, j) * z[j]);
        fac2 += T(alpha[i] * dz);
    }
    fac2 = T(two * fac2);
    const T m2 = T(fac2 + m1);
    d.mean = m1;
    d.scv = T((m2 - m1 * m1) / (m1 * m1));
}

template <class T>
mam::Map<T> dist_to_map(const Distrib<T>& d) {
    if (d.is_prior()) reject_prior("dist_to_map");
    if (d.disabled) throw InputError("dist_to_map: the distribution is disabled");
    if (d.type == ProcessType::NORMAL)
        throw UnsupportedError(
            "dist_to_map: a Normal puts mass below zero, so it is not the law of any duration and "
            "has no Markovian representation. The default arm of this function would hand back the "
            "Erlang fit of its mean, which is a positive law with the same mean and nothing else "
            "in common; refusing instead. A Normal reaches this port only as the parameter density "
            "of a continuous Prior");
    if (d.has_map()) {
        mam::Map<T> m;
        m.D0 = d.D0;
        m.D1 = d.D1;
        return m;
    }
    if (d.type == ProcessType::REPLAYER) {
        // MATLAB fits an APH to the trace's first three moments. The fit is
        // Bobbio-Horvath-Telek, which needs roots and exponentials, so it is
        // gated: without the guard the static_assert inside aph_fit fires for
        // EVERY exact-arithmetic caller of dist_to_map, trace or not, because
        // the branch is instantiated whether or not it is taken.
        if constexpr (num_traits<T>::has_transcendental) {
            const T n = num_traits<T>::from_int(static_cast<long>(d.trace.size()));
            T m1 = num_traits<T>::from_int(0), m2 = num_traits<T>::from_int(0),
              m3 = num_traits<T>::from_int(0);
            for (const T& x : d.trace) {
                m1 += x;
                m2 += T(x * x);
                m3 += T(x * x * x);
            }
            return mam::aph_fit(T(m1 / n), T(m2 / n), T(m3 / n)).aph;
        } else {
            throw UnsupportedError(
                "dist_to_map: fitting a Replayer trace to an acyclic phase-type needs "
                "transcendental arithmetic, which exact rational arithmetic does not provide; "
                "solve in double or Real<n>, or replace the trace by a fitted distribution");
        }
    }
    return mam::map_erlang(d.mean, convert_to_map_phases(num_traits<T>::to_double(d.scv)));
}

/** `sn.pie`: the phase distribution seen by an arriving job. */
template <class T>
std::vector<T> dist_pie(const Distrib<T>& d) {
    return mam::map_pie(dist_to_map(d));
}

/**
 * Fill in the first two moments of a distribution given by its matrices.
 *
 * `Distrib::map_dist` cannot compute them -- they need the stationary vector --
 * so a MAP built directly from (D0,D1) leaves mean and scv at their defaults
 * until this runs. Every builder call that installs such a distribution passes
 * through here.
 */
template <class T>
void dist_refresh_moments(Distrib<T>& d) {
    if (d.disabled || !d.has_map()) return;
    if (d.type == ProcessType::DMAP) {
        dmap_refresh_moments(d);
        return;
    }
    mam::Map<T> m;
    m.D0 = d.D0;
    m.D1 = d.D1;
    d.mean = mam::map_mean(m);
    d.scv = mam::map_scv(m);
}

/**
 * `sn.lst`: the Laplace-Stieltjes transform E[exp(-sX)].
 *
 * The phase-type families evaluate the closed form pie (sI - D0)^-1 (-D0) e;
 * for a MAP that is the transform of its stationary interarrival time, which is
 * the quantity the M/G/1 analyzers want.
 */
template <class T>
T dist_lst(const Distrib<T>& d, const T& s) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (d.is_prior()) reject_prior("dist_lst");
    if (d.disabled) throw InputError("dist_lst: the distribution is disabled");
    if (s == zero) return one;

    switch (d.type) {
        case ProcessType::IMMEDIATE:
            return one;
        case ProcessType::REPLAYER: {
            // The EMPIRICAL transform mean(exp(-s x)) over the trace, which is
            // what MATLAB `Replayer.evalLST` and the JAR/python twins return.
            // Falling through to the phase-type arm would transform the APH fit
            // of the first three moments instead, i.e. a different law: on
            // gallery_replayerm1 that moves the G/M/1 caudal root from 0.332835
            // to 0.332887 and QLen by 1e-4 relative.
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Replayer trace is a mean of exp, which exact "
                    "arithmetic has no representation for; use the double or real backend");
            } else {
                const double sv = num_traits<T>::to_double(s);
                double acc = 0.0;
                for (const T& x : d.trace) acc += std::exp(-sv * num_traits<T>::to_double(x));
                return num_traits<T>::from_double(acc / static_cast<double>(d.trace.size()));
            }
        }
        case ProcessType::DET: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Det is exp(-s m), which exact arithmetic has "
                    "no representation for; use the double or real backend");
            } else {
                return num_traits<T>::from_double(
                    std::exp(-num_traits<T>::to_double(s) * num_traits<T>::to_double(d.mean)));
            }
        }
        case ProcessType::UNIFORM: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Uniform evaluates exp, which exact arithmetic "
                    "has no representation for; use the double or real backend");
            } else {
                const double a = num_traits<T>::to_double(d.params[0]);
                const double b = num_traits<T>::to_double(d.params[1]);
                const double sv = num_traits<T>::to_double(s);
                return num_traits<T>::from_double((std::exp(-sv * a) - std::exp(-sv * b)) /
                                                  (sv * (b - a)));
            }
        }
        case ProcessType::NORMAL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Normal is a value of exp, which exact arithmetic "
                    "has no representation for; use the double or real backend");
            } else {
                // exp(-mu s + sigma^2 s^2 / 2), the Laplace-Stieltjes transform.
                //
                // ONE DELIBERATE DIVERGENCE, and it is a defect on the other
                // side: `Normal.m:96-106` returns exp(+mu s + sigma^2 s^2 / 2),
                // which is the MGF at +s and not E[exp(-sX)] at all -- its own
                // comment says "the moment-generating function evaluated at -s",
                // which that expression is also not. Nothing reads it: a Normal
                // is never a service process, so no analyzer reaches this arm,
                // and reproducing the sign would put a wrong transform in the
                // one place a future caller would trust. Recorded in BUGS.md.
                const double mu = num_traits<T>::to_double(d.params[0]);
                const double sg = num_traits<T>::to_double(d.params[1]);
                const double sv = num_traits<T>::to_double(s);
                return num_traits<T>::from_double(std::exp(-mu * sv + sg * sg * sv * sv / 2.0));
            }
        }
        case ProcessType::PARETO: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Pareto is an integral of exp, which exact "
                    "arithmetic has no representation for; use the double or real backend");
            } else {
                // alpha * int_0^1 u^(alpha-1) exp(-s k / u) du, the substitution
                // x = k/u the reference uses to stay accurate as s -> 0.
                const T alpha = d.params[0], k = d.params[1];
                const T tiny = num_traits<T>::from_double(1e-300);
                auto g = [&alpha, &k, &s, &tiny, &one](const T& u) -> T {
                    if (!(u > num_traits<T>::from_int(0))) return num_traits<T>::from_int(0);
                    const T uu = u > tiny ? u : tiny;
                    const double e = std::exp(-num_traits<T>::to_double(T(s * k / uu)));
                    const double p = std::pow(num_traits<T>::to_double(u),
                                              num_traits<T>::to_double(T(alpha - one)));
                    return num_traits<T>::from_double(p * e);
                };
                const T val = qsys::detail::num_integral<T>(g, zero, one,
                                                    num_traits<T>::from_double(1e-12),
                                                    num_traits<T>::from_double(1e-300), 50);
                return T(alpha * val);
            }
        }
        case ProcessType::WEIBULL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Weibull is a numerical integral, which exact "
                    "arithmetic has no representation for; use the double or real backend");
            } else {
                // The reference's 1000-panel right-endpoint sum, reproduced.
                const double a = num_traits<T>::to_double(d.params[0]);
                const double r = num_traits<T>::to_double(d.params[1]);
                const double sv = num_traits<T>::to_double(s);
                const double upper = a * std::pow(-std::log(1e-10), 1.0 / r);
                const int n = 1000;
                const double dx = upper / n;
                double acc = 0.0;
                for (int i = 1; i <= n; ++i) {
                    const double x = i * dx;
                    const double pdf =
                        (r / a) * std::pow(x / a, r - 1.0) * std::exp(-std::pow(x / a, r));
                    acc += std::exp(-sv * x) * pdf;
                }
                return num_traits<T>::from_double(acc * dx);
            }
        }
        case ProcessType::LOGNORMAL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Lognormal is a numerical integral, which exact "
                    "arithmetic has no representation for; use the double or real backend");
            } else {
                const double mu = num_traits<T>::to_double(d.params[0]);
                const double sg = num_traits<T>::to_double(d.params[1]);
                const double sv = num_traits<T>::to_double(s);
                const double upper = std::exp(mu + 5.0 * sg);
                const int n = 1000;
                const double dx = upper / n;
                double acc = 0.0;
                for (int i = 1; i <= n; ++i) {
                    const double x = i * dx;
                    const double lx = std::log(x);
                    const double pdf = std::exp(-(lx - mu) * (lx - mu) / (2.0 * sg * sg)) /
                                       (x * sg * std::sqrt(2.0 * M_PI));
                    acc += std::exp(-sv * x) * pdf;
                }
                return num_traits<T>::from_double(acc * dx);
            }
        }
        case ProcessType::GAMMA: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_lst: the transform of a Gamma is (1 + s theta)^-k, which exact "
                    "arithmetic has no representation for; use the double or real backend");
            } else {
                const double shape = num_traits<T>::to_double(d.params[0]);
                const double scale = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(
                    std::pow(1.0 + num_traits<T>::to_double(s) * scale, -shape));
            }
        }
        default:
            break;
    }

    // Phase-type / MAP families: pie (sI - D0)^-1 (-D0) e, a rational function
    // of s and therefore exact wherever the arithmetic is.
    const mam::Map<T> m = dist_to_map(d);
    const std::vector<T> pie = mam::map_pie(m);
    const std::size_t n = m.D0.rows();
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = T((i == j ? s : zero) - m.D0(i, j));
    // rhs = (-D0) e, the exit-rate vector
    std::vector<T> rhs(n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T acc = zero;
        for (std::size_t j = 0; j < n; ++j) acc += m.D0(i, j);
        rhs[i] = T(-acc);
    }
    // Solve A x = rhs by Gaussian elimination with partial pivoting.
    std::vector<T> x = rhs;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best = col;
        double bv = std::fabs(num_traits<T>::to_double(A(col, col)));
        for (std::size_t r = col + 1; r < n; ++r) {
            const double v = std::fabs(num_traits<T>::to_double(A(r, col)));
            if (v > bv) {
                bv = v;
                best = r;
            }
        }
        if (best != col) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(col, j), A(best, j));
            std::swap(x[col], x[best]);
        }
        if (A(col, col) == zero) throw NumericError("dist_lst: singular transform matrix");
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) continue;
            const T f = T(A(r, col) / A(col, col));
            if (f == zero) continue;
            for (std::size_t j = 0; j < n; ++j) A(r, j) = T(A(r, j) - f * A(col, j));
            x[r] = T(x[r] - f * x[col]);
        }
    }
    T out = zero;
    for (std::size_t i = 0; i < n; ++i) out += pie[i] * T(x[i] / A(i, i));
    return out;
}

// Forward declarations: the complex transform below reaches the CDF and the raw
// moments, both defined further down, and a dependent call would resolve only by
// ADL at instantiation.
template <class T>
T dist_cdf(const Distrib<T>& d, const T& x);
template <class T>
T dist_moment(const Distrib<T>& d, unsigned k);

/**
 * `sn.lst` at a COMPLEX argument, E[exp(-sX)] with s off the real axis.
 *
 * WHY A SECOND OVERLOAD. A transform is evaluated off the real axis by anything
 * that inverts it or locates its roots: the Abate-Whitt Euler sum walks the line
 * Re(s) = A/(2t), and a matrix transform int exp(Ut) dF(t) is read off the
 * spectrum of U, which is complex in general. `dist_lst(d, s)` above is
 * templated on the arithmetic type T and returns T, so it cannot answer either;
 * this twin fixes the argument and the result at std::complex<double>, since a
 * complex transform is meaningless without transcendental arithmetic anyway.
 * Parity note: the JAR carries the same capability by widening `sn.lst` to
 * SerializableFunction<Complex, Complex>, MATLAB by its own closed forms, and
 * python by Distribution.evalLST accepting a complex argument.
 *
 * THE TIERS mirror the real overload exactly: a closed form where the family has
 * one, the phase-type solve where the law is Markovian, and the CDF-increment
 * sum otherwise -- the last being a proper measure for ANY law, including one
 * with an atom and one with no density.
 */
template <class T>
std::complex<double> dist_lst(const Distrib<T>& d, const std::complex<double>& s) {
    static_assert(num_traits<T>::has_transcendental,
                  "dist_lst at a complex argument requires transcendental arithmetic");
    if (d.is_prior()) reject_prior("dist_lst");
    if (d.disabled) throw InputError("dist_lst: the distribution is disabled");
    if (std::abs(s) == 0.0) return std::complex<double>(1.0, 0.0);

    switch (d.type) {
        case ProcessType::IMMEDIATE:
            return std::complex<double>(1.0, 0.0);
        case ProcessType::REPLAYER: {
            std::complex<double> acc(0.0, 0.0);
            for (const T& x : d.trace) acc += std::exp(-s * num_traits<T>::to_double(x));
            return acc / static_cast<double>(d.trace.size());
        }
        case ProcessType::DET:
            return std::exp(-s * num_traits<T>::to_double(d.mean));
        case ProcessType::UNIFORM: {
            const double a = num_traits<T>::to_double(d.params[0]);
            const double b = num_traits<T>::to_double(d.params[1]);
            return (std::exp(-s * a) - std::exp(-s * b)) / (s * (b - a));
        }
        case ProcessType::NORMAL: {
            const double mu = num_traits<T>::to_double(d.params[0]);
            const double sg = num_traits<T>::to_double(d.params[1]);
            return std::exp(-mu * s + sg * sg * s * s / 2.0);
        }
        case ProcessType::GAMMA: {
            const double shape = num_traits<T>::to_double(d.params[0]);
            const double scale = num_traits<T>::to_double(d.params[1]);
            return std::pow(std::complex<double>(1.0, 0.0) + s * scale, -shape);
        }
        default:
            break;
    }

    if (process_is_markovian(d.type)) {
        // pie (sI - D0)^-1 (-D0) e, the same rational function as the real
        // overload, continued to the complex plane.
        const mam::Map<T> m = dist_to_map(d);
        const std::vector<T> pie = mam::map_pie(m);
        const std::size_t n = m.D0.rows();
        std::vector<std::vector<std::complex<double> > > A(
            n, std::vector<std::complex<double> >(n, std::complex<double>(0.0, 0.0)));
        std::vector<std::complex<double> > x(n, std::complex<double>(0.0, 0.0));
        for (std::size_t i = 0; i < n; ++i) {
            double exit = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                const double d0 = num_traits<T>::to_double(m.D0(i, j));
                A[i][j] = (i == j ? s : std::complex<double>(0.0, 0.0)) - d0;
                exit += d0;
            }
            x[i] = std::complex<double>(-exit, 0.0);
        }
        for (std::size_t col = 0; col < n; ++col) {
            std::size_t piv = col;
            double bv = std::abs(A[col][col]);
            for (std::size_t r = col + 1; r < n; ++r) {
                if (std::abs(A[r][col]) > bv) { bv = std::abs(A[r][col]); piv = r; }
            }
            if (piv != col) { std::swap(A[piv], A[col]); std::swap(x[piv], x[col]); }
            if (std::abs(A[col][col]) == 0.0)
                throw NumericError("dist_lst: singular transform matrix");
            for (std::size_t r = col + 1; r < n; ++r) {
                const std::complex<double> f = A[r][col] / A[col][col];
                for (std::size_t c = col; c < n; ++c) A[r][c] -= f * A[col][c];
                x[r] -= f * x[col];
            }
        }
        for (std::size_t row = n; row-- > 0;) {
            std::complex<double> acc = x[row];
            for (std::size_t c = row + 1; c < n; ++c) acc -= A[row][c] * x[c];
            x[row] = acc / A[row][row];
        }
        std::complex<double> out(0.0, 0.0);
        for (std::size_t i = 0; i < n; ++i) out += num_traits<T>::to_double(pie[i]) * x[i];
        return out;
    }

    // Riemann-Stieltjes sum with true CDF increments, renormalized for the cut
    // tail, so the result is still a transform.
    const std::size_t n_grid = 2400;
    const double mean = num_traits<T>::to_double(dist_moment(d, 1u));
    double hi = mean * 60.0;
    const double m2 = num_traits<T>::to_double(dist_moment(d, 2u));
    const double var = m2 - mean * mean;
    if (std::isfinite(var) && var > 0.0) hi = std::max(hi, mean + 12.0 * std::sqrt(var));
    if (!std::isfinite(hi) || hi <= 0.0) return std::complex<double>(1.0, 0.0);
    const double step = hi / static_cast<double>(n_grid);
    std::complex<double> acc(0.0, 0.0);
    double mass = 0.0;
    double prev = num_traits<T>::to_double(dist_cdf(d, num_traits<T>::from_double(0.0)));
    for (std::size_t i = 0; i < n_grid; ++i) {
        const double right = static_cast<double>(i + 1) * step;
        const double cur = num_traits<T>::to_double(dist_cdf(d, num_traits<T>::from_double(right)));
        const double w = cur - prev;
        prev = cur;
        if (w == 0.0) continue;
        mass += w;
        acc += w * std::exp(-s * ((static_cast<double>(i) + 0.5) * step));
    }
    return mass > 0.0 ? acc / mass : std::complex<double>(1.0, 0.0);
}

/**
 * The k-th raw moment.
 *
 * The closed-form families are evaluated from their parameters, as the MATLAB
 * classes do, rather than from the Erlang approximation of `sn.proc`: the
 * approximation matches only the mean once the SCV exceeds 1.
 */
template <class T>
T dist_moment(const Distrib<T>& d, unsigned k) {
    const T one = num_traits<T>::from_int(1);
    if (k == 0) return one;
    if (d.is_prior()) reject_prior("dist_moment");
    if (d.disabled) throw InputError("dist_moment: the distribution is disabled");
    switch (d.type) {
        case ProcessType::IMMEDIATE:
            return num_traits<T>::from_int(0);
        case ProcessType::DET: {
            T v = one;
            for (unsigned i = 0; i < k; ++i) v = T(v * d.mean);
            return v;
        }
        case ProcessType::UNIFORM: {
            // (b^(k+1) - a^(k+1)) / ((k+1)(b-a))
            const T a = d.params[0], b = d.params[1];
            T pa = one, pb = one;
            for (unsigned i = 0; i <= k; ++i) {
                pa = T(pa * a);
                pb = T(pb * b);
            }
            return T((pb - pa) / (num_traits<T>::from_int(static_cast<long>(k) + 1) * (b - a)));
        }
        case ProcessType::PARETO: {
            // alpha k^m / (alpha - m), finite only for m < alpha
            const T alpha = d.params[0], scale = d.params[1];
            const T m = num_traits<T>::from_int(static_cast<long>(k));
            if (!(alpha > m))
                throw NumericError("dist_moment: the Pareto moment of this order is infinite");
            T ps = one;
            for (unsigned i = 0; i < k; ++i) ps = T(ps * scale);
            return T(alpha * ps / (alpha - m));
        }
        case ProcessType::REPLAYER: {
            const T n = num_traits<T>::from_int(static_cast<long>(d.trace.size()));
            T acc = num_traits<T>::from_int(0);
            for (const T& x : d.trace) {
                T v = one;
                for (unsigned i = 0; i < k; ++i) v = T(v * x);
                acc += v;
            }
            return T(acc / n);
        }
        case ProcessType::GAMMA: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_moment: the moments of a Gamma are values of the gamma function, which "
                    "exact arithmetic has no representation for");
            } else {
                const double shape = num_traits<T>::to_double(d.params[0]);
                const double scale = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(std::tgamma(shape + k) / std::tgamma(shape) *
                                                  std::pow(scale, static_cast<double>(k)));
            }
        }
        case ProcessType::WEIBULL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_moment: the moments of a Weibull are values of the gamma function, "
                    "which exact arithmetic has no representation for");
            } else {
                const double a = num_traits<T>::to_double(d.params[0]);
                const double r = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(std::pow(a, static_cast<double>(k)) *
                                                 std::tgamma(1.0 + k / r));
            }
        }
        case ProcessType::LOGNORMAL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_moment: the moments of a Lognormal are values of exp, which exact "
                    "arithmetic has no representation for");
            } else {
                const double mu = num_traits<T>::to_double(d.params[0]);
                const double sg = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(std::exp(k * mu + k * k * sg * sg / 2.0));
            }
        }
        case ProcessType::NORMAL: {
            // The raw moments from the recurrence m_k = mu m_{k-1} + (k-1)
            // sigma^2 m_{k-2}, which is exact in any arithmetic and needs no
            // double factorial: m_1 = mu and m_2 = mu^2 + sigma^2 seed it.
            const T mu = d.params[0], sg = d.params[1];
            T prev2 = one, prev1 = mu;
            if (k == 1) return prev1;
            for (unsigned i = 2; i <= k; ++i) {
                const T next = T(mu * prev1 +
                                 num_traits<T>::from_int(static_cast<long>(i) - 1) * sg * sg *
                                     prev2);
                prev2 = prev1;
                prev1 = next;
            }
            return prev1;
        }
        default:
            break;
    }
    return mam::map_moment(dist_to_map(d), k);
}

/**
 * F(x) = P{X <= x}, MATLAB's `Distribution.evalCDF`.
 *
 * WHY IT EXISTS AT ALL in a port whose solvers read moments and transforms: it
 * is the only thing `Prior.discretize` needs. The quadrature design of SolverUQ
 * places its nodes at the conditional medians of equal-mass strata of the
 * parameter density, which is an inverse CDF and nothing else, so a parameter
 * law of any family can be discretized with no per-family quantile.
 *
 * PER FAMILY, from the closed form the reference's own class uses -- the Erlang
 * from its Poisson sum, the Gamma from the regularized incomplete gamma, the
 * Pareto from `gpcdf` reduced to 1 - (k/x)^alpha -- rather than from the Erlang
 * approximation of `dist_to_map`, for the same reason `dist_moment` does: the
 * approximation matches only the mean once the SCV exceeds one. The phase-type
 * and MAP families fall through to 1 - pie exp(D0 x) e, which is `map_cdf`.
 *
 * ONE DELIBERATE DIVERGENCE FROM THE REFERENCE, and it is a defect on the other
 * side: `Uniform.evalCDF` in MATLAB returns the constant DENSITY 1/(b-a) inside
 * the support and 0 above it, so it is neither a CDF nor monotone. Reproducing
 * that would make the bisection below fail to bracket rather than return a
 * matching wrong number, and no ported quantity reads it, so the correct
 * (x-a)/(b-a) is computed here. Recorded in BUGS.md.
 */
template <class T>
T dist_cdf(const Distrib<T>& d, const T& x) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (d.is_prior())
        throw UnsupportedError(
            "dist_cdf: a Prior's law is the mixture over its alternatives, which lives with the "
            "Prior; call prior_cdf (lang/prior.h)");
    if (d.disabled) throw InputError("dist_cdf: the distribution is disabled");
    switch (d.type) {
        case ProcessType::IMMEDIATE:
            // Immediate.evalCDF returns 1 everywhere, the point mass at zero.
            return one;
        case ProcessType::DET:
            return x < d.params[0] ? zero : one;
        case ProcessType::UNIFORM: {
            const T a = d.params[0], b = d.params[1];
            if (x <= a) return zero;
            if (x >= b) return one;
            return T((x - a) / (b - a));
        }
        case ProcessType::NORMAL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_cdf: the Gaussian CDF is a value of erf, which exact arithmetic has no "
                    "representation for");
            } else {
                // 0.5 (1 + erf((x - mu) / (sigma sqrt 2))), `Normal.m:87-95`.
                const double mu = num_traits<T>::to_double(d.params[0]);
                const double sg = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(
                    0.5 * (1.0 + std::erf((num_traits<T>::to_double(x) - mu) /
                                          (sg * std::sqrt(2.0)))));
            }
        }
        case ProcessType::REPLAYER: {
            // The empirical CDF of the trace: the fraction of samples <= x.
            if (d.trace.empty()) throw InputError("dist_cdf: the Replayer carries no samples");
            std::size_t below = 0;
            for (const T& v : d.trace)
                if (!(v > x)) ++below;
            return T(num_traits<T>::from_int(static_cast<long>(below)) /
                     num_traits<T>::from_int(static_cast<long>(d.trace.size())));
        }
        default:
            break;
    }
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "dist_cdf: the law of this family is an exponential, which exact arithmetic has no "
            "representation for; use the double or real backend");
    } else {
        if (!(x > zero)) return zero;
        const double xv = num_traits<T>::to_double(x);
        switch (d.type) {
            case ProcessType::EXP: {
                const double lam = num_traits<T>::to_double(d.params[0]);
                return num_traits<T>::from_double(1.0 - std::exp(-lam * xv));
            }
            case ProcessType::ERLANG: {
                // 1 - sum_{j<r} exp(-alpha x) (alpha x)^j / j!, the reference's form.
                const double alpha = num_traits<T>::to_double(d.params[0]);
                const long r = std::lround(num_traits<T>::to_double(d.params[1]));
                const double z = alpha * xv;
                double term = std::exp(-z), acc = term;
                for (long j = 1; j < r; ++j) {
                    term *= z / static_cast<double>(j);
                    acc += term;
                }
                return num_traits<T>::from_double(1.0 - acc);
            }
            case ProcessType::HYPEREXP: {
                const double p = num_traits<T>::to_double(d.params[0]);
                const double m1 = num_traits<T>::to_double(d.params[1]);
                const double m2 = num_traits<T>::to_double(d.params[2]);
                return num_traits<T>::from_double(p * (1.0 - std::exp(-m1 * xv)) +
                                                  (1.0 - p) * (1.0 - std::exp(-m2 * xv)));
            }
            case ProcessType::PARETO: {
                const double alpha = num_traits<T>::to_double(d.params[0]);
                const double scale = num_traits<T>::to_double(d.params[1]);
                if (xv <= scale) return zero;
                return num_traits<T>::from_double(1.0 - std::pow(scale / xv, alpha));
            }
            case ProcessType::GAMMA: {
                const double shape = num_traits<T>::to_double(d.params[0]);
                const double scale = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(mam::gammainc_lower(shape, xv / scale));
            }
            case ProcessType::WEIBULL: {
                const double a = num_traits<T>::to_double(d.params[0]);
                const double r = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(1.0 - std::exp(-std::pow(xv / a, r)));
            }
            case ProcessType::LOGNORMAL: {
                const double mu = num_traits<T>::to_double(d.params[0]);
                const double sg = num_traits<T>::to_double(d.params[1]);
                return num_traits<T>::from_double(
                    0.5 * std::erfc(-(std::log(xv) - mu) / (sg * std::sqrt(2.0))));
            }
            default:
                break;
        }
        // The phase-type and MAP families: map_cdf, 1 - pie exp(D0 x) e.
        const mam::Map<T> m = dist_to_map(d);
        const std::vector<T> pie = mam::map_pie(m);
        Matrix<T> A(m.D0.rows(), m.D0.cols(), zero);
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) A(i, j) = T(m.D0(i, j) * x);
        const Matrix<T> E = expm(A);
        T acc = zero;
        for (std::size_t i = 0; i < E.rows(); ++i)
            for (std::size_t j = 0; j < E.cols(); ++j) acc += T(pie[i] * E(i, j));
        return T(one - acc);
    }
}

/**
 * The p-quantile, by bisection on `dist_cdf`.
 *
 * Port of `Prior.quantile`: bracketing starts at the mean and doubles outward,
 * which terminates for any law with a finite mean, and the search then halves
 * 200 times or until the bracket is within FineTol of its own width. Using only
 * the CDF is what makes it applicable to every family at once, which is the
 * reason `Prior.discretize` is written in probability space rather than in
 * parameter space.
 */
template <class T>
T dist_quantile(const Distrib<T>& d, const T& p) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (!(p > zero) || !(p < one))
        throw InputError("dist_quantile: p must lie strictly between 0 and 1");
    T lo = zero;
    T hi = d.mean > num_traits<T>::from_double(GlobalConstants::FineTol)
               ? d.mean
               : num_traits<T>::from_double(GlobalConstants::FineTol);
    const int max_expand = 200;
    // A law with mass BELOW ZERO needs the lower end walked down as well. Every
    // ProcessType family but one is a duration and starts at 0, which is why the
    // bracket did; a Normal parameter density does not, and leaving lo at 0
    // would have returned a non-negative "quantile" for any p under F(0) -- a
    // wrong stratum median for `Prior.discretize`, with no diagnostic.
    if (d.type == ProcessType::NORMAL) {
        // Walk both ends out from the mean in doubling multiples of sigma.
        T step = d.params.size() > 1 ? d.params[1] : one;
        lo = T(d.mean - step);
        hi = T(d.mean + step);
        int j = 0;
        for (; j < max_expand; ++j) {
            const bool low_ok = dist_cdf(d, lo) <= p, high_ok = dist_cdf(d, hi) >= p;
            if (low_ok && high_ok) break;
            step = T(step * two);
            if (!low_ok) lo = T(d.mean - step);
            if (!high_ok) hi = T(d.mean + step);
        }
        if (j == max_expand) throw NumericError("dist_quantile: failed to bracket the quantile");
    } else {
        int i = 0;
        for (; i < max_expand; ++i) {
            if (dist_cdf(d, hi) >= p) break;
            hi = T(hi * two);
        }
        if (i == max_expand) throw NumericError("dist_quantile: failed to bracket the quantile");
    }
    const T tol = num_traits<T>::from_double(GlobalConstants::FineTol);
    for (int k = 0; k < 200; ++k) {
        const T mid = T((lo + hi) / two);
        if (dist_cdf(d, mid) < p)
            lo = mid;
        else
            hi = mid;
        const T scale = hi > one ? hi : one;
        if (T(hi - lo) <= T(tol * scale)) break;
    }
    return T((lo + hi) / two);
}

/**
 * Port of `Replayer.isNHPP`: test whether a trace is a sample path of a
 * NON-HOMOGENEOUS POISSON process, by the conditional-uniform KS test with the
 * Lewis refinement (`infer_nhpp_ks`).
 *
 * WHY THE QUESTION IS WORTH ASKING. A Replayer is used wherever a measured
 * stream is fed to a solver, and every analytical method that consumes it as an
 * arrival process assumes SOMETHING about its dependence structure. This test
 * says whether the Poisson assumption -- independent increments, whatever the
 * rate does with time -- survives contact with the data, which is the
 * assumption a time-varying analysis (`mtginf`, `mol`, `tvms`) rests on. A
 * small p-value says the stream is not Poisson at any rate function, so those
 * methods are answering a different process.
 *
 * The trace holds INTER-ARRIVAL times, so the arrival epochs are their
 * cumulative sum and the horizon is the last of them.
 */
template <class T>
infer::NhppKsResult<T> dist_is_nhpp(const Distrib<T>& d) {
    if (d.trace.size() < 2)
        throw InputError("dist_is_nhpp: the trace needs at least two inter-arrival times to test");
    std::vector<T> epochs;
    epochs.reserve(d.trace.size());
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < d.trace.size(); ++i) {
        acc = T(acc + d.trace[i]);
        epochs.push_back(acc);
    }
    return infer::infer_nhpp_ks<T>(epochs, epochs.back());
}

}  // namespace lang
}  // namespace line

#endif  // LINE_LANG_DISTRIBUTION_H
