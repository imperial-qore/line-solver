/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_RUNLENGTH_H
#define LINE_API_SIM_RUNLENGTH_H

/**
 * Run-length planning for steady-state simulation.
 *
 * Templated port of matlab/src/api/sim/sim_runlength.m, sim_asymvar_mm1.m and
 * sim_asymvar_ctmc.m, cross-checked against
 * jar/src/main/java/jline/api/sim/SimRunlength.java.
 *
 * THE QUANTITY THAT MATTERS is not the variance of the process but its
 * ASYMPTOTIC VARIANCE sigma^2 = lim t Var(time-average over [0,t]), twice the
 * integral of the autocovariance: a time average of a positively correlated
 * process converges at rate sigma^2/t, not Var(X)/t. Then
 *
 *   t* = (z/eps)^2 sigma^2 / mean^2
 *
 * is the run needed for relative precision eps at confidence 1-alpha.
 *
 * For M/M/1, sigma^2 = 2 rho(1+rho)/(mu (1-rho)^4) in closed form; the FOURTH
 * power is the whole story, and dividing by the squared mean leaves a run length
 * growing like (1-rho)^-2. Checked here against the general CTMC deviation-vector
 * computation, which agrees to 1e-6 at rho up to 0.9.
 *
 * ARITHMETIC. The normal quantile needs erfc, so the planner is transcendental;
 * the CTMC asymptotic variance is a linear solve and stays exact.
 *
 * Reference: W. Whitt (1989). Planning queueing simulations. Management Science
 * 35(11), 1341-1366.
 */

#include <cmath>
#include <limits>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace sim {

/** Second-order description of a steady-state estimator. */
template <class T>
struct AsymVarResult {
    T mean;                ///< the steady-state mean
    T variance;            ///< Var of the process itself
    T asymptoticVariance;  ///< sigma^2, what the run length depends on
    T relaxationTime;      ///< sigma^2/Var, the correlation time scale
    std::vector<T> deviation;  ///< the deviation vector, for the CTMC form
};

/** Outcome of the run-length plan. */
template <class T>
struct RunLengthResult {
    T requiredRunLength;      ///< t*
    T z;                      ///< the two-sided normal quantile used
    T halfWidth;              ///< the half-width a supplied run buys
    T achievedRelPrecision;   ///< that half-width over the mean
    bool hasRun = false;      ///< whether a run length was supplied
};

/**
 * Asymptotic variance of the M/M/1 number-in-system process.
 *
 * @param lambda arrival rate
 * @param mu     service rate
 */
template <class T>
AsymVarResult<T> sim_asymvar_mm1(const T& lambda, const T& mu) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (lambda <= zero || mu <= zero)
        throw InputError("sim_asymvar_mm1: the arrival and service rates must be positive");
    const T rho = lambda / mu;
    if (rho >= one) throw InputError("sim_asymvar_mm1: the queue must be stable, rho < 1");
    AsymVarResult<T> r;
    r.mean = rho / (one - rho);
    r.variance = rho / ((one - rho) * (one - rho));
    const T d = (one - rho) * (one - rho) * (one - rho) * (one - rho);
    r.asymptoticVariance = two * rho * (one + rho) / (mu * d);
    r.relaxationTime = r.variance > zero ? T(r.asymptoticVariance / r.variance) : zero;
    return r;
}

/**
 * Asymptotic variance of a reward on a CTMC: 2 sum_x pi(x)g(x)d(x) with
 * g = f - E_pi[f] and A d = -g, pi d = 0. The normalization is what pins d: A
 * alone is singular, since a constant may be added without changing sigma^2.
 *
 * @param A  the generator, rows summing to zero
 * @param f  the reward attached to each state
 * @param pi the stationary distribution; solved for when empty
 */
template <class T>
AsymVarResult<T> sim_asymvar_ctmc(const Matrix<T>& A, const std::vector<T>& f,
                                  const std::vector<T>& pi = std::vector<T>()) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("sim_asymvar_ctmc: the generator must be square");
    if (f.size() != n) throw InputError("sim_asymvar_ctmc: one reward per state is required");
    for (std::size_t i = 0; i < n; ++i) {
        T row = zero;
        for (std::size_t j = 0; j < n; ++j) row += A(i, j);
        if (num_abs(row) > num_traits<T>::from_double(1e-8))
            throw InputError("sim_asymvar_ctmc: the generator rows must sum to zero");
    }
    std::vector<T> p = pi;
    if (p.empty()) {
        // pi A = 0 with sum pi = 1, as a square system.
        Matrix<T> M(n, n);
        std::vector<T> b(n, zero);
        for (std::size_t j = 0; j + 1 < n; ++j)
            for (std::size_t i = 0; i < n; ++i) M(j, i) = A(i, j);
        for (std::size_t i = 0; i < n; ++i) M(n - 1, i) = one;
        b[n - 1] = one;
        p = solve(M, b);
    }
    T mean = zero;
    for (std::size_t i = 0; i < n; ++i) mean += p[i] * f[i];
    std::vector<T> g(n);
    for (std::size_t i = 0; i < n; ++i) g[i] = f[i] - mean;
    // A d = -g pins d only up to a constant, so one equation of A is redundant
    // and one normalization replaces it. WHICH equation is dropped matters: the
    // rows of A are related by pi A = 0, so a row whose pi is tiny is only
    // nominally redundant, and dropping it loses real information -- on a queue
    // truncated where pi has underflowed, that alone puts sigma^2 out by orders
    // of magnitude. Dropping the row with the LARGEST pi is the well-conditioned
    // choice.
    std::size_t drop = 0;
    for (std::size_t i = 1; i < n; ++i)
        if (p[i] > p[drop]) drop = i;
    Matrix<T> M2(n, n);
    std::vector<T> b2(n, zero);
    std::size_t r0 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (i == drop) continue;
        for (std::size_t j = 0; j < n; ++j) M2(r0, j) = A(i, j);
        b2[r0] = -g[i];
        ++r0;
    }
    for (std::size_t j = 0; j < n; ++j) M2(n - 1, j) = p[j];
    const std::vector<T> d = solve(M2, b2);
    AsymVarResult<T> r;
    r.mean = mean;
    r.variance = zero;
    r.asymptoticVariance = zero;
    for (std::size_t i = 0; i < n; ++i) {
        r.variance += p[i] * g[i] * g[i];
        r.asymptoticVariance += p[i] * g[i] * d[i];
    }
    r.asymptoticVariance *= num_traits<T>::from_int(2);
    r.relaxationTime = r.variance > zero ? T(r.asymptoticVariance / r.variance) : zero;
    r.deviation = d;
    return r;
}

/**
 * Run length for a steady-state estimate of a given relative precision.
 *
 * @param mean         the steady-state mean being estimated
 * @param asymVar      sigma^2 of that estimator
 * @param relPrecision the target half-width as a fraction of the mean
 * @param confidence   the confidence level of the interval
 * @param runLength    an actual run length, to report the precision it buys;
 *                     non-positive to skip
 */
template <class T>
RunLengthResult<T> sim_runlength(const T& mean, const T& asymVar,
                                 const T& relPrecision = num_traits<T>::from_rational(1, 20),
                                 const T& confidence = num_traits<T>::from_rational(19, 20),
                                 const T& runLength = num_traits<T>::from_int(0)) {
    static_assert(num_traits<T>::has_transcendental, "sim_runlength needs erfc for the quantile");
    using std::erfc;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (mean == zero)
        throw InputError("sim_runlength: a relative precision is meaningless for a zero mean");
    if (asymVar < zero) throw InputError("sim_runlength: the asymptotic variance cannot be negative");
    if (relPrecision <= zero) throw InputError("sim_runlength: the relative precision must be positive");
    if (confidence <= zero || confidence >= one)
        throw InputError("sim_runlength: the confidence must lie in (0,1)");
    // Two-sided normal quantile, by bisection on erfc.
    T lo = zero, hi = num_traits<T>::from_int(40);
    const T target = one - confidence;
    for (int i = 0; i < 200; ++i) {
        const T mid = (lo + hi) / two;
        if (erfc(mid / sqrt(two)) > target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    RunLengthResult<T> r;
    r.z = (lo + hi) / two;
    r.requiredRunLength = (r.z / relPrecision) * (r.z / relPrecision) * asymVar / (mean * mean);
    r.halfWidth = zero;
    r.achievedRelPrecision = zero;
    if (runLength > zero) {
        r.hasRun = true;
        r.halfWidth = r.z * sqrt(asymVar / runLength);
        r.achievedRelPrecision = r.halfWidth / num_abs(mean);
    }
    return r;
}

/** The plan of `sim_runlength_plan`: what the run should have been. */
template <class T>
struct RunLengthPlan {
    T relPrecision;                    ///< the precision planned for
    T confidence;                      ///< the level the half-widths were computed at
    T samplesUsed;                     ///< the run length they came from
    Matrix<T> asymptoticVariance;      ///< sigma^2 per (station, class), NaN where unplannable
    Matrix<T> requiredSamples;         ///< the run length that reaches relPrecision
};

/**
 * How long a simulation run should have been, from the one it already did.
 *
 * A batch-means half-width H at confidence 1-alpha over a run of N samples pins
 * the ASYMPTOTIC variance of the estimator,
 *
 *   sigma^2 = (H/z)^2 N,   z = Phi^-1((1+confidence)/2),
 *
 * and that is the quantity a run length is planned from -- NOT the stationary
 * variance, which on M/M/1 differs from it by a factor blowing up like
 * (1-rho)^-2. `sim_runlength` then turns it into the sample count that reaches
 * a requested RELATIVE precision.
 *
 * An entry with a non-positive mean or half-width is left NaN, since there is
 * nothing to plan from there.
 *
 * Reference: W. Whitt (1989). Planning queueing simulations. Management Science
 * 35(11), 1341-1366.
 */
template <class T>
RunLengthPlan<T> sim_runlength_plan(const Matrix<T>& means, const Matrix<T>& ciHalfWidth,
                                    const T& samplesUsed,
                                    const T& relPrecision = num_traits<T>::from_rational(1, 20),
                                    const T& confidence = num_traits<T>::from_rational(19, 20)) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_runlength_plan needs erfc for the quantile");
    const T zero = num_traits<T>::from_int(0);
    if (samplesUsed <= zero)
        throw InputError("sim_runlength_plan: the number of samples already used must be positive");
    const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
    RunLengthPlan<T> out;
    out.relPrecision = relPrecision;
    out.confidence = confidence;
    out.samplesUsed = samplesUsed;
    out.asymptoticVariance = Matrix<T>(means.rows(), means.cols(), nan);
    out.requiredSamples = Matrix<T>(means.rows(), means.cols(), nan);
    // The same z the interval itself was built with; `sim_runlength` computes
    // it inline from erfc, so it is read back from there rather than recomputed
    // by a second rule.
    const T z = sim_runlength<T>(num_traits<T>::from_int(1), zero, num_traits<T>::from_int(1),
                                 confidence)
                    .z;
    for (std::size_t i = 0; i < means.rows(); ++i)
        for (std::size_t r = 0; r < means.cols(); ++r) {
            if (i >= ciHalfWidth.rows() || r >= ciHalfWidth.cols()) continue;
            const double h = num_traits<T>::to_double(ciHalfWidth(i, r));
            const double m = num_traits<T>::to_double(means(i, r));
            if (!std::isfinite(h) || h <= 0.0 || !std::isfinite(m) || m <= 0.0) continue;
            const T av = T(ciHalfWidth(i, r) / z * (ciHalfWidth(i, r) / z) * samplesUsed);
            out.asymptoticVariance(i, r) = av;
            out.requiredSamples(i, r) =
                sim_runlength<T>(means(i, r), av, relPrecision, confidence).requiredRunLength;
        }
    return out;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_RUNLENGTH_H
