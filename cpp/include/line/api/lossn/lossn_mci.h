/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LOSSN_LOSSN_MCI_H
#define LINE_API_LOSSN_LOSSN_MCI_H

/**
 * Monte Carlo importance-sampling summation for product-form loss networks.
 *
 * Templated port of matlab/src/api/lossn/lossn_mci.m, implementing Ross and
 * Wang, "Monte Carlo Summation Applied to Product-Form Loss Networks",
 * Probability in the Engineering and Informational Sciences 6 (1992), 323-348.
 *
 * Links j = 1..J carry capacity C(j), routes r = 1..R carry offered load nu(r)
 * and need A(j, r) circuits on link j. The state n is feasible iff A n <= C,
 * the equilibrium law is product form with normalizing constant
 * g(C) = sum_{n feasible} prod_r nu_r^{n_r} / n_r!, and the class-r acceptance
 * probability is g(C - A(:, r)) / g(C). States are drawn from the truncated
 * Poisson importance law (Eq. 6) over {0..N_1} x ... x {0..N_R} with
 * N_r = min_j floor(C_j / A_jr), and the ratio estimators (Eq. 8) give g and
 * the blocking probabilities with delta-method confidence intervals.
 *
 * RANDOMNESS. There is no global generator here and no hidden seeding. The
 * caller passes its own engine by reference, or a seed, and the uniforms are
 * drawn as (gen() >> 11) * 2^-53 rather than through a distribution object, so
 * a given engine and seed reproduce the same estimate on every standard
 * library. The draw order mirrors the reference (all S samples of route 1,
 * then of route 2, and so on) so that the two implementations traverse the
 * same importance law in the same order, though their engines differ and their
 * sample paths therefore cannot be compared point by point. What CAN be
 * compared, and is, is convergence to a closed form: a single link with unit
 * circuit requirements is an Erlang loss system and its blocking probability
 * is Erlang B.
 *
 * ARITHMETIC. The importance weights are formed in log space, exactly as in
 * the reference, so log, exp, lgamma and the inverse error function are all
 * unavoidable; the header is gated on num_traits<T>::has_transcendental and is
 * instantiated at double and Real only. An exact instantiation would be
 * meaningless in any case: the estimate is a random variable.
 *
 * REFERENCE NOTE. MATLAB's var normalizes by S - 1 and so does the port, and
 * the covariance is likewise the S - 1 form; with S below 2 the estimator is
 * undefined and the port rejects it rather than dividing by zero, which the
 * reference does silently.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lossn {

/** Options of lossn_mci, MATLAB's options struct. */
template <class T>
struct LossnMciOptions {
    std::size_t samples = 100000;                          ///< S
    double alpha = 0.05;                                   ///< 1 - confidence level
    std::vector<T> gamma;                                  ///< importance parameters, empty = heuristic
};

/** Result of lossn_mci. */
template <class T>
struct LossnMciResult {
    std::vector<T> QLen;      ///< mean carried load per route
    std::vector<T> Loss;      ///< blocking probability per route
    double lG = 0.0;          ///< log of the estimated normalizing constant
    Matrix<T> acceptCI;       ///< (R x 2) acceptance confidence interval
    Matrix<T> lossCI;         ///< (R x 2) blocking confidence interval
    std::vector<T> acceptPoint;
    std::vector<T> lossPoint;
    T level = num_traits<T>::from_int(0);  ///< 1 - alpha
    std::size_t nsamples = 0;
};

namespace detail {

/** log(sum exp(x)), shifted by the maximum. */
template <class T>
T logsumexp(const std::vector<T>& x) {
    using std::exp;
    using std::log;
    if (x.empty()) throw InputError("logsumexp: empty argument");
    T m = x[0];
    for (std::size_t i = 1; i < x.size(); ++i)
        if (x[i] > m) m = x[i];
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < x.size(); ++i) s += exp(T(x[i] - m));
    return T(m + log(s));
}

/**
 * log(l!) accumulated term by term, which is what MATLAB's gammaln(l + 1)
 * evaluates to at the integer arguments this algorithm uses. Summing the logs
 * rather than taking the log of a factorial keeps the working precision and
 * cannot overflow, and it removes any dependency on a special-function
 * library: Boost.Math's lgamma promotes through __float128 for its precision
 * policy, which drags libquadmath into every binary that links this header.
 */
template <class T>
std::vector<T> log_factorials(std::size_t upTo) {
    using std::log;
    std::vector<T> lf(upTo + 1, num_traits<T>::from_int(0));
    for (std::size_t l = 1; l <= upTo; ++l)
        lf[l] = T(lf[l - 1] + log(num_traits<T>::from_int(static_cast<long>(l))));
    return lf;
}

/**
 * The standard-normal 1 - alpha/2 quantile, MATLAB's sqrt(2) erfinv(1 - alpha).
 *
 * Computed in double by bisection on std::erf, and only then converted to T.
 * That is not a precision compromise in disguise: alpha reaches this function
 * as a double, so the quantile it determines carries double information and no
 * more, whatever the working type of the estimator is. It is also a random
 * estimator's confidence half width, whose own statistical error dwarfs any
 * arithmetic one.
 */
inline double normal_quantile(double alpha) {
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("lossn_mci: the confidence level alpha must lie in (0, 1)");
    const double target = 1.0 - alpha;
    double lo = 0.0, hi = 10.0;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (std::erf(mid) < target)
            lo = mid;
        else
            hi = mid;
    }
    return std::sqrt(2.0) * 0.5 * (lo + hi);
}

/** A uniform on [0, 1) built from the engine directly, for reproducibility. */
template <class Rng>
double uniform01(Rng& gen) {
    return static_cast<double>(gen() >> 11) * (1.0 / 9007199254740992.0);
}

}  // namespace detail

/**
 * Estimate the normalizing constant and the blocking probabilities.
 *
 * @param nu  (R) offered load per route
 * @param A   (J x R) circuit requirement of link j for route r
 * @param C   (J) link capacity
 * @param opt options
 * @param gen the caller's random engine, advanced in place
 */
template <class T, class Rng>
LossnMciResult<T> lossn_mci(const std::vector<T>& nu, const Matrix<T>& A, const std::vector<T>& C,
                            const LossnMciOptions<T>& opt, Rng& gen) {
    static_assert(num_traits<T>::has_transcendental,
                  "lossn_mci requires transcendental arithmetic (log-space importance weights)");
    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;

    const std::size_t R = nu.size(), J = C.size();
    if (R == 0 || J == 0) throw InputError("lossn_mci: empty loss network");
    if (A.rows() != J || A.cols() != R) throw InputError("lossn_mci: A must be J x R");
    const std::size_t S = opt.samples;
    if (S < 2) throw InputError("lossn_mci: at least two samples are required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // per-route maximum feasible occupancy
    std::vector<std::size_t> N(R, 0);
    for (std::size_t k = 0; k < R; ++k) {
        bool any = false;
        T best = zero;
        for (std::size_t j = 0; j < J; ++j) {
            if (!(A(j, k) > zero)) continue;
            const T ratio = C[j] / A(j, k);
            if (!any || ratio < best) best = ratio;
            any = true;
        }
        if (!any) continue;
        const double bd = num_traits<T>::to_double(best);
        N[k] = bd <= 0.0 ? 0 : static_cast<std::size_t>(std::floor(bd));
    }

    // importance parameters, Section 3.4 heuristic
    std::vector<T> gamma(opt.gamma);
    if (gamma.empty()) {
        T delta = zero;
        for (std::size_t j = 0; j < J; ++j) {
            if (C[j] == zero) throw InputError("lossn_mci: a link has zero capacity");
            T load = zero;
            for (std::size_t k = 0; k < R; ++k) load += A(j, k) * nu[k];
            load = T(load / C[j]);
            if (j == 0 || load > delta) delta = load;
        }
        T base = T(one - num_traits<T>::from_double(0.15) * T(one - delta));
        const T floorBase = num_traits<T>::from_double(1e-6);
        if (base < floorBase) base = floorBase;
        gamma.assign(R, zero);
        for (std::size_t k = 0; k < R; ++k) {
            T b = zero;
            for (std::size_t j = 0; j < J; ++j)
                if (A(j, k) > b) b = A(j, k);
            gamma[k] = nu[k] * pow(base, b);
        }
    }
    if (gamma.size() != R) throw InputError("lossn_mci: gamma must have one entry per route");
    const T gammaFloor = num_traits<T>::from_double(1e-300);
    for (std::size_t k = 0; k < R; ++k)
        if (gamma[k] < gammaFloor) gamma[k] = gammaFloor;

    // normalization constant of the importance law, and its per-route cdfs
    T log_c = zero;
    std::vector<std::vector<T>> cdf(R);
    for (std::size_t k = 0; k < R; ++k) {
        std::vector<T> logterms(N[k] + 1);
        const std::vector<T> lf = detail::log_factorials<T>(N[k]);
        for (std::size_t l = 0; l <= N[k]; ++l)
            logterms[l] =
                T(num_traits<T>::from_int(static_cast<long>(l)) * log(gamma[k]) - lf[l]);
        const T lse = detail::logsumexp(logterms);
        log_c += lse;
        cdf[k].assign(N[k] + 1, zero);
        T acc = zero;
        for (std::size_t l = 0; l <= N[k]; ++l) {
            acc += exp(T(logterms[l] - lse));
            cdf[k][l] = acc;
        }
        cdf[k][N[k]] = one;  // guard rounding, as in the reference
    }

    // S i.i.d. draws, column by column
    std::vector<std::size_t> V(S * R, 0);
    for (std::size_t k = 0; k < R; ++k)
        for (std::size_t s = 0; s < S; ++s) {
            const T u = num_traits<T>::from_double(detail::uniform01(gen));
            std::size_t v = 0;
            for (std::size_t l = 0; l <= N[k]; ++l)
                if (u > cdf[k][l]) ++v;
            V[s * R + k] = v;
        }

    // feasibility indicators and log likelihood ratios
    std::vector<T> logratio(R);
    for (std::size_t k = 0; k < R; ++k) {
        if (!(nu[k] > zero)) throw InputError("lossn_mci: a route has non-positive offered load");
        logratio[k] = T(log(nu[k]) - log(gamma[k]));
    }
    std::vector<T> log_alpha(S, zero);
    std::vector<char> inOmega(S, 0);
    std::vector<char> inOmegaK(S * R, 0);
    std::vector<T> AV(J, zero);
    for (std::size_t s = 0; s < S; ++s) {
        for (std::size_t j = 0; j < J; ++j) {
            T t = zero;
            for (std::size_t k = 0; k < R; ++k)
                t += A(j, k) * num_traits<T>::from_int(static_cast<long>(V[s * R + k]));
            AV[j] = t;
        }
        bool ok = true;
        for (std::size_t j = 0; j < J && ok; ++j)
            if (AV[j] > C[j]) ok = false;
        inOmega[s] = ok ? 1 : 0;
        for (std::size_t k = 0; k < R; ++k) {
            bool okk = true;
            for (std::size_t j = 0; j < J && okk; ++j)
                if (AV[j] > T(C[j] - A(j, k))) okk = false;
            inOmegaK[s * R + k] = okk ? 1 : 0;
        }
        T la = zero;
        for (std::size_t k = 0; k < R; ++k)
            la += num_traits<T>::from_int(static_cast<long>(V[s * R + k])) * logratio[k];
        log_alpha[s] = la;
    }

    LossnMciResult<T> res;
    res.nsamples = S;
    res.level = T(one - num_traits<T>::from_double(opt.alpha));

    // normalizing constant
    std::vector<T> laO;
    for (std::size_t s = 0; s < S; ++s)
        if (inOmega[s]) laO.push_back(log_alpha[s]);
    if (laO.empty()) {
        res.lG = -std::numeric_limits<double>::infinity();
    } else {
        const T lse = detail::logsumexp(laO);
        res.lG = num_traits<T>::to_double(T(log_c + lse)) -
                 std::log(static_cast<double>(S));
    }

    // ratio estimators with shifted weights
    T M = zero;
    bool anyOmega = false;
    for (std::size_t s = 0; s < S; ++s)
        if (inOmega[s]) {
            if (!anyOmega || log_alpha[s] > M) M = log_alpha[s];
            anyOmega = true;
        }
    std::vector<T> w(S), Z(S);
    T meanZ = zero;
    for (std::size_t s = 0; s < S; ++s) {
        w[s] = exp(T(log_alpha[s] - M));
        Z[s] = inOmega[s] ? w[s] : zero;
        meanZ += Z[s];
    }
    meanZ = T(meanZ / num_traits<T>::from_int(static_cast<long>(S)));

    const T crit = num_traits<T>::from_double(detail::normal_quantile(opt.alpha));
    res.acceptPoint.assign(R, zero);
    res.acceptCI = Matrix<T>(R, 2, zero);
    res.lossCI = Matrix<T>(R, 2, zero);
    res.Loss.assign(R, zero);
    res.QLen.assign(R, zero);
    if (!(meanZ > zero))
        throw NumericError("lossn_mci: no sampled state was feasible, the estimator is undefined");

    T varZ = zero;
    for (std::size_t s = 0; s < S; ++s) varZ += T(Z[s] - meanZ) * T(Z[s] - meanZ);
    varZ = T(varZ / num_traits<T>::from_int(static_cast<long>(S - 1)));

    for (std::size_t k = 0; k < R; ++k) {
        std::vector<T> Y(S);
        T meanY = zero;
        for (std::size_t s = 0; s < S; ++s) {
            Y[s] = inOmegaK[s * R + k] ? w[s] : zero;
            meanY += Y[s];
        }
        meanY = T(meanY / num_traits<T>::from_int(static_cast<long>(S)));
        const T phi = T(meanY / meanZ);
        T varY = zero, covYZ = zero;
        for (std::size_t s = 0; s < S; ++s) {
            varY += T(Y[s] - meanY) * T(Y[s] - meanY);
            covYZ += T(Y[s] - meanY) * T(Z[s] - meanZ);
        }
        const T denom = num_traits<T>::from_int(static_cast<long>(S - 1));
        varY = T(varY / denom);
        covYZ = T(covYZ / denom);
        T sig2 = T((varY - num_traits<T>::from_int(2) * phi * covYZ + phi * phi * varZ) /
                   (num_traits<T>::from_int(static_cast<long>(S)) * meanZ * meanZ));
        if (sig2 < zero) sig2 = zero;
        const T half = T(crit * sqrt(sig2));
        res.acceptPoint[k] = phi;
        res.acceptCI(k, 0) = T(phi - half);
        res.acceptCI(k, 1) = T(phi + half);
        res.Loss[k] = T(one - phi);
        res.QLen[k] = nu[k] * phi;
        res.lossCI(k, 0) = T(one - res.acceptCI(k, 1));
        res.lossCI(k, 1) = T(one - res.acceptCI(k, 0));
    }
    res.lossPoint = res.Loss;
    return res;
}

/**
 * Overload seeding a local engine. Deterministic in the seed and independent
 * of any global state; the engine is created and destroyed here.
 */
template <class T>
LossnMciResult<T> lossn_mci(const std::vector<T>& nu, const Matrix<T>& A, const std::vector<T>& C,
                            const LossnMciOptions<T>& opt, std::uint64_t seed) {
    std::mt19937_64 gen(seed);
    return lossn_mci(nu, A, C, opt, gen);
}

}  // namespace lossn
}  // namespace line

#endif  // LINE_API_LOSSN_LOSSN_MCI_H
