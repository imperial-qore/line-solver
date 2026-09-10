/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MCI_H
#define LINE_API_PFQN_MCI_H

/**
 * Monte Carlo Integration estimate of the normalizing constant of a closed
 * product-form network (Ross, Wang and Yao; MonteQueue 2.0).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mci.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_mci.java.
 *
 * The McKenna-Mitra integral form of the constant,
 *
 *   G(N) = 1/prod_r N_r! int_{R_+^M} e^{-sum_i v_i} prod_r (sum_i v_i D(i,r) + Z_r)^{N_r} dv,
 *
 * is estimated by importance sampling with an independent exponential proposal
 * of rate gamma_i per station. Writing V ~ prod_i Exp(gamma_i), one draw
 * contributes
 *
 *   lZ = -sum_i (1 - gamma_i) V_i - sum_i log gamma_i - sum_r log N_r!
 *        + sum_r N_r log(sum_i V_i D(i,r) + Z_r),
 *
 * and lG = log mean exp(lZ). The proposal rates come from a Bard-Schweitzer
 * pre-solve: gamma_i = max(0.01, 1 - U_i) in the IMCI variant, and the
 * saturation-aware gamma_i = 1/sqrt(max N) when U_i > 0.9 in the plain MCI
 * variant. Overshooting rates make the estimator heavy-tailed, which is why
 * the reference clamps them.
 *
 * Deviations from the reference, both deliberate:
 *
 *  - MATLAB caches the uniform matrix in a PERSISTENT variable VL and reuses
 *    it across calls (slicing VL(1:I,1:M)), so two successive calls with the
 *    same shape return the SAME estimate and a call with a larger I silently
 *    reuses the old columns. That is global hidden state; here the deviates
 *    come from the caller's generator, one draw per station per sample.
 *
 *  - The 'rm' variant is REFUSED for M > 1 rather than reproduced. Its rate
 *    line, tput = N./(sum(D,1)+Z+max(D,1)*(sum(N)-1)), broadcasts a (1 x R)
 *    row against the (M x R) matrix max(D,1) and yields an (M x R) "throughput"
 *    whose product D*tput' is (M x M); the subsequent loop then reads gamma
 *    off the wrong axis. The expression is only dimensionally meaningful for a
 *    single station, which is the repairman model the variant is named for, so
 *    that is the only shape accepted here. See the report accompanying this
 *    port; nothing has been substituted for the M > 1 case.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The value is a random variable, the
 * integrand is evaluated through logs, and the proposal rates come from a
 * tolerance-stopped Bard-Schweitzer solve; all three require transcendental
 * arithmetic.
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** The proposal-rate rules the reference selects between. */
enum class MciVariant {
    Imci,  ///< gamma = max(0.01, 1 - U), the MonteQueue 2.0 recommendation
    Mci,   ///< gamma = 1/sqrt(max N) where U > 0.9, else 1 - U
    Rm     ///< repairman: a single station, rates from the balanced bound
};

/**
 * @param D       (M x R) service demands
 * @param N       (R) population per class
 * @param Z       (R) think times; empty for none
 * @param samples number of Monte Carlo samples
 * @param variant proposal-rate rule
 * @param rng     explicit generator, advanced by the call
 */
template <class T>
NcResult<T> pfqn_mci(const Matrix<T>& D, const std::vector<int>& N, const std::vector<T>& Z,
                     std::size_t samples, MciVariant variant, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_mci requires transcendental arithmetic: it is a Monte Carlo estimator over "
                  "an exponential proposal, its integrand is evaluated in the log domain, and its "
                  "proposal rates come from a tolerance-stopped Bard-Schweitzer solve");

    const std::size_t M = D.empty() ? 0 : D.rows();
    const std::size_t R = N.size();
    if (!D.empty() && D.cols() != R)
        throw InputError("pfqn_mci: D and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_mci: Z has the wrong length");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_mci: negative population");
    if (samples == 0) throw InputError("pfqn_mci: at least one sample is required");

    const T zero = num_traits<T>::from_int(0);
    std::vector<T> Zv(R, zero);
    for (std::size_t r = 0; r < Z.size(); ++r) Zv[r] = Z[r];

    // ---- degenerate model: pure delay ------------------------------------
    T dsum = zero;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) dsum += D(i, r);
    if (M == 0 || dsum < num_traits<T>::from_double(1e-4)) {
        double lG = 0.0;
        for (std::size_t r = 0; r < R; ++r) {
            lG -= mc_log_factorial<T>(N[r]);
            if (N[r] > 0) lG += N[r] * num_traits<T>::log_as_double(Zv[r]);
        }
        return {mc_exp<T>(lG), lG};
    }

    // ---- proposal rates ---------------------------------------------------
    std::vector<double> gamma(M);
    if (variant == MciVariant::Rm) {
        if (M != 1)
            throw UnsupportedError(
                "pfqn_mci: the 'rm' variant is defined only for a single station; for M > 1 the "
                "reference's rate expression broadcasts a (1 x R) row against an (M x R) matrix "
                "and reads gamma off the wrong axis, so nothing faithful can be computed");
        long Ntot = 0;
        int Nmax = 0;
        for (int n : N) {
            Ntot += n;
            if (n > Nmax) Nmax = n;
        }
        T util = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T dmax = D(0, r) > num_traits<T>::from_int(1) ? D(0, r) : num_traits<T>::from_int(1);
            const T den = D(0, r) + Zv[r] + dmax * num_traits<T>::from_int(Ntot - 1);
            if (den == zero) throw NumericError("pfqn_mci: zero denominator in the 'rm' bound");
            util += D(0, r) * (num_traits<T>::from_int(N[r]) / den);
        }
        const double u = num_traits<T>::to_double(util);
        gamma[0] = u > 0.9 ? 1.0 / std::sqrt(static_cast<double>(Nmax)) : 1.0 - u;
    } else {
        std::vector<T> Nt(R), Zt(R);
        for (std::size_t r = 0; r < R; ++r) {
            Nt[r] = num_traits<T>::from_int(N[r]);
            Zt[r] = Zv[r];
        }
        const AmvaResult<T> bs = pfqn_bs(D, Nt, Zt);
        int Nmax = 0;
        for (int n : N)
            if (n > Nmax) Nmax = n;
        for (std::size_t i = 0; i < M; ++i) {
            T util = zero;
            for (std::size_t r = 0; r < R; ++r) util += D(i, r) * bs.XN[r];
            const double u = num_traits<T>::to_double(util);
            if (variant == MciVariant::Imci) {
                gamma[i] = 1.0 - u > 0.01 ? 1.0 - u : 0.01;
            } else {
                gamma[i] = u > 0.9 ? 1.0 / std::sqrt(static_cast<double>(Nmax)) : 1.0 - u;
            }
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        if (!(gamma[i] > 0.0))
            throw NumericError(
                "pfqn_mci: a proposal rate is non-positive, the station is saturated beyond what "
                "the variant's clamp covers");

    // ---- constant part of every log-weight ---------------------------------
    double lconst = 0.0;
    for (std::size_t i = 0; i < M; ++i) lconst -= std::log(gamma[i]);
    for (std::size_t r = 0; r < R; ++r) lconst -= mc_log_factorial<T>(N[r]);

    std::vector<double> lZ(samples);
    std::vector<T> V(M);
    for (std::size_t s = 0; s < samples; ++s) {
        double lz = lconst;
        for (std::size_t i = 0; i < M; ++i) {
            double u = mc_uniform01(rng);
            if (u <= 0.0) u = 1.0 / 9007199254740992.0;
            const double v = -std::log(u) / gamma[i];
            V[i] = num_traits<T>::from_double(v);
            lz -= (1.0 - gamma[i]) * v;
        }
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == 0) continue;
            T inner = Zv[r];
            for (std::size_t i = 0; i < M; ++i) inner += V[i] * D(i, r);
            lz += N[r] * num_traits<T>::log_as_double(inner);
        }
        lZ[s] = lz;
    }

    double lG = mc_logmeanexp(lZ);
    if (!std::isfinite(lG)) {
        // Floating-point range exception: the reference falls back to the
        // largest single log-weight, which is the Laplace-style lower bound.
        double m = -std::numeric_limits<double>::infinity();
        for (double v : lZ)
            if (v > m) m = v;
        lG = m;
    }
    return {mc_exp<T>(lG), lG};
}

/** Reference defaults: 1e5 samples, the IMCI proposal. */
template <class T>
NcResult<T> pfqn_mci(const Matrix<T>& D, const std::vector<int>& N, const std::vector<T>& Z,
                     McRng& rng) {
    return pfqn_mci(D, N, Z, static_cast<std::size_t>(100000), MciVariant::Imci, rng);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MCI_H
