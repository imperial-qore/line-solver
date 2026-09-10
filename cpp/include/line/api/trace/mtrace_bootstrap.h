/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_BOOTSTRAP_H
#define LINE_API_TRACE_MTRACE_BOOTSTRAP_H

/**
 * Block-bootstrap confidence intervals for the descriptors of a marked trace.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_bootstrap.m.
 *
 * The statistic is the whole descriptor vector the m3a fitters consume --
 * [pc, backward moment, forward moment, sigma] -- and the resampling is by
 * BLOCK, not by observation: the trace is cut into BN = floor(N / 50) contiguous
 * blocks and the blocks are resampled with replacement. That is the point of the
 * routine. An inter-arrival trace is autocorrelated, and resampling individual
 * observations would destroy exactly the dependence the descriptors measure,
 * giving intervals far too narrow for sigma and the forward moment.
 *
 * The intervals are BCa (bias-corrected and accelerated), which is what MATLAB's
 * `bootci` computes by default:
 *   z0 = Phi^-1( #{theta* < theta_hat} / R ),
 *   a  = sum(d^3) / (6 (sum d^2)^{3/2}) over the jackknife deviations d,
 *   alpha1,2 = Phi( z0 + (z0 -+ z_alpha) / (1 - a (z0 -+ z_alpha)) ),
 * and the endpoints are the alpha1 and alpha2 quantiles of the replicates. The
 * bias correction z0 and the acceleration a are what make BCa transformation
 * respecting, which matters here because several descriptors are probabilities.
 *
 * THE THREE CODEBASES DISAGREE ON WHAT THIS FUNCTION IS. MATLAB returns
 * confidence intervals as above. The JAR returns a statistics object with
 * configurable block size and seed. Native Python returns a resampled (T, A)
 * pair -- a sampler, not an estimator, and no interval at all. MATLAB is the
 * reference and is what is ported; the divergence is recorded rather than
 * papered over, because a caller reading the Python name expects a trace back.
 *
 * RANDOMNESS. `line::pfqn::McRng` by reference, the tree's convention. The
 * stream is not comparable with MATLAB's, so the oracle is distributional:
 * the interval must cover the point estimate and shrink as R grows.
 *
 * ARITHMETIC: transcendental, for the normal quantiles.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/api/sim/sim_dist.h"
#include "line/api/trace/mtrace_moment.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace trace {

/** Lower and upper BCa endpoints of each descriptor, and the point estimate. */
template <class T>
struct MtraceBootstrapResult {
    std::vector<T> estimate;  ///< the statistic on the whole trace
    std::vector<T> lower;     ///< lower confidence limit, same layout
    std::vector<T> upper;     ///< upper confidence limit, same layout
    std::size_t blocks = 0;   ///< BN, the number of blocks resampled
};

namespace bootdetail {

/**
 * The descriptor vector the reference bootstraps: the class probabilities, the
 * first backward moment, the first forward moment, then sigma flattened.
 */
template <class T>
std::vector<T> stat_vector(const std::vector<T>& Tv, const std::vector<int>& A) {
    const std::vector<unsigned> ord(1, 1u);
    const std::vector<T> pc = mtrace_pc<T>(A);
    const Matrix<T> B = mtrace_moment(Tv, A, ord, false, true);
    const Matrix<T> F = mtrace_moment(Tv, A, ord, true, true);
    const Matrix<T> S = mtrace_sigma<T>(A);
    std::vector<T> out;
    out.reserve(pc.size() + B.rows() + F.rows() + S.rows() * S.cols());
    for (std::size_t i = 0; i < pc.size(); ++i) out.push_back(pc[i]);
    for (std::size_t i = 0; i < B.rows(); ++i) out.push_back(B(i, 0));
    for (std::size_t i = 0; i < F.rows(); ++i) out.push_back(F(i, 0));
    // MATLAB's S(:) is column-major.
    for (std::size_t j = 0; j < S.cols(); ++j)
        for (std::size_t i = 0; i < S.rows(); ++i) out.push_back(S(i, j));
    return out;
}

/** Concatenate the listed blocks, in the order given. */
template <class T>
void assemble(const std::vector<T>& Tv, const std::vector<int>& A,
              const std::vector<std::size_t>& first, const std::vector<std::size_t>& len,
              const std::vector<std::size_t>& pick, std::vector<T>* t, std::vector<int>* a) {
    t->clear();
    a->clear();
    for (std::size_t i = 0; i < pick.size(); ++i) {
        const std::size_t b = pick[i];
        for (std::size_t k = 0; k < len[b]; ++k) {
            t->push_back(Tv[first[b] + k]);
            a->push_back(A[first[b] + k]);
        }
    }
}

/** The p-quantile of a sorted sample, by linear interpolation. */
inline double quantile(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) return 0.0;
    if (p <= 0.0) return sorted.front();
    if (p >= 1.0) return sorted.back();
    const double h = p * static_cast<double>(sorted.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(h));
    const std::size_t hi = std::min(lo + 1, sorted.size() - 1);
    return sorted[lo] + (h - static_cast<double>(lo)) * (sorted[hi] - sorted[lo]);
}

}  // namespace bootdetail

/**
 * @param Tv        inter-arrival times
 * @param A         class labels, one per arrival
 * @param rng       generator, advanced by the call
 * @param resamples number of bootstrap replicates; the reference default is 1000
 * @param alpha     two-sided level; 0.05 gives a 95% interval
 * @param blockLen  target block length; the reference uses 50
 */
template <class T>
MtraceBootstrapResult<T> mtrace_bootstrap(const std::vector<T>& Tv, const std::vector<int>& A,
                                          pfqn::McRng& rng, std::size_t resamples = 1000,
                                          double alpha = 0.05, std::size_t blockLen = 50) {
    static_assert(num_traits<T>::has_transcendental,
                  "mtrace_bootstrap inverts the normal distribution");
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError("mtrace_bootstrap: the trace and its labels must agree in length");
    if (resamples < 2) throw InputError("mtrace_bootstrap: at least two replicates are required");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("mtrace_bootstrap: the level must lie strictly inside (0,1)");
    if (blockLen == 0) throw InputError("mtrace_bootstrap: the block length must be positive");

    const std::size_t N = Tv.size();
    const std::size_t BN = N / blockLen;
    if (BN < 2)
        throw InputError(
            "mtrace_bootstrap: the trace is too short to cut into at least two blocks at this "
            "block length; the block bootstrap has nothing to resample");

    // The reference's block layout: BN blocks of floor(N/BN), the first
    // mod(N, BN) of them one longer, so the blocks tile the trace exactly.
    const std::size_t base = N / BN, extra = N % BN;
    std::vector<std::size_t> len(BN, base), first(BN, 0);
    for (std::size_t b = 0; b < extra; ++b) len[b] += 1;
    for (std::size_t b = 1; b < BN; ++b) first[b] = first[b - 1] + len[b - 1];

    MtraceBootstrapResult<T> out;
    out.blocks = BN;
    out.estimate = bootdetail::stat_vector(Tv, A);
    const std::size_t P = out.estimate.size();

    // ---- the replicates -------------------------------------------------
    std::vector<std::vector<double>> rep(P);
    for (std::size_t i = 0; i < P; ++i) rep[i].reserve(resamples);
    std::vector<std::size_t> pick(BN, 0);
    std::vector<T> bt;
    std::vector<int> ba;
    for (std::size_t r = 0; r < resamples; ++r) {
        for (std::size_t b = 0; b < BN; ++b)
            pick[b] = static_cast<std::size_t>(pfqn::mc_uniform01(rng) * static_cast<double>(BN));
        for (std::size_t b = 0; b < BN; ++b)
            if (pick[b] >= BN) pick[b] = BN - 1;
        bootdetail::assemble(Tv, A, first, len, pick, &bt, &ba);
        std::vector<T> s;
        try {
            s = bootdetail::stat_vector(bt, ba);
        } catch (const Error&) {
            continue;  // a replicate that lost a class has no descriptor vector
        }
        if (s.size() != P) continue;  // ditto: the layout changed, so it is not comparable
        for (std::size_t i = 0; i < P; ++i) rep[i].push_back(num_traits<T>::to_double(s[i]));
    }

    // ---- the jackknife, for the acceleration -----------------------------
    std::vector<std::vector<double>> jack(P);
    std::vector<std::size_t> all;
    for (std::size_t b = 0; b < BN; ++b) all.push_back(b);
    for (std::size_t drop = 0; drop < BN; ++drop) {
        std::vector<std::size_t> keep;
        for (std::size_t b = 0; b < BN; ++b)
            if (b != drop) keep.push_back(b);
        bootdetail::assemble(Tv, A, first, len, keep, &bt, &ba);
        std::vector<T> s;
        try {
            s = bootdetail::stat_vector(bt, ba);
        } catch (const Error&) {
            continue;
        }
        if (s.size() != P) continue;
        for (std::size_t i = 0; i < P; ++i) jack[i].push_back(num_traits<T>::to_double(s[i]));
    }

    out.lower.assign(P, num_traits<T>::from_int(0));
    out.upper.assign(P, num_traits<T>::from_int(0));
    const double za = sim::sim_norminv(alpha / 2.0);
    for (std::size_t i = 0; i < P; ++i) {
        std::vector<double> v = rep[i];
        if (v.size() < 2) {  // nothing usable: report the point estimate twice
            out.lower[i] = out.estimate[i];
            out.upper[i] = out.estimate[i];
            continue;
        }
        std::sort(v.begin(), v.end());
        const double theta = num_traits<T>::to_double(out.estimate[i]);

        // Bias correction: the share of replicates below the point estimate.
        std::size_t below = 0;
        for (std::size_t k = 0; k < v.size(); ++k)
            if (v[k] < theta) ++below;
        double frac = static_cast<double>(below) / static_cast<double>(v.size());
        // Guard the endpoints, where the normal quantile is infinite.
        const double eps = 0.5 / static_cast<double>(v.size());
        if (frac < eps) frac = eps;
        if (frac > 1.0 - eps) frac = 1.0 - eps;
        const double z0 = sim::sim_norminv(frac);

        // Acceleration from the jackknife deviations.
        double acc = 0.0;
        if (jack[i].size() >= 2) {
            double mean = 0.0;
            for (std::size_t k = 0; k < jack[i].size(); ++k) mean += jack[i][k];
            mean /= static_cast<double>(jack[i].size());
            double s2 = 0.0, s3 = 0.0;
            for (std::size_t k = 0; k < jack[i].size(); ++k) {
                const double d = mean - jack[i][k];
                s2 += d * d;
                s3 += d * d * d;
            }
            if (s2 > 0.0) acc = s3 / (6.0 * std::pow(s2, 1.5));
        }

        auto endpoint = [&](double z) {
            const double num = z0 + z;
            const double den = 1.0 - acc * num;
            if (!(std::fabs(den) > 0.0)) return 0.5;
            return sim::sim_normcdf(z0 + num / den);
        };
        double a1 = endpoint(za), a2 = endpoint(-za);
        if (a1 > a2) std::swap(a1, a2);
        out.lower[i] = num_traits<T>::from_double(bootdetail::quantile(v, a1));
        out.upper[i] = num_traits<T>::from_double(bootdetail::quantile(v, a2));
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_BOOTSTRAP_H
