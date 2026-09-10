/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_FQUEST_H
#define LINE_API_SIM_SIM_FQUEST_H

/**
 * Fixed-sample-size confidence interval for a steady-state quantile.
 *
 * Port of matlab/src/api/sim/sim_fquest.m. The sample path Y has arbitrary fixed
 * length; no sequential control of the run length is needed. The procedure is
 * FQUEST and has four blocks:
 *
 *   Warmup. Starting from b = b0 and m = m0 it computes the b signed STS areas
 *   of the batched quantile process and tests them for randomness with von
 *   Neumann's ratio at the DECAYING significance beta*exp(-eta*(l-1)^theta) on
 *   iteration l, growing m by sqrt(2) whenever the test rejects. The decay is
 *   what terminates the loop once m can no longer grow: with a fixed
 *   significance the last iteration would repeat forever at m = floor(N/b).
 *
 *   Truncation. The first batch is deleted, which is the entire warmup
 *   treatment; there is no separate transient detector.
 *
 *   Batch-count selection. With b stepping down through s and m = floor(N* / b),
 *   four tests must pass in order: von Neumann and Shapiro-Wilk on the signed
 *   areas, then von Neumann and Shapiro-Wilk on the batched quantile
 *   estimators. b only ever decreases, and a failure at the last entry of s
 *   ends the stage.
 *
 *   Delivery. When all four tests pass the interval is
 *     ytilde_p(n*) +- t_{1-alpha/2, 2b-1} sqrt(V_p(w;b,m)/n*).
 *   Otherwise the sample was too small, `heuristic` is true, and under
 *   options.force the interval returned is the union of the wider of the two
 *   single-component intervals and Willink's skewness- and correlation-adjusted
 *   asymmetric interval.
 *
 * WHAT THIS MAY BE RUN ON. Applicability is a condition on the output process,
 * not on the model that produced it: geometric moment contraction (Wu 2005), a
 * density positive and differentiable at the quantile of interest, short-range
 * dependence and an FCLT for the indicator process. Two exclusions follow and
 * neither raises an error, so they have to be observed by the caller. Do NOT use
 * this on integer-valued output such as a queue length: the marginal has no
 * density and the batched quantile has no Bahadur representation. And
 * heavy-tailed service, which can break geometric moment contraction, is outside
 * the theory. Use it on continuous output, that is response, waiting and
 * sojourn times.
 *
 * A DELIVERED INTERVAL MAY WELL BE THE HEURISTIC ONE: 22% to 65% of runs took
 * that fallback in the reference's own test bed, rising with p. Read
 * `heuristic` before quoting the half-width as an asymptotically justified one.
 * Coverage measured there on M/M/1 waiting times at N = 200000 over 500
 * replications was 95.2% at p = 0.5, 96.2% at p = 0.9 and 95.6% at p = 0.99
 * against a nominal 95%, with the delivered half-width exceeding the empirically
 * needed one by 1.16, 1.24 and 1.99 respectively.
 *
 * Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec, A. C. Mokashi,
 * J. R. Wilson, "A Fixed-Sample-Size Method for Estimating Steady-State
 * Quantiles", Proc. Winter Simulation Conference, 2023.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/sim/sim_dist.h"
#include "line/api/sim/sim_quest_heuristic_ci.h"
#include "line/api/sim/sim_quest_options.h"
#include "line/api/sim/sim_shapirowilk.h"
#include "line/api/sim/sim_sts_quantile_areas.h"
#include "line/api/sim/sim_types.h"
#include "line/api/sim/sim_vonneumann.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

/** Point estimate and interval delivered by the QUEST procedures. */
template <class T>
struct QuestResult {
    T estimate;                         ///< Full-sample empirical p-quantile of the truncated path
    T lower;                            ///< Lower confidence limit, NaN if refused
    T upper;                            ///< Upper confidence limit, NaN if refused
    T halfwidth;                        ///< (upper-lower)/2, attained only on average when asymmetric
    std::size_t b = 0;                  ///< Final batch count, per replication in sim_firquest
    std::size_t m = 0;                  ///< Final batch size
    std::size_t n = 0;                  ///< Observations the interval rests on
    std::size_t R = 1;                  ///< Replications, 1 for sim_fquest
    std::size_t truncated = 0;          ///< Observations deleted from the front of each path
    T Ap;                               ///< STS area variance-parameter estimator
    T Np;                               ///< NBQ variance-parameter estimator
    T Vp;                               ///< Combined variance-parameter estimator
    bool heuristic = false;             ///< True when a stage test failed
    std::vector<std::string> warnings;  ///< Diagnostics, empty on a clean run
};

namespace detail {

/**
 * One warmup iteration's batch size update, shared by the two QUEST procedures.
 * Returns the next batch size and sets atMax when the sample path can no longer
 * support a larger one.
 */
inline long warmup_next_m(long N, long b, long m, bool& atMax) {
    const long full = N / b;
    const long next = static_cast<long>(std::llround(static_cast<double>(m) * std::sqrt(2.0)));
    if (N < b * next && next != full) return full;
    if (N < b * next) {
        atMax = true;
        return full;
    }
    return next;
}

}  // namespace detail

/**
 * @param Y       one simulation sample path, finite
 * @param p       quantile order in (0,1)
 * @param alpha   nominal non-coverage, 0.05 by default
 * @param options procedure constants, see sim_quest_options
 */
template <class T>
QuestResult<T> sim_fquest(const std::vector<T>& Y, double p, double alpha = 0.05,
                          const QuestOptions& options = QuestOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_fquest: the interval is a t quantile times a square root, so exact "
                  "arithmetic is refused");
    const QuestOptions opt = sim_quest_options(options);

    if (!(p > 0.0) || !(p < 1.0))
        throw InputError("sim_fquest: p must be a real scalar in (0,1)");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("sim_fquest: alpha must be a real scalar in (0,1)");
    const long sLast = opt.s.back();
    if (sLast < 3)
        throw InputError("sim_fquest: the stage tests need at least 3 batches, so min(s) >= 3");

    const long N = static_cast<long>(Y.size());
    for (std::size_t i = 0; i < Y.size(); ++i)
        if (!detail::num_isfinite(Y[i]))
            throw InputError("sim_fquest: the sample path must be finite");
    if (N < 2 * sLast)
        throw InputError("sim_fquest: the sample path is too short for the stage tests");

    QuestResult<T> res;

    // ---- warmup: grow the batch size until the signed areas look random
    long b = opt.b0;
    long m = opt.m0;
    if (N < b * m) m = N / b;
    if (m < 1)
        throw InputError("sim_fquest: the sample path is too short for the initial batch count b0");

    long ell = 1;
    bool atMax = false, passed = false;
    while (true) {
        const std::vector<T> head(Y.begin(), Y.begin() + static_cast<std::ptrdiff_t>(b * m));
        const StsQuantileStats<T> st = sim_sts_quantile_areas<T>(
            head, static_cast<std::size_t>(b), static_cast<std::size_t>(m), p, opt.weight);
        const double sig =
            opt.beta * std::exp(-opt.eta * std::pow(static_cast<double>(ell - 1), opt.theta));
        if (!sim_vonneumann<T>(st.areas, sig).reject) {
            passed = true;
            break;
        }
        if (atMax) break;
        ++ell;
        m = detail::warmup_next_m(N, b, m, atMax);
        if (m < 1) break;
    }
    if (!passed)
        res.warnings.push_back("the warmup randomness test could not be passed at the largest "
                               "admissible batch size, the sample path is too short");

    // ---- truncation: delete the first batch
    const long truncated = m > 0 ? m : 0;
    const std::vector<T> Yt(Y.begin() + static_cast<std::ptrdiff_t>(truncated), Y.end());
    const long Nstar = static_cast<long>(Yt.size());

    // ---- batch-count selection: four tests in order, b only decreases
    std::size_t v = 0;
    b = opt.s[v];
    m = Nstar / b;
    bool ok = true, haveStats = false;
    StsQuantileStats<T> stats;
    for (int stage = 1; stage <= 4; ++stage) {
        while (true) {
            if (m < 1) {
                ok = false;
                break;
            }
            const std::vector<T> kept(Yt.end() - static_cast<std::ptrdiff_t>(b * m), Yt.end());
            stats = sim_sts_quantile_areas<T>(kept, static_cast<std::size_t>(b),
                                              static_cast<std::size_t>(m), p, opt.weight);
            haveStats = true;
            const std::vector<T>& sample = stage <= 2 ? stats.areas : stats.bqe;
            const bool reject = (stage % 2 == 1) ? sim_vonneumann<T>(sample, opt.beta).reject
                                                 : sim_shapirowilk<T>(sample, opt.beta).reject;
            if (!reject) break;
            ++v;
            if (v >= opt.s.size()) {
                ok = false;
                break;
            }
            b = opt.s[v];
            m = Nstar / b;
        }
        if (!ok) break;
    }

    if (!haveStats || m < 1)
        throw InputError("sim_fquest: the sample path is too short to form min(s) batches");

    res.b = static_cast<std::size_t>(b);
    res.m = static_cast<std::size_t>(m);
    res.n = stats.n;
    res.truncated = static_cast<std::size_t>(truncated);
    res.estimate = stats.quantile;
    res.Ap = stats.Ap;
    res.Np = stats.Np;
    res.Vp = stats.Vp;

    const T nst = num_traits<T>::from_int(static_cast<long>(stats.n));
    if (ok) {
        const T t = num_traits<T>::from_double(
            sim_tinv(1.0 - alpha / 2.0, static_cast<double>(2 * b - 1)));
        const T half = T(t * detail::num_sqrt(T(stats.Vp / nst)));
        res.lower = T(res.estimate - half);
        res.upper = T(res.estimate + half);
        res.halfwidth = half;
        res.heuristic = false;
    } else {
        res.warnings.push_back("a randomness or normality test failed at b = " +
                               std::to_string(opt.s.back()) +
                               ", the delivered interval is heuristic");
        res.heuristic = true;
        if (opt.force) {
            const QuestInterval<T> ci = sim_quest_heuristic_ci<T>(
                stats.bqe, res.estimate, stats.Ap, stats.Np, stats.n, alpha, true);
            res.lower = ci.lower;
            res.upper = ci.upper;
            res.halfwidth = T(T(ci.upper - ci.lower) / num_traits<T>::from_int(2));
        } else {
            res.lower = detail::num_nan<T>();
            res.upper = detail::num_nan<T>();
            res.halfwidth = detail::num_nan<T>();
        }
    }
    return res;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_FQUEST_H
