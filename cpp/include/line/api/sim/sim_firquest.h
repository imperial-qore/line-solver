/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_FIRQUEST_H
#define LINE_API_SIM_SIM_FIRQUEST_H

/**
 * Fixed-sample-size quantile interval from independent replications.
 *
 * Port of matlab/src/api/sim/sim_firquest.m. FIRQUEST is the replicated
 * counterpart of FQUEST and differs from sim_fquest in four places:
 *
 *   The warmup randomness test runs independently on each replicate path, and
 *   the batch size it settles on may differ between replications.
 *
 *   Truncation removes the LARGEST of those batch sizes from the front of every
 *   replication, not just from one path. This is more aggressive than FQUEST on
 *   purpose: an untruncated transient common to all replications biases every
 *   replicate estimate the same way, and averaging cannot remove it.
 *
 *   The four stage tests act on the R*b signed areas and R*b replicate batched
 *   quantile estimators POOLED IN REPLICATION-MAJOR ORDER, i.e. all b statistics
 *   of replication 1, then those of replication 2, and so on, which is MATLAB's
 *   column-major areas(:). Both stage tests are order-sensitive, so pooling in
 *   any other order silently changes which batch count is selected.
 *
 *   The delivered interval is
 *     ytilde_p(N*) +- t_{1-alpha/2, 2Rb-1} sqrt(Vtilde_p(w;R,b,m)/N*),
 *   N* = R*b*m, with the pooled combined variance-parameter estimator
 *     A_p(w;R,b,m)    = (Rb)^{-1} sum_j A_p(w;j,m)^2
 *     Ntilde_p(R,b,m) = m (Rb-1)^{-1} sum_j (yhat_p(j,m) - ytilde_p(N*))^2
 *     Vtilde_p        = [Rb A_p + (Rb-1) Ntilde_p] / (2Rb-1).
 *   The heuristic fallback drops FQUEST's residual-autocorrelation correction,
 *   since the pooled batch quantiles come from independent paths.
 *
 * DEFAULTS DIFFER FROM FQUEST and that is why the options argument is a pointer
 * here: b0 = 25 rather than 50, and s is chosen from R rather than fixed,
 * because the stage tests act on the R*b pooled statistics and so need fewer
 * batches per replication. Passing nullptr takes both; passing a default
 * constructed QuestOptions takes the FQUEST values instead, which is a
 * different procedure, not a formality. Start from sim_firquest_options(R) when
 * overriding one field.
 *
 * Independent replications shorten the correlation the estimator has to fight,
 * and they parallelize, but they reintroduce initialization bias in every path,
 * so a short run length per replication is worse here than in sim_fquest. The
 * reference reports slight undercoverage at p = 0.99 when the total sample is
 * under 500000, down to 90.8%.
 *
 * Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec, A. C. Mokashi,
 * J. R. Wilson, "A Fixed-Sample-Size Procedure for Estimating Steady-State
 * Quantiles Based on Independent Replications", Proc. Winter Simulation
 * Conference, 2025.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/sim/sim_dist.h"
#include "line/api/sim/sim_fquest.h"
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

/**
 * The article's batch counts as a function of the replication count, chosen so
 * that R*b pooled statistics remain enough to test while every replication still
 * contributes at least one batch.
 *
 * @param R number of replications
 * @return the descending batch-count ladder
 */
inline std::vector<long> sim_firquest_batchcounts(std::size_t R) {
    if (R == 2) return std::vector<long>{14, 11, 8, 5};
    if (R == 3) return std::vector<long>{10, 8, 6, 4};
    if (R == 4) return std::vector<long>{6, 5, 4, 3};
    if (R < 10) return std::vector<long>{5, 4, 3, 2};
    if (R < 17) return std::vector<long>{4, 3, 2, 1};
    if (R < 23) return std::vector<long>{3, 2, 1};
    if (R < 33) return std::vector<long>{2, 1};
    return std::vector<long>{1};
}

/**
 * The FIRQUEST defaults at R replications: b0 = 25 and the R-dependent ladder,
 * every other constant as in FQUEST.
 *
 * @param R number of replications
 * @return the option set sim_firquest uses when none is supplied
 */
inline QuestOptions sim_firquest_options(std::size_t R) {
    QuestOptions opt;
    opt.b0 = 25;
    opt.s = sim_firquest_batchcounts(R);
    return opt;
}

namespace detail {

/**
 * Pooled statistics over the last b*m observations of every replication. The
 * areas and batch quantiles come back in replication-major order, see the note
 * at the top of this header.
 */
template <class T>
StsQuantileStats<T> firquest_pool(const std::vector<std::vector<T>>& Yt, long b, long m, double p,
                                  double weight) {
    const std::size_t R = Yt.size();
    const std::size_t bm = static_cast<std::size_t>(b) * static_cast<std::size_t>(m);

    StsQuantileStats<T> pooled;
    pooled.b = static_cast<std::size_t>(b);
    pooled.m = static_cast<std::size_t>(m);
    pooled.n = R * bm;
    pooled.areas.reserve(R * static_cast<std::size_t>(b));
    pooled.bqe.reserve(R * static_cast<std::size_t>(b));

    std::vector<T> all;
    all.reserve(pooled.n);
    for (std::size_t r = 0; r < R; ++r) {
        const std::vector<T> kept(Yt[r].end() - static_cast<std::ptrdiff_t>(bm), Yt[r].end());
        const StsQuantileStats<T> st = sim_sts_quantile_areas<T>(
            kept, static_cast<std::size_t>(b), static_cast<std::size_t>(m), p, weight);
        pooled.areas.insert(pooled.areas.end(), st.areas.begin(), st.areas.end());
        pooled.bqe.insert(pooled.bqe.end(), st.bqe.begin(), st.bqe.end());
        all.insert(all.end(), kept.begin(), kept.end());
    }

    std::sort(all.begin(), all.end());
    pooled.quantile =
        all[static_cast<std::size_t>(std::ceil(static_cast<double>(pooled.n) * p)) - 1];

    const std::size_t K = R * static_cast<std::size_t>(b);
    T sumSq = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i) sumSq += T(pooled.areas[i] * pooled.areas[i]);
    pooled.Ap = T(sumSq / num_traits<T>::from_int(static_cast<long>(K)));

    T sd = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < K; ++i) {
        const T d = T(pooled.bqe[i] - pooled.quantile);
        sd += T(d * d);
    }
    pooled.Np = T(num_traits<T>::from_int(m) * sd /
                  num_traits<T>::from_int(static_cast<long>(K - 1)));
    pooled.Vp = T((num_traits<T>::from_int(static_cast<long>(K)) * pooled.Ap +
                   num_traits<T>::from_int(static_cast<long>(K - 1)) * pooled.Np) /
                  num_traits<T>::from_int(static_cast<long>(2 * K - 1)));
    return pooled;
}

}  // namespace detail

/**
 * @param Y       R replicate sample paths of equal length, finite
 * @param p       quantile order in (0,1)
 * @param alpha   nominal non-coverage, 0.05 by default
 * @param options procedure constants, nullptr for the FIRQUEST defaults at this R
 */
template <class T>
QuestResult<T> sim_firquest(const std::vector<std::vector<T>>& Y, double p, double alpha = 0.05,
                            const QuestOptions* options = nullptr) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_firquest: the interval is a t quantile times a square root, so exact "
                  "arithmetic is refused");
    const std::size_t R = Y.size();
    if (R < 2)
        throw InputError("sim_firquest: at least 2 replications are required, use sim_fquest for a "
                         "single path");
    const std::size_t nRep = Y[0].size();
    if (nRep == 0)
        throw InputError("sim_firquest: the replicate paths must be nonempty and of equal length");
    for (std::size_t r = 0; r < R; ++r) {
        if (Y[r].size() != nRep)
            throw InputError("sim_firquest: the replicate paths must be nonempty and of equal "
                             "length");
        for (std::size_t i = 0; i < nRep; ++i)
            if (!detail::num_isfinite(Y[r][i]))
                throw InputError("sim_firquest: the sample paths must be finite");
    }

    const QuestOptions opt =
        sim_quest_options(options == nullptr ? sim_firquest_options(R) : *options);

    if (!(p > 0.0) || !(p < 1.0))
        throw InputError("sim_firquest: p must be a real scalar in (0,1)");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("sim_firquest: alpha must be a real scalar in (0,1)");
    if (static_cast<long>(R) * opt.s.back() < 3)
        throw InputError("sim_firquest: R*min(s) pooled batches is below the 3 the stage tests "
                         "need");

    QuestResult<T> res;
    res.R = R;

    // ---- warmup: one randomness loop per replicate path
    const long n = static_cast<long>(nRep);
    const long b0 = opt.b0;
    long mStart = opt.m0;
    if (n < b0 * mStart) mStart = n / b0;
    if (mStart < 1)
        throw InputError("sim_firquest: each replication is too short for the initial batch count "
                         "b0");

    long mMax = 0;
    bool failed = false;
    for (std::size_t r = 0; r < R; ++r) {
        long m = mStart;
        long ell = 1;
        bool atMax = false, passed = false;
        while (true) {
            const std::vector<T> head(Y[r].begin(),
                                      Y[r].begin() + static_cast<std::ptrdiff_t>(b0 * m));
            const StsQuantileStats<T> st = sim_sts_quantile_areas<T>(
                head, static_cast<std::size_t>(b0), static_cast<std::size_t>(m), p, opt.weight);
            const double sig =
                opt.beta * std::exp(-opt.eta * std::pow(static_cast<double>(ell - 1), opt.theta));
            if (!sim_vonneumann<T>(st.areas, sig).reject) {
                passed = true;
                break;
            }
            if (atMax) break;
            ++ell;
            m = detail::warmup_next_m(n, b0, m, atMax);
            if (m < 1) break;
        }
        failed = failed || !passed;
        mMax = std::max(mMax, m);
    }
    if (failed)
        res.warnings.push_back("the warmup randomness test could not be passed in every "
                               "replication, the replicate paths are too short");

    // ---- truncation: delete the longest warmup batch from every replication
    const long truncated = mMax > 0 ? mMax : 0;
    if (truncated >= n)
        throw InputError("sim_firquest: the warmup batch size exhausts the replication length");
    std::vector<std::vector<T>> Yt(R);
    for (std::size_t r = 0; r < R; ++r)
        Yt[r].assign(Y[r].begin() + static_cast<std::ptrdiff_t>(truncated), Y[r].end());
    const long nstarRep = static_cast<long>(Yt[0].size());

    // ---- batch-count selection on the pooled statistics
    std::size_t v = 0;
    long b = opt.s[v];
    long m = nstarRep / b;
    bool ok = true, havePooled = false;
    StsQuantileStats<T> pooled;
    for (int stage = 1; stage <= 4; ++stage) {
        while (true) {
            if (m < 1) {
                ok = false;
                break;
            }
            pooled = detail::firquest_pool<T>(Yt, b, m, p, opt.weight);
            havePooled = true;
            const std::vector<T>& sample = stage <= 2 ? pooled.areas : pooled.bqe;
            const bool reject = (stage % 2 == 1) ? sim_vonneumann<T>(sample, opt.beta).reject
                                                 : sim_shapirowilk<T>(sample, opt.beta).reject;
            if (!reject) break;
            ++v;
            if (v >= opt.s.size()) {
                ok = false;
                break;
            }
            b = opt.s[v];
            m = nstarRep / b;
        }
        if (!ok) break;
    }

    if (!havePooled || m < 1)
        throw InputError("sim_firquest: the replicate paths are too short to form min(s) batches "
                         "each");

    res.b = static_cast<std::size_t>(b);
    res.m = static_cast<std::size_t>(m);
    res.n = pooled.n;
    res.truncated = static_cast<std::size_t>(truncated);
    res.estimate = pooled.quantile;
    res.Ap = pooled.Ap;
    res.Np = pooled.Np;
    res.Vp = pooled.Vp;

    const T nst = num_traits<T>::from_int(static_cast<long>(pooled.n));
    if (ok) {
        const T t = num_traits<T>::from_double(
            sim_tinv(1.0 - alpha / 2.0, static_cast<double>(2 * static_cast<long>(R) * b - 1)));
        const T half = T(t * detail::num_sqrt(T(pooled.Vp / nst)));
        res.lower = T(res.estimate - half);
        res.upper = T(res.estimate + half);
        res.halfwidth = half;
        res.heuristic = false;
    } else {
        res.warnings.push_back("a randomness or normality test failed at b = " +
                               std::to_string(opt.s.back()) +
                               " per replication, the delivered interval is heuristic");
        res.heuristic = true;
        if (opt.force) {
            const QuestInterval<T> ci = sim_quest_heuristic_ci<T>(
                pooled.bqe, res.estimate, pooled.Ap, pooled.Np, pooled.n, alpha, false);
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

#endif  // LINE_API_SIM_SIM_FIRQUEST_H
