/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_STS_QUANTILE_AREAS_H
#define LINE_API_SIM_SIM_STS_QUANTILE_AREAS_H

/**
 * Standardized time series areas of the batched quantile process.
 *
 * Port of matlab/src/api/sim/sim_sts_quantile_areas.m. The b*m observations are
 * split into b nonoverlapping batches of size m; the function returns the
 * signed standardized time series (STS) areas of the quantile-estimation
 * process, the batched quantile estimators, and the three variance-parameter
 * estimators the QUEST procedures build confidence intervals from.
 *
 * With yhat_p(j,m) the empirical p-quantile of batch j and yhat_p(j,k) that of
 * its first k observations, the STS process of batch j is
 *   T_{j,m}(k/m) = (k/sqrt(m)) (yhat_p(j,m) - yhat_p(j,k)),
 * its signed area is
 *   A_p(w;j,m) = m^{-1} sum_{k=1}^{m} w(k/m) T_{j,m}(k/m),
 * and the three estimators of sigma_p^2 = lim n Var(ytilde_p(n)) are
 *   A_p(w;b,m) = b^{-1} sum_j A_p(w;j,m)^2                       (STS area)
 *   N_p(b,m)   = (b-1)^{-1} m sum_j (yhat_p(j,m)-ytilde_p(n))^2  (NBQ)
 *   V_p(w;b,m) = [b A_p(w;b,m) + (b-1) N_p(b,m)] / (2b-1)        (combined)
 * where ytilde_p(n) is the full-sample empirical p-quantile over all n = b*m
 * observations. The first two have limiting chi-square laws on b and b-1
 * degrees of freedom and are asymptotically independent, so the combined
 * estimator carries 2b-1 degrees of freedom and is about sqrt(2) less variable
 * than either component.
 *
 * WEIGHT FUNCTION. The requirement on w is that int_0^1 w(t)B(t)dt be standard
 * normal for a standard Brownian bridge B; for a constant w = c that variance is
 * c^2/12, so c = sqrt(12) is the normalizing choice and any other constant
 * rescales every area and A_p by (c/sqrt(12))^2. Only constant weights are
 * supported, as in the reference.
 *
 * COST. The prefix quantiles yhat_p(j,k) are exact order statistics, not a
 * running approximation: a Fenwick tree over the within-batch ranks is advanced
 * one observation at a time and searched by binary lifting, so the whole
 * function is O(b m log m). MATLAB advances the b trees in lockstep to keep its
 * k loop vectorized; here the batches are simply looped, which is the same
 * arithmetic in a different order.
 *
 * A SINGLE BATCH carries no between-batch degrees of freedom, so N_p and V_p are
 * NaN at b = 1 while areas and bqe stay valid; sim_firquest pools them across
 * replications instead of reading them per replication.
 *
 * Reference: C. Alexopoulos, D. Goldsman, A. Lolos, K. D. Dingec, J. R. Wilson,
 * "Steady-State Quantile Estimation Using Standardized Time Series", 2020/2023;
 * A. Lolos et al., Proc. Winter Simulation Conference, 2023, theorems 1-3.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

#include "line/api/sim/sim_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

/** Batched-quantile statistics of one sample path. */
template <class T>
struct StsQuantileStats {
    std::vector<T> areas;  ///< b signed STS areas A_p(w;j,m)
    std::vector<T> bqe;    ///< b batched quantile estimators yhat_p(j,m)
    T quantile;            ///< Full-sample empirical p-quantile ytilde_p(n)
    T Ap;                  ///< Batched STS area estimator A_p(w;b,m)
    T Np;                  ///< NBQ variance-parameter estimator N_p(b,m), NaN at b = 1
    T Vp;                  ///< Combined variance-parameter estimator, NaN at b = 1
    std::size_t b = 0;     ///< Batch count
    std::size_t m = 0;     ///< Batch size
    std::size_t n = 0;     ///< Number of observations used, b*m
};

namespace detail {

/** Fenwick prefix-count tree over the ranks 1..m of one batch. */
class RankTree {
public:
    explicit RankTree(std::size_t m) : f_(m + 1, 0), m_(m) {
        step_ = 1;
        while (step_ * 2 <= m) step_ *= 2;
    }

    void add(std::size_t rank) {
        for (std::size_t i = rank; i <= m_; i += i & (~i + 1)) ++f_[i];
    }

    /**
     * Zero-based index of the L-th smallest rank inserted so far, i.e. MATLAB's
     * pos + 1 with pos the last position whose prefix count stays below L.
     */
    std::size_t select(std::size_t L) const {
        std::size_t pos = 0, rem = L;
        for (std::size_t step = step_; step > 0; step >>= 1) {
            const std::size_t cand = pos + step;
            if (cand <= m_ && f_[cand] < rem) {
                rem -= f_[cand];
                pos = cand;
            }
        }
        return pos;
    }

private:
    std::vector<std::size_t> f_;
    std::size_t m_;
    std::size_t step_;
};

}  // namespace detail

/**
 * @param Y      exactly b*m finite observations, in sample-path order
 * @param b      batch count, positive
 * @param m      batch size, positive
 * @param p      quantile order in (0,1)
 * @param weight constant STS weight function, sqrt(12) by default
 */
template <class T>
StsQuantileStats<T> sim_sts_quantile_areas(const std::vector<T>& Y, std::size_t b, std::size_t m,
                                           double p, double weight = std::sqrt(12.0)) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_sts_quantile_areas: the areas carry the irrational normalizing weight "
                  "sqrt(12)/(m sqrt(m)), so exact arithmetic is refused");
    if (b < 1) throw InputError("sim_sts_quantile_areas: the batch count b must be positive");
    if (m < 1) throw InputError("sim_sts_quantile_areas: the batch size m must be positive");
    if (!(p > 0.0) || !(p < 1.0))
        throw InputError("sim_sts_quantile_areas: p must be a real scalar in (0,1)");
    if (weight == 0.0)
        throw InputError("sim_sts_quantile_areas: weight must be a nonzero real scalar");

    const std::size_t n = b * m;
    if (Y.size() != n)
        throw InputError("sim_sts_quantile_areas: Y must hold exactly b*m observations");
    for (std::size_t i = 0; i < n; ++i)
        if (!detail::num_isfinite(Y[i]))
            throw InputError("sim_sts_quantile_areas: the sample path must be finite");

    StsQuantileStats<T> st;
    st.b = b;
    st.m = m;
    st.n = n;
    st.areas.assign(b, num_traits<T>::from_int(0));
    st.bqe.assign(b, num_traits<T>::from_int(0));

    const T wgt = num_traits<T>::from_double(weight);
    const T md = num_traits<T>::from_int(static_cast<long>(m));
    const T scale = T(wgt / T(md * detail::num_sqrt(md)));
    const std::size_t bqeIdx = static_cast<std::size_t>(
                                   std::ceil(static_cast<double>(m) * p)) - 1;

    std::vector<std::size_t> ord(m), rnk(m);
    for (std::size_t j = 0; j < b; ++j) {
        const T* col = &Y[j * m];

        std::iota(ord.begin(), ord.end(), static_cast<std::size_t>(0));
        // stable, so ties keep sample-path order exactly as MATLAB's sort does
        std::stable_sort(ord.begin(), ord.end(),
                         [col](std::size_t a, std::size_t c) { return col[a] < col[c]; });
        std::vector<T> sorted(m);
        for (std::size_t r = 0; r < m; ++r) {
            sorted[r] = col[ord[r]];
            rnk[ord[r]] = r + 1;
        }

        st.bqe[j] = sorted[bqeIdx];

        detail::RankTree tree(m);
        T acc = num_traits<T>::from_int(0);
        for (std::size_t k = 1; k <= m; ++k) {
            tree.add(rnk[k - 1]);
            const std::size_t L =
                static_cast<std::size_t>(std::ceil(p * static_cast<double>(k)));
            const T qk = sorted[tree.select(L)];
            acc += T(num_traits<T>::from_int(static_cast<long>(k)) * T(st.bqe[j] - qk));
        }
        st.areas[j] = T(scale * acc);
    }

    std::vector<T> all(Y);
    std::sort(all.begin(), all.end());
    st.quantile = all[static_cast<std::size_t>(std::ceil(static_cast<double>(n) * p)) - 1];

    T sumSq = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < b; ++j) sumSq += T(st.areas[j] * st.areas[j]);
    st.Ap = T(sumSq / num_traits<T>::from_int(static_cast<long>(b)));

    if (b >= 2) {
        T sd = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < b; ++j) {
            const T d = T(st.bqe[j] - st.quantile);
            sd += T(d * d);
        }
        st.Np = T(md * sd / num_traits<T>::from_int(static_cast<long>(b - 1)));
        st.Vp = T((num_traits<T>::from_int(static_cast<long>(b)) * st.Ap +
                   num_traits<T>::from_int(static_cast<long>(b - 1)) * st.Np) /
                  num_traits<T>::from_int(static_cast<long>(2 * b - 1)));
    } else {
        st.Np = detail::num_nan<T>();
        st.Vp = detail::num_nan<T>();
    }
    return st;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_STS_QUANTILE_AREAS_H
