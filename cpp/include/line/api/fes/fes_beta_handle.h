/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_BETA_HANDLE_H
#define LINE_API_FES_BETA_HANDLE_H

/**
 * Wraps a flow-equivalent-server throughput table as a per-class
 * class-dependence function beta_{i,r}(n).
 *
 * Templated port of matlab/src/api/fes/fes_beta_handle.m. The returned
 * callable takes the per-class population n at the FES station and returns
 *
 *   beta_r(n) = X_r(n) |n| / n_r,
 *
 * the DIMENSIONLESS scaling relative to the nominal rate-1 service of the FES
 * station. The |n|/n_r factor cancels the processor-sharing share that the
 * convolution applies (Sauer 1983, eq. (40), with mu_{r,i}(n) = (n_r/|n|)
 * beta_r(n)), leaving the aggregate completing class r at exactly the
 * subnetwork throughput X_r(n). The population is clamped to the cutoffs, so
 * the scaling saturates beyond the tabulated range as the table intends;
 * entries with n_r = 0 are never consulted by the recurrence and are returned
 * as 1.
 *
 * Arithmetic. One table lookup, one multiplication and one division, so this
 * is exact at T = Rational.
 *
 * Deviation from MATLAB, mechanical: the MATLAB handle is a closure over the
 * table; here it is a std::function that owns copies of the table and the
 * cutoffs, so the returned callable outlives its arguments.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/fes/ljd_linearize.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fes {

/** Class-dependence handle: population vector -> per-class beta. */
template <class T>
using FesBetaFun = std::function<std::vector<T>(const std::vector<int>&)>;

/**
 * @param scalingTable (K) linearized per-class throughput vectors, as produced
 *                     by fes_compute_throughputs
 * @param cutoffs      (K) the population vector the table was tabulated on
 */
template <class T>
FesBetaFun<T> fes_beta_handle(const std::vector<std::vector<T>>& scalingTable,
                              const std::vector<int>& cutoffs) {
    const std::vector<std::vector<T>> table = scalingTable;
    const std::vector<int> cut = cutoffs;
    return [table, cut](const std::vector<int>& nin) -> std::vector<T> {
        const std::size_t K = table.size();
        const T one = num_traits<T>::from_int(1);
        std::vector<T> v(K, one);

        // pad with zeros or truncate to the length of the cutoff vector
        std::vector<int> n(cut.size(), 0);
        for (std::size_t k = 0; k < cut.size() && k < nin.size(); ++k) n[k] = nin[k];

        std::vector<int> nClamped(cut.size(), 0);
        for (std::size_t k = 0; k < cut.size(); ++k) {
            int x = n[k] < cut[k] ? n[k] : cut[k];
            nClamped[k] = x < 0 ? 0 : x;
        }
        const std::size_t idx = ljd_linearize(nClamped, cut);

        int tot = 0;
        for (int x : n) tot += x;

        for (std::size_t r = 0; r < K; ++r) {
            if (r >= n.size()) break;
            const std::vector<T>& tbl = table[r];
            if (!tbl.empty() && idx >= 1 && idx <= tbl.size() && n[r] > 0) {
                // FES station rate-sharing rationale: see _kb/03-api-layer.md (cpp port notes: fes)
                v[r] = tbl[idx - 1] * num_traits<T>::from_int(tot) / num_traits<T>::from_int(n[r]);
            }
        }
        return v;
    };
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_BETA_HANDLE_H
