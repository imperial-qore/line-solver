/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_LJD_LINEARIZE_H
#define LINE_API_FES_LJD_LINEARIZE_H

/**
 * Linearized index of a per-class population vector, for Limited Joint
 * Dependence (LJD) tables.
 *
 * Port of matlab/src/api/ljd/ljd_linearize.m,
 *
 *   idx = 1 + n1 + n2 (N1+1) + n3 (N1+1)(N2+1) + ...
 *
 * with every n_k clamped to its cutoff N_k. The result is 1-BASED, exactly as
 * in MATLAB, because the FES scaling tables it addresses are transcribed from
 * MATLAB indices; the callers in this port subtract one when they index a
 * std::vector. It is otherwise the mixed-radix index of line/util/population.h
 * (pop_index over plane_sizes), with the clamp added.
 *
 * The header lives under api/fes because the FES functions are its only
 * consumers in this port; when the ljd domain itself is ported it should move
 * to api/ljd unchanged.
 */

#include <cstddef>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace fes {

/**
 * @param nvec    per-class populations
 * @param cutoffs per-class cutoffs, same length as nvec
 * @return the 1-based linearized index
 */
inline std::size_t ljd_linearize(const std::vector<int>& nvec, const std::vector<int>& cutoffs) {
    if (nvec.size() != cutoffs.size())
        throw InputError("ljd_linearize: population and cutoff vectors have different lengths");
    std::size_t idx = 1;  // 1-indexed, as in MATLAB
    std::size_t multiplier = 1;
    for (std::size_t k = 0; k < nvec.size(); ++k) {
        const int nk = nvec[k] < cutoffs[k] ? nvec[k] : cutoffs[k];  // clamp to cutoff
        if (nk > 0) idx += static_cast<std::size_t>(nk) * multiplier;
        multiplier *= static_cast<std::size_t>(cutoffs[k] + 1);
    }
    return idx;
}

/**
 * Inverse of `ljd_linearize`: the population vector behind an index.
 *
 * The forward map is a mixed-radix numeral with class k in radix (Nk+1), so the
 * inverse is the digit-by-digit division that reads it back. It is what lets a
 * caller walk a tabulated dependence in index order and still know which
 * population each entry belongs to.
 *
 * @param idx     the 1-based linearized index
 * @param cutoffs per-class cutoffs
 * @return the per-class population vector
 */
inline std::vector<int> ljd_delinearize(std::size_t idx, const std::vector<int>& cutoffs) {
    std::vector<int> nvec(cutoffs.size(), 0);
    std::size_t rem = idx - 1;
    for (std::size_t k = 0; k < cutoffs.size(); ++k) {
        const std::size_t radix = static_cast<std::size_t>(cutoffs[k] + 1);
        nvec[k] = static_cast<int>(rem % radix);
        rem /= radix;
    }
    return nvec;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_LJD_LINEARIZE_H
