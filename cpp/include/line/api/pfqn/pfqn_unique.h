/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_UNIQUE_H
#define LINE_API_PFQN_UNIQUE_H

/**
 * Station consolidation and its inverse: merge identical demand rows into one
 * station with a multiplicity, and expand per-station results back.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_unique.m and pfqn_expand.m, plus
 * the load-dependent rate shift pfqn_mushift.m.
 *
 * One deliberate divergence from MATLAB, the same one already taken in
 * pfqn_recal: MATLAB compares rows within GlobalConstants.Zero() (1e-14), so it
 * merges stations whose demands only nearly agree. A tolerance-based merge
 * perturbs the normalizing constant, which is unacceptable for an
 * exact-capable algorithm and has no meaning at all in the rational field. The
 * port merges on exact equality. Callers that want tolerant merging should
 * round their demands before calling, where the rounding is visible.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct UniqueResult {
    Matrix<T> L;                    ///< consolidated demands (Mu x R)
    Matrix<T> mu;                   ///< consolidated LD rates, empty if none
    Matrix<T> gamma;                ///< consolidated CD scalings, empty if none
    std::vector<int> mi;            ///< multiplicity of each consolidated station
    std::vector<std::size_t> mapping;  ///< original station -> consolidated index
};

/** Merge stations whose (L, mu, gamma) rows are exactly equal. */
template <class T>
UniqueResult<T> pfqn_unique(const Matrix<T>& L, const Matrix<T>& mu, const Matrix<T>& gamma) {
    const std::size_t M = L.rows(), R = L.cols();
    if (!mu.empty() && mu.rows() != M) throw InputError("pfqn_unique: mu has the wrong row count");
    if (!gamma.empty() && gamma.rows() != M)
        throw InputError("pfqn_unique: gamma has the wrong row count");

    UniqueResult<T> r;
    r.mapping.assign(M, 0);
    std::vector<std::size_t> reps;  // representative original index per group

    const auto rows_equal = [&](std::size_t a, std::size_t b) {
        for (std::size_t j = 0; j < R; ++j)
            if (L(a, j) != L(b, j)) return false;
        for (std::size_t j = 0; j < mu.cols(); ++j)
            if (mu(a, j) != mu(b, j)) return false;
        for (std::size_t j = 0; j < gamma.cols(); ++j)
            if (gamma(a, j) != gamma(b, j)) return false;
        return true;
    };

    for (std::size_t i = 0; i < M; ++i) {
        bool merged = false;
        for (std::size_t g = 0; g < reps.size(); ++g) {
            if (rows_equal(i, reps[g])) {
                r.mapping[i] = g;
                r.mi[g] += 1;
                merged = true;
                break;
            }
        }
        if (!merged) {
            r.mapping[i] = reps.size();
            reps.push_back(i);
            r.mi.push_back(1);
        }
    }

    r.L = Matrix<T>(reps.size(), R);
    for (std::size_t g = 0; g < reps.size(); ++g)
        for (std::size_t j = 0; j < R; ++j) r.L(g, j) = L(reps[g], j);
    if (!mu.empty()) {
        r.mu = Matrix<T>(reps.size(), mu.cols());
        for (std::size_t g = 0; g < reps.size(); ++g)
            for (std::size_t j = 0; j < mu.cols(); ++j) r.mu(g, j) = mu(reps[g], j);
    }
    if (!gamma.empty()) {
        r.gamma = Matrix<T>(reps.size(), gamma.cols());
        for (std::size_t g = 0; g < reps.size(); ++g)
            for (std::size_t j = 0; j < gamma.cols(); ++j) r.gamma(g, j) = gamma(reps[g], j);
    }
    return r;
}

template <class T>
UniqueResult<T> pfqn_unique(const Matrix<T>& L) {
    return pfqn_unique(L, Matrix<T>(), Matrix<T>());
}

/**
 * Fold a caller-supplied multiplicity vector along a consolidation mapping.
 *
 * Port of pfqn_combine_mi in
 * jar/src/main/java/jline/api/pfqn/Pfqn_replicas.java, the third member of the
 * replica trio alongside pfqn_unique and pfqn_expand. MATLAB has no
 * counterpart.
 *
 * When the caller already carries its own per-station multiplicities mi (a
 * station standing for mi(i) identical replicas) and pfqn_unique then merges
 * further stations, the two multiplicities must COMPOSE: the consolidated
 * station g represents sum over the original stations mapped to g of mi(i).
 * The result is therefore a SUM along the mapping, not a count of the group,
 * which is what makes it different from the group sizes pfqn_unique itself
 * returns.
 *
 * Arithmetic: EXACT-CAPABLE. Integer addition only.
 *
 * @param mi       (M) multiplicity of each original station
 * @param mapping  (M) original station -> consolidated index, as pfqn_unique
 *                 returns it
 * @param M_unique number of consolidated stations
 * @return (M_unique) combined multiplicities
 */
inline std::vector<int> pfqn_combine_mi(const std::vector<int>& mi,
                                        const std::vector<std::size_t>& mapping,
                                        std::size_t M_unique) {
    if (mi.size() != mapping.size())
        throw InputError("pfqn_combine_mi: mi and mapping have different lengths");
    std::vector<int> combined(M_unique, 0);
    for (std::size_t i = 0; i < mapping.size(); ++i) {
        if (mapping[i] >= M_unique)
            throw InputError("pfqn_combine_mi: mapping index out of range");
        combined[mapping[i]] += mi[i];
    }
    return combined;
}

// pfqn_expand and pfqn_mushift live in their own headers, matching the 1:1
// file-per-function convention: see pfqn_expand.h and pfqn_mushift.h.

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_UNIQUE_H
