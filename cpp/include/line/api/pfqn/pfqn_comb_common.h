/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_COMB_COMMON_H
#define LINE_API_PFQN_COMB_COMMON_H

/**
 * Integer-composition enumeration shared by the CoMoM and MVAC ports.
 *
 * Templated-free helpers mirroring matlab/src/util/multichoose.m and
 * matlab/src/util/matchrow.m. Both are pure index arithmetic on integers, so
 * they carry no number type at all and are usable from any instantiation.
 *
 * ORDERING MATTERS. MATLAB's multichoose(n,k) recurses as
 *
 *   for i = 0:k,  rows = [ i , multichoose(n-1, k-i) ]
 *
 * so the first component ascends slowest. Several callers (pfqn_comom's basis
 * layout, pfqn_mvac's multiplicity list) index into the result by position
 * rather than by content, and a different but equally valid enumeration order
 * would silently permute their bases. The recursion is reproduced verbatim
 * rather than replaced by a "nicer" odometer.
 *
 * `line::multichoose` in util/population.h returns the COUNT C(n+k-1,k);
 * `multichoose_rows` here returns the compositions themselves. The names are
 * kept distinct so that a call site cannot pick up the wrong one.
 */

#include <cstddef>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace pfqn {

/**
 * All n-vectors of nonnegative integers summing to k, in MATLAB
 * multichoose(n,k) order.
 */
inline std::vector<std::vector<int>> multichoose_rows(int n, int k) {
    if (n < 1) throw InputError("multichoose_rows: at least one component is required");
    if (k < 0) throw InputError("multichoose_rows: negative total");
    std::vector<std::vector<int>> out;
    if (n == 1) {
        out.push_back(std::vector<int>(1, k));
        return out;
    }
    if (k == 0) {
        out.push_back(std::vector<int>(static_cast<std::size_t>(n), 0));
        return out;
    }
    for (int i = 0; i <= k; ++i) {
        const std::vector<std::vector<int>> w = multichoose_rows(n - 1, k - i);
        for (std::size_t j = 0; j < w.size(); ++j) {
            std::vector<int> row;
            row.reserve(static_cast<std::size_t>(n));
            row.push_back(i);
            row.insert(row.end(), w[j].begin(), w[j].end());
            out.push_back(row);
        }
    }
    return out;
}

/**
 * Position of `row` in `rows`, or -1 when absent. MATLAB's matchrow checks the
 * LAST row first and otherwise returns the first match; with the distinct row
 * sets used here the two rules coincide, and the first-match rule is used.
 */
inline int matchrow(const std::vector<std::vector<int>>& rows, const std::vector<int>& row) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].size() != row.size()) return -2;
        bool eq = true;
        for (std::size_t j = 0; j < row.size(); ++j)
            if (rows[i][j] != row[j]) {
                eq = false;
                break;
            }
        if (eq) return static_cast<int>(i);
    }
    return -1;
}

/**
 * MATLAB's sortbynnzpos: a stable bubble sort putting the rows with FEWER
 * nonzeros first and, among rows with equally many, the row whose leftmost
 * differing entry is nonzero first. Reproduced exactly, including the O(n^2)
 * shape, because the resulting order is the basis layout that pfqn_comom and
 * pfqn_procomom index into.
 */
inline void sort_by_nnz_pos(std::vector<std::vector<int>>& I) {
    const auto nnz = [](const std::vector<int>& v) {
        int c = 0;
        for (int x : v)
            if (x != 0) ++c;
        return c;
    };
    // returns true when i1 must come AFTER i2 (MATLAB nnzcmp(i1,i2) == 1)
    const auto after = [&](const std::vector<int>& i1, const std::vector<int>& i2) {
        const int n1 = nnz(i1), n2 = nnz(i2);
        if (n1 > n2) return true;
        if (n1 < n2) return false;
        for (std::size_t j = 0; j < i1.size(); ++j) {
            if (i1[j] == 0 && i2[j] > 0) return true;
            if (i1[j] > 0 && i2[j] == 0) return false;
        }
        return false;
    };
    for (std::size_t i = 0; i + 1 < I.size(); ++i)
        for (std::size_t j = i + 1; j < I.size(); ++j)
            if (after(I[i], I[j])) std::swap(I[i], I[j]);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_COMB_COMMON_H
