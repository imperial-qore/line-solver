/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_TRAFFIC_MERGE_CS_H
#define LINE_API_NPFQN_TRAFFIC_MERGE_CS_H

/**
 * Merge of marked arrival flows with class switching.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_merge_cs.m. See
 * npfqn_traffic_merge.h for the shared conventions and for the reference
 * defects of the plain merge.
 *
 * prob((i-1)R + r, s) is the probability that a class-r arrival from flow i
 * switches to class s. Flow i is first re-marked with its own (R x R) block
 * (mmap_mark.m), and the re-marked flows are then superposed class by class.
 *
 * RATE IDENTITY. Marking redistributes an arrival's class without touching
 * (D0, D1), so lambda_s of the merged flow is sum_i sum_r lambda_{i,r} P_i(r,s)
 * exactly, and when every P_i is stochastic the aggregate rate sum_s lambda_s
 * equals the sum of the aggregate rates of the operands. Both are exact
 * identities in the rational instantiation.
 *
 * REFERENCE DEFECTS in npfqn_traffic_merge_cs.m:
 *
 *  1. UNREACHABLE MERGE RULES. The function accepts config.merge but its
 *     switch has only the {'default','super'} case and no otherwise branch, so
 *     any other value (including the 'mixture' and 'interpos' that
 *     npfqn_traffic_merge accepts) silently leaves SMMAP undefined and MATLAB
 *     then raises "Output argument SMMAP not assigned". This port raises
 *     UnsupportedError naming the rule instead.
 *
 *  2. NO COMPRESSION AND NO NORMALIZATION. Unlike npfqn_traffic_merge, the
 *     class-switching variant never looks at config.compress and never calls
 *     mmap_normalize on the result; the normalization it does get comes from
 *     mmap_super itself. Reproduced as written: this port applies no
 *     compression here either.
 *
 *  3. R IS TAKEN FROM prob, NOT FROM THE FLOWS. size(prob, 2) sets the output
 *     class count and the (i-1)R + r indexing assumes size(prob, 1) == n R.
 *     A prob whose row count is not a multiple of the flow count silently
 *     mis-indexes in MATLAB; this port rejects it.
 *
 * ARITHMETIC. Marking and Kronecker sums only, so this is exact-capable and is
 * instantiated at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/npfqn/npfqn_traffic_merge.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/**
 * Merge flows that switch class on departure.
 *
 * @param flows n MMAPs, each carrying R marked classes
 * @param prob  (n R x R), prob((i-1)R + r, s) in MATLAB 1-based terms
 * @param merge merge rule; only Default and Super are ported
 * @return the merged MMAP with R classes
 */
template <class T>
mam::Mmap<T> npfqn_traffic_merge_cs(const std::vector<mam::Mmap<T>>& flows, const Matrix<T>& prob,
                                    Merge merge = Merge::Default) {
    const std::size_t n = flows.size();
    if (n == 0) throw InputError("npfqn_traffic_merge_cs: no flows to merge");
    const std::size_t R = prob.cols();
    if (R == 0 || prob.rows() != n * R)
        throw InputError("npfqn_traffic_merge_cs: prob must be (n R) x R");

    std::vector<mam::Mmap<T>> marked;
    marked.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        Matrix<T> P(R, R);
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) P(r, s) = prob(i * R + r, s);
        marked.push_back(mmap_mark_types(flows[i], P));
    }

    if (n == 1) return marked.front();

    switch (merge) {
        case Merge::Default:
        case Merge::Super:
            break;
        default:
            throw UnsupportedError("npfqn_traffic_merge_cs: only the 'default' and 'super' merge "
                                   "rules are defined for the class-switching merge");
    }

    mam::Mmap<T> s = marked.front();
    for (std::size_t j = 1; j < n; ++j) s = mmap_super_match(s, marked[j]);
    return s;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_TRAFFIC_MERGE_CS_H
