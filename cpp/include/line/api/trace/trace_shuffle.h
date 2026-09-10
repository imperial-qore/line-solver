/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_SHUFFLE_H
#define LINE_API_TRACE_TRACE_SHUFFLE_H

/**
 * Random permutation of a trace.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_shuffle.m.
 *
 * Shuffling destroys the autocorrelation while leaving every marginal moment
 * untouched, which is what makes it the null model for a correlation test: any
 * statistic that moves under a shuffle is reading the ORDER of the trace, not
 * its distribution. Callers that want an uncorrelated trace with the same
 * marginal use this rather than resampling, because resampling would perturb
 * the empirical moments as well.
 */

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace trace {

/** A uniformly random permutation of the samples, drawn with the given engine. */
template <class T, class Gen>
std::vector<T> trace_shuffle(const std::vector<T>& S, Gen& gen) {
    if (S.empty()) throw InputError("trace_shuffle: the trace is empty");
    // Fisher-Yates over the index set, so T needs no swap beyond a copy
    std::vector<std::size_t> idx(S.size());
    std::iota(idx.begin(), idx.end(), static_cast<std::size_t>(0));
    for (std::size_t i = S.size(); i > 1; --i) {
        std::uniform_int_distribution<std::size_t> pick(0, i - 1);
        std::swap(idx[i - 1], idx[pick(gen)]);
    }
    std::vector<T> out(S.size());
    for (std::size_t i = 0; i < S.size(); ++i) out[i] = S[idx[i]];
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_SHUFFLE_H
