/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_UPFACTORIAL_FROM_NEGBINOMIAL_H
#define LINE_API_MOMENT_MOMENT_UPFACTORIAL_FROM_NEGBINOMIAL_H

/**
 * Rising-factorial moments from negative-binomial moments.
 *
 * Templated port of matlab/src/api/moment/moment_upfactorial_from_negbinomial.m. Every operation is integer
 * or rational, so the exact instantiation returns the transform with no
 * rounding at all. That matters more here than almost anywhere else in the API:
 * the alternating binomial sums of the moment conversions cancel
 * catastrophically in double arithmetic once the order grows.
 */

#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace moment {

/** fp_i = i! bm_i. */
template <class T>
std::vector<T> moment_upfactorial_from_negbinomial(const std::vector<T>& bm) {
    std::vector<T> fp(bm.size());
    for (std::size_t i = 0; i < bm.size(); ++i)
        fp[i] = num_factorial<T>(static_cast<unsigned>(i)) * bm[i];
    return fp;
}

}  // namespace moment
}  // namespace line

#endif
