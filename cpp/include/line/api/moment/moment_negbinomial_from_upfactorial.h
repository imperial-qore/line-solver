/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_NEGBINOMIAL_FROM_UPFACTORIAL_H
#define LINE_API_MOMENT_MOMENT_NEGBINOMIAL_FROM_UPFACTORIAL_H

/**
 * Negative-binomial moments from rising-factorial moments.
 *
 * Templated port of matlab/src/api/moment/moment_negbinomial_from_upfactorial.m. Every operation is integer
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

/** bm_i = fp_i / i!. */
template <class T>
std::vector<T> moment_negbinomial_from_upfactorial(const std::vector<T>& fp) {
    std::vector<T> bm(fp.size());
    for (std::size_t i = 0; i < fp.size(); ++i)
        bm[i] = fp[i] / num_factorial<T>(static_cast<unsigned>(i));
    return bm;
}

}  // namespace moment
}  // namespace line

#endif
