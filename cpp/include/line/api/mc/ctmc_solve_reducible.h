/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SOLVE_REDUCIBLE_H
#define LINE_API_MC_CTMC_SOLVE_REDUCIBLE_H

/**
 * Limiting distribution of a CTMC whose generator may be reducible.
 *
 * Templated port of matlab/src/api/mc/ctmc_solve_reducible.m and
 * jar/src/main/java/jline/api/mc/Ctmc_solve_reducible.java: the generator is
 * uniformized and handed to dtmc_solve_reducible, which does all the work. The
 * uniformized chain has the same strongly connected components and the same
 * limiting distribution as the CTMC, so nothing is lost by the detour, and the
 * component logic lives in exactly one place.
 *
 * EXACT. Uniformization at the deterministic rate (21/20) max|Q| is a field
 * operation and the reducible solver carries no tolerance either, so a
 * generator with rational rates yields the exact limiting distribution. This
 * is the routine ctmc_courtois uses for its diagonal blocks, which is how the
 * exactness reaches the aggregation methods' microprobabilities.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * @param Q   generator
 * @param pi0 initial distribution; empty to let the routine pick one
 * @param zeroColTol column-sum threshold below which a state counts as unreachable
 */
template <class T>
ReducibleResult<T> ctmc_solve_reducible(const Matrix<T>& Q, const std::vector<T>& pi0,
                                        double zeroColTol = 1e-12) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_solve_reducible: generator is not square");
    return dtmc_solve_reducible(ctmc_randomization(Q).P, pi0, zeroColTol);
}

/** Overload without an initial vector. */
template <class T>
ReducibleResult<T> ctmc_solve_reducible(const Matrix<T>& Q, double zeroColTol = 1e-12) {
    return ctmc_solve_reducible(Q, std::vector<T>(), zeroColTol);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SOLVE_REDUCIBLE_H
