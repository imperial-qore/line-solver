/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CUB_EVALS_H
#define LINE_API_PFQN_CUB_EVALS_H

/**
 * Integrand-evaluation count of pfqn_cub, and the budget pfqn_nc prices it
 * against.
 *
 * Port of matlab/src/api/pfqn/pfqn_cub_evals.m.
 *
 * The Grundmann-Moeller rule of degree `order` on the (M-1)-simplex evaluates
 *   sum_{d = 0..order} C(M-1+2d, M-1)
 * points, and a non-zero think time makes pfqn_cub repeat the whole rule at
 * each of its v-quadrature steps (a uniform grid of CUB_V_STEPS points, which
 * must match the outer McKenna-Mitra integral in pfqn_cub.h).
 *
 * NOTE the two DIFFERENT binomials in play. The cost model pfqn_nc uses to
 * RAISE the order counts C(M + 2d, M-1); the count here, which pfqn_nc then
 * uses to LOWER it again, counts C(M-1+2d, M-1). Both are reproduced as in the
 * reference: they are not the same expression and folding one into the other
 * would change the selected order.
 *
 * Arithmetic: this is a cost model, not a numerical result. It is a plain
 * double count and carries no number type.
 */

#include <cstddef>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** The v-quadrature grid size of pfqn_cub; must match `steps` in pfqn_cub.h. */
constexpr long CUB_V_STEPS = 10000;

/**
 * GlobalConstants.CubMaxEvals: the integrand-evaluation budget above which
 * pfqn_nc lowers the cubature order (and, at order 0, prefers le over cub).
 */
constexpr double CUB_MAX_EVALS = 1e7;

/**
 * @param M     number of queueing stations
 * @param order Grundmann-Moeller degree
 * @param Zsum  total think time; a positive value costs the v-quadrature
 * @return number of integrand evaluations pfqn_cub performs
 */
inline double pfqn_cub_evals(int M, int order, double Zsum) {
    if (M < 1) throw InputError("pfqn_cub_evals: at least one station is required");
    if (order < 0) throw InputError("pfqn_cub_evals: negative cubature order");
    const int n = M - 1;
    double nodes = 0.0;
    for (int d = 0; d <= order; ++d) nodes += nck(n + 2 * d, n);
    // GlobalConstants.FineTol, the reference's threshold for "has a think time"
    if (Zsum >= 1e-8) return nodes * static_cast<double>(CUB_V_STEPS);
    return nodes;
}

/** Zero think time, i.e. the bare simplex rule. */
inline double pfqn_cub_evals(int M, int order) { return pfqn_cub_evals(M, order, 0.0); }

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CUB_EVALS_H
