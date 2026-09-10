/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_QR_LD_H
#define LINE_API_MAPQN_MAPQN_BND_QR_LD_H

/**
 * Quadratic-reduction bound on a marginal of a load-dependent MAP queueing
 * network.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_qr_ld.m (ground truth) and of
 * jar/src/main/java/jline/api/mapqn/Mapqn_bnd_qr_ld.java.
 *
 * What it computes: the exact stationary distribution of a MAP queueing
 * network is intractable, so the quadratic reduction keeps only the PAIRWISE
 * joint distribution p2(j,nj,k,i,ni,h) -- the probability that queue j holds
 * nj jobs in phase k while queue i holds ni jobs in phase h -- and imposes
 * every linear relation the true distribution must satisfy: normalization
 * (ONE), the states that cannot occur (ZERO1/2/3), symmetry of the pair,
 * consistency of the pairwise law with its own marginal (MARGINALS), Little's
 * law (THM1, THM1c), the second moment of the population (PC2), phase balance
 * (THM2), population flow balance (THM3a, THM3b), queue balance (QBAL), the
 * order-1 correlation cuts (COR1a, COR1b) and the QMIN inequality (THM4). The
 * true distribution is feasible for that polytope, so maximizing or minimizing
 * one coordinate over it BOUNDS the corresponding true probability. The bound
 * is a relaxation, not an approximation: it is valid, not merely close.
 *
 * The objective coordinate is the diagonal entry
 * p2(objective_queue, objective_n, objective_phase, same, same, same), which
 * is P(queue j holds n jobs in phase k). Summing the bound over n >= 1 and k
 * is how the reference obtains a utilization bound.
 *
 * Arithmetic: the whole assembly is +, -, * and / on the model data, and the
 * LP is solved by line::lp::simplex_solve, which uses Bland's rule and no
 * tolerance. At T = line::Rational the returned bound is therefore the EXACT
 * optimum of the exact polytope, with no LP tolerance at all -- unlike the
 * MATLAB reference, which reaches it with interior-point linprog and lands a
 * few digits short (its own comment records ~1e-7 residual error on the
 * paper's BAS instance), and unlike the JAR, which uses a double-precision
 * Apache Commons simplex.
 *
 * Variable bounds: 0 <= p2 <= 1 with the ZERO states pinned to 0 are passed to
 * the solver as bounds. They do NOT become rows the way the JAR needs them to,
 * because Apache SimplexSolver does not box variables and the maximization is
 * unbounded without them; LpModel boxes natively, substitutes out the pinned
 * variables entirely, and materializes at most one row per genuinely
 * upper-bounded variable.
 *
 * Cost: the model has ((N+1) sum_i K(i))^2 variables, so it grows as the
 * fourth power of the population. The tableau here is dense, which confines
 * the port to small and medium instances; the reference's large blocking
 * instances need a sparse revised simplex.
 */

#include "line/api/mapqn/mapqn_params.h"
#include "line/api/mapqn/mapqn_qr_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/simplex.h"

namespace line {
namespace mapqn {

/**
 * Bound P(queue objective_queue holds objective_n jobs in phase
 * objective_phase), in the requested direction.
 *
 * @param p               network parameters; queues and phases are 0-based
 * @param objective_queue queue index, 0..M-1
 * @param objective_phase phase index, 0..K(objective_queue)-1
 * @param objective_n     population level, 0..N
 * @param sense           Max for an upper bound, Min for a lower bound
 */
template <class T>
MapqnQrResult<T> mapqn_bnd_qr_ld(const MapqnParams<T>& p, int objective_queue, int objective_phase,
                                 int objective_n, MapqnSense sense = MapqnSense::Max) {
    p.validate();
    if (p.N < 1) throw InputError("mapqn_bnd_qr_ld: N must be at least 1");
    detail::qr_check_objective(p, objective_queue, objective_phase, objective_n);

    const P2Index idx(p.M, p.N, p.K);
    lp::LpModel<T> m(idx.num_vars());
    const T one = num_traits<T>::from_int(1);
    for (std::size_t j = 0; j < idx.num_vars(); ++j) m.set_bounds(j, T(), one);

    // Families, in the order the reference emits them.
    const std::vector<char> is_zero = qr_zero_bounds(p, idx, m);
    qr_one(p, idx, m);
    qr_symmetry(p, idx, m, is_zero);
    qr_marginals(p, idx, m);
    qr_thm1(p, idx, m);
    qr_thm1c(p, idx, m);
    qr_pc2(p, idx, m);
    qr_thm2(p, idx, m);
    qr_thm3a(p, idx, m);
    qr_thm3b(p, idx, m);
    qr_qbal(p, idx, m);
    qr_cor1a(p, idx, m);
    qr_cor1b(p, idx, m);
    qr_thm4(p, idx, m);

    return detail::qr_finish(p, idx, m, objective_queue, objective_phase, objective_n, sense);
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_BND_QR_LD_H
