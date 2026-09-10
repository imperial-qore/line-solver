/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_BND_QR_DELAY_H
#define LINE_API_MAPQN_MAPQN_BND_QR_DELAY_H

/**
 * Quadratic-reduction bound for a MAP queueing network with a delay (think
 * time) station.
 *
 * Templated port of matlab/lib/qrf/mapqn_bnd_qr_delay.m (ground truth) and of
 * jar/src/main/java/jline/api/mapqn/Mapqn_bnd_qr_delay.java.
 *
 * The polytope is exactly the load-dependent one of mapqn_bnd_qr_ld plus the
 * XZ family:
 *   sum_{n>=1,k} n p2(M,n,k,M,n,k) = (Z/D1) sum_{k,n>=1} p2(1,n,k,1,n,k)
 * which is Little's law across the think time: the mean population at the
 * delay station (queue M) is Z times the throughput of queue 1, whose service
 * demand is D1. The delay station is made infinite-server by the caller
 * through the load-dependent scalings, alpha(M,n) = n; nothing in the
 * assembly special-cases it.
 *
 * Everything the mapqn_bnd_qr_ld header says about exactness, bound handling
 * and cost applies here unchanged. On the ordering of families: the MATLAB
 * delay file emits PC2 before THM1 and the ld file emits it after. Row order
 * has no effect on the optimum -- it is the same polytope either way -- so the
 * shared assembly here uses one order. The families themselves are identical:
 * despite what a summary of this domain might suggest, the MATLAB delay
 * reference does NOT drop THM1 or THM1c, it keeps both and adds XZ.
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
 * objective_phase) for a network with a delay station at queue M.
 *
 * p.Z is the think time and p.D1 the service demand at queue 1; both enter
 * only through the XZ family.
 */
template <class T>
MapqnQrResult<T> mapqn_bnd_qr_delay(const MapqnParams<T>& p, int objective_queue,
                                    int objective_phase, int objective_n,
                                    MapqnSense sense = MapqnSense::Max) {
    p.validate();
    if (p.N < 1) throw InputError("mapqn_bnd_qr_delay: N must be at least 1");
    if (p.M < 2) throw InputError("mapqn_bnd_qr_delay: the delay model needs at least two queues");
    detail::qr_check_objective(p, objective_queue, objective_phase, objective_n);

    const P2Index idx(p.M, p.N, p.K);
    lp::LpModel<T> m(idx.num_vars());
    const T one = num_traits<T>::from_int(1);
    for (std::size_t j = 0; j < idx.num_vars(); ++j) m.set_bounds(j, T(), one);

    const std::vector<char> is_zero = qr_zero_bounds(p, idx, m);
    qr_one(p, idx, m);
    qr_symmetry(p, idx, m, is_zero);
    qr_marginals(p, idx, m);
    qr_pc2(p, idx, m);
    qr_thm1(p, idx, m);
    qr_thm1c(p, idx, m);
    qr_xz(p, idx, m);  // the only family the load-dependent model does not have
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

#endif  // LINE_API_MAPQN_MAPQN_BND_QR_DELAY_H
