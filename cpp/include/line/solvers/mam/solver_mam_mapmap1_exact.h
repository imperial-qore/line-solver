/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_MAPMAP1_EXACT_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_MAPMAP1_EXACT_H

/**
 * Port of `solver_mam_mapmap1_exact.m`: the exact fast path SolverMAM tries
 * BEFORE anything else, for a single-class open Source -> FCFS Queue -> Sink
 * model whose arrival or service is a genuinely CORRELATED MAP.
 *
 * Why it comes first. The decomposition methods approximate exactly this case:
 * `dec.source` hands the service to `MMAPPH1FCFS` as a renewal phase type,
 * which keeps the service-time marginal and discards the correlation between
 * consecutive services. When the process is renewal the two agree and the fast
 * path stands down (the `is_renewal_map` test below); when it is not, the fast
 * path returns the exact matrix-geometric answer instead.
 *
 * The reference reaches `Q_CT_MAP_MAP_1`; this port reaches `qbd_mapmap1`,
 * which solves the same level-independent QBD from the port's own machinery
 * (see `api/mam/qbd_mapmap1.h` for the measured agreement).
 *
 * `ok = false` means "not an exactly-solvable single MAP/MAP/1 queue", and the
 * caller falls through to the decomposition. Every early return in the
 * reference is reproduced, including the stability test lambda < mu: an
 * unstable or degenerate model is left to the fallback rather than answered
 * with a divergent geometric series.
 */

#include <cmath>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mam {

/** Result of the fast path; `ok` false means the model is not in its regime. */
template <class T>
struct MapMap1Exact {
    bool ok = false;
    mva::MvaSolution<T> sol;
};

template <class T>
MapMap1Exact<T> solver_mam_mapmap1_exact(const qn::NetworkStruct<T>& L) {
    MapMap1Exact<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        // The QBD's cyclic reduction is tolerance-terminated; under exact
        // arithmetic the fast path simply does not claim the model, and the
        // dispatch below refuses by name.
        return out;
    } else {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;
    if (K != 1) return out;
    if (!std::isinf(L.classes[0].population)) return out;

    std::size_t src = 0, q = 0;
    std::size_t nsrc = 0, nq = 0;
    for (std::size_t i = 1; i <= M; ++i) {
        if (L.stations[i - 1].sched == SchedStrategy::EXT) {
            src = i;
            ++nsrc;
        } else if (L.stations[i - 1].sched == SchedStrategy::FCFS) {
            q = i;
            ++nq;
        }
    }
    if (nsrc != 1 || nq != 1) return out;
    if (L.stations[q - 1].nservers != 1.0) return out;
    // The arrival the queue sees equals the source MAP only when nothing sits
    // between them, so the two must be the model's only stations.
    if (M != 2) return out;
    if (L.disabled[src - 1][0] || L.disabled[q - 1][0]) return out;

    const Map<T> arv = lang::dist_to_map(L.service[src - 1][0]);
    const Map<T> svc = lang::dist_to_map(L.service[q - 1][0]);
    // A RAP or ME pair is a valid point process but not a MAP, and the QBD
    // takes Markovian blocks; leave those to the fallback.
    if (!basic_detail::is_markovian_map(arv) || !basic_detail::is_markovian_map(svc)) return out;
    // Renewal on both sides means the decomposition is already exact.
    if (basic_detail::is_renewal_map(arv) && basic_detail::is_renewal_map(svc)) return out;

    const T lambda = map_lambda(arv);
    const T mu = map_lambda(svc);
    if (!(lambda < mu)) return out;

    const QbdMapMap1Result<T> r = qbd_mapmap1(arv, svc);

    mva::MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.Tp(src - 1, 0) = lambda;
    s.Tp(q - 1, 0) = lambda;
    s.Q(q - 1, 0) = r.QN;
    s.U(q - 1, 0) = T(lambda / mu);
    s.R(q - 1, 0) = T(r.QN / lambda);
    s.C[0] = s.R(q - 1, 0);
    s.X[0] = lambda;
    s.iter = 1;
    out.ok = true;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_MAPMAP1_EXACT_H
