/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_RT_STATIONS_H
#define LINE_API_SN_SN_RT_STATIONS_H

/**
 * Port of matlab/src/api/sn/sn_rt_stations.m.
 *
 * `sn.rt` is over STATEFUL nodes, and a Router or a Cache is stateful without
 * being a station. Callers that need a station-to-station routing matrix must
 * therefore eliminate the non-station stateful rows, not index around them:
 * the reference forms the stochastic complement
 *
 *   rtst = P(A,A) + P(A,B) (I - P(B,B))^-1 P(B,A),
 *
 * with A the (station, class) rows in STATION order and B the rest. That is
 * exactly `dtmc_stochcomp(rt, A)`, which is what this port calls; the explicit
 * block form is only how the reference spells it.
 *
 * `Vst` is `cellsum(sn.visits)` restricted to the same station rows, in the
 * same order, so the two outputs are index-compatible.
 *
 * ARITHMETIC: field. dtmc_stochcomp is a linear solve.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/dtmc_stochcomp.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** What sn_rt_stations returns: the complemented routing and the station visits. */
template <class T>
struct SnRtStations {
    Matrix<T> rtst;  ///< (nstations*nclasses) square, STATION-major rows
    Matrix<T> Vst;   ///< (nstations x nclasses), cellsum(sn.visits) by station
};

template <class T>
SnRtStations<T> sn_rt_stations(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = sn.nclasses, M = sn.nstations, S = sn.nof_stateful();
    SnRtStations<T> out;
    std::vector<std::size_t> keep;
    keep.reserve(M * K);
    for (std::size_t ist = 0; ist < M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist + 1) - 1;
        for (std::size_t r = 0; r < K; ++r) keep.push_back(isf * K + r);
    }
    if (sn.rt.rows() == S * K) {
        out.rtst = mc::dtmc_stochcomp(sn.rt, keep);
    } else {
        out.rtst = Matrix<T>(M * K, M * K, zero);
    }
    out.Vst = Matrix<T>(M, K, zero);
    for (std::size_t ist = 0; ist < M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist + 1) - 1;
        for (std::size_t c = 0; c < sn.visits.size(); ++c)
            for (std::size_t r = 0; r < K; ++r)
                out.Vst(ist, r) = T(out.Vst(ist, r) + sn.visits[c](isf, r));
    }
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_RT_STATIONS_H
