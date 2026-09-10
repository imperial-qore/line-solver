/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_BUSYP_H
#define LINE_SOLVERS_NC_SOLVER_NC_BUSYP_H

/**
 * Port of `@SolverNC/getAvgBusyPeriod.m` and of the Python-native
 * `SolverNC.getAvgBusyPeriod`: the mean busy period of order n for a set of
 * stations, from `pfqn_busyp` (Daduna, J. ACM 35(3), 1988).
 *
 * WHAT THIS LAYER ADDS OVER THE API. `pfqn_busyp` takes the paper's inputs --
 * relative arrival rates, a rate law, a routing matrix -- and this reads them
 * off a NetworkStruct:
 *
 *   alpha = Vchain(:,1), the chain visit ratios. The CLOSED formula is
 *     homogeneous of degree zero in alpha, so unnormalized visits serve; the
 *     OPEN one is not, and alpha there is scaled to absolute rates by lambda.
 *   mu(j,k) = scaling(j,k) / STchain(j), a RATE and not the dimensionless
 *     lldscaling `pfqn_ncld` takes. The precedence is the one `solver_ncld`
 *     uses: an infinite server first (scaling = k), then a declared
 *     lldscaling row, then the multiserver staircase min(k,c).
 *   P = the station-to-station chain routing, the class blocks of
 *     `sn_rt_stations` weighted by the class visits. That weighting is exact,
 *     being a flow balance and not an approximation.
 *
 * AN OPEN MODEL DROPS THE SOURCE. The Jackson network of the paper has no
 * Source station: its outflow is the external stream gamma, so the Source is
 * removed from the node set and `gamma_j = lambda * P(source, j)`.
 *
 * SINGLE CHAIN ONLY. The paper is written for identical customers; Section 5
 * only sketches the multichain extension, which no codebase implements, so a
 * multichain model is refused rather than answered from a chain aggregate that
 * the theorem does not cover.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_busyp.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"

namespace line {
namespace nc {

/**
 * Mean busy period of order n for a set of stations.
 *
 * @param sn     the network structure, single-chain
 * @param subnet zero-based STATION indexes forming the subnetwork
 * @param orders busy period orders, 1 for the ordinary busy period
 * @return one duration per requested order
 */
template <class T>
std::vector<double> solver_nc_busyp(const qn::NetworkStruct<T>& sn,
                                    const std::vector<std::size_t>& subnet,
                                    const std::vector<std::size_t>& orders) {
    if (sn.nchains > 1)
        throw UnsupportedError(
            "solver_nc_busyp: the busy period of a subnetwork is defined for "
            "single-chain models only; Section 5 of Daduna (1988) sketches the "
            "multichain extension, which is not implemented");

    const std::size_t M = sn.nstations, K = sn.nclasses;
    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
    const api::SnRtStations<T> rt = api::sn_rt_stations(sn);

    // station-to-station routing of the chain, the class blocks weighted by the
    // class visits: exact, being a flow balance
    Matrix<double> Pst(M, M, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        double vtot = 0.0;
        for (std::size_t r = 0; r < K; ++r) vtot += num_traits<T>::to_double(rt.Vst(i, r));
        for (std::size_t j = 0; j < M; ++j) {
            double flow = 0.0;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s)
                    flow += num_traits<T>::to_double(rt.Vst(i, r)) *
                            num_traits<T>::to_double(rt.rtst(i * K + r, j * K + s));
            Pst(i, j) = (vtot > 0) ? flow / vtot : 0.0;
        }
    }

    std::vector<double> st_time(M, 0.0), visits(M, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        st_time[i] = num_traits<T>::to_double(d.STchain(i, 0));
        visits[i] = num_traits<T>::to_double(d.Vchain(i, 0));
    }

    // mu(j,k): the rate of station j holding k jobs, in the solver_ncld order
    const auto rate_of = [&](std::size_t j, std::size_t k) -> double {
        double scaling;
        const double servers = num_traits<T>::to_double(sn.stations[j].nservers);
        const std::vector<T>& lld = sn.stations[j].lldscaling;
        if (!std::isfinite(servers)) {
            scaling = static_cast<double>(k);
        } else if (!lld.empty()) {
            const std::size_t idx = std::min(k, lld.size()) - 1;
            scaling = num_traits<T>::to_double(lld[idx]);
        } else {
            scaling = std::min(static_cast<double>(k), servers);
        }
        return scaling / st_time[j];
    };

    const double N = d.Nchain.empty() ? std::numeric_limits<double>::infinity()
                                      : d.Nchain[0];
    if (!std::isfinite(N)) {
        // the Source is not a node of the Jackson network: its outflow is gamma
        std::size_t source = M;
        for (std::size_t i = 0; i < M; ++i)
            if (sn.stations[i].nodetype == lang::NodeType::Source) {
                source = i;
                break;
            }
        if (source == M)
            throw InputError("solver_nc_busyp: an open model must own a Source station");
        for (std::size_t t = 0; t < subnet.size(); ++t)
            if (subnet[t] == source)
                throw InputError(
                    "solver_nc_busyp: the Source cannot belong to the subnetwork");

        double lambda = 0.0;
        for (std::size_t r = 0; r < K; ++r) {
            const double rate = num_traits<T>::to_double(sn.rates(source, r));
            if (std::isfinite(rate)) lambda += rate;
        }
        std::vector<std::size_t> keep;
        std::vector<std::size_t> remap(M, 0);
        for (std::size_t i = 0; i < M; ++i)
            if (i != source) {
                remap[i] = keep.size();
                keep.push_back(i);
            }
        std::vector<double> alpha(keep.size(), 0.0), gamma(keep.size(), 0.0);
        Matrix<double> P(keep.size(), keep.size(), 0.0);
        for (std::size_t i = 0; i < keep.size(); ++i) {
            alpha[i] = lambda * visits[keep[i]] / visits[source];
            gamma[i] = lambda * Pst(source, keep[i]);
            for (std::size_t j = 0; j < keep.size(); ++j) P(i, j) = Pst(keep[i], keep[j]);
        }
        std::vector<std::size_t> mapped;
        for (std::size_t t = 0; t < subnet.size(); ++t) mapped.push_back(remap[subnet[t]]);
        const std::function<double(std::size_t, std::size_t)> mu =
            [&](std::size_t j, std::size_t k) { return rate_of(keep[j], k); };
        return pfqn::pfqn_busyp(alpha, mu, P, N, mapped, orders, gamma).b;
    }

    const std::function<double(std::size_t, std::size_t)> mu = rate_of;
    return pfqn::pfqn_busyp(visits, mu, Pst, N, subnet, orders).b;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_BUSYP_H
