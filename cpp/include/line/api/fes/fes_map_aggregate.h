/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_AGGREGATE_H
#define LINE_API_FES_MAP_AGGREGATE_H

/**
 * Recursive MAP flow-equivalent server for a station subset.
 *
 * Templated port of matlab/src/api/fes/fes_map_aggregate.m, mirrored by the JAR
 * and native Python. Implements the recursion of Section 5.2.1 of Casale, Mi,
 * Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011.
 *
 * The first station seeds the flow-equivalent server; every further station is
 * folded against the running server by building the inter-departure MAP of the
 * resulting pair at each population level and fitting a MAP(2) to its first three
 * moments and index of dispersion. The result is one MAP per level, the service
 * process of a single load-dependent station replacing the whole subnetwork.
 * Unlike the classic flow-equivalent server, which keeps only the mean throughput
 * of the subnetwork, this one also carries the burstiness of its departure
 * stream, so a bottleneck switch across the aggregated resources stays visible to
 * the rest of the model.
 *
 * Levels are evaluated on a grid and the four descriptors are interpolated
 * between grid points, as MAPs fitted at neighbouring populations are similar.
 * The MAP is refitted at every level from the interpolated descriptors, never
 * interpolated entrywise.
 *
 * ARITHMETIC: transcendental, inherited from map2_fit_idc.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/fes/fes_map_interdeparture.h"
#include "line/api/fes/fes_map_interp.h"
#include "line/api/fes/fes_map_levels.h"
#include "line/api/fes/fes_map_moments.h"
#include "line/api/mam/map2_fit_idc.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fes {

/** The load-dependent MAP that replaces a subnetwork, with its descriptors. */
template <class T>
struct FesMapAggregateResult {
    /** Service process of the flow-equivalent server, index k-1 holding k jobs. */
    std::vector<mam::Map<T>> fes;
    /** Throughput of the subnetwork at each population. */
    std::vector<T> throughput;
    /** Descriptors e1, e2, e3 and the index of dispersion at each population. */
    std::vector<std::vector<T>> moments;
    /** Fit status at each population, see map2_fit_idc. */
    std::vector<int> status;
    /** Populations at which the inter-departure MAP was evaluated. */
    std::vector<std::size_t> grid;
};

/**
 * @param maps    service process of each station, already scaled by its visit ratio
 * @param servers number of servers of each station, infinite for a delay
 * @param n       largest population the flow-equivalent server must serve
 * @param grid    populations at which the inter-departure MAP is evaluated
 * @param method  moment evaluation method, "ssolve" or "euler"
 */
template <class T>
FesMapAggregateResult<T> fes_map_aggregate(const std::vector<mam::Map<T>>& maps,
                                           const std::vector<double>& servers, std::size_t n,
                                           const std::vector<std::size_t>& grid,
                                           const std::string& method = "ssolve") {
    static_assert(num_traits<T>::has_transcendental,
                  "fes_map_aggregate inherits map2_fit_idc's arithmetic");
    const std::size_t M = maps.size();
    if (M < 1) throw InputError("fes_map_aggregate: at least one station is required");
    if (servers.size() != M)
        throw InputError("fes_map_aggregate: one server count per station is required");

    FesMapAggregateResult<T> out;
    out.grid = grid;
    out.status.assign(n, 0);
    out.moments.assign(4, std::vector<T>(n, num_traits<T>::from_int(0)));

    std::vector<mam::Map<T>> fes = fes_map_levels(maps[0], n, servers[0]);
    for (std::size_t k = 0; k < n; ++k) {
        out.moments[0][k] = mam::map_moment(fes[k], 1);
        out.moments[1][k] = mam::map_moment(fes[k], 2);
        out.moments[2][k] = mam::map_moment(fes[k], 3);
        out.moments[3][k] = mam::map_idc(fes[k]);
    }

    for (std::size_t i = 1; i < M; ++i) {
        const std::vector<mam::Map<T>> stationLev = fes_map_levels(maps[i], n, servers[i]);
        std::vector<std::vector<T>> gmom(4, std::vector<T>(grid.size()));
        for (std::size_t g = 0; g < grid.size(); ++g) {
            const mam::Map<T> T01 = fes_map_interdeparture(stationLev, fes, grid[g]);
            const FesMapMoments<T> mom = fes_map_moments(T01, method);
            gmom[0][g] = mom.e1;
            gmom[1][g] = mom.e2;
            gmom[2][g] = mom.e3;
            gmom[3][g] = mom.idc;
        }

        if (grid.size() < n) {
            std::vector<T> xs(grid.size()), xq(n);
            for (std::size_t g = 0; g < grid.size(); ++g)
                xs[g] = num_traits<T>::from_int(static_cast<long>(grid[g]));
            for (std::size_t k = 0; k < n; ++k)
                xq[k] = num_traits<T>::from_int(static_cast<long>(k + 1));
            for (std::size_t r = 0; r < 4; ++r) out.moments[r] = fes_map_interp(xs, gmom[r], xq);
        } else {
            out.moments = gmom;
        }

        std::vector<mam::Map<T>> newFes;
        newFes.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            const mam::Map2FitIdcResult<T> fit = mam::map2_fit_idc(
                out.moments[0][k], out.moments[1][k], out.moments[2][k], out.moments[3][k]);
            newFes.push_back(fit.map);
            out.status[k] = fit.status;
        }
        fes = newFes;
    }

    out.fes = fes;
    out.throughput.assign(n, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < n; ++k)
        out.throughput[k] = num_traits<T>::from_int(1) / out.moments[0][k];
    return out;
}

/** Aggregation on the default population grid. */
template <class T>
FesMapAggregateResult<T> fes_map_aggregate(const std::vector<mam::Map<T>>& maps,
                                           const std::vector<double>& servers, std::size_t n) {
    return fes_map_aggregate(maps, servers, n, fes_map_grid(n), "ssolve");
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_AGGREGATE_H
