/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_CDF_H
#define LINE_SOLVERS_NC_SOLVER_NC_CDF_H

/**
 * Port of `@@SolverNC/getCdfRespT.m`, and of its two aliases `getSjrnT` and
 * `sjrnT`.
 *
 * WHAT IT ADDS OVER THE AVERAGE TABLE. `getAvg` reports E[R]; this reports the
 * whole law Pr[R <= t], which is what a percentile or a service-level target
 * needs. The two are not interchangeable: a product-form response time is not
 * exponential, so E[R] does not determine its tail.
 *
 * TWO KINDS OF STATION, TWO DIFFERENT LAWS. At an FCFS queue the sojourn time
 * is the tagged-job passage time through a closed network and comes from
 * `pfqn_stdf` (exact) or `pfqn_stdf_heur` (the `rd` heuristic). At a delay
 * station there is no queueing at all, so the sojourn time IS the service time
 * and the law is the distribution's own CDF, taken through `map_cdf`. The
 * reference computes both and returns them in one station-indexed table.
 *
 * WHY THE TIME GRID IS WHAT IT IS. The reference builds
 * `logspace(0, 2 log10 T, 100)` with `T = max(sum(N) * mean(1/rates))` over the
 * FCFS stations: the population times the mean service time is the scale on
 * which the FCFS subsystem empties, and squaring it gives a grid that still
 * resolves the tail. It is reproduced exactly, because the grid is part of the
 * answer -- the CDF is reported AT these points. Note which axis that `mean`
 * runs along: `sn.rates(fcfsNodeIds,:)` is (stations x classes) and MATLAB's
 * `mean` takes the FIRST NON-SINGLETON dimension, so it averages over STATIONS
 * and `max` then runs over classes -- except with a single FCFS station, where
 * the row is 1xR and the same call averages over CLASSES instead. Both branches
 * are reproduced below; see `_kb/06-solver-catalog.md` (NC, getCdfRespT).
 *
 * TWO PROPERTIES OF THAT GRID TO KNOW BEFORE READING ITS OUTPUT, neither of
 * them a porting artifact. It DESCENDS when T < 1, because 2 log10 T is then
 * negative, and it COLLAPSES to a hundred identical points at T = 1 exactly.
 * Delay(1)+FCFS(3) with two classes at N=3 hits the second case: T = 3 * 1/3.
 * The law is still correct at whatever points the grid names.
 *
 * REFUSED, as the reference refuses it: a model with no FCFS station. MATLAB
 * warns and returns an empty cell array, which reads as a completed analysis
 * with no data; this throws instead.
 *
 * ARITHMETIC. `pfqn_stdf` and `map_cdf` are transcendental, and the grid itself
 * is a `logspace`, so a non-transcendental backend is refused by name.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/pfqn/pfqn_stdf.h"
#include "line/api/pfqn/pfqn_stdf_heur.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/sn_pf_params.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/**
 * The response-time distributions, station by class.
 *
 * `RD[i][r]` is a two-column matrix whose first column is F(t) and whose second
 * is t, which is the shape `setDistribResults` stores in the reference. An
 * empty entry means the station-class pair has no law (a Queue that is not
 * FCFS, or a class the station does not serve).
 */
template <class T>
struct CdfRespTResult {
    std::vector<std::vector<Matrix<T>>> RD;
    std::vector<T> tset;  ///< the shared evaluation grid
    /**
     * Non-empty when the reference WARNS AND RETURNS EMPTY rather than
     * computing: today only "applies only to FCFS nodes". `RD` is then empty,
     * which a caller can test by size -- unlike the zero table the
     * normalizing-constant path returns when a method declines a model.
     */
    std::string warning;
};

/**
 * Port of `@@SolverNC/getCdfRespT.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls; `opt.cdf_algorithm` selects 'exact' or 'rd'
 */
template <class T>
CdfRespTResult<T> solver_nc_cdf_respt(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    CdfRespTResult<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc_cdf_respt: the sojourn-time law is evaluated on a logarithmic time grid "
            "and inverts a generating function; it needs transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, R = sn.nclasses;

        if (opt.cdf_algorithm != "exact" && opt.cdf_algorithm != "rd")
            throw UnsupportedError("solver_nc_cdf_respt: config.algorithm '" + opt.cdf_algorithm +
                                   "' is unsupported; use 'exact' (pfqn_stdf) or 'rd' "
                                   "(pfqn_stdf_heur)");

        // The reference indexes fcfsNodes into the rows of the NON-DELAY
        // stations and fcfsNodeIds into all stations; both are kept, because
        // pfqn_stdf is given the queueing rows only while the result is
        // reported against the full station list.
        std::vector<std::size_t> fcfsStations, delayStations;
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.stations[i].sched == qn::SchedStrategy::INF) {
                delayStations.push_back(i + 1);
            } else if (sn.stations[i].sched == qn::SchedStrategy::FCFS) {
                fcfsStations.push_back(i + 1);
            }
        }
        if (fcfsStations.empty()) {
            // ALIGNED TO MATLAB (the empty-result ruling of 2026-07-25, register
            // row N1 in a second shape). `getCdfRespT.m:44-45` warns
            // "getCdfRespT applies only to FCFS nodes" and RETURNS with RD = {},
            // never calling setDistribResults. This returns the same empty
            // result rather than throwing. It is less lossy than N1's zero
            // table: a caller can detect an empty result by SIZE, where a table
            // of zeros is indistinguishable from a model that holds no jobs.
            out.warning = "getCdfRespT applies only to FCFS nodes.";
            return out;
        }
        for (double n : sn.njobs())
            if (!std::isfinite(n))
                throw UnsupportedError(
                    "solver_nc_cdf_respt: the tagged-job sojourn-time law is defined on a CLOSED "
                    "network; this model has an open class");

        const PfParams<T> pf = sn_get_product_form_params(sn);

        // The FCFS stations in the index space of the NON-DELAY stations, which
        // is `fcfsNodes` in the reference and is what pfqn_stdf wants: the rows
        // of D are the Queue nodes, and on a closed model the non-delay
        // stations are exactly those.
        std::vector<std::size_t> fcfsNonDelayPos;
        {
            std::size_t pos = 0;
            for (std::size_t i = 0; i < M; ++i) {
                if (sn.stations[i].sched == qn::SchedStrategy::INF) continue;
                ++pos;
                if (sn.stations[i].sched == qn::SchedStrategy::FCFS) fcfsNonDelayPos.push_back(pos);
            }
        }

        // T = max over the FCFS stations of sum(N) * mean(1/rate).
        //
        // The reference USED to index `sn.rates` with `fcfsNodes`, which is in
        // the non-delay station space while sn.rates is indexed by station, so
        // with a Delay declared first the grid was scaled by the DELAY's service
        // time. Found by this port and fixed in `@@SolverNC/getCdfRespT.m` with
        // the user's approval; `fcfsNodeIds` is the station-indexed variable the
        // next line of the reference already computes. `fcfsNodes` is still
        // correct where it is passed to pfqn_stdf, and is left alone there.
        double Nsum = 0.0;
        for (double n : sn.njobs()) Nsum += n;
        const std::size_t K = fcfsStations.size();
        double Tmax = 0.0;
        bool anyRow = false;
        if (K == 1) {
            // 1xR: MATLAB's mean collapses the CLASS axis and max sees a scalar.
            const std::size_t ist = fcfsStations[0];
            double acc = 0.0;
            bool nan_row = false;
            for (std::size_t r = 0; r < R; ++r) {
                if (sn.disabled[ist - 1][r]) {
                    nan_row = true;  // MATLAB: 1/NaN is NaN and mean() propagates it
                    break;
                }
                acc += 1.0 / num_traits<T>::to_double(sn.rates(ist - 1, r));
            }
            if (!nan_row) {
                Tmax = Nsum * acc / static_cast<double>(R);
                anyRow = Tmax > 0.0;
            }
        } else {
            for (std::size_t r = 0; r < R; ++r) {
                double acc = 0.0;
                bool nan_col = false;
                for (std::size_t ist : fcfsStations) {
                    if (sn.disabled[ist - 1][r]) {
                        nan_col = true;  // one disabled pair makes the whole class NaN
                        break;
                    }
                    acc += 1.0 / num_traits<T>::to_double(sn.rates(ist - 1, r));
                }
                if (nan_col) continue;  // MATLAB's max ignores the NaN this class produces
                const double v = Nsum * acc / static_cast<double>(K);
                anyRow = true;
                if (v > Tmax) Tmax = v;
            }
        }
        if (!anyRow || !(Tmax > 0.0))
            throw NumericError(
                "solver_nc_cdf_respt: the time grid has no scale; every FCFS station's service "
                "rate is undefined or non-positive");

        // logspace(0, 2*log10(T), 100)
        const std::size_t npts = 100;
        const double hi = 2.0 * std::log10(Tmax);
        out.tset.assign(npts, zero);
        for (std::size_t j = 0; j < npts; ++j) {
            const double e = hi * static_cast<double>(j) / static_cast<double>(npts - 1);
            out.tset[j] = num_traits<T>::from_double(std::pow(10.0, e));
        }

        // pfqn_stdf wants the FCFS stations as 0-based indices into the ROWS OF
        // D, which are the Queue nodes. On a closed model the non-delay
        // stations ARE the queues, so the reference's `fcfsNodes` is already in
        // that space and only needs rebasing to zero.
        std::vector<std::size_t> fcfsRows;
        for (std::size_t pos : fcfsNonDelayPos) {
            if (pos > pf.queue_stations.size())
                throw UnsupportedError(
                    "solver_nc_cdf_respt: an FCFS station has no row in the product-form demand "
                    "matrix; the model has a non-Queue, non-Delay station between them");
            fcfsRows.push_back(pos - 1);
        }
        Matrix<T> rates(pf.queue_stations.size(), R, zero);
        for (std::size_t i = 0; i < pf.queue_stations.size(); ++i)
            for (std::size_t r = 0; r < R; ++r)
                if (!sn.disabled[pf.queue_stations[i] - 1][r])
                    rates(i, r) = sn.rates(pf.queue_stations[i] - 1, r);

        std::vector<int> N(R, 0);
        for (std::size_t r = 0; r < R; ++r)
            N[r] = static_cast<int>(std::llround(pf.N[r]));

        // pfqn_stdf takes integer server counts; an infinite one cannot appear
        // here because every row of D is a Queue node.
        std::vector<int> S(pf.S.size(), 1);
        for (std::size_t i = 0; i < pf.S.size(); ++i) {
            if (!std::isfinite(pf.S[i]))
                throw UnsupportedError(
                    "solver_nc_cdf_respt: a queueing station has infinitely many servers, which "
                    "the tagged-job passage time cannot represent");
            S[i] = static_cast<int>(std::llround(pf.S[i]));
        }
        const pfqn::StdfResult<T> sd =
            opt.cdf_algorithm == "exact"
                ? pfqn::pfqn_stdf(pf.D, N, pf.Z, S, fcfsRows, rates, out.tset)
                : pfqn::pfqn_stdf_heur(pf.D, N, pf.Z, S, fcfsRows, rates, out.tset);

        out.RD.assign(M, std::vector<Matrix<T>>(R));
        for (std::size_t i = 0; i < fcfsStations.size() && i < sd.RD.size(); ++i)
            for (std::size_t r = 0; r < R && r < sd.RD[i].size(); ++r)
                out.RD[fcfsStations[i] - 1][r] = sd.RD[i][r];

        // A delay station queues for nothing, so its sojourn law is the service
        // distribution itself.
        for (std::size_t ist : delayStations)
            for (std::size_t r = 0; r < R; ++r) {
                if (sn.disabled[ist - 1][r]) continue;
                const std::vector<T> F = mam::map_cdf(
                    lang::dist_to_map(sn.service[ist - 1][r]), out.tset);
                Matrix<T> A(npts, 2, zero);
                for (std::size_t j = 0; j < npts; ++j) {
                    A(j, 0) = F[j];
                    A(j, 1) = out.tset[j];
                }
                out.RD[ist - 1][r] = A;
            }
        return out;
    }
}

/** Port of `@@SolverNC/getSjrnT.m`: an alias of getCdfRespT. */
template <class T>
CdfRespTResult<T> solver_nc_sjrnt(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    return solver_nc_cdf_respt(sn, opt);
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_CDF_H
