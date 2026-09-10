/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_DIFFUSION_H
#define LINE_SOLVERS_FLUID_FLUID_DIFFUSION_H

/**
 * The `diffusion` method: a port of `solver_fluid_diffusion.m`.
 *
 * WHAT IT COMPUTES. The fluid limit follows the MEAN drift and says nothing
 * about fluctuation. The diffusion approximation adds a noise term and
 * integrates the resulting stochastic differential equation by
 * Euler-Maruyama,
 *
 *     x <- max(0, x + drift(x) dt + sqrt(dt) Z),   Z standard normal,
 *
 * renormalising each class back to its population after every step so the
 * closed network stays closed. The reported queue lengths are the TIME AVERAGE
 * over the whole trajectory, not the end state: a single noisy path is
 * meaningless at its last instant and informative in the mean.
 *
 * THE DRIFT HERE IS SIMPLER THAN THE FLUID ONE. It is flow in minus flow out
 * with the service rate applied to the whole queue -- x/mu_inv -- and no
 * server-sharing term at all. That is why the reference restricts the method
 * to single-server and infinite-server stations: with c = 1 the fluid rate
 * min(x, c)/mu_inv and this x/mu_inv differ, and the reference chooses the
 * latter, so the two methods answer slightly different questions. All the
 * restrictions below are the reference's and are enforced by name.
 *
 * PARITY IS STATISTICAL, NOT EXACT. The trajectory is driven by pseudorandom
 * normals, and MATLAB's `randn` and this port's Mersenne Twister produce
 * different streams from the same seed. Two runs therefore agree in
 * distribution and not digit for digit -- the same situation the SSA solvers
 * are in across codebases. Tests must assert on averages and invariants
 * (population conservation, ordering), never on a specific trajectory.
 */

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "line/api/mc/dtmc_stochcomp.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** Controls of the Euler-Maruyama trajectory. */
struct DiffusionOptions {
    std::size_t steps = 10000;      ///< number of steps (the reference's iter_max)
    double dt = 0.01;               ///< step size (the reference's timestep)
    unsigned long seed = 23000;     ///< RNG seed; see the header on parity
};

/** Time-averaged queue lengths of the diffusion trajectory. */
struct DiffusionResult {
    Matrix<double> QN;
};

/**
 * Run the diffusion approximation of `sn`.
 *
 * Refuses by name every model shape the reference refuses: open classes, a
 * Source, a discipline outside {PS, FCFS, INF, SIRO}, and any finite
 * multiserver station.
 */
template <class T>
DiffusionResult fluid_diffusion(const qn::NetworkStruct<T>& sn, const DiffusionOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    for (std::size_t r = 0; r < K; ++r)
        if (!std::isfinite(sn.classes[r].population))
            throw UnsupportedError(
                "fluid diffusion: the method supports closed networks only; class '" +
                sn.classes[r].name + "' is open");
    for (std::size_t i = 0; i < M; ++i) {
        const lang::SchedStrategy sc = sn.stations[i].sched;
        if (sc == lang::SchedStrategy::EXT)
            throw UnsupportedError("fluid diffusion: a Source is not supported (station '" +
                                   sn.stations[i].name + "'); the method is for closed networks");
        if (!(sc == lang::SchedStrategy::PS || sc == lang::SchedStrategy::FCFS ||
              sc == lang::SchedStrategy::INF || sc == lang::SchedStrategy::SIRO))
            throw UnsupportedError(
                "fluid diffusion: scheduling at station '" + sn.stations[i].name +
                "' is outside the supported set {PS, FCFS, INF, SIRO}");
        const double c = sn.stations[i].nservers;
        if (std::isfinite(c) && c > 1.0)
            throw UnsupportedError(
                "fluid diffusion: only single-server or infinite-server stations are supported; "
                "station '" + sn.stations[i].name + "' has more than one server");
    }

    // Mean service time per (station, class); Inf marks "not served here".
    Matrix<double> mu_inv(M, K, std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (sn.disabled[i][r] || sn.service[i][r].D0.rows() == 0) continue;
            const double rate = num_traits<T>::to_double(sn.rates(i, r));
            if (rate > 0.0 && std::isfinite(rate)) mu_inv(i, r) = 1.0 / rate;
        }

    // Station-space routing, as the reference builds it.
    const std::size_t S = sn.nof_stateful();
    std::vector<std::size_t> keep;
    keep.reserve(M * K);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
        for (std::size_t r = 0; r < K; ++r) keep.push_back(isf * K + r);
    }
    Matrix<double> rt_full(S * K, S * K, 0.0);
    if (sn.rt.rows() == S * K)
        for (std::size_t a = 0; a < S * K; ++a)
            for (std::size_t b = 0; b < S * K; ++b)
                rt_full(a, b) = num_traits<T>::to_double(sn.rt(a, b));
    const Matrix<double> P = mc::dtmc_stochcomp(rt_full, keep);

    const std::size_t steps = std::max<std::size_t>(2, opt.steps);
    std::mt19937_64 rng(opt.seed);
    std::normal_distribution<double> gauss(0.0, 1.0);

    // Start with each class spread evenly over the stations, as the reference.
    Matrix<double> x(M, K, 0.0), avg(M, K, 0.0);
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t i = 0; i < M; ++i)
            x(i, r) = sn.classes[r].population / static_cast<double>(M);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) avg(i, r) = x(i, r) / static_cast<double>(steps);

    Matrix<double> xn(M, K, 0.0);
    for (std::size_t step = 1; step < steps; ++step) {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                // drift = inflow - outflow, both at the full queue rate
                const double out =
                    std::isfinite(mu_inv(i, r)) && mu_inv(i, r) > 0.0 ? x(i, r) / mu_inv(i, r) : 0.0;
                double in = 0.0;
                for (std::size_t j = 0; j < M; ++j)
                    for (std::size_t q = 0; q < K; ++q) {
                        if (!(std::isfinite(mu_inv(j, q)) && mu_inv(j, q) > 0.0)) continue;
                        in += (x(j, q) / mu_inv(j, q)) * P(j * K + q, i * K + r);
                    }
                const double dW = std::sqrt(opt.dt) * gauss(rng);
                double v = x(i, r) + (in - out) * opt.dt + dW;
                if (v < 0.0) v = 0.0;
                xn(i, r) = v;
            }
        // Renormalise each class back to its population: the network is closed.
        for (std::size_t r = 0; r < K; ++r) {
            double tot = 0.0;
            for (std::size_t i = 0; i < M; ++i) tot += xn(i, r);
            for (std::size_t i = 0; i < M; ++i)
                xn(i, r) = (tot > 0.0) ? xn(i, r) * sn.classes[r].population / tot
                                       : sn.classes[r].population / static_cast<double>(M);
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                x(i, r) = xn(i, r);
                avg(i, r) += x(i, r) / static_cast<double>(steps);
            }
    }

    DiffusionResult out;
    out.QN = avg;
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_DIFFUSION_H
