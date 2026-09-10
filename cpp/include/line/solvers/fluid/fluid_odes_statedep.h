/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_ODES_STATEDEP_H
#define LINE_SOLVERS_FLUID_FLUID_ODES_STATEDEP_H

/**
 * The state-dependent fluid drifts: ports of `ode_statedep.m`, `ode_softmin.m`
 * and `ode_pnorm.m`.
 *
 * HOW THESE DIFFER FROM `closing`. The closing drift factors every rate into a
 * constant times one state entry, which is why it can precompute `rateBase`
 * and evaluate the drift as a sum over events. These three cannot: the server
 * share a job receives depends on the whole station's occupancy, and for FCFS
 * on the MEAN SERVICE TIME OF THE PHASE the job is in, so the multiplier
 * changes per (class, phase) at every step. They are therefore written the way
 * the reference writes them -- accumulate directly into dx, station by station
 * -- and are correspondingly slower. That is the trade the reference names in
 * its own comment: "slower than ODE_RATES_STATEINDEP, but allows rates that
 * are more complex functions of x".
 *
 * WHAT THE THREE VARY. Only how a saturated station's capacity enters:
 *
 *   statedep  the hard `min(ni, c)`, which has a kink at ni = c
 *   softmin   `softmin(ni, c, alpha)`, the weighted average
 *             (x e^-ax + y e^-ay)/(e^-ax + e^-ay), smooth everywhere
 *   pnorm     `ghat = 1/(1 + (ni/c)^p)^(1/p)`, a smooth stand-in for
 *             min(1, c/ni) from Ruuskanen et al., PEVA 151 (2021)
 *
 * The kink is what makes the closing/statedep drift stiff near saturation, and
 * smoothing it is what lets the integrator take longer steps. The fixed point
 * moves slightly in exchange, which is why these are separate methods rather
 * than a faster way to compute the same answer.
 *
 * OPEN MODELS ARE REFUSED. `ode_statedep.m` errors on an EXT station -- the
 * family has no source term -- and so does this port, by name.
 *
 * A NOTE ON TWO REFERENCE ASYMMETRIES, both reproduced deliberately.
 *   1. At an INF station the completion loop SKIPS j == i, so a delay that
 *      routes to itself contributes no flow; the PS, FCFS and DPS branches do
 *      not skip it. This is `ode_statedep.m`'s `if j~=i`.
 *   2. `ode_pnorm.m`'s DPS branch scales by `ghat * w/wni` where statedep
 *      scales by `nservers * w/wni` -- ghat where the others have the server
 *      count. Faithful to the reference; noted because it does not follow from
 *      the smoothing argument.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/util/error.h"

namespace line {
namespace fluid {

/** Which smoothing the drift applies at a saturated station. */
enum class StateDepKind { StateDep, SoftMin, PNorm };

/** Everything the state-dependent drifts read, lowered once to doubles. */
struct FluidStateDepSystem {
    FluidLayout layout;
    std::vector<std::vector<std::vector<double>>> mu, phi, pie;
    std::vector<std::vector<Matrix<double>>> d0;   ///< the PH generator per (i,c)
    std::vector<std::vector<std::vector<std::vector<double>>>> rt;  ///< rt[i][c][j][l]
    std::vector<lang::SchedStrategy> sched;
    std::vector<double> nservers;
    std::vector<std::vector<double>> weight;
    StateDepKind kind = StateDepKind::StateDep;
    double alpha = 20.0;  ///< softmin sharpness
    double pstar = 20.0;  ///< p-norm exponent
};

namespace detail {

/** Port of `util/softmin.m`, including its overflow guard. */
inline double fluid_softmin(double x, double y, double alpha) {
    const double lo = std::min(x, y), hi = std::max(x, y);
    const double gap = hi - lo;
    // exp(-alpha*gap) underflows past 745/alpha; there the min IS the answer.
    if (!(gap < 745.0 / alpha)) return lo;
    const double w = std::exp(-alpha * gap);
    return lo + gap * w / (1.0 + w);
}

/** Port of `util/pnorm_smooth.m`: a smooth stand-in for min(1, c/x). */
inline double fluid_pnorm_smooth(double x, double c, double p) {
    if (x <= 0.0 || c <= 0.0) return 0.0;
    const double ratio = x / c;
    double g;
    if (p <= 0.0) {
        g = std::min(1.0, c / x);  // the reference's hard-min fallback
    } else {
        g = 1.0 / std::pow(1.0 + std::pow(ratio, p), 1.0 / p);
    }
    return std::isnan(g) ? 0.0 : g;
}

}  // namespace detail

/** Assemble what the state-dependent drifts need from `sn`. */
template <class T>
FluidStateDepSystem fluid_statedep_system(const qn::NetworkStruct<T>& sn, StateDepKind kind,
                                          double alpha = 20.0, double pstar = 20.0) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    FluidStateDepSystem s;
    s.kind = kind;
    s.alpha = alpha;
    s.pstar = pstar;
    s.layout = fluid_layout(sn);
    const FluidLayout& L = s.layout;

    s.mu.assign(M, std::vector<std::vector<double>>(K));
    s.phi.assign(M, std::vector<std::vector<double>>(K));
    s.pie.assign(M, std::vector<std::vector<double>>(K));
    s.d0.assign(M, std::vector<Matrix<double>>(K));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) {
                s.pie[i][r] = std::vector<double>{1.0};
                continue;
            }
            detail::fluid_mu_phi(sn.service[i][r], s.mu[i][r], s.phi[i][r]);
            s.pie[i][r] = detail::fluid_pie(sn.service[i][r]);
            const std::size_t n = sn.service[i][r].D0.rows();
            s.d0[i][r] = Matrix<double>(n, n, 0.0);
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b)
                    s.d0[i][r](a, b) = num_traits<T>::to_double(sn.service[i][r].D0(a, b));
        }

    const std::size_t S = sn.nof_stateful();
    const bool have_rt = sn.rt.rows() == S * K;
    std::vector<std::size_t> sf(M, 0);
    for (std::size_t i = 0; i < M; ++i) sf[i] = sn.stateful_of_station(i + 1) - 1;
    s.rt.assign(M, std::vector<std::vector<std::vector<double>>>(
                       K, std::vector<std::vector<double>>(M, std::vector<double>(K, 0.0))));
    for (std::size_t i = 0; i < M && have_rt; ++i)
        for (std::size_t c = 0; c < K; ++c)
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t l = 0; l < K; ++l)
                    s.rt[i][c][j][l] = num_traits<T>::to_double(sn.rt(sf[i] * K + c, sf[j] * K + l));

    double closed_pop = 0.0;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population)) closed_pop += sn.classes[r].population;
    s.sched.resize(M);
    s.nservers.resize(M);
    s.weight.assign(M, std::vector<double>(K, 1.0));
    for (std::size_t i = 0; i < M; ++i) {
        s.sched[i] = sn.stations[i].sched;
        const double c = sn.stations[i].nservers;
        s.nservers[i] = std::isfinite(c) ? c : closed_pop;
        if (sn.stations[i].sched == lang::SchedStrategy::DPS)
            for (std::size_t r = 0; r < K && r < sn.stations[i].schedparam.size(); ++r)
                s.weight[i][r] = num_traits<T>::to_double(sn.stations[i].schedparam[r]);
        if (sn.stations[i].sched == lang::SchedStrategy::EXT)
            throw UnsupportedError(
                "fluid: the 'statedep', 'softmin' and 'pnorm' methods have no source term and are "
                "refused on an open model (station '" + sn.stations[i].name +
                "' is a Source); use method 'closing' or 'matrix'");
    }
    return s;
}

/**
 * The drift dx/dt for the state-dependent family.
 *
 * Written as the reference writes it: for each station, first the phase
 * changes, then the completions, each moving `x[from] * rate` of mass.
 */
inline std::function<void(double, const double*, double*)> fluid_drift_statedep(
    const FluidStateDepSystem& s) {
    const std::size_t n = s.layout.nstates;
    return [s, n](double, const double* x, double* dx) {
        const FluidLayout& L = s.layout;
        const std::size_t M = L.qidx.size();
        const std::size_t K = M ? L.qidx[0].size() : 0;
        for (std::size_t i = 0; i < n; ++i) dx[i] = 0.0;

        for (std::size_t i = 0; i < M; ++i) {
            const lang::SchedStrategy sc = s.sched[i];
            const double c = s.nservers[i];

            // The station's total occupancy, which every branch but INF needs.
            double ni = 0.0;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t k = 0; k < L.kic[i][r]; ++k) ni += x[L.qidx[i][r] + k];

            // FCFS weights the share by the mean duration of the phase a job
            // is in: w = -1/D0(k,k). DPS weights it by the class weight.
            std::vector<std::vector<double>> wfcfs;
            std::vector<double> wdps;
            double wni = lang::GlobalConstants::FineTol;
            if (sc == lang::SchedStrategy::FCFS) {
                wfcfs.assign(K, std::vector<double>());
                for (std::size_t r = 0; r < K; ++r) {
                    wfcfs[r].assign(L.kic[i][r], 0.0);
                    if (!L.enabled[i][r]) continue;
                    for (std::size_t k = 0; k < L.kic[i][r]; ++k) {
                        const double d = s.d0[i][r](k, k);
                        wfcfs[r][k] = (d != 0.0) ? -1.0 / d : 0.0;
                        wni += wfcfs[r][k] * x[L.qidx[i][r] + k];
                    }
                }
            } else if (sc == lang::SchedStrategy::DPS) {
                double wsum = 0.0;
                for (std::size_t r = 0; r < K; ++r) wsum += s.weight[i][r];
                wdps.assign(K, 0.0);
                for (std::size_t r = 0; r < K; ++r)
                    wdps[r] = (wsum > 0.0) ? s.weight[i][r] / wsum : 0.0;
                for (std::size_t r = 0; r < K; ++r) {
                    double blk = 0.0;
                    for (std::size_t k = 0; k < L.kic[i][r]; ++k) blk += x[L.qidx[i][r] + k];
                    wni += wdps[r] * blk;
                }
            }

            // The smoothed capacity, shared by every rate at this station.
            double ghat = 1.0;
            if (s.kind == StateDepKind::PNorm)
                ghat = (ni > 0.0 && c > 0.0) ? detail::fluid_pnorm_smooth(ni, c, s.pstar) : 1.0;
            const double capped = (s.kind == StateDepKind::SoftMin)
                                      ? detail::fluid_softmin(ni, c, s.alpha)
                                      : std::min(ni, c);

            // The multiplier applied to a rate of class r in phase k.
            const auto factor = [&](std::size_t r, std::size_t k) -> double {
                switch (sc) {
                    case lang::SchedStrategy::INF:
                        return 1.0;
                    case lang::SchedStrategy::FCFS:
                        // statedep/softmin: min-or-softmin times the phase share.
                        // pnorm: ghat in place of the capacity.
                        if (s.kind == StateDepKind::PNorm) return ghat * wfcfs[r][k] / wni;
                        return capped * wfcfs[r][k] / wni;
                    case lang::SchedStrategy::DPS:
                        if (s.kind == StateDepKind::PNorm) return ghat * wdps[r] / wni;
                        return (ni > c) ? c * wdps[r] / wni : 1.0;
                    default:  // PS and anything else the reference treats as PS
                        if (s.kind == StateDepKind::PNorm) return ghat;
                        return (ni > c && ni > 0.0) ? c / ni : 1.0;
                }
            };

            // ---- phase changes within (i,r) --------------------------------
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.enabled[i][r]) continue;
                const std::size_t b = L.qidx[i][r];
                for (std::size_t k = 0; k + 1 < L.kic[i][r]; ++k)
                    for (std::size_t kp = 0; kp < L.kic[i][r]; ++kp) {
                        if (kp == k) continue;
                        const double flow = x[b + k] * s.d0[i][r](k, kp) * factor(r, k);
                        dx[b + k] -= flow;
                        dx[b + kp] += flow;
                    }
            }

            // ---- service completions ---------------------------------------
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.enabled[i][r]) continue;
                const std::size_t b = L.qidx[i][r];
                for (std::size_t j = 0; j < M; ++j) {
                    // The reference's `if j~=i`, INF only: see the header.
                    if (sc == lang::SchedStrategy::INF && j == i) continue;
                    for (std::size_t l = 0; l < K; ++l) {
                        if (!L.enabled[j][l]) continue;
                        const double p = s.rt[i][r][j][l];
                        if (!(p > 0.0)) continue;
                        const std::size_t bj = L.qidx[j][l];
                        for (std::size_t k = 0; k < L.kic[i][r]; ++k) {
                            const double base = s.phi[i][r][k] * s.mu[i][r][k] * p * factor(r, k);
                            for (std::size_t kj = 0; kj < L.kic[j][l]; ++kj) {
                                const double flow = x[b + k] * base * s.pie[j][l][kj];
                                dx[b + k] -= flow;
                                dx[bj + kj] += flow;
                            }
                        }
                    }
                }
            }
        }
    };
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_ODES_STATEDEP_H
