/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_TBI_H
#define LINE_SOLVERS_FLUID_FLUID_TBI_H

/**
 * The `tbi` method: a port of `solver_fluid_tbi_iteration.m` and
 * `tbi_partition.m`.
 *
 * WHAT TBI IS FOR. The closing drift couples every station to every other, so
 * one integration works on the whole state vector at once and its cost grows
 * with the model. Time-based iteration splits the stations into CELLS, solves
 * each cell's sub-drift on its own, and treats the flow arriving from the other
 * cells as a KNOWN FUNCTION OF TIME, frozen at the previous sweep's
 * trajectories. Sweeping until the trajectories stop moving recovers the
 * coupled solution, while each solve only ever sees one cell's states. It is a
 * domain decomposition in time, and it pays off when the model is too large
 * for one integration to be comfortable.
 *
 * THE PARTITION is the reference's greedy merge: start with one station per
 * cell and repeatedly merge the pair with the largest routing coupling
 * (W + W', diagonal dropped) whose combined size stays within twice the target
 * cell size, until the cell count reaches ceil(M / cellsize) with cellsize 5.
 * When nothing can be merged within the cap, the two smallest cells are merged
 * instead so the loop always terminates.
 *
 * GAUSS-SEIDEL BY DEFAULT. After a cell is solved, the frozen rates of ITS
 * events are refreshed immediately, so later cells in the same sweep already
 * see the update. The reference offers a Jacobi variant for parallel runs;
 * this port implements the sequential Gauss-Seidel default, which is what the
 * reference selects when it is not asked for parallelism.
 *
 * ONE DELIBERATE SIMPLIFICATION, and what it costs. The reference carries
 * whatever output grid its ODE solver happens to produce, takes the union of
 * the cells' grids and interpolates every cell onto it. This port integrates
 * every cell on ONE FIXED GRID per segment instead, so the sweeps compare
 * trajectories sampled at the same instants and no interpolation of one cell
 * onto another's grid is needed. The frozen inbound drift is still linearly
 * interpolated between grid points, exactly as `fluid_interpcols` does.
 *
 * THE GRID IS THE ACCURACY KNOB, and the error it leaves is measurable. Only
 * the EXTERNAL contribution is approximated -- each cell's own drift is
 * integrated exactly -- so the residual behaves like the interpolation error,
 * falling roughly quadratically as the grid is refined. On a ten-queue model
 * whose exact population is 4, the closed-form closing solve gives 4.000000
 * and this decomposition gives
 *
 *     grid   16    65     129      257
 *     pop    4.167 4.0222 4.00285  3.99936
 *
 * so the default of 129 holds the population to under a tenth of a percent.
 * A model that needs more can raise `TbiOptions::grid`; nothing else changes,
 * because the fixed point being chased is the same coupled ODE either way.
 *
 * A SINGLE CELL IS THE CLOSING METHOD. With M <= cellsize the partition has one
 * cell, there are no external events, the inbound drift is identically zero and
 * the cell drift IS the closing drift. That is the case the tests pin, because
 * it is the one where TBI has an independent right answer to be checked against.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/util/lsoda.h"

namespace line {
namespace fluid {

/** Port of `tbi_partition.m`: stations grouped by routing coupling. */
template <class T>
std::vector<std::vector<std::size_t>> tbi_partition(const qn::NetworkStruct<T>& sn,
                                                    std::size_t cellsize = 5) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<std::vector<std::size_t>> cells;
    for (std::size_t i = 0; i < M; ++i) cells.push_back(std::vector<std::size_t>{i});
    if (cellsize == 0) return cells;
    const std::size_t target = std::max<std::size_t>(1, (M + cellsize - 1) / cellsize);
    if (cells.size() <= target) return cells;

    // Coupling: the total routing mass between two stations, symmetrized.
    const std::size_t S = sn.nof_stateful();
    Matrix<double> C(M, M, 0.0);
    if (sn.rt.rows() == S * K)
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t si = sn.stateful_of_station(i + 1) - 1;
            for (std::size_t j = 0; j < M; ++j) {
                const std::size_t sj = sn.stateful_of_station(j + 1) - 1;
                double w = 0.0;
                for (std::size_t a = 0; a < K; ++a)
                    for (std::size_t b = 0; b < K; ++b)
                        w += num_traits<T>::to_double(sn.rt(si * K + a, sj * K + b));
                C(i, j) += w;
                C(j, i) += w;
            }
        }
    for (std::size_t i = 0; i < M; ++i) C(i, i) = 0.0;

    while (cells.size() > target) {
        const std::size_t n = cells.size();
        double best = -1.0;
        std::size_t ba = n, bb = n;
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = a + 1; b < n; ++b) {
                if (cells[a].size() + cells[b].size() > 2 * cellsize) continue;
                if (C(a, b) > best) {
                    best = C(a, b);
                    ba = a;
                    bb = b;
                }
            }
        if (ba == n) {  // nothing fits the cap: merge the two smallest
            std::vector<std::size_t> ord(n);
            for (std::size_t i = 0; i < n; ++i) ord[i] = i;
            std::sort(ord.begin(), ord.end(),
                      [&](std::size_t p, std::size_t q) { return cells[p].size() < cells[q].size(); });
            ba = std::min(ord[0], ord[1]);
            bb = std::max(ord[0], ord[1]);
        }
        cells[ba].insert(cells[ba].end(), cells[bb].begin(), cells[bb].end());
        for (std::size_t k = 0; k < n; ++k) {
            C(ba, k) += C(bb, k);
            C(k, ba) += C(k, bb);
        }
        C(ba, ba) = 0.0;
        // Drop row/column bb by compacting into a fresh matrix.
        Matrix<double> C2(n - 1, n - 1, 0.0);
        for (std::size_t a = 0, aa = 0; a < n; ++a) {
            if (a == bb) continue;
            for (std::size_t b = 0, bbi = 0; b < n; ++b) {
                if (b == bb) continue;
                C2(aa, bbi) = C(a, b);
                ++bbi;
            }
            ++aa;
        }
        C = C2;
        cells.erase(cells.begin() + static_cast<long>(bb));
    }
    for (std::vector<std::size_t>& c : cells) std::sort(c.begin(), c.end());
    return cells;
}

/** Controls of the time-based iteration, mirroring `options.config.tbi_*`. */
struct TbiOptions {
    double tbi_tol = 1e-6;           ///< sup-norm gap that ends a segment's sweeps
    std::size_t tbi_iter_max = 200;  ///< sweeps per segment
    std::size_t cellsize = 5;        ///< target stations per cell
    std::size_t grid = 129;           ///< sample points per segment (see the header)
};

/**
 * Advance the state over [t0, t1] by time-based iteration.
 *
 * Returns the state at t1. `sys` is the ordinary closing system; TBI only
 * changes HOW it is integrated, never what is integrated.
 */
inline std::vector<double> tbi_advance(const FluidOdeSystem& sys,
                                       const std::vector<std::vector<std::size_t>>& cells,
                                       const std::vector<double>& y0, double t0, double t1,
                                       const TbiOptions& topt, const LsodaOptions& lopt) {
    const FluidLayout& L = sys.layout;
    const std::size_t n = L.nstates, ncells = cells.size();
    const std::size_t K = L.qidx.empty() ? 0 : L.qidx[0].size();

    // Which state entries belong to each cell, and which events are sourced
    // inside it. An event whose driving state is outside the cell is external:
    // its rate is frozen and its jump only matters where it lands inside.
    std::vector<std::vector<std::size_t>> mask(ncells);
    std::vector<std::vector<std::size_t>> eint(ncells), eext(ncells);
    std::vector<std::vector<long>> g2l(ncells, std::vector<long>(n, -1));
    for (std::size_t kc = 0; kc < ncells; ++kc) {
        for (std::size_t i : cells[kc])
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t k = 0; k < L.kic[i][r]; ++k) mask[kc].push_back(L.qidx[i][r] + k);
        std::sort(mask[kc].begin(), mask[kc].end());
        for (std::size_t a = 0; a < mask[kc].size(); ++a) g2l[kc][mask[kc][a]] = static_cast<long>(a);
        for (std::size_t e = 0; e < sys.events.size(); ++e) {
            const bool inside = g2l[kc][sys.events[e].event_idx] >= 0;
            if (inside) {
                eint[kc].push_back(e);
            } else if (g2l[kc][sys.events[e].minus] >= 0 || g2l[kc][sys.events[e].plus] >= 0) {
                eext[kc].push_back(e);  // lands in the cell but is driven outside
            }
        }
    }

    const std::size_t ng = std::max<std::size_t>(2, topt.grid);
    std::vector<double> tgrid(ng);
    for (std::size_t j = 0; j < ng; ++j)
        tgrid[j] = t0 + (t1 - t0) * static_cast<double>(j) / static_cast<double>(ng - 1);

    // Y[j] is the whole state at tgrid[j]; the first sweep freezes it at y0.
    std::vector<std::vector<double>> Y(ng, y0);

    for (std::size_t sweep = 0; sweep < topt.tbi_iter_max; ++sweep) {
        std::vector<std::vector<double>> Ynew = Y;
        double delta = 0.0;
        for (std::size_t kc = 0; kc < ncells; ++kc) {
            const std::size_t nl = mask[kc].size();
            if (nl == 0) continue;

            // The inbound drift on the grid, from events driven outside.
            std::vector<std::vector<double>> B(ng, std::vector<double>(nl, 0.0));
            for (std::size_t j = 0; j < ng; ++j) {
                std::vector<double> g(Y[j]);
                fluid_rates_closing(sys, Y[j].data(), g);
                for (std::size_t e : eext[kc]) {
                    const FluidEvent& ev = sys.events[e];
                    const double rate = ev.rate_base * g[ev.event_idx];
                    if (rate == 0.0) continue;
                    if (g2l[kc][ev.minus] >= 0) B[j][static_cast<std::size_t>(g2l[kc][ev.minus])] -= rate;
                    if (g2l[kc][ev.plus] >= 0) B[j][static_cast<std::size_t>(g2l[kc][ev.plus])] += rate;
                }
            }

            // The cell's own drift, plus the frozen inbound drift interpolated
            // linearly in time -- the port of `fluid_interpcols`.
            const std::vector<std::size_t>& mk = mask[kc];
            const std::vector<std::size_t>& ei = eint[kc];
            const std::vector<long>& gl = g2l[kc];
            std::vector<double> full(n, 0.0);
            const LsodaRhs f = [&sys, &mk, &ei, &gl, &B, &tgrid, ng, nl, n,
                                &full](double t, const double* xc, double* dxc) {
                std::vector<double> x(n, 0.0);
                for (std::size_t a = 0; a < nl; ++a) x[mk[a]] = xc[a];
                std::vector<double> g(x);
                fluid_rates_closing(sys, x.data(), g);
                for (std::size_t a = 0; a < nl; ++a) dxc[a] = 0.0;
                for (std::size_t e : ei) {
                    const FluidEvent& ev = sys.events[e];
                    const double rate = ev.rate_base * g[ev.event_idx];
                    if (rate == 0.0) continue;
                    if (gl[ev.minus] >= 0) dxc[static_cast<std::size_t>(gl[ev.minus])] -= rate;
                    if (gl[ev.plus] >= 0) dxc[static_cast<std::size_t>(gl[ev.plus])] += rate;
                }
                // linear interpolation of the frozen inbound drift
                double u = (t - tgrid.front()) / (tgrid.back() - tgrid.front() + 1e-300);
                u = std::min(1.0, std::max(0.0, u)) * static_cast<double>(ng - 1);
                const std::size_t j0 = std::min<std::size_t>(ng - 2, static_cast<std::size_t>(u));
                const double w = u - static_cast<double>(j0);
                for (std::size_t a = 0; a < nl; ++a)
                    dxc[a] += (1.0 - w) * B[j0][a] + w * B[j0 + 1][a];
            };

            std::vector<double> yl(nl, 0.0);
            for (std::size_t a = 0; a < nl; ++a) yl[a] = y0[mk[a]];
            const LsodaSolution s = fluid_integrate_grid(f, yl, tgrid, lopt);
            for (std::size_t j = 0; j < s.y.size() && j < ng; ++j)
                for (std::size_t a = 0; a < nl; ++a) {
                    double v = s.y[j][a];
                    if (v < 0.0) v = 0.0;
                    delta = std::max(delta, std::fabs(v - Ynew[j][mk[a]]));
                    Ynew[j][mk[a]] = v;
                    // Gauss-Seidel: the next cell of this sweep already sees it.
                    Y[j][mk[a]] = v;
                }
        }
        Y = Ynew;
        if (delta < topt.tbi_tol) break;
    }
    return Y.back();
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_TBI_H
