/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/advanced/largeScale/`, `python/examples/advanced/largeScale/`:
 * transient fluid analysis of a model too large for the monolithic arm.
 *
 * THE PARAMETERS ARE DRAWN, AND EACH CODEBASE DRAWS ITS OWN. The reference
 * seeds MATLAB's twister, the Python twin seeds numpy's PCG64, and neither
 * stream is the other's -- so the three editions already build three different
 * networks of the same SHAPE, and this one draws from `std::mt19937` rather
 * than reproducing a stream that is not shared anyway. Nothing here is
 * goldened for that reason: what the example demonstrates is how TBI scales on
 * a network of this size, and that is a property of the shape.
 *
 * The reference's explicit cell partition (`options.config.tbi_cells`: one
 * queueing station together with its outbound transit delays) is passed as
 * `FluidOptions::tbi_cells`, in 0-based station indices.
 */

#include <cstdio>
#include <ctime>
#include <random>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/fluid/fluid_runner.h"

namespace line {
namespace examples {

namespace {

const std::size_t kM = 10;     ///< queueing stations; the model holds M + M*(M-1)
const double kN = 250.0;       ///< closed population (vehicles)
const std::size_t kKph = 16;   ///< Erlang phases per transit delay
const double kTend = 20.0;     ///< transient horizon

}  // namespace

/**
 * Large-scale transient fluid analysis with trajectory-based iteration (TBI).
 *
 * A vehicle-sharing-style closed network with M queueing stations and one
 * Erlang transit delay per ordered station pair, giving O(M^2) stations and,
 * with 16 Erlang phases per delay, a fluid ODE with about 1500 state variables.
 * At this size the monolithic stiff fluid solution (method 'closing') takes
 * several minutes, dominated by the Jacobian factorizations, while
 * trajectory-based iteration solves each station cell separately against frozen
 * inbound trajectories and completes in seconds.
 *
 * Reference: M. Sheldon, D. Tuncer, G. Casale, "TBI: Transient Hierarchical
 * Modeling of Large-Scale Vehicle Sharing Systems", IEEE Transactions on
 * Intelligent Transportation Systems.
 */
void largescale_tbi() {
    Net m("tbi_largescale");
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<std::size_t> Q;
    for (std::size_t i = 0; i < kM; ++i)
        Q.push_back(m.add_queue("Q" + std::to_string(i + 1), SchedStrategy::PS));
    std::vector<std::vector<std::size_t> > Dl(kM, std::vector<std::size_t>(kM, 0));
    for (std::size_t i = 0; i < kM; ++i)
        for (std::size_t j = 0; j < kM; ++j)
            if (i != j)
                Dl[i][j] = m.add_delay("D" + std::to_string(i + 1) + "_" + std::to_string(j + 1));

    const std::size_t job = m.add_closed_class("C1", kN, Q[0]);
    Matrix<double> lambda(kM, kM, 0.0);
    for (std::size_t i = 0; i < kM; ++i) {
        m.set_service(Q[i], job,
                      lang::erlang_fit_mean_order<double>(1.0 / (1.0 + 3.0 * u01(rng)), 4));
        for (std::size_t j = 0; j < kM; ++j)
            if (i != j) {
                m.set_service(Dl[i][j], job,
                              lang::erlang_fit_mean_order<double>(0.2 + 2.0 * u01(rng), kKph));
                lambda(i, j) = u01(rng);
            }
    }

    Routing P;
    for (std::size_t i = 0; i < kM; ++i) {
        double tot = 0.0;
        for (std::size_t j = 0; j < kM; ++j)
            if (i != j) tot += lambda(i, j);
        for (std::size_t j = 0; j < kM; ++j)
            if (i != j) {
                P.set(job, job, Q[i], Dl[i][j], lambda(i, j) / tot);
                P.set(job, job, Dl[i][j], Q[j], 1.0);
            }
    }
    m.link(P);
    const Sn& sn = m.get_struct();

    // One cell per queueing station plus its outbound transit delays, as 0-based
    // station indices.
    std::vector<std::vector<std::size_t> > cells(kM);
    for (std::size_t i = 0; i < kM; ++i) {
        cells[i].push_back(m.station_index(Q[i]) - 1);
        for (std::size_t j = 0; j < kM; ++j)
            if (i != j) cells[i].push_back(m.station_index(Dl[i][j]) - 1);
    }

    fluid::FluidOptions opt;
    opt.method = "tbi";
    opt.tbi_cells = cells;
    opt.timespan_end = kTend;
    opt.stiff = true;
    const std::clock_t t0 = std::clock();
    const fluid::FluidSolution sol = fluid::solver_fluid_run_analyzer(sn, opt);
    const double secs = double(std::clock() - t0) / double(CLOCKS_PER_SEC);
    std::printf("TBI solved %zu stations (%zu ODE variables) in %.1f seconds.\n", sn.nstations,
                kM * (kM - 1) * kKph + kM * 4, secs);

    // The queue-length summary at the queueing stations, which are the first M
    // rows: the transit delays carry the rest of the population and are not
    // what the paper's tables report.
    std::printf("%-10s %12s %12s %12s\n", "Station", "QLen", "Util", "Tput");
    for (std::size_t i = 0; i < kM; ++i)
        std::printf("%-10s %12.5g %12.5g %12.5g\n", sn.stations[i].name.c_str(), sol.QN(i, 0),
                    sol.UN(i, 0), sol.TN(i, 0));
}

LINE_EXAMPLE("advanced/largeScale", largescale_tbi);

}  // namespace examples
}  // namespace line
