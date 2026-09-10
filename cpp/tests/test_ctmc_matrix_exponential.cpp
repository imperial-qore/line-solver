/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverCTMC on an EXPLICITLY DECLARED matrix exponential, which is a different
 * case from `test_ctmc_fluid_nonmarkov.cpp`. There the model carries a Gamma and
 * `sn_nonmarkov_toph` chooses the surrogate, so the port forces a phase-type fit
 * and the ME never reaches the generator. Here the caller declares the ME
 * itself: it embeds in the generator exactly as a phase-type does, negative
 * off-diagonal entries and all, and the stationary vector is then a genuinely
 * SIGNED measure whose aggregates over each phase block are still the exact
 * probabilities. Mean measures are linear in that vector and stay exact.
 *
 * WHAT THIS FILE PINS is the one step that is not linear: the clamp
 * `solver_ctmc_avg_from_pi` applies before renormalizing. Deleting the negative
 * entries deletes real mass, and it does so QUIETLY -- M/CME/1 at rho 0.3 came
 * back with QLen 0.4469 against the Pollaczek-Khinchine 0.3772, Util 0.7108 for
 * a station whose true utilization is rho, and a departure rate that did not
 * match the arrival rate. Pollaczek-Khinchine is the oracle because it fixes the
 * mean queue length from the SCV alone, so a fit that is not the declared law
 * cannot pass it.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/cme.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace qn = line::qn;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> FCFS Queue -> Sink, Poisson(rho) arrivals, unit-mean CME service. */
qn::Network<double> mme1(double rho, std::size_t order) {
    const line::mam::CmeRepresentation<double> rep = line::mam::cme_representation<double>(order);
    qn::Network<double> m("M/CME/1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, Dist::exp_rate(rho));
    m.set_service(q, c, Dist::me(rep.alpha, rep.A));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("SolverCTMC reproduces Pollaczek-Khinchine on an M/CME/1") {
    // Order 7 at cutoff 30 is 211 states; the smaller order keeps the case cheap
    // while still carrying an oscillating ME (order 3 already has a rotation
    // block, hence the negative off-diagonal entries this test is about).
    struct Case {
        double rho;
        std::size_t order;
        std::size_t cutoff;
    };
    const Case cases[] = {{0.3, 3, 30}, {0.5, 7, 30}};

    for (std::size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        const Case& tc = cases[k];
        qn::Network<double> m = mme1(tc.rho, tc.order);
        line::ctmc::CtmcOptions o;
        o.cutoff = static_cast<double>(tc.cutoff);
        const line::ctmc::CtmcSolution<double> s =
            line::ctmc::solver_ctmc_analyzer(m.get_struct(), o);

        const double scv = line::mam::cme_min_scv(tc.order);
        const double pk =
            tc.rho + tc.rho * tc.rho * (1.0 + scv) / (2.0 * (1.0 - tc.rho));

        // Station 1 is the Queue; station 0 is the Source, whose queue length is
        // an infinite reservoir rather than a number.
        CHECK(s.avg.QN(1, 0) == doctest::Approx(pk).epsilon(1e-6));

        // The two properties the clamp broke, and which no tolerance on QLen
        // alone would have caught: utilization IS rho for a unit-mean service,
        // and what leaves the queue is what entered the system.
        CHECK(s.avg.UN(1, 0) == doctest::Approx(tc.rho).epsilon(1e-6));
        CHECK(s.avg.TN(1, 0) == doctest::Approx(tc.rho).epsilon(1e-6));
    }
}

TEST_CASE("the phase-type predicate separates an ME from a Coxian") {
    // ctmc_all_phasetype is what decides both the clamp and the refusal of the
    // per-state getters, so it is pinned on both sides rather than only where it
    // returns false.
    const line::mam::CmeRepresentation<double> rep = line::mam::cme_representation<double>(3);
    qn::Network<double> me = mme1(0.3, 3);
    CHECK_FALSE(line::ctmc::ctmc_all_phasetype(me.get_struct()));

    qn::Network<double> ph("M/Cox2/1");
    const std::size_t src = ph.add_source("Source");
    const std::size_t q = ph.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = ph.add_sink("Sink");
    const std::size_t c = ph.add_open_class("Class1");
    ph.set_arrival(src, c, Dist::exp_rate(0.3));
    ph.set_service(q, c, Dist::cox2(2.0, 4.0, 0.5));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    ph.link(P);
    CHECK(line::ctmc::ctmc_all_phasetype(ph.get_struct()));

    // Guards against the predicate degenerating to "the model has an ME
    // somewhere": the representation must be the reason, not the model name.
    CHECK(rep.A.rows() == 3);
}
