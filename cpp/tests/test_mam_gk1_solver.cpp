/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverMAM answers a station with NON-PHASE-TYPE, class-dependent service
 * through MMAP[K]/G[K]/1 rather than through the phase-type fit.
 *
 * The generic path reads the service law out of the PH fit, which for a Uniform
 * is a 20-phase approximation: it matches the mean and, above SCV 1, nothing
 * else. He's transform analysis takes the original law. THE ORACLE IS THE MATLAB
 * REFERENCE, which reaches the same numbers through sn.lst, and JMT at 8e6
 * samples agrees with both to 5e-4 per class.
 */
#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mam/solver_mam_runner.h"

namespace qn = line::qn;
namespace mam = line::mam;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

TEST_CASE("mam open: general class-dependent service takes the transform analysis") {
    qn::Network<double> m("gk1");
    const std::size_t src = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("K");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dist::exp_rate(0.35));
    {
        // MMPP2(0.9, 0.2, 0.25, 0.35)
        line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
        D0(0, 0) = -(0.9 + 0.25);
        D0(0, 1) = 0.25;
        D0(1, 0) = 0.35;
        D0(1, 1) = -(0.2 + 0.35);
        D1(0, 0) = 0.9;
        D1(1, 1) = 0.2;
        m.set_arrival(src, c2, Dist::map_dist(D0, D1, line::lang::ProcessType::MMPP2));
    }
    m.set_service(q, c1, Dist::det(0.5));
    m.set_service(q, c2, Dist::uniform(0.2, 1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);

    const line::mva::MvaSolution<double> s = mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
    CHECK(s.R(1, 0) == doctest::Approx(0.936113).epsilon(1e-4));
    CHECK(s.R(1, 1) == doctest::Approx(1.121427).epsilon(1e-4));
}

TEST_CASE("mam open: the declared law survives the phase-type conversion") {
    // Through the RUNNER, not solver_mam_basic: sn_nonmarkov_toph runs there and
    // replaces the Gamma by a phase-type fit, retagging it PH. The case above
    // never sees that, so it held while a Gamma on its own fell to the fit and
    // read 1.0785538590. Pollaczek-Khinchine needs only the first two moments of
    // an M/G/1 service, so it is an oracle independent of every path here.
    qn::Network<double> m("declared");
    const std::size_t src = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("K");
    const std::size_t c1 = m.add_open_class("C1");
    m.set_arrival(src, c1, Dist::exp_rate(0.5));
    m.set_service(q, c1, Dist::gamma_dist(1.0 / 3.0, 3.0));  // mean 1, SCV 3
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    m.link(P);

    const mam::MamSolution<double> s = mam::solver_mam_solve(m.get_struct(), mam::MamOptions());
    CHECK(s.sol.Q(1, 0) == doctest::Approx(1.5).epsilon(1e-6));
}
