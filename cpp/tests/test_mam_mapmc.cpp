/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverMAM answers a multiserver station fed by a CORRELATED arrival stream
 * with the exact MAP/M/c queue, not with the single-fast-server surrogate.
 *
 * The PH/M/c fast path gates on the arrival being RENEWAL, because qsys_phmc
 * reads the arrival through its (pie, D0) marginal alone. An MMPP2 source is
 * not renewal, so it fell through to the generic path, which divides the
 * service time by the server count and adds a surrogate delay -- an
 * approximation that discards the arrival correlation entirely. On this model
 * that read 1.757271 against the exact 1.521875, 15.5% high.
 *
 * THE ORACLE IS AN INDEPENDENT ALGORITHM, not recorded output: qsys_mapmc
 * solves the level-dependent QBD directly, and the MATLAB reference reports the
 * same 1.521875 through its own Q-MAM path.
 */
#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_basic.h"

namespace qn = line::qn;
namespace mam = line::mam;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> Queue(c servers, Exp) -> Sink, with MMPP2(1.8, 0.4, 0.15, 0.25) arrivals. */
qn::Network<double> mapmc_model(double c) {
    qn::Network<double> m("mapmc");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Snk");
    const std::size_t cls = m.add_open_class("C1");
    line::Matrix<double> D0(2, 2), D1(2, 2);
    D0(0, 0) = -(1.8 + 0.15);
    D0(0, 1) = 0.15;
    D0(1, 0) = 0.25;
    D0(1, 1) = -(0.4 + 0.25);
    D1(0, 0) = 1.8;
    D1(0, 1) = 0.0;
    D1(1, 0) = 0.0;
    D1(1, 1) = 0.4;
    m.set_arrival(src, cls, Dist::map_dist(D0, D1, line::lang::ProcessType::MMPP2));
    m.set_service(q, cls, Dist::exp_rate(1.0));
    m.set_number_of_servers(q, c);
    qn::RoutingMatrix<double> P;
    P.set(cls, cls, src, q, 1.0);
    P.set(cls, cls, q, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("mam open: a correlated arrival at c > 1 servers takes the exact MAP/M/c path") {
    qn::Network<double> m = mapmc_model(3.0);
    const line::mva::MvaSolution<double> s = mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
    // station 1 is the queue; the MATLAB reference and the JAR both report this
    CHECK(s.Q(1, 0) == doctest::Approx(1.521875).epsilon(1e-6));
}

TEST_CASE("mam open: one server is unaffected, the exact MAP/MAP/1 path still owns it") {
    qn::Network<double> m = mapmc_model(1.0);
    const line::mva::MvaSolution<double> s = mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
    // MAP/M/1 with the same stream: strictly more work than the 3-server case
    CHECK(s.Q(1, 0) > 1.521875);
}
