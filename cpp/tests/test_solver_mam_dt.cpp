/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverMAM on a discrete (slotted) time scale, reached through the ordinary
 * dispatch: the model is recognized as discrete-time from its distributions
 * alone. Oracles are the Geo/Geo/1 closed form under LAS_DA and the MATLAB
 * SolverMAM numbers recorded in _kb/06-solver-catalog.md, which the LDES
 * slotted engine independently corroborated.
 */
#include <cmath>
#include <string>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_runner.h"

using line::Matrix;
using line::lang::Distrib;
using line::lang::SchedStrategy;
using line::mam::MamOptions;
using line::mam::MamSolution;
using line::mam::solver_mam_solve;
using line::qn::Network;
using line::qn::RoutingMatrix;

namespace {

/** Source -> FCFS Queue -> Sink with the given arrival and service laws. */
line::qn::NetworkStruct<double> single_queue(const Distrib<double>& arv,
                                             const Distrib<double>& svc) {
    Network<double> b("dt1");
    const std::size_t src = b.add_source("Source");
    const std::size_t q = b.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = b.add_sink("Sink");
    const std::size_t c = b.add_open_class("Class1");
    b.set_arrival(src, c, arv);
    b.set_service(q, c, svc);
    RoutingMatrix<double> P = b.init_routing_matrix();
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    b.link(P);
    return b.get_struct();
}

MamSolution<double> solve(const line::qn::NetworkStruct<double>& sn) {
    return solver_mam_solve(sn, MamOptions());
}

}  // namespace

TEST_CASE("SolverMAM recognizes a slotted model and solves Geo/Geo/1 exactly") {
    const double a = 0.2, s = 0.5;
    MamSolution<double> r =
        solve(single_queue(Distrib<double>::geometric(a), Distrib<double>::geometric(s)));

    CHECK(r.actualmethod == "dt.qmam");
    // E[N] = A(1-A)/(S-A), U = A/S under LAS-DA
    CHECK(r.sol.Q(1, 0) == doctest::Approx(a * (1 - a) / (s - a)).epsilon(1e-6));
    CHECK(r.sol.U(1, 0) == doctest::Approx(a / s).epsilon(1e-6));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(a).epsilon(1e-6));
    CHECK(r.sol.Tp(0, 0) == doctest::Approx(a).epsilon(1e-6));
}

TEST_CASE("a Det service on the slot lattice stays on the discrete path") {
    // MATLAB SolverMAM reports 0.466667 here, LDES slotted 0.466615
    MamSolution<double> r =
        solve(single_queue(Distrib<double>::geometric(0.2), Distrib<double>::det(2.0)));
    CHECK(r.actualmethod == "dt.qmam");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(0.466667).epsilon(1e-5));
    CHECK(r.sol.U(1, 0) == doctest::Approx(0.4).epsilon(1e-6));
}

TEST_CASE("the slotted tandem reproduces the discrete-time Burke result") {
    // The stationary departure stream of a Geo/Geo/1 queue is Bernoulli, so the
    // second queue is EXACT: 0.8 = A(1-A)/(S-A) with A = 0.2, S = 0.4.
    // MATLAB reports 0.5333/0.8000, LDES slotted 0.5363/0.8038.
    Network<double> b("dttandem");
    const std::size_t src = b.add_source("Source");
    const std::size_t q1 = b.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = b.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t snk = b.add_sink("Sink");
    const std::size_t c = b.add_open_class("Class1");
    b.set_arrival(src, c, Distrib<double>::geometric(0.2));
    b.set_service(q1, c, Distrib<double>::geometric(0.5));
    b.set_service(q2, c, Distrib<double>::geometric(0.4));
    RoutingMatrix<double> P = b.init_routing_matrix();
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, snk, 1.0);
    b.link(P);

    MamSolution<double> r = solve(b.get_struct());
    CHECK(r.actualmethod == "dt.dec");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(0.533333).epsilon(1e-4));
    CHECK(r.sol.Q(2, 0) == doctest::Approx(0.800000).epsilon(1e-3));
    CHECK(r.sol.U(1, 0) == doctest::Approx(0.4).epsilon(1e-6));
    CHECK(r.sol.U(2, 0) == doctest::Approx(0.5).epsilon(1e-6));
}

TEST_CASE("a continuous model is not swallowed by the discrete-time detection") {
    // An exponential Source and service must not be read as slotted; M/M/1 with
    // rho = 0.4 has E[N] = rho/(1-rho)
    MamSolution<double> r =
        solve(single_queue(Distrib<double>::exp_rate(0.2), Distrib<double>::exp_rate(0.5)));
    CHECK(r.actualmethod != "dt.qmam");
    CHECK(r.actualmethod != "dt.dec");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(0.4 / 0.6).epsilon(1e-4));
}
