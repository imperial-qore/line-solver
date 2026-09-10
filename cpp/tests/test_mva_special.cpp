/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The two analyzers the dispatch reaches on a model no product-form path can
 * solve: the closed LCFS + LCFS-PR pair, and blocking-after-service.
 *
 * Both are reached WITHOUT naming a method, which is the part worth testing:
 * LCFS-QN sits inside `solver_mva` ahead of its product-form test, and SQD is a
 * branch of the `default` ladder keyed on the drop rule. A model that reaches
 * the wrong one still returns numbers, so the reported method is asserted too.
 *
 * Numbers are MATLAB's `SolverMVA(model).getAvg()`.
 */

#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::DropStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

TEST_CASE("a closed LCFS + LCFS-PR pair takes the specialized recursion") {
    qn::Network<double> m("lcfsqn");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::LCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::LCFSPR);
    const std::size_t c = m.add_closed_class("C1", 3.0, q1);
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    // the reference reports this as exact: the gate lives inside solver_mva
    CHECK(r.actualmethod == "exact");
    CHECK(r.QN(0, 0) == doctest::Approx(1.984615384615).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(1.015384615385).epsilon(1e-9));
    CHECK(r.UN(0, 0) == doctest::Approx(0.876923076923).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.584615384615).epsilon(1e-9));
    CHECK(r.RN(0, 0) == doctest::Approx(1.131578947368).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.578947368421).epsilon(1e-9));
    CHECK(r.TN(0, 0) == doctest::Approx(1.753846153846).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(1.753846153846).epsilon(1e-9));
    // the population is conserved across the two stations
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("an LCFS station without its LCFS-PR partner is refused, not approximated") {
    qn::Network<double> m("lcfs-alone");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::LCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    opt.method = "exact";
    Matrix<double> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
}

TEST_CASE("blocking after service reaches SQD from the default ladder") {
    qn::Network<double> m("bas");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t b1 = m.add_queue("B1", SchedStrategy::FCFS);
    const std::size_t b2 = m.add_queue("B2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(b1, c, D::exp_rate(2.0));
    m.set_service(b2, c, D::exp_rate(3.0));
    m.set_capacity(b1, 2.0);
    m.set_drop_rule(b1, c, DropStrategy::BAS);
    m.set_capacity(b2, 2.0);
    m.set_drop_rule(b2, c, DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(d, b1, 1.0);
    P.set(b1, b2, 1.0);
    P.set(b2, d, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    CHECK(r.actualmethod == "sqd");
    CHECK(r.QN(0, 0) == doctest::Approx(1.503416595735).epsilon(1e-8));
    CHECK(r.QN(1, 0) == doctest::Approx(1.572254991721).epsilon(1e-8));
    CHECK(r.QN(2, 0) == doctest::Approx(0.924328412544).epsilon(1e-8));
    CHECK(r.UN(1, 0) == doctest::Approx(0.751708297868).epsilon(1e-8));
    CHECK(r.UN(2, 0) == doctest::Approx(0.501138865245).epsilon(1e-8));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(r.RN(1, 0) == doctest::Approx(1.045787971332).epsilon(1e-8));
    CHECK(r.RN(2, 0) == doctest::Approx(0.614818550737).epsilon(1e-8));
    CHECK(r.TN(0, 0) == doctest::Approx(1.503416595735).epsilon(1e-8));
    CHECK(r.TN(2, 0) == doctest::Approx(1.503416595735).epsilon(1e-8));
    // the blocked jobs are still in the network, so the population is conserved
    CHECK(r.QN(0, 0) + r.QN(1, 0) + r.QN(2, 0) == doctest::Approx(4.0).epsilon(1e-8));

    // and the method is reachable by name on this model
    CHECK(mva::list_valid_methods(m.get_struct()).size() > 0);
    mva::MvaOptions o2;
    o2.method = "sqd";
    const mva::AvgResult<double> r2 = mva::solver_mva_run_analyzer(m.get_struct(), o2, init);
    CHECK(r2.QN(1, 0) == doctest::Approx(r.QN(1, 0)).epsilon(1e-12));
}

}  // namespace
