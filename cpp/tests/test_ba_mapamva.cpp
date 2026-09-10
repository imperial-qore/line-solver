/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * MAP-AMVA LP bounds (Casale-Smirni, IEEE/IFIP DSN 2009) in SolverBA.
 *
 * THIS IS THE ONLY SolverBA FAMILY DERIVED FOR A CORRELATED SERVICE PROCESS.
 * Every other family is a function of the demands D = V./rates, so a MAP is
 * indistinguishable there from the exponential of the same mean and the bracket
 * returned is for a DIFFERENT system; MAP-AMVA is written over the exact
 * mean-value balances of the MAP network itself, in the per-phase variables
 * QN(i,k) and UN(i,k). These tests pin the two consequences: the bracket really
 * does contain the exact solution of a MAP model, and the feature gate lets a
 * MAP reach this family and no other.
 *
 * THE ORACLE IS CONTAINMENT, NOT CLOSENESS. An LP relaxation bound is only
 * required to bracket, so every check here is an inequality against SolverCTMC
 * rather than an equality; asserting a value would pin the simplex vertex a
 * particular solver happens to reach.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
const double kTol = 1e-9;

/** D0 + D1 is a proper generator; D1 has an off-diagonal entry, so services correlate. */
Distrib<double> corr_map() {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -3.0;
    D0(0, 1) = 0.5;
    D0(1, 0) = 0.2;
    D0(1, 1) = -0.4;
    D1(0, 0) = 2.5;
    D1(1, 0) = 0.2;
    return D::map_dist(D0, D1, lang::ProcessType::MAP);
}

qn::Network<double> tandem(bool map_first, double N) {
    qn::Network<double> m("mapamvaTandem");
    const std::size_t q0 = m.add_queue("Q0", SchedStrategy::FCFS);
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", N, q0);
    if (map_first) {
        m.set_service(q0, c, corr_map());
        m.set_service(q1, c, D::exp_rate(2.0));
    } else {
        m.set_service(q0, c, D::exp_rate(2.0));
        m.set_service(q1, c, corr_map());
    }
    qn::RoutingMatrix<double> P;
    P.set(q0, q1, 1.0);
    P.set(q1, q0, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run_ba(qn::Network<double>& m, const std::string& method) {
    ba::BaOptions opt;
    opt.method = method;
    return ba::solver_ba_run_analyzer(m.get_struct(), opt);
}

void check_brackets(qn::Network<double>& m) {
    const mva::AvgResult<double> lo = run_ba(m, "mapamva.lower");
    const mva::AvgResult<double> up = run_ba(m, "mapamva.upper");
    const mva::AvgResult<double> ex =
        ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());
    for (std::size_t i = 0; i < ex.QN.rows(); ++i) {
        for (std::size_t r = 0; r < ex.QN.cols(); ++r) {
            CHECK(lo.QN(i, r) <= ex.QN(i, r) + kTol);
            CHECK(ex.QN(i, r) <= up.QN(i, r) + kTol);
            CHECK(lo.UN(i, r) <= ex.UN(i, r) + kTol);
            CHECK(ex.UN(i, r) <= up.UN(i, r) + kTol);
            CHECK(lo.RN(i, r) <= ex.RN(i, r) + kTol);
            CHECK(ex.RN(i, r) <= up.RN(i, r) + kTol);
            CHECK(lo.TN(i, r) <= ex.TN(i, r) + kTol);
            CHECK(ex.TN(i, r) <= up.TN(i, r) + kTol);
        }
    }
}

}  // namespace

TEST_CASE("ba mapamva brackets a map at the last station") {
    qn::Network<double> m = tandem(false, 8.0);
    check_brackets(m);
}

TEST_CASE("ba mapamva brackets a map at the first station") {
    // The LP requires the phase-carrying queue to be index M; the analyzer
    // PERMUTES it there and inverts the permutation before reporting, so the
    // bracket must hold with the stations the other way round too.
    qn::Network<double> m = tandem(true, 8.0);
    check_brackets(m);
}

TEST_CASE("ba mapamva degenerates to one level with no phase-carrying station") {
    // K collapses to 1 and the balances become the product-form ones. Still a
    // bracket, and the exact solution is in it.
    qn::Network<double> m("mapamvaExp");
    const std::size_t q0 = m.add_queue("Q0", SchedStrategy::FCFS);
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 5.0, q0);
    m.set_service(q0, c, D::exp_rate(2.0));
    m.set_service(q1, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(q0, q1, 1.0);
    P.set(q1, q0, 1.0);
    m.link(P);
    check_brackets(m);
}

TEST_CASE("ba mapamva sides are ordered") {
    qn::Network<double> m = tandem(false, 8.0);
    const mva::AvgResult<double> lo = run_ba(m, "mapamva.lower");
    const mva::AvgResult<double> up = run_ba(m, "mapamva.upper");
    for (std::size_t i = 0; i < lo.QN.rows(); ++i) {
        CHECK(lo.QN(i, 0) <= up.QN(i, 0) + kTol);
        CHECK(lo.UN(i, 0) <= up.UN(i, 0) + kTol);
        CHECK(lo.RN(i, 0) <= up.RN(i, 0) + kTol);
        CHECK(lo.TN(i, 0) <= up.TN(i, 0) + kTol);
    }
}

TEST_CASE("ba mapamva is offered on a closed single-server model") {
    qn::Network<double> m = tandem(false, 4.0);
    const std::vector<std::string> listed = ba::list_valid_methods(m.get_struct());
    CHECK(std::find(listed.begin(), listed.end(), "mapamva.upper") != listed.end());
    CHECK(std::find(listed.begin(), listed.end(), "mapamva.lower") != listed.end());
}

TEST_CASE("ba mapamva refuses a delay station") {
    // The LP is a network of queues, and Casale-Smirni name the delay extension
    // as open work.
    qn::Network<double> m("mapamvaDelay");
    const std::size_t d0 = m.add_delay("D0");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 4.0, d0);
    m.set_service(d0, c, D::exp_rate(1.0));
    m.set_service(q1, c, corr_map());
    qn::RoutingMatrix<double> P;
    P.set(d0, q1, 1.0);
    P.set(q1, d0, 1.0);
    m.link(P);
    const std::vector<std::string> listed = ba::list_valid_methods(m.get_struct());
    CHECK(std::find(listed.begin(), listed.end(), "mapamva.upper") == listed.end());
    CHECK_THROWS(run_ba(m, "mapamva.upper"));
}

TEST_CASE("ba mapamva refuses two phase-carrying stations") {
    // The LP gives queue M the (D0,D1) pair and every other queue a SCALAR rate,
    // so it has nowhere to put a second phase-type station.
    qn::Network<double> m("mapamvaTwoPh");
    const std::size_t q0 = m.add_queue("Q0", SchedStrategy::FCFS);
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 4.0, q0);
    m.set_service(q0, c, D::erlang(4.0, 2));
    m.set_service(q1, c, corr_map());
    qn::RoutingMatrix<double> P;
    P.set(q0, q1, 1.0);
    P.set(q1, q0, 1.0);
    m.link(P);
    CHECK_THROWS(run_ba(m, "mapamva.upper"));
}
