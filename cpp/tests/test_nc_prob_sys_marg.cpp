/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `@@SolverNC/getProbSysMarg.m`: the joint law of the per-station TOTAL queue
 * lengths, evaluated through the Calame permanent identity.
 *
 * WHY THE CHECKS RUN IN THIS ORDER. Normalization alone would pass on a law
 * that is uniformly wrong by a constant, so the FIRST MOMENT of the same sweep
 * is compared against the exact CTMC queue lengths: a normalizing bug and a
 * per-state bug are distinguishable only when both are tested. The zero-bearing
 * demand matrix and the two-delay model are then run through the same pair,
 * because each breaks a different part of the identity -- a structural zero
 * kills the approximate estimators' support, and a second infinite server needs
 * its OWN 1/n_j! that a single aggregated delay row would not supply.
 *
 * The refusals are asserted by message, not merely by "something was raised":
 * the point of refusing a multiserver station and a zero-bearing approximate
 * engine is that the caller is told which one, and a path that throws anything
 * at all is indistinguishable from a numerical accident.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/nc/solver_nc_prob.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Delay + two PS queues, 2 classes, N = (2,1). Dense demand matrix. */
qn::Network<double> dense_model() {
    qn::Network<double> m("marg_dense");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("Class1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("Class2", 1.0, d);
    m.set_service(d, c1, D::exp_rate(1.0 / 0.7));
    m.set_service(d, c2, D::exp_rate(1.0 / 1.3));
    m.set_service(q1, c1, D::exp_rate(1.0 / 1.5));
    m.set_service(q1, c2, D::exp_rate(1.0 / 0.8));
    m.set_service(q2, c1, D::exp_rate(1.0 / 0.9));
    m.set_service(q2, c2, D::exp_rate(1.0 / 1.1));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, d, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

/** The same three stations, but class 2 never visits Queue2. */
qn::Network<double> zero_demand_model() {
    qn::Network<double> m("marg_zero");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("Class1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("Class2", 1.0, d);
    m.set_service(d, c1, D::exp_rate(1.0 / 0.7));
    m.set_service(d, c2, D::exp_rate(1.0 / 1.3));
    m.set_service(q1, c1, D::exp_rate(1.0 / 1.5));
    m.set_service(q1, c2, D::exp_rate(1.0 / 0.8));
    m.set_service(q2, c1, D::exp_rate(1.0 / 0.9));
    m.set_service(q2, c2, D::exp_rate(1.0 / 1.1));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, d, 1.0);
    // Class 2 bypasses Queue2 entirely, so its demand there is a structural zero.
    P.set(c2, c2, d, q1, 1.0);
    P.set(c2, c2, q1, d, 1.0);
    m.link(P);
    return m;
}

/** TWO infinite servers plus one queue: each delay carries its own 1/n_j!. */
qn::Network<double> two_delay_model() {
    qn::Network<double> m("marg_twodelay");
    const std::size_t d1 = m.add_delay("Delay1");
    const std::size_t d2 = m.add_delay("Delay2");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("Class1", 2.0, d1);
    const std::size_t c2 = m.add_closed_class("Class2", 1.0, d1);
    m.set_service(d1, c1, D::exp_rate(1.0 / 0.7));
    m.set_service(d1, c2, D::exp_rate(1.0 / 1.3));
    m.set_service(d2, c1, D::exp_rate(1.0 / 1.1));
    m.set_service(d2, c2, D::exp_rate(1.0 / 0.6));
    m.set_service(q, c1, D::exp_rate(1.0 / 1.5));
    m.set_service(q, c2, D::exp_rate(1.0 / 0.8));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, d1, d2, 1.0);
        P.set(c, c, d2, q, 1.0);
        P.set(c, c, q, d1, 1.0);
    }
    m.link(P);
    return m;
}

/** Sum of the law and its induced mean queue length per station. */
void law_moments(const qn::NetworkStruct<double>& sn, std::size_t M, int total, double* mass,
                 std::vector<double>* mean) {
    const nc::NcSolverOptions opt;
    const std::vector<std::vector<int> > states =
        pfqn::multichoose_rows(static_cast<int>(M), total);
    *mass = 0.0;
    mean->assign(M, 0.0);
    for (std::size_t j = 0; j < states.size(); ++j) {
        const double p = nc::solver_nc_getprob_sys_marg(sn, opt, states[j]);
        CHECK(p >= -1e-12);
        *mass += p;
        for (std::size_t i = 0; i < M; ++i) (*mean)[i] += p * states[j][i];
    }
}

/** Station totals from the exact CTMC solution, classes summed out. */
std::vector<double> ctmc_totals(qn::Network<double>& m) {
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());
    std::vector<double> out(r.QN.rows(), 0.0);
    for (std::size_t i = 0; i < r.QN.rows(); ++i)
        for (std::size_t c = 0; c < r.QN.cols(); ++c) out[i] += r.QN(i, c);
    return out;
}

/** Rows as plain integers, so an expected set can be written literally. */
std::vector<std::vector<int> > as_ints(const std::vector<std::vector<double> >& rows) {
    std::vector<std::vector<int> > out;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        std::vector<int> r;
        for (std::size_t j = 0; j < rows[i].size(); ++j)
            r.push_back(static_cast<int>(std::llround(rows[i][j])));
        out.push_back(r);
    }
    return out;
}

/** Delay + PS queue + FCFS queue, 2 classes: one shared server, one buffered. */
qn::Network<double> started_model() {
    qn::Network<double> m("started");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, c1, D::exp_rate(1.0 / 0.7));
    m.set_service(d, c2, D::exp_rate(1.0 / 1.3));
    m.set_service(q1, c1, D::exp_rate(1.0 / 1.5));
    m.set_service(q1, c2, D::exp_rate(1.0 / 0.8));
    m.set_service(q2, c1, D::exp_rate(1.0 / 0.9));
    m.set_service(q2, c2, D::exp_rate(1.0 / 1.1));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, d, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace


TEST_CASE("getProbSysMarg: the dense model normalizes and its mean is the CTMC one") {
    qn::Network<double> m = dense_model();
    double mass = 0.0;
    std::vector<double> mean;
    law_moments(m.get_struct(), 3, 3, &mass, &mean);
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));

    qn::Network<double> mc = dense_model();
    const std::vector<double> exact = ctmc_totals(mc);
    for (std::size_t i = 0; i < 3; ++i) CHECK(mean[i] == doctest::Approx(exact[i]).epsilon(1e-9));
}

TEST_CASE("getProbSysMarg: a structural zero disturbs neither the mass nor the mean") {
    qn::Network<double> m = zero_demand_model();
    double mass = 0.0;
    std::vector<double> mean;
    law_moments(m.get_struct(), 3, 3, &mass, &mean);
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));

    qn::Network<double> mc = zero_demand_model();
    const std::vector<double> exact = ctmc_totals(mc);
    for (std::size_t i = 0; i < 3; ++i) CHECK(mean[i] == doctest::Approx(exact[i]).epsilon(1e-9));
}

TEST_CASE("getProbSysMarg: each infinite server contributes its own 1/n_j!") {
    qn::Network<double> m = two_delay_model();
    double mass = 0.0;
    std::vector<double> mean;
    law_moments(m.get_struct(), 3, 3, &mass, &mean);
    // Dividing ONCE, as a single aggregated delay row would, leaves this short of 1.
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));

    qn::Network<double> mc = two_delay_model();
    const std::vector<double> exact = ctmc_totals(mc);
    for (std::size_t i = 0; i < 3; ++i) CHECK(mean[i] == doctest::Approx(exact[i]).epsilon(1e-9));
}

TEST_CASE("getProbSysMarg: a state whose total differs from sum(N) has probability 0") {
    qn::Network<double> m = dense_model();
    const nc::NcSolverOptions opt;
    const std::vector<int> bad{1, 0, 0};
    CHECK(nc::solver_nc_getprob_sys_marg(m.get_struct(), opt, bad) == 0.0);
}

TEST_CASE("getProbSysMarg: the approximate engines refuse a structural zero") {
    qn::Network<double> m = zero_demand_model();
    const nc::NcSolverOptions opt;
    const std::vector<int> n{1, 1, 1};
    const char* engines[] = {"bethe", "heur", "huberlaw", "adapart"};
    for (std::size_t e = 0; e < 4; ++e) {
        bool refused = false;
        std::string msg;
        try {
            nc::solver_nc_getprob_sys_marg(m.get_struct(), opt, n, engines[e]);
        } catch (const std::exception& ex) {
            refused = true;
            msg = ex.what();
        }
        CHECK(refused);
        CHECK(msg.find("zero") != std::string::npos);
        CHECK(msg.find(engines[e]) != std::string::npos);
    }
    // The exact engine is unaffected on the same state.
    CHECK(nc::solver_nc_getprob_sys_marg(m.get_struct(), opt, n) > 0.0);
}

TEST_CASE("getProbSysMarg: a multiserver station is refused by name") {
    qn::Network<double> m("marg_multiserver");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::FCFS);
    m.set_number_of_servers(q, 2);
    const std::size_t c1 = m.add_closed_class("Class1", 3.0, d);
    m.set_service(d, c1, D::exp_rate(1.0 / 0.7));
    m.set_service(q, c1, D::exp_rate(1.0 / 1.5));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    m.link(P);

    const nc::NcSolverOptions opt;
    const std::vector<int> n{1, 2};
    bool refused = false;
    std::string msg;
    try {
        nc::solver_nc_getprob_sys_marg(m.get_struct(), opt, n);
    } catch (const std::exception& ex) {
        refused = true;
        msg = ex.what();
    }
    CHECK(refused);
    CHECK(msg.find("multiserver") != std::string::npos);
}


TEST_CASE("from_marginal_and_started: pinned against the MATLAB reference") {
    qn::Network<double> m = started_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::size_t> ph(sn.nclasses, 1);

    // A SHARED SERVER holds every job present, so the block carries n and not
    // s. Writing s here is the defect that was fixed in the Python twin.
    CHECK(as_ints(qn::from_marginal_node_and_started(sn, 2, {2, 0}, {1, 0}, ph)) ==
          std::vector<std::vector<int> >{{2, 0}});
    CHECK(as_ints(qn::from_marginal_node_and_started(sn, 2, {1, 1}, {0, 1}, ph)) ==
          std::vector<std::vector<int> >{{1, 1}});

    // An ORDERED buffer carries the waiting class tags, the server block the
    // started counts in phase one: [buffer | srv(C1) | srv(C2)].
    CHECK(as_ints(qn::from_marginal_node_and_started(sn, 3, {2, 0}, {1, 0}, ph)) ==
          std::vector<std::vector<int> >{{1, 1, 0}});
    CHECK(as_ints(qn::from_marginal_node_and_started(sn, 3, {1, 1}, {0, 1}, ph)) ==
          std::vector<std::vector<int> >{{1, 0, 1}});
    // The empty station keeps the decoder's width: one buffer column here.
    CHECK(as_ints(qn::from_marginal_node_and_started(sn, 3, {0, 0}, {0, 0}, ph)) ==
          std::vector<std::vector<int> >{{0, 0, 0}});

    // More started than present is no state at all, not an error.
    CHECK(qn::from_marginal_node_and_started(sn, 3, {1, 0}, {1, 1}, ph).empty());
}

TEST_CASE("from_marg_node_started: the union over both totals, pinned against MATLAB") {
    qn::Network<double> m = started_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::size_t> ph(sn.nclasses, 1);

    // PS: the started total is immaterial, so this is from_marg_node's answer.
    // Class 2 has population 1, so the split [0,2] is not enumerated at all.
    CHECK(as_ints(qn::from_marg_node_started(sn, 2, 2, 1, ph)) ==
          std::vector<std::vector<int> >{{2, 0}, {1, 1}});

    // FCFS: three (n,s) pairs survive, and the rows come out sorted descending
    // as the reference's trailing unique/flip leaves them.
    CHECK(as_ints(qn::from_marg_node_started(sn, 3, 2, 1, ph)) ==
          std::vector<std::vector<int> >{{2, 1, 0}, {1, 1, 0}, {1, 0, 1}});

    // Both totals zero is the one empty state, at the discipline's own width.
    CHECK(as_ints(qn::from_marg_node_started(sn, 3, 0, 0, ph)) ==
          std::vector<std::vector<int> >{{0, 0, 0}});

    // stot > ntot cannot be realized.
    CHECK(qn::from_marg_node_started(sn, 3, 1, 2, ph).empty());
}

TEST_CASE("from_marg_node: the class splits of a total, capped by classcap") {
    qn::Network<double> m = dense_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::size_t> phases(sn.nclasses, 1);

    // Queue1 is node 2. Class 2 has population 1, so [0,2] is not a state it
    // can hold and the split enumeration must drop it up front.
    const std::vector<std::vector<double> > s2 = qn::from_marg_node(sn, 2, 2, phases);
    CHECK(s2.size() == 2);
    for (std::size_t r = 0; r < s2.size(); ++r) {
        double tot = 0.0;
        for (std::size_t c = 0; c < s2[r].size(); ++c) tot += s2[r][c];
        CHECK(tot == doctest::Approx(2.0));
        CHECK(s2[r][1] <= 1.0);
    }

    // ntot == 0 is the single empty state, at the discipline's own width.
    const std::vector<std::vector<double> > s0 = qn::from_marg_node(sn, 2, 0, phases);
    CHECK(s0.size() == 1);
    for (std::size_t c = 0; c < s0[0].size(); ++c) CHECK(s0[0][c] == doctest::Approx(0.0));
}
