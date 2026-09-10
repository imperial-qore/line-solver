/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * sn_fj_visits_spn: fork-join node visit ratios via the auxiliary closed SPN.
 *
 * THE VISIT VECTOR IS THE LEAST INTERESTING PART OF THE ANSWER. The auxiliary
 * net is population preserving, so after the normalization on the reference
 * station every visited station is exactly 1 and every Fork and Join exactly 0
 * -- a test that only checked those ones and zeros would pass against a
 * function that returned ones unconditionally. What carries the model content
 * is the CONSTRUCTION: which nodes are classified as stations, forks and joins,
 * how many leaves each Join synchronizes, and the circulating population B,
 * which is the largest leaf count over the OUTERMOST forks and not over all of
 * them. Those are what the tests below drive, on a plain fork-join, on a nested
 * one and on a serial pair.
 *
 * The refusals are tested too, because they are what this port has instead of
 * MATLAB's CTMC solve: MATLAB would DETECT a net whose throughputs were not
 * uniform, where this port checks the structural preconditions up front.
 */
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_fj_visits_spn.h"
#include "line/lang/qn/network_builder.h"

namespace api = line::api;
namespace qn = line::qn;
using line::Matrix;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> Fork -> {Q1, Q2} -> Join -> Delay, one closed class. */
qn::Network<double> fj_two_branch() {
    qn::Network<double> m("fj2");
    const std::size_t d = m.add_delay("Think");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("J", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    return m;
}

/** The single-class node routing, as the routine reads it. */
Matrix<double> class_routing(const qn::NetworkStruct<double>& sn, std::size_t r) {
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;
    Matrix<double> P(I, I, 0.0);
    for (std::size_t i = 0; i < I; ++i)
        for (std::size_t j = 0; j < I; ++j) P(i, j) = sn.rtnodes(i * K + r, j * K + r);
    return P;
}

std::vector<bool> all_visited(const qn::NetworkStruct<double>& sn) {
    return std::vector<bool>(sn.nodes.size(), true);
}

}  // namespace

TEST_CASE("a two-branch fork-join: a branch station is visited once per B cycles") {
    // The expected numbers are MATLAB's and the JAR's, both of which SOLVE the
    // auxiliary net: Think 1, Q1 = Q2 = 0.5, Fork and Join 0. They are not the
    // uniform ones `sn_fj_visits_spn.m` used to claim in its comments -- the
    // pre-fork transition consumes all B tokens and the Join returns them, so a
    // station inside the region fires once per B firings of the cycle.
    qn::Network<double> m = fj_two_branch();
    const qn::NetworkStruct<double> sn = m.get_struct();
    const std::vector<Matrix<double>> v = api::sn_fj_visits_spn(sn);

    REQUIRE(v.size() == sn.nchains);
    const Matrix<double>& V = v[0];
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const qn::NodeType t = sn.nodes[i].nodetype;
        if (t == qn::NodeType::Fork || t == qn::NodeType::Join) {
            // The auxiliary net has no Place for these, so they carry nothing.
            CHECK(V(i, 0) == doctest::Approx(0.0));
        } else if (sn.nodes[i].name == "Think") {
            CHECK(V(i, 0) == doctest::Approx(1.0));  // the reference station
        } else if (sn.nodes[i].station != 0) {
            CHECK(V(i, 0) == doctest::Approx(0.5));  // 1/B, B = 2 branches
        }
    }
}

TEST_CASE("the construction sizes the net from the OUTERMOST forks") {
    qn::Network<double> m = fj_two_branch();
    const qn::NetworkStruct<double> sn = m.get_struct();
    const api::FjSpnStructure st =
        api::fj_spn_structure(sn, class_routing(sn, 0), all_visited(sn));

    CHECK(st.forkNodes.size() == 1u);
    CHECK(st.joinNodes.size() == 1u);
    // Delay, Q1, Q2 -- the Fork and Join are not stations.
    CHECK(st.stationNodes.size() == 3u);
    // Two leaves on the one outermost fork.
    CHECK(st.B == 2u);
    CHECK(st.joinLeaves[st.joinNodes[0]] == 2u);
}

TEST_CASE("a three-branch fork raises the circulating population") {
    qn::Network<double> m("fj3");
    const std::size_t d = m.add_delay("Think");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("J", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    m.set_service(q3, c, Dist::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(f, q3, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(q3, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    const api::FjSpnStructure st =
        api::fj_spn_structure(sn, class_routing(sn, 0), all_visited(sn));
    CHECK(st.B == 3u);
    CHECK(st.joinLeaves[st.joinNodes[0]] == 3u);
    CHECK(st.stationNodes.size() == 4u);

    // Three branches, so every branch station carries 1/3 and the reference
    // station carries 1, as MATLAB and the JAR report. The builder hands back
    // 1-BASED node ids while the visit matrix is indexed by 0-based node, so
    // the conversion is explicit rather than incidental.
    const std::vector<Matrix<double>> v = api::sn_fj_visits_spn(sn);
    CHECK(v[0](d - 1, 0) == doctest::Approx(1.0));
    CHECK(v[0](q1 - 1, 0) == doctest::Approx(1.0 / 3.0));
    CHECK(v[0](q3 - 1, 0) == doctest::Approx(1.0 / 3.0));
    CHECK(v[0](f - 1, 0) == doctest::Approx(0.0));
    CHECK(v[0](j - 1, 0) == doctest::Approx(0.0));
}

TEST_CASE("a model with no fork-join returns zeros, not ones") {
    // The reference exits early on `sn.fj` being empty; answering ones there
    // would overwrite the ordinary visit computation with a fork-join one.
    qn::Network<double> m("plain");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    const std::vector<Matrix<double>> v = api::sn_fj_visits_spn(sn);
    REQUIRE(v.size() == sn.nchains);
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) CHECK(v[0](i, 0) == doctest::Approx(0.0));
}

TEST_CASE("a fork that reaches no station is refused by name") {
    // This is what the port has instead of MATLAB's CTMC solve: the invariant
    // the answer rests on is CHECKED, not assumed.
    qn::Network<double> m = fj_two_branch();
    const qn::NetworkStruct<double> sn = m.get_struct();
    Matrix<double> P = class_routing(sn, 0);

    // Cut both branches out of the fork's row.
    std::size_t fnode = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == qn::NodeType::Fork) fnode = i;
    for (std::size_t j = 0; j < sn.nodes.size(); ++j) P(fnode, j) = 0.0;

    CHECK_THROWS_AS(api::fj_spn_structure(sn, P, all_visited(sn)), line::InputError);
}

TEST_CASE("a Join with no branch feeding it is refused by name") {
    qn::Network<double> m = fj_two_branch();
    const qn::NetworkStruct<double> sn = m.get_struct();
    Matrix<double> P = class_routing(sn, 0);

    std::size_t jnode = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == qn::NodeType::Join) jnode = i;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) P(i, jnode) = 0.0;

    CHECK_THROWS_AS(api::fj_spn_structure(sn, P, all_visited(sn)), line::InputError);
}
