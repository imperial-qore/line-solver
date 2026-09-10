/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Routing reducibility and its repair, against the MATLAB reference
 * `@MNetwork/isRoutingErgodic`, `getReducibilityInfo`, `getAbsorbingStations`
 * and `makeErgodic`.
 *
 * Fixture: D -> Q1 -> Q2 -> Q2. Q2 absorbs, D and Q1 are transient, and the
 * routing graph has three strongly connected components. MATLAB reports
 * numSCCs 3, absorbing {Q2}, transient {Q1, D}, one suggested fix, and
 * makeErgodic returns [0 1 0; 0 0 1; 1 0 0].
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_routing_ergodic.h"
#include "line/lang/qn/network_struct.h"

using line::Matrix;
using line::qn::NetworkStruct;

namespace {

/** D (station 1), Q1 (station 2), Q2 (station 3), one closed class. */
NetworkStruct<double> reducible() {
    NetworkStruct<double> sn;
    sn.nstations = 3;
    sn.nclasses = 1;
    sn.nchains = 1;
    sn.nodes.resize(3);
    const char* names[3] = {"D", "Q1", "Q2"};
    for (std::size_t i = 0; i < 3; ++i) {
        sn.nodes[i].name = names[i];
        sn.nodes[i].station = i + 1;
        sn.nodes[i].nodetype =
            (i == 0) ? line::lang::NodeType::Delay : line::lang::NodeType::Queue;
    }
    Matrix<double> M(3, 3, 0.0);
    M(0, 1) = 1.0;
    M(1, 2) = 1.0;
    M(2, 2) = 1.0;
    sn.P[std::make_pair(std::size_t(1), std::size_t(1))] = M;
    return sn;
}

/** D <-> Q, a two-station cycle. */
NetworkStruct<double> ergodic() {
    NetworkStruct<double> sn;
    sn.nstations = 2;
    sn.nclasses = 1;
    sn.nchains = 1;
    sn.nodes.resize(2);
    sn.nodes[0].name = "D";
    sn.nodes[0].station = 1;
    sn.nodes[0].nodetype = line::lang::NodeType::Delay;
    sn.nodes[1].name = "Q";
    sn.nodes[1].station = 2;
    sn.nodes[1].nodetype = line::lang::NodeType::Queue;
    Matrix<double> M(2, 2, 0.0);
    M(0, 1) = 1.0;
    M(1, 0) = 1.0;
    sn.P[std::make_pair(std::size_t(1), std::size_t(1))] = M;
    return sn;
}

bool contains(const std::vector<std::string>& v, const std::string& s) {
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] == s) return true;
    return false;
}

}  // namespace

TEST_CASE("sn_routing_ergodic: reducible structure matches MATLAB") {
    const NetworkStruct<double> sn = reducible();
    const line::sn::RoutingErgodicityInfo info = line::sn::sn_is_routing_ergodic(sn);
    CHECK(info.isRoutingErgodic == false);
    CHECK(info.isReducible == true);
    CHECK(info.numSCCs == 3u);
    REQUIRE(info.absorbingStations.size() == 1u);
    CHECK(info.absorbingStations[0] == "Q2");
    CHECK(contains(info.transientStations, "D"));
    CHECK(contains(info.transientStations, "Q1"));
}

TEST_CASE("sn_routing_ergodic: one suggested fix per absorbing station") {
    const line::sn::RoutingErgodicityInfo info = line::sn::sn_reducibility_info(reducible());
    REQUIRE(info.suggestedFixes.size() == 1u);
    CHECK(info.suggestedFixes[0] ==
          "Route jobs from Q2 back to D (e.g., P{class}(Q2, D) = 1.0)");
}

TEST_CASE("sn_routing_ergodic: absorbing stations by 1-based node index") {
    const std::vector<std::size_t> idx = line::sn::sn_absorbing_stations(reducible());
    REQUIRE(idx.size() == 1u);
    CHECK(idx[0] == 3u);
}

TEST_CASE("sn_routing_ergodic: make_ergodic redirects the absorbing row") {
    const NetworkStruct<double> sn = reducible();
    const std::map<std::pair<std::size_t, std::size_t>, Matrix<double>> P =
        line::sn::sn_make_ergodic(sn);
    const Matrix<double>& M = P.at(std::make_pair(std::size_t(1), std::size_t(1)));
    const double want[3][3] = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(M(i, j) == doctest::Approx(want[i][j]));
    // the original is untouched: the repair is returned, never applied
    CHECK(sn.P.at(std::make_pair(std::size_t(1), std::size_t(1)))(2, 2) == doctest::Approx(1.0));
}

TEST_CASE("sn_routing_ergodic: an ergodic routing is left alone") {
    const NetworkStruct<double> sn = ergodic();
    const line::sn::RoutingErgodicityInfo info = line::sn::sn_reducibility_info(sn);
    CHECK(info.isRoutingErgodic == true);
    CHECK(info.isReducible == false);
    CHECK(info.absorbingStations.empty());
    CHECK(info.suggestedFixes.empty());
    CHECK(line::sn::sn_absorbing_stations(sn).empty());
}
