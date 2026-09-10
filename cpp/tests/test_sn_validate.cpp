/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * sn_validate and the two predicates that had no C++ twin until 2026-08-01.
 *
 * The oracle is the reference's own contract: a consistent struct yields the
 * empty list, and each defect injected in isolation yields exactly one message
 * naming the field it belongs to. Checking the COUNT matters -- the reference
 * accumulates rather than throwing, so a check that silently never fires and
 * one that fires correctly are indistinguishable from a pass/fail alone.
 */
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_predicates.h"
#include "line/api/sn/sn_validate.h"
#include "line/lang/qn/network_builder.h"

namespace qn = line::qn;
namespace api = line::api;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> PS Queue -> Delay, one closed class. */
qn::Network<double> closed_pf() {
    qn::Network<double> m("pf");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** How many messages mention `needle`. */
std::size_t count_with(const std::vector<std::string>& msgs, const std::string& needle) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < msgs.size(); ++i)
        if (msgs[i].find(needle) != std::string::npos) ++n;
    return n;
}

}  // namespace

TEST_CASE("sn_validate accepts a well-formed struct") {
    qn::Network<double> m = closed_pf();
    CHECK(api::sn_validate(m.get_struct()).empty());
    CHECK(api::sn_validate(m.get_struct(), api::ValidationLevel::Minimal).empty());
    CHECK(api::sn_validate(m.get_struct(), api::ValidationLevel::None).empty());
}

TEST_CASE("sn_validate reports a negative rate and a negative SCV") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.rates(1, 0) = -2.0;
    sn.scv(0, 0) = -1.0;

    const std::vector<std::string> errs = api::sn_validate(sn);
    CHECK(count_with(errs, "rates[1,0]") == 1);
    CHECK(count_with(errs, "scv[0,0]") == 1);
    CHECK(count_with(errs, "is negative") == 2);

    // Minimal only checks shapes, so neither value defect is reported.
    CHECK(api::sn_validate(sn, api::ValidationLevel::Minimal).empty());
}

TEST_CASE("sn_validate skips a disabled station-class pair") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.rates(1, 0) = -2.0;
    sn.disabled[1][0] = true;
    CHECK(api::sn_validate(sn).empty());
}

TEST_CASE("sn_validate reports a routing row that does not sum to one") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    // Halve a row: it stays non-empty, so the sum test must fire.
    for (std::size_t j = 0; j < sn.rt.cols(); ++j) sn.rt(0, j) *= 0.5;
    const std::vector<std::string> errs = api::sn_validate(sn);
    CHECK(count_with(errs, "rt row 0 sum") == 1);
}

TEST_CASE("sn_validate reports a non-positive server count but spares a Source") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.stations[1].nservers = 0.0;
    CHECK(count_with(api::sn_validate(sn), "nservers[1]") == 1);

    // A Source carries nservers 0 legitimately; the reference exempts it.
    sn.stations[1].nodetype = qn::NodeType::Source;
    CHECK(count_with(api::sn_validate(sn), "nservers[1]") == 0);
}

TEST_CASE("sn_validate reports a negative closed population") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.classes[0].population = -1.0;
    CHECK(count_with(api::sn_validate(sn), "njobs for class 0") == 1);
}

TEST_CASE("sn_validate reports a rates matrix of the wrong shape") {
    qn::Network<double> m = closed_pf();
    qn::NetworkStruct<double> sn = m.get_struct();
    sn.rates = line::Matrix<double>(3, 3, 0.0);
    const std::vector<std::string> errs = api::sn_validate(sn, api::ValidationLevel::Minimal);
    CHECK(count_with(errs, "rates matrix dimensions") == 1);
}

TEST_CASE("sn_validate index guards return a message only out of bounds") {
    qn::Network<double> m = closed_pf();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(api::sn_validate_station_index(sn, 0).empty());
    CHECK(!api::sn_validate_station_index(sn, 2).empty());
    CHECK(api::sn_validate_class_index(sn, 0).empty());
    CHECK(!api::sn_validate_class_index(sn, 1).empty());
    CHECK(api::sn_validate_node_index(sn, 0).empty());
    CHECK(!api::sn_validate_node_index(sn, 99).empty());
    CHECK(api::sn_validate_node_type(sn, 1, qn::NodeType::Queue).empty());
    CHECK(api::sn_validate_node_type(sn, 1, qn::NodeType::Delay).find("expected Delay") !=
          std::string::npos);
}

TEST_CASE("sn_has_polling and sn_has_srpt read the station disciplines") {
    qn::Network<double> m = closed_pf();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(api::sn_has_polling(sn) == false);
    CHECK(api::sn_has_srpt(sn) == false);

    qn::NetworkStruct<double> s2 = sn;
    s2.stations[1].sched = SchedStrategy::POLLING;
    CHECK(api::sn_has_polling(s2) == true);
    s2.stations[1].sched = SchedStrategy::SRPT;
    CHECK(api::sn_has_srpt(s2) == true);
    CHECK(api::sn_has_polling(s2) == false);
}
