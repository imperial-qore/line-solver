/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The interchange keys that had no round trip: the marking-dependent firing
 * rate of an SPN mode, the memory budget of a finite capacity region, and the
 * class dimension of a transition arc.
 *
 * WHY EACH IS A WRONG NUMBER AND NOT A MISSING FEATURE.
 *
 *  - `firingRateDependence` was REFUSED by the reader, so every MATLAB or
 *    python net whose transition scales its rate with the marking was
 *    unreachable under `lang='cpp'`; and it was dropped by the writer, so a net
 *    built in C++ wrote a document that reloads at the NOMINAL rate. Same net
 *    on the wire, different marking process on either side of it.
 *  - `globalMaxMemory` reaches three consumers (`solver_ctmc_fcr.h`,
 *    `solver_nc_lossn.h`, `jmt_writer.h`) and was never written back, so a
 *    region reloaded with unbounded memory and all three silently dropped a
 *    declared constraint.
 *  - a transition arc written with no `class` key is not merely lossy: both
 *    reference readers dereference it unguarded (`linemodel_load.m:852`,
 *    `linemodel_io.py:2399`), so the document does not load at all.
 */

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "json.hpp"
#include "line/io/network_reader.h"
#include "line/io/network_writer.h"

using namespace line;

namespace {

using json = nlohmann::json;

json roundtrip(const std::string& text) {
    const json root = json::parse(text);
    qn::Network<double> net = io::build_network_from_json<double>(root);
    return io::network_to_json(net.get_struct());
}

const json& node_of(const json& model, const std::string& nm) {
    for (const json& nd : model.at("nodes"))
        if (nd.at("name") == nm) return nd;
    FAIL("no node named ", nm);
    return model;
}

/**
 * A 2-place cyclic net whose first transition scales with the marking of P0.
 *
 * The `firingRateDependence` block is in the reference writers' own spelling:
 * `slots` names the enabling (place, class) pairs the multiplier reads,
 * `cutoffs` saturates each of them, and `scaling` is keyed by the comma-joined
 * 0-based counts in slot order. The tabulated points are deliberately NOT the
 * neutral 1, so a dropped table is a different answer rather than the same one,
 * and the count 1 is deliberately ABSENT so the sparse-table default is
 * exercised beside the two tabulated points.
 */
const char* kFiringDepJson = R"JSON({
  "format": "line-model", "version": "1.0",
  "model": {
    "type": "Network", "name": "spn_firingdep",
    "nodes": [
      {"name": "P0", "type": "Place"},
      {"name": "P1", "type": "Place"},
      {"name": "T0", "type": "Transition",
       "modes": [{"name": "M0", "timingStrategy": "TIMED", "firingPriority": 1,
                  "distribution": {"type": "Exp", "params": {"lambda": 1.0}},
                  "enablingConditions": [{"node": "P0", "class": "Class1", "count": 1}],
                  "firingOutcomes": [{"node": "P1", "class": "Class1", "count": 1}],
                  "firingRateDependence": {
                    "slots": [{"node": "P0", "class": "Class1"}],
                    "cutoffs": [2],
                    "scaling": {"0": 0.0, "2": 1.5}}}]},
      {"name": "T1", "type": "Transition",
       "modes": [{"name": "M1", "timingStrategy": "TIMED", "firingPriority": 1,
                  "distribution": {"type": "Exp", "params": {"lambda": 2.0}},
                  "enablingConditions": [{"node": "P1", "class": "Class1", "count": 1}],
                  "firingOutcomes": [{"node": "P0", "class": "Class1", "count": 1}]}]}
    ],
    "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "P0"}],
    "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
      "P0": {"T0": 1.0}, "T0": {"P1": 1.0}, "P1": {"T1": 1.0}, "T1": {"P0": 1.0}}}}
  }
})JSON";

}  // namespace

TEST_CASE("a marking-dependent firing rate rebuilds into the handle the CTMC calls") {
    const json root = json::parse(kFiringDepJson);
    qn::Network<double> net = io::build_network_from_json<double>(root);
    const qn::NetworkStruct<double>& sn = net.get_struct();

    // T0 is the third declared node, so its 1-based node index is 3.
    const std::map<std::size_t, qn::TransitionParam<double> >::const_iterator tp =
        sn.transparam.find(3);
    REQUIRE(tp != sn.transparam.end());
    REQUIRE(tp->second.firingdep.size() == 1);
    REQUIRE(bool(tp->second.firingdep[0]));

    // The argument is the node-indexed marking summed over classes, which is
    // what `state_events.h` hands the handle; P0 is node 1, so slot 0 is mk[0].
    std::vector<double> mk(sn.nodes.size(), 0.0);
    CHECK(tp->second.firingdep[0](mk) == doctest::Approx(0.0));
    // The count 1 is not in the table: a marking the table does not name is
    // NEUTRAL, not zero, which is the default all three reference handles apply.
    mk[0] = 1.0;
    CHECK(tp->second.firingdep[0](mk) == doctest::Approx(1.0));
    mk[0] = 2.0;
    CHECK(tp->second.firingdep[0](mk) == doctest::Approx(1.5));
    // Above the cutoff the multiplier SATURATES rather than falling to the
    // neutral 1: the count is clamped before the lookup, as the reference
    // handles clamp it.
    mk[0] = 7.0;
    CHECK(tp->second.firingdep[0](mk) == doctest::Approx(1.5));
    // A negative count cannot arise from a marking, but the reference clamps it
    // to 0 rather than failing the lookup, and so does this.
    mk[0] = -3.0;
    CHECK(tp->second.firingdep[0](mk) == doctest::Approx(0.0));

    // The mode that declares no dependence keeps an EMPTY slot, so the vector
    // stays aligned with the mode index; a shifted vector would apply one
    // mode's multiplier to another.
    const std::map<std::size_t, qn::TransitionParam<double> >::const_iterator tp1 =
        sn.transparam.find(4);
    REQUIRE(tp1 != sn.transparam.end());
    REQUIRE(tp1->second.firingdep.size() == 1);
    CHECK_FALSE(bool(tp1->second.firingdep[0]));
}

TEST_CASE("the firing-rate table survives being written back") {
    const json out = roundtrip(kFiringDepJson);
    const json& t0 = node_of(out, "T0");
    REQUIRE(t0.contains("modes"));
    const json& m0 = t0.at("modes")[0];
    REQUIRE(m0.contains("firingRateDependence"));
    const json& frm = m0.at("firingRateDependence");

    REQUIRE(frm.at("slots").size() == 1);
    CHECK(frm.at("slots")[0].at("node") == "P0");
    CHECK(frm.at("slots")[0].at("class") == "Class1");
    // P0 declares no buffer, so it is an open place and the writer tabulates it
    // to the saturation cutoff of 10 the reference writers use -- an unbounded
    // place would otherwise make the lattice infinite.
    CHECK(frm.at("cutoffs")[0] == 10);
    CHECK(frm.at("scaling").at("0").get<double>() == doctest::Approx(0.0));
    // The gap in the input table reappears as the neutral multiplier, because
    // that is what the handle answers there; writing 0 would invent a barrier.
    CHECK(frm.at("scaling").at("1").get<double>() == doctest::Approx(1.0));
    CHECK(frm.at("scaling").at("2").get<double>() == doctest::Approx(1.5));
    // ABOVE the input's cutoff of 2 the handle SATURATES, so the wider lattice
    // the writer tabulates carries the cutoff value and not the neutral one --
    // the saturation the input declared, made explicit rather than lost.
    CHECK(frm.at("scaling").at("7").get<double>() == doctest::Approx(1.5));
    CHECK(frm.at("scaling").size() == 11);

    // Reading the written document back must land the SAME function of the
    // marking, now tabulated to 10 instead of clamped at 2.
    const json again = roundtrip(out.dump());
    const json& frm2 = node_of(again, "T0").at("modes")[0].at("firingRateDependence");
    CHECK(frm2.at("scaling").at("0").get<double>() == doctest::Approx(0.0));
    CHECK(frm2.at("scaling").at("1").get<double>() == doctest::Approx(1.0));
    CHECK(frm2.at("scaling").at("2").get<double>() == doctest::Approx(1.5));
    CHECK(frm2.at("scaling").at("9").get<double>() == doctest::Approx(1.5));
    CHECK(frm2.at("cutoffs")[0] == 10);
    // Idempotent from here on: the second document is a fixed point of the pair.
    CHECK(frm2.at("scaling") == frm.at("scaling"));

    // The mode with no dependence must not acquire one.
    CHECK_FALSE(node_of(out, "T1").at("modes")[0].contains("firingRateDependence"));
}

TEST_CASE("a transition arc is written with the class the reference readers require") {
    const json out = roundtrip(kFiringDepJson);
    const json& m0 = node_of(out, "T0").at("modes")[0];
    REQUIRE(m0.at("enablingConditions").size() == 1);
    CHECK(m0.at("enablingConditions")[0].at("node") == "P0");
    CHECK(m0.at("enablingConditions")[0].at("class") == "Class1");
    CHECK(m0.at("enablingConditions")[0].at("count") == 1.0);
    REQUIRE(m0.at("firingOutcomes").size() == 1);
    CHECK(m0.at("firingOutcomes")[0].at("class") == "Class1");
}

TEST_CASE("the memory budget of a finite capacity region survives the round trip") {
    // `globalMaxMemory` is enforced BESIDE the job cap rather than folded into
    // it, so unlike `classMaxMemory` -- which the reader converts into a job cap
    // through the class size -- it has no other spelling to survive under.
    static const char* kJson = R"JSON({
      "model": {
        "type": "Network", "name": "fcr_mem",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}},
        "finiteCapacityRegions": [
          {"name": "FCR1",
           "stations": [{"node": "Queue1", "classSize": {"Class1": 2.0}}],
           "globalMaxJobs": 4, "globalMaxMemory": 12,
           "dropRule": {"Class1": "drop"}}
        ]
      }
    })JSON";
    const json out = roundtrip(kJson);
    REQUIRE(out.contains("finiteCapacityRegions"));
    const json& fcr = out.at("finiteCapacityRegions")[0];
    CHECK(fcr.at("globalMaxJobs") == 4.0);
    CHECK(fcr.at("globalMaxMemory") == 12.0);
    CHECK(fcr.at("stations")[0].at("classSize").at("Class1") == 2.0);

    // Idempotent: a second pass through the reader must not move the budget.
    const json again = roundtrip(out.dump());
    CHECK(again.at("finiteCapacityRegions")[0].at("globalMaxMemory") == 12.0);
}

/**
 * The four wire node types no interchange test reached: Fork, Join, Router and
 * the `LogTunnel` spelling of a Logger.
 *
 * `node_type_str` in linemodel_save.m emits exactly twelve type names. The other
 * eight are exercised by `test_network_reader.cpp`, `test_network_roundtrip.cpp`
 * and `test_network_reader_spn.cpp`; with these four the reader is covered on
 * all twelve, which is the denominator that matters -- a type with no branch is
 * a model that cannot cross from MATLAB or python at all.
 */
TEST_CASE("fork, join, router and the LogTunnel spelling all cross the wire") {
    static const char* kJson = R"JSON({
      "model": {
        "type": "Network", "name": "forkjoin",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Router1", "type": "Router"},
          {"name": "Fork1", "type": "Fork", "tasksPerLink": 2},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 2.0}}}},
          {"name": "Join1", "type": "Join", "forkNode": "Fork1",
           "joinStrategy": "PARTIAL", "joinQuorum": 1},
          {"name": "Log1", "type": "LogTunnel", "fileName": "trace.csv"},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Router1": 1.0}, "Router1": {"Fork1": 1.0},
          "Fork1": {"Queue1": 1.0}, "Queue1": {"Join1": 1.0},
          "Join1": {"Log1": 1.0}, "Log1": {"Sink": 1.0}}}}
      }
    })JSON";
    const json root = json::parse(kJson);
    qn::Network<double> net = io::build_network_from_json<double>(root);
    const qn::NetworkStruct<double>& sn = net.get_struct();
    // BY NAME, not by position: the reader declares the nodes in dependency
    // order (a Join needs its Fork, a ClassSwitch the class count), so
    // `sn.nodes` is in creation order and a Join declared in the middle of the
    // document lands at the end of it.
    std::map<std::string, lang::NodeType> kind;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) kind[sn.nodes[i].name] = sn.nodes[i].nodetype;
    CHECK(kind["Router1"] == lang::NodeType::Router);
    CHECK(kind["Fork1"] == lang::NodeType::Fork);
    CHECK(kind["Join1"] == lang::NodeType::Join);
    // `LogTunnel` is the second spelling of a Logger, not a node of its own, so
    // it lands on NodeType::Logger and is written back under the one name.
    CHECK(kind["Log1"] == lang::NodeType::Logger);

    const json out = io::network_to_json(sn);
    CHECK(node_of(out, "Fork1").at("tasksPerLink") == 2.0);
    CHECK(node_of(out, "Join1").at("forkNode") == "Fork1");
    // Both are SCALAR on the wire, not per-class maps: linemodel_save.m:501-524
    // scans the classes and writes the last non-default one as a bare value.
    CHECK(node_of(out, "Join1").at("joinStrategy") == "PARTIAL");
    CHECK(node_of(out, "Join1").at("joinQuorum") == 1.0);
    CHECK(node_of(out, "Log1").at("type") == "Logger");
    CHECK(node_of(out, "Log1").at("fileName") == "trace.csv");
    CHECK(node_of(out, "Router1").at("type") == "Router");
}

TEST_CASE("a region with no memory budget does not acquire one") {
    static const char* kJson = R"JSON({
      "model": {
        "type": "Network", "name": "fcr_nomem",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}},
        "finiteCapacityRegions": [
          {"name": "FCR1", "stations": [{"node": "Queue1"}], "globalMaxJobs": 3}
        ]
      }
    })JSON";
    const json out = roundtrip(kJson);
    CHECK_FALSE(out.at("finiteCapacityRegions")[0].contains("globalMaxMemory"));
}
