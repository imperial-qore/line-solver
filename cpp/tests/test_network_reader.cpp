/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The model.json front end (Stage 5): a Network model reaches SolverMVA through
 * the same builder + finalize the programmatic API uses. The reference numbers
 * are MATLAB's SolverMVA getAvgTable on the same model, exported to JSON by
 * `linemodel_save` -- so this asserts the reader reconstructs the model the
 * reference solved, not merely that it parses.
 */

#include <string>

#include "doctest.h"
#include "json.hpp"
#include "line/io/network_reader.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;

namespace {

using json = nlohmann::json;

mva::AvgResult<double> solve(const std::string& text) {
    const json root = json::parse(text);
    qn::Network<double> net = io::build_network_from_json<double>(root);
    mva::MvaOptions opt;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
}

// Station row lookup by name is what the AvgTable comparison keys on.
std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

TEST_CASE("open M/M/1 with a class switch matches the MATLAB getAvgTable") {
    // The cs_implicit example (matlab/examples/basic/classSwitching), as its
    // linemodel_save JSON: Source -> ClassSwitch -> Queue(FCFS) -> Sink, two
    // open classes, the switch matrix [[0.3,0.7],[1,0]]. Distributions arrive as
    // fit blocks (fitMean), the form the MATLAB writer emits.
    const std::string text = R"({
      "format":"line-model","version":"1.0","model":{
       "type":"Network","name":"mm1cs",
       "nodes":[
        {"name":"Source 1","type":"Source","service":{
           "Class1":{"type":"Exp","fit":{"method":"fitMean","mean":10.0}},
           "Class2":{"type":"Exp","fit":{"method":"fitMean","mean":2.0}}}},
        {"name":"Queue 1","type":"Queue","scheduling":"FCFS","service":{
           "Class1":{"type":"Exp","fit":{"method":"fitMean","mean":1.0}},
           "Class2":{"type":"Exp","fit":{"method":"fitMean","mean":1.0}}}},
        {"name":"Sink 1","type":"Sink"},
        {"name":"ClassSwitch 1","type":"ClassSwitch","classSwitchMatrix":{
           "Class1":{"Class1":0.3,"Class2":0.7},"Class2":{"Class1":1.0,"Class2":0.0}}}
       ],
       "classes":[{"name":"Class1","type":"Open"},{"name":"Class2","type":"Open"}],
       "routing":{"type":"matrix","matrix":{
         "Class1,Class1":{"Source 1":{"ClassSwitch 1":1.0},"ClassSwitch 1":{"Queue 1":1.0},"Queue 1":{"Sink 1":1.0}},
         "Class2,Class2":{"Source 1":{"ClassSwitch 1":1.0},"ClassSwitch 1":{"Queue 1":1.0},"Queue 1":{"Sink 1":1.0}}}}
      }})";
    const mva::AvgResult<double> r = solve(text);
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const std::size_t q = station_of(sn, "Queue 1");
    // class order is Class1, Class2 (declaration order)
    CHECK(r.QN(q, 0) == doctest::Approx(1.325).epsilon(1e-6));
    CHECK(r.QN(q, 1) == doctest::Approx(0.175).epsilon(1e-6));
    CHECK(r.UN(q, 0) == doctest::Approx(0.53).epsilon(1e-6));
    CHECK(r.UN(q, 1) == doctest::Approx(0.07).epsilon(1e-6));
    CHECK(r.RN(q, 0) == doctest::Approx(2.5).epsilon(1e-6));
    CHECK(r.WN(q, 0) == doctest::Approx(2.20833333).epsilon(1e-6));
    CHECK(r.WN(q, 1) == doctest::Approx(0.29166667).epsilon(1e-6));
    CHECK(r.TN(q, 0) == doctest::Approx(0.53).epsilon(1e-6));
    CHECK(r.TN(q, 1) == doctest::Approx(0.07).epsilon(1e-6));
}

TEST_CASE("closed two-class network solves through the JSON path") {
    // Delay + PS Queue, populations 2 and 1; the explicit-params Exp form.
    const std::string text = R"({
      "model":{"type":"Network","name":"cqn",
       "nodes":[
        {"name":"Delay","type":"Delay","service":{
           "C1":{"type":"Exp","params":{"lambda":1.0}},"C2":{"type":"Exp","params":{"lambda":1.0}}}},
        {"name":"Q","type":"Queue","scheduling":"PS","service":{
           "C1":{"type":"Exp","params":{"lambda":2.0}},"C2":{"type":"Exp","params":{"lambda":2.0}}}}
       ],
       "classes":[
        {"name":"C1","type":"Closed","population":2,"refNode":"Delay"},
        {"name":"C2","type":"Closed","population":1,"refNode":"Delay"}],
       "routing":{"type":"matrix","matrix":{
         "C1,C1":{"Delay":{"Q":1.0},"Q":{"Delay":1.0}},
         "C2,C2":{"Delay":{"Q":1.0},"Q":{"Delay":1.0}}}}}})";
    const mva::AvgResult<double> r = solve(text);
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const std::size_t q = station_of(sn, "Q");
    // Closed product-form reference (exact MVA): utilisations sum below one and
    // the throughputs balance the delay visits.
    CHECK(r.QN(q, 0) == doctest::Approx(0.94736842).epsilon(1e-6));
    CHECK(r.QN(q, 1) == doctest::Approx(0.47368421).epsilon(1e-6));
    CHECK(r.TN(q, 0) == doctest::Approx(1.05263158).epsilon(1e-6));
}

TEST_CASE("a non-Network model is refused by name") {
    const std::string lqn = R"({"model":{"type":"LayeredNetwork","name":"x"}})";
    CHECK_THROWS_AS(io::build_network_from_json<double>(json::parse(lqn)), UnsupportedError);
}

TEST_CASE("an unsupported node type is refused by name") {
    // The exemplar was "Place" until Place and Transition became readable, then
    // "Logger" until the Logger did. EVERY node class the reference declares
    // (matlab/src/lang/nodes) now has a branch, so the exemplar can no longer be
    // a real one and the case pins what is left: a type this reader does not
    // recognise is NAMED, not skipped past into a model missing a node.
    const std::string text = R"({"model":{"type":"Network","name":"x",
      "nodes":[{"name":"P","type":"Junction"}],
      "classes":[{"name":"C","type":"Open"}],
      "routing":{"type":"matrix","matrix":{}}}})";
    CHECK_THROWS_AS(io::build_network_from_json<double>(json::parse(text)), UnsupportedError);
}

TEST_CASE("an unreconstructable fitMeanAndSCV family is refused, not degraded") {
    // A Gamma asked for by moments has no exact reconstruction here; the reader
    // must refuse rather than silently fit an exponential of the wrong SCV.
    const std::string text = R"({"model":{"type":"Network","name":"x",
      "nodes":[
       {"name":"S","type":"Source","service":{"C":{"type":"Gamma","fit":{"method":"fitMeanAndSCV","mean":1.0,"scv":4.0}}}},
       {"name":"K","type":"Queue","scheduling":"FCFS","service":{"C":{"type":"Exp","params":{"lambda":2.0}}}},
       {"name":"T","type":"Sink"}],
      "classes":[{"name":"C","type":"Open"}],
      "routing":{"type":"matrix","matrix":{"C,C":{"S":{"K":1.0},"K":{"T":1.0}}}}}})";
    CHECK_THROWS_AS(io::build_network_from_json<double>(json::parse(text)), UnsupportedError);
}

TEST_CASE("an embedded Cache node (cacheqn) solves through the JSON path") {
    // gallery_cache_routing as model.json: Source -> Cache(LRU,4,cap2) with
    // hit/miss routed to distinct FCFS queues. Exercises the Cache reader
    // (numItems, itemLevelCap, replacementStrategy, popularity, hit/missClass).
    const std::string text = R"({"model":{"type":"Network","name":"Cache-Routing",
      "nodes":[
       {"name":"Source","type":"Source","service":{"InitClass":{"type":"Exp","params":{"lambda":1.0}}}},
       {"name":"Cache","type":"Cache","numItems":4,"itemLevelCap":[2],"replacementStrategy":"LRU",
        "hitClass":{"InitClass":"HitClass"},"missClass":{"InitClass":"MissClass"},
        "popularity":{"InitClass":{"type":"DiscreteSampler","params":{"p":[0.25,0.25,0.25,0.25],"x":[1,2,3,4]}}}},
       {"name":"HitQueue","type":"Queue","scheduling":"FCFS","service":{"HitClass":{"type":"Exp","params":{"lambda":2.0}}}},
       {"name":"MissQueue","type":"Queue","scheduling":"FCFS","service":{"MissClass":{"type":"Exp","params":{"lambda":1.0}}}},
       {"name":"Sink","type":"Sink"}],
      "classes":[{"name":"InitClass","type":"Open"},{"name":"HitClass","type":"Open"},{"name":"MissClass","type":"Open"}],
      "routing":{"type":"matrix","matrix":{
        "InitClass,InitClass":{"Source":{"Cache":1.0}},
        "HitClass,HitClass":{"Cache":{"HitQueue":1.0},"HitQueue":{"Sink":1.0}},
        "MissClass,MissClass":{"Cache":{"MissQueue":1.0},"MissQueue":{"Sink":1.0}}}}}})";
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
    const qn::NetworkStruct<double>& sn = net.get_struct();
    CHECK(r.QN(station_of(sn, "HitQueue"), 1) == doctest::Approx(0.14286).epsilon(1e-4));
    CHECK(r.QN(station_of(sn, "MissQueue"), 2) == doctest::Approx(0.33333).epsilon(1e-4));
    CHECK(r.actualmethod == "default");
}

// A closed two-station model carrying all three rate-scaling blocks. The
// classDependence table is beta(n) = min(n1, 2) broadcast to both classes over a
// 2x2 lattice, the form `cd_scaling_table` in linemodel_save.m emits.
const std::string kDepModel = R"({"model":{"type":"Network","name":"dep",
  "nodes":[
   {"name":"Delay","type":"Delay","service":{
      "C1":{"type":"Exp","params":{"lambda":1.0}},"C2":{"type":"Exp","params":{"lambda":1.0}}}},
   {"name":"Q","type":"Queue","scheduling":"PS","servers":1,
    "service":{"C1":{"type":"Exp","params":{"lambda":2.0}},"C2":{"type":"Exp","params":{"lambda":2.0}}},
    "loadDependence":{"type":"loadDependent","scaling":[1.0,2.0]},
    "classDependence":{"type":"classDependent","cutoffs":[2,2],"peak":[2.0,2.0],
      "scaling":{"0,0":[1.0,1.0],"1,0":[1.0,1.0],"2,0":[2.0,2.0],
                 "0,1":[1.0,1.0],"1,1":[1.0,1.0],"2,1":[2.0,2.0],
                 "0,2":[1.0,1.0],"1,2":[1.0,1.0],"2,2":[2.0,2.0]}},
    "jointDependence":{"type":"jointDependent","cutoffs":[2,2],"peak":[3.0,3.0],
      "scaling":{"0,0":[1.0,1.0],"1,0":[1.5,1.5],"2,0":[3.0,3.0],
                 "0,1":[1.0,1.0],"1,1":[1.5,1.5],"2,1":[3.0,3.0],
                 "0,2":[1.0,1.0],"1,2":[1.5,1.5],"2,2":[3.0,3.0]}}}],
  "classes":[{"name":"C1","type":"Closed","population":1,"refNode":"Delay"},
             {"name":"C2","type":"Closed","population":1,"refNode":"Delay"}],
  "routing":{"type":"matrix","matrix":{
    "C1,C1":{"Delay":{"Q":1.0},"Q":{"Delay":1.0}},
    "C2,C2":{"Delay":{"Q":1.0},"Q":{"Delay":1.0}}}}}})";

TEST_CASE("loadDependence, classDependence and jointDependence survive the JSON path") {
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(kDepModel));
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const std::size_t q = station_of(sn, "Q");
    const qn::Station<double>& st = sn.stations[q];

    // lldscaling is indexed FROM POPULATION ONE, so entry 0 is the one-job rate
    REQUIRE(st.lldscaling.size() == 2);
    CHECK(st.lldscaling[0] == doctest::Approx(1.0));
    CHECK(st.lldscaling[1] == doctest::Approx(2.0));

    REQUIRE(static_cast<bool>(st.cdscaling));
    REQUIRE(static_cast<bool>(st.jdscaling));
    std::vector<double> n(2, 0.0);
    n[0] = 1;
    CHECK(st.cdscaling(n)[0] == doctest::Approx(1.0));
    CHECK(st.jdscaling(n)[0] == doctest::Approx(1.5));
    n[0] = 2;
    CHECK(st.cdscaling(n)[0] == doctest::Approx(2.0));
    CHECK(st.jdscaling(n)[1] == doctest::Approx(3.0));
    // Beyond the tabulated cutoff the scaling SATURATES rather than falling back
    // to neutral, which is what clamping the population to `cutoffs` buys.
    n[0] = 7;
    CHECK(st.cdscaling(n)[0] == doctest::Approx(2.0));
    CHECK(st.jdscaling(n)[0] == doctest::Approx(3.0));

    REQUIRE(st.cdscalingpeak.size() == 2);
    CHECK(st.cdscalingpeak[0] == doctest::Approx(2.0));
    REQUIRE(st.jdscalingpeak.size() == 2);
    CHECK(st.jdscalingpeak[1] == doctest::Approx(3.0));
}

TEST_CASE("a legacy dependence block with no declared peak takes the table maximum") {
    // `cd_peak_scaling` skips the all-zero composition and every non-finite
    // entry, so the peak here is 2 (the Inf is not the normalizer).
    std::string text = kDepModel;
    const std::string from = "\"peak\":[2.0,2.0],\n      \"scaling\":{\"0,0\":[1.0,1.0]";
    const std::string to = "\"scaling\":{\"0,0\":[9.0,9.0]";
    const std::size_t at = text.find(from);
    REQUIRE(at != std::string::npos);
    text = text.substr(0, at) + to + text.substr(at + from.size());
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const qn::Station<double>& st = sn.stations[station_of(sn, "Q")];
    REQUIRE(st.cdscalingpeak.size() == 2);
    CHECK(st.cdscalingpeak[0] == doctest::Approx(2.0));
    CHECK(st.cdscalingpeak[1] == doctest::Approx(2.0));
}

TEST_CASE("a rate-scaling block outside a Queue or Delay is refused by name") {
    const std::string text = R"({"model":{"type":"Network","name":"bad",
      "nodes":[
       {"name":"Source","type":"Source","service":{"C":{"type":"Exp","params":{"lambda":1.0}}},
        "classDependence":{"type":"classDependent","cutoffs":[1],"peak":[1.0],
          "scaling":{"0":[1.0],"1":[1.0]}}},
       {"name":"Sink","type":"Sink"}],
      "classes":[{"name":"C","type":"Open"}],
      "routing":{"type":"matrix","matrix":{"C,C":{"Source":{"Sink":1.0}}}}}})";
    CHECK_THROWS_AS(io::build_network_from_json<double>(json::parse(text)), UnsupportedError);
}

TEST_CASE("a dependence block declaring a foreign type is refused, not ignored") {
    std::string text = kDepModel;
    const std::string from = "\"type\":\"jointDependent\"";
    const std::size_t at = text.find(from);
    REQUIRE(at != std::string::npos);
    text = text.substr(0, at) + "\"type\":\"marginDependent\"" + text.substr(at + from.size());
    CHECK_THROWS_AS(io::build_network_from_json<double>(json::parse(text)), UnsupportedError);
}

TEST_CASE("the three dependence blocks a writer actually emits are read verbatim") {
    // Copied BYTE FOR BYTE out of `save_model` on a one-class model with
    // beta(n)=min(n,2) peak 2, eta(n)=3 peak 3 and lld [1,2,2,2]. The fixture
    // above is hand-written and 2-class, so it never exercises what a real
    // writer produces at K=1: keys carry NO comma, and `cutoffs` comes from the
    // class population rather than being chosen by the test.
    const std::string text = R"({"model":{"type":"Network","name":"cdjd",
      "nodes":[
       {"name":"Think","type":"Delay","scheduling":"INF",
        "service":{"C1":{"type":"Exp","params":{"lambda":1.0}}}},
       {"name":"Q","type":"Queue","scheduling":"PS",
        "service":{"C1":{"type":"Exp","params":{"lambda":1.5}}},
        "loadDependence":{"type":"loadDependent","scaling":[1.0,2.0,2.0,2.0]},
        "classDependence":{"type":"classDependent","cutoffs":[4],
          "scaling":{"0":[0.0],"1":[1.0],"2":[2.0],"3":[2.0],"4":[2.0]},"peak":[2.0]},
        "jointDependence":{"type":"jointDependent","cutoffs":[4],
          "scaling":{"0":[3.0],"1":[3.0],"2":[3.0],"3":[3.0],"4":[3.0]},"peak":[3.0]}}],
      "classes":[{"name":"C1","type":"Closed","population":4,"refNode":"Think"}],
      "routing":{"type":"matrix","matrix":{
        "C1,C1":{"Think":{"Q":1.0},"Q":{"Think":1.0}}}}}})";
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const qn::Station<double>& st = sn.stations[station_of(sn, "Q")];

    REQUIRE(st.lldscaling.size() == 4);
    CHECK(st.lldscaling[0] == doctest::Approx(1.0));
    CHECK(st.lldscaling[3] == doctest::Approx(2.0));

    REQUIRE(static_cast<bool>(st.cdscaling));
    REQUIRE(static_cast<bool>(st.jdscaling));
    std::vector<double> n(1, 1.0);
    CHECK(st.cdscaling(n)[0] == doctest::Approx(1.0));
    n[0] = 2;
    CHECK(st.cdscaling(n)[0] == doctest::Approx(2.0));
    n[0] = 7;  // clamped to the population-4 cutoff, so still saturated
    CHECK(st.cdscaling(n)[0] == doctest::Approx(2.0));
    CHECK(st.jdscaling(n)[0] == doctest::Approx(3.0));

    REQUIRE(st.cdscalingpeak.size() == 1);
    CHECK(st.cdscalingpeak[0] == doctest::Approx(2.0));
    REQUIRE(st.jdscalingpeak.size() == 1);
    CHECK(st.jdscalingpeak[0] == doctest::Approx(3.0));
}

}  // namespace
