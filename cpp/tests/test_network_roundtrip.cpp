/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * model.json -> `qn::Network` -> model.json, on the keys the reader gained
 * beyond the product-form core.
 *
 * WHY A ROUND TRIP AND NOT A PARSE TEST. A reader that DROPS a key parses
 * every file it is given and answers for a different model; that is the exact
 * failure `reject_unconsumed_model_keys` exists to catch, and the exact failure
 * a writer paired with it turns into a visible diff. Reading a document and
 * writing it back is the only check that covers both directions at once: a key
 * the reader ignores never reappears, and a key the writer invents was never
 * read.
 *
 * The documents below are in the reference writers' own spelling
 * (`linemodel_save.m`, `save_model`), not a convenient subset of it.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "doctest.h"
#include "json.hpp"
#include "line/io/lqn_json_reader.h"
#include "line/io/network_reader.h"
#include "line/io/network_writer.h"

// A directory the test may write a probe file into. The build defines
// LINE_MP_REPO_ROOT; the probe is removed as soon as it is read back.
#ifndef LINE_MP_SCRATCH
#define LINE_MP_SCRATCH "/tmp"
#endif

using namespace line;

namespace {

using json = nlohmann::json;

/** Read a document and write the model back out. */
json roundtrip(const std::string& text) {
    const json root = json::parse(text);
    qn::Network<double> net = io::build_network_from_json<double>(root);
    return io::network_to_json(net.get_struct());
}

/** The node object of a given name in a written model. */
const json& node_of(const json& model, const std::string& nm) {
    for (const json& nd : model.at("nodes"))
        if (nd.at("name") == nm) return nd;
    FAIL("no node named ", nm);
    return model;
}

TEST_CASE("a finite capacity region survives the round trip") {
    // The block `fcr_mm1kdrop` exports: one region over one station, a global
    // cap of 3 and a drop rule. Before the region was read, this model solved
    // as the UNBOUNDED M/M/1 and reported rho/(1-rho) = 4 against the exact
    // M/M/1/K value.
    static const char* kJson = R"({
      "format": "line-model", "version": "1.0",
      "model": {
        "type": "Network", "name": "fcr",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.8}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}},
        "finiteCapacityRegions": [
          {"name": "FCR1", "stations": [{"node": "Queue1"}],
           "globalMaxJobs": 3,
           "dropRule": {"Class1": "drop"}}
        ]
      }
    })";
    const json out = roundtrip(kJson);
    REQUIRE(out.contains("finiteCapacityRegions"));
    const json& fcr = out.at("finiteCapacityRegions")[0];
    CHECK(fcr.at("globalMaxJobs") == 3.0);
    CHECK(fcr.at("stations")[0].at("node") == "Queue1");
    CHECK(fcr.at("dropRule").at("Class1") == "drop");
}

TEST_CASE("a non-probabilistic dispatcher is not flattened into the matrix") {
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "rrobin",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 2.0}}}}
        ],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2,
                     "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}}}},
        "routingStrategies": {"Delay": {"Class1": "RROBIN"}}
      }
    })";
    const json out = roundtrip(kJson);
    REQUIRE(out.contains("routingStrategies"));
    CHECK(out.at("routingStrategies").at("Delay").at("Class1") == "RROBIN");
}

TEST_CASE("impatience, balking and the retrial orbit stay distinct") {
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "impatience",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}},
           "patience": {"Class1": {"distribution": {"type": "Exp", "params": {"lambda": 0.25}},
                                    "impatienceType": "RENEGING"}},
           "orbitImpatience": {"Class1": {"type": "Exp", "params": {"lambda": 0.1}}},
           "batchRejectProb": {"Class1": 0.2},
           "balking": {"Class1": {"strategy": "QUEUE_LENGTH",
                                   "thresholds": [{"minJobs": 3, "maxJobs": -1,
                                                   "probability": 0.5}]}},
           "retrial": {"Class1": {"delay": {"type": "Exp", "params": {"lambda": 2.0}},
                                   "maxAttempts": 4}}},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}}
      }
    })";
    const json out = roundtrip(kJson);
    const json& q = node_of(out, "Queue1");
    CHECK(q.at("patience").at("Class1").at("impatienceType") == "RENEGING");
    CHECK(q.at("orbitImpatience").at("Class1").at("params").at("lambda") == doctest::Approx(0.1));
    CHECK(q.at("batchRejectProb").at("Class1") == doctest::Approx(0.2));
    CHECK(q.at("balking").at("Class1").at("strategy") == "QUEUE_LENGTH");
    CHECK(q.at("balking").at("Class1").at("thresholds")[0].at("maxJobs") == -1.0);
    CHECK(q.at("retrial").at("Class1").at("maxAttempts") == 4);
}

TEST_CASE("setup, delay-off and the heterogeneous pools survive") {
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "hetero",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS", "servers": 3,
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}},
           "setupTime": {"Class1": {"type": "Exp", "params": {"lambda": 4.0}}},
           "delayOffTime": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}},
           "serverTypes": [
             {"name": "fast", "count": 1, "compatibleClasses": ["Class1"],
              "service": {"Class1": {"type": "Exp", "params": {"lambda": 2.0}}}},
             {"name": "slow", "count": 2,
              "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}}],
           "heteroSchedPolicy": "ALIS"},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}}
      }
    })";
    const json out = roundtrip(kJson);
    const json& q = node_of(out, "Queue1");
    CHECK(q.at("servers") == 3.0);
    CHECK(q.at("setupTime").at("Class1").at("params").at("lambda") == doctest::Approx(4.0));
    CHECK(q.at("delayOffTime").at("Class1").at("params").at("lambda") == doctest::Approx(0.5));
    REQUIRE(q.at("serverTypes").size() == 2);
    CHECK(q.at("serverTypes")[0].at("name") == "fast");
    CHECK(q.at("serverTypes")[1].at("count") == 2.0);
    CHECK(q.at("heteroSchedPolicy") == "ALIS");
}

TEST_CASE("the distribution families the writers emit all reconstruct") {
    // One service law per family, at its own class, in the exact spelling the
    // reference writers use -- including the two `alpha` keys that mean the
    // SHAPE for Pareto and the SCALE for Weibull.
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "dists",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {
             "C1": {"type": "Gamma", "params": {"alpha": 2.0, "beta": 3.0}},
             "C2": {"type": "Pareto", "params": {"alpha": 3.0, "scale": 1.0}},
             "C3": {"type": "Weibull", "params": {"alpha": 1.0, "beta": 2.0}},
             "C4": {"type": "Coxian", "params": {"mu": [2.0, 3.0], "phi": [0.4, 1.0]}},
             "C5": {"type": "Geometric", "params": {"p": 0.25}},
             "C6": {"type": "HyperExp", "params": {"p": [0.2, 0.3, 0.5],
                                                    "lambda": [1.0, 2.0, 4.0]}}}},
          {"name": "Queue1", "type": "Queue", "scheduling": "PS",
           "service": {
             "C1": {"type": "Exp", "params": {"lambda": 1.0}},
             "C2": {"type": "Exp", "params": {"lambda": 1.0}},
             "C3": {"type": "Exp", "params": {"lambda": 1.0}},
             "C4": {"type": "Exp", "params": {"lambda": 1.0}},
             "C5": {"type": "Exp", "params": {"lambda": 1.0}},
             "C6": {"type": "Exp", "params": {"lambda": 1.0}}}}
        ],
        "classes": [
          {"name": "C1", "type": "Closed", "population": 1, "refNode": "Delay"},
          {"name": "C2", "type": "Closed", "population": 1, "refNode": "Delay"},
          {"name": "C3", "type": "Closed", "population": 1, "refNode": "Delay"},
          {"name": "C4", "type": "Closed", "population": 1, "refNode": "Delay"},
          {"name": "C5", "type": "Closed", "population": 1, "refNode": "Delay"},
          {"name": "C6", "type": "Closed", "population": 1, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {
          "C1,C1": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}},
          "C2,C2": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}},
          "C3,C3": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}},
          "C4,C4": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}},
          "C5,C5": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}},
          "C6,C6": {"Delay": {"Queue1": 1.0}, "Queue1": {"Delay": 1.0}}}}
      }
    })";
    const json out = roundtrip(kJson);
    const json& svc = node_of(out, "Delay").at("service");
    CHECK(svc.at("C1").at("type") == "Gamma");
    CHECK(svc.at("C1").at("params").at("alpha") == doctest::Approx(2.0));
    CHECK(svc.at("C1").at("params").at("beta") == doctest::Approx(3.0));
    CHECK(svc.at("C2").at("type") == "Pareto");
    CHECK(svc.at("C2").at("params").at("alpha") == doctest::Approx(3.0));
    CHECK(svc.at("C3").at("type") == "Weibull");
    CHECK(svc.at("C3").at("params").at("alpha") == doctest::Approx(1.0));
    CHECK(svc.at("C3").at("params").at("beta") == doctest::Approx(2.0));
    CHECK(svc.at("C4").at("type") == "Coxian");
    CHECK(svc.at("C4").at("params").at("mu")[1] == doctest::Approx(3.0));
    CHECK(svc.at("C5").at("type") == "Geometric");
    CHECK(svc.at("C5").at("params").at("p") == doctest::Approx(0.25));
    CHECK(svc.at("C6").at("type") == "HyperExp");
    REQUIRE(svc.at("C6").at("params").at("p").size() == 3);
    CHECK(svc.at("C6").at("params").at("lambda")[2] == doctest::Approx(4.0));
}

TEST_CASE("the moments of the discrete families are MATLAB's") {
    // Geometric on the NUMBER OF TRIALS (mean 1/p, SCV 1-p), Poisson's SCV of
    // 1/lambda rather than the exponential's 1, and the DiscreteUniform's
    // ((b-a+1)^2-1)/12 variance: three places where the continuous analogue
    // would give a different number and no error.
    const lang::Distrib<double> g = lang::Distrib<double>::geometric(0.25);
    CHECK(g.mean == doctest::Approx(4.0));
    CHECK(g.scv == doctest::Approx(0.75));
    const lang::Distrib<double> p = lang::Distrib<double>::poisson(4.0);
    CHECK(p.mean == doctest::Approx(4.0));
    CHECK(p.scv == doctest::Approx(0.25));
    const lang::Distrib<double> du = lang::Distrib<double>::discrete_uniform(1.0, 6.0);
    CHECK(du.mean == doctest::Approx(3.5));
    CHECK(du.scv == doctest::Approx((36.0 - 1.0) / 12.0 / (3.5 * 3.5)));
    const lang::Distrib<double> b = lang::Distrib<double>::binomial(10.0, 0.3);
    CHECK(b.mean == doctest::Approx(3.0));
    CHECK(b.scv == doctest::Approx(0.7 / 3.0));
}

TEST_CASE("a DMAP takes the discrete moments, not the continuous ones") {
    // D0 + D1 is STOCHASTIC, so the continuous map_scv the MATLAB class
    // inherited reads it as a generator and solves a singular system. The
    // discrete law here is a geometric interarrival of parameter 0.4: one
    // phase, D0 = [0.6], D1 = [0.4].
    Matrix<double> D0(1, 1), D1(1, 1);
    D0(0, 0) = 0.6;
    D1(0, 0) = 0.4;
    lang::Distrib<double> d = lang::Distrib<double>::dmap(D0, D1);
    lang::dist_refresh_moments(d);
    CHECK(d.mean == doctest::Approx(1.0 / 0.4));
    CHECK(d.scv == doctest::Approx(1.0 - 0.4));
}

TEST_CASE("a Logger keeps the file it writes") {
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "logged",
        "logPath": "/tmp",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 0.5}}}},
          {"name": "Log", "type": "Logger", "fileName": "trace.csv"},
          {"name": "Queue1", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}},
          {"name": "Sink", "type": "Sink"}
        ],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Log": 1.0}, "Log": {"Queue1": 1.0}, "Queue1": {"Sink": 1.0}}}}
      }
    })";
    const json out = roundtrip(kJson);
    const json& lg = node_of(out, "Log");
    CHECK(lg.at("type") == "Logger");
    CHECK(lg.at("fileName") == "trace.csv");
}

TEST_CASE("a LayeredNetwork model.json reaches LqnStruct") {
    // The interchange a `save_model` on a LayeredNetwork writes. It is NOT
    // .lqnx: a task think time, a CacheTask, an ItemEntry and a SetupTask all
    // cross here and none of them can cross the XML, so a layered model that
    // used any of them could previously reach this port only by losing it.
    static const char* kJson = R"({
      "format": "line-model", "version": "1.0",
      "model": {
        "type": "LayeredNetwork", "name": "lqn",
        "hosts": [
          {"name": "P1", "multiplicity": 2147483647, "scheduling": "INF"},
          {"name": "P2", "scheduling": "PS"}
        ],
        "tasks": [
          {"name": "T1", "host": "P1", "scheduling": "REF", "multiplicity": 5,
           "thinkTime": {"type": "Exp", "params": {"lambda": 0.1}}},
          {"name": "T2", "host": "P2", "scheduling": "FCFS"}
        ],
        "entries": [
          {"name": "E1", "task": "T1"},
          {"name": "E2", "task": "T2"}
        ],
        "activities": [
          {"name": "A1", "task": "T1", "boundToEntry": "E1", "repliesTo": "E1",
           "hostDemand": {"type": "Exp", "params": {"lambda": 1.0}},
           "synchCalls": [{"dest": "E2", "mean": 2.0}]},
          {"name": "A2", "task": "T2", "boundToEntry": "E2", "repliesTo": "E2",
           "hostDemand": {"type": "Exp", "params": {"lambda": 0.5}}}
        ]
      }
    })";
    const json root = json::parse(kJson);
    const lqn::LqnStruct<double> m = io::build_lqn_from_json<double>(root);
    CHECK(m.nhosts == 2);
    CHECK(m.ntasks == 2);
    CHECK(m.nentries == 2);
    CHECK(m.nacts == 2);
    // `2147483647` is the wire's INFINITY, not a count of two billion servers.
    CHECK(std::isinf(m.mult[1]));
}

TEST_CASE("a layered model.json is told apart from an .lqnx by its content") {
    // `linemodel_save` writes `.json` for BOTH model kinds, so the extension
    // cannot decide which reader applies and the first character plus the
    // `type` field has to.
    const std::string path =
        std::string(LINE_MP_SCRATCH) + "/roundtrip_layered_probe.json";
    {
        std::ofstream out(path.c_str());
        out << R"({"model": {"type": "LayeredNetwork", "name": "p", "hosts": [],
                              "tasks": [], "entries": [], "activities": []}})";
    }
    CHECK(io::is_layered_json(path));
    {
        std::ofstream out(path.c_str());
        out << R"({"model": {"type": "Network", "name": "p"}})";
    }
    CHECK_FALSE(io::is_layered_json(path));
    std::remove(path.c_str());
}

TEST_CASE("an unknown model or node key is named, not dropped") {
    static const char* kJson = R"({
      "model": {
        "type": "Network", "name": "unknown-key",
        "nodes": [
          {"name": "Delay", "type": "Delay", "quantumOfSolace": 7,
           "service": {"Class1": {"type": "Exp", "params": {"lambda": 1.0}}}}
        ],
        "classes": [{"name": "Class1", "type": "Closed", "population": 1,
                     "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Delay": 1.0}}}}
      }
    })";
    CHECK_THROWS_AS(roundtrip(kJson), UnsupportedError);
}

}  // namespace
