/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Place and Transition nodes over the model.json wire.
 *
 * THE JSON BELOW IS MATLAB'S OWN, byte for byte: `linemodel_save` was run on
 * the 3-place cyclic net at N = 4 with firing rates {1, 1.5, 2} and its output
 * pasted here. A hand-written fixture would only prove the reader agrees with
 * this file's author about the schema; the writer's output proves it agrees
 * with the writer, which is what makes a C++ SPN reachable from a model the
 * other codebases produced.
 *
 * WHAT IT PINS. Reading must land the same net the programmatic builder does in
 * `test_spn_mdd.cpp`: |S| = 15 and the queue lengths MATLAB and python report,
 * [2.249874392899, 1.069837548149, 0.680288058952]. That closes the loop
 * writer -> reader -> descriptor -> measures without a single hand-entered
 * parameter.
 */

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mdd/mdd_mcd.h"
#include "line/api/spn/spn_mdd.h"
#include "line/io/network_reader.h"
#include "line/num/number.h"

using namespace line;

namespace {

const double QLEN_REF[3] = {2.249874392899, 1.069837548149, 0.680288058952};

const char* kCyclicSpnJson = R"JSON({"format":"line-model","model":{"classes":[{"name":"Class1","population":4,"refNode":"P0","type":"Closed"}],"name":"spn","nodes":[{"dropRule":{"Class1":"waitingQueue"},"name":"P0","type":"Place"},{"dropRule":{"Class1":"waitingQueue"},"name":"P1","type":"Place"},{"dropRule":{"Class1":"waitingQueue"},"name":"P2","type":"Place"},{"modes":[{"distribution":{"params":{"lambda":1},"type":"Exp"},"enablingConditions":[{"class":"Class1","count":1,"node":"P0"}],"firingOutcomes":[{"class":"Class1","count":1,"node":"P1"}],"firingPriority":1,"name":"M0","timingStrategy":"TIMED"}],"name":"T0","type":"Transition"},{"modes":[{"distribution":{"params":{"lambda":1.5},"type":"Exp"},"enablingConditions":[{"class":"Class1","count":1,"node":"P1"}],"firingOutcomes":[{"class":"Class1","count":1,"node":"P2"}],"firingPriority":1,"name":"M1","timingStrategy":"TIMED"}],"name":"T1","type":"Transition"},{"modes":[{"distribution":{"params":{"lambda":2},"type":"Exp"},"enablingConditions":[{"class":"Class1","count":1,"node":"P2"}],"firingOutcomes":[{"class":"Class1","count":1,"node":"P0"}],"firingPriority":1,"name":"M2","timingStrategy":"TIMED"}],"name":"T2","type":"Transition"}],"routing":{"matrix":{"Class1,Class1":{"P0":{"T0":1},"P1":{"T1":1},"P2":{"T2":1},"T0":{"P1":1},"T1":{"P2":1},"T2":{"P0":1}}},"type":"matrix"},"type":"Network"},"version":"1.0"})JSON";

/** The reader takes a path, so the fixture is materialized next to the test. */
std::string write_temp_json(const char* text) {
    const std::string path = std::string(LINE_MP_REPO_ROOT) + "/cpp/tests/.spn_reader_fixture.json";
    std::ofstream out(path.c_str());
    if (!out) throw line::InputError("test: cannot write the fixture");
    out << text;
    out.close();
    return path;
}

}  // namespace

TEST_CASE("a Petri net written by MATLAB reads back and solves") {
    const std::string path = write_temp_json(kCyclicSpnJson);
    qn::Network<double> m = io::read_network_json<double>(path);
    std::remove(path.c_str());

    const spn::SpnResult<double> r = spn::spn_mdd(m.get_struct());
    CHECK(r.info.nplacelevels == 3);
    CHECK(r.info.diagram.cardinality() == 15);
    CHECK(r.info.placenames[0] == "P0");

    const mdd::MddMcdResult<double> out = mdd::mdd_mcd(r.mdds, r.desc);
    for (std::size_t l = 0; l < 3; ++l)
        CHECK(out.QLen[l] == doctest::Approx(QLEN_REF[l]).epsilon(1e-9));
}

/**
 * A mode with NO `firingPriority` key reads back at the builders' default.
 *
 * THE FIXTURE ABOVE CANNOT CATCH THIS, and that is why the defect survived: the
 * writers emit the key whenever the priority is positive, so MATLAB's own output
 * carries `"firingPriority":1` on every mode and the reader's default is never
 * exercised by it. The key IS omitted by a producer that leaves it out, and the
 * reader then decides what the model means.
 *
 * The number is 1 because `Transition.addMode` appends 1 in all four codebases
 * (MATLAB Transition.m:88, verified in a live session; python nodes.py:2660;
 * JAR Transition.java:114) and the MATLAB, python and JAR readers all leave that
 * default standing when the key is absent. This reader used to answer 0, which
 * is a different net: firing priority arbitrates which IMMEDIATE mode fires, so
 * the marking process changes rather than a reported decimal.
 */
TEST_CASE("a mode with no firingPriority key reads back at the builder default") {
    const char* kNoPrioJson =
        R"JSON({"format":"line-model","model":{"classes":[{"name":"Class1","population":1,)JSON"
        R"JSON("refNode":"P0","type":"Closed"}],"name":"spn","nodes":[{"name":"P0","type":"Place"},)JSON"
        R"JSON({"name":"P1","type":"Place"},{"modes":[{"distribution":{"params":{"lambda":1},)JSON"
        R"JSON("type":"Exp"},"enablingConditions":[{"class":"Class1","count":1,"node":"P0"}],)JSON"
        R"JSON("firingOutcomes":[{"class":"Class1","count":1,"node":"P1"}],"name":"M0",)JSON"
        R"JSON("timingStrategy":"TIMED"}],"name":"T0","type":"Transition"},{"modes":[{)JSON"
        R"JSON("distribution":{"params":{"lambda":2},"type":"Exp"},"enablingConditions":)JSON"
        R"JSON([{"class":"Class1","count":1,"node":"P1"}],"firingOutcomes":[{"class":"Class1",)JSON"
        R"JSON("count":1,"node":"P0"}],"firingPriority":3,"name":"M1","timingStrategy":"TIMED"}],)JSON"
        R"JSON("name":"T1","type":"Transition"}],"routing":{"matrix":{"Class1,Class1":)JSON"
        R"JSON({"P0":{"T0":1},"P1":{"T1":1},"T0":{"P1":1},"T1":{"P0":1}}},"type":"matrix"},)JSON"
        R"JSON("type":"Network"},"version":"1.0"})JSON";

    const std::string path = write_temp_json(kNoPrioJson);
    qn::Network<double> m = io::read_network_json<double>(path);
    std::remove(path.c_str());
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // T0 omits the key and must land on 1; T1 states 3 and must keep it, so the
    // default cannot be restored by clamping every priority to a constant.
    std::size_t t0 = 0, t1 = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        if (sn.nodes[i].name == "T0") t0 = i + 1;
        if (sn.nodes[i].name == "T1") t1 = i + 1;
    }
    REQUIRE(t0 != 0);
    REQUIRE(t1 != 0);
    REQUIRE(sn.transparam.count(t0) == 1);
    REQUIRE(sn.transparam.count(t1) == 1);
    CHECK(sn.transparam.at(t0).firingprio[0] == doctest::Approx(1.0));
    CHECK(sn.transparam.at(t1).firingprio[0] == doctest::Approx(3.0));

    // The two neighbouring defaults travel with it and are asserted here rather
    // than in their own case: all three are read off the same absent-key path.
    CHECK(sn.transparam.at(t0).nmodeservers[0] == doctest::Approx(1.0));
    CHECK(num_traits<double>::to_double(sn.transparam.at(t0).fireweight[0]) ==
          doctest::Approx(1.0));
}
