/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#include <sstream>
#include <string>

#include "doctest.h"
#include "line/lang/distributions.h"
#include "line/lang/qn/nodes.h"
#include "line/solvers/solver.h"

using namespace line;
using namespace line::lang;

TEST_CASE("public syntax: serial routing returns a linkable class block") {
    Network model("routing");
    Source source(model, "Source");
    Queue queue(model, "Queue", SchedStrategy::FCFS);
    Sink sink(model, "Sink");
    OpenClass first(model, "First");
    OpenClass second(model, "Second");

    const Routing block = model.serial_routing({source, queue, sink});
    CHECK(block.get(1, 1, source, queue) == doctest::Approx(1.0));
    CHECK(block.get(1, 1, queue, sink) == doctest::Approx(1.0));
    CHECK(block.get(1, 1, sink, source) == doctest::Approx(0.0));

    Routing P = model.init_routing_matrix();
    P.set(first, block);
    P.set(second, block);
    CHECK(P.get(first, first, source, queue) == doctest::Approx(1.0));
    CHECK(P.get(second, second, queue, sink) == doctest::Approx(1.0));
    model.link(P);
}

TEST_CASE("public syntax: serial routing closes a path without a Sink") {
    Network model("cycle");
    Delay delay(model, "Delay");
    Queue queue(model, "Queue", SchedStrategy::PS);
    ClosedClass jobs(model, "Jobs", 2.0, delay);

    const Routing P = model.serial_routing({delay, queue});
    CHECK(P.get(jobs, jobs, delay, queue) == doctest::Approx(1.0));
    CHECK(P.get(jobs, jobs, queue, delay) == doctest::Approx(1.0));
}

TEST_CASE("public syntax: Replayer accepts the same trace path as MATLAB") {
    const std::string path =
        std::string(LINE_MP_REPO_ROOT) + "/matlab/examples/gettingstarted/example_trace.txt";
    const Replayer trace(path);
    CHECK(!trace.trace.empty());
    CHECK(trace.trace_file == path);
    CHECK(trace.mean > 0.0);
}

TEST_CASE("public syntax: AvgTable filters and prints by model handles") {
    Network model("table");
    Source source(model, "Source");
    Queue queue(model, "Queue", SchedStrategy::FCFS);
    OpenClass jobs(model, "Jobs");

    AvgTable table;
    table.Station = {"Source", "Queue"};
    table.JobClass = {"Jobs", "Jobs"};
    table.QLen = {0.0, 0.5};
    table.Util = {0.0, 0.5};
    table.RespT = {0.0, 1.0};
    table.ResidT = {0.0, 1.0};
    table.ArvR = {0.0, 0.5};
    table.Tput = {0.5, 0.5};

    CHECK(table(queue, jobs).size() == 1);
    CHECK(table(jobs, queue).get("QLen", "Queue", "Jobs") == doctest::Approx(0.5));
    CHECK(table.filter_by(source).Station[0] == "Source");
    CHECK(table.get(jobs).size() == 2);
    CHECK(table.tget(queue, jobs).Tput[0] == doctest::Approx(0.5));

    std::ostringstream out;
    out << table(queue, jobs);
    CHECK(out.str().find("Station") != std::string::npos);
    CHECK(out.str().find("Queue") != std::string::npos);
}

