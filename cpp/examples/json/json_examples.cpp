/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/json/` and `python/examples/json/`: the six models the
 * cross-codebase JSON round-trip is measured on.
 *
 * A REFERENCE FILE HERE IS A MODEL FACTORY AND NOTHING ELSE. Neither
 * `json_balking_retrial.py` nor its MATLAB twin constructs a solver, prints a
 * table or reads a file: each returns a Network whose only purpose is to carry
 * one feature group across `linemodel_save` / `linemodel_load`, and the harness
 * (`line-test.git/parity/test_json.sh`) is what saves, reloads and diffs it.
 * So there is no `model.json` in the tree for these, and there is nothing for
 * this port to READ: the C++ counterparts build the same models and report
 * `model.getStruct()`, which is the quantity the round-trip compares. Running a
 * solver here would answer a question the reference never asks.
 *
 * The C++ engine carries BOTH halves of the wire format: a reader
 * (`io::read_network_json`, `include/line/io/network_reader.h`) and a writer
 * (`io::write_network_json`, `include/line/io/network_writer.h`), so the round
 * trip can be closed in this port. The reader side is additionally exercised by
 * `line-cli -i json` on the file the Python exporter writes, which is what makes
 * the round trip a CROSS-codebase check rather than a self-consistency one.
 */

#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"

namespace line {
namespace examples {

namespace {

const char* kGroup = "json";

/**
 * `model.getStruct()`, the quantity a JSON round trip is compared on.
 *
 * The station and class tables are the part every one of the six models has;
 * the feature each example exists for is printed by the example itself, since
 * only it knows which field it set.
 */
void struct_summary(Net& m) {
    const Sn& sn = m.get_struct();
    kv("Model", sn.name);
    kv("Nodes", static_cast<double>(sn.nof_nodes()));
    kv("Stations", static_cast<double>(sn.nstations));
    kv("Classes", static_cast<double>(sn.nclasses));
    kv("Chains", static_cast<double>(sn.nchains));
    std::printf("%-16s %-10s %10s %10s\n", "Station", "Sched", "Servers", "Capacity");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const double c = sn.stations[i].cap;
        char cbuf[24];
        if (std::isinf(c)) std::snprintf(cbuf, sizeof cbuf, "%10s", "Inf");
        else std::snprintf(cbuf, sizeof cbuf, "%10g", c);
        char sbuf[24];
        if (std::isinf(sn.stations[i].nservers)) std::snprintf(sbuf, sizeof sbuf, "%10s", "Inf");
        else std::snprintf(sbuf, sizeof sbuf, "%10g", sn.stations[i].nservers);
        std::printf("%-16s %-10s %s %s\n", sn.stations[i].name.c_str(),
                    lang::sched_to_text(sn.stations[i].sched), sbuf, cbuf);
    }
    std::printf("%-16s %-10s %12s\n", "JobClass", "Type", "Population");
    for (std::size_t c = 0; c < sn.nclasses; ++c) {
        const bool open = sn.classes[c].type == lang::JobClassType::OPEN;
        char pbuf[24];
        if (open) std::snprintf(pbuf, sizeof pbuf, "%12s", "Inf");
        else std::snprintf(pbuf, sizeof pbuf, "%12g", sn.classes[c].population);
        std::printf("%-16s %-10s %s\n", sn.classes[c].name.c_str(), open ? "open" : "closed",
                    pbuf);
    }
}

// ---------------------------------------------------------------------------
// json_balking_retrial
// ---------------------------------------------------------------------------

/** Open 2-class network with balking, retrial, and patience. */
void json_balking_retrial() {
    Net m("Balking_Retrial");

    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");

    queue.set_capacity(15.0);

    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");

    source.set_arrival(c1, Exp(1.0));
    source.set_arrival(c2, Exp(0.5));
    queue.set_service(c1, Exp(2.0));
    queue.set_service(c2, Exp(3.0));

    // Class1: balking based on queue length
    std::vector<qn::Station<double>::BalkingThreshold> thr(2);
    thr[0].min_jobs = 5.0;
    thr[0].max_jobs = 10.0;
    thr[0].probability = 0.3;
    thr[1].min_jobs = 11.0;
    thr[1].max_jobs = std::numeric_limits<double>::infinity();
    thr[1].probability = 1.0;
    m.set_balking(queue, c1, lang::BalkingStrategy::QUEUE_LENGTH, thr);

    // Class1: patience (reneging)
    m.set_patience(queue, c1, Exp(0.1));

    // Class2: retrial with max attempts
    m.set_retrial(queue, c2, Exp(0.5), 0.5, 3);

    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);

    struct_summary(m);
    kv("Balking (Class1)", std::string("QUEUE_LENGTH, [5,10]->0.3, [11,Inf)->1.0"));
    kv("Patience (Class1)", std::string("Exp(0.1)"));
    kv("Retrial (Class2)", std::string("Exp(0.5), at most 3 attempts"));
}
LINE_EXAMPLE(kGroup, json_balking_retrial);

// ---------------------------------------------------------------------------
// json_classcap_droprule
// ---------------------------------------------------------------------------

/** Open 2-class network with per-class capacity, drop rules, and load-dependence. */
void json_classcap_droprule() {
    Net m("ClassCap_DropRule");

    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::PS);
    Sink sink(m, "Sink");

    queue.set_number_of_servers(2.0);
    queue.set_capacity(20.0);

    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");

    source.set_arrival(c1, Exp(1.0));
    source.set_arrival(c2, Exp(0.5));
    queue.set_service(c1, Exp(3.0));
    queue.set_service(c2, Exp(2.0));

    // Per-class capacity
    queue.set_class_capacity(c1, 8.0);
    queue.set_class_capacity(c2, 15.0);

    // Per-class drop rules. WAITQ is not usable here: no solver honours "wait
    // upstream" for an open class at a plain finite capacity, so the capacity
    // refresh rejects that combination. BAS is the honoured blocking policy.
    queue.set_drop_rule(c1, lang::DropStrategy::DROP);
    queue.set_drop_rule(c2, lang::DropStrategy::BAS);

    // Load-dependent scaling
    queue.set_load_dependence(std::vector<double>{1.0, 0.9, 0.8, 0.7});

    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);

    struct_summary(m);
    kv("Class capacity", std::string("Class1=8, Class2=15"));
    kv("Drop rule", std::string("Class1=DROP, Class2=BAS"));
    kv("Load dependence", std::string("[1, 0.9, 0.8, 0.7]"));
}
LINE_EXAMPLE(kGroup, json_classcap_droprule);

// ---------------------------------------------------------------------------
// json_fcr_details
// ---------------------------------------------------------------------------

/** Open 2-class network with a Finite Capacity Region. */
void json_fcr_details() {
    Net m("FCR_Details");

    Source source(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    Sink sink(m, "Sink");

    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");

    source.set_arrival(c1, Exp(1.0));
    source.set_arrival(c2, Exp(0.5));
    q1.set_service(c1, Exp(3.0));
    q1.set_service(c2, Exp(2.0));
    q2.set_service(c1, Exp(4.0));
    q2.set_service(c2, Exp(3.0));

    Routing P;
    serial(P, c1, {source, q1, q2, sink});
    serial(P, c2, {source, q1, q2, sink});
    m.link(P);

    // Finite Capacity Region: per-class job caps, a global cap, per-class
    // weights against that global cap, and a per-job footprint.
    const std::size_t fcr =
        m.add_region({q1, q2}, std::vector<double>{8.0, 10.0}, 15.0,
                     std::vector<lang::DropStrategy>(), std::vector<double>(),
                     std::vector<double>{1.0, 3.0});
    m.set_region_weights(fcr, std::vector<double>{1.0, 2.0});

    struct_summary(m);
    kv("Region members", std::string("Queue1, Queue2"));
    kv("Region global max jobs", 15.0);
    kv("Region class max jobs", std::string("Class1=8, Class2=10"));
    kv("Region class weight", std::string("Class1=1, Class2=2"));
    kv("Region class size", std::string("Class1=1, Class2=3"));
}
LINE_EXAMPLE(kGroup, json_fcr_details);

// ---------------------------------------------------------------------------
// json_hetero_servers
// ---------------------------------------------------------------------------

/** Closed 2-class network with heterogeneous servers. */
void json_hetero_servers() {
    Net m("HeteroServers");

    Delay delay(m, "Delay");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    queue.set_number_of_servers(3.0);

    ClosedClass c1(m, "Class1", 3.0, delay);
    ClosedClass c2(m, "Class2", 2.0, delay);

    delay.set_service(c1, Exp(1.0));
    delay.set_service(c2, Exp(2.0));

    // Default service (required before hetero setup)
    queue.set_service(c1, Exp(3.0));
    queue.set_service(c2, Exp(2.0));

    // Heterogeneous server types: Fast serves both classes, Slow only Class1.
    qn::Station<double>::ServerType fast;
    fast.name = "Fast";
    fast.count = 2.0;
    fast.compatible = std::vector<bool>{true, true};
    fast.service = std::vector<D>{Exp(5.0), Exp(3.0)};
    qn::Station<double>::ServerType slow;
    slow.name = "Slow";
    slow.count = 1.0;
    slow.compatible = std::vector<bool>{true, false};
    slow.service = std::vector<D>{Exp(1.0), Disabled()};

    m.add_server_type(queue, fast);
    m.add_server_type(queue, slow);
    m.set_hetero_sched_policy(queue, lang::HeteroSchedPolicy::ORDER);

    Routing P;
    cyclic(P, c1, {delay, queue});
    cyclic(P, c2, {delay, queue});
    m.link(P);

    struct_summary(m);
    kv("Server type Fast", std::string("2 servers, Class1 Exp(5), Class2 Exp(3)"));
    kv("Server type Slow", std::string("1 server, Class1 Exp(1)"));
    kv("Hetero policy", std::string("ORDER"));
}
LINE_EXAMPLE(kGroup, json_hetero_servers);

// ---------------------------------------------------------------------------
// json_join_deadline
// ---------------------------------------------------------------------------

/**
 * Closed 1-class fork-join with quorum join and deadline.
 *
 * `class1.setDeadline(5.0)` HAS NO C++ COUNTERPART: no field of
 * `qn::JobClass` carries a per-class deadline and no reader key sets one, so
 * the deadline is reported as absent rather than dropped in silence. The join
 * rule itself is expressible: QUORUM is the JAR's spelling of PARTIAL, and the
 * quorum count is `set_join_strategy`'s second argument.
 */
void json_join_deadline() {
    Net m("JoinDeadline");

    Delay delay(m, "Delay");
    const std::size_t fork = m.add_fork("Fork");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Queue q3(m, "Queue3", SchedStrategy::FCFS);
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "Class1", 5.0, delay);

    delay.set_service(c1, Exp(1.0));
    q1.set_service(c1, Exp(3.0));
    q2.set_service(c1, Exp(4.0));
    q3.set_service(c1, Exp(2.0));

    m.set_join_strategy(join, lang::JoinStrategy::PARTIAL, 2.0);

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, q1, 1.0);
    P.set(c1, c1, fork, q2, 1.0);
    P.set(c1, c1, fork, q3, 1.0);
    P.set(c1, c1, q1, join, 1.0);
    P.set(c1, c1, q2, join, 1.0);
    P.set(c1, c1, q3, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);

    struct_summary(m);
    kv("Join strategy", std::string("PARTIAL (QUORUM), 2 of 3 siblings"));
    na("Class deadline", "ClosedClass.setDeadline(5.0) has no field in the C++ NetworkStruct");
}
LINE_EXAMPLE(kGroup, json_join_deadline);

// ---------------------------------------------------------------------------
// json_signal_classes
// ---------------------------------------------------------------------------

/** Open network with G-network negative signal. */
void json_signal_classes() {
    Net m("SignalClasses");

    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");

    OpenClass c1(m, "Class1");
    OpenClass s1(m, "Signal1");
    // `OpenSignal(model, 'Signal1', NEGATIVE).forJobClass(Class1)`: the signal is
    // an ordinary open class whose target is the class it removes a job from.
    m.set_signal(s1, lang::SignalType::NEGATIVE, lang::RemovalPolicy::RANDOM, c1);

    source.set_arrival(c1, Exp(2.0));
    source.set_arrival(s1, Exp(0.5));
    queue.set_service(c1, Exp(5.0));
    queue.set_service(s1, Immediate());

    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, s1, {source, queue, sink});
    m.link(P);

    struct_summary(m);
    kv("Signal1 type", std::string("NEGATIVE"));
    kv("Signal1 target", std::string("Class1"));
    kv("Signal1 removal", std::string("RANDOM"));
}
LINE_EXAMPLE(kGroup, json_signal_classes);

}  // namespace

}  // namespace examples
}  // namespace line
