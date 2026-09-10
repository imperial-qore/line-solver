/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverLQNS: the .lqnx writer, and the wrapper around the external binary.
 *
 * TWO HALVES, AND ONLY ONE OF THEM CAN ALWAYS RUN. The writer and the reply
 * inference are this port's own code and are checked unconditionally, by
 * round-tripping a model through the document and comparing the struct it comes
 * back as. The wrapper needs an lqns installation, which LINE does not ship and
 * may not redistribute, so those cases are SKIPPED by name where the binary is
 * absent -- never silently passed, since a green suite that exercised nothing
 * is the failure mode that matters here.
 *
 * THE NUMERIC FIXTURE IS ANALYTIC, not read back out of lqns. A closed layered
 * model with ONE customer has cycle time Z + sum of demands exactly, whatever
 * the layering, so throughput and every utilization follow in closed form; a
 * wrapper that mis-mapped a column or lost the think time would still produce
 * plausible numbers and would fail this.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/solvers/wrappers/lqns/solver_lqns.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

const double kInf = std::numeric_limits<double>::infinity();

/** Z = 2, one activity of 0.3 calling an entry of 1.0, a single customer. */
lqn::LqnBuilder<double> tandem_builder() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", kInf, SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::FCFS);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", E(2.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", E(0.3), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    return b;
}

std::string temp_path(const char* stem) {
    const char* base = std::getenv("TMPDIR");
    return (base && *base ? std::string(base) : std::string("/tmp")) + "/" + stem;
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& name) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.names[i] == name) return i;
    return 0;
}

bool lqns_here() { return lqns::lqns_is_available(); }

}  // namespace

TEST_CASE("lqnx writer: a built model round-trips through the document") {
    const lqn::LqnModel<double> m = tandem_builder().model();
    const std::string path = temp_path("line_lqns_roundtrip.lqnx");
    const lqn::LqnWriteReport rep = lqn::write_lqnx(m, path, "roundtrip");
    CHECK(rep.dropped.empty());

    const lqn::LqnStruct<double> before = lqn::lqn_finalize(m);
    const lqn::LqnStruct<double> after = lqn::read_lqnx<double>(path);
    REQUIRE(after.nidx == before.nidx);
    REQUIRE(after.ncalls == before.ncalls);
    REQUIRE(after.nhosts == before.nhosts);
    REQUIRE(after.ntasks == before.ntasks);
    for (std::size_t i = 1; i <= before.nidx; ++i) {
        CHECK(after.names[i] == before.names[i]);
        CHECK(after.type[i] == before.type[i]);
        CHECK(after.parent[i] == before.parent[i]);
        CHECK(after.hostdem[i].mean == doctest::Approx(before.hostdem[i].mean));
    }
    for (std::size_t t = before.tshift + 1; t <= before.tshift + before.ntasks; ++t) {
        CHECK(after.sched[t] == before.sched[t]);
        CHECK(after.mult[t] == doctest::Approx(before.mult[t]));
        CHECK(after.think[t].mean == doctest::Approx(before.think[t].mean));
    }
    for (std::size_t c = 1; c <= before.ncalls; ++c) {
        CHECK(after.callnames[c] == before.callnames[c]);
        CHECK(after.calltype[c] == before.calltype[c]);
        CHECK(after.callproc_mean[c] == doctest::Approx(before.callproc_mean[c]));
    }
    std::remove(path.c_str());
}

TEST_CASE("lqnx writer: the reply lqns requires is inferred when none is declared") {
    // The builder above never calls replies_to, so a writer that only copied
    // declarations would emit a document lqns rejects.
    const lqn::LqnModel<double> m = tandem_builder().model();
    const std::string path = temp_path("line_lqns_reply.lqnx");
    lqn::write_lqnx(m, path, "reply");
    std::string text;
    {
        std::ifstream in(path.c_str());
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    CHECK(text.find("<reply-entry name=\"E2\">") != std::string::npos);
    CHECK(text.find("<reply-activity name=\"A2\"") != std::string::npos);
    // E1 belongs to the REFERENCE task, which replies to nobody.
    CHECK(text.find("<reply-entry name=\"E1\"") == std::string::npos);
    std::remove(path.c_str());
}

TEST_CASE("lqnx writer: a think time the schema cannot carry is reported, not dropped") {
    lqn::LqnBuilder<double> b = tandem_builder();
    b.think_time("T2", E(4.0));  // legal in the API, illegal on a non-ref task in .lqnx
    const std::string path = temp_path("line_lqns_think.lqnx");
    const lqn::LqnWriteReport rep = lqn::write_lqnx(b.model(), path, "think");
    REQUIRE(rep.dropped.size() == 1);
    CHECK(rep.dropped[0].find("T2") != std::string::npos);
    const lqn::LqnStruct<double> after = lqn::read_lqnx<double>(path);
    CHECK(after.think[idx_of(after, "T2")].mean == doctest::Approx(0.0));
    std::remove(path.c_str());
}

TEST_CASE("lqnx writer: a construct with no element in the schema is refused by name") {
    lqn::LqnBuilder<double> cache;
    cache.processor("P1", kInf, SchedStrategy::INF);
    cache.processor("P2", 1, SchedStrategy::FCFS);
    cache.task("T1", 1, SchedStrategy::REF, "P1");
    std::vector<int> cap(1, 1);
    cache.cache_task("C1", 1, SchedStrategy::FCFS, "P2", 3, cap, lang::ReplacementStrategy::LRU);
    cache.entry("E1", "T1");
    cache.activity("A1", E(1.0), "T1");
    cache.bound_to("A1", "E1");
    const std::string path = temp_path("line_lqns_cache.lqnx");
    CHECK_THROWS_AS(lqn::write_lqnx(cache.model(), path, "cache"), UnsupportedError);
}

TEST_CASE("lqnx writer: speed-factor survives the round trip") {
    // getStruct reads neither speed-factor nor quantum, so nothing downstream
    // would notice them being lost -- but lqns honours both, and a processor
    // silently written back at speed 1 answers a different model.
    lqn::LqnModel<double> m = tandem_builder().model();
    m.procs[1].speed_factor = 2.5;
    const std::string path = temp_path("line_lqns_speed.lqnx");
    lqn::write_lqnx(m, path, "speed");
    const lqn::LqnModel<double> back = lqn::read_lqnx_model<double>(path);
    CHECK(back.procs[1].speed_factor == doctest::Approx(2.5));
    std::remove(path.c_str());
}

TEST_CASE("lqnx writer: abstract names keep a shared name apart by kind") {
    // lqngen names a processor, its task and its entry alike, so a SINGLE name
    // map -- which is what the reference keeps -- collapses all three onto the
    // last one written and emits a document naming elements that do not exist.
    lqn::LqnBuilder<double> b;
    b.processor("c0", kInf, SchedStrategy::INF);
    b.processor("p0", 1, SchedStrategy::FCFS);
    b.task("c0", 1, SchedStrategy::REF, "c0");
    b.task("t0", 1, SchedStrategy::FCFS, "p0");
    b.entry("c0", "c0");
    b.entry("e0", "t0");
    b.activity("c0_1", E(1.0), "c0");
    b.bound_to("c0_1", "c0");
    b.sync_call("c0_1", "e0", 1.0);
    b.activity("e0_1", E(1.0), "t0");
    b.bound_to("e0_1", "e0");

    const std::string path = temp_path("line_lqns_abstract.lqnx");
    lqn::write_lqnx(b.model(), path, "abstract", true);
    std::string text;
    {
        std::ifstream in(path.c_str());
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    CHECK(text.find("<processor name=\"P1\"") != std::string::npos);
    CHECK(text.find("<task name=\"T1\"") != std::string::npos);
    CHECK(text.find("<entry name=\"E1\"") != std::string::npos);
    CHECK(text.find("bound-to-entry=\"E1\"") != std::string::npos);
    // The call names the ENTRY e0 (E2), not the task or processor named alike.
    CHECK(text.find("<synch-call dest=\"E2\"") != std::string::npos);
    // Renamed consistently, so the document still reads back as the same model.
    const lqn::LqnStruct<double> after = lqn::read_lqnx<double>(path);
    const lqn::LqnStruct<double> before = b.build();
    REQUIRE(after.nidx == before.nidx);
    REQUIRE(after.ncalls == before.ncalls);
    for (std::size_t i = 1; i <= before.nidx; ++i) {
        CHECK(after.type[i] == before.type[i]);
        CHECK(after.parent[i] == before.parent[i]);
        CHECK(after.hostdem[i].mean == doctest::Approx(before.hostdem[i].mean));
    }
    std::remove(path.c_str());
}

TEST_CASE("SolverLQNS: an unknown method or multiserver policy is refused") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    const lqn::LqnModel<double> m = tandem_builder().model();
    lqns::LqnsOptions bad;
    bad.method = "amva";
    CHECK_THROWS_AS(lqns::SolverLQNS<double>(m, bad), InputError);
    lqns::LqnsOptions badms;
    badms.multiserver = "linearizer";
    CHECK_THROWS_AS(lqns::SolverLQNS<double>(m, badms), InputError);
}

TEST_CASE("SolverLQNS: one customer gives the closed-form cycle time") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    const lqn::LqnModel<double> m = tandem_builder().model();
    lqns::SolverLQNS<double> s(m, lqns::LqnsOptions());
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();

    // N = 1, so the cycle is Z + 0.3 + 1.0 and nothing queues.
    const double X = 1.0 / (2.0 + 0.3 + 1.0);
    const std::size_t t1 = idx_of(sn, "T1"), e2 = idx_of(sn, "E2"), p2 = idx_of(sn, "P2");
    REQUIRE(t1);
    REQUIRE(e2);
    REQUIRE(p2);
    REQUIRE(sol.defined_T[t1]);
    CHECK(sol.TN[t1] == doctest::Approx(X).epsilon(1e-4));
    CHECK(sol.TN[e2] == doctest::Approx(X).epsilon(1e-4));
    // RespT is the phase-1 service time of the entry: its own 1.0, unqueued.
    REQUIRE(sol.defined_R[e2]);
    CHECK(sol.RN[e2] == doctest::Approx(1.0).epsilon(1e-4));
    // Util is the processor utilization: X jobs/s times 1.0 s of demand.
    REQUIRE(sol.defined_U[p2]);
    CHECK(sol.UN[p2] == doctest::Approx(X).epsilon(1e-4));
    // lqns reports neither a residence time nor an arrival rate.
    CHECK_FALSE(sol.defined_W[e2]);
    CHECK_FALSE(sol.defined_A[e2]);
    CHECK(s.raw_avg().iterations >= 1);
}

TEST_CASE("SolverLQNS: the multiserver pragma reaches a multi-server task") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    // Two customers against a two-server task. The utilization column is the
    // proc-utilization SUMMED over the host's servers, which is what SolverLN
    // reports too, so it is bounded by the multiplicity and not by one, and the
    // wrapper passes it through verbatim. Dividing it by the multiplicity was a
    // defect: see _kb/06-solver-catalog.md (LQNS wrapper) and the 2026-08-04
    // log entry, adjudicated against the raw .lqxo on test_LQN_18.
    lqn::LqnBuilder<double> b;
    b.processor("P1", kInf, SchedStrategy::INF);
    b.processor("P2", 2, SchedStrategy::FCFS);
    b.task("T1", 4, SchedStrategy::REF, "P1");
    b.think_time("T1", E(0.1));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", E(0.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");

    lqns::LqnsOptions opt;
    opt.multiserver = "rolia";
    lqns::SolverLQNS<double> s(b.model(), opt);
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();
    const std::size_t p2 = idx_of(sn, "P2");
    REQUIRE(sol.defined_U[p2]);
    CHECK(sol.UN[p2] > 1.0);
    CHECK(sol.UN[p2] <= 2.0 + 1e-9);
    // Both servers near saturation, and the column is the raw one.
    CHECK(s.raw_avg().procutil[p2] == doctest::Approx(sol.UN[p2]).epsilon(1e-12));
}

TEST_CASE("SolverLQNS: the simulator answers the same model within its noise") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    const lqn::LqnModel<double> m = tandem_builder().model();
    lqns::LqnsOptions opt;
    opt.method = "lqsim";
    opt.samples = 20000;
    CHECK(lqns::SolverLQNS<double>::is_stochastic_method(opt.method));
    lqns::SolverLQNS<double> s(m, opt);
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const std::size_t t1 = idx_of(s.get_struct(), "T1");
    const double X = 1.0 / (2.0 + 0.3 + 1.0);
    REQUIRE(sol.defined_T[t1]);
    CHECK(sol.TN[t1] == doctest::Approx(X).epsilon(0.05));
}

TEST_CASE("SolverLQNS: a processor, task and entry sharing a name get their own rows") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    // lqngen names a client processor, its task and its entry all `c0`, so
    // `lqn.names` holds the name three times. The reference resolves the .lqxo
    // by that flat list and lands one element's numbers on all three; the kind
    // is in the tag being read, so this port keeps them apart.
    // Ground truth, lqns 6.2.28 on this file: result-processor c0 utilization
    // 0.710457, result-task c0 utilization 4 throughput 0.340088.
    const std::string model = std::string(LINE_MP_REPO_ROOT) +
                              "/jar/src/test/resources/lqn/randomLQN/"
                              "model_C1_L2_T2_P2_c4_z1_s1_t1_p1_y1_4.lqnx";
    lqns::SolverLQNS<double> s(lqn::read_lqnx_model<double>(model), lqns::LqnsOptions());
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();

    std::size_t host = 0, task = 0, entry = 0;
    for (std::size_t i = 1; i <= sn.nidx; ++i) {
        if (sn.names[i] != "c0") continue;
        if (sn.type[i] == lang::LqnElement::HOST) host = i;
        if (sn.type[i] == lang::LqnElement::TASK) task = i;
        if (sn.type[i] == lang::LqnElement::ENTRY) entry = i;
    }
    REQUIRE(host);
    REQUIRE(task);
    REQUIRE(entry);
    REQUIRE(sol.defined_U[host]);
    CHECK(sol.UN[host] == doctest::Approx(0.710457).epsilon(1e-4));
    // The task row is the one a flat-name lookup loses.
    REQUIRE(sol.defined_Q[task]);
    CHECK(sol.QN[task] == doctest::Approx(4.0).epsilon(1e-4));
    REQUIRE(sol.defined_T[task]);
    CHECK(sol.TN[task] == doctest::Approx(0.340088).epsilon(1e-4));
    // A processor has no utilization of its own in the .lqxo, only a
    // proc-utilization, so QLen there stays undefined rather than borrowing one.
    CHECK_FALSE(sol.defined_Q[host]);
    REQUIRE(sol.defined_Q[entry]);
    CHECK(sol.QN[entry] == doctest::Approx(4.0).epsilon(1e-4));
}

TEST_CASE("SolverLQNS CLI: -s lqns solves the layered path and refuses the rest") {
    const std::string model = std::string(LINE_MP_REPO_ROOT) +
                              "/jar/src/test/resources/lqn/randomLQN/"
                              "model_C1_L2_T2_P2_c4_z1_s1_t1_p1_y1_4.lqnx";
    auto run = [](const std::string& args) {
        const std::string cmd = std::string(LINE_MP_CLI_BINARY) + " " + args + " 2>&1";
        std::string out;
        FILE* p = popen(cmd.c_str(), "r");
        REQUIRE(p != nullptr);
        char buf[512];
        while (std::fgets(buf, sizeof(buf), p) != nullptr) out += buf;
        pclose(p);
        return out;
    };

    // A layer engine is not a knob of a solver that solves no layers.
    const std::string clash = run("-f " + model + " -s lqns --layer-solver mva");
    CHECK(clash.find("solves no layers") != std::string::npos);
    // Neither is an analysis the binary does not compute.
    const std::string tran = run("-f " + model + " -s lqns -a tran --tspan 10");
    CHECK(tran.find("-a avg") != std::string::npos);
    // And the wrapper's own flags do not leak onto the in-process solvers.
    const std::string leak = run("-f " + model + " -s ln --keep");
    CHECK(leak.find("-s lqns only") != std::string::npos);

    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the end-to-end CLI solve");
        return;
    }
    const std::string out = run("-f " + model + " -s lqns");
    CHECK(out.find("SolverLQNS") != std::string::npos);
    CHECK(out.find("Version 6") != std::string::npos);
    CHECK(out.find("NaN") != std::string::npos);  // ResidT and ArvR, per element
}

TEST_CASE("SolverLQNS: an entry's Util sums ITS OWN activities") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    // lqns credits host work to whichever level declares the host demand, and in
    // the activity-graph form -- the only form lqn_writer.h emits -- an entry
    // declares none, so every `result-entry` carries a literal
    // proc-utilization="0" and the work sits on the `result-activity` rows.
    // T2 hosts TWO entries here, which is what pins the rule down: the entry
    // sums are NOT the task's proc-utilization, so a fix that copied the task
    // value, or that summed every activity of the task into each of its
    // entries, would pass a one-entry task and fail this.
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 100, SchedStrategy::REF, "P1");
    b.think_time("T1", E(10.0));
    b.entry("E1", "T1");
    b.activity("A1", E(1.0), "T1");
    b.bound_to("A1", "E1");

    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E2", "T2");
    b.entry("E3", "T2");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 1.0);

    b.activity("A20", E(1.0), "T2");
    b.bound_to("A20", "E2");
    b.activity("A21", E(1.0), "T2");
    b.activity("A22", E(1.0), "T2");
    b.serial("A20", "A21");
    b.serial("A21", "A22");
    b.replies_to("A22", "E2");

    b.activity("A3", E(1.0), "T2");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");

    lqns::SolverLQNS<double> s(b.model(), lqns::LqnsOptions());
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();

    const std::size_t t2 = idx_of(sn, "T2"), e1 = idx_of(sn, "E1");
    const std::size_t e2 = idx_of(sn, "E2"), e3 = idx_of(sn, "E3");
    const std::size_t a1 = idx_of(sn, "A1"), a3 = idx_of(sn, "A3");
    const std::size_t a20 = idx_of(sn, "A20"), a21 = idx_of(sn, "A21");
    const std::size_t a22 = idx_of(sn, "A22");
    REQUIRE(t2);
    REQUIRE(e2);
    REQUIRE(e3);
    REQUIRE(a22);

    REQUIRE(sol.defined_U[e2]);
    REQUIRE(sol.defined_U[e3]);
    const double chain = sol.UN[a20] + sol.UN[a21] + sol.UN[a22];
    CHECK(sol.UN[e2] == doctest::Approx(chain).epsilon(1e-6));
    CHECK(sol.UN[e3] == doctest::Approx(sol.UN[a3]).epsilon(1e-6));
    CHECK(sol.UN[e1] == doctest::Approx(sol.UN[a1]).epsilon(1e-6));

    // Three unit-demand activities against one, at the same throughput.
    CHECK(sol.UN[e2] == doctest::Approx(3.0 * sol.UN[e3]).epsilon(1e-3));
    // Together they are the task, so neither entry can be carrying its value.
    CHECK(sol.UN[e2] + sol.UN[e3] == doctest::Approx(sol.UN[t2]).epsilon(1e-6));
    CHECK(sol.UN[e2] > 1.1 * sol.UN[e3]);

    // The regression itself: read verbatim, every entry comes back 0.
    CHECK(sol.UN[e1] > 0.0);
    CHECK(sol.UN[e2] > 0.0);
    CHECK(sol.UN[e3] > 0.0);
}

/**
 * An entry lqns never invoked has a service time of ZERO, not an absent one.
 *
 * lqns omits `phase1-service-time` from `result-entry` exactly when the entry's
 * throughput is zero -- nothing was served, so there is no per-invocation mean
 * to report. Read verbatim that is a NaN, and it landed in the RespT column
 * where the LQN table says an entry HAS a response time and every other solver
 * reports one. The NaN mask is part of the answer (see
 * `_kb/06-solver-catalog.md`), and a tolerance-based parity comparison cannot
 * see a break in it: `compare_values` passes any cell where either side is NaN.
 *
 * The value is derived from the activity rows and ONLY where they are unanimous:
 * if every activity reachable from the entry reports a zero service time then
 * every aggregation law agrees on zero. It is deliberately not generalised the
 * way the proc-utilization case above is; the second case below pins that down.
 *
 * Twins: `python/tests/test_lqns_entry_svct.py`,
 * `jar/src/test/java/jline/solvers/wrappers/lqns/SolverLQNSEntryServiceTimeTest.java`,
 * `line-test.git/test/testsMisc/test_lqns_entry_svct.m`.
 */
TEST_CASE("SolverLQNS: an entry lqns never invoked reports zero, not NaN") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    // A working two-tier model plus a component no reference task reaches: T3 is
    // not a reference task and nobody calls E3, so lqns solves it at zero
    // throughput and writes a result-entry with no phase1-service-time. This is
    // the shape lqn_ofbiz carries in its USAGE_DELAY component. E3 runs TWO
    // activities in series so the writer emits the activity-graph form; a single
    // bound activity goes out as entry-phase-activities, where the reader has a
    // fallback of its own and the omission never surfaces.
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", E(1.0));
    b.entry("E1", "T1");
    b.activity("A1", E(1.0), "T1");
    b.bound_to("A1", "E1");

    b.processor("P2", 1, SchedStrategy::INF);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E2", "T2");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");

    b.processor("P3", 1, SchedStrategy::INF);
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E3", "T3");
    b.activity("A3", E(1.0), "T3");
    b.bound_to("A3", "E3");
    b.activity("A3b", E(1.0), "T3");
    b.serial("A3", "A3b");
    b.replies_to("A3b", "E3");

    lqns::SolverLQNS<double> s(b.model(), lqns::LqnsOptions());
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();

    const std::size_t e1 = idx_of(sn, "E1"), e2 = idx_of(sn, "E2");
    const std::size_t e3 = idx_of(sn, "E3"), a3 = idx_of(sn, "A3");
    const std::size_t a3b = idx_of(sn, "A3b");
    REQUIRE(e3);
    REQUIRE(a3b);

    // the unreachable component solves at zero throughput
    CHECK(sol.TN[e3] == doctest::Approx(0.0));
    CHECK(sol.TN[a3] == doctest::Approx(0.0));
    CHECK(sol.TN[a3b] == doctest::Approx(0.0));

    // the regression itself: RespT was undefined, because lqns omits the attribute
    CHECK(sol.defined_R[e3]);
    CHECK(sol.RN[e3] == doctest::Approx(0.0));

    // the reachable entries are untouched and still carry lqns' own numbers
    CHECK(sol.defined_R[e1]);
    CHECK(sol.defined_R[e2]);
    CHECK(sol.RN[e1] == doctest::Approx(2.0).epsilon(2e-3));
    CHECK(sol.RN[e2] == doctest::Approx(1.0).epsilon(1e-3));
}

/**
 * Where lqns DOES report the attribute, its value survives verbatim.
 *
 * Utilizations add over an activity graph and response times do not, so the
 * fallback must never become a sum: on this OrFork the sum over the entry's
 * activities exceeds what lqns reports, and the reported value is the one that
 * must come back.
 */
TEST_CASE("SolverLQNS: the service-time fallback does not touch a branching entry") {
    if (!lqns_here()) {
        MESSAGE("lqns is not installed: skipping the cases that run it");
        return;
    }
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 10, SchedStrategy::REF, "P1");
    b.think_time("T1", E(0.1));
    b.entry("E1", "T1");
    b.activity("A1", E(1.0), "T1");
    b.bound_to("A1", "E1");

    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E2", "T2");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A20", E(1.0), "T2");
    b.bound_to("A20", "E2");
    b.activity("A21", E(1.0), "T2");
    b.activity("A22", E(1.0), "T2");
    b.or_fork("A20", std::vector<std::string>{"A21", "A22"}, std::vector<double>{0.5, 0.5});
    b.replies_to("A21", "E2");
    b.replies_to("A22", "E2");

    lqns::SolverLQNS<double> s(b.model(), lqns::LqnsOptions());
    const lqns::LqnsSolution<double> sol = s.get_ensemble_avg();
    const lqn::LqnStruct<double>& sn = s.get_struct();

    const std::size_t e2 = idx_of(sn, "E2"), a20 = idx_of(sn, "A20");
    const std::size_t a21 = idx_of(sn, "A21"), a22 = idx_of(sn, "A22");
    REQUIRE(e2);
    REQUIRE(a22);
    REQUIRE(sol.defined_R[e2]);

    const double summed = sol.RN[a20] + sol.RN[a21] + sol.RN[a22];
    // one arm of the fork is not taken, so the sum over-counts
    INFO("sum over the entry's activities " << summed << " against reported " << sol.RN[e2]);
    CHECK(summed > sol.RN[e2] + 1e-6);
}
