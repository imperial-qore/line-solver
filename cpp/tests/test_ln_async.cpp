/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Asynchronous calls: an open chain out of the layer's own Source.
 *
 * Model mirrors line-test.git/test/testsAdvFeatures/lqn/test_asynchCall.m,
 * Test 1 ("Simple Pure Async Model"): a reference task with 5 clients and a
 * 50s think time fires one async call per firing at AsyncEntry, whose task is
 * FCFS on its own PS processor. The caller does not block, so the two tasks
 * couple in ONE direction only -- through the caller's throughput, which
 * becomes the arrival rate of the async stream.
 *
 * Reference numbers are MATLAB SolverLN.getAvgTable() on the identical model
 * (defaultOptions, 2026-07-29).
 */

#include <limits>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

lqn::LqnStruct<double> build_pure_async() {
    lqn::LqnBuilder<double> b;
    b.processor("ClientProc", 1, SchedStrategy::PS);
    b.processor("AsyncServerProc", 1, SchedStrategy::PS);

    b.task("Client", 5, SchedStrategy::REF, "ClientProc");
    b.think_time("Client", E(50.0));
    b.task("AsyncServer", 1, SchedStrategy::FCFS, "AsyncServerProc");

    b.entry("ClientEntry", "Client");
    b.entry("AsyncEntry", "AsyncServer");

    b.activity("ClientAct", E(0.5), "Client");
    b.bound_to("ClientAct", "ClientEntry");
    b.async_call("ClientAct", "AsyncEntry", 1.0);

    b.activity("AsyncAct", E(1.0), "AsyncServer");
    b.bound_to("AsyncAct", "AsyncEntry");

    return b.build();
}

/**
 * One task reached BOTH ways: SE synchronously, AE asynchronously.
 *
 * An entry may not be called both ways at once (getStruct.m rejects it), so
 * the two land on separate entries of the same task. Its layer then carries
 * closed classes and an open one side by side, which is the case that would
 * break refresh_chains if the two ever shared a chain: they carry different
 * reference stations (the client Delay against the Source).
 */
lqn::LqnStruct<double> build_mixed() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);

    b.task("Client", 2, SchedStrategy::REF, "P1");
    b.think_time("Client", E(5.0));
    b.task("Server", 1, SchedStrategy::FCFS, "P2");

    b.entry("CE", "Client");
    b.entry("SE", "Server");
    b.entry("AE", "Server");

    b.activity("CA", E(0.5), "Client");
    b.bound_to("CA", "CE");
    b.sync_call("CA", "SE", 1.0);
    b.async_call("CA", "AE", 1.0);

    b.activity("SA", E(1.0), "Server");
    b.bound_to("SA", "SE");
    b.replies_to("SA", "SE");

    b.activity("AA", E(0.8), "Server");
    b.bound_to("AA", "AE");

    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("async: closed and open classes coexist in one layer, matching MATLAB") {
    const lqn::LqnStruct<double> l = build_mixed();
    ln::LnOptions opt;
    // one class per task, entry, activity and call is the ROUTING encoding; the
    // 'srvn' alias resolves to 'srvn.ph' here, whose layers carry one class per
    // caller task instead, so the counts below are read under the encoding they
    // describe
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);

    // MATLAB: T:Server has 5 classes and 3 stations (the Source is the third).
    // P:P2 has 5 classes but only 2 stations: an async class is declared in the
    // layer whose SERVER owns the called entry, and P2 is a host, not the task.
    const qn::Layer<double>* tserver = nullptr;
    for (const qn::Layer<double>& L : s.layers()) {
        if (L.name == "T:Server") tserver = &L;
        if (L.name == "P:P2") {
            CHECK(L.nstations == 2);
            CHECK(L.sourceIdx == 0);
        }
    }
    REQUIRE(tserver != nullptr);
    CHECK(tserver->nstations == 3);
    CHECK(tserver->nclasses == 5);
    CHECK(tserver->sourceIdx != 0);
    // the mix itself: at least one finite-population class and one infinite
    bool any_closed = false, any_open = false;
    for (const qn::JobClass& c : tserver->classes)
        (std::isinf(c.population) ? any_open : any_closed) = true;
    CHECK(any_closed);
    CHECK(any_open);

    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:P1", 0.14109, 0.0, 0.0},
        {"P:P2", 0.50791, 0.0, 0.0},
        {"R:Client", 0.14109, 0.28217, 0.0},
        {"T:Server", 0.50791, 0.56434, 0.0},
        {"E:CE", 0.14109, 0.28217, 2.0879},
        {"E:SE", 0.28217, 0.28217, 1.0},
        {"E:AE", 0.22574, 0.28217, 0.8},
        {"A:CA", 0.14109, 0.28217, 2.0879},
        {"A:SA", 0.28217, 0.28217, 1.0},
        {"A:AA", 0.22574, 0.28217, 0.8},
    };
    for (const Row& row : rows) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(row.hn);
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
        if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
    }
}

TEST_CASE("async: the target task's layer carries a Source and one open class") {
    const lqn::LqnStruct<double> l = build_pure_async();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);

    // MATLAB reports 3 layers, the third being T:AsyncServer with 3 stations
    // (Clients, the server, the Source) and exactly ONE class: the async call.
    // The caller blocks on nothing, so it contributes no closed class here.
    CHECK(s.nlayers() == 3);
    const qn::Layer<double>* async_layer = nullptr;
    for (const qn::Layer<double>& L : s.layers())
        if (L.name == "T:AsyncServer") async_layer = &L;
    REQUIRE(async_layer != nullptr);
    CHECK(async_layer->nstations == 3);
    CHECK(async_layer->nclasses == 1);
    CHECK(async_layer->sourceIdx != 0);
    CHECK(async_layer->sinkNode != 0);
    CHECK(async_layer->has_open_classes());
    CHECK(std::isinf(async_layer->classes[0].population));
    CHECK(async_layer->classes[0].refstat == async_layer->sourceIdx);
}

TEST_CASE("async: a pure async model matches MATLAB SolverLN") {
    const lqn::LqnStruct<double> l = build_pure_async();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // MATLAB AvgTable:
    //   ClientProc       Util 0.049485
    //   AsyncServerProc  Util 0.09897
    //   Client   RefTask QLen 0.051504 Util 0.049485 ResidT 0.5204 Tput 0.09897
    //   AsyncServer Task QLen 0.09897  Util 0.09897  ResidT 1      Tput 0.09897
    //   ClientEntry Entry RespT 0.5204 Tput 0.09897
    //   AsyncEntry  Entry RespT 1      Tput 0.09897
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:ClientProc", 0.049485, 0.0, 0.0},
        {"P:AsyncServerProc", 0.09897, 0.0, 0.0},
        {"R:Client", 0.049485, 0.09897, 0.0},
        {"T:AsyncServer", 0.09897, 0.09897, 0.0},
        {"E:ClientEntry", 0.049485, 0.09897, 0.5204},
        {"E:AsyncEntry", 0.09897, 0.09897, 1.0},
        {"A:ClientAct", 0.049485, 0.09897, 0.5204},
        {"A:AsyncAct", 0.09897, 0.09897, 1.0},
    };
    for (const Row& row : rows) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(row.hn);
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
        if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
    }

    // The async server is driven entirely by the caller's firing rate, so the
    // two throughputs must agree: every firing releases exactly one job.
    CHECK(sol.TN[idx_of(l, "E:AsyncEntry")] ==
          doctest::Approx(sol.TN[idx_of(l, "E:ClientEntry")]).epsilon(1e-6));
    // ResidT of an async-only entry is its own service, NOT a caller-side
    // residence: nobody waits on it (updateMetricsDefault.m:63-78).
    CHECK(sol.RN[idx_of(l, "E:AsyncEntry")] == doctest::Approx(1.0).epsilon(1e-3));
}
