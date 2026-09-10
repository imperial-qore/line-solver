/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * An open arrival at an entry: an exogenous stream, not a call.
 *
 * A closed client calls SE synchronously while SE ALSO receives its own
 * Poisson arrivals. The arrival is not a visit of any caller, so it enters the
 * layer of the entry's own host through a Source and leaves at the Sink after
 * one pass -- it does not walk the activity graph the way a call-fed visit
 * does (buildLayersRecursive.m:474-495).
 *
 * Unlike an async call's stream, this one is EXOGENOUS: its rate is the
 * user's distribution and never moves with the fixed point, so it carries no
 * row in arv_call_map.
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

lqn::LqnStruct<double> build_open_arrival() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);

    b.task("Client", 2, SchedStrategy::REF, "P1");
    b.think_time("Client", E(1.0));
    b.task("Server", 1, SchedStrategy::FCFS, "P2");

    b.entry("CE", "Client");
    b.entry("SE", "Server");
    b.open_arrival("SE", E(10.0));

    b.activity("CA", E(0.5), "Client");
    b.bound_to("CA", "CE");
    b.sync_call("CA", "SE", 1.0);

    b.activity("SA", E(1.0), "Server");
    b.bound_to("SA", "SE");
    b.replies_to("SA", "SE");

    return b.build();
}

/**
 * The PURE-OPEN model: nothing calls T1, its single entry only receives a stream.
 * With no caller there is no task layer, so the caller chain on the host layer would
 * never be given a surrogate delay; carrying the arrival as an open class ON TOP of
 * that unthrottled chain loaded the processor twice (0.68 against 0.32). The stream
 * is represented by the thread pool instead -- see open_arrival_rate_of.
 */
lqn::LqnStruct<double> build_pure_open_arrival() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 1, SchedStrategy::FCFS, "P1");
    b.think_time("T1", D::immediate());
    b.entry("E1", "T1");
    b.open_arrival("E1", E(5.0));  // rate 0.2
    b.activity("A1", E(1.6), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("open arrival: only the entry's host layer gains a Source") {
    const lqn::LqnStruct<double> l = build_open_arrival();
    ln::LnOptions opt;
    // the ROUTING encoding, which is what carries an open class in a layer at
    // all; the 'srvn' alias now resolves to 'srvn.ph' on this model, where the
    // stream is a term of the composed law instead. See the case below.
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);

    // MATLAB: 3 layers; P:P2 has 3 stations (Clients, P2, Source) and 4
    // classes -- one more than the other layers, the open stream. The entry
    // belongs to Server, whose host is P2, so no other layer sees the arrival.
    CHECK(s.nlayers() == 3);
    for (const qn::Layer<double>& L : s.layers()) {
        CAPTURE(L.name);
        if (L.name == "P:P2") {
            CHECK(L.nstations == 3);
            CHECK(L.nclasses == 4);
            CHECK(L.sourceIdx != 0);
            CHECK(L.sinkNode != 0);
            CHECK(L.has_open_classes());
        } else {
            CHECK(L.sourceIdx == 0);
            CHECK(!L.has_open_classes());
        }
    }
}

TEST_CASE("open arrival: the model matches MATLAB SolverLN") {
    const lqn::LqnStruct<double> l = build_open_arrival();
    ln::LnOptions opt;
    // the encoding this golden was recorded under; the two encodings reach
    // DIFFERENT fixed points on this model, and both are asserted here
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // MATLAB AvgTable:
    //   P1 Util 0.32426, P2 Util 0.64852
    //   Client RefTask QLen 1.3515  Util 0.32426 ResidT 0.5    Tput 0.64852
    //   Server Task    QLen 0.72058 Util 0.64852 ResidT 1.1111 Tput 0.64852
    //   CE Entry RespT 2.0839 Tput 0.64852
    //   SE Entry RespT 1.1111 Tput 0.64852
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:P1", 0.32426, 0.0, 0.0},
        {"P:P2", 0.64852, 0.0, 0.0},
        {"R:Client", 0.32426, 0.64852, 0.0},
        {"T:Server", 0.64852, 0.64852, 0.0},
        {"E:CE", 0.32426, 0.64852, 2.0839},
        {"E:SE", 0.64852, 0.64852, 1.1111},
        {"A:CA", 0.32426, 0.64852, 2.0839},
        {"A:SA", 0.64852, 0.64852, 1.1111},
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

TEST_CASE("open arrival: the default encoding matches MATLAB srvn.ph") {
    const lqn::LqnStruct<double> l = build_open_arrival();
    ln::LnOptions opt;  // the 'srvn' alias, which resolves to 'srvn.ph' here
    ln::SolverLN<double> s(l, opt);
    CHECK(s.state_lnmethod() == "srvn.ph");
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // MATLAB AvgTable with method='srvn.ph', recorded 2026-08-11. The composed
    // entry law charges the stream at the server's own law rather than through a
    // class of its own, so SE answers in its bare 1 and the layer runs hotter
    // than the routing encoding above: a DIFFERENT fixed point, not a rounding.
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:P1", 0.34483, 0.0, 0.0},
        {"P:P2", 0.76484, 0.0, 0.0},
        {"R:Client", 0.34483, 0.68966, 0.0},
        {"T:Server", 0.76484, 0.76484, 0.0},
        {"E:CE", 0.34483, 0.68966, 1.9},
        {"E:SE", 0.76484, 0.76484, 1.0},
        {"A:CA", 0.34483, 0.68966, 1.9},
        {"A:SA", 0.76484, 0.76484, 1.0},
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

TEST_CASE("pure open arrival: the stream drives the thread pool, not a class") {
    const lqn::LqnStruct<double> l = build_pure_open_arrival();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);

    // Nothing calls T1, so the arrival is carried by its caller chain and the layer
    // needs no Source at all.
    for (const qn::Layer<double>& L : s.layers()) {
        CAPTURE(L.name);
        CHECK(L.sourceIdx == 0);
        CHECK(!L.has_open_classes());
    }

    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // lqns 6.2.28 gives exactly 0.2 / 0.32 / 1.6; lqsim 0.1996 / 0.319 / 1.580 and
    // LDES 0.19986 / 0.31957 / 1.599 confirm it. MATLAB, the JAR and native Python all
    // report the same after the 2026-08-11 fix.
    struct PRow { const char* hn; double util; double tput; double respt; };
    const PRow prows[] = {
        {"P:P1", 0.32, 0.0, 0.0},
        {"T:T1", 0.32, 0.2, 0.0},
        {"E:E1", 0.32, 0.2, 1.6},
        {"A:A1", 0.32, 0.2, 1.6},
    };
    for (const PRow& row : prows) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(row.hn);
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
        if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
    }
}
