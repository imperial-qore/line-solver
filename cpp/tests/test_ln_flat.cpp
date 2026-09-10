/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The squashed (`flat`) layering: one submodel holding every server.
 *
 * Under `srvn` each server gets a submodel of its own and sees the others
 * through a surrogate delay; under `flat` they are stations of ONE network, so
 * two tasks sharing a processor contend inside a single AMVA solve instead of
 * through the outer fixed point. The two therefore answer differently, and the
 * point of this file is that C++ answers what MATLAB answers, not merely
 * something plausible.
 *
 * The model is the smallest one on which the distinction bites: T2 and T3 share
 * processor P2, so squashing changes which contention the layer solver sees.
 * The reference values are MATLAB SolverLN(SolverMVA) with
 * options.config.layering='flat' on exactly this model.
 *
 * A second case pins the refusals down: squashing is REFUSED, by name, for a
 * replicated element, a cache task or a setup task, each of which carries state
 * that only a submodel of its own can hold. A silent fall-back to `srvn` there
 * would answer a question nobody asked.
 */

#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using namespace line::lang;
using D = Distrib<double>;

namespace {

/** T1 calls E2 once and E3 twice; T2 and T3 share processor P2. */
lqn::LqnStruct<double> build_shared() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", D::exp_mean(1.0));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.task("T3", 2, SchedStrategy::FCFS, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", D::exp_mean(0.2), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 2.0);
    b.activity("A2", D::exp_mean(0.4), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", D::exp_mean(0.3), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("flat: one layer holds every server, and it is the MATLAB answer") {
    const lqn::LqnStruct<double> l = build_shared();
    ln::LnOptions opt;
    opt.method = "flat";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    // 'flat' is an alias, and unconditional: it does not probe flat.ph
    CHECK(s.state_lnmethod() == "flat.cs");
    // ONE submodel, not one per server
    CHECK(s.nlayers() == 1);
    // and it carries a station for each of the two processors and two called
    // tasks (T1 is a reference task, so it is a caller and not a server)
    const qn::Layer<double>& L = s.layers()[0];
    CHECK(L.flat);
    CHECK(L.host_stations.size() == 2);
    CHECK(L.task_stations.size() == 2);
    // every served element resolves to a station of its own, none to the
    // scalar fallback shared by all of them
    CHECK(L.server_idx_of[idx_of(l, "P:P1")] != L.server_idx_of[idx_of(l, "P:P2")]);
    CHECK(L.server_idx_of[idx_of(l, "T:T2")] != L.server_idx_of[idx_of(l, "T:T3")]);

    // MATLAB SolverLN(SolverMVA), options.config.layering='flat'
    struct Row {
        const char* name;
        double qlen, util, respt, tput;
    };
    const Row rows[] = {
        {"R:T1", 3.8146656, 0.20592732, 0.0, 1.0296366},
        {"T:T2", 1.1370420, 0.39487322, 0.0, 0.98718304},
        {"T:T3", 1.5210085, 0.58568036, 0.0, 1.9522679},
        {"E:E1", 3.8146656, 0.20592732, 3.7048659, 1.0296366},
        {"E:E2", 1.1370420, 0.39487322, 1.1518046, 0.98718304},
        {"E:E3", 1.5210085, 0.58568036, 0.77909823, 1.9522679},
        {"A:A1", 3.9703633, 0.20592732, 3.8560821, 1.0296366},
        {"A:A2", 1.2195052, 0.39487322, 1.2353385, 0.98718304},
        {"A:A3", 1.6249338, 0.58568036, 0.83233137, 1.9522679},
    };
    for (const Row& r : rows) {
        const std::size_t i = idx_of(l, r.name);
        REQUIRE(i != 0);
        CHECK(sol.QN[i] == doctest::Approx(r.qlen).epsilon(1e-6));
        CHECK(sol.UN[i] == doctest::Approx(r.util).epsilon(1e-6));
        CHECK(sol.TN[i] == doctest::Approx(r.tput).epsilon(1e-6));
        if (r.respt > 0.0) CHECK(sol.RN[i] == doctest::Approx(r.respt).epsilon(1e-5));
    }
    // the processors report utilization only, as under srvn
    CHECK(sol.UN[idx_of(l, "P:P1")] == doctest::Approx(0.20592732).epsilon(1e-6));
    CHECK(sol.UN[idx_of(l, "P:P2")] == doctest::Approx(0.98055358).epsilon(1e-6));
}

TEST_CASE("flat and srvn are different answers, not the same one twice") {
    const lqn::LqnStruct<double> l = build_shared();
    ln::LnOptions fo;
    fo.method = "flat";
    ln::LnOptions so;
    so.method = "srvn.cs";
    const ln::LnSolution<double> f = ln::SolverLN<double>(l, fo).get_ensemble_avg();
    const ln::LnSolution<double> r = ln::SolverLN<double>(l, so).get_ensemble_avg();
    // P2 is the shared processor: squashing resolves its contention inside one
    // solve and reads it 2.5% lower than the decomposition does
    CHECK(f.UN[idx_of(l, "P:P2")] < r.UN[idx_of(l, "P:P2")] - 0.01);
    CHECK(r.UN[idx_of(l, "P:P2")] == doctest::Approx(1.0058462).epsilon(1e-6));
}

TEST_CASE("flat refuses what a single submodel cannot carry") {
    // a replicated task: the routing addresses it once, the copies need a
    // station each, and only a submodel of its own can hold both readings
    {
        lqn::LqnStruct<double> l = build_shared();
        l.repl[idx_of(l, "T:T2")] = 2.0;
        ln::LnOptions opt;
        opt.method = "flat";
        CHECK_THROWS_AS(ln::SolverLN<double>(l, opt).get_ensemble_avg(), UnsupportedError);
    }
    // a setup task: the delay-off belongs to the station that powers down
    {
        lqn::LqnStruct<double> l = build_shared();
        l.hassetup[idx_of(l, "T:T2")] = true;
        ln::LnOptions opt;
        opt.method = "flat";
        CHECK_THROWS_AS(ln::SolverLN<double>(l, opt).get_ensemble_avg(), UnsupportedError);
    }
    // a cache task: the Cache node lives in the host layer its reads queue at
    {
        lqn::LqnStruct<double> l = build_shared();
        l.iscache[idx_of(l, "T:T2")] = true;
        ln::LnOptions opt;
        opt.method = "flat";
        CHECK_THROWS_AS(ln::SolverLN<double>(l, opt).get_ensemble_avg(), UnsupportedError);
    }
}

TEST_CASE("the method name carries the layering and the encoding") {
    const lqn::LqnStruct<double> l = build_shared();
    struct Case {
        const char* asked;
        const char* got;
    };
    // 'flat' and 'squashed' are aliases of flat.cs; 'srvn' probes srvn.ph and
    // takes it here, this model being series-parallel throughout
    const Case cases[] = {{"flat", "flat.cs"},   {"flat.cs", "flat.cs"}, {"squashed", "flat.cs"},
                          {"srvn.cs", "srvn.cs"}, {"srvn", "srvn.ph"},   {"default", "srvn.ph"}};
    for (const Case& c : cases) {
        ln::LnOptions opt;
        opt.method = c.asked;
        ln::SolverLN<double> s(l, opt);
        s.get_ensemble_avg();
        CHECK(s.state_lnmethod() == std::string(c.got));
    }
}
