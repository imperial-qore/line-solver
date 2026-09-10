/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * A setup task: a server that powers down when idle and pays to restart.
 *
 * T2 has a 0.5 setup and a 1.0 delay-off timer. Its own service is 1.0, so the
 * setup is the whole difference: without it E2 answers in exactly 1.0, with it
 * in 1.2227 under the default encoding and 1.2147 under 'srvn.cs'. That gap is
 * what the test is really pinning -- a wiring that dropped the setup would
 * still converge, and would still look plausible.
 *
 * The restart is charged to the ENTRY, not wired onto the layer station: no
 * layer carries a SetupDelayOffParam and none is diverted to SolverMAM, so the
 * layer solver the user asked for serves it.
 *
 * Reference numbers are MATLAB SolverLN.getAvgTable() on the identical model,
 * RE-RECORDED 2026-08-12, one table per srvn encoding. The charge is the
 * analytical one of lqn_setup_charge (Gandhi, Harchol-Balter and Adan,
 * Performance Evaluation 67(11), 2010), which moved E2 from the earlier 1.724
 * to its present value. Both encodings call that same closed form, but they do
 * not converge to the same fixed point -- the alias resolves to 'srvn.ph',
 * whose composed entry law is a different layer from the routing one -- so the
 * table is asserted per encoding rather than once under the alias.
 */

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

lqn::LqnBuilder<double> setup_builder(bool with_setup) {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);

    b.task("T1", 2, SchedStrategy::REF, "P1");
    b.think_time("T1", E(2.0));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    if (with_setup) b.setup_time("T2", E(0.5), E(1.0));

    b.entry("E1", "T1");
    b.entry("E2", "T2");

    b.activity("A1", E(0.3), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);

    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");

    return b;
}

lqn::LqnStruct<double> build_setup_model(bool with_setup) {
    return setup_builder(with_setup).build();
}

std::string temp_path(const char* stem) {
    const char* base = std::getenv("TMPDIR");
    return (base && *base ? std::string(base) : std::string("/tmp")) + "/" + stem;
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("setup task: the setup reaches the struct") {
    const lqn::LqnStruct<double> l = build_setup_model(true);
    // a setup task keeps the plain T: hashname of any task
    const std::size_t t2 = idx_of(l, "T:T2");
    REQUIRE(t2 > 0);
    CHECK(l.hassetup[t2]);
    CHECK(l.setuptime[t2].mean == doctest::Approx(0.5));
    CHECK(l.delayofftime[t2].mean == doctest::Approx(1.0));
    // without one it is a plain task again
    const lqn::LqnStruct<double> plain = build_setup_model(false);
    CHECK(idx_of(plain, "T:T2") > 0);
    CHECK(!plain.hassetup[idx_of(plain, "T:T2")]);
}

TEST_CASE("setup task: the .lqnx dialect round-trips the cold start") {
    // THE ONLY TEST THAT OPENS A DOCUMENT. Every other case here builds through
    // LqnBuilder, which fills setuptime directly, so all of them passed while
    // the reader ignored <setup> and <delay-off> outright: the CPP parity row,
    // which drives line-cli on the .lqnx an example exports, was the sole place
    // the gap was visible, and it read as a 3x throughput defect on lqn_setup.
    // <setup>/<delay-off> are a LINE extension of the schema that MATLAB, the
    // JAR and Python all write and read (writeXML.m:170-183); the writer used
    // to refuse a setup task by name instead of emitting them.
    const std::string path = temp_path("line_setup_roundtrip.lqnx");
    lqn::write_lqnx(setup_builder(true).model(), path, "setup");

    const lqn::LqnStruct<double> after = lqn::read_lqnx<double>(path);
    const std::size_t t2 = idx_of(after, "T:T2");
    REQUIRE(t2 > 0);
    CHECK(after.hassetup[t2]);
    CHECK(after.setuptime[t2].mean == doctest::Approx(0.5));
    CHECK(after.delayofftime[t2].mean == doctest::Approx(1.0));

    // and the charge it buys survives with it, to the digit the in-memory
    // model answers: a document that lost the elements would come back with
    // E2 at its bare 1.0 host demand.
    ln::SolverLN<double> from_doc(after, ln::LnOptions());
    ln::SolverLN<double> in_mem(build_setup_model(true), ln::LnOptions());
    const ln::LnSolution<double> a = from_doc.get_ensemble_avg();
    const ln::LnSolution<double> b = in_mem.get_ensemble_avg();
    const std::size_t e2a = idx_of(after, "E:E2");
    const std::size_t e2b = idx_of(build_setup_model(true), "E:E2");
    REQUIRE(e2a > 0);
    REQUIRE(e2b > 0);
    CHECK(a.RN[e2a] == doctest::Approx(b.RN[e2b]).epsilon(1e-9));
    CHECK(a.RN[e2a] > 1.0);

    // a task with no setup reads back as a plain one, not as a zero-mean setup
    lqn::write_lqnx(setup_builder(false).model(), path, "plain");
    const lqn::LqnStruct<double> plain = lqn::read_lqnx<double>(path);
    CHECK(!plain.hassetup[idx_of(plain, "T:T2")]);
    std::remove(path.c_str());
}

TEST_CASE("setup task: no layer station carries the setup") {
    // A SETUP NO LONGER CHANGES HOW A LAYER IS BUILT. Wiring it onto the server
    // station routed the layer through the open M/G/1-with-setup QBD, which reads
    // the idle period off the Poisson rate 1/X and powers the thread down far
    // more often than a CLOSED layer does, and it charged the restart to the
    // ACTIVITY, where it is not host demand. The cold start is charged to the
    // ENTRY instead, with the probability that the thread was really found down
    // (setup_charge), so EVERY layer is an ordinary one under BOTH encodings and
    // the user's own layer solver serves it -- no SolverMAM/dec.poisson arm any
    // more. Twin of buildLayersRecursive.m and SolverLN.m, 2026-08-11.
    const lqn::LqnStruct<double> l = build_setup_model(true);
    for (const char* mth : {"srvn.cs", "srvn.ph"}) {
        CAPTURE(mth);
        ln::LnOptions opt;
        opt.method = mth;
        ln::SolverLN<double> s(l, opt);
        for (const qn::Layer<double>& L : s.layers()) {
            CAPTURE(L.name);
            CHECK(L.setupparam.empty());
        }
    }
    // the routing encoding still builds the three layers it always did
    ln::LnOptions cs;
    cs.method = "srvn.cs";
    ln::SolverLN<double> scs(l, cs);
    CHECK(scs.nlayers() == 3);
}

TEST_CASE("setup task: the model matches MATLAB SolverLN") {
    // ONE TABLE PER ENCODING. The two charge the restart from the SAME closed
    // form -- lqn_setup_charge, one on the layer station and one inside the
    // composed entry law -- but they do not reach the same fixed point, because
    // 'srvn.ph' composes the activity graph into a phase-type entry law and
    // 'srvn.cs' encodes it as routing. E2 answers in 1.2227 under the first and
    // 1.2147 under the second, and the whole table moves with it by the same
    // 2.9%. Both are MATLAB SolverLN AvgTable on the identical model, RE-RECORDED
    // 2026-08-12 with options.method named explicitly; the 'srvn.cs' column is
    // the table this case carried while the alias still resolved there.
    const lqn::LqnStruct<double> l = build_setup_model(true);
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row ph_rows[] = {
        {"P:P1", 0.14810452276, 0.0, 0.0},
        {"P:P2", 0.493681750993, 0.0, 0.0},
        {"R:T1", 0.14810452276, 0.493681742533, 0.0},
        {"T:T2", 0.493681750993, 0.493681750993, 0.0},
        {"E:E1", 0.14810452276, 0.493681742533, 2.05119295034},
        {"E:E2", 0.493681750993, 0.493681750993, 1.22267268771},
        {"A:A1", 0.14810452276, 0.493681742533, 2.05119295034},
        {"A:A2", 0.493681750993, 0.493681750993, 1.0},
    };
    const Row cs_rows[] = {
        {"P:P1", 0.152496745675, 0.0, 0.0},
        {"P:P2", 0.508322481709, 0.0, 0.0},
        {"R:T1", 0.152496745675, 0.508322485584, 0.0},
        {"T:T2", 0.508322481709, 0.508322481709, 0.0},
        {"E:E1", 0.152496745675, 0.508322485584, 1.9345101813},
        {"E:E2", 0.508322481709, 0.508322481709, 1.21470194823},
        {"A:A1", 0.152496745675, 0.508322485584, 1.9345101813},
        {"A:A2", 0.508322481709, 0.508322481709, 1.0},
    };

    ln::LnSolution<double> ph_sol, cs_sol;
    for (int enc = 0; enc < 2; ++enc) {
        const bool ph = (enc == 0);
        ln::LnOptions opt;
        // The first pass leaves the method at its default ON PURPOSE: what is
        // under test there is where the 'srvn' alias lands on a setup task.
        if (!ph) opt.method = "srvn.cs";
        CAPTURE(ph ? "default (srvn.ph)" : "srvn.cs");
        ln::SolverLN<double> s(l, opt);
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        CHECK(sol.converged);

        const Row* rows = ph ? ph_rows : cs_rows;
        for (std::size_t k = 0; k < 8; ++k) {
            const Row& row = rows[k];
            const std::size_t i = idx_of(l, row.hn);
            REQUIRE(i > 0);
            CAPTURE(std::string(row.hn));
            if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
            if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
            if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
        }
        (ph ? ph_sol : cs_sol) = sol;
    }

    // The charge itself is what must not drift, so each entry response time is
    // read to six digits against its own encoding rather than to the 1e-3 the
    // table above is read to. They differ by 0.8% -- close enough that a
    // tolerance loose enough to cover both would pin neither.
    const std::size_t e2 = idx_of(l, "E:E2");
    CHECK(ph_sol.RN[e2] == doctest::Approx(1.22267268771).epsilon(1e-6));
    CHECK(cs_sol.RN[e2] == doctest::Approx(1.21470194823).epsilon(1e-6));
    CHECK(ph_sol.RN[e2] > cs_sol.RN[e2]);
}

TEST_CASE("setup task: the setup is what slows the entry down") {
    // The same model without the setup answers in exactly its service time.
    // If the setup were silently dropped the two would agree, which is the
    // failure mode this pins.
    const lqn::LqnStruct<double> plain = build_setup_model(false);
    ln::LnOptions opt;
    ln::SolverLN<double> sp(plain, opt);
    const ln::LnSolution<double> a = sp.get_ensemble_avg();
    const double respt_plain = a.RN[idx_of(plain, "E:E2")];
    CHECK(respt_plain == doctest::Approx(1.0).epsilon(1e-3));

    const lqn::LqnStruct<double> l = build_setup_model(true);
    ln::SolverLN<double> ss(l, opt);
    const ln::LnSolution<double> b = ss.get_ensemble_avg();
    const double respt_setup = b.RN[idx_of(l, "E:E2")];
    CHECK(respt_setup > respt_plain + 0.1);
    // and the throughput must fall, the server being busy restarting
    CHECK(b.TN[idx_of(l, "E:E2")] < a.TN[idx_of(plain, "E:E2")]);
}
