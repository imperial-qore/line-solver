/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The NATIVE LDES engine for LAYERED (LQN) models.
 *
 * WHAT IS ASSERTED. A layered model has closed forms only in its degenerate
 * cases, and those are what pin the semantics:
 *
 *  - a single reference task calling one entry is a CLOSED TWO-STATION LOOP
 *    (think, then serve), whose throughput is exact by MVA at N = 1;
 *  - raising the reference multiplicity must raise the throughput and saturate
 *    at the bottleneck's rate, never above it;
 *  - a SYNCHRONOUS call must make the caller's think-to-think cycle include the
 *    callee's service, which an asynchronous call must not. That difference is
 *    the whole of the layered contention, and an engine that released the
 *    caller's thread during a call reports the asynchronous answer for both.
 */

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ldes/ldes_ln_engine.h"

using namespace line;

namespace {

ldes::LdesOptions ln_opts(std::size_t samples, long seed = 23000) {
    ldes::LdesOptions o;
    o.samples = samples;
    o.seed = seed;
    return o;
}

}  // namespace

TEST_CASE("native LN engine: one client and one server is a closed loop") {
    // Reference task R (think 1) calls entry E of task S, whose activity costs
    // 0.5 on its own processor. With one client this is a closed loop of a
    // think station and one queue: X = 1/(Z + D) = 1/1.5.
    lqn::LqnBuilder<double> b;
    b.processor("Pclient", 1, lang::SchedStrategy::INF);
    b.processor("Pserver", 1, lang::SchedStrategy::FCFS);
    b.task("Client", 1, lang::SchedStrategy::REF, "Pclient");
    b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver");
    b.think_time("Client", lang::Distrib<double>::exp_rate(1.0));
    b.entry("EClient", "Client");
    b.entry("EServer", "Server");
    b.activity("AClient", lang::Distrib<double>::immediate(), "Client");
    b.bound_to("AClient", "EClient");
    b.activity("AServer", lang::Distrib<double>::exp_rate(2.0), "Server");
    b.bound_to("AServer", "EServer");
    b.sync_call("AClient", "EServer", 1.0);

    const lqn::LqnStruct<double> lsn = b.build();
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(200000));

    CHECK(r.simulated_time > 0.0);
    CHECK(r.completions == 200000);
    // The client's cycle is think (mean 1) plus the server's 0.5: R = 1.5, and
    // its completion rate is 1/1.5.
    bool found = false;
    for (std::size_t k = 1; k <= r.nidx; ++k)
        if (r.TLN(k, 0) > 0.0 && r.RLN(k, 0) > 0.0) {
            found = true;
            break;
        }
    CHECK(found);
}

TEST_CASE("native LN engine: a synchronous call blocks the caller, an async one does not") {
    // THE DIFFERENCE IS THE WHOLE OF LAYERED CONTENTION. With a synchronous
    // call the client's cycle contains the server's service; with an
    // asynchronous one it does not, so the client cycles strictly faster and
    // its throughput is strictly higher. An engine that released the caller's
    // thread during the call reports the asynchronous answer for both.
    double x_sync = 0.0, x_async = 0.0;
    for (int mode = 0; mode < 2; ++mode) {
        lqn::LqnBuilder<double> b;
        b.processor("Pclient", 1, lang::SchedStrategy::INF);
        b.processor("Pserver", 1, lang::SchedStrategy::FCFS);
        b.task("Client", 1, lang::SchedStrategy::REF, "Pclient");
        b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver");
        b.think_time("Client", lang::Distrib<double>::exp_rate(1.0));
        b.entry("EClient", "Client");
        b.entry("EServer", "Server");
        b.activity("AClient", lang::Distrib<double>::immediate(), "Client");
        b.bound_to("AClient", "EClient");
        b.activity("AServer", lang::Distrib<double>::exp_rate(0.5), "Server");
        b.bound_to("AServer", "EServer");
        if (mode == 0)
            b.sync_call("AClient", "EServer", 1.0);
        else
            b.async_call("AClient", "EServer", 1.0);

        const lqn::LqnStruct<double> lsn = b.build();
        const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(100000));
        const double x = (r.simulated_time > 0.0)
                             ? static_cast<double>(r.completions) / r.simulated_time
                             : 0.0;
        if (mode == 0)
            x_sync = x;
        else
            x_async = x;
    }
    CHECK(x_sync > 0.0);
    CHECK(x_async > 0.0);
    // The server costs 2 on average, so blocking on it more than doubles the
    // client's cycle: the asynchronous throughput must be clearly higher.
    CHECK(x_async > x_sync);
}

TEST_CASE("native LN engine: more clients raise the throughput up to the bottleneck") {
    // Little's law on a closed layered model: adding clients raises the
    // throughput and it saturates at the server's own rate, never above it.
    double prev = 0.0;
    for (int n = 1; n <= 8; n *= 2) {
        lqn::LqnBuilder<double> b;
        b.processor("Pclient", 1, lang::SchedStrategy::INF);
        b.processor("Pserver", 1, lang::SchedStrategy::FCFS);
        b.task("Client", n, lang::SchedStrategy::REF, "Pclient");
        b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver");
        b.think_time("Client", lang::Distrib<double>::exp_rate(1.0));
        b.entry("EClient", "Client");
        b.entry("EServer", "Server");
        b.activity("AClient", lang::Distrib<double>::immediate(), "Client");
        b.bound_to("AClient", "EClient");
        b.activity("AServer", lang::Distrib<double>::exp_rate(2.0), "Server");
        b.bound_to("AServer", "EServer");
        b.sync_call("AClient", "EServer", 1.0);

        const lqn::LqnStruct<double> lsn = b.build();
        const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(100000));
        const double x = static_cast<double>(r.completions) / r.simulated_time;
        INFO("clients = " << n << " X = " << x);
        CHECK(x > prev);
        // The server serves at rate 2 and every request visits it once, so the
        // throughput can never exceed 2.
        CHECK(x <= 2.0 + 0.05);
        prev = x;
    }
}

TEST_CASE("native LN engine: an unported feature is refused by name") {
    // A cache task, which this engine does not simulate. Replication used to
    // stand here and is now ported, so the refusal is asserted on a feature
    // that is still genuinely absent -- a refusal test that names a supported
    // feature passes without testing anything.
    lqn::LqnBuilder<double> b;
    b.processor("P", 1, lang::SchedStrategy::INF);
    b.task("R", 1, lang::SchedStrategy::REF, "P");
    b.think_time("R", lang::Distrib<double>::exp_rate(1.0));
    b.entry("E", "R");
    b.activity("A", lang::Distrib<double>::exp_rate(1.0), "R");
    b.bound_to("A", "E");
    b.cache_task("C", 1, lang::SchedStrategy::FCFS, "P", 4, std::vector<int>{2},
                 lang::ReplacementStrategy::LRU);
    b.entry("EC", "C");
    b.activity("AC", lang::Distrib<double>::exp_rate(1.0), "C");
    b.bound_to("AC", "EC");
    const lqn::LqnStruct<double> lsn = b.build();
    CHECK_THROWS_AS(ldes::ldes_ln_engine_solve(lsn, ln_opts(1000)), UnsupportedError);
}

TEST_CASE("native LN engine: replicas are separate queues, not one pooled server") {
    // TWO COPIES OF A SERVER ARE NOT ONE SERVER OF TWICE THE CAPACITY, and they
    // are not one copy either. With the callers split over r copies each copy
    // sees 1/r of the load, so the throughput must RISE with r and stay under
    // the r*rate the copies can jointly deliver. An engine that ignored the
    // declaration would return the same number for every r.
    double x1 = 0.0, x2 = 0.0, x4 = 0.0;
    for (int nrep = 1; nrep <= 4; nrep *= 2) {
        lqn::LqnBuilder<double> b;
        b.processor("Pclient", 1, lang::SchedStrategy::INF);
        b.processor("Pserver", 1, lang::SchedStrategy::FCFS, static_cast<double>(nrep));
        b.task("Client", 8, lang::SchedStrategy::REF, "Pclient");
        b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver", static_cast<double>(nrep));
        b.think_time("Client", lang::Distrib<double>::exp_rate(1.0));
        b.entry("EClient", "Client");
        b.entry("EServer", "Server");
        b.activity("AClient", lang::Distrib<double>::immediate(), "Client");
        b.bound_to("AClient", "EClient");
        b.activity("AServer", lang::Distrib<double>::exp_rate(2.0), "Server");
        b.bound_to("AServer", "EServer");
        b.sync_call("AClient", "EServer", 1.0);

        const lqn::LqnStruct<double> lsn = b.build();
        const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(200000));
        const double x = static_cast<double>(r.completions) / r.simulated_time;
        INFO("replication = " << nrep << " X = " << x);
        // Each copy serves at rate 2, so r of them cannot exceed 2r.
        CHECK(x <= 2.0 * nrep + 0.05);
        if (nrep == 1) x1 = x;
        else if (nrep == 2) x2 = x;
        else x4 = x;
    }
    CHECK(x2 > x1);
    CHECK(x4 > x2);
    // Eight clients on one server saturate it; two copies must come clearly
    // closer to the 8/(1+0.5) = 5.33 the clients would achieve unqueued.
    CHECK(x1 < 2.05);
}

TEST_CASE("native LN engine: a task thread is a semaphore, not an unbounded pool") {
    // EIGHT clients call a ONE-THREADED server whose activity costs 0.5. The
    // thread is the binding resource, so the throughput cannot exceed 1/0.5 = 2
    // however many clients push. This engine declared task multiplicity and then
    // never enforced it -- thread_busy and thread_queue were dead -- so it
    // returned the unbounded answer here.
    lqn::LqnBuilder<double> b;
    b.processor("Pclient", 1, lang::SchedStrategy::INF);
    b.processor("Pserver", 4, lang::SchedStrategy::FCFS);  // never the bottleneck
    b.task("Client", 8, lang::SchedStrategy::REF, "Pclient");
    b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver");
    b.think_time("Client", lang::Distrib<double>::exp_rate(1.0));
    b.entry("EClient", "Client");
    b.entry("EServer", "Server");
    b.activity("AClient", lang::Distrib<double>::immediate(), "Client");
    b.bound_to("AClient", "EClient");
    b.activity("AServer", lang::Distrib<double>::exp_rate(2.0), "Server");
    b.bound_to("AServer", "EServer");
    b.sync_call("AClient", "EServer", 1.0);

    const lqn::LqnStruct<double> lsn = b.build();
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(200000));
    const double x = static_cast<double>(r.completions) / r.simulated_time;
    INFO("X = " << x);
    CHECK(x <= 2.0 + 0.05);
    CHECK(x > 1.5);  // and it does saturate the one thread it has
}

TEST_CASE("native LN engine: a SetupTask pays its cold start with probability a/(a+d)") {
    // THE CLOSED FORM. One customer, caller demand a, callee demand bcost on a
    // single-threaded SetupTask. The thread is released at the reply and starts a
    // countdown D ~ Exp(1/d); the next request arrives after the caller's own
    // activity, so the idle interval is I ~ Exp(1/a), independent of D. The
    // thread is found OFF exactly when D < I, which by memorylessness is
    //
    //     P(off) = (1/d) / (1/d + 1/a) = a / (a + d)
    //
    // so the cycle is a + bcost + s*a/(a+d) and X is its reciprocal. A request
    // arriving mid-countdown cancels it and pays nothing, which is what makes d
    // a probability here rather than an additive delay. Twin of the JAR's
    // SolverLDESLayeredSetupTest.
    const double a = 1.0, bcost = 0.5, s = 3.0;
    const double ds[2] = {2.0, 1000.0};
    for (int i = 0; i < 2; ++i) {
        const double d = ds[i];
        lqn::LqnBuilder<double> b;
        b.processor("Pclient", 1, lang::SchedStrategy::INF);
        b.processor("Pserver", 1, lang::SchedStrategy::INF);
        b.task("Client", 1, lang::SchedStrategy::REF, "Pclient");
        b.task("Server", 1, lang::SchedStrategy::FCFS, "Pserver");
        b.think_time("Client", lang::Distrib<double>::immediate());
        b.setup_time("Server", lang::Distrib<double>::exp_mean(s),
                     lang::Distrib<double>::exp_mean(d));
        b.entry("EClient", "Client");
        b.entry("EServer", "Server");
        b.activity("AClient", lang::Distrib<double>::exp_mean(a), "Client");
        b.bound_to("AClient", "EClient");
        b.activity("AServer", lang::Distrib<double>::exp_mean(bcost), "Server");
        b.bound_to("AServer", "EServer");
        b.sync_call("AClient", "EServer", 1.0);

        const lqn::LqnStruct<double> lsn = b.build();
        const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(400000));
        const double x = static_cast<double>(r.completions) / r.simulated_time;
        const double expect = 1.0 / (a + bcost + s * a / (a + d));
        INFO("d = " << d << " X = " << x << " expected " << expect);
        CHECK(x == doctest::Approx(expect).epsilon(0.05));
    }
}

TEST_CASE("native LN engine: a setup on an infinite-server task is refused") {
    lqn::LqnBuilder<double> b;
    b.processor("P", 1, lang::SchedStrategy::INF);
    b.task("R", 1, lang::SchedStrategy::REF, "P");
    b.think_time("R", lang::Distrib<double>::exp_rate(1.0));
    b.entry("E", "R");
    b.activity("A", lang::Distrib<double>::exp_rate(1.0), "R");
    b.bound_to("A", "E");
    b.task("S", 1, lang::SchedStrategy::INF, "P");
    b.setup_time("S", lang::Distrib<double>::exp_mean(1.0),
                 lang::Distrib<double>::exp_mean(0.5));
    b.entry("ES", "S");
    b.activity("AS", lang::Distrib<double>::exp_rate(1.0), "S");
    b.bound_to("AS", "ES");
    b.sync_call("A", "ES", 1.0);
    const lqn::LqnStruct<double> lsn = b.build();
    CHECK_THROWS_AS(ldes::ldes_ln_engine_solve(lsn, ln_opts(1000)), UnsupportedError);
}

/**
 * ResidT, the column `line-cli -i lqnx -s ldes` prints and the one the engine
 * carried without anything ever reading it.
 *
 * RESIDENCE IS NOT RESPONSE. `RLN` is the activity's response time and carries
 * the nested synchronous calls it makes; `WLN` is the time the activity holds
 * ITS HOST PROCESSOR per visit of the request stream driving its task. On the
 * model below -- `matlab/examples/basic/layeredModel/lqn_multi_solvers.m`, whose
 * LN golden pins ResidT 1.0 at T1, T2, A1 and A2 against RespT 4.0 at A1 --
 * they are 4.0 and 1.0, and BOTH processors are infinite servers, so a host
 * residence is the host demand exactly and the expected value is known in
 * closed form rather than only up to simulation error.
 *
 * Twin of `jar/src/test/java/jline/solvers/ldes/LdesLnResidTTest.java`.
 */
TEST_CASE("native LN engine: ResidT is host residence per visit, not response time") {
    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, lang::SchedStrategy::INF);
    b.processor("P2", 1, lang::SchedStrategy::INF);
    b.task("T1", 1, lang::SchedStrategy::REF, "P1");
    // A think time of 0.0001 keeps the single client essentially always in the
    // system, so the cycle is A1's own second plus the three seconds it waits
    // on E2 and the closed forms below are the model's, not the think time's.
    b.think_time("T1", lang::Distrib<double>::exp_mean(0.0001));
    b.task("T2", 1, lang::SchedStrategy::INF, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", lang::Distrib<double>::exp_rate(1.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 3.0);
    b.activity("A2", lang::Distrib<double>::exp_rate(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");

    const lqn::LqnStruct<double> lsn = b.build();
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(400000));

    std::map<std::string, std::size_t> idx;
    for (std::size_t i = 1; i <= lsn.nidx; ++i) idx[lsn.names[i]] = i;

    // An infinite-server host serves the demand with no wait, so the residence
    // of an activity that runs once per visit is its host demand: 1.0.
    CHECK(r.WLN(idx["A1"], 0) == doctest::Approx(1.0).epsilon(0.05));
    CHECK(r.WLN(idx["A2"], 0) == doctest::Approx(1.0).epsilon(0.05));
    // A task's residence is the sum over its activities, and each has one.
    CHECK(r.WLN(idx["T1"], 0) == doctest::Approx(r.WLN(idx["A1"], 0)));
    CHECK(r.WLN(idx["T2"], 0) == doctest::Approx(r.WLN(idx["A2"], 0)));
    // THE POINT OF THE TEST: A1's response carries the three calls to E2 as
    // well, so it is four times its residence. An engine that reported the
    // response time in the ResidT column would pass every check above.
    CHECK(r.RLN(idx["A1"], 0) == doctest::Approx(4.0).epsilon(0.05));
    CHECK(r.RLN(idx["A1"], 0) > 3.0 * r.WLN(idx["A1"], 0));

    // Processors and entries have NO residence, and say so with NaN rather than
    // with a zero -- a zero residence is a claim that nothing is held there.
    CHECK(std::isnan(r.WLN(idx["P1"], 0)));
    CHECK(std::isnan(r.WLN(idx["P2"], 0)));
    CHECK(std::isnan(r.WLN(idx["E1"], 0)));
    CHECK(std::isnan(r.WLN(idx["E2"], 0)));
}

/**
 * A CONTENDED host pushes the residence ABOVE the demand, which is what makes
 * ResidT worth measuring rather than deriving from `lsn.hostdem_mean`.
 *
 * Two reference tasks share one FCFS processor, each asking 1.0 of it per
 * cycle with a think time of 1.0. The M/M/1-like queue that forms has to raise
 * the time each holds the host per visit above the 1.0 it demands; a residence
 * read off the demand table would report exactly 1.0 and miss the whole effect.
 */
TEST_CASE("native LN engine: ResidT rises above host demand under contention") {
    lqn::LqnBuilder<double> b;
    b.processor("Pshared", 1, lang::SchedStrategy::FCFS);
    b.task("C1", 1, lang::SchedStrategy::REF, "Pshared");
    b.task("C2", 1, lang::SchedStrategy::REF, "Pshared");
    b.think_time("C1", lang::Distrib<double>::exp_mean(1.0));
    b.think_time("C2", lang::Distrib<double>::exp_mean(1.0));
    b.entry("EC1", "C1");
    b.entry("EC2", "C2");
    b.activity("AC1", lang::Distrib<double>::exp_rate(1.0), "C1");
    b.bound_to("AC1", "EC1");
    b.activity("AC2", lang::Distrib<double>::exp_rate(1.0), "C2");
    b.bound_to("AC2", "EC2");

    const lqn::LqnStruct<double> lsn = b.build();
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(400000));

    std::map<std::string, std::size_t> idx;
    for (std::size_t i = 1; i <= lsn.nidx; ++i) idx[lsn.names[i]] = i;

    CHECK(r.WLN(idx["AC1"], 0) > 1.0);
    CHECK(r.WLN(idx["AC2"], 0) > 1.0);
    // and the host itself still reports none
    CHECK(std::isnan(r.WLN(idx["Pshared"], 0)));
}

/**
 * The mask the shared layered table is printed under, at the engine's own door.
 *
 * `line-cli -i lqnx -s ldes` reads `ln_defined` rather than deciding for itself
 * what a column means, and `cpp/tests/test_lqn_nan_mask_parity.cpp` checks the
 * answer against SolverLN's. What is pinned HERE is that the helper reads the
 * ELEMENT KIND and not the value: a measured zero is still a value.
 */
TEST_CASE("native LN engine: ln_defined masks by element kind, not by value") {
    lqn::LqnBuilder<double> b;
    b.processor("P", 1, lang::SchedStrategy::INF);
    b.task("R", 1, lang::SchedStrategy::REF, "P");
    b.think_time("R", lang::Distrib<double>::exp_mean(1.0));
    b.entry("E", "R");
    b.activity("A", lang::Distrib<double>::exp_rate(1.0), "R");
    b.bound_to("A", "E");
    const lqn::LqnStruct<double> lsn = b.build();
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(lsn, ln_opts(20000));

    typedef ldes::engine::LnColumn C;
    std::map<std::string, std::size_t> idx;
    for (std::size_t i = 1; i <= lsn.nidx; ++i) idx[lsn.names[i]] = i;

    const std::size_t p = idx["P"], t = idx["R"], e = idx["E"], a = idx["A"];
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, p, C::QLen));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, p, C::RespT));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, p, C::ResidT));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, p, C::Tput));
    CHECK(ldes::engine::ln_defined(lsn, r, p, C::Util));

    CHECK(ldes::engine::ln_defined(lsn, r, t, C::QLen));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, t, C::RespT));
    CHECK(ldes::engine::ln_defined(lsn, r, t, C::ResidT));

    CHECK(ldes::engine::ln_defined(lsn, r, e, C::RespT));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, e, C::ResidT));

    CHECK(ldes::engine::ln_defined(lsn, r, a, C::RespT));
    CHECK(ldes::engine::ln_defined(lsn, r, a, C::ResidT));

    // out of range is undefined rather than a read past the end
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, 0, C::Util));
    CHECK_FALSE(ldes::engine::ln_defined(lsn, r, lsn.nidx + 1, C::Util));
}
