/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `solver_mam_bgchain`: the mixed-network analyzer that solves the closed
 * population vector as a background modulating chain and each open station as a
 * level-dependent QBD driven by it.
 *
 * THE ORACLES, in decreasing order of independence.
 *
 *  1. EXACT MVA. Every PS model here is product-form (class-independent service
 *     rates at the shared station), so MVA is exact and the method has no way to
 *     know it: the background chain and the QBD are built from scratch. The
 *     reference measures agreement to five significant digits on exactly these
 *     shapes, so they are asserted at 1e-4 relative.
 *  2. IDENTITIES the analyzer does not impose: flow balance at the open source,
 *     population conservation of each closed chain, and the utilization law
 *     U = T S / c on the open classes.
 *  3. MATLAB, for the FCFS cases with no closed form. Every such
 *     value was produced by
 *     `SolverMAM(model,'method','bgchain').getAvgTable()` against
 *     matlab/src/solvers/MAM/solver_mam_bgchain.m, and the JAR twin reproduces
 *     them to 3e-16.
 *
 * Each analyzer is also called directly on its own, so a lumping or a
 * tagged-class defect is attributed to the piece that owns it rather than to the
 * fixed point that consumes it.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/mam_dispatch.h"
#include "line/solvers/mam/solver_mam_bgchain.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/mva_dispatch.h"

using namespace line;
using Dd = lang::Distrib<double>;
using lang::SchedStrategy;

namespace {

/** Source -> Q -> Sink for the open class, Delay -> Q -> Delay for the closed one. */
qn::Network<double> bg_sdq(const std::string& name, SchedStrategy sched, double lambda,
                           const Dd& svcOpen, const Dd& svcClosed, double N, double think,
                           double servers = 1.0) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", sched);
    const std::size_t k = m.add_sink("Sink");
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    const std::size_t o = m.add_open_class("Open");
    const std::size_t c = m.add_closed_class("Closed", N, d);
    m.set_arrival(s, o, Dd::exp_rate(lambda));
    m.set_service(q, o, svcOpen);
    m.set_service(d, c, Dd::exp_rate(think));
    m.set_service(q, c, svcClosed);
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q, 1.0);
    P.set(o, o, q, k, 1.0);
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

mva::MvaSolution<double> run_bg(qn::Network<double>& m) {
    mam::MamOptions opt;
    opt.method = "bgchain";
    return mam::solver_mam_bgchain(m.get_struct(), opt);
}

const double kRel = 1e-4;

}  // namespace

// ---------------------------------------------------------------------------
// Against exact MVA, on the product-form shapes
// ---------------------------------------------------------------------------

TEST_CASE("bgchain: one open + one closed class at a PS station matches exact MVA") {
    // Source(0.4) -> PS(mu=2) -> Sink, and a closed class of 2 cycling
    // Delay(1) -> PS(mu=2). Class-independent rates under PS, so the network is
    // product-form and MVA gives Delay QLen 1.0722362, Queue Open 0.4819594,
    // Queue Closed 0.9277638.
    qn::Network<double> m =
        bg_sdq("bgA", SchedStrategy::PS, 0.4, Dd::exp_rate(2.0), Dd::exp_rate(2.0), 2.0, 1.0);
    const mva::MvaSolution<double> r = run_bg(m);

    CHECK(r.Q(1, 1) == doctest::Approx(1.0722362).epsilon(kRel));  // Delay, closed
    CHECK(r.Q(2, 0) == doctest::Approx(0.4819594).epsilon(kRel));  // Queue, open
    CHECK(r.Q(2, 1) == doctest::Approx(0.9277638).epsilon(kRel));  // Queue, closed

    // Population conservation of the closed chain, which the analyzer never
    // imposes: the two stations are the only places its 2 jobs can be.
    CHECK(r.Q(1, 1) + r.Q(2, 1) == doctest::Approx(2.0).epsilon(1e-6));
    // Flow balance: the open class leaves at the rate it arrives.
    CHECK(r.Tp(2, 0) == doctest::Approx(0.4).epsilon(1e-9));
    // Utilization law on the open class: U = T S / c.
    CHECK(r.U(2, 0) == doctest::Approx(0.4 * 0.5).epsilon(1e-9));
}

TEST_CASE("bgchain: a multiserver PS station matches exact MVA") {
    // Two servers, class-independent rates: still product-form.
    qn::Network<double> m = bg_sdq("bgMS", SchedStrategy::PS, 0.8, Dd::exp_rate(1.5),
                                   Dd::exp_rate(1.5), 3.0, 1.0, 2.0);
    const mva::MvaSolution<double> r = run_bg(m);
    CHECK(r.Q(1, 1) + r.Q(2, 1) == doctest::Approx(3.0).epsilon(1e-6));
    CHECK(r.Tp(2, 0) == doctest::Approx(0.8).epsilon(1e-9));
    CHECK(r.U(2, 0) == doctest::Approx(0.8 / 1.5 / 2.0).epsilon(1e-9));
    // MVA: Delay 1.5438843, Queue open 0.9282210, Queue closed 1.4561157.
    CHECK(r.Q(1, 1) == doctest::Approx(1.5438843).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.9282210).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(1.4561157).epsilon(kRel));
}

TEST_CASE("bgchain: FCFS with class-independent rates matches exact MVA") {
    // FCFS with one shared rate is distributionally the same queue as PS here,
    // so the product-form answer still applies and the random-order surrogate
    // costs nothing.
    qn::Network<double> m =
        bg_sdq("bgFC", SchedStrategy::FCFS, 0.5, Dd::exp_rate(2.0), Dd::exp_rate(2.0), 3.0, 1.0);
    const mva::MvaSolution<double> r = run_bg(m);
    CHECK(r.Q(1, 1) == doctest::Approx(1.2985075).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.9004975).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(1.7014925).epsilon(kRel));
    CHECK(r.Q(1, 1) + r.Q(2, 1) == doctest::Approx(3.0).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// Against MATLAB, on the shapes with no closed form
// ---------------------------------------------------------------------------

TEST_CASE("bgchain: FCFS with class-dependent rates and two servers matches MATLAB") {
    // SolverMAM(model,'method','bgchain') in MATLAB, reproduced by the JAR to
    // 3e-16: Delay 1.923295, Queue open 0.210412, Queue closed 1.076705.
    qn::Network<double> m = bg_sdq("bgFCd", SchedStrategy::FCFS, 0.5, Dd::exp_rate(3.0),
                                   Dd::exp_rate(2.0), 3.0, 1.0, 2.0);
    const mva::MvaSolution<double> r = run_bg(m);
    CHECK(r.Q(1, 1) == doctest::Approx(1.923295).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.210412).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(1.076705).epsilon(kRel));
    CHECK(r.Q(1, 1) + r.Q(2, 1) == doctest::Approx(3.0).epsilon(1e-6));
}

TEST_CASE("bgchain: a PS station is insensitive to the service law beyond its mean") {
    // Erlang-2 on both classes at a PS station. PS is INSENSITIVE, so the exact
    // answer is the one the same means give with exponential service, and
    // SolverCTMC confirms it: Delay 1.09589, Queue open 0.38082, Queue closed
    // 0.90411 for both service laws. The method honours that by exponentializing
    // the open service at a PS station -- carrying the phase-type there instead
    // made the queue length inherit the SCV-sensitivity of an M/PH/1 FCFS queue,
    // reading 4% low here and 21% high on a HyperExp of SCV 4.
    qn::Network<double> m =
        bg_sdq("bgPH", SchedStrategy::PS, 0.5, Dd::erlang(6.0, 2),
               Dd::erlang(4.0, 2), 2.0, 1.0);
    const mva::MvaSolution<double> r = run_bg(m);
    CHECK(r.Q(1, 1) == doctest::Approx(1.095947).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.380822).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(0.904053).epsilon(kRel));
    CHECK(r.Q(1, 1) + r.Q(2, 1) == doctest::Approx(2.0).epsilon(1e-6));

    // The same model with EXPONENTIAL service of the same means must give the
    // same answer. That is the insensitivity claim itself, rather than a
    // restatement of the numbers above.
    qn::Network<double> mexp =
        bg_sdq("bgPHexp", SchedStrategy::PS, 0.5, Dd::exp_rate(3.0),
               Dd::exp_rate(2.0), 2.0, 1.0);
    const mva::MvaSolution<double> rexp = run_bg(mexp);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(r.Q(i, c) == doctest::Approx(rexp.Q(i, c)).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// The tagged-class iteration
// ---------------------------------------------------------------------------

TEST_CASE("bgchain: two closed chains are carried by the tagged-class iteration") {
    // Source(0.3) -> PS -> PS -> Sink for the open class; two closed chains of
    // 2 and 1 jobs cycling Delay -> PS -> PS. Product-form, so exact MVA gives
    // Delay A 0.6941265, Delay B 0.3470633, Queue1 open 0.3493137,
    // Queue1 A 0.6529316, Queue1 B 0.3264658, and Queue2 the same as Queue1.
    qn::Network<double> m("bgTag");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("Open");
    const std::size_t ca = m.add_closed_class("ClosedA", 2.0, d);
    const std::size_t cb = m.add_closed_class("ClosedB", 1.0, d);
    m.set_arrival(s, o, Dd::exp_rate(0.3));
    m.set_service(q1, o, Dd::exp_rate(2.0));
    m.set_service(q2, o, Dd::exp_rate(2.0));
    m.set_service(d, ca, Dd::exp_rate(1.0));
    m.set_service(q1, ca, Dd::exp_rate(2.0));
    m.set_service(q2, ca, Dd::exp_rate(2.0));
    m.set_service(d, cb, Dd::exp_rate(1.0));
    m.set_service(q1, cb, Dd::exp_rate(2.0));
    m.set_service(q2, cb, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q1, 1.0);
    P.set(o, o, q1, q2, 1.0);
    P.set(o, o, q2, k, 1.0);
    P.set(ca, ca, d, q1, 1.0);
    P.set(ca, ca, q1, q2, 1.0);
    P.set(ca, ca, q2, d, 1.0);
    P.set(cb, cb, d, q1, 1.0);
    P.set(cb, cb, q1, q2, 1.0);
    P.set(cb, cb, q2, d, 1.0);
    m.link(P);

    const mva::MvaSolution<double> r = run_bg(m);
    // Station order is Source(0), Delay(1), Queue1(2), Queue2(3); class order
    // Open(0), ClosedA(1), ClosedB(2).
    CHECK(r.Q(1, 1) == doctest::Approx(0.6941265).epsilon(kRel));
    CHECK(r.Q(1, 2) == doctest::Approx(0.3470633).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.3493137).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(0.6529316).epsilon(kRel));
    CHECK(r.Q(2, 2) == doctest::Approx(0.3264658).epsilon(kRel));
    CHECK(r.Q(3, 0) == doctest::Approx(0.3493137).epsilon(kRel));
    // Each closed chain conserves its own population across the three stations.
    CHECK(r.Q(1, 1) + r.Q(2, 1) + r.Q(3, 1) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(r.Q(1, 2) + r.Q(2, 2) + r.Q(3, 2) == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("bgchain: three closed chains aggregate two at a time") {
    // R = 3 is the first population where the tagged iteration has to aggregate
    // MORE THAN ONE chain, so the flow-equivalent weights actually mix. Rates
    // are class-DEPENDENT at the queue, so the model is not product-form and
    // MVA is no oracle: the reference is SolverCTMC at cutoff 12, which gives
    // Delay A 0.49582, B 0.43573, C 0.58149, Q1 open 0.27633, A 0.50418,
    // B 0.56427, C 0.41851. MATLAB bgchain lands within 0.11% of it and the JAR
    // reproduces MATLAB to 6e-16, so this asserts the MATLAB values at 1e-4.
    qn::Network<double> m("bgR3");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("Open");
    const std::size_t ca = m.add_closed_class("ClosedA", 1.0, d);
    const std::size_t cb = m.add_closed_class("ClosedB", 1.0, d);
    const std::size_t cc = m.add_closed_class("ClosedC", 1.0, d);
    m.set_arrival(s, o, Dd::exp_rate(0.4));
    m.set_service(q, o, Dd::exp_rate(4.0));
    m.set_service(d, ca, Dd::exp_rate(1.0));
    m.set_service(q, ca, Dd::exp_rate(2.0));
    m.set_service(d, cb, Dd::exp_rate(2.0));
    m.set_service(q, cb, Dd::exp_rate(3.0));
    m.set_service(d, cc, Dd::exp_rate(0.5));
    m.set_service(q, cc, Dd::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q, 1.0);
    P.set(o, o, q, k, 1.0);
    P.set(ca, ca, d, q, 1.0);
    P.set(ca, ca, q, d, 1.0);
    P.set(cb, cb, d, q, 1.0);
    P.set(cb, cb, q, d, 1.0);
    P.set(cc, cc, d, q, 1.0);
    P.set(cc, cc, q, d, 1.0);
    m.link(P);

    const mva::MvaSolution<double> r = run_bg(m);
    // Station order Source(0), Delay(1), Q1(2); class order Open(0), A(1), B(2), C(3).
    CHECK(r.Q(1, 1) == doctest::Approx(0.4956555).epsilon(kRel));
    CHECK(r.Q(1, 2) == doctest::Approx(0.4357193).epsilon(kRel));
    CHECK(r.Q(1, 3) == doctest::Approx(0.5814728).epsilon(kRel));
    CHECK(r.Q(2, 0) == doctest::Approx(0.2766369).epsilon(kRel));
    CHECK(r.Q(2, 1) == doctest::Approx(0.5043445).epsilon(kRel));
    CHECK(r.Q(2, 2) == doctest::Approx(0.5642807).epsilon(kRel));
    CHECK(r.Q(2, 3) == doctest::Approx(0.4185272).epsilon(kRel));
    // Each of the three chains conserves its own single job.
    for (std::size_t c = 1; c <= 3; ++c)
        CHECK(r.Q(1, c) + r.Q(2, c) == doctest::Approx(1.0).epsilon(1e-6));
    // Flow balance on the open class, and its utilization law.
    CHECK(r.Tp(2, 0) == doctest::Approx(0.4).epsilon(1e-9));
    CHECK(r.U(2, 0) == doctest::Approx(0.4 / 4.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// config.bgaggr, and the per-class station support it stresses
// ---------------------------------------------------------------------------

namespace {

/**
 * Source -> Q1 -> Q2 -> Sink for the open class; three closed chains with
 * DIFFERENT routes -- one confined to Q1, one to Q2, one crossing both. This is
 * the shape where aggregating a group really distorts, because the aggregate
 * carries a flow-weighted MEAN of its members' routing matrices, and it is also
 * the shape that requires each class to be enumerated over its OWN stations.
 */
qn::Network<double> bg_routes() {
    qn::Network<double> m("bgRoutes");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("Open");
    const std::size_t e1 = m.add_closed_class("OnQ1", 1.0, d);
    const std::size_t e2 = m.add_closed_class("OnQ2", 1.0, d);
    const std::size_t e3 = m.add_closed_class("Both", 1.0, d);
    m.set_arrival(s, o, Dd::exp_rate(0.4));
    m.set_service(q1, o, Dd::exp_rate(4.0));
    m.set_service(q2, o, Dd::exp_rate(4.0));
    m.set_service(d, e1, Dd::exp_rate(1.0));
    m.set_service(q1, e1, Dd::exp_rate(1.5));
    m.set_service(d, e2, Dd::exp_rate(1.0));
    m.set_service(q2, e2, Dd::exp_rate(1.5));
    m.set_service(d, e3, Dd::exp_rate(1.0));
    m.set_service(q1, e3, Dd::exp_rate(3.0));
    m.set_service(q2, e3, Dd::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q1, 1.0);
    P.set(o, o, q1, q2, 1.0);
    P.set(o, o, q2, k, 1.0);
    P.set(e1, e1, d, q1, 1.0);
    P.set(e1, e1, q1, d, 1.0);
    P.set(e2, e2, d, q2, 1.0);
    P.set(e2, e2, q2, d, 1.0);
    P.set(e3, e3, d, q1, 1.0);
    P.set(e3, e3, q1, q2, 1.0);
    P.set(e3, e3, q2, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("bgchain: a chain confined to some stations conserves its population") {
    // Each closed chain holds ONE job. Enumerating a chain over the union of
    // every chain's stations would place probability on stations it never
    // reaches, and those states ABSORB (a zero routing row normalizes to a
    // self-loop), so the mass would leak: measured, HALF of it. Population
    // conservation is therefore the test that the per-class station support is
    // being honoured, and it holds whatever bgaggr does.
    qn::Network<double> m = bg_routes();
    for (std::size_t g = 1; g <= 3; ++g) {
        mam::MamOptions opt;
        opt.method = "bgchain";
        opt.bgaggr = g;
        const mva::MvaSolution<double> r = mam::solver_mam_bgchain(m.get_struct(), opt);
        // Stations Source(0), Delay(1), Q1(2), Q2(3); classes Open(0), OnQ1(1),
        // OnQ2(2), Both(3).
        for (std::size_t c = 1; c <= 3; ++c) {
            const double held = r.Q(1, c) + r.Q(2, c) + r.Q(3, c);
            CHECK(held == doctest::Approx(1.0).epsilon(1e-6));
        }
        // A chain confined to Q1 holds nothing at Q2, and conversely.
        CHECK(r.Q(3, 1) == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
        CHECK(r.Q(2, 2) == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
    }
}

TEST_CASE("bgchain: raising bgaggr separates the chains an aggregate misrepresents") {
    // SolverCTMC at cutoff 8 on the same model: Delay OnQ1 0.53043, Delay Both
    // 0.48639, Q1 OnQ1 0.46957, Q1 Both 0.25680. With ONE aggregate the two
    // untagged chains share a mean routing matrix they have no business sharing
    // and the worst queue length is ~1.7% out; with two, every chain is carried
    // exactly and it is ~0.01%. That gap IS the feature.
    qn::Network<double> m = bg_routes();
    mam::MamOptions opt1;
    opt1.method = "bgchain";
    opt1.bgaggr = 1;
    const mva::MvaSolution<double> r1 = mam::solver_mam_bgchain(m.get_struct(), opt1);
    mam::MamOptions opt2;
    opt2.method = "bgchain";
    opt2.bgaggr = 2;
    const mva::MvaSolution<double> r2 = mam::solver_mam_bgchain(m.get_struct(), opt2);

    const double refDelayBoth = 0.48639;
    const double err1 = std::fabs(r1.Q(1, 3) - refDelayBoth) / refDelayBoth;
    const double err2 = std::fabs(r2.Q(1, 3) - refDelayBoth) / refDelayBoth;
    CHECK(err1 > 5e-3);   // one aggregate cannot represent both routes
    CHECK(err2 < 5e-4);   // two carry them exactly
    CHECK(err2 < err1);

    // bgaggr above R-1 clamps rather than erroring, so "no aggregation" needs no
    // magic value: 3 and 2 must answer identically here (R = 3).
    mam::MamOptions opt3;
    opt3.method = "bgchain";
    opt3.bgaggr = 3;
    const mva::MvaSolution<double> r3 = mam::solver_mam_bgchain(m.get_struct(), opt3);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t c = 0; c < 4; ++c)
            CHECK(r3.Q(i, c) == doctest::Approx(r2.Q(i, c)).epsilon(1e-12));
}

TEST_CASE("mam_bgchain_groups: demand similarity drives the grouping") {
    // Four chains over two stations: two near-twins at a light demand, two at a
    // heavy one. Asking for two groups must split them that way, and the labels
    // are canonical (group 0 owns the smallest member).
    std::vector<std::vector<double>> D(2, std::vector<double>(4, 0.0));
    D[0][0] = 1.00; D[1][0] = 0.20;
    D[0][1] = 1.02; D[1][1] = 0.21;
    D[0][2] = 4.00; D[1][2] = 3.90;
    D[0][3] = 4.10; D[1][3] = 3.95;
    const std::vector<std::size_t> g2 = mam::bgchain_detail::mam_bgchain_groups(D, 2);
    CHECK(g2[0] == g2[1]);
    CHECK(g2[2] == g2[3]);
    CHECK(g2[0] != g2[2]);
    CHECK(g2[0] == 0u);  // relabelled by smallest member

    // One group swallows everything; four leave every chain alone.
    const std::vector<std::size_t> g1 = mam::bgchain_detail::mam_bgchain_groups(D, 1);
    for (std::size_t c = 0; c < 4; ++c) CHECK(g1[c] == 0u);
    const std::vector<std::size_t> g4 = mam::bgchain_detail::mam_bgchain_groups(D, 4);
    for (std::size_t c = 0; c < 4; ++c) CHECK(g4[c] == c);

    // The distance is scale-RELATIVE, so a chain with the same profile at ten
    // times the magnitude is not a twin: with two groups it stands apart from
    // the pair it is proportional to.
    std::vector<std::vector<double>> E(2, std::vector<double>(3, 0.0));
    E[0][0] = 1.0;  E[1][0] = 0.2;
    E[0][1] = 1.02; E[1][1] = 0.21;
    E[0][2] = 10.0; E[1][2] = 2.0;
    const std::vector<std::size_t> ge = mam::bgchain_detail::mam_bgchain_groups(E, 2);
    CHECK(ge[0] == ge[1]);
    CHECK(ge[2] != ge[0]);
}

// ---------------------------------------------------------------------------
// The pieces on their own
// ---------------------------------------------------------------------------

TEST_CASE("mam_bgchain_ctmc: a two-station closed chain with no open work is the exact M/M/1//N") {
    // Delay(mu=1) -> Queue(mu=2), N = 3, no open class: cshare is min(e,1) and
    // the chain is the machine-repairman birth-death process, whose stationary
    // law is known in closed form.
    const int N = 3;
    std::vector<int> Nb(1, N);
    std::vector<std::vector<double>> STb(2, std::vector<double>(1, 0.0));
    STb[0][0] = 1.0;   // Delay, mean 1
    STb[1][0] = 0.5;   // Queue, mean 1/2
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;
    std::vector<Matrix<double>> Pb(1, P);
    std::vector<bool> isinf_i(2, false);
    isinf_i[0] = true;
    std::vector<double> nsrv(2, 1.0);
    std::vector<std::vector<double>> cshare(2, std::vector<double>(N + 1, 0.0));
    for (int e = 0; e <= N; ++e) {
        cshare[0][e] = e;              // INF, unused
        cshare[1][e] = std::min(e, 1); // single server, no open work
    }
    const std::vector<std::vector<bool>> supp(2, std::vector<bool>(1, true));
    const mam::BgchainCtmc<double> bg =
        mam::mam_bgchain_ctmc(Nb, STb, Pb, isinf_i, nsrv, cshare, supp, 20000);

    // Birth-death: p(n) at the queue is proportional to prod_{j<n} (N-j)*1 / 2.
    std::vector<double> w(N + 1, 0.0);
    w[0] = 1.0;
    for (int n = 1; n <= N; ++n) w[n] = w[n - 1] * (N - (n - 1)) * 1.0 / 2.0;
    double tot = 0.0, mean = 0.0;
    for (int n = 0; n <= N; ++n) tot += w[n];
    for (int n = 0; n <= N; ++n) mean += n * w[n] / tot;
    CHECK(bg.QLen[1][0] == doctest::Approx(mean).epsilon(1e-9));
    CHECK(bg.QLen[0][0] + bg.QLen[1][0] == doctest::Approx(N).epsilon(1e-9));
    // The queue is busy exactly when it holds a job.
    CHECK(bg.Ubusy[1][0] == doctest::Approx(1.0 - w[0] / tot).epsilon(1e-9));
}

TEST_CASE("mam_bgchain_env: lumping a two-station chain is lossless") {
    // With two stations the closed occupancy at one determines the occupancy at
    // the other, so the level sets are singletons and the lumped generator is
    // the chain itself: the exact-aggregation approximation costs nothing here,
    // which is what makes it the right control on the code path.
    const int N = 2;
    std::vector<int> Nb(1, N);
    std::vector<std::vector<double>> STb(2, std::vector<double>(1, 0.0));
    STb[0][0] = 1.0;
    STb[1][0] = 0.5;
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;
    std::vector<Matrix<double>> Pb(1, P);
    std::vector<bool> isinf_i(2, false);
    isinf_i[0] = true;
    std::vector<double> nsrv(2, 1.0);
    std::vector<std::vector<double>> cshare(2, std::vector<double>(N + 1, 0.0));
    for (int e = 0; e <= N; ++e) {
        cshare[0][e] = e;
        cshare[1][e] = std::min(e, 1);
    }
    const std::vector<std::vector<bool>> supp(2, std::vector<bool>(1, true));
    const mam::BgchainCtmc<double> bg =
        mam::mam_bgchain_ctmc(Nb, STb, Pb, isinf_i, nsrv, cshare, supp, 20000);
    const mam::BgchainEnv<double> env = mam::mam_bgchain_env(bg, 1);

    CHECK(env.esup.size() == static_cast<std::size_t>(N + 1));
    // A generator: every row sums to zero.
    for (std::size_t e = 0; e < env.esup.size(); ++e) {
        double rowsum = 0.0;
        for (std::size_t ep = 0; ep < env.esup.size(); ++ep) rowsum += env.A(e, ep);
        CHECK(rowsum == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    }
    // phi is the marginal of the chain at that station.
    double mean = 0.0;
    for (std::size_t e = 0; e < env.esup.size(); ++e) mean += env.esup[e] * env.phi[e];
    CHECK(mean == doctest::Approx(bg.QLen[1][0]).epsilon(1e-9));
    // The birth rate out of level e is the delay's (N-e) * 1, exactly, because
    // the partition is a bijection here.
    for (std::size_t e = 0; e + 1 < env.esup.size(); ++e)
        CHECK(env.A(e, e + 1) ==
              doctest::Approx(static_cast<double>(N - env.esup[e])).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// The PURELY CLOSED model: the degenerate case of the same construction
// ---------------------------------------------------------------------------

namespace {

/** Delay -> Q1 -> ... -> Qmq -> Delay, one closed class of N jobs. */
qn::Network<double> closed_tandem(const std::string& name, std::size_t mq, double N,
                                  SchedStrategy sched, bool erlang) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("Think");
    std::vector<std::size_t> q;
    for (std::size_t i = 0; i < mq; ++i)
        q.push_back(m.add_queue("Q" + std::to_string(i + 1), sched));
    const std::size_t c = m.add_closed_class("C", N, d);
    m.set_service(d, c, Dd::exp_rate(1.0));
    for (std::size_t i = 0; i < mq; ++i) {
        const double rate = 1.0 + 0.3 * static_cast<double>(i + 1);
        m.set_service(q[i], c, erlang ? Dd::erlang(3.0 * rate, 3) : Dd::exp_rate(rate));
    }
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q[0], 1.0);
    for (std::size_t i = 0; i + 1 < mq; ++i) P.set(c, c, q[i], q[i + 1], 1.0);
    P.set(c, c, q[mq - 1], d, 1.0);
    m.link(P);
    return m;
}

/** Exact MVA on the same model, the independent oracle for a product-form net. */
mva::MvaSolution<double> exact_mva(qn::Network<double>& m) {
    mva::MvaOptions opt;
    opt.method = "exact";
    Matrix<double> init;
    return mva::mva_dispatch(m.get_struct(), opt, init).sol;
}

}  // namespace

TEST_CASE("bgchain: a purely closed tandem is the EXACT closed CTMC") {
    // With no open class cshare never leaves its initial min(e,c), so the
    // background chain alone answers and it answers exactly. The tandem is
    // product-form under PS, so exact MVA is an independent oracle and the
    // agreement is machine precision, not the 1e-4 the mixed cases assert.
    for (std::size_t mq : {std::size_t(2), std::size_t(4), std::size_t(8)}) {
        qn::Network<double> m = closed_tandem("bgClosed", mq, 5.0, SchedStrategy::PS, false);
        const mva::MvaSolution<double> r = run_bg(m);
        qn::Network<double> mref = closed_tandem("ref", mq, 5.0, SchedStrategy::PS, false);
        const mva::MvaSolution<double> ref = exact_mva(mref);
        double held = 0.0;
        for (std::size_t i = 0; i <= mq; ++i) {
            held += r.Q(i, 0);
            CHECK(r.Q(i, 0) == doctest::Approx(ref.Q(i, 0)).epsilon(1e-9));
        }
        CHECK(held == doctest::Approx(5.0).epsilon(1e-9));
    }
}

TEST_CASE("bgchain: a purely closed FCFS tandem with exponential service is the same chain") {
    qn::Network<double> m = closed_tandem("bgClosedFcfs", 3, 4.0, SchedStrategy::FCFS, false);
    const mva::MvaSolution<double> r = run_bg(m);
    qn::Network<double> mref = closed_tandem("ref", 3, 4.0, SchedStrategy::FCFS, false);
    const mva::MvaSolution<double> ref = exact_mva(mref);
    double held = 0.0;
    for (std::size_t i = 0; i <= 3; ++i) {
        held += r.Q(i, 0);
        CHECK(r.Q(i, 0) == doctest::Approx(ref.Q(i, 0)).epsilon(1e-9));
    }
    CHECK(held == doctest::Approx(4.0).epsilon(1e-9));
}

TEST_CASE("the default takes bgchain on a closed model whose chain fits") {
    qn::Network<double> m = closed_tandem("bgClosedDef", 4, 5.0, SchedStrategy::PS, false);
    mam::MamOptions opt;
    opt.method = "default";
    CHECK(mam::solver_mam_run_analyzer(m.get_struct(), opt).actualmethod == "default/bgchain");
}

namespace {

/** Delay -> Q -> Delay for two closed chains of 2 jobs each. */
qn::Network<double> two_chain_cycle(double rate1, double rate2, SchedStrategy sched) {
    qn::Network<double> m("twoChain");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q1", sched);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
    m.set_service(d, c1, Dd::exp_rate(1.0));
    m.set_service(d, c2, Dd::exp_rate(1.0));
    m.set_service(q, c1, Dd::exp_rate(rate1));
    m.set_service(q, c2, Dd::exp_rate(rate2));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the default leaves class-dependent FCFS rates alone") {
    // The chain splits a station's capacity in proportion to the job COUNTS, i.e.
    // service in random order. Under FCFS that is exact only at one shared rate:
    // with Exp(3) against Exp(0.8) it reads 25.2% off SolverCTMC, so the default
    // must stand aside. Under PS, which splits by DEMAND, the same rates cost
    // nothing (measured 0.0e+00), and one shared FCFS rate is exact to 3.5e-16.
    mam::MamOptions opt;
    opt.method = "default";
    qn::Network<double> dep = two_chain_cycle(3.0, 0.8, SchedStrategy::FCFS);
    CHECK(mam::solver_mam_run_analyzer(dep.get_struct(), opt).actualmethod != "default/bgchain");
    qn::Network<double> indep = two_chain_cycle(1.5, 1.5, SchedStrategy::FCFS);
    CHECK(mam::solver_mam_run_analyzer(indep.get_struct(), opt).actualmethod == "default/bgchain");
    qn::Network<double> ps = two_chain_cycle(3.0, 0.8, SchedStrategy::PS);
    CHECK(mam::solver_mam_run_analyzer(ps.get_struct(), opt).actualmethod == "default/bgchain");
}

TEST_CASE("the default leaves a non-exponential closed FCFS model to mna") {
    // The background chain carries the MEAN service time only, which FCFS is not
    // insensitive to, so the mean-only surrogate must not be the DEFAULT here.
    qn::Network<double> m = closed_tandem("bgClosedErl", 2, 4.0, SchedStrategy::FCFS, true);
    mam::MamOptions opt;
    opt.method = "default";
    CHECK(mam::solver_mam_run_analyzer(m.get_struct(), opt).actualmethod != "default/bgchain");
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

TEST_CASE("bgchain refuses a purely OPEN model, which has no closed population") {
    qn::Network<double> m("bgOpen");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("Open");
    m.set_arrival(s, o, Dd::exp_rate(0.5));
    m.set_service(q, o, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    mam::MamOptions opt;
    opt.method = "bgchain";
    CHECK_THROWS_AS(mam::solver_mam_bgchain(m.get_struct(), opt), UnsupportedError);
}

TEST_CASE("bgchain refuses a background chain above bgstates_max") {
    qn::Network<double> m =
        bg_sdq("bgBig", SchedStrategy::PS, 0.4, Dd::exp_rate(2.0), Dd::exp_rate(2.0), 4.0, 1.0);
    mam::MamOptions opt;
    opt.method = "bgchain";
    opt.bgstates_max = 2;  // the chain needs 5 states
    CHECK_THROWS_AS(mam::solver_mam_bgchain(m.get_struct(), opt), UnsupportedError);
}

// ---------------------------------------------------------------------------
// The DEFAULT sizes the chain before choosing it
// ---------------------------------------------------------------------------

TEST_CASE("the default leaves an oversized background chain to dec.source") {
    // `bgchain` by name is refused above bgstates_max, which is right: the user
    // asked for a chain this model cannot afford. `default` must not inherit
    // that refusal -- it has to land on a method that answers. Before the guard
    // (2026-08-14) MATLAB and the JAR both THREW here on mqn_singleserver_fcfs
    // (one closed chain, N = 100 over 4 stations, nchoosek(103,3) = 176851
    // states). bgchain_states sizes the chain up front, so the ladder falls
    // through.
    qn::Network<double> m =
        bg_sdq("bgGuard", SchedStrategy::PS, 0.4, Dd::exp_rate(2.0), Dd::exp_rate(2.0), 4.0, 1.0);

    // The chain this model needs: 4 jobs over the 2 stations the closed class
    // visits, i.e. nchoosek(5,1) = 5 states.
    mam::MamOptions sized;
    CHECK(mam::bgchain_states(m.get_struct(), sized) == doctest::Approx(5.0));

    mam::MamOptions wide;
    wide.method = "default";
    CHECK(mam::solver_mam_run_analyzer(m.get_struct(), wide).actualmethod == "default/bgchain");

    mam::MamOptions narrow;
    narrow.method = "default";
    narrow.bgstates_max = 2;   // below the 5 the chain needs
    CHECK(mam::solver_mam_run_analyzer(m.get_struct(), narrow).actualmethod == "default/dec.source");

    // By name the refusal stands, whatever the default does.
    mam::MamOptions byname;
    byname.method = "bgchain";
    byname.bgstates_max = 2;
    CHECK_THROWS_AS(mam::solver_mam_run_analyzer(m.get_struct(), byname), UnsupportedError);
}

TEST_CASE("bgchain_states counts the AGGREGATED partition, not just the tagged chain") {
    // Three closed chains and one open class, so bgaggr = 1 groups the two
    // untagged chains into one aggregate: B = 2 per pass and the count is the
    // LARGEST pass, not the first. MATLAB, the JAR and python all report 315
    // here and 900 at bgaggr = 2 (mam_bgchain_states.m; the whole AvgTable
    // agrees to 8 digits across the four).
    qn::Network<double> m("bgMulti");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("Open");
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 3.0, d);
    const std::size_t c3 = m.add_closed_class("C3", 4.0, d);
    m.set_arrival(s, o, Dd::exp_rate(0.3));
    m.set_service(q1, o, Dd::exp_rate(2.0));
    m.set_service(q2, o, Dd::exp_rate(3.0));
    m.set_service(d, c1, Dd::exp_rate(1.0));
    m.set_service(q1, c1, Dd::exp_rate(2.0));
    m.set_service(q2, c1, Dd::exp_rate(4.0));
    m.set_service(d, c2, Dd::exp_rate(0.5));
    m.set_service(q1, c2, Dd::exp_rate(1.0));
    m.set_service(q2, c2, Dd::exp_rate(3.0));
    m.set_service(d, c3, Dd::exp_rate(2.0));
    m.set_service(q1, c3, Dd::exp_rate(5.0));
    m.set_service(q2, c3, Dd::exp_rate(6.0));
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q1, 1.0);
    P.set(o, o, q1, q2, 1.0);
    P.set(o, o, q2, k, 1.0);
    const std::size_t cls[3] = {c1, c2, c3};
    for (std::size_t ci = 0; ci < 3; ++ci) {
        P.set(cls[ci], cls[ci], d, q1, 1.0);
        P.set(cls[ci], cls[ci], q1, q2, 1.0);
        P.set(cls[ci], cls[ci], q2, d, 1.0);
    }
    m.link(P);

    mam::MamOptions one;
    CHECK(mam::bgchain_states(m.get_struct(), one) == doctest::Approx(315.0));
    mam::MamOptions two;
    two.bgaggr = 2;
    CHECK(mam::bgchain_states(m.get_struct(), two) == doctest::Approx(900.0));

    // and the default still takes bgchain, since 315 is well under bgstates_max
    mam::MamOptions run;
    run.method = "default";
    const mva::AvgResult<double> r = mam::solver_mam_run_analyzer(m.get_struct(), run);
    CHECK(r.actualmethod == "default/bgchain");
    CHECK(r.QN(1, 1) == doctest::Approx(0.38964452).epsilon(1e-6));   // Delay, C1
    CHECK(r.QN(2, 3) == doctest::Approx(2.5139685).epsilon(1e-6));    // Queue1, C3
}
