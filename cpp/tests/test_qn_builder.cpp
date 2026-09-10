/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The general Network builder and the refresh chain it feeds.
 *
 * Until now the only model this port could build was a SolverLN layer, whose
 * shape is fixed by the layer construction; these are the first models built
 * from the constructor API (`Queue`, `Delay`, `Source`, `ClassSwitch`, `link`)
 * and refreshed through `NetworkStruct::refresh_struct`.
 *
 * THE ORACLES ARE INDEPENDENT OF THE PORT, deliberately: an M/M/1 in closed
 * form, an exact MVA recursion carried out by hand on a two-station closed
 * model, and the stationary vector of the class-switch chain solved by hand.
 * A regression against numbers this code produced would only pin its behaviour
 * in place, not check it.
 *
 * Every value asserted here was afterwards confirmed against MATLAB, on the
 * same three models built through the MATLAB constructor API and solved by
 * SolverMVA: struct fields (nstations, nclasses, nchains, rates, visits, cap,
 * classcap, droprule) and the AvgTable agree to the digits printed.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::DropStrategy;
using lang::RoutingStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Source -> Queue -> Sink, one open class: the M/M/1 of the manual. */
TEST_CASE("the builder produces an M/M/1 whose struct and solution are the closed form") {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("Class1");
    m.set_arrival(src, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nstations == 2);   // the Source is a station, the Sink is not
    CHECK(sn.nof_nodes() == 3);
    CHECK(sn.nclasses == 1);
    CHECK(sn.nchains == 1);
    CHECK(sn.rates(0, 0) == doctest::Approx(1.0));  // arrival rate at the Source
    CHECK(sn.rates(1, 0) == doctest::Approx(2.0));
    CHECK(sn.visits[0](1, 0) == doctest::Approx(1.0));
    CHECK(std::isinf(sn.classes[0].population));

    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::MvaSolution<double> r = mva::solver_mva_analyzer(sn, opt, init);
    // rho = 1/2: Q = rho/(1-rho) = 1, R = 1/(mu-lambda) = 1, U = rho, X = lambda
    CHECK(r.Q(1, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(r.R(1, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(r.U(1, 0) == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(r.X[0] == doctest::Approx(1.0).epsilon(1e-6));
}

/**
 * Delay(Z=1) + FCFS Queue(D=0.5), one closed class of 3 jobs.
 *
 * The exact MVA recursion, carried out by hand:
 *   n=1  R=0.5      X=1/1.5     =0.666667  Q=0.333333
 *   n=2  R=0.666667 X=2/1.666667=1.2       Q=0.8
 *   n=3  R=0.9      X=3/1.9     =1.578947  Q=1.421053
 * so U = X D = 0.789474 and the delay holds X Z = 1.578947 jobs.
 */
TEST_CASE("a closed two-station model reproduces the exact MVA recursion") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nstations == 2);
    CHECK(sn.nchains == 1);
    CHECK(sn.visits[0](0, 0) == doctest::Approx(1.0));
    CHECK(sn.visits[0](1, 0) == doctest::Approx(1.0));
    // no explicit capacity: the bound is the chain population, at BOTH
    // stations, which is what MATLAB reports here (cap = [3 3])
    CHECK(sn.cap[0] == doctest::Approx(3.0));
    CHECK(sn.cap[1] == doctest::Approx(3.0));

    mva::MvaOptions opt;
    opt.method = "exact";
    Matrix<double> init;
    const mva::MvaSolution<double> r = mva::solver_mva_analyzer(sn, opt, init);
    CHECK(r.method == "exact");
    CHECK(r.X[0] == doctest::Approx(1.578947368).epsilon(1e-9));
    CHECK(r.Q(1, 0) == doctest::Approx(1.421052632).epsilon(1e-9));
    CHECK(r.Q(0, 0) == doctest::Approx(1.578947368).epsilon(1e-9));
    CHECK(r.R(1, 0) == doctest::Approx(0.9).epsilon(1e-9));
    CHECK(r.U(1, 0) == doctest::Approx(0.789473684).epsilon(1e-9));
}

/**
 * Delay -> Queue -> ClassSwitch -> Delay, two classes switching by
 *   C = [0.3 0.7; 0.6 0.4].
 *
 * The two classes are one chain, and the class process on a cycle is the chain
 * C itself: pi C = pi gives 0.7 pi1 = 0.6 pi2, so pi = (6/13, 7/13). The visits
 * are normalised at the reference station, whose row must therefore sum to 1.
 */
TEST_CASE("a ClassSwitch node merges the classes into one chain and splits the visits") {
    qn::Network<double> m("cs");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    Matrix<double> C(2, 2, 0.0);
    C(0, 0) = 0.3;
    C(0, 1) = 0.7;
    C(1, 0) = 0.6;
    C(1, 1) = 0.4;
    const std::size_t cs = m.add_class_switch("CS", C);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(2.0));
    m.set_service(q, c1, D::exp_rate(4.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, cs, 1.0);
    P.set(c1, c1, cs, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, cs, 1.0);
    P.set(c2, c2, cs, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nclasses == 2);
    CHECK(sn.nchains == 1);  // the switch couples them
    const double p1 = 6.0 / 13.0, p2 = 7.0 / 13.0;
    CHECK(sn.visits[0](0, 0) == doctest::Approx(p1).epsilon(1e-10));
    CHECK(sn.visits[0](0, 1) == doctest::Approx(p2).epsilon(1e-10));
    CHECK(sn.visits[0](1, 0) == doctest::Approx(p1).epsilon(1e-10));
    CHECK(sn.visits[0](1, 1) == doctest::Approx(p2).epsilon(1e-10));
    // the ClassSwitch is not stateful, so it carries no visits row of its own
    CHECK(sn.nof_stateful() == 2);
}

TEST_CASE("RAND routing is expanded to the uniform split over the connected nodes") {
    qn::Network<double> m("rand");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    // the weights given are connections only: RAND overrides them
    P.set(d, q1, 0.9);
    P.set(d, q2, 0.1);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    m.set_routing(d, c, RoutingStrategy::RAND);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.route_eff(1, 1, d, q1) == doctest::Approx(0.5));
    CHECK(sn.route_eff(1, 1, d, q2) == doctest::Approx(0.5));
    CHECK(sn.visits[0](1, 0) == doctest::Approx(0.5));
    CHECK(sn.visits[0](2, 0) == doctest::Approx(0.5));
}

TEST_CASE("a finite buffer derives DROP for an open class and the capacity fields") {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("O1");
    m.set_arrival(src, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_capacity(q, 3.0);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t iq = sn.nodes[q - 1].station;
    CHECK(sn.cap[iq - 1] == doctest::Approx(3.0));
    CHECK(sn.classcap[iq - 1][0] == doctest::Approx(3.0));
    CHECK(sn.droprule[iq - 1][0] == DropStrategy::DROP);
    // the Source is unbounded and its rule is never consulted
    CHECK(std::isinf(sn.cap[sn.sourceIdx - 1]));
}

TEST_CASE("setNumberOfServers is a no-op on an inf-scheduled station, as in MATLAB") {
    qn::Network<double> m("inf");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 1.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.0));
    m.set_number_of_servers(d, 4.0);  // ignored: the Delay stays infinite
    m.set_number_of_servers(q, 4.0);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(std::isinf(sn.stations[sn.nodes[d - 1].station - 1].nservers));
    CHECK(sn.stations[sn.nodes[q - 1].station - 1].nservers == doctest::Approx(4.0));
}

TEST_CASE("the builder and the refresh refuse by name what they cannot represent") {
    SUBCASE("a routing strategy with no routing-matrix expansion") {
        // FIRING is the SPN transition rule and has no expansion into `rt`, so
        // it is refused by name. JSQ and SQ are NOT: `getRoutingMatrix.m:117`
        // spreads them uniformly over the connected destinations and
        // refresh_routing does the same, so refusing them would refuse a model
        // the reference solves. The refusal is about the expansion existing,
        // not about the dispatcher being stateless.
        qn::Network<double> m("firing");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t c = m.add_closed_class("C1", 1.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        m.set_routing(d, c, RoutingStrategy::FIRING);
        CHECK_THROWS_AS(m.get_struct(), UnsupportedError);
    }
    SUBCASE("a state-dependent dispatcher the reference expands uniformly") {
        qn::Network<double> m("jsq");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
        const std::size_t c = m.add_closed_class("C1", 1.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q1, c, D::exp_rate(1.0));
        m.set_service(q2, c, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q1, 1.0);
        P.set(d, q2, 1.0);
        P.set(q1, d, 1.0);
        P.set(q2, d, 1.0);
        m.link(P);
        m.set_routing(d, c, RoutingStrategy::JSQ);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        // Half each, as getRoutingMatrix.m:117 gives it. The determinism the
        // dispatcher adds is a higher moment the routing matrix cannot carry.
        CHECK(sn.route_eff(c, c, d, q1) == doctest::Approx(0.5));
        CHECK(sn.route_eff(c, c, d, q2) == doctest::Approx(0.5));
    }
    SUBCASE("a class with no service anywhere") {
        qn::Network<double> m("nosvc");
        const std::size_t d = m.add_delay("Delay");
        m.add_closed_class("C1", 1.0, d);
        CHECK_THROWS_AS(m.get_struct(), InputError);
    }
    SUBCASE("a Fork with no Join") {
        qn::Network<double> m("nojoin");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t f = m.add_fork("Fork");
        const std::size_t c = m.add_closed_class("C1", 1.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(d, f, 1.0);
        P.set(f, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(m.get_struct(), InputError);
    }
    SUBCASE("a setting that names a node which serves no jobs") {
        qn::Network<double> m("sink");
        m.add_source("Source");
        const std::size_t snk = m.add_sink("Sink");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
        const std::size_t c = m.add_open_class("O1");
        m.set_service(q, c, D::exp_rate(1.0));
        CHECK_THROWS_AS(m.set_service(snk, c, D::exp_rate(1.0)), InputError);
    }
}

TEST_CASE("the distribution descriptor carries the moments and the representation") {
    const D e = D::exp_mean(2.0);
    CHECK(e.mean == doctest::Approx(2.0));
    CHECK(e.scv == doctest::Approx(1.0));
    CHECK(e.phases() == 1);

    const D er = D::erlang(3.0, 3);
    CHECK(er.mean == doctest::Approx(1.0));
    CHECK(er.scv == doctest::Approx(1.0 / 3.0));
    CHECK(er.phases() == 3);
    CHECK(er.mu_vec()[0] == doctest::Approx(3.0));
    CHECK(er.phi_vec()[0] == doctest::Approx(0.0));
    CHECK(er.phi_vec()[2] == doctest::Approx(1.0));

    const D h = D::hyperexp(0.3, 1.0, 4.0);
    CHECK(h.mean == doctest::Approx(0.3 * 1.0 + 0.7 * 0.25));
    // E[X^2] = 2 (0.3/1 + 0.7/16) = 0.6875, SCV = 0.6875/0.475^2 - 1
    CHECK(h.scv == doctest::Approx(0.6875 / (0.475 * 0.475) - 1.0));

    const D c2 = D::cox2(2.0, 5.0, 0.4);
    CHECK(c2.mean == doctest::Approx(0.5 + 0.6 * 0.2));

    const D p = D::pareto(3.0, 1.0);
    CHECK(p.mean == doctest::Approx(1.5));
    CHECK(p.scv == doctest::Approx(1.0 / 3.0));

    const D u = D::uniform(1.0, 3.0);
    CHECK(u.mean == doctest::Approx(2.0));
    CHECK(u.scv == doctest::Approx((4.0 / 12.0) / 4.0));
}

}  // namespace

TEST_CASE("network: an open class reached only by CLASS SWITCHING is accepted") {
    // `cs_multi_diamond`: Class1 carries the only Source arrival, Class2 and
    // Class3 are fed by routing blocks with r != s. MATLAB and python both
    // solve it; the C++ validator used to refuse it, because it asked "has a
    // Source arrival" where the question is "can a job ever enter this class".
    qn::Network<double> m("cs_diamond");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", lang::SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Class1");
    const std::size_t c2 = m.add_open_class("Class2");
    const std::size_t c3 = m.add_open_class("Class3");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(1.0));
    m.set_service(q1, c1, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q2, c2, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q3, c3, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q1, 1.0);
    P.set(c1, c2, q1, q2, 0.5);   // Class1 -> Class2 by class switching
    P.set(c1, c3, q1, q3, 0.5);   // Class1 -> Class3 by class switching
    P.set(c2, c2, q2, k, 1.0);
    P.set(c3, c3, q3, k, 1.0);
    m.link(P);
    CHECK_NOTHROW(m.get_struct());
}

TEST_CASE("network: an UNREACHABLE open class is ACCEPTED and reports zero") {
    // THIS MODEL IS `gallery_erlerl1`, which ships in all three reference
    // suites: an open class whose Source arrival is Disabled, which is served
    // at a Queue and routed on to the Sink, and which nothing ever switches
    // into. No job can enter it.
    //
    // Until 2026-08-05 this port REFUSED it, on the reasoning that a class no
    // job can enter is a modelling error. MATLAB and native Python both SOLVE
    // it and the class simply carries zero of every metric. The refusal made
    // this port the only one that could not read its own gallery, which is what
    // the parity-static cpp row measured. What is asserted here is therefore
    // the acceptance AND the zero: accepting the model and then reporting a
    // nonzero row for a class with no arrivals would be the worse defect, and
    // only the second half catches it.
    //
    // The flow-sink guard added on 2026-08-31 broke this again in all four
    // codebases, and for a reason worth keeping: it read the node visits of the
    // Orphan CHAIN, which has no inflow to normalise on and so came back
    // carrying Fed's route -- Source = Q1 = Sink = 1 -- and accused Q1 of not
    // serving a class never routed there. An unfed chain is now skipped; see
    // _kb/07-cross-language-parity.md.
    qn::Network<double> m("orphan");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Fed");
    const std::size_t c2 = m.add_open_class("Orphan");
    m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(1.0));
    m.set_arrival(src, c2, lang::Distrib<double>::disabled_dist());
    m.set_service(q1, c1, lang::Distrib<double>::exp_rate(2.0));
    m.set_service(q2, c2, lang::Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q1, 1.0);
    P.set(c1, c1, q1, k, 1.0);
    P.set(c2, c2, q2, k, 1.0);   // Orphan circulates but nothing ever creates it
    m.link(P);
    CHECK_NOTHROW(m.get_struct());
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nclasses == 2);

    mva::MvaOptions o;
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(sn, o, init);
    // Stations: 0 Source, 1 Q1, 2 Q2. The fed class is the plain M/M/1 with
    // rho = 1/2; the orphan carries nothing anywhere.
    CHECK(r.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r.QN(2, 1) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(r.TN(2, 1) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(r.UN(2, 1) == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("network: a class-switched open model SOLVES and its answer MOVES") {
    // Passing validation is not the bar. A validator that accepts a model it
    // cannot actually solve is the same family of defect as a reader that
    // accepts a key it does not consume, so the acceptance side is held to the
    // same evidence: the model must SOLVE, and the metrics of the switched
    // classes must RESPOND to the switching probabilities. If Class2 and Class3
    // were merely tolerated and carried no traffic, their throughputs would be
    // zero and would not move with the split.
    auto build = [](double p2) {
        qn::Network<double> m("cs_diamond");
        const std::size_t src = m.add_source("Source");
        const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
        const std::size_t q3 = m.add_queue("Q3", lang::SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c1 = m.add_open_class("Class1");
        const std::size_t c2 = m.add_open_class("Class2");
        const std::size_t c3 = m.add_open_class("Class3");
        m.set_arrival(src, c1, lang::Distrib<double>::exp_rate(1.0));
        m.set_service(q1, c1, lang::Distrib<double>::exp_rate(4.0));
        m.set_service(q2, c2, lang::Distrib<double>::exp_rate(4.0));
        m.set_service(q3, c3, lang::Distrib<double>::exp_rate(4.0));
        qn::RoutingMatrix<double> P;
        P.set(c1, c1, src, q1, 1.0);
        P.set(c1, c2, q1, q2, p2);
        P.set(c1, c3, q1, q3, 1.0 - p2);
        P.set(c2, c2, q2, k, 1.0);
        P.set(c3, c3, q3, k, 1.0);
        m.link(P);
        return m;
    };
    auto tput = [](qn::Network<double>& m, std::size_t st, std::size_t cls) {
        mva::MvaOptions o;
        Matrix<double> init;
        return mva::solver_mva_run_analyzer(m.get_struct(), o, init).TN(st, cls);
    };
    qn::Network<double> a = build(0.25), b = build(0.75);
    // Stations: 0 Source, 1 Q1, 2 Q2, 3 Q3. Class indices 0,1,2.
    const double a2 = tput(a, 2, 1), a3 = tput(a, 3, 2);
    const double b2 = tput(b, 2, 1), b3 = tput(b, 3, 2);

    // It SOLVES: the switched classes carry the split of the arrival rate, so
    // their throughputs are the exact routing probabilities times lambda = 1.
    CHECK(a2 == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(a3 == doctest::Approx(0.75).epsilon(1e-9));
    CHECK(b2 == doctest::Approx(0.75).epsilon(1e-9));
    CHECK(b3 == doctest::Approx(0.25).epsilon(1e-9));
    // And the answer MOVES with the class-switching route, in the right
    // direction, which a merely-tolerated class could not do.
    CHECK(b2 > a2);
    CHECK(b3 < a3);
    // Flow is conserved across the switch: everything that arrives leaves.
    CHECK(a2 + a3 == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(b2 + b3 == doctest::Approx(1.0).epsilon(1e-9));
}
