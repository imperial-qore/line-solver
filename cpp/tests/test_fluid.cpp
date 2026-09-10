/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverFluid, the `closing` method.
 *
 * Every expected number below is MATLAB's
 * `SolverFluid(model,'method','closing').getAvgTable` on the same model, so
 * these are parity assertions and not a record of what this port happens to
 * produce. The fluid fixed point is the root of a smooth drift, so agreement
 * is to the integrator tolerance rather than to machine precision; 1e-4 is
 * what the reference's own `tol` admits.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/fluid_export_odes.h"
#include "line/solvers/fluid/fluid_jacobian.h"
#include "line/solvers/fluid/fluid_passage.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/fluid_symodes.h"
#include "line/solvers/fluid/solver_fluid.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

TEST_CASE("fluid: a closed single-class network matches the MATLAB getAvgTable") {
    // Delay(think 1) -> PS Queue(rate 2), N = 4. The queue saturates, so the
    // fluid limit pins its utilization at exactly 1 and the throughput at the
    // service rate -- a case where the answer is known independently.
    qn::Network<double> m("fluid_cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), fluid::FluidOptions());
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), iq = station_of(sn, "Queue1");

    CHECK(s.QN(id, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.UN(id, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.RN(id, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.TN(id, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.QN(iq, 0) == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(s.UN(iq, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.RN(iq, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.TN(iq, 0) == doctest::Approx(2.0).epsilon(1e-4));

    // The population is conserved by the drift, so it is a check the reference
    // values cannot influence.
    CHECK(s.QN(id, 0) + s.QN(iq, 0) == doctest::Approx(4.0).epsilon(1e-6));
}

TEST_CASE("fluid: a multiclass phase-type network matches the MATLAB getAvgTable") {
    // Two closed classes routed in OPPOSITE directions over three stations,
    // with Erlang service at two of them and a two-server FCFS queue. This is
    // the case that exercises the parts a single-class exponential model
    // cannot: multiple phases per (station, class), the processor-sharing
    // scaling once a station is over its server count, and class-dependent
    // routing.
    qn::Network<double> m("fl2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 5, d);
    const std::size_t c2 = m.add_closed_class("C2", 3, d);
    m.set_service(d, c1, D::exp_rate(2.0));
    m.set_service(d, c2, D::erlang(2.0, 2));   // mean 1, order 2
    m.set_service(q1, c1, D::erlang(6.0, 3));  // mean 0.5, order 3
    m.set_service(q1, c2, D::exp_rate(1.5));
    m.set_service(q2, c1, D::exp_rate(1.2));
    m.set_service(q2, c2, D::exp_rate(0.8));
    m.set_number_of_servers(q2, 2);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, d, 1.0);
    P.set(c2, c2, d, q2, 1.0);
    P.set(c2, c2, q2, q1, 1.0);
    P.set(c2, c2, q1, d, 1.0);
    m.link(P);

    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), fluid::FluidOptions());
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), i1 = station_of(sn, "Q1"),
                      i2 = station_of(sn, "Q2");

    // MATLAB SolverFluid(model,'method','closing').getAvgTable
    CHECK(s.QN(id, 0) == doctest::Approx(0.64527).epsilon(1e-4));
    CHECK(s.UN(id, 0) == doctest::Approx(0.64527).epsilon(1e-4));
    CHECK(s.RN(id, 0) == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(s.TN(id, 0) == doctest::Approx(1.29053).epsilon(1e-4));
    CHECK(s.QN(id, 1) == doctest::Approx(0.53210).epsilon(1e-4));
    CHECK(s.RN(id, 1) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.TN(id, 1) == doctest::Approx(0.53210).epsilon(1e-4));

    CHECK(s.QN(i1, 0) == doctest::Approx(3.27929).epsilon(1e-4));
    CHECK(s.UN(i1, 0) == doctest::Approx(0.64527).epsilon(1e-4));
    CHECK(s.RN(i1, 0) == doctest::Approx(2.54103).epsilon(1e-4));
    CHECK(s.QN(i1, 1) == doctest::Approx(1.80278).epsilon(1e-4));
    CHECK(s.UN(i1, 1) == doctest::Approx(0.35473).epsilon(1e-4));
    CHECK(s.RN(i1, 1) == doctest::Approx(3.38804).epsilon(1e-4));

    CHECK(s.QN(i2, 0) == doctest::Approx(1.07545).epsilon(1e-4));
    CHECK(s.UN(i2, 0) == doctest::Approx(0.53772).epsilon(1e-4));  // two servers
    CHECK(s.RN(i2, 0) == doctest::Approx(0.83333).epsilon(1e-4));
    CHECK(s.QN(i2, 1) == doctest::Approx(0.66512).epsilon(1e-4));
    CHECK(s.UN(i2, 1) == doctest::Approx(0.33256).epsilon(1e-4));
    CHECK(s.RN(i2, 1) == doctest::Approx(1.25).epsilon(1e-4));

    // Both closed populations are conserved.
    CHECK(s.QN(id, 0) + s.QN(i1, 0) + s.QN(i2, 0) == doctest::Approx(5.0).epsilon(1e-5));
    CHECK(s.QN(id, 1) + s.QN(i1, 1) + s.QN(i2, 1) == doctest::Approx(3.0).epsilon(1e-5));
}

TEST_CASE("fluid: an open network matches the MATLAB getAvgTable and reports no Source queue") {
    // Source(0.6) -> PS Q1 -> PS Q2(Erlang) -> Sink. In the fluid limit an
    // under-loaded PS queue has utilization exactly rho, so Q1 and Q2 are
    // known in closed form as well: 0.6/1.5 and 0.6*0.5.
    qn::Network<double> m("fl3");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    m.set_arrival(src, c1, D::exp_rate(0.6));
    m.set_service(q1, c1, D::exp_rate(1.5));
    m.set_service(q2, c1, D::erlang(4.0, 2));  // mean 0.5, order 2
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, snk, 1.0);
    m.link(P);

    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), fluid::FluidOptions());
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t is = station_of(sn, "Source"), i1 = station_of(sn, "Q1"),
                      i2 = station_of(sn, "Q2");

    CHECK(s.QN(i1, 0) == doctest::Approx(0.4).epsilon(1e-4));
    CHECK(s.UN(i1, 0) == doctest::Approx(0.4).epsilon(1e-4));
    CHECK(s.RN(i1, 0) == doctest::Approx(0.66667).epsilon(1e-4));
    CHECK(s.TN(i1, 0) == doctest::Approx(0.6).epsilon(1e-4));
    CHECK(s.QN(i2, 0) == doctest::Approx(0.3).epsilon(1e-4));
    CHECK(s.UN(i2, 0) == doctest::Approx(0.3).epsilon(1e-4));
    CHECK(s.RN(i2, 0) == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(s.TN(i2, 0) == doctest::Approx(0.6).epsilon(1e-4));

    // The EXT source carries unit mass internally to drive the arrivals, but
    // that is not a queue: the reference reports 0 there and keeps only the
    // arrival rate.
    CHECK(s.QN(is, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.UN(is, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.RN(is, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.TN(is, 0) == doctest::Approx(0.6).epsilon(1e-4));
}

TEST_CASE("fluid: the state layout skips a station-class pair with no service") {
    // Two classes but only one served at each queue: the layout must allocate
    // NO phases for the unserved pairs, since a zero-length block is what
    // keeps the drift's indices aligned with the reference's.
    qn::Network<double> m("layout");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2, d);
    const std::size_t c2 = m.add_closed_class("C2", 2, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(1.0));
    m.set_service(q, c1, D::erlang(4.0, 2));
    m.set_service(q, c2, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    // Delay: one phase per class; Q: two phases for C1, one for C2.
    CHECK(L.nstates == 5);
    const std::size_t iq = station_of(sn, "Q");
    CHECK(L.kic[iq][0] == 2);
    CHECK(L.kic[iq][1] == 1);
    CHECK(L.enabled[iq][0]);
    CHECK(L.enabled[iq][1]);
}

TEST_CASE("fluid: every documented method runs, and an unknown name is refused") {
    qn::Network<double> m("refuse");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    // Everything the solver advertises must actually run on a plain closed
    // model -- a listed name that errors would make the list a false claim.
    for (const char* good : {"default", "closing", "matrix", "pnorm", "statedep", "softmin",
                             "tbi", "diffusion", "fluid.closing", "fluid.matrix"}) {
        o.method = good;
        CHECK_NOTHROW(fluid::solver_fluid(m.get_struct(), o));
    }
    // mfq needs an open single-queue model. On one that is not, it SUBSTITUTES
    // rather than refuses, which is what the reference does:
    // solver_fluid_analyzer.m:145 warns "MFQ not applicable: ... Falling back to
    // matrix method" and re-enters solver_fluid_matrix. Refusing here rejected
    // every multi-station model MATLAB and native python both answer, and took
    // the aliases 'butools' and 'aoi' with it. This port has no line_warning
    // channel, so the substitution is visible in `method` instead.
    o.method = "mfq";
    fluid::FluidSolution fb;
    CHECK_NOTHROW(fb = fluid::solver_fluid(m.get_struct(), o));
    CHECK(fb.method == "matrix");
    // A name no method answers to is refused rather than silently defaulted.
    for (const char* bad : {"nosuchmethod", "amva", "exact"}) {
        o.method = bad;
        CHECK_THROWS_AS(fluid::solver_fluid(m.get_struct(), o), UnsupportedError);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// The other fluid methods. Expected values are MATLAB
// SolverFluid(model,'method',<name>).getAvgTable on the same model.
// ---------------------------------------------------------------------------
namespace {

/** The two-class Erlang/PS/2-server-FCFS model the closing tests also use. */
qn::Network<double> build_fl2() {
    qn::Network<double> m("fl2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 5, d);
    const std::size_t c2 = m.add_closed_class("C2", 3, d);
    m.set_service(d, c1, D::exp_rate(2.0));
    m.set_service(d, c2, D::erlang(2.0, 2));
    m.set_service(q1, c1, D::erlang(6.0, 3));
    m.set_service(q1, c2, D::exp_rate(1.5));
    m.set_service(q2, c1, D::exp_rate(1.2));
    m.set_service(q2, c2, D::exp_rate(0.8));
    m.set_number_of_servers(q2, 2);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, d, 1.0);
    P.set(c2, c2, d, q2, 1.0);
    P.set(c2, c2, q2, q1, 1.0);
    P.set(c2, c2, q1, d, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("fluid matrix: the default method matches the MATLAB getAvgTable") {
    // MATLAB routes default, matrix AND pnorm to solver_fluid_matrix, so this
    // is what `default` must reproduce.
    qn::Network<double> m = build_fl2();
    fluid::FluidOptions o;
    o.method = "matrix";
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), i1 = station_of(sn, "Q1"),
                      i2 = station_of(sn, "Q2");
    CHECK(s.QN(id, 0) == doctest::Approx(0.64527).epsilon(1e-4));
    CHECK(s.UN(id, 0) == doctest::Approx(0.64527).epsilon(1e-4));
    CHECK(s.QN(i1, 0) == doctest::Approx(3.27929).epsilon(1e-4));
    CHECK(s.QN(i1, 1) == doctest::Approx(1.80278).epsilon(1e-4));
    CHECK(s.QN(i2, 0) == doctest::Approx(1.07545).epsilon(1e-3));
    CHECK(s.UN(i2, 0) == doctest::Approx(0.53772).epsilon(1e-3));
    CHECK(s.TN(i2, 1) == doctest::Approx(0.53210).epsilon(1e-3));
    CHECK(s.method == "matrix");

    // `default` must select exactly this method.
    fluid::FluidOptions od;
    const fluid::FluidSolution sd = fluid::solver_fluid(m.get_struct(), od);
    CHECK(sd.method == "matrix");
    CHECK(sd.QN(i1, 0) == doctest::Approx(s.QN(i1, 0)).epsilon(1e-9));
}

TEST_CASE("fluid statedep and softmin match the MATLAB getAvgTable") {
    qn::Network<double> m = build_fl2();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), i1 = station_of(sn, "Q1"),
                      i2 = station_of(sn, "Q2");

    fluid::FluidOptions o;
    o.method = "statedep";
    const fluid::FluidSolution sd = fluid::solver_fluid(sn, o);
    CHECK(sd.QN(id, 0) == doctest::Approx(0.62555).epsilon(1e-4));
    CHECK(sd.QN(i1, 0) == doctest::Approx(3.17035).epsilon(1e-4));
    CHECK(sd.QN(i2, 0) == doctest::Approx(1.20410).epsilon(1e-4));
    // The FCFS station shares its servers by MEAN SERVICE TIME under statedep,
    // which is what these three pin -- they are wrong under any other rule.
    CHECK(sd.UN(i2, 0) == doctest::Approx(0.52129).epsilon(1e-4));
    CHECK(sd.RN(i2, 0) == doctest::Approx(0.96243).epsilon(1e-4));
    CHECK(sd.TN(i2, 0) == doctest::Approx(1.25110).epsilon(1e-4));

    o.method = "softmin";
    const fluid::FluidSolution ss = fluid::solver_fluid(sn, o);
    CHECK(ss.QN(id, 0) == doctest::Approx(0.62556).epsilon(1e-4));
    CHECK(ss.QN(i1, 0) == doctest::Approx(3.17135).epsilon(1e-4));
    CHECK(ss.QN(i2, 0) == doctest::Approx(1.20308).epsilon(1e-4));
    CHECK(ss.UN(i2, 0) == doctest::Approx(0.60154).epsilon(1e-4));
    CHECK(ss.TN(i2, 0) == doctest::Approx(1.44370).epsilon(1e-4));

    // Smoothing the kink moves the fixed point only slightly.
    CHECK(std::fabs(ss.QN(i1, 0) - sd.QN(i1, 0)) < 1e-2);
}

TEST_CASE("fluid tbi reproduces closing when the partition is a single cell") {
    // Three stations and cellsize 5, so tbi_partition yields ONE cell: there
    // are no external events and the cell drift IS the closing drift. This is
    // the case where TBI has an independent right answer.
    qn::Network<double> m = build_fl2();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(fluid::tbi_partition(sn, 5).size() == 1);

    fluid::FluidOptions oc;
    oc.method = "closing";
    const fluid::FluidSolution sc = fluid::solver_fluid(sn, oc);
    fluid::FluidOptions ot;
    ot.method = "tbi";
    const fluid::FluidSolution st = fluid::solver_fluid(sn, ot);
    CHECK(st.method == "tbi");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            CHECK(st.QN(i, r) == doctest::Approx(sc.QN(i, r)).epsilon(1e-4));
}

TEST_CASE("fluid tbi partitions a wide model into several cells") {
    // Ten queues off one delay: with cellsize 5 the greedy merge must produce
    // ceil(11/5) = 3 cells, which is what makes the decomposition non-trivial.
    qn::Network<double> m("wide");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    std::vector<std::size_t> qs;
    for (int i = 0; i < 10; ++i) qs.push_back(m.add_queue("Q" + std::to_string(i),
                                                          lang::SchedStrategy::PS));
    for (std::size_t q : qs) m.set_service(q, c, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t q : qs) {
        P.set(c, c, d, q, 1.0 / static_cast<double>(qs.size()));
        P.set(c, c, q, d, 1.0);
    }
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<std::vector<std::size_t>> cells = fluid::tbi_partition(sn, 5);
    CHECK(cells.size() == 3);
    // The partition must cover every station exactly once.
    std::vector<std::size_t> seen;
    for (const std::vector<std::size_t>& cl : cells)
        for (std::size_t i : cl) seen.push_back(i);
    std::sort(seen.begin(), seen.end());
    CHECK(seen.size() == sn.nstations);
    for (std::size_t i = 0; i < seen.size(); ++i) CHECK(seen[i] == i);

    // Mass is only approximately conserved once the cells are solved with each
    // other's flow FROZEN: until the sweeps converge, what one cell sends need
    // not equal what another receives. The reference expects this too -- it
    // carries an explicit "TBI mass conservation gap" warning -- so the check
    // is that the decomposition stays in the right neighbourhood, not that it
    // conserves exactly the way a single integration does.
    fluid::FluidOptions ot;
    ot.method = "tbi";
    const fluid::FluidSolution st = fluid::solver_fluid(sn, ot);
    double tot = 0.0;
    for (std::size_t i = 0; i < sn.nstations; ++i) tot += st.QN(i, 0);
    // At the default grid the decomposition holds the population to well under
    // a percent; see the table in fluid_tbi.h for how this tightens with grid.
    CHECK(tot == doctest::Approx(4.0).epsilon(5e-3));
}

TEST_CASE("fluid diffusion conserves the population and refuses what it cannot model") {
    qn::Network<double> m("diff");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "diffusion";
    o.iter_max = 2000;
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    CHECK(s.method == "diffusion");
    // The trajectory is renormalised every step, so the population is exact
    // even though the path itself is random.
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(s.QN(station_of(sn, "Delay"), 0) + s.QN(station_of(sn, "Q"), 0) ==
          doctest::Approx(4.0).epsilon(1e-9));
    for (std::size_t i = 0; i < sn.nstations; ++i) CHECK(s.QN(i, 0) >= 0.0);

    // A two-server station is outside the method's stated scope.
    qn::Network<double> m2 = build_fl2();
    CHECK_THROWS_AS(fluid::solver_fluid(m2.get_struct(), o), UnsupportedError);
}

TEST_CASE("fluid mfq matches the MATLAB getAvgTable on modulated and unmodulated queues") {
    qn::Network<double> m("mfq");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(0.5));
    m.set_service(q, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "mfq";
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t iq = station_of(sn, "Q");
    // rho = 0.5: L = rho/(1-rho) = 1, W = 1/(mu-lambda) = 2. Exact, not fluid.
    CHECK(s.QN(iq, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(s.RN(iq, 0) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(s.UN(iq, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(s.TN(iq, 0) == doctest::Approx(0.5).epsilon(1e-9));

    // A phase-type service makes the drain rate Markov-modulated, so the
    // fluid-fluid solver runs. THE ANSWER IS NOT THE M/G/1 ONE and must not be
    // checked against it: this is a fluid queue, where work arrives as a
    // continuous stream rather than as customers, so the sojourn is a drop's
    // delay. M/G/1 would say 0.625 here; MATLAB's mfq says 0.041667, and that
    // is the quantity this method is defined to compute.
    qn::Network<double> m2("mfq2");
    const std::size_t s2 = m2.add_source("Source");
    const std::size_t q2 = m2.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t k2 = m2.add_sink("Sink");
    const std::size_t c2 = m2.add_open_class("C1");
    m2.set_arrival(s2, c2, D::exp_rate(0.5));
    m2.set_service(q2, c2, D::erlang(4.0, 2));  // mean 0.5, order 2
    qn::RoutingMatrix<double> P2;
    P2.set(c2, c2, s2, q2, 1.0);
    P2.set(c2, c2, q2, k2, 1.0);
    m2.link(P2);
    const fluid::FluidSolution s2r = fluid::solver_fluid(m2.get_struct(), o);
    const std::size_t iq2 = station_of(m2.get_struct(), "Q");
    CHECK(s2r.QN(iq2, 0) == doctest::Approx(0.020833).epsilon(1e-4));
    CHECK(s2r.RN(iq2, 0) == doctest::Approx(0.041667).epsilon(1e-4));
    CHECK(s2r.UN(iq2, 0) == doctest::Approx(0.020833).epsilon(1e-4));
    CHECK(s2r.TN(iq2, 0) == doctest::Approx(0.5).epsilon(1e-9));

    // A Markov-modulated ARRIVAL as well: MMPP2(0.6, 0.4, 0.5, 0.3).
    qn::Network<double> m3("mfq3");
    const std::size_t s3 = m3.add_source("Source");
    const std::size_t q3 = m3.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t k3 = m3.add_sink("Sink");
    const std::size_t c3 = m3.add_open_class("C1");
    {
        Matrix<double> D0(2, 2), D1(2, 2);
        const double l0 = 0.6, l1 = 0.4, sg0 = 0.5, sg1 = 0.3;
        D1(0, 0) = l0;
        D1(1, 1) = l1;
        D0(0, 0) = -(l0 + sg0);
        D0(0, 1) = sg0;
        D0(1, 0) = sg1;
        D0(1, 1) = -(l1 + sg1);
        m3.set_arrival(s3, c3, D::map_dist(D0, D1, lang::ProcessType::MMPP2));
    }
    m3.set_service(q3, c3, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P3;
    P3.set(c3, c3, s3, q3, 1.0);
    P3.set(c3, c3, q3, k3, 1.0);
    m3.link(P3);
    const fluid::FluidSolution s3r = fluid::solver_fluid(m3.get_struct(), o);
    const std::size_t iq3 = station_of(m3.get_struct(), "Q");
    CHECK(s3r.QN(iq3, 0) == doctest::Approx(0.019508).epsilon(1e-4));
    CHECK(s3r.RN(iq3, 0) == doctest::Approx(0.041070).epsilon(1e-4));
    CHECK(s3r.TN(iq3, 0) == doctest::Approx(0.475).epsilon(1e-4));
}

TEST_CASE("fluid: a DPS model resolves to closing, and matrix refuses it") {
    // The reference resolves the method from the MODEL: a DPS station forces
    // `closing` at the ANALYZER level, because the matrix method's theta has no
    // per-class weight and cannot express DPS. Asking for matrix explicitly is an
    // error, not a silent downgrade. Expected values are MATLAB
    // SolverFLD(model,'method','closing') at LINE 3.0.7.
    //
    // THESE NUMBERS CHANGED WHEN THE DPS BRANCH WAS CORRECTED, and the old ones
    // were a stale MATLAB. The reference's DPS share used to carry an ADDITIVE
    // mean(w) term in its denominator and to scale by the full server count rather
    // than by psi(xi) = min(xi,c)*alpha(xi); this port had copied that. With both
    // fixed, equal weights reduce DPS to PS identically.
    //
    // Re-recorded 2026-09-08 against MATLAB R2026a: 34ca22581 replaced the moment
    // closure with the joint share-capacity one, moving both rows. Every codebase
    // meets them -- the JAR and native python to 2e-7, C++ to the digits it
    // prints -- so these are a four-way consensus and not a MATLAB-only row.
    qn::Network<double> m("dps");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", lang::SchedStrategy::DPS);
    const std::size_t c1 = m.add_closed_class("Class1", 2, d);
    const std::size_t c2 = m.add_closed_class("Class2", 1, d);
    m.set_service(d, c1, D::exp_rate(3.0));
    m.set_service(d, c2, D::exp_rate(0.5));
    m.set_service(q1, c1, D::exp_rate(0.1));
    m.set_service(q1, c2, D::exp_rate(1.0));
    m.set_service(q2, c1, D::exp_rate(0.1));
    m.set_service(q2, c2, D::exp_rate(1.0));
    m.set_sched_param(q2, c1, 1.0);
    m.set_sched_param(q2, c2, 5.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q1, 0.3);
    P.set(c1, c1, d, q2, 0.7);
    P.set(c1, c1, q1, d, 1.0);
    P.set(c1, c1, q2, d, 1.0);
    P.set(c2, c2, d, q1, 0.7);
    P.set(c2, c2, d, q2, 0.3);
    P.set(c2, c2, q1, d, 1.0);
    P.set(c2, c2, q2, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const fluid::FluidSolution s = fluid::solver_fluid(sn, fluid::FluidOptions());
    CHECK(s.method == "closing");  // resolved from the model, not the name
    const std::size_t id = station_of(sn, "Delay"), i1 = station_of(sn, "Queue1"),
                      i2 = station_of(sn, "Queue2");
    CHECK(s.QN(id, 0) == doctest::Approx(0.042857171).epsilon(1e-3));
    CHECK(s.QN(i1, 0) == doctest::Approx(0.385715121).epsilon(1e-3));
    CHECK(s.QN(i2, 0) == doctest::Approx(1.571427708).epsilon(1e-3));
    CHECK(s.UN(i2, 0) == doctest::Approx(0.758620589).epsilon(1e-3));
    CHECK(s.TN(i2, 0) == doctest::Approx(0.0758620589).epsilon(1e-3));

    // Explicitly asking for the method that cannot express DPS is refused.
    fluid::FluidOptions o;
    o.method = "matrix";
    CHECK_THROWS_AS(fluid::solver_fluid(sn, o), UnsupportedError);
    o.method = "pnorm";
    CHECK_THROWS_AS(fluid::solver_fluid(sn, o), UnsupportedError);

    // AT THE RUNNER LEVEL the same model resolves to `minnormal` instead, which is
    // `runAnalyzer`'s resolution and not the analyzer's: the second-order closure
    // is preferred wherever it applies, and `closing` for DPS is only the
    // fallback. `solver_fluid` stays the port of solver_fluid_analyzer.m alone,
    // which is why the two answers differ here by design.
    // MATLAB SolverFLD(model) -> 'default/minnormal'.
    const fluid::FluidSolution sr = fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions());
    CHECK(sr.method == "minnormal");
    CHECK(sr.QN(id, 0) == doctest::Approx(0.038569432).epsilon(1e-3));
    CHECK(sr.QN(i1, 0) == doctest::Approx(0.515215613).epsilon(1e-3));
    CHECK(sr.QN(i2, 0) == doctest::Approx(1.446214955).epsilon(1e-3));
    CHECK(sr.UN(i2, 0) == doctest::Approx(0.809956966).epsilon(1e-3));
    CHECK(sr.TN(i2, 0) == doctest::Approx(0.0809956966).epsilon(1e-3));
}

TEST_CASE("fluid transient matches the MATLAB getTranAvg trajectory") {
    // Delay(1) <-> PS Queue(2), N = 4, started with every job at the reference
    // station. The trajectory, not just its limit, is the thing under test.
    qn::Network<double> m("tr");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const std::vector<fluid::FluidTranPoint> tr =
        fluid::solver_fluid_transient(m.get_struct(), fluid::FluidOptions(), 5.0, 501);
    REQUIRE(tr.size() == 501);

    // The initial state is the reference station holding the whole population;
    // spreading it would reach the same fixed point by the wrong path.
    CHECK(tr.front().t == doctest::Approx(0.0));
    CHECK(tr.front().QN(0, 0) == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(tr.front().QN(1, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // MATLAB getTranAvg on the same model, read at t = 5.
    CHECK(tr.back().t == doctest::Approx(5.0));
    CHECK(tr.back().QN(0, 0) == doctest::Approx(2.010717).epsilon(1e-4));
    CHECK(tr.back().QN(1, 0) == doctest::Approx(1.989283).epsilon(1e-4));

    // Population is conserved at every point, and the queue fills monotonically
    // here, which no reference value can fake.
    for (const fluid::FluidTranPoint& pt : tr)
        CHECK(pt.QN(0, 0) + pt.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-5));
    for (std::size_t j = 1; j < tr.size(); ++j)
        CHECK(tr[j].QN(1, 0) >= tr[j - 1].QN(1, 0) - 1e-9);

    // The long-run limit of the transient IS the steady state.
    const fluid::FluidSolution ss = fluid::solver_fluid(m.get_struct(), fluid::FluidOptions());
    const std::vector<fluid::FluidTranPoint> lg =
        fluid::solver_fluid_transient(m.get_struct(), fluid::FluidOptions(), 200.0, 51);
    CHECK(lg.back().QN(0, 0) == doctest::Approx(ss.QN(0, 0)).epsilon(1e-4));
    CHECK(lg.back().QN(1, 0) == doctest::Approx(ss.QN(1, 0)).epsilon(1e-4));

    CHECK_THROWS_AS(fluid::solver_fluid_transient(m.get_struct(), fluid::FluidOptions(), -1.0),
                    InputError);
}

TEST_CASE("fluid Jacobian is the linearization of the drift") {
    // The reference builds this symbolically and ships it to a SAGE service;
    // this is the numerical counterpart. The check that does not depend on the
    // method: a closed model conserves mass, so every COLUMN of the Jacobian
    // must sum to zero -- perturbing one state cannot create or destroy any.
    qn::Network<double> m("jac");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), fluid::FluidOptions());
    const Matrix<double> J = fluid::fluid_jacobian(m.get_struct(), s.xvec);
    REQUIRE(J.rows() == s.xvec.size());
    for (std::size_t j = 0; j < J.cols(); ++j) {
        double col = 0.0;
        for (std::size_t i = 0; i < J.rows(); ++i) col += J(i, j);
        CHECK(col == doctest::Approx(0.0).epsilon(1e-6).scale(1.0));
    }
    CHECK_THROWS_AS(fluid::fluid_jacobian(m.get_struct(), std::vector<double>{1.0}), InputError);
}

TEST_CASE("fluid getProbAggr matches the MATLAB value") {
    // Delay <-> PS Queue, N = 4. A FIRST-ORDER method carries no second moment,
    // so the probability of the default marginal is FITTED to the mean queue
    // lengths: a Schmidt binomial per closed class. MATLAB reports 0.0625 at
    // both stations under `closing`, which is also what the symmetric split 2/2
    // out of 4 gives: C(4,4) p^4 with p = 1/2 at the reference station. The
    // method is pinned because `default` resolves to `minnormal`, which answers
    // the same question from its covariance instead (next case).
    qn::Network<double> m("pa");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "closing";
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    double lp = 0.0;
    const double p1 = fluid::fluid_prob_aggr(m.get_struct(), s, 1, &lp);
    const double p2 = fluid::fluid_prob_aggr(m.get_struct(), s, 2);
    CHECK(p1 == doctest::Approx(0.0625).epsilon(1e-6));
    CHECK(p2 == doctest::Approx(0.0625).epsilon(1e-6));
    CHECK(std::exp(lp) == doctest::Approx(p1).epsilon(1e-12));
    CHECK_THROWS_AS(fluid::fluid_prob_aggr(m.get_struct(), s, 99), InputError);
}

TEST_CASE("fluid getProbAggr reads the cell of the moment closure") {
    // The same model under `minnormal`, which DOES carry a covariance: the
    // answer is the probability its multivariate normal puts on the unit cell
    // around the state, evaluated by the Genz transformation with a
    // deterministic lattice. MATLAB returns 0.0750552761 at both stations
    // (symmetric here: the state is (4,0), so one cell is [3.5, Inf) and the
    // other (-Inf, 0.5] on the same one-dimensional fluctuation).
    //
    // RE-RECORDED 2026-09-01 against MATLAB R2026a. The previous 0.0750694807
    // was measured before the closure alternation and its inner mean solve were
    // tightened from CoarseTol to mom_tol = 1e-6, which moves the converged
    // answer by 1.9e-4 -- far past the 1e-6 asserted here. These are MATLAB's
    // numbers: this port reproduces them to 1.8e-7 and native python to 1.8e-7,
    // which is what licenses the re-record. Cf. the same re-pin of
    // test_fld_probaggr.py and of MinNormalTest.testOpenMm1 in the JAR.
    qn::Network<double> m("pa");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "minnormal";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(m.get_struct(), o);
    REQUIRE(s.has_moments);
    double lp = 0.0;
    const double p1 = fluid::fluid_prob_aggr(m.get_struct(), s, 1, &lp);
    const double p2 = fluid::fluid_prob_aggr(m.get_struct(), s, 2);
    CHECK(p1 == doctest::Approx(0.0750552761).epsilon(1e-6));
    CHECK(p2 == doctest::Approx(0.0750552761).epsilon(1e-6));
    CHECK(std::exp(lp) == doctest::Approx(p1).epsilon(1e-12));
}

TEST_CASE("fluid passage time gives the response-time CDF") {
    // Delay(1) <-> PS Queue(2), N = 4. The queue is saturated, so it drains at
    // a CONSTANT rate and the marked fluid leaves exponentially: the CDF must
    // be 1 - exp(-t/R) with R the mean response time the solver reports. That
    // is an analytic check the reference values cannot influence.
    qn::Network<double> m("pt");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "closing";
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    const fluid::FluidPassage pt = fluid::fluid_passage_time(m.get_struct(), s.xvec, 2, 1);

    REQUIRE(pt.t.size() > 10);
    CHECK(pt.fluid0 == doctest::Approx(2.0).epsilon(1e-4));  // the queue's mass
    CHECK(pt.cdf.front() == doctest::Approx(0.0).epsilon(1e-9));

    // A CDF: starts at 0, never decreases, ends at 1.
    for (std::size_t j = 1; j < pt.cdf.size(); ++j) CHECK(pt.cdf[j] >= pt.cdf[j - 1] - 1e-9);
    CHECK(pt.cdf.back() > 0.99);

    const double R = s.RN(1, 0);
    CHECK(R == doctest::Approx(1.0).epsilon(1e-4));
    for (double x : {0.5, 1.0, 2.0, 5.0}) {
        std::size_t b = 0;
        double bd = 1e9;
        for (std::size_t j = 0; j < pt.t.size(); ++j) {
            const double dd = std::fabs(pt.t[j] - x);
            if (dd < bd) { bd = dd; b = j; }
        }
        CHECK(pt.cdf[b] == doctest::Approx(1.0 - std::exp(-x / R)).epsilon(2e-3));
    }
    // MATLAB getCdfRespT on the same model, for the record: 0.393056, 0.631468,
    // 0.864318, 0.993277 -- this port agrees to ~6e-4 and sits closer to the
    // analytic curve, MATLAB reading its coarser adaptive grid through interp1.

    CHECK_THROWS_AS(fluid::fluid_passage_time(m.get_struct(), s.xvec, 99, 1), InputError);
    CHECK_THROWS_AS(fluid::fluid_passage_time(m.get_struct(), s.xvec, 2, 99), InputError);
}

TEST_CASE("getCdfRespT solves first and pairs the state with the struct it was solved on") {
    // `solver_fluid_cdf_respt` is the reference's getCdfRespT: it runs the solve
    // itself, because the state vector to mark a job in only exists afterwards
    // and is laid out by the phase counts of the struct the FCFS refit produced,
    // not of the model's. Two models, one of each kind.
    qn::Network<double> m("cdfrt");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "closing";
    const std::vector<std::vector<fluid::FluidPassage> > RD =
        fluid::solver_fluid_cdf_respt(m.get_struct(), o);
    REQUIRE(RD.size() == 2);
    REQUIRE(RD[0].size() == 1);
    // Every pair the model serves carries a law, and it is the SAME law the
    // per-pair entry point returns from the same solve.
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    const fluid::FluidPassage one = fluid::fluid_passage_time(m.get_struct(), s.xvec, 2, 1);
    REQUIRE(RD[1][0].cdf.size() == one.cdf.size());
    for (std::size_t j = 0; j < one.cdf.size(); ++j)
        CHECK(RD[1][0].cdf[j] == doctest::Approx(one.cdf[j]).epsilon(1e-9));
    for (std::size_t j = 1; j < RD[0][0].cdf.size(); ++j)
        CHECK(RD[0][0].cdf[j] >= RD[0][0].cdf[j - 1] - 1e-9);

    // An FCFS station with SCV != 1 is REFITTED, so the solved state is longer
    // than the model's own layout: pairing it with the input struct would be a
    // length error, and this entry point is what prevents a caller from doing so.
    qn::Network<double> mf("cdfrt-fcfs");
    const std::size_t d2 = mf.add_delay("Delay");
    const std::size_t q2 = mf.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t c2 = mf.add_closed_class("C1", 3, d2);
    mf.set_service(d2, c2, D::exp_rate(1.0));
    mf.set_service(q2, c2, D::hyperexp(0.3, 4.0, 0.5));
    qn::RoutingMatrix<double> P2;
    P2.set(c2, c2, d2, q2, 1.0);
    P2.set(c2, c2, q2, d2, 1.0);
    mf.link(P2);
    const std::vector<std::vector<fluid::FluidPassage> > RF =
        fluid::solver_fluid_cdf_respt(mf.get_struct(), o);
    REQUIRE(RF.size() == 2);
    REQUIRE(!RF[1][0].cdf.empty());
    CHECK(RF[1][0].cdf.front() == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(RF[1][0].cdf.back() > 0.99);
    for (std::size_t j = 1; j < RF[1][0].cdf.size(); ++j)
        CHECK(RF[1][0].cdf[j] >= RF[1][0].cdf[j - 1] - 1e-9);
}

TEST_CASE("getTranAvg resolves its own horizon when the timespan is unbounded") {
    // `options.timespan = [0, Inf]` does not mean "integrate forever":
    // `@NetworkSolver/getTranAvg.m` resolves an unspecified end time to
    // 30/minrate, minrate being the slowest finite entry of sn.rates. Service
    // is Exp(1) at the Delay and Exp(2) at the Queue, so MATLAB warns
    // "setting the timespan option to [0,30]" on this model and the port must
    // land on the same 30. A caller who names a horizon gets exactly that.
    qn::Network<double> m("tranauto");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    const std::vector<fluid::FluidTranPoint> au =
        fluid::solver_fluid_tran_avg(m.get_struct(), fluid::FluidOptions());
    REQUIRE(au.size() == 101);
    CHECK(au.front().t == doctest::Approx(0.0));
    CHECK(au.front().QN(0, 0) == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(au.back().t == doctest::Approx(30.0));
    // The horizon is long enough to have converged, which is what makes it the
    // reference's own stopping rule and not an arbitrary number.
    fluid::FluidOptions oc;
    oc.method = "closing";
    const fluid::FluidSolution ss = fluid::solver_fluid(m.get_struct(), oc);
    CHECK(au.back().QN(0, 0) == doctest::Approx(ss.QN(0, 0)).epsilon(1e-4));
    CHECK(au.back().QN(1, 0) == doctest::Approx(ss.QN(1, 0)).epsilon(1e-4));

    fluid::FluidOptions ot;
    ot.timespan_end = 3.0;
    const std::vector<fluid::FluidTranPoint> gi =
        fluid::solver_fluid_tran_avg(m.get_struct(), ot, 61);
    REQUIRE(gi.size() == 61);
    CHECK(gi.back().t == doctest::Approx(3.0));
}

TEST_CASE("fluid mfq priority matches the MATLAB getAvg on two modulated classes") {
    // Source -> Queue -> Sink with TWO open classes at different priorities, so
    // `mfq` takes the fluid priority branch. Both arrivals are MMPP2 and the
    // service is class-independent, which is what the branch requires.
    //
    // A FLUID LEVEL IS NOT A JOB COUNT: the high-priority class peaks at rate
    // 0.6 against a drain rate of 0.8, so it never builds up and its level is
    // ZERO exactly. Only the low-priority class, drained on what is left, holds
    // fluid. Expected values are MATLAB SolverFLD(model, method 'mfq').getAvg.
    qn::Network<double> m("mfqprio");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1", 0);  // smaller prio = higher
    const std::size_t c2 = m.add_open_class("C2", 1);
    auto mmpp2 = [](double l0, double l1, double s0, double s1) {
        Matrix<double> D0(2, 2), D1(2, 2);
        D1(0, 0) = l0;
        D1(1, 1) = l1;
        D0(0, 0) = -(l0 + s0);
        D0(0, 1) = s0;
        D0(1, 0) = s1;
        D0(1, 1) = -(l1 + s1);
        return D::map_dist(D0, D1, lang::ProcessType::MMPP2);
    };
    m.set_arrival(src, c1, mmpp2(0.6, 0.4, 0.5, 0.3));
    m.set_arrival(src, c2, mmpp2(0.3, 0.2, 0.4, 0.6));
    m.set_service(q, c1, D::exp_rate(0.8));
    m.set_service(q, c2, D::exp_rate(0.8));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "mfq";
    const fluid::FluidSolution s = fluid::solver_fluid(m.get_struct(), o);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t iq = station_of(sn, "Q"), is = station_of(sn, "Source");

    CHECK(s.QN(iq, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.QN(iq, 1) == doctest::Approx(0.1024710990).epsilon(1e-4));
    CHECK(s.RN(iq, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.RN(iq, 1) == doctest::Approx(0.3941196117).epsilon(1e-4));
    CHECK(s.UN(iq, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(s.UN(iq, 1) == doctest::Approx(0.1024710990).epsilon(1e-4));
    CHECK(s.TN(iq, 0) == doctest::Approx(0.475).epsilon(1e-9));
    CHECK(s.TN(iq, 1) == doctest::Approx(0.26).epsilon(1e-9));
    CHECK(s.TN(is, 0) == doctest::Approx(0.475).epsilon(1e-9));
    CHECK(s.TN(is, 1) == doctest::Approx(0.26).epsilon(1e-9));
    CHECK_FALSE(s.has_aoi);

    // Unmodulated arrivals leave the priority model degenerate, and the
    // reference falls back to the matrix method rather than reporting zeros.
    qn::Network<double> m2("mfqprio2");
    const std::size_t s2 = m2.add_source("Source");
    const std::size_t q2 = m2.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t k2 = m2.add_sink("Sink");
    const std::size_t d1 = m2.add_open_class("C1", 0);
    const std::size_t d2 = m2.add_open_class("C2", 1);
    m2.set_arrival(s2, d1, D::exp_rate(0.2));
    m2.set_arrival(s2, d2, D::exp_rate(0.1));
    m2.set_service(q2, d1, D::exp_rate(0.8));
    m2.set_service(q2, d2, D::exp_rate(0.8));
    qn::RoutingMatrix<double> P2;
    P2.set(d1, d1, s2, q2, 1.0);
    P2.set(d1, d1, q2, k2, 1.0);
    P2.set(d2, d2, s2, q2, 1.0);
    P2.set(d2, d2, q2, k2, 1.0);
    m2.link(P2);
    const fluid::FluidSolution s2r = fluid::solver_fluid(m2.get_struct(), o);
    const std::size_t iq2 = station_of(m2.get_struct(), "Q");
    CHECK(s2r.TN(iq2, 0) == doctest::Approx(0.2).epsilon(1e-4));
    CHECK(s2r.TN(iq2, 1) == doctest::Approx(0.1).epsilon(1e-4));
    CHECK(s2r.QN(iq2, 0) > 0.0);  // the matrix method, not a degenerate zero
}

TEST_CASE("fluid mfq AoI matches the MATLAB getAvgAoI on bufferless and buffered queues") {
    // Source(Exp 0.5) -> Queue(capacity 1 or 2, one server) -> Sink. A capacity
    // of one or two is what makes `mfq` take the AGE branch instead of the
    // ordinary fluid queue. Expected values are MATLAB
    // SolverFLD(model, method 'mfq').getAvgAoI().
    auto build = [](const std::string& nm, double cap, lang::SchedStrategy sched,
                    const D& svc) {
        qn::Network<double>* m = new qn::Network<double>(nm);
        const std::size_t src = m->add_source("Source");
        const std::size_t q = m->add_queue("Q", sched);
        const std::size_t snk = m->add_sink("Sink");
        const std::size_t c = m->add_open_class("C1");
        m->set_arrival(src, c, D::exp_rate(0.5));
        m->set_service(q, c, svc);
        m->set_number_of_servers(q, 1);
        m->set_capacity(q, cap);
        qn::RoutingMatrix<double> P;
        P.set(c, c, src, q, 1.0);
        P.set(c, c, q, snk, 1.0);
        m->link(P);
        return m;
    };
    fluid::FluidOptions o;
    o.method = "mfq";

    struct Case {
        const char* name;
        double cap;
        lang::SchedStrategy sched;
        bool erlang;
        double aoi_mean, aoi_var, paoi_mean, paoi_var, QN, UN, RN;
    };
    // Erlang(rate 2, 2 phases) has mean 1, but the reference's aoi_dist2ph
    // reads it as an age process of mean 0.5, and the M/M/1 metrics beside the
    // age laws follow that reading. Both codebases agree because both do it.
    const Case cases[] = {
        {"bufferless FCFS Exp", 1.0, lang::SchedStrategy::FCFS, false, 3.3333333333, 5.5555555556,
         4.0000000000, 6.0000000000, 1.0, 0.5, 2.0},
        {"bufferless LCFSPR Exp", 1.0, lang::SchedStrategy::LCFSPR, false, 3.0000000000,
         5.0000000000, 3.6666666667, 5.4444444444, 1.0, 0.5, 2.0},
        {"bufferless FCFS Erlang", 1.0, lang::SchedStrategy::FCFS, true, 2.6000000000, 4.3400000000,
         3.0000000000, 4.5000000000, 0.3333333333, 0.25, 0.6666666667},
        {"single buffer FCFS Exp", 2.0, lang::SchedStrategy::FCFS, false, 3.2857142857,
         5.2040816327, 3.6666666667, 5.2222222222, 1.0, 0.5, 2.0},
        {"single buffer LCFS Exp", 2.0, lang::SchedStrategy::LCFS, false, 3.1746031746,
         4.8954396573, 3.5555555556, 4.9135802469, 1.0, 0.5, 2.0},
        {"single buffer LCFS Erlang", 2.0, lang::SchedStrategy::LCFS, true, 2.5276190476,
         4.1915229025, 2.6800000000, 4.0776000000, 0.3333333333, 0.25, 0.6666666667},
    };
    for (const Case& cs : cases) {
        CAPTURE(cs.name);
        qn::Network<double>* m =
            build(cs.name, cs.cap, cs.sched, cs.erlang ? D::erlang(2.0, 2) : D::exp_rate(1.0));
        const fluid::FluidSolution s = fluid::solver_fluid(m->get_struct(), o);
        REQUIRE(s.has_aoi);
        CHECK(s.aoi.system_type == (cs.cap == 1.0 ? "bufferless" : "singlebuffer"));
        CHECK(s.aoi.aoi.mean == doctest::Approx(cs.aoi_mean).epsilon(1e-6));
        CHECK(s.aoi.aoi.var == doctest::Approx(cs.aoi_var).epsilon(1e-6));
        CHECK(s.aoi.paoi.mean == doctest::Approx(cs.paoi_mean).epsilon(1e-6));
        CHECK(s.aoi.paoi.var == doctest::Approx(cs.paoi_var).epsilon(1e-6));
        const std::size_t iq = station_of(m->get_struct(), "Q");
        CHECK(s.QN(iq, 0) == doctest::Approx(cs.QN).epsilon(1e-6));
        CHECK(s.UN(iq, 0) == doctest::Approx(cs.UN).epsilon(1e-6));
        CHECK(s.RN(iq, 0) == doctest::Approx(cs.RN).epsilon(1e-6));
        CHECK(s.TN(iq, 0) == doctest::Approx(0.5).epsilon(1e-6));

        CHECK(fluid::aoi_cdf(s.aoi.aoi, 0.0) == doctest::Approx(0.0).epsilon(1e-12));
        delete m;
    }

    // The preemption probability is overridable, as options.config.aoi_preemption
    // is in the reference: forcing p = 1 on an FCFS model must reproduce the
    // LCFS-PR answer.
    qn::Network<double>* mf = build("override", 1.0, lang::SchedStrategy::FCFS, D::exp_rate(1.0));
    fluid::FluidOptions ov = o;
    ov.aoi_preemption = 1.0;
    const fluid::FluidSolution so = fluid::solver_fluid(mf->get_struct(), ov);
    CHECK(so.aoi.preemption == doctest::Approx(1.0));
    CHECK(so.aoi.aoi.mean == doctest::Approx(3.0000000000).epsilon(1e-6));
    delete mf;

    // The age CDF, against MATLAB getCdfAoI on the same two models (re-derived
    // 2026-07-31; MATLAB, native Python and this port agree to ten digits).
    //
    // The law is F(t) = 1 + g exp(A t) A^-1 h: (g, A, h) is a DENSITY triple
    // (g is normalized by -g A^-1 h), so g exp(A t) h is the density and the
    // survival function carries the extra A^-1. All four codebases used to
    // report 1 - g exp(A t) h, which FELL before it rose. Integrating the
    // corrected survival function returns the mean the same solve reports
    // (3.333333 and 3.174603 here), which is what pins these numbers down.
    struct CdfCase {
        double cap;
        lang::SchedStrategy sched;
        double t[5], aoi[5], paoi[5];
    };
    const CdfCase cdfs[] = {
        {1.0,
         lang::SchedStrategy::FCFS,
         {0.5, 1.0, 2.0, 5.0, 10.0},
         {0.0351707879, 0.1183437898, 0.3347704844, 0.8035664937, 0.9822591410},
         {0.0076541767, 0.0453951258, 0.2051586515, 0.7255635815, 0.9736384111}},
        {2.0,
         lang::SchedStrategy::LCFS,
         {0.5, 1.0, 2.0, 5.0, 10.0},
         {0.0331582787, 0.1177941801, 0.3488869584, 0.8300834581, 0.9863345348},
         {0.0101434992, 0.0592510793, 0.2561431516, 0.7935848406, 0.9838344652}},
    };
    for (const CdfCase& cc : cdfs) {
        qn::Network<double>* mc = build("cdf", cc.cap, cc.sched, D::exp_rate(1.0));
        const fluid::FluidSolution sc = fluid::solver_fluid(mc->get_struct(), o);
        REQUIRE(sc.has_aoi);
        for (int j = 0; j < 5; ++j) {
            CHECK(fluid::aoi_cdf(sc.aoi.aoi, cc.t[j]) == doctest::Approx(cc.aoi[j]).epsilon(1e-6));
            CHECK(fluid::aoi_cdf(sc.aoi.paoi, cc.t[j]) ==
                  doctest::Approx(cc.paoi[j]).epsilon(1e-6));
        }
        delete mc;
    }

    // A capacity the age model cannot express is not an AoI model at all, and
    // falls through to the ordinary single fluid queue.
    qn::Network<double>* mb = build("cap3", 3.0, lang::SchedStrategy::FCFS, D::exp_rate(1.0));
    const fluid::FluidSolution sb = fluid::solver_fluid(mb->get_struct(), o);
    CHECK_FALSE(sb.has_aoi);
    delete mb;
}

TEST_CASE("fluid exportODEs writes the system MATLAB writes") {
    // Delay <-> PS Queue, N = 2, exponential service: small enough that every
    // line of the export can be compared by eye against MATLAB's, which is how
    // the expected strings below were obtained.
    qn::Network<double> m("fluid_cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "closing";
    const std::string tex =
        fluid::fluid_export_odes(m.get_struct(), o, "scalar", "fluid_cqn");
    auto has = [&tex](const std::string& s) { return tex.find(s) != std::string::npos; };

    // The machine-readable header is the contract a parser reads.
    CHECK(has("% form: dx/dt = J*r(x)"));
    CHECK(has("% method: closing"));
    CHECK(has("% nstates: 2"));
    CHECK(has("% nevents: 2"));
    CHECK(has("% STATE 1 station=Delay class=C1 phase=1"));
    CHECK(has("% STATE 2 station=Queue1 class=C1 phase=1"));
    CHECK(has("% EVENT 1 var=1 type=lin coeff=1"));
    CHECK(has("% EVENT 2 var=2 type=min coeff=2"));
    // The Delay is a linear term, the queue a min-scaled one, and the two
    // equations are each other's negative because the population is conserved.
    CHECK(has("\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -x_{1} + "
              "2\\,x_{2}\\,g_{2}(\\mathbf{x})"));
    CHECK(has("\\frac{\\mathrm{d}x_{2}}{\\mathrm{d}t} &= x_{1} - "
              "2\\,x_{2}\\,g_{2}(\\mathbf{x})"));
    CHECK(has("g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 1)}{n_{2}(\\mathbf{x})}"));
    CHECK(has("n_{2}(\\mathbf{x}) &= x_{2}"));
    CHECK(has("\\mathbf{x}(0) = \\begin{pmatrix} 2 & 0 \\end{pmatrix}^{\\top}"));
    CHECK(has("\\documentclass{article}"));
    CHECK(has("\\end{document}"));

    // The matrix notation prints J and the rate functions instead.
    const std::string texm =
        fluid::fluid_export_odes(m.get_struct(), o, "matrix", "fluid_cqn");
    CHECK(texm.find("J = \\begin{bmatrix} -1 & 1 \\\\ 1 & -1 \\end{bmatrix}") !=
          std::string::npos);
    CHECK(texm.find("r_{1}(\\mathbf{x}) &= x_{1}") != std::string::npos);
    CHECK(texm.find("r_{2}(\\mathbf{x}) &= 2\\,x_{2}\\,g_{2}(\\mathbf{x})") != std::string::npos);
    CHECK(texm.find("% notation: matrix") != std::string::npos);

    // The W form is the other solver path, and prints W' and theta.
    fluid::FluidOptions ow;
    ow.method = "default";
    const std::string texw =
        fluid::fluid_export_odes(m.get_struct(), ow, "matrix", "fluid_cqn");
    CHECK(texw.find("% form: dx/dt = W^T*theta(x) + lambda") != std::string::npos);
    CHECK(texw.find("% method: matrix") != std::string::npos);
    CHECK(texw.find("W^{\\top} = \\begin{bmatrix} -1 & 2 \\\\ 1 & -2 \\end{bmatrix}") !=
          std::string::npos);
    CHECK(texw.find("\\theta(\\mathbf{x}) = \\begin{bmatrix} x_{1}\\,g_{1}(\\mathbf{x}) \\\\ "
                    "x_{2}\\,g_{2}(\\mathbf{x}) \\end{bmatrix}") != std::string::npos);

    // An open model: the Source contributes the `ext1` factor under closing,
    // and a two-phase service splits the queue block into two states.
    qn::Network<double> m2("fluid_oqn");
    const std::size_t src = m2.add_source("Source");
    const std::size_t q2 = m2.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t snk = m2.add_sink("Sink");
    const std::size_t c2 = m2.add_open_class("C1");
    m2.set_arrival(src, c2, D::exp_rate(0.5));
    m2.set_service(q2, c2, D::erlang(2.0, 2));
    qn::RoutingMatrix<double> P2;
    P2.set(c2, c2, src, q2, 1.0);
    P2.set(c2, c2, q2, snk, 1.0);
    m2.link(P2);
    const std::string texo =
        fluid::fluid_export_odes(m2.get_struct(), o, "scalar", "fluid_oqn");
    CHECK(texo.find("% EVENT 1 var=1 type=ext1 coeff=0.5") != std::string::npos);
    CHECK(texo.find("% STATE 2 station=Q class=C1 phase=1") != std::string::npos);
    CHECK(texo.find("% STATE 3 station=Q class=C1 phase=2") != std::string::npos);
    // A single-phase source class carries unit mass, so its `ext1` factor is
    // empty and the arrival appears as the bare rate.
    CHECK(texo.find("\\frac{\\mathrm{d}x_{2}}{\\mathrm{d}t} &= 0.5 - "
                    "2\\,x_{2}\\,g_{2}(\\mathbf{x})") != std::string::npos);

    // Methods without an ODE system of either shape are refused by name.
    fluid::FluidOptions bad;
    bad.method = "mfq";
    CHECK_THROWS_AS(fluid::fluid_export_odes(m.get_struct(), bad, "scalar", "x"),
                    UnsupportedError);
    bad.method = "tbi";
    CHECK_THROWS_AS(fluid::fluid_export_odes(m.get_struct(), bad, "scalar", "x"),
                    UnsupportedError);
    CHECK_THROWS_AS(fluid::fluid_export_odes(m.get_struct(), o, "graph", "x"), InputError);
    // statedep has no open-model branch, as in the reference.
    fluid::FluidOptions sd;
    sd.method = "statedep";
    CHECK_THROWS_AS(fluid::fluid_export_odes(m2.get_struct(), sd, "scalar", "x"),
                    UnsupportedError);
}

TEST_CASE("fluid exportODEs typesets the DPS drift the solver integrates") {
    // The exported DPS factor used to be x_v/(mean(w) + ntilde_i) with S_i*w_ir
    // in the coefficient, which is what the closing rates computed BEFORE the
    // additive seed was dropped and the full server count replaced by
    // min(n_i,S_i). The export kept the old form, so it typeset a drift the
    // solver had stopped integrating -- and because every codebase copied the
    // same stale factor, no cross-codebase comparison could see it.
    //
    // TWO SERVERS ON PURPOSE. At S_i = 1 the two coefficient conventions,
    // S_i*w_ir and w_ir, are numerically identical, so a single-server model
    // pins the numerator and misses the regression entirely.
    qn::Network<double> m("dps_export");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("QueueDPS", lang::SchedStrategy::DPS);
    m.set_number_of_servers(q, 2);
    const std::size_t c1 = m.add_closed_class("C1", 2, d);
    const std::size_t c2 = m.add_closed_class("C2", 1, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(1.0));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(1.0));
    m.set_sched_param(q, c1, 2.0);
    m.set_sched_param(q, c2, 1.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);

    fluid::FluidOptions o;
    o.method = "closing";
    const std::string tex = fluid::fluid_export_odes(m.get_struct(), o, "scalar", "dps_export");
    auto has = [&tex](const std::string& s) { return tex.find(s) != std::string::npos; };

    // The share divides the CAPACITY min(n_i, S_i), with no additive seed.
    CHECK(has("g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 2)}{\\tilde{n}_{2}(\\mathbf{x})}"));
    CHECK_FALSE(has("\\frac{1}{0.5 + \\tilde{n}_{2}(\\mathbf{x})}"));
    // Weights normalized to (2/3, 1/3); the denominator needs BOTH n_2 and
    // ntilde_2 defined, where the old smooth form needed only ntilde_2.
    CHECK(has("\\tilde{n}_{2}(\\mathbf{x}) &= 0.66666667\\,(x_{3}) + 0.33333333\\,(x_{4})"));
    CHECK(has("n_{2}(\\mathbf{x}) &= "));
    // The coefficient carries w_ir alone: 2 * 2/3, not 2 * 2 * 2/3.
    CHECK(has("1.3333333\\,x_{3}\\,g_{2}(\\mathbf{x})"));
    CHECK_FALSE(has("2.6666667\\,x_{3}\\,g_{2}(\\mathbf{x})"));

    // min() has no derivative at n_i = S_i, so the symbolic drift must refuse
    // DPS under the closing family rather than export the old smooth ratio.
    const fluid::FluidSymSystem sys =
        fluid::fluid_symodes(m.get_struct(), "closing", 0.0, std::vector<double>());
    CHECK_THROWS_AS(fluid::fluid_symbolic_drift(sys), UnsupportedError);
    std::string msg;
    try {
        fluid::fluid_symbolic_drift(sys);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("dpsmin") != std::string::npos);
}

TEST_CASE("fluid: earlystop=false runs to iter_max, and the default stops on the tail") {
    // `solver_fluid_iteration.m` USED to have no working break, so it always ran
    // iter_max passes and this case asserted that. It now stops on the geometric
    // tail r*rho/(1-rho) of the window iteration, which bounds the distance still
    // left to the fixed point where the mass moved by ONE window understates it;
    // `options.config.fluid_earlystop = false` is what restores the unconditional
    // sweep. The invariant that survives the change is the one below: the two
    // answers agree to the tolerance the stop is sized against.
    qn::Network<double> m = build_fl2();
    const qn::NetworkStruct<double>& sn = m.get_struct();

    fluid::FluidOptions full;
    full.method = "closing";
    full.iter_max = 12;  // the reference's 200 is the same behaviour, slower
    full.earlystop = false;
    const fluid::FluidSolution sl = fluid::solver_fluid(sn, full);

    fluid::FluidOptions fast = full;
    fast.earlystop = true;
    fast.iter_tol = 1e-4;
    const fluid::FluidSolution sf = fluid::solver_fluid(sn, fast);

    CHECK(sl.iters == 12);       // every pass, as the reference does with the stop off
    CHECK(sf.iters < sl.iters);  // the default early stop
    // this model settles inside the window; the models where it does not are
    // why the stop is off by default (M/M/1 at rho = 0.9 lands 0.1% short)
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            CHECK(sl.QN(i, r) == doctest::Approx(sf.QN(i, r)).epsilon(1e-4));
}

}  // namespace

TEST_CASE("fluid: the FCFS non-exponential refit reproduces the MATLAB getAvgTable") {
    // Delay(Exp, Z = 1) -> FCFS Q1(HyperExp, mean 0.5, SCV 8) -> FCFS Q2(Erlang,
    // mean 0.4, SCV 1/3), N = 6. The fluid drift of an FCFS station is the PS
    // drift, so without `solver_fluid_analyzer.m`'s refit loop both queues would
    // be integrated at their declared phase structure and neither number below
    // is reachable. The oracle is MATLAB 3.0.7 `SolverFluid(model).getAvgTable`,
    // which resolves to `minnormal` on this model, as the runner does here.
    qn::Network<double> m("fluid_fcfs_nonexp");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 6, d);
    m.set_service(d, c, D::exp_rate(1.0));
    // HyperExp.fitMeanAndSCV(0.5, 8) in the balanced-means convention.
    const double scv = 8.0, mean = 0.5;
    const double p = 0.5 * (1.0 + std::sqrt((scv - 1.0) / (scv + 1.0)));
    m.set_service(q1, c, D::hyperexp(p, 2.0 * p / mean, 2.0 * (1.0 - p) / mean));
    m.set_service(q2, c, D::erlang(3.0 / 0.4, 3));  // mean 0.4, SCV 1/3
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions());
    const std::size_t id = station_of(sn, "Think"), i1 = station_of(sn, "Q1"),
                      i2 = station_of(sn, "Q2");
    CHECK(s.QN(id, 0) == doctest::Approx(1.7779).epsilon(1e-4));
    CHECK(s.QN(i1, 0) == doctest::Approx(2.6760).epsilon(1e-4));
    CHECK(s.QN(i2, 0) == doctest::Approx(1.5461).epsilon(1e-4));
    CHECK(s.UN(i1, 0) == doctest::Approx(0.88896).epsilon(1e-4));
    CHECK(s.UN(i2, 0) == doctest::Approx(0.71117).epsilon(1e-4));
    CHECK(s.RN(i1, 0) == doctest::Approx(1.5051).epsilon(1e-4));
    CHECK(s.RN(i2, 0) == doctest::Approx(0.86962).epsilon(1e-4));
    CHECK(s.TN(id, 0) == doctest::Approx(1.7779).epsilon(1e-4));
}

TEST_CASE("fluid kp: a cyclic MAPt is averaged over its period, not read at the horizon") {
    // The (MAP_t/Ph_t/inf)^N model of `matlab/examples/basic/openQN/oqn_mapt.m`:
    // a two-segment 2-phase MAPt repeating with period 2.5, feeding a Delay of
    // rate 2. A cyclic schedule has NO fixed point, so the steady table is the
    // time average over the last full period; reading the horizon instead gives
    // an arbitrary point of the cycle at which the Source and the Delay
    // throughputs do not even agree, which the last check here would catch.
    qn::Network<double> m("fluid_kp_mapt");
    const std::size_t src = m.add_source("Source");
    const std::size_t del = m.add_delay("Delay");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("OpenClass", src);
    std::vector<double> bp;
    bp.push_back(0.0);
    bp.push_back(1.0);
    bp.push_back(2.5);
    std::vector<Matrix<double> > D0s, D1s;
    Matrix<double> a0(2, 2), a1(2, 2), b0(2, 2), b1(2, 2);
    a0(0, 0) = -5.0; a0(0, 1) = 1.0;  a0(1, 0) = 2.0; a0(1, 1) = -4.0;
    a1(0, 0) = -12.0; a1(0, 1) = 3.0; a1(1, 0) = 5.0; a1(1, 1) = -9.0;
    b0(0, 0) = 3.0; b0(0, 1) = 1.0;   b0(1, 0) = 1.0; b0(1, 1) = 1.0;
    b1(0, 0) = 7.0; b1(0, 1) = 2.0;   b1(1, 0) = 2.0; b1(1, 1) = 2.0;
    D0s.push_back(a0); D0s.push_back(a1);
    D1s.push_back(b0); D1s.push_back(b1);
    m.set_service(src, c, D::mapt(bp, D0s, D1s, true));
    m.set_service(del, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, del, 1.0);
    P.set(c, c, del, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    fluid::FluidOptions o;
    o.method = "kp";
    o.tol = 1e-9;
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    const std::size_t is = station_of(sn, "Source"), idl = station_of(sn, "Delay");
    CHECK(s.TN(is, 0) == doctest::Approx(5.4301).epsilon(1e-4));
    CHECK(s.QN(idl, 0) == doctest::Approx(2.7151).epsilon(1e-4));
    CHECK(s.RN(idl, 0) == doctest::Approx(0.5).epsilon(1e-4));
    // Over a full period the flow in equals the flow out; at any single instant
    // of the cycle it does not.
    CHECK(s.TN(idl, 0) == doctest::Approx(s.TN(is, 0)).epsilon(1e-4));
}

TEST_CASE("fluid kp: the mean and covariance seeds are honoured, and a misfit is refused") {
    // M/M/inf started from a POISSON(x0) queue stays Poisson at every t, so a run
    // seeded with mean x0 and covariance x0 must return the exact transient mean
    // AND Var == Mean along the whole trajectory. Honouring only one of the two
    // seeds breaks the identity, so this pins both at once.
    const double lambda = 3.0, mu = 2.0, x0 = 5.0;
    qn::Network<double> m("fluid_kp_seed");
    const std::size_t src = m.add_source("Source");
    const std::size_t del = m.add_delay("Delay");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("OpenClass", src);
    std::vector<double> bp;
    bp.push_back(0.0);
    bp.push_back(1.0);
    std::vector<Matrix<double> > D0s, D1s;
    Matrix<double> a0(1, 1), b0(1, 1);
    a0(0, 0) = -lambda;
    b0(0, 0) = lambda;
    D0s.push_back(a0);
    D1s.push_back(b0);
    m.set_service(src, c, D::mapt(bp, D0s, D1s, true));
    m.set_service(del, c, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, del, 1.0);
    P.set(c, c, del, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t idl = station_of(sn, "Delay");
    fluid::FluidOptions o;
    o.method = "kp";
    o.tol = 1e-9;
    o.timespan_end = 4.0;
    // dim = 2: the arrival phase (u-block) then the Delay service phase.
    o.kp_init_sol.push_back(1.0);
    o.kp_init_sol.push_back(x0);
    o.init_cov = Matrix<double>(2, 2, 0.0);
    o.init_cov(1, 1) = x0;
    const fluid::FluidKpTransient tr = fluid::solver_fluid_tran_avg_var(sn, o);
    REQUIRE(tr.t.size() > 2);
    REQUIRE(tr.q.front().size() == 2);
    // The seed reached the integrator in the layout it was documented in.
    CHECK(tr.q.front()[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(tr.q.front()[1] == doctest::Approx(x0).epsilon(1e-12));
    for (std::size_t k = 0; k < tr.t.size(); ++k) {
        const double t = tr.t[k];
        const double exact =
            x0 * std::exp(-mu * t) + (lambda / mu) * (1.0 - std::exp(-mu * t));
        CHECK(tr.q[k][1] == doctest::Approx(exact).epsilon(1e-6));
        CHECK(tr.QVar[k](idl, 0) == doctest::Approx(tr.q[k][1]).epsilon(1e-6));
    }

    // A seed that does not fit is REFUSED, not dropped: dropping it would
    // integrate from the default initial condition under the caller's name and
    // hand back a plausible wrong trajectory.
    fluid::FluidOptions bad = o;
    bad.init_cov = Matrix<double>();
    bad.kp_init_sol.push_back(0.0);  // three entries against dim = 2
    CHECK_THROWS_AS(fluid::solver_fluid_kp(sn, bad), InputError);
    bad = o;
    bad.init_cov = Matrix<double>(1, 2, 0.0);
    CHECK_THROWS_AS(fluid::solver_fluid_kp(sn, bad), InputError);
    bad = o;
    bad.init_cov(0, 1) = 2.0;  // no longer symmetric
    CHECK_THROWS_AS(fluid::solver_fluid_kp(sn, bad), InputError);
}

TEST_CASE("fluid: a schedule is solved at the nominal, and exactly under a rate schedule") {
    // Three claims, and they are the reference's split between the model layer
    // and the options:
    //
    //   1. every method ACCEPTS a schedule. `solver_fluid.m:20-35` installs the
    //      width-weighted nominal pair and integrates it deliberately, so a
    //      steady-state request is answered at the time-averaged rate.
    //   2. asking for the schedule ITSELF is what `nhpp_sched` does, and it
    //      changes the answer -- otherwise the multiplier never reached the drift.
    //   3. the moment closures refuse the time-varying drift by name, because a
    //      non-autonomous system has no stationary covariance to converge to.
    qn::Network<double> m("fluid_kp_reject");
    const std::size_t src = m.add_source("Source");
    const std::size_t del = m.add_delay("Delay");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("OpenClass", src);
    std::vector<double> bp;
    bp.push_back(0.0);
    bp.push_back(1.0);
    bp.push_back(3.0);
    std::vector<double> rates;
    rates.push_back(1.0);
    rates.push_back(4.0);
    m.set_service(src, c, D::nhpp(bp, rates, true));
    m.set_service(del, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, del, 1.0);
    P.set(c, c, del, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t idl = station_of(sn, "Delay");

    // 1. The first-order methods answer, at the nominal. NHPP([0,1,3],[1,4])
    // averages to (1*1 + 4*2)/3 = 3, and the Delay of rate 2 then holds
    // lambda/mu = 1.5 jobs -- a number the nominal fixes exactly.
    fluid::FluidOptions o;
    o.method = "closing";
    o.tol = 1e-9;
    fluid::FluidSolution s1;
    CHECK_NOTHROW(s1 = fluid::solver_fluid_run_analyzer(sn, o));
    CHECK(s1.QN(idl, 0) == doctest::Approx(1.5).epsilon(1e-6));
    o.method = "matrix";
    CHECK_NOTHROW(fluid::solver_fluid_run_analyzer(sn, o));
    o.method = "kp";
    CHECK_NOTHROW(fluid::solver_fluid_run_analyzer(sn, o));

    // 2. `nhpp_sched` makes the drift follow the intensity. The transient
    // detects it itself, as `getTranAvg.m:76` does, so its trajectory must
    // OSCILLATE with the period-3 schedule rather than sit at the nominal.
    o.method = "closing";
    // `timespan_end` is the horizon the cyclic schedule is EXPANDED over, so it
    // is set to the same 6 the MATLAB run used; leaving it unbounded expands
    // three periods instead and puts the step grid at different points.
    o.timespan_end = 6.0;
    // A dense output grid: the pinned values below are read by LINEAR
    // interpolation, and the trajectory is steep where the schedule steps, so a
    // coarse grid would be testing the interpolation rather than the drift.
    const std::vector<fluid::FluidTranPoint> tr =
        fluid::solver_fluid_transient(sn, o, 6.0, 1201);
    REQUIRE(tr.size() > 2);
    double lo = tr.front().QN(idl, 0), hi = lo;
    for (std::size_t k = 0; k < tr.size(); ++k) {
        lo = std::min(lo, tr[k].QN(idl, 0));
        hi = std::max(hi, tr[k].QN(idl, 0));
    }
    CHECK(hi - lo > 0.5);  // the nominal drift would settle, not swing
    // Pinned to MATLAB `solver_fluid_closing` with options.config.nhpp_sched set,
    // read at points of ITS OWN integrator grid so no interpolation is involved.
    const double want[3] = {0.43234136, 1.97127833, 0.69912487};
    const double at[3] = {1.0, 3.0, 4.0};
    for (std::size_t a = 0; a < 3; ++a) {
        double got = tr.back().QN(idl, 0);
        for (std::size_t k = 0; k + 1 < tr.size(); ++k)
            if (tr[k].t <= at[a] && at[a] <= tr[k + 1].t) {
                const double w = (tr[k + 1].t > tr[k].t)
                                     ? (at[a] - tr[k].t) / (tr[k + 1].t - tr[k].t) : 0.0;
                got = (1.0 - w) * tr[k].QN(idl, 0) + w * tr[k + 1].QN(idl, 0);
                break;
            }
        CHECK(got == doctest::Approx(want[a]).epsilon(1e-6));
    }

    // 3. A time-varying drift has no stationary covariance, so the closures say so.
    fluid::FluidOptions mm;
    mm.method = "minnormal";
    mm.nhpp_sched = fluid::fluid_detect_nhpp(sn);
    REQUIRE(mm.nhpp_sched.size() == 1);
    CHECK_THROWS_AS(fluid::solver_fluid_run_analyzer(sn, mm), UnsupportedError);
}

TEST_CASE("fluid: the symbolic Jacobian is the derivative of the symbolic drift") {
    // `getJacobian` in the reference posts the drift to the line-sage-rest
    // backend; this port differentiates the STRUCTURED system instead. The oracle
    // is the definition: J(i,j) evaluated at a point must equal a central
    // difference of the printed drift there. Both sides are built from the same
    // FluidSymSystem, so this catches a wrong derivative rule, not a wrong drift.
    //
    // Only the smooth methods have a Jacobian at all, which is the same gate
    // `fluid_symbolic_drift` applies: `min(n,S)` has no derivative at n = S.
    // A closed Delay -> FCFS(Erlang) pair. Every event of this model carries a
    // factor the softmin smooths; `build_fl2`'s multiserver station does not, and
    // is used below as the refusal case.
    qn::Network<double> m("symjac");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 3, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidSymSystem sys =
        fluid::fluid_symodes(sn, "softmin", 20.0, std::vector<double>());
    const fluid::FluidSymbolicJacobian jj = fluid::fluid_symbolic_jacobian(sys);
    REQUIRE(jj.J.size() == sys.nstates);
    for (std::size_t i = 0; i < sys.nstates; ++i) REQUIRE(jj.J[i].size() == sys.nstates);
    // Every entry is an expression, never an empty string: a structurally zero
    // entry is the literal "0" so a consumer never has to guess.
    for (std::size_t i = 0; i < sys.nstates; ++i)
        for (std::size_t j = 0; j < sys.nstates; ++j) CHECK(!jj.J[i][j].empty());
    // A non-smooth method is refused by FACTOR TYPE rather than by name: the
    // `closing` drift scales by min(n_i, S_i), which has no derivative at n = S.
    const fluid::FluidSymSystem cl =
        fluid::fluid_symodes(sn, "closing", 20.0, std::vector<double>());
    CHECK_THROWS_AS(fluid::fluid_symbolic_jacobian(cl), UnsupportedError);
}

TEST_CASE("fluid: getJacobian answers locally and refuses the equilibria without a backend") {
    // `fluid_jacobian` is the whole of `@SolverFLD/getJacobian`: rhs and vars
    // from the drift, J from the backend when there is one, and the solutions of
    // f(x) = 0 only when asked for. This case pins the NO-BACKEND contract,
    // which is the half that must never reach the network: 'none' is honoured
    // rather than read as "not given", so the test cannot start a container.
    qn::Network<double> m("symjac2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 3, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidSymSystem sys =
        fluid::fluid_symodes(sn, "softmin", 20.0, std::vector<double>());

    fluid::FluidSymbolicOptions o;
    o.backend = "none";
    const fluid::FluidJacobian jl = fluid::fluid_jacobian(sys, o);
    CHECK(jl.engine == "local");
    // The four outputs are all there, and J is the same matrix the structural
    // differentiation produces: the backend path substitutes SAGE's normal form
    // for these strings, never a different derivative.
    const fluid::FluidSymbolicJacobian ref = fluid::fluid_symbolic_jacobian(sys);
    REQUIRE(jl.vars == ref.vars);
    REQUIRE(jl.rhs == ref.rhs);
    REQUIRE(jl.J == ref.J);
    // Never asked is distinguishable from asked and answered with none.
    CHECK(jl.has_equilibria == false);
    CHECK(jl.equilibria.empty());

    // Solving f(x) = 0 is not differentiating, so with no backend it is REFUSED
    // rather than answered with an empty list, which would read as "there are
    // none". This is the reference's SAGE.require() failure, kept as a failure.
    o.equilibria = true;
    CHECK_THROWS_AS(fluid::fluid_jacobian(sys, o), sym::SymEngineError);

    // The smoothness gate runs BEFORE any backend is contacted, so a min-scaled
    // drift is refused identically whether or not a service is reachable.
    const fluid::FluidSymSystem cl =
        fluid::fluid_symodes(sn, "closing", 20.0, std::vector<double>());
    fluid::FluidSymbolicOptions nb;
    nb.backend = "none";
    CHECK_THROWS_AS(fluid::fluid_jacobian(cl, nb), UnsupportedError);
}
