/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The Whitt family as SOLVER METHODS, not as API calls: the abandonment
 * dispatch of SolverMVA, the single-station limits of SolverFLD, the run-length
 * plan and the NHPP test.
 *
 * Every reference below is the API function the solver is supposed to be
 * calling, evaluated on the same inputs, so the test measures the WIRING --
 * which station the metrics land on, whether the carried rate reaches TN and
 * the offered rate AN, whether `default` resolves -- and not the formula, which
 * test_qsys_abandonment.cpp and test_qsys_tvfluid.cpp already cover.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mgisrgi_whitt.h"
#include "line/api/qsys/qsys_ggisgi_fluid.h"
#include "line/api/sim/sim_runlength.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::ImpatienceType;
using lang::ProcessType;
using lang::SchedStrategy;

namespace {

/** Source(Exp(lambda)) -> Queue(Exp(mu), s servers, Exp(theta) patience) -> Sink. */
qn::Network<double> abandoning(double lambda, double mu, double theta, double s) {
    qn::Network<double> m("abandon");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(src, c, Distrib<double>::exp_rate(lambda));
    m.set_service(q, c, Distrib<double>::exp_rate(mu));
    m.set_number_of_servers(q, s);
    m.set_patience(q, c, Distrib<double>::exp_rate(theta), ImpatienceType::RENEGING);
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    return m;
}

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

}  // namespace

TEST_CASE("SolverMVA resolves an exponential-patience station to erlanga") {
    qn::Network<double> m = abandoning(10.0, 1.0, 0.5, 8);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(mva::resolve_method(sn, std::string("default")) == "erlanga");
    // And the method list offers both abandonment names on this shape.
    const std::vector<std::string> valid = mva::list_valid_methods(sn);
    CHECK(std::find(valid.begin(), valid.end(), "erlanga") != valid.end());
    CHECK(std::find(valid.begin(), valid.end(), "mgisrgi") != valid.end());
}

TEST_CASE("SolverMVA erlanga reports the exact Erlang A row") {
    qn::Network<double> m = abandoning(10.0, 1.0, 0.5, 8);
    mva::MvaOptions opt;  // method = "default", which resolves to erlanga
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t q = station_of(sn, "Q");
    const qsys::QsysAbandonResult<double> ref = qsys::qsys_erlanga<double>(10.0, 1.0, 0.5, 8);
    CHECK(r.actualmethod == "erlanga");
    CHECK(r.QN(q, 0) == doctest::Approx(ref.meanNumber).epsilon(1e-9));
    CHECK(r.UN(q, 0) == doctest::Approx(ref.utilization).epsilon(1e-9));
    // The CARRIED rate reaches TN and the OFFERED one reaches AN, which is what
    // makes the loss table read the abandonment as ArvR - Tput.
    CHECK(r.TN(q, 0) == doctest::Approx(ref.throughput).epsilon(1e-9));
    CHECK(r.AN(q, 0) == doctest::Approx(10.0).epsilon(1e-9));
    // R is Little's law on the CARRIED rate, not on lambda.
    CHECK(r.RN(q, 0) == doctest::Approx(ref.meanNumber / ref.throughput).epsilon(1e-9));
}

TEST_CASE("A reneging model is refused by a method that does not declare Reneging") {
    qn::Network<double> m = abandoning(10.0, 1.0, 0.5, 8);
    mva::MvaOptions opt;
    opt.method = "mmk";
    Matrix<double> init;
    CHECK_THROWS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init));
}

TEST_CASE("SolverFLD ggisgi.fluid reports the fluid stationary point") {
    qn::Network<double> m = abandoning(120.0, 1.0, 0.5, 100);
    fluid::FluidOptions opt;
    opt.method = "ggisgi.fluid";
    const fluid::FluidSolution r = fluid::solver_fluid_run_analyzer(m.get_struct(), opt);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t q = station_of(sn, "Q");
    const qsys::QsysFluidAbandonResult<double> ref = qsys::qsys_ggisgi_fluid<double>(
        120.0, 1.0, 100, [](const double& t) { return std::exp(-0.5 * t); });
    CHECK(r.method == "ggisgi.fluid");
    CHECK(r.QN(q, 0) == doctest::Approx(ref.meanNumber).epsilon(1e-9));
    CHECK(r.TN(q, 0) == doctest::Approx(ref.throughput).epsilon(1e-9));
    CHECK(r.UN(q, 0) == doctest::Approx(ref.utilization).epsilon(1e-9));
    // rho = 1.2 > 1, so the fluid model saturates the servers and still has a
    // FINITE queue: 140 waiting-plus-serving against a carried rate of 100.
    CHECK(r.QN(q, 0) == doctest::Approx(140.0).epsilon(1e-9));
    CHECK(r.TN(q, 0) == doctest::Approx(100.0).epsilon(1e-9));
}

TEST_CASE("A fluid limit that needs patience refuses a model without it") {
    qn::Network<double> m("nopatience");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(src, c, Distrib<double>::exp_rate(1.0));
    m.set_service(q, c, Distrib<double>::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    fluid::FluidOptions opt;
    opt.method = "ggisgi.fluid";
    CHECK_THROWS(fluid::solver_fluid_run_analyzer(m.get_struct(), opt));
}

TEST_CASE("sim_runlength_plan recovers the asymptotic variance it was given") {
    // A half-width H at confidence 1-alpha over N samples must give back
    // sigma^2 = (H/z)^2 N, and tightening the relative precision fivefold must
    // cost 25 times the samples.
    Matrix<double> means(1, 1, 2.0);
    Matrix<double> half(1, 1, 0.1);
    const sim::RunLengthPlan<double> p =
        sim::sim_runlength_plan<double>(means, half, 10000.0, 0.01, 0.95);
    const double z = sim::sim_runlength<double>(1.0, 0.0, 1.0, 0.95).z;
    CHECK(p.asymptoticVariance(0, 0) == doctest::Approx((0.1 / z) * (0.1 / z) * 10000.0));
    const sim::RunLengthPlan<double> p5 =
        sim::sim_runlength_plan<double>(means, half, 10000.0, 0.05, 0.95);
    CHECK(p.requiredSamples(0, 0) / p5.requiredSamples(0, 0) == doctest::Approx(25.0).epsilon(1e-9));
    // A non-positive mean or half-width has nothing to plan from.
    Matrix<double> zero(1, 1, 0.0);
    const sim::RunLengthPlan<double> pz =
        sim::sim_runlength_plan<double>(means, zero, 10000.0, 0.05, 0.95);
    CHECK(std::isnan(pz.requiredSamples(0, 0)));
}

TEST_CASE("dist_is_nhpp keeps a Poisson trace and rejects an Erlang renewal one") {
    // Deterministic streams rather than sampled ones: the point is the wiring
    // (inter-arrivals to epochs to the horizon), and a fixed input makes the
    // verdict reproducible. A CONSTANT inter-arrival stream is as far from
    // Poisson as a renewal process gets, and its uniforms are the regular grid,
    // so the KS statistic is bounded well away from zero.
    std::vector<double> det(200, 1.0);
    Distrib<double> d = Distrib<double>::replayer(det);
    const infer::NhppKsResult<double> r = lang::dist_is_nhpp(d);
    CHECK(r.n == 200);
    CHECK(r.pvalue < 0.05);
    // A trace of fewer than two inter-arrivals has no epoch structure to test.
    std::vector<double> one(1, 1.0);
    Distrib<double> d1 = Distrib<double>::replayer(one);
    CHECK_THROWS(lang::dist_is_nhpp(d1));
}
