/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Wiring of the p-norm smoothing in the `matrix` fluid method.
 *
 * Ruuskanen et al., PEVA 151 (2021) smooth the PS min() of eq. (12) with the
 * inverse p-norm of eq. (26). Two consequences a mean-only comparison hides:
 *
 *  1. An INF station has k = infinity, so min(k, sum x) = sum x identically and
 *     there is nothing to smooth. `sa` carries the total population at a delay,
 *     so smoothing it would run the delay as a k = N queue.
 *  2. Eq. (23) reads the utilization off the SAME share the drift integrated,
 *     k rho / E[sum X] = ghat, so the metrics may not revert to the hard min().
 *
 * The M/M/1 case is the paper's own Example 1: at p = 1 the smoothed model IS
 * the Tipper/PSFFA model, whose fixed point is the exact mean rho/(1-rho),
 * where the unsmoothed mean-field model returns lambda.
 */

#include <cmath>
#include <string>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/solver_fluid.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

std::size_t pnorm_station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

/** Delay(rate 1) + PS Queue(rate 2), one closed class of 10. */
qn::Network<double> pnorm_cqn() {
    qn::Network<double> m("pnorm_cqn");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 10, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

qn::Network<double> pnorm_mm1(double rho) {
    qn::Network<double> m("pnorm_mm1");
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::PS);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, D::exp_rate(rho));
    m.set_service(q, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, s, q, 1.0);
    P.set(c, c, q, k, 1.0);
    m.link(P);
    return m;
}

fluid::FluidOptions pnorm_opts(double pstar) {
    fluid::FluidOptions o;
    o.method = "matrix";
    if (pstar > 0.0) {
        o.pstar = pstar;
        o.pstar_set = true;  // `matrix` keeps the hard min() until asked
    }
    return o;
}

}  // namespace

TEST_CASE("fluid pnorm: the smoothing leaves an INF station alone") {
    // every job at a delay is in service, so its departure rate is Q * mu with
    // no share to apply; smoothing it as a k = N queue returns Q * ghat instead
    qn::Network<double> m = pnorm_cqn();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = pnorm_station_of(sn, "D");
    for (double p : {1.0, 4.0, 20.0}) {
        const fluid::FluidSolution s = fluid::solver_fluid(sn, pnorm_opts(p));
        CHECK(s.TN(id, 0) == doctest::Approx(s.QN(id, 0)).epsilon(1e-9));
    }
    // a large exponent recovers the hard min(), i.e. the unsmoothed fixed point
    const fluid::FluidSolution hard = fluid::solver_fluid(sn, pnorm_opts(0.0));
    const fluid::FluidSolution soft = fluid::solver_fluid(sn, pnorm_opts(20.0));
    CHECK(soft.QN(id, 0) == doctest::Approx(hard.QN(id, 0)).epsilon(1e-6));
}

TEST_CASE("fluid pnorm: the metrics read the smoothed share") {
    // T is read off theta, so a metric taken from the hard min() while x came
    // from the smoothed drift shows up as a flow imbalance around the cycle
    qn::Network<double> m = pnorm_cqn();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = pnorm_station_of(sn, "D"), iq = pnorm_station_of(sn, "Q");
    for (double p : {1.0, 4.0}) {
        const fluid::FluidSolution s = fluid::solver_fluid(sn, pnorm_opts(p));
        CHECK(s.TN(id, 0) == doctest::Approx(s.TN(iq, 0)).epsilon(1e-6));
    }
}

TEST_CASE("fluid pnorm: p = 1 recovers the exact M/M/1 mean") {
    // Ruuskanen et al., Example 1: p = 1 is the Tipper/PSFFA model
    for (double rho : {0.3, 0.5, 0.7}) {
        qn::Network<double> m = pnorm_mm1(rho);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        const std::size_t iq = pnorm_station_of(sn, "Q");
        const fluid::FluidSolution s = fluid::solver_fluid(sn, pnorm_opts(1.0));
        CHECK(s.QN(iq, 0) == doctest::Approx(rho / (1.0 - rho)).epsilon(1e-5));
        CHECK(s.UN(iq, 0) == doctest::Approx(rho).epsilon(1e-6));  // ODE tol
        // the unsmoothed mean-field model returns lambda, the failure the
        // smoothing exists to repair
        const fluid::FluidSolution h = fluid::solver_fluid(sn, pnorm_opts(0.0));
        CHECK(h.QN(iq, 0) == doctest::Approx(rho).epsilon(1e-5));
    }
}
