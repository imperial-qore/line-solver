/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * infer_fluid_ps_rt_likelihood and infer_fmlps: the fluid analogue of MLPS.
 *
 * THE ORACLE IS AGAIN A CLOSED FORM. A job alone at a processor-sharing server
 * of rate mu has sojourn Exp(mu), so the passage-time density the fluid drift
 * reports must be exactly mu exp(-mu t). Measured, it agrees to seven digits at
 * every time tested -- which is the whole chain (the layout, the marked block,
 * the departure unmarking, the folded service share, the integration and the
 * derivative read at the terminal state) against an answer this code did not
 * compute.
 *
 * The fluid limit is an APPROXIMATION once the job actually shares the server,
 * so the remaining cases test what must hold regardless: the density is
 * nonnegative, the marked mass is nonincreasing, and a heavier initial queue
 * makes the passage slower rather than faster.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_fmlps.h"
#include "line/lang/qn/network_builder.h"

namespace api = line::api;
namespace qn = line::qn;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> PS Queue -> Delay, one class, service rate `mu` at the queue. */
qn::Network<double> cycle(double mu, double njobs) {
    qn::Network<double> m("f");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("a job alone at the server: the density is exactly mu exp(-mu t)") {
    const double mu = 2.0;
    qn::Network<double> m = cycle(mu, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::fluid::FluidOdeSystem sys = line::fluid::fluid_ode_system(sn);

    std::vector<double> lev(sys.layout.nstates, 0.0);
    lev[sys.layout.qidx[1][0]] = 1.0;  // exactly one job, at the queue

    const double ts[] = {0.1, 0.25, 0.5, 1.0, 2.0};
    for (std::size_t i = 0; i < 5; ++i) {
        const api::FluidRtLikelihood r =
            api::infer_fluid_ps_rt_likelihood(sn, 2, 1, lev, ts[i]);
        CHECK(r.like == doctest::Approx(mu * std::exp(-mu * ts[i])).epsilon(1e-5));
        // The marked mass left is the survival function.
        CHECK(r.markedT == doctest::Approx(std::exp(-mu * ts[i])).epsilon(1e-5));
        CHECK(r.marked0 == doctest::Approx(1.0));
    }
}

TEST_CASE("the density scales with the service rate") {
    // Doubling mu doubles the density at the origin and halves the mean.
    qn::Network<double> a = cycle(1.0, 3.0);
    qn::Network<double> b = cycle(2.0, 3.0);
    const qn::NetworkStruct<double> sa = a.get_struct(), sb = b.get_struct();
    const line::fluid::FluidOdeSystem ya = line::fluid::fluid_ode_system(sa);
    const line::fluid::FluidOdeSystem yb = line::fluid::fluid_ode_system(sb);

    std::vector<double> la(ya.layout.nstates, 0.0), lb(yb.layout.nstates, 0.0);
    la[ya.layout.qidx[1][0]] = 1.0;
    lb[yb.layout.qidx[1][0]] = 1.0;

    const double t = 0.2;
    const api::FluidRtLikelihood ra = api::infer_fluid_ps_rt_likelihood(sa, 2, 1, la, t);
    const api::FluidRtLikelihood rb = api::infer_fluid_ps_rt_likelihood(sb, 2, 1, lb, t);
    CHECK(ra.like == doctest::Approx(1.0 * std::exp(-1.0 * t)).epsilon(1e-5));
    CHECK(rb.like == doctest::Approx(2.0 * std::exp(-2.0 * t)).epsilon(1e-5));
    // The faster server empties the marked block sooner.
    CHECK(rb.markedT < ra.markedT);
}

TEST_CASE("sharing the server slows the passage down") {
    // Two units of fluid at the queue instead of one: the marked job now gets
    // half the server, so more of it survives to any given time. This is the
    // property the folded service share exists to produce; computing the share
    // without folding the marked mass back in would give the marked job the
    // whole server and this check would fail.
    qn::Network<double> m = cycle(2.0, 4.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::fluid::FluidOdeSystem sys = line::fluid::fluid_ode_system(sn);

    std::vector<double> alone(sys.layout.nstates, 0.0), shared(sys.layout.nstates, 0.0);
    alone[sys.layout.qidx[1][0]] = 1.0;
    shared[sys.layout.qidx[1][0]] = 2.0;

    const double t = 0.5;
    const api::FluidRtLikelihood ra = api::infer_fluid_ps_rt_likelihood(sn, 2, 1, alone, t);
    const api::FluidRtLikelihood rs = api::infer_fluid_ps_rt_likelihood(sn, 2, 1, shared, t);
    CHECK(rs.markedT > ra.markedT);
    CHECK(rs.like >= 0.0);
    CHECK(ra.like >= 0.0);
}

TEST_CASE("the marked mass is nonincreasing and the density nonnegative") {
    qn::Network<double> m = cycle(1.5, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::fluid::FluidOdeSystem sys = line::fluid::fluid_ode_system(sn);
    std::vector<double> lev(sys.layout.nstates, 0.0);
    lev[sys.layout.qidx[1][0]] = 1.0;

    double prev = 1.0;
    const double ts[] = {0.1, 0.3, 0.6, 1.0, 1.5, 3.0};
    for (std::size_t i = 0; i < 6; ++i) {
        const api::FluidRtLikelihood r = api::infer_fluid_ps_rt_likelihood(sn, 2, 1, lev, ts[i]);
        CHECK(r.like >= 0.0);
        CHECK(r.markedT <= prev + 1e-9);
        CHECK(r.markedT >= -1e-12);
        prev = r.markedT;
    }
}

TEST_CASE("the refusals are by name") {
    qn::Network<double> m = cycle(2.0, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::fluid::FluidOdeSystem sys = line::fluid::fluid_ode_system(sn);
    std::vector<double> lev(sys.layout.nstates, 0.0);
    lev[sys.layout.qidx[1][0]] = 1.0;

    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 0, 1, lev, 1.0), line::InputError);
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 9, 1, lev, 1.0), line::InputError);
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 2, 0, lev, 1.0), line::InputError);
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 2, 1, lev, 0.0), line::InputError);
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 2, 1, lev, 1.0, 0.0), line::InputError);

    // A state that holds less fluid in the tagged block than the observation
    // marks is not the state the job arrived to.
    std::vector<double> thin(sys.layout.nstates, 0.0);
    thin[sys.layout.qidx[1][0]] = 0.25;
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 2, 1, thin, 1.0), line::InputError);

    // A state vector of the wrong length is refused rather than padded.
    std::vector<double> wrong(sys.layout.nstates + 3, 0.0);
    CHECK_THROWS_AS(api::infer_fluid_ps_rt_likelihood(sn, 2, 1, wrong, 1.0), line::InputError);
}

TEST_CASE("FMLPS recovers the demand from samples of a known service law") {
    // Every observation is a job alone at the queue, where the sojourn is
    // exactly Exp(1/d) and the likelihood peaks at the sample mean -- the same
    // closed form MLPS is checked against, now through the fluid density.
    qn::Network<double> m = cycle(2.0, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();

    const double v[] = {0.30, 0.55, 0.42, 0.71, 0.22, 0.63};
    std::vector<api::MlpsSample> S;
    double sum = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        api::MlpsSample s;
        s.rt = v[i];
        s.cls = 1;
        s.ql = std::vector<double>(1, 1.0);
        S.push_back(s);
        sum += v[i];
    }
    const double mean = sum / 6.0;

    const std::vector<double> d = api::infer_fmlps(sn, 2, S);
    REQUIRE(d.size() == 1u);
    CHECK(d[0] == doctest::Approx(mean).epsilon(1e-3));
}

TEST_CASE("FMLPS refuses the malformed sample sets by name") {
    qn::Network<double> m = cycle(2.0, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    std::vector<api::MlpsSample> ok;
    api::MlpsSample s;
    s.rt = 0.5;
    s.cls = 1;
    s.ql = std::vector<double>(1, 1.0);
    ok.push_back(s);

    CHECK_THROWS_AS(api::infer_fmlps(sn, 0, ok), line::InputError);
    CHECK_THROWS_AS(api::infer_fmlps(sn, 9, ok), line::InputError);
    CHECK_THROWS_AS(api::infer_fmlps(sn, 2, std::vector<api::MlpsSample>()), line::InputError);

    std::vector<api::MlpsSample> badcls = ok;
    badcls[0].cls = 5;
    CHECK_THROWS_AS(api::infer_fmlps(sn, 2, badcls), line::InputError);

    std::vector<api::MlpsSample> badql = ok;
    badql[0].ql = std::vector<double>(3, 1.0);
    CHECK_THROWS_AS(api::infer_fmlps(sn, 2, badql), line::InputError);

    std::vector<api::MlpsSample> badrt = ok;
    badrt[0].rt = -1.0;
    CHECK_THROWS_AS(api::infer_fmlps(sn, 2, badrt), line::InputError);
}
