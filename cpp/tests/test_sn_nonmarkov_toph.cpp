/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * sn_nonmarkov_toph: the solver-side conversion of a non-Markovian service law
 * into a Markovian surrogate. Ported 2026-08-01; before it, SolverMAM refused a
 * Gamma, Weibull, Lognormal, Pareto or Uniform service by name.
 *
 * ORACLES.
 *  (a) Moment fidelity, which is the whole point of the default fit: with
 *      phfit = Cme and an SCV inside (0,1) both moments must survive exactly.
 *  (b) The type tag, which decides what every downstream solver will accept: a
 *      concentrated ME is not a phase-type and must be tagged ME, an Erlang is
 *      and must be tagged PH.
 *  (c) What must NOT change: the Markovian families, and the schedule families
 *      whose whole content would be erased by a homogeneous surrogate.
 */
#include <cmath>
#include <cstddef>

#include "doctest.h"
#include "line/api/sn/sn_nonmarkov_toph.h"
#include "line/lang/qn/network_builder.h"

namespace qn = line::qn;
namespace api = line::api;
using line::lang::ProcessType;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> Queue -> Sink, with the service law supplied by the caller. */
qn::Network<double> open_with_service(const Dist& svc) {
    qn::Network<double> m("nonmkv");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1", src);
    m.set_service(src, c, Dist::exp_rate(0.2));
    m.set_service(q, c, svc);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the default fit keeps both moments of a Gamma service law") {
    // Gamma(4, 0.5): mean 2, SCV 1/4, so it sits inside the range a
    // concentrated ME covers exactly.
    qn::Network<double> m = open_with_service(Dist::gamma_dist(4.0, 0.5));
    qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.service[1][0].type == ProcessType::GAMMA);

    api::sn_nonmarkov_toph(sn);

    const line::lang::Distrib<double>& d = sn.service[1][0];
    CHECK(d.type == ProcessType::ME);  // a CME is not a phase-type
    CHECK(d.has_map());
    CHECK(d.mean == doctest::Approx(2.0).epsilon(1e-8));
    CHECK(d.scv == doctest::Approx(0.25).epsilon(1e-8));
    CHECK(sn.phases_of(2, 1) == d.D0.rows());
    CHECK(d.D0.rows() <= 21);  // within the default budget of 20 plus the tail
}

TEST_CASE("the Bernstein arm returns a genuine phase-type") {
    qn::Network<double> m = open_with_service(Dist::gamma_dist(4.0, 0.5));
    qn::NetworkStruct<double> sn = m.get_struct();

    api::NonmarkovOptions o;
    o.phfit = api::PhFit::Ph;
    api::sn_nonmarkov_toph(sn, o);

    const line::lang::Distrib<double>& d = sn.service[1][0];
    CHECK(d.type == ProcessType::PH);
    CHECK(d.D0.rows() == 20);
    // The density fit pins the mean but not the SCV; that is the trade the
    // reference documents against the two-moment arm.
    CHECK(d.mean == doctest::Approx(2.0).epsilon(1e-8));
    CHECK(d.scv > 0.0);
}

TEST_CASE("a Uniform is fitted despite its density vanishing on the grid") {
    // Uniform(1,3): mean 2, SCV = (4/12)/4 = 1/12.
    qn::Network<double> m = open_with_service(Dist::uniform(1.0, 3.0));
    qn::NetworkStruct<double> sn = m.get_struct();
    api::sn_nonmarkov_toph(sn);
    const line::lang::Distrib<double>& d = sn.service[1][0];
    CHECK(d.has_map());
    CHECK(d.mean == doctest::Approx(2.0).epsilon(1e-8));
    CHECK(d.scv == doctest::Approx(1.0 / 12.0).epsilon(1e-8));

    // Under the density arm the same law must still give a usable fit.
    qn::NetworkStruct<double> s2 = m.get_struct();
    api::NonmarkovOptions o;
    o.phfit = api::PhFit::Ph;
    api::sn_nonmarkov_toph(s2, o);
    CHECK(s2.service[1][0].type == ProcessType::PH);
    CHECK(s2.service[1][0].mean == doctest::Approx(2.0).epsilon(1e-8));
}

TEST_CASE("a Pareto above SCV one takes the density fit even under Cme") {
    // A Pareto has SCV 1/(alpha(alpha-2)), so alpha = 2.2 gives 2.27, outside
    // the range dist_fit_me covers; the Cme arm must fall through to the
    // Bernstein fit rather than refuse.
    qn::Network<double> m = open_with_service(Dist::pareto(2.2, 1.2));
    qn::NetworkStruct<double> sn = m.get_struct();
    const double mean0 = sn.service[1][0].mean;
    REQUIRE(sn.service[1][0].scv > 1.0);

    api::sn_nonmarkov_toph(sn);
    const line::lang::Distrib<double>& d = sn.service[1][0];
    CHECK(d.type == ProcessType::PH);
    CHECK(d.D0.rows() == 20);
    CHECK(d.mean == doctest::Approx(mean0).epsilon(1e-8));
}

TEST_CASE("Det becomes an Erlang of the full budget, and preserveDet keeps it") {
    qn::Network<double> m = open_with_service(Dist::det(1.5));
    qn::NetworkStruct<double> sn = m.get_struct();
    api::sn_nonmarkov_toph(sn);
    CHECK(sn.service[1][0].type == ProcessType::PH);  // a generator, not an ME
    CHECK(sn.service[1][0].D0.rows() == 20);
    CHECK(sn.service[1][0].mean == doctest::Approx(1.5).epsilon(1e-9));
    CHECK(sn.service[1][0].scv == doctest::Approx(1.0 / 20.0).epsilon(1e-9));

    qn::NetworkStruct<double> s2 = m.get_struct();
    api::NonmarkovOptions o;
    o.preserve_det = true;
    api::sn_nonmarkov_toph(s2, o);
    CHECK(s2.service[1][0].type == ProcessType::DET);
}

TEST_CASE("the phase budget is honoured") {
    qn::Network<double> m = open_with_service(Dist::gamma_dist(4.0, 0.5));
    qn::NetworkStruct<double> sn = m.get_struct();
    api::NonmarkovOptions o;
    o.order = 6;
    o.phfit = api::PhFit::Ph;
    api::sn_nonmarkov_toph(sn, o);
    CHECK(sn.service[1][0].D0.rows() == 6);

    CHECK_THROWS_AS(
        [&] {
            qn::NetworkStruct<double> s3 = m.get_struct();
            api::NonmarkovOptions bad;
            bad.order = 0;
            api::sn_nonmarkov_toph(s3, bad);
        }(),
        line::InputError);
}

TEST_CASE("nonmkv = none leaves the struct untouched") {
    qn::Network<double> m = open_with_service(Dist::gamma_dist(4.0, 0.5));
    qn::NetworkStruct<double> sn = m.get_struct();
    api::NonmarkovOptions o;
    o.enabled = false;
    api::sn_nonmarkov_toph(sn, o);
    CHECK(sn.service[1][0].type == ProcessType::GAMMA);
}

TEST_CASE("Markovian and schedule families are left alone") {
    qn::Network<double> m = open_with_service(Dist::erlang_fit(2.0, 0.25));
    qn::NetworkStruct<double> sn = m.get_struct();
    const ProcessType before = sn.service[1][0].type;
    api::sn_nonmarkov_toph(sn);
    CHECK(sn.service[1][0].type == before);
    // The arrival side is an exponential and must survive as declared.
    CHECK(sn.service[0][0].type == ProcessType::EXP);
}
