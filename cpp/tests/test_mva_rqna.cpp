/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * RQNA (WS-B): the robust queueing-network analyzer and the resolveMethod
 * default->rqna upgrade. The model is a single-class open tandem
 * Source -> Q1(FCFS) -> Q2(FCFS) -> Sink with an MMPP2(0.4,0.2,2.0,0.1) bursty
 * arrival. Reference numbers are MATLAB SolverMVA(model,'rqna').getAvgTable
 * (which the plain SolverMVA(model) reproduces, since resolveMethod upgrades a
 * bursty single-class open network to rqna).
 */

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::ProcessType;
using lang::SchedStrategy;

namespace {

// MMPP2(lambda0, lambda1, sigma0, sigma1) as its (D0, D1) MAP.
template <class T>
Distrib<T> mmpp2(double l0, double l1, double s0, double s1) {
    Matrix<T> D0(2, 2), D1(2, 2);
    D1(0, 0) = num_traits<T>::from_double(l0);
    D1(1, 1) = num_traits<T>::from_double(l1);
    D0(0, 0) = num_traits<T>::from_double(-(l0 + s0));
    D0(0, 1) = num_traits<T>::from_double(s0);
    D0(1, 0) = num_traits<T>::from_double(s1);
    D0(1, 1) = num_traits<T>::from_double(-(l1 + s1));
    return Distrib<T>::map_dist(D0, D1, ProcessType::MMPP2);
}

template <class T>
qn::Network<T> tandem_bursty() {
    qn::Network<T> m("rqna");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(src, c, mmpp2<T>(0.4, 0.2, 2.0, 0.1));
    m.set_service(q1, c, Distrib<T>::exp_rate(num_traits<T>::from_int(2)));
    m.set_service(q2, c, Distrib<T>::exp_mean(num_traits<T>::from_double(1.0 / 1.5)));
    qn::RoutingMatrix<T> P;
    P.set(c, c, src, q1, num_traits<T>::from_int(1));
    P.set(c, c, q1, q2, num_traits<T>::from_int(1));
    P.set(c, c, q2, snk, num_traits<T>::from_int(1));
    m.link(P);
    return m;
}

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

TEST_CASE("RQNA matches MATLAB on a bursty-arrival open tandem") {
    qn::Network<double> m = tandem_bursty<double>();
    mva::MvaOptions opt;
    opt.method = "rqna";
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    CHECK(r.actualmethod == "rqna");
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t q1 = station_of(sn, "Q1"), q2 = station_of(sn, "Q2");
    CHECK(r.QN(q1, 0) == doctest::Approx(0.117053).epsilon(1e-5));
    CHECK(r.UN(q1, 0) == doctest::Approx(0.104762).epsilon(1e-5));
    CHECK(r.RN(q1, 0) == doctest::Approx(0.558662).epsilon(1e-5));
    CHECK(r.QN(q2, 0) == doctest::Approx(0.162415).epsilon(1e-5));
    CHECK(r.UN(q2, 0) == doctest::Approx(0.139683).epsilon(1e-5));
    CHECK(r.RN(q2, 0) == doctest::Approx(0.775163).epsilon(1e-5));
    CHECK(r.TN(q1, 0) == doctest::Approx(0.209524).epsilon(1e-5));
}

TEST_CASE("resolveMethod upgrades a bursty single-class open default to rqna") {
    qn::Network<double> m = tandem_bursty<double>();
    mva::MvaOptions opt;  // method = "default"
    Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    CHECK(r.actualmethod == "rqna");
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(r.QN(station_of(sn, "Q2"), 0) == doctest::Approx(0.162415).epsilon(1e-5));
}

TEST_CASE("RQNA refuses by name under exact arithmetic") {
    qn::Network<Rational> m = tandem_bursty<Rational>();
    mva::MvaOptions opt;
    opt.method = "rqna";
    Matrix<Rational> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
}

}  // namespace
