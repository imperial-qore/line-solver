/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverMAM answers a multiserver station with PHASE-TYPE service exactly.
 *
 * The multiserver fast paths that existed before covered exponential service
 * only. With any other phase-type law the generic path divided the service time
 * by the server count and added a surrogate delay, which reproduces the mean
 * service but not its shape: on the Erlang-2 model below that read 1.176136
 * against the exact 1.075875.
 *
 * THE ORACLE IS THE MATLAB REFERENCE, which reaches these numbers through the
 * same multiset QBD, and independently JMT, which agrees with them to 0.07% and
 * 0.11% at 8e6 samples.
 */
#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_basic.h"

namespace qn = line::qn;
namespace mam = line::mam;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

qn::Network<double> erlang_two_servers() {
    qn::Network<double> m("mphc");
    const std::size_t src = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("K");
    const std::size_t cls = m.add_open_class("C");
    m.set_arrival(src, cls, Dist::exp_rate(0.9));
    // Erlang.fitMeanAndSCV(1.0, 0.5): two phases of rate 2
    m.set_service(q, cls, Dist::erlang(2.0, 2));
    m.set_number_of_servers(q, 2.0);
    qn::RoutingMatrix<double> P;
    P.set(cls, cls, src, q, 1.0);
    P.set(cls, cls, q, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("mam open: phase-type service at c > 1 takes the exact MAP/PH/c path") {
    qn::Network<double> m = erlang_two_servers();
    const line::mva::MvaSolution<double> s = mam::solver_mam_basic(m.get_struct(), mam::MamOptions());
    CHECK(s.Q(1, 0) == doctest::Approx(1.075875).epsilon(1e-5));
}
