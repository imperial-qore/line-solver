/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `map2renv` / `mapqn2renv`: the random-environment image of a MAP-modulated
 * network.
 *
 * THE ORACLE IS THE CONSTRUCTION ITSELF, checked field by field against the
 * definition: an MMPP2 service process with intensities (mu0, mu1) and
 * switching rates (s0, s1) must yield exactly 2 stages whose frozen service
 * rates are mu0 and mu1 and whose stage transitions carry Exp(s0) and
 * Exp(s1) -- for an MMPP the modulating chain is preserved exactly, so every
 * number is read off the input, not off another solver.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/io/map2renv.h"
#include "line/lang/qn/network_builder.h"

namespace qn = line::qn;
using line::Matrix;
using line::lang::Distrib;
using line::lang::ProcessType;
using line::lang::SchedStrategy;
using line::num_traits;
using D = Distrib<double>;

namespace {
D mmpp2(double l0, double l1, double s0, double s1) {
    Matrix<double> D0(2, 2), D1(2, 2);
    D1(0, 0) = l0;
    D1(1, 1) = l1;
    D0(0, 0) = -(l0 + s0);
    D0(0, 1) = s0;
    D0(1, 0) = s1;
    D0(1, 1) = -(l1 + s1);
    return D::map_dist(D0, D1, ProcessType::MMPP2);
}
}  // namespace

TEST_CASE("map2renv freezes an MMPP2 service into a 2-stage random environment") {
    qn::Network<double> m("mmpp2_cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    const double mu0 = 2.0, mu1 = 5.0, s0 = 0.3, s1 = 0.7;
    m.set_service(q, c, mmpp2(mu0, mu1, s0, s1));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    line::io::Map2RenvInfo<double> info;
    line::env::Environment<double> env = line::io::map2renv(m.get_struct(), &info);

    CHECK(env.nstages() == 2);
    CHECK(info.nstages == 2);
    CHECK(info.is_mmpp);
    REQUIRE(info.orders.size() == 1);
    CHECK(info.orders[0] == 2);

    // stage service rates are the phase-conditional intensities
    const qn::NetworkStruct<double>& st0 = env.stage(0).model;
    const qn::NetworkStruct<double>& st1 = env.stage(1).model;
    CHECK(st0.rates(1, 0) == doctest::Approx(mu0).epsilon(1e-12));
    CHECK(st1.rates(1, 0) == doctest::Approx(mu1).epsilon(1e-12));
    // the delay is untouched in both stages
    CHECK(st0.rates(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(st1.rates(0, 0) == doctest::Approx(1.0).epsilon(1e-12));

    // the environment switches at the MMPP's own phase rates
    REQUIRE(env.arc(0, 1).enabled);
    REQUIRE(env.arc(1, 0).enabled);
    CHECK(num_traits<double>::to_double(env.arc(0, 1).dist.mean) ==
          doctest::Approx(1.0 / s0).epsilon(1e-12));
    CHECK(num_traits<double>::to_double(env.arc(1, 0).dist.mean) ==
          doctest::Approx(1.0 / s1).epsilon(1e-12));

    // the retained name is the same transformation
    line::env::Environment<double> env2 = line::io::mapqn2renv(m.get_struct());
    CHECK(env2.nstages() == 2);
}

TEST_CASE("map2renv refuses a model with no modulated process") {
    qn::Network<double> m("expcqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    CHECK_THROWS(line::io::map2renv(m.get_struct()));
}
