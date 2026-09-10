/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The arithmetic gate on the model-solve path (WS-A). A closed product-form
 * network is field arithmetic, so it solves exactly under Rational; the
 * open-queue closed forms, DPS-exact, size-based and Marie analyzers evaluate
 * transcendentals, so under Rational they must refuse BY NAME rather than fail
 * to compile or silently degrade. This pins that mva_dispatch<Rational> both
 * compiles and gives the exact rational answer where the algorithm is exact.
 */

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

// A Delay(Z=1) + Queue(mu=2) closed network, N=2: the exact MVA fixed point is
// rational, so Rational must reproduce it bit-for-bit against double.
template <class T>
mva::AvgResult<T> closed_pf(std::size_t N) {
    qn::Network<T> m("pf");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", double(N), d);
    m.set_service(d, c, Distrib<T>::exp_rate(num_traits<T>::from_int(1)));
    m.set_service(q, c, Distrib<T>::exp_rate(num_traits<T>::from_int(2)));
    qn::RoutingMatrix<T> P;
    P.set(c, c, d, q, num_traits<T>::from_int(1));
    P.set(c, c, q, d, num_traits<T>::from_int(1));
    m.link(P);
    mva::MvaOptions opt;
    Matrix<T> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

TEST_CASE("closed product-form solves identically in double and exact") {
    const mva::AvgResult<double> rd = closed_pf<double>(2);
    const mva::AvgResult<Rational> re = closed_pf<Rational>(2);
    // station 1 is the Queue; compare the exact rational to the double answer
    CHECK(num_traits<Rational>::to_double(re.QN(1, 0)) ==
          doctest::Approx(rd.QN(1, 0)).epsilon(1e-12));
    CHECK(num_traits<Rational>::to_double(re.TN(1, 0)) ==
          doctest::Approx(rd.TN(1, 0)).epsilon(1e-12));
    // and that the exact backend really is exact: N=2, Z=1, mu=2 gives
    // X = 6/5 at the queue by hand (rational, no rounding)
    CHECK(re.TN(1, 0) == Rational(6) / Rational(5));
}

TEST_CASE("an open M/M/1 refuses by name under Rational, solves under double") {
    auto build = []() {
        qn::Network<Rational> m("mm1");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C");
        m.set_arrival(s, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
        m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
        qn::RoutingMatrix<Rational> P;
        P.set(c, c, s, q, num_traits<Rational>::from_int(1));
        P.set(c, c, q, k, num_traits<Rational>::from_int(1));
        m.link(P);
        return m;
    };
    qn::Network<Rational> m = build();
    mva::MvaOptions opt;
    Matrix<Rational> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
}

}  // namespace
