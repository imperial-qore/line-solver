/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The order-independent (OI) analyzer, the first branch of the dispatch ladder.
 *
 * The single-class model is self-checking: an OI station whose total rate is the
 * CONSTANT 2 is an ordinary M/M/1 queue, so the whole model must give exactly
 * what exact MVA gives for the same Delay(1) + Queue(mu=2) network. MATLAB
 * agrees on both, which is what makes the two-class numbers below trustworthy
 * rather than merely reproducible.
 */

#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

mva::AvgResult<double> run(qn::Network<double>& m) {
    mva::MvaOptions opt;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

TEST_CASE("an OI station with a constant rate IS an M/M/1 queue") {
    qn::Network<double> m("oi1");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service_rate_function(o, [](const std::vector<std::size_t>&) { return 2.0; });
    qn::RoutingMatrix<double> P;
    P.set(d, o, 1.0);
    P.set(o, d, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "oi");
    CHECK(r.QN(0, 0) == doctest::Approx(1.578947368421).epsilon(1e-10));
    CHECK(r.QN(1, 0) == doctest::Approx(1.421052631579).epsilon(1e-10));
    CHECK(r.UN(1, 0) == doctest::Approx(0.789473684211).epsilon(1e-10));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-10));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-10));
    CHECK(r.TN(0, 0) == doctest::Approx(1.578947368421).epsilon(1e-10));
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-10));

    // and the SAME model built as an ordinary FCFS queue must agree exactly,
    // which is the independent check on the whole OI path
    qn::Network<double> ref("ref1");
    const std::size_t dr = ref.add_delay("Delay");
    const std::size_t qr = ref.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t cr = ref.add_closed_class("C1", 3.0, dr);
    ref.set_service(dr, cr, D::exp_rate(1.0));
    ref.set_service(qr, cr, D::exp_rate(2.0));
    qn::RoutingMatrix<double> Pr;
    Pr.set(dr, qr, 1.0);
    Pr.set(qr, dr, 1.0);
    ref.link(Pr);
    const mva::AvgResult<double> rr = run(ref);
    CHECK(rr.actualmethod == "exact");
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.QN(i, 0) == doctest::Approx(rr.QN(i, 0)).epsilon(1e-10));
        CHECK(r.UN(i, 0) == doctest::Approx(rr.UN(i, 0)).epsilon(1e-10));
        CHECK(r.RN(i, 0) == doctest::Approx(rr.RN(i, 0)).epsilon(1e-10));
        CHECK(r.TN(i, 0) == doctest::Approx(rr.TN(i, 0)).epsilon(1e-10));
    }
}

TEST_CASE("a class-dependent OI rate is solved exactly by the CMVA") {
    qn::Network<double> m("oi2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t a = m.add_closed_class("C1", 2.0, d);
    const std::size_t b = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, a, D::exp_rate(1.0));
    m.set_service(d, b, D::exp_rate(2.0));
    // the mean over the queue of the per-position class rate (1+r): permutation
    // invariant, so the station is order-independent
    m.set_service_rate_function(o, [](const std::vector<std::size_t>& cc) {
        double s = 0.0;
        for (std::size_t x : cc) s += 1.0 + static_cast<double>(x);
        return s / static_cast<double>(cc.size());
    });
    qn::RoutingMatrix<double> P;
    // per class: the two-argument shorthand routes class 1 only
    P.set(a, a, d, o, 1.0);
    P.set(a, a, o, d, 1.0);
    P.set(b, b, d, o, 1.0);
    P.set(b, b, o, d, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "oi");
    CHECK(r.QN(0, 0) == doctest::Approx(1.030303030303).epsilon(1e-10));
    CHECK(r.QN(0, 1) == doctest::Approx(0.454545454545).epsilon(1e-10));
    CHECK(r.QN(1, 0) == doctest::Approx(0.969696969697).epsilon(1e-10));
    CHECK(r.QN(1, 1) == doctest::Approx(0.545454545455).epsilon(1e-10));
    // the in-service utilization convention: E[sir_r]/nservers, which can sum
    // above one at an OI station because more than one job may hold a
    // strictly positive rank rate
    CHECK(r.UN(1, 0) == doctest::Approx(0.558441558442).epsilon(1e-10));
    CHECK(r.UN(1, 1) == doctest::Approx(0.545454545455).epsilon(1e-10));
    CHECK(r.RN(1, 0) == doctest::Approx(0.941176470588).epsilon(1e-10));
    CHECK(r.RN(1, 1) == doctest::Approx(0.6).epsilon(1e-10));
    CHECK(r.TN(0, 0) == doctest::Approx(1.030303030303).epsilon(1e-10));
    CHECK(r.TN(0, 1) == doctest::Approx(0.909090909091).epsilon(1e-10));
    // both populations are conserved
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-10));
    CHECK(r.QN(0, 1) + r.QN(1, 1) == doctest::Approx(1.0).epsilon(1e-10));
}

TEST_CASE("a rate function is refused on a discipline that does not take one") {
    qn::Network<double> m("bad");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    m.add_closed_class("C1", 1.0, d);
    CHECK_THROWS_AS(
        m.set_service_rate_function(q, [](const std::vector<std::size_t>&) { return 1.0; }),
        InputError);
}

TEST_CASE("MVA reaches an OI station only through the exact analyzer") {
    qn::Network<double> m("oi1");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service_rate_function(o, [](const std::vector<std::size_t>&) { return 2.0; });
    qn::RoutingMatrix<double> P;
    P.set(d, o, 1.0);
    P.set(o, d, 1.0);
    m.link(P);

    // An AMVA method only ever sees the single-job rates, so it cannot
    // represent the rate function at all; the reference errors rather than
    // reporting a zero queue length at the OI station.
    mva::MvaOptions opt;
    opt.method = "bs";
    Matrix<double> init;
    CHECK_THROWS_AS(mva::mva_dispatch(m.get_struct(), opt, init), UnsupportedError);
    opt.method = "exact";
    CHECK(mva::mva_dispatch(m.get_struct(), opt, init).actualmethod == "oi");
}

}  // namespace
