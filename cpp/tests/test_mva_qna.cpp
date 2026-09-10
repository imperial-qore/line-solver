/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * QNA, the two-moment open-network decomposition.
 *
 * Numbers are MATLAB's `SolverMVA(model,'method','qna').getAvg()`. The
 * exponential tandem is additionally checkable in closed form: QNA is exact on
 * an M/M/1 chain, so Q1 must be rho/(1-rho) = 1 and Q2 = 0.5.
 *
 * THE GENERAL CASE PINS THE ONE-SWEEP BEHAVIOUR. The reference's open-chain
 * renormalisation produces a NaN convergence measure that stops the iteration
 * after a single sweep, with every flow SCV still at its initial 1 -- so the
 * Erlang/HyperExp tandem reports 0.875 and 0.75, not the 0.75 and (lower)
 * values the converged fixed point would give. Reproducing that is the point of
 * these two rows: they fail if the port ever "fixes" the loop.
 */

#include <vector>

#include "doctest.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/**
 * `HyperExp.fitMeanAndSCV(mean, scv)`, which is `map_hyperexp` at the MATLAB
 * default branching probability p = 0.99 -- NOT the balanced-means fit. The
 * phase rates differ between the two conventions, and although QNA reads only
 * the mean and the SCV, building it the reference's way keeps the model
 * identical rather than merely equivalent for this one analyzer.
 */
D hyperexp_fit(double mean, double scv) {
    const mam::Map<double> m = mam::map_hyperexp<double>(mean, scv);
    const double mu1 = -m.D0(0, 0), mu2 = -m.D0(1, 1);
    return D::hyperexp(m.D1(0, 0) / mu1, mu1, mu2);
}

qn::Network<double> tandem(const D& arrival, const D& s1, const D& s2) {
    qn::Network<double> m("oqn");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, arrival);
    m.set_service(q1, o, s1);
    m.set_service(q2, o, s2);
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, snk, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run_qna(qn::Network<double>& m) {
    mva::MvaOptions opt;
    opt.method = "qna";
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

TEST_CASE("QNA is exact on an exponential open tandem") {
    qn::Network<double> m = tandem(D::exp_rate(1.0), D::exp_rate(2.0), D::exp_rate(3.0));
    const mva::AvgResult<double> r = run_qna(m);
    CHECK(r.actualmethod == "qna");
    // stations are Source, Q1, Q2; the Sink is not a station
    CHECK(r.QN(1, 0) == doctest::Approx(0.999999997500).epsilon(1e-9));
    CHECK(r.QN(2, 0) == doctest::Approx(0.499999999722).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.UN(2, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.999999997500).epsilon(1e-9));
    CHECK(r.RN(2, 0) == doctest::Approx(0.499999999722).epsilon(1e-9));
    CHECK(r.TN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.TN(2, 0) == doctest::Approx(1.0).epsilon(1e-9));
    // a Source holds no jobs
    CHECK(r.QN(0, 0) == doctest::Approx(0.0));
    // the M/M/1 closed forms, which QNA reproduces exactly here
    CHECK(r.QN(1, 0) == doctest::Approx(0.5 / (1.0 - 0.5)).epsilon(1e-8));
    CHECK(r.QN(2, 0) == doctest::Approx((1.0 / 3.0) / (1.0 - 1.0 / 3.0)).epsilon(1e-8));
}

TEST_CASE("QNA on a non-exponential tandem reports its first sweep, as the reference does") {
    qn::Network<double> m = tandem(D::erlang_fit(1.0, 0.5), D::erlang_fit(0.5, 0.5),
                                   hyperexp_fit(1.0 / 3.0, 4.0));
    const mva::AvgResult<double> r = run_qna(m);
    CHECK(r.actualmethod == "qna");
    CHECK(r.QN(1, 0) == doctest::Approx(0.874999998125).epsilon(1e-9));
    CHECK(r.QN(2, 0) == doctest::Approx(0.749999999306).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.UN(2, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.874999998125).epsilon(1e-9));
    CHECK(r.RN(2, 0) == doctest::Approx(0.749999999306).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("QNA refuses a closed chain rather than seeding it from a state it cannot encode") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    opt.method = "qna";
    Matrix<double> init;
    // the whitelist withholds qna on a closed model, so the gate refuses first
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
    CHECK_THROWS_AS(mva::solver_qna(m.get_struct(), opt), UnsupportedError);
}

}  // namespace
