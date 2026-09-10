/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The SolverNC class surface: the method gate, the multiserver-to-load-
 * dependence conversion, and the two analyzers behind them.
 *
 * The numbers are MATLAB's `[Q,U,R,T,A,W] = SolverNC(model, ...).getAvg()` and
 * `getProbNormConstAggr()` on the same models. That is the right comparison for
 * this layer: the arrival rate and the residence time are computed by the
 * runner and not by any analyzer, so an analyzer-level test cannot see them.
 *
 * A NOTE ON THE TOLERANCES. The default route is not one algorithm. A model
 * that reduces to a single queueing station goes to CoMoM and is EXACT, so it
 * is checked to 1e-9; a multi-station one goes to the Grundmann-Moeller
 * cubature, which is a quadrature stopped on a cost budget, so it is checked to
 * the accuracy MATLAB itself reports there (about 1e-6 against the exact
 * convolution). Checking cub to 1e-9 would be testing the quadrature's
 * rounding, not the port.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/api/mam/map_transform.h"
#include "line/api/me/me_gegec_mql.h"
#include "line/api/me/me_gegecn.h"
#include "line/api/me/me_oqn_blk.h"
#include "line/solvers/nc/solver_nc_cache.h"
#include "line/solvers/nc/solver_nc_cacheqn.h"
#include "line/solvers/nc/solver_nc_retrieval.h"
#include "line/solvers/nc/solver_nc_cdf.h"
#include "line/solvers/nc/solver_nc_prob.h"
#include "line/solvers/nc/solver_nc_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

nc::NcSolverOptions options(const std::string& method) {
    nc::NcSolverOptions opt;
    opt.method = method;
    return opt;
}

mva::AvgResult<double> run(qn::Network<double>& m, const std::string& method = "default") {
    return nc::solver_nc_run_analyzer(m.get_struct(), options(method));
}

/** Delay + single-server FCFS Queue, one closed class of 3 jobs. */
qn::Network<double> closed_single_class() {
    qn::Network<double> m("cqn1");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Delay + two PS queues, two closed classes of 2 and 3 jobs. */
qn::Network<double> closed_multiclass() {
    qn::Network<double> m("cqn2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 3.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(2.0));
    m.set_service(q1, c1, D::exp_rate(2.0));
    m.set_service(q1, c2, D::exp_rate(4.0));
    m.set_service(q2, c1, D::exp_rate(3.0));
    m.set_service(q2, c2, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, d, q1, 0.5);
        P.set(c, c, d, q2, 0.5);
        P.set(c, c, q1, d, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

/** Delay + a 2-server FCFS queue, 4 jobs: the load-dependent route. */
qn::Network<double> closed_multiserver() {
    qn::Network<double> m("cqnms");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    m.set_number_of_servers(q, 2);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue1 -> Queue2 -> Sink, one open class: the open branch. */
qn::Network<double> open_tandem() {
    qn::Network<double> m("oqn");
    const std::size_t s = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, D::exp_rate(0.5));
    m.set_service(q1, o, D::exp_rate(2.0));
    m.set_service(q2, o, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);
    return m;
}

/** A closed chain and an open chain sharing one PS queue: the mixed branch. */
qn::Network<double> mixed_model() {
    qn::Network<double> m("mixed");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    const std::size_t o = m.add_open_class("O1");
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_arrival(s, o, D::exp_rate(0.3));
    m.set_service(q, o, D::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    P.set(o, o, s, q, 1.0);
    P.set(o, o, q, k, 1.0);
    m.link(P);
    return m;
}

/** Two closed classes forming ONE chain through a class switch on the arcs. */
qn::Network<double> class_switching() {
    qn::Network<double> m("cs");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(2.0));
    m.set_service(q, c1, D::exp_rate(3.0));
    m.set_service(q, c2, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c2, d, q, 1.0);
    P.set(c2, c1, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("nc solves a closed single-class product-form network exactly") {
    qn::Network<double> m = closed_single_class();
    const mva::AvgResult<double> r = run(m);
    // The model reduces to one queueing station plus a delay, which is CoMoM's
    // exact repairman branch.
    CHECK(r.actualmethod == "default/comom");
    CHECK(r.QN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(1.421052631578948).epsilon(1e-9));
    CHECK(r.UN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.7894736842105261).epsilon(1e-9));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9000000000000004).epsilon(1e-9));
    CHECK(r.TN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    // arrival rate and residence time, both computed by the runner
    CHECK(r.AN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    CHECK(r.AN(1, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
    CHECK(r.WN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.WN(1, 0) == doctest::Approx(0.9000000000000004).epsilon(1e-9));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(-0.233614851181505).epsilon(1e-9));
}

TEST_CASE("nc dispatches 'ble' to the bias-corrected logistic expansion") {
    // Delay Exp(1) -> PS Exp(2) -> PS Exp(3), one class of 12 jobs. The numbers are
    // MATLAB SolverNC(m,'le'|'ble').getAvgTput / getProbNormConstAggr, and agree
    // with the JAR and Python to 10 digits. The think time puts this on the Z > 0
    // branch, which Laplaces M = 2 directions, so the correction is
    // 2*(1 - log(2 pi)/2); it cancels in G(N-e_r)/G(N), hence the identical
    // throughput.
    qn::Network<double> m("cqn_ble");
    const std::size_t d = m.add_delay("think");
    const std::size_t q1 = m.add_queue("q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("c", 12.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);

    const double lGle = nc::solver_nc_lognormconst(m.get_struct(), options("le"));
    const double lGble = nc::solver_nc_lognormconst(m.get_struct(), options("ble"));
    CHECK(lGle == doctest::Approx(-5.3741456592).epsilon(1e-9));
    CHECK(lGble == doctest::Approx(-5.2120227256).epsilon(1e-9));
    CHECK(lGble - lGle == doctest::Approx(2 * (1.0 - std::log(2 * M_PI) / 2)).epsilon(1e-12));

    const mva::AvgResult<double> rle = run(m, "le");
    const mva::AvgResult<double> rble = run(m, "ble");
    CHECK(rble.actualmethod == "ble");
    CHECK(rle.TN(0, 0) == doctest::Approx(1.9375590305).epsilon(1e-9));
    CHECK(rble.TN(0, 0) == doctest::Approx(rle.TN(0, 0)).epsilon(1e-12));
}

TEST_CASE("nc dispatches 'bkt' to the Stirling-corrected Knessl-Tier expansion") {
    // MATLAB SolverNC(cqn2,'method',m).getProbNormConstAggr():
    //   kt -0.270654206757, bkt -0.339672828397, ble -0.339672828425.
    // With a think time the two corrected expansions are ONE estimator (the R- and
    // (M-1)-dimensional saddle points are one point in dual coordinates), so bkt
    // must land on ble to the accuracy of the two saddle-point solvers.
    qn::Network<double> m = closed_multiclass();
    const double lGkt = nc::solver_nc_lognormconst(m.get_struct(), options("kt"));
    const double lGpp = nc::solver_nc_lognormconst(m.get_struct(), options("bkt"));
    const double lGble = nc::solver_nc_lognormconst(m.get_struct(), options("ble"));
    CHECK(lGkt == doctest::Approx(-0.270654206757).epsilon(1e-8));
    CHECK(lGpp == doctest::Approx(-0.339672828397).epsilon(1e-8));
    CHECK(lGpp == doctest::Approx(lGble).epsilon(1e-8));
    const mva::AvgResult<double> r = run(m, "bkt");
    CHECK(r.actualmethod == "bkt");
    // 'lekt' takes the KT side here (R = 2 <= M = 3) and so IS bkt
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("lekt")) == doctest::Approx(lGpp).epsilon(1e-14));
    CHECK(run(m, "lekt").actualmethod == "lekt");
}

TEST_CASE("the exact methods agree with the default on the single-class model") {
    qn::Network<double> m = closed_single_class();
    for (const std::string& method : {std::string("exact"), std::string("ca"),
                                      std::string("comom")}) {
        INFO("method ", method);
        const mva::AvgResult<double> r = run(m, method);
        CHECK(r.QN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
        CHECK(r.QN(1, 0) == doctest::Approx(1.421052631578947).epsilon(1e-9));
        CHECK(r.TN(1, 0) == doctest::Approx(1.578947368421052).epsilon(1e-9));
        CHECK(nc::solver_nc_lognormconst(m.get_struct(), options(method)) ==
              doctest::Approx(-0.233614851181505).epsilon(1e-9));
    }
}

TEST_CASE("nc solves a multiclass closed network by cubature on the default route") {
    qn::Network<double> m = closed_multiclass();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/cub");
    const double Q[3][2] = {{1.056869207529217, 0.7518565201955865},
                            {0.38705114360974, 0.3092302996974814},
                            {0.5560796488610427, 1.938913180106932}};
    const double U[3][2] = {{1.056869207529217, 0.7518565201955865},
                            {0.2642173018823043, 0.1879641300488966},
                            {0.1761448679215362, 0.7518565201955865}};
    const double R[3][2] = {{1.0, 0.5},
                            {0.7324485203133141, 0.4112889778717871},
                            {1.052314978806249, 2.578834030198405}};
    const double Tp[3][2] = {{1.056869207529217, 1.503713040391173},
                             {0.5284346037646086, 0.7518565201955865},
                             {0.5284346037646086, 0.7518565201955865}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
            INFO("station ", i, " class ", k);
            CHECK(r.QN(i, k) == doctest::Approx(Q[i][k]).epsilon(1e-6));
            CHECK(r.UN(i, k) == doctest::Approx(U[i][k]).epsilon(1e-6));
            CHECK(r.RN(i, k) == doctest::Approx(R[i][k]).epsilon(1e-6));
            CHECK(r.TN(i, k) == doctest::Approx(Tp[i][k]).epsilon(1e-6));
            CHECK(r.AN(i, k) == doctest::Approx(Tp[i][k]).epsilon(1e-6));
        }
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(-0.3191222656741051).epsilon(1e-6));
}

TEST_CASE("the exact convolution on the multiclass model") {
    qn::Network<double> m = closed_multiclass();
    for (const std::string& method : {std::string("exact"), std::string("ca")}) {
        INFO("method ", method);
        const mva::AvgResult<double> r = run(m, method);
        const double Q[3][2] = {{1.056869502469612, 0.7518568744790177},
                                {0.3870510220585489, 0.3092302524353982},
                                {0.5560794754718389, 1.938912873085584}};
        const double Tp[3][2] = {{1.056869502469612, 1.503713748958035},
                                 {0.5284347512348061, 0.7518568744790177},
                                 {0.5284347512348061, 0.7518568744790177}};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t k = 0; k < 2; ++k) {
                INFO("station ", i, " class ", k);
                CHECK(r.QN(i, k) == doctest::Approx(Q[i][k]).epsilon(1e-9));
                CHECK(r.TN(i, k) == doctest::Approx(Tp[i][k]).epsilon(1e-9));
            }
        CHECK(nc::solver_nc_lognormconst(m.get_struct(), options(method)) ==
              doctest::Approx(-0.3191196881340845).epsilon(1e-9));
    }
}

TEST_CASE("a multiserver station is rewritten as load dependence and solved exactly") {
    qn::Network<double> m = closed_multiserver();
    for (const std::string& method : {std::string("default"), std::string("exact")}) {
        INFO("method ", method);
        const mva::AvgResult<double> r = run(m, method);
        // mu(n) = min(n, 2) sends the model to the load-dependent convolution.
        CHECK(r.actualmethod == (method == "default" ? "default/exact/gld" : "exact/gld"));
        CHECK(r.QN(0, 0) == doctest::Approx(2.195744680851064).epsilon(1e-9));
        CHECK(r.QN(1, 0) == doctest::Approx(1.804255319148936).epsilon(1e-9));
        CHECK(r.UN(0, 0) == doctest::Approx(2.195744680851064).epsilon(1e-9));
        // the fraction of the TWO servers busy, not of one
        CHECK(r.UN(1, 0) == doctest::Approx(0.7319148936170212).epsilon(1e-9));
        CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.RN(1, 0) == doctest::Approx(0.8217054263565891).epsilon(1e-9));
        CHECK(r.TN(0, 0) == doctest::Approx(2.195744680851064).epsilon(1e-9));
        CHECK(r.TN(1, 0) == doctest::Approx(2.195744680851064).epsilon(1e-9));
        CHECK(nc::solver_nc_lognormconst(m.get_struct(), options(method)) ==
              doctest::Approx(-1.014305182208116).epsilon(1e-9));
    }
}

TEST_CASE("nc solves an open tandem by the open-class formulas") {
    qn::Network<double> m = open_tandem();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/exact");
    // Source row: no queue length, no utilization, no response time, but it
    // carries the flow. Exact comparisons, not Approx: filter_metric skips the
    // Source for these kinds and `out` starts at zero, so the row is zero BY
    // CONSTRUCTION and any nonzero value is a defect rather than a rounding
    // difference. This is the C++ counterpart of MATLAB zeroSourceMetrics and
    // of the JAR zeroSourceMetrics, both of which had to zero the row after the
    // fact because their solvers wrote into it (measured on an open M/M/1:
    // NC U=1, MAM Q=1, CTMC Q=Inf and R=Inf).
    CHECK(r.QN(0, 0) == 0.0);
    CHECK(r.UN(0, 0) == 0.0);
    CHECK(r.RN(0, 0) == 0.0);
    CHECK(r.WN(0, 0) == 0.0);
    CHECK(r.QN(1, 0) == doctest::Approx(0.3333333333333333).epsilon(1e-9));
    CHECK(r.QN(2, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(r.UN(2, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.6666666666666666).epsilon(1e-9));
    CHECK(r.RN(2, 0) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(r.TN(0, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.TN(2, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.AN(0, 0) == doctest::Approx(0.0));
    CHECK(r.AN(1, 0) == doctest::Approx(0.5).epsilon(1e-9));
    // A purely open model has no closed population to normalize over.
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(0.0));
}

TEST_CASE("nc solves a mixed open and closed network") {
    qn::Network<double> m = mixed_model();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/comom");
    const double Q[3][2] = {{0.0, 0.0},
                            {1.072164948453608, 0.0},
                            {0.9278350515463918, 0.4819587628865978}};
    const double U[3][2] = {{0.0, 0.0},
                            {1.072164948453608, 0.0},
                            {0.5360824742268042, 0.2}};
    const double R[3][2] = {{0.0, 0.0},
                            {1.0, 0.0},
                            {0.8653846153846153, 1.606529209621993}};
    const double Tp[3][2] = {{0.0, 0.3},
                             {1.072164948453608, 0.0},
                             {1.072164948453608, 0.3}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
            INFO("station ", i, " class ", k);
            CHECK(r.QN(i, k) == doctest::Approx(Q[i][k]).epsilon(1e-9));
            CHECK(r.UN(i, k) == doctest::Approx(U[i][k]).epsilon(1e-9));
            CHECK(r.RN(i, k) == doctest::Approx(R[i][k]).epsilon(1e-9));
            CHECK(r.TN(i, k) == doctest::Approx(Tp[i][k]).epsilon(1e-9));
        }
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(0.4158278951437109).epsilon(1e-9));
}

TEST_CASE("nc solves a class-switching chain") {
    qn::Network<double> m = class_switching();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/comom");
    // Each class holds the population only where its own arc leaves it.
    CHECK(r.QN(0, 0) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(r.QN(0, 1) == doctest::Approx(0.0));
    CHECK(r.QN(1, 0) == doctest::Approx(0.0));
    CHECK(r.QN(1, 1) == doctest::Approx(2.0625).epsilon(1e-9));
    CHECK(r.UN(0, 0) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(r.UN(1, 1) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.RN(1, 1) == doctest::Approx(2.2).epsilon(1e-9));
    CHECK(r.TN(0, 0) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(r.TN(1, 1) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(0.980829253011726).epsilon(1e-9));
}

TEST_CASE("a non-Markovian FCFS station drives the outer eta loop") {
    // SolverOptions('NC') sets config.highvar = 'interp', so an FCFS station
    // whose SCV is not 1 has its service time rescaled after each pass and the
    // analyzer re-solves. MATLAB reports "Iterations: 2" on this model, and
    // this is the branch that breaks if the disabled (station, class) pairs are
    // not excluded from the interpolation the way isfinite(SCV) excludes them.
    qn::Network<double> m("mem_check");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, D::exp_rate(0.5));
    m.set_service(q, o, D::erlang(4.0, 2));  // mean 0.5, SCV 0.5
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/exact");
    CHECK(r.QN(1, 0) == doctest::Approx(0.3333333332858712).epsilon(1e-9));
    CHECK(r.RN(1, 0) == doctest::Approx(0.6666666665717425).epsilon(1e-9));
}

TEST_CASE("a fork-join model runs the shared fixed point with nc_dispatch inside") {
    // The MMT transform turns the Fork into a router, the Join into a delay
    // carrying E[max] - mean, and the branches into auxiliary open classes; the
    // NC analyzer then sees a plain mixed network. Numbers are MATLAB's
    // SolverNC getAvg.
    qn::Network<double> m("fj");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    const mva::AvgResult<double> r = run(m);
    REQUIRE(r.QN.rows() == 4);
    REQUIRE(r.QN.cols() == 1);
    // A forked job is present on every branch at once, so the queue lengths of
    // a fork-join model do not partition the population.
    CHECK(r.QN(0, 0) == doctest::Approx(1.011350273353725).epsilon(1e-6));
    CHECK(r.QN(1, 0) == doctest::Approx(0.8716946350488413).epsilon(1e-6));
    CHECK(r.QN(2, 0) == doctest::Approx(0.4725955493011986).epsilon(1e-6));
    CHECK(r.QN(3, 0) == doctest::Approx(0.7313824892199967).epsilon(1e-6));
    CHECK(r.UN(1, 0) == doctest::Approx(0.5056714595317254).epsilon(1e-6));
    CHECK(r.UN(2, 0) == doctest::Approx(0.3371143063544836).epsilon(1e-6));
    CHECK(r.UN(3, 0) == doctest::Approx(0.0));
    // the Join's response time IS the synchronisation delay the fixed point set
    CHECK(r.RN(3, 0) == doctest::Approx(0.3615897612143712).epsilon(1e-6));
    CHECK(r.RN(1, 0) == doctest::Approx(0.8619179692839198).epsilon(1e-6));
    CHECK(r.TN(0, 0) == doctest::Approx(1.011350273353725).epsilon(1e-6));
    CHECK(r.TN(1, 0) == doctest::Approx(1.011342919063451).epsilon(1e-6));
}

TEST_CASE("the deterministic estimators reproduce MATLAB on the multiclass model") {
    // Every method here is an asymptotic expansion or a quadrature, so the
    // numbers below are NOT the exact answer and are not meant to be: they are
    // what MATLAB's own estimator produces, which is what a port has to match.
    // The Monte Carlo family (ls, is, mci, imci, sampling) is left out because
    // its estimate is a function of the generator, not of the model.
    struct Case {
        const char* method;
        const char* actual;
        double Q0, Q1, Q2, X0, X1;
    };
    const Case cases[] = {
        {"cub", "cub", 1.056869, 0.387051, 0.556080, 1.056869, 1.503713},
        {"gm", "gm", 1.056869, 0.387051, 0.556080, 1.056869, 1.503713},
        {"le", "le", 1.034023, 0.390120, 0.575858, 1.034023, 1.436385},
        {"kt", "kt", 1.053844, 0.382114, 0.564042, 1.053844, 1.480229},
        // the Stirling correction is additive in lG and cancels in every ratio, so the
        // mean values are kt's (MATLAB SolverNC(cqn2,'method','bkt') prints the same row)
        {"bkt", "bkt", 1.053844, 0.382114, 0.564042, 1.053844, 1.480229},
        {"lekt", "lekt", 1.053844, 0.382114, 0.564042, 1.053844, 1.480229},
        {"clw", "clw", 1.056869, 0.387052, 0.556079, 1.056869, 1.503714},
        {"propfair", "propfair", 0.646127, 0.278621, 1.075253, 0.646127, 0.619391},
    };
    qn::Network<double> m = closed_multiclass();
    for (const Case& c : cases) {
        // std::string, not the raw const char*: doctest stringifies a pointer as its address
        INFO("method ", std::string(c.method));
        const mva::AvgResult<double> r = run(m, c.method);
        CHECK(r.actualmethod == c.actual);
        CHECK(r.QN(0, 0) == doctest::Approx(c.Q0).epsilon(1e-5));
        CHECK(r.QN(1, 0) == doctest::Approx(c.Q1).epsilon(1e-5));
        CHECK(r.QN(2, 0) == doctest::Approx(c.Q2).epsilon(1e-5));
        CHECK(r.TN(0, 0) == doctest::Approx(c.X0).epsilon(1e-5));
        CHECK(r.TN(0, 1) == doctest::Approx(c.X1).epsilon(1e-5));
    }
}

TEST_CASE("a method outside its own domain returns the reference's zero table") {
    // DELIBERATE REPRODUCTION OF THE REFERENCE CONVENTION, ruled by the user on
    // 2026-07-25 (register row N1). Where a method declines a model -- pana
    // outside normal usage, mmint2/gleint on more than one queueing station --
    // MATLAB warns, returns lG = [], and getAvg renders a table of ZEROS while
    // reporting a completed analysis. The port used to throw; it now matches.
    // The consequence the ruling accepts, and it is a real one: a caller cannot
    // distinguish "the model holds no jobs" from "the solver declined".
    // Do not restore the refusal here without a new ruling.
    qn::Network<double> m = closed_multiclass();
    for (const std::string& method : {std::string("pana"), std::string("mmint2")}) {
        INFO("method ", method);
        const mva::AvgResult<double> r = run(m, method);
        CHECK(r.actualmethod == method);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t k = 0; k < 2; ++k) {
                CHECK(r.QN(i, k) == doctest::Approx(0.0));
                CHECK(r.UN(i, k) == doctest::Approx(0.0));
                CHECK(r.TN(i, k) == doctest::Approx(0.0));
            }
    }
    // comom on a multi-queue model still THROWS, and that is the reference's
    // own gate rather than a port addition: pfqn_nc.m calls line_error there,
    // added deliberately with the comment that "a silently zeroed result is
    // worse than no result". The N1 ruling covers the empty-lG convention only.
    CHECK_THROWS_AS(run(m, "comom"), InputError);
}

// ---------------------------------------------------------------------------
// Tier A: the sojourn-time distribution
// ---------------------------------------------------------------------------

namespace {

/** Delay(a) + FCFS(b), one closed class of N jobs. */
qn::Network<double> delay_fcfs(double a, double b, double N) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", N, d);
    m.set_service(d, c, D::exp_rate(a));
    m.set_service(q, c, D::exp_rate(b));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the sojourn-time CDF reproduces MATLAB, grid and values") {
    // Both columns are MATLAB's: the CDF is reported AT the grid points, so the
    // grid is part of the answer and is asserted alongside the values. These
    // are the values AFTER the getCdfRespT.m grid fix; the grid now ends at
    // 2.25 = T^2 with T = N/rate = 3/2, where it used to end at 9 because the
    // horizon was scaled by the Delay's rate.
    qn::Network<double> m = delay_fcfs(1.0, 2.0, 3.0);
    const nc::CdfRespTResult<double> r =
        nc::solver_nc_cdf_respt(m.get_struct(), options("default"));
    REQUIRE(r.RD.size() == 2u);
    REQUIRE(r.RD[0][0].rows() == 100u);
    REQUIRE(r.RD[1][0].rows() == 100u);
    // grid
    CHECK(r.RD[0][0](0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.RD[0][0](49, 1) == doctest::Approx(1.49386915261219).epsilon(1e-12));
    CHECK(r.RD[0][0](99, 1) == doctest::Approx(2.25).epsilon(1e-12));
    // the Delay's law is its own service distribution
    CHECK(r.RD[0][0](0, 0) == doctest::Approx(0.632120558828558).epsilon(1e-12));
    CHECK(r.RD[0][0](49, 0) == doctest::Approx(0.77549766088019).epsilon(1e-12));
    CHECK(r.RD[0][0](99, 0) == doctest::Approx(0.894600775438136).epsilon(1e-12));
    // the FCFS station's law is the tagged-job passage time
    CHECK(r.RD[1][0](0, 0) == doctest::Approx(0.648128263584807).epsilon(1e-12));
    CHECK(r.RD[1][0](49, 0) == doctest::Approx(0.814256037652956).epsilon(1e-12));
    CHECK(r.RD[1][0](99, 0) == doctest::Approx(0.936400994818563).epsilon(1e-12));
    // a CDF, on both rows
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 1; j < 100; ++j) {
            INFO("station ", i, " point ", j);
            CHECK(r.RD[i][0](j, 0) >= r.RD[i][0](j - 1, 0) - 1e-12);
            CHECK(r.RD[i][0](j, 0) <= 1.0 + 1e-12);
        }
}

TEST_CASE("at one job the passage time IS the service law, exactly") {
    // With a single job there is nothing to queue behind, so the FCFS sojourn
    // time must be Exp(2) at every grid point. This is the check the expected
    // values cannot fake: MATLAB's own deviation from 1-exp(-2t) is 0.
    qn::Network<double> m = delay_fcfs(0.5, 2.0, 1.0);
    const nc::CdfRespTResult<double> r =
        nc::solver_nc_cdf_respt(m.get_struct(), options("default"));
    const Matrix<double>& A = r.RD[1][0];
    REQUIRE(A.rows() == 100u);
    // The corrected horizon is T = N/rate = 1/2 < 1, so the grid DESCENDS:
    // 2*log10(0.5) is negative. That is a property of the reference's
    // logspace(0, 2 log10 T, 100), not of the fix.
    CHECK(A(0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(A(49, 1) == doctest::Approx(0.50351302719175).epsilon(1e-12));
    CHECK(A(99, 1) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(A(0, 0) == doctest::Approx(0.864664716763387).epsilon(1e-12));
    CHECK(A(49, 0) == doctest::Approx(0.634696240752367).epsilon(1e-12));
    CHECK(A(99, 0) == doctest::Approx(0.393469340287367).epsilon(1e-12));
    double worst = 0.0;
    for (std::size_t j = 0; j < 100; ++j)
        worst = std::max(worst, std::fabs(A(j, 0) - (1.0 - std::exp(-2.0 * A(j, 1)))));
    CHECK(worst < 1e-12);
}

TEST_CASE("the rd heuristic agrees with the exact passage time where it must") {
    // On a single-class model the heuristic and the exact algorithm solve the
    // same problem; MATLAB's difference between them is 0 here.
    qn::Network<double> m = delay_fcfs(0.5, 2.0, 1.0);
    nc::NcSolverOptions opt = options("default");
    opt.cdf_algorithm = "rd";
    const nc::CdfRespTResult<double> r = nc::solver_nc_cdf_respt(m.get_struct(), opt);
    CHECK(r.RD[1][0](0, 0) == doctest::Approx(0.864664716763387).epsilon(1e-12));
    CHECK(r.RD[1][0](49, 0) == doctest::Approx(0.634696240752367).epsilon(1e-12));
    CHECK(r.RD[1][0](99, 0) == doctest::Approx(0.393469340287367).epsilon(1e-12));
    CHECK_THROWS_AS(
        [&] {
            nc::NcSolverOptions bad = options("default");
            bad.cdf_algorithm = "nosuch";
            return nc::solver_nc_cdf_respt(m.get_struct(), bad);
        }(),
        UnsupportedError);
}

TEST_CASE("the multiclass sojourn law is reported per class") {
    qn::Network<double> m("cqn2c");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(2.0));
    m.set_service(q, c1, D::exp_rate(3.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
    }
    m.link(P);
    const nc::CdfRespTResult<double> r =
        nc::solver_nc_cdf_respt(m.get_struct(), options("default"));
    // Both classes are served at rate 3 and N = 3, so the corrected horizon is
    // T = 3 * (1/3) = 1 EXACTLY and logspace(0, 0, 100) is a hundred identical
    // points. That is a property of the reference's grid design at T = 1, not
    // of the indexing fix, and the law itself is unaffected: F is the same
    // value at every point because t is.
    CHECK(r.RD[1][0](0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.RD[1][0](99, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.RD[1][0](0, 0) == doctest::Approx(0.821219163588125).epsilon(1e-12));
    CHECK(r.RD[1][0](99, 0) == doctest::Approx(0.821219163588125).epsilon(1e-12));
    CHECK(r.RD[1][1](0, 0) == doctest::Approx(0.853567445976871).epsilon(1e-12));
    CHECK(r.RD[1][1](99, 0) == doctest::Approx(0.853567445976871).epsilon(1e-12));
}

TEST_CASE("the sojourn law declines or refuses, each as the reference does") {
    // NO FCFS STATION: ALIGNED TO MATLAB (empty-result ruling, 2026-07-25).
    // getCdfRespT.m:44-45 warns and returns RD = {} without calling
    // setDistribResults, so this returns an empty result carrying the same
    // warning text rather than throwing. This shape is LESS LOSSY than the
    // zero table the normalizing-constant path returns when a method declines:
    // a caller can detect an empty result by SIZE, whereas zeros cannot be told
    // apart from a model that holds no jobs.
    qn::Network<double> m("ps_only");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const nc::CdfRespTResult<double> none =
        nc::solver_nc_cdf_respt(m.get_struct(), options("default"));
    CHECK(none.RD.empty());
    CHECK(none.warning == "getCdfRespT applies only to FCFS nodes.");
    // AN OPEN CLASS still REFUSES, in both codebases. MATLAB had no gate and
    // failed inside pfqn_stdf with the unspecific MATLAB:nonaninf ("NaN and Inf
    // not allowed"), which names neither the model nor the requirement; an
    // unhandled internal failure is not a specification, so the gate was added
    // to @SolverNC/getCdfRespT.m to match this message rather than this port
    // being aligned down to it.
    qn::Network<double> o = open_tandem();
    CHECK_THROWS_AS(nc::solver_nc_cdf_respt(o.get_struct(), options("default")),
                    UnsupportedError);
}

// ---------------------------------------------------------------------------
// The MEM finite-buffer building blocks (me_gegec_mql, me_gegecn, me_gegecn_pb)
// ---------------------------------------------------------------------------

/**
 * TWO ORACLES WITH DIFFERENT JOBS, and the distinction is the point.
 *
 * The ANALYTIC values pin the ANSWER: at Ca = Cs = 1 the GE blocks degenerate
 * to M/M/1 and M/M/1/N, whose laws are known in closed form, so a match says
 * the algorithm is right on this model. MATLAB's values pin the PORT: they say
 * the transcription is right, and they would agree with a wrong algorithm
 * faithfully copied. Neither alone is sufficient.
 */
TEST_CASE("the GE/GE/c blocks degenerate to the exact M/M/1 and M/M/1/N laws") {
    // ANSWER ORACLE. me_gegec_mql at Ca = Cs = 1 is M/M/1: L = rho/(1-rho).
    CHECK(me::me_gegec_mql<double>(0.5, 1.0, 1.0, 1.0, 1) ==
          doctest::Approx(1.0).epsilon(1e-14));
    CHECK(me::me_gegec_mql<double>(0.25, 1.0, 1.0, 1.0, 1) ==
          doctest::Approx(0.25 / 0.75).epsilon(1e-14));

    // me_gegecn at Ca = Cs = 1 is M/M/1/N: p(n) = rho^n (1-rho)/(1-rho^(N+1)).
    const double rho = 0.5;
    for (long N : {2L, 4L, 6L}) {
        INFO("N = ", N);
        const me::GegecnResult<double> r = me::me_gegecn<double>(0.5, 1.0, 1.0, 1.0, 1, 0, N);
        const double den = 1.0 - std::pow(rho, static_cast<double>(N + 1));
        double L = 0.0, U = 0.0;
        for (long n = 0; n <= N; ++n) {
            const double pn = std::pow(rho, static_cast<double>(n)) * (1.0 - rho) / den;
            CHECK(r.p[static_cast<std::size_t>(n)] == doctest::Approx(pn).epsilon(1e-12));
            L += n * pn;
            if (n > 0) U += pn;
        }
        CHECK(r.L == doctest::Approx(L).epsilon(1e-12));
        CHECK(r.U == doctest::Approx(U).epsilon(1e-12));
        CHECK(r.Lq == doctest::Approx(L - U).epsilon(1e-12));
        // PASTA: a Poisson stream has tau = 1, so every term of (4.3) but n = N
        // vanishes and the blocking probability collapses to p(N). That the
        // formula does this is a check on the batch-overflow term itself.
        CHECK(r.PB == doctest::Approx(r.p[static_cast<std::size_t>(N)]).epsilon(1e-12));
        // and a law, not a collection of numbers
        double tot = 0.0;
        for (double v : r.p) tot += v;
        CHECK(tot == doctest::Approx(1.0).epsilon(1e-14));
    }
}

TEST_CASE("the censored queue obeys the properties it must, reference or not") {
    // These hold for ANY correct implementation and are checked against no
    // reference at all, so they survive a reference being wrong -- which in
    // this one function family has today happened in MATLAB, the JAR and
    // python.

    // (a) L is the mean of the law the function ITSELF returned. This catches a
    // law and a mean drifting apart, which every external comparison misses
    // because it checks them separately against the same source.
    for (double lam : {0.3, 0.8, 1.5}) {
        for (long N : {2L, 5L, 9L}) {
            INFO("lambda = ", lam, ", N = ", N);
            const me::GegecnResult<double> r =
                me::me_gegecn<double>(lam, 2.0, 1.0, 3.0, 1, 0, N);
            double tot = 0.0, L = 0.0, Ebusy = 0.0;
            for (long n = 0; n <= N; ++n) {
                const double pn = r.p[static_cast<std::size_t>(n)];
                CHECK(pn >= -1e-15);          // a probability
                CHECK(pn <= 1.0 + 1e-12);
                tot += pn;
                L += n * pn;
                Ebusy += std::min<long>(n, 1) * pn;
            }
            CHECK(tot == doctest::Approx(1.0).epsilon(1e-14));
            CHECK(r.L == doctest::Approx(L).epsilon(1e-12));
            CHECK(r.U == doctest::Approx(Ebusy).epsilon(1e-12));
            CHECK(r.Lq == doctest::Approx(L - Ebusy).epsilon(1e-12));
            CHECK(r.PB >= 0.0);
            CHECK(r.PB <= 1.0);
        }
    }

    // (b) Blocking is MONOTONE IN THE LOAD at a fixed buffer. More offered work
    // cannot make a full buffer less likely.
    {
        double prev = -1.0;
        for (double lam : {0.2, 0.4, 0.6, 0.9, 1.3, 2.0}) {
            const me::GegecnResult<double> r =
                me::me_gegecn<double>(lam, 2.0, 1.0, 3.0, 1, 0, 5);
            INFO("lambda = ", lam, " PB = ", r.PB);
            CHECK(r.PB > prev);
            prev = r.PB;
        }
    }

    // (c) Blocking VANISHES AS THE BUFFER GROWS, monotonically. A queue with
    // room to spare turns nobody away, and the limit identifies the censoring
    // as censoring rather than as a rescaling of the law.
    {
        double prev = 2.0;
        for (long N : {2L, 4L, 8L, 16L, 32L}) {
            const me::GegecnResult<double> r =
                me::me_gegecn<double>(0.4, 2.0, 1.0, 3.0, 1, 0, N);
            INFO("N = ", N, " PB = ", r.PB);
            CHECK(r.PB < prev);
            prev = r.PB;
        }
        CHECK(me::me_gegecn<double>(0.4, 2.0, 1.0, 3.0, 1, 0, 64).PB < 1e-6);
        // and the mean converges on the UNBOUNDED queue's, which is the other
        // building block -- two independent code paths meeting in the limit.
        const double Linf = me::me_gegec_mql<double>(0.4, 2.0, 1.0, 3.0, 1);
        CHECK(me::me_gegecn<double>(0.4, 2.0, 1.0, 3.0, 1, 0, 64).L ==
              doctest::Approx(Linf).epsilon(1e-6));
    }
}

TEST_CASE("the GE/GE/c blocks reproduce MATLAB away from the exponential case") {
    // PORT ORACLE: genuine GE inputs, where no closed form is available.
    CHECK(me::me_gegec_mql<double>(0.5, 2.0, 1.0, 3.0, 1) ==
          doctest::Approx(2.0).epsilon(1e-12));
    CHECK(me::me_gegec_mql<double>(1.2, 2.0, 1.0, 3.0, 2) ==
          doctest::Approx(3.19864864864865).epsilon(1e-12));
    {
        const me::GegecnResult<double> r = me::me_gegecn<double>(0.8, 2.0, 1.0, 3.0, 2, 0, 5);
        CHECK(r.L == doctest::Approx(1.00550409792819).epsilon(1e-12));
        CHECK(r.U == doctest::Approx(0.373435366959977).epsilon(1e-12));
        CHECK(r.PB == doctest::Approx(0.0664115826000571).epsilon(1e-12));
    }
    {
        // A POSITIVE K, the closed-network censoring the open case never
        // exercises: departures are forbidden from state K, so the law starts
        // at K rather than 0.
        const me::GegecnResult<double> r = me::me_gegecn<double>(0.9, 2.0, 1.0, 2.0, 2, 1, 5);
        CHECK(r.L == doctest::Approx(1.84382019262052).epsilon(1e-12));
        CHECK(r.U == doctest::Approx(0.699927799673999).epsilon(1e-12));
        CHECK(r.PB == doctest::Approx(0.111432001448893).epsilon(1e-12));
        CHECK(r.Lq == doctest::Approx(0.443964593272524).epsilon(1e-12));
        REQUIRE(r.p.size() == 5u);
        const double ref[5] = {0.600144400652, 0.156559408866, 0.106641916184,
                               0.0726401458064, 0.0640141284919};
        double tot = 0.0;
        for (std::size_t i = 0; i < 5; ++i) {
            CHECK(r.p[i] == doctest::Approx(ref[i]).epsilon(1e-11));
            tot += r.p[i];
        }
        CHECK(tot == doctest::Approx(1.0).epsilon(1e-14));
    }
}

TEST_CASE("the GE/GE/c blocks refuse what the GE distribution cannot express") {
    // scv below 1 is not a GE at all; the reference errors rather than
    // approximating, and so does this.
    CHECK_THROWS_AS(me::me_gegecn<double>(0.5, 0.5, 1.0, 1.0, 1, 0, 4), InputError);
    CHECK_THROWS_AS(me::me_gegecn<double>(0.5, 1.0, 1.0, 0.5, 1, 0, 4), InputError);
    CHECK_THROWS_AS(me::me_gegecn<double>(0.5, 1.0, 1.0, 1.0, 1, 4, 4), InputError);
    CHECK_THROWS_AS(me::me_gegecn<double>(0.5, 1.0, 0.0, 1.0, 1, 0, 4), InputError);
    CHECK_THROWS_AS(me::me_gegecn<double>(0.5, 1.0, 1.0, 1.0, 0, 0, 4), InputError);
}

// ---------------------------------------------------------------------------
// Tier 3: the non-reentrant cache (Source-Cache-Sink)
// ---------------------------------------------------------------------------

namespace {

using lang::ReplacementStrategy;

/**
 * Source -> Cache -> Sink, one read class switching to a hit and a miss class.
 *
 * Every cache parameter is a constructor argument so the tests can MOVE ONE and
 * watch a metric respond; see the "honoured, not merely read" case below.
 */
qn::Network<double> cache_model(std::size_t nitems = 5, int cap = 2,
                                std::vector<double> pread = {0.4, 0.25, 0.15, 0.12, 0.08},
                                ReplacementStrategy rs = ReplacementStrategy::RR,
                                bool swap_hit_miss = false) {
    qn::Network<double> m("model");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = nitems;
    ch.itemcap = std::vector<int>{cap};
    ch.replacestrat = rs;
    ch.pread = std::vector<std::vector<double>>{pread, {}, {}};
    ch.hitclass = std::vector<std::size_t>{swap_hit_miss ? 3u : 2u, 0, 0};
    ch.missclass = std::vector<std::size_t>{swap_hit_miss ? 2u : 3u, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t job = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t mis = m.add_open_class("MissClass");
    m.set_arrival(src, job, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(job, job, src, cn, 1.0);
    P.set(hit, hit, cn, snk, 1.0);
    P.set(mis, mis, cn, snk, 1.0);
    m.link(P);
    return m;
}

/** cache_model() with per-item storage costs, and optionally per-list cost caps. */
qn::Network<double> sized_cache_model(std::vector<int> sizes = {1, 1, 2, 2, 3},
                                      std::vector<int> caps = {}) {
    qn::Network<double> m("model");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 5;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double>>{{0.4, 0.25, 0.15, 0.12, 0.08}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2u, 0, 0};
    ch.missclass = std::vector<std::size_t>{3u, 0, 0};
    ch.itemsize = std::move(sizes);
    ch.costcap = std::move(caps);
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t job = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t mis = m.add_open_class("MissClass");
    m.set_arrival(src, job, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(job, job, src, cn, 1.0);
    P.set(hit, hit, cn, snk, 1.0);
    P.set(mis, mis, cn, snk, 1.0);
    m.link(P);
    return m;
}

/** The hit ratio implied by a run: the flow that left as a hit, over the total. */
double hit_ratio(const nc::NcCacheSolution<double>& r) {
    return r.sol.sol.X[1] / (r.sol.sol.X[1] + r.sol.sol.X[2]);
}

}  // namespace

TEST_CASE("the non-reentrant cache reproduces MATLAB on both branches") {
    // Zipf-like popularity over 5 items, capacity 2, RR replacement.
    const double spm_pij0[5] = {0.320433556233686, 0.475385923172859, 0.640589133941376,
                                0.701809392870969, 0.791941248021822};
    const double exa_pij0[5] = {0.346227186052847, 0.489239989103786, 0.652683192590575,
                                0.712339961863252, 0.79950967038954};
    {
        qn::Network<double> m = cache_model();
        const nc::NcCacheSolution<double> r =
            nc::solver_nc_cache_analyzer(m.get_struct(), options("default"));
        CHECK(r.sol.actualmethod == "spm");
        CHECK(r.sol.sol.X[1] == doctest::Approx(0.975232724166409).epsilon(1e-12));
        CHECK(r.sol.sol.X[2] == doctest::Approx(1.02476727583359).epsilon(1e-12));
        for (std::size_t k = 0; k < 5; ++k)
            CHECK(r.pij(k, 0) == doctest::Approx(spm_pij0[k]).epsilon(1e-12));
        // The per-list breakdown is NaN outside the exact branch, and that is a
        // deliberate refusal to report a number: the approximate algorithms
        // derive the miss column and the per-list columns from DIFFERENT
        // expansions, so their breakdown is not a distribution.
        CHECK(std::isnan(r.hitproblist(0, 0)));
        // The per-item table is still the EXACT recursion, because the cache is
        // product form -- so it does NOT equal this branch's own pij.
        for (std::size_t k = 0; k < 5; ++k)
            CHECK(r.itemprob(k, 0) == doctest::Approx(exa_pij0[k]).epsilon(1e-12));
        CHECK(r.itemprob(0, 0) != doctest::Approx(r.pij(0, 0)).epsilon(1e-6));
    }
    {
        qn::Network<double> m = cache_model();
        const nc::NcCacheSolution<double> r =
            nc::solver_nc_cache_analyzer(m.get_struct(), options("exact"));
        CHECK(r.sol.actualmethod == "exact");
        CHECK(r.sol.sol.X[1] == doctest::Approx(0.98371016071915).epsilon(1e-12));
        CHECK(r.sol.sol.X[2] == doctest::Approx(1.01628983928085).epsilon(1e-12));
        for (std::size_t k = 0; k < 5; ++k)
            CHECK(r.pij(k, 0) == doctest::Approx(exa_pij0[k]).epsilon(1e-12));
        CHECK(r.hitproblist(0, 0) == doctest::Approx(0.491855080359575).epsilon(1e-12));
        CHECK(std::isnan(r.hitproblist(1, 0)));  // a class that reads nothing
    }
}

TEST_CASE("'rayint' and 'spm' are aliases, and item sizes pick the tilted kernel") {
    // One method name set -- default / spm / rayint -- choosing the kernel from the
    // model rather than from the name. See _kb/09-ldes-and-cache.md.
    SUBCASE("no item sizes: all three serve the size-free SPM") {
        double ref = 0;
        for (const std::string& method :
             {std::string("default"), std::string("spm"), std::string("rayint")}) {
            INFO("method ", method);
            qn::Network<double> m = cache_model();
            const nc::NcCacheSolution<double> r =
                nc::solver_nc_cache_analyzer(m.get_struct(), options(method));
            CHECK(r.sol.actualmethod == "spm");
            if (method == "default") ref = hit_ratio(r);
            else CHECK(hit_ratio(r) == doctest::Approx(ref).epsilon(1e-14));
        }
    }
    SUBCASE("item sizes: all three serve cache_spm_size and say spm.size") {
        double ref = 0;
        for (const std::string& method :
             {std::string("default"), std::string("spm"), std::string("rayint")}) {
            INFO("method ", method);
            qn::Network<double> m = sized_cache_model();
            const nc::NcCacheSolution<double> r =
                nc::solver_nc_cache_analyzer(m.get_struct(), options(method));
            CHECK(r.sol.actualmethod == "spm.size");
            if (method == "default") ref = hit_ratio(r);
            else CHECK(hit_ratio(r) == doctest::Approx(ref).epsilon(1e-14));
        }
    }
    SUBCASE("cost caps no longer force the exact/sampling switch on a sized cache") {
        // cache_spm_size carries the caps itself, so the default branch keeps them
        // instead of being auto-switched away.
        qn::Network<double> m = sized_cache_model({1, 1, 2, 2, 3}, {4});
        const nc::NcCacheSolution<double> r =
            nc::solver_nc_cache_analyzer(m.get_struct(), options("default"));
        CHECK(r.sol.actualmethod == "spm.size");
        CHECK_FALSE(r.costcap_method_switched);
        // Flow conservation still holds: every read leaves as a hit or a miss.
        CHECK(r.sol.sol.X[1] + r.sol.sol.X[2] == doctest::Approx(2.0).epsilon(1e-12));
    }
}

TEST_CASE("the cache result obeys the identities it must, reference or not") {
    for (const std::string& method : {std::string("default"), std::string("exact")}) {
        INFO("method ", method);
        qn::Network<double> m = cache_model();
        const nc::NcCacheSolution<double> r =
            nc::solver_nc_cache_analyzer(m.get_struct(), options(method));
        // FLOW CONSERVATION: every read leaves either as a hit or as a miss.
        CHECK(r.sol.sol.X[1] + r.sol.sol.X[2] == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(r.sol.sol.Tp(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
        // Each item is somewhere: miss plus the per-list columns sum to 1.
        for (std::size_t k = 0; k < 5; ++k) {
            double tot = 0.0;
            for (std::size_t l = 0; l < r.pij.cols(); ++l) tot += r.pij(k, l);
            INFO("item ", k);
            CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));
        }
        // The miss RATE ought to be the read rate weighted by the miss
        // probability. IT ONLY IS IN THE EXACT BRANCH, and that is a real
        // property of the method rather than a defect: `cache_miss_spm` and
        // `cache_prob_spm` are DIFFERENT expansions, so the approximate branch
        // reports a miss rate and a per-item occupancy that disagree -- here by
        // 4.4% (1.0248 reported against 0.9814 implied). MATLAB does the same
        // to the digit, so this is the method's, not the port's. It is asserted
        // in BOTH directions so that a future change which accidentally made
        // spm self-consistent would be noticed rather than silently welcomed.
        const double pread[5] = {0.4, 0.25, 0.15, 0.12, 0.08};
        double mr = 0.0;
        for (std::size_t k = 0; k < 5; ++k) mr += 2.0 * pread[k] * r.pij(k, 0);
        if (method == "exact") {
            CHECK(r.missrate[0] == doctest::Approx(mr).epsilon(1e-9));
            CHECK(r.sol.sol.X[2] == doctest::Approx(mr).epsilon(1e-9));
        } else {
            CHECK(mr == doctest::Approx(0.9813614007283151).epsilon(1e-9));
            CHECK(r.sol.sol.X[2] == doctest::Approx(1.02476727583359).epsilon(1e-9));
        }
        // A popular item is more likely resident than an unpopular one.
        for (std::size_t k = 1; k < 5; ++k) CHECK(r.pij(k, 0) > r.pij(k - 1, 0));
    }
}

TEST_CASE("every cache key is HONOURED, not merely read") {
    // The reader's unknown-key guard proves a key is READ. It cannot prove the
    // key is USED: the failure this port has already hit once is a key read, a
    // route built, and the split never applied. Each check below MOVES ONE KEY
    // and requires a metric to respond, which is what proves consumption.
    const double base = hit_ratio(
        nc::solver_nc_cache_analyzer(cache_model().get_struct(), options("exact")));
    CHECK(base == doctest::Approx(0.491855080359575).epsilon(1e-9));

    // itemLevelCap: a bigger cache must hit more often.
    const double bigger = hit_ratio(nc::solver_nc_cache_analyzer(
        cache_model(5, 3).get_struct(), options("exact")));
    CHECK(bigger > base + 0.05);
    const double smaller = hit_ratio(nc::solver_nc_cache_analyzer(
        cache_model(5, 1).get_struct(), options("exact")));
    CHECK(smaller < base - 0.05);

    // popularity: concentrating the reads on fewer items must hit more often,
    // at the SAME capacity and item count.
    const double skewed = hit_ratio(nc::solver_nc_cache_analyzer(
        cache_model(5, 2, {0.8, 0.1, 0.05, 0.03, 0.02}).get_struct(), options("exact")));
    const double flat = hit_ratio(nc::solver_nc_cache_analyzer(
        cache_model(5, 2, {0.2, 0.2, 0.2, 0.2, 0.2}).get_struct(), options("exact")));
    CHECK(skewed > base);
    CHECK(flat < base);

    // numItems: spreading the same reads over more items must hit less often.
    const double more_items = hit_ratio(nc::solver_nc_cache_analyzer(
        cache_model(8, 2, {0.3, 0.2, 0.15, 0.1, 0.1, 0.06, 0.05, 0.04}).get_struct(),
        options("exact")));
    CHECK(more_items < base);

    // hitClass / missClass: swapping them must swap WHICH COLUMN carries the
    // flow. A port that read the keys and ignored them would leave this fixed.
    const nc::NcCacheSolution<double> swapped = nc::solver_nc_cache_analyzer(
        cache_model(5, 2, {0.4, 0.25, 0.15, 0.12, 0.08}, ReplacementStrategy::RR, true)
            .get_struct(),
        options("exact"));
    CHECK(swapped.sol.sol.X[2] == doctest::Approx(0.98371016071915).epsilon(1e-9));
    CHECK(swapped.sol.sol.X[1] == doctest::Approx(1.01628983928085).epsilon(1e-9));

    // replacementStrategy: it changes the ADMISSIBILITY of the exact branch,
    // which is the only way this key can show itself here -- RR and FIFO are
    // exchangeable, LRU is not, and the exact recursion must refuse it rather
    // than return the exchangeable answer.
    CHECK_NOTHROW(nc::solver_nc_cache_analyzer(
        cache_model(5, 2, {0.4, 0.25, 0.15, 0.12, 0.08}, ReplacementStrategy::FIFO).get_struct(),
        options("exact")));
    CHECK_THROWS_AS(
        nc::solver_nc_cache_analyzer(
            cache_model(5, 2, {0.4, 0.25, 0.15, 0.12, 0.08}, ReplacementStrategy::LRU)
                .get_struct(),
            options("exact")),
        UnsupportedError);
}

TEST_CASE("the cache analyzer refuses a model outside its domain, by name") {
    // Fewer than capacity+2 items: the reference's own gate. The comparison is
    // against the CAPACITY VALUE, not the list count -- 3 items with capacity 2
    // is short (3 < 4) and must be refused, while 4 items with capacity 2 is
    // not. The port originally compared against the number of lists, which let
    // this model through into a recursion with no headroom; a property test
    // over the boundary is what found it.
    CHECK_THROWS_AS(nc::solver_nc_cache_analyzer(
                        cache_model(3, 2, {0.5, 0.3, 0.2}).get_struct(), options("default")),
                    UnsupportedError);
    CHECK_NOTHROW(nc::solver_nc_cache_analyzer(
        cache_model(4, 2, {0.4, 0.3, 0.2, 0.1}).get_struct(), options("default")));
    // The shape predicate, which the dispatch uses to select this analyzer.
    CHECK(nc::nc_is_noreentrant_cache(cache_model().get_struct()));
    CHECK_FALSE(nc::nc_is_noreentrant_cache(closed_single_class().get_struct()));
    CHECK_FALSE(nc::nc_is_noreentrant_cache(open_tandem().get_struct()));
}

namespace {

/** Delay + Cache in a closed loop: the INTEGRATED caching-queueing shape. */
qn::Network<double> cacheqn_model(std::size_t nitems = 5, int cap = 2,
                                  std::vector<double> pread = {0.2, 0.2, 0.2, 0.2, 0.2}) {
    qn::Network<double> m("cacheqn");
    const std::size_t d = m.add_delay("Delay");
    qn::CacheParam<double> ch;
    ch.nitems = nitems;
    ch.itemcap = std::vector<int>{cap};
    ch.replacestrat = ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double>>{pread, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t job = m.add_closed_class("JobClass", 1.0, d);
    const std::size_t hit = m.add_closed_class("HitClass", 0.0, d);
    const std::size_t mis = m.add_closed_class("MissClass", 0.0, d);
    m.set_service(d, job, D::exp_rate(1.0));
    m.set_service(d, hit, D::exp_rate(1.0));
    m.set_service(d, mis, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(job, job, d, cn, 1.0);
    P.set(hit, job, cn, d, 1.0);
    P.set(mis, job, cn, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the integrated caching-queueing network reproduces MATLAB") {
    {
        qn::Network<double> m = cacheqn_model();
        const nc::NcCacheqnSolution<double> r =
            nc::solver_nc_cacheqn_analyzer(m.get_struct(), options("default"));
        CHECK(r.sol.actualmethod == "spm");
        CHECK(r.sol.sol.iter == 2);
        CHECK(r.hitprob(0, 0) == doctest::Approx(0.394262669203887).epsilon(1e-12));
        CHECK(r.missprob(0, 0) == doctest::Approx(0.605737330796113).epsilon(1e-12));
    }
    {
        qn::Network<double> m = cacheqn_model();
        const nc::NcCacheqnSolution<double> r =
            nc::solver_nc_cacheqn_analyzer(m.get_struct(), options("exact"));
        CHECK(r.sol.actualmethod == "exact");
        CHECK(r.sol.sol.iter == 2);
        // THE ANSWER ORACLE, and it needs no reference: with UNIFORM popularity
        // and RR replacement, a 2-of-5 cache holds each item with probability
        // 2/5 by symmetry, so the hit ratio is EXACTLY 0.4. That the exact
        // branch lands on it to machine precision -- while the spm branch is
        // 1.5% off, as an approximation should be -- says the two branches are
        // doing what their names claim.
        CHECK(r.hitprob(0, 0) == doctest::Approx(0.4).epsilon(1e-12));
        CHECK(r.missprob(0, 0) == doctest::Approx(0.6).epsilon(1e-12));
    }
    // Hit and miss must partition the flow, in both branches.
    for (const std::string& method : {std::string("default"), std::string("exact")}) {
        INFO("method ", method);
        qn::Network<double> m = cacheqn_model();
        const nc::NcCacheqnSolution<double> r =
            nc::solver_nc_cacheqn_analyzer(m.get_struct(), options(method));
        CHECK(r.hitprob(0, 0) + r.missprob(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("the cacheqn hit ratio responds to the keys that determine it") {
    // The moving-metric discipline again, on the integrated shape: the exact
    // branch's hit ratio is cap/nitems under uniform popularity, so it is
    // PREDICTABLE and not merely monotone. That is a stronger check than a
    // direction, and it is available only because the uniform case has a closed
    // form.
    struct Case { std::size_t n; int cap; double expect; };
    const Case cases[] = {{5, 1, 0.2}, {5, 2, 0.4}, {5, 3, 0.6}, {8, 2, 0.25}, {10, 5, 0.5}};
    for (const Case& c : cases) {
        INFO("nitems = ", c.n, " cap = ", c.cap);
        std::vector<double> pr(c.n, 1.0 / static_cast<double>(c.n));
        qn::Network<double> m = cacheqn_model(c.n, c.cap, pr);
        const nc::NcCacheqnSolution<double> r =
            nc::solver_nc_cacheqn_analyzer(m.get_struct(), options("exact"));
        CHECK(r.hitprob(0, 0) == doctest::Approx(c.expect).epsilon(1e-9));
    }
    // Skewing the popularity at a fixed capacity must RAISE the hit ratio above
    // the uniform value, since the cache keeps what is asked for most.
    qn::Network<double> skew = cacheqn_model(5, 2, {0.6, 0.2, 0.1, 0.07, 0.03});
    const nc::NcCacheqnSolution<double> rs =
        nc::solver_nc_cacheqn_analyzer(skew.get_struct(), options("exact"));
    CHECK(rs.hitprob(0, 0) > 0.4);
}

TEST_CASE("the cacheqn analyzer refuses an inadmissible cache, by name") {
    // Here the gate compares against the SUM of the list capacities, where the
    // non-reentrant analyzer compares entrywise. Both are the reference's.
    CHECK_THROWS_AS(nc::solver_nc_cacheqn_analyzer(
                        cacheqn_model(3, 2, {0.4, 0.3, 0.3}).get_struct(), options("default")),
                    UnsupportedError);
    CHECK_NOTHROW(nc::solver_nc_cacheqn_analyzer(
        cacheqn_model(4, 2, {0.25, 0.25, 0.25, 0.25}).get_struct(), options("default")));
}

namespace {

/** Source -> Cache -> {Q1, Q2} -> Cache, an OPEN delayed-hit retrieval cache. */
qn::Network<double> retrieval_model(int cap = 1, double rate1 = 2.0, double rate2 = 3.0) {
    qn::Network<double> m("DelayedHits");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 3;
    ch.itemcap = std::vector<int>{cap};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double>>{{0.6, 0.3, 0.1}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t q1 = m.add_queue("Queue_1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue_2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t job = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t mis = m.add_open_class("MissClass");
    m.set_arrival(src, job, D::exp_rate(1.0));
    m.set_service(q1, job, D::exp_rate(rate1));
    m.set_service(q2, job, D::exp_rate(rate2));
    m.set_retrieval_system(cn, job, mis, std::vector<std::size_t>{q1, q2});
    qn::RoutingMatrix<double> P;
    P.set(job, job, src, cn, 1.0);
    P.set(job, job, cn, q1, 0.5);
    P.set(job, job, cn, q2, 0.5);
    P.set(job, job, q1, cn, 1.0);
    P.set(job, job, q2, cn, 1.0);
    P.set(hit, hit, cn, snk, 1.0);
    P.set(mis, mis, cn, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the open delayed-hit cache reproduces MATLAB") {
    qn::Network<double> m = retrieval_model();
    const nc::NcRetrievalSolution<double> r =
        nc::solver_nc_retrieval_analyzer(m.get_struct(), options("default"));
    CHECK(r.sol.actualmethod == "exact");
    CHECK(r.sol.sol.lG == doctest::Approx(0.214506357918546).epsilon(1e-12));
    CHECK(r.hitprob[0] == doctest::Approx(0.447649788178334).epsilon(1e-12));
    CHECK(r.missprob[0] == doctest::Approx(0.472059713536413).epsilon(1e-12));
    CHECK(r.delayedprob[0] == doctest::Approx(0.0802904982852532).epsilon(1e-12));
    CHECK(r.sol.sol.X[1] == doctest::Approx(0.527940286463587).epsilon(1e-12));
    CHECK(r.sol.sol.X[2] == doctest::Approx(0.472059713536413).epsilon(1e-12));
    // The retrieval stations, read off phi rather than solved.
    CHECK(r.sol.sol.Q(1, 0) == doctest::Approx(0.1234617712326).epsilon(1e-12));
    CHECK(r.sol.sol.Q(2, 0) == doctest::Approx(0.081097437966512).epsilon(1e-12));
    CHECK(r.sol.sol.Tp(1, 0) == doctest::Approx(0.236029856768207).epsilon(1e-12));
    CHECK(r.sol.sol.Tp(2, 0) == doctest::Approx(0.236029856768207).epsilon(1e-12));
    const double ip[3] = {0.342949364535001, 0.605204760944119, 0.847286665321767};
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(r.itemprob(i, 0) == doctest::Approx(ip[i]).epsilon(1e-12));
    CHECK(r.hitproblist(0, 0) == doctest::Approx(0.447649788178334).epsilon(1e-12));
    // The latency deliberately belongs to SolverMVA and is reported as NaN.
    CHECK(std::isnan(r.latency[0]));
}

TEST_CASE("the delayed-hit split obeys its defining identities") {
    qn::Network<double> m = retrieval_model();
    const nc::NcRetrievalSolution<double> r =
        nc::solver_nc_retrieval_analyzer(m.get_struct(), options("default"));
    // THE DEFINING IDENTITY of a delayed-hit cache: every read is a true hit, a
    // delayed hit, or a miss. Nothing else can happen to it.
    CHECK(r.hitprob[0] + r.missprob[0] + r.delayedprob[0] ==
          doctest::Approx(1.0).epsilon(1e-14));
    // A DELAYED HIT LEAVES THROUGH THE HIT CLASS, because the fetch it waits on
    // was already paid for by the miss that started it. So the hit class must
    // carry hit + delayed, not hit alone -- a port that sent delayed hits to
    // the miss class would still conserve flow and still sum to 1.
    CHECK(r.sol.sol.X[1] ==
          doctest::Approx(1.0 * (r.hitprob[0] + r.delayedprob[0])).epsilon(1e-12));
    CHECK(r.sol.sol.X[2] == doctest::Approx(1.0 * r.missprob[0]).epsilon(1e-12));
    // Flow conservation over the two exit classes.
    CHECK(r.sol.sol.X[1] + r.sol.sol.X[2] == doctest::Approx(1.0).epsilon(1e-12));
    // Little's law at each retrieval station.
    for (std::size_t i = 1; i <= 2; ++i)
        CHECK(r.sol.sol.Q(i, 0) ==
              doctest::Approx(r.sol.sol.Tp(i, 0) * r.sol.sol.R(i, 0)).epsilon(1e-12));
}

TEST_CASE("the delayed-hit fractions respond to the fetch speed and the capacity") {
    // Moving-metric checks again, and here they exercise the delayed-hit
    // mechanism specifically rather than the cache.
    const nc::NcRetrievalSolution<double> base =
        nc::solver_nc_retrieval_analyzer(retrieval_model().get_struct(), options("default"));
    // A SLOWER FETCH leaves items in flight longer, so more reads arrive during
    // one and the DELAYED fraction must rise. Nothing else in the model moves.
    const nc::NcRetrievalSolution<double> slow = nc::solver_nc_retrieval_analyzer(
        retrieval_model(1, 0.5, 0.75).get_struct(), options("default"));
    CHECK(slow.delayedprob[0] > base.delayedprob[0]);
    // and a faster one lowers it
    const nc::NcRetrievalSolution<double> fast = nc::solver_nc_retrieval_analyzer(
        retrieval_model(1, 20.0, 30.0).get_struct(), options("default"));
    CHECK(fast.delayedprob[0] < base.delayedprob[0]);
    // A BIGGER CACHE raises the true hit fraction and lowers the miss.
    const nc::NcRetrievalSolution<double> big = nc::solver_nc_retrieval_analyzer(
        retrieval_model(2).get_struct(), options("default"));
    CHECK(big.hitprob[0] > base.hitprob[0]);
    CHECK(big.missprob[0] < base.missprob[0]);
    // The identity survives every one of them.
    for (const nc::NcRetrievalSolution<double>* p : {&slow, &fast, &big})
        CHECK(p->hitprob[0] + p->missprob[0] + p->delayedprob[0] ==
              doctest::Approx(1.0).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// Tier B: exact state probabilities
// ---------------------------------------------------------------------------

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/** Delay(1) + PS(2), one closed class of 3: the smallest probability model. */
qn::Network<double> delay_ps() {
    qn::Network<double> m("cqnP");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Delay(1) -> PS(2) -> PS(4) -> Delay, one closed class of 4. */
qn::Network<double> delay_two_ps() {
    qn::Network<double> m("cqn3");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

/** Delay + two PS queues, two closed classes of 2 and 1. */
qn::Network<double> prob_multiclass() {
    qn::Network<double> m("cqn3s");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t ca = m.add_closed_class("C1", 2.0, d);
    const std::size_t cb = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, ca, D::exp_rate(1.0));
    m.set_service(d, cb, D::exp_rate(2.0));
    m.set_service(q1, ca, D::exp_rate(2.0));
    m.set_service(q1, cb, D::exp_rate(4.0));
    m.set_service(q2, ca, D::exp_rate(3.0));
    m.set_service(q2, cb, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {ca, cb}) {
        P.set(c, c, d, q1, 0.5);
        P.set(c, c, d, q2, 0.5);
        P.set(c, c, q1, d, 1.0);
        P.set(c, c, q2, d, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the aggregate state probability reproduces MATLAB and sums to one") {
    qn::Network<double> m = delay_ps();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const double ref[4] = {0.210526315789474, 0.315789473684211, 0.315789473684211,
                           0.157894736842105};
    double sum = 0.0;
    for (int n = 0; n <= 3; ++n) {
        const double p =
            nc::solver_nc_getprob_aggr(sn, options("default"), 2, std::vector<int>{n}, kNaN);
        INFO("n = ", n);
        CHECK(p == doctest::Approx(ref[n]).epsilon(1e-12));
        sum += p;
    }
    // A distribution, not a collection of numbers that happen to match.
    CHECK(sum == doctest::Approx(1.0).epsilon(1e-12));
    const nc::NcQueueLengthDist<double> md =
        nc::solver_nc_getprob_marg(sn, options("default"), 2);
    REQUIRE(md.P.size() == 4u);
    for (int n = 0; n <= 3; ++n) CHECK(md.P[n] == doctest::Approx(ref[n]).epsilon(1e-12));
}

TEST_CASE("the exact marginal matches an INDEPENDENT product-form enumeration") {
    // Delay(Z=1) -> PS(L1=0.5) -> PS(L2=0.25), N = 4. The product form gives
    //   Pr[n] ~ Z^{n0}/n0! * L1^{n1} * L2^{n2},  n0+n1+n2 = N
    // which is summed here directly, touching no normalizing-constant code at
    // all. This is the check MATLAB's own numbers cannot fake: it is a
    // different formula reaching the same distribution.
    qn::Network<double> m = delay_two_ps();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const int N = 4;
    const double Z = 1.0, L1 = 0.5, L2 = 0.25;
    std::vector<double> w(N + 1, 0.0);
    double G = 0.0;
    for (int n1 = 0; n1 <= N; ++n1)
        for (int n2 = 0; n2 + n1 <= N; ++n2) {
            const int n0 = N - n1 - n2;
            double f = std::pow(L1, n1) * std::pow(L2, n2) * std::pow(Z, n0);
            for (int k = 2; k <= n0; ++k) f /= static_cast<double>(k);
            w[n1] += f;
            G += f;
        }
    const nc::NcQueueLengthDist<double> md =
        nc::solver_nc_getprob_marg(sn, options("default"), 2);
    REQUIRE(md.P.size() == static_cast<std::size_t>(N + 1));
    double sum = 0.0;
    for (int n = 0; n <= N; ++n) {
        INFO("n = ", n);
        CHECK(md.P[n] == doctest::Approx(w[n] / G).epsilon(1e-10));
        sum += md.P[n];
    }
    CHECK(sum == doctest::Approx(1.0).epsilon(1e-12));
    // and against MATLAB's getProbMarg on the same model
    const double ref[5] = {0.181019332161687, 0.249560632688928, 0.274165202108963,
                           0.210896309314587, 0.0843585237258348};
    for (int n = 0; n <= N; ++n) CHECK(md.P[n] == doctest::Approx(ref[n]).epsilon(1e-9));
    // the second queue, whose demand is half the first's
    const double ref2[5] = {0.590509666080844, 0.26713532513181, 0.105448154657293,
                            0.0316344463971881, 0.00527240773286468};
    const nc::NcQueueLengthDist<double> md2 =
        nc::solver_nc_getprob_marg(sn, options("default"), 3);
    for (int n = 0; n <= N; ++n) CHECK(md2.P[n] == doctest::Approx(ref2[n]).epsilon(1e-9));
}

TEST_CASE("the joint and the detailed marginal reproduce MATLAB") {
    qn::Network<double> m = delay_ps();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // The default state: every closed job at its reference station.
    const nc::MarginalState st{{3}, {0}};
    CHECK(nc::solver_nc_getprob_sys(sn, options("default"), st) ==
          doctest::Approx(0.210526315789474).epsilon(1e-12));
    CHECK(nc::solver_nc_getprob_sys_aggr(sn, options("default"), st) ==
          doctest::Approx(0.210526315789474).epsilon(1e-12));
    // solver_nc_marg reports one probability per station; MATLAB's lPr on this
    // model is [-1.55814461804655, -1.15267950993839].
    const nc::NcMargResult<double> mg =
        nc::solver_nc_marg(sn, options("default"), nc::MarginalState{{3}, {2}}, kNaN);
    CHECK(mg.logP[0] == doctest::Approx(-1.55814461804655).epsilon(1e-12));
    CHECK(mg.logP[1] == doctest::Approx(-1.15267950993839).epsilon(1e-12));
    CHECK(mg.P[0] == doctest::Approx(0.210526315789474).epsilon(1e-12));
    CHECK(mg.P[1] == doctest::Approx(0.315789473684211).epsilon(1e-12));
    // getProb returns the LOG probability. THIS IS DELIBERATE, not an
    // oversight, and must not be "fixed": MATLAB's @SolverNC/getProb.m passes
    // solver_nc_marg's first output (lPr) through under the name Pnir without
    // exponentiating, and the user ruled for strict bug-for-bug parity over the
    // safer API (register row N8). exp() of this is 0.31578947368, which is
    // also what getProbAggr returns for the same state -- so the two accessors
    // genuinely disagree in KIND, and that is the reference's behaviour.
    CHECK(nc::solver_nc_getprob(sn, options("default"), 2, std::vector<int>{2}) ==
          doctest::Approx(-1.15267950993839).epsilon(1e-12));
    CHECK(std::exp(nc::solver_nc_getprob(sn, options("default"), 2, std::vector<int>{2})) ==
          doctest::Approx(0.315789473684211).epsilon(1e-12));
}

TEST_CASE("the multiclass state probabilities reproduce MATLAB") {
    qn::Network<double> m = prob_multiclass();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // DO NOT TIGHTEN THESE TO 1e-12. On this model the default NC route is
    // `cub`, whose lG (0.31020865315028889, the same number MATLAB reports) is
    // an APPROXIMATION -- the exact ladder gives 0.31021279699978 -- and the
    // per-partition constants inherit its error. getProbMarg therefore
    // RENORMALIZES, in the reference (`getProbMarg.m`, last block) and here,
    // which is why its last entry (0.0171865054191192) is not its own
    // getProbAggr for the only partition of n=3 (0.0171865054523948): the two
    // differ by exactly the total the vector was divided by. ~1e-9 is the
    // ceiling on any comparison against MATLAB here.
    CHECK(nc::solver_nc_getprob_aggr(sn, options("default"), 2, std::vector<int>{1, 0}, kNaN) ==
          doctest::Approx(0.229153405173061).epsilon(1e-8));
    CHECK(nc::solver_nc_getprob_aggr(sn, options("default"), 2, std::vector<int>{0, 1}, kNaN) ==
          doctest::Approx(0.0636537238977584).epsilon(1e-8));
    CHECK(nc::solver_nc_getprob_aggr(sn, options("default"), 2, std::vector<int>{2, 1}, kNaN) ==
          doctest::Approx(0.0171865054523948).epsilon(1e-8));
    const nc::NcQueueLengthDist<double> md =
        nc::solver_nc_getprob_marg(sn, options("default"), 2);
    const double ref[4] = {0.590706556988734, 0.292807128503902, 0.0992998090882443,
                           0.0171865054191192};
    double sum = 0.0;
    for (int n = 0; n <= 3; ++n) {
        CHECK(md.P[n] == doctest::Approx(ref[n]).epsilon(1e-8));
        sum += md.P[n];
    }
    // The renormalization is what makes this hold to the last bit rather than
    // to the constant's own accuracy.
    CHECK(sum == doctest::Approx(1.0).epsilon(1e-12));
    // The only partition of n = 3 is [2,1], so the marginal entry is that
    // aggregate probability DIVIDED BY THE TOTAL -- equal to within the
    // renormalization and not beyond it, exactly as in the reference.
    CHECK(md.P[3] ==
          doctest::Approx(nc::solver_nc_getprob_aggr(sn, options("default"), 2,
                                                     std::vector<int>{2, 1}, kNaN))
              .epsilon(1e-8));
    const nc::MarginalState st{{2, 1}, {0, 0}, {0, 0}};
    CHECK(nc::solver_nc_getprob_sys(sn, options("default"), st) ==
          doctest::Approx(0.183322724825544).epsilon(1e-8));
    CHECK(nc::solver_nc_getprob_sys_aggr(sn, options("default"), st) ==
          doctest::Approx(0.183323484042729).epsilon(1e-12));
}

TEST_CASE("the state probability refuses what a marginal cannot carry, by name") {
    const qn::NetworkStruct<double>& sn = delay_ps().get_struct();
    // A malformed marginal is an input error, not a wrong number.
    CHECK_THROWS_AS(nc::solver_nc_margaggr(sn, options("default"),
                                           nc::MarginalState{{3}}, kNaN),
                    InputError);
    CHECK_THROWS_AS(
        nc::solver_nc_getprob_aggr(sn, options("default"), 9, std::vector<int>{1}, kNaN),
        InputError);
    // An open class has no closed state space to normalize over.
    qn::Network<double> o = open_tandem();
    CHECK_THROWS_AS(nc::solver_nc_getprob_aggr(o.get_struct(), options("default"), 2,
                                               std::vector<int>{1}, kNaN),
                    UnsupportedError);
    // SIRO's detailed marginal needs the class of the job IN SERVICE, which a
    // per-class marginal does not carry. The aggregate route has no such
    // dependency and is not refused.
    qn::Network<double> s("siro");
    const std::size_t d = s.add_delay("Delay");
    const std::size_t q = s.add_queue("Queue", SchedStrategy::SIRO);
    const std::size_t c = s.add_closed_class("C1", 2.0, d);
    s.set_service(d, c, D::exp_rate(1.0));
    s.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    s.link(P);
    CHECK_THROWS_AS(nc::solver_nc_getprob(s.get_struct(), options("default"), 2,
                                          std::vector<int>{1}),
                    UnsupportedError);
    CHECK_NOTHROW(nc::solver_nc_getprob_aggr(s.get_struct(), options("default"), 2,
                                             std::vector<int>{1}, kNaN));
}

TEST_CASE("MEM with a finite buffer, against CTMC for the answer and MATLAB for the port") {
    // BOTH ORACLES, IN THEIR STATED ROLES. The CTMC values are the EXACT
    // M/M/1/N answer and say the algorithm is right on this model; MATLAB's MEM
    // values say the transcription is right and would agree with a faithfully
    // copied wrong algorithm. The tolerances are the ones the reference test
    // (`line-test.git .../test_me_oqn_blk.m` Test 6) applies to its own MEM
    // result: 1e-6 on the queue length, 1e-5 on the throughput.
    struct Row {
        double Nb, ctmcQ, ctmcT, memQ, memT;
    };
    const Row rows[3] = {
        {2.0, 0.85245902, 0.59016393, 0.852459016393443, 0.590164734887295},
        {4.0, 1.56306521, 0.70252261, 1.56306520704426, 0.702524095668729},
        {6.0, 2.14243372, 0.74692668, 2.14243371503539, 0.746927487407071},
    };
    for (const Row& row : rows) {
        INFO("buffer N = ", row.Nb);
        qn::Network<double> m("mm1n");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue1", SchedStrategy::FCFS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t o = m.add_open_class("Class1");
        m.set_arrival(s, o, D::exp_rate(0.8));
        m.set_service(q, o, D::exp_rate(1.0));
        m.set_capacity(q, row.Nb);
        qn::RoutingMatrix<double> P;
        P.set(s, q, 1.0);
        P.set(q, k, 1.0);
        m.link(P);

        const nc::NcSolution<double> r = nc::solver_nc_mem(m.get_struct(), options("mem"));
        // The finite buffer selects the censored blocks, which the reported
        // method name has to say -- a caller cannot otherwise tell which of the
        // four MEM algorithms ran.
        CHECK(r.actualmethod == "mem.blocking");
        // ANSWER: against the exact CTMC.
        CHECK(std::fabs(r.sol.Q(1, 0) - row.ctmcQ) < 1e-6);
        CHECK(std::fabs(r.sol.Tp(1, 0) - row.ctmcT) < 1e-5);
        // PORT: against MATLAB's own MEM, to every printed digit.
        CHECK(r.sol.Q(1, 0) == doctest::Approx(row.memQ).epsilon(1e-12));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(row.memT).epsilon(1e-12));
        // Utilization is the carried flow over the unit service rate.
        CHECK(r.sol.U(1, 0) == doctest::Approx(r.sol.Tp(1, 0)).epsilon(1e-12));
    }
}

TEST_CASE("the blocking network reproduces MATLAB on loss and on BAS") {
    // Two stations in tandem, the second one finite. Loss and transfer blocking
    // differ in what happens to a job that finds it full, so the two runs must
    // NOT agree: under BAS the job is held in station 1's server, which is why
    // its queue length rises from 3.15 to 8.44 and the carried flow rises to
    // the full 0.6 rather than being thinned to 0.496.
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    const std::vector<double> l0{0.6, 0.0}, ca0{2.0, 1.0}, mu{1.0, 1.0}, cs{3.0, 2.0};
    {
        const me::MeBlkResult<double> r =
            me::me_oqn_blk<double>(2, l0, ca0, mu, cs, P, {1L, 1L}, {0L, 4L}, {0, 0});
        CHECK(r.Q[0] == doctest::Approx(3.15).epsilon(1e-12));
        CHECK(r.Q[1] == doctest::Approx(1.1433205664861).epsilon(1e-12));
        CHECK(r.T_[1] == doctest::Approx(0.496436304409967).epsilon(1e-12));
        CHECK(r.PBa[1] == doctest::Approx(0.172606159316721).epsilon(1e-12));
        // Under LOSS the carried flow is thinned by the blocking probability.
        CHECK(r.T_[1] == doctest::Approx(r.T_[0] * (1.0 - r.PBa[1])).epsilon(1e-9));
    }
    {
        const me::MeBlkResult<double> r =
            me::me_oqn_blk<double>(2, l0, ca0, mu, cs, P, {1L, 1L}, {0L, 4L}, {0, 1});
        CHECK(r.Q[0] == doctest::Approx(8.44040595418807).epsilon(1e-12));
        CHECK(r.Q[1] == doctest::Approx(1.49927623037146).epsilon(1e-12));
        CHECK(r.PBa[1] == doctest::Approx(0.245763735293554).epsilon(1e-12));
        // Under BAS NOTHING IS LOST: every job eventually enters, so the second
        // station carries the whole offered flow. That is the property which
        // distinguishes the two policies, and it holds exactly rather than to
        // the fixed point's tolerance.
        CHECK(r.T_[1] == doctest::Approx(0.6).epsilon(1e-12));
        CHECK(r.T_[0] == doctest::Approx(r.T_[1]).epsilon(1e-12));
    }
}

TEST_CASE("MEM solves the mixed model, whose gate used to make it unreachable") {
    // This branch was DEAD until 2026-07-25: `sn_get_buffer_size` compared the
    // capacity against sum(njobs), which a mixed model's open class makes
    // infinite, so the Delay's classcap (the closed chain's population, 2) read
    // as a binding buffer and solver_nc_mem_supports rejected every mixed model
    // with "finite station buffers only in open models". Register row N10; the
    // MATLAB gate now compares against the population that can REACH the
    // station. A branch that has never run is not an oracle, so the checks
    // below are INVARIANTS as well as MATLAB's numbers.
    qn::Network<double> m("mem_mixed");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t cc = m.add_closed_class("C1", 2.0, d);
    const std::size_t oo = m.add_open_class("O1");
    m.set_service(d, cc, D::exp_rate(1.0));
    m.set_service(q, cc, D::exp_rate(2.0));
    m.set_arrival(s, oo, D::exp_rate(0.3));
    // HyperExp of mean 0.5, scv 3, built through map_hyperexp as MATLAB's
    // HyperExp.fitMeanAndSCV does.
    const mam::Map<double> hyp = mam::map_hyperexp(0.5, 3.0);
    m.set_service(q, oo, D::map_dist(hyp.D0, hyp.D1, lang::ProcessType::HYPEREXP));
    qn::RoutingMatrix<double> P;
    P.set(cc, cc, d, q, 1.0);
    P.set(cc, cc, q, d, 1.0);
    P.set(oo, oo, s, q, 1.0);
    P.set(oo, oo, q, k, 1.0);
    m.link(P);

    nc::NcSolverOptions opt = options("mem");
    const nc::NcSolution<double> r = nc::solver_nc_mem(m.get_struct(), opt);
    CHECK(r.actualmethod == "mem");
    // MATLAB's values, every printed digit.
    CHECK(r.sol.Q(1, 0) == doctest::Approx(1.10735826296743).epsilon(1e-12));
    CHECK(r.sol.Q(2, 0) == doctest::Approx(0.892641737032569).epsilon(1e-12));
    CHECK(r.sol.Q(2, 1) == doctest::Approx(0.384094940750727).epsilon(1e-12));
    CHECK(r.sol.U(2, 0) == doctest::Approx(0.553679131483715).epsilon(1e-12));
    CHECK(r.sol.U(2, 1) == doctest::Approx(0.15).epsilon(1e-12));
    CHECK(r.sol.R(1, 0) == doctest::Approx(1.00000073802581).epsilon(1e-12));
    CHECK(r.sol.R(2, 1) == doctest::Approx(1.28031646916909).epsilon(1e-12));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(1.10735744570905).epsilon(1e-12));
    CHECK(r.sol.Tp(0, 1) == doctest::Approx(0.3).epsilon(1e-12));

    // INVARIANTS, which is what makes a first-ever output trustworthy.
    // The closed class conserves its population across the two stations.
    CHECK(r.sol.Q(1, 0) + r.sol.Q(2, 0) == doctest::Approx(2.0).epsilon(1e-9));
    // The open class carries exactly its arrival rate through the queue.
    CHECK(r.sol.Tp(2, 1) == doctest::Approx(0.3).epsilon(1e-12));
    // Little's law at the queue, per class.
    CHECK(r.sol.Q(2, 0) ==
          doctest::Approx(r.sol.Tp(2, 0) * r.sol.R(2, 0)).epsilon(1e-9));
    CHECK(r.sol.Q(2, 1) ==
          doctest::Approx(r.sol.Tp(2, 1) * r.sol.R(2, 1)).epsilon(1e-9));
    // Utilization is the carried rate times the mean service time -- but only
    // to 7.4e-7 for the CLOSED class, because U comes from the converged rho
    // and Tp from lambda, and the two differ by the fixed point's residual.
    // That is the SAME 7e-7 the product-form recovery shows on the exponential
    // model above, which is what identifies it as convergence rather than a
    // defect. The open class carries an exact arrival rate and needs no slack.
    CHECK(r.sol.U(2, 0) == doctest::Approx(r.sol.Tp(2, 0) * 0.5).epsilon(1e-5));
    CHECK(std::fabs(r.sol.U(2, 0) - r.sol.Tp(2, 0) * 0.5) / r.sol.U(2, 0) < 1e-6);
    CHECK(r.sol.U(2, 1) == doctest::Approx(r.sol.Tp(2, 1) * 0.5).epsilon(1e-12));
}

TEST_CASE("the method whitelist is the solver's, and what is unported is named") {
    qn::Network<double> m = closed_single_class();
    const std::vector<std::string> valid = nc::list_valid_methods();
    auto has = [&](const std::string& s) {
        return std::find(valid.begin(), valid.end(), s) != valid.end();
    };
    CHECK(has("default"));
    CHECK(has("exact"));
    CHECK(has("comom"));
    CHECK(has("cub"));
    CHECK(has("le"));
    // 'ms' names the lossn_manjunath transform, reachable only on a loss-network
    // shape; see test_nc_lossn.cpp for it actually running.
    CHECK(has("ms"));
    // ANOTHER SITE THAT ONCE PINNED AN ABSENCE: 'rgf' and 'panald' had no
    // port anywhere in this tree and were asserted absent. Both landed, so the
    // absence is now the wrong claim -- 'rgf' answers this model by its own
    // recursion and reports 'rgf', and 'panald' refuses BY NAME when the model
    // is outside normal usage, which is the reference's own refusal.
    CHECK(has("rgf"));
    CHECK(has("panald"));
    CHECK_FALSE(has("nosuchmethod"));
    CHECK_THROWS_AS(run(m, "nosuchmethod"), UnsupportedError);
    // THIS SITE ONCE PINNED AN ABSENCE. Until MEM landed it asserted
    // CHECK_THROWS_AS(run(m,"mem"), UnsupportedError); the refusal was correct
    // then and wrong the moment solver_nc_mem.h existed, so it was replaced
    // with the reference's numbers under approval on 2026-07-25.
    //
    // The values are MATLAB's SolverNC(m,'method','mem').getAvg() on this model.
    // What makes them worth asserting: the model is ENTIRELY EXPONENTIAL, so
    // MEM's fixed point must converge on the product form, and it does -- the
    // exact queue lengths are [1.578947368421052, 1.421052631578947] and the
    // Delay's response time is exactly 1. The residual 7e-7 is the fixed point
    // converging, not a transcription, which is what a pure golden could not
    // tell you.
    {
        const mva::AvgResult<double> mem = run(m, "mem");
        CHECK(mem.actualmethod == "mem");
        CHECK(mem.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
        CHECK(mem.QN(1, 0) == doctest::Approx(1.42105263157895).epsilon(1e-12));
        CHECK(mem.UN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
        CHECK(mem.UN(1, 0) == doctest::Approx(0.789473684210526).epsilon(1e-12));
        CHECK(mem.RN(0, 0) == doctest::Approx(1.00000075136613).epsilon(1e-12));
        CHECK(mem.RN(1, 0) == doctest::Approx(0.900000676229517).epsilon(1e-12));
        CHECK(mem.TN(0, 0) == doctest::Approx(1.57894618205437).epsilon(1e-12));
        CHECK(mem.TN(1, 0) == doctest::Approx(1.57894618205437).epsilon(1e-12));
        // the product-form recovery, to the accuracy the fixed point reaches
        CHECK(mem.QN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-6));
        CHECK(mem.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-5));
    }
    // PERMANENT for THIS model, unlike the mem site above. `erlangfp` is the
    // Erlang fixed point of a LOSS NETWORK, reachable only through a Finite
    // Capacity Region under a DROP rule, and `m` declares no region -- so the
    // refusal is a property of the model, not a missing port. The analyzer and
    // all three of its methods are ported; see test_nc_lossn.cpp.
    CHECK_THROWS_AS(run(m, "erlangfp"), UnsupportedError);
    // THE THIRD SITE THAT PINNED AN ABSENCE, and the last: the load-dependent
    // Norlund-Rice family used to be whitelisted-and-refused. All three run
    // now, so what is asserted is the reference's numbers, each method naming
    // ITSELF rather than silently answering as 'default'. Values are MATLAB
    // SolverNC(m,'method',<name>).getAvg() on this model, whose exact queue
    // lengths are [1.578947368421052, 1.421052631578947].
    {
        const mva::AvgResult<double> rd = run(m, "rd");
        CHECK(rd.actualmethod == "rd");
        CHECK(rd.QN(0, 0) == doctest::Approx(1.65425582052974).epsilon(1e-11));
        CHECK(rd.QN(1, 0) == doctest::Approx(1.34574417947026).epsilon(1e-11));

        const mva::AvgResult<double> nrp = run(m, "nrp");
        CHECK(nrp.actualmethod == "nrp");
        CHECK(nrp.QN(0, 0) == doctest::Approx(1.57894780666713).epsilon(1e-11));
        CHECK(nrp.QN(1, 0) == doctest::Approx(1.42105219333287).epsilon(1e-11));

        // NRL IS LOOSE AGAINST MATLAB BY CONSTRUCTION, not by a porting defect.
        // laplaceapprox differentiates log h TWICE by central differences at
        // step 1e-5, so a last-bit disagreement in the integrand is amplified
        // by 1/(2h)^2 ~ 2.5e9: infradius_h agrees with MATLAB's to 2 ulp at
        // every point checked, and det(-H) still parts company in the 7th
        // digit. That is the method's own conditioning at this step, and it
        // sets the tolerance here. pfqn_nrp is unaffected because its Hessian
        // comes out of the normal substitution bit-for-bit.
        const mva::AvgResult<double> nrl = run(m, "nrl");
        CHECK(nrl.actualmethod == "nrl");
        CHECK(nrl.QN(0, 0) == doctest::Approx(1.57894824491454).epsilon(1e-5));
        CHECK(nrl.QN(1, 0) == doctest::Approx(1.42105175508546).epsilon(1e-5));

        // 'rgf' USED TO RESOLVE TO THE EXACT CONVOLUTION HERE and reported
        // 'rgf/ca'. The name a single-class model reports is decided by the
        // arrival-theorem call, which replicates the tagged station into an
        // auxiliary class and so asks pfqn_nc a TWO-class question; until the
        // multiclass residue recursion landed there was nothing to answer it
        // with. pfqn_rgfmc answers it now, and agrees with pfqn_ca to the last
        // bit on this input (lG = 0.8109302162163288 from both), so the name
        // is 'rgf' -- as in pfqn_nc.m, which falls back only on the refusal.
        const mva::AvgResult<double> rgf = run(m, "rgf");
        CHECK(rgf.actualmethod == "rgf");
        CHECK(rgf.QN(0, 0) == doctest::Approx(1.578947368421052).epsilon(1e-12));
        CHECK(rgf.QN(1, 0) == doctest::Approx(1.421052631578947).epsilon(1e-12));
    }
    // 'panald' is whitelisted and REFUSES BY NAME outside normal usage
    // (1 - lambda_i/mu_i(Ntot) <= 0 at a queueing centre), which is what
    // pfqn_ncld.m does on this model too.
    {
        std::string what;
        CHECK_THROWS_AS(
            [&] {
                try {
                    run(m, "panald");
                } catch (const UnsupportedError& e) {
                    what = e.what();
                    throw;
                }
            }(),
            UnsupportedError);
        CHECK(what.find("panald") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// Tier D: the specialised analyzers -- order-independent, pass-and-swap
// importance sampling, class-dependent convolution, and the LCFS closed form
// ---------------------------------------------------------------------------

namespace {

/** Delay + OI queue in a closed loop, one class, a CONSTANT rank rate. */
qn::Network<double> oi_delay_queue(double n = 3.0, double mu = 2.0) {
    qn::Network<double> m("oi1");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service_rate_function(o, [mu](const std::vector<std::size_t>&) { return mu; });
    qn::RoutingMatrix<double> P;
    P.set(d, o, 1.0);
    P.set(o, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the order-independent analyzer reproduces MATLAB on a Delay+OI loop") {
    // MATLAB SolverNC(m).getAvgTable / getProbNormConstAggr; the reported method
    // is 'default/oi', so the exact balanced-fairness route is what ran.
    qn::Network<double> m = oi_delay_queue();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/oi");
    CHECK(r.QN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.42105263157895).epsilon(1e-12));
    CHECK(r.UN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.789473684210526).epsilon(1e-12));
    CHECK(r.RN(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(0.9).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(r.TN(1, 0) == doctest::Approx(1.57894736842105).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(-0.233614851181505).epsilon(1e-12));

    // ANSWER ORACLE, independent of MATLAB: a constant rank rate mu makes the OI
    // station an ordinary M/M/1, so the whole model is the exact closed
    // Delay+queue product form G(N) = sum_k Z^{N-k} (1/mu)^k / (N-k)!, and the
    // throughput is G(N-1)/G(N).
    double G3 = 0.0, G2 = 0.0;
    for (int k = 0; k <= 3; ++k) G3 += std::pow(0.5, k) / std::tgamma(3.0 - k + 1.0);
    for (int k = 0; k <= 2; ++k) G2 += std::pow(0.5, k) / std::tgamma(2.0 - k + 1.0);
    CHECK(r.TN(0, 0) == doctest::Approx(G2 / G3).epsilon(1e-12));
    CHECK(std::log(G3) == doctest::Approx(-0.233614851181505).epsilon(1e-12));
}

TEST_CASE("the OI analyzer carries a count-dependent rank rate and two classes") {
    qn::Network<double> m("oi2");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t a = m.add_closed_class("C1", 2.0, d);
    const std::size_t b = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, a, D::exp_rate(1.0));
    m.set_service(d, b, D::exp_rate(2.0));
    // the mean over the queue of the per-position class rate (1+r): permutation
    // invariant, hence order independent, but NOT constant in the counts
    m.set_service_rate_function(o, [](const std::vector<std::size_t>& cc) {
        double s = 0.0;
        for (std::size_t x : cc) s += 1.0 + static_cast<double>(x);
        return s / static_cast<double>(cc.size());
    });
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {a, b}) {
        P.set(c, c, d, o, 1.0);
        P.set(c, c, o, d, 1.0);
    }
    m.link(P);

    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/oi");
    CHECK(r.QN(0, 0) == doctest::Approx(1.03030303030303).epsilon(1e-12));
    CHECK(r.QN(0, 1) == doctest::Approx(0.454545454545455).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(0.96969696969697).epsilon(1e-12));
    CHECK(r.QN(1, 1) == doctest::Approx(0.545454545454545).epsilon(1e-12));
    // the OI row is the IN-SERVICE utilization, which is NOT T/mu(e_r) here
    CHECK(r.UN(1, 0) == doctest::Approx(0.558441558441558).epsilon(1e-12));
    CHECK(r.UN(1, 1) == doctest::Approx(0.545454545454545).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(1.03030303030303).epsilon(1e-12));
    CHECK(r.TN(0, 1) == doctest::Approx(0.909090909090909).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(0.318453731118535).epsilon(1e-12));
    // conservation, which no tolerance can hide
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.QN(0, 1) + r.QN(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("the OI analyzer folds an ordinary BCMP queue in by lattice convolution") {
    qn::Network<double> m("oi3");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t o = m.add_queue("OIQ", SchedStrategy::OI);
    const std::size_t p = m.add_queue("PSQ", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(p, c, D::exp_rate(2.5));
    m.set_service_rate_function(o, [](const std::vector<std::size_t>&) { return 2.0; });
    qn::RoutingMatrix<double> P;
    P.set(d, o, 1.0);
    P.set(o, p, 1.0);
    P.set(p, d, 1.0);
    m.link(P);

    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/oi");
    CHECK(r.QN(0, 0) == doctest::Approx(1.25966158345519).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.00584917484855).epsilon(1e-12));
    CHECK(r.QN(2, 0) == doctest::Approx(0.734489241696261).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.629830791727595).epsilon(1e-12));
    CHECK(r.UN(2, 0) == doctest::Approx(0.503864633382076).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(1.25966158345519).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(0.467291621742262).epsilon(1e-12));
}

TEST_CASE("the OI and P&S predicates separate the two analyzers") {
    qn::Network<double> oi = oi_delay_queue();
    CHECK(nc::nc_is_oi_model(oi.get_struct()));
    // two stations, but one is a Delay, and the P&S predicate needs BOTH to be
    // OI/PAS with a rate function
    CHECK_FALSE(nc::nc_is_pas_model(oi.get_struct()));
    CHECK_FALSE(nc::nc_is_oi_model(closed_single_class().get_struct()));
    CHECK_FALSE(nc::nc_is_oi_model(open_tandem().get_struct()));

    // A PAS station built with NO explicit swap graph takes MATLAB's default,
    // the complete compatibility graph ones(R,R)-eye(R), so it is a GENUINE
    // pass-and-swap station and NOT order independent. Reading the raw empty
    // matrix instead would have made it OI and routed it to the exact analyzer.
    qn::Network<double> pas("pasdef");
    const std::size_t p1 = pas.add_queue("PASQueue1", SchedStrategy::PAS);
    const std::size_t p2 = pas.add_queue("PASQueue2", SchedStrategy::PAS);
    const std::size_t c1 = pas.add_closed_class("C1", 1.0, p1);
    const std::size_t c2 = pas.add_closed_class("C2", 1.0, p1);
    pas.set_service_rate_function(p1, [](const std::vector<std::size_t>&) { return 1.0; });
    pas.set_service_rate_function(p2, [](const std::vector<std::size_t>&) { return 1.3; });
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {
        P.set(c, c, p1, p2, 1.0);
        P.set(c, c, p2, p1, 1.0);
    }
    pas.link(P);
    CHECK_FALSE(nc::nc_is_oi_model(pas.get_struct()));
    CHECK(nc::nc_is_pas_model(pas.get_struct()));
}

namespace {

/** The Fig. 6 closed P&S tandem of Comte and Dorsman, trimmed to three classes. */
qn::Network<double> pas_tandem() {
    qn::Network<double> m("pas4");
    const std::size_t q1 = m.add_queue("PASQueue1", SchedStrategy::PAS);
    const std::size_t q2 = m.add_queue("PASQueue2", SchedStrategy::PAS);
    std::vector<std::size_t> jc;
    for (int r = 1; r <= 3; ++r)
        jc.push_back(m.add_closed_class("Class" + std::to_string(r), 1.0, q1));
    Matrix<double> G(3, 3, 0.0);
    G(0, 1) = 1.0;
    G(1, 0) = 1.0;
    m.set_service_rate_function(q1, [](const std::vector<std::size_t>&) { return 1.0; }, G);
    m.set_service_rate_function(q2, [](const std::vector<std::size_t>&) { return 1.3; }, G);
    m.set_number_of_servers(q1, 1.0);
    m.set_number_of_servers(q2, 1.0);
    qn::RoutingMatrix<double> P;
    for (std::size_t c : jc) {
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, q1, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the pass-and-swap importance sampler agrees with MATLAB up to its own noise") {
    // THIS IS A MONTE CARLO ESTIMATOR AND THE STREAMS DIFFER. MATLAB draws from
    // its own `rand`, this port from mt19937_64, so a bit-for-bit match is not
    // available at any sample count and pinning one would be pinning the RNG.
    // What IS comparable is the estimate: the values below are MATLAB's
    // SolverNC(m,'method','is','samples',2e4,'seed',23456).getAvgTable, and the
    // tolerance is the measured agreement, reported by the MESSAGE below.
    qn::Network<double> m = pas_tandem();
    nc::NcSolverOptions opt = options("is");
    opt.samples = 20000;
    opt.seed = 23456;
    const mva::AvgResult<double> r = nc::solver_nc_run_analyzer(m.get_struct(), opt);
    CHECK(r.actualmethod == "is");

    const double mQ[2][3] = {{0.768350123568727, 0.447382673814972, 0.605990170128877},
                             {0.231649876431273, 0.552617326185028, 0.394009829871123}};
    const double mT[2][3] = {{0.558876593184683, 0.558876593184683, 0.279438296592341},
                             {0.558876593184683, 0.558876593184683, 0.279438296592341}};
    double worstQ = 0.0, worstT = 0.0;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            worstQ = std::max(worstQ, std::fabs(r.QN(i, k) - mQ[i][k]));
            worstT = std::max(worstT, std::fabs(r.TN(i, k) - mT[i][k]));
        }
    INFO("pas_is vs MATLAB at 2e4 samples: max|dQ| = ", worstQ, ", max|dT| = ", worstT);
    CHECK(worstQ < 2e-2);
    CHECK(worstT < 2e-2);

    // PROPERTIES THE ESTIMATOR MUST SATISFY EXACTLY, whatever the stream. The
    // per-class queue lengths are complementary by construction inside
    // pfqn_pas_is, and the closed tandem carries the same throughput at both
    // stations because the visits are unit.
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.QN(0, k) + r.QN(1, k) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(r.TN(0, k) == doctest::Approx(r.TN(1, k)).epsilon(1e-12));
        // utilization at station 2 is T/1.3, the offered-load form the analyzer
        // uses for a P&S station
        CHECK(r.UN(1, k) == doctest::Approx(r.TN(1, k) / 1.3).epsilon(1e-12));
    }
}

namespace {

/** Delay + PS queue whose class-1 count scales the rate, as `ld_class_dependence`. */
qn::Network<double> cd_model(double n1 = 4.0, double n2 = 2.0, double cap = 2.0) {
    qn::Network<double> m("cd5");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t a = m.add_closed_class("Class1", n1, d);
    const std::size_t b = m.add_closed_class("Class2", n2, d);
    m.set_service(d, a, D::exp_rate(1.0));
    m.set_service(d, b, D::exp_rate(0.5));
    m.set_service(q, a, D::exp_rate(1.0 / 1.5));
    m.set_service(q, b, D::exp_rate(1.0 / 2.5));
    m.set_class_dependence(
        q,
        [cap](const std::vector<double>& ni) {
            return std::vector<double>{std::min(ni[0], cap)};
        },
        std::vector<double>{cap});
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {a, b}) {
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the convolution analyzer reproduces MATLAB on a class-dependent station") {
    qn::Network<double> m = cd_model();
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/conv");
    CHECK(r.QN(0, 0) == doctest::Approx(1.23798882681564).epsilon(1e-12));
    CHECK(r.QN(0, 1) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(2.76201117318436).epsilon(1e-12));
    CHECK(r.UN(0, 0) == doctest::Approx(1.23798882681564).epsilon(1e-12));
    CHECK(r.UN(0, 1) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.928491620111732).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(2.23104693140794).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(1.23798882681564).epsilon(1e-12));
    CHECK(r.TN(0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(1.53932834624707).epsilon(1e-12));

    // WHY CLASS 2 HOLDS NOTHING AT THE QUEUE, and it is not a port defect. The
    // reference builds the station factor X_m(n) by picking the FIRST class with
    // n_r > 0 and asserting the recurrence is path independent. It is not when
    // beta can vanish: beta(n) = min(n_1, c) is zero at every n with n_1 = 0, so
    // X_m(0, k) = 0, and every state with class-2 jobs is reached through one of
    // those. Reproduced rather than corrected -- MATLAB is ground truth and the
    // conservation sums below still hold exactly.
    CHECK(r.QN(1, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-12));
    CHECK(r.QN(0, 1) + r.QN(1, 1) == doctest::Approx(2.0).epsilon(1e-12));

    // The utilization normalization is what makes a beta emulating c servers
    // report E[busy]/c: without the division by the declared peak, U at the
    // queue would be T*S = 2.785, well above one.
    CHECK(r.TN(1, 0) * 1.5 == doctest::Approx(r.UN(1, 0) * 2.0).epsilon(1e-12));
    CHECK(r.UN(1, 0) < 1.0);
}

TEST_CASE("a class-dependent station with no declared peak is refused, by name") {
    // MATLAB's getLimitedClassDependencePeak errors on exactly this, since the
    // peak cannot be recovered from beta without sweeping the whole lattice.
    qn::Network<double> m("cdnopeak");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t a = m.add_closed_class("Class1", 2.0, d);
    m.set_service(d, a, D::exp_rate(1.0));
    m.set_service(q, a, D::exp_rate(1.0));
    m.set_class_dependence(q, [](const std::vector<double>& ni) {
        return std::vector<double>{std::min(ni[0], 2.0)};
    });
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    CHECK_THROWS_AS(run(m), UnsupportedError);
}

namespace {

/** LCFS + LCFS-PR closed tandem; `mu(i,r)` is the service RATE at station i. */
qn::Network<double> lcfs_tandem(const std::vector<std::vector<double>>& mu,
                                const std::vector<double>& pop) {
    qn::Network<double> m("lcfs");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::LCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::LCFSPR);
    for (std::size_t r = 0; r < pop.size(); ++r) {
        const std::size_t c = m.add_closed_class("Class" + std::to_string(r + 1), pop[r], q1);
        m.set_service(q1, c, D::exp_rate(mu[0][r]));
        m.set_service(q2, c, D::exp_rate(mu[1][r]));
    }
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 1; r <= pop.size(); ++r) {
        P.set(r, r, q1, q2, 1.0);
        P.set(r, r, q2, q1, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the LCFS closed form reproduces MATLAB with one job per class") {
    qn::Network<double> m =
        lcfs_tandem({{1.0, 1.0 / 3.0, 0.2}, {0.5, 0.25, 1.0 / 6.0}}, {1.0, 1.0, 1.0});
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/lcfsqn.ca");
    CHECK(r.QN(0, 0) == doctest::Approx(0.460992907801418).epsilon(1e-12));
    CHECK(r.QN(0, 1) == doctest::Approx(0.312394461330632).epsilon(1e-12));
    CHECK(r.QN(0, 2) == doctest::Approx(0.238095238095238).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(0.539007092198582).epsilon(1e-12));
    CHECK(r.QN(1, 1) == doctest::Approx(0.687605538669368).epsilon(1e-12));
    CHECK(r.QN(1, 2) == doctest::Approx(0.761904761904762).epsilon(1e-12));
    CHECK(r.UN(0, 0) == doctest::Approx(0.204322863897332).epsilon(1e-12));
    CHECK(r.UN(1, 0) == doctest::Approx(0.408645727794664).epsilon(1e-12));
    CHECK(r.RN(0, 0) == doctest::Approx(2.25619834710744).epsilon(1e-12));
    CHECK(r.RN(1, 2) == doctest::Approx(22.3366336633663).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(0.204322863897332).epsilon(1e-12));
    CHECK(r.TN(0, 2) == doctest::Approx(0.0341100979398852).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(9.78504179732965).epsilon(1e-12));
}

TEST_CASE("the LCFS closed form expands a class of multiplicity two") {
    // The permanent formulas assume one job per class; a class with N_r > 1 is
    // expanded into N_r exchangeable copies against G_exp = G prod_r N_r!. This
    // is the case that distinguishes a correct expansion from an omitted one.
    qn::Network<double> m = lcfs_tandem({{1.0, 1.0 / 3.0}, {0.5, 0.25}}, {2.0, 1.0});
    const mva::AvgResult<double> r = run(m);
    CHECK(r.actualmethod == "default/lcfsqn.ca");
    CHECK(r.QN(0, 0) == doctest::Approx(0.650602409638554).epsilon(1e-12));
    CHECK(r.QN(0, 1) == doctest::Approx(0.180722891566265).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(1.34939759036145).epsilon(1e-12));
    CHECK(r.QN(1, 1) == doctest::Approx(0.819277108433735).epsilon(1e-12));
    CHECK(r.UN(0, 0) == doctest::Approx(0.36144578313253).epsilon(1e-12));
    CHECK(r.UN(1, 1) == doctest::Approx(0.183132530120482).epsilon(1e-12));
    CHECK(r.TN(0, 0) == doctest::Approx(0.36144578313253).epsilon(1e-12));
    CHECK(r.TN(0, 1) == doctest::Approx(0.0457831325301205).epsilon(1e-12));
    CHECK(nc::solver_nc_lognormconst(m.get_struct(), options("default")) ==
          doctest::Approx(6.0282785202307).epsilon(1e-12));
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.QN(0, 1) + r.QN(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("an unpaired LCFS station is refused, by name") {
    qn::Network<double> m("lcfsalone");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", SchedStrategy::LCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    // A non-preemptive LCFS station is not a BCMP station; solving it as an
    // ordinary queue would return the wrong network's numbers.
    CHECK_THROWS_AS(run(m), UnsupportedError);
}

namespace {

/** Delay -> Cache -> Fetch -> Cache, CLOSED, with a delayed-hit retrieval system. */
qn::Network<double> closed_retrieval_model() {
    qn::Network<double> m("ClosedDelayedHits");
    const std::size_t d = m.add_delay("Delay");
    qn::CacheParam<double> ch;
    ch.nitems = 3;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double>>{{0.6, 0.3, 0.1}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t q1 = m.add_queue("Fetch", SchedStrategy::PS);
    const std::size_t job = m.add_closed_class("InitClass", 2.0, d);
    const std::size_t hit = m.add_closed_class("HitClass", 0.0, d);
    const std::size_t mis = m.add_closed_class("MissClass", 0.0, d);
    m.set_service(d, job, D::exp_rate(1.0));
    m.set_service(d, hit, D::exp_rate(1.0));
    m.set_service(d, mis, D::exp_rate(1.0));
    m.set_service(q1, job, D::exp_rate(2.0));
    m.set_retrieval_system(cn, job, mis, {q1});
    qn::RoutingMatrix<double> P;
    P.set(job, job, d, cn, 1.0);
    P.set(job, job, cn, q1, 1.0);
    P.set(job, job, q1, cn, 1.0);
    P.set(hit, job, cn, d, 1.0);
    P.set(mis, job, cn, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the closed delayed-hit retrieval analyzer reproduces MATLAB") {
    // MATLAB's solver_nc_cacheqn_retrieval_analyzer on the same model. Reaching
    // these numbers took two builder fixes: `link()` now consumes the read
    // class's template edges over the retrieval queue set, and `refresh_chains`
    // couples the retrieval classes to the read class's chain. Without them the
    // model carried four singleton chains where MATLAB carries one, and the
    // whole table was a different network's.
    qn::Network<double> m = closed_retrieval_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nchains == 1);
    CHECK(sn.nclasses == 6);

    const nc::NcCacheqnRetrievalSolution<double> r =
        nc::solver_nc_cacheqn_retrieval_analyzer(sn, options("default"));
    CHECK(r.sol.sol.Q(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(r.sol.sol.Tp(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    const double fq[3] = {0.247568598687817, 0.166148415818464, 0.0717543454721912};
    const double ft[3] = {0.445312216396303, 0.298858254605033, 0.12906760706991};
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.sol.sol.Q(1, 3 + k) == doctest::Approx(fq[k]).epsilon(1e-9));
        CHECK(r.sol.sol.Tp(1, 3 + k) == doctest::Approx(ft[k]).epsilon(1e-9));
    }
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.sol.sol.Q(1, k) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(r.sol.sol.Tp(1, k) == doctest::Approx(0.0).epsilon(1e-12));
    }
    CHECK(r.hitprob[0] == doctest::Approx(0.423425840227865).epsilon(1e-9));
    CHECK(r.missprob[0] == doctest::Approx(0.576574159772135).epsilon(1e-9));
    CHECK(r.hitprob[0] + r.missprob[0] == doctest::Approx(1.0).epsilon(1e-12));

    // THE WHOLE RUNNER TABLE, which is what a caller sees. MATLAB's
    // `SolverNC(model).getAvg()` on this model, measured rather than assumed.
    //
    // ONE ANALYZER-LEVEL DIFFERENCE SURVIVES AND IT IS INVISIBLE HERE. MATLAB's
    // analyzer reports Tput 5.41682399813583 for each retrieval class in the
    // DELAY row -- where those classes are disabled -- and this port reports
    // zero. `filterMetric` zeroes it in MATLAB too: a disabled pair keeps its
    // throughput only when the class is a cache HIT or MISS class, and a minted
    // retrieval class is neither. Measured, not reasoned: MATLAB's analyzer TN
    // row is [1.5145 0 0 5.4168 5.4168 5.4168] and its getAvg TN row is
    // [1.5145 0 0 0 0 0].
    const mva::AvgResult<double> a = run(m);
    CHECK(a.TN(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(a.QN(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(a.WN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    for (std::size_t k = 1; k < 6; ++k) {
        CHECK(a.TN(0, k) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(a.QN(0, k) == doctest::Approx(0.0).epsilon(1e-12));
    }
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(a.TN(1, 3 + k) == doctest::Approx(ft[k]).epsilon(1e-9));
        CHECK(a.QN(1, 3 + k) == doctest::Approx(fq[k]).epsilon(1e-9));
        CHECK(a.WN(1, 3 + k) == doctest::Approx(0.83391581088187).epsilon(1e-9));
        CHECK(a.TN(1, k) == doctest::Approx(0.0).epsilon(1e-12));
    }
    // ARRIVAL RATE IS AN OPEN DIVERGENCE AND IS DELIBERATELY NOT PINNED HERE.
    // MATLAB reports an ALL-ZERO AN on this model; this port reports
    // AN(Delay, InitClass) = 1.51452864 and AN(Fetch, retrieval r) equal to the
    // corresponding throughputs. Both are measured, not inferred. The likely
    // cause is `sn_get_arvr_from_tput.m`, which fills the Cache row of
    // TN_stateful only when the node carries `actualhitprob` and otherwise
    // `continue`s with the comment "Arrival rates will be computed later when
    // probabilities are available" -- on this path they never are, so every
    // station fed through the Cache inherits a zero. That reading is NOT yet
    // verified. Pinning either value would settle a question that belongs to
    // the divergence register, so the assertion is withheld rather than guessed.
}

namespace {

/** mu(n) = min(n, 3) on `q`, the only form the load-dependent route admits. */
void three_servers_as_load_dependence(qn::Network<double>& m, std::size_t q) {
    std::vector<double> alpha(40);
    for (std::size_t n = 1; n <= alpha.size(); ++n) alpha[n - 1] = std::min<double>(n, 3.0);
    m.set_load_dependence(q, alpha);
}

/** Source -> 3-server FCFS queue -> Sink: an M/M/3 and nothing closed at all. */
qn::Network<double> open_ld_multiserver() {
    qn::Network<double> m("open_ld");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, D::exp_rate(1.5));
    m.set_service(q, o, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(o, o, s, q, 1.0);
    P.set(o, o, q, k, 1.0);
    m.link(P);
    three_servers_as_load_dependence(m, q);
    return m;
}

/** The same queue shared with a closed chain of ONE job, i.e. fewer than c. */
qn::Network<double> mixed_ld_multiserver() {
    qn::Network<double> m("mixed_ld");
    const std::size_t s = m.add_source("Source");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_closed_class("C1", 1.0, d);
    const std::size_t o = m.add_open_class("O1");
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.0));
    m.set_arrival(s, o, D::exp_rate(0.4));
    m.set_service(q, o, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    P.set(o, o, s, q, 1.0);
    P.set(o, o, q, k, 1.0);
    m.link(P);
    three_servers_as_load_dependence(m, q);
    return m;
}

}  // namespace

TEST_CASE("the mixed load-dependent route reads the whole rate lattice") {
    // `pfqn_ldmx_ec` infers the limited-load-dependence level b_i as the first
    // column of mu equal to the LAST one and treats every rate past it as
    // saturated. The lattice handed to it used to be cut at the closed
    // population, so a c-server station was read as saturated at min(n,c) with
    // n < c whenever c exceeded that population -- and with no closed class at
    // all the row collapsed to mu(1), one server. Both models below have c = 3
    // and at most one closed job, which is exactly that case.
    SUBCASE("a purely open multiserver is the M/M/c queue") {
        qn::Network<double> m = open_ld_multiserver();
        const mva::AvgResult<double> r = run(m, "exact");
        CHECK(r.actualmethod.find("ncldmx") != std::string::npos);
        // exact M/M/3 at lambda = 1.5, mu = 1: L = Lq + a = 0.236842 + 1.5
        CHECK(r.QN(1, 0) == doctest::Approx(1.736842105263158).epsilon(1e-9));
        CHECK(r.UN(1, 0) == doctest::Approx(0.5).epsilon(1e-9));   // busy fraction of 3
        CHECK(r.TN(1, 0) == doctest::Approx(1.5).epsilon(1e-9));
    }
    SUBCASE("c above the closed population is still c") {
        // SolverCTMC (cutoff 18) on the same model, agreeing to 1e-12
        qn::Network<double> m = mixed_ld_multiserver();
        const mva::AvgResult<double> r = run(m, "exact");
        CHECK(r.actualmethod.find("ncldmx") != std::string::npos);
        CHECK(r.QN(2, 0) == doctest::Approx(0.500791765637).epsilon(1e-9));
        CHECK(r.QN(2, 1) == doctest::Approx(0.405871246726).epsilon(1e-9));
        CHECK(r.QN(1, 0) == doctest::Approx(0.499208234363).epsilon(1e-9));
        CHECK(r.UN(2, 0) == doctest::Approx(0.166402744788).epsilon(1e-9));
        CHECK(r.UN(2, 1) == doctest::Approx(0.133333333333).epsilon(1e-9));
        CHECK(r.TN(2, 0) == doctest::Approx(0.499208234363).epsilon(1e-9));
        CHECK(r.TN(2, 1) == doctest::Approx(0.4).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_rgfmc matches the convolution with think times") {
    // Multiclass RGF: iterated residues (Harrison-Coury 2002 Thm 1) with the delay
    // carried by the Bertozzi-McKenna truncation, which neither RGF paper has.
    using namespace line;
    using namespace line::pfqn;
    const double L2d[3][2] = {{0.9, 0.4}, {0.6, 0.7}, {0.3, 1.1}};
    Matrix<double> L2(3, 2);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j) L2(i, j) = L2d[i][j];
    const std::vector<int> N2{3, 4};
    const double Z2d[3][2] = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.5}};
    for (int t = 0; t < 3; ++t) {
        const std::vector<double> Zv{Z2d[t][0], Z2d[t][1]};
        Matrix<double> Zm(1, 2);
        Zm(0, 0) = Z2d[t][0];
        Zm(0, 1) = Z2d[t][1];
        CHECK(pfqn_rgfmc<double>(L2, N2, Zv).lG ==
              doctest::Approx(pfqn_ca<double>(L2, N2, Zm).lG).epsilon(1e-10));
    }
    // A class removed by residues enters only as a pole ORDER, so its population is
    // nearly free while it carries no think time (Harrison-Lee sec. 4).
    const double L4d[4][2] = {{0.9, 0.4}, {0.6, 0.7}, {0.3, 1.1}, {0.5, 0.2}};
    Matrix<double> L4(4, 2);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 2; ++j) L4(i, j) = L4d[i][j];
    const std::vector<double> Zfree{1.5, 0.0};
    Matrix<double> Zfm(1, 2);
    Zfm(0, 0) = 1.5;
    Zfm(0, 1) = 0.0;
    const int n2[3] = {20, 200, 2000};
    for (int a = 0; a < 3; ++a) {
        const std::vector<int> N{6, n2[a]};
        CHECK(pfqn_rgfmc<double>(L4, N, Zfree).lG ==
              doctest::Approx(pfqn_ca<double>(L4, N, Zfm).lG).epsilon(1e-10));
    }
}

TEST_CASE("pfqn_rgfmc refuses rather than returning a wrong constant") {
    // Near-coincident loads over the eliminated class make the alternating residue
    // sum cancel; the guard must refuse, not answer.
    using namespace line;
    using namespace line::pfqn;
    const double Ld[2][3] = {{0.81733524, 1.24343101, 0.86870538},
                             {1.4732791, 0.38631325, 0.87522251}};
    Matrix<double> L(2, 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j) L(i, j) = Ld[i][j];
    const std::vector<int> N{2, 3, 2};
    const std::vector<double> Z{0.0, 0.0, 0.0};
    CHECK_THROWS_AS(pfqn_rgfmc<double>(L, N, Z, 1e-12, std::size_t(1000000), 1.0),
                    InputError);
}
