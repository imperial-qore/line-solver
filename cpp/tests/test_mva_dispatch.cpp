/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The mvaDispatch ladder: which analyzer a model gets, and what it returns.
 *
 * Two things are asserted here and they are different in kind. The NUMBERS come
 * from MATLAB (`SolverMVA(model).getAvg()`) and pin the ported analyzers to the
 * reference. The `actualmethod` strings pin the ROUTING -- the branches are not
 * disjoint, so a model reaching the right numbers through the wrong branch is a
 * latent defect that only shows up on the next model.
 *
 * The refusals matter as much: a cache or a polling model must FAIL rather than
 * fall through to the generic analyzer, which would return a plausible answer
 * for a different model.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/mva_dispatch.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

mva::DispatchResult<double> go(qn::Network<double>& m, const std::string& method = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::mva_dispatch(m.get_struct(), opt, init);
}

/** Source -> Queue -> Sink with one open class, the shape of every qsys case. */
qn::Network<double> sqs(const std::string& name, SchedStrategy sched, const D& arrival,
                        const D& service, double servers = 1.0) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(s, o, arrival);
    m.set_service(q, o, service);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("a single-class open Source-Queue-Sink takes the queueing-system closed forms") {
    SUBCASE("M/M/1") {
        qn::Network<double> m = sqs("mm1", SchedStrategy::FCFS, D::exp_rate(1.0), D::exp_rate(2.0));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "mm1");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }
    SUBCASE("M/M/2") {
        qn::Network<double> m =
            sqs("mmk", SchedStrategy::FCFS, D::exp_rate(1.5), D::exp_rate(1.0), 2.0);
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "mmk");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(3.428571429).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(2.285714286).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.75).epsilon(1e-9));
    }
    SUBCASE("M/G/1, Erlang service") {
        qn::Network<double> m =
            sqs("mg1", SchedStrategy::FCFS, D::exp_rate(1.0), D::erlang_fit(0.5, 0.25));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "mg1");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(0.8125).epsilon(1e-9));
    }
    SUBCASE("G/M/1, Erlang arrivals: the PH sigma-root, exact for a Markovian arrival") {
        qn::Network<double> m =
            sqs("gm1", SchedStrategy::FCFS, D::erlang_fit(1.0, 0.5), D::exp_rate(2.0));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "gm1");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(0.809016994).epsilon(1e-9));
    }
    SUBCASE("G/G/1 falls back to the KLB approximation") {
        qn::Network<double> m =
            sqs("gig1", SchedStrategy::FCFS, D::erlang_fit(1.0, 0.5), D::erlang_fit(0.5, 0.25));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "gig1.klb");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(0.650138263).epsilon(1e-9));
    }
    SUBCASE("Trace/M/1 takes the EMPIRICAL transform, not the APH fit of the trace") {
        // MATLAB `Replayer.evalLST` (and the JAR/python twins) return
        // mean(exp(-s x)) over the trace, so the G/M/1 caudal root is the root
        // of the empirical transform. Reading it off the APH fit of the first
        // three moments solves a different arrival law and moves QLen by ~1e-4
        // relative on gallery_replayerm1.
        std::vector<double> trace;
        trace.push_back(1.0);
        trace.push_back(3.0);
        const D arv = D::replayer(trace);
        CHECK(lang::dist_lst(arv, 0.5) ==
              doctest::Approx(0.5 * (std::exp(-0.5) + std::exp(-1.5))).epsilon(1e-12));

        const double mu = 1.0;
        qn::Network<double> m = sqs("gm1trace", SchedStrategy::FCFS, arv, D::exp_rate(mu));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "gm1");
        // sigma = 0.5 (exp(-mu(1-sigma)) + exp(-3 mu(1-sigma))), by bisection
        double lo = 1e-12, hi = 1.0 - 1e-12;
        for (int it = 0; it < 200; ++it) {
            const double x = 0.5 * (lo + hi);
            const double s = mu * (1.0 - x);
            const double f = 0.5 * (std::exp(-s) + std::exp(-3.0 * s)) - x;
            if (f > 0.0) lo = x; else hi = x;
        }
        const double sigma = 0.5 * (lo + hi);
        CHECK(r.sol.R(1, 0) == doctest::Approx(1.0 / (mu * (1.0 - sigma))).epsilon(1e-9));
    }
    SUBCASE("exact refuses what no closed form covers") {
        qn::Network<double> m =
            sqs("gig1x", SchedStrategy::FCFS, D::erlang_fit(1.0, 0.5), D::erlang_fit(0.5, 0.25));
        CHECK_THROWS_AS(go(m, "exact"), UnsupportedError);
    }
}

TEST_CASE("a multiclass open HOL queue takes the exact Cobham formula") {
    qn::Network<double> m("hol");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::HOL);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t hi = m.add_open_class("Hi", 0);
    const std::size_t lo = m.add_open_class("Lo", 1);
    m.set_arrival(s, hi, D::exp_rate(0.4));
    m.set_arrival(s, lo, D::exp_rate(0.4));
    m.set_service(q, hi, D::exp_rate(2.0));
    m.set_service(q, lo, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(hi, hi, s, q, 1.0);
    P.set(hi, hi, q, k, 1.0);
    P.set(lo, lo, s, q, 1.0);
    P.set(lo, lo, q, k, 1.0);
    m.link(P);
    const mva::DispatchResult<double> r = go(m);
    CHECK(r.actualmethod == "mg1.prio");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(0.3).epsilon(1e-9));
    CHECK(r.sol.Q(1, 1) == doctest::Approx(0.366666667).epsilon(1e-9));
    CHECK(r.sol.R(1, 0) == doctest::Approx(0.75).epsilon(1e-9));
    CHECK(r.sol.R(1, 1) == doctest::Approx(0.916666667).epsilon(1e-9));
}

TEST_CASE("a multiclass open DPS queue takes the exact truncated chain") {
    qn::Network<double> m("dpsopen");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::DPS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t a = m.add_open_class("A");
    const std::size_t b = m.add_open_class("B");
    m.set_arrival(s, a, D::exp_rate(0.4));
    m.set_arrival(s, b, D::exp_rate(0.4));
    m.set_service(q, a, D::exp_rate(2.0));
    m.set_service(q, b, D::exp_rate(2.0));
    m.set_sched_param(q, a, 1.0);
    m.set_sched_param(q, b, 3.0);
    qn::RoutingMatrix<double> P;
    P.set(a, a, s, q, 1.0);
    P.set(a, a, q, k, 1.0);
    P.set(b, b, s, q, 1.0);
    P.set(b, b, q, k, 1.0);
    m.link(P);
    const mva::DispatchResult<double> r = go(m);
    CHECK(r.actualmethod == "mm1.dps");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(0.375).epsilon(1e-9));
    CHECK(r.sol.Q(1, 1) == doctest::Approx(0.291666667).epsilon(1e-9));
    CHECK(r.sol.R(1, 0) == doctest::Approx(0.9375).epsilon(1e-9));
    CHECK(r.sol.R(1, 1) == doctest::Approx(0.729166667).epsilon(1e-9));
}

TEST_CASE("a load-dependent closed model takes the exact load-dependent recursion") {
    qn::Network<double> m("ld");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 4.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_load_dependence(q, std::vector<double>{1.0, 1.5, 1.8, 2.0});
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const mva::DispatchResult<double> r = go(m);
    CHECK(r.actualmethod == "exact");
    CHECK(r.sol.Q(0, 0) == doctest::Approx(2.365217391).epsilon(1e-9));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(1.634782609).epsilon(1e-9));
    // the NC convention: the carried load over max(nservers, max lldscaling)
    CHECK(r.sol.U(1, 0) == doctest::Approx(0.591304348).epsilon(1e-9));
    CHECK(r.sol.R(1, 0) == doctest::Approx(0.691176471).epsilon(1e-9));
}

TEST_CASE("an ordinary closed network still reaches the generic analyzer") {
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
    const mva::DispatchResult<double> r = go(m);
    CHECK(r.actualmethod == "exact");
    CHECK(r.sol.Q(1, 0) == doctest::Approx(1.421052632).epsilon(1e-9));
}

TEST_CASE("the branches whose analyzers are not ported refuse by name") {
    SUBCASE("an order-independent station") {
        qn::Network<double> m("oi");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::OI);
        const std::size_t c = m.add_closed_class("C1", 2.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(go(m), UnsupportedError);
    }
    SUBCASE("an open polling system now reaches the polling analyzer") {
        // polling was refused before its analyzer landed; a POLLING queue with
        // a discipline set now dispatches to solver_mva_polling_analyzer, which
        // the dedicated test_mva_polling.cpp checks numerically
        qn::Network<double> m("poll");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::POLLING);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t a = m.add_open_class("A");
        const std::size_t b = m.add_open_class("B");
        m.set_arrival(s, a, D::exp_rate(0.2));
        m.set_arrival(s, b, D::exp_rate(0.2));
        m.set_service(q, a, D::exp_rate(2.0));
        m.set_service(q, b, D::exp_rate(2.0));
        m.set_switchover(q, a, D::exp_rate(1e8));
        m.set_switchover(q, b, D::exp_rate(1e8));
        m.set_polling_type(q, lang::PollingType::EXHAUSTIVE);
        qn::RoutingMatrix<double> P;
        P.set(a, a, s, q, 1.0);
        P.set(a, a, q, k, 1.0);
        P.set(b, b, s, q, 1.0);
        P.set(b, b, q, k, 1.0);
        m.link(P);
        CHECK_NOTHROW(go(m));
    }
    SUBCASE("the bound family, which moved to SolverBA") {
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
        CHECK_THROWS_AS(go(m, "aba.upper"), UnsupportedError);
    }
}

TEST_CASE("the size-based disciplines take the M/G/1 formulas") {
    // Two open classes, arrival 0.3 each, service rates 2 and 1. The reference
    // values come from solver_mva_qsys_sizebased_analyzer called directly:
    // SolverMVA's own featset does not list these disciplines, so the branch is
    // unreachable through the MATLAB solver even though the analyzer exists.
    // 2026-07-29: mva_feature_set no longer lists them either, so the branch is
    // now unreachable through the C++ gate as well. This case still exercises
    // it because go() calls mva_dispatch directly and never passes the gate,
    // which is only reached from solver_mva_runner.h.
    struct Case {
        SchedStrategy sched;
        const char* actual;
        double Qa, Qb, Ra, Rb;
    };
    const Case cases[] = {
        {SchedStrategy::SRPT, "mg1.srpt", 0.186877828, 0.422319447, 0.622926092, 1.407731492},
        {SchedStrategy::PSJF, "mg1.psjf", 0.190565371, 0.456842793, 0.635217903, 1.522809311},
        {SchedStrategy::FB, "mg1.fb", 0.245640547, 0.558997878, 0.818801822, 1.863326259},
        {SchedStrategy::LRPT, "mg1.lrpt", 0.644628087, 0.917355348, 2.148760289, 3.057851158},
        {SchedStrategy::SETF, "mg1.setf", 0.499106083, 0.851483112, 1.663686944, 2.838277041},
    };
    for (const Case& c : cases) {
        qn::Network<double> m("sb");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", c.sched);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t a = m.add_open_class("A");
        const std::size_t b = m.add_open_class("B");
        m.set_arrival(s, a, D::exp_rate(0.3));
        m.set_arrival(s, b, D::exp_rate(0.3));
        m.set_service(q, a, D::exp_rate(2.0));
        m.set_service(q, b, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(a, a, s, q, 1.0);
        P.set(a, a, q, k, 1.0);
        P.set(b, b, s, q, 1.0);
        P.set(b, b, q, k, 1.0);
        m.link(P);
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == c.actual);
        CHECK(r.sol.Q(1, 0) == doctest::Approx(c.Qa).epsilon(1e-9));
        CHECK(r.sol.Q(1, 1) == doctest::Approx(c.Qb).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(c.Ra).epsilon(1e-9));
        CHECK(r.sol.R(1, 1) == doctest::Approx(c.Rb).epsilon(1e-9));
    }
}

TEST_CASE("Marie's aggregation-decomposition solves a closed non-exponential FCFS model") {
    qn::Network<double> m("marie");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 6.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::erlang_fit(0.4, 0.25));  // FCFS is service-sensitive
    m.set_service(q2, c, D::exp_rate(3.0));          // PS is insensitive: SCV forced to 1
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    const mva::DispatchResult<double> r = go(m, "marie");
    CHECK(r.actualmethod == "marie");
    CHECK(r.sol.Q(0, 0) == doctest::Approx(2.179188329).epsilon(1e-9));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(2.198679953).epsilon(1e-9));
    CHECK(r.sol.Q(2, 0) == doctest::Approx(1.622131718).epsilon(1e-9));
    CHECK(r.sol.U(1, 0) == doctest::Approx(0.871675332).epsilon(1e-9));
    CHECK(r.sol.R(1, 0) == doctest::Approx(1.008944442).epsilon(1e-9));
}

TEST_CASE("Marie refuses the model shapes its decomposition does not cover") {
    SUBCASE("an open model") {
        qn::Network<double> m =
            sqs("marieopen", SchedStrategy::FCFS, D::exp_rate(1.0), D::exp_rate(2.0));
        CHECK_THROWS_AS(go(m, "marie"), UnsupportedError);
    }
    SUBCASE("a discipline outside FCFS / PS / LCFSPR / Delay") {
        qn::Network<double> m("mariedps");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::DPS);
        const std::size_t c = m.add_closed_class("C1", 2.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        m.set_sched_param(q, c, 1.0);
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(go(m, "marie"), UnsupportedError);
    }
}

TEST_CASE("the summation method solves closed and open models") {
    SUBCASE("closed, two classes, FCFS with Erlang service and PS") {
        qn::Network<double> m("sum");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
        const std::size_t c1 = m.add_closed_class("C1", 4.0, d);
        const std::size_t c2 = m.add_closed_class("C2", 3.0, d);
        m.set_service(d, c1, D::exp_rate(1.0));
        m.set_service(d, c2, D::exp_rate(2.0));
        m.set_service(q1, c1, D::erlang_fit(0.4, 0.25));
        m.set_service(q1, c2, D::exp_rate(2.5));
        m.set_service(q2, c1, D::exp_rate(3.0));
        m.set_service(q2, c2, D::exp_rate(4.0));
        qn::RoutingMatrix<double> P;
        P.set(c1, c1, d, q1, 1.0);
        P.set(c1, c1, q1, q2, 1.0);
        P.set(c1, c1, q2, d, 1.0);
        P.set(c2, c2, d, q1, 1.0);
        P.set(c2, c2, q1, q2, 1.0);
        P.set(c2, c2, q2, d, 1.0);
        m.link(P);
        // SUM and ESUM coincide here: only FCFS and SIRO are service-sensitive,
        // and the correction is applied by both under the same name.
        for (const char* meth : {"sum", "esum"}) {
            const mva::DispatchResult<double> r = go(m, meth);
            CHECK(r.sol.Q(0, 0) == doctest::Approx(1.157322539).epsilon(1e-9));
            CHECK(r.sol.Q(0, 1) == doctest::Approx(0.542395739).epsilon(1e-9));
            CHECK(r.sol.Q(1, 0) == doctest::Approx(1.959659974).epsilon(1e-9));
            CHECK(r.sol.Q(1, 1) == doctest::Approx(1.836845276).epsilon(1e-9));
            CHECK(r.sol.Q(2, 0) == doctest::Approx(0.883018696).epsilon(1e-9));
            CHECK(r.sol.Q(2, 1) == doctest::Approx(0.620758987).epsilon(1e-9));
            // the closed populations, to the method's own tolerance
            double n1 = 0.0, n2 = 0.0;
            for (std::size_t i = 0; i < r.sol.Q.rows(); ++i) {
                n1 += r.sol.Q(i, 0);
                n2 += r.sol.Q(i, 1);
            }
            CHECK(n1 == doctest::Approx(4.0).epsilon(1e-6));
            CHECK(n2 == doctest::Approx(3.0).epsilon(1e-6));
        }
    }
    SUBCASE("open, one class through two queues: the closing method") {
        qn::Network<double> m("sumopen");
        const std::size_t s = m.add_source("Source");
        const std::size_t qa = m.add_queue("Qa", SchedStrategy::FCFS);
        const std::size_t qb = m.add_queue("Qb", SchedStrategy::PS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t o = m.add_open_class("O1");
        m.set_arrival(s, o, D::exp_rate(0.5));
        m.set_service(qa, o, D::exp_rate(2.0));
        m.set_service(qb, o, D::exp_rate(1.5));
        qn::RoutingMatrix<double> P;
        P.set(s, qa, 1.0);
        P.set(qa, qb, 1.0);
        P.set(qb, k, 1.0);
        m.link(P);
        const mva::DispatchResult<double> r = go(m, "sum");
        CHECK(r.sol.Q(1, 0) == doctest::Approx(0.333310901).epsilon(1e-9));
        CHECK(r.sol.Q(2, 0) == doctest::Approx(0.499949647).epsilon(1e-9));
        // Both stations are exponential single servers, so the product form
        // gives rho/(1-rho) = 1/3 and 1/2; the gap is the Kclosed = 5000
        // truncation of the closing method, not an error in the port.
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-3));
        CHECK(r.sol.Q(2, 0) == doctest::Approx(0.5).epsilon(1e-3));
    }
}

// ---------------------------------------------------------------------------
// Branch 0: shortest-job-next
//
// Every number below is MATLAB's, from `SolverMVA(model).getAvg()` on the same
// model, and the `actualmethod` strings pin the lattice/fixed-point choice that
// mvaDispatch.m delegates to solver_mva_sjn_analyzer.m. SJN is NOT product form
// and the branch sits ABOVE every other one, so a model that reached the
// generic AMVA path instead would still return plausible numbers -- these cases
// are what separates the two.
// ---------------------------------------------------------------------------

/** Delay -> SJF queue -> Delay, one closed class, the shape of every SJN case. */
qn::Network<double> sjn_closed(const std::string& name, double njobs, const D& think,
                               const D& service, int prio = 0) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
    const std::size_t c = m.add_closed_class("C1", njobs, d, prio);
    m.set_service(d, c, think);
    m.set_service(q, c, service);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("a closed model with a shortest-job-next station takes the Kant recursion") {
    SUBCASE("exponential service, the lattice recursion") {
        // N = 3 gives a lattice of 4 states, far below sjn_lattice_max = 1e5,
        // so `default` runs the exact population recursion.
        qn::Network<double> m = sjn_closed("sjnA", 3.0, D::exp_rate(1.0), D::exp_rate(2.0));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "sjn.mva");
        CHECK(r.sol.Q(0, 0) == doctest::Approx(1.57374156251).epsilon(1e-9));
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.42625843749).epsilon(1e-9));
        CHECK(r.sol.U(0, 0) == doctest::Approx(1.57374156251).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.786870781255).epsilon(1e-9));
        CHECK(r.sol.R(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(0.906285041627).epsilon(1e-9));
        CHECK(r.sol.Tp(0, 0) == doctest::Approx(1.57374156251).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(1.57374156251).epsilon(1e-9));
        // The two populations must add up to N, which no single expression above
        // enforces: Q at the delay and Q at the queue are written from different
        // branches of the assembly.
        CHECK(r.sol.Q(0, 0) + r.sol.Q(1, 0) == doctest::Approx(3.0).epsilon(1e-12));
    }
    SUBCASE("hyperexponential service, the SCV > 1 branch of the two-moment fit") {
        // The size density is not an input: it is reconstructed from the mean
        // and the SCV, and HyperExp(0.9, 2, 0.5) (mean 0.65, SCV 1.9585798817)
        // is what sends sjn_fit down its balanced-means hyperexponential arm.
        qn::Network<double> m =
            sjn_closed("sjnE", 3.0, D::exp_rate(1.0), D::hyperexp(0.9, 2.0, 0.5));
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "sjn.mva");
        CHECK(r.sol.Q(0, 0) == doctest::Approx(1.22067663896).epsilon(1e-9));
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.77932336104).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.793439815327).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(1.45765332459).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(1.22067663896).epsilon(1e-9));
    }
    SUBCASE("an SJF station alongside a PS station") {
        // The remaining stations take the ordinary single-server MVA equation,
        // which is the other half of the recursion and is silent when wrong:
        // dropping the PS station's own queue-length term still leaves every
        // conservation law below satisfied.
        qn::Network<double> m("sjnD");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::SJF);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
        const std::size_t c = m.add_closed_class("C1", 4.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q1, c, D::exp_rate(3.0));
        m.set_service(q2, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q1, 1.0);
        P.set(q1, q2, 1.0);
        P.set(q2, d, 1.0);
        m.link(P);
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "sjn.mva");
        CHECK(r.sol.Q(0, 0) == doctest::Approx(1.56461003043).epsilon(1e-9));
        CHECK(r.sol.Q(1, 0) == doctest::Approx(0.81363969946).epsilon(1e-9));
        CHECK(r.sol.Q(2, 0) == doctest::Approx(1.62175027011).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.521536676811).epsilon(1e-9));
        CHECK(r.sol.U(2, 0) == doctest::Approx(0.782305015217).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(0.520027152858).epsilon(1e-9));
        CHECK(r.sol.R(2, 0) == doctest::Approx(1.03652043548).epsilon(1e-9));
        CHECK(r.sol.Q(0, 0) + r.sol.Q(1, 0) + r.sol.Q(2, 0) ==
              doctest::Approx(4.0).epsilon(1e-12));
    }
}

TEST_CASE("the SJN fixed point is what 'amva' selects, and the utilization cap binds") {
    // N = 6 at rho -> 1 drives the SJN station into the starvation regime, where
    // the conditional waiting time equation underestimates the residence time
    // and the implied throughput would exceed the station capacity. sjn_cap
    // bisects the waiting time until U = umax = 0.999 EXACTLY, which is the
    // number asserted here: it is the cap itself, not a solved value, and a
    // port that capped the throughput instead would break Little's law below.
    qn::Network<double> m = sjn_closed("sjnB", 6.0, D::exp_rate(1.0), D::exp_rate(2.0));
    const mva::DispatchResult<double> r = go(m, "amva");
    CHECK(r.actualmethod == "sjn.amva");
    CHECK(r.sol.Q(0, 0) == doctest::Approx(1.998).epsilon(1e-9));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(4.002).epsilon(1e-9));
    CHECK(r.sol.U(1, 0) == doctest::Approx(0.999).epsilon(1e-9));
    CHECK(r.sol.R(1, 0) == doctest::Approx(2.003003003).epsilon(1e-9));
    CHECK(r.sol.Tp(1, 0) == doctest::Approx(1.998).epsilon(1e-9));
    // X (Z + sum_m C) = N still holds exactly, which is the property the cap is
    // written to preserve.
    CHECK(r.sol.Q(0, 0) + r.sol.Q(1, 0) == doctest::Approx(6.0).epsilon(1e-12));
    CHECK(r.sol.Q(1, 0) == doctest::Approx(r.sol.Tp(1, 0) * r.sol.R(1, 0)).epsilon(1e-9));
    // MATLAB calls line_warning here and still returns. The numbers above are
    // stable and the population law holds exactly, but their ACCURACY is not
    // warranted, so the text must reach the caller rather than be dropped: a
    // silent 1.998 carries an accuracy claim the reference declines to make.
    CHECK_FALSE(r.warning.empty());
    CHECK(r.warning.find("starvation regime") != std::string::npos);
    CHECK(r.warning.find("SolverLDES") != std::string::npos);
}

TEST_CASE("an SJN solve that does not hit the cap carries no warning") {
    // The other side of the boundary: the warning must be set BY the condition,
    // not by taking the SJN branch at all, or it becomes noise that stops
    // meaning anything.
    qn::Network<double> m = sjn_closed("sjnquiet", 3.0, D::exp_rate(1.0), D::exp_rate(2.0));
    const mva::DispatchResult<double> r = go(m);
    CHECK(r.actualmethod == "sjn.mva");
    CHECK(r.warning.empty());
}

TEST_CASE("the SJN multiclass readings: pooled by default, priority when the levels differ") {
    SUBCASE("pooled: equal priorities compare the jobs of every class by size") {
        // Two classes, the second Erlang-2 (SCV 1/2) at the SJF station, which
        // is the branching-Erlang arm of the fit. Equal class priorities, so
        // options.prio stays empty and the sums run over both classes.
        qn::Network<double> m("sjnC");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
        const std::size_t a = m.add_closed_class("CA", 2.0, d);
        const std::size_t b = m.add_closed_class("CB", 2.0, d);
        m.set_service(d, a, D::exp_rate(1.0));
        m.set_service(d, b, D::exp_rate(2.0));
        m.set_service(q, a, D::exp_rate(2.0));
        m.set_service(q, b, D::erlang(2.0, 2));
        qn::RoutingMatrix<double> P;
        P.set(a, a, d, q, 1.0);
        P.set(a, a, q, d, 1.0);
        P.set(b, b, d, q, 1.0);
        P.set(b, b, q, d, 1.0);
        m.link(P);
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "sjn.mva");
        CHECK(r.sol.Q(0, 0) == doctest::Approx(0.738632992435).epsilon(1e-9));
        CHECK(r.sol.Q(0, 1) == doctest::Approx(0.271434874967).epsilon(1e-9));
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.26136700757).epsilon(1e-9));
        CHECK(r.sol.Q(1, 1) == doctest::Approx(1.72856512503).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.369316496217).epsilon(1e-9));
        CHECK(r.sol.U(1, 1) == doctest::Approx(0.542869749934).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(1.7077046659).epsilon(1e-9));
        CHECK(r.sol.R(1, 1) == doctest::Approx(3.18412496781).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 0) == doctest::Approx(0.738632992435).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 1) == doctest::Approx(0.542869749934).epsilon(1e-9));
    }
    SUBCASE("priority: distinct levels select Kant's method A") {
        // Same shape with DISTINCT class priorities, which is what makes
        // solver_mva_sjn_analyzer.m fill options.prio: SJN then applies only
        // WITHIN a class and the classes are non-preemptively prioritised. The
        // numbers differ from the pooled reading on the same demands, so this
        // case is what proves the branch is reached rather than defaulted.
        qn::Network<double> m("sjnF");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
        const std::size_t a = m.add_closed_class("FA", 2.0, d, 1);
        const std::size_t b = m.add_closed_class("FB", 2.0, d, 2);
        m.set_service(d, a, D::exp_rate(1.0));
        m.set_service(d, b, D::exp_rate(2.0));
        m.set_service(q, a, D::exp_rate(2.0));
        m.set_service(q, b, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(a, a, d, q, 1.0);
        P.set(a, a, q, d, 1.0);
        P.set(b, b, d, q, 1.0);
        P.set(b, b, q, d, 1.0);
        m.link(P);
        const mva::DispatchResult<double> r = go(m);
        CHECK(r.actualmethod == "sjn.mva");
        CHECK(r.sol.Q(0, 0) == doctest::Approx(0.847641172963).epsilon(1e-9));
        CHECK(r.sol.Q(0, 1) == doctest::Approx(0.246200704637).epsilon(1e-9));
        CHECK(r.sol.Q(1, 0) == doctest::Approx(1.15235882704).epsilon(1e-9));
        CHECK(r.sol.Q(1, 1) == doctest::Approx(1.75379929536).epsilon(1e-9));
        CHECK(r.sol.U(1, 0) == doctest::Approx(0.423820586481).epsilon(1e-9));
        CHECK(r.sol.U(1, 1) == doctest::Approx(0.492401409274).epsilon(1e-9));
        CHECK(r.sol.R(1, 0) == doctest::Approx(1.35948897221).epsilon(1e-9));
        CHECK(r.sol.R(1, 1) == doctest::Approx(3.56172679917).epsilon(1e-9));
        // The higher-priority class is served first, so its response time is
        // the smaller of the two despite the larger mean job size.
        CHECK(r.sol.R(1, 0) < r.sol.R(1, 1));
    }
}

TEST_CASE("SJN refuses what its recursion does not cover") {
    SUBCASE("an OPEN model is refused by name, never solved size-blind") {
        // The conditional waiting time equation is a POPULATION recursion, so
        // the open case has nothing to recur over. MATLAB's getFeatureSet
        // declares SchedStrategy_SJF unconditionally and mvaDispatch.m carries
        // the restriction as an imperative refusal; a boolean feature cannot
        // say "this discipline, but only in a closed model".
        qn::Network<double> m =
            sqs("sjnopen", SchedStrategy::SJF, D::exp_rate(0.5), D::exp_rate(2.0));
        CHECK_THROWS_AS(go(m), UnsupportedError);
    }
    SUBCASE("a multi-server SJF station is refused: the equation is single server") {
        qn::Network<double> m("sjnms");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
        const std::size_t c = m.add_closed_class("C1", 3.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        m.set_number_of_servers(q, 2.0);
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(go(m), UnsupportedError);
    }
    SUBCASE("a discipline the SJN analyzer cannot solve at the OTHER stations") {
        // DPS is not in the reference's switch, so the whole model is refused
        // rather than the DPS station being quietly solved as PS.
        qn::Network<double> m("sjndps");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::SJF);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::DPS);
        const std::size_t c = m.add_closed_class("C1", 3.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q1, c, D::exp_rate(3.0));
        m.set_service(q2, c, D::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q1, 1.0);
        P.set(q1, q2, 1.0);
        P.set(q2, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(go(m), UnsupportedError);
    }
    SUBCASE("exact arithmetic, which the quadrature-based recursion cannot serve") {
        // The recursion evaluates regularized incomplete gammas on a Simpson
        // grid, so it is refused by name under Rational rather than returning a
        // rounded double lifted back into an exact type.
        qn::Network<Rational> m("sjnexact");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::SJF);
        const std::size_t c = m.add_closed_class("C1", 3.0, d);
        m.set_service(d, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
        m.set_service(q, c, Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
        qn::RoutingMatrix<Rational> P;
        P.set(d, q, num_traits<Rational>::from_int(1));
        P.set(q, d, num_traits<Rational>::from_int(1));
        m.link(P);
        mva::MvaOptions opt;
        Matrix<Rational> init;
        CHECK_THROWS_AS(mva::mva_dispatch(m.get_struct(), opt, init), UnsupportedError);
    }
}

}  // namespace
