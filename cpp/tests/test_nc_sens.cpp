/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `@@NetworkSolver/getSensitivityTable.m` under SolverNC: the analytic branch
 * SolverNC is one of the two engines entitled to, and the finite differences
 * every engine falls back to.
 *
 * WHY THE TWO BRANCHES ARE CHECKED AGAINST EACH OTHER. The exact branch
 * differentiates the product-form recursion in closed form, so it shares no code
 * with the difference quotient, which re-solves rate-perturbed copies of the
 * model through the normalizing-constant path itself. Agreement to the order of
 * the step is therefore evidence about the derivative and not about one
 * implementation reproducing its own error. The scope conditions are checked by
 * the refusals they cause, since an out-of-scope model that answered anyway
 * would report a number no branch computed.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/sens/solver_sens_table.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

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

/** Source -> Queue1 -> Queue2 -> Sink, one open class: the open exact branch. */
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

/** Delay + a 2-server FCFS queue: out of the analytic branch's scope. */
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

/** The table SolverNC produces, i.e. `getAvg` re-run on the perturbed struct. */
sens::SensTable<double> nc_sens(qn::Network<double>& model, const std::string& method,
                                const std::string& scheme = "forward") {
    qn::NetworkStruct<double> sn = model.get_struct();
    nc::NcSolverOptions opt;
    sens::SensOptions so;
    so.method = method;
    so.scheme = scheme;
    return sens::solver_sensitivity_table<double>(sn, so, /*exact_available=*/true, [&sn, &opt]() {
        const mva::AvgResult<double> a = nc::solver_nc_run_analyzer(sn, opt);
        mva::MvaSolution<double> s;
        s.Q = a.QN;
        s.U = a.UN;
        s.R = a.RN;
        s.Tp = a.TN;
        s.C = a.CN;
        s.X = a.XN;
        s.method = a.actualmethod;
        return s;
    });
}

const sens::SensRow<double>& row_of(const sens::SensTable<double>& t, const std::string& station,
                                    const std::string& jobclass) {
    for (const sens::SensRow<double>& r : t.rows)
        if (r.station == station && r.jobclass == jobclass) return r;
    FAIL("no row for " << station << " / " << jobclass);
    return t.rows[0];
}

}  // namespace

TEST_CASE("NC sensitivity: the closed analytic branch agrees with central differences") {
    qn::Network<double> m = closed_multiclass();
    const sens::SensTable<double> exact = nc_sens(m, "exact");
    const sens::SensTable<double> fd = nc_sens(m, "fd", "central");

    CHECK(exact.method == "exact");
    CHECK(fd.method == "fd");
    // The analytic branch of a closed model carries the pfqn_sens Jacobian; the
    // difference quotient has none to carry, which is the reference's contract.
    CHECK(exact.has_jacobian);
    CHECK_FALSE(fd.has_jacobian);
    REQUIRE(exact.rows.size() == fd.rows.size());
    REQUIRE(exact.rows.size() == 4);  // two queues x two classes, the delay excluded

    for (std::size_t i = 0; i < exact.rows.size(); ++i) {
        const sens::SensRow<double>& a = exact.rows[i];
        const sens::SensRow<double>& b = fd.rows[i];
        CAPTURE(a.station);
        CAPTURE(a.jobclass);
        CHECK(a.station == b.station);
        CHECK(a.jobclass == b.jobclass);
        CHECK(a.dTput == doctest::Approx(b.dTput).epsilon(1e-3));
        CHECK(a.dRespT == doctest::Approx(b.dRespT).epsilon(1e-3));
        CHECK(a.dQLen == doctest::Approx(b.dQLen).epsilon(1e-3));
        CHECK(a.dUtil == doctest::Approx(b.dUtil).epsilon(1e-3));
    }

    // Raising a station's own service rate cannot raise its own queue length,
    // its response time or its utilization, and cannot lower its throughput.
    for (const sens::SensRow<double>& r : exact.rows) {
        CAPTURE(r.station);
        CAPTURE(r.jobclass);
        CHECK(r.dQLen <= 0.0);
        CHECK(r.dRespT <= 0.0);
        CHECK(r.dUtil <= 0.0);
        CHECK(r.dTput >= 0.0);
    }
}

TEST_CASE("NC sensitivity: the open analytic branch agrees with central differences") {
    qn::Network<double> m = open_tandem();
    const sens::SensTable<double> exact = nc_sens(m, "exact");
    const sens::SensTable<double> fd = nc_sens(m, "fd", "central");

    CHECK(exact.method == "exact");
    // The open branch differentiates a closed form station by station and forms
    // no Jacobian, the reference returning an empty `sens` there.
    CHECK_FALSE(exact.has_jacobian);
    REQUIRE(exact.rows.size() == 2);
    REQUIRE(fd.rows.size() == 2);

    for (std::size_t i = 0; i < exact.rows.size(); ++i) {
        const sens::SensRow<double>& a = exact.rows[i];
        const sens::SensRow<double>& b = fd.rows[i];
        CAPTURE(a.station);
        CHECK(a.dRespT == doctest::Approx(b.dRespT).epsilon(1e-3));
        CHECK(a.dQLen == doctest::Approx(b.dQLen).epsilon(1e-3));
        CHECK(a.dUtil == doctest::Approx(b.dUtil).epsilon(1e-3));
    }

    // An open throughput is the arrival rate times the visits and does not move
    // with a service rate, so the column is EXACTLY zero rather than small.
    for (const sens::SensRow<double>& r : exact.rows) CHECK(r.dTput == 0.0);

    // M/M/1 in closed form: with lambda = 0.5 and mu = 2, U = 1/4 and
    // dU/dmu = -rho/mu = -1/8, dR/dmu = -(1-rho+rho)/(mu^2 (1-rho)^2) is
    // -1/(mu^2(1-rho)^2) + ... which the branch composes; check the two whose
    // algebra is unambiguous.
    const sens::SensRow<double>& r1 = row_of(exact, "Queue1", "O1");
    CHECK(r1.dUtil == doctest::Approx(-0.125).epsilon(1e-12));
    const sens::SensRow<double>& r2 = row_of(exact, "Queue2", "O1");
    CHECK(r2.dUtil == doctest::Approx(-0.5).epsilon(1e-12));
}

TEST_CASE("NC sensitivity: out of scope refuses rather than downgrading") {
    qn::Network<double> ms = closed_multiserver();
    // A multiserver station has no single-server recursion to differentiate.
    CHECK_THROWS_AS(nc_sens(ms, "exact"), UnsupportedError);
    // ... and `auto` falls back to the differences instead of refusing.
    const sens::SensTable<double> a = nc_sens(ms, "auto");
    CHECK(a.method == "fd");
    CHECK(a.rows.size() == 1);

    qn::Network<double> ok = closed_multiclass();
    CHECK_THROWS_AS(nc_sens(ok, "bisection"), InputError);
    CHECK_THROWS_AS(nc_sens(ok, "fd", "backward"), InputError);
}
