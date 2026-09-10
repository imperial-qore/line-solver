/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * One load-dependent fixture, every solver, ONE utilization convention.
 *
 * Utilization has two readings that COINCIDE for a load-independent single
 * server and part company the moment alpha(n) != 1:
 *
 *   work-based   U = X * E[S] / max(c, max alpha)   the fraction of the
 *                                                   station's PEAK capacity
 *                                                   actually being delivered
 *   time-based   U = P(at least one server busy)    the fraction of time the
 *                                                   station is occupied
 *
 * A server running alpha(n) times faster does the same work in less time, so
 * the time-based reading calls it no busier than one at its nominal rate. On
 * the fixture below that is 0.9587 against 0.6612 -- a 45% spread on the same
 * model, with QLen and Tput agreeing to every printed digit.
 *
 * LINE reports the WORK-BASED number everywhere. It used to be split: MAM's
 * LD-QBD, the NRM SSA engine and the C++ LDES engine measured busy time while
 * CTMC, MVA, NC, serial SSA and the Java LDES engine measured work, so
 * SolverSSA disagreed with ITSELF depending on which engine ran.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mam/solver_mam_ldqbd.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ldes/ldes_engine.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace mam = line::mam;
namespace ctmc = line::ctmc;
namespace mva = line::mva;
namespace nc = line::nc;
namespace ssa = line::ssa;
namespace ldes = line::ldes;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

const std::size_t QI = 1;   // 0 = Delay, 1 = Queue

/**
 * Closed form for the fixture: the queue is a birth-death chain on n = 0..4
 * with birth (N-n)*1.0 and death alpha(n)*1.0.
 */
const double X_EXACT = 1.6528925619834711;
const double U_EXACT = X_EXACT / 2.5;    // 0.6611570247933884
const double Q_EXACT = 2.3471074380165293;

qn::Network<double> lld_model(int N, double think_rate, const std::vector<double>& alpha) {
    qn::Network<double> m("lld_util");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", static_cast<double>(N), d);
    m.set_service(d, c, D::exp_rate(think_rate));
    m.set_service(q, c, D::exp_rate(1.0));
    m.set_load_dependence(q, alpha);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

const std::vector<double> ALPHA{1.0, 1.5, 2.0, 2.5};

}  // namespace

TEST_CASE("lld util: the closed form is what CTMC computes") {
    // anchors the other cases to arithmetic rather than to a solver
    qn::Network<double> m = lld_model(4, 1.0, ALPHA);
    const line::mva::AvgResult<double> e =
        ctmc::solver_ctmc_run_analyzer(m.get_struct(), ctmc::CtmcOptions());
    CHECK(e.QN(QI, 0) == doctest::Approx(Q_EXACT).epsilon(1e-9));
    CHECK(e.UN(QI, 0) == doctest::Approx(U_EXACT).epsilon(1e-9));
    CHECK(e.TN(QI, 0) == doctest::Approx(X_EXACT).epsilon(1e-9));
    // and it is emphatically NOT the time-based reading
    CHECK(std::abs(e.UN(QI, 0) - 0.9586776860) > 0.29);
}

TEST_CASE("lld util: MAM's LD-QBD is an exact chain and agrees exactly") {
    // the regression: it used to report 1 - p(0) = 0.9587 here
    qn::Network<double> m = lld_model(4, 1.0, ALPHA);
    mam::MamOptions mo;
    mo.method = "ldqbd";
    const mam::LdqbdSolution<double> s = mam::solver_mam_ldqbd(m.get_struct(), mo);
    CHECK(s.sol.Q(QI, 0) == doctest::Approx(Q_EXACT).epsilon(1e-9));
    CHECK(s.sol.U(QI, 0) == doctest::Approx(U_EXACT).epsilon(1e-9));
    CHECK(s.sol.Tp(QI, 0) == doctest::Approx(X_EXACT).epsilon(1e-9));
}

TEST_CASE("lld util: MVA and NC report the same number") {
    qn::Network<double> mm = lld_model(4, 1.0, ALPHA);
    mva::MvaOptions mopt;
    line::Matrix<double> init;
    const mva::AvgResult<double> a = mva::solver_mva_run_analyzer(mm.get_struct(), mopt, init);
    CHECK(a.QN(QI, 0) == doctest::Approx(Q_EXACT).epsilon(1e-9));
    CHECK(a.UN(QI, 0) == doctest::Approx(U_EXACT).epsilon(1e-9));

    qn::Network<double> mn = lld_model(4, 1.0, ALPHA);
    const mva::AvgResult<double> b = nc::solver_nc_run_analyzer(mn.get_struct(), nc::NcSolverOptions());
    CHECK(b.QN(QI, 0) == doctest::Approx(Q_EXACT).epsilon(1e-9));
    CHECK(b.UN(QI, 0) == doctest::Approx(U_EXACT).epsilon(1e-9));
}

TEST_CASE("lld util: both SSA engines report the same convention") {
    // the two measure it differently -- NRM integrates busy time and serial
    // applies the utilization law -- and used to differ by 0.30 here
    for (const char* method : {"nrm", "serial"}) {
        CAPTURE(method);
        qn::Network<double> m = lld_model(4, 1.0, ALPHA);
        ssa::SsaOptions opt;
        opt.method = method;
        opt.samples = 200000;
        opt.seed = 23000;
        const ssa::SsaSolution sim = ssa::solver_ssa(m.get_struct(), opt);
        CHECK(sim.UN(QI, 0) == doctest::Approx(U_EXACT).epsilon(3e-2));
        CHECK(sim.QN(QI, 0) == doctest::Approx(Q_EXACT).epsilon(3e-2));
    }
}

TEST_CASE("lld util: utilization stays below one while capacity remains") {
    // the time-based reading saturates at 1 long before the work-based one does
    const std::vector<double> ramp{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    qn::Network<double> mc = lld_model(6, 100.0, ramp);   // a near-zero think time
    const line::mva::AvgResult<double> e =
        ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());

    qn::Network<double> mm = lld_model(6, 100.0, ramp);
    mam::MamOptions mo;
    mo.method = "ldqbd";
    const mam::LdqbdSolution<double> s = mam::solver_mam_ldqbd(mm.get_struct(), mo);

    CHECK(s.sol.U(QI, 0) == doctest::Approx(e.UN(QI, 0)).epsilon(1e-9));
    CHECK(s.sol.U(QI, 0) > 0.0);
    CHECK(s.sol.U(QI, 0) < 1.0);
}

TEST_CASE("lld util: the native LDES engine measures work, not busy time") {
    // this engine integrates the BUSY TIME of each server along the sample
    // path, so a load-dependent station needed the integrand scaled by the
    // alpha(n) in force over each interval and the divisor raised from the
    // server count to max(c, max alpha). Without both it reported 0.9598 here,
    // and at a load-dependent PS station -- whose shares already carried the
    // alpha -- it could report a utilization ABOVE ONE.
    qn::Network<double> m = lld_model(4, 1.0, ALPHA);
    ldes::LdesOptions o;
    o.samples = 200000;
    o.seed = 23000;
    const ldes::LdesResult r = ldes::ldes_engine_solve(m.get_struct(), o);
    CHECK(r.UN(QI, 0) == doctest::Approx(U_EXACT).epsilon(3e-2));
    CHECK(r.QN(QI, 0) == doctest::Approx(Q_EXACT).epsilon(3e-2));
    CHECK(r.TN(QI, 0) == doctest::Approx(X_EXACT).epsilon(3e-2));
    CHECK(r.UN(QI, 0) < 1.0);
}
