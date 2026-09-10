/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The CLOSED setup/delay-off queue, BUG-78.
 *
 * A closed queue carrying a setup and a delay-off was solved by the per-instance
 * cold-start race `R = p_cold*E[setup] + S`: it raced the delay-off against the
 * per-instance idle time and carried NO queueing term, so it described a
 * serverless instance pool rather than a single-server vacation queue. The
 * measured symptom was that the reported response time was byte-identical across
 * a tenfold change in the setup mean while two independent simulators moved from
 * 0.90 to 3.08.
 *
 * `qbd_setupdelayoff_closed` solves the finite level-dependent chain the
 * simulator actually walks: lambda(n) = (N-n)/Z with the level bounded by N, the
 * setup phases above level 0 and the delay-off phases at level 0. An arrival
 * during the delay-off finds the server still warm and resumes without setup,
 * which is Solver_ssj's cancelDelayoff and the discipline the analysis has to
 * reproduce.
 *
 * THE REFERENCE IS THE SIMULATORS, since no analytical row is available: the
 * model is Delay(Z=1) + Queue(FCFS, D=0.5), one closed class with N=3 and a
 * delay-off of Exp(4), and SolverCTMC refuses it by featset. LDES reads
 * 0.9009 / 1.1157 / 3.0842 and JMT 0.9002 / 1.1137 / 3.0703 at setup means
 * none / 0.5 / 5.0. The bar below is the bracket between them, not a tolerance
 * around one of them: the two are 0.2 per cent apart and the analysis has to sit
 * inside that rather than merely near it.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/qbd_setupdelayoff.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_runner.h"

namespace qn = line::qn;
namespace mam = line::mam;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** Delay(Exp(1)) -> Queue(FCFS, Exp(2)), one closed class of 3 jobs. */
qn::Network<double> closed_setup(const D* setup, const D* delayoff) {
    qn::Network<double> m("closed_setup");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    if (setup != 0) m.set_setup_delayoff(q, c, *setup, *delayoff);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

double resp_at_queue(const D* setup, const D* delayoff) {
    qn::Network<double> m = closed_setup(setup, delayoff);
    mam::MamOptions mo;
    const line::mva::AvgResult<double> r = mam::solver_mam_run_analyzer(m.get_struct(), mo);
    return r.RN(1, 0);
}

double tput_at_queue(const D* setup, const D* delayoff) {
    qn::Network<double> m = closed_setup(setup, delayoff);
    mam::MamOptions mo;
    const line::mva::AvgResult<double> r = mam::solver_mam_run_analyzer(m.get_struct(), mo);
    return r.TN(1, 0);
}

}  // namespace

TEST_CASE("qbd_setupdelayoff_closed lands between the two simulators") {
    const mam::SetupDelayoffClosed<double> fast =
        mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 2.0, 1.0, 4.0, 1.0);
    const double rFast = fast.QN / fast.XN;
    CHECK(rFast >= 1.1137 * 0.995);
    CHECK(rFast <= 1.1157 * 1.005);

    const mam::SetupDelayoffClosed<double> slow =
        mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 0.2, 1.0, 4.0, 1.0);
    const double rSlow = slow.QN / slow.XN;
    CHECK(rSlow >= 3.0703 * 0.995);
    CHECK(rSlow <= 3.0842 * 1.005);
}

TEST_CASE("qbd_setupdelayoff_closed degenerates to the plain closed queue") {
    // An instantaneous setup is a server that is never cold, so the vacation
    // chain has to collapse onto the exact closed queue: Q = 1.421053,
    // X = 1.578947, R = 0.9, which is what MVA, NC and CTMC all return.
    const mam::SetupDelayoffClosed<double> instant =
        mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 1e8, 1.0, 4.0, 1.0);
    CHECK(instant.QN == doctest::Approx(1.42105263).epsilon(1e-8));
    CHECK(instant.XN == doctest::Approx(1.57894737).epsilon(1e-8));

    // A delay-off that never expires is the same statement from the other side:
    // the server is always warm, so no arrival ever pays a setup.
    const mam::SetupDelayoffClosed<double> warm =
        mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 0.2, 1.0, 1e-8, 1.0);
    CHECK(warm.QN == doctest::Approx(1.42105263).epsilon(1e-6));
}

TEST_CASE("qbd_setupdelayoff_closed reads the setup SCV") {
    // A non-exponential setup must not be read as exponential: the Coxian form
    // carries the SCV, and a lower SCV is a more predictable -- and so cheaper --
    // cold start.
    const double rExp =
        [] {
            const mam::SetupDelayoffClosed<double> r =
                mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 0.2, 1.0, 4.0, 1.0);
            return r.QN / r.XN;
        }();
    const double rLow =
        [] {
            const mam::SetupDelayoffClosed<double> r =
                mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 0.2, 0.5, 4.0, 1.0);
            return r.QN / r.XN;
        }();
    const double rHigh =
        [] {
            const mam::SetupDelayoffClosed<double> r =
                mam::qbd_setupdelayoff_closed<double>(3.0, 1.0, 2.0, 0.2, 2.0, 4.0, 1.0);
            return r.QN / r.XN;
        }();
    CHECK(rLow < rExp);
    CHECK(rExp < rHigh);
}

TEST_CASE("SolverMAM does not ignore a closed setup") {
    const D fast = D::exp_rate(2.0);      // setup mean 0.5
    const D slow = D::exp_rate(0.2);      // setup mean 5.0
    const D doff = D::exp_rate(4.0);

    const double rPlain = resp_at_queue(0, 0);
    const double rFast = resp_at_queue(&fast, &doff);
    const double rSlow = resp_at_queue(&slow, &doff);

    // The defect this pins: every analytical row was byte-identical across a
    // tenfold change in the setup mean while the simulators moved 0.90 -> 3.08.
    CHECK(rPlain < rFast);
    CHECK(rFast < rSlow);
    CHECK(rSlow / rPlain > 3.0);

    CHECK(rPlain == doctest::Approx(0.9).epsilon(1e-6));
    CHECK(rFast >= 1.1137 * 0.995);
    CHECK(rFast <= 1.1157 * 1.005);
    CHECK(rSlow >= 3.0703 * 0.995);
    CHECK(rSlow <= 3.0842 * 1.005);
}

TEST_CASE("a closed setup slows the chain and conserves its population") {
    const D slow = D::exp_rate(0.2);
    const D doff = D::exp_rate(4.0);

    // A slower setup is a slower server, so the closed chain's throughput has to
    // drop with it; the p_cold formula left it at the setup-free 1.578947.
    CHECK(tput_at_queue(0, 0) == doctest::Approx(1.578947).epsilon(1e-5));
    CHECK(tput_at_queue(&slow, &doff) < tput_at_queue(0, 0));

    // Little across the two stations: a job is thinking or it is at the queue.
    qn::Network<double> m = closed_setup(&slow, &doff);
    mam::MamOptions mo;
    const line::mva::AvgResult<double> r = mam::solver_mam_run_analyzer(m.get_struct(), mo);
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-9));
}
