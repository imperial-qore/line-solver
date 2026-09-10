/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Offered-vs-carried arrival rate (ArvR) at a bound in solver_ctmc.
 *
 * A finite bound on an open class is EITHER a physical capacity (set_capacity)
 * OR a state-space cutoff imposed only to keep the CTMC enumeration finite. They
 * must report ArvR differently, and the discriminator is the station drop rule
 * (is_physical_capacity / arrival_is_lost in lang/qn/state_events.h, applied by
 * the state-events generator the CTMC analyzer calls):
 *
 *   PHYSICAL cap  -> a job meeting a full buffer is really lost, so the loss is a
 *                    real event: ArvR reports the OFFERED rate (lambda), Tput<lambda.
 *   state-space CUTOFF -> the refused job never existed; a truncation artifact
 *                    that counts nowhere: ArvR reports the CARRIED rate (== Tput).
 *
 * This is the USER DECISION of 2026-07-17 (a cutoff is not a physical capacity),
 * documented in _kb/06-solver-catalog.md, and mirrors the SSA convention (BUG-85).
 * test_ctmc_analyzer.cpp already pins the PHYSICAL-cap side (M/M/1/K reports the
 * offered 0.6); this file adds the PURE-CUTOFF side and the decisive control that
 * separates the two, which no test covered before.
 *
 * The BINDING control is that K=2 (physical) and cutoff=2 (truncation) yield an
 * IDENTICAL stationary distribution (M/M/1/2 == M/M/1 truncated at 2, both Tput
 * 0.630996) yet must report DIFFERENT ArvR: 0.9 (cap) vs 0.630996 (cutoff).
 */
#include <algorithm>
#include <string>

#include "doctest.h"
#include "line/num/number.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

constexpr double LAM = 0.9;
constexpr double MU = 1.0;

/** Source(Exp lambda) -> FCFS Queue -> Sink, one open class; K>0 sets a cap. */
qn::Network<double> mm1(int K) {
    qn::Network<double> m("mm1");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(LAM));
    m.set_service(q, c, Dist::exp_rate(MU));
    if (K > 0) m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** {ArvR, Tput} at the Queue (station row 1, class 0). */
std::pair<double, double> solve(int K, int cutoff) {
    qn::Network<double> m = mm1(K);
    ctmc::CtmcOptions opt;
    opt.cutoff = cutoff;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);
    return std::make_pair(r.AN(1, 0), r.TN(1, 0));
}

}  // namespace

TEST_CASE("a state-space cutoff reports the carried ArvR") {
    // infinite capacity + tight cutoff: ArvR == carried Tput, both < lambda.
    const int cutoffs[] = {2, 3, 5, 8};
    for (int cutoff : cutoffs) {
        const std::pair<double, double> r = solve(-1, cutoff);
        CHECK(r.first == doctest::Approx(r.second).epsilon(1e-9));
        CHECK(r.first < LAM - 1e-6);
    }
}

TEST_CASE("a physical capacity reports the offered ArvR") {
    // finite physical capacity K + loose cutoff: ArvR == offered lambda, Tput < lambda.
    const int Ks[] = {2, 3, 5};
    for (int K : Ks) {
        const std::pair<double, double> r = solve(K, std::max(K + 5, 20));
        CHECK(r.first == doctest::Approx(LAM).epsilon(1e-9));
        CHECK(r.second < LAM - 1e-6);
    }
}

TEST_CASE("cutoff and cap share a distribution but differ in ArvR") {
    // Decisive control: K=2 (physical) and cutoff=2 (truncation) share the same
    // stationary throughput but must report DIFFERENT ArvR.
    const std::pair<double, double> cut = solve(-1, 2);
    const std::pair<double, double> cap = solve(2, 20);
    CHECK(cut.second == doctest::Approx(cap.second).epsilon(1e-9));  // identical carried rate
    CHECK(cut.first == doctest::Approx(cut.second).epsilon(1e-9));   // cutoff -> carried
    CHECK(cap.first == doctest::Approx(LAM).epsilon(1e-9));          // cap -> offered
    CHECK(cap.first > cut.first + 1e-3);                             // and they differ
}
