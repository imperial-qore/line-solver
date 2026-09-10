/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The SPN firing-phase progression: `after_event_transition` and the per-mode
 * PHASE synchronizations of `refresh_sync`.
 *
 * THE ORACLE IS MATLAB SolverCTMC on the identical net, run 2026-08-05:
 * P1 -> T1 -> P2 -> T2 -> P1, one token, T1 firing Erlang(4,2) (mean 0.5),
 * T2 firing Exp(3): QLen = {0.6, 0.4}, Tput = 1.2 at both places, which is
 * also the analytic renewal split 0.5 : 1/3 of the 5/6 cycle. Before the
 * PHASE port a multi-phase firing process could START in a phase (the ENABLE
 * seeding) but never ADVANCE, so the Erlang net either deadlocked or solved
 * with a one-phase firing law.
 */

#include <limits>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::Distrib;
using line::lang::EventType;
using Dist = Distrib<double>;

namespace {
qn::Network<double> erlang_spn() {
    qn::Network<double> m("spn_erl");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, Dist::exp_rate(1.0));
    m.set_service(p2, c, Dist::exp_rate(1.0));

    // node-indexed arc tables sized over {P1, P2, T1, T2}
    qn::TransitionParam<double> t1;
    t1.nmodes = 1;
    t1.modenames.push_back("M1");
    t1.enabling.assign(1, line::Matrix<double>(4, 1, 0.0));
    t1.inhibiting.assign(1, line::Matrix<double>(4, 1, std::numeric_limits<double>::infinity()));
    t1.firing.assign(1, line::Matrix<double>(4, 1, 0.0));
    t1.enabling[0](p1 - 1, 0) = 1.0;
    t1.firing[0](p2 - 1, 0) = 1.0;
    t1.nmodeservers.push_back(1.0);
    t1.firingphases.push_back(2);
    t1.timing.push_back(line::lang::TimingStrategy::TIMED);
    t1.fireweight.push_back(1.0);
    t1.firingproc.push_back(Dist::erlang(4.0, 2));
    m.add_transition("T1", t1);

    qn::TransitionParam<double> t2;
    t2.nmodes = 1;
    t2.modenames.push_back("M2");
    t2.enabling.assign(1, line::Matrix<double>(4, 1, 0.0));
    t2.inhibiting.assign(1, line::Matrix<double>(4, 1, std::numeric_limits<double>::infinity()));
    t2.firing.assign(1, line::Matrix<double>(4, 1, 0.0));
    t2.enabling[0](p2 - 1, 0) = 1.0;
    t2.firing[0](p1 - 1, 0) = 1.0;
    t2.nmodeservers.push_back(1.0);
    t2.firingphases.push_back(1);
    t2.timing.push_back(line::lang::TimingStrategy::TIMED);
    t2.fireweight.push_back(1.0);
    t2.firingproc.push_back(Dist::exp_rate(3.0));
    m.add_transition("T2", t2);
    return m;
}
}  // namespace

TEST_CASE("refresh_sync emits one PHASE action per transition mode") {
    qn::Network<double> m = erlang_spn();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::Sync<double>> sync = qn::refresh_sync(sn);
    std::size_t nphase = 0;
    for (const qn::Sync<double>& s : sync)
        if (s.active.event == EventType::PHASE &&
            sn.nodes[s.active.node - 1].nodetype == qn::NodeType::Transition)
            ++nphase;
    CHECK(nphase == 2);  // one per mode of T1 and of T2
}

TEST_CASE("the Erlang firing SPN reproduces the MATLAB CTMC means") {
    qn::Network<double> m = erlang_spn();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> r = ctmc::solver_ctmc_analyzer(sn, opt);
    // MATLAB SolverCTMC: QLen 0.6 / 0.4 (Tput 1.2 there; this port does not
    // yet attribute firing completions to place throughput on ANY SPN, also
    // with single-phase modes, so TN is not asserted here -- see _kb).
    CHECK(r.avg.QN(0, 0) == doctest::Approx(0.6).epsilon(1e-8));
    CHECK(r.avg.QN(1, 0) == doctest::Approx(0.4).epsilon(1e-8));
}
