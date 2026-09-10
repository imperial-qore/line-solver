/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `advanced/agentModels`: SolverAG's agent-based INAP arm.
 *
 * SolverAG decomposes a network into
 * interacting agents and closes the loop with a fixed point over the
 * synchronization rates. `method = 'inap'` selects `solver_ag`, which is the
 * analyzer these five examples exercise; each then re-solves the same model
 * with the solver the reference compares it against, so the two numbers are
 * printed side by side under their own banners.
 *
 * These methods used to be reached as `MAM(model, 'inap')`; they are SolverAG's
 * since those methods moved there, and `mam::check_method` now refuses them by name.
 *
 * `exact` (AutoCAT) is NOT one of them: it needs an LP/NLP solver and is
 * refused by name in `solver_ag_autocat.h`, exactly as the reference's own
 * note in `ag_tandem_open.py` says.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/ag/ag_dispatch.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

/** `AG(model, 'inap')`, the agent-based arm every example here runs first. */
void run_ag_inap(const Sn& sn) {
    ag::AgOptions o;
    o.method = "inap";
    const mva::AvgResult<double> r = ag::solver_ag_run_analyzer(sn, o);
    std::printf("AG (method=%s):\n", r.actualmethod.c_str());
    print_avg(sn, r);
}

void run_mva(const Sn& sn) {
    mva::MvaOptions o;
    Matrix<double> init;
    print_avg(sn, mva::solver_mva_run_analyzer(sn, o, init));
}

void run_ctmc(const Sn& sn) {
    ctmc::CtmcOptions o;
    print_avg(sn, ctmc::solver_ctmc_run_analyzer(sn, o));
}

}  // namespace

// ---------------------------------------------------------------------------
// ag_closed_network
// ---------------------------------------------------------------------------

/** A closed two-queue PS loop under the agent-based INAP method, against MVA and exact CTMC. */
void ag_closed_network() {
    const double N = 10.0, mu1 = 2.0, mu2 = 1.0;

    Net m("Closed-2Q");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass c(m, "Class1", N, q1);
    q1.set_service(c, Exp(mu1));
    q2.set_service(c, Exp(mu2));
    Routing P;
    cyclic(P, c, {q1, q2});
    m.link(P);

    note("=== Closed Network ===\n");
    const Sn& sn = m.get_struct();
    section("AG");
    run_ag_inap(sn);
    section("MVA");
    run_mva(sn);
    // The reference guards the exact solve with N <= 20, which holds here.
    section("CTMC");
    run_ctmc(sn);
}

// ---------------------------------------------------------------------------
// ag_multiclass_closed
// ---------------------------------------------------------------------------

/** Two closed classes over two PS queues: the solver builds one agent per (station, class). */
void ag_multiclass_closed() {
    const double N1 = 5.0, N2 = 3.0;
    const double mu1_q1 = 2.0, mu2_q1 = 1.5, mu1_q2 = 1.0, mu2_q2 = 0.8;

    Net m("Multiclass-Closed");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    ClosedClass c1(m, "Class1", N1, q1);
    ClosedClass c2(m, "Class2", N2, q1);
    q1.set_service(c1, Exp(mu1_q1));
    q1.set_service(c2, Exp(mu2_q1));
    q2.set_service(c1, Exp(mu1_q2));
    q2.set_service(c2, Exp(mu2_q2));
    Routing P;
    cyclic(P, c1, {q1, q2});
    cyclic(P, c2, {q1, q2});
    m.link(P);

    note("=== Multiclass Closed Network ===\n");
    const Sn& sn = m.get_struct();
    section("AG");
    run_ag_inap(sn);
    section("MVA");
    run_mva(sn);
}

// ---------------------------------------------------------------------------
// ag_tandem_open
// ---------------------------------------------------------------------------

/** An open M/M/1 -> M/M/1 tandem, whose exact answer is closed form. */
void ag_tandem_open() {
    const double lam = 0.5, s1 = 1.0, s2 = 1.5;

    Net m("Tandem-MM1");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c(m, "Class1");
    src.set_arrival(c, Exp(lam));
    q1.set_service(c, Exp(s1));
    q2.set_service(c, Exp(s2));
    Routing P;
    serial(P, c, {src, q1, q2, snk});
    m.link(P);

    const double U1 = lam / s1, U2 = lam / s2;
    note("=== Open Tandem Queue (M/M/1 -> M/M/1) ===\n");
    note("Analytical (M/M/1):");
    std::printf("  Queue1: U=%.4f, Q=%.4f, R=%.4f\n", U1, U1 / (1.0 - U1), 1.0 / (s1 - lam));
    std::printf("  Queue2: U=%.4f, Q=%.4f, R=%.4f\n\n", U2, U2 / (1.0 - U2), 1.0 / (s2 - lam));

    const Sn& sn = m.get_struct();
    section("AG");
    run_ag_inap(sn);
    // NOT a gap in this port: the reference never calls 'exact', it only notes
    // why, so there is no call to record and nothing to refuse. The C++ side
    // declines it for the same reason (solver_ag_autocat.h needs an LP/NLP
    // solver), which is the reference's situation and not a divergence from it.
    note("Note: 'exact' method (AutoCAT) is not yet implemented in native Python");
    note("It falls back to INAP with a warning in the JAR version");
    section("MVA");
    run_mva(sn);
}

// ---------------------------------------------------------------------------
// ag_tandem_phasetype
// ---------------------------------------------------------------------------

/**
 * The same open tandem with PHASE-TYPE service, which the agent-based methods represent exactly
 * rather than collapsing to a mean rate.
 *
 * Each component is a QBD over (queue length, phase), so Queue1 -- an isolated
 * M/PH/1, since it sees the Poisson source directly -- comes out at the
 * Pollaczek-Khinchine mean whatever the reversed-rate iteration does. Both
 * stations carry the SAME mean service time and differ only in variability,
 * which is exactly what the earlier scalar birth-death construction could not
 * see: it returned the M/M/1 answer for either.
 */
void ag_tandem_phasetype() {
    const double lam = 0.5, mean_s = 1.0, scv1 = 0.5, scv2 = 4.0;

    Net m("Tandem-MPH1");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c(m, "Class1");
    src.set_arrival(c, Exp(lam));
    q1.set_service(c, D::erlang_fit(mean_s, scv1));
    q2.set_service(c, lang::hyperexp_fit_mean_scv<double>(mean_s, scv2));
    Routing P;
    serial(P, c, {src, q1, q2, snk});
    m.link(P);

    const double rho = lam * mean_s;
    const double pk1 = rho + rho * rho * (1.0 + scv1) / (2.0 * (1.0 - rho));
    note("=== Open Tandem with Phase-Type Service ===\n");
    note("Queue1 is an isolated M/Er2/1, so its exact mean queue length is the");
    std::printf("Pollaczek-Khinchine value %.6f. The M/M/1 reading would be %.6f.\n\n",
                pk1, rho / (1.0 - rho));

    const Sn& sn = m.get_struct();
    section("AG");
    run_ag_inap(sn);
    // 'inapinf' additionally drops the maxStates truncation, solving each open
    // component on its infinite state space through Neuts' rate matrix R. That
    // matters most at Queue2, whose service law has the heavier tail.
    section("AG inapinf");
    ag::AgOptions oinf;
    oinf.method = "inapinf";
    print_avg(sn, ag::solver_ag_run_analyzer(sn, oinf));
}

// ---------------------------------------------------------------------------
// ag_jackson_network
// ---------------------------------------------------------------------------

/** An open three-queue Jackson network with probabilistic routing. */
void ag_jackson_network() {
    const double lam = 1.0;
    const double mu[3] = {2.0, 3.0, 2.5};
    const double p12 = 0.4, p13 = 0.3, p1s = 0.3;
    const double p21 = 0.2, p23 = 0.3, p2s = 0.5;
    const double p3s = 1.0;

    Net m("Jackson-3Q");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Queue q3(m, "Queue3", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c(m, "Class1");
    src.set_arrival(c, Exp(lam));
    q1.set_service(c, Exp(mu[0]));
    q2.set_service(c, Exp(mu[1]));
    q3.set_service(c, Exp(mu[2]));
    Routing P;
    P.set(c, c, src, q1, 1.0);
    P.set(c, c, q1, q2, p12);
    P.set(c, c, q1, q3, p13);
    P.set(c, c, q1, snk, p1s);
    P.set(c, c, q2, q1, p21);
    P.set(c, c, q2, q3, p23);
    P.set(c, c, q2, snk, p2s);
    P.set(c, c, q3, snk, p3s);
    m.link(P);

    note("=== Jackson Network ===\n");
    const Sn& sn = m.get_struct();
    section("AG");
    run_ag_inap(sn);
    section("MVA");
    run_mva(sn);
}

// ---------------------------------------------------------------------------
// ag_gnetwork
// ---------------------------------------------------------------------------

/**
 * A G-network: a second class of NEGATIVE customers that removes jobs instead
 * of joining a queue.
 *
 * A signal is an ordinary open class carrying `set_signal`, which is where the
 * removal semantics live; the class still needs a service distribution at every
 * station it passes, exactly as the reference gives it one.
 */
void ag_gnetwork() {
    const double lam_pos = 1.0, lam_neg = 0.3, mu1 = 2.0, mu2 = 3.0;

    Net m("GNetwork-Example");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");

    OpenClass pos(m, "Positive");
    src.set_arrival(pos, Exp(lam_pos));
    q1.set_service(pos, Exp(mu1));
    q2.set_service(pos, Exp(mu2));

    OpenClass neg(m, "Negative");
    m.set_signal(neg, lang::SignalType::NEGATIVE);
    src.set_arrival(neg, Exp(lam_neg));
    q1.set_service(neg, Exp(mu1));
    q2.set_service(neg, Exp(mu2));

    Routing P;
    serial(P, pos, {src, q1, q2, snk});
    serial(P, neg, {src, q1, q2, snk});
    m.link(P);

    note("=== G-Network with Negative Customers ===\n");
    section("AG");
    run_ag_inap(m.get_struct());
}

LINE_EXAMPLE("advanced/agentModels", ag_closed_network);
LINE_EXAMPLE("advanced/agentModels", ag_multiclass_closed);
LINE_EXAMPLE("advanced/agentModels", ag_tandem_open);
LINE_EXAMPLE("advanced/agentModels", ag_tandem_phasetype);
LINE_EXAMPLE("advanced/agentModels", ag_jackson_network);
LINE_EXAMPLE("advanced/agentModels", ag_gnetwork);

}  // namespace examples
}  // namespace line
