/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/stochPetriNet/`: the ten stochastic Petri nets.
 *
 * NINE OF THE TEN ARE SOLVED BY SolverJMT ALONE, and this port drives the same
 * JMT engine through `jmt_avg`, so each model is built -- which is where the
 * port's own SPN encoding gets exercised -- printed, and then handed to the
 * simulator the reference used. Answering a JMT block with SolverCTMC would
 * report a number the reference never asked for, under a solver name that did
 * not produce it.
 *
 * THE TENTH, `spn_queueing_place`, is a QUEUEING Petri net: its places carry an
 * embedded queue with a scheduling discipline of their own, which JMT cannot
 * express. The reference solves it with SolverLDES and cross-checks it against
 * the equivalent Delay + M/M/1 finite-population network under SolverMVA, and
 * both blocks run here.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "line/api/spn/spn_metrics.h"
#include "line/api/spn/spn_lpbnd.h"
#include "line/api/spn/spn_pf.h"
#include "line/io/line_citations.h"
#include "line/solvers/ba/solver_ba_runner.h"

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/distribution.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/fluid/fluid_petri.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

typedef qn::TransitionParam<double> TransParam;
using lang::TimingStrategy;

/** `Transition.setNumberOfServers(mode, GlobalConstants.MaxInt)`. */
const double kMaxInt = lang::GlobalConstants::MaxInt;

void print_arcs(const char* label, const Matrix<double>& arc, const Sn& sn, double skip) {
    std::string line;
    for (std::size_t q = 0; q < arc.rows(); ++q)
        for (std::size_t r = 0; r < arc.cols(); ++r) {
            if (arc(q, r) == skip) continue;
            if (!line.empty()) line += "  ";
            line += sn.nodes[q].name + " x" + std::to_string(static_cast<long>(arc(q, r)));
            // The class is named only when there is a choice to make, so a
            // single-class net prints exactly what it printed before.
            if (arc.cols() > 1) line += " " + sn.classes[r].name;
        }
    if (!line.empty()) std::printf("    %-11s %s\n", label, line.c_str());
}

/**
 * The net a JMT-only example built, printed instead of solved.
 *
 * It reads the REFRESHED struct, so printing it also runs the whole refresh
 * chain over the Places and Transitions: an example that prints its net has
 * proved the port can encode it, which is the part of the reference script that
 * survives the missing simulator.
 */
void print_spn(Net& m) {
    const Sn& sn = m.get_struct();
    std::printf("MODEL: %s (%zu nodes, %zu classes)\n", sn.name.c_str(), sn.nodes.size(),
                sn.classes.size());
    for (std::size_t c = 0; c < sn.classes.size(); ++c) {
        const qn::JobClass& cl = sn.classes[c];
        std::printf("  %-12s %-7s population=%-8g priority=%d\n", cl.name.c_str(),
                    cl.type == qn::JobClassType::CLOSED ? "Closed" : "Open", cl.population,
                    cl.prio);
    }
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        if (nd.nodetype == lang::NodeType::Place) {
            std::string marking;
            const typename std::map<std::size_t, std::vector<double> >::const_iterator mk =
                sn.initmarking.find(i + 1);
            if (mk != sn.initmarking.end())
                for (std::size_t c = 0; c < mk->second.size(); ++c)
                    marking += (c ? " " : "") + sn.classes[c].name + "=" +
                               std::to_string(static_cast<long>(mk->second[c]));
            std::printf("  Place      %-12s marking: %s\n", nd.name.c_str(),
                        marking.empty() ? "(from the reference station)" : marking.c_str());
            continue;
        }
        if (nd.nodetype != lang::NodeType::Transition) continue;
        const typename std::map<std::size_t, TransParam>::const_iterator it =
            sn.transparam.find(i + 1);
        if (it == sn.transparam.end()) continue;
        const TransParam& tp = it->second;
        std::printf("  Transition %-12s %zu mode(s)\n", nd.name.c_str(), tp.nmodes);
        for (std::size_t md = 0; md < tp.nmodes; ++md) {
            char servers[32], proc[64];
            if (tp.nmodeservers[md] >= kMaxInt) std::snprintf(servers, sizeof servers, "MaxInt");
            else std::snprintf(servers, sizeof servers, "%g", tp.nmodeservers[md]);
            const bool immediate = tp.timing[md] == TimingStrategy::IMMEDIATE;
            if (immediate || tp.firingproc[md].disabled) proc[0] = '\0';
            else std::snprintf(proc, sizeof proc, "%s mean=%g",
                               lang::process_to_text(tp.firingproc[md].type),
                               tp.firingproc[md].mean);
            std::printf("    %-11s %-9s servers=%-8s weight=%-6g prio=%-4g %s\n",
                        tp.modenames[md].c_str(), immediate ? "IMMEDIATE" : "TIMED", servers,
                        tp.fireweight[md], tp.firingprio[md], proc);
            print_arcs("enabling:", tp.enabling[md], sn, 0.0);
            print_arcs("inhibiting:", tp.inhibiting[md], sn, std::numeric_limits<double>::infinity());
            print_arcs("firing:", tp.firing[md], sn, 0.0);
        }
    }
}

mva::AvgResult<double> mva_run(Net& m) {
    mva::MvaOptions opt;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

/** P1 -> T1 -> P1, one token, and T1 racing an Exp, an Erlang and a HyperExp mode. */
Net spn_basic_closed_model() {
    Net m("model");
    Place p1(m, "P1");
    ClosedClass c1(m, "Class1", 1.0, p1);

    Transition t1(m, "T1");
    const char* names[3] = {"Mode1", "Mode2", "Mode3"};
    const D procs[3] = {D::exp_mean(1.0), lang::erlang_fit_mean_order<double>(1.0, 2),
                        lang::hyperexp_fit_mean_scv<double>(1.0, 4.0)};
    for (int k = 0; k < 3; ++k) {
        const std::size_t md = t1.add_mode(names[k]);
        t1.set_distribution(md, procs[k]);
        t1.set_enabling_conditions(md, c1, p1, 1);
        t1.set_firing_outcome(md, c1, p1, 1);
    }

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, p1, 1.0);
    m.link(R);
    p1.set_initial_marking(std::vector<double>(1, 1.0));
    return m;
}

/** Source -> P1 -> T1 -> Sink, T1 an infinite-server exponential transition. */
Net spn_basic_open_model() {
    Net m("model");
    Source source(m, "Source");
    Sink sink(m, "Sink");
    Place p1(m, "P1");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Exp(1.0));

    Transition t1(m, "T1");
    const std::size_t md = t1.add_mode("Mode1");
    t1.set_number_of_servers(md, kMaxInt);
    t1.set_distribution(md, Exp(4.0));
    t1.set_enabling_conditions(md, c1, p1, 1);
    t1.set_firing_outcome(md, c1, sink, 1);

    Routing R;
    R.set(c1, c1, source, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, sink, 1.0);
    m.link(R);
    return m;
}

/** P1 -> P2 -> P3 -> P4 -> P1 in batches of two, one distribution family each. */
Net spn_closed_fourplaces_model() {
    Net m("model");
    std::vector<Place> p;
    for (int i = 0; i < 4; ++i) p.push_back(Place(m, "P" + std::to_string(i + 1)));
    ClosedClass c1(m, "Class1", 2.0, p[0], 0);

    std::vector<double> mu0, phi0;
    mu0.push_back(1.0);
    mu0.push_back(2.0);
    phi0.push_back(0.6);
    phi0.push_back(1.0);
    const D procs[4] = {Exp(2.0), Erlang(3.0, 4), HyperExp(0.7, 3.0, 1.5),
                        Coxian(mu0, phi0)};

    std::vector<Transition> t;
    for (int i = 0; i < 4; ++i) {
        t.push_back(Transition(m, "T" + std::to_string(i + 1)));
        const std::size_t md = t[i].add_mode("Mode" + std::to_string(i + 1));
        t[i].set_distribution(md, procs[i]);
        t[i].set_enabling_conditions(md, c1, p[i], 2);
        t[i].set_firing_outcome(md, c1, p[(i + 1) % 4], 2);
    }

    Routing R;
    for (int i = 0; i < 4; ++i) {
        R.set(c1, c1, p[i], t[i], 1.0);
        R.set(c1, c1, t[i], p[(i + 1) % 4], 1.0);
    }
    m.link(R);
    for (int i = 0; i < 4; ++i)
        p[i].set_initial_marking(std::vector<double>(1, i == 0 ? 2.0 : 0.0));
    return m;
}

/** P1 branches to P2 (2 tokens) or P3 (3 tokens); both return to P1. */
Net spn_fourmodes_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    Place p3(m, "P3");
    ClosedClass c1(m, "Class1", 8.0, p1, 0);

    std::vector<Place> pl;
    pl.push_back(p1);
    pl.push_back(p2);
    pl.push_back(p3);
    const int from[4] = {0, 0, 1, 2};
    const int to[4] = {1, 2, 0, 0};
    const double count[4] = {2.0, 3.0, 1.0, 2.0};
    const double rate[4] = {2.0, 1.0, 4.0, 2.0};
    std::vector<Transition> t;
    for (int i = 0; i < 4; ++i) {
        t.push_back(Transition(m, "T" + std::to_string(i + 1)));
        const std::size_t md = t[i].add_mode("Mode" + std::to_string(i + 1));
        t[i].set_distribution(md, Exp(rate[i]));
        t[i].set_enabling_conditions(md, c1, pl[from[i]], count[i]);
        t[i].set_firing_outcome(md, c1, pl[to[i]], count[i]);
    }

    Routing R;
    for (int i = 0; i < 4; ++i) {
        R.set(c1, c1, pl[from[i]], t[i], 1.0);
        R.set(c1, c1, t[i], pl[to[i]], 1.0);
    }
    m.link(R);
    p1.set_initial_marking(std::vector<double>(1, 8.0));
    p2.set_initial_marking(std::vector<double>(1, 0.0));
    p3.set_initial_marking(std::vector<double>(1, 0.0));
    return m;
}

/** T3 returns P3 to P1 only while P2 is empty: the inhibitor arc. */
Net spn_inhibiting_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    Place p3(m, "P3");
    ClosedClass c1(m, "Class1", 4.0, p1, 0);

    Transition t1(m, "T1");
    const std::size_t mode1 = t1.add_mode("Mode1");
    t1.set_distribution(mode1, Exp(2.0));
    t1.set_enabling_conditions(mode1, c1, p1, 2);
    t1.set_firing_outcome(mode1, c1, p2, 2);
    const std::size_t mode2 = t1.add_mode("Mode2");
    t1.set_distribution(mode2, Exp(1.0));
    t1.set_enabling_conditions(mode2, c1, p1, 1);
    t1.set_firing_outcome(mode2, c1, p3, 1);

    Transition t2(m, "T2");
    const std::size_t mode3 = t2.add_mode("Mode3");
    t2.set_distribution(mode3, Exp(4.0));
    t2.set_enabling_conditions(mode3, c1, p2, 1);
    t2.set_firing_outcome(mode3, c1, p1, 1);

    Transition t3(m, "T3");
    const std::size_t mode4 = t3.add_mode("Mode4");
    t3.set_distribution(mode4, Exp(1.0));
    t3.set_enabling_conditions(mode4, c1, p3, 3);
    t3.set_inhibiting_conditions(mode4, c1, p2, 1);
    t3.set_firing_outcome(mode4, c1, p1, 3);

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c1, c1, p2, t3, 1.0);
    R.set(c1, c1, p3, t3, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c1, c1, t1, p3, 1.0);
    R.set(c1, c1, t2, p1, 1.0);
    R.set(c1, c1, t3, p1, 1.0);
    m.link(R);
    p1.set_initial_marking(std::vector<double>(1, 4.0));
    p2.set_initial_marking(std::vector<double>(1, 0.0));
    p3.set_initial_marking(std::vector<double>(1, 0.0));
    return m;
}

/** Seven places and eight transitions, four of them immediate, one inhibited. */
Net spn_open_sevenplaces_model() {
    Net m("model");
    Source source(m, "Source");
    Sink sink(m, "Sink");
    std::vector<Place> p;
    for (int i = 0; i < 7; ++i) p.push_back(Place(m, "P" + std::to_string(i + 1)));
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, D::exp_mean(1.0));

    std::vector<Transition> t;
    std::size_t md = 0;

    t.push_back(Transition(m, "T1"));
    md = t[0].add_mode("Mode1");
    t[0].set_number_of_servers(md, kMaxInt);
    t[0].set_distribution(md, Exp(4.0));
    t[0].set_enabling_conditions(md, c1, p[0], 1);
    t[0].set_firing_outcome(md, c1, p[1], 1);

    t.push_back(Transition(m, "T2"));
    md = t[1].add_mode("Mode1");
    t[1].set_number_of_servers(md, kMaxInt);
    t[1].set_timing_strategy(md, TimingStrategy::IMMEDIATE);
    t[1].set_enabling_conditions(md, c1, p[1], 1);
    t[1].set_firing_outcome(md, c1, p[2], 1);

    t.push_back(Transition(m, "T3"));
    md = t[2].add_mode("Mode1");
    t[2].set_number_of_servers(md, kMaxInt);
    t[2].set_timing_strategy(md, TimingStrategy::IMMEDIATE);
    t[2].set_enabling_conditions(md, c1, p[1], 1);
    t[2].set_firing_outcome(md, c1, p[3], 1);

    t.push_back(Transition(m, "T4"));
    md = t[3].add_mode("Mode1");
    t[3].set_number_of_servers(md, kMaxInt);
    t[3].set_timing_strategy(md, TimingStrategy::IMMEDIATE);
    t[3].set_enabling_conditions(md, c1, p[2], 1);
    t[3].set_enabling_conditions(md, c1, p[4], 1);
    t[3].set_firing_outcome(md, c1, p[4], 1);
    t[3].set_firing_outcome(md, c1, p[5], 1);

    t.push_back(Transition(m, "T5"));
    md = t[4].add_mode("Mode1");
    t[4].set_number_of_servers(md, kMaxInt);
    t[4].set_timing_strategy(md, TimingStrategy::IMMEDIATE);
    t[4].set_enabling_conditions(md, c1, p[3], 1);
    t[4].set_enabling_conditions(md, c1, p[4], 1);
    t[4].set_inhibiting_conditions(md, c1, p[5], 1);
    t[4].set_firing_outcome(md, c1, p[6], 1);

    t.push_back(Transition(m, "T6"));
    md = t[5].add_mode("Mode1");
    t[5].set_number_of_servers(md, kMaxInt);
    t[5].set_distribution(md, Erlang(2.0, 2));
    t[5].set_enabling_conditions(md, c1, p[5], 1);
    t[5].set_firing_outcome(md, c1, p[0], 1);

    t.push_back(Transition(m, "T7"));
    md = t[6].add_mode("Mode1");
    t[6].set_number_of_servers(md, kMaxInt);
    t[6].set_distribution(md, Exp(2.0));
    t[6].set_enabling_conditions(md, c1, p[6], 1);
    t[6].set_firing_outcome(md, c1, p[0], 1);
    t[6].set_firing_outcome(md, c1, p[4], 1);

    t.push_back(Transition(m, "T8"));
    md = t[7].add_mode("Mode1");
    t[7].set_number_of_servers(md, kMaxInt);
    t[7].set_distribution(md, Exp(2.0));
    t[7].set_enabling_conditions(md, c1, p[3], 1);
    t[7].set_firing_outcome(md, c1, sink, 1);

    Routing R;
    R.set(c1, c1, source, p[0], 1.0);
    R.set(c1, c1, p[0], t[0], 1.0);
    R.set(c1, c1, p[1], t[1], 1.0);
    R.set(c1, c1, p[1], t[2], 1.0);
    R.set(c1, c1, p[2], t[3], 1.0);
    R.set(c1, c1, p[3], t[4], 1.0);
    R.set(c1, c1, p[4], t[3], 1.0);
    R.set(c1, c1, p[4], t[4], 1.0);
    R.set(c1, c1, p[5], t[4], 1.0);
    R.set(c1, c1, p[5], t[5], 1.0);
    R.set(c1, c1, p[6], t[6], 1.0);
    R.set(c1, c1, p[3], t[7], 1.0);
    R.set(c1, c1, t[0], p[1], 1.0);
    R.set(c1, c1, t[1], p[2], 1.0);
    R.set(c1, c1, t[2], p[3], 1.0);
    R.set(c1, c1, t[3], p[4], 1.0);
    R.set(c1, c1, t[3], p[5], 1.0);
    R.set(c1, c1, t[4], p[6], 1.0);
    R.set(c1, c1, t[5], p[0], 1.0);
    R.set(c1, c1, t[6], sink, 1.0);
    R.set(c1, c1, t[6], p[0], 1.0);
    R.set(c1, c1, t[6], p[4], 1.0);
    R.set(c1, c1, t[7], sink, 1.0);
    m.link(R);

    const double init[7] = {2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    for (int i = 0; i < 7; ++i) p[i].set_initial_marking(std::vector<double>(1, init[i]));
    return m;
}

/** Source -> P1 -> T1 -> Sink with a single-server Pareto(3, 1) firing time. */
Net spn_pareto_service_model() {
    Net m("model");
    Source source(m, "Source");
    Sink sink(m, "Sink");
    Place p1(m, "P1");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Exp(0.5));

    Transition t1(m, "T1");
    const std::size_t md = t1.add_mode("Mode1");
    t1.set_number_of_servers(md, 1.0);
    t1.set_distribution(md, Pareto(3.0, 1.0));
    t1.set_enabling_conditions(md, c1, p1, 1);
    t1.set_firing_outcome(md, c1, sink, 1);

    Routing R;
    R.set(c1, c1, source, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, sink, 1.0);
    m.link(R);
    return m;
}

/**
 * A closed QUEUEING Petri net: two queueing places exchanging N tokens.
 *
 * A queueing place embeds a scheduling station in the place itself. A token
 * arriving there is served by that embedded queue and only then reaches the
 * DEPOSITORY the output transitions consume from, so the place's marking counts
 * waiting, in-service and deposited tokens while its output arcs see the last
 * of the three. `Place(model, name, sched)` declares the discipline and
 * `set_service` is what turns the place into a queueing one.
 *
 * CPU is a single-server FCFS place and Think an infinite-server one; the two
 * IMMEDIATE transitions move one token per firing. That makes the net the
 * queueing-Petri-net rendering of a machine-repairman model, which is why
 * `spn_queueing_place_ref` can state its exact answer.
 */
Net spn_queueing_place_model(double N) {
    Net m("QueueingPetriNet");
    Place cpu(m, "CPU", SchedStrategy::FCFS);
    Place think(m, "Think", SchedStrategy::INF);
    ClosedClass jobs(m, "Jobs", N, think, 0);
    cpu.set_service(jobs, Exp(1.5));
    think.set_service(jobs, Exp(0.5));
    Transition to_cpu(m, "toCPU");
    const std::size_t m1 = to_cpu.add_mode("m1");
    to_cpu.set_timing_strategy(m1, TimingStrategy::IMMEDIATE);
    to_cpu.set_enabling_conditions(m1, jobs, think, 1);
    to_cpu.set_firing_outcome(m1, jobs, cpu, 1);

    Transition to_think(m, "toThink");
    const std::size_t m2 = to_think.add_mode("m2");
    to_think.set_timing_strategy(m2, TimingStrategy::IMMEDIATE);
    to_think.set_enabling_conditions(m2, jobs, cpu, 1);
    to_think.set_firing_outcome(m2, jobs, think, 1);

    Routing P;
    P.set(jobs, jobs, think, to_cpu, 1.0);
    P.set(jobs, jobs, to_cpu, cpu, 1.0);
    P.set(jobs, jobs, cpu, to_think, 1.0);
    P.set(jobs, jobs, to_think, think, 1.0);
    m.link(P);
    think.set_initial_marking(std::vector<double>{N});
    cpu.set_initial_marking(std::vector<double>{0.0});
    return m;
}

/**
 * The exact cross-check `spn_queueing_place.py` runs beside the QPN: the
 * equivalent finite-population Delay + M/M/1 network, whose SolverMVA answer is
 * the reference the queueing-Petri-net simulation is validated against.
 */
Net spn_queueing_place_ref(double N) {
    Net m("ref");
    Delay think(m, "Think");
    Queue cpu(m, "CPU", SchedStrategy::FCFS);
    ClosedClass jobs(m, "Jobs", N, think, 0);
    think.set_service(jobs, Exp(0.5));
    cpu.set_service(jobs, Exp(1.5));
    Routing P;
    cyclic(P, jobs, {think, cpu});
    m.link(P);
    return m;
}

/**
 * A COLOURED net: two closed classes whose tokens move along their own arcs.
 *
 * T1 has one mode per class -- two Class1 tokens at a time and one Class2 --
 * while T2 returns Class1 singly under an Erlang law and T3 returns Class2 four
 * at a time. The two colours share both places and never substitute for one
 * another, which is what the per-(place, class) arcs of `TransitionParam` say.
 */
Net spn_closed_twoplaces_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    ClosedClass c1(m, "Class1", 10.0, p1, 0);
    ClosedClass c2(m, "Class2", 7.0, p1, 0);
    Transition t1(m, "T1");
    const std::size_t m1 = t1.add_mode("Mode1");
    t1.set_distribution(m1, Exp(2.0));
    t1.set_enabling_conditions(m1, c1, p1, 2);
    t1.set_firing_outcome(m1, c1, p2, 2);
    const std::size_t m2 = t1.add_mode("Mode2");
    t1.set_distribution(m2, Exp(3.0));
    t1.set_enabling_conditions(m2, c2, p1, 1);
    t1.set_firing_outcome(m2, c2, p2, 1);

    Transition t2(m, "T2");
    const std::size_t m3 = t2.add_mode("Mode3");
    t2.set_distribution(m3, Erlang(1.5, 2));  // rate 1.5 per phase, mean 2/1.5
    t2.set_enabling_conditions(m3, c1, p2, 1);
    t2.set_firing_outcome(m3, c1, p1, 1);

    Transition t3(m, "T3");
    const std::size_t m4 = t3.add_mode("Mode4");
    t3.set_distribution(m4, Exp(0.5));
    t3.set_enabling_conditions(m4, c2, p2, 4);
    t3.set_firing_outcome(m4, c2, p1, 4);

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c2, c2, p1, t1, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c2, c2, p2, t3, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c2, c2, t1, p2, 1.0);
    R.set(c1, c1, t2, p1, 1.0);
    R.set(c2, c2, t3, p1, 1.0);
    m.link(R);
    p1.set_initial_marking(std::vector<double>{10.0, 7.0});
    p2.set_initial_marking(std::vector<double>{0.0, 0.0});
    return m;
}

/**
 * A COLOURED GENERALIZED net: two colours, an immediate transition and a timed
 * one.
 *
 * `admit` is IMMEDIATE and carries one mode per colour, each held off by
 * inhibitor arcs on BOTH colours, so the Server place is a mutual-exclusion
 * resource and the firing weights arbitrate whenever both colours are waiting.
 * `serve` is TIMED and returns the token at a colour-dependent rate. Admission
 * costs no time, so the server is never idle and the result is exact by hand:
 * the weights split the completions 2:1 in favour of Gold, hence
 * X_Gold = 2*X_Silver, U_c = X_c/mu_c and U_Gold + U_Silver = 1.
 */
Net spn_colored_gspn_model() {
    Net m("ColoredGSPN");
    Place buf(m, "Buffer");
    Place srv(m, "Server");
    ClosedClass gold(m, "Gold", 2.0, buf, 0);
    ClosedClass silver(m, "Silver", 2.0, buf, 0);
    Transition admit(m, "admit");
    const std::size_t ag = admit.add_mode("gold");
    admit.set_timing_strategy(ag, TimingStrategy::IMMEDIATE);
    admit.set_firing_weights(ag, 2.0);
    admit.set_enabling_conditions(ag, gold, buf, 1);
    admit.set_inhibiting_conditions(ag, gold, srv, 1);
    admit.set_inhibiting_conditions(ag, silver, srv, 1);
    admit.set_firing_outcome(ag, gold, srv, 1);
    const std::size_t as = admit.add_mode("silver");
    admit.set_timing_strategy(as, TimingStrategy::IMMEDIATE);
    admit.set_firing_weights(as, 1.0);
    admit.set_enabling_conditions(as, silver, buf, 1);
    admit.set_inhibiting_conditions(as, gold, srv, 1);
    admit.set_inhibiting_conditions(as, silver, srv, 1);
    admit.set_firing_outcome(as, silver, srv, 1);

    Transition serve(m, "serve");
    const std::size_t sg = serve.add_mode("gold");
    serve.set_distribution(sg, Exp(3.0));
    serve.set_enabling_conditions(sg, gold, srv, 1);
    serve.set_firing_outcome(sg, gold, buf, 1);
    const std::size_t ss = serve.add_mode("silver");
    serve.set_distribution(ss, Exp(1.5));
    serve.set_enabling_conditions(ss, silver, srv, 1);
    serve.set_firing_outcome(ss, silver, buf, 1);

    Routing R;
    R.set(gold, gold, buf, admit, 1.0);
    R.set(silver, silver, buf, admit, 1.0);
    R.set(gold, gold, admit, srv, 1.0);
    R.set(silver, silver, admit, srv, 1.0);
    R.set(gold, gold, srv, serve, 1.0);
    R.set(silver, silver, srv, serve, 1.0);
    R.set(gold, gold, serve, buf, 1.0);
    R.set(silver, silver, serve, buf, 1.0);
    m.link(R);
    buf.set_initial_marking(std::vector<double>{2.0, 2.0});
    srv.set_initial_marking(std::vector<double>{0.0, 0.0});
    return m;
}

/** P1 -> P2 in batches of four, P2 -> P1 in batches of two. */
Net spn_twomodes_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    ClosedClass c1(m, "Class1", 10.0, p1, 0);

    Transition t1(m, "T1");
    const std::size_t m1 = t1.add_mode("Mode1");
    t1.set_distribution(m1, Exp(2.0));
    t1.set_enabling_conditions(m1, c1, p1, 4);
    t1.set_firing_outcome(m1, c1, p2, 4);

    Transition t2(m, "T2");
    const std::size_t m2 = t2.add_mode("Mode2");
    t2.set_distribution(m2, Exp(3.0));
    t2.set_enabling_conditions(m2, c1, p2, 2);
    t2.set_firing_outcome(m2, c1, p1, 2);

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c1, c1, t2, p1, 1.0);
    m.link(R);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Examples
// ---------------------------------------------------------------------------

void spn_basic_closed() {
    Net m = spn_basic_closed_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void spn_basic_open() {
    Net m = spn_basic_open_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void spn_closed_fourplaces() {
    Net m = spn_closed_fourplaces_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

/**
 * The COLOURED net of this directory: two classes with their own arcs.
 *
 * T1/Mode1 takes 2 Class1 tokens from P1, T1/Mode2 takes 1 Class2 token from the
 * same place, T2 returns Class1 one at a time and T3 returns Class2 four at a
 * time. Until 2026-08-12 `TransitionParam` stored one arc per (mode, node) and
 * this net was refused rather than approximated; the arcs now carry the class,
 * so it is built and simulated like the other nine.
 */
void spn_closed_twoplaces() {
    Net m = spn_closed_twoplaces_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

/**
 * The coloured GSPN, solved exactly and cross-checked by simulation.
 *
 * Unlike the JMT-only nets above this one is answered by SolverCTMC: the
 * immediate transition is eliminated by stochastic complementation, leaving a
 * chain small enough to solve exactly, and SolverSSA re-runs it as a sample
 * path to confirm the same numbers.
 */
void spn_colored_gspn() {
    Net m = spn_colored_gspn_model();
    print_spn(m);
    section("CTMC");
    SolverOpts ctmc;
    ctmc.cutoff = 4;
    print_avg(solve_avg("CTMC", m, ctmc));
    section("SSA");
    SolverOpts ssa;
    ssa.samples = 200000;
    ssa.seed = 23000;
    print_avg(solve_avg("SSA", m, ssa));
}

// ---------------------------------------------------------------------------
// test_spn_nrm_open
// ---------------------------------------------------------------------------

namespace {

/** Source Exp(lambda) -> P1 -> T1 Exp(mu), one server -> Sink: an M/M/1 at P1. */
Net mm1spn(double lambda, double mu) {
    Net m("mm1spn");
    Source source(m, "Source");
    Sink sink(m, "Sink");
    Place p1(m, "P1");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Exp(lambda));

    Transition t1(m, "T1");
    const std::size_t md = t1.add_mode("Mode1");
    t1.set_distribution(md, Exp(mu));
    t1.set_enabling_conditions(md, c1, p1, 1);
    t1.set_firing_outcome(md, c1, sink, 1);

    Routing R;
    R.set(c1, c1, source, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, sink, 1.0);
    m.link(R);
    return m;
}

/** The same, in series: Source -> P1 -> T1 -> P2 -> T2 -> Sink. */
Net tandemspn(double lambda, double mu1, double mu2) {
    Net m("tandemspn");
    Source source(m, "Source");
    Sink sink(m, "Sink");
    Place p1(m, "P1");
    Place p2(m, "P2");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Exp(lambda));

    Transition t1(m, "T1");
    const std::size_t m1 = t1.add_mode("Mode1");
    t1.set_distribution(m1, Exp(mu1));
    t1.set_enabling_conditions(m1, c1, p1, 1);
    t1.set_firing_outcome(m1, c1, p2, 1);

    Transition t2(m, "T2");
    const std::size_t m2 = t2.add_mode("Mode1");
    t2.set_distribution(m2, Exp(mu2));
    t2.set_enabling_conditions(m2, c1, p2, 1);
    t2.set_firing_outcome(m2, c1, sink, 1);

    Routing R;
    R.set(c1, c1, source, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c1, c1, t2, sink, 1.0);
    m.link(R);
    return m;
}

/** The relative error the reference asserts on, and the assertion itself. */
void check_close(const std::string& what, double got, double want, double rtol) {
    const double rel = std::fabs(got - want) / std::fabs(want);
    std::printf("  %-28s %10.5f  vs %10.5f   rel %.4f  %s\n", what.c_str(), got, want, rel,
                rel < rtol ? "OK" : "FAIL");
    if (rel >= rtol)
        throw NumericError(what + ": " + std::to_string(got) + " vs " + std::to_string(want));
}

}  // namespace

/**
 * The SSA Next-Reaction-Method path on OPEN Petri nets.
 *
 * A Source feeds a Place whose tokens drain through a Transition to a Sink.
 * Before the Source-arrival reaction was added the fed Place stayed empty and
 * the run threw "Deadlock: no transition is enabled". Each net asserts that the
 * solver really ran 'nrm', and that the simulated marking mean and throughput
 * match the analytic M/M/1 result at the Place: mean tokens = rho/(1-rho) and
 * throughput = lambda.
 */
void test_spn_nrm_open() {
    const double RTOL = 0.04;
    const std::size_t SAMPLES = 300000;
    SolverOpts ssa;
    ssa.method = "nrm";
    ssa.samples = SAMPLES;
    ssa.seed = 23000;

    // Net 1: M/M/1 SPN, Source Exp(0.5) -> P1 -> T1 Exp(1.0) -> Sink.
    {
        const double lambda = 0.5, mu = 1.0;
        const double rho = lambda / mu, qExact = rho / (1.0 - rho);
        Net m = mm1spn(lambda, mu);
        const AvgTable t = solve_avg("SSA", m, ssa);
        std::printf("\n--- Net 1: M/M/1 SPN, lambda=%g mu=%g ---\n", lambda, mu);
        print_avg(t);
        if (t.method.find("nrm") == std::string::npos)
            throw NumericError("net1 did not run NRM, it ran '" + t.method + "'");
        check_close("net1 P1 tokens", t.get("QLen", "P1", "Class1"), qExact, RTOL);
        check_close("net1 Source tput", t.get("Tput", "Source", "Class1"), lambda, RTOL);
        check_close("net1 P1 tput", t.get("Tput", "P1", "Class1"), lambda, RTOL);

        // The reference's second solver, and it runs here: the SPN half of
        // `io/jmt_writer.h` exports places, transitions, modes, enabling,
        // inhibiting and firing, so this net reaches JSIM. The block used to
        // say the port had no JMT, which stopped being true.
        Net mj = mm1spn(lambda, mu);
        const AvgTable tj = solve_avg("JMT", mj, sim_opts(23000, SAMPLES));
        section("JMT");
        print_avg(tj);
        check_close("net1 P1 tokens vs JMT", t.get("QLen", "P1", "Class1"),
                    tj.get("QLen", "P1", "Class1"), RTOL);
    }

    // Net 2: open tandem, two places in series.
    {
        const double lambda = 0.5, mu1 = 1.0, mu2 = 2.0;
        const double q1 = (lambda / mu1) / (1.0 - lambda / mu1);
        const double q2 = (lambda / mu2) / (1.0 - lambda / mu2);
        Net m = tandemspn(lambda, mu1, mu2);
        const AvgTable t = solve_avg("SSA", m, ssa);
        std::printf("\n--- Net 2: open tandem, lambda=%g mu1=%g mu2=%g ---\n", lambda, mu1, mu2);
        print_avg(t);
        if (t.method.find("nrm") == std::string::npos)
            throw NumericError("net2 did not run NRM, it ran '" + t.method + "'");
        check_close("net2 P1 tokens", t.get("QLen", "P1", "Class1"), q1, RTOL);
        check_close("net2 P2 tokens", t.get("QLen", "P2", "Class1"), q2, RTOL);
        check_close("net2 P1 tput", t.get("Tput", "P1", "Class1"), lambda, RTOL);
        check_close("net2 P2 tput", t.get("Tput", "P2", "Class1"), lambda, RTOL);
    }

    note("test_spn_nrm_open passed");
}

void spn_fourmodes() {
    Net m = spn_fourmodes_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void spn_inhibiting() {
    Net m = spn_inhibiting_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void spn_open_sevenplaces() {
    Net m = spn_open_sevenplaces_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void spn_pareto_service() {
    Net m = spn_pareto_service_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
}

/**
 * The queueing Petri net and its exact cross-check.
 *
 * LDES is the solver here because a queueing place is a simulated object: the
 * engine holds the embedded queue, serves it under the place's own discipline,
 * and lets the output transitions see the depository alone. The second block is
 * the equivalent Delay + M/M/1 finite-population network under SolverMVA, whose
 * exact answer the simulated one is read against.
 */
void spn_queueing_place() {
    const double N = 4.0;
    Net m = spn_queueing_place_model(N);
    print_spn(m);
    section("LDES");
    print_avg(m.get_struct(), ldes_avg(m, sim_opts(23000, 200000)));
    Net ref = spn_queueing_place_ref(N);
    section("MVA");
    print_avg(ref.get_struct(), mva_run(ref));
}

void spn_twomodes() {
    Net m = spn_twomodes_model();
    print_spn(m);
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

// ---------------------------------------------------------------------------
// spn_productform_nc
// ---------------------------------------------------------------------------

namespace {

/** P0 -> T0 -> P1 -> T1 -> P2 -> T2 -> P0, one class, single-server modes. */
Net pf_cyclic_model(double ntokens) {
    const double rates[3] = {1.0, 1.5, 2.0};
    Net m("spn");
    std::vector<Place> p;
    for (int i = 0; i < 3; ++i) p.push_back(Place(m, "P" + std::to_string(i)));
    ClosedClass c1(m, "Class1", ntokens, p[0], 0);
    std::vector<Transition> t;
    for (int i = 0; i < 3; ++i) {
        t.push_back(Transition(m, "T" + std::to_string(i)));
        const std::size_t md = t[i].add_mode("fire");
        t[i].set_distribution(md, Exp(rates[i]));
        t[i].set_enabling_conditions(md, c1, p[i], 1);
        t[i].set_firing_outcome(md, c1, p[(i + 1) % 3], 1);
    }
    Routing R;
    for (int i = 0; i < 3; ++i) {
        R.set(c1, c1, p[i], t[i], 1.0);
        R.set(c1, c1, t[i], p[(i + 1) % 3], 1.0);
    }
    m.link(R);
    for (int i = 0; i < 3; ++i)
        p[i].set_initial_marking(std::vector<double>(1, i == 0 ? ntokens : 0.0));
    return m;
}

/**
 * P0 -(Tf)-> P1 + P2 -(Tj)-> P3 -(Tb)-> P0.
 *
 * Tf consumes ONE token and produces TWO, Tj the reverse, so the marking is not
 * a conserved job population and the net has no queueing-network counterpart.
 * Its place invariant is 2*m0 + m1 + m2 + 2*m3.
 */
Net pf_forkjoin_model(double ntokens) {
    Net m("fj");
    std::vector<Place> p;
    for (int i = 0; i < 4; ++i) p.push_back(Place(m, "P" + std::to_string(i)));
    ClosedClass c1(m, "C", ntokens, p[0], 0);

    Transition tf(m, "Tf");
    const std::size_t mf = tf.add_mode("f");
    tf.set_distribution(mf, Exp(1.3));
    tf.set_enabling_conditions(mf, c1, p[0], 1);
    tf.set_firing_outcome(mf, c1, p[1], 1);
    tf.set_firing_outcome(mf, c1, p[2], 1);

    Transition tj(m, "Tj");
    const std::size_t mj = tj.add_mode("j");
    tj.set_distribution(mj, Exp(0.7));
    tj.set_enabling_conditions(mj, c1, p[1], 1);
    tj.set_enabling_conditions(mj, c1, p[2], 1);
    tj.set_firing_outcome(mj, c1, p[3], 1);

    Transition tb(m, "Tb");
    const std::size_t mb = tb.add_mode("b");
    tb.set_distribution(mb, Exp(1.9));
    tb.set_enabling_conditions(mb, c1, p[3], 1);
    tb.set_firing_outcome(mb, c1, p[0], 1);

    Routing R;
    R.set(c1, c1, p[0], tf, 1.0);
    R.set(c1, c1, tf, p[1], 1.0);
    R.set(c1, c1, tf, p[2], 1.0);
    R.set(c1, c1, p[1], tj, 1.0);
    R.set(c1, c1, p[2], tj, 1.0);
    R.set(c1, c1, tj, p[3], 1.0);
    R.set(c1, c1, p[3], tb, 1.0);
    R.set(c1, c1, tb, p[0], 1.0);
    m.link(R);
    for (int i = 0; i < 4; ++i)
        p[i].set_initial_marking(std::vector<double>(1, i == 0 ? ntokens : 0.0));
    return m;
}

/** Print what spn_pf derived and what the measures come to. */
void pf_report(Net& m, const char* label) {
    section(label);
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    const spn::SpnMetrics<double> met = spn::spn_metrics<double>(pf.spn.mdds, pf.g, pf.spn.info);
    std::printf("  product form: %s, %d complexes, %zu linkage classes, rank %zu, deficiency %d, "
                "%s\n",
                pf.kind.c_str(), static_cast<int>(pf.complexes.size()), pf.linkage, pf.srank,
                pf.deficiency, pf.weakly_reversible ? "weakly reversible" : "not weakly reversible");
    std::printf("  G = %.9f (complex-balance residual %.2e)\n", met.G, pf.residual);
    for (std::size_t l = 0; l < met.tokens.size(); ++l)
        std::printf("  %-8s tokens %10.6f  util %10.6f  tput %10.6f\n",
                    pf.spn.info.placenames[l].c_str(), met.tokens[l], met.place_util[l],
                    met.place_tput[l]);
    for (std::size_t e = 0; e < met.mode_tput.size(); ++e)
        std::printf("  mode %zu throughput %10.6f\n", e, met.mode_tput[e]);
}

}  // namespace

/**
 * Solve a stochastic Petri net analytically, the way SolverNC's `rec` method
 * does: `spn_pf` derives the product form by complex balance
 * (Coleman-Henderson-Taylor, Perform. Eval. 26(3), 1996), `mdd_rec` evaluates
 * the normalising constant by one memoised walk of the decision diagram holding
 * the reachable set (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490), and
 * `spn_metrics` reads the measures off masked walks of the same diagram.
 *
 * The second net FORKS: Tf consumes one token and produces two, so the marking
 * is not a conserved job population. That is the case the MDD-rec paper opens
 * with, and no queueing network expresses it.
 */
void spn_productform_nc() {
    Net cyc = pf_cyclic_model(4.0);
    print_spn(cyc);
    pf_report(cyc, "NC (rec)");

    Net fj = pf_forkjoin_model(3.0);
    print_spn(fj);
    pf_report(fj, "NC (rec), fork-join");
}

// ---------------------------------------------------------------------------
// spn_fluid_dae
// ---------------------------------------------------------------------------

/**
 * `P1 <-> P2`, every mode infinite-server with one input arc.
 *
 * The drift is then LINEAR, so the fluid mean is the EXACT mean and the fluid
 * covariance the exact covariance: the marking is Binomial(4, 3/5), whose
 * variance is 4*0.6*0.4 = 0.96.
 */
Net spn_fluid_exact_model() {
    Net m("spn_fluid_exact");
    Place p1(m, "P1");
    Place p2(m, "P2");
    ClosedClass c1(m, "Class1", 4.0, p1, 0);

    Transition t1(m, "T1");
    const std::size_t md1 = t1.add_mode("Mode1");
    t1.set_number_of_servers(md1, kMaxInt);
    t1.set_distribution(md1, Exp(2.0));
    t1.set_enabling_conditions(md1, c1, p1, 1);
    t1.set_firing_outcome(md1, c1, p2, 1);

    Transition t2(m, "T2");
    const std::size_t md2 = t2.add_mode("Mode2");
    t2.set_number_of_servers(md2, kMaxInt);
    t2.set_distribution(md2, Exp(3.0));
    t2.set_enabling_conditions(md2, c1, p2, 1);
    t2.set_firing_outcome(md2, c1, p1, 1);

    Routing R;
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, p2, 1.0);
    R.set(c1, c1, p2, t2, 1.0);
    R.set(c1, c1, t2, p1, 1.0);
    m.link(R);
    p1.set_initial_marking({4.0});
    p2.set_initial_marking({0.0});
    return m;
}

/**
 * `P1 -T1-> P2 -(immediate)-> P3 -T3-> P1`.
 *
 * The vanishing place P2 holds exactly zero mass, pinned by the algebraic
 * constraint rather than approached by a large finite firing rate, so the net
 * answers as the reduced two-place net below does.
 */
Net spn_fluid_immediate_model() {
    Net m("spn_fluid_immediate");
    Place q1(m, "P1");
    Place q2(m, "P2");
    Place q3(m, "P3");
    ClosedClass c1(m, "Class1", 4.0, q1, 0);

    Transition u1(m, "T1");
    const std::size_t a1 = u1.add_mode("M1");
    u1.set_distribution(a1, Exp(2.0));
    u1.set_enabling_conditions(a1, c1, q1, 1);
    u1.set_firing_outcome(a1, c1, q2, 1);

    Transition ui(m, "Ti");
    const std::size_t ai = ui.add_mode("Mi");
    ui.set_timing_strategy(ai, TimingStrategy::IMMEDIATE);
    ui.set_distribution(ai, Immediate());
    ui.set_enabling_conditions(ai, c1, q2, 1);
    ui.set_firing_outcome(ai, c1, q3, 1);

    Transition u3(m, "T3");
    const std::size_t a3 = u3.add_mode("M3");
    u3.set_distribution(a3, Exp(3.0));
    u3.set_enabling_conditions(a3, c1, q3, 1);
    u3.set_firing_outcome(a3, c1, q1, 1);

    Routing R;
    R.set(c1, c1, q1, u1, 1.0);
    R.set(c1, c1, u1, q2, 1.0);
    R.set(c1, c1, q2, ui, 1.0);
    R.set(c1, c1, ui, q3, 1.0);
    R.set(c1, c1, q3, u3, 1.0);
    R.set(c1, c1, u3, q1, 1.0);
    m.link(R);
    q1.set_initial_marking({4.0});
    q2.set_initial_marking({0.0});
    q3.set_initial_marking({0.0});
    return m;
}

/** The same net with the immediate transition eliminated by hand. */
Net spn_fluid_reduced_model() {
    Net m("spn_fluid_reduced");
    Place w1(m, "P1");
    Place w3(m, "P3");
    ClosedClass c1(m, "Class1", 4.0, w1, 0);

    Transition v1(m, "T1");
    const std::size_t b1 = v1.add_mode("M1");
    v1.set_distribution(b1, Exp(2.0));
    v1.set_enabling_conditions(b1, c1, w1, 1);
    v1.set_firing_outcome(b1, c1, w3, 1);

    Transition v3(m, "T3");
    const std::size_t b3 = v3.add_mode("M3");
    v3.set_distribution(b3, Exp(3.0));
    v3.set_enabling_conditions(b3, c1, w3, 1);
    v3.set_firing_outcome(b3, c1, w1, 1);

    Routing R;
    R.set(c1, c1, w1, v1, 1.0);
    R.set(c1, c1, v1, w3, 1.0);
    R.set(c1, c1, w3, v3, 1.0);
    R.set(c1, c1, v3, w1, 1.0);
    m.link(R);
    w1.set_initial_marking({4.0});
    w3.set_initial_marking({0.0});
    return m;
}

/**
 * Fluid (mean-field) analysis of a stochastic Petri net, with SolverFLD.
 *
 * A GSPN is a density-dependent Markov population process: the marking is the
 * population, a transition mode is a reaction, and the firing rate
 * lambda*min(enabling degree, servers) is the same min() non-linearity the
 * min-normal closure of SolverFLD exists to smooth. The 'dae' method is the one
 * that can carry it, because a Petri net needs three things stated as EQUATIONS
 * rather than integrated: the P-invariants, which hold to solver tolerance
 * instead of integrator tolerance and supply the rank the drift Jacobian is
 * missing; the firing FLOW of an immediate transition, an algebraic unknown
 * pinned by the constraint that its input place holds no mass; and a bounded
 * place, a linear inequality on the marking.
 *
 * `solver_fluid_run_analyzer` resolves to 'dae' on any model holding a
 * Transition node, so no method has to be named. Unlike every other solver of a
 * Petri net in LINE it also returns a SECOND MOMENT: the marking covariance of
 * the linear noise approximation, which is reached here through
 * `fluid::petri::solver_fluid_petri` because the `PetriReport` carrying it -- the marking
 * variance, the invariants and their residuals -- is wider than the
 * `FluidMomentReport` the station-table wrapper keeps.
 */
void spn_fluid_dae() {
    note("-- a closed net whose fluid answer is EXACT --");
    Net exact = spn_fluid_exact_model();
    const Sn& sn_exact = exact.get_struct();
    section("FLD");
    print_avg_sim(sn_exact, fluid::solver_fluid_run_analyzer(sn_exact, fluid::FluidOptions()));
    ctmc::CtmcOptions copt;
    copt.cutoff = 6.0;
    section("CTMC (cutoff = 6)");
    print_avg(sn_exact, ctmc::solver_ctmc_run_analyzer(sn_exact, copt));

    const fluid::petri::PetriSolution ps =
        fluid::petri::solver_fluid_petri(sn_exact, fluid::petri::PetriOptions());
    std::printf("marking variance: %.6f %.6f  (exact 0.96)\n", ps.petri.marking_var(0, 0),
                ps.petri.marking_var(1, 0));
    std::printf("invariant \"%s\" = %g, error %.2e\n", ps.petri.invariant_label[0].c_str(),
                ps.petri.invariant_value[0], ps.petri.invariant_error[0]);

    note("\n-- an immediate transition, as an algebraic flow --");
    Net imm = spn_fluid_immediate_model();
    const Sn& sn_imm = imm.get_struct();
    section("FLD");
    print_avg_sim(sn_imm, fluid::solver_fluid_run_analyzer(sn_imm, fluid::FluidOptions()));

    note("\n-- the same net with the immediate transition eliminated by hand --");
    Net red = spn_fluid_reduced_model();
    const Sn& sn_red = red.get_struct();
    section("FLD");
    print_avg_sim(sn_red, fluid::solver_fluid_run_analyzer(sn_red, fluid::FluidOptions()));
}

// ---------------------------------------------------------------------------
// spn_lpbounds
// ---------------------------------------------------------------------------

/**
 * Fig. 2b of Liu (1998): four servers in a line, blocking before service.
 *
 * Server i cannot start until the downstream buffer has a free slot. The
 * buffers hold 3, 2 and 4, and each is a conserved pair of places --
 * (p5,p2), (p4,p1), (p3,p0) -- so the net is a strongly connected marked graph
 * and all four transitions carry the same throughput.
 */
Net spn_lpbounds_prodline_model(const double mu[4]) {
    Net m("liu98");
    Place p5(m, "p5"), p4(m, "p4"), p3(m, "p3");
    Place p2(m, "p2"), p1(m, "p1"), p0(m, "p0");
    ClosedClass c1(m, "Class1", 9.0, p2, 0);

    Transition t1(m, "t1");
    const std::size_t m1 = t1.add_mode("m1");
    t1.set_distribution(m1, Exp(mu[0]));
    t1.set_enabling_conditions(m1, c1, p2, 1);
    t1.set_firing_outcome(m1, c1, p5, 1);

    Transition t2(m, "t2");
    const std::size_t m2 = t2.add_mode("m2");
    t2.set_distribution(m2, Exp(mu[1]));
    t2.set_enabling_conditions(m2, c1, p5, 1);
    t2.set_enabling_conditions(m2, c1, p1, 1);
    t2.set_firing_outcome(m2, c1, p4, 1);
    t2.set_firing_outcome(m2, c1, p2, 1);

    Transition t3(m, "t3");
    const std::size_t m3 = t3.add_mode("m3");
    t3.set_distribution(m3, Exp(mu[2]));
    t3.set_enabling_conditions(m3, c1, p4, 1);
    t3.set_enabling_conditions(m3, c1, p0, 1);
    t3.set_firing_outcome(m3, c1, p3, 1);
    t3.set_firing_outcome(m3, c1, p1, 1);

    Transition t4(m, "t4");
    const std::size_t m4 = t4.add_mode("m4");
    t4.set_distribution(m4, Exp(mu[3]));
    t4.set_enabling_conditions(m4, c1, p3, 1);
    t4.set_firing_outcome(m4, c1, p0, 1);

    Routing R;
    R.set(c1, c1, p2, t1, 1.0);
    R.set(c1, c1, t1, p5, 1.0);
    R.set(c1, c1, p5, t2, 1.0);
    R.set(c1, c1, p1, t2, 1.0);
    R.set(c1, c1, t2, p4, 1.0);
    R.set(c1, c1, t2, p2, 1.0);
    R.set(c1, c1, p4, t3, 1.0);
    R.set(c1, c1, p0, t3, 1.0);
    R.set(c1, c1, t3, p3, 1.0);
    R.set(c1, c1, t3, p1, 1.0);
    R.set(c1, c1, p3, t4, 1.0);
    R.set(c1, c1, t4, p0, 1.0);
    m.link(R);
    p5.set_initial_marking({0.0});
    p4.set_initial_marking({0.0});
    p3.set_initial_marking({0.0});
    p2.set_initial_marking({3.0});
    p1.set_initial_marking({2.0});
    p0.set_initial_marking({4.0});
    return m;
}

/** Three places, four modes, one inhibitor arc; the token count is conserved. */
Net spn_lpbounds_inhibiting_model(double n) {
    Net m("spn");
    Place P1(m, "P1"), P2(m, "P2"), P3(m, "P3");
    ClosedClass c1(m, "Class1", n, P1, 0);

    Transition T1(m, "T1");
    const std::size_t md1 = T1.add_mode("Mode1");
    T1.set_distribution(md1, Exp(2.0));
    T1.set_enabling_conditions(md1, c1, P1, 2);
    T1.set_firing_outcome(md1, c1, P2, 2);
    const std::size_t md2 = T1.add_mode("Mode2");
    T1.set_distribution(md2, Exp(1.0));
    T1.set_enabling_conditions(md2, c1, P1, 1);
    T1.set_firing_outcome(md2, c1, P3, 1);

    Transition T2(m, "T2");
    const std::size_t md3 = T2.add_mode("Mode3");
    T2.set_distribution(md3, Exp(4.0));
    T2.set_enabling_conditions(md3, c1, P2, 1);
    T2.set_firing_outcome(md3, c1, P1, 1);

    Transition T3(m, "T3");
    const std::size_t md4 = T3.add_mode("Mode4");
    T3.set_distribution(md4, Exp(1.0));
    T3.set_enabling_conditions(md4, c1, P3, 3);
    T3.set_inhibiting_conditions(md4, c1, P2, 1);
    T3.set_firing_outcome(md4, c1, P1, 3);

    Routing R;
    R.set(c1, c1, P1, T1, 1.0);
    R.set(c1, c1, P2, T2, 1.0);
    R.set(c1, c1, P2, T3, 1.0);
    R.set(c1, c1, P3, T3, 1.0);
    R.set(c1, c1, T1, P2, 1.0);
    R.set(c1, c1, T1, P3, 1.0);
    R.set(c1, c1, T2, P1, 1.0);
    R.set(c1, c1, T3, P1, 1.0);
    m.link(R);
    P1.set_initial_marking({n});
    P2.set_initial_marking({0.0});
    P3.set_initial_marking({0.0});
    return m;
}

/**
 * Bound a stochastic Petri net by linear programming.
 *
 * SolverBA's 'spnlp' family is the first BOUNDING route LINE offers for a Petri
 * net. SolverCTMC builds the explicit generator, SolverSSA and SolverLDES
 * simulate, SolverFLD fluidises, and SolverNC 'rec' needs a product form; this
 * one needs none of that. It relaxes the stationary chain to a MOMENT POLYTOPE
 * -- the uniformized evolution equation written for E[X_p], E[X_p^2] and
 * E[X_p1 X_p2], plus behavioural and probabilistic inequalities -- and then
 * minimises and maximises each reported measure over it. Every stationary point
 * of the true chain satisfies every row, so the two optima bracket the exact
 * value.
 *
 * Table 2 of the paper reports four bound columns on five rate vectors, all
 * four reproduced below. The one column not reproduced is its u.b.1, the upper
 * side further tightened by the subnet-throughput theorems (its Thms 1 and 2),
 * which are not implemented; u.b.2 is the column to compare against.
 *
 * Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
 * Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
 * 1014-1030.
 */
void spn_lpbounds() {
    const double mus[5][4] = {{1.0, 1.25, 2.0, 0.5},
                              {1.0, 1.25, 2.0, 2.5},
                              {1.0, 1.25, 1.25, 2.5},
                              {1.0, 1.25, 1.25, 1.0},
                              {1.111, 1.111, 1.111, 1.111}};
    const double pub[5][3] = {{1.165, 1.951, 2.000},
                              {1.829, 2.978, 3.529},
                              {1.581, 2.873, 3.333},
                              {1.359, 2.757, 3.333},
                              {1.350, 2.667, 2.963}};

    std::printf("\nLiu (1998) Table 2: total throughput of the production line\n");
    std::printf("%-5s %-32s %-19s %-19s\n", "case", "Markovian LP", "published",
                "operational LP");
    std::printf("%-5s %9s %9s %9s   %9s %9s   %9s %9s\n", "", "lower", "simul", "upper", "l.b.",
                "u.b.2", "o.l.b.", "o.u.b.");
    for (std::size_t c = 0; c < 5; ++c) {
        Net m = spn_lpbounds_prodline_model(mus[c]);
        const Sn& sn = m.get_struct();
        // The liveness rows of the reference's Table 1 are OPT-IN, because they
        // hold only on a live net and `spn_lpbnd` cannot certify liveness. This
        // one is live: a strongly connected marked graph with a token on every
        // cycle. They are the whole of the lower side, so the published l.b.
        // needs them.
        // A NetworkStruct carries no per-place state, so `spn_lpbnd` takes the
        // initial marking from the reference station of the closed class unless
        // `init` says otherwise. This net starts its 9 tokens spread over the
        // three free-slot places, so it must be passed: places are levelled in
        // node order, p5 p4 p3 p2 p1 p0.
        const std::vector<double> m0 = {0.0, 0.0, 0.0, 3.0, 2.0, 4.0};
        spn::SpnLpOptions lo;
        lo.markovian = true;
        lo.assumelive = true;
        lo.init = m0;
        spn::SpnLpOptions up;
        up.markovian = true;
        up.init = m0;
        spn::SpnLpOptions op;
        op.markovian = false;
        op.assumelive = true;
        op.init = m0;
        const spn::SpnLpBounds bLo = spn::spn_lpbnd(sn, lo);
        const spn::SpnLpBounds bUp = spn::spn_lpbnd(sn, up);
        const spn::SpnLpBounds bOp = spn::spn_lpbnd(sn, op);
        double sumLo = 0.0, sumUp = 0.0, sumOpLo = 0.0, sumOpUp = 0.0;
        for (std::size_t e = 0; e < bLo.mode_tput_lo.size(); ++e) sumLo += bLo.mode_tput_lo[e];
        for (std::size_t e = 0; e < bUp.mode_tput_hi.size(); ++e) sumUp += bUp.mode_tput_hi[e];
        for (std::size_t e = 0; e < bOp.mode_tput_lo.size(); ++e) sumOpLo += bOp.mode_tput_lo[e];
        for (std::size_t e = 0; e < bOp.mode_tput_hi.size(); ++e) sumOpUp += bOp.mode_tput_hi[e];
        std::printf("%-5zu %9.4f %9.4f %9.4f   %9.3f %9.3f   %9.4f %9.4f\n", c + 1, sumLo,
                    pub[c][1], sumUp, pub[c][0], pub[c][2], sumOpLo, sumOpUp);
    }

    // The four method names are spnlp2.upper/lower and spnlp1.upper/lower.
    // counterparts, which drop the second-moment, covariance and Little's-law
    // families and so need only a mean firing time rather than an exponential
    // one. They are the only family SolverBA offers on a Petri net, and the only
    // one it withholds off a Petri net: every other family is parameterized by
    // demands and a population, which a marking is not.
    Net inh = spn_lpbounds_inhibiting_model(4.0);
    const Sn& sn_inh = inh.get_struct();
    std::string methods;
    const std::vector<std::string> valid = ba::list_valid_methods(sn_inh);
    for (std::size_t i = 0; i < valid.size(); ++i) methods += (i ? ", " : "") + valid[i];
    std::printf("\nmethods offered on this net: %s\n", methods.c_str());

    ba::BaOptions blo;
    blo.method = "spnlp2.lower";
    ba::BaOptions bup;
    bup.method = "spnlp2.upper";
    const mva::AvgResult<double> lo = ba::solver_ba_run_analyzer(sn_inh, blo);
    const mva::AvgResult<double> up = ba::solver_ba_run_analyzer(sn_inh, bup);
    const mva::AvgResult<double> ex =
        ctmc::solver_ctmc_run_analyzer(sn_inh, ctmc::CtmcOptions());

    std::printf("\nMean tokens per place, exact between the two sides:\n");
    std::printf("%-8s %10s %10s %10s\n", "place", "lower", "exact", "upper");
    for (std::size_t i = 0; i < sn_inh.nstations; ++i)
        std::printf("%-8s %10.5f %10.5f %10.5f\n", sn_inh.stations[i].name.c_str(), lo.QN(i, 0),
                    ex.QN(i, 0), up.QN(i, 0));

    // U = Q at a Place, which LINE models as an INF station -- the same
    // convention SolverCTMC and SolverNC report. The paper's place utilization
    // 1 - P(m = 0) is a different quantity and is not this column.
    std::printf("\nsolver.citations():\n");
    const std::vector<io::Citation> cits = io::line_citations_for_method("spnlp2.upper");
    for (std::size_t i = 0; i < cits.size(); ++i)
        std::printf("  [%s] %s\n      covers: %s\n", cits[i].key.c_str(), cits[i].ref.c_str(),
                    cits[i].covers.c_str());
}

LINE_EXAMPLE("basic/stochPetriNet", spn_basic_closed);
LINE_EXAMPLE("basic/stochPetriNet", spn_basic_open);
LINE_EXAMPLE("basic/stochPetriNet", spn_closed_fourplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_closed_twoplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_colored_gspn);
LINE_EXAMPLE("basic/stochPetriNet", spn_fluid_dae);
LINE_EXAMPLE("basic/stochPetriNet", spn_fourmodes);
LINE_EXAMPLE("basic/stochPetriNet", spn_inhibiting);
LINE_EXAMPLE("basic/stochPetriNet", spn_lpbounds);
LINE_EXAMPLE("basic/stochPetriNet", spn_open_sevenplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_pareto_service);
LINE_EXAMPLE("basic/stochPetriNet", spn_queueing_place);
LINE_EXAMPLE("basic/stochPetriNet", spn_twomodes);
LINE_EXAMPLE("basic/stochPetriNet", spn_productform_nc);
LINE_EXAMPLE("basic/stochPetriNet", test_spn_nrm_open);

}  // namespace examples
}  // namespace line
