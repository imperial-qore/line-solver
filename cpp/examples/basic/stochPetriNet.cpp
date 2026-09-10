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
 * ONE OF THEM CANNOT BE BUILT IN FULL, and is refused rather than approximated:
 *
 *  - `spn_queueing_place` needs QUEUEING places, whose embedded queue has its
 *    own scheduling discipline. `add_place(name)` builds the INF pass-through
 *    place only, so the QPN half is refused; its exact cross-check, the
 *    equivalent Delay + M/M/1 finite-population network under SolverMVA, IS in
 *    the reference's `__main__` and does run here.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "line/api/spn/spn_metrics.h"
#include "line/api/spn/spn_pf.h"

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/distribution.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

typedef qn::TransitionParam<double> TransParam;
using lang::TimingStrategy;

/** The reason every JMT block in this directory carries. */
const char* kNoJmt =
    "SolverJMT drives the Java Modelling Tools simulator, which has no C++ port; the reference "
    "runs no second solver on this net, so the model is built and printed rather than answered "
    "by an engine it never named";

/**
 * One arc of a mode: the node it touches, its multiplicity, and the CLASS whose
 * tokens it moves.
 *
 * `cls` is 1-based and defaults to the first class, which is what a single-class
 * net means. It is not decoration: a mode requiring two Class1 tokens at a place
 * must not be enabled by Class2 tokens sitting there.
 */
struct Arc {
    std::size_t node;
    double count;
    std::size_t cls = 1;
};

/**
 * One mode of a Transition: its firing law and its pre, post and inhibitor arcs.
 *
 * The defaults are `Transition.add_mode`'s own -- one server, TIMED, priority
 * ONE and weight one. The priority is 1.0 and not 0.0: `nodes.py:2660` appends
 * `1.0`, and MATLAB's writer emits `"firingPriority": 1` for a transition that
 * never had one set, so a mode built here without an explicit priority must
 * carry the same number a mode built there does.
 */
struct Mode {
    std::string name;
    D proc = Disabled();
    TimingStrategy timing = TimingStrategy::TIMED;
    double servers = 1.0;
    double prio = 1.0;
    double weight = 1.0;
    std::vector<Arc> enab, inhib, fire;
};

/**
 * `Transition.addMode` and its setters, collected into the struct's own shape.
 *
 * The arc matrices are indexed by (1-based NODE, 1-based CLASS) and sized
 * against the FULL node count of the finished model, which is what
 * `network_reader.h` does: a transition is itself a node, so the count cannot be
 * read off the model while the transitions are still being added. An inhibitor
 * defaults to infinity, the encoding for "this place never blocks the mode".
 */
TransParam trans(std::size_t nnodes, const std::vector<Mode>& modes, std::size_t nclasses = 1) {
    const double inf = std::numeric_limits<double>::infinity();
    TransParam tp;
    tp.nmodes = modes.size();
    for (std::size_t m = 0; m < modes.size(); ++m) {
        const Mode& md = modes[m];
        const bool immediate = md.timing == TimingStrategy::IMMEDIATE;
        tp.modenames.push_back(md.name);
        tp.timing.push_back(md.timing);
        tp.firingproc.push_back(md.proc);
        tp.firingphases.push_back(immediate || md.proc.disabled ? 0
                                                                : lang::dist_to_map(md.proc).order());
        tp.nmodeservers.push_back(md.servers);
        tp.firingprio.push_back(md.prio);
        tp.fireweight.push_back(md.weight);
        Matrix<double> enab(nnodes, nclasses, 0.0), inhib(nnodes, nclasses, inf),
            fire(nnodes, nclasses, 0.0);
        for (std::size_t a = 0; a < md.enab.size(); ++a)
            enab(md.enab[a].node - 1, md.enab[a].cls - 1) += md.enab[a].count;
        for (std::size_t a = 0; a < md.inhib.size(); ++a)
            inhib(md.inhib[a].node - 1, md.inhib[a].cls - 1) = md.inhib[a].count;
        for (std::size_t a = 0; a < md.fire.size(); ++a)
            fire(md.fire[a].node - 1, md.fire[a].cls - 1) += md.fire[a].count;
        tp.enabling.push_back(enab);
        tp.inhibiting.push_back(inhib);
        tp.firing.push_back(fire);
    }
    return tp;
}

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

    std::vector<Mode> modes;
    const char* names[3] = {"Mode1", "Mode2", "Mode3"};
    const D procs[3] = {D::exp_mean(1.0), lang::erlang_fit_mean_order<double>(1.0, 2),
                        lang::hyperexp_fit_mean_scv<double>(1.0, 4.0)};
    for (int k = 0; k < 3; ++k) {
        Mode md;
        md.name = names[k];
        md.proc = procs[k];
        md.enab.push_back(Arc{p1, 1.0});
        md.fire.push_back(Arc{p1, 1.0});
        modes.push_back(md);
    }
    const std::size_t t1 = m.add_transition("T1", trans(2, modes));

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

    Mode md;
    md.name = "Mode1";
    md.servers = kMaxInt;
    md.proc = Exp(4.0);
    md.enab.push_back(Arc{p1, 1.0});
    md.fire.push_back(Arc{sink, 1.0});
    const std::size_t t1 = m.add_transition("T1", trans(4, std::vector<Mode>(1, md)));

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
    std::vector<std::size_t> p;
    for (int i = 0; i < 4; ++i) p.push_back(m.add_place("P" + std::to_string(i + 1)));
    ClosedClass c1(m, "Class1", 2.0, p[0], 0);

    std::vector<double> mu0, phi0;
    mu0.push_back(1.0);
    mu0.push_back(2.0);
    phi0.push_back(0.6);
    phi0.push_back(1.0);
    const D procs[4] = {Exp(2.0), Erlang(3.0, 4), HyperExp(0.7, 3.0, 1.5),
                        Coxian(mu0, phi0)};

    std::vector<std::size_t> t;
    for (int i = 0; i < 4; ++i) {
        Mode md;
        md.name = "Mode" + std::to_string(i + 1);
        md.proc = procs[i];
        md.enab.push_back(Arc{p[i], 2.0});
        md.fire.push_back(Arc{p[(i + 1) % 4], 2.0});
        t.push_back(m.add_transition("T" + std::to_string(i + 1), trans(8, std::vector<Mode>(1, md))));
    }

    Routing R;
    for (int i = 0; i < 4; ++i) {
        R.set(c1, c1, p[i], t[i], 1.0);
        R.set(c1, c1, t[i], p[(i + 1) % 4], 1.0);
    }
    m.link(R);
    for (int i = 0; i < 4; ++i)
        m.set_initial_marking(p[i], std::vector<double>(1, i == 0 ? 2.0 : 0.0));
    return m;
}

/** P1 branches to P2 (2 tokens) or P3 (3 tokens); both return to P1. */
Net spn_fourmodes_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    Place p3(m, "P3");
    ClosedClass c1(m, "Class1", 8.0, p1, 0);

    const std::size_t from[4] = {p1, p1, p2, p3};
    const std::size_t to[4] = {p2, p3, p1, p1};
    const double count[4] = {2.0, 3.0, 1.0, 2.0};
    const double rate[4] = {2.0, 1.0, 4.0, 2.0};
    std::vector<std::size_t> t;
    for (int i = 0; i < 4; ++i) {
        Mode md;
        md.name = "Mode" + std::to_string(i + 1);
        md.proc = Exp(rate[i]);
        md.enab.push_back(Arc{from[i], count[i]});
        md.fire.push_back(Arc{to[i], count[i]});
        t.push_back(m.add_transition("T" + std::to_string(i + 1), trans(7, std::vector<Mode>(1, md))));
    }

    Routing R;
    for (int i = 0; i < 4; ++i) {
        R.set(c1, c1, from[i], t[i], 1.0);
        R.set(c1, c1, t[i], to[i], 1.0);
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

    std::vector<Mode> t1modes;
    Mode m1;
    m1.name = "Mode1";
    m1.proc = Exp(2.0);
    m1.enab.push_back(Arc{p1, 2.0});
    m1.fire.push_back(Arc{p2, 2.0});
    t1modes.push_back(m1);
    Mode m2;
    m2.name = "Mode2";
    m2.proc = Exp(1.0);
    m2.enab.push_back(Arc{p1, 1.0});
    m2.fire.push_back(Arc{p3, 1.0});
    t1modes.push_back(m2);
    const std::size_t t1 = m.add_transition("T1", trans(6, t1modes));

    Mode m3;
    m3.name = "Mode3";
    m3.proc = Exp(4.0);
    m3.enab.push_back(Arc{p2, 1.0});
    m3.fire.push_back(Arc{p1, 1.0});
    const std::size_t t2 = m.add_transition("T2", trans(6, std::vector<Mode>(1, m3)));

    Mode m4;
    m4.name = "Mode4";
    m4.proc = Exp(1.0);
    m4.enab.push_back(Arc{p3, 3.0});
    m4.inhib.push_back(Arc{p2, 1.0});
    m4.fire.push_back(Arc{p1, 3.0});
    const std::size_t t3 = m.add_transition("T3", trans(6, std::vector<Mode>(1, m4)));

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
    std::vector<std::size_t> p;
    for (int i = 0; i < 7; ++i) p.push_back(m.add_place("P" + std::to_string(i + 1)));
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, D::exp_mean(1.0));

    // Source, Sink, seven places and eight transitions: the arc vectors are
    // sized against the finished model, which the transitions are part of.
    const std::size_t nn = 17;
    std::vector<std::size_t> t;

    Mode m1;
    m1.name = "Mode1";
    m1.servers = kMaxInt;
    m1.proc = Exp(4.0);
    m1.enab.push_back(Arc{p[0], 1.0});
    m1.fire.push_back(Arc{p[1], 1.0});
    t.push_back(m.add_transition("T1", trans(nn, std::vector<Mode>(1, m1))));

    Mode m2;
    m2.name = "Mode1";
    m2.servers = kMaxInt;
    m2.timing = TimingStrategy::IMMEDIATE;
    m2.prio = 1.0;
    m2.weight = 1.0;
    m2.enab.push_back(Arc{p[1], 1.0});
    m2.fire.push_back(Arc{p[2], 1.0});
    t.push_back(m.add_transition("T2", trans(nn, std::vector<Mode>(1, m2))));

    Mode m3;
    m3.name = "Mode1";
    m3.servers = kMaxInt;
    m3.timing = TimingStrategy::IMMEDIATE;
    m3.prio = 1.0;
    m3.enab.push_back(Arc{p[1], 1.0});
    m3.fire.push_back(Arc{p[3], 1.0});
    t.push_back(m.add_transition("T3", trans(nn, std::vector<Mode>(1, m3))));

    Mode m4;
    m4.name = "Mode1";
    m4.servers = kMaxInt;
    m4.timing = TimingStrategy::IMMEDIATE;
    m4.prio = 1.0;
    m4.enab.push_back(Arc{p[2], 1.0});
    m4.enab.push_back(Arc{p[4], 1.0});
    m4.fire.push_back(Arc{p[4], 1.0});
    m4.fire.push_back(Arc{p[5], 1.0});
    t.push_back(m.add_transition("T4", trans(nn, std::vector<Mode>(1, m4))));

    Mode m5;
    m5.name = "Mode1";
    m5.servers = kMaxInt;
    m5.timing = TimingStrategy::IMMEDIATE;
    m5.prio = 1.0;
    m5.enab.push_back(Arc{p[3], 1.0});
    m5.enab.push_back(Arc{p[4], 1.0});
    m5.inhib.push_back(Arc{p[5], 1.0});
    m5.fire.push_back(Arc{p[6], 1.0});
    t.push_back(m.add_transition("T5", trans(nn, std::vector<Mode>(1, m5))));

    Mode m6;
    m6.name = "Mode1";
    m6.servers = kMaxInt;
    m6.proc = Erlang(2.0, 2);
    m6.enab.push_back(Arc{p[5], 1.0});
    m6.fire.push_back(Arc{p[0], 1.0});
    t.push_back(m.add_transition("T6", trans(nn, std::vector<Mode>(1, m6))));

    Mode m7;
    m7.name = "Mode1";
    m7.servers = kMaxInt;
    m7.proc = Exp(2.0);
    m7.enab.push_back(Arc{p[6], 1.0});
    m7.fire.push_back(Arc{p[0], 1.0});
    m7.fire.push_back(Arc{p[4], 1.0});
    t.push_back(m.add_transition("T7", trans(nn, std::vector<Mode>(1, m7))));

    Mode m8;
    m8.name = "Mode1";
    m8.servers = kMaxInt;
    m8.proc = Exp(2.0);
    m8.enab.push_back(Arc{p[3], 1.0});
    m8.fire.push_back(Arc{sink, 1.0});
    t.push_back(m.add_transition("T8", trans(nn, std::vector<Mode>(1, m8))));

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
    for (int i = 0; i < 7; ++i) m.set_initial_marking(p[i], std::vector<double>(1, init[i]));
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

    Mode md;
    md.name = "Mode1";
    md.servers = 1.0;
    md.proc = Pareto(3.0, 1.0);
    md.enab.push_back(Arc{p1, 1.0});
    md.fire.push_back(Arc{sink, 1.0});
    const std::size_t t1 = m.add_transition("T1", trans(4, std::vector<Mode>(1, md)));

    Routing R;
    R.set(c1, c1, source, p1, 1.0);
    R.set(c1, c1, p1, t1, 1.0);
    R.set(c1, c1, t1, sink, 1.0);
    m.link(R);
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
    const std::size_t nn = 5;  // P1, P2, T1, T2, T3

    Mode m1;
    m1.name = "Mode1";
    m1.proc = Exp(2.0);
    m1.enab.push_back(Arc{p1, 2.0, c1});
    m1.fire.push_back(Arc{p2, 2.0, c1});
    Mode m2;
    m2.name = "Mode2";
    m2.proc = Exp(3.0);
    m2.enab.push_back(Arc{p1, 1.0, c2});
    m2.fire.push_back(Arc{p2, 1.0, c2});
    std::vector<Mode> t1modes;
    t1modes.push_back(m1);
    t1modes.push_back(m2);
    const std::size_t t1 = m.add_transition("T1", trans(nn, t1modes, 2));

    Mode m3;
    m3.name = "Mode3";
    m3.proc = Erlang(1.5, 2);  // Erlang(1.5, 2): rate 1.5 per phase, mean 2/1.5
    m3.enab.push_back(Arc{p2, 1.0, c1});
    m3.fire.push_back(Arc{p1, 1.0, c1});
    const std::size_t t2 = m.add_transition("T2", trans(nn, std::vector<Mode>(1, m3), 2));

    Mode m4;
    m4.name = "Mode4";
    m4.proc = Exp(0.5);
    m4.enab.push_back(Arc{p2, 4.0, c2});
    m4.fire.push_back(Arc{p1, 4.0, c2});
    const std::size_t t3 = m.add_transition("T3", trans(nn, std::vector<Mode>(1, m4), 2));

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

/** P1 -> P2 in batches of four, P2 -> P1 in batches of two. */
Net spn_twomodes_model() {
    Net m("model");
    Place p1(m, "P1");
    Place p2(m, "P2");
    ClosedClass c1(m, "Class1", 10.0, p1, 0);

    Mode m1;
    m1.name = "Mode1";
    m1.proc = Exp(2.0);
    m1.enab.push_back(Arc{p1, 4.0});
    m1.fire.push_back(Arc{p2, 4.0});
    const std::size_t t1 = m.add_transition("T1", trans(4, std::vector<Mode>(1, m1)));

    Mode m2;
    m2.name = "Mode2";
    m2.proc = Exp(3.0);
    m2.enab.push_back(Arc{p2, 2.0});
    m2.fire.push_back(Arc{p1, 2.0});
    const std::size_t t2 = m.add_transition("T2", trans(4, std::vector<Mode>(1, m2)));

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

    Mode md;
    md.name = "Mode1";
    md.proc = Exp(mu);
    md.enab.push_back(Arc{p1, 1.0});
    md.fire.push_back(Arc{sink, 1.0});
    const std::size_t t1 = m.add_transition("T1", trans(4, std::vector<Mode>(1, md)));

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

    Mode m1;
    m1.name = "Mode1";
    m1.proc = Exp(mu1);
    m1.enab.push_back(Arc{p1, 1.0});
    m1.fire.push_back(Arc{p2, 1.0});
    const std::size_t t1 = m.add_transition("T1", trans(6, std::vector<Mode>(1, m1)));

    Mode m2;
    m2.name = "Mode1";
    m2.proc = Exp(mu2);
    m2.enab.push_back(Arc{p2, 1.0});
    m2.fire.push_back(Arc{sink, 1.0});
    const std::size_t t2 = m.add_transition("T2", trans(6, std::vector<Mode>(1, m2)));

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

        // TODO(cpp): jmt = SolverJMT(mm1spn(lambda, mu), samples=SAMPLES, seed=SEED)
        // TODO(cpp): assert abs(Qn(2) - Qj(2)) / Qj(2) < RTOL
        na("JMT", kNoJmt);
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
 * The QPN itself is refused: its CPU and Think places are QUEUEING places, whose
 * embedded queue is served under a scheduling discipline of its own, and
 * `add_place(name)` builds the INF pass-through place only. The reference's
 * second block, the equivalent Delay + M/M/1 finite-population network under
 * SolverMVA, is a plain queueing network and runs.
 */
void spn_queueing_place() {
    const double N = 4.0;
    // TODO(cpp): print(SolverLDES(model, seed=23000, samples=200000).avgTable())
    na("LDES",
       "SolverLDES has no C++ port, and the model it is given here cannot be built either: its "
       "CPU and Think places are QUEUEING places with an embedded FCFS and INF queue, which the "
       "port's add_place(name) -- an INF pass-through place with no service process -- cannot "
       "express. Building it as an ordinary place would simulate a different net");
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
    std::vector<std::size_t> p;
    for (int i = 0; i < 3; ++i) p.push_back(m.add_place("P" + std::to_string(i)));
    ClosedClass c1(m, "Class1", ntokens, p[0], 0);
    std::vector<std::size_t> t;
    for (int i = 0; i < 3; ++i) {
        Mode md;
        md.name = "fire";
        md.proc = Exp(rates[i]);
        md.enab.push_back(Arc{p[i], 1.0});
        md.fire.push_back(Arc{p[(i + 1) % 3], 1.0});
        t.push_back(m.add_transition("T" + std::to_string(i), trans(6, std::vector<Mode>(1, md))));
    }
    Routing R;
    for (int i = 0; i < 3; ++i) {
        R.set(c1, c1, p[i], t[i], 1.0);
        R.set(c1, c1, t[i], p[(i + 1) % 3], 1.0);
    }
    m.link(R);
    for (int i = 0; i < 3; ++i)
        m.set_initial_marking(p[i], std::vector<double>(1, i == 0 ? ntokens : 0.0));
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
    std::vector<std::size_t> p;
    for (int i = 0; i < 4; ++i) p.push_back(m.add_place("P" + std::to_string(i)));
    ClosedClass c1(m, "C", ntokens, p[0], 0);

    Mode mf;
    mf.name = "f";
    mf.proc = Exp(1.3);
    mf.enab.push_back(Arc{p[0], 1.0});
    mf.fire.push_back(Arc{p[1], 1.0});
    mf.fire.push_back(Arc{p[2], 1.0});
    const std::size_t tf = m.add_transition("Tf", trans(7, std::vector<Mode>(1, mf)));

    Mode mj;
    mj.name = "j";
    mj.proc = Exp(0.7);
    mj.enab.push_back(Arc{p[1], 1.0});
    mj.enab.push_back(Arc{p[2], 1.0});
    mj.fire.push_back(Arc{p[3], 1.0});
    const std::size_t tj = m.add_transition("Tj", trans(7, std::vector<Mode>(1, mj)));

    Mode mb;
    mb.name = "b";
    mb.proc = Exp(1.9);
    mb.enab.push_back(Arc{p[3], 1.0});
    mb.fire.push_back(Arc{p[0], 1.0});
    const std::size_t tb = m.add_transition("Tb", trans(7, std::vector<Mode>(1, mb)));

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
        m.set_initial_marking(p[i], std::vector<double>(1, i == 0 ? ntokens : 0.0));
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

LINE_EXAMPLE("basic/stochPetriNet", spn_basic_closed);
LINE_EXAMPLE("basic/stochPetriNet", spn_basic_open);
LINE_EXAMPLE("basic/stochPetriNet", spn_closed_fourplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_closed_twoplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_fourmodes);
LINE_EXAMPLE("basic/stochPetriNet", spn_inhibiting);
LINE_EXAMPLE("basic/stochPetriNet", spn_open_sevenplaces);
LINE_EXAMPLE("basic/stochPetriNet", spn_pareto_service);
LINE_EXAMPLE("basic/stochPetriNet", spn_queueing_place);
LINE_EXAMPLE("basic/stochPetriNet", spn_twomodes);
LINE_EXAMPLE("basic/stochPetriNet", spn_productform_nc);
LINE_EXAMPLE("basic/stochPetriNet", test_spn_nrm_open);

}  // namespace examples
}  // namespace line
