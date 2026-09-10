/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/prioModel/`: the four class-priority models.
 *
 * Each one puts a priority discipline -- HOL, PSPRIO, GPSPRIO -- beside at least
 * one exact or simulated answer for the same network, so the priority
 * approximation can be read against something. JMT and LDES have no C++ port and
 * are refused by name; CTMC, MVA and SSA are the reference's own blocks and run.
 *
 * WHAT A PRIORITY IS HERE. The class priority is the `prio` argument of
 * `add_open_class` / `add_closed_class` -- a LARGER number is a HIGHER priority,
 * as in `ClosedClass(model, name, N, refnode, prio)` -- and the station decides
 * what to do with it through its scheduling strategy. A class carrying a
 * priority at a station whose discipline ignores one (FCFS, SIRO, PS) is not an
 * error and changes nothing there, which is exactly what `prio_hol_open` and
 * `prio_hol_closed` are built to show.
 */

#include <cstddef>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace line {
namespace examples {

namespace {

/** The reason every JMT block in this directory carries. */
const char* kNoJmt =
    "SolverJMT drives the Java Modelling Tools simulator, which has no C++ port; the native "
    "blocks above answer the same model, and reporting one of them under the simulator's name "
    "would attribute a number to an engine that never produced it";

/** The reason `prio_hol_closed` refuses its fourth block. */
const char* kNoLdes =
    "SolverLDES is the SSJ-based discrete-event engine of the JAR, driven over model.json by the "
    "MATLAB and Python clients; the C++ port carries no LDES engine and no client for it, and SSA "
    "above is a different simulator, not this one";

mva::AvgResult<double> mva_run(Net& m, double iter_tol = -1.0) {
    mva::MvaOptions opt;
    if (iter_tol > 0.0) opt.iter_tol = iter_tol;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

mva::AvgResult<double> ctmc_run(Net& m, double cutoff = -1.0) {
    ctmc::CtmcOptions opt;
    opt.cutoff = cutoff;
    return ctmc::solver_ctmc_run_analyzer_any(m.get_struct(), opt);
}

ssa::SsaSolution ssa_run(Net& m, std::size_t samples, unsigned long seed) {
    ssa::SsaOptions opt;
    opt.samples = samples;
    opt.seed = seed;
    return ssa::solver_ssa(m.get_struct(), opt);
}

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

/**
 * The closed six-node network: a slow Delay, four queues under FCFS, SIRO, PS
 * and HOL, and a fast Delay. Class2 carries priority 1, Class1 and Class3
 * priority 0, so only the HOL queue distinguishes them.
 */
Net prio_hol_closed_model() {
    Net m("MyNetwork");
    Delay slow(m, "SlowDelay");
    Queue fcfs(m, "FCFSQueue", SchedStrategy::FCFS);
    Queue siro(m, "SIROQueue", SchedStrategy::SIRO);
    Queue ps(m, "PSQueue", SchedStrategy::PS);
    Queue hol(m, "HOLQueue", SchedStrategy::HOL);
    Delay fast(m, "FastDelay");

    std::vector<std::size_t> cls;
    cls.push_back(m.add_closed_class("Class1", 18.0, slow, 0));
    cls.push_back(m.add_closed_class("Class2", 18.0, slow, 1));
    cls.push_back(m.add_closed_class("Class3", 18.0, slow, 0));

    const double slow_mean[3] = {10.0, 10.0, 10.0};
    const double fcfs_mean[3] = {0.3, 0.5, 0.6};
    const double siro_mean[3] = {1.1, 1.3, 1.5};
    const double ps_mean[3] = {1.0, 1.1, 1.9};
    const double hol_mean[3] = {2.5, 1.9, 4.3};
    const double fast_mean[3] = {1.0, 1.0, 1.0};
    for (std::size_t r = 0; r < 3; ++r) {
        slow.set_service(cls[r], D::exp_mean(slow_mean[r]));
        fcfs.set_service(cls[r], D::exp_mean(fcfs_mean[r]));
        siro.set_service(cls[r], D::exp_mean(siro_mean[r]));
        ps.set_service(cls[r], D::exp_mean(ps_mean[r]));
        hol.set_service(cls[r], D::exp_mean(hol_mean[r]));
        fast.set_service(cls[r], D::exp_mean(fast_mean[r]));
    }

    Routing P;
    for (std::size_t r = 0; r < 3; ++r) {
        P.set(cls[r], cls[r], slow, fcfs, 1.0);
        P.set(cls[r], cls[r], fcfs, siro, 0.25);
        P.set(cls[r], cls[r], fcfs, ps, 0.25);
        P.set(cls[r], cls[r], fcfs, hol, 0.25);
        P.set(cls[r], cls[r], fcfs, fast, 0.25);
        P.set(cls[r], cls[r], siro, fcfs, 1.0);
        P.set(cls[r], cls[r], ps, fcfs, 1.0);
        P.set(cls[r], cls[r], hol, fcfs, 1.0);
        P.set(cls[r], cls[r], fast, slow, 1.0);
    }
    m.link(P);
    return m;
}

/** The open counterpart: a Source, the same four disciplines, and a Sink. */
Net prio_hol_open_model() {
    Net m("MyNetwork");
    Source source(m, "Source");
    Queue web(m, "WebServer", SchedStrategy::FCFS);
    Queue st1(m, "Storage1", SchedStrategy::SIRO);
    Queue st2(m, "Storage2", SchedStrategy::PS);
    Queue st3(m, "Storage3", SchedStrategy::HOL);
    Sink sink(m, "Sink");

    std::vector<std::size_t> cls;
    cls.push_back(m.add_open_class("Class1", 0));
    cls.push_back(m.add_open_class("Class2", 1));
    cls.push_back(m.add_open_class("Class3", 0));

    const double web_mean[3] = {0.3, 0.5, 0.6};
    const double st1_mean[3] = {1.1, 1.3, 1.5};
    const double st2_mean[3] = {2.0, 2.1, 1.9};
    const double st3_mean[3] = {2.5, 1.9, 4.3};
    for (std::size_t r = 0; r < 3; ++r) {
        source.set_arrival(cls[r], D::exp_mean(10.0));
        web.set_service(cls[r], D::exp_mean(web_mean[r]));
        st1.set_service(cls[r], D::exp_mean(st1_mean[r]));
        st2.set_service(cls[r], D::exp_mean(st2_mean[r]));
        st3.set_service(cls[r], D::exp_mean(st3_mean[r]));
    }

    Routing P;
    for (std::size_t r = 0; r < 3; ++r) {
        P.set(cls[r], cls[r], source, web, 1.0);
        P.set(cls[r], cls[r], web, st1, 0.25);
        P.set(cls[r], cls[r], web, st2, 0.25);
        P.set(cls[r], cls[r], web, st3, 0.25);
        P.set(cls[r], cls[r], web, sink, 0.25);
        P.set(cls[r], cls[r], st1, web, 1.0);
        P.set(cls[r], cls[r], st2, web, 1.0);
        P.set(cls[r], cls[r], st3, web, 1.0);
    }
    m.link(P);
    return m;
}

/**
 * GPSPRIO with a TIE: Class2 and Class3 both carry priority 1, so the queue
 * resolves between them by their GPS weights rather than by their priorities.
 */
Net prio_identical_model() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue queue(m, "Queue1", SchedStrategy::GPSPRIO);

    std::vector<std::size_t> cls;
    cls.push_back(m.add_closed_class("Class1", 6.0, delay, 0));
    cls.push_back(m.add_closed_class("Class2", 4.0, delay, 1));
    cls.push_back(m.add_closed_class("Class3", 4.0, delay, 1));
    cls.push_back(m.add_closed_class("Class4", 1.0, delay, 2));

    delay.set_service(cls[0], Erlang(3.0, 2));
    delay.set_service(cls[1], Exp(1.0));
    delay.set_service(cls[2], Exp(1.0));
    delay.set_service(cls[3], Exp(2.0));

    const double rate[4] = {30.0, 2.0, 12.0, 1.0};
    const double weight[4] = {12.0, 3.0, 5.0, 1.0};
    for (std::size_t r = 0; r < 4; ++r) {
        queue.set_service(cls[r], Exp(rate[r]));
        queue.set_sched_param(cls[r], weight[r]);
    }

    Routing P;
    for (std::size_t r = 0; r < 4; ++r) cyclic(P, cls[r], {delay, queue});
    m.link(P);
    return m;
}

/** PSPRIO: processor sharing inside each priority level, strict between them. */
Net prio_psprio_model() {
    Net m("MyNetwork");
    Delay slow(m, "SlowDelay");
    Queue queue(m, "PSPRIOQueue", SchedStrategy::PSPRIO);

    ClosedClass c1(m, "Class1", 2.0, slow, 0);
    ClosedClass c2(m, "Class2", 2.0, slow, 1);

    slow.set_service(c1, Erlang(3.0, 2));
    slow.set_service(c2, HyperExp(0.5, 3.0, 10.0));
    queue.set_service(c1, HyperExp(0.1, 1.0, 10.0));
    queue.set_service(c2, Exp(1.0));

    Routing P;
    cyclic(P, c1, {slow, queue});
    cyclic(P, c2, {slow, queue});
    m.link(P);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Examples
// ---------------------------------------------------------------------------

void prio_hol_closed() {
    Net m = prio_hol_closed_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
    section("SSA");
    print_avg_sim(m.get_struct(), ssa_run(m, 10000, 23000));
    // TODO(cpp): print(LDES(model, seed=23000, samples=10000).avg_table())
    na("LDES", kNoLdes);
}

/**
 * `cutoff=1` on the CTMC block is the reference's own truncation: the chain of
 * an open network is infinite, so its answer is the answer of the truncated
 * chain and is reported as such rather than as the model's.
 */
void prio_hol_open() {
    Net m = prio_hol_open_model();
    section("CTMC");
    print_avg(m.get_struct(), ctmc_run(m, 1.0));
    // The MVA-specific iter_tol (1e-6), not the generic 1e-4: the looser one
    // stopped the AMVA fixed point 2.9e-5 short on Storage3 and split this row
    // from the engines that solve at the solver's own tolerance.
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
    section("SSA");
    print_avg_sim(m.get_struct(), ssa_run(m, 10000, 23000));
}

void prio_identical() {
    Net m = prio_identical_model();
    section("CTMC");
    print_avg(m.get_struct(), ctmc_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
}

void prio_psprio() {
    Net m = prio_psprio_model();
    section("CTMC");
    print_avg(m.get_struct(), ctmc_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 5000)));
    section("SSA");
    print_avg_sim(m.get_struct(), ssa_run(m, 5000, 23000));
}

LINE_EXAMPLE("basic/prioModel", prio_hol_closed);
LINE_EXAMPLE("basic/prioModel", prio_hol_open);
LINE_EXAMPLE("basic/prioModel", prio_identical);
LINE_EXAMPLE("basic/prioModel", prio_psprio);

}  // namespace examples
}  // namespace line
