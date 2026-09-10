/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/mixedQN/`: mixed networks, one closed class and one
 * open class sharing the same stations.
 *
 * The six entries differ in the discipline (FCFS or PS), the server
 * multiplicity and the closed population, and `mqn_multichain_cs` adds class
 * switching inside both an open and a closed chain. Every one of them runs a
 * `SolverJMT` and a `SolverLDES` block the reference prints; neither engine has
 * a C++ port, so those blocks are refused BY NAME and no native solver is put
 * in their place.
 *
 * `mqn_multichain_cs` prints `getAvgChainTable`, which this port has no facade
 * for; the aggregation is `matlab/src/solvers/@NetworkSolver/getAvg*Chain.m`
 * transcribed below -- a sum over `sn.inchain` for every column but RespT,
 * which is the visit-share weighted sum with the `alpha` of
 * `sn_get_demands_chain`.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
#include "line/lang/dist_fitters.h"
#include "line/solvers/mva/sn_chain.h"

namespace line {
namespace examples {

namespace {

/** One cell of a solved table; a (station, class) pair with no row is a zero. */
double cell(const AvgTable& t, const std::string& col, const std::string& station,
            const std::string& jobclass) {
    const double v = t.get(col, station, jobclass);
    return std::isnan(v) ? 0.0 : v;
}

/**
 * `getAvgChainTable`: the class rows of `t` aggregated over `sn.inchain`.
 *
 * QLen, Util, ResidT, ArvR and Tput are summed over the classes of the chain;
 * RespT is weighted by `alpha`, the class share of the chain's visits at the
 * station, because response times do not add across classes that a job moves
 * between while circulating.
 */
void print_avg_chain(const Sn& sn, const AvgTable& t) {
    const Matrix<double> alpha = mva::sn_get_demands_chain(sn).alpha;
    std::printf("Solver%s method=%s\n", t.solver.c_str(), t.method.c_str());
    std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Station", "Chain", "QLen", "Util",
                "RespT", "ResidT", "ArvR", "Tput");
    // A CHAIN TABLE DOES NOT GO THROUGH `avg_rows`, so nothing records it unless
    // this loop does -- and it is the only table these examples print, so the
    // whole row went unmeasured. The aggregation above is the same one
    // `solvers::solver_get_avg_chain` performs, and the chain names are the same
    // `Chain1..ChainC` the goldens key on.
    namespace parity = line::examples::parity;
    if (parity::enabled()) {
        parity::set_solver_from_table(t.solver);
        parity::begin_table("chain", "Station", "JobClass");
    }
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nchains; ++c) {
            double q = 0.0, u = 0.0, r = 0.0, w = 0.0, a = 0.0, x = 0.0;
            for (std::size_t k : sn.inchain[c]) {
                const std::string& st = sn.stations[i].name;
                const std::string& cl = sn.classes[k - 1].name;
                q += cell(t, "QLen", st, cl);
                u += cell(t, "Util", st, cl);
                r += cell(t, "RespT", st, cl) * alpha(i, k - 1);
                w += cell(t, "ResidT", st, cl);
                a += cell(t, "ArvR", st, cl);
                x += cell(t, "Tput", st, cl);
            }
            std::printf("%-16s Chain%-9zu %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                        sn.stations[i].name.c_str(), c + 1, q, u, r, w, a, x);
            if (!parity::enabled()) continue;
            std::vector<parity::Cell> cells;
            cells.push_back(parity::Cell{"QLen", q});
            cells.push_back(parity::Cell{"Util", u});
            cells.push_back(parity::Cell{"RespT", r});
            cells.push_back(parity::Cell{"ResidT", w});
            cells.push_back(parity::Cell{"ArvR", a});
            cells.push_back(parity::Cell{"Tput", x});
            parity::add_row(sn.stations[i].name, "Chain" + std::to_string(c + 1), cells);
        }
}

/** The five-queue mixed model of the two multiserver entries. */
Net multiserver_model(SchedStrategy sched) {
    Net model("model");
    std::vector<std::size_t> node;
    for (std::size_t i = 0; i < 5; ++i) {
        node.push_back(model.add_queue("Queue" + std::to_string(i + 1), sched));
        model.set_number_of_servers(node[i], static_cast<double>(i + 1));
    }
    Source source(model, "Source");
    Sink sink(model, "Sink");

    ClosedClass closed(model, "ClosedClass", 3, node[0], 0);
    OpenClass open(model, "OpenClass", 0);
    for (std::size_t i = 0; i < 5; ++i) {
        model.set_service(node[i], closed, Exp(static_cast<double>(i + 1)));
        model.set_service(node[i], open, Exp(std::sqrt(static_cast<double>(i + 1))));
    }
    source.set_arrival(open, Exp(0.3));

    Routing P;
    cyclic(P, closed, {node[0], node[1], node[2], node[3]});
    serial(P, open, {source, node[0], node[1], node[2], node[4], sink});
    model.link(P);
    return model;
}

/** The four-queue single-server mixed model, with the arrival law as a parameter. */
Net singleserver_model(SchedStrategy sched, const D& arrival) {
    Net model("model");
    std::vector<std::size_t> node;
    for (std::size_t i = 0; i < 4; ++i)
        node.push_back(model.add_queue("Queue" + std::to_string(i + 1), sched));
    Source source(model, "Source");
    Sink sink(model, "Sink");

    ClosedClass closed(model, "ClosedClass", 100, node[0], 0);
    OpenClass open(model, "OpenClass", 0);
    for (std::size_t i = 0; i < 4; ++i) {
        model.set_service(node[i], closed, Exp(static_cast<double>(i + 1)));
        model.set_service(node[i], open, Exp(std::sqrt(static_cast<double>(i + 1))));
    }
    source.set_arrival(open, arrival);

    Routing P;
    cyclic(P, closed, {node[0], node[1], node[2], node[3]});
    serial(P, open, {source, node[0], node[1], node[2], sink});
    model.link(P);
    return model;
}

}  // namespace

/**
 * A Delay and a PS Queue shared by two closed jobs and a Poisson open stream,
 * with Erlang and hyperexponential service.
 */
void mqn_basic() {
    Net model("model");
    Delay delay(model, "Delay");
    Queue queue1(model, "Queue1", SchedStrategy::PS);
    Source source(model, "Source");
    Sink sink(model, "Sink");

    ClosedClass closed(model, "ClosedClass", 2, delay, 0);
    OpenClass open(model, "OpenClass", 0);

    delay.set_service(closed, Erlang(3.0, 2));
    delay.set_service(open, HyperExp(0.5, 3.0, 10.0));
    queue1.set_service(closed, HyperExp(0.1, 1.0, 10.0));
    queue1.set_service(open, Exp(1.0));
    source.set_arrival(open, Exp(0.1));

    Routing P;
    P.set(closed, closed, delay, queue1, 1.0);
    P.set(closed, closed, queue1, delay, 1.0);
    P.set(open, open, source, delay, 1.0);
    P.set(open, open, delay, queue1, 1.0);
    P.set(open, open, queue1, sink, 1.0);
    model.link(P);

    SolverOpts ctmc;
    ctmc.cutoff = 3.0;
    ctmc.seed = 23000;
    section("CTMC");
    print_avg(solve_avg("CTMC", model, ctmc));

    section("JMT");
    print_avg(model.get_struct(), jmt_avg(model, sim_opts(23000)));

    SolverOpts ssa;
    ssa.seed = 23000;
    section("SSA");
    print_avg(solve_avg("SSA", model, ssa));

    section("MVA");
    print_avg(solve_avg("MVA", model));

    // TODO(cpp): LDES(model, keep=True, verbose=True, seed=23000).avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

/**
 * Three chains over four stations: two closed classes that switch into each
 * other at the WebServer, two open classes that do the same, and a third open
 * class that does not switch at all.
 */
void mqn_multichain_cs() {
    Net model("mqn_multichain");
    Delay think(model, "ThinkingTime");
    Queue web(model, "WebServer", SchedStrategy::FCFS);
    Queue app(model, "AppServer", SchedStrategy::PS);
    Queue data(model, "DataServer", SchedStrategy::FCFS);
    Source source(model, "Source");
    Sink sink(model, "Sink");

    ClosedClass ia(model, "InteractiveA", 3, think);
    ClosedClass ib(model, "InteractiveB", 2, think);
    OpenClass ba(model, "BatchA", 0);
    OpenClass bb(model, "BatchB", 0);
    OpenClass ext(model, "ExternalLoad", 0);

    think.set_service(ia, D::exp_mean(1.5));
    think.set_service(ib, D::exp_mean(2.0));

    web.set_service(ia, D::exp_mean(0.4));
    web.set_service(ib, D::exp_mean(0.5));
    web.set_service(ba, D::exp_mean(0.3));
    web.set_service(bb, D::exp_mean(0.35));
    web.set_service(ext, D::exp_mean(0.2));

    app.set_service(ia, D::exp_mean(0.8));
    app.set_service(ib, D::exp_mean(1.0));
    app.set_service(ba, D::exp_mean(0.6));
    app.set_service(bb, D::exp_mean(0.7));
    app.set_service(ext, D::exp_mean(0.5));

    data.set_service(ia, D::exp_mean(0.5));
    data.set_service(ib, D::exp_mean(0.6));
    data.set_service(ba, D::exp_mean(1.2));
    data.set_service(bb, D::exp_mean(1.5));
    data.set_service(ext, D::exp_mean(0.8));

    source.set_arrival(ba, Exp(0.1));
    source.set_arrival(bb, Exp(0.05));
    source.set_arrival(ext, Exp(0.05));

    Routing P;
    P.set(ia, ia, think, web, 1.0);
    P.set(ia, ia, web, app, 0.7);
    P.set(ia, ib, web, app, 0.3);
    P.set(ia, ia, app, data, 1.0);
    P.set(ia, ia, data, think, 1.0);

    P.set(ib, ib, think, web, 1.0);
    P.set(ib, ib, web, app, 0.8);
    P.set(ib, ia, web, app, 0.2);
    P.set(ib, ib, app, data, 1.0);
    P.set(ib, ib, data, think, 1.0);

    P.set(ba, ba, source, web, 1.0);
    P.set(ba, ba, web, app, 0.6);
    P.set(ba, bb, web, app, 0.4);
    P.set(ba, ba, app, data, 1.0);
    P.set(ba, ba, data, sink, 1.0);

    P.set(bb, bb, source, web, 1.0);
    P.set(bb, bb, web, app, 0.7);
    P.set(bb, ba, web, app, 0.3);
    P.set(bb, bb, app, data, 1.0);
    P.set(bb, bb, data, sink, 1.0);

    P.set(ext, ext, source, web, 1.0);
    P.set(ext, ext, web, app, 1.0);
    P.set(ext, ext, app, data, 1.0);
    P.set(ext, ext, data, sink, 1.0);
    model.link(P);

    note("Mixed QN with 3 chains, 5 classes, and class switching:");
    note("- Chain 1 (Closed): InteractiveA and InteractiveB (switch at WebServer)");
    note("- Chain 2 (Open): BatchA and BatchB (switch at WebServer)");
    note("- Chain 3 (Open): ExternalLoad");
    note("");

    const Sn& sn = model.get_struct();

    SolverOpts mva;
    mva.seed = 23000;
    section("MVA");
    print_avg_chain(sn, solve_avg("MVA", model, mva));

    SolverOpts ssa;
    ssa.seed = 23000;
    // The reference also passes cutoff=5 to SSA; neither codebase's simulator
    // reads it (Python carries it into HandlerOptions and no engine consumes it).
    section("SSA");
    print_avg_chain(sn, solve_avg("SSA", model, ssa));
}

/** Five FCFS queues of multiplicity 1 to 5: the sparse LDMX model. */
void mqn_multiserver_fcfs() {
    Net model = multiserver_model(SchedStrategy::FCFS);

    note(
        "This example shows the execution of the solver on a 2-class mixed model with 5 "
        "multi-server FCFS nodes.");

    SolverOpts ctmc;
    ctmc.cutoff = 3.0;
    ctmc.seed = 23000;
    section("CTMC");
    print_avg(solve_avg("CTMC", model, ctmc));

    section("JMT");
    print_avg(model.get_struct(), jmt_avg(model, sim_opts(23000, 100000)));

    SolverOpts ssa;
    ssa.seed = 23000;
    section("SSA");
    print_avg(solve_avg("SSA", model, ssa));

    section("MVA");
    print_avg(solve_avg("MVA", model));

    // TODO(cpp): LDES(model, samples=100000, seed=23000).avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

/** The same five stations under processor sharing, with exact MVA. */
void mqn_multiserver_ps() {
    Net model = multiserver_model(SchedStrategy::PS);

    SolverOpts ctmc;
    ctmc.cutoff = 3.0;
    ctmc.seed = 23000;
    section("CTMC");
    print_avg(solve_avg("CTMC", model, ctmc));

    section("JMT");
    print_avg(model.get_struct(), jmt_avg(model, sim_opts(23000, 100000)));

    SolverOpts mva;
    mva.method = "exact";
    section("MVA");
    print_avg(solve_avg("MVA", model, mva));

    // TODO(cpp): LDES(model, samples=100000, seed=23000).avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

/** Four single-server FCFS queues, 100 closed jobs and an APH arrival stream. */
void mqn_singleserver_fcfs() {
    Net model = singleserver_model(SchedStrategy::FCFS, lang::aph_fit_mean_scv<double>(3.0, 64.0));

    section("JMT");
    print_avg(model.get_struct(), jmt_avg(model, sim_opts(23000, 20000)));

    // The reference leaves its SolverFLD block commented out: the ODEs are too
    // stiff with APH arrivals in a mixed network.

    SolverOpts mva;
    mva.method = "lin";
    section("MVA");
    print_avg(solve_avg("MVA", model, mva));

    SolverOpts nc;
    nc.cutoff = 3.0;
    nc.seed = 23000;
    section("NC");
    print_avg(solve_avg("NC", model, nc));

    // TODO(cpp): LDES(model, keep=False, verbose=True, cutoff=3, seed=23000, samples=20000).avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

/** The same four stations under processor sharing, with Poisson arrivals. */
void mqn_singleserver_ps() {
    Net model = singleserver_model(SchedStrategy::PS, Exp(0.3));

    section("JMT");
    print_avg(model.get_struct(), jmt_avg(model, sim_opts(23000, 100000)));

    SolverOpts ssa;
    ssa.seed = 23000;
    ssa.samples = 100000;
    section("SSA");
    print_avg(solve_avg("SSA", model, ssa));

    SolverOpts fld;
    fld.seed = 23000;
    // The reference also passes cutoff=3 to FLD; a fluid limit enumerates no
    // state space, and neither codebase's fluid options carry the field.
    section("FLD");
    print_avg(solve_avg("FLD", model, fld));

    SolverOpts mva;
    mva.method = "lin";
    section("MVA");
    print_avg(solve_avg("MVA", model, mva));

    SolverOpts nc;
    nc.cutoff = 3.0;
    nc.seed = 23000;
    section("NC");
    print_avg(solve_avg("NC", model, nc));

    // TODO(cpp): LDES(model, keep=False, verbose=True, cutoff=3, seed=23000, samples=100000).avg_table()
    na("LDES", "SolverLDES is the SSJ simulation engine of the JAR and has no C++ port");
}

LINE_EXAMPLE("basic/mixedQN", mqn_basic);
LINE_EXAMPLE("basic/mixedQN", mqn_multichain_cs);
LINE_EXAMPLE("basic/mixedQN", mqn_multiserver_fcfs);
LINE_EXAMPLE("basic/mixedQN", mqn_multiserver_ps);
LINE_EXAMPLE("basic/mixedQN", mqn_singleserver_fcfs);
LINE_EXAMPLE("basic/mixedQN", mqn_singleserver_ps);

}  // namespace examples
}  // namespace line
