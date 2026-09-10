/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/classSwitching/`: the four ways a job changes class.
 *
 * `cs_implicit` puts the switch on a ClassSwitch node and its matrix;
 * the other three write it straight into the routing block `P{r,s}`, which is
 * the implicit form. `cs_multi_diamond` and `cs_transient_class` are reducible:
 * Classes 2 and 3 carry no arrival and no population and exist only because
 * Class 1 switches into them.
 *
 * Each reference script prints `model.print_routing_matrix()` and then
 * `getAvgChainTable`. Neither has a facade in this port, so both are
 * transcribed here from their references -- `sn_print_routing_matrix.m` and
 * `matlab/src/solvers/@NetworkSolver/getAvg*Chain.m`.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"
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
 * Every column but RespT is a plain sum over the classes of the chain; RespT is
 * weighted by `alpha`, the class share of the chain's visits at the station,
 * because the response times of classes a job moves between do not add.
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

/**
 * `model.print_routing_matrix()`: every positive entry of `sn.rtnodes`.
 *
 * A Sink is skipped as a source of routing and a class whose strategy at the
 * node is DISABLED is skipped too, exactly as `sn_print_routing_matrix.m` does.
 */
void print_routing_matrix(const Sn& sn) {
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;
    if (sn.rtnodes.rows() < I * K) {
        note("No routing matrix available.");
        return;
    }
    for (std::size_t i = 0; i < I; ++i) {
        if (sn.nodes[i].nodetype == lang::NodeType::Sink) continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (r < sn.nodes[i].routing.size() &&
                sn.nodes[i].routing[r] == RoutingStrategy::DISABLED)
                continue;
            for (std::size_t j = 0; j < I; ++j)
                for (std::size_t s = 0; s < K; ++s) {
                    const double p = sn.rtnodes(i * K + r, j * K + s);
                    if (p > 0.0)
                        std::printf("%s [%s] => %s [%s] : Pr=%.6f\n", sn.nodes[i].name.c_str(),
                                    sn.classes[r].name.c_str(), sn.nodes[j].name.c_str(),
                                    sn.classes[s].name.c_str(), p);
                }
        }
    }
}

}  // namespace

/**
 * An M/M/1 whose class switch is an explicit ClassSwitch node: Class1 stays
 * Class1 with probability 0.3 and becomes Class2 with 0.7, and Class2 always
 * becomes Class1.
 */
void cs_implicit() {
    Net model("mm1cs");
    Source source(model, "Source 1");
    Queue queue(model, "Queue 1", SchedStrategy::FCFS);
    Sink sink(model, "Sink 1");

    OpenClass class1(model, "Class1", 0);
    OpenClass class2(model, "Class2", 0);

    source.set_arrival(class1, D::exp_mean(10.0));
    source.set_arrival(class2, D::exp_mean(2.0));
    queue.set_service(class1, D::exp_mean(1.0));
    queue.set_service(class2, D::exp_mean(1.0));

    // The matrix is declared after the classes, so the node can be too.
    Matrix<double> csmatrix(2, 2, 0.0);
    csmatrix(class1 - 1, class1 - 1) = 0.3;
    csmatrix(class1 - 1, class2 - 1) = 0.7;
    csmatrix(class2 - 1, class1 - 1) = 1.0;
    const std::size_t cs = model.add_class_switch("ClassSwitch 1", csmatrix);

    Routing P;
    P.set(class1, class1, source, cs, 1.0);
    P.set(class1, class1, queue, sink, 1.0);
    P.set(class1, class1, cs, queue, 1.0);
    P.set(class2, class2, source, cs, 1.0);
    P.set(class2, class2, queue, sink, 1.0);
    P.set(class2, class2, cs, queue, 1.0);
    model.link(P);

    const Sn& sn = model.get_struct();
    print_routing_matrix(sn);

    section("MVA");
    print_avg_chain(sn, solve_avg("MVA", model));
}

/**
 * Class switching driven by a reducible chain: Class1 arrives, loops at Queue 0
 * and then leaves as Class2 or Class3, each of which has its own exit queue.
 */
void cs_multi_diamond() {
    Net model("mm1cs");
    Source source(model, "Source 1");
    Queue queue0(model, "Queue 0", SchedStrategy::FCFS);
    Queue queue1(model, "Queue 1", SchedStrategy::FCFS);
    Queue queue2(model, "Queue 2", SchedStrategy::FCFS);
    Sink sink(model, "Sink 1");

    OpenClass class1(model, "Class1", 0);
    OpenClass class2(model, "Class2", 0);
    OpenClass class3(model, "Class3", 0);

    source.set_arrival(class1, Exp(1.0));
    queue0.set_service(class1, Exp(10.0));
    queue1.set_service(class2, Exp(20.0));
    queue2.set_service(class3, Exp(30.0));

    Routing P;
    P.set(class1, class1, source, queue0, 1.0);
    P.set(class1, class1, queue0, queue0, 0.2);
    P.set(class1, class2, queue0, queue1, 0.3);
    P.set(class1, class3, queue0, queue2, 0.5);
    P.set(class2, class2, queue1, sink, 1.0);
    P.set(class3, class3, queue2, sink, 1.0);
    model.link(P);

    const Sn& sn = model.get_struct();
    print_routing_matrix(sn);

    section("MVA");
    print_avg_chain(sn, solve_avg("MVA", model));
}

/**
 * The closed counterpart: one job circulating between three Delays, changing
 * class on every hop away from Queue 0 and back to Class1 on the return.
 */
void cs_single_diamond() {
    Net model("mm1cs");
    Delay queue0(model, "Queue 0");
    Delay queue1(model, "Queue 1");
    Delay queue2(model, "Queue 2");

    ClosedClass class1(model, "Class1", 1, queue0);
    ClosedClass class2(model, "Class2", 0, queue0);
    ClosedClass class3(model, "Class3", 0, queue0);

    queue0.set_service(class1, D::exp_mean(1.0));
    queue1.set_service(class2, D::exp_mean(2.0));
    queue2.set_service(class3, D::exp_mean(3.0));

    Routing P;
    P.set(class1, class1, queue0, queue0, 0.2);
    P.set(class1, class2, queue0, queue1, 0.3);
    P.set(class1, class3, queue0, queue2, 0.5);
    P.set(class2, class1, queue1, queue0, 1.0);
    P.set(class3, class1, queue2, queue0, 1.0);
    model.link(P);

    const Sn& sn = model.get_struct();
    print_routing_matrix(sn);

    section("MVA");
    print_avg_chain(sn, solve_avg("MVA", model));
}

/**
 * The transient-class variant: once Class1 has switched, the job stays in
 * Class2 or Class3 forever, so Class1 is transient in the switching chain.
 */
void cs_transient_class() {
    Net model("reducible_cs");
    Delay queue0(model, "Queue 0");
    Delay queue1(model, "Queue 1");
    Delay queue2(model, "Queue 2");

    ClosedClass class1(model, "Class1", 1, queue0);
    ClosedClass class2(model, "Class2", 0, queue0);
    ClosedClass class3(model, "Class3", 0, queue0);

    queue0.set_service(class1, D::exp_mean(1.0));
    queue0.set_service(class2, D::exp_mean(1.0));
    queue0.set_service(class3, D::exp_mean(1.0));
    queue1.set_service(class2, D::exp_mean(1.0));
    queue2.set_service(class3, D::exp_mean(1.0));

    Routing P;
    P.set(class1, class1, queue0, queue0, 0.2);
    P.set(class1, class2, queue0, queue1, 0.3);
    P.set(class1, class3, queue0, queue2, 0.5);
    P.set(class2, class2, queue0, queue1, 1.0);
    P.set(class2, class2, queue1, queue0, 1.0);
    P.set(class3, class3, queue0, queue2, 1.0);
    P.set(class3, class3, queue2, queue0, 1.0);
    model.link(P);

    const Sn& sn = model.get_struct();
    print_routing_matrix(sn);

    section("MVA");
    print_avg_chain(sn, solve_avg("MVA", model));
}

LINE_EXAMPLE("basic/classSwitching", cs_implicit);
LINE_EXAMPLE("basic/classSwitching", cs_multi_diamond);
LINE_EXAMPLE("basic/classSwitching", cs_single_diamond);
LINE_EXAMPLE("basic/classSwitching", cs_transient_class);

}  // namespace examples
}  // namespace line
