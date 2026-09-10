/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/fcRegion/`: finite capacity regions, the cap a set
 * of stations shares.
 *
 * SIX OF THE SEVEN ARE SIMULATION SCRIPTS. The reference solves them with
 * SolverJMT (or SolverLDES for the linear-constraint one) because blocking and
 * dropping across a region are sample-path properties no analytical solver in
 * the suite reproduces; `fcr_lossn` is the exception and runs SolverNC's Erlang
 * fixed point. The models are therefore built here in full -- the region, its
 * global and per-class caps and its drop rules are what the example is about --
 * and the simulation block is refused by name.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "example_util.h"
#include "example_node_table.h"
#include "examples_common.h"

namespace line {
namespace examples {

namespace {

// A node index is its 1-based creation order, so the helpers below fix the
// numbering their callers name a region over.

/** The `Source -> Queue -> Sink` single-class M/M/1 the two comparison scripts share. */
Net mm1_open(const std::string& nm, double lambda, double mu) {
    Net m(nm);
    Source src(m, "Source");
    Queue q(m, "Queue", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    src.set_arrival(c1, Exp(lambda));
    q.set_service(c1, Exp(mu));
    Routing P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    m.link(P);
    return m;
}

/** The two-queue multiclass open network of `fcr_oqn*`, region caps aside. */
Net oqn_twoqueues(const std::string& nm) {
    Net m(nm);
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    OpenClass c2(m, "Class2", 1);
    src.set_arrival(c1, Exp(0.4));
    src.set_arrival(c2, Exp(0.3));
    m.set_service(q1, c1, Exp(1.0));
    m.set_service(q1, c2, Exp(0.9));
    q2.set_service(c1, Exp(1.1));
    q2.set_service(c2, Exp(1.0));
    Routing P;
    P.set(c1, c1, src, q1, 0.5);
    P.set(c1, c1, src, q2, 0.5);
    P.set(c1, c1, q1, q2, 0.3);
    P.set(c1, c1, q1, snk, 0.7);
    P.set(c1, c1, q2, snk, 1.0);
    P.set(c2, c2, src, q1, 0.6);
    P.set(c2, c2, src, q2, 0.4);
    P.set(c2, c2, q1, q2, 0.5);
    P.set(c2, c2, q1, snk, 0.5);
    P.set(c2, c2, q2, snk, 1.0);
    m.link(P);
    return m;
}

/** The Source -> Delay -> Sink two-class model both region-on-a-delay scripts use. */
Net delay_open(const std::string& nm, double lambda1, double lambda2, double mu1, double mu2,
               int prio2) {
    Net m(nm);
    Source src(m, "Source");
    Delay dly(m, "Delay");
    Sink snk(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    OpenClass c2(m, "Class2", prio2);
    src.set_arrival(c1, Exp(lambda1));
    src.set_arrival(c2, Exp(lambda2));
    m.set_service(dly, c1, Exp(mu1));
    m.set_service(dly, c2, Exp(mu2));
    Routing P;
    P.set(c1, c1, src, dly, 1.0);
    P.set(c1, c1, dly, snk, 1.0);
    P.set(c2, c2, src, dly, 1.0);
    P.set(c2, c2, dly, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

void fcr_constraints() {
    Net m("FCR Constraints Demo");
    Source src(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink snk(m, "Sink");
    OpenClass hi(m, "HighPriority", 0);
    OpenClass lo(m, "LowPriority", 1);
    src.set_arrival(hi, Exp(0.3));
    src.set_arrival(lo, Exp(0.5));
    m.set_service(q1, hi, Exp(1.0));
    m.set_service(q1, lo, Exp(0.8));
    q2.set_service(hi, Exp(1.2));
    q2.set_service(lo, Exp(1.0));
    Routing P;
    P.set(hi, hi, src, q1, 1.0);
    P.set(hi, hi, q1, q2, 1.0);
    P.set(hi, hi, q2, snk, 1.0);
    P.set(lo, lo, src, q1, 1.0);
    P.set(lo, lo, q1, q2, 1.0);
    P.set(lo, lo, q2, snk, 1.0);
    m.link(P);

    // Global cap 2 over both queues, per-class cap 2 each, both classes dropping.
    m.add_region(std::vector<std::size_t>{q1, q2}, std::vector<double>{2.0, 2.0}, 2.0,
                 std::vector<DropStrategy>{DropStrategy::DROP, DropStrategy::DROP});

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 100000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_constraints);

void fcr_lincon() {
    Net m = delay_open("FCR LinCon", 0.3, 0.2, 1.0, 0.8, 0);
    const std::size_t dly = 2;  // Source, Delay, Sink
    const std::size_t region =
        m.add_region(std::vector<std::size_t>{dly}, std::vector<double>{1000.0, 1000.0}, 1000.0,
                     std::vector<DropStrategy>{DropStrategy::DROP, DropStrategy::DROP});
    // 2 n_1 + 1 n_2 <= 5 and 1 n_1 + 3 n_2 <= 7, the weighted cross-class coupling.
    Matrix<double> A(2, 2, 0.0);
    A(0, 0) = 2.0;
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    A(1, 1) = 3.0;
    m.set_region_constraint(region, A, std::vector<double>{5.0, 7.0});

    // TODO(cpp): solver = LDES(model, seed=23000, samples=500000); print(solver.avg_table())
    // The engine RUNS; what it does not carry is the constraint. `A n <= b` is a
    // sample-path admission rule, and neither the LDES region model nor the JSIM
    // document has a place to put it, so the run would report an UNCONSTRAINED
    // region under a constrained model's name.
    na("LDES",
       "the general linear admission constraint A n <= b has no wire form: the LDES region model "
       "carries per-class and global job caps and a memory budget, not a coupling matrix, so a "
       "run here would answer the unconstrained region");
}

LINE_EXAMPLE("advanced/fcRegion", fcr_lincon);

void fcr_lossn() {
    Net m = delay_open("FCR Loss Network", 0.3, 0.2, 1.0, 0.8, 1);
    const std::size_t dly = 2;  // Source, Delay, Sink
    m.add_region(std::vector<std::size_t>{dly}, std::vector<double>{3.0, 3.0}, 5.0,
                 std::vector<DropStrategy>{DropStrategy::DROP, DropStrategy::DROP});

    note("Running NC solver (lossn method)...");
    section("NC");
    print_avg(solve_avg("NC", m), "NC Results:");

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 500000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_lossn);

void fcr_mm1kdrop() {
    const double arrival_rate = 0.8, service_rate = 1.0, K = 3.0;

    Net m1 = mm1_open("FCR Dropping", arrival_rate, service_rate);
    const std::size_t q1 = 2;  // Source, Queue, Sink
    m1.add_region(std::vector<std::size_t>{q1}, std::vector<double>{-1.0}, K,
                  std::vector<DropStrategy>{DropStrategy::DROP});

    Net m2 = mm1_open("M/M/1/K", arrival_rate, service_rate);
    m2.set_number_of_servers(2, 1.0);  // the Queue node
    m2.set_capacity(2, K);

    section("JMT");
    print_avg(m1.get_struct(), jmt_avg(m1, sim_opts(23000, 100000)));
    section("JMT");
    print_avg(m2.get_struct(), jmt_avg(m2, sim_opts(23000, 100000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_mm1kdrop);

void fcr_mm1waitq() {
    const double arrival_rate = 0.5, service_rate = 1.0;

    Net m1 = mm1_open("FCR Blocking", arrival_rate, service_rate);
    const std::size_t q1 = 2;  // Source, Queue, Sink
    m1.add_region(std::vector<std::size_t>{q1}, std::vector<double>{-1.0}, 10.0,
                  std::vector<DropStrategy>{DropStrategy::WAITQ});

    Net m2 = mm1_open("M/M/1", arrival_rate, service_rate);

    // BOTH MODELS ARE SOLVED, which is the whole point of the example: the
    // FCR-blocking net must answer what the plain M/M/1 answers, because a
    // region that makes jobs WAIT bounds nothing.
    //
    // THE NODE TABLE, which is the table this reference prints (`avg_node_table`,
    // where `fcr_mm1kdrop` beside it prints `avg_table`). The difference is the
    // point of the example: a WAITQ region holds jobs that are in no station's
    // QLen, and its own row -- the one carrying the parked jobs -- is what makes
    // the population balance. The Sink is not a station either.
    section("JMT");
    print_avg_node(m1.get_struct(), jmt_avg(m1, sim_opts(23000, 100000)));
    section("JMT");
    print_avg_node(m2.get_struct(), jmt_avg(m2, sim_opts(23000, 100000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_mm1waitq);

void fcr_oqndrop() {
    Net m = oqn_twoqueues("FCR Dropping Example");
    const std::size_t q1 = 2, q2 = 3;  // Source, Queue1, Queue2, Sink
    m.add_region(std::vector<std::size_t>{q1, q2}, std::vector<double>{3.0, 2.0}, 4.0,
                 std::vector<DropStrategy>{DropStrategy::DROP, DropStrategy::DROP});

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 50000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_oqndrop);

void fcr_oqnwaitq() {
    Net m = oqn_twoqueues("FCR Blocking Example");
    const std::size_t q1 = 2, q2 = 3;  // Source, Queue1, Queue2, Sink
    m.add_region(std::vector<std::size_t>{q1, q2}, std::vector<double>{5.0, 4.0}, 8.0,
                 std::vector<DropStrategy>{DropStrategy::WAITQ, DropStrategy::WAITQ});

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
}

LINE_EXAMPLE("advanced/fcRegion", fcr_oqnwaitq);

}  // namespace examples
}  // namespace line
