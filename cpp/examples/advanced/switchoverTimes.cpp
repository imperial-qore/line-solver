/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/switchoverTimes/`: the time a server needs to
 * switch from serving one class to another.
 *
 * A SWITCHOVER DECLARED AT AN ORDINARY QUEUE IS CARRIED, NOT CONSUMED -- in
 * every codebase, not only here. JMT's `Server` has no switchover parameter, so
 * `writeJSIM` warns and drops the (from, to) times; no analytical solver reads
 * them either, since the MVA polling analyzer and the state machinery reach a
 * switchover through a POLLING station's buffers. The model below is therefore
 * the M/M/1 the same two classes make without it, and the reference reports
 * exactly that. Declaring the times still matters: it is what the example is
 * about, and it is what the document round-trip carries.
 *
 * The queue is deliberately SATURATED (rho = 2.53 from class 1 alone at 2.0),
 * so the table below is a transient of a growing queue rather than a steady
 * state, and it is reproducible only because the seed and the sample count are
 * the reference's.
 */

#include "example_util.h"
#include "examples_common.h"

namespace line {
namespace examples {

/** `switchover_basic.py`: pairwise switchover times at an FCFS queue. */
void switchover_basic() {
    Net m("M[2]/M[2]/1-Gated");
    Source src(m, "mySource");
    Queue q(m, "myQueue", SchedStrategy::FCFS);
    Sink snk(m, "mySink");

    OpenClass c1(m, "myClass1", 0);
    src.set_arrival(c1, Exp(0.2));
    q.set_service(c1, Exp(0.1));

    OpenClass c2(m, "myClass2", 0);
    src.set_arrival(c2, Exp(0.8));
    q.set_service(c2, Exp(1.5));

    // Switching from class 1 to class 2 costs an Exp(1); the other way round an
    // Erlang of two phases of rate 1.
    q.set_switchover(c1, c2, Exp(1.0));
    q.set_switchover(c2, c1, Erlang(1.0, 2));

    Routing P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, snk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, snk, 1.0);
    m.link(P);

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 10000)));
}

LINE_EXAMPLE("advanced/switchoverTimes", switchover_basic);

}  // namespace examples
}  // namespace line
