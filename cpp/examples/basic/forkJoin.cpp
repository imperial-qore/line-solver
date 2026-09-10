/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/forkJoin/`: the eighteen fork-join models.
 *
 * Almost every one of them exists to put SolverMVA's fork-join transform beside
 * SolverJMT's simulation of the same topology, and both blocks run: `jmt_avg`
 * of example_util.h drives the same `common/JMT.jar` the reference drives, at
 * the reference's own seed, so the comparison the example was written for is
 * the comparison it prints. A host without the jar or a JVM makes the JMT call
 * THROW by name, which is a different statement from a missing port.
 *
 * BOTH FORK-JOIN ARMS RUN. `fj_basic_closed` pins `fork_join='ht'`, the
 * Heidelberger-Trivedi response-time method of `solvers/mva/fj_ht.h`; every
 * other example takes the MMT transform of `solvers/mva/fj_mmt.h`, which is the
 * reference's `default`/`mmt`/`fjt` arm.
 */

#include <cstddef>

#include "example_util.h"
#include "examples_common.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace line {
namespace examples {

namespace {

mva::AvgResult<double> mva_run(Net& m, const std::string& method = "default",
                               const std::string& fork_join = "default") {
    mva::MvaOptions opt;
    opt.method = method;
    opt.fork_join = fork_join;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

/** Asymmetric branches: Queue1 alone against the Queue2 -> Queue3 chain. */
Net fj_asymm_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    Queue queue3(m, "Queue3", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(2.0));
    queue3.set_service(c1, Exp(1.0));
    delay.set_service(c1, Exp(0.5));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, queue3, 1.0);
    P.set(c1, c1, queue3, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/** Two symmetric PS branches, Queue1 and Queue2, closed through the Delay. */
Net fj_basic_closed_model() {
    Net m("model");
    Delay delay(m, "Delay");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 5.0, delay);

    delay.set_service(c1, Exp(1.0));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}


/**
 * A fork whose degree is NOT one fixed number, in the four ways it can vary.
 *
 *   0 fixed   the classic fork, one task on each of two links (the baseline)
 *   1 vector  three tasks towards Queue2 and one towards Queue1
 *   2 random  one or three tasks per link, each with probability 1/2
 *   3 prob    the branch towards Queue2 fires only half the time
 *
 * The overrides are set AFTER `link()`, because a per-destination override
 * names a link and the routing is what decides which links exist.
 */
Net fj_variable_fanout_model(int mode) {
    Net m("model");
    Delay delay(m, "Delay");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 4.0, delay);

    delay.set_service(c1, Exp(1.0));
    queue1.set_service(c1, Exp(2.0));
    queue2.set_service(c1, Exp(2.0));

    // Declared BEFORE the routing, exactly as the MATLAB and Python twins do.
    // The overrides are recorded and replayed by `link()`, so either order
    // builds the same model.
    if (mode == 1) {
        m.set_fork_tasks_per_link(fork, c1, 3.0, queue2);
    } else if (mode == 2) {
        std::vector<double> pmf(2, 0.5);
        std::vector<double> pts;
        pts.push_back(1.0);
        pts.push_back(3.0);
        m.set_fork_tasks_per_link_dist(fork, c1, D::discrete_sampler(pmf, pts));
    } else if (mode == 3) {
        // an uncertain branch needs a Join that does not wait for it
        m.set_join_strategy(join, lang::JoinStrategy::PARTIAL, 1.0);
        m.set_fork_branch_probability(fork, c1, queue2, 0.5);
    }

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * Two forks, one enclosing the other on class2 only.
 *
 * class1 takes Fork1 -{Queue1, Queue2}- Join1 and closes at the Delay. class2
 * enters Fork1_1 (two tasks per link) first, then the same Fork1 pair, and its
 * Join1 feeds Join1_1, so on class2 Fork1 is nested inside Fork1_1.
 */
Net fj_basic_nesting_model() {
    Net m("model");
    Delay delay(m, "Delay");
    const std::size_t fork1 = m.add_fork("Fork1", 1.0);
    const std::size_t fork11 = m.add_fork("Fork1_1", 2.0);
    const std::size_t join1 = m.add_join("Join1", fork1);
    const std::size_t join11 = m.add_join("Join1_1", fork11);
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);

    ClosedClass c1(m, "class1", 5.0, delay, 0);
    ClosedClass c2(m, "class2", 2.0, delay, 0);

    delay.set_service(c1, Exp(0.25));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(0.75));
    delay.set_service(c2, Exp(0.25));
    queue1.set_service(c2, Exp(2.0));
    queue2.set_service(c2, Exp(2.0));

    Routing P;
    P.set(c1, c1, delay, fork1, 1.0);
    P.set(c1, c1, fork1, queue1, 1.0);
    P.set(c1, c1, fork1, queue2, 1.0);
    P.set(c1, c1, queue1, join1, 1.0);
    P.set(c1, c1, queue2, join1, 1.0);
    P.set(c1, c1, join1, delay, 1.0);
    P.set(c2, c2, delay, fork11, 1.0);
    P.set(c2, c2, fork11, fork1, 1.0);
    P.set(c2, c2, fork1, queue1, 1.0);
    P.set(c2, c2, fork1, queue2, 1.0);
    P.set(c2, c2, queue1, join1, 1.0);
    P.set(c2, c2, queue2, join1, 1.0);
    P.set(c2, c2, join1, join11, 1.0);
    P.set(c2, c2, join11, delay, 1.0);
    m.link(P);
    return m;
}

/** Source -> Fork -{Queue1, Queue2}- Join -> Sink, the open two-branch shape. */
Net fj_basic_open_model() {
    Net m("model");
    Source source(m, "Source");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1");

    source.set_arrival(c1, Exp(0.05));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(2.0));

    Routing P;
    P.set(c1, c1, source, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, sink, 1.0);
    m.link(P);
    return m;
}

/**
 * Two branches, each a serial chain of its own: Queue1 -> Queue4 -> Queue5 and
 * Queue2 -> Queue3, both ending at the Join.
 */
Net fj_complex_serial_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    Queue queue3(m, "Queue3", SchedStrategy::FCFS);
    Queue queue4(m, "Queue4", SchedStrategy::FCFS);
    Queue queue5(m, "Queue5", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(2.0));
    queue3.set_service(c1, Exp(1.0));
    queue4.set_service(c1, Exp(3.0));
    queue5.set_service(c1, Exp(0.8));
    delay.set_service(c1, Exp(0.5));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, queue4, 1.0);
    P.set(c1, c1, queue4, queue5, 1.0);
    P.set(c1, c1, queue5, join, 1.0);
    P.set(c1, c1, queue2, queue3, 1.0);
    P.set(c1, c1, queue3, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * The same fork traversed twice per job: class1 forks, joins, switches to
 * class2 at the Join and forks again before leaving.
 */
Net fj_cs_multi_visits_model() {
    Net m("model");
    Source source(m, "Source");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1");
    OpenClass c2(m, "class2");

    source.set_arrival(c1, Exp(0.1));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    queue1.set_service(c2, Exp(1.0));
    queue2.set_service(c2, Exp(1.0));

    Routing P;
    P.set(c1, c1, source, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c2, join, fork, 1.0);
    P.set(c2, c2, fork, queue1, 1.0);
    P.set(c2, c2, fork, queue2, 1.0);
    P.set(c2, c2, queue1, join, 1.0);
    P.set(c2, c2, queue2, join, 1.0);
    P.set(c2, c2, join, sink, 1.0);
    m.link(P);
    return m;
}

/** class2 switches to class1 ON the fork's output arcs, so only class1 forks. */
Net fj_cs_postfork_model() {
    Net m("model");
    Delay delay(m, "Delay");
    const std::size_t fork1 = m.add_fork("Fork1");
    const std::size_t join1 = m.add_join("Join1", fork1);
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);

    ClosedClass c1(m, "class1", 1.0, delay, 0);
    ClosedClass c2(m, "class2", 1.0, delay, 0);

    delay.set_service(c1, Exp(0.25));
    queue1.set_service(c1, Exp(2.0));
    queue2.set_service(c1, Exp(2.0));
    delay.set_service(c2, Exp(0.25));
    queue1.set_service(c2, Exp(2.0));
    queue2.set_service(c2, Exp(2.0));

    Routing P;
    P.set(c1, c1, delay, fork1, 1.0);
    P.set(c1, c1, fork1, queue1, 1.0);
    P.set(c1, c1, fork1, queue2, 1.0);
    P.set(c1, c1, queue1, join1, 1.0);
    P.set(c1, c1, queue2, join1, 1.0);
    P.set(c1, c1, join1, delay, 1.0);
    P.set(c2, c2, delay, fork1, 1.0);
    P.set(c2, c1, fork1, queue1, 1.0);
    P.set(c2, c1, fork1, queue2, 1.0);
    m.link(P);
    return m;
}

/**
 * The class switch happens BEFORE the fork: class1 becomes class2 at Delay1,
 * class2 becomes class1 again at Delay2, and only class1 reaches the Fork.
 */
Net fj_cs_prefork_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Delay delay2(m, "Delay2");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);
    ClosedClass c2(m, "class2", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    delay.set_service(c1, Exp(0.5));
    delay.set_service(c2, Exp(0.5));
    delay2.set_service(c2, Exp(2.0));

    Routing P;
    P.set(c1, c2, delay, delay2, 1.0);
    P.set(c2, c1, delay2, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * Fork nested inside a branch of Fork: the outer branch Queue1 forks again into
 * Queue3 and Queue4, whose Join2 feeds the outer Join; Queue2 is the other
 * outer branch and reaches the outer Join directly.
 */
Net fj_deep_nesting_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Queue queue3(m, "Queue3", SchedStrategy::FCFS);
    Queue queue4(m, "Queue4", SchedStrategy::FCFS);
    const std::size_t fork2 = m.add_fork("Fork2");
    const std::size_t join2 = m.add_join("Join2", fork2);

    ClosedClass c1(m, "class1", 1.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    delay.set_service(c1, Exp(0.5));
    queue3.set_service(c1, Exp(2.0));
    queue4.set_service(c1, Exp(2.0));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, fork2, 1.0);
    P.set(c1, c1, fork2, queue3, 1.0);
    P.set(c1, c1, fork2, queue4, 1.0);
    P.set(c1, c1, queue3, join2, 1.0);
    P.set(c1, c1, queue4, join2, 1.0);
    P.set(c1, c1, join2, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/** Two delays in series ahead of the fork: Delay1 -> Delay2 -> Fork. */
Net fj_delays_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Delay delay2(m, "Delay2");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    delay.set_service(c1, Exp(0.5));
    delay2.set_service(c1, Exp(2.0));

    Routing P;
    P.set(c1, c1, delay, delay2, 1.0);
    P.set(c1, c1, delay2, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * One fork traversed by an open AND a closed class at once.
 *
 * Both classes run PreFork -> Fork -{Branch1, Branch2}- Join -> PostJoin; the
 * open one enters at the Source and leaves at the Sink, the closed one
 * circulates through the Client delay.
 */
Net fj_mixed_openclosed_model() {
    Net m("ForkJoinOpenClosed");
    Source source(m, "Source");
    Delay client(m, "Client");
    Queue prefork(m, "PreFork", SchedStrategy::FCFS);
    Queue branch1(m, "Branch1", SchedStrategy::FCFS);
    Queue branch2(m, "Branch2", SchedStrategy::FCFS);
    Queue postjoin(m, "PostJoin", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Sink sink(m, "Sink");

    OpenClass oclass(m, "Open");
    ClosedClass cclass(m, "Closed", 1.0, client);

    source.set_arrival(oclass, Exp(0.1));
    client.set_service(cclass, D::exp_mean(1.0));
    client.set_service(oclass, Disabled());
    const std::size_t station[4] = {prefork, branch1, branch2, postjoin};
    const double mean[4] = {0.2, 0.3, 0.4, 0.1};
    for (std::size_t i = 0; i < 4; ++i) {
        m.set_service(station[i], oclass, D::exp_mean(mean[i]));
        m.set_service(station[i], cclass, D::exp_mean(mean[i]));
    }

    Routing P;
    const std::size_t entry[2] = {source, client};
    const std::size_t exit_node[2] = {sink, client};
    const std::size_t cls[2] = {oclass, cclass};
    for (std::size_t k = 0; k < 2; ++k) {
        P.set(cls[k], cls[k], entry[k], prefork, 1.0);
        P.set(cls[k], cls[k], prefork, fork, 1.0);
        P.set(cls[k], cls[k], fork, branch1, 1.0);
        P.set(cls[k], cls[k], fork, branch2, 1.0);
        P.set(cls[k], cls[k], branch1, join, 1.0);
        P.set(cls[k], cls[k], branch2, join, 1.0);
        P.set(cls[k], cls[k], join, postjoin, 1.0);
        P.set(cls[k], cls[k], postjoin, exit_node[k], 1.0);
    }
    m.link(P);
    return m;
}

/** A Fork with no Join: three branches that each leave through the Sink. */
Net fj_nojoin_model() {
    Net m("model");
    Source source(m, "Source");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    Queue queue3(m, "Queue3", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1");

    source.set_arrival(c1, Exp(0.5));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(2.0));
    queue3.set_service(c1, Exp(3.0));

    Routing P;
    P.set(c1, c1, source, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, fork, queue3, 1.0);
    P.set(c1, c1, queue1, sink, 1.0);
    P.set(c1, c1, queue2, sink, 1.0);
    P.set(c1, c1, queue3, sink, 1.0);
    m.link(P);
    return m;
}

/** Two closed classes over the SAME two branches, with different think times. */
Net fj_route_overlap_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);
    ClosedClass c2(m, "class2", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(2.0));
    delay.set_service(c1, Exp(0.5));
    queue1.set_service(c2, Exp(1.0));
    queue2.set_service(c2, Exp(2.0));
    delay.set_service(c2, Exp(0.2));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    P.set(c2, c2, delay, fork, 1.0);
    P.set(c2, c2, fork, queue1, 1.0);
    P.set(c2, c2, fork, queue2, 1.0);
    P.set(c2, c2, queue1, join, 1.0);
    P.set(c2, c2, queue2, join, 1.0);
    P.set(c2, c2, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * Two fork-join stages in SERIES, neither inside the other: Fork -{Queue1,
 * Queue2}- Join -> Fork2 -{Queue3, Queue4}- Join2, closed through the Delay.
 */
Net fj_serialfjs_closed_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Queue queue3(m, "Queue3", SchedStrategy::PS);
    Queue queue4(m, "Queue4", SchedStrategy::PS);
    const std::size_t fork2 = m.add_fork("Fork2");
    const std::size_t join2 = m.add_join("Join2", fork2);

    ClosedClass c1(m, "class1", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    delay.set_service(c1, Exp(0.5));
    queue3.set_service(c1, Exp(1.0));
    queue4.set_service(c1, Exp(1.0));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, fork2, 1.0);
    P.set(c1, c1, fork2, queue3, 1.0);
    P.set(c1, c1, fork2, queue4, 1.0);
    P.set(c1, c1, queue3, join2, 1.0);
    P.set(c1, c1, queue4, join2, 1.0);
    P.set(c1, c1, join2, delay, 1.0);
    m.link(P);
    return m;
}

/** The same two stages in series, open: Source -> Fork1 ... Join2 -> Sink. */
Net fj_serialfjs_open_model() {
    Net m("model");
    Source source(m, "Source");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork1 = m.add_fork("Fork1");
    const std::size_t join1 = m.add_join("Join1", fork1);
    Queue queue3(m, "Queue3", SchedStrategy::FCFS);
    Queue queue4(m, "Queue4", SchedStrategy::FCFS);
    const std::size_t fork2 = m.add_fork("Fork2");
    const std::size_t join2 = m.add_join("Join2", fork2);
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1");

    source.set_arrival(c1, Exp(0.4));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(1.0));
    queue3.set_service(c1, Exp(1.0));
    queue4.set_service(c1, Exp(1.0));

    Routing P;
    P.set(c1, c1, source, fork1, 1.0);
    P.set(c1, c1, fork1, queue1, 1.0);
    P.set(c1, c1, fork1, queue2, 1.0);
    P.set(c1, c1, queue1, join1, 1.0);
    P.set(c1, c1, queue2, join1, 1.0);
    P.set(c1, c1, join1, fork2, 1.0);
    P.set(c1, c1, fork2, queue3, 1.0);
    P.set(c1, c1, fork2, queue4, 1.0);
    P.set(c1, c1, queue3, join2, 1.0);
    P.set(c1, c1, queue4, join2, 1.0);
    P.set(c1, c1, join2, sink, 1.0);
    m.link(P);
    return m;
}

/**
 * Three fork outputs over two branches: Queue1 alone, and Queue2 -> Queue3 in
 * series, both reaching the Join; two closed classes traverse them.
 */
Net fj_threebranches_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    Queue queue3(m, "Queue3", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 10.0, delay, 0);
    ClosedClass c2(m, "class2", 10.0, delay, 0);

    queue1.set_service(c1, Exp(1.5));
    queue2.set_service(c1, Exp(1.1));
    queue3.set_service(c1, Exp(2.5));
    delay.set_service(c1, Exp(0.5));
    queue1.set_service(c2, Exp(2.8));
    queue2.set_service(c2, Exp(3.0));
    queue3.set_service(c2, Exp(1.0));
    delay.set_service(c2, Exp(0.8));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue2, queue3, 1.0);
    P.set(c1, c1, queue3, join, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    P.set(c2, c2, delay, fork, 1.0);
    P.set(c2, c2, fork, queue1, 1.0);
    P.set(c2, c2, fork, queue2, 1.0);
    P.set(c2, c2, queue2, queue3, 1.0);
    P.set(c2, c2, queue3, join, 1.0);
    P.set(c2, c2, queue1, join, 1.0);
    P.set(c2, c2, join, delay, 1.0);
    m.link(P);
    return m;
}

/** One circulating job over two FCFS branches: X = 1/(Z + E[max(S1,S2)]). */
Net fj_tiny_closed_model() {
    Net m("model");
    Delay delay(m, "Delay1");
    Queue queue1(m, "Queue1", SchedStrategy::FCFS);
    Queue queue2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork1");
    const std::size_t join = m.add_join("Join1", fork);

    ClosedClass c1(m, "Class1", 1.0, delay);

    delay.set_service(c1, Exp(1.0));
    queue1.set_service(c1, Exp(2.0));
    queue2.set_service(c1, Exp(3.0));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * Two open classes over the same two branches, with a fan-out of two tasks per
 * link; class2 is served immediately at Queue1, so only Queue2 delays it.
 */
Net fj_twoclasses_forked_model() {
    Net m("model");
    Source source(m, "Source");
    Queue queue1(m, "Queue1", SchedStrategy::PS);
    Queue queue2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork", 2.0);
    const std::size_t join = m.add_join("Join", fork);
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1");
    OpenClass c2(m, "class2");

    source.set_arrival(c1, Exp(0.25));
    queue1.set_service(c1, Exp(1.0));
    queue2.set_service(c1, Exp(0.75));
    source.set_arrival(c2, Exp(0.25));
    queue1.set_service(c2, Immediate());
    queue2.set_service(c2, Exp(2.0));

    Routing P;
    P.set(c1, c1, source, fork, 1.0);
    P.set(c1, c1, fork, queue1, 1.0);
    P.set(c1, c1, fork, queue2, 1.0);
    P.set(c1, c1, queue1, join, 1.0);
    P.set(c1, c1, queue2, join, 1.0);
    P.set(c1, c1, join, sink, 1.0);
    P.set(c2, c2, source, fork, 1.0);
    P.set(c2, c2, fork, queue1, 1.0);
    P.set(c2, c2, fork, queue2, 1.0);
    P.set(c2, c2, queue1, join, 1.0);
    P.set(c2, c2, queue2, join, 1.0);
    P.set(c2, c2, join, sink, 1.0);
    m.link(P);
    return m;
}

// ---------------------------------------------------------------------------
// Examples
// ---------------------------------------------------------------------------

void fj_asymm() {
    Net m = fj_asymm_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

/**
 * `MVA(model, method='amva', fork_join='ht')`.
 *
 * The Heidelberger-Trivedi response-time method is the second fork-join arm the
 * reference offers beside the MMT transform, and the example pins it together
 * with `amva`, whose fixed point the reference's own comment calls the smoother
 * of the two on this model.
 */
void fj_basic_closed() {
    Net m = fj_basic_closed_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m, "amva", "ht"));
}

void fj_basic_nesting() {
    Net m = fj_basic_nesting_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_basic_open() {
    Net m = fj_basic_open_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    // `LDES(model, seed=23000)`, the reference's third solver. It was left out,
    // so the golden's whole LDES row went unmeasured; the engine is driven the
    // same way the other simulation blocks in this corpus drive it.
    section("LDES");
    print_avg(m.get_struct(), ldes_avg(m, sim_opts(23000)));
}

void fj_complex_serial() {
    Net m = fj_complex_serial_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_cs_multi_visits() {
    Net m = fj_cs_multi_visits_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void fj_cs_postfork() {
    Net m = fj_cs_postfork_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void fj_cs_prefork() {
    Net m = fj_cs_prefork_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_deep_nesting() {
    Net m = fj_deep_nesting_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_delays() {
    Net m = fj_delays_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_mixed_openclosed() {
    Net m = fj_mixed_openclosed_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, 200000)));
}

/**
 * A Fork with NO Join, which this port cannot yet build.
 *
 * `network_builder.h`'s link-time check refuses an unjoined Fork outright, and
 * `fj_driver.h` refuses one again in the transform, both on the ground that
 * MATLAB's `sortForks` recurses without end on such a model. That holds for a
 * CLOSED model; this one is OPEN, and native Python solves it today (MVA gives
 * Queue1 QLen 1, Queue2 0.33314, Queue3 0.19913, agreeing with the JMT run
 * stored in fj_nojoin.ipynb), so the blanket refusal is a port gap and not a
 * reference limitation. The model is built exactly as the reference builds it
 * and the library's own message is what the driver reports; nothing here works
 * around it, because a worked-around gap is one nobody fixes.
 */
void fj_nojoin() {
    Net m = fj_nojoin_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void fj_route_overlap() {
    Net m = fj_route_overlap_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_serialfjs_closed() {
    Net m = fj_serialfjs_closed_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_serialfjs_open() {
    Net m = fj_serialfjs_open_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

void fj_threebranches() {
    Net m = fj_threebranches_model();
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}

/** The one model here small enough for the exact tag-augmented state space. */
void fj_tiny_closed() {
    Net m = fj_tiny_closed_model();
    section("CTMC");
    ctmc::CtmcOptions copt;
    print_avg(m.get_struct(), ctmc::solver_ctmc_run_analyzer(m.get_struct(), copt));

    section("SSA");
    ssa::SsaOptions sopt;
    sopt.seed = 23000;
    sopt.samples = 50000;
    print_avg_sim(m.get_struct(), ssa::solver_ssa(m.get_struct(), sopt));

    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
}

void fj_twoclasses_forked() {
    Net m = fj_twoclasses_forked_model();
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    section("MVA");
    print_avg(m.get_struct(), mva_run(m));
}


void fj_variable_fanout() {
    static const char* names[4] = {"fixed", "vector", "random", "prob"};
    for (int mode = 0; mode < 4; ++mode) {
        Net m = fj_variable_fanout_model(mode);
        section(std::string("mode = ") + names[mode] + " / JMT");
        print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000)));
    }
}

LINE_EXAMPLE("basic/forkJoin", fj_asymm);
LINE_EXAMPLE("basic/forkJoin", fj_basic_closed);
LINE_EXAMPLE("basic/forkJoin", fj_basic_nesting);
LINE_EXAMPLE("basic/forkJoin", fj_basic_open);
LINE_EXAMPLE("basic/forkJoin", fj_complex_serial);
LINE_EXAMPLE("basic/forkJoin", fj_cs_multi_visits);
LINE_EXAMPLE("basic/forkJoin", fj_cs_postfork);
LINE_EXAMPLE("basic/forkJoin", fj_cs_prefork);
LINE_EXAMPLE("basic/forkJoin", fj_deep_nesting);
LINE_EXAMPLE("basic/forkJoin", fj_delays);
LINE_EXAMPLE("basic/forkJoin", fj_mixed_openclosed);
LINE_EXAMPLE("basic/forkJoin", fj_nojoin);
LINE_EXAMPLE("basic/forkJoin", fj_route_overlap);
LINE_EXAMPLE("basic/forkJoin", fj_serialfjs_closed);
LINE_EXAMPLE("basic/forkJoin", fj_serialfjs_open);
LINE_EXAMPLE("basic/forkJoin", fj_threebranches);
LINE_EXAMPLE("basic/forkJoin", fj_tiny_closed);
LINE_EXAMPLE("basic/forkJoin", fj_twoclasses_forked);
LINE_EXAMPLE("basic/forkJoin", fj_variable_fanout);

}  // namespace

}  // namespace examples
}  // namespace line
