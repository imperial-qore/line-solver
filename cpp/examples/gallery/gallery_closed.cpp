/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The closed half of the gallery, plus the caches, the fork-join pair, the
 * finite capacity region and the random environment.
 *
 * FOUR OF THESE FACTORIES DRAW THEIR PARAMETERS. Three of them --
 * `gallery_cqn`, `gallery_cqn_multiclass` and `gallery_repairmen` -- draw from
 * `random.random()` after `random.seed(seed)`, and they draw a NUMBER OF VALUES
 * THAT DEPENDS ON THEIR ARGUMENTS, so `py_random.h` replays that stream bit for
 * bit rather than freezing the default model's draws as literals. They keep the
 * reference's `seed` argument for the same reason: the seed is part of the
 * model. The realized values at the default arguments were checked against the
 * reference (2026-07-31): population 3 and means 1.2090597453293104,
 * 2.282549310366711 for `gallery_cqn`; population 8 and mean 1.301152946901345
 * for `gallery_repairmen`; means 31, 21, 14, 12 for `gallery_cqn_multiclass`.
 *
 * EACH DRAW IS SEQUENCED INTO ITS OWN VARIABLE. Two `rng.random()` calls in one
 * expression are evaluated in an unspecified order in C++, so the model would
 * get the draws transposed on some compilers and not others.
 *
 * `gallery_qn_random` is the fourth and is the one exception: it seeds NUMPY as
 * well and runs `NetworkGenerator`, whose draws interleave two streams, so its
 * realized topology is transcribed as literals instead. It takes no arguments,
 * so there is no second model to get wrong.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "gallery.h"
#include "py_random.h"

namespace line {
namespace examples {

namespace {

/** `DiscreteSampler([1/n] * n)`: the uniform popularity of an n-item cache. */
std::vector<double> uniform_pread(std::size_t n) {
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

}  // namespace

// ---------------------------------------------------------------------------
// Closed networks
// ---------------------------------------------------------------------------

/**
 * Single-class CQN: M PS queues in a cycle, optionally behind a Delay.
 *
 * `random.seed(23000)`, then M + 1 draws: the population is
 * `round(random.random() * M + 2)` and each queue's mean is
 * `random.random() + i + 1`. The Delay's mean is the literal 2.0 the reference
 * writes and consumes no draw.
 */
Net gallery_cqn(std::size_t M, bool use_delay, unsigned long long seed) {
    PyRandom rng(seed);
    Net m("Single-class CQN");
    std::vector<std::size_t> station;
    for (std::size_t i = 0; i < M; ++i)
        station.push_back(m.add_queue("Queue" + std::to_string(i + 1), SchedStrategy::PS));
    if (use_delay) station.push_back(m.add_delay("Delay1"));

    const double u_pop = rng.random();
    ClosedClass c1(m, "Class1", py_round(u_pop * static_cast<double>(M) + 2.0), station[0], 0);
    for (std::size_t i = 0; i < M; ++i) {
        const double u = rng.random();
        m.set_service(station[i], c1, D::exp_mean(u + static_cast<double>(i) + 1.0));
    }
    if (use_delay) m.set_service(station.back(), c1, D::exp_mean(2.0));

    Routing P;
    cyclic(P, c1, station);
    m.link(P);
    return m;
}

/**
 * Multi-class CQN: m PS queues and a Delay, r classes of five jobs each.
 *
 * `random.seed(23000)`, then per class `round(50 * random.random())` at each
 * queue and `round(100 * random.random())` at the Delay, in that order. The
 * draws are interleaved with the classes, so the stream has to be replayed
 * rather than tabulated.
 */
Net gallery_cqn_multiclass(std::size_t m_queues, std::size_t r, bool wantdelay,
                           unsigned long long seed) {
    PyRandom rng(seed);
    Net m("Multi-class CQN");
    std::vector<std::size_t> node;
    for (std::size_t i = 0; i < m_queues; ++i)
        node.push_back(m.add_queue("Queue " + std::to_string(i + 1), SchedStrategy::PS));
    if (wantdelay) node.push_back(m.add_delay("Delay 1"));

    std::vector<std::size_t> jobclass;
    for (std::size_t s = 0; s < r; ++s)
        jobclass.push_back(m.add_closed_class("Class" + std::to_string(s + 1), 5, node[0], 0));

    for (std::size_t s = 0; s < r; ++s) {
        for (std::size_t i = 0; i < m_queues; ++i) {
            const double u = rng.random();
            m.set_service(node[i], jobclass[s], D::exp_mean(py_round(50.0 * u)));
        }
        if (wantdelay) {
            const double u = rng.random();
            m.set_service(node.back(), jobclass[s], D::exp_mean(py_round(100.0 * u)));
        }
    }

    Routing P;
    for (std::size_t s = 0; s < r; ++s) cyclic(P, jobclass[s], node);
    m.link(P);
    return m;
}

/**
 * The finite-repairmen model: a PS repair queue and a Delay holding the
 * up time.
 *
 * `random.seed(2300)`, then `round(random.random() * 10 * M + 3)` for the
 * population and `random.random() + 1` for the repair mean, with the
 * reference's M fixed at 1.
 */
Net gallery_repairmen(double nservers, unsigned long long seed) {
    const double M = 1.0;
    PyRandom rng(seed);
    Net m("Finite repairmen CQN");
    Queue queue(m, "Queue1", SchedStrategy::PS);
    queue.set_number_of_servers(nservers);
    Delay delay(m, "Delay1");

    const double u_pop = rng.random();
    ClosedClass c1(m, "Class1", py_round(u_pop * 10.0 * M + 3.0), queue, 0);
    const double u_svc = rng.random();
    queue.set_service(c1, D::exp_mean(u_svc + 1.0));
    delay.set_service(c1, D::exp_mean(2.0));

    Routing P;
    cyclic(P, c1, {queue, delay});
    m.link(P);
    return m;
}

/**
 * `NetworkGenerator(...).generate(3, 1, 0, 2)` at seed 23000: three FCFS
 * queues, one Delay and two closed classes on a cyclic topology.
 *
 * The generator's ranges are all degenerate, so every service rate is exactly
 * 1 and the only drawn quantities are the two populations (14 and 11) and the
 * reference station (queue2 for both classes). Extracted from the reference's
 * `getStruct()` -- rates, njobs, refstat and rt.
 */
Net gallery_qn_random() {
    Net m("nw");
    Queue q1(m, "queue1", SchedStrategy::FCFS);
    Queue q2(m, "queue2", SchedStrategy::FCFS);
    Queue q3(m, "queue3", SchedStrategy::FCFS);
    Delay d1(m, "delay1");

    ClosedClass c1(m, "CClass1", 14, q2, 0);
    ClosedClass c2(m, "CClass2", 11, q2, 0);
    const std::size_t nodes[4] = {q1, q2, q3, d1};
    for (std::size_t i = 0; i < 4; ++i) {
        m.set_service(nodes[i], c1, Exp(1.0));
        m.set_service(nodes[i], c2, Exp(1.0));
    }

    Routing P;
    cyclic(P, c1, {q1, q2, q3, d1});
    cyclic(P, c2, {q1, q2, q3, d1});
    m.link(P);
    return m;
}

// ---------------------------------------------------------------------------
// Caches
// ---------------------------------------------------------------------------

/** A closed cache: one job circulating between a Delay and an LRU cache. */
Net gallery_cache_lru() {
    const std::size_t n = 5, mslots = 2;
    Net m("Cache-LRU");
    Delay delay(m, "Delay");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{static_cast<int>(mslots)};
    ch.replacestrat = lang::ReplacementStrategy::LRU;
    // JobClass reads uniformly; the hit and miss classes never read.
    ch.pread = std::vector<std::vector<double> >{uniform_pread(n), {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = m.add_cache("Cache", ch);

    ClosedClass job(m, "JobClass", 1, delay, 0);
    ClosedClass hit(m, "HitClass", 0, delay, 0);
    ClosedClass miss(m, "MissClass", 0, delay, 0);
    delay.set_service(job, Exp(1.0));

    Routing P;
    P.set(job, job, delay, cache, 1.0);
    P.set(hit, job, cache, delay, 1.0);
    P.set(miss, job, cache, delay, 1.0);
    m.link(P);
    return m;
}

/** An open cache whose hits and misses are served by distinct queues. */
Net gallery_cache_routing() {
    const std::size_t n = 4, mslots = 2;
    Net m("Cache-Routing");
    Source source(m, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{static_cast<int>(mslots)};
    ch.replacestrat = lang::ReplacementStrategy::LRU;
    ch.pread = std::vector<std::vector<double> >{uniform_pread(n), {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = m.add_cache("Cache", ch);
    Queue hq(m, "HitQueue", SchedStrategy::FCFS);
    Queue mq(m, "MissQueue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");

    OpenClass init(m, "InitClass", 0);
    OpenClass hit(m, "HitClass", 0);
    OpenClass miss(m, "MissClass", 0);
    source.set_arrival(init, Exp(1.0));
    hq.set_service(hit, Exp(2.0));
    mq.set_service(miss, Exp(1.0));

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(hit, hit, cache, hq, 1.0);
    P.set(hit, hit, hq, sink, 1.0);
    P.set(miss, miss, cache, mq, 1.0);
    P.set(miss, miss, mq, sink, 1.0);
    m.link(P);
    return m;
}

// ---------------------------------------------------------------------------
// Fork-join and finite capacity
// ---------------------------------------------------------------------------

/** A closed fork-join: five jobs, each split over two PS queues. */
Net gallery_fj_closed() {
    Net m("Fork-Join-Closed");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 5, delay, 0);
    delay.set_service(c1, Exp(1.0));
    q1.set_service(c1, Exp(1.0));
    q2.set_service(c1, Exp(1.0));

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, q1, 1.0);
    P.set(c1, c1, fork, q2, 1.0);
    P.set(c1, c1, q1, join, 1.0);
    P.set(c1, c1, q2, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/**
 * A closed fork-join whose join fires on a 2-of-3 QUORUM: the third sibling is
 * discarded when it arrives. SolverLDES simulates it exactly; SolverMVA and
 * SolverNC charge the second order statistic of the branch completion times
 * (fj_ordstat_exp), which they can only floor at zero here -- see
 * _kb/05-solvers-overview.md.
 */
Net gallery_fj_quorum() {
    Net m("Fork-Join-Quorum");
    Delay delay(m, "Delay");
    Queue q1(m, "Queue1", SchedStrategy::PS);
    Queue q2(m, "Queue2", SchedStrategy::PS);
    Queue q3(m, "Queue3", SchedStrategy::PS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);

    ClosedClass c1(m, "class1", 5, delay, 0);
    delay.set_service(c1, Exp(1.0));
    q1.set_service(c1, Exp(2.0));
    q2.set_service(c1, Exp(2.0));
    q3.set_service(c1, Exp(2.0));
    m.set_join_strategy(join, lang::JoinStrategy::PARTIAL, 2.0);

    Routing P;
    P.set(c1, c1, delay, fork, 1.0);
    P.set(c1, c1, fork, q1, 1.0);
    P.set(c1, c1, fork, q2, 1.0);
    P.set(c1, c1, fork, q3, 1.0);
    P.set(c1, c1, q1, join, 1.0);
    P.set(c1, c1, q2, join, 1.0);
    P.set(c1, c1, q3, join, 1.0);
    P.set(c1, c1, join, delay, 1.0);
    m.link(P);
    return m;
}

/** The open counterpart: a Poisson stream split over two FCFS queues. */
Net gallery_fj_open() {
    Net m("Fork-Join-Open");
    Source source(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t join = m.add_join("Join", fork);
    Sink sink(m, "Sink");

    OpenClass c1(m, "class1", 0);
    source.set_arrival(c1, Exp(0.05));
    q1.set_service(c1, Exp(1.0));
    q2.set_service(c1, Exp(2.0));

    Routing P;
    P.set(c1, c1, source, fork, 1.0);
    P.set(c1, c1, fork, q1, 1.0);
    P.set(c1, c1, fork, q2, 1.0);
    P.set(c1, c1, q1, join, 1.0);
    P.set(c1, c1, q2, join, 1.0);
    P.set(c1, c1, join, sink, 1.0);
    m.link(P);
    return m;
}

/**
 * A finite capacity region of K jobs around one queue, dropping the overflow.
 *
 * The region carries no per-class cap, so the class cap is the reference's
 * unbounded sentinel and the K is the region's GLOBAL maximum, which is what
 * `setGlobalMaxJobs` sets.
 */
Net gallery_fcr(double K) {
    Net m("FCR-Dropping");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1", 0);
    source.set_arrival(c1, Exp(0.8));
    queue.set_service(c1, Exp(1.0));

    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, sink, 1.0);
    m.link(P);

    m.add_region({queue}, std::vector<double>{-1.0}, K,
                 std::vector<lang::DropStrategy>{lang::DropStrategy::DROP});
    return m;
}

// ---------------------------------------------------------------------------
// Random environment
// ---------------------------------------------------------------------------

/**
 * A single server alternating between an UP stage and a DOWN stage.
 *
 * `addNodeFailureRepair` is a macro over two stages and two arcs: the UP stage
 * holds the base network, the DOWN stage the same network with the server's
 * service replaced by the degraded one, and the arcs carry the time to failure
 * and the time to repair.
 */
env::Environment<double> gallery_renv_breakdown() {
    Net model("ServerWithFailures");
    Source source(model, "Arrivals");
    Queue queue(model, "Server", SchedStrategy::FCFS);
    Sink sink(model, "Departures");
    OpenClass jobs(model, "Jobs", 0);
    source.set_arrival(jobs, Exp(0.8));
    queue.set_service(jobs, Exp(2.0));
    queue.set_number_of_servers(1.0);
    Routing P;
    P.set(jobs, jobs, source, queue, 1.0);
    P.set(jobs, jobs, queue, sink, 1.0);
    model.link(P);

    env::Environment<double> e("ServerEnv", 2);
    e.add_node_failure_repair(0, 1, model.get_struct(), "Server", Exp(0.1),
                              Exp(1.0), Exp(0.5));
    e.init();
    return e;
}

}  // namespace examples
}  // namespace line
