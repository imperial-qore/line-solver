/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * One fork-join traversed by an open AND a closed class at once.
 *
 * This is the flat counterpart of test_ln_fork_openarrival.cpp: the same fork
 * carries an exogenous Poisson stream (rate 0.1) and a closed chain (1 job,
 * think time 1). The MVA fork-join transform serves it -- the auxiliary open
 * class it mints for the branches the circulating job did not take coexists with
 * the model's OWN open class, which is exactly what the layered builder cannot
 * do, because there the transform's Source would displace the layer's one.
 *
 * The numbers are the JAR's SolverMVA on the identical model (2026-07-30), which
 * agrees with this port to its printed precision. SolverJMT at seed 23000 / 2e5
 * samples brackets them: Branch1 open 0.3657 vs 0.36184 simulated, Branch1
 * closed 0.35253 vs 0.30955, Join open 0.22961 vs 0.21537. The closed-branch gap
 * is the AMVA fork-join approximation, not a transform defect: the open class,
 * which sees no population constraint, agrees with the simulation to 1%.
 *
 * Native Python returned 0.519566 for the closed class until 2026-07-30, up to
 * 0.9% away: its open/mixed fork-join loop set a `_force_method` attribute that
 * nothing reads, so the inner solve ran exact mixed MVA where mvaDispatch.m
 * :211-215 forces AMVA. All four codebases now agree on the numbers below.
 */

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Station order: Client(delay), PreFork, Branch1, Branch2, PostJoin, Join. */
qn::Network<double> mixed_fj() {
    qn::Network<double> m("ForkJoinOpenClosed");
    const std::size_t src = m.add_source("Source");
    const std::size_t client = m.add_delay("Client");
    const std::size_t prefork = m.add_queue("PreFork", SchedStrategy::FCFS);
    const std::size_t fork = m.add_fork("Fork");
    const std::size_t branch1 = m.add_queue("Branch1", SchedStrategy::FCFS);
    const std::size_t branch2 = m.add_queue("Branch2", SchedStrategy::FCFS);
    const std::size_t join = m.add_join("Join", fork);
    const std::size_t postjoin = m.add_queue("PostJoin", SchedStrategy::FCFS);
    const std::size_t sink = m.add_sink("Sink");

    const std::size_t o = m.add_open_class("Open");
    const std::size_t c = m.add_closed_class("Closed", 1.0, client);

    m.set_arrival(src, o, D::exp_rate(0.1));
    m.set_service(client, c, D::exp_mean(1.0));
    m.set_service(client, o, D::disabled_dist());
    m.set_service(prefork, o, D::exp_mean(0.2));
    m.set_service(prefork, c, D::exp_mean(0.2));
    m.set_service(branch1, o, D::exp_mean(0.3));
    m.set_service(branch1, c, D::exp_mean(0.3));
    m.set_service(branch2, o, D::exp_mean(0.4));
    m.set_service(branch2, c, D::exp_mean(0.4));
    m.set_service(postjoin, o, D::exp_mean(0.1));
    m.set_service(postjoin, c, D::exp_mean(0.1));

    qn::RoutingMatrix<double> P;
    P.set(o, o, src, prefork, 1.0);
    P.set(o, o, prefork, fork, 1.0);
    P.set(o, o, fork, branch1, 1.0);
    P.set(o, o, fork, branch2, 1.0);
    P.set(o, o, branch1, join, 1.0);
    P.set(o, o, branch2, join, 1.0);
    P.set(o, o, join, postjoin, 1.0);
    P.set(o, o, postjoin, sink, 1.0);
    P.set(c, c, client, prefork, 1.0);
    P.set(c, c, prefork, fork, 1.0);
    P.set(c, c, fork, branch1, 1.0);
    P.set(c, c, fork, branch2, 1.0);
    P.set(c, c, branch1, join, 1.0);
    P.set(c, c, branch2, join, 1.0);
    P.set(c, c, join, postjoin, 1.0);
    P.set(c, c, postjoin, client, 1.0);
    m.link(P);
    return m;
}

mva::AvgResult<double> run(qn::Network<double>& m) {
    mva::MvaOptions opt;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

}  // namespace

TEST_CASE("a fork-join carrying an open and a closed class at once is solved") {
    qn::Network<double> m = mixed_fj();
    const mva::AvgResult<double> r = run(m);
    // Source, Client, PreFork, Branch1, Branch2, PostJoin, Join: the Fork is a
    // router and not a station, the transform's auxiliary columns are merged back
    CHECK(r.QN.cols() == 2);

    // station rows: 0 Source, 1 Client, 2 PreFork, 3 Branch1, 4 Branch2,
    // 5 Join, 6 PostJoin -- the Fork is a router and carries no row
    // the open class is carried at its offered rate everywhere on the path
    CHECK(r.TN(2, 0) == doctest::Approx(0.1).epsilon(1e-6));
    CHECK(r.TN(3, 0) == doctest::Approx(0.1).epsilon(1e-6));
    CHECK(r.TN(4, 0) == doctest::Approx(0.1).epsilon(1e-6));
    // and the closed class circulates at the think-time-limited rate
    CHECK(r.TN(1, 1) == doctest::Approx(0.517506).epsilon(1e-5));

    // per-class response times, open then closed, on the two branches
    CHECK(r.RN(3, 0) == doctest::Approx(0.365702).epsilon(1e-5));
    CHECK(r.RN(3, 1) == doctest::Approx(0.352530).epsilon(1e-5));
    CHECK(r.RN(4, 0) == doctest::Approx(0.524448).epsilon(1e-5));
    CHECK(r.RN(4, 1) == doctest::Approx(0.499853).epsilon(1e-5));

    // the Join holds the synchronisation delay for BOTH classes and is never busy
    CHECK(r.RN(5, 0) == doctest::Approx(0.229615).epsilon(1e-5));
    CHECK(r.RN(5, 1) == doctest::Approx(0.219462).epsilon(1e-5));
    CHECK(r.UN(5, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(r.UN(5, 1) == doctest::Approx(0.0).epsilon(1e-9));
    // it sees both branches arrive, so twice the class throughput
    CHECK(r.AN(5, 0) == doctest::Approx(0.2).epsilon(1e-6));
    // and the post-join station sees each job once again
    CHECK(r.RN(6, 0) == doctest::Approx(0.106293).epsilon(1e-5));
    CHECK(r.RN(6, 1) == doctest::Approx(0.101063).epsilon(1e-5));
}
