/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The SolverSSA queries that are not the average table: the four probabilities
 * and the four samplers of `solver_ssa_getters.h`.
 *
 * THE ORACLE IS THE EXACT CHAIN, ON THE SAME STATE. Both solvers answer
 * `getProb` at the model's default initial state, so the pair measures the
 * simulation error directly rather than inferring it -- which is the reason the
 * C++ CLI reports the SSA probabilities at that state and not at one of the
 * simulator's choosing.
 *
 * TWO IDENTITIES HOLD EXACTLY AND ARE ASSERTED FIRST. A joint probability can
 * never exceed the aggregate one, since the aggregate sums over every
 * arrangement realizing the same per-class counts; and the trajectory's time
 * column is strictly increasing, because every holding time is `-log(u)/R` with
 * `u` strictly inside (0,1). Neither can be perturbed by Monte Carlo error.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ssa/solver_ssa_getters.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace ssa = line::ssa;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Think(mean 2) -> FCFS Queue(rate 3), `n` closed jobs: a four-state chain. */
qn::Network<double> cqn(double n) {
    qn::Network<double> m("ssa-getters-cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("SolverSSA getProb: the time fraction at the initial state, against the chain") {
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;  // stated with the seed: a simulated number needs both
    opt.seed = 23000;
    const ssa::SsaProbReport r = ssa::solver_ssa_prob<double>(sn, opt);

    REQUIRE(r.samples == 200000);
    REQUIRE(r.seed == 23000);
    REQUIRE(r.marg.size() == sn.nstations);

    // EXACT: the joint probability of one arrangement cannot exceed the
    // aggregate over every arrangement with the same counts.
    CHECK(r.sys.prob <= r.sys_aggr.prob + 1e-12);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        CHECK(r.marg[i].prob <= r.aggr[i].prob + 1e-12);
    // The initial state is where the path starts, so it cannot be unseen.
    CHECK(r.sys.seen);

    // Against the exact stationary law at the SAME state. The band is wide
    // because a probability estimated from one path carries the same Monte
    // Carlo error a mean does: 10 per cent here is roughly ten standard errors
    // at 2e5 firings on a four-state chain, so a violation is a bias.
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, ctmc::CtmcOptions());
    qn::NetState<double> init;
    REQUIRE(ctmc::analyzer_detail::default_init_state(sn, init));
    const double psys = ctmc::solver_ctmc_joint(sn, d, init);
    const double psysaggr = ctmc::solver_ctmc_jointaggr(sn, d, init);
    CHECK(r.sys.prob == doctest::Approx(psys).epsilon(0.10));
    CHECK(r.sys_aggr.prob == doctest::Approx(psysaggr).epsilon(0.10));
}

TEST_CASE("SolverSSA sampleSys: one trajectory, its aggregate and one node's block") {
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ssa::SsaSerialOptions opt;
    opt.samples = 5000;
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);
    const ssa::SsaSamplePath<double> sys = ssa::ssa_sample_sys(sn, sim.run);

    REQUIRE(sys.t.size() == 5000);
    REQUIRE(sys.state.rows() == 5000);
    REQUIRE(sys.aggr.rows() == 5000);
    CHECK(sys.seed == 23000);

    // EXACT: time is strictly increasing, and every sampled state holds the
    // whole closed population. The second is the sharpest available check that
    // the trajectory indexes the states it claims to.
    bool increasing = true, conserved = true;
    for (std::size_t i = 0; i < sys.t.size(); ++i) {
        if (i > 0 && !(sys.t[i] > sys.t[i - 1])) increasing = false;
        double tot = 0;
        for (std::size_t c = 0; c < sys.aggr.cols(); ++c) tot += sys.aggr(i, c);
        if (std::fabs(tot - 2.0) > 1e-9) conserved = false;
    }
    CHECK(increasing);
    CHECK(conserved);

    // The per-node view: the same trajectory, one node's row and its counts.
    const std::size_t qnode = sn.node_of_station(2);
    const ssa::SsaSamplePath<double> nd = ssa::ssa_sample_node(sn, sim.run, qnode);
    REQUIRE(nd.aggr.rows() == sys.aggr.rows());
    bool matches = true;
    for (std::size_t i = 0; i < nd.aggr.rows(); ++i)
        if (std::fabs(nd.aggr(i, 0) - sys.aggr(i, (2 - 1) * sn.nclasses)) > 1e-9) matches = false;
    CHECK(matches);

    // A stateless node has no state to sample, and says so.
    CHECK_THROWS_AS(ssa::ssa_sample_node(sn, sim.run, sn.nodes.size() + 1), line::Error);
}

TEST_CASE("SolverSSA getCdfRespT is refused, and the refusal is the reference's answer") {
    // `@@SolverSSA/getCdfRespT.m` raises in MATLAB too: SSA samples state
    // trajectories, not per-job sojourn times, so the inherited exponential fit
    // would look measured while carrying no information about the tail.
    CHECK_THROWS_AS(ssa::ssa_cdf_respt_refuse(), line::UnsupportedError);
}
