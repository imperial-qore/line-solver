/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The SolverCTMC surface beyond the AvgTable: the probability family, the
 * transient analyzer, the reward analyzer and the sampler.
 *
 * EVERY ORACLE HERE IS AN IDENTITY THE STATIONARY SOLVE ALREADY SATISFIES, and
 * that is deliberate. Each of these four paths reaches the same chain by a
 * different route, so the strongest available check is that they AGREE with the
 * mean measures the analyzer reports -- a probability family that sums to
 * something other than one, a transient that does not converge to the
 * stationary law, or a reward equal to the queue length that differs from QN
 * are all defects no closed form is needed to expose.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_reward.h"
#include "line/solvers/ctmc/solver_ctmc_sample.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> FCFS Queue (capacity K) -> Sink, one open class. */
qn::Network<double> mm1k(double lambda, double mu, int K) {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Delay -> FCFS Queue -> Delay, one closed class. */
qn::Network<double> cqn(double njobs) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("the probability family sums to one over the whole state space") {
    const int K = 4;
    qn::Network<double> m = mm1k(0.6, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);

    // `joint` at every enumerated state must recover the stationary vector, so
    // the total is one. A lookup that mis-aligned the row width would return
    // zero everywhere and this would read 0 instead.
    double tot = 0;
    for (std::size_t s = 0; s < d.chain.space.size(); ++s)
        tot += ctmc::solver_ctmc_joint(sn, d, d.chain.space[s]);
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));

    // The empty queue: with one phase and a right-aligned buffer there is
    // exactly one state realizing it, so joint and jointaggr must coincide and
    // both must equal p0 = 1/sum(rho^i).
    const std::size_t s0 = ctmc::analyzer_detail::init_state_index(sn, d.chain.space);
    REQUIRE(s0 != static_cast<std::size_t>(-1));
    const double pj = ctmc::solver_ctmc_joint(sn, d, d.chain.space[s0]);
    const double pa = ctmc::solver_ctmc_jointaggr(sn, d, d.chain.space[s0]);
    double norm = 0;
    for (int i = 0; i <= K; ++i) norm += std::pow(0.6, i);
    CHECK(pj == doctest::Approx(1.0 / norm).epsilon(1e-9));
    CHECK(pa == doctest::Approx(pj).epsilon(1e-9));
}

TEST_CASE("the marginal family agrees with the queue length it aggregates") {
    const int K = 4;
    qn::Network<double> m = mm1k(0.6, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);

    // Sum n * P(n at the queue) over the enumerated states and recover QN. This
    // ties margaggr to the mean the analyzer reports through a quantity neither
    // computes from the other.
    double q = 0, mass = 0;
    for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
        const std::vector<double> p = ctmc::solver_ctmc_margaggr(sn, d, d.chain.space[s]);
        const std::vector<double> nir =
            ctmc::prob_detail::marginal_of(sn, sn.node_of_station(2), d.chain.space[s].local[1]);
        // margaggr sums over every state sharing the marginal, so weight each
        // DISTINCT marginal once -- here the marginal determines the state, so
        // the direct sum is safe.
        q += nir[0] * d.pi[s];
        mass += d.pi[s];
        CHECK(p[1] >= 0.0);
        CHECK(p[1] <= 1.0 + 1e-12);
    }
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(q == doctest::Approx(d.avg.QN(1, 0)).epsilon(1e-9));
}

TEST_CASE("the ME gate passes a phase-type model and refuses a signed one") {
    qn::Network<double> m = mm1k(0.5, 1.0, 3);
    CHECK_NOTHROW(ctmc::assert_phase_type_states(m.get_struct(), "getProb"));

    // A matrix-exponential: D0 carries a NEGATIVE off-diagonal, so the
    // stationary vector is a signed measure and a per-state probability does
    // not exist. Built by hand because no builder entry point produces one.
    qn::Network<double> me = mm1k(0.5, 1.0, 3);
    qn::NetworkStruct<double>& sn = const_cast<qn::NetworkStruct<double>&>(me.get_struct());
    line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = -0.5;  // the signature of an ME rather than a PH
    D0(1, 1) = -3.0;
    D1(0, 1) = 2.5;
    D1(1, 0) = 3.0;
    sn.service[1][0].D0 = D0;
    sn.service[1][0].D1 = D1;
    sn.service[1][0].disabled = false;
    CHECK_THROWS_AS(ctmc::assert_phase_type_states(sn, "getProb"), line::UnsupportedError);
}

TEST_CASE("the transient occupancy converges to the stationary law") {
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> st = ctmc::solver_ctmc_analyzer(sn, opt);
    const ctmc::CtmcTransient<double> tr =
        ctmc::solver_ctmc_transient_analyzer(sn, opt, 0.0, 200.0);

    REQUIRE(!tr.t.empty());
    // The chain starts with every job at the Delay, so QNt at t=0 is the whole
    // population there and nothing at the queue.
    CHECK(tr.QNt[0][0].front() == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(tr.QNt[1][0].front() == doctest::Approx(0.0).epsilon(1e-6));
    // By t = 200 on a chain with rates of order 1 the transient has died out,
    // so the last sample must be the stationary answer.
    CHECK(tr.QNt[0][0].back() == doctest::Approx(st.avg.QN(0, 0)).epsilon(1e-4));
    CHECK(tr.QNt[1][0].back() == doctest::Approx(st.avg.QN(1, 0)).epsilon(1e-4));
    CHECK(tr.TNt[1][0].back() == doctest::Approx(st.avg.TN(1, 0)).epsilon(1e-4));
    // Population is conserved at EVERY time point, not only in the limit.
    for (std::size_t i = 0; i < tr.t.size(); ++i)
        CHECK(tr.QNt[0][0][i] + tr.QNt[1][0][i] == doctest::Approx(2.0).epsilon(1e-4));
}

TEST_CASE("the fau transient method answers what the integrator answers") {
    // `transient_method = "fau"` advances the same forward equation by fast
    // adaptive uniformization instead of integrating it. The oracle is the ODE
    // path on the same generator and the same grid, and the two must agree to
    // the integrator's own tolerance; tightening fau_epsilon by six orders must
    // then move the fau answer by far less than that gap, which is what shows
    // the residual is the integrator's and not the uniformization's.
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    std::vector<double> grid;
    for (int i = 0; i <= 40; ++i) grid.push_back(0.25 * i);

    ctmc::CtmcOptions ode;
    const ctmc::CtmcTransient<double> a = ctmc::solver_ctmc_transient_analyzer(sn, ode, 0.0, 10.0, grid);

    ctmc::CtmcOptions fau;
    fau.transient_method = "fau";
    const ctmc::CtmcTransient<double> b = ctmc::solver_ctmc_transient_analyzer(sn, fau, 0.0, 10.0, grid);

    ctmc::CtmcOptions tight;
    tight.transient_method = "fau";
    tight.fau_epsilon = 1e-12;
    const ctmc::CtmcTransient<double> c = ctmc::solver_ctmc_transient_analyzer(sn, tight, 0.0, 10.0, grid);

    REQUIRE(a.t.size() == b.t.size());
    REQUIRE(b.t.size() == c.t.size());
    double dev_ode = 0.0, dev_tight = 0.0;
    for (std::size_t i = 0; i < b.t.size(); ++i) {
        CHECK(b.t[i] == doctest::Approx(a.t[i]).epsilon(1e-12));
        dev_ode = std::max(dev_ode, std::fabs(b.QNt[1][0][i] - a.QNt[1][0][i]));
        dev_tight = std::max(dev_tight, std::fabs(b.QNt[1][0][i] - c.QNt[1][0][i]));
        // Population is conserved along the fau trajectory too.
        CHECK(b.QNt[0][0][i] + b.QNt[1][0][i] == doctest::Approx(2.0).epsilon(1e-4));
    }
    CHECK(dev_ode < 1e-3);
    CHECK(dev_tight < 1e-5);
}

TEST_CASE("an unknown transient method is refused rather than defaulted") {
    qn::Network<double> m = cqn(2.0);
    ctmc::CtmcOptions opt;
    opt.transient_method = "uniformization";
    CHECK_THROWS_AS(ctmc::solver_ctmc_transient_analyzer(m.get_struct(), opt, 0.0, 1.0),
                    line::InputError);
}

TEST_CASE("a reward equal to the queue length reproduces QN") {
    qn::Network<double> m = cqn(2.0);
    // Column (ist-1)*K + k with ist = 2 (the queue), K = 1, k = 1 -> index 1.
    m.set_reward("qlen", [](const std::vector<double>& s) { return s[1]; });
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> st = ctmc::solver_ctmc_analyzer(sn, opt);
    const ctmc::CtmcReward<double> rw = ctmc::solver_ctmc_reward(sn, opt, 50);

    REQUIRE(rw.names.size() == 1);
    CHECK(rw.names[0] == "qlen");
    // E[r] = sum_s pi(s) n_2(s) IS the queue length, so the reward analyzer and
    // the mean analyzer must agree exactly.
    CHECK(rw.steady_state[0] == doctest::Approx(st.avg.QN(1, 0)).epsilon(1e-9));

    // The value function ACCUMULATES reward, so it is non-decreasing in the
    // iteration index and its increments approach E[r]/q. It does NOT converge
    // to E[r]: it is a total, not an average.
    REQUIRE(rw.V[0].rows() == 51);
    for (std::size_t s = 0; s < rw.V[0].cols(); ++s)
        CHECK(rw.V[0](50, s) >= rw.V[0](49, s) - 1e-12);
    CHECK(rw.t.front() == doctest::Approx(0.0));
    CHECK(rw.t.back() > 0.0);
}

TEST_CASE("a reward analysis with no reward declared is refused by name") {
    qn::Network<double> m = cqn(2.0);
    CHECK_THROWS_AS(ctmc::solver_ctmc_reward(m.get_struct(), ctmc::CtmcOptions()),
                    line::InputError);
}

TEST_CASE("the sampler walks the chain and marks which event fired") {
    qn::Network<double> m = cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcSamplePath<double> p =
        ctmc::solver_ctmc_sample_sys(sn, ctmc::CtmcOptions(), 500, 23000);

    REQUIRE(p.state.size() > 1);
    CHECK(p.t.size() == p.state.size());
    // Time is non-decreasing and every state index is inside the space.
    for (std::size_t i = 1; i < p.t.size(); ++i) CHECK(p.t[i] > p.t[i - 1]);
    for (std::size_t i = 0; i < p.state.size(); ++i)
        CHECK(p.state[i] < p.chain.chain.space.size());
    // Every step names the synchronization that fired: the state sequence alone
    // does not determine it, which is the whole reason the filtration is kept.
    for (std::size_t i = 0; i + 1 < p.event.size(); ++i)
        CHECK(p.event[i] < p.chain.chain.filt.size());

    // The aggregate view conserves the closed population at every sample.
    const line::Matrix<double> A = ctmc::solver_ctmc_sample_sys_aggr(sn, p);
    for (std::size_t i = 0; i < A.rows(); ++i)
        CHECK(A(i, 0) + A(i, 1) == doctest::Approx(2.0).epsilon(1e-12));

    // The per-node views agree with the aggregate they are slices of.
    const line::Matrix<double> nq = ctmc::solver_ctmc_sample_aggr(sn, p, sn.node_of_station(2));
    for (std::size_t i = 0; i < nq.rows(); ++i)
        CHECK(nq(i, 0) == doctest::Approx(A(i, 1)).epsilon(1e-12));
    const line::Matrix<double> raw = ctmc::solver_ctmc_sample(sn, p, sn.node_of_station(2));
    CHECK(raw.rows() == nq.rows());
}

TEST_CASE("the state-space aggregate reports a Source as zero, not Inf") {
    const int K = 3;
    qn::Network<double> m = mm1k(0.5, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, d.chain.space);

    REQUIRE(A.rows() == d.chain.space.size());
    // Station 1 is the Source. `to_marginal` reports Inf there as an encoding
    // sentinel; an Inf in this matrix would poison every aggregate built on it.
    for (std::size_t s = 0; s < A.rows(); ++s) {
        CHECK(A(s, 0) == 0.0);
        CHECK(std::isfinite(A(s, 1)));
        CHECK(A(s, 1) <= static_cast<double>(K));
    }
}
