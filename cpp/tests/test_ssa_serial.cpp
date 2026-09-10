/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The serial SolverSSA engine: `solver_ssa_reachability.m`, the run loop of
 * `solver_ssa.m` and `solver_ssa_analyzer_serial.m`.
 *
 * WHAT A SIMULATION TEST CAN AND CANNOT ASSERT. A sample path is a random
 * object, so nothing here compares a simulated number against an exact one at a
 * tight tolerance. Three kinds of oracle are used instead, in increasing order
 * of strength:
 *
 *   EXACT AND DETERMINISTIC. A seeded run is a pure function of its seed, so two
 *   runs at the same seed must produce the IDENTICAL trace, bit for bit. This
 *   catches any hidden dependence on iteration order or uninitialized state.
 *
 *   EXACT AND STRUCTURAL. The states the path visits must lie in the reachable
 *   space the CTMC walk enumerates, and in a closed network the queue lengths
 *   must sum to the population for EVERY sample, hence for the time average --
 *   an identity no amount of Monte Carlo error can perturb, which is what makes
 *   it the sharpest bias detector available here.
 *
 *   STATISTICAL, WITH A STATED BAND. A time-averaged mean is compared against
 *   the exact CTMC answer inside a band chosen an order of magnitude above the
 *   run's standard error, so that a failure means a BIAS and never noise. The
 *   band is stated with its reasoning at the assertion.
 *
 * Every simulated result below is reported with its sample count and its seed,
 * because a simulated number without them is not a measurement.
 */
#include <cmath>
#include <map>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state_events.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ssa/solver_ssa_serial.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace ssa = line::ssa;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Think -> FCFS Queue -> Think, one closed class of `n` jobs. */
qn::Network<double> cqn(double n, double think_mean, double mu) {
    qn::Network<double> m("ssa-serial-cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(1.0 / think_mean));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The flattened key of a network state, as the state tables index it. */
std::vector<double> key_of(const qn::NetState<double>& s) {
    return ctmc::ctmc_detail::state_key(s);
}

}  // namespace

TEST_CASE("SolverSSA serial: a seeded run reproduces its trace exactly") {
    qn::Network<double> m = cqn(2.0, 2.0, 3.0);
    ssa::SsaSerialOptions opt;
    opt.samples = 2000;
    opt.seed = 12345;

    ssa::SsaSerialEngine<double> e1(m.get_struct(), opt);
    ssa::SsaSerialEngine<double> e2(m.get_struct(), opt);
    const ssa::SsaSerialRun<double> r1 = e1.run();
    const ssa::SsaSerialRun<double> r2 = e2.run();

    REQUIRE(r1.samples == 2000);
    REQUIRE(r1.seed == 12345);
    REQUIRE(r2.tran_sync.size() == r1.tran_sync.size());
    // Bit-for-bit, not to a tolerance: the same seed drives the same draws
    // through the same arithmetic, so any difference at all is a defect and not
    // a rounding artefact.
    CHECK(r1.tran_sync == r2.tran_sync);
    bool same_times = true;
    for (std::size_t i = 0; i < r1.tran_time.size(); ++i)
        same_times = same_times && r1.tran_time[i] == r2.tran_time[i];
    CHECK(same_times);
    CHECK(r1.simulated_time == r2.simulated_time);

    // A DIFFERENT seed is a different stream. With 2000 firings over a chain
    // with several enabled transitions per state, two streams agreeing on the
    // whole trace has probability below 2^-2000; a failure here means the seed
    // is not reaching the generator.
    ssa::SsaSerialOptions other = opt;
    other.seed = 999;
    ssa::SsaSerialEngine<double> e3(m.get_struct(), other);
    const ssa::SsaSerialRun<double> r3 = e3.run();
    CHECK(r3.tran_sync != r1.tran_sync);
    CHECK(r3.seed == 999);
}

TEST_CASE("SolverSSA serial: reachability is the CTMC walk, decomposed per node") {
    qn::Network<double> m = cqn(2.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const ssa::SsaReachability<double> rr = ssa::solver_ssa_reachability(sn);

    // The same walk, reached independently through the CTMC entry point and
    // seeded with the same widened initial state.
    const std::vector<std::size_t> cut(sn.nclasses, 0);  // closed model: no cutoff
    const qn::NetState<double> init = ssa::serial_detail::wide_init_state(sn, cut);
    const std::vector<qn::Sync<double>> sync = qn::refresh_sync(sn);
    const std::vector<qn::GlobalSync<double>> gsync = qn::refresh_global_sync(sn);
    const std::vector<qn::NetState<double>> walk =
        ctmc::reachable_space_generator(sn, init, sync, gsync);

    REQUIRE(!rr.space.empty());
    REQUIRE(rr.space.size() == walk.size());
    std::map<std::vector<double>, int> seen;
    for (std::size_t s = 0; s < walk.size(); ++s) seen[key_of(walk[s])] += 1;
    for (std::size_t s = 0; s < rr.space.size(); ++s) seen[key_of(rr.space[s])] += 2;
    for (std::map<std::vector<double>, int>::const_iterator it = seen.begin(); it != seen.end();
         ++it)
        CHECK(it->second == 3);  // in both, exactly once each

    // With two jobs and a single-server FCFS queue every encoded state is also
    // reachable, so the walk must equal the full lattice enumeration too. This
    // is what makes the check above more than a tautology on the shared walk.
    const std::vector<qn::NetState<double>> full = qn::space_generator(sn, cut);
    CHECK(full.size() == rr.space.size());

    // `SSh` and `space{i}` must reconstruct the state they were split from --
    // the decomposition is the only part of the reference's reachability this
    // port adds on top of the shared walk.
    REQUIRE(rr.hash.size() == rr.space.size());
    for (std::size_t s = 0; s < rr.space.size(); ++s) {
        REQUIRE(rr.hash[s].size() == sn.stateful_nodes.size());
        std::size_t col = 0;
        for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f) {
            const std::size_t h = rr.hash[s][f];
            REQUIRE(h >= 1);
            REQUIRE(h <= rr.node_space[f].size());
            CHECK(rr.node_space[f][h - 1] == rr.space[s].local[f]);
            for (std::size_t j = 0; j < rr.space[s].local[f].size(); ++j) {
                CHECK(rr.ssq(s, col) == rr.space[s].local[f][j]);
                ++col;
            }
        }
        CHECK(col == rr.ssq.cols());
    }
}

TEST_CASE("SolverSSA serial: every state the path visits is a reachable state") {
    qn::Network<double> m = cqn(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ssa::SsaSerialOptions opt;
    opt.samples = 5000;
    opt.seed = 4242;
    ssa::SsaSerialEngine<double> eng(sn, opt);
    const ssa::SsaSerialRun<double> run = eng.run();

    const ssa::SsaReachability<double> rr = ssa::solver_ssa_reachability(sn, opt);
    std::map<std::vector<double>, bool> reachable;
    for (std::size_t s = 0; s < rr.space.size(); ++s) reachable[key_of(rr.space[s])] = true;

    REQUIRE(!run.space.empty());
    for (std::size_t s = 0; s < run.space.size(); ++s)
        CHECK(reachable.find(key_of(run.space[s])) != reachable.end());
    // The path is ergodic enough at 5000 firings to have covered the whole
    // 4-state chain; a path that visits fewer states than it can reach is a
    // stuck sample path, which is the failure mode this catches.
    CHECK(run.space.size() == rr.space.size());

    // The time weights are a probability distribution over the visited states.
    double tot = 0;
    for (std::size_t s = 0; s < run.pi.size(); ++s) {
        CHECK(run.pi[s] >= 0.0);
        tot += run.pi[s];
    }
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("SolverSSA serial: means lie in a Monte Carlo band around the exact CTMC") {
    // Think(mean 2) -> FCFS Queue(rate 3), three jobs: the same closed network
    // SolverCTMC and SolverMVA agree on exactly, so the CTMC answer below is an
    // EXACT oracle and every deviation is the simulation's.
    qn::Network<double> m = cqn(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const line::mva::AvgResult<double> exact =
        ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    ssa::SsaSerialOptions opt;
    opt.samples = 100000;  // reported with the seed below: neither is optional
    opt.seed = 23000;      // LINE's own default seed
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

    REQUIRE(sim.avg.method == "serial");
    REQUIRE(sim.avg.samples == 100000);
    REQUIRE(sim.seed == 23000);
    REQUIRE(sim.avg.simulated_time > 0.0);

    // EXACT, not statistical: every state of this chain holds exactly three
    // jobs, so the time average of the total queue length is three whatever the
    // path did. A bias in the state decode or in the time weighting breaks this
    // identity while leaving the individual means plausible.
    CHECK(sim.avg.QN(0, 0) + sim.avg.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-9));

    // THE BAND. The chain has four states and mixes within a few firings, so at
    // 1e5 firings the effective sample size is of order 1e4 and the relative
    // standard error of a time-averaged mean is around one per cent. The 10 per
    // cent band below is roughly ten standard errors: a violation is a BIAS in
    // the estimator, not an unlucky path. It is deliberately NOT tightened to
    // the observed error, because a band at the noise level turns an honest
    // simulation into a flaky test.
    const double band = 0.10;
    CHECK(sim.avg.QN(0, 0) == doctest::Approx(exact.QN(0, 0)).epsilon(band));
    CHECK(sim.avg.QN(1, 0) == doctest::Approx(exact.QN(1, 0)).epsilon(band));
    CHECK(sim.avg.TN(0, 0) == doctest::Approx(exact.TN(0, 0)).epsilon(band));
    CHECK(sim.avg.TN(1, 0) == doctest::Approx(exact.TN(1, 0)).epsilon(band));
    CHECK(sim.avg.UN(1, 0) == doctest::Approx(exact.UN(1, 0)).epsilon(band));
    CHECK(sim.avg.XN[0] == doctest::Approx(exact.XN[0]).epsilon(band));
    // Little's law holds on the SIMULATED numbers by construction, which is a
    // check on the reduction rather than on the path.
    CHECK(sim.avg.RN(1, 0) * sim.avg.TN(1, 0) ==
          doctest::Approx(sim.avg.QN(1, 0)).epsilon(1e-12));
    // Response time over one cycle: N = X * (R_think + R_queue), the closed
    // network's own Little's law, which the serial reduction must satisfy to
    // simulation accuracy.
    CHECK(sim.avg.XN[0] * (sim.avg.RN(0, 0) + sim.avg.RN(1, 0)) ==
          doctest::Approx(3.0).epsilon(band));
}

TEST_CASE("SolverSSA serial: an open model truncates where SolverCTMC truncates") {
    // Source -> FCFS Queue (capacity 4) -> Sink. The buffer is a PHYSICAL
    // capacity, so both solvers see the same finite chain and the simulated
    // means are comparable to the exact ones inside the same band as above.
    const double lambda = 0.6, mu = 1.0;
    qn::Network<double> m("ssa-serial-mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, 4);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, sk, 1.0);
    m.link(P);

    ctmc::CtmcOptions copt;
    copt.cutoff = 4;
    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(m.get_struct(), copt);

    ssa::SsaSerialOptions opt;
    opt.samples = 100000;
    opt.seed = 7;
    opt.cutoff = 4;
    const ssa::SsaSerialSolution<double> sim =
        ssa::solver_ssa_serial_analyzer(m.get_struct(), opt);

    // M/M/1/4 at rho = 0.6: the exact mean queue length is about 1.24 and the
    // carried throughput about 0.566. The same 10 per cent band applies, and for
    // the same reason: five states, fast mixing, 1e5 firings at seed 7.
    const double band = 0.10;
    CHECK(sim.avg.QN(1, 0) == doctest::Approx(exact.QN(1, 0)).epsilon(band));
    CHECK(sim.avg.TN(1, 0) == doctest::Approx(exact.TN(1, 0)).epsilon(band));
    CHECK(sim.avg.UN(1, 0) == doctest::Approx(exact.UN(1, 0)).epsilon(band));
    // A Source holds no jobs, so its queue length is zero by definition and not
    // by measurement -- the same rule the NRM engine and the CTMC apply.
    CHECK(sim.avg.QN(0, 0) == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("SolverSSA serial: the warmup discard drops only the time average") {
    qn::Network<double> m = cqn(3.0, 2.0, 3.0);
    ssa::SsaSerialOptions opt;
    opt.samples = 4000;
    opt.seed = 31337;
    opt.warmupfrac = 0.25;
    ssa::SsaSerialEngine<double> eng(m.get_struct(), opt);
    const ssa::SsaSerialRun<double> run = eng.run();

    CHECK(run.warmup == 1000);
    CHECK(run.samples == 4000);
    // The trace keeps every firing; only the weights lose the transient, so the
    // simulated time the metrics average over is shorter than the elapsed time.
    CHECK(run.tran_time.size() == 4000);
    CHECK(run.simulated_time < run.tran_time.back());
    double tot = 0;
    for (std::size_t s = 0; s < run.pi.size(); ++s) tot += run.pi[s];
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));
}

/** Delay -> Fork -> {Q1, Q2} -> Join -> Delay, one closed class of `n` jobs. */
qn::Network<double> fjmodel(double n) {
    qn::Network<double> fj("ssa-serial-fj");
    const std::size_t d = fj.add_delay("Delay");
    const std::size_t f = fj.add_fork("Fork");
    const std::size_t q1 = fj.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = fj.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = fj.add_join("Join", f);
    const std::size_t c = fj.add_closed_class("C1", n, d);
    fj.set_service(d, c, Dist::exp_rate(1.0));
    fj.set_service(q1, c, Dist::exp_rate(2.0));
    fj.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    fj.link(P);
    return fj;
}

TEST_CASE("SolverSSA serial: refusals name the construct that is missing") {
    // A fork-join model reaches the engine only through the TAG-AUGMENTED copy,
    // which is what carries the `fjsync` firing list; the analyzer builds it, so
    // the refusal below is what a caller driving the engine by hand reads.
    qn::Network<double> fj("ssa-serial-fj");
    const std::size_t d = fj.add_delay("Delay");
    const std::size_t f = fj.add_fork("Fork");
    const std::size_t q1 = fj.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = fj.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = fj.add_join("Join", f);
    const std::size_t c = fj.add_closed_class("C1", 2.0, d);
    fj.set_service(d, c, Dist::exp_rate(1.0));
    fj.set_service(q1, c, Dist::exp_rate(2.0));
    fj.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    fj.link(P);
    ssa::SsaSerialOptions fopt;
    fopt.samples = 100;
    // The analyzer augments and runs; the engine and the reachability walk, both
    // of which take the struct they are given, refuse the un-augmented one by
    // name rather than simulating a fork that never fires.
    CHECK_NOTHROW(ssa::solver_ssa_serial_analyzer(fj.get_struct(), fopt));
    CHECK_THROWS_AS(ssa::SsaSerialEngine<double>(fj.get_struct(), fopt), line::UnsupportedError);
    CHECK_THROWS_AS(ssa::solver_ssa_reachability(fj.get_struct(), fopt), line::UnsupportedError);

    // The parallel method replicates the engine and averages; one replica has a
    // different variance from the average of many, so it refuses rather than
    // answering with a number at the wrong precision.
    qn::Network<double> m = cqn(2.0, 2.0, 3.0);
    ssa::SsaSerialOptions popt;
    popt.samples = 100;
    popt.method = "parallel";
    CHECK_THROWS_AS(ssa::solver_ssa_serial(m.get_struct(), popt), line::UnsupportedError);
    popt.method = "amva";
    CHECK_THROWS_AS(ssa::solver_ssa_serial(m.get_struct(), popt), line::UnsupportedError);
    // 'serial', 'default' and 'ssa' all reach the engine.
    popt.method = "serial";
    CHECK_NOTHROW(ssa::solver_ssa_serial(m.get_struct(), popt));
}

TEST_CASE("SolverSSA serial: exact arithmetic is refused, not narrowed") {
    // An exponential clock is -log(u)/rate. There is no exact value for a
    // rational backend to compute and no information a wider float would carry,
    // since the error is the Monte Carlo error and not the rounding.
    qn::Network<line::Rational> m("ssa-serial-exact");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, line::lang::Distrib<line::Rational>::exp_rate(
                            line::num_traits<line::Rational>::from_double(0.5)));
    m.set_service(q, c, line::lang::Distrib<line::Rational>::exp_rate(
                            line::num_traits<line::Rational>::from_int(3)));
    qn::RoutingMatrix<line::Rational> P;
    P.set(d, q, line::num_traits<line::Rational>::from_int(1));
    P.set(q, d, line::num_traits<line::Rational>::from_int(1));
    m.link(P);

    ssa::SsaSerialOptions opt;
    opt.samples = 10;
    CHECK_THROWS_AS(ssa::solver_ssa_serial_analyzer(m.get_struct(), opt), line::UnsupportedError);
}

namespace {

/**
 * Source -> POLLING queue -> Sink, two open classes, declared the way
 * `Queue.setPollingType` / `Queue.setSwitchover` declare it: the MATLAB-faithful
 * pair, which is also what the JSON reader emits. `set_polling` writes the same
 * controller into a different field, and the two must resolve identically.
 */
qn::Network<double> polling_open(bool via_set_polling, bool timed_switchover,
                                 line::lang::PollingType ty = line::lang::PollingType::EXHAUSTIVE) {
    const double l1 = 0.2, l2 = 0.3, mu = 2.0;
    qn::Network<double> m("ssa-serial-poll");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("P", SchedStrategy::POLLING);
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dist::exp_rate(l1));
    m.set_arrival(src, c2, Dist::exp_rate(l2));
    m.set_service(q, c1, Dist::exp_rate(mu));
    m.set_service(q, c2, Dist::exp_rate(mu));
    std::vector<Dist> sw(2);
    sw[0] = timed_switchover ? Dist::exp_rate(4.0) : Dist::immediate();
    sw[1] = sw[0];
    if (via_set_polling) {
        m.set_polling(q, ty, sw, 2);
    } else {
        m.set_polling_type(q, ty, ty == line::lang::PollingType::KLIMITED ? 2 : 0);
        m.set_switchover(q, c1, sw[0]);
        m.set_switchover(q, c2, sw[1]);
    }
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, sk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, sk, 1.0);
    m.link(P);
    m.set_capacity(q, 4);
    return m;
}

}  // namespace

TEST_CASE("SolverSSA serial: a polling controller exists however it was declared") {
    // THE SEGMENTATION FAULT THIS PINS. `polling_info` used to build a
    // controller only when `sn.pollingparam` carried one, i.e. only for a model
    // built through `set_polling`. A station declared through setPollingType /
    // setSwitchover got `valid = false` and, with it, an EMPTY `polled` -- and
    // `polling_next` indexes `polled[pos-1]` with no size test. An empty
    // std::vector<bool> holds a null word pointer, so that index is a SIGSEGV
    // and not a garbage read, which is why the failure was a crash inside
    // after_event_station_arv rather than a wrong number.
    //
    // The reference settles it: State.pollingInfo returns [] only for a node
    // that is not a polling station, and otherwise defaults to EXHAUSTIVE with
    // every switchover immediate. Every polling station HAS a controller.
    for (int api = 0; api < 2; ++api) {
        qn::Network<double> m = polling_open(api == 1, false);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
        CHECK(pi.valid);
        REQUIRE(pi.polled.size() == sn.nclasses);
        CHECK(static_cast<bool>(pi.polled[0]));
        CHECK(static_cast<bool>(pi.polled[1]));
        CHECK((pi.ptype == line::lang::PollingType::EXHAUSTIVE));
        // Immediate switchovers are folded, so the controller is zero-width --
        // and the width nvars reserves must agree with it, or every read of the
        // local-variable block is shifted by the difference.
        CHECK(pi.width == 0);
        CHECK(sn.nvars_of(2) == 0);
    }

    // A timed switchover materializes pos and swk, again through either API.
    for (int api = 0; api < 2; ++api) {
        qn::Network<double> m = polling_open(api == 1, true);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
        REQUIRE(pi.valid);
        CHECK(pi.width == 2);
        CHECK(sn.nvars_of(2) == 2);
    }
}

TEST_CASE("SolverSSA serial: a polling station simulates against the exact CTMC") {
    // The model that crashed: the serial engine hands every state to the same
    // handlers the CTMC generator uses, so the polling controller is exercised
    // on the simulated path rather than only on the enumerated space.
    for (int sw = 0; sw < 2; ++sw) {
        qn::Network<double> m = polling_open(false, sw == 1);
        const qn::NetworkStruct<double>& sn = m.get_struct();

        ctmc::CtmcOptions copt;
        copt.cutoff = 4;
        const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, copt);

        ssa::SsaSerialOptions opt;
        opt.samples = 60000;
        opt.seed = 4242;
        opt.cutoff = 4;
        const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

        // 6e4 firings at seed 4242 on a chain of a few dozen states. The band is
        // 12 per cent, an order of magnitude above the run's standard error, so
        // a failure is a bias in the controller dynamics and not noise.
        const double band = 0.12;
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(std::isfinite(sim.avg.QN(1, r)));
            CHECK(sim.avg.QN(1, r) == doctest::Approx(exact.QN(1, r)).epsilon(band));
            CHECK(sim.avg.TN(1, r) == doctest::Approx(exact.TN(1, r)).epsilon(band));
            CHECK(sim.avg.UN(1, r) == doctest::Approx(exact.UN(1, r)).epsilon(band));
        }
        // Physics, independent of the polling order and of the reference table:
        // at rho = 0.25 with capacity 4 the loss is a fraction of a per cent, so
        // each class carries essentially its whole arrival rate.
        CHECK(sim.avg.TN(1, 0) == doctest::Approx(0.2).epsilon(band));
        CHECK(sim.avg.TN(1, 1) == doctest::Approx(0.3).epsilon(band));
    }
}

TEST_CASE("SolverSSA serial: the closed single-buffer polling model the NRM refuses") {
    // The shape the dispatcher would hand the serial engine on a fallback, and
    // the one that faulted: one closed class, so the cyclic order has a single
    // buffer and the server pays a timed switchover before every visit. The
    // controller is declared through setPollingType / setSwitchover, which is
    // what made polling_info build nothing at all.
    qn::Network<double> m("ssa-serial-poll-closed");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("PollQ", SchedStrategy::POLLING);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    m.set_polling_type(q, line::lang::PollingType::EXHAUSTIVE);
    m.set_switchover(q, c, Dist::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // A single non-immediate switchover materializes pos and swk.
    const qn::PollingInfo<double> pi = qn::polling_info(sn, 2);
    REQUIRE(pi.valid);
    REQUIRE(pi.polled.size() == 1);
    CHECK(pi.width == 2);
    CHECK(sn.nvars_of(2) == 2);

    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    ssa::SsaSerialOptions opt;
    opt.samples = 40000;
    opt.seed = 20260728;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

    // The population identity holds for EVERY sample, so it holds exactly for
    // the time average: no Monte Carlo error can perturb it, which makes it the
    // sharpest available check that the controller conserves jobs.
    CHECK(sim.avg.QN(0, 0) + sim.avg.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
    const double band = 0.10;
    CHECK(sim.avg.QN(1, 0) == doctest::Approx(exact.QN(1, 0)).epsilon(band));
    CHECK(sim.avg.TN(1, 0) == doctest::Approx(exact.TN(1, 0)).epsilon(band));
    CHECK(sim.avg.UN(1, 0) == doctest::Approx(exact.UN(1, 0)).epsilon(band));
}

TEST_CASE("SolverSSA serial: every polling discipline reaches the engine") {
    // The controller columns differ by discipline -- EXHAUSTIVE carries none,
    // the other three carry the visit budget ctr -- so each is a distinct path
    // through polling_get / polling_next / polling_land. The oracle is the exact
    // CTMC on the same struct, which drives the same handlers over the
    // enumerated space rather than over a sample path.
    const line::lang::PollingType tys[] = {
        line::lang::PollingType::EXHAUSTIVE, line::lang::PollingType::GATED,
        line::lang::PollingType::KLIMITED, line::lang::PollingType::DECREMENTING};
    for (int t = 0; t < 4; ++t) {
        qn::Network<double> m = polling_open(false, false, tys[t]);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        const std::size_t want = tys[t] == line::lang::PollingType::EXHAUSTIVE ? 0 : 1;
        CHECK(sn.nvars_of(2) == want);

        ctmc::CtmcOptions copt;
        copt.cutoff = 4;
        const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, copt);

        ssa::SsaSerialOptions opt;
        opt.samples = 20000;
        opt.seed = 8080;
        opt.cutoff = 4;
        const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);
        // 2e4 firings: a wider band than the case above, for the same reason a
        // shorter run has a larger standard error.
        const double band = 0.15;
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(std::isfinite(sim.avg.QN(1, r)));
            CHECK(sim.avg.QN(1, r) == doctest::Approx(exact.QN(1, r)).epsilon(band));
            CHECK(sim.avg.TN(1, r) == doctest::Approx(exact.TN(1, r)).epsilon(band));
        }
    }
}

namespace {

/** Think -> preemptive Queue -> Think, one closed class, Erlang-2 service. */
qn::Network<double> preempt_cqn(SchedStrategy s, double njobs) {
    qn::Network<double> m("ssa-serial-preempt");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", s);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::erlang(4.0, 2));  // mean 0.5, two phases
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

/** The same, two closed classes at DISTINCT priorities, for the PRIO variants. */
qn::Network<double> preempt_cqn2(SchedStrategy s) {
    qn::Network<double> m("ssa-serial-preempt2");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", s);
    // A LOWER classprio value is the more urgent group, so C1 preempts C2. The
    // populations put TWO jobs in the buffer at once, without which the group
    // scan has a single candidate and the priority rule is never exercised.
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d, 0);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d, 1);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::erlang(4.0, 2));
    m.set_service(q, c2, Dist::erlang(6.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("SolverSSA serial: the preempt buffer enumerates [class, phase] pairs") {
    // The preempt family records the phase every WAITING job was interrupted
    // in, so a waiting position costs TWO columns and one marginal yields one
    // state per assignment of a phase to each waiting job -- this is what
    // `fromMarginal` builds by interleaving `mi_buf` with `bkstate`
    // (fromMarginal.m:310-325), and what `to_marginal`'s paired branch and the
    // departure handler's `colfirstnnz` both already assume.
    qn::Network<double> m = preempt_cqn(SchedStrategy::LCFSPR, 2.0);
    const std::vector<std::vector<double>> rows =
        qn::from_marginal_node(m.get_struct(), 2, std::vector<std::size_t>(1, 2),
                               std::vector<std::size_t>(1, 2));
    REQUIRE(!rows.empty());
    // Two jobs, one server, Erlang-2: [class, phase | phase1, phase2].
    CHECK(rows[0].size() == 4);
    // One waiting job in either phase, times the phase of the job in service.
    CHECK(rows.size() == 4);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        CHECK(rows[i][0] == 1.0);                        // the waiting class tag
        CHECK((rows[i][1] == 1.0 || rows[i][1] == 2.0));  // its stored phase
        CHECK(rows[i][2] + rows[i][3] == 1.0);           // exactly one in service
    }
    // The empty station keeps the buffer EVEN, so a left pad preserves parity.
    const std::vector<std::vector<double>> idle =
        qn::from_marginal_node(m.get_struct(), 2, std::vector<std::size_t>(1, 0),
                               std::vector<std::size_t>(1, 2));
    REQUIRE(idle.size() == 1);
    CHECK(idle[0].size() == 4);
}

TEST_CASE("DEP at LCFSPR resumes the stored phase, at LCFSPI restarts from pie") {
    // THE PR-VS-PI DISTINCTION, pinned on the DEPARTURE HANDLER itself, which is
    // exact and deterministic and so a sharper oracle than any sample path.
    //
    // Two class-1 jobs wait, the newest interrupted in phase 2 and the oldest in
    // phase 1: buffer [1,2 | 1,1], with a third job in service in phase 2. Only
    // phase 2 completes, since Erlang-2 has phi = [0,1].
    const std::vector<double> in{1.0, 2.0, 1.0, 1.0, 0.0, 1.0};

    qn::Network<double> mpr = preempt_cqn(SchedStrategy::LCFSPR, 3.0);
    const qn::EventOutcome<double> pr =
        qn::after_event_station_dep(mpr.get_struct(), 2, in, 1);
    // LCFS promotes the NEWEST, the FIRST nonzero pair (colfirstnnz, :1052), so
    // the OLDEST pair [1,1] is what stays, shifted right to keep the buffer
    // right-aligned. It resumes in phase 2, the phase it was stored in.
    REQUIRE(pr.space.size() == 1);
    CHECK(pr.space[0][0] == doctest::Approx(0.0));
    CHECK(pr.space[0][1] == doctest::Approx(0.0));
    CHECK(pr.space[0][2] == doctest::Approx(1.0));
    CHECK(pr.space[0][3] == doctest::Approx(1.0));
    CHECK(pr.space[0][4] == doctest::Approx(0.0));  // NOT phase 1
    CHECK(pr.space[0][5] == doctest::Approx(1.0));  // resumed in phase 2
    CHECK(pr.rate[0] == doctest::Approx(4.0));      // mu*phi*kir, phase rate 4

    qn::Network<double> mpi = preempt_cqn(SchedStrategy::LCFSPI, 3.0);
    const qn::EventOutcome<double> pi =
        qn::after_event_station_dep(mpi.get_struct(), 2, in, 1);
    // PI discards that phase and restarts from pie, which for Erlang-2 is
    // [1,0]: the same promoted job, the same buffer, but phase 1 and not 2.
    REQUIRE(pi.space.size() == 2);
    CHECK(pi.space[0][2] == doctest::Approx(1.0));
    CHECK(pi.space[0][3] == doctest::Approx(1.0));
    CHECK(pi.space[0][4] == doctest::Approx(1.0));  // restarted in phase 1
    CHECK(pi.space[0][5] == doctest::Approx(0.0));
    CHECK(pi.rate[0] == doctest::Approx(4.0));
    // A phase pie cannot reach carries rate 0 rather than being dropped.
    CHECK(pi.space[1][5] == doctest::Approx(1.0));
    CHECK(pi.rate[1] == doctest::Approx(0.0));
}

TEST_CASE("DEP at FCFSPR takes the longest waiting job, not the newest") {
    // The mirror image: the newest is stored in phase 1 and the oldest in phase
    // 2, so promoting the OLDEST (colLastNnz, :1121) both leaves the NEWEST pair
    // [1,1] in the buffer and resumes in phase 2. LCFS would do the opposite on
    // this same row, which is what makes the two rules distinguishable here.
    const std::vector<double> in{1.0, 1.0, 1.0, 2.0, 0.0, 1.0};

    qn::Network<double> mpr = preempt_cqn(SchedStrategy::FCFSPR, 3.0);
    const qn::EventOutcome<double> pr =
        qn::after_event_station_dep(mpr.get_struct(), 2, in, 1);
    REQUIRE(pr.space.size() == 1);
    CHECK(pr.space[0][0] == doctest::Approx(0.0));
    CHECK(pr.space[0][1] == doctest::Approx(0.0));
    CHECK(pr.space[0][2] == doctest::Approx(1.0));
    CHECK(pr.space[0][3] == doctest::Approx(1.0));  // the NEWEST pair survives
    CHECK(pr.space[0][4] == doctest::Approx(0.0));
    CHECK(pr.space[0][5] == doctest::Approx(1.0));  // resumed in phase 2
    CHECK(pr.rate[0] == doctest::Approx(4.0));

    qn::Network<double> mpi = preempt_cqn(SchedStrategy::FCFSPI, 3.0);
    const qn::EventOutcome<double> pi =
        qn::after_event_station_dep(mpi.get_struct(), 2, in, 1);
    REQUIRE(pi.space.size() == 2);
    CHECK(pi.space[0][3] == doctest::Approx(1.0));
    CHECK(pi.space[0][4] == doctest::Approx(1.0));  // restarted in phase 1
    CHECK(pi.rate[0] == doctest::Approx(4.0));
    CHECK(pi.rate[1] == doctest::Approx(0.0));

    // Same row under LCFSPR: the OTHER pair survives and the server takes the
    // OTHER phase, so the two selection rules cannot be collapsed.
    qn::Network<double> ml = preempt_cqn(SchedStrategy::LCFSPR, 3.0);
    const qn::EventOutcome<double> lc =
        qn::after_event_station_dep(ml.get_struct(), 2, in, 1);
    REQUIRE(lc.space.size() == 1);
    CHECK(lc.space[0][3] == doctest::Approx(2.0));  // the OLDEST pair survives
    CHECK(lc.space[0][4] == doctest::Approx(1.0));  // resumed in phase 1
}

TEST_CASE("DEP at the PRIO variants promotes the most urgent group") {
    // Two waiting jobs, the less urgent C2 NEWEST and the more urgent C1 OLDEST:
    // buffer [2,1 | 1,1], one C1 job in service in phase 2. A LOWER classprio
    // value is more urgent, so all four PRIO variants must skip the newest and
    // promote C1, where plain LCFS would take C2. The scan must read only the
    // CLASS columns; over all columns a stored phase of 1 would read as class 1
    // and win the group, which is what the reference's PI-PRIO arms do.
    const std::vector<double> in{2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0};
    const SchedStrategy prio[] = {SchedStrategy::LCFSPRPRIO, SchedStrategy::LCFSPIPRIO,
                                  SchedStrategy::FCFSPRPRIO, SchedStrategy::FCFSPIPRIO};
    for (int t = 0; t < 4; ++t) {
        qn::Network<double> m = preempt_cqn2(prio[t]);
        const qn::EventOutcome<double> o =
            qn::after_event_station_dep(m.get_struct(), 2, in, 1);
        INFO("discipline index ", t);
        REQUIRE(!o.space.empty());
        // The C2 pair is what remains, shifted right by a whole empty pair.
        CHECK(o.space[0][0] == doctest::Approx(0.0));
        CHECK(o.space[0][1] == doctest::Approx(0.0));
        CHECK(o.space[0][2] == doctest::Approx(2.0));
        CHECK(o.space[0][3] == doctest::Approx(1.0));
        // The server block is [C1 ph1, C1 ph2, C2 ph1, C2 ph2]: C1 completed and
        // the promoted job is C1 again, so C2 must NOT be in service.
        CHECK(o.space[0][4] + o.space[0][5] == doctest::Approx(1.0));
        CHECK(o.space[0][6] == doctest::Approx(0.0));
        CHECK(o.space[0][7] == doctest::Approx(0.0));
        CHECK(o.rate[0] == doctest::Approx(4.0));
    }
}

TEST_CASE("SolverSSA serial: every preempt discipline reaches the engine") {
    // All eight members, each on the model its own selection rule exercises: the
    // four PRIO variants need two classes at distinct priorities, since with one
    // class the group scan is vacuous.
    //
    // The exact CTMC is the oracle for the three the CTMC feature set admits
    // (SolverCTMC.m:159 lists LCFSPR, LCFSPRPRIO and FCFSPRPRIO and no other
    // member, a gap in the reference's own declaration since afterEventStation
    // handles all eight). For the remaining five the oracle is the population
    // identity, which holds for EVERY sample and so exactly for the time
    // average: a promotion rule that lost or duplicated a job would break it.
    const SchedStrategy plain[] = {SchedStrategy::LCFSPR, SchedStrategy::LCFSPI,
                                   SchedStrategy::FCFSPR, SchedStrategy::FCFSPI};
    const SchedStrategy prio[] = {SchedStrategy::LCFSPRPRIO, SchedStrategy::LCFSPIPRIO,
                                  SchedStrategy::FCFSPRPRIO, SchedStrategy::FCFSPIPRIO};
    for (int t = 0; t < 8; ++t) {
        const bool isprio = t >= 4;
        const SchedStrategy s = isprio ? prio[t - 4] : plain[t];
        qn::Network<double> m = isprio ? preempt_cqn2(s) : preempt_cqn(s, 3.0);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        const std::size_t R = isprio ? 2 : 1;
        const bool ctmc_ok = s == SchedStrategy::LCFSPR || s == SchedStrategy::LCFSPRPRIO ||
                             s == SchedStrategy::FCFSPRPRIO;

        ssa::SsaSerialOptions opt;
        opt.samples = 100000;
        opt.seed = 7 + t;
        const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);
        INFO("discipline index ", t);
        for (std::size_t r = 0; r < R; ++r) {
            const double pop = isprio ? (r == 0 ? 2.0 : 1.0) : 3.0;
            CHECK(sim.avg.QN(0, r) + sim.avg.QN(1, r) == doctest::Approx(pop).epsilon(1e-9));
            CHECK(sim.avg.TN(1, r) > 0.0);
        }
        if (!ctmc_ok) continue;
        const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());
        const double band = 0.12;  // 1e5 firings
        for (std::size_t r = 0; r < R; ++r) {
            CHECK(std::isfinite(exact.QN(1, r)));
            CHECK(sim.avg.QN(1, r) == doctest::Approx(exact.QN(1, r)).epsilon(band));
            CHECK(sim.avg.TN(1, r) == doctest::Approx(exact.TN(1, r)).epsilon(band));
        }
    }
}

/*
 * A CACHE READ IS READ FROM THE CACHE, not drawn from a coin.
 *
 * `refresh_routing` resolves the cache's unresolved hit/miss split to a uniform
 * half-half self-loop so the visit equations have a number; until 2026-08-03
 * `refresh_sync` turned that resolved value into a DEPARTURE synchronization and
 * emitted no READ at all, so the sample path decided hit against miss by a coin
 * and every cache model reported one half whatever its popularity was. Two
 * items, room for one, popularity (0.9, 0.1): a read hits exactly when it
 * repeats the previous one, so the stationary hit probability is
 * 0.9^2 + 0.1^2 = 0.82 and NOT 0.5. The band is chosen far below the gap
 * between the two, which is what makes this a structural assertion rather than
 * a tolerance.
 */
TEST_CASE("ssa serial: a cache reads its contents, it does not flip a coin") {
    using line::lang::ReplacementStrategy;
    typedef line::lang::Distrib<double> D;
    qn::Network<double> m("cache-read");
    const std::size_t src = m.add_source("Source");
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t mis = m.add_open_class("Miss");

    qn::CacheParam<double> cp;
    cp.nitems = 2;
    cp.itemcap.push_back(1);
    cp.replacestrat = ReplacementStrategy::LRU;
    cp.pread.assign(3, std::vector<double>());
    cp.pread[rd - 1] = std::vector<double>{0.9, 0.1};
    cp.hitclass.assign(3, 0);
    cp.missclass.assign(3, 0);
    cp.hitclass[rd - 1] = hit;
    cp.missclass[rd - 1] = mis;
    const std::size_t ca = m.add_cache("C", cp);

    m.set_arrival(src, rd, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(rd, rd, src, ca, 1.0);
    P.set(hit, hit, ca, sk, 1.0);
    P.set(mis, mis, ca, sk, 1.0);
    m.link(P);

    ssa::SsaSerialOptions opt;
    opt.samples = 20000;
    opt.seed = 11;
    const ssa::SsaSerialSolution<double> sim =
        ssa::solver_ssa_serial_analyzer(m.get_struct(), opt);
    REQUIRE(sim.cache.size() == 1);
    CHECK(sim.cache[0].node == ca);
    CHECK(sim.cache[0].hitprob[rd - 1] + sim.cache[0].missprob[rd - 1] ==
          doctest::Approx(1.0).epsilon(1e-9));
    CHECK(sim.cache[0].hitprob[rd - 1] == doctest::Approx(0.82).epsilon(0.05));
    // A cache with no retrieval system merges nothing, so the field stays EMPTY
    // rather than reporting a zero share that was never measured.
    CHECK(sim.cache[0].delayedprob.empty());
}

/*
 * The delayed hit, measured on the sample path.
 *
 * Source -> Cache -> (miss) an infinite-server fetch station -> back to the
 * Cache. A read for an item already being fetched MERGES onto that fetch and is
 * released in the hit class when it completes, so the hit-class departure rate
 * carries true hits AND delayed hits and only the merge transition separates
 * them. `retrieval_simple` of the example suite, whose exact answer SolverNC
 * gives as hit 0.44605, delayed 0.091157, miss 0.46280.
 *
 * The band is 15% of the delayed share at 20000 firings, an order above the
 * run's standard error: what is asserted is that the share is MEASURED and
 * near the exact value, not that a sample path reproduces it.
 */
TEST_CASE("ssa serial: a delayed hit is measured, not folded into the hit rate") {
    using line::lang::ReplacementStrategy;
    typedef line::lang::Distrib<double> D;
    qn::Network<double> m("retrieval");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> cp;
    cp.nitems = 3;
    cp.itemcap.push_back(1);
    cp.replacestrat = ReplacementStrategy::FIFO;
    cp.pread = std::vector<std::vector<double> >{std::vector<double>{0.6, 0.3, 0.1}, {}, {}};
    cp.hitclass = std::vector<std::size_t>{2, 0, 0};
    cp.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t ca = m.add_cache("Cache", cp);
    const std::size_t q = m.add_queue("Queue", line::lang::SchedStrategy::INF);
    const std::size_t sk = m.add_sink("Sink");

    const std::size_t rd = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t mis = m.add_open_class("MissClass");
    m.set_arrival(src, rd, D::exp_rate(1.0));
    m.set_service(q, rd, D::exp_rate(2.0));
    m.set_retrieval_system(ca, rd, mis, std::vector<std::size_t>{q});

    qn::RoutingMatrix<double> P;
    P.set(rd, rd, src, ca, 1.0);
    P.set(rd, rd, ca, q, 1.0);
    P.set(rd, rd, q, ca, 1.0);
    P.set(hit, hit, ca, sk, 1.0);
    P.set(mis, mis, ca, sk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    ssa::SsaSerialOptions opt;
    opt.samples = 20000;
    opt.seed = 1;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);
    REQUIRE(sim.cache.size() == 1);
    REQUIRE(sim.cache[0].delayedprob.size() == sn.nclasses);
    const double h = sim.cache[0].hitprob[rd - 1];
    const double d = sim.cache[0].delayedprob[rd - 1];
    const double mi = sim.cache[0].missprob[rd - 1];
    // The three shares PARTITION every read: this is what makes ArvR the
    // retrieval flow `arvr*(missprob + delayedprob)` rather than `arvr*missprob`.
    CHECK(h + d + mi == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(d == doctest::Approx(0.091157).epsilon(0.15));
    CHECK(h == doctest::Approx(0.44605).epsilon(0.05));
    CHECK(mi == doctest::Approx(0.46280).epsilon(0.05));
    // The fetch station is ENTERED. Before the retrieval branches were ported a
    // miss switched straight into the miss class, so the retrieval sub-network
    // carried no jobs at all and this read zero.
    double fetch = 0.0;
    for (std::size_t k = 0; k < sn.nclasses; ++k) fetch += sim.avg.TN(sn.nodes[q - 1].station - 1, k);
    CHECK(fetch > 0.4);
}
