/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The closed-form environment limits, `SolverENV.solveEnvLimit`.
 *
 * Model: renv_basic again (a think-time Delay and an FCFS queue, five jobs, a
 * server that runs fast in one stage and slow in the other), so the environment
 * process is the one `test_env.cpp` already pinned against MATLAB: probEnv is
 * 2/3 and 1/3.
 *
 * WHAT IS ASSERTED, and why it is not a tautology. Each limit is DEFINED as a
 * different model solved in steady state, so the test builds that model itself
 * and solves it with the same fluid analyzer:
 *
 *   `dec` must equal the probEnv blend of the two stages solved separately;
 *   `avg` must equal ONE network whose server rate is 2/3 * fast + 1/3 * slow,
 *   built here from scratch rather than derived from the environment.
 *
 * The second is the sharper of the two: it checks that
 * `buildRateAveragedModel` averaged the RATES (and left the unmodulated think
 * time alone) rather than averaging the answers, which is the whole difference
 * between the fast-environment limit and the slow one.
 *
 * On a model whose stages both saturate the server, the two limits coincide
 * exactly -- the fluid throughput is linear in the rate there -- so the
 * separating case uses a fast stage that does NOT saturate, where they must
 * disagree.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/env/env_dispatch.h"
#include "line/solvers/env/solver_env_limit.h"
#include "line/solvers/fluid/solver_fluid.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

qn::Network<double> make_stage(const std::string& nm, double svc) {
    qn::Network<double> m(nm);
    const std::size_t d = m.add_delay("ThinkTime");
    const std::size_t q = m.add_queue("Fast/Slow Server", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Jobs", 5, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(svc));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

/** The two-stage environment of renv_basic, with the server rates given. */
env::Environment<double> make_env(double fast_rate, double slow_rate) {
    qn::Network<double> fast = make_stage("Fast", fast_rate);
    qn::Network<double> slow = make_stage("Slow", slow_rate);
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, D::exp_rate(0.5));
    e.add_transition(1, 0, D::exp_rate(1.0));
    return e;
}

/** One network in steady state, with the analyzer the limits use. */
fluid::FluidSolution steady(const qn::NetworkStruct<double>& sn) {
    fluid::FluidOptions fo;
    return fluid::solver_fluid(sn, fo);
}

}  // namespace

TEST_CASE("the dec limit is the probEnv blend of the stages solved separately") {
    env::Environment<double> e = make_env(4.0, 1.0);
    env::EnvOptions o;
    o.method = "dec";
    const env::EnvLimitSolution r = env::solver_env_limit(e, o);

    REQUIRE(r.method == "dec");
    REQUIRE(r.prob_env.size() == 2);
    CHECK(r.prob_env[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-9));

    const fluid::FluidSolution fast = steady(e.stage(0).model);
    const fluid::FluidSolution slow = steady(e.stage(1).model);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.QStage[0](i, 0) == doctest::Approx(fast.QN(i, 0)).epsilon(1e-9));
        CHECK(r.QStage[1](i, 0) == doctest::Approx(slow.QN(i, 0)).epsilon(1e-9));
        CHECK(r.QN(i, 0) == doctest::Approx(2.0 / 3.0 * fast.QN(i, 0) + slow.QN(i, 0) / 3.0)
                                .epsilon(1e-9));
        CHECK(r.TN(i, 0) == doctest::Approx(2.0 / 3.0 * fast.TN(i, 0) + slow.TN(i, 0) / 3.0)
                                .epsilon(1e-9));
    }

    // Both stages saturate their server (5 jobs against a think time of 1), so
    // each stage clears jobs at its own service rate and the blend is exact:
    // 2/3 * 4 + 1/3 * 1 = 3, with one job at the queue in Fast and four in Slow.
    CHECK(r.TN(1, 0) == doctest::Approx(3.0).epsilon(1e-6));
    CHECK(r.UN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("the avg limit is one network at the probEnv-weighted rate") {
    env::Environment<double> e = make_env(4.0, 1.0);
    env::EnvOptions o;
    o.method = "avg";
    const env::EnvLimitSolution r = env::solver_env_limit(e, o);
    REQUIRE(r.method == "avg");

    // The model the limit claims to solve, built from scratch: the server runs
    // at 2/3 * 4 + 1/3 * 1 = 3 and the think time, which no stage modulates,
    // keeps its own distribution.
    qn::Network<double> avg = make_stage("Averaged", 3.0);
    const fluid::FluidSolution ref = steady(avg.get_struct());
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.QN(i, 0) == doctest::Approx(ref.QN(i, 0)).epsilon(1e-9));
        CHECK(r.UN(i, 0) == doctest::Approx(ref.UN(i, 0)).epsilon(1e-9));
        CHECK(r.TN(i, 0) == doctest::Approx(ref.TN(i, 0)).epsilon(1e-9));
    }
    CHECK(r.QStage.empty());  // one model was solved, not one per stage
}

TEST_CASE("the two limits disagree when a stage does not saturate") {
    // A fast stage of rate 10 leaves the server idle part of the time, so the
    // throughput is no longer linear in the rate and averaging the RATES is no
    // longer the same as averaging the ANSWERS. If these two agreed, the avg
    // arm would be silently running the dec one.
    env::Environment<double> e = make_env(10.0, 1.0);
    env::EnvOptions oa, od;
    oa.method = "avg";
    od.method = "dec";
    const env::EnvLimitSolution a = env::solver_env_limit(e, oa);
    env::Environment<double> e2 = make_env(10.0, 1.0);
    const env::EnvLimitSolution d = env::solver_env_limit(e2, od);

    qn::Network<double> avg = make_stage("Averaged", 2.0 / 3.0 * 10.0 + 1.0 / 3.0);
    const fluid::FluidSolution ref = steady(avg.get_struct());
    CHECK(a.QN(1, 0) == doctest::Approx(ref.QN(1, 0)).epsilon(1e-9));
    CHECK(a.QN(1, 0) != doctest::Approx(d.QN(1, 0)).epsilon(1e-3));
    // Both still conserve the population, which is a property of neither limit
    // but of the model they solve.
    CHECK(a.QN(0, 0) + a.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
    CHECK(d.QN(0, 0) + d.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("a rate that is not modulated keeps its own distribution under avg") {
    // The think time is an Erlang here and the same in both stages. The
    // reference rewrites only MODULATED rates as exponentials, so this one must
    // survive as an Erlang -- an exponential of the same mean would change the
    // stage transient of every solver that reads the phase structure.
    qn::Network<double> fast("Fast"), slow("Slow");
    for (int pass = 0; pass < 2; ++pass) {
        qn::Network<double>& m = (pass == 0) ? fast : slow;
        const std::size_t d = m.add_delay("ThinkTime");
        const std::size_t q = m.add_queue("Server", lang::SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("Jobs", 5, d);
        m.set_service(d, c, D::erlang(2.0, 2));
        m.set_service(q, c, D::exp_rate(pass == 0 ? 4.0 : 1.0));
        qn::RoutingMatrix<double> P;
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
        m.link(P);
    }
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, D::exp_rate(0.5));
    e.add_transition(1, 0, D::exp_rate(1.0));

    env::EnvOptions o;
    o.method = "avg";
    const env::EnvLimitSolution r = env::solver_env_limit(e, o);
    // The Erlang think time has squared coefficient of variation 1/2; the
    // averaged model must still show it, and the modulated server must not.
    CHECK(e.stage(0).model.scv(0, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("the limits are reached through the ENV dispatch, and refuse what they cannot do") {
    env::Environment<double> e = make_env(4.0, 1.0);
    env::EnvOptions o;
    o.method = "dec";
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
    CHECK(r.method == "dec");
    CHECK(r.converged);      // a closed form has converged by construction
    CHECK(r.iterations == 0);  // and iterated nothing
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));

    // The mean-field coupling refuses a limit BY NAME rather than serving it:
    // the two carry different things across a switch, and a limit carries none.
    env::Environment<double> e2 = make_env(4.0, 1.0);
    CHECK_THROWS_AS(env::SolverEnv<double>(e2, o), UnsupportedError);

    // And the limit refuses a method that is not one.
    env::EnvOptions bad = o;
    bad.method = "meanfield";
    env::Environment<double> e3 = make_env(4.0, 1.0);
    CHECK_THROWS_AS(env::solver_env_limit(e3, bad), UnsupportedError);
}
