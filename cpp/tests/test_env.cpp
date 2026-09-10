/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverENV: a network in a random environment, with fluid stages.
 *
 * Model: renv_basic (matlab/examples/advanced/randomEnv/renv_basic.m). A closed
 * network of a think-time Delay and an FCFS queue, five jobs, whose server runs
 * at rate 4 in the Fast stage and rate 1 in the Slow stage. Fast -> Slow fires
 * at rate 0.5 and Slow -> Fast at rate 1, so the environment sits in Fast two
 * thirds of the time. The reference spelling is
 * `ENV(env, @(m) FLD(m, fldOptions), options)` with `fldOptions.timespan =
 * [0, 100]`.
 *
 * WHAT IS ASSERTED, in three kinds:
 *
 *  1. The environment process itself -- probEnv, probOrig and the holding-time
 *     means -- against MATLAB. These are exact and carry no quadrature.
 *
 *  2. The blended metrics against MATLAB's, at a tolerance that reflects the
 *     one real numerical difference left between the two ports: the exit
 *     metrics are a Riemann-Stieltjes sum, and since 2026-08-11 all four
 *     codebases sum it on the SAME sojourn-scaled grid, leaving only the spread
 *     between their integrators' own trajectories on that grid.
 *
 *  3. Properties neither codebase can influence: the population is conserved
 *     exactly, the saturated Slow stage pins the server at utilization one, and
 *     the answer does NOT move with the caller's `tran_points`, which is what
 *     says the sum is taken on the sojourn grid rather than on whatever the
 *     integrator happened to report.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/env/solver_env.h"

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

}  // namespace

TEST_CASE("renv_basic under SolverENV with fluid stages, against MATLAB ENV(FLD)") {
    qn::Network<double> fast = make_stage("Fast", 4.0);
    qn::Network<double> slow = make_stage("Slow", 1.0);
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, D::exp_rate(0.5));
    e.add_transition(1, 0, D::exp_rate(1.0));

    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    o.tran_points = 2001;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();

    // 1. The environment process, exactly as MATLAB's env.init() reports it.
    CHECK(e.prob_env[0] == doctest::Approx(0.666666666667).epsilon(1e-9));
    CHECK(e.prob_env[1] == doctest::Approx(0.333333333333).epsilon(1e-9));
    CHECK(e.prob_orig(0, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(e.prob_orig(0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(e.prob_orig(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(e.prob_orig(1, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(mam::map_mean(e.hold_time[0].map()) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(mam::map_mean(e.hold_time[1].map()) == doctest::Approx(1.0).epsilon(1e-9));

    // 2. MATLAB ENV(FLD).getAvg on this model: QN = [2.99611792893, 2.00388207107],
    // UN = [2.99611792893, 1]. The 1e-2 is the quadrature gap described above.
    REQUIRE(sol.converged);
    CHECK(sol.QN(0, 0) == doctest::Approx(2.99611792893).epsilon(1e-2));
    CHECK(sol.QN(1, 0) == doctest::Approx(2.00388207107).epsilon(1e-2));
    CHECK(sol.UN(0, 0) == doctest::Approx(2.99611792893).epsilon(1e-2));
    CHECK(sol.UN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));

    // 3a. The population is conserved exactly: no job is created or lost at a
    // switch, whatever the quadrature does.
    CHECK(sol.QN(0, 0) + sol.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
    // The Delay is an infinite server, so its "utilization" is its queue length.
    CHECK(sol.UN(0, 0) == doctest::Approx(sol.QN(0, 0)).epsilon(1e-9));
    // The throughput blend is exact here and MATLAB does not report it (its FLD
    // transient leaves the T handle empty, so ENV blends NaN): in Fast the
    // server clears 4 per unit time and in Slow 1, both saturated, so the
    // environment average is 2/3 * 4 + 1/3 * 1 = 3.
    CHECK(sol.TN(1, 0) == doctest::Approx(3.0).epsilon(1e-6));

    // 3b. THE CALLER'S GRID NO LONGER MOVES THE ANSWER, and that is the point of
    // the 2026-08-11 alignment: the stage transient is asked for on the SOJOURN
    // grid (`stage_grid`, 5000 points scaled to E[S]), not on the caller's
    // output grid, so `tran_points` seeds `pre()` and nothing else. All four
    // codebases rebuild that grid the same way, which is what makes them sum the
    // same points; a run that still drifted with `tran_points` would mean the
    // refined grid had been dropped somewhere and the exit average had gone back
    // to reading the integrator's own steps.
    double prev = 0.0;
    const std::size_t grids[] = {101, 501, 2001, 10001};
    for (std::size_t np : grids) {
        env::EnvOptions og = o;
        og.tran_points = np;
        env::SolverEnv<double> sg(e, og);
        const env::EnvSolution rg = sg.solve();
        CAPTURE(np);
        if (prev > 0.0) CHECK(rg.QN(0, 0) == doctest::Approx(prev).epsilon(1e-9));
        CHECK(rg.QN(0, 0) + rg.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
        prev = rg.QN(0, 0);
    }
    // and the value that grid converges to, against MATLAB's 2.99611792893
    CHECK(prev == doctest::Approx(2.9994).epsilon(1e-3));

    // The deterministic sojourn evaluates each stage at its MEAN holding time
    // instead of averaging over it, so it must give a different answer while
    // still conserving the population.
    env::EnvOptions od = o;
    od.sojourn = "deterministic";
    env::SolverEnv<double> sd(e, od);
    const env::EnvSolution rd = sd.solve();
    CHECK(rd.QN(0, 0) + rd.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
    CHECK(rd.QN(0, 0) != doctest::Approx(sol.QN(0, 0)).epsilon(1e-4));
}

TEST_CASE("SolverENV refuses by name what it does not implement") {
    qn::Network<double> fast = make_stage("Fast", 4.0);
    qn::Network<double> slow = make_stage("Slow", 1.0);
    auto build = [&]() {
        env::Environment<double>* e = new env::Environment<double>("ServerModes", 2);
        e->set_stage(0, "Fast", "operational", fast.get_struct());
        e->set_stage(1, "Slow", "degraded", slow.get_struct());
        e->add_transition(0, 1, D::exp_rate(0.5));
        e->add_transition(1, 0, D::exp_rate(1.0));
        return e;
    };
    // Every branch of the reference's ladder this class does not serve, by the
    // name the reference gives it. `smp` and `statedep` are NOT here: this
    // class serves both, and a `statedep` run with no arc hook is refused as an
    // InputError further down rather than as an unported method.
    const char* methods[] = {"statevec", "blend", "avg", "dec", "nonsense"};
    for (const char* mth : methods) {
        CAPTURE(mth);
        env::Environment<double>* e = build();
        env::EnvOptions o;
        o.method = mth;
        CHECK_THROWS_AS(env::SolverEnv<double>(*e, o), UnsupportedError);
        delete e;
    }
    {   // Only the fluid analyzer produces the transient the coupling needs.
        env::Environment<double>* e = build();
        env::EnvOptions o;
        o.stage_solver = "mva";
        CHECK_THROWS_AS(env::SolverEnv<double>(*e, o), UnsupportedError);
        delete e;
    }
    {   // A stage transient with no horizon is an input error, as in the
        // reference's own check on options.timespan(2).
        env::Environment<double>* e = build();
        env::EnvOptions o;
        o.timespan_end = 0.0;
        CHECK_THROWS_AS(env::SolverEnv<double>(*e, o), InputError);
        o.timespan_end = 100.0;
        o.tran_points = 1;
        CHECK_THROWS_AS(env::SolverEnv<double>(*e, o), InputError);
        delete e;
    }
    {   // A stage with no model at all cannot be solved.
        env::Environment<double> e("half", 2);
        e.set_stage(0, "Fast", "operational", fast.get_struct());
        e.add_transition(0, 1, D::exp_rate(0.5));
        env::EnvOptions o;
        CHECK_THROWS_AS(env::SolverEnv<double>(e, o), InputError);
    }
}

TEST_CASE("SolverENV reset policies are applied at the switch") {
    // A reset that FLUSHES the queue on every switch cannot leave the
    // population at five: the jobs that were in the network are discarded, so
    // the entry vector, and with it the blend, must fall.
    qn::Network<double> fast = make_stage("Fast", 4.0);
    qn::Network<double> slow = make_stage("Slow", 1.0);
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    const env::ResetMarginal flush_queue = [](const Matrix<double>& Q) {
        Matrix<double> R = Q;
        for (std::size_t r = 0; r < R.cols(); ++r) R(1, r) = 0.0;  // the queue empties
        return R;
    };
    e.add_transition(0, 1, D::exp_rate(0.5), flush_queue);
    e.add_transition(1, 0, D::exp_rate(1.0), flush_queue);

    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();
    for (std::size_t st = 0; st < 2; ++st)
        CHECK(sol.Qentry[st](1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sol.QN(0, 0) + sol.QN(1, 0) < 5.0);
}

TEST_CASE("SolverENV statedep rewrites the environment from the state it is left in") {
    // `resetEnvRates` in the reference: the breakdown becomes likelier the more
    // work is queued at the server. The environment process is therefore not an
    // input but part of the fixed point -- probEnv is whatever the converged
    // load makes it.
    qn::Network<double> fast = make_stage("Fast", 4.0);
    qn::Network<double> slow = make_stage("Slow", 1.0);
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, D::exp_rate(0.5));
    e.add_transition(1, 0, D::exp_rate(1.0));
    e.set_env_rate_reset(0, 1,
                         [](const D&, const Matrix<double>& Q, const Matrix<double>&,
                            const Matrix<double>&) { return D::exp_rate(0.05 * (1.0 + Q(1, 0))); });

    env::EnvOptions o;
    o.method = "statedep";
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();

    // The arc it was given is not the arc it ended with, and the stage
    // probabilities followed: the Fast stage now lasts as long as its own load
    // says it does, which is longer than the rate 0.5 it started from.
    CHECK(e.arc(0, 1).dist.mean != doctest::Approx(2.0).epsilon(1e-6));
    CHECK(e.arc(1, 0).dist.mean == doctest::Approx(1.0).epsilon(1e-12));  // untouched
    CHECK(e.prob_env[0] != doctest::Approx(0.666666666667).epsilon(1e-6));
    CHECK(e.prob_env[0] + e.prob_env[1] == doctest::Approx(1.0).epsilon(1e-9));

    // The rate is 0.05 * (1 + Qexit), so the converged arc and the converged
    // queue length are each other's: this is the fixed point, read back.
    const double q = sol.QExit[0](1, 0);
    CHECK(q > 0.0);
    CHECK(sol.QN(0, 0) + sol.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));

    // Without a hook the method is a mean-field run under another name, and
    // saying so is better than reporting one as the other.
    qn::Network<double> f2 = make_stage("Fast", 4.0);
    qn::Network<double> s2 = make_stage("Slow", 1.0);
    env::Environment<double> plain("ServerModes", 2);
    plain.set_stage(0, "Fast", "operational", f2.get_struct());
    plain.set_stage(1, "Slow", "degraded", s2.get_struct());
    plain.add_transition(0, 1, D::exp_rate(0.5));
    plain.add_transition(1, 0, D::exp_rate(1.0));
    CHECK_THROWS_AS(env::SolverEnv<double>(plain, o), InputError);
}

TEST_CASE("SolverENV accepts smp, which selects no analyzer of its own") {
    // In the reference `smp` only lifts the constructor's check that every arc
    // is Markovian; the analyzer that runs is the mean-field one. An arc here is
    // a (D0, D1) pair by construction, so the two methods must agree exactly.
    qn::Network<double> fast = make_stage("Fast", 4.0);
    qn::Network<double> slow = make_stage("Slow", 1.0);
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "Fast", "operational", fast.get_struct());
    e.set_stage(1, "Slow", "degraded", slow.get_struct());
    e.add_transition(0, 1, D::erlang(1.0, 2));  // a non-exponential holding time
    e.add_transition(1, 0, D::exp_rate(1.0));

    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    env::SolverEnv<double> mf(e, o);
    const env::EnvSolution a = mf.solve();
    env::EnvOptions osmp = o;
    osmp.method = "smp";
    env::SolverEnv<double> smp(e, osmp);
    const env::EnvSolution b = smp.solve();
    CHECK(a.QN(0, 0) == doctest::Approx(b.QN(0, 0)).epsilon(1e-12));
    CHECK(a.QN(1, 0) == doctest::Approx(b.QN(1, 0)).epsilon(1e-12));
}
