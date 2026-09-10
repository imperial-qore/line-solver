/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverENV over LAYERED stages: a LayeredNetwork in a random environment.
 *
 * Model: renv_lqn_twostages (matlab/examples/advanced/randomEnv/renv_lqn_twostages.m
 * and its python twin). A five-job client task with a 5.0 think time calls a
 * database entry 2.5 times per request; the environment alternates between an
 * UP stage whose DB activity has host demand 0.8 and a DOWN stage where it has
 * 3.0. UP -> DOWN fires at rate 0.2 and DOWN -> UP at rate 1.0, so the
 * environment sits in UP five sixths of the time.
 *
 * WHAT AN LQN STAGE CHANGES, and it is only this: the stage has no stations of
 * its own. SolverLN builds one network per layer, and the (station, class) view
 * the coupling blends over is the BLOCK-DIAGONAL UNION of those layers --
 * `LayeredNetwork.layerBlocks` in the reference. The handoff across a switch is
 * the same marginal mean queue lengths as ever, split into per-layer blocks on
 * the way in and reassembled on the way out.
 *
 * WHAT IS ASSERTED, and why properties rather than a golden row:
 *
 *  1. THE HANDOFF IS NOT INERT. This is the sharp one, and it is the defect
 *     both the JAR and python hit: the layered fixed point RESETS its layers as
 *     it converges, so a warm start installed before the solve is gone by the
 *     time the transients run. A silently inert handoff still converges and
 *     still prints a table -- of warm-up averages, each stage restarted empty --
 *     so nothing but reading the transient's own first point catches it.
 *
 *  2. THE THROUGHPUT IS BRACKETED BY THE TWO STAGES, and RISES with the
 *     stationary probability of the fast one. That is the reference's own
 *     oracle, and it is what an inert handoff breaks: warm-up averages fall
 *     BELOW the slowest single-stage value rather than between the two.
 *
 *  3. THE POPULATION IS CONSERVED. Every layer carries the same five client
 *     jobs, so the aggregate sums five per layer whatever the environment does.
 *
 *  4. THE COUPLINGS THAT CANNOT TAKE A LAYERED STAGE REFUSE ONE BY NAME. The
 *     state-vector coupling needs one generator per stage and the closed-form
 *     limits need a station rate table; an LQN has neither, and the reference
 *     refuses the first in exactly the same place.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/env/env_dispatch.h"

using namespace line;
using D = lang::Distrib<double>;
using SS = lang::SchedStrategy;

namespace {

/** The stage LQN; `db_mean` is the only thing that differs between stages. */
lqn::LqnStruct<double> lqn_stage(double db_mean) {
    lqn::LqnBuilder<double> b;
    b.processor("ClientProcessor", 1, SS::PS);
    b.processor("DBProcessor", 1, SS::PS);
    b.task("ClientTask", 5, SS::REF, "ClientProcessor");
    b.think_time("ClientTask", D::exp_mean(5.0));
    b.task("DBTask", std::numeric_limits<double>::infinity(), SS::INF, "DBProcessor");
    b.entry("ClientEntry", "ClientTask");
    b.entry("DBEntry", "DBTask");
    b.activity("ClientActivity", D::exp_mean(1.0), "ClientTask");
    b.bound_to("ClientActivity", "ClientEntry");
    b.sync_call("ClientActivity", "DBEntry", 2.5);
    b.activity("DBActivity", D::exp_mean(db_mean), "DBTask");
    b.bound_to("DBActivity", "DBEntry");
    b.replies_to("DBActivity", "DBEntry");
    return b.build();
}

/** `EnvOptions` as the reference's `lnFactory` spells it: fluid layers over [0, T]. */
env::EnvOptions lqn_env_options(double T, int iter_max, double iter_tol) {
    env::EnvOptions o;
    o.method = "meanfield";
    o.timespan_end = T;
    o.iter_max = iter_max;
    o.iter_tol = iter_tol;
    o.lqn.layer_solver = "fluid";
    return o;
}

/** The two-stage UP/DOWN environment at the given switch rates. */
env::Environment<double> lqn_env(double a, double b) {
    env::Environment<double> e("DBReliability", 2);
    e.set_lqn_stage(0, "UP", "operational", lqn_stage(0.8));
    e.set_lqn_stage(1, "DOWN", "degraded", lqn_stage(3.0));
    e.add_transition(0, 1, D::exp_rate(a));
    e.add_transition(1, 0, D::exp_rate(b));
    return e;
}

/** `sum(TN(isfinite(TN)))`, the reference's total over the aggregate. */
double finite_sum(const Matrix<double>& A) {
    double s = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            if (std::isfinite(A(i, j))) s += A(i, j);
    return s;
}

/** `stageAggregate`: one LQN alone, at the END of the same layered transient. */
void stage_alone(double db_mean, double T, double& sumQ, double& sumT) {
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = T;
    ln::SolverLN<double> s(lqn_stage(db_mean), lo);
    const ln::LnTranSolution tr = s.get_tran_avg();
    sumQ = 0.0;
    sumT = 0.0;
    for (std::size_t l = 0; l < tr.layers.size(); ++l) {
        const ln::LnTranLayer& L = tr.layers[l];
        if (L.t.empty()) continue;
        const std::size_t last = L.t.size() - 1;
        for (std::size_t i = 0; i < L.QN.size(); ++i)
            for (std::size_t r = 0; r < L.QN[i].size(); ++r) {
                sumQ += L.QN[i][r][last];
                sumT += L.TN[i][r][last];
            }
    }
}

/** The environment-averaged total throughput of one two-stage instance. */
double env_tput(double a, double b, double T) {
    env::Environment<double> e = lqn_env(a, b);
    return finite_sum(env::solver_env(e, lqn_env_options(T, 10, 0.03)).TN);
}

}  // namespace

TEST_CASE("layerBlocks tiles the aggregate view of a layered stage") {
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = 50.0;
    ln::SolverLN<double> s(lqn_stage(0.8), lo);
    const ln::LnLayerBlocks b = s.layer_blocks();

    REQUIRE(s.nlayers() > 0);
    REQUIRE(b.roff.size() == s.nlayers());
    // The blocks tile: each starts where the previous ended, and the totals are
    // the sums. Nothing overlaps and nothing is left out, which is what lets an
    // aggregate marginal be split and reassembled without loss.
    std::size_t M = 0, K = 0;
    for (std::size_t l = 0; l < s.nlayers(); ++l) {
        CHECK(b.roff[l] == M);
        CHECK(b.coff[l] == K);
        CHECK(b.msz[l] == s.layers()[l].nstations);
        CHECK(b.ksz[l] == s.layers()[l].nclasses);
        M += b.msz[l];
        K += b.ksz[l];
    }
    CHECK(b.M == M);
    CHECK(b.K == K);
}

TEST_CASE("the LQN warm start reaches the layer transients rather than being reset away") {
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = 20.0;
    ln::SolverLN<double> s(lqn_stage(0.8), lo);
    const ln::LnLayerBlocks b = s.layer_blocks();

    // Put a job on the SECOND station of every layer -- the server, not the
    // client Delay -- and read the trajectory back at t = 0. A fluid layer
    // starts exactly where it is told to, so this is an equality and not a
    // tendency; anything that drops the warm start reports the layer's default
    // marking instead, which puts the whole population on the reference station.
    Matrix<double> n(b.M, b.K, 0.0);
    std::size_t seeded = 0;
    for (std::size_t l = 0; l < b.msz.size(); ++l)
        if (b.msz[l] > 1 && b.ksz[l] > 0) {
            n(b.roff[l] + 1, b.coff[l]) = 2.0;
            ++seeded;
        }
    REQUIRE(seeded > 0);
    s.init_from_marginal(n);

    const ln::LnTranSolution tr = s.get_tran_avg();
    REQUIRE(tr.layers.size() == s.nlayers());
    std::size_t checked = 0;
    for (std::size_t l = 0; l < tr.layers.size(); ++l) {
        if (b.msz[l] <= 1 || b.ksz[l] == 0) continue;
        const ln::LnTranLayer& L = tr.layers[l];
        REQUIRE(!L.t.empty());
        CHECK(L.QN[1][0][0] == doctest::Approx(2.0).epsilon(1e-9));
        ++checked;
    }
    CHECK(checked == seeded);

    // And the marginal is cleared by an empty matrix, which puts the layers back
    // on their own default state rather than on a remembered one.
    s.init_from_marginal(Matrix<double>());
    const ln::LnTranSolution back = s.get_tran_avg();
    bool moved = false;
    for (std::size_t l = 0; l < back.layers.size(); ++l) {
        if (b.msz[l] <= 1 || b.ksz[l] == 0) continue;
        if (std::fabs(back.layers[l].QN[1][0][0] - 2.0) > 1e-6) moved = true;
    }
    CHECK(moved);
}

TEST_CASE("a layered stage is refused a marginal of the wrong shape") {
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = 20.0;
    ln::SolverLN<double> s(lqn_stage(0.8), lo);
    const ln::LnLayerBlocks b = s.layer_blocks();
    CHECK_THROWS_AS(s.init_from_marginal(Matrix<double>(b.M + 1, b.K, 0.0)), InputError);
    CHECK_THROWS_AS(s.init_from_marginal(Matrix<double>(b.M, b.K + 1, 0.0)), InputError);
}

TEST_CASE("renv_lqn_twostages: ENV over LQN stages runs, conserves and brackets") {
    const double T = 50.0;
    env::Environment<double> e = lqn_env(0.2, 1.0);
    const env::EnvAnalyzerSolution<double> r = env::solver_env(e, lqn_env_options(T, 20, 0.02));

    CHECK(r.method == "meanfield");
    CHECK(r.converged);

    // The aggregate has the shape of the layer blocks, not of the LQN elements.
    ln::LnOptions lo;
    lo.layer_solver = "fluid";
    lo.timespan_end = T;
    ln::SolverLN<double> probe(lqn_stage(0.8), lo);
    const ln::LnLayerBlocks b = probe.layer_blocks();
    CHECK(r.QN.rows() == b.M);
    CHECK(r.QN.cols() == b.K);

    for (std::size_t i = 0; i < r.QN.rows(); ++i)
        for (std::size_t j = 0; j < r.QN.cols(); ++j) REQUIRE(std::isfinite(r.QN(i, j)));

    double upQ = 0.0, upT = 0.0, dnQ = 0.0, dnT = 0.0;
    stage_alone(0.8, T, upQ, upT);
    stage_alone(3.0, T, dnQ, dnT);
    const double envQ = finite_sum(r.QN), envT = finite_sum(r.TN);

    // Population: every layer holds the same five client jobs whatever the
    // environment does, so the aggregate agrees with either stage taken alone.
    CHECK(std::fabs(envQ - upQ) < 1e-2);
    CHECK(std::fabs(envQ - dnQ) < 1e-2);

    // The bracket, the reference's own oracle.
    const double lo_x = std::min(upT, dnT), hi_x = std::max(upT, dnT);
    const double tol = 1e-2 * std::max(1.0, hi_x);
    CHECK(upT > dnT);  // a faster database does more work; the bracket is not degenerate
    CHECK(envT >= lo_x - tol);
    CHECK(envT <= hi_x + tol);
}

TEST_CASE("ENV over LQN is monotone in the stationary probability of the fast stage") {
    const double T = 50.0;
    double upQ = 0.0, upT = 0.0, dnQ = 0.0, dnT = 0.0;
    stage_alone(0.8, T, upQ, upT);
    stage_alone(3.0, T, dnQ, dnT);
    const double lo_x = std::min(upT, dnT), hi_x = std::max(upT, dnT);
    const double tol = 1e-2 * std::max(1.0, hi_x);

    const double x_low = env_tput(1.0, 0.2, T);   // P(UP) = 0.167
    const double x_high = env_tput(0.2, 1.0, T);  // P(UP) = 0.833
    CHECK(x_high > x_low + 1e-3);
    CHECK(x_low >= lo_x - tol);
    CHECK(x_high <= hi_x + tol);
}

TEST_CASE("a three-stage layered environment stays bracketed by the extreme stages") {
    const double T = 50.0;
    double upQ = 0.0, upT = 0.0, dnQ = 0.0, dnT = 0.0;
    stage_alone(0.8, T, upQ, upT);
    stage_alone(3.0, T, dnQ, dnT);
    const double lo_x = std::min(upT, dnT), hi_x = std::max(upT, dnT);
    const double tol = 1e-2 * std::max(1.0, hi_x);

    // MID's entry marginal mixes the exits of BOTH its neighbours by probOrig,
    // which is the part of the coupling a two-stage environment cannot exercise.
    env::Environment<double> e("R3", 3);
    e.set_lqn_stage(0, "UP", "operational", lqn_stage(0.8));
    e.set_lqn_stage(1, "MID", "degraded", lqn_stage(1.6));
    e.set_lqn_stage(2, "DOWN", "failed", lqn_stage(3.0));
    e.add_transition(0, 1, D::exp_rate(0.3));
    e.add_transition(1, 2, D::exp_rate(0.3));
    e.add_transition(2, 1, D::exp_rate(0.6));
    e.add_transition(1, 0, D::exp_rate(0.6));
    const double x3 = finite_sum(env::solver_env(e, lqn_env_options(T, 10, 0.03)).TN);
    CHECK(x3 >= lo_x - tol);
    CHECK(x3 <= hi_x + tol);
}

TEST_CASE("ENV refuses a layered stage where the coupling cannot take one, by name") {
    const double T = 50.0;

    // The state-vector coupling propagates a joint law across ONE stage
    // generator, which an LQN does not have. This is the reference's own
    // refusal (@@SolverENV/SolverENV.m), not a limit of this port.
    {
        env::Environment<double> e = lqn_env(0.2, 1.0);
        env::EnvOptions o = lqn_env_options(T, 5, 0.05);
        o.method = "statevec";
        o.stage_solver = "ctmc";
        CHECK_THROWS_AS(env::solver_env(e, o), UnsupportedError);
    }
    // The closed-form limits read a stage's station rate table; a layered model
    // has none, only the layers SolverLN derives from it.
    for (const std::string& m : {std::string("avg"), std::string("dec")}) {
        env::Environment<double> e = lqn_env(0.2, 1.0);
        env::EnvOptions o = lqn_env_options(T, 5, 0.05);
        o.method = m;
        CHECK_THROWS_AS(env::solver_env(e, o), UnsupportedError);
    }
    // A macro-state is one network at averaged rates, which a layered stage
    // cannot be aggregated into.
    {
        env::Environment<double> e = lqn_env(0.2, 1.0);
        env::EnvCompressOptions c;
        CHECK_THROWS_AS(env::solver_env(e, lqn_env_options(T, 5, 0.05), c), UnsupportedError);
    }
    // An enumerated CTMC stage solver has nothing to enumerate: the LQN has no
    // single generator, and `stage_solver` names the FLAT stage engine.
    {
        env::Environment<double> e = lqn_env(0.2, 1.0);
        env::EnvOptions o = lqn_env_options(T, 5, 0.05);
        o.stage_solver = "ctmc";
        o.stage_cutoff = 10.0;
        CHECK_THROWS_AS(env::solver_env(e, o), UnsupportedError);
    }
    // Only the fluid layer engine produces the stage transient the coupling
    // carries queue lengths across; an MVA layer has a fixed point and no
    // trajectory.
    {
        env::Environment<double> e = lqn_env(0.2, 1.0);
        env::EnvOptions o = lqn_env_options(T, 5, 0.05);
        o.lqn.layer_solver = "mva";
        CHECK_THROWS_AS(env::solver_env(e, o), UnsupportedError);
    }
    // Stages whose LAYER SHAPES differ cannot be blended entrywise, and the LQN
    // element count is not what decides it -- the layers SolverLN derives are.
    // The second stage is the first with one more served task behind the
    // database, so it is the same construct throughout and only the layering
    // differs.
    {
        env::Environment<double> e("Mixed", 2);
        lqn::LqnBuilder<double> b;
        b.processor("ClientProcessor", 1, SS::PS);
        b.processor("DBProcessor", 1, SS::PS);
        b.processor("DiskProcessor", 1, SS::PS);
        b.task("ClientTask", 5, SS::REF, "ClientProcessor");
        b.think_time("ClientTask", D::exp_mean(5.0));
        b.task("DBTask", std::numeric_limits<double>::infinity(), SS::INF, "DBProcessor");
        b.task("DiskTask", std::numeric_limits<double>::infinity(), SS::INF, "DiskProcessor");
        b.entry("ClientEntry", "ClientTask");
        b.entry("DBEntry", "DBTask");
        b.entry("DiskEntry", "DiskTask");
        b.activity("ClientActivity", D::exp_mean(1.0), "ClientTask");
        b.bound_to("ClientActivity", "ClientEntry");
        b.sync_call("ClientActivity", "DBEntry", 2.5);
        b.activity("DBActivity", D::exp_mean(0.8), "DBTask");
        b.bound_to("DBActivity", "DBEntry");
        b.sync_call("DBActivity", "DiskEntry", 1.0);
        b.replies_to("DBActivity", "DBEntry");
        b.activity("DiskActivity", D::exp_mean(0.5), "DiskTask");
        b.bound_to("DiskActivity", "DiskEntry");
        b.replies_to("DiskActivity", "DiskEntry");
        e.set_lqn_stage(0, "UP", "operational", lqn_stage(0.8));
        e.set_lqn_stage(1, "DOWN", "degraded", b.build());
        e.add_transition(0, 1, D::exp_rate(0.2));
        e.add_transition(1, 0, D::exp_rate(1.0));
        CHECK_THROWS_AS(env::solver_env(e, lqn_env_options(T, 5, 0.05)), InputError);
    }
}
