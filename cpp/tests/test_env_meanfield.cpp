/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverENV, the default mean-field path, and the environment COMPRESSION that
 * `@SolverENV/SolverENV.m` wraps around it.
 *
 * THE ORACLE EVERY TEST HERE IS BUILT ON is the same one that governs
 * `test_env_statevec.cpp`: an environment whose stages are IDENTICAL cannot
 * change anything the network does, so the coupled answer must be the plain
 * single-model answer of that one stage. For a fluid stage the plain answer is
 * its own steady state, and the identity is exact for a reason worth stating:
 * the steady state is a stationary point of the fluid drift, so a transient
 * started there does not move, and a CDF-weighted average of a constant is that
 * constant. It holds only because every service here is EXPONENTIAL -- the
 * handoff reconstructs a stage's state from marginal means by putting all mass
 * in phase one, which is lossless at one phase per class and lossy beyond it.
 *
 * THE SECOND ORACLE IS SPECIFIC TO COMPRESSION, and it is sharper than a
 * tolerance. With SINGLETON macro-states every block is 1 x 1, so the
 * conditional distribution inside a block is identically one, the aggregated
 * macro chain is the uniformized chain itself, and the whole construction has
 * to be the identity: `p` is then the EXACT stationary vector rather than an
 * approximation of it, `macro_rate` is `E0` entry for entry, and the compressed
 * solve must reproduce the uncompressed one. Any indexing slip in the
 * permutation, in pmicro, or in the macro-rate double sum breaks that identity
 * immediately, which is what makes it the right wiring test for four kernels
 * that are individually already tested.
 *
 * Two things this file does NOT assert. It never compares the mean-field answer
 * against the state-vector one at a tight tolerance: mean-field carries only
 * marginal means across a switch and is an approximation exactly where the joint
 * distribution matters. And it never asserts that a particular pair of stages
 * gets merged: which merge minimizes eps is the search's business, so what is
 * checked is that the search returns a valid partition no worse than the
 * singletons it started from.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/env/solver_env_meanfield.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

/** Closed: think-time Delay and an FCFS queue, `n` jobs, server rate `svc`. */
qn::Network<double> mf_stage(const std::string& nm, double svc, int n) {
    qn::Network<double> m(nm);
    const std::size_t d = m.add_delay("ThinkTime");
    const std::size_t q = m.add_queue("Server", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Jobs", n, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(svc));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

env::EnvOptions mf_opt() {
    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    // The NESTED integrator's tolerance, not the fixed point's. FluidOptions
    // defaults it to 1e-4, and the identical-stages case below compares against
    // the plain fluid solution at 1e-6: an assertion ten times tighter than the
    // integrator it is reading. Measured relative error against that reference
    // is 2.1e-6 at the 1e-4 default and 2.8e-9 from 1e-8 onwards, where it
    // plateaus -- so this buys three orders of margin and nothing beyond it.
    // Tightening iter_tol instead does NOTHING, because the fixed point is not
    // what carries the error.
    o.stage.tol = 1e-8;
    o.timespan_end = 100.0;
    o.tran_points = 2001;
    return o;
}

/** Two stages over the given models, Exp(0.5) forward and Exp(1.0) back. */
env::Environment<double> mf_env2(qn::Network<double>& a, qn::Network<double>& b) {
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "A", "operational", a.get_struct());
    e.set_stage(1, "B", "degraded", b.get_struct());
    e.add_transition(0, 1, D::exp_rate(0.5));
    e.add_transition(1, 0, D::exp_rate(1.0));
    return e;
}

/**
 * Four stages in two tightly coupled pairs, {0,1} and {2,3}, joined by rare
 * transitions: the textbook nearly completely decomposable environment, and the
 * shape the whole compression path exists for.
 */
env::Environment<double> mf_env_ncd(qn::Network<double>& m, double fast, double slow) {
    env::Environment<double> e("NCD", 4);
    const char* nms[4] = {"A1", "A2", "B1", "B2"};
    for (std::size_t i = 0; i < 4; ++i) e.set_stage(i, nms[i], "operational", m.get_struct());
    e.add_transition(0, 1, D::exp_rate(fast));
    e.add_transition(1, 0, D::exp_rate(fast));
    e.add_transition(2, 3, D::exp_rate(fast));
    e.add_transition(3, 2, D::exp_rate(fast));
    e.add_transition(1, 2, D::exp_rate(slow));
    e.add_transition(3, 0, D::exp_rate(slow));
    return e;
}

env::MacroPartition mf_singletons(std::size_t E) {
    env::MacroPartition MS(E);
    for (std::size_t i = 0; i < E; ++i) MS[i] = std::vector<std::size_t>{i};
    return MS;
}

}  // namespace

TEST_CASE("ENV meanfield: identical stages reproduce the plain fluid solution") {
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env2(m, m);
    const env::EnvMeanfieldSolution<double> s = env::solver_env_meanfield(e, mf_opt());
    REQUIRE(s.avg.converged);
    CHECK_FALSE(s.compressed);

    fluid::FluidOptions fo;
    const fluid::FluidSolution ref = fluid::solver_fluid(m.get_struct(), fo);
    for (std::size_t i = 0; i < 2; ++i) {
        CAPTURE(i);
        CHECK(s.avg.QN(i, 0) == doctest::Approx(ref.QN(i, 0)).epsilon(1e-6));
    }
    // No job is created or lost at a switch, whatever the coupling does.
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("ENV meanfield: singleton macro-states make the decomposition exact") {
    // Every block is 1 x 1, so the aggregated chain IS the uniformized chain and
    // the four kernels must all return the environment's exact stationary law,
    // not an approximation of it. The environment is Exp(0.5) one way and
    // Exp(1.0) the other, so it sits in stage one two thirds of the time.
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env2(m, m);
    const Matrix<double> E0 = env::env_rate_matrix(e);
    CHECK(E0(0, 1) == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(E0(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(E0(0, 0) == doctest::Approx(0.0).epsilon(1e-12));

    const Matrix<double> Eutil = mc::ctmc_makeinfgen(E0);
    const env::MacroPartition MS = mf_singletons(2);
    const char* kernels[] = {"courtois", "kms", "takahashi", "multi"};
    for (const char* k : kernels) {
        CAPTURE(k);
        env::EnvCompressOptions c;
        c.da = k;
        const env::EnvDecomp<double> d = env::env_ctmc_decompose(Eutil, MS, c);
        REQUIRE(d.p.size() == 2);
        CHECK(d.p[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
        CHECK(d.p[1] == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
        CHECK(d.p[0] + d.p[1] == doctest::Approx(1.0).epsilon(1e-12));
        // The uniformization rate is 1.05 max|Q|, and max|Q| here is the 1.0 of
        // the second row's diagonal.
        CHECK(d.q == doctest::Approx(1.05).epsilon(1e-9));
    }
}

TEST_CASE("ENV meanfield: singleton compression is the identity, end to end") {
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env2(m, m);

    env::EnvCompressOptions c;
    c.partition = mf_singletons(2);  // skip the search: this is the wiring test
    const env::EnvCompression<double> comp = env::env_compress(e, c);

    REQUIRE(comp.MS.size() == 2);
    CHECK(comp.MS[0].size() == 1);
    CHECK(comp.MS[1].size() == 1);
    REQUIRE(comp.env);
    CHECK(comp.env->nstages() == 2);

    // pmicro is identically one inside a singleton, so the aggregated rates are
    // the original ones and nothing has been averaged away.
    CHECK(comp.pmicro[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(comp.pmicro[1] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(comp.macro_rate(0, 1) == doctest::Approx(comp.E0(0, 1)).epsilon(1e-12));
    CHECK(comp.macro_rate(1, 0) == doctest::Approx(comp.E0(1, 0)).epsilon(1e-12));
    CHECK(comp.pmacro[0] == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    CHECK(comp.pmacro[1] == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(comp.pmacro[0] + comp.pmacro[1] == doctest::Approx(1.0).epsilon(1e-12));
    // A two-stage alternating environment always comes from the other stage.
    CHECK(comp.prob_orig(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(comp.prob_orig(0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(comp.prob_orig(0, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // And the whole solve is unchanged, which is the property that catches an
    // indexing slip anywhere between the permutation and the macro networks.
    env::Environment<double> plain = mf_env2(m, m);
    const env::EnvMeanfieldSolution<double> a = env::solver_env_meanfield(plain, mf_opt());
    env::Environment<double> again = mf_env2(m, m);
    const env::EnvMeanfieldSolution<double> b = env::solver_env_meanfield(again, mf_opt(), c);
    REQUIRE(b.compressed);
    for (std::size_t i = 0; i < 2; ++i) {
        CAPTURE(i);
        CHECK(b.avg.QN(i, 0) == doctest::Approx(a.avg.QN(i, 0)).epsilon(1e-9));
        CHECK(b.avg.UN(i, 0) == doctest::Approx(a.avg.UN(i, 0)).epsilon(1e-9));
    }
    CHECK(b.avg.QN(0, 0) + b.avg.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("ENV meanfield: the macro probabilities agree with the environment's own") {
    // env_apply_macro_probabilities overwrites probEnv with the decomposition's
    // block sums, and Environment::init derives probEnv from the macro arcs.
    // For singletons both are exact, so the overwrite must be a no-op -- which
    // is the check that the two halves of the compression agree rather than
    // merely coexisting.
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env2(m, m);
    env::EnvCompressOptions c;
    c.partition = mf_singletons(2);
    const env::EnvCompression<double> comp = env::env_compress(e, c);

    comp.env->init();
    const std::vector<double> from_init = comp.env->prob_env;
    REQUIRE(from_init.size() == 2);
    env::env_apply_macro_probabilities(*comp.env, comp);
    REQUIRE(comp.env->prob_env.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        CAPTURE(i);
        CHECK(comp.env->prob_env[i] == doctest::Approx(from_init[i]).epsilon(1e-9));
    }
}

TEST_CASE("ENV meanfield: a nearly decomposable environment compresses") {
    // Two fast pairs joined by rare transitions. The point of the model is that
    // the coupling eps is small against the admissible epsMAX, which is the
    // condition under which aggregating the pairs is meaningful at all.
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env_ncd(m, 100.0, 0.1);
    const Matrix<double> Eutil = mc::ctmc_makeinfgen(env::env_rate_matrix(e));

    env::EnvCompressOptions c;
    c.partition = env::MacroPartition{{0, 1}, {2, 3}};
    const env::EnvDecomp<double> d = env::env_ctmc_decompose(Eutil, c.partition, c);

    // The aggregation is only meaningful when eps stays under epsMAX, and on a
    // model built to be NCD it must.
    CHECK(d.eps < d.epsMAX);
    double tot = 0;
    for (double v : d.p) tot += v;
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));

    // Against the exact stationary law of the same generator: Courtois is an
    // approximation, and on a genuinely NCD chain it is a close one.
    const std::vector<double> exact = mc::ctmc_solve(Eutil);
    REQUIRE(exact.size() == d.p.size());
    double l1 = 0;
    for (std::size_t i = 0; i < exact.size(); ++i) l1 += std::fabs(exact[i] - d.p[i]);
    CHECK(l1 < 1e-2);

    // The search must never return something worse than where it started, and
    // must always return a partition of all four stages.
    const env::MacroPartition found = env::env_find_best_partition(Eutil, c);
    std::vector<bool> seen(4, false);
    std::size_t covered = 0;
    for (const std::vector<std::size_t>& blk : found) {
        CHECK(!blk.empty());
        for (std::size_t s : blk) {
            REQUIRE(s < 4);
            CHECK_FALSE(seen[s]);
            seen[s] = true;
            ++covered;
        }
    }
    CHECK(covered == 4);
    const env::EnvDecomp<double> dfound = env::env_ctmc_decompose(Eutil, found, c);
    const env::EnvDecomp<double> dsingle = env::env_ctmc_decompose(Eutil, mf_singletons(4), c);
    CHECK(dfound.eps <= dsingle.eps);

    // The beam search is a different search, not a budgeted one, so it is only
    // required to return a valid partition of the same stages.
    env::EnvCompressOptions cb = c;
    cb.partition.clear();
    const env::MacroPartition beam = env::env_beam_search_partition(Eutil, cb);
    std::size_t bcov = 0;
    for (const std::vector<std::size_t>& blk : beam) bcov += blk.size();
    CHECK(bcov == 4);
}

TEST_CASE("ENV meanfield: compressing an NCD environment conserves the population") {
    // The compressed model is a different model, so its queue lengths may move;
    // what cannot move is the closed population, which is a property of the
    // network rather than of the environment.
    qn::Network<double> m = mf_stage("Server", 2.0, 5);
    env::Environment<double> e = mf_env_ncd(m, 100.0, 0.1);
    env::EnvCompressOptions c;
    c.partition = env::MacroPartition{{0, 1}, {2, 3}};
    const env::EnvMeanfieldSolution<double> s = env::solver_env_meanfield(e, mf_opt(), c);
    REQUIRE(s.compressed);
    CHECK(s.compression.MS.size() == 2);
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
    // Identical stages throughout, so the aggregated rates are the originals and
    // the compressed answer is still the plain fluid one.
    fluid::FluidOptions fo;
    const fluid::FluidSolution ref = fluid::solver_fluid(m.get_struct(), fo);
    CHECK(s.avg.QN(1, 0) == doctest::Approx(ref.QN(1, 0)).epsilon(1e-6));
}

TEST_CASE("ENV meanfield: refuses by name what it does not implement") {
    qn::Network<double> m = mf_stage("Server", 2.0, 3);

    {   // The NCD construction reads the environment as a CTMC through
        // E0 = getRate(), so a non-exponential arc cannot be aggregated without
        // silently replacing it by an exponential of the same mean.
        env::Environment<double> e("ph", 2);
        e.set_stage(0, "A", "op", m.get_struct());
        e.set_stage(1, "B", "op", m.get_struct());
        e.add_transition(0, 1, D::erlang(2.0, 2));
        e.add_transition(1, 0, D::exp_rate(1.0));
        env::EnvCompressOptions c;
        c.partition = mf_singletons(2);
        CHECK_THROWS_AS(env::env_compress(e, c), UnsupportedError);
    }
    {   // options.config.da names one of four kernels and nothing else.
        env::Environment<double> e = mf_env2(m, m);
        const Matrix<double> Eutil = mc::ctmc_makeinfgen(env::env_rate_matrix(e));
        env::EnvCompressOptions c;
        c.da = "nonsense";
        CHECK_THROWS_AS(env::env_ctmc_decompose(Eutil, mf_singletons(2), c), UnsupportedError);
    }
    {   // Every kernel assumes MS partitions the stages and none of them checks.
        env::Environment<double> e = mf_env2(m, m);
        const Matrix<double> Eutil = mc::ctmc_makeinfgen(env::env_rate_matrix(e));
        env::EnvCompressOptions c;
        CHECK_THROWS_AS(  // a stage in two blocks
            env::env_ctmc_decompose(Eutil, env::MacroPartition{{0, 1}, {1}}, c), InputError);
        CHECK_THROWS_AS(  // a stage in none
            env::env_ctmc_decompose(Eutil, env::MacroPartition{{0}}, c), InputError);
        CHECK_THROWS_AS(  // a stage that does not exist
            env::env_ctmc_decompose(Eutil, env::MacroPartition{{0}, {1}, {2}}, c), InputError);
        CHECK_THROWS_AS(  // an empty block
            env::env_ctmc_decompose(Eutil, env::MacroPartition{{0, 1}, {}}, c), InputError);
    }
    {   // Merging every stage into one block leaves nothing for the environment
        // to switch to, and an absorbing environment has no stage probabilities.
        env::Environment<double> e = mf_env2(m, m);
        env::EnvCompressOptions c;
        c.partition = env::MacroPartition{{0, 1}};
        CHECK_THROWS_AS(env::env_compress(e, c), InputError);
    }
    {   // A Cache under a NON-FLUID stage solver: only the fluid stage exposes
        // the RMF transient the blend integrates, so the other backends are
        // refused by name rather than reported without hit and miss ratios.
        qn::Network<double> cm("cacheqn");
        const std::size_t src = cm.add_source("Src");
        const std::size_t sk = cm.add_sink("Sink");
        const std::size_t rd = cm.add_open_class("Read");
        const std::size_t hit = cm.add_open_class("Hit");
        const std::size_t mis = cm.add_open_class("Miss");
        qn::CacheParam<double> cp;
        cp.nitems = 2;
        cp.itemcap.push_back(1);
        cp.replacestrat = lang::ReplacementStrategy::LRU;
        cp.pread.assign(3, std::vector<double>());
        cp.pread[rd - 1] = std::vector<double>{0.5, 0.5};
        cp.hitclass.assign(3, 0);
        cp.missclass.assign(3, 0);
        cp.hitclass[rd - 1] = hit;
        cp.missclass[rd - 1] = mis;
        const std::size_t ca = cm.add_cache("C", cp);
        cm.set_arrival(src, rd, D::exp_rate(1.0));
        qn::RoutingMatrix<double> P;
        P.set(rd, rd, src, ca, 1.0);
        P.set(hit, hit, ca, sk, 1.0);
        P.set(mis, mis, ca, sk, 1.0);
        cm.link(P);
        env::Environment<double> e = mf_env2(cm, cm);
        env::EnvOptions o = mf_opt();
        o.stage_solver = "ctmc";
        CHECK_THROWS_AS(env::solver_env_meanfield(e, o), UnsupportedError);
    }
}

TEST_CASE("the mean-field cache blend runs and reports a probability") {
    // Source -> Cache -{HitQueue, MissQueue}- Sink, in a two-stage environment
    // whose stages differ only in the miss-path service rate, so the cache sees
    // the same read stream in both and the blend has something to average.
    auto cache_model = [](double miss_rate) {
        qn::Network<double> m("cacheqn");
        const std::size_t src = m.add_source("Src");
        qn::CacheParam<double> cp;
        cp.nitems = 4;
        cp.itemcap.push_back(2);
        cp.replacestrat = lang::ReplacementStrategy::RR;
        cp.pread.assign(3, std::vector<double>());
        cp.pread[0] = std::vector<double>(4, 0.25);
        cp.hitclass.assign(3, 0);
        cp.missclass.assign(3, 0);
        cp.hitclass[0] = 2;
        cp.missclass[0] = 3;
        const std::size_t ca = m.add_cache("C", cp);
        const std::size_t hq = m.add_queue("HitQueue", lang::SchedStrategy::PS);
        const std::size_t mq = m.add_queue("MissQueue", lang::SchedStrategy::PS);
        const std::size_t sk = m.add_sink("Sink");
        const std::size_t rd = m.add_open_class("Read");
        const std::size_t hit = m.add_open_class("Hit");
        const std::size_t mis = m.add_open_class("Miss");
        m.set_arrival(src, rd, D::exp_rate(1.0));
        m.set_service(hq, hit, D::exp_rate(4.0));
        m.set_service(mq, mis, D::exp_rate(miss_rate));
        qn::RoutingMatrix<double> P;
        P.set(rd, rd, src, ca, 1.0);
        P.set(hit, hit, ca, hq, 1.0);
        P.set(hit, hit, hq, sk, 1.0);
        P.set(mis, mis, ca, mq, 1.0);
        P.set(mis, mis, mq, sk, 1.0);
        m.link(P);
        return m;
    };

    // A SHORT SWEEP BUDGET ON PURPOSE. Every sweep of the cache fixed point is
    // a full solver_fld_cacheqn_tran per stage -- a fluid network solve plus an
    // adaptive Rosenbrock on the cache drift -- so the 50 sweeps mf_opt() gives
    // the queue-length fixed point would cost minutes here for a quantity that
    // stops moving after a handful. The occupancy handoff converges much faster
    // than the queue-length one because a cache's exit occupancy depends on its
    // entry occupancy only through the sojourn.
    env::EnvOptions co = mf_opt();
    co.iter_max = 3;
    qn::Network<double> a = cache_model(2.0), b = cache_model(0.5);
    env::Environment<double> e = mf_env2(a, b);
    const env::EnvMeanfieldSolution<double> s = env::solver_env_meanfield(e, co);

    REQUIRE(s.cache.nodes.size() == 1u);
    const std::vector<double>& hp = s.cache.hitprob[0];
    const std::vector<double>& mp = s.cache.missprob[0];
    // The READ class is the only one that accesses the cache, so it is the only
    // one carrying a ratio; the switched classes never reach the cache and stay
    // NaN rather than being reported as zero, which would read as "never hits".
    REQUIRE(hp.size() >= 1u);
    CHECK(std::isfinite(hp[0]));
    CHECK(hp[0] >= 0.0);
    CHECK(hp[0] <= 1.0);
    CHECK(hp[0] + mp[0] == doctest::Approx(1.0).epsilon(1e-12));
    // Two items held of four read uniformly: the hit ratio cannot exceed the
    // share of the catalogue the cache can hold under any replacement policy.
    CHECK(hp[0] <= 0.5 + 1e-9);
    for (std::size_t k = 1; k < hp.size(); ++k) {
        CAPTURE(k);
        CHECK(std::isnan(hp[k]));
    }

    // The blend must sit between the two stages' own hit ratios: it is a
    // weighted average of them, so a value outside that band means the
    // per-request-to-real-time rescaling or the prob_env weighting is wrong.
    double lo = 1.0, hi = 0.0;
    for (int which = 0; which < 2; ++which) {
        qn::Network<double> one = cache_model(which == 0 ? 2.0 : 0.5);
        env::Environment<double> se = mf_env2(one, one);
        const env::EnvMeanfieldSolution<double> ss = env::solver_env_meanfield(se, co);
        REQUIRE(ss.cache.nodes.size() == 1u);
        const double h = ss.cache.hitprob[0][0];
        lo = std::min(lo, h);
        hi = std::max(hi, h);
    }
    CHECK(hp[0] >= lo - 1e-6);
    CHECK(hp[0] <= hi + 1e-6);
}
