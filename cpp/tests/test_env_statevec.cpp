/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverENV, `method = "statevec"`: the full state-vector coupling.
 *
 * THE ORACLE EVERY TEST HERE IS BUILT ON, and the one the reference's own
 * comments name: an environment whose stages are IDENTICAL cannot change
 * anything the network does, so SolverENV must reproduce the plain SolverCTMC
 * solution of that one model. It needs no external reference, it is exact rather
 * than approximate, and it is what isolated the recorded Util drift documented
 * in `_kb/06-solver-catalog.md` ("ENV state-vector Util had drifted from the
 * CTMC analyzer"): QLen already matched to 2e-16 there while Util did not, which
 * placed the defect in the Util branch and nowhere else. So the identity is
 * asserted on Q, U AND T -- checking QLen alone is exactly the check that missed
 * the defect in three codebases at once.
 *
 * WHY IT IS EXACT AND NOT APPROXIMATE. The stationary distribution is a fixed
 * point of all three sojourn regimes: pi Q = 0 makes the resolvent
 * s pi (sI - Q)^{-1} = pi, makes the uniformized time-average over any interval
 * pi, and makes the transient pi(t) = pi for every t. `pre()` seeds from pi, so
 * the fixed point is reached at the first iteration and no quadrature error can
 * enter. That in turn lets each of the three regimes -- exponential, Erlang,
 * deterministic -- be exercised against the same oracle, which is the only way
 * to cover the transient branch without an external reference.
 *
 * The second family of checks is population conservation: no job is created or
 * destroyed at an environment switch, so a closed model's stage metrics and its
 * blend must all carry the same total, whatever the coupling does in between.
 *
 * NOT COVERED HERE: the cache hit/miss blend of `aggregateCacheBlend_`. The same
 * identity oracle applies to it and the test was written, but the CTMC cache
 * path it needs currently refuses the model outright ("the state space is empty;
 * no state satisfies the model's capacities") -- the same failure the existing
 * `test_ctmc_analyzer.cpp:285` shows on the identical network. The blend is
 * ported and unverified until that is resolved.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/env/solver_env_statevec.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

/** Closed: think-time Delay and an FCFS queue, `n` jobs, server rate `svc`. */
qn::Network<double> sv_stage(const std::string& nm, double svc, int n) {
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

/** The same closed cycle with a second queue, so it has three stations. */
qn::Network<double> sv_stage3(int n) {
    qn::Network<double> m("three");
    const std::size_t d = m.add_delay("ThinkTime");
    const std::size_t q1 = m.add_queue("Server1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Server2", lang::SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Jobs", n, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(2.0));
    m.set_service(q2, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, d, 1.0);
    m.link(P);
    return m;
}

/** Open: Source -> FCFS Queue -> Sink, one class, arrival `lam`, service `mu`. */
qn::Network<double> sv_mm1(const std::string& nm, double lam, double mu) {
    qn::Network<double> m(nm);
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(src, c, D::exp_rate(lam));
    m.set_service(q, c, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);
    return m;
}

/** A two-stage environment over two given models, with given transitions. */
env::Environment<double> sv_env(qn::Network<double>& a, qn::Network<double>& b, const D& ab,
                                const D& ba) {
    env::Environment<double> e("ServerModes", 2);
    e.set_stage(0, "A", "operational", a.get_struct());
    e.set_stage(1, "B", "degraded", b.get_struct());
    e.add_transition(0, 1, ab);
    e.add_transition(1, 0, ba);
    return e;
}

env::EnvStatevecOptions<double> sv_opt() {
    env::EnvStatevecOptions<double> o;
    o.iter_max = 50;
    o.iter_tol = 1e-10;
    // A few mean holding times, not a large safe number: the general branch
    // integrates with ode23 at a maximum step of one tenth of the span, and an
    // explicit pair started AT an equilibrium cannot detect its own instability
    // -- the local error estimate vanishes with the derivative, so an oversized
    // step is accepted and roundoff is amplified. See EnvStatevecOptions.
    o.timespan_end = 5.0;
    return o;
}

/** Every entry of the two metric matrices, at the given tolerance. */
void check_same(const Matrix<double>& got, const Matrix<double>& want, double tol) {
    REQUIRE(got.rows() == want.rows());
    REQUIRE(got.cols() == want.cols());
    for (std::size_t i = 0; i < got.rows(); ++i)
        for (std::size_t r = 0; r < got.cols(); ++r) {
            CAPTURE(i);
            CAPTURE(r);
            CHECK(got(i, r) == doctest::Approx(want(i, r)).epsilon(tol));
        }
}

}  // namespace

TEST_CASE("ENV statevec: identical stages reproduce SolverCTMC, all three sojourn regimes") {
    const int N = 3;
    qn::Network<double> m = sv_stage("Server", 2.0, N);
    ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> ref = ctmc::solver_ctmc_analyzer(m.get_struct(), copt);

    SUBCASE("exponential transitions take the exact resolvent, with no quadrature") {
        env::Environment<double> e = sv_env(m, m, D::exp_rate(0.5), D::exp_rate(1.0));
        const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, sv_opt());
        REQUIRE(s.converged);
        // The seed is already the fixed point, so it is reached immediately.
        CHECK(s.iterations == 1);
        check_same(s.QN, ref.avg.QN, 1e-12);
        check_same(s.UN, ref.avg.UN, 1e-12);
        check_same(s.TN, ref.avg.TN, 1e-12);
    }

    SUBCASE("an Erlang transition forces the general transient branch") {
        // Not exponential, so neither the resolvent nor the deterministic
        // shortcut applies and the answer goes through ctmc_transient and the
        // Riemann-Stieltjes sum against map_cdf. The identity must survive that,
        // and it does so to the last bit: the mean holding times are 0.5 and
        // 0.25 against a span of 5, so the CDFs have decayed well inside the
        // horizon AND ode23's maximum step stays inside its stability region.
        env::Environment<double> e = sv_env(m, m, D::erlang(4.0, 2), D::erlang(8.0, 2));
        const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, sv_opt());
        REQUIRE(s.converged);
        CHECK(s.iterations == 1);
        check_same(s.QN, ref.avg.QN, 1e-12);
        check_same(s.UN, ref.avg.UN, 1e-12);
        check_same(s.TN, ref.avg.TN, 1e-12);
    }

    SUBCASE("a deterministic sojourn goes through uniformization") {
        env::EnvStatevecOptions<double> o = sv_opt();
        o.sojourn = "deterministic";
        env::Environment<double> e = sv_env(m, m, D::exp_rate(0.5), D::exp_rate(1.0));
        const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, o);
        REQUIRE(s.converged);
        check_same(s.QN, ref.avg.QN, 1e-10);
        check_same(s.UN, ref.avg.UN, 1e-10);
        check_same(s.TN, ref.avg.TN, 1e-10);
    }
}

TEST_CASE("ENV statevec: distinct stages conserve population and stay between the stages") {
    // Fast and Slow differ only in the server rate, so each stage taken alone
    // has a SolverCTMC answer and the coupled one must sit between them: a stage
    // is entered from the other stage's exit and relaxes toward its own
    // stationary law without ever reaching it.
    const int N = 4;
    qn::Network<double> fast = sv_stage("Fast", 4.0, N);
    qn::Network<double> slow = sv_stage("Slow", 1.0, N);
    ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> rf = ctmc::solver_ctmc_analyzer(fast.get_struct(), copt);
    const ctmc::CtmcSolution<double> rs = ctmc::solver_ctmc_analyzer(slow.get_struct(), copt);

    env::Environment<double> e = sv_env(fast, slow, D::exp_rate(0.5), D::exp_rate(1.0));
    const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, sv_opt());
    REQUIRE(s.converged);

    // No job is created or lost at a switch, so every stage and the blend all
    // carry the same closed population exactly.
    for (std::size_t st = 0; st < 2; ++st) {
        CAPTURE(st);
        CHECK(s.QStage[st](0, 0) + s.QStage[st](1, 0) ==
              doctest::Approx(static_cast<double>(N)).epsilon(1e-9));
    }
    CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(static_cast<double>(N)).epsilon(1e-9));

    // The entry distributions are proper distributions over their state spaces.
    for (std::size_t st = 0; st < 2; ++st) {
        double tot = 0;
        for (double p : s.pi_enter[st]) tot += p;
        CAPTURE(st);
        CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));
    }

    // The queue is emptier in Fast than in Slow, and the coupled stages sit
    // strictly inside that interval -- the coupling is doing something, and it
    // is not overshooting either stage's own steady state.
    REQUIRE(rf.avg.QN(1, 0) < rs.avg.QN(1, 0));
    CHECK(s.QStage[0](1, 0) > rf.avg.QN(1, 0));
    CHECK(s.QStage[0](1, 0) < rs.avg.QN(1, 0));
    CHECK(s.QStage[1](1, 0) > rf.avg.QN(1, 0));
    CHECK(s.QStage[1](1, 0) < rs.avg.QN(1, 0));
    CHECK(s.QN(1, 0) > rf.avg.QN(1, 0));
    CHECK(s.QN(1, 0) < rs.avg.QN(1, 0));
}

TEST_CASE("ENV statevec: a state reset policy is applied at the switch") {
    // A reset that sends every switch into ONE state makes both entry
    // distributions that degenerate point mass, whatever the transient did.
    const int N = 3;
    qn::Network<double> m = sv_stage("Server", 2.0, N);
    env::Environment<double> e = sv_env(m, m, D::exp_rate(0.5), D::exp_rate(1.0));

    env::EnvStatevecOptions<double> o = sv_opt();
    const env::ResetStateVec<double> to_first = [](const std::vector<double>& p) {
        std::vector<double> q(p.size(), 0.0);
        q[0] = 1.0;
        return q;
    };
    o.reset_state.assign(2, std::vector<env::ResetStateVec<double>>(2));
    o.reset_state[0][1] = to_first;
    o.reset_state[1][0] = to_first;
    const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, o);
    REQUIRE(s.converged);
    for (std::size_t st = 0; st < 2; ++st) {
        CAPTURE(st);
        REQUIRE(!s.pi_enter[st].empty());
        CHECK(s.pi_enter[st][0] == doctest::Approx(1.0).epsilon(1e-12));
    }
    // The population is a property of the state space, so pinning the entry
    // state cannot break it.
    CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(static_cast<double>(N)).epsilon(1e-9));
    // And the answer is no longer the plain SolverCTMC one: the reset has
    // genuinely displaced the fixed point.
    ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> ref = ctmc::solver_ctmc_analyzer(m.get_struct(), copt);
    CHECK(s.QN(1, 0) != doctest::Approx(ref.avg.QN(1, 0)).epsilon(1e-4));
}

TEST_CASE("ENV statevec: refuses by name what it does not implement") {
    qn::Network<double> m = sv_stage("Server", 2.0, 2);
    auto build = [&]() { return sv_env(m, m, D::exp_rate(0.5), D::exp_rate(1.0)); };

    {   // The MAM backend accepts only what solver_mam_ldqbd accepts: one class
        // and two stations. A three-station stage is refused by that gate, not
        // by the coupling, which is why the message names the LDQBD method.
        qn::Network<double> m3 = sv_stage3(2);
        env::Environment<double> e = sv_env(m3, m3, D::exp_rate(0.5), D::exp_rate(1.0));
        env::EnvStatevecOptions<double> o = sv_opt();
        o.stage_solver = "mam";
        env::SolverEnvStatevec<double> s(e, o);
        CHECK_THROWS_AS(s.solve(), UnsupportedError);
    }
    {   // Anything else: the coupling needs an enumerated generator.
        env::Environment<double> e = build();
        env::EnvStatevecOptions<double> o = sv_opt();
        o.stage_solver = "fluid";
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), UnsupportedError);
    }
    {   // The horizon is REQUIRED, not defaulted: the default is infinite so
        // that a caller who never set one is refused, as the reference is.
        env::Environment<double> e = build();
        env::EnvStatevecOptions<double> o;  // timespan_end left at infinity
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), InputError);
    }
    {   // A per-stage horizon overrides the global one, and is checked too.
        env::Environment<double> e = build();
        env::EnvStatevecOptions<double> o = sv_opt();
        o.stage_timespan_end = {50.0, std::numeric_limits<double>::infinity()};
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), InputError);
        o.stage_timespan_end = {50.0, 0.0};
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), InputError);
    }
    {   // An unknown sojourn regime is an input error, not a silent default.
        env::Environment<double> e = build();
        env::EnvStatevecOptions<double> o = sv_opt();
        o.sojourn = "nonsense";
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), InputError);
    }
    {   // A reset that does not land in the destination's state space is named
        // together with the two sizes, as the reference's message is.
        env::Environment<double> e = build();
        env::EnvStatevecOptions<double> o = sv_opt();
        o.reset_state.assign(2, std::vector<env::ResetStateVec<double>>(2));
        o.reset_state[0][1] = [](const std::vector<double>& p) {
            return std::vector<double>(p.size() + 1, 0.0);
        };
        env::SolverEnvStatevec<double> s(e, o);
        CHECK_THROWS_AS(s.solve(), InputError);
    }
    {   // Stages that do not agree on stations and classes cannot be blended.
        qn::Network<double> wide = sv_stage3(2);
        env::Environment<double> e("mixed", 2);
        e.set_stage(0, "A", "operational", m.get_struct());
        e.set_stage(1, "B", "operational", wide.get_struct());
        e.add_transition(0, 1, D::exp_rate(0.5));
        e.add_transition(1, 0, D::exp_rate(1.0));
        env::EnvStatevecOptions<double> o = sv_opt();
        CHECK_THROWS_AS(env::SolverEnvStatevec<double>(e, o), InputError);
    }
}

TEST_CASE("ENV statevec: the MAM/LDQBD backend reproduces the plain LDQBD solve") {
    // THE SAME IDENTICAL-STAGE ORACLE as the CTMC backend, applied to the other
    // one: an environment whose two stages are the same model cannot change what
    // the network does, so the coupling must return `solver_mam_ldqbd`'s own
    // answer for that model. It needs no external reference and it is EXACT --
    // pi Q = 0 is a fixed point of all three sojourn regimes and `pre()` seeds
    // from it -- so any drift here is the flatten/reduce pair, which is the
    // only thing this backend adds.
    qn::Network<double> m = sv_stage("Server", 2.0, 3);
    env::Environment<double> e = sv_env(m, m, D::exp_rate(0.5), D::exp_rate(1.0));
    env::EnvStatevecOptions<double> o = sv_opt();
    o.stage_solver = "mam";
    const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, o);
    REQUIRE(s.converged);
    CHECK(s.iterations == 1);

    mam::MamOptions mo;
    mo.method = "default";
    const mam::LdqbdSolution<double> ref = mam::solver_mam_ldqbd(m.get_struct(), mo);
    REQUIRE(s.QN.rows() == ref.sol.Q.rows());
    for (std::size_t i = 0; i < s.QN.rows(); ++i) {
        CAPTURE(i);
        CHECK(s.QN(i, 0) == doctest::Approx(ref.sol.Q(i, 0)).epsilon(1e-9));
        CHECK(s.UN(i, 0) == doctest::Approx(ref.sol.U(i, 0)).epsilon(1e-9));
        CHECK(s.TN(i, 0) == doctest::Approx(ref.sol.Tp(i, 0)).epsilon(1e-9));
    }

    // Population conservation, the second family of checks in this file: the
    // closed model holds three jobs and the blend must still hold three.
    double total = 0.0;
    for (std::size_t i = 0; i < s.QN.rows(); ++i) total += s.QN(i, 0);
    CHECK(total == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("the LDQBD flatten/reduce pair is self-consistent") {
    qn::Network<double> m = sv_stage("Server", 2.0, 3);
    mam::MamOptions mo;
    mo.method = "default";
    const mam::LdqbdSolution<double> sol = mam::solver_mam_ldqbd(m.get_struct(), mo);
    const mam::LdqbdFlat<double> flat = mam::solver_mam_ldqbd_flatten(sol.ld);

    // The flat matrix is a generator: rows sum to zero, off-diagonals are
    // non-negative. A block placed at the wrong offset breaks this immediately,
    // which is why it is checked before any number is read off it.
    REQUIRE(flat.Q.rows() == flat.levelOf.size());
    for (std::size_t i = 0; i < flat.Q.rows(); ++i) {
        double rs = 0.0;
        for (std::size_t j = 0; j < flat.Q.cols(); ++j) {
            if (i != j) CHECK(flat.Q(i, j) >= -1e-14);
            rs += flat.Q(i, j);
        }
        CAPTURE(i);
        CHECK(std::fabs(rs) <= 1e-12);
    }
    // Level 0 is ONE state and the rest are nPhases wide, so the level map is
    // not a multiple of the index; it must still be non-decreasing and reach Nlev.
    CHECK(flat.levelOf[0] == 0u);
    CHECK(flat.levelOf.back() == sol.ld.Nlev);
    for (std::size_t i = 1; i < flat.levelOf.size(); ++i)
        CHECK(flat.levelOf[i] >= flat.levelOf[i - 1]);

    // Reducing the flat chain's OWN stationary law must give back what the
    // block recursion reported, which is the check that ties the two halves.
    const std::vector<double> pi = mc::ctmc_solve_reducible(flat.Q).pi;
    const mam::LdqbdAvg<double> a = mam::solver_mam_ldqbd_avg(sol.ld, pi, flat.levelOf);
    for (std::size_t i = 0; i < a.QN.rows(); ++i) {
        CAPTURE(i);
        CHECK(a.QN(i, 0) == doctest::Approx(sol.sol.Q(i, 0)).epsilon(1e-9));
        CHECK(a.UN(i, 0) == doctest::Approx(sol.sol.U(i, 0)).epsilon(1e-9));
        CHECK(a.TN(i, 0) == doctest::Approx(sol.sol.Tp(i, 0)).epsilon(1e-9));
    }
}

/*
 * REGRESSION, `_kb/06-solver-catalog.md` "ENV statevec was warm-started from a
 * law that does not exist".
 *
 * The seed used to be each stage's OWN stationary law. A stage that is
 * individually critical or unstable has none, so the reducible solve returned
 * the stationary law of the TRUNCATED generator -- mass against the truncation
 * wall, mean N/2 -- and the sweeps needed to drain it grew with the cutoff. The
 * fixed point was always correct; at a finite `iter_max` the REPORTED number
 * drifted away from it as the cutoff was RAISED.
 *
 * Here the Slow stage is critical (lambda = mu = 0.8) while the system is stable
 * on average (mu_bar = 1.4), which is the regime that exposed it. Raising the
 * cutoff must not move the answer: the truncation is far past the mass either
 * way, so anything that moves is the seed, not the model. Before the fix the
 * same sweep budget gave 1.478346 / 1.494480 / 3.771534 / 28.608647 at cutoffs
 * 40 / 80 / 150 / 300 against an exact 1.478318.
 *
 * THE SWEEP BUDGET IS PART OF THE TEST. The fixed seed reaches the fixed point
 * in ~220 sweeps at EVERY cutoff; the old per-stage seed needs ~350 at cutoff 40
 * but ~900 at cutoff 300, which is why the LARGE cutoff is what discriminates
 * here. Raising `iter_max` far enough lets the old seed converge too and
 * SILENTLY DEFANGS this test -- it would still pass while testing nothing.
 */
TEST_CASE("ENV statevec: an individually critical stage is cutoff-invariant") {
    // Exact joint (queue,stage) CTMC answer for this environment.
    const double kExactQ = 1.478318;
    const double kLambda = 0.8;

    qn::Network<double> slow = sv_mm1("slow", kLambda, 0.8);  // critical on its own
    qn::Network<double> fast = sv_mm1("fast", kLambda, 2.0);

    std::vector<double> qlen;
    for (double cutoff : {40.0, 300.0}) {
        env::Environment<double> e = sv_env(slow, fast, D::exp_rate(1.0), D::exp_rate(1.0));
        env::EnvStatevecOptions<double> o = sv_opt();
        o.iter_max = 400;
        o.iter_tol = 1e-10;
        o.stage.cutoff = cutoff;
        const env::EnvStatevecSolution<double> s = env::solver_env_statevec(e, o);
        CAPTURE(cutoff);
        // A run that stopped on iter_max says nothing about the fixed point, so
        // assert convergence FIRST -- that is the half that used to be silent.
        CHECK(s.converged);

        double Q = 0.0, T = 0.0;
        for (std::size_t i = 0; i < s.QN.rows(); ++i)
            for (std::size_t r = 0; r < s.QN.cols(); ++r) Q += s.QN(i, r);
        // Flow balance at the queue is the cheap tell the defect violated by up
        // to 33%: in steady state the queue must clear exactly what arrives.
        for (std::size_t r = 0; r < s.TN.cols(); ++r) T += s.TN(1, r);
        CHECK(T == doctest::Approx(kLambda).epsilon(1e-6));
        CHECK(Q == doctest::Approx(kExactQ).epsilon(1e-4));
        qlen.push_back(Q);
    }
    // The invariance itself, independent of how close either is to the truth.
    CHECK(qlen[0] == doctest::Approx(qlen[1]).epsilon(1e-6));
}
