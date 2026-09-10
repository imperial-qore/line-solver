/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The per-method support gates of SolverCTMC and SolverFLD.
 *
 * WHAT THIS PINS. `model.help()` / `findSolver()` reports one row per (solver,
 * method) pair and asks the per-method gate -- the same one the AUTO ranking
 * consults before delegating. Until the rules below reached it, that gate was
 * much weaker than what the ANALYZERS enforce at run time, so the report
 * offered pairs that then threw:
 *
 *   ctmc cftp / cftp.approx   single-class and closed-only
 *   ctmc mdd                  single-class and closed-only, and no fork-join
 *   ctmc and ssa, every one   the fork-join model class sn_fj_validate admits
 *   fluid diffusion           closed-only, and no fork-join
 *   fluid kp                  no fork-join
 *   fluid tbi                 closed-only, and no cache
 *   fluid refined             closed-only
 *   fluid dae                 no OPEN fork-join model
 *   fluid mol / mtginf        a finite horizon
 *
 * TWO ROUTES REACH THE REPORT, and both are asserted here. What the feature
 * registry can NAME rides in `qn::ctmc_feature_set(method)` and
 * `qn::fluid_feature_set(method)`; what it cannot -- a class count, a station
 * count, a horizon -- is a structural predicate that the ANALYZER ITSELF calls,
 * so the gate and the run are one body of rules rather than two copies.
 *
 * Each case asserts the refusal AND its converse: a model the method IS derived
 * for must keep it. Over-tightening a gate hides a method the user could have
 * run, which is the same defect with the sign flipped.
 *
 * The MATLAB twin is @SolverCTMC/supportsModelMethod.m and
 * @SolverFLD/supportsModelMethod.m; the Python twin is
 * python/tests/test_gate_ctmc_fluid.py and the JAR twin
 * jar/src/test/java/jline/solvers/ctmc/CtmcFluidGateTest.java.
 */
#include <limits>
#include <string>

#include "doctest.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cftp.h"
#include "line/solvers/ctmc/solver_ctmc_mdd_analyzer.h"
#include "line/solvers/fluid/fluid_qsys.h"
#include "line/solvers/fluid/fluid_runner.h"

using line::UnsupportedError;
using line::qn::Feature;
using line::qn::feature_set_supports;
using line::qn::FeatureSet;
using line::qn::NetworkStruct;
using line::qn::used_lang_features;
using line::lang::SchedStrategy;

namespace {

using NetD = line::qn::Network<double>;
using Dist = line::lang::Distrib<double>;
using SN = NetworkStruct<double>;

/** Source -> Queue -> Sink, one open class: the smallest open model. */
NetD mm1() {
    NetD m("mm1");
    const std::size_t s = m.add_source("S");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("K");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(s, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Queue -> Delay, ONE closed class: what cftp and mdd are for. */
NetD repairmen() {
    NetD m("rep");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The same shape with TWO closed classes: closed, but not single-class. */
NetD cqn2() {
    NetD m("cqn2");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(2.0));
    m.set_service(q, c1, Dist::exp_rate(2.0));
    m.set_service(q, c2, Dist::exp_rate(3.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Fork -> two FCFS queues -> Join.
 *
 * `paired` names the Fork on the Join, which is what DECLARES the fork-join
 * pairing: `sn.fj` is read off that declaration and not derived from the
 * routing, because a nested model (fj_basic_nesting) has two forks and two
 * joins the routing alone does not pair. A Join built without it leaves
 * `sn.fj` empty, which the exact construction refuses.
 */
NetD forkjoin(bool paired, bool closed) {
    NetD m("fj");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t jn = paired ? m.add_join("J", f) : m.add_join_unbound("J");
    std::size_t entry = 0, exit_node = 0, c = 0;
    if (closed) {
        const std::size_t d = m.add_delay("D");
        c = m.add_closed_class("C", 2.0, d);
        m.set_service(d, c, Dist::exp_rate(1.0));
        entry = d;
        exit_node = d;
    } else {
        const std::size_t s = m.add_source("S");
        const std::size_t k = m.add_sink("K");
        c = m.add_open_class("C");
        m.set_arrival(s, c, Dist::exp_rate(0.5));
        entry = s;
        exit_node = k;
    }
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c, c, entry, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, jn, 1.0);
    P.set(c, c, q2, jn, 1.0);
    P.set(c, c, jn, exit_node, 1.0);
    m.link(P);
    return m;
}

/** Does `declared` accept everything `sn` uses? */
bool accepts(const FeatureSet& declared, const SN& sn) {
    return feature_set_supports("gate", declared, used_lang_features(sn)).ok;
}

}  // namespace

// ---------------------------------------------------------------------------
// SolverCTMC: what the registry CAN name
// ---------------------------------------------------------------------------

TEST_CASE("ctmc_feature_set: cftp and mdd drop the open-model names") {
    const FeatureSet base = line::qn::ctmc_feature_set("default");
    CHECK(base.has(Feature::OpenClass));
    CHECK(base.has(Feature::Source));

    const std::string closed_only[] = {"cftp", "cftp.approx", "mdd"};
    for (const std::string& method : closed_only) {
        const FeatureSet f = line::qn::ctmc_feature_set(method);
        CHECK_MESSAGE(!f.has(Feature::OpenClass), method);
        CHECK_MESSAGE(!f.has(Feature::Source), method);
        CHECK_MESSAGE(!f.has(Feature::Sink), method);
        // ... and each keeps the closed core it does serve.
        CHECK_MESSAGE(f.has(Feature::ClosedClass), method);
        CHECK_MESSAGE(f.has(Feature::Queue), method);
        CHECK_MESSAGE(f.has(Feature::Delay), method);
        CHECK_MESSAGE(f.has(Feature::Exp), method);
    }
}

TEST_CASE("ctmc_feature_set: cftp keeps only the product-form envelope") {
    const FeatureSet f = line::qn::ctmc_feature_set("cftp");
    // The sampler walks Queue, Delay and Router alone, on a product-form
    // discipline, with Markovian routing and no rate scaling.
    CHECK_FALSE(f.has(Feature::Cache));
    CHECK_FALSE(f.has(Feature::ClassSwitch));
    CHECK_FALSE(f.has(Feature::Fork));
    CHECK_FALSE(f.has(Feature::Place));
    CHECK_FALSE(f.has(Feature::Region));
    CHECK_FALSE(f.has(Feature::LoadDependence));
    CHECK_FALSE(f.has(Feature::SchedStrategy_HOL));
    CHECK_FALSE(f.has(Feature::SchedStrategy_DPS));
    CHECK_FALSE(f.has(Feature::RoutingStrategy_JSQ));
    CHECK_FALSE(f.has(Feature::RoutingStrategy_RROBIN));
    CHECK(f.has(Feature::SchedStrategy_FCFS));
    CHECK(f.has(Feature::SchedStrategy_PS));
    CHECK(f.has(Feature::SchedStrategy_INF));
    CHECK(f.has(Feature::SchedStrategy_SIRO));
    CHECK(f.has(Feature::SchedStrategy_LCFSPR));
    CHECK(f.has(Feature::RoutingStrategy_PROB));
    // mdd narrows the open names ONLY: it holds a Petri net's marking too, so
    // its Place/Transition names have to survive.
    const FeatureSet mdd = line::qn::ctmc_feature_set("mdd");
    CHECK(mdd.has(Feature::Place));
    CHECK(mdd.has(Feature::Transition));
    CHECK(mdd.has(Feature::Cache));
}

TEST_CASE("ctmc_feature_set: the open model is refused, the closed one is kept") {
    const SN open_sn = mm1().get_struct();
    const SN closed_sn = repairmen().get_struct();
    const std::string closed_only[] = {"cftp", "cftp.approx", "mdd"};
    for (const std::string& method : closed_only) {
        const FeatureSet f = line::qn::ctmc_feature_set(method);
        CHECK_MESSAGE(!accepts(f, open_sn), method);
        CHECK_MESSAGE(accepts(f, closed_sn), method);
    }
    // The state-space methods take both; only their memory budget stops them.
    CHECK(accepts(line::qn::ctmc_feature_set("default"), open_sn));
    CHECK(accepts(line::qn::ctmc_feature_set("default"), closed_sn));
}

// ---------------------------------------------------------------------------
// SolverCTMC: what the registry CANNOT name
// ---------------------------------------------------------------------------

TEST_CASE("solver_ctmc_cftp_supports: the class count is structural") {
    CHECK(line::ctmc::solver_ctmc_cftp_supports(repairmen().get_struct()).empty());

    const std::string reason = line::ctmc::solver_ctmc_cftp_supports(cqn2().get_struct());
    CHECK_FALSE(reason.empty());
    CHECK(reason.find("single-class models only") != std::string::npos);
    CHECK(reason.find("2 classes") != std::string::npos);
}

TEST_CASE("solver_ctmc_cftp_supports: the gate and the run are one predicate") {
    // Whatever the predicate refuses, the solve refuses with the SAME sentence.
    const SN sn = cqn2().get_struct();
    const std::string reason = line::ctmc::solver_ctmc_cftp_supports(sn);
    REQUIRE_FALSE(reason.empty());
    line::ctmc::CtmcOptions opt;
    opt.method = "cftp";
    line::ctmc::CtmcCftpOptions cftpopt;
    cftpopt.samples = 100;
    try {
        line::ctmc::solver_ctmc_cftp(sn, opt, cftpopt);
        FAIL("solver_ctmc_cftp accepted a model its own gate refuses");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()) == reason);
    }
}

TEST_CASE("solver_ctmc_mdd_supports: the class count is structural") {
    CHECK(line::ctmc::solver_ctmc_mdd_supports(repairmen().get_struct()).empty());

    const std::string reason = line::ctmc::solver_ctmc_mdd_supports(cqn2().get_struct());
    CHECK_FALSE(reason.empty());
    CHECK(reason.find("single-class networks") != std::string::npos);
    CHECK(reason.find("2 classes") != std::string::npos);

    const std::string open_reason = line::ctmc::solver_ctmc_mdd_supports(mm1().get_struct());
    CHECK_FALSE(open_reason.empty());
    CHECK(open_reason.find("CLOSED networks") != std::string::npos);
}

TEST_CASE("solver_ctmc_mdd_supports: the gate and the run are one predicate") {
    const SN sn = cqn2().get_struct();
    const std::string reason = line::ctmc::solver_ctmc_mdd_supports(sn);
    REQUIRE_FALSE(reason.empty());
    line::ctmc::CtmcOptions opt;
    opt.method = "mdd";
    try {
        line::ctmc::solver_ctmc_mdd_analyzer(sn, opt);
        FAIL("solver_ctmc_mdd_analyzer accepted a model its own gate refuses");
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()) == reason);
    }
}

// ---------------------------------------------------------------------------
// SolverCTMC: the fork-join model class, which EVERY method has to clear
// ---------------------------------------------------------------------------

TEST_CASE("an undeclared pairing never reaches the gate in this port") {
    // THE PORT CLOSES THIS HOLE ONE LEVEL EARLIER than the reference does. The
    // pairing is a declaration on the Join in all four codebases, and a Join
    // built without it leaves the Fork unmatched; here `link` refuses to build
    // such a model at all, where MATLAB and native python accept it and leave
    // `sn.fj` empty for `sn_fj_validate` to refuse. Asserted rather than worked
    // around: it is why this port has no undeclared-pairing row to withdraw.
    NetD m("fj");
    const std::size_t f = m.add_fork("F");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t jn = m.add_join_unbound("J");
    const std::size_t d = m.add_delay("D");
    const std::size_t c = m.add_closed_class("C", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, jn, 1.0);
    P.set(c, c, q2, jn, 1.0);
    P.set(c, c, jn, d, 1.0);
    CHECK_THROWS_AS(
        {
            m.link(P);
            (void)m.get_struct();
        },
        line::Error);
}

TEST_CASE("sn_fj_supports: an open class through a Fork is refused") {
    // The pairing is declared here, so this is the second rule of the same
    // validator: the exact construction is stated for closed chains.
    const std::string reason =
        line::qn::sn_fj_supports(forkjoin(true, false).get_struct());
    CHECK_FALSE(reason.empty());

    // The converse: the exact construction IS derived for a declared closed
    // fork-join model, and refusing it would hide the only exact answer such a
    // model has.
    CHECK(line::qn::sn_fj_supports(forkjoin(true, true).get_struct()).empty());
    // A model with no fork at all is never asked about.
    CHECK(line::qn::sn_fj_supports(repairmen().get_struct()).empty());
    CHECK(line::qn::sn_fj_supports(mm1().get_struct()).empty());
}

TEST_CASE("sn_fj_supports: the gate and the run are one predicate") {
    const SN sn = forkjoin(true, false).get_struct();
    const std::string reason = line::qn::sn_fj_supports(sn);
    REQUIRE_FALSE(reason.empty());
    line::ctmc::CtmcOptions opt;
    opt.method = "default";
    try {
        line::ctmc::solver_ctmc_analyzer(sn, opt);
        FAIL("solver_ctmc_analyzer accepted a fork-join model its own validator refuses");
    } catch (const line::Error& e) {
        CHECK(std::string(e.what()) == reason);
    }
}

TEST_CASE("ctmc mdd: a fork-join model is neither of the two shapes it serves") {
    // The tag augmentation adds one auxiliary class per branch, so the struct
    // that reaches the analyzer is never single-class however the model was
    // written. Stated twice on purpose: as a feature, so the report names the
    // construct, and structurally, so the analyzer refuses before it builds.
    const FeatureSet mdd = line::qn::ctmc_feature_set("mdd");
    CHECK_FALSE(mdd.has(Feature::Fork));
    CHECK_FALSE(mdd.has(Feature::Join));
    CHECK_FALSE(mdd.has(Feature::Forker));
    CHECK_FALSE(mdd.has(Feature::Joiner));
    CHECK(line::qn::ctmc_feature_set("default").has(Feature::Fork));

    const SN sn = forkjoin(true, true).get_struct();
    CHECK_FALSE(accepts(mdd, sn));
    const std::string reason = line::ctmc::solver_ctmc_mdd_supports(sn);
    CHECK_FALSE(reason.empty());
    CHECK(reason.find("fork-join") != std::string::npos);
    // The state-space methods keep it: they are the exact answer here.
    CHECK(accepts(line::qn::ctmc_feature_set("default"), sn));
}

// ---------------------------------------------------------------------------
// SolverFLD
// ---------------------------------------------------------------------------

TEST_CASE("fluid_feature_set: refined is closed-only wherever the model is open") {
    // The reference's runAnalyzer has always refused this by name and this set
    // never said so; on an open fork-join model the restriction surfaced as a
    // failure inside the transform instead of as a refusal.
    const FeatureSet f = line::qn::fluid_feature_set("refined");
    CHECK_FALSE(f.has(Feature::OpenClass));
    CHECK_FALSE(f.has(Feature::Source));
    CHECK(f.has(Feature::ClosedClass));
    CHECK(line::qn::fluid_feature_set("minnormal").has(Feature::OpenClass));
    CHECK_FALSE(accepts(f, mm1().get_struct()));
    CHECK_FALSE(accepts(f, forkjoin(true, false).get_struct()));
    // ... and every closed model keeps it.
    CHECK(accepts(f, repairmen().get_struct()));
    CHECK(accepts(f, cqn2().get_struct()));
    // A CLOSED FORK-JOIN MODEL keeps it, and that is the half of the old
    // seven-name Fork rule that was RESOLVED IN THE OTHER DIRECTION: measured on
    // a symmetric closed fork-join, `refined` answers it symmetrically
    // (0.2715 0.2715 0.9713 0.4857 against the exact 0.664 0.664 0.624 1.024),
    // as MATLAB and native python have always done, so the name came back.
    CHECK(f.has(Feature::Fork));
    CHECK(accepts(f, forkjoin(true, true).get_struct()));
}

TEST_CASE("fluid_feature_set: the Fork names, measured method by method") {
    // THE RULE THAT CHANGED, and the evidence that changed it. This set used to
    // withhold Fork from seven methods on the argument that the MMT transform
    // hands the inner solve a MIXED network. Measured on a SYMMETRIC closed
    // fork-join against the exact chain (Q1 = Q2 = 0.664, J = 0.624, D = 1.024),
    // five of the seven answer it and answer it symmetrically, so their names
    // came back; MATLAB and native python have always run them there.
    const std::string integrates[] = {"statedep", "refined", "tbi", "mfq", "rmf"};
    const SN closed_fj = forkjoin(true, true).get_struct();
    for (const std::string& method : integrates) {
        const FeatureSet f = line::qn::fluid_feature_set(method);
        CHECK_MESSAGE(f.has(Feature::Fork), method);
        CHECK_MESSAGE(f.has(Feature::Join), method);
        CHECK_MESSAGE(accepts(f, closed_fj), method);
    }
    // The other two stay out, and for the opposite reason: each ANSWERS a
    // fork-join model rather than refusing one. `diffusion` put the whole
    // population on ONE station and zero elsewhere (a different station on a
    // rerun); `kp` returned an all-zero table on a symmetric OPEN fork-join fed
    // at rate 0.5. A silent wrong answer is the one outcome worth refusing, so
    // MATLAB and native python now withhold the names too.
    const std::string mis_answers[] = {"diffusion", "kp"};
    for (const std::string& method : mis_answers) {
        const FeatureSet f = line::qn::fluid_feature_set(method);
        CHECK_MESSAGE(!f.has(Feature::Fork), method);
        CHECK_MESSAGE(!f.has(Feature::Join), method);
        CHECK_MESSAGE(!accepts(f, closed_fj), method);
    }
    // The base set still declares them: this is a per-method narrowing.
    CHECK(line::qn::fluid_feature_set("default").has(Feature::Fork));
    CHECK(accepts(line::qn::fluid_feature_set("default"), closed_fj));
}

TEST_CASE("ssa_feature_set: the fork-join names are declared and the wiring gated") {
    // All four codebases declare these; native python alone omitted them, so its
    // ssa family never appeared on a fork-join model it solves. Declaring them
    // without the structural gate would move the defect rather than fix it.
    const FeatureSet f = line::qn::ssa_feature_set("default");
    CHECK(f.has(Feature::Fork));
    CHECK(f.has(Feature::Join));
    CHECK(f.has(Feature::Forker));
    CHECK(f.has(Feature::Joiner));
    CHECK(accepts(f, forkjoin(true, true).get_struct()));
    // ... and the wiring rule the names cannot state still refuses.
    CHECK_FALSE(line::qn::sn_fj_supports(forkjoin(true, false).get_struct()).empty());
}

TEST_CASE("fluid_forkjoin_supports: dae refuses an OPEN fork-join model only") {
    // The conjunction the feature set cannot state: Fork and OpenClass are both
    // declared names, and it is having BOTH that the DAE form cannot take --
    // the transform's auxiliary open classes have no unknown in it.
    const std::string reason =
        line::fluid::detail::fluid_forkjoin_supports(forkjoin(true, false).get_struct(), "dae");
    CHECK_FALSE(reason.empty());
    CHECK(reason.find("fork-join fixed point") != std::string::npos);
    CHECK(reason.find("minnormal") != std::string::npos);

    // Each half on its own stays runnable, and so does every other method.
    CHECK(line::fluid::detail::fluid_forkjoin_supports(mm1().get_struct(), "dae").empty());
    CHECK(line::fluid::detail::fluid_forkjoin_supports(forkjoin(true, true).get_struct(), "dae")
              .empty());
    const std::string others[] = {"default", "minnormal", "closing", "matrix", "statedep"};
    for (const std::string& method : others) {
        CHECK_MESSAGE(
            line::fluid::detail::fluid_forkjoin_supports(forkjoin(true, false).get_struct(), method)
                .empty(),
            method);
    }
}

TEST_CASE("fluid_feature_set: diffusion and tbi are closed-only") {
    const SN open_sn = mm1().get_struct();
    const SN closed_sn = repairmen().get_struct();
    const std::string closed_only[] = {"diffusion", "fluid.diffusion", "tbi", "fluid.tbi"};
    for (const std::string& method : closed_only) {
        const FeatureSet f = line::qn::fluid_feature_set(method);
        CHECK_MESSAGE(!f.has(Feature::OpenClass), method);
        CHECK_MESSAGE(!f.has(Feature::Source), method);
        CHECK_MESSAGE(f.has(Feature::ClosedClass), method);
        CHECK_MESSAGE(!accepts(f, open_sn), method);
        CHECK_MESSAGE(accepts(f, closed_sn), method);
    }
    // The methods that DO take an open model must not have moved.
    CHECK(line::qn::fluid_feature_set("default").has(Feature::OpenClass));
    CHECK(accepts(line::qn::fluid_feature_set("default"), open_sn));
}

TEST_CASE("fluid_feature_set: tbi refuses a cache station") {
    CHECK(line::qn::fluid_feature_set("default").has(Feature::Cache));
    CHECK_FALSE(line::qn::fluid_feature_set("tbi").has(Feature::Cache));
    // diffusion never declared a cache in the first place; the point of the
    // delta is the open names, so its closed core has to survive intact.
    CHECK(line::qn::fluid_feature_set("diffusion").has(Feature::SchedStrategy_PS));
    CHECK(line::qn::fluid_feature_set("diffusion").has(Feature::SchedStrategy_FCFS));
    CHECK(line::qn::fluid_feature_set("diffusion").has(Feature::SchedStrategy_INF));
}

TEST_CASE("fluid_qsys_horizon_supports: a trajectory needs a finite horizon") {
    line::fluid::FluidOptions unset;
    unset.timespan_end = std::numeric_limits<double>::infinity();
    const std::string time_varying[] = {"mol", "fluid.mol", "mtginf", "fluid.mtginf", "tvms"};
    for (const std::string& method : time_varying) {
        const std::string reason = line::fluid::fluid_qsys_horizon_supports(method, unset);
        CHECK_MESSAGE(!reason.empty(), method);
        CHECK_MESSAGE(reason.find("finite horizon") != std::string::npos, reason);
    }

    // The converse, and the whole point of the rule being on the OPTIONS: the
    // model never changed, only the horizon did.
    line::fluid::FluidOptions finite;
    finite.timespan_end = 10.0;
    for (const std::string& method : time_varying) {
        CHECK_MESSAGE(line::fluid::fluid_qsys_horizon_supports(method, finite).empty(), method);
    }

    // The stationary limits report a point, not a trajectory, so no horizon
    // rule may touch them; nor may it touch a method of the network family.
    const std::string stationary[] = {"ggisgi.fluid", "ggisgi", "ggingi.tga", "tga",
                                      "default", "closing", "minnormal"};
    for (const std::string& method : stationary) {
        CHECK_MESSAGE(line::fluid::fluid_qsys_horizon_supports(method, unset).empty(), method);
    }
}

TEST_CASE("fluid_qsys_horizon: the gate and the run are one predicate") {
    line::fluid::FluidOptions unset;
    unset.method = "mtginf";
    unset.timespan_end = std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(line::fluid::fluid_qsys_horizon_supports("mtginf", unset).empty());
    const SN sn = mm1().get_struct();
    double t0 = 0.0, t1 = 0.0;
    CHECK_THROWS_AS(line::fluid::detail::fluid_qsys_horizon(sn, unset, t0, t1), UnsupportedError);
    unset.timespan_end = 10.0;
    CHECK_NOTHROW(line::fluid::detail::fluid_qsys_horizon(sn, unset, t0, t1));
    CHECK(t1 == doctest::Approx(10.0));
}
