/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The SolverAUTO chooser: which engine a model is routed to, which slots are
 * skipped for want of an engine, and which selectors refuse.
 *
 * EVERY CASE PINS THE ARM THAT CLAIMS IT, because the arms of avgOrder are not
 * disjoint: a model that reaches the right solver through the wrong arm is a
 * latent defect that surfaces on the next model. Each case below therefore
 * carries the feature that makes its own arm the first match.
 *
 * THE FEATURE-SET GATE IS ASSERTED, not assumed. A ranking is a preference, and
 * the choice must move when the leading candidate cannot take the model -- the
 * MMPP2 and the two priority cases are there for that.
 *
 * THE FINDSTRING DEFECT IS ASSERTED ON PURPOSE. The reference's
 * hasHomogeneousScheduling degenerates to `nstations == 1` for every discipline
 * (see solver_auto.h), so the response-time-CDF branch fires on a single-station
 * model whatever its discipline. If someone repairs the predicate in
 * NetworkStruct, that test must fail: the repair changes which solver real
 * models get, and it should not pass silently.
 *
 * A SKIPPED SLOT IS PART OF THE ANSWER. LDES leads three of the reference's
 * rankings, so where its engine is absent the choice reports it in `skipped`
 * rather than presenting the runner-up as if it had been the first choice.
 *
 * THE LDES SLOT IS MACHINE-DEPENDENT, exactly like the LQNS one at the foot of
 * this file. `auto_solver_is_available(LDES)` is `ldes_is_available()`, which
 * answers whether this machine carries `common/ldes` or `common/ldes.jar` --
 * and cpp/make.sh installs the native engine beside line-cli, so a developer
 * tree normally HAS one. Pinning either outcome would make these cases pass or
 * fail by what is installed, so every ranking LDES leads is asserted against
 * the probe `ldes_here` instead.
 */

#include <algorithm>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/auto/solver_auto.h"

using namespace line;
using autosolver::AutoMode;
using autosolver::AutoSolver;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Delay -> Queue -> Delay, one closed class per chain: two stations, K chains. */
qn::Network<double> closed_chain(const std::string& name, SchedStrategy qsched,
                                 std::size_t nclasses, double njobs, double servers = 1.0) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", qsched);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 0; r < nclasses; ++r) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(r), njobs, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::exp_rate(2.0));
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
    }
    m.link(P);
    return m;
}

/** One Queue routing back to itself: a model whose only station is that queue. */
qn::Network<double> closed_single_station(const std::string& name, SchedStrategy sched,
                                          std::size_t nclasses, double njobs) {
    qn::Network<double> m(name);
    const std::size_t q = m.add_queue("Queue", sched);
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 0; r < nclasses; ++r) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(r), njobs, q);
        m.set_service(q, c, D::exp_rate(2.0));
        P.set(c, c, q, q, 1.0);
    }
    m.link(P);
    return m;
}

/** Source -> Queue -> Sink, one open class: the shape the qsys closed forms serve. */
qn::Network<double> open_mm1(const std::string& name, SchedStrategy sched = SchedStrategy::FCFS) {
    qn::Network<double> m(name);
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("O1");
    m.set_arrival(s, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Is `token` in the list? */
bool lists(const std::vector<std::string>& v, const std::string& token) {
    return std::find(v.begin(), v.end(), token) != v.end();
}

/** The bare family method names the list carries, in the order they appear. */
std::vector<std::string> families_of(const std::vector<std::string>& v) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i].find('.') == std::string::npos && lists(autosolver::auto_network_family_names(), v[i]))
            out.push_back(v[i]);
    return out;
}

std::string avg_choice(qn::Network<double>& m) {
    return autosolver::auto_solver_name(autosolver::auto_choose_avg_solver(m.get_struct()));
}

std::string choice(qn::Network<double>& m, const std::string& method) {
    return autosolver::auto_solver_name(autosolver::auto_choose_solver(m.get_struct(), method));
}

std::string mode_choice(qn::Network<double>& m, const std::string& method, AutoMode mode) {
    return autosolver::auto_solver_name(
        autosolver::auto_choose_solver_mode(m.get_struct(), method, mode).solver);
}

/** Does this machine carry an LDES engine? See the file header. */
bool ldes_here() { return line::ldes::ldes_is_available(); }

/**
 * The outcome of a ranking LDES leads: `ldes` and nothing skipped where the
 * engine is installed, the runner-up plus a reported skip where it is not.
 */
void check_ldes_led(const autosolver::AutoChoice& c, const char* runner_up) {
    if (ldes_here()) {
        CHECK(autosolver::auto_solver_name(c.solver) == std::string("ldes"));
        CHECK(c.skipped.empty());
    } else {
        CHECK(autosolver::auto_solver_name(c.solver) == std::string(runner_up));
        REQUIRE(c.skipped.size() == 1);
        CHECK(autosolver::auto_solver_name(c.skipped[0]) == std::string("ldes"));
    }
}

}  // namespace

TEST_CASE("solver-auto: the traits the rankings are keyed on") {
    qn::Network<double> pf = closed_chain("tr", SchedStrategy::PS, 2, 4.0);
    const autosolver::AutoTraits t = autosolver::auto_traits(pf.get_struct());
    CHECK(t.is_closed);
    CHECK_FALSE(t.is_open);
    CHECK(t.is_product_form);
    CHECK_FALSE(t.has_multi_server);
    CHECK(t.total_jobs == doctest::Approx(8.0));
    CHECK(t.pop_per_chain == doctest::Approx(4.0));
    CHECK(t.prio == "none");
    CHECK_FALSE(t.has_map);

    qn::Network<double> hol = closed_chain("holtr", SchedStrategy::HOL, 2, 4.0);
    CHECK(autosolver::auto_traits(hol.get_struct()).prio == "hol");
    qn::Network<double> pp = closed_chain("pptr", SchedStrategy::FCFSPRPRIO, 2, 4.0);
    CHECK(autosolver::auto_traits(pp.get_struct()).prio == "preempt");
    qn::Network<double> ms = closed_chain("mstr", SchedStrategy::PS, 2, 4.0, 3.0);
    CHECK(autosolver::auto_traits(ms.get_struct()).has_multi_server);
}

TEST_CASE("solver-auto: chooseAvgSolverHeur, arm by arm") {
    SUBCASE("a small population takes an exact solver, MVA first") {
        // totalJobs = 3 <= EXACT_POPULATION_MAX and the model is product form,
        // so method 'exact' is admissible and the exact-first pass claims it
        // before avgOrder is consulted at all.
        qn::Network<double> m = closed_chain("small", SchedStrategy::PS, 1, 3.0);
        const autosolver::AutoChoice c = autosolver::auto_choose_avg_solver_ex(m.get_struct());
        CHECK(autosolver::auto_solver_name(c.solver) == std::string("mva"));
        CHECK(c.method == "exact");
    }
    SUBCASE("a large population takes the fluid arm") {
        // popPerChain = 40 > 30, and no earlier arm claims the model.
        qn::Network<double> m = closed_chain("big", SchedStrategy::PS, 2, 40.0);
        CHECK(avg_choice(m) == "fld");
    }
    SUBCASE("a mid population falls to the global order, MVA") {
        qn::Network<double> m = closed_chain("mid", SchedStrategy::PS, 2, 12.0);
        CHECK(avg_choice(m) == "mva");
    }
    SUBCASE("autocorrelated service promotes MAM over MVA") {
        qn::Network<double> m("mapmodel");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("C0", 12.0, d);
        m.set_service(d, c, D::exp_rate(1.0));
        {
            Matrix<double> D0(2, 2), D1(2, 2);
            const double l0 = 2.0, l1 = 1.0, sg0 = 0.5, sg1 = 0.3;
            D1(0, 0) = l0;
            D1(1, 1) = l1;
            D0(0, 0) = -(l0 + sg0);
            D0(0, 1) = sg0;
            D0(1, 0) = sg1;
            D0(1, 1) = -(l1 + sg1);
            m.set_service(q, c, D::map_dist(D0, D1, lang::ProcessType::MMPP2));
        }
        qn::RoutingMatrix<double> P;
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
        m.link(P);
        CHECK(autosolver::auto_traits(m.get_struct()).has_map);
        CHECK(avg_choice(m) == "mam");
    }
    SUBCASE("head-of-line priority keeps MVA, which declares HOL") {
        qn::Network<double> m = closed_chain("hol", SchedStrategy::HOL, 2, 12.0);
        CHECK(avg_choice(m) == "mva");
    }
    SUBCASE("preemptive priority leads with MAM, and reports what it skipped") {
        // avgOrder is [MAM, LDES, CTMC, SSA]. MAM does not declare FCFSPRPRIO,
        // so the feature gate rejects it; LDES answers where its engine is
        // installed, and is skipped for CTMC where it is not.
        qn::Network<double> m = closed_chain("preempt", SchedStrategy::FCFSPRPRIO, 2, 3.0);
        const autosolver::AutoChoice c = autosolver::auto_choose_avg_solver_ex(m.get_struct());
        check_ldes_led(c, "ctmc");
    }
}

TEST_CASE("solver-auto: the ranked gate is a feature test, not a preference") {
    // A model MVA cannot take must not be handed to MVA whatever the ranking
    // says: PSPRIO is not in mva_feature_set, so the ps arm's [CTMC, LDES, SSA]
    // is what answers.
    qn::Network<double> m = closed_chain("psprio", SchedStrategy::PSPRIO, 2, 3.0);
    CHECK_FALSE(autosolver::auto_supports(AutoSolver::MVA, m.get_struct(), "default"));
    CHECK(avg_choice(m) == "ctmc");
}

TEST_CASE("solver-auto: method 'exact' needs a product-form solution") {
    qn::Network<double> pf = closed_chain("pfx", SchedStrategy::PS, 2, 3.0);
    CHECK(autosolver::auto_supports(AutoSolver::MVA, pf.get_struct(), "exact"));
    CHECK(autosolver::auto_supports(AutoSolver::NC, pf.get_struct(), "exact"));

    qn::Network<double> npf = closed_chain("npfx", SchedStrategy::HOL, 2, 3.0);
    CHECK_FALSE(autosolver::auto_supports(AutoSolver::MVA, npf.get_struct(), "exact"));
    CHECK_FALSE(autosolver::auto_supports(AutoSolver::NC, npf.get_struct(), "exact"));
    // The exact SELECTION mode still answers, from the chain.
    CHECK(mode_choice(npf, "getAvgTable", AutoMode::EXACT) == "ctmc");
}

TEST_CASE("solver-auto: the CTMC slot is screened for a chain that fits") {
    qn::Network<double> small = closed_chain("fits", SchedStrategy::PS, 1, 3.0);
    CHECK(autosolver::auto_ctmc_is_tractable(small.get_struct()));
    CHECK(autosolver::auto_supports(AutoSolver::CTMC, small.get_struct(), "default"));

    // Three classes of 4000 jobs over two stations: the placement factor alone
    // is far past the generator's cap.
    qn::Network<double> huge = closed_chain("huge", SchedStrategy::PS, 3, 4000.0);
    CHECK_FALSE(autosolver::auto_ctmc_is_tractable(huge.get_struct()));
    CHECK_FALSE(autosolver::auto_supports(AutoSolver::CTMC, huge.get_struct(), "default"));
}

TEST_CASE("solver-auto: the getter name selects the family before the tree runs") {
    qn::Network<double> pf = closed_chain("fampf", SchedStrategy::PS, 2, 12.0);

    SUBCASE("the average-metric family defers to the feature tree") {
        CHECK(choice(pf, "getAvgTable") == "mva");
        CHECK(choice(pf, "getAvg") == "mva");
        CHECK(choice(pf, "getAvgSysRespT") == "mva");
    }
    SUBCASE("transient probabilities are CTMC unconditionally") {
        CHECK(choice(pf, "getTranProb") == "ctmc");
        CHECK(choice(pf, "getTranProbSysAggr") == "ctmc");
    }
    SUBCASE("per-job sampling is SSA, the native sample-path engine") {
        CHECK(choice(pf, "sample") == "ssa");
        CHECK(choice(pf, "sampleSys") == "ssa");
    }
    SUBCASE("aggregate sampling prefers LDES and falls to SSA without it") {
        const autosolver::AutoChoice c =
            autosolver::auto_choose_solver_heur(pf.get_struct(), "sampleAggr");
        check_ldes_led(c, "ssa");
    }
    SUBCASE("transient averages take Fluid") {
        CHECK(choice(pf, "getTranAvg") == "fld");
        CHECK(choice(pf, "getTranCdfRespT") == "fld");
    }
    SUBCASE("state probabilities take NC while the model is product form") {
        CHECK(choice(pf, "getProb") == "nc");
        CHECK(choice(pf, "getProbNormConstAggr") == "nc");
    }
    SUBCASE("the response-time CDF needs homogeneous scheduling as well") {
        // Two stations, so the homogeneous conjunct fails and Fluid takes it
        // even though the model is product form and FCFS-solvable.
        CHECK(choice(pf, "getCdfRespT") == "fld");
        qn::Network<double> one = closed_single_station("cdf1", SchedStrategy::FCFS, 1, 12.0);
        CHECK(choice(one, "getCdfRespT") == "nc");
    }
    SUBCASE("the sensitivity family wants a differentiable model") {
        CHECK(choice(pf, "getSensitivityTable") == "fld");
    }
    SUBCASE("moments rank the analytical solvers first") {
        CHECK(choice(pf, "getMomentTable") == "mva");
    }
    SUBCASE("drop counts come from a sample path or an exact chain, never AMVA") {
        // [LDES, CTMC, SSA]: the sample path where it exists, the chain where
        // it does not. Never AMVA, which cannot count a drop.
        const autosolver::AutoChoice c =
            autosolver::auto_choose_solver_heur(pf.get_struct(), "getAvgLossTable");
        check_ldes_led(c, "ctmc");
    }
    SUBCASE("an ensemble getter is a LayeredNetwork method and has no Network arm") {
        CHECK_THROWS_AS(choice(pf, "getEnsembleAvg"), InputError);
        CHECK_THROWS_AS(choice(pf, "getNumberOfModels"), InputError);
    }
    SUBCASE("a getter with no ranking of its own falls to the average heuristic") {
        // The reference's `otherwise` arm, which is a fallback and not an error.
        CHECK(choice(pf, "getSomethingElse") == "mva");
    }
}

TEST_CASE("solver-auto: the selection modes") {
    qn::Network<double> pf = closed_chain("modes", SchedStrategy::PS, 2, 12.0);
    CHECK(mode_choice(pf, "getAvgTable", AutoMode::HEUR) == "mva");
    CHECK(mode_choice(pf, "getAvgTable", AutoMode::FAST) == "mva");
    CHECK(mode_choice(pf, "getAvgTable", AutoMode::ACCURATE) == "fld");
    CHECK(mode_choice(pf, "getAvgTable", AutoMode::EXACT) == "nc");
    // The simulator ranking is [LDES, SSA]; which one answers is which engine
    // this machine carries.
    CHECK(mode_choice(pf, "getAvgTable", AutoMode::SIM) == (ldes_here() ? "ldes" : "ssa"));
}

TEST_CASE("solver-auto: method tokens are intents or families") {
    CHECK(autosolver::auto_resolve_token("").is_intent);
    CHECK(autosolver::auto_resolve_token("exact").mode == AutoMode::EXACT);
    CHECK(autosolver::auto_resolve_token("accurate").mode == AutoMode::ACCURATE);

    const autosolver::AutoToken t = autosolver::auto_resolve_token("nc.comom");
    CHECK_FALSE(t.is_intent);
    CHECK(t.family == "nc");
    CHECK(t.submethod == "comom");

    const autosolver::AutoToken bare = autosolver::auto_resolve_token("fld");
    CHECK(bare.family == "fluid");
    CHECK(bare.submethod == "default");

    // A bare algorithm name needs the per-family registry this port lacks.
    CHECK_THROWS_AS(autosolver::auto_resolve_token("comom"), InputError);
}

TEST_CASE("solver-auto: the candidate fallback order drops what cannot run") {
    qn::Network<double> m = closed_chain("cand", SchedStrategy::PS, 2, 12.0);
    const std::vector<AutoSolver> c = autosolver::auto_candidates(m.get_struct());
    REQUIRE(c.size() >= 4);
    CHECK(autosolver::auto_solver_name(c[0]) == std::string("mva"));
    CHECK(autosolver::auto_solver_name(c[1]) == std::string("nc"));
    for (AutoSolver s : c) CHECK(autosolver::auto_solver_is_available(s));
    // The LDES slot appears exactly when its engine does -- the list drops what
    // cannot run, and nothing else.
    CHECK((std::find(c.begin(), c.end(), AutoSolver::LDES) != c.end()) == ldes_here());

    // The delegate's order is the chosen solver, then the candidates.
    const std::vector<AutoSolver> prop =
        autosolver::auto_proposed_solvers(m.get_struct(), "getTranAvg", AutoMode::HEUR);
    CHECK(autosolver::auto_solver_name(prop[0]) == std::string("fld"));
    CHECK(prop.size() >= c.size());
}

TEST_CASE("solver-auto: JMT is not a candidate, and LDES is skipped not refused") {
    qn::Network<double> m = closed_chain("skip", SchedStrategy::PS, 2, 12.0);
    // The LDES slot answers the engine probe and nothing else -- it is never
    // hard-wired to "unported", which would silently disagree with the
    // reference on every machine that has the engine.
    CHECK(autosolver::auto_solver_is_available(AutoSolver::LDES) == ldes_here());
    CHECK(autosolver::auto_solver_is_available(AutoSolver::MVA));
    // Every ranking that leads with LDES answers either way: from LDES itself,
    // or from the next slot with the skip reported.
    CHECK_NOTHROW(choice(m, "getAvgLossTable"));
    CHECK_NOTHROW(choice(m, "sampleAggr"));
}

TEST_CASE("solver-auto: the layered arm ranks the LN layer engines") {
    // THE LQNS SLOT IS MACHINE-DEPENDENT, in the reference too: LINE ships no
    // lqns binary and `chooseAvgSolverHeur.m` gates the candidate on
    // `SolverLQNS.isAvailable()`. The two rankings LQNS leads are therefore
    // asserted against the probe, not against a fixed name -- pinning either
    // outcome alone would make this test pass or fail by what is installed.
    const bool lqns_here = line::lqns::lqns_is_available();
    CHECK(autosolver::auto_layered_name(
              autosolver::auto_choose_layered_solver("getAvgTable", false).solver) ==
          std::string(lqns_here ? "lqns" : "ln.mva"));
    // A cache task inverts the order: the cache layer is where NC beats MVA.
    CHECK(autosolver::auto_layered_name(
              autosolver::auto_choose_layered_solver("getAvgTable", true).solver) ==
          std::string("ln.comom"));
    CHECK(autosolver::auto_layered_name(
              autosolver::auto_choose_layered_solver("getTranAvg", false).solver) ==
          std::string("ln.fluid"));
    CHECK(autosolver::auto_layered_name(
              autosolver::auto_choose_layered_solver("getCdfRespT", false).solver) ==
          std::string("ln.fluid"));
    CHECK(autosolver::auto_layered_name(
              autosolver::auto_choose_layered_solver("getAvgTable", false, AutoMode::EXACT)
                  .solver) == std::string("ln.comom"));
    // LQNS leads the sampling ranking: lqsim is the layered simulator. Where the
    // binary is absent the slot is skipped, and the skip is REPORTED rather than
    // swallowed, so the caller can see the second choice answered.
    const autosolver::AutoLayeredChoice sim =
        autosolver::auto_choose_layered_solver("sampleSys", false);
    if (lqns_here) {
        CHECK(autosolver::auto_layered_name(sim.solver) == std::string("lqns"));
        CHECK(sim.skipped.empty());
    } else {
        CHECK(autosolver::auto_layered_name(sim.solver) == std::string("ln.mva"));
        REQUIRE(sim.skipped.size() == 1);
        CHECK(autosolver::auto_layered_name(sim.skipped[0]) == std::string("lqns"));
    }
}

TEST_CASE("solver-auto: the environment arm leads with the transient stage solver") {
    const autosolver::AutoEnvChoice c = autosolver::auto_choose_env_solver("getEnsembleAvg");
    CHECK(autosolver::auto_env_name(c.solver) == std::string("env.fluid"));
    CHECK(c.skipped.empty());
    // The exact ranking is [NC, MVA] and neither stage engine exists here.
    CHECK_THROWS_AS(autosolver::auto_choose_env_solver("getEnsembleAvg", AutoMode::EXACT),
                    UnsupportedError);
}

/**
 * `listValidMethods` answers about THE MODEL, and the feature set is what makes
 * it do so.
 *
 * The cases below are chosen so that each isolates one of the three gates the
 * list runs: the flat feature set (a discipline no analytical family declares),
 * the model-aware registry (the queueing-system closed forms, which exist for a
 * two-station open model and nowhere else) and the structural rules a feature
 * name cannot express (a binding buffer). A list that widened to the token
 * universe would pass none of them.
 */
TEST_CASE("solver-auto: listValidMethods is filtered by the feature set") {
    SUBCASE("the selection intents are always offered") {
        // They name a RANKING, not an algorithm: the ranking's own job is to
        // find a family that supports the model, so an intent is never dropped.
        qn::Network<double> m = closed_chain("pf", SchedStrategy::PS, 2, 3.0);
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(m.get_struct());
        for (const char* intent :
             {"default", "heur", "sim", "exact", "fast", "accurate", "bound", "auto"})
            CHECK(lists(v, intent));
    }

    SUBCASE("a family whose feature set refuses the model contributes nothing") {
        // PSPRIO is declared by no analytical family: MVA, NC, MAM and the
        // bounds all drop out, and what remains is the state-space and
        // simulation side. The same gate `auto_supports` applies to the ranked
        // choice is what removes them here.
        qn::Network<double> m = closed_chain("psprio", SchedStrategy::PSPRIO, 2, 3.0);
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(m.get_struct());
        CHECK_FALSE(lists(v, "mva"));
        CHECK_FALSE(lists(v, "mva.default"));
        CHECK_FALSE(lists(v, "nc"));
        CHECK_FALSE(lists(v, "ba"));
        CHECK(lists(v, "ctmc.default"));

        // ... and it comes back the moment the discipline does.
        qn::Network<double> ps = closed_chain("ps", SchedStrategy::PS, 2, 3.0);
        const std::vector<std::string> w = autosolver::auto_list_valid_methods(ps.get_struct());
        CHECK(lists(w, "mva"));
        CHECK(lists(w, "mva.default"));
        CHECK(lists(w, "nc.comom"));
    }

    SUBCASE("a method the family withholds on this model is not offered") {
        // The queueing-system closed forms are advertised for a single-class
        // two-station OPEN model, which is the shape their analyzer serves.
        qn::Network<double> open = open_mm1("qsys");
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(open.get_struct());
        CHECK(lists(v, "mva.mm1"));
        CHECK(lists(v, "mva.qna"));
        // 'marie' and 'mvac' are the closed-model counterparts, and are not.
        CHECK_FALSE(lists(v, "mva.marie"));
        CHECK_FALSE(lists(v, "mva.mvac"));

        qn::Network<double> closed = closed_chain("noqsys", SchedStrategy::PS, 1, 3.0);
        const std::vector<std::string> w = autosolver::auto_list_valid_methods(closed.get_struct());
        CHECK_FALSE(lists(w, "mva.mm1"));
        CHECK_FALSE(lists(w, "mva.qna"));
        CHECK(lists(w, "mva.mvac"));
    }

    SUBCASE("the bounds narrow to the side of the model they are derived for") {
        // 'bpt', 'bgt' and 'snc' are the OPEN-network bounds; every other family
        // in SolverBA rules the open model out, and the reverse holds on a
        // closed one. The registry is asked for the model, so the list follows.
        qn::Network<double> open = open_mm1("obound");
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(open.get_struct());
        CHECK(lists(v, "ba.bgt.upper"));
        CHECK(lists(v, "ba.bpt.lower"));
        CHECK_FALSE(lists(v, "ba.aba.upper"));

        qn::Network<double> closed = closed_chain("cbound", SchedStrategy::PS, 1, 3.0);
        const std::vector<std::string> w = autosolver::auto_list_valid_methods(closed.get_struct());
        CHECK(lists(w, "ba.aba.upper"));
        CHECK_FALSE(lists(w, "ba.bgt.upper"));
    }

    SUBCASE("a binding buffer is a refusal no feature name describes") {
        // setCapacity is not a feature: MVA, NC and FLD refuse a capped station
        // through check_binding_capacity, because none of them represents a
        // finite buffer and all three would answer the UNCONSTRAINED model.
        // The list has to run that gate itself or it advertises exactly that.
        qn::Network<double> m = closed_chain("capped", SchedStrategy::FCFS, 1, 4.0);
        const std::vector<std::string> before =
            autosolver::auto_list_valid_methods(m.get_struct());
        CHECK(lists(before, "mva.default"));

        qn::Network<double> c("capped2");
        const std::size_t d = c.add_delay("Delay");
        const std::size_t q = c.add_queue("Queue", SchedStrategy::FCFS);
        const std::size_t k = c.add_closed_class("C0", 4.0, d);
        c.set_service(d, k, D::exp_rate(1.0));
        c.set_service(q, k, D::exp_rate(2.0));
        c.set_capacity(q, 2.0);
        qn::RoutingMatrix<double> P;
        P.set(k, k, d, q, 1.0);
        P.set(k, k, q, d, 1.0);
        c.link(P);
        const std::vector<std::string> after =
            autosolver::auto_list_valid_methods(c.get_struct());
        CHECK_FALSE(lists(after, "mva"));
        CHECK_FALSE(lists(after, "mva.default"));
        CHECK_FALSE(lists(after, "fluid"));
        // CTMC carries the buffer in its state space and keeps its tokens.
        CHECK(lists(after, "ctmc.default"));
    }

    SUBCASE("every listed token is a family this port can dispatch") {
        // A listed name must actually run: the list may not name a family the
        // AUTO method name resolver refuses, nor one that takes another model class.
        qn::Network<double> m = closed_chain("dispatch", SchedStrategy::PS, 2, 3.0);
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(m.get_struct());
        for (std::size_t i = 0; i < v.size(); ++i) {
            const std::string& tok = v[i];
            const autosolver::AutoToken t = autosolver::auto_resolve_token(tok);
            if (t.is_intent) continue;
            CHECK(lists(autosolver::auto_network_family_names(), t.family));
        }
        for (const char* absent : {"ln", "env", "lqns", "uq"}) CHECK_FALSE(lists(v, absent));
    }

    SUBCASE("a bare family token stands only where one of its methods does") {
        qn::Network<double> m = closed_chain("bare", SchedStrategy::PS, 2, 3.0);
        const std::vector<std::string> v = autosolver::auto_list_valid_methods(m.get_struct());
        const std::vector<std::string> fams = families_of(v);
        REQUIRE_FALSE(fams.empty());
        for (std::size_t i = 0; i < fams.size(); ++i) {
            bool any = false;
            for (std::size_t j = 0; j < v.size() && !any; ++j)
                any = v[j].compare(0, fams[i].size() + 1, fams[i] + ".") == 0;
            CHECK(any);
        }
    }
}
