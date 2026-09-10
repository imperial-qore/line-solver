/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Every LQN layer engine must agree about WHICH CELLS OF THE TABLE EXIST.
 *
 * A NaN in an LQN average table is not a failed computation: it says the
 * quantity is not defined for that element kind. A processor has no queue
 * length, response time, residence or throughput of its own; an entry has no
 * residence; nothing reports an arrival rate. Two engines that disagree about
 * the mask disagree about the MODEL rather than about arithmetic, and a
 * tolerance-based comparison cannot see it -- `parity-static/_compare_vendor.py`
 * passes any cell where either side is NaN, which is why the goldens carry the
 * string "NaN" beside the numbers.
 *
 * In this port the mask is carried by the `defined_*` flags of `LnSolution`
 * rather than by a NaN payload; `line-cli` prints "NaN" wherever a flag is
 * false. Twin of `jar/src/test/java/jline/solvers/ln/LqnNaNMaskParityTest.java`,
 * of `python/tests/test_lqn_nan_mask_parity.py` and of
 * `line-test.git/test/testsMisc/test_lqn_nan_mask_parity.m`. LQNS is
 * deliberately NOT in the panel: it reports no residence time and no arrival
 * rate at all, and leaves both undefined by decision rather than by omission
 * (see `_kb/06-solver-catalog.md`, LQNS wrapper) -- and its binary is not always
 * present.
 */

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ldes/ldes_ln_engine.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

const double INF = std::numeric_limits<double>::infinity();

/** The model of `matlab/examples/basic/layeredModel/lqn_multi_solvers.m`. */
lqn::LqnStruct<double> build_multi_solvers() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", INF, SchedStrategy::INF);
    b.processor("P2", INF, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", E(0.0001));
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", E(1.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 3.0);
    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    return b.build();
}

/**
 * The same model plus a component NOTHING REACHES: no reference task, no
 * caller.
 *
 * SolverLN marks such a component `ignore` and reports zero for it -- which is
 * right for the measures its elements HAVE, and was being written over the ones
 * they do not, flattening the mask. The `lqn_ofbiz` golden carries exactly this
 * shape (its USAGE_DELAY component), which is why the defect survived the
 * 2026-08-21 audit: no golden pairs LQNS with an LN row on a model that has one.
 */
lqn::LqnStruct<double> build_ignored_component() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", INF, SchedStrategy::INF);
    b.processor("P2", INF, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", E(1.0));
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", E(1.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", E(1.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    // unreachable: T3 is not a reference task and nobody calls E3
    b.processor("P3", INF, SchedStrategy::INF);
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E3", "T3");
    b.activity("A3", E(1.0), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** "P1.QLen" -> true when the cell is undefined. The mask, spelled so a failure reads. */
std::map<std::string, bool> mask_of(const lqn::LqnStruct<double>& l,
                                    const ln::LnSolution<double>& s) {
    std::map<std::string, bool> m;
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        const std::string& n = l.names[i];
        m[n + ".QLen"] = !s.defined_Q[i];
        m[n + ".Util"] = !s.defined_U[i];
        m[n + ".RespT"] = !s.defined_R[i];
        m[n + ".ResidT"] = !s.defined_W[i];
        m[n + ".ArvR"] = !s.defined_A[i];
        m[n + ".Tput"] = !s.defined_T[i];
    }
    return m;
}

ln::LnSolution<double> solve(const lqn::LqnStruct<double>& l, const std::string& engine,
                             const std::string& method = "") {
    ln::LnOptions opt;
    opt.layer_solver = engine;
    if (!method.empty()) opt.method = method;
    ln::SolverLN<double> s(l, opt);
    return s.get_ensemble_avg();
}

void check_same_mask(const std::string& ref_name, const std::map<std::string, bool>& ref,
                     const std::string& other_name, const std::map<std::string, bool>& other) {
    std::string disagree;
    for (std::map<std::string, bool>::const_iterator it = ref.begin(); it != ref.end(); ++it) {
        std::map<std::string, bool>::const_iterator o = other.find(it->first);
        if (o == other.end()) continue;  // a solver that omits a column is not a mask defect
        if (o->second != it->second) {
            if (!disagree.empty()) disagree += "; ";
            disagree += it->first + " (" + ref_name + (it->second ? "=NaN" : "=value") + ", " +
                        other_name + (o->second ? "=NaN" : "=value") + ")";
        }
    }
    INFO("LQN table NaN mask differs between " << ref_name << " and " << other_name << ": "
                                               << disagree);
    CHECK(disagree.empty());
}

}  // namespace

TEST_CASE("lqn mask: every layer engine agrees on which cells exist") {
    const lqn::LqnStruct<double> l = build_multi_solvers();
    const std::map<std::string, bool> mva = mask_of(l, solve(l, "mva"));
    const std::map<std::string, bool> nc = mask_of(l, solve(l, "nc"));
    check_same_mask("LN(MVA)", mva, "LN(NC)", nc);
}

/**
 * The mask itself, pinned. Without the test above this would still pass if
 * every engine regressed the same way.
 */
TEST_CASE("lqn mask: the mask is the documented one") {
    const lqn::LqnStruct<double> l = build_multi_solvers();
    const std::map<std::string, bool> mask = mask_of(l, solve(l, "mva"));

    // A processor has no queue, no response time, no residence, no throughput
    // of its own -- only a utilization. TN is left undefined explicitly, "for
    // consistency with LQNS".
    const char* procs[] = {"P1", "P2"};
    for (int i = 0; i < 2; ++i) {
        const std::string p(procs[i]);
        CHECK(mask.at(p + ".QLen"));
        CHECK(mask.at(p + ".RespT"));
        CHECK(mask.at(p + ".ResidT"));
        CHECK(mask.at(p + ".Tput"));
        CHECK_FALSE(mask.at(p + ".Util"));
    }
    // An entry has a response time but no residence.
    const char* entries[] = {"E1", "E2"};
    for (int i = 0; i < 2; ++i) {
        const std::string e(entries[i]);
        CHECK(mask.at(e + ".ResidT"));
        CHECK_FALSE(mask.at(e + ".RespT"));
    }
    // Tasks and activities carry both.
    const char* both[] = {"T1", "T2", "A1", "A2"};
    for (int i = 0; i < 4; ++i) CHECK_FALSE(mask.at(std::string(both[i]) + ".ResidT"));
    // Nobody reports an arrival rate on an LQN.
    const char* any[] = {"P1", "T1", "E1", "A1"};
    for (int i = 0; i < 4; ++i) CHECK(mask.at(std::string(any[i]) + ".ArvR"));
}

/**
 * A component no reference task reaches is IDLE, not undefined.
 *
 * Its elements report zero for every measure their kind has, and NaN for the
 * ones it does not -- the same mask a reachable element carries. Until
 * 2026-08-25 the ignore branch marked all six columns defined and wrote a flat
 * zero, so an unreachable processor claimed a queue length of 0 where every
 * solver, LQNS included, reports none.
 */
TEST_CASE("lqn mask: an ignored component keeps the mask") {
    const lqn::LqnStruct<double> li = build_ignored_component();
    const lqn::LqnStruct<double> lr = build_multi_solvers();
    const ln::LnSolution<double> si = solve(li, "mva");
    const std::map<std::string, bool> mask = mask_of(li, si);
    const std::map<std::string, bool> reachable = mask_of(lr, solve(lr, "mva"));

    // the unreachable rows carry the SAME mask as the reachable ones
    const char* pairs[][3] = {{"P3", "P1", "QLen"},  {"P3", "P1", "RespT"},
                              {"P3", "P1", "ResidT"}, {"P3", "P1", "Tput"},
                              {"P3", "P1", "Util"},   {"T3", "T2", "RespT"},
                              {"T3", "T2", "ResidT"}, {"E3", "E2", "RespT"},
                              {"E3", "E2", "ResidT"}, {"A3", "A2", "ResidT"}};
    for (int i = 0; i < 10; ++i) {
        const std::string unreached = std::string(pairs[i][0]) + "." + pairs[i][2];
        const std::string reached = std::string(pairs[i][1]) + "." + pairs[i][2];
        INFO(unreached << " must carry the mask of " << reached);
        CHECK(mask.at(unreached) == reachable.at(reached));
    }

    // spelled out, so the loop above cannot pass by both sides regressing
    CHECK(mask.at("P3.QLen"));
    CHECK(mask.at("P3.RespT"));
    CHECK(mask.at("P3.Tput"));
    CHECK_FALSE(mask.at("P3.Util"));
    CHECK(mask.at("T3.RespT"));
    CHECK_FALSE(mask.at("T3.ResidT"));
    CHECK(mask.at("E3.ResidT"));
    const char* unreachable[] = {"P3", "T3", "E3", "A3"};
    for (int i = 0; i < 4; ++i) CHECK(mask.at(std::string(unreachable[i]) + ".ArvR"));

    // idle, not undefined: the cells that DO exist there read zero
    for (std::size_t i = 1; i <= li.nidx; ++i) {
        const std::string& n = li.names[i];
        if (n != "P3" && n != "T3" && n != "E3" && n != "A3") continue;
        if (si.defined_Q[i]) CHECK(si.QN[i] == doctest::Approx(0.0));
        if (si.defined_U[i]) CHECK(si.UN[i] == doctest::Approx(0.0));
        if (si.defined_R[i]) CHECK(si.RN[i] == doctest::Approx(0.0));
        if (si.defined_W[i]) CHECK(si.WN[i] == doctest::Approx(0.0));
        if (si.defined_T[i]) CHECK(si.TN[i] == doctest::Approx(0.0));
    }

    // and every layer engine still agrees about it
    check_same_mask("LN(MVA)", mask, "LN(NC)", mask_of(li, solve(li, "nc")));
}

/**
 * `srvn.ph` assembles the table in its own routine (`aggregate_ph`), and had its
 * own copy of the flat-zero branch.
 *
 * The two encodings rebuild every figure differently -- one reads class rows off
 * the ensemble, the other composes a phase-type law -- so agreeing on the mask is
 * a claim about the table, not about shared code.
 */
TEST_CASE("lqn mask: the ph encoding masks an ignored component the same way") {
    const lqn::LqnStruct<double> l = build_ignored_component();
    const std::map<std::string, bool> dflt = mask_of(l, solve(l, "mva"));
    const std::map<std::string, bool> ph = mask_of(l, solve(l, "mva", "srvn.ph"));
    check_same_mask("LN(MVA)", dflt, "LN(MVA, srvn.ph)", ph);
    CHECK(ph.at("P3.QLen"));
    CHECK(ph.at("P3.ArvR"));
    CHECK_FALSE(ph.at("P3.Util"));
}

/**
 * The NATIVE LDES LAYERED ENGINE joins the panel.
 *
 * It is the one engine here that does not decompose the model into layers -- it
 * simulates entries, activities, task threads and calls directly -- and that is
 * exactly why its mask has to be pinned rather than assumed: it MEASURES more
 * than the table reports. A processor's completion rate and an entry's
 * occupancy are perfectly well defined along a sample path and sit in `TLN` and
 * `QLN`, so nothing inside the engine stops them reaching a column where every
 * other solver prints NaN. The JAR shipped exactly that divergence until
 * 2026-08-21. The mask therefore lives beside the result, in
 * `ldes::engine::ln_defined`, and is checked here against SolverLN's.
 *
 * ArvR is absent from the LDES mask rather than false: the layered table of
 * this port has no such column, and `check_same_mask` skips a key the other
 * side does not carry, which is the difference between "no column" and "a
 * column that disagrees".
 */
TEST_CASE("lqn mask: the native LDES layered engine agrees with SolverLN") {
    const lqn::LqnStruct<double> l = build_multi_solvers();

    ldes::LdesOptions o;
    // The mask is a property of the element kinds, not of the estimates, so a
    // short run pins it: every activity gets its WLN written on the first pass
    // through the result assembly regardless of how long the run was.
    o.samples = 20000;
    o.seed = 23000;
    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(l, o);

    typedef ldes::engine::LnColumn C;
    std::map<std::string, bool> ldes_mask;
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        const std::string& n = l.names[i];
        ldes_mask[n + ".QLen"] = !ldes::engine::ln_defined(l, r, i, C::QLen);
        ldes_mask[n + ".Util"] = !ldes::engine::ln_defined(l, r, i, C::Util);
        ldes_mask[n + ".RespT"] = !ldes::engine::ln_defined(l, r, i, C::RespT);
        ldes_mask[n + ".ResidT"] = !ldes::engine::ln_defined(l, r, i, C::ResidT);
        ldes_mask[n + ".Tput"] = !ldes::engine::ln_defined(l, r, i, C::Tput);
    }

    check_same_mask("LN(MVA)", mask_of(l, solve(l, "mva")), "LDES(native LN)", ldes_mask);

    // spelled out, so the comparison above cannot pass by both sides regressing
    CHECK(ldes_mask.at("P1.QLen"));
    CHECK(ldes_mask.at("P1.Tput"));
    CHECK_FALSE(ldes_mask.at("P1.Util"));
    CHECK(ldes_mask.at("E1.ResidT"));
    CHECK_FALSE(ldes_mask.at("E1.RespT"));
    CHECK_FALSE(ldes_mask.at("T1.ResidT"));
    CHECK_FALSE(ldes_mask.at("A1.ResidT"));
    CHECK(ldes_mask.at("T1.RespT"));
}
