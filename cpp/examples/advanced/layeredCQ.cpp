/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `advanced/layeredCQ`: layered cache-queueing models.
 *
 * A CacheTask is an ordinary task whose entry is an ITEM ENTRY -- a read over a
 * population of items with a popularity pmf -- and whose read activity is
 * followed by a POST_CACHE precedence splitting the job onto a hit branch and a
 * miss branch. The hit ratio is not declared: it is what the cache replacement
 * policy produces at the offered read rate, so the whole model turns on it.
 *
 * The Cache NODE lives in the cache task's HOST layer rather than in the task's
 * own layer, because the lookup is what the task does on its processor; that is
 * `buildLayersRecursive`'s `iscachelayer` host-layer test and it is why the
 * three-host example has the hit and the miss queueing at the same server the
 * read arrived at.
 *
 * `DiscreteSampler([1/n]*n)` is a discrete distribution the reference reads the
 * pmf out of; this port has no discrete-distribution type, so `item_entry`
 * takes the pmf directly, which is the same object with one fewer wrapper.
 */

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

namespace line {
namespace examples {

namespace {

using Lqn = lqn::LqnStruct<double>;
using lang::ReplacementStrategy;

/** `ln_element_kind` of the CLI: what an LQN element is, for the NodeType column. */
const char* lqn_kind(const Lqn& l, std::size_t i) {
    switch (l.type[i]) {
        case lang::LqnElement::HOST: return "Processor";
        case lang::LqnElement::TASK: return l.isref[i] ? "RefTask" : "Task";
        case lang::LqnElement::ENTRY: return "Entry";
        default: return "Activity";
    }
}

/** The CLI's `ln_sanitize`: snap to a tenth, and report a negligible value as zero. */
double ln_round(double x) {
    const double r = std::round(x * 10.0);
    if (std::fabs(x * 10.0 - r) < lang::GlobalConstants::CoarseTol * x * 10.0) x = r / 10.0;
    if (x <= lang::GlobalConstants::FineTol) x = 0.0;
    return x;
}

/** The golden's key for the layer engine `layer_solver` names. */
const char* lcq_layer_golden_key(const std::string& layer_solver) {
    if (layer_solver == "fluid") return "FLD";
    if (layer_solver == "nc") return "NC";
    if (layer_solver == "ssa") return "SSA";
    return "MVA";
}

/** `SolverLN.getAvgTable()`: one row per LQN element, undefined metrics as NaN. */
void print_ln_avg(const Lqn& l, const ln::LnSolution<double>& sol,
                  const std::string& solver_key) {
    std::printf("%-24s %-10s %12s %12s %12s %12s %12s\n", "Node", "NodeType", "QLen", "Util",
                "RespT", "ResidT", "Tput");
    // The recorded row is the printed row: the golden was generated from this
    // table, so it holds `ln_round`'s value and not the raw solution's.
    std::vector<LnRow> recorded;
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        char q[24], u[24], rr[24], w[24], t[24];
        const std::vector<double>* vals[5] = {&sol.QN, &sol.UN, &sol.RN, &sol.WN, &sol.TN};
        const std::vector<bool>* defs[5] = {&sol.defined_Q, &sol.defined_U, &sol.defined_R,
                                            &sol.defined_W, &sol.defined_T};
        char* bufs[5] = {q, u, rr, w, t};
        double cell[5];
        for (std::size_t k = 0; k < 5; ++k) {
            cell[k] = (*defs[k])[i] ? ln_round((*vals[k])[i])
                                    : std::numeric_limits<double>::quiet_NaN();
            if (!(*defs[k])[i]) std::snprintf(bufs[k], 24, "%12s", "NaN");
            else std::snprintf(bufs[k], 24, "%12.6g", cell[k]);
        }
        std::printf("%-24s %-10s %s %s %s %s %s\n", l.names[i].c_str(), lqn_kind(l, i), q, u, rr,
                    w, t);
        LnRow row;
        row.name = l.names[i];
        row.q = cell[0];
        row.u = cell[1];
        row.r = cell[2];
        row.w = cell[3];
        row.t = cell[4];
        recorded.push_back(row);
    }
    record_ln(solver_key, recorded);
}

/** Run the ensemble with the named layer engine and print its table. */
void run_ln(const Lqn& l, const std::string& layer_solver) {
    ln::LnOptions opt;
    opt.layer_solver = layer_solver;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    std::printf("SolverLN(layers=%s) layers=%zu iterations=%d converged=%d\n",
                layer_solver.c_str(), s.nlayers(), sol.iterations, int(sol.converged));
    print_ln_avg(l, sol, ln_solver_key(lcq_layer_golden_key(layer_solver), opt.method));
}

/** A uniform read popularity over `n` items, `DiscreteSampler([1/n]*n)`. */
std::vector<double> uniform_access(std::size_t n) {
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

}  // namespace

// ---------------------------------------------------------------------------
// lcq_singlehost
// ---------------------------------------------------------------------------

/**
 * A client whose single activity is a SYNCHRONOUS read of a cache holding 2 of
 * 4 items under random replacement; the read splits into a fast hit and a slow
 * miss on the same host.
 */
void lcq_singlehost() {
    const std::size_t totalitems = 4;

    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");

    b.processor("PC", 1, SchedStrategy::PS);
    b.cache_task("C2", 1, SchedStrategy::FCFS, "PC", totalitems, {2}, ReplacementStrategy::RR);
    b.item_entry("I2", "C2", totalitems, uniform_access(totalitems));

    b.activity("A1", Immediate(), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "I2", 1.0);

    b.activity("AC2", Immediate(), "C2");
    b.bound_to("AC2", "I2");
    b.activity("AC2h", Exp(1.0), "C2");
    b.activity("AC2m", Exp(0.5), "C2");
    b.cache_access("AC2", "AC2h", "AC2m");
    b.replies_to("AC2h", "I2");
    b.replies_to("AC2m", "I2");

    const Lqn l = b.build();
    section("LN(MVA)");
    run_ln(l, "mva");
}

// ---------------------------------------------------------------------------
// lcq_async_cache
// ---------------------------------------------------------------------------

/**
 * The same shape with an ASYNCHRONOUS read: the client fires the cache call and
 * continues without waiting for the reply, so it never blocks on the cache.
 * The POST_CACHE hit/miss split is unchanged, and the replacement policy here
 * is LRU rather than random.
 */
void lcq_async_cache() {
    const std::size_t totalitems = 4;

    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");

    b.processor("PC", 1, SchedStrategy::PS);
    b.cache_task("C2", 1, SchedStrategy::FCFS, "PC", totalitems, {2}, ReplacementStrategy::LRU);
    b.item_entry("I2", "C2", totalitems, uniform_access(totalitems));

    b.activity("A1", Immediate(), "T1");
    b.bound_to("A1", "E1");
    b.async_call("A1", "I2", 1.0);

    b.activity("AC2", Immediate(), "C2");
    b.bound_to("AC2", "I2");
    b.activity("AC2h", Exp(1.0), "C2");
    b.activity("AC2m", Exp(0.5), "C2");
    b.cache_access("AC2", "AC2h", "AC2m");
    b.replies_to("AC2h", "I2");
    b.replies_to("AC2m", "I2");

    const Lqn l = b.build();
    note("=== Solving Async Cache Model ===");
    section("LN(MVA)");
    run_ln(l, "mva");

    note("\n=== Compare with Synchronous Version ===");
    note("To compare async vs sync cache access, run lcq_singlehost.");
    note("");
    note("Expected differences:");
    note("- Async: Lower client response time (no blocking)");
    note("- Async: Higher client throughput (no cache bottleneck)");
    note("- Similar cache hit/miss ratios (same cache logic)");
}

// ---------------------------------------------------------------------------
// lcq_threehosts
// ---------------------------------------------------------------------------

/**
 * A two-level cache in front of a back-end server: the miss branch makes its
 * own synchronous call to a second task, so a miss pays the server's time on
 * top of its own and the hit ratio decides how often that happens.
 */
void lcq_threehosts() {
    const std::size_t totalitems = 4;
    const double nusers = 1.0, ntokens = 1.0;

    lqn::LqnBuilder<double> b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", nusers, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");

    // Two cache lists of one item each: the multi-level cache of the reference.
    b.processor("Pc", 1, SchedStrategy::PS);
    b.cache_task("CT", ntokens, SchedStrategy::FCFS, "Pc", totalitems, {1, 1},
                 ReplacementStrategy::RR);
    b.item_entry("IE", "CT", totalitems, uniform_access(totalitems));

    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E2", "T2");
    b.activity("A2", Exp(5.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");

    b.activity("A1", Immediate(), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "IE", 1.0);

    b.activity("Ac", Immediate(), "CT");
    b.bound_to("Ac", "IE");
    b.activity("Ac_hit", Exp(1.0), "CT");
    b.activity("Ac_miss", Exp(0.5), "CT");
    b.cache_access("Ac", "Ac_hit", "Ac_miss");
    b.sync_call("Ac_miss", "E2", 1.0);
    b.replies_to("Ac_hit", "IE");
    b.replies_to("Ac_miss", "IE");

    const Lqn l = b.build();
    section("LN(NC)");
    run_ln(l, "nc");
    section("LN(MVA)");
    run_ln(l, "mva");
}

LINE_EXAMPLE("advanced/layeredCQ", lcq_singlehost);
LINE_EXAMPLE("advanced/layeredCQ", lcq_async_cache);
LINE_EXAMPLE("advanced/layeredCQ", lcq_threehosts);

}  // namespace examples
}  // namespace line
