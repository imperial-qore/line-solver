/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The layered half of the gallery: `gallery_lqn_*`, `gallery_multitier` and
 * `gallery_multitier_storage`.
 *
 * DECLARATION ORDER IS THE MODEL. `lqn_finalize` assigns element indices by
 * kind and, within a kind, by declaration order, exactly as
 * `@@LayeredNetwork/getStruct.m` does, so each factory declares its processors,
 * tasks, entries and activities in the reference file's order and the two
 * codebases index the same element the same way.
 *
 * `ActivityPrecedence.Serial(a, b, c)` is n-1 pairwise sequences in the
 * reference (`ActivityPrecedence.m:79-83`) and is written out as such here,
 * because the builder's `serial` is the pair.
 */

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "gallery.h"

namespace line {
namespace examples {

namespace {

const double INF_MULT = std::numeric_limits<double>::infinity();

/** `Exp(r)`: the reference's constructor takes the RATE, not the mean. */
D E(double rate) { return Exp(rate); }

/** `DiscreteSampler([1/n] * n)`, the uniform popularity of an ItemEntry. */
std::vector<double> uniform_pmf(std::size_t n) {
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

/**
 * The three tiers `gallery_multitier` and `gallery_multitier_storage` share:
 * the client T0, the application T1 with its four entries, and the database T2
 * with its four. Everything the two models differ in is declared by the
 * callers, after this returns.
 */
void multitier_tiers(Lqn& b) {
    b.processor("P0", 1, SchedStrategy::PS);
    b.task("T0", 1, SchedStrategy::REF, "P0");
    b.entry("E0", "T0");

    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 1, SchedStrategy::FCFS, "P1");
    b.entry("E10", "T1");
    b.entry("E11", "T1");
    b.entry("E12", "T1");
    b.entry("E13", "T1");

    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.entry("E20", "T2");
    b.entry("E21", "T2");
    b.entry("E22", "T2");
    b.entry("E23", "T2");
}

/** The client activities A0..A3, identical in both multitier models. */
void multitier_client_acts(Lqn& b) {
    b.activity("A0", E(1.0), "T0");
    b.bound_to("A0", "E0");
    b.sync_call("A0", "E12", 1.0);
    b.activity("A1", E(1.0), "T0");
    b.sync_call("A1", "E10", 1.0);
    b.activity("A2", E(1.0), "T0");
    b.sync_call("A2", "E11", 1.0);
    b.activity("A3", E(1.0), "T0");
    b.sync_call("A3", "E13", 1.0);
}

/** The application activities B0..B7 that the two multitier models share. */
void multitier_app_acts(Lqn& b) {
    b.activity("B0", E(1.0), "T1");
    b.bound_to("B0", "E10");
    b.activity("B1", E(1.0), "T1");
    b.replies_to("B1", "E10");
    b.activity("B2", E(1.0), "T1");
    b.bound_to("B2", "E11");
    b.activity("B3", E(1.0), "T1");
    b.sync_call("B3", "E21", 1.0);
    b.replies_to("B3", "E11");
    b.activity("B4", E(1.0), "T1");
    b.bound_to("B4", "E12");
    b.sync_call("B4", "E20", 1.0);
    b.replies_to("B4", "E12");
    b.activity("B5", E(1.0), "T1");
    b.bound_to("B5", "E13");
    b.activity("B6", E(1.0), "T1");
    b.activity("B7", E(1.0), "T1");
    b.sync_call("B7", "E22", 1.0);
}

/** The precedences of the two shared tiers, minus the or-join of B7a/B7b. */
void multitier_shared_precedences(Lqn& b) {
    b.serial("A0", "A1");
    b.serial("A1", "A2");
    b.serial("A2", "A3");
    b.serial("B0", "B1");
    b.serial("B2", "B3");
    b.serial("B5", "B6");
    b.serial("B6", "B7");
    b.or_fork("B7", {"B7a", "B7b"}, {0.7, 0.3});
    b.serial("C0", "C1");
    b.serial("C3", "C4");
    b.serial("C4", "C5");
}

}  // namespace

// ---------------------------------------------------------------------------
// Basic layered networks
// ---------------------------------------------------------------------------

/** A three-task client/server chain: T1 calls T2, which calls T3 five times. */
Lqn gallery_lqn_basic() {
    Lqn b;
    b.processor("P1", 2, SchedStrategy::PS);
    b.processor("P2", 3, SchedStrategy::PS);

    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", E(1.0 / 2.0));
    b.task("T2", 50, SchedStrategy::FCFS, "P1");
    b.think_time("T2", E(1.0 / 3.0));
    b.task("T3", 25, SchedStrategy::FCFS, "P2");
    b.think_time("T3", E(1.0 / 4.0));

    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");

    b.activity("AS1", E(10.0), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", E(20.0), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5.0);
    b.replies_to("AS2", "E2");
    b.activity("AS3", E(50.0), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b;
}

/**
 * `LayeredNetworkGenerator().generate(1, 2, 4, 2)` at seed 23000.
 *
 * Every range of the default generator is degenerate, so the population, the
 * think time, both multiplicities, every host demand and every call mean come
 * out as exactly 1 and only the TOPOLOGY is drawn. The realized topology,
 * read off the reference's `getStruct()` (parent, callpair), is: tasks 1-3 on
 * processor_1 and task_4 on processor_2; the client calls entry_1 and entry_2
 * (the two first-level tasks); activity_1 calls entry_4 and activity_2 calls
 * entry_3.
 *
 * The client's think time is carried as the HOST DEMAND of `c_activity_1`,
 * which is where `_create_clients` puts it; the reference task's own think time
 * stays Immediate.
 */
Lqn gallery_lqn_random() {
    Lqn b;
    b.processor("c_processor_1", INF_MULT, SchedStrategy::INF);
    b.processor("processor_1", 1, SchedStrategy::PS);
    b.processor("processor_2", 1, SchedStrategy::PS);

    b.task("c_task_1", 1, SchedStrategy::REF, "c_processor_1");
    b.task("task_1", 1, SchedStrategy::FCFS, "processor_1");
    b.task("task_2", 1, SchedStrategy::FCFS, "processor_1");
    b.task("task_3", 1, SchedStrategy::FCFS, "processor_1");
    b.task("task_4", 1, SchedStrategy::FCFS, "processor_2");

    b.entry("c_entry_1", "c_task_1");
    for (std::size_t t = 1; t <= 4; ++t)
        b.entry("entry_" + std::to_string(t), "task_" + std::to_string(t));

    b.activity("c_activity_1", E(1.0), "c_task_1");
    b.bound_to("c_activity_1", "c_entry_1");
    for (std::size_t t = 1; t <= 4; ++t) {
        const std::string a = "activity_" + std::to_string(t);
        const std::string e = "entry_" + std::to_string(t);
        b.activity(a, E(1.0), "task_" + std::to_string(t));
        b.bound_to(a, e);
        b.replies_to(a, e);
    }

    b.sync_call("c_activity_1", "entry_1", 1.0);
    b.sync_call("c_activity_1", "entry_2", 1.0);
    b.sync_call("activity_1", "entry_4", 1.0);
    b.sync_call("activity_2", "entry_3", 1.0);
    return b;
}

/**
 * A loop on T1, an and-fork/join on T2 and an or-fork/join on T3.
 *
 * T3 carries BOTH an infinite-server discipline and a think time of 10, which
 * the .lqnx writer cannot express, so this model exists only through the
 * builder.
 */
Lqn gallery_lqn_workflows() {
    Lqn b;
    b.processor("P1", INF_MULT, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", Immediate());
    b.entry("Entry", "T1");

    b.processor("P2", INF_MULT, SchedStrategy::INF);
    b.task("T2", INF_MULT, SchedStrategy::INF, "P2");
    b.think_time("T2", Immediate());
    b.entry("E2", "T2");

    b.processor("P3", 5, SchedStrategy::PS);
    b.task("T3", INF_MULT, SchedStrategy::INF, "P3");
    b.think_time("T3", D::exp_mean(10.0));
    b.entry("E3", "T3");

    b.activity("A1", D::exp_mean(1.0), "T1");
    b.bound_to("A1", "Entry");
    b.activity("A2", D::exp_mean(2.0), "T1");
    b.activity("A3", D::exp_mean(3.0), "T1");
    b.sync_call("A3", "E2", 1.0);

    b.activity("B1", D::exp_mean(0.1), "T2");
    b.bound_to("B1", "E2");
    b.activity("B2", D::exp_mean(0.2), "T2");
    b.activity("B3", D::exp_mean(0.3), "T2");
    b.activity("B4", D::exp_mean(0.4), "T2");
    b.activity("B5", D::exp_mean(0.5), "T2");
    b.activity("B6", D::exp_mean(0.6), "T2");
    b.sync_call("B6", "E3", 1.0);
    b.replies_to("B6", "E2");

    b.activity("C1", D::exp_mean(0.1), "T3");
    b.bound_to("C1", "E3");
    b.activity("C2", D::exp_mean(0.2), "T3");
    b.activity("C3", D::exp_mean(0.3), "T3");
    b.activity("C4", D::exp_mean(0.4), "T3");
    b.activity("C5", D::exp_mean(0.5), "T3");
    b.replies_to("C5", "E3");

    // Loop(A1, {A2, A3}, 3) runs A2 three times and exits to A3.
    b.loop("A1", {"A2"}, "A3", 3.0);
    b.serial("B4", "B5");
    b.and_fork("B1", {"B2", "B3", "B4"});
    b.and_join({"B2", "B3", "B5"}, "B6");
    b.or_fork("C1", {"C2", "C3", "C4"}, {0.3, 0.3, 0.4});
    b.or_join({"C2", "C3", "C4"}, "C5");
    return b;
}

// ---------------------------------------------------------------------------
// The multitier J2EE models
// ---------------------------------------------------------------------------

/** The 3-tier reference LQN: client, application server and database. */
Lqn gallery_multitier() {
    Lqn b;
    multitier_tiers(b);
    multitier_client_acts(b);
    multitier_app_acts(b);
    b.activity("B7a", E(1.0), "T1");
    b.activity("B7b", E(1.0), "T1");
    b.sync_call("B7b", "E23", 1.0);
    b.activity("B8", E(1.0), "T1");
    b.replies_to("B8", "E13");

    b.activity("C0", E(1.0), "T2");
    b.bound_to("C0", "E20");
    b.activity("C1", E(1.0), "T2");
    b.replies_to("C1", "E20");
    b.activity("C2", E(1.0), "T2");
    b.bound_to("C2", "E21");
    b.replies_to("C2", "E21");
    b.activity("C3", E(1.0), "T2");
    b.bound_to("C3", "E22");
    b.activity("C4", E(1.0), "T2");
    b.activity("C5", E(1.0), "T2");
    b.replies_to("C5", "E22");
    b.activity("C6", E(1.0), "T2");
    b.bound_to("C6", "E23");
    b.replies_to("C6", "E23");

    multitier_shared_precedences(b);
    b.or_join({"B7a", "B7b"}, "B8");
    return b;
}

/**
 * The same four tiers with a dedicated cache layer.
 *
 * The database entries E20 and E21 look the item up in the CacheTask T3, whose
 * read activity D0 branches to a fast hit (D1a) or a slow miss (D1b) through
 * POST_CACHE. There is no B8 here: B7a and B7b reply to E13 themselves, so the
 * or-join of `gallery_multitier` is absent rather than replaced.
 */
Lqn gallery_multitier_storage() {
    const std::size_t totalitems = 10;
    const int cachecapacity = 2;

    Lqn b;
    multitier_tiers(b);
    b.processor("P3", 1, SchedStrategy::PS);
    b.cache_task("T3", 1, SchedStrategy::FCFS, "P3", totalitems, {cachecapacity},
                 lang::ReplacementStrategy::LRU);
    b.item_entry("E3", "T3", totalitems, uniform_pmf(totalitems));

    multitier_client_acts(b);
    multitier_app_acts(b);
    b.activity("B7a", E(1.0), "T1");
    b.replies_to("B7a", "E13");
    b.activity("B7b", E(1.0), "T1");
    b.sync_call("B7b", "E23", 1.0);
    b.replies_to("B7b", "E13");

    b.activity("C0", E(1.0), "T2");
    b.bound_to("C0", "E20");
    b.activity("C1", E(1.0), "T2");
    b.sync_call("C1", "E3", 1.0);
    b.replies_to("C1", "E20");
    b.activity("C2", E(1.0), "T2");
    b.bound_to("C2", "E21");
    b.sync_call("C2", "E3", 1.0);
    b.replies_to("C2", "E21");
    b.activity("C3", E(1.0), "T2");
    b.bound_to("C3", "E22");
    b.activity("C4", E(1.0), "T2");
    b.activity("C5", E(1.0), "T2");
    b.replies_to("C5", "E22");
    b.activity("C6", E(1.0), "T2");
    b.bound_to("C6", "E23");
    b.replies_to("C6", "E23");

    b.activity("D0", Immediate(), "T3");
    b.bound_to("D0", "E3");
    b.activity("D1a", E(1.0), "T3");
    b.replies_to("D1a", "E3");
    b.activity("D1b", E(0.5), "T3");
    b.replies_to("D1b", "E3");

    multitier_shared_precedences(b);
    b.cache_access("D0", "D1a", "D1b");
    return b;
}

}  // namespace examples
}  // namespace line
