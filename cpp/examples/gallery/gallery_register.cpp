/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The runnable form of every gallery entry.
 *
 * A gallery file's `__main__` builds the model and prints its name, and nothing
 * else: the entry is a FACTORY, and solving it is what `test_gallery` does. So
 * each runner here is that `__main__` -- one call and one line -- and the
 * factories themselves stay in `gallery_open.cpp`, `gallery_closed.cpp` and
 * `gallery_lqn.cpp` where a test can reuse them.
 *
 * THE RUNNER MUST CARRY THE REFERENCE'S NAME, because `LINE_EXAMPLE` registers
 * a function under its own identifier and `line-examples gallery_mm1` has to
 * find `gallery_mm1.py`'s model. A runner therefore cannot live beside the
 * factory of the same name, and sits in a nested namespace instead; the macro
 * writes the pair so the 70 registrations read as 70 lines rather than 210.
 */

#include <cstdio>
#include <string>

#include "gallery.h"

namespace line {
namespace examples {
namespace gallery_runners {

namespace {

/** The `print(f"Model: {model.getName()}")` of a Network entry. */
void show(Net m) { std::printf("Model: %s\n", m.raw_struct().name.c_str()); }

/** The same for an Environment, whose name is on the environment itself. */
void show_env(const env::Environment<double>& e) {
    std::printf("Model: %s\n", e.name().c_str());
}

/**
 * The same for a layered network. `LqnBuilder` carries no model name -- the
 * reference's `LayeredNetwork('...')` label survives nowhere in `LqnStruct` --
 * so the name is passed in as the literal the reference file declares, and the
 * builder is constructed and discarded exactly as the reference's `__main__`
 * does.
 */
void show_lqn(const Lqn&, const std::string& nm) {
    std::printf("Model: %s\n", nm.c_str());
}

}  // namespace

#define GALLERY_RUNNER(fn)                       \
    void fn() { show(::line::examples::fn()); }  \
    LINE_EXAMPLE("gallery", fn)

#define GALLERY_RUNNER_ENV(fn)                       \
    void fn() { show_env(::line::examples::fn()); }  \
    LINE_EXAMPLE("gallery", fn)

#define GALLERY_RUNNER_LQN(fn, model_name)                       \
    void fn() { show_lqn(::line::examples::fn(), model_name); }  \
    LINE_EXAMPLE("gallery", fn)


// --- Open queues, single class ---
GALLERY_RUNNER(gallery_mm1);
GALLERY_RUNNER(gallery_mm1_ps);
GALLERY_RUNNER(gallery_mmk);
GALLERY_RUNNER(gallery_mm1k);
GALLERY_RUNNER(gallery_mdk);
GALLERY_RUNNER(gallery_merl1);
GALLERY_RUNNER(gallery_merlk);
GALLERY_RUNNER(gallery_mhyp1);
GALLERY_RUNNER(gallery_mhypk);
GALLERY_RUNNER(gallery_mpar1);
GALLERY_RUNNER(gallery_erlm1);
GALLERY_RUNNER(gallery_erlm1_ps);
GALLERY_RUNNER(gallery_erlm1ps);
GALLERY_RUNNER(gallery_erldk);
GALLERY_RUNNER(gallery_hyperlk);
GALLERY_RUNNER(gallery_hypm1);
GALLERY_RUNNER(gallery_detm1);
GALLERY_RUNNER(gallery_dm1);
GALLERY_RUNNER(gallery_gamm1);
GALLERY_RUNNER(gallery_parm1);
GALLERY_RUNNER(gallery_um1);
GALLERY_RUNNER(gallery_aphm1);
GALLERY_RUNNER(gallery_coxm1);
GALLERY_RUNNER(gallery_mapm1);
GALLERY_RUNNER(gallery_mapmk);
GALLERY_RUNNER(gallery_mmap1);
GALLERY_RUNNER(gallery_mmapk);
GALLERY_RUNNER(gallery_replayerm1);

// --- Open queues, multiclass, feedback, reentrant, tandem ---
GALLERY_RUNNER(gallery_mm1_multiclass);
GALLERY_RUNNER(gallery_mm1_ps_multiclass);
GALLERY_RUNNER(gallery_mm1_prio);
GALLERY_RUNNER(gallery_mm1_feedback);
GALLERY_RUNNER(gallery_mm1_ps_feedback);
GALLERY_RUNNER(gallery_mm1_reentrant);
GALLERY_RUNNER(gallery_mm1_ps_reentrant);
GALLERY_RUNNER(gallery_mm1_linear);
GALLERY_RUNNER(gallery_mm1_tandem);
GALLERY_RUNNER(gallery_mm1_tandem_multiclass);
GALLERY_RUNNER(gallery_merl1_linear);
GALLERY_RUNNER(gallery_merl1_tandem);
GALLERY_RUNNER(gallery_merl1_reentrant);
GALLERY_RUNNER(gallery_mhyp1_linear);
GALLERY_RUNNER(gallery_mhyp1_tandem);
GALLERY_RUNNER(gallery_mhyp1_reentrant);
GALLERY_RUNNER(gallery_hyphyp1_linear);
GALLERY_RUNNER(gallery_hyphyp1_tandem);
GALLERY_RUNNER(gallery_hyphyp1_reentrant);
GALLERY_RUNNER(gallery_hyperl1_feedback);
GALLERY_RUNNER(gallery_hyperl1_reentrant);
GALLERY_RUNNER(gallery_hypm1_reentrant);
GALLERY_RUNNER(gallery_erlerl1);
GALLERY_RUNNER(gallery_erlerl1_reentrant);
GALLERY_RUNNER(gallery_erlm1_reentrant);
GALLERY_RUNNER(gallery_mmap1_multiclass);
GALLERY_RUNNER(gallery_lukumar_reentrant);

// --- Closed networks ---
GALLERY_RUNNER(gallery_cqn);
GALLERY_RUNNER(gallery_cqn_multiclass);
GALLERY_RUNNER(gallery_repairmen);
GALLERY_RUNNER(gallery_qn_random);

// --- Caches, fork-join, finite capacity ---
GALLERY_RUNNER(gallery_cache_lru);
GALLERY_RUNNER(gallery_cache_routing);
GALLERY_RUNNER(gallery_fj_closed);
GALLERY_RUNNER(gallery_fj_open);
GALLERY_RUNNER(gallery_fj_quorum);
GALLERY_RUNNER(gallery_fcr);

// --- Random environment ---
GALLERY_RUNNER_ENV(gallery_renv_breakdown);

// --- Layered networks ---
GALLERY_RUNNER_LQN(gallery_lqn_basic, "LQN-Basic");
GALLERY_RUNNER_LQN(gallery_lqn_random, "lnw");
GALLERY_RUNNER_LQN(gallery_lqn_workflows, "LQN-Workflows");
GALLERY_RUNNER_LQN(gallery_multitier, "testLQN3");
GALLERY_RUNNER_LQN(gallery_multitier_storage, "testLQN3_Cache");

#undef GALLERY_RUNNER
#undef GALLERY_RUNNER_ENV
#undef GALLERY_RUNNER_LQN

}  // namespace gallery_runners
}  // namespace examples
}  // namespace line
