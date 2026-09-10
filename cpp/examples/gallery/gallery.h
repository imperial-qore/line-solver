/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_EXAMPLES_GALLERY_H
#define LINE_EXAMPLES_GALLERY_H

/**
 * The gallery: `matlab/examples/gallery/*.m` and `python/examples/gallery/*.py`.
 *
 * A gallery entry is a MODEL FACTORY and nothing else -- the MATLAB file is one
 * `function model = gallery_x`, and the Python file's `__main__` prints the
 * model name. They are declared here so the `test_gallery` runners, and any
 * future test, consume the same models the reference does rather than a second
 * transcription of them.
 *
 * DEFAULT ARGUMENTS ARE THE REFERENCE'S. `gallery_mmk(k = 2)` keeps its
 * default, so `gallery_mmk()` is the model `test_gallery_mmk.py` solves.
 */

#include <string>
#include <vector>

#include "examples_common.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/qn/environment.h"

namespace line {
namespace examples {

using Lqn = lqn::LqnBuilder<double>;

// --- Open queues, single class ---------------------------------------------
Net gallery_mm1();
Net gallery_mm1_ps();
Net gallery_mmk(double k = 2);
Net gallery_mm1k(double K = 3);
Net gallery_mdk(double k = 2);
Net gallery_merl1();
Net gallery_merlk(double k = 2);
Net gallery_mhyp1();
Net gallery_mhypk(double k = 2);
Net gallery_mpar1();
Net gallery_erlm1();
Net gallery_erlm1_ps();
Net gallery_erlm1ps();
Net gallery_erldk(double k = 2);
Net gallery_hyperlk(double k = 2);
Net gallery_hypm1();
Net gallery_detm1();
Net gallery_dm1();
Net gallery_gamm1();
Net gallery_parm1();
Net gallery_um1();
Net gallery_aphm1();
Net gallery_coxm1();
Net gallery_mapm1();
Net gallery_mapmk(double k = 2);
Net gallery_mmap1();
Net gallery_mmapk(double k = 2);
Net gallery_replayerm1();

// --- Open queues, multiclass, feedback, reentrant, tandem ------------------
Net gallery_mm1_multiclass();
Net gallery_mm1_ps_multiclass();
Net gallery_mm1_prio();
Net gallery_mm1_feedback();
Net gallery_mm1_ps_feedback();
Net gallery_mm1_reentrant();
Net gallery_mm1_ps_reentrant();
Net gallery_mm1_linear(std::size_t n = 2, double umax = 0.9);
Net gallery_mm1_tandem();
Net gallery_mm1_tandem_multiclass();
Net gallery_merl1_linear(std::size_t n = 2, double umax = 0.9);
Net gallery_merl1_tandem();
Net gallery_merl1_reentrant();
Net gallery_mhyp1_linear(std::size_t n = 2, double umax = 0.9);
Net gallery_mhyp1_tandem();
Net gallery_mhyp1_reentrant();
Net gallery_hyphyp1_linear(std::size_t n = 2, double umax = 0.9);
Net gallery_hyphyp1_tandem();
Net gallery_hyphyp1_reentrant();
Net gallery_hyperl1_feedback();
Net gallery_hyperl1_reentrant();
Net gallery_hypm1_reentrant();
Net gallery_erlerl1(std::size_t n = 5);
Net gallery_erlerl1_reentrant();
Net gallery_erlm1_reentrant();
Net gallery_mmap1_multiclass();
Net gallery_lukumar_reentrant(const std::string& sched = "FCFS");

// --- Closed networks --------------------------------------------------------
// The three drawn factories keep the reference's `seed` argument: their
// parameters come from `random.random()`, which `py_random.h` replays, so the
// seed is part of the model and not a detail of how it was built.
Net gallery_cqn(std::size_t M = 2, bool use_delay = false, unsigned long long seed = 23000);
Net gallery_cqn_multiclass(std::size_t m = 1, std::size_t r = 2, bool wantdelay = true,
                           unsigned long long seed = 23000);
Net gallery_repairmen(double nservers = 1, unsigned long long seed = 2300);
Net gallery_qn_random();

// --- Caches, fork-join, finite capacity ------------------------------------
Net gallery_cache_lru();
Net gallery_cache_routing();
Net gallery_fj_closed();
Net gallery_fj_open();
Net gallery_fj_quorum();
Net gallery_fcr(double K = 3);

// --- Random environment -----------------------------------------------------
env::Environment<double> gallery_renv_breakdown();

// --- Layered networks -------------------------------------------------------
Lqn gallery_lqn_basic();
Lqn gallery_lqn_random();
Lqn gallery_lqn_workflows();
Lqn gallery_multitier();
Lqn gallery_multitier_storage();

}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_GALLERY_H
