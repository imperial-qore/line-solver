/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/test_gallery/`: every gallery model that MVA can solve,
 * solved by MVA and printed as its AvgTable.
 *
 * Each reference file is four lines -- import the factory, `MVA(model)`,
 * `getAvgTable()`, print -- so each runner here is the same four, and the model
 * itself comes from `gallery.h` rather than being transcribed a second time.
 * That is the point of the split: if a factory drifts, these move with it.
 *
 * THE 57 ENTRIES ARE THE REFERENCE'S 57, no more. The gallery also carries
 * caches, fork-join, a finite capacity region, a random environment and five
 * layered networks, and `test_gallery` covers none of them -- they need a
 * solver other than plain MVA, and inventing a runner for them here would
 * report numbers the reference never produced.
 */

#include <cstdio>

#include "gallery/gallery.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

/** `MVA(model).getAvgTable()` at the default options, printed as it prints. */
void solve_and_print(Net m) {
    std::printf("Model: %s\n", m.raw_struct().name.c_str());
    const Sn& sn = m.get_struct();
    mva::MvaOptions opt;
    const Matrix<double> init_sol;
    // NO BANNER, BUT STILL ATTRIBUTED. The reference is four lines that print
    // the table and nothing else, so `section("MVA")` would put a line in this
    // twin's output that its reference never had. `attribute` makes the same
    // declaration to the recorder alone -- without it every one of these 57
    // tables is recorded under no solver, dropped, and read downstream as
    // "solver MVA missing from the recorded results".
    attribute("MVA");
    print_avg(sn, mva::solver_mva_run_analyzer(sn, opt, init_sol));
}

}  // namespace

#define TEST_GALLERY(fn)                                           \
    void test_gallery_##fn() { solve_and_print(gallery_##fn()); }  \
    LINE_EXAMPLE("test_gallery", test_gallery_##fn)

TEST_GALLERY(aphm1);
TEST_GALLERY(coxm1);
TEST_GALLERY(cqn);
TEST_GALLERY(cqn_multiclass);
TEST_GALLERY(detm1);
TEST_GALLERY(dm1);
TEST_GALLERY(erldk);
TEST_GALLERY(erlerl1);
TEST_GALLERY(erlerl1_reentrant);
TEST_GALLERY(erlm1);
TEST_GALLERY(erlm1_ps);
TEST_GALLERY(erlm1_reentrant);
TEST_GALLERY(erlm1ps);
TEST_GALLERY(gamm1);
TEST_GALLERY(hyperl1_feedback);
TEST_GALLERY(hyperl1_reentrant);
TEST_GALLERY(hyperlk);
TEST_GALLERY(hyphyp1_linear);
TEST_GALLERY(hyphyp1_reentrant);
TEST_GALLERY(hyphyp1_tandem);
TEST_GALLERY(hypm1);
TEST_GALLERY(hypm1_reentrant);
TEST_GALLERY(lukumar_reentrant);
TEST_GALLERY(mapm1);
TEST_GALLERY(mapmk);
TEST_GALLERY(mdk);
TEST_GALLERY(merl1);
TEST_GALLERY(merl1_linear);
TEST_GALLERY(merl1_reentrant);
TEST_GALLERY(merl1_tandem);
TEST_GALLERY(merlk);
TEST_GALLERY(mhyp1);
TEST_GALLERY(mhyp1_linear);
TEST_GALLERY(mhyp1_reentrant);
TEST_GALLERY(mhyp1_tandem);
TEST_GALLERY(mhypk);
TEST_GALLERY(mm1);
TEST_GALLERY(mm1_feedback);
TEST_GALLERY(mm1_linear);
TEST_GALLERY(mm1_multiclass);
TEST_GALLERY(mm1_prio);
TEST_GALLERY(mm1_ps);
TEST_GALLERY(mm1_ps_feedback);
TEST_GALLERY(mm1_ps_multiclass);
TEST_GALLERY(mm1_ps_reentrant);
TEST_GALLERY(mm1_reentrant);
TEST_GALLERY(mm1_tandem);
TEST_GALLERY(mm1_tandem_multiclass);
TEST_GALLERY(mmap1);
TEST_GALLERY(mmap1_multiclass);
TEST_GALLERY(mmapk);
TEST_GALLERY(mmk);
TEST_GALLERY(mpar1);
TEST_GALLERY(parm1);
TEST_GALLERY(repairmen);
TEST_GALLERY(replayerm1);
TEST_GALLERY(um1);

#undef TEST_GALLERY

}  // namespace examples
}  // namespace line
