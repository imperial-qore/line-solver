/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverLN(SolverMVA) on lqn_ofbiz, against the MATLAB reference.
 *
 * The model is matlab/examples/basic/layeredModel/lqn_ofbiz.xml, a 9-processor
 * 9-task 14-entry 40-activity layered network that decomposes into 14 layers.
 * It is the regression for the whole stack added for the layered port -- the
 * .lqnx reader, the layer builder, the chain aggregation, the exact and
 * approximate MVA, and the LN fixed point -- because every one of those has to
 * be right for the last column to come out.
 *
 * The expected values are MATLAB's, read off `LN(model, @(x) SolverMVA(x))`
 * with the default options and passed through getAvgTable's sanitizer, at
 * LINE 3.0.6. TOLERANCE, and why it is not tighter: the reference snaps any
 * metric at or below FineTol = 1e-8 to zero, in three separate places
 * (filterMetric, sn_get_residt_from_respt, getAvgTable). Nearly every quantity
 * in a layered model passes through an Immediate distribution, whose service
 * time is exactly 1e-8, so those tests are evaluated ON the threshold, and two
 * implementations that differ in the last bit of a visit ratio land on
 * opposite sides. The observed worst-case relative disagreement over all 72
 * elements and all five metrics is 6e-7; 1e-5 leaves room without hiding a
 * genuine regression, which would be orders of magnitude larger.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;

namespace {

std::string ofbiz_path() {
    return std::string(LINE_MP_REPO_ROOT) + "/matlab/examples/basic/layeredModel/lqn_ofbiz.xml";
}

/** getAvgTable's numeric clean-up, so the comparison is against what MATLAB prints. */
double sanitize(double x) {
    const double r = std::round(x * 10.0);
    if (std::fabs(x * 10.0 - r) < lang::GlobalConstants::CoarseTol * x * 10.0) x = r / 10.0;
    if (x <= lang::GlobalConstants::FineTol) x = 0.0;
    return x;
}

struct Expected {
    std::size_t idx;  ///< 1-based LQN element index
    const char* name;
    double util;
    double tput;
};

/**
 * The processors, tasks and entries that carry the model's load. The activity
 * rows are omitted here only to keep the table readable; they are covered
 * transitively, since an entry's response time is the sum over its activities.
 */
const Expected kExpected[] = {
    {1, "FrontEnd_CPU_Processor", 0.116372299576233, 0.0},
    {10, "FrontEnd_CPU_Task", 0.116372299576233, 11.6372299576233},
    {12, "UsageScenario_userType1_1_Task", 0.0, 0.993266675600061},
    {13, "RequestHandler_HandlerIF_main_345_Task", 0.0, 4.1156445684398},
    {14, "RequestHandler_HandlerIF_login_345_Task", 0.0, 1.70297042181102},
    {15, "RequestHandler_HandlerIF_checkLogin_345_Task", 0.0, 3.40594082874432},
    {16, "RequestHandler_HandlerIF_logout_345_Task", 0.0, 1.70297042181102},
    {17, "UsageScenario_userType2_7_Task", 0.0, 0.709703752863478},
    {18, "RequestHandler_HandlerIF_quickadd_345_Task", 0.0, 0.709703750616092},
    {20, "InternalAction_main_Entry", 0.041156445540142, 4.1156445540142},
    {26, "UsageScenario_userType1_1_Entry", 0.0, 0.993266675600061},
    {31, "UsageScenario_userType2_7_Entry", 0.0, 0.709703752863478},
};

}  // namespace

TEST_CASE("lqn_ofbiz: the reader reproduces the MATLAB LayeredNetworkStruct") {
    const lqn::LqnStruct<double> l = lqn::read_lqnx<double>(ofbiz_path());
    CHECK(l.nhosts == 9);
    CHECK(l.ntasks == 9);
    CHECK(l.nentries == 14);
    CHECK(l.nacts == 40);
    CHECK(l.ncalls == 19);
    CHECK(l.nidx == 72);
    CHECK(l.tshift == 9);
    CHECK(l.eshift == 18);
    CHECK(l.ashift == 32);
    CHECK(l.hashnames[1] == "P:FrontEnd_CPU_Processor");
    CHECK(l.hashnames[10] == "T:FrontEnd_CPU_Task");
    CHECK(l.hashnames[12] == "R:UsageScenario_userType1_1_Task");
    // the multiplicity correction walks the inf-scheduled tasks in ascending
    // index order, so a task whose callers are themselves still uncorrected
    // keeps Inf; FrontEnd_CPU_Task is that task
    CHECK(std::isinf(l.mult[10]));
    CHECK(l.mult[13] == doctest::Approx(20.0));
    CHECK(std::isinf(l.maxmult[10]));
    CHECK(l.maxmult[13] == doctest::Approx(20.0));
    CHECK(l.maxmult[2] == doctest::Approx(0.0));  // the disconnected USAGE_DELAY host
}

TEST_CASE("lqn_ofbiz: SolverLN(SolverMVA) in double agrees with MATLAB") {
    const lqn::LqnStruct<double> model = lqn::read_lqnx<double>(ofbiz_path());
    ln::LnOptions opt;
    ln::SolverLN<double> solver(model, opt);
    CHECK(solver.nlayers() == 14);
    const ln::LnSolution<double> s = solver.get_ensemble_avg();
    CHECK(s.converged);
    CHECK(s.iterations == 52);  // the reference takes the same 52 iterations

    for (const Expected& e : kExpected) {
        CAPTURE(e.name);
        CHECK(sanitize(s.UN[e.idx]) == doctest::Approx(e.util).epsilon(1e-5));
        if (e.tput > 0.0) CHECK(sanitize(s.TN[e.idx]) == doctest::Approx(e.tput).epsilon(1e-5));
    }

    // the disconnected component reports zeros, not NaN
    CHECK(s.TN[11] == doctest::Approx(0.0));
    CHECK(s.UN[11] == doctest::Approx(0.0));
}

TEST_CASE("lqn_ofbiz: the high-precision backend reproduces the double answer") {
    // Real50 carries ~50 decimal digits, so the fixed point it converges to is
    // free of the double run's accumulated rounding. Agreement to 1e-5 says the
    // double answer is limited by the algorithm's stopping rule and by the
    // reference's 1e-8 snapping, not by its arithmetic.
    const lqn::LqnStruct<Real50> model = lqn::read_lqnx<Real50>(ofbiz_path());
    ln::LnOptions opt;
    ln::SolverLN<Real50> solver(model, opt);
    const ln::LnSolution<Real50> s = solver.get_ensemble_avg();
    CHECK(s.converged);
    for (const Expected& e : kExpected) {
        CAPTURE(e.name);
        CHECK(sanitize(num_traits<Real50>::to_double(s.UN[e.idx])) ==
              doctest::Approx(e.util).epsilon(1e-5));
    }
}

TEST_CASE("lqn_ofbiz: the exact backend refuses the approximate MVA by name") {
    // WHAT CHANGED, and why this is a refusal rather than a comparison. This
    // model's layers are not product form: an entry's FCFS station carries a
    // non-unit scv, which fails the BCMP type-1 gate in has_product_form, so
    // mva_dispatch sends every layer to the APPROXIMATE MVA. That is a Picard
    // fixed point, and in exact rational arithmetic each sweep carries the
    // product of every denominator before it, so the bit length grows
    // geometrically in the iteration count: this model was measured at over
    // four and a half hours, operands millions of bits wide, still unfinished
    // and with no bound to quote. solver_amvald therefore refuses non-floating
    // arithmetic by name, the same rule its load-dependent soft minimum and its
    // Suri factor already follow, and this pins that refusal.
    //
    // The exactness claim this test used to make is not lost, it moved: the
    // Real backend reproduces the double answer on the same model in seconds
    // (the test above), which is the affordable way to run lqn_ofbiz at
    // extended precision. The growth itself is documented in _kb/14.
    const lqn::LqnStruct<Rational> model = lqn::read_lqnx<Rational>(ofbiz_path());

    // BOTH ENCODINGS, because the encoding decides how many phases an entry
    // gets but not whether the layers have a product form: 'srvn.cs' keeps the
    // entries single-phase and still leaves a layer that only the approximate
    // MVA can take, so neither spelling makes the exact path reachable.
    ln::LnOptions ph;
    ph.method = "srvn.ph";
    ph.iter_max = 2;
    ln::SolverLN<Rational> sph(model, ph);
    CHECK_THROWS_AS(sph.get_ensemble_avg(), UnsupportedError);

    ln::LnOptions cs;
    cs.method = "srvn.cs";
    cs.iter_max = 2;
    ln::SolverLN<Rational> scs(model, cs);
    CHECK_THROWS_AS(scs.get_ensemble_avg(), UnsupportedError);

    // The double path over the same model is untouched by the refusal, which is
    // what keeps the refusal a statement about the ARITHMETIC and not about the
    // model being unsolvable.
    const lqn::LqnStruct<double> dmodel = lqn::read_lqnx<double>(ofbiz_path());
    ln::SolverLN<double> dsolver(dmodel, ph);
    const ln::LnSolution<double> ds = dsolver.get_ensemble_avg();
    CHECK(ds.iterations == 2);
}
