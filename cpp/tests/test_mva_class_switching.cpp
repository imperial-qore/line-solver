/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The closed-population AMVA family on a CLASS-SWITCHING model.
 *
 * THE DEFECT THIS PINS, and why this port is the reference for it. These sixteen
 * algorithms recur on a population vector, and the recursion presumes each entry
 * of it is CONSERVED: the arrival-instant estimate E[Q(N - 1_r)] is only
 * meaningful if removing a customer of r leaves a network of the same shape.
 * Under class switching a job CHANGES CLASS as it moves, so no per-class
 * population is conserved -- the conserved quantity is the CHAIN population.
 * `solver_amva` builds its whole product-form branch out of
 * `sn_get_product_form_chain_params` and deaggregates through
 * `sn_deaggregate_chain_results` at the end, which is exactly why this port and
 * MATLAB were right about the model below while native python and the JAR's two
 * Schmidt arms were not: they handed the kernels the CLASS vector and so solved
 * a different network, returning [2, 0] -- both jobs parked at the delay, none
 * at the queue -- from every one of the sixteen names.
 *
 * These cases therefore assert the reference behaviour rather than a fix: they
 * are what the other two codebases are now measured against, and they fail if
 * this branch ever stops aggregating.
 *
 * WHAT IS ASSERTED: the exact answer comes from SolverCTMC, which solves the
 * generator as written and never goes near these kernels; every approximation
 * must land near it; and they must not all land on the SAME number, because
 * sixteen differently-derived estimators agreeing bit for bit is the signature
 * of the defect rather than of accuracy.
 */

#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** The closed-population algorithms, less 'sqni'; see the note on it below. */
const char* kFamily[] = {"bs",   "aql",   "qsa",  "tay",     "scat", "lcp",
                         "chow", "pamb",  "pami", "pamt",    "clust", "dmlin",
                         "ab",   "schmidt", "schmidt-ext"};

/**
 * Delay -> Queue -> Delay with the class relabelled on each hop.
 *
 * One chain of 2 jobs spread over two classes: C1 is served at the delay and
 * switches to C2 on the way to the queue, C2 switches back on the way out. C2
 * therefore has a population of 0 of its own while the chain holds 2.
 */
qn::Network<double> class_switching(SchedStrategy sched) {
    qn::Network<double> m("cs");
    const std::size_t d = m.add_delay("D");
    const std::size_t q = m.add_queue("Q", sched);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(1.0));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c2, d, q, 1.0);
    P.set(c2, c1, q, d, 1.0);
    m.link(P);
    return m;
}

/** Total mean queue length at station `i`, summed over the classes. */
double qlen_at(const Matrix<double>& QN, std::size_t i) {
    double s = 0.0;
    for (std::size_t r = 0; r < QN.cols(); ++r) s += QN(i, r);
    return s;
}

Matrix<double> mva_qlen(const qn::NetworkStruct<double>& sn, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(sn, opt, init).QN;
}

}  // namespace

TEST_CASE("the closed-population family solves a class-switching model") {
    qn::Network<double> m = class_switching(SchedStrategy::PS);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // The reference, from a solver that never touches these kernels.
    const Matrix<double> exact =
        ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions()).QN;
    const double eq0 = qlen_at(exact, 0), eq1 = qlen_at(exact, 1);
    CHECK(eq0 == doctest::Approx(1.4117647058823530).epsilon(1e-9));
    CHECK(eq1 == doctest::Approx(0.5882352941176471).epsilon(1e-9));

    std::set<long long> distinct;
    for (const char* name : kFamily) {
        CAPTURE(name);
        const Matrix<double> q = mva_qlen(sn, name);
        const double q0 = qlen_at(q, 0), q1 = qlen_at(q, 1);
        // N = 2 jobs are somewhere, whatever approximation is used.
        CHECK(q0 + q1 == doctest::Approx(2.0).epsilon(1e-6));
        // None may be off by the whole queue: [2, 0] is 0.59 out at both
        // stations, which is the entire content of the queue row.
        CHECK(q1 > 0.4);
        CHECK(std::fabs(q0 - eq0) < 0.25);
        CHECK(std::fabs(q1 - eq1) < 0.25);
        distinct.insert(static_cast<long long>(std::llround(q1 * 1e9)));
    }
    // Sixteen differently-derived estimators cannot all land on one number.
    CHECK(distinct.size() > 3);
}

TEST_CASE("the same class-switching model served FCFS") {
    // FCFS is where the extended Schmidt correction is formed, so the two
    // Schmidt arms take a different route through the kernel here.
    qn::Network<double> m = class_switching(SchedStrategy::FCFS);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const Matrix<double> exact =
        ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions()).QN;
    const double eq0 = qlen_at(exact, 0), eq1 = qlen_at(exact, 1);
    const char* names[] = {"schmidt", "schmidt-ext", "ab", "bs"};
    for (const char* name : names) {
        CAPTURE(name);
        const Matrix<double> q = mva_qlen(sn, name);
        CHECK(std::fabs(qlen_at(q, 0) - eq0) < 0.25);
        CHECK(std::fabs(qlen_at(q, 1) - eq1) < 0.25);
    }
}

TEST_CASE("sqni reports its own closed form on the class-switching model") {
    // 'sqni' is held out of the accuracy case above on purpose. Its closed form
    // reports Q = N - X Z off a square-root estimate of X, and on this model that
    // estimate is the saturation value X = 2, which drives the queue term to
    // exactly zero. That is the algorithm and not the defect -- the reference
    // formula in pfqn_sqni.m gives the same 2 and 0 -- so what is pinned here is
    // that it still conserves the population and still runs.
    qn::Network<double> m = class_switching(SchedStrategy::PS);
    const Matrix<double> q = mva_qlen(m.get_struct(), "sqni");
    CHECK(qlen_at(q, 0) + qlen_at(q, 1) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(qlen_at(q, 0) == doctest::Approx(2.0).epsilon(1e-6));
}
