/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The AMVA-LD peak utilization post-pass (solver_amvald.m:255-291).
 *
 * A station with limited class dependence reports utilization as T*S/peak
 * from the declared peak rate scaling, the same convention the NC, CTMC and
 * SSA paths already apply. THE ORACLE IS MATLAB on the same model
 * (SolverMVA default -> egflin -> amvald): QN=1.737033963, TN=2.262965855,
 * UN=0.7543219517, run 2026-08-05. The CTMC bound is looser because amvald's
 * class-dependence term is a fixed-point approximation, not the exact
 * recursion. Before the dispatch fixes this model ran through the
 * PRODUCT-FORM kernels, which silently ignored the class dependence
 * (QN 2.5719, the unscaled recursion, 43% off the CTMC).
 */

#include <algorithm>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace mva = line::mva;
using line::Matrix;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {
/** Think -> PS Queue -> Think, one closed class; the queue speeds up with
 * occupancy as beta(n) = min(n, 2), declared peak 2. */
qn::Network<double> cd_cqn(double njobs) {
    qn::Network<double> m("cd_cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(1.5));
    m.set_class_dependence(q, [](const std::vector<double>& n) {
        return std::vector<double>(1, std::min(n[0], 2.0));
    }, std::vector<double>(1, 2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}
}  // namespace

TEST_CASE("amvald reports class-dependent utilization as T*S/peak, matching the CTMC") {
    qn::Network<double> m = cd_cqn(4.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    mva::MvaOptions opt;
    Matrix<double> init;
    const mva::MvaSolution<double> r = mva::solver_mva_analyzer(sn, opt, init);

    const ctmc::CtmcOptions copt;
    const ctmc::CtmcSolution<double> e = ctmc::solver_ctmc_analyzer(sn, copt);

    // station index 1 is the scaled PS queue (0 is the Think delay)
    const double bmax = 2.0, rate = 1.5;
    CHECK(r.U(1, 0) == doctest::Approx(r.Tp(1, 0) / rate / bmax).epsilon(1e-12));
    // MATLAB SolverMVA (default) on the identical model
    CHECK(r.Q(1, 0) == doctest::Approx(1.737033963).epsilon(1e-8));
    CHECK(r.Tp(1, 0) == doctest::Approx(2.262965855).epsilon(1e-8));
    CHECK(r.U(1, 0) == doctest::Approx(0.7543219517).epsilon(1e-8));
    // the CTMC is exact; amvald's cd term is approximate, so the bound is loose
    CHECK(r.Q(1, 0) == doctest::Approx(e.avg.QN(1, 0)).epsilon(0.05));
    CHECK(r.U(1, 0) == doctest::Approx(e.avg.UN(1, 0)).epsilon(0.05));
}
