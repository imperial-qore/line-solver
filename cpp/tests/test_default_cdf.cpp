/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The base-class response-time CDF fallback (`solver_default_cdf.h`) and the
 * per-entry LDES LN ecdf (`ldes_ln_cdf_respt`), the two pieces that close the
 * get*Cdf* support matrix against the MATLAB reference.
 *
 * THE FALLBACK IS A CONTRACT, NOT AN APPROXIMATION TO TEST LOOSELY: it must be
 * exactly `@@NetworkSolver/getCdfRespT.m` -- 100 quantile points 0.001..0.999,
 * t = -log(1-q) * RN(i,r), a Source left empty, a served cell with no finite
 * positive mean reduced to the point mass [1, 0] -- because a caller reading
 * it beside MATLAB's must see the same table.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ldes/ldes_ln_engine.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/solver_default_cdf.h"

namespace qn = line::qn;
namespace mva = line::mva;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source(rate 1) -> FCFS Queue(rate 2) -> Sink, one open class. */
qn::Network<double> mm1() {
    qn::Network<double> m("default-cdf-mm1");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("solver_default_cdf_respt is the reference exponential fallback") {
    qn::Network<double> net = mm1();
    const qn::NetworkStruct<double>& sn = net.get_struct();
    mva::MvaOptions opt;
    line::Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(sn, opt, init);

    const std::vector<std::vector<line::solvers::DefaultCdfCurve>> RD =
        line::solvers::solver_default_cdf_respt<double>(sn, r.RN);
    REQUIRE(RD.size() == sn.nstations);

    // The Source row is left empty, as the reference leaves its cells empty
    std::size_t src = 0, q = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].nodetype == line::lang::NodeType::Source) src = i;
        else q = i;
    }
    CHECK(RD[src][0].t.empty());

    // The served cell is the 100-point exponential law with the MVA mean
    const line::solvers::DefaultCdfCurve& cell = RD[q][0];
    REQUIRE(cell.t.size() == 100);
    const double rn = r.RN(q, 0);
    CHECK(rn > 0.0);
    CHECK(cell.F.front() == doctest::Approx(0.001));
    CHECK(cell.F.back() == doctest::Approx(0.999));
    for (std::size_t j = 0; j < cell.t.size(); j += 20)
        CHECK(cell.t[j] == doctest::Approx(-std::log(1.0 - cell.F[j]) * rn).epsilon(1e-12));
    // M/M/1: RN = 1/(mu - lambda) = 1, so the median sits at ln 2
    CHECK(rn == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("ctmc_cdf_firstpasst reproduces the exact passage mean") {
    // Think(rate 1) <-> FCFS Queue(rate 2), 2 closed jobs: a three-state chain
    qn::Network<double> m("firstpasst-cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    line::ctmc::CtmcOptions opt;
    const line::ctmc::CtmcSolution<double> sol = line::ctmc::solver_ctmc_analyzer(sn, opt);
    const std::size_t n = sol.chain.Q.rows();
    REQUIRE(n > 1);

    // Into the LAST state, from the conditional stationary law on the rest
    line::Matrix<double> B(1, 1);
    B(0, 0) = static_cast<double>(n);
    const line::ctmc::CtmcFirstPassage fp =
        line::ctmc::ctmc_cdf_firstpasst<double>(sn, sol, line::Matrix<double>(0, 0), B);
    REQUIRE(fp.t.size() == 1000);
    REQUIRE(fp.target.size() == 1);
    CHECK(fp.source.empty());
    double prev = 0.0;
    for (std::size_t i = 0; i < fp.F.size(); ++i) {
        CHECK(fp.F[i] >= prev - 1e-12);
        CHECK(fp.F[i] >= 0.0);
        CHECK(fp.F[i] <= 1.0);
        prev = fp.F[i];
    }
    CHECK(fp.F.back() > 0.99);

    // The curve's mean integrates to the first moment ctmc_passage_moments
    // solves for exactly -- the check that curve and moments answer ONE question
    std::vector<double> pi0;  // empty: same conditional stationary law
    const line::mc::PassageMoments<double> pm =
        line::mc::ctmc_passage_moments(sol.chain.Q, pi0, fp.target, 1);
    // Trapezoid, not a right-endpoint sum: the grid is thor/999 = 0.1 wide
    // against a mean of 0.92, so the endpoint rule's -dt/2 bias is 5%.
    double mean_from_curve = 0.0;
    for (std::size_t i = 1; i < fp.t.size(); ++i)
        mean_from_curve +=
            0.5 * ((1.0 - fp.F[i]) + (1.0 - fp.F[i - 1])) * (fp.t[i] - fp.t[i - 1]);
    CHECK(mean_from_curve == doctest::Approx(pm.m[0]).epsilon(1e-2));
}

TEST_CASE("ldes_ln_cdf_respt collapses ties keeping the largest F") {
    line::ldes::engine::LnResult r;
    r.entry_resp_samples.resize(2);
    r.entry_resp_samples[0] = {2.0, 1.0, 2.0, 3.0};
    // entry 1 observed nothing

    const std::vector<line::ldes::engine::LnEntryCdf> cdf =
        line::ldes::engine::ldes_ln_cdf_respt(r, 2);
    REQUIRE(cdf.size() == 2);
    REQUIRE(cdf[0].t.size() == 3);
    CHECK(cdf[0].t[0] == 1.0);
    CHECK(cdf[0].t[1] == 2.0);
    CHECK(cdf[0].t[2] == 3.0);
    CHECK(cdf[0].F[0] == doctest::Approx(0.25));
    CHECK(cdf[0].F[1] == doctest::Approx(0.75));
    CHECK(cdf[0].F[2] == doctest::Approx(1.0));
    CHECK(cdf[1].t.empty());
}
