/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The peak rate that normalizes utilization is a MODEL INPUT, and EVERY reader
 * refuses a model that omits it.
 *
 * Utilization at a rate-dependent station is reported as U = T*E[S]/peak, the
 * same fraction-of-capacity as the T*S/c of an ordinary multiserver station.
 * Load dependence always knows its peak, max(c, max alpha), because the user
 * supplied every alpha. CLASS and JOINT dependence do not: the scaling is a
 * HANDLE, so recovering max_n beta(n) means sweeping the population lattice --
 * which needs a bound the handle does not carry, and an open class has no bound
 * at all. So the peak is declared, and a declaration without one is a defect.
 *
 * `set_class_dependence` accepts an empty peak and the READER refuses, which is
 * what lets the negative cases below exist at all. That contract is only worth
 * anything if every reader honours it: SolverMVA used to be the one that did
 * not, silently writing a column of ZEROS into U, and the LDES engine used to
 * ignore the declared peak and normalize by the server count instead.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ldes/ldes_engine.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace mva = line::mva;
namespace nc = line::nc;
namespace ssa = line::ssa;
namespace ldes = line::ldes;
using line::InputError;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** beta(n) = min(n_1, 2): a handle whose lattice peak is 2. */
std::vector<double> beta_two(const std::vector<double>& n) {
    return std::vector<double>(1, std::min(n[0], 2.0));
}

/** Closed Delay -> Queue, one class of N jobs, class dependence with `peak`. */
qn::Network<double> cd_model(const std::vector<double>& peak,
                             SchedStrategy sched = SchedStrategy::PS) {
    qn::Network<double> m("cdpeak");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", sched);
    const std::size_t c = m.add_closed_class("Class1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.0));
    m.set_class_dependence(q, beta_two, peak);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

const std::vector<double> NO_PEAK;
const std::vector<double> PEAK_TWO(1, 2.0);

}  // namespace

TEST_CASE("peak required: SolverMVA refuses a class dependence with no peak") {
    // it used to write U = 0 into the column instead, which reads as a station
    // that is never busy rather than as the model defect it is
    qn::Network<double> m = cd_model(NO_PEAK);
    mva::MvaOptions opt;
    line::Matrix<double> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), InputError);
}

TEST_CASE("peak required: the refusal names the setter that should carry it") {
    qn::Network<double> m = cd_model(NO_PEAK);
    mva::MvaOptions opt;
    line::Matrix<double> init;
    std::string msg;
    try {
        mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
    } catch (const InputError& e) {
        msg = e.what();
    }
    CHECK(msg.find("setClassDependence") != std::string::npos);
    CHECK(msg.find("Queue") != std::string::npos);
}

TEST_CASE("peak required: the LDES engine refuses it too") {
    // this engine used to normalize a class-dependent station by its SERVER
    // COUNT, so it neither used the declared peak nor noticed its absence
    qn::Network<double> m = cd_model(NO_PEAK);
    ldes::LdesOptions o;
    o.samples = 2000;
    o.seed = 23000;
    CHECK_THROWS_AS(ldes::ldes_engine_solve(m.get_struct(), o), InputError);
}

TEST_CASE("peak required: CTMC, NC and the serial SSA engine agree") {
    // pinned together so that a reader added later cannot quietly opt out
    qn::Network<double> mc = cd_model(NO_PEAK);
    CHECK_THROWS(ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions()));

    qn::Network<double> mn = cd_model(NO_PEAK);
    CHECK_THROWS(nc::solver_nc_run_analyzer(mn.get_struct(), nc::NcSolverOptions()));

    qn::Network<double> ms = cd_model(NO_PEAK);
    ssa::SsaOptions so;
    so.method = "serial";
    so.samples = 2000;
    so.seed = 23000;
    CHECK_THROWS(ssa::solver_ssa(ms.get_struct(), so));
}

TEST_CASE("peak required: a declared peak is accepted and used, not re-derived") {
    // beta maxes at 2 over the lattice but is DECLARED as 4, and every solver
    // must divide by 4: the declaration is the normalizer. A reader that swept
    // the handle instead would report twice this station's utilization -- and
    // that sweep is exactly what the Java LDES engine used to do.
    const std::vector<double> declared_four(1, 4.0);

    qn::Network<double> two = cd_model(PEAK_TWO);
    const mva::AvgResult<double> a =
        ctmc::solver_ctmc_run_analyzer(two.get_struct(), ctmc::CtmcOptions());

    qn::Network<double> four = cd_model(declared_four);
    const mva::AvgResult<double> b =
        ctmc::solver_ctmc_run_analyzer(four.get_struct(), ctmc::CtmcOptions());

    // the CHAIN is identical -- the peak touches only the utilization column
    CHECK(b.QN(1, 0) == doctest::Approx(a.QN(1, 0)).epsilon(1e-12));
    CHECK(b.TN(1, 0) == doctest::Approx(a.TN(1, 0)).epsilon(1e-12));
    CHECK(a.UN(1, 0) > 0.0);
    CHECK(b.UN(1, 0) == doctest::Approx(a.UN(1, 0) * 2.0 / 4.0).epsilon(1e-9));
}

TEST_CASE("peak required: joint dependence refuses an empty peak at DECLARATION") {
    // its twin takes the peak as a required argument rather than a defaulted
    // one, so the refusal lands on the setter and never reaches a solver
    qn::Network<double> m("jdpeak");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("Class1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(1.0));
    CHECK_THROWS_AS(m.set_joint_dependence(q, beta_two, NO_PEAK), InputError);
}

TEST_CASE("peak required: the LDES engine reports the same U as CTMC at a cd station") {
    // Once the peak is declared the two must agree on the CONVENTION, not just
    // on the refusal. beta(n) = min(n,2) emulates a second server, so a station
    // read as a fraction of TIME occupied would report up to 2x what one read as
    // a fraction of WORK delivered does -- and this engine used to read time and
    // then divide it by the server count, which is neither.
    qn::Network<double> mc = cd_model(PEAK_TWO, SchedStrategy::FCFS);
    const mva::AvgResult<double> e =
        ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());

    qn::Network<double> ml = cd_model(PEAK_TWO, SchedStrategy::FCFS);
    ldes::LdesOptions o;
    o.samples = 300000;
    o.seed = 23000;
    const ldes::LdesResult r = ldes::ldes_engine_solve(ml.get_struct(), o);

    CHECK(e.UN(1, 0) > 0.0);
    CHECK(r.QN(1, 0) == doctest::Approx(e.QN(1, 0)).epsilon(3e-2));
    CHECK(r.TN(1, 0) == doctest::Approx(e.TN(1, 0)).epsilon(3e-2));
    CHECK(r.UN(1, 0) == doctest::Approx(e.UN(1, 0)).epsilon(3e-2));
}

TEST_CASE("peak required: a SHARING station applies the class dependence too") {
    // `ps_advance` and `ps_reschedule` built their shares from `ps_shares` and
    // the breakdown clock and never read the handle, so a PS station carrying a
    // `cdscaling` ran its whole sample path at the UNSCALED rates. THE TELL IS
    // THROUGHPUT, not utilization: on this model the exact answer is 1.4118 and
    // the engine returned 0.9348, a 34% error in a quantity no utilization
    // convention can explain.
    qn::Network<double> mc = cd_model(PEAK_TWO, SchedStrategy::PS);
    const mva::AvgResult<double> e =
        ctmc::solver_ctmc_run_analyzer(mc.get_struct(), ctmc::CtmcOptions());
    // the closed form: birth (3-n), death min(n,2), so pi ~ [1, 3, 3, 1.5]
    CHECK(e.TN(1, 0) == doctest::Approx(12.0 / 8.5).epsilon(1e-9));
    CHECK(e.QN(1, 0) == doctest::Approx(13.5 / 8.5).epsilon(1e-9));

    qn::Network<double> ml = cd_model(PEAK_TWO, SchedStrategy::PS);
    ldes::LdesOptions o;
    o.samples = 300000;
    o.seed = 23000;
    const ldes::LdesResult r = ldes::ldes_engine_solve(ml.get_struct(), o);
    CHECK(r.TN(1, 0) == doctest::Approx(e.TN(1, 0)).epsilon(3e-2));
    CHECK(r.QN(1, 0) == doctest::Approx(e.QN(1, 0)).epsilon(3e-2));
    CHECK(r.UN(1, 0) == doctest::Approx(e.UN(1, 0)).epsilon(3e-2));
}
