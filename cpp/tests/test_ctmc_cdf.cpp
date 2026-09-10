/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The response-time DISTRIBUTION: `tag_chain` and the two CDF entry points it
 * feeds, `solver_ctmc_cdf_respt` and `solver_ctmc_cdf_sys_respt`.
 *
 * WHY THE ORACLES ARE WHAT THEY ARE. A CDF has no scalar to compare against a
 * reference table, so every check here is either a closed form or an identity
 * the construction must satisfy whatever the numbers turn out to be:
 *
 *   - TAGGING RELABELS A JOB, it does not change the dynamics. The twin class
 *     carries the same service process as the class it copies, so the tagged
 *     model is a LUMPING of the original: every aggregate the analyzer reports,
 *     summed back over the twin and its original, must be unchanged. That is the
 *     strongest available check on the transform, and it fails loudly if the
 *     population bookkeeping, the routing copy or the capacity split is wrong.
 *   - A CDF MUST BE A CDF: F(0) = 0, non-decreasing, bounded by 1, reaching 1.
 *     A sub-generator that is not sub-stochastic -- the usual symptom of a
 *     filtration summed on the wrong half of the synchronization -- breaks all
 *     four at once.
 *   - AT POPULATION 1 THERE IS NO QUEUEING, so the response time at a station is
 *     its service time exactly and the closed form is known.
 *   - ON A SINGLE-STATION CHAIN the system response time IS the station's, since
 *     the cycle is the visit.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/tag_chain.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cdf.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/**
 * Delay -> FCFS Queue -> Delay, one closed class.
 *
 * `completes` is set explicitly although it is now the default in this port as in
 * the other three: the CDF has no completion event to absorb at without it, so
 * stating it here keeps the fixture readable next to `cdf_cqn_incomplete`, which
 * clears it.
 */
qn::Network<double> cdf_cqn(double njobs, double rate_think, double rate_queue) {
    qn::Network<double> m("cdf_cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(rate_think));
    m.set_service(q, c, Dist::exp_rate(rate_queue));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    m.raw_struct().classes[c - 1].completes = true;
    return m;
}

/** The same model with `completes` CLEARED on its only class. */
qn::Network<double> cdf_cqn_incomplete(double njobs) {
    qn::Network<double> m("cdf_cqn_incomplete");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    m.raw_struct().classes[c - 1].completes = false;
    return m;
}

/** One FCFS Queue routing back to itself: a chain whose only station is that queue. */
qn::Network<double> cdf_single(double njobs, double rate) {
    qn::Network<double> m("cdf_single");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, q);
    m.set_service(q, c, Dist::exp_rate(rate));
    qn::RoutingMatrix<double> P;
    P.set(q, q, 1.0);
    m.link(P);
    m.raw_struct().classes[c - 1].completes = true;
    return m;
}

/** Source -> FCFS Queue -> Sink, one open class. */
qn::Network<double> cdf_open(double lambda, double mu, int K) {
    qn::Network<double> m("cdf_open");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/**
 * Every property that makes a curve a distribution function.
 *
 * The monotonicity slack is one step's rounding: `F = 1 - sum(v)` with v carried
 * by repeated multiplication by a sub-stochastic matrix, so a genuine decrease
 * means the sub-generator was not sub-stochastic and not that the arithmetic
 * drifted. `tol` is the truncation tolerance the entry point uses, and the last
 * point is where the curve crossed it.
 */
void check_is_cdf(const ctmc::CdfCurve<double>& c, double tol) {
    REQUIRE(!c.empty());
    REQUIRE(c.t.size() == c.F.size());
    CHECK(c.t[0] == 0.0);
    CHECK(std::fabs(c.F[0]) < 1e-12);
    // Scanned in plain C++ and reported as three verdicts rather than three
    // assertions per grid point: a curve runs to thousands of points and the
    // per-point form would bury the rest of the file in its assertion count.
    bool grid_ok = true, monotone = true, bounded = true;
    for (std::size_t k = 1; k < c.F.size(); ++k) {
        if (!(c.t[k] > c.t[k - 1])) grid_ok = false;
        if (c.F[k] < c.F[k - 1] - 1e-12) monotone = false;
        if (c.F[k] > 1.0 + 1e-12) bounded = false;
    }
    CHECK(grid_ok);
    CHECK(monotone);
    CHECK(bounded);
    CHECK(c.F.back() > 1.0 - tol);
}

/** Total population of a struct, over every class. */
double total_pop(const qn::NetworkStruct<double>& sn) {
    double s = 0;
    for (std::size_t r = 0; r < sn.classes.size(); ++r) s += sn.classes[r].population;
    return s;
}

}  // namespace

TEST_CASE("ctmc cdf: tagging moves one job into a twin class and conserves the population") {
    qn::Network<double> m = cdf_cqn(3.0, 0.5, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.nchains == 1);

    const qn::TaggedChain<double> tg = qn::tag_chain(sn, 1, 1);

    // One twin per class of the chain, so the tagged job stays tagged across
    // every class switch it can make.
    REQUIRE(tg.tagged.size() == sn.inchain[0].size());
    CHECK(tg.V.classes.size() == sn.classes.size() + sn.inchain[0].size());
    REQUIRE(tg.taggedjob != 0);
    for (std::size_t a = 0; a < tg.tagged.size(); ++a) {
        // Every twin is a NEW class carrying the name of the one it copies, and
        // the pairing is what a caller reads its results back through.
        CHECK(tg.tagged[a] > sn.classes.size());
        CHECK(tg.V.classes[tg.tagged[a] - 1].name ==
              sn.classes[tg.orig[a] - 1].name + ".tagged");
        CHECK(tg.V.classes[tg.tagged[a] - 1].refstat == sn.classes[tg.orig[a] - 1].refstat);
    }

    // The job was MOVED, not minted: one out of the tagged class, one into its
    // twin, and every other twin empty.
    CHECK(tg.V.classes[0].population == doctest::Approx(2.0));
    CHECK(tg.V.classes[tg.taggedjob - 1].population == doctest::Approx(1.0));
    CHECK(total_pop(tg.V) == doctest::Approx(total_pop(sn)));

    // The twin is not a fresh class with default parameters: it inherits the
    // service of the class it copies, which is what makes the tagged model a
    // lumping of the original rather than a different model.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(!tg.V.service[i][tg.taggedjob - 1].disabled);
        CHECK(tg.V.service[i][tg.taggedjob - 1].mean ==
              doctest::Approx(sn.service[i][0].mean).epsilon(1e-12));
    }
    // The twins form their own chain, since nothing routes between the tagged
    // and the untagged block.
    CHECK(tg.V.nchains == sn.nchains + 1);
}

TEST_CASE("ctmc cdf: the tagged model reports the same aggregates as the untagged one") {
    // Both classes of the tagged model carry the SAME exponential service, so
    // relabelling one job cannot change any aggregate the analyzer reports. The
    // tagged chain is a strictly larger state space reaching identical numbers,
    // which no coincidence of parameters could produce.
    qn::Network<double> m = cdf_cqn(3.0, 0.5, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const line::mva::AvgResult<double> base = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    const qn::TaggedChain<double> tg = qn::tag_chain(sn, 1, 1);
    const line::mva::AvgResult<double> tagd = ctmc::solver_ctmc_run_analyzer(tg.V, ctmc::CtmcOptions());

    for (std::size_t i = 0; i < sn.nstations; ++i) {
        double q0 = 0, u0 = 0, x0 = 0;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            q0 += base.QN(i, r);
            u0 += base.UN(i, r);
            x0 += base.TN(i, r);
        }
        double q1 = 0, u1 = 0, x1 = 0;
        for (std::size_t r = 0; r < tg.V.nclasses; ++r) {
            q1 += tagd.QN(i, r);
            u1 += tagd.UN(i, r);
            x1 += tagd.TN(i, r);
        }
        CHECK(q1 == doctest::Approx(q0).epsilon(1e-8));
        CHECK(u1 == doctest::Approx(u0).epsilon(1e-8));
        CHECK(x1 == doctest::Approx(x0).epsilon(1e-8));
    }
    // Population is conserved exactly on both sides, which is the invariant the
    // capacity split is most likely to break: taking a buffer slot away from the
    // original class without giving one to the twin loses a job outright.
    double n0 = 0, n1 = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) n0 += base.QN(i, r);
    for (std::size_t i = 0; i < tg.V.nstations; ++i)
        for (std::size_t r = 0; r < tg.V.nclasses; ++r) n1 += tagd.QN(i, r);
    CHECK(n0 == doctest::Approx(3.0).epsilon(1e-9));
    CHECK(n1 == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("ctmc cdf: at population one the response time is the service time exactly") {
    // WITH ONE JOB IN THE NETWORK THERE IS NOTHING TO QUEUE BEHIND. The job's
    // sojourn at a station is one exponential service, so the response-time CDF
    // is that exponential's, at BOTH stations:
    //
    //     F_i(t) = 1 - exp(-mu_i t)
    //
    // This is the only place in the file where a closed form pins the ABSOLUTE
    // scale of the curve rather than an identity between two of them.
    const double mu_d = 0.5, mu_q = 3.0;
    qn::Network<double> m = cdf_cqn(1.0, mu_d, mu_q);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<std::vector<ctmc::CdfCurve<double>>> RD =
        ctmc::solver_ctmc_cdf_respt(sn, ctmc::CtmcOptions());
    REQUIRE(sn.nstations == 2);
    REQUIRE(RD.size() == sn.nstations);
    REQUIRE(RD[0].size() == sn.nclasses);

    const double rate[2] = {mu_d, mu_q};
    for (std::size_t i = 0; i < 2; ++i) {
        const ctmc::CdfCurve<double>& c = RD[i][0];
        check_is_cdf(c, line::lang::GlobalConstants::CoarseTol);
        for (std::size_t k = 0; k < c.t.size(); k += 137) {
            const double want = 1.0 - std::exp(-rate[i] * c.t[k]);
            // Below 1e-6 the relative comparison is meaningless: the curve is
            // still at the origin and any absolute error dominates.
            if (want < 1e-6) continue;
            CHECK(c.F[k] == doctest::Approx(want).epsilon(1e-7));
        }
    }
}

TEST_CASE("ctmc cdf: the system response time at population one is the hypoexponential cycle") {
    // The split is taken at the tagged job's arrival at its own REFERENCE
    // station, so the passage is a full cycle: think, then queue. At population
    // one those are two independent exponentials in series, and the cycle time
    // is hypoexponential,
    //
    //     F(t) = 1 - (mu_q exp(-mu_d t) - mu_d exp(-mu_q t)) / (mu_q - mu_d)
    //
    // which requires mu_d != mu_q to be non-degenerate.
    const double mu_d = 0.5, mu_q = 3.0;
    qn::Network<double> m = cdf_cqn(1.0, mu_d, mu_q);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const std::vector<ctmc::CdfCurve<double>> RD =
        ctmc::solver_ctmc_cdf_sys_respt(sn, ctmc::CtmcOptions());
    REQUIRE(RD.size() == sn.nchains);
    const ctmc::CdfCurve<double>& c = RD[0];
    check_is_cdf(c, line::lang::GlobalConstants::FineTol);

    for (std::size_t k = 0; k < c.t.size(); k += 37) {
        const double t = c.t[k];
        const double want =
            1.0 - (mu_q * std::exp(-mu_d * t) - mu_d * std::exp(-mu_q * t)) / (mu_q - mu_d);
        if (want < 1e-6) continue;
        CHECK(c.F[k] == doctest::Approx(want).epsilon(1e-7));
        // The cycle is stochastically LARGER than either leg alone, so its CDF
        // lies strictly below the think leg's. A split taken at the wrong
        // station would return one leg and land on that curve instead.
        CHECK(c.F[k] < 1.0 - std::exp(-mu_d * t));
    }
}

TEST_CASE("ctmc cdf: on a single-station chain the system response time is the station's") {
    // The chain's only station is its reference station, so the cycle IS the
    // visit: the arrival at the reference station that ends one passage is the
    // same event as the departure that ends the station's. The two entry points
    // reach that from opposite halves of the synchronization -- one filters the
    // PASSIVE arrival, the other the ACTIVE departure -- so agreement here is a
    // real check that both halves name the same transition.
    const double mu = 2.0;
    qn::Network<double> m = cdf_single(2.0, mu);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.nstations == 1);
    REQUIRE(sn.nchains == 1);

    const std::vector<std::vector<ctmc::CdfCurve<double>>> RD =
        ctmc::solver_ctmc_cdf_respt(sn, ctmc::CtmcOptions());
    const std::vector<ctmc::CdfCurve<double>> SD =
        ctmc::solver_ctmc_cdf_sys_respt(sn, ctmc::CtmcOptions());
    const ctmc::CdfCurve<double>& st = RD[0][0];
    const ctmc::CdfCurve<double>& sy = SD[0];
    check_is_cdf(st, line::lang::GlobalConstants::CoarseTol);
    check_is_cdf(sy, line::lang::GlobalConstants::FineTol);

    // The two grids share a horizon and differ only in how finely they cut it,
    // 100000 intervals against 10000, so every system point lands on a station
    // point exactly ten apart.
    std::size_t compared = 0;
    double worst_t = 0, worst_F = 0;
    for (std::size_t k = 0; k < sy.t.size() && 10 * k < st.t.size(); ++k) {
        worst_t = std::max(worst_t, std::fabs(sy.t[k] - st.t[10 * k]));
        worst_F = std::max(worst_F, std::fabs(sy.F[k] - st.F[10 * k]));
        ++compared;
    }
    CHECK(worst_t < 1e-9);
    CHECK(worst_F < 1e-9);
    // A vacuous loop would pass silently; the station curve runs to F = 1-1e-3
    // and the system curve is sampled ten times more coarsely, so there is real
    // overlap to compare.
    CHECK(compared > 100);
}

TEST_CASE("ctmc cdf: every class of a chain shares the chain's curve at a station") {
    // The tagged job carries its class with it, so "the tagged job is at station
    // i" is one event per chain and not one per class. The reference reaches the
    // same conclusion through a loop variable that shadows its own index; here
    // it is computed once and written to every column of the chain.
    qn::Network<double> m("cdf_twoclass");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("A", 1.0, d);
    const std::size_t c2 = m.add_closed_class("B", 0.0, d);
    m.set_service(d, c1, Dist::exp_rate(0.5));
    m.set_service(q, c1, Dist::exp_rate(3.0));
    m.set_service(d, c2, Dist::exp_rate(0.5));
    m.set_service(q, c2, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    // A switches to B on the way to the queue and back to A on the way home, so
    // the two classes are one chain and the tagged job changes class every hop.
    P.set(c1, c2, d, q, 1.0);
    P.set(c2, c1, q, d, 1.0);
    m.link(P);
    m.raw_struct().classes[c1 - 1].completes = true;
    m.raw_struct().classes[c2 - 1].completes = true;

    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.nchains == 1);
    REQUIRE(sn.inchain[0].size() == 2);

    const std::vector<std::vector<ctmc::CdfCurve<double>>> RD =
        ctmc::solver_ctmc_cdf_respt(sn, ctmc::CtmcOptions());
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        REQUIRE(RD[i][0].t.size() == RD[i][1].t.size());
        double worst = 0;
        for (std::size_t k = 0; k < RD[i][0].t.size(); ++k)
            worst = std::max(worst, std::fabs(RD[i][0].F[k] - RD[i][1].F[k]));
        // Bit-for-bit: the two columns are copies of one computed curve, so any
        // difference at all would mean they were computed separately.
        CHECK(worst == 0.0);
        check_is_cdf(RD[i][0], line::lang::GlobalConstants::CoarseTol);
    }
}

TEST_CASE("ctmc cdf: an open model and an open chain are refused by name") {
    qn::Network<double> m = cdf_open(0.6, 1.0, 4);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = 4;
    // getCdfRespT's own guard: a truncated open chain has no population to move
    // one job out of.
    CHECK_THROWS_AS(ctmc::solver_ctmc_cdf_respt(sn, opt), line::UnsupportedError);
    // The reference's getCdfSysRespT has no such guard; tag_chain supplies one,
    // so the refusal arrives from the transform rather than from the entry point.
    CHECK_THROWS_AS(ctmc::solver_ctmc_cdf_sys_respt(sn, opt), line::UnsupportedError);
    CHECK_THROWS_AS(qn::tag_chain(sn, 1, 1), line::UnsupportedError);
}

TEST_CASE("ctmc cdf: a chain with no completing class is refused rather than answered") {
    // The fixture CLEARS `completes`, which now defaults to true here as it does
    // in MATLAB, the JAR and native Python. An empty selection is refused by name
    // instead of returning a curve built from an all-zero filtration, which
    // would be a CDF that never reaches 1.
    qn::Network<double> m = cdf_cqn_incomplete(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK_THROWS_AS(ctmc::solver_ctmc_cdf_respt(sn, ctmc::CtmcOptions()), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::solver_ctmc_cdf_sys_respt(sn, ctmc::CtmcOptions()),
                    line::UnsupportedError);
}

TEST_CASE("ctmc cdf: tag_chain refuses every feature whose per-class tables it cannot widen") {
    // The fixtures are built by mutating a COPY of a solved struct rather than
    // through the builder: the point is the refusal, and going through the
    // builder would drag in the rest of each feature's setup for no added
    // coverage. Each of these declares something indexed by class that the
    // transform would have to invent an entry for.
    qn::Network<double> m = cdf_cqn(2.0, 0.5, 3.0);
    const qn::NetworkStruct<double> good = m.get_struct();
    REQUIRE_NOTHROW(qn::tag_chain(good, 1, 1));

    {
        qn::NetworkStruct<double> bad = good;
        bad.issignal.assign(bad.nclasses, false);
        bad.issignal[0] = true;
        CHECK_THROWS_AS(qn::tag_chain(bad, 1, 1), line::UnsupportedError);
    }
    {
        qn::NetworkStruct<double> bad = good;
        bad.regions.resize(1);
        CHECK_THROWS_AS(qn::tag_chain(bad, 1, 1), line::UnsupportedError);
    }
    {
        qn::NetworkStruct<double> bad = good;
        bad.syncreply.assign(bad.nclasses, 0);
        bad.syncreply[0] = 1;
        CHECK_THROWS_AS(qn::tag_chain(bad, 1, 1), line::UnsupportedError);
    }
    {
        qn::NetworkStruct<double> bad = good;
        bad.stations[0].cdscaling = [](const std::vector<double>&) {
            return std::vector<double>();
        };
        CHECK_THROWS_AS(qn::tag_chain(bad, 1, 1), line::UnsupportedError);
    }
    {
        qn::NetworkStruct<double> bad = good;
        bad.stations[0].svc_rate_fun = [](const std::vector<std::size_t>&) { return 1.0; };
        CHECK_THROWS_AS(qn::tag_chain(bad, 1, 1), line::UnsupportedError);
    }
}

TEST_CASE("ctmc cdf: tag_chain rejects a chain or class it cannot name") {
    qn::Network<double> m = cdf_cqn(2.0, 0.5, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK_THROWS_AS(qn::tag_chain(sn, 0, 1), line::InputError);
    CHECK_THROWS_AS(qn::tag_chain(sn, sn.nchains + 1, 1), line::InputError);
    CHECK_THROWS_AS(qn::tag_chain(sn, 1, sn.nclasses + 1), line::InputError);
    // A class with no job cannot give one up, and the reference picks the class
    // to tag by exactly this predicate.
    qn::Network<double> empty = cdf_cqn(0.0, 0.5, 3.0);
    CHECK_THROWS_AS(qn::tag_chain(empty.get_struct(), 1, 1), line::InputError);
}
