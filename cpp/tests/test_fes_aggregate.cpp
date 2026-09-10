/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * fes_aggregate: Norton aggregation of a station subset into one FES.
 *
 * THE ORACLE IS THE THEOREM, NOT A STORED NUMBER. For a closed product-form
 * network the flow-equivalent replacement is EXACT at every population, so the
 * aggregated model must reproduce the original model's complement queue
 * lengths and throughputs to solver precision, and the FES must hold exactly
 * the sum of the queue lengths it replaced. That is a far stronger check than
 * comparing against a golden: a golden only says the code did not change,
 * while conservation says it is right.
 *
 * Each case below breaks a different piece if it is wrong:
 *  - the plain cycle checks the wiring;
 *  - INTERNAL FEEDBACK inside the subset checks the Norton escape factor, and
 *    is the case that fails loudly without it -- the isolated throughput counts
 *    every internal hop, so the FES would complete jobs 1/escape times too
 *    fast (2.5x at the 0.4 escape used here);
 *  - two classes check that the per-class tables and the class dependence are
 *    not collapsed onto a single scaling;
 *  - a Delay inside the subset checks the infinite-server branch of the
 *    isolated solve.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/fes/fes_aggregate.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/nc/solver_nc_runner.h"

namespace fes = line::fes;
namespace qn = line::qn;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

line::mva::AvgResult<double> solve(const qn::NetworkStruct<double>& sn) {
    line::nc::NcSolverOptions o;
    return line::nc::solver_nc_run_analyzer(sn, o);
}

/** Delay -> Q1 -> Q2 -> Delay, one closed class. */
qn::Network<double> cycle3(double njobs) {
    qn::Network<double> m("cycle3");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("aggregating a two-queue subset conserves the whole solution") {
    qn::Network<double> m = cycle3(3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::mva::AvgResult<double> ref = solve(sn);

    std::vector<std::size_t> sub;
    sub.push_back(2);  // Q1, 1-based station index
    sub.push_back(3);  // Q2
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub);

    const qn::NetworkStruct<double> sn2 = r.model.get_struct();
    CHECK(sn2.nstations == 2);  // the Delay plus the FES
    const line::mva::AvgResult<double> agg = solve(sn2);

    // The complement station is untouched, to solver precision.
    CHECK(agg.QN(0, 0) == doctest::Approx(ref.QN(0, 0)).epsilon(1e-6));
    CHECK(agg.TN(0, 0) == doctest::Approx(ref.TN(0, 0)).epsilon(1e-6));
    // The FES holds exactly what the subset held.
    CHECK(agg.QN(1, 0) == doctest::Approx(ref.QN(1, 0) + ref.QN(2, 0)).epsilon(1e-6));
    // The population is conserved.
    CHECK(agg.QN(0, 0) + agg.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-6));
    // Nothing left the subset, so the escape factor is one.
    CHECK(r.deagg.escape[0] == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("internal feedback inside the subset needs the escape factor") {
    // Q2 returns to Q1 with probability 0.6, so a job makes 1/0.4 = 2.5 subset
    // visits per escape. Without the escape correction the FES would serve 2.5
    // times too fast and the conservation checks below would all fail.
    qn::Network<double> m("fb");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 0.6);
    P.set(q2, d, 0.4);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::mva::AvgResult<double> ref = solve(sn);

    std::vector<std::size_t> sub;
    sub.push_back(2);
    sub.push_back(3);
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub);
    CHECK(r.deagg.escape[0] == doctest::Approx(0.4).epsilon(1e-9));

    const line::mva::AvgResult<double> agg = solve(r.model.get_struct());
    CHECK(agg.QN(0, 0) == doctest::Approx(ref.QN(0, 0)).epsilon(1e-5));
    CHECK(agg.TN(0, 0) == doctest::Approx(ref.TN(0, 0)).epsilon(1e-5));
    CHECK(agg.QN(1, 0) == doctest::Approx(ref.QN(1, 0) + ref.QN(2, 0)).epsilon(1e-5));
    // The FES completes at the EXTERNAL rate, which is the internal subset
    // throughput times the escape probability.
    CHECK(agg.TN(1, 0) == doctest::Approx(ref.TN(1, 0) * 0.4).epsilon(1e-5));
}

TEST_CASE("two classes keep their own throughput tables") {
    qn::Network<double> m("mc");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t a = m.add_closed_class("A", 2.0, d);
    const std::size_t b = m.add_closed_class("B", 2.0, d);
    m.set_service(d, a, Dist::exp_rate(1.0));
    m.set_service(q1, a, Dist::exp_rate(2.0));
    m.set_service(q2, a, Dist::exp_rate(3.0));
    m.set_service(d, b, Dist::exp_rate(0.5));
    m.set_service(q1, b, Dist::exp_rate(4.0));
    m.set_service(q2, b, Dist::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(a, a, d, q1, 1.0);
    P.set(a, a, q1, q2, 1.0);
    P.set(a, a, q2, d, 1.0);
    P.set(b, b, d, q1, 1.0);
    P.set(b, b, q1, q2, 1.0);
    P.set(b, b, q2, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::mva::AvgResult<double> ref = solve(sn);

    std::vector<std::size_t> sub;
    sub.push_back(2);
    sub.push_back(3);
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub);
    CHECK(r.deagg.throughputTable.size() == 2);
    // The two classes have different demands, so their tables must differ; a
    // single shared scaling is the failure this catches.
    bool differ = false;
    for (std::size_t i = 0; i < r.deagg.throughputTable[0].size(); ++i)
        if (std::fabs(r.deagg.throughputTable[0][i] - r.deagg.throughputTable[1][i]) > 1e-9)
            differ = true;
    CHECK(differ);

    const line::mva::AvgResult<double> agg = solve(r.model.get_struct());
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(agg.QN(0, c) == doctest::Approx(ref.QN(0, c)).epsilon(1e-5));
        CHECK(agg.TN(0, c) == doctest::Approx(ref.TN(0, c)).epsilon(1e-5));
        CHECK(agg.QN(1, c) == doctest::Approx(ref.QN(1, c) + ref.QN(2, c)).epsilon(1e-5));
    }
}

TEST_CASE("a Delay inside the subset takes the infinite-server branch") {
    // Q1 -> Think2 -> Q1 aggregated, leaving Q0 outside.
    qn::Network<double> m("withdelay");
    const std::size_t q0 = m.add_queue("Q0", SchedStrategy::PS);
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t d2 = m.add_delay("Think2");
    const std::size_t c = m.add_closed_class("C1", 3.0, q0);
    m.set_service(q0, c, Dist::exp_rate(1.5));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(d2, c, Dist::exp_rate(0.8));
    qn::RoutingMatrix<double> P;
    P.set(q0, q1, 1.0);
    P.set(q1, d2, 1.0);
    P.set(d2, q0, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    const line::mva::AvgResult<double> ref = solve(sn);

    std::vector<std::size_t> sub;
    sub.push_back(2);  // Q1
    sub.push_back(3);  // Think2
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub);
    CHECK(r.deagg.isolatedIsDelay.size() == 2);
    CHECK(r.deagg.isolatedIsDelay[0] == false);
    CHECK(r.deagg.isolatedIsDelay[1] == true);

    const line::mva::AvgResult<double> agg = solve(r.model.get_struct());
    CHECK(agg.QN(0, 0) == doctest::Approx(ref.QN(0, 0)).epsilon(1e-5));
    CHECK(agg.QN(1, 0) == doctest::Approx(ref.QN(1, 0) + ref.QN(2, 0)).epsilon(1e-5));
    CHECK(agg.TN(0, 0) == doctest::Approx(ref.TN(0, 0)).epsilon(1e-5));
}

TEST_CASE("the deaggregation record describes the transform it performed") {
    qn::Network<double> m = cycle3(2.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    std::vector<std::size_t> sub;
    sub.push_back(2);
    sub.push_back(3);
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub);

    CHECK(r.deagg.subsetIndices == sub);
    CHECK(r.deagg.complementIndices.size() == 1);
    CHECK(r.deagg.complementIndices[0] == 1);
    // The default cutoff is the class population.
    CHECK(r.deagg.cutoffs.size() == 1);
    CHECK(r.deagg.cutoffs[0] == 2);
    // The table is tabulated over the whole lattice, i.e. prod(cutoffs+1).
    CHECK(r.deagg.throughputTable[0].size() == 3u);
    // The isolated visit ratios are normalized within the subset.
    CHECK(r.deagg.isolatedVisits(0, 0) + r.deagg.isolatedVisits(1, 0) ==
          doctest::Approx(1.0).epsilon(1e-9));
    // The two stochastic complements are square on their own index sets.
    CHECK(r.deagg.stochCompSubset.rows() == 2);
    CHECK(r.deagg.stochCompComplement.rows() == 1);
    // The FES is the last node added, after the single complement station.
    CHECK(r.fesNode == 2);
    CHECK(r.deagg.fesNode == r.fesNode);
}

TEST_CASE("a smaller cutoff truncates the table without breaking the model") {
    qn::Network<double> m = cycle3(4.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    std::vector<std::size_t> sub;
    sub.push_back(2);
    sub.push_back(3);

    fes::FesOptions opt;
    opt.cutoffs.push_back(2);  // below the population of four
    fes::FesAggregateResult<double> r = fes::fes_aggregate(sn, sub, opt);
    CHECK(r.deagg.cutoffs[0] == 2);
    CHECK(r.deagg.throughputTable[0].size() == 3u);

    // The model still solves and conserves the population; the scaling
    // saturates past the cutoff, which is what the truncation means.
    const line::mva::AvgResult<double> agg = solve(r.model.get_struct());
    CHECK(agg.QN(0, 0) + agg.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-6));
}

TEST_CASE("the refusals are the reference's, by name") {
    qn::Network<double> m = cycle3(2.0);
    const qn::NetworkStruct<double> sn = m.get_struct();

    CHECK_THROWS_AS(fes::fes_aggregate(sn, std::vector<std::size_t>()), line::InputError);

    std::vector<std::size_t> all;
    all.push_back(1);
    all.push_back(2);
    all.push_back(3);
    CHECK_THROWS_AS(fes::fes_aggregate(sn, all), line::InputError);  // not a proper subset

    std::vector<std::size_t> bad;
    bad.push_back(9);
    CHECK_THROWS_AS(fes::fes_aggregate(sn, bad), line::InputError);

    // A cutoff vector of the wrong length is refused rather than padded.
    std::vector<std::size_t> sub;
    sub.push_back(2);
    fes::FesOptions opt;
    opt.cutoffs.push_back(1);
    opt.cutoffs.push_back(1);
    CHECK_THROWS_AS(fes::fes_aggregate(sn, sub, opt), line::InputError);
}

TEST_CASE("an open model is refused: Norton aggregation is a closed-network theorem") {
    qn::Network<double> m("open");
    const std::size_t src = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("O1");
    m.set_arrival(src, c, Dist::exp_rate(0.4));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, snk, 1.0);
    m.link(P);

    const qn::NetworkStruct<double> sn = m.get_struct();
    std::vector<std::size_t> sub;
    sub.push_back(2);
    CHECK_THROWS_AS(fes::fes_aggregate(sn, sub), line::InputError);
}
