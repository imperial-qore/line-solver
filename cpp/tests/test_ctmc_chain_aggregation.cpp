/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `CtmcOptions::chain_aggregation`, the first solver consumer of the two
 * class-level ModelAdapter transforms.
 *
 * `api::sn_aggregate_chains` collapses every chain onto a single class and
 * `mva::sn_deaggregate_chain_results` maps chain-level metrics back through
 * alpha. Both existed in all four codebases with no solver consumer at all: the
 * transform was exercised by examples and tests only, so nothing in the solver
 * stack depended on it and a defect in it could not surface as a wrong answer.
 *
 * THE ORACLE IS AN IDENTITY, not a golden. On a PRODUCT-FORM model the chain is
 * the unit MVA and convolution already solve in, so the aggregation is exact and
 * the aggregated solve must reproduce the exact multiclass CTMC table to machine
 * precision, station by station and class by class. That is a statement about
 * the transform and the deaggregation together, and it is what fails if either
 * mis-derives alpha.
 *
 * The trade the option makes is the state space: the aggregated model carries
 * one class per chain, so a model with several classes in one chain is solved
 * over a strictly smaller chain. The exactness is what makes the trade free
 * here; on a non-product-form model one aggregate service law replaces the
 * per-class ones and the answer becomes an approximation, which is why the
 * option is off by default.
 */
#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace {

using DD = line::lang::Distrib<double>;
using Net = line::qn::Network<double>;
using Routing = line::qn::RoutingMatrix<double>;
using line::lang::SchedStrategy;

/**
 * Delay -> Q1 -> Q2 -> Delay with the class switch on the Q1 -> Q2 link, so the
 * two classes are ONE chain and the aggregation is not the identity.
 */
Net two_class_one_chain(double n = 3.0) {
    Net m("agg");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", n, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, DD::exp_rate(1.0));
    m.set_service(d, c2, DD::exp_rate(1.0));
    m.set_service(q1, c1, DD::exp_rate(2.0));
    m.set_service(q1, c2, DD::exp_rate(2.0));
    m.set_service(q2, c1, DD::exp_rate(3.0));
    m.set_service(q2, c2, DD::exp_rate(3.0));
    Routing P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c1, c2, q1, q2, 1.0);
    P.set(c2, c1, q2, d, 1.0);
    m.link(P);
    return m;
}

/** One class per chain: the transform is the identity and the guard declines it. */
Net one_class() {
    Net m("plain");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, DD::exp_rate(1.0));
    m.set_service(q, c, DD::exp_rate(2.0));
    Routing P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

line::mva::AvgResult<double> solve(Net& m, bool aggregate) {
    line::ctmc::CtmcOptions o;
    o.chain_aggregation = aggregate;
    return line::ctmc::solver_ctmc_run_analyzer(m.get_struct(), o);
}

}  // namespace

TEST_CASE("chain aggregation reproduces the exact table on a product-form model") {
    Net a = two_class_one_chain();
    Net b = two_class_one_chain();
    const line::mva::AvgResult<double> exact = solve(a, false);
    const line::mva::AvgResult<double> agg = solve(b, true);

    REQUIRE(agg.QN.rows() == exact.QN.rows());
    REQUIRE(agg.QN.cols() == exact.QN.cols());
    for (std::size_t i = 0; i < exact.QN.rows(); ++i)
        for (std::size_t k = 0; k < exact.QN.cols(); ++k) {
            CHECK(agg.QN(i, k) == doctest::Approx(exact.QN(i, k)).epsilon(1e-9));
            CHECK(agg.UN(i, k) == doctest::Approx(exact.UN(i, k)).epsilon(1e-9));
            CHECK(agg.TN(i, k) == doctest::Approx(exact.TN(i, k)).epsilon(1e-9));
        }
}

TEST_CASE("the aggregated solve names itself") {
    Net b = two_class_one_chain();
    const line::mva::AvgResult<double> agg = solve(b, true);
    // A caller reading the table has to be able to tell that the number came
    // from the reduced model rather than from the exact chain.
    CHECK(agg.actualmethod == "default/chainaggr");
}

TEST_CASE("the option is declined where the transform is the identity") {
    Net a = one_class();
    Net b = one_class();
    const line::mva::AvgResult<double> exact = solve(a, false);
    const line::mva::AvgResult<double> agg = solve(b, true);
    // nchains == nclasses, so the guard sends it down the ordinary path and the
    // method is NOT renamed.
    CHECK(agg.actualmethod != "default/chainaggr");
    for (std::size_t i = 0; i < exact.QN.rows(); ++i)
        CHECK(agg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-12));
}

TEST_CASE("flow is conserved through the deaggregation") {
    Net b = two_class_one_chain();
    const line::mva::AvgResult<double> agg = solve(b, true);
    // Every station on the cycle is visited once per chain traversal, so its
    // throughput summed over the classes must be the chain's.
    double ref = 0.0;
    for (std::size_t k = 0; k < agg.TN.cols(); ++k) ref += agg.TN(0, k);
    REQUIRE(ref > 0.0);
    for (std::size_t i = 1; i < agg.TN.rows(); ++i) {
        double t = 0.0;
        for (std::size_t k = 0; k < agg.TN.cols(); ++k) t += agg.TN(i, k);
        CHECK(t == doctest::Approx(ref).epsilon(1e-9));
    }
}

TEST_CASE("the population is conserved through the deaggregation") {
    Net b = two_class_one_chain(5.0);
    const line::mva::AvgResult<double> agg = solve(b, true);
    double total = 0.0;
    for (std::size_t i = 0; i < agg.QN.rows(); ++i)
        for (std::size_t k = 0; k < agg.QN.cols(); ++k) total += agg.QN(i, k);
    CHECK(total == doctest::Approx(5.0).epsilon(1e-9));
}
