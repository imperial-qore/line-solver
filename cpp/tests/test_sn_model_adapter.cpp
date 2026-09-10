/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The two class-level ModelAdapter transforms: chain aggregation and class
 * removal.
 *
 * THE ORACLE IS AN IDENTITY, NOT A STORED NUMBER, for both of them.
 *
 *  - Chain aggregation is EXACT on a product-form model, because MVA and
 *    convolution already solve such a model at the CHAIN level: the chain's
 *    per-station demand is all either of them reads. So the aggregate's
 *    per-station metrics must equal the ORIGINAL's metrics summed over the
 *    classes of that chain, to solver precision, and the aggregate must carry
 *    one class per chain rather than one per class. A golden would only say
 *    the code had not changed; this says it is right, and it is the same check
 *    the MATLAB and python wiring tests make.
 *
 *  - Class removal has an even cleaner oracle: the model without class r,
 *    built directly, is a model this transform must reproduce. Metrics from
 *    the two must agree exactly, since they are the same network.
 *
 * Each case below breaks a different piece if it is wrong:
 *  - a station visited by BOTH classes of a chain checks the alpha weighting;
 *    a chain whose classes have different rates there makes SCVchain leave 1,
 *    so the Erlang / HyperExp rungs of the refit ladder are exercised rather
 *    than the exponential shortcut;
 *  - class switching along the cycle checks that the synthesized ClassSwitch
 *    nodes are folded away by `rt` rather than carried into the aggregate;
 *  - the open chain checks that the Sink, which is not stateful, is folded by
 *    the same complement and the chain still routes;
 *  - C == K checks the copy branch, which must NOT refit the service laws;
 *  - removing a MIDDLE class checks the index shift in every table at once --
 *    a transform that only handled the last class would pass on class 3 and
 *    fail here;
 *  - the refusals check that a construct whose class indexing cannot be
 *    rewritten is named rather than silently mis-sliced.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_aggregate_chains.h"
#include "line/api/sn/sn_remove_class.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"

namespace qn = line::qn;
namespace api = line::api;
using line::InputError;
using line::UnsupportedError;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;
using Mat = line::Matrix<double>;

namespace {

line::mva::AvgResult<double> solve_nc(const qn::NetworkStruct<double>& sn) {
    line::nc::NcSolverOptions o;
    o.method = "exact";
    return line::nc::solver_nc_run_analyzer(sn, o);
}

line::mva::AvgResult<double> solve_mva(const qn::NetworkStruct<double>& sn) {
    line::mva::MvaOptions o;
    return line::mva::solver_mva_run_analyzer(sn, o, Mat());
}

/** The original's metric at station `ist` (1-based), summed over chain `c`. */
double chain_sum(const Mat& X, const qn::NetworkStruct<double>& sn, std::size_t c,
                 std::size_t ist) {
    double s = 0.0;
    for (std::size_t a = 0; a < sn.inchain[c].size(); ++a)
        s += X(ist - 1, sn.inchain[c][a] - 1);
    return s;
}

/**
 * A closed cycle Delay -> Q1 -> Q2 -> Delay whose two classes form ONE chain.
 *
 * The switch happens on the way out of Q1: half the departures leave as C2.
 * Both classes are therefore served at Q2 and at the Delay, which is what makes
 * the alpha weighting observable -- a model where each station saw one class
 * would pass with alpha ignored.
 */
qn::Network<double> cycle_two_class_one_chain(double njobs, double q1r1, double q1r2,
                                              double q2r1, double q2r2) {
    qn::Network<double> m("cyc2");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", njobs, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(q1, c1, Dist::exp_rate(q1r1));
    m.set_service(q1, c2, Dist::exp_rate(q1r2));
    m.set_service(q2, c1, Dist::exp_rate(q2r1));
    m.set_service(q2, c2, Dist::exp_rate(q2r2));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c2, c2, d, q1, 1.0);
    // Out of Q1 the chain splits its label; both halves go on to Q2.
    P.set(c1, c1, q1, q2, 0.5);
    P.set(c1, c2, q1, q2, 0.5);
    P.set(c2, c2, q1, q2, 1.0);
    P.set(c1, c1, q2, d, 1.0);
    P.set(c2, c1, q2, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("sn_aggregate_chains: one chain of two classes, exact on product form") {
    qn::Network<double> m = cycle_two_class_one_chain(6.0, 2.0, 3.0, 4.0, 5.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.nclasses == 2);
    REQUIRE(sn.nchains == 1);

    api::ChainAggregationResult<double> agg = api::sn_aggregate_chains(sn);
    CHECK(agg.deagg.isaggregated);
    CHECK(agg.deagg.nchains == 1);
    CHECK(agg.deagg.nclasses == 2);

    const qn::NetworkStruct<double>& sa = agg.model.get_struct();
    CHECK(sa.nclasses == 1);
    CHECK(sa.nchains == 1);
    // The three stations survive; the ClassSwitch `link` synthesized does not.
    CHECK(sa.nstations == 3);
    for (std::size_t i = 0; i < sa.nodes.size(); ++i)
        CHECK(sa.nodes[i].nodetype != qn::NodeType::ClassSwitch);

    const line::mva::AvgResult<double> ro = solve_nc(sn);
    const line::mva::AvgResult<double> ra = solve_nc(sa);

    for (std::size_t i = 1; i <= sn.nstations; ++i) {
        const std::size_t nd = agg.stationnode[i - 1];
        REQUIRE(nd != 0);
        const std::size_t ia = agg.model.station_index(nd);
        CHECK(ra.QN(ia - 1, 0) == doctest::Approx(chain_sum(ro.QN, sn, 0, i)).epsilon(1e-9));
        CHECK(ra.UN(ia - 1, 0) == doctest::Approx(chain_sum(ro.UN, sn, 0, i)).epsilon(1e-9));
        CHECK(ra.TN(ia - 1, 0) == doctest::Approx(chain_sum(ro.TN, sn, 0, i)).epsilon(1e-9));
    }
    // The population is conserved: the chain still holds every job.
    double qtot = 0.0;
    for (std::size_t i = 0; i < sa.nstations; ++i) qtot += ra.QN(i, 0);
    CHECK(qtot == doctest::Approx(6.0).epsilon(1e-9));
}

TEST_CASE("sn_aggregate_chains: the alpha shares are the class split of the chain visits") {
    qn::Network<double> m = cycle_two_class_one_chain(4.0, 2.0, 2.0, 3.0, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const api::ChainAggregationResult<double> agg = api::sn_aggregate_chains(sn);
    // Every station's shares sum to one over the chain's classes, which is what
    // makes a chain result splittable back over them.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        double s = 0.0;
        for (std::size_t a = 0; a < sn.inchain[0].size(); ++a)
            s += agg.alpha(i, sn.inchain[0][a] - 1);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-9));
    }
    // Q1 is entered only in C1, so its share there is the whole chain.
    const std::size_t q1 = 2;  // Think, Q1, Q2 in creation order
    CHECK(agg.alpha(q1 - 1, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(agg.alpha(q1 - 1, 1) == doctest::Approx(0.0).epsilon(1e-9));
    // Q2 is entered half in each, since the switch is on the Q1 -> Q2 link.
    const std::size_t q2 = 3;
    CHECK(agg.alpha(q2 - 1, 0) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(agg.alpha(q2 - 1, 1) == doctest::Approx(0.5).epsilon(1e-9));
}

TEST_CASE("sn_aggregate_chains: C == K copies the model instead of refitting it") {
    qn::Network<double> m("indep");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    // An Erlang, so a refit to (mean, SCV) would be visible in the struct.
    m.set_service(q, c1, Dist::erlang(4.0, 2));
    m.set_service(q, c2, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.nchains == sn.nclasses);

    api::ChainAggregationResult<double> agg = api::sn_aggregate_chains(sn);
    CHECK_FALSE(agg.deagg.isaggregated);
    const qn::NetworkStruct<double>& sa = agg.model.get_struct();
    CHECK(sa.nclasses == 2);
    CHECK(sa.nchains == 2);
    // The Erlang is still an Erlang: the copy branch must not two-moment it.
    const std::size_t iq = agg.model.station_index(agg.stationnode[m.station_index(q) - 1]);
    CHECK(sa.service[iq - 1][0].type == line::lang::ProcessType::ERLANG);
    CHECK(sa.scv(iq - 1, 0) == doctest::Approx(0.5).epsilon(1e-9));
}

TEST_CASE("sn_aggregate_chains: an open chain still routes once the Sink is folded") {
    qn::Network<double> m("open2");
    const std::size_t src = m.add_source("Source");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dist::exp_rate(0.5));
    m.set_arrival(src, c2, Dist::disabled_dist());
    m.set_service(q1, c1, Dist::exp_rate(2.0));
    m.set_service(q1, c2, Dist::exp_rate(2.0));
    m.set_service(q2, c1, Dist::exp_rate(3.0));
    m.set_service(q2, c2, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q1, 1.0);
    P.set(c1, c2, q1, q2, 1.0);
    P.set(c2, c2, q2, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.nclasses == 2);
    REQUIRE(sn.nchains == 1);

    api::ChainAggregationResult<double> agg = api::sn_aggregate_chains(sn);
    const qn::NetworkStruct<double>& sa = agg.model.get_struct();
    CHECK(sa.nclasses == 1);

    const line::mva::AvgResult<double> ra = solve_mva(sa);
    // lambda = 0.5 through both queues: U = lambda / mu, exactly.
    const std::size_t i1 = agg.model.station_index(agg.stationnode[m.station_index(q1) - 1]);
    const std::size_t i2 = agg.model.station_index(agg.stationnode[m.station_index(q2) - 1]);
    CHECK(ra.UN(i1 - 1, 0) == doctest::Approx(0.25).epsilon(1e-6));
    CHECK(ra.UN(i2 - 1, 0) == doctest::Approx(0.5 / 3.0).epsilon(1e-6));
    CHECK(ra.TN(i1 - 1, 0) == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(ra.TN(i2 - 1, 0) == doctest::Approx(0.5).epsilon(1e-6));
}

TEST_CASE("sn_aggregate_chains: a node kind it cannot carry is named, not dropped") {
    // A Logger is the reference's `otherwise` branch: it warns and SKIPS,
    // returning a model with the node gone and the routing rewired around it.
    // This port refuses instead, and the refusal has to be reached, so the two
    // classes here form ONE chain -- at C == K the copy branch runs and no node
    // sweep happens at all.
    qn::Network<double> m("withlogger");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t lg = m.add_logger("Log");
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::exp_rate(2.0));
    m.set_service(q, c2, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c1, c1, q, lg, 1.0);
    P.set(c2, c2, q, lg, 1.0);
    // The switch is on the way out of the Logger, which is what merges the two
    // classes into one chain.
    P.set(c1, c2, lg, d, 1.0);
    P.set(c2, c1, lg, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.nclasses == 2);
    REQUIRE(sn.nchains == 1);
    CHECK_THROWS_AS(api::sn_aggregate_chains(sn), UnsupportedError);
}

TEST_CASE("sn_remove_class: the transform reproduces the model built without the class") {
    // Three classes on two stations, and the MIDDLE one is removed: a transform
    // that only handled the last class would pass on C3 and fail here.
    qn::Network<double> m3("three");
    const std::size_t d3 = m3.add_delay("Think");
    const std::size_t q3 = m3.add_queue("Q1", SchedStrategy::PS);
    const std::size_t a1 = m3.add_closed_class("C1", 3.0, d3);
    const std::size_t a2 = m3.add_closed_class("C2", 2.0, d3);
    const std::size_t a3 = m3.add_closed_class("C3", 4.0, d3);
    m3.set_service(d3, a1, Dist::exp_rate(1.0));
    m3.set_service(d3, a2, Dist::exp_rate(2.0));
    m3.set_service(d3, a3, Dist::exp_rate(3.0));
    m3.set_service(q3, a1, Dist::exp_rate(4.0));
    m3.set_service(q3, a2, Dist::exp_rate(5.0));
    m3.set_service(q3, a3, Dist::exp_rate(6.0));
    qn::RoutingMatrix<double> P3;
    P3.set(a1, a1, d3, q3, 1.0);
    P3.set(a1, a1, q3, d3, 1.0);
    P3.set(a2, a2, d3, q3, 1.0);
    P3.set(a2, a2, q3, d3, 1.0);
    P3.set(a3, a3, d3, q3, 1.0);
    P3.set(a3, a3, q3, d3, 1.0);
    m3.link(P3);
    const qn::NetworkStruct<double> sn3 = m3.get_struct();
    REQUIRE(sn3.nclasses == 3);

    const qn::NetworkStruct<double> cut = api::sn_remove_class(sn3, std::size_t(2));
    CHECK(cut.nclasses == 2);
    CHECK(cut.classes[0].name == "C1");
    CHECK(cut.classes[1].name == "C3");
    // The tables shifted with the class list, not around it.
    CHECK(cut.classes[0].population == doctest::Approx(3.0));
    CHECK(cut.classes[1].population == doctest::Approx(4.0));

    // The same model, built directly.
    qn::Network<double> m2("two");
    const std::size_t d2 = m2.add_delay("Think");
    const std::size_t q2 = m2.add_queue("Q1", SchedStrategy::PS);
    const std::size_t b1 = m2.add_closed_class("C1", 3.0, d2);
    const std::size_t b3 = m2.add_closed_class("C3", 4.0, d2);
    m2.set_service(d2, b1, Dist::exp_rate(1.0));
    m2.set_service(d2, b3, Dist::exp_rate(3.0));
    m2.set_service(q2, b1, Dist::exp_rate(4.0));
    m2.set_service(q2, b3, Dist::exp_rate(6.0));
    qn::RoutingMatrix<double> P2;
    P2.set(b1, b1, d2, q2, 1.0);
    P2.set(b1, b1, q2, d2, 1.0);
    P2.set(b3, b3, d2, q2, 1.0);
    P2.set(b3, b3, q2, d2, 1.0);
    m2.link(P2);
    const qn::NetworkStruct<double> sn2 = m2.get_struct();

    const line::mva::AvgResult<double> rc = solve_nc(cut);
    const line::mva::AvgResult<double> rd = solve_nc(sn2);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(rc.QN(i, r) == doctest::Approx(rd.QN(i, r)).epsilon(1e-9));
            CHECK(rc.UN(i, r) == doctest::Approx(rd.UN(i, r)).epsilon(1e-9));
            CHECK(rc.TN(i, r) == doctest::Approx(rd.TN(i, r)).epsilon(1e-9));
        }
    // The rates really moved: C3's row is now at index 1, not 2.
    CHECK(cut.rates(1, 1) == doctest::Approx(6.0).epsilon(1e-12));
}

TEST_CASE("sn_remove_class: the input is left untouched, and the name resolves") {
    qn::Network<double> m("three");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
    const std::size_t c3 = m.add_closed_class("C3", 4.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(2.0));
    m.set_service(d, c3, Dist::exp_rate(3.0));
    m.set_service(q, c1, Dist::exp_rate(4.0));
    m.set_service(q, c2, Dist::exp_rate(5.0));
    m.set_service(q, c3, Dist::exp_rate(6.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    P.set(c3, c3, d, q, 1.0);
    P.set(c3, c3, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();

    const qn::NetworkStruct<double> byidx = api::sn_remove_class(sn, std::size_t(2));
    const qn::NetworkStruct<double> byname = api::sn_remove_class(sn, std::string("C2"));
    CHECK(byidx.nclasses == byname.nclasses);
    for (std::size_t r = 0; r < byidx.nclasses; ++r)
        CHECK(byidx.classes[r].name == byname.classes[r].name);
    // The non-mutating contract: the source struct still has all three.
    CHECK(sn.nclasses == 3);
    CHECK(sn.classes.size() == 3);

    CHECK_THROWS_AS(api::sn_remove_class(sn, std::size_t(0)), InputError);
    CHECK_THROWS_AS(api::sn_remove_class(sn, std::size_t(4)), InputError);
    CHECK_THROWS_AS(api::sn_remove_class(sn, std::string("nosuch")), InputError);
}

TEST_CASE("sn_remove_class: the last class of a model cannot be removed") {
    qn::Network<double> m("one");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    CHECK_THROWS_AS(api::sn_remove_class(sn, std::size_t(1)), InputError);
}

TEST_CASE("sn_remove_class: the class-switching matrix is sliced with the class list") {
    // C1 -> C2 at a ClassSwitch node, plus a third class that ignores it.
    qn::Network<double> m("cs3");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 0.0, d);
    const std::size_t c3 = m.add_closed_class("C3", 3.0, d);
    Mat C(3, 3, 0.0);
    C(0, 1) = 1.0;  // C1 -> C2
    C(1, 0) = 1.0;  // C2 -> C1
    C(2, 2) = 1.0;  // C3 stays
    const std::size_t cs = m.add_class_switch("CS", C);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(d, c3, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::exp_rate(4.0));
    m.set_service(q, c2, Dist::exp_rate(4.0));
    m.set_service(q, c3, Dist::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 1; r <= 3; ++r) {
        P.set(r, r, d, q, 1.0);
        P.set(r, r, q, cs, 1.0);
        P.set(r, r, cs, d, 1.0);
    }
    m.link(P);
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.nclasses == 3);

    // Removing C3, which the switch leaves alone, must leave the C1 <-> C2
    // block intact at its new width.
    const qn::NetworkStruct<double> cut = api::sn_remove_class(sn, std::size_t(3));
    CHECK(cut.nclasses == 2);
    const Mat& Cc = cut.csmatrix.at(cs);
    CHECK(Cc.rows() == 2);
    CHECK(Cc.cols() == 2);
    CHECK(Cc(0, 1) == doctest::Approx(1.0));
    CHECK(Cc(1, 0) == doctest::Approx(1.0));
    CHECK(cut.nchains == 1);
}

TEST_CASE("sn_remove_class: a construct whose class indexing cannot be sliced is refused") {
    SUBCASE("a Fork") {
        qn::Network<double> m("fj");
        const std::size_t d = m.add_delay("Think");
        const std::size_t f = m.add_fork("F");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
        const std::size_t j = m.add_join("J", f);
        const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
        const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
        m.set_service(d, c1, Dist::exp_rate(1.0));
        m.set_service(d, c2, Dist::exp_rate(1.0));
        m.set_service(q1, c1, Dist::exp_rate(2.0));
        m.set_service(q1, c2, Dist::exp_rate(2.0));
        m.set_service(q2, c1, Dist::exp_rate(2.0));
        m.set_service(q2, c2, Dist::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        for (std::size_t r = 1; r <= 2; ++r) {
            P.set(r, r, d, f, 1.0);
            P.set(r, r, f, q1, 1.0);
            P.set(r, r, f, q2, 1.0);
            P.set(r, r, q1, j, 1.0);
            P.set(r, r, q2, j, 1.0);
            P.set(r, r, j, d, 1.0);
        }
        m.link(P);
        const qn::NetworkStruct<double> sn = m.get_struct();
        CHECK_THROWS_AS(api::sn_remove_class(sn, std::size_t(1)), UnsupportedError);
    }
    SUBCASE("a class-dependent scaling function") {
        qn::Network<double> m("cd");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
        const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
        m.set_service(d, c1, Dist::exp_rate(1.0));
        m.set_service(d, c2, Dist::exp_rate(1.0));
        m.set_service(q, c1, Dist::exp_rate(2.0));
        m.set_service(q, c2, Dist::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        for (std::size_t r = 1; r <= 2; ++r) {
            P.set(r, r, d, q, 1.0);
            P.set(r, r, q, d, 1.0);
        }
        m.link(P);
        const line::lang::CdScaling<double> beta =
            [](const std::vector<double>& n) { return std::vector<double>(n.size(), 1.0); };
        m.set_class_dependence(q, beta, std::vector<double>(2, 1.0));
        const qn::NetworkStruct<double> sn = m.get_struct();
        CHECK_THROWS_AS(api::sn_remove_class(sn, std::size_t(1)), UnsupportedError);
    }
}
