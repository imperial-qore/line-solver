/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The api/sn layer: the sn_has_* / sn_is_* predicate family, the node-level
 * arrival and throughput tables, the station-to-station routing complement and
 * the algorithm-selection feature map.
 *
 * ORACLES.
 *  (a) The reference .m files read line by line; every assertion below names
 *      the MATLAB expression it reproduces.
 *  (b) Conservation, which is stronger than agreement on a number: the node
 *      throughput of a station row must equal its station throughput exactly,
 *      and a Router's node throughput must equal the flow the routing carries
 *      into it.
 *  (c) The stochastic complement identity: eliminating a Router from the
 *      stateful routing must leave the station-to-station routing a stochastic
 *      matrix whose entries are the path probabilities through it.
 *
 * The predicate assertions include the two that NetworkStruct's own member
 * functions get WRONG (see sn_predicates.h): a round-robin model is not
 * product form, and neither is a model with a live Fork.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_has_classdep_routing.h"
#include "line/api/sn/sn_node_metrics.h"
#include "line/api/sn/sn_predicates.h"
#include "line/api/sn/sn_print_routing_matrix.h"
#include "line/api/sn/sn_region_members.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/api/sn/sn_setters.h"
#include "line/lang/qn/network_builder.h"

namespace qn = line::qn;
namespace api = line::api;
using line::Matrix;
using line::lang::RoutingStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> PS Queue -> Delay, one closed class: the product-form baseline. */
qn::Network<double> closed_pf(double njobs) {
    qn::Network<double> m("pf");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The same but the Delay dispatches round-robin over two queues. */
qn::Network<double> closed_rrobin(double njobs) {
    qn::Network<double> m("rr");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 0.5);
    P.set(d, q2, 0.5);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    m.set_routing(d, c, RoutingStrategy::RROBIN);
    return m;
}

/** FCFS station serving two classes at DIFFERENT rates: heterogeneous FCFS. */
qn::Network<double> closed_het_fcfs() {
    qn::Network<double> m("het");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::exp_rate(2.0));
    m.set_service(q, c2, Dist::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Router -> Q1 (0.25) / Q2 (0.75) -> Delay: a stateful non-station. */
qn::Network<double> closed_router(double njobs) {
    qn::Network<double> m("router");
    const std::size_t d = m.add_delay("Think");
    const std::size_t rt = m.add_router("R");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(d, rt, 1.0);
    P.set(rt, q1, 0.25);
    P.set(rt, q2, 0.75);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("sn predicate family reproduces the reference expressions") {
    qn::Network<double> a = closed_pf(3.0);
    const qn::NetworkStruct<double>& sa = a.get_struct();

    // sn_has_ps / sn_has_inf / sn_has_fcfs: any(sn.sched == <strategy>)
    CHECK(api::sn_has_ps(sa));
    CHECK(api::sn_has_inf(sa));
    CHECK_FALSE(api::sn_has_fcfs(sa));
    CHECK_FALSE(api::sn_has_lcfs(sa));

    // the class and chain counts
    CHECK(api::sn_has_closed_classes(sa));
    CHECK_FALSE(api::sn_has_open_classes(sa));
    CHECK(api::sn_is_closed_model(sa));
    CHECK_FALSE(api::sn_is_open_model(sa));
    CHECK_FALSE(api::sn_is_mixed_model(sa));
    CHECK(api::sn_has_single_class(sa));
    CHECK(api::sn_has_single_chain(sa));
    CHECK_FALSE(api::sn_has_multiple_closed_classes(sa));
    CHECK_FALSE(api::sn_has_class_switching(sa));
    CHECK_FALSE(api::sn_has_priorities(sa));
    CHECK_FALSE(api::sn_has_multi_server(sa));
    CHECK_FALSE(api::sn_has_fractional_populations(sa));
    CHECK_FALSE(api::sn_has_load_dependence(sa));
    CHECK_FALSE(api::sn_has_fork_join(sa));
    CHECK_FALSE(api::sn_has_sd_routing(sa));

    // PS and INF only, no heterogeneous FCFS, no priorities, no fork, no
    // state-dependent routing: product form both ways
    CHECK(api::sn_has_product_form(sa));
    CHECK(api::sn_has_product_form_not_het_fcfs(sa));
    CHECK(api::sn_is_population_model(sa));
    CHECK_FALSE(api::sn_is_bas_model(sa));
    CHECK_FALSE(api::sn_is_mm1k_loss(sa));

    // sn_has_homogeneous_scheduling reduces to nstations == 1 in the reference
    CHECK_FALSE(api::sn_has_homogeneous_scheduling(sa, SchedStrategy::INF));
}

TEST_CASE("round-robin routing takes product form away") {
    qn::Network<double> b = closed_rrobin(3.0);
    const qn::NetworkStruct<double>& sb = b.get_struct();

    // sn_has_sd_routing: RROBIN is state dependent even though the refresh
    // spreads its probabilities uniformly
    CHECK(api::sn_has_sd_routing(sb));
    CHECK_FALSE(api::sn_has_product_form(sb));
    CHECK_FALSE(api::sn_has_product_form_not_het_fcfs(sb));

    // and this is precisely where NetworkStruct's own member DISAGREES with
    // the reference: it omits the state-dependent-routing conjunct, so it
    // reports a product form the model does not have
    CHECK(sb.has_product_form());
}

TEST_CASE("a binding buffer takes product form away, a non-binding one does not") {
    // cqn_bas_blocking: Queue1 blocks after service, Queue2 holds one job of the
    // two circulating, so the buffer BINDS and the truncation couples the two
    // station occupancies. Before sn_has_blocking existed, no conjunct of
    // sn_has_product_form read sn.cap/sn.classcap/sn.droprule and this model
    // reported hasProductForm=1.
    qn::Network<double> m("bas");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("Class1", 2.0, q1);
    m.set_service(q1, c, Dist::exp_rate(1.0));
    m.set_service(q2, c, Dist::exp_rate(0.8));
    m.set_capacity(q2, 1.0);
    m.set_drop_rule(q1, c, line::lang::DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sm = m.get_struct();
    CHECK(api::sn_has_blocking(sm));
    CHECK_FALSE(api::sn_has_product_form(sm));
    CHECK_FALSE(sm.has_product_form());  // the member carries the conjunct too
    CHECK(api::sn_is_bas_model(sm));

    // A capacity at least as large as the population that can reach the station
    // can never refuse a job, so the declaration is a no-op: the common
    // setCap(N) idiom on an N-job closed model stays product form.
    qn::Network<double> n = closed_pf(3.0);
    n.set_capacity(1, 3.0);  // the PS queue, node index 1
    const qn::NetworkStruct<double>& snn = n.get_struct();
    CHECK_FALSE(api::sn_has_blocking(snn));
    CHECK(api::sn_has_product_form(snn));
}

TEST_CASE("heterogeneous FCFS is detected and excludes product form") {
    qn::Network<double> d = closed_het_fcfs();
    const qn::NetworkStruct<double>& sd = d.get_struct();
    CHECK(api::sn_has_fcfs(sd));
    CHECK(api::sn_has_multi_class(sd));
    CHECK(api::sn_has_multi_class_fcfs(sd));
    CHECK(api::sn_has_multi_class_heter_fcfs(sd));
    // both rates exponential, so the "heterogeneous EXPONENTIAL FCFS" arm fires
    CHECK(api::sn_has_multi_class_heter_exp_fcfs(sd));
    CHECK_FALSE(api::sn_has_product_form(sd));
    // sn_has_product_form_not_het_fcfs asks BCMP type 1 in full: the FCFS
    // service must be exponential AND class-independent. Two exponentials
    // satisfy the SCV half, but the means here are 1/2 against 1/4 in TWO
    // separate chains, so the between-chain service-time spread takes the
    // product form away (the mean test is visit-weighted per chain; a spread
    // confined within one chain, or on zero-visit classes, does not divert).
    CHECK_FALSE(api::sn_has_product_form_not_het_fcfs(sd));
    // check_means=false is the bypass the class-dependent FCFS algorithms
    // (ab, schmidt, schmidt-ext) take, for which the exclusion is the point:
    // it restores the SCV-only test, which this model passes.
    CHECK(api::sn_has_product_form_not_het_fcfs(sd, false));
    CHECK(api::sn_has_multiple_closed_classes(sd));
    CHECK_FALSE(api::sn_is_population_model(sd));  // FCFS is not population only
}

TEST_CASE("sn_has_classdep_routing separates shared from per-class routing") {
    // one class: the reference returns false without inspecting anything
    qn::Network<double> a = closed_pf(3.0);
    { const qn::NetworkStruct<double>& s0 = a.get_struct(); CHECK_FALSE(api::sn_has_classdep_routing(s0)); }

    // two classes routed identically: every (i,j) block has one shared value
    qn::Network<double> d = closed_het_fcfs();
    { const qn::NetworkStruct<double>& s1 = d.get_struct(); CHECK_FALSE(api::sn_has_classdep_routing(s1)); }

    // two classes routed differently at the same station
    qn::Network<double> m("cdr");
    const std::size_t th = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 1.0, th);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, th);
    m.set_service(th, c1, Dist::exp_rate(1.0));
    m.set_service(th, c2, Dist::exp_rate(1.0));
    m.set_service(q1, c1, Dist::exp_rate(2.0));
    m.set_service(q1, c2, Dist::exp_rate(2.0));
    m.set_service(q2, c1, Dist::exp_rate(2.0));
    m.set_service(q2, c2, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, th, q1, 1.0);   // class 1 always to Q1
    P.set(c1, c1, q1, th, 1.0);
    P.set(c2, c2, th, q2, 1.0);   // class 2 always to Q2
    P.set(c2, c2, q2, th, 1.0);
    m.link(P);
    { const qn::NetworkStruct<double>& s2 = m.get_struct(); CHECK(api::sn_has_classdep_routing(s2)); }
}

TEST_CASE("sn_rt_stations eliminates a Router by stochastic complement") {
    qn::Network<double> m = closed_router(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.nstations == 3);
    // MATLAB's refreshStruct marks a Router STATEFUL, so there its `rt` is 4
    // stateful rows and sn_rt_stations genuinely eliminates one; this port's
    // builder creates a Router with stateful = false and folds it into `P` at
    // refresh time instead, so `rt` is already station-to-station and the
    // complement is the identity. Either way the station-to-station routing
    // below must be the path probabilities THROUGH the Router, which is what
    // the assertions check and what any caller of sn_rt_stations relies on.
    REQUIRE(sn.nof_stateful() == 3);

    const api::SnRtStations<double> r = api::sn_rt_stations(sn);
    REQUIRE(r.rtst.rows() == 3);  // 3 stations x 1 class

    // stations are Think(0), Q1(1), Q2(2) in creation order, minus the Router
    // Oracle (c): the path Think -> Router -> Q1 has probability 0.25 and
    // Think -> Router -> Q2 has 0.75, which is what eliminating R must give
    CHECK(r.rtst(0, 1) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(r.rtst(0, 2) == doctest::Approx(0.75).epsilon(1e-12));
    CHECK(r.rtst(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.rtst(2, 0) == doctest::Approx(1.0).epsilon(1e-12));
    for (std::size_t a = 0; a < 3; ++a) {
        double s = 0.0;
        for (std::size_t b = 0; b < 3; ++b) s += r.rtst(a, b);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
    // Vst is cellsum(sn.visits) at the station rows, normalised to 1 at the
    // reference station
    CHECK(r.Vst(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.Vst(1, 0) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(r.Vst(2, 0) == doctest::Approx(0.75).epsilon(1e-12));
}

TEST_CASE("node arrival and throughput tables agree with the station tables") {
    qn::Network<double> m = closed_router(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t I = sn.nodes.size(), R = sn.nclasses, M = sn.nstations;

    // a consistent station throughput: X = 1 at the reference station, scaled
    // by the visit ratios everywhere else
    Matrix<double> TN(M, R, 0.0);
    const double X = 0.8;
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
        double v = 0.0;
        for (std::size_t c = 0; c < sn.visits.size(); ++c) v += sn.visits[c](isf, 0);
        TN(i, 0) = X * v;
    }
    Matrix<double> AN(M, R, 0.0);
    for (std::size_t i = 0; i < M; ++i) AN(i, 0) = TN(i, 0);  // closed model, flow balance

    const Matrix<double> ANn = api::sn_get_node_arvr_from_tput(sn, TN, AN);
    REQUIRE(ANn.rows() == I);
    // Oracle (b): each station's node row repeats its station row
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t nd = sn.station_to_node[i];
        CHECK(ANn(nd - 1, 0) == doctest::Approx(AN(i, 0)).epsilon(1e-12));
    }
    // the Router is node 2 (1-based), and everything that leaves Think enters it
    CHECK(ANn(1, 0) == doctest::Approx(TN(0, 0)).epsilon(1e-12));

    const Matrix<double> TNn = api::sn_get_node_tput_from_tput(sn, TN, ANn);
    REQUIRE(TNn.rows() == I);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t nd = sn.station_to_node[i];
        CHECK(TNn(nd - 1, 0) == doctest::Approx(TN(i, 0)).epsilon(1e-12));
    }
    // the Router passes its whole arriving flow onward
    CHECK(TNn(1, 0) == doctest::Approx(TN(0, 0)).epsilon(1e-12));
}

TEST_CASE("sn_region_members prefers the declared membership") {
    qn::Network<double> m = closed_pf(2.0);
    m.get_struct();
    qn::NetworkStruct<double>& sn = m.raw_struct();
    CHECK(api::sn_region_members(sn, 0).empty());  // no such region

    typename qn::NetworkStruct<double>::Region rg;
    rg.name = "FCR";
    rg.members.assign(sn.nstations, false);
    rg.members[1] = true;
    rg.cap.assign(sn.nstations, std::vector<double>(sn.nclasses + 1, -1.0));
    rg.maxmem.assign(sn.nstations, -1.0);
    sn.regions.push_back(rg);
    const std::vector<bool> mask = api::sn_region_members(sn, 0);
    REQUIRE(mask.size() == sn.nstations);
    CHECK_FALSE(mask[0]);
    CHECK(mask[1]);
    const std::vector<std::size_t> lst = api::sn_region_member_list(sn, 0);
    REQUIRE(lst.size() == 1);
    CHECK(lst[0] == 2);

    // with no declared flags the -1 sentinel decides instead
    sn.regions[0].members.assign(sn.nstations, false);
    sn.regions[0].cap[0][0] = 3.0;
    const std::vector<bool> m2 = api::sn_region_members(sn, 0);
    CHECK(m2[0]);
    CHECK_FALSE(m2[1]);
}

TEST_CASE("sn_set_service and sn_refresh_process_fields rebuild the representation") {
    qn::Network<double> m = closed_pf(2.0);
    m.get_struct();
    qn::NetworkStruct<double>& sn = m.raw_struct();

    // exponential stays one phase
    api::sn_set_service(sn, 2, 1, 5.0, 1.0, true);
    CHECK(sn.rates(1, 0) == doctest::Approx(5.0));
    CHECK(sn.service[1][0].type == line::lang::ProcessType::EXP);
    CHECK(sn.phases_of(2, 1) == 1);

    // SCV 1/4 becomes Erlang-4, the reference's ceil(1/scv)
    api::sn_set_service(sn, 2, 1, 5.0, 0.25, true);
    CHECK(sn.service[1][0].type == line::lang::ProcessType::ERLANG);
    CHECK(sn.phases_of(2, 1) == 4);
    CHECK(line::num_traits<double>::to_double(sn.service[1][0].mean) ==
          doctest::Approx(0.2).epsilon(1e-12));

    // SCV 4 becomes a two-phase hyperexponential
    api::sn_set_service(sn, 2, 1, 5.0, 4.0, true);
    CHECK(sn.service[1][0].type == line::lang::ProcessType::HYPEREXP);
    CHECK(sn.phases_of(2, 1) == 2);

    api::sn_set_servers(sn, 2, 3.0);
    CHECK(api::sn_has_multi_server(sn));
    api::sn_set_priority(sn, 1, 2);
    CHECK(api::sn_has_priorities(sn));
    api::sn_set_population(sn, 1, 5.0, false);
    CHECK(sn.classes[0].population == doctest::Approx(5.0));
}

TEST_CASE("sn_print_routing_matrix renders every positive node edge") {
    qn::Network<double> m = closed_router(2.0);
    const qn::NetworkStruct<double>& sp = m.get_struct();
    const std::string s = api::sn_print_routing_matrix(sp);
    CHECK(s.find("Think [C1] => R [C1]") != std::string::npos);
    CHECK(s.find("R [C1] => Q1 [C1]") != std::string::npos);
    CHECK(s.find("R [C1] => Q2 [C1]") != std::string::npos);
    CHECK(s.find("Q1 [C1] => Think [C1]") != std::string::npos);
}
