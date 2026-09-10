/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `solver_mam_traffic` and `solver_mam_traffic_mmap`: the per-link traffic
 * descriptors of an open network.
 *
 * WHAT CAN AND CANNOT BE ASSERTED. The descriptors are APPROXIMATIONS of the
 * true superposed processes -- merging n marked MAPs exactly costs the product
 * of their orders, so the reference compresses back to an APH(2) mixture -- so
 * nothing here is compared against an exact solver, and no correlation or
 * higher moment is checked. What IS exact, and is what these tests pin, is the
 * first-moment algebra around the approximation:
 *
 *  1. SPLITTING IS LINEAR in the routing probabilities. A flow split p / (1-p)
 *     must produce two links whose rates are exactly p lambda and (1-p) lambda,
 *     and their sum must be exactly lambda. `npfqn_traffic_split_cs` performs
 *     only additions and multiplications by those probabilities, so this holds
 *     entry by entry and not to within a tolerance.
 *  2. A CLASS SWITCH IS A PERMUTATION OF THE MARKINGS. Hiding a ClassSwitch
 *     node by the stochastic complement must move the class-r rate onto the
 *     class the switch sends it to, exactly.
 *  3. MERGING POISSON GIVES POISSON of the summed rate. The Kronecker sum of
 *     two exponential MAPs IS a Poisson process of the summed rate, and the
 *     APH(2) fit of an exponential returns the exponential, so split-then-merge
 *     must recover the input rate. This one goes through the compression, so it
 *     is asserted to floating tolerance rather than exactly.
 *  4. THE JOIN SYNCHRONIZATION HAS A CLOSED FORM for Poisson branches. With
 *     both branches Poisson, `mmap_max`'s phase process is the birth-death walk
 *     of the lead over -k..k, whose stationary law is geometric in the rate
 *     ratio, so the join's throughput is known in closed form and is strictly
 *     below the slower branch: the sync queue blocks. Two cases are computed by
 *     hand below and neither comes from running this code.
 *
 * Every model is built through `network_builder`, so the routing, the chains
 * and `rtnodes` are the ones the rest of the port sees.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_traffic.h"

using namespace line;
using Dd = lang::Distrib<double>;
using lang::SchedStrategy;
using mam::DepTable;
using mam::Mmap;
using mam::TrafficConfig;

namespace {

/** A DEP table of the right shape, every cell empty (MATLAB's `[]`). */
DepTable<double> empty_dep(std::size_t rows, std::size_t nclasses) {
    return DepTable<double>(rows, std::vector<mam::Map<double>>(nclasses));
}

/** Poisson departures of class `cls` (1-based) from the station behind `node`. */
void dep_at_station(DepTable<double>& DEP, const qn::NetworkStruct<double>& sn, std::size_t node,
                    std::size_t cls, double rate) {
    const std::size_t ist = sn.nodes[node - 1].station;
    REQUIRE(ist > 0);
    DEP[ist - 1][cls - 1] = mam::map_exponential(rate);
}

/** The same, for the node-indexed table the fork-join variant takes. */
void dep_at_node(DepTable<double>& DEP, std::size_t node, std::size_t cls, double rate) {
    DEP[node - 1][cls - 1] = mam::map_exponential(rate);
}

double rate_of(const Mmap<double>& m) {
    double s = 0.0;
    for (double v : mam::mmap_lambda(m)) s += v;
    return s;
}

/** A single-class MMAP for Poisson(rate), the shape mmap_max is fed. */
Mmap<double> poisson_mmap(double rate) {
    Mmap<double> m;
    m.D0 = Matrix<double>(1, 1, -rate);
    m.D1 = Matrix<double>(1, 1, rate);
    m.Dc.assign(1, Matrix<double>(1, 1, rate));
    return m;
}

template <class F>
void refuses(F f, const std::string& needle) {
    try {
        f();
        FAIL("expected a refusal naming: ", needle);
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find(needle) != std::string::npos);
    }
}

}  // namespace

TEST_CASE("mam traffic: a tandem link carries the departure rate unchanged") {
    const double lambda = 0.75;
    qn::Network<double> m("trafA");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(lambda));
    m.set_service(q, c, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c, lambda);
    dep_at_station(DEP, sn, q, c, lambda);  // flow balance at the queue

    const std::vector<Mmap<double>> ARV = mam::solver_mam_traffic(sn, DEP, TrafficConfig());
    REQUIRE(ARV.size() == sn.nof_nodes());
    // a Source has no arrival stream to describe; MATLAB leaves it []
    CHECK(ARV[src - 1].order() == 0);
    // a single incoming link is passed through, so the rate is the source's
    CHECK(rate_of(ARV[q - 1]) == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(rate_of(ARV[k - 1]) == doctest::Approx(lambda).epsilon(1e-12));
}

TEST_CASE("mam traffic: a probabilistic split conserves the first moment exactly") {
    const double lambda = 2.0, p = 0.25;
    qn::Network<double> m("trafB");
    const std::size_t src = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(lambda));
    m.set_service(q1, c, Dd::exp_rate(4.0));
    m.set_service(q2, c, Dd::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, p);
    P.set(src, q2, 1.0 - p);
    P.set(q1, k, 1.0);
    P.set(q2, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c, lambda);
    dep_at_station(DEP, sn, q1, c, p * lambda);
    dep_at_station(DEP, sn, q2, c, (1.0 - p) * lambda);

    const std::vector<Mmap<double>> ARV = mam::solver_mam_traffic(sn, DEP, TrafficConfig());
    const double r1 = rate_of(ARV[q1 - 1]), r2 = rate_of(ARV[q2 - 1]);
    CHECK(r1 == doctest::Approx(p * lambda).epsilon(1e-12));
    CHECK(r2 == doctest::Approx((1.0 - p) * lambda).epsilon(1e-12));
    CHECK(r1 + r2 == doctest::Approx(lambda).epsilon(1e-12));
}

TEST_CASE("mam traffic: splitting a Poisson stream and merging it back recovers its rate") {
    const double lambda = 2.0, p = 0.25;
    qn::Network<double> m("trafC");
    const std::size_t src = m.add_source("Src");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(lambda));
    m.set_service(q1, c, Dd::exp_rate(4.0));
    m.set_service(q2, c, Dd::exp_rate(4.0));
    m.set_service(q3, c, Dd::exp_rate(8.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q1, p);
    P.set(src, q2, 1.0 - p);
    P.set(q1, q3, 1.0);
    P.set(q2, q3, 1.0);
    P.set(q3, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c, lambda);
    dep_at_station(DEP, sn, q1, c, p * lambda);
    dep_at_station(DEP, sn, q2, c, (1.0 - p) * lambda);
    dep_at_station(DEP, sn, q3, c, lambda);

    const std::vector<Mmap<double>> ARV = mam::solver_mam_traffic(sn, DEP, TrafficConfig());
    // Two Poisson flows merge into a Poisson of the summed rate; the APH(2) fit
    // the compression applies is exact for an exponential, so the rate survives
    // to floating tolerance. The ORDER does not survive and is not checked.
    CHECK(rate_of(ARV[q3 - 1]) == doctest::Approx(lambda).epsilon(1e-9));
}

TEST_CASE("mam traffic: a ClassSwitch is hidden and its permutation rides the links") {
    const double l1 = 1.0, l2 = 3.0;
    qn::Network<double> m("trafD");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    Matrix<double> C(2, 2, 0.0);
    C(0, 1) = 1.0;  // class 1 leaves as class 2
    C(1, 0) = 1.0;  // and class 2 as class 1
    // the switching matrix is indexed by the classes, so the node comes after them
    const std::size_t cs = m.add_class_switch("CS", C);
    m.set_arrival(src, c1, Dd::exp_rate(l1));
    m.set_arrival(src, c2, Dd::exp_rate(l2));
    m.set_service(q, c1, Dd::exp_rate(8.0));
    m.set_service(q, c2, Dd::exp_rate(8.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, cs, 1.0);
    P.set(c2, c2, src, cs, 1.0);
    P.set(c1, c1, cs, q, 1.0);
    P.set(c2, c2, cs, q, 1.0);
    P.set(c1, c1, q, k, 1.0);
    P.set(c2, c2, q, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c1, l1);
    dep_at_station(DEP, sn, src, c2, l2);
    dep_at_station(DEP, sn, q, c1, l2);
    dep_at_station(DEP, sn, q, c2, l1);

    const std::vector<Mmap<double>> ARV = mam::solver_mam_traffic(sn, DEP, TrafficConfig());
    // the switch node itself carries no arrival descriptor: it is eliminated
    CHECK(ARV[cs - 1].order() == 0);
    const std::vector<double> lk = mam::mmap_lambda(ARV[q - 1]);
    REQUIRE(lk.size() == 2);
    CHECK(lk[0] == doctest::Approx(l2).epsilon(1e-12));
    CHECK(lk[1] == doctest::Approx(l1).epsilon(1e-12));
}

TEST_CASE("mam traffic: a node the reference's switch cannot describe is refused by name") {
    qn::Network<double> m("trafE");
    const std::size_t src = m.add_source("Src");
    const std::size_t rt = m.add_router("Router");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(1.0));
    m.set_service(q, c, Dd::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(src, rt, 1.0);
    P.set(rt, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c, 1.0);
    dep_at_station(DEP, sn, q, c, 1.0);
    refuses([&] { mam::solver_mam_traffic(sn, DEP, TrafficConfig()); }, "Router");
}

TEST_CASE("mam traffic: a Fork-Join model is refused by the plain traffic step") {
    qn::Network<double> m("trafF");
    const std::size_t src = m.add_source("Src");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(0.5));
    m.set_service(q1, c, Dd::exp_rate(1.0));
    m.set_service(q2, c, Dd::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    DepTable<double> DEP = empty_dep(sn.nstations, sn.nclasses);
    dep_at_station(DEP, sn, src, c, 0.5);
    refuses([&] { mam::solver_mam_traffic(sn, DEP, TrafficConfig()); }, "solver_mam_traffic_mmap");
}

/**
 * The oracles below are derived, not measured. With both branches Poisson the
 * phase process of mmap_max is the lead of stream a over stream b, a birth-death
 * walk on -k..k with up-rate mu_a and down-rate mu_b, blocked at both ends.
 * Detailed balance gives pi_s proportional to (mu_a/mu_b)^s, and an arrival is
 * emitted on every transition towards 0, so
 *   throughput = mu_b P(s > 0) + mu_a P(s < 0).
 * Symmetric case mu_a = mu_b = mu: pi is uniform on 2k+1 states, so the
 * throughput is mu 2k/(2k+1), which is 4 mu / 5 at k = 2.
 * Asymmetric case mu_a = 1, mu_b = 2, k = 1: pi is proportional to (2, 1, 1/2)
 * over s = -1, 0, 1, so the throughput is (2 * 1 + 0.5 * 2) / 3.5 = 6/7.
 */
TEST_CASE("mam traffic: mmap_max reaches the closed-form throughput of a Poisson join") {
    const double mu = 1.5;
    const Mmap<double> s2 =
        mam::mmap_normalize(mam::mmap_max(poisson_mmap(mu), poisson_mmap(mu), 2));
    CHECK(s2.order() == 5);  // na * nb * (1 + 2k)
    CHECK(rate_of(s2) == doctest::Approx(mu * 4.0 / 5.0).epsilon(1e-10));
    // the synchronization queue blocks, so the join is strictly slower
    CHECK(rate_of(s2) < mu);

    const Mmap<double> s1 =
        mam::mmap_max(poisson_mmap(1.0), poisson_mmap(2.0), 1);
    CHECK(s1.order() == 3);
    CHECK(rate_of(s1) == doctest::Approx(6.0 / 7.0).epsilon(1e-10));

    // D0 + D1 must be a generator without any normalization, and the class
    // matrices must partition D1; both are identities of the construction.
    for (std::size_t i = 0; i < s1.order(); ++i) {
        double row = 0.0;
        for (std::size_t jj = 0; jj < s1.order(); ++jj) row += s1.D0(i, jj) + s1.D1(i, jj);
        CHECK(row == doctest::Approx(0.0).epsilon(1e-12));
        for (std::size_t jj = 0; jj < s1.order(); ++jj) CHECK(s1.D1(i, jj) == s1.Dc[0](i, jj));
    }
}

TEST_CASE("mam traffic: the sync map groups the parallel branches of a fork-join") {
    qn::Network<double> m("trafG");
    const std::size_t src = m.add_source("Src");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(0.5));
    m.set_service(q1, c, Dd::exp_rate(1.0));
    m.set_service(q2, c, Dd::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mam::FjSyncMap fs = mam::sn_build_fj_sync_map(sn);
    REQUIRE(fs.ngroups == 1);
    CHECK(fs.fork_of_group[0] == f);
    CHECK(fs.join_of_group[0] == j);
    // both branches sit on the parallel path, the source does not
    CHECK(fs.node_sync[j - 1][q1 - 1] == 1);
    CHECK(fs.node_sync[j - 1][q2 - 1] == 1);
    CHECK(fs.node_sync[j - 1][src - 1] == 0);
    CHECK(fs.node_sync[q1 - 1][f - 1] == 0);
}

TEST_CASE("mam traffic: the fork-join variant synchronizes the branches at the join") {
    const double mu = 1.0, lambda = 0.4;
    qn::Network<double> m("trafH");
    const std::size_t src = m.add_source("Src");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dd::exp_rate(lambda));
    m.set_service(q1, c, Dd::exp_rate(mu));
    m.set_service(q2, c, Dd::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(src, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, k, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    // the fork-join variant indexes DEP by NODE, not by station
    DepTable<double> DEP = empty_dep(sn.nof_nodes(), sn.nclasses);
    dep_at_node(DEP, src, c, lambda);
    dep_at_node(DEP, f, c, lambda);
    dep_at_node(DEP, q1, c, mu);
    dep_at_node(DEP, q2, c, mu);
    dep_at_node(DEP, j, c, lambda);

    const mam::FjSyncMap fs = mam::sn_build_fj_sync_map(sn);
    TrafficConfig cfg;
    cfg.fj_sync_q_len = 2;
    const std::vector<Mmap<double>> ARV = mam::solver_mam_traffic_mmap(sn, DEP, cfg, fs);
    REQUIRE(ARV.size() == sn.nof_nodes());

    // Both branches carry Poisson(mu) into the join along a single edge each, so
    // the arrival descriptor at the join is exactly the synchronization of two
    // Poisson streams and its rate is the closed form above: mu 2k/(2k+1).
    CHECK(rate_of(ARV[j - 1]) == doctest::Approx(mu * 4.0 / 5.0).epsilon(1e-10));
    // superposing the branches instead would have given 2 mu, twice the join's
    // true throughput; the point of the FJ variant is that it does not
    CHECK(rate_of(ARV[j - 1]) < 2.0 * mu);
}
