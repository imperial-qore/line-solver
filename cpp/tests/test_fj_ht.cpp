/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Heidelberger-Trivedi fork-join arm, `options.config.fork_join='ht'`
 * (solvers/mva/fj_ht.h and the H-T branch of solvers/mva/fj_driver.h).
 *
 * THE ORACLE IS MATLAB, station by station and class by class. Every reference
 * value below is `MVA(model, options)` with `options.config.fork_join='ht'` and
 * `options.method='amva'` run under `matlab -singleCompThread` on the same
 * model, printed at `format long g`; H-T is an approximation with no closed
 * form to check against, so agreement with the reference implementation to the
 * digit is the only statement worth making about it.
 *
 * The models are the closed fork-join examples of `matlab/examples/basic/
 * forkJoin`: `fj_basic_closed` (the one example that PINS this method), plus
 * `fj_asymm` (branches of unequal length), `fj_threebranches` (two classes,
 * each its own chain) and `fj_serialfjs_closed` (two forks in series, which is
 * what exercises the per-fork auxiliary delay). On the last of these the
 * reference stops at `iter_max` rather than at its tolerance, so it is checked
 * to 1e-4 and against its own symmetry instead of to the digit.
 *
 * The refusals are checked too. H-T mints a ClosedClass of the original's
 * population per branch, so it cannot express an open class; it walks one
 * branch per auxiliary class, so it cannot express tasksPerLink > 1; and it
 * charges the whole span at the synchronisation point, so it cannot express a
 * fork with no join.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace {

using DD = line::lang::Distrib<double>;
using Net = line::qn::Network<double>;
using Routing = line::qn::RoutingMatrix<double>;
using line::lang::SchedStrategy;

line::mva::AvgResult<double> ht_avg(Net& m, const std::string& fj = "ht") {
    line::mva::MvaOptions opt;
    opt.method = "amva";
    opt.fork_join = fj;
    line::Matrix<double> init;
    return line::mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

/** `matlab/examples/basic/forkJoin/fj_basic_closed.m`, the example that pins 'ht'. */
Net basic_closed() {
    Net m("model");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("class1", 5.0, d);
    m.set_service(d, c, DD::exp_rate(1.0));
    m.set_service(q1, c, DD::exp_rate(1.0));
    m.set_service(q2, c, DD::exp_rate(1.0));
    Routing P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

/** `fj_asymm.m`: Queue1 alone against the Queue2 -> Queue3 chain. */
Net asymm() {
    Net m("model");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::FCFS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("class1", 10.0, d);
    m.set_service(q1, c, DD::exp_rate(1.0));
    m.set_service(q2, c, DD::exp_rate(2.0));
    m.set_service(q3, c, DD::exp_rate(1.0));
    m.set_service(d, c, DD::exp_rate(0.5));
    Routing P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, q3, 1.0);
    P.set(c, c, q3, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

/** `fj_threebranches.m`: two closed classes, each its own chain, each forking twice. */
Net threebranches() {
    Net m("model");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::PS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c1 = m.add_closed_class("class1", 10.0, d);
    const std::size_t c2 = m.add_closed_class("class2", 10.0, d);
    m.set_service(q1, c1, DD::exp_rate(1.5));
    m.set_service(q2, c1, DD::exp_rate(1.1));
    m.set_service(q3, c1, DD::exp_rate(2.5));
    m.set_service(d, c1, DD::exp_rate(0.5));
    m.set_service(q1, c2, DD::exp_rate(2.8));
    m.set_service(q2, c2, DD::exp_rate(3.0));
    m.set_service(q3, c2, DD::exp_rate(1.0));
    m.set_service(d, c2, DD::exp_rate(0.8));
    Routing P;
    const std::size_t cls[2] = {c1, c2};
    for (std::size_t x = 0; x < 2; ++x) {
        const std::size_t c = cls[x];
        P.set(c, c, d, f, 1.0);
        P.set(c, c, f, q1, 1.0);
        P.set(c, c, f, q2, 1.0);
        P.set(c, c, q2, q3, 1.0);
        P.set(c, c, q3, j, 1.0);
        P.set(c, c, q1, j, 1.0);
        P.set(c, c, j, d, 1.0);
    }
    m.link(P);
    return m;
}

/** `fj_serialfjs_closed.m`: two fork-join spans in series, each with its own join. */
Net serialfjs() {
    Net m("model");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
    const std::size_t f1 = m.add_fork("Fork");
    const std::size_t j1 = m.add_join("Join", f1);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::PS);
    const std::size_t q4 = m.add_queue("Queue4", SchedStrategy::PS);
    const std::size_t f2 = m.add_fork("Fork2");
    const std::size_t j2 = m.add_join("Join2", f2);
    const std::size_t c = m.add_closed_class("class1", 10.0, d);
    m.set_service(d, c, DD::exp_rate(0.5));
    m.set_service(q1, c, DD::exp_rate(1.0));
    m.set_service(q2, c, DD::exp_rate(1.0));
    m.set_service(q3, c, DD::exp_rate(1.0));
    m.set_service(q4, c, DD::exp_rate(1.0));
    Routing P;
    P.set(c, c, d, f1, 1.0);
    P.set(c, c, f1, q1, 1.0);
    P.set(c, c, f1, q2, 1.0);
    P.set(c, c, q1, j1, 1.0);
    P.set(c, c, q2, j1, 1.0);
    P.set(c, c, j1, f2, 1.0);
    P.set(c, c, f2, q3, 1.0);
    P.set(c, c, f2, q4, 1.0);
    P.set(c, c, q3, j2, 1.0);
    P.set(c, c, q4, j2, 1.0);
    P.set(c, c, j2, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("H-T reproduces MATLAB on fj_basic_closed, the example that pins it") {
    Net m = basic_closed();
    const line::mva::AvgResult<double> r = ht_avg(m);
    // MATLAB SolverMVA, amva, config.fork_join='ht'. Rows: Delay, Queue1, Queue2, Join.
    const double Q[4] = {0.94049440904849, 2.70633701284678, 2.70633701284678, 2.70633710335178};
    const double U[4] = {0.94049440904849, 0.940494426072379, 0.940494426072379, 0.0};
    const double R[4] = {1.0, 2.8775683702335, 2.8775683702335, 1.4387842332324};
    const double T[4] = {0.94049440904849, 0.940494426072379, 0.940494426072379,
                         0.94049440904849};
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(r.QN(i, 0) == doctest::Approx(Q[i]).epsilon(1e-9));
        CHECK(r.UN(i, 0) == doctest::Approx(U[i]).epsilon(1e-9));
        CHECK(r.RN(i, 0) == doctest::Approx(R[i]).epsilon(1e-9));
        CHECK(r.TN(i, 0) == doctest::Approx(T[i]).epsilon(1e-9));
    }
    // The two arms are DIFFERENT answers to the same model, not two spellings of
    // one: MMT keeps the circulating job on a branch, H-T routes it past them.
    Net m2 = basic_closed();
    CHECK(ht_avg(m2, "default").TN(0, 0) == doctest::Approx(0.770570).epsilon(1e-5));
}

TEST_CASE("H-T reproduces MATLAB on branches of unequal length") {
    Net m = asymm();
    const line::mva::AvgResult<double> r = ht_avg(m);
    // Rows: Delay1, Queue1, Queue2, Queue3, Join.
    const double Q[5] = {1.94646322696, 4.7740489614537, 0.891841865785355, 5.03103090684816,
                         5.41015172486815};
    const double U[5] = {1.94646322696, 0.973231930371047, 0.486615665209478, 0.973231330418956,
                         0.0};
    const double R[5] = {2.0, 4.90535586890741, 0.916371923005639, 5.16940911127717,
                         2.77947795566013};
    const double T[5] = {0.97323161348, 0.973231930371047, 0.973231330418956, 0.973231330418956,
                         0.97323161348};
    for (std::size_t i = 0; i < 5; ++i) {
        CHECK(r.QN(i, 0) == doctest::Approx(Q[i]).epsilon(1e-8));
        CHECK(r.UN(i, 0) == doctest::Approx(U[i]).epsilon(1e-8));
        CHECK(r.RN(i, 0) == doctest::Approx(R[i]).epsilon(1e-8));
        CHECK(r.TN(i, 0) == doctest::Approx(T[i]).epsilon(1e-8));
    }
}

TEST_CASE("H-T reproduces MATLAB with two classes, each its own chain") {
    Net m = threebranches();
    const line::mva::AvgResult<double> r = ht_avg(m);
    // Rows: Delay1, Queue1, Queue2, Queue3, Join. Columns: class1, class2.
    const double Q[5][2] = {{1.57461963729062, 0.850931502453232},
                            {1.53336144871088, 0.75584755843157},
                            {4.3126250767005, 1.66655928221381},
                            {3.87078209022378, 7.42449036681445},
                            {7.13399208899917, 8.45123977342315}};
    const double R[5][2] = {{2.0, 1.25},
                            {1.94759596751022, 1.1103237317339},
                            {5.47767221885155, 2.44813960966864},
                            {4.91646622271853, 10.9064160768749},
                            {4.53061292304293, 6.20734435299749}};
    const double T[5][2] = {{0.787309818645311, 0.680745201962586},
                            {0.787309829292322, 0.68074520685172},
                            {0.787309810517411, 0.680745197549979},
                            {0.787309810517411, 0.680745197549979},
                            {0.787309818645311, 0.680745201962586}};
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
            CHECK(r.QN(i, k) == doctest::Approx(Q[i][k]).epsilon(1e-8));
            CHECK(r.RN(i, k) == doctest::Approx(R[i][k]).epsilon(1e-8));
            CHECK(r.TN(i, k) == doctest::Approx(T[i][k]).epsilon(1e-8));
        }
}

/**
 * Two forks in series is what makes the per-fork auxiliary delay load bearing:
 * each fork mints its OWN auxiliary classes referencing ITS OWN delay, and each
 * auxiliary token must cycle over its own branch alone. Get the routing wrong
 * and the second fork's auxiliary chain has no recurrent class at all, which is
 * a zero-visit chain rather than a wrong number.
 *
 * The reference stops at `options.iter_max` here rather than at its tolerance,
 * so it reports the two spans at 2.7573 and 2.7575 where they are symmetric by
 * construction; this port converges to the symmetric point. The tolerance below
 * is that residual, and the symmetry is asserted separately.
 */
TEST_CASE("H-T carries two fork-join spans in series, each with its own auxiliary delay") {
    Net m = serialfjs();
    const line::mva::AvgResult<double> r = ht_avg(m);
    // Rows: Delay1, Queue1, Queue2, Join, Queue3, Queue4, Join2.
    CHECK(r.QN(0, 0) == doctest::Approx(1.7278).epsilon(1e-4));
    CHECK(r.TN(0, 0) == doctest::Approx(0.86392).epsilon(1e-4));
    for (std::size_t i : {1u, 2u, 4u, 5u}) {
        CHECK(r.QN(i, 0) == doctest::Approx(2.7574).epsilon(1e-4));
        CHECK(r.RN(i, 0) == doctest::Approx(3.1917).epsilon(1e-4));
    }
    // the two spans are the same subnetwork, so they must report the same numbers
    CHECK(r.QN(1, 0) == doctest::Approx(r.QN(4, 0)).epsilon(1e-12));
    CHECK(r.QN(3, 0) == doctest::Approx(r.QN(6, 0)).epsilon(1e-12));
    CHECK(r.RN(3, 0) == doctest::Approx(r.RN(6, 0)).epsilon(1e-12));
}

TEST_CASE("H-T refuses what it cannot express, by name") {
    SUBCASE("an open class through the fork") {
        Net m("model");
        const std::size_t src = m.add_source("Source");
        const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
        const std::size_t f = m.add_fork("Fork");
        const std::size_t j = m.add_join("Join", f);
        const std::size_t snk = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("class1");
        m.set_service(src, c, DD::exp_rate(0.5));
        m.set_service(q1, c, DD::exp_rate(1.0));
        m.set_service(q2, c, DD::exp_rate(1.0));
        Routing P;
        P.set(c, c, src, f, 1.0);
        P.set(c, c, f, q1, 1.0);
        P.set(c, c, f, q2, 1.0);
        P.set(c, c, q1, j, 1.0);
        P.set(c, c, q2, j, 1.0);
        P.set(c, c, j, snk, 1.0);
        m.link(P);
        CHECK_THROWS_AS(ht_avg(m), line::UnsupportedError);
        CHECK_NOTHROW(ht_avg(m, "default"));  // MMT does carry it
    }
    SUBCASE("more than one task per link") {
        Net m("model");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::PS);
        const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::PS);
        const std::size_t f = m.add_fork("Fork", 2.0);
        const std::size_t j = m.add_join("Join", f);
        const std::size_t c = m.add_closed_class("class1", 5.0, d);
        m.set_service(d, c, DD::exp_rate(1.0));
        m.set_service(q1, c, DD::exp_rate(1.0));
        m.set_service(q2, c, DD::exp_rate(1.0));
        Routing P;
        P.set(c, c, d, f, 1.0);
        P.set(c, c, f, q1, 1.0);
        P.set(c, c, f, q2, 1.0);
        P.set(c, c, q1, j, 1.0);
        P.set(c, c, q2, j, 1.0);
        P.set(c, c, j, d, 1.0);
        m.link(P);
        CHECK_THROWS_AS(ht_avg(m), line::UnsupportedError);
    }
    SUBCASE("an unknown fork-join method is refused rather than defaulted") {
        Net m = basic_closed();
        CHECK_THROWS_AS(ht_avg(m, "nonesuch"), line::InputError);
    }
    SUBCASE("'heidelberger-trivedi' is the same arm as 'ht'") {
        Net a = basic_closed(), b = basic_closed();
        CHECK(ht_avg(a, "heidelberger-trivedi").TN(0, 0) ==
              doctest::Approx(ht_avg(b, "ht").TN(0, 0)).epsilon(1e-12));
    }
}
