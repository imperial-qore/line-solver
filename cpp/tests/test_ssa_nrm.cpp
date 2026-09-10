/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverSSA, the `nrm` method.
 *
 * HOW EACH ASSERTION IS JUSTIFIED. A simulator cannot be pinned to a golden the
 * way a deterministic algorithm can: the answer is a function of the random
 * stream, and no two LINE codebases consume draws in the same order, so an
 * exact cross-codebase comparison is meaningless here (see `ssa_types.h`).
 * Every numeric assertion below is therefore in one of three modes, and each
 * case says which:
 *
 *   ANALYTICAL + t-TEST  the model has a closed form (M/M/1, a closed
 *                        product-form network, M/M/1/K). Eight independent
 *                        seeds give eight replicate means; the two-sided
 *                        one-sample t statistic against the closed form must
 *                        fall inside the 1% critical value on 7 degrees of
 *                        freedom, 3.4995. This is the same discipline the LDES
 *                        tests use against JMT.
 *   STRUCTURAL           an invariant the engine never uses to compute
 *                        anything: closed populations sum to N, a capped
 *                        station is never over its cap, utilization lies in
 *                        [0,1], throughput balances around a cycle. These hold
 *                        on EVERY path, not just on average, so they are
 *                        asserted at full strength.
 *   REFUSAL              an unported construct must throw UnsupportedError and
 *                        the message must NAME it.
 *
 * There is no "exact agreement" mode anywhere in this file, because this port
 * shares its random stream with nothing.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"

using namespace line;
using D = lang::Distrib<double>;
using lang::SchedStrategy;

namespace {

/** Student t at 1% two-sided on 7 degrees of freedom. */
constexpr double kT7 = 3.4995;
constexpr std::size_t kReps = 8;

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

/** Sample mean and standard error of a replicate vector. */
struct Repl {
    double mean = 0.0, se = 0.0;
};
Repl summarize(const std::vector<double>& x) {
    Repl r;
    for (double v : x) r.mean += v;
    r.mean /= static_cast<double>(x.size());
    double s2 = 0.0;
    for (double v : x) s2 += (v - r.mean) * (v - r.mean);
    s2 /= static_cast<double>(x.size() - 1);
    r.se = std::sqrt(s2 / static_cast<double>(x.size()));
    return r;
}

/**
 * The t statistic of the replicates against `ref`.
 *
 * A zero standard error means every replicate produced the same number, which
 * for a simulator happens only when the quantity is deterministic along the
 * path (a saturated utilization, say); the statistic is then 0 if the value is
 * right and infinite if it is not.
 */
double tstat(const std::vector<double>& x, double ref) {
    const Repl r = summarize(x);
    if (r.se <= 0.0)
        return std::fabs(r.mean - ref) < 1e-12 ? 0.0 : std::numeric_limits<double>::infinity();
    return std::fabs(r.mean - ref) / r.se;
}

/**
 * Welch's two-sample t statistic between this port's replicates and a MATLAB
 * replicate summary (its mean and the standard error of that mean).
 *
 * This is the ONLY defensible cross-codebase comparison for a simulator whose
 * random stream differs: it asks whether the two estimators have the same
 * expectation, which is the property a correct port must have, rather than
 * whether they produced the same path, which it must not.
 */
double tstat2(const std::vector<double>& x, double ref_mean, double ref_se) {
    const Repl r = summarize(x);
    const double sd = std::sqrt(r.se * r.se + ref_se * ref_se);
    if (sd <= 0.0) return 0.0;
    return std::fabs(r.mean - ref_mean) / sd;
}

// ---------------------------------------------------------------------------
// Model builders
// ---------------------------------------------------------------------------

/** Source(lambda) -> FCFS Queue(mu) -> Sink, the textbook M/M/1. */
qn::Network<double> build_mm1(double lambda, double mu, double cap = -1.0) {
    qn::Network<double> m("mm1");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, D::exp_rate(lambda));
    m.set_service(q, c, D::exp_rate(mu));
    if (cap > 0.0) {
        m.set_capacity(q, cap);
        m.set_drop_rule(q, c, lang::DropStrategy::DROP);
    }
    qn::RoutingMatrix<double> P;
    P.set(c, c, s, q, 1.0);
    P.set(c, c, q, k, 1.0);
    m.link(P);
    return m;
}

/** Delay(Z) + PS Queue(mu), one closed class of N jobs: exact product form. */
qn::Network<double> build_cqn(double N, double z_rate, double mu, SchedStrategy qs) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", qs);
    const std::size_t c = m.add_closed_class("C1", N, d);
    m.set_service(d, c, D::exp_rate(z_rate));
    m.set_service(q, c, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// ANALYTICAL + t-TEST
// ---------------------------------------------------------------------------

TEST_CASE("ssa nrm: an open M/M/1 matches its closed form (t-test)") {
    // lambda = 0.6, mu = 1: rho = 0.6, Q = 1.5, U = 0.6, T = 0.6, R = 2.5.
    const double lambda = 0.6, mu = 1.0;
    const qsys::QsysResult<double> ref = qsys::qsys_mm1(lambda, mu);

    qn::Network<double> m0 = build_mm1(lambda, mu);
    const std::size_t iq = station_of(m0.get_struct(), "Queue");

    std::vector<double> Q, U, Tp;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 400000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build_mm1(lambda, mu);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Q.push_back(s.QN(iq, 0));
        U.push_back(s.UN(iq, 0));
        Tp.push_back(s.TN(iq, 0));
    }

    // Closed form: W = 2.5, so Q = lambda*W = 1.5 by Little; U = rho = 0.6 and
    // T = lambda = 0.6.
    CHECK(ref.W == doctest::Approx(2.5).epsilon(1e-12));
    CHECK(lambda * ref.W == doctest::Approx(1.5).epsilon(1e-12));
    CHECK(ref.rhohat == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(tstat(Q, 1.5) < kT7);
    CHECK(tstat(U, 0.6) < kT7);
    CHECK(tstat(Tp, 0.6) < kT7);

    // Little's law on the replicate means is an invariant the engine never
    // uses: R is computed as Q/T and T as a propensity integral.
    const double qm = summarize(Q).mean, tm = summarize(Tp).mean;
    CHECK(qm / tm == doctest::Approx(1.0 / (mu - lambda)).epsilon(0.05));
}

TEST_CASE("ssa nrm: a closed product-form network matches exact MVA (t-test)") {
    // Delay(Z=1) + PS Queue(mu=2), N=4. PS is symmetric, so the network is
    // product form and MVA is EXACT -- an oracle whose derivation shares
    // nothing with the sample path.
    const double N = 4, zr = 1.0, mu = 2.0;
    qn::Network<double> mref = build_cqn(N, zr, mu, SchedStrategy::PS);
    mva::MvaOptions mo;
    Matrix<double> init;
    const mva::AvgResult<double> exact = mva::solver_mva_run_analyzer(mref.get_struct(), mo, init);
    const std::size_t id = station_of(mref.get_struct(), "Delay");
    const std::size_t iq = station_of(mref.get_struct(), "Queue");

    std::vector<double> Qq, Qd, X;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 200000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build_cqn(N, zr, mu, SchedStrategy::PS);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Qq.push_back(s.QN(iq, 0));
        Qd.push_back(s.QN(id, 0));
        X.push_back(s.XN[0]);
        // STRUCTURAL: the closed population is conserved on every path.
        CHECK(s.QN(id, 0) + s.QN(iq, 0) == doctest::Approx(N).epsilon(1e-9));
    }
    CHECK(tstat(Qq, exact.QN(iq, 0)) < kT7);
    CHECK(tstat(Qd, exact.QN(id, 0)) < kT7);
    CHECK(tstat(X, exact.XN[0]) < kT7);
}

TEST_CASE("ssa nrm: a closed FCFS network matches exact MVA (t-test)") {
    // Single class, so FCFS and PS have the same product-form solution; this is
    // the case that exercises the ORDERED BUFFER (promotion on a departure,
    // join on an arrival) against an exact reference.
    const double N = 3, zr = 1.0, mu = 1.5;
    qn::Network<double> mref = build_cqn(N, zr, mu, SchedStrategy::FCFS);
    mva::MvaOptions mo;
    Matrix<double> init;
    const mva::AvgResult<double> exact = mva::solver_mva_run_analyzer(mref.get_struct(), mo, init);
    const std::size_t iq = station_of(mref.get_struct(), "Queue");

    std::vector<double> Qq, Uq;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 200000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build_cqn(N, zr, mu, SchedStrategy::FCFS);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Qq.push_back(s.QN(iq, 0));
        Uq.push_back(s.UN(iq, 0));
        // STRUCTURAL: a single-server utilization is a probability.
        CHECK(s.UN(iq, 0) >= 0.0);
        CHECK(s.UN(iq, 0) <= 1.0 + 1e-12);
    }
    CHECK(tstat(Qq, exact.QN(iq, 0)) < kT7);
    CHECK(tstat(Uq, exact.UN(iq, 0)) < kT7);
}

TEST_CASE("ssa nrm: phase-type service at a PS station matches exact MVA (t-test)") {
    // Erlang-2 service at a PS station in a closed network. PS is symmetric, so
    // the stationary law is insensitive to the service distribution beyond its
    // mean and the product form still holds EXACTLY -- which makes this the
    // strongest available check on the phase expansion: the phase machinery
    // must reproduce a number it cannot influence.
    const double N = 3, zr = 1.0;
    qn::Network<double> mref = build_cqn(N, zr, 1.5, SchedStrategy::PS);
    mva::MvaOptions mo;
    Matrix<double> init;
    const mva::AvgResult<double> exact = mva::solver_mva_run_analyzer(mref.get_struct(), mo, init);
    const std::size_t iq = station_of(mref.get_struct(), "Queue");

    auto build = []() {
        qn::Network<double> m("cqn_ph");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
        const std::size_t c = m.add_closed_class("C1", 3, d);
        m.set_service(d, c, D::exp_rate(1.0));
        m.set_service(q, c, D::erlang(3.0, 2));  // mean 2/3, i.e. rate 1.5
        qn::RoutingMatrix<double> P;
        P.set(c, c, d, q, 1.0);
        P.set(c, c, q, d, 1.0);
        m.link(P);
        return m;
    };

    // The state vector must actually have grown: two phases at the queue.
    {
        qn::Network<double> m = build();
        ssa::SsaOptions o;
        o.samples = 1;
        ssa::NrmEngine<double> eng(m.get_struct(), o);
        CHECK(eng.nstates() == 3);  // Delay 1 phase + Queue 2 phases, one class
    }

    std::vector<double> Qq;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 200000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build();
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Qq.push_back(s.QN(iq, 0));
        CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(N).epsilon(1e-9));
    }
    CHECK(tstat(Qq, exact.QN(iq, 0)) < kT7);
}

TEST_CASE("ssa nrm: a multiserver PS station matches exact MVA (t-test)") {
    const double N = 5, zr = 1.0, mu = 1.0;
    qn::Network<double> mref = build_cqn(N, zr, mu, SchedStrategy::PS);
    mref.set_number_of_servers(2, 2);
    mva::MvaOptions mo;
    Matrix<double> init;
    const mva::AvgResult<double> exact = mva::solver_mva_run_analyzer(mref.get_struct(), mo, init);
    const std::size_t iq = station_of(mref.get_struct(), "Queue");

    std::vector<double> Qq;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 200000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build_cqn(N, zr, mu, SchedStrategy::PS);
        m.set_number_of_servers(2, 2);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Qq.push_back(s.QN(iq, 0));
    }
    CHECK(tstat(Qq, exact.QN(iq, 0)) < kT7);
}

TEST_CASE("ssa nrm: a finite-capacity M/M/1/K loses arrivals at the cap") {
    // lambda = mu = 1, K = 3. The truncated birth-death chain is uniform on
    // 0..K, so Q = K/2 = 1.5 and the carried throughput is lambda*(1 - p_K) =
    // 1 - 1/(K+1) = 0.75.
    const double lambda = 1.0, mu = 1.0, K = 3.0;
    qn::Network<double> m0 = build_mm1(lambda, mu, K);
    const std::size_t iq = station_of(m0.get_struct(), "Queue");

    std::vector<double> Qq, Tq;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 400000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build_mm1(lambda, mu, K);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Qq.push_back(s.QN(iq, 0));
        Tq.push_back(s.TN(iq, 0));
        // STRUCTURAL: the time-average population cannot exceed the cap. This
        // is the invariant the reference's NRM violated before the loss branch
        // existed, when a capped queue overflowed well past sn.cap.
        CHECK(s.QN(iq, 0) <= K + 1e-9);
    }
    CHECK(tstat(Qq, 1.5) < kT7);
    CHECK(tstat(Tq, 0.75) < kT7);
}

TEST_CASE("ssa nrm: a closed job that finds no room BLOCKS, it is not admitted") {
    // BUG-81. `capacity_loss` above declines to DROP a closed job -- a closed
    // network's population is an invariant -- but until 2026-08-19 nothing then
    // stopped the reaction, so the firing went into the full station anyway and
    // every codebase's NRM returned the UNCONSTRAINED product-form answer:
    // QLen [1.96 2.07 1.97] on this model, against the exact
    // [3.60897726 0.97107998 1.41994276]. A mean of 2.07 at a station capped at
    // 2 is not attainable, which is the assertion that catches it; population
    // conservation does NOT, because the broken engine conserved it too.
    //
    // The exact figures are the stationary solution of the constrained chain,
    // hard-coded rather than solved here so the assertion cannot drift with a
    // reference solver; MATLAB and Python SolverCTMC agree on them to the digits
    // shown, and both serial SSA engines reproduce them.
    auto build = [](double cap) {
        qn::Network<double> m("tandem");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
        const std::size_t q3 = m.add_queue("Q3", SchedStrategy::FCFS);
        const std::size_t c = m.add_closed_class("C1", 6, q1);
        m.set_service(q1, c, D::exp_rate(1.0));
        m.set_service(q2, c, D::exp_rate(1.0));
        m.set_service(q3, c, D::exp_rate(1.0));
        m.set_class_capacity(q2, c, cap);
        qn::RoutingMatrix<double> P;
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, q3, 1.0);
        P.set(c, c, q3, q1, 1.0);
        m.link(P);
        return m;
    };
    const double exactQ[3] = {3.60897726, 0.97107998, 1.41994276};
    const double exactX = 0.65220666;

    std::vector<double> Q1r, Q2r, Q3r, Xr;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 200000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build(2.0);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        // STRUCTURAL, and the assertion the defect fails: a station that holds
        // two jobs cannot average more than two.
        CHECK(s.QN(1, 0) <= 2.0 + 1e-9);
        // The six jobs cannot leave a closed network, blocked or not.
        CHECK(s.QN(0, 0) + s.QN(1, 0) + s.QN(2, 0) == doctest::Approx(6.0).epsilon(1e-9));
        Q1r.push_back(s.QN(0, 0));
        Q2r.push_back(s.QN(1, 0));
        Q3r.push_back(s.QN(2, 0));
        Xr.push_back(s.TN(0, 0));
    }
    CHECK(tstat(Q1r, exactQ[0]) < kT7);
    CHECK(tstat(Q2r, exactQ[1]) < kT7);
    CHECK(tstat(Q3r, exactQ[2]) < kT7);
    // The departure-rate integral must net out the firings that were blocked:
    // uncorrected it reads the unconstrained 0.744 rather than 0.652.
    CHECK(tstat(Xr, exactX) < kT7);
}

TEST_CASE("ssa nrm: a two-class closed network conserves both populations") {
    // Two classes routed in OPPOSITE directions over three stations, one of
    // them a two-server FCFS queue. No closed form is claimed here; what is
    // asserted is the pair of conservation laws and the flow balance around
    // each cycle, none of which the engine computes from the other.
    auto build = []() {
        qn::Network<double> m("cqn2");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
        const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
        const std::size_t c1 = m.add_closed_class("C1", 5, d);
        const std::size_t c2 = m.add_closed_class("C2", 3, d);
        m.set_service(d, c1, D::exp_rate(2.0));
        m.set_service(d, c2, D::exp_rate(1.0));
        m.set_service(q1, c1, D::exp_rate(3.0));
        m.set_service(q1, c2, D::exp_rate(1.5));
        m.set_service(q2, c1, D::exp_rate(1.2));
        m.set_service(q2, c2, D::exp_rate(0.8));
        m.set_number_of_servers(q2, 2);
        qn::RoutingMatrix<double> P;
        P.set(c1, c1, d, q1, 1.0);
        P.set(c1, c1, q1, q2, 1.0);
        P.set(c1, c1, q2, d, 1.0);
        P.set(c2, c2, d, q2, 1.0);
        P.set(c2, c2, q2, q1, 1.0);
        P.set(c2, c2, q1, d, 1.0);
        m.link(P);
        return m;
    };
    ssa::SsaOptions o;
    o.samples = 300000;
    o.seed = 23000;
    qn::Network<double> m = build();
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
    double n1 = 0.0, n2 = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        n1 += s.QN(i, 0);
        n2 += s.QN(i, 1);
    }
    CHECK(n1 == doctest::Approx(5.0).epsilon(1e-9));
    CHECK(n2 == doctest::Approx(3.0).epsilon(1e-9));
    // Flow balance on a tandem cycle: every station on class r's cycle carries
    // the same throughput. 2% is the Monte Carlo spread at this sample count.
    for (std::size_t i = 1; i < 3; ++i) {
        CHECK(s.TN(i, 0) == doctest::Approx(s.TN(0, 0)).epsilon(0.02));
        CHECK(s.TN(i, 1) == doctest::Approx(s.TN(0, 1)).epsilon(0.02));
    }
    // Little's law per station, which the engine derives but never checks.
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t r = 0; r < 2; ++r)
            if (s.TN(i, r) > 0.0)
                CHECK(s.RN(i, r) == doctest::Approx(s.QN(i, r) / s.TN(i, r)).epsilon(1e-9));
}

TEST_CASE("ssa nrm: the buffered disciplines all run and conserve population") {
    // FCFS, LCFS, SIRO, HOL, SEPT, LEPT and LCFSPR share the rate law and
    // differ only in which waiting job is promoted, so on a SINGLE-CLASS closed
    // network they must all reproduce the same product-form marginal. That is
    // an oracle none of the promotion rules can influence.
    const double N = 3, zr = 1.0, mu = 1.5;
    qn::Network<double> mref = build_cqn(N, zr, mu, SchedStrategy::PS);
    mva::MvaOptions mo;
    Matrix<double> init;
    const mva::AvgResult<double> exact = mva::solver_mva_run_analyzer(mref.get_struct(), mo, init);
    const std::size_t iq = station_of(mref.get_struct(), "Queue");

    const SchedStrategy pol[] = {SchedStrategy::FCFS, SchedStrategy::LCFS, SchedStrategy::SIRO,
                                 SchedStrategy::HOL,  SchedStrategy::SEPT, SchedStrategy::LEPT,
                                 SchedStrategy::LCFSPR};
    for (SchedStrategy p : pol) {
        std::vector<double> Qq;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 200000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_cqn(N, zr, mu, p);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Qq.push_back(s.QN(iq, 0));
            CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(N).epsilon(1e-9));
        }
        INFO("policy ", lang::sched_to_text(p));
        CHECK(tstat(Qq, exact.QN(iq, 0)) < kT7);
    }
}

TEST_CASE("ssa nrm: DPS with equal weights reduces to PS (t-test)") {
    // Two classes at one DPS station with EQUAL weights is exactly PS, which is
    // product form, so MVA is the oracle. Unequal weights have no closed form,
    // so the DPS-specific branch is checked by the degenerate case plus the
    // monotone response of the class shares below.
    auto build = [](double w1, double w2) {
        qn::Network<double> m("dps");
        const std::size_t d = m.add_delay("Delay");
        const std::size_t q = m.add_queue("Q", SchedStrategy::DPS);
        const std::size_t c1 = m.add_closed_class("C1", 2, d);
        const std::size_t c2 = m.add_closed_class("C2", 2, d);
        m.set_service(d, c1, D::exp_rate(1.0));
        m.set_service(d, c2, D::exp_rate(1.0));
        m.set_service(q, c1, D::exp_rate(2.0));
        m.set_service(q, c2, D::exp_rate(2.0));
        m.set_sched_param(q, c1, w1);
        m.set_sched_param(q, c2, w2);
        qn::RoutingMatrix<double> P;
        P.set(c1, c1, d, q, 1.0);
        P.set(c1, c1, q, d, 1.0);
        P.set(c2, c2, d, q, 1.0);
        P.set(c2, c2, q, d, 1.0);
        m.link(P);
        return m;
    };
    ssa::SsaOptions o;
    o.samples = 300000;
    o.seed = 23000;
    // Equal weights: the two symmetric classes must get the same queue length.
    {
        qn::Network<double> m = build(1.0, 1.0);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        CHECK(s.QN(1, 0) == doctest::Approx(s.QN(1, 1)).epsilon(0.03));
    }
    // Four times the weight must leave class 1 with the shorter queue; the
    // direction is an invariant of the DPS share and needs no reference value.
    {
        qn::Network<double> m = build(4.0, 1.0);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        CHECK(s.QN(1, 0) < s.QN(1, 1));
        CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(s.QN(0, 1) + s.QN(1, 1) == doctest::Approx(2.0).epsilon(1e-9));
    }
}

TEST_CASE("ssa nrm: GPS shares by active classes, not by jobs") {
    // GPS gives class r the share w_r / sum over ACTIVE classes, so a class
    // with one job gets the same capacity as one with many. The check is the
    // structural one that distinguishes GPS from DPS: with equal weights and
    // very unequal populations the per-class utilizations must be nearly equal
    // under GPS, where DPS would split them in proportion to the populations.
    qn::Network<double> m("gps");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::GPS);
    const std::size_t c1 = m.add_closed_class("C1", 8, d);
    const std::size_t c2 = m.add_closed_class("C2", 1, d);
    m.set_service(d, c1, D::exp_rate(10.0));
    m.set_service(d, c2, D::exp_rate(10.0));
    m.set_service(q, c1, D::exp_rate(1.0));
    m.set_service(q, c2, D::exp_rate(1.0));
    m.set_sched_param(q, c1, 1.0);
    m.set_sched_param(q, c2, 1.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    o.samples = 300000;
    o.seed = 23000;
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
    // Both classes are essentially always present, so both hold half the
    // server; the throughputs, not the queue lengths, are what GPS equalizes.
    CHECK(s.TN(1, 0) == doctest::Approx(s.TN(1, 1)).epsilon(0.10));
    CHECK(s.UN(1, 0) + s.UN(1, 1) <= 1.0 + 1e-9);
}

TEST_CASE("ssa nrm: an INF station is its own utilization") {
    // At an infinite server every job present is in service, so U == Q exactly,
    // on every path and at every sample count.
    qn::Network<double> m = build_cqn(4, 1.0, 2.0, SchedStrategy::PS);
    ssa::SsaOptions o;
    o.samples = 50000;
    o.seed = 23000;
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
    CHECK(s.UN(0, 0) == doctest::Approx(s.QN(0, 0)).epsilon(1e-12));
}

TEST_CASE("ssa nrm: the reaction network has the shape the layout predicts") {
    // A structural pin on the builder: Delay + Queue + one class with Erlang-3
    // service at the queue gives 1 + 3 state slots, 4 departure reactions and 2
    // phase transitions (the Erlang chain 1->2->3).
    qn::Network<double> m("shape");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(3.0, 3));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    o.samples = 1;
    ssa::NrmEngine<double> eng(m.get_struct(), o);
    CHECK(eng.nstates() == 4);
    CHECK(eng.nreactions() == 6);
}

TEST_CASE("ssa nrm: a phase change at a full capped station is not a loss") {
    // M/Er2/1/K under PS, lambda = 1, mean service 1, K = 3. PS is insensitive,
    // so the stationary law is that of M/M/1/K with the same mean: uniform on
    // 0..3, hence Q = 1.5 and carried throughput lambda*(1 - p_K) = 0.75 --
    // EXACTLY the exponential case, which is what makes this a clean test.
    //
    // This is the model on which MATLAB, the JAR and native Python all destroy
    // jobs. Their finite-capacity loss test is applied to the destination of
    // every single-destination firing, and a PHASE CHANGE is one, so once the
    // queue holds K jobs the Erlang phase transition of the job in service is
    // declared a loss and the job disappears. This port excludes phase changes
    // from the test, so it must reproduce the insensitive answer.
    auto build = []() {
        qn::Network<double> m("mphk");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(s, c, D::exp_rate(1.0));
        m.set_service(q, c, D::erlang(2.0, 2));  // mean 1, order 2
        m.set_capacity(q, 3.0);
        m.set_drop_rule(q, c, lang::DropStrategy::DROP);
        qn::RoutingMatrix<double> P;
        P.set(c, c, s, q, 1.0);
        P.set(c, c, q, k, 1.0);
        m.link(P);
        return m;
    };
    std::vector<double> Q, Tp;
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaOptions o;
        o.samples = 400000;
        o.seed = 23000 + 7919 * rep;
        qn::Network<double> m = build();
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        Q.push_back(s.QN(1, 0));
        Tp.push_back(s.TN(1, 0));
        // STRUCTURAL: the cap is never exceeded, in either reading of the loss
        // rule, so this cannot be what the t-tests below are detecting.
        CHECK(s.QN(1, 0) <= 3.0 + 1e-9);
    }
    CHECK(tstat(Q, 1.5) < kT7);
    CHECK(tstat(Tp, 0.75) < kT7);
}

TEST_CASE("ssa nrm: HyperExp separates multi-phase from phase-CHANGE firings") {
    // The discriminator for the defect above, and the case that would catch a
    // future "fix" special-casing all multi-phase service. HyperExp2 has TWO
    // phases but NO phase-change firings: the phase is drawn at the start of
    // service and the job completes from it, so its D0 has no off-diagonal
    // mass. Erlang2 and Erlang3 do have such firings, and their D0 off-diagonal
    // mass grows with the order -- which is why the defect was monotone in the
    // phase count where it existed.
    //
    // All four service laws have mean 1 and sit at the same PS station with the
    // same cap, so PS insensitivity says all four must give the SAME answer,
    // the M/M/1/3 one. A capacity gate that mistook a phase change for an
    // arrival would leave HyperExp2 clean and drag the Erlangs down.
    auto build = [](int variant) {
        qn::Network<double> m("mphk");
        const std::size_t s = m.add_source("Source");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
        const std::size_t k = m.add_sink("Sink");
        const std::size_t c = m.add_open_class("C1");
        m.set_arrival(s, c, D::exp_rate(1.0));
        // Every variant has MEAN 1. The HyperExp is the balanced-means form
        // with SCV 4: p/l1 = (1-p)/l2 = 1/2 forces l1+l2 = 2 and l1*l2 = 2/5,
        // so l = 1 +- sqrt(3/5) and p = l1/2. Built with hyperexp() and NOT
        // erlang_fit(1,4), which follows MATLAB in taking r = ceil(1/SCV) = 1
        // and so collapses to an exponential -- no phases at all, and the
        // discriminator would test nothing.
        if (variant == 0) m.set_service(q, c, D::exp_rate(1.0));
        else if (variant == 1)
            m.set_service(q, c, D::hyperexp(0.8872983346207417, 1.7745966692414834,
                                            0.22540333075851657));
        else if (variant == 2) m.set_service(q, c, D::erlang(2.0, 2));
        else m.set_service(q, c, D::erlang(3.0, 3));
        m.set_capacity(q, 3.0);
        m.set_drop_rule(q, c, lang::DropStrategy::DROP);
        qn::RoutingMatrix<double> P;
        P.set(c, c, s, q, 1.0);
        P.set(c, c, q, k, 1.0);
        m.link(P);
        return m;
    };
    // The premise, asserted rather than assumed: HyperExp2 must have TWO state
    // slots at the queue like Erlang2, but the SAME reaction count as Exp,
    // because it contributes no phase-change reaction.
    {
        ssa::SsaOptions o;
        o.samples = 1;
        qn::Network<double> m0 = build(0), m1 = build(1), m2 = build(2);
        ssa::NrmEngine<double> e0(m0.get_struct(), o), e1(m1.get_struct(), o),
            e2(m2.get_struct(), o);
        CHECK(e1.nstates() == e2.nstates());        // both carry two phases
        CHECK(e1.nstates() == e0.nstates() + 1);
        CHECK(e1.nreactions() == e0.nreactions() + 1);   // one extra departure only
        CHECK(e2.nreactions() == e0.nreactions() + 2);   // departure AND phase change
    }
    const char* nm[] = {"Exp", "HyperExp2", "Erlang2", "Erlang3"};
    for (int v = 0; v < 4; ++v) {
        std::vector<double> Q, Tp;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 400000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build(v);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Q.push_back(s.QN(1, 0));
            Tp.push_back(s.TN(1, 0));
        }
        INFO("service law ", nm[v]);
        CHECK(tstat(Q, 1.5) < kT7);
        CHECK(tstat(Tp, 0.75) < kT7);
        // STRUCTURAL: the gate must still be LIVE. A lazy fix that simply
        // stopped dropping would push the carried throughput toward the offered
        // 1.0 and the queue past the cap-limited 1.5, so the loss rate is
        // asserted from below as well as the metrics from above.
        const double loss = 1.0 - summarize(Tp).mean;
        CHECK(loss > 0.20);
        CHECK(loss < 0.30);
    }
}

TEST_CASE("ssa nrm: a Source reports a zero row, never a negative one") {
    // The Source's state slot is a FICTITIOUS TOKEN seeded at 1, decremented by
    // every arrival and replenished only when a job reaches the Sink, so its
    // time average is 1 - E[jobs in system]. Reported as a queue length it goes
    // NEGATIVE on any model holding more than one job -- observed as
    // QLen = -0.0294588 against a Queue at 1.02946, the two summing to exactly
    // 1. A Source holds no jobs, so its QLen and Util are zero BY DEFINITION
    // and only its throughput, the arrival rate, is a quantity.
    //
    // The assertions are exact equalities, not tolerances: this is a definition
    // and not an estimate, so any nonzero value is a defect however small.
    const double lambda = 0.5, mu = 1.0;
    qn::Network<double> m = build_mm1(lambda, mu);
    const std::size_t isrc = station_of(m.get_struct(), "Source");
    const std::size_t iq = station_of(m.get_struct(), "Queue");
    ssa::SsaOptions o;
    o.samples = 200000;
    o.seed = 23000;
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);

    CHECK(s.QN(isrc, 0) == 0.0);
    CHECK(s.UN(isrc, 0) == 0.0);
    CHECK(s.RN(isrc, 0) == 0.0);
    // The Source's throughput IS meaningful and must survive: it is the
    // arrival rate. Zeroing the whole row would pass a naive "not negative"
    // test while destroying the one quantity a Source has.
    CHECK(s.TN(isrc, 0) == doctest::Approx(lambda).epsilon(0.05));

    // And the fix must not have touched the Queue, which is where the real
    // metrics live. rho = 0.5 gives QLen 1 and Util 0.5 exactly.
    CHECK(s.QN(iq, 0) == doctest::Approx(1.0).epsilon(0.10));
    CHECK(s.UN(iq, 0) == doctest::Approx(0.5).epsilon(0.10));
    CHECK(s.QN(iq, 0) > 0.0);
}

TEST_CASE("ssa nrm: no station ever reports a negative metric") {
    // The generalisation of the case above, over every discipline this port
    // supports and both open and closed models. A queue length, a utilization
    // and a throughput are all nonnegative quantities; the Source defect was
    // the only way this port produced one, and this is what would catch a
    // second route to it.
    const SchedStrategy pol[] = {SchedStrategy::FCFS, SchedStrategy::PS,
                                 SchedStrategy::LCFS, SchedStrategy::SIRO,
                                 SchedStrategy::HOL,  SchedStrategy::SEPT,
                                 SchedStrategy::LEPT, SchedStrategy::LCFSPR};
    for (SchedStrategy p : pol) {
        INFO("policy ", lang::sched_to_text(p));
        // open
        {
            qn::Network<double> m("mm1p");
            const std::size_t src = m.add_source("Source");
            const std::size_t q = m.add_queue("Queue", p);
            const std::size_t k = m.add_sink("Sink");
            const std::size_t c = m.add_open_class("C1");
            m.set_arrival(src, c, D::exp_rate(0.7));
            m.set_service(q, c, D::exp_rate(1.0));
            qn::RoutingMatrix<double> P;
            P.set(c, c, src, q, 1.0);
            P.set(c, c, q, k, 1.0);
            m.link(P);
            ssa::SsaOptions o;
            o.samples = 100000;
            o.seed = 23000;
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            for (std::size_t i = 0; i < m.get_struct().nstations; ++i) {
                CHECK(s.QN(i, 0) >= 0.0);
                CHECK(s.UN(i, 0) >= 0.0);
                CHECK(s.RN(i, 0) >= 0.0);
                CHECK(s.TN(i, 0) >= 0.0);
            }
        }
        // closed
        {
            qn::Network<double> m = build_cqn(3, 1.0, 1.5, p);
            ssa::SsaOptions o;
            o.samples = 100000;
            o.seed = 23000;
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            for (std::size_t i = 0; i < m.get_struct().nstations; ++i) {
                CHECK(s.QN(i, 0) >= 0.0);
                CHECK(s.UN(i, 0) >= 0.0);
                CHECK(s.RN(i, 0) >= 0.0);
                CHECK(s.TN(i, 0) >= 0.0);
            }
        }
    }
}

TEST_CASE("ssa nrm: open station sums equal the model total, not one") {
    // THE DISCRIMINATING TEST between two hypotheses for the "station QLen sums
    // to exactly 1" observation on open models:
    //   (A) the engine NORMALIZES station queue lengths to a unit job pool and
    //       gives the remainder to the Source, which would bias every station;
    //   (B) the Source slot holds `1 - N(t)` by construction, so the sum is
    //       1 - N + N == 1 as an ARITHMETIC IDENTITY and no station is biased.
    // They are separated by the SUM AGAINST THE MODEL'S ACTUAL TOTAL once the
    // Source is out of the table. Under (A) the sum stays pinned at 1 whatever
    // the model; under (B) it equals the true total.
    //
    // Both directions are covered, because a test keyed on the negative row
    // catches only half the parameter space: at rho = 0.3 the same defect would
    // produce an INFLATION to 1 from a true 0.43, with every printed number
    // looking plausible and no negative anywhere.
    //
    // References are the M/M/1 closed form Q = rho/(1-rho), exact, not another
    // implementation.
    struct Case { double lambda, mu, exact; const char* name; };
    const Case cases[] = {
        {0.3, 1.0, 0.3 / 0.7, "rho=0.3, INFLATION case: true 0.4286, a unit pool would say 1"},
        {0.5, 1.0, 1.0, "rho=0.5, the boundary where a unit pool is invisible"},
        {0.8, 1.0, 4.0, "rho=0.8, DEFLATION case: true 4, a unit pool would say 1"},
    };
    for (const Case& cs : cases) {
        INFO(cs.name);
        std::vector<double> sums;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 2000000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_mm1(cs.lambda, cs.mu);
            const qn::NetworkStruct<double>& sn = m.get_struct();
            const ssa::SsaSolution s = ssa::solver_ssa(sn, o);
            double tot = 0.0;
            for (std::size_t i = 0; i < sn.nstations; ++i) tot += s.QN(i, 0);
            sums.push_back(tot);
        }
        CHECK(tstat(sums, cs.exact) < kT7);
        // And the sum must be nowhere near 1 when the model's total is not 1,
        // which is the assertion a unit-pool normalization would fail outright.
        if (std::fabs(cs.exact - 1.0) > 0.3) {
            const double m = summarize(sums).mean;
            CHECK(std::fabs(m - 1.0) > 0.3);
        }
    }
}

TEST_CASE("ssa nrm: a closed N=1 model sums to one LEGITIMATELY") {
    // The control for the test above. A closed model of population 1 has a
    // station sum of exactly 1 because that IS its population, so "the sum is
    // 1" is not on its own evidence of anything. A regression keyed on the
    // naive predicate would fire here and be disabled by the next reader.
    qn::Network<double> m = build_cqn(1, 1.0, 2.0, SchedStrategy::PS);
    ssa::SsaOptions o;
    o.samples = 200000;
    o.seed = 23000;
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ssa::SsaSolution s = ssa::solver_ssa(sn, o);
    double tot = 0.0;
    for (std::size_t i = 0; i < sn.nstations; ++i) tot += s.QN(i, 0);
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// CROSS-CODEBASE, STATISTICAL
// ---------------------------------------------------------------------------

/**
 * MATLAB `SolverSSA(model, 'method','nrm')` on the same six models, five seeds
 * each at 200000 samples, reported as (mean, standard error of the mean).
 * Produced by the batch in this port's scratchpad and transcribed here.
 *
 * These are NOT goldens. Two simulators with different random streams cannot
 * agree on a number; what they must agree on is the EXPECTATION, so each is
 * compared by Welch's two-sample t against this port's own replicates. The
 * threshold is the same 3.4995 the one-sample cases use, which is conservative
 * for the Welch degrees of freedom here (about ten).
 */
TEST_CASE("ssa nrm: agrees statistically with MATLAB SolverSSA(method='nrm')") {
    SUBCASE("open M/M/1") {
        std::vector<double> Q, U, Tp;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 400000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_mm1(0.6, 1.0);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Q.push_back(s.QN(1, 0));
            U.push_back(s.UN(1, 0));
            Tp.push_back(s.TN(1, 0));
        }
        CHECK(tstat2(Q, 1.515708, 0.009971) < kT7);
        CHECK(tstat2(U, 0.601748, 0.001164) < kT7);
        CHECK(tstat2(Tp, 0.601748, 0.001164) < kT7);
    }
    SUBCASE("closed PS, N=4") {
        std::vector<double> Qd, Qq, X;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 200000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_cqn(4, 1.0, 2.0, SchedStrategy::PS);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Qd.push_back(s.QN(0, 0));
            Qq.push_back(s.QN(1, 0));
            X.push_back(s.XN[0]);
        }
        CHECK(tstat2(Qd, 1.809369, 0.003077) < kT7);
        CHECK(tstat2(Qq, 2.190631, 0.003077) < kT7);
        CHECK(tstat2(X, 1.809276, 0.000692) < kT7);
    }
    SUBCASE("closed FCFS, N=3") {
        std::vector<double> Qq, Uq;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 200000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_cqn(3, 1.0, 1.5, SchedStrategy::FCFS);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Qq.push_back(s.QN(1, 0));
            Uq.push_back(s.UN(1, 0));
        }
        CHECK(tstat2(Qq, 1.695433, 0.001629) < kT7);
        CHECK(tstat2(Uq, 0.863803, 0.000198) < kT7);
    }
    SUBCASE("M/M/1/K, K=3") {
        std::vector<double> Qq, Tq;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 400000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m = build_mm1(1.0, 1.0, 3.0);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Qq.push_back(s.QN(1, 0));
            Tq.push_back(s.TN(1, 0));
        }
        CHECK(tstat2(Qq, 1.498711, 0.003876) < kT7);
        CHECK(tstat2(Tq, 0.749284, 0.001395) < kT7);
    }
    SUBCASE("closed PS with Erlang-2 service, N=3") {
        std::vector<double> Qq, X;
        for (std::size_t rep = 0; rep < kReps; ++rep) {
            ssa::SsaOptions o;
            o.samples = 200000;
            o.seed = 23000 + 7919 * rep;
            qn::Network<double> m("cqn_ph");
            const std::size_t d = m.add_delay("Delay");
            const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
            const std::size_t c = m.add_closed_class("C1", 3, d);
            m.set_service(d, c, D::exp_rate(1.0));
            m.set_service(q, c, D::erlang(3.0, 2));
            qn::RoutingMatrix<double> P;
            P.set(c, c, d, q, 1.0);
            P.set(c, c, q, d, 1.0);
            m.link(P);
            const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
            Qq.push_back(s.QN(1, 0));
            X.push_back(s.XN[0]);
        }
        CHECK(tstat2(Qq, 1.703814, 0.001088) < kT7);
        CHECK(tstat2(X, 1.297996, 0.001739) < kT7);
    }
}

// ---------------------------------------------------------------------------
// REFUSALS, each asserted BY NAME
// ---------------------------------------------------------------------------

// 2026-07-28: the refusal this case used to assert was wrong. SolverSSA.m:55
// lists 'serial', 'para' and 'parallel' as valid methods and
// solver_ssa_analyzer.m:135-176 dispatches all three, so with both engines now
// ported the dispatcher answers them instead of refusing.
TEST_CASE("ssa nrm: the serial and parallel methods are dispatched, not refused") {
    qn::Network<double> m = build_cqn(2, 1.0, 2.0, SchedStrategy::PS);
    // 'serial' and its 'ssa' alias reach the event-driven engine
    // (solver_ssa_analyzer.m:128-131, 136-140).
    for (const std::string& meth : {std::string("serial"), std::string("ssa")}) {
        ssa::SsaOptions o;
        o.method = meth;
        o.samples = 2000;
        o.seed = 4242;
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        INFO("method ", meth);
        CHECK(s.method == "serial");
        CHECK(s.QN(1, 0) > 0.0);
    }
    // 'para' / 'parallel' prefer the NRM on an eligible model, exactly as
    // solver_ssa_analyzer.m:143-157 does.
    for (const std::string& meth : {std::string("parallel"), std::string("para")}) {
        ssa::SsaOptions o;
        o.method = meth;
        o.samples = 2000;
        o.seed = 4242;
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        INFO("method ", meth);
        CHECK(s.method == "nrm");
        CHECK(s.QN(1, 0) > 0.0);
    }
}

// The other side of the same dispatch: on a model the NRM cannot run,
// 'parallel' does NOT stop at the NRM gate. LCFSPR with phase-type service
// fails phaseNrmOK (solver_ssa_analyzer.m:402-429), so the dispatcher hands it
// to the replicated serial engine, which RUNS it: the serial engine's LCFSPR
// arm is ported, so what says the model got past the gate is that an answer
// comes back and no message mentions the gate. 'default' is deliberately not
// asserted here: it keeps the NRM's message, see ssa_dispatch.h.
TEST_CASE("ssa nrm: an NRM-ineligible model reaches the replicated engine") {
    qn::Network<double> m("cqn_lcfspr_ph");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::LCFSPR);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(3.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);

    ssa::SsaOptions o;
    o.samples = 200;
    o.seed = 909;

    o.method = "nrm";  // spelled out, so the NRM gate refuses by name
    std::string gate;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        gate = e.what();
    }
    CHECK(gate.find("SolverSSA(method='nrm')") != std::string::npos);
    CHECK(gate.find("lcfspr") != std::string::npos);

    // Same model, method 'parallel': the replicated engine runs it to a number,
    // so the gate is behind it. No message at all is a STRONGER statement than
    // the engine naming its own gap, which is what this asserted while the
    // serial engine's LCFSPR arm was missing.
    o.method = "parallel";
    std::string eng;
    double qn = -1.0;
    try {
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), o);
        qn = s.QN(1, 0);
    } catch (const UnsupportedError& e) {
        eng = e.what();
    }
    CHECK(eng.empty());
    CHECK(qn > 0.0);
}

TEST_CASE("ssa nrm: an unknown method refuses by name") {
    qn::Network<double> m = build_cqn(2, 1.0, 2.0, SchedStrategy::PS);
    ssa::SsaOptions o;
    o.method = "gillespie";
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("gillespie") != std::string::npos);
    CHECK(msg.find("'nrm'") != std::string::npos);
}

TEST_CASE("ssa nrm: a non-double arithmetic refuses by name") {
    qn::Network<Rational> m("pf");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C", 2, d);
    m.set_service(d, c, lang::Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(1)));
    m.set_service(q, c, lang::Distrib<Rational>::exp_rate(num_traits<Rational>::from_int(2)));
    qn::RoutingMatrix<Rational> P;
    P.set(c, c, d, q, num_traits<Rational>::from_int(1));
    P.set(c, c, q, d, num_traits<Rational>::from_int(1));
    m.link(P);
    ssa::SsaOptions o;
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("exponential clocks") != std::string::npos);
    CHECK(msg.find("--arith double") != std::string::npos);
}

TEST_CASE("ssa nrm: a POLLING station refuses by name") {
    qn::Network<double> m("poll");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("PollQ", SchedStrategy::POLLING);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_polling_type(q, lang::PollingType::EXHAUSTIVE);
    m.set_switchover(q, c, D::exp_rate(5.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("PollQ") != std::string::npos);
    CHECK(msg.find("POLLING") != std::string::npos);
}

TEST_CASE("ssa nrm: a pass-and-swap station refuses by name") {
    qn::Network<double> m("pas");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("OiQ", SchedStrategy::PAS);
    const std::size_t c = m.add_closed_class("C1", 2, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service_rate_function(q, [](const std::vector<std::size_t>& cs) {
        return static_cast<double>(cs.size());
    });
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("OiQ") != std::string::npos);
    CHECK(msg.find("pass-and-swap") != std::string::npos);
}

TEST_CASE("ssa nrm: a Cache node refuses by name") {
    qn::Network<double> m("cache");
    const std::size_t s = m.add_source("Source");
    qn::CacheParam<double> cp;
    cp.nitems = 4;
    cp.itemcap.push_back(2);
    cp.replacestrat = lang::ReplacementStrategy::LRU;
    cp.pread.assign(2, std::vector<double>());
    cp.pread[0] = std::vector<double>(4, 0.25);
    cp.hitclass.assign(2, 0);
    cp.missclass.assign(2, 0);
    cp.hitclass[0] = 2;
    cp.missclass[0] = 2;
    const std::size_t ca = m.add_cache("Cache1", cp);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("Read");
    const std::size_t c2 = m.add_open_class("Out");
    m.set_arrival(s, c1, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, s, ca, 1.0);
    P.set(c2, c2, ca, k, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("Cache1") != std::string::npos);
    CHECK(msg.find("Cache") != std::string::npos);
}

TEST_CASE("ssa nrm: class-dependent scaling refuses by name") {
    // The refusal must name the construct AND the remedy. The reason moved on
    // 2026-07-30: the NRM used to be refused for want of a declared peak, but it
    // now carries one, so what disqualifies it is that it builds its reaction
    // rates without ever evaluating the handle -- the sample path would run at
    // the unscaled rates while the utilization still looked consistent.
    qn::Network<double> m = build_cqn(3, 1.0, 2.0, SchedStrategy::PS);
    m.set_class_dependence(2, [](const std::vector<double>& n) {
        return std::vector<double>(n.size(), 1.0);
    });
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("class- or joint-dependent scaling") != std::string::npos);
    CHECK(msg.find("without evaluating the dependence handle") != std::string::npos);
    CHECK(msg.find("method='serial'") != std::string::npos);
}

TEST_CASE("ssa nrm: joint-dependent scaling refuses on the same gate") {
    // The gate was widened to eta_i(n) in the same change that reworded it, and
    // the NRM evaluates that handle no more than it evaluates beta_r(n).
    qn::Network<double> m = build_cqn(3, 1.0, 2.0, SchedStrategy::PS);
    m.set_joint_dependence(2, [](const std::vector<double>&) {
        return std::vector<double>(1, 1.0);
    }, std::vector<double>(1, 1.0));
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("class- or joint-dependent scaling") != std::string::npos);
    CHECK(msg.find("method='serial'") != std::string::npos);
}

TEST_CASE("ssa nrm: phase-type service at LCFSPR refuses by name") {
    // The gate mirrors phaseNrmOK: preempt-resume would need the preempted
    // job's phase remembered in the buffer, which the buffer does not record.
    qn::Network<double> m("lcfspr_ph");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::LCFSPR);
    const std::size_t c = m.add_closed_class("C1", 3, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    // NAMED, not `default`: since the reference's fallback was restored
    // (2026-07-31) `default` answers this model on the serial engine, and the
    // refusal below belongs to the caller who asked for the NRM by name.
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("non-exponential") != std::string::npos);
    CHECK(msg.find("lcfspr") != std::string::npos);
}

TEST_CASE("ssa nrm: a finite capacity region refuses under nrm and is answered under default") {
    // Before this gate, a region reaching the NRM was silently simulated as absent.
    qn::Network<double> mm = build_cqn(3, 1.0, 2.0, SchedStrategy::FCFS);
    mm.add_region(std::vector<std::size_t>{2}, std::vector<double>{1.0});
    ssa::SsaOptions o;
    o.method = "nrm";
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(mm.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("region") != std::string::npos);
    CHECK(msg.find("fcrBuf") != std::string::npos);

    // `default` ANSWERS IT, on the serial engine, which censors the path under
    // DROP exactly as SolverCTMC censors its space. That is the reference's own
    // ladder and the reason the refusal above had to be asked for by name.
    ssa::SsaOptions dflt;
    dflt.samples = 200;
    const ssa::SsaSolution r = ssa::solver_ssa(mm.get_struct(), dflt);
    CHECK(r.method == "serial");
}

TEST_CASE("ssa nrm: multi-server DPS refuses by name") {
    // State.afterEventStation rejects it in the reference too, so this is a
    // shared refusal and not a gap in the port.
    qn::Network<double> m("dpsms");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Q", SchedStrategy::DPS);
    const std::size_t c = m.add_closed_class("C1", 3, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    m.set_number_of_servers(q, 3);
    m.set_sched_param(q, c, 1.0);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    m.link(P);
    ssa::SsaOptions o;
    o.samples = 10;
    std::string msg;
    try {
        ssa::solver_ssa(m.get_struct(), o);
    } catch (const UnsupportedError& e) {
        msg = e.what();
    }
    CHECK(msg.find("multi-server") != std::string::npos);
    CHECK(msg.find("dps") != std::string::npos);
}
