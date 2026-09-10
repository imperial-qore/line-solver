/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `solver_mva_marie_analyzer` (mva_dispatch.h), the SolverMVA analyzer for
 * Marie's iterative aggregation-decomposition of a closed network with
 * non-exponential FCFS service.
 *
 * MARIE'S METHOD IS AN APPROXIMATION, so almost nothing here may be checked
 * against an exact solver at a tight tolerance. The exceptions are the two
 * regimes where the decomposition is provably exact, and they are the strongest
 * oracles available:
 *
 *   scv == 1 everywhere. The isolation chain is the M/M/1(-m) queue, its
 *   conditional departure rate is the initial multiplier lattice min(n,m)/L,
 *   so the aggregate load-dependent solve reproduces exact product form on the
 *   first pass (pfqn_marie.h documents this and exits on the second). Marie
 *   must then equal exact MVA to solver tolerance.
 *
 *   A non-FCFS product-form discipline. PS and LCFSPR are insensitive, so the
 *   analyzer hands them an SCV of one; two models differing only in the service
 *   LAW at such a station therefore produce identical input to pfqn_marie and
 *   must agree to rounding.
 *
 * Everything else is checked as an identity that binds any correct answer --
 * population conservation, Little's law, the sign of the variability
 * correction -- never against a value read back out of this implementation.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/mva_dispatch.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace mva = line::mva;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay(Z) -> Q1 -> Q2 -> Delay, one closed class, every service exponential. */
qn::Network<double> marie_two_queue(double njobs) {
    qn::Network<double> m("marie2q");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
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

/**
 * Delay(Z=1) -> Q -> Delay, one closed class, N jobs, with the queue's service
 * law supplied by the caller. Every law used below has mean 2/3, so the models
 * differ ONLY in the SCV and their demands are identical.
 */
qn::Network<double> marie_one_queue(double njobs, SchedStrategy sched, const Dist& svc) {
    qn::Network<double> m("marie1q");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", sched);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, svc);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Erlang-2, exponential and HyperExp-2 laws, all of mean 2/3. */
Dist svc_lowvar() { return Dist::erlang(3.0, 2); }              // scv 1/2
Dist svc_exp() { return Dist::exp_rate(1.5); }                  // scv 1
Dist svc_highvar() { return Dist::hyperexp(0.5, 1.0, 3.0); }    // scv 3/2

mva::MvaOptions marie_opt() {
    mva::MvaOptions o;
    o.method = "marie";
    return o;
}

/**
 * The analyzer entered below the ladder. It reports the algorithm alongside the
 * metrics, and only the metrics are of interest here.
 */
template <class T>
mva::MvaSolution<T> marie_direct(const qn::NetworkStruct<T>& L, const mva::MvaOptions& opt) {
    return mva::solver_mva_marie_analyzer(L, opt).sol;
}

}  // namespace

TEST_CASE("mva marie: exponential service reproduces exact MVA") {
    qn::Network<double> m = marie_two_queue(3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const mva::MvaSolution<double> r = marie_direct(sn, marie_opt());

    mva::MvaOptions ex;
    ex.method = "exact";
    line::Matrix<double> init;
    const mva::AvgResult<double> e = mva::solver_mva_run_analyzer(sn, ex, init);

    // With every SCV at one the isolation chains are M/M/1 queues whose
    // conditional departure rate IS the initial multiplier lattice, so the
    // aggregate is the exact load-dependent product form and the iteration adds
    // nothing. This is the one place a tight tolerance against an exact solver
    // is legitimate; everywhere else Marie is an approximation.
    CHECK(r.method == "marie");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(r.Q(i, 0) == doctest::Approx(e.QN(i, 0)).epsilon(1e-7));
        CHECK(r.U(i, 0) == doctest::Approx(e.UN(i, 0)).epsilon(1e-7));
        CHECK(r.Tp(i, 0) == doctest::Approx(e.TN(i, 0)).epsilon(1e-7));
        CHECK(r.R(i, 0) == doctest::Approx(e.RN(i, 0)).epsilon(1e-7));
    }
    CHECK(r.X[0] == doctest::Approx(e.XN[0]).epsilon(1e-7));
}

TEST_CASE("mva marie: the closed population is conserved and Little's law holds") {
    const double N = 4.0;
    qn::Network<double> m = marie_one_queue(N, SchedStrategy::FCFS, svc_lowvar());
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mva::MvaSolution<double> r = marie_direct(sn, marie_opt());

    // No job is created or destroyed in a closed network, so the queue lengths
    // must add up to the population whatever the decomposition did. An
    // approximation is free to misplace jobs but not to lose them.
    double tot = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i) tot += r.Q(i, 0);
    CHECK(tot == doctest::Approx(N).epsilon(1e-6));

    // Little's law at each station, and per class over the whole cycle.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(r.Tp(i, 0) > 0.0);
        CHECK(r.Q(i, 0) == doctest::Approx(r.R(i, 0) * r.Tp(i, 0)).epsilon(1e-9));
    }
    CHECK(r.C[0] == doctest::Approx(N / r.X[0]).epsilon(1e-9));

    // The single server cannot be busy more than all of the time, and the delay
    // station's "utilization" is its own population.
    CHECK(r.U(1, 0) > 0.0);
    CHECK(r.U(1, 0) < 1.0);
    CHECK(r.U(0, 0) == doctest::Approx(r.Q(0, 0)).epsilon(1e-9));
    CHECK(r.iter >= 1);
}

TEST_CASE("mva marie: FCFS service variability moves the queue in the right direction") {
    const double N = 4.0;
    qn::Network<double> lo = marie_one_queue(N, SchedStrategy::FCFS, svc_lowvar());
    qn::Network<double> md = marie_one_queue(N, SchedStrategy::FCFS, svc_exp());
    qn::Network<double> hi = marie_one_queue(N, SchedStrategy::FCFS, svc_highvar());

    // The three laws share a mean, so the demands handed to Marie are identical
    // and the SCV is the only thing that differs. Reading the ordering off the
    // struct rather than asserting the fitted moments keeps this a statement
    // about the analyzer, not about the distribution constructors.
    REQUIRE(lo.get_struct().scv(1, 0) < md.get_struct().scv(1, 0));
    REQUIRE(md.get_struct().scv(1, 0) < hi.get_struct().scv(1, 0));

    const mva::MvaSolution<double> rl = marie_direct(lo.get_struct(), marie_opt());
    const mva::MvaSolution<double> rm = marie_direct(md.get_struct(), marie_opt());
    const mva::MvaSolution<double> rh = marie_direct(hi.get_struct(), marie_opt());

    // More variable FCFS service congests the station: a long service holds
    // every job behind it, which is exactly the effect the Coxian isolation
    // chain exists to capture. If the SCV never reached pfqn_marie the three
    // would coincide, so this is the check that the plumbing is live AND
    // correctly signed. Only the ORDER is asserted; the sizes of the gaps are
    // the approximation's own and no closed form pins them.
    CHECK(rl.Q(1, 0) < rm.Q(1, 0));
    CHECK(rm.Q(1, 0) < rh.Q(1, 0));
    // Throughput moves the other way, the population being fixed.
    CHECK(rl.X[0] > rm.X[0]);
    CHECK(rm.X[0] > rh.X[0]);
    // Conservation survives in all three.
    for (const mva::MvaSolution<double>* r : {&rl, &rm, &rh})
        CHECK(r->Q(0, 0) + r->Q(1, 0) == doctest::Approx(N).epsilon(1e-6));
}

TEST_CASE("mva marie: a PS station is insensitive to its service variability") {
    const double N = 4.0;
    qn::Network<double> ps_e = marie_one_queue(N, SchedStrategy::PS, svc_exp());
    qn::Network<double> ps_h = marie_one_queue(N, SchedStrategy::PS, svc_highvar());
    const mva::MvaSolution<double> re = marie_direct(ps_e.get_struct(), marie_opt());
    const mva::MvaSolution<double> rh = marie_direct(ps_h.get_struct(), marie_opt());

    // PS is a product-form insensitive discipline, so the analyzer overrides its
    // SCV with one and the two models present pfqn_marie with literally the same
    // demands. Anything but agreement here would mean the FCFS-only guard leaked
    // a service law into a station whose answer cannot depend on it.
    REQUIRE(ps_h.get_struct().scv(1, 0) > 1.0);
    CHECK(rh.Q(1, 0) == doctest::Approx(re.Q(1, 0)).epsilon(1e-12));
    CHECK(rh.X[0] == doctest::Approx(re.X[0]).epsilon(1e-12));

    // And, being exponential-equivalent and single chain, that common answer is
    // the exact product-form one.
    line::Matrix<double> init;
    mva::MvaOptions ex;
    ex.method = "exact";
    const mva::AvgResult<double> e = mva::solver_mva_run_analyzer(ps_e.get_struct(), ex, init);
    CHECK(re.Q(1, 0) == doctest::Approx(e.QN(1, 0)).epsilon(1e-7));
}

TEST_CASE("mva marie: the dispatch and the runner both route method 'marie' here") {
    qn::Network<double> m = marie_two_queue(3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const mva::MvaSolution<double> direct = marie_direct(sn, marie_opt());

    line::Matrix<double> init;
    const mva::DispatchResult<double> viaDispatch = mva::mva_dispatch(sn, marie_opt(), init);
    CHECK(viaDispatch.actualmethod == "marie");

    // The ladder must land on THIS analyzer and not on an AMVA that would also
    // return plausible numbers for a closed network. The direct call above is
    // the same function the dispatch reaches, so what the comparison pins is the
    // ROUTING: a branch that chose differently would not reproduce it.
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(viaDispatch.sol.Q(i, 0) == doctest::Approx(direct.Q(i, 0)).epsilon(1e-12));
        CHECK(viaDispatch.sol.U(i, 0) == doctest::Approx(direct.U(i, 0)).epsilon(1e-12));
        CHECK(viaDispatch.sol.R(i, 0) == doctest::Approx(direct.R(i, 0)).epsilon(1e-12));
    }
    CHECK(viaDispatch.sol.X[0] == doctest::Approx(direct.X[0]).epsilon(1e-12));
    CHECK(viaDispatch.sol.iter == direct.iter);

    // The runner adds the metric filter and the derived ArvR/ResidT columns on
    // top of the analyzer, so its table must still carry the analyzer's own
    // numbers -- a branch that reached a different algorithm would not.
    const mva::AvgResult<double> table = mva::solver_mva_run_analyzer(sn, marie_opt(), init);
    CHECK(table.actualmethod == "marie");
    CHECK(table.method == "marie");
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(table.QN(i, 0) == doctest::Approx(direct.Q(i, 0)).epsilon(1e-12));
        CHECK(table.TN(i, 0) == doctest::Approx(direct.Tp(i, 0)).epsilon(1e-12));
    }
    // One visit per cycle at each station, so residence and response coincide.
    CHECK(table.WN(1, 0) == doctest::Approx(table.RN(1, 0)).epsilon(1e-9));
}

TEST_CASE("mva marie: an open model is refused by name") {
    // Marie's aggregate is a closed load-dependent product form over a finite
    // population lattice, so an unbounded chain has nothing to iterate over.
    qn::Network<double> m("open");
    const std::size_t s = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(s, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(s, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);

    CHECK_THROWS_AS(marie_direct(m.get_struct(), marie_opt()),
                    line::UnsupportedError);
}

TEST_CASE("mva marie: an unsupported discipline is refused by name") {
    // SIRO is neither service-sensitive in the way the Coxian isolation assumes
    // nor insensitive the way PS is, so it is refused rather than silently
    // solved as one of the two.
    qn::Network<double> m = marie_one_queue(3.0, SchedStrategy::SIRO, svc_exp());
    CHECK_THROWS_AS(marie_direct(m.get_struct(), marie_opt()),
                    line::UnsupportedError);
}

TEST_CASE("mva marie: a multichain multiserver station is refused by name") {
    // The multiclass isolation chain has no multiserver form, in the reference
    // or here, so the combination is named rather than approximated by one
    // server. The same station with a single chain is accepted, which is what
    // makes this a refusal of the COMBINATION.
    qn::Network<double> m("marie-ms");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 2.0, d);
    m.set_service(d, c1, Dist::exp_rate(1.0));
    m.set_service(d, c2, Dist::exp_rate(1.0));
    m.set_service(q, c1, Dist::exp_rate(2.0));
    m.set_service(q, c2, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    m.set_number_of_servers(q, 2.0);
    REQUIRE(m.get_struct().nchains > 1);
    CHECK_THROWS_AS(marie_direct(m.get_struct(), marie_opt()),
                    line::UnsupportedError);

    qn::Network<double> one = marie_one_queue(3.0, SchedStrategy::FCFS, svc_exp());
    one.set_number_of_servers(2, 2.0);
    REQUIRE(one.get_struct().nchains == 1);
    CHECK_NOTHROW(marie_direct(one.get_struct(), marie_opt()));
}

TEST_CASE("mva marie: exact rational arithmetic is refused by name") {
    // The Coxian fit needs a square root and the outer loop stops on a
    // tolerance, so there is no exact-arithmetic answer to give. The branch is
    // discarded at compile time and the model refused, which is what lets
    // mva_dispatch<Rational> compile at all.
    qn::Network<line::Rational> m("marie-exact");
    using RDist = line::lang::Distrib<line::Rational>;
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, RDist::exp_rate(line::num_traits<line::Rational>::from_int(1)));
    m.set_service(q, c, RDist::exp_rate(line::num_traits<line::Rational>::from_int(2)));
    qn::RoutingMatrix<line::Rational> P;
    P.set(d, q, line::num_traits<line::Rational>::from_int(1));
    P.set(q, d, line::num_traits<line::Rational>::from_int(1));
    m.link(P);

    CHECK_THROWS_AS(marie_direct(m.get_struct(), marie_opt()),
                    line::UnsupportedError);
}
