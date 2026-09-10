/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * The second-order fluid methods (`minnormal`, `refined`), the Ko-Pender joint
 * mean/covariance limit (`kp`), and the two closing-drift branches they brought
 * with them: load dependence and GPS.
 *
 * EVERY EXPECTED VALUE IS MATLAB `SolverFLD(model,'method',...)` at LINE 3.0.7,
 * printed to nine significant digits. The `minnormal` rows were RE-RECORDED on
 * 2026-09-01 against MATLAB R2026a: the closure alternation and its inner mean
 * solve were tightened from CoarseTol to mom_tol = 1e-6 (solver_fluid_moments.m,
 * minnormal.py, fluid_moments.h), which moves every converged second-order
 * answer by 2e-5 to 3e-5 -- past the 1e-5 asserted here. This port and native
 * python reproduce the new MATLAB values to 4e-8 and 2.3e-6 respectively, which
 * is what licenses the re-record; the `closing` rows did not move. Where the exact answer is also known --
 * a birth-death chain or an M/G/inf -- it is asserted as well, because agreeing
 * with the reference to 1e-6 and with the truth to 2% are different claims and
 * the second one is what says the closure is doing its job.
 */
#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/fluid/fluid_kp.h"
#include "line/solvers/fluid/fluid_moments.h"
#include "line/solvers/fluid/fluid_runner.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& name) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        if (sn.stations[i].name == name) return i;
    throw line::InputError("no station named " + name);
}

/** Delay(Z=1) -> Queue(PS, mu=1, c) with N jobs, the reference's own sweep model. */
qn::Network<double> closed_ps(std::size_t N, double c) {
    qn::Network<double> m("mn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    m.set_number_of_servers(q, c);
    const std::size_t k = m.add_closed_class("Class1", static_cast<double>(N), d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q, k, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q, 1.0);
    P.set(k, k, q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fluid minnormal: the closure corrects the mean where min() bends") {
    // Delay(Z=1) -> Queue(PS, mu=1, c=2), N=6. The first-order closure replaces
    // E[min(X,2)] by min(E[X],2) and saturates the server exactly, returning
    // Q1 = 2 and Q2 = 4; the exact stationary mean of the birth-death chain with
    // birth (N-n) and death min(n,2) is 1.95137. The Gaussian closure recovers
    // most of that gap, which is the whole point of the method.
    qn::Network<double> m = closed_ps(6, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t d = station_of(sn, "Delay"), q = station_of(sn, "Queue1");

    fluid::FluidOptions o;
    o.method = "minnormal";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(s.method == "minnormal");
    CHECK(s.QN(d, 0) == doctest::Approx(1.96063265).epsilon(1e-5));
    CHECK(s.QN(q, 0) == doctest::Approx(4.03936735).epsilon(1e-5));
    CHECK(s.UN(d, 0) == doctest::Approx(1.96063265).epsilon(1e-5));
    CHECK(s.UN(q, 0) == doctest::Approx(0.980316327).epsilon(1e-5));
    CHECK(s.TN(d, 0) == doctest::Approx(1.96063265).epsilon(1e-5));
    CHECK(s.TN(q, 0) == doctest::Approx(1.96063265).epsilon(1e-5));

    // The first-order closure on the same model, for contrast: it saturates.
    // Read at 1e-5, not 1e-6: the fixed point is exactly (2,4) and the reference
    // returns it to nine digits, but this port restarts LSODA once per pass and
    // at the default tol = 1e-4 the restart error accumulates to 7e-6 over the
    // 200 passes the reference also runs. See the header of solver_fluid.h.
    fluid::FluidOptions oc;
    oc.method = "closing";
    const fluid::FluidSolution sc = fluid::solver_fluid_run_analyzer(sn, oc);
    CHECK(sc.QN(d, 0) == doctest::Approx(2.0).epsilon(1e-5));
    CHECK(sc.QN(q, 0) == doctest::Approx(4.0).epsilon(1e-5));

    // And against the exact chain: the closure is 5x closer than the first-order
    // method, which is the claim the method exists to make.
    const double exact = 1.95137;
    CHECK(std::fabs(s.QN(d, 0) - exact) < std::fabs(sc.QN(d, 0) - exact) / 4.0);
}

TEST_CASE("fluid minnormal: the covariance is reported, and only by the closures") {
    // getMoments: the state-level covariance, the per-station-class queue-length
    // variance and its square root. QVar = sigma2 here because the station carries
    // one class of one phase.
    qn::Network<double> m = closed_ps(6, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t d = station_of(sn, "Delay"), q = station_of(sn, "Queue1");

    fluid::FluidOptions o;
    o.method = "minnormal";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    REQUIRE(s.has_moments);
    CHECK(s.moments.QVar(d, 0) == doctest::Approx(1.83872243).epsilon(1e-5));
    CHECK(s.moments.QVar(q, 0) == doctest::Approx(1.83872243).epsilon(1e-5));
    CHECK(s.moments.QStd(q, 0) == doctest::Approx(std::sqrt(1.83872243)).epsilon(1e-5));
    CHECK(s.moments.sigma2[q] == doctest::Approx(1.83872243).epsilon(1e-5));
    // Population is conserved, so the two coordinates are perfectly negatively
    // correlated and the total variance is zero: Sigma is not diagonal.
    CHECK(s.moments.Sigma.rows() == 2);
    CHECK(s.moments.Sigma(0, 1) == doctest::Approx(-s.moments.Sigma(0, 0)).epsilon(1e-6));

    // A first-order method computes no second moment at all, and says so rather
    // than reporting a variance of zero.
    fluid::FluidOptions oc;
    oc.method = "closing";
    CHECK_FALSE(fluid::solver_fluid_run_analyzer(sn, oc).has_moments);
}

TEST_CASE("fluid refined: the 1/N correction is added to the MEAN-FIELD point") {
    // 'refined' recomputes its base point with the first-order closure and adds
    // Gast's O(1/N) term to THAT, because the Gaussian fixed point already resums
    // it. The two therefore do not agree, and must not.
    qn::Network<double> m = closed_ps(6, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t d = station_of(sn, "Delay"), q = station_of(sn, "Queue1");

    fluid::FluidOptions o;
    o.method = "refined";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(s.method == "refined");
    CHECK(s.QN(d, 0) == doctest::Approx(1.91433883).epsilon(1e-5));
    CHECK(s.QN(q, 0) == doctest::Approx(4.08566117).epsilon(1e-5));
    CHECK(s.UN(q, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(s.TN(d, 0) == doctest::Approx(1.91433883).epsilon(1e-5));
    REQUIRE(s.has_moments);
    CHECK(s.moments.refinement.size() == 2);
}

TEST_CASE("fluid default: it resolves to minnormal, not to matrix") {
    // The reference's `default` prefers the second-order closure wherever
    // `fluid_minnormal_applicable` accepts the model, which is what makes the
    // ported default the same METHOD as the reference's and not merely the same
    // name. Before this the port answered `matrix` here.
    qn::Network<double> m = closed_ps(6, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    std::string why;
    CHECK(fluid::detail::fluid_resolve_method(sn, "default", fluid::FluidOptions(), why) ==
          "minnormal");
    CHECK(why.empty());
    CHECK(fluid::solver_fluid_run_analyzer(sn, fluid::FluidOptions()).method == "minnormal");

    // And it declines with a REASON when the model is outside the closure's reach.
    // A phase-type source is the canonical case: those coordinates track the phase
    // of one arrival process, not a population.
    qn::Network<double> mo("open");
    const std::size_t src = mo.add_source("Source");
    const std::size_t qq = mo.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t snk = mo.add_sink("Sink");
    const std::size_t oc = mo.add_open_class("Class1");
    mo.set_service(src, oc, D::erlang(2.0, 2));
    mo.set_service(qq, oc, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(oc, oc, src, qq, 1.0);
    P.set(oc, oc, qq, snk, 1.0);
    mo.link(P);
    const qn::NetworkStruct<double>& sno = mo.get_struct();
    CHECK(fluid::detail::fluid_resolve_method(sno, "default", fluid::FluidOptions(), why) ==
          "matrix");
    CHECK(why.find("non-Poisson") != std::string::npos);
}

TEST_CASE("fluid minnormal: an open model projects the source pool out") {
    // Source(Exp(1)) -> Queue(PS, mu=0.8, c=2) -> Sink. The covariance lives on
    // the queue coordinate only: the EXT coordinate is a normalisation constant,
    // not a job count, and building diffusion over it would invent noise for a
    // direction with no population.
    qn::Network<double> m("mnopen");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t oc = m.add_open_class("Class1");
    m.set_number_of_servers(q, 2.0);
    m.set_service(src, oc, D::exp_rate(1.0));
    m.set_service(q, oc, D::exp_rate(0.8));
    qn::RoutingMatrix<double> P;
    P.set(oc, oc, src, q, 1.0);
    P.set(oc, oc, q, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t iq = station_of(sn, "Queue1"), is = station_of(sn, "Source");

    fluid::FluidOptions o;
    o.method = "minnormal";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(s.QN(iq, 0) == doctest::Approx(1.69397477).epsilon(1e-5));
    CHECK(s.TN(is, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(s.UN(iq, 0) == doctest::Approx(0.625).epsilon(1e-5));

    // `refined` keeps the closed-model restriction: its correction is solved over
    // the FULL state and would perturb the source pool mass.
    fluid::FluidOptions orf;
    orf.method = "refined";
    CHECK_THROWS_AS(fluid::solver_fluid_run_analyzer(sn, orf), line::Error);
}

TEST_CASE("fluid: load dependence reaches the drift, and only the closing family") {
    // Delay(Z=1) -> Queue(PS, mu=1, alpha=[1 1.5 2 2]), N=4. alpha multiplies the
    // scheduling share in the drift AND the capacity term in the metric reader,
    // and the utilization normalises by the PEAK scaling. Exact CTMC: Q = [1.6 2.4].
    qn::Network<double> m("lld");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t k = m.add_closed_class("Class1", 4.0, d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q, k, D::exp_rate(1.0));
    std::vector<double> alpha;
    alpha.push_back(1.0);
    alpha.push_back(1.5);
    alpha.push_back(2.0);
    alpha.push_back(2.0);
    m.set_load_dependence(q, alpha);
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q, 1.0);
    P.set(k, k, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), iq = station_of(sn, "Queue1");

    fluid::FluidOptions oc;
    oc.method = "closing";
    const fluid::FluidSolution sc = fluid::solver_fluid_run_analyzer(sn, oc);
    CHECK(sc.QN(id, 0) == doctest::Approx(1.66666667).epsilon(1e-5));
    CHECK(sc.QN(iq, 0) == doctest::Approx(2.33333333).epsilon(1e-5));
    CHECK(sc.TN(iq, 0) == doctest::Approx(1.66666667).epsilon(1e-5));
    CHECK(sc.UN(iq, 0) == doctest::Approx(0.833333333).epsilon(1e-5));

    fluid::FluidOptions om;
    om.method = "minnormal";
    const fluid::FluidSolution sm = fluid::solver_fluid_run_analyzer(sn, om);
    CHECK(sm.QN(id, 0) == doctest::Approx(1.58897870).epsilon(1e-5));
    CHECK(sm.QN(iq, 0) == doctest::Approx(2.41102130).epsilon(1e-5));
    CHECK(sm.UN(iq, 0) == doctest::Approx(0.794489294).epsilon(1e-5));

    // Every other method builds its drift independently and would silently return
    // the answer for alpha == 1, so it is refused by name.
    fluid::FluidOptions ox;
    ox.method = "matrix";
    CHECK_THROWS_AS(fluid::solver_fluid_run_analyzer(sn, ox), line::Error);
}

TEST_CASE("fluid GPS: the backlog closure is the mechanism, and minnormal only") {
    // Delay -> Queue(GPS, weights [1 3]), one job per class. GPS splits the server
    // by weight among the BACKLOGGED classes: with continuous mass every class is
    // always backlogged, so a first-order closure prices the share at the constant
    // w_k/sum_j w_j regardless of load, which is why the featset gate refuses GPS
    // for every method except `minnormal`.
    qn::Network<double> m("gps");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::GPS);
    const std::size_t a1 = m.add_closed_class("Class1", 1.0, d);
    const std::size_t a2 = m.add_closed_class("Class2", 1.0, d);
    m.set_service(d, a1, D::exp_rate(1.0));
    m.set_service(d, a2, D::exp_rate(1.0));
    m.set_service(q, a1, D::exp_rate(2.0));
    m.set_service(q, a2, D::exp_rate(1.0));
    m.set_sched_param(q, a1, 1.0);
    m.set_sched_param(q, a2, 3.0);
    qn::RoutingMatrix<double> P;
    P.set(a1, a1, d, q, 1.0);
    P.set(a1, a1, q, d, 1.0);
    P.set(a2, a2, d, q, 1.0);
    P.set(a2, a2, q, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t id = station_of(sn, "Delay"), iq = station_of(sn, "Queue1");

    fluid::FluidOptions o;
    o.method = "minnormal";
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(s.QN(id, 0) == doctest::Approx(0.556724182).epsilon(1e-4));
    CHECK(s.QN(id, 1) == doctest::Approx(0.465455957).epsilon(1e-4));
    CHECK(s.QN(iq, 0) == doctest::Approx(0.443275818).epsilon(1e-4));
    CHECK(s.QN(iq, 1) == doctest::Approx(0.534544043).epsilon(1e-4));
    CHECK(s.TN(iq, 0) == doctest::Approx(0.556724182).epsilon(1e-4));
    CHECK(s.TN(iq, 1) == doctest::Approx(0.465455957).epsilon(1e-4));

    // The gate refuses GPS for the first-order methods by naming the discipline.
    fluid::FluidOptions oc;
    oc.method = "closing";
    CHECK_THROWS_AS(fluid::solver_fluid_run_analyzer(sn, oc), line::Error);
    oc.method = "matrix";
    CHECK_THROWS_AS(fluid::solver_fluid_run_analyzer(sn, oc), line::Error);
    // But the initial condition DOES decode a GPS station: the refusal is the
    // share's, not the state encoding's, and `solver_fluid_initsol.m` lists GPS.
    CHECK_NOTHROW(fluid::fluid_initsol(sn));
}

TEST_CASE("fluid kp: the mean and the covariance of a MAP/Erlang/inf network") {
    // Source(MAP) -> Delay(Erlang-2, mean 0.5) -> Sink. The arrival stream is
    // non-renewal, which is the point of the method: the closing drift would
    // return mass through the stationary arrival-instant vector and lose the
    // autocorrelation. At an infinite server the rate functions are affine, so BOTH
    // equations close exactly: E[N] = lambda*E[S] = 0.5 and, the arrivals being
    // Poisson-like at this order, Var[N] = 0.5 as well.
    qn::Network<double> m("kp");
    const std::size_t src = m.add_source("Source");
    const std::size_t inf = m.add_delay("Delay1");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t oc = m.add_open_class("Class1");
    Matrix<double> A0(2, 2, 0.0), A1(2, 2, 0.0);
    A0(0, 0) = -2.0;
    A0(0, 1) = 1.0;
    A0(1, 0) = 0.5;
    A0(1, 1) = -1.5;
    A1(0, 0) = 0.8;
    A1(0, 1) = 0.2;
    A1(1, 0) = 0.6;
    A1(1, 1) = 0.4;
    m.set_service(src, oc, D::map_dist(A0, A1, lang::ProcessType::MAP));
    m.set_service(inf, oc, D::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(oc, oc, src, inf, 1.0);
    P.set(oc, oc, inf, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t is = station_of(sn, "Source"), id = station_of(sn, "Delay1");

    fluid::FluidOptions o;
    o.method = "kp";
    o.timespan_end = 20.0;
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(s.method == "kp");
    CHECK(s.QN(id, 0) == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(s.TN(id, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.TN(is, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(s.UN(id, 0) == doctest::Approx(0.5).epsilon(1e-4));
    REQUIRE(s.has_moments);
    CHECK(s.moments.QVar(id, 0) == doctest::Approx(0.5).epsilon(1e-4));
    // The state covariance carries the arrival phase indicator too, whose diagonal
    // is theta*(1-theta) at the stationary phase distribution: 4x4 here, arrival
    // phases before service phases.
    CHECK(s.moments.Sigma.rows() == 4);

    // getTranAvgVar: the same covariance along the trajectory.
    const fluid::FluidKpTransient tr = fluid::solver_fluid_tran_avg_var(sn, o);
    REQUIRE(tr.t.size() > 2);
    CHECK(tr.t.back() == doctest::Approx(20.0).epsilon(1e-9));
    CHECK(tr.QVar.back()(id, 0) == doctest::Approx(0.5).epsilon(1e-4));
    // It starts empty, so the variance grows from zero at the served station.
    CHECK(tr.QVar.front()(id, 0) == doctest::Approx(0.0).epsilon(1e-9));

    // Only `kp` has a second moment along the trajectory.
    fluid::FluidOptions ocl;
    ocl.method = "closing";
    CHECK_THROWS_AS(fluid::solver_fluid_tran_avg_var(sn, ocl), line::UnsupportedError);
}

TEST_CASE("fluid kp: a closed class is refused, a finite server is not") {
    // The method analyses an OPEN network: a closed class has no arrival process to
    // modulate, so it is refused by name rather than answered with an empty
    // u-block. A finite-server station IS admitted, through the fluid min(x,c),
    // and only degrades the covariance to a linear noise approximation.
    qn::Network<double> mc = closed_ps(4, 1.0);
    fluid::FluidOptions o;
    o.method = "kp";
    CHECK_THROWS_AS(fluid::solver_fluid_kp(mc.get_struct(), o), line::UnsupportedError);

    qn::Network<double> m("kp2");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t oc = m.add_open_class("Class1");
    m.set_number_of_servers(q, 2.0);
    m.set_service(src, oc, D::exp_rate(1.0));
    m.set_service(q, oc, D::exp_rate(1.5));
    qn::RoutingMatrix<double> P;
    P.set(oc, oc, src, q, 1.0);
    P.set(oc, oc, q, snk, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t iq = station_of(sn, "Queue1");
    o.timespan_end = 40.0;
    const fluid::FluidSolution s = fluid::solver_fluid_kp(sn, o);
    CHECK(s.QN(iq, 0) == doctest::Approx(0.666666667).epsilon(1e-3));
    CHECK(s.TN(iq, 0) == doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("fluid closures: the invariants each one exists to preserve") {
    // The share closure is exactly capacity conserving, because
    // sum_j Cov(u_j,v) = Var(v) makes the two correction terms cancel in the sum.
    // A work-conserving discipline that lost that identity would leak capacity.
    std::vector<double> x, w;
    x.push_back(2.0);
    x.push_back(3.0);
    w.push_back(1.0);
    w.push_back(4.0);
    Matrix<double> C(2, 2, 0.0);
    C(0, 0) = 0.7;
    C(1, 1) = 1.3;
    C(0, 1) = C(1, 0) = -0.4;
    const fluid::ShareValue sh = fluid::fluid_share_closure(x, w, C, true);
    CHECK(sh.s[0] + sh.s[1] == doctest::Approx(1.0).epsilon(1e-12));

    // With no covariance it is the plug-in share, and the Jacobian rows sum to
    // zero: moving mass between coordinates cannot change the total share.
    const fluid::ShareValue s0 = fluid::fluid_share_closure(x, w, Matrix<double>(0, 0, 0.0), true);
    CHECK(s0.s[0] == doctest::Approx(2.0 / 14.0).epsilon(1e-12));
    CHECK(s0.ds(0, 0) + s0.ds(1, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // The min closure collapses to min(n,c) at zero variance, and its derivative
    // at the kink is the RIGHT derivative of min, i.e. zero.
    CHECK(fluid::fluid_min_closure(3.0, 2.0, 0.0).h == doctest::Approx(2.0));
    CHECK(fluid::fluid_min_closure(3.0, 2.0, 0.0).dh == doctest::Approx(0.0));
    CHECK(fluid::fluid_min_closure(2.0, 2.0, 0.0).dh == doctest::Approx(0.0));
    CHECK(fluid::fluid_min_closure(1.0, 2.0, 0.0).dh == doctest::Approx(1.0));
    // At the kink with variance the closure is strictly below the kink value, by
    // theta/sqrt(2 pi): that gap IS the correction the method applies.
    const double s2 = 4.0;
    CHECK(fluid::fluid_min_closure(2.0, 2.0, s2).h ==
          doctest::Approx(2.0 - std::sqrt(s2) / std::sqrt(2.0 * 3.14159265358979323846))
              .epsilon(1e-9));

    // The GPS shares sum to 1 - P(station empty), NOT to one: that is how the idle
    // single server is represented. With both classes surely backlogged they sum
    // to one and reduce to the weights.
    std::vector<double> xk, wk, vk;
    xk.push_back(3.0);
    xk.push_back(1.0);
    wk.push_back(1.0);
    wk.push_back(3.0);
    vk.push_back(0.0);
    vk.push_back(0.0);
    const fluid::ShareValue g0 = fluid::fluid_gps_share(xk, wk, vk, false);
    CHECK(g0.s[0] == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(g0.s[1] == doctest::Approx(0.75).epsilon(1e-12));
    vk[0] = 1.0;
    vk[1] = 1.0;
    const fluid::ShareValue g1 = fluid::fluid_gps_share(xk, wk, vk, false);
    CHECK(g1.s[0] + g1.s[1] < 1.0);
    CHECK(g1.s[0] + g1.s[1] > 0.5);

    // The load-dependent scaling interpolates the integer table and clamps both
    // tails, matching the clamping the CTMC applies.
    std::vector<double> alpha;
    alpha.push_back(1.0);
    alpha.push_back(1.5);
    alpha.push_back(2.0);
    CHECK(fluid::fluid_lld_scaling(alpha, 0.4).h == doctest::Approx(1.0));
    CHECK(fluid::fluid_lld_scaling(alpha, 1.5).h == doctest::Approx(1.25));
    CHECK(fluid::fluid_lld_scaling(alpha, 1.5).dh == doctest::Approx(0.5));
    CHECK(fluid::fluid_lld_scaling(alpha, 9.0).h == doctest::Approx(2.0));
    CHECK(fluid::fluid_lld_scaling(std::vector<double>(), 3.0).h == doctest::Approx(1.0));

    // The load-dependent capacity closure integrates psi EXACTLY against the
    // normal density, so at zero variance it is psi itself: min(n,c)*alpha(n).
    CHECK(fluid::fluid_capacity_closure(1.5, 4.0, 0.0, alpha, false).h ==
          doctest::Approx(1.5 * 1.25).epsilon(1e-12));
    // And it is a genuine expectation: with variance it differs from psi(n).
    CHECK(fluid::fluid_capacity_closure(1.5, 4.0, 0.25, alpha, false).h !=
          doctest::Approx(1.5 * 1.25).epsilon(1e-6));
}
