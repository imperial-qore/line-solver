/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Immediate-transition elimination and the stiff arm of the fluid solver.
 *
 * WHY ELIMINATION CAN BE CHECKED SHARPLY AND INTEGRATION CANNOT. The stochastic
 * complement is an identity: the generator it returns is the one an observer
 * who watches only the retained states sees, so its stationary law is the
 * original's conditioned on those states and renormalized. That is a theorem
 * about the construction, it holds at any rate ratio, and it can be checked
 * against a stationary solve of the original generator without any reference
 * value at all. Every elimination test below is that identity in one form or
 * another.
 *
 * The integrator, by contrast, can only be checked against ANOTHER integrator
 * or against a problem with a closed-form solution, and only to the tolerance
 * both were asked for. Both are used: a linear drift whose fixed point is known
 * exactly, and LSODA, which reaches the same trajectory by a completely
 * different method.
 *
 * THE MODELS ARE ALL INFINITE SERVERS ON PURPOSE. The elimination reads the
 * drift as `dx/dt = x W`, which is what the fluid drift IS wherever the
 * state-dependent factor is x itself. At an infinite server it always is, so
 * these models are exactly the regime in which the reconstructed generator is
 * the drift and the complement is exact rather than nearly so. A saturated
 * queue would test the approximation instead of the algebra.
 */
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/util/lsoda.h"

using namespace line;
using D = lang::Distrib<double>;

namespace {

/**
 * D1 -> Dmid -> D2 -> D1, three infinite servers, one closed class.
 *
 * `mu_mid` at `GlobalConstants::Immediate` is what an Immediate distribution
 * contributes to the drift; the rate is used rather than the process tag
 * because that is what the detector keys on and what makes the system stiff.
 */
qn::Network<double> delay_cycle3(double n, double mu_mid) {
    qn::Network<double> m("fluid_stiff_cycle3");
    const std::size_t d1 = m.add_delay("D1");
    const std::size_t dm = m.add_delay("Dmid");
    const std::size_t d2 = m.add_delay("D2");
    const std::size_t c = m.add_closed_class("C1", n, d1);
    m.set_service(d1, c, D::exp_rate(1.0));
    m.set_service(dm, c, D::exp_rate(mu_mid));
    m.set_service(d2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d1, dm, 1.0);
    P.set(c, c, dm, d2, 1.0);
    P.set(c, c, d2, d1, 1.0);
    m.link(P);
    return m;
}

/** D1 <-> D2, two infinite servers, one closed class. */
qn::Network<double> delay_pair(double n) {
    qn::Network<double> m("fluid_stiff_pair");
    const std::size_t d1 = m.add_delay("D1");
    const std::size_t d2 = m.add_delay("D2");
    const std::size_t c = m.add_closed_class("C1", n, d1);
    m.set_service(d1, c, D::exp_rate(1.0));
    m.set_service(d2, c, D::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d1, d2, 1.0);
    P.set(c, c, d2, d1, 1.0);
    m.link(P);
    return m;
}

/** Close a rate matrix into a generator by putting minus the row sum on it. */
line::Matrix<double> as_generator(line::Matrix<double> W) {
    for (std::size_t i = 0; i < W.rows(); ++i) {
        double s = 0;
        for (std::size_t j = 0; j < W.cols(); ++j)
            if (j != i) s += W(i, j);
        W(i, i) = -s;
    }
    return W;
}

/** Integrate a fluid drift to `t_end` with LSODA, from `x0`. */
std::vector<double> settle(const fluid::FluidOdeSystem& sys, const std::vector<double>& x0,
                           double t_end) {
    LsodaOptions lopt;
    lopt.rtol = 1e-10;
    lopt.atol = 1e-12;
    return lsoda_integrate(fluid::fluid_drift(sys), x0, std::vector<double>{0.0, t_end}, lopt)
        .final_state();
}

}  // namespace

TEST_CASE("fluid immediate: eliminating an immediate transition leaves the fixed point alone") {
    // D1(rate 1) -> Dmid(rate 1e8) -> D2(rate 3) -> D1, six jobs. Every station
    // is an infinite server, so the drift is linear and its fixed point is the
    // stationary mean exactly: mass proportional to the service demand 1/mu.
    const double N = 6.0, imm = lang::GlobalConstants::Immediate;
    qn::Network<double> m = delay_cycle3(N, imm);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(sn);
    REQUIRE(sys.layout.nstates == 3);

    const fluid::FluidImmediateResult r = fluid::fluid_eliminate_immediate(sys);
    REQUIRE(r.eliminated);
    CHECK(r.fallback.empty());
    CHECK(r.n_immediate == 1);
    // Exactly the coordinate that SOURCES the immediate transition is removed.
    // The state it leads into is entered instantaneously but left at rate 3, so
    // it is timed and must survive: that is where the mass actually sits.
    REQUIRE(r.state_map.size() == 2);
    CHECK(r.state_map[0] == sys.layout.qidx[0][0]);
    CHECK(r.state_map[1] == sys.layout.qidx[2][0]);

    // The reduced drift on the two surviving coordinates is the two-delay cycle
    // that remains once the instantaneous hop is contracted away: mass splits
    // in proportion to 1/mu, which is the exact stationary mean of a linear
    // drift and owes nothing to the integrator beyond its tolerance.
    std::vector<double> x0(sys.layout.nstates, 0.0);
    x0[sys.layout.qidx[0][0]] = N;
    const std::vector<double> xr = settle(r.sys, x0, 60.0);
    const double d1 = 1.0, d2 = 1.0 / 3.0;
    CHECK(xr[sys.layout.qidx[0][0]] == doctest::Approx(N * d1 / (d1 + d2)).epsilon(1e-6));
    CHECK(xr[sys.layout.qidx[2][0]] == doctest::Approx(N * d2 / (d1 + d2)).epsilon(1e-6));
    // The eliminated coordinate receives no transitions at all, so a state that
    // starts empty there stays empty; that is what makes the reduced system a
    // drop-in on the original layout.
    CHECK(xr[sys.layout.qidx[1][0]] == doctest::Approx(0.0).epsilon(1e-12));

    // And against the UNREDUCED system, whose fixed point is the same law with
    // the immediate station included. Its demand is 1/1e8, so the two answers
    // differ by a relative 1e-8: the mass the immediate state holds and the
    // complement conditions away. This is the only gap elimination introduces,
    // and it shrinks as the immediate rate grows.
    const double dm = 1.0 / imm, tot = d1 + dm + d2;
    CHECK(xr[sys.layout.qidx[0][0]] == doctest::Approx(N * d1 / tot).epsilon(1e-6));
    CHECK(std::fabs(N * d1 / (d1 + d2) - N * d1 / tot) < 1e-6);
}

TEST_CASE("fluid immediate: the reduced drift still conserves mass") {
    // Every event moves one unit of mass from one coordinate to another, so the
    // drift sums to zero at EVERY state, not only at the fixed point. A
    // reconstruction that dropped the -1 of a jump, or that pointed it at the
    // wrong coordinate, breaks this at the first state tried and would
    // otherwise show up only as a slow leak over a long integration.
    qn::Network<double> m = delay_cycle3(6.0, lang::GlobalConstants::Immediate);
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(m.get_struct());
    const fluid::FluidImmediateResult r = fluid::fluid_eliminate_immediate(sys);
    REQUIRE(r.eliminated);

    const std::function<void(double, const double*, double*)> f = fluid::fluid_drift(r.sys);
    const std::size_t n = sys.layout.nstates;
    const double probes[3][3] = {{6.0, 0.0, 0.0}, {1.5, 0.0, 4.5}, {2.0, 1.0, 3.0}};
    for (std::size_t p = 0; p < 3; ++p) {
        std::vector<double> x(probes[p], probes[p] + 3), dx(n, 0.0);
        f(0.0, x.data(), dx.data());
        double s = 0;
        for (std::size_t i = 0; i < n; ++i) s += dx[i];
        CHECK(s == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
    }
}

TEST_CASE("fluid immediate: elimination is inert when no transition is immediate") {
    // The strongest check available on the detector: a model whose rates are
    // all of order one must come back untouched, event for event. A threshold
    // that keyed on the model's own rates rather than on the Immediate constant
    // would start eliminating the fastest ordinary transition here.
    qn::Network<double> m = delay_pair(4.0);
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(m.get_struct());
    const fluid::FluidImmediateResult r = fluid::fluid_eliminate_immediate(sys);

    CHECK(!r.eliminated);
    CHECK(r.n_immediate == 0);
    CHECK(r.fallback.empty());  // nothing was attempted, so there is nothing to warn about
    REQUIRE(r.state_map.size() == sys.layout.nstates);
    for (std::size_t i = 0; i < r.state_map.size(); ++i) CHECK(r.state_map[i] == i);
    REQUIRE(r.sys.events.size() == sys.events.size());
    for (std::size_t e = 0; e < sys.events.size(); ++e) {
        CHECK(r.sys.events[e].minus == sys.events[e].minus);
        CHECK(r.sys.events[e].plus == sys.events[e].plus);
        CHECK(r.sys.events[e].rate_base == doctest::Approx(sys.events[e].rate_base).epsilon(1e-15));
    }
    // The departure count survives an inert pass; it is only a reduced system
    // that cannot say which of its edges are departures.
    CHECK(r.sys.n_departures == sys.n_departures);
}

TEST_CASE("fluid immediate: the matrix elimination is the censored chain") {
    // A four-state ring, one of whose states is left at 1e9 while the others
    // run at rates of order one. The stochastic complement's defining property
    // is that its stationary law is the original's, conditioned on the retained
    // states and renormalized; both sides are computed here by an independent
    // stationary solve, so nothing in the check comes from the elimination.
    line::Matrix<double> W(4, 4, 0.0);
    W(0, 1) = 2.0;
    W(1, 2) = 1e9;
    W(2, 3) = 3.0;
    W(3, 0) = 4.0;
    W = as_generator(W);

    const fluid::FluidImmediateMatrix<double> r = fluid::fluid_eliminate_immediate_matrix(W);
    REQUIRE(r.eliminated);
    CHECK(r.fallback.empty());
    REQUIRE(r.state_map.size() == 3);
    CHECK(r.state_map[0] == 0);
    CHECK(r.state_map[1] == 2);
    CHECK(r.state_map[2] == 3);

    // The complement of a generator is a generator: its rows sum to zero. A
    // sign error in the correction term breaks this before it breaks anything
    // a stationary solve would notice.
    for (std::size_t i = 0; i < r.W.rows(); ++i) {
        double s = 0;
        for (std::size_t j = 0; j < r.W.cols(); ++j) s += r.W(i, j);
        CHECK(s == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
    }

    const std::vector<double> pi_full = mc::ctmc_solve(W);
    const std::vector<double> pi_red = mc::ctmc_solve(r.W);
    double mass = 0;
    for (std::size_t k = 0; k < r.state_map.size(); ++k) mass += pi_full[r.state_map[k]];
    REQUIRE(mass > 0.0);
    for (std::size_t k = 0; k < r.state_map.size(); ++k)
        CHECK(pi_red[k] == doctest::Approx(pi_full[r.state_map[k]] / mass).epsilon(1e-6));

    // The mass that conditioning divided out is the immediate state's own, and
    // it is of order 1/rate: that is the entire difference between eliminating
    // and not eliminating.
    CHECK(pi_full[1] < 1e-8);
}

TEST_CASE("fluid immediate: the matrix elimination falls back rather than trivializing") {
    // No state is immediate: the generator comes back untouched. A rule that
    // compared each row against the model's own largest rate would eliminate
    // the fastest state of any model at all, which is the failure this guards.
    line::Matrix<double> slow(3, 3, 0.0);
    slow(0, 1) = 1.0;
    slow(1, 2) = 2.0;
    slow(2, 0) = 3.0;
    slow = as_generator(slow);
    const fluid::FluidImmediateMatrix<double> a = fluid::fluid_eliminate_immediate_matrix(slow);
    CHECK(!a.eliminated);
    CHECK(a.fallback.empty());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(a.W(i, j) == doctest::Approx(slow(i, j)));

    // Everything is immediate: there is no timed process left to watch, so the
    // reference keeps the original rather than returning a system with nothing
    // in it. The fallback must SAY so -- a caller that cannot tell elimination
    // failed would integrate a stiff drift on the non-stiff path.
    line::Matrix<double> fast(2, 2, 0.0);
    fast(0, 1) = 1e9;
    fast(1, 0) = 2e9;
    fast = as_generator(fast);
    const fluid::FluidImmediateMatrix<double> b = fluid::fluid_eliminate_immediate_matrix(fast);
    CHECK(!b.eliminated);
    CHECK(!b.fallback.empty());
    CHECK(b.W.rows() == 2);

    CHECK_THROWS_AS(fluid::fluid_eliminate_immediate_matrix(line::Matrix<double>(2, 3, 0.0)),
                    line::InputError);
}

TEST_CASE("fluid stiff: a stiff and a non-stiff integration of the same drift agree") {
    // A two-delay cycle is not stiff, so both integrators are inside their
    // comfortable range and the only thing being compared is whether they solve
    // the same equations. They share no code: one is a four-stage Rosenbrock
    // with a numeric Jacobian, the other is LSODA running Adams.
    const double N = 4.0;
    qn::Network<double> m = delay_pair(N);
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(m.get_struct());
    std::vector<double> x0(sys.layout.nstates, 0.0);
    x0[sys.layout.qidx[0][0]] = N;

    fluid::FluidStiffOptions so;
    so.rtol = 1e-10;
    so.atol = 1e-12;
    const std::vector<double> xs =
        fluid::fluid_ode_solve_stiff(fluid::fluid_drift(sys), 0.0, 20.0, x0, so).final_state();
    const std::vector<double> xl = settle(sys, x0, 20.0);
    REQUIRE(xs.size() == xl.size());
    for (std::size_t i = 0; i < xs.size(); ++i)
        CHECK(xs[i] == doctest::Approx(xl[i]).epsilon(1e-6));

    // And both land on the exact fixed point of the linear drift, which neither
    // of them was told about.
    const double d1 = 1.0, d2 = 1.0 / 3.0;
    CHECK(xs[sys.layout.qidx[0][0]] == doctest::Approx(N * d1 / (d1 + d2)).epsilon(1e-6));
    CHECK(xs[sys.layout.qidx[0][0]] + xs[sys.layout.qidx[1][0]] ==
          doctest::Approx(N).epsilon(1e-8));
}

TEST_CASE("fluid stiff: the stiff arm integrates what elimination could not remove") {
    // The same three-station model, integrated WITHOUT eliminating anything.
    // This is the fallback path, and the point of it is that an L-stable method
    // damps the 1e-8 mode instead of being pinned to it: the answer is the same
    // fixed point, reached at the integrator's tolerance rather than exactly.
    const double N = 6.0, imm = lang::GlobalConstants::Immediate;
    qn::Network<double> m = delay_cycle3(N, imm);
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(m.get_struct());
    std::vector<double> x0(sys.layout.nstates, 0.0);
    x0[sys.layout.qidx[0][0]] = N;

    // Loose tolerances on purpose: the assertion below is at 1e-3, and asking a
    // step controller for more accuracy than the check needs across a 1e-8 mode
    // buys nothing but steps.
    fluid::FluidStiffOptions so;
    so.rtol = 1e-6;
    so.atol = 1e-8;
    const std::vector<double> xs =
        fluid::fluid_ode_solve_stiff(fluid::fluid_drift(sys), 0.0, 20.0, x0, so).final_state();

    const double d1 = 1.0, dm = 1.0 / imm, d2 = 1.0 / 3.0, tot = d1 + dm + d2;
    CHECK(xs[sys.layout.qidx[0][0]] == doctest::Approx(N * d1 / tot).epsilon(1e-3));
    CHECK(xs[sys.layout.qidx[2][0]] == doctest::Approx(N * d2 / tot).epsilon(1e-3));
    // The immediate station holds the mass the elimination conditions away.
    CHECK(xs[sys.layout.qidx[1][0]] < 1e-6);
}

TEST_CASE("fluid stiff: NonNegative is kept on the accurate arm and cleared on the fast one") {
    // y' = -1 from y0 = 1 has the exact solution 1 - t, so at t = 2 the true
    // value is -1. The reference keeps MATLAB's NonNegative for its accurate
    // stiff solver and CLEARS it for the fast one, which cannot honour it; that
    // difference is the only one this port can reproduce, so it is the one
    // asserted. The analytic solution is what makes the two arms
    // distinguishable at all.
    const std::function<void(double, const double*, double*)> decay =
        [](double, const double*, double* dy) { dy[0] = -1.0; };
    const std::vector<double> y0(1, 1.0);

    fluid::FluidStiffOptions keep;
    keep.stiff = true;
    const std::vector<double> a =
        fluid::fluid_ode_solve_stiff(decay, 0.0, 2.0, y0, keep).final_state();
    CHECK(a[0] == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));

    fluid::FluidStiffOptions clear;
    clear.stiff = false;
    const std::vector<double> b =
        fluid::fluid_ode_solve_stiff(decay, 0.0, 2.0, y0, clear).final_state();
    CHECK(b[0] == doctest::Approx(-1.0).epsilon(1e-8));
}

TEST_CASE("fluid stiff: a named MATLAB integrator is refused rather than impersonated") {
    // `options.odesolvers.*StiffOdeSolver` are function handles in the
    // reference, so a caller can ask for any MATLAB solver. This port has one
    // stiff method, and serving ode15s or ode23s under its own name would
    // report a method that did not run.
    const std::function<void(double, const double*, double*)> decay =
        [](double, const double*, double* dy) { dy[0] = -1.0; };
    const std::vector<double> y0(1, 1.0);

    const char* named[] = {"ode15s", "ode23s", "ode23t", "radau", "sundials"};
    for (std::size_t i = 0; i < sizeof(named) / sizeof(named[0]); ++i) {
        fluid::FluidStiffOptions o;
        o.solver = named[i];
        CHECK_THROWS_AS(fluid::fluid_ode_solve_stiff(decay, 0.0, 1.0, y0, o),
                        line::UnsupportedError);
    }
    fluid::FluidStiffOptions ok;
    ok.solver = "rosenbrock4";
    CHECK_NOTHROW(fluid::fluid_ode_solve_stiff(decay, 0.0, 1.0, y0, ok));

    // A reversed or empty range is a caller error, not a backwards solve.
    CHECK_THROWS_AS(fluid::fluid_ode_solve_stiff(decay, 1.0, 1.0, y0), line::InputError);
    CHECK_THROWS_AS(fluid::fluid_ode_solve_stiff(decay, 0.0, 1.0, std::vector<double>()),
                    line::InputError);
}

TEST_CASE("fluid stiff: the fluid path stays double, and says so by name") {
    // The elimination itself is field arithmetic and would be exact over the
    // rationals, but everything downstream of it integrates, so the entry point
    // that knows what the model is made of refuses a non-double one rather than
    // letting it reach an integrator whose coefficients are irrational.
    qn::Network<line::Rational> mr("fluid_stiff_rational");
    const std::size_t d1 = mr.add_delay("D1");
    const std::size_t d2 = mr.add_delay("D2");
    const std::size_t c = mr.add_closed_class("C1", 2.0, d1);
    mr.set_service(d1, c, lang::Distrib<line::Rational>::exp_rate(
                              line::num_traits<line::Rational>::from_int(1)));
    mr.set_service(d2, c, lang::Distrib<line::Rational>::exp_rate(
                              line::num_traits<line::Rational>::from_int(3)));
    qn::RoutingMatrix<line::Rational> P;
    P.set(c, c, d1, d2, line::num_traits<line::Rational>::from_int(1));
    P.set(c, c, d2, d1, line::num_traits<line::Rational>::from_int(1));
    mr.link(P);

    // The drift handed over is the DOUBLE model's: the gate is a statement
    // about the arithmetic the model was built in, and it has to fire before
    // anything in the fluid path is asked to run at that arithmetic.
    qn::Network<double> md = delay_pair(2.0);
    const fluid::FluidOdeSystem sd = fluid::fluid_ode_system(md.get_struct());
    CHECK_THROWS_AS(fluid::fluid_eliminate_immediate(mr.get_struct(), sd),
                    line::UnsupportedError);

    // The same call on a double model reaches the elimination and finds nothing
    // to do, so the gate is what refused above and not the model's shape.
    CHECK_NOTHROW(fluid::fluid_eliminate_immediate(md.get_struct(), sd));
}
