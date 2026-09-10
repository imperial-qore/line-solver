/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * The min-normal closure solved as a differential-algebraic system: the port of
 * `solver_fluid_dae.m`, its dispatch through the runner, and the RODAS leg.
 *
 * THE ORACLE IS MATLAB `SolverFLD(model,'method','dae')` AT LINE 3.0.7, printed
 * to six significant digits, which is what `t_dae.m` and `t_fcr.m` report. The
 * finite capacity region cases carry a SECOND oracle: SolverLDES on the same
 * model at 200000 samples, which is what says the cap means the same thing to a
 * simulator as it does to the constraint. The LDES throughputs are quoted in the
 * comments rather than asserted -- they are a sample mean with its own error and
 * the fluid answer is an approximation, so agreeing to 5% is the claim, not
 * agreeing to a tolerance.
 *
 * WHAT IS ASSERTED ABOUT THE TRAJECTORY IS STRUCTURAL, and deliberately so.
 * MATLAB integrates the covariance ALONGSIDE the mean when the closable state is
 * small (nc <= dae_maxcov); this port always holds the variance at its
 * stationary value, which is MATLAB's own fallback above that cap but not what
 * it does on a model this size. A point-by-point trajectory comparison would
 * therefore fail for a documented reason and prove nothing. What can be asserted
 * without an oracle is stronger anyway:
 *
 *   - the path STARTS at the model's initial condition, not at the answer
 *   - population conservation holds AT EVERY POINT to rounding, which is the
 *     whole reason for the singular mass matrix
 *   - the path ENDS at the fixed point the Newton solve found independently,
 *     so the algebraic route and the integrated route agree
 *   - the path AGREES with the LSODA closing trajectory while t is small, where
 *     the two closures have not yet separated
 *
 * The first of those is a regression test with a history: the transient used to
 * start from `seed.xvec`, the seed's CONVERGED state, so every horizon reported
 * the fixed point and the trajectory was invisible.
 *
 * @see fluid_dae.h - the solver
 * @see test_rodas.cpp - that the integrator under it matches its Fortran
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/fluid_dae.h"
#include "line/solvers/fluid/fluid_runner.h"

using namespace line;
using D = lang::Distrib<double>;
using lang::DropStrategy;

namespace {

/** Delay(Z=1) -> Q1(PS, c=2, mu=1.5) -> Q2(FCFS, c=1, mu=0.9), the `t_dae.m` model. */
qn::Network<double> tandem(std::size_t N) {
    qn::Network<double> m("tandem");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::PS);
    m.set_number_of_servers(q1, 2);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::FCFS);
    m.set_number_of_servers(q2, 1);
    const std::size_t k = m.add_closed_class("C1", static_cast<double>(N), d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q1, k, D::exp_rate(1.5));
    m.set_service(q2, k, D::exp_rate(0.9));
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q1, 1.0);
    P.set(k, k, q1, q2, 1.0);
    P.set(k, k, q2, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Delay(Z=1) -> Q1(PS,c=1,mu=2) -> Q2(PS,c=1,mu=2), the `t_fcr.m` model, with
 * {Q1,Q2} held to `B` jobs by a waiting-queue region when B > 0.
 *
 * Q1 AND Q2 ARE THE SAME STATION TWICE, which is what makes it a useful FCR
 * case: their means must be equal, so any asymmetry in the answer is the
 * solver's and not the model's.
 */
qn::Network<double> fcr(std::size_t N, double B) {
    qn::Network<double> m("fcr");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::PS);
    m.set_number_of_servers(q1, 1);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::PS);
    m.set_number_of_servers(q2, 1);
    const std::size_t k = m.add_closed_class("C1", static_cast<double>(N), d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q1, k, D::exp_rate(2.0));
    m.set_service(q2, k, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q1, 1.0);
    P.set(k, k, q1, q2, 1.0);
    P.set(k, k, q2, d, 1.0);
    m.link(P);
    if (B > 0.0)
        m.add_region(std::vector<std::size_t>{q1, q2}, std::vector<double>{-1.0}, B,
                     std::vector<DropStrategy>{DropStrategy::WAITQ});
    return m;
}

/** Source -> Cache -> HitQueue / MissQueue -> Sink: a decomposition, not one drift. */
qn::Network<double> cache_model() {
    qn::Network<double> m("cache");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = lang::ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double> >{std::vector<double>(4, 0.25), {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = m.add_cache("Cache", ch);
    const std::size_t hq = m.add_queue("HitQueue", lang::SchedStrategy::PS);
    const std::size_t mq = m.add_queue("MissQueue", lang::SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t init = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t miss = m.add_open_class("MissClass");
    m.set_arrival(src, init, D::exp_rate(1.0));
    m.set_service(hq, hit, D::exp_rate(4.0));
    m.set_service(mq, miss, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(init, init, src, cache, 1.0);
    P.set(hit, hit, cache, hq, 1.0);
    P.set(hit, hit, hq, snk, 1.0);
    P.set(miss, miss, cache, mq, 1.0);
    P.set(miss, miss, mq, snk, 1.0);
    m.link(P);
    return m;
}

/** Delay -> Queue(DPS), whose closure state is a matrix block per station. */
qn::Network<double> dps_model() {
    qn::Network<double> m("dps");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::DPS);
    m.set_number_of_servers(q, 1);
    const std::size_t c1 = m.add_closed_class("C1", 3.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 3.0, d);
    m.set_service(d, c1, D::exp_rate(1.0));
    m.set_service(d, c2, D::exp_rate(1.0));
    m.set_service(q, c1, D::exp_rate(2.0));
    m.set_service(q, c2, D::exp_rate(1.0));
    m.set_sched_param(q, c1, 1.0);
    m.set_sched_param(q, c2, 3.0);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 1.0);
    P.set(c1, c1, q, d, 1.0);
    P.set(c2, c2, d, q, 1.0);
    P.set(c2, c2, q, d, 1.0);
    m.link(P);
    return m;
}

fluid::FluidOptions dae_opts() {
    fluid::FluidOptions o;
    o.method = "dae";
    o.tol = 1e-8;
    return o;
}

/** The message of whatever `f` throws, or "" if it returns. */
std::string refusal(const std::function<void()>& f) {
    try {
        f();
    } catch (const std::exception& e) {
        return e.what();
    }
    return std::string();
}

}  // namespace

TEST_CASE("fluid dae: the steady state is MATLAB's, and it is the min-normal closure's") {
    qn::Network<double> model = tandem(6);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_dae(sn, dae_opts());

    // MATLAB SolverFLD(model,'method','dae'): [0.899261 0.611942 4.4888].
    CHECK(s.QN(0, 0) == doctest::Approx(0.899261).epsilon(1e-4));
    CHECK(s.QN(1, 0) == doctest::Approx(0.611942).epsilon(1e-4));
    CHECK(s.QN(2, 0) == doctest::Approx(4.488800).epsilon(1e-4));

    // Same closure as `minnormal`, so the same answer: what differs is that the
    // mean and the variance are solved together rather than alternated.
    fluid::FluidOptions om = dae_opts();
    om.method = "minnormal";
    const fluid::FluidSolution mn = fluid::solver_fluid_moments(sn, om);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(s.QN(i, 0) == doctest::Approx(mn.QN(i, 0)).epsilon(1e-4));

    // The second moment is reported, which is what `-a var` needs.
    CHECK(s.has_moments);
    CHECK(s.moments.QVar.rows() == 3);
    for (std::size_t i = 0; i < 3; ++i) CHECK(s.moments.QVar(i, 0) >= 0.0);
}

TEST_CASE("fluid dae: conservation is an equation, not a consequence of the drift") {
    for (std::size_t N : {4, 6, 12}) {
        qn::Network<double> model = tandem(N);
        const qn::NetworkStruct<double>& sn = model.get_struct();
        const fluid::FluidSolution s = fluid::solver_fluid_dae(sn, dae_opts());
        double sum = 0.0;
        for (std::size_t i = 0; i < s.QN.rows(); ++i) sum += s.QN(i, 0);
        // Held to the Newton tolerance because it IS one of the equations, not
        // to integrator tolerance because it happened to be preserved.
        CHECK(sum == doctest::Approx(static_cast<double>(N)).epsilon(1e-9));
    }
}

TEST_CASE("fluid dae: reachable by options.method, through the runner") {
    qn::Network<double> model = tandem(6);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const fluid::FluidSolution direct = fluid::solver_fluid_dae(sn, dae_opts());

    for (const char* spelling : {"dae", "fluid.dae"}) {
        fluid::FluidOptions o = dae_opts();
        o.method = spelling;
        const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, o);
        CHECK(s.method == "dae");
        CHECK(s.has_moments);
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(s.QN(i, 0) == doctest::Approx(direct.QN(i, 0)).epsilon(1e-12));
    }

    // `solver_fluid` is the port of solver_fluid_analyzer.m alone and cannot
    // call the DAE without a cyclic include, so it names where the method lives
    // instead of reporting it as unknown.
    const std::string e = refusal([&] { fluid::solver_fluid(sn, dae_opts()); });
    CHECK(e.find("solver_fluid_run_analyzer") != std::string::npos);
}

TEST_CASE("fluid dae: the FCFS non-exponential refit runs for it, as it does for minnormal") {
    // `fluid_method_refits_fcfs` now names `dae`, matching the reference's second
    // switch (`solver_fluid_analyzer.m:125`). It refits because it IS the
    // min-normal closure -- same drift, same rate factors -- so the phase count
    // its answer depends on has to be refitted on the same schedule. An Erlang
    // FCFS station is what makes the loop actually iterate: with exponential
    // service everywhere the loop is entered and immediately satisfied.
    qn::Network<double> m("erl");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    m.set_number_of_servers(q1, 1);
    const std::size_t q2 = m.add_queue("Q2", lang::SchedStrategy::PS);
    m.set_number_of_servers(q2, 1);
    const std::size_t k = m.add_closed_class("C1", 6.0, d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q1, k, D::erlang(3.0, 3));  // mean 1, SCV 1/3
    m.set_service(q2, k, D::exp_rate(1.2));
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q1, 1.0);
    P.set(k, k, q1, q2, 1.0);
    P.set(k, k, q2, d, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, dae_opts());
    fluid::FluidOptions omn = dae_opts();
    omn.method = "minnormal";
    const fluid::FluidSolution mn = fluid::solver_fluid_run_analyzer(sn, omn);

    // The same number of sweeps as `minnormal`, which is the claim: one schedule,
    // not two.
    CHECK(s.refit_sweeps > 0);
    CHECK(s.refit_sweeps == mn.refit_sweeps);
    // ... and the refitted answer still agrees with the closure it shares.
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(s.QN(i, 0) == doctest::Approx(mn.QN(i, 0)).epsilon(1e-3));
    double sum = 0.0;
    for (std::size_t i = 0; i < s.QN.rows(); ++i) sum += s.QN(i, 0);
    CHECK(sum == doctest::Approx(6.0).epsilon(1e-9));
}

TEST_CASE("fluid dae: 'default' resolves to it only where nothing else can run") {
    // The DAE route is a request everywhere the model has an alternative:
    // `fluid_resolve_default_method.m` chooses among rmf, minnormal, closing and
    // matrix, and a method that WIDENS the feature set must not be reached on a
    // model that does not need the widening. The one exception is a binding
    // buffer or a capacity region, where every other method is refused outright
    // -- see the capacity and region cases below.
    for (std::size_t N : {4, 8}) {
        qn::Network<double> model = tandem(N);
        const qn::NetworkStruct<double>& sn = model.get_struct();
        fluid::FluidOptions o = dae_opts();
        o.method = "default";
        CHECK(fluid::detail::fluid_resolve_method(sn, "default", o) != "dae");
        CHECK(fluid::solver_fluid_run_analyzer(sn, o).method != "dae");
    }
}

TEST_CASE("fluid dae: a finite capacity region binds, and only for this method") {
    // MATLAB SolverFLD(model,'method','dae') on the same five configurations,
    // with SolverLDES at 200000 samples beside it:
    //   N=8  B=2   X 1.26677  (ldes 1.33151)
    //   N=8  B=4   X 1.60178  (ldes 1.59691)
    //   N=16 B=4   X 1.60178  (ldes 1.59809)
    //   N=16 B=8   X 1.79953  (ldes 1.77546)
    //   N=30 B=6   X 1.73272  (ldes 1.71256)
    // LDES also reports the region population itself at the cap to within its
    // own noise (1.99994, 3.97705, 4.00000, 7.99996, 6.00000), which is what
    // says the constraint means what the simulator means.
    struct Case {
        std::size_t N;
        double B, X;
    };
    const Case cases[] = {{8, 2.0, 1.26677},  {8, 4.0, 1.60178}, {16, 4.0, 1.60178},
                          {16, 8.0, 1.79953}, {30, 6.0, 1.73272}};
    for (const Case& c : cases) {
        CAPTURE(c.N);
        CAPTURE(c.B);
        qn::Network<double> model = fcr(c.N, c.B);
        const qn::NetworkStruct<double>& sn = model.get_struct();
        const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, dae_opts());

        // The cap is an equation once it binds, so it holds exactly and not to
        // within a penalty.
        const double region = s.QN(1, 0) + s.QN(2, 0);
        CHECK(region == doctest::Approx(c.B).epsilon(1e-6));
        CHECK(s.TN(1, 0) == doctest::Approx(c.X).epsilon(1e-4));

        // Little's law at the delay: Q = X * Z with Z = 1. The blocked mass is
        // held in the staging coordinate and NOT counted at the delay, which is
        // what this checks -- throttling admission by slowing the upstream
        // station instead would break exactly this identity.
        CHECK(s.QN(0, 0) / s.TN(0, 0) == doctest::Approx(1.0).epsilon(1e-6));

        // Q1 and Q2 are the same station twice, so the answer must be symmetric.
        CHECK(s.QN(1, 0) == doctest::Approx(s.QN(2, 0)).epsilon(1e-6));
    }

    // Every other fluid method still refuses a region outright, which for them
    // is still true: an ODE integrated through a cap it cannot see returns the
    // UNCONSTRAINED population, with no warning.
    for (const char* meth : {"closing", "matrix", "minnormal", "statedep"}) {
        CAPTURE(meth);
        fluid::FluidOptions o = dae_opts();
        o.method = meth;
        qn::Network<double> model = fcr(8, 2.0);
        const qn::NetworkStruct<double>& sn = model.get_struct();
        const std::string e = refusal([&] { fluid::solver_fluid_run_analyzer(sn, o); });
        CHECK(e.find("Finite Capacity Region") != std::string::npos);
    }

    // `default` is NOT among them: it stands for the one method that carries the
    // region, so the resolution sends the model here rather than to a refusal
    // whose advice was to type `dae`.
    {
        qn::Network<double> model = fcr(8, 2.0);
        const qn::NetworkStruct<double>& sn = model.get_struct();
        CHECK(fluid::detail::fluid_resolve_method(sn, "default", dae_opts()) == "dae");
        fluid::FluidOptions o = dae_opts();
        o.method = "default";
        CHECK(fluid::solver_fluid_run_analyzer(sn, o).method == "dae");
    }

    // With no region the constraint machinery is inert, not merely inactive.
    qn::Network<double> plain_model = fcr(8, 0.0);
    const qn::NetworkStruct<double>& plain = plain_model.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(plain, dae_opts());
    double sum = 0.0;
    for (std::size_t i = 0; i < s.QN.rows(); ++i) sum += s.QN(i, 0);
    CHECK(sum == doctest::Approx(8.0).epsilon(1e-9));
    CHECK(s.QN(1, 0) == doctest::Approx(s.QN(2, 0)).epsilon(1e-6));
}

TEST_CASE("fluid dae: the refusals name the model feature") {
    // A cache model is a decomposition, so there is no single drift to attach
    // the constraint to; `minnormal` carries the closure inside the cache
    // analyzer's network step and still answers it.
    qn::Network<double> cache_net = cache_model();
    const qn::NetworkStruct<double>& cm = cache_net.get_struct();
    const std::string ec = refusal([&] { fluid::solver_fluid_run_analyzer(cm, dae_opts()); });
    CHECK(ec.find("caching") != std::string::npos);
    fluid::FluidOptions omn = dae_opts();
    omn.method = "minnormal";
    CHECK(refusal([&] { fluid::solver_fluid_run_analyzer(cm, omn); }).empty());

    // DPS closes on the covariance BETWEEN a station's class coordinates, where
    // this carries one scalar variance per station. The per-method feature set
    // refuses it before the solve, so the message names the discipline.
    qn::Network<double> dps_net = dps_model();
    const qn::NetworkStruct<double>& dm = dps_net.get_struct();
    const std::string ed = refusal([&] { fluid::solver_fluid_run_analyzer(dm, dae_opts()); });
    CHECK(ed.find("DPS") != std::string::npos);
    // ... and the direct call refuses it too, since `solver_fluid_dae` is
    // reachable without the gate.
    CHECK(refusal([&] { fluid::solver_fluid_dae(dm, dae_opts()); }).find("DPS") !=
          std::string::npos);
}

TEST_CASE("fluid dae: the RODAS transient starts at the model and conserves along it") {
    const std::size_t N = 6;
    qn::Network<double> model = tandem(N);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const fluid::FluidSolution ss = fluid::solver_fluid_dae(sn, dae_opts());

    const std::vector<fluid::FluidTranPoint> tr =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 40.0, 41);
    REQUIRE(tr.size() == 41);

    // t = 0 is the model's own initial condition: `initDefault` puts the whole
    // closed population at the reference station. THIS IS THE REGRESSION: the
    // integration used to start from the seed's CONVERGED state, so every
    // horizon reported the fixed point and no trajectory existed at all.
    CHECK(tr.front().t == doctest::Approx(0.0));
    CHECK(tr.front().QN(0, 0) == doctest::Approx(static_cast<double>(N)).epsilon(1e-12));
    CHECK(tr.front().QN(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(tr.front().QN(2, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // The singular mass matrix carries `C x = N` as an equation, so the
    // population is conserved AT EVERY reported point rather than drifting with
    // the integrator. This is the property the whole DAE form exists for.
    for (std::size_t j = 0; j < tr.size(); ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < tr[j].QN.rows(); ++i) sum += tr[j].QN(i, 0);
        CAPTURE(tr[j].t);
        CHECK(std::fabs(sum - static_cast<double>(N)) < 1e-9);
    }

    // The grid is the caller's, monotone, and ends at the horizon.
    for (std::size_t j = 1; j < tr.size(); ++j) CHECK(tr[j].t > tr[j - 1].t);
    CHECK(tr.back().t == doctest::Approx(40.0));

    // THE TWO ROUTES MEET. The last point is reached by integrating the drift;
    // the table is reached by solving the same equations algebraically with
    // Newton. Neither is used to compute the other, so agreeing here says both
    // are solving the same system.
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(tr.back().QN(i, 0) == doctest::Approx(ss.QN(i, 0)).epsilon(1e-6));

    // The trajectory actually moves: a path that started at the answer would
    // pass every check above.
    double travel = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        travel = std::max(travel, std::fabs(tr.front().QN(i, 0) - tr.back().QN(i, 0)));
    CHECK(travel > 1.0);

    // A caller's own grid is honoured, including one that skips t=0.
    const std::vector<double> own{0.5, 1.0, 2.0, 40.0};
    const std::vector<fluid::FluidTranPoint> pick =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 40.0, 101, own);
    REQUIRE(pick.size() == own.size());
    for (std::size_t j = 0; j < own.size(); ++j) CHECK(pick[j].t == doctest::Approx(own[j]));
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(pick.back().QN(i, 0) == doctest::Approx(ss.QN(i, 0)).epsilon(1e-6));
}

TEST_CASE("fluid dae: early in the transient the two closures have not yet separated") {
    // There is no MATLAB oracle for this trajectory, so the check is
    // against the OTHER trajectory the tree can produce: the first-order closing
    // drift under LSODA. The two start from the same state and are driven by the
    // same rates until the variance correction has had time to matter, so a
    // disagreement at small t would mean the drift itself was wired wrong --
    // while a disagreement at large t is the closure doing its job.
    qn::Network<double> model = tandem(6);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    // THE GRID MUST START AT ZERO for `solver_fluid_transient`: it hands the
    // vector to LSODA as a tspan, so the first entry is read as the INTEGRATION
    // START and the trajectory comes back shifted if it is not t0. That is the
    // LSODA path's convention and not this method's -- `solver_fluid_dae_transient`
    // integrates from t=0 whatever the grid asks for -- but a comparison between
    // the two has to be made on ground both agree on.
    const std::vector<double> grid{0.0, 0.002, 0.02, 0.1};

    fluid::FluidOptions oc = dae_opts();
    oc.method = "closing";
    const std::vector<fluid::FluidTranPoint> cls =
        fluid::solver_fluid_transient(sn, oc, 0.1, 101, grid);
    const std::vector<fluid::FluidTranPoint> dae =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 0.1, 101, grid);
    REQUIRE(cls.size() == grid.size());
    REQUIRE(dae.size() == grid.size());
    for (std::size_t j = 0; j < grid.size(); ++j) {
        CAPTURE(grid[j]);
        CHECK(dae[j].t == doctest::Approx(cls[j].t));
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(std::fabs(dae[j].QN(i, 0) - cls[j].QN(i, 0)) < 5e-3);
    }
    // ... and by the time the horizon is long the two HAVE separated, which is
    // the correction the closure exists to make. A test that only showed
    // agreement would be satisfied by two copies of the same drift.
    const std::vector<fluid::FluidTranPoint> late_cls =
        fluid::solver_fluid_transient(sn, oc, 2.0, 3);
    const std::vector<fluid::FluidTranPoint> late_dae =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 2.0, 3);
    double sep = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        sep = std::max(sep, std::fabs(late_dae.back().QN(i, 0) - late_cls.back().QN(i, 0)));
    CHECK(sep > 0.01);
}

TEST_CASE("fluid dae: the transient covariance is integrated, not held") {
    // THE PROPERTY THAT SEPARATES THIS FROM `minnormal`, which evaluates its
    // whole transient at the single STATIONARY variance. Here the nc x nc
    // covariance block joins the state as ordinary differential rows -- only the
    // one row per closed chain is algebraic -- so the drift is read at the
    // variance the trajectory actually HAS at each instant.
    const std::size_t N = 6;
    qn::Network<double> model = tandem(N);
    const qn::NetworkStruct<double>& sn = model.get_struct();

    const std::vector<fluid::FluidTranPoint> tr =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 40.0, 41);
    REQUIRE(tr.size() == 41);

    // Sigma(0) = 0 is the consistent AND the physically right initialisation:
    // the population at t=0 is a known deterministic state, so it has no
    // variance.
    REQUIRE(tr.front().QVar.rows() > 0);
    for (std::size_t i = 0; i < tr.front().QVar.rows(); ++i)
        CHECK(tr.front().QVar(i, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // and it must actually move, or the rows were carried and never integrated
    double grew = 0.0;
    for (std::size_t i = 0; i < tr.back().QVar.rows(); ++i)
        grew = std::max(grew, tr.back().QVar(i, 0));
    CHECK(grew > 1e-3);

    // A variance is a variance: never negative, at any reported point.
    for (std::size_t j = 0; j < tr.size(); ++j)
        for (std::size_t i = 0; i < tr[j].QVar.rows(); ++i) {
            CAPTURE(tr[j].t);
            CHECK(tr[j].QVar(i, 0) >= 0.0);
        }

    // The population is STILL conserved with the covariance rows present: they
    // are differential, so the mass matrix leaves them at 1 and only the chain
    // row stays algebraic. Getting that backwards would break this.
    for (std::size_t j = 0; j < tr.size(); ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < tr[j].QN.rows(); ++i) sum += tr[j].QN(i, 0);
        CAPTURE(tr[j].t);
        CHECK(std::fabs(sum - static_cast<double>(N)) < 1e-9);
    }

    // A long horizon must reach the STATIONARY covariance the algebraic solve
    // converged to, which is the only external check available: the trajectory
    // and the table are then demonstrably the same closure.
    const fluid::FluidSolution ss = fluid::solver_fluid_dae(sn, dae_opts());
    for (std::size_t i = 0; i < tr.back().QVar.rows(); ++i)
        CHECK(tr.back().QVar(i, 0) == doctest::Approx(ss.moments.QVar(i, 0)).epsilon(1e-3));
}

TEST_CASE("fluid dae: above dae_maxcov the variance is held rather than integrated") {
    // The cap is a cost decision, not a correctness one: nc^2 extra states
    // differenced numerically is nc^4 work. Below it the covariance is carried;
    // above it the mean is STILL a DAE -- conservation stays an equation -- and
    // the variance falls back to the stationary one, which is what `minnormal`
    // uses for the whole of its transient anyway.
    qn::Network<double> model = tandem(6);
    const qn::NetworkStruct<double>& sn = model.get_struct();

    fluid::FluidDaeOptions held;
    held.maxcov = 0;   // nothing is small enough, so nothing is integrated
    const std::vector<fluid::FluidTranPoint> tr =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 40.0, 41,
                                          std::vector<double>(), held);
    REQUIRE(tr.size() == 41);
    // no per-point second moment is reported, rather than a repeated constant
    CHECK(tr.front().QVar.rows() == 0);
    // and the mean is unaffected in what it must satisfy
    for (std::size_t j = 0; j < tr.size(); ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < tr[j].QN.rows(); ++i) sum += tr[j].QN(i, 0);
        CHECK(std::fabs(sum - 6.0) < 1e-9);
    }
}

TEST_CASE("fluid dae: '-a tran' honours the method instead of substituting closing") {
    // `solver_fluid_tran_avg` integrates the FIRST-ORDER closing drift whatever
    // the method says, which is the reference's rule for every steady-state
    // device. `dae` is the exception the reference itself makes
    // (`getTranAvg.m:155-166`), so the router has to send it to its own
    // trajectory -- otherwise a caller asking for `-a tran --method dae` is
    // handed a different method's answer under this method's name.
    qn::Network<double> model = tandem(6);
    const qn::NetworkStruct<double>& sn = model.get_struct();

    fluid::FluidOptions od = dae_opts();
    od.timespan_end = 30.0;
    const std::vector<fluid::FluidTranPoint> viadae =
        fluid::solver_fluid_run_transient(sn, od, 31);
    const std::vector<fluid::FluidTranPoint> direct =
        fluid::solver_fluid_dae_transient(sn, od, 30.0, 31);
    REQUIRE(viadae.size() == direct.size());
    for (std::size_t j = 0; j < viadae.size(); ++j)
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(viadae[j].QN(i, 0) == doctest::Approx(direct[j].QN(i, 0)).epsilon(1e-12));

    // Every other method keeps the path it had, byte for byte.
    fluid::FluidOptions oc = od;
    oc.method = "closing";
    const std::vector<fluid::FluidTranPoint> viacls =
        fluid::solver_fluid_run_transient(sn, oc, 31);
    const std::vector<fluid::FluidTranPoint> oldcls = fluid::solver_fluid_tran_avg(sn, oc, 31);
    REQUIRE(viacls.size() == oldcls.size());
    for (std::size_t j = 0; j < viacls.size(); ++j)
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(viacls[j].QN(i, 0) == doctest::Approx(oldcls[j].QN(i, 0)).epsilon(1e-12));

    // With no horizon on the options both arms fall back to the same rule --
    // thirty mean events of the slowest rate in `sn.rates` -- so the two
    // trajectories are read at the same times and only the drift differs.
    fluid::FluidOptions on = dae_opts();
    const std::vector<fluid::FluidTranPoint> nod = fluid::solver_fluid_run_transient(sn, on, 11);
    fluid::FluidOptions onc = on;
    onc.method = "closing";
    const std::vector<fluid::FluidTranPoint> noc = fluid::solver_fluid_run_transient(sn, onc, 11);
    REQUIRE(nod.size() == 11);
    REQUIRE(noc.size() == 11);
    for (std::size_t j = 0; j < nod.size(); ++j) CHECK(nod[j].t == doctest::Approx(noc[j].t));
}

TEST_CASE("fluid dae: a fork model reaches the DAE through the MMT fixed point") {
    // THE FEATURE SET LEAVES Fork DECLARED for `dae`, as the reference does, so
    // `solver_fluid_run_analyzer` routes a fork model into `fluid_fork_join_run` and the
    // inner solve is the DAE on the MIXED network the transform emits -- open
    // auxiliary classes over a closed model. What is asserted here is only that
    // the path is REACHABLE and self-consistent, because there is no oracle for
    // it: the transform's state is not the model's population, so the value
    // cannot be checked against MATLAB the way the others are. It is worth
    // pinning even so, since the alternative to a weak test on this path is no
    // test at all, and the failure mode of an undeclared-but-reached method is a
    // number rather than a refusal.
    qn::Network<double> m("fj");
    const std::size_t d = m.add_delay("Delay1");
    const std::size_t q1 = m.add_queue("Queue1", lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", lang::SchedStrategy::PS);
    m.set_number_of_servers(q1, 1);
    m.set_number_of_servers(q2, 1);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("class1", 5.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);

    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(m.get_struct(), dae_opts());
    CHECK(s.method == "dae");
    for (std::size_t i = 0; i < s.QN.rows(); ++i) {
        CAPTURE(i);
        CHECK(std::isfinite(s.QN(i, 0)));
        CHECK(s.QN(i, 0) >= 0.0);
    }
    // The two branches are the same queue twice, so the answer must be symmetric
    // whatever its level.
    REQUIRE(s.QN.rows() >= 3);
    CHECK(s.QN(1, 0) == doctest::Approx(s.QN(2, 0)).epsilon(1e-6));
}

TEST_CASE("fluid dae: the FCR transient holds the cap and ends where the table says") {
    // The transient under a cap is a HYBRID DAE: the system switches every time the
    // region fills or drains, so the horizon is covered by segments, each ending at
    // a LOCATED crossing. Integrating with the binding set frozen would report the
    // unconstrained path through a cap the model declares -- so what is tested is
    // that the path never exceeds the cap and that it ends where the steady-state
    // solve says it should.
    qn::Network<double> model = fcr(8, 2.0);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const std::vector<fluid::FluidTranPoint> path =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 30.0, 31);
    REQUIRE(path.size() >= 2);
    for (std::size_t s = 0; s < path.size(); ++s) {
        const double region = path[s].QN(1, 0) + path[s].QN(2, 0);
        CHECK(region <= 2.0 + 1e-6);
    }
    const fluid::FluidSolution ss = fluid::solver_fluid_run_analyzer(sn, dae_opts());
    CHECK(path.back().QN(1, 0) == doctest::Approx(ss.QN(1, 0)).epsilon(1e-3));
    CHECK(path.back().QN(2, 0) == doctest::Approx(ss.QN(2, 0)).epsilon(1e-3));

    // The same path reaches a caller who asks for a horizon on the options instead
    // of calling the transient by name: the state at that horizon, on the cap.
    fluid::FluidOptions o = dae_opts();
    o.timespan_end = 30.0;
    const fluid::FluidSolution at = fluid::solver_fluid_run_analyzer(sn, o);
    CHECK(at.QN(1, 0) + at.QN(2, 0) <= 2.0 + 1e-6);

    // With no region that same horizon is integrated, as it always was.
    qn::Network<double> plain_model = fcr(8, 0.0);
    const qn::NetworkStruct<double>& plain = plain_model.get_struct();
    CHECK(refusal([&] { fluid::solver_fluid_run_analyzer(plain, o); }).empty());
}

// --------------------------------------------------------------- station buffers

namespace {

/** Delay -> Queue(cap), one closed class: the blocked job is HELD upstream. */
qn::Network<double> capped(std::size_t N, double cap, double mu) {
    qn::Network<double> m("capped");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    m.set_number_of_servers(q, 1);
    if (cap > 0.0) m.set_capacity(q, cap);
    const std::size_t k = m.add_closed_class("C1", static_cast<double>(N), d);
    m.set_service(d, k, D::exp_rate(1.0));
    m.set_service(q, k, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(k, k, d, q, 1.0);
    P.set(k, k, q, d, 1.0);
    m.link(P);
    return m;
}

/** Source -> Queue(cap) -> Sink, one open class: the arrival is LOST. */
qn::Network<double> loss_model(double lambda, double mu, double cap) {
    qn::Network<double> m("loss");
    const std::size_t s = m.add_source("Source");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t k0 = m.add_sink("Sink");
    m.set_number_of_servers(q, 1);
    m.set_capacity(q, cap);
    const std::size_t k = m.add_open_class("C1");
    m.set_arrival(s, k, D::exp_rate(lambda));
    m.set_service(q, k, D::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(k, k, s, q, 1.0);
    P.set(k, k, q, k0, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fluid dae: a station buffer binds and holds the blocked job upstream") {
    // A station buffer is the one-station case of the same row -- but NOT of the same
    // model. LINE refuses to lose a closed job (State.arrivalIsLost) and disables the
    // upstream departure instead, so the blocked mass is still AT the upstream
    // station and still counted there: the station queues sum to the whole
    // population and nothing is staged. That is the opposite of a region, whose
    // blocked jobs are reported separately.
    qn::Network<double> model = capped(8, 2.0, 2.0);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, dae_opts());
    CHECK(s.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(8.0).epsilon(1e-6));

    // A cap the population cannot reach is pruned, and the answer is the
    // unconstrained one: refresh_capacity derives a finite classcap for every closed
    // model, so that is the common case rather than a corner one.
    qn::Network<double> wide = capped(8, 50.0, 2.0);
    qn::Network<double> free_model = capped(8, 0.0, 2.0);
    const fluid::FluidSolution sw = fluid::solver_fluid_run_analyzer(wide.get_struct(), dae_opts());
    const fluid::FluidSolution sf = fluid::solver_fluid_run_analyzer(free_model.get_struct(), dae_opts());
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(sw.QN(i, 0) == doctest::Approx(sf.QN(i, 0)).epsilon(1e-6));
}

TEST_CASE("fluid dae: an arrival is lost at a full buffer of an open class") {
    // The same predicate that HOLDS a closed job LOSES an open one: the external
    // stream is memoryless, so a job that finds the buffer full never enters. In
    // overload the fluid loss rate is exact -- what gets through is the server -- and
    // the answer sits on the cap.
    qn::Network<double> model = loss_model(4.0, 2.0, 5.0);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(sn, dae_opts());
    CHECK(s.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-6));
    CHECK(s.TN(1, 0) == doctest::Approx(2.0).epsilon(1e-4));

    // Below saturation the same buffer is inert: the fluid never reaches it, so the
    // answer is the uncapped one and the loss it reports is zero.
    qn::Network<double> light = loss_model(1.0, 2.0, 5.0);
    const fluid::FluidSolution sl = fluid::solver_fluid_run_analyzer(light.get_struct(), dae_opts());
    CHECK(sl.TN(1, 0) == doctest::Approx(1.0).epsilon(1e-4));
    CHECK(sl.QN(1, 0) < 5.0);
}

TEST_CASE("fluid dae: a station buffer transient never leaves the feasible set") {
    qn::Network<double> model = capped(8, 2.0, 2.0);
    const qn::NetworkStruct<double>& sn = model.get_struct();
    const std::vector<fluid::FluidTranPoint> path =
        fluid::solver_fluid_dae_transient(sn, dae_opts(), 30.0, 31);
    REQUIRE(path.size() >= 2);
    for (std::size_t s = 0; s < path.size(); ++s) CHECK(path[s].QN(1, 0) <= 2.0 + 1e-6);
    const fluid::FluidSolution ss = fluid::solver_fluid_run_analyzer(sn, dae_opts());
    CHECK(path.back().QN(1, 0) == doctest::Approx(ss.QN(1, 0)).epsilon(1e-3));
}

// ---------------------------------------------------------------------------
// The non-hyperbolic fallback ladder: minnormal -> dae -> first order
// ---------------------------------------------------------------------------
//
// `fluid_lyapunov` raises when the drift Jacobian at the converged mean is not
// exponentially stable on range(D). The dominant case is NEUTRAL rather than
// unstable, and it is an artifact of the alternation: `solver_fluid_moments`
// must start at sigma2 = 0, where min(n,c) has no derivative, so a saturated or
// balanced model's first-order fixed point lands on the kink and sits on a
// continuum of equilibria. `solver_fluid_dae` seeds the variance positive and
// never adopts sigma2 = 0, so the same closure has an isolated, hyperbolic fixed
// point there.

/**
 * Two identical stations in a closed cycle: the fluid drift is DEGENERATE. Every
 * state with both populations at or above the server count is a fluid
 * equilibrium, so a first-order method has no reason to prefer one point of that
 * continuum over another and lands wherever the integrator stopped.
 */
qn::Network<double> balanced_cycle(double N, std::size_t nservers,
                                   lang::SchedStrategy sched) {
    qn::Network<double> m("balanced");
    const std::size_t q1 = m.add_queue("Q1", sched);
    const std::size_t q2 = m.add_queue("Q2", sched);
    m.set_number_of_servers(q1, nservers);
    m.set_number_of_servers(q2, nservers);
    const std::size_t c = m.add_closed_class("C", N, q1);
    m.set_service(q1, c, D::exp_rate(1.0));
    m.set_service(q2, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, q1, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("fluid ladder: a declined minnormal lands on dae and on the exact answer") {
    // By symmetry the exact answer splits the population evenly, which is what the
    // closure gives once the degeneracy is broken. The first-order fallback this
    // replaces returned [9 1] at N=10 -- a point of the continuum, not the mean.
    const double pop[] = {4.0, 6.0, 6.0, 10.0, 10.0, 6.0, 6.0};
    const std::size_t servers[] = {1, 1, 2, 1, 2, 1, 2};
    const lang::SchedStrategy sched[] = {
        lang::SchedStrategy::PS,   lang::SchedStrategy::PS,   lang::SchedStrategy::PS,
        lang::SchedStrategy::PS,   lang::SchedStrategy::PS,   lang::SchedStrategy::FCFS,
        lang::SchedStrategy::FCFS};
    for (std::size_t k = 0; k < 7; ++k) {
        qn::Network<double> m = balanced_cycle(pop[k], servers[k], sched[k]);
        fluid::FluidOptions o = dae_opts();
        o.method = "minnormal";
        const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(m.get_struct(), o);
        CHECK(s.QN(0, 0) == doctest::Approx(pop[k] / 2.0).epsilon(1e-6));
        CHECK(s.QN(1, 0) == doctest::Approx(pop[k] / 2.0).epsilon(1e-6));
    }
}

TEST_CASE("fluid ladder: a genuinely unstable fixed point walks past dae to first order") {
    // An overloaded open station has no stationary distribution at all, so no
    // closure has a stationary covariance there and `dae` declines on the same
    // error `minnormal` did. The mean is still reported, by the first-order
    // method, which is the whole reason the last rung exists.
    qn::Network<double> m("overloaded");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Q1", lang::SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C");
    m.set_arrival(src, c, D::exp_rate(1.1));
    m.set_service(q, c, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(c, c, src, q, 1.0);
    P.set(c, c, q, snk, 1.0);
    m.link(P);

    fluid::FluidOptions o = dae_opts();
    o.method = "minnormal";
    // reaching a finite answer at all IS the assertion: both closure rungs
    // decline, and without the last rung this throws
    const fluid::FluidSolution s = fluid::solver_fluid_run_analyzer(m.get_struct(), o);
    CHECK(std::isfinite(s.QN(1, 0)));
}

TEST_CASE("fluid ladder: the dae rung is declined in advance for what it cannot take") {
    // `fluid_dae_applicable` is the STATIC difference set between the closures. A
    // rung entered only to be refused a moment later would spend a whole seed
    // integration to learn what the model already declares.
    std::string why;

    // DPS closes on the covariance BETWEEN class coordinates, not the station
    // total, so its closure state is a matrix block rather than the scalar the
    // Newton vector carries.
    qn::Network<double> dm = dps_model();
    CHECK(!fluid::detail::fluid_dae_applicable(dm.get_struct(), dae_opts(), why));
    CHECK(why.find("covariance between") != std::string::npos);

    // A cache model is a decomposition, and `dae` has no arm for it.
    qn::Network<double> cm = cache_model();
    CHECK(!fluid::detail::fluid_dae_applicable(cm.get_struct(), dae_opts(), why));
    CHECK(why.find("decomposition") != std::string::npos);

    // The simultaneous solve is quartic where one Lyapunov solve is cubic, so it
    // carries its own cap, lower than moment_maxstate. It is read off the
    // OPTIONS, as `options.config.dae_maxstate` is in the other three codebases.
    qn::Network<double> bc = balanced_cycle(6.0, 1, lang::SchedStrategy::PS);
    fluid::FluidOptions capped = dae_opts();
    capped.dae_maxstate = 1;
    CHECK(!fluid::detail::fluid_dae_applicable(bc.get_struct(), capped, why));
    CHECK(why.find("dae_maxstate") != std::string::npos);
    // and the same cap pinned on the struct a caller may hand in directly
    fluid::FluidDaeOptions capped_opt;
    capped_opt.maxstate = 1;
    CHECK(!fluid::detail::fluid_dae_applicable(bc.get_struct(), dae_opts(), why, capped_opt));
    CHECK(why.find("dae_maxstate") != std::string::npos);

    // and the model the ladder does take
    CHECK(fluid::detail::fluid_dae_applicable(bc.get_struct(), dae_opts(), why));
    CHECK(why.empty());
}

TEST_CASE("fluid dae: the caps are reachable from the options, and the Newton cap follows iter_max") {
    // `dae_maxstate` and `dae_maxcov` are `options.config` entries in MATLAB,
    // the JAR and native Python; here they arrive on FluidOptions and have to
    // reach the struct the route reads. Zero means NOT SET, so a caller pinning
    // the struct still wins where the options are silent -- and the Newton cap
    // is max(50, iter_max) in all three references, not a flat 50.
    fluid::FluidOptions o;
    fluid::FluidDaeOptions d;
    o.iter_max = 200;
    fluid::FluidDaeOptions m = fluid::fluid_dae_options(o, d);
    CHECK(m.maxstate == 100);
    CHECK(m.maxcov == 25);
    CHECK(m.newton_max == 200);

    o.dae_maxstate = 7;
    o.dae_maxcov = 3;
    o.iter_max = 10;   // below the floor, so the floor stands
    m = fluid::fluid_dae_options(o, d);
    CHECK(m.maxstate == 7);
    CHECK(m.maxcov == 3);
    CHECK(m.newton_max == 50);

    // a pinned struct against silent options
    fluid::FluidOptions q;
    q.iter_max = 0;
    d.maxstate = 4;
    d.maxcov = 0;
    m = fluid::fluid_dae_options(q, d);
    CHECK(m.maxstate == 4);
    CHECK(m.maxcov == 0);

    // and the refusal the solve itself raises names the option
    qn::Network<double> bc = balanced_cycle(6.0, 1, lang::SchedStrategy::PS);
    fluid::FluidOptions capped = dae_opts();
    capped.dae_maxstate = 1;
    CHECK(refusal([&] { fluid::solver_fluid_dae(bc.get_struct(), capped); })
              .find("dae_maxstate") != std::string::npos);
}

TEST_CASE("fluid: a binding station buffer is refused, except by the method that carries it") {
    // Nothing in the fluid tree reads `sn.cap` or `sn.classcap`, so every method
    // but `dae` -- which carries the buffer as an algebraic constraint on the
    // drift -- integrated a capped station as an UNBOUNDED one and reported more
    // jobs in the buffer than the buffer holds. `@@SolverFLD/runAnalyzer.m` has
    // refused that since 2026-08-21; this port carried no capacity gate at all,
    // while SolverMVA, SolverNC and SolverAG already called the same
    // `qn::check_binding_capacity`.
    //
    // The model is `cqn_bas_blocking`: two FCFS queues, one closed class of 2,
    // `setCap(1)` on the second and the BAS drop rule on the first.
    qn::Network<double> m("cqn_bas_blocking");
    const std::size_t q1 = m.add_queue("Queue1", lang::SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", lang::SchedStrategy::FCFS);
    const std::size_t k = m.add_closed_class("Class1", 2.0, q1);
    m.set_service(q1, k, D::exp_rate(1.0));
    m.set_service(q2, k, D::exp_rate(0.8));
    m.set_capacity(q2, 1.0);
    m.set_drop_rule(q1, k, DropStrategy::BAS);
    qn::RoutingMatrix<double> P;
    P.set(k, k, q1, q2, 1.0);
    P.set(k, k, q2, q1, 1.0);
    m.link(P);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    // `default` IS NOT IN THIS LIST, and that is the point: it stands for the
    // method that carries the buffer, so `fluid_resolve_method` sends a blocked
    // model to `dae` and the gate never sees a blind method. Asked by NAME, each
    // blind method is still refused.
    CHECK(fluid::detail::fluid_resolve_method(sn, "default") == "dae");
    for (const char* blind : {"minnormal", "closing", "matrix"}) {
        fluid::FluidOptions o;
        o.method = blind;
        const std::string e = refusal([&] { fluid::solver_fluid_run_analyzer(sn, o); });
        CHECK(e.find("SolverFLD") != std::string::npos);
        CHECK(e.find("finite station capacity") != std::string::npos);
        CHECK(e.find("Queue2") != std::string::npos);
    }

    // `dae` states the buffer as a constraint, so the capacity gate must let it
    // through. What is asserted is that it is NOT REFUSED FOR THE CAPACITY --
    // not that it returns, because whether the DAE leg runs at all depends on
    // the build (`svd_full` needs LINE_MP_USE_LAPACK, which the ad-hoc link of a
    // single TU does not carry) and on the drop rule, which DaeAnalyzer refuses
    // by name for BAS. Either of those is a different refusal, and either is
    // fine here.
    for (const char* carried : {"dae", "default"}) {
        fluid::FluidOptions dae;
        dae.method = carried;
        const std::string daeErr = refusal([&] { fluid::solver_fluid_run_analyzer(sn, dae); });
        CHECK(daeErr.find("finite station capacity") == std::string::npos);
    }

    // An UNCAPPED twin must resolve exactly as it did before: the branch fires on
    // a buffer that binds, not on every model.
    qn::Network<double> free_("cqn_plain");
    const std::size_t f1 = free_.add_queue("Queue1", lang::SchedStrategy::FCFS);
    const std::size_t f2 = free_.add_queue("Queue2", lang::SchedStrategy::FCFS);
    const std::size_t fk = free_.add_closed_class("Class1", 2.0, f1);
    free_.set_service(f1, fk, D::exp_rate(1.0));
    free_.set_service(f2, fk, D::exp_rate(0.8));
    qn::RoutingMatrix<double> Pf;
    Pf.set(fk, fk, f1, f2, 1.0);
    Pf.set(fk, fk, f2, f1, 1.0);
    free_.link(Pf);
    CHECK(fluid::detail::fluid_resolve_method(free_.get_struct(), "default") == "minnormal");
}
