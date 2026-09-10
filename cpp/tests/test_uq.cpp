/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverUQ and the Prior it expands.
 *
 * Three kinds of assertion, and they fail for different reasons:
 *
 *   THE EXPECTATION IS THE WEIGHTED SOLVE. Every aggregate here is checked
 *   against the same models solved one at a time and averaged by hand, so a
 *   defect in the substitution -- a Prior replaced at the wrong class, a design
 *   weight applied to the wrong point -- shows up as a number rather than as a
 *   crash. Solving the alternatives directly is the whole specification of UQ.
 *
 *   THE DESIGN IS A CONSTRUCTION, not a sample. The tensor product's SIZE, its
 *   weights and its ORDER are pinned, because a permuted design leaves the
 *   expectation unchanged and silently renames every row of the per-point
 *   table.
 *
 *   THE REFUSALS. A Prior that reaches any other solver, a Prior at a node the
 *   reference does not expand, a design above the cap: each must fail rather
 *   than fall through to an answer computed from the epistemic mixture, which
 *   would be a confident number for a model nobody described.
 */

#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "json.hpp"
#include "line/io/network_reader.h"
#include "line/lang/prior.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/uq/uq_dispatch.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
using json = nlohmann::json;

/** Delay -> Queue closed loop, one class of population N; the repairmen shape. */
qn::Network<double> loop(const std::string& name, const D& service, double N = 2.0,
                         double think_rate = 1.0) {
    qn::Network<double> m(name);
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", N, d);
    m.set_service(d, c, D::exp_rate(think_rate));
    m.set_service(q, c, service);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** The queue-length of the loop above under plain SolverMVA. */
double loop_qlen(const D& service, double N = 2.0) {
    qn::Network<double> m = loop("ref", service, N);
    mva::MvaOptions opt;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init).QN(1, 0);
}

uq::UqStageSolver<double> mva_stage() {
    uq::UqStageOptions o;
    o.solver = "mva";
    return uq::uq_stage_solver<double>(o);
}

}  // namespace

// ---------------------------------------------------------------------------
// The Prior itself
// ---------------------------------------------------------------------------

TEST_CASE("a discrete Prior validates its alternatives and weights") {
    std::vector<D> alts{D::exp_rate(2.0), D::exp_rate(4.0)};
    SUBCASE("well formed") {
        const D p = lang::prior_discrete<double>(alts, {0.25, 0.75});
        CHECK(p.is_prior());
        CHECK(p.type == lang::ProcessType::PRIOR);
        CHECK_FALSE(p.disabled);
        CHECK(p.prior->alternatives.size() == 2);
        CHECK_FALSE(p.prior->continuous);
    }
    SUBCASE("the weights must sum to one") {
        CHECK_THROWS_AS(lang::prior_discrete<double>(alts, {0.25, 0.25}), InputError);
    }
    SUBCASE("a negative weight is not a probability") {
        CHECK_THROWS_AS(lang::prior_discrete<double>(alts, {-0.5, 1.5}), InputError);
    }
    SUBCASE("the two vectors must agree in length") {
        CHECK_THROWS_AS(lang::prior_discrete<double>(alts, {1.0}), InputError);
    }
    SUBCASE("an empty alternative set is not a prior") {
        CHECK_THROWS_AS(lang::prior_discrete<double>(std::vector<D>(), std::vector<double>()),
                        InputError);
    }
    SUBCASE("a Prior cannot be an alternative of a Prior") {
        std::vector<D> nested{lang::prior_discrete<double>(alts, {0.5, 0.5}), D::exp_rate(1.0)};
        CHECK_THROWS_AS(lang::prior_discrete<double>(nested, {0.5, 0.5}), InputError);
    }
}

TEST_CASE("the CDF of each family matches its closed form") {
    // The values are the analytic laws, not another implementation of them.
    CHECK(lang::dist_cdf(D::exp_rate(2.0), 0.5) == doctest::Approx(1.0 - std::exp(-1.0)));
    CHECK(lang::dist_cdf(D::erlang(3.0, 2), 1.0) ==
          doctest::Approx(1.0 - std::exp(-3.0) * (1.0 + 3.0)));
    CHECK(lang::dist_cdf(D::hyperexp(0.3, 1.0, 5.0), 2.0) ==
          doctest::Approx(0.3 * (1.0 - std::exp(-2.0)) + 0.7 * (1.0 - std::exp(-10.0))));
    CHECK(lang::dist_cdf(D::det(2.0), 1.9) == doctest::Approx(0.0));
    CHECK(lang::dist_cdf(D::det(2.0), 2.0) == doctest::Approx(1.0));
    CHECK(lang::dist_cdf(D::immediate(), 0.0) == doctest::Approx(1.0));
    SUBCASE("the Uniform is the RAMP, against the reference's own defect") {
        // MATLAB's Uniform.evalCDF returns the density 1/(b-a) inside the
        // support and 0 above it; the port computes the CDF. See BUGS.md.
        const D u = D::uniform(1.0, 3.0);
        CHECK(lang::dist_cdf(u, 0.5) == doctest::Approx(0.0));
        CHECK(lang::dist_cdf(u, 2.0) == doctest::Approx(0.5));
        CHECK(lang::dist_cdf(u, 4.0) == doctest::Approx(1.0));
    }
    SUBCASE("the phase-type path agrees with the Erlang closed form") {
        // An Erlang built as a PH representation takes the matrix branch, so
        // the two agree only if 1 - pie exp(D0 x) e is right.
        const D e = D::erlang(2.0, 3);
        const D ph = D::phase_type({1.0, 0.0, 0.0},
                                   Matrix<double>{{-2.0, 2.0, 0.0}, {0.0, -2.0, 2.0},
                                                  {0.0, 0.0, -2.0}},
                                   false);
        for (double x : {0.25, 1.0, 3.0})
            CHECK(lang::dist_cdf(ph, x) == doctest::Approx(lang::dist_cdf(e, x)).epsilon(1e-10));
    }
}

TEST_CASE("the quantile inverts the CDF") {
    const D e = D::exp_rate(2.0);
    CHECK(lang::dist_quantile(e, 0.5) == doctest::Approx(std::log(2.0) / 2.0).epsilon(1e-7));
    CHECK(lang::dist_quantile(e, 0.9) == doctest::Approx(-std::log(0.1) / 2.0).epsilon(1e-7));
    const D g = D::erlang(1.0, 4);
    const double q = lang::dist_quantile(g, 0.75);
    CHECK(lang::dist_cdf(g, q) == doctest::Approx(0.75).epsilon(1e-6));
    CHECK_THROWS_AS(lang::dist_quantile(e, 0.0), InputError);
    CHECK_THROWS_AS(lang::dist_quantile(e, 1.0), InputError);
}

TEST_CASE("a Prior reports the mixture moments of its alternatives") {
    // Exp(2) and Exp(4) with equal weight: means 0.5 and 0.25.
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    CHECK(lang::prior_mean(p) == doctest::Approx(0.375));
    // Var = E[Var|D] + Var(E[X|D]) = (0.25 + 0.0625)/2 + (0.5*0.25 + 0.5*0.0625 - 0.375^2)
    const double e_var = (0.25 + 0.0625) / 2.0;
    const double var_mean = 0.5 * 0.25 + 0.5 * 0.0625 - 0.375 * 0.375;
    CHECK(lang::prior_scv(p) == doctest::Approx((e_var + var_mean) / (0.375 * 0.375)));
    CHECK(lang::prior_cdf(p, 0.5) ==
          doctest::Approx(0.5 * (1.0 - std::exp(-1.0)) + 0.5 * (1.0 - std::exp(-2.0))));
    SUBCASE("and the model layer writes them onto the struct") {
        qn::Network<double> m = loop("moments", p);
        const qn::NetworkStruct<double>& sn = m.get_struct();
        CHECK(num_traits<double>::to_double(sn.rates(1, 0)) ==
              doctest::Approx(1.0 / 0.375).epsilon(1e-9));
    }
}

// ---------------------------------------------------------------------------
// The design
// ---------------------------------------------------------------------------

TEST_CASE("a discrete Prior is expanded as given by the quadrature design") {
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.4, 0.6});
    qn::Network<double> m = loop("discrete", p);
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());

    CHECK(s.method == "quadrature");
    CHECK(s.sites.size() == 1);
    CHECK(s.sites[0].cls == 1);
    CHECK_FALSE(s.sites[0].arrival);
    CHECK(s.points.size() == 2);
    CHECK(s.weights[0] == doctest::Approx(0.4));
    CHECK(s.weights[1] == doctest::Approx(0.6));

    const double q1 = loop_qlen(D::exp_rate(2.0)), q2 = loop_qlen(D::exp_rate(4.0));
    CHECK(s.points[0].QN(1, 0) == doctest::Approx(q1));
    CHECK(s.points[1].QN(1, 0) == doctest::Approx(q2));
    CHECK(s.avg.QN(1, 0) == doctest::Approx(0.4 * q1 + 0.6 * q2));
    // The system-level vectors ride along with the station matrices.
    CHECK(s.avg.XN.size() == 1);
}

TEST_CASE("the samples knob does not touch a discrete Prior") {
    // `n` is ignored by the quadrature design of a discrete prior: the set is
    // already exact, and 50 nodes would be 50 copies of two models.
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    qn::Network<double> m = loop("nodes", p);
    uq::UqOptions o;
    o.samples = 50;
    CHECK(uq::solver_uq_run_analyzer<double>(m, mva_stage(), o).points.size() == 2);
}

TEST_CASE("a continuous Prior is placed at the medians of equal-mass strata") {
    // The rate is Erlang(4, 2)-distributed and the service is Exp of that rate,
    // i.e. Prior.fromSample(2, 4): two observations summing to 4.
    const D p = lang::prior_from_sample<double>(2, 4.0);
    CHECK(p.prior->continuous);
    qn::Network<double> m = loop("continuous", p);
    uq::UqOptions o;
    o.samples = 5;
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage(), o);

    REQUIRE(s.points.size() == 5);
    const D param = D::erlang(4.0, 2);
    for (std::size_t i = 0; i < 5; ++i) {
        CHECK(s.weights[i] == doctest::Approx(0.2));
        // Node i is the ((i+0.5)/n)-quantile of the parameter law, and the
        // alternative is an Exp of THAT rate.
        const double theta = lang::dist_quantile(param, (static_cast<double>(i) + 0.5) / 5.0);
        CHECK(s.design[i].dists[0].params[0] == doctest::Approx(theta).epsilon(1e-6));
        CHECK(s.points[i].QN(1, 0) == doctest::Approx(loop_qlen(D::exp_rate(theta))).epsilon(1e-9));
    }
    double acc = 0.0;
    for (std::size_t i = 0; i < 5; ++i) acc += 0.2 * s.points[i].QN(1, 0);
    CHECK(s.avg.QN(1, 0) == doctest::Approx(acc));
}

TEST_CASE("several Priors take the tensor product, first coordinate fastest") {
    const D ps = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.4, 0.6});
    const D pz = lang::prior_discrete<double>(
        {D::exp_rate(1.0), D::exp_rate(2.0), D::exp_rate(4.0)}, {0.2, 0.3, 0.5});

    qn::Network<double> m("tensor");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, pz);
    m.set_service(q, c, ps);
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());
    REQUIRE(s.sites.size() == 2);
    // Detection is in NODE order: the Delay is node 1, the Queue node 2.
    CHECK(s.sites[0].node == d);
    CHECK(s.sites[1].node == q);
    REQUIRE(s.points.size() == 6);
    // Site 0 (the Delay's 3 alternatives) varies FASTEST.
    const double wz[3] = {0.2, 0.3, 0.5};
    const double ws[2] = {0.4, 0.6};
    double tot = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        CHECK(s.weights[i] == doctest::Approx(wz[i % 3] * ws[i / 3]));
        tot += s.weights[i];
    }
    CHECK(tot == doctest::Approx(1.0));
}

TEST_CASE("a design above the cap is refused, not truncated") {
    const D p = lang::prior_from_sample<double>(2, 4.0);
    qn::Network<double> m("cap");
    const std::size_t dl = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 2.0, dl);
    m.set_service(dl, c, p);
    m.set_service(q, c, p);
    qn::RoutingMatrix<double> P;
    P.set(dl, q, 1.0);
    P.set(q, dl, 1.0);
    m.link(P);

    uq::UqOptions o;
    o.samples = 100;  // 100 x 100 = 10000 > 4096
    CHECK_THROWS_AS(uq::solver_uq_run_analyzer<double>(m, mva_stage(), o), UnsupportedError);
    o.method = "montecarlo";  // the same request, at a cost that does not grow with L
    CHECK(uq::solver_uq_run_analyzer<double>(m, mva_stage(), o).points.size() == 100);
}

TEST_CASE("the Monte Carlo design draws n points of equal weight, reproducibly") {
    const D p = lang::prior_from_sample<double>(3, 6.0);
    qn::Network<double> m = loop("mc", p);
    uq::UqOptions o;
    o.method = "montecarlo";
    o.samples = 8;
    o.seed = 12345;
    const uq::UqSolution<double> a = uq::solver_uq_run_analyzer<double>(m, mva_stage(), o);
    const uq::UqSolution<double> b = uq::solver_uq_run_analyzer<double>(m, mva_stage(), o);
    REQUIRE(a.points.size() == 8);
    CHECK(a.method == "montecarlo");
    for (std::size_t i = 0; i < 8; ++i) {
        CHECK(a.weights[i] == doctest::Approx(0.125));
        CHECK(a.design[i].dists[0].params[0] == doctest::Approx(b.design[i].dists[0].params[0]));
    }
    o.seed = 999;
    const uq::UqSolution<double> cdiff = uq::solver_uq_run_analyzer<double>(m, mva_stage(), o);
    bool moved = false;
    for (std::size_t i = 0; i < 8; ++i)
        if (std::fabs(cdiff.design[i].dists[0].params[0] - a.design[i].dists[0].params[0]) > 1e-12)
            moved = true;
    CHECK(moved);
}

TEST_CASE("a model with no Prior is solved once, and the expectation is that solve") {
    qn::Network<double> m = loop("noprior", D::exp_rate(3.0));
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());
    CHECK(s.sites.empty());
    CHECK(s.points.size() == 1);
    CHECK(s.weights[0] == doctest::Approx(1.0));
    CHECK(s.avg.QN(1, 0) == doctest::Approx(loop_qlen(D::exp_rate(3.0))));
}

// ---------------------------------------------------------------------------
// The arrival side, the wire, and the other arithmetics
// ---------------------------------------------------------------------------

TEST_CASE("a Prior on a Source is expanded as an ARRIVAL process") {
    const D p = lang::prior_discrete<double>({D::exp_rate(1.0), D::exp_rate(1.5)}, {0.5, 0.5});
    qn::Network<double> m("openprior");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, p);
    m.set_service(q, o, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);

    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());
    REQUIRE(s.sites.size() == 1);
    CHECK(s.sites[0].arrival);
    REQUIRE(s.points.size() == 2);
    // M/M/1 queue lengths at rho = 1/4 and 3/8.
    CHECK(s.points[0].QN(1, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(s.points[1].QN(1, 0) == doctest::Approx(0.6).epsilon(1e-9));
    CHECK(s.avg.QN(1, 0) == doctest::Approx(0.5 * (1.0 / 3.0) + 0.5 * 0.6).epsilon(1e-9));
}

TEST_CASE("a Prior crosses model.json in its discrete form") {
    const std::string text = R"({
      "format":"line-model","version":"1.0","model":{
       "type":"Network","name":"jsonprior",
       "nodes":[
        {"name":"Source 1","type":"Source","service":{
           "Class1":{"type":"Exp","params":{"lambda":1.0}}}},
        {"name":"Queue 1","type":"Queue","scheduling":"FCFS","service":{
           "Class1":{"type":"Prior",
                     "distributions":[{"type":"Exp","params":{"lambda":4.0}},
                                      {"type":"Exp","params":{"lambda":2.0}}],
                     "probabilities":[0.25,0.75]}}},
        {"name":"Sink 1","type":"Sink"}
       ],
       "classes":[{"name":"Class1","type":"Open"}],
       "routing":{"type":"matrix","matrix":{
         "Class1,Class1":{"Source 1":{"Queue 1":1.0},"Queue 1":{"Sink 1":1.0}}}}
      }})";
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(net, mva_stage());
    REQUIRE(s.sites.size() == 1);
    REQUIRE(s.points.size() == 2);
    // M/M/1 at rho = 1/4 and 1/2.
    CHECK(s.avg.QN(1, 0) == doctest::Approx(0.25 * (1.0 / 3.0) + 0.75 * 1.0).epsilon(1e-9));
    SUBCASE("and a Prior missing its alternative set is refused, not half read") {
        const std::string bad = R"({"type":"Prior","probabilities":[1.0]})";
        CHECK_THROWS_AS(io::detail::dist_from_json<double>(json::parse(bad)), InputError);
    }
}

TEST_CASE("a continuous Prior crosses model.json as a factory template") {
    // The wire form of `Prior(Erlang(4,2), @(lambda) Exp(lambda))`: the factory
    // handle cannot cross, the distribution it BUILDS can, with the parameter
    // left in the slot the handle fills.
    const std::string text = R"({
      "format":"line-model","version":"1.0","model":{
       "type":"Network","name":"jsoncontinuous",
       "nodes":[
        {"name":"Source 1","type":"Source","service":{
           "Class1":{"type":"Exp","params":{"lambda":1.0}}}},
        {"name":"Queue 1","type":"Queue","scheduling":"FCFS","service":{
           "Class1":{"type":"Prior","kind":"continuous",
                     "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}},
                     "factory":{"template":{"type":"Exp","params":{"lambda":1.0}},
                                "slots":["lambda"]}}}},
        {"name":"Sink 1","type":"Sink"}
       ],
       "classes":[{"name":"Class1","type":"Open"}],
       "routing":{"type":"matrix","matrix":{
         "Class1,Class1":{"Source 1":{"Queue 1":1.0},"Queue 1":{"Sink 1":1.0}}}}
      }})";
    qn::Network<double> net = io::build_network_from_json<double>(json::parse(text));
    uq::UqOptions o;
    o.samples = 5;
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(net, mva_stage(), o);

    REQUIRE(s.sites.size() == 1);
    CHECK(s.sites[0].prior.prior->continuous);
    REQUIRE(s.points.size() == 5);
    // THE WIRE FORM AND THE API FORM ARE THE SAME PRIOR, node for node: the
    // parameter reaches the template's slot, so the design is the one
    // `prior_from_sample(2, 4)` builds in process and options.samples still
    // chooses the node count -- which a materialized alternative set could not.
    const D api = lang::prior_from_sample<double>(2, 4.0);
    lang::PriorRng rng(o.seed);
    const lang::PriorDesign<double> want = lang::prior_discretize(api, 5, "quadrature", rng);
    for (std::size_t i = 0; i < 5; ++i) {
        CHECK(s.weights[i] == doctest::Approx(0.2));
        CHECK(s.design[i].dists[0].type == lang::ProcessType::EXP);
        CHECK(s.design[i].dists[0].params[0] ==
              doctest::Approx(want.dists[i].params[0]).epsilon(1e-12));
    }

    SUBCASE("a second node count re-discretizes rather than replaying a table") {
        uq::UqOptions o2;
        o2.samples = 9;
        const uq::UqSolution<double> s9 = uq::solver_uq_run_analyzer<double>(net, mva_stage(), o2);
        CHECK(s9.points.size() == 9);
        for (std::size_t i = 0; i < 9; ++i) CHECK(s9.weights[i] == doctest::Approx(1.0 / 9.0));
    }

    SUBCASE("every way the encoding can be incomplete is named") {
        const char* bad[] = {
            R"({"type":"Prior","kind":"sampled","distributions":[],"probabilities":[]})",
            R"({"type":"Prior","kind":"continuous",
                "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}}})",
            R"({"type":"Prior","kind":"continuous",
                "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}},
                "factory":{"template":{"type":"Exp","params":{"lambda":1.0}}}})",
            R"({"type":"Prior","kind":"continuous",
                "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}},
                "factory":{"template":{"type":"Exp","params":{"lambda":1.0}},"slots":[]}})",
            R"({"type":"Prior","kind":"continuous",
                "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}},
                "factory":{"template":{"type":"Exp","params":{"lambda":1.0}},"slots":["rate"]}})",
            R"({"type":"Prior","kind":"continuous",
                "paramDist":{"type":"Erlang","params":{"lambda":4.0,"k":2}},
                "factory":{"template":{"type":"APH","ph":{"alpha":[1.0],"T":[[-1.0]]}},
                           "slots":["lambda"]}})"};
        for (std::size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i)
            CHECK_THROWS_AS(io::detail::dist_from_json<double>(json::parse(bad[i])), InputError);
    }
}

// ---------------------------------------------------------------------------
// The ensemble lifecycle
// ---------------------------------------------------------------------------

TEST_CASE("the ensemble can be driven point by point, and reports the same aggregate") {
    const D p = lang::prior_discrete<double>(
        {D::exp_rate(2.0), D::exp_rate(3.0), D::exp_rate(6.0)}, {0.2, 0.3, 0.5});
    qn::Network<double> m = loop("lifecycle", p);
    uq::SolverUq<double> s(m, mva_stage());

    CHECK(s.has_prior_distribution());
    CHECK(s.get_num_alternatives() == 3);
    CHECK(s.get_number_of_models() == 3);
    CHECK(s.get_uq_method() == "quadrature");
    CHECK(s.get_uq_nodes() == lang::kPriorDefaultNodes);
    REQUIRE(s.get_probabilities().size() == 3);
    CHECK(s.get_probabilities()[1] == doctest::Approx(0.3));

    // THE DESIGN IS KNOWN BEFORE ANYTHING IS SOLVED, which is the point of the
    // lifecycle: `init()` has run, and no stage solve has.
    CHECK(s.get_ensemble_avg().empty());

    const int it = 1;
    s.pre(it);
    for (std::size_t e = 1; e <= s.get_number_of_models(); ++e) {
        const mva::AvgResult<double>& r = s.analyze(it, e);
        CHECK(r.QN(1, 0) == doctest::Approx(loop_qlen(s.get_design()[e - 1].dists[0])));
    }
    s.post(it);
    CHECK(s.converged(it));
    s.finish();

    const double want = 0.2 * loop_qlen(D::exp_rate(2.0)) + 0.3 * loop_qlen(D::exp_rate(3.0)) +
                        0.5 * loop_qlen(D::exp_rate(6.0));
    CHECK(s.get_avg().QN(1, 0) == doctest::Approx(want).epsilon(1e-12));
    // The one-shot spelling is the same lifecycle, so it is the same numbers.
    CHECK(uq::solver_uq_run_analyzer<double>(m, mva_stage()).avg.QN(1, 0) ==
          doctest::Approx(s.get_avg().QN(1, 0)).epsilon(1e-15));

    SUBCASE("a point outside the design is refused") {
        CHECK_THROWS_AS(s.analyze(1, 0), InputError);
        CHECK_THROWS_AS(s.analyze(1, 4), InputError);
        CHECK_THROWS_AS(s.expand(4), InputError);
    }
}

TEST_CASE("the expanded model of a design point carries the alternative, not the Prior") {
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(5.0)}, {0.5, 0.5});
    qn::Network<double> m = loop("expand", p);
    uq::SolverUq<double> s(m, mva_stage());
    REQUIRE(s.get_number_of_models() == 2);
    for (std::size_t e = 1; e <= 2; ++e) {
        qn::Network<double> alt = s.expand(e);
        const qn::NetworkStruct<double>& sn = alt.get_struct();
        // Station 2 (the Queue) serves the class at the alternative's rate, and
        // the model that reaches the stage solver has no Prior left in it.
        CHECK(sn.rates(1, 0) == doctest::Approx(e == 1 ? 2.0 : 5.0));
        CHECK_FALSE(sn.service[1][0].is_prior());
    }
    // The original model is UNTOUCHED: the expansion copies.
    CHECK(m.get_struct().service[1][0].is_prior());
}

TEST_CASE("a model with no Prior has an ensemble of one") {
    qn::Network<double> m = loop("noprior-lifecycle", D::exp_rate(3.0));
    uq::SolverUq<double> s(m, mva_stage());
    CHECK_FALSE(s.has_prior_distribution());
    CHECK(s.get_number_of_models() == 1);
    CHECK(s.get_probabilities()[0] == doctest::Approx(1.0));
    s.iterate();
    CHECK(s.get_avg().QN(1, 0) == doctest::Approx(loop_qlen(D::exp_rate(3.0))));
}

TEST_CASE("a discrete Prior expands under exact arithmetic") {
    // The quadrature design of a discrete prior draws nothing and evaluates no
    // quantile, so it carries through to the rational backend unchanged.
    using R = Rational;
    qn::Network<R> m("exactprior");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, Distrib<R>::exp_rate(num_traits<R>::from_int(1)));
    std::vector<Distrib<R>> alts{Distrib<R>::exp_rate(num_traits<R>::from_int(2)),
                                 Distrib<R>::exp_rate(num_traits<R>::from_int(4))};
    std::vector<R> pr{num_traits<R>::from_double(0.5), num_traits<R>::from_double(0.5)};
    m.set_service(q, c, lang::prior_discrete<R>(alts, pr));
    qn::RoutingMatrix<R> P;
    P.set(d, q, num_traits<R>::from_int(1));
    P.set(q, d, num_traits<R>::from_int(1));
    m.link(P);

    uq::UqStageOptions so;
    so.solver = "mva";
    const uq::UqSolution<R> s = uq::solver_uq_run_analyzer<R>(m, uq::uq_stage_solver<R>(so));
    REQUIRE(s.points.size() == 2);
    // 0.8 and 6/13, averaged: the same numbers the double run reports.
    CHECK(num_traits<R>::to_double(s.avg.QN(1, 0)) ==
          doctest::Approx(0.5 * (0.8 + 6.0 / 13.0)).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// Posterior summaries
// ---------------------------------------------------------------------------

TEST_CASE("the posterior summaries read the design, not the expectation") {
    const D p = lang::prior_discrete<double>(
        {D::exp_rate(2.0), D::exp_rate(3.0), D::exp_rate(6.0)}, {0.2, 0.3, 0.5});
    qn::Network<double> m = loop("posterior", p);
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());

    const std::vector<double> q = uq::uq_samples(s, "Q", 2, 1);
    REQUIRE(q.size() == 3);
    CHECK(q[0] == doctest::Approx(loop_qlen(D::exp_rate(2.0))));
    CHECK(q[2] == doctest::Approx(loop_qlen(D::exp_rate(6.0))));

    const uq::UqMoments<double> mo = uq::uq_moments(s, "Q", 2, 1);
    double mean = 0.0, var = 0.0;
    const double w[3] = {0.2, 0.3, 0.5};
    for (std::size_t i = 0; i < 3; ++i) mean += w[i] * q[i];
    for (std::size_t i = 0; i < 3; ++i) var += w[i] * (q[i] - mean) * (q[i] - mean);
    CHECK(mo.mean == doctest::Approx(mean));
    CHECK(mo.var == doctest::Approx(var));
    CHECK(mo.mean == doctest::Approx(s.avg.QN(1, 0)));

    SUBCASE("the empirical CDF is sorted and its weights are the design's") {
        const uq::UqEmpiricalCdf<double> ec = uq::uq_posterior_cdf(s, "Q", 2, 1);
        REQUIRE(ec.values.size() == 3);
        CHECK(ec.values[0] <= ec.values[1]);
        CHECK(ec.values[1] <= ec.values[2]);
        CHECK(ec.cdf[2] == doctest::Approx(1.0));
        // The smallest queue length is the fastest server, weight 0.5.
        CHECK(ec.probabilities[0] == doctest::Approx(0.5));
    }
    SUBCASE("the credible interval names two design points") {
        const std::pair<double, double> ci = uq::uq_credible_interval(s, "Q", 2, 1, 0.95);
        const uq::UqEmpiricalCdf<double> ec = uq::uq_posterior_cdf(s, "Q", 2, 1);
        CHECK(ci.first == doctest::Approx(ec.values.front()));
        CHECK(ci.second == doctest::Approx(ec.values.back()));
        CHECK_THROWS_AS(uq::uq_credible_interval(s, "Q", 2, 1, 1.0), InputError);
    }
    SUBCASE("an unknown metric is refused") {
        CHECK_THROWS_AS(uq::uq_samples(s, "Z", 2, 1), InputError);
    }
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

TEST_CASE("a Prior that reaches any other solver is refused at the feature gate") {
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    qn::Network<double> m = loop("gated", p);
    mva::MvaOptions opt;
    Matrix<double> init;
    CHECK_THROWS_AS(mva::solver_mva_run_analyzer(m.get_struct(), opt, init), UnsupportedError);
    SUBCASE("and the used-feature set names it") {
        const qn::FeatureSet u = qn::used_lang_features(m.get_struct());
        CHECK(u.has(qn::Feature::Prior));
    }
    SUBCASE("while SolverUQ declares it") {
        CHECK(uq::uq_feature_set().has(qn::Feature::Prior));
    }
}

TEST_CASE("a Prior has no representation to lower onto") {
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    CHECK_THROWS_AS(lang::dist_to_map(p), UnsupportedError);
    CHECK_THROWS_AS(lang::dist_lst(p, 1.0), UnsupportedError);
    CHECK_THROWS_AS(lang::dist_moment(p, 2), UnsupportedError);
    CHECK_THROWS_AS(lang::dist_cdf(p, 1.0), UnsupportedError);
}

TEST_CASE("SolverUQ refuses what it cannot expand") {
    SUBCASE("no stage solver") {
        const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
        qn::Network<double> m = loop("nostage", p);
        CHECK_THROWS_AS(uq::solver_uq_run_analyzer<double>(m, uq::UqStageSolver<double>()), InputError);
    }
    SUBCASE("an unknown discretization method") {
        const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
        qn::Network<double> m = loop("badmethod", p);
        uq::UqOptions o;
        o.method = "latin";
        CHECK_THROWS_AS(uq::solver_uq_run_analyzer<double>(m, mva_stage(), o), InputError);
    }
    SUBCASE("an unknown stage solver") {
        uq::UqStageOptions o;
        o.solver = "ldes";
        CHECK_THROWS_AS(uq::uq_stage_solver<double>(o), UnsupportedError);
        o.solver = "";
        CHECK_THROWS_AS(uq::uq_stage_solver<double>(o), InputError);
    }
    SUBCASE("a Prior at a node the reference does not expand") {
        // The detection walks the service table, so a Prior parked on a station
        // that is neither a Queue, a Delay nor a Source must be NAMED rather
        // than left in place for the stage solver to lower onto a rate.
        const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
        qn::Network<double> m = loop("wrongnode", p);
        m.get_struct();
        m.raw_struct().nodes[1].nodetype = qn::NodeType::Cache;
        CHECK_THROWS_AS(uq::uq_detect_priors(m.raw_struct()), UnsupportedError);
    }
}

// ---------------------------------------------------------------------------
// Support-only (interval) uncertainty
// ---------------------------------------------------------------------------

TEST_CASE("the exact interval is the attained hull of MVA over the demand box") {
    // Think Exp(1), service uncertain between mean 0.25 and mean 0.5, N = 2.
    // The whole design lies inside the box, so the hull must CONTAIN every
    // design point and be attained at the two endpoints -- which for a discrete
    // Prior on two alternatives means it equals the sampled range exactly.
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    qn::Network<double> m = loop("interval", p);
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());
    const uq::UqInterval<double> iv = uq::uq_interval(s, m.get_struct());

    CHECK(iv.exact);
    CHECK(iv.method == "mvainterval");
    CHECK(iv.why.empty());
    CHECK(iv.has_totals);
    const double qlo = loop_qlen(D::exp_rate(4.0)), qup = loop_qlen(D::exp_rate(2.0));
    CHECK(iv.Qlo(1, 0) == doctest::Approx(qlo).epsilon(1e-9));
    CHECK(iv.Qup(1, 0) == doctest::Approx(qup).epsilon(1e-9));
    // The throughput moves the other way: a slower server carries fewer jobs.
    CHECK(iv.Xlo <= iv.Xup);
    CHECK(iv.Tlo(1, 0) <= iv.Tup(1, 0));
    // Every design point lies inside the hull, which is what an enclosure means.
    for (const mva::AvgResult<double>& r : s.points) {
        CHECK(r.QN(1, 0) >= iv.Qlo(1, 0) - 1e-9);
        CHECK(r.QN(1, 0) <= iv.Qup(1, 0) + 1e-9);
    }
    SUBCASE("and the utilization enclosure never exceeds one") {
        CHECK(iv.Uup(1, 0) <= 1.0);
        CHECK(iv.Ulo(1, 0) <= iv.Uup(1, 0));
    }
}

TEST_CASE("a continuous Prior's exact interval is wider than its sampled range") {
    // The quadrature nodes are stratum MEDIANS, so the sampled range stops
    // short of the discretized support's endpoints only when the design is
    // coarser than the range used for the box -- here both use the same node
    // count, so the two agree and the point is that the exact path needs no
    // ensemble at all.
    const D p = lang::prior_from_sample<double>(2, 4.0);
    qn::Network<double> m = loop("contint", p);
    uq::UqOptions o;
    o.samples = 5;
    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage(), o);
    const uq::UqInterval<double> iv = uq::uq_interval(s, m.get_struct());
    REQUIRE(iv.exact);
    for (const mva::AvgResult<double>& r : s.points) {
        CHECK(r.QN(1, 0) >= iv.Qlo(1, 0) - 1e-9);
        CHECK(r.QN(1, 0) <= iv.Qup(1, 0) + 1e-9);
    }
}

TEST_CASE("a model outside the monotonicity theorems falls back and says why") {
    // An open model: the Prior sits on a Source arrival, which is the first
    // condition qualifiesForIntervalMVA tests.
    const D p = lang::prior_discrete<double>({D::exp_rate(1.0), D::exp_rate(1.5)}, {0.5, 0.5});
    qn::Network<double> m("openint");
    const std::size_t src = m.add_source("Source");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t o = m.add_open_class("O1");
    m.set_arrival(src, o, p);
    m.set_service(q, o, D::exp_rate(4.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);

    const uq::UqSolution<double> s = uq::solver_uq_run_analyzer<double>(m, mva_stage());
    const uq::UqInterval<double> iv = uq::uq_interval(s, m.get_struct());
    CHECK_FALSE(iv.exact);
    CHECK(iv.method == "sampled");
    CHECK(iv.why.find("arrival") != std::string::npos);
    // The sampled range spans the two solved points and nothing more.
    CHECK(iv.Qlo(1, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(iv.Qup(1, 0) == doctest::Approx(0.6).epsilon(1e-9));
    CHECK_FALSE(iv.has_totals);
}

TEST_CASE("each disqualifying condition is named by the one that fires first") {
    const D p = lang::prior_discrete<double>({D::exp_rate(2.0), D::exp_rate(4.0)}, {0.5, 0.5});
    SUBCASE("two classes") {
        qn::Network<double> m("twoclass");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Queue", SchedStrategy::PS);
        const std::size_t c1 = m.add_closed_class("C1", 2.0, d);
        const std::size_t c2 = m.add_closed_class("C2", 1.0, d);
        m.set_service(d, c1, D::exp_rate(1.0));
        m.set_service(d, c2, D::exp_rate(1.0));
        m.set_service(q, c1, p);
        m.set_service(q, c2, D::exp_rate(3.0));
        qn::RoutingMatrix<double> P;
        P.set(c1, c1, d, q, 1.0);
        P.set(c1, c1, q, d, 1.0);
        P.set(c2, c2, d, q, 1.0);
        P.set(c2, c2, q, d, 1.0);
        m.link(P);
        const std::pair<bool, std::string> r =
            uq::uq_qualifies_for_interval_mva(m.get_struct(), uq::uq_detect_priors(m.get_struct()));
        CHECK_FALSE(r.first);
        CHECK(r.second.find("single class") != std::string::npos);
    }
    SUBCASE("a multiserver queue") {
        qn::Network<double> m = loop("multiserver", p);
        m.set_number_of_servers(2, 3.0);
        const std::pair<bool, std::string> r =
            uq::uq_qualifies_for_interval_mva(m.get_struct(), uq::uq_detect_priors(m.get_struct()));
        CHECK_FALSE(r.first);
        CHECK(r.second.find("server") != std::string::npos);
    }
    SUBCASE("a discipline the theorems do not cover") {
        qn::Network<double> m = loop("siro", p);
        m.get_struct();
        m.raw_struct().stations[1].sched = SchedStrategy::SIRO;
        const std::pair<bool, std::string> r =
            uq::uq_qualifies_for_interval_mva(m.raw_struct(), uq::uq_detect_priors(m.raw_struct()));
        CHECK_FALSE(r.first);
        CHECK(r.second.find("delay, PS nor FCFS") != std::string::npos);
    }
}

TEST_CASE("pfqn_mva_interval reduces to pfqn_mva on a thin box") {
    // A box of zero width is one model, so every endpoint must coincide with
    // the ordinary MVA answer; without that the corner selection is untestable.
    Matrix<double> L(2, 2);
    L(0, 0) = 0.6; L(0, 1) = 0.6;
    L(1, 0) = 0.4; L(1, 1) = 0.4;
    const pfqn::MvaIntervalResult<double> iv = pfqn::pfqn_mva_interval(L, 3, 3, 1.0, 1.0);
    Matrix<double> D2(2, 1);
    D2(0, 0) = 0.6;
    D2(1, 0) = 0.4;
    const pfqn::MvaResult<double> r =
        pfqn::pfqn_mva(D2, std::vector<int>(1, 3), Matrix<double>(1, 1, 1.0));
    CHECK(iv.Xlo == doctest::Approx(r.XN[0]).epsilon(1e-12));
    CHECK(iv.Xup == doctest::Approx(r.XN[0]).epsilon(1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(iv.Q(i, 0) == doctest::Approx(r.QN(i, 0)).epsilon(1e-12));
        CHECK(iv.Q(i, 1) == doctest::Approx(r.QN(i, 0)).epsilon(1e-12));
        CHECK(iv.R(i, 0) == doctest::Approx(r.CN(i, 0)).epsilon(1e-12));
    }
    SUBCASE("and a thick box encloses every corner of it") {
        Matrix<double> W(2, 2);
        W(0, 0) = 0.4; W(0, 1) = 0.8;
        W(1, 0) = 0.4; W(1, 1) = 0.4;
        const pfqn::MvaIntervalResult<double> wv = pfqn::pfqn_mva_interval(W, 3, 3, 1.0, 1.0);
        CHECK(wv.Xlo <= iv.Xlo);
        CHECK(wv.Xup >= iv.Xup);
        CHECK(wv.Q(0, 0) <= iv.Q(0, 0));
        CHECK(wv.Q(0, 1) >= iv.Q(0, 1));
    }
    SUBCASE("an inverted or empty box is refused") {
        Matrix<double> bad(1, 2);
        bad(0, 0) = 0.8; bad(0, 1) = 0.4;
        CHECK_THROWS_AS(pfqn::pfqn_mva_interval(bad, 2, 2, 0.0, 0.0), InputError);
        CHECK_THROWS_AS(pfqn::pfqn_mva_interval(Matrix<double>(), 2, 2, 0.0, 0.0), InputError);
        Matrix<double> ok(1, 2);
        ok(0, 0) = 0.4; ok(0, 1) = 0.8;
        CHECK_THROWS_AS(pfqn::pfqn_mva_interval(ok, 0, 2, 0.0, 0.0), InputError);
    }
}

TEST_CASE("pfqn_mva_interval reproduces the published figures of Luthi and Haring") {
    // Section 3.2 of the paper, and the same box the MATLAB, JAR and Python
    // ports are pinned to: D_cpu = [12,16], D_disk = 10, Z = [15,20], n = 10.
    // The paper prints the throughput endpoints; the remaining numbers are the
    // cross-codebase golden, so a corner picked on the wrong side of a
    // monotonicity sign shows up here rather than as a plausible-looking hull.
    Matrix<double> L(2, 2);
    L(0, 0) = 12.0; L(0, 1) = 16.0;
    L(1, 0) = 10.0; L(1, 1) = 10.0;
    const pfqn::MvaIntervalResult<double> iv = pfqn::pfqn_mva_interval(L, 10, 10, 15.0, 20.0);

    // the paper prints four decimals, the golden carries the full precision
    CHECK(iv.Xlo == doctest::Approx(0.0620).epsilon(1e-3));
    CHECK(iv.Xup == doctest::Approx(0.0799).epsilon(1e-3));
    CHECK(iv.Xlo == doctest::Approx(0.06204324).epsilon(1e-7));
    CHECK(iv.Xup == doctest::Approx(0.07985136).epsilon(1e-7));
    CHECK(iv.Q(0, 0) == doctest::Approx(5.49188426).epsilon(1e-7));
    CHECK(iv.Q(0, 1) == doctest::Approx(7.49723174).epsilon(1e-7));
    CHECK(iv.Q(1, 0) == doctest::Approx(1.55704371).epsilon(1e-7));
    CHECK(iv.Q(1, 1) == doctest::Approx(3.01527403).epsilon(1e-7));
    CHECK(iv.R(0, 0) == doctest::Approx(69.09872305).epsilon(1e-7));
    CHECK(iv.R(0, 1) == doctest::Approx(120.68538538).epsilon(1e-7));
    CHECK(iv.R(1, 0) == doctest::Approx(25.09610503).epsilon(1e-7));
    CHECK(iv.R(1, 1) == doctest::Approx(37.76108548).epsilon(1e-7));
    CHECK(iv.Rtot_lo == doctest::Approx(105.81970016).epsilon(1e-7));
    CHECK(iv.Rtot_up == doctest::Approx(145.97326256).epsilon(1e-7));
    CHECK(iv.Qtot_lo == doctest::Approx(8.41042381).epsilon(1e-7));
    CHECK(iv.Qtot_up == doctest::Approx(9.06816823).epsilon(1e-7));
    // The CPU's demand is thick, so its utilization enclosure saturates.
    CHECK(iv.U(0, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(iv.U(1, 1) == doctest::Approx(0.7985136).epsilon(1e-7));

    SUBCASE("the hull contains and attains the grid over the box") {
        double xlo = 1e300, xup = -1e300, rtlo = 1e300, rtup = -1e300;
        std::vector<double> qlo(2, 1e300), qup(2, -1e300);
        for (int a = 0; a <= 8; ++a) {
            for (int b = 0; b <= 5; ++b) {
                Matrix<double> d(2, 1);
                d(0, 0) = 12.0 + 0.5 * a;
                d(1, 0) = 10.0;
                const pfqn::MvaResult<double> r = pfqn::pfqn_mva(
                    d, std::vector<int>(1, 10), Matrix<double>(1, 1, 15.0 + b));
                xlo = std::min(xlo, r.XN[0]);
                xup = std::max(xup, r.XN[0]);
                double rt = 0.0;
                for (std::size_t i = 0; i < 2; ++i) {
                    qlo[i] = std::min(qlo[i], r.QN(i, 0));
                    qup[i] = std::max(qup[i], r.QN(i, 0));
                    rt += r.CN(i, 0);
                }
                rtlo = std::min(rtlo, rt);
                rtup = std::max(rtup, rt);
            }
        }
        CHECK(xlo == doctest::Approx(iv.Xlo).epsilon(1e-12));
        CHECK(xup == doctest::Approx(iv.Xup).epsilon(1e-12));
        CHECK(rtlo == doctest::Approx(iv.Rtot_lo).epsilon(1e-12));
        CHECK(rtup == doctest::Approx(iv.Rtot_up).epsilon(1e-12));
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(qlo[i] == doctest::Approx(iv.Q(i, 0)).epsilon(1e-12));
            CHECK(qup[i] == doctest::Approx(iv.Q(i, 1)).epsilon(1e-12));
        }
    }

    SUBCASE("a population interval widens the throughput on the lower side only") {
        const pfqn::MvaIntervalResult<double> wv = pfqn::pfqn_mva_interval(L, 5, 10, 15.0, 20.0);
        Matrix<double> d(2, 1);
        d(0, 0) = 16.0;
        d(1, 0) = 10.0;
        const pfqn::MvaResult<double> r =
            pfqn::pfqn_mva(d, std::vector<int>(1, 5), Matrix<double>(1, 1, 20.0));
        CHECK(wv.Xlo == doctest::Approx(r.XN[0]).epsilon(1e-12));
        CHECK(wv.Xup == doctest::Approx(iv.Xup).epsilon(1e-12));
        CHECK(wv.Rtot_lo <= iv.Rtot_lo + 1e-12);
    }
}
