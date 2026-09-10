/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The SSJ variate layer against SSJ itself.
 *
 * Every expected value below was printed at 17 significant digits by a Java
 * program built on `common/jline.jar`, which bundles SSJ, driving the same
 * `*Gen` classes `Solver_ssj` builds from an MRG32k3a seeded the way the engine
 * seeds it (`{seed+offset .. seed+offset+5}`, 23000 + 1000). The C++ side draws
 * from `rng::Mrg32k3a` with the same seed, so the two are compared draw by
 * draw on the same uniforms.
 *
 * WHAT THE TOLERANCES MEAN, and why they are not uniform. The closed-form
 * quantiles (Exponential, Uniform, Pareto, Weibull, Bernoulli) are required to
 * agree EXACTLY: they are the same arithmetic on the same uniform, and a
 * difference would be a transcription error, not rounding. Lognormal, Gamma and
 * Erlang go through AS 241 and a Newton inversion, so they are held to 1e-13
 * relative -- far below the resolution at which a service time could reorder
 * two events, and the stream cannot desynchronize regardless because each draw
 * costs exactly one uniform. Poisson and Binomial return integers and are
 * required to be equal.
 *
 * THE UNIFORM COUNT IS ITSELF PINNED at the end: if any family ever consumed
 * two uniforms where SSJ consumes one, every later draw in the run would shift
 * and the seeded goldens would silently move. That test is the load-bearing
 * one.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/solvers/ldes/ldes_ssj_variates.h"
#include "line/util/rng_ssj.h"

using namespace line;
namespace v = line::ldes::ssj;

namespace {

rng::Mrg32k3a engine_stream() {
    rng::Mrg32k3a s;
    s.set_seed_offset(23000, 1000);
    return s;
}

/** The first five uniforms of that stream, which every table below is keyed to. */
std::vector<double> five_uniforms() {
    rng::Mrg32k3a s = engine_stream();
    std::vector<double> u(5);
    for (int i = 0; i < 5; ++i) u[i] = s.next_double();
    return u;
}

}  // namespace

TEST_CASE("ldes ssj variates: the closed-form quantiles match SSJ exactly") {
    const std::vector<double> u = five_uniforms();

    // ExponentialGen(stream, 2.0)
    const std::vector<double> exp_expected = {0.012154767488118706, 0.063972018301778890,
                                              0.089224158962993170, 2.0670789262834948,
                                              0.032452964419235480};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::exponential_inverse(2.0, u[i]) == exp_expected[i]);

    // UniformGen(stream, 1.0, 4.0)
    const std::vector<double> uni_expected = {1.0720493141529750, 1.3602921396821657,
                                              1.4902981261215198, 3.9519515577717510,
                                              1.1885331443545628};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::uniform_inverse(1.0, 4.0, u[i]) == uni_expected[i]);

    // ParetoGen(stream, 2.5, 1.0)
    const std::vector<double> par_expected = {1.0097712438782378, 1.0525098178363330,
                                              1.0739885416450412, 5.2260887065025345,
                                              1.0263023295659728};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::pareto_inverse(2.5, 1.0, u[i]) == par_expected[i]);

    // WeibullGen(stream, 1.5, 2.0, 0.0)
    const std::vector<double> wei_expected = {0.041958611219593720, 0.12695506631201467,
                                              0.15848140402631666, 1.2879372605973654,
                                              0.080753170763643540};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::weibull_inverse(1.5, 2.0, 0.0, u[i]) == wei_expected[i]);

    // BernoulliGen(stream, 0.4)
    const std::vector<double> ber_expected = {0.0, 0.0, 0.0, 1.0, 0.0};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::bernoulli_inverse(0.4, u[i]) == ber_expected[i]);
}

TEST_CASE("ldes ssj variates: the special-function quantiles match SSJ to 1e-13") {
    const std::vector<double> u = five_uniforms();

    // LognormalGen(stream, 0.5, 0.75)
    const std::vector<double> logn_expected = {0.37425867033054850, 0.68326322543066220,
                                               0.79030573277788830, 8.2318069000492230,
                                               0.52283218172642640};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::lognormal_inverse(0.5, 0.75, u[i]) ==
              doctest::Approx(logn_expected[i]).epsilon(1e-13));

    // GammaGen(stream, 2.5, 1.5)
    const std::vector<double> gam_expected = {0.27209226428451766, 0.59016244857977200,
                                              0.69664221131633850, 4.6455870286705350,
                                              0.42607717230478980};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::gamma_inverse(2.5, 1.5, u[i]) == doctest::Approx(gam_expected[i]).epsilon(1e-13));

    // ErlangGen(stream, 3, 2.0): the integer-shape Gamma, ONE uniform per draw
    const std::vector<double> erl_expected = {0.30451359118916360, 0.59904134553448140,
                                              0.69363163180164900, 3.9021559669903790,
                                              0.44996188091528666};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::erlang_inverse(3, 2.0, u[i]) == doctest::Approx(erl_expected[i]).epsilon(1e-13));
}

TEST_CASE("ldes ssj variates: the discrete quantiles match SSJ exactly") {
    const std::vector<double> u = five_uniforms();

    // PoissonGen(stream, 3.7)
    const std::vector<double> poi_expected = {0.0, 2.0, 2.0, 8.0, 1.0};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::poisson_inverse(3.7, u[i]) == poi_expected[i]);

    // BinomialGen(stream, 10, 0.3)
    const std::vector<double> bin_expected = {0.0, 1.0, 2.0, 6.0, 1.0};
    for (std::size_t i = 0; i < u.size(); ++i)
        CHECK(v::binomial_inverse(10, 0.3, u[i]) == bin_expected[i]);
}

TEST_CASE("ldes ssj variates: one uniform per draw, which is what keeps the stream in step") {
    // Each family is drawn five times from its own copy of the stream; the
    // stream must then be exactly five draws along, matching the measurement
    // taken against SSJ. A family that consumed two would shift every later
    // event in the run.
    rng::Mrg32k3a probe = engine_stream();
    for (int i = 0; i < 5; ++i) probe.next_double();
    const double after_five = probe.next_double();

    struct Case {
        const char* name;
        double (*draw)(double);
    };
    const std::vector<Case> cases = {
        {"exponential", [](double x) { return v::exponential_inverse(2.0, x); }},
        {"erlang", [](double x) { return v::erlang_inverse(3, 2.0, x); }},
        {"uniform", [](double x) { return v::uniform_inverse(1.0, 4.0, x); }},
        {"weibull", [](double x) { return v::weibull_inverse(1.5, 2.0, 0.0, x); }},
        {"pareto", [](double x) { return v::pareto_inverse(2.5, 1.0, x); }},
        {"lognormal", [](double x) { return v::lognormal_inverse(0.5, 0.75, x); }},
        {"gamma", [](double x) { return v::gamma_inverse(2.5, 1.5, x); }},
        {"poisson", [](double x) { return v::poisson_inverse(3.7, x); }},
        {"binomial", [](double x) { return v::binomial_inverse(10, 0.3, x); }},
        {"bernoulli", [](double x) { return v::bernoulli_inverse(0.4, x); }},
    };
    for (std::size_t c = 0; c < cases.size(); ++c) {
        rng::Mrg32k3a s = engine_stream();
        for (int i = 0; i < 5; ++i) cases[c].draw(s.next_double());
        CHECK_MESSAGE(s.next_double() == after_five,
                      "family consumed the wrong number of uniforms: ", cases[c].name);
    }
}

TEST_CASE("ldes ssj variates: the quantiles are monotone and land inside their support") {
    rng::Mrg32k3a s = engine_stream();
    double prev_exp = -1.0;
    for (int i = 1; i < 100; ++i) {
        const double u = static_cast<double>(i) / 100.0;
        const double x = v::exponential_inverse(2.0, u);
        CHECK(x > prev_exp);
        prev_exp = x;
        CHECK(v::pareto_inverse(2.5, 1.0, u) >= 1.0);
        CHECK(v::weibull_inverse(1.5, 2.0, 0.5, u) >= 0.5);
        const double g = v::gamma_inverse(2.5, 1.5, u);
        CHECK(g > 0.0);
        CHECK(std::isfinite(g));
    }
    (void)s;
    CHECK_THROWS_AS(v::gamma_inverse(-1.0, 1.0, 0.5), InputError);
    CHECK_THROWS_AS(v::exponential_inverse(0.0, 0.5), InputError);
}
