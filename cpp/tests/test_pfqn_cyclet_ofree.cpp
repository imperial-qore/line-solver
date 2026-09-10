/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Exact passage-time law along an overtake-free path of a closed tree-like
 * product-form network.
 *
 * THE ORACLE IS THE PAPER ITSELF. Harrison and Knottenbelt 2002 publishes the
 * first three moments for four models, computed on their own passage-time
 * analyser over a 1.2-million-state chain; those numbers are asserted here
 * directly. The exact (Theorem 2) and transform (Theorem 1) routes are then
 * checked against each other, which they share no code with beyond the Buzen
 * convolution, and the density is checked to integrate to one.
 *
 * The MATLAB and native Python twins report the same numbers.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_cyclet_ofree.h"

namespace {
namespace pfqn = line::pfqn;

std::vector<double> grid(double a, double b, double step) {
    std::vector<double> t;
    for (double x = a; x <= b + 1e-9; x += step) t.push_back(x);
    return t;
}

}  // namespace

TEST_CASE("the tree network of Fig. 6 reproduces the paper's published moments") {
    const std::vector<double> mu{3.0, 5.0, 4.0, 6.0, 2.0, 1.0};
    const double p12 = 0.2, p13 = 0.5, p14 = 0.3;
    const std::vector<double> v{1.0, p12, p13, p14, p12, p14};
    const std::vector<std::vector<std::size_t>> paths{{0, 2}, {0, 1, 4}, {0, 3, 5}};
    const std::vector<double> pp{p13, p12, p14};
    const std::vector<double> t = grid(0.0, 40.0, 0.25);

    const pfqn::CycletResult r =
        pfqn::pfqn_cyclet_ofree(v, mu, 18, paths, t, "auto", 3, pp);

    // Sec. 7.2: "6.12717, 53.3067 and 612.887".
    CHECK(r.mom[0] == doctest::Approx(6.12717).epsilon(1e-6));
    CHECK(r.mom[1] == doctest::Approx(53.3067).epsilon(1e-6));
    CHECK(r.mom[2] == doctest::Approx(612.887).epsilon(1e-6));
    for (const pfqn::CycletPathInfo& i : r.info) CHECK(i.method == "exact");

    // The density integrates to one, and its first moment is the one above.
    double m0 = 0.0, m1 = 0.0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        const double dt = t[i] - t[i - 1];
        m0 += 0.5 * dt * (r.f[i] + r.f[i - 1]);
        m1 += 0.5 * dt * (t[i] * r.f[i] + t[i - 1] * r.f[i - 1]);
    }
    CHECK(m0 == doctest::Approx(1.0).epsilon(1e-3));
    CHECK(m1 == doctest::Approx(r.mom[0]).epsilon(1e-4));
    CHECK(r.F.back() == doctest::Approx(1.0).epsilon(1e-5));
}

TEST_CASE("the Erlang and branching-Erlang examples of Sec. 6") {
    // Sec. 6.1: a 3-stage Erlang with lambda = 2, moments 1.5, 3 and 7.5. Equal
    // rates, so 'auto' MUST leave Theorem 2 alone: its partial fractions divide
    // by the rate differences.
    const pfqn::CycletResult e3 = pfqn::pfqn_cyclet_ofree(
        {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 1, {{0, 1, 2}}, {1.0});
    CHECK(e3.mom[0] == doctest::Approx(1.5).epsilon(1e-10));
    CHECK(e3.mom[1] == doctest::Approx(3.0).epsilon(1e-10));
    CHECK(e3.mom[2] == doctest::Approx(7.5).epsilon(1e-10));
    CHECK(e3.info[0].method == "lt");

    // Sec. 6.2: a 3-stage Erlang(1) branch and a 12-stage Erlang(2) branch, each
    // with probability 1/2. Moments 4.5, 25.5 and 166.5.
    std::vector<double> vb, mub;
    for (std::size_t i = 0; i < 3; ++i) {
        vb.push_back(0.5);
        mub.push_back(1.0);
    }
    for (std::size_t i = 0; i < 12; ++i) {
        vb.push_back(0.5);
        mub.push_back(2.0);
    }
    std::vector<std::size_t> lower, upper;
    for (std::size_t i = 0; i < 3; ++i) lower.push_back(i);
    for (std::size_t i = 3; i < 15; ++i) upper.push_back(i);
    const pfqn::CycletResult br =
        pfqn::pfqn_cyclet_ofree(vb, mub, 1, {lower, upper}, {1.0}, "auto", 3, {0.5, 0.5});
    CHECK(br.mom[0] == doctest::Approx(4.5).epsilon(1e-10));
    CHECK(br.mom[1] == doctest::Approx(25.5).epsilon(1e-10));
    CHECK(br.mom[2] == doctest::Approx(166.5).epsilon(1e-10));

    // Sec. 6.2, the first-to-last passage of the upper branch: an Erlang-11,
    // moments 5.5, 33 and 214.5.
    std::vector<std::size_t> u11;
    for (std::size_t i = 3; i < 14; ++i) u11.push_back(i);
    const pfqn::CycletResult up = pfqn::pfqn_cyclet_ofree(vb, mub, 1, {u11}, {1.0});
    CHECK(up.mom[0] == doctest::Approx(5.5).epsilon(1e-10));
    CHECK(up.mom[1] == doctest::Approx(33.0).epsilon(1e-10));
    CHECK(up.mom[2] == doctest::Approx(214.5).epsilon(1e-10));
}

TEST_CASE("Theorem 2 and Theorem 1 agree where both apply") {
    const std::vector<double> mu{3.0, 5.0, 4.0, 6.0, 2.0, 1.0};
    const std::vector<double> v{1.0, 0.2, 0.5, 0.3, 0.2, 0.3};
    const std::vector<std::vector<std::size_t>> z{{0, 1, 4}};
    const std::vector<double> t = grid(1.0, 15.0, 0.5);
    const pfqn::CycletResult a = pfqn::pfqn_cyclet_ofree(v, mu, 18, z, t, "exact");
    const pfqn::CycletResult b = pfqn::pfqn_cyclet_ofree(v, mu, 18, z, t, "lt");
    for (std::size_t i = 0; i < t.size(); ++i) {
        CHECK(b.f[i] == doctest::Approx(a.f[i]).epsilon(1e-6));
        CHECK(b.F[i] == doctest::Approx(a.F[i]).epsilon(1e-6));
    }
    for (std::size_t q = 0; q < 3; ++q) CHECK(a.mom[q] == doctest::Approx(b.mom[q]));
}

TEST_CASE("the refusals are by name") {
    const std::vector<double> v{1.0, 0.5}, mu{2.0, 3.0};
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, {2.0}, 3, {{0}}, {1.0}), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, {2.0, -1.0}, 3, {{0}}, {1.0}), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, mu, 3, {{0, 9}}, {1.0}), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, mu, 3, {{0, 0}}, {1.0}), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, mu, 3, {{0}}, {1.0}, "nope"), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_cyclet_ofree(v, mu, 3, {{0}}, {1.0}, "auto", 3, {0.5, 0.5}),
                    line::InputError);
    // Equal rates on the path are refused by the exact route by name, not
    // silently divided by zero.
    CHECK_THROWS_AS(
        pfqn::pfqn_cyclet_ofree({1.0, 1.0}, {2.0, 2.0}, 2, {{0, 1}}, {1.0}, "exact"),
        line::InputError);
}
