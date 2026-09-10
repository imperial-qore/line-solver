/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Cost-capped cache normalizing constant, per-item occupancy and mean per-list
 * storage cost. Oracles, in order of strength:
 *   1. The worked example of Casale-Gast, IEEE/ACM Trans. Networking 29(2),
 *      2021, Sec. IX (Model 7): E(m) = 238.7982 and E(m,k) = 7.7963 for the
 *      Fig. 2 tree with n = 10, m = (2,1,1,2), sizes 1 and 2, caps (2,1,2,4).
 *   2. Reduction: unit sizes with caps equal to the capacities must reproduce
 *      the unconstrained constant exactly.
 *   3. Conservation: sum_i pi_ij = m_j still holds under the caps, and an item
 *      too large for a list has probability exactly zero of sitting in it --
 *      both asserted with == in the exact (rational) instantiation.
 *   4. Agreement between the importance sampler and the exact recursion.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_cost.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_is.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_spm_size.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/nc/solver_nc_cache.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::cache::cache_cost;
using line::cache::cache_cost_pathcheck;
using line::cache::cache_erec;
using line::cache::cache_is;
using line::cache::cache_prob_erec;
using line::cache::CacheCostMode;
using line::cache::cache_spm_size;

namespace {

/**
 * Fig. 2 tree of the paper: p(1) = p(2) = miss, p(3) = p(4) = list 2. Stream 1
 * requests items 1..5 at rate 0.9 with c0 = 0.025 and c = c' = 0; stream 2
 * requests items 6..10 at rate 1.0 with c0 = 1.0 and c = c' = 0.5.
 */
template <class T>
Matrix<T> paper_gamma() {
    Matrix<T> g(10, 4, line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < 5; ++i) {
        g(i, 0) = line::num_traits<T>::from_double(0.9 * (1 - 0.025));
        g(i, 1) = line::num_traits<T>::from_double(0.9 * 0.025);
    }
    for (std::size_t i = 5; i < 10; ++i) {
        g(i, 1) = line::num_traits<T>::from_int(1);
        g(i, 2) = line::num_traits<T>::from_double(0.5);
        g(i, 3) = line::num_traits<T>::from_double(0.5);
    }
    return g;
}

std::vector<int> paper_m() { return std::vector<int>{2, 1, 1, 2}; }
std::vector<int> paper_sigma() { return std::vector<int>{1, 1, 1, 1, 1, 2, 2, 2, 2, 2}; }
std::vector<int> paper_k() { return std::vector<int>{2, 1, 2, 4}; }

}  // namespace

TEST_CASE("cache_erec reproduces the paper's constrained normalizing constant") {
    const Matrix<double> g = paper_gamma<double>();
    CHECK(cache_erec(g, paper_m()) == doctest::Approx(238.7982).epsilon(1e-6));
    CHECK(cache_erec(g, paper_m(), paper_sigma(), paper_k()) ==
          doctest::Approx(7.7963).epsilon(1e-5));
}

TEST_CASE("unit sizes with saturated caps reduce to the unconstrained constant") {
    const Matrix<Rational> g = paper_gamma<Rational>();
    const std::vector<int> unit(10, 1);
    CHECK(cache_erec(g, paper_m(), unit, paper_m()) == cache_erec(g, paper_m()));
}

TEST_CASE("constrained marginals fill every list and exclude oversized items") {
    const Matrix<Rational> g = paper_gamma<Rational>();
    const Matrix<Rational> pij = cache_prob_erec(g, paper_m(), paper_sigma(), paper_k());
    const std::vector<int> m = paper_m();
    for (std::size_t j = 0; j < 4; ++j) {
        Rational col = line::num_traits<Rational>::from_int(0);
        for (std::size_t i = 0; i < 10; ++i) col += pij(i, j + 1);
        CHECK(col == line::num_traits<Rational>::from_int(m[j]));
    }
    // list 2 has cap 1, so the size-2 items cannot sit in it
    for (std::size_t i = 5; i < 10; ++i)
        CHECK(pij(i, 2) == line::num_traits<Rational>::from_int(0));
    // in this instance every cap is saturated
    const std::vector<Rational> K =
        cache_cost(g, paper_m(), paper_sigma(), paper_k(), pij);
    const std::vector<int> k = paper_k();
    for (std::size_t j = 0; j < 4; ++j)
        CHECK(K[j] == line::num_traits<Rational>::from_int(k[j]));
}

TEST_CASE("the importance sampler agrees with the exact constrained constant") {
    // small linear cache where every size-feasible state is reachable
    Matrix<double> g(6, 2);
    for (std::size_t i = 0; i < 6; ++i) {
        g(i, 0) = 1.0 / 3.0;
        g(i, 1) = 1.0 / 9.0;
    }
    const std::vector<int> m{1, 1};
    const std::vector<int> sigma{1, 1, 1, 2, 2, 2};
    const std::vector<int> k{2, 1};
    const double exact = cache_erec(g, m, sigma, k);
    const double sampled = cache_is(g, m, static_cast<std::size_t>(200000),
                                    static_cast<std::uint64_t>(7), sigma, k)
                               .E;
    CHECK(sampled == doctest::Approx(exact).epsilon(0.05));
}

TEST_CASE("the path check flags the blocked promotion of the paper instance") {
    const Matrix<double> g = paper_gamma<double>();
    const std::vector<int> parent{-1, -1, 1, 1};
    const std::vector<line::cache::CacheBlockedPair> viol =
        cache_cost_pathcheck(g, paper_sigma(), paper_k(), parent);
    // the size-2 items may sit in lists 3 and 4 but can never traverse list 2
    CHECK(!viol.empty());
    for (std::size_t v = 0; v < viol.size(); ++v) {
        CHECK(viol[v].item >= 5);
        CHECK((viol[v].list == 2 || viol[v].list == 3));
        CHECK(viol[v].blocking_list == 1);
    }
    // the same cache without a cap on the intermediate list is clean
    const std::vector<int> kopen{2, 2, 2, 4};
    CHECK(cache_cost_pathcheck(g, paper_sigma(), kopen, parent).empty());
}

namespace {

/**
 * Two-list RR cache over six items of sizes 1 and 2, read uniformly. The list
 * capacities and the cost caps are the knobs: together they decide whether the
 * SPM family can carry the caps itself or the analyzer has to leave it.
 */
line::nc::NcCacheSolution<double> costcap_cache(const std::vector<int>& itemcap,
                                                const std::vector<int>& costcap) {
    line::qn::Network<double> m("cacheCostCap");
    const std::size_t src = m.add_source("Source");
    line::qn::CacheParam<double> ch;
    ch.nitems = 6;
    ch.itemcap = itemcap;
    ch.itemsize = std::vector<int>{1, 1, 1, 2, 2, 2};
    ch.costcap = costcap;
    ch.replacestrat = line::lang::ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double> >{
        std::vector<double>(6, 1.0 / 6.0), std::vector<double>(), std::vector<double>()};
    ch.hitclass = std::vector<std::size_t>{2u, 0, 0};
    ch.missclass = std::vector<std::size_t>{3u, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t job = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t mis = m.add_open_class("MissClass");
    m.set_arrival(src, job, line::lang::Distrib<double>::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(job, job, src, cn, 1.0);
    P.set(hit, hit, cn, snk, 1.0);
    P.set(mis, mis, cn, snk, 1.0);
    m.link(P);

    line::nc::NcSolverOptions opt;
    opt.method = "default";  // the SPM branch
    return line::nc::solver_nc_cache_analyzer(m.get_struct(), opt);
}

}  // namespace

/**
 * Sizes and caps with room under the saddle -- sum(m) = 2 against six items --
 * so the SPM family carries the caps ITSELF, through its size-tilted form
 * cache_spm_size, and there is nothing to switch away from. List 2 has cap 1,
 * so it admits the size-1 items only.
 */
TEST_CASE("the size-tilted SPM carries the cost caps with no method switch") {
    const line::nc::NcCacheSolution<double> r = costcap_cache(std::vector<int>{1, 1},
                                                              std::vector<int>{2, 1});

    CHECK_FALSE(r.costcap_method_switched);
    CHECK(r.sol.warning.find("does not support storage cost caps") == std::string::npos);
    CHECK(r.sol.actualmethod == "spm.size");
    // the per-list cost survives on the inner solution the runner returns
    CHECK(r.listcost.size() == 2);
    CHECK(r.sol.listcost.size() == 2);
    CHECK(line::num_traits<double>::to_double(r.listcost[0]) == doctest::Approx(1.6).epsilon(1e-9));
    CHECK(line::num_traits<double>::to_double(r.listcost[1]) == doctest::Approx(1.0).epsilon(1e-9));
    // caps honoured: the size-2 items never enter list 2. The saddle point is an
    // expansion, so it returns that zero to its own accuracy rather than the hard
    // zero the exact recursion below returns.
    for (std::size_t i = 3; i < 6; ++i) CHECK(std::abs(r.pij(i, 2)) < 1e-9);
}

/**
 * The saddle escapes to infinity at sum(m) = n: with the cache holding every
 * item the size-tilted expansion has no room, so the analyzer leaves the SPM
 * family for the exact recursion and SAYS SO on the solution's warning channel
 * rather than only setting a flag.
 *
 * The caps bind exactly here, which fixes the occupancy by hand. Total cost is
 * 3*1 + 3*2 = 9, so caps (4,5) admit only the partitions costing 4 in list 1
 * and 5 in list 2 -- two small items and one large against one small and two
 * large. Reads are uniform, so the nine such partitions are equally likely and
 * each small item sits in list 2 in three of them.
 */
TEST_CASE("the cost-cap conditions reach the solution's warning channel") {
    const line::nc::NcCacheSolution<double> r = costcap_cache(std::vector<int>{3, 3},
                                                              std::vector<int>{4, 5});

    // the switch is REPORTED, not merely flagged
    CHECK(r.costcap_method_switched);
    CHECK(r.sol.warning.find("does not support storage cost caps") != std::string::npos);
    CHECK(r.sol.actualmethod == "exact");
    // no slack is left under either cap, so the mean cost is the cap itself
    REQUIRE(r.listcost.size() == 2);
    REQUIRE(r.sol.listcost.size() == 2);
    CHECK(line::num_traits<double>::to_double(r.listcost[0]) == doctest::Approx(4.0).epsilon(1e-9));
    CHECK(line::num_traits<double>::to_double(r.listcost[1]) == doctest::Approx(5.0).epsilon(1e-9));
    // every item is held, and the exact recursion places it by size
    for (std::size_t i = 0; i < 6; ++i) {
        CAPTURE(i);
        const double inlist1 = i < 3 ? 2.0 / 3.0 : 1.0 / 3.0;
        CHECK(std::abs(r.pij(i, 0)) < 1e-12);
        CHECK(r.pij(i, 1) == doctest::Approx(inlist1).epsilon(1e-9));
        CHECK(r.pij(i, 2) == doctest::Approx(1.0 - inlist1).epsilon(1e-9));
    }
}

/**
 * Ray (WKB) expansion of the cost-capped constant, cache_spm_size. The
 * profile is the one the derivation note validates against its own dynamic
 * program: gamma_i = 0.3 + 2.7 i/n, sizes {1,2,3} in equal thirds of the item
 * index (so gcd = 1 and the sizes are genuinely diverse), m = n/4, k = 1.8 m.
 */
TEST_CASE("cache_spm_size: exact-cost mode reproduces the note's table") {
    // log H with the cost resolved exactly at k, i.e. the raw Laplace formula.
    const int ns[4] = {50, 100, 200, 400};
    const int ms[4] = {12, 25, 50, 100};
    const double want[4] = {27.234323, 57.620291, 118.438233, 240.843731};
    for (int t = 0; t < 4; ++t) {
        const int n = ns[t], m = ms[t], k = static_cast<int>(std::ceil(1.8 * m));
        Matrix<double> gamma(static_cast<std::size_t>(n), 1);
        std::vector<int> sigma(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            gamma(static_cast<std::size_t>(i), 0) = 0.3 + 2.7 * ((i + 1.0) / n);
            sigma[static_cast<std::size_t>(i)] = 1 + (3 * i) / n;
        }
        const auto r = cache_spm_size<double>(gamma, {m}, sigma, {k}, CacheCostMode::Exact);
        double logfact = 0.0;
        for (int i = 2; i <= m; ++i) logfact += std::log(static_cast<double>(i));
        CHECK(r.log_e - logfact == doctest::Approx(want[t]).epsilon(1e-7));
    }
}

TEST_CASE("cache_spm_size: cumulative caps track cache_erec and degenerate to size-free") {
    const int n = 120;
    Matrix<double> gamma(static_cast<std::size_t>(n), 1);
    std::vector<int> sigma(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        gamma(static_cast<std::size_t>(i), 0) = 0.3 + 2.7 * ((i + 1.0) / n);
        sigma[static_cast<std::size_t>(i)] = 1 + (3 * i) / n;
    }
    const std::vector<int> m{30};

    // binding cap: the expansion sits close to the exact recursion
    {
        const auto r = cache_spm_size<double>(gamma, m, sigma, {54});
        const double exact = std::log(cache_erec<double>(gamma, m, sigma, {54}));
        CHECK(r.binding[0]);
        CHECK(r.zeta[0] < 1.0);
        CHECK(std::abs(r.log_e - exact) < 0.15);
        CHECK(r.method == "spm-size");
        CHECK(r.iterations < 30);
    }
    // slack cap: zeta returns to 1, the cost coordinate leaves the saddle
    {
        const auto r = cache_spm_size<double>(gamma, m, sigma, {200});
        const double exact = std::log(cache_erec<double>(gamma, m, sigma, {200}));
        CHECK_FALSE(r.binding[0]);
        CHECK(r.zeta[0] == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(std::abs(r.log_e - exact) < 0.05);
    }
    // a cap below the cheapest m items admits no state at all
    {
        const auto r = cache_spm_size<double>(gamma, m, sigma, {10});
        CHECK(r.e == 0.0);
        CHECK(r.method == "boundary");
    }
}

TEST_CASE("cache_spm_size: the size lattice is divided out, uniform sizes are exact") {
    const int n = 120;
    Matrix<double> gamma(static_cast<std::size_t>(n), 1);
    std::vector<int> sigma(static_cast<std::size_t>(n)), doubled(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        gamma(static_cast<std::size_t>(i), 0) = 0.3 + 2.7 * ((i + 1.0) / n);
        sigma[static_cast<std::size_t>(i)] = 1 + (3 * i) / n;
        doubled[static_cast<std::size_t>(i)] = 2 * sigma[static_cast<std::size_t>(i)];
    }
    const std::vector<int> m{30};
    const auto base = cache_spm_size<double>(gamma, m, sigma, {54});
    // sizes 2/4/6 with cap 108 is the same problem; 109 snaps down to the same lattice point
    const auto twice = cache_spm_size<double>(gamma, m, doubled, {108});
    const auto odd = cache_spm_size<double>(gamma, m, doubled, {109});
    CHECK(twice.span == 2);
    CHECK(twice.log_e == doctest::Approx(base.log_e).epsilon(1e-12));
    CHECK(odd.log_e == doctest::Approx(base.log_e).epsilon(1e-12));

    // a single item size: the cost of the list is sigma*m identically, so the cap
    // carries no information and the 2h saddle would be singular (note Sec. 5)
    const std::vector<int> flat(static_cast<std::size_t>(n), 2);
    const auto slack = cache_spm_size<double>(gamma, m, flat, {60});
    CHECK(slack.method == "uniform-size");
    CHECK(slack.zeta[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(slack.log_e == doctest::Approx(std::log(cache_erec<double>(gamma, m, flat, {60})))
                             .epsilon(1e-3));
    const auto tight = cache_spm_size<double>(gamma, m, flat, {59});
    CHECK(tight.e == 0.0);   // 2*30 = 60 > 59, no state fits
}

TEST_CASE("cache_spm_size: occupancy and mean cost are consistent") {
    const int n = 80;
    Matrix<double> gamma(static_cast<std::size_t>(n), 1);
    std::vector<int> sigma(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        gamma(static_cast<std::size_t>(i), 0) = 0.3 + 2.7 * ((i + 1.0) / n);
        sigma[static_cast<std::size_t>(i)] = 1 + (3 * i) / n;
    }
    const auto r = cache_spm_size<double>(gamma, {20}, sigma, {36});
    double occ = 0.0, cost = 0.0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
        CHECK(r.pij(i, 1) >= 0.0);
        CHECK(r.pij(i, 0) >= 0.0);
        CHECK(r.pij(i, 0) + r.pij(i, 1) == doctest::Approx(1.0).epsilon(1e-12));
        occ += r.pij(i, 1);
        cost += sigma[i] * r.pij(i, 1);
    }
    // the saddle conditions are exactly sum_i pi_i = m and sum_i sigma_i pi_i = k
    CHECK(occ == doctest::Approx(20.0).epsilon(1e-9));
    CHECK(cost == doctest::Approx(36.0).epsilon(1e-9));
    CHECK(r.k_mean[0] == doctest::Approx(cost).epsilon(1e-12));
}
