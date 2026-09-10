/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * TTL cache approximations that need a matrix exponential and a scalar root
 * finder: cache_lrum_map_levelstats, cache_t_lrum_map, cache_ttl_lrum_map
 * (LRU(m) with MAP request streams) and cache_ttl_lrua (LRU over an access
 * graph).
 *
 * Oracles:
 *   - the conservation laws already used across test_cache_*.cpp: every item's
 *     level probabilities sum to one, and each list column sums to its
 *     capacity (that second law IS the equation the characteristic times
 *     solve, so it is the residual of the fixed point);
 *   - the h = 1 access graph of cache_ttl_lrua reduces analytically to the Che
 *     approximation, prob(k,1) = 1 - exp(-lambda_k T), which is checked in
 *     closed form and against cache_t_hlru, an independently ported function;
 *   - a MAP degenerating to a Poisson process makes the LRU(m)-MAP model
 *     coincide with the Poisson TTL model of cache_ttl_hlru;
 *   - MATLAB reference values (R2025a). MATLAB solves both fixed points with
 *     fsolve and stops well short of convergence, leaving capacity residuals
 *     of 1.4e-7 (cache_ttl_lrua) and 1e-10 (cache_t_lrum_map); the comparisons
 *     are asserted at the accuracy MATLAB delivers, and the conservation laws
 *     at the far tighter accuracy this port delivers.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_lrum_map_levelstats.h"
#include "line/api/cache/cache_t_hlru.h"
#include "line/api/cache/cache_t_lrum_map.h"
#include "line/api/cache/cache_ttl_hlru.h"
#include "line/api/cache/cache_ttl_lrua.h"
#include "line/api/cache/cache_ttl_lrum_map.h"

using line::Matrix;
using line::mam::Map;
using namespace line::cache;

namespace {

constexpr double MATLAB_TOL = 1e-9;
constexpr double FSOLVE_TOL = 1e-5;  // what MATLAB's fsolve actually delivers

/** Item k as a two-phase MMPP with the given per-phase request rates. */
Map<double> item(double r0, double r1) {
    Map<double> m;
    const double a = 0.1, b = 0.2;
    m.D0 = Matrix<double>{{-(r0 + a), a}, {b, -(r1 + b)}};
    m.D1 = Matrix<double>{{r0, 0.0}, {0.0, r1}};
    return m;
}

std::vector<Map<double>> catalogue() {
    std::vector<Map<double>> v;
    v.push_back(item(0.5, 1.5));
    v.push_back(item(0.3, 0.9));
    v.push_back(item(0.8, 2.4));
    return v;
}

/** Linear access graph: list l promotes to l+1, the last list self-loops. */
std::vector<Matrix<double>> linear_graph(std::size_t n, std::size_t h) {
    std::vector<Matrix<double>> R;
    for (std::size_t k = 0; k < n; ++k) {
        Matrix<double> Rk(h + 1, h + 1, 0.0);
        for (std::size_t l = 0; l < h; ++l) Rk(l, l + 1) = 1.0;
        Rk(h, h) = 1.0;
        R.push_back(Rk);
    }
    return R;
}

}  // namespace

TEST_CASE("cache_lrum_map_levelstats matches MATLAB and conserves probability") {
    const std::vector<Map<double>> items = catalogue();
    const std::vector<double> T{1.5, 3.0};

    const CacheLrumMapLevelStats<double> s1 =
        cache_lrum_map_levelstats(items[0].D0, items[0].D1, T);
    CHECK(s1.prob[0] == doctest::Approx(0.082111388651669365).epsilon(MATLAB_TOL));
    CHECK(s1.prob[1] == doctest::Approx(0.11836826223776498).epsilon(MATLAB_TOL));
    CHECK(s1.prob[2] == doctest::Approx(0.79952034911056569).epsilon(MATLAB_TOL));
    CHECK(s1.occ[0] == s1.prob[1]);
    CHECK(s1.occ[1] == s1.prob[2]);
    CHECK(s1.hitfrac[0] == doctest::Approx(0.088635874825935546).epsilon(MATLAB_TOL));
    CHECK(s1.hitfrac[1] == doctest::Approx(0.85351712014183467).epsilon(MATLAB_TOL));

    const CacheLrumMapLevelStats<double> s3 =
        cache_lrum_map_levelstats(items[2].D0, items[2].D1, T);
    CHECK(s3.prob[0] == doctest::Approx(0.017567910076823012).epsilon(MATLAB_TOL));
    CHECK(s3.prob[2] == doctest::Approx(0.9328240721267177).epsilon(MATLAB_TOL));
    CHECK(s3.hitfrac[0] == doctest::Approx(0.03385411452377915).epsilon(MATLAB_TOL));
    CHECK(s3.hitfrac[1] == doctest::Approx(0.95442519424103645).epsilon(MATLAB_TOL));

    // The level probabilities are a distribution, and the hit fractions of an
    // item never exceed one (they are a fraction of its own requests).
    for (const CacheLrumMapLevelStats<double>* s : {&s1, &s3}) {
        double tot = 0.0, hits = 0.0;
        for (std::size_t l = 0; l < s->prob.size(); ++l) tot += s->prob[l];
        for (std::size_t l = 0; l < s->hitfrac.size(); ++l) {
            CHECK(s->hitfrac[l] >= 0.0);
            hits += s->hitfrac[l];
        }
        CHECK(tot == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(hits <= 1.0 + 1e-12);
    }

    SUBCASE("single list") {
        const CacheLrumMapLevelStats<double> s =
            cache_lrum_map_levelstats(items[1].D0, items[1].D1, std::vector<double>{2.0});
        CHECK(s.prob[0] == doctest::Approx(0.41227081726012982).epsilon(MATLAB_TOL));
        CHECK(s.prob[1] == doctest::Approx(0.58772918273987018).epsilon(MATLAB_TOL));
        CHECK(s.hitfrac[0] == doctest::Approx(0.66610995533321948).epsilon(MATLAB_TOL));
    }

    SUBCASE("a Poisson item reduces to the Che approximation") {
        // For a one-phase MAP (a Poisson process of rate lam) and one list, the
        // embedded chain is the two-state Che chain and the time-stationary
        // probability of being cached is 1 - exp(-lam T).
        const double lam = 0.75, T1 = 1.6;
        Map<double> p;
        p.D0 = Matrix<double>(1, 1, -lam);
        p.D1 = Matrix<double>(1, 1, lam);
        const CacheLrumMapLevelStats<double> s =
            cache_lrum_map_levelstats(p.D0, p.D1, std::vector<double>{T1});
        CHECK(s.prob[1] == doctest::Approx(1.0 - std::exp(-lam * T1)).epsilon(1e-12));
        CHECK(s.hitfrac[0] == doctest::Approx(1.0 - std::exp(-lam * T1)).epsilon(1e-12));
    }

    SUBCASE("malformed input is rejected") {
        CHECK_THROWS_AS(cache_lrum_map_levelstats(items[0].D0, items[0].D1, std::vector<double>{}),
                        line::InputError);
        CHECK_THROWS_AS(
            cache_lrum_map_levelstats(items[0].D0, items[0].D1, std::vector<double>{-1.0}),
            line::InputError);
    }
}

TEST_CASE("cache_t_lrum_map solves the capacity equations to full accuracy") {
    const std::vector<Map<double>> items = catalogue();
    const std::vector<double> m{1.0, 1.0};
    const std::vector<double> t = cache_t_lrum_map(items, m, 1e-13);

    // MATLAB (fsolve on log T): [1.1271170497683802, 0.75964586208290563]
    CHECK(t[0] == doctest::Approx(1.1271170497683802).epsilon(1e-8));
    CHECK(t[1] == doctest::Approx(0.75964586208290563).epsilon(1e-8));

    // The defining equations, at this port's accuracy rather than fsolve's.
    std::vector<double> occ(2, 0.0);
    for (std::size_t k = 0; k < items.size(); ++k) {
        const CacheLrumMapLevelStats<double> s =
            cache_lrum_map_levelstats(items[k].D0, items[k].D1, t);
        occ[0] += s.occ[0];
        occ[1] += s.occ[1];
    }
    CHECK(occ[0] == doctest::Approx(1.0).epsilon(1e-11));
    CHECK(occ[1] == doctest::Approx(1.0).epsilon(1e-11));

    SUBCASE("a cache that cannot be smaller than the catalogue is rejected") {
        CHECK_THROWS_AS(cache_t_lrum_map(items, std::vector<double>{2.0, 1.0}, 1e-12),
                        line::InputError);
        CHECK_THROWS_AS(cache_t_lrum_map(items, std::vector<double>{0.0}, 1e-12), line::InputError);
    }
}

TEST_CASE("cache_ttl_lrum_map matches MATLAB and conserves probability and capacity") {
    const std::vector<Map<double>> items = catalogue();
    const std::vector<double> m{1.0, 1.0};
    const CacheTtlLrumMapResult<double> r = cache_ttl_lrum_map(items, m, 1e-13);

    // MATLAB pij, column major: miss, hit in list 1, hit in list 2.
    const double miss[3] = {0.2402243348684272, 0.41996019060673462, 0.11504819831016455};
    const double hit1[3] = {0.33531417467700525, 0.34350740933567808, 0.26769043026046363};
    const double hit2[3] = {0.42446149045456755, 0.23653240005758727, 0.61726137142937176};
    const double tim0[3] = {0.32261540005148792, 0.50894059082371434, 0.1684440090868376};
    const double tim1[3] = {0.3537543460142184, 0.32070315752119505, 0.32554249654847456};
    const double tim2[3] = {0.32363025393429368, 0.17035625165509058, 0.5060134943646879};
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.pij(k, 0) == doctest::Approx(miss[k]).epsilon(1e-8));
        CHECK(r.pij(k, 1) == doctest::Approx(hit1[k]).epsilon(1e-8));
        CHECK(r.pij(k, 2) == doctest::Approx(hit2[k]).epsilon(1e-8));
        CHECK(r.pijtime(k, 0) == doctest::Approx(tim0[k]).epsilon(1e-8));
        CHECK(r.pijtime(k, 1) == doctest::Approx(tim1[k]).epsilon(1e-8));
        CHECK(r.pijtime(k, 2) == doctest::Approx(tim2[k]).epsilon(1e-8));
    }

    // Conservation: each row of both matrices is a distribution over the
    // levels, and each list column of the time-stationary matrix sums to the
    // capacity of that list.
    double occ1 = 0.0, occ2 = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        double rowp = 0.0, rowt = 0.0;
        for (std::size_t l = 0; l < 3; ++l) {
            rowp += r.pij(k, l);
            rowt += r.pijtime(k, l);
            CHECK(r.pij(k, l) >= 0.0);
            CHECK(r.pijtime(k, l) >= 0.0);
        }
        CHECK(rowp == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(rowt == doctest::Approx(1.0).epsilon(1e-14));
        occ1 += r.pijtime(k, 1);
        occ2 += r.pijtime(k, 2);
    }
    CHECK(occ1 == doctest::Approx(m[0]).epsilon(1e-11));
    CHECK(occ2 == doctest::Approx(m[1]).epsilon(1e-11));

    // The more popular the item, the more of its requests hit.
    CHECK(r.pij(2, 0) < r.pij(0, 0));  // item 3 is the most popular
    CHECK(r.pij(1, 0) > r.pij(0, 0));  // item 2 the least
}

TEST_CASE("cache_ttl_lrum_map with Poisson items agrees with cache_ttl_hlru") {
    // One list, four items whose MAPs are Poisson processes: both models then
    // reduce to the same Che fixed point, though they were ported from
    // different MATLAB files and solve different equations.
    const double lam[4] = {0.4, 0.3, 0.2, 0.1};
    std::vector<Map<double>> items;
    Matrix<double> rates(1, 4);
    for (std::size_t k = 0; k < 4; ++k) {
        Map<double> p;
        p.D0 = Matrix<double>(1, 1, -lam[k]);
        p.D1 = Matrix<double>(1, 1, lam[k]);
        items.push_back(p);
        rates(0, k) = lam[k];
    }
    const CacheTtlLrumMapResult<double> r =
        cache_ttl_lrum_map(items, std::vector<double>{2.0}, 1e-13);
    const Matrix<double> p = cache_ttl_hlru(rates, std::vector<int>{2});
    const std::vector<double> t = cache_t_hlru(Matrix<double>(rates.transpose()),
                                               std::vector<int>{2});
    CHECK(r.t[0] == doctest::Approx(t[0]).epsilon(1e-8));
    for (std::size_t k = 0; k < 4; ++k) {
        CHECK(r.pijtime(k, 1) == doctest::Approx(p(k, 1)).epsilon(1e-8));
        CHECK(r.pijtime(k, 1) == doctest::Approx(1.0 - std::exp(-lam[k] * r.t[0])).epsilon(1e-10));
    }
}

TEST_CASE("cache_ttl_lrua on a single list is the Che approximation") {
    // Five items, rates 1.0 .. 0.2, one list of capacity 2, rate unchanged by
    // caching: prob(k,1) = 1 - exp(-lam_k T) and sum_k prob(k,1) = 2.
    const double base[5] = {1.0, 0.8, 0.6, 0.4, 0.2};
    Matrix<double> lam(5, 2);
    for (std::size_t k = 0; k < 5; ++k) {
        lam(k, 0) = base[k];
        lam(k, 1) = base[k];
    }
    const Matrix<double> P =
        cache_ttl_lrua(lam, linear_graph(5, 1), std::vector<double>{2.0}, 1e-13);

    // MATLAB reference (fsolve, residual 5e-8 on the capacity)
    const double ref0[5] = {0.4042687142692426, 0.48454722041868115, 0.58076719896534623,
                            0.69609426115913431, 0.83432263613013302};
    const double ref1[5] = {0.5957312857307574, 0.51545277958131897, 0.41923280103465382,
                            0.30390573884086575, 0.16567736386986698};
    double occ = 0.0;
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(P(k, 0) == doctest::Approx(ref0[k]).epsilon(FSOLVE_TOL));
        CHECK(P(k, 1) == doctest::Approx(ref1[k]).epsilon(FSOLVE_TOL));
        CHECK(P(k, 0) + P(k, 1) == doctest::Approx(1.0).epsilon(1e-14));
        occ += P(k, 1);
    }
    CHECK(occ == doctest::Approx(2.0).epsilon(1e-11));

    // Closed form: the characteristic time recovered from item 1 reproduces
    // every other item, and matches the independently ported cache_t_hlru.
    const double T = -std::log(1.0 - P(0, 1)) / base[0];
    for (std::size_t k = 0; k < 5; ++k)
        CHECK(P(k, 1) == doctest::Approx(1.0 - std::exp(-base[k] * T)).epsilon(1e-10));
    Matrix<double> g(5, 1);
    for (std::size_t k = 0; k < 5; ++k) g(k, 0) = base[k];
    const std::vector<double> t = cache_t_hlru(g, std::vector<int>{2});
    CHECK(T == doctest::Approx(t[0]).epsilon(1e-8));
}

TEST_CASE("cache_ttl_lrua on a two-list access graph matches MATLAB") {
    const double base[5] = {1.0, 0.8, 0.6, 0.4, 0.2};
    Matrix<double> lam(5, 3);
    for (std::size_t k = 0; k < 5; ++k)
        for (std::size_t l = 0; l < 3; ++l) lam(k, l) = base[k];
    const std::vector<double> m{2.0, 1.0};
    const Matrix<double> P = cache_ttl_lrua(lam, linear_graph(5, 2), m, 1e-13);

    const double ref0[5] = {0.15671452658097063, 0.2360139676514984, 0.35255416529257544,
                            0.5176660129215469, 0.73705146500903151};
    const double ref1[5] = {0.46284505334831261, 0.47277469654429549, 0.45173115670574743,
                            0.37943052610218797, 0.23321847810377194};
    const double ref2[5] = {0.38044042007071682, 0.29121133580420611, 0.19571467800167705,
                            0.10290346097626511, 0.029730056887196644};
    double occ1 = 0.0, occ2 = 0.0;
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(P(k, 0) == doctest::Approx(ref0[k]).epsilon(FSOLVE_TOL));
        CHECK(P(k, 1) == doctest::Approx(ref1[k]).epsilon(FSOLVE_TOL));
        CHECK(P(k, 2) == doctest::Approx(ref2[k]).epsilon(FSOLVE_TOL));
        CHECK(P(k, 0) + P(k, 1) + P(k, 2) == doctest::Approx(1.0).epsilon(1e-14));
        occ1 += P(k, 1);
        occ2 += P(k, 2);
    }
    // MATLAB leaves 1.4e-7 here; the bisection fixed point closes it to 1e-11.
    CHECK(occ1 == doctest::Approx(m[0]).epsilon(1e-11));
    CHECK(occ2 == doctest::Approx(m[1]).epsilon(1e-11));

    SUBCASE("rates that fall once an item is cached") {
        Matrix<double> lam2(5, 3);
        for (std::size_t k = 0; k < 5; ++k) {
            lam2(k, 0) = base[k];
            lam2(k, 1) = base[k] / 2.0;
            lam2(k, 2) = base[k] / 4.0;
        }
        const Matrix<double> Q = cache_ttl_lrua(lam2, linear_graph(5, 2), m, 1e-13);
        const double q0[5] = {0.18052414308308898, 0.25139374011344356, 0.35374060649001621,
                              0.50171615756179688, 0.71262535277558292};
        const double q1[5] = {0.45407512591248711, 0.46173181112542644, 0.44573239744841425,
                              0.38636450904736891, 0.25209615644928551};
        const double q2[5] = {0.36540073100442388, 0.28687444876112989, 0.20052699606156946,
                              0.11191933339083415, 0.035278490775131593};
        double o1 = 0.0, o2 = 0.0;
        for (std::size_t k = 0; k < 5; ++k) {
            CHECK(Q(k, 0) == doctest::Approx(q0[k]).epsilon(FSOLVE_TOL));
            CHECK(Q(k, 1) == doctest::Approx(q1[k]).epsilon(FSOLVE_TOL));
            CHECK(Q(k, 2) == doctest::Approx(q2[k]).epsilon(FSOLVE_TOL));
            o1 += Q(k, 1);
            o2 += Q(k, 2);
        }
        CHECK(o1 == doctest::Approx(m[0]).epsilon(1e-11));
        CHECK(o2 == doctest::Approx(m[1]).epsilon(1e-11));
    }
}

TEST_CASE("cache_ttl_lrua handles never-requested items and rejects bad input") {
    // An item with zero request rate stays out of the cache with probability
    // one and takes no capacity from the others.
    const double base[4] = {1.0, 0.5, 0.0, 0.25};
    Matrix<double> lam(4, 2);
    for (std::size_t k = 0; k < 4; ++k) {
        lam(k, 0) = base[k];
        lam(k, 1) = base[k];
    }
    const Matrix<double> P =
        cache_ttl_lrua(lam, linear_graph(4, 1), std::vector<double>{1.0}, 1e-13);
    CHECK(P(2, 0) == 1.0);
    CHECK(P(2, 1) == 0.0);
    double occ = 0.0;
    for (std::size_t k = 0; k < 4; ++k) occ += P(k, 1);
    CHECK(occ == doctest::Approx(1.0).epsilon(1e-11));

    // An access graph that routes back into the not-cached node has no
    // characteristic time for that node; MATLAB indexes x(0) and errors out.
    std::vector<Matrix<double>> bad = linear_graph(4, 1);
    for (std::size_t k = 0; k < 4; ++k) bad[k](1, 0) = 1.0;
    CHECK_THROWS_AS(cache_ttl_lrua(lam, bad, std::vector<double>{1.0}, 1e-12), line::InputError);

    CHECK_THROWS_AS(cache_ttl_lrua(lam, linear_graph(4, 1), std::vector<double>{}, 1e-12),
                    line::InputError);
    CHECK_THROWS_AS(cache_ttl_lrua(lam, linear_graph(4, 2), std::vector<double>{1.0, 1.0}, 1e-12),
                    line::InputError);
}

TEST_CASE("the TTL cache approximations instantiate at high precision") {
    using line::Real50;
    // One list, two Poisson items of rate 1 and capacity 1: by symmetry the
    // characteristic time solves 2(1 - exp(-T)) = 1, i.e. T = log 2, and the
    // bisection at 50 digits must reproduce it far below double epsilon.
    std::vector<Map<Real50>> items;
    for (int k = 0; k < 2; ++k) {
        Map<Real50> p;
        p.D0 = Matrix<Real50>(1, 1, Real50(-1));
        p.D1 = Matrix<Real50>(1, 1, Real50(1));
        items.push_back(p);
    }
    const std::vector<Real50> m{Real50(1)};
    const std::vector<Real50> t = cache_t_lrum_map(items, m, Real50("1e-40"));
    const Real50 log2 = log(Real50(2));
    CHECK(static_cast<double>(abs(t[0] - log2)) < 1e-30);

    const CacheTtlLrumMapResult<Real50> r = cache_ttl_lrum_map(items, m, Real50("1e-40"));
    CHECK(static_cast<double>(abs(r.pijtime(0, 1) - Real50("0.5"))) < 1e-30);

    // cache_ttl_lrua on the same instance: one list, two items of rate 1,
    // capacity 1, so prob(k,1) = 1 - exp(-T) = 1/2.
    Matrix<Real50> lam(2, 2, Real50(1));
    std::vector<Matrix<Real50>> R;
    for (int k = 0; k < 2; ++k) {
        Matrix<Real50> Rk(2, 2, Real50(0));
        Rk(0, 1) = Real50(1);
        Rk(1, 1) = Real50(1);
        R.push_back(Rk);
    }
    const Matrix<Real50> P = cache_ttl_lrua(lam, R, m, Real50("1e-40"));
    CHECK(static_cast<double>(abs(P(0, 1) - Real50("0.5"))) < 1e-30);
    CHECK(static_cast<double>(abs(P(1, 0) - Real50("0.5"))) < 1e-30);
}
