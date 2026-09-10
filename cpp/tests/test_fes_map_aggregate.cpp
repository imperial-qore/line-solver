/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP flow-equivalent server (Casale, Mi, Cherkasova and Smirni, IEEE Trans.
 * Soft. Eng. 37(5), 2011, Section 5.2). Oracles, in order of strength:
 *   1. The mean inter-departure time of the aggregated subnetwork is the
 *      reciprocal of its throughput, so for exponential service the recursion
 *      must reproduce EXACT MVA at every population: the MAP flow-equivalent
 *      server degrades to the classic Norton one, which is exact for
 *      product-form. Checked against pfqn_mva.
 *   2. A closed birth-death reference for a delay plus a two-server queue.
 *   3. Round trip through map2_fit_idc: the fitted MAP(2) must return the four
 *      descriptors it was fitted from.
 *   4. First-order convergence of the Euler quadrature to the linear solve.
 * The moment machinery uses field operations only, so case 1 also runs at
 * Rational, where MVA and the inter-departure mean must agree exactly.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/fes/fes_map_aggregate.h"
#include "line/api/fes/fes_map_interdeparture.h"
#include "line/api/fes/fes_map_interp.h"
#include "line/api/fes/fes_map_levels.h"
#include "line/api/fes/fes_map_moments.h"
#include "line/api/fes/fes_map_solve.h"
#include "line/api/mam/map2_fit_idc.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/pfqn/pfqn_mva.h"

using line::Matrix;
using line::Rational;
using line::fes::fes_map_aggregate;
using line::fes::fes_map_grid;
using line::fes::fes_map_interdeparture;
using line::fes::fes_map_interp;
using line::fes::fes_map_levels;
using line::fes::fes_map_deaggregate;
using line::fes::fes_map_moments;
using line::fes::fes_map_solve;
using line::mam::Map;
using line::mam::map2_fit_idc;
using line::mam::map_exponential_mean;
using line::mam::map_idc;
using line::mam::map_moment;
using line::mam::map_scv;

namespace {

constexpr double TOL = 1e-9;

/** Throughput of a single-class closed network with unit-server stations. */
template <class T>
T mva_throughput(const std::vector<T>& demands, std::size_t n, const std::vector<int>& servers) {
    Matrix<T> L(demands.size(), 1);
    for (std::size_t i = 0; i < demands.size(); ++i) L(i, 0) = demands[i];
    const std::vector<int> N(1, static_cast<int>(n));
    Matrix<T> Z(1, 1, line::num_traits<T>::from_int(0));
    return line::pfqn::pfqn_mva(L, N, Z, servers).XN[0];
}

/** Stationary throughput of a delay of rate 2 feeding a two-server queue of rate 1. */
double birth_death_throughput(std::size_t k) {
    std::vector<double> p(k + 1, 1.0);
    for (std::size_t j = 1; j <= k; ++j)
        p[j] = p[j - 1] * (static_cast<double>(k - j + 1) * 2.0) / std::min<double>(j, 2.0);
    double s = 0;
    for (std::size_t j = 0; j <= k; ++j) s += p[j];
    double x = 0;
    for (std::size_t j = 0; j <= k; ++j) x += (p[j] / s) * std::min<double>(j, 2.0);
    return x;
}

}  // namespace

TEST_CASE("fes_map interdeparture mean is the inverse throughput") {
    const std::vector<Map<double>> station = fes_map_levels(map_exponential_mean(1.0 / 1.5), 10);
    const std::vector<Map<double>> fes = fes_map_levels(map_exponential_mean(1.0), 10);
    for (std::size_t n : {1u, 2u, 5u, 10u}) {
        const Map<double> T = fes_map_interdeparture(station, fes, n);
        const line::fes::FesMapMoments<double> mom = fes_map_moments(T);
        const double x = mva_throughput<double>({1.0 / 1.5, 1.0}, n, {1, 1});
        CHECK(std::abs(mom.e1 - 1.0 / x) <= TOL / x);
    }
}

TEST_CASE("fes_map interdeparture mean is exact at Rational") {
    const std::vector<Map<Rational>> station =
        fes_map_levels(map_exponential_mean(Rational(2, 3)), 4);
    const std::vector<Map<Rational>> fes = fes_map_levels(map_exponential_mean(Rational(1)), 4);
    for (std::size_t n = 1; n <= 4; ++n) {
        const Map<Rational> T = fes_map_interdeparture(station, fes, n);
        const line::fes::FesMapMoments<Rational> mom = fes_map_moments(T);
        const Rational x = mva_throughput<Rational>({Rational(2, 3), Rational(1)}, n, {1, 1});
        CHECK(mom.e1 == Rational(1) / x);
    }
}

TEST_CASE("fes_map euler quadrature converges with the step") {
    const std::vector<Map<double>> station = fes_map_levels(map_exponential_mean(0.8), 6);
    const std::vector<Map<double>> fes = fes_map_levels(map_exponential_mean(1.3), 6);
    const Map<double> T = fes_map_interdeparture(station, fes, 6);
    const double exact = fes_map_moments(T, "ssolve").e1;
    const double coarse = std::abs(fes_map_moments(T, "euler", 0.1).e1 - exact);
    const double fine = std::abs(fes_map_moments(T, "euler", 0.01).e1 - exact);
    CHECK(coarse / exact < 5e-2);
    CHECK(fine < 0.2 * coarse);
}

TEST_CASE("fes_map aggregation of exponential stations reproduces exact MVA") {
    const std::vector<double> rates = {2.0, 1.5, 1.0};
    std::vector<Map<double>> maps;
    std::vector<double> demands;
    for (double r : rates) {
        maps.push_back(map_exponential_mean(1.0 / r));
        demands.push_back(1.0 / r);
    }
    const std::size_t nmax = 8;
    const line::fes::FesMapAggregateResult<double> res =
        fes_map_aggregate(maps, {1, 1, 1}, nmax);
    for (std::size_t k = 1; k <= nmax; ++k) {
        const double x = mva_throughput<double>(demands, k, {1, 1, 1});
        CHECK(std::abs(res.throughput[k - 1] - x) <= TOL * x);
    }
}

TEST_CASE("fes_map aggregation of a delay and a multiserver station") {
    const std::vector<Map<double>> maps = {map_exponential_mean(0.5), map_exponential_mean(1.0)};
    const std::size_t nmax = 6;
    const line::fes::FesMapAggregateResult<double> res = fes_map_aggregate(
        maps, {std::numeric_limits<double>::infinity(), 2.0}, nmax);
    for (std::size_t k = 1; k <= nmax; ++k) {
        const double x = birth_death_throughput(k);
        CHECK(std::abs(res.throughput[k - 1] - x) <= 1e-8 * x);
    }
}

TEST_CASE("fes_map aggregation keeps the burstiness of the subnetwork") {
    Map<double> bursty;
    bursty.D0 = Matrix<double>(2, 2, 0.0);
    bursty.D1 = Matrix<double>(2, 2, 0.0);
    bursty.D0(0, 0) = -1.9;
    bursty.D0(1, 1) = -0.1;
    bursty.D1(0, 0) = 1.71;
    bursty.D1(0, 1) = 0.19;
    bursty.D1(1, 0) = 0.01;
    bursty.D1(1, 1) = 0.09;
    const std::vector<Map<double>> maps = {bursty, map_exponential_mean(1.2)};
    const line::fes::FesMapAggregateResult<double> res = fes_map_aggregate(maps, {1, 1}, 5);
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(res.status[k] == 0);
        CHECK(std::abs(map_idc(res.fes[k]) - res.moments[3][k]) <= 1e-6 * res.moments[3][k]);
        CHECK(std::abs(map_moment(res.fes[k], 1) - res.moments[0][k]) <= 1e-9);
        CHECK(map_scv(res.fes[k]) > 1.0);
    }
}

TEST_CASE("fes_map fit falls back to an exponential when the aggregate is not bursty") {
    const line::mam::Map2FitIdcResult<double> fit = map2_fit_idc(1.0, 1.5, 4.0, 0.7);
    CHECK(fit.status == 1);
    CHECK(std::abs(map_moment(fit.map, 1) - 1.0) <= TOL);
}

TEST_CASE("fes_map grid skips levels only above twenty") {
    CHECK(fes_map_grid(20).size() == 20);
    const std::vector<std::size_t> grid = fes_map_grid(100);
    CHECK(grid.size() < 100);
    CHECK(grid.front() == 1);
    CHECK(grid.back() == 100);
}

TEST_CASE("fes_map interpolation is shape preserving") {
    const std::vector<double> x = {1, 2, 5, 9};
    const std::vector<double> y = {1.0, 2.0, 2.5, 2.6};
    std::vector<double> xq;
    for (int i = 0; i <= 40; ++i) xq.push_back(1.0 + 8.0 * i / 40.0);
    const std::vector<double> yq = fes_map_interp(x, y, xq);
    for (std::size_t i = 1; i < yq.size(); ++i) CHECK(yq[i] >= yq[i - 1] - 1e-12);
    for (std::size_t i = 0; i < yq.size(); ++i) CHECK(yq[i] <= 2.6 + 1e-12);
}

TEST_CASE("fes_map reduced model solve matches exact MVA") {
    // the closed model of a delay and an aggregated exponential subnetwork is
    // product-form, so the reduced solve must return exact MVA
    const std::vector<double> rates = {2.0, 1.5, 1.0};
    const double Z = 0.8;
    std::vector<Map<double>> maps;
    std::vector<double> demands;
    for (double r : rates) {
        maps.push_back(map_exponential_mean(1.0 / r));
        demands.push_back(1.0 / r);
    }
    Matrix<double> L(3, 1);
    for (std::size_t i = 0; i < 3; ++i) L(i, 0) = demands[i];
    Matrix<double> Zm(1, 1, Z);

    for (std::size_t n : {1u, 3u, 6u, 10u}) {
        const line::fes::FesMapAggregateResult<double> agg = fes_map_aggregate(maps, {1, 1, 1}, n);
        const line::fes::FesMapSolveResult<double> sol =
            fes_map_solve(agg.fes, map_exponential_mean(Z), n);
        const std::vector<int> N(1, static_cast<int>(n));
        const line::pfqn::MvaResult<double> ref =
            line::pfqn::pfqn_mva(L, N, Zm, std::vector<int>{1, 1, 1});
        CHECK(std::abs(sol.X - ref.XN[0]) <= TOL * ref.XN[0]);
        double R = 0;
        for (std::size_t i = 0; i < 3; ++i) R += ref.CN(i, 0);
        CHECK(std::abs(sol.R - R) <= 1e-8 * R);
        double mass = 0;
        for (double p : sol.pk) mass += p;
        CHECK(std::abs(mass - 1.0) <= 1e-12);
    }
}

TEST_CASE("fes_map deaggregation matches exact MVA per station") {
    const std::vector<double> rates = {2.0, 1.5, 1.0};
    const double Z = 0.8;
    const std::size_t n = 8;
    std::vector<Map<double>> maps;
    std::vector<double> demands;
    for (double r : rates) {
        maps.push_back(map_exponential_mean(1.0 / r));
        demands.push_back(1.0 / r);
    }
    const line::fes::FesMapAggregateResult<double> agg = fes_map_aggregate(maps, {1, 1, 1}, n);
    const line::fes::FesMapSolveResult<double> sol =
        fes_map_solve(agg.fes, map_exponential_mean(Z), n);
    const line::fes::FesMapDeaggregateResult<double> de =
        fes_map_deaggregate(sol.pk, demands, std::vector<int>{1, 1, 1},
                            std::vector<bool>{false, false, false});

    Matrix<double> L(3, 1);
    for (std::size_t i = 0; i < 3; ++i) L(i, 0) = demands[i];
    Matrix<double> Zm(1, 1, Z);
    const std::vector<int> N(1, static_cast<int>(n));
    const line::pfqn::MvaResult<double> ref =
        line::pfqn::pfqn_mva(L, N, Zm, std::vector<int>{1, 1, 1});
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(std::abs(de.Q[i] - ref.QN(i, 0)) <= 1e-8 * ref.QN(i, 0));
        CHECK(std::abs(de.U[i] - ref.UN(i, 0)) <= 1e-8 * ref.UN(i, 0));
        CHECK(std::abs(de.X[i] - ref.XN[0]) <= 1e-8 * ref.XN[0]);
    }
}
