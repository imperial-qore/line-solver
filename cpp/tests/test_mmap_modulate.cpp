/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mmap_modulate and mmap_mixture_order2.
 *
 * ORACLES. A modulated process has closed-form aggregate rates that the
 * construction must reproduce and that share no code with it: the stationary
 * environment law is the invariant vector of P, and the per-class arrival rate
 * of the modulated process is the P-weighted average of the components' rates
 * weighted by their mean holding times,
 *
 *   lambda_k = sum_j w_j E[HT_j] lambda_k(j) / sum_j w_j E[HT_j],
 *
 * with w the invariant law of P. The structural invariants (D0 + sum_c Dc a
 * generator, Dc non-negative, sum_c Dc = D1) hold for any marked MAP.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_modulate.h"
#include "line/api/mam/mmap_stats.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

/** A one-class MMAP from a MAP. */
mam::Mmap<double> one_class(const mam::Map<double>& m) {
    mam::Mmap<double> x;
    x.D0 = m.D0;
    x.D1 = m.D1;
    x.Dc.assign(1, m.D1);
    return x;
}

void check_is_mmap(const mam::Mmap<double>& m) {
    const std::size_t n = m.order();
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-12);
            double dc = 0.0;
            for (std::size_t c = 0; c < m.classes(); ++c) {
                CHECK(m.Dc[c](i, j) >= -1e-12);
                dc += m.Dc[c](i, j);
            }
            CHECK(m.D1(i, j) == doctest::Approx(dc).epsilon(1e-12));
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-9);
    }
}

}  // namespace

TEST_CASE("mmap_modulate builds a marked MAP of the right order and shape") {
    // Two environments, exponential holding times of means 1 and 2, Poisson
    // arrivals of rates 3 and 0.5.
    std::vector<mam::Map<double>> HT;
    HT.push_back(mam::map_exponential_mean(1.0));
    HT.push_back(mam::map_exponential_mean(2.0));
    std::vector<mam::Mmap<double>> C;
    C.push_back(one_class(mam::map_exponential_mean(1.0 / 3.0)));
    C.push_back(one_class(mam::map_exponential_mean(2.0)));

    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;

    const mam::Mmap<double> m = mam::mmap_modulate(P, HT, C);
    CHECK(m.order() == 2);  // 1*1 + 1*1
    CHECK(m.classes() == 1);
    check_is_mmap(m);

    // Alternating environments: the invariant law is (1/2, 1/2), so the time
    // shares are E[HT]/sum = 1/3 and 2/3 and the arrival rate is the weighted
    // average 3*(1/3) + 0.5*(2/3) = 4/3.
    CHECK(mam::map_lambda(m.map()) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
}

TEST_CASE("a modulated process reduces to its component when the environment is one") {
    // A single environment cannot switch, so the result is that component.
    std::vector<mam::Map<double>> HT(1, mam::map_exponential_mean(1.0));
    std::vector<mam::Mmap<double>> C(1, one_class(mam::map_erlang(0.5, 2)));
    Matrix<double> P(1, 1, 1.0);

    const mam::Mmap<double> m = mam::mmap_modulate(P, HT, C);
    check_is_mmap(m);
    CHECK(mam::map_lambda(m.map()) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(mam::map_scv(m.map()) == doctest::Approx(0.5).epsilon(1e-9));
}

TEST_CASE("the per-class rates follow the environment time shares") {
    // Two environments, two classes each, with different class splits.
    std::vector<mam::Map<double>> HT;
    HT.push_back(mam::map_exponential_mean(1.0));
    HT.push_back(mam::map_exponential_mean(3.0));

    std::vector<mam::Mmap<double>> C;
    for (int e = 0; e < 2; ++e) {
        const double lam = (e == 0) ? 4.0 : 1.0;
        const double split = (e == 0) ? 0.25 : 0.75;  // share of class 0
        mam::Mmap<double> x;
        x.D0 = Matrix<double>(1, 1, -lam);
        x.D1 = Matrix<double>(1, 1, lam);
        x.Dc.assign(2, Matrix<double>(1, 1, 0.0));
        x.Dc[0](0, 0) = lam * split;
        x.Dc[1](0, 0) = lam * (1.0 - split);
        C.push_back(x);
    }
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;

    const mam::Mmap<double> m = mam::mmap_modulate(P, HT, C);
    check_is_mmap(m);
    CHECK(m.classes() == 2);

    // Time shares 1/4 and 3/4; class 0 rate = 4*0.25*(1/4) + 1*0.75*(3/4).
    const std::vector<double> lk = mam::mmap_lambda(m);
    CHECK(lk[0] == doctest::Approx(4.0 * 0.25 * 0.25 + 1.0 * 0.75 * 0.75).epsilon(1e-9));
    CHECK(lk[1] == doctest::Approx(4.0 * 0.75 * 0.25 + 1.0 * 0.25 * 0.75).epsilon(1e-9));
    CHECK(lk[0] + lk[1] == doctest::Approx(mam::map_lambda(m.map())).epsilon(1e-10));
}

TEST_CASE("mmap_modulate refuses mismatched inputs by name") {
    std::vector<mam::Map<double>> HT(2, mam::map_exponential_mean(1.0));
    std::vector<mam::Mmap<double>> C(1, one_class(mam::map_exponential_mean(1.0)));
    Matrix<double> P(2, 2, 0.5);
    CHECK_THROWS_AS(mam::mmap_modulate(P, HT, C), line::InputError);

    // Disagreeing class counts.
    std::vector<mam::Mmap<double>> C2;
    C2.push_back(one_class(mam::map_exponential_mean(1.0)));
    mam::Mmap<double> two = one_class(mam::map_exponential_mean(1.0));
    two.Dc.push_back(Matrix<double>(1, 1, 0.0));
    C2.push_back(two);
    CHECK_THROWS_AS(mam::mmap_modulate(P, HT, C2), line::InputError);
}

TEST_CASE("mmap_mixture_order2 records the previous and the current class") {
    const std::size_t m = 2;
    std::vector<std::vector<mam::Map<double>>> PHs(m, std::vector<mam::Map<double>>(m));
    // Distinct two-phase sojourns so the pair structure is observable.
    PHs[0][0] = mam::map_erlang(1.0, 2);
    PHs[0][1] = mam::map_erlang(2.0, 2);
    PHs[1][0] = mam::map_erlang(0.5, 2);
    PHs[1][1] = mam::map_erlang(1.5, 2);

    Matrix<double> P2(2, 2, 0.0);
    P2(0, 0) = 0.6; P2(0, 1) = 0.4;
    P2(1, 0) = 0.3; P2(1, 1) = 0.7;

    const mam::Mmap<double> x = mam::mmap_mixture_order2(PHs, P2);
    CHECK(x.order() == 2 * m * m);
    CHECK(x.classes() == m);
    check_is_mmap(x);
    CHECK(mam::map_lambda(x.map()) > 0.0);
}

TEST_CASE("mmap_mixture_order2 refuses a non-square table or a wrong order") {
    std::vector<std::vector<mam::Map<double>>> bad(2, std::vector<mam::Map<double>>(1));
    bad[0][0] = mam::map_erlang(1.0, 2);
    bad[1][0] = mam::map_erlang(1.0, 2);
    Matrix<double> P2(2, 2, 0.5);
    CHECK_THROWS_AS(mam::mmap_mixture_order2(bad, P2), line::InputError);

    std::vector<std::vector<mam::Map<double>>> one(1, std::vector<mam::Map<double>>(1));
    one[0][0] = mam::map_erlang(1.0, 3);  // order three, not two
    Matrix<double> P1(1, 1, 1.0);
    CHECK_THROWS_AS(mam::mmap_mixture_order2(one, P1), line::InputError);
}

TEST_CASE("mmap_mixture_fit builds a marked MAP over the class pairs") {
    // Two classes, so four ordered pairs and a state space of 2*2^2 = 8.
    const std::size_t m = 2;
    line::Matrix<double> M1(m, m, 0.0), M2(m, m, 0.0), M3(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            const double mu = 1.0 + 0.5 * static_cast<double>(i) + 0.25 * static_cast<double>(j);
            const double scv = 2.0;
            M1(i, j) = mu;
            M2(i, j) = (1.0 + scv) * mu * mu;
            M3(i, j) = 3.0 * 1.5 * M2(i, j) * M2(i, j) / mu;
        }
    // A flattened triple sigma, (m x m*m), entry (i, j*m + h).
    line::Matrix<double> P2(m, m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            for (std::size_t h = 0; h < m; ++h)
                P2(i, j * m + h) = 0.1 + 0.05 * static_cast<double>(i + j + h);

    const mam::Mmap<double> x = mam::mmap_mixture_fit(P2, M1, M2, M3);
    CHECK(x.order() == 2 * m * m);
    CHECK(x.classes() == m);
    check_is_mmap(x);
    CHECK(mam::map_lambda(x.map()) > 0.0);
}

TEST_CASE("mmap_mixture_fit refuses tables of the wrong shape by name") {
    line::Matrix<double> M(2, 2, 1.0), bad(2, 3, 1.0), P2(2, 4, 0.25);
    CHECK_THROWS_AS(mam::mmap_mixture_fit(P2, bad, M, M), line::InputError);
    line::Matrix<double> badP2(2, 3, 0.25);
    CHECK_THROWS_AS(mam::mmap_mixture_fit(badP2, M, M, M), line::InputError);
}

TEST_CASE("mmap_mixture_fit_trace measures its inputs from the trace") {
    // Every ORDERED PAIR must occur, or the pair's component has no cross
    // moment to fit; a strictly alternating label sequence never produces the
    // (1,1) or (2,2) pairs and is refused for exactly that reason.
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 2000; ++i) {
        Tv.push_back(0.5 + 0.5 * static_cast<double>((i * 7919) % 11) / 11.0);
        A.push_back(static_cast<int>((i % 5) < 2) + 1);  // 2,2,1,1,1: all four pairs occur
    }
    const mam::Mmap<double> x = mam::mmap_mixture_fit_trace(Tv, A);
    check_is_mmap(x);
    CHECK(x.classes() == 2);
    CHECK(mam::map_lambda(x.map()) > 0.0);

    CHECK_THROWS_AS(mam::mmap_mixture_fit_trace(std::vector<double>(), std::vector<int>()),
                    line::InputError);
}

TEST_CASE("a class pair that never occurs is refused by name") {
    // Strictly alternating labels: the (1,1) and (2,2) pairs never happen.
    std::vector<double> Tv;
    std::vector<int> A;
    for (std::size_t i = 0; i < 500; ++i) {
        Tv.push_back(1.0);
        A.push_back(static_cast<int>(i % 2) + 1);
    }
    CHECK_THROWS_AS(mam::mmap_mixture_fit_trace(Tv, A), line::InputError);
}
