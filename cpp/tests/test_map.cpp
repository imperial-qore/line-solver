/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP descriptors. The oracles are the closed forms of processes whose
 * descriptors are known independently: a Poisson process, a two-phase Erlang,
 * and a two-state MMPP. Every identity that must hold exactly (SCV of Erlang-2
 * is 1/2, a renewal MAP has zero autocorrelation, IDC of Poisson is 1) is
 * asserted as an exact rational equality.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"

using line::Matrix;
using line::Rational;
using line::mam::Map;
using namespace line::mam;

namespace {

/** Erlang-2 with rate mu per phase: mean 2/mu, SCV 1/2. */
template <class T>
Map<T> erlang2(long mu_num, long mu_den) {
    const T mu = line::num_traits<T>::from_rational(mu_num, mu_den);
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D1 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D0(0, 0) = -mu;
    m.D0(0, 1) = mu;
    m.D0(1, 1) = -mu;
    m.D1(1, 0) = mu;
    return m;
}

/** Two-state MMPP: phase process with rates q01, q10, arrival rates l0, l1. */
template <class T>
Map<T> mmpp2(const T& l0, const T& l1, const T& q01, const T& q10) {
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D1 = Matrix<T>(2, 2, line::num_traits<T>::from_int(0));
    m.D0(0, 0) = -(l0 + q01);
    m.D0(0, 1) = q01;
    m.D0(1, 0) = q10;
    m.D0(1, 1) = -(l1 + q10);
    m.D1(0, 0) = l0;
    m.D1(1, 1) = l1;
    return m;
}

}  // namespace

TEST_CASE("Poisson process descriptors are exact") {
    const Rational lam(3, 2);
    Map<Rational> m = map_exponential(lam);
    CHECK(map_lambda(m) == lam);
    CHECK(map_mean(m) == Rational(2, 3));
    CHECK(map_moment(m, 2) == Rational(2) / (lam * lam));  // 2/lambda^2
    CHECK(map_scv(m) == Rational(1));                      // exponential
    CHECK(map_idc(m) == Rational(1));                      // Poisson counts
    CHECK(map_acf(m, std::vector<unsigned>{1, 2})[0] == Rational(0));  // renewal, no correlation
}

TEST_CASE("Erlang-2 descriptors match the closed form, exactly") {
    Map<Rational> m = erlang2<Rational>(2, 1);  // mu = 2 per phase
    CHECK(map_mean(m) == Rational(1));          // 2/mu
    CHECK(map_lambda(m) == Rational(1));
    CHECK(map_moment(m, 2) == Rational(3, 2));  // (k)(k+1)/mu^2 = 6/4
    CHECK(map_var(m) == Rational(1, 2));
    CHECK(map_scv(m) == Rational(1, 2));  // 1/k for Erlang-k
}

TEST_CASE("MMPP2 is more variable than Poisson and positively correlated") {
    Map<Rational> m = mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    CHECK(map_scv(m) > Rational(1));  // burstiness
    CHECK(map_idc(m) > Rational(1));  // over-dispersed counts

    // Autocorrelation decays: the lag-2 coefficient is closer to the renewal
    // value than lag-1 is (both computed exactly). map_acf is normalized as in
    // map_acf.m, so the renewal value is 0 and the coefficients lie in [-1, 1].
    const std::vector<Rational> acf = map_acf(m, std::vector<unsigned>{1, 2, 3});
    CHECK(acf[0] > Rational(0));
    CHECK(acf[0] < Rational(1));
    CHECK(line::num_abs(Rational(acf[1])) < line::num_abs(Rational(acf[0])));
    CHECK(line::num_abs(Rational(acf[2])) < line::num_abs(Rational(acf[1])));
}

TEST_CASE("map_prob and map_pie are proper distributions") {
    Map<Rational> m = mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    const std::vector<Rational> p = map_prob(m);
    const std::vector<Rational> pie = map_pie(m);
    Rational sp(0), spie(0);
    for (const Rational& v : p) sp += v;
    for (const Rational& v : pie) spie += v;
    CHECK(sp == Rational(1));
    CHECK(spie == Rational(1));
    // pi (D0 + D1) = 0 identically.
    const Matrix<Rational> Q = map_infgen(m);
    const std::vector<Rational> r = line::vecmul(p, Q);
    for (const Rational& v : r) CHECK(v == Rational(0));
}

TEST_CASE("the embedded chain is stochastic, exactly") {
    Map<Rational> m = mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    const Matrix<Rational> P = map_embedded(m);
    for (std::size_t i = 0; i < P.rows(); ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < P.cols(); ++j) s += P(i, j);
        CHECK(s == Rational(1));
    }
}

TEST_CASE("MAP descriptors agree between double and exact arithmetic") {
    Map<Rational> mq = mmpp2<Rational>(Rational(4), Rational(1), Rational(1, 2), Rational(1, 3));
    Map<double> md = mmpp2<double>(4.0, 1.0, 0.5, 1.0 / 3.0);
    CHECK(static_cast<double>(map_lambda(mq)) == doctest::Approx(map_lambda(md)).epsilon(1e-12));
    CHECK(static_cast<double>(map_scv(mq)) == doctest::Approx(map_scv(md)).epsilon(1e-12));
    CHECK(static_cast<double>(map_idc(mq)) == doctest::Approx(map_idc(md)).epsilon(1e-12));
}
