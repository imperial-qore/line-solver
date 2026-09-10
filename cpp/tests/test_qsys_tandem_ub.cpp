/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Ciucu-Mehri tail bounds for a two-station tandem. The load-bearing check is
 * the exact reduction: for M/M/1 -> ./M/1 the five inequalities of Lemma 4 hold
 * as equalities, so the bound must return the exact tails (1+theta x)e^{-theta x}
 * and the Kraemer form for W. The remaining cases pin the digits produced by
 * matlab/src/api/qsys/qsys_tandem_ub_ciucu.m, which were themselves checked
 * against an exact CTMC reference for the Erlang(2)/M/1 tandem.
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_tandem_ub_ciucu.h"

using line::qsys::qsys_tandem_ub_ciucu;
using line::qsys::TandemUbResult;

namespace {
constexpr double REL_TOL = 1e-9;

std::function<double(const double&)> expLst(double rate) {
    return [rate](const double& s) { return rate / (rate + s); };
}
std::function<double(const double&)> expDlst(double rate) {
    return [rate](const double& s) { return rate / ((rate + s) * (rate + s)); };
}
std::function<double(const double&)> detLst(double d) {
    return [d](const double& s) { return std::exp(-s * d); };
}
std::function<double(const double&)> detDlst(double d) {
    return [d](const double& s) { return d * std::exp(-s * d); };
}
std::function<double(const double&)> erlang2Lst(double rate) {
    return [rate](const double& s) { return (rate / (rate + s)) * (rate / (rate + s)); };
}
std::function<double(const double&)> erlang2Dlst(double rate) {
    return [rate](const double& s) {
        return 2.0 * rate * rate / ((rate + s) * (rate + s) * (rate + s));
    };
}
}  // namespace

TEST_CASE("qsys_tandem_ub_ciucu is exact for M/M/1 -> ./M/1") {
    const double mu = 1.0, lambda = 0.5, theta = 0.5;
    const std::vector<double> x = {1.0, 5.0, 10.0};
    const std::vector<double> p = {1.0}, mus = {mu};
    const TandemUbResult<double> r =
        qsys_tandem_ub_ciucu<double>(x, expLst(lambda), p, mus, expDlst(lambda));

    CHECK(r.theta == doctest::Approx(theta).epsilon(1e-12));
    CHECK(r.A == doctest::Approx(1.0));
    CHECK(r.B == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.D == 0.0);
    for (std::size_t k = 0; k < x.size(); ++k) {
        const double S = (1.0 + theta * x[k]) * std::exp(-theta * x[k]);
        const double W = (1.0 - 2.0 * theta * theta / (mu * (mu + theta)) +
                          x[k] * (mu - theta) * theta / (mu + theta)) *
                         std::exp(-theta * x[k]);
        CHECK(r.S[k] == doctest::Approx(S).epsilon(REL_TOL));
        CHECK(r.W[k] == doctest::Approx(W).epsilon(REL_TOL));
    }
}

TEST_CASE("qsys_tandem_ub_ciucu matches MATLAB on D/M/1 -> ./M/1") {
    const std::vector<double> x = {5.0, 10.0, 20.0};
    const std::vector<double> p = {1.0}, mus = {1.0};
    const double d = 4.0 / 3.0;  // utilization 3/4
    const TandemUbResult<double> r = qsys_tandem_ub_ciucu<double>(x, detLst(d), p, mus, detDlst(d));

    CHECK(r.theta == doctest::Approx(0.454394983439251).epsilon(1e-10));
    CHECK(r.D > 0.0);  // arrivals less variable than Poisson put the bound on the D>0 branch
    CHECK(r.S[0] == doctest::Approx(0.493699069414895).epsilon(REL_TOL));
    CHECK(r.S[1] == doctest::Approx(0.09117765981667).epsilon(REL_TOL));
    CHECK(r.S[2] == doctest::Approx(0.00182565464670795).epsilon(REL_TOL));
    CHECK(r.W[0] == doctest::Approx(0.249922122549679).epsilon(REL_TOL));
}

TEST_CASE("qsys_tandem_ub_ciucu matches MATLAB on Erlang(2)/M/1 -> ./M/1") {
    const std::vector<double> x = {10.0, 20.0};
    const std::vector<double> p = {1.0}, mus = {1.0};
    const TandemUbResult<double> r =
        qsys_tandem_ub_ciucu<double>(x, erlang2Lst(1.0), p, mus, erlang2Dlst(1.0));

    CHECK(r.theta == doctest::Approx(0.618033988749894).epsilon(1e-10));
    CHECK(r.S[0] == doctest::Approx(0.0170463897138194).epsilon(REL_TOL));
    CHECK(r.S[1] == doctest::Approx(6.62788942098923e-05).epsilon(REL_TOL));
}

TEST_CASE("qsys_tandem_ub_ciucu matches MATLAB on Erlang(2)/H2/1 -> ./H2/1") {
    const double p2 = 0.9, p1 = 0.1, m2 = 1.69;
    const double m1 = p1 * m2 / (m2 - p2);          // CV(Y) = 2 with E[Y] = 1
    const double meanY = p1 / m1 + p2 / m2;
    const double rate = 2.0 / (meanY / 0.5);        // Erlang(2) arrivals at rho = 1/2
    const std::vector<double> x = {10.0, 25.0, 50.0};
    const std::vector<double> p = {p1, p2}, mus = {m1, m2};
    const TandemUbResult<double> r =
        qsys_tandem_ub_ciucu<double>(x, erlang2Lst(rate), p, mus, erlang2Dlst(rate));

    CHECK(r.theta == doctest::Approx(0.1500514416369).epsilon(1e-10));
    CHECK(r.S[0] == doctest::Approx(0.524484692983612).epsilon(REL_TOL));
    CHECK(r.S[1] == doctest::Approx(0.0877479694886629).epsilon(REL_TOL));
    CHECK(r.S[2] == doctest::Approx(0.0033336271735292).epsilon(REL_TOL));
    for (std::size_t k = 0; k < x.size(); ++k) CHECK(std::isnan(r.W[k]));  // Exp service only
}

TEST_CASE("qsys_tandem_ub_ciucu differentiates the transform when no derivative is given") {
    const std::vector<double> x = {5.0, 10.0};
    const std::vector<double> p = {1.0}, mus = {1.0};
    const double d = 4.0 / 3.0;
    const TandemUbResult<double> ref =
        qsys_tandem_ub_ciucu<double>(x, detLst(d), p, mus, detDlst(d));
    const TandemUbResult<double> num = qsys_tandem_ub_ciucu<double>(x, detLst(d), p, mus);
    CHECK(num.alpha == doctest::Approx(ref.alpha).epsilon(1e-10));
    for (std::size_t k = 0; k < x.size(); ++k)
        CHECK(num.S[k] == doctest::Approx(ref.S[k]).epsilon(1e-10));
}

TEST_CASE("qsys_tandem_ub_ciucu rejects an unstable tandem and malformed input") {
    const std::vector<double> x = {1.0};
    const std::vector<double> p = {1.0}, mus = {1.0};
    CHECK_THROWS_AS(qsys_tandem_ub_ciucu<double>(x, expLst(2.0), p, mus, expDlst(2.0)),
                    line::InputError);          // arrival rate above the service rate
    const std::vector<double> xneg = {-1.0};
    CHECK_THROWS_AS(qsys_tandem_ub_ciucu<double>(xneg, expLst(0.5), p, mus, expDlst(0.5)),
                    line::InputError);
    const std::vector<double> pbad = {0.3, 0.3}, mubad = {1.0, 2.0};
    CHECK_THROWS_AS(qsys_tandem_ub_ciucu<double>(x, expLst(0.5), pbad, mubad),
                    line::InputError);
}
