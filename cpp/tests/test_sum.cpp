/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Summation method (SUM/ESUM) and the closing method. Oracles:
 *   1. The population constraint the method is defined by: the returned
 *      throughputs must satisfy sum_i Q(i,r) + X_r Z_r = N_r, per class, to
 *      the bisection tolerance. That is the equation being solved, so it is
 *      the only law available exactly.
 *   2. The operational identities U = X L / m and R = Q / X, asserted as
 *      value identities.
 *   3. The throughput never exceeds the saturation bound min_i m_i / L(i,r),
 *      which the bisection brackets from above by construction.
 *   4. The Erlang-C waiting probability against its hand value, asserted
 *      EXACTLY in rational arithmetic (it is a finite field computation even
 *      though the bisection around it is not).
 *   5. MATLAB reference values from matlab/src/api/sum/*.m.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/sum/sum_closed.h"
#include "line/api/sum/sum_closing.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::sum::ClosingOptions;
using line::sum::Servers;
using line::sum::sum_closed;
using line::sum::sum_closing;
using line::sum::SumOptions;

namespace {

double rel(double got, double want) {
    const double d = std::fabs(got - want);
    return std::fabs(want) > 1e-30 ? d / std::fabs(want) : d;
}

SumOptions tight() {
    SumOptions o;
    o.tol = 1e-9;
    o.maxiter = 100000;
    return o;
}

}  // namespace

TEST_CASE("sum_erlangc hand value, exactly") {
    // M/M/2 at rho = 1/2: a = 1, the tail term is 1/(2!(1-rho)) = 1, the head
    // sum is 1 + 1 = 2, so C = 1/3. A finite field computation: the a^k/k!
    // terms are built incrementally, no factorial and no real power.
    CHECK(line::sum::detail::sum_erlangc<Rational>(2, Rational(1, 2)) == Rational(1, 3));
    // M/M/1 collapses to C = rho.
    CHECK(line::sum::detail::sum_erlangc<Rational>(1, Rational(3, 8)) == Rational(3, 8));
    // Saturation is answered with certainty of waiting.
    CHECK(line::sum::detail::sum_erlangc<Rational>(3, Rational(1)) == Rational(1));
    CHECK(line::sum::detail::sum_erlangc<Rational>(3, Rational(5, 4)) == Rational(1));
}

TEST_CASE("sum_closed single class matches MATLAB sum_closed.m") {
    // L = [1;2], N = 3, Z = 1, single servers, exponential service. The
    // bisection is run to tol = 1e-9 on both sides; the comparison is held at
    // 1e-8 relative, which is the accuracy the stopping rule guarantees.
    const Matrix<double> L{{1.0}, {2.0}};
    const std::vector<long> N{3};
    const std::vector<double> Z{1.0};
    const std::vector<Servers> mi{Servers::of(1), Servers::of(1)};
    const Matrix<double> scv{{1.0}, {1.0}};

    const auto r = sum_closed(L, N, Z, mi, scv, tight());
    CHECK(rel(r.XN[0], 0.426497321110219) < 1e-8);
    CHECK(rel(r.QN(0, 0), 0.595942603820006) < 1e-8);
    CHECK(rel(r.QN(1, 0), 1.97756007418811) < 1e-8);
    CHECK(rel(r.UN(0, 0), 0.426497321110219) < 1e-8);
    CHECK(rel(r.UN(1, 0), 0.852994642220438) < 1e-8);
    CHECK(rel(r.RN(0, 0), 1.39729506921334) < 1e-8);
    CHECK(rel(r.RN(1, 0), 4.63674676558415) < 1e-8);

    // the population constraint the method solves
    const double pop = r.QN(0, 0) + r.QN(1, 0) + r.XN[0] * Z[0];
    CHECK(std::fabs(pop - 3.0) < 1e-8);
}

TEST_CASE("sum_closed ESUM node functions match MATLAB sum_closed.m") {
    // A 2-server station with scv = 4 (so Eq. 10.89 with the Erlang-C factor)
    // and a single-server station with scv = 1/4 (Eq. 10.88).
    const Matrix<double> L{{1.0}, {2.0}};
    const std::vector<long> N{5};
    const std::vector<double> Z{0.5};
    const std::vector<Servers> mi{Servers::of(2), Servers::of(1)};
    const Matrix<double> scv{{4.0}, {0.25}};

    const auto r = sum_closed(L, N, Z, mi, scv, tight());
    CHECK(rel(r.XN[0], 0.4850802232977) < 1e-8);
    CHECK(rel(r.QN(0, 0), 0.544911825904915) < 1e-8);
    CHECK(rel(r.QN(1, 0), 4.21254806350151) < 1e-8);
    CHECK(rel(r.UN(0, 0), 0.24254011164885) < 1e-8);
    CHECK(rel(r.UN(1, 0), 0.970160446595401) < 1e-8);

    // U is per server at a queueing station
    CHECK(std::fabs(r.UN(0, 0) - r.XN[0] * L(0, 0) / 2.0) < 1e-15);
    CHECK(std::fabs(r.UN(1, 0) - r.XN[0] * L(1, 0) / 1.0) < 1e-15);
    // and the saturation bound brackets the answer
    CHECK(r.XN[0] <= 2.0 / L(0, 0) + 1e-12);
    CHECK(r.XN[0] <= 1.0 / L(1, 0) + 1e-12);
}

TEST_CASE("sum_closed multiclass matches MATLAB sum_closed.m") {
    // Two classes, one single-server and one two-server station, mixed SCVs.
    Matrix<double> L(2, 2);
    L(0, 0) = 1.0;
    L(0, 1) = 2.0;
    L(1, 0) = 3.0;
    L(1, 1) = 1.0;
    const std::vector<long> N{2, 3};
    const std::vector<double> Z{1.0, 0.0};
    const std::vector<Servers> mi{Servers::of(1), Servers::of(2)};
    Matrix<double> scv(2, 2);
    scv(0, 0) = 1.0;
    scv(0, 1) = 1.0;
    scv(1, 0) = 2.0;
    scv(1, 1) = 0.5;

    const auto r = sum_closed(L, N, Z, mi, scv, tight());
    CHECK(rel(r.XN[0], 0.230570135727421) < 1e-7);
    CHECK(rel(r.XN[1], 0.342046962352058) < 1e-7);
    CHECK(rel(r.QN(0, 0), 0.859474442346299) < 1e-7);
    CHECK(rel(r.QN(1, 0), 0.909955422668779) < 1e-7);
    CHECK(rel(r.QN(0, 1), 2.55003208716781) < 1e-7);
    CHECK(rel(r.QN(1, 1), 0.44996791282573) < 1e-7);
    CHECK(rel(r.UN(0, 0), 0.230570135727421) < 1e-7);
    CHECK(rel(r.UN(1, 0), 0.345855203591132) < 1e-7);
    CHECK(rel(r.RN(0, 0), 3.72760522361129) < 1e-7);
    CHECK(rel(r.RN(1, 1), 1.31551500920097) < 1e-7);

    // per-class population constraints, the equations the sweeps solve
    for (std::size_t c = 0; c < 2; ++c) {
        const double pop = r.QN(0, c) + r.QN(1, c) + r.XN[c] * Z[c];
        CHECK(std::fabs(pop - static_cast<double>(N[c])) < 1e-6);
    }
    // residence times are queue lengths over throughputs, identically
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(std::fabs(r.RN(i, c) * r.XN[c] - r.QN(i, c)) < 1e-14);
}

TEST_CASE("sum_closed with an infinite-server station matches MATLAB") {
    Matrix<double> L(2, 2);
    L(0, 0) = 1.0;
    L(0, 1) = 2.0;
    L(1, 0) = 3.0;
    L(1, 1) = 1.0;
    const std::vector<Servers> mi{Servers::of(1), Servers::inf()};
    const auto r = sum_closed(L, std::vector<long>{2, 3}, std::vector<double>{0.0, 0.0}, mi,
                              Matrix<double>(2, 2, 1.0), tight());
    CHECK(rel(r.XN[0], 0.283430614359664) < 1e-7);
    CHECK(rel(r.XN[1], 0.329207214347883) < 1e-7);
    CHECK(rel(r.QN(0, 0), 1.14970815807226) < 1e-7);
    CHECK(rel(r.QN(1, 0), 0.850291843078992) < 1e-7);
    CHECK(rel(r.QN(0, 1), 2.67079278564955) < 1e-7);
    CHECK(rel(r.QN(1, 1), 0.329207214347883) < 1e-7);

    // an infinite-server station holds exactly X L jobs, and its utilization
    // is reported as that same offered load rather than per server
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(std::fabs(r.QN(1, c) - r.XN[c] * L(1, c)) < 1e-14);
        CHECK(std::fabs(r.UN(1, c) - r.XN[c] * L(1, c)) < 1e-14);
    }
}

TEST_CASE("sum_closed edge cases") {
    const Matrix<double> L{{1.0}, {2.0}};
    // an empty population leaves everything at zero and never iterates
    const auto z = sum_closed(L, std::vector<long>{0}, std::vector<double>{0.0},
                              std::vector<Servers>{Servers::of(1), Servers::of(1)},
                              Matrix<double>(2, 1, 1.0));
    CHECK(z.XN[0] == 0.0);
    CHECK(z.it == 0);
    // all demands zero and no think time is not a model
    CHECK_THROWS_AS(sum_closed(Matrix<double>(2, 1, 0.0), std::vector<long>{2},
                               std::vector<double>{0.0},
                               std::vector<Servers>{Servers::of(1), Servers::of(1)},
                               Matrix<double>(2, 1, 1.0)),
                    line::InputError);
    // a negative population is rejected
    CHECK_THROWS_AS(sum_closed(L, std::vector<long>{-1}, std::vector<double>{0.0},
                               std::vector<Servers>{Servers::of(1), Servers::of(1)},
                               Matrix<double>(2, 1, 1.0)),
                    line::InputError);
}

TEST_CASE("sum_closing matches MATLAB sum_closing.m") {
    // One open class at rate 1/2 and one closed class of population 2, over
    // two single-server stations. Kclosed = 200.
    const std::vector<double> lambda0{0.5, 0.0};
    const std::vector<double> scva{1.0, 1.0};
    Matrix<double> L(2, 2);
    L(0, 0) = 1.0;
    L(0, 1) = 2.0;
    L(1, 0) = 0.5;
    L(1, 1) = 1.0;
    const std::vector<Servers> mi{Servers::of(1), Servers::of(1)};
    const Matrix<double> scv(2, 2, 1.0);
    const std::vector<long> N{0, 2};  // the open class entry is ignored
    const std::vector<double> Z{0.0, 1.0};

    ClosingOptions opt;
    opt.Kclosed = 200;
    opt.sum = tight();
    const auto r = sum_closing(lambda0, scva, L, mi, scv, N, Z, opt);

    CHECK(rel(r.XN[0], 0.499937392478842) < 1e-7);
    CHECK(rel(r.XN[1], 0.154597075828142) < 1e-7);
    CHECK(rel(r.QN(0, 0), 2.56543841041066) < 1e-6);
    CHECK(rel(r.QN(1, 0), 0.418401746453664) < 1e-6);
    CHECK(rel(r.QN(0, 1), 1.58663577653264) < 1e-6);
    CHECK(rel(r.QN(1, 1), 0.258767147631838) < 1e-6);
    CHECK(rel(r.UN(0, 0), 0.499937392478842) < 1e-7);
    CHECK(rel(r.UN(1, 1), 0.154597075828142) < 1e-7);
    CHECK(rel(r.RN(0, 0), 5.13151936423566) < 1e-6);
    CHECK(rel(r.TN[0], 5.96842765064949) < 1e-6);
    CHECK(rel(r.TN[1], 11.936855301299) < 1e-6);

    // the open-class throughput approaches lambda0 from below
    CHECK(r.XN[0] < 0.5);
    CHECK(r.XN[0] > 0.499);
    // and the closed class still satisfies its population constraint
    const double pop = r.QN(0, 1) + r.QN(1, 1) + r.XN[1] * Z[1];
    CHECK(std::fabs(pop - 2.0) < 1e-5);
    // Little's law on the original stations
    for (std::size_t c = 0; c < 2; ++c)
        CHECK(std::fabs(r.TN[c] * r.XN[c] - (r.QN(0, c) + r.QN(1, c))) < 1e-12);
}

TEST_CASE("sum_closing raises the open throughput toward lambda0 with Kclosed") {
    // The closing station's demand is 1/(Ropen lambda0), so a larger closing
    // population pushes the open class closer to its offered rate. Monotone by
    // construction of the method, and the test asserts that direction.
    const std::vector<double> lambda0{0.4};
    const std::vector<double> scva{1.0};
    const Matrix<double> L{{1.0}};
    const std::vector<Servers> mi{Servers::of(1)};
    const Matrix<double> scv{{1.0}};

    double prev = 0.0;
    for (long K : {50L, 200L, 1000L}) {
        ClosingOptions opt;
        opt.Kclosed = K;
        opt.sum = tight();
        const auto r =
            sum_closing(lambda0, scva, L, mi, scv, std::vector<long>{0}, std::vector<double>{0.0},
                        opt);
        CHECK(r.XN[0] < 0.4);
        CHECK(r.XN[0] > prev);
        prev = r.XN[0];
    }
    CHECK(prev > 0.399);
}

TEST_CASE("sum_closing rejects a network with no open class") {
    const std::vector<double> lambda0{0.0};
    CHECK_THROWS_AS(sum_closing(lambda0, std::vector<double>{1.0}, Matrix<double>{{1.0}},
                                std::vector<Servers>{Servers::of(1)}, Matrix<double>{{1.0}},
                                std::vector<long>{2}, std::vector<double>{0.0}),
                    line::InputError);
}

TEST_CASE("sum_closed at high precision agrees with the double path") {
    // The method is arithmetic-agnostic: at Real50 the same bisection lands on
    // the same root, so the two paths agree far beyond the stopping tolerance
    // once both are driven to it.
    const Matrix<double> Ld{{1.0}, {2.0}};
    Matrix<Real50> Lr(2, 1);
    Lr(0, 0) = Real50(1);
    Lr(1, 0) = Real50(2);
    const auto rd = sum_closed(Ld, std::vector<long>{3}, std::vector<double>{1.0},
                               std::vector<Servers>{Servers::of(1), Servers::of(1)},
                               Matrix<double>(2, 1, 1.0), tight());
    const auto rr = sum_closed(Lr, std::vector<long>{3}, std::vector<Real50>{Real50(1)},
                               std::vector<Servers>{Servers::of(1), Servers::of(1)},
                               Matrix<Real50>(2, 1, Real50(1)), tight());
    CHECK(std::fabs(static_cast<double>(rr.XN[0]) - rd.XN[0]) < 1e-12);
    CHECK(std::fabs(static_cast<double>(rr.QN(1, 0)) - rd.QN(1, 0)) < 1e-11);
}
