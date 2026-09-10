/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Non-product-form traffic approximations. Oracles, in order of strength:
 *   1. Conservation identities the traffic algebra must satisfy as identities,
 *      checked for exact equality in the rational instantiation: the split
 *      flows sum back to the flow that was split, the total arrival rates
 *      satisfy the traffic equations, and Xi is a true inverse.
 *   2. Fixed points of the variability equations that can be read off by hand
 *      (a Poisson feed-forward network has c2a = 1 everywhere).
 *   3. Closed-form values of w*(t) at its two limits.
 * Everything that stays in the field also runs at T = Rational and must agree
 * with double to rounding.
 */
#include <cmath>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_nonexp_approx.h"
#include "line/api/npfqn/npfqn_rqna_weight.h"
#include "line/api/npfqn/npfqn_traffic_idc.h"
#include "line/api/npfqn/npfqn_traffic_split_cs.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::npfqn::Mmap;
using line::npfqn::npfqn_nonexp_approx;
using line::npfqn::npfqn_rqna_weight;
using line::npfqn::npfqn_traffic_idc;
using line::npfqn::npfqn_traffic_idc_at;
using line::npfqn::npfqn_traffic_split_cs;

namespace {

constexpr double TOL = 1e-9;

/** A two-phase MMAP with two marked classes, exact in every arithmetic. */
template <class T>
Mmap<T> mmap_2class() {
    using nt = line::num_traits<T>;
    Mmap<T> M(4, Matrix<T>(2, 2, nt::from_int(0)));
    // D1^(1) and D1^(2), then D1 = sum, then D0 with zero row sums
    M[2](0, 0) = nt::from_rational(1, 2);
    M[2](0, 1) = nt::from_rational(1, 4);
    M[2](1, 1) = nt::from_rational(1, 5);
    M[3](0, 1) = nt::from_rational(1, 10);
    M[3](1, 0) = nt::from_rational(3, 10);
    M[3](1, 1) = nt::from_rational(1, 20);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) M[1](i, j) = M[2](i, j) + M[3](i, j);
    M[0](0, 1) = nt::from_rational(1, 3);
    M[0](1, 0) = nt::from_rational(1, 6);
    for (std::size_t k = 0; k < 2; ++k) {
        T s = nt::from_int(0);
        for (std::size_t j = 0; j < 2; ++j) s += M[0](k, j) + M[1](k, j);
        M[0](k, k) = -s;
    }
    return M;
}

/** Routing with class switching to two destinations, rows summing to one. */
template <class T>
Matrix<T> split_probs() {
    using nt = line::num_traits<T>;
    Matrix<T> P(2, 4, nt::from_int(0));
    // destination 1: columns 0,1 ; destination 2: columns 2,3
    P(0, 0) = nt::from_rational(1, 4);
    P(0, 1) = nt::from_rational(1, 4);
    P(0, 2) = nt::from_rational(1, 8);
    P(0, 3) = nt::from_rational(3, 8);
    P(1, 0) = nt::from_rational(1, 5);
    P(1, 1) = nt::from_rational(1, 10);
    P(1, 2) = nt::from_rational(1, 2);
    P(1, 3) = nt::from_rational(1, 5);
    return P;
}

/** Three-queue open tandem with feedback, used for the IDC equations. */
template <class T>
Matrix<T> routing_3() {
    using nt = line::num_traits<T>;
    Matrix<T> P(3, 3, nt::from_int(0));
    P(0, 1) = nt::from_rational(1, 2);
    P(0, 2) = nt::from_rational(1, 4);
    P(1, 2) = nt::from_rational(1, 2);
    P(2, 0) = nt::from_rational(1, 5);
    return P;
}

template <class T>
line::npfqn::TrafficIdcContext<T> idc_3() {
    using nt = line::num_traits<T>;
    const std::vector<T> lambda0{nt::from_rational(1, 2), nt::from_rational(1, 10),
                                 nt::from_rational(1, 20)};
    const std::vector<T> c2a0{nt::from_int(1), nt::from_rational(3, 2), nt::from_rational(1, 2)};
    const std::vector<T> mu{nt::from_int(2), nt::from_int(3), nt::from_int(4)};
    const std::vector<T> cs2{nt::from_int(1), nt::from_rational(1, 4), nt::from_int(2)};
    return npfqn_traffic_idc(lambda0, routing_3<T>(), c2a0, mu, cs2);
}

}  // namespace

// ---------------------------------------------------------------------------
// npfqn_traffic_split_cs
// ---------------------------------------------------------------------------

TEST_CASE("npfqn_traffic_split_cs conserves the split flow exactly") {
    const Mmap<Rational> M = mmap_2class<Rational>();
    const Matrix<Rational> P = split_probs<Rational>();
    const std::vector<Mmap<Rational>> S = npfqn_traffic_split_cs(M, P);
    REQUIRE(S.size() == 2);

    // The routing rows sum to one, so summing the destination arrival matrices
    // must return the arrival matrix that was split, entry by entry and with
    // no rounding at all.
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 2; ++b) {
            Rational sum = 0;
            for (std::size_t j = 0; j < 2; ++j) sum += S[j][1](a, b);
            CHECK(sum == M[1](a, b));
        }

    // Per class as well: sum_j sum over destinations of D1^(s)_j equals
    // sum_r D1^(r) weighted by the class-switching probabilities, whose column
    // sums over (j,s) are one for every r.
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 2; ++b) {
            Rational sum = 0;
            for (std::size_t j = 0; j < 2; ++j)
                for (std::size_t s = 0; s < 2; ++s) sum += S[j][2 + s](a, b);
            CHECK(sum == M[1](a, b));
        }
}

TEST_CASE("npfqn_traffic_split_cs returns proper MMAPs") {
    const Mmap<Rational> M = mmap_2class<Rational>();
    const std::vector<Mmap<Rational>> S = npfqn_traffic_split_cs(M, split_probs<Rational>());
    for (const Mmap<Rational>& s : S) {
        REQUIRE(s.size() == 4);
        for (std::size_t a = 0; a < 2; ++a) {
            // D0 + D1 is a generator: every row sums to zero exactly
            Rational rowsum = 0;
            for (std::size_t b = 0; b < 2; ++b) rowsum += s[0](a, b) + s[1](a, b);
            CHECK(rowsum == 0);
            // D1 is the sum of the markings
            for (std::size_t b = 0; b < 2; ++b) CHECK(s[1](a, b) == s[2](a, b) + s[3](a, b));
            // off-diagonal of D0 is nonnegative, diagonal is not
            for (std::size_t b = 0; b < 2; ++b)
                if (a != b) CHECK(s[0](a, b) >= 0);
            CHECK(s[0](a, a) <= 0);
        }
    }
}

TEST_CASE("npfqn_traffic_split_cs matches by hand on a one-class one-phase flow") {
    // D0 = -1, D1 = D1^(1) = 1, split 1/4 to destination 1 and 3/4 to 2.
    Mmap<Rational> M(3, Matrix<Rational>(1, 1, Rational(0)));
    M[0](0, 0) = -1;
    M[1](0, 0) = 1;
    M[2](0, 0) = 1;
    Matrix<Rational> P(1, 2, Rational(0));
    P(0, 0) = Rational(1, 4);
    P(0, 1) = Rational(3, 4);

    const std::vector<Mmap<Rational>> S = npfqn_traffic_split_cs(M, P);
    CHECK(S[0][1](0, 0) == Rational(1, 4));
    CHECK(S[0][0](0, 0) == Rational(-1, 4));
    CHECK(S[1][1](0, 0) == Rational(3, 4));
    CHECK(S[1][0](0, 0) == Rational(-3, 4));
}

TEST_CASE("npfqn_traffic_split_cs exact equals double to rounding") {
    const std::vector<Mmap<double>> d = npfqn_traffic_split_cs(mmap_2class<double>(), split_probs<double>());
    const std::vector<Mmap<Rational>> q =
        npfqn_traffic_split_cs(mmap_2class<Rational>(), split_probs<Rational>());
    for (std::size_t j = 0; j < 2; ++j)
        for (std::size_t c = 0; c < 4; ++c)
            for (std::size_t a = 0; a < 2; ++a)
                for (std::size_t b = 0; b < 2; ++b)
                    CHECK(static_cast<double>(q[j][c](a, b)) ==
                          doctest::Approx(d[j][c](a, b)).epsilon(TOL));
}

TEST_CASE("npfqn_traffic_split_cs rejects a malformed class-switching matrix") {
    const Mmap<double> M = mmap_2class<double>();
    CHECK_THROWS_AS(npfqn_traffic_split_cs(M, Matrix<double>(2, 3, 0.0)), line::InputError);
}

// ---------------------------------------------------------------------------
// npfqn_rqna_weight
// ---------------------------------------------------------------------------

TEST_CASE("npfqn_rqna_weight hits its two limits and stays monotone in [0,1]") {
    CHECK(npfqn_rqna_weight(0.0) == 0.0);
    CHECK(npfqn_rqna_weight(-1.0) == 0.0);
    CHECK(npfqn_rqna_weight(std::numeric_limits<double>::infinity()) == 1.0);
    CHECK(npfqn_rqna_weight(1e-9) == 0.0);  // series branch

    double prev = 0.0;
    for (double t = 0.01; t < 60.0; t *= 1.5) {
        const double w = npfqn_rqna_weight(t);
        CHECK(w >= prev);
        CHECK(w <= 1.0);
        prev = w;
    }
    CHECK(npfqn_rqna_weight(1e4) == doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("npfqn_rqna_weight agrees between double and high precision") {
    for (double t : {0.25, 1.0, 2.5, 10.0}) {
        const Real50 wr = npfqn_rqna_weight(Real50(t));
        CHECK(static_cast<double>(wr) == doctest::Approx(npfqn_rqna_weight(t)).epsilon(TOL));
    }
}

TEST_CASE("npfqn_rqna_weight vector form is the elementwise scalar form") {
    const std::vector<double> t{0.0, 0.5, 3.0};
    const std::vector<double> w = npfqn_rqna_weight(t);
    REQUIRE(w.size() == 3);
    for (std::size_t i = 0; i < 3; ++i) CHECK(w[i] == npfqn_rqna_weight(t[i]));
}

// ---------------------------------------------------------------------------
// npfqn_traffic_idc
// ---------------------------------------------------------------------------

TEST_CASE("npfqn_traffic_idc solves the traffic rate equations exactly") {
    const line::npfqn::TrafficIdcContext<Rational> ctx = idc_3<Rational>();
    const Matrix<Rational> P = routing_3<Rational>();
    // lambda = lambda0 + P' lambda, as an identity in exact arithmetic
    for (std::size_t i = 0; i < 3; ++i) {
        Rational rhs = ctx.lambda0[i];
        for (std::size_t j = 0; j < 3; ++j) rhs += P(j, i) * ctx.lambda[j];
        CHECK(ctx.lambda[i] == rhs);
        CHECK(ctx.rho[i] == ctx.lambda[i] / ctx.mu[i]);
    }
    // Xi is a true inverse of I - P'
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            Rational s = 0;
            for (std::size_t k = 0; k < 3; ++k) s += ((i == k ? Rational(1) : Rational(0)) - P(k, i)) * ctx.Xi(k, j);
            CHECK(s == (i == j ? Rational(1) : Rational(0)));
        }
}

TEST_CASE("npfqn_traffic_idc satisfies its own variability equations exactly") {
    const line::npfqn::TrafficIdcContext<Rational> ctx = idc_3<Rational>();
    const Matrix<Rational> P = routing_3<Rational>();
    for (std::size_t i = 0; i < 3; ++i) {
        // c2d_i = c2a_i, and c2x_i = c2a_i + cs2_i
        CHECK(ctx.c2d[i] == ctx.c2a[i]);
        CHECK(ctx.c2x[i] == ctx.c2a[i] + ctx.cs2[i]);
        // c2aij_{i,j} = p_{i,j} c2d_i + (1 - p_{i,j}) + c2alpha_{i,j}
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(ctx.c2aij(i, j) == P(i, j) * ctx.c2d[i] + (Rational(1) - P(i, j)) + ctx.c2alpha(i, j));
    }
}

TEST_CASE("npfqn_traffic_idc reproduces Poisson superposition on a feed-forward network") {
    // No internal routing at all: every arrival stream is the external one, so
    // the asymptotic variability parameter is exactly the external SCV and the
    // corrections vanish (P = 0 makes both alpha and beta zero).
    using nt = line::num_traits<Rational>;
    const std::vector<Rational> lambda0{nt::from_int(1), nt::from_rational(1, 2)};
    const std::vector<Rational> c2a0{nt::from_int(1), nt::from_int(3)};
    const std::vector<Rational> mu{nt::from_int(4), nt::from_int(4)};
    const std::vector<Rational> cs2{nt::from_int(1), nt::from_int(1)};
    const Matrix<Rational> P(2, 2, Rational(0));
    const line::npfqn::TrafficIdcContext<Rational> ctx = npfqn_traffic_idc(lambda0, P, c2a0, mu, cs2);
    CHECK(ctx.c2a[0] == Rational(1));
    CHECK(ctx.c2a[1] == Rational(3));
    CHECK(ctx.c2d[0] == Rational(1));
    CHECK(ctx.c2x[1] == Rational(4));
}

TEST_CASE("npfqn_traffic_idc exact equals double to rounding") {
    const line::npfqn::TrafficIdcContext<double> d = idc_3<double>();
    const line::npfqn::TrafficIdcContext<Rational> q = idc_3<Rational>();
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(static_cast<double>(q.lambda[i]) == doctest::Approx(d.lambda[i]).epsilon(TOL));
        CHECK(static_cast<double>(q.c2a[i]) == doctest::Approx(d.c2a[i]).epsilon(TOL));
        CHECK(static_cast<double>(q.c2d[i]) == doctest::Approx(d.c2d[i]).epsilon(TOL));
        CHECK(static_cast<double>(q.c2x[i]) == doctest::Approx(d.c2x[i]).epsilon(TOL));
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(static_cast<double>(q.c2aij(i, j)) == doctest::Approx(d.c2aij(i, j)).epsilon(TOL));
    }
}

TEST_CASE("npfqn_traffic_idc_at converges to the limiting variability parameters") {
    const line::npfqn::TrafficIdcContext<double> ctx = idc_3<double>();
    // Renewal inputs: the IDC of a renewal process is its SCV at every scale,
    // so the time-dependent system must reduce to the limiting one as w -> 1.
    const std::function<std::vector<double>(const double&)> a0 = [&](const double&) {
        return ctx.c2a0;
    };
    const std::function<std::vector<double>(const std::vector<double>&)> sfun =
        [&](const std::vector<double>&) { return ctx.cs2; };

    const std::vector<double> near = npfqn_traffic_idc_at(ctx, 1e7, a0, sfun);
    for (std::size_t i = 0; i < 3; ++i) CHECK(near[i] == doctest::Approx(ctx.c2a[i]).epsilon(1e-4));

    // At t -> 0 the weights vanish, the departure IDC collapses onto the
    // service IDC and the answer stays finite and positive.
    const std::vector<double> early = npfqn_traffic_idc_at(ctx, 1e-9, a0, sfun);
    for (std::size_t i = 0; i < 3; ++i) CHECK(early[i] > 0.0);
}

TEST_CASE("npfqn_traffic_idc corrections can be switched off independently") {
    using nt = line::num_traits<Rational>;
    const std::vector<Rational> lambda0{nt::from_rational(1, 2), nt::from_rational(1, 10),
                                        nt::from_rational(1, 20)};
    const std::vector<Rational> c2a0{nt::from_int(1), nt::from_rational(3, 2), nt::from_rational(1, 2)};
    const std::vector<Rational> mu{nt::from_int(2), nt::from_int(3), nt::from_int(4)};
    const std::vector<Rational> cs2{nt::from_int(1), nt::from_rational(1, 4), nt::from_int(2)};
    line::npfqn::TrafficIdcCorrections off;
    off.alpha = false;
    off.beta = false;
    const line::npfqn::TrafficIdcContext<Rational> plain =
        npfqn_traffic_idc(lambda0, routing_3<Rational>(), c2a0, mu, cs2, off);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK(plain.c2alpha(i, j) == 0);
            CHECK(plain.zetaAll[i](0, j) == 0);
        }
    // Without corrections, c2aij_{i,j} = p_{i,j} c2d_i + (1 - p_{i,j}) exactly.
    const Matrix<Rational> P = routing_3<Rational>();
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(plain.c2aij(i, j) == P(i, j) * plain.c2d[i] + (Rational(1) - P(i, j)));
}

// ---------------------------------------------------------------------------
// npfqn_nonexp_approx
// ---------------------------------------------------------------------------

namespace {

struct NonexpFixture {
    Matrix<double> ST{{0.5, 0.2}, {0.4, 0.4}};
    Matrix<double> V{{1.0, 1.0}, {1.0, 1.0}};
    Matrix<double> SCV{{4.0, 0.5}, {1.0, 1.0}};
    Matrix<double> T{{1.0, 2.0}, {1.5, 0.5}};
    Matrix<double> U{{0.5, 0.4}, {0.6, 0.2}};
    Matrix<double> rates{{2.0, 5.0}, {2.5, 2.5}};
    std::vector<bool> isFCFS{true, true};
    std::vector<double> gamma{0.0, 0.0};
    std::vector<double> nservers{1.0, 1.0};
};

}  // namespace

TEST_CASE("npfqn_nonexp_approx no-op methods return their inputs") {
    NonexpFixture f;
    for (const char* m : {"default", "none", "hvmva"}) {
        const auto r = npfqn_nonexp_approx<double>(m, f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T, f.U,
                                                   f.gamma, f.nservers);
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(r.nservers[i] == f.nservers[i]);
            CHECK(r.rho[i] == 0.0);
            CHECK(r.scvs[i] == 1.0);
            for (std::size_t k = 0; k < 2; ++k) CHECK(r.ST(i, k) == f.ST(i, k));
        }
    }
    CHECK_THROWS_AS(npfqn_nonexp_approx<double>("nosuch", f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T,
                                                f.U, f.gamma, f.nservers),
                    line::InputError);
}

TEST_CASE("npfqn_nonexp_approx interp rescales only the non-product-form station") {
    NonexpFixture f;
    // Station 2 is product form under this test: both classes have SCV 1 and
    // the same demand, so its row must come back untouched with one server.
    f.ST(1, 0) = 0.4;
    f.ST(1, 1) = 0.4;
    const auto r = npfqn_nonexp_approx<double>("interp", f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T,
                                               f.U, f.gamma, f.nservers);
    CHECK(r.ST(1, 0) == f.ST(1, 0));
    CHECK(r.ST(1, 1) == f.ST(1, 1));
    CHECK(r.scvs[1] == 1.0);

    // Station 1 mixes SCV 4 and SCV 0.5, so it is rescaled.
    CHECK(r.rho[0] == doctest::Approx(0.9).epsilon(TOL));
    // scvs = (SCV . T)/sum(T) = (4*1 + 0.5*2)/3
    CHECK(r.scvs[0] == doctest::Approx(5.0 / 3.0).epsilon(TOL));
    CHECK(r.scva[0] == 1.0);
    // gamma = (rho^c + rho)/2 with c = 1, i.e. rho
    CHECK(r.gamma[0] == doctest::Approx(0.9).epsilon(TOL));
    // eta = exp(-2(1-rho)/(scvs + scva rho))
    const double eta = std::exp(-2.0 * (1.0 - 0.9) / (5.0 / 3.0 + 0.9));
    CHECK(r.eta[0] == doctest::Approx(eta).epsilon(TOL));
    // ST(1,k) = (1-a) ST + a (b eta + (1-b) gamma) (c / sum T), a = b = rho^8
    const double a = std::pow(0.9, 8);
    for (std::size_t k = 0; k < 2; ++k) {
        const double want = (1 - a) * f.ST(0, k) + a * (a * eta + (1 - a) * 0.9) * (1.0 / 3.0);
        CHECK(r.ST(0, k) == doctest::Approx(want).epsilon(TOL));
    }
    CHECK(r.nservers[0] == 1.0);
}

TEST_CASE("npfqn_nonexp_approx leaves non-FCFS stations alone") {
    NonexpFixture f;
    f.isFCFS[0] = false;
    const auto r = npfqn_nonexp_approx<double>("interp", f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T,
                                               f.U, f.gamma, f.nservers);
    for (std::size_t k = 0; k < 2; ++k) CHECK(r.ST(0, k) == f.ST(0, k));
    // rho is still accumulated for every station, FCFS or not
    CHECK(r.rho[0] == doctest::Approx(0.9).epsilon(TOL));
}

TEST_CASE("npfqn_nonexp_approx uses the M/M/1 branch at unit SCV with one server") {
    NonexpFixture f;
    // Unequal demands but both SCVs at one: non-product-form by the demand
    // test, and eta falls back to rho rather than the diffusion formula.
    f.SCV(0, 0) = 1.0;
    f.SCV(0, 1) = 1.0;
    const auto r = npfqn_nonexp_approx<double>("interp", f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T,
                                               f.U, f.gamma, f.nservers);
    CHECK(r.scvs[0] == doctest::Approx(1.0).epsilon(TOL));
    CHECK(r.eta[0] == doctest::Approx(r.rho[0]).epsilon(TOL));
}

TEST_CASE("npfqn_nonexp_approx agrees between double and high precision") {
    NonexpFixture f;
    const auto d = npfqn_nonexp_approx<double>("interp", f.isFCFS, f.rates, f.ST, f.V, f.SCV, f.T,
                                               f.U, f.gamma, f.nservers);
    Matrix<Real50> ST(2, 2), V(2, 2), SCV(2, 2), T(2, 2), U(2, 2), rates(2, 2);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
            ST(i, k) = Real50(f.ST(i, k));
            V(i, k) = Real50(f.V(i, k));
            SCV(i, k) = Real50(f.SCV(i, k));
            T(i, k) = Real50(f.T(i, k));
            U(i, k) = Real50(f.U(i, k));
            rates(i, k) = Real50(f.rates(i, k));
        }
    const std::vector<Real50> gamma(2, Real50(0)), nservers(2, Real50(1));
    const auto r = npfqn_nonexp_approx<Real50>("interp", f.isFCFS, rates, ST, V, SCV, T, U, gamma,
                                               nservers);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(static_cast<double>(r.eta[i]) == doctest::Approx(d.eta[i]).epsilon(TOL));
        CHECK(static_cast<double>(r.scvs[i]) == doctest::Approx(d.scvs[i]).epsilon(TOL));
        for (std::size_t k = 0; k < 2; ++k)
            CHECK(static_cast<double>(r.ST(i, k)) == doctest::Approx(d.ST(i, k)).epsilon(TOL));
    }
}
