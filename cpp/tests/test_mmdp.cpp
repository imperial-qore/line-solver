/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * api/mmdp: the Markov-modulated DETERMINISTIC process.
 *
 * The reference's docstrings carry a worked example -- Q = [[-0.5,0.5],
 * [0.3,-0.3]], R = diag(2,5) -- and claim a mean rate of 3.125. THE DOCSTRING IS
 * WRONG AND ITS OWN CODE DISAGREES WITH IT. The stationary law of that generator
 * is (sigma1, sigma0)/(sigma0+sigma1) = (0.375, 0.625), so the mean is
 * 0.375*2 + 0.625*5 = 3.875, which is what native Python actually returns
 * (measured: `mmdp_mean_rate` and `MMDP2.get_mean_rate` both give 3.875). The
 * docstring appears to have dropped the first term. 3.875 is asserted here,
 * together with the two things a test of this module has to separate:
 *
 *  - the SCV is that of the RATE, not of any holding time. Equal rates give SCV
 *    zero however fast the chain switches, and that is what distinguishes this
 *    from an MMPP's interarrival SCV;
 *  - the two-state closed forms and the general stationary solve must agree,
 *    since they are two routes to one answer and the closed form is the one
 *    that would silently drift.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mmdp/mmdp.h"

namespace mmdp = line::mmdp;
using line::Matrix;

namespace {

Matrix<double> gen2(double s0, double s1) {
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -s0;
    Q(0, 1) = s0;
    Q(1, 0) = s1;
    Q(1, 1) = -s1;
    return Q;
}

Matrix<double> diag2(double a, double b) {
    Matrix<double> R(2, 2, 0.0);
    R(0, 0) = a;
    R(1, 1) = b;
    return R;
}

}  // namespace

TEST_CASE("the reference's worked example, at the value its CODE returns") {
    const Matrix<double> Q = gen2(0.5, 0.3), R = diag2(2.0, 5.0);
    CHECK(mmdp::mmdp_isfeasible(Q, R));
    // pi = (sigma1, sigma0)/(sigma0+sigma1) = (0.3, 0.5)/0.8 = (0.375, 0.625),
    // so E[r] = 0.375*2 + 0.625*5 = 3.875. Native Python returns 3.875; only
    // its docstring says 3.125. See the file header.
    CHECK(mmdp::mmdp_mean_rate(Q, R) == doctest::Approx(3.875).epsilon(1e-12));
    CHECK(mmdp::mmdp2_mean_rate(2.0, 5.0, 0.5, 0.3) == doctest::Approx(3.875).epsilon(1e-12));
    // Native Python, same inputs.
    CHECK(mmdp::mmdp_scv(Q, R) == doctest::Approx(0.1404786680541103).epsilon(1e-10));
    CHECK(mmdp::mmdp2_scv(2.0, 5.0, 0.5, 0.3) == doctest::Approx(0.1404786680541103).epsilon(1e-10));
}

TEST_CASE("the closed forms agree with the general stationary solve") {
    const double r0 = 2.0, r1 = 5.0, s0 = 0.5, s1 = 0.3;
    const mmdp::MmdpPair<double> p = mmdp::mmdp2(r0, r1, s0, s1);
    CHECK(mmdp::mmdp_mean_rate(p.Q, p.R) ==
          doctest::Approx(mmdp::mmdp2_mean_rate(r0, r1, s0, s1)).epsilon(1e-12));
    CHECK(mmdp::mmdp_scv(p.Q, p.R) ==
          doctest::Approx(mmdp::mmdp2_scv(r0, r1, s0, s1)).epsilon(1e-10));

    // The chain spends LONGER in the state it leaves more slowly, so raising
    // sigma0 (the rate OUT of state 0) shifts the mean toward r1.
    const double a = mmdp::mmdp2_mean_rate(r0, r1, 0.1, 0.3);
    const double b = mmdp::mmdp2_mean_rate(r0, r1, 2.0, 0.3);
    CHECK(a < b);
    CHECK(a > r0);
    CHECK(b < r1);
}

TEST_CASE("the SCV is the RATE's dispersion, not a holding time's") {
    // Equal rates: no variability at all, however fast the chain switches.
    for (double s = 0.01; s < 100.0; s *= 10.0) {
        const Matrix<double> Q = gen2(s, 2.0 * s), R = diag2(3.0, 3.0);
        CHECK(mmdp::mmdp_scv(Q, R) == doctest::Approx(0.0).epsilon(1e-12));
    }
    // A one-state process likewise.
    Matrix<double> Q1(1, 1, 0.0), R1(1, 1, 7.0);
    CHECK(mmdp::mmdp_scv(Q1, R1) == doctest::Approx(0.0));
    CHECK(mmdp::mmdp_mean_rate(Q1, R1) == doctest::Approx(7.0));

    // Unequal rates: the SCV is Var/mean^2 under the stationary law.
    const Matrix<double> Q = gen2(0.5, 0.3), R = diag2(2.0, 5.0);
    const double p0 = 0.3 / 0.8, p1 = 0.5 / 0.8;
    const double m = p0 * 2.0 + p1 * 5.0;
    const double v = p0 * 4.0 + p1 * 25.0 - m * m;
    CHECK(mmdp::mmdp_scv(Q, R) == doctest::Approx(v / (m * m)).epsilon(1e-10));
    CHECK(mmdp::mmdp_scv(Q, R) > 0.0);
}

TEST_CASE("a zero mean rate reports infinity rather than dividing") {
    // A caller can test infinity; a NaN it cannot.
    const Matrix<double> Q = gen2(0.5, 0.3), R = diag2(0.0, 0.0);
    CHECK(std::isinf(mmdp::mmdp_scv(Q, R)));
    CHECK(std::isinf(mmdp::mmdp2_scv(0.0, 0.0, 0.5, 0.3)));
}

TEST_CASE("feasibility rejects what is not an MMDP") {
    CHECK(mmdp::mmdp_isfeasible(gen2(0.5, 0.3), diag2(2.0, 5.0)));

    // Rows that do not sum to zero are not a generator.
    Matrix<double> bad = gen2(0.5, 0.3);
    bad(0, 1) = 0.7;
    CHECK(!mmdp::mmdp_isfeasible(bad, diag2(2.0, 5.0)));

    // A negative off-diagonal is not a rate.
    Matrix<double> neg(2, 2, 0.0);
    neg(0, 0) = 0.5;
    neg(0, 1) = -0.5;
    neg(1, 0) = 0.3;
    neg(1, 1) = -0.3;
    CHECK(!mmdp::mmdp_isfeasible(neg, diag2(2.0, 5.0)));

    // A NON-DIAGONAL R is refused, and that is a modelling test rather than a
    // storage one: it would make the flow rate depend on a transition.
    Matrix<double> offdiag = diag2(2.0, 5.0);
    offdiag(0, 1) = 0.1;
    CHECK(!mmdp::mmdp_isfeasible(gen2(0.5, 0.3), offdiag));

    // A negative rate is not a flow.
    CHECK(!mmdp::mmdp_isfeasible(gen2(0.5, 0.3), diag2(-1.0, 5.0)));

    // Mismatched sizes.
    Matrix<double> R3(3, 3, 1.0);
    CHECK(!mmdp::mmdp_isfeasible(gen2(0.5, 0.3), R3));
}

TEST_CASE("the MMDP of a MAP keeps the phase process and takes D1's row sums") {
    Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = 0.5;
    D0(1, 0) = 0.3;
    D0(1, 1) = -1.5;
    D1(0, 0) = 1.0;
    D1(0, 1) = 0.5;
    D1(1, 0) = 0.4;
    D1(1, 1) = 0.8;

    const mmdp::MmdpPair<double> p = mmdp::mmdp_from_map(D0, D1);
    // Q is the MAP's full generator, so its rows sum to zero.
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) {
            s += p.Q(i, j);
            CHECK(p.Q(i, j) == doctest::Approx(D0(i, j) + D1(i, j)));
        }
        CHECK(s == doctest::Approx(0.0).epsilon(1e-12));
    }
    // R holds D1's row sums: the rate at which that phase generates arrivals.
    CHECK(p.R(0, 0) == doctest::Approx(1.5));
    CHECK(p.R(1, 1) == doctest::Approx(1.2));
    CHECK(p.R(0, 1) == doctest::Approx(0.0));
    CHECK(mmdp::mmdp_isfeasible(p.Q, p.R));

    // The mean rate of that MMDP is the MAP's own arrival rate, since both are
    // pi times the D1 row sums.
    CHECK(mmdp::mmdp_mean_rate(p.Q, p.R) > 0.0);
}

TEST_CASE("the refusals are by name") {
    Matrix<double> empty(0, 0, 0.0);
    CHECK_THROWS_AS(mmdp::mmdp_mean_rate(empty, empty), line::InputError);
    CHECK_THROWS_AS(mmdp::mmdp_scv(empty, empty), line::InputError);
    CHECK_THROWS_AS(mmdp::mmdp_mean_rate(gen2(0.5, 0.3), Matrix<double>(3, 3, 1.0)),
                    line::InputError);
    CHECK_THROWS_AS(mmdp::mmdp2(2.0, 5.0, 0.0, 0.3), line::InputError);
    CHECK_THROWS_AS(mmdp::mmdp2(2.0, 5.0, 0.5, -1.0), line::InputError);
    CHECK_THROWS_AS(mmdp::mmdp_from_map(Matrix<double>(2, 3, 0.0), Matrix<double>(2, 3, 0.0)),
                    line::InputError);
}
