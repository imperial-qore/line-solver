/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Covers api/mc/ctmc_isfeasible.h, ctmc_relsolve.h, ctmc_timereverse.h and
 * dtmc_transient.h. Every literal here is reproducible in closed form on the
 * two three-state chains used throughout, so the test states the closed form
 * next to the number rather than only the number.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_isfeasible.h"
#include "line/api/mc/ctmc_relsolve.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_timereverse.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/dtmc_transient.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"

using line::Matrix;
using namespace line::mc;

namespace {

Matrix<double> mk(std::size_t n, const std::vector<double>& v) {
    Matrix<double> A(n, n);
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = v[k++];
    return A;
}

/** Stationary measure (1, 1, 4), so the relative solution is exactly integral. */
Matrix<double> ref_generator() { return mk(3, {-3, 2, 1, 1, -4, 3, 0.5, 0.5, -1}); }

Matrix<double> ref_stochastic() { return mk(3, {0.2, 0.5, 0.3, 0.1, 0.6, 0.3, 0.4, 0.4, 0.2}); }

}  // namespace

TEST_CASE("ctmc_relsolve normalises on the reference state, not on the sum") {
    const Matrix<double> Q = ref_generator();
    const std::vector<double> p = ctmc_relsolve(Q);
    REQUIRE(p.size() == 3);
    CHECK(p[0] == doctest::Approx(1.0).epsilon(1e-13));
    CHECK(p[1] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(p[2] == doctest::Approx(4.0).epsilon(1e-12));
    // The measure still solves the balance equations
    const std::vector<double> r = line::vecmul(p, Q);
    for (std::size_t j = 0; j < 3; ++j) CHECK(r[j] == doctest::Approx(0.0).epsilon(1e-11).scale(1.0));
    // Renormalising recovers ctmc_solve
    const std::vector<double> pi = ctmc_solve(Q);
    double s = p[0] + p[1] + p[2];
    for (std::size_t j = 0; j < 3; ++j) CHECK(p[j] / s == doctest::Approx(pi[j]).epsilon(1e-11));
    // A different reference state rescales the same measure
    const std::vector<double> p2 = ctmc_relsolve(Q, static_cast<std::size_t>(2));
    CHECK(p2[2] == doctest::Approx(1.0).epsilon(1e-12));
    for (std::size_t j = 0; j < 3; ++j) CHECK(p2[j] * 4.0 == doctest::Approx(p[j]).epsilon(1e-11));
}

TEST_CASE("ctmc_relsolve renormalises globally on a reducible generator") {
    // Two disconnected two-state chains: the reference-state closure is
    // meaningless across components, so the reference returns a probability
    // vector there instead of a relative measure
    Matrix<double> Q(4, 4, 0.0);
    Q(0, 1) = 1.0;
    Q(0, 0) = -1.0;
    Q(1, 0) = 3.0;
    Q(1, 1) = -3.0;
    Q(2, 3) = 2.0;
    Q(2, 2) = -2.0;
    Q(3, 2) = 2.0;
    Q(3, 3) = -2.0;
    const std::vector<double> p = ctmc_relsolve(Q);
    double s = 0;
    for (std::size_t i = 0; i < 4; ++i) s += p[i];
    CHECK(s == doctest::Approx(1.0).epsilon(1e-11));
    // Within a component the ratio is the component's own stationary law
    CHECK(p[0] / p[1] == doctest::Approx(3.0).epsilon(1e-11));
    CHECK(p[2] / p[3] == doctest::Approx(1.0).epsilon(1e-11));
}

TEST_CASE("ctmc_timereverse and dtmc_timereverse are involutions on the stationary law") {
    const Matrix<double> Q = ref_generator();
    const Matrix<double> R = ctmc_timereverse(Q);
    // Qrev(i,j) = Q(j,i) pi(j) / pi(i) with pi proportional to (1, 1, 4)
    CHECK(R(0, 0) == doctest::Approx(-3.0).epsilon(1e-11));
    CHECK(R(0, 1) == doctest::Approx(1.0).epsilon(1e-11));
    CHECK(R(0, 2) == doctest::Approx(2.0).epsilon(1e-11));
    CHECK(R(2, 0) == doctest::Approx(0.25).epsilon(1e-11));
    CHECK(R(2, 1) == doctest::Approx(0.75).epsilon(1e-11));
    // The reversed chain is a generator with the same stationary law
    for (std::size_t i = 0; i < 3; ++i) {
        double row = 0;
        for (std::size_t j = 0; j < 3; ++j) row += R(i, j);
        CHECK(row == doctest::Approx(0.0).epsilon(1e-11).scale(1.0));
    }
    const std::vector<double> pi = ctmc_solve(Q), pr = ctmc_solve(R);
    for (std::size_t j = 0; j < 3; ++j) CHECK(pr[j] == doctest::Approx(pi[j]).epsilon(1e-10));
    const Matrix<double> RR = ctmc_timereverse(R);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(RR(i, j) == doctest::Approx(Q(i, j)).epsilon(1e-9));

    const Matrix<double> P = ref_stochastic();
    const Matrix<double> S = dtmc_timereverse(P);
    CHECK(S(0, 0) == doctest::Approx(0.2).epsilon(1e-11));
    CHECK(S(0, 1) == doctest::Approx(0.26).epsilon(1e-11));
    CHECK(S(0, 2) == doctest::Approx(0.54).epsilon(1e-11));
    CHECK(S(1, 0) == doctest::Approx(0.19230769230769).epsilon(1e-11));
    CHECK(S(2, 0) == doctest::Approx(0.22222222222222).epsilon(1e-11));
    for (std::size_t i = 0; i < 3; ++i) {
        double row = 0;
        for (std::size_t j = 0; j < 3; ++j) row += S(i, j);
        CHECK(row == doctest::Approx(1.0).epsilon(1e-11));
    }
}

TEST_CASE("ctmc_isfeasible and dtmc_isfeasible separate generators from kernels") {
    const Matrix<double> Q = ref_generator(), P = ref_stochastic();
    CHECK(ctmc_isfeasible(Q));
    CHECK_FALSE(ctmc_isfeasible(P));
    // dtmc_isfeasible returns a precision level, not a flag
    CHECK(dtmc_isfeasible(P) == 15);
    CHECK(dtmc_isfeasible(Q) == 0);
    // A row that misses one by 1e-7 is feasible only to six digits
    Matrix<double> Pl = P;
    Pl(0, 0) += 1e-7;
    CHECK(dtmc_isfeasible(Pl) == 6);
    // A negative off-diagonal is not a generator
    Matrix<double> Qb = Q;
    Qb(0, 1) = -2.0;
    Qb(0, 0) = 1.0;
    CHECK_FALSE(ctmc_isfeasible(Qb));
}

TEST_CASE("dtmc_transient walks the chain and keeps every row a distribution") {
    const Matrix<double> P = ref_stochastic();
    const Matrix<double> pit = dtmc_transient(P, std::vector<double>{1.0, 0.0, 0.0}, 4);
    REQUIRE(pit.rows() == 5);
    REQUIRE(pit.cols() == 3);
    // Row zero is the initial condition, not the first step
    CHECK(pit(0, 0) == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(pit(1, 0) == doctest::Approx(0.2).epsilon(1e-13));
    CHECK(pit(1, 1) == doctest::Approx(0.5).epsilon(1e-13));
    CHECK(pit(2, 0) == doctest::Approx(0.21).epsilon(1e-12));
    CHECK(pit(4, 0) == doctest::Approx(0.2021).epsilon(1e-11));
    CHECK(pit(4, 1) == doctest::Approx(0.5252).epsilon(1e-11));
    CHECK(pit(4, 2) == doctest::Approx(0.2727).epsilon(1e-11));
    for (std::size_t k = 0; k < 5; ++k) {
        double s = 0;
        for (std::size_t j = 0; j < 3; ++j) s += pit(k, j);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
    // The walk converges to the stationary law
    const Matrix<double> longrun = dtmc_transient(P, std::vector<double>{1.0, 0.0, 0.0}, 200);
    const std::vector<double> pi = dtmc_solve(P);
    for (std::size_t j = 0; j < 3; ++j)
        CHECK(longrun(200, j) == doctest::Approx(pi[j]).epsilon(1e-10));
}

TEST_CASE("dtmc_hitting_time solves the first-passage system") {
    const Matrix<double> P = ref_stochastic();
    const std::vector<double> h = dtmc_hitting_time(P, std::vector<std::size_t>{2});
    REQUIRE(h.size() == 3);
    // Both non-target states reach state 3 with probability 0.3 in one step,
    // so the mean number of steps is 1/0.3 from either
    CHECK(h[0] == doctest::Approx(10.0 / 3.0).epsilon(1e-11));
    CHECK(h[1] == doctest::Approx(10.0 / 3.0).epsilon(1e-11));
    CHECK(h[2] == doctest::Approx(0.0).epsilon(1e-14));
    // A state that cannot reach the target has infinite hitting time
    Matrix<double> Pa(3, 3, 0.0);
    Pa(0, 0) = 1.0;
    Pa(1, 2) = 1.0;
    Pa(2, 2) = 1.0;
    const std::vector<double> ha = dtmc_hitting_time(Pa, std::vector<std::size_t>{2});
    CHECK(ha[1] == doctest::Approx(1.0).epsilon(1e-13));
    CHECK(ha[2] == doctest::Approx(0.0).epsilon(1e-14));
    CHECK(!std::isfinite(ha[0]));
}

TEST_CASE("dtmc_uniformization is the matrix exponential of the uniformized generator") {
    const Matrix<double> P = ref_stochastic();
    const std::vector<double> pi0{1.0, 0.0, 0.0};
    const UniformizationResult<double> u = dtmc_uniformization(pi0, P, 2.0, 1e-12, -1);
    REQUIRE(u.pi.size() == 3);
    CHECK(u.pi[0] == doctest::Approx(0.3116316677797).epsilon(1e-10));
    CHECK(u.pi[1] == doctest::Approx(0.44586010268219).epsilon(1e-10));
    CHECK(u.pi[2] == doctest::Approx(0.24250822953733).epsilon(1e-10));
    double s = u.pi[0] + u.pi[1] + u.pi[2];
    CHECK(s == doctest::Approx(1.0).epsilon(1e-11));
    // The reference reads the DTMC through ctmc_makeinfgen, so the answer must
    // be pi0 * expm(Q t) with Q that generator, NOT pi0 * P^2
    Matrix<double> Qe = line::expm(ctmc_makeinfgen(P), 2.0);
    const std::vector<double> direct = line::vecmul(pi0, Qe);
    for (std::size_t j = 0; j < 3; ++j) CHECK(u.pi[j] == doctest::Approx(direct[j]).epsilon(1e-9));
}
