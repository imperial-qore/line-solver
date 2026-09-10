/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_qli / pfqn_fli: the two Wang-Sevcik approximate MVA schemes.
 *
 * Every golden below was measured in native Python on the same inputs, and the
 * port reproduces them to ten digits INCLUDING THE ITERATION COUNTS -- which is
 * the sharper check of the two, because an iteration count agrees only if the
 * convergence test, the initial guess and every intermediate iterate agree too.
 *
 * The schemes are APPROXIMATIONS, so the second oracle is exact single-class
 * MVA: with one class the exact arrival queue is known in closed form, and both
 * corrections overestimate it in a fixed order (QLI then FLI). That ordering is
 * a property of the corrections rather than of these numbers, so it is asserted
 * separately.
 *
 * pfqn_qdlin was tested here until it moved to line/api/pfqn/pfqn_qdlin.h; its
 * Wang-Sevcik arm was Bard-Schweitzer under another name. See
 * test_pfqn_qdlin.cpp.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_wangsevcik.h"

namespace pfqn = line::pfqn;
using line::Matrix;

namespace {

/** Two stations, one class, demands 0.5 and 0.3, three jobs, think time one. */
void one_class(Matrix<double>* L, std::vector<double>* N, std::vector<double>* Z) {
    *L = Matrix<double>(2, 1, 0.0);
    (*L)(0, 0) = 0.5;
    (*L)(1, 0) = 0.3;
    N->assign(1, 3.0);
    Z->assign(1, 1.0);
}

/** Exact single-class MVA, the closed recursion. */
void exact_mva(const Matrix<double>& L, double N, double Z, double* Q) {
    const std::size_t M = L.rows();
    for (std::size_t k = 0; k < M; ++k) Q[k] = 0.0;
    for (int n = 1; n <= static_cast<int>(N); ++n) {
        std::vector<double> R(M, 0.0);
        double tot = 0.0;
        for (std::size_t k = 0; k < M; ++k) {
            R[k] = L(k, 0) * (1.0 + Q[k]);
            tot += R[k];
        }
        const double X = static_cast<double>(n) / (Z + tot);
        for (std::size_t k = 0; k < M; ++k) Q[k] = X * R[k];
    }
}

}  // namespace

TEST_CASE("both schemes match native Python, iteration counts included") {
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);

    const pfqn::WsResult<double> b = pfqn::pfqn_qli(L, N, Z);
    CHECK(b.Q(0, 0) == doctest::Approx(1.197328922).epsilon(1e-8));
    CHECK(b.Q(1, 0) == doctest::Approx(0.5496876719).epsilon(1e-8));
    CHECK(b.X[0] == doctest::Approx(1.252983406).epsilon(1e-8));
    CHECK(b.iterations == 9u);

    const pfqn::WsResult<double> c = pfqn::pfqn_fli(L, N, Z);
    CHECK(c.Q(0, 0) == doctest::Approx(1.227250062).epsilon(1e-8));
    CHECK(c.Q(1, 0) == doctest::Approx(0.5660640186).epsilon(1e-8));
    CHECK(c.X[0] == doctest::Approx(1.20668592).epsilon(1e-8));
    CHECK(c.iterations == 7u);
}

TEST_CASE("two classes match native Python too") {
    Matrix<double> L(2, 2, 0.0);
    L(0, 0) = 0.5;
    L(1, 0) = 0.3;
    L(0, 1) = 0.2;
    L(1, 1) = 0.7;
    std::vector<double> N, Z;
    N.push_back(3.0);
    N.push_back(2.0);
    Z.push_back(1.0);
    Z.push_back(0.5);

    const pfqn::WsResult<double> b = pfqn::pfqn_qli(L, N, Z);
    CHECK(b.Q(0, 0) == doctest::Approx(1.06227).epsilon(1e-5));
    CHECK(b.X[1] == doctest::Approx(0.70355).epsilon(1e-5));

    const pfqn::WsResult<double> c = pfqn::pfqn_fli(L, N, Z);
    CHECK(c.Q(0, 0) == doctest::Approx(1.07356).epsilon(1e-5));
    CHECK(c.X[1] == doctest::Approx(0.701047).epsilon(1e-5));
}

TEST_CASE("both bracket exact single-class MVA in a fixed order") {
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);
    double Qe[2];
    exact_mva(L, N[0], Z[0], Qe);
    CHECK(Qe[0] == doctest::Approx(1.106372303).epsilon(1e-8));

    const double qb = pfqn::pfqn_qli(L, N, Z).Q(0, 0);
    const double qc = pfqn::pfqn_fli(L, N, Z).Q(0, 0);
    // Both overestimate the busiest station's queue, and they do so in this
    // order -- a property of the corrections, not of these numbers.
    CHECK(qb > Qe[0]);
    CHECK(qc > qb);
    // None is wildly off: within 12 per cent of exact on this instance.
    CHECK((qc - Qe[0]) / Qe[0] < 0.12);
}

TEST_CASE("the conservation laws hold for every scheme") {
    Matrix<double> L(3, 2, 0.0);
    L(0, 0) = 0.4;
    L(1, 0) = 0.6;
    L(2, 0) = 0.2;
    L(0, 1) = 0.9;
    L(1, 1) = 0.1;
    L(2, 1) = 0.5;
    std::vector<double> N, Z;
    N.push_back(4.0);
    N.push_back(3.0);
    Z.push_back(2.0);
    Z.push_back(1.0);

    const pfqn::WsScheme all[2] = {pfqn::WsScheme::Qli, pfqn::WsScheme::Fli};
    for (std::size_t s = 0; s < 2; ++s) {
        const pfqn::WsResult<double> r = pfqn::pfqn_wangsevcik(L, N, Z, all[s]);
        for (std::size_t cls = 0; cls < 2; ++cls) {
            // Little's law over the whole network: the jobs in the queues plus
            // those thinking add up to the population.
            double q = 0.0;
            for (std::size_t k = 0; k < 3; ++k) q += r.Q(k, cls);
            CHECK(q + r.X[cls] * Z[cls] == doctest::Approx(N[cls]).epsilon(1e-6));
            // Utilization is throughput times demand.
            for (std::size_t k = 0; k < 3; ++k) {
                CHECK(r.U(k, cls) == doctest::Approx(r.X[cls] * L(k, cls)).epsilon(1e-9));
                CHECK(r.Q(k, cls) >= 0.0);
            }
            // The cycle time is the residence times summed.
            double c = 0.0;
            for (std::size_t k = 0; k < 3; ++k) c += r.R(k, cls);
            CHECK(r.C[cls] == doctest::Approx(c).epsilon(1e-9));
        }
    }
}

TEST_CASE("an empty network and an absent class are answered, not refused") {
    Matrix<double> L(2, 2, 0.5);
    std::vector<double> Zero(2, 0.0), Z(2, 0.0);
    // No jobs at all: every metric is zero and no iteration runs.
    const pfqn::WsResult<double> e = pfqn::pfqn_qli(L, Zero, Z);
    CHECK(e.iterations == 0u);
    for (std::size_t k = 0; k < 2; ++k)
        for (std::size_t r = 0; r < 2; ++r) CHECK(e.Q(k, r) == doctest::Approx(0.0));

    // One class present, one empty: the empty one stays at zero throughout.
    std::vector<double> N;
    N.push_back(2.0);
    N.push_back(0.0);
    const pfqn::WsResult<double> r = pfqn::pfqn_qli(L, N, Z);
    CHECK(r.X[1] == doctest::Approx(0.0));
    for (std::size_t k = 0; k < 2; ++k) CHECK(r.Q(k, 1) == doctest::Approx(0.0));
    CHECK(r.X[0] > 0.0);

    // No think time at all is a valid model, and an empty Z means exactly that.
    const pfqn::WsResult<double> nz = pfqn::pfqn_qli(L, N, std::vector<double>());
    CHECK(nz.X[0] > 0.0);
}

TEST_CASE("the refusals are by name") {
    Matrix<double> L(2, 2, 0.5);
    std::vector<double> N(2, 1.0);
    CHECK_THROWS_AS(pfqn::pfqn_qli(Matrix<double>(0, 0, 0.0), N, N), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_qli(L, std::vector<double>(3, 1.0), N), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_qli(L, N, std::vector<double>(5, 1.0)), line::InputError);
}
