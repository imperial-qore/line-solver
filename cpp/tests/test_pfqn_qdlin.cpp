/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * pfqn_qdlin: the Linearizer arm of AMVA-LD on a plain demand matrix.
 *
 * Every golden below was measured with the native-Python pfqn_qdlin, which was
 * itself validated against SolverMVA(model,'qdlin') on 640 random closed models
 * (single server, multiserver, with and without think time, load dependent):
 * every metric agreed to 3e-16 relative, so these numbers are the solver's, not
 * an independent approximation of it. The iteration counts are pinned wherever
 * the fixed point converges, because a count agrees only if the convergence
 * test, the initial guess and every intermediate iterate agree too.
 *
 * WHAT THIS FUNCTION IS NOT: it is not the Wang-Sevcik QDLIN that used to carry
 * this name, which was Bard-Schweitzer written out; and it is not a textbook
 * Linearizer, because the reference's gamma is class-aggregate and lands in
 * slice 0 of a per-class array. That property is asserted here rather than
 * merely documented, so a future "fix" to the solver cannot pass silently.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_qdlin.h"

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

}  // namespace

TEST_CASE("one class matches the solver, iteration count included") {
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);

    const pfqn::QdLinResult<double> a = pfqn::pfqn_qdlin(L, N, Z);
    CHECK(a.Q(0, 0) == doctest::Approx(1.1055055966).epsilon(1e-9));
    CHECK(a.Q(1, 0) == doctest::Approx(0.5464847591).epsilon(1e-9));
    CHECK(a.U(0, 0) == doctest::Approx(0.6740047306).epsilon(1e-9));
    CHECK(a.U(1, 0) == doctest::Approx(0.4044028384).epsilon(1e-9));
    CHECK(a.R(0, 0) == doctest::Approx(0.8201022533).epsilon(1e-9));
    CHECK(a.X(0, 0) == doctest::Approx(1.3480094612).epsilon(1e-9));
    CHECK(a.C(0, 0) == doctest::Approx(2.2255034473).epsilon(1e-9));
    CHECK(a.iter == 397u);
}

TEST_CASE("two classes and three stations match the solver") {
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

    const pfqn::QdLinResult<double> a = pfqn::pfqn_qdlin(L, N, Z);
    CHECK(a.Q(0, 0) == doctest::Approx(1.2023471976).epsilon(1e-9));
    CHECK(a.Q(1, 0) == doctest::Approx(0.8455406444).epsilon(1e-9));
    CHECK(a.Q(2, 0) == doctest::Approx(0.2939223883).epsilon(1e-9));
    CHECK(a.Q(0, 1) == doctest::Approx(1.7138257835).epsilon(1e-9));
    CHECK(a.Q(1, 1) == doctest::Approx(0.1321216124).epsilon(1e-9));
    CHECK(a.Q(2, 1) == doctest::Approx(0.5172388252).epsilon(1e-9));
    CHECK(a.X(0, 0) == doctest::Approx(0.8290948878).epsilon(1e-9));
    CHECK(a.X(0, 1) == doctest::Approx(0.6368140023).epsilon(1e-9));
    CHECK(a.iter == 853u);

    // Little's law over the whole network, to the tolerance the fixed point
    // was actually converged to rather than to machine precision.
    for (std::size_t r = 0; r < 2; ++r) {
        double q = 0.0;
        for (std::size_t k = 0; k < 3; ++k) q += a.Q(k, r);
        CHECK(q + a.X(0, r) * Z[r] == doctest::Approx(N[r]).epsilon(1e-6));
    }
    // The cycle time is the residence times summed, plus the think time.
    for (std::size_t r = 0; r < 2; ++r) {
        double c = 0.0;
        for (std::size_t k = 0; k < 3; ++k) c += a.R(k, r);
        CHECK(a.C(0, r) == doctest::Approx(c + Z[r]).epsilon(1e-6));
    }
}

TEST_CASE("the reported utilization is the analytic T*S/c, not the iterate") {
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);
    std::vector<double> srv;
    srv.push_back(2.0);
    srv.push_back(1.0);

    const pfqn::QdLinResult<double> a =
        pfqn::pfqn_qdlin(L, N, Z, Matrix<double>(), srv);
    CHECK(a.Q(0, 0) == doctest::Approx(0.7844085389).epsilon(1e-9));
    CHECK(a.Q(1, 0) == doctest::Approx(0.6467257168).epsilon(1e-9));
    CHECK(a.X(0, 0) == doctest::Approx(1.568865852).epsilon(1e-9));
    CHECK(a.iter == 300u);
    // SolverMVA recomputes utilization from the NOMINAL demand and the server
    // count, with no lld or softmin scaling folded in.
    CHECK(a.U(0, 0) == doctest::Approx(a.X(0, 0) * L(0, 0) / 2.0).epsilon(1e-12));
    CHECK(a.U(1, 0) == doctest::Approx(a.X(0, 0) * L(1, 0)).epsilon(1e-12));
}

TEST_CASE("mu and nservers are different mechanisms, not two spellings of one") {
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);

    std::vector<double> srv;
    srv.push_back(2.0);
    srv.push_back(1.0);
    const pfqn::QdLinResult<double> viaServers =
        pfqn::pfqn_qdlin(L, N, Z, Matrix<double>(), srv);

    // The same two-server station spelled as a load-dependent rate lattice.
    Matrix<double> mu(2, 4, 1.0);
    for (std::size_t j = 0; j < 4; ++j) {
        mu(0, j) = std::min<double>(static_cast<double>(j + 1), 2.0);
        mu(1, j) = 1.0;
    }
    const pfqn::QdLinResult<double> viaMu =
        pfqn::pfqn_qdlin(L, N, Z, mu, std::vector<double>());
    CHECK(viaMu.Q(0, 0) == doctest::Approx(0.7844030387).epsilon(1e-9));
    CHECK(viaMu.iter == 300u);

    // A load-dependent model reports the ITERATED utilization, because the
    // analyzer forwards Uchain to the deaggregation only under lld/cd/jd
    // scaling. It is NOT X*L here: the lld term is folded into it.
    CHECK(viaMu.U(0, 0) == doctest::Approx(0.5088714356).epsilon(1e-9));
    CHECK(viaMu.U(1, 0) == doctest::Approx(0.4705485024).epsilon(1e-9));
    CHECK(viaMu.U(0, 0) != doctest::Approx(viaMu.X(0, 0) * L(0, 0)).epsilon(1e-6));
    // The server-count spelling has no lld scaling, so it reports the analytic
    // T*S/c instead -- the same station, two different reported utilizations.
    CHECK(viaServers.U(0, 0) == doctest::Approx(viaServers.X(0, 0) * L(0, 0) / 2.0).epsilon(1e-12));

    // Close, because both describe the same station, but NOT equal: the server
    // count goes through the softmin and the lattice through the interpolation.
    CHECK(viaMu.Q(0, 0) != doctest::Approx(viaServers.Q(0, 0)).epsilon(1e-12));
    CHECK(std::fabs(viaMu.Q(0, 0) - viaServers.Q(0, 0)) < 1e-4);
}

TEST_CASE("a single-server station still carries a softmin term") {
    // The multiserver factor is evaluated even at c = 1, where the softmin is
    // not exactly one, so qdlin does not reduce to a textbook single-server
    // AMVA. Spelling the same model with an explicit c = 1 must therefore give
    // the identical answer, and it does.
    Matrix<double> L;
    std::vector<double> N, Z;
    one_class(&L, &N, &Z);
    const pfqn::QdLinResult<double> implicit = pfqn::pfqn_qdlin(L, N, Z);
    const pfqn::QdLinResult<double> explicitOne =
        pfqn::pfqn_qdlin(L, N, Z, Matrix<double>(), std::vector<double>(2, 1.0));
    CHECK(implicit.Q(0, 0) == doctest::Approx(explicitOne.Q(0, 0)).epsilon(1e-15));
    CHECK(implicit.iter == explicitOne.iter);
    // A textbook single-server AMVA would put the utilization at X*L exactly,
    // and so does the REPORTED one; the softmin lives in the residence time.
    CHECK(implicit.U(0, 0) == doctest::Approx(implicit.X(0, 0) * L(0, 0)).epsilon(1e-12));
}

TEST_CASE("think time may be absent and a class may be empty") {
    Matrix<double> L(2, 2, 0.0);
    L(0, 0) = 0.5;
    L(1, 0) = 0.3;
    L(0, 1) = 0.2;
    L(1, 1) = 0.7;
    std::vector<double> N;
    N.push_back(3.0);
    N.push_back(2.0);

    // No delay station is appended when there is no think time.
    const pfqn::QdLinResult<double> nz = pfqn::pfqn_qdlin(L, N, std::vector<double>());
    CHECK(nz.Q(0, 0) == doctest::Approx(1.6035096548).epsilon(1e-9));
    CHECK(nz.Q(1, 1) == doctest::Approx(1.4989846125).epsilon(1e-9));
    CHECK(nz.X(0, 0) == doctest::Approx(1.3103309302).epsilon(1e-9));
    for (std::size_t r = 0; r < 2; ++r) {
        double q = 0.0;
        for (std::size_t k = 0; k < 2; ++k) q += nz.Q(k, r);
        CHECK(q == doctest::Approx(N[r]).epsilon(1e-6));
    }

    // An empty population is answered with zeros, not refused, and no sweep runs.
    const pfqn::QdLinResult<double> e =
        pfqn::pfqn_qdlin(L, std::vector<double>(2, 0.0), std::vector<double>());
    CHECK(e.iter == 0u);
    for (std::size_t k = 0; k < 2; ++k)
        for (std::size_t r = 0; r < 2; ++r) CHECK(e.Q(k, r) == doctest::Approx(0.0));

    // One class present, one empty: the empty one stays at zero throughout.
    std::vector<double> N1;
    N1.push_back(2.0);
    N1.push_back(0.0);
    const pfqn::QdLinResult<double> p = pfqn::pfqn_qdlin(L, N1, std::vector<double>());
    CHECK(p.X(0, 1) == doctest::Approx(0.0));
    for (std::size_t k = 0; k < 2; ++k) CHECK(p.Q(k, 1) == doctest::Approx(0.0));
    CHECK(p.X(0, 0) > 0.0);
}

TEST_CASE("the refusals are by name") {
    Matrix<double> L(2, 2, 0.5);
    std::vector<double> N(2, 1.0);
    CHECK_THROWS_AS(pfqn::pfqn_qdlin(L, std::vector<double>(3, 1.0), N), line::InputError);
    CHECK_THROWS_AS(pfqn::pfqn_qdlin(L, N, std::vector<double>(5, 1.0)), line::InputError);
    CHECK_THROWS_AS(
        pfqn::pfqn_qdlin(L, N, N, Matrix<double>(), std::vector<double>(4, 1.0)),
        line::InputError);
    std::vector<double> inf(2, 1.0);
    inf[0] = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(pfqn::pfqn_qdlin(L, inf, N), line::InputError);
}
