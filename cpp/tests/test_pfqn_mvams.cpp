/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Tests for the multiserver MVA dispatcher and the three load-dependent and
 * mixed variants it dispatches to. The oracles are:
 *
 *  1. pfqn_mva. With one server everywhere the multiserver rates degenerate to
 *     mu(k) = 1, so pfqn_mvald must reproduce pfqn_mva as a rational, with no
 *     rounding, and pfqn_mvams must take the pfqn_mva branch outright.
 *  2. pfqn_gld. The normalizing constant that pfqn_mvald accumulates along the
 *     lattice is the constant of the same load-dependent model, so it must
 *     equal pfqn_gld's to the last bit of the rational.
 *  3. The laws. Little's law station by station, the utilization law
 *     U = X L / S, and the population constraint sum_i Q(i,r) + X(r) Z_r = N_r,
 *     all checked as exact identities rather than to a tolerance.
 *  4. MATLAB. Reference values produced by pfqn_mvams and pfqn_mvald on the
 *     same models are embedded for the four branches.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mvams.h"

using line::InputError;
using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::INF_SERVERS;
using line::pfqn::MvaResult;
using line::pfqn::OPEN_CLASS;
using line::pfqn::pfqn_gld;
using line::pfqn::pfqn_mva;
using line::pfqn::pfqn_mvald;
using line::pfqn::pfqn_mvaldms;
using line::pfqn::pfqn_mvamx;
using line::pfqn::pfqn_mvams;

namespace {

constexpr double MATLAB_TOL = 1e-12;

/** Two stations, two classes; the demand matrix shared by most cases. */
template <class T>
Matrix<T> demands2x2() {
    Matrix<T> D(2, 2);
    D(0, 0) = line::num_traits<T>::from_rational(1, 5);   // 0.2
    D(0, 1) = line::num_traits<T>::from_rational(1, 2);   // 0.5
    D(1, 0) = line::num_traits<T>::from_rational(3, 10);  // 0.3
    D(1, 1) = line::num_traits<T>::from_rational(2, 5);   // 0.4
    return D;
}

/** Think times (0, 1) as a one-row matrix. */
template <class T>
Matrix<T> think01() {
    Matrix<T> Z(1, 2);
    Z(0, 0) = line::num_traits<T>::from_int(0);
    Z(0, 1) = line::num_traits<T>::from_int(1);
    return Z;
}

/** Exact-arithmetic closed model with rational demands. */
template <class T>
Matrix<T> demands_exact() {
    Matrix<T> L(2, 2);
    L(0, 0) = line::num_traits<T>::from_rational(1, 2);
    L(0, 1) = line::num_traits<T>::from_rational(3, 10);
    L(1, 0) = line::num_traits<T>::from_rational(2, 5);
    L(1, 1) = line::num_traits<T>::from_rational(3, 5);
    return L;
}

/** (M x Nt) rate matrix with mu(i,k) = min(k, S(i)). */
template <class T>
Matrix<T> ms_rates(const std::vector<int>& S, std::size_t Nt) {
    Matrix<T> mu(S.size(), Nt);
    for (std::size_t i = 0; i < S.size(); ++i)
        for (std::size_t k = 1; k <= Nt; ++k) {
            const long c = static_cast<long>(k) < S[i] ? static_cast<long>(k) : S[i];
            mu(i, k - 1) = line::num_traits<T>::from_int(c);
        }
    return mu;
}

}  // namespace

TEST_CASE("pfqn_mvald with unit rates is pfqn_mva, exactly in rational arithmetic") {
    const Matrix<Rational> L = demands_exact<Rational>();
    const std::vector<int> N{2, 1};
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(3, 10);
    Z(0, 1) = Rational(1, 5);

    const Matrix<Rational> mu(2, 3, Rational(1));
    const auto ld = pfqn_mvald(L, N, Z, mu);
    const auto mv = pfqn_mva(L, N, Z);

    REQUIRE(ld.isNumStable);
    for (std::size_t r = 0; r < 2; ++r) CHECK(ld.XN[r] == mv.XN[r]);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(ld.QN(i, r) == mv.QN(i, r));
            // pfqn_mvald's WN is pfqn_mva's per-station residence time.
            CHECK(ld.WN(i, r) == mv.CN(i, r));
        }
    CHECK(ld.G == mv.G);

    // pfqn_mvald reports the aggregate utilization 1 - P(empty), which for a
    // single server is the sum of pfqn_mva's per-class utilizations.
    for (std::size_t i = 0; i < 2; ++i) {
        Rational rowsum(0);
        for (std::size_t r = 0; r < 2; ++r) rowsum += mv.UN(i, r);
        CHECK(ld.UN[i] == rowsum);
    }
}

TEST_CASE("pfqn_mvams with unit server counts is pfqn_mva, exactly") {
    const Matrix<Rational> L = demands_exact<Rational>();
    const std::vector<int> N{2, 1};
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(3, 10);
    Z(0, 1) = Rational(1, 5);

    const auto ms = pfqn_mvams(std::vector<Rational>{}, L, N, Z, std::vector<int>{1, 1},
                               std::vector<int>{1, 1});
    const auto mv = pfqn_mva(L, N, Z);

    for (std::size_t r = 0; r < 2; ++r) CHECK(ms.XN[r] == mv.XN[r]);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(ms.QN(i, r) == mv.QN(i, r));
            CHECK(ms.UN(i, r) == mv.UN(i, r));
            CHECK(ms.CN(i, r) == mv.CN(i, r));
        }
    CHECK(ms.G == mv.G);

    // The overloads that default the multiplicities and the server counts must
    // reach the same branch.
    CHECK(pfqn_mvams(std::vector<Rational>{}, L, N, Z).G == mv.G);
    CHECK(pfqn_mvams(std::vector<Rational>{}, L, N, Z, std::vector<int>{1, 1}).G == mv.G);
}

TEST_CASE("pfqn_mvams multiserver constant equals pfqn_gld on the same rates") {
    const Matrix<Rational> L = demands_exact<Rational>();
    const std::vector<int> N{2, 1};
    const std::vector<int> S{2, 1};
    const Matrix<Rational> mu = ms_rates<Rational>(S, 3);

    // No think time: the two algorithms see exactly the same model.
    const auto ms = pfqn_mvams(std::vector<Rational>{}, L, N, Matrix<Rational>(), std::vector<int>{},
                               S);
    CHECK(ms.G == pfqn_gld(L, N, mu).G);

    // With think time, pfqn_gld carries the delay as an infinite-server row.
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(3, 10);
    Z(0, 1) = Rational(1, 5);
    Matrix<Rational> Lz(3, 2);
    Matrix<Rational> muz(3, 3);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) Lz(i, r) = L(i, r);
    for (std::size_t r = 0; r < 2; ++r) Lz(2, r) = Z(0, r);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 3; ++k) muz(i, k) = mu(i, k);
    for (std::size_t k = 0; k < 3; ++k) muz(2, k) = Rational(static_cast<long>(k) + 1);

    const auto msz = pfqn_mvams(std::vector<Rational>{}, L, N, Z, std::vector<int>{}, S);
    CHECK(msz.G == pfqn_gld(Lz, N, muz).G);
}

TEST_CASE("pfqn_mvams closed metrics satisfy Little, the utilization law and the population constraint") {
    const Matrix<Rational> L = demands_exact<Rational>();
    const std::vector<int> N{2, 1};
    const std::vector<int> S{2, 1};
    Matrix<Rational> Z(1, 2);
    Z(0, 0) = Rational(3, 10);
    Z(0, 1) = Rational(1, 5);

    const auto r = pfqn_mvams(std::vector<Rational>{}, L, N, Z, std::vector<int>{}, S);

    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 2; ++c) {
            CHECK(r.QN(i, c) == r.XN[c] * r.CN(i, c));                        // Little's law
            CHECK(r.UN(i, c) == r.XN[c] * L(i, c) / Rational(S[i]));          // utilization law
        }

    for (std::size_t c = 0; c < 2; ++c) {
        Rational pop = r.XN[c] * Z(0, c);
        for (std::size_t i = 0; i < 2; ++i) pop += r.QN(i, c);
        CHECK(pop == Rational(N[c]));  // every job is somewhere
    }

    // No station may be more than fully utilized.
    for (std::size_t i = 0; i < 2; ++i) {
        Rational rho(0);
        for (std::size_t c = 0; c < 2; ++c) rho += r.UN(i, c);
        CHECK(rho <= Rational(1));
    }
}

TEST_CASE("pfqn_mvams closed multiserver matches MATLAB") {
    // MATLAB: [X,Q,U,C,lG] = pfqn_mvams([0 0], [0.2 0.5; 0.3 0.4], [2 3],
    //                                   [0 1], [1;1], [2;1])
    const Matrix<double> D = demands2x2<double>();
    const auto r =
        pfqn_mvams(std::vector<double>{0.0, 0.0}, D, std::vector<int>{2, 3}, think01<double>(),
                   std::vector<int>{1, 1}, std::vector<int>{2, 1});

    CHECK(r.XN[0] == doctest::Approx(1.83754541253049).epsilon(MATLAB_TOL));
    CHECK(r.XN[1] == doctest::Approx(1.02555665404813).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 0) == doctest::Approx(0.427590477315796).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 1) == doctest::Approx(0.590605944124452).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 0) == doctest::Approx(1.5724095226842).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 1) == doctest::Approx(1.38383740182742).epsilon(MATLAB_TOL));
    CHECK(r.CN(0, 0) == doctest::Approx(0.232696549647151).epsilon(MATLAB_TOL));
    CHECK(r.CN(0, 1) == doctest::Approx(0.57588816940847).epsilon(MATLAB_TOL));
    CHECK(r.CN(1, 0) == doctest::Approx(0.855711925246425).epsilon(MATLAB_TOL));
    CHECK(r.CN(1, 1) == doctest::Approx(1.34935246762338).epsilon(MATLAB_TOL));
    CHECK(r.lG == doctest::Approx(-0.58490133148295).epsilon(MATLAB_TOL));

    // MATLAB's UN on this branch is pfqn_mvald's (M x 1) aggregate 1 - P(empty)
    // rather than the documented (M x R) matrix; see the contract note in
    // pfqn_mvams.h. The single-server station, where the two definitions must
    // coincide, is checked against the MATLAB value.
    const auto ld = pfqn_mvald(D, std::vector<int>{2, 3}, think01<double>(),
                               ms_rates<double>(std::vector<int>{2, 1}, 5));
    CHECK(ld.UN[0] == doctest::Approx(0.617704642435682).epsilon(MATLAB_TOL));
    CHECK(ld.UN[1] == doctest::Approx(0.961486285378399).epsilon(MATLAB_TOL));
    CHECK(r.UN(1, 0) + r.UN(1, 1) == doctest::Approx(ld.UN[1]).epsilon(MATLAB_TOL));
}

TEST_CASE("pfqn_mvald matches MATLAB on a three-station load-dependent model") {
    // MATLAB: L = [0.5 0.3; 0.4 0.6; 0.2 0.9], N = [2 2], Z = [0.3 0.2],
    // mu rows min(1:4,2), min(1:4,3), 1:4 (a two-server, a three-server and an
    // infinite-server station).
    Matrix<double> L(3, 2);
    L(0, 0) = 0.5;
    L(0, 1) = 0.3;
    L(1, 0) = 0.4;
    L(1, 1) = 0.6;
    L(2, 0) = 0.2;
    L(2, 1) = 0.9;
    Matrix<double> Z(1, 2);
    Z(0, 0) = 0.3;
    Z(0, 1) = 0.2;
    Matrix<double> mu = ms_rates<double>(std::vector<int>{2, 3, INF_SERVERS}, 4);
    for (std::size_t k = 1; k <= 4; ++k) mu(2, k - 1) = static_cast<double>(k);

    const auto r = pfqn_mvald(L, std::vector<int>{2, 2}, Z, mu);
    CHECK(r.XN[0] == doctest::Approx(1.39305269938802).epsilon(MATLAB_TOL));
    CHECK(r.XN[1] == doctest::Approx(0.980525755060578).epsilon(MATLAB_TOL));
    const double Q[6] = {0.741495503084661, 0.328349148930899, 0.561978147221327,
                         0.593072520502465, 0.278610539877605, 0.88247317955452};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(r.QN(i, c) == doctest::Approx(Q[i * 2 + c]).epsilon(MATLAB_TOL));
    CHECK(r.UN[0] == doctest::Approx(0.707626173781621).epsilon(MATLAB_TOL));
    CHECK(r.UN[1] == doctest::Approx(0.738609053294021).epsilon(MATLAB_TOL));
    CHECK(r.UN[2] == doctest::Approx(0.764859641733356).epsilon(MATLAB_TOL));
    CHECK(r.CN[0] == doctest::Approx(1.13569586482881).epsilon(MATLAB_TOL));
    CHECK(r.CN[1] == doctest::Approx(1.83972204674668).epsilon(MATLAB_TOL));
    CHECK(r.lG == doctest::Approx(0.70213169863141).epsilon(MATLAB_TOL));

    // The marginal distribution the recursion carries must be a distribution.
    for (std::size_t i = 0; i < 3; ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < r.PI.cols(); ++k) s += r.PI(i, k);
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_mvams mixed single-server matches MATLAB and pfqn_mvamx") {
    // MATLAB: pfqn_mvams([0.5 0], [0.2 0.5; 0.3 0.4], [Inf 3], [0 1], [1;1], [1;1])
    const Matrix<double> D = demands2x2<double>();
    const std::vector<double> lambda{0.5, 0.0};
    const std::vector<int> N{OPEN_CLASS, 3};
    const auto r =
        pfqn_mvams(lambda, D, N, think01<double>(), std::vector<int>{1, 1}, std::vector<int>{1, 1});

    CHECK(r.XN[0] == doctest::Approx(0.5).epsilon(MATLAB_TOL));
    CHECK(r.XN[1] == doctest::Approx(1.15008307371642).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 0) == doctest::Approx(0.226037256877681).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 1) == doctest::Approx(1.03433531189913).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 0) == doctest::Approx(0.320396755479608).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 1) == doctest::Approx(0.815581614384444).epsilon(MATLAB_TOL));
    CHECK(r.UN(0, 0) == doctest::Approx(0.1).epsilon(MATLAB_TOL));
    CHECK(r.UN(0, 1) == doctest::Approx(0.575041536858212).epsilon(MATLAB_TOL));
    CHECK(r.UN(1, 0) == doctest::Approx(0.15).epsilon(MATLAB_TOL));
    CHECK(r.UN(1, 1) == doctest::Approx(0.460033229486569).epsilon(MATLAB_TOL));
    CHECK(r.CN(0, 0) == doctest::Approx(0.452074513755363).epsilon(MATLAB_TOL));
    CHECK(r.CN(1, 1) == doctest::Approx(0.709150176212003).epsilon(MATLAB_TOL));
    CHECK(r.lG == doctest::Approx(0.700731196706087).epsilon(MATLAB_TOL));

    // The dispatcher must have taken the pfqn_mvamx branch verbatim.
    const auto mx = pfqn_mvamx(lambda, D, N, think01<double>(), std::vector<int>{1, 1});
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 2; ++c) CHECK(r.QN(i, c) == mx.QN(i, c));

    // The open utilization is the offered load, and Little's law holds.
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(r.UN(i, 0) == doctest::Approx(lambda[0] * D(i, 0)).epsilon(MATLAB_TOL));
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(r.QN(i, c) == doctest::Approx(r.XN[c] * r.CN(i, c)).epsilon(MATLAB_TOL));
    }
}

TEST_CASE("pfqn_mvams mixed multiserver matches MATLAB and has no normalizing constant") {
    // MATLAB: pfqn_mvams([0.5 0], [0.2 0.5; 0.3 0.4], [Inf 3], [0 1], [1;1], [2;1])
    const Matrix<double> D = demands2x2<double>();
    const std::vector<double> lambda{0.5, 0.0};
    const std::vector<int> N{OPEN_CLASS, 3};
    const auto r =
        pfqn_mvams(lambda, D, N, think01<double>(), std::vector<int>{1, 1}, std::vector<int>{2, 1});

    CHECK(r.XN[0] == doctest::Approx(0.5).epsilon(MATLAB_TOL));
    CHECK(r.XN[1] == doctest::Approx(1.32644908070282).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 0) == doctest::Approx(0.111838089262978).epsilon(MATLAB_TOL));
    CHECK(r.QN(0, 1) == doctest::Approx(0.701244092052267).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 0) == doctest::Approx(0.348054145984397).epsilon(MATLAB_TOL));
    CHECK(r.QN(1, 1) == doctest::Approx(0.972306827244915).epsilon(MATLAB_TOL));
    CHECK(r.UN(0, 0) == doctest::Approx(0.05).epsilon(MATLAB_TOL));
    CHECK(r.UN(0, 1) == doctest::Approx(0.331612270175705).epsilon(MATLAB_TOL));
    CHECK(r.UN(1, 0) == doctest::Approx(0.15).epsilon(MATLAB_TOL));
    CHECK(r.UN(1, 1) == doctest::Approx(0.530579632281128).epsilon(MATLAB_TOL));
    CHECK(r.CN(0, 0) == doctest::Approx(0.223676178525955).epsilon(MATLAB_TOL));
    CHECK(r.CN(1, 1) == doctest::Approx(0.733014814808977).epsilon(MATLAB_TOL));

    // MATLAB reports lG = NaN on this branch.
    CHECK(std::isnan(r.lG));

    // Utilization law with the server count in the denominator.
    for (std::size_t i = 0; i < 2; ++i) {
        const double S = i == 0 ? 2.0 : 1.0;
        CHECK(r.UN(i, 0) == doctest::Approx(lambda[0] * D(i, 0) / S).epsilon(MATLAB_TOL));
        CHECK(r.UN(i, 1) == doctest::Approx(r.XN[1] * D(i, 1) / S).epsilon(MATLAB_TOL));
    }
}

TEST_CASE("pfqn_mvaldms with one server everywhere reproduces pfqn_mvamx") {
    // Two exact algorithms for the same product-form model, reached by
    // different recursions: the limited-load-dependent one must collapse onto
    // the single-server one when every rate is 1.
    const Matrix<double> D = demands2x2<double>();
    const std::vector<double> lambda{0.5, 0.0};
    const std::vector<int> N{OPEN_CLASS, 3};
    const auto ldms = pfqn_mvaldms(lambda, D, N, think01<double>(), std::vector<int>{1, 1});
    const auto mx = pfqn_mvamx(lambda, D, N, think01<double>(), std::vector<int>{1, 1});

    for (std::size_t c = 0; c < 2; ++c) CHECK(ldms.XN[c] == doctest::Approx(mx.XN[c]).epsilon(1e-12));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 2; ++c) {
            CHECK(ldms.QN(i, c) == doctest::Approx(mx.QN(i, c)).epsilon(1e-12));
            CHECK(ldms.UN(i, c) == doctest::Approx(mx.UN(i, c)).epsilon(1e-12));
            CHECK(ldms.CN(i, c) == doctest::Approx(mx.CN(i, c)).epsilon(1e-12));
        }
    CHECK(ldms.lG == doctest::Approx(mx.lG).epsilon(1e-12));
}

TEST_CASE("pfqn_mvams runs in high precision as well as in double") {
    const Matrix<Real50> L = demands_exact<Real50>();
    const std::vector<int> N{2, 1};
    const auto r =
        pfqn_mvams(std::vector<Real50>{}, L, N, Matrix<Real50>(), std::vector<int>{}, std::vector<int>{2, 1});
    const Matrix<double> Ld = demands_exact<double>();
    const auto rd =
        pfqn_mvams(std::vector<double>{}, Ld, N, Matrix<double>(), std::vector<int>{}, std::vector<int>{2, 1});
    for (std::size_t c = 0; c < 2; ++c)
        CHECK(static_cast<double>(r.XN[c]) == doctest::Approx(rd.XN[c]).epsilon(1e-12));
    CHECK(r.lG == doctest::Approx(rd.lG).epsilon(1e-12));
}

TEST_CASE("pfqn_mvams edge cases follow the MATLAB contract") {
    const Matrix<Rational> L = demands_exact<Rational>();

    // Empty population: no throughput, and the absent-class residence time is
    // L times the multiplicity, as in pfqn_mva.
    const auto none = pfqn_mvams(std::vector<Rational>{}, L, std::vector<int>{0, 0},
                                 Matrix<Rational>(), std::vector<int>{}, std::vector<int>{1, 1});
    for (std::size_t c = 0; c < 2; ++c) CHECK(none.XN[c] == Rational(0));
    CHECK(none.G == Rational(1));

    // One class absent while the other is populated, on the multiserver branch.
    const auto part = pfqn_mvams(std::vector<Rational>{}, L, std::vector<int>{2, 0},
                                 Matrix<Rational>(), std::vector<int>{}, std::vector<int>{2, 1});
    CHECK(part.XN[1] == Rational(0));
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(part.CN(i, 1) == L(i, 1));      // L times mi, with mi = 1
        CHECK(part.QN(i, 1) == Rational(0));  // no jobs of an absent class
    }

    // An infinite-server station alone does not make the model multiserver, so
    // the dispatcher stays on the pfqn_mva branch, exactly as MATLAB's isfinite
    // guard does.
    const auto inf = pfqn_mvams(std::vector<Rational>{}, L, std::vector<int>{2, 1},
                                Matrix<Rational>(), std::vector<int>{}, std::vector<int>{INF_SERVERS, 1});
    CHECK(inf.G == pfqn_mva(L, std::vector<int>{2, 1}).G);
}

TEST_CASE("pfqn_mvams rejects inputs the exact algorithms cannot honour") {
    const Matrix<double> D = demands2x2<double>();

    // Queue replicas together with multiservers: MATLAB rejects this for mixed
    // models and silently drops mi for closed ones; the port rejects both.
    CHECK_THROWS_AS(pfqn_mvams(std::vector<double>{0.5, 0.0}, D, std::vector<int>{OPEN_CLASS, 3},
                               think01<double>(), std::vector<int>{2, 1}, std::vector<int>{2, 1}),
                    InputError);
    CHECK_THROWS_AS(pfqn_mvams(std::vector<double>{0.0, 0.0}, D, std::vector<int>{2, 3},
                               think01<double>(), std::vector<int>{2, 1}, std::vector<int>{2, 1}),
                    InputError);

    // An arrival rate on a closed class.
    CHECK_THROWS_AS(pfqn_mvamx(std::vector<double>{0.5, 0.5}, D, std::vector<int>{OPEN_CLASS, 3},
                               think01<double>(), std::vector<int>{}),
                    InputError);

    // Dimension mismatches.
    CHECK_THROWS_AS(pfqn_mvams(std::vector<double>{0.0}, D, std::vector<int>{2},
                               Matrix<double>(), std::vector<int>{}, std::vector<int>{}),
                    InputError);
    CHECK_THROWS_AS(pfqn_mvams(std::vector<double>{0.0, 0.0}, D, std::vector<int>{2, 1},
                               Matrix<double>(), std::vector<int>{1, 1, 1}, std::vector<int>{}),
                    InputError);

    // Rates given for fewer jobs than the population.
    CHECK_THROWS_AS(pfqn_mvald(D, std::vector<int>{2, 3}, think01<double>(), Matrix<double>(2, 2, 1.0)),
                    InputError);
}
