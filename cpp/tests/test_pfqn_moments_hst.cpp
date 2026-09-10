/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Queue-length joint moments, PS sojourn moments, the HST robustness
 * certificate, and the FES / access-graph helpers. Oracles, in order of
 * strength:
 *   1. Cross-algorithm identity. The mean of the joint moment array IS the mean
 *      queue length, so pfqn_qlen_joint_moments must reproduce pfqn_mva term for
 *      term on any product-form model; and the single-class exact PS sojourn of
 *      Mitra-Morrison must equal the MVA RESIDENCE TIME at that station, the two
 *      being the same quantity computed from different objects.
 *   2. Lemma 3.1 of Suri (1983): the unconstrained HST certificate sum_n |c_n|
 *      equals Q_i(N) - Q_i(N-1). Checked against pfqn_mva at two populations.
 *   3. Published values of the reference implementation, quoted to 12 digits.
 *   4. Definitional invariants: a survival array is non-increasing and starts at
 *      one; a covariance matrix is symmetric; the (P1) optimum is feasible and
 *      no larger than the unconstrained certificate.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_gamma.h"
#include "line/api/cache/cache_gamma_lp.h"
#include "line/api/pfqn/pfqn_hst.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_pff_delay.h"
#include "line/api/pfqn/pfqn_qlen_joint_moments.h"
#include "line/api/pfqn/pfqn_respt_ps_moments.h"
#include "line/api/pfqn/pfqn_xia.h"

using line::Matrix;
using line::Rational;
using line::pfqn::pfqn_hst;
using line::pfqn::pfqn_pff_delay;
using line::pfqn::pfqn_qlen_joint_moments;
using line::pfqn::pfqn_respt_ps_moments;
using line::pfqn::pfqn_xia;
using line::pfqn::QlenJointRoute;
using line::pfqn::ResptPsMethod;
using line::pfqn::ResptPsRoute;

namespace {

constexpr double TOL = 1e-9;

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_hst
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_hst reproduces the reference certificate term for term") {
    // matlab: pfqn_hst([0.9;0.4;0.2],5,1.0), every value quoted to 12 digits
    const std::vector<double> L = {0.9, 0.4, 0.2};
    const auto s = pfqn_hst<double>(L, 5, 1.0);
    CHECK(s.station == 0);  // the bottleneck is picked by default
    CHECK(s.X == doctest::Approx(1.05231780339).epsilon(TOL));
    CHECK(s.U == doctest::Approx(0.947086023051).epsilon(TOL));
    CHECK(s.Q == doctest::Approx(3.02199199328).epsilon(TOL));
    CHECK(s.total == doctest::Approx(0.83116041122).epsilon(TOL));
    CHECK(s.worst == doctest::Approx(0.088391089577).epsilon(TOL));

    const double Pgeq[6] = {1, 0.947086023051, 0.84463330281, 0.667547052725, 0.413768834334,
                            0.14895678036};
    const double p[6] = {0.0529139769485, 0.102452720242, 0.177086250085,
                         0.253778218391,  0.264812053973, 0.14895678036};
    const double c[5] = {-0.0552628071534, -0.139790198301, -0.230660830922, -0.256489794484,
                         -0.14895678036};
    const double astar[5] = {1, 1, 0.528925866285, -1, -1};
    REQUIRE(s.Pgeq.size() == 6);
    REQUIRE(s.c.size() == 5);
    for (int k = 0; k < 6; ++k) {
        CHECK(s.Pgeq[k] == doctest::Approx(Pgeq[k]).epsilon(TOL));
        CHECK(s.p[k] == doctest::Approx(p[k]).epsilon(TOL));
    }
    for (int k = 0; k < 5; ++k) {
        CHECK(s.c[k] == doctest::Approx(c[k]).epsilon(TOL));
        CHECK(s.astar[k] == doctest::Approx(astar[k]).epsilon(TOL));
    }
}

TEST_CASE("pfqn_hst on an explicitly chosen non-bottleneck station") {
    // matlab: pfqn_hst([1;0.5],4,0,2)
    const std::vector<double> L = {1.0, 0.5};
    const auto s = pfqn_hst<double>(L, 4, 0.0, 1);
    CHECK(s.station == 1);
    CHECK(s.X == doctest::Approx(0.967741935484).epsilon(TOL));
    CHECK(s.U == doctest::Approx(0.483870967742).epsilon(TOL));
    CHECK(s.Q == doctest::Approx(0.838709677419).epsilon(TOL));
    CHECK(s.total == doctest::Approx(0.105376344086).epsilon(TOL));
    CHECK(s.worst == doctest::Approx(0.0731182795699).epsilon(TOL));
    const double c[4] = {-0.0172043010753, -0.0258064516129, -0.0301075268817, -0.0322580645161};
    const double astar[4] = {0.875, -1, -1, -1};
    for (int k = 0; k < 4; ++k) {
        CHECK(s.c[k] == doctest::Approx(c[k]).epsilon(TOL));
        CHECK(s.astar[k] == doctest::Approx(astar[k]).epsilon(TOL));
    }
}

TEST_CASE("pfqn_hst: Lemma 3.1, the certificate is Q_i(N) - Q_i(N-1)") {
    const std::vector<double> L = {0.9, 0.4, 0.2};
    Matrix<double> Lm(3, 1);
    for (std::size_t i = 0; i < 3; ++i) Lm(i, 0) = L[i];
    for (int N = 2; N <= 7; ++N) {
        for (std::size_t ist = 0; ist < 3; ++ist) {
            const auto s = pfqn_hst<double>(L, N, 0.7, ist);
            Matrix<double> Z(1, 1, 0.7);
            const std::vector<int> Nv(1, N), Nv1(1, N - 1);
            const double QN = line::pfqn::pfqn_mva(Lm, Nv, Z).QN(ist, 0);
            const double QN1 = line::pfqn::pfqn_mva(Lm, Nv1, Z).QN(ist, 0);
            CHECK(s.total == doctest::Approx(QN - QN1).epsilon(1e-9));
            CHECK(s.Q == doctest::Approx(QN).epsilon(1e-9));
        }
    }
}

TEST_CASE("pfqn_hst: the constrained optimum is feasible and tighter") {
    const std::vector<double> L = {0.9, 0.4, 0.2};
    for (int N = 2; N <= 6; ++N) {
        const auto s = pfqn_hst<double>(L, N, 1.0);
        // the (P1) optimum can never exceed the unconstrained certificate
        CHECK(s.worst <= s.total + 1e-12);
        CHECK(s.worst >= 0.0);
        // the box constraint |a_n| <= 1
        double feas = 0.0;
        for (std::size_t k = 0; k < s.astar.size(); ++k) {
            CHECK(std::abs(s.astar[k]) <= 1.0 + 1e-12);
            feas += s.p[k + 1] * s.astar[k];
        }
        // and the operational-consistency equality sum_n p_n a_n = 0
        CHECK(feas == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
        // the objective attained by astar is the reported worst case
        double obj = 0.0;
        for (std::size_t k = 0; k < s.c.size(); ++k) obj += s.c[k] * s.astar[k];
        CHECK(std::abs(obj) == doctest::Approx(s.worst).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_hst: the marginals are a proper distribution") {
    const std::vector<double> L = {0.6, 0.5, 0.4};
    const auto s = pfqn_hst<double>(L, 6, 2.0, 1);
    CHECK(s.Pgeq[0] == doctest::Approx(1.0).epsilon(1e-12));
    double tot = 0.0;
    for (std::size_t k = 0; k < s.p.size(); ++k) {
        CHECK(s.p[k] >= -1e-12);
        tot += s.p[k];
    }
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));
    for (std::size_t k = 1; k < s.Pgeq.size(); ++k) CHECK(s.Pgeq[k] <= s.Pgeq[k - 1] + 1e-12);
}

TEST_CASE("pfqn_hst rejects a degenerate station and a sub-unit population") {
    const std::vector<double> L = {0.9, 0.0};
    CHECK_THROWS(pfqn_hst<double>(L, 3, 0.0, 1));   // zero demand: no marginals
    CHECK_THROWS(pfqn_hst<double>(L, 0, 0.0, 0));   // no job to perturb
    CHECK_THROWS(pfqn_hst<double>(L, 3, 0.0, 5));   // station out of range
}

// ---------------------------------------------------------------------------
// pfqn_respt_ps_moments
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_respt_ps_moments exact route is the MVA residence time") {
    // One class: the closed terminal/PS system is product form, so the exact
    // Mitra-Morrison mean sojourn must be the MVA residence time at the PS
    // station. Two independent objects, one number.
    for (double S : {0.5, 0.25, 0.8}) {
        for (double Z : {3.0, 1.0}) {
            for (long N = 1; N <= 5; ++N) {
                const std::vector<double> Sv(1, S), Zv(1, Z);
                const std::vector<long> Nv(1, N);
                const auto r = pfqn_respt_ps_moments(Sv, Nv, Zv, ResptPsRoute::Exact);
                Matrix<double> Lm(1, 1, S), Zm(1, 1, Z);
                const std::vector<int> Nm(1, static_cast<int>(N));
                const double W = line::pfqn::pfqn_mva(Lm, Nm, Zm).CN(0, 0);
                CHECK(r.W[0] == doctest::Approx(W).epsilon(1e-9));
                CHECK(r.method[0] == ResptPsMethod::Exact);
            }
        }
    }
}

TEST_CASE("pfqn_respt_ps_moments: the second moment dominates the square of the first") {
    const std::vector<double> S = {0.2, 0.4, 0.3}, Z = {5, 7, 4};
    const std::vector<long> N = {2, 3, 1};
    const auto r = pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Exact);
    for (std::size_t j = 0; j < 3; ++j) {
        CHECK(r.W[j] > 0.0);
        CHECK(r.W2[j] > 0.0);
        // Var = E[W^2] - E[W]^2 must be non-negative
        CHECK(r.W2[j] >= r.W[j] * r.W[j] - 1e-9);
        CHECK(r.method[j] == ResptPsMethod::Exact);
    }
}

TEST_CASE("pfqn_respt_ps_moments asymptotic route approaches the exact one") {
    // Nexp = max_r Z_r/S_r is the expansion parameter, so driving the think
    // times up must close the gap. Three decades of Z on a two-class model.
    const std::vector<double> S = {0.2, 0.4};
    const std::vector<long> N = {2, 2};
    double prev = 1e9;
    for (double zscale : {10.0, 100.0, 1000.0}) {
        const std::vector<double> Z = {5 * zscale, 7 * zscale};
        const auto e = pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Exact);
        const auto a = pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Asymptotic);
        double err = 0.0;
        for (std::size_t j = 0; j < 2; ++j)
            err = std::max(err, std::abs(a.W[j] - e.W[j]) / e.W[j]);
        CHECK(err < prev);
        prev = err;
    }
    CHECK(prev < 1e-4);
}

TEST_CASE("pfqn_respt_ps_moments: auto picks exact below the state-space bound") {
    const std::vector<double> S = {0.5, 0.8}, Z = {4, 6};
    const std::vector<long> N = {3, 2};
    const auto au = pfqn_respt_ps_moments(S, N, Z);
    const auto ex = pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Exact);
    for (std::size_t j = 0; j < 2; ++j) {
        CHECK(au.method[j] == ResptPsMethod::Exact);
        CHECK(au.W[j] == doctest::Approx(ex.W[j]).epsilon(1e-12));
        CHECK(au.W2[j] == doctest::Approx(ex.W2[j]).epsilon(1e-12));
    }
    // the state space of the tagged class is prod_r (K_r + 1) with K = N - e_J
    CHECK(au.nstates[0] == doctest::Approx(3.0 * 3.0));
    CHECK(au.nstates[1] == doctest::Approx(4.0 * 2.0));
    CHECK(au.expansionParam == doctest::Approx(8.0).epsilon(1e-12));  // max(Z/S) = 6/0.8
}

TEST_CASE("pfqn_respt_ps_moments: an unpopulated class has no sojourn time") {
    const std::vector<double> S = {0.5, 0.8}, Z = {4, 6};
    const std::vector<long> N = {3, 0};
    const auto r = pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Exact);
    CHECK(std::isnan(r.W[1]));
    CHECK(std::isnan(r.W2[1]));
    CHECK(r.method[1] == ResptPsMethod::None);
    CHECK(!std::isnan(r.W[0]));
    CHECK(r.method[0] == ResptPsMethod::Exact);
}

TEST_CASE("pfqn_respt_ps_moments rejects a malformed model") {
    const std::vector<double> S = {0.5, 0.0}, Z = {4, 6};
    const std::vector<long> N = {3, 2};
    CHECK_THROWS(pfqn_respt_ps_moments(S, N, Z));                    // zero service time
    const std::vector<double> Sg = {0.5, 0.8}, Zbad = {4, 0};
    CHECK_THROWS(pfqn_respt_ps_moments(Sg, N, Zbad));                // no think time, N > 0
    const std::vector<long> Nbad = {3};
    CHECK_THROWS(pfqn_respt_ps_moments(Sg, Nbad, Z));                // arity mismatch
}

TEST_CASE("pfqn_respt_ps_moments: the asymptotic route refuses heavy usage") {
    // alpha = 1 - sum_r lambda_r/q_r <= 0 is outside normal usage. Short think
    // times relative to the service times drive it there.
    const std::vector<double> S = {0.9, 0.9}, Z = {1.0, 1.0};
    const std::vector<long> N = {3, 3};
    CHECK_THROWS(pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Asymptotic));
}

// ---------------------------------------------------------------------------
// pfqn_qlen_joint_moments
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_qlen_joint_moments tail route: the mean IS the MVA queue length") {
    Matrix<double> L(3, 1);
    L(0, 0) = 0.9; L(1, 0) = 0.4; L(2, 0) = 0.2;
    const std::vector<int> N(1, 5);
    const std::vector<double> Z(1, 1.0);
    const auto out = pfqn_qlen_joint_moments(L, N, Z);
    CHECK(out.route == "tail");
    Matrix<double> Zm(1, 1, 1.0);
    const auto mva = line::pfqn::pfqn_mva(L, N, Zm);
    REQUIRE(out.mean.size() == 3);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(out.mean[i] == doctest::Approx(mva.QN(i, 0)).epsilon(1e-8));
}

TEST_CASE("pfqn_qlen_joint_moments pmf route: the mean IS the MVA queue length") {
    Matrix<double> L(3, 2);
    L(0, 0) = 1.0; L(0, 1) = 0.5;
    L(1, 0) = 0.7; L(1, 1) = 0.9;
    L(2, 0) = 0.3; L(2, 1) = 0.4;
    const std::vector<int> N = {3, 2};
    const std::vector<double> Z = {0.0, 0.0};
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    pairs.push_back(std::make_pair(0, 0));
    pairs.push_back(std::make_pair(1, 1));
    const auto out = pfqn_qlen_joint_moments(L, N, Z, pairs);
    CHECK(out.route == "pmf");
    Matrix<double> Zm(1, 2, 0.0);
    const auto mva = line::pfqn::pfqn_mva(L, N, Zm);
    CHECK(out.mean[0] == doctest::Approx(mva.QN(0, 0)).epsilon(1e-8));
    CHECK(out.mean[1] == doctest::Approx(mva.QN(1, 1)).epsilon(1e-8));
}

TEST_CASE("pfqn_qlen_joint_moments pmf route with a think time") {
    Matrix<double> L(2, 2);
    L(0, 0) = 2.0; L(0, 1) = 1.0;
    L(1, 0) = 0.5; L(1, 1) = 1.5;
    const std::vector<int> N = {4, 3};
    const std::vector<double> Z = {0.5, 0.8};
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    pairs.push_back(std::make_pair(0, 0));
    pairs.push_back(std::make_pair(0, 1));
    const auto out = pfqn_qlen_joint_moments(L, N, Z, pairs);
    Matrix<double> Zm(1, 2);
    Zm(0, 0) = 0.5; Zm(0, 1) = 0.8;
    const auto mva = line::pfqn::pfqn_mva(L, N, Zm);
    CHECK(out.mean[0] == doctest::Approx(mva.QN(0, 0)).epsilon(1e-8));
    CHECK(out.mean[1] == doctest::Approx(mva.QN(0, 1)).epsilon(1e-8));
    // two coordinates at the SAME station: the cross-class covariance
    CHECK(out.cov.rows() == 2);
    CHECK(out.cov(0, 1) == doctest::Approx(out.cov(1, 0)).epsilon(1e-8));
}

TEST_CASE("pfqn_qlen_joint_moments: the survival array is a survival array") {
    Matrix<double> L(2, 1);
    L(0, 0) = 0.8; L(1, 0) = 0.3;
    const std::vector<int> N(1, 4);
    const std::vector<double> Z(1, 0.0);
    const auto out = pfqn_qlen_joint_moments(L, N, Z);
    // P(n >= 0 for all) = 1
    std::vector<std::size_t> origin(out.tail.order(), 0);
    CHECK(out.tail.at(origin) == doctest::Approx(1.0).epsilon(1e-9));
    // non-increasing along every mode, and never outside [0,1]
    for (std::size_t k = 0; k < out.tail.data.size(); ++k) {
        CHECK(out.tail.data[k] >= -1e-12);
        CHECK(out.tail.data[k] <= 1.0 + 1e-9);
    }
    for (std::size_t i = 0; i + 1 < out.tail.sz[0]; ++i) {
        std::vector<std::size_t> a(out.tail.order(), 0), b(out.tail.order(), 0);
        a[0] = i;
        b[0] = i + 1;
        CHECK(out.tail.at(b) <= out.tail.at(a) + 1e-12);
    }
}

TEST_CASE("pfqn_qlen_joint_moments: the covariance matrix is symmetric") {
    Matrix<double> L(3, 1);
    L(0, 0) = 0.9; L(1, 0) = 0.4; L(2, 0) = 0.2;
    const std::vector<int> N(1, 4);
    const std::vector<double> Z(1, 0.5);
    const auto out = pfqn_qlen_joint_moments(L, N, Z);
    for (std::size_t j = 0; j < 3; ++j) {
        CHECK(out.cov(j, j) >= -1e-12);  // a variance
        for (std::size_t l = 0; l < 3; ++l)
            CHECK(out.cov(j, l) == doctest::Approx(out.cov(l, j)).epsilon(1e-8));
    }
    // a closed single-class network with no delay conserves the population, so
    // the queue lengths are negatively correlated across stations
    const std::vector<double> Z0(1, 0.0);
    const auto tight = pfqn_qlen_joint_moments(L, N, Z0);
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t l = 0; l < 3; ++l)
            if (j != l) CHECK(tight.cov(j, l) <= 1e-9);
}

TEST_CASE("pfqn_qlen_joint_moments: the tail route is refused for several classes") {
    Matrix<double> L(2, 2, 0.5);
    const std::vector<int> N = {2, 2};
    const std::vector<double> Z = {0.0, 0.0};
    CHECK_THROWS(pfqn_qlen_joint_moments(L, N, Z,
                                         std::vector<std::pair<std::size_t, std::size_t>>(),
                                         QlenJointRoute::Tail, line::pfqn::QlenJointLgSource<double>(),
                                         line::pfqn::NcMethod::Exact, line::pfqn::NcOptions()));
}

TEST_CASE("pfqn_qlen_joint_moments rejects duplicate and out-of-range coordinates") {
    Matrix<double> L(2, 1, 0.5);
    const std::vector<int> N(1, 3);
    const std::vector<double> Z(1, 0.0);
    std::vector<std::pair<std::size_t, std::size_t>> dup;
    dup.push_back(std::make_pair(0, 0));
    dup.push_back(std::make_pair(0, 0));
    CHECK_THROWS(pfqn_qlen_joint_moments(L, N, Z, dup));
    std::vector<std::pair<std::size_t, std::size_t>> oor;
    oor.push_back(std::make_pair(5, 0));
    CHECK_THROWS(pfqn_qlen_joint_moments(L, N, Z, oor));
}

// ---------------------------------------------------------------------------
// pfqn_pff_delay
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_pff_delay is prod_r Z_r^n_r / n_r!") {
    std::vector<double> Z = {0.5, 0.8};
    std::vector<int> n = {2, 3};
    // 0.5^2/2 * 0.8^3/6
    CHECK(pfqn_pff_delay(Z, n) == doctest::Approx(0.125 * 0.512 / 6.0).epsilon(1e-12));
    n[0] = 0;
    n[1] = 0;
    CHECK(pfqn_pff_delay(Z, n) == doctest::Approx(1.0).epsilon(1e-15));
    n[0] = 1;
    CHECK(pfqn_pff_delay(Z, n) == doctest::Approx(0.5).epsilon(1e-12));
}

TEST_CASE("pfqn_pff_delay is zero when a populated class cannot think") {
    std::vector<double> Z = {0.5, 0.0};
    std::vector<int> n = {1, 2};
    CHECK(pfqn_pff_delay(Z, n) == doctest::Approx(0.0));
    // but a class with no jobs does not need a think time
    n[1] = 0;
    CHECK(pfqn_pff_delay(Z, n) == doctest::Approx(0.5).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// pfqn_xia
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_xia is the saturated bottleneck term on a one-station model") {
    // One station, s servers: B = {0}, no non-bottleneck term, scalefactor = s/L
    // log G ~ -log(0!) - N log(s/L) + s log(L s / L) - log(s!)
    const std::vector<double> L(1, 2.0), s(1, 3.0);
    const int N = 10;
    const double sf = 3.0 / 2.0;
    const double want = -static_cast<double>(N) * std::log(sf) + 3.0 * std::log(2.0 * sf) -
                        std::log(6.0);
    CHECK(pfqn_xia(L, N, s) == doctest::Approx(want).epsilon(1e-10));
}

TEST_CASE("pfqn_xia grows linearly in N at the bottleneck rate") {
    const std::vector<double> L = {1.0, 0.5, 0.25}, s = {1, 2, 1};
    const double a = pfqn_xia(L, 10, s);
    const double b = pfqn_xia(L, 20, s);
    const double c = pfqn_xia(L, 30, s);
    // -N log(scalefactor) is the only N dependence, so the increments are equal
    CHECK(b - a == doctest::Approx(c - b).epsilon(1e-10));
}

TEST_CASE("pfqn_xia rejects a degenerate model") {
    const std::vector<double> Lz = {0.0, 1.0}, s = {1, 1};
    CHECK_THROWS(pfqn_xia(Lz, 5, s));
    const std::vector<double> L = {1.0, 1.0}, sz = {1, 0};
    CHECK_THROWS(pfqn_xia(L, 5, sz));
    const std::vector<double> smis = {1.0};
    CHECK_THROWS(pfqn_xia(L, 5, smis));
}

// ---------------------------------------------------------------------------
// cache_gamma (the general access GRAPH)
// ---------------------------------------------------------------------------

TEST_CASE("cache_gamma follows the breadth-first path of the access graph") {
    // One user, one item, two lists. The graph is the linear chain
    // 0 -> 1 -> 2 with probability one on each edge.
    const std::size_t h = 2;
    Matrix<double> lam(1, h + 1);
    lam(0, 0) = 3.0;
    lam(0, 1) = 2.0;
    lam(0, 2) = 1.0;
    std::vector<Matrix<double>> lambda(1, lam);
    Matrix<double> Rm(h + 1, h + 1, 0.0);
    Rm(0, 1) = 1.0;
    Rm(1, 2) = 1.0;
    std::vector<std::vector<Matrix<double>>> R(1, std::vector<Matrix<double>>(1, Rm));
    const auto g = line::cache::cache_gamma(lambda, R);
    CHECK(g.u == 1);
    CHECK(g.n == 1);
    CHECK(g.h == 2);
    // column 0 walks to NODE 0, the trivial path: just the miss-node rate
    CHECK(g.gamma(0, 0) == doctest::Approx(3.0).epsilon(1e-12));
    // column 1 walks 0 -> 1: 3 * lambda(0)*R(0,1) = 3 * 3
    CHECK(g.gamma(0, 1) == doctest::Approx(9.0).epsilon(1e-12));
}

TEST_CASE("cache_gamma gives zero for an unreachable list") {
    const std::size_t h = 2;
    Matrix<double> lam(1, h + 1, 1.0);
    std::vector<Matrix<double>> lambda(1, lam);
    Matrix<double> Rm(h + 1, h + 1, 0.0);
    Rm(0, 1) = 1.0;  // node 2 is not reachable from node 0
    std::vector<std::vector<Matrix<double>>> R(1, std::vector<Matrix<double>>(1, Rm));
    const auto g = line::cache::cache_gamma(lambda, R);
    CHECK(g.gamma(0, 1) == doctest::Approx(1.0).epsilon(1e-12));  // node 1, reachable
    // and with h = 3 the unreachable node 2 yields exactly zero
    Matrix<double> lam3(1, 4, 1.0);
    std::vector<Matrix<double>> lambda3(1, lam3);
    Matrix<double> R3(4, 4, 0.0);
    R3(0, 1) = 1.0;
    std::vector<std::vector<Matrix<double>>> RR(1, std::vector<Matrix<double>>(1, R3));
    const auto g3 = line::cache::cache_gamma(lambda3, RR);
    CHECK(g3.gamma(0, 2) == doctest::Approx(0.0));
}

TEST_CASE("cache_gamma takes a shortest path where a list has two parents") {
    // 0 -> 1, 0 -> 2 and 1 -> 2: node 2 has two parents, which cache_gamma_lp
    // rejects outright and this resolves by breadth-first order (the direct
    // edge wins, being shorter).
    const std::size_t h = 2;
    Matrix<double> lam(1, h + 1);
    lam(0, 0) = 2.0; lam(0, 1) = 5.0; lam(0, 2) = 7.0;
    std::vector<Matrix<double>> lambda(1, lam);
    Matrix<double> Rm(h + 1, h + 1, 0.0);
    Rm(0, 1) = 0.5;
    Rm(0, 2) = 0.5;
    Rm(1, 2) = 1.0;
    std::vector<std::vector<Matrix<double>>> R(1, std::vector<Matrix<double>>(1, Rm));
    const auto g = line::cache::cache_gamma(lambda, R);
    // 0 -> 2 directly: 2 * (lambda(0) * R(0,2)) = 2 * 1 = 2
    CHECK(g.gamma(0, 1) == doctest::Approx(2.0 * (2.0 * 0.5)).epsilon(1e-12));
    // cache_gamma_lp refuses the same structure
    CHECK_THROWS(line::cache::cache_gamma_lp(lambda, R));
}

TEST_CASE("cache_gamma is exact in rational arithmetic") {
    Matrix<Rational> lam(1, 3);
    lam(0, 0) = Rational(1, 3);
    lam(0, 1) = Rational(1, 7);
    lam(0, 2) = Rational(0);
    std::vector<Matrix<Rational>> lambda(1, lam);
    Matrix<Rational> Rm(3, 3, Rational(0));
    Rm(0, 1) = Rational(1);
    Rm(1, 2) = Rational(1);
    std::vector<std::vector<Matrix<Rational>>> R(1, std::vector<Matrix<Rational>>(1, Rm));
    const auto g = line::cache::cache_gamma(lambda, R);
    CHECK(g.gamma(0, 0) == Rational(1, 3));
    CHECK(g.gamma(0, 1) == Rational(1, 9));
}
