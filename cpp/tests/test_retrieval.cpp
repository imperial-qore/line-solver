/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Delayed-hit (list-based) cache retrieval. Oracles:
 *   1. The balance law pi_{i,0} + sum_s phi_{s,i} + sum_j pi_{i,j} = 1, asserted
 *      EXACTLY (as an identity between rationals) for retrieval_metrics and
 *      retrieval_mva, not merely to a tolerance.
 *   2. retrieval_mva and retrieval_metrics are two independent routes to the
 *      same quantities -- an MVA-style recursion over item subsets against a
 *      ratio of normalizing constants -- so in exact arithmetic they must
 *      agree as values.
 *   3. With no retrieval system at all (eta = 0) the model degenerates to the
 *      plain list-based cache, so retrieval_nc must equal cache_erec and
 *      retrieval_metrics must equal cache_prob_erec, again exactly.
 *   4. MATLAB reference values from matlab/src/api/retrieval/*.m.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/retrieval/retrieval_fpi.h"
#include "line/api/retrieval/retrieval_fpi_latency.h"
#include "line/api/retrieval/retrieval_metrics.h"
#include "line/api/retrieval/retrieval_mva.h"
#include "line/api/retrieval/retrieval_nc.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::cache::cache_erec;
using line::cache::cache_prob_erec;
using line::retrieval::FpiOptions;
using line::retrieval::retrieval_fpi;
using line::retrieval::retrieval_fpi_latency;
using line::retrieval::retrieval_metrics;
using line::retrieval::retrieval_mva;
using line::retrieval::retrieval_nc;
using line::retrieval::RetrievalStationPH;
using line::retrieval::RetrievalStationType;

namespace {

// The MATLAB reference instance: 3 items, one list of capacity 1, one PS
// station. lambda = [1 2 1], eta = [1/2 1/3; 1/4 1/5; 1/6 1/7],
// gamma = [1/2; 1/3; 1/4].
template <class T>
std::vector<T> lambda_3() {
    return {line::num_traits<T>::from_int(1), line::num_traits<T>::from_int(2),
            line::num_traits<T>::from_int(1)};
}

template <class T>
Matrix<T> eta_3x2() {
    Matrix<T> e(3, 2);
    e(0, 0) = line::num_traits<T>::from_rational(1, 2);
    e(0, 1) = line::num_traits<T>::from_rational(1, 3);
    e(1, 0) = line::num_traits<T>::from_rational(1, 4);
    e(1, 1) = line::num_traits<T>::from_rational(1, 5);
    e(2, 0) = line::num_traits<T>::from_rational(1, 6);
    e(2, 1) = line::num_traits<T>::from_rational(1, 7);
    return e;
}

template <class T>
Matrix<T> gamma_3x1() {
    Matrix<T> g(3, 1);
    g(0, 0) = line::num_traits<T>::from_rational(1, 2);
    g(1, 0) = line::num_traits<T>::from_rational(1, 3);
    g(2, 0) = line::num_traits<T>::from_rational(1, 4);
    return g;
}

// The second reference instance: 4 items, two lists of capacities [2 1], two
// PS stations.
template <class T>
std::vector<T> lambda_4() {
    return {line::num_traits<T>::from_int(1), line::num_traits<T>::from_int(2),
            line::num_traits<T>::from_int(1), line::num_traits<T>::from_int(3)};
}

template <class T>
Matrix<T> eta_4x3() {
    Matrix<T> e(4, 3);
    const int den[4][3] = {{2, 3, 4}, {5, 6, 7}, {8, 9, 10}, {11, 12, 13}};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j) e(i, j) = line::num_traits<T>::from_rational(1, den[i][j]);
    return e;
}

template <class T>
Matrix<T> gamma_4x2() {
    Matrix<T> g(4, 2);
    const int den[4][2] = {{2, 3}, {4, 5}, {6, 7}, {8, 9}};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j) g(i, j) = line::num_traits<T>::from_rational(1, den[i][j]);
    return g;
}

/** Relative discrepancy against a MATLAB-printed double. */
double rel(double got, double want) {
    const double d = std::fabs(got - want);
    return std::fabs(want) > 1e-30 ? d / std::fabs(want) : d;
}

}  // namespace

TEST_CASE("retrieval_metrics satisfies the balance law exactly") {
    // pi_{i,0} + sum_s phi_{s,i} + sum_j pi_{i,j} = 1 for every item. This is
    // an identity between rationals, so == is the right assertion.
    const std::vector<int> m{1};
    const auto r = retrieval_metrics(m, lambda_3<Rational>(), eta_3x2<Rational>(),
                                     gamma_3x1<Rational>());
    for (std::size_t i = 0; i < 3; ++i) {
        Rational tot = r.pmiss[i];
        for (std::size_t s = 0; s < r.pdh.rows(); ++s) tot += r.pdh(s, i);
        for (std::size_t j = 0; j < r.phit.rows(); ++j) tot += r.phit(j, i);
        CHECK(tot == Rational(1));
    }

    const std::vector<int> m2{2, 1};
    const auto r2 = retrieval_metrics(m2, lambda_4<Rational>(), eta_4x3<Rational>(),
                                      gamma_4x2<Rational>());
    for (std::size_t i = 0; i < 4; ++i) {
        Rational tot = r2.pmiss[i];
        for (std::size_t s = 0; s < r2.pdh.rows(); ++s) tot += r2.pdh(s, i);
        for (std::size_t j = 0; j < r2.phit.rows(); ++j) tot += r2.phit(j, i);
        CHECK(tot == Rational(1));
    }
}

TEST_CASE("retrieval_mva reproduces retrieval_metrics exactly") {
    // Two independent exact routes to the same metrics: the subset recursion
    // and the ratio of normalizing constants. In rationals they must agree as
    // values, which is the strongest available check on either one.
    const std::vector<int> m{1};
    const auto a = retrieval_metrics(m, lambda_3<Rational>(), eta_3x2<Rational>(),
                                     gamma_3x1<Rational>());
    const auto b = retrieval_mva(m, lambda_3<Rational>(), eta_3x2<Rational>(),
                                 gamma_3x1<Rational>());
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(a.pmiss[i] == b.pmiss[i]);
        for (std::size_t j = 0; j < a.phit.rows(); ++j) CHECK(a.phit(j, i) == b.phit(j, i));
        for (std::size_t s = 0; s < a.pdh.rows(); ++s) CHECK(a.pdh(s, i) == b.pdh(s, i));
    }

    const std::vector<int> m2{2, 1};
    const auto c = retrieval_metrics(m2, lambda_4<Rational>(), eta_4x3<Rational>(),
                                     gamma_4x2<Rational>());
    const auto d = retrieval_mva(m2, lambda_4<Rational>(), eta_4x3<Rational>(),
                                 gamma_4x2<Rational>());
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(c.pmiss[i] == d.pmiss[i]);
        for (std::size_t j = 0; j < c.phit.rows(); ++j) CHECK(c.phit(j, i) == d.phit(j, i));
        for (std::size_t s = 0; s < c.pdh.rows(); ++s) CHECK(c.pdh(s, i) == d.pdh(s, i));
    }

    // retrieval_mva also satisfies the balance law exactly.
    for (std::size_t i = 0; i < 4; ++i) {
        Rational tot = d.pmiss[i];
        for (std::size_t s = 0; s < d.pdh.rows(); ++s) tot += d.pdh(s, i);
        for (std::size_t j = 0; j < d.phit.rows(); ++j) tot += d.phit(j, i);
        CHECK(tot == Rational(1));
    }
}

TEST_CASE("retrieval_nc hand values on the smallest instances") {
    // One item, one list of capacity one. The only surviving branch of the
    // recurrence is "item 1 stored in list 1", so E = gamma(1,1).
    Matrix<Rational> g(1, 1);
    g(0, 0) = Rational(3, 7);
    Matrix<Rational> e(1, 2);
    e(0, 0) = Rational(1, 2);
    e(0, 1) = Rational(1, 5);
    const std::vector<Rational> lam{Rational(2)};
    CHECK(retrieval_nc(std::vector<int>{0}, std::vector<int>{1}, lam, e, g) == Rational(3, 7));

    // Empty cache: the item is either outside or being fetched, so
    // E(0,[0]) = (1 + lambda eta_0) + lambda eta_1 (v_1+1) = 2 + 2/5.
    CHECK(retrieval_nc(std::vector<int>{0}, std::vector<int>{0}, lam, e, g) == Rational(12, 5));

    // A capacity larger than the item count admits no placement at all.
    CHECK(retrieval_nc(std::vector<int>{0}, std::vector<int>{2}, lam, e, g) == Rational(0));

    // The first PS moment doubles the fetching branch, (v_1+1) = 2:
    // E(1_1,[0]) = 2 + 2 * 2/5.
    CHECK(retrieval_nc(std::vector<int>{1}, std::vector<int>{0}, lam, e, g) == Rational(14, 5));
}

TEST_CASE("retrieval with no retrieval system degenerates to the plain cache") {
    // eta = 0 kills every fetching branch of the recurrence, leaving exactly
    // the cache_erec recursion, and the metrics must then collapse onto
    // cache_prob_erec. Both are exact identities in rationals.
    Matrix<Rational> gamma(3, 2);
    const int den[3][2] = {{2, 3}, {4, 5}, {6, 7}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) gamma(i, j) = Rational(1, den[i][j]);
    const Matrix<Rational> eta(3, 3, Rational(0));  // one IS + two PS columns, all zero
    const std::vector<Rational> lam{Rational(1), Rational(2), Rational(3)};

    for (const std::vector<int>& m : std::vector<std::vector<int>>{{1, 1}, {2, 1}, {1, 0}}) {
        CHECK(retrieval_nc(std::vector<int>{0, 0}, m, lam, eta, gamma) == cache_erec(gamma, m));
        const auto r = retrieval_metrics(m, lam, eta, gamma);
        const Matrix<Rational> p = cache_prob_erec(gamma, m);
        for (std::size_t i = 0; i < 3; ++i) {
            CHECK(r.pmiss[i] == p(i, 0));
            for (std::size_t j = 0; j < 2; ++j) CHECK(r.phit(j, i) == p(i, 1 + j));
            // no fetching at all, so every delayed-hit probability vanishes
            for (std::size_t s = 0; s < r.pdh.rows(); ++s) CHECK(r.pdh(s, i) == Rational(0));
        }
    }
}

TEST_CASE("retrieval_metrics matches MATLAB retrieval_metrics.m") {
    // Reference from matlab -batch retrieval_metrics/retrieval_nc. Both sides
    // evaluate the same finite recursion, so the only discrepancy is the
    // double rounding of a few dozen operations: 1e-12 relative is already
    // three orders looser than that.
    const std::vector<int> m{1};
    const double E0 = 2.99292328042328, E1 = 4.12056878306878;
    CHECK(rel(retrieval_nc(std::vector<int>{0}, m, lambda_3<double>(), eta_3x2<double>(),
                           gamma_3x1<double>()),
              E0) < 1e-12);
    CHECK(rel(retrieval_nc(std::vector<int>{1}, m, lambda_3<double>(), eta_3x2<double>(),
                           gamma_3x1<double>()),
              E1) < 1e-12);

    const auto r = retrieval_metrics(m, lambda_3<double>(), eta_3x2<double>(), gamma_3x1<double>());
    const double pmiss[3] = {0.304554394183811, 0.371909044704219, 0.521600777848982};
    const double phit[3] = {0.425209378383754, 0.272689103484852, 0.302101518131395};
    const double pdh[2][3] = {{0.152277197091906, 0.185954522352109, 0.0869334629748304},
                              {0.11795903034053, 0.16944732945882, 0.0893642410447926}};
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(rel(r.pmiss[i], pmiss[i]) < 1e-12);
        CHECK(rel(r.phit(0, i), phit[i]) < 1e-12);
        for (std::size_t s = 0; s < 2; ++s) CHECK(rel(r.pdh(s, i), pdh[s][i]) < 1e-12);
    }

    // Two lists, two PS stations.
    const std::vector<int> m2{2, 1};
    const auto r2 = retrieval_metrics(m2, lambda_4<double>(), eta_4x3<double>(), gamma_4x2<double>());
    const double pmiss2[4] = {0.0624211621021485, 0.11830444687689, 0.173253749176287,
                              0.227891723983904};
    const double phit2[2][4] = {
        {0.623798054638859, 0.504088239645337, 0.492478182527952, 0.379635523187852},
        {0.246157857648332, 0.257049448565133, 0.276035558155953, 0.220757135630582}};
    const double pdh2[3][4] = {
        {0.0312105810510742, 0.047321778750756, 0.0216567186470359, 0.0621522883592465},
        {0.0208070540340495, 0.03943481562563, 0.0192504165751431, 0.0569729309959759},
        {0.0156052905255371, 0.0338012705362543, 0.0173253749176287, 0.0525903978424393}};
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(rel(r2.pmiss[i], pmiss2[i]) < 1e-12);
        for (std::size_t j = 0; j < 2; ++j) CHECK(rel(r2.phit(j, i), phit2[j][i]) < 1e-12);
        for (std::size_t s = 0; s < 3; ++s) CHECK(rel(r2.pdh(s, i), pdh2[s][i]) < 1e-12);
    }
}

TEST_CASE("retrieval_mva matches MATLAB retrieval_mva.m") {
    // Same reference values: MATLAB's own retrieval_mva agrees with its
    // retrieval_metrics to 5.6e-17 on this instance.
    const std::vector<int> m{1};
    const auto r = retrieval_mva(m, lambda_3<double>(), eta_3x2<double>(), gamma_3x1<double>());
    const double pmiss[3] = {0.304554394183811, 0.371909044704219, 0.521600777848982};
    const double phit[3] = {0.425209378383754, 0.272689103484852, 0.302101518131395};
    const double pdh0[3] = {0.152277197091905, 0.185954522352109, 0.0869334629748304};
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(rel(r.pmiss[i], pmiss[i]) < 1e-12);
        CHECK(rel(r.phit(0, i), phit[i]) < 1e-12);
        CHECK(rel(r.pdh(0, i), pdh0[i]) < 1e-12);
    }
}

TEST_CASE("retrieval_mva leaves the saturated cache with an unattributed hit") {
    // DOCUMENTED MATLAB DEFECT, reproduced deliberately. When the lists can
    // hold every item, retrieval_mva takes its boundary branch and records
    // pihit = 1 without attributing it to any list, so phit comes back all
    // zeros; retrieval_metrics attributes it and returns 1. MATLAB does
    // exactly this (BND_MVA phit = 0 0, BND_MET phit = 1 1), and the port is
    // faithful to it rather than inventing an attribution rule.
    const std::vector<int> m{2};
    const std::vector<Rational> lam{Rational(1), Rational(2)};
    Matrix<Rational> eta(2, 2);
    eta(0, 0) = Rational(1, 2);
    eta(0, 1) = Rational(1, 3);
    eta(1, 0) = Rational(1, 4);
    eta(1, 1) = Rational(1, 5);
    Matrix<Rational> g(2, 1);
    g(0, 0) = Rational(1, 2);
    g(1, 0) = Rational(1, 3);

    const auto met = retrieval_metrics(m, lam, eta, g);
    const auto mva = retrieval_mva(m, lam, eta, g);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(met.pmiss[i] == Rational(0));
        CHECK(met.phit(0, i) == Rational(1));
        CHECK(mva.pmiss[i] == Rational(0));
        CHECK(mva.phit(0, i) == Rational(0));  // the defect
    }
}

TEST_CASE("retrieval_fpi matches MATLAB retrieval_fpi.m") {
    // Driven to tol = 1e-12 with the same iterate ordering as MATLAB, so the
    // two fixed points agree to the stopping tolerance; 1e-9 absolute is a
    // decade of slack on that.
    FpiOptions opt;
    opt.max_iter = 100000;
    opt.tol = 1e-12;

    const std::vector<int> m{1};
    const auto r = retrieval_fpi(m, lambda_3<double>(), eta_3x2<double>(), gamma_3x1<double>(), opt);
    REQUIRE(r.converged);
    const double pmiss[3] = {0.316070368254084, 0.354788366459017, 0.506597121070067};
    const double phit[3] = {0.392198191863786, 0.293494465780216, 0.314307342356368};
    const double pdh[2][3] = {{0.158035184127042, 0.177394183229509, 0.0844328535116778},
                              {0.133696255755088, 0.174322984531258, 0.0946626830618872}};
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(std::fabs(r.pmiss[i] - pmiss[i]) < 1e-9);
        CHECK(std::fabs(r.phit(0, i) - phit[i]) < 1e-9);
        for (std::size_t s = 0; s < 2; ++s) CHECK(std::fabs(r.pdh(s, i) - pdh[s][i]) < 1e-9);
    }

    const std::vector<int> m2{2, 1};
    const auto r2 =
        retrieval_fpi(m2, lambda_4<double>(), eta_4x3<double>(), gamma_4x2<double>(), opt);
    REQUIRE(r2.converged);
    const double pmiss2[4] = {0.0715267021157106, 0.119841666305603, 0.177375315534569,
                              0.198256355153419};
    for (std::size_t i = 0; i < 4; ++i) CHECK(std::fabs(r2.pmiss[i] - pmiss2[i]) < 1e-9);
    const double phit2[2][4] = {
        {0.595716468226224, 0.499055961523485, 0.492428670670093, 0.412798899578572},
        {0.250500518484157, 0.251825391001482, 0.266229896528635, 0.231444193984921}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(std::fabs(r2.phit(j, i) - phit2[j][i]) < 1e-9);
}

TEST_CASE("retrieval_fpi obeys the balance law at its fixed point") {
    // The heuristic is not exact, but its fixed point does satisfy the same
    // balance identity, so the residual measures convergence rather than the
    // approximation error. Held to the stopping tolerance, not to the distance
    // from the exact solution (which is 1e-2 here and is the point of the
    // exact routines).
    FpiOptions opt;
    opt.max_iter = 100000;
    opt.tol = 1e-14;
    const std::vector<int> m{2, 1};
    const auto r = retrieval_fpi(m, lambda_4<Real50>(), eta_4x3<Real50>(), gamma_4x2<Real50>(), opt);
    REQUIRE(r.converged);
    for (std::size_t i = 0; i < 4; ++i) {
        Real50 tot = r.pmiss[i];
        for (std::size_t s = 0; s < r.pdh.rows(); ++s) tot += r.pdh(s, i);
        for (std::size_t j = 0; j < r.phit.rows(); ++j) tot += r.phit(j, i);
        CHECK(std::fabs(static_cast<double>(tot) - 1.0) < 1e-12);
    }

    // The heuristic stays within a few percent of the exact metrics.
    const auto exact = retrieval_metrics(m, lambda_4<double>(), eta_4x3<double>(),
                                         gamma_4x2<double>());
    const auto approx = retrieval_fpi(m, lambda_4<double>(), eta_4x3<double>(),
                                      gamma_4x2<double>(), opt);
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(std::fabs(approx.pmiss[i] - exact.pmiss[i]) < 0.05);
}

TEST_CASE("retrieval_fpi_latency matches MATLAB retrieval_fpi_latency.m") {
    // Two items, one list of capacity one, one PS station with exponential
    // service of rates 3 and 4, routed outside -> station -> outside. Run at
    // the MATLAB defaults (1000 iterations, tol 1e-6), so the comparison is
    // held at 1e-9 absolute: the fixed point itself is only located to 1e-6,
    // but both sides take the identical iterate sequence to it.
    const std::vector<int> m{1};
    const std::vector<double> lambda{1.0, 2.0};
    Matrix<double> gamma(2, 1);
    gamma(0, 0) = 0.5;
    gamma(1, 0) = 1.0 / 3.0;

    RetrievalStationPH<double> st;
    st.type = RetrievalStationType::PS;
    st.alpha = Matrix<double>(2, 1, 1.0);
    st.sub.push_back(Matrix<double>{{-3.0}});
    st.sub.push_back(Matrix<double>{{-4.0}});

    Matrix<double> Ri(2, 2, 0.0);
    Ri(0, 1) = 1.0;
    Ri(1, 0) = 1.0;
    const std::vector<Matrix<double>> R{Ri, Ri};

    const auto r = retrieval_fpi_latency(m, lambda, gamma, std::vector<RetrievalStationPH<double>>{st}, R);
    CHECK(std::fabs(r.Z - 0.301127289784926) < 1e-9);
    CHECK(std::fabs(r.d[0] - 0.0416156308318846) < 1e-9);
    CHECK(std::fabs(r.d[1] - 0.10152703841542) < 1e-9);
    CHECK(std::fabs(r.phi[0] - 0.124805373741914) < 1e-9);
    CHECK(std::fabs(r.phi[1] - 0.203003427273137) < 1e-9);
    CHECK(std::fabs(r.pi0[0] - 0.31123442510907) < 1e-9);
    CHECK(std::fabs(r.pi0[1] - 0.36095739414325) < 1e-9);
}

TEST_CASE("retrieval_fpi_latency rejects the unsupported station configurations") {
    const std::vector<int> m{1};
    const std::vector<double> lambda{1.0, 2.0};
    Matrix<double> gamma(2, 1, 0.5);
    Matrix<double> Ri(2, 2, 0.0);
    Ri(0, 1) = 1.0;
    Ri(1, 0) = 1.0;
    const std::vector<Matrix<double>> R{Ri, Ri};

    // FCFS with class-dependent rates has no single-exponential sojourn.
    RetrievalStationPH<double> fcfs;
    fcfs.type = RetrievalStationType::FCFS;
    fcfs.alpha = Matrix<double>(2, 1, 1.0);
    fcfs.sub.push_back(Matrix<double>{{-3.0}});
    fcfs.sub.push_back(Matrix<double>{{-4.0}});
    CHECK_THROWS_AS(
        retrieval_fpi_latency(m, lambda, gamma, std::vector<RetrievalStationPH<double>>{fcfs}, R),
        line::UnsupportedError);

    // SIRO with two phases is rejected before the rates are even looked at.
    RetrievalStationPH<double> siro;
    siro.type = RetrievalStationType::SIRO;
    siro.alpha = Matrix<double>(2, 2, 0.5);
    siro.sub.push_back(Matrix<double>{{-3.0, 1.0}, {0.0, -3.0}});
    siro.sub.push_back(Matrix<double>{{-3.0, 1.0}, {0.0, -3.0}});
    CHECK_THROWS_AS(
        retrieval_fpi_latency(m, lambda, gamma, std::vector<RetrievalStationPH<double>>{siro}, R),
        line::UnsupportedError);
}

TEST_CASE("retrieval input validation") {
    Matrix<double> g(2, 1, 0.5);
    Matrix<double> e(2, 2, 0.1);
    const std::vector<double> lam{1.0};
    CHECK_THROWS_AS(retrieval_nc(std::vector<int>{0}, std::vector<int>{1}, lam, e, g),
                    line::InputError);
    CHECK_THROWS_AS(retrieval_metrics(std::vector<int>{1}, lam, e, g), line::InputError);
    CHECK_THROWS_AS(retrieval_mva(std::vector<int>{1}, lam, e, g), line::InputError);
}
