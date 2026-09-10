/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The inexact half of the cache family: the fixed-point multipliers
 * (cache_xi_fp, cache_xi_iter), the saddle-point approximation (cache_spm and
 * its probability/miss front-ends) and the h-LRU TTL approximation. All of
 * these are gated on transcendental arithmetic, so only double and Real50 are
 * instantiated -- an attempt at line::Rational is a compile error, which is
 * the intended behaviour and cannot be exercised from a runtime test.
 *
 * Oracles, in order of strength:
 *   1. Values from the MATLAB reference, pinned at 1e-12 relative. Every
 *      number below was produced by matlab/src/api/cache on the same model.
 *   2. cache_xi_fp and cache_xi_iter are two unrelated iterations for the same
 *      fixed point (Gauss-Seidel bisection versus alternating substitution),
 *      so they must agree; and the fixed point itself is checked by
 *      substituting back into the capacity constraint.
 *   3. Exact coincidence points: on a symmetric single-list model the
 *      decoupling approximation is exact, so cache_prob_fpi must reproduce
 *      cache_prob_erec; and when the item count equals the total capacity
 *      cache_spm falls back to cache_erec.
 *   4. Closed forms: for h = 1 the TTL fixed point is the Che approximation
 *      1 - exp(-lam T), which on a symmetric model gives T in closed form.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_miss.h"
#include "line/api/cache/cache_miss_fpi.h"
#include "line/api/cache/cache_miss_spm.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_prob_fpi.h"
#include "line/api/cache/cache_prob_spm.h"
#include "line/api/cache/cache_spm.h"
#include "line/api/cache/cache_t_hlru.h"
#include "line/api/cache/cache_ttl_hlru.h"
#include "line/api/cache/cache_xi_fp.h"
#include "line/api/cache/cache_xi_iter.h"

using line::Matrix;
using line::Real50;
using line::cache::cache_erec;
using line::cache::cache_miss;
using line::cache::cache_miss_fpi;
using line::cache::cache_miss_spm;
using line::cache::cache_prob_erec;
using line::cache::cache_prob_fpi;
using line::cache::cache_prob_spm;
using line::cache::cache_spm;
using line::cache::cache_t_hlru;
using line::cache::cache_ttl_hlru;
using line::cache::cache_xi_fp;
using line::cache::cache_xi_iter;

namespace {

constexpr double MATLAB_TOL = 1e-12;
constexpr double TOL = 1e-9;

/** Six items, two lists, second column exactly half the first. */
template <class T>
Matrix<T> gamma_6x2() {
    static const int num[6] = {50, 40, 30, 20, 10, 5};
    Matrix<T> g(6, 2);
    for (int i = 0; i < 6; ++i) {
        g(i, 0) = line::num_traits<T>::from_rational(num[i], 100);
        g(i, 1) = line::num_traits<T>::from_rational(num[i], 200);
    }
    return g;
}

/** First-list request rates of the same six items, as a (1 x 6) row. */
template <class T>
Matrix<T> lambda_1x6() {
    static const int num[6] = {50, 40, 30, 20, 10, 5};
    Matrix<T> l(1, 6);
    for (int i = 0; i < 6; ++i) l(0, i) = line::num_traits<T>::from_rational(num[i], 100);
    return l;
}

const std::vector<int> M21{2, 1};

}  // namespace

TEST_CASE("cache_xi_fp and cache_xi_iter solve the same capacity fixed point") {
    const Matrix<double> g = gamma_6x2<double>();
    const auto fp = cache_xi_fp(g, M21);
    const std::vector<double> it = cache_xi_iter(g, M21);

    // Two unrelated iterations, same root.
    for (std::size_t l = 0; l < 2; ++l)
        CHECK(fp.xi[l] == doctest::Approx(it[l]).epsilon(TOL));

    // Substitute back: sum_k gamma(k,l) xi(l) / (1 + sum_s gamma(k,s) xi(s)) = m(l).
    for (const std::vector<double>* xi : {&fp.xi, &it}) {
        for (std::size_t l = 0; l < 2; ++l) {
            double occ = 0.0;
            for (std::size_t k = 0; k < 6; ++k) {
                double S = 0.0;
                for (std::size_t s = 0; s < 2; ++s) S += g(k, s) * (*xi)[s];
                occ += g(k, l) * (*xi)[l] / (1.0 + S);
            }
            CHECK(occ == doctest::Approx(static_cast<double>(M21[l])).epsilon(1e-9));
        }
    }
}

TEST_CASE("cache_xi_fp and cache_xi_iter match the MATLAB reference") {
    const auto fp = cache_xi_fp(gamma_6x2<double>(), M21);
    const std::vector<double> it = cache_xi_iter(gamma_6x2<double>(), M21);
    // MATLAB cache_xi_fp / cache_xi_iter on this model.
    CHECK(fp.xi[0] == doctest::Approx(3.2968685280906236).epsilon(MATLAB_TOL));
    CHECK(fp.xi[1] == doctest::Approx(3.2968685280906236).epsilon(MATLAB_TOL));
    CHECK(it[0] == doctest::Approx(3.296868528089659).epsilon(MATLAB_TOL));
    CHECK(it[1] == doctest::Approx(3.2968685280891297).epsilon(MATLAB_TOL));

    // Structural identity: the second column is half the first and m(2) is
    // half m(1), so the two capacity constraints coincide and xi(1) = xi(2).
    CHECK(fp.xi[0] == doctest::Approx(fp.xi[1]).epsilon(1e-14));

    // Real50 must land on the same root.
    const auto fpr = cache_xi_fp(gamma_6x2<Real50>(), M21);
    CHECK(static_cast<double>(fpr.xi[0]) == doctest::Approx(fp.xi[0]).epsilon(TOL));
}

TEST_CASE("cache_prob_fpi is exact on a symmetric single-list cache") {
    // Five identical items, two slots: by symmetry every item is cached with
    // probability m/n = 2/5, and the decoupling approximation reproduces it.
    Matrix<double> g(5, 1, 0.3);
    const std::vector<int> m{2};
    const Matrix<double> approx = cache_prob_fpi(g, m);
    const Matrix<double> exact = cache_prob_erec(g, m);
    for (std::size_t i = 0; i < 5; ++i) {
        CHECK(exact(i, 1) == doctest::Approx(0.4).epsilon(1e-14));
        CHECK(approx(i, 1) == doctest::Approx(exact(i, 1)).epsilon(1e-12));
        CHECK(approx(i, 0) + approx(i, 1) == doctest::Approx(1.0).epsilon(1e-14));
    }
    // MATLAB cache_prob_fpi on the same model.
    CHECK(approx(0, 1) == doctest::Approx(0.40000000000000058).epsilon(MATLAB_TOL));
}

TEST_CASE("cache_prob_fpi reproduces the duplicated-column reference defect") {
    // Documented in the header: for h > 1 the references write the aggregate
    // hit probability into every list column, so the columns are equal and the
    // row does not sum to one. MATLAB gives row 1 = [0.28796440700390641,
    // 0.71203559299609365, 0.71203559299609365].
    const Matrix<double> p = cache_prob_fpi(gamma_6x2<double>(), M21);
    CHECK(p(0, 0) == doctest::Approx(0.28796440700390641).epsilon(MATLAB_TOL));
    CHECK(p(0, 1) == doctest::Approx(0.71203559299609365).epsilon(MATLAB_TOL));
    CHECK(p(0, 2) == p(0, 1));
    CHECK(p(0, 0) + p(0, 1) == doctest::Approx(1.0).epsilon(1e-14));
}

TEST_CASE("cache_spm answers the full cache exactly instead of iterating") {
    // Three items, three slots: every item is cached, so the capacity equations
    // force every multiplier to infinity and there is no interior saddle to
    // find. `cache_spm.m` used to call cache_xi_iter here anyway, whose bracket
    // doubling and Gauss-Seidel sweep never terminate --
    // `cache_spm([.5 .25;.4 .2;.3 .15],[2 1])` did not return within two
    // minutes -- and this port capped the sweep and raised NumericError. The
    // reference now short-circuits the case in all four codebases: n == sum(m)
    // takes Z from the exact recursion and reports the multipliers' limit,
    // `Z=cache_erec(gamma,m); lZ=log(Z); xi=inf(1,h)`. So the degenerate cache
    // is answered EXACTLY rather than approximated or refused.
    Matrix<double> g(3, 2);
    g(0, 0) = 0.5;
    g(0, 1) = 0.25;
    g(1, 0) = 0.4;
    g(1, 1) = 0.2;
    g(2, 0) = 0.3;
    g(2, 1) = 0.15;
    const std::vector<int> m{2, 1};
    const auto s = cache_spm(g, m);
    // Hand computation: three placements of {1,2,3} into a 2-slot and a
    // 1-slot list, each worth 0.03 (.5*.4*.15 = .5*.3*.2 = .4*.3*.25), times
    // the 2!*1! list multiplicity, so E = 0.18.
    CHECK(s.Z == doctest::Approx(0.18).epsilon(1e-14));
    CHECK(s.lZ == doctest::Approx(std::log(0.18)).epsilon(1e-14));
    for (std::size_t l = 0; l < 2; ++l) CHECK(std::isinf(s.xi[l]));
    // The same constant the short circuit took, reached directly.
    CHECK(cache_erec(g, m) == doctest::Approx(0.18).epsilon(1e-14));
}

TEST_CASE("cache_spm matches the MATLAB reference and is consistent with cache_erec") {
    const auto s = cache_spm(gamma_6x2<double>(), M21);
    // MATLAB cache_spm: Z = 0.90905945820570455, lZ = -0.095344776376504514.
    CHECK(s.Z == doctest::Approx(0.90905945820570455).epsilon(MATLAB_TOL));
    CHECK(s.lZ == doctest::Approx(-0.095344776376504514).epsilon(MATLAB_TOL));
    CHECK(s.lZ == doctest::Approx(std::log(s.Z)).epsilon(1e-13));

    // The exact constant is 0.8025; the saddle point is 13% high on a cache
    // this small, which is the documented accuracy of the approximation and
    // not an implementation error (the value agrees with MATLAB to 1e-15).
    const double exact = cache_erec(gamma_6x2<double>(), M21);
    CHECK(exact == doctest::Approx(0.8025).epsilon(1e-15));
    CHECK(std::fabs(s.Z / exact - 1.0) < 0.14);

    const auto sr = cache_spm(gamma_6x2<Real50>(), M21);
    CHECK(static_cast<double>(sr.lZ) == doctest::Approx(s.lZ).epsilon(TOL));
}

TEST_CASE("cache_prob_spm matches MATLAB and stays a distribution at a zero-capacity list") {
    // Column 2 asks for cache_spm at capacity [2,0]. A list with no capacity
    // has xi = 0, which is a BOUNDARY of the Laplace integral rather than a
    // direction of it: kept in the expansion its -sum_l log(sqrt(xi_l))
    // prefactor diverges, and the ratio of constants used to come back at
    // 5.97e6, taking the miss column with it through 1 - sum(hits). The
    // reference now drops an empty list before expanding -- exactly, since
    // setting z_l = 0 in the generating function removes list l from E(m) and
    // 0! leaves prod_l m_l! unchanged -- so the row is a distribution again.
    const Matrix<double> p = cache_prob_spm(gamma_6x2<double>(), M21);
    CHECK(p(0, 0) == doctest::Approx(0.25323703432606559).epsilon(MATLAB_TOL));
    CHECK(p(0, 1) == doctest::Approx(0.51746656445373584).epsilon(MATLAB_TOL));
    CHECK(p(0, 2) == doctest::Approx(0.22929640122019859).epsilon(MATLAB_TOL));
    CHECK(p(0, 0) + p(0, 1) + p(0, 2) == doctest::Approx(1.0).epsilon(1e-14));

    // Against the exact cache_prob_erec probabilities the approximation is
    // within 4% on list 1 and a tenth on the zero-capacity list 2, which is the
    // accuracy of the saddle point on a cache this small and no longer a
    // divergence.
    const Matrix<double> e = cache_prob_erec(gamma_6x2<double>(), M21);
    CHECK(std::fabs(p(0, 1) / e(0, 1) - 1.0) < 0.04);
    CHECK(std::fabs(p(0, 2) / e(0, 2) - 1.0) < 0.10);
}

TEST_CASE("cache_miss, cache_miss_fpi and cache_miss_spm agree with MATLAB") {
    const Matrix<double> g = gamma_6x2<double>();
    const Matrix<double> lam = lambda_1x6<double>();

    const auto exact = cache_miss(g, M21, lam);
    CHECK(exact.M == doctest::Approx(0.57794392523364502).epsilon(MATLAB_TOL));
    CHECK(exact.pi0[0] == doctest::Approx(0.25233644859813087).epsilon(MATLAB_TOL));
    CHECK(exact.MU[0] == doctest::Approx(0.57794392523364491).epsilon(MATLAB_TOL));

    const auto fpi = cache_miss_fpi(g, M21, lam);
    CHECK(fpi.M == doctest::Approx(0.60663626194348108).epsilon(MATLAB_TOL));
    CHECK(fpi.pi0[0] == doctest::Approx(0.28796440700390641).epsilon(MATLAB_TOL));

    const auto spm = cache_miss_spm(g, M21, lam);
    CHECK(spm.M == doctest::Approx(0.57544119106685176).epsilon(MATLAB_TOL));
    CHECK(spm.pi0[0] == doctest::Approx(0.25386083144429716).epsilon(MATLAB_TOL));
    CHECK(spm.lE == doctest::Approx(-0.095344776376504514).epsilon(MATLAB_TOL));

    // The saddle point tracks the exact rate to 0.5%, the decoupling
    // approximation to 5%, on this six-item cache.
    CHECK(std::fabs(spm.M / exact.M - 1.0) < 0.005);
    CHECK(std::fabs(fpi.M / exact.M - 1.0) < 0.06);

    // MI must be the request rate times the miss probability, identically.
    for (std::size_t k = 0; k < 6; ++k) {
        CHECK(exact.MI[k] == doctest::Approx(lam(0, k) * exact.pi0[k]).epsilon(1e-15));
        CHECK(fpi.MI[k] == doctest::Approx(lam(0, k) * fpi.pi0[k]).epsilon(1e-15));
    }
}

TEST_CASE("cache_t_hlru reduces to the Che closed form on a symmetric cache") {
    // Four identical items of rate 1/4 and two slots: the fixed point is
    // 1 - exp(-lam T) = m/n = 1/2, so T = log(2)/lam = 4 log 2.
    Matrix<double> g(4, 1, 0.25);
    const std::vector<double> t = cache_t_hlru(g, std::vector<int>{2});
    CHECK(t[0] == doctest::Approx(4.0 * std::log(2.0)).epsilon(1e-10));

    const Matrix<double> p = cache_ttl_hlru(Matrix<double>(1, 4, 0.25), std::vector<int>{2});
    for (std::size_t k = 0; k < 4; ++k) CHECK(p(k, 1) == doctest::Approx(0.5).epsilon(1e-10));
}

TEST_CASE("cache_ttl_hlru matches MATLAB and conserves probability and capacity") {
    Matrix<double> lam(1, 4);
    lam(0, 0) = 0.4;
    lam(0, 1) = 0.3;
    lam(0, 2) = 0.2;
    lam(0, 3) = 0.1;
    Matrix<double> g(4, 1);
    for (std::size_t k = 0; k < 4; ++k) g(k, 0) = lam(0, k);

    SUBCASE("single list is the Che approximation") {
        const std::vector<double> t = cache_t_hlru(g, std::vector<int>{2});
        CHECK(t[0] == doctest::Approx(2.9938911919827547).epsilon(MATLAB_TOL));
        const Matrix<double> p = cache_ttl_hlru(lam, std::vector<int>{2});
        CHECK(p(0, 0) == doctest::Approx(0.30193108687768755).epsilon(MATLAB_TOL));
        CHECK(p(0, 1) == doctest::Approx(0.69806891312231245).epsilon(MATLAB_TOL));
        CHECK(p(3, 1) == doctest::Approx(0.25872908943399797).epsilon(MATLAB_TOL));
        // pi_1(k) = 1 - exp(-lam_k T), the closed form of the h = 1 fixed point
        double occ = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            CHECK(p(k, 1) == doctest::Approx(1.0 - std::exp(-lam(0, k) * t[0])).epsilon(1e-12));
            CHECK(p(k, 0) + p(k, 1) == doctest::Approx(1.0).epsilon(1e-14));
            occ += p(k, 1);
        }
        CHECK(occ == doctest::Approx(2.0).epsilon(1e-9));
    }

    SUBCASE("two lists") {
        const std::vector<int> m{2, 1};
        const std::vector<double> t = cache_t_hlru(g, m);
        CHECK(t[0] == doctest::Approx(5.5191305663272701).epsilon(MATLAB_TOL));
        CHECK(t[1] == doctest::Approx(1.5157121232396675).epsilon(MATLAB_TOL));

        const Matrix<double> p = cache_ttl_hlru(lam, m);
        CHECK(p(0, 0) == doctest::Approx(0.063124000380167841).epsilon(MATLAB_TOL));
        CHECK(p(0, 1) == doctest::Approx(0.51094711236129609).epsilon(MATLAB_TOL));
        CHECK(p(0, 2) == doctest::Approx(0.42592888725853612).epsilon(MATLAB_TOL));

        double occ0 = 0.0, occ1 = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            CHECK(p(k, 0) + p(k, 1) + p(k, 2) == doctest::Approx(1.0).epsilon(1e-14));
            occ0 += p(k, 1);
            occ1 += p(k, 2);
        }
        CHECK(occ0 == doctest::Approx(2.0).epsilon(1e-8));
        CHECK(occ1 == doctest::Approx(1.0).epsilon(1e-8));

        // Real50 must find the same fixed point.
        Matrix<Real50> gr(4, 1);
        for (std::size_t k = 0; k < 4; ++k) gr(k, 0) = Real50(lam(0, k));
        const std::vector<Real50> tr = cache_t_hlru(gr, m);
        CHECK(static_cast<double>(tr[0]) == doctest::Approx(t[0]).epsilon(TOL));
        CHECK(static_cast<double>(tr[1]) == doctest::Approx(t[1]).epsilon(TOL));
    }
}
