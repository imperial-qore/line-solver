/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * THE CROSS-ARITHMETIC CONTRACT.
 *
 * Every function the registry declares at more than one arithmetic promises,
 * implicitly, that the arithmetic is an implementation detail: the Double
 * result must be the Real50 result rounded, and where an Exact instantiation
 * is declared it must be that same number computed without rounding at all.
 * This file asserts that promise directly.
 *
 * WHY NOTHING ELSE IN THE SUITE CATCHES THIS. Every other test fixes one
 * arithmetic and compares against MATLAB or a closed form. An algorithm that
 * takes a DIFFERENT PATH at a different precision still passes all of them:
 *   - a convergence loop with an absolute tolerance stops at a different
 *     iterate when the residual is computed at 50 digits;
 *   - a hardcoded epsilon (1e-12, GlobalConstants.Zero) is a meaningless
 *     threshold in a field with no rounding, so an exact instantiation can
 *     take a branch the double one never takes, or vice versa;
 *   - a branch on floating-point equality (x == 1) fires at one precision and
 *     not at another, which is exactly the qsys_mapg1 defect already on
 *     record.
 * The failure is silent: each instantiation is self-consistent and only the
 * COMPARISON between them shows the algorithm is not arithmetic-agnostic.
 *
 * METHOD. Fixtures are taken from the tests that already exercise each
 * function, not invented, because the object under test is the contract and
 * not new coverage. Each case reduces the function's outputs to a vector of
 * doubles and hands the same input to every declared arithmetic. The
 * comparison is relative, 1e-9, against the higher-precision result.
 *
 * READING A FAILURE. A disagreement here is NOT necessarily a bug in the port.
 * It says the algorithm is precision-dependent and that the registry entry
 * promises more than the algorithm delivers. The resolution may be to narrow
 * the registry entry rather than to change the code; both are decisions for
 * the caller, so this file reports numbers and does not paper over them.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"

// mc
#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_multi.h"
#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_sens.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/ctmc_solve_reducible_blkdecomp.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/api/mc/ctmc_uniformization.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/mc/dtmc_solve_reducible.h"
// cache
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_gamma_lp.h"
#include "line/api/cache/cache_miss.h"
#include "line/api/cache/cache_miss_fpi.h"
#include "line/api/cache/cache_miss_spm.h"
#include "line/api/cache/cache_mva.h"
#include "line/api/cache/cache_mva_miss.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_prob_fpi.h"
#include "line/api/cache/cache_prob_spm.h"
#include "line/api/cache/cache_spm.h"
#include "line/api/cache/cache_t_hlru.h"
#include "line/api/cache/cache_ttl_hlru.h"
#include "line/api/cache/cache_xi_fp.h"
#include "line/api/cache/cache_xi_iter.h"
// mam
#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/aph2_adjust.h"
#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/aph2_fitall.h"
#include "line/api/mam/aph_fit.h"
#include "line/api/mam/map2_fit.h"
#include "line/api/mam/map2mmpp.h"
#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmpp2_fit3.h"
#include "line/api/mam/mmpp2_fitc.h"
#include "line/api/mam/mmpp2_fitc_approx.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/mam/qbd_r.h"
// pfqn
#include "line/api/pfqn/pfqn_aql.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_cub.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_linearizer.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_propfair.h"
#include "line/api/pfqn/pfqn_recal.h"
#include "line/api/pfqn/pas_placement.h"
#include "line/api/pfqn/pas_swap2order.h"
// qsys and the small domains
#include "line/api/aoi/aoi_dist2ph.h"
#include "line/api/fj/fj_dist2fj.h"
#include "line/api/da/da_fpi.h"
#include "line/api/lossn/lossn_erlangfp.h"
#include "line/api/lossn/lossn_mci.h"
#include "line/api/lossn/lossn_manjunath.h"
#include "line/api/npfqn/npfqn_traffic_merge.h"
#include "line/api/npfqn/npfqn_traffic_split_cs.h"
#include "line/api/qsys/qsys_bmapphnn_retrial.h"
#include "line/api/dqsys/dqsys_geogeo1.h"
#include "line/api/qsys/qsys_gg1.h"
#include "line/api/qsys/qsys_gig1_approx_klb.h"
#include "line/api/qsys/qsys_mapg1.h"
#include "line/api/qsys/qsys_mapm1.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/api/qsys/qsys_mmk.h"

using line::Matrix;
using line::Rational;
using line::Real50;

namespace {

// ---------------------------------------------------------------------------
// The harness
// ---------------------------------------------------------------------------

/** Outputs of one instantiation, flattened. */
using Out = std::vector<double>;

template <class T>
double as_double(const T& v) {
    return line::num_traits<T>::to_double(v);
}

/** Append helpers, so a case body reads as a list of what it compares. */
template <class T>
void put(Out& o, const T& v) {
    o.push_back(as_double(v));
}
template <class T>
void put(Out& o, const std::vector<T>& v) {
    for (const T& x : v) o.push_back(as_double(x));
}
template <class T>
void put(Out& o, const Matrix<T>& m) {
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) o.push_back(as_double(m(i, j)));
}

/**
 * Compare two instantiations elementwise, relative to the higher-precision
 * one. The tolerance is the contract's 1e-9; it is NOT widened per function.
 */
void agree(const char* name, const char* lhs, const char* rhs, const Out& a, const Out& b,
           double tol = 1e-9) {
    const std::string tag = std::string(name) + ": " + lhs + " vs " + rhs;
    INFO(tag);
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double scale = std::abs(b[i]) > 1.0 ? std::abs(b[i]) : 1.0;
        INFO("output ", i, ": ", std::string(lhs), " = ", a[i], ", ", std::string(rhs), " = ",
             b[i], ", |d| = ", std::abs(a[i] - b[i]));
        CHECK(std::abs(a[i] - b[i]) <= tol * scale);
    }
}

/**
 * Run one generic lambda at double and Real50 and compare. F must be callable
 * as f(T{}) and return an Out.
 */
template <class F>
void cross_dr(const char* name, F f) {
    const Out d = f(double(0));
    const Out r = f(Real50(0));
    agree(name, "double", "real50", d, r);
}

/** Run at double, Real50 and Rational, and compare all three. */
template <class F>
void cross_dre(const char* name, F f) {
    const Out d = f(double(0));
    const Out r = f(Real50(0));
    const Out e = f(Rational(0));
    agree(name, "double", "real50", d, r);
    agree(name, "double", "exact", d, e);
    agree(name, "real50", "exact", r, e);
}

// ---------------------------------------------------------------------------
// Fixtures, copied from the tests that already exercise these functions
// ---------------------------------------------------------------------------

/** test_ctmc_solve.cpp: birth-death chain with K+1 states. */
template <class T>
Matrix<T> birth_death(int K, long ln, long ld, long mn, long md) {
    Matrix<T> Q(K + 1, K + 1, line::num_traits<T>::from_int(0));
    for (int i = 0; i < K; ++i) {
        Q(i, i + 1) = line::num_traits<T>::from_rational(ln, ld);
        Q(i + 1, i) = line::num_traits<T>::from_rational(mn, md);
    }
    return Q;
}

/** test_mc_aggregation.cpp: nearly completely decomposable chain. */
template <class T>
Matrix<T> ncd_chain_A() {
    using nt = line::num_traits<T>;
    Matrix<T> Q(6, 6, nt::from_int(0));
    Q(0, 1) = nt::from_int(3);
    Q(1, 0) = nt::from_int(2);
    Q(2, 3) = nt::from_int(5);
    Q(3, 2) = nt::from_int(1);
    Q(4, 5) = nt::from_int(2);
    Q(5, 4) = nt::from_int(4);
    Q(1, 2) = nt::from_rational(1, 100);
    Q(2, 1) = nt::from_rational(1, 100);
    Q(3, 4) = nt::from_rational(1, 100);
    Q(4, 3) = nt::from_rational(1, 100);
    Q(5, 0) = nt::from_rational(1, 100);
    Q(0, 5) = nt::from_rational(1, 100);
    return line::mc::ctmc_makeinfgen(Q);
}

/** test_cache_asymptotic.cpp: six items, two lists. */
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

/** test_cache_asymptotic.cpp: the matching arrival rates. */
template <class T>
Matrix<T> lambda_1x6() {
    static const int num[6] = {50, 40, 30, 20, 10, 5};
    Matrix<T> l(1, 6);
    for (int i = 0; i < 6; ++i) l(0, i) = line::num_traits<T>::from_rational(num[i], 100);
    return l;
}

}  // namespace

// ===========================================================================
// mc
// ===========================================================================

TEST_CASE("cross-arithmetic: mc, the exact-capable stationary solvers") {
    cross_dre("ctmc_makeinfgen", [](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_makeinfgen(birth_death<T>(6, 3, 2, 5, 2)));
        return o;
    });
    cross_dre("ctmc_solve", [](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_solve(line::mc::ctmc_makeinfgen(birth_death<T>(6, 3, 2, 5, 2))));
        return o;
    });
    cross_dre("ctmc_stochcomp", [](auto proto) {
        using T = decltype(proto);
        Out o;
        const line::mc::StochCompResult<T> r =
            line::mc::ctmc_stochcomp(ncd_chain_A<T>(), std::vector<std::size_t>{0, 1, 2});
        put(o, r.S);
        return o;
    });
    cross_dre("dtmc_solve", [](auto proto) {
        using T = decltype(proto);
        Out o;
        const Matrix<T> P = line::mc::dtmc_makestochastic(
            line::mc::ctmc_randomization(ncd_chain_A<T>()).P);
        put(o, line::mc::dtmc_solve(P));
        return o;
    });
    cross_dre("dtmc_makestochastic", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> P(3, 3, nt::from_int(0));
        P(0, 0) = nt::from_rational(1, 2);
        P(0, 1) = nt::from_rational(1, 4);
        P(1, 2) = nt::from_int(2);
        P(2, 0) = nt::from_rational(3, 5);
        P(2, 2) = nt::from_rational(1, 5);
        Out o;
        put(o, line::mc::dtmc_makestochastic(P));
        return o;
    });
    cross_dre("ctmc_randomization", [](auto proto) {
        using T = decltype(proto);
        const line::mc::RandomizationResult<T> r = line::mc::ctmc_randomization(ncd_chain_A<T>());
        Out o;
        put(o, r.q);
        put(o, r.P);
        return o;
    });
    cross_dre("ctmc_sens", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const Matrix<T> Q = line::mc::ctmc_makeinfgen(birth_death<T>(4, 3, 2, 5, 2));
        Matrix<T> dQ(5, 5, nt::from_int(0));
        for (std::size_t i = 0; i + 1 < 5; ++i) {
            dQ(i, i + 1) = nt::from_int(1);
            dQ(i, i) = nt::from_int(-1);
        }
        Out o;
        put(o, line::mc::ctmc_sens(Q, dQ));
        return o;
    });
    cross_dre("ctmc_solve_reducible", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        // two absorbing states reachable from a transient pair
        Matrix<T> Q(4, 4, nt::from_int(0));
        Q(0, 1) = nt::from_int(1);
        Q(1, 0) = nt::from_int(1);
        Q(0, 2) = nt::from_int(2);
        Q(1, 3) = nt::from_int(3);
        Q = line::mc::ctmc_makeinfgen(Q);
        std::vector<T> pi0(4, nt::from_int(0));
        pi0[0] = nt::from_int(1);
        Out o;
        put(o, line::mc::ctmc_solve_reducible(Q, pi0).pi);
        return o;
    });
    cross_dre("dtmc_solve_reducible", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> P(4, 4, nt::from_int(0));
        P(0, 1) = nt::from_rational(1, 2);
        P(0, 2) = nt::from_rational(1, 2);
        P(1, 0) = nt::from_rational(1, 2);
        P(1, 3) = nt::from_rational(1, 2);
        P(2, 2) = nt::from_int(1);
        P(3, 3) = nt::from_int(1);
        std::vector<T> pi0(4, nt::from_int(0));
        pi0[0] = nt::from_int(1);
        Out o;
        put(o, line::mc::dtmc_solve_reducible(P, pi0).pi);
        return o;
    });
    cross_dre("ctmc_solve_reducible_blkdecomp", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> Q(4, 4, nt::from_int(0));
        Q(0, 1) = nt::from_int(1);
        Q(1, 0) = nt::from_int(1);
        Q(2, 3) = nt::from_int(2);
        Q(3, 2) = nt::from_int(4);
        Q = line::mc::ctmc_makeinfgen(Q);
        std::vector<T> pi0(4, nt::from_rational(1, 4));
        Out o;
        put(o, line::mc::ctmc_solve_reducible_blkdecomp(Q, pi0).pi);
        return o;
    });
}

TEST_CASE("cross-arithmetic: mc, the iterative and transient solvers") {
    using Blocks = std::vector<std::vector<std::size_t>>;
    const Blocks MS{{0, 1}, {2, 3}, {4, 5}};
    cross_dr("ctmc_courtois", [&MS](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_courtois(ncd_chain_A<T>(), MS).p);
        return o;
    });
    cross_dr("ctmc_kms", [&MS](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_kms(ncd_chain_A<T>(), MS, 5).p);
        return o;
    });
    cross_dr("ctmc_takahashi", [&MS](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_takahashi(ncd_chain_A<T>(), MS, 5).p);
        return o;
    });
    const Blocks MSS{{0, 1}, {2}};
    cross_dr("ctmc_multi", [&MS, &MSS](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::mc::ctmc_multi(ncd_chain_A<T>(), MS, MSS).p);
        return o;
    });
    cross_dr("ctmc_uniformization", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const Matrix<T> Q = line::mc::ctmc_makeinfgen(birth_death<T>(4, 3, 2, 5, 2));
        std::vector<T> pi0(5, nt::from_int(0));
        pi0[0] = nt::from_int(1);
        Out o;
        put(o, line::mc::ctmc_uniformization(pi0, Q, nt::from_int(2)).pi);
        return o;
    });
    cross_dr("ctmc_transient", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const Matrix<T> Q = line::mc::ctmc_makeinfgen(birth_death<T>(3, 3, 2, 5, 2));
        std::vector<T> pi0(4, nt::from_int(0));
        pi0[0] = nt::from_int(1);
        Out o;
        const line::mc::TransientResult<T> r =
            line::mc::ctmc_transient(Q, pi0, nt::from_int(0), nt::from_int(1), 1e-8, 1e-10);
        for (std::size_t j = 0; j < r.pi.cols(); ++j) put(o, r.pi(r.pi.rows() - 1, j));
        return o;
    });
}

// ===========================================================================
// cache
// ===========================================================================

TEST_CASE("cross-arithmetic: cache, the exact-capable recursions") {
    const std::vector<int> M21{2, 1};
    cross_dre("cache_erec", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_erec(gamma_6x2<T>(), M21));
        return o;
    });
    cross_dre("cache_prob_erec", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_prob_erec(gamma_6x2<T>(), M21));
        return o;
    });
    cross_dre("cache_miss", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        const line::cache::CacheMissResult<T> r =
            line::cache::cache_miss(gamma_6x2<T>(), M21, lambda_1x6<T>());
        put(o, r.M);
        put(o, r.MU);
        put(o, r.MI);
        put(o, r.pi0);
        return o;
    });
    cross_dre("cache_mva", [&M21](auto proto) {
        using T = decltype(proto);
        const line::cache::CacheMvaResult<T> r = line::cache::cache_mva(gamma_6x2<T>(), M21);
        Out o;
        put(o, r.pi);
        put(o, r.pi0);
        return o;
    });
    cross_dre("cache_mva_miss", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> p{nt::from_rational(2, 5), nt::from_rational(3, 10),
                         nt::from_rational(1, 5), nt::from_rational(1, 10)};
        Matrix<T> R(2, 4, nt::from_rational(9, 10));
        const line::cache::CacheMvaMissResult<T> r =
            line::cache::cache_mva_miss(p, std::vector<int>{2, 1}, R);
        Out o;
        put(o, r.M);
        put(o, r.Mk);
        return o;
    });
    cross_dre("cache_gamma_lp", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<Matrix<T>> lambda(1, Matrix<T>(1, 3));
        lambda[0](0, 0) = nt::from_rational(1, 2);
        lambda[0](0, 1) = nt::from_rational(1, 3);
        lambda[0](0, 2) = nt::from_rational(1, 5);
        std::vector<std::vector<Matrix<T>>> R(1, std::vector<Matrix<T>>(1));
        R[0][0] = Matrix<T>(3, 3, nt::from_int(0));
        R[0][0](0, 1) = nt::from_int(1);
        R[0][0](1, 2) = nt::from_int(1);
        Out o;
        put(o, line::cache::cache_gamma_lp(lambda, R).gamma);
        return o;
    });
}

TEST_CASE("cross-arithmetic: cache, the fixed points and asymptotics") {
    const std::vector<int> M21{2, 1};
    cross_dr("cache_xi_fp", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_xi_fp(gamma_6x2<T>(), M21).xi);
        return o;
    });
    cross_dr("cache_xi_iter", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_xi_iter(gamma_6x2<T>(), M21));
        return o;
    });
    cross_dr("cache_prob_fpi", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_prob_fpi(gamma_6x2<T>(), M21));
        return o;
    });
    cross_dr("cache_miss_fpi", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        const line::cache::CacheMissResult<T> r =
            line::cache::cache_miss_fpi(gamma_6x2<T>(), M21, lambda_1x6<T>());
        put(o, r.M);
        put(o, r.MU);
        return o;
    });
    cross_dr("cache_spm", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_spm(gamma_6x2<T>(), M21).lZ);
        return o;
    });
    cross_dr("cache_prob_spm", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_prob_spm(gamma_6x2<T>(), M21));
        return o;
    });
    cross_dr("cache_miss_spm", [&M21](auto proto) {
        using T = decltype(proto);
        Out o;
        const line::cache::CacheMissSpmResult<T> r =
            line::cache::cache_miss_spm(gamma_6x2<T>(), M21, lambda_1x6<T>());
        put(o, r.M);
        put(o, r.MU);
        return o;
    });
    cross_dr("cache_ttl_hlru", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> lam(1, 4, nt::from_rational(1, 4));
        Out o;
        put(o, line::cache::cache_ttl_hlru(lam, std::vector<int>{2}));
        return o;
    });
    cross_dr("cache_t_hlru", [](auto proto) {
        using T = decltype(proto);
        Out o;
        put(o, line::cache::cache_t_hlru(gamma_6x2<T>(), std::vector<int>{2}));
        return o;
    });
}

// ===========================================================================
// mam
// ===========================================================================

namespace {

/** test_mam_fit_optim.cpp: the MMPP(2) every optimizer-based fit is fed. */
template <class T>
line::mam::Map<T> reference_mmpp2() {
    using nt = line::num_traits<T>;
    line::mam::Map<T> m;
    m.D0 = Matrix<T>(2, 2, nt::from_int(0));
    m.D1 = Matrix<T>(2, 2, nt::from_int(0));
    const T l1 = nt::from_int(2), l2 = nt::from_rational(1, 2);
    const T r1 = nt::from_rational(3, 10), r2 = nt::from_rational(7, 10);
    m.D0(0, 0) = -T(l1 + r1);
    m.D0(0, 1) = r1;
    m.D0(1, 0) = r2;
    m.D0(1, 1) = -T(l2 + r2);
    m.D1(0, 0) = l1;
    m.D1(1, 1) = l2;
    return m;
}

template <class T>
T idc_at(const line::mam::Map<T>& m, const T& t) {
    const std::vector<T> tv(1, t);
    const std::vector<T> v = line::mam::map_count_var(m, tv);
    const std::vector<T> mu = line::mam::map_count_mean(m, tv);
    return T(v[0] / mu[0]);
}

template <class T>
T m3_counts_at(const line::mam::Map<T>& m, const T& t) {
    std::vector<unsigned> ord{1u, 2u, 3u};
    const std::vector<T> mt = line::mam::map_count_moment(m, t, ord);
    return T(mt[2] - line::num_traits<T>::from_int(3) * mt[1] * mt[0] +
             line::num_traits<T>::from_int(2) * mt[0] * mt[0] * mt[0]);
}

/** test_qbd_family.cpp: the M/M/1 QBD blocks, rho = lambda/mu. */
template <class T>
void mm1_qbd(Matrix<T>& B, Matrix<T>& L, Matrix<T>& F, long ln, long ld, long mn, long md) {
    using nt = line::num_traits<T>;
    const T lam = nt::from_rational(ln, ld), mu = nt::from_rational(mn, md);
    B = Matrix<T>(1, 1, mu);
    L = Matrix<T>(1, 1, T(-T(lam + mu)));
    F = Matrix<T>(1, 1, lam);
}

}  // namespace

TEST_CASE("cross-arithmetic: mam, the closed-form moment fits") {
    cross_dr("aph_fit", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::AphFitResult<T> r =
            line::mam::aph_fit(nt::from_int(1), nt::from_int(3), nt::from_int(20));
        Out o;
        put(o, static_cast<double>(r.order));
        put(o, r.aph.D0);
        put(o, r.aph.D1);
        return o;
    });
    cross_dr("aph2_fitall", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const std::vector<line::mam::Map<T>> v =
            line::mam::aph2_fitall(nt::from_int(1), nt::from_int(3), nt::from_int(20));
        Out o;
        put(o, static_cast<double>(v.size()));
        for (const line::mam::Map<T>& m : v) {
            put(o, m.D0);
            put(o, m.D1);
        }
        return o;
    });
    cross_dr("aph2_fit", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Aph2FitResult<T> r =
            line::mam::aph2_fit(nt::from_int(1), nt::from_int(3), nt::from_int(20));
        Out o;
        put(o, r.aph.D0);
        put(o, r.aph.D1);
        put(o, static_cast<double>(r.adjusted));
        return o;
    });
    cross_dr("aph2_adjust", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        // an infeasible triple, so the adjustment really runs
        const line::mam::Aph2AdjustResult<T> r = line::mam::aph2_adjust(
            nt::from_int(1), nt::from_rational(3, 2), nt::from_int(2));
        Out o;
        put(o, r.M2a);
        put(o, r.M3a);
        return o;
    });
    cross_dr("map2_fit", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Map2FitResult<T> r = line::mam::map2_fit(
            nt::from_int(1), nt::from_int(4), nt::from_int(30), nt::from_rational(3, 10));
        Out o;
        put(o, static_cast<double>(r.has_map));
        put(o, static_cast<double>(r.err));
        put(o, r.map.D0);
        put(o, r.map.D1);
        put(o, r.e3_used);
        return o;
    });
    cross_dr("mmpp2_fit3", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Map<T> m = line::mam::mmpp2_fit3(
            nt::from_int(1), nt::from_int(4), nt::from_int(30), nt::from_rational(1, 2));
        Out o;
        put(o, m.D0);
        put(o, m.D1);
        return o;
    });
    cross_dr("amap2_fit_gamma", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Amap2FitGammaResult<T> r = line::mam::amap2_fit_gamma(
            nt::from_int(1), nt::from_int(4), nt::from_int(30), nt::from_rational(3, 10));
        Out o;
        put(o, static_cast<double>(r.poisson_fallback));
        put(o, r.amap.D0);
        put(o, r.amap.D1);
        return o;
    });
    cross_dre("map2mmpp", [](auto proto) {
        using T = decltype(proto);
        const line::mam::Map2mmppResult<T> r = line::mam::map2mmpp(reference_mmpp2<T>());
        Out o;
        put(o, r.Q);
        put(o, r.LAMBDA);
        put(o, r.offdiag_norm);
        return o;
    });
    cross_dre("map_mark", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> prob(2, 2);
        prob(0, 0) = nt::from_rational(3, 10);
        prob(0, 1) = nt::from_rational(7, 10);
        prob(1, 0) = nt::from_rational(1, 2);
        prob(1, 1) = nt::from_rational(1, 2);
        const line::mam::Mmap<T> r = line::mam::mmap_mark(reference_mmpp2<T>(), prob);
        Out o;
        put(o, r.D0);
        put(o, r.D1);
        for (const Matrix<T>& d : r.Dc) put(o, d);
        return o;
    });
}

TEST_CASE("cross-arithmetic: mam, the counting statistics and the optimizer fits") {
    cross_dr("map_count_mean/var/moment, idc", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Map<T> m = reference_mmpp2<T>();
        Out o;
        put(o, idc_at(m, nt::from_int(1)));
        put(o, idc_at(m, nt::from_int(5)));
        put(o, m3_counts_at(m, nt::from_int(5)));
        return o;
    });
    cross_dr("mmpp2_fitc", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Map<T> src = reference_mmpp2<T>();
        const T t1 = nt::from_int(1), t2 = nt::from_int(5);
        const line::mam::Mmpp2FitcResult<T> r = line::mam::mmpp2_fitc(
            line::mam::map_lambda(src), idc_at(src, t1), idc_at(src, t2),
            idc_at(src, nt::from_double(1e6)), m3_counts_at(src, t2), t1, t2);
        Out o;
        put(o, static_cast<double>(r.degenerate));
        put(o, static_cast<double>(r.third_moment_ok));
        put(o, r.map.D0);
        put(o, r.map.D1);
        return o;
    });
    cross_dr("mmpp2_fitc_approx", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::mam::Map<T> src = reference_mmpp2<T>();
        const T t1 = nt::from_int(1), t2 = nt::from_int(5);
        const line::mam::Mmpp2FitcApproxResult<T> r = line::mam::mmpp2_fitc_approx(
            line::mam::map_lambda(src), idc_at(src, t1), idc_at(src, t2),
            idc_at(src, nt::from_double(1e6)), m3_counts_at(src, t2), t1, t2);
        Out o;
        put(o, r.map.D0);
        put(o, r.map.D1);
        put(o, r.objective);
        return o;
    });
}

TEST_CASE("cross-arithmetic: mam, the QBD family") {
    cross_dr("qbd_R", [](auto proto) {
        using T = decltype(proto);
        Matrix<T> B, L, F;
        mm1_qbd(B, L, F, 1, 2, 1, 1);
        Out o;
        put(o, line::mam::qbd_R(B, L, F));
        return o;
    });
    cross_dr("qbd_fundmat", [](auto proto) {
        using T = decltype(proto);
        Matrix<T> B, L, F;
        mm1_qbd(B, L, F, 1, 2, 1, 1);
        const line::mam::QbdFundMat<T> r = line::mam::qbd_fundmat(B, L, F);
        Out o;
        put(o, r.R);
        put(o, r.G);
        return o;
    });
    cross_dr("qbd_caudal", [](auto proto) {
        using T = decltype(proto);
        Matrix<T> B, L, F;
        mm1_qbd(B, L, F, 1, 2, 1, 1);
        Out o;
        put(o, line::mam::qbd_caudal(line::mam::qbd_R(B, L, F)));
        return o;
    });
    cross_dre("qbd_pi", [](auto proto) {
        using T = decltype(proto);
        Matrix<T> B, L, F;
        mm1_qbd(B, L, F, 1, 2, 1, 1);
        // For the M/M/1 QBD the minimal solution is exactly rho, so R can be
        // written down and qbd_pi is reached without qbd_R's cyclic reduction
        // (which is double/real only, and correctly declared so).
        const Matrix<T> R(1, 1, line::num_traits<T>::from_rational(1, 2));
        Matrix<T> Lbar = L;
        for (std::size_t i = 0; i < Lbar.rows(); ++i) Lbar(i, i) = Lbar(i, i) + B(i, i);
        Out o;
        put(o, line::mam::qbd_pi(B, Lbar, R, static_cast<std::size_t>(12),
                                 line::num_traits<T>::from_int(0)));
        return o;
    });
    cross_dr("qbd_mapmap1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> a;
        a.D0 = Matrix<T>(1, 1, nt::from_int(-1));
        a.D1 = Matrix<T>(1, 1, nt::from_int(1));
        line::mam::Map<T> s;
        s.D0 = Matrix<T>(1, 1, nt::from_int(-2));
        s.D1 = Matrix<T>(1, 1, nt::from_int(2));
        const line::mam::QbdMapMap1Result<T> r = line::mam::qbd_mapmap1(a, s);
        Out o;
        put(o, r.QN);
        put(o, r.UN);
        put(o, r.RN);
        put(o, r.XN);
        return o;
    });
}

// ===========================================================================
// pfqn
// ===========================================================================

namespace {

/**
 * test_pfqn_linearizer.cpp / test_pfqn_amva2.cpp: two stations, two classes,
 * demands 0.5/0.3 and 0.4/0.6. Written as rationals so the exact
 * instantiation gets the SAME model and not a rounded one.
 */
template <class T>
Matrix<T> pfqn_demands() {
    using nt = line::num_traits<T>;
    Matrix<T> L(2, 2);
    L(0, 0) = nt::from_rational(1, 2);
    L(0, 1) = nt::from_rational(3, 10);
    L(1, 0) = nt::from_rational(2, 5);
    L(1, 1) = nt::from_rational(3, 5);
    return L;
}

/** Think times of the same fixture. */
template <class T>
Matrix<T> pfqn_think() {
    using nt = line::num_traits<T>;
    Matrix<T> Z(1, 2);
    Z(0, 0) = nt::from_int(1);
    Z(0, 1) = nt::from_rational(1, 2);
    return Z;
}

}  // namespace

TEST_CASE("cross-arithmetic: pfqn, the exact normalizing constants") {
    const std::vector<int> N{4, 3};
    cross_dre("pfqn_ca", [&N](auto proto) {
        using T = decltype(proto);
        const line::pfqn::NcResult<T> r = line::pfqn::pfqn_ca(pfqn_demands<T>(), N, pfqn_think<T>());
        Out o;
        put(o, r.G);
        put(o, r.lG);
        return o;
    });
    cross_dre("pfqn_recal", [&N](auto proto) {
        using T = decltype(proto);
        const line::pfqn::NcResult<T> r =
            line::pfqn::pfqn_recal(pfqn_demands<T>(), N, pfqn_think<T>());
        Out o;
        put(o, r.G);
        return o;
    });
    cross_dre("pfqn_mva", [&N](auto proto) {
        using T = decltype(proto);
        const line::pfqn::MvaResult<T> r =
            line::pfqn::pfqn_mva(pfqn_demands<T>(), N, pfqn_think<T>());
        Out o;
        put(o, r.XN);
        put(o, r.QN);
        put(o, r.UN);
        put(o, r.CN);
        return o;
    });
    cross_dre("pfqn_comomrm", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> L(1, 2);
        L(0, 0) = nt::from_rational(1, 2);
        L(0, 1) = nt::from_rational(3, 10);
        Matrix<T> Z(1, 2);
        Z(0, 0) = nt::from_int(1);
        Z(0, 1) = nt::from_rational(1, 2);
        const line::pfqn::ComomResult<T> r =
            line::pfqn::pfqn_comomrm(L, std::vector<int>{3, 2}, Z);
        Out o;
        put(o, r.G);
        put(o, r.lG);
        return o;
    });
    cross_dre("pfqn_gld", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> L(2, 1);
        L(0, 0) = nt::from_rational(1, 2);
        L(1, 0) = nt::from_rational(2, 5);
        Matrix<T> mu(2, 4);
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t n = 0; n < 4; ++n)
                mu(i, n) = nt::from_int(static_cast<long>(n) + 1);
        const line::pfqn::NcResult<T> r = line::pfqn::pfqn_gld(L, std::vector<int>{4}, mu);
        Out o;
        put(o, r.G);
        return o;
    });
}

TEST_CASE("cross-arithmetic: pfqn, the approximate mean value analyses") {
    const std::vector<int> N{4, 3};
    cross_dr("pfqn_linearizer", [&N](auto proto) {
        using T = decltype(proto);
        const line::pfqn::LinearizerResult<T> r =
            line::pfqn::pfqn_linearizer(pfqn_demands<T>(), N, pfqn_think<T>());
        Out o;
        put(o, r.X);
        put(o, r.Q);
        put(o, r.U);
        put(o, r.C);
        return o;
    });
    cross_dr("pfqn_bs", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> N{nt::from_int(4), nt::from_int(3)};
        std::vector<T> Z{nt::from_int(1), nt::from_rational(1, 2)};
        const line::pfqn::AmvaResult<T> r = line::pfqn::pfqn_bs(pfqn_demands<T>(), N, Z);
        Out o;
        put(o, r.XN);
        put(o, r.QN);
        put(o, r.UN);
        return o;
    });
    cross_dr("pfqn_aql", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> N{nt::from_int(4), nt::from_int(3)};
        std::vector<T> Z{nt::from_int(1), nt::from_rational(1, 2)};
        const line::pfqn::AmvaResult<T> r = line::pfqn::pfqn_aql(pfqn_demands<T>(), N, Z);
        Out o;
        put(o, r.XN);
        put(o, r.QN);
        return o;
    });
    cross_dr("pfqn_propfair", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> N{nt::from_int(4), nt::from_int(3)};
        std::vector<T> Z{nt::from_int(1), nt::from_rational(1, 2)};
        const line::pfqn::PropfairResult<T> r =
            line::pfqn::pfqn_propfair(pfqn_demands<T>(), N, Z);
        Out o;
        put(o, r.G);
        put(o, r.lG);
        put(o, r.Xasy);
        return o;
    });
    cross_dr("pfqn_cub", [&N](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> Z{nt::from_int(0), nt::from_int(0)};
        const line::pfqn::CubResult<T> r = line::pfqn::pfqn_cub(pfqn_demands<T>(), N, Z);
        Out o;
        put(o, r.lG);
        return o;
    });
}

// ===========================================================================
// qsys
// ===========================================================================

TEST_CASE("cross-arithmetic: qsys, the closed-form queues") {
    cross_dre("qsys_mm1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::QsysResult<T> r =
            line::qsys::qsys_mm1(nt::from_int(2), nt::from_int(3));
        Out o;
        put(o, r.W);
        put(o, r.rhohat);
        return o;
    });
    cross_dre("qsys_mmk", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::QsysResult<T> r =
            line::qsys::qsys_mmk(nt::from_int(5), nt::from_int(2), 3u);
        Out o;
        put(o, r.W);
        put(o, r.rhohat);
        return o;
    });
    cross_dre("qsys_mg1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::QsysResult<T> r =
            line::qsys::qsys_mg1(nt::from_int(2), nt::from_int(3), nt::from_rational(1, 2));
        Out o;
        put(o, r.W);
        put(o, r.rhohat);
        return o;
    });
    cross_dre("qsys_mmck", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::MmckResult<T> r =
            line::qsys::qsys_mmck(nt::from_int(5), nt::from_int(2), 2u, 6u);
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.meanWaitingTime);
        put(o, r.utilization);
        put(o, r.throughput);
        return o;
    });
    cross_dr("qsys_gig1_approx_klb", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::QsysResult<T> r = line::qsys::qsys_gig1_approx_klb(
            nt::from_int(2), nt::from_int(3), nt::from_rational(3, 2), nt::from_rational(1, 2));
        Out o;
        put(o, r.W);
        put(o, r.rhohat);
        return o;
    });
    cross_dr("qsys_gg1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::qsys::QsysResult<T> r = line::qsys::qsys_gg1(
            nt::from_int(2), nt::from_int(3), nt::from_rational(3, 2), nt::from_rational(1, 2));
        Out o;
        put(o, r.W);
        put(o, r.rhohat);
        return o;
    });
    cross_dre("dqsys_geogeo1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        const line::dqsys::GeoGeo1Result<T> r =
            line::dqsys::dqsys_geogeo1(nt::from_rational(1, 4), nt::from_rational(1, 2));
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.utilization);
        put(o, r.emptyProb);
        put(o, r.ratio);
        return o;
    });
    cross_dr("qsys_mapm1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> a;
        a.D0 = Matrix<T>(1, 1, nt::from_int(-2));
        a.D1 = Matrix<T>(1, 1, nt::from_int(2));
        const line::qsys::MapMcResult<T> r = line::qsys::qsys_mapm1(a, nt::from_int(3));
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.meanWaitingTime);
        put(o, r.utilization);
        return o;
    });
    cross_dr("qsys_mapph1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> a;
        a.D0 = Matrix<T>(1, 1, nt::from_int(-2));
        a.D1 = Matrix<T>(1, 1, nt::from_int(2));
        std::vector<T> sigma{nt::from_int(1), nt::from_int(0)};
        Matrix<T> S(2, 2, nt::from_int(0));
        S(0, 0) = nt::from_int(-6);
        S(0, 1) = nt::from_int(6);
        S(1, 1) = nt::from_int(-6);
        const line::qsys::MapMap1Result<T> r = line::qsys::qsys_mapph1(a, sigma, S, 20);
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.meanWaitingTime);
        put(o, r.utilization);
        return o;
    });
    cross_dr("qsys_mapg1", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> a;
        a.D0 = Matrix<T>(1, 1, nt::from_int(-2));
        a.D1 = Matrix<T>(1, 1, nt::from_int(2));
        // interior moment set, so the fit is not on the APH(2) boundary
        std::vector<T> mom{nt::from_rational(1, 3), nt::from_rational(8, 45),
                           nt::from_rational(2, 15)};
        const line::qsys::MapG1Result<T> r = line::qsys::qsys_mapg1(a, mom, 20);
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.meanWaitingTime);
        put(o, r.utilization);
        // servicePhases IS compared: it used to diverge (2 at double, 3 at
        // Real50) and no longer does, see the regression case at the end.
        put(o, static_cast<double>(r.servicePhases));
        return o;
    });
    cross_dre("qsys_bmapphnn_retrial", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<Matrix<T>> D;
        D.push_back(Matrix<T>(1, 1, nt::from_int(-1)));
        D.push_back(Matrix<T>(1, 1, nt::from_int(1)));
        std::vector<T> beta(1, nt::from_int(1));
        Matrix<T> S(1, 1, nt::from_int(-1));
        line::qsys::BmapPhNnRetrialOptions opt;
        opt.maxLevel = 6;
        const line::qsys::BmapPhNnRetrialResult<T> r = line::qsys::qsys_bmapphnn_retrial(
            D, beta, S, 2, nt::from_rational(1, 2), nt::from_rational(1, 5),
            nt::from_rational(1, 4), std::vector<long>(1, 1L), opt);
        Out o;
        put(o, r.L_orbit);
        put(o, r.N_server);
        put(o, r.throughput);
        put(o, r.P_idle);
        return o;
    });
}

// ===========================================================================
// the small domains: npfqn, lossn, fj, aoi, pfqn/pas
// ===========================================================================

TEST_CASE("cross-arithmetic: npfqn, lossn, fj, aoi and the pass-and-swap order") {
    // cross_dr, NOT cross_dre: `npfqn_traffic_merge` cannot be instantiated at
    // Rational at all. Its Merge::Interpos branch calls
    // `m3pp2m_fitc_theoretical`, whose `static_assert(has_transcendental)` fires
    // when the template is instantiated, not when the branch is taken -- so even
    // the Merge::Super config used below fails to compile under exact
    // arithmetic. The counting-characteristic fit needs logarithms; the claim was
    // the defect, not the fitter.
    cross_dr("npfqn_traffic_merge", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Mmap<T> A;
        A.D0 = Matrix<T>{{nt::from_int(-3), nt::from_int(1)},
                         {nt::from_int(1), nt::from_int(-4)}};
        Matrix<T> a1{{nt::from_int(1), nt::from_rational(1, 2)},
                     {nt::from_int(1), nt::from_int(1)}};
        Matrix<T> a2{{nt::from_rational(1, 2), nt::from_int(0)},
                     {nt::from_rational(1, 2), nt::from_rational(1, 2)}};
        A.D1 = Matrix<T>(2, 2, nt::from_int(0));
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) A.D1(i, j) = a1(i, j) + a2(i, j);
        A.Dc.push_back(a1);
        A.Dc.push_back(a2);
        const line::mam::Mmap<T> S = line::npfqn::npfqn_traffic_merge<T>(
            {A, A}, line::npfqn::MergeConfig{line::npfqn::Merge::Super,
                                             line::npfqn::Compress::None});
        Out o;
        put(o, S.D0);
        put(o, line::mam::mmap_lambda(S));
        return o;
    });
    cross_dre("npfqn_traffic_split_cs", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::npfqn::Mmap<T> M(4, Matrix<T>(2, 2, nt::from_int(0)));
        M[2](0, 0) = nt::from_rational(1, 2);
        M[2](0, 1) = nt::from_rational(1, 4);
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
        Matrix<T> P(2, 4, nt::from_rational(1, 4));
        const std::vector<line::npfqn::Mmap<T>> out = line::npfqn::npfqn_traffic_split_cs(M, P);
        Out o;
        for (const line::npfqn::Mmap<T>& s : out)
            for (const Matrix<T>& b : s) put(o, b);
        return o;
    });
    cross_dr("lossn_erlangfp", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> nu{nt::from_int(3), nt::from_int(2)};
        Matrix<T> A(2, 2, nt::from_int(0));
        A(0, 0) = nt::from_int(1);
        A(0, 1) = nt::from_int(1);
        A(1, 1) = nt::from_int(1);
        const std::vector<int> C{6, 4};
        const line::lossn::ErlangFpResult<T> r =
            line::lossn::lossn_erlangfp(nu, A, C, line::da::FpiOptions());
        Out o;
        put(o, r.Loss);
        put(o, r.QLen);
        return o;
    });
    // The one loss-network analyzer that reaches exact arithmetic: the series is
    // rational and the metrics are ratios of series values, so unlike its two
    // siblings it has no tolerance and no sampling error to lose precision to.
    cross_dre("lossn_manjunath", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        std::vector<T> nu{nt::from_int(3), nt::from_rational(3, 2)};
        Matrix<T> A(2, 2, nt::from_int(0));
        A(0, 0) = nt::from_int(1);
        A(0, 1) = nt::from_int(1);
        A(1, 0) = nt::from_int(1);
        A(1, 1) = nt::from_int(2);
        const std::vector<T> C{nt::from_int(4), nt::from_int(5)};
        const line::lossn::LossnManjunathResult<T> r = line::lossn::lossn_manjunath<T>(nu, A, C);
        Out o;
        put(o, r.Loss);
        put(o, r.QLen);
        put(o, r.lG);
        return o;
    });
    cross_dr("lossn_mci", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::lossn::LossnMciOptions<T> opt;
        opt.samples = 2000;
        const line::lossn::LossnMciResult<T> r =
            line::lossn::lossn_mci<T>({nt::from_int(4)}, Matrix<T>(1, 1, nt::from_int(1)),
                                      {nt::from_int(5)}, opt, static_cast<std::uint64_t>(7));
        Out o;
        put(o, r.Loss);
        put(o, r.lG);
        return o;
    });
    cross_dre("fj_dist2fj", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> m;
        m.D0 = Matrix<T>{{nt::from_int(-2), nt::from_int(0)},
                         {nt::from_int(0), nt::from_int(-5)}};
        m.D1 = Matrix<T>{{nt::from_rational(6, 5), nt::from_rational(4, 5)},
                         {nt::from_int(2), nt::from_int(3)}};
        const line::fj::FjDist<T> a =
            line::fj::fj_dist2fj(m, line::fj::FjDistKind::Service, line::fj::FjProcType::Erlang);
        Out o;
        put(o, a.mu);
        put(o, a.St);
        put(o, a.tau_st);
        return o;
    });
    cross_dre("aoi_dist2ph", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> m;
        m.D0 = Matrix<T>{{nt::from_int(-3), nt::from_int(1)},
                         {nt::from_int(1), nt::from_int(-4)}};
        m.D1 = Matrix<T>{{nt::from_rational(3, 2), nt::from_rational(1, 2)},
                         {nt::from_rational(3, 2), nt::from_rational(3, 2)}};
        const line::aoi::AoiPh<T> ph = line::aoi::aoi_dist2ph(m);
        Out o;
        put(o, ph.alpha);
        put(o, ph.Tmat);
        return o;
    });
    cross_dre("pas_placement", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> H(3, 3, nt::from_int(0));
        H(0, 1) = nt::from_int(1);
        H(1, 2) = nt::from_int(1);
        const line::pfqn::PasPlacement<T> p = line::pfqn::pas_placement(H);
        Out o;
        put(o, p.P);
        const std::vector<T> x(3, nt::from_int(1));
        for (std::size_t i : p.placeable(x)) put(o, static_cast<double>(i));
        return o;
    });
    cross_dre("pas_swap2order", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        Matrix<T> G(3, 3, nt::from_int(0));
        G(0, 1) = nt::from_int(1);
        G(1, 0) = nt::from_int(1);
        const std::vector<T> mu{nt::from_int(1), nt::from_int(2), nt::from_int(3)};
        line::pfqn::PasRateFun<T> rate = [mu](const std::vector<int>& c) -> T {
            if (c.empty()) return line::num_traits<T>::from_int(0);
            return mu[static_cast<std::size_t>(c[0]) - 1];
        };
        Out o;
        put(o, line::pfqn::pas_swap2order<T>({G, G}, {rate, rate}, {1, 1, 1}));
        return o;
    });
}

// ===========================================================================
// Findings: divergences the sweep uncovered, pinned with their numbers
// ===========================================================================

TEST_CASE("regression: aph_fit picks the same order at every precision") {
    // HISTORY. This was a FINDING of the cross-arithmetic sweep before it was a
    // regression test. On (e1, e2, e3) = (1/3, 8/45, 2/15) -- cv2 = 0.6, the
    // interior set the MATLAB comparison in test_qsys_mapg1.cpp uses -- the
    // fitted ORDER used to depend on the working precision, and not
    // monotonically: 2 at double, 3 at Real50, 2 at Real100.
    //
    // The cause was NOT a hardcoded tolerance in the order search; the search
    // had no tolerance at all. n3 = 9/4 and the APH(2) upper bound un(2) = 9/4
    // are EQUAL for this moment set, and the exact comparison n3 <= un was
    // decided by the sign of a rounding residual: +8.88e-16 at double,
    // -1.07e-50 at Real50, exactly 0 at Real100. A frozen constant cannot fix
    // that; a precision-derived one can, and fitdetail::aph_boundary_slack now
    // supplies it (16 eps of the working type, zero in an exact field).
    //
    // What this test pins: every arithmetic returns MATLAB's answer, order 2,
    // and the fit still reproduces the moments it was given.
    using line::mam::aph_fit;
    using nt_d = line::num_traits<double>;
    const line::mam::AphFitResult<double> d =
        aph_fit(nt_d::from_rational(1, 3), nt_d::from_rational(8, 45), nt_d::from_rational(2, 15));
    using nt_r = line::num_traits<Real50>;
    const line::mam::AphFitResult<Real50> r =
        aph_fit(nt_r::from_rational(1, 3), nt_r::from_rational(8, 45), nt_r::from_rational(2, 15));
    using R100 = line::Real<100>;
    using nt_h = line::num_traits<R100>;
    const line::mam::AphFitResult<R100> h =
        aph_fit(nt_h::from_rational(1, 3), nt_h::from_rational(8, 45), nt_h::from_rational(2, 15));
    CHECK(d.order == 2);  // MATLAB aph_fit(1/3, 8/45, 2/15) -> order 2
    CHECK(r.order == 2);
    CHECK(h.order == 2);
    // the rates agree across the three, which they could not while the orders
    // differed
    CHECK(static_cast<double>(r.aph.D0(0, 0)) ==
          doctest::Approx(d.aph.D0(0, 0)).epsilon(1e-12));
    CHECK(static_cast<double>(h.aph.D0(1, 1)) ==
          doctest::Approx(d.aph.D0(1, 1)).epsilon(1e-12));
    // and each still carries the moments it was asked for
    CHECK(line::mam::map_moment(d.aph, 1) == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(line::mam::map_moment(d.aph, 2) == doctest::Approx(8.0 / 45.0).epsilon(1e-12));
    CHECK(line::mam::map_moment(d.aph, 3) == doctest::Approx(2.0 / 15.0).epsilon(1e-12));
    CHECK(static_cast<double>(line::mam::map_moment(r.aph, 3)) ==
          doctest::Approx(2.0 / 15.0).epsilon(1e-12));

    // The Erlang(2) moment set (1/3, 1/6, 1/9) is the other attained boundary,
    // and it bites one level below the comparisons: n2 = 3/2 sits ON
    // n2 = (n+1)/n, where the radicand of the APH(2) upper bound un is exactly
    // zero. It came out as 0 at double, -1.07e-50 at Real50 and 0 at Real100,
    // and the placeholder that guards the genuinely infeasible region turned
    // that noise into un = 0, so the n3 <= un test failed and Real50 returned
    // order 3 for a distribution that IS an APH(2).
    //
    // fitdetail::read_nonneg_radicand now reads a provably nonnegative radicand
    // that rounds negative as zero, at the three sites where one vanishes on an
    // attained boundary. MATLAB aph_fit(1/3, 1/6, 1/9) returns order 2 with
    // isexact true, and every arithmetic now agrees with it.
    const line::mam::AphFitResult<double> ed =
        aph_fit(nt_d::from_rational(1, 3), nt_d::from_rational(1, 6), nt_d::from_rational(1, 9));
    const line::mam::AphFitResult<Real50> er =
        aph_fit(nt_r::from_rational(1, 3), nt_r::from_rational(1, 6), nt_r::from_rational(1, 9));
    const line::mam::AphFitResult<R100> eh =
        aph_fit(nt_h::from_rational(1, 3), nt_h::from_rational(1, 6), nt_h::from_rational(1, 9));
    CHECK(ed.order == 2);
    CHECK(er.order == 2);
    CHECK(eh.order == 2);
    CHECK(ed.isexact);
    CHECK(er.isexact);
    CHECK(eh.isexact);
    // an order alone is not a fit: each still carries the Erlang(2) moments
    CHECK(line::mam::map_moment(ed.aph, 1) == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(line::mam::map_moment(ed.aph, 2) == doctest::Approx(1.0 / 6.0).epsilon(1e-12));
    CHECK(line::mam::map_moment(ed.aph, 3) == doctest::Approx(1.0 / 9.0).epsilon(1e-12));
    CHECK(static_cast<double>(line::mam::map_moment(er.aph, 2)) ==
          doctest::Approx(1.0 / 6.0).epsilon(1e-12));
    CHECK(static_cast<double>(line::mam::map_moment(er.aph, 3)) ==
          doctest::Approx(1.0 / 9.0).epsilon(1e-12));
    CHECK(static_cast<double>(line::mam::map_moment(eh.aph, 3)) ==
          doctest::Approx(1.0 / 9.0).epsilon(1e-12));
}

TEST_CASE("cross-arithmetic: qsys_mapg1 under a MAP arrival, where the fitted order matters") {
    // THE HAZARD THIS PINS. With Poisson arrivals the M/G/1 mean wait is
    // Pollaczek-Khinchine and depends on the first two service moments alone,
    // so two different three-moment fits give the SAME answer and a divergent
    // fitted order is invisible in the metrics. Under a MAP arrival the queue
    // depends on the whole service law, so an order that differed between
    // arithmetics would move meanQueueLength and this case would fail.
    //
    // The moment set is the one that used to flip (1/3, 8/45, 2/15), and the
    // arrival is the correlated MMPP(2) of the MATLAB reference runs, so the
    // test exercises exactly the combination that would have been wrong.
    cross_dr("qsys_mapg1 (MAP arrival, boundary moments)", [](auto proto) {
        using T = decltype(proto);
        using nt = line::num_traits<T>;
        line::mam::Map<T> a;
        a.D0 = Matrix<T>{{nt::from_double(-2.5), nt::from_double(0.2)},
                         {nt::from_double(0.1), nt::from_double(-0.7)}};
        a.D1 = Matrix<T>{{nt::from_double(2.3), nt::from_int(0)},
                         {nt::from_int(0), nt::from_double(0.6)}};
        const std::vector<T> mom{nt::from_rational(1, 3), nt::from_rational(8, 45),
                                 nt::from_rational(2, 15)};
        const line::qsys::MapG1Result<T> r = line::qsys::qsys_mapg1(a, mom, 40);
        Out o;
        put(o, r.meanQueueLength);
        put(o, r.meanWaitingTime);
        put(o, r.meanSojournTime);
        put(o, r.utilization);
        put(o, static_cast<double>(r.servicePhases));
        for (std::size_t i = 0; i < 8 && i < r.queueLengthDist.size(); ++i)
            put(o, r.queueLengthDist[i]);
        return o;
    });
    // and the double instantiation is still the one MATLAB agrees with:
    // qsys_mapg1(A0, A1, [1/3 8/45 2/15]) -> meanQueueLength 0.858359337155326
    using nt = line::num_traits<double>;
    line::mam::Map<double> a;
    a.D0 = Matrix<double>{{-2.5, 0.2}, {0.1, -0.7}};
    a.D1 = Matrix<double>{{2.3, 0.0}, {0.0, 0.6}};
    const std::vector<double> mom{nt::from_rational(1, 3), nt::from_rational(8, 45),
                                  nt::from_rational(2, 15)};
    const line::qsys::MapG1Result<double> r = line::qsys::qsys_mapg1(a, mom, 40);
    CHECK(r.meanQueueLength == doctest::Approx(0.858359337155326).epsilon(1e-10));
    CHECK(r.meanWaitingTime == doctest::Approx(0.402403241371232).epsilon(1e-10));
    CHECK(r.servicePhases == 2);
}
