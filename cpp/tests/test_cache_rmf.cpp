/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The refined mean-field cache miss rates (cache_miss_rmf.h) and the RANDOM(m)
 * multi-list mean-field steady state (cache_rrm_meanfield.h), both of which the
 * stiff integrator in line/util/ode.h unblocked.
 *
 * The oracle is MATLAB R2025a; the exact commands that produced each number are
 * quoted at the assertion. Beyond the oracle the tests check what the reference
 * does not: that the fixed point really is one (the drift vanishes there), that
 * the occupancies stay on the simplex, that the reference's rmf_jacobian
 * departs from the true derivative of its rmf_drift in exactly the way the
 * header documents, and that the null-space assumption the dimension reduction
 * rests on holds.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_miss_rmf.h"
#include "line/api/cache/cache_rrm_meanfield.h"
#include "line/util/eig.h"

using line::Matrix;
using line::Real50;
using line::cache::cache_miss_rmf;
using line::cache::cache_miss_rmf_index;
using line::cache::CacheMissRmfResult;

namespace {

Matrix<double> two_user_four_item() {
    Matrix<double> lam(2, 4, 0.0);
    const double a[2][4] = {{0.5, 0.3, 0.15, 0.05}, {0.1, 0.2, 0.3, 0.4}};
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 4; ++j) lam(i, j) = a[i][j];
    return lam;
}

Matrix<double> one_user_seven_item() {
    Matrix<double> lam(1, 7, 0.0);
    const double b[7] = {49, 49, 49, 49, 7, 1, 1};
    for (std::size_t j = 0; j < 7; ++j) lam(0, j) = b[j] / 205.0;
    return lam;
}

}  // namespace

TEST_CASE("cache_miss_rmf: against MATLAB") {
    // MATLAB, for each case:
    //   [M,MU,MI,pi0] = cache_miss_rmf(ones(1,n), m, lambda)
    // The tolerance is 1e-6 relative, which is loose against the 1e-8 the two
    // codes actually agree to; the slack is there because both sides locate the
    // fixed point by integration and the reference's is the less converged of
    // the two (see the header).
    SUBCASE("two users, four items, one list of size 2") {
        const CacheMissRmfResult<double> r =
            cache_miss_rmf(std::vector<double>(4, 1.0), std::vector<int>{2}, two_user_four_item());
        CHECK(r.M == doctest::Approx(0.992377021789716).epsilon(1e-6));
        CHECK(r.MU[0] == doctest::Approx(0.479544182002491).epsilon(1e-6));
        CHECK(r.MU[1] == doctest::Approx(0.512832839787225).epsilon(1e-6));
        CHECK(r.MI[0] == doctest::Approx(0.269891715469127).epsilon(1e-6));
        CHECK(r.MI[1] == doctest::Approx(0.249040929224349).epsilon(1e-6));
        CHECK(r.pi0[0] == doctest::Approx(0.449819525781878).epsilon(1e-6));
        CHECK(r.pi0[1] == doctest::Approx(0.498081858448699).epsilon(1e-6));
        CHECK(r.pi0[2] == doctest::Approx(0.526049307884710).epsilon(1e-6));
        CHECK(r.pi0[3] == doctest::Approx(0.526049307884710).epsilon(1e-6));
        CHECK(r.refined);
    }

    SUBCASE("two users, four items, two lists") {
        const CacheMissRmfResult<double> r = cache_miss_rmf(
            std::vector<double>(4, 1.0), std::vector<int>{1, 2}, two_user_four_item());
        CHECK(r.M == doctest::Approx(0.490863272756674).epsilon(1e-6));
        CHECK(r.MU[0] == doctest::Approx(0.225357496042656).epsilon(1e-6));
        CHECK(r.MU[1] == doctest::Approx(0.265505776714018).epsilon(1e-6));
        CHECK(r.pi0[0] == doctest::Approx(0.190478156602372).epsilon(1e-6));
        CHECK(r.pi0[1] == doctest::Approx(0.245831177173324).epsilon(1e-6));
        CHECK(r.pi0[2] == doctest::Approx(0.281845323440741).epsilon(1e-6));
        CHECK(r.refined);
    }

    SUBCASE("one user, seven items, three lists") {
        const CacheMissRmfResult<double> r = cache_miss_rmf(
            std::vector<double>(7, 1.0), std::vector<int>{1, 1, 3}, one_user_seven_item());
        CHECK(r.M == doctest::Approx(0.0240177201681542).epsilon(1e-6));
        CHECK(r.pi0[0] == doctest::Approx(0.00581319251768532).epsilon(1e-6));
        CHECK(r.pi0[4] == doctest::Approx(0.301249947485579).epsilon(1e-6));
        CHECK(r.pi0[5] == doctest::Approx(0.837748630899249).epsilon(1e-6));
        CHECK(r.refined);
    }

    SUBCASE("one user, seven items, one list of size 3: the reference falls back") {
        // Here cache_miss_rmf.m's 1/N refinement produces non-finite entries and
        // the reference silently keeps the plain mean-field fixed point,
        // M = 0.348362022568998. The port's rank decision keeps the reduced
        // Jacobian non-singular, so it applies the correction and returns
        // 0.345660745093153, 0.8 percent lower. Both numbers are asserted: the
        // port's own, and that the port's UNREFINED value is the reference's,
        // which is what makes the difference attributable to the refinement
        // rather than to the mean field.
        const CacheMissRmfResult<double> r = cache_miss_rmf(
            std::vector<double>(7, 1.0), std::vector<int>{3}, one_user_seven_item());
        CHECK(r.refined);
        CHECK(r.M == doctest::Approx(0.345660745093153).epsilon(1e-6));

        // The plain mean field, recomputed here without the refinement.
        const Matrix<double> lam = one_user_seven_item();
        const std::size_t n = 7, h = 1;
        std::vector<double> p(n, 0.0), m{3.0};
        double tot = 0.0;
        for (std::size_t j = 0; j < n; ++j) tot += lam(0, j);
        for (std::size_t j = 0; j < n; ++j) p[j] = lam(0, j) / tot;
        std::vector<double> x0(n * (h + 1), 0.0);
        x0[cache_miss_rmf_index(0, 1, n)] = 1.0;
        x0[cache_miss_rmf_index(1, 1, n)] = 1.0;
        x0[cache_miss_rmf_index(2, 1, n)] = 1.0;
        for (std::size_t i = 3; i < n; ++i) x0[cache_miss_rmf_index(i, 0, n)] = 1.0;
        const std::vector<double> xss =
            line::cache::rmf_detail::fixed_point(x0, p, m, n, h, 10000.0);
        double M_mf = 0.0;
        for (std::size_t i = 0; i < n; ++i) M_mf += lam(0, i) * xss[cache_miss_rmf_index(i, 0, n)];
        CHECK(M_mf == doctest::Approx(0.348362022568998).epsilon(1e-6));
    }
}

TEST_CASE("cache_miss_rmf: structural properties the reference does not check") {
    const Matrix<double> lam = two_user_four_item();
    const std::size_t n = 4, h = 2;
    const CacheMissRmfResult<double> r =
        cache_miss_rmf(std::vector<double>(n, 1.0), std::vector<int>{1, 2}, lam);

    // Each item's occupancies are a distribution over the h+1 lists.
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k <= h; ++k) s += r.xss[cache_miss_rmf_index(i, k, n)];
        CHECK(s == doctest::Approx(1.0).epsilon(1e-9));
    }
    // Miss probabilities are probabilities, and the global miss rate is the sum
    // of the per-item ones and also the sum of the per-user ones.
    double sum_MI = 0.0, sum_MU = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(r.pi0[i] >= 0.0);
        CHECK(r.pi0[i] <= 1.0);
        sum_MI += r.MI[i];
    }
    for (double v : r.MU) sum_MU += v;
    CHECK(sum_MI == doctest::Approx(r.M).epsilon(1e-12));
    CHECK(sum_MU == doctest::Approx(r.M).epsilon(1e-12));

    // Items with equal aggregate request rates must come out with equal miss
    // probabilities. Items 3 and 4 both have rate 0.45 here, and the asymmetry
    // is a sensitive probe of the refinement: an ill-conditioned reduction
    // breaks it long before it breaks anything else.
    CHECK(std::fabs(r.pi0[2] - r.pi0[3]) < 1e-10);
}

TEST_CASE("cache_miss_rmf: the drift, its Jacobian and its null space") {
    const std::size_t n = 4, h = 2;
    std::vector<double> p{0.6, 0.5, 0.45, 0.45}, m{1.0, 2.0};
    double tot = 0.0;
    for (double v : p) tot += v;
    for (double& v : p) v /= tot;

    std::vector<double> x0(n * (h + 1), 0.0);
    x0[cache_miss_rmf_index(0, 1, n)] = 1.0;
    x0[cache_miss_rmf_index(1, 2, n)] = 1.0;
    x0[cache_miss_rmf_index(2, 2, n)] = 1.0;
    x0[cache_miss_rmf_index(3, 0, n)] = 1.0;
    std::vector<double> xss = line::cache::rmf_detail::fixed_point(x0, p, m, n, h, 10000.0);
    xss = line::cache::rmf_detail::fixed_point(xss, p, m, n, h, 10000.0, 1e-13, 1e-16);

    // It is a fixed point: the drift vanishes there.
    const std::vector<double> F = line::cache::rmf_detail::drift(xss, p, m, n, h);
    for (double v : F) CHECK(std::fabs(v) < 1e-13);

    // The drift conserves each item's total occupancy, everywhere, not just at
    // the fixed point.
    std::vector<double> xr(n * (h + 1), 0.0);
    for (std::size_t i = 0; i < n * (h + 1); ++i) xr[i] = 0.1 + 0.03 * static_cast<double>(i % 7);
    const std::vector<double> Fr = line::cache::rmf_detail::drift(xr, p, m, n, h);
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k <= h; ++k) s += Fr[cache_miss_rmf_index(i, k, n)];
        CHECK(std::fabs(s) < 1e-14);
    }

    // The reference's rmf_jacobian is NOT the derivative of its rmf_drift, and
    // the port reproduces the reference (see the defect note on jacobian()).
    // This pins the discrepancy exactly rather than papering over it: build the
    // true derivative by central differences, and check that jacobian() equals
    // it EXCEPT in the entries d F(i,k) / d x(j,k+1), for EVERY j, where the
    // reference adds -p(i) x(i,k) / m(k+1) (and +p(i) x(i,k)/m(k+1) in the
    // paired row F(i,k+1)). If the transcription were wrong in any other way,
    // the difference would not have exactly this form.
    const std::size_t md = n * (h + 1);
    const Matrix<double> J = line::cache::rmf_detail::jacobian(xr, p, m, n, h);
    Matrix<double> D(md, md, 0.0);  // the true derivative
    const double d = 1e-6;
    for (std::size_t j = 0; j < md; ++j) {
        std::vector<double> xp = xr, xm = xr;
        xp[j] += d;
        xm[j] -= d;
        const std::vector<double> fp = line::cache::rmf_detail::drift(xp, p, m, n, h);
        const std::vector<double> fm = line::cache::rmf_detail::drift(xm, p, m, n, h);
        for (std::size_t i = 0; i < md; ++i) D(i, j) = (fp[i] - fm[i]) / (2 * d);
    }
    Matrix<double> extra(md, md, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k + 1 <= h; ++k) {
            const std::size_t ik = cache_miss_rmf_index(i, k, n);
            const std::size_t ik1 = cache_miss_rmf_index(i, k + 1, n);
            for (std::size_t j = 0; j < n; ++j) {
                const std::size_t jk1 = cache_miss_rmf_index(j, k + 1, n);
                extra(ik, jk1) -= p[i] * xr[ik] / m[k];
                extra(ik1, jk1) += p[i] * xr[ik] / m[k];
            }
        }
    bool any_extra = false;
    for (std::size_t i = 0; i < md; ++i)
        for (std::size_t j = 0; j < md; ++j) {
            CHECK(J(i, j) == doctest::Approx(D(i, j) + extra(i, j)).epsilon(1e-6));
            if (std::fabs(extra(i, j)) > 1e-9) any_extra = true;
        }
    CHECK(any_extra);  // the defect is present on this input, not vacuously absent

    // The null space the dimension reduction rests on. The elimination path and
    // the SVD path must agree on the RANK, which is the only thing about them
    // that affects the answer.
    const Matrix<double> Fp = line::cache::rmf_detail::jacobian(xss, p, m, n, h);
    const line::cache::rmf_detail::NullSpace<double> a =
        line::cache::rmf_detail::left_null_space(Fp);
    const line::cache::rmf_detail::NullSpace<double> b =
        line::cache::rmf_detail::left_null_space_svd(Fp);
    CHECK(a.rank == b.rank);
    CHECK(a.rank < n * (h + 1));
    CHECK(a.basis.size() == n * (h + 1) - a.rank);
    // Every returned vector really is a left null vector.
    for (const std::vector<double>& w : a.basis)
        for (std::size_t j = 0; j < n * (h + 1); ++j) {
            double s = 0.0;
            for (std::size_t i = 0; i < n * (h + 1); ++i) s += w[i] * Fp(i, j);
            CHECK(std::fabs(s) < 1e-10);
        }
    // ... and the item-conservation vectors are in the span, which is what makes
    // the rank at most model_dim - n_items.
    CHECK(a.rank <= n * (h + 1) - n);
}

TEST_CASE("cache_miss_rmf: transient trajectory") {
    const Matrix<double> lam = two_user_four_item();
    const std::size_t n = 4, h = 2;
    std::vector<double> x0(n * (h + 1), 0.0);
    for (std::size_t i = 0; i < n; ++i) x0[cache_miss_rmf_index(i, 0, n)] = 1.0;
    const CacheMissRmfResult<double> tr =
        line::cache::cache_miss_rmf_transient(std::vector<int>{1, 2}, lam, 0.0, 50.0, x0);
    REQUIRE(tr.tout.size() > 2);
    CHECK(tr.tout.front() == 0.0);
    CHECK(tr.tout.back() == doctest::Approx(50.0));
    CHECK(tr.pi0_t.rows() == n);
    CHECK(tr.pi0_t.cols() == tr.tout.size());
    // Everything starts outside the cache, so the initial miss probability is 1
    // for every item and decreasing thereafter for the popular ones.
    for (std::size_t i = 0; i < n; ++i) CHECK(tr.pi0_t(i, 0) == doctest::Approx(1.0));
    CHECK(tr.pi0_t(0, tr.tout.size() - 1) < 0.5);
    // The trajectory stays on the simplex throughout.
    for (std::size_t j = 0; j < tr.tout.size(); ++j)
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t k = 0; k <= h; ++k) s += tr.xtraj(cache_miss_rmf_index(i, k, n), j);
            CHECK(s == doctest::Approx(1.0).epsilon(1e-8));
        }
    // The long-run limit of the transient is the PLAIN mean-field fixed point.
    // Not the steady-state result of cache_miss_rmf, which carries the 1/N
    // refinement on top of it: the two differ here by about 1 percent, and
    // confusing them would make the refinement look like an integration error.
    const CacheMissRmfResult<double> longrun =
        line::cache::cache_miss_rmf_transient(std::vector<int>{1, 2}, lam, 0.0, 10000.0, x0);
    std::vector<double> p(n, 0.0), m{1.0, 2.0};
    double tot = 0.0;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < lam.rows(); ++i) tot += lam(i, j);
    for (std::size_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < lam.rows(); ++i) s += lam(i, j);
        p[j] = s / tot;
    }
    const std::vector<double> mf = line::cache::rmf_detail::fixed_point(x0, p, m, n, h, 10000.0);
    for (std::size_t i = 0; i < n; ++i)
        CHECK(longrun.pi0_t(i, longrun.tout.size() - 1) ==
              doctest::Approx(mf[cache_miss_rmf_index(i, 0, n)]).epsilon(1e-6));
}

TEST_CASE("cache_miss_rmf: input validation") {
    const Matrix<double> lam = two_user_four_item();
    CHECK_THROWS_AS(cache_miss_rmf(std::vector<double>(4, 1.0), std::vector<int>{}, lam),
                    line::InputError);
    CHECK_THROWS_AS(cache_miss_rmf(std::vector<double>(4, 1.0), std::vector<int>{0}, lam),
                    line::InputError);
    CHECK_THROWS_AS(
        cache_miss_rmf(std::vector<double>(0), std::vector<int>{1}, Matrix<double>(2, 4, 0.0)),
        line::InputError);
}

TEST_CASE("cache_rrm_meanfield: against MATLAB") {
    // MATLAB, the body of cache_rrm_meanfield.m with its own data:
    //   [t,x] = ode23s(@(t,x) cache_rrm_meanfield_ode(t,x,lambda,m,n,h), [0 10000], x0(:));
    //   xe = reshape(x(end,:),[n,1+h]); lambda*xe(:,1)
    SUBCASE("the reference's own seven-item three-list case") {
        std::vector<double> lambda(7, 0.0);
        const double b[7] = {49, 49, 49, 49, 7, 1, 1};
        for (std::size_t i = 0; i < 7; ++i) lambda[i] = b[i] / 205.0;
        const line::cache::CacheRrmMeanfieldResult<double> r =
            line::cache::cache_rrm_meanfield(lambda, std::vector<int>{1, 1, 3});
        CHECK(r.missrate == doctest::Approx(0.0254753216627099).epsilon(1e-7));
        CHECK(r.missratio == doctest::Approx(0.0254753216627099).epsilon(1e-7));
        CHECK(r.x[0] == doctest::Approx(0.00707147033548331).epsilon(1e-6));
        CHECK(r.x[4] == doctest::Approx(0.310786439499164).epsilon(1e-6));
        CHECK(r.x[5] == doctest::Approx(0.830463839303598).epsilon(1e-6));
        // The four equally popular items must stay equal.
        for (std::size_t i = 1; i < 4; ++i) CHECK(std::fabs(r.x[i] - r.x[0]) < 1e-9);
    }

    SUBCASE("five items, a single list of size 2") {
        const std::vector<double> lambda{0.4, 0.3, 0.15, 0.1, 0.05};
        const line::cache::CacheRrmMeanfieldResult<double> r =
            line::cache::cache_rrm_meanfield(lambda, std::vector<int>{2});
        CHECK(r.missrate == doctest::Approx(0.494488965618057).epsilon(1e-7));
        CHECK(r.missratio == doctest::Approx(0.494488965618057).epsilon(1e-7));
        CHECK(r.x[0] == doctest::Approx(0.381995504598189).epsilon(1e-6));
        CHECK(r.x[1] == doctest::Approx(0.451798950150968).epsilon(1e-6));
        CHECK(r.x[4] == doctest::Approx(0.831788299222618).epsilon(1e-6));
        // Miss probability must increase as popularity decreases.
        for (std::size_t i = 0; i + 1 < 5; ++i) CHECK(r.x[i] < r.x[i + 1]);
        // The occupancies of each item sum to one, and the total occupancy of
        // the cache list equals its capacity.
        double occupied = 0.0;
        for (std::size_t i = 0; i < 5; ++i) {
            CHECK(r.x[i] + r.x[i + 5] == doctest::Approx(1.0).epsilon(1e-8));
            occupied += r.x[i + 5];
        }
        CHECK(occupied == doctest::Approx(2.0).epsilon(1e-6));
    }
}

TEST_CASE("cache_rrm_meanfield: Real50 agrees with double") {
    // The whole pipeline instantiated at 50 digits: the integrator, the drift
    // and (in cache_miss_rmf) the elimination null space rather than LAPACK.
    // The answer must be the same to the accuracy the double run can carry.
    const std::vector<Real50> lambda{Real50("0.4"), Real50("0.3"), Real50("0.15"), Real50("0.1"),
                                     Real50("0.05")};
    const line::cache::CacheRrmMeanfieldResult<Real50> r =
        line::cache::cache_rrm_meanfield(lambda, std::vector<int>{2});
    CHECK(static_cast<double>(r.missrate) == doctest::Approx(0.494488965618057).epsilon(1e-7));

    const std::vector<double> lambda_d{0.4, 0.3, 0.15, 0.1, 0.05};
    const line::cache::CacheRrmMeanfieldResult<double> rd =
        line::cache::cache_rrm_meanfield(lambda_d, std::vector<int>{2});
    CHECK(static_cast<double>(r.missrate) == doctest::Approx(rd.missrate).epsilon(1e-7));
}

TEST_CASE("cache_miss_rmf: Real50 agrees with double") {
    Matrix<Real50> lam(1, 4, Real50(0));
    const char* v[4] = {"0.5", "0.3", "0.15", "0.05"};
    for (std::size_t j = 0; j < 4; ++j) lam(0, j) = Real50(v[j]);
    const line::cache::CacheMissRmfResult<Real50> r =
        cache_miss_rmf(std::vector<Real50>(4, Real50(1)), std::vector<int>{2}, lam);

    Matrix<double> lamd(1, 4, 0.0);
    const double vd[4] = {0.5, 0.3, 0.15, 0.05};
    for (std::size_t j = 0; j < 4; ++j) lamd(0, j) = vd[j];
    const CacheMissRmfResult<double> rd =
        cache_miss_rmf(std::vector<double>(4, 1.0), std::vector<int>{2}, lamd);

    CHECK(r.refined == rd.refined);
    CHECK(static_cast<double>(r.M) == doctest::Approx(rd.M).epsilon(1e-6));
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(static_cast<double>(r.pi0[i]) == doctest::Approx(rd.pi0[i]).epsilon(1e-6));
}
