/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_M1PS_H
#define LINE_API_MAM_MAP_M1PS_H

/**
 * Sojourn time distribution in a MAP/M/1 processor-sharing queue.
 *
 * Templated port of matlab/src/api/mam/map_compute_R.m,
 * map_m1ps_h_recursive.m, map_m1ps_sojourn.m and map_m1ps_cdfrespt.m, which
 * implement Theorem 1 of H. Masuyama and T. Takine, "Sojourn time distribution
 * in a MAP/M/1 processor-sharing queue", Operations Research Letters 31(6),
 * 2003, 406-412.
 *
 * The arrival process is the MAP (C, D) -- C carries the hidden transitions,
 * D the arrivals -- and service is exponential of rate mu shared equally, so
 * with n jobs present each is served at rate mu/n. The queue length is a QBD
 * whose rate matrix R is the minimal nonnegative solution of
 *
 *     D + R (C - mu I) + mu R^2 = 0,
 *
 * and the complementary sojourn time distribution is, by uniformization at
 * theta + mu with theta = max_i |C_ii|,
 *
 *     W^c(x) = (1/lambda) sum_n pi_0 R^n D sum_k e^-(theta+mu)x
 *              ((theta+mu) x)^k / k! h_{n,k},
 *
 * with pi_0 = pi (I - R) and the vectors h_{n,k} from the recursion
 *
 *     h_{n,0}   = e
 *     h_{n,k+1} = [ n mu/(n+1) h_{n-1,k} + (theta I + C) h_{n,k}
 *                   + D h_{n+1,k} ] / (theta + mu),   h_{-1,k} = 0.
 *
 * ARITHMETIC.
 *   - map_m1ps_h_recursive is a FINITE recursion in the entries of C and D,
 *     with theta a maximum of absolute diagonal entries, so it instantiates at
 *     every arithmetic including Rational and returns exact fractions.
 *   - map_compute_R is a fixed-point iteration driven to a tolerance and is
 *     gated on num_traits<T>::has_transcendental, for the same reason as
 *     qbd_R (see qbd_r.h).
 *   - map_m1ps_sojourn and map_m1ps_cdfrespt need the Poisson weights
 *     e^-a a^k / k! and are gated as well.
 *
 * THE TWO ENTRY POINTS ARE NOT THE SAME FUNCTION, despite identical
 * signatures and identical documentation in the reference. They differ in
 * three ways, all reproduced here:
 *   1. The queue-length truncation. map_m1ps_sojourn finds the smallest N with
 *      (1/lambda) sum_{n<=N} pi_0 R^n D e > 1 - epsilon, scanning n = 0..1000,
 *      and falls back to N = 100 when the scan never gets there.
 *      map_m1ps_cdfrespt instead estimates N from the spectral radius of R,
 *      N = ceil(log(epsilon (1 - sp)) / log(sp)), clamps it to [10, 10000] and
 *      then truncates AGAIN at run time as soon as ||pi_0 R^n D||_inf drops
 *      below epsilon/100.
 *   2. The stationary vector. map_m1ps_sojourn solves the (M+1) x M
 *      overdetermined system [Q; e'] x = [0; 1] in the least-squares sense,
 *      whereas map_m1ps_cdfrespt replaces the LAST ROW of Q' by e' and solves
 *      the resulting square system. Both give the stationary vector of an
 *      irreducible generator; the port uses the square solve of ctmc_solve for
 *      both, which is the same vector.
 *   3. R itself. map_m1ps_sojourn calls map_compute_R, iterating
 *      R <- -D (C - mu I + mu R)^-1. map_m1ps_cdfrespt has a PRIVATE
 *      compute_R_matrix that iterates the different splitting
 *      R <- (D + mu R^2) (mu I - C)^-1 from the warm start -D (C - mu I)^-1,
 *      with 5000 rather than 1000 iterations, and for M = 1 solves the scalar
 *      quadratic in closed form. The two fixed points coincide -- both are the
 *      minimal nonnegative solution -- so the two are exposed here as
 *      map_compute_R and map_compute_R_quadratic and the tests check they
 *      agree to the residual of the defining equation.
 *
 * REFERENCE DEFECT. The second output of both functions, W_bar_n, is
 * documented as "the conditional complementary distribution for customers
 * finding n customers in the system", but both compute it as
 * sum(sum_k)/M -- the arithmetic MEAN OVER PHASES of the uniformized
 * h-weighted sum, with no reference to the phase distribution at an arrival
 * and no normalization by the probability of finding n customers. It is not a
 * conditional distribution: on the M/M/1-PS instance of the tests (M = 1,
 * lambda = 0.8, mu = 1) the n = 0 curve at x = 0 is 1.0 and DECREASES
 * correctly, but for M > 1 the phases are weighted uniformly rather than by
 * pie, so it is not a probability of anything. The port returns it under the
 * name w_bar_n_unweighted to make the meaning explicit, and computes it
 * identically so a caller comparing against MATLAB sees the same numbers.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_r.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Residual ||D + R (C - mu I) + mu R^2||_inf of the MAP/M/1 rate equation. */
template <class T>
T map_compute_R_residual(const Matrix<T>& C, const Matrix<T>& D, const T& mu, const Matrix<T>& R) {
    const std::size_t m = C.rows();
    Matrix<T> res = D;
    const Matrix<T> CmuI = [&]() {
        Matrix<T> X = C;
        for (std::size_t i = 0; i < m; ++i) X(i, i) -= mu;
        return X;
    }();
    const Matrix<T> t1 = matmul(R, CmuI);
    const Matrix<T> t2 = matmul(R, R);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) res(i, j) += t1(i, j) + mu * t2(i, j);
    T worst = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < m; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < m; ++j) s += num_abs(T(res(i, j)));
        if (s > worst) worst = s;
    }
    return worst;
}

/**
 * Rate matrix R of a MAP/M/1 queue, the minimal nonnegative solution of
 * D + R (C - mu I) + mu R^2 = 0, by the iteration
 * R <- -D (C - mu I + mu R)^-1 (map_compute_R.m).
 *
 * Negative entries produced by rounding are clamped to zero, as in the
 * reference; the clamp is a no-op whenever the iteration has converged.
 */
template <class T>
Matrix<T> map_compute_R(const Matrix<T>& C, const Matrix<T>& D, const T& mu, unsigned max_iter,
                        const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_compute_R requires transcendental arithmetic");
    const std::size_t m = C.rows();
    if (C.cols() != m || D.rows() != m || D.cols() != m)
        throw InputError("map_compute_R: C and D must be square and of equal order");
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> R(m, m, zero);
    for (unsigned it = 0; it < max_iter; ++it) {
        Matrix<T> X(m, m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) X(i, j) = C(i, j) + mu * R(i, j);
        for (std::size_t i = 0; i < m; ++i) X(i, i) -= mu;
        Matrix<T> Rn = matmul(D, inverse(X));
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) Rn(i, j) = -Rn(i, j);
        T diff = zero;
        for (std::size_t i = 0; i < m; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < m; ++j) s += num_abs(T(Rn(i, j) - R(i, j)));
            if (s > diff) diff = s;
        }
        R = Rn;
        if (diff < tol) break;
    }
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            if (R(i, j) < zero) R(i, j) = zero;
    return R;
}

/** map_compute_R with the reference defaults, 1000 iterations and tolerance 1e-10. */
template <class T>
Matrix<T> map_compute_R(const Matrix<T>& C, const Matrix<T>& D, const T& mu) {
    return map_compute_R(C, D, mu, 1000u, T(num_traits<T>::from_double(1e-10)));
}

/**
 * The same R by the other splitting, R <- (D + mu R^2) (mu I - C)^-1, warm
 * started at -D (C - mu I)^-1 and with the scalar case solved in closed form
 * (the private compute_R_matrix of map_m1ps_cdfrespt.m).
 *
 * For M = 1 the equation is mu R^2 + (C - mu) R + D = 0 and the root in [0, 1)
 * is selected, which for Poisson arrivals (C = -lambda, D = lambda) is exactly
 * the utilization rho = lambda/mu. That closed form is the sharpest available
 * oracle for the matrix iteration.
 */
template <class T>
Matrix<T> map_compute_R_quadratic(const Matrix<T>& C, const Matrix<T>& D, const T& mu,
                                  unsigned max_iter, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_compute_R_quadratic requires transcendental arithmetic");
    const std::size_t m = C.rows();
    if (C.cols() != m || D.rows() != m || D.cols() != m)
        throw InputError("map_compute_R_quadratic: C and D must be square and of equal order");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> R(m, m, zero);

    if (m == 1) {
        const T a = mu, b = T(C(0, 0) - mu), c = D(0, 0);
        const T disc = b * b - num_traits<T>::from_int(4) * a * c;
        if (disc < zero) throw NumericError("map_compute_R: no real solution for R");
        using std::sqrt;
        const T sd = T(sqrt(disc));
        const T r1 = T((-b - sd) / (num_traits<T>::from_int(2) * a));
        const T r2 = T((-b + sd) / (num_traits<T>::from_int(2) * a));
        if (r1 >= zero && r1 < one)
            R(0, 0) = r1;
        else if (r2 >= zero && r2 < one)
            R(0, 0) = r2;
        else
            throw NumericError("map_compute_R: no valid solution in [0,1) for R");
        return R;
    }

    Matrix<T> CmuI = C;
    for (std::size_t i = 0; i < m; ++i) CmuI(i, i) -= mu;
    R = matmul(D, inverse(CmuI));
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) R(i, j) = -R(i, j);

    Matrix<T> muImC(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) muImC(i, j) = -C(i, j);
    for (std::size_t i = 0; i < m; ++i) muImC(i, i) += mu;
    const Matrix<T> muImCinv = inverse(muImC);

    for (unsigned it = 0; it < max_iter; ++it) {
        const Matrix<T> R2 = matmul(R, R);
        Matrix<T> X(m, m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) X(i, j) = D(i, j) + mu * R2(i, j);
        const Matrix<T> Rn = matmul(X, muImCinv);
        T diff = zero;
        for (std::size_t i = 0; i < m; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < m; ++j) s += num_abs(T(Rn(i, j) - R(i, j)));
            if (s > diff) diff = s;
        }
        R = Rn;
        if (diff < tol) break;
    }
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
            if (R(i, j) < zero) R(i, j) = zero;
    return R;
}

/** map_compute_R_quadratic with the reference defaults, 5000 iterations, 1e-10. */
template <class T>
Matrix<T> map_compute_R_quadratic(const Matrix<T>& C, const Matrix<T>& D, const T& mu) {
    return map_compute_R_quadratic(C, D, mu, 5000u, T(num_traits<T>::from_double(1e-10)));
}

/**
 * The vectors h_{n,k} of Theorem 1 (map_m1ps_h_recursive.m).
 *
 * @return h[n][k], each of length M, for n = 0..N and k = 0..K
 */
template <class T>
std::vector<std::vector<std::vector<T>>> map_m1ps_h_recursive(const Matrix<T>& C,
                                                             const Matrix<T>& D, const T& mu,
                                                             std::size_t N, std::size_t K) {
    const std::size_t M = C.rows();
    if (C.cols() != M || D.rows() != M || D.cols() != M)
        throw InputError("map_m1ps_h_recursive: C and D must be square and of equal order");
    const T zero = num_traits<T>::from_int(0);

    T theta = zero;
    for (std::size_t i = 0; i < M; ++i) {
        const T a = num_abs(T(C(i, i)));
        if (a > theta) theta = a;
    }
    const T theta_plus_mu = theta + mu;
    if (theta_plus_mu == zero) throw NumericError("map_m1ps_h_recursive: theta + mu is zero");
    Matrix<T> thetaIplusC = C;
    for (std::size_t i = 0; i < M; ++i) thetaIplusC(i, i) += theta;

    std::vector<std::vector<std::vector<T>>> h(N + 1,
                                               std::vector<std::vector<T>>(K + 1, std::vector<T>()));
    for (std::size_t n = 0; n <= N; ++n) h[n][0] = ones<T>(M);

    for (std::size_t k = 0; k + 1 <= K; ++k) {
        for (std::size_t n = 0; n <= N; ++n) {
            std::vector<T> acc = mulvec(thetaIplusC, h[n][k]);
            if (n > 0) {
                const T c = num_traits<T>::from_int(static_cast<long>(n)) * mu /
                            num_traits<T>::from_int(static_cast<long>(n + 1));
                for (std::size_t i = 0; i < M; ++i) acc[i] += c * h[n - 1][k][i];
            }
            if (n < N) {
                const std::vector<T> t3 = mulvec(D, h[n + 1][k]);
                for (std::size_t i = 0; i < M; ++i) acc[i] += t3[i];
            }
            for (std::size_t i = 0; i < M; ++i) acc[i] /= theta_plus_mu;
            h[n][k + 1] = acc;
        }
    }
    return h;
}

/** What the two MAP/M/1-PS sojourn entry points return. */
template <class T>
struct MapM1psResult {
    std::vector<T> w_bar;  ///< Pr[W > x] at each requested point
    /// The reference's second output, sum_i (sum_k)_i / M, per level n and
    /// per point. NOT a conditional distribution; see the header note.
    std::vector<std::vector<T>> w_bar_n_unweighted;
    std::size_t n_levels;   ///< levels actually summed
    std::size_t k_max;      ///< largest uniformization index used
};

namespace m1ps_detail {

/** Poisson pmf a^k e^-a / k!, evaluated stably through logs. */
template <class T>
T poisson_pmf(const T& a, std::size_t k) {
    using std::exp;
    using std::log;
    if (a == num_traits<T>::from_int(0)) return k == 0 ? num_traits<T>::from_int(1)
                                                       : num_traits<T>::from_int(0);
    T lg = num_traits<T>::from_int(0);
    for (std::size_t j = 2; j <= k; ++j) lg += T(log(num_traits<T>::from_int(static_cast<long>(j))));
    const T lp = num_traits<T>::from_int(static_cast<long>(k)) * T(log(a)) - a - lg;
    return T(exp(lp));
}

/** The reference's L and R window for the Poisson weights at mean a. */
template <class T>
void poisson_window(const T& a, const T& eps_prime, std::size_t& L, std::size_t& K) {
    if (!(a > num_traits<T>::from_int(0))) {
        L = 0;
        K = 0;
        return;
    }
    const double ad = num_traits<T>::to_double(a);
    const double lo = ad - 10.0 * std::sqrt(ad);
    L = lo > 0.0 ? static_cast<std::size_t>(std::floor(lo)) : 0;
    std::size_t hi = static_cast<std::size_t>(std::ceil(ad + 10.0 * std::sqrt(ad)));
    const T one = num_traits<T>::from_int(1);
    for (;;) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t k = L; k <= hi; ++k) s += poisson_pmf(a, k);
        if (s >= one - eps_prime || hi >= 10000) break;
        hi += 10;
    }
    K = hi;
}

/** Stationary vector of the MAP phase process, pi (C + D) = 0. */
template <class T>
std::vector<T> m1ps_pi(const Matrix<T>& C, const Matrix<T>& D) {
    Matrix<T> Q(C.rows(), C.cols());
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < C.cols(); ++j) Q(i, j) = C(i, j) + D(i, j);
    return mc::ctmc_solve(Q);
}

}  // namespace m1ps_detail

/**
 * Complementary sojourn time distribution of a MAP/M/1-PS queue
 * (map_m1ps_sojourn.m).
 *
 * The queue-length truncation is the reference's: the smallest N with
 * (1/lambda) sum_{n<=N} pi_0 R^n D e > 1 - epsilon over n = 0..1000, falling
 * back to N = 100. The uniformization window is recomputed per evaluation
 * point, as in the reference, so the h recursion is rebuilt for each point.
 */
template <class T>
MapM1psResult<T> map_m1ps_sojourn(const Matrix<T>& C, const Matrix<T>& D, const T& mu,
                                  const std::vector<T>& x, const T& epsilon,
                                  const T& epsilon_prime) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_m1ps_sojourn requires transcendental arithmetic");
    const std::size_t M = C.rows();
    if (C.cols() != M || D.rows() != M || D.cols() != M)
        throw InputError("map_m1ps_sojourn: C and D must be square and of equal order");
    if (!(mu > num_traits<T>::from_int(0)))
        throw InputError("map_m1ps_sojourn: the service rate must be positive");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const std::vector<T> pi = m1ps_detail::m1ps_pi(C, D);
    const std::vector<T> e = ones<T>(M);
    T lambda = zero;
    {
        const std::vector<T> t = vecmul(pi, D);
        for (const T& v : t) lambda += v;
    }
    if (!(lambda > zero)) throw NumericError("map_m1ps_sojourn: zero arrival rate");
    if (!(T(lambda / mu) < one))
        throw NumericError("map_m1ps_sojourn: the system is unstable, rho >= 1");

    const Matrix<T> R = map_compute_R(C, D, mu);
    std::vector<T> pi0 = pi;
    {
        const std::vector<T> t = vecmul(pi, R);
        for (std::size_t i = 0; i < M; ++i) pi0[i] = pi[i] - t[i];
    }

    // N(epsilon) by the reference's cumulative scan.
    std::size_t N_epsilon = 0;
    {
        T cum = zero;
        std::vector<T> row = pi0;
        bool found = false;
        for (std::size_t n = 0; n <= 1000; ++n) {
            const std::vector<T> t = vecmul(row, D);
            T add = zero;
            for (const T& v : t) add += v;
            cum += add / lambda;
            if (cum > one - epsilon) {
                N_epsilon = n;
                found = true;
                break;
            }
            row = vecmul(row, R);
        }
        if (!found || N_epsilon == 0) N_epsilon = 100;
    }

    T theta = zero;
    for (std::size_t i = 0; i < M; ++i) {
        const T a = num_abs(T(C(i, i)));
        if (a > theta) theta = a;
    }
    const T theta_plus_mu = theta + mu;

    // Weights pi_0 R^n D, one row per level.
    std::vector<std::vector<T>> weights(N_epsilon + 1);
    {
        std::vector<T> row = pi0;
        for (std::size_t n = 0; n <= N_epsilon; ++n) {
            weights[n] = vecmul(row, D);
            row = vecmul(row, R);
        }
    }

    MapM1psResult<T> out;
    out.n_levels = N_epsilon + 1;
    out.k_max = 0;
    out.w_bar.assign(x.size(), zero);
    out.w_bar_n_unweighted.assign(N_epsilon + 1, std::vector<T>(x.size(), zero));

    for (std::size_t idx = 0; idx < x.size(); ++idx) {
        const T xv = x[idx];
        if (xv < zero) throw InputError("map_m1ps_sojourn: the evaluation points must be >= 0");
        std::size_t L = 0, K = 0;
        m1ps_detail::poisson_window(T(theta_plus_mu * xv), epsilon_prime, L, K);
        if (K > out.k_max) out.k_max = K;
        const std::vector<std::vector<std::vector<T>>> h =
            map_m1ps_h_recursive(C, D, mu, N_epsilon, K);

        std::vector<T> pois(K + 1 - L);
        for (std::size_t k = L; k <= K; ++k)
            pois[k - L] = m1ps_detail::poisson_pmf(T(theta_plus_mu * xv), k);

        for (std::size_t n = 0; n <= N_epsilon; ++n) {
            std::vector<T> sum_k(M, zero);
            for (std::size_t k = L; k <= K; ++k)
                for (std::size_t i = 0; i < M; ++i) sum_k[i] += pois[k - L] * h[n][k][i];
            T term = zero;
            for (std::size_t i = 0; i < M; ++i) term += weights[n][i] * sum_k[i];
            out.w_bar[idx] += term / lambda;
            T s = zero;
            for (const T& v : sum_k) s += v;
            out.w_bar_n_unweighted[n][idx] = s / num_traits<T>::from_int(static_cast<long>(M));
        }
    }
    (void)e;
    return out;
}

/** map_m1ps_sojourn with the reference defaults, epsilon 1e-11 and 1e-10. */
template <class T>
MapM1psResult<T> map_m1ps_sojourn(const Matrix<T>& C, const Matrix<T>& D, const T& mu,
                                  const std::vector<T>& x) {
    return map_m1ps_sojourn(C, D, mu, x, T(num_traits<T>::from_double(1e-11)),
                            T(num_traits<T>::from_double(1e-10)));
}

/**
 * Complementary sojourn time distribution of a MAP/M/1-PS queue by the
 * spectral-radius truncation (map_m1ps_cdfrespt.m).
 *
 * Differences from map_m1ps_sojourn are listed in the header note: R comes
 * from the other splitting, the level truncation is estimated from Sp(R) and
 * then cut again at ||pi_0 R^n D||_inf < epsilon/100, and the h recursion is
 * built ONCE at the largest uniformization index over all evaluation points
 * instead of once per point.
 *
 * Sp(R) is obtained from the caudal-characteristic bracket on the nonnegative
 * R (qbd_caudal) rather than from a double-precision eigensolve, so the
 * truncation estimate is computed at the working precision.
 */
template <class T>
MapM1psResult<T> map_m1ps_cdfrespt(const Matrix<T>& C, const Matrix<T>& D, const T& mu,
                                   const std::vector<T>& x, const T& epsilon,
                                   const T& epsilon_prime) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_m1ps_cdfrespt requires transcendental arithmetic");
    const std::size_t M = C.rows();
    if (C.cols() != M || D.rows() != M || D.cols() != M)
        throw InputError("map_m1ps_cdfrespt: C and D must be square and of equal order");
    if (!(mu > num_traits<T>::from_int(0)))
        throw InputError("map_m1ps_cdfrespt: the service rate must be positive");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    const std::vector<T> pi = m1ps_detail::m1ps_pi(C, D);
    T lambda = zero;
    {
        const std::vector<T> t = vecmul(pi, D);
        for (const T& v : t) lambda += v;
    }
    if (!(lambda > zero)) throw NumericError("map_m1ps_cdfrespt: zero arrival rate");
    if (!(T(lambda / mu) < one))
        throw NumericError("map_m1ps_cdfrespt: the system is unstable, rho >= 1");

    const Matrix<T> R = map_compute_R_quadratic(C, D, mu);
    std::vector<T> pi0(M);
    {
        const std::vector<T> t = vecmul(pi, R);
        for (std::size_t i = 0; i < M; ++i) pi0[i] = pi[i] - t[i];
    }

    const T rho_R = qbd_caudal(R);
    if (!(rho_R < one))
        throw NumericError("map_m1ps_cdfrespt: R has spectral radius >= 1");

    std::size_t N_epsilon;
    if (rho_R > zero) {
        const double e_ = num_traits<T>::to_double(epsilon);
        const double s_ = num_traits<T>::to_double(rho_R);
        const double est = std::ceil(std::log(e_ * (1.0 - s_)) / std::log(s_));
        const long n_est = static_cast<long>(est);
        N_epsilon = static_cast<std::size_t>(n_est < 10 ? 10 : (n_est > 10000 ? 10000 : n_est));
    } else {
        N_epsilon = 10;
    }

    T theta = zero;
    for (std::size_t i = 0; i < M; ++i) {
        const T a = num_abs(T(C(i, i)));
        if (a > theta) theta = a;
    }
    const T theta_plus_mu = theta + mu;

    // One uniformization window per point, and the global maximum over them.
    std::vector<std::size_t> Lp(x.size(), 0), Kp(x.size(), 0);
    std::size_t K_global = 0;
    for (std::size_t idx = 0; idx < x.size(); ++idx) {
        if (x[idx] < zero) throw InputError("map_m1ps_cdfrespt: the points must be >= 0");
        m1ps_detail::poisson_window(T(theta_plus_mu * x[idx]), epsilon_prime, Lp[idx], Kp[idx]);
        if (Kp[idx] > K_global) K_global = Kp[idx];
    }

    const std::vector<std::vector<std::vector<T>>> h =
        map_m1ps_h_recursive(C, D, mu, N_epsilon, K_global);

    // Weights with the reference's early truncation.
    const T weight_tol = epsilon * num_traits<T>::from_double(1e-2);
    std::vector<std::vector<T>> weights(N_epsilon + 1);
    std::size_t N_actual = N_epsilon;
    {
        std::vector<T> row = pi0;
        for (std::size_t n = 0; n <= N_epsilon; ++n) {
            weights[n] = vecmul(row, D);
            T wn = zero;
            for (const T& v : weights[n]) {
                const T a = num_abs(T(v));
                if (a > wn) wn = a;
            }
            if (n > 0 && wn < weight_tol) {
                N_actual = n;
                break;
            }
            row = vecmul(row, R);
        }
    }

    MapM1psResult<T> out;
    out.n_levels = N_actual + 1;
    out.k_max = K_global;
    out.w_bar.assign(x.size(), zero);
    out.w_bar_n_unweighted.assign(N_actual + 1, std::vector<T>(x.size(), zero));

    for (std::size_t idx = 0; idx < x.size(); ++idx) {
        const std::size_t L = Lp[idx], K = Kp[idx];
        std::vector<T> pois(K + 1 - L);
        for (std::size_t k = L; k <= K; ++k)
            pois[k - L] = m1ps_detail::poisson_pmf(T(theta_plus_mu * x[idx]), k);
        for (std::size_t n = 0; n <= N_actual; ++n) {
            std::vector<T> sum_k(M, zero);
            for (std::size_t k = L; k <= K; ++k)
                for (std::size_t i = 0; i < M; ++i) sum_k[i] += pois[k - L] * h[n][k][i];
            T term = zero;
            for (std::size_t i = 0; i < M; ++i) term += weights[n][i] * sum_k[i];
            out.w_bar[idx] += term / lambda;
            T s = zero;
            for (const T& v : sum_k) s += v;
            out.w_bar_n_unweighted[n][idx] = s / num_traits<T>::from_int(static_cast<long>(M));
        }
    }
    return out;
}

/** map_m1ps_cdfrespt with the reference defaults, epsilon 1e-11 and 1e-10. */
template <class T>
MapM1psResult<T> map_m1ps_cdfrespt(const Matrix<T>& C, const Matrix<T>& D, const T& mu,
                                   const std::vector<T>& x) {
    return map_m1ps_cdfrespt(C, D, mu, x, T(num_traits<T>::from_double(1e-11)),
                             T(num_traits<T>::from_double(1e-10)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_M1PS_H
