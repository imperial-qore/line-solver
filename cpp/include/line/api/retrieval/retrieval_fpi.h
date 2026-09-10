/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_FPI_H
#define LINE_API_RETRIEVAL_RETRIEVAL_FPI_H

/**
 * Fixed-point heuristic for a delayed-hit (list-based) cache.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_fpi.m, cross-checked
 * against jar/src/main/java/jline/api/retrieval/Retrieval_fpi.java.
 *
 * The exact recursion of retrieval_metrics costs O(2^n); the heuristic
 * truncates the perturbation expansion at zeroth order,
 * pi_{i,l}(m-1_j) ~ pi_{i,l}(m), and solves the resulting nonlinear system by
 * successive substitution. One sweep is
 *
 *   F_{s,i}   = 1 + sum_{k!=i} phi_{s,k}
 *   D_i       = 1 + lambda_i eta_{0,i} + sum_s lambda_i eta_{s,i} F_{s,i}
 *   theta_ij  = gamma_ij / D_i
 *   xi_j      = m_j / sum_k theta_kj (1 - sum_l pi_kl)
 *   pi_ij     = theta_ij xi_j / (1 + sum_l theta_il xi_l)
 *   pi_i0     = (1 - sum_j pi_ij) / D_i
 *   phi_{0,i} = lambda_i eta_{0,i} pi_i0,  phi_{s,i} = lambda_i eta_{s,i} F_{s,i} pi_i0
 *
 * stopping when the largest relative change of the miss, hit and delayed-hit
 * ratios falls below tol. Note that xi and the "1 - sum_l pi_kl" factor read
 * the PREVIOUS sweep's pi, exactly as MATLAB and the JAR do, so the port is
 * iterate-for-iterate identical to them and not merely fixed-point identical.
 *
 * ARITHMETIC: the iteration stops on a tolerance, so its answer is the fixed
 * point only to within tol whatever the arithmetic; it is gated on
 * has_transcendental so nobody instantiates it at exact arithmetic expecting
 * an exact result. Use retrieval_metrics or retrieval_mva for that.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** Options mirroring the trailing (max_iter, tol) arguments of the MATLAB function. */
struct FpiOptions {
    std::size_t max_iter = 1000;
    double tol = 1e-6;
};

/** Mirrors the [pmiss, phit, pdh] return list, plus the iteration diagnostics. */
template <class T>
struct RetrievalFpiResult {
    std::vector<T> pmiss;  ///< (n) miss ratios pi_{i,0}
    Matrix<T> phit;        ///< (h x n) hit ratios pi_{i,j}
    Matrix<T> pdh;         ///< ((r+1) x n) delayed-hit probabilities phi_{s,i}, s = 0..r
    std::size_t iterations = 0;
    bool converged = false;
    /// True when the iterate left the finite range; MATLAB warns and breaks here.
    bool diverged = false;
};

namespace detail {

/** MATLAB's local reldiff: max|a-b| over max|b|, with an all-zero b read as 1. */
template <class T>
double fpi_reldiff(const std::vector<T>& a, const std::vector<T>& b) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = num_traits<T>::to_double(num_abs(T(a[i] - b[i])));
        const double e = num_traits<T>::to_double(num_abs(T(b[i])));
        if (d > num) num = d;
        if (e > den) den = e;
    }
    if (den == 0.0) den = 1.0;
    return num / den;
}

}  // namespace detail

/**
 * @param m      (h) cache list capacities
 * @param lambda (n) per-item arrival rates
 * @param eta    (n x (r+1)) fetching demands, column 0 = IS station, columns 1..r = PS
 * @param gamma  (n x h) access factors
 * @param options fixed-point options (tolerance, iteration cap, damping)
 */
template <class T>
RetrievalFpiResult<T> retrieval_fpi(const std::vector<int>& m, const std::vector<T>& lambda,
                                    const Matrix<T>& eta, const Matrix<T>& gamma,
                                    const FpiOptions& options = FpiOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "retrieval_fpi requires transcendental arithmetic: it is a successive "
                  "substitution stopped on a relative tolerance, so its answer is the fixed "
                  "point only to within tol whatever the arithmetic");
    const std::size_t n = lambda.size();
    const std::size_t h = m.size();
    if (eta.rows() != n || gamma.rows() != n)
        throw InputError("retrieval_fpi: eta/gamma and lambda disagree on the item count");
    if (gamma.cols() != h)
        throw InputError("retrieval_fpi: gamma and m disagree on the number of lists");
    if (eta.cols() == 0) throw InputError("retrieval_fpi: eta has no columns");
    const std::size_t r = eta.cols() - 1;

    const T one = num_traits<T>::from_int(1);

    // initial guess of the paper: the h+1 cache states and the r+2 retrieval
    // states are each given the same mass.
    const T init_phi = one / num_traits<T>::from_int(static_cast<long>((h + 1) * (r + 2)));
    const T init_pij = one / num_traits<T>::from_int(static_cast<long>(h + 1));

    // phi is stored (r+1) x n and pij is stored h x n, as MATLAB returns them.
    Matrix<T> phi(r + 1, n, init_phi);
    Matrix<T> pij(h, n, init_pij);
    std::vector<T> pi0(n, init_phi);

    RetrievalFpiResult<T> out;
    for (std::size_t t = 1; t <= options.max_iter; ++t) {
        out.iterations = t;

        // F(s,i) = 1 + sum_{k != i} phi_{s,k}
        Matrix<T> F(r, n, one);
        for (std::size_t s = 0; s < r; ++s) {
            T tot = num_traits<T>::from_int(0);
            for (std::size_t i = 0; i < n; ++i) tot += phi(s + 1, i);
            for (std::size_t i = 0; i < n; ++i) F(s, i) = one + (tot - phi(s + 1, i));
        }

        std::vector<T> D(n);
        for (std::size_t i = 0; i < n; ++i) {
            D[i] = one + lambda[i] * eta(i, 0);
            for (std::size_t s = 0; s < r; ++s) D[i] += lambda[i] * eta(i, s + 1) * F(s, i);
        }

        Matrix<T> theta(n, h);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < h; ++j) theta(i, j) = gamma(i, j) / D[i];

        // 1 - sum_l pi_kl, read off the previous sweep as MATLAB does
        std::vector<T> oneminus(n);
        for (std::size_t k = 0; k < n; ++k) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < h; ++j) s += pij(j, k);
            oneminus[k] = one - s;
        }

        std::vector<T> xi(h);
        for (std::size_t j = 0; j < h; ++j) {
            T den = num_traits<T>::from_int(0);
            for (std::size_t k = 0; k < n; ++k) den += theta(k, j) * oneminus[k];
            xi[j] = num_traits<T>::from_int(static_cast<long>(m[j])) / den;
        }

        Matrix<T> pij_new(h, n);
        std::vector<T> pi0_new(n);
        Matrix<T> phi_new(r + 1, n);
        for (std::size_t i = 0; i < n; ++i) {
            T denom = one;
            for (std::size_t l = 0; l < h; ++l) denom += theta(i, l) * xi[l];
            T sum_pij = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < h; ++j) {
                pij_new(j, i) = theta(i, j) * xi[j] / denom;
                sum_pij += pij_new(j, i);
            }
            pi0_new[i] = (one - sum_pij) / D[i];
            phi_new(0, i) = lambda[i] * eta(i, 0) * pi0_new[i];
            for (std::size_t s = 0; s < r; ++s)
                phi_new(s + 1, i) = lambda[i] * eta(i, s + 1) * F(s, i) * pi0_new[i];
        }

        const std::vector<T> pij_flat(pij_new.data(), pij_new.data() + pij_new.size());
        const std::vector<T> pij_old(pij.data(), pij.data() + pij.size());
        const std::vector<T> phi_flat(phi_new.data(), phi_new.data() + phi_new.size());
        const std::vector<T> phi_old(phi.data(), phi.data() + phi.size());
        double delta = detail::fpi_reldiff(pi0_new, pi0);
        const double d2 = detail::fpi_reldiff(pij_flat, pij_old);
        const double d3 = detail::fpi_reldiff(phi_flat, phi_old);
        if (d2 > delta) delta = d2;
        if (d3 > delta) delta = d3;

        pij = pij_new;
        pi0 = pi0_new;
        phi = phi_new;

        if (!std::isfinite(delta)) {
            out.diverged = true;
            break;
        }
        if (delta < options.tol) {
            out.converged = true;
            break;
        }
    }

    out.pmiss = pi0;
    out.phit = pij;
    out.pdh = phi;
    return out;
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_FPI_H
