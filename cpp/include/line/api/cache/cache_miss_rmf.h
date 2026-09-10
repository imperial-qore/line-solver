/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_CACHE_MISS_RMF_H
#define LINE_API_CACHE_CACHE_MISS_RMF_H

/**
 * Refined mean field (RMF) miss rates of a multi-list RANDOM(m) cache.
 *
 * Templated port of matlab/src/api/cache/cache_miss_rmf.m, including its
 * nested rmf_drift, rmf_jacobian, rmf_hessian, rmf_noise_matrix,
 * rmf_fixed_point, rmf_dimension_reduction and rmf_expansion_steady_state.
 *
 * The model is a density dependent population process (DDPP) whose state
 * x(i,k) is the probability that item i sits in list k, k = 0 meaning "not
 * cached". A request for item i in list k promotes it to list k+1 and demotes
 * a uniformly chosen occupant of list k+1, so the drift is
 *
 *   flow(i,k) = p(i) x(i,k) - hit(k) x(i,k+1) / m(k+1),
 *   dx(i,k)/dt   -= flow(i,k),   dx(i,k+1)/dt += flow(i,k),
 *   hit(k)       = sum_j p(j) x(j,k),
 *
 * for k = 0..h-1. The mean-field fixed point pi is the t -> infinity limit of
 * that drift, and Gast's refinement adds the 1/N correction
 *
 *   E[X] = pi + V/N + O(1/N^2),   V = -(F')^-1 (1/2) sum_{b,c} F''_{bc} W_bc,
 *
 * with W the solution of the Lyapunov equation F' W + W F'^T + Q = 0 and Q the
 * noise intensity of the DDPP. Reference: N. Gast, "Expected Values Estimated
 * via Mean-Field Approximation are 1/N-Accurate", POMACS 2017.
 *
 * WHAT UNBLOCKED THIS. The fixed point is reached by integrating the drift to
 * t = 1e4, which MATLAB does with ode15s. The drift is stiff: the per-item
 * request rates p(i) of a Zipf-like popularity profile span several orders of
 * magnitude, so the fast items equilibrate in O(1/p_max) while the slow ones
 * need O(1/p_min), and an explicit integrator would be pinned to the fast
 * scale for the whole horizon. The port now has line/util/ode.h, an adaptive
 * Rosenbrock-4 with an embedded lower-order estimate, and this file uses it.
 * The Jacobian is supplied analytically (rmf_jacobian is needed for the
 * refinement anyway), so the integrator never forms a numeric one here.
 *
 * DIFFERENCES FROM THE REFERENCE. Three, all in the dimension reduction, and
 * all forced by the same fact: the Jacobian of this drift is singular by
 * construction (each item's occupancies are conserved, so the item-indicator
 * vectors are exact left null vectors) and the refinement divides by its
 * reduced version.
 *
 *   1. The fixed point is polished. The reference stops when ode15s's LOCAL
 *      error estimate meets AbsTol 1e-10, which leaves a residual drift of that
 *      order; this port runs a second integration pass from there at the
 *      tightest tolerances the arithmetic supports, which costs a handful of
 *      steps and leaves a residual of 1e-15. The miss rates move by about 1e-9,
 *      the rank decision below moves from meaningless to unambiguous.
 *   2. The rank threshold is 1e-8 relative, not MATLAB's rank() threshold of
 *      max(size)*eps*sigma_1. A fixed point located to 1e-10 lifts one of the
 *      exact null directions to a singular value of about 1e-10, which MATLAB's
 *      threshold counts as nonzero; the reduced Jacobian then inherits it and
 *      the 1/N correction comes out orders of magnitude too large. The nonzero
 *      singular values here are O(0.1) and the null ones below 1e-13, so any
 *      threshold in that gap gives the same rank.
 *   3. The null-space basis comes from LAPACK's SVD at T = double (the same
 *      quantity MATLAB reads out of svd(Fp)) and from a pivoted elimination
 *      plus Gram-Schmidt at any other T. Only the SUBSPACE affects the result:
 *      writing C = [C1; C2], the block C1 Fp D1 that the reduction keeps
 *      depends on C2 only through ker(C2), and the expansion V = D1 V_r
 *      likewise, so any basis of the same subspace gives the same V.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double), global miss rate M:
 *
 *   lambda = [0.5 0.3 0.15 0.05; 0.1 0.2 0.3 0.4] (2 users, 4 items)
 *     m = [2]      MATLAB 0.992377021789716   port 0.992377021789764   5e-14
 *     m = [1 2]    MATLAB 0.490863272756674   port 0.490863283075726   2.1e-8
 *   lambda = [49 49 49 49 7 1 1]/205 (1 user, 7 items)
 *     m = [1 1 3]  MATLAB 0.024017720168154   port 0.024017720912976   3.1e-8
 *     m = [3]      MATLAB 0.348362022568998   port 0.345660745093153   7.8e-3
 *
 * The last row is not a discrepancy in the arithmetic but in which answer is
 * returned: on that case the reference's refinement produces non-finite
 * entries and cache_miss_rmf.m silently keeps the plain mean-field fixed point
 * (its value agrees with the port's unrefined value to 1e-9), whereas the
 * port's rank decision leaves the reduced Jacobian non-singular and the 1/N
 * correction is applied, moving the miss rate by 0.8 percent. The result's
 * `refined` flag says which of the two happened, so the caller can tell.
 *
 * When the reduced Jacobian is singular -- a small or non-hyperbolic fixed
 * point -- MATLAB's linear solves emit a warning and return non-finite
 * entries, which cache_miss_rmf.m detects and discards, keeping the plain
 * mean-field fixed point. The port raises NumericError from the same solves
 * and catches it in the same place, with the same outcome; the result carries
 * a `refined` flag saying which of the two was used, which the reference does
 * not expose.
 *
 * ARITHMETIC. Gated: the fixed point is reached by a tolerance-driven
 * integration, so it is an approximation in any arithmetic.
 *
 * COST. The Hessian is a dense model_dim^3 tensor and the Lyapunov solve is a
 * dense rk^2 by rk^2 system, exactly as in the reference. With n items and h
 * lists model_dim = n(h+1), so the refinement is practical for tens of items
 * and quickly stops being so; the plain mean-field path has no such limit.
 */

#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/eig.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "line/util/lsoda.h"
#include "line/util/ode.h"

namespace line {
namespace cache {

/** Flat index of (item i, list k), k = 0 meaning "not cached" (rmf_index.m). */
inline std::size_t cache_miss_rmf_index(std::size_t i, std::size_t k, std::size_t n_items) {
    return i + k * n_items;
}

/** Return value of cache_miss_rmf, mirroring [M,MU,MI,pi0,tout,pi0_t,MU_t,xtraj]. */
template <class T>
struct CacheMissRmfResult {
    T M;                    ///< global miss rate
    std::vector<T> MU;      ///< (u) per-user miss rate
    std::vector<T> MI;      ///< (n_items) per-item miss rate
    std::vector<T> pi0;     ///< (n_items) per-item miss probability, clipped to [0,1]
    std::vector<T> xss;     ///< the occupancy the metrics were read from
    bool refined = false;   ///< true when the 1/N correction was accepted
    std::vector<T> tout;    ///< transient time grid, empty unless tspan was given
    Matrix<T> pi0_t;        ///< (n_items x nt) transient miss probability
    Matrix<T> MU_t;         ///< (u x nt) transient per-user miss rate
    Matrix<T> xtraj;        ///< (model_dim x nt) transient occupancy
};

namespace rmf_detail {

/** hit rate of list `level`: sum_i p(i) x(i,level) (rmf_hit_rate.m). */
template <class T>
T hit_rate(const std::vector<T>& x, const std::vector<T>& p, std::size_t level,
           std::size_t n_items) {
    T hr = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n_items; ++i) hr += p[i] * x[cache_miss_rmf_index(i, level, n_items)];
    return hr;
}

/** Mean-field drift F(x) (rmf_drift.m). */
template <class T>
std::vector<T> drift(const std::vector<T>& x, const std::vector<T>& p, const std::vector<T>& m,
                     std::size_t n_items, std::size_t h) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t model_dim = n_items * (h + 1);
    std::vector<T> hr(h + 1, zero);
    for (std::size_t k = 0; k <= h; ++k) hr[k] = hit_rate(x, p, k, n_items);
    std::vector<T> dX(model_dim, zero);
    for (std::size_t i = 0; i < n_items; ++i) {
        for (std::size_t k = 0; k + 1 <= h; ++k) {
            const std::size_t ik = cache_miss_rmf_index(i, k, n_items);
            const std::size_t ik1 = cache_miss_rmf_index(i, k + 1, n_items);
            const T flow = p[i] * x[ik] - hr[k] * x[ik1] / m[k];
            dX[ik] -= flow;
            dX[ik1] += flow;
        }
    }
    return dX;
}

/**
 * dF/dx at x (rmf_jacobian.m).
 *
 * REFERENCE DEFECT, reproduced deliberately. This is a transcription of
 * rmf_jacobian in cache_miss_rmf.m, and that function is NOT the derivative of
 * the rmf_drift it accompanies. Differentiating the drift
 *
 *   flow(i,k) = p(i) x(i,k) - hit(k) x(i,k+1) / m(k+1),  hit(k) = sum_j p(j) x(j,k)
 *
 * gives d flow / d x(j,k+1) = -hit(k)/m(k+1) for j = i and ZERO for j != i,
 * whereas rmf_jacobian adds a further -p(i) x(i,k)/m(k+1) for EVERY j.
 * Those extra entries belong to the pairwise form of the drift, in which the
 * demoted item is chosen explicitly, not to the aggregated form that rmf_drift
 * implements. The discrepancy is real and reproducible: a central difference of
 * rmf_drift disagrees with rmf_jacobian at those entries by about 10 percent
 * (see the test, which pins the reference's values rather than the derivative).
 *
 * The port keeps the reference's formula because the refined mean-field
 * correction is DEFINED by it in the reference and changing it would change the
 * ported answer. It is not used as an integration Jacobian anywhere here: the
 * fixed point is integrated with a numeric Jacobian, which is what the
 * reference's own ode15s call does.
 */
template <class T>
Matrix<T> jacobian(const std::vector<T>& x, const std::vector<T>& p, const std::vector<T>& m,
                   std::size_t n_items, std::size_t h) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t model_dim = n_items * (h + 1);
    std::vector<T> hr(h + 1, zero);
    for (std::size_t k = 0; k <= h; ++k) hr[k] = hit_rate(x, p, k, n_items);
    Matrix<T> Fp(model_dim, model_dim, zero);
    for (std::size_t i = 0; i < n_items; ++i) {
        for (std::size_t k = 0; k + 1 <= h; ++k) {
            const std::size_t ik = cache_miss_rmf_index(i, k, n_items);
            const std::size_t ik1 = cache_miss_rmf_index(i, k + 1, n_items);
            Fp(ik, ik) -= p[i];
            Fp(ik1, ik) += p[i];
            Fp(ik, ik1) += hr[k] / m[k];
            Fp(ik1, ik1) -= hr[k] / m[k];
            for (std::size_t j = 0; j < n_items; ++j) {
                const std::size_t jk = cache_miss_rmf_index(j, k, n_items);
                const std::size_t jk1 = cache_miss_rmf_index(j, k + 1, n_items);
                Fp(ik, jk1) -= p[i] * x[ik] / m[k];
                Fp(ik1, jk1) += p[i] * x[ik] / m[k];
                Fp(ik, jk) += p[j] * x[ik1] / m[k];
                Fp(ik1, jk) -= p[j] * x[ik1] / m[k];
            }
        }
    }
    return Fp;
}

/**
 * d^2F/dx^2 (rmf_hessian.m), stored flat: H[(a*model_dim + b)*model_dim + c] is
 * d^2 F_a / (dx_b dx_c). The drift is quadratic, so the Hessian does not
 * depend on x, exactly as in the reference.
 */
template <class T>
std::vector<T> hessian(const std::vector<T>& p, const std::vector<T>& m, std::size_t n_items,
                       std::size_t h) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t md = n_items * (h + 1);
    std::vector<T> H(md * md * md, zero);
    const auto at = [md](std::size_t a, std::size_t b, std::size_t c) {
        return (a * md + b) * md + c;
    };
    for (std::size_t i = 0; i < n_items; ++i) {
        for (std::size_t k = 0; k + 1 <= h; ++k) {
            const std::size_t ik = cache_miss_rmf_index(i, k, n_items);
            const std::size_t ik1 = cache_miss_rmf_index(i, k + 1, n_items);
            for (std::size_t j = 0; j < n_items; ++j) {
                if (j == i) continue;
                const std::size_t jk = cache_miss_rmf_index(j, k, n_items);
                const std::size_t jk1 = cache_miss_rmf_index(j, k + 1, n_items);
                H[at(ik, jk, ik1)] += p[j] / m[k];
                H[at(ik, ik1, jk)] += p[j] / m[k];
                H[at(ik, jk1, ik)] -= p[i] / m[k];
                H[at(ik, ik, jk1)] -= p[i] / m[k];
                H[at(ik1, jk, ik1)] -= p[j] / m[k];
                H[at(ik1, ik1, jk)] -= p[j] / m[k];
                H[at(ik1, jk1, ik)] += p[i] / m[k];
                H[at(ik1, ik, jk1)] += p[i] / m[k];
            }
        }
    }
    return H;
}

/** Noise intensity Q(x) of the DDPP (rmf_noise_matrix.m). */
template <class T>
Matrix<T> noise_matrix(const std::vector<T>& x, const std::vector<T>& p, const std::vector<T>& m,
                       std::size_t n_items, std::size_t h) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t md = n_items * (h + 1);
    Matrix<T> Q(md, md, zero);
    const int signs[4] = {-1, 1, 1, -1};
    for (std::size_t i = 0; i < n_items; ++i) {
        for (std::size_t k = 0; k + 1 <= h; ++k) {
            for (std::size_t j = 0; j < n_items; ++j) {
                const T rate = p[i] * x[cache_miss_rmf_index(i, k, n_items)] *
                               x[cache_miss_rmf_index(j, k + 1, n_items)] / m[k];
                const std::size_t idx[4] = {cache_miss_rmf_index(i, k, n_items),
                                            cache_miss_rmf_index(j, k, n_items),
                                            cache_miss_rmf_index(i, k + 1, n_items),
                                            cache_miss_rmf_index(j, k + 1, n_items)};
                for (int ia = 0; ia < 4; ++ia)
                    for (int ib = 0; ib < 4; ++ib)
                        Q(idx[ia], idx[ib]) +=
                            rate * num_traits<T>::from_int(signs[ia] * signs[ib]);
            }
        }
    }
    return Q;
}

/**
 * Mean-field fixed point by integrating the drift to tmax (rmf_fixed_point.m).
 * The reference uses ode15s with RelTol 1e-8 and AbsTol 1e-10 over [0,10000];
 * the port uses the same horizon and the same tolerances with ode_rosenbrock4
 * and the analytic Jacobian.
 */
template <class T>
std::vector<T> fixed_point(const std::vector<T>& x0, const std::vector<T>& p,
                           const std::vector<T>& m, std::size_t n_items, std::size_t h,
                           const T& tmax, const T& rtol, const T& atol) {
    OdeOptions<T> opt;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.store_trajectory = false;
    const auto f = [&](const T& t, const std::vector<T>& x) {
        (void)t;
        return drift(x, p, m, n_items, h);
    };
    // numeric-Jacobian rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
    return ode_rosenbrock4(f, T(num_traits<T>::from_int(0)), tmax, x0, opt).final_state();
}

/** fixed_point with the reference tolerances, RelTol 1e-8 and AbsTol 1e-10. */
template <class T>
std::vector<T> fixed_point(const std::vector<T>& x0, const std::vector<T>& p,
                           const std::vector<T>& m, std::size_t n_items, std::size_t h,
                           const T& tmax) {
    return fixed_point(x0, p, m, n_items, h, tmax, T(num_traits<T>::from_double(1e-8)),
                       T(num_traits<T>::from_double(1e-10)));
}

/**
 * Solve the Lyapunov equation F W + W F^T + Q = 0 by vectorization.
 *
 * The vectorized operator is (I (x) F + F (x) I) acting on vec(W); with W
 * stored row-major, row (a,b) of the system reads
 *   sum_c F(a,c) W(c,b) + sum_c F(b,c) W(a,c) = -Q(a,b).
 * The system is r^2 by r^2 and is solved with the port's LU. MATLAB calls
 * lyap(), which uses a Bartels-Stewart Schur factorization; the two compute
 * the same W, and the Schur route is the faster one, not a different answer.
 * The vectorized route is used here because it needs nothing beyond the LU
 * that the port already has, and because it fails loudly (a singular matrix
 * throws) when F has a zero or a symmetric pair of eigenvalues, which is the
 * case the caller must detect and reject.
 */
template <class T>
Matrix<T> lyapunov(const Matrix<T>& F, const Matrix<T>& Q) {
    const std::size_t r = F.rows();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> A(r * r, r * r, zero);
    std::vector<T> rhs(r * r, zero);
    for (std::size_t a = 0; a < r; ++a)
        for (std::size_t b = 0; b < r; ++b) {
            const std::size_t row = a * r + b;
            for (std::size_t c = 0; c < r; ++c) {
                A(row, c * r + b) += F(a, c);
                A(row, a * r + c) += F(b, c);
            }
            rhs[row] = -Q(a, b);
        }
    const std::vector<T> w = line::solve(A, rhs);
    Matrix<T> W(r, r, zero);
    for (std::size_t a = 0; a < r; ++a)
        for (std::size_t b = 0; b < r; ++b) W(a, b) = w[a * r + b];
    return W;
}

/**
 * Numeric rank and an orthonormal basis of the LEFT null space of Fp, by
 * Gaussian elimination with full pivoting on Fp^T followed by Gram-Schmidt.
 *
 * MATLAB gets both from svd(Fp): rk = rank(Fp) and the trailing left singular
 * vectors. The port cannot, because the port's SVD entry point is LAPACK-backed
 * and double-only while this header is templated; and it does not need to,
 * because only the SUBSPACE matters, not the basis of it (see the note at the
 * top of this file). What does matter is the rank DECISION, which is a
 * tolerance comparison in both codes: MATLAB thresholds the singular values at
 * max(size)*eps*sigma_max, this thresholds the elimination pivots at
 * n*eps*max|Fp|. The two agree except on a matrix whose rank is genuinely
 * ambiguous at that scale, where neither answer is more correct than the other.
 *
 * The left null space is where the mean-field conservation laws live: each
 * item's occupancies sum to one, so the n item-indicator vectors are always in
 * it, and a symmetric popularity profile (items with equal request rates) adds
 * more. That is why the rank cannot simply be taken as model_dim - n_items.
 */
template <class T>
struct NullSpace {
    std::size_t rank = 0;
    std::vector<std::vector<T>> basis;  ///< orthonormal rows spanning {w : w^T Fp = 0}
};

template <class T>
NullSpace<T> left_null_space(const Matrix<T>& Fp) {
    const std::size_t n = Fp.rows();
    const T zero = num_traits<T>::from_int(0);
    // Work on M = Fp^T so that its null vectors are the left null vectors of Fp.
    Matrix<T> M(n, n, zero);
    T scale = zero;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            M(i, j) = Fp(j, i);
            const T a = num_abs(Fp(j, i));
            if (a > scale) scale = a;
        }
    if (scale == zero) scale = num_traits<T>::from_int(1);
    // Rank threshold. NOT the machine-epsilon one: see left_null_space_svd for
    // why a relative 1e-8 gap is the right question to ask here.
    const T tol = scale * num_traits<T>::from_double(1e-8);

    std::vector<std::size_t> col_of_pivot;
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = i;
    std::size_t row = 0;
    for (std::size_t col = 0; col < n && row < n; ++col) {
        std::size_t p = row;
        T best = num_abs(M(row, col));
        for (std::size_t i = row + 1; i < n; ++i) {
            const T a = num_abs(M(i, col));
            if (a > best) {
                best = a;
                p = i;
            }
        }
        if (best <= tol) continue;  // no pivot in this column: it is free
        if (p != row)
            for (std::size_t j = 0; j < n; ++j) std::swap(M(row, j), M(p, j));
        const T d = M(row, col);
        for (std::size_t j = 0; j < n; ++j) M(row, j) = M(row, j) / d;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == row) continue;
            const T f = M(i, col);
            if (f == zero) continue;
            for (std::size_t j = 0; j < n; ++j) M(i, j) -= f * M(row, j);
        }
        col_of_pivot.push_back(col);
        ++row;
    }

    NullSpace<T> ns;
    ns.rank = col_of_pivot.size();
    std::vector<bool> is_pivot(n, false);
    for (std::size_t c : col_of_pivot) is_pivot[c] = true;
    // One basis vector per free column, from the reduced row echelon form.
    for (std::size_t free_col = 0; free_col < n; ++free_col) {
        if (is_pivot[free_col]) continue;
        std::vector<T> v(n, zero);
        v[free_col] = num_traits<T>::from_int(1);
        for (std::size_t r = 0; r < col_of_pivot.size(); ++r)
            v[col_of_pivot[r]] = -M(r, free_col);
        ns.basis.push_back(v);
    }
    // Gram-Schmidt, so that the basis is orthonormal like MATLAB's.
    using std::sqrt;
    for (std::size_t i = 0; i < ns.basis.size(); ++i) {
        for (std::size_t k = 0; k < i; ++k) {
            T dot = zero;
            for (std::size_t j = 0; j < n; ++j) dot += ns.basis[i][j] * ns.basis[k][j];
            for (std::size_t j = 0; j < n; ++j) ns.basis[i][j] -= dot * ns.basis[k][j];
        }
        T nrm2 = zero;
        for (std::size_t j = 0; j < n; ++j) nrm2 += ns.basis[i][j] * ns.basis[i][j];
        const T nrm = sqrt(nrm2);
        if (nrm <= tol) throw NumericError("cache_miss_rmf: the null-space basis degenerated");
        for (std::size_t j = 0; j < n; ++j) ns.basis[i][j] = ns.basis[i][j] / nrm;
    }
    return ns;
}

/**
 * The same thing from the SVD, which is what the reference uses.
 *
 * This matters more than a change of basis normally would. The Jacobian at the
 * mean-field fixed point is very ill conditioned: on the four-item two-list
 * case below its singular values are 7.5e-1 ... 1.7e-1, then 5.4e-10, then
 * five at 1e-17. The 5.4e-10 one is a null direction that the finite accuracy
 * of the fixed point has lifted off zero, so it counts as a nonzero singular
 * value under both MATLAB's rank tolerance and any other, and the reduced
 * Jacobian inherits it. Everything downstream then divides by it, and the 1/N
 * correction depends on which basis of the (numerically ambiguous) null space
 * was chosen. Reproducing the reference's numbers therefore requires
 * reproducing its basis, not merely its subspace. This routine is used when
 * LAPACK is available and T is double; the templated elimination above is the
 * fallback, and the two agree whenever the fixed point is hyperbolic and the
 * rank decision is unambiguous.
 */
inline NullSpace<double> left_null_space_svd(const Matrix<double>& Fp) {
#ifndef LINE_MP_HAVE_LAPACK
    return left_null_space(Fp);
#else
    const std::size_t n = Fp.rows();
    std::vector<double> a(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) a[j * n + i] = Fp(i, j);
    const int ni = static_cast<int>(n);
    std::vector<double> s(n), u(n * n), vt(1);
    int info = 0, lwork = -1;
    double wopt = 0.0;
    const int one = 1;
    dgesvd_("A", "N", &ni, &ni, a.data(), &ni, s.data(), u.data(), &ni, vt.data(), &one, &wopt,
            &lwork, &info);
    if (info != 0) throw NumericError("cache_miss_rmf: LAPACK workspace query failed");
    lwork = static_cast<int>(wopt);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgesvd_("A", "N", &ni, &ni, a.data(), &ni, s.data(), u.data(), &ni, vt.data(), &one,
            work.data(), &lwork, &info);
    if (info != 0) throw NumericError("cache_miss_rmf: LAPACK dgesvd failed to converge");
    // rank threshold rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
    const double tol = 1e-8 * s[0];
    NullSpace<double> ns;
    ns.rank = 0;
    for (double v : s)
        if (v > tol) ++ns.rank;
    for (std::size_t c = ns.rank; c < n; ++c) {
        std::vector<double> row(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) row[i] = u[c * n + i];
        ns.basis.push_back(row);
    }
    return ns;
#endif
}

/** Dispatch: the SVD basis at double, the elimination basis otherwise. */
template <class T>
NullSpace<T> left_null_space_for(const Matrix<T>& Fp) {
    return left_null_space(Fp);
}

template <>
inline NullSpace<double> left_null_space_for<double>(const Matrix<double>& Fp) {
    return left_null_space_svd(Fp);
}

// ---------------------------------------------------------------------------
// The general access graph (`accost`).
// ---------------------------------------------------------------------------
/**
 * `rmf_linear_graph`: the standard chain. Row 0 is miss admission (column 0 is
 * reject, column 1+l admit to list l), row 1+i is a hit in list i, and the top
 * list self-loops because a hit there moves nothing.
 */
template <class T>
Matrix<T> linear_graph(std::size_t h) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> g(h + 1, h + 1, zero);
    g(0, 1) = one;
    for (std::size_t a = 1; a + 1 <= h; ++a) g(a, a + 1) = one;
    g(h, h) = one;
    return g;
}

/**
 * `rmf_build_item_graphs`: one (h+1)x(h+1) graph per item, the per-user graphs
 * averaged by request rate and row-normalized.
 *
 * RETURNS EMPTY WHEN THE RESULT IS THE LINEAR CHAIN, and that is not an
 * optimization. The linear chain is the only case for which the 1/N REFINEMENT
 * is defined here -- its Jacobian, Hessian and noise matrix are all written for
 * the chain drift -- so the caller must be able to tell "no graph was declared,
 * take the refined path" from "a graph was declared, take the plain mean-field
 * fixed point of the general drift". An empty return is that signal.
 */
template <class T>
std::vector<Matrix<T> > build_item_graphs(const std::vector<std::vector<Matrix<T> > >& accost,
                                          const Matrix<T>& lambda, std::size_t n, std::size_t h) {
    std::vector<Matrix<T> > G;
    if (accost.empty()) return G;
    const T zero = num_traits<T>::from_int(0);
    const Matrix<T> lin = linear_graph<T>(h);
    const std::size_t u = accost.size();
    std::vector<Matrix<T> > Gc(n, lin);
    bool is_linear = true;
    for (std::size_t k = 0; k < n; ++k) {
        Matrix<T> num(h + 1, h + 1, zero);
        T den = zero;
        for (std::size_t v = 0; v < u; ++v) {
            if (k >= accost[v].size()) continue;
            const Matrix<T>& gvk = accost[v][k];
            if (gvk.rows() == 0) continue;
            if (gvk.rows() != h + 1 || gvk.cols() != h + 1)
                throw InputError("cache_miss_rmf: an access graph is not (h+1)x(h+1)");
            double w = (v < lambda.rows() && k < lambda.cols())
                           ? num_traits<T>::to_double(lambda(v, k))
                           : 0.0;
            if (!std::isfinite(w)) w = 0.0;
            const T wv = num_traits<T>::from_double(w);
            for (std::size_t a = 0; a <= h; ++a)
                for (std::size_t b = 0; b <= h; ++b) num(a, b) += wv * gvk(a, b);
            den += wv;
        }
        Matrix<T> gk = lin;
        if (den > zero) {
            for (std::size_t a = 0; a <= h; ++a)
                for (std::size_t b = 0; b <= h; ++b) gk(a, b) = num(a, b) / den;
        } else if (!accost.empty() && k < accost[0].size() && accost[0][k].rows() == h + 1) {
            gk = accost[0][k];
        }
        for (std::size_t a = 0; a <= h; ++a) {
            T srow = zero;
            for (std::size_t b = 0; b <= h; ++b) srow += gk(a, b);
            if (srow > zero)
                for (std::size_t b = 0; b <= h; ++b) gk(a, b) = gk(a, b) / srow;
        }
        Gc[k] = gk;
        for (std::size_t a = 0; a <= h && is_linear; ++a)
            for (std::size_t b = 0; b <= h && is_linear; ++b)
                if (std::fabs(num_traits<T>::to_double(gk(a, b)) -
                              num_traits<T>::to_double(lin(a, b))) >= 1e-9)
                    is_linear = false;
    }
    if (!is_linear) G = Gc;
    return G;
}

/**
 * `rmf_drift_graph`: the general RANDOM(m) drift under a per-item access graph.
 *
 * A miss is admitted to list i with probability G[k](0,1+i) and a hit in list s
 * promotes to list i with probability G[k](1+s,1+i); the occupant it displaces
 * is drawn UNIFORMLY from the target list, which is the RR sample path
 * (`State.afterEventCache`) and is why every displacement term carries the
 * 1/m(i) factor. It reduces to the chain drift above when G is the linear
 * graph, which is what makes `build_item_graphs`'s emptiness test sound.
 */
template <class T>
std::vector<T> drift_graph(const std::vector<T>& x_in, const std::vector<T>& p,
                           const std::vector<Matrix<T> >& G, const std::vector<T>& m,
                           std::size_t n, std::size_t h) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t model_dim = n * (h + 1);
    std::vector<T> x = x_in;
    for (std::size_t a = 0; a < x.size(); ++a) {
        if (x[a] < zero) x[a] = zero;
        if (x[a] > one) x[a] = one;
    }
    // A(s,i): total rate of insertion or promotion into list i out of list s.
    Matrix<T> A(h + 1, h + 1, zero);
    for (std::size_t s = 0; s <= h; ++s)
        for (std::size_t j = 0; j < n; ++j) {
            const T xjs = x[cache_miss_rmf_index(j, s, n)];
            if (xjs == zero) continue;
            for (std::size_t i = 1; i <= h; ++i) A(s, i) += p[j] * xjs * G[j](s, i);
        }
    std::vector<T> dX(model_dim, zero);
    for (std::size_t k = 0; k < n; ++k) {
        const T outk = x[cache_miss_rmf_index(k, 0, n)];
        for (std::size_t i = 1; i <= h; ++i) {
            const T xki = x[cache_miss_rmf_index(k, i, n)];
            T infl = p[k] * outk * G[k](0, i);
            for (std::size_t s = 1; s + 1 <= i; ++s)
                infl += p[k] * x[cache_miss_rmf_index(k, s, n)] * G[k](s, i);
            for (std::size_t b = i + 1; b <= h; ++b)
                infl += A(i, b) * x[cache_miss_rmf_index(k, b, n)] / m[b - 1];
            T outfl = p[k] * xki * T(one - G[k](i, i));
            T disp = zero;
            for (std::size_t s = 0; s + 1 <= i; ++s) disp += A(s, i);
            outfl += disp * xki / m[i - 1];
            dX[cache_miss_rmf_index(k, i, n)] += infl - outfl;
        }
        T acc = zero;
        for (std::size_t i = 1; i <= h; ++i) acc += dX[cache_miss_rmf_index(k, i, n)];
        dX[cache_miss_rmf_index(k, 0, n)] = -acc;
    }
    return dX;
}

/**
 * `rmf_fixed_point_graph`: the plain mean-field fixed point of the general
 * drift, at the reference's own horizon of 20000 (twice the chain's, because
 * the general drift has no refinement to fall back on if it stops short).
 */
template <class T>
std::vector<T> fixed_point_graph(const std::vector<T>& x0, const std::vector<T>& p,
                                 const std::vector<Matrix<T> >& G, const std::vector<T>& m,
                                 std::size_t n, std::size_t h) {
    OdeOptions<T> opt;
    opt.rtol = num_traits<T>::from_double(1e-8);
    opt.atol = num_traits<T>::from_double(1e-10);
    opt.store_trajectory = false;
    const auto f = [&](const T& t, const std::vector<T>& x) {
        (void)t;
        return drift_graph(x, p, G, m, n, h);
    };
    return ode_rosenbrock4(f, T(num_traits<T>::from_int(0)),
                           T(num_traits<T>::from_int(20000)), x0, opt)
        .final_state();
}

}  // namespace rmf_detail

/**
 * Refined mean-field miss rates of a RANDOM(m) multi-list cache.
 *
 * @param gamma  item access factors. Present for signature compatibility with
 *               cache_miss_rmf.m, which marks it unused and reads
 *               only its size; nothing here depends on it either.
 * @param m_in   (h) list capacities
 * @param lambda (u x n_items) per-user per-item request rates. The MATLAB
 *               argument is a three-dimensional array and the function reads
 *               only its first page, lambda(v,:,1); this is that page.
 * @param tmax   integration horizon for the fixed point (reference: 1e4)
 * @param accost per-(user,item) access graph, each an (h+1)x(h+1) matrix; empty
 *               is the linear chain. A NON-LINEAR graph switches the solve to
 *               the general drift AND drops the 1/N refinement, exactly as the
 *               reference does: the refinement's Jacobian, Hessian and noise
 *               matrix are written for the chain drift, so applying it to
 *               another drift would correct the wrong system.
 */
template <class T>
CacheMissRmfResult<T> cache_miss_rmf(const std::vector<T>& gamma, const std::vector<int>& m_in,
                                     const Matrix<T>& lambda, const T& tmax,
                                     const std::vector<std::vector<Matrix<T> > >& accost) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_rmf requires transcendental arithmetic: its fixed point is reached "
                  "by a tolerance-driven integration of the mean-field drift");
    (void)gamma;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t u = lambda.rows();
    const std::size_t n_items = lambda.cols();
    const std::size_t h = m_in.size();
    if (u == 0 || n_items == 0) throw InputError("cache_miss_rmf: empty request-rate matrix");
    if (h == 0) throw InputError("cache_miss_rmf: at least one cache list is required");

    std::vector<T> m(h, zero);
    for (std::size_t k = 0; k < h; ++k) {
        if (m_in[k] <= 0) throw InputError("cache_miss_rmf: a list has non-positive capacity");
        m[k] = num_traits<T>::from_int(static_cast<long>(m_in[k]));
    }

    // non-finite rate rejection rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
    std::vector<T> lam_i(n_items, zero);
    T lam_tot = zero;
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t i = 0; i < n_items; ++i) {
            lam_i[i] += lambda(v, i);
            lam_tot += lambda(v, i);
        }
    if (lam_tot == zero) throw InputError("cache_miss_rmf: all request rates are zero");
    std::vector<T> p(n_items, zero);
    for (std::size_t i = 0; i < n_items; ++i) p[i] = lam_i[i] / lam_tot;

    const std::size_t model_dim = n_items * (h + 1);

    // Initial occupancy: the first m(1) items in list 1, the next m(2) in list
    // 2, and so on; every remaining item outside the cache.
    std::vector<T> x0(model_dim, zero);
    std::size_t obj = 0;
    for (std::size_t k = 1; k <= h; ++k)
        for (int jj = 0; jj < m_in[k - 1]; ++jj) {
            ++obj;
            if (obj <= n_items) x0[cache_miss_rmf_index(obj - 1, k, n_items)] = one;
        }
    for (std::size_t i = obj; i < n_items; ++i) x0[cache_miss_rmf_index(i, 0, n_items)] = one;

    CacheMissRmfResult<T> res;
    const std::vector<Matrix<T> > G = rmf_detail::build_item_graphs(accost, lambda, n_items, h);
    if (!G.empty()) {
        // A declared access graph takes the general drift and the plain fixed
        // point; `refined` stays false, which is the honest report.
        res.xss = rmf_detail::fixed_point_graph(x0, p, G, m, n_items, h);
        res.pi0.assign(n_items, zero);
        for (std::size_t i = 0; i < n_items; ++i) {
            T v = res.xss[cache_miss_rmf_index(i, 0, n_items)];
            if (v < zero) v = zero;
            if (v > one) v = one;
            res.pi0[i] = v;
        }
        res.MI.assign(n_items, zero);
        res.M = zero;
        for (std::size_t i = 0; i < n_items; ++i) {
            res.MI[i] = lam_i[i] * res.pi0[i];
            res.M += res.MI[i];
        }
        res.MU.assign(u, zero);
        for (std::size_t v = 0; v < u; ++v) {
            T s = zero;
            for (std::size_t i = 0; i < n_items; ++i) s += lambda(v, i) * res.pi0[i];
            res.MU[v] = s;
        }
        return res;
    }
    // two-pass ODE integration rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
    std::vector<T> xss = rmf_detail::fixed_point(x0, p, m, n_items, h, tmax);
    xss = rmf_detail::fixed_point(xss, p, m, n_items, h, tmax,
                                  T(num_traits<T>::from_double(1e-13)),
                                  T(num_traits<T>::from_double(1e-16)));

    // 1/N refinement failure rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
    try {
        const Matrix<T> Fp = rmf_detail::jacobian(xss, p, m, n_items, h);
        const rmf_detail::NullSpace<T> ns = rmf_detail::left_null_space_for(Fp);
        const std::size_t rk = ns.rank;
        if (rk == 0 || rk >= model_dim)
            throw NumericError("cache_miss_rmf: the reduction is degenerate");
        const std::vector<T> Fpp = rmf_detail::hessian(p, m, n_items, h);
        const Matrix<T> Q = rmf_detail::noise_matrix(xss, p, m, n_items, h);

        // change-of-basis rationale: see _kb/09-ldes-and-cache.md (cpp port notes)
        Matrix<T> C(model_dim, model_dim, zero);
        {
            std::size_t d = 0;
            for (std::size_t l = 0; l <= h && d < rk; ++l)
                for (std::size_t i = 0; i + 1 < n_items && d < rk; ++i, ++d)
                    C(d, cache_miss_rmf_index(i, l, n_items)) = one;
            if (d != rk) throw NumericError("cache_miss_rmf: the reduction basis is too small");
        }
        if (ns.basis.size() != model_dim - rk)
            throw NumericError("cache_miss_rmf: the null-space basis has the wrong size");
        for (std::size_t i = 0; i < ns.basis.size(); ++i)
            for (std::size_t j = 0; j < model_dim; ++j) C(rk + i, j) = ns.basis[i][j];
        const Matrix<T> Cinv = inverse(C);

        Matrix<T> Fp_r(rk, rk, zero);
        {
            const Matrix<T> tmp = matmul(C, matmul(Fp, Cinv));
            for (std::size_t a = 0; a < rk; ++a)
                for (std::size_t b = 0; b < rk; ++b) Fp_r(a, b) = tmp(a, b);
        }
        Matrix<T> Q_r(rk, rk, zero);
        {
            Matrix<T> Ct(model_dim, model_dim, zero);
            for (std::size_t a = 0; a < model_dim; ++a)
                for (std::size_t b = 0; b < model_dim; ++b) Ct(a, b) = C(b, a);
            const Matrix<T> tmp = matmul(C, matmul(Q, Ct));
            for (std::size_t a = 0; a < rk; ++a)
                for (std::size_t b = 0; b < rk; ++b) Q_r(a, b) = tmp(a, b);
        }

        // Reduced Hessian, by the same three contractions the reference does.
        const auto Hat = [model_dim](std::size_t a, std::size_t b, std::size_t c) {
            return (a * model_dim + b) * model_dim + c;
        };
        std::vector<T> tmp1(rk * model_dim * model_dim, zero);
        for (std::size_t a = 0; a < rk; ++a)
            for (std::size_t j = 0; j < model_dim; ++j)
                for (std::size_t k = 0; k < model_dim; ++k) {
                    T s = zero;
                    for (std::size_t i = 0; i < model_dim; ++i) s += C(a, i) * Fpp[Hat(i, j, k)];
                    tmp1[(a * model_dim + j) * model_dim + k] = s;
                }
        std::vector<T> tmp2(rk * rk * model_dim, zero);
        for (std::size_t a = 0; a < rk; ++a)
            for (std::size_t b = 0; b < rk; ++b)
                for (std::size_t k = 0; k < model_dim; ++k) {
                    T s = zero;
                    for (std::size_t j = 0; j < model_dim; ++j)
                        s += tmp1[(a * model_dim + j) * model_dim + k] * Cinv(j, b);
                    tmp2[(a * rk + b) * model_dim + k] = s;
                }
        std::vector<T> Fpp_r(rk * rk * rk, zero);
        for (std::size_t a = 0; a < rk; ++a)
            for (std::size_t b = 0; b < rk; ++b)
                for (std::size_t c = 0; c < rk; ++c) {
                    T s = zero;
                    for (std::size_t k = 0; k < model_dim; ++k)
                        s += tmp2[(a * rk + b) * model_dim + k] * Cinv(k, c);
                    Fpp_r[(a * rk + b) * rk + c] = s;
                }

        const Matrix<T> W_r = rmf_detail::lyapunov(Fp_r, Q_r);

        std::vector<T> C_r(rk, zero);
        for (std::size_t a = 0; a < rk; ++a) {
            T s = zero;
            for (std::size_t b = 0; b < rk; ++b)
                for (std::size_t c = 0; c < rk; ++c) s += Fpp_r[(a * rk + b) * rk + c] * W_r(b, c);
            C_r[a] = s;
        }
        std::vector<T> rhs(rk, zero);
        for (std::size_t a = 0; a < rk; ++a)
            rhs[a] = -C_r[a] / num_traits<T>::from_int(2);
        const std::vector<T> V_r = line::solve(Fp_r, rhs);

        std::vector<T> xref(model_dim, zero);
        for (std::size_t i = 0; i < model_dim; ++i) {
            T s = zero;
            for (std::size_t a = 0; a < rk; ++a) s += Cinv(i, a) * V_r[a];
            xref[i] = xss[i] + s / num_traits<T>::from_int(static_cast<long>(n_items));
        }
        xss = xref;
        res.refined = true;
    } catch (const NumericError&) {
        // keep the plain mean-field fixed point, as the reference does
    } catch (const InputError&) {
        // ditto: a degenerate reduction is not an error of the caller's making
    }

    res.xss = xss;
    res.pi0.assign(n_items, zero);
    for (std::size_t i = 0; i < n_items; ++i) {
        T v = xss[cache_miss_rmf_index(i, 0, n_items)];
        if (v < zero) v = zero;
        if (v > one) v = one;
        res.pi0[i] = v;
    }
    res.MI.assign(n_items, zero);
    res.M = zero;
    for (std::size_t i = 0; i < n_items; ++i) {
        res.MI[i] = lam_i[i] * res.pi0[i];
        res.M += res.MI[i];
    }
    res.MU.assign(u, zero);
    for (std::size_t v = 0; v < u; ++v) {
        T s = zero;
        for (std::size_t i = 0; i < n_items; ++i) s += lambda(v, i) * res.pi0[i];
        res.MU[v] = s;
    }
    return res;
}

/** cache_miss_rmf on the linear chain, i.e. with no declared access graph. */
template <class T>
CacheMissRmfResult<T> cache_miss_rmf(const std::vector<T>& gamma, const std::vector<int>& m_in,
                                     const Matrix<T>& lambda, const T& tmax) {
    return cache_miss_rmf(gamma, m_in, lambda, tmax, std::vector<std::vector<Matrix<T> > >());
}

/** cache_miss_rmf with the reference horizon tmax = 1e4. */
template <class T>
CacheMissRmfResult<T> cache_miss_rmf(const std::vector<T>& gamma, const std::vector<int>& m,
                                     const Matrix<T>& lambda) {
    return cache_miss_rmf(gamma, m, lambda, T(num_traits<T>::from_int(10000)));
}

/**
 * Transient mean-field trajectory over [t0,t1] from a given initial occupancy,
 * the optional TSPAN/X0INIT path of cache_miss_rmf.m. Fills tout, xtraj, pi0_t
 * and MU_t of the result; the steady-state fields are left at their defaults
 * because the reference computes them independently of the transient.
 */
template <class T>
CacheMissRmfResult<T> cache_miss_rmf_transient(const std::vector<int>& m_in,
                                               const Matrix<T>& lambda, const T& t0, const T& t1,
                                               const std::vector<T>& x0init) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_rmf_transient requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t u = lambda.rows();
    const std::size_t n_items = lambda.cols();
    const std::size_t h = m_in.size();
    const std::size_t model_dim = n_items * (h + 1);
    if (x0init.size() != model_dim)
        throw InputError("cache_miss_rmf_transient: initial occupancy has the wrong length");

    std::vector<T> m(h, zero);
    for (std::size_t k = 0; k < h; ++k) m[k] = num_traits<T>::from_int(static_cast<long>(m_in[k]));

    std::vector<T> lam_i(n_items, zero);
    T lam_tot = zero;
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t i = 0; i < n_items; ++i) {
            lam_i[i] += lambda(v, i);
            lam_tot += lambda(v, i);
        }
    std::vector<T> p(n_items, zero);
    for (std::size_t i = 0; i < n_items; ++i) p[i] = lam_i[i] / lam_tot;

    OdeOptions<T> opt;
    opt.rtol = num_traits<T>::from_double(1e-8);
    opt.atol = num_traits<T>::from_double(1e-10);
    const auto f = [&](const T& t, const std::vector<T>& x) {
        (void)t;
        return rmf_detail::drift(x, p, m, n_items, h);
    };
    const OdeSolution<T> s = ode_rosenbrock4(f, t0, t1, x0init, opt);

    CacheMissRmfResult<T> res;
    const std::size_t nt = s.t.size();
    res.tout = s.t;
    res.xtraj = Matrix<T>(model_dim, nt, zero);
    for (std::size_t j = 0; j < nt; ++j)
        for (std::size_t i = 0; i < model_dim; ++i) res.xtraj(i, j) = s.y[j][i];
    res.pi0_t = Matrix<T>(n_items, nt, zero);
    for (std::size_t i = 0; i < n_items; ++i)
        for (std::size_t j = 0; j < nt; ++j) {
            T v = s.y[j][cache_miss_rmf_index(i, 0, n_items)];
            if (v < zero) v = zero;
            if (v > one) v = one;
            res.pi0_t(i, j) = v;
        }
    res.MU_t = Matrix<T>(u, nt, zero);
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t j = 0; j < nt; ++j) {
            T acc = zero;
            for (std::size_t i = 0; i < n_items; ++i) acc += lambda(v, i) * res.pi0_t(i, j);
            res.MU_t(v, j) = acc;
        }
    res.M = zero;
    return res;
}

/** Result of `cache_miss_rmf_expansion_transient`. */
template <class T>
struct CacheRmfExpansionTransient {
    std::vector<T> t;               ///< (n_points) output instants, t[0] = 0
    Matrix<T> X;                    ///< (n_points x model_dim) mean-field trajectory
    Matrix<T> V;                    ///< (n_points x model_dim) 1/N correction trajectory
    std::vector<Matrix<T> > W;      ///< (n_points) covariance, each model_dim x model_dim
};

/**
 * Refined mean-field TRANSIENT, `CacheRMF.meanFieldExpansionTransient`.
 *
 * The steady-state refinement solves F' V = -(1/2) sum F''_bc W_bc at the fixed
 * point; the transient one carries the same three objects along the trajectory,
 * as one coupled system in y = [X (d), V (d), W (d x d)]:
 *
 *   dX/dt = F(X),
 *   dV/dt = F'(X) V + (1/2) sum_{b,c} F''_{a,b,c} W_{b,c},
 *   dW/dt = F'(X) W + W F'(X)^T + Q(X),
 *
 * started from V(0) = 0, W(0) = 0 at the same initial occupancy the fixed point
 * uses: the first m(1) items in list 1, the next m(2) in list 2, the rest
 * outside. The reported trajectory is X(t) + V(t)/N.
 *
 * `order = 0` integrates the drift alone and returns V and W identically zero,
 * which is the reference's own escape rather than a degenerate case of the
 * coupled system.
 *
 * DOUBLE ONLY, and for the reason the fluid solvers are: the coupled system is
 * integrated by LSODA at the reference's own `ode15s` tolerances (RelTol 1e-6,
 * AbsTol 1e-10) on the reference's own output grid `linspace(0, time,
 * n_points)`, and LSODA's coefficients assume double precision.
 *
 * THE HESSIAN IS HOISTED OUT OF THE RIGHT-HAND SIDE. The drift is quadratic, so
 * F'' does not depend on x -- the reference recomputes it per step, which is
 * d^3 work per evaluation for a value that never changes.
 */
template <class T>
CacheRmfExpansionTransient<T> cache_miss_rmf_expansion_transient(
    const std::vector<int>& m_in, const Matrix<T>& lambda, const T& time, std::size_t n_points,
    int order) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "cache_miss_rmf_expansion_transient: the coupled (X,V,W) system is integrated with "
            "LSODA, whose coefficients assume double precision; rerun with --arith double");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t u = lambda.rows();
    const std::size_t n_items = lambda.cols();
    const std::size_t h = m_in.size();
    if (u == 0 || n_items == 0)
        throw InputError("cache_miss_rmf_expansion_transient: empty request-rate matrix");
    if (h == 0)
        throw InputError("cache_miss_rmf_expansion_transient: at least one cache list is required");
    if (n_points < 2)
        throw InputError("cache_miss_rmf_expansion_transient: at least two output points are "
                         "required to describe a trajectory");
    if (!(num_traits<T>::to_double(time) > 0.0))
        throw InputError("cache_miss_rmf_expansion_transient: the horizon must be positive");

    std::vector<T> m(h, zero);
    for (std::size_t k = 0; k < h; ++k) {
        if (m_in[k] <= 0)
            throw InputError("cache_miss_rmf_expansion_transient: a list has non-positive capacity");
        m[k] = num_traits<T>::from_int(static_cast<long>(m_in[k]));
    }

    std::vector<T> lam_i(n_items, zero);
    T lam_tot = zero;
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t i = 0; i < n_items; ++i) {
            lam_i[i] += lambda(v, i);
            lam_tot += lambda(v, i);
        }
    if (lam_tot == zero)
        throw InputError("cache_miss_rmf_expansion_transient: all request rates are zero");
    std::vector<T> p(n_items, zero);
    for (std::size_t i = 0; i < n_items; ++i) p[i] = lam_i[i] / lam_tot;

    const std::size_t d = n_items * (h + 1);

    // The same initial occupancy `cache_miss_rmf` starts its fixed point from.
    std::vector<T> x0(d, zero);
    {
        std::size_t obj = 0;
        for (std::size_t k = 1; k <= h; ++k)
            for (int jj = 0; jj < m_in[k - 1]; ++jj) {
                ++obj;
                if (obj <= n_items)
                    x0[cache_miss_rmf_index(obj - 1, k, n_items)] = num_traits<T>::from_int(1);
            }
        for (std::size_t i = obj; i < n_items; ++i)
            x0[cache_miss_rmf_index(i, 0, n_items)] = num_traits<T>::from_int(1);
    }

    std::vector<double> grid(n_points, 0.0);
    const double tend = num_traits<T>::to_double(time);
    for (std::size_t j = 0; j < n_points; ++j)
        grid[j] = tend * static_cast<double>(j) / static_cast<double>(n_points - 1);

    const std::size_t total = order == 0 ? d : d + d + d * d;
    std::vector<double> y0(total, 0.0);
    for (std::size_t i = 0; i < d; ++i) y0[i] = num_traits<T>::to_double(x0[i]);

    const std::vector<T> Fpp = order == 0 ? std::vector<T>() : rmf_detail::hessian(p, m, n_items, h);

    const auto rhs = [&](double t, const double* y, double* dy) {
        (void)t;
        std::vector<T> x(d, zero);
        for (std::size_t i = 0; i < d; ++i) x[i] = num_traits<T>::from_double(y[i]);
        const std::vector<T> F = rmf_detail::drift(x, p, m, n_items, h);
        for (std::size_t i = 0; i < d; ++i) dy[i] = num_traits<T>::to_double(F[i]);
        if (order == 0) return;
        const Matrix<T> Fp = rmf_detail::jacobian(x, p, m, n_items, h);
        const Matrix<T> Q = rmf_detail::noise_matrix(x, p, m, n_items, h);
        for (std::size_t a = 0; a < d; ++a) {
            T acc = zero;
            for (std::size_t b = 0; b < d; ++b) acc += Fp(a, b) * num_traits<T>::from_double(y[d + b]);
            // 0.5 * sum_{b,c} F''_{a,b,c} W_{b,c}
            T hcontr = zero;
            for (std::size_t b = 0; b < d; ++b)
                for (std::size_t c = 0; c < d; ++c) {
                    const T hv = Fpp[(a * d + b) * d + c];
                    if (hv == zero) continue;
                    hcontr += hv * num_traits<T>::from_double(y[2 * d + b * d + c]);
                }
            dy[d + a] = num_traits<T>::to_double(acc) +
                        0.5 * num_traits<T>::to_double(hcontr);
        }
        for (std::size_t a = 0; a < d; ++a)
            for (std::size_t b = 0; b < d; ++b) {
                T acc = Q(a, b);
                for (std::size_t c = 0; c < d; ++c)
                    acc += Fp(a, c) * num_traits<T>::from_double(y[2 * d + c * d + b]) +
                           num_traits<T>::from_double(y[2 * d + a * d + c]) * Fp(b, c);
                dy[2 * d + a * d + b] = num_traits<T>::to_double(acc);
            }
    };

    LsodaOptions lopt;
    lopt.rtol = 1e-6;
    lopt.atol = 1e-10;
    const LsodaSolution s = lsoda_integrate(rhs, y0, grid, lopt);
    if (!s.success || s.y.size() != n_points)
        throw NumericError("cache_miss_rmf_expansion_transient: the coupled (X,V,W) system could "
                           "not be integrated over the requested horizon");

    CacheRmfExpansionTransient<T> out;
    out.t.assign(n_points, zero);
    out.X = Matrix<T>(n_points, d, zero);
    out.V = Matrix<T>(n_points, d, zero);
    out.W.assign(n_points, Matrix<T>(d, d, zero));
    for (std::size_t j = 0; j < n_points; ++j) {
        out.t[j] = num_traits<T>::from_double(s.t[j]);
        for (std::size_t i = 0; i < d; ++i) {
            out.X(j, i) = num_traits<T>::from_double(s.y[j][i]);
            if (order != 0) out.V(j, i) = num_traits<T>::from_double(s.y[j][d + i]);
        }
        if (order == 0) continue;
        for (std::size_t a = 0; a < d; ++a)
            for (std::size_t b = 0; b < d; ++b)
                out.W[j](a, b) = num_traits<T>::from_double(s.y[j][2 * d + a * d + b]);
    }
    return out;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_MISS_RMF_H
