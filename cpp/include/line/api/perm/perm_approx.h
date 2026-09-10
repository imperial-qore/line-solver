/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PERM_PERM_APPROX_H
#define LINE_API_PERM_PERM_APPROX_H

/**
 * APPROXIMATE permanents: the Sinkhorn heuristic, the Bethe estimate and the
 * saddle-point expansion.
 *
 * Port of python/line_solver/api/perm/approx.py, itself a twin of MATLAB's
 * `perm_heur.m` and of `jline.lib.perm.BethePermanent`.
 *
 * WHY APPROXIMATE AT ALL. The exact algorithms in `permanent.h` cost 2^n or
 * n!, and the multiplicity method only escapes that when columns repeat. On a
 * dense matrix with distinct columns neither is usable past about thirty, and
 * these two are what remain.
 *
 * THE THREE ARE NOT INTERCHANGEABLE, and a caller has to know which guarantee it
 * is getting:
 *
 *  - `perm_heur` is a HEURISTIC with no error bound in either direction. It
 *    Sinkhorn-scales the matrix toward double stochasticity, averages a
 *    mean-field van der Waerden estimate with a Gurvits capacity bound on the
 *    scaled matrix, and undoes the scaling. The average of a lower bound and a
 *    mean-field estimate is neither.
 *  - `perm_bethe` is a LOWER BOUND on a STRICTLY POSITIVE matrix, which is a
 *    real guarantee and the reason to prefer it when one is needed. It runs
 *    sum-product message passing on the square root of the matrix to a fixed
 *    point and exponentiates the Bethe free energy there.
 *  - `perm_spm` is the saddle-point (Laplace) expansion of the coefficient
 *    integral, the HOMOGENEOUS variant of `cache_spm`. It is exact in the limit
 *    of large column multiplicities and, at unit multiplicities, overestimates
 *    by a factor near (e/sqrt(2 pi))^n with a tight spread across matrices. No
 *    bound in either direction, but far closer than the bare capacity it
 *    corrects, and it is the one of the three that takes repeated columns.
 *
 * ALL REQUIRE A STRICTLY POSITIVE MATRIX and refuse otherwise, by name and
 * with the offending position. Two separate reasons:
 *
 *  - a NEGATIVE entry: Sinkhorn scaling diverges rather than failing, and the
 *    Bethe bound is simply not a bound off the nonnegative orthant;
 *  - a ZERO entry: the matrix has no full support. Both routines used to floor
 *    a zero to a small eps first, and that substitution is NOT INVERTIBLE.
 *    Every permutation takes one entry per row, so the floored matrix has
 *    permanent n! eps times the permanent of the rest against a true permanent
 *    that may be 0; n! outruns eps by n = 18, and the order of the replicated
 *    demand matrix in pfqn_jointmarg is sum(N). The floored answers were also
 *    simply wrong: on a 3x3 of ones with two zeroed entries, whose permanent is
 *    3, perm_bethe returned 2748880111.1018 -- the lower-bound property gone by
 *    nine orders of magnitude -- and perm_heur returned the van der Waerden
 *    bound of the Sinkhorn limit of the FLOORED matrix.
 *
 * Positivity is sufficient but not necessary. The sharp precondition of the
 * Sinkhorn scaling is TOTAL SUPPORT: every positive entry lies on a positive
 * permutation. A matrix with a strictly positive permanent can still fail it --
 * [[J3, 0], [J3, J3]] has permanent 36 and no total support, and the scaling
 * then stalls on its tolerance instead of converging. Positivity is the test
 * used because it is O(n^2). Use the exact engine for a matrix with zeros.
 *
 * ARITHMETIC: double. All are iterative floating-point schemes.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace perm {

namespace approxdetail {

/**
 * n! in double.
 *
 * Exact to n = 170, the largest factorial a double holds; past that Stirling
 * avoids the overflow, which is the reference's own rule.
 */
inline double factorial_d(std::size_t n) {
    if (n <= 1) return 1.0;
    if (n <= 170) {
        double v = 1.0;
        for (std::size_t i = 2; i <= n; ++i) v *= static_cast<double>(i);
        return v;
    }
    const double d = static_cast<double>(n);
    return std::sqrt(2.0 * M_PI * d) * std::pow(d / std::exp(1.0), d);
}

/** Refuse a matrix that is not square and nonnegative, naming the entry. */
inline void require_nonnegative_square(const Matrix<double>& m, const char* who) {
    const std::size_t n = m.rows();
    if (n == 0 || m.cols() != n)
        throw InputError(std::string(who) + ": the matrix must be square and non-empty");
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (m(i, j) < 0.0)
                throw InputError(std::string(who) + ": the matrix must be non-negative; entry (" +
                                 std::to_string(i) + ", " + std::to_string(j) + ") is " +
                                 std::to_string(m(i, j)));
}

/**
 * Log determinant of a symmetric positive definite k x k matrix, by Cholesky.
 *
 * A failed factorisation is a degenerate saddle point, not a rounding accident,
 * so it is reported rather than nudged.
 */
inline double log_det_cholesky(const std::vector<double>& h, std::size_t k, const char* who) {
    std::vector<double> l(k * k, 0.0);
    double logdet = 0.0;
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double acc = h[i * k + j];
            for (std::size_t t = 0; t < j; ++t) acc -= l[i * k + t] * l[j * k + t];
            if (i == j) {
                if (!(acc > 0.0))
                    throw InputError(std::string(who) + ": the reduced Hessian is not positive"
                                     " definite (leading minor " + std::to_string(i + 1) +
                                     "), so the saddle point is degenerate and the Gaussian"
                                     " factor does not exist");
                l[i * k + i] = std::sqrt(acc);
                logdet += 2.0 * std::log(l[i * k + i]);
            } else {
                l[i * k + j] = acc / l[j * k + j];
            }
        }
    }
    return logdet;
}

}  // namespace approxdetail

/**
 * Sinkhorn scaling toward double stochasticity.
 *
 * Returns the scaled matrix and the two diagonal scalings, since undoing them
 * is what turns an estimate on the scaled matrix back into one on the original.
 */
/**
 * Refuse a matrix the permanent approximations cannot take.
 *
 * The approximations rest, directly or through the Sinkhorn scaling they
 * share, on a strictly positive matrix. A zero used to be floored to a small
 * eps first, and that substitution is not invertible: a matrix with an
 * identically zero row has permanent 0 while the floored matrix has permanent
 * n! eps times the permanent of the rest, which is O(1) by n = 18. Positivity
 * is sufficient but not necessary -- the sharp precondition is TOTAL SUPPORT,
 * which a matrix with a positive permanent can still fail -- but it is O(n^2)
 * and is the contract this header already states.
 */
inline void require_full_support(const Matrix<double>& m, const char* who) {
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (!(m(i, j) > 0.0))
                throw InputError(std::string(who) +
                                 ": requires a strictly positive matrix, but entry (" +
                                 std::to_string(i + 1) + "," + std::to_string(j + 1) +
                                 ") is " + std::to_string(m(i, j)) +
                                 ", so the matrix has no full support. Flooring it would change"
                                 " the permanent by n!*eps, which is O(1) by n=18. Use the exact"
                                 " engine.");
}

inline void sinkhorn_scaling(const Matrix<double>& in, Matrix<double>* B, std::vector<double>* r,
                             std::vector<double>* c, double tolerance = 1e-10,
                             std::size_t max_iterations = 1000) {
    const std::size_t n = in.rows();
    r->assign(n, 1.0);
    c->assign(n, 1.0);
    bool converged = false;
    double last_error = std::numeric_limits<double>::infinity();
    for (std::size_t it = 0; it < max_iterations; ++it) {
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < n; ++j) s += in(i, j) * (*c)[j];
            (*r)[i] = (s != 0.0) ? 1.0 / s : 0.0;
        }
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) s += in(i, j) * (*r)[i];
            (*c)[j] = (s != 0.0) ? 1.0 / s : 0.0;
        }
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < n; ++j) s += in(i, j) * (*c)[j];
            worst = std::max(worst, std::fabs((*r)[i] * s - 1.0));
        }
        if (worst < tolerance) {
            converged = true;
            break;
        }
        last_error = worst;
    }
    if (!converged)
        throw InputError("sinkhorn_scaling: did not converge to a doubly stochastic matrix in " +
                         std::to_string(max_iterations) + " sweeps (margin error " +
                         std::to_string(last_error) + " against a tolerance of " +
                         std::to_string(tolerance) +
                         "). The usual cause is a matrix without total support.");
    *B = Matrix<double>(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) (*B)(i, j) = (*r)[i] * in(i, j) * (*c)[j];
}

/**
 * The Sinkhorn heuristic. NO error bound in either direction -- see the header.
 *
 * A zero entry used to be nudged to 1e-15 before scaling. That is not
 * invertible -- it changes the permanent by n!*eps -- so a non-positive entry
 * is now REFUSED by `require_full_support` instead.
 */
inline double perm_heur(const Matrix<double>& m, double tolerance = 1e-10,
                        std::size_t max_iterations = 1000) {
    approxdetail::require_nonnegative_square(m, "perm_heur");
    const std::size_t n = m.rows();
    if (n > 0) require_full_support(m, "perm_heur");

    Matrix<double> work = m;

    Matrix<double> B;
    std::vector<double> r, c;
    sinkhorn_scaling(work, &B, &r, &c, tolerance, max_iterations);

    std::vector<double> rowsum(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) rowsum[i] += B(i, j);

    const double nd = static_cast<double>(n);
    double rowprod = 1.0, logsum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        rowprod *= rowsum[i];
        logsum += std::log(rowsum[i]);
    }
    const double fact = approxdetail::factorial_d(n);
    const double p_meanfield = fact * (rowprod / std::pow(nd, nd));
    const double cap = std::exp(logsum / nd);
    const double p_gurvits = fact * std::pow(cap / nd, nd);
    const double p_est = 0.5 * (p_meanfield + p_gurvits);

    double scale = 1.0;
    for (std::size_t i = 0; i < n; ++i) scale *= 1.0 / r[i];
    for (std::size_t j = 0; j < n; ++j) scale *= 1.0 / c[j];
    return p_est * scale;
}

/**
 * The Bethe permanent, by sum-product message passing.
 *
 * A LOWER BOUND of the permanent for a nonnegative matrix. Two message
 * families -- right-going `r` and left-going `l` -- are iterated to a fixed
 * point on the SQUARE ROOT of the matrix, and the Bethe free energy there is
 * exponentiated.
 *
 * THE DENOMINATOR EXCLUDES THE DIAGONAL, NOT THE ENTRY ITSELF. A textbook
 * sum-product message from (i,j) omits column j of row i; this scheme omits
 * the DIAGONAL element of the row instead, so one denominator serves the whole
 * row and the numerator carries `s(i,j)`. The two agree on a symmetric matrix
 * and disagree otherwise -- measured on a 4x4 dense instance, 386.3 for the
 * textbook form against the reference's 325.7, both below the exact 1092 and
 * so both bounds, but only one of them the reference's. Transcribed as written.
 *
 * @param epsilon       squared message change below which the iteration stops
 * @param max_iteration cap on the message passing
 */
inline double perm_bethe(const Matrix<double>& m, double epsilon = 0.001,
                         std::size_t max_iteration = 200000) {
    approxdetail::require_nonnegative_square(m, "perm_bethe");
    const std::size_t n = m.rows();
    if (n > 0) require_full_support(m, "perm_bethe");
    // kMin guards the LOGARITHMS of message products below against underflow.
    // It is deliberately NOT applied to the input: flooring the input is what
    // fabricates a permanent of n!*eps where the truth is zero.
    const double kMin = 2.220446049250314e-16;
    (void)kMin;

    Matrix<double> s(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) s(i, j) = std::sqrt(m(i, j));

    // One update: right-going from the current left-going, then left-going
    // from the right-going just produced. The second half reads the FRESH r,
    // which is why the two cannot be swapped.
    auto update = [&s, n](const Matrix<double>& l, Matrix<double>* r1, Matrix<double>* l1) {
        *r1 = Matrix<double>(n, n, 0.0);
        *l1 = Matrix<double>(n, n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double d = 0.0;
            for (std::size_t j = 0; j < n; ++j) d += s(i, j) * l(i, j);
            d -= s(i, i) * l(i, i);  // the DIAGONAL, not the (i,j) term
            for (std::size_t j = 0; j < n; ++j) (*r1)(i, j) = (d != 0.0) ? s(i, j) / d : 0.0;
        }
        for (std::size_t j = 0; j < n; ++j) {
            double d = 0.0;
            for (std::size_t i = 0; i < n; ++i) d += s(i, j) * (*r1)(i, j);
            d -= s(j, j) * (*r1)(j, j);
            for (std::size_t i = 0; i < n; ++i) (*l1)(i, j) = (d != 0.0) ? s(i, j) / d : 0.0;
        }
    };

    Matrix<double> r_past(n, n, 1.0), l_past(n, n, 1.0), r, l;
    update(l_past, &r, &l);
    for (std::size_t it = 0; it < max_iteration; ++it) {
        double change = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const double a = r_past(i, j) - r(i, j), b = l_past(i, j) - l(i, j);
                change += a * a + b * b;
            }
        if (change <= epsilon) break;
        r_past = r;
        l_past = l;
        update(l_past, &r, &l);
    }

    // The Bethe free energy at the fixed point, in logs so the products do not
    // overflow: row sums of s*l, column sums of s*r, less the edge term r*l+1.
    double logv = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double t = 0.0;
        for (std::size_t j = 0; j < n; ++j) t += s(i, j) * l(i, j);
        logv += std::log(std::max(t, kMin));
    }
    for (std::size_t j = 0; j < n; ++j) {
        double t = 0.0;
        for (std::size_t i = 0; i < n; ++i) t += s(i, j) * r(i, j);
        logv += std::log(std::max(t, kMin));
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            logv -= std::log(std::max(r(i, j) * l(i, j) + 1.0, kMin));

    const double out = std::exp(logv);
    // A non-finite free energy is not a bound; the reference reports zero.
    return std::isfinite(out) ? out : 0.0;
}

/** Outcome of the saddle-point expansion: the estimate and what produced it. */
struct PermSpmResult {
    double value = 1.0;          ///< the approximate permanent
    double log_value = 0.0;      ///< its logarithm, correct even when `value` overflows
    double log_capacity = 0.0;   ///< log Gurvits capacity, an upper bound on the log permanent
    std::vector<double> xi;      ///< saddle point, unit geometric mean, 0 on a dropped column
};

/**
 * Saddle-point (SPM) approximation of the permanent. THE HOMOGENEOUS CACHE_SPM.
 *
 * Approximates the permanent of the matrix built from the n x h matrix `a` by
 * repeating column l exactly `mult[l]` times, so sum(mult) must equal n; an
 * empty `mult` means all ones, which requires a square matrix.
 *
 * cache_spm and this routine evaluate the SAME Cauchy integral by Laplace's
 * method and differ only in the generating function whose coefficient they
 * extract:
 *
 *   cache_spm  E(m) = prod_l m_l! [prod_l z_l^m_l] prod_k (1 + sum_l g_kl z_l)
 *   perm_spm   P    = prod_l m_l! [prod_l z_l^m_l] prod_k (    sum_l a_kl z_l)
 *
 * The cache factor carries a "+1" because an item may stay out of the cache, so
 * what it extracts is a RECTANGULAR permanent over n items and sum(m) < n
 * slots. Dropping the "+1" forces every row to be matched, which is exactly the
 * permanent and requires sum(m) == n -- the one case cache_spm cannot serve,
 * since at n == sum(m) its multipliers diverge and it falls back on cache_erec.
 * Here the integrand is homogeneous and the saddle point is interior in the
 * h-1 directions that survive.
 *
 * METHOD. With z_l = xi_l exp(i th_l) the saddle point in xi solves
 *
 *   sum_k a_kl xi_l / (sum_j a_kj xi_j) = m_l,    l = 1..h,
 *
 * so p_kl = a_kl xi_l / s_k with s = a*xi is the diagonal scaling of `a` to row
 * sums 1 and column sums m (Sinkhorn; doubly stochastic when m is all ones).
 * There phi = sum_k log s_k - sum_l m_l log xi_l is the log Gurvits capacity.
 * The Gaussian correction uses H = diag(m) - p'p, a weighted graph Laplacian on
 * the columns: H*ones = 0, which is the invariance of the integrand under
 * th -> th + c*ones that homogeneity creates. That direction is a full period
 * rather than a Gaussian, so it contributes 2*pi and leaves an (h-1)
 * dimensional Laplace integral; any principal (h-1) submatrix serves, because
 * all cofactors of a Laplacian are equal. The estimate is
 *
 *   log P = sum_l log(m_l!) - (h-1)/2 log(2 pi) + phi - 1/2 log det(H_red).
 *
 * ACCURACY, AND WHAT IT IS NOT. Exact for h == 1, where the permanent is
 * n! prod_k a(k,0). It is a genuine asymptotic expansion as min(m) grows with h
 * fixed, the ratio to the exact permanent falling from 1.11 at m = (2,2,2) to
 * 1.02 at m = (3,3). At m = ones the dimension of the integral grows with the
 * expansion parameter and the leading term keeps a systematic BIAS: on the n x n
 * matrix of ones it returns (2 pi)^(-(n-1)/2) n^(n+1/2) against the exact n!, a
 * ratio tending to (e/sqrt(2 pi))^n = 1.084^n, and random positive matrices
 * track that closely (1.31 at n = 4, 1.87 at n = 8). So at m = ones it
 * OVERESTIMATES, with a spread across matrices far tighter than the bias
 * itself, and it is NOT a bound in either direction; `perm_bethe` is the
 * routine to use when a bound is needed.
 *
 * REQUIRES A STRICTLY POSITIVE MATRIX, for the reasons the header states: the
 * scaling is what needs it, and flooring a zero is not invertible. Positivity
 * also makes the column graph complete, hence H_red positive definite.
 */
inline PermSpmResult perm_spm_expand(const Matrix<double>& a,
                                     const std::vector<std::size_t>& mult,
                                     double tolerance = 1e-11,
                                     std::size_t max_iterations = 10000) {
    PermSpmResult out;
    const std::size_t n = a.rows();
    const std::size_t h = a.cols();
    if (n == 0 || h == 0) return out;   // the permanent of the empty matrix is 1

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < h; ++j)
            if (a(i, j) < 0.0)
                throw InputError("perm_spm: the matrix must be non-negative; entry (" +
                                 std::to_string(i) + ", " + std::to_string(j) + ") is " +
                                 std::to_string(a(i, j)));

    std::vector<std::size_t> m = mult;
    if (m.empty()) {
        if (h != n)
            throw InputError("perm_spm: without column multiplicities the matrix must be square;"
                             " it is " + std::to_string(n) + "x" + std::to_string(h));
        m.assign(n, 1);
    }
    if (m.size() != h)
        throw InputError("perm_spm: the multiplicity vector has " + std::to_string(m.size()) +
                         " entries against " + std::to_string(h) + " columns");
    std::size_t total = 0;
    for (std::size_t j = 0; j < h; ++j) total += m[j];
    if (total != n)
        throw InputError("perm_spm: the column multiplicities must sum to the number of rows"
                         " (sum(m) = " + std::to_string(total) + " against " + std::to_string(n) +
                         " rows). The integrand is homogeneous of degree " + std::to_string(n) +
                         ", so every other coefficient of it is exactly zero");
    require_full_support(a, "perm_spm");

    // A column repeated zero times leaves the permanent unchanged, and its xi is a
    // boundary of the Laplace integral rather than a direction of it, so it must leave
    // the expansion. Dropping it is exact: setting z_l = 0 removes the column, and
    // prod_l m_l! is unchanged because 0! = 1.
    std::vector<std::size_t> keep;
    for (std::size_t j = 0; j < h; ++j)
        if (m[j] > 0) keep.push_back(j);
    const std::size_t hk = keep.size();
    std::vector<double> mk(hk);
    for (std::size_t l = 0; l < hk; ++l) mk[l] = static_cast<double>(m[keep[l]]);

    std::vector<double> ak(n * hk);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t l = 0; l < hk; ++l) ak[i * hk + l] = a(i, keep[l]);

    // Saddle point: scale to row sums 1 and column sums mk. The row sums are 1 by
    // construction of p, so only the column sums are iterated on.
    std::vector<double> xik(hk, 1.0), s(n, 0.0), colsum(hk, 0.0);
    bool converged = false;
    double margin = std::numeric_limits<double>::infinity();
    for (std::size_t it = 0; it < max_iterations && !converged; ++it) {
        for (std::size_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (std::size_t l = 0; l < hk; ++l) acc += ak[i * hk + l] * xik[l];
            s[i] = acc;
        }
        margin = 0.0;
        for (std::size_t l = 0; l < hk; ++l) {
            double acc = 0.0;
            for (std::size_t i = 0; i < n; ++i) acc += ak[i * hk + l] / s[i];
            colsum[l] = xik[l] * acc;
            margin = std::max(margin, std::fabs(colsum[l] - mk[l]));
        }
        if (margin < tolerance) {
            converged = true;
            break;
        }
        double logmean = 0.0;
        for (std::size_t l = 0; l < hk; ++l) {
            xik[l] *= mk[l] / colsum[l];
            logmean += std::log(xik[l]);
        }
        logmean /= static_cast<double>(hk);
        const double scale = std::exp(logmean);
        for (std::size_t l = 0; l < hk; ++l) xik[l] /= scale;   // the saddle is a ray; pin its scale
    }
    if (!converged)
        throw InputError("perm_spm: the scaling to row sums 1 and column sums m did not converge"
                         " in " + std::to_string(max_iterations) + " sweeps (margin error " +
                         std::to_string(margin) + " against a tolerance of " +
                         std::to_string(tolerance) + "). The expansion assumes the saddle point,"
                         " so no value is returned. The usual cause is a matrix without total"
                         " support.");

    out.xi.assign(h, 0.0);
    for (std::size_t l = 0; l < hk; ++l) out.xi[keep[l]] = xik[l];

    for (std::size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t l = 0; l < hk; ++l) acc += ak[i * hk + l] * xik[l];
        s[i] = acc;
    }
    std::vector<double> p(n * hk);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t l = 0; l < hk; ++l) p[i * hk + l] = ak[i * hk + l] * xik[l] / s[i];

    double log_capacity = 0.0;
    for (std::size_t i = 0; i < n; ++i) log_capacity += std::log(s[i]);
    for (std::size_t l = 0; l < hk; ++l) log_capacity -= mk[l] * std::log(xik[l]);

    // H = diag(mk) - p'p is a Laplacian, so it is singular along ones and all of its
    // principal cofactors are equal; the last index is dropped only because one has to
    // be. Strict positivity makes the column graph complete, hence H_red positive
    // definite and Cholesky the right factor. At hk == 1 no direction survives the
    // homogeneity, and the determinant of the empty matrix is 1.
    double log_det = 0.0;
    if (hk > 1) {
        const std::size_t k = hk - 1;
        std::vector<double> hred(k * k);
        for (std::size_t l = 0; l < k; ++l)
            for (std::size_t j = 0; j < k; ++j) {
                double dot = 0.0;
                for (std::size_t i = 0; i < n; ++i) dot += p[i * hk + l] * p[i * hk + j];
                hred[l * k + j] = (l == j ? mk[l] : 0.0) - dot;
            }
        log_det = approxdetail::log_det_cholesky(hred, k, "perm_spm");
    }

    double log_fact = 0.0;
    for (std::size_t l = 0; l < hk; ++l) log_fact += std::lgamma(mk[l] + 1.0);

    out.log_capacity = log_capacity;
    out.log_value = log_fact - 0.5 * static_cast<double>(hk - 1) * std::log(2.0 * M_PI) +
                    log_capacity - 0.5 * log_det;
    out.value = std::exp(out.log_value);
    return out;
}

/** The saddle-point estimate of the permanent of a square strictly positive matrix. */
inline double perm_spm(const Matrix<double>& a, double tolerance = 1e-11,
                       std::size_t max_iterations = 10000) {
    return perm_spm_expand(a, std::vector<std::size_t>(), tolerance, max_iterations).value;
}

/** The saddle-point estimate with column l of `a` repeated `mult[l]` times. */
inline double perm_spm(const Matrix<double>& a, const std::vector<std::size_t>& mult,
                       double tolerance = 1e-11, std::size_t max_iterations = 10000) {
    return perm_spm_expand(a, mult, tolerance, max_iterations).value;
}

}  // namespace perm
}  // namespace line

#endif  // LINE_API_PERM_PERM_APPROX_H
