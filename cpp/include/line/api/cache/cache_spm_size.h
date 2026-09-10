/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_CACHE_SPM_SIZE_H
#define LINE_API_CACHE_CACHE_SPM_SIZE_H

/**
 * Ray (WKB) asymptotic expansion of the cost-capped cache normalizing constant.
 *
 * Templated port of matlab/src/api/cache/cache_spm_size.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_spm_size.java and
 * python/line_solver/api/cache/spm_size.py.
 *
 * Approximates what cache_erec(gamma, m, sigma, k) computes exactly, in the SAME
 * normalization, so the two are interchangeable. This is the item-size extension
 * of retrieval_rayint, which carries the size-free expansion; call that one when
 * there are no storage costs.
 *
 * Writing E = prod_j m_j! * H, the size-free recursion
 *
 *   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),  E(0,0)=1
 *
 * relaxes to H ~ exp(phi/eps) with n = y/eps, m_j = x_j/eps, whose eikonal
 * e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j} carries the ray constants
 * xi_j = e^{-phi_j}. With per-item storage costs sigma_i and per-list cost caps
 * k_j the recursion gains the cost coordinate,
 *
 *   E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j),
 *
 * so the shift 1_j becomes e_j(y) = (1_j, s(y) 1_j) in the enlarged space
 * X = (x,kappa) and the eikonal picks up the size tilt
 *
 *   e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_{x_j} - s(y) phi_{kappa_j}},
 *
 * with the second family of ray constants zeta_j = e^{-phi_{kappa_j}}. The rays
 * integrate to the discrete saddle point of the product generating function
 *
 *   sum_{m,k} H(m,k) prod_j z_j^{m_j} w_j^{k_j}
 *       = prod_i ( 1 + sum_j gamma_ij z_j w_j^{sigma_i} ),
 *
 * namely, with D_i = 1 + sum_j gamma_ij xi_j zeta_j^{sigma_i} and
 * Psi = sum_i log D_i,
 *
 *   m_j = sum_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 *   k_j = sum_i sigma_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 *   log H(m,k) ~ Psi - sum_j m_j log xi_j - sum_j k_j log zeta_j
 *                - (d/2) log(2 pi) - (1/2) log det grad^2 Psi,
 *
 * where d is the number of saddle coordinates and, with
 * pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i and
 * Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il,
 *
 *   grad^2 Psi = sum_i [1; sigma_i] [1; sigma_i]' (x) Q^i .
 *
 * Setting zeta_j = 1 recovers the size-free expansion exactly.
 *
 * CAPS ARE CUMULATIVE. cache_erec sums over the states of cost AT MOST k_j, so
 * this function does the same by default (CacheCostMode::AtMost). The shadow
 * price eta_j = log zeta_j <= 0 obeys complementary slackness: a list whose
 * unconstrained mean cost already meets its cap is SLACK, keeps zeta_j = 1, and
 * drops out of the saddle, which then degenerates continuously to the size-free
 * expansion; a list whose cap BINDS sits at eta_j < 0, and the states below the
 * boundary decay geometrically with ratio zeta_j, contributing the amplitude
 * factor 1/(1-zeta_j). CacheCostMode::Exact gives instead the constant resolving
 * the cost exactly at k_j, the raw Laplace formula above with no such factor.
 *
 * SIZE DIVERSITY IS REQUIRED. The Hessian integrand
 * [1;sigma_i][1;sigma_i]' (x) Q^i has rank h, not 2h, so grad^2 Psi is
 * nonsingular only if the sizes actually vary. This is not an artefact: with a
 * single item size the cost of list j is sigma*m_j identically and the cap
 * carries no information. That case is detected and answered exactly rather than
 * passed to a singular saddle. If the sizes share a common divisor the cost
 * lives on a sublattice; the sizes and caps are divided through by their gcd,
 * which is an exact reduction and removes the corresponding lattice factor.
 *
 * OCCUPANCY. The result's pij is the saddle occupancy
 * pi_il = gamma_il xi_l zeta_l^{sigma_i} / D_i and k_mean its per-list cost.
 * These are EXACT-COST quantities: the saddle conditions are sum_i pi_ij = m_j
 * and sum_i sigma_i pi_ij = k_j, so k_mean equals the cap exactly on every
 * binding list. Under cumulative caps the true mean cost is strictly below the
 * cap; use cache_cost and cache_prob_erec for that.
 *
 * ACCURACY. The expansion is O(1/n) at fixed occupancy. With a well-separated
 * cap the observed error in log E is around 1e-2 at n = 200 and halves at each
 * doubling of n. It degrades as a binding zeta_j approaches 1, i.e. in the
 * transition between the binding and slack regimes, where the geometric
 * resummation 1/(1-zeta_j) is no longer sharp; the result's zeta and binding
 * report where the saddle sits, and relerr_est is the size-free O(1/n) baseline
 * that does NOT cover that transition.
 *
 * ARITHMETIC: the expansion is a Laplace approximation built out of logs, exps
 * and a square root, so it is meaningless at exact arithmetic and is gated on
 * has_transcendental. It is also an APPROXIMATION whatever the arithmetic --
 * widening the type sharpens the saddle solve, never the O(1/n) model error.
 * Use cache_erec when the exact constant is wanted.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Whether the caps bound the cost from above (matching cache_erec) or resolve it exactly. */
enum class CacheCostMode { AtMost, Exact };

/** Outcome of the expansion. */
template <class T>
struct CacheSpmSizeResult {
    /** Normalizing constant, same normalization as cache_erec (may overflow; use log_e). */
    T e;
    /** Natural logarithm of e, safe for large n. */
    T log_e;
    /** Saddle point xi_j, one entry per list (0 for a list of zero capacity). */
    std::vector<T> xi;
    /** Cost tilt zeta_j on the original size lattice (1 for a slack or absent list). */
    std::vector<T> zeta;
    /** Whether each list's cost cap binds. */
    std::vector<bool> binding;
    /** Occupancy pi, n x (h+1), column 0 the miss probability. */
    Matrix<T> pij;
    /** Mean storage cost held by each list. */
    std::vector<T> k_mean;
    /** The exponent Psi - m.log xi - k.log zeta. */
    T phi;
    /** log det of the Hessian in (log xi, log zeta), restricted to the free coordinates. */
    T logdet_sigma;
    /** gcd of the item sizes, divided out as an exact lattice reduction. */
    int span = 1;
    /** "spm-size", "spm", "uniform-size", "lattice" or "boundary". */
    std::string method;
    /** Size-free error baseline 0.14*(1/min_j m_j + 1/(n - sum_j m_j)); see ACCURACY. */
    double relerr_est = 0.0;
    /** Newton iterations used. */
    std::size_t iterations = 0;
};

namespace detail {

/**
 * pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i, D_i = 1 + sum_j (that numerator).
 * `d_out`, when non-null, receives the D_i so the caller need not recompute them.
 */
template <class T>
Matrix<T> spm_size_occupancy(const Matrix<T>& g, const std::vector<T>& sg,
                                const std::vector<T>& th, const std::vector<T>& et,
                                std::vector<T>* d_out = nullptr) {
    using std::exp;
    const std::size_t n = g.rows();
    const std::size_t h = th.size();
    const T one = num_traits<T>::from_int(1);
    Matrix<T> p(n, h, num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < n; ++a) {
        T s = one;
        for (std::size_t b = 0; b < h; ++b) {
            p(a, b) = g(a, b) * exp(th[b] + sg[a] * et[b]);
            s += p(a, b);
        }
        if (d_out != nullptr) (*d_out)[a] = s;
        for (std::size_t b = 0; b < h; ++b) p(a, b) = p(a, b) / s;
    }
    return p;
}

/**
 * grad^2 Psi in (theta, eta), restricted to the free eta coordinates `ix`. Coordinate
 * u < h is theta_u with weight 1; coordinate h+b is eta_{ix[b]} with weight sigma_i.
 * The block is sum_i w_u(i) w_v(i) Q^i_{ju,jv} with Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il.
 */
template <class T>
Matrix<T> spm_size_hessian(const Matrix<T>& p, const std::vector<T>& sg,
                              const std::vector<std::size_t>& ix) {
    const std::size_t n = p.rows();
    const std::size_t h = p.cols();
    const std::size_t nb = ix.size();
    const std::size_t dim = h + nb;
    std::vector<std::size_t> coord(dim);
    for (std::size_t u = 0; u < h; ++u) coord[u] = u;
    for (std::size_t b = 0; b < nb; ++b) coord[h + b] = ix[b];
    Matrix<T> hess(dim, dim, num_traits<T>::from_int(0));
    std::vector<T> w(dim, num_traits<T>::from_int(1));
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = 0; b < nb; ++b) w[h + b] = sg[a];
        for (std::size_t u = 0; u < dim; ++u) {
            const std::size_t ju = coord[u];
            for (std::size_t v = 0; v < dim; ++v) {
                const std::size_t jv = coord[v];
                const T q = (ju == jv ? p(a, ju) : zero) - p(a, ju) * p(a, jv);
                hess(u, v) += w[u] * w[v] * q;
            }
        }
    }
    return hess;
}

/** Gaussian elimination with partial pivoting; the order is at most 2h, so tiny. */
template <class T>
std::vector<T> spm_size_solve(const Matrix<T>& a, const std::vector<T>& b) {
    using std::abs;
    const std::size_t dim = b.size();
    Matrix<T> m(dim, dim + 1);
    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) m(i, j) = a(i, j);
        m(i, dim) = b[i];
    }
    for (std::size_t c = 0; c < dim; ++c) {
        std::size_t piv = c;
        for (std::size_t i = c + 1; i < dim; ++i)
            if (abs(m(i, c)) > abs(m(piv, c))) piv = i;
        if (num_traits<T>::to_double(abs(m(piv, c))) == 0.0)
            throw NumericError("cache_spm_size: the saddle-point Newton step is not finite. With "
                               "item sizes this is the rank-h degeneracy of the size-tilted Hessian: "
                               "the sizes must genuinely vary for the cost coordinate to carry "
                               "information");
        if (piv != c)
            for (std::size_t j = 0; j <= dim; ++j) {
                const T tmp = m(c, j);
                m(c, j) = m(piv, j);
                m(piv, j) = tmp;
            }
        for (std::size_t i = c + 1; i < dim; ++i) {
            const T f = m(i, c) / m(c, c);
            for (std::size_t j = c; j <= dim; ++j) m(i, j) -= f * m(c, j);
        }
    }
    std::vector<T> x(dim);
    for (std::size_t ii = dim; ii-- > 0;) {
        T s = m(ii, dim);
        for (std::size_t j = ii + 1; j < dim; ++j) s -= m(ii, j) * x[j];
        x[ii] = s / m(ii, ii);
    }
    return x;
}

/** log det via Cholesky; the Hessian of a strictly convex objective is positive definite. */
template <class T>
T spm_size_logdet(const Matrix<T>& a) {
    using std::log;
    using std::sqrt;
    const std::size_t dim = a.rows();
    const T half = num_traits<T>::from_rational(1, 2);
    const T two = num_traits<T>::from_int(2);
    Matrix<T> l(dim, dim, num_traits<T>::from_int(0));
    T ld = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T s = half * (a(i, j) + a(j, i));
            for (std::size_t c = 0; c < j; ++c) s -= l(i, c) * l(j, c);
            if (i == j) {
                if (num_traits<T>::to_double(s) <= 0.0)
                    throw NumericError("cache_spm_size: the saddle-point Hessian is not positive "
                                       "definite; the ray map is singular here. With item sizes this "
                                       "happens when the sizes do not vary over the items the cache "
                                       "can hold, in which case the cost cap carries no information");
                l(i, j) = sqrt(s);
                ld += two * log(l(i, j));
            } else {
                l(i, j) = s / l(j, j);
            }
        }
    }
    return ld;
}

/** Keeps the tilts within a factor e^2 per iteration. */
template <class T>
double spm_size_damp(const std::vector<T>& d) {
    double dmax = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i)
        dmax = std::max(dmax, std::abs(num_traits<T>::to_double(d[i])));
    double step = 1.0;
    while (step * dmax > 2.0) step /= 2.0;
    return step;
}

template <class T>
std::vector<T> spm_size_theta0(const Matrix<T>& g, const std::vector<T>& tgt) {
    using std::log;
    const std::size_t n = g.rows();
    const std::size_t h = tgt.size();
    double tsum = 0.0;
    for (std::size_t b = 0; b < h; ++b) tsum += num_traits<T>::to_double(tgt[b]);
    const double slack = std::max(1.0 - tsum / static_cast<double>(n), 1e-9);
    const T slackT = num_traits<T>::from_double(slack);
    const T tiny = num_traits<T>::from_double(1e-12);
    std::vector<T> th(h);
    for (std::size_t b = 0; b < h; ++b) {
        T gb = num_traits<T>::from_int(0);
        for (std::size_t a = 0; a < n; ++a) gb += g(a, b);
        T den = gb * slackT;
        if (num_traits<T>::to_double(den) < 1e-12) den = tiny;
        T num = tgt[b];
        if (num_traits<T>::to_double(num) < 1e-12) num = tiny;
        th[b] = log(num / den);
    }
    return th;
}

/** The convex dual f = sum_i log D_i - m.theta - k.eta, over the binding eta only. */
template <class T>
T spm_size_obj(const Matrix<T>& g, const std::vector<T>& sg, const std::vector<T>& th,
                  const std::vector<T>& et, const std::vector<T>& tgtm, const std::vector<T>& tgtk,
                  const std::vector<bool>& bind) {
    using std::log;
    std::vector<T> d(g.rows());
    spm_size_occupancy(g, sg, th, et, &d);
    T f = num_traits<T>::from_int(0);
    for (std::size_t a = 0; a < g.rows(); ++a) f += log(d[a]);
    for (std::size_t b = 0; b < th.size(); ++b) {
        f -= tgtm[b] * th[b];
        if (bind[b]) f -= tgtk[b] * et[b];
    }
    return f;
}

/** Size-free saddle: Newton on theta = log xi for sum_i gamma_ij xi_j / D_i = m_j. */
template <class T>
std::vector<T> spm_size_saddle_free(const Matrix<T>& g, const std::vector<T>& tgt,
                                       std::size_t& iterations) {
    const std::size_t n = g.rows();
    const std::size_t h = tgt.size();
    const std::vector<T> zeros_n(n, num_traits<T>::from_int(0));
    const std::vector<T> zeros_h(h, num_traits<T>::from_int(0));
    const std::vector<std::size_t> noix;
    std::vector<T> th = spm_size_theta0(g, tgt);
    double tmax = 1.0;
    for (std::size_t b = 0; b < h; ++b)
        tmax = std::max(tmax, std::abs(num_traits<T>::to_double(tgt[b])));
    std::size_t it = 0;
    for (it = 1; it <= 200; ++it) {
        const Matrix<T> p = spm_size_occupancy(g, zeros_n, th, zeros_h);
        std::vector<T> grad(h);
        double gmax = 0.0;
        for (std::size_t b = 0; b < h; ++b) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t a = 0; a < n; ++a) s += p(a, b);
            grad[b] = s - tgt[b];
            gmax = std::max(gmax, std::abs(num_traits<T>::to_double(grad[b])));
        }
        if (gmax <= 1e-12 * tmax) break;
        std::vector<T> d = spm_size_solve(spm_size_hessian(p, zeros_n, noix), grad);
        for (std::size_t b = 0; b < h; ++b) d[b] = -d[b];
        const T step = num_traits<T>::from_double(spm_size_damp(d));
        for (std::size_t b = 0; b < h; ++b) th[b] += step * d[b];
    }
    iterations = it;
    return th;
}

/**
 * Cost-constrained saddle. Minimises the convex dual
 * f(theta,eta) = sum_i log D_i - m.theta - k.eta over eta <= 0 when the caps are
 * cumulative, so that complementary slackness selects the binding lists; over all
 * of R^{2h} when the cost is resolved exactly.
 */
template <class T>
void spm_size_saddle(const Matrix<T>& g, const std::vector<T>& tgtm, const std::vector<T>& sg,
                        const std::vector<T>& tgtk, bool cumulative, std::vector<T>& th,
                        std::vector<T>& et, std::vector<bool>& bind, std::size_t& iterations) {
    const std::size_t n = g.rows();
    const std::size_t h = tgtm.size();
    const T zero = num_traits<T>::from_int(0);
    th = spm_size_theta0(g, tgtm);
    et.assign(h, zero);
    bind.assign(h, true);
    double tol = 1.0;
    for (std::size_t b = 0; b < h; ++b) {
        tol = std::max(tol, std::abs(num_traits<T>::to_double(tgtm[b])));
        tol = std::max(tol, std::abs(num_traits<T>::to_double(tgtk[b])));
    }
    tol *= 1e-12;
    std::vector<T> thn(h), etn(h);
    std::size_t it = 0;
    for (it = 1; it <= 200; ++it) {
        const Matrix<T> p = spm_size_occupancy(g, sg, th, et);
        std::vector<T> gth(h), get(h);
        for (std::size_t b = 0; b < h; ++b) {
            T s1 = zero;
            T s2 = zero;
            for (std::size_t a = 0; a < n; ++a) {
                s1 += p(a, b);
                s2 += sg[a] * p(a, b);
            }
            gth[b] = s1 - tgtm[b];
            get[b] = s2 - tgtk[b];
        }
        if (cumulative)
            for (std::size_t b = 0; b < h; ++b)   // at eta_j = 0 the cap binds when the cost exceeds it
                bind[b] = num_traits<T>::to_double(et[b]) < 0.0 ||
                          num_traits<T>::to_double(get[b]) > 0.0;
        std::vector<std::size_t> ix;
        for (std::size_t b = 0; b < h; ++b)
            if (bind[b]) ix.push_back(b);
        const std::size_t nb = ix.size();
        std::vector<T> grad(h + nb);
        double gmax = 0.0;
        for (std::size_t b = 0; b < h; ++b) {
            grad[b] = gth[b];
            gmax = std::max(gmax, std::abs(num_traits<T>::to_double(grad[b])));
        }
        for (std::size_t b = 0; b < nb; ++b) {
            grad[h + b] = get[ix[b]];
            gmax = std::max(gmax, std::abs(num_traits<T>::to_double(grad[h + b])));
        }
        if (gmax <= tol) break;
        std::vector<T> d = spm_size_solve(spm_size_hessian(p, sg, ix), grad);
        for (std::size_t b = 0; b < d.size(); ++b) d[b] = -d[b];
        double step = spm_size_damp(d);
        const T fcur = spm_size_obj(g, sg, th, et, tgtm, tgtk, bind);
        // Backtrack until the dual decreases. The slack is essential, not cosmetic:
        // Newton reaches the floating-point floor of f in a handful of steps, and a
        // strict test then rejects every step and halves to zero without converging.
        const double ftol = 1e-12 * (1.0 + std::abs(num_traits<T>::to_double(fcur)));
        for (int ls = 0; ls < 40; ++ls) {
            const T stepT = num_traits<T>::from_double(step);
            for (std::size_t b = 0; b < h; ++b) {
                thn[b] = th[b] + stepT * d[b];
                etn[b] = et[b];
            }
            for (std::size_t b = 0; b < nb; ++b) etn[ix[b]] = et[ix[b]] + stepT * d[h + b];
            if (cumulative)
                for (std::size_t b = 0; b < h; ++b)
                    if (num_traits<T>::to_double(etn[b]) > 0.0) etn[b] = zero;
            const T fn = spm_size_obj(g, sg, thn, etn, tgtm, tgtk, bind);
            if (num_traits<T>::to_double(fn) <= num_traits<T>::to_double(fcur) + ftol) break;
            step /= 2.0;
        }
        double moved = 0.0;
        for (std::size_t b = 0; b < h; ++b) {
            moved = std::max(moved, std::abs(num_traits<T>::to_double(thn[b] - th[b])));
            moved = std::max(moved, std::abs(num_traits<T>::to_double(etn[b] - et[b])));
        }
        th = thn;
        et = etn;
        if (moved <= 1e-13) break;   // the iterate can no longer move: at the floor
    }
    iterations = it;
    if (cumulative)
        for (std::size_t b = 0; b < h; ++b) bind[b] = num_traits<T>::to_double(et[b]) < 0.0;
}

/** log(x!) for a non-negative integer x, by summing logs; the counts here are small. */
template <class T>
T spm_size_logfact(int x) {
    using std::log;
    T s = num_traits<T>::from_int(0);
    for (int i = 2; i <= x; ++i) s += log(num_traits<T>::from_int(i));
    return s;
}

inline long spm_size_gcd(long a, long b) {
    long x = a < 0 ? -a : a;
    long y = b < 0 ? -b : b;
    while (y != 0) {
        const long r = x % y;
        x = y;
        y = r;
    }
    return x;
}

}  // namespace detail

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities, non-negative
 * @param sigma (n) per-item storage costs, positive
 * @param k     (h) per-list storage cost caps
 * @param mode  AtMost (matches cache_erec) or Exact
 */
template <class T>
CacheSpmSizeResult<T> cache_spm_size(const Matrix<T>& gamma, const std::vector<int>& m,
                                           const std::vector<int>& sigma, const std::vector<int>& k,
                                           CacheCostMode mode = CacheCostMode::AtMost) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_spm_size is a Laplace approximation built out of logs, exps and a "
                  "square root: it is meaningless at exact arithmetic, and widening the type "
                  "sharpens the saddle solve but never the O(1/n) model error. Use cache_erec "
                  "for the exact constant");
    using std::exp;
    using std::log;

    if (gamma.rows() == 0 || gamma.cols() == 0)
        throw InputError("cache_spm_size: gamma must be a non-empty n x h matrix");
    const std::size_t n0 = gamma.rows();
    const std::size_t h0 = gamma.cols();
    if (m.size() != h0)
        throw InputError("cache_spm_size: the capacity vector must have one entry per cache list");
    for (std::size_t j = 0; j < h0; ++j)
        if (m[j] < 0) throw InputError("cache_spm_size: list capacities must be non-negative");
    if (sigma.empty() || k.empty())
        throw InputError("cache_spm_size: the item sizes and the cost caps are both required; use "
                         "retrieval_rayint for the size-free expansion");
    if (sigma.size() != n0)
        throw InputError("cache_spm_size: the item size vector must have one entry per item");
    if (k.size() != h0)
        throw InputError("cache_spm_size: the cost cap vector must have one entry per cache list");
    for (std::size_t i = 0; i < n0; ++i)
        if (sigma[i] <= 0) throw InputError("cache_spm_size: item sizes must be positive integers");
    const bool exact_mode = (mode == CacheCostMode::Exact);
    bool capped = true;   // cleared below when a single item size makes the cap uninformative

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    CacheSpmSizeResult<T> out;
    out.e = zero;
    out.log_e = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
    out.xi.assign(h0, zero);
    out.zeta.assign(h0, one);
    out.binding.assign(h0, false);
    out.k_mean.assign(h0, zero);
    out.phi = zero;
    out.logdet_sigma = zero;
    out.pij = Matrix<T>(n0, h0 + 1, zero);
    for (std::size_t i = 0; i < n0; ++i) out.pij(i, 0) = one;

    // --- boundaries, matching cache_erec ---
    long msum = 0;
    for (std::size_t j = 0; j < h0; ++j) msum += m[j];
    bool negcap = false;
    for (std::size_t j = 0; j < h0; ++j)
        if (k[j] < 0) negcap = true;
    if (msum > static_cast<long>(n0) || negcap) {
        out.method = "boundary";
        return out;
    }
    if (msum == 0) {
        bool poscap = false;
        for (std::size_t j = 0; j < h0; ++j)
            if (k[j] > 0) poscap = true;
        out.method = "boundary";
        if (!(exact_mode && poscap)) {
            out.e = one;
            out.log_e = zero;
        }
        return out;
    }

    // --- items that can never be cached and lists of zero capacity drop out ---
    std::vector<std::size_t> alive;
    for (std::size_t i = 0; i < n0; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < h0; ++j) s += gamma(i, j);
        if (num_traits<T>::to_double(s) > 0.0) alive.push_back(i);
    }
    std::vector<std::size_t> keep;
    for (std::size_t j = 0; j < h0; ++j)
        if (m[j] > 0) keep.push_back(j);
    const std::size_t n = alive.size();
    const std::size_t hk = keep.size();
    long mksum = 0;
    for (std::size_t b = 0; b < hk; ++b) mksum += m[keep[b]];
    if (mksum > static_cast<long>(n)) {
        out.method = "boundary";
        return out;
    }
    if (mksum == static_cast<long>(n))
        throw InputError("cache_spm_size: the expansion requires sum(m) < n; at sum(m) = n the "
                         "saddle point escapes to infinity, use cache_erec for a full cache");

    Matrix<T> g(n, hk);
    std::vector<long> sgi(n);
    for (std::size_t a = 0; a < n; ++a) {
        sgi[a] = sigma[alive[a]];
        for (std::size_t b = 0; b < hk; ++b) g(a, b) = gamma(alive[a], keep[b]);
    }
    std::vector<long> kki(hk);
    std::vector<long> mki(hk);
    for (std::size_t b = 0; b < hk; ++b) {
        mki[b] = m[keep[b]];
        kki[b] = k[keep[b]];
    }

    // --- exact reductions on the cost lattice ---
    long span = 0;
    for (std::size_t a = 0; a < n; ++a) span = detail::spm_size_gcd(span, sgi[a]);
    out.span = static_cast<int>(span);
    if (exact_mode)
        for (std::size_t b = 0; b < hk; ++b)
            if (kki[b] % span != 0) {   // unreachable off the sublattice
                out.method = "lattice";
                return out;
            }
    for (std::size_t a = 0; a < n; ++a) sgi[a] /= span;
    for (std::size_t b = 0; b < hk; ++b) kki[b] = kki[b] / span;
    // per-list feasibility: the m_j cheapest (dearest) reachable items bound the cost
    for (std::size_t b = 0; b < hk; ++b) {
        std::vector<long> srt;
        for (std::size_t a = 0; a < n; ++a)
            if (num_traits<T>::to_double(g(a, b)) > 0.0) srt.push_back(sgi[a]);
        if (static_cast<long>(srt.size()) < mki[b]) {
            out.method = "boundary";
            return out;
        }
        std::sort(srt.begin(), srt.end());
        long lo = 0;
        for (long a = 0; a < mki[b]; ++a) lo += srt[static_cast<std::size_t>(a)];
        if (lo > kki[b]) {
            out.method = "boundary";
            return out;
        }
        if (exact_mode) {
            long hi = 0;
            for (std::size_t a = srt.size() - static_cast<std::size_t>(mki[b]); a < srt.size(); ++a)
                hi += srt[a];
            if (hi < kki[b]) {
                out.method = "boundary";
                return out;
            }
        }
    }
    // a single item size makes the cost of list j equal to sigma*m_j identically,
    // so the cap carries no information and the 2h saddle is singular (rank h)
    bool uniform = true;
    for (std::size_t a = 1; a < n; ++a)
        if (sgi[a] != sgi[0]) uniform = false;
    if (uniform) {
        bool feasible = true;
        for (std::size_t b = 0; b < hk; ++b) {
            const long cost = sgi[0] * mki[b];
            if (exact_mode ? (cost != kki[b]) : (cost > kki[b])) feasible = false;
        }
        if (!feasible) {
            out.method = "uniform-size";
            return out;
        }
        capped = false;                 // fall through to the size-free expansion
        out.method = "uniform-size";
    }

    std::vector<T> sg(n);
    std::vector<T> mk(hk), kk(hk);
    for (std::size_t a = 0; a < n; ++a) sg[a] = num_traits<T>::from_int(sgi[a]);
    for (std::size_t b = 0; b < hk; ++b) {
        mk[b] = num_traits<T>::from_int(mki[b]);
        kk[b] = num_traits<T>::from_int(kki[b]);
    }

    // --- saddle point ---
    std::vector<T> th, et;
    std::vector<bool> bind;
    std::size_t iters = 0;
    std::vector<T> dvec(n);
    Matrix<T> p;
    T phi = zero;
    T logdet = zero;
    T log_h = zero;
    const T half = num_traits<T>::from_rational(1, 2);
    const T log2pi = num_traits<T>::from_double(std::log(2.0 * M_PI));
    if (capped) {
        detail::spm_size_saddle(g, mk, sg, kk, !exact_mode, th, et, bind, iters);
        p = detail::spm_size_occupancy(g, sg, th, et, &dvec);
        std::vector<std::size_t> ix;
        for (std::size_t b = 0; b < hk; ++b)
            if (bind[b]) ix.push_back(b);
        for (std::size_t a = 0; a < n; ++a) phi += log(dvec[a]);
        for (std::size_t b = 0; b < hk; ++b) phi -= mk[b] * th[b];
        for (std::size_t b = 0; b < ix.size(); ++b) phi -= kk[ix[b]] * et[ix[b]];
        logdet = detail::spm_size_logdet(detail::spm_size_hessian(p, sg, ix));
        const T dof = num_traits<T>::from_int(static_cast<long>(hk + ix.size()));
        log_h = phi - half * dof * log2pi - half * logdet;
        if (!exact_mode)
            for (std::size_t b = 0; b < ix.size(); ++b)   // geometric resummation below the cap
                log_h -= log(one - exp(et[ix[b]]));
        if (out.method.empty()) out.method = "spm-size";
    } else {
        const std::vector<T> zeros_n(n, zero);
        th = detail::spm_size_saddle_free(g, mk, iters);
        et.assign(hk, zero);
        bind.assign(hk, false);
        p = detail::spm_size_occupancy(g, zeros_n, th, et, &dvec);
        for (std::size_t a = 0; a < n; ++a) phi += log(dvec[a]);
        for (std::size_t b = 0; b < hk; ++b) phi -= mk[b] * th[b];
        const std::vector<std::size_t> noix;
        logdet = detail::spm_size_logdet(detail::spm_size_hessian(p, zeros_n, noix));
        const T dof = num_traits<T>::from_int(static_cast<long>(hk));
        log_h = phi - half * dof * log2pi - half * logdet;
        if (out.method.empty()) out.method = "spm";
    }

    T logfact = zero;
    for (std::size_t j = 0; j < h0; ++j) logfact += detail::spm_size_logfact<T>(m[j]);
    out.log_e = log_h + logfact;         // back to the cache_erec normalization
    out.e = exp(out.log_e);

    // --- ray quantities, reported on the original item and list indexing ---
    const T spanT = num_traits<T>::from_int(span);
    for (std::size_t b = 0; b < hk; ++b) {
        out.xi[keep[b]] = exp(th[b]);
        out.zeta[keep[b]] = exp(et[b] / spanT);
        out.binding[keep[b]] = bind[b];
    }
    for (std::size_t a = 0; a < n; ++a) {
        T hit = zero;
        for (std::size_t b = 0; b < hk; ++b) {
            out.pij(alive[a], 1 + keep[b]) = p(a, b);
            hit += p(a, b);
        }
        out.pij(alive[a], 0) = one - hit;
    }
    for (std::size_t j = 0; j < h0; ++j) {
        T s = zero;
        for (std::size_t i = 0; i < n0; ++i) s += num_traits<T>::from_int(sigma[i]) * out.pij(i, 1 + j);
        out.k_mean[j] = s;
    }
    out.phi = phi;
    out.logdet_sigma = logdet;
    out.iterations = iters;
    long mmin = mki[0];
    for (std::size_t b = 1; b < hk; ++b) mmin = std::min(mmin, mki[b]);
    out.relerr_est = 0.14 * (1.0 / static_cast<double>(mmin) +
                             1.0 / static_cast<double>(n - static_cast<std::size_t>(mksum)));
    return out;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_SPM_SIZE_H
