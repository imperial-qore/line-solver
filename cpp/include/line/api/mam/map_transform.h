/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_TRANSFORM_H
#define LINE_API_MAM_MAP_TRANSFORM_H

/**
 * MAP constructors and structural transformations.
 *
 * Templated port of the kpctoolbox MAP algebra that the QBD solvers consume
 * (matlab/lib/kpctoolbox/map/map_normalize.m, map_scale.m, map_erlang.m,
 * map_exponential.m, map_hyperexp.m, map_sum.m, map_sumind.m, map_mixture.m,
 * map_renewal.m, map_stochcomp.m, map2ph.m, map_skew.m, map_kurt.m,
 * map_joint.m, map_isfeasible.m).
 *
 * Everything here except map_hyperexp and map_skew is a finite sequence of
 * field operations on the entries of (D0, D1) -- block assembly, a Kronecker
 * product, one linear solve -- so the exact instantiation carries the MAP
 * identities (row sums of D0 + D1 vanish, the embedded chain is stochastic)
 * with no residual at all. map_hyperexp needs a square root of the moment
 * discriminant and map_skew a square root of the SCV, so both are gated on
 * num_traits<T>::has_transcendental.
 *
 * Naming note: map_exponential in map_moment.h is rate-parameterized,
 * map_exponential(lambda), whereas the MATLAB map_exponential(MEAN) is
 * mean-parameterized. map_exponential_mean below is the MATLAB spelling; the
 * two differ by the reciprocal and are otherwise identical.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Clamp negative off-diagonal entries of D0 and negative entries of D1 to
 * zero, then rebuild the diagonal of D0 so that every row of D0 + D1 sums to
 * zero (map_normalize.m).
 */
template <class T>
Map<T> map_normalize(const Map<T>& in) {
    const T zero = num_traits<T>::from_int(0);
    Map<T> m = in;
    const std::size_t n = m.order();
    if (m.D1.rows() != n || m.D0.cols() != n || m.D1.cols() != n)
        throw InputError("map_normalize: D0 and D1 must be square and of equal order");
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j && m.D0(i, j) < zero) m.D0(i, j) = zero;
            if (m.D1(i, j) < zero) m.D1(i, j) = zero;
        }
    for (std::size_t i = 0; i < n; ++i) {
        m.D0(i, i) = zero;
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += m.D0(i, j) + m.D1(i, j);
        m.D0(i, i) = -s;
    }
    return m;
}

/** Rescale time so that the mean inter-arrival time becomes new_mean. */
template <class T>
Map<T> map_scale(const Map<T>& in, const T& new_mean) {
    if (new_mean == num_traits<T>::from_int(0)) throw InputError("map_scale: zero target mean");
    const T ratio = map_mean(in) / new_mean;
    Map<T> m = in;
    for (std::size_t i = 0; i < m.order(); ++i)
        for (std::size_t j = 0; j < m.order(); ++j) {
            m.D0(i, j) *= ratio;
            m.D1(i, j) *= ratio;
        }
    return map_normalize(m);
}

/**
 * Rescale to a target mean WITHOUT the feasibility repair, for a matrix
 * exponential.
 *
 * `map_scale` finishes with `map_normalize`, which zeroes every negative entry;
 * that is a repair for a MAP whose blocks drifted, and it is DESTRUCTION for an
 * ME or a RAP, whose negative off-diagonals are the representation. Scaling both
 * blocks by one positive factor already preserves every normalized moment and
 * the zero row sums of D0 + D1, so nothing needs repairing.
 *
 * The distinction is the reference's: `solver_mam_basic.m:82-88` branches on
 * `procid == ME || procid == RAP` and rescales by the rate alone. Without the
 * branch the clamp turns a two-moment CME fit into a different process
 * altogether -- measured on M/Pareto/1 at rho 0.5, which reported the queue
 * length of an infinite server (0.5) against the exact 0.95.
 */
template <class T>
Map<T> map_scale_rate(const Map<T>& in, const T& new_mean) {
    if (new_mean == num_traits<T>::from_int(0)) throw InputError("map_scale_rate: zero target mean");
    const T ratio = map_mean(in) / new_mean;
    Map<T> m = in;
    for (std::size_t i = 0; i < m.order(); ++i)
        for (std::size_t j = 0; j < m.order(); ++j) {
            m.D0(i, j) *= ratio;
            m.D1(i, j) *= ratio;
        }
    return m;
}

/** Poisson process with the given mean inter-arrival time (map_exponential.m). */
template <class T>
Map<T> map_exponential_mean(const T& mean) {
    if (mean == num_traits<T>::from_int(0)) throw InputError("map_exponential_mean: zero mean");
    return map_exponential(T(num_traits<T>::from_int(1) / mean));
}

/** Erlang-k renewal MAP with the given mean (map_erlang.m). */
template <class T>
Map<T> map_erlang(const T& mean, unsigned k) {
    if (k == 0) throw InputError("map_erlang: k must be positive");
    if (mean == num_traits<T>::from_int(0)) throw InputError("map_erlang: zero mean");
    const T mu = num_traits<T>::from_int(static_cast<long>(k)) / mean;
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = Matrix<T>(k, k, zero);
    m.D1 = Matrix<T>(k, k, zero);
    for (unsigned i = 0; i + 1 < k; ++i) m.D0(i, i + 1) = mu;
    m.D1(k - 1, 0) = mu;
    return map_normalize(m);
}

/**
 * Two-phase hyperexponential renewal MAP matching a mean and an SCV >= 1,
 * with branching probability p (map_hyperexp.m, default p = 0.99).
 *
 * Gated on transcendental arithmetic: the phase rates come from the root of a
 * quadratic moment-matching condition and need a square root, which has no
 * exact rational counterpart. MATLAB falls back to the second root and then to
 * a smaller p when the first root is infeasible; the same two fallbacks are
 * reproduced here, and an infeasible result throws rather than returning the
 * empty MAP that MATLAB returns.
 */
template <class T>
Map<T> map_hyperexp(const T& mean, const T& scv, const T& p_in) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_hyperexp requires transcendental arithmetic");
    using std::sqrt;
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T four = num_traits<T>::from_int(4);
    T p = p_in;
    for (unsigned attempt = 0; attempt < 8; ++attempt) {
        const T E2 = T((one + scv) * mean * mean);
        const T Delta = T(-four * p * mean * mean + four * p * p * mean * mean + two * E2 * p -
                          two * E2 * p * p);
        if (Delta >= num_traits<T>::from_int(0)) {
            const T sD = T(sqrt(Delta));
            const T den = T(E2 * p - two * mean * mean);
            if (den != num_traits<T>::from_int(0)) {
                for (int root = 0; root < 2; ++root) {
                    const T mu2 = T((root == 0 ? T(-two * mean + two * p * mean + sD)
                                               : T(-two * mean + two * p * mean - sD)) /
                                    den);
                    const T dd = T(p - one + mean * mu2);
                    if (dd == num_traits<T>::from_int(0)) continue;
                    const T mu1 = T(mu2 * p / dd);
                    Map<T> m;
                    m.D0 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
                    m.D1 = Matrix<T>(2, 2, num_traits<T>::from_int(0));
                    m.D0(0, 0) = -mu1;
                    m.D0(1, 1) = -mu2;
                    m.D1(0, 0) = mu1 * p;
                    m.D1(0, 1) = mu1 * (one - p);
                    m.D1(1, 0) = mu2 * p;
                    m.D1(1, 1) = mu2 * (one - p);
                    if (mu1 > num_traits<T>::from_int(0) && mu2 > num_traits<T>::from_int(0))
                        return m;
                }
            }
        }
        p /= num_traits<T>::from_int(10);
    }
    throw NumericError("map_hyperexp: no feasible two-phase fit for this (mean, scv)");
}

/** map_hyperexp with the MATLAB default branching probability p = 0.99. */
template <class T>
Map<T> map_hyperexp(const T& mean, const T& scv) {
    return map_hyperexp(mean, scv, T(num_traits<T>::from_rational(99, 100)));
}

/**
 * n-fold convolution of a MAP with itself: the inter-arrival time of the
 * result is the sum of n consecutive inter-arrival times (map_sum.m).
 */
template <class T>
Map<T> map_sum(const Map<T>& in, unsigned n) {
    if (n == 0) throw InputError("map_sum: n must be positive");
    const std::size_t order = in.order();
    const std::size_t N = order * n;
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = Matrix<T>(N, N, zero);
    m.D1 = Matrix<T>(N, N, zero);
    std::size_t cur = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (std::size_t a = 0; a < order; ++a)
            for (std::size_t b = 0; b < order; ++b) m.D0(cur + a, cur + b) = in.D0(a, b);
        if (i + 1 < n) {
            for (std::size_t a = 0; a < order; ++a)
                for (std::size_t b = 0; b < order; ++b)
                    m.D0(cur + a, cur + order + b) = in.D1(a, b);
        } else {
            for (std::size_t a = 0; a < order; ++a)
                for (std::size_t b = 0; b < order; ++b) m.D1(cur + a, b) = in.D1(a, b);
        }
        cur += order;
    }
    return m;
}

/**
 * Sum of independent, not necessarily identical MAPs: after each component
 * completes, the next one restarts from its own stationary arrival phase
 * distribution pie (map_sumind.m).
 */
template <class T>
Map<T> map_sumind(const std::vector<Map<T>>& maps) {
    if (maps.empty()) throw InputError("map_sumind: empty list");
    const std::size_t n = maps.size();
    std::vector<std::size_t> order(n), off(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = maps[i].order();
        off[i + 1] = off[i] + order[i];
    }
    const std::size_t N = off[n];
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = Matrix<T>(N, N, zero);
    m.D1 = Matrix<T>(N, N, zero);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t a = 0; a < order[i]; ++a)
            for (std::size_t b = 0; b < order[i]; ++b)
                m.D0(off[i] + a, off[i] + b) = maps[i].D0(a, b);
        // D1 of component i, collapsed to a column, times the entry vector of
        // the next component (or of the first one, closing the cycle).
        const std::size_t nxt = (i + 1 < n) ? i + 1 : 0;
        const std::vector<T> pie = map_pie(maps[nxt]);
        std::vector<T> rowsum(order[i], zero);
        for (std::size_t a = 0; a < order[i]; ++a)
            for (std::size_t b = 0; b < order[i]; ++b) rowsum[a] += maps[i].D1(a, b);
        for (std::size_t a = 0; a < order[i]; ++a)
            for (std::size_t b = 0; b < order[nxt]; ++b) {
                const T val = rowsum[a] * pie[b];
                if (i + 1 < n)
                    m.D0(off[i] + a, off[nxt] + b) = val;
                else
                    m.D1(off[i] + a, off[nxt] + b) = val;
            }
    }
    return m;
}

/**
 * Probabilistic mixture of MAPs with weights alpha: after an arrival from
 * component i the process jumps to component j with probability alpha(j),
 * entering at its stationary arrival phase (map_mixture.m).
 */
template <class T>
Map<T> map_mixture(const std::vector<T>& alpha, const std::vector<Map<T>>& maps) {
    if (maps.empty()) throw InputError("map_mixture: empty list");
    if (alpha.size() != maps.size()) throw InputError("map_mixture: weight/list size mismatch");
    const std::size_t n = maps.size();
    std::vector<std::size_t> order(n), off(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = maps[i].order();
        off[i + 1] = off[i] + order[i];
    }
    const std::size_t N = off[n];
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = Matrix<T>(N, N, zero);
    m.D1 = Matrix<T>(N, N, zero);
    std::vector<std::vector<T>> pies(n);
    for (std::size_t j = 0; j < n; ++j) pies[j] = map_pie(maps[j]);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t a = 0; a < order[i]; ++a)
            for (std::size_t b = 0; b < order[i]; ++b)
                m.D0(off[i] + a, off[i] + b) = maps[i].D0(a, b);
        std::vector<T> rowsum(order[i], zero);
        for (std::size_t a = 0; a < order[i]; ++a)
            for (std::size_t b = 0; b < order[i]; ++b) rowsum[a] += maps[i].D1(a, b);
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t a = 0; a < order[i]; ++a)
                for (std::size_t b = 0; b < order[j]; ++b)
                    m.D1(off[i] + a, off[j] + b) = rowsum[a] * alpha[j] * pies[j][b];
    }
    return map_normalize(m);
}

/**
 * Renewal process with the same inter-arrival distribution: D1 is replaced by
 * (D1 e) pie, which destroys the correlation but preserves every marginal
 * moment (map_renewal.m).
 */
template <class T>
Map<T> map_renewal(const Map<T>& in) {
    const std::size_t n = in.order();
    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> pie = map_pie(in);
    Map<T> m;
    m.D0 = in.D0;
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t a = 0; a < n; ++a) {
        T rs = zero;
        for (std::size_t b = 0; b < n; ++b) rs += in.D1(a, b);
        for (std::size_t b = 0; b < n; ++b) m.D1(a, b) = rs * pie[b];
    }
    return m;
}

/** A phase-type representation (alpha, T) of a MAP's inter-arrival time. */
template <class T>
struct PhType {
    std::vector<T> alpha;  ///< entry distribution, = pie
    Matrix<T> subgen;      ///< sub-generator, = D0
};

/** (alpha, T) of the inter-arrival distribution: alpha = pie, T = D0 (map2ph.m). */
template <class T>
PhType<T> map2ph(const Map<T>& in) {
    PhType<T> ph;
    ph.alpha = map_pie(in);
    ph.subgen = in.D0;
    return ph;
}

/** MAP whose inter-arrival time is the PH (alpha, T): the renewal MAP with D1 = (-T e) alpha. */
template <class T>
Map<T> ph2map(const PhType<T>& ph) {
    const std::size_t n = ph.subgen.rows();
    if (ph.alpha.size() != n) throw InputError("ph2map: alpha and T have different orders");
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = ph.subgen;
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t a = 0; a < n; ++a) {
        T rs = zero;
        for (std::size_t b = 0; b < n; ++b) rs += ph.subgen(a, b);
        for (std::size_t b = 0; b < n; ++b) m.D1(a, b) = -rs * ph.alpha[b];
    }
    return m;
}

/**
 * Stochastic complement of a MAP on the retained phases (map_stochcomp.m):
 * the eliminated phases are censored out of the generator and of D1, giving a
 * smaller MAP with the same behaviour observed on the retained phases.
 */
template <class T>
Map<T> map_stochcomp(const Map<T>& in, const std::vector<std::size_t>& retain) {
    const std::size_t n = in.order();
    std::vector<bool> keep(n, false);
    for (std::size_t r : retain) {
        if (r >= n) throw InputError("map_stochcomp: retained index out of range");
        keep[r] = true;
    }
    std::vector<std::size_t> elim;
    for (std::size_t i = 0; i < n; ++i)
        if (!keep[i]) elim.push_back(i);
    const std::size_t nr = retain.size(), ne = elim.size();
    if (nr == 0) throw InputError("map_stochcomp: no phase retained");
    const Matrix<T> Q = map_infgen(in);
    Map<T> out;
    out.D0 = Matrix<T>(nr, nr);
    out.D1 = Matrix<T>(nr, nr);
    if (ne == 0) {
        for (std::size_t a = 0; a < nr; ++a)
            for (std::size_t b = 0; b < nr; ++b) {
                out.D0(a, b) = in.D0(retain[a], retain[b]);
                out.D1(a, b) = in.D1(retain[a], retain[b]);
            }
        return map_normalize(out);
    }
    Matrix<T> Q_EE(ne, ne), Q_RE(nr, ne), Q_ER(ne, nr), D1_ER(ne, nr);
    for (std::size_t a = 0; a < ne; ++a)
        for (std::size_t b = 0; b < ne; ++b) Q_EE(a, b) = Q(elim[a], elim[b]);
    for (std::size_t a = 0; a < nr; ++a)
        for (std::size_t b = 0; b < ne; ++b) Q_RE(a, b) = Q(retain[a], elim[b]);
    for (std::size_t a = 0; a < ne; ++a)
        for (std::size_t b = 0; b < nr; ++b) {
            Q_ER(a, b) = Q(elim[a], retain[b]);
            D1_ER(a, b) = in.D1(elim[a], retain[b]);
        }
    Matrix<T> negQEE(ne, ne);
    for (std::size_t a = 0; a < ne; ++a)
        for (std::size_t b = 0; b < ne; ++b) negQEE(a, b) = -Q_EE(a, b);
    const Matrix<T> inv = inverse(negQEE);
    const Matrix<T> corrQ = matmul(Q_RE, matmul(inv, Q_ER));
    const Matrix<T> corrD1 = matmul(Q_RE, matmul(inv, D1_ER));
    for (std::size_t a = 0; a < nr; ++a)
        for (std::size_t b = 0; b < nr; ++b) {
            const T Qnew = Q(retain[a], retain[b]) + corrQ(a, b);
            out.D1(a, b) = in.D1(retain[a], retain[b]) + corrD1(a, b);
            out.D0(a, b) = Qnew - out.D1(a, b);
        }
    return map_normalize(out);
}

/** Kurtosis of the inter-arrival time (map_kurt.m); rational in the entries. */
template <class T>
T map_kurt(const Map<T>& m) {
    const T m1 = map_moment(m, 1);
    const T m2 = map_moment(m, 2);
    const T m3 = map_moment(m, 3);
    const T m4 = map_moment(m, 4);
    const T v = map_var(m);
    const T num = T(m4 - num_traits<T>::from_int(4) * m3 * m1 +
                    num_traits<T>::from_int(6) * m2 * m1 * m1 -
                    num_traits<T>::from_int(3) * m1 * m1 * m1 * m1);
    return num / (v * v);
}

/**
 * Skewness of the inter-arrival time (map_skew.m). Gated on transcendental
 * arithmetic: the denominator is (sqrt(SCV) * mean)^3, and the square root of
 * a rational SCV is irrational in general.
 */
template <class T>
T map_skew(const Map<T>& m) {
    static_assert(num_traits<T>::has_transcendental, "map_skew requires transcendental arithmetic");
    using std::sqrt;
    const T m1 = map_moment(m, 1);
    const T m2 = map_moment(m, 2);
    const T m3 = map_moment(m, 3);
    const T M3 = T(m3 - num_traits<T>::from_int(3) * m2 * m1 +
                   num_traits<T>::from_int(2) * m1 * m1 * m1);
    const T s = T(sqrt(map_scv(m)));
    const T den = T(s * m1);
    return M3 / (den * den * den);
}

/**
 * Joint moment of K consecutive inter-arrival times observed at the cumulative
 * lags a, with orders i (map_joint.m). a and i have the same length K; a[0] is
 * ignored as a base point exactly as MATLAB's cumsum makes it.
 */
template <class T>
T map_joint(const Map<T>& m, const std::vector<unsigned>& a, const std::vector<unsigned>& i) {
    if (a.size() != i.size() || a.empty()) throw InputError("map_joint: a and i size mismatch");
    const std::size_t n = m.order();
    const std::size_t K = a.size();
    std::vector<unsigned> ca(K);
    unsigned acc = 0;
    for (std::size_t k = 0; k < K; ++k) {
        acc += a[k];
        ca[k] = acc;
    }
    Matrix<T> negD0 = m.D0;
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c) negD0(r, c) = -negD0(r, c);
    const Matrix<T> invD0 = inverse(negD0);
    const Matrix<T> P = matmul(invD0, m.D1);
    Matrix<T> JM = eye<T>(n);
    for (std::size_t k = 0; k + 1 < K; ++k) {
        Matrix<T> blk = matpow(invD0, i[k]);
        const T f = num_factorial<T>(i[k]);
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = 0; c < n; ++c) blk(r, c) *= f;
        JM = matmul(JM, matmul(blk, matpow(P, ca[k + 1] - ca[k])));
    }
    Matrix<T> last = matpow(invD0, i[K - 1]);
    const T fl = num_factorial<T>(i[K - 1]);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c) last(r, c) *= fl;
    JM = matmul(JM, last);
    const std::vector<T> v = vecmul(map_pie(m), JM);
    T s = num_traits<T>::from_int(0);
    for (const T& x : v) s += x;
    return s;
}

/**
 * Structural feasibility of a MAP within a tolerance (map_isfeasible.m):
 * off-diagonal D0 and all of D1 non-negative, diagonal of D0 non-positive,
 * D0 + D1 a generator, and the embedded chain P = (-D0)^-1 D1 non-negative and
 * stochastic. The eigenvalue-multiplicity screen of the MATLAB routine is not
 * reproduced -- it needs a full eigendecomposition, which is not in this port;
 * a MAP that passes here can in principle still have a defective generator.
 */
template <class T>
bool map_isfeasible(const Map<T>& m, const T& tol) {
    const std::size_t n = m.order();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T ntol = T(num_traits<T>::from_int(static_cast<long>(n)) * tol);
    // The reference's first guard: a NaN or an infinity anywhere makes every
    // test below meaningless, and a NaN comparison is false in BOTH directions,
    // so without this an unusable pair reads as feasible.
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double a = num_traits<T>::to_double(m.D0(i, j));
            const double b = num_traits<T>::to_double(m.D1(i, j));
            if (std::isnan(a) || std::isnan(b) || std::isinf(a) || std::isinf(b)) return false;
        }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j && m.D0(i, j) < -tol) return false;
            if (i == j && m.D0(i, j) > tol) return false;
            if (m.D1(i, j) < -tol) return false;
        }
    const Matrix<T> Q = map_infgen(m);
    for (std::size_t i = 0; i < n; ++i) {
        T rs = zero;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j && Q(i, j) < -tol) return false;
            rs += Q(i, j);
        }
        if (num_abs(T(rs)) > ntol) return false;
    }
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negD0(i, j) = -negD0(i, j);
    Matrix<T> P;
    try {
        P = matmul(inverse(negD0), m.D1);
    } catch (const NumericError&) {
        return false;  // -D0 singular: no embedded chain, hence infeasible
    }
    for (std::size_t i = 0; i < n; ++i) {
        T rs = zero;
        for (std::size_t j = 0; j < n; ++j) {
            if (P(i, j) < -tol) return false;
            rs += P(i, j);
        }
        if (num_abs(T(rs - one)) > ntol) return false;
    }
    return true;
}

/**
 * The reference's `map_checkfeasible`, i.e. `map_isfeasible` at a given
 * tolerance. It is a nested function of map_isfeasible.m and a top-level class
 * in the JAR; both spellings name the same test.
 */
template <class T>
bool map_checkfeasible(const Map<T>& m, const T& tol) {
    return map_isfeasible(m, tol);
}

/**
 * `map_isfeasible(MAP)` with no tolerance, which is NOT the zero-tolerance test.
 *
 * The reference scans k from 15 downwards, i.e. from the tightest tolerance
 * 1e-15 to the loosest 1e-1, stops at the first k whose `map_checkfeasible`
 * passes, and reports feasible iff that tightest passing tolerance is tighter
 * than `map_feastol` (1e-8). A MAP whose row sums are exact only to rounding is
 * therefore feasible to the reference and INFEASIBLE to a zero-tolerance test.
 *
 * This overload used to call the zero-tolerance form, which no assembled MAP
 * can pass in floating point: `map_block` was sent to its fallback on moment
 * sets MATLAB fits exactly, purely on the row sums' last bits.
 */
template <class T>
bool map_isfeasible(const Map<T>& m) {
    for (int k = 15; k >= 1; --k) {
        const T tol = num_traits<T>::from_double(std::pow(10.0, -static_cast<double>(k)));
        if (map_isfeasible(m, tol)) return k > map_feastol();
    }
    return false;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_TRANSFORM_H
