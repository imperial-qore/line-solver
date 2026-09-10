/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CLW_H
#define LINE_API_PFQN_CLW_H

/**
 * Choudhury-Leung-Whitt normalization constant by numerical inversion of the
 * generating function (JACM 42(5):935-970, 1995), and its limited
 * load-dependent extension through the per-center transforms of Bertozzi and
 * McKenna (SIAM Review 35(2):239-268, 1993).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_clw.m and pfqn_clw_lld.m.
 *
 * The generating function of g(K) is (CLW eq. 4.5)
 *
 *   G(z) = exp( sum_j rho_{j0} z_j ) / prod_i ( 1 - sum_j rho_{ji} z_j )^{m_i}
 *
 * and, with limited load-dependent stations (Bertozzi-McKenna 2.17/2.23),
 *
 *   G(z) = exp( sum_j rho_{j0} z_j ) prod_i F_i( sum_j rho_{ji} z_j ),
 *   F_i(x) = [ c_i + sum_{n=1}^{l_i-1} (c_i - S_i(n)) / prod_{k<=n} S_i(k) x^n ]
 *            / (c_i - x),
 *
 * analytic except for a simple pole at x = c_i. g(K) is the coefficient of
 * prod_j z_j^{K_j}, recovered by p NESTED one-dimensional lattice-Poisson
 * inversions (eq. 2.3) on contours of radius r_j = 10^{-gamma_j/(2 l_j K_j)},
 * with the restrictive static scaling of eqs. 5.41-5.46 and log-domain
 * recovery (eq. 7.1). pfqn_clw applies both of the paper's speed-ups: dimension
 * reduction by decomposition (Sec. 3, Sec. 5.4), which inverts the subset D
 * minimizing |D| + max_i |S_i(D)| (eq. 3.3) and then each connected component of
 * the remainder separately, and Euler summation of the inner sums (Sec. 2.4,
 * eq. 2.22), which replaces 2 l_j K_j contour points by 2 l_j (n+m+1) wherever
 * K_j > n+m, refining m until |E(m,n) - E(m,n+1)| settles. pfqn_clw_lld keeps
 * the plain nested inversion of cost prod_j 2 l_j K_j.
 *
 * COMPLEX ARITHMETIC WITHOUT std::complex. The contour integrand is genuinely
 * complex, and std::complex is specified only for float, double and long
 * double; instantiating it on a Boost.Multiprecision number is unspecified
 * behavior. detail::Cx<T> below is a two-field complex with the six operations
 * this routine needs, so the Real backends get real high precision on the
 * inversion rather than silently falling back to double.
 *
 * Arithmetic: TRANSCENDENTAL, double and Real only. exp, log, atan2, sqrt and
 * a fractional power all appear; the contour radius alone is 10^{-gamma/(2 l
 * K)}, which is not in the field of the inputs. This is the one member of the
 * normalizing-constant family that CANNOT be instantiated exactly, and that is
 * intrinsic to inverting a generating function numerically, not an artifact of
 * the port. Accuracy against the exact pfqn_ca on the models tested is ~4e-10
 * for p = 2, ~7e-7 for p = 3, matching the reference's own claim of about 1e-9
 * and confirming that the port reproduces the reference's aliasing rather than
 * adding error of its own.
 *
 * REFERENCE DEFECTS
 *
 *  1. pfqn_clw returns NaN when SOME chain has zero population. The contour
 *     count is 2 l_j K_j, so a chain with K_j = 0 makes the final division
 *     acc / (2 l_j K_j r_j^{K_j}) a 0/0. pfqn_clw_lld guards against exactly
 *     this by dropping the zero-population chains up front ("the coefficient of
 *     z_j^0 equals the pgf restricted to z_j = 0, so chain j is removed
 *     exactly"); pfqn_clw never received that guard. Reproduce with
 *     pfqn_clw([0.1 0.2; 0.3 0.05], [2 0], [1.0 0.5]), which returns NaN where
 *     the answer is the p = 1 constant. THIS PORT APPLIES THE GUARD to both
 *     routines, so pfqn_clw here returns the finite value; that is the only
 *     input on which the two disagree.
 *  2. Dead code in both routines: alpha0_j = exp(-alpha_j rho_{j0}) is computed
 *     in the scaling loop and never read. The recovery (eq. 7.1) is written in
 *     terms of sum_j alpha_j rho_{j0}, i.e. -sum_j log alpha0_j, so the array
 *     is redundant rather than wrong. It is not carried here.
 *  3. `denom(denom <= 0) = eps` in the scaling loop silently substitutes
 *     2.2e-16 for a nonpositive deflated denominator, which turns a chain whose
 *     predecessors have already saturated a queue into an enormous effective
 *     intensity rather than reporting the saturation. Reproduced, because the
 *     scaling only has to keep the contour inside the disc of analyticity and
 *     the recovery divides the choice back out, so the result is unaffected;
 *     but it is a silent branch, not a designed one.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_clw and pfqn_clw_lld, mirroring [G, lG]. */
template <class T>
struct ClwResult {
    T G;    ///< normalization constant, +infinity when it overflows the range of T
    T lG;   ///< its natural logarithm, always finite
};

/** Optional lattice and aliasing parameters; empty means "use the CLW defaults". */
struct ClwOptions {
    std::vector<int> l;       ///< inner lattice parameters l_j (roundoff control)
    std::vector<double> gamma;  ///< aliasing parameters, aliasing ~ 10^-gamma_j
    bool euler = true;        ///< Euler-sum the inner sums where K_j > eulerN + eulerM
    int eulerN = 11;          ///< terms summed exactly before averaging (n in eq. 2.22)
    int eulerM = 20;          ///< starting order of the averaging (m in eq. 2.22)
    double eulerTol = 1e-10;  ///< relative tolerance on |E(m,n) - E(m,n+1)|
    int eulerMaxM = 160;      ///< largest Euler order reached by doubling
    bool dimred = true;       ///< dimension reduction by decomposition (Section 3)
    int dimredMaxD = 4;       ///< largest |D| examined when minimizing (3.3)
    std::vector<double> beta;  ///< multipliers on alpha_j, the manual tuning of page 956
};

namespace detail {

/**
 * Cx<T>, cx_add, cx_mul, cx_scale, cx_div and cx_expi come from
 * pfqn_asympt_common.h, which introduced the same two-field complex for the
 * Norlund-Rice integrands. Only the three operations that inversion needs and
 * that file does not have are added here.
 */
template <class T>
Cx<T> cx_sub(const Cx<T>& a, const Cx<T>& b) {
    return Cx<T>(T(a.re - b.re), T(a.im - b.im));
}

/** Principal complex logarithm, log|z| + i arg z. */
template <class T>
Cx<T> cx_log(const Cx<T>& a) {
    using std::atan2;
    using std::log;
    using std::sqrt;
    const T mod = sqrt(T(a.re * a.re + a.im * a.im));
    return Cx<T>(T(log(mod)), T(atan2(a.im, a.re)));
}

template <class T>
Cx<T> cx_exp(const Cx<T>& a) {
    using std::cos;
    using std::exp;
    using std::sin;
    const T e = exp(a.re);
    return Cx<T>(T(e * cos(a.im)), T(e * sin(a.im)));
}

/**
 * x^y for positive x. std::pow, not exp(y log x): the latter costs a couple of
 * ulps on the contour radius r_j = 10^{-gamma_j/(2 l_j K_j)}, and the inversion
 * is an alternating sum that amplifies exactly that kind of error.
 */
template <class T>
T clw_pow_real(const T& x, const T& y) {
    using std::pow;
    return T(pow(x, y));
}

/** The page-956 scale multipliers beta_j, in the retained chain index space. */
template <class T>
std::vector<T> clw_beta(const std::vector<std::size_t>& keep, const ClwOptions& opt) {
    std::vector<T> beta(keep.size(), num_traits<T>::from_int(1));
    if (opt.beta.empty()) return beta;
    for (std::size_t j = 0; j < keep.size(); ++j) {
        if (keep[j] >= opt.beta.size()) throw InputError("pfqn_clw: options.beta has the wrong length");
        beta[j] = num_traits<T>::from_double(opt.beta[keep[j]]);
    }
    return beta;
}

/**
 * Inversion order of the variables: the subset D of the interdependence graph
 * that is inverted first, then the connected components of what is left
 * (CLW Section 3). Indices are into the retained chains.
 */
struct ClwPlan {
    std::vector<std::size_t> D;
    std::vector<std::vector<std::size_t>> comps;
};

/** Connected components of adj with the nodes in mask removed. */
inline std::vector<std::vector<std::size_t>> clw_components(
    const std::vector<std::vector<char>>& adj, const std::vector<char>& mask) {
    const std::size_t p = adj.size();
    std::vector<long> lab(p, -1);
    long nc = 0;
    for (std::size_t s = 0; s < p; ++s) {
        if (mask[s] || lab[s] >= 0) continue;
        lab[s] = nc;
        std::vector<std::size_t> stack(1, s);
        while (!stack.empty()) {
            const std::size_t v = stack.back();
            stack.pop_back();
            for (std::size_t u = 0; u < p; ++u)
                if (adj[v][u] && !mask[u] && lab[u] < 0) {
                    lab[u] = nc;
                    stack.push_back(u);
                }
        }
        ++nc;
    }
    std::vector<std::vector<std::size_t>> out(static_cast<std::size_t>(nc));
    for (std::size_t j = 0; j < p; ++j)
        if (lab[j] >= 0) out[static_cast<std::size_t>(lab[j])].push_back(j);
    return out;
}

inline double clw_binom(std::size_t n, std::size_t k) {
    double v = 1.0;
    for (std::size_t i = 1; i <= k; ++i) v = v * static_cast<double>(n - k + i) / static_cast<double>(i);
    return v;
}

/**
 * The interdependence graph of the factors of (4.5) and the subset D minimizing
 * the inversion dimension |D| + max_i |S_i(D)| (eqs. 3.1-3.3), by enumeration in
 * increasing cardinality. A plan of dimension p is no reduction at all, and is
 * returned as the single all-chain component so the recursion is unchanged.
 */
template <class T>
ClwPlan clw_plan(const Matrix<T>& L, const ClwOptions& opt) {
    const std::size_t qd = L.rows(), p = L.cols();
    ClwPlan trivial;
    trivial.comps.resize(1);
    for (std::size_t j = 0; j < p; ++j) trivial.comps[0].push_back(j);
    if (!opt.dimred || p <= 2) return trivial;
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::vector<char>> adj(p, std::vector<char>(p, 0));
    for (std::size_t i = 0; i < qd; ++i)
        for (std::size_t a = 0; a < p; ++a) {
            if (L(i, a) == zero) continue;
            for (std::size_t b = 0; b < p; ++b)
                if (b != a && L(i, b) != zero) adj[a][b] = 1;   // each factor is a clique
        }
    const std::vector<char> none(p, 0);
    std::vector<std::vector<std::size_t>> bestC = clw_components(adj, none);
    std::size_t best = 0;
    for (std::size_t c = 0; c < bestC.size(); ++c) best = std::max(best, bestC[c].size());
    std::vector<std::size_t> bestD;
    const std::size_t maxd = std::min<std::size_t>(
        (opt.dimredMaxD > 0) ? static_cast<std::size_t>(opt.dimredMaxD) : 0, p - 1);
    for (std::size_t dd = 1; dd <= maxd; ++dd) {
        if (dd >= best) break;                       // dimension is at least |D|
        if (clw_binom(p, dd) > 2e5) break;           // (3.3) is solved by enumeration only
        std::vector<std::size_t> sub(dd);
        for (std::size_t t = 0; t < dd; ++t) sub[t] = t;
        while (true) {
            std::vector<char> mask(p, 0);
            for (std::size_t t = 0; t < dd; ++t) mask[sub[t]] = 1;
            const std::vector<std::vector<std::size_t>> cc = clw_components(adj, mask);
            std::size_t mx = 0;
            for (std::size_t c = 0; c < cc.size(); ++c) mx = std::max(mx, cc[c].size());
            if (dd + mx < best) {
                best = dd + mx;
                bestD = sub;
                bestC = cc;
            }
            std::size_t t = dd;
            while (t-- > 0 && sub[t] == p - dd + t) {
            }
            if (t >= dd) break;                      // wrapped: enumeration exhausted
            ++sub[t];
            for (std::size_t u = t + 1; u < dd; ++u) sub[u] = sub[u - 1] + 1;
        }
    }
    if (best >= p) return trivial;
    ClwPlan plan;
    plan.D = bestD;
    plan.comps = bestC;
    return plan;
}

/**
 * The CLW lattice and aliasing parameters, defaulted by INVERSION DEPTH: the
 * dimension reduction sets the order of the variables (Section 5.4), and every
 * component restarts at depth |D|+1 because the components are inverted in
 * parallel. depth is 1-based and indexed by retained chain; keep maps retained
 * chains back to the caller's index space, which is where opt.l lives.
 */
inline void clw_defaults(std::size_t pfull, const std::vector<std::size_t>& keep,
                         const std::vector<std::size_t>& depth, const ClwOptions& opt,
                         std::vector<int>& l, std::vector<double>& gam) {
    const std::size_t p = keep.size();
    l.assign(p, 3);
    gam.assign(p, 15.0);
    for (std::size_t j = 0; j < p; ++j) {
        if (depth[j] == 1) {
            l[j] = 1;
            gam[j] = 11.0;
        } else if (depth[j] <= 3) {
            l[j] = 2;
            gam[j] = 13.0;
        }
    }
    if (!opt.l.empty()) {
        if (opt.l.size() != pfull) throw InputError("pfqn_clw: options.l has the wrong length");
        for (std::size_t j = 0; j < p; ++j) l[j] = opt.l[keep[j]];
    }
    if (!opt.gamma.empty()) {
        if (opt.gamma.size() != pfull)
            throw InputError("pfqn_clw: options.gamma has the wrong length");
        for (std::size_t j = 0; j < p; ++j) gam[j] = opt.gamma[keep[j]];
    }
    for (std::size_t j = 0; j < p; ++j)
        if (l[j] < 1) throw InputError("pfqn_clw: the lattice parameters must be positive");
}

/**
 * Positional form of the defaults, for the members of the family that invert in
 * chain order and take no plan: pfqn_clwoi, pfqn_clwjd and the pfqn_ncld cost
 * estimate.
 */
inline void clw_defaults(std::size_t p, const ClwOptions& opt, std::vector<int>& l,
                         std::vector<double>& gam) {
    std::vector<std::size_t> keep(p, 0), depth(p, 0);
    for (std::size_t j = 0; j < p; ++j) {
        keep[j] = j;
        depth[j] = j + 1;
    }
    clw_defaults(p, keep, depth, opt, l, gam);
}

/**
 * The restrictive static scaling of CLW eqs. 5.41-5.46.
 *
 * @param Lsc  intensities the constraint is phrased on (rho for pfqn_clw, the
 *             unit-pole rho/c_i for pfqn_clw_lld)
 * @param mult per-queue multiplicity entering N_{ij}; all ones in the LLD form,
 *             where every F_i contributes a single pole
 * @param Lraw the unscaled demands, before the chain-level scaling
 * @param N (R) population per chain
 * @param Z (R) think times
 * @param l per-chain station index list
 * @param r per-chain scaling factors
 */
template <class T>
std::vector<T> clw_scaling(const Matrix<T>& Lsc, const Matrix<T>& Lraw, const std::vector<int>& N,
                           const std::vector<T>& Z, const std::vector<int>& l,
                           const std::vector<T>& r, const std::vector<long>& mult,
                           const std::vector<std::size_t>& order, const std::vector<T>& beta) {
    const std::size_t qd = Lsc.rows(), p = Lsc.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T eps = num_traits<T>::from_double(std::numeric_limits<double>::epsilon());
    std::vector<T> alpha(p, one), used(qd, zero);
    // in inversion order, which is what the dimension reduction changes (Sec. 5.4)
    for (std::size_t t = 0; t < p; ++t) {
        const std::size_t j = order[t];
        const long Kj = N[j], lj = l[j];
        if (Kj == 0) continue;  // empty chain: no lattice, and 2*lj*Kj = 0 below
        std::vector<std::size_t> posq;
        std::vector<T> e(qd, zero);
        for (std::size_t i = 0; i < qd; ++i) {
            T den = one - used[i];
            if (!(den > zero)) den = eps;
            e[i] = Lsc(i, j) / den;
            if (Lsc(i, j) > zero) posq.push_back(i);
        }
        bool have = false;
        T aj = zero;
        if (!posq.empty()) {
            std::stable_sort(posq.begin(), posq.end(),
                             [&](std::size_t a, std::size_t b) { return e[b] < e[a]; });
            T cum = zero;
            long cummb = 0;
            for (std::size_t n = 0; n < posq.size(); ++n) {
                const std::size_t qi = posq[n];
                cum += e[qi];
                cummb += mult[qi];
                const T rhobar = cum / num_traits<T>::from_int(static_cast<long>(n) + 1);
                long Nn = cummb - 1;
                for (std::size_t u = t + 1; u < p; ++u)
                    if (Lraw(qi, order[u]) != zero) Nn += N[order[u]];
                T an = one;
                if (Nn > 0) {
                    // in the log domain: the product runs over N_{ij} factors below
                    // one and underflows to zero at a few hundred of them, which
                    // would silently set alpha_j = 0 and lG = NaN
                    using std::exp;
                    using std::log;
                    T lp = zero;
                    for (long ll = 1; ll <= Nn; ++ll)
                        lp += log(T(num_traits<T>::from_int(Kj + ll) /
                                    num_traits<T>::from_int(Kj + 2 * lj * Kj + ll)));
                    an = exp(T(lp / num_traits<T>::from_int(2 * lj * Kj)));
                }
                const T cand = an / rhobar;
                if (!have || cand < aj) {
                    aj = cand;
                    have = true;
                }
            }
        }
        if (Z[j] > zero) {
            const T cand = num_traits<T>::from_int(Kj) / Z[j];
            if (!have || cand < aj) {
                aj = cand;
                have = true;
            }
        }
        if (!have) aj = one;  // chain with no demand anywhere
        alpha[j] = T(beta[j] * aj);
        for (std::size_t i = 0; i < qd; ++i) used[i] += alpha[j] * Lsc(i, j) * r[j];
    }
    return alpha;
}

/** Scaling in chain order with no page-956 tuning, for the same three callers. */
template <class T>
std::vector<T> clw_scaling(const Matrix<T>& Lsc, const Matrix<T>& Lraw, const std::vector<int>& N,
                           const std::vector<T>& Z, const std::vector<int>& l,
                           const std::vector<T>& r, const std::vector<long>& mult) {
    std::vector<std::size_t> order(Lsc.cols(), 0);
    for (std::size_t j = 0; j < order.size(); ++j) order[j] = j;
    return clw_scaling(Lsc, Lraw, N, Z, l, r, mult, order,
                       std::vector<T>(order.size(), num_traits<T>::from_int(1)));
}

/**
 * Euler weights of eq. (2.22): E(m,n) = sum_i w_i (-1)^i a_i.
 *
 * E(m,n) = 2^-m sum_{k=0}^{m} C(m,k) S_{n+k} with S_t = sum_{i<=t} (-1)^i a_i,
 * so a_i carries the mass of every partial sum that contains it.
 */
template <class T>
std::vector<T> clw_euler_weights(int n, int mm) {
    const T one = num_traits<T>::from_int(1);
    std::vector<T> b(static_cast<std::size_t>(mm) + 1, one);
    for (int k = 1; k <= mm; ++k)
        b[static_cast<std::size_t>(k)] =
            T(b[static_cast<std::size_t>(k) - 1] * num_traits<T>::from_int(mm - k + 1) /
              num_traits<T>::from_int(k));                       // C(mm,k)
    const T scale = T(one / num_pow_int(num_traits<T>::from_int(2), static_cast<unsigned>(mm)));
    for (int k = 0; k <= mm; ++k) b[static_cast<std::size_t>(k)] *= scale;
    std::vector<T> tail(static_cast<std::size_t>(mm) + 1, num_traits<T>::from_int(0));
    T acc = num_traits<T>::from_int(0);
    for (int k = mm; k >= 0; --k) {
        acc += b[static_cast<std::size_t>(k)];
        tail[static_cast<std::size_t>(k)] = acc;                 // 2^-mm sum_{j>=k} C(mm,j)
    }
    std::vector<T> w(static_cast<std::size_t>(n + mm) + 1, one);
    for (int i = n + 1; i <= n + mm; ++i)
        w[static_cast<std::size_t>(i)] = tail[static_cast<std::size_t>(i - n)];
    return w;
}

/**
 * The inner sum of (2.3) over the lattice index k, with Euler summation.
 *
 * The sum splits at k = 0 into two nearly alternating series (Section 2.4);
 * each is replaced by its Euler sum (eq. 2.22). The order m is doubled until the
 * paper's own estimate |E(m,n) - E(m,n+1)| falls under the tolerance, and the
 * exact sum is taken once n+m reaches K_j, so accuracy is not traded away.
 * `ev(k)` returns the inverted function at the lattice point of index k.
 */
template <class T, class EV>
Cx<T> clw_inner(long Kj, const ClwOptions& opt, const EV& ev) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    using std::sqrt;
    int mCur = opt.eulerM;
    while (true) {
        const long TT = static_cast<long>(opt.eulerN) + mCur;
        if (!opt.euler || Kj <= TT + 1) {
            Cx<T> s(zero, zero);
            for (long k = -Kj; k <= Kj - 1; ++k)
                s = cx_add(s, cx_scale(ev(k), (k % 2 == 0) ? one : T(-one)));
            return s;
        }
        const std::vector<T> w1 = clw_euler_weights<T>(opt.eulerN, mCur);
        const std::vector<T> w2 = clw_euler_weights<T>(opt.eulerN + 1, mCur);
        Cx<T> e1(zero, zero), e2(zero, zero);
        for (long s = 0; s <= TT + 1; ++s) {
            const Cx<T> dv = cx_sub(ev(s), ev(-(s + 1)));
            const T sgn = (s % 2 == 0) ? one : T(-one);
            if (s <= TT) e1 = cx_add(e1, cx_scale(dv, T(sgn * w1[static_cast<std::size_t>(s)])));
            e2 = cx_add(e2, cx_scale(dv, T(sgn * w2[static_cast<std::size_t>(s)])));
        }
        const Cx<T> df = cx_sub(e1, e2);
        const T dn = sqrt(T(df.re * df.re + df.im * df.im));
        const T en = sqrt(T(e2.re * e2.re + e2.im * e2.im));
        if (dn <= num_traits<T>::from_double(opt.eulerTol) * en || mCur >= opt.eulerMaxM) return e2;
        mCur *= 2;
    }
}

/** Everything the nested inversion needs besides the contour point itself. */
template <class T>
struct ClwCtx {
    const std::vector<int>* N;
    const std::vector<int>* l;
    const std::vector<T>* r;
    const ClwPlan* plan;
    const ClwOptions* opt;
    T pi;
};

template <class T, class F>
Cx<T> clw_invert_d(std::size_t t, std::vector<Cx<T>>& w, const ClwCtx<T>& cx, const F& gbar);

/**
 * One lattice-Poisson inversion (CLW eq. 2.3), scaled: extracts the coefficient
 * of w_j^{K_j} from whatever `next` evaluates at the contour points.
 */
template <class T, class NEXT>
Cx<T> clw_lattice(std::size_t j, std::vector<Cx<T>>& w, const ClwCtx<T>& cx, const NEXT& next) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const long Kj = (*cx.N)[j], lj = (*cx.l)[j];
    const T rj = (*cx.r)[j];
    if (Kj == 0) {
        // [w_j^0] Gbar = Gbar(w_j=0): the K=0 lattice is the single point 0, and
        // exp(-arho0_j) there cancels the +arho0_j added back into lG
        w[j] = Cx<T>(zero, zero);
        return next();
    }
    Cx<T> acc(zero, zero);
    for (long k1 = 0; k1 < lj; ++k1) {
        const Cx<T> ph =
            cx_expi(T(-cx.pi * num_traits<T>::from_int(k1) / num_traits<T>::from_int(lj)));
        const Cx<T> inner = clw_inner<T>(Kj, *cx.opt, [&](long k) {
            const T theta = T(cx.pi * num_traits<T>::from_int(k1 + lj * k) /
                              num_traits<T>::from_int(lj * Kj));
            w[j] = cx_scale(cx_expi(theta), rj);
            return next();
        });
        acc = cx_add(acc, cx_mul(ph, inner));
    }
    const T den =
        T(num_traits<T>::from_int(2 * lj * Kj) * num_pow_int(rj, static_cast<unsigned>(Kj)));
    return cx_scale(acc, T(one / den));
}

/** Inversion of one component of the interdependence graph minus D. */
template <class T, class F>
Cx<T> clw_invert_c(std::size_t c, std::size_t s, std::vector<Cx<T>>& w, const ClwCtx<T>& cx,
                   const F& gbar) {
    const std::vector<std::size_t>& vars = cx.plan->comps[c];
    if (s >= vars.size()) return gbar(w, static_cast<long>(c));
    Cx<T> val = clw_lattice(vars[s], w, cx, [&]() { return clw_invert_c(c, s + 1, w, cx, gbar); });
    if (cx.plan->D.empty() && s == 0) {
        // with D empty every component is an independent subnetwork, so its
        // coefficient is real
        val.im = num_traits<T>::from_int(0);
    }
    return val;
}

/** Outer inversion over the committed variables D (Section 3). */
template <class T, class F>
Cx<T> clw_invert_d(std::size_t t, std::vector<Cx<T>>& w, const ClwCtx<T>& cx, const F& gbar) {
    if (t >= cx.plan->D.size()) {
        // D fixed: the remaining factors have no variable in common, so the
        // coefficient of the inner monomial is the product of the components'
        Cx<T> val = gbar(w, -1);
        for (std::size_t c = 0; c < cx.plan->comps.size(); ++c)
            val = cx_mul(val, clw_invert_c(c, 0, w, cx, gbar));
        return val;
    }
    Cx<T> val =
        clw_lattice(cx.plan->D[t], w, cx, [&]() { return clw_invert_d(t + 1, w, cx, gbar); });
    if (t == 0) val.im = num_traits<T>::from_int(0);
    return val;
}

/**
 * The plain nested inversion of eq. (2.3) over all p variables in order, with
 * neither acceleration: the entry point pfqn_clwoi and pfqn_clwjd use, whose
 * `gbar` takes the contour point alone.
 */
template <class T, class F>
Cx<T> clw_invert(std::size_t j, std::vector<Cx<T>>& w, const std::vector<int>& N,
                 const std::vector<int>& l, const std::vector<T>& r, std::size_t p,
                 const F& gbar) {
    using std::acos;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    ClwPlan plan;
    plan.comps.resize(1);
    for (std::size_t u = j; u < p; ++u) plan.comps[0].push_back(u);
    ClwOptions opt;
    opt.euler = false;
    opt.dimred = false;
    ClwCtx<T> cx;
    cx.N = &N;
    cx.l = &l;
    cx.r = &r;
    cx.plan = &plan;
    cx.opt = &opt;
    cx.pi = T(num_traits<T>::from_int(2) * acos(zero));
    return clw_invert_d(0, w, cx, [&](const std::vector<Cx<T>>& wv, long bucket) {
        return (bucket < 0) ? Cx<T>(one, zero) : gbar(wv);
    });
}

}  // namespace detail

/**
 * @param L   (q' x p) single-server relative traffic intensities, L(i,j) = rho_{ji}
 * @param N   (p) closed-chain population vector
 * @param Z   (p) aggregate infinite-server relative intensities rho_{j0}
 * @param m   (q') queue multiplicities; empty for all ones
 * @param opt lattice and aliasing parameters
 */
template <class T>
ClwResult<T> pfqn_clw(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                      const std::vector<long>& m, const ClwOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_clw requires transcendental arithmetic (contour integration of a "
                  "generating function)");
    using std::exp;
    using std::log;
    const std::size_t qd = L.rows();
    if (L.cols() != N.size()) throw InputError("pfqn_clw: L and N disagree on the chain count");
    if (Z.size() != N.size()) throw InputError("pfqn_clw: Z and N disagree on the chain count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    ClwResult<T> res;
    for (std::size_t j = 0; j < N.size(); ++j)
        if (N[j] < 0) {
            res.G = zero;
            res.lG = T(-std::numeric_limits<T>::infinity());
            return res;
        }
    bool allzero = true;
    for (std::size_t j = 0; j < N.size(); ++j)
        if (N[j] > 0) allzero = false;
    if (allzero) {
        res.G = one;
        res.lG = zero;
        return res;
    }

    // zero-population chain drop rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<std::size_t> keep;
    for (std::size_t j = 0; j < N.size(); ++j)
        if (N[j] > 0) keep.push_back(j);
    const std::size_t p = keep.size();
    Matrix<T> Lk(qd, p, zero);
    std::vector<int> Nk(p, 0), l(p, 1);
    std::vector<T> Zk(p, zero);
    std::vector<double> gam(p, 0.0);
    for (std::size_t j = 0; j < p; ++j) {
        for (std::size_t i = 0; i < qd; ++i) Lk(i, j) = L(i, keep[j]);
        Nk[j] = N[keep[j]];
        Zk[j] = Z[keep[j]];
    }

    // dimension reduction (Section 3): D is inverted first, then each connected
    // component of the interdependence graph minus D, independently
    const detail::ClwPlan plan = detail::clw_plan(Lk, opt);
    std::vector<std::size_t> order, depth(p, 0);
    for (std::size_t t = 0; t < plan.D.size(); ++t) {
        order.push_back(plan.D[t]);
        depth[plan.D[t]] = t + 1;
    }
    for (std::size_t c = 0; c < plan.comps.size(); ++c)
        for (std::size_t s = 0; s < plan.comps[c].size(); ++s) {
            order.push_back(plan.comps[c][s]);
            depth[plan.comps[c][s]] = plan.D.size() + s + 1;
        }
    detail::clw_defaults(N.size(), keep, depth, opt, l, gam);
    std::vector<long> mult(qd, 1);
    if (!m.empty()) {
        if (m.size() != qd) throw InputError("pfqn_clw: m must have one entry per queue");
        for (std::size_t i = 0; i < qd; ++i) {
            if (m[i] < 1) throw InputError("pfqn_clw: the queue multiplicities must be positive");
            mult[i] = m[i];
        }
    }

    // contour radii r_j = 10^{-gamma_j / (2 l_j K_j)} (eq. 2.7)
    std::vector<T> r(p, one);
    for (std::size_t j = 0; j < p; ++j)
        r[j] = detail::clw_pow_real(
            num_traits<T>::from_int(10),
            T(num_traits<T>::from_double(-gam[j]) /
              num_traits<T>::from_int(2 * static_cast<long>(l[j]) * Nk[j])));

    const std::vector<T> alpha =
        detail::clw_scaling(Lk, Lk, Nk, Zk, l, r, mult, order, detail::clw_beta<T>(keep, opt));

    std::vector<T> arho0(p, zero);
    Matrix<T> rhoS(qd, p, zero);
    for (std::size_t j = 0; j < p; ++j) {
        arho0[j] = alpha[j] * Zk[j];
        for (std::size_t i = 0; i < qd; ++i) rhoS(i, j) = Lk(i, j) * alpha[j];
    }

    // split the factors of (4.5) over D and the components: a queue whose chains
    // all lie in D is constant during the component inversions, and every other
    // queue has all of its non-D chains inside a single component
    std::vector<char> inD(p, 0);
    for (std::size_t t = 0; t < plan.D.size(); ++t) inD[plan.D[t]] = 1;
    std::vector<long> compOf(p, 0);
    for (std::size_t c = 0; c < plan.comps.size(); ++c)
        for (std::size_t s = 0; s < plan.comps[c].size(); ++s)
            compOf[plan.comps[c][s]] = static_cast<long>(c);
    std::vector<long> qBucket(qd, -1);
    for (std::size_t i = 0; i < qd; ++i)
        for (std::size_t j = 0; j < p; ++j)
            if (Lk(i, j) != zero && !inD[j]) {
                qBucket[i] = compOf[j];
                break;
            }

    // per-group normalization (Section 2.2, page 944): the scaling normalizes the
    // whole generating function, not each group, and a decomposition multiplies
    // the groups together. Every factor has nonnegative coefficients, so the
    // group modulus is maximized at w = r; that constant cancels in the recovery
    // and is left at zero on the undecomposed path, which is thus unchanged.
    const bool decomposed = !plan.D.empty() || plan.comps.size() > 1;
    std::vector<T> off(plan.comps.size() + 1, zero);   // [0] is the D group
    if (decomposed) {
        for (std::size_t b = 0; b <= plan.comps.size(); ++b) {
            const long tag = static_cast<long>(b) - 1;
            T o = zero;
            const std::vector<std::size_t>& ch =
                (tag < 0) ? plan.D : plan.comps[static_cast<std::size_t>(tag)];
            for (std::size_t t = 0; t < ch.size(); ++t) o += arho0[ch[t]] * T(r[ch[t]] - one);
            for (std::size_t i = 0; i < qd; ++i) {
                if (qBucket[i] != tag) continue;
                T x = zero;
                for (std::size_t j = 0; j < p; ++j) x += rhoS(i, j) * r[j];
                T pole = T(one - x);
                if (!(pole > zero)) pole = num_traits<T>::from_double(std::numeric_limits<double>::min());
                o -= num_traits<T>::from_int(mult[i]) * log(pole);
            }
            off[b] = o;
        }
    }

    const auto gbar = [&](const std::vector<detail::Cx<T>>& w, long bucket) {
        const std::vector<std::size_t>& ch =
            (bucket < 0) ? plan.D : plan.comps[static_cast<std::size_t>(bucket)];
        detail::Cx<T> expo(zero, zero);
        for (std::size_t t = 0; t < ch.size(); ++t) {
            const std::size_t j = ch[t];
            expo = detail::cx_add(
                expo, detail::cx_scale(detail::Cx<T>(T(w[j].re - one), w[j].im), arho0[j]));
        }
        detail::Cx<T> logden(zero, zero);
        for (std::size_t i = 0; i < qd; ++i) {
            if (qBucket[i] != bucket) continue;
            detail::Cx<T> a(zero, zero);
            for (std::size_t j = 0; j < p; ++j)
                a = detail::cx_add(a, detail::cx_scale(w[j], rhoS(i, j)));
            const detail::Cx<T> lg = detail::cx_log(detail::Cx<T>(T(one - a.re), T(-a.im)));
            logden = detail::cx_add(
                logden, detail::cx_scale(lg, num_traits<T>::from_int(mult[i])));
        }
        const T o = off[static_cast<std::size_t>(bucket + 1)];
        return detail::cx_exp(detail::Cx<T>(T(expo.re - logden.re - o), T(expo.im - logden.im)));
    };

    detail::ClwCtx<T> cx;
    cx.N = &Nk;
    cx.l = &l;
    cx.r = &r;
    cx.plan = &plan;
    cx.opt = &opt;
    {
        using std::acos;
        cx.pi = T(num_traits<T>::from_int(2) * acos(zero));
    }
    std::vector<detail::Cx<T>> w(p, detail::Cx<T>(zero, zero));
    const detail::Cx<T> gv = detail::clw_invert_d(0, w, cx, gbar);
    if (!(gv.re > zero))
        throw NumericError("pfqn_clw: the inverted generating function is not positive");

    T lG = log(gv.re);
    for (std::size_t b = 0; b < off.size(); ++b) lG += off[b];
    for (std::size_t j = 0; j < p; ++j)
        lG += arho0[j] - num_traits<T>::from_int(Nk[j]) * log(alpha[j]);
    res.lG = lG;
    res.G = (lG > num_traits<T>::from_int(709)) ? T(std::numeric_limits<T>::infinity()) : T(exp(lG));
    return res;
}

/** Overload with unit multiplicities and the CLW default parameters. */
template <class T>
ClwResult<T> pfqn_clw(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_clw(L, N, Z, std::vector<long>(), ClwOptions());
}

/** Overload with the CLW default parameters. */
template <class T>
ClwResult<T> pfqn_clw(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                      const std::vector<long>& m) {
    return pfqn_clw(L, N, Z, m, ClwOptions());
}

/**
 * Limited load-dependent form (matlab pfqn_clw_lld.m).
 *
 * @param L   (q' x p) relative traffic intensities
 * @param N   (p) populations
 * @param Z   (p) infinite-server intensities
 * @param mu  (q' x n) load-dependent rate scalings S_i(k); the last column is
 *            extended when fewer than sum(N) are supplied (the LLD assumption),
 *            and empty means all queues are load independent
 * @param opt lattice and aliasing parameters
 */
template <class T>
ClwResult<T> pfqn_clw_lld(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                          const Matrix<T>& mu, const ClwOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_clw_lld requires transcendental arithmetic (contour integration of a "
                  "generating function)");
    using std::exp;
    using std::log;
    const std::size_t qd = L.rows();
    if (L.cols() != N.size()) throw InputError("pfqn_clw_lld: L and N disagree on the chain count");
    if (Z.size() != N.size()) throw InputError("pfqn_clw_lld: Z and N disagree on the chain count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    ClwResult<T> res;
    long Ntot = 0;
    for (std::size_t j = 0; j < N.size(); ++j) {
        if (N[j] < 0) {
            res.G = zero;
            res.lG = T(-std::numeric_limits<T>::infinity());
            return res;
        }
        Ntot += N[j];
    }
    if (Ntot == 0) {
        res.G = one;
        res.lG = zero;
        return res;
    }

    // extend or truncate mu to sum(N) columns (LLD extension of the last column)
    Matrix<T> S(qd, static_cast<std::size_t>(Ntot), one);
    if (!mu.empty()) {
        if (mu.rows() != qd) throw InputError("pfqn_clw_lld: mu must have one row per queue");
        for (std::size_t i = 0; i < qd; ++i)
            for (long k = 0; k < Ntot; ++k) {
                const std::size_t src = (static_cast<std::size_t>(k) < mu.cols())
                                            ? static_cast<std::size_t>(k)
                                            : mu.cols() - 1;
                if (!(mu(i, src) > zero))
                    throw InputError("pfqn_clw_lld: the load-dependent rates must be positive");
                S(i, static_cast<std::size_t>(k)) = mu(i, src);
            }
    }

    std::vector<std::size_t> keep;
    for (std::size_t j = 0; j < N.size(); ++j)
        if (N[j] > 0) keep.push_back(j);
    const std::size_t p = keep.size();
    Matrix<T> Lk(qd, p, zero);
    std::vector<int> Nk(p, 0), l(p, 1);
    std::vector<T> Zk(p, zero);
    std::vector<double> gam(p, 0.0);
    std::vector<std::size_t> order(p, 0), depth(p, 0);
    for (std::size_t j = 0; j < p; ++j) {
        for (std::size_t i = 0; i < qd; ++i) Lk(i, j) = L(i, keep[j]);
        Nk[j] = N[keep[j]];
        Zk[j] = Z[keep[j]];
        order[j] = j;
        depth[j] = j + 1;
    }
    // the accelerations are wired for pfqn_clw only: extending them to the LLD
    // form means porting the same change to all four codebases, so this routine
    // keeps the plain nested inversion and its numerics
    ClwOptions lopt = opt;
    lopt.euler = false;
    lopt.dimred = false;
    detail::clw_defaults(N.size(), keep, depth, lopt, l, gam);
    detail::ClwPlan lplan;
    lplan.comps.assign(1, order);

    // pole c_i and LLD cutoff l_i: S_i(k) = c_i for k >= l_i
    std::vector<T> cpole(qd, one);
    std::vector<std::vector<T>> numc(qd);
    for (std::size_t i = 0; i < qd; ++i) {
        cpole[i] = S(i, static_cast<std::size_t>(Ntot) - 1);
        long last = -1;
        for (long k = 0; k < Ntot; ++k)
            if (S(i, static_cast<std::size_t>(k)) != cpole[i]) last = k;
        const long li = (last < 0) ? 1 : last + 2;  // MATLAB last is 1-based
        std::vector<T> a(static_cast<std::size_t>(li), zero);
        a[0] = cpole[i];
        T cp = one;
        for (long n = 1; n < li; ++n) {
            cp *= S(i, static_cast<std::size_t>(n) - 1);
            a[static_cast<std::size_t>(n)] = (cpole[i] - S(i, static_cast<std::size_t>(n) - 1)) / cp;
        }
        numc[i] = a;
    }

    std::vector<T> r(p, one);
    for (std::size_t j = 0; j < p; ++j)
        r[j] = detail::clw_pow_real(
            num_traits<T>::from_int(10),
            T(num_traits<T>::from_double(-gam[j]) /
              num_traits<T>::from_int(2 * static_cast<long>(l[j]) * Nk[j])));

    // unit-pole intensities: each F_i behaves as a simple pole at 1
    Matrix<T> Lt(qd, p, zero);
    for (std::size_t i = 0; i < qd; ++i)
        for (std::size_t j = 0; j < p; ++j) Lt(i, j) = Lk(i, j) / cpole[i];
    const std::vector<long> mult(qd, 1);
    const std::vector<T> alpha =
        detail::clw_scaling(Lt, Lk, Nk, Zk, l, r, mult, order, detail::clw_beta<T>(keep, lopt));

    std::vector<T> arho0(p, zero);
    Matrix<T> rhoS(qd, p, zero);
    for (std::size_t j = 0; j < p; ++j) {
        arho0[j] = alpha[j] * Zk[j];
        for (std::size_t i = 0; i < qd; ++i) rhoS(i, j) = Lk(i, j) * alpha[j];
    }

    const auto gbar = [&](const std::vector<detail::Cx<T>>& w, long bucket) {
        if (bucket < 0) return detail::Cx<T>(one, zero);   // no D group in the trivial plan
        detail::Cx<T> expo(zero, zero);
        for (std::size_t j = 0; j < p; ++j)
            expo = detail::cx_add(
                expo, detail::cx_scale(detail::Cx<T>(T(w[j].re - one), w[j].im), arho0[j]));
        detail::Cx<T> logF(zero, zero);
        for (std::size_t i = 0; i < qd; ++i) {
            detail::Cx<T> x(zero, zero);
            for (std::size_t j = 0; j < p; ++j)
                x = detail::cx_add(x, detail::cx_scale(w[j], rhoS(i, j)));
            const std::vector<T>& a = numc[i];
            detail::Cx<T> num(a.back(), zero);  // Horner on N_i(x)
            for (std::size_t k = a.size() - 1; k-- > 0;)
                num = detail::cx_add(detail::cx_mul(num, x), detail::Cx<T>(a[k], zero));
            logF = detail::cx_add(logF, detail::cx_log(num));
            logF = detail::cx_sub(
                logF, detail::cx_log(detail::Cx<T>(T(cpole[i] - x.re), T(-x.im))));
        }
        return detail::cx_exp(detail::cx_add(expo, logF));
    };

    detail::ClwCtx<T> cx;
    cx.N = &Nk;
    cx.l = &l;
    cx.r = &r;
    cx.plan = &lplan;
    cx.opt = &lopt;
    {
        using std::acos;
        cx.pi = T(num_traits<T>::from_int(2) * acos(zero));
    }
    std::vector<detail::Cx<T>> w(p, detail::Cx<T>(zero, zero));
    const detail::Cx<T> gv = detail::clw_invert_d(0, w, cx, gbar);
    if (!(gv.re > zero))
        throw NumericError("pfqn_clw_lld: the inverted generating function is not positive");

    T lG = log(gv.re);
    for (std::size_t j = 0; j < p; ++j)
        lG += arho0[j] - num_traits<T>::from_int(Nk[j]) * log(alpha[j]);
    res.lG = lG;
    res.G = (lG > num_traits<T>::from_int(709)) ? T(std::numeric_limits<T>::infinity()) : T(exp(lG));
    return res;
}

/** Overload with the CLW default parameters. */
template <class T>
ClwResult<T> pfqn_clw_lld(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                          const Matrix<T>& mu) {
    return pfqn_clw_lld(L, N, Z, mu, ClwOptions());
}

/** Overload with all queues load independent. */
template <class T>
ClwResult<T> pfqn_clw_lld(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_clw_lld(L, N, Z, Matrix<T>(), ClwOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CLW_H
