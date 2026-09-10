/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_RESPT_PS_MOMENTS_H
#define LINE_API_PFQN_RESPT_PS_MOMENTS_H

/**
 * Sojourn-time moments at the processor-sharing station of a closed
 * terminal-driven system (Mitra and Morrison 1983).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_respt_ps_moments.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/Pfqn_respt_ps_moments.java.
 *
 * The model is a bank of terminals in series with a single processor-sharing
 * CPU, with class-dependent exponential think times (mean Z_r) and
 * class-dependent exponential service times (mean S_r), and N_r jobs of class r
 * cycling between the two. Two routes to the moments are implemented, both from
 * that paper:
 *
 *   Exact       solves the linear system c'[A - q_J I] = -pi'B of Proposition 3
 *               on the state space {n : 0 <= n <= K}, K being the population
 *               vector with the tagged class decremented by one. The moments are
 *               then E[W_J] = sum_n c(n) and (q_J/2) E[W_J^2] =
 *               sum_n (n'1 + 1) c(n). Exact to solver precision, at the cost of
 *               a linear solve of dimension prod_r (K_r + 1).
 *
 *   Asymptotic  evaluates the two leading terms of the expansion in inverse
 *               powers of the large parameter Nexp = max_r Z_r / S_r,
 *               E[W_J^2] ~ c0 + c1/Nexp, of Proposition 6. The cost is a linear
 *               system of dimension R and is therefore independent of the
 *               populations. NOTE the expansion parameter is the
 *               THINK-TO-SERVICE RATIO and NOT the population, so a model with
 *               short think times is expanded in a small parameter no matter how
 *               many jobs it holds.
 *
 *   Auto        (default) takes the exact route when the state space has at most
 *               4096 states and the asymptotic route otherwise.
 *
 * The asymptotic route requires the normal-usage condition alpha > 0, where
 * alpha = 1 - sum_r lambda_r / q_r with lambda_r = K_r / Z_r and q_r = 1 / S_r,
 * is the unutilized fraction of the CPU in the corresponding open system. Where
 * it fails and the exact route is not affordable, the entry of W and W2 is NaN
 * and the per-class method records Unavailable; asking for Asymptotic explicitly
 * in that regime is an error rather than a blank.
 *
 * A class with N_r = 0 has no sojourn time and its entries are NaN.
 *
 * Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of the
 * Waiting Time in Closed and Open Processor-Sharing Systems with Multiple Job
 * Classes", Adv. Appl. Prob. 15(4), 1983, Propositions 3 and 6.
 *
 * Arithmetic: TRANSCENDENTAL. The exact route builds its stationary law in the
 * log domain so that large populations do not overflow.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Which route produced the moments of a given class. */
enum class ResptPsMethod {
    None,        ///< the class is unpopulated
    Exact,       ///< Proposition 3, the linear solve
    Asymptotic,  ///< Proposition 6, the two-term expansion
    Unavailable  ///< normal usage fails and the exact route was not affordable
};

/** Requested route. */
enum class ResptPsRoute { Auto, Exact, Asymptotic };

/** Sojourn-time moments at the PS station, per class. */
template <class T>
struct ResptPsMomentsResult {
    std::vector<T> W;   ///< (R) mean sojourn times; NaN for an unpopulated class
    std::vector<T> W2;  ///< (R) second moments
    std::vector<ResptPsMethod> method;
    std::vector<T> c0;      ///< (R) leading asymptotic term, NaN off that route
    std::vector<T> c1;      ///< (R) first correction, NaN off that route
    std::vector<T> alpha;   ///< (R) unutilized CPU fraction of the open system
    std::vector<double> nstates;  ///< (R) size of the exact state space
    T expansionParam;       ///< Nexp = max_r q_r / p_r
};

namespace detail {

/** State-space size below which Auto goes exact. */
constexpr double RESPT_PS_AUTO_MAX = 4096.0;
/** Hard bound on an explicitly requested exact solve. */
constexpr double RESPT_PS_EXACT_MAX = 65536.0;

/**
 * Proposition 3: the moments follow from c, the solution of
 * c'[A - q_J I] = -pi'B, with A the generator-like operator of equation (26) and
 * B the diagonal operator B(n,n) = n'1 + 1.
 */
template <class T>
void respt_ps_exact(const std::vector<T>& p, const std::vector<T>& q, const std::vector<long>& K,
                    std::size_t J, T& W, T& W2) {
    using std::exp;
    using std::log;
    const std::size_t R = K.size();
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::size_t> dims(R);
    std::size_t ns = 1;
    for (std::size_t j = 0; j < R; ++j) {
        dims[j] = static_cast<std::size_t>(K[j]) + 1;
        ns *= dims[j];
    }
    std::vector<std::size_t> stride(R, 1);
    for (std::size_t j = 1; j < R; ++j) stride[j] = stride[j - 1] * dims[j - 1];

    Matrix<long> states(ns, R);
    std::vector<long> tot(ns, 0);
    for (std::size_t lin = 0; lin < ns; ++lin) {
        std::size_t res = lin;
        for (std::size_t j = 0; j < R; ++j) {
            states(lin, j) = static_cast<long>(res % dims[j]);
            res /= dims[j];
            tot[lin] += states(lin, j);
        }
    }

    // stationary law (15), in logs so that large populations do not overflow
    std::vector<T> r(R);
    for (std::size_t j = 0; j < R; ++j) r[j] = T(p[j] / q[j]);
    const double ninf = -std::numeric_limits<double>::infinity();
    std::vector<double> logpi(ns, 0.0);
    for (std::size_t lin = 0; lin < ns; ++lin) {
        double v = num_traits<T>::to_double(
            num_lgamma<T>(num_traits<T>::from_int(tot[lin] + 1)));
        for (std::size_t j = 0; j < R; ++j) {
            const long nj = states(lin, j);
            v += num_traits<T>::to_double(num_lgamma<T>(num_traits<T>::from_int(K[j] + 1))) -
                 num_traits<T>::to_double(num_lgamma<T>(num_traits<T>::from_int(nj + 1))) -
                 num_traits<T>::to_double(num_lgamma<T>(num_traits<T>::from_int(K[j] - nj + 1)));
            if (r[j] > num_traits<T>::from_int(0)) {
                v += static_cast<double>(nj) * num_traits<T>::log_as_double(r[j]);
            } else if (nj > 0) {
                v = ninf;
            }
        }
        logpi[lin] = v;
    }
    double mx = ninf;
    for (std::size_t lin = 0; lin < ns; ++lin) mx = std::max(mx, logpi[lin]);
    std::vector<T> pin(ns, zero);
    T psum = zero;
    for (std::size_t lin = 0; lin < ns; ++lin) {
        pin[lin] = num_traits<T>::from_double(std::exp(logpi[lin] - mx));
        psum += pin[lin];
    }
    for (std::size_t lin = 0; lin < ns; ++lin) pin[lin] = T(pin[lin] / psum);

    Matrix<T> A(ns, ns, zero);
    for (std::size_t lin = 0; lin < ns; ++lin) {
        T diagv = zero;
        for (std::size_t j = 0; j < R; ++j) {
            const long nj = states(lin, j);
            if (nj >= 1)
                A(lin - stride[j], lin) +=
                    p[j] * num_traits<T>::from_int(K[j] - nj + 1) *
                    num_traits<T>::from_int(tot[lin]);
            if (nj <= K[j] - 1)
                A(lin + stride[j], lin) += num_traits<T>::from_int(nj + 1) * q[j];
            diagv -= p[j] * num_traits<T>::from_int(K[j] - nj) *
                         num_traits<T>::from_int(tot[lin] + 1) +
                     num_traits<T>::from_int(nj) * q[j];
        }
        A(lin, lin) += diagv;
    }

    // c = ((A - q_J I)')^{-1} (-(tot+1) pi)
    Matrix<T> Mt(ns, ns, zero);
    for (std::size_t i = 0; i < ns; ++i)
        for (std::size_t k = 0; k < ns; ++k) Mt(i, k) = A(k, i);
    for (std::size_t i = 0; i < ns; ++i) Mt(i, i) -= q[J];
    std::vector<T> rhs(ns, zero);
    for (std::size_t lin = 0; lin < ns; ++lin)
        rhs[lin] = T(-num_traits<T>::from_int(tot[lin] + 1) * pin[lin]);
    const std::vector<T> c = solve(Mt, rhs);

    W = zero;
    T acc = zero;
    for (std::size_t lin = 0; lin < ns; ++lin) {
        W += c[lin];
        acc += num_traits<T>::from_int(tot[lin] + 1) * c[lin];
    }
    W2 = T(T(num_traits<T>::from_int(2) / q[J]) * acc);
}

/**
 * Proposition 6: the two leading terms of the expansion in 1/Nexp. Equation
 * numbers are those of Mitra and Morrison (1983).
 */
template <class T>
void respt_ps_asymptotic(const std::vector<T>& p, const std::vector<T>& q,
                         const std::vector<long>& K, std::size_t J, T& W, T& W2, T& c0, T& c1) {
    const std::size_t R = K.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> lambda(R);
    for (std::size_t j = 0; j < R; ++j) lambda[j] = T(p[j] * num_traits<T>::from_int(K[j]));
    T alpha = one;
    for (std::size_t j = 0; j < R; ++j) alpha -= T(lambda[j] / q[j]);
    T Nexp = T(q[0] / p[0]);  // (50)
    for (std::size_t j = 1; j < R; ++j) {
        const T v = T(q[j] / p[j]);
        if (v > Nexp) Nexp = v;
    }
    std::vector<T> Gam(R), beta(R);
    for (std::size_t j = 0; j < R; ++j) {
        Gam[j] = T(Nexp * p[j] / q[j]);                             // (51)
        beta[j] = T(num_traits<T>::from_int(K[j]) / Nexp);          // (51)
    }
    const T qJ = q[J];

    T den = one;
    for (std::size_t j = 0; j < R; ++j) den -= T(lambda[j] / T(q[j] + qJ));
    T num = one;
    for (std::size_t j = 0; j < R; ++j)
        num -= T(lambda[j] * T(q[j] - qJ) / T(q[j] * T(q[j] + qJ)));
    const T alpha2 = T(alpha * alpha);
    const T F10 = T(T(-one / T(alpha2 * qJ)) * num / den);          // (110)
    c0 = T(T(-num_traits<T>::from_int(2) / qJ) * F10);

    T bg2 = zero;
    for (std::size_t j = 0; j < R; ++j) bg2 += beta[j] * Gam[j] * Gam[j];
    std::vector<T> f1(R);
    for (std::size_t j = 0; j < R; ++j)                             // (113iii)
        f1[j] = T(T(lambda[j] / T(q[j] + qJ)) *
                  T(F10 - T(num_traits<T>::from_int(2) / T(alpha2 * q[j]))));
    const T alpha3 = T(alpha2 * alpha);
    const T alpha4 = T(alpha3 * alpha);
    std::vector<T> S2j(R);
    for (std::size_t j = 0; j < R; ++j)                             // (113i)
        S2j[j] = T(T(num_traits<T>::from_int(6) / alpha4) *
                   T(alpha * beta[j] * Gam[j] * Gam[j] +
                     num_traits<T>::from_int(2) * bg2 * beta[j] * Gam[j]));
    Matrix<T> S2js(R, R, zero);                                     // (113ii)
    for (std::size_t j = 0; j < R; ++j)
        for (std::size_t s = 0; s < R; ++s)
            S2js(j, s) = T(T(num_traits<T>::from_int(3) / alpha3) * beta[j] * Gam[j] * beta[s] *
                           Gam[s]);

    Matrix<T> Amat(R, R, zero);
    for (std::size_t j = 0; j < R; ++j) Amat(j, j) = one;
    std::vector<T> rhs(R, zero);
    for (std::size_t j = 0; j < R; ++j) {
        for (std::size_t s = 0; s < R; ++s) {
            const T d = T(q[j] + q[s] + qJ);
            Amat(j, s) -= T(lambda[j] / d);
            Amat(j, j) -= T(lambda[s] / d);
            rhs[j] += T(S2js(j, s) / d);
        }
        rhs[j] -= f1[j];
    }
    const std::vector<T> F2 = solve(Amat, rhs);                     // (112)

    const T f10 = T(T(-num_traits<T>::from_int(3) / T(alpha3 * qJ)) * bg2);  // (98)
    T acc = zero;
    for (std::size_t j = 0; j < R; ++j)
        acc += T(T(num_traits<T>::from_int(2) * Gam[j] * q[j] * F2[j] + S2j[j]) / T(q[j] + qJ));
    const T F20 = T(T(acc - f10) / den);                            // (111)
    c1 = T(T(T(-num_traits<T>::from_int(2) / qJ) * F20) + T(c0 / alpha2 * bg2));  // (114ii)

    W = T(T(one / T(alpha * qJ)) *
          T(one - T(num_traits<T>::from_int(2) / Nexp * bg2 / alpha2)));  // (68)
    W2 = T(c0 + T(c1 / Nexp));                                            // (114i)
}

}  // namespace detail

/**
 * @param S     (R) mean service times at the PS station, positive
 * @param N     (R) populations, non-negative integers
 * @param Z     (R) mean think times, positive where N > 0
 * @param route which route to take
 */
template <class T>
ResptPsMomentsResult<T> pfqn_respt_ps_moments(const std::vector<T>& S, const std::vector<long>& N,
                                              const std::vector<T>& Z, ResptPsRoute route) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_respt_ps_moments builds its stationary law in the log domain and needs "
                  "transcendental arithmetic");
    const std::size_t R = S.size();
    const T zero = num_traits<T>::from_int(0);
    if (N.size() != R || Z.size() != R)
        throw InputError("pfqn_respt_ps_moments: S, N and Z must have the same number of classes");
    for (std::size_t r = 0; r < R; ++r)
        if (S[r] <= zero) throw InputError("pfqn_respt_ps_moments: S must be finite and positive");
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0)
            throw InputError("pfqn_respt_ps_moments: N must contain non-negative integers");

    std::vector<std::size_t> act;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0) act.push_back(r);
    for (std::size_t k = 0; k < act.size(); ++k)
        if (Z[act[k]] <= zero)
            throw InputError(
                "pfqn_respt_ps_moments: Z must be finite and positive for every populated class");

    const T nan = std::numeric_limits<T>::quiet_NaN();
    ResptPsMomentsResult<T> res;
    res.W.assign(R, nan);
    res.W2.assign(R, nan);
    res.method.assign(R, ResptPsMethod::None);
    res.c0.assign(R, nan);
    res.c1.assign(R, nan);
    res.alpha.assign(R, nan);
    res.nstates.assign(R, std::numeric_limits<double>::quiet_NaN());
    res.expansionParam = nan;
    if (act.empty()) return res;

    const std::size_t A = act.size();
    std::vector<T> qa(A), pa(A);
    for (std::size_t k = 0; k < A; ++k) {
        qa[k] = T(num_traits<T>::from_int(1) / S[act[k]]);
        pa[k] = T(num_traits<T>::from_int(1) / Z[act[k]]);
    }
    res.expansionParam = T(qa[0] / pa[0]);
    for (std::size_t k = 1; k < A; ++k) {
        const T v = T(qa[k] / pa[k]);
        if (v > res.expansionParam) res.expansionParam = v;
    }

    for (std::size_t jj = 0; jj < A; ++jj) {
        const std::size_t J = act[jj];
        std::vector<long> K(A);
        for (std::size_t k = 0; k < A; ++k) K[k] = N[act[k]];
        K[jj] -= 1;
        double ns = 1.0;
        for (std::size_t k = 0; k < A; ++k) ns *= static_cast<double>(K[k] + 1);
        T alpha = num_traits<T>::from_int(1);
        for (std::size_t k = 0; k < A; ++k)
            alpha -= T(T(pa[k] * num_traits<T>::from_int(K[k])) / qa[k]);
        res.alpha[J] = alpha;
        res.nstates[J] = ns;

        const bool useExact =
            route == ResptPsRoute::Exact ||
            (route == ResptPsRoute::Auto && ns <= detail::RESPT_PS_AUTO_MAX);
        if (useExact) {
            if (ns > detail::RESPT_PS_EXACT_MAX)
                throw InputError(
                    "pfqn_respt_ps_moments: the exact route needs a linear solve above the "
                    "supported dimension; use the asymptotic route");
            detail::respt_ps_exact(pa, qa, K, jj, res.W[J], res.W2[J]);
            res.method[J] = ResptPsMethod::Exact;
            continue;
        }
        if (alpha <= zero) {
            if (route == ResptPsRoute::Asymptotic)
                throw InputError(
                    "pfqn_respt_ps_moments: the asymptotic expansion needs normal usage alpha > 0, "
                    "which this model violates");
            res.method[J] = ResptPsMethod::Unavailable;
            continue;
        }
        detail::respt_ps_asymptotic(pa, qa, K, jj, res.W[J], res.W2[J], res.c0[J], res.c1[J]);
        res.method[J] = ResptPsMethod::Asymptotic;
    }
    return res;
}

/** MATLAB default: the automatic route. */
template <class T>
ResptPsMomentsResult<T> pfqn_respt_ps_moments(const std::vector<T>& S, const std::vector<long>& N,
                                              const std::vector<T>& Z) {
    return pfqn_respt_ps_moments(S, N, Z, ResptPsRoute::Auto);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_RESPT_PS_MOMENTS_H
