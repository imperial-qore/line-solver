/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BUSYP_H
#define LINE_API_PFQN_PFQN_BUSYP_H

/**
 * Mean busy period of order n for a subnetwork of a product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_busyp.m,
 * jar/src/main/java/jline/api/pfqn/Pfqn_busyp.java and
 * python/line_solver/api/pfqn/busyp.py. Implements H. Daduna, "Busy Periods for
 * Subnetworks in Stochastic Networks: Mean Value Analysis", J. ACM 35(3), 1988:
 * Theorem 1 for a closed Gordon-Newell network and Theorem 3 for an open
 * Jackson network.
 *
 * The busy period of order n for a set of nodes I is the interval from the
 * instant a job entering I finds n-1 jobs in it until fewer than n remain. With
 * G(m,I) the normalizing constant of I at population m, H(m,I) that of the
 * complement, and A(I) the total rate at which jobs enter I from outside it,
 *
 *   closed:  b(n,I) = sum_{m=n}^{N} G(m,I) H(N-m,I) / [G(n-1,I) H(N-n,I) A(I)]
 *   open:    b(n,I) = sum_{m>=n}   G(m,I)           / [G(n-1,I) A(I)]
 *
 * The paper is single-chain: alpha solves x*P = x for a closed network and
 * x = gamma + x*P for an open one, and every node is a state-dependent
 * single-server FCFS station. By the insensitivity noted in Section 5 the
 * result depends on the service processes only through the rates mu.
 *
 * ARITHMETIC: the sums are accumulated in the log domain in double, exactly as
 * the three reference implementations do, which is what keeps G(m,I) from
 * overflowing on its own well before the ratio does. The templated inputs are
 * therefore read through num_traits<T>::to_double; there is no exact-rational
 * path, since a logarithm has none.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Default relative tolerance of the open-network tail truncation. */
inline constexpr double PFQN_BUSYP_DEFAULT_TOL = 1e-12;

namespace detail {

/** log-sum-exp, stable when every entry is -infinity. */
inline double busyp_lse(const std::vector<double>& v) {
    double m = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] > m) m = v[i];
    if (!std::isfinite(m)) return m;
    double s = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) s += std::exp(v[i] - m);
    return m + std::log(s);
}

/**
 * Rates of the selected nodes for populations 1..K. A rate table shorter than K
 * keeps its last rate, the saturated-server convention; a node whose rate does
 * not saturate (an infinite server) must be supplied through the callable
 * overload instead.
 */
template <class T>
std::vector<std::vector<double>> busyp_rates(const Matrix<T>& mu,
                                             const std::vector<std::size_t>& idx,
                                             std::size_t K) {
    std::vector<std::vector<double>> out(idx.size(), std::vector<double>(K, 0.0));
    const std::size_t cols = mu.cols();
    for (std::size_t i = 0; i < idx.size(); ++i)
        for (std::size_t k = 0; k < K; ++k)
            out[i][k] = num_traits<T>::to_double(mu(idx[i], std::min(k, cols - 1)));
    return out;
}

/** Same, for a rate law given as a callable rate(node, population). */
inline std::vector<std::vector<double>> busyp_rates(
    const std::function<double(std::size_t, std::size_t)>& mu,
    const std::vector<std::size_t>& idx, std::size_t K) {
    std::vector<std::vector<double>> out(idx.size(), std::vector<double>(K, 0.0));
    for (std::size_t i = 0; i < idx.size(); ++i)
        for (std::size_t k = 1; k <= K; ++k) out[i][k - 1] = mu(idx[i], k);
    return out;
}

/**
 * Log normalizing constants of orders 0..K of a set of nodes: lg[m] is the log
 * of the sum over the compositions n_1+...+n_L = m of the product over the
 * nodes of prod_{k=1}^{n_i} alpha_i/mu_i(k), the G(m,I) and H(m,I) of the
 * paper. The nodes are convolved one at a time in the log domain.
 */
inline std::vector<double> busyp_lgvec(const std::vector<double>& alpha,
                                       const std::vector<std::vector<double>>& mu,
                                       std::size_t K) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    std::vector<double> lg(K + 1, neg_inf);
    lg[0] = 0.0;
    for (std::size_t i = 0; i < alpha.size(); ++i) {
        std::vector<double> li(K + 1, 0.0);
        double acc = 0.0;
        for (std::size_t m = 1; m <= K; ++m) {
            acc += std::log(alpha[i]) - std::log(mu[i][m - 1]);
            li[m] = acc;
        }
        std::vector<double> lgnew(K + 1, neg_inf);
        std::vector<double> terms;
        for (std::size_t m = 0; m <= K; ++m) {
            terms.assign(m + 1, neg_inf);
            for (std::size_t k = 0; k <= m; ++k) terms[k] = lg[m - k] + li[k];
            lgnew[m] = busyp_lse(terms);
        }
        lg.swap(lgnew);
    }
    return lg;
}

/**
 * Truncation order of the open-network sum sum_{m>=n} G(m,I). The subnetwork
 * terms decay geometrically once the rates saturate, so the truncation grows
 * until the geometric tail estimate is negligible against the partial sum.
 */
template <class RateSource>
std::size_t busyp_trunc(const std::vector<double>& alpha, const RateSource& mu,
                        const std::vector<std::size_t>& subnet, std::size_t nmax,
                        double tol) {
    std::size_t K = std::max<std::size_t>(nmax + 8, 16);
    std::vector<std::vector<double>> rows = busyp_rates(mu, subnet, K);
    double rho_max = 0.0;
    for (std::size_t i = 0; i < alpha.size(); ++i)
        rho_max = std::max(rho_max, alpha[i] / rows[i][K - 1]);
    if (rho_max >= 1)
        throw InputError("pfqn_busyp: the subnetwork is not stable, its busy period is infinite");
    while (true) {
        const std::vector<double> lg = busyp_lgvec(alpha, busyp_rates(mu, subnet, K), K);
        // decay rate read off the last two orders, the exact ratio for a
        // saturated single-server subnetwork and an upper estimate otherwise
        double r = std::exp(lg[K] - lg[K - 1]);
        if (!(r < 1)) r = rho_max;
        const double ltail = lg[K] + std::log(r) - std::log1p(-r);
        std::vector<double> partial(lg.begin() + static_cast<std::ptrdiff_t>(nmax), lg.end());
        if (ltail - busyp_lse(partial) < std::log(tol)) return K;
        K = 2 * K;
        if (K > 1000000)
            throw NumericError(
                "pfqn_busyp: the open busy period sum did not converge, the subnetwork "
                "is nearly saturated");
    }
}

}  // namespace detail

/** What pfqn_busyp returns: the durations and the two constant sequences. */
struct BusyPeriodResult {
    std::vector<double> b;   ///< mean duration per requested order
    std::vector<double> lG;  ///< log normalizing constants of the subnetwork
    std::vector<double> lH;  ///< log constants of the complement, empty when open
};

/**
 * Mean busy period of order n for the subnetwork.
 *
 * @param alpha  relative arrival rates, one per node
 * @param mu     load-dependent rates, either a (J x K) matrix mu(j,k-1) with k
 *               jobs at node j, or a callable mu(j, k) when the rates do not
 *               saturate (an infinite server)
 * @param P      (J x J) routing matrix
 * @param N      population, infinity for an open network
 * @param subnet zero-based node indexes forming the subnetwork
 * @param n      busy period orders, 1 <= n <= N
 * @param gamma  external arrival rates, empty for a closed network
 * @param tol    relative tolerance of the open-network tail truncation
 */
template <class T, class RateSource>
BusyPeriodResult pfqn_busyp(const std::vector<double>& alpha, const RateSource& mu,
                            const Matrix<T>& P, double N,
                            const std::vector<std::size_t>& subnet,
                            const std::vector<std::size_t>& n,
                            const std::vector<double>& gamma = {},
                            double tol = PFQN_BUSYP_DEFAULT_TOL) {
    const std::size_t J = alpha.size();
    const bool is_closed = std::isfinite(N);

    std::vector<std::size_t> target = subnet;
    std::sort(target.begin(), target.end());
    target.erase(std::unique(target.begin(), target.end()), target.end());
    if (target.empty()) throw InputError("pfqn_busyp: the subnetwork must be non-empty");
    if (is_closed && target.size() >= J)
        // a closed network needs jobs outside the subnetwork to start a busy period
        throw InputError(
            "pfqn_busyp: in a closed network the subnetwork must be a proper subset of "
            "the nodes");
    if (target.back() >= J)
        throw InputError("pfqn_busyp: the subnetwork indexes are out of range");

    std::vector<bool> in_subnet(J, false);
    for (std::size_t i = 0; i < target.size(); ++i) in_subnet[target[i]] = true;
    std::vector<std::size_t> compl_nodes;
    for (std::size_t j = 0; j < J; ++j)
        if (!in_subnet[j]) compl_nodes.push_back(j);

    std::size_t nmax = 0;
    for (std::size_t t = 0; t < n.size(); ++t) {
        if (n[t] < 1)
            throw InputError("pfqn_busyp: the busy period order must be a positive integer");
        if (is_closed && static_cast<double>(n[t]) > N)
            throw InputError("pfqn_busyp: the busy period order must be an integer in 1..N");
        nmax = std::max(nmax, n[t]);
    }
    if (!is_closed && gamma.empty())
        throw InputError("pfqn_busyp: an open network requires the external arrival rates gamma");

    // A(I) for a closed network, C(I) for an open one: both are the total rate at
    // which jobs enter the subnetwork from outside it, which is what starts a busy
    // period. The closed network has no external stream.
    double inflow = 0.0;
    for (std::size_t i = 0; i < compl_nodes.size(); ++i)
        for (std::size_t j = 0; j < target.size(); ++j)
            inflow += alpha[compl_nodes[i]] *
                      num_traits<T>::to_double(P(compl_nodes[i], target[j]));
    if (!gamma.empty())
        for (std::size_t j = 0; j < target.size(); ++j) inflow += gamma[target[j]];
    if (inflow <= 0)
        throw InputError(
            "pfqn_busyp: no job ever enters the subnetwork, its busy period is undefined");

    std::vector<double> alpha_sub, alpha_compl;
    for (std::size_t i = 0; i < target.size(); ++i) alpha_sub.push_back(alpha[target[i]]);
    for (std::size_t i = 0; i < compl_nodes.size(); ++i)
        alpha_compl.push_back(alpha[compl_nodes[i]]);

    BusyPeriodResult out;
    out.b.assign(n.size(), 0.0);
    if (is_closed) {
        const std::size_t pop = static_cast<std::size_t>(std::llround(N));
        out.lG = detail::busyp_lgvec(alpha_sub, detail::busyp_rates(mu, target, pop), pop);
        out.lH = detail::busyp_lgvec(alpha_compl,
                                     detail::busyp_rates(mu, compl_nodes, pop), pop);
        for (std::size_t t = 0; t < n.size(); ++t) {
            // Theorem 1
            std::vector<double> terms;
            for (std::size_t m = n[t]; m <= pop; ++m)
                terms.push_back(out.lG[m] + out.lH[pop - m]);
            out.b[t] = std::exp(detail::busyp_lse(terms) - out.lG[n[t] - 1] -
                                out.lH[pop - n[t]] - std::log(inflow));
        }
    } else {
        const std::size_t K = detail::busyp_trunc(alpha_sub, mu, target, nmax, tol);
        out.lG = detail::busyp_lgvec(alpha_sub, detail::busyp_rates(mu, target, K), K);
        for (std::size_t t = 0; t < n.size(); ++t) {
            // Theorem 3, the tail summed up to the truncation order
            std::vector<double> terms;
            for (std::size_t m = n[t]; m <= K; ++m) terms.push_back(out.lG[m]);
            out.b[t] = std::exp(detail::busyp_lse(terms) - out.lG[n[t] - 1] -
                                std::log(inflow));
        }
    }
    return out;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BUSYP_H
