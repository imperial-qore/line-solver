/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BUSYP_CLW_H
#define LINE_API_PFQN_PFQN_BUSYP_CLW_H

/**
 * Busy period of a subnetwork from point evaluations of the normalizing constant.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_busyp_clw.m,
 * jar/src/main/java/jline/api/pfqn/Pfqn_busyp_clw.java and
 * python/line_solver/api/pfqn/busyp_clw.py.
 *
 * WHAT THIS BUYS OVER `pfqn_busyp` / `pfqn_busyp_multiclass`. Those walk the
 * whole population ladder (the whole lattice, multichain) because the numerator
 * sums over {|m| >= n}. The complement of that set is the SHELLS |m| <= n-1, and
 * summing the product form over the WHOLE lattice is the full network's own
 * normalizing constant, G_I and H convolving to it:
 *
 *   sum_{m : |m| >= n} G_I(m) H(N-m) = G(N) - sum_{m : |m| <= n-1} G_I(m) H(N-m)
 *
 * so order n needs only the n lowest shells plus ONE evaluation of G(N). The
 * ordinary busy period n=1 collapses to three constants,
 *
 *   b(1,I) = [G(N) - H(N)] / sum_r A_r(I) H(N-e_r)
 *
 * all point evaluations at or near the full population, which is what the
 * normalizing-constant methods are built for. This routine calls CLW
 * (Choudhury-Leung-Whitt, J. ACM 42, 1995, numerical inversion of the generating
 * function); any method returning lG(N) can take its place. The cost stops
 * depending on N.
 *
 * THE OPEN CASE NEEDS NO INVERSION. The subnetwork's constant sequence has
 * generating function g(z) = prod_{i in I} f_i(z) and the tail is g(1) minus a
 * partial sum, with f_i(1) = 1/(1-rho_i) at a single server and exp(rho_i) at an
 * infinite one. That removes the tail TRUNCATION of the ladder routine, not just
 * its cost: the tail is exact, and the C++ test checks it at == against the
 * closed form.
 *
 * ACCURACY: the numerator is a difference of two nearly equal quantities when the
 * level set is unlikely, so the relative error grows with n -- 6e-12 at n=1
 * against 2.7e-08 at n=N on a three-station closed model at N=20. Cost grows with
 * n too, so the routine is most accurate where it is fastest.
 *
 * SCOPE: CLW's generating function covers single-server and infinite-server
 * stations, so a general load-dependent scaling belongs to
 * `pfqn_busyp_multiclass`. The identity is for the AGGREGATE level set: a
 * per-class one has complement {m_r <= n-1}, the whole lattice in the other
 * chains, which buys nothing.
 *
 * ARITHMETIC: log domain in double, as in the three reference implementations.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_busyp.h"
// busyp_factln and the station function live with the multichain routine
#include "line/api/pfqn/pfqn_busyp_multiclass.h"
#include "line/api/pfqn/pfqn_clw.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * log G(k) of a set of nodes: the single servers go to the normalizing-constant
 * method and the infinite servers into its aggregate think time.
 */
inline double busyp_clw_lognc(const std::vector<std::vector<double>>& L,
                              const std::vector<bool>& isdelay,
                              const std::vector<std::size_t>& nodes,
                              const std::vector<int>& k, const std::string& method) {
    const std::size_t R = k.size();
    bool all_zero = true;
    for (std::size_t r = 0; r < R; ++r) {
        if (k[r] < 0) return -std::numeric_limits<double>::infinity();
        if (k[r] != 0) all_zero = false;
    }
    if (all_zero) return 0.0;
    std::vector<std::size_t> queues;
    std::vector<double> Z(R, 0.0);
    for (std::size_t t = 0; t < nodes.size(); ++t) {
        if (isdelay[nodes[t]]) {
            for (std::size_t r = 0; r < R; ++r) Z[r] += L[nodes[t]][r];
        } else {
            queues.push_back(nodes[t]);
        }
    }
    if (queues.empty()) {
        // only infinite servers left: G(k) = prod_r Z_r^k_r / k_r!
        double out = 0.0;
        for (std::size_t r = 0; r < R; ++r) {
            if (k[r] == 0) continue;
            if (Z[r] <= 0) return -std::numeric_limits<double>::infinity();
            out += k[r] * std::log(Z[r]) - busyp_factln(static_cast<std::size_t>(k[r]));
        }
        return out;
    }
    if (method != "clw")
        throw InputError(
            "pfqn_busyp_clw: only the clw method is wired here; the point evaluation is a "
            "plug-in, so another token needs its own call rather than a silent substitution");
    Matrix<double> Lq(queues.size(), R, 0.0);
    for (std::size_t i = 0; i < queues.size(); ++i)
        for (std::size_t r = 0; r < R; ++r) Lq(i, r) = L[queues[i]][r];
    return static_cast<double>(pfqn_clw(Lq, k, Z).lG);
}

/** log X_i(k): multinomial at a single server, 1/prod k_r! at an infinite one. */
inline double busyp_clw_node(const std::vector<double>& Li, bool isdelay,
                             const std::vector<std::size_t>& k) {
    std::size_t tot = 0;
    for (std::size_t r = 0; r < k.size(); ++r) tot += k[r];
    if (tot == 0) return 0.0;
    double v = isdelay ? 0.0 : busyp_factln(tot);
    for (std::size_t r = 0; r < k.size(); ++r) {
        if (k[r] == 0) continue;
        if (Li[r] <= 0) return -std::numeric_limits<double>::infinity();
        v += -busyp_factln(k[r]) + static_cast<double>(k[r]) * std::log(Li[r]);
    }
    return v;
}

/**
 * log G_I(m) by direct enumeration of the splits of m across the subnetwork
 * nodes, cheap because m is bounded by n-1 and n is small wherever this wins.
 */
inline double busyp_clw_station(const std::vector<std::vector<double>>& L,
                                const std::vector<bool>& isdelay,
                                const std::vector<std::size_t>& nodes, std::size_t from,
                                const std::vector<std::size_t>& m) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const std::size_t R = m.size();
    if (from >= nodes.size()) {
        for (std::size_t r = 0; r < R; ++r)
            if (m[r] > 0) return neg_inf;
        return 0.0;
    }
    std::size_t total = 1;
    for (std::size_t r = 0; r < R; ++r) total *= m[r] + 1;
    std::vector<double> acc;
    std::vector<std::size_t> head(R, 0), rest(R, 0);
    for (std::size_t idx = 0; idx < total; ++idx) {
        std::size_t t = idx;
        for (std::size_t r = 0; r < R; ++r) {
            head[r] = t % (m[r] + 1);
            t /= m[r] + 1;
        }
        const double lterm = busyp_clw_node(L[nodes[from]], isdelay[nodes[from]], head);
        if (lterm == neg_inf) continue;
        for (std::size_t r = 0; r < R; ++r) rest[r] = m[r] - head[r];
        const double lrest = busyp_clw_station(L, isdelay, nodes, from + 1, rest);
        if (lrest == neg_inf) continue;
        acc.push_back(lterm + lrest);
    }
    return busyp_lse(acc);
}

/** log G_I(0..kmax) of an OPEN subnetwork, convolving the per-node series. */
inline std::vector<double> busyp_clw_open(const std::vector<double>& rho,
                                          const std::vector<bool>& isdelay,
                                          std::size_t kmax) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    std::vector<double> lg(kmax + 1, neg_inf);
    lg[0] = 0.0;
    for (std::size_t i = 0; i < rho.size(); ++i) {
        std::vector<double> li(kmax + 1, 0.0);
        double acc = 0.0;
        for (std::size_t k = 1; k <= kmax; ++k) {
            // 1/(1-rho z) has coefficients rho^k; exp(rho z) has rho^k/k!
            acc += std::log(rho[i]) - (isdelay[i] ? std::log(static_cast<double>(k)) : 0.0);
            li[k] = acc;
        }
        std::vector<double> lgnew(kmax + 1, neg_inf);
        for (std::size_t m = 0; m <= kmax; ++m) {
            std::vector<double> terms(m + 1, neg_inf);
            for (std::size_t k = 0; k <= m; ++k) terms[k] = lg[m - k] + li[k];
            lgnew[m] = busyp_lse(terms);
        }
        lg.swap(lgnew);
    }
    return lg;
}

}  // namespace detail

/**
 * Mean busy period of order n for the subnetwork, via NC point evaluations.
 *
 * @param alpha   (J x R) relative arrival rates, one column per chain
 * @param mu      (J x R) service rates, the chain-r rate at node j
 * @param P       routing matrices, one per chain (size 1 = shared by all chains)
 * @param N       population per chain, infinite entries for an open chain
 * @param subnet  zero-based node indexes forming the subnetwork
 * @param n       busy period orders, counting the jobs of every chain
 * @param gamma   (J x R) external arrival rates, empty for a closed network
 * @param isdelay infinite-server nodes, empty meaning all single servers
 * @param method  method name of the normalizing-constant method ("clw")
 */
template <class T>
std::vector<double> pfqn_busyp_clw(const Matrix<T>& alpha, const Matrix<T>& mu,
                                   const std::vector<Matrix<T>>& P,
                                   const std::vector<double>& N,
                                   const std::vector<std::size_t>& subnet,
                                   const std::vector<std::size_t>& n,
                                   const Matrix<T>& gamma = Matrix<T>(),
                                   const std::vector<bool>& isdelay = std::vector<bool>(),
                                   const std::string& method = "clw") {
    const std::size_t J = alpha.rows(), R = alpha.cols();
    bool is_closed = true, is_open = true;
    for (std::size_t r = 0; r < R; ++r) {
        if (std::isinf(N[r]))
            is_closed = false;
        else
            is_open = false;
    }
    if (!is_closed && !is_open)
        throw InputError(
            "pfqn_busyp_clw: a mixed model needs the lattice routine pfqn_busyp_multiclass");
    std::vector<bool> delay = isdelay.empty() ? std::vector<bool>(J, false) : isdelay;

    std::vector<std::size_t> target = subnet;
    std::sort(target.begin(), target.end());
    target.erase(std::unique(target.begin(), target.end()), target.end());
    if (target.empty()) throw InputError("pfqn_busyp_clw: the subnetwork must be non-empty");
    if (is_closed && target.size() >= J)
        throw InputError(
            "pfqn_busyp_clw: in a closed network the subnetwork must be a proper subset");
    std::vector<bool> in_subnet(J, false);
    for (std::size_t i = 0; i < target.size(); ++i) in_subnet[target[i]] = true;
    std::vector<std::size_t> compl_nodes, all_nodes;
    for (std::size_t j = 0; j < J; ++j) {
        all_nodes.push_back(j);
        if (!in_subnet[j]) compl_nodes.push_back(j);
    }

    // demands L(i,r) = alpha(i,r)/mu(i,r), zero where chain r does not visit node i
    std::vector<std::vector<double>> L(J, std::vector<double>(R, 0.0));
    for (std::size_t i = 0; i < J; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            const double a = num_traits<T>::to_double(alpha(i, r));
            const double m = num_traits<T>::to_double(mu(i, r));
            L[i][r] = (a > 0 && m > 0) ? a / m : 0.0;
        }

    // A_r(I): the chain-r rate at which jobs enter the subnetwork from outside it
    std::vector<double> A(R, 0.0);
    double inflow = 0.0;
    for (std::size_t r = 0; r < R; ++r) {
        const Matrix<T>& Pr = (P.size() == 1) ? P[0] : P[r];
        for (std::size_t i = 0; i < compl_nodes.size(); ++i)
            for (std::size_t j = 0; j < target.size(); ++j)
                A[r] += num_traits<T>::to_double(alpha(compl_nodes[i], r)) *
                        num_traits<T>::to_double(Pr(compl_nodes[i], target[j]));
        if (gamma.rows() == J)
            for (std::size_t j = 0; j < target.size(); ++j)
                A[r] += num_traits<T>::to_double(gamma(target[j], r));
        inflow += A[r];
    }
    if (inflow <= 0)
        throw InputError(
            "pfqn_busyp_clw: no job ever enters the subnetwork, its busy period is undefined");

    std::size_t nmax = 1;
    for (std::size_t t = 0; t < n.size(); ++t) nmax = std::max(nmax, n[t]);

    if (is_open) {
        // g_I(1) in closed form, so the tail is exact rather than truncated
        std::vector<double> rho(target.size(), 0.0);
        std::vector<bool> delayI(target.size(), false);
        double lg1 = 0.0;
        for (std::size_t t = 0; t < target.size(); ++t) {
            for (std::size_t r = 0; r < R; ++r) rho[t] += L[target[t]][r];
            delayI[t] = delay[target[t]];
            if (delayI[t]) {
                lg1 += rho[t];
            } else {
                if (rho[t] >= 1)
                    throw InputError(
                        "pfqn_busyp_clw: the subnetwork is not stable, its busy period is infinite");
                lg1 += -std::log1p(-rho[t]);
            }
        }
        const std::vector<double> lseq = detail::busyp_clw_open(rho, delayI, nmax);
        std::vector<double> b(n.size(), 0.0);
        for (std::size_t t = 0; t < n.size(); ++t) {
            const std::vector<double> head(lseq.begin(),
                                           lseq.begin() + static_cast<std::ptrdiff_t>(n[t]));
            const double tail = lg1 + std::log1p(-std::exp(detail::busyp_lse(head) - lg1));
            b[t] = std::exp(tail - lseq[n[t] - 1] - std::log(inflow));
        }
        return b;
    }

    std::vector<int> pop(R, 0);
    std::vector<std::size_t> bound(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        pop[r] = static_cast<int>(std::llround(N[r]));
        bound[r] = std::min<std::size_t>(static_cast<std::size_t>(pop[r]), nmax - 1);
    }
    std::vector<std::size_t> stride(R, 1);
    for (std::size_t r = 1; r < R; ++r) stride[r] = stride[r - 1] * (bound[r - 1] + 1);
    std::size_t size = 1;
    for (std::size_t r = 0; r < R; ++r) size *= bound[r] + 1;
    std::vector<std::vector<std::size_t>> mvec(size, std::vector<std::size_t>(R, 0));
    for (std::size_t idx = 0; idx < size; ++idx)
        for (std::size_t r = 0; r < R; ++r) mvec[idx][r] = (idx / stride[r]) % (bound[r] + 1);

    // the low shells of the subnetwork, the only lattice this routine walks
    std::vector<double> lGlow(size, 0.0);
    for (std::size_t idx = 0; idx < size; ++idx)
        lGlow[idx] = detail::busyp_clw_station(L, delay, target, 0, mvec[idx]);

    const double lGfull = detail::busyp_clw_lognc(L, delay, all_nodes, pop, method);

    std::vector<double> b(n.size(), 0.0);
    for (std::size_t t = 0; t < n.size(); ++t) {
        std::vector<double> corr, den;
        for (std::size_t idx = 0; idx < size; ++idx) {
            std::size_t level = 0;
            for (std::size_t r = 0; r < R; ++r) level += mvec[idx][r];
            if (level + 1 <= n[t]) {
                // numerator: the full constant minus the shells the level set excludes
                std::vector<int> left(R, 0);
                for (std::size_t r = 0; r < R; ++r)
                    left[r] = pop[r] - static_cast<int>(mvec[idx][r]);
                corr.push_back(lGlow[idx] +
                               detail::busyp_clw_lognc(L, delay, compl_nodes, left, method));
            }
            if (level + 1 != n[t]) continue;
            // denominator: the flow out of the shell |m| = n-1
            std::vector<double> terms;
            for (std::size_t r = 0; r < R; ++r) {
                if (A[r] <= 0) continue;
                std::vector<int> left(R, 0);
                bool ok = true;
                for (std::size_t s = 0; s < R; ++s) {
                    left[s] = pop[s] - static_cast<int>(mvec[idx][s]) - (s == r ? 1 : 0);
                    if (left[s] < 0) ok = false;
                }
                if (!ok) continue;
                terms.push_back(std::log(A[r]) +
                                detail::busyp_clw_lognc(L, delay, compl_nodes, left, method));
            }
            if (!terms.empty()) den.push_back(lGlow[idx] + detail::busyp_lse(terms));
        }
        const double num = lGfull + std::log1p(-std::exp(detail::busyp_lse(corr) - lGfull));
        b[t] = std::exp(num - detail::busyp_lse(den));
    }
    return b;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BUSYP_CLW_H
