/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BUSYP_MULTICLASS_H
#define LINE_API_PFQN_PFQN_BUSYP_MULTICLASS_H

/**
 * Multichain generalization of `pfqn_busyp`.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_busyp_multiclass.m,
 * jar/src/main/java/jline/api/pfqn/Pfqn_busyp_multiclass.java and
 * python/line_solver/api/pfqn/busyp_multiclass.py.
 *
 * Daduna (J. ACM 35(3), 1988) states Theorems 1 and 3 for a single chain and
 * notes in Section 5 that they carry over to the whole product-form class. The
 * proof uses only that the stationary law is product form and that the busy
 * period is Keilson's mean ergodic sojourn time on a level set, neither of which
 * is single-chain, so replacing the scalar population by a per-chain vector m
 * gives, for a closed network,
 *
 *             sum_{m : |m| >= n}   G_I(m) H(N-m)
 *   b(n,I) = --------------------------------------------------
 *             sum_{m : |m| = n-1}  G_I(m) sum_r A_r(I) H(N-m-e_r)
 *
 * with G_I and H the normalizing constants of the subnetwork and of its
 * complement at a population VECTOR and A_r(I) the chain-r arrival flow into I.
 * The denominator is the exact chain-r flow across the cut: a chain-r departure
 * from the complement at population k occurs at rate alpha_ir H(k-e_r)/H(k),
 * and the H(k) cancels the state weight. At R=1 the inner sum holds the single
 * term m=n-1 and H(N-m-e_1)=H(N-n), so it collapses to Theorem 1 exactly, which
 * is the regression the test runs.
 *
 * THE OPEN CASE NEEDS NO LATTICE. In an open product-form network the stations
 * are independent and the total occupancy of a node depends on the AGGREGATE
 * load sum_r alpha_ir/mu_ir alone, since summing the station function over the
 * compositions of t collapses the multinomial to (sum_r rho_ir)^t. It is
 * therefore reduced here to the single-chain routine on aggregated demands.
 *
 * PER CLASS: with jobclass = r the level set becomes {m_r >= n}, the jobs of
 * chain r alone. Only chain-r arrivals move that level, so the flow sum loses its
 * sum over r and the same two lattices serve every class.
 *
 * A MIXED MODEL keeps the closed lattice with its OPEN dimensions TRUNCATED. The
 * closed chains are conserved between the subnetwork and its complement, the open
 * ones are not: the complement's open count is free, so its open dimensions are
 * summed out and no e_r shift applies to an open chain, removing one job from an
 * unbounded dimension leaving the same sum. The truncation grows until the answer
 * stops moving, and is the only approximation in that branch. For a per-class
 * query on a CLOSED chain there is an exact shortcut: marginalizing the open
 * chains leaves a closed network with the demands deflated by 1/(1-rho_i^open),
 * where the RATES and not the visits must be deflated, since A_r is built from
 * the visit ratios.
 *
 * ARITHMETIC: log domain in double, as in the three reference implementations.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_busyp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** log(k!). */
inline double busyp_factln(std::size_t k) {
    return (k < 2) ? 0.0 : std::lgamma(static_cast<double>(k) + 1.0);
}

/**
 * log X_i(m) over the lattice for one node.
 *
 * X_i(m) = multinomial(|m|; m) prod_r L(i,r)^m_r / prod_{k=1}^{|m|} phi_i(k),
 * which at R=1 is the prod_k alpha_i/mu_i(k) of the single-chain routine and at
 * phi(k)=k the infinite-server form prod_r L^m_r/m_r!.
 */
inline std::vector<double> busyp_station(const std::vector<double>& Li,
                                         const std::vector<double>& phii,
                                         const std::vector<std::vector<std::size_t>>& mvec) {
    const std::size_t size = mvec.size(), R = mvec.empty() ? 0 : mvec[0].size();
    std::vector<double> out(size, 0.0);
    for (std::size_t idx = 0; idx < size; ++idx) {
        std::size_t tot = 0;
        for (std::size_t r = 0; r < R; ++r) tot += mvec[idx][r];
        double v = busyp_factln(tot);
        bool ok = true;
        for (std::size_t r = 0; r < R && ok; ++r) {
            if (mvec[idx][r] == 0) continue;
            if (Li[r] <= 0)
                ok = false;
            else
                v += -busyp_factln(mvec[idx][r]) +
                     static_cast<double>(mvec[idx][r]) * std::log(Li[r]);
        }
        if (!ok) {
            out[idx] = -std::numeric_limits<double>::infinity();
            continue;
        }
        for (std::size_t k = 1; k <= tot; ++k)
            v -= std::log(phii[std::min(k, phii.size()) - 1]);
        out[idx] = v;
    }
    return out;
}

/**
 * Log normalizing constants over the whole lattice of a set of nodes. A node
 * whose scaling row is all ones takes the Buzen recursion, O(R) per lattice
 * point; any other node needs the full sub-lattice convolution.
 */
inline std::vector<double> busyp_lgvec_multi(
    const std::vector<std::vector<double>>& L, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<std::size_t>>& mvec, const std::vector<std::size_t>& stride,
    const std::vector<std::size_t>& pop) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const std::size_t nodes = L.size(), R = pop.size(), size = mvec.size();
    std::size_t totalN = 0;
    for (std::size_t r = 0; r < R; ++r) totalN += pop[r];
    std::vector<double> lg(size, neg_inf);
    lg[0] = 0.0;
    for (std::size_t i = 0; i < nodes; ++i) {
        bool is_li = true;
        for (std::size_t k = 0; k < std::min(phi[i].size(), std::max<std::size_t>(1, totalN)); ++k)
            if (phi[i][k] != 1.0) is_li = false;
        if (is_li) {
            std::vector<double> lgnew = lg;
            for (std::size_t idx = 0; idx < size; ++idx) {
                double acc = lgnew[idx];
                for (std::size_t r = 0; r < R; ++r)
                    if (mvec[idx][r] > 0 && L[i][r] > 0) {
                        const double alt = std::log(L[i][r]) + lgnew[idx - stride[r]];
                        acc = busyp_lse(std::vector<double>{acc, alt});
                    }
                lgnew[idx] = acc;
            }
            lg.swap(lgnew);
        } else {
            const std::vector<double> lX = busyp_station(L[i], phi[i], mvec);
            std::vector<double> lgnew(size, neg_inf);
            for (std::size_t a = 0; a < size; ++a) {
                if (lg[a] == neg_inf) continue;
                for (std::size_t c = 0; c < size; ++c) {
                    if (lX[c] == neg_inf) continue;
                    bool fits = true;
                    std::size_t j = 0;
                    for (std::size_t r = 0; r < R && fits; ++r) {
                        const std::size_t s = mvec[a][r] + mvec[c][r];
                        if (s > pop[r])
                            fits = false;
                        else
                            j += s * stride[r];
                    }
                    if (fits) lgnew[j] = busyp_lse(std::vector<double>{lgnew[j], lg[a] + lX[c]});
                }
            }
            lg.swap(lgnew);
        }
    }
    return lg;
}

/**
 * The lattice evaluation shared by the closed and the mixed branch. `pop` bounds
 * every chain: the population of a closed one, the truncation of an open one;
 * `open_chain` names the dimensions that are NOT conserved, whose complement
 * counts are summed out rather than read at N-m.
 */
inline std::vector<double> busyp_lattice(const std::vector<std::vector<double>>& L,
                                         const std::vector<std::vector<double>>& scaling,
                                         const std::vector<std::size_t>& target,
                                         const std::vector<std::size_t>& compl_nodes,
                                         const std::vector<std::size_t>& pop,
                                         const std::vector<std::size_t>& n,
                                         const std::vector<double>& A, int jobclass,
                                         const std::vector<bool>& open_chain) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    const std::size_t R = pop.size();
    std::vector<std::size_t> stride(R, 1);
    for (std::size_t r = 1; r < R; ++r) stride[r] = stride[r - 1] * (pop[r - 1] + 1);
    std::size_t size = 1;
    for (std::size_t r = 0; r < R; ++r) size *= pop[r] + 1;
    std::vector<std::vector<std::size_t>> mvec(size, std::vector<std::size_t>(R, 0));
    for (std::size_t idx = 0; idx < size; ++idx)
        for (std::size_t r = 0; r < R; ++r) mvec[idx][r] = (idx / stride[r]) % (pop[r] + 1);

    std::vector<std::vector<double>> Lsub, Lcompl, Psub, Pcompl;
    for (std::size_t i = 0; i < target.size(); ++i) {
        Lsub.push_back(L[target[i]]);
        Psub.push_back(scaling[target[i]]);
    }
    for (std::size_t i = 0; i < compl_nodes.size(); ++i) {
        Lcompl.push_back(L[compl_nodes[i]]);
        Pcompl.push_back(scaling[compl_nodes[i]]);
    }
    const std::vector<double> lG = busyp_lgvec_multi(Lsub, Psub, mvec, stride, pop);
    const std::vector<double> lH = busyp_lgvec_multi(Lcompl, Pcompl, mvec, stride, pop);

    // Hbar sums the complement over its unconserved dimensions, so it is indexed by
    // the CLOSED components alone; with no open chain it is lH itself.
    bool any_open = false;
    for (std::size_t r = 0; r < R; ++r) any_open = any_open || open_chain[r];
    std::vector<double> lHbar = lH;
    if (any_open) {
        lHbar.assign(size, neg_inf);
        for (std::size_t idx = 0; idx < size; ++idx) {
            std::size_t j = 0;
            for (std::size_t r = 0; r < R; ++r)
                if (!open_chain[r]) j += mvec[idx][r] * stride[r];
            lHbar[j] = busyp_lse(std::vector<double>{lHbar[j], lH[idx]});
        }
    }

    std::vector<double> b(n.size(), 0.0);
    for (std::size_t t = 0; t < n.size(); ++t) {
        const std::size_t nt = n[t];
        std::vector<double> num, den;
        for (std::size_t idx = 0; idx < size; ++idx) {
            // the level set is |m| for the aggregate busy period and m_r for the
            // class-r one; only chain-r arrivals move m_r
            std::size_t level = 0;
            if (jobclass < 0)
                for (std::size_t r = 0; r < R; ++r) level += mvec[idx][r];
            else
                level = mvec[idx][static_cast<std::size_t>(jobclass)];
            if (level >= nt) {
                std::size_t j = 0;
                for (std::size_t r = 0; r < R; ++r)
                    if (!open_chain[r]) j += (pop[r] - mvec[idx][r]) * stride[r];
                num.push_back(lG[idx] + lHbar[j]);
            }
            if (level + 1 != nt) continue;
            std::vector<double> terms;
            for (std::size_t r = 0; r < R; ++r) {
                if (jobclass >= 0 && r != static_cast<std::size_t>(jobclass)) continue;
                if (A[r] <= 0) continue;
                if (!open_chain[r] && pop[r] == mvec[idx][r]) continue;
                std::size_t j = 0;
                for (std::size_t s = 0; s < R; ++s) {
                    if (open_chain[s]) continue;
                    // a closed chain conserves jobs, so the departing one is removed
                    j += (pop[s] - mvec[idx][s] - (s == r ? 1 : 0)) * stride[s];
                }
                terms.push_back(std::log(A[r]) + lHbar[j]);
            }
            if (!terms.empty()) den.push_back(lG[idx] + busyp_lse(terms));
        }
        b[t] = std::exp(busyp_lse(num) - busyp_lse(den));
    }
    return b;
}

}  // namespace detail

/**
 * Mean busy period of order n for the subnetwork, multichain.
 *
 * @param alpha  (J x R) relative arrival rates, one column per chain
 * @param mu     (J x R) service rates, the chain-r rate at node j
 * @param P      routing matrices, one per chain (size 1 = shared by all chains)
 * @param N      population per chain, infinite entries for an open chain
 * @param subnet zero-based node indexes forming the subnetwork
 * @param n      busy period orders, counting the jobs of every chain
 * @param gamma  (J x R) external arrival rates, empty for a closed network
 * @param phi    (J x K) dimensionless load-dependent scaling, empty = single server
 * @param tol    relative tolerance of the open-network tail truncation
 * @param jobclass zero-based chain whose own jobs are counted, -1 for every chain
 */
template <class T>
std::vector<double> pfqn_busyp_multiclass(const Matrix<T>& alpha, const Matrix<T>& mu,
                                          const std::vector<Matrix<T>>& P,
                                          const std::vector<double>& N,
                                          const std::vector<std::size_t>& subnet,
                                          const std::vector<std::size_t>& n,
                                          const Matrix<T>& gamma = Matrix<T>(),
                                          const Matrix<T>& phi = Matrix<T>(),
                                          double tol = PFQN_BUSYP_DEFAULT_TOL,
                                          int jobclass = -1) {
    const std::size_t J = alpha.rows(), R = alpha.cols();
    if (N.size() != R)
        throw InputError(
            "pfqn_busyp_multiclass: the population vector must have one entry per chain");
    bool is_closed = true, is_open = true;
    for (std::size_t r = 0; r < R; ++r) {
        if (std::isinf(N[r]))
            is_closed = false;
        else
            is_open = false;
    }
    const bool is_mixed = !is_closed && !is_open;

    std::vector<std::size_t> target = subnet;
    std::sort(target.begin(), target.end());
    target.erase(std::unique(target.begin(), target.end()), target.end());
    if (target.empty())
        throw InputError("pfqn_busyp_multiclass: the subnetwork must be non-empty");
    if (is_closed && target.size() >= J)
        throw InputError(
            "pfqn_busyp_multiclass: in a closed network the subnetwork must be a proper "
            "subset of the nodes");
    if (target.back() >= J)
        throw InputError("pfqn_busyp_multiclass: the subnetwork indexes are out of range");
    std::vector<bool> in_subnet(J, false);
    for (std::size_t i = 0; i < target.size(); ++i) in_subnet[target[i]] = true;
    std::vector<std::size_t> compl_nodes;
    for (std::size_t j = 0; j < J; ++j)
        if (!in_subnet[j]) compl_nodes.push_back(j);

    std::size_t totalN = 0;
    for (std::size_t r = 0; r < R; ++r)
        if (!std::isinf(N[r])) totalN += static_cast<std::size_t>(std::llround(N[r]));

    std::vector<std::vector<double>> scaling(J);
    for (std::size_t i = 0; i < J; ++i) {
        if (phi.rows() == J && phi.cols() > 0) {
            scaling[i].resize(phi.cols());
            for (std::size_t k = 0; k < phi.cols(); ++k)
                scaling[i][k] = num_traits<T>::to_double(phi(i, k));
        } else {
            scaling[i].assign(std::max<std::size_t>(1, totalN), 1.0);
        }
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
            "pfqn_busyp_multiclass: no job ever enters the subnetwork, its busy period is "
            "undefined");
    if (jobclass >= static_cast<int>(R))
        throw InputError("pfqn_busyp_multiclass: the job class index is out of range");
    if (jobclass >= 0 && A[static_cast<std::size_t>(jobclass)] <= 0)
        throw InputError(
            "pfqn_busyp_multiclass: no job of that class ever enters the subnetwork");

    if (is_open && !is_mixed) {
        // exact reduction to the single-chain routine on a per-station scalar; the
        // synthetic problem carries no routing, the whole inflow riding on gamma
        std::vector<double> rho(J, 0.0);
        for (std::size_t i = 0; i < J; ++i)
            for (std::size_t r = 0; r < R; ++r) rho[i] += L[i][r];
        Matrix<double> zeroP(J, J, 0.0);
        std::vector<double> gsyn(J, 0.0);
        if (jobclass < 0) {
            gsyn[target[0]] = inflow;
        } else {
            // the class-r marginal is geometric in rho_ir/(1-rho_i+rho_ir), NOT in
            // rho_ir: the other classes inflate the queue the class-r jobs sit in.
            // That collapse assumes a load-INDEPENDENT station.
            const std::size_t rr = static_cast<std::size_t>(jobclass);
            for (std::size_t t = 0; t < target.size(); ++t)
                for (std::size_t k = 0; k < scaling[target[t]].size(); ++k)
                    if (scaling[target[t]][k] != 1.0)
                        throw InputError(
                            "pfqn_busyp_multiclass: a per-class busy period of an open "
                            "subnetwork requires load-independent stations");
            for (std::size_t i = 0; i < J; ++i) {
                const double den = 1.0 - rho[i] + L[i][rr];
                rho[i] = (den > 0) ? L[i][rr] / den : 0.0;
            }
            gsyn[target[0]] = A[rr];
        }
        const std::function<double(std::size_t, std::size_t)> rates =
            [&](std::size_t j, std::size_t k) {
                return scaling[j][std::min(k, scaling[j].size()) - 1];
            };
        return pfqn_busyp(rho, rates, zeroP, std::numeric_limits<double>::infinity(), target,
                          n, gsyn, tol)
            .b;
    }

    if (!is_mixed) {
        const std::size_t bound = (jobclass < 0)
            ? totalN
            : static_cast<std::size_t>(std::llround(N[static_cast<std::size_t>(jobclass)]));
        for (std::size_t t = 0; t < n.size(); ++t)
            if (n[t] < 1 || n[t] > bound)
                throw InputError(
                    "pfqn_busyp_multiclass: the busy period order must be an integer in "
                    "1..sum(N), or in 1..N(r) for the busy period of class r alone");
        std::vector<std::size_t> pop(R, 0);
        for (std::size_t r = 0; r < R; ++r)
            pop[r] = static_cast<std::size_t>(std::llround(N[r]));
        const std::vector<bool> open_chain(R, false);
        return detail::busyp_lattice(L, scaling, target, compl_nodes, pop, n, A, jobclass,
                                     open_chain);
    }

    // A mixed model grows the truncation of its open dimensions until the answer stops
    // moving; that truncation is the only approximation in this branch.
    std::vector<bool> open_chain(R, false);
    for (std::size_t r = 0; r < R; ++r) open_chain[r] = std::isinf(N[r]);
    std::size_t nmax = 1;
    for (std::size_t t = 0; t < n.size(); ++t) nmax = std::max(nmax, n[t]);
    std::size_t trunc = 8 + 2 * nmax;
    std::vector<double> prev;
    while (true) {
        std::vector<std::size_t> pop(R, 0);
        for (std::size_t r = 0; r < R; ++r)
            pop[r] = open_chain[r] ? trunc : static_cast<std::size_t>(std::llround(N[r]));
        const std::vector<double> b = detail::busyp_lattice(L, scaling, target, compl_nodes,
                                                           pop, n, A, jobclass, open_chain);
        if (!prev.empty()) {
            bool settled = true;
            for (std::size_t t = 0; t < b.size(); ++t)
                if (std::fabs(b[t] - prev[t]) > 1e-10 * std::fabs(b[t])) settled = false;
            if (settled) return b;
        }
        prev = b;
        trunc *= 2;
        if (trunc > 4096)
            throw NumericError(
                "pfqn_busyp_multiclass: the mixed busy period did not converge, so some "
                "station of the subnetwork is nearly saturated");
    }
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BUSYP_MULTICLASS_H
