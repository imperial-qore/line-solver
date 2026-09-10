/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DPFQN_DPFQN_NC_H
#define LINE_API_DPFQN_DPFQN_NC_H

/**
 * Normalizing constants of a discrete-time closed cycle of Bernoulli servers.
 *
 * Templated port of matlab/src/api/dpfqn/dpfqn_nc.m and dpfqn_ncld.m. With
 * q_j = 1 - p_j the queue length vector has the product form of Daduna (2001),
 * corollary 3.4,
 *
 *   pi(n_1,...,n_J) = prod_j (q_j/p_j)^n_j (1/q_j)^{1{n_j>0}} / G(N,J),
 *
 * whose extra factor on the busy nodes is what separates it from the
 * continuous-time Gordon-Newell form: a homogeneous cycle is uniform on the
 * state space in continuous time and is not here.
 *
 * `dpfqn_nc` runs the three-term recursion of proposition 3.18,
 *
 *   G(k,j) = G(k,j-1) + (q_j/p_j) G(k-1,j) + G(k-1,j-1),
 *
 * together with the arrival constants of proposition 3.19,
 *
 *   G1(k,J) = G(k-1,J-1) + (q_J/p_J) G1(k-1,J),   k >= 3.
 *
 * Unlike the continuous-time convolution algorithm neither recursion is
 * invariant to the numbering of the nodes. By lemma 7.3 the arrival constant is
 * the same at every node, so one family suffices and
 *
 *   throughput per slot   X   = G1(N,J) / G(N,J)   (equal at every node),
 *   utilization           U_j = X / p_j,
 *   tail probability      P(X_j >= k)
 *                         = (q_j/p_j)^k (1/q_j) G1(N-k+1,J) / G(N,J).
 *
 * The last identity is corollary 3.20(a) with its index corrected: as printed
 * there the right-hand side evaluates to P(X_j >= k+1). Corollary 3.20(c),
 * which transfers the tail from node 1 to node j, holds for k >= 1 only; at
 * k = 0 both tails are 1 while the stated ratio is q_1/q_j.
 *
 * `dpfqn_ncld` takes the state dependent case of theorem 3.2 by truncated
 * convolution of the per-node weights
 *
 *   w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h),
 *
 * with the complement constants (the cycle without node j) from a prefix/suffix
 * pass, so that P(X_j = n) = W[j][n] * Gc[j][N-n] / G[N]. Deconvolution is
 * never used, so a node with a near-unit service probability does not spoil the
 * accuracy of the other marginals.
 *
 * Everything here is a rational function of the service probabilities, so the
 * exact instantiation returns the constants with no rounding. That matters in
 * the interesting regime: a slow node makes (q/p)^N large, and the double
 * evaluation of the ratio G1(N-k+1)/G loses digits exactly there. No
 * transcendental appears in any metric: every quantity the analyzer reads is a
 * ratio of constants in the same common scale, and `lG` is a double-valued
 * diagnostic obtained through `to_double`, as in `solver_nc.h`.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace dpfqn {

/** Constants of the state independent cycle, in one common scale. */
template <class T>
struct DtNcResult {
    double lG;          ///< true log of G(N,J), in double whatever T is
    T G;                ///< G(N,J), in the common scale
    std::vector<T> G1;  ///< G1[k] = G_1(k,J), k = 0..N, same scale

    /** Throughput per slot, equal at every node of the cycle. */
    T throughput() const { return G1[G1.size() - 1] / G; }
};

/** Constants of the state dependent cycle, in one common scale. */
template <class T>
struct DtNcLdResult {
    double lG;                       ///< true log of G(N,J), in double whatever T is
    std::vector<T> G;                ///< G(k), k = 0..N
    std::vector<std::vector<T> > W;   ///< time-stationary node weights w_j(n)
    std::vector<std::vector<T> > Gc;  ///< cycle without node j, at population k
    std::vector<std::vector<T> > Wa;  ///< arrival weights v_j(n) = w_j(n) q_j(n)

    /** Marginal queue length law of node j, P(X_j = n) for n = 0..N. */
    std::vector<T> marginal(std::size_t j) const {
        const std::size_t N = G.size() - 1;
        std::vector<T> marg(N + 1);
        for (std::size_t n = 0; n <= N; ++n) {
            marg[n] = W[j][n] * Gc[j][N - n] / G[N];
        }
        return marg;
    }
};

/**
 * Propositions 3.18 and 3.19 for a cycle with state independent service
 * probabilities.
 *
 * @param p per-slot service completion probabilities p_j in (0,1)
 * @param N number of customers cycling in the J nodes
 */
template <class T>
DtNcResult<T> dpfqn_nc(const std::vector<T>& p, std::size_t N) {
    if (p.empty()) {
        throw InputError("dpfqn_nc: the cycle must contain at least one node");
    }
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t j = 0; j < p.size(); ++j) {
        if (!(p[j] > zero) || !(p[j] < one)) {
            throw InputError("dpfqn_nc: service probabilities must be in the open interval (0,1)");
        }
    }
    const std::size_t J = p.size();
    std::vector<T> x(J);
    for (std::size_t j = 0; j < J; ++j) {
        x[j] = (one - p[j]) / p[j];
    }

    // The recursion is linear and homogeneous in the whole table, so rescaling
    // every entry at once preserves it; that is what keeps (q/p)^N from
    // overflowing for a slow node and a large population. G1 is rescaled in
    // step so the two families stay in one common scale.
    std::vector<std::vector<T> > Gt(N + 1, std::vector<T>(J + 1, zero));
    for (std::size_t j = 0; j <= J; ++j) {
        Gt[0][j] = one;
    }
    std::vector<T> G1(N + 1, zero);
    double lscale = 0.0;
    const T limit = num_traits<T>::from_double(1e250);
    for (std::size_t k = 1; k <= N; ++k) {
        for (std::size_t j = 1; j <= J; ++j) {
            Gt[k][j] = Gt[k][j - 1] + x[j - 1] * Gt[k - 1][j] + Gt[k - 1][j - 1];
        }
        if (k == 1) {
            G1[1] = Gt[0][0];
        } else if (k == 2) {
            T s = x[0];
            for (std::size_t j = 1; j < J; ++j) {
                s = s + one / p[j];
            }
            G1[2] = s * Gt[0][0];
        } else {
            G1[k] = Gt[k - 1][J - 1] + x[J - 1] * G1[k - 1];
        }
        T mx = zero;
        for (std::size_t j = 0; j <= J; ++j) {
            if (Gt[k][j] > mx) {
                mx = Gt[k][j];
            }
        }
        if (mx > limit) {
            for (std::size_t a = 0; a <= N; ++a) {
                for (std::size_t b = 0; b <= J; ++b) {
                    Gt[a][b] = Gt[a][b] / mx;
                }
                G1[a] = G1[a] / mx;
            }
            lscale += std::log(num_traits<T>::to_double(mx));
        }
    }
    DtNcResult<T> out;
    out.G = Gt[N][J];
    out.lG = std::log(num_traits<T>::to_double(out.G)) + lscale;
    out.G1 = G1;
    return out;
}

namespace detail {

/** Convolution of two population tables truncated at N. */
template <class T>
std::vector<T> dt_conv(const std::vector<T>& a, const std::vector<T>& b, std::size_t N) {
    std::vector<T> c(N + 1, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k <= N; ++k) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t m = 0; m <= k; ++m) {
            s = s + a[m] * b[k - m];
        }
        c[k] = s;
    }
    return c;
}

}  // namespace detail

/**
 * Theorem 3.2 for a cycle whose service probabilities depend on the local
 * queue length.
 *
 * @param P service probabilities, P[j][n-1] = p_j(n) in (0,1]
 * @param N number of customers cycling in the J nodes
 */
template <class T>
DtNcLdResult<T> dpfqn_ncld(const std::vector<std::vector<T> >& P, std::size_t N) {
    if (P.empty()) {
        throw InputError("dpfqn_ncld: P must be a non-empty matrix of service probabilities");
    }
    const std::size_t J = P.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    DtNcLdResult<T> out;
    if (N == 0) {
        out.lG = 0.0;
        out.G.assign(1, one);
        out.W.assign(J, std::vector<T>(1, one));
        out.Gc = out.W;
        out.Wa = out.W;
        return out;
    }
    for (std::size_t j = 0; j < J; ++j) {
        if (P[j].size() < N) {
            throw InputError("dpfqn_ncld: P must supply p_j(n) for every n = 1..N");
        }
        for (std::size_t n = 0; n < N; ++n) {
            if (!(P[j][n] > zero) || P[j][n] > one) {
                throw InputError("dpfqn_ncld: service probabilities must be in the interval (0,1]");
            }
            // p_j(n)=1 is admissible only at the last reachable population,
            // where the missing q_j(n) never multiplies any weight.
            if (n + 1 < N && !(P[j][n] < one)) {
                throw InputError("dpfqn_ncld: service probabilities below the population bound must be "
                           "strictly less than 1");
            }
        }
    }

    // Per-node weights of theorem 3.2 by their exact recurrence,
    // w_j(0) = 1 and w_j(n) = w_j(n-1) q_j(n-1) / p_j(n) with q_j(0) := 1, then
    // a per-node division by the largest entry. The recurrence keeps every
    // weight in the field, so nothing here is transcendental; the divisions
    // cancel in every ratio below and are tracked only for lG.
    std::vector<double> shift(J, 0.0);
    out.W.assign(J, std::vector<T>(N + 1, zero));
    out.Wa.assign(J, std::vector<T>(N + 1, zero));
    for (std::size_t j = 0; j < J; ++j) {
        out.W[j][0] = one;
        for (std::size_t n = 1; n <= N; ++n) {
            const T q = (n >= 2) ? (one - P[j][n - 2]) : one;
            out.W[j][n] = out.W[j][n - 1] * q / P[j][n - 1];
        }
        T mx = out.W[j][0];
        for (std::size_t n = 1; n <= N; ++n) {
            if (out.W[j][n] > mx) {
                mx = out.W[j][n];
            }
        }
        shift[j] = std::log(num_traits<T>::to_double(mx));
        for (std::size_t n = 0; n <= N; ++n) {
            out.W[j][n] = out.W[j][n] / mx;
        }
        out.Wa[j][0] = out.W[j][0];
        for (std::size_t n = 1; n <= N; ++n) {
            out.Wa[j][n] = out.W[j][n] * (one - P[j][n - 1]);
        }
    }

    // Prefix/suffix convolution: pre[j] covers nodes 0..j-1 and suf[j] covers
    // nodes j..J-1, so the complement of node j is their convolution.
    std::vector<std::vector<T> > pre(J + 1, std::vector<T>(N + 1, zero));
    pre[0][0] = one;
    for (std::size_t j = 0; j < J; ++j) {
        pre[j + 1] = detail::dt_conv(pre[j], out.W[j], N);
    }
    std::vector<std::vector<T> > suf(J + 1, std::vector<T>(N + 1, zero));
    suf[J][0] = one;
    for (std::size_t j = J; j-- > 0;) {
        suf[j] = detail::dt_conv(suf[j + 1], out.W[j], N);
    }
    out.G = pre[J];
    out.Gc.assign(J, std::vector<T>(N + 1, zero));
    for (std::size_t j = 0; j < J; ++j) {
        out.Gc[j] = detail::dt_conv(pre[j], suf[j + 1], N);
    }

    double total = 0.0;
    for (std::size_t j = 0; j < J; ++j) {
        total += shift[j];
    }
    out.lG = std::log(num_traits<T>::to_double(out.G[N])) + total;
    return out;
}

}  // namespace dpfqn
}  // namespace line

#endif  // LINE_API_DPFQN_DPFQN_NC_H
