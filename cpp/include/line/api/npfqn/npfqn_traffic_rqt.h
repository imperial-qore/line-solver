/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_NPFQN_TRAFFIC_RQT_H
#define LINE_API_NPFQN_NPFQN_TRAFFIC_RQT_H

/**
 * Effective arrival processes of a network under the Robust Queueing calculus.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_traffic_rqt.m, cross-checked
 * against jar/src/main/java/jline/api/npfqn/Npfqn_traffic_rqt.java.
 *
 * The network characterization composes three operators: passage through a
 * queue with adversarial servers leaves the uncertainty set unchanged (robust
 * Burke, Theorem 4), superposition merges sets by Theorem 5, and thinning by a
 * fraction f scales the rate by f and the variability by f^(-1/alpha)
 * (Theorem 6). The resulting equations are
 *
 *   lambda_j = lambda0_j + sum_i lambda_i f_ij,
 *   Gamma_j  = (1/lambda_j) [ 1{a0_j=ab_j} (lambda0_j Gamma0_j)^(p_j)
 *                             + sum_i 1{ab_i=ab_j} (lambda_i Gamma_i)^(p_i) f_ij ]^(1/p_j),
 *
 * with p_j = ab_j/(ab_j-1) and ab_j the minimum tail coefficient among the
 * streams feeding j: the heaviest tail upstream dominates. Both are solved
 * exactly rather than iteratively. The rate equations are the usual traffic
 * equations, and in the variables z_j = (lambda_j Gamma_j)^(p_j) the variability
 * equations are linear as well, so each is one linear system; ab is obtained by
 * propagating the minimum to a fixed point.
 *
 * ARITHMETIC. Real exponents make this transcendental.
 *
 * Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
 * Operations Research 63(3), 676-700, Theorems 4-7 and 10.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

template <class T>
struct TrafficRqt {
    std::vector<T> lambda;  ///< effective arrival rate at each node
    std::vector<T> Gamma;   ///< effective variability parameter at each node
    std::vector<T> alpha;   ///< effective tail coefficient at each node
};

/**
 * @param lambda0 external arrival rate at each node, 0 where there is none
 * @param Gamma0  variability parameter of each external arrival process, which
 *                for a renewal stream is the interarrival standard deviation
 * @param alpha0  tail coefficient in (1,2] of each external arrival process
 * @param F       routing probabilities, F(i,j) = fraction of the jobs leaving
 *                node i that go to node j (row sums <= 1)
 */
template <class T>
TrafficRqt<T> npfqn_traffic_rqt(const std::vector<T>& lambda0, const std::vector<T>& Gamma0,
                                const std::vector<T>& alpha0, const Matrix<T>& F) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_traffic_rqt requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const std::size_t J = lambda0.size();
    TrafficRqt<T> out;
    out.lambda.assign(J, zero);
    out.Gamma.assign(J, zero);
    out.alpha.assign(J, two);

    // Traffic equations, lambda = (I - F')^{-1} lambda0.
    Matrix<T> ImFt(J, J, zero);
    for (std::size_t i = 0; i < J; ++i)
        for (std::size_t j = 0; j < J; ++j)
            ImFt(i, j) = T((i == j ? one : zero) - F(j, i));
    const Matrix<T> Xi = inverse(ImFt);
    for (std::size_t j = 0; j < J; ++j) {
        T s = zero;
        for (std::size_t i = 0; i < J; ++i) s = T(s + Xi(j, i) * lambda0[i]);
        out.lambda[j] = s;
    }

    // Effective tail coefficient: the minimum propagated along the routing graph.
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> alpha(J);
    for (std::size_t j = 0; j < J; ++j)
        alpha[j] = lambda0[j] > zero ? num_traits<T>::to_double(alpha0[j]) : inf;
    for (std::size_t it = 0; it < J; ++it) {
        bool changed = false;
        for (std::size_t j = 0; j < J; ++j)
            for (std::size_t i = 0; i < J; ++i)
                if (F(i, j) > zero && out.lambda[i] > zero && alpha[i] < alpha[j]) {
                    alpha[j] = alpha[i];
                    changed = true;
                }
        if (!changed) break;
    }
    for (std::size_t j = 0; j < J; ++j) {
        // an unreachable node keeps the light-tailed default
        if (!(alpha[j] < inf)) alpha[j] = 2.0;
        out.alpha[j] = num_traits<T>::from_double(alpha[j]);
    }

    // Variability equations, linear in z_j = (lambda_j Gamma_j)^(p_j).
    std::vector<T> p(J, two);
    for (std::size_t j = 0; j < J; ++j) p[j] = out.alpha[j] / (out.alpha[j] - one);
    Matrix<T> z0(J, 1, zero);
    for (std::size_t j = 0; j < J; ++j) {
        const double d = num_traits<T>::to_double(alpha0[j]) - alpha[j];
        if (lambda0[j] > zero && d < 1e-12 && d > -1e-12)
            z0(j, 0) = qsys::detail::num_pow(T(lambda0[j] * Gamma0[j]), p[j]);
    }
    Matrix<T> ImAt(J, J, zero);
    for (std::size_t i = 0; i < J; ++i)
        for (std::size_t j = 0; j < J; ++j) {
            T a = zero;
            const double d = alpha[j] - alpha[i];
            if (F(j, i) > zero && d < 1e-12 && d > -1e-12) a = F(j, i);
            ImAt(i, j) = T((i == j ? one : zero) - a);
        }
    const Matrix<T> z = matmul(inverse(ImAt), z0);

    for (std::size_t j = 0; j < J; ++j) {
        T zj = z(j, 0);
        if (zj < zero) zj = zero;
        if (out.lambda[j] > zero)
            out.Gamma[j] = qsys::detail::num_pow(zj, T(one / p[j])) / out.lambda[j];
    }
    return out;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_NPFQN_TRAFFIC_RQT_H
