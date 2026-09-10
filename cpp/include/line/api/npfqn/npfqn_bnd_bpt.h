/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_NPFQN_BND_BPT_H
#define LINE_API_NPFQN_NPFQN_BND_BPT_H

/**
 * First-order linear-programming relaxation of the achievable region of a
 * multiclass open Markovian queueing network.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_bnd_bpt.m, cross-checked against
 * jar/src/main/java/jline/api/npfqn/Npfqn_bnd_bpt.java.
 *
 * Returns a LOWER bound on sum_r c_r x_r, where x_r is the mean sojourn time of
 * class r, valid for EVERY non-idling scheduling policy. A "class" here is a
 * buffer with its own exponential service rate and its own Markovian routing,
 * so a station serving several customer types owns one class per type. The
 * network is open: class r receives external Poisson arrivals at rate
 * lambda0(r) and, on completing service, becomes class r' with probability
 * P(r,r') or leaves with the row deficit.
 *
 * METHOD. Uniformize the chain and let R(t) = sum_r f(r) n_r(t) for an
 * arbitrary vector f. The steady-state balance of E[R^2] is an identity
 * quadratic in f; since it holds for every f, the two sides' coefficient
 * matrices agree entrywise. Diagonal entries give one equation per class,
 * off-diagonal entries one per unordered pair, in the variables
 *
 *   x_r    = E[T_r],  the mean sojourn time of class r,
 *   I(r,l) = E[1{server sigma(r) busy with class r} n_l],
 *   N(i,l) = E[1{server i idle} n_l].
 *
 * A third block states that the events "station i serves class r" and "station
 * i idle" are mutually exclusive and exhaustive, so their terms sum to
 * E[n_l] = lambda_l x_l. Minimizing over this polyhedron is a relaxation of the
 * achievable region, hence a lower bound.
 *
 * EXACT ON M/M/1. The LP reduces to mu*I11 - lambda^2*x = lambda and
 * I11 + N11 = lambda*x with N11 >= 0, whence x >= 1/(mu-lambda) with equality.
 *
 * NOT INCLUDED, DELIBERATELY. The valid inequality I(r,r) >= rho_r would
 * tighten the relaxation but is not part of the reference's characterization,
 * and reproducing the reference's published bounds is the acceptance test.
 *
 * ARITHMETIC. Rational-clean: the equations are polynomial in the data and the
 * dense simplex is exact, so at T = line::Rational the returned value is the
 * exact optimum of the exact polytope. Nothing here is transcendental.
 *
 * Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994). Optimization
 * of multiclass queueing networks: polyhedral and nonlinear characterizations
 * of achievable performance. Annals of Applied Probability 4(1), 43-75. See
 * also D. Bertsimas (1995), Queueing Systems 21, 337-389, Theorem 9, which
 * restates the same characterization.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lp_highs.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"

namespace line {
namespace npfqn {

template <class T>
struct BndBpt {
    T zlb = T();                 ///< lower bound on sum_r c_r x_r
    std::vector<T> x;            ///< the x block of the LP optimizer
    std::vector<T> lambda;       ///< effective arrival rate of each class
    std::vector<T> rho;          ///< per-class utilization lambda_r/mu_r
    std::vector<T> rhoStation;   ///< per-station utilization
    std::size_t nvars = 0;
    std::size_t nrows = 0;
};

/**
 * @param lambda0   external Poisson arrival rate into each class (0 if none)
 * @param mu        exponential service rate of each class
 * @param P         K x K routing, P(r,r') = P(class r becomes r' after service)
 * @param stationOf zero-based station index of each class
 * @param c         objective weights; empty means all ones
 */
template <class T>
BndBpt<T> npfqn_bnd_bpt(const std::vector<T>& lambda0, const std::vector<T>& mu,
                        const Matrix<T>& P, const std::vector<std::size_t>& stationOf,
                        const std::vector<T>& c) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const std::size_t K = lambda0.size();
    if (mu.size() != K || stationOf.size() != K)
        throw InputError("npfqn_bnd_bpt: lambda0, mu and stationOf must have the same length");
    if (P.rows() != K || P.cols() != K)
        throw InputError("npfqn_bnd_bpt: P must be K x K");
    std::vector<T> cost = c;
    if (cost.empty()) cost.assign(K, one);
    if (cost.size() != K) throw InputError("npfqn_bnd_bpt: c must have K entries");
    for (std::size_t r = 0; r < K; ++r) {
        if (!(mu[r] > zero))
            throw InputError("npfqn_bnd_bpt: every class needs a strictly positive service rate");
        T rowSum = zero;
        for (std::size_t s = 0; s < K; ++s) rowSum = T(rowSum + P(r, s));
        if (rowSum > T(one + num_traits<T>::from_double(1e-9)))
            throw InputError("npfqn_bnd_bpt: the routing matrix has a row summing above one");
    }
    std::size_t M = 0;
    for (std::size_t r = 0; r < K; ++r) M = std::max(M, stationOf[r] + 1);

    // ---- traffic equations, lambda = lambda0 + P' lambda ----
    Matrix<T> ImPt(K, K, zero);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) ImPt(i, j) = T((i == j ? one : zero) - P(j, i));
    Matrix<T> rhs0(K, 1, zero);
    for (std::size_t r = 0; r < K; ++r) rhs0(r, 0) = lambda0[r];
    const Matrix<T> lamM = matmul(inverse(ImPt), rhs0);

    BndBpt<T> out;
    out.lambda.assign(K, zero);
    out.rho.assign(K, zero);
    out.rhoStation.assign(M, zero);
    for (std::size_t r = 0; r < K; ++r) {
        out.lambda[r] = lamM(r, 0) < zero ? zero : lamM(r, 0);
        out.rho[r] = T(out.lambda[r] / mu[r]);
        out.rhoStation[stationOf[r]] = T(out.rhoStation[stationOf[r]] + out.rho[r]);
    }
    for (std::size_t i = 0; i < M; ++i)
        if (!(out.rhoStation[i] < one))
            throw UnsupportedError("npfqn_bnd_bpt: station " + std::to_string(i + 1) +
                                   " is saturated: no policy stabilizes the network");

    // ---- variable layout: x(r) | I(r,l) | N(i,l) ----
    const std::size_t oI = K, oN = K + K * K, nv = K + K * K + M * K;
    lp::LpModel<T> m(nv);
    m.set_maximize(false);

    // (a) diagonal equations, test function n_r^2
    for (std::size_t r = 0; r < K; ++r) {
        m.row_clear();
        m.row_add(oI + r * K + r, T(two * mu[r]));
        for (std::size_t w = 0; w < K; ++w)
            if (P(w, r) != zero) m.row_add(oI + w * K + r, T(-(two * mu[w] * P(w, r))));
        m.row_add(r, T(-(two * lambda0[r] * out.lambda[r])));
        m.emit(lp::LpSense::EQ, T(two * out.lambda[r] * (one - P(r, r))));
    }

    // (b) off-diagonal equations, test function n_r*n_s
    for (std::size_t r = 1; r < K; ++r) {
        for (std::size_t s = 0; s < r; ++s) {
            m.row_clear();
            m.row_add(oI + r * K + s, mu[r]);
            m.row_add(oI + s * K + r, mu[s]);
            for (std::size_t w = 0; w < K; ++w) {
                if (P(w, r) != zero) m.row_add(oI + w * K + s, T(-(mu[w] * P(w, r))));
                if (P(w, s) != zero) m.row_add(oI + w * K + r, T(-(mu[w] * P(w, s))));
            }
            m.row_add(s, T(-(lambda0[r] * out.lambda[s])));
            m.row_add(r, T(-(lambda0[s] * out.lambda[r])));
            m.emit(lp::LpSense::EQ,
                   T(-(out.lambda[r] * P(r, s)) - out.lambda[s] * P(s, r)));
        }
    }

    // (c) exhaustiveness at each station
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t l = 0; l < K; ++l) {
            m.row_clear();
            for (std::size_t r = 0; r < K; ++r)
                if (stationOf[r] == i) m.row_add(oI + r * K + l, one);
            m.row_add(oN + i * K + l, one);
            m.row_add(l, T(-out.lambda[l]));
            m.emit(lp::LpSense::EQ, zero);
        }
    }

    for (std::size_t r = 0; r < K; ++r) m.set_cost(r, cost[r]);

    out.nvars = m.num_vars();
    out.nrows = m.num_rows();
    const lp::LpSolution<T> s = lp::lp_solve(m);
    if (!s.ok())
        throw UnsupportedError(std::string("npfqn_bnd_bpt: the achievable-region LP did not "
                                           "solve to optimality (") +
                               lp::lp_status_name(s.status) + ")");
    out.zlb = s.objective;
    out.x.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r) out.x[r] = s.x[r];
    return out;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_NPFQN_BND_BPT_H
