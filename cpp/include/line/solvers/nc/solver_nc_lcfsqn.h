/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_LCFSQN_H
#define LINE_SOLVERS_NC_SOLVER_NC_LCFSQN_H

/**
 * The two-station LCFS + LCFS-PR closed network. Port of `solver_nc_lcfsqn.m`.
 *
 * A last-come-first-served station WITHOUT preemption is not a BCMP station:
 * the queue-length distribution depends on the order of arrival, so no
 * product form and no ordinary normalizing constant exist for it. The pair
 * LCFS / LCFS-PR does admit a closed form (Casale, QUESTA 2026) whose constant
 * `pfqn_lcfsqn_ca` computes by recursion, and whose per-class measures are sums
 * of PERMANENTS over the boundary position between the two stations.
 *
 * THE PERMANENT FORMULAS ASSUME ONE JOB PER CLASS. A class with N_r > 1 is
 * expanded into N_r exchangeable single-job copies, the constant is rescaled to
 * the distinguishable-jobs one G_exp = G prod_r N_r!, and the per-copy measures
 * are scaled back by N_r. This is what makes the routine cost
 * (sum N)! -- exponential in the population -- so it is a special-case analyzer,
 * not a general one.
 *
 * Arithmetic: the permanents are exact-capable (`pfqn_perm` uses the integer
 * Pascal recurrence), but lG is a log, so the analyzer is guarded on
 * `has_transcendental` for that alone.
 */

#include <cmath>
#include <vector>

#include "line/api/pfqn/pfqn_lcfsqn_ca.h"
#include "line/api/pfqn/pfqn_perm.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

namespace detail {

/**
 * `make_Tx`: the (K-1 x K-1) throughput matrix at boundary position xt, with
 * the expanded copy r removed. Row i is the copy, column j the position.
 */
template <class T>
Matrix<T> lcfs_make_Tx(const std::vector<T>& alpha, const std::vector<T>& beta, std::size_t xt,
                       std::size_t K, std::size_t r) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Tx(K - 1, K - 1, zero);
    std::size_t idx = 0;
    for (std::size_t i = 0; i < K; ++i) {
        if (i + 1 == r) continue;
        for (std::size_t j = 1; j <= xt - 1; ++j) Tx(idx, j - 1) = num_pow_int(alpha[i], j);
        for (std::size_t j = 1; j <= K - xt; ++j)
            Tx(idx, xt - 2 + j) = T(num_pow_int(alpha[i], xt + j - 1) * beta[i]);
        ++idx;
    }
    return Tx;
}

/**
 * `make_Yx`: the (K x K) queue-length matrix at boundary position xt. The
 * LCFS-PR side of copy r is left at zero, which is what makes the permanent
 * count the states holding that copy at station 1.
 */
template <class T>
Matrix<T> lcfs_make_Yx(const std::vector<T>& alpha, const std::vector<T>& beta, std::size_t xt,
                       std::size_t K, std::size_t r) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> Y(K, K, zero);
    for (std::size_t i = 0; i < K; ++i) {
        for (std::size_t j = 1; j <= xt; ++j) Y(i, j - 1) = num_pow_int(alpha[i], j);
        for (std::size_t j = 1; j <= K - xt; ++j)
            if (i + 1 != r) Y(i, xt + j - 1) = T(num_pow_int(alpha[i], xt + j - 1) * beta[i]);
    }
    return Y;
}

}  // namespace detail

/**
 * Port of `solver_nc_lcfsqn.m`.
 *
 * @param sn         the refreshed struct
 * @param opt        solver controls; unused, the closed form has no tuning
 * @param lcfsStat   1-based index of the LCFS station
 * @param lcfsprStat 1-based index of the LCFS-PR station
 */
template <class T>
NcSolution<T> solver_nc_lcfsqn(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt,
                               std::size_t lcfsStat, std::size_t lcfsprStat) {
    (void)opt;
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)lcfsStat;
        (void)lcfsprStat;
        throw UnsupportedError(
            "solver_nc_lcfsqn: the LCFS closed form reports lG = log(G); this backend has no "
            "transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, R = sn.nclasses;

        std::vector<T> alpha(R, zero), beta(R, zero);
        std::vector<int> N(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            if (std::isinf(sn.classes[r].population))
                throw UnsupportedError("solver_nc_lcfsqn: requires a closed queueing network");
            N[r] = static_cast<int>(std::llround(sn.classes[r].population));
            if (N[r] <= 0) continue;
            const double mu_l = num_traits<T>::to_double(sn.rates(lcfsStat - 1, r));
            const double mu_p = num_traits<T>::to_double(sn.rates(lcfsprStat - 1, r));
            if (!(mu_l > 0.0) || !std::isfinite(mu_l))
                throw UnsupportedError(
                    "solver_nc_lcfsqn: invalid service rate at the LCFS station");
            if (!(mu_p > 0.0) || !std::isfinite(mu_p))
                throw UnsupportedError(
                    "solver_nc_lcfsqn: invalid service rate at the LCFS-PR station");
            alpha[r] = T(num_traits<T>::from_int(1) / sn.rates(lcfsStat - 1, r));
            beta[r] = T(num_traits<T>::from_int(1) / sn.rates(lcfsprStat - 1, r));
        }

        const T G = pfqn::pfqn_lcfsqn_ca(alpha, beta, N).G;
        out.sol.lG = G > zero ? num_traits<T>::log_as_double(G)
                              : -std::numeric_limits<double>::infinity();

        std::size_t K = 0;
        for (int v : N) K += static_cast<std::size_t>(v);

        // expand each class into N_r exchangeable single-job copies
        std::vector<T> alphaE, betaE;
        std::vector<std::size_t> ecls(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            ecls[r] = alphaE.size() + 1;
            for (int a = 0; a < N[r]; ++a) {
                alphaE.push_back(alpha[r]);
                betaE.push_back(beta[r]);
            }
        }
        T Gexp = G;
        for (std::size_t r = 0; r < R; ++r)
            Gexp = T(Gexp * num_factorial<T>(static_cast<unsigned>(N[r])));

        const std::vector<int> ones(K, 1);
        Matrix<T> Q_l(2, R, zero), U_l(2, R, zero), T_l(2, R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] <= 0) continue;
            const std::size_t e = ecls[r];
            T Tcopy = zero, Qcopy = zero;
            for (std::size_t xt = 1; xt <= K; ++xt) {
                const Matrix<T> Tx = detail::lcfs_make_Tx(alphaE, betaE, xt, K, e);
                const std::vector<int> onesTx(K - 1, 1);
                Tcopy = T(Tcopy + num_pow_int(alphaE[e - 1], xt - 1) *
                                      pfqn::pfqn_perm(Tx, onesTx) / Gexp);
                const Matrix<T> Yx = detail::lcfs_make_Yx(alphaE, betaE, xt, K, e);
                Qcopy = T(Qcopy + pfqn::pfqn_perm(Yx, ones) / Gexp);
            }
            const T nr = num_traits<T>::from_int(N[r]);
            T_l(0, r) = T(nr * Tcopy);
            T_l(1, r) = T(nr * Tcopy);
            Q_l(0, r) = T(nr * Qcopy);
            Q_l(1, r) = T(nr - Q_l(0, r));  // conservation at the LCFS-PR station
        }
        for (std::size_t r = 0; r < R; ++r) {
            U_l(0, r) = T(T_l(0, r) * alpha[r]);
            U_l(1, r) = T(T_l(1, r) * beta[r]);
        }

        Matrix<T> Q(M, R, zero), U(M, R, zero), Tp(M, R, zero), Rt(M, R, zero);
        std::vector<T> X(R, zero), C(R, zero);
        for (std::size_t r = 0; r < R; ++r) {
            Q(lcfsStat - 1, r) = Q_l(0, r);
            Q(lcfsprStat - 1, r) = Q_l(1, r);
            U(lcfsStat - 1, r) = U_l(0, r);
            U(lcfsprStat - 1, r) = U_l(1, r);
            if (N[r] <= 0) continue;
            X[r] = T_l(0, r);
            Tp(lcfsStat - 1, r) = T_l(0, r);
            Tp(lcfsprStat - 1, r) = T_l(1, r);
        }
        for (std::size_t i : {lcfsStat, lcfsprStat})
            for (std::size_t r = 0; r < R; ++r)
                if (Tp(i - 1, r) > zero) Rt(i - 1, r) = T(Q(i - 1, r) / Tp(i - 1, r));
        for (std::size_t r = 0; r < R; ++r)
            if (N[r] > 0) C[r] = T(Rt(lcfsStat - 1, r) + Rt(lcfsprStat - 1, r));

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = Rt;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C = C;
        out.sol.iter = 1;
        out.sol.method = "lcfsqn.ca";
        out.actualmethod = "lcfsqn.ca";
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_LCFSQN_H
