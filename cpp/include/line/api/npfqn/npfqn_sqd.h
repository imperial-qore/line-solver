/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_SQD_H
#define LINE_API_NPFQN_SQD_H

/**
 * Smith Queue Decomposition (SQD): approximate MVA for closed networks under
 * Blocking-After-Service (manufacturing / transfer blocking).
 *
 * Templated port of matlab/src/api/npfqn/npfqn_sqd.m. Each finite-buffer
 * station is described by a load-dependent EFFECTIVE service rate calibrated
 * from an M/M/1/K blocking probability, and downstream blocking is propagated
 * through the effective routing between service stations. The recursion is an
 * AMVA-style sweep over the population 1..N:
 *
 *   mu_i(n)  = n V1_i exp( -[ max(0, (n-1)/beta_i) ]^{gamma_i} )        (13)
 *   W^buf_i  = (1 + n) / mu_i(n)                                        (18)
 *   W^svr_i  = ST_i (1 + L^svr_i) + (server blocking time)              (17)
 *   X        = pop / sum_i V_i (W^buf_i + W^svr_i)
 *   L^buf_i  = X V_i W^buf_i,   L^svr_i = X V_i W^svr_i
 *
 * with V1_i deflated by the blocking probability after every population step.
 * Delay (INF / EXT) stations are infinite-capacity pure delays. The method is
 * single-chain: it consumes the CHAIN-AGGREGATED demands, i.e. the first
 * column of MATLAB's sn_get_demands_chain output.
 *
 * SIGNATURE. MATLAB's entry point is `npfqn_sqd(sn, N, ...)` and its first act
 * is to unpack `sn` into six plain arrays; after that line the routine is pure
 * numerics with no reference to the model layer at all. The port therefore
 * takes those arrays directly, exactly as the whole pfqn_* family takes L, N, Z
 * rather than a Network. THE NAME IS UNCHANGED so that the registry keeps its
 * one-to-one mapping to the MATLAB source. To reconstruct the MATLAB call, a
 * caller unpacks:
 *
 *   ST                 = STchain(:,1)   from sn_get_demands_chain(sn)
 *   V                  = Vchain(:,1)    from sn_get_demands_chain(sn)
 *   isDelay(i)         = sn.sched(i) is SchedStrategy.INF or .EXT
 *   cap(i)             = infinity when isDelay(i), else sn.cap(i), with any
 *                        sn.cap(i) > 1e14 also treated as infinite
 *   rt                 = sn.rt          (stateful-indexed routing matrix)
 *   stationToStateful  = sn.stationToStateful  (ONE-based, as MATLAB stores it)
 *   nclasses           = sn.nclasses
 *   N                  = sn.nclosedjobs when the caller passes nothing
 *
 * `computeEffectiveRouting` is NOT hoisted out: collapsing pass-through delay
 * nodes into a station-to-station chain is algorithm content, not unpacking, so
 * it is ported here and reads `rt` in the same class-1 slice the reference does,
 * rt((sf_i-1) nclasses + 1, (sf_j-1) nclasses + 1).
 *
 * NO ARITY-COLLIDING OVERLOAD. MATLAB accepts one to seven positional
 * arguments; every overload here takes eight or nine, and of unrelated types,
 * so a call transcribed from MATLAB cannot bind to one of them and silently
 * solve a different model.
 *
 * Arithmetic: TRANSCENDENTAL, double and Real only. `exp` in (13), `log` in the
 * (beta, gamma) calibration, and real powers in both the calibration and the
 * M/M/1/K blocking probability rho^K / (1 - rho^{K+1}).
 *
 * UNEXPLAINED CONSTANTS, reproduced verbatim and NOT rationalized here.
 * `INITIAL_V1 = 692.192` is the default initial value of the effective-rate
 * scale V1_i, and `CALIBRATION_EPSILON = 0.05` is the target V_b of the fixed
 * heuristic calibration (mode 1). Neither is derived anywhere in the reference;
 * both arrive from the original SolverDBT contribution (Avinash Bommareddy,
 * Imperial College London FYP, 2026). 692.192 in particular has no stated
 * units or provenance, and the results DO depend on it whenever the population
 * sweep is short enough that the deflation V1 <- V1 (1 - pBlock) has not washed
 * the initial value out. Treat it as an inherited magic number pending an
 * answer from the contributor, not as a tuned parameter.
 *
 * REFERENCE DEFECTS: none found. Reproduces MATLAB to 1e-15 relative on every
 * combination of the four option switches tried.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace npfqn {

/** Return value of npfqn_sqd, mirroring [X, Q, U, R]. */
template <class T>
struct SqdResult {
    std::vector<T> X;  ///< (M) per-station throughput
    std::vector<T> Q;  ///< (M) per-station queue length
    std::vector<T> U;  ///< (M) per-station utilization, capped at one
    std::vector<T> R;  ///< (M) per-station residence time
};

/** The four option switches of the reference, with its defaults. */
template <class T>
struct SqdOptions {
    /// 0 = base (beta = K, gamma = 1), 1 = fixed heuristic, 2 = blocking-aware
    int calibrationMode = 0;
    /// include the server-blocking-time term in W^svr (MATLAB serverBlockingTime)
    bool serverBlockingTime = true;
    /// true = MATLAB neighborMode 'ownserver', false = 'downstream'
    bool ownServerNeighbor = false;
    /// true = MATLAB v1Policy 'fresh', false = 'compound'
    bool freshV1 = false;
    /// per-station initial V1; empty selects INITIAL_V1 for every station
    std::vector<T> initialV1;
};

namespace detail {

/** Infinite-capacity test; the routine is gated to inexact T, so this is safe. */
template <class T>
bool sqd_is_inf(const T& v) {
    return !std::isfinite(num_traits<T>::to_double(v));
}

/** Steady-state blocking probability of an M/M/1/K queue at load rho. */
template <class T>
T sqd_mm1k_blocking(const T& K, const T& rho) {
    using std::pow;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (rho <= num_traits<T>::from_double(1e-15)) return zero;
    if (num_abs(T(rho - one)) < num_traits<T>::from_double(1e-9))
        return T(one / T(K + one));
    return T(T(one - rho) * pow(rho, K) / T(one - pow(rho, T(K + one))));
}

/** Calibrate (beta, gamma) of the load-dependent effective service rate. */
template <class T>
void sqd_beta_gamma(const T& K, const T& pBlockDown, int mode, const T& calEps, T& beta,
                    T& gamma) {
    using std::log;
    using std::pow;
    const T one = num_traits<T>::from_int(1);
    const T a = num_traits<T>::from_int(2);
    if (K <= a || sqd_is_inf(K)) {
        beta = K;
        gamma = one;
        return;
    }
    const T b = K;
    T Va, Vb;
    if (mode == 2) {  // blocking-aware
        T pK = pBlockDown;
        const T lo = num_traits<T>::from_double(1e-6);
        const T hi = T(one - lo);
        if (pK > hi) pK = hi;
        if (pK < lo) pK = lo;
        Va = T(one - pK * T(a - one) / T(b - one));
        Vb = T(one - pK);
    } else if (mode == 1) {  // fixed heuristic
        Va = T(T(b - a) / b);
        Vb = calEps;
    } else {  // mode 0: base
        beta = K;
        gamma = one;
        return;
    }
    const T va_hi = num_traits<T>::from_double(0.999);
    const T va_lo = num_traits<T>::from_double(0.01);
    const T vb_lo = num_traits<T>::from_double(1e-6);
    if (Va > va_hi) Va = va_hi;
    if (Va < va_lo) Va = va_lo;
    if (Vb < vb_lo) Vb = vb_lo;
    if (Vb >= Va) {
        beta = K;
        gamma = one;
        return;
    }
    const T lnVa = log(Va);
    const T lnVb = log(Vb);
    gamma = T(log(T(lnVa / lnVb)) / log(T(T(a - one) / T(b - one))));
    const T g_lo = num_traits<T>::from_double(0.5);
    const T g_hi = num_traits<T>::from_int(10);
    if (gamma < g_lo) gamma = g_lo;
    if (gamma > g_hi) gamma = g_hi;
    beta = T(T(a - one) / pow(T(-lnVa), T(one / gamma)));
    if (!std::isfinite(num_traits<T>::to_double(gamma)) ||
        !std::isfinite(num_traits<T>::to_double(beta)) ||
        beta <= num_traits<T>::from_int(0)) {
        beta = K;
        gamma = one;
    }
}

/**
 * Station-to-station effective routing, collapsing pass-through delay nodes.
 * Reads the class-1 slice of the stateful-indexed routing matrix, exactly as
 * MATLAB's computeEffectiveRouting does.
 */
template <class T>
Matrix<T> sqd_effective_routing(const Matrix<T>& rt,
                                const std::vector<std::size_t>& stationToStateful,
                                std::size_t nclasses, const std::vector<bool>& isDelay) {
    const std::size_t M = isDelay.size();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> p(M, M, zero);
    const auto rtat = [&](std::size_t i, std::size_t j) {
        // stationToStateful is one-based, as MATLAB stores it
        const std::size_t ri = (stationToStateful[i] - 1) * nclasses;
        const std::size_t rj = (stationToStateful[j] - 1) * nclasses;
        if (ri >= rt.rows() || rj >= rt.cols())
            throw InputError("npfqn_sqd: stationToStateful indexes outside the routing matrix");
        return rt(ri, rj);
    };
    for (std::size_t i = 0; i < M; ++i) {
        if (isDelay[i]) continue;
        for (std::size_t j = 0; j < M; ++j) {
            const T p_ij = rtat(i, j);
            if (p_ij <= zero) continue;
            if (!isDelay[j]) {
                p(i, j) += p_ij;
            } else {
                for (std::size_t k = 0; k < M; ++k) {
                    if (isDelay[k]) continue;
                    const T p_jk = rtat(j, k);
                    if (p_jk > zero) p(i, k) += p_ij * p_jk;
                }
            }
        }
    }
    return p;
}

}  // namespace detail

/**
 * @param ST                (M) chain-aggregated service times, STchain(:,1)
 * @param V                 (M) chain-aggregated visit ratios, Vchain(:,1)
 * @param cap               (M) buffer capacities; use infinity for unbounded
 * @param isDelay           (M) true at the INF / EXT stations
 * @param rt                stateful-indexed routing matrix (sn.rt)
 * @param stationToStateful (M) ONE-based station-to-stateful map
 * @param nclasses          class count, the stride of rt
 * @param N                 total closed population
 * @param opt               the four option switches and the initial V1
 */
template <class T>
SqdResult<T> npfqn_sqd(const std::vector<T>& ST, const std::vector<T>& V,
                       const std::vector<T>& cap, const std::vector<bool>& isDelay,
                       const Matrix<T>& rt, const std::vector<std::size_t>& stationToStateful,
                       std::size_t nclasses, int N, const SqdOptions<T>& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_sqd requires transcendental arithmetic (the effective-rate "
                  "calibration uses exp, log and real powers)");
    using std::exp;
    using std::pow;
    const std::size_t M = ST.size();
    if (M == 0) throw InputError("npfqn_sqd: empty station list");
    if (V.size() != M || cap.size() != M || isDelay.size() != M ||
        stationToStateful.size() != M)
        throw InputError("npfqn_sqd: the per-station inputs disagree on the station count");
    if (nclasses == 0) throw InputError("npfqn_sqd: the class count must be positive");
    if (N < 0) throw InputError("npfqn_sqd: the population must be nonnegative");
    if (!opt.initialV1.empty() && opt.initialV1.size() != M)
        throw InputError("npfqn_sqd: initialV1 must have one entry per station");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T INITIAL_V1 = num_traits<T>::from_double(692.192);
    const T CALIBRATION_EPSILON = num_traits<T>::from_double(0.05);
    const T tiny = num_traits<T>::from_double(1e-10);
    const T eps15 = num_traits<T>::from_double(1e-15);

    const Matrix<T> pEff =
        detail::sqd_effective_routing(rt, stationToStateful, nclasses, isDelay);

    std::vector<T> V1(M, zero), v1init(M, zero), L_buf(M, zero), L_svr(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        v1init[i] = opt.initialV1.empty() ? INITIAL_V1 : opt.initialV1[i];
        V1[i] = v1init[i];
    }
    T X = zero;
    std::vector<T> W_buf(M, zero), W_svr(M, zero);

    for (int pop = 1; pop <= N; ++pop) {
        // ---- wait times ------------------------------------------------------
        for (std::size_t i = 0; i < M; ++i) {
            if (isDelay[i]) {
                W_buf[i] = zero;
                W_svr[i] = ST[i];
                continue;
            }
            T pBlockDown = zero;
            if (opt.ownServerNeighbor) {
                if (!detail::sqd_is_inf(cap[i]))
                    pBlockDown = detail::sqd_mm1k_blocking(cap[i], T(X * V[i] * ST[i]));
            } else {
                for (std::size_t j = 0; j < M; ++j)
                    if (!isDelay[j] && pEff(i, j) > zero && !detail::sqd_is_inf(cap[j]))
                        pBlockDown += pEff(i, j) *
                                      detail::sqd_mm1k_blocking(cap[j], T(X * V[j] * ST[j]));
            }

            T n = zero;
            if (opt.ownServerNeighbor) {
                n = L_svr[i];
            } else {
                for (std::size_t j = 0; j < M; ++j) n += pEff(i, j) * L_svr[j];
            }

            T mu_n;
            if (n < tiny) {
                mu_n = V1[i];
            } else {
                T beta, gamma;
                detail::sqd_beta_gamma(cap[i], pBlockDown, opt.calibrationMode,
                                       CALIBRATION_EPSILON, beta, gamma);
                T base = T(T(n - one) / beta);
                if (base < zero) base = zero;
                const T expArg = pow(base, gamma);
                mu_n = T(n * V1[i] * exp(T(-expArg)));  // Eq. 13
            }
            if (mu_n < tiny) mu_n = tiny;

            W_buf[i] = T(one / mu_n * T(one + n));   // Eq. 18
            W_svr[i] = T(ST[i] * T(one + L_svr[i]));  // Eq. 17

            if (opt.serverBlockingTime) {
                T bt = zero;
                for (std::size_t j = 0; j < M; ++j) {
                    if (isDelay[j] || !(pEff(i, j) > zero) || detail::sqd_is_inf(cap[j])) continue;
                    const T pBj = detail::sqd_mm1k_blocking(cap[j], T(X * V[j] * ST[j]));
                    const T denom = T(ST[i] + ST[j]);
                    T theta = zero;
                    if (denom > eps15) theta = T(ST[j] / denom);
                    bt += pEff(i, j) * pBj * ST[j] * theta;
                }
                W_svr[i] += bt;
            }
        }

        // ---- throughput ------------------------------------------------------
        T sumVW = zero;
        for (std::size_t i = 0; i < M; ++i) sumVW += V[i] * T(W_buf[i] + W_svr[i]);
        X = (sumVW > eps15) ? T(num_traits<T>::from_int(pop) / sumVW) : zero;

        // ---- queue lengths ---------------------------------------------------
        for (std::size_t i = 0; i < M; ++i) {
            L_buf[i] = X * V[i] * W_buf[i];
            L_svr[i] = X * V[i] * W_svr[i];
        }

        // ---- deflate V1 ------------------------------------------------------
        if (pop < N) {
            for (std::size_t i = 0; i < M; ++i) {
                if (isDelay[i]) continue;
                T pBlock = zero;
                if (opt.ownServerNeighbor) {
                    if (!detail::sqd_is_inf(cap[i]))
                        pBlock = detail::sqd_mm1k_blocking(cap[i], T(X * V[i] * ST[i]));
                } else {
                    for (std::size_t j = 0; j < M; ++j)
                        if (!isDelay[j] && pEff(i, j) > zero && !detail::sqd_is_inf(cap[j]))
                            pBlock += pEff(i, j) *
                                      detail::sqd_mm1k_blocking(cap[j], T(X * V[j] * ST[j]));
                }
                V1[i] = opt.freshV1 ? T(v1init[i] * T(one - pBlock))
                                    : T(V1[i] * T(one - pBlock));
                if (V1[i] < tiny) V1[i] = tiny;
            }
        }
    }

    SqdResult<T> res;
    res.X.assign(M, zero);
    res.Q.assign(M, zero);
    res.U.assign(M, zero);
    res.R.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const T T_i = X * V[i];
        const T Q_i = T(L_buf[i] + L_svr[i]);
        const T W_tot = T(W_buf[i] + W_svr[i]);
        res.X[i] = T_i;
        res.Q[i] = Q_i;
        res.R[i] = (T_i > eps15) ? T(Q_i / T_i) : W_tot;
        const T u = T(T_i * ST[i]);
        res.U[i] = (u < one) ? u : one;
    }
    return res;
}

/**
 * Overload with the reference's default option switches. Eight arguments, so it
 * cannot collide with a MATLAB call transcribed positionally (the reference
 * accepts at most seven, starting with an sn).
 */
template <class T>
SqdResult<T> npfqn_sqd(const std::vector<T>& ST, const std::vector<T>& V,
                       const std::vector<T>& cap, const std::vector<bool>& isDelay,
                       const Matrix<T>& rt, const std::vector<std::size_t>& stationToStateful,
                       std::size_t nclasses, int N) {
    return npfqn_sqd(ST, V, cap, isDelay, rt, stationToStateful, nclasses, N, SqdOptions<T>());
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_SQD_H
