/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CFTP_H
#define LINE_API_PFQN_CFTP_H

/**
 * Perfect stationary state sampling for closed single-class multiserver
 * product-form networks, by monotone Coupling From The Past.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_cftp.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/Pfqn_cftp.java. Reference: S. Kijima and
 * T. Matsui, "Approximate/Perfect Samplers for Closed Jackson Networks",
 * Winter Simulation Conference 2005.
 *
 * The chain moves one adjacent station pair at a time. A single uniform u
 * encodes both the pair, through lam = 1 + u (M-1) and j = floor(lam), and the
 * split of the pair's combined occupancy k = x_j + x_{j+1}, through the
 * fractional part Lambda used as an inverse-CDF argument against
 *
 *   w(s) proportional to alpha_j(s) alpha_{j+1}(k-s),
 *   log alpha_i(m) = m log L_i - sum_{t=1}^{m} log min(t, S_i).
 *
 * That update is monotone with respect to the componentwise partial order on
 * the population simplex, so running the top state (K,0,...,0) and the bottom
 * state (0,...,0,K) from -T with a FIXED randomness tape and doubling T until
 * they coalesce returns a draw from the exact stationary distribution
 * (Propp-Wilson). The randomness for the steps already simulated must be
 * REUSED as T doubles, which is why the tape grows at its far end and is
 * replayed oldest first; re-drawing it would destroy the perfection guarantee.
 *
 * The 'approx' method is the rapidly-mixing sampler M_A of the same paper: a
 * fixed number ceil(M(M-1)/2 log(K/eps)) of updates on uniformly random
 * DISTINCT (not necessarily adjacent) pairs, from an arbitrary feasible start.
 * It is not exact; the reference offers it because its running time is
 * deterministic.
 *
 * Log-space weights are kept in double, as in the reference: the split CDF is
 * compared against a double uniform, so carrying the weights at a higher
 * precision cannot change the sampled state. The reported mean queue length is
 * accumulated in the working arithmetic.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The output is a random state; the
 * balance functions are formed in the log domain.
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Sentinel for an infinite-server (delay) station, the reference's S = Inf. */
constexpr int cftp_inf_servers = -1;

/** Which sampler to run. */
enum class CftpMethod {
    Cftp,   ///< exact, monotone coupling from the past
    Approx  ///< the rapidly-mixing approximate sampler M_A
};

/** Return value of pfqn_cftp, mirroring [Q, X, T]. */
template <class T>
struct CftpResult {
    std::vector<T> Q;              ///< (M) empirical mean queue length
    Matrix<int> X;                 ///< (nsamples x M) sampled states, rows sum to K
    std::vector<long> horizon;     ///< (nsamples) coalescence horizon or step count
};

namespace detail {

/**
 * Inverse-CDF split of k jobs between stations i and j: the smallest s with
 * Lambda <= cdf(s). Weights are normalized by their maximum before
 * exponentiating, exactly as the reference does.
 */
inline int cftp_split_index(const std::vector<double>& logL, const Matrix<double>& logfac,
                            std::size_t i, std::size_t j, int k, double Lambda) {
    std::vector<double> lw(static_cast<std::size_t>(k) + 1);
    double m = -std::numeric_limits<double>::infinity();
    for (int s = 0; s <= k; ++s) {
        const double val = s * logL[i] - logfac(i, static_cast<std::size_t>(s)) +
                           (k - s) * logL[j] - logfac(j, static_cast<std::size_t>(k - s));
        lw[static_cast<std::size_t>(s)] = val;
        if (val > m) m = val;
    }
    double tot = 0.0;
    for (int s = 0; s <= k; ++s) {
        lw[static_cast<std::size_t>(s)] = std::exp(lw[static_cast<std::size_t>(s)] - m);
        tot += lw[static_cast<std::size_t>(s)];
    }
    double c = 0.0;
    for (int s = 0; s <= k; ++s) {
        c += lw[static_cast<std::size_t>(s)] / tot;
        if (Lambda <= c) return s;
    }
    return k;
}

/** One monotone update of an adjacent pair, driven by a single uniform. */
inline void cftp_monotone_update(std::vector<int>& x, double u, const std::vector<double>& logL,
                                 const Matrix<double>& logfac, std::size_t M) {
    const double lam = 1.0 + u * static_cast<double>(M - 1);
    std::size_t j = static_cast<std::size_t>(std::floor(lam));
    if (j > M - 1) j = M - 1;
    if (j < 1) j = 1;
    const double Lambda = lam - static_cast<double>(j);
    const std::size_t a = j - 1, b = j;  // 0-based pair (j, j+1)
    const int k = x[a] + x[b];
    const int l = cftp_split_index(logL, logfac, a, b, k, Lambda);
    x[a] = l;
    x[b] = k - l;
}

}  // namespace detail

/**
 * @param L        (M) demands L_i = theta_i / mu_i, strictly positive
 * @param N        total closed population K
 * @param S        (M) servers per station; cftp_inf_servers for a delay.
 *                 Empty for all single-server.
 * @param nsamples number of independent draws
 * @param method   exact CFTP or the approximate sampler
 * @param rng      explicit generator, advanced by the call
 */
template <class T>
CftpResult<T> pfqn_cftp(const std::vector<T>& L, int N, const std::vector<int>& S,
                        std::size_t nsamples, CftpMethod method, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_cftp requires transcendental arithmetic: it draws random states and forms "
                  "the station balance functions in the log domain");

    const std::size_t M = L.size();
    if (M < 2) throw InputError("pfqn_cftp: at least two stations are required");
    if (N < 0) throw InputError("pfqn_cftp: negative population");
    if (nsamples == 0) throw InputError("pfqn_cftp: at least one sample is required");
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < M; ++i)
        if (!(L[i] > zero)) throw InputError("pfqn_cftp: all demands L must be strictly positive");
    std::vector<int> Sv(M, 1);
    if (!S.empty()) {
        if (S.size() != M) throw InputError("pfqn_cftp: S has the wrong station count");
        Sv = S;
        for (std::size_t i = 0; i < M; ++i)
            if (Sv[i] == 0 || Sv[i] < cftp_inf_servers)
                throw InputError("pfqn_cftp: server counts must be positive or cftp_inf_servers");
    }
    const int K = N;

    // logfac(i, m) = sum_{t=1}^{m} log min(t, S_i), m = 0 .. K
    Matrix<double> logfac(M, static_cast<std::size_t>(K) + 1, 0.0);
    std::vector<double> logL(M);
    for (std::size_t i = 0; i < M; ++i) {
        logL[i] = num_traits<T>::log_as_double(L[i]);
        double acc = 0.0;
        for (int m = 1; m <= K; ++m) {
            const double cap = Sv[i] == cftp_inf_servers
                                   ? static_cast<double>(m)
                                   : static_cast<double>(m < Sv[i] ? m : Sv[i]);
            acc += std::log(cap);
            logfac(i, static_cast<std::size_t>(m)) = acc;
        }
    }

    CftpResult<T> res;
    res.X = Matrix<int>(nsamples, M, 0);
    res.horizon.assign(nsamples, 0);
    res.Q.assign(M, zero);

    std::vector<double> tape;
    std::vector<int> xU(M), xL(M), x(M);
    for (std::size_t smp = 0; smp < nsamples; ++smp) {
        if (method == CftpMethod::Cftp) {
            tape.clear();
            long Tback = 1;
            while (true) {
                // Extend the tape at its far end; the recent steps keep the
                // randomness they were already simulated with.
                while (static_cast<long>(tape.size()) < Tback) tape.push_back(mc_uniform01(rng));
                xU.assign(M, 0);
                xU[0] = K;
                xL.assign(M, 0);
                xL[M - 1] = K;
                for (long t = Tback; t >= 1; --t) {
                    const double u = tape[static_cast<std::size_t>(t - 1)];
                    detail::cftp_monotone_update(xU, u, logL, logfac, M);
                    detail::cftp_monotone_update(xL, u, logL, logfac, M);
                }
                if (xU == xL) {
                    x = xU;
                    res.horizon[smp] = Tback;
                    break;
                }
                Tback *= 2;
            }
        } else {
            const double eps = 1e-2;
            const long steps = static_cast<long>(std::ceil(
                static_cast<double>(M * (M - 1)) / 2.0 *
                std::log(static_cast<double>(K > 0 ? K : 1) / eps)));
            res.horizon[smp] = steps;
            x.assign(M, 0);
            x[0] = K;
            for (long t = 0; t < steps; ++t) {
                const std::size_t i = static_cast<std::size_t>(mc_uniform_int(rng, M));
                std::size_t j = static_cast<std::size_t>(mc_uniform_int(rng, M - 1));
                if (j >= i) ++j;  // uniform over the M-1 stations other than i
                const int k = x[i] + x[j];
                const int l = detail::cftp_split_index(logL, logfac, i, j, k, mc_uniform01(rng));
                x[i] = l;
                x[j] = k - l;
            }
        }
        for (std::size_t i = 0; i < M; ++i) {
            res.X(smp, i) = x[i];
            res.Q[i] += num_traits<T>::from_int(x[i]);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        res.Q[i] /= num_traits<T>::from_int(static_cast<long>(nsamples));
    return res;
}

/** Reference defaults: single servers, one sample, exact CFTP. */
template <class T>
CftpResult<T> pfqn_cftp(const std::vector<T>& L, int N, McRng& rng) {
    return pfqn_cftp(L, N, std::vector<int>(), static_cast<std::size_t>(1), CftpMethod::Cftp, rng);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CFTP_H
