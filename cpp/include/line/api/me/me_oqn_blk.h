/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_ME_ME_OQN_BLK_H
#define LINE_API_ME_ME_OQN_BLK_H

/**
 * Maximum Entropy algorithm for single-class OPEN networks with FINITE BUFFERS,
 * under loss or transfer blocking.
 *
 * Templated port of `matlab/src/api/me/me_oqn_blk.m`: Kouvatsos (1994) Section
 * 4 for the loss case, and Tahilramani, Manjunath and Bose (1999) for transfer
 * blocking.
 *
 * TWO POLICIES, AND ONLY ONE OF THEM IS EASY.
 *
 *   LOSS (blockrule 0). A job finding the destination full is discarded. Every
 *   station is then a censored GE/GE/c/0;N queue and the network is a
 *   straightforward ME decomposition.
 *
 *   TRANSFER BLOCKING (blockrule 1, BAS). A job that completes at i and finds j
 *   full is held in i's SERVER, which can serve nobody else until j has room.
 *   That is NOT WORK CONSERVING, so no product-form approximation applies to
 *   the network as it stands. The reference first makes it work conserving by
 *   inserting a GE/GE/inf HOLDING NODE h_ij on every routing pair with a
 *   finite-buffer destination: the holding node absorbs the blocked job so i's
 *   server is released, and the delay it charges is the residual life of the
 *   minimum of the c_j services in progress at j, inflated geometrically
 *   because the released job may find j full again. Station i's own service is
 *   inflated by the same blocking probability, so the jobs queued behind the
 *   blocked one still see the server busy. THE HELD JOBS ARE ADDED BACK TO
 *   STATION i at the end -- they are physically in i's servers, and reporting
 *   them at the holding node would lose them from the station table.
 *
 * WHY THE SELF-LOOP IS ELIMINATED FIRST. The flow decomposition assumes RENEWAL
 * arrival streams, which immediate feedback breaks. A job returning straight to
 * i receives a geometric number of passes, so the composite service has rate
 * mu(1-p_ii) and scv p_ii + (1-p_ii)Cs; the loop is removed and the residual
 * routing renormalized.
 *
 * TWO CONVENTIONS THAT ARE EASY TO GET WRONG, both carried deliberately.
 * An EXTERNAL arrival finding the buffer full is LOST whatever the drop rule --
 * there is no upstream server to hold it in -- which is both the source's
 * convention and what SolverCTMC does, returning the same answer for DROP and
 * BAS on a source-fed finite queue. And the reported UTILIZATION is LINE's, the
 * carried flow times the nominal mean service time per server, so a server held
 * blocked after service does NOT count as busy; the ME solution's own
 * E[min(n,c)]/c is computed under the INFLATED service and would include the
 * blocking time, so it is recomputed at the end.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/me/me_gegec_mql.h"
#include "line/api/me/me_gegecn.h"
#include "line/api/me/me_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace me {

/** Controls of the blocking fixed point, `me_oqn_blk`'s options struct. */
struct MeBlkOptions {
    double tol = 1e-6;
    long maxiter = 1000;
    /** Relaxation weight on the blocking probabilities, the source's scheme. */
    double damping = 0.5;
};

/** What `me_oqn_blk` returns, per station. */
template <class T>
struct MeBlkResult {
    std::vector<T> Q;       ///< mean number present, the jobs held blocked included
    std::vector<T> W;       ///< mean response time, Q / T
    std::vector<T> T_;      ///< throughput, the CARRIED flow
    std::vector<T> U;       ///< utilization on LINE's convention (see the header)
    std::vector<T> Ca;      ///< interarrival scv of the offered flow
    std::vector<T> Cd;      ///< interdeparture scv
    std::vector<T> PBa;     ///< probability an arrival finds the station full
    std::vector<T> lambda;  ///< OFFERED arrival rate, the attempts included
    long iter = 0;
    bool converged = false;
};

/**
 * Port of `me_oqn_blk`.
 *
 * @param M         number of stations
 * @param lambda0   external arrival rates (M)
 * @param Ca0       external interarrival scv (M), at least 1 where lambda0 > 0
 * @param mu        service rates (M)
 * @param Cs        service scv (M), at least 1 at every finite-buffer station
 * @param P         routing (M x M), `P(i,j) = p_ij`; a row sum below one sends
 *                  the residual flow out of the network
 * @param c         servers per station (M); 0 marks an infinite server, as
 *                  everywhere else in `api/me`
 * @param N         buffer capacity per station (M) in jobs, service included;
 *                  0 marks an unbounded buffer
 * @param blockrule per station (M): 0 = loss, 1 = transfer blocking
 * @param opt       tolerance, iteration budget and relaxation weight
 */
template <class T>
MeBlkResult<T> me_oqn_blk(std::size_t M, const std::vector<T>& lambda0,
                          const std::vector<T>& Ca0, const std::vector<T>& mu,
                          const std::vector<T>& Cs, const Matrix<T>& P,
                          const std::vector<long>& c, const std::vector<long>& N,
                          const std::vector<int>& blockrule,
                          const MeBlkOptions& opt = MeBlkOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "me_oqn_blk requires transcendental arithmetic: its building blocks assemble "
                  "the censored state law in logarithms");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    if (lambda0.size() != M || Ca0.size() != M || mu.size() != M || Cs.size() != M ||
        c.size() != M || N.size() != M || blockrule.size() != M)
        throw InputError("me_oqn_blk: every per-station vector must have M entries");
    if (P.rows() != M || P.cols() != M)
        throw InputError("me_oqn_blk: the routing matrix must be M x M");

    // 0 marks "unbounded" for N and "infinite server" for c, the convention the
    // rest of api/me uses in place of MATLAB's Inf.
    std::vector<char> finiteBuf(M, 0), bas(M, 0);
    for (std::size_t i = 0; i < M; ++i) {
        finiteBuf[i] = (N[i] > 0 && c[i] > 0) ? 1 : 0;
        bas[i] = (finiteBuf[i] && blockrule[i] == 1) ? 1 : 0;
        if (finiteBuf[i] && num_traits<T>::to_double(Cs[i]) < 1.0 - 1e-12)
            throw InputError(
                "me_oqn_blk: a finite buffer requires a service scv of at least 1 at station " +
                std::to_string(i + 1) + ": the GE distribution is not defined for scv < 1.");
        if (lambda0[i] > zero && num_traits<T>::to_double(Ca0[i]) < 1.0 - 1e-12)
            throw InputError(
                "me_oqn_blk: a finite buffer requires an external interarrival scv of at least 1 "
                "at station " + std::to_string(i + 1) + ".");
    }

    // Immediate feedback elimination: the flow decomposition assumes renewal
    // streams, which a self loop breaks.
    Matrix<T> Pf = P;
    std::vector<T> muf = mu, Csf = Cs;
    for (std::size_t i = 0; i < M; ++i) {
        const T pii = P(i, i);
        if (pii > zero) {
            const T q = T(one - pii);
            muf[i] = T(mu[i] * q);
            Csf[i] = T(pii + q * Cs[i]);
            for (std::size_t j = 0; j < M; ++j) Pf(i, j) = T(P(i, j) / q);
            Pf(i, i) = zero;
        }
    }

    // Residual-life rate of the minimum of the c_j services in progress at j. A
    // GE service is zero with probability 1-sigma and exponential with rate
    // mu*sigma otherwise, so its equilibrium residual life is exponential with
    // rate mu*sigma and the minimum over c_j busy servers has rate c*mu*sigma.
    std::vector<T> sigmaS(M, one), muRes(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        sigmaS[i] = T(two / (Csf[i] + one));
        muRes[i] = T(num_traits<T>::from_int(c[i]) * muf[i] * sigmaS[i]);
    }

    MeBlkResult<T> res;
    res.Ca.assign(M, one);
    res.Cd = Csf;
    res.Q.assign(M, zero);
    res.U.assign(M, zero);
    res.T_.assign(M, zero);
    res.lambda.assign(M, zero);
    res.PBa.assign(M, zero);
    Matrix<T> PBs(M, M, zero), PBh(M, M, zero);
    std::vector<T> PBe(M, zero);

    double delta = std::numeric_limits<double>::infinity();
    long it = 0;
    for (it = 1; it <= opt.maxiter; ++it) {
        const std::vector<T> Ca_old = res.Ca;
        const Matrix<T> PBs_old = PBs;
        const std::vector<T> PBe_old = PBe;

        // Service inflation at the blocking stations: the fraction PBf(i) of
        // i's completions is followed by a blocking period.
        std::vector<T> PBf(M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j)
                if (Pf(i, j) > zero && bas[j]) PBf[i] += T(Pf(i, j) * PBs(i, j));
        for (std::size_t i = 0; i < M; ++i)
            if (num_traits<T>::to_double(PBf[i]) >= 1.0 - 1e-9)
                throw NumericError(
                    "me_oqn_blk: the transfer-blocking fixed point saturates, a station is blocked "
                    "with probability one. The network has no stable operating point under BAS.");
        std::vector<T> muEff(M, zero), CsEff(M, one);
        for (std::size_t i = 0; i < M; ++i) {
            muEff[i] = T(muf[i] * (one - PBf[i]));
            CsEff[i] = T(PBf[i] + Csf[i] * (one - PBf[i]));
        }

        // Flow balance on the CARRIED flow. Under loss a fraction PB of a
        // stream is discarded; under transfer blocking every job eventually
        // enters, the delay being charged to the holding node.
        Matrix<T> A(M, M, zero);
        std::vector<T> b(M, zero);
        for (std::size_t j = 0; j < M; ++j) {
            b[j] = T(lambda0[j] * (one - PBe[j]));
            for (std::size_t i = 0; i < M; ++i)
                if (Pf(i, j) > zero)
                    A(i, j) = (finiteBuf[j] && !bas[j]) ? T(Pf(i, j) * (one - PBs(i, j)))
                                                        : Pf(i, j);
        }
        // T = (I - A')^{-1} b
        Matrix<T> S(M, M, zero);
        for (std::size_t r = 0; r < M; ++r) {
            for (std::size_t k = 0; k < M; ++k) S(r, k) = T(-A(k, r));
            S(r, r) = T(S(r, r) + one);
        }
        res.T_ = solve(S, b);
        for (T& v : res.T_)
            if (v < zero) v = zero;

        // Offered (attempt) rates. A job blocked under BAS re-attempts from the
        // holding node, so its stream contributes carried/(1-PB) attempts.
        std::vector<T> attExt(M, zero);
        Matrix<T> attInt(M, M, zero);
        const T tiny = num_traits<T>::from_double(1e-12);
        for (std::size_t j = 0; j < M; ++j) {
            attExt[j] = lambda0[j];
            for (std::size_t i = 0; i < M; ++i)
                if (Pf(i, j) > zero) {
                    if (finiteBuf[j] && bas[j]) {
                        T den = T(one - PBs(i, j));
                        if (den < tiny) den = tiny;
                        attInt(i, j) = T(res.T_[i] * Pf(i, j) / den);
                    } else {
                        attInt(i, j) = T(res.T_[i] * Pf(i, j));
                    }
                }
        }
        for (std::size_t j = 0; j < M; ++j) {
            T tot = attExt[j];
            for (std::size_t i = 0; i < M; ++i) tot += attInt(i, j);
            res.lambda[j] = tot;
            if (tot > zero) {
                T acc = T(attExt[j] * PBe[j]);
                for (std::size_t i = 0; i < M; ++i) acc += T(attInt(i, j) * PBs(i, j));
                res.PBa[j] = T(acc / tot);
            } else {
                res.PBa[j] = zero;
            }
        }

        // Interarrival scv of the offered flow, by GE splitting and merging: a
        // stream thinned with probability p has scv 1-p+p*Cd, and the merge
        // satisfies 1/(Cm+1) = sum_s (lam_s/lam)/(Cs+1).
        Matrix<T> CaStreamInt(M, M, one);
        for (std::size_t j = 0; j < M; ++j) {
            if (!(res.lambda[j] > zero)) continue;
            T sum_inv = zero;
            if (attExt[j] > zero)
                sum_inv += T((attExt[j] / res.lambda[j]) / (Ca0[j] + one));
            for (std::size_t i = 0; i < M; ++i)
                if (attInt(i, j) > zero) {
                    CaStreamInt(i, j) = T(one - Pf(i, j) + Pf(i, j) * res.Cd[i]);
                    sum_inv += T((attInt(i, j) / res.lambda[j]) / (CaStreamInt(i, j) + one));
                }
            if (sum_inv > zero) res.Ca[j] = T(-one + one / sum_inv);
        }

        // Station solution in isolation, and the per-STREAM blocking
        // probabilities: eq. (4.3) evaluated with each stream's own scv.
        std::vector<T> PBe_new(M, zero);
        Matrix<T> PBs_new(M, M, zero), PBh_new(M, M, zero);
        for (std::size_t j = 0; j < M; ++j) {
            if (c[j] <= 0) {
                // Infinite server: no queueing and no blocking.
                res.Q[j] = (muEff[j] > zero) ? T(res.lambda[j] / muEff[j]) : zero;
                res.U[j] = res.Q[j];
                res.Cd[j] = res.Ca[j];
                continue;
            }
            const T cT = num_traits<T>::from_int(c[j]);
            if (!finiteBuf[j]) {
                // Unbounded buffer: the infinite-capacity GE blocks.
                T rho = zero;
                if (muEff[j] > zero) rho = T(res.lambda[j] / (cT * muEff[j]));
                if (num_traits<T>::to_double(rho) >= 1.0) {
                    res.Q[j] = num_traits<T>::from_double(
                        std::numeric_limits<double>::infinity());
                    res.U[j] = one;
                    res.Cd[j] = CsEff[j];
                } else if (c[j] == 1) {
                    res.Q[j] = T(rho * (res.Ca[j] + one) / two +
                                 rho * rho * (CsEff[j] + res.Ca[j]) / (two * (one - rho)));
                    res.U[j] = rho;
                    res.Cd[j] = T(rho * rho * CsEff[j] + (one - rho) * res.Ca[j] +
                                  rho * (one - rho));
                } else {
                    res.Q[j] = me_gegec_mql(res.lambda[j], res.Ca[j], muEff[j], CsEff[j], c[j]);
                    res.U[j] = rho;
                    res.Cd[j] = T(rho * rho * CsEff[j] + (one - rho) * res.Ca[j] +
                                  rho * (one - rho));
                }
                continue;
            }
            const GegecnResult<T> gj =
                me_gegecn(res.lambda[j], res.Ca[j], muEff[j], CsEff[j], c[j], 0, N[j]);
            res.Q[j] = gj.L;
            res.U[j] = gj.U;
            // Interdeparture scv at the CENSORED utilization: with losses the
            // offered load can exceed one while the busy fraction cannot.
            res.Cd[j] = T(gj.U * gj.U * CsEff[j] + (one - gj.U) * res.Ca[j] +
                          gj.U * (one - gj.U));
            if (attExt[j] > zero)
                PBe_new[j] = me_gegecn_pb(gj.p, 0L, N[j], c[j], CsEff[j], Ca0[j]);
            for (std::size_t i = 0; i < M; ++i)
                if (attInt(i, j) > zero) {
                    PBs_new(i, j) =
                        me_gegecn_pb(gj.p, 0L, N[j], c[j], CsEff[j], CaStreamInt(i, j));
                    if (bas[j]) {
                        // The flow released by h_ij is a Bernoulli sample of i's
                        // departures with probability p_ij PB^i_j, carried
                        // through a GE/GE/inf queue whose interdeparture scv
                        // equals its interarrival scv.
                        const T q = T(Pf(i, j) * PBs(i, j));
                        const T CaH = T(one - q + q * res.Cd[i]);
                        PBh_new(i, j) = me_gegecn_pb(gj.p, 0L, N[j], c[j], CsEff[j], CaH);
                    }
                }
        }

        // Relaxation on the blocking probabilities.
        const T w = num_traits<T>::from_double(opt.damping);
        const T w1 = T(one - w);
        for (std::size_t j = 0; j < M; ++j) PBe[j] = T(w1 * PBe[j] + w * PBe_new[j]);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) {
                PBs(i, j) = T(w1 * PBs(i, j) + w * PBs_new(i, j));
                PBh(i, j) = T(w1 * PBh(i, j) + w * PBh_new(i, j));
            }

        delta = 0.0;
        for (std::size_t j = 0; j < M; ++j) {
            delta = std::max(delta, std::fabs(num_traits<T>::to_double(res.Ca[j]) -
                                              num_traits<T>::to_double(Ca_old[j])));
            delta = std::max(delta, std::fabs(num_traits<T>::to_double(PBe[j]) -
                                              num_traits<T>::to_double(PBe_old[j])));
            for (std::size_t i = 0; i < M; ++i)
                delta = std::max(delta, std::fabs(num_traits<T>::to_double(PBs(i, j)) -
                                                  num_traits<T>::to_double(PBs_old(i, j))));
        }
        if (delta < opt.tol) {
            res.converged = true;
            break;
        }
    }
    res.iter = std::min(it, opt.maxiter);

    // Holding-node occupancy. The jobs held in h_ij are physically blocked in
    // the SERVERS OF STATION i, so they are added back to i rather than
    // reported at a node the station table does not have.
    for (std::size_t i = 0; i < M; ++i) {
        T held = zero;
        for (std::size_t j = 0; j < M; ++j)
            if (bas[j] && Pf(i, j) > zero && PBs(i, j) > zero) {
                const T rateH = T(res.T_[i] * Pf(i, j) * PBs(i, j));
                const T muH = T(muRes[j] * (one - PBh(i, j)));
                if (muH > zero) held += T(rateH / muH);
            }
        res.Q[i] += held;
    }

    // Utilization on LINE's convention: the carried flow times the NOMINAL mean
    // service time per server, so a server held blocked after service does not
    // count as busy. The ME solution reports E[min(n,c)]/c under the INFLATED
    // service, which would include the blocking time.
    for (std::size_t i = 0; i < M; ++i) {
        if (muf[i] > zero) {
            res.U[i] = (c[i] <= 0) ? T(res.T_[i] / muf[i])
                                   : T(res.T_[i] / (num_traits<T>::from_int(c[i]) * muf[i]));
        } else {
            res.U[i] = zero;
        }
    }
    res.W.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        if (res.T_[i] > zero) res.W[i] = T(res.Q[i] / res.T_[i]);
    return res;
}

}  // namespace me
}  // namespace line

#endif  // LINE_API_ME_ME_OQN_BLK_H
