/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_CONV_H
#define LINE_SOLVERS_NC_SOLVER_NC_CONV_H

/**
 * Exact convolution analysis of a closed network with class-dependent service
 * rates. Port of `solver_nc_conv.m`.
 *
 * A class-dependent station scales its nominal demand by beta_{i,r}(n), a
 * function of the per-class population AT THAT STATION. None of the
 * normalizing-constant algorithms in `solver_ncld` can apply it -- they read a
 * rate lattice mu(n) and nothing else -- so a cdscaling model is routed here
 * UNCONDITIONALLY on the method, not only on 'exact'. The multichain convolution
 * of Sauer (1983), Sect. 5.2, eq. (40) is the exact algorithm for beta, and it
 * is what `pfqn_conv` implements.
 *
 * QUEUE LENGTHS COME FROM THE MARGINAL, NOT FROM MVA. The station factor
 * X_m(n) is built by the chain-dependent recurrence, and the marginal
 * P_m(n | N) = X_m(n) G_{-m}(N-n) / G(N) is summed against n_k. G_{-m} is a
 * fresh convolution over every station but m, so the cost is one convolution per
 * lattice point per class-dependent station: this analyzer is exact and slow,
 * which is the trade the reference makes.
 *
 * UTILIZATION IS NORMALIZED BY THE DECLARED PEAK. T*ST measures capacity used in
 * units of the NOMINAL rate and reaches max_n beta, not 1, at saturation, so a
 * beta emulating two servers would report U = 2 * (true utilization). Dividing
 * by `sn.cdscalingpeak` restores the T*S/c convention of an ordinary multiserver
 * station. No cap is needed: sum_r U(i,r) is a convex combination of the
 * beta_r(n) over max_n beta, hence at most one by construction.
 *
 * Arithmetic: the convolution itself is field arithmetic, but lG is a log and
 * the analyzer is guarded on `has_transcendental` for it.
 */

#include <cmath>
#include <vector>

#include "line/api/pfqn/pfqn_conv.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

namespace detail {

/** Advance a population vector over 0 <= n <= N, as `pprod_next`; false at end. */
inline bool conv_pprod_next(std::vector<int>& n, const std::vector<int>& N) {
    std::size_t s = n.size();
    while (s > 0 && n[s - 1] == N[s - 1]) {
        n[s - 1] = 0;
        --s;
    }
    if (s == 0) return false;
    n[s - 1] += 1;
    return true;
}

/** Column-major linear index of a population vector, as `hashpop`. */
inline std::size_t conv_hashpop(const std::vector<int>& n, const std::vector<int>& N) {
    std::size_t idx = 0, stride = 1;
    for (std::size_t r = 0; r < N.size(); ++r) {
        idx += stride * static_cast<std::size_t>(n[r]);
        stride *= static_cast<std::size_t>(N[r]) + 1;
    }
    return idx;
}

}  // namespace detail

/**
 * Port of `solver_nc_conv.m`.
 *
 * @param sn  the refreshed struct; at least one station carries a cdscaling
 * @param opt solver controls; unused, the convolution is exact and has no tuning
 */
template <class T>
NcSolution<T> solver_nc_conv(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    (void)opt;
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        throw UnsupportedError(
            "solver_nc_conv: the convolution analyzer reports lG = log(G); this backend has no "
            "transcendental arithmetic");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations;

        // The convolution runs on CHAINS, not classes: a class-switching model
        // splits one circulating population across several classes, so the
        // per-class populations are zero in the classes that hold no reference
        // jobs and the class-level recursion charges those stations nothing at
        // all. Every other closed NC and MVA path aggregates the same way and
        // deaggregates at the end.
        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const std::size_t K = sn.nchains;

        std::vector<int> NK(K, 0);
        for (std::size_t c = 0; c < K; ++c) {
            if (std::isinf(d.Nchain[c]))
                throw UnsupportedError(
                    "solver_nc_conv: the multichain convolution requires a closed queueing "
                    "network; this model has an open class");
            NK[c] = static_cast<int>(std::llround(d.Nchain[c]));
        }

        const Matrix<T>& V = d.Vchain;
        const Matrix<T>& ST = d.STchain;
        const Matrix<T>& Ldemand = d.Lchain;

        std::vector<std::size_t> delayIdx, queueIdx;
        for (std::size_t i = 0; i < M; ++i)
            (std::isinf(sn.stations[i].nservers) ? delayIdx : queueIdx).push_back(i + 1);
        const std::size_t nQueues = queueIdx.size();

        Matrix<T> Z_conv(1, K, zero);
        for (std::size_t i : delayIdx)
            for (std::size_t r = 0; r < K; ++r) Z_conv(0, r) = T(Z_conv(0, r) + Ldemand(i - 1, r));

        Matrix<T> L_conv(nQueues, K, zero);
        for (std::size_t q = 0; q < nQueues; ++q)
            for (std::size_t r = 0; r < K; ++r) L_conv(q, r) = Ldemand(queueIdx[q] - 1, r);

        std::vector<lang::CdScaling<T>> cd_conv(nQueues);
        for (std::size_t q = 0; q < nQueues; ++q)
            cd_conv[q] = sn.stations[queueIdx[q] - 1].cdscaling;

        const T G_N = pfqn::pfqn_conv(L_conv, NK, Z_conv, cd_conv).G;
        out.sol.lG = num_traits<T>::log_as_double(G_N);

        std::vector<T> XN(K, zero);
        for (std::size_t r = 0; r < K; ++r)
            if (NK[r] > 0) {
                std::vector<int> Nm = NK;
                --Nm[r];
                const T G_Nk = pfqn::pfqn_conv(L_conv, Nm, Z_conv, cd_conv).G;
                XN[r] = T(G_Nk / G_N);
            }

        Matrix<T> TN(M, K, zero), QN(M, K, zero), RN(M, K, zero), UN(M, K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) TN(i, r) = T(V(i, r) * XN[r]);

        for (std::size_t i : delayIdx)
            for (std::size_t r = 0; r < K; ++r)
                QN(i - 1, r) = T(Ldemand(i - 1, r) * XN[r]);

        std::size_t stateSpaceSize = 1;
        for (int v : NK) stateSpaceSize *= static_cast<std::size_t>(v) + 1;

        for (std::size_t q = 0; q < nQueues; ++q) {
            const std::size_t ist = queueIdx[q];
            const bool isCd = static_cast<bool>(cd_conv[q]);

            // the station factor X_m(n), by the chain-dependent recurrence
            std::vector<T> Xm(stateSpaceSize, zero);
            Xm[0] = num_traits<T>::from_int(1);
            std::vector<int> n(K, 0);
            do {
                int tot = 0;
                for (int v : n) tot += v;
                if (tot == 0) continue;
                const std::size_t idx = detail::conv_hashpop(n, NK);
                if (isCd) {
                    // X_m(n) = (|n|/n_r) (L/beta_r(n)) X_m(n - e_r) for ANY r with
                    // n_r > 0; the result is path independent, so the first is used
                    for (std::size_t r = 0; r < K; ++r) {
                        if (n[r] == 0) continue;
                        std::vector<T> row(K, zero);
                        for (std::size_t j = 0; j < K; ++j)
                            row[j] = num_traits<T>::from_int(n[j]);
                        const std::vector<T> bval = cd_conv[q](row);
                        if (bval.empty())
                            throw UnsupportedError(
                                "solver_nc_conv: a class-dependence map returned no value");
                        const T beta = bval.size() > 1 ? bval[r] : bval[0];
                        const int nr = n[r];
                        n[r] -= 1;
                        const std::size_t idx_prev = detail::conv_hashpop(n, NK);
                        n[r] += 1;
                        if (beta > zero)
                            Xm[idx] = T(num_traits<T>::from_int(tot) /
                                        num_traits<T>::from_int(nr) * (L_conv(q, r) / beta) *
                                        Xm[idx_prev]);
                        break;
                    }
                } else {
                    for (std::size_t r = 0; r < K; ++r) {
                        if (n[r] == 0) continue;
                        n[r] -= 1;
                        const std::size_t idx_prev = detail::conv_hashpop(n, NK);
                        n[r] += 1;
                        Xm[idx] = T(Xm[idx] + L_conv(q, r) * Xm[idx_prev]);
                    }
                }
            } while (detail::conv_pprod_next(n, NK));

            // the complement network, every station but this one
            Matrix<T> L_comp(nQueues > 0 ? nQueues - 1 : 0, K, zero);
            std::vector<lang::CdScaling<T>> cd_comp;
            for (std::size_t j = 0, w = 0; j < nQueues; ++j) {
                if (j == q) continue;
                for (std::size_t r = 0; r < K; ++r) L_comp(w, r) = L_conv(j, r);
                cd_comp.push_back(cd_conv[j]);
                ++w;
            }

            std::vector<int> m(K, 0);
            do {
                bool anyPos = false;
                for (int v : m)
                    if (v > 0) anyPos = true;
                if (!anyPos) continue;
                const std::size_t idx = detail::conv_hashpop(m, NK);
                std::vector<int> nmi(K, 0);
                for (std::size_t r = 0; r < K; ++r) nmi[r] = NK[r] - m[r];
                const T G_comp = pfqn::pfqn_conv(L_comp, nmi, Z_conv, cd_comp).G;
                const T prob = T(Xm[idx] * G_comp / G_N);
                for (std::size_t r = 0; r < K; ++r)
                    QN(ist - 1, r) = T(QN(ist - 1, r) + num_traits<T>::from_int(m[r]) * prob);
            } while (detail::conv_pprod_next(m, NK));
        }

        // RN is the PER-VISIT response time Qchain/Tchain: the deaggregation
        // below multiplies the visit ratio back in, so dividing by Xchain would
        // count it twice
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                if (TN(i, r) != zero) RN(i, r) = T(QN(i, r) / TN(i, r));
                UN(i, r) = T(TN(i, r) * ST(i, r));
            }

        for (std::size_t q = 0; q < nQueues; ++q) {
            if (!cd_conv[q]) continue;
            const std::size_t ist = queueIdx[q];
            const std::vector<T>& peak = sn.stations[ist - 1].cdscalingpeak;
            if (peak.empty())
                throw UnsupportedError(
                    "solver_nc_conv: a class-dependent station has no declared peak rate; pass "
                    "peakRatePerClass to setClassDependence");
            for (std::size_t c = 0; c < K; ++c) {
                // The peaks are declared per class, so the chain takes the
                // largest peak among its classes: utilization is a per-station
                // quantity with one normalizer.
                T bmax = zero;
                for (std::size_t kk : sn.inchain[c])
                    if (peak[kk - 1] > bmax) bmax = peak[kk - 1];
                if (bmax > zero) UN(ist - 1, c) = T(UN(ist - 1, c) / bmax);
            }
        }

        const mva::ClassResults<T> cls =
            mva::sn_deaggregate_chain_results(sn, d, Matrix<T>(), UN, RN, TN, XN);

        out.sol.Q = cls.Q;
        out.sol.U = cls.U;
        out.sol.R = cls.R;
        out.sol.Tp = cls.Tp;
        out.sol.X = cls.X;
        out.sol.C = cls.C;
        out.sol.iter = 1;
        out.sol.method = "conv";
        out.actualmethod = "conv";
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_CONV_H
