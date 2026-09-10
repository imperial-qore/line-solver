/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_STDF_HEUR_H
#define LINE_API_PFQN_STDF_HEUR_H

/**
 * Heuristic sojourn-time distribution at multiserver FCFS stations, a variant
 * of J. McKenna, JACM 1987.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_stdf_heur.m. See pfqn_stdf.h for
 * the method; this file documents only where the heuristic departs from it.
 *
 * 1. The level CDF is built PER CLASS. Below the server count it uses that
 *    class's own rate, Exp(rates(k,r)), instead of the common station rate, so
 *    unlike pfqn_stdf the heuristic accepts an FCFS station whose per-class
 *    rates differ and performs no rate-agreement check.
 * 2. At and above the server count the Erlang residual of pfqn_stdf, which
 *    assumes one common rate, is replaced by a sum of R independent
 *    exponentials, one per class, whose means split the n - S(k) + 1 waiting
 *    jobs in proportion to the mean queue lengths Q1(k,s) of the aggregate
 *    solve:  Exp( Q1(k,s) (n - S(k) + 1) / sum_s Q1(k,s) / rates(k,s) ).
 *    Q1 comes from pfqn_mvald, or, when there is a single station, from the
 *    ratio of two pfqn_comomrm_ld constants, the second evaluated on the
 *    aggregate rate lattice pfqn_mu_ms(sum(N), 2, S(k)) of two S(k)-server
 *    stations.
 * 3. The inner constants Y_ks(t) come from the reduction heuristic pfqn_rd, not
 *    from an exact load-dependent solve, and the tilted lattice is NOT
 *    truncated before being handed over.
 * 4. The outer constant lGk always comes from pfqn_mvald, even when there is a
 *    single station and the other constants come from pfqn_comomrm_ld.
 * 5. There is NO min(1, .) clamp on the result. The heuristic therefore reports
 *    values above one where the approximation overshoots -- this is observed,
 *    not hypothetical, and is reproduced deliberately rather than repaired.
 *
 * TWO DEFECTS OF THE REFERENCE, reproduced as failures rather than as wrong
 * numbers. Neither is repaired here, since MATLAB is the ground truth and both
 * are failures in it too.
 *
 *   (a) Small t. h_k(t | n) underflows to zero for the higher levels, so the
 *       tilted rate gamma_k(t, n) becomes infinite, every entry of the beta
 *       transform inside pfqn_rd becomes infinite, and pfqn_rd's
 *       `lastfinite = max(find(isfinite(...)))` is empty, whereupon
 *       `s(ist) = lastfinite` raises "Unable to perform assignment because the
 *       left and right sides have a different number of elements". The ported
 *       pfqn_rd raises NumericError("a station has no finite load-dependent
 *       rate") at the same point, which is the same condition with a diagnosis.
 *       CORRECTED, MEASURED: this note previously claimed the failure fires at
 *       the FineTol substitute for t = 0 "in every model tried", so that the
 *       heuristic could not be evaluated at t = 0 at all. That is false. For
 *       the single-delay family MATLAB returns 1.9999999767e-08 through
 *       1.9047620835e-09 at N = 1 to 5 there and the port matches every one.
 *       The condition is real but model-dependent, not universal.
 *
 *   (b) A class with population zero at N - e_r. Then Q1(k,s) = 0 for that
 *       class and item 2 above asks for Exp with mean 0, i.e. an infinite rate.
 *       MATLAB builds D0 = -Inf, map_cdf returns NaN, the tilt is NaN and
 *       pfqn_rd fails as in (a). This makes the heuristic unusable for every t
 *       on any model with a class whose population is one. The port throws
 *       InputError from map_exponential_mean at the point the degenerate MAP is
 *       requested, which is the earliest place the defect is detectable.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION, gated on has_transcendental, for the
 * reasons given in pfqn_stdf.h and additionally because pfqn_rd is itself a
 * truncated correction series.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_transform.h"
#include "line/api/pfqn/pfqn_comomrm_ld.h"
#include "line/api/pfqn/pfqn_mu_ms.h"
#include "line/api/pfqn/pfqn_mushift.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/api/pfqn/pfqn_rd.h"
#include "line/api/pfqn/pfqn_stdf.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * Heuristic sojourn-time distribution at the listed FCFS stations.
 *
 * Arguments as pfqn_stdf, except that `rates` may differ across classes at an
 * analyzed station: that is the point of the heuristic.
 */
template <class T>
StdfResult<T> pfqn_stdf_heur(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                             const std::vector<int>& S,
                             const std::vector<std::size_t>& fcfsNodes, const Matrix<T>& rates,
                             const std::vector<T>& tset) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_stdf_heur requires transcendental arithmetic: the level CDFs are matrix "
                  "exponentials, the inner constants come from a truncated correction series, and "
                  "the constants are combined in the log domain");

    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    if (R != N.size()) throw InputError("pfqn_stdf_heur: L and N disagree on the class count");
    if (S.size() != M) throw InputError("pfqn_stdf_heur: S has the wrong length");
    if (rates.rows() != M || rates.cols() != R)
        throw InputError("pfqn_stdf_heur: rates has the wrong shape");
    for (std::size_t k = 0; k < M; ++k)
        if (S[k] < 1) throw InputError("pfqn_stdf_heur: the server count must be at least one");

    int Nt = 0;
    for (int n : N) {
        if (n < 0) throw InputError("pfqn_stdf_heur: negative population");
        Nt += n;
    }
    if (Nt < 1) throw InputError("pfqn_stdf_heur: the population must be at least one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T fine = num_traits<T>::from_double(kStdfFineTol);
    const std::size_t nt = tset.size();

    const Matrix<T> mu = detail::stdf_mu<T>(S, Nt);
    const std::vector<T> tv = detail::stdf_guard_tset(tset);

    StdfResult<T> res;
    res.tset = tv;
    res.RD.assign(M, std::vector<Matrix<T>>(R));

    const bool singleStation = (M == 1);

    for (std::size_t ki = 0; ki < fcfsNodes.size(); ++ki) {
        const std::size_t k = fcfsNodes[ki];
        if (k >= M) throw InputError("pfqn_stdf_heur: FCFS station index out of range");
        if (!(rates(k, 0) > zero))
            throw InputError("pfqn_stdf_heur: the FCFS service rate must be strictly positive");

        const Matrix<T> Lk = detail::stdf_drop_row(L, k);
        const Matrix<T> muk = detail::stdf_drop_row(mu, k);

        for (std::size_t r = 0; r < R; ++r) {
            if (!(L(k, r) > fine)) continue;
            std::vector<int> Nr = N;
            Nr[r] -= 1;
            int sumNr = 0;
            for (int n : Nr) sumNr += n;

            // ---- aggregate solve: the constant and the mean queue lengths ----
            double lGr;
            Matrix<T> Q1(M, R, zero);
            if (singleStation) {
                lGr = pfqn_comomrm_ld(L, Nr, Z, mu).lG;
                // Aggregate of two S(k)-server stations, as the reference does.
                const std::vector<T> ms = pfqn_mu_ms<T>(Nt, 2, S[k]);
                Matrix<T> muA(1, ms.size());
                for (std::size_t j = 0; j < ms.size(); ++j) muA(0, j) = ms[j];
                const double lGra = pfqn_comomrm_ld(L, Nr, Z, muA).lG;
                const T ratio = num_traits<T>::from_double(std::exp(lGra - lGr));
                for (std::size_t s = 0; s < R; ++s) Q1(k, s) = L(k, s) * ratio;
            } else {
                const MvaLdResult<T> agg = pfqn_mvald(L, Nr, Z, mu);
                if (!agg.isNumStable) res.isNumStable = false;
                lGr = agg.lG;
                Q1 = agg.QN;
            }

            // ---- per-class level CDFs ---------------------------------------
            T Q1sum = zero;
            for (std::size_t s = 0; s < R; ++s) Q1sum += Q1(k, s);
            Matrix<T> hkc(nt, static_cast<std::size_t>(Nt) + 1);
            for (int n = 0; n <= Nt; ++n) {
                mam::Map<T> h;
                if (!(rates(k, r) > zero))
                    throw InputError("pfqn_stdf_heur: non-positive per-class service rate");
                const T mean = T(one / rates(k, r));
                if (n < S[k]) {
                    h = mam::map_exponential_mean(mean);
                } else {
                    std::vector<mam::Map<T>> parts;
                    parts.push_back(mam::map_exponential_mean(mean));
                    const T lvl = num_traits<T>::from_int(n - S[k] + 1);
                    for (std::size_t s = 0; s < R; ++s)
                        parts.push_back(mam::map_exponential_mean(
                            T(Q1(k, s) * lvl / Q1sum / rates(k, s))));
                    h = mam::map_sumind(parts);
                }
                const std::vector<T> F = mam::map_cdf(h, tv);
                for (std::size_t t = 0; t < nt; ++t) hkc(t, static_cast<std::size_t>(n)) = F[t];
            }

            // ---- outer constant: always the load-dependent MVA ---------------
            const MvaLdResult<T> outer = pfqn_mvald(Lk, Nr, Z, muk);
            const double lGk = outer.lG;

            Matrix<T> RD(nt, 2);
            for (std::size_t t = 0; t < nt; ++t) RD(t, 1) = tv[t];

            Matrix<T> gammat, gammak;
            for (std::size_t t = 0; t < nt; ++t) {
                detail::stdf_gamma(mu, hkc, t, k, sumNr, false, gammat, gammak);
                T H = hkc(t, 0) * num_traits<T>::from_double(std::exp(lGk));
                for (std::size_t s = 0; s < R; ++s) {
                    if (Nr[s] <= 0) continue;
                    std::vector<int> Nrs = Nr;
                    Nrs[s] -= 1;
                    const double lY = pfqn_rd(L, Nrs, Z, gammak).lGN;
                    H += L(k, s) * hkc(t, 0) / gammat(k, 0) *
                         num_traits<T>::from_double(std::exp(lY));
                }
                if (!(H == H)) H = fine;  // the reference's isnan guard
                const double lH = num_traits<T>::log_as_double(H);
                // No min(1, .): the heuristic can and does overshoot one.
                RD(t, 0) = num_traits<T>::from_double(std::exp(lH - lGr));
            }
            res.RD[k][r] = RD;
        }
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_STDF_HEUR_H
