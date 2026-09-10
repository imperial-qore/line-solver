/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_AB_AMVA_H
#define LINE_API_PFQN_AB_AMVA_H

/**
 * Akyildiz-Bolch approximate MVA for multi-server BCMP networks.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ab_amva.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_ab_amva.java.
 *
 * The driver is a Linearizer-shaped three-pass scheme around a fixed-point
 * core:
 *
 *   step 1  run the core at the full population N, from the flat start
 *           L(i,r) = N_r / M;
 *   step 2  run the core once per class at N - e_r, seeding it with the
 *           per-class queue lengths of step 1;
 *   step 3  form the fractional-change tensor
 *           Delta(i,r,t) = Q(i,t | N - e_r)/den - Q(i,r | N)/N_r,
 *           with den = N_r - 1 when r == t and N_r otherwise;
 *   step 4  rerun the core at N with those Delta held fixed.
 *
 * Inside the core the arrival-theorem queue length seen by a class-r job is
 * L(i,c | N - e_r) = scalar (F(i,c) + Delta(i,c,r)), and the residence time is
 * formed per station kind:
 *
 *   INF                 W = S(i,r);
 *   single server       W = S(i,r) (1 + sum_c L(i,c | N - e_r));
 *   multiserver         W = S(i,r)/c (1 + Qtot + sum_{j<c} (c-j) Pr(j)),
 *                       with Pr(j) the Akyildiz-Bolch marginal weight (or the
 *                       two-point 'scat' scatter) of the queue length;
 *   FCFS, fcfsSchmidt   W = sum_{n <= N, n_r > 0} B_r(n) Pr(n - e_r),
 *                       with Pr a per-class binomial product and B_r the
 *                       queue-composition-weighted mean service time.
 *
 * Reference behaviour preserved verbatim, including the parts that look like
 * defects but define the numbers the reference produces:
 *
 *  - the throughput is read off STATION 1, XN(r) = Q(1,r)/W(1,r), so it is the
 *    class-r throughput AT that station, i.e. v(1,r) times the system
 *    throughput, not the system throughput itself;
 *  - the convergence tolerance is the reference's 1/(4000 + 16 sum(N)) and the
 *    iteration cap is 100 passes, with no error on non-convergence;
 *  - a Schmidt-FCFS wait below 1e-3 is snapped to zero;
 *  - the AB marginal weights use the reference's ALPHA = 45, BETA = 0.7 and
 *    its distance cutoff of 25.
 *
 * One sizing divergence, deliberate. weightFun builds its weight table with
 * max(N) + 1 rows, but the row index used later is floor(Qtot) where Qtot is
 * the TOTAL queue length over all classes, which exceeds max(N) as soon as two
 * classes are both loaded at a multiserver station; MATLAB then raises an
 * index-out-of-bounds error. The recursion defining the table is index-generic
 * (row l is built from row l-1 and the geometric scaling sequence), so the port
 * sizes the table by the largest row actually required. That evaluates the SAME
 * function at a larger argument; it changes no in-range entry and it removes an
 * error the reference cannot otherwise avoid.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The core is a tolerance-stopped fixed
 * point, so the answer depends on where the iteration is cut; the marginal
 * weights additionally need floor and integer powers of a non-integer
 * fraction. It is gated on has_transcendental accordingly.
 */

#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** Which marginal-probability rule the multiserver correction uses. */
enum class AbMarginalMethod {
    Ab,   ///< the Akyildiz-Bolch weight function
    Scat  ///< two-point scatter around floor(Qtot)
};

/** Return value of pfqn_ab_amva, mirroring [QN,UN,RN,CN,XN,totiter]. */
template <class T>
struct AbAmvaResult {
    Matrix<T> QN;        ///< (M x R) mean queue length
    Matrix<T> UN;        ///< (M x R) utilization
    Matrix<T> RN;        ///< (M x R) residence (wait) time per visit
    std::vector<T> CN;   ///< (R) cycle time
    std::vector<T> XN;   ///< (R) class throughput AT STATION 1
    std::size_t totiter; ///< iterations of the final core pass
};

namespace detail {

/** Akyildiz-Bolch weight table, rows 0..lmax, w(l,j) defined for j <= l. */
template <class T>
Matrix<T> ab_weight_fun(const std::vector<int>& population, int lmax, double alpha, double beta) {
    int maxPop = 0;
    for (int n : population)
        if (n > maxPop) maxPop = n;
    if (lmax > maxPop) maxPop = lmax;
    if (maxPop < 0) maxPop = 0;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t n1 = static_cast<std::size_t>(maxPop) + 1;

    std::vector<T> scaling(n1, zero);
    if (maxPop >= 1) {
        scaling[1] = num_traits<T>::from_double(alpha);
        for (int n = 2; n <= maxPop; ++n)
            scaling[static_cast<std::size_t>(n)] =
                num_traits<T>::from_double(beta) * scaling[static_cast<std::size_t>(n) - 1];
    }
    Matrix<T> w(n1, n1, zero);
    w(0, 0) = one;
    const T hundred = num_traits<T>::from_int(100);
    for (int l = 1; l <= maxPop; ++l) {
        const std::size_t lu = static_cast<std::size_t>(l);
        T sum = zero;
        for (int j = 0; j <= l - 1; ++j) {
            const std::size_t ju = static_cast<std::size_t>(j);
            w(lu, ju) = w(lu - 1, ju) - w(lu - 1, ju) * scaling[lu] / hundred;
            sum += w(lu, ju);
        }
        w(lu, lu) = one - sum;
    }
    return w;
}

/** matlab findMarginalProbs, as an integer-keyed map (keys may be negative). */
template <class T>
std::map<int, T> ab_marginal_probs(const T& avgJobs, int numServers,
                                   const std::vector<int>& population, std::size_t classIdx,
                                   AbMarginalMethod method) {
    const T zero = num_traits<T>::from_int(0);
    std::map<int, T> mp;
    const double aj = num_traits<T>::to_double(avgJobs);
    const int floorVal = static_cast<int>(std::floor(aj));

    if (method == AbMarginalMethod::Scat) {
        const int ceilVal = floorVal + 1;
        mp[floorVal] = num_traits<T>::from_int(ceilVal) - avgJobs;
        mp[ceilVal] = avgJobs - num_traits<T>::from_int(floorVal);
        return mp;
    }

    const int ceiling = floorVal + 1;
    const int maxVal = (2 * floorVal + 1) < (numServers - 2) ? (2 * floorVal + 1) : (numServers - 2);
    const Matrix<T> w = ab_weight_fun<T>(population, floorVal < 0 ? 0 : floorVal, 45.0, 0.7);
    const std::size_t wn = w.rows();
    const auto wat = [&](int l, int j) -> T {
        if (l < 0 || j < 0 || static_cast<std::size_t>(l) >= wn ||
            static_cast<std::size_t>(j) >= wn)
            return num_traits<T>::from_int(0);
        return w(static_cast<std::size_t>(l), static_cast<std::size_t>(j));
    };
    const int popc = population[classIdx];

    for (int j = 0; j <= maxVal; ++j) {
        if (j <= floorVal) {
            const int lDist = floorVal - j;
            const int lowerVal = floorVal - lDist;
            const int upperVal = ceiling + lDist;
            T prob = zero;
            if (lDist <= 25 && floorVal < popc && upperVal != lowerVal)
                prob = wat(floorVal, lDist) *
                       ((num_traits<T>::from_int(upperVal) - avgJobs) /
                        num_traits<T>::from_int(upperVal - lowerVal));
            mp[j] = prob;
        } else {
            const int uDist = j - ceiling;
            if (uDist > 25) {
                mp[j] = zero;
            } else if (j > popc - 1 && uDist < 25) {
                // uDist == 25 falls through to the plain branch in the
                // reference, which keys the entry on j rather than on popc-1
                const auto ite = mp.find(popc - 1);
                const T existing = ite == mp.end() ? zero : ite->second;
                const auto itf = mp.find(floorVal - uDist);
                const T mfu = itf == mp.end() ? zero : itf->second;
                mp[popc - 1] = existing + (wat(floorVal, uDist) - mfu);
            } else {
                const auto itf = mp.find(floorVal - uDist);
                const T mfu = itf == mp.end() ? zero : itf->second;
                mp[j] = wat(floorVal, uDist) - mfu;
            }
        }
    }
    return mp;
}

/** matlab getBcnForAB: queue-composition-weighted mean service time. */
template <class T>
T ab_bcn(const Matrix<T>& S, std::size_t i, std::size_t c, const std::vector<int>& nvec, int ns) {
    T bcn = S(i, c);
    long nsum = 0;
    for (int t : nvec) nsum += t;
    if (nsum > 1) {
        const T eps = num_traits<T>::from_double(1e-12);
        T sumVal = num_traits<T>::from_int(0);
        for (std::size_t t = 0; t < nvec.size(); ++t)
            sumVal += num_traits<T>::from_int(nvec[t]) * S(i, t);
        const T num = num_traits<T>::from_int(nsum - ns > 0 ? nsum - ns : 0);
        const T den0 = num_traits<T>::from_int(ns * (nsum - 1));
        const T den = den0 > eps ? den0 : eps;
        bcn += num / den * (sumVal - S(i, c));
    }
    return bcn;
}

/** matlab getMarginalProb: per-class binomial product. */
template <class T>
T ab_binomial_prob(const std::vector<int>& n, const std::vector<int>& Kpop, const T& Ljr) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    T prob = one;
    for (std::size_t r = 0; r < Kpop.size(); ++r) {
        if (Kpop[r] <= 0) continue;
        const T frac = Ljr / num_traits<T>::from_int(Kpop[r]);
        if (frac == zero) continue;
        if (n[r] < 0 || n[r] > Kpop[r]) continue;
        const T t1 = num_nck<T>(Kpop[r], n[r]);
        const T t2 = num_pow_int(frac, static_cast<unsigned>(n[r]));
        const T t3 = num_pow_int(T(one - frac), static_cast<unsigned>(Kpop[r] - n[r]));
        prob *= t1 * t2 * t3;
    }
    return prob;
}

/** The Akyildiz-Bolch fixed-point core, matlab pfqn_ab_core. */
template <class T>
AbAmvaResult<T> ab_core(const std::vector<int>& population, const std::vector<int>& nservers,
                        const std::vector<SchedStrategy>& type, const Matrix<T>& v,
                        const Matrix<T>& S, std::size_t maxiter, const std::vector<T>& Delta,
                        const Matrix<T>& lIn, bool fcfsSchmidt, AbMarginalMethod method) {
    const std::size_t M = S.rows(), K = S.cols();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const auto dat = [&](std::size_t i, std::size_t r, std::size_t t) -> const T& {
        return Delta[(i * K + r) * K + t];
    };

    long Npop = 0;
    for (int n : population) Npop += n;
    const T tol = one / num_traits<T>::from_int(4000 + 16 * Npop);

    AbAmvaResult<T> res;
    res.QN = lIn;
    res.RN = Matrix<T>(M, K, zero);
    res.UN = Matrix<T>(M, K, zero);
    res.CN.assign(K, zero);
    res.XN.assign(K, zero);
    res.totiter = 0;

    Matrix<T>& L = res.QN;
    Matrix<T>& W = res.RN;
    Matrix<T> F(M, K, zero);
    std::vector<T> lWJ(M * K * K, zero);
    const auto lwj = [&](std::size_t i, std::size_t r, std::size_t t) -> T& {
        return lWJ[(i * K + r) * K + t];
    };
    const T milli = num_traits<T>::from_double(1e-3);

    while (res.totiter < maxiter) {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                F(i, r) = population[r] > 0 ? T(L(i, r) / num_traits<T>::from_int(population[r]))
                                            : zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t t = 0; t < K; ++t) {
                    const long scalar = r == t ? population[r] - 1 : population[r];
                    lwj(i, r, t) = num_traits<T>::from_int(scalar) * (F(i, r) + dat(i, r, t));
                }

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                if (type[i] == SchedStrategy::INF) {
                    W(i, r) = S(i, r);
                } else if (nservers[i] == 1) {
                    T qtot = zero;
                    for (std::size_t c = 0; c < K; ++c) qtot += lwj(i, c, r);
                    W(i, r) = S(i, r) * (one + qtot);
                } else if (fcfsSchmidt && type[i] == SchedStrategy::FCFS) {
                    T wait = zero;
                    std::vector<int> nvec(K, 0);
                    bool more = true;
                    while (more) {
                        if (nvec[r] > 0) {
                            const T bcn = ab_bcn(S, i, r, nvec, nservers[i]);
                            std::vector<int> nm = nvec;
                            nm[r] -= 1;
                            // seeded-queue-length binomial rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                            wait += bcn * ab_binomial_prob<T>(nm, population, lIn(i, r));
                        }
                        more = next_pop(nvec, population);
                    }
                    W(i, r) = wait <= milli ? zero : wait;
                } else {
                    T qtot = zero;
                    for (std::size_t j = 0; j < K; ++j) qtot += lwj(i, j, r);
                    const int c = nservers[i];
                    T corr = zero;
                    if (c > 1) {
                        std::vector<int> popWithoutR = population;
                        popWithoutR[r] -= 1;
                        const std::map<int, T> mp =
                            ab_marginal_probs<T>(qtot, c, popWithoutR, r, method);
                        for (int j = 1; j <= c - 1; ++j) {
                            const auto it = mp.find(j);
                            if (it != mp.end())
                                corr += it->second * num_traits<T>::from_int(c - j);
                        }
                    }
                    W(i, r) = S(i, r) / num_traits<T>::from_int(c) * (one + qtot + corr);
                }
            }

        for (std::size_t r = 0; r < K; ++r) {
            T cyc = zero;
            for (std::size_t i = 0; i < M; ++i) cyc += v(i, r) * W(i, r);
            res.CN[r] = cyc;
        }

        Matrix<T> itQ(M, K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                itQ(i, r) = res.CN[r] > zero
                                ? T(num_traits<T>::from_int(population[r]) *
                                    (v(i, r) * W(i, r) / res.CN[r]))
                                : zero;

        T maxDiff = zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                if (population[r] <= 0) continue;
                const T diff =
                    num_abs(T(L(i, r) - itQ(i, r))) / num_traits<T>::from_int(population[r]);
                if (diff > maxDiff) maxDiff = diff;
            }

        res.totiter += 1;
        L = itQ;
        if (maxDiff < tol) break;
    }

    for (std::size_t r = 0; r < K; ++r)
        res.XN[r] = W(0, r) > zero ? T(L(0, r) / W(0, r)) : zero;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!(S(i, r) > zero)) continue;
            res.UN(i, r) = type[i] == SchedStrategy::INF
                               ? T(res.XN[r] * S(i, r))
                               : T(res.XN[r] * S(i, r) / num_traits<T>::from_int(nservers[i]));
        }
    return res;
}

}  // namespace detail

/**
 * @param S         (M x R) service demands
 * @param N         (R) population per class
 * @param v         (M x R) visit ratios
 * @param nservers  (M) server counts
 * @param sched     (M) scheduling discipline
 * @param fcfsSchmidt use the Schmidt state-sum wait at FCFS stations
 * @param method    marginal-probability rule for the multiserver correction
 */
template <class T>
AbAmvaResult<T> pfqn_ab_amva(const Matrix<T>& S, const std::vector<int>& N, const Matrix<T>& v,
                             const std::vector<int>& nservers,
                             const std::vector<SchedStrategy>& sched, bool fcfsSchmidt,
                             AbMarginalMethod method) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_ab_amva requires transcendental arithmetic: its core is a "
                  "tolerance-stopped fixed point, so the answer depends on where the iteration "
                  "stops, and its marginal weights use floor and non-integer fractions");

    const std::size_t M = S.rows(), K = S.cols();
    if (N.size() != K) throw InputError("pfqn_ab_amva: S and N disagree on the class count");
    if (v.rows() != M || v.cols() != K)
        throw InputError("pfqn_ab_amva: visit-ratio matrix has the wrong shape");
    if (nservers.size() != M) throw InputError("pfqn_ab_amva: nservers has the wrong length");
    if (sched.size() != M) throw InputError("pfqn_ab_amva: sched has the wrong length");
    for (std::size_t i = 0; i < M; ++i)
        if (nservers[i] < 1) throw InputError("pfqn_ab_amva: server count below one");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_ab_amva: negative population");

    const T zero = num_traits<T>::from_int(0);
    const std::size_t maxiter = 100;

    // Flat start L(i,r) = N_r / M, and the per-class seeds of step 2.
    Matrix<T> L(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            L(i, r) = num_traits<T>::from_int(N[r]) / num_traits<T>::from_int(static_cast<long>(M));

    std::vector<Matrix<T>> lWithoutR(K, Matrix<T>(M, K, zero));
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t t = 0; t < K; ++t)
                lWithoutR[r](i, t) =
                    r == t ? T(num_traits<T>::from_int(N[r] - 1) /
                               num_traits<T>::from_int(static_cast<long>(M)))
                           : L(i, r);

    std::vector<T> Delta(M * K * K, zero);

    // STEP 1: the core at the full population.
    const AbAmvaResult<T> step1 =
        detail::ab_core(N, nservers, sched, v, S, maxiter, Delta, L, fcfsSchmidt, method);
    const Matrix<T> LUpdated = step1.QN;

    // STEP 2: the core at N - e_r, one class at a time.
    for (std::size_t r = 0; r < K; ++r) {
        std::vector<int> popWithout = N;
        popWithout[r] -= 1;
        Matrix<T> lWithoutC(M, K, zero);
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t c = 0; c < K; ++c) lWithoutC(j, c) = lWithoutR[c](j, c);
        const AbAmvaResult<T> ret = detail::ab_core(popWithout, nservers, sched, v, S, maxiter,
                                                    Delta, lWithoutC, fcfsSchmidt, method);
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t c = 0; c < K; ++c) lWithoutR[c](j, r) = ret.QN(j, c);
    }

    // STEP 3: the fractional-change tensor.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            const T F_ir = N[r] != 0 ? T(LUpdated(i, r) / num_traits<T>::from_int(N[r])) : zero;
            for (std::size_t t = 0; t < K; ++t) {
                const long divisor = r == t ? N[r] - 1 : N[r];
                const T F_irt =
                    divisor != 0 ? T(lWithoutR[r](i, t) / num_traits<T>::from_int(divisor)) : zero;
                Delta[(i * K + r) * K + t] = F_irt - F_ir;
            }
        }

    // STEP 4: the core again at N, with the step-1 queue lengths and step-3
    // fractional changes.
    return detail::ab_core(N, nservers, sched, v, S, maxiter, Delta, LUpdated, fcfsSchmidt, method);
}

/** Reference defaults: no Schmidt FCFS wait, the AB marginal rule. */
template <class T>
AbAmvaResult<T> pfqn_ab_amva(const Matrix<T>& S, const std::vector<int>& N, const Matrix<T>& v,
                             const std::vector<int>& nservers,
                             const std::vector<SchedStrategy>& sched) {
    return pfqn_ab_amva(S, N, v, nservers, sched, false, AbMarginalMethod::Ab);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_AB_AMVA_H
