/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SCHMIDT_H
#define LINE_API_PFQN_SCHMIDT_H

/**
 * Schmidt's MVA for closed networks with general scheduling disciplines and
 * class-dependent multiserver FCFS stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_schmidt.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/ld/Pfqn_schmidt.java.
 *
 * The recursion walks the population lattice 0 <= k <= N and, at each k,
 * computes the residence time by the arrival theorem. Three station kinds are
 * distinguished, as in the reference:
 *
 *  - INF: w = D(i,c).
 *  - PS, and FCFS with class-independent demands: the standard
 *    w = D(i,c)/s (1 + sum_r L(i,r | k - e_c)) plus, for s > 1, the
 *    idle-server correction sum_{j=1}^{s-1} (s-j) Pr(j-1 busy | k - e_c) D/s.
 *    Both need the scalar busy-server distribution Pr(j | k).
 *  - FCFS with class-dependent demands and s > 1: the full per-class state
 *    distribution Pr(nvec | k) is carried, and
 *    w = sum_{nvec <= k, nvec_c > 0} B_c(nvec) Pr(nvec - e_c | k - e_c),
 *    with B_c the queue-composition-weighted mean service time.
 *
 * Pure service times. The B_c terms need per-visit service times S = D/v, not
 * demands, so the visit ratios enter explicitly; with v == 1 the two coincide.
 *
 * Arithmetic: EXACT-CAPABLE, no transcendental gate. Every step is a finite
 * sum, product or quotient over the population lattice, so the whole recursion
 * stays in the field of the inputs and is exact in rational arithmetic. The
 * reference's three floating-point guards (max(v, 1e-12) on the visit ratios,
 * max(s (sum(nvec) - 1), 1e-12) on the B_c denominator, and the max(eps, .)
 * and max(1e-12, .) floors on the idle-state probability) are carried over
 * verbatim as constants of the algorithm, since removing them would change the
 * numbers the reference produces; they are the only place a magic constant
 * enters, and none of them makes the arithmetic inexact.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_schmidt, mirroring [XN,QN,UN,CN]. */
template <class T>
struct SchmidtResult {
    std::vector<T> XN;  ///< (R) per-class throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization, D X / s
    Matrix<T> CN;       ///< (M x R) residence time
};

namespace detail {

/** Which marginal distribution a station needs, if any. */
enum class SchmidtPc { None, Scalar, Vector };

}  // namespace detail

/**
 * @param D     (M x R) service demands
 * @param N     (R) population per class
 * @param S     (M x R) or (M x 1) server counts
 * @param sched (M) scheduling discipline per station
 * @param v     (M x R) visit ratios; empty for all ones
 */
template <class T>
SchmidtResult<T> pfqn_schmidt(const Matrix<T>& D, const std::vector<int>& N, const Matrix<int>& S,
                              const std::vector<SchedStrategy>& sched, const Matrix<T>& v) {
    const std::size_t M = D.rows();
    const std::size_t R = N.size();
    if (!D.empty() && D.cols() != R)
        throw InputError("pfqn_schmidt: demand matrix and population vector disagree on the class count");
    if (sched.size() != M)
        throw InputError("pfqn_schmidt: scheduling vector has the wrong station count");
    if (S.rows() != M || (S.cols() != R && S.cols() != 1))
        throw InputError("pfqn_schmidt: server-count matrix has the wrong shape");
    if (!v.empty() && (v.rows() != M || v.cols() != R))
        throw InputError("pfqn_schmidt: visit-ratio matrix has the wrong shape");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_schmidt: negative population");
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < S.cols(); ++c)
            if (S(i, c) < 1) throw InputError("pfqn_schmidt: server count below one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T tiny = num_traits<T>::from_double(1e-12);
    const T epsT = num_traits<T>::from_double(2.220446049250313e-16);

    SchmidtResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    if (M == 0) return res;

    Matrix<T> vis(M, R, one);
    if (!v.empty()) vis = v;
    // Pure per-visit service times, S_pure = D / max(v, 1e-12).
    Matrix<T> Sp(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < R; ++c) {
            const T den = vis(i, c) > tiny ? vis(i, c) : tiny;
            Sp(i, c) = D(i, c) / den;
        }

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);
    long Ntot = 0;
    for (int n : N) Ntot += n;

    auto nserv = [&](std::size_t i, std::size_t c) {
        return S.cols() == 1 ? S(i, 0) : S(i, static_cast<std::size_t>(c));
    };
    // A station has class-independent demands when every class shares D(i,0).
    std::vector<bool> classIndep(M, true);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 1; c < R; ++c)
            if (D(i, c) != D(i, 0)) classIndep[i] = false;

    std::vector<detail::SchmidtPc> kind(M, detail::SchmidtPc::None);
    for (std::size_t i = 0; i < M; ++i) {
        bool single = true;
        for (std::size_t c = 0; c < (S.cols() == 1 ? std::size_t(1) : R); ++c)
            if (nserv(i, c) != 1) single = false;
        switch (sched[i]) {
            case SchedStrategy::INF:
                break;
            case SchedStrategy::PS:
                if (!single) kind[i] = detail::SchmidtPc::Scalar;
                break;
            case SchedStrategy::FCFS:
                if (classIndep[i]) {
                    if (!single) kind[i] = detail::SchmidtPc::Scalar;
                } else {
                    kind[i] = detail::SchmidtPc::Vector;
                }
                break;
        }
    }

    std::vector<Matrix<T>> Lq(M, Matrix<T>(R, total, zero));
    std::vector<Matrix<T>> Pc(M);
    for (std::size_t i = 0; i < M; ++i) {
        if (kind[i] == detail::SchmidtPc::Scalar)
            Pc[i] = Matrix<T>(static_cast<std::size_t>(1 + Ntot), total, zero);
        else if (kind[i] == detail::SchmidtPc::Vector)
            Pc[i] = Matrix<T>(total, total, zero);
        if (kind[i] != detail::SchmidtPc::None) Pc[i](0, 0) = one;  // Pr(0 | 0) = 1
    }

    // x[(i*R + c)*total + h] and w likewise.
    std::vector<T> x(M * R * total, zero), w(M * R * total, zero);

    std::vector<int> kvec(R, 0);
    std::size_t hlast = 0;
    bool more = true;
    while (more) {
        const std::size_t hk = pop_index(kvec, prods);
        hlast = hk;
        long kpop = 0;
        for (int t : kvec) kpop += t;

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < R; ++c) {
                if (kvec[c] <= 0) continue;
                const std::size_t hkc = hk - prods[c];
                const int ns = nserv(i, c);
                T& wi = w[(i * R + c) * total + hk];
                if (sched[i] == SchedStrategy::INF) {
                    wi = D(i, c);
                    continue;
                }
                const bool vectorPc = kind[i] == detail::SchmidtPc::Vector;
                if (!vectorPc || ns == 1) {
                    T qtot = zero;
                    for (std::size_t r = 0; r < R; ++r) qtot += Lq[i](r, hkc);
                    if (ns == 1) {
                        wi = D(i, c) * (one + qtot);
                    } else {
                        const T nsT = num_traits<T>::from_int(ns);
                        wi = D(i, c) / nsT * (one + qtot);
                        for (int j = 1; j <= ns - 1; ++j)
                            wi += num_traits<T>::from_int(ns - j) *
                                  Pc[i](static_cast<std::size_t>(j - 1), hkc) * (D(i, c) / nsT);
                    }
                } else {
                    // Class-dependent multiserver FCFS: sum over the joint state.
                    const T nsT = num_traits<T>::from_int(ns);
                    std::vector<int> nvec(R, 0);
                    bool more_n = true;
                    while (more_n) {
                        if (nvec[c] > 0) {
                            long nsum = 0;
                            for (int t : nvec) nsum += t;
                            const std::size_t hnc = pop_index(nvec, prods) - prods[c];
                            T Bcn = Sp(i, c);
                            if (nsum > ns) {
                                T sumVal = zero;
                                for (std::size_t r = 0; r < R; ++r)
                                    sumVal += num_traits<T>::from_int(nvec[r]) * Sp(i, r);
                                const T den0 = nsT * num_traits<T>::from_int(nsum - 1);
                                const T den = den0 > tiny ? den0 : tiny;
                                Bcn += num_traits<T>::from_int(nsum - ns) / den *
                                       (sumVal - Sp(i, c));
                            }
                            wi += Bcn * Pc[i](hnc, hkc);
                        }
                        more_n = next_pop(nvec, kvec);
                    }
                }
            }

        for (std::size_t c = 0; c < R; ++c) {
            T denom = zero;
            for (std::size_t i = 0; i < M; ++i) denom += vis(i, c) * w[(i * R + c) * total + hk];
            for (std::size_t i = 0; i < M; ++i)
                x[(i * R + c) * total + hk] =
                    denom > zero ? vis(i, c) * num_traits<T>::from_int(kvec[c]) / denom : zero;
        }

        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t c = 0; c < R; ++c)
                Lq[i](c, hk) = x[(i * R + c) * total + hk] * w[(i * R + c) * total + hk];

            if (kind[i] == detail::SchmidtPc::Scalar) {
                // Pr(j busy | k) from Pr(j-1 busy | k - e_c).
                const int s0 = nserv(i, 0);
                const long jmax = sched[i] == SchedStrategy::PS
                                      ? (s0 < kpop ? s0 : kpop)
                                      : (s0 < kpop ? s0 : kpop) - 1;
                for (long n = 1; n <= jmax; ++n)
                    for (std::size_t c = 0; c < R; ++c) {
                        if (kvec[c] <= 0) continue;
                        const std::size_t hkc = hk - prods[c];
                        Pc[i](static_cast<std::size_t>(n), hk) +=
                            D(i, c) / num_traits<T>::from_int(n) * x[(i * R + c) * total + hk] *
                            Pc[i](static_cast<std::size_t>(n - 1), hkc);
                    }
                if (jmax >= 1) {
                    const long top = sched[i] == SchedStrategy::PS ? (s0 < kpop ? s0 : kpop)
                                                                   : (s0 < kpop ? s0 : kpop);
                    T acc = zero;
                    for (long n = 1; n <= top; ++n) acc += Pc[i](static_cast<std::size_t>(n), hk);
                    const T p0 = one - acc;
                    Pc[i](0, hk) = p0 > epsT ? p0 : epsT;
                }
            } else if (kind[i] == detail::SchmidtPc::Vector) {
                const int ns = nserv(i, 0);
                const T nsT = num_traits<T>::from_int(ns);
                T sumAll = zero;
                std::vector<int> nvec(R, 0);
                bool more_n = next_pop(nvec, kvec);  // skip the zero vector
                while (more_n) {
                    const std::size_t hn = pop_index(nvec, prods);
                    long nsum = 0;
                    for (int t : nvec) nsum += t;
                    T prob = zero;
                    for (std::size_t c = 0; c < R; ++c) {
                        if (nvec[c] <= 0 || kvec[c] <= 0) continue;
                        const std::size_t hnc = hn - prods[c];
                        const std::size_t hkc = hk - prods[c];
                        T Bcn = Sp(i, c);
                        if (nsum > 1) {
                            T sumVal = zero;
                            for (std::size_t r = 0; r < R; ++r)
                                sumVal += num_traits<T>::from_int(nvec[r]) * Sp(i, r);
                            const T den0 = nsT * num_traits<T>::from_int(nsum - 1);
                            const T den = den0 > tiny ? den0 : tiny;
                            Bcn += num_traits<T>::from_int(nsum - ns > 0 ? nsum - ns : 0) / den *
                                   (sumVal - Sp(i, c));
                        }
                        prob += Bcn / num_traits<T>::from_int(nsum) *
                                x[(i * R + c) * total + hk] * Pc[i](hnc, hkc);
                    }
                    Pc[i](hn, hk) = prob;
                    sumAll += prob;
                    more_n = next_pop(nvec, kvec);
                }
                const T p0 = one - sumAll;
                Pc[i](0, hk) = p0 > tiny ? p0 : tiny;
            }
        }

        more = next_pop(kvec, N);
    }

    for (std::size_t c = 0; c < R; ++c) {
        T tot = zero;
        for (std::size_t i = 0; i < M; ++i) tot += w[(i * R + c) * total + hlast];
        res.XN[c] = tot > zero ? num_traits<T>::from_int(N[c]) / tot : zero;
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < R; ++c) {
            res.UN(i, c) = D(i, c) * res.XN[c] / num_traits<T>::from_int(nserv(i, c));
            res.CN(i, c) = w[(i * R + c) * total + hlast];
            res.QN(i, c) = Lq[i](c, hlast);
        }
    return res;
}

/** Unit visit ratios, the MATLAB default. */
template <class T>
SchmidtResult<T> pfqn_schmidt(const Matrix<T>& D, const std::vector<int>& N, const Matrix<int>& S,
                              const std::vector<SchedStrategy>& sched) {
    return pfqn_schmidt(D, N, S, sched, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SCHMIDT_H
