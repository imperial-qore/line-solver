/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVAMS_H
#define LINE_API_PFQN_MVAMS_H

/**
 * Exact Mean Value Analysis for mixed open/closed networks with multiserver
 * stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mvams.m. That function is a
 * dispatcher over four exact MVA variants, so the port carries them all:
 *
 *   no multiserver, closed   -> pfqn_mva      (pfqn_mva.h, already ported)
 *   no multiserver, mixed    -> pfqn_mvamx    (below)
 *   multiserver,    closed   -> pfqn_mvald    (below)
 *   multiserver,    mixed    -> pfqn_mvaldms  (below, via pfqn_mvaldmx)
 *
 * pfqn_mvald is the load-dependent MVA of Reiser and Lavenberg, carrying the
 * marginal queue-length distribution pi(k|n) along the population lattice;
 * pfqn_mvaldmx is the Bruell-Balbo-Ashfari mixed extension for limited load
 * dependence, whose effective-capacity terms come from pfqn_ldmx_ec. A
 * multiserver station with S servers is the load-dependent station with rates
 * mu(i,k) = min(k, S), which is where the "ms" in the name comes from.
 *
 * Arithmetic. Every step of all four variants is an addition, a subtraction, a
 * multiplication, a division or an integer power in the field of the inputs:
 * the family is exact-capable end to end and needs no transcendental function.
 * The only non-field step in MATLAB is the log used to accumulate lG, and the
 * port does what pfqn_mva.h already does, accumulating the product of the
 * reciprocal throughputs along the lattice path and taking the log once at the
 * end of a value that is still exact.
 *
 * Contract notes, where this port deliberately differs from MATLAB:
 *
 *  1. UN in the closed multiserver branch. MATLAB documents pfqn_mvams as
 *     returning an (M x R) utilization, and three of its four branches do, but
 *     the closed multiserver branch forwards pfqn_mvald's UN, which is the
 *     (M x 1) aggregate 1 - P(station i empty). This port returns the (M x R)
 *     per-class utilization law U(i,r) = X(r) L(i,r) / S(i), the same formula
 *     pfqn_mvaldms already uses in the mixed multiserver branch, so the metric
 *     set matches MvaResult and is consistent across all four branches. The
 *     aggregate form is still available from pfqn_mvald directly, whose result
 *     carries both UN and the full marginal distribution.
 *
 *  2. Queue replicas with multiservers. MATLAB rejects mi > 1 for the mixed
 *     multiserver branch but silently drops mi in the closed multiserver
 *     branch, because pfqn_mvald has no mi argument: the replicas then affect
 *     only the residence time reported for absent classes. This port rejects
 *     mi > 1 in both multiserver branches rather than returning numbers that
 *     ignore an argument the caller supplied.
 *
 * Index conventions. An open class is marked by OPEN_CLASS in the population
 * vector (MATLAB uses Inf, which an int vector cannot hold) and an infinite
 * server station by INF_SERVERS in the server-count vector (MATLAB uses Inf).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** Marks an open (infinite-population) class in a population vector. */
constexpr int OPEN_CLASS = -1;

/** Marks an infinite-server station in a server-count vector. */
constexpr int INF_SERVERS = -1;

/** True when the population entry denotes an open class. */
inline bool is_open_class(int n) { return n < 0; }

// ---------------------------------------------------------------------------
// pfqn_mvald: exact MVA for closed load-dependent networks
// ---------------------------------------------------------------------------

/**
 * Result of pfqn_mvald, mirroring the seven MATLAB outputs. It is deliberately
 * not MvaResult: the load-dependent family reports an aggregate utilization
 * and a per-class cycle time rather than per-station, per-class matrices, and
 * it additionally exposes the marginal queue-length distribution.
 */
template <class T>
struct MvaLdResult {
    std::vector<T> XN;   ///< (R) per-class throughput
    Matrix<T> QN;        ///< (M x R) mean queue length
    std::vector<T> UN;   ///< (M) utilization, 1 - P(station empty)
    std::vector<T> CN;   ///< (R) cycle time, exclusive of think time
    Matrix<T> WN;        ///< (M x R) residence time at the full population
    Matrix<T> PI;        ///< (M x (Nt+1)) marginal queue-length distribution at N
    T G;                 ///< normalizing constant
    double lG;           ///< log of the normalizing constant
    bool isNumStable;    ///< false once a marginal probability had to be clamped
};

/**
 * Exact MVA for a closed network of load-dependent stations.
 *
 * Port of matlab/src/api/pfqn/pfqn_mvald.m. The recursion over the population
 * lattice is
 *
 *   W(i,r|n)  = sum_{k=1}^{|n|} L(i,r)/mu(i,k) * k * pi(i,k-1|n - e_r)
 *   X(r|n)    = n_r / (Z_r + sum_i W(i,r|n))
 *   pi(i,k|n) = sum_r L(i,r)/mu(i,k) * X(r|n) * pi(i,k-1|n - e_r)
 *   pi(i,0|n) = 1 - sum_{k>=1} pi(i,k|n)
 *
 * @param L         (M x R) service demands
 * @param N         (R) population per class, all finite and non-negative
 * @param Z         (K x R) think times, summed over rows; may be empty
 * @param mu        (M x Nt') service rates, Nt' >= sum(N); mu(i,k-1) is the
 *                  rate of station i while it holds k jobs
 * @param stabilize when true (the MATLAB default) a marginal probability that
 *                  comes out negative is clamped to the double epsilon rather
 *                  than propagated. The clamp is a floating-point guard: in
 *                  exact arithmetic pi(i,0|n) of a well-posed product-form
 *                  model is non-negative and the branch is never taken. The
 *                  clamp value is a dyadic rational, so it is representable
 *                  without rounding in every supported arithmetic.
 */
template <class T>
MvaLdResult<T> pfqn_mvald(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                          const Matrix<T>& mu, bool stabilize = true) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mvald: demand matrix and population vector disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    MvaLdResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN.assign(M, zero);
    res.CN.assign(R, zero);
    res.WN = Matrix<T>(M, R, zero);
    res.G = one;
    res.lG = 0.0;
    res.isNumStable = true;

    long Nt = 0;
    for (int v : N) {
        // MATLAB returns all-zero metrics and lGN = -Inf for a negative population.
        if (v < 0) {
            res.PI = Matrix<T>(M, 1, zero);
            res.G = zero;
            res.lG = -std::numeric_limits<double>::infinity();
            return res;
        }
        Nt += v;
    }

    if (mu.rows() != M)
        throw InputError("pfqn_mvald: rate matrix and demand matrix disagree on the station count");
    if (static_cast<long>(mu.cols()) < Nt)
        throw InputError("pfqn_mvald: rate matrix needs one column per job in the total population");

    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_mvald: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    const std::size_t K = static_cast<std::size_t>(Nt) + 1;  // queue lengths 0 ... Nt
    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // pi[(idx * M + i) * K + k] is P(station i holds k jobs | population idx).
    std::vector<T> pi(total * M * K, zero);
    const T tiny = num_traits<T>::from_double(std::numeric_limits<double>::epsilon());

    std::vector<int> n(R, 0);
    std::vector<T> x(R, zero);
    bool more = true;
    while (more) {
        const std::size_t idx = pop_index(n, prods);
        long nsum = 0;
        for (int v : n) nsum += v;

        res.WN.fill(zero);
        std::fill(x.begin(), x.end(), zero);

        for (std::size_t s = 0; s < R; ++s) {
            if (n[s] == 0) continue;
            const std::size_t idx_s = idx - prods[s];
            T ctot = Zsum[s];
            for (std::size_t i = 0; i < M; ++i) {
                T w = zero;
                for (long k = 1; k <= nsum; ++k)
                    w += (L(i, s) / mu(i, static_cast<std::size_t>(k - 1))) *
                         num_traits<T>::from_int(k) *
                         pi[(idx_s * M + i) * K + static_cast<std::size_t>(k - 1)];
                res.WN(i, s) = w;
                ctot += w;
            }
            if (ctot == zero) throw NumericError("pfqn_mvald: zero total residence time");
            x[s] = num_traits<T>::from_int(n[s]) / ctot;
        }

        for (long k = 1; k <= nsum; ++k) {
            for (std::size_t i = 0; i < M; ++i) {
                T acc = zero;
                for (std::size_t s = 0; s < R; ++s) {
                    if (n[s] == 0) continue;
                    const std::size_t idx_s = idx - prods[s];
                    acc += (L(i, s) / mu(i, static_cast<std::size_t>(k - 1))) * x[s] *
                           pi[(idx_s * M + i) * K + static_cast<std::size_t>(k - 1)];
                }
                pi[(idx * M + i) * K + static_cast<std::size_t>(k)] = acc;
            }
        }

        for (std::size_t i = 0; i < M; ++i) {
            T acc = zero;
            for (long k = 1; k <= nsum; ++k) acc += pi[(idx * M + i) * K + static_cast<std::size_t>(k)];
            T p0 = one - acc;
            if (p0 < zero) {
                res.isNumStable = false;
                if (stabilize) p0 = tiny;
            }
            pi[(idx * M + i) * K + 0] = p0;
        }

        // Normalizing constant along the lattice path that fills class 0, then
        // class 1, and so on, exactly as in pfqn_mva.
        long last_nnz = -1;
        for (long r = static_cast<long>(R) - 1; r >= 0; --r)
            if (n[r] != 0) {
                last_nnz = r;
                break;
            }
        if (last_nnz >= 0) {
            bool prefixFull = true;
            for (long r = 0; r < last_nnz; ++r)
                if (n[r] != N[r]) {
                    prefixFull = false;
                    break;
                }
            bool suffixEmpty = true;
            for (std::size_t r = static_cast<std::size_t>(last_nnz) + 1; r < R; ++r)
                if (n[r] != 0) {
                    suffixEmpty = false;
                    break;
                }
            if (prefixFull && suffixEmpty) {
                const T& xr = x[static_cast<std::size_t>(last_nnz)];
                if (xr == zero) throw NumericError("pfqn_mvald: zero throughput on the G path");
                res.G /= xr;
            }
        }

        more = next_pop(n, N);
    }

    // The loop leaves x, WN and the pi slice at the full population N.
    const std::size_t idxN = total - 1;
    res.XN = x;
    res.PI = Matrix<T>(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) res.PI(i, k) = pi[(idxN * M + i) * K + k];

    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) res.QN(i, r) = res.WN(i, r) * res.XN[r];
        res.UN[i] = one - res.PI(i, 0);
    }
    for (std::size_t r = 0; r < R; ++r) {
        // An absent class has zero throughput, so N/X would be 0/0: MATLAB
        // reports no cycle time for it rather than NaN.
        if (N[r] > 0) res.CN[r] = num_traits<T>::from_int(N[r]) / res.XN[r] - Zsum[r];
    }

    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

// ---------------------------------------------------------------------------
// pfqn_mvamx: exact MVA for mixed single-server networks
// ---------------------------------------------------------------------------

/**
 * Exact MVA for a mixed open/closed network of single-server stations.
 *
 * Port of matlab/src/api/pfqn/pfqn_mvamx.m. The open classes are absorbed by
 * inflating the closed demands, D_c(i,r) / (1 - sum_{open} lambda D), after
 * which the closed subnetwork is solved by pfqn_mva; the open metrics then
 * follow from the closed queue lengths.
 *
 * @param lambda (R) arrival rates, zero on closed classes
 * @param D      (M x R) service demands
 * @param N      (R) population, OPEN_CLASS on open classes
 * @param Z      (K x R) think times, summed over rows; may be empty
 * @param mi     (M) station multiplicities; empty for all ones
 *
 * G and lG describe the closed subnetwork on the inflated demands, and are set
 * to 0 and NaN respectively when there is no closed class, as in MATLAB.
 */
template <class T>
MvaResult<T> pfqn_mvamx(const std::vector<T>& lambda, const Matrix<T>& D,
                        const std::vector<int>& N, const Matrix<T>& Z,
                        const std::vector<int>& mi) {
    const std::size_t M = D.rows();
    const std::size_t R = N.size();
    if (!D.empty() && D.cols() != R)
        throw InputError("pfqn_mvamx: demand matrix and population vector disagree on the class count");
    if (lambda.size() != R)
        throw InputError("pfqn_mvamx: arrival-rate vector and population vector disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    MvaResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.G = zero;
    res.lG = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::size_t> openClasses, closedClasses;
    for (std::size_t r = 0; r < R; ++r) {
        if (is_open_class(N[r]))
            openClasses.push_back(r);
        else
            closedClasses.push_back(r);
    }
    for (std::size_t r = 0; r < R; ++r)
        if (lambda[r] > zero && !is_open_class(N[r]) && N[r] > 0)
            throw InputError("pfqn_mvamx: arrival rate cannot be specified on a closed class");

    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_mvamx: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    for (std::size_t r : openClasses) {
        for (std::size_t i = 0; i < M; ++i) res.UN(i, r) = lambda[r] * D(i, r);
        res.XN[r] = lambda[r];
    }

    std::vector<T> UNt(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) UNt[i] += res.UN(i, r);

    const std::size_t nClosed = closedClasses.size();
    Matrix<T> QNc;
    if (nClosed > 0) {
        Matrix<T> Dc(M, nClosed, zero);
        for (std::size_t i = 0; i < M; ++i) {
            const T slack = one - UNt[i];
            if (!(slack > zero))
                throw NumericError("pfqn_mvamx: open classes saturate a station, no closed solution exists");
            for (std::size_t c = 0; c < nClosed; ++c) Dc(i, c) = D(i, closedClasses[c]) / slack;
        }
        std::vector<int> Nc(nClosed);
        Matrix<T> Zc(1, nClosed, zero);
        for (std::size_t c = 0; c < nClosed; ++c) {
            Nc[c] = N[closedClasses[c]];
            Zc(0, c) = Zsum[closedClasses[c]];
        }
        MvaResult<T> closed = pfqn_mva(Dc, Nc, Zc, mi);
        QNc = closed.QN;
        for (std::size_t c = 0; c < nClosed; ++c) {
            const std::size_t r = closedClasses[c];
            res.XN[r] = closed.XN[c];
            for (std::size_t i = 0; i < M; ++i) {
                res.QN(i, r) = closed.QN(i, c);
                res.CN(i, r) = closed.CN(i, c);
            }
        }
        res.G = closed.G;
        res.lG = closed.lG;
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < nClosed; ++c) {
            const std::size_t r = closedClasses[c];
            res.UN(i, r) = res.XN[r] * D(i, r);
        }

    for (std::size_t i = 0; i < M; ++i) {
        T qc = zero;
        for (std::size_t c = 0; c < nClosed; ++c) qc += QNc(i, c);
        const T slack = one - UNt[i];
        for (std::size_t r : openClasses) {
            if (!(slack > zero))
                throw NumericError("pfqn_mvamx: open classes saturate a station");
            res.CN(i, r) = D(i, r) * (one + qc) / slack;
            res.QN(i, r) = res.CN(i, r) * res.XN[r];
        }
    }
    return res;
}

// ---------------------------------------------------------------------------
// pfqn_ldmx_ec: effective-capacity terms of the mixed load-dependent MVA
// ---------------------------------------------------------------------------

namespace detail {

/** Outputs of pfqn_ldmx_ec, with E and Eprime indexed by queue length 0 ... Nt. */
template <class T>
struct LdmxEc {
    Matrix<T> EC;      ///< (M x Nt) effective capacity, EC(i,k-1) for k jobs
    Matrix<T> E;       ///< (M x (Nt+1))
    Matrix<T> Eprime;  ///< (M x (Nt+1))
    std::vector<T> Lo; ///< (M) open-class load at each station
};

/**
 * Port of matlab/src/api/pfqn/pfqn_ldmx_ec.m. Computes the terms that let the
 * mixed load-dependent MVA absorb the open classes into a station-dependent
 * effective capacity, under the limited-load-dependence assumption that the
 * rate saturates from the b(i)-th job onwards.
 *
 * Every step is a field operation or an integer power, so the routine is exact.
 *
 * @param lambda (R) arrival rates
 * @param D      (M x R) service demands
 * @param mu     (M x Nt) service rates
 */
template <class T>
LdmxEc<T> pfqn_ldmx_ec(const std::vector<T>& lambda, const Matrix<T>& D, const Matrix<T>& mu) {
    const std::size_t M = mu.rows();
    const std::size_t Nt = mu.cols();
    const std::size_t R = D.cols();
    if (D.rows() != M)
        throw InputError("pfqn_ldmx_ec: demand matrix and rate matrix disagree on the station count");
    if (lambda.size() != R)
        throw InputError("pfqn_ldmx_ec: arrival-rate vector and demand matrix disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    LdmxEc<T> out;
    out.Lo.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) out.Lo[i] += lambda[r] * D(i, r);

    // b(i): first queue length at which the rate has reached its final value.
    std::vector<std::size_t> b(M, 1);
    std::size_t maxb = 1;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 1; j <= Nt; ++j) {
            if (mu(i, j - 1) == mu(i, Nt - 1)) {
                b[i] = j;
                break;
            }
        }
        if (b[i] > maxb) maxb = b[i];
    }

    // C(i,j) = 1/mu(i,j), extended past Nt by repeating the saturated rate, so
    // the F-recursions below can index up to Nt + 1 + maxb as MATLAB does.
    const std::size_t Cn = Nt + 1 + maxb;
    Matrix<T> C(M, Cn, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 1; j <= Cn; ++j) {
            const T& m = j <= Nt ? mu(i, j - 1) : mu(i, Nt - 1);
            if (m == zero) throw InputError("pfqn_ldmx_ec: service rate must be nonzero");
            C(i, j - 1) = one / m;
        }

    out.EC = Matrix<T>(M, Nt, zero);
    out.E = Matrix<T>(M, Nt + 1, zero);
    out.Eprime = Matrix<T>(M, Nt + 1, zero);

    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t bi = b[i];
        const T Cb = C(i, bi - 1);
        const T slack = one - out.Lo[i] * Cb;
        if (slack == zero)
            throw NumericError("pfqn_ldmx_ec: open-class load saturates a station, effective capacity diverges");
        const T invSlack = one / slack;

        // prod_{j=1}^{b-1} C(j)/C(b), shared by E1(0) and F3(0,0).
        T ratioProd = one;
        for (std::size_t j = 1; j + 1 <= bi; ++j) ratioProd *= C(i, j - 1) / Cb;

        const long nb = static_cast<long>(bi) - 2;  // last F-column index, -1 when b == 1
        std::vector<T> E1(Nt + 1, zero);
        T F3prev0 = zero;

        for (std::size_t n = 0; n <= Nt; ++n) {
            if (n >= bi) {
                out.E(i, n) = num_pow_int(invSlack, static_cast<unsigned>(n + 1));
                out.Eprime(i, n) = Cb * out.E(i, n);
                continue;
            }

            if (n == 0)
                E1[0] = invSlack * ratioProd;
            else
                E1[n] = invSlack * Cb / C(i, n - 1) * E1[n - 1];

            T E2 = zero, E3 = zero, E2prime = zero;

            T F2 = zero, F3 = zero, F2p = zero;
            for (long n0 = 0; n0 <= nb; ++n0) {
                if (n0 == 0) {
                    F2 = one;
                    F3 = (n == 0) ? ratioProd : Cb / C(i, n - 1) * F3prev0;
                    F2p = C(i, n);  // C(n+1) in 1-based terms
                } else {
                    const T fac = num_traits<T>::from_int(static_cast<long>(n) + n0) /
                                  num_traits<T>::from_int(n0);
                    F2 = fac * out.Lo[i] * C(i, static_cast<std::size_t>(static_cast<long>(n) + n0) - 1) * F2;
                    F3 = fac * out.Lo[i] * Cb * F3;
                    F2p = fac * out.Lo[i] * C(i, static_cast<std::size_t>(static_cast<long>(n) + n0)) * F2p;
                }
                if (n0 == 0) F3prev0 = F3;
                E2 += F2;
                E3 += F3;
                E2prime += F2p;
            }
            if (nb < 0) {
                // b == 1: the F-sums are empty and F3(n,0) is never formed, so
                // the carry-over stays at its unused zero, as in MATLAB.
                F3prev0 = zero;
            }

            out.E(i, n) = E1[n] + E2 - E3;
            if (n + 1 < bi)
                out.Eprime(i, n) = Cb * E1[n] + E2prime - Cb * E3;
            else
                out.Eprime(i, n) = Cb * out.E(i, n);
        }

        for (std::size_t n = 1; n <= Nt; ++n) {
            if (out.E(i, n - 1) == zero)
                throw NumericError("pfqn_ldmx_ec: vanishing effective capacity");
            out.EC(i, n - 1) = C(i, n - 1) * out.E(i, n) / out.E(i, n - 1);
        }
    }
    return out;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// pfqn_mvaldmx / pfqn_mvaldms: mixed networks with load-dependent stations
// ---------------------------------------------------------------------------

/**
 * Exact MVA for mixed open/closed networks with limited load dependence.
 *
 * Port of matlab/src/api/pfqn/pfqn_mvaldmx.m (Bruell, Balbo and Ashfari). The
 * MATLAB signature carries a trailing server-count argument S that its body
 * never reads; the port drops it, since pfqn_mvaldms is the caller that turns
 * server counts into rates.
 *
 * @param lambda (R) arrival rates, zero on closed classes
 * @param D      (M x R) service demands
 * @param N      (R) population, OPEN_CLASS on open classes
 * @param Z      (K x R) think times, summed over rows; may be empty
 * @param mu     (M x Nc') rates, Nc' >= the total closed population
 */
template <class T>
MvaResult<T> pfqn_mvaldmx(const std::vector<T>& lambda, const Matrix<T>& D,
                          const std::vector<int>& N, const Matrix<T>& Z, const Matrix<T>& mu) {
    const std::size_t M = D.rows();
    const std::size_t R = N.size();
    if (!D.empty() && D.cols() != R)
        throw InputError("pfqn_mvaldmx: demand matrix and population vector disagree on the class count");
    if (lambda.size() != R)
        throw InputError("pfqn_mvaldmx: arrival-rate vector and population vector disagree on the class count");
    if (mu.rows() != M)
        throw InputError("pfqn_mvaldmx: rate matrix and demand matrix disagree on the station count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> openClasses, closedClasses;
    long Nct = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (is_open_class(N[r])) {
            openClasses.push_back(r);
        } else {
            closedClasses.push_back(r);
            Nct += N[r];
        }
    }
    for (std::size_t r = 0; r < R; ++r)
        if (lambda[r] > zero && !is_open_class(N[r]) && N[r] > 0)
            throw InputError("pfqn_mvaldmx: arrival rate cannot be specified on a closed class");
    if (static_cast<long>(mu.cols()) < Nct)
        throw InputError(
            "pfqn_mvaldmx: the load-dependent rates must be given for at least the maximum closed population");

    // MATLAB appends one extra column so the recursion can look one job ahead.
    Matrix<T> mux(M, mu.cols() + 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < mu.cols(); ++j) mux(i, j) = mu(i, j);
        mux(i, mu.cols()) = mu(i, mu.cols() - 1);
    }
    const detail::LdmxEc<T> ec = detail::pfqn_ldmx_ec(lambda, D, mux);

    MvaResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.G = one;
    res.lG = 0.0;

    const std::size_t Cc = closedClasses.size();
    Matrix<T> Dc(M, Cc, zero);
    std::vector<int> Nc(Cc, 0);
    std::vector<T> Zc(Cc, zero);
    for (std::size_t c = 0; c < Cc; ++c) {
        const std::size_t r = closedClasses[c];
        Nc[c] = N[r];
        for (std::size_t i = 0; i < M; ++i) Dc(i, c) = D(i, r);
    }
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_mvaldmx: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t c = 0; c < Cc; ++c) Zc[c] += Z(k, closedClasses[c]);
    }

    const std::vector<std::size_t> prods = plane_sizes(Nc);
    const std::size_t total = population_count(Nc);
    const std::size_t K = static_cast<std::size_t>(Nct) + 1;

    // Pc[(idx * M + i) * K + k] is P(station i holds k closed jobs | idx).
    std::vector<T> Pc(total * M * K, zero);
    for (std::size_t i = 0; i < M; ++i) Pc[(0 * M + i) * K + 0] = one;
    const T tiny = num_traits<T>::from_double(std::numeric_limits<double>::epsilon());

    Matrix<T> w(M, Cc, zero);
    std::vector<T> x(Cc, zero);
    std::vector<int> nvec(Cc, 0);
    bool more = true;
    while (more) {
        const std::size_t idx = pop_index(nvec, prods);
        long nc = 0;
        for (int v : nvec) nc += v;

        w.fill(zero);
        std::fill(x.begin(), x.end(), zero);

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < Cc; ++c) {
                if (nvec[c] == 0) continue;
                const std::size_t idx_c = idx - prods[c];
                T acc = zero;
                for (long k = 1; k <= nc; ++k)
                    acc += Dc(i, c) * num_traits<T>::from_int(k) *
                           ec.EC(i, static_cast<std::size_t>(k - 1)) *
                           Pc[(idx_c * M + i) * K + static_cast<std::size_t>(k - 1)];
                w(i, c) = acc;
            }

        for (std::size_t c = 0; c < Cc; ++c) {
            if (nvec[c] == 0) continue;  // MATLAB forms 0/0 here and never reads it
            T ctot = Zc[c];
            for (std::size_t i = 0; i < M; ++i) ctot += w(i, c);
            if (ctot == zero) throw NumericError("pfqn_mvaldmx: zero total residence time");
            x[c] = num_traits<T>::from_int(nvec[c]) / ctot;
        }

        for (std::size_t i = 0; i < M; ++i) {
            for (long k = 1; k <= nc; ++k) {
                T acc = zero;
                for (std::size_t c = 0; c < Cc; ++c) {
                    if (nvec[c] == 0) continue;
                    const std::size_t idx_c = idx - prods[c];
                    acc += Dc(i, c) * ec.EC(i, static_cast<std::size_t>(k - 1)) * x[c] *
                           Pc[(idx_c * M + i) * K + static_cast<std::size_t>(k - 1)];
                }
                Pc[(idx * M + i) * K + static_cast<std::size_t>(k)] = acc;
            }
            T acc = zero;
            for (long k = 1; k <= nc; ++k) acc += Pc[(idx * M + i) * K + static_cast<std::size_t>(k)];
            T p0 = one - acc;
            if (p0 < tiny) p0 = tiny;  // MATLAB: max(eps, 1 - sum)
            Pc[(idx * M + i) * K + 0] = p0;
        }

        long last_nnz = -1;
        for (long c = static_cast<long>(Cc) - 1; c >= 0; --c)
            if (nvec[c] != 0) {
                last_nnz = c;
                break;
            }
        if (last_nnz >= 0) {
            bool prefixFull = true;
            for (long c = 0; c < last_nnz; ++c)
                if (nvec[c] != Nc[c]) {
                    prefixFull = false;
                    break;
                }
            bool suffixEmpty = true;
            for (std::size_t c = static_cast<std::size_t>(last_nnz) + 1; c < Cc; ++c)
                if (nvec[c] != 0) {
                    suffixEmpty = false;
                    break;
                }
            if (prefixFull && suffixEmpty) {
                const T& xr = x[static_cast<std::size_t>(last_nnz)];
                if (xr > zero) res.G /= xr;
            }
        }

        more = next_pop(nvec, Nc);
    }

    const std::size_t idxN = total - 1;

    for (std::size_t c = 0; c < Cc; ++c) {
        const std::size_t r = closedClasses[c];
        res.XN[r] = x[c];
        for (std::size_t i = 0; i < M; ++i) {
            res.CN(i, r) = w(i, c);
            res.QN(i, r) = x[c] * w(i, c);
            if (Nc[c] > 0) {
                const std::size_t idx_c = idxN - prods[c];
                T u = zero;
                for (long k = 1; k <= Nct; ++k) {
                    const std::size_t kk = static_cast<std::size_t>(k - 1);
                    if (ec.E(i, kk) == zero) throw NumericError("pfqn_mvaldmx: vanishing effective capacity");
                    u += Dc(i, c) * x[c] * ec.Eprime(i, kk) / ec.E(i, kk) * Pc[(idx_c * M + i) * K + kk];
                }
                res.UN(i, r) = u;
            }
        }
    }

    for (std::size_t r : openClasses) {
        res.XN[r] = lambda[r];
        for (std::size_t i = 0; i < M; ++i) {
            T q = zero, u = zero;
            for (long k = 0; k <= Nct; ++k) {
                const std::size_t kk = static_cast<std::size_t>(k);
                q += lambda[r] * D(i, r) * num_traits<T>::from_int(k + 1) * ec.EC(i, kk) *
                     Pc[(idxN * M + i) * K + kk];
                if (ec.E(i, kk + 1) == zero) throw NumericError("pfqn_mvaldmx: vanishing effective capacity");
                u += lambda[r] * ec.Eprime(i, kk + 1) / ec.E(i, kk + 1) * Pc[(idxN * M + i) * K + kk];
            }
            res.QN(i, r) = q;
            if (lambda[r] == zero)
                throw NumericError("pfqn_mvaldmx: an open class must have a positive arrival rate");
            res.CN(i, r) = q / lambda[r];
            res.UN(i, r) = u;
        }
    }

    res.lG = Cc > 0 ? num_traits<T>::log_as_double(res.G) : std::numeric_limits<double>::quiet_NaN();
    return res;
}

/**
 * Exact MVA for mixed open/closed networks with multiserver stations.
 *
 * Port of matlab/src/api/pfqn/pfqn_mvaldms.m: builds the multiserver rates
 * mu(i,k) = min(k, S(i)), calls pfqn_mvaldmx and replaces its utilizations by
 * the per-server utilization law U(i,r) = X(r) D(i,r) / S(i).
 *
 * @param S (M) servers per station, INF_SERVERS for an infinite server
 * @param lambda (R) arrival rates, zero on the closed classes
 * @param D (M x R) service demands
 * @param N (R) populations, negative on the open classes
 * @param Z (K x R) think times
 */
template <class T>
MvaResult<T> pfqn_mvaldms(const std::vector<T>& lambda, const Matrix<T>& D,
                          const std::vector<int>& N, const Matrix<T>& Z, const std::vector<int>& S) {
    const std::size_t M = D.rows();
    const std::size_t R = N.size();
    if (S.size() != M) throw InputError("pfqn_mvaldms: server-count vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    long Nct = 0;
    for (std::size_t r = 0; r < R; ++r)
        if (!is_open_class(N[r])) Nct += N[r];

    Matrix<T> mu(M, static_cast<std::size_t>(Nct > 0 ? Nct : 1), num_traits<T>::from_int(1));
    for (std::size_t i = 0; i < M; ++i)
        for (long k = 1; k <= Nct; ++k) {
            const long c = S[i] == INF_SERVERS ? k : (k < S[i] ? k : S[i]);
            mu(i, static_cast<std::size_t>(k - 1)) = num_traits<T>::from_int(c);
        }

    MvaResult<T> res = pfqn_mvaldmx(lambda, D, N, Z, mu);

    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) {
            const T rate = is_open_class(N[r]) ? lambda[r] : res.XN[r];
            // An infinite server has no server count to divide by: the law
            // degenerates to the mean number of busy servers, X * D.
            res.UN(i, r) = S[i] == INF_SERVERS
                               ? T(rate * D(i, r))
                               : T(rate * D(i, r) / num_traits<T>::from_int(S[i]));
            if (rate == zero) res.UN(i, r) = zero;
        }
    }
    return res;
}

// ---------------------------------------------------------------------------
// pfqn_mvams
// ---------------------------------------------------------------------------

/**
 * General-purpose exact MVA for mixed networks with multiserver stations.
 *
 * @param lambda (R) arrival rates, zero on closed classes; may be empty when
 *               the model has no open class
 * @param L      (M x R) service demands
 * @param N      (R) population, OPEN_CLASS on open classes
 * @param Z      (K x R) think times, summed over rows; may be empty
 * @param mi     (M) station multiplicities; empty for all ones
 * @param S      (M) servers per station, INF_SERVERS for an infinite server;
 *               empty for all ones
 *
 * The returned CN is the (M x R) per-station residence time of the pfqn_mva
 * contract in every branch, and UN the (M x R) per-class utilization; see the
 * contract notes at the top of this header for the two points at which that
 * differs from MATLAB. In the mixed multiserver branch the normalizing
 * constant is not available: G is 0 and lG is NaN, as in MATLAB.
 *
 * Standard arrival theorem throughout. For the interlocked-flow correction of
 * Franks (1999), Ch. 4, Eq. (4.7), call pfqn_mvams_ilock instead.
 */
template <class T>
MvaResult<T> pfqn_mvams(const std::vector<T>& lambda, const Matrix<T>& L,
                        const std::vector<int>& N, const Matrix<T>& Z, const std::vector<int>& mi,
                        const std::vector<int>& S) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mvams: demand matrix and population vector disagree on the class count");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_mvams: multiplicity vector has the wrong length");
    if (!S.empty() && S.size() != M)
        throw InputError("pfqn_mvams: server-count vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<int> Sv = S.empty() ? std::vector<int>(M, 1) : S;
    std::vector<int> miv = mi.empty() ? std::vector<int>(M, 1) : mi;
    std::vector<T> lam = lambda.empty() ? std::vector<T>(R, zero) : lambda;
    if (lam.size() != R)
        throw InputError("pfqn_mvams: arrival-rate vector and population vector disagree on the class count");

    long Ntot = 0;
    bool hasOpenClasses = false;
    for (std::size_t r = 0; r < R; ++r) {
        if (is_open_class(N[r]))
            hasOpenClasses = true;
        else
            Ntot += N[r];
    }
    // An infinite server is not counted here, mirroring MATLAB's isfinite guard.
    bool hasMultiServer = false;
    for (std::size_t i = 0; i < M; ++i)
        if (Sv[i] != INF_SERVERS && Sv[i] > 1) {
            hasMultiServer = true;
            break;
        }
    bool hasReplicas = false;
    for (std::size_t i = 0; i < M; ++i)
        if (miv[i] != 1) {
            hasReplicas = true;
            break;
        }

    if (!hasMultiServer) {
        if (hasOpenClasses) return pfqn_mvamx(lam, L, N, Z, miv);
        return pfqn_mva(L, N, Z, miv);
    }

    if (hasReplicas)
        throw InputError("pfqn_mvams: queue replicas are not available in exact MVA with multiserver stations");

    if (hasOpenClasses) {
        MvaResult<T> res = pfqn_mvaldms(lam, L, N, Z, Sv);
        // MATLAB discards the normalizing constant on this branch.
        res.G = zero;
        res.lG = std::numeric_limits<double>::quiet_NaN();
        return res;
    }

    Matrix<T> mu(M, static_cast<std::size_t>(Ntot > 0 ? Ntot : 1), one);
    for (std::size_t i = 0; i < M; ++i)
        for (long k = 1; k <= Ntot; ++k) {
            const long c = Sv[i] == INF_SERVERS ? k : (k < Sv[i] ? k : Sv[i]);
            mu(i, static_cast<std::size_t>(k - 1)) = num_traits<T>::from_int(c);
        }

    const MvaLdResult<T> ld = pfqn_mvald(L, N, Z, mu);

    MvaResult<T> res;
    res.XN = ld.XN;
    res.QN = ld.QN;
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.G = ld.G;
    res.lG = ld.lG;

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] > 0) {
                res.CN(i, r) = res.QN(i, r) / res.XN[r];
                res.UN(i, r) = Sv[i] == INF_SERVERS
                                   ? T(res.XN[r] * L(i, r))
                                   : T(res.XN[r] * L(i, r) / num_traits<T>::from_int(Sv[i]));
            } else {
                // An absent class, as in pfqn_mva.
                res.CN(i, r) = L(i, r) * num_traits<T>::from_int(miv[i]);
            }
        }
    return res;
}

/** Overload with unit multiplicities. */
template <class T>
MvaResult<T> pfqn_mvams(const std::vector<T>& lambda, const Matrix<T>& L,
                        const std::vector<int>& N, const Matrix<T>& Z, const std::vector<int>& S) {
    return pfqn_mvams(lambda, L, N, Z, std::vector<int>(), S);
}

/** Overload with unit multiplicities and a single server everywhere. */
template <class T>
MvaResult<T> pfqn_mvams(const std::vector<T>& lambda, const Matrix<T>& L,
                        const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_mvams(lambda, L, N, Z, std::vector<int>(), std::vector<int>());
}

/**
 * MVA entry point for models carrying the interlocked-flow correction.
 *
 * The interlock of Franks (1999), Ch. 4, Eq. (4.7) is defined only for closed
 * single-server models, so that is the one shape accepted here; anything else is
 * refused rather than served without the correction. Models with no interlock go to
 * pfqn_mvams.
 *
 * @param IL (R x R) interlock matrix, see pfqn_mva_ilock. Required.
 */
template <class T>
MvaResult<T> pfqn_mvams_ilock(const std::vector<T>& lambda, const Matrix<T>& L,
                              const std::vector<int>& N, const Matrix<T>& Z,
                              const std::vector<int>& mi, const std::vector<int>& S,
                              const Matrix<T>& IL) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (IL.empty())
        throw InputError("pfqn_mvams_ilock: an interlock matrix is required; use pfqn_mvams for the standard arrival theorem");
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mvams_ilock: demand matrix and population vector disagree on the class count");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_mvams_ilock: multiplicity vector has the wrong length");
    if (!S.empty() && S.size() != M)
        throw InputError("pfqn_mvams_ilock: server-count vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    std::vector<int> Sv = S.empty() ? std::vector<int>(M, 1) : S;
    std::vector<int> miv = mi.empty() ? std::vector<int>(M, 1) : mi;
    std::vector<T> lam = lambda.empty() ? std::vector<T>(R, zero) : lambda;
    if (lam.size() != R)
        throw InputError("pfqn_mvams_ilock: arrival-rate vector and population vector disagree on the class count");

    for (std::size_t r = 0; r < R; ++r)
        if (is_open_class(N[r]) || lam[r] != zero)
            throw InputError("pfqn_mvams_ilock: the interlock correction is available in exact MVA "
                             "for closed single-server models only; use an AMVA method");
    for (std::size_t i = 0; i < M; ++i)
        if (Sv[i] != INF_SERVERS && Sv[i] > 1)
            throw InputError("pfqn_mvams_ilock: the interlock correction is available in exact MVA "
                             "for closed single-server models only; use an AMVA method");

    return pfqn_mva_ilock(L, N, Z, miv, IL);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVAMS_H
