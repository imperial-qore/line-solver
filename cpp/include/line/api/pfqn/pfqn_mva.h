/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVA_H
#define LINE_API_PFQN_MVA_H

/**
 * Exact Mean Value Analysis for closed product-form networks
 * (Reiser and Lavenberg 1980).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mva.m, cross-checked against
 * mp_pfqn's mva/mva-multi.c for the exact path. The recursion over the
 * population lattice is
 *
 *   C(i,r|n) = L(i,r) (mi(i) + Q(i|n - e_r))
 *   X(r|n)   = n_r / (Z_r + sum_i C(i,r|n))
 *   Q(i,r|n) = X(r|n) C(i,r|n)
 *
 * every operation of which stays in the field of the inputs, so the algorithm
 * is exact in rational arithmetic with no reformulation.
 *
 * Normalizing constant: MATLAB accumulates lG as -sum log X along the lattice
 * path (0 -> N) that fills one class at a time. Logs do not exist in an exact
 * field, so the port accumulates the product of the reciprocals instead and
 * takes the log once at the end, of a value that is still exact. The two agree
 * to rounding in double.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct MvaResult {
    std::vector<T> XN;  ///< (R) per-class throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> CN;       ///< (M x R) residence time
    T G;                ///< normalizing constant
    double lG;          ///< log of the normalizing constant
};

/**
 * @param L  (M x R) service demands
 * @param N  (R) population per class
 * @param Z  (K x R) think times, summed over rows; empty for no delay
 * @param mi (M) additive term of the residence-time recursion
 *        C(i,s)=L(i,s)*(mi[i]+Qarv), 1 for a queueing station; empty for all ones.
 *        THIS IS NOT A SERVER COUNT: mi[i]=c inflates the residence time by c
 *        rather than adding c servers. For multiserver stations call
 *        pfqn_mvams(lambda, L, N, Z, mi, S), which passes S to the load-dependent
 *        recursion with mu(i,n)=min(n,S(i)).
 *
 * Standard arrival theorem. For the interlocked-flow correction of Franks (1999),
 * Ch. 4, Eq. (4.7), call pfqn_mva_ilock instead.
 */
template <class T>
MvaResult<T> pfqn_mva(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                      const std::vector<int>& mi) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mva: demand matrix and population vector disagree on the class count");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_mva: multiplicity vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    MvaResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.G = one;
    res.lG = 0.0;

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_mva: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0) return res;  // empty closed population: nothing to compute

    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_mva: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    std::vector<T> multi(M, one);
    for (std::size_t i = 0; i < mi.size(); ++i) multi[i] = num_traits<T>::from_int(mi[i]);

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // lattice indexing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<T> Q(total * M, zero);

    std::vector<int> n(R, 0);
    bool more = true;
    while (more) {
        int npop = 0;
        for (int v : n) npop += v;
        if (npop > 0) {
            const std::size_t idx = pop_index(n, prods);
            for (std::size_t s = 0; s < R; ++s) {
                // empty-class zero-population-row rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                const std::size_t idx_1s = n[s] > 0 ? idx - prods[s] : 0;
                T ctot = Zsum[s];
                for (std::size_t i = 0; i < M; ++i) {
                    const T qarv = Q[idx_1s * M + i];
                    res.CN(i, s) = L(i, s) * (multi[i] + qarv);
                    ctot += res.CN(i, s);
                }
                if (ctot == zero) throw NumericError("pfqn_mva: zero total residence time");
                res.XN[s] = num_traits<T>::from_int(n[s]) / ctot;
                for (std::size_t i = 0; i < M; ++i) {
                    res.QN(i, s) = res.XN[s] * res.CN(i, s);
                    Q[idx * M + i] += res.QN(i, s);
                }
            }

            // normalizing-constant accumulation rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
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
                    const T& x = res.XN[static_cast<std::size_t>(last_nnz)];
                    if (x == zero) throw NumericError("pfqn_mva: zero throughput on the G path");
                    res.G /= x;
                }
            }
        }
        more = next_pop(n, N);
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.UN(i, r) = res.XN[r] * L(i, r);

    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

/**
 * Exact MVA recursion carrying the interlocked-flow correction.
 *
 * The correction replaces the arrival theorem term Q(n-1_s,i) by a per-class weighted
 * sum, so the recursion has to carry per-class queue lengths that pfqn_mva does not
 * need. Closed single-server models only.
 *
 * The discounted arrival-instant queue is floored at the in-service component, as in
 * lqns MVA::queueOnly_adjusted, so the correction damps itself out as a station
 * saturates. That is a self-limiting guard, NOT a hard capacity test:
 * sum_s XN[s]*L(i,s) <= mi[i] is still asserted nowhere.
 * See git show 449847e7b:_kb/log.md.
 *
 * @param L  (M x R) service demands
 * @param N  (R) population per class
 * @param Z  (K x R) think times, summed over rows; empty for no delay
 * @param mi (M) additive term of the residence-time recursion
 *        C(i,s)=L(i,s)*(mi[i]+Qarv), 1 for a queueing station; empty for all ones.
 *        THIS IS NOT A SERVER COUNT: mi[i]=c inflates the residence time by c
 *        rather than adding c servers. For multiserver stations call
 *        pfqn_mvams(lambda, L, N, Z, mi, S), which passes S to the load-dependent
 *        recursion with mu(i,n)=min(n,S(i)).
 * @param IL (R x R) interlock matrix of Franks (1999), Eq. (4.7): IL(r,s) is the
 *        share of the class-s queue that a class-r arrival cannot see, because that
 *        work was itself caused by the class-r request.
 *        Required. The model is outside product form under it, so G and lG are not
 *        meaningful and come back as one and zero.
 */
template <class T>
MvaResult<T> pfqn_mva_ilock(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                            const std::vector<int>& mi, const Matrix<T>& IL) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_mva_ilock: demand matrix and population vector disagree on the class count");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_mva_ilock: multiplicity vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    MvaResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.G = one;
    res.lG = 0.0;

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_mva_ilock: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0) return res;  // empty closed population: nothing to compute

    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_mva_ilock: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    std::vector<T> multi(M, one);
    for (std::size_t i = 0; i < mi.size(); ++i) multi[i] = num_traits<T>::from_int(mi[i]);

    if (IL.empty())
        throw InputError("pfqn_mva_ilock: an interlock matrix is required; use pfqn_mva for the standard arrival theorem");
    Matrix<T> ILw;
    {
        if (IL.rows() != R || IL.cols() != R)
            throw InputError("pfqn_mva_ilock: the interlock matrix must be nclasses x nclasses");
        ILw = Matrix<T>(R, R, zero);
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) {
                T w = (r == s) ? one : T(one - IL(r, s));
                if (w < zero) w = zero;
                if (w > one) w = one;
                ILw(r, s) = w;
            }
    }

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // lattice indexing rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<T> Q(total * M, zero);
    // per-class queue lengths, needed by the interlock
    std::vector<T> Qc(total * M * R, zero);
    // per-class in-service component, the interlock's floor
    std::vector<T> Uc(total * M * R, zero);

    std::vector<int> n(R, 0);
    bool more = true;
    while (more) {
        int npop = 0;
        for (int v : n) npop += v;
        if (npop > 0) {
            const std::size_t idx = pop_index(n, prods);
            for (std::size_t s = 0; s < R; ++s) {
                // empty-class zero-population-row rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                const std::size_t idx_1s = n[s] > 0 ? idx - prods[s] : 0;
                T ctot = Zsum[s];
                for (std::size_t i = 0; i < M; ++i) {
                    T qarv = zero;
                    for (std::size_t r = 0; r < R; ++r) {
                        // In-service protection, as in lqns MVA::queueOnly_adjusted: the
                        // discount bites on the WAITING part only, never on the job already
                        // in service, so it damps itself out as the station saturates.
                        const T disc = ILw(s, r) * Qc[(idx_1s * M + i) * R + r];
                        const T inSvc = Uc[(idx_1s * M + i) * R + r];
                        qarv += disc > inSvc ? disc : inSvc;
                    }
                    res.CN(i, s) = L(i, s) * (multi[i] + qarv);
                    ctot += res.CN(i, s);
                }
                if (ctot == zero) throw NumericError("pfqn_mva_ilock: zero total residence time");
                res.XN[s] = num_traits<T>::from_int(n[s]) / ctot;
                for (std::size_t i = 0; i < M; ++i) {
                    res.QN(i, s) = res.XN[s] * res.CN(i, s);
                    Q[idx * M + i] += res.QN(i, s);
                    Qc[(idx * M + i) * R + s] = res.QN(i, s);
                    Uc[(idx * M + i) * R + s] = res.XN[s] * L(i, s);
                }
            }

            // the interlock leaves the model outside product form, so G is never
            // accumulated and res.G/res.lG keep their neutral initial values
        }
        more = next_pop(n, N);
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.UN(i, r) = res.XN[r] * L(i, r);

    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

template <class T>
MvaResult<T> pfqn_mva(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_mva(L, N, Z, std::vector<int>());
}

template <class T>
MvaResult<T> pfqn_mva(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_mva(L, N, Matrix<T>(), std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVA_H
