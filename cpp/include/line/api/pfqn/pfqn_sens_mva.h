/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_MVA_H
#define LINE_API_PFQN_SENS_MVA_H

/**
 * Exact per-station queue-length variances and covariances of a closed
 * product-form network, by the MVA-like moment recursion of de Souza e Silva
 * and Muntz (IEEE TC 37(9):1125-1129, 1988, Corollary 1).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_mva.m. Writing
 * W(k,i;v,j|n) = Cov[n(i,k),n(j,v)], differentiating the Reiser-Lavenberg
 * equation and rescaling gives
 *
 *   W(k,i;v,j|n) = Q(j,v|n) (Q(i,k|n-e_v) - Q(i,k|n))
 *                + [i==j && k==v] Q(j,v|n)
 *                + X(v|n) L(j,v) sum_t W(k,i;t,j|n-e_v)
 *
 * with W(.|0) = 0. Only the same-station case i==j is evaluated, because only
 * then does the inner sum stay at the station and close on the single scalar
 * Ssum(j,k|n) = sum_t W(k,j;t,j|n) carried along the lattice. The
 * cross-station blocks are not self-contained and are left to pfqn_sens.
 *
 * Arithmetic. Every operation is an addition, a multiplication or a division
 * by a lattice quantity, so the recursion stays in the field of the inputs and
 * instantiates at line::Rational with no reformulation: the covariances of a
 * rational model are exact rationals. There is deliberately no
 * static_assert(has_transcendental) here.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/**
 * Radix weights of the MVA population lattice, class R-1 varying fastest.
 * This is the transpose of line::plane_sizes and matches the `prods` vector of
 * pfqn_mva.m, whose lattice index the sensitivity routines must reproduce
 * entry by entry so that their base measures agree with pfqn_mva's.
 */
inline std::vector<std::size_t> sens_lattice_radix(const std::vector<int>& N) {
    const std::size_t R = N.size();
    std::vector<std::size_t> radix(R, 1);
    for (std::size_t w = R; w-- > 0;) {
        if (w + 1 == R)
            radix[w] = 1;
        else
            radix[w] = radix[w + 1] * static_cast<std::size_t>(N[w + 1] + 1);
    }
    return radix;
}

/** Decode a lattice index back into a population vector. */
inline std::vector<int> sens_lattice_decode(std::size_t k, const std::vector<int>& N,
                                            const std::vector<std::size_t>& radix) {
    const std::size_t R = N.size();
    std::vector<int> n(R, 0);
    for (std::size_t w = 0; w < R; ++w)
        n[w] = static_cast<int>((k / radix[w]) % static_cast<std::size_t>(N[w] + 1));
    return n;
}

template <class T>
struct SensMvaResult {
    std::vector<T> XN;          ///< (R) throughput, identical to pfqn_mva
    Matrix<T> QN;               ///< (M x R) mean queue length
    Matrix<T> UN;               ///< (M x R) utilization
    Matrix<T> CN;               ///< (M x R) residence time (MATLAB field .R)
    std::vector<Matrix<T>> QCov;  ///< (M) matrices R x R, QCov[i](r,s) = Cov[n(i,r),n(i,s)]
    Matrix<T> QVar;             ///< (M x R) Var[n(i,r)]
    std::vector<T> QTotVar;     ///< (M) Var[sum_r n(i,r)]
    T QCovAsym;                 ///< raw asymmetry of QCov before symmetrization
};

/**
 * @param L  (M x R) service demands
 * @param N  (R) population per class, closed only
 * @param Z  (R) think times, empty for none
 * @param mi (M) station multiplicities, empty for all ones
 */
template <class T>
SensMvaResult<T> pfqn_sens_mva(const Matrix<T>& L, const std::vector<int>& N,
                               const std::vector<T>& Z, const std::vector<int>& mi) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_sens_mva: demand matrix and population vector disagree on the class count");
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_sens_mva: think-time vector has the wrong length");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_sens_mva: multiplicity vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);

    SensMvaResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.QCov.assign(M, Matrix<T>(R, R, zero));
    res.QVar = Matrix<T>(M, R, zero);
    res.QTotVar.assign(M, zero);
    res.QCovAsym = zero;

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens_mva: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0 || R == 0) return res;

    const auto Zr = [&](std::size_t r) -> T { return Z.empty() ? zero : Z[r]; };
    const auto miT = [&](std::size_t i) -> T {
        return num_traits<T>::from_int(mi.empty() ? 1 : mi[i]);
    };

    const std::vector<std::size_t> radix = sens_lattice_radix(N);
    const std::size_t totpop = population_count(N);

    // Lattice state. Qcls and Ssum are the only per-population tables the
    // recursion reads back; Qtot and Xall are needed by the MVA step itself.
    Matrix<T> Qtot(totpop, M, zero);
    std::vector<Matrix<T>> Qcls(totpop, Matrix<T>(M, R, zero));
    Matrix<T> Xall(totpop, R, zero);
    std::vector<Matrix<T>> Ssum(totpop, Matrix<T>(M, R, zero));

    std::vector<std::size_t> rows(R, 0);

    for (std::size_t k = 1; k < totpop; ++k) {
        const std::vector<int> n = sens_lattice_decode(k, N, radix);

        // ---- mean value analysis step -------------------------------------
        for (std::size_t s = 0; s < R; ++s) {
            // empty-population index collapse rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
            const std::size_t row = n[s] > 0 ? k - radix[s] : 0;
            rows[s] = row;
            T CNtot = zero;
            for (std::size_t i = 0; i < M; ++i) {
                res.CN(i, s) = L(i, s) * (miT(i) + Qtot(row, i));
                CNtot += res.CN(i, s);
            }
            const T den = Zr(s) + CNtot;
            if (den == zero) {
                res.XN[s] = zero;
            } else {
                res.XN[s] = num_traits<T>::from_int(n[s]) / den;
            }
            Xall(k, s) = res.XN[s];
            for (std::size_t i = 0; i < M; ++i) {
                res.QN(i, s) = res.XN[s] * res.CN(i, s);
                Qcls[k](i, s) = res.QN(i, s);
                Qtot(k, i) += res.QN(i, s);
            }
        }

        // ---- moment step ---------------------------------------------------
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t kk = 0; kk < R; ++kk) {
                const T Qjk = Qcls[k](j, kk);
                T sk = zero;
                for (std::size_t t = 0; t < R; ++t) {
                    const T Qjt = Qcls[k](j, t);
                    T wkt = Qjt * (Qcls[rows[t]](j, kk) - Qjk);
                    if (kk == t) wkt += Qjt;
                    wkt += Xall(k, t) * L(j, t) * Ssum[rows[t]](j, kk);
                    res.QCov[j](kk, t) = wkt;
                    sk += wkt;
                }
                Ssum[k](j, kk) = sk;
            }
        }
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.UN(i, r) = res.XN[r] * L(i, r);

    // covariance-triangle symmetrization rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t s = r + 1; s < R; ++s) {
                const T d = num_abs(T(res.QCov[i](r, s) - res.QCov[i](s, r)));
                if (d > res.QCovAsym) res.QCovAsym = d;
                const T avg = (res.QCov[i](r, s) + res.QCov[i](s, r)) / num_traits<T>::from_int(2);
                res.QCov[i](r, s) = avg;
                res.QCov[i](s, r) = avg;
            }
        }
        T tot = zero;
        for (std::size_t r = 0; r < R; ++r) {
            res.QVar(i, r) = res.QCov[i](r, r);
            for (std::size_t s = 0; s < R; ++s) tot += res.QCov[i](r, s);
        }
        res.QTotVar[i] = tot;
    }

    return res;
}

/** pfqn_sens_mva with unit multiplicities. */
template <class T>
SensMvaResult<T> pfqn_sens_mva(const Matrix<T>& L, const std::vector<int>& N,
                               const std::vector<T>& Z) {
    return pfqn_sens_mva(L, N, Z, std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_MVA_H
