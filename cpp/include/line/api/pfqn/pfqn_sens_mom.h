/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_MOM_H
#define LINE_API_PFQN_SENS_MOM_H

/**
 * Exact moments E[Q], Var[Q], E[Q^2] and E[Q^3] of the grouped queue lengths of
 * a closed product-form network, by second-order differentiation of the MVA
 * recursion.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_mom.m. Theorem 3.1 of
 * Strelen (Performance Evaluation 11:127-142, 1990) states that one further
 * factor Q_i in a moment costs one differentiation with respect to x_i, the
 * reciprocal capacity of station i, so with m_i = E[Q_i] and equation (3.2)
 *
 *   Var[Q_i]     = x_i dm_i/dx_i
 *   Cov[Q_i,Q_j] = x_j dm_i/dx_j
 *   E[Q_i^2]     = x_i dm_i/dx_i + m_i^2
 *   E[Q_i^3]     = x_i^2 d2m_i/dx_i^2 + (x_i + 3 x_i m_i) dm_i/dx_i + m_i^3
 *
 * The third moment therefore needs the SECOND derivative of the recursion,
 * which is what this routine adds over pfqn_sens_mva and pfqn_sens. The
 * parameter y_(h,g) rescales the demands of the classes of group g at station
 * h; at y = 1 the y-derivatives are exactly the scaled x-derivatives that
 * (3.2) asks for, because a pure rescaling gives d/dy = x d/dx and
 * d2/dy2 = x^2 d2/dx2. The grouping is the class subset T of Theorem 1 of
 * Akyildiz and Strelen (IEEE TComm 39(6):828-832, 1991); groups all equal to 0
 * gives the per-station totals of Strelen's x_i, one group per class gives the
 * per-class moments.
 *
 * Arithmetic. The recursion is field-only, so the moments instantiate at
 * line::Rational and are exact rationals. The skewness is the one derived
 * quantity that leaves the field (it divides by Var^1.5), so it is returned as
 * a double rather than as a T, and the header carries no transcendental gate.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_sens_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct SensMomResult {
    std::vector<T> XN;  ///< (R) throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> CN;       ///< (M x R) residence time (MATLAB field .R)

    Matrix<T> m;    ///< (M x G) E[Q_(i,g)]
    Matrix<T> d2m;  ///< (M x G) the scaled pure second derivative
    /// (M*G x M*G) Cov[Q_(i,g),Q_(j,g')] at row i*G+g, column j*G+g'.
    Matrix<T> Cov;
    Matrix<T> dm;  ///< (M*G x M*G) the raw scaled first derivative, same layout
    Matrix<T> Var;   ///< (M x G)
    Matrix<T> M2;    ///< (M x G) E[Q^2]
    Matrix<T> M3;    ///< (M x G) E[Q^3]
    Matrix<double> Skew;  ///< (M x G) skewness, NaN where the variance vanishes
    T CovAsym;       ///< raw asymmetry of Cov before symmetrization
};

/**
 * @param L      (M x R) service demands
 * @param N      (R) population per class, closed only
 * @param Z      (R) think times, empty for none
 * @param mi     (M) station multiplicities, empty for all ones
 * @param groups (R) 0-based group label of each class; empty means one group
 *               holding every class, i.e. the per-station totals
 */
template <class T>
SensMomResult<T> pfqn_sens_mom(const Matrix<T>& L, const std::vector<int>& N,
                               const std::vector<T>& Z, const std::vector<int>& mi,
                               const std::vector<int>& groups) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_sens_mom: demand matrix and population vector disagree on the class count");
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_sens_mom: think-time vector has the wrong length");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_sens_mom: multiplicity vector has the wrong length");

    std::vector<int> grp = groups.empty() ? std::vector<int>(R, 0) : groups;
    if (grp.size() != R)
        throw InputError("pfqn_sens_mom: groups must have one label per class");
    std::size_t G = 0;
    for (int g : grp) {
        if (g < 0) throw InputError("pfqn_sens_mom: group labels start at zero");
        if (static_cast<std::size_t>(g) + 1 > G) G = static_cast<std::size_t>(g) + 1;
    }
    {
        std::vector<bool> seen(G, false);
        for (int g : grp) seen[static_cast<std::size_t>(g)] = true;
        for (std::size_t g = 0; g < G; ++g)
            if (!seen[g])
                throw InputError("pfqn_sens_mom: groups must label the classes consecutively with no empty group");
    }

    const T zero = num_traits<T>::from_int(0);
    const std::size_t P = M * G;

    SensMomResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    res.m = Matrix<T>(M, G, zero);
    res.d2m = Matrix<T>(M, G, zero);
    res.Cov = Matrix<T>(P, P, zero);
    res.dm = Matrix<T>(P, P, zero);
    res.Var = Matrix<T>(M, G, zero);
    res.M2 = Matrix<T>(M, G, zero);
    res.M3 = Matrix<T>(M, G, zero);
    res.Skew = Matrix<double>(M, G, 0.0);
    res.CovAsym = zero;

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens_mom: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0 || R == 0) return res;

    const auto Zr = [&](std::size_t r) -> T { return Z.empty() ? zero : Z[r]; };
    const auto miT = [&](std::size_t i) -> T {
        return num_traits<T>::from_int(mi.empty() ? 1 : mi[i]);
    };
    const auto pidx = [&](std::size_t h, std::size_t g) { return h * G + g; };

    const std::vector<std::size_t> radix = sens_lattice_radix(N);
    const std::size_t totpop = population_count(N);

    Matrix<T> Qtot(totpop, M, zero);
    std::vector<Matrix<T>> D1(totpop, Matrix<T>(M, P, zero));
    std::vector<Matrix<T>> D2(totpop, Matrix<T>(M, P, zero));

    // The group accumulators describe one population only and are overwritten
    // on every sweep, so at the end of the walk they hold the values at N.
    Matrix<T> Qg(M, G, zero);
    std::vector<Matrix<T>> D1Qg(P, Matrix<T>(M, G, zero));
    std::vector<Matrix<T>> D2Qg(P, Matrix<T>(M, G, zero));

    std::vector<T> Cs(M, zero), dCNtot(P, zero), d2CNtot(P, zero), dX(P, zero), d2X(P, zero);
    Matrix<T> dCs(M, P, zero), d2Cs(M, P, zero);

    for (std::size_t k = 1; k < totpop; ++k) {
        const std::vector<int> n = sens_lattice_decode(k, N, radix);
        Qg.fill(zero);
        for (std::size_t p = 0; p < P; ++p) {
            D1Qg[p].fill(zero);
            D2Qg[p].fill(zero);
        }
        for (std::size_t s = 0; s < R; ++s) {
            const std::size_t row = n[s] > 0 ? k - radix[s] : 0;
            const std::size_t gs = static_cast<std::size_t>(grp[s]);

            // ---- residence times and their first two derivatives ----------
            T CNtot = zero;
            for (std::size_t p = 0; p < P; ++p) {
                dCNtot[p] = zero;
                d2CNtot[p] = zero;
            }
            for (std::size_t i = 0; i < M; ++i) {
                const T A = miT(i) + Qtot(row, i);
                Cs[i] = L(i, s) * A;
                res.CN(i, s) = Cs[i];
                CNtot += Cs[i];
                for (std::size_t p = 0; p < P; ++p) {
                    const T dA = D1[row](i, p);
                    const T d2A = D2[row](i, p);
                    if (p == pidx(i, gs)) {
                        dCs(i, p) = L(i, s) * (A + dA);
                        d2Cs(i, p) = L(i, s) * (num_traits<T>::from_int(2) * dA + d2A);
                    } else {
                        dCs(i, p) = L(i, s) * dA;
                        d2Cs(i, p) = L(i, s) * d2A;
                    }
                    dCNtot[p] += dCs(i, p);
                    d2CNtot[p] += d2Cs(i, p);
                }
            }

            // ---- throughput and its first two derivatives ------------------
            const T den = Zr(s) + CNtot;
            const T nsT = num_traits<T>::from_int(n[s]);
            if (den == zero) {
                res.XN[s] = zero;
                for (std::size_t p = 0; p < P; ++p) {
                    dX[p] = zero;
                    d2X[p] = zero;
                }
            } else {
                const T den2 = den * den;
                const T den3 = den2 * den;
                res.XN[s] = nsT / den;
                for (std::size_t p = 0; p < P; ++p) {
                    dX[p] = -nsT * dCNtot[p] / den2;
                    d2X[p] = -nsT * d2CNtot[p] / den2 +
                             num_traits<T>::from_int(2) * nsT * dCNtot[p] * dCNtot[p] / den3;
                }
            }

            // ---- queue lengths ---------------------------------------------
            for (std::size_t i = 0; i < M; ++i) {
                res.QN(i, s) = res.XN[s] * Cs[i];
                Qtot(k, i) += res.QN(i, s);
                Qg(i, gs) += res.QN(i, s);
                for (std::size_t p = 0; p < P; ++p) {
                    const T dQ = dX[p] * Cs[i] + res.XN[s] * dCs(i, p);
                    const T d2Q = d2X[p] * Cs[i] +
                                  num_traits<T>::from_int(2) * dX[p] * dCs(i, p) +
                                  res.XN[s] * d2Cs(i, p);
                    D1[k](i, p) += dQ;
                    D2[k](i, p) += d2Q;
                    D1Qg[p](i, gs) += dQ;
                    D2Qg[p](i, gs) += d2Q;
                }
            }
        }
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.UN(i, r) = res.XN[r] * L(i, r);

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t g = 0; g < G; ++g) {
            res.m(i, g) = Qg(i, g);
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t g2 = 0; g2 < G; ++g2)
                    res.dm(pidx(i, g), pidx(j, g2)) = D1Qg[pidx(j, g2)](i, g);
            res.d2m(i, g) = D2Qg[pidx(i, g)](i, g);
        }

    // Cov((i,g),(j,g')) and its transpose are distinct expressions for the same
    // quantity; report the raw disagreement, then symmetrize.
    for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b < P; ++b) {
            const T d = num_abs(T(res.dm(a, b) - res.dm(b, a)));
            if (d > res.CovAsym) res.CovAsym = d;
        }
    for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b < P; ++b)
            res.Cov(a, b) = (res.dm(a, b) + res.dm(b, a)) / num_traits<T>::from_int(2);

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t g = 0; g < G; ++g) {
            const T d1 = res.dm(pidx(i, g), pidx(i, g));
            const T mi_g = res.m(i, g);
            res.Var(i, g) = d1;
            res.M2(i, g) = d1 + mi_g * mi_g;
            res.M3(i, g) = res.d2m(i, g) +
                           (num_traits<T>::from_int(1) + num_traits<T>::from_int(3) * mi_g) * d1 +
                           mi_g * mi_g * mi_g;
            const T mu3 = res.M3(i, g) - num_traits<T>::from_int(3) * mi_g * res.M2(i, g) +
                          num_traits<T>::from_int(2) * mi_g * mi_g * mi_g;
            const double var = num_traits<T>::to_double(res.Var(i, g));
            res.Skew(i, g) = var > 0.0 ? num_traits<T>::to_double(mu3) / std::pow(var, 1.5)
                                       : std::nan("");
        }

    return res;
}

/** pfqn_sens_mom with unit multiplicities and per-station totals. */
template <class T>
SensMomResult<T> pfqn_sens_mom(const Matrix<T>& L, const std::vector<int>& N,
                               const std::vector<T>& Z) {
    return pfqn_sens_mom(L, N, Z, std::vector<int>(), std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_MOM_H
