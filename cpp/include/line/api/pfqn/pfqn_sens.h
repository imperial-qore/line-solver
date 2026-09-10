/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_H
#define LINE_API_PFQN_SENS_H

/**
 * Exact analytic derivatives of the mean performance measures {X,Q,U,R} of a
 * closed product-form (BCMP) network with respect to the demands L(i,r) and
 * the think times Z(r).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens.m. Two exact kernels are
 * dispatched transparently, exactly as the reference does:
 *
 *  - pfqn_sens_dmva : forward-mode differentiation of the exact
 *    Reiser-Lavenberg MVA recursion. Handles every model.
 *  - pfqn_sens_comom : the Class-Oriented Method of Moments specialization for
 *    the repairman model (one single-server queue plus a delay, every
 *    populated class having a strictly positive demand and think time).
 *
 * Both return the identical layout, so the choice is invisible to callers.
 *
 * Arithmetic. The derivatives are analytic, not finite differences, and every
 * step of both kernels is a field operation, so both instantiate at
 * line::Rational and return exact rationals. Where the reference is forced
 * through logarithms this port is not: MATLAB's sens_comom evaluates the
 * normalizing-constant ratios as exp(lgm(m,n-e_s) - lgm(m,n)) because
 * pfqn_comomrm returns only a log, whereas this port takes the ratio of the
 * ComomResult::G values directly. That is the same identity evaluated in the
 * field instead of through a transcendental round trip, so it is exact at
 * Rational and agrees with MATLAB to rounding at double. This is the same
 * substitution pfqn_mva.h already makes for the normalizing constant.
 *
 * Reference: Z. Liu and P. Nain, INRIA RR-1144, 1989; X.-R. Cao and D.-J. Ma,
 * Performance Evaluation 26:181-199, 1996; G. Casale, IEEE TSE 2011.
 */

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_sens_mva.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** One differentiation parameter. */
struct SensParam {
    char type;         ///< 'L' for a demand, 'Z' for a think time
    int station;       ///< 0-based station index, -1 for a 'Z' parameter
    std::size_t cls;   ///< 0-based class index
};

template <class T>
struct SensResult {
    std::vector<T> XN;  ///< (R) throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> CN;       ///< (M x R) residence time (MATLAB field .R)

    std::vector<SensParam> params;  ///< (P) the differentiation parameters

    Matrix<T> dX;               ///< (R x P) dX(r)/dparam(p)
    std::vector<Matrix<T>> dQ;  ///< (P) matrices M x R, dQ[p](i,r)
    std::vector<Matrix<T>> dU;  ///< (P) matrices M x R
    std::vector<Matrix<T>> dR;  ///< (P) matrices M x R

    /// (M*R x M*R) Cov[n(i,r),n(j,s)] at row i*R+r, column j*R+s.
    Matrix<T> QCov;
    Matrix<T> QVar;          ///< (M x R)
    std::vector<T> QTotVar;  ///< (M)
    T QCovAsym;              ///< residual of the moment recursion
};

namespace detail {

template <class T>
SensResult<T> sens_pack(std::size_t M, std::size_t R, std::size_t P) {
    const T zero = num_traits<T>::from_int(0);
    SensResult<T> s;
    s.XN.assign(R, zero);
    s.QN = Matrix<T>(M, R, zero);
    s.UN = Matrix<T>(M, R, zero);
    s.CN = Matrix<T>(M, R, zero);
    s.dX = Matrix<T>(R, P, zero);
    s.dQ.assign(P, Matrix<T>(M, R, zero));
    s.dU.assign(P, Matrix<T>(M, R, zero));
    s.dR.assign(P, Matrix<T>(M, R, zero));
    s.QCov = Matrix<T>(M * R, M * R, zero);
    s.QVar = Matrix<T>(M, R, zero);
    s.QTotVar.assign(M, zero);
    s.QCovAsym = zero;
    return s;
}

}  // namespace detail

/**
 * Forward-mode differentiation of the exact MVA recursion. Parameters are
 * ordered L(0,0), L(0,1), ..., L(M-1,R-1), then Z(0), ..., Z(R-1).
 */
template <class T>
SensResult<T> pfqn_sens_dmva(const Matrix<T>& L, const std::vector<int>& N,
                             const std::vector<T>& Z, const std::vector<int>& mi) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_sens: demand matrix and population vector disagree on the class count");
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_sens: think-time vector has the wrong length");
    if (!mi.empty() && mi.size() != M)
        throw InputError("pfqn_sens: multiplicity vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const std::size_t P = M * R + R;

    SensResult<T> res = detail::sens_pack<T>(M, R, P);
    res.params.resize(P);
    std::vector<std::vector<std::size_t>> pL(M, std::vector<std::size_t>(R, 0));
    std::vector<std::size_t> pZ(R, 0);
    {
        std::size_t p = 0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                res.params[p].type = 'L';
                res.params[p].station = static_cast<int>(i);
                res.params[p].cls = r;
                pL[i][r] = p;
                ++p;
            }
        for (std::size_t r = 0; r < R; ++r) {
            res.params[p].type = 'Z';
            res.params[p].station = -1;
            res.params[p].cls = r;
            pZ[r] = p;
            ++p;
        }
    }

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0 || R == 0) return res;

    const auto Zr = [&](std::size_t r) -> T { return Z.empty() ? zero : Z[r]; };
    const auto miT = [&](std::size_t i) -> T {
        return num_traits<T>::from_int(mi.empty() ? 1 : mi[i]);
    };

    const std::vector<std::size_t> radix = sens_lattice_radix(N);
    const std::size_t totpop = population_count(N);

    Matrix<T> Qtot(totpop, M, zero);
    std::vector<Matrix<T>> Qtotd(totpop, Matrix<T>(M, P, zero));

    std::vector<T> CNtotd(P, zero), Xd(P, zero), Qd(P, zero);
    Matrix<T> Cd(M, P, zero);

    for (std::size_t k = 1; k < totpop; ++k) {
        const std::vector<int> n = sens_lattice_decode(k, N, radix);
        for (std::size_t s = 0; s < R; ++s) {
            const std::size_t row = n[s] > 0 ? k - radix[s] : 0;
            T CNtot = zero;
            for (std::size_t p = 0; p < P; ++p) CNtotd[p] = zero;
            for (std::size_t i = 0; i < M; ++i) {
                const T base = miT(i) + Qtot(row, i);
                res.CN(i, s) = L(i, s) * base;
                for (std::size_t p = 0; p < P; ++p) Cd(i, p) = L(i, s) * Qtotd[row](i, p);
                Cd(i, pL[i][s]) += base;
                CNtot += res.CN(i, s);
                for (std::size_t p = 0; p < P; ++p) CNtotd[p] += Cd(i, p);
            }
            const T den = Zr(s) + CNtot;
            const T nsT = num_traits<T>::from_int(n[s]);
            if (den == zero) {
                res.XN[s] = zero;
                for (std::size_t p = 0; p < P; ++p) Xd[p] = zero;
            } else {
                const T den2 = den * den;
                res.XN[s] = nsT / den;
                for (std::size_t p = 0; p < P; ++p) Xd[p] = -nsT * CNtotd[p] / den2;
                Xd[pZ[s]] -= nsT / den2;
            }
            for (std::size_t p = 0; p < P; ++p) res.dX(s, p) = Xd[p];
            for (std::size_t i = 0; i < M; ++i) {
                res.QN(i, s) = res.XN[s] * res.CN(i, s);
                for (std::size_t p = 0; p < P; ++p) {
                    Qd[p] = Xd[p] * res.CN(i, s) + res.XN[s] * Cd(i, p);
                    res.dQ[p](i, s) = Qd[p];
                    res.dR[p](i, s) = Cd(i, p);
                    Qtotd[k](i, p) += Qd[p];
                }
                Qtot(k, i) += res.QN(i, s);
            }
        }
    }

    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) {
            res.UN(i, r) = res.XN[r] * L(i, r);
            for (std::size_t p = 0; p < P; ++p) res.dU[p](i, r) = res.dX(r, p) * L(i, r);
            res.dU[pL[i][r]](i, r) += res.XN[r];
        }
    }
    return res;
}

/**
 * CoMoM-backed kernel for the repairman model (M = 1). Parameters are ordered
 * L(0,0), ..., L(0,R-1), then Z(0), ..., Z(R-1), matching pfqn_sens_dmva at
 * M = 1.
 */
template <class T>
SensResult<T> pfqn_sens_comom(const Matrix<T>& L, const std::vector<int>& N,
                              const std::vector<T>& Z) {
    const std::size_t R = N.size();
    if (L.rows() != 1) throw InputError("pfqn_sens_comom: the kernel takes a single station");
    if (L.cols() != R || Z.size() != R)
        throw InputError("pfqn_sens_comom: L, N and Z disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const std::size_t P = 2 * R;

    SensResult<T> res = detail::sens_pack<T>(1, R, P);
    res.params.resize(P);
    for (std::size_t r = 0; r < R; ++r) {
        res.params[r].type = 'L';
        res.params[r].station = 0;
        res.params[r].cls = r;
        res.params[R + r].type = 'Z';
        res.params[R + r].station = -1;
        res.params[R + r].cls = r;
    }

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens_comom: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || R == 0) return res;

    // memoized-constant class-stripping rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::map<std::pair<int, std::vector<int>>, T> cache;
    const T one = num_traits<T>::from_int(1);
    const auto Gm = [&](int m, const std::vector<int>& n) -> T {
        std::vector<int> nz;
        std::vector<std::size_t> keep;
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0) {
                nz.push_back(n[r]);
                keep.push_back(r);
            }
        if (nz.empty()) return one;
        const std::pair<int, std::vector<int>> key(m, n);
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        Matrix<T> Ls(1, keep.size(), zero);
        Matrix<T> Zs(1, keep.size(), zero);
        for (std::size_t j = 0; j < keep.size(); ++j) {
            Ls(0, j) = L(0, keep[j]);
            Zs(0, j) = Z[keep[j]];
        }
        const T g = pfqn_comomrm(Ls, nz, Zs, m).G;
        cache.emplace(key, g);
        return g;
    };
    const auto minus = [&](const std::vector<int>& n, std::size_t s) {
        std::vector<int> m = n;
        m[s] -= 1;
        return m;
    };
    // Mean class-s queue at population n of the m = 1 model.
    const auto qmean = [&](const std::vector<int>& n, std::size_t s) -> T {
        if (n[s] < 1) return zero;
        return L(0, s) * Gm(2, minus(n, s)) / Gm(1, n);
    };
    // Class-s throughput of the m-replica model.
    const auto xput = [&](int m, const std::vector<int>& n, std::size_t s) -> T {
        if (n[s] < 1) return zero;
        return Gm(m, minus(n, s)) / Gm(m, n);
    };
    // Class-s queue at one replica of the doubled station.
    const auto qplus = [&](const std::vector<int>& n, std::size_t s) -> T {
        if (n[s] < 1) return zero;
        return L(0, s) * Gm(3, minus(n, s)) / Gm(2, n);
    };

    for (std::size_t r = 0; r < R; ++r)
        if (N[r] >= 1) res.XN[r] = xput(1, N, r);
    for (std::size_t s = 0; s < R; ++s) res.QN(0, s) = qmean(N, s);
    for (std::size_t r = 0; r < R; ++r) {
        res.UN(0, r) = res.XN[r] * L(0, r);
        if (res.XN[r] != zero) res.CN(0, r) = res.QN(0, r) / res.XN[r];
    }

    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 1) continue;  // empty class: X = Q = 0 and every derivative vanishes
        const std::vector<int> Nr = minus(N, r);
        const T Xr = res.XN[r];
        const T Qr = res.QN(0, r);
        for (std::size_t s = 0; s < R; ++s) {
            const T dlt = num_traits<T>::from_int(r == s ? 1 : 0);
            // L(0,s): D_s dQ_r/dD_s = Cov[n_r,n_s]
            const T Vrs = Qr * (dlt + num_traits<T>::from_int(2) * qplus(Nr, s) - res.QN(0, s));
            const T dQ_L = Vrs / L(0, s);
            const T dX_L = Xr * (qmean(Nr, s) - res.QN(0, s)) / L(0, s);
            const T dU_L = dX_L * L(0, r) + Xr * dlt;
            res.dQ[s](0, r) = dQ_L;
            res.dX(r, s) = dX_L;
            res.dU[s](0, r) = dU_L;
            if (Xr != zero) res.dR[s](0, r) = (dQ_L * Xr - Qr * dX_L) / (Xr * Xr);

            // Z(s): d log G_m(n)/dZ_s = G_m(n - e_s)/G_m(n), exact for any Z_s >= 0
            const std::size_t p = R + s;
            const T dX_Z = Xr * (xput(1, Nr, s) - res.XN[s]);
            const T dQ_Z = Qr * (xput(2, Nr, s) - res.XN[s]);
            res.dQ[p](0, r) = dQ_Z;
            res.dX(r, p) = dX_Z;
            res.dU[p](0, r) = dX_Z * L(0, r);
            if (Xr != zero) res.dR[p](0, r) = (dQ_Z * Xr - Qr * dX_Z) / (Xr * Xr);
        }
    }
    return res;
}

/**
 * @param L  (M x R) service demands
 * @param N  (R) population per class
 * @param Z  (R) think times, empty for none
 * @param mi (M) station multiplicities, empty for all ones
 */
template <class T>
SensResult<T> pfqn_sens(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                        const std::vector<int>& mi) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);

    // Dispatch: the repairman model goes to CoMoM, everything else to the
    // differentiated MVA, exactly as the reference decides.
    const T fineTol = num_traits<T>::from_double(1e-8);  // GlobalConstants.FineTol
    bool anyPopulated = false, useComom = (M == 1);
    if (!mi.empty())
        for (std::size_t i = 0; i < M; ++i)
            if (mi[i] != 1) useComom = false;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] <= 0) continue;
        anyPopulated = true;
        if (Zv[r] <= fineTol) useComom = false;
        if (M == 1 && L(0, r) <= fineTol) useComom = false;
    }
    useComom = useComom && anyPopulated;

    SensResult<T> res = useComom ? pfqn_sens_comom(L, N, Zv) : pfqn_sens_dmva(L, N, Zv, mi);

    // queue-length covariance identity rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const SensMvaResult<T> mom = pfqn_sens_mva(L, N, Zv, mi);
    std::vector<std::vector<std::size_t>> pL(M, std::vector<std::size_t>(R, 0));
    std::vector<std::vector<bool>> haveL(M, std::vector<bool>(R, false));
    for (std::size_t p = 0; p < res.params.size(); ++p)
        if (res.params[p].type == 'L') {
            pL[static_cast<std::size_t>(res.params[p].station)][res.params[p].cls] = p;
            haveL[static_cast<std::size_t>(res.params[p].station)][res.params[p].cls] = true;
        }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t s = 0; s < R; ++s) {
                    T v = zero;
                    if (i == j) {
                        v = mom.QCov[i](r, s);
                    } else if (haveL[j][s]) {
                        v = L(j, s) * res.dQ[pL[j][s]](i, r);
                    }
                    res.QCov(i * R + r, j * R + s) = v;
                }
    res.QVar = mom.QVar;
    res.QTotVar = mom.QTotVar;
    res.QCovAsym = mom.QCovAsym;
    return res;
}

/** pfqn_sens with unit multiplicities. */
template <class T>
SensResult<T> pfqn_sens(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z) {
    return pfqn_sens(L, N, Z, std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_H
