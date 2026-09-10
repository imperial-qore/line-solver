/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_QSA_H
#define LINE_API_PFQN_QSA_H

/**
 * Queue-Shift Approximation (QSA) for closed product-form networks.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_qsa.m. Schweitzer, Serazzi and
 * Broglia, "A Queue-Shift Approximation Technique for Product-Form Queueing
 * Networks", Tools'98, LNCS 1469, pp. 267-279.
 *
 * Where Linearizer extrapolates the fractional deviation D_rit, QSA
 * extrapolates the ABSOLUTE shift of the aggregate queue length,
 *   Y_ri(K) = 1 + Q_i(K - e_r) - Q_i(K)      i in QC
 * so the unknowns are one per station rather than one per station-class. The
 * core equation (13a),
 *   Q_i(K) = sum_r K_r L_ri [Q_i(K) + Y_ri(K)] / C_r(K)
 * is imposed at K, at every K - e_s and, in the three-level variant of eq. (16),
 * at every K - e_s - e_t through the affine extrapolation of eq. (15).
 *
 * The quintuple (16) is solved as ONE system by damped Newton, as Sect. 4 of
 * the paper prescribes. The decomposed successive substitution that works for
 * Linearizer must NOT be used here: it drifts to the degenerate root in which
 * the bottleneck absorbs the whole population, and does so even when seeded at
 * the exact solution, because the instability is a positive real eigenvalue
 * rather than an oscillation that under-relaxation could damp.
 *
 * Iterates to a residual tolerance, so exact arithmetic buys nothing: the
 * static_assert records that.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/// Shift matrix (M x R) of (16d)-(16e), or (15) when both s and t are set.
template <class T>
Matrix<T> qsa_shift(const Matrix<T>& q, const Matrix<T>& pops, const std::vector<long>& sIdx,
                    const Matrix<long>& pIdx, long s, long t, int levels) {
    const std::size_t M = q.rows(), R = pops.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> Y(M, R, zero);
    if (s < 0) {
        for (std::size_t r = 0; r < R; ++r)
            if (sIdx[r] >= 0)
                for (std::size_t i = 0; i < M; ++i)
                    Y(i, r) = one + q(i, static_cast<std::size_t>(sIdx[r])) - q(i, 0);
    } else if (t < 0) {
        if (levels < 3) return qsa_shift(q, pops, sIdx, pIdx, -1, -1, levels);  // (14)
        const std::size_t ps = static_cast<std::size_t>(sIdx[static_cast<std::size_t>(s)]);
        for (std::size_t r = 0; r < R; ++r) {
            const long pr = pIdx(static_cast<std::size_t>(s), r);
            if (pr >= 0 && pops(ps, r) >= one)
                for (std::size_t i = 0; i < M; ++i)
                    Y(i, r) = one + q(i, static_cast<std::size_t>(pr)) - q(i, ps);
        }
    } else {
        const Matrix<T> Ys = qsa_shift(q, pops, sIdx, pIdx, s, -1, levels);
        const Matrix<T> Yt = qsa_shift(q, pops, sIdx, pIdx, t, -1, levels);
        const Matrix<T> Y0 = qsa_shift(q, pops, sIdx, pIdx, -1, -1, levels);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) Y(i, r) = Ys(i, r) + Yt(i, r) - Y0(i, r);
    }
    return Y;
}

/// Decode a population index into the removed classes.
inline void qsa_which(std::size_t p, const std::vector<long>& sIdx, const Matrix<long>& pIdx,
                      long& s, long& t) {
    s = -1;
    t = -1;
    if (p == 0) return;
    for (std::size_t k = 0; k < sIdx.size(); ++k)
        if (sIdx[k] == static_cast<long>(p)) {
            s = static_cast<long>(k);
            return;
        }
    for (std::size_t a = 0; a < pIdx.rows(); ++a)
        for (std::size_t b = 0; b < pIdx.cols(); ++b)
            if (pIdx(a, b) == static_cast<long>(p)) {
                s = static_cast<long>(a);
                t = static_cast<long>(b);
                return;
            }
}

/**
 * Residual of (13) imposed simultaneously at every population of (16). adm
 * carries the side conditions of Remark 2: non-negative queue lengths and
 * positive cycle times.
 */
template <class T>
std::vector<T> qsa_resid(const std::vector<T>& x, const Matrix<T>& L, const std::vector<T>& Z,
                         const Matrix<T>& pops, const std::vector<long>& sIdx,
                         const Matrix<long>& pIdx, const std::vector<std::size_t>& qc,
                         const std::vector<T>& Ldc, int levels, bool& adm) {
    const std::size_t M = L.rows(), R = L.cols(), nP = pops.rows(), mq = qc.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    adm = true;
    Matrix<T> q(M, nP, zero);
    for (std::size_t a = 0; a < mq; ++a)
        for (std::size_t p = 0; p < nP; ++p) {
            q(qc[a], p) = x[a * nP + p];
            if (x[a * nP + p] < zero) adm = false;
        }
    std::vector<T> F(mq * nP, zero);
    for (std::size_t p = 0; p < nP; ++p) {
        long s = 0, t = 0;
        qsa_which(p, sIdx, pIdx, s, t);
        const Matrix<T> Y = qsa_shift(q, pops, sIdx, pIdx, s, t, levels);
        std::vector<T> acc(mq, zero);
        for (std::size_t r = 0; r < R; ++r) {
            if (pops(p, r) < one) continue;
            T c = Z.empty() ? Ldc[r] : T(Z[r] + Ldc[r]);
            for (std::size_t a = 0; a < mq; ++a) c += L(qc[a], r) * (q(qc[a], p) + Y(qc[a], r));
            if (!(c > zero) || !std::isfinite(num_traits<T>::to_double(c))) {
                adm = false;
                c = num_traits<T>::from_double(1e-300);
            }
            const T xr = pops(p, r) / c;
            for (std::size_t a = 0; a < mq; ++a)
                acc[a] += xr * L(qc[a], r) * (q(qc[a], p) + Y(qc[a], r));
        }
        for (std::size_t a = 0; a < mq; ++a) F[a * nP + p] = q(qc[a], p) - acc[a];
    }
    return F;
}

/**
 * Aggregate Bard-Schweitzer queue lengths at population n, with the
 * delay-centre demands folded into the think time.
 */
template <class T>
std::vector<T> qsa_aggbs(const Matrix<T>& L, const std::vector<T>& n, const std::vector<T>& Z,
                         const std::vector<bool>& isQC) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> q(M, zero);
    std::vector<T> nn(R, zero);
    bool empty = true;
    for (std::size_t r = 0; r < R; ++r) {
        nn[r] = (n[r] > zero) ? n[r] : zero;
        if (nn[r] > zero) empty = false;
    }
    if (empty) return q;

    std::vector<T> Zeff(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        Zeff[r] = Z.empty() ? zero : Z[r];
        for (std::size_t i = 0; i < M; ++i)
            if (!isQC[i]) Zeff[r] += L(i, r);
    }
    std::size_t mq = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (isQC[i]) ++mq;

    std::vector<T> X(R, zero);
    if (mq > 0) {
        Matrix<T> Lq(mq, R, zero);
        std::size_t a = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (isQC[i]) {
                for (std::size_t r = 0; r < R; ++r) Lq(a, r) = L(i, r);
                ++a;
            }
        const AmvaResult<T> bs = pfqn_bs(Lq, nn, Zeff, std::vector<AmvaSched>());
        a = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (isQC[i]) {
                T s = zero;
                for (std::size_t r = 0; r < R; ++r) s += bs.QN(a, r);
                q[i] = s;
                ++a;
            }
        X = bs.XN;
    } else {
        for (std::size_t r = 0; r < R; ++r)
            if (nn[r] >= one && Zeff[r] > zero) X[r] = nn[r] / Zeff[r];
    }
    for (std::size_t i = 0; i < M; ++i)
        if (!isQC[i]) {
            T s = zero;
            for (std::size_t r = 0; r < R; ++r) s += X[r] * L(i, r);
            q[i] = s;
        }
    return q;
}

}  // namespace detail

/**
 * @param L    (M x R) demands
 * @param N    (R) populations
 * @param Z    (R) think times, empty for none
 * @param type (M) per-station scheduling; AmvaSched::INF marks a delay centre,
 *             whose demand enters the cycle time without a queueing term (the
 *             paper's DC set). Empty means every station queues.
 * @param tol  residual tolerance of the Newton iteration
 * @param maxiter maximum Newton iterations
 * @param levels 2 for the two-level QSA of eq. (14), 3 for eq. (16)
 */
template <class T>
AmvaResult<T> pfqn_qsa(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                       const std::vector<AmvaSched>& type, double tol = 1e-10,
                       std::size_t maxiter = 100, int levels = 3) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_qsa requires transcendental arithmetic: the Newton iteration stops on a "
                  "residual tolerance, so its answer is a fixed point only to within tol");
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_qsa: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_qsa: Z has the wrong length");
    if (!type.empty() && type.size() != M) throw InputError("pfqn_qsa: type has the wrong length");
    if (levels != 2 && levels != 3) throw InputError("pfqn_qsa: levels must be 2 or 3");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<bool> isQC(M, true);
    for (std::size_t i = 0; i < M && !type.empty(); ++i) isQC[i] = type[i] != AmvaSched::INF;

    AmvaResult<T> out;
    out.XN.assign(R, zero);
    out.QN = Matrix<T>(M, R, zero);
    out.UN = Matrix<T>(M, R, zero);
    out.RN = Matrix<T>(M, R, zero);

    bool emptyDemands = true;
    for (std::size_t i = 0; i < M && emptyDemands; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) != zero) {
                emptyDemands = false;
                break;
            }
    bool emptyPop = true;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > zero) emptyPop = false;
    if (M == 0 || emptyDemands || emptyPop) {
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] > zero && !Z.empty() && Z[r] > zero) out.XN[r] = N[r] / Z[r];
            for (std::size_t i = 0; i < M; ++i) out.UN(i, r) = out.XN[r] * L(i, r);
        }
        return out;
    }

    // Populations touched by (16): K, every K - e_s, every K - e_s - e_t.
    std::vector<std::vector<T>> popList;
    popList.push_back(N);
    std::vector<long> sIdx(R, -1);
    Matrix<long> pIdx(R, R, -1);
    for (std::size_t s = 0; s < R; ++s) {
        std::vector<T> n = N;
        n[s] -= one;
        bool ok = true;
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] < zero) ok = false;
        if (ok) {
            popList.push_back(n);
            sIdx[s] = static_cast<long>(popList.size()) - 1;
        }
    }
    if (levels >= 3) {
        for (std::size_t s = 0; s < R; ++s)
            for (std::size_t t = s; t < R; ++t) {
                std::vector<T> n = N;
                n[s] -= one;
                n[t] -= one;
                bool ok = true;
                for (std::size_t r = 0; r < R; ++r)
                    if (n[r] < zero) ok = false;
                if (ok) {
                    popList.push_back(n);
                    pIdx(s, t) = static_cast<long>(popList.size()) - 1;
                    pIdx(t, s) = pIdx(s, t);
                }
            }
    }
    const std::size_t nP = popList.size();
    Matrix<T> pops(nP, R, zero);
    for (std::size_t p = 0; p < nP; ++p)
        for (std::size_t r = 0; r < R; ++r) pops(p, r) = popList[p][r];

    std::vector<std::size_t> qc;
    for (std::size_t i = 0; i < M; ++i)
        if (isQC[i]) qc.push_back(i);
    const std::size_t mq = qc.size();
    std::vector<T> Ldc(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i)
            if (!isQC[i]) Ldc[r] += L(i, r);

    // Bard-Schweitzer at every population supplies the Newton starting point.
    Matrix<T> q(M, nP, zero);
    for (std::size_t p = 0; p < nP; ++p) {
        const std::vector<T> qp = detail::qsa_aggbs(L, popList[p], Z, isQC);
        for (std::size_t i = 0; i < M; ++i) q(i, p) = qp[i];
    }

    const std::size_t nUnk = mq * nP;
    std::vector<T> x(nUnk, zero);
    for (std::size_t a = 0; a < mq; ++a)
        for (std::size_t p = 0; p < nP; ++p) x[a * nP + p] = q(qc[a], p);

    bool adm = true;
    std::vector<T> F = detail::qsa_resid(x, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, adm);
    auto nrm = [](const std::vector<T>& v) {
        double s = 0.0;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const double d = num_traits<T>::to_double(v[i]);
            s += d * d;
        }
        return std::sqrt(s);
    };
    double fnrm = nrm(F);
    for (std::size_t it = 1; it <= maxiter; ++it) {
        if (fnrm < tol) break;
        out.iterations = it;
        Matrix<T> J(nUnk, nUnk, zero);
        for (std::size_t col = 0; col < nUnk; ++col) {
            const double xc = num_traits<T>::to_double(x[col]);
            const double hd = 1e-7 * std::max(1.0, std::fabs(xc));
            const T h = num_traits<T>::from_double(hd);
            std::vector<T> xp = x;
            xp[col] += h;
            bool dummy = true;
            const std::vector<T> Fp =
                detail::qsa_resid(xp, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, dummy);
            for (std::size_t row = 0; row < nUnk; ++row) J(row, col) = (Fp[row] - F[row]) / h;
        }
        std::vector<T> rhs(nUnk, zero);
        for (std::size_t row = 0; row < nUnk; ++row) rhs[row] = zero - F[row];
        std::vector<T> step;
        try {
            step = solve(J, rhs);
        } catch (const std::exception&) {
            break;   // singular Jacobian: keep the best iterate rather than guessing
        }
        bool accepted = false;
        double lambda = 1.0;
        for (int ls = 0; ls < 40; ++ls) {
            std::vector<T> xn(nUnk, zero);
            bool finite = true;
            const T lam = num_traits<T>::from_double(lambda);
            for (std::size_t row = 0; row < nUnk; ++row) {
                xn[row] = x[row] + lam * step[row];
                if (!std::isfinite(num_traits<T>::to_double(xn[row]))) finite = false;
            }
            if (finite) {
                bool admn = true;
                const std::vector<T> Fn =
                    detail::qsa_resid(xn, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, admn);
                const double nn = nrm(Fn);
                if (admn && nn < fnrm) {
                    x = xn;
                    F = Fn;
                    fnrm = nn;
                    adm = admn;
                    accepted = true;
                    break;
                }
            }
            lambda /= 2.0;
        }
        if (!accepted) break;
    }
    out.converged = adm && fnrm < tol;

    // Disaggregate (13) at K into the per-class measures
    for (std::size_t a = 0; a < mq; ++a)
        for (std::size_t p = 0; p < nP; ++p) q(qc[a], p) = x[a * nP + p];
    const Matrix<T> Y0 = detail::qsa_shift(q, pops, sIdx, pIdx, -1, -1, levels);
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < one) continue;
        T sumW = zero;
        for (std::size_t i = 0; i < M; ++i) {
            out.RN(i, r) = isQC[i] ? T(L(i, r) * (q(i, 0) + Y0(i, r))) : L(i, r);
            sumW += out.RN(i, r);
        }
        const T denom = (Z.empty() ? zero : Z[r]) + sumW;
        out.XN[r] = N[r] / denom;
        for (std::size_t i = 0; i < M; ++i) {
            out.QN(i, r) = out.XN[r] * out.RN(i, r);
            out.UN(i, r) = out.XN[r] * L(i, r);
        }
    }
    return out;
}

template <class T>
AmvaResult<T> pfqn_qsa(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_qsa(L, N, Z, std::vector<AmvaSched>());
}

template <class T>
AmvaResult<T> pfqn_qsa(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_qsa(L, N, std::vector<T>(), std::vector<AmvaSched>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_QSA_H
