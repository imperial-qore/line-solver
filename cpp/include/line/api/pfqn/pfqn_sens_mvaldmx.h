/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_MVALDMX_H
#define LINE_API_PFQN_SENS_MVALDMX_H

/**
 * Exact queue-length variances and covariances of a mixed open/closed
 * product-form network with limited load dependence, the load-dependent and
 * mixed counterpart of pfqn_sens_mva.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_mvaldmx.m. Theorem 1,
 * equation (11), of Akyildiz and Strelen states that one further factor Q_jT
 * in a moment costs one derivative with respect to a parameter y_j scaling the
 * demands of the classes of T at station j; taking k = 2 and T = {s} gives
 *
 *   Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s) |_{y=1}
 *
 * which is evaluated here by forward-mode differentiation of the mixed
 * load-dependent MVA of Bruell-Balbo-Afshari, i.e. of the recursion
 * pfqn_mvaldmx implements. The differentiated equations are (13) for the
 * residence times, (15)-(17) for the conditional marginals, (18) for the
 * throughputs, (19) and (24)-(31) for the effective capacities (delegated to
 * pfqn_sens_ldmx_ec), (32) for the closed-class queue lengths and (33) for the
 * open-class ones.
 *
 * A demand-scaling parameter perturbs the whole network through the
 * closed-class throughputs, so the derivatives must be propagated for every
 * parameter and the cross-station covariances come out at no extra cost; they
 * are returned in QCovFull. That is unlike pfqn_sens_mva, whose cheaper
 * same-station recursion cannot reach them.
 *
 * Arithmetic. Field operations only, so the routine instantiates at
 * line::Rational. The one non-field constant is the max(eps, .) floor the
 * reference keeps on P(0) so that the base measures agree with pfqn_mvaldmx
 * entry by entry; it is a numerical guard, not part of the model, and the
 * derivative is the exact -sum of the derivatives, as in the reference.
 *
 * Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
 * Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
 * Communications 39(6):828-832, 1991.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_sens_ldmx_ec.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct SensMvaldmxResult {
    std::vector<T> XN;  ///< (R) throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> CN;       ///< (M x R) residence time (MATLAB field .R)

    std::vector<Matrix<T>> QCov;  ///< (M) matrices R x R, same-station covariance
    /// (M*R x M*R) Cov[n(i,r),n(j,s)] at row i*R+r, column j*R+s.
    Matrix<T> QCovFull;
    Matrix<T> QVar;         ///< (M x R)
    std::vector<T> QTotVar; ///< (M)
    T QCovAsym;             ///< residual before symmetrization
};

/**
 * @param lambda (R) arrival rates, zero on the closed classes
 * @param D      (M x R) service demands
 * @param N      (R) population, negative marks an open class (MATLAB uses Inf)
 * @param Z      (R) think times
 * @param mu     (M x >= sum of the closed populations) load-dependent rates
 */
template <class T>
SensMvaldmxResult<T> pfqn_sens_mvaldmx(const std::vector<T>& lambda, const Matrix<T>& D,
                                       const std::vector<int>& N, const std::vector<T>& Z,
                                       const Matrix<T>& mu) {
    const std::size_t M = D.rows();
    const std::size_t R = D.cols();
    if (lambda.size() != R || N.size() != R || Z.size() != R)
        throw InputError("pfqn_sens_mvaldmx: lambda, N, Z and D disagree on the class count");
    if (mu.rows() != M)
        throw InputError("pfqn_sens_mvaldmx: mu and D disagree on the station count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> openClasses, closedClasses;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) {
            openClasses.push_back(r);
        } else {
            if (lambda[r] != zero && N[r] > 0)
                throw InputError("pfqn_sens_mvaldmx: an arrival rate cannot be set on a closed class");
            closedClasses.push_back(r);
        }
    }
    const std::size_t Cn = closedClasses.size();
    if (Cn == 0)
        throw InputError(
            "pfqn_sens_mvaldmx: at least one closed class is required; use the open-class formulas "
            "directly otherwise");

    std::vector<int> Nc(Cn, 0);
    std::vector<T> Zc(Cn, zero);
    Matrix<T> Dc(M, Cn, zero);
    int NCtot = 0;
    for (std::size_t c = 0; c < Cn; ++c) {
        Nc[c] = N[closedClasses[c]];
        Zc[c] = Z[closedClasses[c]];
        NCtot += Nc[c];
        for (std::size_t i = 0; i < M; ++i) Dc(i, c) = D(i, closedClasses[c]);
    }
    if (static_cast<int>(mu.cols()) < NCtot)
        throw InputError(
            "pfqn_sens_mvaldmx: the load-dependent rates must be given up to the maximum closed "
            "population");

    // The reference appends one more saturation column so that the effective
    // capacities reach sum(N)+1, which the open-class equation (33) reads.
    Matrix<T> mux(M, mu.cols() + 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t k = 0; k < mu.cols(); ++k) mux(i, k) = mu(i, k);
        mux(i, mu.cols()) = mu(i, mu.cols() - 1);
    }
    const SensLdmxEcResult<T> ec = pfqn_sens_ldmx_ec(lambda, D, mux);

    // ---- parameter list: y(j,r) multiplies D(j,r) --------------------------
    const std::size_t P = M * R;
    const auto pidx = [&](std::size_t j, std::size_t r) { return j * R + r; };
    // eq. (21): only an open-class parameter at station i perturbs Lo(i).
    Matrix<T> dLo(M, P, zero);
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t r = 0; r < R; ++r)
            if (N[r] < 0) dLo(j, pidx(j, r)) = lambda[r] * D(j, r);

    // ---- population lattice over the closed classes, first class fastest ----
    const std::vector<std::size_t> prods = plane_sizes(Nc);
    const std::size_t NT = population_count(Nc);
    const std::size_t NL = static_cast<std::size_t>(NCtot) + 1;

    std::vector<Matrix<T>> Pc(NT, Matrix<T>(M, NL, zero));
    std::vector<std::vector<Matrix<T>>> dPc(NT, std::vector<Matrix<T>>(M, Matrix<T>(NL, P, zero)));
    Matrix<T> x(NT, Cn, zero);
    std::vector<Matrix<T>> dx(NT, Matrix<T>(Cn, P, zero));
    std::vector<Matrix<T>> w(NT, Matrix<T>(M, Cn, zero));
    std::vector<std::vector<Matrix<T>>> dw(NT, std::vector<Matrix<T>>(M, Matrix<T>(Cn, P, zero)));

    for (std::size_t i = 0; i < M; ++i) Pc[0](i, 0) = one;  // eq. (16)

    std::vector<T> dacc(P, zero);
    for (std::size_t k = 0; k < NT; ++k) {
        std::vector<int> nvec(Cn, 0);
        int nc = 0;
        for (std::size_t c = 0; c < Cn; ++c) {
            nvec[c] = static_cast<int>((k / prods[c]) % static_cast<std::size_t>(Nc[c] + 1));
            nc += nvec[c];
        }

        // ---- residence times, eq. (12) and its derivative eq. (13) ---------
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t c = 0; c < Cn; ++c) {
                if (nvec[c] <= 0) continue;
                const std::size_t kc = k - prods[c];
                const std::size_t cls = closedClasses[c];
                T acc = zero;
                for (std::size_t p = 0; p < P; ++p) dacc[p] = zero;
                for (int n = 1; n <= nc; ++n) {
                    const T Pprev = Pc[kc](i, static_cast<std::size_t>(n - 1));
                    const T nT = num_traits<T>::from_int(n);
                    acc += nT * ec.EC(i, static_cast<std::size_t>(n - 1)) * Pprev;
                    for (std::size_t p = 0; p < P; ++p)
                        dacc[p] += nT * (ec.dEC(i, static_cast<std::size_t>(n - 1)) * dLo(i, p) * Pprev +
                                         ec.EC(i, static_cast<std::size_t>(n - 1)) *
                                             dPc[kc][i](static_cast<std::size_t>(n - 1), p));
                }
                w[k](i, c) = Dc(i, c) * acc;
                for (std::size_t p = 0; p < P; ++p) {
                    T dwq = Dc(i, c) * dacc[p];
                    if (p == pidx(i, cls)) dwq += Dc(i, c) * acc;  // d(D y)/dy = D
                    dw[k][i](c, p) = dwq;
                }
            }
        }

        // ---- throughputs, eq. (18) -------------------------------------------
        for (std::size_t c = 0; c < Cn; ++c) {
            T sw = zero;
            for (std::size_t i = 0; i < M; ++i) sw += w[k](i, c);
            const T den = Zc[c] + sw;
            if (den == zero) {
                x(k, c) = zero;
                continue;
            }
            const T nT = num_traits<T>::from_int(nvec[c]);
            x(k, c) = nT / den;
            if (nvec[c] > 0) {
                const T den2 = den * den;
                for (std::size_t p = 0; p < P; ++p) {
                    T sdw = zero;
                    for (std::size_t i = 0; i < M; ++i) sdw += dw[k][i](c, p);
                    dx[k](c, p) = -nT / den2 * sdw;
                }
            }
        }

        // ---- conditional marginals, eq. (14)-(15) -----------------------------
        for (std::size_t i = 0; i < M; ++i) {
            for (int n = 1; n <= nc; ++n) {
                const std::size_t nu = static_cast<std::size_t>(n);
                for (std::size_t c = 0; c < Cn; ++c) {
                    if (nvec[c] <= 0) continue;
                    const std::size_t kc = k - prods[c];
                    const std::size_t cls = closedClasses[c];
                    const T Pprev = Pc[kc](i, nu - 1);
                    const T ECn = ec.EC(i, nu - 1);
                    Pc[k](i, nu) += Dc(i, c) * ECn * x(k, c) * Pprev;
                    for (std::size_t p = 0; p < P; ++p) {
                        T dt = Dc(i, c) * (ec.dEC(i, nu - 1) * dLo(i, p) * x(k, c) * Pprev +
                                           ECn * dx[k](c, p) * Pprev +
                                           ECn * x(k, c) * dPc[kc][i](nu - 1, p));
                        if (p == pidx(i, cls)) dt += Dc(i, c) * ECn * x(k, c) * Pprev;
                        dPc[k][i](nu, p) += dt;
                    }
                }
            }
            // eq. (17), with the reference's floor on the primal only
            T s1 = zero;
            for (int n = 1; n <= nc; ++n) s1 += Pc[k](i, static_cast<std::size_t>(n));
            const T epsT = num_traits<T>::from_double(2.220446049250313e-16);
            const T p0 = one - s1;
            Pc[k](i, 0) = p0 > epsT ? p0 : epsT;
            for (std::size_t p = 0; p < P; ++p) {
                T ds = zero;
                for (int n = 1; n <= nc; ++n) ds += dPc[k][i](static_cast<std::size_t>(n), p);
                dPc[k][i](0, p) = -ds;
            }
        }
    }

    // ---- measures at the full population -----------------------------------
    const std::size_t kN = NT - 1;
    SensMvaldmxResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    std::vector<Matrix<T>> dQN(P, Matrix<T>(M, R, zero));

    // closed classes, eq. (32)
    for (std::size_t c = 0; c < Cn; ++c) {
        const std::size_t cls = closedClasses[c];
        res.XN[cls] = x(kN, c);
        const std::size_t kc = Nc[c] > 0 ? kN - prods[c] : kN;
        for (std::size_t i = 0; i < M; ++i) {
            res.CN(i, cls) = w[kN](i, c);
            res.QN(i, cls) = res.XN[cls] * res.CN(i, cls);
            for (std::size_t p = 0; p < P; ++p)
                dQN[p](i, cls) = dx[kN](c, p) * w[kN](i, c) + x(kN, c) * dw[kN][i](c, p);
            T uacc = zero;
            for (int n = 1; n <= NCtot; ++n) {
                const std::size_t nu = static_cast<std::size_t>(n);
                uacc += Dc(i, c) * x(kN, c) * ec.Eprime(i, nu - 1) / ec.E(i, nu - 1) *
                        Pc[kc](i, nu - 1);
            }
            res.UN(i, cls) = uacc;
        }
    }

    // open classes, eq. (33)
    for (std::size_t oi = 0; oi < openClasses.size(); ++oi) {
        const std::size_t r = openClasses[oi];
        res.XN[r] = lambda[r];
        for (std::size_t i = 0; i < M; ++i) {
            T acc = zero;
            for (std::size_t p = 0; p < P; ++p) dacc[p] = zero;
            for (int n = 0; n <= NCtot; ++n) {
                const std::size_t nu = static_cast<std::size_t>(n);
                const T Pn = Pc[kN](i, nu);
                const T nT = num_traits<T>::from_int(n + 1);
                acc += nT * ec.EC(i, nu) * Pn;
                for (std::size_t p = 0; p < P; ++p)
                    dacc[p] += nT * (ec.dEC(i, nu) * dLo(i, p) * Pn +
                                     ec.EC(i, nu) * dPc[kN][i](nu, p));
            }
            res.QN(i, r) = lambda[r] * D(i, r) * acc;
            if (lambda[r] != zero) res.CN(i, r) = res.QN(i, r) / lambda[r];
            for (std::size_t p = 0; p < P; ++p) {
                T dq = lambda[r] * D(i, r) * dacc[p];
                if (p == pidx(i, r)) dq += lambda[r] * D(i, r) * acc;
                dQN[p](i, r) = dq;
            }
            T uacc = zero;
            for (int n = 0; n <= NCtot; ++n) {
                const std::size_t nu = static_cast<std::size_t>(n);
                uacc += lambda[r] * ec.Eprime(i, nu + 1) / ec.E(i, nu + 1) * Pc[kN](i, nu);
            }
            res.UN(i, r) = uacc;
        }
    }

    // ---- moments -------------------------------------------------------------
    res.QCovFull = Matrix<T>(M * R, M * R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t s = 0; s < R; ++s)
                    res.QCovFull(i * R + r, j * R + s) = dQN[pidx(j, s)](i, r);
    res.QCovAsym = zero;
    for (std::size_t a = 0; a < M * R; ++a)
        for (std::size_t b = 0; b < M * R; ++b) {
            const T d = num_abs(T(res.QCovFull(a, b) - res.QCovFull(b, a)));
            if (d > res.QCovAsym) res.QCovAsym = d;
        }
    for (std::size_t a = 0; a < M * R; ++a)
        for (std::size_t b = a + 1; b < M * R; ++b) {
            const T avg = (res.QCovFull(a, b) + res.QCovFull(b, a)) / num_traits<T>::from_int(2);
            res.QCovFull(a, b) = avg;
            res.QCovFull(b, a) = avg;
        }

    res.QCov.assign(M, Matrix<T>(R, R, zero));
    res.QVar = Matrix<T>(M, R, zero);
    res.QTotVar.assign(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        T tot = zero;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t s = 0; s < R; ++s) {
                res.QCov[i](r, s) = res.QCovFull(i * R + r, i * R + s);
                tot += res.QCov[i](r, s);
            }
            res.QVar(i, r) = res.QCov[i](r, r);
        }
        res.QTotVar[i] = tot;
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_MVALDMX_H
