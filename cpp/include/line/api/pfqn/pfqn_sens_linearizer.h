/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SENS_LINEARIZER_H
#define LINE_API_PFQN_SENS_LINEARIZER_H

/**
 * Approximate moments E[Q_i], Var[Q_i], Cov[Q_i,Q_j], E[Q_i^2] and E[Q_i^3] of
 * the per-station total queue lengths of a closed product-form network, by the
 * LINEARIZER-2 / LINEARIZER-3 algorithms of Strelen (Performance Evaluation
 * 11:127-142, 1990, Section 5).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_sens_linearizer.m. The exact
 * moment analysis of pfqn_sens_mom walks the whole population lattice and so
 * costs O(prod(N+1)); the Linearizer replaces that lattice by a fixed point
 * over R+1 populations, and the reference observes that the same trick applies
 * to the derivatives: differentiate the Linearizer equations, append them to
 * the originals and iterate everything together. Carrying the first derivative
 * is LINEARIZER-2, carrying the second as well is LINEARIZER-3; both are done
 * here, so the third moment is available.
 *
 * CORE (equations (5.1)-(5.2)) estimates the queue lengths one job down by
 *   v_i(l) = m_i(l)/n(l),  m_i^(n-e_l')(l) = (n - e_l')_l (v_i(l) + delta_i(l',l))
 * and substitutes them into the exact MVA equations; delta comes from (5.3),
 * delta_i(l',l) = v_i^(N-e_l')(l) - v_i^(N)(l), and is held fixed across
 * populations by the heuristic (5.4). Differentiating (5.1)-(5.3) gives
 * (5.5)-(5.8), carried alongside.
 *
 * Accuracy. The reference reports, over 51 networks, relative errors below
 * 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3]. This is an approximation
 * and is expected to disagree with pfqn_sens_mom by about that much. CovAsym
 * is a genuine error indicator here, not a roundoff residual: the product form
 * makes Cov symmetric but the Linearizer fixed point does not enforce it.
 *
 * Arithmetic. static_assert(has_transcendental) -- unlike every other member
 * of the sensitivity family, this routine reaches a fixed point only to within
 * a stopping tolerance, so its output is a function of the termination test
 * rather than of the model alone and exact arithmetic buys nothing. Same
 * rationale as pfqn_linearizer.h, on top of which this is built.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

template <class T>
struct SensLinearizerResult {
    std::vector<T> XN;  ///< (R) throughput
    Matrix<T> QN;       ///< (M x R) queue length per class
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> WN;       ///< (M x R) residence time

    std::vector<T> m;       ///< (M) E[Q_i], the total queue at station i
    Matrix<T> dm;           ///< (M x M) dm(i,h) = x_h dm_i/dx_h
    std::vector<T> d2m;     ///< (M) x_i^2 d2m_i/dx_i^2
    Matrix<T> Cov;          ///< (M x M) symmetrized dm
    std::vector<T> Var;     ///< (M)
    std::vector<T> M2;      ///< (M) E[Q^2]
    std::vector<T> M3;      ///< (M) E[Q^3]
    std::vector<double> Skew;  ///< (M) skewness, NaN where the variance vanishes
    T CovAsym;              ///< raw asymmetry of Cov, an error indicator here
    unsigned iter;          ///< total CORE iterations performed
};

namespace detail {

/** (5.1) and (5.5): v_i(l) = m_i(l)/n(l), likewise for the derivatives. */
template <class T>
void sens_lin_fractions(const Matrix<T>& m, const std::vector<Matrix<T>>& dm,
                        const std::vector<Matrix<T>>& d2m, const std::vector<int>& n,
                        Matrix<T>& v, std::vector<Matrix<T>>& dv, std::vector<Matrix<T>>& d2v) {
    const std::size_t M = m.rows(), R = m.cols(), MM = dm.size();
    const T zero = num_traits<T>::from_int(0);
    v = Matrix<T>(M, R, zero);
    dv.assign(MM, Matrix<T>(M, R, zero));
    d2v.assign(MM, Matrix<T>(M, R, zero));
    for (std::size_t l = 0; l < R; ++l) {
        if (n[l] <= 0) continue;
        const T nT = num_traits<T>::from_int(n[l]);
        for (std::size_t i = 0; i < M; ++i) {
            v(i, l) = m(i, l) / nT;
            for (std::size_t h = 0; h < MM; ++h) {
                dv[h](i, l) = dm[h](i, l) / nT;
                d2v[h](i, l) = d2m[h](i, l) / nT;
            }
        }
    }
}

/**
 * CORE-2 of the reference extended to second derivatives: iterates (5.1),
 * (5.2), the MVA equations (1.4)-(1.5) and their derivatives (5.5)-(5.6),
 * (3.4) until the mean queue lengths and the variances both stop moving.
 */
template <class T>
unsigned sens_lin_core2(const Matrix<T>& L, const std::vector<T>& Z, const std::vector<int>& n,
                        const std::vector<Matrix<T>>& delta,
                        const std::vector<std::vector<Matrix<T>>>& ddelta,
                        const std::vector<std::vector<Matrix<T>>>& d2delta, Matrix<T>& m,
                        std::vector<Matrix<T>>& dm, std::vector<Matrix<T>>& d2m,
                        std::vector<T>& lam, Matrix<T>& w, const T* tol, unsigned maxiter) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    int nc = 0;
    for (int v : n) nc += v;
    const T tolm = tol ? *tol
                       : num_traits<T>::from_int(1) /
                             num_traits<T>::from_int(4000 + 16 * static_cast<long>(nc));
    const T tolv = num_traits<T>::from_double(1e-3);

    lam.assign(R, zero);
    w = Matrix<T>(M, R, zero);
    std::vector<T> varprev(M, zero);
    unsigned it = 0;
    for (unsigned iter = 1; iter <= maxiter; ++iter) {
        it = iter;
        const Matrix<T> mprev = m;

        // ---- (5.1)-(5.2) and (5.5)-(5.6) ---------------------------------
        Matrix<T> v;
        std::vector<Matrix<T>> dv, d2v;
        sens_lin_fractions(m, dm, d2m, n, v, dv, d2v);
        Matrix<T> mtot(M, R, zero);
        std::vector<Matrix<T>> dmtot(M, Matrix<T>(M, R, zero));
        std::vector<Matrix<T>> d2mtot(M, Matrix<T>(M, R, zero));
        for (std::size_t l = 0; l < R; ++l) {
            if (n[l] <= 0) continue;
            for (std::size_t i = 0; i < M; ++i) {
                T acc = zero;
                std::vector<T> dacc(M, zero), d2acc(M, zero);
                for (std::size_t l2 = 0; l2 < R; ++l2) {
                    const int cnt = n[l2] - (l2 == l ? 1 : 0);
                    if (cnt <= 0) continue;
                    const T cT = num_traits<T>::from_int(cnt);
                    acc += cT * (v(i, l2) + delta[i](l, l2));
                    for (std::size_t h = 0; h < M; ++h) {
                        dacc[h] += cT * (dv[h](i, l2) + ddelta[i][h](l, l2));
                        d2acc[h] += cT * (d2v[h](i, l2) + d2delta[i][h](l, l2));
                    }
                }
                mtot(i, l) = acc;
                for (std::size_t h = 0; h < M; ++h) {
                    dmtot[h](i, l) = dacc[h];
                    d2mtot[h](i, l) = d2acc[h];
                }
            }
        }

        // ---- MVA (1.4)-(1.5) and its derivatives (3.4) --------------------
        w = Matrix<T>(M, R, zero);
        std::vector<Matrix<T>> dw(M, Matrix<T>(M, R, zero)), d2w(M, Matrix<T>(M, R, zero));
        const T one = num_traits<T>::from_int(1);
        for (std::size_t l = 0; l < R; ++l) {
            if (n[l] <= 0) continue;
            for (std::size_t i = 0; i < M; ++i) {
                const T A = one + mtot(i, l);
                w(i, l) = L(i, l) * A;
                for (std::size_t h = 0; h < M; ++h) {
                    const T dA = dmtot[h](i, l);
                    const T d2A = d2mtot[h](i, l);
                    if (i == h) {
                        dw[h](i, l) = L(i, l) * (A + dA);
                        d2w[h](i, l) = L(i, l) * (num_traits<T>::from_int(2) * dA + d2A);
                    } else {
                        dw[h](i, l) = L(i, l) * dA;
                        d2w[h](i, l) = L(i, l) * d2A;
                    }
                }
            }
        }
        lam.assign(R, zero);
        Matrix<T> dlam(R, M, zero), d2lam(R, M, zero);
        for (std::size_t l = 0; l < R; ++l) {
            if (n[l] <= 0) continue;
            T sw = zero;
            for (std::size_t i = 0; i < M; ++i) sw += w(i, l);
            const T den = (Z.empty() ? zero : Z[l]) + sw;
            const T nT = num_traits<T>::from_int(n[l]);
            lam[l] = nT / den;
            const T den2 = den * den;
            const T den3 = den2 * den;
            for (std::size_t h = 0; h < M; ++h) {
                T dden = zero, d2den = zero;
                for (std::size_t i = 0; i < M; ++i) {
                    dden += dw[h](i, l);
                    d2den += d2w[h](i, l);
                }
                dlam(l, h) = -nT * dden / den2;
                d2lam(l, h) = -nT * d2den / den2 + num_traits<T>::from_int(2) * nT * dden * dden / den3;
            }
        }
        m = Matrix<T>(M, R, zero);
        dm.assign(M, Matrix<T>(M, R, zero));
        d2m.assign(M, Matrix<T>(M, R, zero));
        for (std::size_t l = 0; l < R; ++l) {
            if (n[l] <= 0) continue;
            for (std::size_t i = 0; i < M; ++i) {
                m(i, l) = lam[l] * w(i, l);
                for (std::size_t h = 0; h < M; ++h) {
                    dm[h](i, l) = dlam(l, h) * w(i, l) + lam[l] * dw[h](i, l);
                    d2m[h](i, l) = d2lam(l, h) * w(i, l) +
                                   num_traits<T>::from_int(2) * dlam(l, h) * dw[h](i, l) +
                                   lam[l] * d2w[h](i, l);
                }
            }
        }

        // ---- termination test of the reference ----------------------------
        T dev = zero;
        for (std::size_t l = 0; l < R; ++l) {
            if (n[l] <= 0) continue;
            const T nT = num_traits<T>::from_int(n[l]);
            for (std::size_t i = 0; i < M; ++i) {
                const T d = num_abs(T(m(i, l) - mprev(i, l))) / nT;
                if (d > dev) dev = d;
            }
        }
        std::vector<T> varnow(M, zero);
        T sv = zero;
        for (std::size_t i = 0; i < M; ++i) {
            T acc = zero;
            for (std::size_t l = 0; l < R; ++l) acc += dm[i](i, l);
            varnow[i] = acc;
            sv += acc;
        }
        T vdev = zero;
        if (sv > zero)
            for (std::size_t i = 0; i < M; ++i) {
                const T d = num_abs(T(varnow[i] - varprev[i])) / sv;
                if (d > vdev) vdev = d;
            }
        varprev = varnow;
        if (dev <= tolm && vdev <= tolv) break;
    }
    return it;
}

}  // namespace detail

/**
 * @param L       (M x R) service demands
 * @param N       (R) population per class, closed only
 * @param Z       (R) think times, empty for none
 * @param tol     stopping tolerance on the mean queue lengths; null for the
 *                reference's own test 1/(4000 + 16 sum(n))
 * @param maxiter maximum CORE iterations, 200 in the reference
 */
template <class T>
SensLinearizerResult<T> pfqn_sens_linearizer(const Matrix<T>& L, const std::vector<int>& N,
                                             const std::vector<T>& Z, const T* tol,
                                             unsigned maxiter) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_sens_linearizer requires transcendental arithmetic: it is a fixed point "
                  "stopped on a tolerance, so its result is a property of the stopping test rather "
                  "than of the model, and exact arithmetic buys nothing");

    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_sens_linearizer: demand matrix and population vector disagree on the class count");
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_sens_linearizer: think-time vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);

    SensLinearizerResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.WN = Matrix<T>(M, R, zero);
    res.m.assign(M, zero);
    res.dm = Matrix<T>(M, M, zero);
    res.d2m.assign(M, zero);
    res.Cov = Matrix<T>(M, M, zero);
    res.Var.assign(M, zero);
    res.M2.assign(M, zero);
    res.M3.assign(M, zero);
    res.Skew.assign(M, 0.0);
    res.CovAsym = zero;
    res.iter = 0;

    bool anyPositive = false;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_sens_linearizer: negative population");
        if (v > 0) anyPositive = true;
    }
    if (!anyPositive || M == 0 || R == 0) return res;

    // pops: index 0 = N, index 1+l = N - e_l
    const std::size_t npops = 1 + R;
    std::vector<std::vector<int>> pv(npops, N);
    for (std::size_t l = 0; l < R; ++l)
        if (N[l] > 0) pv[1 + l][l] = N[l] - 1;

    std::vector<Matrix<T>> mE(npops, Matrix<T>(M, R, zero));
    std::vector<std::vector<Matrix<T>>> dmE(npops, std::vector<Matrix<T>>(M, Matrix<T>(M, R, zero)));
    std::vector<std::vector<Matrix<T>>> d2mE(npops, std::vector<Matrix<T>>(M, Matrix<T>(M, R, zero)));
    for (std::size_t p = 0; p < npops; ++p)
        for (std::size_t l = 0; l < R; ++l)
            for (std::size_t i = 0; i < M; ++i)
                mE[p](i, l) = num_traits<T>::from_int(pv[p][l]) / num_traits<T>::from_int(static_cast<long>(M));

    std::vector<Matrix<T>> delta(M, Matrix<T>(R, R, zero));
    std::vector<std::vector<Matrix<T>>> ddelta(M, std::vector<Matrix<T>>(M, Matrix<T>(R, R, zero)));
    std::vector<std::vector<Matrix<T>>> d2delta(M, std::vector<Matrix<T>>(M, Matrix<T>(R, R, zero)));

    std::vector<T> lam;
    Matrix<T> wmat;
    for (int outer = 0; outer < 3; ++outer) {
        res.iter += detail::sens_lin_core2(L, Z, N, delta, ddelta, d2delta, mE[0], dmE[0], d2mE[0],
                                           lam, wmat, tol, maxiter);
        if (outer == 2) break;

        for (std::size_t l = 0; l < R; ++l) {
            if (N[l] == 0) continue;
            std::vector<T> lam2;
            Matrix<T> w2;
            res.iter += detail::sens_lin_core2(L, Z, pv[1 + l], delta, ddelta, d2delta, mE[1 + l],
                                               dmE[1 + l], d2mE[1 + l], lam2, w2, tol, maxiter);
        }

        // ---- refresh delta from (5.1) and (5.3) ---------------------------
        Matrix<T> vN;
        std::vector<Matrix<T>> dvN, d2vN;
        detail::sens_lin_fractions(mE[0], dmE[0], d2mE[0], N, vN, dvN, d2vN);
        for (std::size_t lp = 0; lp < R; ++lp) {
            if (N[lp] == 0) continue;
            Matrix<T> vL;
            std::vector<Matrix<T>> dvL, d2vL;
            detail::sens_lin_fractions(mE[1 + lp], dmE[1 + lp], d2mE[1 + lp], pv[1 + lp], vL, dvL,
                                       d2vL);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t l = 0; l < R; ++l) {
                    delta[i](lp, l) = vL(i, l) - vN(i, l);
                    for (std::size_t h = 0; h < M; ++h) {
                        ddelta[i][h](lp, l) = dvL[h](i, l) - dvN[h](i, l);
                        d2delta[i][h](lp, l) = d2vL[h](i, l) - d2vN[h](i, l);
                    }
                }
        }
    }

    res.QN = mE[0];
    res.WN = wmat;
    res.XN = lam;
    for (std::size_t i = 0; i < M; ++i) {
        T acc = zero;
        for (std::size_t l = 0; l < R; ++l) acc += res.QN(i, l);
        res.m[i] = acc;
        for (std::size_t h = 0; h < M; ++h) {
            T d = zero;
            for (std::size_t l = 0; l < R; ++l) d += dmE[0][h](i, l);
            res.dm(i, h) = d;
        }
        T d2 = zero;
        for (std::size_t l = 0; l < R; ++l) d2 += d2mE[0][i](i, l);
        res.d2m[i] = d2;
        for (std::size_t r = 0; r < R; ++r) res.UN(i, r) = res.XN[r] * L(i, r);
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            const T d = num_abs(T(res.dm(i, j) - res.dm(j, i)));
            if (d > res.CovAsym) res.CovAsym = d;
        }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j)
            res.Cov(i, j) = (res.dm(i, j) + res.dm(j, i)) / num_traits<T>::from_int(2);

    for (std::size_t i = 0; i < M; ++i) {
        const T d1 = res.dm(i, i);
        const T mi = res.m[i];
        res.Var[i] = d1;
        res.M2[i] = d1 + mi * mi;
        res.M3[i] = res.d2m[i] +
                    (num_traits<T>::from_int(1) + num_traits<T>::from_int(3) * mi) * d1 + mi * mi * mi;
        const T mu3 = res.M3[i] - num_traits<T>::from_int(3) * mi * res.M2[i] +
                      num_traits<T>::from_int(2) * mi * mi * mi;
        const double var = num_traits<T>::to_double(res.Var[i]);
        res.Skew[i] = var > 0.0 ? num_traits<T>::to_double(mu3) / std::pow(var, 1.5) : std::nan("");
    }
    return res;
}

/** pfqn_sens_linearizer with the defaults of the reference. */
template <class T>
SensLinearizerResult<T> pfqn_sens_linearizer(const Matrix<T>& L, const std::vector<int>& N,
                                             const std::vector<T>& Z) {
    return pfqn_sens_linearizer(L, N, Z, static_cast<const T*>(nullptr), 200u);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SENS_LINEARIZER_H
