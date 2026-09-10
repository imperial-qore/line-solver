/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_EKF_H
#define LINE_API_INFER_INFER_LQN_EKF_H

/**
 * Extended Kalman Filter for LQN parameter identification.
 *
 * Templated port of matlab/src/api/infer/infer_lqn_ekf.m. The JAR carries the
 * same recursion inside jline/api/infer/InferLqn.java.
 *
 * Tracks a hidden parameter vector across a measurement sequence, following
 * Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying Parameters in
 * Software Systems with Extended Kalman Filters", CASCON 2005, equations 1-9.
 * The parameter follows a zero-mean random walk a_k = a_{k-1} + w and the
 * measurement is z_k = h(a_k) + v, with h the nonlinear performance model.
 * Per step:
 *
 *   a_pred = a,  P_pred = P + Q                       (prediction, eq 5)
 *   H, zpred = jacobian of h at a_pred                (eq 4)
 *   e = z_k - zpred                                   (innovation)
 *   S = H P_pred H' + R,  K = P_pred H' S^-1          (gain, eq 6)
 *   a = a_pred + K e,  P = (I - K H) P_pred           (update, eq 7)
 *
 * h enters as a std::function, exactly as in infer_lqn_jacobian, so the model
 * layer that evaluates an LQN stays outside this header: any caller that can
 * produce z = h(a) can run the filter. That is the whole reason the MATLAB
 * splits infer_lqn into infer_lqn_ekf plus a model-evaluating closure.
 *
 * MATLAB forms the gain as (Ppred*H')/S, a right division, i.e. it SOLVES
 * K S = Ppred H' rather than inverting S. The port does the same, row by row
 * through one LU factorization of S', which matters: S is the innovation
 * covariance and is ill-conditioned exactly when a parameter is unobservable,
 * which is the regime the filter is used in. The covariance is symmetrized
 * after the update, as in MATLAB, because (I - K H) P_pred is symmetric only
 * up to roundoff and an asymmetric P feeds back into every later gain.
 *
 * MATLAB's OPTIONS struct is read through infer_lqn_optget, a
 * field-present-and-non-empty default helper. That helper has no counterpart
 * here and needs none: EkfOptions below carries the same five defaults as
 * member initializers, which is what a struct with defaulted fields IS in
 * C++. The one option not carried is 'verbose', a console trace: the port
 * writes nothing to any stream (an API-layer invariant), and the caller has
 * the same information in the returned per-step innovations.
 *
 * ARITHMETIC: the recursion itself is rational, but Er and Ea are root mean
 * squares, so the routine is gated on transcendental arithmetic and registered
 * for Double and Real only. That gate is not merely about the two summary
 * numbers: the filter is a fixed number of steps of a linear-algebraic
 * recursion, so in exact arithmetic the numerators of P would grow without
 * bound while the estimate it produces stays an approximation of a nonlinear
 * problem, and the finite differences inside H carry a truncation error that
 * exact arithmetic cannot remove.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/infer/infer_lqn_jacobian.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/** MATLAB's OPTIONS struct, with the same defaults infer_lqn_optget supplies. */
template <class T>
struct EkfOptions {
    double fd_step = 1e-3;         ///< relative finite-difference step
    double fd_floor = 1e-6;        ///< minimum absolute perturbation scale
    bool clamp_positive = true;    ///< clamp each estimate to at least fd_floor
    std::vector<T> a_true;         ///< ground truth for Ea, empty for none
};

/** MATLAB's [ahat, info] return list. */
template <class T>
struct EkfResult {
    Matrix<T> ahat;                 ///< (np x nsteps) estimate trajectory
    Matrix<T> e;                    ///< (no x nsteps) prediction errors
    Matrix<T> zpred;                ///< (no x nsteps) predicted measurements
    std::vector<Matrix<T>> Phist;   ///< per-step posterior covariances
    Matrix<T> P;                    ///< final covariance
    T Er;                           ///< RMS prediction error
    T Ea;                           ///< RMS parameter tracking error, if a_true given
    bool has_Ea = false;            ///< false mirrors MATLAB's info.Ea = []
};

/**
 * @param hfun observation map, a parameter vector to an observation vector
 * @param a0   (np) initial parameter estimate
 * @param P0   (np x np) initial estimation-error covariance
 * @param Z    (no x nsteps) measurements, one column per step
 * @param Q    (np x np) parameter-drift covariance
 * @param R    (no x no) measurement-error covariance
 * @param opts filter options
 */
template <class T>
EkfResult<T> infer_lqn_ekf(const std::function<std::vector<T>(const std::vector<T>&)>& hfun,
                           const std::vector<T>& a0, const Matrix<T>& P0, const Matrix<T>& Z,
                           const Matrix<T>& Q, const Matrix<T>& R,
                           const EkfOptions<T>& opts = EkfOptions<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "infer_lqn_ekf requires transcendental arithmetic: its Er and Ea summaries "
                  "are root mean squares, and the sensitivity matrix it corrects on is a "
                  "finite difference whose truncation error no arithmetic can remove");

    const std::size_t np = a0.size();
    const std::size_t no = Z.rows();
    const std::size_t nsteps = Z.cols();
    if (np == 0) throw InputError("infer_lqn_ekf: empty parameter vector");
    if (nsteps == 0) throw InputError("infer_lqn_ekf: no measurements");
    if (P0.rows() != np || P0.cols() != np) throw InputError("infer_lqn_ekf: P0 is not np x np");
    if (Q.rows() != np || Q.cols() != np) throw InputError("infer_lqn_ekf: Q is not np x np");
    if (R.rows() != no || R.cols() != no) throw InputError("infer_lqn_ekf: R is not no x no");
    if (!opts.a_true.empty() && opts.a_true.size() != np)
        throw InputError("infer_lqn_ekf: aTrue has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T half = num_traits<T>::from_double(0.5);
    const T floor = num_traits<T>::from_double(opts.fd_floor);

    EkfResult<T> out;
    out.ahat = Matrix<T>(np, nsteps, zero);
    out.e = Matrix<T>(no, nsteps, zero);
    out.zpred = Matrix<T>(no, nsteps, zero);
    out.Phist.reserve(nsteps);

    std::vector<T> a = a0;
    Matrix<T> P = P0;

    for (std::size_t k = 0; k < nsteps; ++k) {
        // (1) predict: the drift is zero mean, so only the covariance moves
        const std::vector<T> aPred = a;
        Matrix<T> Ppred(np, np, zero);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) Ppred(i, j) = P(i, j) + Q(i, j);

        // (2,4) predicted measurement and sensitivity matrix at a_pred
        const JacobianResult<T> jac =
            infer_lqn_jacobian(hfun, aPred, opts.fd_step, opts.fd_floor);
        if (jac.h0.size() != no)
            throw InputError("infer_lqn_ekf: the observation map and Z disagree on the "
                             "observation count");
        const Matrix<T>& H = jac.H;

        // (3) prediction error
        std::vector<T> e(no, zero);
        for (std::size_t i = 0; i < no; ++i) e[i] = Z(i, k) - jac.h0[i];

        // (6) gain: B = Ppred H' (np x no), S = H B + R (no x no), K S = B
        Matrix<T> B(np, no, zero);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < no; ++j) {
                T s = zero;
                for (std::size_t l = 0; l < np; ++l) s += Ppred(i, l) * H(j, l);
                B(i, j) = s;
            }
        Matrix<T> St(no, no, zero);  // S transposed, the matrix of the row solves
        for (std::size_t i = 0; i < no; ++i)
            for (std::size_t j = 0; j < no; ++j) {
                T s = R(i, j);
                for (std::size_t l = 0; l < np; ++l) s += H(i, l) * B(l, j);
                St(j, i) = s;
            }
        Matrix<T> LU = St;
        const std::vector<std::size_t> piv = lu_factor(LU);
        Matrix<T> K(np, no, zero);
        for (std::size_t i = 0; i < np; ++i) {
            std::vector<T> rhs(no, zero);
            for (std::size_t j = 0; j < no; ++j) rhs[j] = B(i, j);
            lu_solve(LU, piv, rhs);
            for (std::size_t j = 0; j < no; ++j) K(i, j) = rhs[j];
        }

        // (4-update) improved parameter estimate
        for (std::size_t i = 0; i < np; ++i) {
            T s = aPred[i];
            for (std::size_t j = 0; j < no; ++j) s += K(i, j) * e[j];
            a[i] = s;
        }
        if (opts.clamp_positive)
            for (std::size_t i = 0; i < np; ++i)
                if (a[i] < floor) a[i] = floor;

        // (7) covariance update, symmetrized as in MATLAB
        Matrix<T> KH(np, np, zero);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) {
                T s = zero;
                for (std::size_t l = 0; l < no; ++l) s += K(i, l) * H(l, j);
                KH(i, j) = s;
            }
        Matrix<T> Pnew(np, np, zero);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j) {
                T s = Ppred(i, j);
                for (std::size_t l = 0; l < np; ++l) s -= KH(i, l) * Ppred(l, j);
                Pnew(i, j) = s;
            }
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = i; j < np; ++j) {
                const T s = half * (Pnew(i, j) + Pnew(j, i));
                Pnew(i, j) = s;
                Pnew(j, i) = s;
            }
        P = Pnew;

        for (std::size_t i = 0; i < np; ++i) out.ahat(i, k) = a[i];
        for (std::size_t i = 0; i < no; ++i) {
            out.e(i, k) = e[i];
            out.zpred(i, k) = jac.h0[i];
        }
        out.Phist.push_back(P);
    }
    out.P = P;

    using std::sqrt;
    T se = zero;
    for (std::size_t i = 0; i < no; ++i)
        for (std::size_t k = 0; k < nsteps; ++k) se += out.e(i, k) * out.e(i, k);
    const std::size_t ne = no * nsteps;
    out.Er = ne == 0 ? zero
                     : T(sqrt(T(se / num_traits<T>::from_int(static_cast<long>(ne)))));
    out.Ea = zero;
    if (!opts.a_true.empty()) {
        T sa = zero;
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < nsteps; ++k) {
                const T d = out.ahat(i, k) - opts.a_true[i];
                sa += d * d;
            }
        const std::size_t na = np * nsteps;
        out.Ea = T(sqrt(T(sa / num_traits<T>::from_int(static_cast<long>(na)))));
        out.has_Ea = true;
    }
    return out;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_EKF_H
