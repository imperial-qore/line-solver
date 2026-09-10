/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_H
#define LINE_API_INFER_INFER_LQN_H

/**
 * Identify hidden LQN parameters from measured performance data.
 *
 * Templated port of matlab/src/api/infer/infer_lqn.m. No JAR counterpart.
 * Estimates the LQN parameters named in `spec` (activity host demands and
 * task think times) from the sequence of measurements Z, using the Extended
 * Kalman Filter of infer_lqn_ekf.h over the observation model defined by
 * `obs`. It implements Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking
 * Time-Varying Parameters in Software Systems with Extended Kalman Filters",
 * CASCON 2005. A single measurement column with QFac = 0 reduces to one-shot
 * least-squares calibration.
 *
 * The DEFAULT COVARIANCES are the paper's equations (9a) and (9b):
 *   Q_ii = (QFac a0_i cvA)^2, taking mean(a_i) ~ a0_i
 *   R_ii = ((RFac zbar_i)/1.96)^2 / gammaT, with zbar_i the mean of the i-th
 *          measured row over the steps
 * and P0 = diag((0.5 a0)^2). Each is overridable, and each is floored at the
 * machine epsilon exactly as MATLAB floors it: a parameter whose initial
 * estimate is zero would otherwise give a singular Q and stall the filter on
 * that coordinate forever.
 *
 * gammaT is T/Tstar when both are supplied and 1 otherwise, which is what
 * MATLAB's cascade of infer_lqn_optget calls resolves to.
 *
 * THE OBSERVATION MODEL is a required argument here rather than the defaulted
 * @SolverLN of MATLAB. It is MATLAB's options.solver, and in C++ the caller
 * passes the solve step directly: injecting a default would make api/infer
 * depend on the whole layered solver, and the algorithm is identical either
 * way. The callable receives the struct with the candidate parameters already
 * injected and returns the per-element metric vectors that infer_lqn_getobs
 * selects from.
 *
 * The final estimate is applied to the struct before returning, so the caller
 * gets back a model carrying the identified parameters, as in MATLAB.
 *
 * ARITHMETIC: the filter needs transcendental arithmetic (see
 * infer_lqn_ekf.h), so the exact instantiation is refused.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/infer/infer_lqn_ekf.h"
#include "line/api/infer/infer_lqn_getobs.h"
#include "line/api/infer/infer_lqn_setparams.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * MATLAB's OPTIONS struct for infer_lqn, with the same defaults the
 * infer_lqn_optget calls supply. An empty override matrix means "use the
 * default construction"; a non-empty one is used verbatim.
 */
template <class T>
struct InferLqnOptions {
    double QFac = 0.1;          ///< drift-noise factor
    double RFac = 0.2;          ///< measurement-noise factor
    double cvA = 1.0;           ///< parameter drift coefficient of variation
    bool has_gammaT = false;    ///< true when gammaT is given outright
    double gammaT = 1.0;        ///< T/Tstar ratio used in R
    bool has_T = false;         ///< true when the measurement interval is given
    double T_interval = 0.0;    ///< measurement interval
    bool has_Tstar = false;     ///< true when the system constant is given
    double Tstar = 0.0;         ///< system constant
    Matrix<T> Q;                ///< explicit drift covariance, empty for the default
    Matrix<T> R;                ///< explicit measurement covariance, empty for the default
    Matrix<T> P0;               ///< explicit initial covariance, empty for the default
    std::vector<T> a0;          ///< initial estimate, empty to read it off the model
    std::vector<T> a_true;      ///< ground truth, enables the Ea metric
    double fdStep = 1e-3;       ///< relative finite-difference step
    double fdFloor = 1e-6;      ///< minimum absolute perturbation scale
    bool clampPositive = true;  ///< clamp each estimate to at least fdFloor
};

/** MATLAB's INFO struct: the EKF result plus what the driver constructed. */
template <class T>
struct InferLqnResult {
    EkfResult<T> ekf;   ///< the filter output, ahat included
    std::vector<T> a0;  ///< initial estimate actually used
    Matrix<T> Q;        ///< drift covariance actually used
    Matrix<T> R;        ///< measurement covariance actually used
    Matrix<T> P0;       ///< initial covariance actually used
};

/**
 * @param lsn   layered struct, mutated in place to carry the final estimate
 * @param spec  parameters to identify
 * @param obs   observations to read at each step, numel(obs) == Z.rows()
 * @param Z     (no x nsteps) measurements, one column per step
 * @param solve observation model: given the struct with candidate parameters
 *              injected, return the per-element metric vectors
 * @param opt   filter and covariance options
 */
template <class T>
InferLqnResult<T> infer_lqn(
    lqn::LqnStruct<T>& lsn, const std::vector<LqnParamSpec>& spec,
    const std::vector<LqnObsSpec>& obs, const Matrix<T>& Z,
    const std::function<LqnMetrics<T>(const lqn::LqnStruct<T>&)>& solve,
    const InferLqnOptions<T>& opt = InferLqnOptions<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "infer_lqn requires transcendental arithmetic: its filter corrects on a "
                  "finite-difference sensitivity matrix and reports root mean square errors");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t no = obs.size();
    if (Z.rows() != no) throw InputError("infer_lqn: row count of Z must equal numel(obsSpec)");
    if (spec.empty()) throw InputError("infer_lqn: no parameters to identify");
    if (!solve) throw InputError("infer_lqn: an observation model is required");

    // initial estimate: the option, or the values currently on the model
    std::vector<T> a0 = opt.a0;
    if (a0.empty()) a0 = infer_lqn_getparams(lsn, spec);
    if (a0.size() != spec.size())
        throw InputError("infer_lqn: a0 has the wrong length");
    const std::size_t np = a0.size();

    // MATLAB's eps floor, applied to every constructed diagonal
    const T eps = num_traits<T>::from_double(2.220446049250313e-16);

    double gammaT = 1.0;
    if (opt.has_gammaT)
        gammaT = opt.gammaT;
    else if (opt.has_T && opt.has_Tstar)
        gammaT = opt.T_interval / opt.Tstar;

    InferLqnResult<T> out;
    out.a0 = a0;

    if (opt.Q.rows() != 0) {
        out.Q = opt.Q;
    } else {
        out.Q = Matrix<T>(np, np, zero);
        const T qf = num_traits<T>::from_double(opt.QFac * opt.cvA);
        for (std::size_t i = 0; i < np; ++i) {
            const T qi = qf * num_abs(a0[i]);
            const T v = qi * qi;
            out.Q(i, i) = v > eps ? v : eps;
        }
    }

    if (opt.R.rows() != 0) {
        out.R = opt.R;
    } else {
        out.R = Matrix<T>(no, no, zero);
        const T nsteps = num_traits<T>::from_int(static_cast<long>(Z.cols()));
        const T rf = num_traits<T>::from_double(opt.RFac / 1.96);
        const T gt = num_traits<T>::from_double(gammaT);
        for (std::size_t i = 0; i < no; ++i) {
            T zbar = zero;
            for (std::size_t k = 0; k < Z.cols(); ++k) zbar += Z(i, k);
            zbar = zbar / nsteps;
            const T ri = rf * num_abs(zbar);
            const T v = ri * ri / gt;
            out.R(i, i) = v > eps ? v : eps;
        }
    }

    if (opt.P0.rows() != 0) {
        out.P0 = opt.P0;
    } else {
        out.P0 = Matrix<T>(np, np, zero);
        const T half = num_traits<T>::from_double(0.5);
        for (std::size_t i = 0; i < np; ++i) {
            const T pi_ = half * num_abs(a0[i]);
            const T v = pi_ * pi_;
            out.P0(i, i) = v > eps ? v : eps;
        }
    }

    // observation model h(a): inject the parameters, solve, read the metrics
    lqn::LqnStruct<T>* model = &lsn;
    const std::vector<LqnParamSpec>* sp = &spec;
    const std::vector<LqnObsSpec>* ob = &obs;
    std::function<std::vector<T>(const std::vector<T>&)> hfun =
        [model, sp, ob, &solve](const std::vector<T>& a) {
            infer_lqn_setparams(*model, *sp, a);
            const LqnMetrics<T> m = solve(*model);
            return infer_lqn_getobs(model->names, m, *ob);
        };

    EkfOptions<T> ek;
    ek.fd_step = opt.fdStep;
    ek.fd_floor = opt.fdFloor;
    ek.clamp_positive = opt.clampPositive;
    ek.a_true = opt.a_true;
    out.ekf = infer_lqn_ekf<T>(hfun, a0, out.P0, Z, out.Q, out.R, ek);

    // apply the final estimate to the returned model
    std::vector<T> afinal(np, zero);
    const std::size_t last = out.ekf.ahat.cols() - 1;
    for (std::size_t i = 0; i < np; ++i) afinal[i] = out.ekf.ahat(i, last);
    infer_lqn_setparams(lsn, spec, afinal);
    return out;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_H
