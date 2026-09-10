/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_TRANSIENT_QBD_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_TRANSIENT_QBD_H

/**
 * Transient analysis of a single-class open queue by the Laplace-domain
 * transient QBD plus numerical inverse Laplace, the port of
 * `matlab/src/solvers/MAM/solver_mam_transient_qbd.m`.
 *
 * Covers MAP/MAP/1 (infinite buffer) and MAP/MAP/1/N (finite buffer): arrival
 * and service are both read as (D0,D1) MAPs, so M/M/1, M/PH/1 and correlated
 * arrival or service all go through the same construction. That is exactly the
 * set the expm / libQBD fast path in `solver_mam_ldqbd_transient` cannot
 * represent, which is why `mam_transient_qbd_applicable` routes here.
 *
 * The level-to-level transform V(s,0,m) comes from `mam_transient2_open`
 * (infinite) or `mam_transient2` (finite) and is inverted per metric with the
 * CME-based numerical inverse Laplace transform.
 *
 * WHY THE INFINITE BRANCH IS EXACT RATHER THAN TRUNCATED. Above the boundary
 * the chain is homogeneous, so V(s,0,m) = V(s,0,1) R^(m-1) for m >= 1 and the
 * level sums close in closed form:
 *
 *     E[N](s) = pi0 V(s,0,1) (I-R)^-2 e
 *     Tput(s) = pi0 V(s,0,1) (I-R)^-1 wDep
 *     P0(s)   = pi0 V(s,0,0) e
 *
 * so no level truncation enters the infinite-buffer answer at all.
 *
 * ARITHMETIC. Double only. The transform is evaluated at complex quadrature
 * nodes and the inversion weights are transcendental, so an exact
 * instantiation refuses by name.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/api/mam/mam_transient2.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/matlab_ilt.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/complex_number.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_ldqbd_transient.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace transient_qbd_detail {

/** Real matrix to complex, entrywise. */
template <class T>
inline CMat to_complex(const Matrix<T>& A) {
    CMat C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            C(i, j) = Complex(num_traits<T>::to_double(A(i, j)), 0.0);
    return C;
}

/** Kronecker product of complex matrices. */
inline CMat ckron(const CMat& A, const CMat& B) {
    CMat C(A.rows() * B.rows(), A.cols() * B.cols(), Complex(0.0, 0.0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            if (A(i, j) == Complex(0.0, 0.0)) continue;
            for (std::size_t k = 0; k < B.rows(); ++k)
                for (std::size_t l = 0; l < B.cols(); ++l)
                    C(i * B.rows() + k, j * B.cols() + l) = A(i, j) * B(k, l);
        }
    return C;
}

/** Kronecker sum, MATLAB's krons. */
inline CMat ckrons(const CMat& A, const CMat& B) {
    const CMat left = ckron(A, eye<Complex>(B.rows()));
    const CMat right = ckron(eye<Complex>(A.rows()), B);
    return transient_detail::madd(left, right);
}

/** pi0 V w, the scalar the inversion actually needs. */
inline Complex quad(const std::vector<Complex>& pi0, const CMat& V, const std::vector<Complex>& w) {
    if (pi0.size() != V.rows() || V.cols() != w.size())
        throw InputError("solver_mam_transient_qbd: quadratic form dimension mismatch");
    Complex acc(0.0, 0.0);
    for (std::size_t i = 0; i < V.rows(); ++i) {
        if (pi0[i] == Complex(0.0, 0.0)) continue;
        Complex row(0.0, 0.0);
        for (std::size_t j = 0; j < V.cols(); ++j) row += V(i, j) * w[j];
        acc += pi0[i] * row;
    }
    return acc;
}

}  // namespace transient_qbd_detail

/**
 * Port of `solver_mam_transient_qbd.m`.
 *
 * @param opt `timespan` bounds the horizon; `iter_max` sets the inversion
 *            budget exactly as the reference does
 * @param L the refreshed struct
 */
template <class T>
TranResult<T> solver_mam_transient_qbd(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_transient_qbd: the level-to-level transform is evaluated at complex "
            "Laplace nodes and inverted with transcendental CME weights; rerun with --arith "
            "double or --arith real");
    } else {
        using namespace transient_qbd_detail;
        using lang::SchedStrategy;
        const std::size_t M = L.nstations, K = L.nclasses;
        if (K != 1)
            throw UnsupportedError(
                "solver_mam_transient_qbd: the transient QBD method requires a single-class model");
        if (!std::isinf(L.classes[0].population))
            throw UnsupportedError(
                "solver_mam_transient_qbd: the transient QBD method requires an open model");

        std::size_t src = 0, q = 0, nsrc = 0, nq = 0;
        for (std::size_t i = 1; i <= M; ++i) {
            if (L.stations[i - 1].sched == SchedStrategy::EXT) { src = i; ++nsrc; }
            else if (L.stations[i - 1].sched == SchedStrategy::FCFS) { q = i; ++nq; }
        }
        if (nsrc != 1 || nq != 1)
            throw UnsupportedError(
                "solver_mam_transient_qbd: the transient QBD method requires exactly one Source "
                "and one FCFS Queue");
        if (L.stations[q - 1].nservers != 1.0)
            throw UnsupportedError(
                "solver_mam_transient_qbd: the Laplace transient QBD supports single-server "
                "queues only");

        // ---- arrival and service MAPs, uniform (D0,D1) ----------------------
        const Map<T> arr = lang::dist_to_map(L.service[src - 1][0]);
        const Map<T> svc = lang::dist_to_map(L.service[q - 1][0]);
        const CMat Da0 = to_complex(arr.D0), Da1 = to_complex(arr.D1);
        const CMat Ds0 = to_complex(svc.D0), Ds1 = to_complex(svc.D1);
        const std::size_t na = Da0.rows(), ns = Ds0.rows();
        const CMat Ina = eye<Complex>(na), Ins = eye<Complex>(ns);

        // Repeating-level blocks (levels >= 1), the qbd_rg convention.
        const CMat Lrep = ckrons(Da0, Ds0);
        const CMat Frep = ckron(Da1, Ins);
        const CMat Brep = ckron(Ina, Ds1);
        // Level 0: the service phase is frozen, only the arrival evolves.
        const CMat Lv0 = ckron(Da0, Ins);
        const CMat F0 = ckron(Da1, Ins);
        const CMat B0 = ckron(Ina, Ds1);

        // Empty system, both phase processes stationary.
        const std::vector<T> piArr = map_prob(arr);
        const std::vector<T> piSvc = map_prob(svc);
        std::vector<Complex> pi0(na * ns);
        for (std::size_t i = 0; i < na; ++i)
            for (std::size_t j = 0; j < ns; ++j)
                pi0[i * ns + j] = Complex(num_traits<T>::to_double(piArr[i]) *
                                              num_traits<T>::to_double(piSvc[j]),
                                          0.0);

        // Service-completion rate over the (arrival, service) phase pairs.
        std::vector<Complex> wDep(na * ns, Complex(0.0, 0.0));
        for (std::size_t i = 0; i < na; ++i)
            for (std::size_t j = 0; j < ns; ++j) {
                Complex s(0.0, 0.0);
                for (std::size_t l = 0; l < ns; ++l) s += Ds1(j, l);
                wDep[i * ns + j] = s;
            }
        const std::vector<Complex> wOne(na * ns, Complex(1.0, 0.0));

        const double bufCap = L.cap[q - 1];
        const bool isFinite = std::isfinite(bufCap);

        // ---- time grid and inversion budget ---------------------------------
        const double T_start = opt.timespan_start;
        const double T_end = opt.timespan_end;
        if (!(T_end > T_start) || !std::isfinite(T_end))
            throw InputError(
                "solver_mam_transient_qbd: the timespan must be a finite interval with a positive "
                "duration");
        const double dur = T_end - T_start;
        const std::size_t nTimePoints = static_cast<std::size_t>(std::min(
            101.0, std::max(11.0, static_cast<double>(std::llround(dur * 10.0)))));
        std::vector<double> times(nTimePoints);
        for (std::size_t i = 0; i < nTimePoints; ++i)
            times[i] = nTimePoints == 1 ? T_start
                                        : T_start + dur * static_cast<double>(i) /
                                                        static_cast<double>(nTimePoints - 1);
        // The inversion is singular at t = 0 (s = beta/t -> Inf). The system
        // starts empty there, so E[N] = U = Tput = 0 exactly; invert only on
        // strictly positive times and write the initial condition directly.
        std::vector<double> tpos;
        std::vector<std::size_t> posIdx;
        for (std::size_t i = 0; i < times.size(); ++i)
            if (times[i] > 0.0) { tpos.push_back(times[i]); posIdx.push_back(i); }

        std::size_t maxFnEvals = 100;
        if (opt.iter_max > 1)
            maxFnEvals = static_cast<std::size_t>(
                std::min(1000L, std::max(11L, static_cast<long>(std::llround(
                                                  static_cast<double>(opt.iter_max))))));

        std::vector<double> EN(times.size(), 0.0), DEP(times.size(), 0.0), Uval(times.size(), 0.0);

        if (isFinite) {
            // Closed piecewise QBD on T = [0, Ncap], one regime.
            const long Ncap = static_cast<long>(std::llround(bufCap));
            if (Ncap < 1)
                throw InputError(
                    "solver_mam_transient_qbd: a finite buffer must hold at least one job");
            // Arrivals are lost at a full buffer, so the top level keeps the
            // arrival phase moving without a level change.
            const CMat LvTop = transient_detail::madd(
                ckron(transient_detail::madd(Da0, Da1), Ins), ckron(Ina, Ds0));
            const TransientQbd qbd = make_transient_qbd({Brep}, {Lrep}, {Frep}, {Lv0, LvTop},
                                                        {0, Ncap});

            const std::vector<double> en = matlab_ilt(
                [&](const Complex& s) {
                    Complex acc(0.0, 0.0);
                    for (long m = 1; m <= Ncap; ++m)
                        acc += static_cast<double>(m) *
                               quad(pi0, mam_transient2(qbd, 0, m, s), wOne);
                    return acc;
                },
                tpos, maxFnEvals);
            const std::vector<double> p0 = matlab_ilt(
                [&](const Complex& s) { return quad(pi0, mam_transient2(qbd, 0, 0, s), wOne); },
                tpos, maxFnEvals);
            const std::vector<double> dep = matlab_ilt(
                [&](const Complex& s) {
                    Complex acc(0.0, 0.0);
                    for (long m = 1; m <= Ncap; ++m)
                        acc += quad(pi0, mam_transient2(qbd, 0, m, s), wDep);
                    return acc;
                },
                tpos, maxFnEvals);
            for (std::size_t i = 0; i < posIdx.size(); ++i) {
                EN[posIdx[i]] = en[i];
                DEP[posIdx[i]] = dep[i];
                Uval[posIdx[i]] = 1.0 - p0[i];
            }
        } else {
            // Open piecewise QBD on T = [0, 1]: regime 1 is the boundary level,
            // regime 2 repeats. L{1} is never referenced from level 0, which is
            // why the reference leaves it empty.
            const TransientQbd qbd = make_transient_qbd({B0, Brep}, {CMat(), Lrep}, {F0, Frep},
                                                        {Lv0, Lrep}, {0, 1});
            const std::size_t nrep = Lrep.rows();

            // (I-R)^-1 and (I-R)^-2 close the geometric level sums exactly.
            const auto geom = [&](const Complex& s, unsigned power,
                                  const std::vector<Complex>& w) {
                const CMat Lk = transient_detail::msub(Lrep, transient_detail::sI(s, nrep));
                const CQbdFundMat gr = qbd_fundmat_laplace(Brep, Lk, Frep);
                const CMat ImR = transient_detail::msub(eye<Complex>(nrep), gr.R);
                CMat rhs(nrep, 1);
                for (std::size_t i = 0; i < nrep; ++i) rhs(i, 0) = w[i];
                for (unsigned p = 0; p < power; ++p) rhs = transient_detail::mldivide(ImR, rhs);
                const CMat V1 = mam_transient2_open(qbd, 0, 1, s);
                std::vector<Complex> col(nrep);
                for (std::size_t i = 0; i < nrep; ++i) col[i] = rhs(i, 0);
                return quad(pi0, V1, col);
            };

            const std::vector<double> en =
                matlab_ilt([&](const Complex& s) { return geom(s, 2, wOne); }, tpos, maxFnEvals);
            const std::vector<double> p0 = matlab_ilt(
                [&](const Complex& s) { return quad(pi0, mam_transient2_open(qbd, 0, 0, s), wOne); },
                tpos, maxFnEvals);
            const std::vector<double> dep =
                matlab_ilt([&](const Complex& s) { return geom(s, 1, wDep); }, tpos, maxFnEvals);
            for (std::size_t i = 0; i < posIdx.size(); ++i) {
                EN[posIdx[i]] = en[i];
                DEP[posIdx[i]] = dep[i];
                Uval[posIdx[i]] = 1.0 - p0[i];
            }
        }

        TranResult<T> out;
        out.Qt.assign(M, std::vector<TranCurve<T>>(K));
        out.Ut.assign(M, std::vector<TranCurve<T>>(K));
        out.Tt.assign(M, std::vector<TranCurve<T>>(K));
        std::vector<T> qv(times.size()), uv(times.size()), tv(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            qv[i] = num_traits<T>::from_double(EN[i]);
            uv[i] = num_traits<T>::from_double(Uval[i]);
            tv[i] = num_traits<T>::from_double(DEP[i]);
        }
        out.Qt[q - 1][0].values = qv;
        out.Qt[q - 1][0].times = times;
        out.Ut[q - 1][0].values = uv;
        out.Ut[q - 1][0].times = times;
        out.Tt[q - 1][0].values = tv;
        out.Tt[q - 1][0].times = times;
        return out;
    }
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_TRANSIENT_QBD_H
