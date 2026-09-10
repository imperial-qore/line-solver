/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_RQNA_H
#define LINE_SOLVERS_MVA_SOLVER_RQNA_H

/**
 * Robust Queueing Network Analyzer (RQNA), a port of
 * matlab/src/solvers/MVA/solver_rqna.m.
 *
 * Single-class open network of single-server FCFS queues with Markovian routing
 * and general (non-renewal) external arrival and (non-exponential) service.
 * Implements Whitt & You (2018) Algorithm 1: traffic-rate + limiting-variability
 * equations, time-dependent IDC equations (via npfqn_traffic_idc /
 * npfqn_traffic_idc_at), and the robust-queueing workload approximation
 * (qsys_gig1_rq), plus the near-immediate feedback elimination of Algorithm 2 /
 * Section 4.1-4.2, applied by default.
 *
 * REFERENCE INDEXING, REPRODUCED. Like solver_qna, the reference indexes the
 * stateful-indexed sn.rt with STATION indices, which is exact only when every
 * stateful node is a station; that precondition is asserted. Single class, so
 * sn.rt(i,j) is the station-to-station routing directly.
 *
 * ARITHMETIC: transcendental (counting-process IDC needs expm, the RQ workload
 * a square root), so under Rational the whole body is discarded and RQNA is
 * refused by name, matching solver_qna.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/mam/map_count_idc.h"
#include "line/api/mam/map_moment.h"
#include "line/api/npfqn/npfqn_feedback_elim.h"
#include "line/api/npfqn/npfqn_traffic_idc.h"
#include "line/api/qsys/qsys_gig1_rq.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mva {

namespace rqna_detail {

/** Submatrix A(rows, cols) by index lists (0-based). */
template <class T>
Matrix<T> submat(const Matrix<T>& A, const std::vector<std::size_t>& r,
                 const std::vector<std::size_t>& c) {
    Matrix<T> B(r.size(), c.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < r.size(); ++i)
        for (std::size_t j = 0; j < c.size(); ++j) B(i, j) = A(r[i], c[j]);
    return B;
}

/**
 * Near-immediate feedback probability at station a (Whitt-You eq. 3.8/3.9,
 * H={a}): the probability of returning to a before visiting any strictly
 * higher-rho station.
 *
 * Delegated to `npfqn_feedback_elim` so that the solver and the API function
 * cannot drift apart: they answer the same question, and a private copy of the
 * rule here is how the two came to differ on ties in the first place.
 */
template <class T>
T rqna_phat(const Matrix<T>& P, const std::vector<T>& rho, std::size_t a) {
    return npfqn::npfqn_feedback_elim(P, rho).feedbackProb[a];
}

/**
 * Geometric random sum of i.i.d. PH service with success prob (1-p):
 * PH(alpha, T) -> D0 + p t0 alpha, D1 = (1-p) t0 alpha, t0 = -D0 e.
 */
template <class T>
mam::Map<T> rqna_geom_map(const mam::Map<T>& map, const T& p) {
    const std::size_t s = map.D0.rows();
    std::vector<T> e(s, num_traits<T>::from_int(1));
    std::vector<T> t0(s, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < s; ++i) {
        T acc = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < s; ++j) acc = T(acc - map.D0(i, j) * e[j]);
        t0[i] = acc;
    }
    const std::vector<T> al = mam::map_pie(map);
    mam::Map<T> out;
    out.D0 = Matrix<T>(s, s, num_traits<T>::from_int(0));
    out.D1 = Matrix<T>(s, s, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < s; ++i)
        for (std::size_t j = 0; j < s; ++j) {
            out.D0(i, j) = T(map.D0(i, j) + p * (t0[i] * al[j]));
            out.D1(i, j) = T((num_traits<T>::from_int(1) - p) * (t0[i] * al[j]));
        }
    return out;
}

}  // namespace rqna_detail

/**
 * @param L   the single-class open network struct
 * @param opt options (opt.tol is the feedback-elimination threshold). The
 *            reference's config knobs rqna_feedback_elim / rqna_alpha /
 *            rqna_beta are not exposed by MvaOptions, so this port runs their
 *            MATLAB defaults: feedback elimination ON, no IDC corrections.
 * @return the AvgResult-style MvaSolution
 */
template <class T>
MvaSolution<T> solver_rqna(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_rqna: the robust queueing-network analyzer needs transcendental arithmetic "
            "(counting-process indices of dispersion); rerun with --arith double or --arith real");
    } else {
    using rqna_detail::rqna_phat;
    using rqna_detail::rqna_geom_map;
    using rqna_detail::submat;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // One predicate for the gate and the run: list_valid_methods asks the same
    // question before it advertises 'rqna', so the sentence a caller reads here is
    // the sentence that kept the row off the report.
    {
        const std::string rqna_reason = mva_single_class_open_reason(L, "rqna");
        if (!rqna_reason.empty()) throw UnsupportedError(rqna_reason);
    }
    const std::size_t M = L.nstations;
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (std::isfinite(L.classes[r].population))
            throw UnsupportedError("solver_rqna: RQNA supports open networks only (no closed classes)");
    if (L.nof_stateful() != M)
        throw UnsupportedError(
            "solver_rqna: the reference indexes the stateful-indexed sn.rt with station indices, "
            "which is only correct when every stateful node is a station");

    MvaSolution<T> s;
    s.Q = Matrix<T>(M, 1, zero);
    s.U = Matrix<T>(M, 1, zero);
    s.R = Matrix<T>(M, 1, zero);
    s.Tp = Matrix<T>(M, 1, zero);
    s.X.assign(1, zero);
    s.C.assign(1, zero);
    s.method = opt.method;
    s.iter = 1;

    // Source and queueing (incl. delay) stations.
    std::size_t src = M;  // sentinel
    std::vector<bool> isSource(M, false), schedInf(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        isSource[i] = (L.stations[i].nodetype == qn::NodeType::Source);
        schedInf[i] = (L.stations[i].sched == qn::SchedStrategy::INF);
        if (isSource[i] && src == M) src = i;
    }
    if (src == M)
        throw UnsupportedError("solver_rqna: RQNA requires an open network with a Source station");
    std::vector<std::size_t> qstat;
    for (std::size_t i = 0; i < M; ++i)
        if (!isSource[i]) qstat.push_back(i);
    const std::size_t nq = qstat.size();

    // Station-to-station single-class routing sn.rt (K == 1).
    auto rtS = [&](std::size_t i, std::size_t j) -> T { return L.rt(i, j); };

    // External arrival process from the Source.
    const mam::Map<T> arvMAP = lang::dist_to_map(L.service[src][0]);
    const T lambda_src = mam::map_lambda(arvMAP);
    const T c2_src = mam::map_idc(arvMAP);

    // Per-queue data.
    std::vector<T> mu(nq, zero), cs2(nq, zero), lambda0(nq, zero), qsplit(nq, zero);
    std::vector<mam::Map<T> > svcMAP;
    svcMAP.reserve(nq);
    Matrix<T> P(nq, nq, zero);
    for (std::size_t a = 0; a < nq; ++a) {
        const std::size_t ia = qstat[a];
        mu[a] = L.rates(ia, 0);
        cs2[a] = L.scv(ia, 0);
        svcMAP.push_back(lang::dist_to_map(L.service[ia][0]));
        qsplit[a] = rtS(src, ia);
        lambda0[a] = T(lambda_src * qsplit[a]);
        for (std::size_t b = 0; b < nq; ++b) P(a, b) = rtS(ia, qstat[b]);
    }

    // External-arrival IDC seen by each queue (eq. 31): split of the source.
    std::vector<T> c2a0(nq, zero);
    for (std::size_t a = 0; a < nq; ++a) c2a0[a] = T(qsplit[a] * c2_src + (one - qsplit[a]));

    // IDC handles for the time-dependent solve.
    std::function<std::vector<T>(const T&)> a0IdcFun = [&](const T& t) {
        std::vector<T> out(nq, zero);
        const T Iarv = mam::map_count_idc(arvMAP, t);
        for (std::size_t a = 0; a < nq; ++a) out[a] = T(qsplit[a] * Iarv + (one - qsplit[a]));
        return out;
    };
    std::function<std::vector<T>(const std::vector<T>&)> sIdcFun =
        [&](const std::vector<T>& tt) {
            std::vector<T> out(nq, zero);
            for (std::size_t a = 0; a < nq; ++a) out[a] = mam::map_count_idc(svcMAP[a], tt[a]);
            return out;
        };

    // Corrections (eqs. 34, 38-39), configurable.
    npfqn::TrafficIdcCorrections corr;
    // Defaults follow npfqn_traffic_idc's own defaults; config overrides not
    // wired through MvaOptions here (parity with the analytic default path).

    const npfqn::TrafficIdcContext<T> ctx =
        npfqn::npfqn_traffic_idc(lambda0, P, c2a0, mu, cs2, corr);
    const std::vector<T>& lambda = ctx.lambda;
    const std::vector<T>& rho = ctx.rho;

    // component IDC I_a(x) reached through ctx (npfqn_traffic_idc_at).
    auto component = [&](const T& x, std::size_t a) -> T {
        return npfqn::npfqn_traffic_idc_at(ctx, x, a0IdcFun, sIdcFun)[a];
    };

    for (std::size_t a = 0; a < nq; ++a) {
        const std::size_t ia = qstat[a];
        s.Tp(ia, 0) = lambda[a];
        if (!(lambda[a] > zero)) continue;
        if (schedInf[ia]) {
            s.U(ia, 0) = T(lambda[a] / mu[a]);
            s.Q(ia, 0) = T(lambda[a] / mu[a]);
            s.R(ia, 0) = T(one / mu[a]);
            continue;
        }
        T phat = zero;
        phat = rqna_phat(P, rho, a);
        if (phat > num_traits<T>::from_double(opt.tol)) {
            // Near-immediate feedback elimination at station a.
            const T tol = num_traits<T>::from_double(1e-9);
            std::vector<std::size_t> Hc, Hi;
            for (std::size_t i = 0; i < nq; ++i) {
                if (i == a) continue;
                if (rho[i] <= T(rho[a] + tol)) Hc.push_back(i);
                else Hi.push_back(i);
            }
            std::vector<std::size_t> R;  // retained, a first
            R.push_back(a);
            for (std::size_t i : Hi) R.push_back(i);
            const std::size_t m = R.size();

            Matrix<T> Fhc, G, Pred;
            if (Hc.empty()) {
                G = Matrix<T>(0, m, zero);
                Pred = submat(P, R, R);
            } else {
                const std::size_t nc = Hc.size();
                Matrix<T> ImP(nc, nc, zero);
                for (std::size_t i = 0; i < nc; ++i) {
                    ImP(i, i) = one;
                    for (std::size_t j = 0; j < nc; ++j) ImP(i, j) = T(ImP(i, j) - P(Hc[i], Hc[j]));
                }
                Fhc = inverse(ImP);
                G = matmul(Fhc, submat(P, Hc, R));                  // nc x m
                Pred = submat(P, R, R);
                const Matrix<T> corr2 = matmul(matmul(submat(P, R, Hc), Fhc), submat(P, Hc, R));
                for (std::size_t i = 0; i < m; ++i)
                    for (std::size_t j = 0; j < m; ++j) Pred(i, j) = T(Pred(i, j) + corr2(i, j));
            }
            T ph = Pred(0, 0);
            if (ph < zero) ph = zero;
            if (ph > num_traits<T>::from_double(1.0 - 1e-9)) ph = num_traits<T>::from_double(1.0 - 1e-9);
            if (ph > zero)
                for (std::size_t j = 0; j < m; ++j) Pred(0, j) = T(Pred(0, j) / (one - ph));
            Pred(0, 0) = zero;

            // First-passage external rate and IDC into each retained station.
            const T idc_arv_inf = mam::map_idc(arvMAP);
            std::vector<T> lam0R(m, zero);
            for (std::size_t rr = 0; rr < m; ++rr) {
                lam0R[rr] = lambda0[R[rr]];
                for (std::size_t ii = 0; ii < Hc.size(); ++ii)
                    lam0R[rr] = T(lam0R[rr] + lambda0[Hc[ii]] * G(ii, rr));
            }
            std::vector<T> c2a0R(m, zero);
            for (std::size_t rr = 0; rr < m; ++rr) {
                if (!(lam0R[rr] > zero)) continue;
                T acc = T(lambda0[R[rr]] * (qsplit[R[rr]] * idc_arv_inf + (one - qsplit[R[rr]])));
                for (std::size_t ii = 0; ii < Hc.size(); ++ii) {
                    const T g = G(ii, rr);
                    const T ci = T(qsplit[Hc[ii]] * idc_arv_inf + (one - qsplit[Hc[ii]]));
                    acc = T(acc + lambda0[Hc[ii]] * g * (g * ci + (one - g)));
                }
                c2a0R[rr] = T(acc / lam0R[rr]);
            }
            std::function<std::vector<T>(const T&)> a0IdcR = [&, R, Hc, G, lam0R, m](const T& t) {
                std::vector<T> Ir(m, one);
                const T Iarv = mam::map_count_idc(arvMAP, t);
                auto c2fun_i = [&](std::size_t i) -> T {
                    return T(qsplit[i] * Iarv + (one - qsplit[i]));
                };
                for (std::size_t rr = 0; rr < m; ++rr) {
                    if (!(lam0R[rr] > zero)) continue;
                    T num = T(lambda0[R[rr]] * c2fun_i(R[rr]));
                    for (std::size_t ii = 0; ii < Hc.size(); ++ii) {
                        const T g = G(ii, rr);
                        num = T(num + lambda0[Hc[ii]] * g * (g * c2fun_i(Hc[ii]) + (one - g)));
                    }
                    Ir[rr] = T(num / lam0R[rr]);
                }
                return Ir;
            };

            // Service data on R; a gets the geometric-sum (folded) service.
            std::vector<T> muR(m, zero), cs2R(m, zero);
            std::vector<mam::Map<T> > svcR;
            svcR.reserve(m);
            for (std::size_t rr = 0; rr < m; ++rr) {
                muR[rr] = mu[R[rr]];
                cs2R[rr] = cs2[R[rr]];
                svcR.push_back(svcMAP[R[rr]]);
            }
            svcR[0] = rqna_geom_map(svcMAP[a], ph);
            muR[0] = T((one - ph) * mu[a]);
            cs2R[0] = T(ph + (one - ph) * cs2[a]);
            std::function<std::vector<T>(const std::vector<T>&)> sIdcR =
                [&, m](const std::vector<T>& tt) {
                    std::vector<T> out(m, zero);
                    for (std::size_t rr = 0; rr < m; ++rr)
                        out[rr] = mam::map_count_idc(svcR[rr], tt[rr]);
                    return out;
                };

            const npfqn::TrafficIdcContext<T> ctxR =
                npfqn::npfqn_traffic_idc(lam0R, Pred, c2a0R, muR, cs2R, corr);
            auto IaFunA = [&](const T& x) -> T {
                return npfqn::npfqn_traffic_idc_at(ctxR, x, a0IdcR, sIdcR)[0];
            };
            const T Wt = qsys::qsys_gig1_rq(rho[a], muR[0], cs2R[0], IaFunA).W;
            s.R(ia, 0) = T((one - ph) * Wt + one / mu[a]);
        } else {
            auto IaFun_a = [&](const T& x) -> T { return component(x, a); };
            const T Wa = qsys::qsys_gig1_rq(rho[a], mu[a], cs2[a], IaFun_a).W;
            s.R(ia, 0) = T(Wa + one / mu[a]);
        }
        s.U(ia, 0) = rho[a];
        s.Q(ia, 0) = T(lambda[a] * s.R(ia, 0));
    }

    s.Tp(src, 0) = lambda_src;
    T Csum = zero;
    for (std::size_t i = 0; i < M; ++i) Csum = T(Csum + s.R(i, 0));
    s.C[0] = Csum;
    s.X[0] = lambda_src;
    return s;
    }  // if constexpr has_transcendental
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_RQNA_H
