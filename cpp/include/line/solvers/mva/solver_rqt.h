/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_RQT_H
#define LINE_SOLVERS_MVA_SOLVER_RQT_H

/**
 * Robust Queueing Network Analyzer (RQNA) of Robust Queueing Theory, a port of
 * matlab/src/solvers/MVA/solver_rqt.m, cross-checked against
 * jar/src/main/java/jline/solvers/mva/handlers/Solver_rqt.java.
 *
 * Single-class open network of FCFS queues with Markovian routing. The
 * stochastic primitives are replaced by polyhedral uncertainty sets and each
 * node is analysed in isolation under a worst case. The algorithm is Section 7.2
 * of the reference: the external streams get Gamma_a = sigma_a, the effective
 * arrival process at each node follows from the network characterization of
 * Theorem 10 (npfqn_traffic_rqt), the service variability parameter from the
 * adaptation of Section 7.1 (qsys_gigk_rqt_gamma), and the system time at each
 * node is the worst-case bound of Theorem 3 (qsys_gigk_rqt). The published step
 * 3, path enumeration, is not needed here: LINE aggregates per-node system times
 * into per-class response times through the visit ratios.
 *
 * The adaptation is regressed against simulation in heavy traffic, so accuracy
 * degrades at low utilization: on M/M/1 the error is about 5% at rho=0.9 but over
 * 50% at rho=0.5.
 *
 * REFERENCE INDEXING, REPRODUCED. Like solver_qna and solver_rqna, the reference
 * indexes the stateful-indexed sn.rt with STATION indices, which is exact only
 * when every stateful node is a station; that precondition is asserted.
 *
 * ARITHMETIC: transcendental (real exponents), so under Rational the body is
 * discarded and RQT is refused by name, matching solver_qna and solver_rqna.
 *
 * Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
 * Operations Research 63(3), 676-700.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/npfqn/npfqn_traffic_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt_gamma.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mva {

template <class T>
MvaSolution<T> solver_rqt(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_rqt: robust queueing theory needs transcendental arithmetic (real exponents); "
            "rerun with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    // One predicate for the gate and the run: list_valid_methods asks the same
    // question before it advertises 'rqt', so the sentence a caller reads here is
    // the sentence that kept the row off the report.
    {
        const std::string rqt_reason = mva_single_class_open_reason(L, "rqt");
        if (!rqt_reason.empty()) throw UnsupportedError(rqt_reason);
    }
    const std::size_t M = L.nstations;
    for (std::size_t r = 0; r < L.nclasses; ++r)
        if (std::isfinite(L.classes[r].population))
            throw UnsupportedError("solver_rqt: RQT supports open networks only (no closed classes)");
    if (L.nof_stateful() != M)
        throw UnsupportedError(
            "solver_rqt: the reference indexes the stateful-indexed sn.rt with station indices, "
            "which is only correct when every stateful node is a station");

    MvaSolution<T> s;
    s.Q = Matrix<T>(M, 1, zero);
    s.U = Matrix<T>(M, 1, zero);
    s.R = Matrix<T>(M, 1, zero);
    s.Tp = Matrix<T>(M, 1, zero);
    s.X.assign(1, zero);
    s.C.assign(1, zero);
    s.method = "rqt";
    s.iter = 1;

    // Configuration: adaptation regime, exact worst case, tail coefficients.
    std::string regime = opt.rqt_regime.empty() ? std::string("independent") : opt.rqt_regime;
    const bool useExact = opt.rqt_exact;
    const T alpha_a_src = num_traits<T>::from_double(opt.rqt_alpha_a > 0 ? opt.rqt_alpha_a : 2.0);
    const T alpha_s_cfg = num_traits<T>::from_double(opt.rqt_alpha_s > 0 ? opt.rqt_alpha_s : 2.0);

    // Source and queueing (incl. delay) stations.
    std::size_t src = M;
    std::vector<bool> isSource(M, false), schedInf(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        isSource[i] = (L.stations[i].nodetype == qn::NodeType::Source);
        schedInf[i] = (L.stations[i].sched == qn::SchedStrategy::INF);
        if (isSource[i] && src == M) src = i;
    }
    if (src == M)
        throw UnsupportedError("solver_rqt: RQT requires an open network with a Source station");
    std::vector<std::size_t> qstat;
    for (std::size_t i = 0; i < M; ++i)
        if (!isSource[i]) qstat.push_back(i);
    const std::size_t nq = qstat.size();

    auto rtS = [&](std::size_t i, std::size_t j) -> T { return L.rt(i, j); };

    // External arrival process from the Source.
    const mam::Map<T> arvMAP = lang::dist_to_map(L.service[src][0]);
    const T lambda_src = mam::map_lambda(arvMAP);
    const T sigma_a_src = qsys::detail::num_sqrt(mam::map_scv(arvMAP)) / lambda_src;

    // Per-node primitives.
    std::vector<T> mu(nq, zero), sigma_s(nq, zero), lambda0(nq, zero), Gamma0(nq, zero);
    std::vector<T> alpha0(nq, alpha_a_src), alpha_s(nq, alpha_s_cfg);
    std::vector<std::size_t> nserv(nq, 1);
    Matrix<T> F(nq, nq, zero);
    for (std::size_t a = 0; a < nq; ++a) {
        const std::size_t ia = qstat[a];
        mu[a] = L.rates(ia, 0);
        sigma_s[a] = qsys::detail::num_sqrt(L.scv(ia, 0)) / mu[a];
        const double ns = num_traits<T>::to_double(L.stations[ia].nservers);
        if (std::isfinite(ns) && ns > 0) nserv[a] = static_cast<std::size_t>(ns);
        // the source stream reaches node a thinned by q, Theorem 6
        const T q = rtS(src, ia);
        lambda0[a] = T(lambda_src * q);
        if (q > zero)
            Gamma0[a] = sigma_a_src * qsys::detail::num_pow(T(one / q), T(one / alpha_a_src));
        for (std::size_t b = 0; b < nq; ++b) F(a, b) = rtS(ia, qstat[b]);
    }

    // Effective arrival processes, Theorem 10.
    const npfqn::TrafficRqt<T> eff = npfqn::npfqn_traffic_rqt(lambda0, Gamma0, alpha0, F);

    // Per-node worst-case analysis.
    for (std::size_t a = 0; a < nq; ++a) {
        const std::size_t ia = qstat[a];
        const T lam = eff.lambda[a];
        s.Tp(ia, 0) = lam;
        if (lam <= zero) continue;
        if (schedInf[ia]) {
            s.U(ia, 0) = T(lam / mu[a]);
            s.Q(ia, 0) = T(lam / mu[a]);
            s.R(ia, 0) = T(one / mu[a]);
            continue;
        }
        const T kk = num_traits<T>::from_int(static_cast<int>(nserv[a]));
        const T rho = lam / (kk * mu[a]);
        const T Gamma_s = qsys::qsys_gigk_rqt_gamma(rho, mu[a], eff.Gamma[a], sigma_s[a], nserv[a],
                                                    eff.alpha[a], regime);
        const qsys::GigkRqtResult<T> w =
            qsys::qsys_gigk_rqt(lam, mu[a], eff.Gamma[a], Gamma_s, nserv[a], eff.alpha[a], alpha_s[a]);
        const T Ra = useExact ? w.Sworst : w.W;
        s.R(ia, 0) = Ra;
        s.U(ia, 0) = rho;
        s.Q(ia, 0) = T(lam * Ra);  // Little's law, number in system
    }

    s.Tp(src, 0) = lambda_src;
    T csum = zero;
    for (std::size_t i = 0; i < M; ++i) csum = T(csum + s.R(i, 0));
    s.C[0] = csum;
    s.X[0] = lambda_src;
    (void)two;
    return s;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_RQT_H
