/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_BPT_H
#define LINE_SOLVERS_BA_SOLVER_BA_BPT_H

/**
 * Achievable-region LOWER bound on the mean response times of a multiclass
 * open Markovian network, valid for EVERY non-idling scheduling policy at
 * every station.
 *
 * Templated port of matlab/src/solvers/BA/solver_ba_bpt_analyzer.m,
 * cross-checked against
 * jar/src/main/java/jline/solvers/ba/analyzers/Solver_ba_bpt_analyzer.java.
 * The polyhedron is `npfqn_bnd_bpt`; this analyzer maps the LINE model onto it
 * and reads the bound back per station and class.
 *
 * CLASS SPACE. The reference's "class" is a buffer: one exponential service
 * rate, one Markovian routing law. LINE's (station, job class) pair is exactly
 * that, so a pair carrying traffic becomes one LP class, the Source is absorbed
 * into the external arrival vector, and class switching needs no special
 * treatment because sn.rt already carries it.
 *
 * BOUND CONVENTION. R(i,r) minimizes x over the polyhedron with the objective
 * set to that pair's unit vector, so each entry is a valid lower bound on its
 * own. Q follows by Little's law from the bounded R and the EXACT throughput T
 * (an open network's per-class rates are fixed by the traffic equations, not by
 * the policy), and so does C. U is exact for the same reason.
 *
 * TIGHTNESS. Exact on M/M/1 and tight on the externally fed classes, but weak
 * on a class whose arrivals are all internal: the only term coupling x_r to the
 * second-moment block carries the factor lambda0_r, so an internally fed class
 * can fall back to its own mean service time.
 *
 * ARITHMETIC. Rational-clean, like the API core it calls.
 *
 * Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994). Optimization
 * of multiclass queueing networks: polyhedral and nonlinear characterizations
 * of achievable performance. Annals of Applied Probability 4(1), 43-75.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_bnd_bpt.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ba {

/**
 * @param L   the model
 * @param out the (Q,U,R,Tp,C,X) block to fill; shapes are set here
 */
template <class T, class Solution>
void solver_ba_bpt(const qn::NetworkStruct<T>& L, Solution& out) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;

    out.Q = Matrix<T>(M, K, zero);
    out.U = Matrix<T>(M, K, zero);
    out.R = Matrix<T>(M, K, zero);
    out.Tp = Matrix<T>(M, K, zero);
    out.C.assign(K, zero);
    out.X.assign(K, zero);
    out.lG = std::numeric_limits<double>::quiet_NaN();
    out.iter = 1;

    // ---- model gates ----
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(L.classes[r].population))
            throw UnsupportedError(
                "solver_ba_bpt: method 'bpt.lower' supports fully open networks only "
                "(no closed classes)");
    std::vector<std::size_t> srcList, qstat;
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].nodetype == qn::NodeType::Source)
            srcList.push_back(i);
        else
            qstat.push_back(i);
    }
    if (srcList.empty())
        throw UnsupportedError(
            "solver_ba_bpt: method 'bpt.lower' requires an open network with a Source station");
    for (std::size_t a = 0; a < qstat.size(); ++a) {
        const std::size_t i = qstat[a];
        if (L.stations[i].sched == lang::SchedStrategy::INF)
            throw UnsupportedError(
                "solver_ba_bpt: method 'bpt.lower' does not support delay (infinite-server) "
                "stations: the achievable region is derived for one server per station");
        const double ns = num_traits<T>::to_double(L.stations[i].nservers);
        if (std::isfinite(ns) && ns > 1)
            throw UnsupportedError(
                "solver_ba_bpt: method 'bpt.lower' does not support multi-server stations");
    }

    // ---- station-space routing, with the Source absorbed into lambda0 ----
    const Matrix<T> rtst = api::sn_rt_stations(L).rtst;

    std::vector<std::size_t> pairStation, pairClass, pairFlat;
    for (std::size_t a = 0; a < qstat.size(); ++a) {
        const std::size_t i = qstat[a];
        for (std::size_t r = 0; r < K; ++r) {
            pairStation.push_back(i);
            pairClass.push_back(r);
            pairFlat.push_back(i * K + r);
        }
    }
    const std::size_t np = pairFlat.size();

    std::vector<T> lambda0(np, zero);
    for (std::size_t si = 0; si < srcList.size(); ++si) {
        const std::size_t s = srcList[si];
        for (std::size_t r0 = 0; r0 < K; ++r0) {
            const T arr = L.rates(s, r0);
            const double ad = num_traits<T>::to_double(arr);
            if (!std::isfinite(ad) || !(arr > zero)) continue;
            for (std::size_t p = 0; p < np; ++p)
                lambda0[p] = T(lambda0[p] + arr * rtst(s * K + r0, pairFlat[p]));
        }
    }

    // Flow to the Sink or back to a Source is the exit probability, i.e. the
    // row deficit, and needs no column.
    Matrix<T> P(np, np, zero);
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t q = 0; q < np; ++q) P(p, q) = rtst(pairFlat[p], pairFlat[q]);

    // ---- restrict to the pairs that actually carry traffic ----
    Matrix<T> ImPt(np, np, zero);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t j = 0; j < np; ++j) ImPt(i, j) = T((i == j ? one : zero) - P(j, i));
    Matrix<T> rhs0(np, 1, zero);
    for (std::size_t p = 0; p < np; ++p) rhs0(p, 0) = lambda0[p];
    const Matrix<T> lamAll = matmul(inverse(ImPt), rhs0);
    double lamMax = 0.0;
    for (std::size_t p = 0; p < np; ++p)
        lamMax = std::max(lamMax, num_traits<T>::to_double(lamAll(p, 0)));
    const T lamTol = num_traits<T>::from_double(1e-12 * std::max(1.0, lamMax));
    std::vector<std::size_t> keep;
    for (std::size_t p = 0; p < np; ++p)
        if (lamAll(p, 0) > lamTol) keep.push_back(p);
    if (keep.empty()) throw UnsupportedError("solver_ba_bpt: the model carries no open traffic");

    const std::size_t nk = keep.size();
    std::vector<T> lam0k(nk, zero), muk(nk, zero);
    std::vector<std::size_t> statk(nk), clsk(nk), ustat;
    Matrix<T> Pk(nk, nk, zero);
    for (std::size_t a = 0; a < nk; ++a) {
        const std::size_t p = keep[a];
        lam0k[a] = lambda0[p];
        statk[a] = pairStation[p];
        clsk[a] = pairClass[p];
        for (std::size_t b = 0; b < nk; ++b) Pk(a, b) = P(p, keep[b]);
        muk[a] = L.rates(statk[a], clsk[a]);
        const double md = num_traits<T>::to_double(muk[a]);
        if (!std::isfinite(md) || !(muk[a] > zero))
            throw UnsupportedError("solver_ba_bpt: station " + std::to_string(statk[a] + 1) +
                                   " has no service rate for class " + std::to_string(clsk[a] + 1) +
                                   " but carries its traffic");
        if (L.procid(statk[a] + 1, clsk[a] + 1) != lang::ProcessType::EXP)
            throw UnsupportedError(
                "solver_ba_bpt: method 'bpt.lower' requires exponential service: station " +
                std::to_string(statk[a] + 1) + " class " + std::to_string(clsk[a] + 1) +
                " is not exponential");
        bool seen = false;
        for (std::size_t u = 0; u < ustat.size(); ++u)
            if (ustat[u] == statk[a]) seen = true;
        if (!seen) ustat.push_back(statk[a]);
    }
    // Dense station index space for the LP.
    std::vector<std::size_t> stationOf(nk, 0);
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t u = 0; u < ustat.size(); ++u)
            if (ustat[u] == statk[a]) stationOf[a] = u;

    // ---- one LP per pair, objective = that pair's unit vector ----
    for (std::size_t a = 0; a < nk; ++a) {
        std::vector<T> e(nk, zero);
        e[a] = one;
        const npfqn::BndBpt<T> info = npfqn::npfqn_bnd_bpt(lam0k, muk, Pk, stationOf, e);
        out.R(statk[a], clsk[a]) = info.zlb;
        out.Tp(statk[a], clsk[a]) = info.lambda[a];
        out.U(statk[a], clsk[a]) = info.rho[a];
    }

    // ---- exact open-network quantities ----
    for (std::size_t si = 0; si < srcList.size(); ++si) {
        const std::size_t s = srcList[si];
        for (std::size_t r = 0; r < K; ++r) {
            const T arr = L.rates(s, r);
            const double ad = num_traits<T>::to_double(arr);
            if (std::isfinite(ad) && arr > zero) {
                out.Tp(s, r) = T(out.Tp(s, r) + arr);
                out.X[r] = T(out.X[r] + arr);
            }
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) out.Q(i, r) = T(out.Tp(i, r) * out.R(i, r));
    for (std::size_t r = 0; r < K; ++r) {
        if (out.X[r] > zero) {
            T sum = zero;
            for (std::size_t i = 0; i < M; ++i) sum = T(sum + out.Q(i, r));
            out.C[r] = T(sum / out.X[r]);
        }
    }
}

}  // namespace ba
}  // namespace line

#endif  // LINE_SOLVERS_BA_SOLVER_BA_BPT_H
