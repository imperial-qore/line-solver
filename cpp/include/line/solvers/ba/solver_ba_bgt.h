/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_BGT_H
#define LINE_SOLVERS_BA_SOLVER_BA_BGT_H

/**
 * Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
 * multitype open Markovian network, valid for EVERY work-conserving Markovian
 * policy.
 *
 * Templated port of matlab/src/solvers/BA/solver_ba_bgt_analyzer.m,
 * cross-checked against
 * jar/src/main/java/jline/solvers/ba/analyzers/Solver_ba_bgt_analyzer.java.
 * The polyhedron and the bound are `npfqn_bnd_bgt`; this analyzer maps the LINE
 * model onto them and reads the bound back per station and class.
 *
 * CLASS SPACE. The reference's network is a MULTITYPE one: each type follows a
 * FIXED sequence of stages, and stage k of type i is its own buffer. LINE's
 * (station, job class) pair is that buffer, so the analyzer walks the routing
 * matrix from the Source and turns each open class into one type whose stages
 * are the pairs it visits. Two gates follow and are enforced by name rather
 * than approximated: routing must be DETERMINISTIC (a pair sends everything to
 * one successor, or everything to the Sink) and routes must NOT MERGE (a pair
 * belongs to exactly one type, else the reference's class index (i,k) is not
 * defined). A re-entrant line is expressible by giving the revisits distinct
 * LINE classes.
 *
 * THE BOUND IS LOOSE, and knowingly so: the exception parameter of the smoothed
 * Lyapunov function carries (Lmax+gamma)^3/gamma^2 and dominates as soon as
 * there is more than one station. What is sharp is the STABILITY CERTIFICATE --
 * a feasible gamma > 0 proves every work-conserving policy stable, and the LP
 * correctly refuses the Lu-Kumar network at per-station loads of 0.7, where
 * global stability genuinely fails -- and the geometric tail RATE.
 *
 * ARITHMETIC. Rational-clean, like the API core it calls.
 *
 * Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance of
 * multiclass Markovian queueing networks via piecewise linear Lyapunov
 * functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_bnd_bgt.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ba {

namespace detail {

/**
 * The single successor of a routing row, as a pair index, or `npairs` when
 * everything leaves the network. A probabilistic split is refused by name: the
 * reference's network has deterministic routing, and a split is a different
 * model rather than an approximation of this one.
 */
template <class T>
std::size_t bgt_single_successor(const Matrix<T>& rtst, std::size_t row,
                                 const std::vector<std::size_t>& pairFlat,
                                 const std::string& who) {
    const T tol = num_traits<T>::from_double(1e-9);
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T mass = zero, best = zero;
    std::size_t p = pairFlat.size();
    for (std::size_t q = 0; q < pairFlat.size(); ++q) {
        const T v = rtst(row, pairFlat[q]);
        if (v > tol) {
            mass = T(mass + v);
            if (v > best) {
                best = v;
                p = q;
            }
        }
    }
    if (!(mass > tol)) return pairFlat.size();
    const T dm = T(mass - one), db = T(best - one);
    const bool ok = (dm < tol && T(zero - dm) < tol) && (db < tol && T(zero - db) < tol);
    if (!ok)
        throw UnsupportedError(
            "solver_ba_bgt: method 'bgt.upper' needs deterministic routing: " + who +
            " splits its departures. The reference's network routes each type along a fixed "
            "sequence of stages");
    return p;
}

}  // namespace detail

/**
 * @param L   the model
 * @param out the (Q,U,R,Tp,C,X) block to fill; shapes are set here
 */
template <class T, class Solution>
void solver_ba_bgt(const qn::NetworkStruct<T>& L, Solution& out) {
    const T zero = num_traits<T>::from_int(0);
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
                "solver_ba_bgt: method 'bgt.upper' supports fully open networks only "
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
            "solver_ba_bgt: method 'bgt.upper' requires an open network with a Source station");
    for (std::size_t a = 0; a < qstat.size(); ++a) {
        const std::size_t i = qstat[a];
        if (L.stations[i].sched == lang::SchedStrategy::INF)
            throw UnsupportedError(
                "solver_ba_bgt: method 'bgt.upper' does not support delay (infinite-server) "
                "stations: the reference's network has one server per station");
        const double ns = num_traits<T>::to_double(L.stations[i].nservers);
        if (std::isfinite(ns) && ns > 1)
            throw UnsupportedError(
                "solver_ba_bgt: method 'bgt.upper' does not support multi-server stations");
    }

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

    // ---- walk one deterministic route per source class ----
    std::vector<T> lambda;
    std::vector<std::vector<std::size_t> > routes;
    std::vector<bool> used(np, false);
    for (std::size_t si = 0; si < srcList.size(); ++si) {
        const std::size_t s = srcList[si];
        for (std::size_t r0 = 0; r0 < K; ++r0) {
            const T arr = L.rates(s, r0);
            const double ad = num_traits<T>::to_double(arr);
            if (!std::isfinite(ad) || !(arr > zero)) continue;
            std::size_t cur = detail::bgt_single_successor(
                rtst, s * K + r0, pairFlat, "the Source for class " + std::to_string(r0 + 1));
            std::vector<std::size_t> route;
            while (cur < np) {
                if (used[cur])
                    throw UnsupportedError(
                        "solver_ba_bgt: method 'bgt.upper' needs routes that do not merge: "
                        "station " +
                        std::to_string(pairStation[cur] + 1) + " class " +
                        std::to_string(pairClass[cur] + 1) +
                        " is visited by more than one type. Give the visits distinct job classes");
                used[cur] = true;
                route.push_back(cur);
                cur = detail::bgt_single_successor(
                    rtst, pairFlat[route.back()], pairFlat,
                    "station " + std::to_string(pairStation[route.back()] + 1) + " class " +
                        std::to_string(pairClass[route.back()] + 1));
            }
            if (route.empty())
                throw UnsupportedError("solver_ba_bgt: class " + std::to_string(r0 + 1) +
                                       " leaves the Source and reaches no station");
            routes.push_back(route);
            lambda.push_back(arr);
        }
    }
    if (routes.empty()) throw UnsupportedError("solver_ba_bgt: the model carries no open traffic");

    const std::size_t I = routes.size();
    std::vector<std::vector<T> > mu(I);
    std::vector<std::vector<std::size_t> > sigma(I);
    std::vector<std::size_t> ustat;
    for (std::size_t i = 0; i < I; ++i) {
        for (std::size_t k = 0; k < routes[i].size(); ++k) {
            const std::size_t q = routes[i][k];
            const T m = L.rates(pairStation[q], pairClass[q]);
            const double md = num_traits<T>::to_double(m);
            if (!std::isfinite(md) || !(m > zero))
                throw UnsupportedError("solver_ba_bgt: station " +
                                       std::to_string(pairStation[q] + 1) +
                                       " has no service rate for class " +
                                       std::to_string(pairClass[q] + 1) +
                                       " but carries its traffic");
            if (L.procid(pairStation[q] + 1, pairClass[q] + 1) != lang::ProcessType::EXP)
                throw UnsupportedError(
                    "solver_ba_bgt: method 'bgt.upper' requires exponential service: station " +
                    std::to_string(pairStation[q] + 1) + " class " +
                    std::to_string(pairClass[q] + 1) + " is not exponential");
            mu[i].push_back(m);
            bool seen = false;
            for (std::size_t u = 0; u < ustat.size(); ++u)
                if (ustat[u] == pairStation[q]) seen = true;
            if (!seen) ustat.push_back(pairStation[q]);
        }
    }
    // Dense station index space for the LP.
    for (std::size_t i = 0; i < I; ++i) {
        sigma[i].assign(routes[i].size(), 0);
        for (std::size_t k = 0; k < routes[i].size(); ++k)
            for (std::size_t u = 0; u < ustat.size(); ++u)
                if (ustat[u] == pairStation[routes[i][k]]) sigma[i][k] = u;
    }

    const npfqn::BndBgt<T> info = npfqn::npfqn_bnd_bgt(lambda, mu, sigma, ustat.size());

    // ---- read the bound back per station and class ----
    for (std::size_t i = 0; i < I; ++i) {
        for (std::size_t k = 0; k < routes[i].size(); ++k) {
            const std::size_t q = routes[i][k];
            const std::size_t ist = pairStation[q], r = pairClass[q];
            out.Q(ist, r) = T(out.Q(ist, r) + info.Qub[i][k]);
            out.Tp(ist, r) = T(out.Tp(ist, r) + lambda[i]);
            out.U(ist, r) = T(out.U(ist, r) + lambda[i] / L.rates(ist, r));
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            if (out.Tp(i, r) > zero) out.R(i, r) = T(out.Q(i, r) / out.Tp(i, r));

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

#endif  // LINE_SOLVERS_BA_SOLVER_BA_BGT_H
