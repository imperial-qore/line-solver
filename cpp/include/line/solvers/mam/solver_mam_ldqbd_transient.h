/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_TRANSIENT_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_TRANSIENT_H

/**
 * Port of `solver_mam_ldqbd_transient.m`: transient queue length, utilization
 * and throughput of a single-class OPEN queue, and the fast path behind
 * `getTranAvg`.
 *
 * TWO ENGINES, chosen by the buffer, and the split is forced rather than
 * stylistic:
 *  - FINITE capacity: the generator is a finite matrix, so the law is
 *    pi(t) = pi(0) exp(Qt), stepped on a uniform grid by ONE matrix exponential
 *    `expm(Q dt)` reused at every point. That is the reference's construction
 *    and it is exact up to `expm`.
 *  - INFINITE capacity: there is no finite generator to exponentiate, so the
 *    reference calls libQBD's adaptive Taylor series, which grows the
 *    represented level depth as mass reaches it (`api/mam/libqbd_taylor.h`).
 *    The reference grid is libQBD's own, one point per 1/|min diagonal|.
 *
 * WHAT THE PORT DOES NOT DO. The reference's time grid in the finite branch is
 * `min(101, max(11, round(10 T)))` points, which is reproduced exactly, because
 * a transient result read off a different grid cannot be compared point for
 * point with the reference at all.
 *
 * PH SERVICE IS SINGLE-SERVER ONLY, in both branches, and the reference says so
 * (`Transient QBD with PH service supports single-server only`): the level
 * phase would have to carry the multiset of in-service phases. Refused by name.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/libqbd_taylor.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** One station-class transient curve, the reference's `[metric, time]` pair. */
template <class T>
struct TranCurve {
    std::vector<T> values;
    std::vector<double> times;
};

/** What `getTranAvg` returns: queue length, utilization and throughput curves. */
template <class T>
struct TranResult {
    /** Indexed [station][class]; only the queue station is populated. */
    std::vector<std::vector<TranCurve<T>>> Qt, Ut, Tt;
};

/**
 * Port of `mam_transient_qbd_applicable.m`: true when the Laplace-domain
 * transient QBD should run instead of this fast path.
 *
 * The Laplace solver is for single-server open queues whose ARRIVAL is
 * non-Poisson or whose SERVICE is a correlated (non-renewal) MAP -- exactly
 * what the level structure here cannot represent. Poisson arrival with PH or
 * exponential service, and M/M/c, stay on the fast path.
 */
template <class T>
bool mam_transient_qbd_applicable(const qn::NetworkStruct<T>& L) {
    using lang::SchedStrategy;
    if (L.nclasses != 1) return false;
    if (!std::isinf(L.classes[0].population)) return false;
    std::size_t src = 0, q = 0, nsrc = 0, nq = 0;
    for (std::size_t i = 1; i <= L.nstations; ++i) {
        if (L.stations[i - 1].sched == SchedStrategy::EXT) { src = i; ++nsrc; }
        else if (L.stations[i - 1].sched == SchedStrategy::FCFS) { q = i; ++nq; }
    }
    if (nsrc != 1 || nq != 1) return false;
    if (L.stations[q - 1].nservers != 1.0) return false;
    const Map<T> arv = lang::dist_to_map(L.service[src - 1][0]);
    const Map<T> svc = lang::dist_to_map(L.service[q - 1][0]);
    const bool arrivalIsPoisson = (arv.D0.rows() == 1);
    const bool serviceIsRenewal = basic_detail::is_renewal_map(svc);
    return !arrivalIsPoisson || !serviceIsRenewal;
}

/**
 * Port of `solver_mam_ldqbd_transient.m`.
 *
 * @param opt `timespan` bounds the horizon; `tol` is the Taylor truncation
 *            target in the infinite branch
 * @param L the refreshed struct
 */
template <class T>
TranResult<T> solver_mam_ldqbd_transient(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_ldqbd_transient: the transient law is a matrix exponential (finite "
            "buffer) or a tolerance-truncated Taylor series with an incomplete-gamma error bound "
            "(infinite buffer); rerun with --arith double or --arith real");
    } else {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;
    if (K != 1)
        throw UnsupportedError(
            "solver_mam_ldqbd_transient: the transient QBD method requires a single-class model");
    if (!std::isinf(L.classes[0].population))
        throw UnsupportedError(
            "solver_mam_ldqbd_transient: the transient QBD method requires an open model");

    std::size_t src = 0, q = 0, nsrc = 0, nq = 0;
    for (std::size_t i = 1; i <= M; ++i) {
        if (L.stations[i - 1].sched == SchedStrategy::EXT) { src = i; ++nsrc; }
        else if (L.stations[i - 1].sched == SchedStrategy::FCFS) { q = i; ++nq; }
    }
    if (nsrc != 1 || nq != 1)
        throw UnsupportedError(
            "solver_mam_ldqbd_transient: the transient QBD method requires exactly one Source and "
            "one FCFS Queue");

    const T lambda = L.rates(src - 1, 0);
    const Map<T> PHq = lang::dist_to_map(L.service[q - 1][0]);
    const double nServers = L.stations[q - 1].nservers;
    const double bufCap = L.cap[q - 1];
    const std::size_t nPhases = PHq.D0.rows();
    const bool isPH = nPhases > 1;
    T mu = zero;
    std::vector<T> alphaV;
    Matrix<T> texit;
    if (!isPH) {
        mu = T(-PHq.D0(0, 0));
    } else {
        alphaV = map_pie(PHq);
        texit = Matrix<T>(nPhases, 1, zero);
        for (std::size_t i = 0; i < nPhases; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < nPhases; ++j) s += PHq.D0(i, j);
            texit(i, 0) = -s;
        }
    }
    if (isPH && nServers > 1.0)
        throw UnsupportedError(
            "solver_mam_ldqbd_transient: transient QBD with PH service supports single-server "
            "queues only; the level phase would have to carry the multiset of in-service phases");

    const double T_start = opt.timespan_start;
    const double T_end = opt.timespan_end;
    if (!(T_end > T_start) || !std::isfinite(T_end))
        throw InputError(
            "solver_mam_ldqbd_transient: the timespan must be a finite interval with a positive "
            "duration");
    const double T_duration = T_end - T_start;
    const unsigned c = static_cast<unsigned>(
        std::isfinite(nServers) ? std::llround(nServers) : 1);

    std::vector<double> times;
    std::vector<T> qlen, util, tput;

    if (std::isfinite(bufCap)) {
        // ---- finite capacity: one expm, stepped ---------------------------
        const std::size_t Cap = static_cast<std::size_t>(std::llround(bufCap));
        std::size_t dim;
        Matrix<T> Q;
        if (!isPH) {
            dim = Cap + 1;
            Q = Matrix<T>(dim, dim, zero);
            for (std::size_t n = 0; n <= Cap; ++n) {
                const T dep = num_traits<T>::from_double(
                                  std::min(static_cast<double>(n), static_cast<double>(c))) * mu;
                const T arr = (n < Cap) ? lambda : zero;
                if (n > 0) Q(n, n - 1) = dep;
                if (n < Cap) Q(n, n + 1) = arr;
                Q(n, n) = T(-(dep + arr));
            }
        } else {
            dim = 1 + Cap * nPhases;
            Q = Matrix<T>(dim, dim, zero);
            Q(0, 0) = -lambda;
            for (std::size_t j = 0; j < nPhases; ++j) Q(0, 1 + j) = T(lambda * alphaV[j]);
            for (std::size_t n = 1; n <= Cap; ++n) {
                const std::size_t r0 = 1 + (n - 1) * nPhases;
                for (std::size_t i = 0; i < nPhases; ++i)
                    for (std::size_t j = 0; j < nPhases; ++j) {
                        Q(r0 + i, r0 + j) = PHq.D0(i, j);
                        if (i == j && n < Cap) Q(r0 + i, r0 + j) -= lambda;
                    }
                if (n < Cap)
                    for (std::size_t i = 0; i < nPhases; ++i)
                        Q(r0 + i, r0 + nPhases + i) = lambda;
                if (n == 1) {
                    for (std::size_t i = 0; i < nPhases; ++i) Q(r0 + i, 0) = texit(i, 0);
                } else {
                    const std::size_t p0 = r0 - nPhases;
                    for (std::size_t i = 0; i < nPhases; ++i)
                        for (std::size_t j = 0; j < nPhases; ++j) Q(r0 + i, p0 + j) = PHq.D1(i, j);
                }
            }
        }
        const std::size_t nT = static_cast<std::size_t>(std::min<double>(
            101.0, std::max<double>(11.0, std::round(T_duration * 10.0))));
        const double dt = T_duration / static_cast<double>(nT - 1);
        const Matrix<T> eQdt = expm(Q, num_traits<T>::from_double(dt));

        std::vector<T> pi(dim, zero);
        pi[0] = one;
        for (std::size_t t = 0; t < nT; ++t) {
            times.push_back(T_start + static_cast<double>(t) * dt);
            T qv = zero, uv = zero, tv = zero;
            for (std::size_t n = 0; n <= Cap; ++n) {
                T pn = zero;
                if (!isPH) {
                    pn = pi[n];
                } else if (n == 0) {
                    pn = pi[0];
                } else {
                    for (std::size_t i = 0; i < nPhases; ++i) pn += pi[1 + (n - 1) * nPhases + i];
                }
                qv += num_traits<T>::from_int(static_cast<int>(n)) * pn;
                if (n >= 1) {
                    uv += num_traits<T>::from_double(
                              std::min(static_cast<double>(n), static_cast<double>(c)) /
                              static_cast<double>(c)) * pn;
                    if (!isPH) {
                        tv += num_traits<T>::from_double(
                                  std::min(static_cast<double>(n), static_cast<double>(c))) * mu *
                              pn;
                    } else {
                        for (std::size_t i = 0; i < nPhases; ++i)
                            tv += pi[1 + (n - 1) * nPhases + i] * texit(i, 0);
                    }
                }
            }
            qlen.push_back(qv);
            util.push_back(uv);
            tput.push_back(tv);
            if (t + 1 < nT) pi = vecmul(pi, eQdt);
        }
    } else {
        // ---- infinite capacity: libQBD's adaptive Taylor series -----------
        LibQbdProcess<T> proc;
        if (!isPH) {
            proc.add_zero_level(Matrix<T>(1, 1, T(-lambda)), Matrix<T>(1, 1, lambda));
            for (unsigned n = 1; n + 1 <= c; ++n) {
                const T dep = num_traits<T>::from_int(static_cast<int>(n)) * mu;
                proc.add_level(Matrix<T>(1, 1, dep), Matrix<T>(1, 1, T(-(lambda + dep))),
                               Matrix<T>(1, 1, lambda));
            }
            const T dep = num_traits<T>::from_int(static_cast<int>(c)) * mu;
            proc.add_final_level(Matrix<T>(1, 1, dep), Matrix<T>(1, 1, T(-(lambda + dep))));
        } else {
            Matrix<T> up0(1, nPhases, zero);
            for (std::size_t j = 0; j < nPhases; ++j) up0(0, j) = T(lambda * alphaV[j]);
            proc.add_zero_level(Matrix<T>(1, 1, T(-lambda)), up0);
            Matrix<T> A10(nPhases, nPhases, zero), A1p(nPhases, nPhases, zero);
            for (std::size_t i = 0; i < nPhases; ++i) {
                for (std::size_t j = 0; j < nPhases; ++j) A10(i, j) = PHq.D0(i, j);
                A10(i, i) -= lambda;
                A1p(i, i) = lambda;
            }
            proc.add_level(texit, A10, A1p);
            proc.add_final_level(PHq.D1, A10);
        }
        std::vector<std::vector<T>> pi0(1, std::vector<T>(1, one));
        const TaylorSeriesResult<T> ts =
            taylor_series_adaptive(proc, pi0, opt.tol, T_duration);
        for (std::size_t t = 0; t < ts.times.size(); ++t) {
            times.push_back(ts.times[t] + T_start);
            const std::vector<std::vector<T>>& d = ts.dists[t];
            T qv = zero, uv = zero, tv = zero;
            for (std::size_t n = 1; n < d.size(); ++n) {
                T pn = zero;
                for (const T& v : d[n]) pn += v;
                qv += num_traits<T>::from_int(static_cast<int>(n)) * pn;
                uv += num_traits<T>::from_double(
                          std::min(static_cast<double>(n), static_cast<double>(c)) /
                          static_cast<double>(c)) * pn;
                if (!isPH) {
                    tv += num_traits<T>::from_double(
                              std::min(static_cast<double>(n), static_cast<double>(c))) * mu * pn;
                } else {
                    for (std::size_t i = 0; i < d[n].size() && i < nPhases; ++i)
                        tv += d[n][i] * texit(i, 0);
                }
            }
            qlen.push_back(qv);
            util.push_back(uv);
            tput.push_back(tv);
        }
    }

    TranResult<T> out;
    out.Qt.assign(M, std::vector<TranCurve<T>>(K));
    out.Ut.assign(M, std::vector<TranCurve<T>>(K));
    out.Tt.assign(M, std::vector<TranCurve<T>>(K));
    out.Qt[q - 1][0] = TranCurve<T>{qlen, times};
    out.Ut[q - 1][0] = TranCurve<T>{util, times};
    out.Tt[q - 1][0] = TranCurve<T>{tput, times};
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_TRANSIENT_H
