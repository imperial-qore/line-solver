/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_DT_H
#define LINE_SOLVERS_NC_SOLVER_NC_DT_H

/**
 * Exact normalizing-constant analysis of a discrete-time (slotted) model.
 *
 * Port of matlab/src/solvers/NC/nc_is_dt_model.m and
 * solver_nc_dt_analyzer.m. The route is requested with
 * `NcSolverOptions::slotted`, the same switch SolverLDES uses to run on a
 * discrete time scale, and is never auto-detected: a Geometric service time is
 * a perfectly ordinary continuous-time model unless the caller says the model
 * lives on a slot lattice.
 *
 * Two families are covered, both from Daduna (2001):
 *
 *   chapter 2  a Bernoulli server fed by a Bernoulli arrival stream, with an
 *              unbounded buffer (theorem 2.3, corollary 2.7), a finite buffer
 *              (corollary 2.8) or a load-dependent service probability
 *              (example 2.10), evaluated by `dqsys_bernoulli1`;
 *   chapter 3  a closed cycle of Bernoulli servers (theorem 3.2, corollary
 *              3.4), evaluated by `dpfqn_nc` when the service probabilities
 *              are state independent and by `dpfqn_ncld` otherwise.
 *
 * THE ADMISSIBLE FEATURE SET IS NARROW BECAUSE THE PRODUCT FORM IS NARROW.
 * Beyond the geometric service requirement:
 *
 *   - a cycle is the only topology; section 4.1 of the reference records that
 *     general discrete-time topologies of FCFS Bernoulli servers have no
 *     product form;
 *   - every station must be a single server. Pestien and Ramakrishnan, quoted
 *     before example 2.10, proved that a multiserver node inside a cycle of
 *     geometrical queues destroys the product form for ANY finite server
 *     count, so the multiserver case is refused rather than approximated;
 *   - class switching is rejected, and a multichain cycle is admitted only
 *     through its aggregate population.
 *
 * A model that fails any of these is an ERROR, not a fallback to the
 * continuous-time analyzer: the continuous-time answer to a slotted question is
 * a different number, not a worse one.
 *
 * Every metric is on the slot lattice: a rate is a per-slot probability and a
 * time is a number of slots. `NcSolverOptions::slotlength` rescales both to
 * model time units.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/dpfqn/dpfqn_nc.h"
#include "line/api/dqsys/dqsys_bernoulli1.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"

namespace line {
namespace nc {

/** Classification of a model against the discrete-time product form. */
template <class T>
struct DtModel {
    std::string kind = "none";  ///< "bernoulli1", "cycle" or "none"
    std::string reason;         ///< why not, when kind is "none"
    std::size_t station = 0;    ///< 0-based queueing station (bernoulli1)
    std::size_t source = 0;     ///< 0-based Source station (bernoulli1)
    T arrival_prob;             ///< offered b (bernoulli1)
    std::vector<T> service_single;  ///< p(n), n = 1..L (bernoulli1)
    std::size_t capacity = 0;   ///< 0 means unbounded
    std::vector<std::size_t> order;      ///< stations in cycle order
    std::vector<std::vector<T> > service;  ///< p_j(n), [cycle position][n-1]
    std::size_t population = 0;
};

namespace detail {

/** Load-dependent scaling of station `ist` expanded to alpha(1..N). */
template <class T>
std::vector<T> dt_lld_vector(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t N) {
    const T one = num_traits<T>::from_int(1);
    std::vector<T> alpha(N, one);
    const std::vector<T>& lld = sn.stations[ist].lldscaling;
    if (lld.empty()) return alpha;
    bool varies = false;
    for (std::size_t k = 0; k < lld.size(); ++k) {
        if (!(lld[k] == one)) { varies = true; break; }
    }
    if (!varies) return alpha;
    // Hold the last declared value beyond the tabulated range, as the
    // load-dependent normalizing-constant analyzers do.
    for (std::size_t n = 0; n < N; ++n) {
        alpha[n] = (n < lld.size()) ? lld[n] : lld[lld.size() - 1];
    }
    return alpha;
}

/** True when station `ist` declares a load dependence. */
template <class T>
bool dt_has_lld(const qn::NetworkStruct<T>& sn, std::size_t ist) {
    const T one = num_traits<T>::from_int(1);
    const std::vector<T>& lld = sn.stations[ist].lldscaling;
    for (std::size_t k = 0; k < lld.size(); ++k) {
        if (!(lld[k] == one)) return true;
    }
    return false;
}

/** Effective buffer capacity of station `ist`, 0 when unbounded. */
template <class T>
std::size_t dt_capacity(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r) {
    double cap = sn.stations[ist].cap;
    if (!sn.classcap.empty() && ist < sn.classcap.size() && r < sn.classcap[ist].size()) {
        const double cc = sn.classcap[ist][r];
        if (cc < cap) cap = cc;
    }
    if (!std::isfinite(cap) || cap <= 0) return 0;
    return static_cast<std::size_t>(cap);
}

/**
 * Station order along the cycle, or empty when the routing is not a single
 * deterministic cycle visiting every station exactly once.
 */
template <class T>
std::vector<std::size_t> dt_cycle_order(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, R = sn.nclasses;
    const api::SnRtStations<T> rt = api::sn_rt_stations(sn);
    const Matrix<T>& rtst = rt.rtst;
    const double tol = 1e-8, coarse = 1e-3;

    std::vector<std::size_t> succ(M, 0);
    for (std::size_t i = 0; i < M; ++i) {
        long tgt = -1;
        for (std::size_t j = 0; j < M; ++j) {
            double w = 0;
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t s = 0; s < R; ++s) {
                    w += num_traits<T>::to_double(rtst(i * R + r, j * R + s));
                }
            }
            if (w > tol) {
                if (std::fabs(w - static_cast<double>(R)) > coarse && std::fabs(w - 1.0) > coarse) {
                    return std::vector<std::size_t>();  // fractional routing out of i
                }
                if (tgt >= 0) return std::vector<std::size_t>();  // more than one successor
                tgt = static_cast<long>(j);
            }
        }
        if (tgt < 0 || static_cast<std::size_t>(tgt) == i) return std::vector<std::size_t>();
        succ[i] = static_cast<std::size_t>(tgt);
    }

    std::vector<bool> visited(M, false);
    std::vector<std::size_t> order;
    std::size_t cur = 0;
    for (std::size_t k = 0; k < M; ++k) {
        if (visited[cur]) return std::vector<std::size_t>();
        visited[cur] = true;
        order.push_back(cur);
        cur = succ[cur];
    }
    if (cur != 0) return std::vector<std::size_t>();
    for (std::size_t i = 0; i < M; ++i) {
        if (!visited[i]) return std::vector<std::size_t>();
    }
    return order;
}

}  // namespace detail

/** Classify `sn` against the two discrete-time product-form families. */
template <class T>
DtModel<T> nc_is_dt_model(const qn::NetworkStruct<T>& sn) {
    DtModel<T> dt;
    const std::size_t M = sn.nstations, R = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) {
            const T rate = sn.rates(i, r);
            if (rate > zero) {
                if (sn.procid(i + 1, r + 1) != qn::ProcessType::GEOMETRIC) {
                    dt.reason = "a station serves a class with a non-Geometric process; a "
                                "discrete-time model needs Geometric service and interarrival times";
                    return dt;
                }
            }
        }
    }
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].cdscaling || sn.stations[i].jdscaling) {
            dt.reason = "class- or joint-dependent scaling is not covered by the discrete-time "
                        "product form";
            return dt;
        }
    }

    const std::vector<double> njobs = sn.njobs();
    bool is_open = false;
    for (std::size_t r = 0; r < njobs.size(); ++r) {
        if (std::isinf(njobs[r])) { is_open = true; break; }
    }

    if (is_open) {
        if (R != 1) {
            dt.reason = "the discrete-time single-node route handles one open class";
            return dt;
        }
        std::size_t src = 0, nsrc = 0, ist = 0, nq = 0;
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.stations[i].sched == qn::SchedStrategy::EXT) { src = i; ++nsrc; }
            else { ist = i; ++nq; }
        }
        if (nsrc != 1) {
            dt.reason = "an open discrete-time model needs exactly one Source";
            return dt;
        }
        if (nq != 1) {
            dt.reason = "the discrete-time single-node route handles one queueing station";
            return dt;
        }
        if (sn.stations[ist].sched != qn::SchedStrategy::FCFS) {
            dt.reason = "a Bernoulli server is a FCFS station";
            return dt;
        }
        if (std::isfinite(sn.stations[ist].nservers) && sn.stations[ist].nservers != 1.0) {
            dt.reason = "a Bernoulli server is a single-server station; use load dependence for "
                        "the multiserver approximation of example 2.10";
            return dt;
        }
        const T b = sn.rates(src, 0), p = sn.rates(ist, 0);
        if (!(b > zero) || b > one) {
            dt.reason = "the source arrival probability must lie in (0,1]";
            return dt;
        }
        if (!(p > zero) || p > one) {
            dt.reason = "the service probability must lie in (0,1]";
            return dt;
        }
        const std::size_t cap = detail::dt_capacity(sn, ist, 0);
        const bool has_lld = detail::dt_has_lld(sn, ist);
        if (cap == 0 && has_lld) {
            dt.reason = "a load-dependent Bernoulli server needs a finite capacity to bound the "
                        "state space";
            return dt;
        }
        if (cap == 0 && !(b < p)) {
            dt.reason = "an unbounded discrete-time queue needs an arrival probability below the "
                        "service probability";
            return dt;
        }
        if (cap == 0) {
            dt.service_single.assign(1, p);
        } else {
            const std::vector<T> alpha = detail::dt_lld_vector(sn, ist, cap);
            dt.service_single.resize(cap);
            for (std::size_t n = 0; n < cap; ++n) {
                dt.service_single[n] = p * alpha[n];
                if (!(dt.service_single[n] > zero) || dt.service_single[n] > one) {
                    dt.reason = "load dependence must keep the service probability inside (0,1]";
                    return dt;
                }
            }
        }
        dt.kind = "bernoulli1";
        dt.station = ist;
        dt.source = src;
        dt.arrival_prob = b;
        dt.capacity = cap;
        return dt;
    }

    double Ntot = 0;
    for (std::size_t r = 0; r < njobs.size(); ++r) {
        if (std::isfinite(njobs[r])) Ntot += njobs[r];
    }
    if (Ntot <= 0 || Ntot != std::floor(Ntot)) {
        dt.reason = "the closed population must be a positive integer";
        return dt;
    }
    const std::size_t N = static_cast<std::size_t>(Ntot);

    std::vector<T> p(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched != qn::SchedStrategy::FCFS) {
            dt.reason = "a cycle of Bernoulli servers is FCFS throughout";
            return dt;
        }
        if (std::isfinite(sn.stations[i].nservers) && sn.stations[i].nservers != 1.0) {
            dt.reason = "a multiserver node inside a cycle of geometrical queues has no product "
                        "form";
            return dt;
        }
        bool seen = false;
        T lo = zero, hi = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T v = sn.rates(i, r);
            if (v > zero) {
                if (!seen) { lo = v; hi = v; seen = true; }
                else { if (v < lo) lo = v; if (v > hi) hi = v; }
            }
        }
        if (!seen) {
            dt.reason = "a station of the cycle serves no class";
            return dt;
        }
        if (num_traits<T>::to_double(hi - lo) > 1e-8 * num_traits<T>::to_double(hi)) {
            dt.reason = "a station has a class-dependent service probability; the discrete-time "
                        "cycle needs one Bernoulli server per node";
            return dt;
        }
        p[i] = lo;
        if (!(p[i] > zero) || !(p[i] < one)) {
            dt.reason = "the product form of theorem 3.2 needs every service probability in (0,1)";
            return dt;
        }
    }

    const std::vector<std::size_t> order = detail::dt_cycle_order(sn);
    if (order.empty()) {
        dt.reason = "the stations do not form a single deterministic cycle; discrete-time FCFS "
                    "networks of other topologies have no product form";
        return dt;
    }

    dt.service.assign(M, std::vector<T>(N, zero));
    for (std::size_t k = 0; k < M; ++k) {
        const std::vector<T> alpha = detail::dt_lld_vector(sn, order[k], N);
        for (std::size_t n = 0; n < N; ++n) {
            dt.service[k][n] = p[order[k]] * alpha[n];
            if (!(dt.service[k][n] > zero) || dt.service[k][n] > one) {
                dt.reason = "load dependence must keep every service probability inside (0,1]";
                return dt;
            }
        }
    }
    dt.kind = "cycle";
    dt.order = order;
    dt.population = N;
    return dt;
}

/**
 * Exact discrete-time analysis of `sn`.
 *
 * On the cycle route the per-class split is proportional to the per-class
 * population. Service in the cycle is type independent and FCFS forbids
 * overtaking, so the cyclic order of the jobs is frozen; the marginal law of
 * the queue lengths carries no class information, and the long-run share of
 * station j held by chain g is its population share N_g/N. That is the sense in
 * which section 3.2 of the reference calls the multichain case a direct
 * adaptation of the unichain one.
 */
template <class T>
NcSolution<T> solver_nc_dt(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    const DtModel<T> dt = nc_is_dt_model(sn);
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    NcSolution<T> out;
    out.sol.Q = Matrix<T>(M, K, zero);
    out.sol.U = Matrix<T>(M, K, zero);
    out.sol.R = Matrix<T>(M, K, zero);
    out.sol.Tp = Matrix<T>(M, K, zero);
    out.sol.C.assign(K, zero);
    out.sol.X.assign(K, zero);
    out.sol.iter = 1;

    if (dt.kind == "bernoulli1") {
        out.sol.method = "dt.bernoulli1";
        std::vector<T> b(1, dt.arrival_prob);
        dqsys::Bernoulli1Result<T> r;
        if (dt.capacity == 0) {
            // An unbounded buffer is the finite chain taken far enough out that
            // the geometric tail is below the working precision; the closed form
            // of corollary 2.7 agrees with it term for term.
            const double ratio = num_traits<T>::to_double(dt.arrival_prob)
                    * (1.0 - num_traits<T>::to_double(dt.service_single[0]))
                    / (num_traits<T>::to_double(dt.service_single[0])
                       * (1.0 - num_traits<T>::to_double(dt.arrival_prob)));
            std::size_t L = 64;
            if (ratio > 0 && ratio < 1) {
                const std::size_t need = static_cast<std::size_t>(std::ceil(std::log(1e-18)
                        / std::log(ratio)));
                if (need > L) L = need;
            }
            if (L > 100000) L = 100000;
            r = dqsys::dqsys_bernoulli1(b, dt.service_single, L);
        } else {
            r = dqsys::dqsys_bernoulli1(b, dt.service_single, dt.capacity);
        }
        out.sol.Q(dt.station, 0) = r.meanQueueLength;
        out.sol.U(dt.station, 0) = r.utilization;
        out.sol.Tp(dt.station, 0) = r.throughput;
        out.sol.R(dt.station, 0) = r.meanSojournTime;
        out.sol.X[0] = r.throughput;
        out.sol.C[0] = r.meanSojournTime;
        // The Source row carries the offered stream, as on the continuous-time route.
        out.sol.Tp(dt.source, 0) = dt.arrival_prob;
        out.sol.lG = std::log(num_traits<T>::to_double(r.normConst));
    } else if (dt.kind == "cycle") {
        const std::size_t N = dt.population;
        bool state_independent = true;
        for (std::size_t k = 0; k < dt.service.size() && state_independent; ++k) {
            for (std::size_t n = 1; n < dt.service[k].size(); ++n) {
                if (!(dt.service[k][n] == dt.service[k][0])) { state_independent = false; break; }
            }
        }
        std::vector<T> Qs(M, zero), Us(M, zero), Ts(M, zero);
        if (state_independent) {
            // Propositions 3.18 and 3.19 end to end, with the index of corollary
            // 3.20(a) corrected (see dpfqn_nc.h).
            out.sol.method = "dt.cycle";
            std::vector<T> p(M);
            for (std::size_t k = 0; k < M; ++k) p[k] = dt.service[k][0];
            const dpfqn::DtNcResult<T> nc = dpfqn::dpfqn_nc(p, N);
            const T x = nc.throughput();
            for (std::size_t k = 0; k < M; ++k) {
                const T q = one - p[k];
                T ratio = one, tail = zero;
                for (std::size_t n = 1; n <= N; ++n) {
                    ratio = ratio * q / p[k];
                    tail = tail + ratio / q * nc.G1[N - n + 1] / nc.G;
                }
                Qs[dt.order[k]] = tail;      // E[X_j] = sum_{n>=1} P(X_j>=n)
                Us[dt.order[k]] = x / p[k];  // P(X_j >= 1)
                Ts[dt.order[k]] = x;
            }
            out.sol.lG = nc.lG;
        } else {
            out.sol.method = "dt.cycleld";
            const dpfqn::DtNcLdResult<T> nc = dpfqn::dpfqn_ncld(dt.service, N);
            for (std::size_t k = 0; k < M; ++k) {
                const std::vector<T> marg = nc.marginal(k);
                T q = zero, t = zero;
                for (std::size_t n = 0; n <= N; ++n) {
                    q = q + marg[n] * num_traits<T>::from_int(static_cast<long>(n));
                    if (n >= 1) t = t + marg[n] * dt.service[k][n - 1];
                }
                Qs[dt.order[k]] = q;
                Us[dt.order[k]] = one - marg[0];
                Ts[dt.order[k]] = t;
            }
            out.sol.lG = nc.lG;
        }

        const std::vector<double> njobs = sn.njobs();
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < K; ++r) {
                const T share = num_traits<T>::from_double(
                        std::isfinite(njobs[r]) ? njobs[r] / static_cast<double>(N) : 0.0);
                out.sol.Q(i, r) = Qs[i] * share;
                out.sol.U(i, r) = Us[i] * share;
                out.sol.Tp(i, r) = Ts[i] * share;
                if (out.sol.Tp(i, r) > zero) {
                    out.sol.R(i, r) = out.sol.Q(i, r) / out.sol.Tp(i, r);
                }
            }
        }
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t ref = sn.classes[r].refstat;
            if (ref >= 1 && ref <= M) out.sol.X[r] = out.sol.Tp(ref - 1, r);
            if (out.sol.X[r] > zero && std::isfinite(njobs[r])) {
                out.sol.C[r] = num_traits<T>::from_double(njobs[r]) / out.sol.X[r];
            }
        }
    } else {
        throw InputError("solver_nc_dt: the slotted option was requested but the model is not a "
                         "discrete-time product-form model: " + dt.reason);
    }

    out.actualmethod = out.sol.method;
    if (opt.slotlength != 1.0) {
        const T d = num_traits<T>::from_double(opt.slotlength);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < K; ++r) {
                out.sol.Tp(i, r) = out.sol.Tp(i, r) / d;
                out.sol.R(i, r) = out.sol.R(i, r) * d;
            }
        }
        for (std::size_t r = 0; r < K; ++r) {
            out.sol.X[r] = out.sol.X[r] / d;
            out.sol.C[r] = out.sol.C[r] * d;
        }
    }
    return out;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_DT_H
