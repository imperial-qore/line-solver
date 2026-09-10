/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_OI_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_OI_H

/**
 * Exact mean-value analysis for networks with order-independent stations.
 *
 * Templated port of `matlab/src/solvers/MVA/solver_mva_oi_analyzer.m`.
 *
 * An ORDER-INDEPENDENT (OI) station is a class-dependent load-dependent server
 * whose total service rate mu(n) is a permutation-invariant function of the
 * per-class count vector. A closed network of delay stations, load-independent
 * single-server queues and ANY number of OI stations is product-form, and
 * `pfqn_mvaoi` -- the Conditional MVA that carries one rate-shift vector per OI
 * station -- solves it exactly, with no normalizing constant and no joint
 * marginal.
 *
 * References: Reiser and Lavenberg, JACM 27(2), 1980; the load-dependent
 * extension of Bruell, Balbo and Afshari, 1984; OI stations and CMVA in
 * Casale 2009 and Casale, Comte and Dorsman 2026.
 *
 * A MULTISERVER BCMP QUEUE IS PROMOTED TO AN OI STATION rather than
 * approximated. `pfqn_mvaoi` models a load-independent queue as a SINGLE
 * server, so a c-server station cannot take that path; but its BCMP weight
 * satisfies the OI balance recursion exactly for the permutation-invariant rate
 *
 *     mu(n) = (min(|n|, c) / |n|) * sum_{r : n_r > 0} n_r / D_r,
 *
 * so it goes down the OI path instead and the answer stays exact. Reporting
 * still classifies it as a queue, so its utilization keeps the per-server
 * offered-load convention.
 *
 * DETECTION mirrors `nc_is_oi_model`: PAS or OI scheduling, a service-rate
 * function, and an ALL-ZERO swap graph. A nonzero swap graph is a genuine
 * pass-and-swap station, which is not order-independent and is not product
 * form, so it is excluded from the list -- and if that leaves no OI station at
 * all, the analyzer refuses rather than solving a different model.
 *
 * Arithmetic: whatever `pfqn_mvaoi` accepts; the analyzer itself is field
 * arithmetic.
 */

#include <cmath>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_mvaoi.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mva {

/**
 * The OI stations of a model, 1-based station indices.
 *
 * Port of `find_oi_stations`. Empty when the model has none, which is what the
 * dispatch tests before entering this analyzer.
 */
template <class T>
std::vector<std::size_t> find_oi_stations(const qn::NetworkStruct<T>& L) {
    std::vector<std::size_t> oi;
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const qn::Station<T>& st = L.stations[i];
        if (st.sched != qn::SchedStrategy::PAS && st.sched != qn::SchedStrategy::OI) continue;
        if (!st.svc_rate_fun) continue;
        // The graph is read through `station_swap_graph`, which applies the
        // defaults `refreshLocalVars.m` installs: an OI station is always the
        // zero graph, and a PAS station with no explicit graph is the complete
        // compatibility graph, hence a genuine pass-and-swap station and NOT
        // order-independent. MATLAB never stores an empty graph, so the
        // `isempty(sg)` arm of `nc_is_oi_model` is unreachable there.
        if (!qn::station_swap_graph_is_zero(L, i + 1)) continue;
        oi.push_back(i + 1);
    }
    return oi;
}

/**
 * Port of `nc_is_oi_model`: whether the model as a WHOLE is order-independent,
 * which is a strictly stronger condition than having an OI station.
 *
 * The whole network must be product-form for the CMVA to be exact, so every
 * station has to be a delay, an OI station with an all-zero swap graph, a PS or
 * LCFS-PR queue, or an FCFS / SIRO queue whose service rate does NOT depend on
 * the class. A closed population is required as well: the recursion is over a
 * finite lattice.
 */
template <class T>
bool nc_is_oi_model(const qn::NetworkStruct<T>& L) {
    using qn::SchedStrategy;
    for (const auto& c : L.classes)
        if (std::isinf(c.population)) return false;
    bool hasOI = false;
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const qn::Station<T>& st = L.stations[i];
        if (st.sched == SchedStrategy::INF) continue;
        if (st.sched == SchedStrategy::PAS || st.sched == SchedStrategy::OI) {
            if (!st.svc_rate_fun) return false;
            // see the note in find_oi_stations on the swap-graph defaults
            if (!qn::station_swap_graph_is_zero(L, i + 1))
                return false;  // genuine pass-and-swap: not order-independent
            hasOI = true;
        } else if (st.sched == SchedStrategy::PS || st.sched == SchedStrategy::LCFSPR ||
                   st.sched == SchedStrategy::SIRO || st.sched == SchedStrategy::FCFS) {
            if (st.sched != SchedStrategy::FCFS && st.sched != SchedStrategy::SIRO) continue;
            // a class-dependent FCFS or SIRO rate breaks product form
            double lo = 0.0, hi = 0.0;
            bool any = false;
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                if (!(L.classes[r].population > 0.0)) continue;
                const double v = num_traits<T>::to_double(L.rates(i, r));
                if (!std::isfinite(v)) continue;
                if (!any) {
                    lo = hi = v;
                    any = true;
                } else {
                    lo = std::min(lo, v);
                    hi = std::max(hi, v);
                }
            }
            if (any && (hi - lo) > 1e-9 * hi) return false;
        } else {
            return false;  // an unsupported (non-product-form) station
        }
    }
    return hasOI;
}

namespace detail {

/**
 * The OI rate at a count vector, evaluated through the station's microstate
 * function: `repelem(1:R, n)` is the canonical 1-based ordered state with n_r
 * jobs of class r, and permutation invariance makes any other order give the
 * same rate.
 */
template <class T>
T oi_rate(const std::function<T(const std::vector<std::size_t>&)>& f, const std::vector<int>& n,
          std::size_t R) {
    int tot = 0;
    for (int v : n) tot += v;
    if (tot == 0) return num_traits<T>::from_int(0);
    std::vector<std::size_t> micro;
    micro.reserve(static_cast<std::size_t>(tot));
    for (std::size_t r = 0; r < R && r < n.size(); ++r)
        for (int a = 0; a < n[r]; ++a) micro.push_back(r + 1);
    return f(micro);
}

/**
 * The OI rate that reproduces a c-server BCMP queue with per-class demands Dq:
 * `mu(n) = (min(|n|,c)/|n|) sum_{r : n_r>0} n_r/Dq_r`. At c = 1 this is the
 * ordinary total completion rate of a multiclass single server, and for one
 * class it reduces to min(n,c)/Dq, the M/M/c rate.
 */
template <class T>
T ms_oi_rate(const std::vector<int>& n, const std::vector<T>& Dq, double c) {
    const T zero = num_traits<T>::from_int(0);
    int tot = 0;
    for (int v : n) tot += v;
    if (tot == 0) return zero;
    T acc = zero;
    for (std::size_t r = 0; r < n.size() && r < Dq.size(); ++r)
        if (n[r] > 0 && Dq[r] > zero) acc += T(num_traits<T>::from_int(n[r]) / Dq[r]);
    const double cc = (std::isfinite(c) && c > 0.0) ? c : 1.0;
    const double m = std::min<double>(static_cast<double>(tot), cc);
    return T(num_traits<T>::from_double(m / static_cast<double>(tot)) * acc);
}

}  // namespace detail

/** Port of `solver_mva_oi_analyzer.m`. */
template <class T>
MvaSolution<T> solver_mva_oi_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::SchedStrategy;
    (void)opt;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, R = L.nclasses;

    const std::vector<std::size_t> oi_list = find_oi_stations(L);
    if (oi_list.empty())
        throw UnsupportedError(
            "solver_mva_oi_analyzer: the OI solver requires at least one order-independent "
            "station (PAS or OI scheduling, a service-rate function, and an all-zero swap graph)");

    // The OI rank rates are indexed by RAW class, so a chain that merges several
    // classes has no rate to evaluate: the recursion is driven by the per-class
    // population, which class switching makes meaningless (a class appearing only
    // mid-chain carries population 0, so the OI station reads as empty).
    for (std::size_t c = 0; c < L.nchains; ++c)
        if (L.inchain[c].size() > 1)
            throw UnsupportedError(
                "solver_mva_oi: requires one class per chain (no class switching)");

    std::vector<int> N(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        const double nr = L.classes[r].population;
        if (!std::isfinite(nr))
            throw UnsupportedError(
                "solver_mva_oi_analyzer: the CMVA recursion is over a closed lattice; this model "
                "has an open class");
        N[r] = static_cast<int>(std::llround(nr));
    }

    std::vector<bool> isOI(M, false), isDelay(M, false);
    for (std::size_t i : oi_list) isOI[i - 1] = true;
    for (std::size_t i = 0; i < M; ++i) isDelay[i] = (L.stations[i].sched == SchedStrategy::INF);

    // per-class demand D = V / rate at every station
    Matrix<T> V(M, R, zero), D(M, R, zero);
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                V(i, r) = T(V(i, r) + L.visits[c](L.stateful_of_station(i + 1) - 1, r));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            const T mu = L.rates(i, r);
            if (std::isfinite(num_traits<T>::to_double(mu)) && mu > zero)
                D(i, r) = T(V(i, r) / mu);
        }

    // aggregate the delays into Z; split the rest into single-server queues and
    // multiserver queues, the latter promoted to OI stations
    std::vector<T> Z(R, zero);
    std::vector<std::size_t> li_list, ms_list;
    for (std::size_t i = 0; i < M; ++i) {
        if (isOI[i]) continue;
        if (isDelay[i]) {
            for (std::size_t r = 0; r < R; ++r) Z[r] += D(i, r);
        } else if (std::isfinite(L.stations[i].nservers) && L.stations[i].nservers > 1.0) {
            ms_list.push_back(i + 1);
        } else {
            li_list.push_back(i + 1);
        }
    }
    Matrix<T> Dli(li_list.size(), R, zero);
    for (std::size_t a = 0; a < li_list.size(); ++a)
        for (std::size_t r = 0; r < R; ++r) Dli(a, r) = D(li_list[a] - 1, r);

    std::vector<std::function<T(const std::vector<int>&)>> mu;
    mu.reserve(oi_list.size() + ms_list.size());
    for (std::size_t o = 0; o < oi_list.size(); ++o) {
        const std::function<T(const std::vector<std::size_t>&)> f =
            L.stations[oi_list[o] - 1].svc_rate_fun;
        mu.push_back([f, R](const std::vector<int>& n) { return detail::oi_rate<T>(f, n, R); });
    }
    for (std::size_t j = 0; j < ms_list.size(); ++j) {
        std::vector<T> Dq(R, zero);
        for (std::size_t r = 0; r < R; ++r) Dq[r] = D(ms_list[j] - 1, r);
        const double c = L.stations[ms_list[j] - 1].nservers;
        mu.push_back([Dq, c](const std::vector<int>& n) { return detail::ms_oi_rate<T>(n, Dq, c); });
    }

    const pfqn::MvaoiResult<T> res = pfqn::pfqn_mvaoi(Z, N, mu, Dli, /*want_soi=*/true);

    Matrix<T> QN(M, R, zero), TN(M, R, zero), RN(M, R, zero), UN(M, R, zero);
    for (std::size_t o = 0; o < oi_list.size(); ++o)
        for (std::size_t r = 0; r < R; ++r) QN(oi_list[o] - 1, r) = res.Qoi(o, r);
    for (std::size_t j = 0; j < ms_list.size(); ++j)
        for (std::size_t r = 0; r < R; ++r)
            QN(ms_list[j] - 1, r) = res.Qoi(oi_list.size() + j, r);
    for (std::size_t a = 0; a < li_list.size(); ++a)
        for (std::size_t r = 0; r < R; ++r) QN(li_list[a] - 1, r) = res.Qli(a, r);
    // a delay station takes its exact product-form share, X * D
    for (std::size_t i = 0; i < M; ++i)
        if (isDelay[i])
            for (std::size_t r = 0; r < R; ++r) QN(i, r) = T(res.X[r] * D(i, r));

    // which row of Soi / Qoi holds each OI station
    std::vector<std::size_t> oiRow(M, 0);
    for (std::size_t o = 0; o < oi_list.size(); ++o) oiRow[oi_list[o] - 1] = o + 1;

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            TN(i, r) = T(res.X[r] * V(i, r));
            if (res.X[r] > zero) RN(i, r) = T(QN(i, r) / res.X[r]);
            const double sv_raw = L.stations[i].nservers;
            const double sv = (std::isfinite(sv_raw) && sv_raw > 0.0) ? sv_raw : 1.0;
            if (isOI[i]) {
                // IN-SERVICE utilization, E[sir_r]/nservers, the convention
                // solver_nc_oi_analyzer and the exact CTMC/LDES also use: sir_r
                // counts the class-r jobs receiving a strictly positive rank
                // rate, which Soi holds.
                UN(i, r) = T(res.Soi(oiRow[i] - 1, r) / num_traits<T>::from_double(sv));
            } else if (isDelay[i]) {
                UN(i, r) = QN(i, r);
            } else {
                UN(i, r) = T(res.X[r] * D(i, r) / num_traits<T>::from_double(sv));
            }
        }

    MvaSolution<T> out;
    out.Q = QN;
    out.U = UN;
    out.R = RN;
    out.Tp = TN;
    out.X = res.X;
    // the reference reports CN as the per-station response time, not as a
    // per-class system time; it returns the RN matrix itself
    out.C.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) out.C[r] += RN(i, r);
    out.method = "oi";
    int iter = 0;
    for (int v : N) iter += v;
    out.iter = iter;
    out.lG = 0.0;
    return out;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_OI_H
