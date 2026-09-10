/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_DT_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_DT_H

/**
 * Port of `solver_mam_dt.m`: discrete-time (slotted) analysis of an open
 * network whose interarrival and service laws all live on the slot lattice.
 *
 * A single queueing station is solved EXACTLY by the Q-MAM discrete-time
 * queues, `q_dt_ph_ph_1` when both laws are renewal discrete phase-type and
 * `q_dt_map_map_1` when either side is a D-MAP. Several stations are solved by
 * a discrete-time parametric decomposition, which is an approximation: the
 * departure process is truncated at a finite level and compressed back to a
 * bounded phase dimension, and correlations between the streams entering a
 * station are not preserved.
 *
 * Time is measured in slots internally and converted back on exit, so Q and U
 * are dimensionless, T is per time unit and R is in time units.
 *
 * Convention: late arrival system with delayed access (LAS-DA), matching the
 * Q-MAM discrete-time queues and the LDES slotted engine. `A_1 = kron(C1,D0)`
 * is the claim that a job arriving at the end of a slot is not served in it.
 *
 * ARITHMETIC: field with transcendentals. The queue solve iterates to a
 * tolerance, so the analyzer stands down under exact arithmetic exactly as the
 * MAP/MAP/1 fast path does.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/dtime.h"
#include "line/api/sn/sn_is_discrete_time.h"
#include "line/api/sn/sn_predicates.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"

namespace line {
namespace mam {

namespace dt_detail {

/** Rejects the model features the discrete-time path cannot represent. */
template <class T>
void assert_scope(const qn::NetworkStruct<T>& sn) {
    using lang::SchedStrategy;
    if (!api::sn_is_open_model(sn))
        throw UnsupportedError(
            "SolverMAM: the discrete-time path supports open models only. A closed slotted model "
            "needs a level-dependent discrete chain, which the Q-MAM discrete-time catalogue does "
            "not cover.");
    if (sn.nclasses > 1)
        throw UnsupportedError(
            "SolverMAM: the discrete-time path supports one class only. Independent per-class "
            "lattice sources fire in the same slot with positive probability, and a batch of "
            "simultaneous arrivals of different classes is not an MMAP[K], which is what "
            "Q_DT_MMAPK_PHK_1 consumes.");
    for (std::size_t ist = 0; ist < sn.nstations; ++ist) {
        const SchedStrategy sched = sn.stations[ist].sched;
        if (sched == SchedStrategy::EXT) continue;
        if (sched != SchedStrategy::FCFS)
            throw UnsupportedError(
                "SolverMAM: the discrete-time path supports FCFS single-server stations and the "
                "Source only.");
        if (sn.stations[ist].nservers > 1)
            throw UnsupportedError(
                "SolverMAM: the discrete-time path models one server per station; a slotted "
                "multiserver queue needs the level-dependent boundary of Geo/Geo/c.");
    }
}

/** Discrete-time law of a station, expressed in slots. */
template <class T>
Dmap<T> station_law(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r,
                    double slot_length) {
    using lang::ProcessType;
    const ProcessType pt = sn.procid(ist + 1, r + 1);
    if (pt == ProcessType::DMAP) {
        if (slot_length != 1.0)
            throw InputError(
                "SolverMAM: a DMAP is defined on its own slot, so it cannot be combined with a "
                "slotlength other than one.");
        Dmap<T> d;
        d.D0 = sn.service[ist][r].D0;
        d.D1 = sn.service[ist][r].D1;
        return d;
    }
    const T mean_slots =
        num_traits<T>::from_double(1.0 / (num_traits<T>::to_double(sn.rates(ist, r)) * slot_length));
    return dph_to_dmap(dph_from_dist<T>(pt, mean_slots, sn.scv(ist, r)));
}

/** Exact single-station analysis through the Q-MAM discrete-time queues. */
template <class T>
std::vector<T> solve_single_station(const Dmap<T>& arv, const Dmap<T>& svc,
                                    std::size_t max_num_comp) {
    if (dmap_is_renewal(arv) && dmap_is_renewal(svc))
        return q_dt_ph_ph_1(dmap_to_dph(arv), dmap_to_dph(svc), max_num_comp);
    return q_dt_map_map_1(arv, svc, max_num_comp);
}

/**
 * Station-to-station routing probabilities, source first and sink stripped.
 *
 * The Sink feeds back into the Source to keep the routing matrix stochastic; an
 * open traffic equation must not see that edge. Routers and class switches
 * carry no service, so their transit is censored into (I - Pnn)^-1.
 */
template <class T>
Matrix<T> routing(const qn::NetworkStruct<T>& sn, std::size_t source_idx,
                  const std::vector<std::size_t>& queue_idx) {
    using lang::NodeType;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t I = sn.nof_nodes(), K = sn.nclasses;
    Matrix<T> Pn(I, I, zero);
    for (std::size_t a = 0; a < I; ++a)
        for (std::size_t b = 0; b < I; ++b) Pn(a, b) = sn.rtnodes(a * K, b * K);

    std::vector<std::size_t> sink_nodes;
    for (std::size_t ind = 0; ind < I; ++ind)
        if (sn.nodes[ind].nodetype == NodeType::Sink) {
            sink_nodes.push_back(ind);
            for (std::size_t b = 0; b < I; ++b) Pn(ind, b) = zero;
        }

    std::vector<std::size_t> station_nodes;
    station_nodes.push_back(sn.station_to_node[source_idx] - 1);
    for (std::size_t idx = 0; idx < queue_idx.size(); ++idx)
        station_nodes.push_back(sn.station_to_node[queue_idx[idx]] - 1);

    std::vector<std::size_t> inter_nodes;
    for (std::size_t ind = 0; ind < I; ++ind) {
        const bool is_station =
            std::find(station_nodes.begin(), station_nodes.end(), ind) != station_nodes.end();
        const bool is_sink = std::find(sink_nodes.begin(), sink_nodes.end(), ind) != sink_nodes.end();
        if (!is_station && !is_sink) inter_nodes.push_back(ind);
    }

    auto sub = [&](const std::vector<std::size_t>& rows, const std::vector<std::size_t>& cols) {
        Matrix<T> out(rows.size(), cols.size(), zero);
        for (std::size_t a = 0; a < rows.size(); ++a)
            for (std::size_t b = 0; b < cols.size(); ++b) out(a, b) = Pn(rows[a], cols[b]);
        return out;
    };

    const std::size_t ns = station_nodes.size();
    Matrix<T> full = sub(station_nodes, station_nodes);
    if (!inter_nodes.empty()) {
        Matrix<T> Psn = sub(station_nodes, inter_nodes);
        Matrix<T> Pnn = sub(inter_nodes, inter_nodes);
        Matrix<T> Pns = sub(inter_nodes, station_nodes);
        Matrix<T> Inn = detail::dt_eye<T>(inter_nodes.size());
        Matrix<T> inv = detail::dt_left_solve(detail::dt_sub(Inn, Pnn), Inn);
        full = detail::dt_add(full, detail::dt_mul(detail::dt_mul(Psn, inv), Pns));
    }
    // column 0 is the Source, which receives nothing
    Matrix<T> P(ns, ns - 1, zero);
    for (std::size_t a = 0; a < ns; ++a)
        for (std::size_t b = 1; b < ns; ++b) P(a, b - 1) = full(a, b);
    return P;
}

/** Arrival stream of a queue: source share plus thinned upstream departures. */
template <class T>
DBatch<T> arrivals(const DBatch<T>& src, const std::vector<DBatch<T>>& dep, const Matrix<T>& P,
                   std::size_t idx, std::size_t space_max) {
    DBatch<T> arv;
    bool have = false;
    if (num_traits<T>::to_double(P(0, idx)) > 0) {
        arv = dmap_thin(src, P(0, idx));
        have = true;
    }
    for (std::size_t j = 0; j < dep.size(); ++j) {
        const T p = P(1 + j, idx);
        if (num_traits<T>::to_double(p) <= 0) continue;
        DBatch<T> contrib = dmap_thin(dep[j], p);
        if (!have) {
            arv = contrib;
            have = true;
        } else {
            arv = dmap_super(arv, contrib);
            arv = dmap_compress_batch(arv, space_max);
        }
    }
    if (!have)
        throw InputError(
            "SolverMAM: a queue receives no arrivals in the discrete-time routing matrix.");
    return arv;
}

/** Total visit ratio of a station across the chains. */
template <class T>
T visit_ratio(const qn::NetworkStruct<T>& sn, std::size_t ist) {
    T v = num_traits<T>::from_int(0);
    for (std::size_t c = 0; c < sn.visits.size(); ++c) {
        const Matrix<T>& vm = sn.visits[c];
        if (ist >= vm.rows()) continue;
        for (std::size_t r = 0; r < vm.cols(); ++r) v = v + vm(ist, r);
    }
    return v;
}

}  // namespace dt_detail

/**
 * Discrete-time analysis of `sn`, returning the same metric tuple as every
 * other MAM analyzer. `slot_length` is the slot in model time units.
 */
template <class T>
MamSolution<T> solver_mam_dt(const qn::NetworkStruct<T>& sn, const MamOptions& opt,
                             double slot_length) {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0);
    dt_detail::assert_scope(sn);

    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<Dmap<T>> law(M);
    for (std::size_t ist = 0; ist < M; ++ist) law[ist] = dt_detail::station_law(sn, ist, 0, slot_length);

    std::size_t source_idx = M;
    std::vector<std::size_t> queue_idx;
    for (std::size_t ist = 0; ist < M; ++ist) {
        if (sn.stations[ist].sched == SchedStrategy::EXT)
            source_idx = ist;
        else
            queue_idx.push_back(ist);
    }
    if (source_idx == M)
        throw UnsupportedError("SolverMAM: the discrete-time path requires an open model with a Source.");

    const std::size_t nq = queue_idx.size();
    std::vector<T> QN(nq, zero), UN(nq, zero), TN(nq, zero);
    std::string method;
    int totiter = 1;

    DBatch<T> src_batch;
    src_batch.push_back(law[source_idx].D0);
    src_batch.push_back(law[source_idx].D1);
    const T lambda_slot = dmap_lambda_batch(src_batch);

    if (nq == 1) {
        const std::vector<T> ql =
            dt_detail::solve_single_station(law[source_idx], law[queue_idx[0]], 1000);
        for (std::size_t i = 0; i < ql.size(); ++i)
            QN[0] = QN[0] + num_traits<T>::from_int(static_cast<int>(i)) * ql[i];
        UN[0] = num_traits<T>::from_int(1) - ql[0];
        TN[0] = lambda_slot;
        method = "dt.qmam";
    } else {
        // Initial departure streams: Bernoulli of the exact station throughput,
        // the slotted counterpart of seeding a decomposition with Poisson streams
        Matrix<T> P = dt_detail::routing(sn, source_idx, queue_idx);
        std::vector<DBatch<T>> dep(nq);
        for (std::size_t idx = 0; idx < nq; ++idx) {
            const T p = lambda_slot * dt_detail::visit_ratio(sn, queue_idx[idx]);
            DBatch<T> b;
            b.push_back(Matrix<T>(1, 1, num_traits<T>::from_int(1) - p));
            b.push_back(Matrix<T>(1, 1, p));
            dep[idx] = b;
        }

        std::vector<T> QNprev(nq, zero), QNprev2(nq, zero), UNprev(nq, zero), TNprev(nq, zero);
        const int iter_max = opt.iter_max > 0 ? opt.iter_max : 100;
        const double iter_tol = opt.tol > 0 ? opt.tol : 1e-3;
        auto max_rel = [&](const std::vector<T>& a, const std::vector<T>& b) {
            double worst = 0;
            for (std::size_t i = 0; i < a.size(); ++i) {
                const double d = std::max(std::abs(num_traits<T>::to_double(b[i])), 1e-14);
                worst = std::max(worst, std::abs(num_traits<T>::to_double(a[i] - b[i])) / d);
            }
            return worst;
        };

        for (int it = 1; it <= iter_max; ++it) {
            totiter = it;
            for (std::size_t idx = 0; idx < nq; ++idx) {
                DBatch<T> arv = dt_detail::arrivals(src_batch, dep, P, idx, opt.space_max);
                DtQueueResult<T> r = mg1_dt_queue(arv, law[queue_idx[idx]], 1000, true);
                QN[idx] = r.QN;
                UN[idx] = r.UN;
                TN[idx] = r.TN;
                DBatch<T> d;
                d.push_back(r.dep.D0);
                d.push_back(r.dep.D1);
                dep[idx] = dmap_compress_batch(d, opt.space_max);
            }
            if (it > 1 && max_rel(QN, QNprev) < iter_tol) break;
            if (it > 2 && max_rel(QN, QNprev2) < iter_tol) {
                // Feedback loops settle into a period-two cycle rather than a
                // point: re-solving a station with the departure process it just
                // produced moves it back. The cycle amplitude sits far below the
                // error of the decomposition itself, so the midpoint is reported
                // instead of burning iter_max sweeps on an orbit that will not close.
                const T half = num_traits<T>::from_double(0.5);
                for (std::size_t idx = 0; idx < nq; ++idx) {
                    QN[idx] = (QN[idx] + QNprev[idx]) * half;
                    UN[idx] = (UN[idx] + UNprev[idx]) * half;
                    TN[idx] = (TN[idx] + TNprev[idx]) * half;
                }
                break;
            }
            QNprev2 = QNprev;
            QNprev = QN;
            UNprev = UN;
            TNprev = TN;
        }
        method = "dt.dec";
    }

    MamSolution<T> out;
    mva::MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);

    const T slot = num_traits<T>::from_double(slot_length);
    for (std::size_t idx = 0; idx < nq; ++idx) {
        const std::size_t ist = queue_idx[idx];
        s.Q(ist, 0) = QN[idx];
        s.U(ist, 0) = UN[idx];
        s.Tp(ist, 0) = TN[idx] / slot;
        if (num_traits<T>::to_double(TN[idx]) > 0) s.R(ist, 0) = QN[idx] / TN[idx] * slot;
    }
    const T lambda = lambda_slot / slot;
    s.Tp(source_idx, 0) = lambda;
    s.X[0] = lambda;
    T qtot = zero;
    for (std::size_t idx = 0; idx < nq; ++idx) qtot = qtot + QN[idx];
    s.C[0] = qtot / lambda;
    s.iter = totiter;
    out.actualmethod = method;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_DT_H
