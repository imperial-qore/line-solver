/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_QUICK_MODEL_H
#define LINE_API_INFER_INFER_QUICK_MODEL_H

/**
 * @file
 * @ingroup api_infer
 * Convenience factory for the single-layer networks the inference estimators fit.
 *
 * Port of matlab/src/api/infer/infer_quick_model.m and
 * python/line_solver/inference/api/infer_quick_model.py.
 *
 * The model is named `quickModel` and has one Queue `QueueStation<i>` per entry
 * of `stations`, with that discipline, `servers[i]` servers and infinite
 * capacity, and one class `Class<c>` per row of `demands`, served at station i
 * by `Exp.fitMean(demands(c,i))`.
 *
 *   - OPEN: a Source `mySource` and a Sink `mySink`, every class routed serially
 *     Source -> QueueStation1 -> ... -> QueueStationM -> Sink. No arrival
 *     process is set, as in both references: the estimators supply it.
 *   - CLOSED: class c has `jobs[c]` jobs referenced at QueueStation1 and follows
 *     `routing[c]` (an M x M station-to-station matrix) when given, and the
 *     cycle 1 -> 2 -> ... -> M -> 1 otherwise. `routing` is ignored when open.
 *
 * The open branch follows the Python reference. The MATLAB open branch stored
 * the Source in `node{1}` and then overwrote it with the first queue, so the
 * serial path started at a queue and the Source was never routed.
 */

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * @brief Build a simple open or closed queueing network.
 *
 * @param is_open  true for an open network, false for a closed one
 * @param stations (M) scheduling discipline of each queue
 * @param demands  (K x M) mean service demand of class c at station i
 * @param servers  (M) server count per station; empty means one each
 * @param jobs     (K) population per closed class; empty means one each
 * @param routing  (K) per-class M x M station routing for a closed model; empty means serial
 * @return         the linked network
 */
template <class T>
qn::Network<T> infer_quick_model(bool is_open, const std::vector<lang::SchedStrategy>& stations,
                                 const Matrix<T>& demands,
                                 const std::vector<double>& servers = std::vector<double>(),
                                 const std::vector<double>& jobs = std::vector<double>(),
                                 const std::vector<Matrix<T>>& routing = std::vector<Matrix<T>>()) {
    const std::size_t M = stations.size();
    const std::size_t K = demands.rows();
    if (M == 0) throw InputError("infer_quick_model: at least one station is required");
    if (K == 0) throw InputError("infer_quick_model: at least one class is required");
    if (demands.cols() != M)
        throw InputError("infer_quick_model: demands must have one column per station");
    if (!servers.empty() && servers.size() != M)
        throw InputError("infer_quick_model: servers must have one entry per station");
    if (!jobs.empty() && jobs.size() != K)
        throw InputError("infer_quick_model: jobs must have one entry per class");
    if (!is_open && !routing.empty()) {
        if (routing.size() != K)
            throw InputError("infer_quick_model: routing must have one matrix per class");
        for (const Matrix<T>& P : routing)
            if (P.rows() != M || P.cols() != M)
                throw InputError("infer_quick_model: each routing matrix must be M x M");
    }

    qn::Network<T> model("quickModel");
    std::size_t source = 0, sink = 0;
    if (is_open) {
        source = model.add_source("mySource");
        sink = model.add_sink("mySink");
    }

    std::vector<std::size_t> queue(M);
    for (std::size_t i = 0; i < M; ++i) {
        queue[i] = model.add_queue("QueueStation" + std::to_string(i + 1), stations[i]);
        model.set_number_of_servers(queue[i], servers.empty() ? 1.0 : servers[i]);
        model.set_capacity(queue[i], std::numeric_limits<double>::infinity());
    }

    std::vector<std::size_t> cls(K);
    for (std::size_t c = 0; c < K; ++c) {
        const std::string nm = "Class" + std::to_string(c + 1);
        cls[c] = is_open ? model.add_open_class(nm)
                         : model.add_closed_class(nm, jobs.empty() ? 1.0 : jobs[c], queue[0]);
        for (std::size_t i = 0; i < M; ++i)
            model.set_service(queue[i], cls[c], lang::Distrib<T>::exp_mean(demands(c, i)));
    }

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    qn::RoutingMatrix<T> P = model.init_routing_matrix();
    for (std::size_t c = 0; c < K; ++c) {
        const std::size_t r = cls[c];
        if (is_open) {
            P.set(r, r, source, queue[0], one);
            for (std::size_t i = 0; i + 1 < M; ++i) P.set(r, r, queue[i], queue[i + 1], one);
            P.set(r, r, queue[M - 1], sink, one);
        } else if (!routing.empty()) {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j)
                    if (routing[c](i, j) > zero) P.set(r, r, queue[i], queue[j], routing[c](i, j));
        } else {
            for (std::size_t i = 0; i < M; ++i) P.set(r, r, queue[i], queue[(i + 1) % M], one);
        }
    }
    model.link(P);
    return model;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_QUICK_MODEL_H
