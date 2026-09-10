/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_COMPUTE_QL_AT_ARRIVAL_H
#define LINE_API_INFER_INFER_COMPUTE_QL_AT_ARRIVAL_H

/**
 * Per-class queue lengths seen by each arriving job, reconstructed from
 * arrival and response time samples.
 *
 * Templated port of matlab/src/api/infer/infer_compute_ql_at_arrival.m (the
 * shipped .mexa64 is a compiled copy of that same .m). No JAR counterpart.
 *
 * Arrivals and response times are matched by job id, so the two sample sets
 * need not be in the same order nor come from the same source. Each job is
 * turned into an arrival event at t and a departure event at t + rt; the
 * events are replayed in time order, departures before arrivals at a tie, and
 * the state recorded at each arrival is the queue seen by that job INCLUDING
 * itself (the arriving job is counted first, then the state is read).
 *
 * ARITHMETIC: only comparisons and one addition per job, so a finite field
 * computation, exact in the exact instantiation. Exactness matters here: with
 * rational timestamps a tie between a departure and an arrival is decided by
 * the documented rule rather than by rounding.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * @param at       (n) arrival times
 * @param at_jobid (n) job id of each arrival
 * @param rt       (m) response times, m >= n
 * @param rt_jobid (m) job id of each response time
 * @param cls      (n) class of each arrival, 0-based in [0,R)
 * @param R        number of classes
 * @return         (n x R) queue length at each arrival, in the input order
 */
template <class T>
Matrix<T> infer_compute_ql_at_arrival(const std::vector<T>& at, const std::vector<long>& at_jobid,
                                      const std::vector<T>& rt, const std::vector<long>& rt_jobid,
                                      const std::vector<std::size_t>& cls, std::size_t R) {
    const std::size_t n = at.size();
    if (at_jobid.size() != n || cls.size() != n)
        throw InputError("infer_compute_ql_at_arrival: arrival arrays have different lengths");
    if (rt.size() != rt_jobid.size())
        throw InputError("infer_compute_ql_at_arrival: response-time arrays have different lengths");
    for (std::size_t c : cls)
        if (c >= R) throw InputError("infer_compute_ql_at_arrival: class index out of range");

    // match response times to arrivals by job id
    std::vector<T> rt_matched(n);
    for (std::size_t i = 0; i < n; ++i) {
        bool found = false;
        for (std::size_t k = 0; k < rt_jobid.size(); ++k)
            if (rt_jobid[k] == at_jobid[i]) {
                rt_matched[i] = rt[k];
                found = true;
                break;
            }
        if (!found)
            throw InputError(
                "infer_compute_ql_at_arrival: not all arrival job ids found among the response "
                "time job ids");
    }

    // stable sort of the arrivals by time, as MATLAB's sort is stable
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return at[a] < at[b]; });

    struct Event {
        T time;
        int type;  // -1 departure, +1 arrival
        std::size_t slot;
        std::size_t cls;
    };
    std::vector<Event> events;
    events.reserve(2 * n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = order[k];
        Event a;
        a.time = at[i];
        a.type = 1;
        a.slot = k;
        a.cls = cls[i];
        events.push_back(a);
    }
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = order[k];
        Event d;
        d.time = at[i] + rt_matched[i];
        d.type = -1;
        d.slot = k;
        d.cls = cls[i];
        events.push_back(d);
    }
    // sort by time, departures before arrivals at a tie
    std::stable_sort(events.begin(), events.end(), [](const Event& x, const Event& y) {
        if (x.time < y.time) return true;
        if (y.time < x.time) return false;
        return x.type < y.type;
    });

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> state(R, zero);
    Matrix<T> ql_sorted(n, R, zero);
    for (const Event& e : events) {
        if (e.type == 1) {
            state[e.cls] += one;
            for (std::size_t c = 0; c < R; ++c) ql_sorted(e.slot, c) = state[c];
        } else {
            state[e.cls] -= one;
        }
    }

    Matrix<T> ql(n, R, zero);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t c = 0; c < R; ++c) ql(order[k], c) = ql_sorted(k, c);
    return ql;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_COMPUTE_QL_AT_ARRIVAL_H
