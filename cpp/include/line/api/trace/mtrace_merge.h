/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_MTRACE_MERGE_H
#define LINE_API_TRACE_MTRACE_MERGE_H

/**
 * Superposes two single-class traces into one marked trace, labelling the
 * events of the first stream 1 and those of the second 2.
 *
 * Templated port of matlab/lib/m3a/m3a/mtrace/mtrace_merge.m, cross-checked
 * against jar/src/main/java/jline/api/trace/Mtrace_merge.java (identical; the
 * JAR carries an explicit comment that the single shared time origin
 * `sort([0; cumsum(t1); cumsum(t2)])` matters, since merging the two
 * cumulative-sum vectors with separate origins injects a spurious zero-length
 * interval).
 *
 * Ties are broken in favour of the first stream, which is MATLAB's stable
 * sort and is reproduced with std::stable_sort.
 *
 * ARITHMETIC: a merge and pairwise differences, exact in Rational.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of mtrace_merge, mirroring [T, L]. */
template <class T>
struct MtraceMergeResult {
    std::vector<T> times;   ///< inter-arrival times of the merged process
    std::vector<int> labels;  ///< 1 for the first stream, 2 for the second
};

/**
 * @param t1 inter-arrival times of the first trace
 * @param t2 inter-arrival times of the second trace
 */
template <class T>
MtraceMergeResult<T> mtrace_merge(const std::vector<T>& t1, const std::vector<T>& t2) {
    if (t1.empty() && t2.empty()) throw InputError("mtrace_merge: both traces are empty");
    struct Event {
        T time;
        int label;
    };
    std::vector<Event> ev;
    ev.reserve(t1.size() + t2.size() + 1);
    Event origin;
    origin.time = num_traits<T>::from_int(0);
    origin.label = 0;
    ev.push_back(origin);
    T acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < t1.size(); ++i) {
        acc += t1[i];
        Event e;
        e.time = acc;
        e.label = 1;
        ev.push_back(e);
    }
    acc = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < t2.size(); ++i) {
        acc += t2[i];
        Event e;
        e.time = acc;
        e.label = 2;
        ev.push_back(e);
    }
    std::stable_sort(ev.begin(), ev.end(),
                     [](const Event& a, const Event& b) { return a.time < b.time; });

    MtraceMergeResult<T> out;
    out.times.reserve(ev.size() - 1);
    out.labels.reserve(ev.size() - 1);
    for (std::size_t k = 1; k < ev.size(); ++k) {
        out.times.push_back(ev[k].time - ev[k - 1].time);
        out.labels.push_back(ev[k].label);
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_MTRACE_MERGE_H
