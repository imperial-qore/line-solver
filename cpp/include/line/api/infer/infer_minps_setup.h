/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_MINPS_SETUP_H
#define LINE_API_INFER_INFER_MINPS_SETUP_H

/**
 * Turn a raw per-class trace into the sample set MINPS estimates from.
 *
 * Port of matlab/src/api/infer/infer_minps_setup.m, which is MATLAB-ONLY --
 * there is no JAR or native-Python twin.
 *
 * WHAT IT DOES, in the reference's order, because each step changes what the
 * estimator sees:
 *
 *  1. DROP the classes with no samples. A class that never appears carries no
 *     information, and leaving it in shifts every later class's index, so the
 *     drop is a RELABELLING and the caller is told the surviving order.
 *  2. Derive the per-class queue lengths at arrival (`infer_get_qlen_arrival`).
 *  3. Merge the classes into one stream SORTED BY ARRIVAL TIME. That is the
 *     whole point of the sort: the estimator conditions on the state a job
 *     found, so the samples have to be in the order the system produced them,
 *     not grouped by class.
 *  4. Take the contiguous window `[initSample, initSample+sampleSize)`. A
 *     window, not a random subset, again because the state a job finds is only
 *     meaningful within a contiguous stretch of the trace.
 *  5. Drop samples whose response time is not positive; they have no
 *     phase-type density.
 *  6. Estimate the think-time rates, capped at 1e6, which is also the value a
 *     negative estimate is replaced by.
 *
 * ARITHMETIC: double, following the estimator it feeds.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/infer/infer_get_qlen_arrival.h"
#include "line/api/infer/infer_mlps.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/** One class's raw trace: when its jobs arrived, and how long they took. */
struct MinpsClassTrace {
    std::vector<double> arrival_ms;  ///< arrival times, MILLISECONDS
    std::vector<double> rt;          ///< response times, seconds
};

/** The prepared sample set, plus what the preparation had to decide. */
struct MinpsSetup {
    std::vector<MlpsSample> samples;   ///< sorted by arrival time, window applied
    std::vector<double> lambda;        ///< per-surviving-class think-time rate
    /** Original 1-based class label of each surviving class, in the new order. */
    std::vector<std::size_t> classMap;
    double threads = 0.0;              ///< Wexp, the largest total queue length seen
};

/**
 * @param traces     per-class raw trace, one entry per ORIGINAL class
 * @param initSample 1-based index of the first sample of the window
 * @param sampleSize window length; 0 means every sample
 */
inline MinpsSetup infer_minps_setup(const std::vector<MinpsClassTrace>& traces,
                                    std::size_t initSample, std::size_t sampleSize) {
    if (traces.empty()) throw InputError("infer_minps_setup: no classes");
    if (initSample < 1) throw InputError("infer_minps_setup: initSample is 1-based");

    // ---- 1. drop the classes with no samples, and record the relabelling ---
    std::vector<std::vector<double>> at_ms, rt;
    MinpsSetup out;
    for (std::size_t k = 0; k < traces.size(); ++k) {
        if (traces[k].arrival_ms.size() != traces[k].rt.size())
            throw InputError(
                "infer_minps_setup: a class's arrival and response vectors disagree in length");
        if (traces[k].rt.empty()) continue;
        at_ms.push_back(traces[k].arrival_ms);
        rt.push_back(traces[k].rt);
        out.classMap.push_back(k + 1);
    }
    const std::size_t R = at_ms.size();
    if (R == 0) throw InputError("infer_minps_setup: every class is empty");

    // ---- 2. the queue length each job found ------------------------------
    const std::vector<Matrix<double>> qls = infer::infer_get_qlen_arrival<double>(at_ms, rt);

    // ---- 3. one stream, sorted by arrival time ---------------------------
    struct Row {
        double at, rt;
        std::size_t cls;
        std::vector<double> ql;
    };
    std::vector<Row> rows;
    for (std::size_t k = 0; k < R; ++k)
        for (std::size_t i = 0; i < rt[k].size(); ++i) {
            Row r;
            r.at = at_ms[k][i] / 1000.0;  // the reference works in seconds here
            r.rt = rt[k][i];
            r.cls = k + 1;
            r.ql.assign(R, 0.0);
            for (std::size_t c = 0; c < R && c < qls[k].cols(); ++c) r.ql[c] = qls[k](i, c);
            rows.push_back(r);
        }
    std::stable_sort(rows.begin(), rows.end(),
                     [](const Row& a, const Row& b) { return a.at < b.at; });

    // Wexp and the busy-processor estimate are formed over the WHOLE stream,
    // not over the window, exactly as the reference does.
    double Wexp = 0.0, qlTotal = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        double s = 0.0;
        for (std::size_t c = 0; c < R; ++c) s += rows[i].ql[c];
        Wexp = std::max(Wexp, s);
        qlTotal += s;
    }
    out.threads = Wexp;
    double numNotProc = Wexp - qlTotal / static_cast<double>(rows.size());
    numNotProc /= static_cast<double>(R);

    // ---- 4. the contiguous window ----------------------------------------
    if (sampleSize == 0) sampleSize = rows.size();
    if (initSample - 1 + sampleSize > rows.size())
        throw InputError(
            "infer_minps_setup: the requested window runs past the end of the trace; a window is "
            "contiguous by construction and cannot be shortened silently");
    const std::size_t first = initSample - 1, last = first + sampleSize - 1;

    std::vector<std::size_t> perClass(R, 0);
    for (std::size_t i = first; i <= last; ++i) ++perClass[rows[i].cls - 1];

    // ---- 6. the think-time rates -----------------------------------------
    const double span = rows[last].at + rows[last].rt - rows[first].at;
    out.lambda.assign(R, 1e6);
    for (std::size_t k = 0; k < R; ++k) {
        double v = 1e6;
        if (span > 0.0 && numNotProc != 0.0)
            v = (static_cast<double>(perClass[k]) / span) / numNotProc;
        // A negative rate is not a slow class, it is an unusable estimate, and
        // the reference replaces it with the same cap.
        if (!(v >= 0.0)) v = 1e6;
        out.lambda[k] = std::min(1e6, v);
    }

    // ---- 5. drop the non-positive response times -------------------------
    for (std::size_t i = first; i <= last; ++i) {
        if (!(rows[i].rt > 0.0)) continue;
        MlpsSample s;
        s.rt = rows[i].rt;
        s.cls = rows[i].cls;
        s.ql = rows[i].ql;
        out.samples.push_back(s);
    }
    if (out.samples.empty())
        throw InputError("infer_minps_setup: the window holds no usable sample");
    return out;
}

/** Prepare the trace and run MINPS on it, as the reference's last line does. */
inline std::vector<double> infer_minps_from_trace(const std::vector<MinpsClassTrace>& traces,
                                                  std::size_t initSample, std::size_t sampleSize,
                                                  double nCores) {
    const MinpsSetup st = infer_minps_setup(traces, initSample, sampleSize);
    return infer_minps(st.lambda, nCores, st.samples);
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_INFER_INFER_MINPS_SETUP_H
