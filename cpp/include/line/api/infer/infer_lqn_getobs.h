/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_GETOBS_H
#define LINE_API_INFER_INFER_LQN_GETOBS_H

/**
 * Observation vector of a solved LQN, z = h(a).
 *
 * Templated port of matlab/src/api/infer/infer_lqn_getobs.m. The JAR carries
 * the same selection inline inside jline/api/infer/InferLqn.java.
 *
 * This is the measurement half of the LQN parameter identification method of
 * Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying Parameters in
 * Software Systems with Extended Kalman Filters", CASCON 2005: given the
 * per-element metric vectors of a solved model, it selects the entries the
 * filter observes, in the order the specification lists them. It is the map
 * whose sensitivity infer_lqn_jacobian differentiates and whose residual
 * infer_lqn_ekf corrects on.
 *
 * The metric name is an enum here rather than MATLAB's case-insensitive
 * string, so an unknown metric is a compile error instead of a runtime one;
 * the only runtime rejection left is the one MATLAB also makes, an element
 * name absent from the LQN. That rejection is an InputError and NOT a
 * not-a-number or a zero: an observation silently read off the wrong element
 * would be indistinguishable from a badly fitting model.
 *
 * ARITHMETIC: selection and copying only, no operation on the values at all,
 * so the port is exact in every instantiation.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/infer/infer_lqn_findbyname.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace infer {

/** The four per-element metrics of a solved LQN, as in getEnsembleAvg. */
enum class LqnMetric { QLen, Util, RespT, Tput };

/**
 * Per-element metric vectors, each aligned with the element name list.
 * Mirrors the struct MATLAB passes as METRICS.
 */
template <class T>
struct LqnMetrics {
    std::vector<T> QLen;
    std::vector<T> Util;
    std::vector<T> RespT;
    std::vector<T> Tput;
};

/** One row of MATLAB's OBSSPEC struct array. */
struct LqnObsSpec {
    LqnMetric metric;
    std::string name;
};

/**
 * @param names   element names, as in LayeredNetworkStruct.names
 * @param metrics per-element metric vectors, indexed like names
 * @param spec    observations to extract, in order
 * @return        (numel(spec)) observation vector
 */
template <class T>
std::vector<T> infer_lqn_getobs(const std::vector<std::string>& names,
                                const LqnMetrics<T>& metrics,
                                const std::vector<LqnObsSpec>& spec) {
    std::vector<T> z;
    z.reserve(spec.size());
    for (std::size_t i = 0; i < spec.size(); ++i) {
        const std::size_t idx = infer_lqn_findbyname(names, spec[i].name);
        if (idx == INFER_LQN_NOT_FOUND)
            throw InputError("infer_lqn_getobs: element '" + spec[i].name +
                             "' not found in the LQN");
        const std::vector<T>* v = nullptr;
        switch (spec[i].metric) {
            case LqnMetric::QLen:
                v = &metrics.QLen;
                break;
            case LqnMetric::Util:
                v = &metrics.Util;
                break;
            case LqnMetric::RespT:
                v = &metrics.RespT;
                break;
            case LqnMetric::Tput:
                v = &metrics.Tput;
                break;
        }
        if (idx >= v->size())
            throw InputError("infer_lqn_getobs: metric vector shorter than the element list");
        z.push_back((*v)[idx]);
    }
    return z;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_GETOBS_H
