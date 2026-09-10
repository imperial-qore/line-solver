/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_SETPARAMS_H
#define LINE_API_INFER_INFER_LQN_SETPARAMS_H

/**
 * Apply a parameter vector to a LayeredNetworkStruct, and read it back.
 *
 * Templated port of matlab/src/api/infer/infer_lqn_setparams.m, plus the
 * getparams half that MATLAB keeps as a local function inside infer_lqn.m.
 * No JAR counterpart. This is the parameter-injection side of the LQN
 * parameter identification method; the filter itself is in infer_lqn_ekf.h
 * and the driver in infer_lqn.h.
 *
 * A parameter is named by an element name and a kind:
 *   HOSTDEM  the mean host demand of an ACTIVITY
 *   THINK    the mean think time of a TASK
 *
 * Means are injected as EXPONENTIAL laws (SCV 1), which is the assumption of
 * the CASCON 2005 tracking method. Injecting the mean while keeping the
 * previous shape would make the observation map depend on a shape the filter
 * never estimates, so the reference does not do it and neither does this port.
 *
 * MATLAB mutates the model objects and then invalidates model.lsn so the next
 * solve re-reads them. The C++ port operates on the STRUCT directly, which is
 * what the solvers consume, so there is no cache to invalidate; the caller is
 * responsible for handing the mutated struct back to the solver.
 *
 * The lookup is by declared name over the whole element list, so a name that
 * denotes both a task and an activity is ambiguous. The kind resolves it:
 * HOSTDEM searches ACTIVITY elements only and THINK searches TASK elements
 * only, which is exactly what MATLAB's model.activities and model.tasks lists
 * do and is stricter than a bare name search over `names`.
 *
 * ARITHMETIC: assignment only, so any T works.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace infer {

/** The two parameter kinds MATLAB's paramSpec(i).type names. */
enum class LqnParamType { HOSTDEM, THINK };

/** One row of MATLAB's PARAMSPEC struct array. */
struct LqnParamSpec {
    LqnParamType type;
    std::string name;
};

namespace detail {

/** Index of the element of the given kind and name, or 0 when absent. */
template <class T>
std::size_t lqn_find_element(const lqn::LqnStruct<T>& lsn, lang::LqnElement kind,
                             const std::string& name) {
    for (std::size_t i = 1; i < lsn.names.size(); ++i)
        if (lsn.type[i] == kind && lsn.names[i] == name) return i;
    return 0;
}

/** Index of the element a parameter row names, throwing when it is absent. */
template <class T>
std::size_t lqn_param_index(const lqn::LqnStruct<T>& lsn, const LqnParamSpec& p) {
    const lang::LqnElement kind =
        (p.type == LqnParamType::HOSTDEM) ? lang::LqnElement::ACTIVITY : lang::LqnElement::TASK;
    const std::size_t idx = lqn_find_element(lsn, kind, p.name);
    if (idx == 0) {
        if (p.type == LqnParamType::HOSTDEM)
            throw InputError("infer_lqn_setparams: activity '" + p.name + "' not found");
        throw InputError("infer_lqn_setparams: task '" + p.name + "' not found");
    }
    return idx;
}

}  // namespace detail

/**
 * Read the current values of the parameters named in `spec`.
 *
 * @param lsn  layered struct to read
 * @param spec parameters to read, in order
 * @return     (numel(spec)) current values
 */
template <class T>
std::vector<T> infer_lqn_getparams(const lqn::LqnStruct<T>& lsn,
                                   const std::vector<LqnParamSpec>& spec) {
    std::vector<T> a;
    a.reserve(spec.size());
    for (std::size_t i = 0; i < spec.size(); ++i) {
        const std::size_t idx = detail::lqn_param_index(lsn, spec[i]);
        a.push_back(spec[i].type == LqnParamType::HOSTDEM ? lsn.hostdem[idx].mean
                                                          : lsn.think[idx].mean);
    }
    return a;
}

/**
 * Set the parameters named in `spec` to the values in `a`, in place.
 *
 * @param lsn  layered struct to mutate
 * @param spec parameters to set, in order
 * @param a    (numel(spec)) values, injected as exponential means
 */
template <class T>
void infer_lqn_setparams(lqn::LqnStruct<T>& lsn, const std::vector<LqnParamSpec>& spec,
                         const std::vector<T>& a) {
    if (a.size() != spec.size())
        throw InputError("infer_lqn_setparams: length of parameter vector does not match spec");
    for (std::size_t i = 0; i < spec.size(); ++i) {
        const std::size_t idx = detail::lqn_param_index(lsn, spec[i]);
        const lang::Distrib<T> d = lang::Distrib<T>::exp_mean(a[i]);
        if (spec[i].type == LqnParamType::HOSTDEM)
            lsn.hostdem[idx] = d;
        else
            lsn.think[idx] = d;
    }
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_SETPARAMS_H
