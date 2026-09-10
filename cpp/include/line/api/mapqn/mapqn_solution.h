/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_SOLUTION_H
#define LINE_API_MAPQN_MAPQN_SOLUTION_H

/**
 * The result of a MAPQN linear or nonlinear bound program.
 *
 * Templated port of jar/src/main/java/jline/api/mapqn/Mapqn_solution.java and
 * the native Python `MapqnSolution`. The QRF bound models write their answers
 * into a NAME-KEYED variable table rather than into a fixed struct, because the
 * variable set depends on the model that was built -- a linear-reduction bound
 * carries per-phase utilizations `U_i_k`, a queue-recursive one carries the
 * aggregate `U_i`, and the MVA-shaped variant writes `UN_i_k` / `QN_i_k`. This
 * class is the accessor layer over that table, so a caller does not have to know
 * which spelling the model it ran happens to use.
 *
 * THE `e_i_k` FALLBACK IS NOT A GUESS. Some solvers name the per-phase
 * utilization variable `e_i_k` rather than `U_i_k`; the JAR's phase-indexed
 * getter falls back to it and Python's does not, which means the same solution
 * table reads as zero in one codebase and as the utilization in the other. The
 * JAR's behaviour is kept here, since a zero utilization is indistinguishable
 * from "the model did not write one" and silently understates every bound
 * derived from it.
 *
 * A MISSING VARIABLE READS AS ZERO, in both references. That is a deliberate
 * convention of the bound models -- a variable the model did not create is a
 * quantity it does not constrain -- and not an error path.
 *
 * ARITHMETIC: field. Pure table lookup.
 */

#include <map>
#include <string>

#include "line/num/number.h"

namespace line {
namespace mapqn {

/** Objective value and the name-keyed variable table of one bound program. */
template <class T>
struct MapqnSolution {
    T objectiveValue = num_traits<T>::from_int(0);
    std::map<std::string, T> variables;

    /** The named variable, or zero when the model did not create it. */
    T get_variable(const std::string& name) const {
        typename std::map<std::string, T>::const_iterator it = variables.find(name);
        return it == variables.end() ? num_traits<T>::from_int(0) : it->second;
    }

    /** Aggregate utilization of queue i (1-based), as the QR models write it. */
    T get_utilization(int i) const { return get_variable("U_" + std::to_string(i)); }

    /** Per-phase utilization, with the `e_i_k` spelling as the fallback. */
    T get_utilization(int i, int k) const {
        const T u = get_variable("U_" + std::to_string(i) + "_" + std::to_string(k));
        if (!(u == num_traits<T>::from_int(0))) return u;
        return get_variable("e_" + std::to_string(i) + "_" + std::to_string(k));
    }

    /** Aggregate mean queue length of queue i (1-based). */
    T get_queue_length(int i) const { return get_variable("Q_" + std::to_string(i)); }

    /** Per-phase mean queue length, as the LR model writes it. */
    T get_queue_length(int i, int k) const {
        return get_variable("Q_" + std::to_string(i) + "_" + std::to_string(k));
    }

    /** The MVA-shaped variant's spellings (Mapqn_bnd_lr_mva). */
    T get_utilization_mva(int i, int k) const {
        return get_variable("UN_" + std::to_string(i) + "_" + std::to_string(k));
    }
    T get_queue_length_mva(int i, int k) const {
        return get_variable("QN_" + std::to_string(i) + "_" + std::to_string(k));
    }
};

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_SOLUTION_H
