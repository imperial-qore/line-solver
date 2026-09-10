/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_HAS_CLASSDEP_ROUTING_H
#define LINE_API_SN_SN_HAS_CLASSDEP_ROUTING_H

/**
 * Port of matlab/src/api/sn/sn_has_classdep_routing.m.
 *
 * True when the routing depends on the class: either some edge switches class
 * (an off-diagonal (r,s) block entry) or two classes leave the same pair (i,j)
 * with different probabilities. SolverMVA reads it to decide whether a
 * per-chain aggregate routing is faithful, so a false negative aggregates a
 * class-dependent model and returns a chain answer for a class question.
 *
 * INDEX SPACE. The reference walks i and j over STATIONS while indexing `sn.rt`,
 * which is over STATEFUL NODES. The two coincide unless the model has a
 * stateful non-station node (a Router, a Cache), and this port reproduces the
 * reference walk rather than correcting it: the predicate gates a dispatch
 * branch, and answering it on a different index space would send models down a
 * different path here than in every other codebase.
 *
 * ARITHMETIC: field. Comparisons only.
 */

#include <cmath>
#include <cstddef>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

template <class T>
bool sn_has_classdep_routing(const qn::NetworkStruct<T>& sn) {
    const std::size_t K = sn.nclasses, M = sn.nstations;
    if (K <= 1) return false;
    const double tol = qn::GlobalConstants::FineTol;
    const std::size_t dim = sn.rt.rows();
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            bool have_shared = false;
            double shared = 0.0;
            for (std::size_t r = 0; r < K; ++r) {
                for (std::size_t s = 0; s < K; ++s) {
                    if (r == s) continue;
                    const std::size_t a = i * K + r, b = j * K + s;
                    if (a >= dim || b >= dim) continue;
                    if (num_traits<T>::to_double(sn.rt(a, b)) > tol) return true;
                }
                const std::size_t a = i * K + r, b = j * K + r;
                if (a >= dim || b >= dim) continue;
                const double p = num_traits<T>::to_double(sn.rt(a, b));
                if (!have_shared) {
                    shared = p;
                    have_shared = true;
                } else if (std::fabs(p - shared) > tol) {
                    return true;
                }
            }
        }
    return false;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_HAS_CLASSDEP_ROUTING_H
