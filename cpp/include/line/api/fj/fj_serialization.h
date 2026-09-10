/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_SERIALIZATION_H
#define LINE_API_FJ_SERIALIZATION_H

/**
 * Blocking probability and pseudoserver delay of serialization phases.
 *
 * Templated port of matlab/src/api/fj/fj_serialization.m.
 *
 * A serialization phase is a stretch of execution protected by an exclusive
 * lock, so at most one of the M circulating jobs may occupy it. Treating the
 * other M-1 jobs as independently placed in proportion to the residence times,
 *
 *   P_s(M) = 1 - [ 1 - R_s(M)/R(M) ]^(M-1),   R(M) = R_0 + sum_s R_s(M),
 *
 * and the delay charged at the pseudoserver is alpha R_s(M), with alpha = 1/2
 * for an arrival uniform in a lightly utilized phase, the regime in which the
 * approximation is stated.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [P, delay, Rtot] of fj_serialization. */
template <class T>
struct FJSerializationResult {
    std::vector<T> P;
    std::vector<T> delay;
    T Rtot;
};

/**
 * @param Rs    mean residence time inside each serialization phase
 * @param R0    mean residence time in the nonserialized phase
 * @param M     number of circulating jobs, M >= 1
 * @param alpha fraction of the phase charged to a blocked job, in [0,1]
 * @return      the blocking probabilities, the pseudoserver delays and the cycle time
 */
template <class T>
FJSerializationResult<T> fj_serialization(const std::vector<T>& Rs, const T& R0, unsigned M,
                                          const T& alpha = num_traits<T>::from_double(0.5)) {
    const std::size_t S = Rs.size();
    if (S < 1) throw InputError("fj_serialization: at least one serialization phase is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t s = 0; s < S; ++s)
        if (Rs[s] < zero)
            throw InputError("fj_serialization: the phase residence times must be non-negative");
    if (R0 < zero)
        throw InputError("fj_serialization: the nonserialized residence time must be non-negative");
    if (M < 1) throw InputError("fj_serialization: M must be a positive integer");
    if (alpha < zero || alpha > one)
        throw InputError("fj_serialization: alpha must lie in [0,1]");

    T R = R0;
    for (std::size_t s = 0; s < S; ++s) R += Rs[s];
    if (!(R > zero))
        throw NumericError("fj_serialization: the total residence time vanished");

    FJSerializationResult<T> out;
    out.P.resize(S);
    out.delay.resize(S);
    out.Rtot = R;
    for (std::size_t s = 0; s < S; ++s) {
        T pw = one;
        for (unsigned e = 0; e + 1 < M; ++e) pw *= (one - Rs[s] / R);
        out.P[s] = one - pw;
        out.delay[s] = out.P[s] * (alpha * Rs[s]);
        out.Rtot += out.delay[s];
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_SERIALIZATION_H
