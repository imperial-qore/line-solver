/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_MVA_INTERVAL_H
#define LINE_API_PFQN_PFQN_MVA_INTERVAL_H

/**
 * Exact interval-valued MVA for single-class closed product-form networks.
 *
 * Port of `matlab/src/api/pfqn/pfqn_mva_interval.m`, the algorithm of J. Luthi
 * and G. Haring, "Mean value analysis for queueing network models with
 * intervals as input parameters", Performance Evaluation 32(3):185-215, 1998.
 *
 * WHY CORNERS AND NOT INTERVAL ARITHMETIC. Single-class MVA is monotone in
 * every input -- the throughput decreases in each demand and in the think time
 * and increases in the population, a station's own queue length and residence
 * time increase in its own demand and in the population and decrease in the
 * other demands and in the think time, and the totals increase in every demand
 * and in the population and decrease in the think time (their Theorems 2-5,
 * Table 1). By their Theorem 1 the exact range of a function monotone in each
 * argument is attained AT A CORNER of the input box, so every bound below is
 * one ordinary MVA call at the corner the sign pattern selects: 2*(m+2) calls,
 * m being the number of thick demand intervals. Running the recursion in
 * interval arithmetic instead is valid but far wider, because every input
 * recurs at each step -- the dependency problem, 14x too wide on the paper's
 * own example.
 *
 * WHAT THE INTERVAL IS, AND WHAT IT IS NOT. It is the exact hull of MVA over
 * the input box, conditional on the demands lying in that box; it says nothing
 * about the accuracy of MVA itself. It must therefore NOT be composed with the
 * brackets of SolverBA, which bracket the exact solution of a model whose
 * demands are known. The two answer different questions and intersecting them
 * would claim a guarantee neither provides.
 *
 * Delay stations are folded into Z exactly as `pfqn_mva` folds them: a delay
 * demand interval enters as a term of the think-time interval, and the hull of
 * the sum is the sum of the hulls when the delays vary independently.
 * Load-independent single-server queueing stations only, one class only; the
 * monotonicity theorems cover no other case.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_mva.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Every output of `pfqn_mva_interval`, each a [lower, upper] pair. */
template <class T>
struct MvaIntervalResult {
    T Xlo, Xup;                ///< throughput interval
    Matrix<T> Q;               ///< (M x 2) mean queue length per station
    Matrix<T> U;               ///< (M x 2) utilization enclosure per station
    Matrix<T> R;               ///< (M x 2) residence time per station
    T Rtot_lo, Rtot_up;        ///< total response time
    T Qtot_lo, Qtot_up;        ///< total number of jobs at the stations
};

/**
 * @param L   (M x 2) demand intervals, column 0 lower and column 1 upper.
 * @param nlo lower population endpoint (integer, at least 1).
 * @param nup upper population endpoint.
 * @param zlo lower think-time endpoint.
 * @param zup upper think-time endpoint.
 */
template <class T>
MvaIntervalResult<T> pfqn_mva_interval(const Matrix<T>& L, int nlo, int nup, const T& zlo,
                                       const T& zup) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (L.empty()) throw InputError("pfqn_mva_interval requires at least one queueing station");
    if (L.cols() != 2)
        throw InputError(
            "pfqn_mva_interval is a single-class method: L must be M x 2, one [lower upper] demand "
            "interval per station");
    const std::size_t M = L.rows();
    for (std::size_t i = 0; i < M; ++i) {
        if (L(i, 0) < zero || L(i, 1) < zero)
            throw InputError("pfqn_mva_interval: demands must be nonnegative");
        if (L(i, 0) > L(i, 1))
            throw InputError(
                "pfqn_mva_interval: an interval lower endpoint exceeds its upper endpoint");
    }
    if (zlo < zero || zup < zero)
        throw InputError("pfqn_mva_interval: think times must be nonnegative");
    if (zlo > zup)
        throw InputError("pfqn_mva_interval: the think-time interval is inverted");
    if (nlo > nup) throw InputError("pfqn_mva_interval: the population interval is inverted");
    if (nlo < 1)
        throw InputError(
            "pfqn_mva_interval requires a population interval with at least one job; the "
            "monotonicity theorems assume n >= 1");

    const T tiny = num_traits<T>::from_double(lang::GlobalConstants::Zero);
    std::vector<bool> thick(M, false);
    Matrix<T> Llo(M, 1, zero), Lup(M, 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        Llo(i, 0) = L(i, 0);
        Lup(i, 0) = L(i, 1);
        thick[i] = L(i, 1) > T(L(i, 0) + tiny);
    }

    auto call = [&](const Matrix<T>& d, int n, const T& z) {
        return pfqn_mva(d, std::vector<int>(1, n), Matrix<T>(1, 1, z));
    };

    MvaIntervalResult<T> out;
    out.Q = Matrix<T>(M, 2, zero);
    out.U = Matrix<T>(M, 2, zero);
    out.R = Matrix<T>(M, 2, zero);

    // S1: the throughput upper bound, and the upper bounds of the stations whose
    // demand is thin -- their own demand is fixed, so lowering the others
    // maximizes them.
    const MvaResult<T> s1 = call(Llo, nup, zlo);
    out.Xup = s1.XN[0];
    // S2: the same quantities at the opposite corner, giving the lower bounds.
    const MvaResult<T> s2 = call(Lup, nlo, zup);
    out.Xlo = s2.XN[0];
    for (std::size_t i = 0; i < M; ++i) {
        if (thick[i]) continue;
        out.Q(i, 1) = s1.QN(i, 0);
        out.R(i, 1) = s1.CN(i, 0);
        out.Q(i, 0) = s2.QN(i, 0);
        out.R(i, 0) = s2.CN(i, 0);
    }

    // S3/S4: the totals increase in the demands and the population and decrease
    // in the think time, so their corners differ from the throughput's.
    const MvaResult<T> s3 = call(Llo, nlo, zup);
    const MvaResult<T> s4 = call(Lup, nup, zlo);
    out.Rtot_lo = zero;
    out.Qtot_lo = zero;
    out.Rtot_up = zero;
    out.Qtot_up = zero;
    for (std::size_t i = 0; i < M; ++i) {
        out.Rtot_lo += s3.CN(i, 0);
        out.Qtot_lo += s3.QN(i, 0);
        out.Rtot_up += s4.CN(i, 0);
        out.Qtot_up += s4.QN(i, 0);
    }

    // S5/S6: one pair of calls per thick station, its own demand at the endpoint
    // that maximizes (minimizes) it and the others at the opposite endpoint.
    for (std::size_t k = 0; k < M; ++k) {
        if (!thick[k]) continue;
        Matrix<T> d = Llo;
        d(k, 0) = Lup(k, 0);
        const MvaResult<T> s5 = call(d, nup, zlo);
        out.Q(k, 1) = s5.QN(k, 0);
        out.R(k, 1) = s5.CN(k, 0);
        d = Lup;
        d(k, 0) = Llo(k, 0);
        const MvaResult<T> s6 = call(d, nlo, zup);
        out.Q(k, 0) = s6.QN(k, 0);
        out.R(k, 0) = s6.CN(k, 0);
    }

    // U = X D is not covered by the monotonicity table, so it is enclosed by the
    // product of the two intervals, capped at the range of a single-server
    // utilization. Where the demand is thin the product is already exact.
    for (std::size_t i = 0; i < M; ++i) {
        out.U(i, 0) = T(out.Xlo * Llo(i, 0));
        const T up = T(out.Xup * Lup(i, 0));
        out.U(i, 1) = up < one ? up : one;
    }
    return out;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_MVA_INTERVAL_H
