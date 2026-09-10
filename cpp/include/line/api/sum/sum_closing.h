/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SUM_SUM_CLOSING_H
#define LINE_API_SUM_SUM_CLOSING_H

/**
 * Closing method for open and mixed non-product-form queueing networks,
 * solved with the summation method.
 *
 * Templated port of matlab/src/api/sum/sum_closing.m, cross-checked against
 * jar/src/main/java/jline/api/sum/Sum_closing.java.
 *
 * The external world of the open classes is replaced by one extra -/G/1
 * station with demand 1/(Ropen lambda0_r) for open class r, service SCV equal
 * to that class's interarrival SCV, and unit visit ratio. The resulting closed
 * network is solved by sum_closed with a large closing population Kclosed for
 * the open classes (5000 by default, the value recommended for the summation
 * method). Closed classes pass through untouched, which is what makes the
 * method applicable to mixed networks. The open-class throughput approaches
 * lambda0 from below as Kclosed grows.
 *
 * Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing Networks
 * and Markov Chains, 2nd ed., Wiley, 2006, Sec. 10.1.5.
 *
 * ARITHMETIC: it is a wrapper around sum_closed, whose bisection stops on a
 * tolerance, so it carries the same has_transcendental gate. The closing
 * itself is exact: it only builds one extra row of demands.
 *
 * MATLAB marks the open classes by N(r) = Inf; here the class is open exactly
 * when lambda0(r) > 0, and the N entry of an open class is ignored (it is
 * overwritten by Kclosed), so no infinite population ever has to be
 * represented in the number type.
 */

#include <cstddef>
#include <vector>

#include "line/api/sum/sum_closed.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace sum {

/** Mirrors the [XN, QN, UN, RN, TN, it] return list of the MATLAB function. */
template <class T>
struct SumClosingResult {
    std::vector<T> XN;  ///< (R) class throughputs
    Matrix<T> QN;       ///< (M x R) queue lengths at the original stations
    Matrix<T> UN;       ///< (M x R) utilizations at the original stations
    Matrix<T> RN;       ///< (M x R) residence times at the original stations
    std::vector<T> TN;  ///< (R) mean response time in the original network
    std::size_t it = 0;
};

/** Closing controls; Kclosed is the population given to the open classes. */
struct ClosingOptions {
    long Kclosed = 5000;
    SumOptions sum;
};

/**
 * @param lambda0 (R) external arrival rates, zero for a closed class
 * @param scva    (R) interarrival-time SCVs of the open classes, 1 if Poisson
 * @param L       (M x R) service demands, open-class visits per external arrival
 * @param mi      (M) servers per station
 * @param scv     (M x R) service-time SCVs
 * @param N       (R) populations of the closed classes; entries of open classes are ignored
 * @param Z       (R) think times
 * @param options tolerances, iteration caps and the closing method
 */
template <class T>
SumClosingResult<T> sum_closing(const std::vector<T>& lambda0, const std::vector<T>& scva,
                                const Matrix<T>& L, const std::vector<Servers>& mi,
                                const Matrix<T>& scv, const std::vector<long>& N,
                                const std::vector<T>& Z,
                                const ClosingOptions& options = ClosingOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "sum_closing requires transcendental arithmetic: it drives sum_closed, whose "
                  "bisection stops on a tolerance");
    const std::size_t M = L.rows(), R = L.cols();
    if (lambda0.size() != R || scva.size() != R || N.size() != R || Z.size() != R)
        throw InputError("sum_closing: an input disagrees with L on the class count");
    if (mi.size() != M) throw InputError("sum_closing: L and mi disagree on the station count");
    if (scv.rows() != M || scv.cols() != R) throw InputError("sum_closing: scv has the wrong shape");
    if (options.Kclosed <= 0) throw InputError("sum_closing: the closing population must be positive");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> open;
    for (std::size_t r = 0; r < R; ++r)
        if (lambda0[r] > zero) open.push_back(r);
    if (open.empty())
        throw InputError("sum_closing: no open class, use sum_closed for closed networks");
    const long Ropen = static_cast<long>(open.size());

    // augment with the closing -/G/1 station, visited by the open classes only
    Matrix<T> Laug(M + 1, R, zero);
    Matrix<T> scvaug(M + 1, R, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            Laug(i, r) = L(i, r);
            scvaug(i, r) = scv(i, r);
        }
    std::vector<long> Naug = N;
    for (std::size_t r : open) {
        Laug(M, r) = one / (num_traits<T>::from_int(Ropen) * lambda0[r]);
        scvaug(M, r) = scva[r];
        Naug[r] = options.Kclosed;
    }
    std::vector<Servers> miaug = mi;
    miaug.push_back(Servers::of(1));

    const SumClosedResult<T> in = sum_closed(Laug, Naug, Z, miaug, scvaug, options.sum);

    SumClosingResult<T> out;
    out.XN = in.XN;
    out.it = in.it;
    out.QN = Matrix<T>(M, R, zero);
    out.UN = Matrix<T>(M, R, zero);
    out.RN = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            out.QN(i, r) = in.QN(i, r);
            out.UN(i, r) = in.UN(i, r);
            out.RN(i, r) = in.RN(i, r);
        }
    out.TN.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        if (out.XN[r] > zero) {
            T q = zero;
            for (std::size_t i = 0; i < M; ++i) q += out.QN(i, r);
            out.TN[r] = q / out.XN[r];
        }
    }
    return out;
}

}  // namespace sum
}  // namespace line

#endif  // LINE_API_SUM_SUM_CLOSING_H
