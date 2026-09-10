/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_RPS_H
#define LINE_API_INFER_INFER_RPS_H

/**
 * Regression for Processor Sharing (RPS) demand estimator.
 *
 * Templated port of matlab/src/api/infer/infer_rps.m. No JAR counterpart.
 *
 * Mean value analysis of a PS station gives E[R_r] = E[D_r] E[Qbar_A]/V, with
 * Qbar_A the total number of jobs seen on admission INCLUDING the arriving job
 * and V the number of servers. The demand of each class is the non-negative
 * least squares fit of its response times against Qbar_A/V.
 *
 * MATLAB calls lsqnonneg, but the design matrix here has a SINGLE column, so
 * the non-negative least squares problem has the closed form
 * max(0, a.b / a.a): the unconstrained minimizer is the ordinary projection
 * and the active-set method returns 0 exactly when it is negative. The port
 * evaluates that closed form, so it needs no optimizer and reproduces
 * lsqnonneg exactly rather than approximately.
 *
 * MATLAB takes the class count from max(class); the port does the same, so a
 * trailing class with no samples at all is simply not represented in the
 * output, exactly as in MATLAB.
 *
 * ARITHMETIC: two inner products, one division and one comparison against
 * zero, so a finite field computation, exact in the exact instantiation. The
 * response times are non-negative and Qbar_A >= 1, so the estimate can only
 * hit the bound when every response time is zero.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * @param rt  (n) response time samples
 * @param cls (n) class of each sample, 0-based
 * @param ql  (n x R) per-class queue lengths at arrival, excluding the arriving job
 * @param V   number of servers of the PS station
 * @return    (max(cls)+1) estimated mean service demands
 */
template <class T>
std::vector<T> infer_rps(const std::vector<T>& rt, const std::vector<std::size_t>& cls,
                         const Matrix<T>& ql, long V) {
    const std::size_t n = rt.size();
    if (cls.size() != n) throw InputError("infer_rps: rt and class disagree on the sample count");
    if (ql.rows() != n) throw InputError("infer_rps: ql and rt disagree on the sample count");
    if (V <= 0) throw InputError("infer_rps: the number of servers must be positive");
    if (n == 0) throw InputError("infer_rps: no samples");

    std::size_t R = 0;
    for (std::size_t c : cls) R = c + 1 > R ? c + 1 : R;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T Vt = num_traits<T>::from_int(V);

    std::vector<T> demand(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        T aa = zero, ab = zero;
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (cls[i] != r) continue;
            ++count;
            T qbar = one;  // the arriving job itself
            for (std::size_t c = 0; c < ql.cols(); ++c) qbar += ql(i, c);
            const T a = qbar / Vt;
            aa += a * a;
            ab += a * rt[i];
        }
        if (count == 0)
            throw InputError("infer_rps: a class below the maximum has no samples at all");
        if (aa == zero) throw NumericError("infer_rps: degenerate regressor for a class");
        const T x = ab / aa;
        demand[r] = x > zero ? x : zero;  // the non-negativity bound of lsqnonneg
    }
    return demand;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_RPS_H
