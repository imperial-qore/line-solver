/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_GAMMA_H
#define LINE_API_TRACE_TRACE_GAMMA_H

/**
 * Autocorrelation decay rate of a trace: the gamma of the geometric model
 * rho(k) = rho0 * gamma^k, with rho0 = (1 - 1/scv)/2 fixed by the second
 * moment and gamma fitted by least squares on the empirical acf.
 *
 * Templated port of `jar/src/main/java/jline/api/trace/Trace_var.java#trace_gamma`,
 * cross-checked against matlab/lib/kpctoolbox/trace/trace_gamma.m.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB fits gamma by nonlinear regression
 * (nlinfit with a fair robust weight, falling back to lsqcurvefit), the JAR
 * by an exhaustive search over the ten grid points 0.990, 0.991, ..., 0.999.
 * The JAR can therefore never report a decay rate outside that window and its
 * resolution is 1e-3, so the two agree only when the true rate happens to sit
 * on the grid. The grid search is what is ported, because it is the only one
 * of the two that is deterministic and free of an optimizer dependency; the
 * grid is exposed as a parameter so a caller can refine it.
 *
 * REFERENCE DEFECT (JAR): trace_gamma builds the lags 1..min(limit, n-1) and
 * then indexes the trace_acf result by those lags, but trace_acf silently
 * DROPS every lag above n-2. For any trace with n <= limit+1 -- with the
 * default limit of 1000, any trace of at most 1001 samples -- the returned
 * array is one element shorter than the loop bound and the method throws
 * ArrayIndexOutOfBoundsException. This port builds the lags as
 * 1..min(limit, n-2) so that the residual is computed on the acf that
 * actually exists.
 *
 * ARITHMETIC: rho0 and the residual sum of squares are rational in the
 * samples, and the grid points are raised to INTEGER lag powers, so the whole
 * fit is a field computation; instantiated for double, Rational and Real50.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_acf.h"
#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_scv.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of trace_gamma, mirroring [GAMMA, RHO0, RESIDUALS]. */
template <class T>
struct TraceGammaResult {
    T gamma;      ///< best decay rate on the grid
    T rho0;       ///< (1 - 1/scv)/2, the GE-type lag-0 amplitude
    T residuals;  ///< sum of squared deviations from the fitted model
};

/**
 * @param S     the trace
 * @param limit largest lag considered (the JAR default is 1000)
 * @param grid  candidate decay rates; the JAR grid 0.990..0.999 by default
 */
template <class T>
TraceGammaResult<T> trace_gamma(const std::vector<T>& S, long limit = 1000,
                                const std::vector<T>& grid = std::vector<T>()) {
    detail::require_nonempty(S, "trace_gamma");
    const long n = static_cast<long>(S.size());
    if (n < 4) throw InputError("trace_gamma: the fit needs at least four samples");

    const long maxlag = std::min<long>(limit, n - 2);
    std::vector<int> lags;
    for (long l = 1; l <= maxlag; ++l) lags.push_back(static_cast<int>(l));
    const std::vector<T> rho = trace_acf(S, lags);

    // rho0 uses the POPULATION scv, as both references do (they form
    // M2 - M1^2 explicitly rather than calling var).
    const T scv = trace_scv(S, false);
    if (scv == num_traits<T>::from_int(0))
        throw NumericError("trace_gamma: the trace has zero scv, rho0 is undefined");
    TraceGammaResult<T> out;
    out.rho0 = num_traits<T>::from_rational(1, 2) *
               (num_traits<T>::from_int(1) - num_traits<T>::from_int(1) / scv);

    std::vector<T> g = grid;
    if (g.empty())
        for (long k = 990; k <= 999; ++k) g.push_back(num_traits<T>::from_rational(k, 1000));

    bool first = true;
    for (std::size_t a = 0; a < g.size(); ++a) {
        T res = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < rho.size(); ++i) {
            const T expected = out.rho0 * num_pow_int(g[a], static_cast<unsigned>(lags[i]));
            const T d = rho[i] - expected;
            res += d * d;
        }
        if (first || res < out.residuals) {
            out.residuals = res;
            out.gamma = g[a];
            first = false;
        }
    }
    if (first) throw InputError("trace_gamma: the candidate grid is empty");
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_GAMMA_H
