/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_ENV_MAP_H
#define LINE_API_SNC_ENV_MAP_H

/**
 * MGF arrival envelope of a MAP/MMPP flow with unit-size jobs.
 *
 * For a Markovian arrival process (D0,D1) counting N(0,t) unit-work jobs,
 * `E[exp(theta*N(0,t))] = pi*expm((D0+D1*exp(theta))*t)*1`. With lstar the
 * eigenvalue of maximal real part of `A(theta)=D0+D1*e^theta` and v > 0 its
 * right Perron eigenvector, bounding `1 <= v/min(v)` entrywise gives
 *
 *   rho(theta)   = lstar/theta,
 *   sigma(theta) = log(max(v)/min(v))/theta,
 *
 * the standard exponential-form envelope of a Markov-modulated source. The burst
 * term is what the modulating chain contributes: 0 for a one-phase MAP, where
 * this reproduces `snc_env_poisson` exactly, and positive for an MMPP.
 *
 * THE PERRON PAIR IS COMPUTED BY POWER ITERATION ON THE SHIFTED MATRIX `A+cI`,
 * not by a general eigensolver. A(theta) is essentially nonnegative, so the
 * shift makes it nonnegative with a positive diagonal, hence primitive whenever
 * the MAP is irreducible, and the iteration converges to the pair the bound
 * needs without asking a general solver which of its eigenvectors is the
 * positive one. The MATLAB reference uses `eig` and agrees to machine
 * precision; the JAR port iterates the same way.
 *
 * Port of matlab/src/api/snc/snc_env_map.m. Reference: C.-S. Chang, Performance
 * Guarantees in Communication Networks, Springer 2000, Ch. 7.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace snc {

/**
 * @param D0    hidden-transition generator block of the MAP
 * @param D1    arrival-transition block of the MAP
 * @param theta Chernoff parameter, theta > 0
 */
inline Env snc_env_map(const Matrix<double>& D0, const Matrix<double>& D1, double theta) {
    if (theta <= 0) throw UnsupportedError("snc_env_map: theta must be positive");
    const std::size_t n = D0.rows();
    if (D0.cols() != n || D1.rows() != n || D1.cols() != n)
        throw UnsupportedError("snc_env_map: D0 and D1 must be square and of equal size");

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    const double etheta = std::exp(theta);
    double shift = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) A[i][j] = D0(i, j) + D1(i, j) * etheta;
        shift = std::max(shift, -A[i][i]);
    }
    shift += 1.0;
    for (std::size_t i = 0; i < n; ++i) A[i][i] += shift;

    std::vector<double> v(n, 1.0), w(n, 0.0);
    double lambdaShift = 0.0;
    for (int it = 0; it < 100000; ++it) {
        for (std::size_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (std::size_t j = 0; j < n; ++j) acc += A[i][j] * v[j];
            w[i] = acc;
        }
        double norm = 0.0;
        for (std::size_t i = 0; i < n; ++i) norm = std::max(norm, std::fabs(w[i]));
        if (!(norm > 0))
            throw UnsupportedError(
                "snc_env_map: MAP is not irreducible: the Perron eigenvector is not positive");
        double delta = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            w[i] /= norm;
            delta = std::max(delta, std::fabs(w[i] - v[i]));
        }
        v = w;
        lambdaShift = norm;
        if (delta < 1e-14) break;
    }
    double vmin = std::numeric_limits<double>::infinity(), vmax = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        vmin = std::min(vmin, v[i]);
        vmax = std::max(vmax, v[i]);
    }
    if (!(vmin > 0))
        throw UnsupportedError(
            "snc_env_map: MAP is not irreducible: the Perron eigenvector is not positive");
    return Env{std::log(vmax / vmin) / theta, (lambdaShift - shift) / theta};
}

/** The same envelope as a function of theta. */
inline Envelope snc_env_map_fn(const Matrix<double>& D0, const Matrix<double>& D1) {
    return [D0, D1](double theta) { return snc_env_map(D0, D1, theta); };
}

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_ENV_MAP_H
