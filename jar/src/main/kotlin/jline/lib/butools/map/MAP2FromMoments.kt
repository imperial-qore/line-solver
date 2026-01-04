/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.ph.ph2From3Moments
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Returns a MAP(2) which has the same 3 marginal moments
 * and lag-1 autocorrelation as given.
 *
 * @param moms First three marginal moments of the inter-arrival times
 * @param corr1 The lag-1 autocorrelation of the inter-arrival times
 * @return Pair of (D0, D1) matrices of the MAP(2)
 */
fun map2FromMoments(moms: DoubleArray, corr1: Double): Pair<Matrix, Matrix> {
    val m1 = moms[0]
    val m2 = moms[1]

    val prec = 1e-12

    // If we have an exponential distribution, we do not allow correlation
    if (abs(m2 - 2.0 * m1 * m1) < prec && abs(corr1) > prec) {
        throw IllegalArgumentException("We do not allow correlation in case of exponentially distributed marginal")
    }

    // Perform PH fitting
    val phResult = ph2From3Moments(moms)
    val tau = phResult.alpha
    val A = phResult.A

    val l1 = -A[0, 0]
    val l2 = -A[1, 1]
    val p = tau[0, 0]
    val alpha = l1 / l2

    // Check the feasibility of the correlation parameter
    val (corrl, corru) = map2CorrelationBounds(moms)
    if (corr1 < corrl) {
        throw IllegalArgumentException("The correlation parameter is too small!")
    }
    if (corr1 > corru) {
        throw IllegalArgumentException("The correlation parameter is too large!")
    }

    val gamma = corr1 * (m2 - m1 * m1) / (m2 / 2.0 - m1 * m1)

    // Perform MAP fitting
    val D0: Matrix
    val D1: Matrix

    if (gamma > 0) {
        val discriminant = (1.0 + alpha * gamma - p * (1.0 - gamma)) * (1.0 + alpha * gamma - p * (1.0 - gamma)) - 4.0 * alpha * gamma
        val a = (1.0 + alpha * gamma - p * (1.0 - gamma) - sqrt(discriminant)) / (2.0 * alpha)
        val b = (1.0 + alpha * gamma - p * (1.0 - gamma) + sqrt(discriminant)) / 2.0

        D0 = Matrix(2, 2)
        D0[0, 0] = -l1
        D0[0, 1] = (1.0 - a) * l1
        D0[1, 0] = 0.0
        D0[1, 1] = -l2

        D1 = Matrix(2, 2)
        D1[0, 0] = a * l1
        D1[0, 1] = 0.0
        D1[1, 0] = (1.0 - b) * l2
        D1[1, 1] = b * l2
    } else if (gamma < 0) {
        val a = gamma / (alpha * gamma - p * (1.0 - gamma))
        val b = p * (1.0 - gamma) - alpha * gamma

        D0 = Matrix(2, 2)
        D0[0, 0] = -l1
        D0[0, 1] = (1.0 - a) * l1
        D0[1, 0] = 0.0
        D0[1, 1] = -l2

        D1 = Matrix(2, 2)
        D1[0, 0] = 0.0
        D1[0, 1] = a * l1
        D1[1, 0] = b * l2
        D1[1, 1] = (1.0 - b) * l2
    } else {
        D0 = Matrix(2, 2)
        D0[0, 0] = -l1
        D0[0, 1] = l1
        D0[1, 0] = 0.0
        D0[1, 1] = -l2

        D1 = Matrix(2, 2)
        D1[0, 0] = 0.0
        D1[0, 1] = 0.0
        D1[1, 0] = p * l2
        D1[1, 1] = (1.0 - p) * l2
    }

    return Pair(D0, D1)
}
