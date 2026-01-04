/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * M. Telek and A. Heindl, "Moment bounds for acyclic discrete and continuous
 * phase-type distributions of second order," in Proc. of UK Performance
 * Evaluation Workshop, UKPEW, 2002
 */
package jline.lib.butools.ph

import jline.lib.butools.APH2ndMomentLowerBound
import jline.lib.butools.APH3rdMomentLowerBound
import jline.lib.butools.APH3rdMomentUpperBound
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Result class for PH2From3Moments.
 */
data class PH2Representation(val alpha: Matrix, val A: Matrix)

/**
 * Returns a PH(2) which has the same 3 moments as given.
 *
 * @param moms The moments to match (length 3).
 * @param prec Numerical precision, default value is 1e-14.
 * @return The PH2Representation containing alpha (initial vector) and A (generator).
 * @throws IllegalArgumentException if the moments are not feasible with a PH(2).
 */
fun ph2From3Moments(moms: DoubleArray, prec: Double = 1e-14): PH2Representation {
    val m1 = moms[0]
    val m2 = moms[1]
    val m3 = moms[2]

    // Check moment bounds
    val m2l = APH2ndMomentLowerBound(m1, 2)
    val m3l = APH3rdMomentLowerBound(m1, m2, 2)
    val m3u = APH3rdMomentUpperBound(m1, m2, 2)

    if (m2 < m2l) {
        throw IllegalArgumentException("The given second moment is not feasible!")
    }
    if (m3 < m3l) {
        throw IllegalArgumentException("The given third moment is not feasible (too small)!")
    }
    if (m3 > m3u) {
        throw IllegalArgumentException("The given third moment is not feasible (too large)!")
    }

    // Check if we have an exponential distribution
    if (abs(m2 / m1 / m1 - 2.0) < prec) {
        val alpha = Matrix(1, 1)
        alpha[0, 0] = 1.0
        val A = Matrix(1, 1)
        A[0, 0] = -1.0 / m1
        return PH2Representation(alpha, A)
    }

    // Calculate parameters
    val b = 3.0 * m1 * m2 - m3
    val c = 3.0 * m2 * m2 - 2.0 * m1 * m3
    val e = -2.0 * m1 * m1 + m2
    var a = b * b + 6.0 * c * e
    if (a < 0) {
        a = 0.0
    }
    a = sqrt(a)

    val lambda1: Double
    val lambda2: Double
    val p: Double

    when {
        c > 0 -> {
            lambda1 = (b - a) / c
            lambda2 = (b + a) / c
            p = (-b - 6.0 * m1 * e + a) / (b + a)
        }
        c < 0 -> {
            lambda1 = (b + a) / c
            lambda2 = (b - a) / c
            p = (b + 6.0 * m1 * e + a) / (-b + a)
        }
        else -> {
            lambda1 = 0.0
            lambda2 = 1.0 / m1
            p = 0.0
        }
    }

    // Return the result
    val alpha = Matrix(1, 2)
    alpha[0, 0] = p
    alpha[0, 1] = 1.0 - p

    val A = Matrix(2, 2)
    A[0, 0] = -lambda1
    A[0, 1] = lambda1
    A[1, 0] = 0.0
    A[1, 1] = -lambda2

    return PH2Representation(alpha, A)
}
