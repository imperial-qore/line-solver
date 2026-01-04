/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the vector is a valid probability vector: the
 * vector has only non-negative elements, the sum of the
 * vector elements is 1.
 *
 * If parameter "sub" is set to true, it checks if the
 * vector is a valid substochastic vector: the vector has
 * only non-negative elements, the sum of the elements are
 * less than or equal to 1.
 *
 * @param pi The vector to check.
 * @param sub If false, the procedure checks for stochastic, if
 *        true, it checks for sub-stochastic property. The
 *        default value is false.
 * @param prec Numerical precision. Entries with absolute value
 *        less than prec are considered to be zeros. The
 *        default value is 1e-14.
 * @return The result of the check.
 */
fun checkProbVector(pi: Matrix, sub: Boolean = false, prec: Double = 1e-14): Boolean {
    val n = pi.length()

    // Check for negative elements
    for (i in 0 until n) {
        if (pi[i] < -prec) {
            return false
        }
    }

    val sum = pi.elementSum()

    if (sub) {
        // Check if sum is less than or equal to 1
        if (sum > 1 + prec * n) {
            return false
        }
    } else {
        // Check if sum equals 1
        if (abs(sum - 1.0) > prec * n) {
            return false
        }
    }

    return true
}

/**
 * Overload for DoubleArray input.
 */
fun checkProbVector(pi: DoubleArray, sub: Boolean = false, prec: Double = 1e-14): Boolean {
    return checkProbVector(Matrix(pi), sub, prec)
}
