/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the given vector and matrix define a valid matrix-
 * exponential representation.
 *
 * @param alpha Initial vector of the matrix-exponential distribution to check
 * @param A Matrix parameter of the matrix-exponential distribution to check
 * @param prec Numerical precision. The default value is 1e-14.
 * @return True, if the matrix is a square matrix, the vector and
 *         the matrix have the same size, the dominant eigenvalue
 *         is negative and real
 *
 * Note: This procedure does not check the positivity of the density!
 * Call 'checkMEPositiveDensity' if it is needed, but keep in mind
 * that it can be time-consuming, while this procedure is fast.
 */
fun checkMERepresentation(alpha: Matrix, A: Matrix, prec: Double = 1e-14): Boolean {
    // Check if matrix is square
    if (A.numRows != A.numCols) {
        return false
    }

    // Check if vector and matrix have same size
    if (alpha.length() != A.numRows) {
        return false
    }

    // Check if sum of alpha is between 0 and 1
    val sumAlpha = alpha.elementSum()
    if (sumAlpha < -prec * alpha.length() || sumAlpha > 1 + prec * alpha.length()) {
        return false
    }

    // Check if all eigenvalues have negative real parts
    val eigenvalues = A.eig()
    for (ev in eigenvalues) {
        if (ev.real >= prec) {
            return false
        }
    }

    // Find dominant eigenvalue (closest to zero)
    var maxEv = eigenvalues[0]
    for (ev in eigenvalues) {
        if (abs(ev.real) < abs(maxEv.real)) {
            maxEv = ev
        }
    }

    // Check if dominant eigenvalue is real
    if (abs(maxEv.imaginary) > prec) {
        return false
    }

    return true
}

/**
 * Overload for DoubleArray alpha.
 */
fun checkMERepresentation(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): Boolean {
    return checkMERepresentation(Matrix(alpha), A, prec)
}
