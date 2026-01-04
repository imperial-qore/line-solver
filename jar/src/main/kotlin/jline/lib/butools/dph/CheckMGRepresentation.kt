/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs

/**
 * Checks if the given vector and matrix define a valid matrix-geometric representation.
 *
 * @param alpha Initial vector of the matrix-geometric distribution to check
 * @param A Matrix parameter of the matrix-geometric distribution to check
 * @param prec Numerical precision. The default value is 1e-14.
 * @return True if the matrix is a square matrix, the vector and the matrix have the same size,
 *         the dominant eigenvalue is positive, less than 1 and real.
 *
 * Note: This procedure does not check the positivity of the density!
 * The discrete counterpart of 'CheckMEPositiveDensity' does not exist yet.
 */
fun checkMGRepresentation(alpha: Matrix, A: Matrix, prec: Double = 1e-14): Boolean {
    // Check if matrix is square
    if (A.numRows != A.numCols) {
        return false
    }

    // Check if vector and matrix have matching sizes
    if (alpha.length() != A.numRows) {
        return false
    }

    // Check if sum of alpha is between 0 and 1
    val alphaSum = alpha.elementSum()
    if (alphaSum < -prec || alphaSum > 1 + prec) {
        return false
    }

    // Get eigenvalues and sort by absolute value (descending)
    val eigenvalues: List<Complex> = A.eig()
    val sortedEigs = eigenvalues.sortedByDescending { ev: Complex -> abs(ev.real * ev.real + ev.imaginary * ev.imaginary) }

    if (sortedEigs.isEmpty()) {
        return false
    }

    val maxEv = sortedEigs[0]

    // Check if largest eigenvalue is real
    if (abs(maxEv.imaginary) > prec) {
        return false
    }

    // Check if largest eigenvalue is not greater than 1
    if (maxEv.real > 1 + prec) {
        return false
    }

    return true
}

/**
 * Overload for DoubleArray alpha.
 */
fun checkMGRepresentation(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): Boolean {
    return checkMGRepresentation(Matrix(alpha), A, prec)
}
