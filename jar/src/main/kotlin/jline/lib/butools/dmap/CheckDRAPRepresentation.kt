/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs

/**
 * Checks if the input matrices define a discrete time RAP.
 *
 * Matrices H0 and H1 must have the same size, the dominant eigenvalue
 * of H0 is real and less than 1, and the rowsum of H0+H1 is 1 (up to
 * the numerical precision).
 *
 * @param D0 The D0 matrix of the DRAP to check
 * @param D1 The D1 matrix of the DRAP to check
 * @param prec Numerical precision, default value is 1e-14
 * @return True if the matrices define a valid DRAP
 */
fun checkDRAPRepresentation(D0: Matrix, D1: Matrix, prec: Double = 1e-14): Boolean {
    // Check if D0 and D1 have the same size and D0 is square
    if (D0.numRows != D1.numRows || D0.numCols != D1.numCols || D0.numRows != D0.numCols) {
        return false
    }

    // Check if rowsums of D0+D1 equal 1
    val D0D1 = D0.add(D1)
    for (i in 0 until D0D1.numRows) {
        var rowSum = 0.0
        for (j in 0 until D0D1.numCols) {
            rowSum += D0D1[i, j]
        }
        if (abs(rowSum - 1.0) > prec) {
            return false
        }
    }

    // Get eigenvalues and sort by absolute value (descending)
    val eigenvalues: List<Complex> = D0.eig()
    val sortedEigs = eigenvalues.sortedByDescending { ev: Complex -> abs(ev.real * ev.real + ev.imaginary * ev.imaginary) }

    if (sortedEigs.isEmpty()) {
        return false
    }

    val maxEv = sortedEigs[0]

    // Check if dominant eigenvalue is real
    if (abs(maxEv.imaginary) > prec) {
        return false
    }

    // Check if dominant eigenvalue is not greater than 1
    if (maxEv.real > 1 + prec) {
        return false
    }

    return true
}
