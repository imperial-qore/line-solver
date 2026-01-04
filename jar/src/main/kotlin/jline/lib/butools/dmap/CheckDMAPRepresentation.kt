/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.mc.checkProbMatrix
import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the input matrices define a discrete time MAP.
 *
 * Matrices D0 and D1 must have the same size, D0 must be a transient
 * probability matrix, D1 has only non-negative elements, and the rowsum
 * of D0+D1 is 1 (up to the numerical precision).
 *
 * @param D0 The D0 matrix of the DMAP to check
 * @param D1 The D1 matrix of the DMAP to check
 * @param prec Numerical precision, default value is 1e-14
 * @return True if the matrices define a valid DMAP
 */
fun checkDMAPRepresentation(D0: Matrix, D1: Matrix, prec: Double = 1e-14): Boolean {
    // Check if D0 is a transient probability matrix
    if (!checkProbMatrix(D0, transient = true, prec = prec)) {
        return false
    }

    // Check if D0 and D1 have the same size
    if (D0.numRows != D1.numRows || D0.numCols != D1.numCols) {
        return false
    }

    // Check for negative elements
    if (D0.elementMin() < -prec || D1.elementMin() < -prec) {
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

    return true
}
