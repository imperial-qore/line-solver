/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the matrix is a valid generator matrix: the
 * matrix is a square matrix, the matrix has positive or
 * zero off-diagonal elements, the diagonal of the matrix
 * is negative, the rowsum of the matrix is 0.
 *
 * If the "transient" parameter is set to true, it checks
 * if the real part of the maximum absolute eigenvalue is
 * less than zero and the rowsum is equal or less than 0.
 *
 * @param Q The generator to check.
 * @param transient If true, the procedure checks if Q is a transient
 *        generator, otherwise it checks if it is a valid
 *        generator. The default value is false.
 * @param prec Entries with absolute value less than prec are
 *        considered to be zeros. The default value is 1e-14.
 * @return The result of the check.
 */
fun checkGenerator(Q: Matrix, transient: Boolean = false, prec: Double = 1e-14): Boolean {
    // Check if matrix is square
    if (Q.numRows != Q.numCols) {
        return false
    }

    val N = Q.numRows

    // Check if diagonal is negative
    for (i in 0 until N) {
        if (Q[i, i] >= prec) {
            return false
        }
    }

    // Check if off-diagonal elements are non-negative
    for (i in 0 until N) {
        for (j in 0 until N) {
            if (i != j && Q[i, j] < -prec) {
                return false
            }
        }
    }

    if (transient) {
        // Check if rowsum is less than or equal to 0
        for (i in 0 until N) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += Q[i, j]
            }
            if (rowSum > prec) {
                return false
            }
        }

        // Check if maximum eigenvalue has negative real part
        val eigenvalues = Q.eig()
        for (ev in eigenvalues) {
            if (ev.real >= prec) {
                return false
            }
        }
    } else {
        // Check if rowsum is 0
        for (i in 0 until N) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += Q[i, j]
            }
            if (abs(rowSum) > prec) {
                return false
            }
        }
    }

    return true
}
