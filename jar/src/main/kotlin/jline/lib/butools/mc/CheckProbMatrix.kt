/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the matrix is a valid probability matrix: the
 * matrix is a square matrix, the matrix has positive or
 * zero elements, the rowsum of the matrix is 1.
 *
 * If "transient" is true, it checks if the matrix is a
 * valid transient probability matrix: the matrix is a square
 * matrix, the matrix has positive or zero elements,
 * the rowsum of the matrix is less than or equal to 1,
 * the maximum absolute eigenvalue is less than 1.
 *
 * @param P The matrix to check.
 * @param transient If true, the procedure checks if P is a transient
 *        probability matrix, otherwise it checks if it is
 *        a valid probability matrix. The default value is false.
 * @param prec Entries with absolute value less than prec are
 *        considered to be zeros. The default value is 1e-14.
 * @return The result of the check.
 */
fun checkProbMatrix(P: Matrix, transient: Boolean = false, prec: Double = 1e-14): Boolean {
    // Check if matrix is square
    if (P.numRows != P.numCols) {
        return false
    }

    val N = P.numRows

    // Check for negative elements
    for (i in 0 until N) {
        for (j in 0 until N) {
            if (P[i, j] < -prec) {
                return false
            }
        }
    }

    if (transient) {
        // Check if rowsum is less than or equal to 1
        for (i in 0 until N) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += P[i, j]
            }
            if (rowSum - 1.0 > N * prec) {
                return false
            }
        }

        // Check if maximum eigenvalue is less than 1
        val eigenvalues = P.eig()
        for (ev in eigenvalues) {
            if (ev.real >= 1.0 - prec) {
                return false
            }
        }
    } else {
        // Check if rowsum is 1
        for (i in 0 until N) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += P[i, j]
            }
            if (abs(rowSum - 1.0) > N * prec) {
                return false
            }
        }
    }

    return true
}
