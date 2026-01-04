/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.lib.butools.mc.checkProbMatrix
import jline.lib.butools.mc.checkProbVector
import jline.util.matrix.Matrix

/**
 * Checks if the given vector and matrix define a valid discrete phase-type representation.
 *
 * @param alpha Initial vector of the phase-type distribution to check
 * @param A Transient generator of the phase-type distribution to check
 * @param prec Numerical precision. The default value is 1e-14.
 * @return True if vector alpha is a probability vector and matrix A is substochastic, and they have the same size.
 */
fun checkDPHRepresentation(alpha: Matrix, A: Matrix, prec: Double = 1e-14): Boolean {
    // Check if alpha and A have the same size
    if (alpha.length() != A.numRows) {
        return false
    }

    // Check if A is a valid transient (substochastic) probability matrix
    if (!checkProbMatrix(A, true, prec)) {
        return false
    }

    // Check if alpha is a valid substochastic probability vector
    if (!checkProbVector(alpha, true, prec)) {
        return false
    }

    return true
}

/**
 * Overload for DoubleArray alpha.
 */
fun checkDPHRepresentation(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): Boolean {
    return checkDPHRepresentation(Matrix(alpha), A, prec)
}
