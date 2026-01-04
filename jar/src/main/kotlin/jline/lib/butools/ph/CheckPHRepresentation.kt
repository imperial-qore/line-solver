/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.lib.butools.mc.checkGenerator
import jline.lib.butools.mc.checkProbVector
import jline.util.matrix.Matrix

/**
 * Checks if the given vector and matrix define a valid phase-
 * type representation.
 *
 * @param alpha Initial vector of the phase-type distribution to check
 * @param A Transient generator of the phase-type distribution to check
 * @param prec Numerical precision. The default value is 1e-14.
 * @return True, if vector alpha is a probability vector and matrix
 *         A is a transient generator, and they have the same size.
 */
fun checkPHRepresentation(alpha: Matrix, A: Matrix, prec: Double = 1e-14): Boolean {
    // Check if vector and matrix have same size
    if (alpha.length() != A.numRows) {
        return false
    }

    // Check if A is a valid transient generator
    if (!checkGenerator(A, transient = true, prec = prec)) {
        return false
    }

    // Check if alpha is a valid substochastic probability vector
    if (!checkProbVector(alpha, sub = true, prec = prec)) {
        return false
    }

    return true
}

/**
 * Overload for DoubleArray alpha.
 */
fun checkPHRepresentation(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): Boolean {
    return checkPHRepresentation(Matrix(alpha), A, prec)
}
