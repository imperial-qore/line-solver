/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.mc.checkGenerator
import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the input matrices define a continuous time MAP.
 *
 * Matrices D0 and D1 must have the same size, D0 must be a
 * transient generator matrix, D1 has only non-negative
 * elements, and the rowsum of D0+D1 is 0 (up to the numerical precision).
 *
 * @param D0 The D0 matrix of the MAP to check
 * @param D1 The D1 matrix of the MAP to check
 * @param prec Numerical precision, the default value is 1e-12
 * @return The result of the check
 */
fun checkMAPRepresentation(D0: Matrix, D1: Matrix, prec: Double = 1e-12): Boolean {
    // Check if D0 is a transient generator
    if (!checkGenerator(D0, transient = true, prec = prec)) {
        return false
    }

    // Check if D0 and D1 have same size
    if (D0.numRows != D1.numRows || D0.numCols != D1.numCols) {
        return false
    }

    val n = D0.numRows

    // Check if D1 has non-negative elements
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (D1[i, j] < -prec) {
                return false
            }
        }
    }

    // Check if rowsum of D0+D1 is 0
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            rowSum += D0[i, j] + D1[i, j]
        }
        if (abs(rowSum) > prec) {
            return false
        }
    }

    return true
}

/**
 * Checks if the input matrices define a valid RAP representation.
 *
 * @param H0 The H0 matrix of the RAP to check
 * @param H1 The H1 matrix of the RAP to check
 * @param prec Numerical precision, the default value is 1e-12
 * @return The result of the check
 */
fun checkRAPRepresentation(H0: Matrix, H1: Matrix, prec: Double = 1e-12): Boolean {
    // Check matrix sizes
    if (H0.numRows != H0.numCols || H1.numRows != H1.numCols) {
        return false
    }
    if (H0.numRows != H1.numRows) {
        return false
    }

    val n = H0.numRows

    // Check if rowsum of H0+H1 is 0
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            rowSum += H0[i, j] + H1[i, j]
        }
        if (abs(rowSum) > prec) {
            return false
        }
    }

    // Check if dominant eigenvalue of H0 is negative and real
    val eigenvalues = H0.eig()
    var maxRealEv = Double.NEGATIVE_INFINITY
    for (ev in eigenvalues) {
        if (ev.real > maxRealEv) {
            maxRealEv = ev.real
        }
    }

    if (maxRealEv >= prec) {
        return false
    }

    return true
}
