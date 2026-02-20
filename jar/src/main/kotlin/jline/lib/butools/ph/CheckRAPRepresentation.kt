/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Checks if the given matrices define a valid Rational Arrival Process (RAP) representation.
 *
 * A RAP is valid if:
 * - H0 and H1 are square matrices of the same size
 * - H0 + H1 forms a valid infinitesimal generator (row sums = 0)
 * - All eigenvalues of H0 have negative real parts
 * - Dominant eigenvalue of H0 is negative and real
 *
 * @param H0 The H0 matrix of the RAP (hidden transition rates)
 * @param H1 The H1 matrix of the RAP (visible transition rates)
 * @param prec Numerical precision. The default value is 1e-14.
 * @return True if the representation is valid, false otherwise
 *
 * Note: RAP is a generalization of MAP. Unlike MAP, RAP allows H1 to have
 * negative entries, but H0 + H1 must still form a valid generator.
 */
fun checkRAPRepresentation(H0: Matrix, H1: Matrix, prec: Double = 1e-14): Boolean {
    // Check if H0 is square
    if (H0.numRows != H0.numCols) {
        return false
    }

    // Check if H1 is square
    if (H1.numRows != H1.numCols) {
        return false
    }

    // Check if H0 and H1 have the same size
    if (H0.numRows != H1.numRows) {
        return false
    }

    val n = H0.numRows

    // Check if H0 + H1 forms a valid infinitesimal generator
    // (row sums should be 0)
    val Q = H0.add(H1)
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            rowSum += Q.get(i, j)
        }
        if (abs(rowSum) > prec * n) {
            return false
        }
    }

    // Check if all eigenvalues of H0 have negative real parts
    val eigenvalues = H0.eig()
    for (ev in eigenvalues) {
        if (ev.real > -prec) {
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
