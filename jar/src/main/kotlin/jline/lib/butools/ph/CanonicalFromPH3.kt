/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix

/**
 * Returns the canonical form of an order-3 phase-type distribution.
 *
 * @param alpha Initial vector of the phase-type distribution.
 * @param A Transient generator of the phase-type distribution.
 * @param prec Numerical precision, default value is 1e-10.
 * @return The PH3Representation containing beta (canonical initial vector) and B (canonical generator).
 *
 * Note: This procedure calculates 5 moments of the input and calls 'ph3From5Moments'.
 */
fun canonicalFromPH3(alpha: Matrix, A: Matrix, prec: Double = 1e-10): PH3Representation {
    if (A.numRows != 3) {
        throw IllegalArgumentException("CanonicalFromPH3: Dimension is not 3!")
    }

    if (!checkMERepresentation(alpha, A)) {
        throw IllegalArgumentException("CanonicalFromPH3: Input isn't a valid ME distribution!")
    }

    val moms = momentsFromME(alpha, A, 5)
    return ph3From5Moments(moms, prec)
}

/**
 * Overload for DoubleArray alpha.
 */
fun canonicalFromPH3(alpha: DoubleArray, A: Matrix, prec: Double = 1e-10): PH3Representation {
    return canonicalFromPH3(Matrix(alpha), A, prec)
}
