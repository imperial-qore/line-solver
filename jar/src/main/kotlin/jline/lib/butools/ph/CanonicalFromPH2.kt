/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix

/**
 * Returns the canonical form of an order-2 phase-type distribution.
 *
 * @param alpha Initial vector of the phase-type distribution.
 * @param A Transient generator of the phase-type distribution.
 * @param prec Numerical precision, default value is 1e-14.
 * @return The PH2Representation containing beta (canonical initial vector) and B (canonical generator).
 *
 * Note: This procedure calculates 3 moments of the input and calls 'ph2From3Moments'.
 */
@JvmOverloads
fun canonicalFromPH2(alpha: Matrix, A: Matrix, prec: Double = 1e-14): PH2Representation {
    if (A.numRows != 2) {
        throw IllegalArgumentException("CanonicalFromPH2: Dimension is not 2!")
    }

    if (!checkMERepresentation(alpha, A)) {
        throw IllegalArgumentException("CanonicalFromPH2: Input isn't a valid ME distribution!")
    }

    val moms = momentsFromME(alpha, A, 3)
    return ph2From3Moments(moms, prec)
}

/**
 * Overload for DoubleArray alpha.
 */
@JvmOverloads
fun canonicalFromPH2(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): PH2Representation {
    return canonicalFromPH2(Matrix(alpha), A, prec)
}
