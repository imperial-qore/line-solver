/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs

/**
 * Result class for CanonicalFromDPH2 containing both beta and B.
 */
data class DPH2Representation(val beta: Matrix, val B: Matrix)

/**
 * Returns the canonical form of an order-2 discrete phase-type distribution.
 *
 * @param alpha Initial vector of the discrete phase-type distribution
 * @param A Transition probability matrix of the discrete phase-type distribution
 * @param prec Numerical precision for checking the input, default value is 1e-14
 * @return The DPH2Representation containing beta (canonical initial vector) and B (canonical transition matrix)
 */
fun canonicalFromDPH2(alpha: Matrix, A: Matrix, prec: Double = 1e-14): DPH2Representation {
    if (A.numRows != 2 || A.numCols != 2) {
        throw IllegalArgumentException("CanonicalFromDPH2: Dimension must be 2!")
    }

    if (!checkMGRepresentation(alpha, A, prec)) {
        throw IllegalArgumentException("CanonicalFromDPH2: Input isn't a valid MG distribution!")
    }

    // Get eigenvalues sorted by absolute value (descending)
    val eigenvalues: List<Complex> = A.eig()
    val lambda = eigenvalues.sortedByDescending { ev: Complex -> abs(ev.real * ev.real + ev.imaginary * ev.imaginary) }
        .map { ev: Complex -> ev.real }

    // e = [1; 1]
    // p1 = alpha * (e - A*e)
    val N = A.numRows
    val e = Matrix(N, 1)
    for (i in 0 until N) {
        e[i, 0] = 1.0
    }
    val Ae = A.mult(e)
    val eMinusAe = e.sub(Ae)
    val p1 = alpha.mult(eMinusAe)[0, 0]

    val beta: Matrix
    val B: Matrix

    if (lambda[0] > 0 && lambda[1] > 0 && abs(lambda[0] - lambda[1]) > prec) {
        // Case: Two distinct positive eigenvalues
        val d1 = (1 - lambda[0]) * (1 - p1 - lambda[1]) / (lambda[0] - lambda[1])
        val d2 = p1 - d1

        beta = Matrix(1, 2)
        beta[0, 0] = d1 * (lambda[0] - lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1]))
        beta[0, 1] = (d1 + d2) / (1 - lambda[1])

        B = Matrix(2, 2)
        B[0, 0] = lambda[0]
        B[0, 1] = 1 - lambda[0]
        B[1, 0] = 0.0
        B[1, 1] = lambda[1]
    } else if (lambda[0] > 0 && abs(lambda[0] - lambda[1]) <= prec) {
        // Case: Two equal positive eigenvalues
        val d2 = p1
        val d1 = (1 - lambda[0]) * (1 - d2 - lambda[0]) / lambda[0]

        beta = Matrix(1, 2)
        beta[0, 0] = d1 * lambda[0] / ((1 - lambda[0]) * (1 - lambda[0]))
        beta[0, 1] = d2 / (1 - lambda[0])

        B = Matrix(2, 2)
        B[0, 0] = lambda[0]
        B[0, 1] = 1 - lambda[0]
        B[1, 0] = 0.0
        B[1, 1] = lambda[0]
    } else if (lambda[0] > 0) {
        // Case: One positive, one non-positive eigenvalue
        val d1 = (1 - lambda[0]) * (1 - p1 - lambda[1]) / (lambda[0] - lambda[1])
        val d2 = p1 - d1

        beta = Matrix(1, 2)
        beta[0, 0] = (d1 * lambda[0] + d2 * lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1]))
        beta[0, 1] = (d1 + d2) * (1 - lambda[0] - lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1]))

        B = Matrix(2, 2)
        B[0, 0] = lambda[0] + lambda[1]
        B[0, 1] = 1 - lambda[0] - lambda[1]
        B[1, 0] = lambda[0] * lambda[1] / (lambda[0] + lambda[1] - 1)
        B[1, 1] = 0.0
    } else {
        throw IllegalArgumentException("CanonicalFromDPH2: Cannot convert to canonical form!")
    }

    return DPH2Representation(beta, B)
}

/**
 * Overload for DoubleArray alpha.
 */
fun canonicalFromDPH2(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): DPH2Representation {
    return canonicalFromDPH2(Matrix(alpha), A, prec)
}
