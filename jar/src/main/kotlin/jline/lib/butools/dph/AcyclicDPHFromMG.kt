/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.lib.butools.reptrans.similarityMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs

/**
 * Transforms a matrix-geometric representation to an acyclic DPH representation
 * of the same size, if possible.
 *
 * @param alpha Initial vector of the distribution
 * @param A Matrix parameter of the distribution
 * @param prec Vector and matrix entries smaller than the precision are considered to be zeros.
 *        The default value is 1e-14.
 * @return The MGRepresentation containing beta (initial probability vector of acyclic DPH) and
 *         B (transition probability matrix of acyclic DPH)
 *
 * Note: Contrary to 'AcyclicPHFromME' of the 'ph' package, this procedure is not able to extend
 *       the size in order to obtain a Markovian initial vector.
 *       Raises an error if A has complex eigenvalues. In this case the transformation to an
 *       acyclic representation is not possible.
 */
fun acyclicDPHFromMG(alpha: Matrix, A: Matrix, prec: Double = 1e-14): MGRepresentation {
    if (!checkMGRepresentation(alpha, A, prec)) {
        throw IllegalArgumentException("AcyclicDPHFromMG: Input isn't a valid MG distribution!")
    }

    val N = A.numRows

    // Get eigenvalues
    val eigenvalues: List<Complex> = A.eig()

    // Check for complex eigenvalues
    for (ev in eigenvalues) {
        if (abs(ev.imaginary) > prec) {
            throw IllegalArgumentException("AcyclicDPHFromMG: The input matrix has complex eigenvalue!")
        }
    }

    // Sort eigenvalues by real part (descending)
    val lambda = eigenvalues.map { ev: Complex -> ev.real }.sortedDescending()

    // Create target acyclic matrix: mx = diag(lambda) + diag(1-lambda(1:end-1), 1)
    val mx = Matrix.zeros(N, N)
    for (i in 0 until N) {
        mx[i, i] = lambda[i]
    }
    for (i in 0 until N - 1) {
        mx[i, i + 1] = 1.0 - lambda[i]
    }

    // Find similarity transformation matrix
    val T = similarityMatrix(A, mx)

    // beta = alpha * T
    val beta = alpha.mult(T)
    val B = mx

    // Verify the result is a valid DPH representation
    if (!checkDPHRepresentation(beta, B, prec)) {
        throw IllegalArgumentException("AcyclicDPHFromMG: No acyclic representation found!")
    }

    return MGRepresentation(beta, B)
}

/**
 * Overload for DoubleArray alpha.
 */
fun acyclicDPHFromMG(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): MGRepresentation {
    return acyclicDPHFromMG(Matrix(alpha), A, prec)
}
