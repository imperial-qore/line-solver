package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Compute waiting time distribution for Fork-Join queue
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result of returnWait
 */
data class ReturnWaitResult(
    val wait_alpha: Matrix,    // Initial probability for waiting time PH
    val wait_Smat: Matrix,      // Sub-generator for waiting time PH
    val prob_wait: Double,      // Probability of waiting
    val alfa: Matrix            // Steady-state distribution
)

/**
 * Compute Phase-Type representation of stationary waiting time
 *
 * Calculates the phase-type representation of the stationary waiting time
 * for a K-node Fork-Join queue.
 *
 * @param En1 Expected number in system
 * @param pi0 Steady-state probability vector
 * @param T T-matrix
 * @param phi Row sum vector
 * @param sum_Ajump Row sum of A_jump matrix
 * @return ReturnWaitResult with waiting time distribution
 */
fun returnWait(
    En1: Double,
    pi0: Matrix,
    T: Matrix,
    phi: Matrix,
    sum_Ajump: Matrix
): ReturnWaitResult {
    // alfa = -pi0 / T
    val alfa = pi0.scale(-1.0).rightMatrixDivide(T)

    val ds = phi.getNumRows()

    // rhos = diag(phi * alfa)' / (alfa * phi)
    // phi is 101x1 column, alfa is 1x101 row
    // phi * alfa is outer product: 101x1 mult 1x101 = 101x101
    val phiAlfa = phi.mult(alfa)  // Outer product, not phi * alfa'
    val alfaPhi = alfa.mult(phi)  // Inner product: 1x101 mult 101x1 = scalar
    val alfaPhiSum = alfaPhi.get(0, 0)

    val rhos = Matrix(1, ds)
    for (i in 0 until ds) {
        rhos.set(0, i, phiAlfa.get(i, 0) / alfaPhiSum)
    }

    // En0 = alfa * sum_Ajump / sum(pi0)
    val alfaSumAjump = alfa.mult(sum_Ajump).get(0, 0)
    val sumPi0 = pi0.elementSum()
    val En0 = alfaSumAjump / sumPi0

    // prob_wait = (En0 - 1) / (En0 - 1 + En1)
    val prob_wait = (En0 - 1.0) / (En0 - 1.0 + En1)

    // wait_alpha = prob_wait * rhos
    val wait_alpha = rhos.scale(prob_wait)

    // wait_Smat: G matrix
    val wait_Smat = Matrix(ds, ds)
    for (i in 0 until ds) {
        for (j in 0 until ds) {
            val gij = alfa.get(0, j) * T.get(j, i) / alfa.get(0, i)
            wait_Smat.set(i, j, gij)
        }
    }

    return ReturnWaitResult(wait_alpha, wait_Smat, prob_wait, alfa)
}

/**
 * Compute element-wise sum of matrix (helper extension)
 */
private fun Matrix.elementSum(): Double {
    var sum = 0.0
    for (i in 0 until getNumRows()) {
        for (j in 0 until getNumCols()) {
            sum += get(i, j)
        }
    }
    return sum
}
