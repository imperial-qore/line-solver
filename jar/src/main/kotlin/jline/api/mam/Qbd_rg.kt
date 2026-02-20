/**
 * @file Quasi-Birth-Death process R and G matrix computation
 *
 * Computes fundamental R and G matrices for QBD analysis of MAP/MAP/1 queues.
 * Essential for solving infinite-dimensional Markov chains with structured transitions.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.lib.smc.QBD_CR

/**
 * Compute R and G matrices for MAP/MAP/1 queue using QBD approach
 *
 * @param MAPa Arrival MAP (MatrixCell with D0 and D1)
 * @param MAPs Service MAP (MatrixCell with D0 and D1)
 * @param util Optional utilization parameter for scaling
 * @return QbdRgResult containing R, G, B, L, F, and U matrices
 */
data class QbdRgResult(
    val R: Matrix,
    val G: Matrix,
    val B: Matrix,
    val L: Matrix,
    val F: Matrix,
    val U: Matrix
)

fun qbd_rg(MAPa: MatrixCell, MAPs: MatrixCell, util: Double? = null): QbdRgResult {
    val na = MAPa[0].numRows
    val ns = MAPs[0].numRows

    var scaledMAPs = MAPs
    if (util != null) {
        val lambdaA = map_lambda(MAPa)
        scaledMAPs = map_scale(MAPs, util / lambdaA)
    }

    // Construct QBD matrices using Kronecker sum for L (matching MATLAB krons)
    val IA = Matrix.eye(na)
    val IS = Matrix.eye(ns)
    val F = MAPa[1].kron(IS)
    val L = MAPa[0].kron(IS).add(IA.kron(scaledMAPs[0]))  // Kronecker sum
    val B = IA.kron(scaledMAPs[1])

    // Solve for G, R, U using SMC library Cyclic Reduction
    val qbdResult = QBD_CR(B, L, F, null, null, null, null)
    val G = qbdResult["G"] ?: throw RuntimeException("QBD_CR failed to compute G matrix")
    val R = qbdResult["R"] ?: throw RuntimeException("QBD_CR failed to compute R matrix")
    val U = qbdResult["U"] ?: throw RuntimeException("QBD_CR failed to compute U matrix")

    return QbdRgResult(R, G, B, L, F, U)
}

/**
 * Result from simplified QBD CR solver (kept for backward compatibility)
 */
data class QbdCrResult(val G: Matrix, val R: Matrix, val U: Matrix)

/**
 * QBD rg algorithms
 */
@Suppress("unused")
class QbdRgAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
