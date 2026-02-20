/**
 * @file Quasi-Birth-Death process MAP/MAP/1 queue analysis
 *
 * Analyzes MAP/MAP/1 queueing systems using QBD matrix analytic methods.
 * Computes performance measures including throughput, utilization, and queue length.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.lib.smc.QBD_CR
import jline.lib.smc.QBD_pi
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Analyze MAP/MAP/1 queue using QBD methods, computing performance measures
 *
 * @param MAPa Arrival MAP
 * @param MAPs Service MAP
 * @param util Optional utilization parameter
 * @return QbdMapMap1Result containing performance measures
 */
data class QbdMapMap1Result(
    val XN: Double,        // Throughput
    val QN: Double,        // Queue length
    val UN: Double,        // Utilization
    val pqueue: Matrix,    // Queue state probabilities
    val R: Matrix,         // Rate matrix
    val eta: Matrix,       // Caudal characteristic
    val G: Matrix,         // G matrix
    val A_1: Matrix,       // Downward transition block
    val A0: Matrix,        // Local transition block
    val A1: Matrix,        // Upward transition block
    val U: Matrix,         // Matrix U
    val MAPs: MatrixCell   // Scaled service process
)

fun qbd_mapmap1(MAPa: MatrixCell, MAPs: MatrixCell, util: Double? = null): QbdMapMap1Result {
    val na = MAPa[0].numRows
    val ns = MAPs[0].numRows

    // Scale service MAP if utilization is specified
    var scaledMAPs = MAPs
    if (util != null) {
        val lambdaA = map_lambda(MAPa)
        scaledMAPs = map_scale(MAPs, util / lambdaA)
    }

    val lambdaA = map_lambda(MAPa)
    val lambdaS = map_lambda(scaledMAPs)
    val actualUtil = lambdaA / lambdaS

    // Construct QBD matrices matching MATLAB qbd_mapmap1.m
    val IA = Matrix.eye(na)
    val IS = Matrix.eye(ns)
    val A1 = MAPa[1].kron(IS)                                          // Forward (arrivals)
    val A0 = MAPa[0].kron(IS).add(IA.kron(scaledMAPs[0]))             // Local (Kronecker sum)
    val A_1 = IA.kron(scaledMAPs[1])                                    // Backward (services)
    val A0bar = MAPa[0].kron(IS)                                        // Boundary local

    // Solve QBD using Cyclic Reduction: QBD_CR(A_1, A0, A1)
    val qbdResult = QBD_CR(A_1, A0, A1, null, null, null, null)
    val G = qbdResult["G"] ?: throw RuntimeException("QBD_CR failed to compute G matrix")
    val R = qbdResult["R"] ?: throw RuntimeException("QBD_CR failed to compute R matrix")
    val U = qbdResult["U"] ?: throw RuntimeException("QBD_CR failed to compute U matrix")

    val n = na * ns

    // Compute queue length distribution using QBD_pi
    // MATLAB: pqueue = QBD_pi(A_1, A0bar, R, 'MaxNumComp', 1e2)
    val pi = QBD_pi(A_1, A0bar, R, 100, 0, null, 0)
    val numLevels = pi.numCols / n

    // Reshape pqueue: each row is a level
    var pqueue = Matrix(numLevels, n)
    for (i in 0 until numLevels) {
        for (j in 0 until n) {
            pqueue[i, j] = pi[0, i * n + j]
        }
    }

    // Retry with more components if needed (MATLAB lines 75-77)
    var busyProb = 0.0
    for (i in 1 until numLevels) {
        for (j in 0 until n) {
            busyProb += pqueue[i, j]
        }
    }
    if (busyProb < actualUtil * 0.99) {
        val pi2 = QBD_pi(A_1, A0bar, R, 20000, 0, null, 0)
        val numLevels2 = pi2.numCols / n
        pqueue = Matrix(numLevels2, n)
        for (i in 0 until numLevels2) {
            for (j in 0 until n) {
                pqueue[i, j] = pi2[0, i * n + j]
            }
        }
    }

    // Compute performance measures (matching MATLAB lines 79-88)
    var UN: Double
    var QN: Double

    if (na == 1 && ns == 1) {
        UN = 1.0 - pqueue[0, 0]
        QN = 0.0
        for (i in 0 until pqueue.numRows) {
            QN += i.toDouble() * pqueue[i, 0]
        }
    } else {
        var sum0 = 0.0
        for (j in 0 until n) {
            sum0 += pqueue[0, j]
        }
        UN = 1.0 - sum0
        QN = 0.0
        for (i in 0 until pqueue.numRows) {
            var levelProb = 0.0
            for (j in 0 until n) {
                levelProb += pqueue[i, j]
            }
            QN += i.toDouble() * levelProb
        }
    }

    val XN = lambdaA

    // Compute eta as spectral radius of R stored in 1x1 matrix
    val eta = Matrix(1, 1)
    eta[0, 0] = spectralRadiusMapmap1(R)

    // Build result MAPs cell
    val resultMAPs = MatrixCell(2)
    resultMAPs[0] = scaledMAPs[0]
    resultMAPs[1] = scaledMAPs[1]

    return QbdMapMap1Result(XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, resultMAPs)
}

/**
 * Compute spectral radius of a matrix
 */
private fun spectralRadiusMapmap1(A: Matrix): Double {
    val n = A.numRows
    val dm = org.ejml.data.DMatrixRMaj(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            dm[i, j] = A[i, j]
        }
    }
    val evd = org.ejml.dense.row.factory.DecompositionFactory_DDRM.eig(n, false)
    evd.decompose(dm)
    var maxAbs = 0.0
    for (i in 0 until evd.numberOfEigenvalues) {
        val ev = evd.getEigenvalue(i)
        val absVal = Math.sqrt(ev.real * ev.real + ev.imaginary * ev.imaginary)
        if (absVal > maxAbs) {
            maxAbs = absVal
        }
    }
    return maxAbs
}

/**
 * QBD mapmap1 algorithms
 */
@Suppress("unused")
class QbdMapmap1 {
    companion object {
        // Class documentation marker for Dokka
    }
}
