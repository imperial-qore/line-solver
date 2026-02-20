/**
 * @file Quasi-Birth-Death process RAP/RAP/1 queue analysis
 *
 * Analyzes RAP/RAP/1 queueing systems using QBD methods with rational arrival processes.
 * Extends QBD analysis to more general arrival and service process models.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.lib.smc.QBD_CR
import jline.lib.smc.QBD_pi
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.ejml.dense.row.factory.DecompositionFactory_DDRM
import org.ejml.data.DMatrixRMaj

/**
 * Analyze RAP/RAP/1 queue using QBD methods
 *
 * @param RAPa Arrival RAP (Rational Arrival Process)
 * @param RAPs Service RAP (Rational Service Process)
 * @param util Optional utilization parameter
 * @return QbdRapRap1Result containing performance measures
 */
data class QbdRapRap1Result(
    val XN: Double,        // Throughput
    val QN: Double,        // Queue length
    val UN: Double,        // Utilization
    val pqueue: Matrix,    // Queue state probabilities
    val R: Matrix,         // Rate matrix
    val eta: Matrix,       // Caudal characteristic (scalar stored as 1x1)
    val G: Matrix,         // G matrix
    val B: Matrix,         // Backward matrix
    val L: Matrix,         // Local matrix
    val F: Matrix          // Forward matrix
)

fun qbd_raprap1(RAPa: MatrixCell, RAPs: MatrixCell, util: Double? = null): QbdRapRap1Result {
    val na = RAPa[0].numRows
    val ns = RAPs[0].numRows

    var scaledRAPs = RAPs
    if (util != null) {
        val lambdaA = map_lambda(RAPa)
        scaledRAPs = map_scale(RAPs, util / lambdaA)
    }

    val lambdaA = map_lambda(RAPa)
    val lambdaS = map_lambda(scaledRAPs)
    val actualUtil = lambdaA / lambdaS

    // Construct QBD matrices for RAP/RAP/1 using Kronecker sum for L
    val IA = Matrix.eye(na)
    val IS = Matrix.eye(ns)
    val F = RAPa[1].kron(IS)
    val L = RAPa[0].kron(IS).add(IA.kron(scaledRAPs[0]))  // Kronecker sum
    val B = IA.kron(scaledRAPs[1])

    // Boundary matrices
    val B0 = B.copy()
    val B1 = RAPa[0].kron(IS)

    // Solve using QBD CR with RAP components (RAPComp=1)
    val qbdResult = QBD_CR(B, L, F, null, null, null, 1)
    val R = qbdResult["R"] ?: throw RuntimeException("QBD_CR failed to compute R matrix")
    val G = qbdResult["G"] ?: throw RuntimeException("QBD_CR failed to compute G matrix")

    // Compute stationary distribution using QBD_pi with RAP components
    val pi = QBD_pi(B0, B1, R, 100, 0, null, 1)
    val n = na * ns
    val numLevels = pi.numCols / n

    // Reshape pqueue: each row is a level's distribution summed over phases
    val pqueue = Matrix(numLevels, n)
    for (i in 0 until numLevels) {
        for (j in 0 until n) {
            pqueue[i, j] = pi[0, i * n + j]
        }
    }

    // Compute eta = max(abs(eig(R)))
    val etaVal = spectralRadius(R)
    val eta = Matrix(1, 1)
    eta[0, 0] = etaVal

    // Compute performance measures
    val XN = lambdaA
    var UN: Double
    var QN: Double

    if (na == 1 && ns == 1) {
        UN = 1.0 - pqueue[0, 0]
    } else {
        var sum0 = 0.0
        for (j in 0 until n) {
            sum0 += pqueue[0, j]
        }
        UN = 1.0 - sum0
    }

    // QN = sum over levels k of k * sum_j(pqueue(k,j))
    QN = 0.0
    for (i in 0 until numLevels) {
        var levelProb = 0.0
        for (j in 0 until n) {
            levelProb += pqueue[i, j]
        }
        QN += i.toDouble() * levelProb
    }

    return QbdRapRap1Result(XN, QN, UN, pqueue, R, eta, G, B, L, F)
}

/**
 * Compute spectral radius (maximum absolute eigenvalue) of a matrix
 */
private fun spectralRadius(A: Matrix): Double {
    val n = A.numRows
    val dm = DMatrixRMaj(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            dm[i, j] = A[i, j]
        }
    }
    val evd = DecompositionFactory_DDRM.eig(n, false)
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
 * QBD raprap1 algorithms
 */
@Suppress("unused")
class QbdRaprap1 {
    companion object {
        // Class documentation marker for Dokka
    }
}
