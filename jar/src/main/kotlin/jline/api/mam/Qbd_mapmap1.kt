/**
 * @file Quasi-Birth-Death process MAP/MAP/1 queue analysis
 * 
 * Analyzes MAP/MAP/1 queueing systems using QBD matrix analytic methods.
 * Computes performance measures including throughput, utilization, and queue length.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.GlobalConstants.NegInf
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
    val eta: Matrix,       // Auxiliary matrix
    val G: Matrix          // G matrix
)

fun qbd_mapmap1(MAPa: MatrixCell, MAPs: MatrixCell, util: Double? = null): QbdMapMap1Result {
    // Following the MATLAB MAPMAP1.m implementation structure
    val na = MAPa[0].numRows
    val ns = MAPs[0].numRows
    val IA = Matrix.eye(na)
    val IS = Matrix.eye(ns)
    
    // Scale service MAP if utilization is specified
    var scaledMAPs = MAPs
    if (util != null) {
        val lambdaA = map_lambda(MAPa)
        scaledMAPs = map_scale(MAPs, util / lambdaA)  
    }
    
    val lambdaA = map_lambda(MAPa)
    val lambdaS = map_lambda(scaledMAPs) 
    val actualUtil = lambdaA / lambdaS
    
    if (actualUtil >= 1.0) {
        throw IllegalArgumentException("System is unstable: utilization >= 1.0 (ρ = $actualUtil)")
    }
    
    // Construct QBD matrices as in MATLAB MAPMAP1.m
    // B = kron(IA, S1) - backward transitions (service completions)
    // L = kron(D0, IS) + kron(IA, S0) - local transitions  
    // F = kron(D1, IS) - forward transitions (arrivals)
    // L0 = kron(D0, IS) - boundary local transitions
    val B = IA.kron(scaledMAPs[1])
    val L = MAPa[0].kron(IS).add(IA.kron(scaledMAPs[0]))
    val F = MAPa[1].kron(IS)
    val L0 = MAPa[0].kron(IS)
    
    // Solve QBD using existing qbd_rg implementation
    // Create temporary MatrixCells for the QBD structure
    val MAPa_temp = MatrixCell(2)
    MAPa_temp[0] = L  // Local transitions correspond to D0
    MAPa_temp[1] = F  // Forward transitions correspond to D1
    
    val MAPs_temp = MatrixCell(2) 
    MAPs_temp[0] = Matrix.zeros(B.numRows, B.numCols)  // No hidden service transitions
    MAPs_temp[1] = B  // Backward transitions correspond to service completions
    
    val rgResult = qbd_rg(MAPa_temp, MAPs_temp, null)
    val R = rgResult.R
    val G = rgResult.G
    val n = L0.numRows
    val I = Matrix.eye(n)
    
    // Convert to discrete time if needed (following QBDSolve.m logic)
    var B_dt = B.copy()
    var L0_dt = L0.copy()
    
    // Check if continuous time (diagonal elements are negative)
    var isContinuousTime = false
    var maxDiag = NegInf
    for (i in 0 until L0.numRows) {
        val diagVal = L0[i, i]
        if (diagVal < 0) {
            isContinuousTime = true
        }
        if (diagVal > maxDiag) {
            maxDiag = diagVal
        }
    }
    
    if (isContinuousTime) { // continuous time case
        val lamb = -maxDiag  // max(-diag(L0)) - most negative diagonal element
        B_dt = B.scale(1.0 / lamb)
        L0_dt = L0.scale(1.0 / lamb).add(I)
    }
    
    // Compute boundary probabilities following MATLAB QBD_pi.m algorithm
    // π₀ is computed using stat(B1 + R*B0) where B1=L0_dt, B0=B_dt
    val boundaryMatrix = L0_dt.add(R.mult(B_dt))
    
    // Find stationary distribution of boundaryMatrix
    // Using nullspace approach: (boundaryMatrix)^T * pi = 0
    val boundaryT = boundaryMatrix.transpose()
    val eyeMinusBT = Matrix.eye(n).sub(boundaryT)
    
    // Solve via pseudo-inverse or add normalization constraint
    val A = Matrix.zeros(n + 1, n)
    val b = Matrix.zeros(n + 1, 1)
    
    // Set up system: [I - boundaryMatrix^T; ones^T] * pi0 = [0; 1]
    for (i in 0 until n) {
        for (j in 0 until n) {
            A[i, j] = eyeMinusBT[i, j]
        }
    }
    // Normalization constraint
    for (j in 0 until n) {
        A[n, j] = 1.0
    }
    b[n, 0] = 1.0
    
    // Solve using least squares
    val AtA = A.transpose().mult(A)
    val Atb = A.transpose().mult(b)
    val pi0Vec = AtA.inv().mult(Atb)
    
    val piBoundary = Matrix.zeros(1, n)
    for (i in 0 until n) {
        piBoundary[0, i] = pi0Vec[i, 0]
    }
    
    // Normalize following MATLAB: pi0 = pi0 / (pi0 * (I-R)^{-1} * ones)
    val eyeMinusR = I.sub(R)
    val eyeMinusRInv = eyeMinusR.inv()
    val nr = piBoundary.mult(eyeMinusRInv).mult(Matrix.ones(n, 1)).toDouble()
    val pi0Normalized = piBoundary.scale(1.0 / nr)
    
    // Performance measures
    val XN = lambdaA  // Throughput equals arrival rate in steady state
    val UN = actualUtil  // Utilization
    
    // Mean queue length using standard QBD formula: sum(pi0 * (I-R)^{-2} * R)
    val QN = pi0Normalized.mult(eyeMinusRInv).mult(eyeMinusRInv).mult(R).mult(Matrix.ones(n, 1)).toDouble()
    
    // Create truncated queue probabilities (first 10 levels)
    val maxLevels = 10
    val pqueue = Matrix.zeros(1, n * maxLevels)
    
    // Level 0 probabilities
    for (i in 0 until n) {
        pqueue[0, i] = pi0Normalized[0, i] 
    }
    
    // Higher level probabilities using π_k = π_0 * R^k
    var Rk = Matrix.eye(n)
    for (k in 1 until maxLevels) {
        Rk = Rk.mult(R)
        val piK = pi0Normalized.mult(Rk)
        for (i in 0 until n) {
            pqueue[0, k * n + i] = piK[0, i]
        }
    }
    
    // Compute auxiliary eta matrix (fundamental matrix)
    val eta = eyeMinusRInv.copy()
    
    return QbdMapMap1Result(XN, QN, UN, pqueue, R, eta, G)
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