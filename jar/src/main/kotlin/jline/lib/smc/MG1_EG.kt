/**
 * @file M/G/1-type Explicit G computation
 *
 * Determines G directly when rank(A0)=1 for M/G/1-type Markov chains.
 * This is an optimization that avoids iterative computation when possible.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Determines G directly if rank(A0)=1 for M/G/1-type Markov chains.
 *
 * @param A The block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns
 * @param verbose When true, prints residual error
 * @return G matrix if explicit solution exists, null otherwise
 */
fun mg1_eg(A: Matrix, verbose: Boolean = false): Matrix? {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Compute sumA and beta
    var sumA = A.extractCols(dega * m, (dega + 1) * m)
    var beta = sumA.sumRows()

    for (i in dega - 1 downTo 1) {
        sumA = sumA.add(A.extractCols(i * m, (i + 1) * m))
        beta = beta.add(sumA.sumRows())
    }
    sumA = sumA.add(A.extractCols(0, m))

    // Compute stationary distribution
    val theta = stat(sumA)
    val drift = theta.mult(beta)[0, 0]

    var G: Matrix? = null

    if (drift < 1) {
        // Positive recurrent case
        val A0 = A.extractCols(0, m)
        if (matrixRank(A0) == 1) {
            // A0 = alpha * beta, find first nonzero row
            val rowSums = A0.sumRows()
            var temp = -1
            for (i in 0 until m) {
                if (rowSums[i, 0] > 0) {
                    temp = i
                    break
                }
            }
            if (temp >= 0) {
                val rowSum = A0.extractRows(temp, temp + 1).elementSum()
                val betaVec = A0.extractRows(temp, temp + 1).scale(1.0 / rowSum)
                G = Matrix.ones(m, 1).mult(betaVec)
            }
        }
    } else if (drift > 1) {
        // Transient case
        val A0 = A.extractCols(0, m)
        if (matrixRank(A0) == 1) {
            // Transform A matrices
            val Atransformed = Matrix(m, A.numCols)
            for (i in 0..dega) {
                val Ai = A.extractCols(i * m, (i + 1) * m)
                // A_i' = diag(theta^-1) * A_i^T * diag(theta)
                val thetaInvDiag = Matrix.diag(*theta.scale(-1.0).exp().getRow(0).toArray1D())
                val thetaDiag = Matrix.diag(*theta.getRow(0).toArray1D())
                val transformedAi = thetaInvDiag.mult(Ai.transpose()).mult(thetaDiag)
                for (r in 0 until m) {
                    for (c in 0 until m) {
                        Atransformed[r, i * m + c] = transformedAi[r, c]
                    }
                }
            }

            // Compute etahat using GIM1_Caudal
            val etahat = gim1_caudal(Atransformed)

            if (etahat != null) {
                // temp = A_max
                var temp = Atransformed.extractCols(dega * m, (dega + 1) * m)
                for (i in dega - 1 downTo 1) {
                    temp = temp.mult(etahat).add(Atransformed.extractCols(i * m, (i + 1) * m))
                }

                // G = diag(theta^-1) * (A0 * (I - temp)^-1)^T * diag(theta)
                val A0t = Atransformed.extractCols(0, m)
                val ImT = Matrix.eye(m).sub(temp)
                val invImT = ImT.inv()
                val product = A0t.mult(invImT).transpose()

                val thetaInvDiag2 = Matrix.diag(*theta.scale(-1.0).exp().getRow(0).toArray1D())
                val thetaDiag2 = Matrix.diag(*theta.getRow(0).toArray1D())
                G = thetaInvDiag2.mult(product).mult(thetaDiag2)
            }
        }
    }

    return G
}

/**
 * Compute matrix rank (simplified estimation)
 */
private fun matrixRank(A: Matrix): Int {
    val svd = A.svd()
    val singularValues = svd.s
    val tol = 1e-10 * singularValues[0, 0]
    var rank = 0
    for (i in 0 until minOf(A.numRows, A.numCols)) {
        if (singularValues[i, 0] > tol) {
            rank++
        }
    }
    return rank
}

/**
 * Placeholder for GIM1_Caudal function
 * TODO: Implement full GIM1_Caudal
 */
private fun gim1_caudal(A: Matrix): Matrix? {
    // This requires full GIM1 implementation
    // For now, return null to fall back to iterative methods
    return null
}
