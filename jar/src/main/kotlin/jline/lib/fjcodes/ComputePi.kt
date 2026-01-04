package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Compute steady-state distribution for 2-replica case
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Compute pseudo-inverse of a matrix using SVD
 * Handles rank-deficient matrices by zeroing out small singular values
 */
private fun pseudoInverse(A: Matrix, tol: Double = 1e-10): Matrix {
    val svdResult = A.svd()
    val U = svdResult.u  // Left singular vectors
    val S = svdResult.s  // Singular values as column vector
    val V = svdResult.v  // Right singular vectors

    // For compact SVD: U is mxk, S is kx1, V is nxk, where k = min(m,n)
    // We need Splus to be kxk diagonal matrix
    val m = A.getNumRows()
    val n = A.getNumCols()
    val k = Math.min(m, n)  // Rank of compact SVD
    val Splus = Matrix(k, k)  // Square diagonal matrix for compact SVD

    val numSingularValues = S.getNumRows()
    for (i in 0 until numSingularValues) {
        val sigma = S.get(i, 0)
        if (kotlin.math.abs(sigma) > tol) {
            Splus.set(i, i, 1.0 / sigma)
        }
        // else: leave as 0 (zero out small singular values)
    }

    // A+ = V * S+ * U^T
    return V.mult(Splus).mult(U.transpose())
}

/**
 * Solve overdetermined/rank-deficient system A * x = b using pseudo-inverse
 */
private fun solveLeastSquares(A: Matrix, b: Matrix): Matrix {
    val Aplus = pseudoInverse(A)
    return Aplus.mult(b)
}

/**
 * Arrival representation for FJ_codes
 */
data class FJArrival(
    val lambda: Double,     // Arrival rate
    val lambda0: Matrix,    // D0 matrix (transitions without arrivals)
    val lambda1: Matrix,    // D1 matrix (transitions with arrivals)
    val ArrChoice: Int      // Arrival type: 1=Exponential, other=MAP
)

/**
 * Result of computePi
 */
data class PiResult(
    val pi0: Matrix,   // Steady-state probability (row vector)
    val En1: Double    // Expected number in system
)

/**
 * Compute steady-state distribution for the case with 2 replicas
 *
 * Solves for the steady-state probability pi0 and expected number En1.
 *
 * @param T T-matrix from NARE solution
 * @param arrival Arrival process
 * @param services Service process for single subtask
 * @param service_h Service representation for 2-node job
 * @param C Capacity parameter
 * @param S State transition matrix
 * @param A_jump Jump matrix
 * @return PiResult with pi0 and En1
 */
fun computePi(
    T: Matrix,
    arrival: FJArrival,
    services: FJService,
    service_h: FJServiceH,
    C: Int,
    S: Matrix,
    A_jump: Matrix
): PiResult {
    if (services.SerChoice == 1) {
        // Exponential service case
        val ms = S.getNumCols()
        val S_notallbusy = constructNotAllBusy(C, services, service_h)
        val Q0 = kronsum(S_notallbusy, arrival.lambda0)

        val da = arrival.lambda0.getNumRows()

        // Solve Lyapunov equation: T*Igral + Igral*(kron(I, lambda0)) = -I
        // This is: A*X + X*B + C = 0 with A=T, B=kron(I, lambda0), C=-I
        val B = Matrix.eye(ms).kron(arrival.lambda0)
        val C = Matrix.eye(da * ms).scale(-1.0)
        val Igral = Matrix.lyap(T, B, C, null)

        // pi0mat = Igral * kron(A_jump, I) * inv(Q0) * kron(I, lambda1)
        // Use pseudo-inverse for Q0 (rate matrix, potentially singular)
        val Q0inv = pseudoInverse(Q0)
        val pi0mat = Igral
            .mult(A_jump.kron(Matrix.eye(da)))
            .mult(Q0inv)
            .mult(Matrix.eye(ms).kron(arrival.lambda1))

        // Solve for pi0: pi0 * (pi0mat - I) = 0, pi0 * 1 = 1
        // MATLAB: pi0 = [zeros(1,ms*da),-1]/[(pi0mat - eye(ms*da)),sum(inv(T),2)]
        // This is: pi0 * lhsMat = rhsVec
        //   where lhsMat = [(pi0mat - I), sum(inv(T),2)] is [ms*da × (ms*da+1)]
        //   and rhsVec = [zeros(1,ms*da), -1] is [1 × (ms*da+1)]

        val rhsVec = Matrix(1, ms * da + 1)
        for (i in 0 until ms * da) {
            rhsVec.set(0, i, 0.0)
        }
        rhsVec.set(0, ms * da, -1.0)

        val lhsMat = Matrix(ms * da, ms * da + 1)

        // First ms*da columns: (pi0mat - I)
        for (i in 0 until ms * da) {
            for (j in 0 until ms * da) {
                val value = if (i == j) pi0mat.get(i, j) - 1.0 else pi0mat.get(i, j)
                lhsMat.set(i, j, value)
            }
        }

        // Last column: sum(inv(T), 2) = row sums of inv(T)
        // Use pseudo-inverse for robustness
        val Tinv = pseudoInverse(T)
        for (i in 0 until ms * da) {
            var rowSum = 0.0
            for (j in 0 until ms * da) {
                rowSum += Tinv.get(i, j)
            }
            lhsMat.set(i, ms * da, rowSum)
        }

        // Solve: pi0 * lhsMat = rhsVec
        // Transpose to get: lhsMat^T * pi0^T = rhsVec^T
        // Use pseudo-inverse for rank-deficient system
        val pi0_T = solveLeastSquares(lhsMat.transpose(), rhsVec.transpose())
        val pi0 = pi0_T.transpose()

        // En1 = (1/sum(pi0)) * sum(pi0 * pi0mat)
        var sumPi0 = 0.0
        for (i in 0 until ms * da) {
            sumPi0 += pi0.get(0, i)
        }

        val pi0_times_pi0mat = pi0.mult(pi0mat)
        var sumResult = 0.0
        for (i in 0 until ms * da) {
            sumResult += pi0_times_pi0mat.get(0, i)
        }

        val En1 = (1.0 / sumPi0) * sumResult

        return PiResult(pi0, En1)

    } else {
        // General PH service case
        val result = constructSRK(C, services, service_h, S)
        val Se = result.Se
        val Sestar = result.Sestar
        val R0 = result.R0
        val Ke = result.Ke
        val Kc = result.Kc

        val dtmat = T.getNumRows() / arrival.lambda0.getNumCols()
        val dsexp = Se.getNumRows()

        // Extract submatrices
        val Sedash = Matrix.getSubMatrix(Se, dtmat, dsexp, dtmat, dsexp)
        val Rbusy = Matrix.getSubMatrix(R0, dtmat, dsexp, 0, R0.getNumCols())

        val da = arrival.lambda0.length()
        val Iidle_rows = dsexp - dtmat
        val Iidle = Matrix(Iidle_rows, da)
        for (i in 0 until Math.min(Iidle_rows, da)) {
            Iidle.set(i, i, 1.0)
        }

        // Q0 = kronsum(Sedash, lambda0)
        val Q0 = kronsum(Sedash, arrival.lambda0)

        // Qbusy = kron(Rbusy, lambda1)
        val Qbusy = Rbusy.kron(arrival.lambda1)

        // Bmap = -(Iidle / Qidle) * Qbusy
        // Use pseudo-inverse for Qidle (rate matrix)
        val Qidle = Q0
        val Qidle_inv = pseudoInverse(Qidle)
        val Bmap = Iidle.mult(Qidle_inv).scale(-1.0).mult(Qbusy)

        val Kemap = Ke.kron(Matrix.eye(da))
        val Kcmap = Kc.kron(Matrix.eye(da))

        // Solve Lyapunov equation for Igral
        // T*Igral + Igral*BB + Kemap*kron(Sestar, I) = 0
        val BB = Matrix.eye(Sestar.getNumCols()).kron(arrival.lambda0)
        val C_lyap = Kemap.mult(Sestar.kron(Matrix.eye(da)))
        val Igral = Matrix.lyap(T, BB, C_lyap, null)

        // pi0mat = Igral * Bmap * Kcmap
        val pi0mat = Igral.mult(Bmap).mult(Kcmap)

        // Solve for pi0 (same structure as above)
        val ms_da = dtmat * da

        val rhsVec = Matrix(1, ms_da + 1)
        for (i in 0 until ms_da) {
            rhsVec.set(0, i, 0.0)
        }
        rhsVec.set(0, ms_da, -1.0)

        val lhsMat = Matrix(ms_da, ms_da + 1)

        for (i in 0 until ms_da) {
            for (j in 0 until ms_da) {
                val value = if (i == j) pi0mat.get(i, j) - 1.0 else pi0mat.get(i, j)
                lhsMat.set(i, j, value)
            }
        }

        // Use pseudo-inverse for robustness
        val Tinv = pseudoInverse(T)
        for (i in 0 until ms_da) {
            var rowSum = 0.0
            for (j in 0 until ms_da) {
                rowSum += Tinv.get(i, j)
            }
            lhsMat.set(i, ms_da, rowSum)
        }

        // Solve: pi0 * lhsMat = rhsVec using pseudo-inverse
        val pi0_T = solveLeastSquares(lhsMat.transpose(), rhsVec.transpose())
        val pi0 = pi0_T.transpose()

        // En1 = (1/sum(pi0)) * sum(-pi0 * Igral * Iidle / Qidle * kron(I, lambda1))
        var sumPi0 = 0.0
        for (i in 0 until ms_da) {
            sumPi0 += pi0.get(0, i)
        }

        // Use pseudo-inverse for Qidle
        val term = pi0.scale(-1.0)
            .mult(Igral)
            .mult(Iidle)
            .mult(Qidle_inv)
            .mult(Matrix.eye(Sedash.getNumRows()).kron(arrival.lambda1))

        var sumTerm = 0.0
        for (i in 0 until term.getNumCols()) {
            sumTerm += term.get(0, i)
        }

        val En1 = (1.0 / sumPi0) * sumTerm

        return PiResult(pi0, En1)
    }
}
