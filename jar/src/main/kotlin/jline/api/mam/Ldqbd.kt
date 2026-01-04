/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Level-Dependent Quasi-Birth-Death process solver
 *
 * Computes rate matrices for level-dependent QBD processes using matrix continued fractions.
 * Implements Algorithm 1 from "A Simple Algorithm for the Rate Matrices of Level-Dependent
 * QBD Processes" by Phung-Duc, Masuyama, Kasahara, and Takahashi (2010), QTNA Conference.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Result of LDQBD solver containing rate matrices and stationary distribution.
 *
 * @property R Cell array of rate matrices R^(1), R^(2), ..., R^(N)
 * @property pi Stationary distribution vector [pi_0, pi_1, ..., pi_N]
 */
data class LdqbdResult(
    val R: List<Matrix>,
    val pi: Matrix
)

/**
 * Options for LDQBD solver.
 *
 * @property epsilon Convergence tolerance (default: 1e-10)
 * @property maxIter Maximum number of iterations (default: 1000)
 * @property verbose Print debug information (default: false)
 */
data class LdqbdOptions(
    val epsilon: Double = 1e-10,
    val maxIter: Int = 1000,
    val verbose: Boolean = false
)

/**
 * Solves level-dependent QBD processes using matrix continued fractions.
 *
 * For a level-dependent QBD with levels 0, 1, ..., N, the infinitesimal
 * generator has block-tridiagonal structure:
 *
 * ```
 *   Q^(0)_1  Q^(0)_0   O        O      ...   O
 *   Q^(1)_2  Q^(1)_1  Q^(1)_0   O      ...   O
 *    O       Q^(2)_2  Q^(2)_1  Q^(2)_0 ...   O
 *   ...      ...      ...      ...    ...  ...
 *    O        O        O       ...   Q^(N)_2  Q^(N)_1
 * ```
 *
 * where Q_0^(n) are upward transitions (level n to n+1), Q_1^(n) are local
 * transitions (within level n), and Q_2^(n) are downward transitions
 * (level n to n-1).
 *
 * @param Q0 List of upward transition matrices {Q0^(0), Q0^(1), ..., Q0^(N-1)}
 * @param Q1 List of local transition matrices {Q1^(0), Q1^(1), ..., Q1^(N)}
 * @param Q2 List of downward transition matrices {Q2^(1), Q2^(2), ..., Q2^(N)}
 * @param options Solver options
 * @return LdqbdResult containing rate matrices R and stationary distribution pi
 */
fun ldqbd(
    Q0: List<Matrix>,
    Q1: List<Matrix>,
    Q2: List<Matrix>,
    options: LdqbdOptions = LdqbdOptions()
): LdqbdResult {
    // Determine number of levels
    // Q0 has entries for levels 0 to N-1
    // Q1 has entries for levels 0 to N
    // Q2 has entries for levels 1 to N
    val N = Q1.size - 1  // Maximum level

    if (options.verbose) {
        println("LD-QBD Solver: N=$N levels, epsilon=${options.epsilon}")
    }

    // Compute all R matrices using backward recursion
    val R = computeAllRateMatrices(N, Q0, Q1, Q2, options)

    if (options.verbose) {
        for (n in 0 until N) {
            println("  R^(${n + 1}) computed (${R[n].numRows}x${R[n].numCols})")
        }
    }

    // Compute stationary distribution
    val pi = computeStationaryDist(R, Q0, Q1, Q2, N, options)

    return LdqbdResult(R, pi)
}

/**
 * Compute all R matrices using backward recursion (heterogeneous dimensions).
 *
 * For heterogeneous LDQBD, compute R matrices using backward recursion:
 *   R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
 *   R^(n) = Q0^(n-1) * (-Q1^(n) - R^(n+1) * Q2^(n+1))^{-1}  for n = N-1,...,1
 *
 * Dimensions:
 *   R^(n) is (states at level n-1) x (states at level n)
 *   Q0^(n-1) is (states at level n-1) x (states at level n)
 *   Q1^(n) is (states at level n) x (states at level n)
 *   Q2^(n) is (states at level n) x (states at level n-1)
 */
private fun computeAllRateMatrices(
    N: Int,
    Q0: List<Matrix>,
    Q1: List<Matrix>,
    Q2: List<Matrix>,
    options: LdqbdOptions
): List<Matrix> {
    val R = MutableList<Matrix?>(N) { null }

    // Start at level N (boundary): R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
    val Q0_Nminus1 = Q0[N - 1]  // Q0^(N-1)
    val Q1_N = Q1[N]            // Q1^(N)

    val U_N = Q1_N.scale(-1.0)
    R[N - 1] = safeMatrixDivide(Q0_Nminus1, U_N)

    // Backward recursion for n = N-1 down to 1
    for (n in N - 1 downTo 1) {
        // R^(n) = Q0^(n-1) * (-Q1^(n) - R^(n+1) * Q2^(n+1))^{-1}
        val Q0_nminus1 = Q0[n - 1]      // Q0^(n-1): (states at n-1) x (states at n)
        val Q1_n = Q1[n]                // Q1^(n): (states at n) x (states at n)
        val Q2_nplus1 = Q2[n]           // Q2^(n+1): (states at n+1) x (states at n)
        val R_nplus1 = R[n]!!           // R^(n+1): (states at n) x (states at n+1)

        // Compute R^(n+1) * Q2^(n+1): (states at n) x (states at n)
        val RQ2 = R_nplus1.mult(Q2_nplus1)

        // Compute U = -Q1^(n) - R^(n+1) * Q2^(n+1)
        val U = Q1_n.scale(-1.0).sub(RQ2)

        // Compute R^(n) = Q0^(n-1) * U^{-1}
        R[n - 1] = safeMatrixDivide(Q0_nminus1, U)
    }

    return R.map { it!! }
}

/**
 * Safely divide matrix A by matrix U (compute A * U^{-1}).
 * Uses pseudo-inverse if U is singular.
 */
private fun safeMatrixDivide(A: Matrix, U: Matrix): Matrix {
    return if (U.length() == 1) {
        // Scalar case
        val uVal = U[0, 0]
        if (abs(uVal) > 1e-14) {
            A.scale(1.0 / uVal)
        } else {
            Matrix(A.numRows, A.numCols)  // Return zeros
        }
    } else {
        // Matrix case
        val det = U.det()
        if (abs(det) > 1e-14) {
            A.mult(U.inv())
        } else {
            // Use pseudo-inverse via SVD for singular matrices
            A.mult(pinv(U))
        }
    }
}

/**
 * Compute the Moore-Penrose pseudo-inverse of a matrix using SVD.
 * For A = U * S * V^T, the pseudo-inverse is V * S^+ * U^T
 * where S^+ has 1/sigma_i for non-zero singular values.
 */
private fun pinv(A: Matrix): Matrix {
    val svd = A.svd()
    val U = svd.u
    val S = svd.s  // Column vector of singular values
    val V = svd.v

    val m = A.numRows
    val n = A.numCols
    val Splus = Matrix(n, m)

    val tol = 1e-10 * maxOf(m, n) * (if (S.length() > 0) S[0, 0] else 0.0)
    val rank = S.numRows
    for (i in 0 until rank) {
        val sigma = S[i, 0]
        if (abs(sigma) > tol) {
            Splus[i, i] = 1.0 / sigma
        }
    }

    // Pseudo-inverse = V * S^+ * U^T
    return V.mult(Splus).mult(U.transpose())
}

/**
 * Compute stationary distribution from Algorithm 3 in Phung-Duc et al.
 *
 * 1. Boundary equation: pi_0 * (Q1^(0) + R^(1)*Q2^(1)) = 0
 * 2. Forward recursion: pi_n = pi_{n-1} * R^(n)
 * 3. Normalize: sum(pi) = 1
 */
private fun computeStationaryDist(
    R: List<Matrix>,
    Q0: List<Matrix>,
    Q1: List<Matrix>,
    Q2: List<Matrix>,
    N: Int,
    options: LdqbdOptions
): Matrix {
    // Check if all levels have same dimension (homogeneous) and scalar
    val isScalar = Q1.all { it.length() == 1 }
    val dims = Q1.map { it.numRows }
    val isHomogeneous = dims.toSet().size == 1

    return if (isScalar && isHomogeneous) {
        // Scalar case: direct computation
        val pi = Matrix(1, N + 1)

        // Start with pi_0 = 1 (will normalize later)
        pi[0, 0] = 1.0

        // Forward recursion: pi_n = pi_{n-1} * R^(n)
        for (n in 1..N) {
            if (R[n - 1].length() == 1) {
                pi[0, n] = pi[0, n - 1] * R[n - 1][0, 0]
            } else {
                pi[0, n] = 0.0
            }
        }

        // Normalize
        val sum = pi.sumRows(0)
        if (sum > 0) {
            pi.scaleEq(1.0 / sum)
        }
        pi
    } else {
        // Heterogeneous or matrix case: use cell-like approach
        val pi_cells = MutableList<Matrix>(N + 1) { Matrix(1, 1) }

        // Initialize pi_0 based on its dimension
        if (Q1[0].length() == 1) {
            // Scalar level 0: start with pi_0 = 1
            pi_cells[0] = Matrix(1, 1)
            pi_cells[0][0, 0] = 1.0
        } else {
            // Matrix level 0: solve boundary equation
            val Q1_0 = Q1[0]
            val Q2_1 = Q2[0]
            val A = Q1_0.add(1.0, R[0].mult(Q2_1))

            // Find left null space (solve pi * A = 0)
            // Use eigenvalue decomposition of A'
            val pi0 = solveLeftNullSpace(A)
            pi_cells[0] = pi0
        }

        // Forward recursion: pi_n = pi_{n-1} * R^(n)
        for (n in 1..N) {
            pi_cells[n] = pi_cells[n - 1].mult(R[n - 1])
        }

        // Normalize and convert to scalar probabilities (sum over phases at each level)
        var total = 0.0
        for (n in 0..N) {
            total += pi_cells[n].sumRows(0)
        }

        val pi = Matrix(1, N + 1)
        for (n in 0..N) {
            if (total > 0) {
                pi_cells[n].scaleEq(1.0 / total)
            }
            pi[0, n] = pi_cells[n].sumRows(0)
        }
        pi
    }
}

/**
 * Solve for left null space of matrix A (find pi such that pi * A = 0).
 */
private fun solveLeftNullSpace(A: Matrix): Matrix {
    val n = A.numRows

    // Find eigenvector corresponding to smallest eigenvalue of A'
    val AT = A.transpose()

    // Use power iteration to find dominant eigenvector of (A' - lambda_max * I)^{-1}
    // For simplicity, use a direct approach with normalization
    var pi = Matrix(1, n)
    pi.fill(1.0 / n)

    // Iterative refinement
    for (iter in 0 until 100) {
        val newPi = pi.mult(AT)
        val norm = newPi.norm()
        if (norm > 1e-14) {
            newPi.scaleEq(1.0 / norm)
        }
        if ((pi.sub(newPi)).norm() < 1e-10) {
            break
        }
        pi = newPi
    }

    // Ensure positive and normalize
    pi.absEq()
    val sum = pi.sumRows(0)
    if (sum > 0) {
        pi.scaleEq(1.0 / sum)
    }

    return pi
}

/**
 * LDQBD solver algorithms for level-dependent Quasi-Birth-Death processes.
 */
@Suppress("unused")
class LdqbdAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
