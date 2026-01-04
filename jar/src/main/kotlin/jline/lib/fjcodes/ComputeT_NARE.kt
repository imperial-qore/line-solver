package jline.lib.fjcodes

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.EigenDecomposition
import kotlin.math.abs
import kotlin.math.min

/**
 * Compute T-matrix using NARE (Nonsymmetric Algebraic Riccati Equation) method
 *
 * Solves the matrix equation for T via Schur decomposition of the Hamiltonian matrix.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Compute T-matrix via NARE method using Schur decomposition
 *
 * The T-matrix satisfies the equation:
 * T*X + X*kron(I_ms, D0) + kron(A_jump, I_ma) = 0
 *
 * where X is the solution to the algebraic Riccati equation associated with
 * the Hamiltonian matrix H = [A, B; -C, -D].
 *
 * @param D0 D0 matrix (transitions without arrivals)
 * @param D1 D1 matrix (transitions with arrivals)
 * @param S State transition matrix
 * @param A_jump Jump matrix
 * @return T-matrix
 */
fun computeT_NARE(D0: Matrix, D1: Matrix, S: Matrix, A_jump: Matrix): Matrix {
    val m = S.getNumRows()
    val ma = D0.getNumRows()
    val ms = m / ma

    // Debug: print input matrices
//    println("=== computeT_NARE inputs ===")
//    println("m = $m, ma = $ma, ms = $ms")
//    println("D0 (${D0.getNumRows()}x${D0.getNumCols()}): ${D0}")
//    println("D1 (${D1.getNumRows()}x${D1.getNumCols()}): ${D1}")
//    println("S (${S.getNumRows()}x${S.getNumCols()}) - first 3x3:")
//    for (i in 0 until minOf(3, S.getNumRows())) {
//        print("  [")
//        for (j in 0 until minOf(3, S.getNumCols())) {
//            print("%.6f ".format(S.get(i, j)))
//        }
//        println("...]")
//    }
//    println("A_jump (${A_jump.getNumRows()}x${A_jump.getNumCols()}): ${A_jump}")

    // Build Hamiltonian matrix H = [A, B; -C, -D]
    // where A = kron(I_ms, D0), B = kron(I_ms, D1), C = kron(A_jump, I_ma), D = S
    val Ims = Matrix.eye(ms)
    val Ima = Matrix.eye(ma)

    val A = Ims.kron(D0)
    val B = Ims.kron(D1)
    val C = A_jump.kron(Ima)
    val D = S

    // H = [A, B; -C, -D]
    val H = Matrix(2 * m, 2 * m)
    setSubMatrix(H, 0, 0, A)
    setSubMatrix(H, 0, m, B)
    setSubMatrix(H, m, 0, C.scale(-1.0))
    setSubMatrix(H, m, m, D.scale(-1.0))

    // Use Schur decomposition via reflection (SchurTransformer is package-private)
    val H_commons = Array2DRowRealMatrix(H.toArray2D())

    // Access package-private SchurTransformer via reflection
    val schurClass = Class.forName("org.apache.commons.math3.linear.SchurTransformer")
    val constructor = schurClass.getDeclaredConstructor(org.apache.commons.math3.linear.RealMatrix::class.java)
    constructor.isAccessible = true
    val schurTransformer = constructor.newInstance(H_commons)

    // Get T and P matrices via reflection
    val getT = schurClass.getDeclaredMethod("getT")
    getT.isAccessible = true
    val T_commons = getT.invoke(schurTransformer) as org.apache.commons.math3.linear.RealMatrix

    val getP = schurClass.getDeclaredMethod("getP")
    getP.isAccessible = true
    val P_commons = getP.invoke(schurTransformer) as org.apache.commons.math3.linear.RealMatrix

    // Convert to LINE Matrix format
    val T = Matrix(T_commons.data)
    val U = Matrix(P_commons.data)

    // Extract eigenvalues from Schur form
    val eigenvalues = extractSchurEigenvalues(T)

    // Build list of (real part of eigenvalue, index) pairs
    val eigenPairs = mutableListOf<Pair<Double, Int>>()
    for (i in 0 until 2 * m) {
        eigenPairs.add(Pair(eigenvalues[i].real, i))
    }

    // Sort by real part (ascending = most stable first)
    eigenPairs.sortBy { it.first }

    // Create selection vector: mark first m eigenvalues for top-left block
    val sel = BooleanArray(2 * m) { false }
    for (i in 0 until m) {
        val eigIndex = eigenPairs[i].second
        sel[eigIndex] = true
    }

    // Reorder Schur decomposition to move selected eigenvalues to top-left
    val Q1 = orderedSchur(U, T, sel)

    // Extract solution X = Q1(m+1:2m, 1:m) / Q1(1:m, 1:m)
    val Q1_11 = Matrix.getSubMatrix(Q1, 0, m, 0, m)
    val Q1_21 = Matrix.getSubMatrix(Q1, m, 2 * m, 0, m)

    val X = Q1_21.rightMatrixDivide(Q1_11)

    // Compute T2 = S + X * kron(I_ms, D1)
    val T2 = S.add(1.0, X.mult(Ims.kron(D1)))

    // Compute residual norm for verification
    // Equation: T*X + X*A + C = 0 where A = kron(I_ms, D0), C = kron(A_jump, I_ma)
    val termTX = T2.mult(X)
    val termXA = X.mult(Ims.kron(D0))
    val termC = A_jump.kron(Ima)

    val residual = termTX.add(1.0, termXA).add(1.0, termC)
    val residualNorm = normInf(residual)

    if (residualNorm > 1e-6) {
//        System.err.println("WARNING: NARE T-matrix residual = $residualNorm (expected < 1e-6)")
    }

    return T2
}

/**
 * Extract eigenvalues from Schur form (real or quasi-triangular)
 *
 * For real Schur form, eigenvalues are either on the diagonal (real) or
 * in 2x2 blocks (complex conjugate pairs).
 *
 * @param Q Schur matrix (upper triangular or quasi-triangular)
 * @return List of eigenvalues as Complex numbers
 */
private fun extractSchurEigenvalues(Q: Matrix): List<Complex> {
    val n = Q.getNumRows()
    val eigenvalues = mutableListOf<Complex>()

    var i = 0
    while (i < n) {
        if (i < n - 1 && abs(Q.get(i + 1, i)) > 1e-10) {
            // 2x2 block - complex conjugate pair
            val a = Q.get(i, i)
            val b = Q.get(i, i + 1)
            val c = Q.get(i + 1, i)
            val d = Q.get(i + 1, i + 1)

            // Eigenvalues of 2x2 matrix: lambda = (a+d)/2 ± sqrt((a+d)^2/4 - (ad-bc))
            val trace = a + d
            val det = a * d - b * c
            val discriminant = trace * trace / 4.0 - det

            if (discriminant >= 0) {
                // Real eigenvalues
                val sqrtDisc = kotlin.math.sqrt(discriminant)
                eigenvalues.add(Complex(trace / 2.0 + sqrtDisc, 0.0))
                eigenvalues.add(Complex(trace / 2.0 - sqrtDisc, 0.0))
            } else {
                // Complex eigenvalues
                val sqrtDisc = kotlin.math.sqrt(-discriminant)
                eigenvalues.add(Complex(trace / 2.0, sqrtDisc))
                eigenvalues.add(Complex(trace / 2.0, -sqrtDisc))
            }
            i += 2
        } else {
            // Real eigenvalue on diagonal
            eigenvalues.add(Complex(Q.get(i, i), 0.0))
            i++
        }
    }

    return eigenvalues
}

/**
 * Reorder Schur decomposition based on selection vector
 *
 * Uses QR-based algorithm to reorder eigenvalues while maintaining orthogonality.
 * Based on the LAPACK DTRSEN/DTREXC algorithms.
 *
 * @param U Orthogonal matrix from Schur decomposition
 * @param T Schur matrix (upper triangular or quasi-triangular)
 * @param sel Selection vector (true = move to top-left block)
 * @return Reordered orthogonal matrix U1
 */
private fun orderedSchur(U: Matrix, T: Matrix, sel: BooleanArray): Matrix {
    val n = U.getNumRows()

    // Make copies to avoid modifying inputs
    val U1 = U.copy()
    val T1 = T.copy()

    // Bubble sort: move selected eigenvalues to the top
    var changed = true
    var iter = 0
    val maxIter = n * n  // Prevent infinite loops

    while (changed && iter < maxIter) {
        changed = false
        iter++

        for (i in 0 until n - 1) {
            // If position i is not selected but position i+1 is selected,
            // swap them to move the selected one upward
            if (!sel[i] && sel[i + 1]) {
                // Swap eigenvalues at positions i and i+1
                swapAdjacentEigenvalues(U1, T1, i)

                // Swap the selection flags
                val temp = sel[i]
                sel[i] = sel[i + 1]
                sel[i + 1] = temp

                changed = true
            }
        }
    }

    return U1
}

/**
 * Swap two adjacent 1x1 blocks in real Schur form
 *
 * Based on LAPACK DLAEXC algorithm. Swaps diagonal elements at positions k and k+1
 * while maintaining the Schur decomposition property.
 *
 * @param U Orthogonal matrix (modified in-place)
 * @param T Schur matrix (modified in-place)
 * @param k Position of first element to swap
 */
private fun swapAdjacentEigenvalues(U: Matrix, T: Matrix, k: Int) {
    val n = T.getNumRows()

    // Extract the 2x2 block
    val t11 = T.get(k, k)
    val t22 = T.get(k + 1, k + 1)

    // Check if subdiagonal is already zero
    if (abs(T.get(k + 1, k)) < 1e-14) {
        // Already upper triangular, just need to swap diagonal elements
        // This requires rotating to swap them
        val t12 = T.get(k, k + 1)

        if (abs(t12) < 1e-14) {
            // Diagonal matrix, just swap the elements
            T.set(k, k, t22)
            T.set(k + 1, k + 1, t11)
            return
        }
    }

    // Compute Givens rotation to perform the swap
    // We want to zero out the (k+1, k) element after rotation
    // The rotation should also swap the eigenvalues

    val t12 = T.get(k, k + 1)

    // Standard Givens rotation to swap elements in a 2x2 upper triangular matrix
    // We solve for the rotation angle that exchanges eigenvalues
    val cs: Double
    val sn: Double

    if (abs(t11 - t22) < 1e-14) {
        // Eigenvalues are equal
        cs = 1.0
        sn = 0.0
    } else {
        // Compute rotation angle
        // For 2x2 matrix [t11, t12; 0, t22], we want rotation that makes
        // the result have t22 in upper-left and t11 in lower-right
        val sign = if (t22 - t11 > 0) 1.0 else -1.0
        val temp = (t22 - t11) / (2.0 * t12)
        val tau = sign / (abs(temp) + kotlin.math.sqrt(1.0 + temp * temp))
        cs = 1.0 / kotlin.math.sqrt(1.0 + tau * tau)
        sn = tau * cs
    }

    // Apply Givens rotation G = [cs, sn; -sn, cs]
    // Update T := G^T * T * G

    // First apply G from the right to columns k:k+1: T := T * G
    for (i in 0 until n) {
        val temp1 = T.get(i, k)
        val temp2 = T.get(i, k + 1)
        T.set(i, k, cs * temp1 + sn * temp2)
        T.set(i, k + 1, -sn * temp1 + cs * temp2)
    }

    // Then apply G^T from the left to rows k:k+1: T := G^T * T
    for (j in 0 until n) {
        val temp1 = T.get(k, j)
        val temp2 = T.get(k + 1, j)
        T.set(k, j, cs * temp1 + sn * temp2)
        T.set(k + 1, j, -sn * temp1 + cs * temp2)
    }

    // Apply same rotation to U from the right: U := U * G
    for (i in 0 until n) {
        val temp1 = U.get(i, k)
        val temp2 = U.get(i, k + 1)
        U.set(i, k, cs * temp1 + sn * temp2)
        U.set(i, k + 1, -sn * temp1 + cs * temp2)
    }
}

/**
 * Compute infinity norm (max row sum of absolute values)
 *
 * @param A Matrix to compute norm of
 * @return Infinity norm
 */
private fun normInf(A: Matrix): Double {
    var maxSum = 0.0
    for (i in 0 until A.getNumRows()) {
        var rowSum = 0.0
        for (j in 0 until A.getNumCols()) {
            rowSum += kotlin.math.abs(A.get(i, j))
        }
        maxSum = kotlin.math.max(maxSum, rowSum)
    }
    return maxSum
}
