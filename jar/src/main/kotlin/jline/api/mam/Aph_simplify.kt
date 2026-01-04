/**
 * @file Absorbing Phase-type distribution simplification and combination
 * 
 * Simplifies and combines APH distributions using structural pattern operations.
 * Used for building complex phase-type distributions from simpler components.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Simplifies and combines two APH distributions using different structural patterns.
 *
 * This function provides various ways to combine two Acyclic Phase-type (APH) distributions
 * represented by their initial probability vectors and generator matrices.
 *
 * @param a1 Initial probability vector for the first APH
 * @param T1 Generator matrix for the first APH
 * @param a2 Initial probability vector for the second APH
 * @param T2 Generator matrix for the second APH
 * @param p1 Branch probability for the first APH (used in branch pattern)
 * @param p2 Branch probability for the second APH (used in branch pattern)
 * @param pattern The combination pattern: 1=sequence, 2=parallel, 3=branch
 * @return A pair of (alpha, T) representing the combined APH distribution
 */
fun aph_simplify(a1: Matrix, T1: Matrix, a2: Matrix, T2: Matrix, p1: Double, p2: Double, pattern: Int): Pair<Matrix, Matrix> {
    return when (pattern) {
        1 -> sequencePattern(a1, T1, a2, T2)
        2 -> parallelPattern(a1, T1, a2, T2)
        3 -> branchPattern(a1, T1, a2, T2, p1, p2)
        else -> throw IllegalArgumentException("Pattern must be 1 (sequence), 2 (parallel), or 3 (branch)")
    }
}

/**
 * Sequence structure: first APH followed by second APH.
 */
private fun sequencePattern(a1: Matrix, T1: Matrix, a2: Matrix, T2: Matrix): Pair<Matrix, Matrix> {
    val order1 = a1.numCols
    val order2 = a2.numCols
    val e1 = Matrix.ones(order1, 1)
    
    // Compute alpha = [a1, (1-a1*e1)*a2]
    val a1_e1 = a1.mult(e1)
    val one_minus_a1_e1 = Matrix.ones(1, 1).sub(1.0, a1_e1)
    val alpha_right = one_minus_a1_e1.mult(a2)
    val alpha = Matrix.concatColumns(a1, alpha_right, null)
    
    // Compute T = [T1, (-T1*e1)*a2; zeros(order2, order1), T2]
    val minusT1_e1 = T1.scale(-1.0).mult(e1)
    val T_top_right = minusT1_e1.mult(a2)
    val T_bottom_left = Matrix.zeros(order2, order1)
    
    val T_top = Matrix.concatColumns(T1, T_top_right, null)
    val T_bottom = Matrix.concatColumns(T_bottom_left, T2, null)
    val T = Matrix.concatRows(T_top, T_bottom, null)
    
    return Pair(alpha, T)
}

/**
 * Parallel structure: both APHs run in parallel.
 */
private fun parallelPattern(a1: Matrix, T1: Matrix, a2: Matrix, T2: Matrix): Pair<Matrix, Matrix> {
    val order1 = a1.numCols
    val order2 = a2.numCols
    val e1 = Matrix.ones(order1, 1)
    val e2 = Matrix.ones(order2, 1)
    
    // Compute alpha = [kron(a1,a2), (1-a2*e2)*a1, (1-a1*e1)*a2]
    val kron_a1_a2 = kroneckerProduct(a1, a2)
    val one_minus_a2_e2 = Matrix.ones(1, 1).sub(1.0, a2.mult(e2))
    val alpha_middle = one_minus_a2_e2.mult(a1)
    val one_minus_a1_e1 = Matrix.ones(1, 1).sub(1.0, a1.mult(e1))
    val alpha_right = one_minus_a1_e1.mult(a2)
    
    val alpha_temp = Matrix.concatColumns(kron_a1_a2, alpha_middle, null)
    val alpha = Matrix.concatColumns(alpha_temp, alpha_right, null)
    
    // Compute T matrix components
    val eye_order1 = Matrix.eye(order1)
    val eye_order2 = Matrix.eye(order2)
    val kron_T1_eye2 = kroneckerProduct(T1, eye_order2)
    val kron_eye1_T2 = kroneckerProduct(eye_order1, T2)
    val kron_eye1_minusT2_e2 = kroneckerProduct(eye_order1, T2.scale(-1.0).mult(e2))
    val kron_minusT1_e1_eye2 = kroneckerProduct(T1.scale(-1.0).mult(e1), eye_order2)
    
    val Tr1_left = kron_T1_eye2.add(1.0, kron_eye1_T2)
    val Tr1_middle = kron_eye1_minusT2_e2
    val Tr1_right = kron_minusT1_e1_eye2
    val Tr1_temp = Matrix.concatColumns(Tr1_left, Tr1_middle, null)
    val Tr1 = Matrix.concatColumns(Tr1_temp, Tr1_right, null)
    
    val Tr2_temp = Matrix.concatColumns(Matrix.zeros(order1, order1 * order2), T1, null)
    val Tr2 = Matrix.concatColumns(Tr2_temp, Matrix.zeros(order1, order2), null)
    
    val Tr3_temp = Matrix.concatColumns(Matrix.zeros(order2, order1 * order2), Matrix.zeros(order2, order1), null)
    val Tr3 = Matrix.concatColumns(Tr3_temp, T2, null)
    
    val T_temp = Matrix.concatRows(Tr1, Tr2, null)
    val T = Matrix.concatRows(T_temp, Tr3, null)
    
    return Pair(alpha, T)
}

/**
 * Branch structure: select between two APHs with given probabilities.
 */
private fun branchPattern(a1: Matrix, T1: Matrix, a2: Matrix, T2: Matrix, p1: Double, p2: Double): Pair<Matrix, Matrix> {
    val order1 = a1.numCols
    val order2 = a2.numCols
    
    // Compute alpha = [p1*a1, p2*a2]
    val alpha_left = a1.scale(p1)
    val alpha_right = a2.scale(p2)
    val alpha = Matrix.concatColumns(alpha_left, alpha_right, null)
    
    // Compute T = [T1, zeros(order1, order2); zeros(order2, order1), T2]
    val T_top = Matrix.concatColumns(T1, Matrix.zeros(order1, order2), null)
    val T_bottom = Matrix.concatColumns(Matrix.zeros(order2, order1), T2, null)
    val T = Matrix.concatRows(T_top, T_bottom, null)
    
    return Pair(alpha, T)
}

/**
 * Computes the Kronecker product of two matrices.
 */
private fun kroneckerProduct(A: Matrix, B: Matrix): Matrix {
    val result = Matrix(A.numRows * B.numRows, A.numCols * B.numCols)
    
    for (i in 0 until A.numRows) {
        for (j in 0 until A.numCols) {
            val aij = A.get(i, j)
            for (k in 0 until B.numRows) {
                for (l in 0 until B.numCols) {
                    result.set(i * B.numRows + k, j * B.numCols + l, aij * B.get(k, l))
                }
            }
        }
    }
    
    return result
}
/**
 * APH  simplify algorithms
 */
@Suppress("unused")
class AphSimplifyAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}